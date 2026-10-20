"""sheaf-core m2-native public parity run (reference-only tooling; numpy + scipy).
Materialises a CLOSED-WORLD native builder dir exactly as SPEC_ADDENDUM R-M2-LOAD lists it -- sheafimpl.py (this repo's
python/m2/sheafimpl.py), sheaf_native.py (the contract-owned pinned loader, copied byte for byte) and libsheaf_cli.so
(target/release, which must hash to the loader's PINNED_SHA256) -- and nothing else; imports the shim the way the m1
runner imports a builder (importlib spec_from_file_location); then:
  1. pin checks: exact file set, library sha256 == PINNED_SHA256, the only newly mapped sheaf object is that library;
  2. loader refusals: a 1-byte-tampered copy and a symlinked library both raise ImportError (no fallback);
  3. the full public parity suite python/m1/rust_parity.py against this shim (values, errors, observable dtype/type).
usage: cargo build --release && <venv>/python python/m2/m2_parity.py [--loader /path/to/pinned/sheaf_native.py]
NOTE (open ambiguity, see the report): the shim imports the loader as a plain top-level module, so the builder dir must
be on sys.path when the shim is imported; this script puts it there. The public addendum does not say how the gate does it."""
import argparse, hashlib, importlib.util, os, pathlib, re, runpy, shutil, sys, tempfile
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.dont_write_bytecode = True   # closed world: importing must not add __pycache__ to the builder dir
ap = argparse.ArgumentParser(); ap.add_argument("--loader", default="/home/box/frontier/contracts/m2-native/pinned/sheaf_native.py")
ap.add_argument("--lib", default=str(ROOT / "target/release/libsheaf_cli.so")); args, rest = ap.parse_known_args()
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
pin = re.search(r'^PINNED_SHA256 = "([0-9a-f]{64})"$', pathlib.Path(args.loader).read_text(), re.M).group(1)
FILES = {"sheafimpl.py": ROOT / "python/m2/sheafimpl.py", "sheaf_native.py": pathlib.Path(args.loader), "libsheaf_cli.so": pathlib.Path(args.lib)}
fails = []

def materialise(d, lib_bytes=None, symlink=False):
    for name, src in FILES.items():
        dst = pathlib.Path(d) / name
        if name == "libsheaf_cli.so" and symlink:
            real = pathlib.Path(d).parent / (pathlib.Path(d).name + ".real.so"); shutil.copyfile(src, real); os.symlink(real, dst)
        elif name == "libsheaf_cli.so" and lib_bytes is not None: dst.write_bytes(lib_bytes)
        else: shutil.copyfile(src, dst)
def load(d, modname):
    sys.path.insert(0, d)
    try:
        sys.modules.pop("sheaf_native", None)
        spec = importlib.util.spec_from_file_location(modname, os.path.join(d, "sheafimpl.py")); m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m); return m
    finally: sys.path.remove(d)

base = tempfile.mkdtemp(prefix="m2-parity-")
# 2. refusals first (each in its own dir; the dynamic loader caches by path, so distinct dirs)
raw = FILES["libsheaf_cli.so"].read_bytes(); flip = bytearray(raw); flip[len(flip) // 2] ^= 1
for tag, kw in (("tampered", dict(lib_bytes=bytes(flip))), ("symlink", dict(symlink=True))):
    d = os.path.join(base, tag); os.mkdir(d); materialise(d, **kw)
    try: load(d, "sheafimpl_" + tag); fails.append(f"{tag} library was NOT refused")
    except ImportError as e: print(f"refused {tag}: {e}")
sys.modules.pop("sheaf_native", None)
# 1. closed world + pin + mapped object
d = os.path.join(base, "builder"); os.mkdir(d); materialise(d)
if sorted(os.listdir(d)) != sorted(FILES): fails.append(f"file set {sorted(os.listdir(d))}")
if sha(os.path.join(d, "libsheaf_cli.so")) != pin: fails.append("library sha256 != PINNED_SHA256 (rebuild sheaf-core at the pinned commit)")
before = set(l.split()[-1] for l in open("/proc/self/maps") if "libsheaf_cli" in l)
R = load(d, "sheafimpl")
mapped = sorted(set(l.split()[-1] for l in open("/proc/self/maps") if "libsheaf_cli" in l) - before)
if mapped != [os.path.join(d, "libsheaf_cli.so")]: fails.append(f"newly mapped sheaf objects: {mapped}")
if sorted(os.listdir(d)) != sorted(FILES): fails.append(f"import changed the builder dir: {sorted(os.listdir(d))}")   # no __pycache__ written
print(f"pin ok: {sha(os.path.join(d, 'libsheaf_cli.so'))[:16]}... mapped {mapped}")
if fails: print("M2-PIN FAIL", fails); sys.exit(1)
# 3. the public parity suite against this shim
sys.modules["sheafimpl_rs"] = R
sys.argv = [str(ROOT / "python/m1/rust_parity.py")] + rest
runpy.run_path(sys.argv[0], run_name="__main__")
