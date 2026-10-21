"""Regression (public review finding, R-PERF-18): component sizes / caps must come from the ORIGINAL dtype's structural
nonzeros. A float64 cast drops imaginary parts (and longdouble values below the float64 range), so purely imaginary
couplings looked like zeros, components looked disconnected, the cap never fired and the reference-lines fallback then
densified the over-cap block. The shim must raise exactly the reference's prefixed MemoryError, and stay in behaviour
parity with the reference on small mixed real / imaginary inputs.
Materialises the same closed-world builder dir as m2_parity.py (shim + pinned loader + library) and imports the shim BY
FILE PATH (its §6(1) self-insert finds the loader); the reference is python/sheafcore.py.
usage: <venv>/python python/m2/test_m2_perf18_dtype.py [--loader /path/to/sheaf_native.py] [--lib target/release/libsheaf_cli.so]
(pytest also collects it; set M2_LOADER / M2_LIB to override the defaults)."""
import argparse, importlib.util, os, pathlib, shutil, sys, tempfile, warnings
import numpy as np, scipy.sparse as sp
ROOT = pathlib.Path(__file__).resolve().parents[2]
LOADER = os.environ.get("M2_LOADER", "/home/box/frontier/contracts/m2-native/pinned/sheaf_native.py")
LIB = os.environ.get("M2_LIB", str(ROOT / "target/release/libsheaf_cli.so"))
SHIM = os.environ.get("M2_SHIM", str(ROOT / "python/m2/sheafimpl.py"))
_mods = {}

def _by_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
def mods():
    if not _mods:
        d = tempfile.mkdtemp(prefix="m2-perf18-")
        for name, src in (("sheafimpl.py", SHIM), ("sheaf_native.py", LOADER), ("libsheaf_cli.so", LIB)):
            shutil.copyfile(src, os.path.join(d, name))
        sys.modules.pop("sheaf_native", None)
        _mods["shim"] = _by_path("sheafimpl_perf18", os.path.join(d, "sheafimpl.py"))
        _mods["ref"] = _by_path("sheafcore_ref_perf18", str(ROOT / "python/sheafcore.py"))
    return _mods["shim"], _mods["ref"]

def _outcome(fn, *a):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try: return ("ok", fn(*a))
        except Exception as e: return ("err", type(e), str(e))
def _same(fname, *a):
    S, R = mods(); r = _outcome(getattr(R, fname), *a); s = _outcome(getattr(S, fname), *a)
    assert s[0] == r[0], f"{fname}: shim {s[:3]} vs reference {r[:3]}"
    if s[0] == "err": assert s[1:] == r[1:], f"{fname}: shim {s[1:]} vs reference {r[1:]}"
    else:
        sv, rv = s[1], r[1]
        assert type(sv) is type(rv), f"{fname}: type {type(sv)} vs {type(rv)}"
        if isinstance(rv, np.ndarray):
            assert sv.dtype == rv.dtype and sv.shape == rv.shape, f"{fname}: {sv.dtype}{sv.shape} vs {rv.dtype}{rv.shape}"
            assert np.allclose(sv, rv, atol=1e-12, rtol=0), f"{fname}: values differ"
        else: assert sv == rv, f"{fname}: {sv} vs {rv}"
    return s

def _imag_row(n=5001, real_every=0):
    v = np.full(n, 1j)
    if real_every: v[::real_every] = 1.0
    return sp.csr_matrix(v.reshape(1, n))

def test_imag_row_over_cap_global_sections():          # the repro: 1 x 5001 purely imaginary row
    s = _same("global_sections", _imag_row())
    assert s[0] == "err" and s[1] is MemoryError and s[2].startswith("R-PERF-18:"), s
def test_mixed_row_over_cap_global_sections():         # alternating real / imaginary entries in one over-cap row
    s = _same("global_sections", _imag_row(real_every=2))
    assert s[0] == "err" and s[1] is MemoryError and s[2].startswith("R-PERF-18:"), s
def test_imag_star_over_cap_h0_dim_from_L():            # L: real diagonal, purely imaginary star coupling 5001 nodes
    n = 5001; r = np.r_[np.arange(n), np.zeros(n - 1, int), np.arange(1, n)]; c = np.r_[np.arange(n), np.arange(1, n), np.zeros(n - 1, int)]
    v = np.r_[np.full(n, 2.0 + 0j), np.full(n - 1, 1j), np.full(n - 1, -1j)]
    s = _same("h0_dim_from_L", sp.csr_matrix((v, (r, c)), shape=(n, n)))
    assert s[0] == "err" and s[1] is MemoryError and s[2].startswith("R-PERF-18:"), s
def test_mixed_small_parity():                         # within the caps: same values / dtypes / types as the reference
    D = sp.csr_matrix(np.array([[1.0, 1j, 0, 0, 0], [0, 0, 2.0, 0.5j, 0], [0, 0, 0, 1j, 0]]))
    _same("global_sections", D)
    _same("global_sections", D.toarray())
    _same("global_sections", sp.csr_matrix(np.array([[1j, -1j, 0], [0, 0, 0]])))   # imaginary-only coupling, below the caps
    L = sp.csr_matrix(np.array([[2.0, 1j, 0], [-1j, 2.0, 0], [0, 0, 0.0]]))
    _same("h0_dim_from_L", L); _same("h0_dim_from_L", L.toarray())
def test_longdouble_below_f64_range():                 # a value a float64 cast turns into 0 must still couple
    if np.finfo(np.longdouble).tiny >= 1e-300: return   # platform without extended precision: nothing to test
    t = np.longdouble("1e-400"); D = sp.csr_matrix(np.array([[t, t, 0], [0, 0, 1]], dtype=np.longdouble))
    _same("global_sections", D)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--loader", default=LOADER); ap.add_argument("--lib", default=LIB)
    ap.add_argument("--shim", default=SHIM); a = ap.parse_args(); LOADER, LIB, SHIM = a.loader, a.lib, a.shim
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]; bad = 0
    for k, f in tests:
        try: f(); print("PASS", k)
        except AssertionError as e: bad += 1; print("FAIL", k, str(e)[:300])
    print(f"{len(tests) - bad}/{len(tests)} passed"); sys.exit(1 if bad else 0)
