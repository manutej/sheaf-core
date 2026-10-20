"""REFERENCE-ONLY Python. Parity: Python reference sheafcore (python/sheafcore.py, m1 SPEC v1.18) vs Rust v1 core (sheaf-cli), on shared v1 spec fixtures.
Generates the fixtures (deterministic), runs `cargo run --example parity_dump`, compares coboundary / Laplacian /
energy / per-edge discord to 1e-9 (max abs, relative to max(1, |ref|max)) and H0 dims exactly.
usage: python python/parity/parity.py   (needs numpy, scipy; run from anywhere)"""
import json, subprocess, sys, os
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "python"))
import sheafcore as sc
FIX = HERE / "fixtures"; FIX.mkdir(exist_ok=True)
rng = np.random.default_rng(20261019)
def rot(t): c, s = np.cos(t), np.sin(t); return [[c, -s], [s, c]]
I2 = [[1.0, 0.0], [0.0, 1.0]]
fx = {}
fx["translation_loop_rot90"] = dict(nodes=[dict(id=i, dim=2) for i in ("en", "es", "fr")],
    edges=[dict(u="en", v="es", Fu=rot(np.pi/2), Fv=I2), dict(u="es", v="fr", Fu=rot(np.pi/2), Fv=I2), dict(u="en", v="fr", Fu=I2, Fv=rot(np.pi/2))],
    x=[[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]])
# config drift: timeout in each file's own unit, per-edge maps = unit factor to seconds / edge scale (as audit-local does)
units = dict(env=1.0, docker=1e-3, readme=60.0, k8s=1e-3); vals = dict(env=60.0, docker=60000.0, readme=1.0, k8s=60.0)
names = list(units); E = []
for a in range(4):
    for b in range(a + 1, 4):
        u, v = names[a], names[b]; sc_ = max(abs(vals[u] * units[u]), abs(vals[v] * units[v]))
        E.append(dict(u=u, v=v, Fu=[[units[u] / sc_]], Fv=[[units[v] / sc_]], label=f"{u}<->{v} timeout"))
fx["config_drift_units"] = dict(nodes=[dict(id=n, dim=1) for n in names], edges=E, x=[[vals[n]] for n in names])
dims = [1, 2, 3, 2, 1, 3]; E = []
for (u, v) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (1, 4), (2, 5)]:
    de = int(rng.integers(1, 4)); E.append(dict(u=u, v=v, Fu=rng.standard_normal((de, dims[u])).tolist(), Fv=rng.standard_normal((de, dims[v])).tolist()))
fx["hetero_random"] = dict(nodes=[dict(id=f"n{i}", dim=d) for i, d in enumerate(dims)], edges=E, x=[rng.standard_normal(d).tolist() for d in dims])
k4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
fx["k4_rot60"] = dict(nodes=[dict(id=f"v{i}", dim=2) for i in range(4)],
    edges=[dict(u=a, v=b, Fu=I2, Fv=(rot(np.pi/3) if (a, b) == (0, 1) else I2)) for a, b in k4], x=[rng.standard_normal(2).tolist() for _ in range(4)])
fx["triangle_w_only"] = dict(nodes=[dict(id=c, dim=1) for c in "abc"], edges=[dict(u="a", v="b"), dict(u="b", v="c"), dict(u="a", v="c")], w=[[0.0], [0.0], [1.0]])
paths = []
for name, f in fx.items():
    f = dict(format="frontier-sheaf/v1", **f, options=dict(null_reps=0)); p = FIX / f"{name}.json"; p.write_text(json.dumps(f, indent=1)); paths.append(str(p))
out = subprocess.run(["cargo", "run", "-q", "--release", "--offline", "-p", "sheaf-cli", "--example", "parity_dump", "--", *paths], cwd=ROOT, capture_output=True, text=True, check=True).stdout
R = {r["name"]: r for r in map(json.loads, out.strip().splitlines())}
TOL = 1e-9; worst = 0.0; ok = True; rows = []
def close(a, b, what, name):
    global worst, ok
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape: print(f"FAIL {name} {what}: shape {a.shape} vs {b.shape}"); ok = False; return float("inf")
    e = float(np.max(np.abs(a - b))) / max(1.0, float(np.max(np.abs(b))) if b.size else 1.0) if a.size else 0.0
    worst = max(worst, e)
    if e > TOL: print(f"FAIL {name} {what}: rel err {e:.3e}"); ok = False
    return e
for name, f in fx.items():
    ix = {n["id"]: i for i, n in enumerate(f["nodes"])}; dv = [n["dim"] for n in f["nodes"]]
    edges = []
    for e in f["edges"]:
        u, v = ix[e["u"]] if isinstance(e["u"], str) else e["u"], ix[e["v"]] if isinstance(e["v"], str) else e["v"]
        Fu = np.array(e.get("Fu", np.eye(dv[u])), float); Fv = np.array(e.get("Fv", np.eye(dv[v])), float); edges.append((u, v, Fu, Fv))
    D = sc.coboundary(dv, edges); L = sc.laplacian(D); r = R[name]
    e1 = close(np.array(r["delta"]).reshape(D.shape), D.toarray(), "coboundary", name)
    e2 = close(np.array(r["laplacian"]).reshape(L.shape), L.toarray(), "laplacian", name)
    h0_py_svd = int(sc.global_sections(D).shape[1]); h0_py_eig = int(sc.h0_dim_from_L(L))
    h0_ok = (h0_py_svd == h0_py_eig == r["h0"]); ok &= h0_ok
    e3 = e4 = 0.0
    if "x" in f:
        x = np.concatenate([np.asarray(z, float) for z in f["x"]]); de = [e[2].shape[0] for e in edges]
        e3 = close([r["energy"]], [sc.energy(L, x)], "energy", name); e4 = close(r["edge_discord"], sc.edge_discord(D, x, de), "edge_discord", name)
    rows.append((name, e1, e2, e3, e4, h0_py_svd, h0_py_eig, r["h0"], h0_ok))
print(f"{'fixture':26} {'delta':>9} {'L':>9} {'energy':>9} {'discord':>9}  h0 py(svd/eig) rust")
for n, a, b, c, d, hs, he, hr, hk in rows: print(f"{n:26} {a:9.1e} {b:9.1e} {c:9.1e} {d:9.1e}  {hs}/{he} {hr} {'OK' if hk else 'MISMATCH'}")
print(f"worst relative error {worst:.3e} (tol {TOL:g}) -> {'PASS' if ok else 'FAIL'}")
json.dump(dict(rows=rows, worst=worst, tol=TOL, passed=bool(ok)), open(HERE / "parity_result.json", "w"), indent=1)
sys.exit(0 if ok else 1)
