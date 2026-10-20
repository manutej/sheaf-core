"""sheaf-core ADR-016: the Rust-backed m1 shim (python/m1/sheafimpl_rs.py + libsheaf_cli.so) vs the Python reference
(python/sheafcore.py) through the PUBLIC function surface, in-process via ctypes. Compares values (stated tolerances as in
crates/sheaf-cli/tests/m1_parity.rs) and, on a battery of invalid inputs, the exception type and message.
usage: cargo build --release && cp target/release/libsheaf_cli.so python/m1/ && <venv>/python python/m1/rust_parity.py"""
import json, pathlib, sys
import numpy as np, scipy.sparse as sp
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "python")); sys.path.insert(0, str(ROOT / "python/m1"))
import sheafcore as P, sheafimpl_rs as R
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
W = dict(delta=0.0, L=0.0, discord=0.0, energy=0.0, proj=0.0); bad = []; n_val = n_err = 0

def err(f, *a):
    try: f(*a); return None
    except Exception as e: return (type(e).__name__, str(e))
def cmp(name, dv, edges, x):
    global n_val
    n_val += 1; a, b = err(P.coboundary, dv, edges), err(R.coboundary, dv, edges)
    if a is not None or b is not None:
        if a != b: bad.append((name, "coboundary error parity", a, b))
        return
    Dp, Dr = P.coboundary(dv, edges), R.coboundary(dv, edges)
    if Dp.shape != Dr.shape or Dp.nnz != Dr.nnz: bad.append((name, "delta shape/nnz", Dp.shape, Dr.shape, Dp.nnz, Dr.nnz)); return
    W["delta"] = max(W["delta"], abs(Dp - Dr).max() if Dp.nnz else 0.0)
    Lp, Lr = P.laplacian(Dp), R.laplacian(Dr); e = abs(Lp - Lr).max() if Lp.nnz or Lr.nnz else 0.0
    W["L"] = max(W["L"], e / max(1.0, abs(Lp).max() if Lp.nnz else 0.0))
    if Lp.nnz != Lr.nnz: bad.append((name, "L nnz", Lp.nnz, Lr.nnz))
    de = [int(np.asarray(E[2], float).shape[0]) for E in edges]
    for op in (Dp, Dp.toarray(), Dp.tocsc(), Dp.tocoo()):
        dp, dr = P.edge_discord(op, x, de), R.edge_discord(op, x, de)
        if not np.all(np.abs(dp - dr) <= 1e-9 * np.abs(dp) + 1e-12): bad.append((name, "discord"))
        W["discord"] = max(W["discord"], float(np.max(np.abs(dp - dr) / np.maximum(np.abs(dp), 1e-300), initial=0.0)))
    tol = 1e-9 * float(dp.sum()) + 1e-12 * (abs(Lp).sum(axis=1).max() if Lp.nnz else 0.0) * float(x @ x)
    for op in (Lp, Lp.toarray(), Lp.tocsc()):
        ep, er = P.energy(op, x), R.energy(op, x); W["energy"] = max(W["energy"], abs(ep - er))
        if abs(ep - er) > tol: bad.append((name, "energy", ep, er))
    Gp, Gr = P.global_sections(Dp), R.global_sections(Dr)
    if Gp.shape != Gr.shape: bad.append((name, "H0 svd", Gp.shape, Gr.shape)); return
    if Gp.size: W["proj"] = max(W["proj"], float(abs(Gp @ Gp.T - Gr @ Gr.T).max()))
    if Gp.size:                                          # column-block order: column j lives on the same component in both
        lab = np.full(Dp.shape[1], -1)
        for i, cols in enumerate(P._components(Dp, Dp.shape[1])): lab[cols] = i
        cp = [lab[np.flatnonzero(np.abs(Gp[:, j]) > 1e-9)[0]] for j in range(Gp.shape[1])]
        cr = [lab[np.flatnonzero(np.abs(Gr[:, j]) > 1e-9)[0]] for j in range(Gr.shape[1])]
        if cp != cr: bad.append((name, "G column order", cp, cr))
    for t in (-1.0, 0.0, 1e-3):                          # caller tolerances, incl. negative (Grok R1 critical)
        if P.global_sections(Dp, t).shape != R.global_sections(Dr, t).shape: bad.append((name, "H0 svd tol", t))
        if t != 0.0 and P.h0_dim_from_L(Lp, t) != R.h0_dim_from_L(Lr, t): bad.append((name, "H0 eig tol", t))   # tol=0 on λ is the sign of round-off (±1e-16), LAPACK- vs faer-dependent: not compared
    if P.h0_dim_from_L(Lp) != R.h0_dim_from_L(Lr): bad.append((name, "H0 eig"))
    if P.h0_dim_from_L(Lp.toarray()) != R.h0_dim_from_L(Lp.toarray()): bad.append((name, "H0 eig dense"))

fx = json.loads((ROOT / "fixtures/m1/rust_parity.json").read_text())
for c in fx["cases"]:
    dv = c["dv"]; E = [(e["u"], e["v"], np.asarray(e["fu"]).reshape(e["de"], dv[e["u"]]), np.asarray(e["fv"]).reshape(e["de"], dv[e["v"]])) for e in c["edges"]]
    cmp(c["name"], dv, E, np.asarray(c["x"])); cmp(c["name"] + "/lists", dv, [(u, v, a.tolist(), b.tolist()) for u, v, a, b in E], np.asarray(c["x"]))
dv, E, x = P.from_spec(json.loads((ROOT / "fixtures/m1/spec.json").read_text())); cmp("spec", dv, E, x)
# hand-built operators: duplicates, explicit zeros, unsorted indices, all-zero rows (pattern semantics)
for k in range(30):
    m, n = int(rng.integers(1, 9)), int(rng.integers(1, 9)); nnz = int(rng.integers(0, 3 * m * n // 2 + 1))
    r, cc = rng.integers(0, m, nnz), rng.integers(0, n, nnz); v = rng.normal(size=nnz); v[rng.random(nnz) < 0.2] = 0.0
    A = sp.csr_matrix((v, (r, cc)), shape=(m, n)); A.has_canonical_format = False   # explicit zeros kept, then rows reversed
    for i in range(m): s, t = A.indptr[i], A.indptr[i + 1]; A.indices[s:t] = A.indices[s:t][::-1]; A.data[s:t] = A.data[s:t][::-1]
    Gp, Gr = P.global_sections(A), R.global_sections(A); n_val += 1
    if Gp.shape != Gr.shape or (Gp.size and abs(Gp @ Gp.T - Gr @ Gr.T).max() > 1e-8): bad.append((f"hand{k}", "H0 svd", Gp.shape, Gr.shape))
    Lp, Lr = P.laplacian(A), R.laplacian(A)
    if abs(Lp - Lr).max() > 1e-12 * max(1, abs(Lp).max()) if Lp.nnz else Lr.nnz: bad.append((f"hand{k}", "L"))
    S = rng.normal(size=(n, n)); S[rng.random((n, n)) < 0.5] = 0.0
    if P.h0_dim_from_L(S) != R.h0_dim_from_L(S): bad.append((f"hand{k}", "h0 nonsym"))
# raw CSR with unmerged duplicate (i, j) entries (incl. pairs that cancel) and L cancellation (Grok R1 minor/major)
for k in range(20):
    m, n = int(rng.integers(1, 6)), int(rng.integers(1, 6)); cnt = rng.integers(0, 2 * n + 1, m)
    ind = np.concatenate([rng.integers(0, n, c) for c in cnt]) if cnt.sum() else np.zeros(0, int)
    dat = rng.choice([-1.0, 1.0, 0.0, 2.5], ind.size)
    A = sp.csr_matrix((dat, ind, np.r_[0, np.cumsum(cnt)]), shape=(m, n)); x = rng.normal(size=n); n_val += 1
    Gp, Gr = P.global_sections(A), R.global_sections(A)
    if Gp.shape != Gr.shape or (Gp.size and abs(Gp @ Gp.T - Gr @ Gr.T).max() > 1e-8): bad.append((f"dup{k}", "H0 svd"))
    Lp, Lr = P.laplacian(A), R.laplacian(A)
    if Lp.nnz != Lr.nnz or (Lp.nnz and abs(Lp - Lr).max() > 1e-12 * max(1, abs(Lp).max())): bad.append((f"dup{k}", "L", Lp.nnz, Lr.nnz))
    if P.h0_dim_from_L(A.T @ A) != R.h0_dim_from_L(A.T @ A): bad.append((f"dup{k}", "h0"))
    dp, dr = P.edge_discord(A, x, [1] * m), R.edge_discord(A, x, [1] * m)
    if not np.all(np.abs(dp - dr) <= 1e-9 * np.abs(dp) + 1e-12): bad.append((f"dup{k}", "discord"))
    if abs(P.energy(A.T @ A, x) - R.energy(A.T @ A, x)) > 1e-9 * max(1.0, float(dp.sum())): bad.append((f"dup{k}", "energy"))
cmp("cancel", [1, 1], [(0, 1, [[1.0]], [[1.0]]), (0, 1, [[1.0]], [[-1.0]])], np.array([0.3, -0.7]))
cmp("self_loop", [2, 1], [(0, 0, [[1.0, 2.0]], [[1.0, 0.5]]), (0, 1, [[0.0, 1.0]], [[2.0]]), (1, 1, [[3.0]], [[3.0]])], np.array([1.0, -2.0, 0.5]))
# error-path parity: same exception type and message
Z = np.zeros((2, 2)); I1 = [[1.0]]
BAD = [("coboundary", ([1, 1], [(0, 1, [[np.inf]], I1)])), ("coboundary", ([1, 1], [(0, 2, I1, I1)])), ("coboundary", ([1, -1], [])),
       ("coboundary", ([1, 1], [(0, 1, [[1.0, 2.0]], I1)])), ("coboundary", ([1, 1], [(0.0, 1, I1, I1)])), ("coboundary", ([2, 2], [(0, 1, Z, np.full((2, 2), np.nan))])),
       ("coboundary", ([1, 1], [(0, 1, sp.csr_matrix(I1), I1)])), ("coboundary", ([1, 1], [(0, 1, I1)])),
       ("laplacian", ([1, 2],)), ("laplacian", (np.zeros(3),)), ("edge_discord", (sp.csr_matrix(Z), np.zeros(2), [1, 2])),
       ("edge_discord", (sp.csr_matrix(Z), np.zeros(3), [1, 1])), ("edge_discord", (sp.csr_matrix(Z), [np.nan, 0], [1, 1])),
       ("edge_discord", (sp.csr_matrix([[np.inf, 0], [0, 0]]), np.zeros(2), [1, 1])), ("edge_discord", (sp.csr_matrix([[1e300, 0]]), [1e10, 0], [1])),
       ("energy", (sp.csr_matrix(Z), np.zeros(3))), ("energy", (sp.csr_matrix([[np.nan, 0], [0, 0]]), np.zeros(2))), ("energy", (np.full((1, 1), 1e300), [1e10])),
       ("energy", (np.zeros((2, 3)), np.zeros(3))), ("global_sections", (sp.csr_matrix([[np.inf]]),)), ("global_sections", (np.zeros((2, 2, 2)),)),
       ("global_sections", (sp.csr_matrix((np.r_[np.ones(5999), -np.ones(5999)], (np.r_[np.arange(5999), np.arange(5999)], np.r_[np.arange(5999), np.arange(1, 6000)])), shape=(5999, 6000)),)),
       ("h0_dim_from_L", (sp.csr_matrix((2, 3)),)), ("h0_dim_from_L", (sp.csr_matrix([[np.nan]]),)),
       ("h0_dim_from_L", (sp.diags(np.ones(5999), 1, shape=(6000, 6000)).tocsr(),))]
for fn, args in BAD:
    n_err += 1; a, b = err(getattr(P, fn), *args), err(getattr(R, fn), *args)
    if a != b or a is None: bad.append((fn, "error parity", a, b))
print(json.dumps(dict(value_cases=n_val, error_cases=n_err, max_abs_delta=W["delta"], max_rel_L=W["L"], max_rel_discord=W["discord"],
                      max_abs_energy=W["energy"], max_abs_projector=W["proj"], mismatches=len(bad))))
for b in bad: print("MISMATCH", b)
sys.exit(1 if bad else 0)
