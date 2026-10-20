"""sheaf-core ADR-016: Python-reference outputs for the Rust m1 port's parity tests (reference-only; numpy + scipy).
Writes fixtures/m1/rust_parity.json:
  * "cases": seeded small heterogeneous sheaves (dim-0 stalks, multi-edges, isolated nodes, self-loops, explicit zero map
    entries, near-singular maps on both sides of the R-H0-08 cut with a wide gap, constant sheaves) with the reference's
    dense delta, dense L, per-edge discord, energy, dim H0 by both routes and the kernel projector G G^T;
  * "spec": the synthetic m1 spec (fixtures/m1/spec.json) as delta triplets + discord / energy / H0;
  * "table": every row of the regenerated shared mutation table (fixtures/m1/mutations.json, 114 rows) with the
    models verdict and, for accepted rows, the reference's discord / energy / H0 / nnz on the mutated spec.
usage: /path/to/venv/python python/m1/make_rust_parity.py   (uses python/sheafcore.py = the gated reference)"""
import copy, json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "python")); sys.path.insert(0, str(ROOT / "python/m1"))
import sheafcore as P, models
rng = np.random.default_rng(20261020)

def ref(dv, edges, x):
    D = P.coboundary(dv, edges); L = P.laplacian(D); de = [int(np.asarray(e[2], float).shape[0]) for e in edges]
    G = P.global_sections(D)
    return D, L, dict(discord=P.edge_discord(D, x, de).tolist(), energy=P.energy(L, x), h0_svd=int(G.shape[1]),
                      h0_eig=int(P.h0_dim_from_L(L)), nnz_delta=int(D.nnz), nnz_L=int(L.nnz)), G

def rand_case(k):
    nv = int(rng.integers(1, 7)); dv = [int(rng.integers(0, 4)) for _ in range(nv)]
    if k % 3 == 0: dv = [2] * nv                                   # uniform fast path
    edges = []
    for _ in range(int(rng.integers(0, 9))):
        u, v = int(rng.integers(0, nv)), int(rng.integers(0, nv))
        if k % 4 != 1 and u == v: continue                         # self-loops only in some cases
        de = dv[u] if k % 3 == 0 else int(rng.integers(0, 4))
        Fu = rng.normal(size=(de, dv[u])); Fv = rng.normal(size=(de, dv[v]))
        Fu[rng.random(Fu.shape) < 0.3] = 0.0; Fv[rng.random(Fv.shape) < 0.3] = 0.0   # explicit zero entries
        if rng.random() < 0.3 and de and dv[u] and dv[v] and dv[u] == dv[v]: Fv = Fu.copy() if u != v else Fv   # restriction agreement -> sections
        edges.append((u, v, Fu.tolist(), Fv.tolist()))
        if rng.random() < 0.25: edges.append((u, v, Fu.tolist(), Fv.tolist()))   # multi-edge
    return dv, edges

cases = []
for k in range(60): cases.append(("rand%02d" % k,) + rand_case(k))
n = 6; I2 = np.eye(2).tolist()
cases.append(("const_cycle", [2] * n, [(i, (i + 1) % n, I2, I2) for i in range(n)]))                       # H0 = 2
cases.append(("path_isolated", [1, 1, 1, 3, 0], [(0, 1, [[1.0]], [[1.0]]), (1, 2, [[2.0]], [[-1.0]])]))    # isolated dims 3 and 0
for eps, nm in ((1e-13, "near_singular_below_cut"), (1e-6, "near_singular_above_cut")):
    F = np.array([[1.0, 2.0], [2.0, 4.0 + eps]])                    # rank 1 up to eps: sigma ~ eps/5 vs cut 1e-9*max(1,sigma_max)
    cases.append((nm, [2, 2], [(0, 1, F.tolist(), (-F).tolist())]))
cases.append(("self_loop_cancel", [2], [(0, 0, [[1.0, 0.0], [0.0, 3.0]], [[1.0, 0.0], [0.0, 1.0]])]))      # explicit 0 after summing
cases.append(("zero_edge_dims", [2, 3], [(0, 1, [], [])]))
cases.append(("empty", [0, 0], []))

out_cases = []
for name, dv, edges in cases:
    nc = sum(dv); x = rng.normal(size=nc)
    E = [(u, v, np.asarray(Fu, float).reshape(len(Fu), dv[u]) if len(Fu) else np.zeros((0, dv[u])),
          np.asarray(Fv, float).reshape(len(Fv), dv[v]) if len(Fv) else np.zeros((0, dv[v]))) for u, v, Fu, Fv in edges]
    D, L, r, G = ref(dv, E, x)
    out_cases.append(dict(name=name, dv=dv, edges=[dict(u=u, v=v, de=int(Fu.shape[0]), fu=Fu.reshape(-1).tolist(), fv=Fv.reshape(-1).tolist()) for u, v, Fu, Fv in E],
                          x=x.tolist(), delta=D.toarray().tolist(), L=L.toarray().tolist(), proj=(G @ G.T).tolist(), **r))

def apply(doc, op):
    def ptr(p):
        cur = doc
        for k in p.split("/")[1:]: cur = cur[int(k)] if isinstance(cur, list) else cur[k]
        return cur
    par, key = op["path"].rsplit("/", 1); tgt = ptr(par) if par else doc
    if op["op"] == "set":
        if isinstance(tgt, list): tgt[int(key)] = copy.deepcopy(op["value"])
        else: tgt[key] = copy.deepcopy(op["value"])
    elif op["op"] == "remove": tgt.pop(int(key) if isinstance(tgt, list) else key)
    elif op["op"] == "append": ptr(op["path"]).append(copy.deepcopy(op["value"]))

S = json.loads((ROOT / "fixtures/m1/spec.json").read_text()); MUT = json.loads((ROOT / "fixtures/m1/mutations.json").read_text())
dv, E, x = P.from_spec(S); D, L, r, G = ref(dv, E, x); C = D.tocoo()
spec = dict(rows=C.row.tolist(), cols=C.col.tolist(), vals=C.data.tolist(), shape=list(D.shape), **r)
table = []
for c in MUT:
    T = copy.deepcopy(S)
    for op in c["ops"]: apply(T, op)
    try: models.validate_spec(T); got = "ok"
    except Exception as e: got = "reject"
    row = dict(name=c["name"], models=got, rust=c["rust"])
    if got == "ok":
        dv, E, x = P.from_spec(T); row.update(ref(dv, E, x)[2])
    table.append(row)
json.dump(dict(cases=out_cases, spec=spec, table=table), open(ROOT / "fixtures/m1/rust_parity.json", "w"))
print(f"cases {len(out_cases)}  spec delta {spec['shape']} nnz {len(spec['vals'])}  table {len(table)} rows ({sum(t['models'] == 'ok' for t in table)} accepted)")
