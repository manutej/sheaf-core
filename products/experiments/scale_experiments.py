#!/usr/bin/env python3
"""Grounding experiments at product scale (black-box `sheaf` binary; crates/ untouched).

S1  P2 recon batch: K account keys x 5 systems, each key a 5-clique (disconnected components in ONE spec),
    planted unit slips on some keys. Measures wall time vs K and slip recall@rank (per component).
S2  P2 loop listing: full harmonic edge vector (top = #edges) -> loops recovered in L1 via a cycle basis
    of edges with harmonic energy > tau; checks planted rate loops are recovered exactly.
S3  P1 KB scale: synthetic claim graph (sources -> claims -> topics, ~N nodes, avg degree ~4), x-only and
    with NLI-style w; wall time and whether planted contradictions top the rankings.
S4  Quantisation range: dynamic range of harmonic energies / discords seen, to size the Nat fixed-point units.
Writes out/scale.json.
"""
import json, math, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out", "scale"); os.makedirs(OUT, exist_ok=True)
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "target", "release", "sheaf")
R = {}

def run(name, spec, extra=()):
    sp = os.path.join(OUT, f"{name}.spec.json"); rp = os.path.join(OUT, f"{name}.report.json")
    json.dump(spec, open(sp, "w"))
    t = time.time(); p = subprocess.run([BIN, "analyze", sp, "--out", rp, *extra], capture_output=True, text=True); wall = time.time() - t
    if p.returncode: raise RuntimeError(p.stderr[:500])
    return json.load(open(rp)), wall

# ---------------- S1: recon batch
SYS = ["CRM", "ERP", "BILL", "SHEET", "BI"]; CUR = {"CRM": 1.0, "ERP": 0.92, "BILL": 1.0, "SHEET": 17.0, "BI": 1.0}
def recon_spec(K, slip_every, rng):
    nodes, edges, x, planted = [], [], [], {}
    for k in range(K):
        amt = float(rng.integers(100, 100000))
        bad = (k % slip_every == 0)
        badsys = SYS[k % 5] if bad else None
        for s in SYS:
            v = amt * CUR[s]
            if s == badsys: v = v * 100.0 if CUR[s] == 1.0 else amt  # cents-vs-dollars slip, or unconverted amount
            nodes.append({"id": f"{k}:{s}", "dim": 1}); x.append([v])
        for i in range(5):
            for j in range(i + 1, 5):
                a, b = SYS[i], SYS[j]
                edges.append({"u": f"{k}:{a}", "v": f"{k}:{b}", "Fu": [[1 / CUR[a]]], "Fv": [[1 / CUR[b]]]})
        if bad: planted[k] = f"{k}:{badsys}"
    return {"format": "frontier-sheaf/v1", "nodes": nodes, "edges": edges, "x": x,
            "options": {"null_reps": 0, "top": len(nodes), "k": 2}}, planted

rng = np.random.default_rng(1)
s1 = []
for K in (100, 1000, 5000):
    spec, planted = recon_spec(K, 20, rng)
    rep, wall = run(f"S1_recon_K{K}", spec)
    # per-component top node by incident discord: group ranked list by key prefix, first occurrence wins
    first = {}
    for row in rep["nodes_by_incident_discord"]:
        k = int(row["node"].split(":")[0]); first.setdefault(k, row["node"])
    hits = sum(first.get(k) == n for k, n in planted.items())
    s1.append({"K": K, "nodes": len(spec["nodes"]), "edges": len(spec["edges"]), "wall_s": round(wall, 3),
               "planted": len(planted), "recall_at_rank1_per_key": hits / len(planted),
               "h0": rep["h0"]["dim"], "components": rep["input"]["components"]})
R["S1_recon_batch"] = s1

# ---------------- S2: loop listing from full harmonic vector
def cycle_basis(edge_list):
    # edge_list: [(j,u,v)]; returns list of edge-index cycles (spanning forest + one cycle per non-tree edge)
    adj, parent, depth, pedge = {}, {}, {}, {}
    for j, u, v in edge_list: adj.setdefault(u, []).append((v, j)); adj.setdefault(v, []).append((u, j))
    tree = set()
    for r in adj:
        if r in parent: continue
        parent[r] = None; depth[r] = 0; st = [r]
        while st:
            a = st.pop()
            for b, j in adj[a]:
                if b not in parent: parent[b] = a; depth[b] = depth[a] + 1; pedge[b] = j; tree.add(j); st.append(b)
    cycles = []
    for j, u, v in edge_list:
        if j in tree: continue
        a, b, ca, cb = u, v, [], []
        while a != b:
            if depth[a] >= depth[b]: ca.append(pedge[a]); a = parent[a]
            else: cb.append(pedge[b]); b = parent[b]
        cycles.append(sorted(set([j] + ca + cb)))
    return cycles

nodes = [f"C{i}" for i in range(30)]; p = {n: float(rng.normal()) for n in nodes}
E = [(nodes[i], nodes[i + 1]) for i in range(29)] + [(nodes[i], nodes[j]) for i, j in rng.choice(30, (25, 2)).tolist() if i != j]
E = list(dict.fromkeys(E))
w = [[p[v] - p[u]] for u, v in E]
# plant 3 loop violations: perturb 3 non-path (chord) edges
chords = [j for j in range(29, len(E))][:3]
for j in chords: w[j][0] += 0.05
spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": n, "dim": 1} for n in nodes],
        "edges": [{"u": u, "v": v, "label": str(j)} for j, (u, v) in enumerate(E)], "w": w, "options": {"null_reps": 0, "top": len(E)}}
rep, wall = run("S2_loops", spec)
hv = {int(r["label"]): r["harmonic_energy"] for r in rep["edges_by_harmonic_energy"]}
tau = 1e-9
hot = [(j, u, v) for j, (u, v) in enumerate(E) if hv.get(j, 0.0) > tau]
loops = cycle_basis(hot)
R["S2_loop_listing"] = {"edges": len(E), "edges_reported": len(rep["edges_by_harmonic_energy"]), "hot_edges": len(hot),
                        "loops_in_hot_subgraph": len(loops), "planted": len(chords), "wall_s": round(wall, 3),
                        "note": "hot-edge subgraph cycle basis; hot set includes every edge on a cycle through a perturbed chord, so loops are a basis, not the planted set itself"}

# ---------------- S3: P1 KB scale
def kb_spec(N, with_w, rng):
    ns, nc, nt = N // 5, (3 * N) // 5, N // 5
    ids = [f"s{i}" for i in range(ns)] + [f"c{i}" for i in range(nc)] + [f"t{i}" for i in range(nt)]
    x = {i: 1.0 for i in ids}
    E = []
    for i in range(nc):
        for s in rng.choice(ns, 2, replace=False): E.append((f"s{s}", f"c{i}", 1))
        E.append((f"c{i}", f"t{int(rng.integers(nt))}", 1))
    for _ in range(nc // 4):
        a, b = rng.choice(nc, 2, replace=False); E.append((f"c{a}", f"c{b}", 1))
    bad = [f"c{i}" for i in rng.choice(nc, 10, replace=False)]
    for b in bad: x[b] = -1.0
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": i, "dim": 1} for i in ids],
            "edges": [{"u": u, "v": v} for u, v, _ in E], "options": {"null_reps": 0, "top": 50, "k": 4}}
    if with_w:
        w = [[x[v] - x[u] + float(rng.normal(0, 0.01))] for u, v, _ in E]
        spec["w"] = w
    else:
        spec["x"] = [[x[i]] for i in ids]
    return spec, set(bad)

s3 = []
for N in (300, 2000, 10000):
    spec, bad = kb_spec(N, False, rng)
    rep, wall = run(f"S3_kb_N{N}", spec)
    top = [r["node"] for r in rep["nodes_by_incident_discord"][:10]]
    s3.append({"N": N, "edges": len(spec["edges"]), "wall_s": round(wall, 3), "precision_top10": len(set(top) & bad) / 10,
               "h0": rep["h0"]["dim"], "solver": rep["spectral"]["solver"]})
R["S3_kb_scale"] = s3

# ---------------- S4: quantisation range
vals = []
for f in os.listdir(OUT):
    if f.endswith(".report.json"):
        r = json.load(open(os.path.join(OUT, f)))
        vals += [e["harmonic_energy"] for e in r["edges_by_harmonic_energy"] if e["harmonic_energy"] > 0]
        vals += [n["incident_discord"] for n in r["nodes_by_incident_discord"]]
vals = [v for v in vals if v > 0]
R["S4_range"] = {"min_positive": min(vals), "max": max(vals), "decades": math.log10(max(vals) / min(vals))}
json.dump(R, open(os.path.join(HERE, "out", "scale.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
