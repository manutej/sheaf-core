#!/usr/bin/env python3
"""Tool-first experiments: sheaf-core's `sheaf` binary used as a black-box tool call.

Nothing under crates/ is touched. Every case states its ground truth BEFORE the run:
either a hand-derived value (noted in the case) or an independent numpy oracle that
re-derives delta0 / Hodge split / harmonic extension from the spec alone.

Usage: python3 run_experiments.py [path/to/sheaf]   -> writes out/<case>.{spec,report}.json + out/results.json
Exit code 0 iff every check passes.
"""
import json, math, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "target", "release", "sheaf")

# ---------------------------------------------------------------- independent oracle (numpy only)
def oracle_d0(spec):
    ids = [n["id"] for n in spec["nodes"]]; dims = [n["dim"] for n in spec["nodes"]]
    on = np.concatenate([[0], np.cumsum(dims)]).astype(int)
    rows = []
    for e in spec["edges"]:
        u, v = ids.index(e["u"]), ids.index(e["v"])
        fu = np.array(e.get("Fu", np.eye(dims[u])), float); fv = np.array(e.get("Fv", np.eye(dims[v])), float)
        blk = np.zeros((fu.shape[0], on[-1])); blk[:, on[u]:on[u+1]] -= fu; blk[:, on[v]:on[v+1]] += fv
        rows.append(blk)
    return np.vstack(rows), on

def oracle_h0(d0, tol=1e-9):
    s = np.linalg.svd(d0, compute_uv=False) if d0.size else np.array([])
    return d0.shape[1] - int((s > tol * max(1.0, s.max() if s.size else 1.0)).sum())

def oracle_hodge(d0, w):
    pot, *_ = np.linalg.lstsq(d0, w, rcond=None)
    ex = d0 @ pot; ha = w - ex
    t = float(w @ w)
    return float(ex @ ex) / t, float(ha @ ha) / t, ha

def oracle_extension(d0, on, bnd):
    n = d0.shape[1]; fixed = np.zeros(n, bool); xb = np.zeros(n)
    for i, val in bnd.items():
        fixed[on[i]:on[i+1]] = True; xb[on[i]:on[i+1]] = val
    A = d0[:, ~fixed]; rhs = -d0[:, fixed] @ xb[fixed]
    xi, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    x = xb.copy(); x[~fixed] = xi
    r = d0 @ x
    unique = (A.shape[1] - np.linalg.matrix_rank(A)) == 0
    return x, float(r @ r), bool(unique)

# ---------------------------------------------------------------- harness
def run(name, spec, boundary=None, extra=()):
    os.makedirs(OUT, exist_ok=True)
    sp = os.path.join(OUT, f"{name}.spec.json"); rp = os.path.join(OUT, f"{name}.report.json")
    json.dump(spec, open(sp, "w"), indent=1)
    cmd = [BIN, "analyze", sp, "--out", rp, *extra]
    if boundary is not None:
        bp = os.path.join(OUT, f"{name}.boundary.json"); json.dump(boundary, open(bp, "w"), indent=1); cmd += ["--boundary", bp]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{name}: sheaf exited {p.returncode}: {p.stderr}")
    return json.load(open(rp))

RESULTS = []
def check(case, what, ok, observed, expected):
    RESULTS.append({"case": case, "check": what, "pass": bool(ok), "observed": observed, "expected": expected})

def close(a, b, rtol=1e-8, atol=1e-10):
    return a is not None and abs(a - b) <= atol + rtol * abs(b)

def ident(n): return np.eye(n).tolist()
def rot(deg):
    t = math.radians(deg); return [[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]]

# ================================================================ E1 memory / KB integrity
def e1():
    # E1a: one corrupted claim. Constant R^1 sheaf; status +1 = "true". c4 wrongly -1.
    # Hand truth: only c4's edges carry discord (each 2^2 = 4); top incident-discord node is c4; H0 = 1; w := dx -> exact share 1.
    nodes = ["s1", "s2", "c1", "c2", "c3", "c4", "c5", "t1"]
    edges = [("s1","c1"),("s1","c2"),("s1","c3"),("s2","c3"),("s2","c4"),("s2","c5"),("c1","t1"),("c4","t1"),("c5","t1"),("c2","c3")]
    x = {n: [1.0] for n in nodes}; x["c4"] = [-1.0]
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": n, "dim": 1} for n in nodes],
            "edges": [{"u": u, "v": v, "label": f"{u}->{v}"} for u, v in edges], "x": [x[n] for n in nodes], "options": {"null_reps": 0}}
    r = run("E1a_corrupted_claim", spec)
    top = r["nodes_by_incident_discord"][0]["node"]
    check("E1a", "top incident-discord node is the corrupted claim", top == "c4", top, "c4")
    check("E1a", "h0 = 1 (one connected KB)", r["h0"]["dim"] == 1, r["h0"]["dim"], 1)
    check("E1a", "exact share = 1 (single bad node is a fixable, exact error)", close(r["hodge"]["exact_share"], 1.0, 1e-9), r["hodge"]["exact_share"], 1.0)

    # E1b: irreducible contradiction. Triangle A,B,C with "contradicts" maps (Fu=-1: B should be -A).
    # Odd cycle of sign flips -> holonomy -1 -> H0 = 0. Pin A = 1: hand solution B = C = -1/3, energy 4/3, unique.
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": n, "dim": 1} for n in "ABC"],
            "edges": [{"u": "A", "v": "B", "Fu": [[-1]], "Fv": [[1]], "label": "A contradicts B"},
                      {"u": "B", "v": "C", "Fu": [[-1]], "Fv": [[1]], "label": "B contradicts C"},
                      {"u": "C", "v": "A", "Fu": [[-1]], "Fv": [[1]], "label": "C contradicts A"}],
            "x": [[1.0], [-1.0], [1.0]], "options": {"null_reps": 0}}
    r = run("E1b_odd_contradiction_cycle", spec, {"boundary": [{"node": "A", "value": [1.0]}]})
    ext = r["extension"]; vals = {v["id"]: v["value"] for v in ext["values"]}
    check("E1b", "h0 = 0 (no globally consistent truth assignment exists)", r["h0"]["dim"] == 0, r["h0"]["dim"], 0)
    check("E1b", "extension energy = 4/3 (irreducible contradiction)", close(ext["energy"], 4/3), ext["energy"], 4/3)
    check("E1b", "B = C = -1/3", close(vals["B"][0], -1/3) and close(vals["C"][0], -1/3), [vals["B"], vals["C"]], [-1/3, -1/3])
    check("E1b", "unique", ext["unique"] is True, ext["unique"], True)

    # E1c: fill unknown statuses. Path s -> c1 -> c2 -> t, pin s=1, t=1 plus unknown orphan claim o (isolated).
    # Hand truth: c1 = c2 = 1, energy 0; orphan is unconstrained -> unique = false (h0_relative = 1).
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": n, "dim": 1} for n in ["s","c1","c2","t","o"]],
            "edges": [{"u": "s", "v": "c1"}, {"u": "c1", "v": "c2"}, {"u": "c2", "v": "t"}],
            "x": [[1.0],[0.0],[0.0],[1.0],[0.0]], "options": {"null_reps": 0}}
    r = run("E1c_fill_unknown", spec, {"boundary": [{"node": "s", "value": [1.0]}, {"node": "t", "value": [1.0]}]})
    ext = r["extension"]; vals = {v["id"]: v["value"] for v in ext["values"]}
    check("E1c", "filled c1 = c2 = 1", close(vals["c1"][0], 1.0) and close(vals["c2"][0], 1.0), [vals["c1"], vals["c2"]], [1, 1])
    check("E1c", "orphan claim flagged non-unique (h0_relative = 1)", ext["unique"] is False and ext["h0_relative"] == 1, [ext["unique"], ext["h0_relative"]], [False, 1])

# ================================================================ E2 multi-system reconciliation
def e2():
    # E2a: triangular FX arbitrage. Stalk = log price; edge observation w_e = log rate. Arbitrage-free => w exact.
    # Two triangles sharing USD-EUR. Inject eps on EUR->MXN. Oracle: numpy Hodge split.
    nodes = ["USD", "EUR", "MXN", "GBP"]; p = {"USD": 0.0, "EUR": 0.08, "MXN": -2.9, "GBP": 0.25}
    edges = [("USD","EUR"),("EUR","MXN"),("MXN","USD"),("USD","GBP"),("GBP","EUR")]
    eps = 0.03
    w = [[p[v] - p[u] + (eps if (u, v) == ("EUR","MXN") else 0.0)] for u, v in edges]
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": n, "dim": 1} for n in nodes],
            "edges": [{"u": u, "v": v, "label": f"{u}->{v}"} for u, v in edges], "w": w, "options": {"null_reps": 50, "seed": 7}}
    r = run("E2a_fx_arbitrage", spec)
    d0, _ = oracle_d0(spec); ex, ha, hv = oracle_hodge(d0, np.array(w).ravel())
    check("E2a", "harmonic share matches numpy oracle", close(r["hodge"]["harmonic_share"], ha, 1e-7), r["hodge"]["harmonic_share"], ha)
    tri = {"USD->EUR","EUR->MXN","MXN->USD"}
    top3 = {e["label"] for e in r["edges_by_harmonic_energy"][:3]}
    check("E2a", "top-3 harmonic edges = the arbitrage triangle", top3 == tri, sorted(top3), sorted(tri))
    # control: no arbitrage -> harmonic share ~ 0
    spec0 = dict(spec); spec0["w"] = [[p[v] - p[u]] for u, v in edges]
    r0 = run("E2a_fx_control", spec0)
    check("E2a", "control (arbitrage-free) harmonic share ~ 0", r0["hodge"]["harmonic_share"] < 1e-12, r0["hodge"]["harmonic_share"], 0)

    # E2b: one system wrong, non-identity (currency) maps. Revenue for one account in CRM(USD), ERP(EUR), Sheet(MXN), BI(USD).
    rate = {"USD": 1.0, "EUR": 0.92, "MXN": 17.0}
    sysc = {"CRM": "USD", "ERP": "EUR", "SHEET": "MXN", "BI": "USD"}
    truth = 1000.0
    x = {s: truth * rate[c] for s, c in sysc.items()}; x["ERP"] = 1000.0  # ERP booked USD amount as EUR (unit slip)
    pairs = [("CRM","ERP"),("CRM","SHEET"),("CRM","BI"),("ERP","SHEET"),("ERP","BI"),("SHEET","BI")]
    # restriction to a common USD edge stalk: (dx)_e = x_v / rate_v - x_u / rate_u
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": s, "dim": 1} for s in sysc],
            "edges": [{"u": u, "v": v, "Fu": [[1/rate[sysc[u]]]], "Fv": [[1/rate[sysc[v]]]], "label": f"{u}~{v}"} for u, v in pairs],
            "x": [[x[s]] for s in sysc], "options": {"null_reps": 0}}
    r = run("E2b_unit_slip", spec)
    top = r["nodes_by_incident_discord"][0]["node"]
    check("E2b", "wrong system ranked first by incident discord", top == "ERP", top, "ERP")
    check("E2b", "h0 = 1 (rates consistent: one true amount up to scale)", r["h0"]["dim"] == 1, r["h0"]["dim"], 1)

# ================================================================ E3 multi-agent consensus
def e3():
    # E3a: 5 agents, belief over 3 hypotheses, complete graph, identity frames; agent a4 is the outlier.
    ag = ["a0","a1","a2","a3","a4"]; b = {a: [0.7, 0.2, 0.1] for a in ag}; b["a4"] = [0.1, 0.2, 0.7]
    b["a1"] = [0.68, 0.22, 0.10]; b["a2"] = [0.72, 0.18, 0.10]
    edges = [(ag[i], ag[j]) for i in range(5) for j in range(i+1, 5)]
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": a, "dim": 3} for a in ag],
            "edges": [{"u": u, "v": v} for u, v in edges], "x": [b[a] for a in ag], "options": {"null_reps": 0}}
    r = run("E3a_outlier_agent", spec)
    top = r["nodes_by_incident_discord"][0]["node"]
    check("E3a", "outlier agent ranked first", top == "a4", top, "a4")
    check("E3a", "h0 = 3 (consensus space = shared belief vector)", r["h0"]["dim"] == 3, r["h0"]["dim"], 3)

    # E3b/c: frame holonomy. 3 agents, 2-d frames, each edge rotates by theta. Hand: theta=120 -> R^3 = I -> H0 = 2;
    # theta=90 -> R^3 = R(270) has no fixed vector -> H0 = 0 (no consensus possible in any frame). H1 = H0 (rank-nullity, n0 = n1 = 6).
    for tag, th, h in (("E3b_frames_consistent", 120, 2), ("E3c_frames_twisted", 90, 0)):
        spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": a, "dim": 2} for a in ["p","q","r"]],
                "edges": [{"u": u, "v": v, "Fu": rot(th), "Fv": ident(2)} for u, v in [("p","q"),("q","r"),("r","p")]],
                "x": [[1,0],[0,1],[1,1]], "options": {"null_reps": 0}}
        rr = run(tag, spec)
        check(tag[:3], f"theta={th}: h0 = {h}", rr["h0"]["dim"] == h, rr["h0"]["dim"], h)
        check(tag[:3], f"theta={th}: h1 = {h}", rr["h1"]["dim"] == h, rr["h1"]["dim"], h)

# ================================================================ E4 studio: harmonic extension + oracle sweep
def e4():
    # E4a: path of 5, pin ends 0 and 4 -> interior 1,2,3 (hand), energy 4, unique.
    spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": f"n{i}", "dim": 1} for i in range(5)],
            "edges": [{"u": f"n{i}", "v": f"n{i+1}"} for i in range(4)], "x": [[0]]*5, "options": {"null_reps": 0}}
    r = run("E4a_path_extension", spec, {"boundary": [{"node": "n0", "value": [0.0]}, {"node": "n4", "value": [4.0]}]})
    vals = {v["id"]: v["value"][0] for v in r["extension"]["values"]}
    check("E4a", "interior = 1,2,3", all(close(vals[f"n{i}"], i) for i in (1,2,3)), [vals[f"n{i}"] for i in (1,2,3)], [1,2,3])
    check("E4a", "energy = 4", close(r["extension"]["energy"], 4.0), r["extension"]["energy"], 4.0)

    # E4b: randomized oracle sweep (the 'studio' upload path): random graph sheaves with random maps vs numpy.
    rng = np.random.default_rng(20261102); worst = {"h0": 0, "harm": 0.0, "ext": 0.0}; n_cases = 40
    for t in range(n_cases):
        n = int(rng.integers(4, 12)); m = int(rng.integers(n - 1, 2 * n + 2)); d = int(rng.integers(1, 4))
        E = [(i, i + 1) for i in range(n - 1)] + [tuple(sorted(rng.choice(n, 2, replace=False).tolist())) for _ in range(m - n + 1)]
        spec = {"format": "frontier-sheaf/v1", "nodes": [{"id": f"v{i}", "dim": d} for i in range(n)],
                "edges": [{"u": f"v{u}", "v": f"v{v}", "Fu": rng.normal(size=(d, d)).round(3).tolist(), "Fv": rng.normal(size=(d, d)).round(3).tolist()} for u, v in E],
                "w": rng.normal(size=(len(E), d)).round(4).tolist(), "options": {"null_reps": 0}}
        bnd = {0: rng.normal(size=d).round(3).tolist(), n - 1: rng.normal(size=d).round(3).tolist()}
        r = run(f"E4b_sweep_{t:02d}", spec, {"boundary": [{"node": f"v{i}", "value": v} for i, v in bnd.items()]})
        d0, on = oracle_d0(spec)
        worst["h0"] = max(worst["h0"], abs(r["h0"]["dim"] - oracle_h0(d0)))
        _, ha, _ = oracle_hodge(d0, np.array(spec["w"]).ravel())
        worst["harm"] = max(worst["harm"], abs(r["hodge"]["harmonic_share"] - ha))
        _, en, _ = oracle_extension(d0, on, bnd)
        worst["ext"] = max(worst["ext"], abs(r["extension"]["energy"] - en) / max(1.0, en))
    check("E4b", f"{n_cases} random sheaves: H0 equals oracle", worst["h0"] == 0, worst["h0"], 0)
    check("E4b", f"{n_cases} random sheaves: harmonic share within 1e-7", worst["harm"] < 1e-7, worst["harm"], "<1e-7")
    check("E4b", f"{n_cases} random sheaves: extension energy within 1e-7 rel", worst["ext"] < 1e-7, worst["ext"], "<1e-7")

if __name__ == "__main__":
    for f in (e1, e2, e3, e4):
        try: f()
        except Exception as ex: RESULTS.append({"case": f.__name__, "check": "ran", "pass": False, "observed": repr(ex), "expected": "no error"})
    json.dump(RESULTS, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=str)
    for r in RESULTS: print(("PASS " if r["pass"] else "FAIL ") + f"{r['case']:<4} {r['check']}  obs={r['observed']}")
    npass = sum(r["pass"] for r in RESULTS); print(f"\n{npass}/{len(RESULTS)} checks pass")
    sys.exit(0 if npass == len(RESULTS) else 1)
