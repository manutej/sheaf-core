"""sheaf-core: deterministic SYNTHETIC m1 fixtures (reference-only Python).
The gated frontier suite runs on a frozen spec.json exported from a private Notion research graph; that export is NOT
shipped here. This script builds a structurally equivalent synthetic spec (same schema, same node kinds / edge types /
ordering conventions, random hex ids, fake titles/urls, seeded x and maps), then regenerates the shared mutation table
with the unchanged generator and the Python reference discord used by sheaf-contract's parity tests.
usage: python python/m1/make_fixtures.py   (needs numpy, scipy)"""
import json, pathlib, subprocess, sys
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; ROOT = HERE.parents[1]; OUT = ROOT / "fixtures/m1"; OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261019); d = 10
hexid = lambda: "".join(rng.choice(list("0123456789abcdef"), 32))
STAT = ["verified", "partly verified", "search-verified", "claimed", "unverified", "refuted", "not found"]
nodes = []
S, T, L, CPS = 24, 4, 2, 2
for i in range(S):
    nodes.append(dict(id=hexid(), kind="source", url=f"https://example.org/paper/{i}", slug=f"src-{i}", title=f"Synthetic source {i}",
                      status=STAT[i % len(STAT)], fetched="None"))
for i in range(T): nodes.append(dict(id=hexid(), kind="topic", slug=f"topic-{i}", title=f"topic {i}"))
for i in range(L): nodes.append(dict(id=hexid(), kind="log", slug=f"Round {i}", title=f"Round {i} log"))
edges = []
for s in range(S):
    for c in range(CPS):
        p = nodes[s]; nodes.append(dict(id=p["id"] + f"c{c}", kind="claim", slug=f"{p['slug']}#c{c}", title=f"synthetic claim {c} of {p['slug']}", parent=p["id"]))
        edges.append(dict(u=s, v=len(nodes) - 1, type="claim_of"))
for s in range(S):
    for t in rng.choice(T, size=2, replace=False): edges.append(dict(u=s, v=S + int(t), type="about"))
    if s % 2 == 0: edges.append(dict(u=s, v=S + T + (s // 2) % L, type="logged"))
for a in range(0, S - 3, 3): edges.append(dict(u=a, v=a + 2, type="related"))
I = np.eye(d); Fv = I.copy(); Fv[d - 1, d - 1] = 0.9   # like the frozen spec: neither side the identity
maps = {}
for t in ["about", "related", "logged", "claim_of"]: maps[t] = [(I + 0.3 * rng.standard_normal((d, d))).round(10).tolist(), Fv.tolist()]
Fu = I.copy(); Fu[d - 1, d - 1] = 0.5; maps["same_slug"] = [Fu.tolist(), Fv.tolist()]
x = rng.standard_normal((len(nodes), d)).round(10).tolist()
spec = dict(d=d, nodes=nodes, x=x, edges=edges, maps=maps)
(OUT / "spec.json").write_text(json.dumps(spec, indent=0) + "\n")
subprocess.run([sys.executable, str(HERE / "gen_mutations.py")], check=True)
sys.path.insert(0, str(ROOT / "python")); import sheafcore as sc
dv, E, X = sc.from_spec(spec); D = sc.coboundary(dv, E); de = [d] * len(E)
(OUT / "py_discord.json").write_text(json.dumps(sc.edge_discord(D, X, de).tolist()) + "\n")
r = np.random.default_rng(7); out = []
for k in range(6):   # same recipe as frontier conform.sh, FIXED seed here (the gate draws a fresh one)
    Xk = r.standard_normal((len(nodes), d)) * (1.0 if k < 4 else 3.0)
    if k == 5: Xk[:, -2:] = 0.0
    out.append({"x": Xk.tolist(), "discord": sc.edge_discord(D, Xk.reshape(-1), de).tolist()})
(OUT / "py_discord_rand.json").write_text(json.dumps(out) + "\n")
print("fixtures written:", sorted(p.name for p in OUT.iterdir()))
