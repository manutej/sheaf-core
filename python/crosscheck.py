"""Independent numpy cross-check of the Rust spike.
(1) Reads out/fixtures.json, rebuilds each coboundary from triplets, recomputes ranks with numpy SVD,
    checks d^2=0 and Betti numbers vs Rust.
(2) Builds torus7 / RP2 / sphere boundary matrices from scratch in numpy (no Rust input) and checks Betti.
"""
import json, itertools, numpy as np
def rank(A):
    if A.size == 0: return 0
    return int(np.linalg.matrix_rank(A))
fx = json.load(open("out/fixtures.json"))
ok = True
for f in fx:
    D = []
    for d in f["delta"]:
        M = np.zeros((d["rows"], d["cols"]))
        for i, j, v in d["t"]: M[i, j] += v
        D.append(M)
    for k in range(len(D) - 1):
        r = np.linalg.norm(D[k+1] @ D[k]); assert r < 1e-10, (f["name"], r)
    dims = f["dims"]; rk = [rank(M) for M in D]
    b = [dims[k] - (rk[k-1] if k else 0) - (rk[k] if k < len(rk) else 0) for k in range(len(dims))]
    match = b == f["betti_rust"]; ok &= match
    print(f"{f['name']:24s} dims={dims} numpy={b} rust={f['betti_rust']} {'OK' if match else 'MISMATCH'}")
def betti_from_faces(n, faces):
    faces = [tuple(sorted(t)) for t in faces]
    edges = sorted({e for t in faces for e in itertools.combinations(t, 2)})
    ei = {e: i for i, e in enumerate(edges)}
    d0 = np.zeros((len(edges), n)); d1 = np.zeros((len(faces), len(edges)))
    for i, (a, b) in enumerate(edges): d0[i, a] -= 1; d0[i, b] += 1
    for j, (a, b, c) in enumerate(faces): d1[j, ei[(b, c)]] += 1; d1[j, ei[(a, c)]] -= 1; d1[j, ei[(a, b)]] += 1
    r0, r1 = rank(d0), rank(d1)
    return [n - r0, len(edges) - r0 - r1, len(faces) - r1]
torus = [f for i in range(7) for f in ([i, (i+1) % 7, (i+3) % 7], [i, (i+2) % 7, (i+3) % 7])]
rp2 = [[a-1, b-1, c-1] for a, b, c in [[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,6],[2,3,6],[2,4,5],[2,5,6],[3,4,5],[3,4,6]]]
for name, n, F, exp in [("torus7", 7, torus, [1,2,1]), ("rp2", 6, rp2, [1,0,0]), ("sphere", 4, list(itertools.combinations(range(4), 3)), [1,0,1])]:
    b = betti_from_faces(n, F); ok &= b == exp
    print(f"independent numpy {name:8s} {b} expected {exp} {'OK' if b == exp else 'MISMATCH'}")
print("ALL_OK" if ok else "FAIL")
