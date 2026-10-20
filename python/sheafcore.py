"""Reference sheaf core (Python), SPEC v1.18. Heterogeneous stalks: vertex dims dv[i] >= 0, edge dims de[j] >= 0.
Convention (Hansen & Ghrist 2019): (delta x)_e = F_v x_v - F_u x_u for oriented edge e=(u,v), F_* : R^{dv} -> R^{de}.
L = delta^T delta. H^0 = ker delta = ker L.
Unified H0 rank rule (R-H0-08): SVD route counts sigma > RANK_TOL*max(1, sigma_max) as nonzero;
eigen route counts lambda > EIG_TOL*max(1, lambda_max), lambda_max := max|lambda|. sigma_max / lambda_max are GLOBAL maxima (largest over all
components = the matrix's own sigma_max / lambda_max); the values themselves are computed per connected component
(sparse block structure, at most one component dense at a time, never a globally densified matrix). A component above
MAX_DENSE_DIM (5000), a dense block above MAX_DENSE_BYTES (2 GiB, upper bound incl. U/Vt/LAPACK workspace), or a
returned basis above MAX_DENSE_BYTES raises MemoryError('R-PERF-18: ...') before allocating (no in-contract fallback: the
LDL^T route is the unbuilt Julia SHOULD; the contract Rust crate implements validate() and discord only).
Fourth cap (v1.8): a sparse operator whose storage measure 8*(m+1) + 16*nnz + 8*(m+n) exceeds MAX_SPARSE_BYTES (4 GiB)
is refused with the same prefixed MemoryError before any allocation -- by the H0 routes on the given operator and (v1.9)
by coboundary before assembly (nnz bounded by sum_j de_j*(dv_u+dv_v)). Total edge dimension <= MAX_TOTAL_DIM (R-DIM-02).
v1.11: an offset index over k stalks or edges (8*(k+1) bytes) counts against the same cap and is refused from len() alone.
v1.12/v1.13: ONE admission (_op_guard) for every operator argument of laplacian / edge_discord / energy / global_sections /
h0_dim_from_L -- type (sparse or ndarray, else R-DIM-02), rank (ndarray.ndim == 2, else R-DIM-02) and the storage measure
with dense logical size m*n -- before any scan; energy admits L before x; laplacian also measures its OUTPUT
(nnz <= min(n*n, sum_i r_i**2)) before forming delta^T delta."""
import math
import numpy as np, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
RANK_TOL = 1e-9; EIG_TOL = 1e-10; MAX_DENSE_DIM = 5000; MAX_DENSE_BYTES = 2 * 2**30

MAX_TOTAL_DIM = 2**31 - 1
MAX_SPARSE_BYTES = 4 * 2**30
def _meta_guard(m, n, nnz):
    """R-PERF-18 sparse-storage cap (v1.8, SPEC R-PERF-18 clause 4): the H0 routes build CSR/CSC index arrays and
    component labels, upper bound 8*(m+1) + 16*nnz + 8*(m+n) bytes; above MAX_SPARSE_BYTES (4 GiB) they raise the prefixed
    MemoryError BEFORE anything is allocated (e.g. a (2**32, 3) operator with 3 nonzeros). Ordinary large sparse sheaves
    (m+n+nnz ~ 1e8) are far below it."""
    m, n, nnz = int(m), int(n), int(nnz)
    if 8 * (m + 1) + 16 * nnz + 8 * (m + n) > MAX_SPARSE_BYTES: raise MemoryError(f"R-PERF-18: sparse operator {m}x{n} (nnz {nnz}) storage exceeds MAX_SPARSE_BYTES")  # M:sparse_storage_cap
def _op_guard(A):
    """v1.12 (R11 security) / v1.13 (R12 security): the SAME admission measure for EVERY operator argument of laplacian /
    edge_discord / energy / global_sections / h0_dim_from_L -- sparse: (m, n, nnz); dense ndarray: (m, n, m*n), i.e. every
    logical entry counts (a zero-stride view of shape (2**31, 1) is 8 bytes real but 2**31 logical entries). v1.13: an
    ndarray of any rank other than 2 (rank-1 / rank-3 zero-stride views included) and any non-matrix object are refused
    with R-DIM-02 from .ndim / type alone, so no isfinite scan, conversion, transpose or matvec ever touches them."""
    if sp.issparse(A): _meta_guard(*A.shape, A.nnz)
    elif isinstance(A, np.ndarray):
        if A.ndim != 2: raise ValueError(f"R-DIM-02: operator must be 2-D, got ndim {A.ndim}")  # M:op_rank
        _meta_guard(A.shape[0], A.shape[1], A.size)
    else: raise ValueError("R-DIM-02: operator must be a scipy.sparse matrix or a 2-D numpy array")   # untagged on purpose: every public non-matrix input dies in the type/shape access of the next statement anyway, so a mutant is behaviourally equivalent except for the message; pinned by the direct-core call energy([[1.0]], ...) in test_R_PERF_18_core_x_checked_before_conversion
def _gram_guard(delta):
    """v1.13 (R12 security): laplacian's OUTPUT delta^T delta is measured too, before the product is formed. Its nnz is at most
    min(n*n, sum_i r_i**2) with r_i the stored entries of row i (each row contributes at most an r_i x r_i block), so a
    (1, n) row -- dense or sparse -- whose input passes the cap but whose n x n Gram would not, is refused (R-PERF-18).
    delta is CSR here (laplacian converts first); r = diff(indptr) counts stored entries (explicit zeros and duplicates
    included, so it can only over-estimate): O(m), within the already-admitted input measure."""
    n = delta.shape[1]; r = np.diff(delta.indptr).astype(np.int64)
    _meta_guard(n, n, min(int(n) * int(n), int(np.dot(r, r))))
def _index_guard(k, what):
    """v1.11 (R10 security): an offset index over k stalks/edges is an int64 vector of k+1 entries (8*(k+1) bytes), part of the
    operator's index storage; above MAX_SPARSE_BYTES it is refused (MemoryError R-PERF-18) from len() alone, BEFORE any
    list()/tolist()/offsets()/bincount materialisation (a zero-stride dims vector of logical length 2**30 costs 8 bytes real)."""
    if 8 * (int(k) + 1) > MAX_SPARSE_BYTES: raise MemoryError(f"R-PERF-18: {what} offset index of {k} entries exceeds MAX_SPARSE_BYTES")
def offsets(dims):
    o = np.zeros(len(dims) + 1, dtype=np.int64); o[1:] = np.cumsum(dims); return o

def _blocks(F, j, de, dn):
    """Strict R-DIM-02: map must be a 2-D array of exactly (de, dn), including dimension-0 cases."""
    F = np.asarray(F, float)
    if F.ndim != 2 or F.shape != (de, dn): raise ValueError(f"R-DIM-02: edge {j}: map shape {F.shape} != ({de},{dn})")  # M:blocks_exact_shape
    return F

def _is_int(a): return isinstance(a, (int, np.integer)) and not isinstance(a, (bool, np.bool_))

def _check_graph(dv, edges):
    """Input hygiene (R-DIM-02, and the endpoint-type/range part of R-ORIENT-01): stalk dims are non-negative integers,
    endpoints are integers (not bool, not float: the vectorised path would otherwise truncate 0.7 -> 0 silently) in [0, n).
    u < v is NOT required here: orientation is a spec-document rule (schema/models/Rust/validate.py); the core computes the
    same delta for either orientation of an edge (sign convention (dx)_e = Fv x_v - Fu x_u)."""
    if not all(_is_int(k) and k >= 0 for k in dv): raise ValueError("R-DIM-02: stalk dims must be non-negative integers")  # M:dims_nonneg_int
    if sum(int(k) for k in dv) > MAX_TOTAL_DIM: raise ValueError(f"R-DIM-02: total stalk dimension exceeds {MAX_TOTAL_DIM} (sparse int32 index range)")  # v1.7: no index wrap # M:total_dim_cap
    for j, e in enumerate(edges):
        if len(e) != 4: raise ValueError(f"R-DIM-02: edge {j} must be (u, v, Fu, Fv)")  # M:edge_arity
        if not (_is_int(e[0]) and _is_int(e[1])): raise ValueError(f"R-ORIENT-01: edge {j} endpoints must be integers")  # M:endpoints_int
        if not (0 <= e[0] < len(dv) and 0 <= e[1] < len(dv)): raise ValueError(f"R-ORIENT-01: edge {j} endpoint out of range")  # M:endpoints_range

def _map_shape_guard(edges):
    """v1.8 (R7 security): ndarray maps may LIE about their size (zero-stride views). Before anything is stacked, copied
    or scanned: total edge dimension <= MAX_TOTAL_DIM (R-DIM-02) and total map storage <= MAX_SPARSE_BYTES (R-PERF-18).
    v1.15 (R14 security): a scipy.sparse map is refused with R-DIM-02 HERE, before any np.asarray / np.stack (in the locked
    numpy 2.5.4 / scipy 1.18.1, np.asarray(csr, float) raises an unprefixed ValueError rather than densifying; the contract
    makes that an explicit, prefixed rejection so no future numpy/scipy __array__ change can densify a 20000^2 map)."""
    rows = elems = 0
    for e in edges:
        for F in e[2:4]:
            if isinstance(F, np.ndarray): elems += int(F.size); rows += int(F.shape[0]) if F.ndim else 0   # v1.17 (R16 security): Python-int accumulation, never an int64 that can wrap
            elif sp.issparse(F): raise ValueError("R-DIM-02: restriction maps must be dense 2-D arrays or nested lists, not scipy.sparse")  # M:sparse_map_rejected
    # untagged on purpose: edge-dimension / map-storage analogues of M:total_dim_cap / M:sparse_storage_cap for ndarray maps
    # that lie about their size; not reachable through the JSON call boundary, exercised by the direct core test
    # test_R_PERF_18_core_lying_shapes_rejected_before_alloc (tests/test_spec_public.py)
    if rows > 2 * MAX_TOTAL_DIM: raise ValueError(f"R-DIM-02: total edge dimension exceeds {MAX_TOTAL_DIM}")
    if 8 * elems > MAX_SPARSE_BYTES: raise MemoryError("R-PERF-18: restriction maps exceed MAX_SPARSE_BYTES")
def coboundary(dv, edges):
    """edges: list of (u, v, Fu, Fv) with Fu.shape == (de, dv[u]), Fv.shape == (de, dv[v]). Raises ValueError on shape mismatch."""
    _index_guard(len(dv), "stalk"); _index_guard(len(edges), "edge")   # v1.11 (R10 security): before list()/offsets() materialise them; untagged on purpose: a 2**29-entry dims vector cannot cross the JSON runner, so the core audit could not kill it -- covered by the direct-core test test_R_PERF_18_core_offset_index_cap_before_materialisation
    dv = list(dv); _check_graph(dv, edges); _map_shape_guard(edges); ov = offsets(dv); m = len(edges)
    fast = _uniform(dv, edges)
    if fast is not None: return fast
    de = []
    for j, (_, _, Fu, _) in enumerate(edges):
        A = np.asarray(Fu, float)
        if A.ndim != 2: raise ValueError(f"R-DIM-02: edge {j}: map must be 2-D, got shape {A.shape}")  # M:generic_map_2d
        de.append(int(A.shape[0]))   # v1.17: Python ints, so sum(de) and the assembly measure cannot wrap
    if sum(de) > MAX_TOTAL_DIM: raise ValueError(f"R-DIM-02: total edge dimension exceeds {MAX_TOTAL_DIM}")   # v1.8 (lists too); untagged on purpose: the EDGE-dimension analogue of M:total_dim_cap (which caps stalk dims); reaching it needs > 2**31 list rows, so no public test can kill a mutant -- documented, not tagged
    _meta_guard(sum(de), sum(dv), sum(de[j] * (dv[e[0]] + dv[e[1]]) for j, e in enumerate(edges)))   # v1.9: assembly cap (tagged in _meta_guard)
    oe = offsets(de); R, C, V = [], [], []
    for j, (u, v, Fu, Fv) in enumerate(edges):
        for F, node, sign in ((_blocks(Fu, j, de[j], dv[u]), u, -1.0), (_blocks(Fv, j, de[j], dv[v]), v, 1.0)):
            if F.size == 0: continue
            if not np.isfinite(F).all(): raise ValueError(f"R-FINITE-01: edge {j}: non-finite map")  # M:generic_finite
            r, c = np.nonzero(F); R.append(oe[j] + r); C.append(ov[node] + c); V.append(sign * F[r, c])
    cat = lambda a, dt: np.concatenate(a).astype(dt) if a else np.zeros(0, dt)
    return sp.csr_matrix((cat(V, float), (cat(R, np.int64), cat(C, np.int64))), shape=(int(oe[-1]), int(ov[-1])))

def _uniform(dv, edges):
    """Vectorised path when every map is a 2-D array of one shape (de, d) and all incident stalks have dim d."""
    if not edges: return None
    try:
        Fu = np.stack([np.asarray(e[2], float) for e in edges]); Fv = np.stack([np.asarray(e[3], float) for e in edges])
    except ValueError: return None
    if Fu.ndim != 3 or Fu.shape != Fv.shape: return None
    m, de, d = Fu.shape; U = np.fromiter((e[0] for e in edges), np.int64, m); V = np.fromiter((e[1] for e in edges), np.int64, m)
    dva = np.asarray(dv, np.int64); n = len(dva)
    if (dva[U] != d).any() or (dva[V] != d).any(): return None   # mismatch -> generic path raises ValueError
    if not (np.isfinite(Fu).all() and np.isfinite(Fv).all()): raise ValueError("R-FINITE-01: non-finite map")  # M:fast_finite
    _meta_guard(m * de, int(dva.sum()), 2 * m * de * d)   # v1.9 (R8 security): before ANY O(d) / O(m*de*d) index allocation
    if de == 0 or d == 0: return sp.csr_matrix((m * de, int(dva.sum())))   # v1.10 (R9 performance): no index scratch for an all-zero-size block
    ov = offsets(dv); a = np.arange(de); b = np.arange(d)
    rows = (np.arange(m)[:, None, None] * de + a[None, :, None]) + np.zeros((1, 1, d), np.int64)
    R, C, Vv = [], [], []
    for F, node, sign in ((Fu, U, -1.0), (Fv, V, 1.0)):
        cols = ov[node][:, None, None] + b[None, None, :] + np.zeros((1, de, 1), np.int64); nz = F != 0
        R.append(rows[nz]); C.append(cols[nz]); Vv.append(sign * F[nz])
    return sp.csr_matrix((np.concatenate(Vv), (np.concatenate(R), np.concatenate(C))), shape=(m * de, int(ov[-1])))

def laplacian(delta):
    """L = delta^T delta (CSR). delta: scipy.sparse matrix or 2-D ndarray (v1.13: a dense delta is converted to CSR AFTER
    the input measure admits it; the Gram output is measured by _gram_guard before it is formed)."""
    _op_guard(delta)   # v1.10 (R9 security): before transpose/indptr allocation; v1.12: dense operators too
    delta = delta.tocsr() if sp.issparse(delta) else sp.csr_matrix(delta)   # storage O(nnz), admitted by _op_guard above
    _gram_guard(delta)   # v1.13 (R12 security); untagged on purpose: an unguarded mutant needs a >= 2**28-nnz Gram to observe, killed only by the direct-core test test_R_PERF_18_core_gram_output_capped (monkeypatched cap) and the builder test test_R_PERF_18_gram_output_capped
    return (delta.T @ delta).tocsr()

def _check_x(delta, x):
    """x: 1-d, length delta.shape[1], finite. v1.13 (R12 security): an ndarray x is shape-checked from .shape BEFORE any
    float conversion or isfinite scan, and every caller admits the operator first (_op_guard), so len(x) is bounded."""
    if isinstance(x, np.ndarray) and (x.ndim != 1 or x.shape[0] != delta.shape[1]): raise ValueError(f"R-DIM-02: x has shape {x.shape}, expected ({delta.shape[1]},)")   # untagged on purpose: the M:x_shape check below raises the SAME error, so a mutant differs only in what it allocates first (direct-core test test_R_PERF_18_core_x_checked_before_conversion)
    x = np.asarray(x, float)
    if x.ndim != 1 or x.shape[0] != delta.shape[1]: raise ValueError(f"R-DIM-02: x has shape {x.shape}, expected ({delta.shape[1]},)")  # M:x_shape
    if not np.isfinite(x).all(): raise ValueError("R-FINITE-01: x must be finite")  # M:x_finite
    return x

def edge_discord(delta, x, de):
    """Vectorised: a sparse matvec plus a bincount of squared residuals over the delta rows (O(rows + nnz))."""
    _op_guard(delta)   # v1.12 (R11 security): sparse AND dense, before the isfinite scan below; v1.13: the ONE admission flow (type, rank, measure) shared by all five operator entry points
    if sp.issparse(delta) and delta.format in ("dok", "lil"): delta = delta.tocsr()   # v1.18 (R17 security minor): DOK/LIL .data is not a numeric array; admitted (bounded) first
    del_data = delta.data if sp.issparse(delta) else np.asarray(delta)
    if not np.isfinite(del_data).all(): raise ValueError("R-FINITE-01: delta must be finite")      # v1.7: hand-built operators too # M:discord_delta_finite
    if isinstance(de, np.ndarray) and de.ndim != 1: raise ValueError("R-DIM-02: edge dims must be a 1-d sequence")   # v1.12 (R11 security): before tolist(), which expands EVERY axis
    _index_guard(len(de), "edge")   # v1.11 (R10 security): before tolist()/bincount(minlength=len(de)); untagged on purpose (same reason as in coboundary)
    de_l = de.tolist() if isinstance(de, np.ndarray) else list(de)                                   # v1.7: exact Python-int sum, no float truncation
    if not all(_is_int(k) and k >= 0 for k in de_l) or sum(int(k) for k in de_l) != delta.shape[0]: raise ValueError("R-DIM-02: edge dims must be non-negative integers summing to delta's row count")  # M:de_consistent
    de = np.asarray(de_l, dtype=np.int64)
    if de.ndim != 1 or (de < 0).any() or int(de.sum()) != delta.shape[0]: raise ValueError("R-DIM-02: edge dims must be >= 0 and sum to delta's row count")  # defence in depth for the SUM only (after int64 coercion it cannot reject 1.7 / True): NOT equivalent to M:de_consistent, whose exact-int check is the one that rejects non-integer dims
    r = delta @ _check_x(delta, x)
    out = np.bincount(np.repeat(np.arange(len(de)), de), weights=r * r, minlength=len(de)).astype(float)
    if not np.isfinite(out).all(): raise ValueError("R-FINITE-01: discord overflows float64 (finite inputs, non-finite residual)")   # v1.10 # M:discord_overflow
    return out

def energy(L, x):
    _op_guard(L)   # v1.12 (R11 security): before the isfinite scan of L (dense too); v1.13 (R12 security/performance): also before x is converted or scanned
    _check_x(L, x)  # M:energy_checks_x
    if sp.issparse(L) and L.format in ("dok", "lil"): L = L.tocsr()   # v1.18 (R17 security minor): as in edge_discord
    if not np.isfinite(L.data if sp.issparse(L) else np.asarray(L)).all(): raise ValueError("R-FINITE-01: L must be finite")   # v1.7 # M:energy_L_finite
    x = np.asarray(x, float); e = float(x @ (L @ x))
    if not math.isfinite(e): raise ValueError("R-FINITE-01: energy overflows float64 (finite inputs, non-finite result)")   # v1.10 # M:energy_overflow
    return e

def _components(A, n):
    """Column groups coupled through A (connected components of the bipartite row/column incidence graph).
    Uses only A's sparsity pattern (no Gram matrix, no fill-in)."""
    if n == 0: return []
    m = A.shape[0]; r, c = A.nonzero()                  # nonzero() drops explicitly stored zeros, so size the data from it
    P = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(m, n)) if m else sp.csr_matrix((0, n))
    B = sp.bmat([[None, P], [P.T, None]], format="csr") if m else sp.csr_matrix((n, n))
    k, lab = connected_components(B, directed=False); lab = lab[m:] if m else lab
    order = np.argsort(lab, kind="stable"); cuts = np.flatnonzero(np.diff(lab[order])) + 1
    return np.split(order, cuts)

def _svd_bytes(m, n):
    M = max(m, n); return 8 * (2 * M * n + 6 * n * n)        # padded block + U + Vt + LAPACK workspace (upper bound)

def _dense_block(delta, cols):
    sub = delta[:, cols]; rows = np.unique(sub.nonzero()[0])
    if _svd_bytes(len(rows), len(cols)) > MAX_DENSE_BYTES: raise MemoryError(f"R-PERF-18: dense SVD of {len(rows)}x{len(cols)} block exceeds MAX_DENSE_BYTES")  # defence in depth (fires only when this component is reached); pass-0 M:svd_byte_cap is the before-any-densification check
    if len(rows) == 0: return None
    D = sub[rows].toarray()
    if D.shape[0] < D.shape[1]: D = np.vstack([D, np.zeros((D.shape[1] - D.shape[0], D.shape[1]))])   # pad so thin SVD yields all n right vectors
    return D

def global_sections(delta, tol=RANK_TOL):
    """Orthonormal basis of ker delta: per-component SVD, rank rule sigma > tol*max(1, sigma_max) (sigma_max global).
    Pass 0 checks every component's limits before any densification. Pass 1: singular values only (one component
    dense at a time) -> global sigma_max and each component's rank -> exact output width; the output basis size is
    checked and allocated once. Pass 2: per component, full SVD again, write its null-space rows into the output."""
    _op_guard(delta)                   # v1.7: before ANY allocation (v1.8: dense input counts every entry)
    delta = sp.csr_matrix(delta); m, n = delta.shape
    if not np.isfinite(delta.data).all(): raise ValueError("R-FINITE-01: delta must be finite")  # M:delta_finite
    if n == 0: return np.zeros((0, 0))
    comps = _components(delta, n)
    for cols in comps:                                   # pass 0: validate limits for every component before any SVD
        if len(cols) > MAX_DENSE_DIM: raise MemoryError(f"R-PERF-18: component of dim {len(cols)} exceeds MAX_DENSE_DIM={MAX_DENSE_DIM}")  # M:svd_dim_cap
        rows = np.unique(delta[:, cols].nonzero()[0])
        if _svd_bytes(len(rows), len(cols)) > MAX_DENSE_BYTES: raise MemoryError(f"R-PERF-18: dense SVD of {len(rows)}x{len(cols)} block exceeds MAX_DENSE_BYTES")  # M:svd_byte_cap
    svs = []
    for cols in comps:                                   # pass 1: singular values only
        D = _dense_block(delta, cols)
        svs.append(np.linalg.svd(D, compute_uv=False) if D is not None else np.zeros(0)); del D
    smax = max((float(v.max()) for v in svs if v.size), default=0.0); thr = tol * max(1.0, smax)
    ranks = []
    for sv in svs: r = int((sv > thr).sum()); ranks.append(r)
    N = sum(len(c) - r for c, r in zip(comps, ranks))
    if 8 * n * N > MAX_DENSE_BYTES: raise MemoryError(f"R-PERF-18: dense basis {n}x{N} exceeds MAX_DENSE_BYTES")  # M:basis_bytes
    G = np.zeros((n, N)); k = 0
    for cols, r in zip(comps, ranks):                    # pass 2: one component dense at a time
        D = _dense_block(delta, cols); w = len(cols) - r
        if D is None:                                    # no incident rows: every coordinate of the component is a section
            G[np.asarray(cols), k + np.arange(len(cols))] = 1.0  # v1.8: scatter, no ambient n-vector per isolated column
        else:
            Vt = np.linalg.svd(D, full_matrices=False)[2]  # v1.7: U/S not kept alive into the next component
            G[cols, k:k + w] = Vt[r:].T; del Vt
        del D                                            # v1.8: component i's dense block is freed BEFORE component i+1 is densified
        k += w
    return G

def h0_dim_from_L(L, tol=EIG_TOL):
    """dim ker L via per-component eigenvalues; lambda counted nonzero iff lambda > tol*max(1, lambda_max), where
    lambda_max := max |lambda| over all components (= the largest eigenvalue for a PSD Laplacian, R-LAP-07)."""
    _op_guard(L)                                   # v1.7: before ANY allocation (v1.8: dense input counts every entry)
    L = sp.csr_matrix(L); n = L.shape[0]
    if L.shape[0] != L.shape[1]: raise ValueError(f"R-DIM-02: L must be square, got {L.shape}")  # M:L_square
    if not np.isfinite(L.data).all(): raise ValueError("R-FINITE-01: L must be finite")  # M:L_finite
    if n == 0: return 0
    r, c = L.nonzero()                                   # pattern only; explicitly stored zeros are not couplings
    k, lab = connected_components(sp.csr_matrix((np.ones(len(r)), (r, c)), shape=L.shape), directed=False); ws = []
    order = np.argsort(lab, kind="stable"); groups = np.split(order, np.flatnonzero(np.diff(lab[order])) + 1)   # O(n log n), not O(k*n)
    for cols in groups:                                  # all limits checked before any densification
        if len(cols) > MAX_DENSE_DIM: raise MemoryError(f"R-PERF-18: component of dim {len(cols)} exceeds MAX_DENSE_DIM={MAX_DENSE_DIM}")  # M:eig_dim_cap
        if 8 * 3 * len(cols) ** 2 > MAX_DENSE_BYTES: raise MemoryError(f"R-PERF-18: dense eig of dim {len(cols)} exceeds MAX_DENSE_BYTES")  # non-binding under the 5000-dim cap (SPEC R-PERF-18 note); kept as defence
    for cols in groups:
        ws.append(np.linalg.eigvalsh(L[cols][:, cols].toarray()))
    w = np.concatenate(ws); thr = tol * max(1.0, float(abs(w).max())); return int((w <= thr).sum())

def _is_num(v): return isinstance(v, (int, float)) and not isinstance(v, bool)
def _f64(v):
    """v1.9 (R8 security): a JSON number must be REPRESENTABLE as float64. An integer beyond float range (e.g. 10**309)
    is rejected as R-DIM-02 'not a valid number' -- the same verdict as models.validate_spec (pydantic) and validate.py;
    serde_json rejects it at parse time ('number out of range') -- instead of escaping as OverflowError or Inf."""
    try: float(v); return True
    except OverflowError: return False
def from_spec(spec):
    """spec.json -> (dv, edges, x). Intended for specs accepted by models.validate_spec (R-SCHEMA-01); v1.7 re-checks the
    typing it relies on instead of coercing (no strings/bools, x is exactly n rows of d numbers, maps are numeric)."""
    d = spec["d"]; n = len(spec["nodes"])
    if not (isinstance(d, int) and not isinstance(d, bool)): raise ValueError("R-DIM-02: d must be a JSON integer")
    X = spec["x"]
    if len(X) != n or any(not isinstance(r, list) or len(r) != d or not all(_is_num(v) for v in r) for r in X): raise ValueError("R-DIM-02: x must be n rows of d numbers")
    if not all(_f64(v) for r in X for v in r): raise ValueError("R-DIM-02: x entries must be float64-representable numbers")   # v1.9
    if not all(math.isfinite(v) for r in X for v in r): raise ValueError("R-FINITE-01: x must be finite")       # v1.8
    for k, pair in spec["maps"].items():
        if not (isinstance(pair, list) and len(pair) == 2 and all(isinstance(F, list) and all(isinstance(r, list) for r in F) for F in pair)): raise ValueError(f"R-DIM-02: map {k} must be a pair of matrices (lists of rows)")   # v1.11 (R10 security): no TypeError on scalar/non-list containers
        if any(not all(_is_num(v) for v in row) for F in pair for row in F): raise ValueError(f"R-DIM-02: map {k} has non-numeric entries")
        if any(len({len(row) for row in F}) > 1 for F in pair): raise ValueError(f"R-DIM-02: map {k} must be rectangular (ragged rows)")   # v1.18 (R17 security minor): prefixed, not NumPy's 'inhomogeneous shape'
        if not all(_f64(v) for F in pair for row in F for v in row): raise ValueError(f"R-DIM-02: map {k} entries must be float64-representable numbers")  # v1.9
        if not all(math.isfinite(v) for F in pair for row in F for v in row): raise ValueError(f"R-FINITE-01: map {k} must be finite")  # v1.8
    M = {k: (np.array(a, float), np.array(b, float)) for k, (a, b) in spec["maps"].items()}
    edges = [(e["u"], e["v"], *M[e["type"]]) for e in spec["edges"]]
    return [d] * n, edges, np.array(spec["x"], float).reshape(-1)
