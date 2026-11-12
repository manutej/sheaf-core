use proptest::prelude::*;
use sheaf_cell::{complex::CellComplex, num::{linalg::*, sparse::Csr, tol::H0Policy}, sheaf::Sheaf};
use sheaf_cohom::cohomology::*;
use sheaf_cohom::fp;
use sheaf_cohom::relative::h0_relative;
use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

/// random simplicial 2-complex on n vertices: random triangles + random extra edges
fn arb_complex() -> impl Strategy<Value = CellComplex> {
    (3usize..8).prop_flat_map(|n| {
        let tri = prop::collection::vec((0..n, 0..n, 0..n), 0..6);
        let ed = prop::collection::vec((0..n, 0..n), 0..8);
        (Just(n), tri, ed)
    }).prop_map(|(n, tri, ed)| {
        let faces: Vec<[usize; 3]> = tri.into_iter().filter(|&(a, b, c)| a != b && b != c && a != c).map(|(a, b, c)| [a, b, c]).collect();
        let edges: Vec<(usize, usize)> = ed.into_iter().filter(|&(a, b)| a != b).collect();
        CellComplex::simplicial(n, &edges, &faces)
    })
}

/// random sheaf that is guaranteed functorial: vertex/edge stalks R^d with random edge maps on graph
/// part, and face stalks = 0 or (if the complex has faces) a constant sheaf. Two families:
fn arb_sheaf() -> impl Strategy<Value = Sheaf> {
    (arb_complex(), 1usize..3, any::<u64>(), any::<bool>()).prop_map(|(cx, d, seed, twist)| {
        let has_faces = cx.n_cells.len() > 2;
        let mut s = Sheaf::constant(cx, d);
        if twist && !has_faces {
            // graph only: arbitrary per-incidence maps are always a valid sheaf (no codim-2 condition)
            let mut st = seed | 1;
            for m in s.maps[0].iter_mut() {
                *m = Mat::from_fn(d, d, |_, _| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 });
            }
        }
        s
    })
}

fn arb_cochain(n: usize, seed: u64) -> Vec<f64> {
    let mut st = seed | 1;
    (0..n).map(|_| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 20001) as f64) / 10000.0 - 1.0 }).collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn d_squared_is_zero(s in arb_sheaf()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        prop_assert!(cc.d2_residual <= D2_TOL);
    }

    // NOTE (Grok r1 #4): rank-nullity and Euler are identities given how dim is computed; the non-trivial
    // assertion here is `consistent()` (third, independent SVD of [d_k; d_{k-1}^T]). Golden values live in known.rs.
    #[test]
    fn rank_nullity_and_euler(s in arb_sheaf()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        let mut chi_c = 0i64; let mut chi_h = 0i64;
        for k in 0..cc.dims.len() {
            let h = cc.cohomology(k).unwrap();
            prop_assert!(h.consistent(), "dim {} vs basis {}", h.dim, h.harmonic_basis.ncols());
            // rank-nullity: dim ker d_k = dim C^k - rank d_k, and dim H^k = dim ker d_k - rank d_{k-1}
            prop_assert_eq!(h.dim + h.rank_in.rank + h.rank_out.rank, cc.dims[k]);
            let sgn = if k % 2 == 0 { 1 } else { -1 };
            chi_c += sgn * cc.dims[k] as i64; chi_h += sgn * h.dim as i64;
        }
        prop_assert_eq!(chi_c, chi_h);
    }

    #[test]
    fn hodge_components_orthogonal_and_sum(s in arb_sheaf(), seed in any::<u64>()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        for k in 0..cc.dims.len() {
            if cc.dims[k] == 0 { continue; }
            let w = arb_cochain(cc.dims[k], seed);
            let h = cc.hodge(k, &w);
            // ADR-015: delta_0 is cut by the H0 policy (sigma^2 <= 1e-10 * ||L0|_c||), so a dropped sigma (up to ~1e-5 *
            // sqrt(bound)) legitimately leaves (closed / coclosed) residuals of that size in the delta_0-facing checks;
            // every other identity must hold to 1e-9.
            let sd = if cc.delta.is_empty() { 0.0 } else { cc.d0_split().rank.sigma_dropped_max.unwrap_or(0.0) };
            let c = &h.cert;
            let d0_tol = 1e-9f64.max(1.01 * sd);
            let (closed_tol, coclosed_tol) = match k { 0 => (d0_tol, 1e-9), 1 => (1e-9, d0_tol), _ => (1e-9, 1e-9) };
            prop_assert!(c.potential_residual.max(c.ortho_exact_harm).max(c.ortho_exact_coex).max(c.ortho_harm_coex) < 1e-9
                && c.harmonic_closed <= closed_tol && c.harmonic_coclosed <= coclosed_tol, "k={} sigma_dropped={} cert={:?}", k, sd, h.cert);
            let (a, b, c) = h.energies();
            prop_assert!((a + b + c - dot(&w, &w)).abs() < 1e-9 * dot(&w, &w).max(1.0)); // Pythagoras
            // harmonic part lies in span of the harmonic basis
            let hb = cc.cohomology(k).unwrap().harmonic_basis;
            let p = project(&hb, &h.harmonic);
            prop_assert!(norm(&p.iter().zip(&h.harmonic).map(|(x, y)| x - y).collect::<Vec<_>>()) < 1e-8);
        }
    }
}

/// Gauge-twisted sheaf on an arbitrary 2-complex: F_{s<=t} = A_t A_s^{-1} with random invertible A_c.
/// Always functorial (d^2 = 0) and isomorphic to the constant sheaf, so its Betti numbers must match.
fn gauge(cx: &CellComplex, d: usize, seed: u64) -> Sheaf {
    let mut st = seed | 1;
    let mut rnd = move || { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 };
    let a: Vec<Vec<Mat<f64>>> = cx.n_cells.iter().map(|&n| (0..n).map(|_| {
        Mat::from_fn(d, d, |i, j| rnd() + if i == j { 3.0 } else { 0.0 })
    }).collect()).collect();
    let mut s = Sheaf::constant(cx.clone(), d);
    for k in 0..cx.inc.len() {
        for (j, inc) in cx.inc[k].iter().enumerate() {
            let inv = a[k][inc.facet].partial_piv_lu().inverse();
            s.maps[k][j] = &a[k + 1][inc.cofacet] * &inv;
        }
    }
    s
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]
    #[test]
    fn gauge_twist_on_2_complex_preserves_complex_and_betti(cx in arb_complex(), d in 1usize..3, seed in any::<u64>()) {
        let base = CochainComplex::from_sheaf(&Sheaf::constant(cx.clone(), d)).unwrap();
        let tw = CochainComplex::from_sheaf(&gauge(&cx, d, seed)).unwrap();
        prop_assert!(tw.d2_residual <= D2_TOL);
        for k in 0..base.dims.len() {
            prop_assert_eq!(base.cohomology(k).unwrap().dim, tw.cohomology(k).unwrap().dim);
            let w = arb_cochain(tw.dims[k], seed ^ 0x9e37);
            if tw.dims[k] > 0 { prop_assert!(tw.hodge(k, &w).cert.max() < 1e-8); }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]
    /// Grok r1 #5: generic random maps on a complex WITH faces violate functoriality and must be rejected.
    #[test]
    fn random_maps_on_faces_are_rejected(seed in any::<u64>()) {
        use sheaf_cell::fixtures::{triangle_filled, sphere};
        for cx in [triangle_filled(), sphere()] {
            let mut s = Sheaf::constant(cx, 2);
            let mut st = seed | 1;
            for m in s.maps[0].iter_mut() {
                *m = Mat::from_fn(2, 2, |_, _| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 });
            }
            let rejected = matches!(CochainComplex::from_sheaf(&s), Err(CohomError::NotAComplex { .. }));
            prop_assert!(rejected);
        }
    }
}

/// Field product and difference modulo a positive prime, reduced after every operation.
/// Independent of `fp::rank_mod_p` (no elimination, no Fermat inverse).
fn mod_mul(x: i64, y: i64, p: i64) -> i64 { ((x as i128 * y as i128).rem_euclid(p as i128)) as i64 }
fn mod_sub(x: i64, y: i64, p: i64) -> i64 { ((x as i128) - (y as i128)).rem_euclid(p as i128) as i64 }
fn mod_add(x: i64, y: i64, p: i64) -> i64 { ((x as i128) + (y as i128)).rem_euclid(p as i128) as i64 }

/// Largest k such that some k×k minor is nonzero in F_p. Entries are reduced with `rem_euclid` first.
fn det_rank(m: &[Vec<i64>], ncols: usize, p: i64) -> usize {
    debug_assert!(m.iter().all(|r| r.len() == ncols), "row length must equal ncols");
    if m.is_empty() || ncols == 0 { return 0; }
    let a: Vec<Vec<i64>> = m.iter().map(|r| r.iter().map(|x| x.rem_euclid(p)).collect()).collect();
    let (n, c) = (a.len(), ncols);
    for k in (1..=n.min(c)).rev() {
        if minor_nonzero(&a, k, p) { return k; }
    }
    0
}

fn minor_nonzero(a: &[Vec<i64>], k: usize, p: i64) -> bool {
    let (n, c) = (a.len(), a[0].len());
    match k {
        1 => a.iter().any(|row| row.iter().any(|&x| x != 0)),
        2 => {
            for i in 0..n { for i2 in (i + 1)..n {
                for j in 0..c { for j2 in (j + 1)..c {
                    let d = mod_sub(mod_mul(a[i][j], a[i2][j2], p), mod_mul(a[i][j2], a[i2][j], p), p);
                    if d != 0 { return true; }
                }}
            }}
            false
        }
        3 => {
            for i in 0..n { for i2 in (i + 1)..n { for i3 in (i2 + 1)..n {
                for j in 0..c { for j2 in (j + 1)..c { for j3 in (j2 + 1)..c {
                    if det3(a, [i, i2, i3], [j, j2, j3], p) != 0 { return true; }
                }}}
            }}}
            false
        }
        _ => false,
    }
}

/// det of the 3×3 submatrix on these rows and columns, via cofactor expansion in F_p.
fn det3(a: &[Vec<i64>], rs: [usize; 3], cs: [usize; 3], p: i64) -> i64 {
    let g = |i: usize, j: usize| a[rs[i]][cs[j]];
    let c0 = mod_sub(mod_mul(g(1, 1), g(2, 2), p), mod_mul(g(1, 2), g(2, 1), p), p);
    let c1 = mod_sub(mod_mul(g(1, 0), g(2, 2), p), mod_mul(g(1, 2), g(2, 0), p), p);
    let c2 = mod_sub(mod_mul(g(1, 0), g(2, 1), p), mod_mul(g(1, 1), g(2, 0), p), p);
    let pos = mod_add(mod_mul(g(0, 0), c0, p), mod_mul(g(0, 2), c2, p), p);
    mod_sub(pos, mod_mul(g(0, 1), c1, p), p)
}

fn for_each_matrix(rows: usize, cols: usize, p: i64, f: &mut impl FnMut(&[Vec<i64>])) {
    let mut m = vec![vec![0i64; cols]; rows];
    fn rec(m: &mut [Vec<i64>], p: i64, cell: usize, f: &mut impl FnMut(&[Vec<i64>])) {
        let cols = m.first().map_or(0, |r| r.len());
        let n = m.len() * cols;
        if cell == n { f(m); return; }
        let (i, j) = (cell / cols, cell % cols);
        for v in 0..p { m[i][j] = v; rec(m, p, cell + 1, f); }
    }
    rec(&mut m, p, 0, f);
}

/// Every matrix over a tiny prime field, plus anchors whose rank is fixed by hand
/// (not by the oracle): zero, identity, a scalar multiple of a row, and out-of-range residues.
#[test]
fn rank_mod_p_equals_determinantal_rank_on_small_fields() {
    assert_eq!(fp::rank_mod_p(&[vec![0, 0], vec![0, 0]], 2, 5), 0);
    assert_eq!(fp::rank_mod_p(&[vec![1, 0], vec![0, 1]], 2, 5), 2);
    assert_eq!(fp::rank_mod_p(&[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]], 3, 5), 3);
    // row2 = 2·row1 over F_5, so the matrix is singular. A Fermat inverse that is not a true
    // inverse reports rank 2 here (the pivot is not cleared).
    assert_eq!(fp::rank_mod_p(&[vec![2, 1], vec![4, 2]], 2, 5), 1);
    assert_eq!(det_rank(&[vec![2, 1], vec![4, 2]], 2, 5), 1);
    assert_eq!(fp::rank_mod_p(&[vec![-1, 0], vec![0, -1]], 2, 5), 2);
    assert_eq!(fp::rank_mod_p(&[vec![i64::MIN]], 1, 2), 0);
    assert_eq!(fp::rank_mod_p(&[vec![i64::MAX]], 1, 2), 1);
    for &p in &[2i64, 3, 5, 7, 11] {
        for rows in 0..=2 { for cols in 0..=2 {
            for_each_matrix(rows, cols, p, &mut |m| {
                assert_eq!(fp::rank_mod_p(m, cols, p), det_rank(m, cols, p), "p={p} m={m:?}");
            });
        }}
    }
    for &p in &[2i64, 3] {
        for rows in 0..=3 { for cols in 0..=3 {
            if rows <= 2 && cols <= 2 { continue; }
            for_each_matrix(rows, cols, p, &mut |m| {
                assert_eq!(fp::rank_mod_p(m, cols, p), det_rank(m, cols, p), "p={p} m={m:?}");
            });
        }}
    }
}

fn arb_prime() -> impl Strategy<Value = i64> {
    prop::sample::select(vec![
        2i64, 3, 5, 7, 11, 13, 17, 19, 31, 97, 251, 257, 65537,
        2_147_483_647, 4_294_967_291, 9_223_372_036_854_775_783,
    ])
}

fn arb_entry() -> impl Strategy<Value = i64> {
    prop_oneof![6 => -40i64..41, 1 => Just(i64::MIN), 1 => Just(i64::MAX), 1 => Just(-1i64)]
}

fn arb_fp_matrix() -> impl Strategy<Value = (i64, Vec<Vec<i64>>, usize)> {
    (arb_prime(), 0usize..=3, 0usize..=3).prop_flat_map(|(p, rows, cols)| {
        prop::collection::vec(arb_entry(), rows * cols).prop_map(move |flat| {
            let m = if cols == 0 { vec![Vec::new(); rows] } else { flat.chunks(cols).map(|row| row.to_vec()).collect() };
            (p, m, cols)
        })
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]
    /// Same claim as the exhaustive test, for primes past the i64-product overflow line
    /// and for residues outside 0..p (including i64::MIN / i64::MAX).
    #[test]
    fn rank_mod_p_equals_determinantal_rank(case in arb_fp_matrix()) {
        let (p, m, cols) = case;
        prop_assert_eq!(fp::rank_mod_p(&m, cols, p), det_rank(&m, cols, p), "p={} m={:?}", p, m);
    }
}

/// Component label of each vertex in the 1-skeleton (edges only; faces do not affect H⁰).
fn vertex_components(s: &Sheaf) -> Vec<usize> {
    let nv = s.complex.n_cells.first().copied().unwrap_or(0);
    let mut parent: Vec<usize> = (0..nv).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x { p[x] = p[p[x]]; x = p[x]; }
        x
    }
    if let Some(edges) = s.complex.verts.get(1) {
        for e in edges {
            if e.len() >= 2 {
                let (a, b) = (find(&mut parent, e[0]), find(&mut parent, e[1]));
                if a != b { parent[a] = b; }
            }
        }
    }
    let mut label = vec![usize::MAX; nv];
    let mut n = 0;
    (0..nv).map(|v| {
        let r = find(&mut parent, v);
        if label[r] == usize::MAX { label[r] = n; n += 1; }
        label[r]
    }).collect()
}

/// Constant sheaf: δ⁰ is the graph incidence matrix, once per stalk coordinate, and those
/// coordinates do not mix. A coordinate on a component survives in H⁰(X,B) exactly when
/// B misses it on every vertex of that component.
fn unpinned_constant_sections(s: &Sheaf, d: usize, boundary: &[usize]) -> usize {
    let comp = vertex_components(s);
    let ncomp = comp.iter().copied().max().map(|m| m + 1).unwrap_or(0);
    let mut hit = vec![vec![false; d]; ncomp];
    for &b in boundary {
        let (v, k) = (b / d, b % d);
        if v < comp.len() && k < d { hit[comp[v]][k] = true; }
    }
    hit.iter().map(|row| row.iter().filter(|h| !**h).count()).sum()
}

fn relative_h0(s: &Sheaf, boundary: &[usize]) -> sheaf_cohom::relative::RelativeH0 {
    let (rows, cols, t) = s.coboundary_triplets(0);
    h0_relative(&Csr::from_triplets(rows, cols, t), boundary, H0Policy::adr012(), 1500).unwrap().value
}

#[test]
fn relative_h0_unpinned_count_anchors() {
    let path = CellComplex::simplicial(4, &[(0, 1), (1, 2), (2, 3)], &[]);
    let s = Sheaf::constant(path, 1);
    assert_eq!(unpinned_constant_sections(&s, 1, &[]), 1);
    assert_eq!(relative_h0(&s, &[]).dim, 1);
    assert_eq!(relative_h0(&s, &[0]).dim, 0);
    assert_eq!(relative_h0(&s, &[0, 0]).dim, 0);
    // two paths, stalk R^2. Dofs are vertex-major: vertex v owns 2v and 2v+1.
    // Pinning coordinate 0 on vertex 0 leaves three sections; pinning it on both components leaves the two coordinate-1 sections.
    let split = CellComplex::simplicial(4, &[(0, 1), (2, 3)], &[]);
    let t = Sheaf::constant(split, 2);
    assert_eq!(unpinned_constant_sections(&t, 2, &[]), 4);
    assert_eq!(relative_h0(&t, &[]).dim, 4);
    assert_eq!(unpinned_constant_sections(&t, 2, &[0]), 3);
    assert_eq!(relative_h0(&t, &[0]).dim, 3);
    assert_eq!(unpinned_constant_sections(&t, 2, &[0, 2]), 3);
    assert_eq!(relative_h0(&t, &[0, 2]).dim, 3);
    assert_eq!(unpinned_constant_sections(&t, 2, &[0, 4]), 2);
    assert_eq!(relative_h0(&t, &[0, 4]).dim, 2);
    assert_eq!(unpinned_constant_sections(&t, 2, &[0, 1]), 2);
    assert_eq!(relative_h0(&t, &[0, 1]).dim, 2);
    assert_eq!(unpinned_constant_sections(&t, 2, &[0, 1, 2, 3]), 2);
    assert_eq!(relative_h0(&t, &[0, 1, 2, 3]).dim, 2);
}

fn arb_constant_relative() -> impl Strategy<Value = (Sheaf, usize, Vec<usize>)> {
    (arb_complex(), 1usize..3).prop_flat_map(|(cx, d)| {
        let n0 = cx.n_cells[0] * d;
        prop_oneof![
            1 => Just(Vec::<usize>::new()),
            1 => Just((0..n0).collect::<Vec<_>>()),
            4 => prop::collection::vec(0..n0, 0..=n0.min(6)),
        ].prop_map(move |boundary| (Sheaf::constant(cx.clone(), d), d, boundary))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]
    /// Constant sheaf: dim H⁰(X,B) equals the number of (1-skeleton component, stalk coordinate)
    /// pairs that B misses. The basis is orthonormal, supported off B, and δ⁰-closed.
    #[test]
    fn relative_h0_counts_unpinned_component_coordinates((s, d, boundary) in arb_constant_relative()) {
        let want = unpinned_constant_sections(&s, d, &boundary);
        let h = relative_h0(&s, &boundary);
        prop_assert_eq!(h.dim, want);
        prop_assert_eq!(h.basis.ncols(), h.dim);
        let n0 = s.cochain_dim(0);
        prop_assert_eq!(h.basis.nrows(), n0);
        let abs = CochainComplex::from_sheaf(&s).unwrap().cohomology(0).unwrap().dim;
        prop_assert!(h.dim <= abs);
        if boundary.is_empty() { prop_assert_eq!(h.dim, abs); }
        let (rows, cols, t) = s.coboundary_triplets(0);
        let d0 = Csr::from_triplets(rows, cols, t);
        let mut on_b = vec![false; n0];
        for &i in &boundary { on_b[i] = true; }
        let mut cols_g = Vec::with_capacity(h.dim);
        for j in 0..h.dim {
            let g: Vec<f64> = (0..n0).map(|i| h.basis[(i, j)]).collect();
            for (i, &on) in on_b.iter().enumerate() { if on { prop_assert_eq!(g[i], 0.0); } }
            let resid = norm(&d0.mul(&g));
            prop_assert!(resid < 1e-8, "||δg||={resid}");
            cols_g.push(g);
        }
        for j in 0..h.dim { for i in 0..=j {
            let want_dot = if i == j { 1.0 } else { 0.0 };
            prop_assert!((dot(&cols_g[i], &cols_g[j]) - want_dot).abs() < 1e-8, "basis columns {i},{j}");
        }}
    }
}
