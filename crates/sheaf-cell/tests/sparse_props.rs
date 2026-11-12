use proptest::prelude::*;
use sheaf_cell::complex::{CellComplex, Incidence};
use sheaf_cell::num::sparse::{transpose, Csr};
use sheaf_cell::num::dot;
use sheaf_cell::ops::{edge_discord, energy, HodgeOp};
use sheaf_cell::sheaf::Sheaf;

/// Random r x c CSR (1 <= r, c <= 5) with small-integer entries; repeated (row, col) triplets are frequent, and
/// `Csr` keeps them as separate entries. Integer values keep every product and sum exact in f64.
fn csr(r: usize, c: usize) -> impl Strategy<Value = Csr> {
    prop::collection::vec((0..r, 0..c, -3i8..=3), 0..20)
        .prop_map(move |t| Csr::from_triplets(r, c, t.into_iter().map(|(i, j, v)| (i, j, v as f64)).collect()))
}
fn arb_csr() -> impl Strategy<Value = Csr> { (1usize..6, 1usize..6).prop_flat_map(|(r, c)| csr(r, c)) }
/// (A, B) with A: n x m (the d_{k-1} of L_k) and B: k x n (the d_k), n = A.nrows.
fn arb_down_up() -> impl Strategy<Value = (Csr, Csr)> {
    (1usize..6, 1usize..6, 1usize..6).prop_flat_map(|(n, m, k)| (csr(n, m), csr(k, n)))
}

proptest! {
    /// transpose(transpose(A)) == A entry for entry, and A.mul_t(y) == transpose(A).mul(y) for every y.
    #[test]
    fn transpose_is_an_involution_and_mul_t_is_transpose_mul(a in arb_csr(), y in prop::collection::vec(-1e3f64..1e3, 5)) {
        let tt = transpose(&transpose(&a));
        prop_assert_eq!((tt.nrows, tt.ncols, &tt.ptr, &tt.idx, &tt.val), (a.nrows, a.ncols, &a.ptr, &a.idx, &a.val));
        let y = &y[..a.nrows];
        prop_assert_eq!(a.mul_t(y), transpose(&a).mul(y));
    }

    /// HodgeOp::diagonal() == diag(A A^T + B^T B) as applied by HodgeOp, including CSRs with duplicate entries.
    #[test]
    fn hodge_diagonal_is_the_diagonal_of_the_operator((a, b) in arb_down_up(), use_a: bool, use_b: bool) {
        let n = a.nrows;
        let op = HodgeOp::general(n, use_a.then_some(&a), use_b.then_some(&b));
        let dense = op.dense();
        prop_assert_eq!(op.diagonal(), (0..n).map(|i| dense[(i, i)]).collect::<Vec<_>>());
    }

    /// For any d0 and any edge_offsets partitioning 0..d0.nrows into contiguous blocks (0 = o_0 <= ... <= o_m = nrows):
    /// edge_discord has m entries, each >= 0, they sum to energy(d0, x), and energy == x^T L0 x for L0 = HodgeOp::l0.
    /// Integer entries keep every square, sum and product exact, so the equalities are bitwise.
    #[test]
    fn edge_discord_sums_to_energy_and_energy_is_the_l0_quadratic_form(
        a in arb_csr(), x in prop::collection::vec(-4i8..=4, 5), cuts in prop::collection::vec(0usize..6, 0..6),
    ) {
        let x: Vec<f64> = x[..a.ncols].iter().map(|&v| v as f64).collect();
        let mut offs: Vec<usize> = cuts.into_iter().map(|c| c.min(a.nrows)).collect();
        offs.extend([0, a.nrows]); offs.sort_unstable();
        let dis = edge_discord(&a, &x, &offs);
        let e = energy(&a, &x);
        prop_assert_eq!(dis.len(), offs.len() - 1);
        prop_assert!(dis.iter().all(|&v| v >= 0.0));
        prop_assert_eq!(dis.iter().sum::<f64>(), e);
        let mut l0x = vec![0.0; a.ncols];
        HodgeOp::l0(&a).apply(&x, &mut l0x);
        prop_assert_eq!(dot(&x, &l0x), e);
    }
}

/// CW circle (one vertex, one loop edge): [v:e] = -1 and +1 give duplicate δ⁰ entries that cancel, so L₀ = 0.
#[test]
fn loop_edge_has_zero_laplacian_diagonal() {
    let inc = vec![vec![Incidence { facet: 0, cofacet: 0, sign: -1 }, Incidence { facet: 0, cofacet: 0, sign: 1 }]];
    let s = Sheaf::constant(CellComplex { n_cells: vec![1, 1], inc, verts: vec![] }, 1);
    s.validate().unwrap();
    let (r, c, t) = s.coboundary_triplets(0);
    let op = HodgeOp::l0(&Csr::from_triplets(r, c, t));
    assert_eq!(op.dense()[(0, 0)], 0.0);
    assert_eq!(op.diagonal(), vec![0.0]);
}
