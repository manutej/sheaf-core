use proptest::prelude::*;
use sheaf_cell::m1::{drop_zeros, laplacian, laplacian_nnz};
use sheaf_cell::num::sparse::Csr;

/// Random r x c CSR (1 <= r, c <= 5) with small-integer entries; repeated (row, col) triplets are kept as separate
/// entries and explicit zeros occur. Integer values keep every product and sum exact in f64.
fn arb_csr() -> impl Strategy<Value = Csr> {
    (1usize..6, 1usize..6).prop_flat_map(|(r, c)| {
        prop::collection::vec((0..r, 0..c, -3i8..=3), 0..20)
            .prop_map(move |t| Csr::from_triplets(r, c, t.into_iter().map(|(i, j, v)| (i, j, v as f64)).collect()))
    })
}

fn dense(a: &Csr) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; a.ncols]; a.nrows];
    for (i, row) in m.iter_mut().enumerate() { for p in a.ptr[i]..a.ptr[i + 1] { row[a.idx[p]] += a.val[p]; } }
    m
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// For every CSR δ: `laplacian(δ)` is n x n (n = δ.ncols) with strictly increasing columns per row, densely equals
    /// δᵀδ; `laplacian_nnz(δ)` is exactly its stored-entry count (the FFI sizes its buffer with it); `drop_zeros`
    /// keeps the dense value and stores no 0.0.
    #[test]
    fn m1_laplacian_is_dt_d_with_the_symbolic_nnz(d in arb_csr()) {
        let l = laplacian(&d);
        let n = d.ncols;
        prop_assert_eq!((l.nrows, l.ncols, l.ptr.len()), (n, n, n + 1));
        prop_assert_eq!(laplacian_nnz(&d), l.idx.len());
        prop_assert_eq!(l.val.len(), l.idx.len());
        for i in 0..n { prop_assert!(l.idx[l.ptr[i]..l.ptr[i + 1]].windows(2).all(|w| w[0] < w[1]), "row {} unsorted", i); }
        let (dd, ld) = (dense(&d), dense(&l));
        for (i, row) in ld.iter().enumerate() { for (j, &got) in row.iter().enumerate() {
            let want: f64 = dd.iter().map(|r| r[i] * r[j]).sum();
            prop_assert_eq!(got, want, "L[{}, {}]", i, j);
        }}
        let z = drop_zeros(l);
        prop_assert!(z.val.iter().all(|&v| v != 0.0));
        prop_assert_eq!(dense(&z), ld);
    }
}
