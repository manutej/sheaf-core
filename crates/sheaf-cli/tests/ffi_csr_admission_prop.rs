//! CSR admission law of the C ABI, stated on `sheaf_m1_energy_csr`.
//! The `ffi` module docs promise that every entry point re-validates the CSR structure and
//! returns -2 on malformed input. Values are small integers, so the f64 arithmetic is exact.
use proptest::prelude::*;
use sheaf_cli::ffi::sheaf_m1_energy_csr;

const E_INPUT: i64 = -2;
const SENTINEL: f64 = 12345.5;

/// Returns (return code, value left in the output slot, which starts as `SENTINEL`).
fn call(m: i64, n: i64, nnz: i64, ptr: &[i64], idx: &[i64], val: &[f64], x: &[f64]) -> (i64, f64) {
    let mut out = SENTINEL;
    // SAFETY: every caller below passes either the true sizes of these slices (ptr has m + 1
    // entries, idx and val have nnz, x has n) or one negative size, which is refused before any read.
    let rc = unsafe { sheaf_m1_energy_csr(m, n, nnz, ptr.as_ptr(), idx.as_ptr(), val.as_ptr(), x.as_ptr(), &mut out) };
    (rc, out)
}

/// (n, rows of (column, value) entries, x). Columns may repeat and are unsorted; zeros may be stored.
type Case = (usize, Vec<Vec<(usize, i8)>>, Vec<i8>);

fn arb_case() -> impl Strategy<Value = Case> {
    (0usize..=5).prop_flat_map(|n| (
        Just(n),
        prop::collection::vec(prop::collection::vec((0..n.max(1), -2i8..=2), 0..=4), n),
        prop::collection::vec(-3i8..=3, n),
    ))
}

proptest! {
    /// ∀ n ≤ 5, ∀ n×n CSR with entries in {-2..2} (≤ 4 per row), ∀ x in {-3..3}^n:
    /// the well-formed call returns 0 and writes exactly xᵀLx;
    /// every single corruption below returns -2 and leaves the output slot untouched:
    /// indptr[0] ≠ 0, indptr[m] ≠ nnz, indptr decreasing at any interior position,
    /// any one column index set to -1 or to n, or a negative m, n or nnz.
    #[test]
    fn energy_csr_admits_wellformed_and_refuses_each_corruption((n, rows, x) in arb_case()) {
        let mut ptr = vec![0i64];
        let (mut idx, mut val) = (Vec::new(), Vec::new());
        let mut want = 0i64;
        for (i, row) in rows.iter().enumerate() {
            for &(j, v) in row {
                idx.push(j as i64);
                val.push(v as f64);
                want += x[i] as i64 * v as i64 * x[j] as i64;
            }
            ptr.push(idx.len() as i64);
        }
        let xf: Vec<f64> = x.iter().map(|&v| v as f64).collect();
        let (m, nnz) = (n as i64, idx.len() as i64);

        let (rc, out) = call(m, m, nnz, &ptr, &idx, &val, &xf);
        prop_assert_eq!(rc, 0);
        prop_assert_eq!(out, want as f64);

        let refused = |m: i64, n: i64, nnz: i64, ptr: &[i64], idx: &[i64]| {
            let (rc, out) = call(m, n, nnz, ptr, idx, &val, &xf);
            rc == E_INPUT && out.to_bits() == SENTINEL.to_bits()
        };
        let mut p = ptr.clone();
        p[0] = 1;
        prop_assert!(refused(m, m, nnz, &p, &idx), "indptr[0] = 1");
        let mut p = ptr.clone();
        p[n] = nnz + 1;
        prop_assert!(refused(m, m, nnz, &p, &idx), "indptr[m] = nnz + 1");
        for k in 1..n {
            let mut p = ptr.clone();
            p[k] = p[k - 1] - 1;
            prop_assert!(refused(m, m, nnz, &p, &idx), "indptr decreasing at {}", k);
        }
        for q in 0..idx.len() {
            for bad in [-1, m] {
                let mut ix = idx.clone();
                ix[q] = bad;
                prop_assert!(refused(m, m, nnz, &ptr, &ix), "index {} set to {}", q, bad);
            }
        }
        prop_assert!(refused(-1, m, nnz, &ptr, &idx), "m = -1");
        prop_assert!(refused(m, -1, nnz, &ptr, &idx), "n = -1");
        prop_assert!(refused(m, m, -1, &ptr, &idx), "nnz = -1");
    }
}
