use proptest::prelude::*;
use sheaf_cell::num::krylov::{cg, SolveTol};
use sheaf_cell::num::norm;
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::ops::HodgeOp;

/// Random r x c CSR (1 <= r, c <= 5) with small-integer entries, plus an integer vector y of length c.
fn arb_system() -> impl Strategy<Value = (Csr, Vec<f64>)> {
    (1usize..6, 1usize..6).prop_flat_map(|(r, c)| {
        let d = prop::collection::vec((0..r, 0..c, -2i8..=2), 0..15)
            .prop_map(move |t| Csr::from_triplets(r, c, t.into_iter().map(|(i, j, v)| (i, j, v as f64)).collect()));
        (d, prop::collection::vec((-4i8..=4).prop_map(f64::from), c))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// For A = D^T D (symmetric PSD, possibly singular) and a consistent right-hand side b = A y, CG from x0 = 0
    /// returns x with A x = b (true residual, not the recurrence's) and ||x|| <= ||y||: x is the minimum-norm
    /// solution, so it is no longer than the solution y the system was built from. Two runs agree bitwise.
    #[test]
    fn cg_from_zero_solves_consistent_psd_systems_with_minimum_norm((d, y) in arb_system()) {
        let apply = |v: &[f64]| d.mul_t(&d.mul(v));
        let b = apply(&y);
        let op_scale = HodgeOp::l0(&d).norm_bound();
        let run = || cg(&apply, &b, vec![0.0; d.ncols], norm(&b), op_scale, SolveTol::default());
        let out = run();
        let r: Vec<f64> = b.iter().zip(apply(&out.x)).map(|(p, q)| p - q).collect();
        prop_assert!(norm(&r) <= 1e-9 * norm(&b).max(1.0), "true residual {} (b {}, {:?})", norm(&r), norm(&b), out);
        prop_assert!(norm(&out.x) <= norm(&y) * (1.0 + 1e-9) + 1e-9, "||x|| {} > ||y|| {}", norm(&out.x), norm(&y));
        prop_assert_eq!(run().x, out.x);
    }
}
