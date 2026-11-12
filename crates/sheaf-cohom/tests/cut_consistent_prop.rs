//! `map_or(true, f)` and `is_none_or(f)` agree on the zero side of the residual cut.
//! None stays consistent. Some still evaluates `value + residual <= zero_tol` once.
use proptest::prelude::*;
use sheaf_cohom::spectral::{EigPair, Spectral};

fn arb_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        8 => -1.0e6f64..1.0e6,
        1 => Just(0.0f64),
        1 => Just(-0.0f64),
        1 => Just(f64::INFINITY),
        1 => Just(f64::NEG_INFINITY),
        1 => Just(f64::NAN),
        1 => Just(f64::MIN),
        1 => Just(f64::MAX),
    ]
}

fn arb_case() -> impl Strategy<Value = (Vec<EigPair>, usize, f64, bool)> {
    prop::collection::vec(arb_f64().prop_flat_map(|value| arb_f64().prop_map(move |residual| EigPair { value, residual })), 0..6)
        .prop_flat_map(|pairs| {
            let n = pairs.len();
            (Just(pairs), 0..=n, arb_f64(), any::<bool>())
        })
}

/// Match form of the predicate before the `is_none_or` rewrite.
/// `None` on the zero side is the old `map_or(true, _)` default.
fn match_oracle(pairs: &[EigPair], found_zero: usize, zero_tol: f64, h0_lower_bound_only: bool) -> bool {
    let zero_ok = match pairs[..found_zero].last() {
        None => true,
        Some(p) => p.value + p.residual <= zero_tol,
    };
    let nonzero_ok = match pairs.get(found_zero) {
        None => !h0_lower_bound_only,
        Some(p) => p.value - p.residual > zero_tol,
    };
    zero_ok && nonzero_ok
}

fn spectral(pairs: Vec<EigPair>, found_zero: usize, zero_tol: f64, h0_lower_bound_only: bool) -> Spectral {
    Spectral {
        method: "prop",
        n: 0,
        known_kernel_dim: 0,
        pairs,
        vectors: Vec::new(),
        c_shift: 0.0,
        lambda_max_bound: 0.0,
        zero_tol,
        found_zero,
        rounds: 0,
        inner: None,
        sigma: 0.0,
        converged: false,
        h0_lower_bound_only,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// ∀ pairs with length < 6, ∀ found_zero ≤ pairs.len(), ∀ zero_tol ∈ f64
    /// (finite, ±0, ±∞, NaN), ∀ h0_lower_bound_only:
    /// `cut_consistent_with_residuals` equals the match-form of
    /// `map_or(true, |p| p.value + p.residual <= zero_tol)` on the zero side
    /// and the unchanged nonzero-side `map_or`.
    #[test]
    fn cut_consistent_matches_option_match((pairs, found_zero, zero_tol, lower) in arb_case()) {
        let want = match_oracle(&pairs, found_zero, zero_tol, lower);
        let got = spectral(pairs, found_zero, zero_tol, lower).cut_consistent_with_residuals();
        prop_assert_eq!(got, want);
    }
}
