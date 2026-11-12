use proptest::prelude::*;
use sheaf_cell::num::linalg::rank_report;

/// Singular-value lists: each entry is 0 or a normal f64 in [1e-20, 1e3], so no subnormals reach `gap_ratio`.
fn arb_sigmas() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(prop_oneof![Just(0.0), 1e-20f64..1e3], 0..8)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]
    /// For every s (entries 0 or normal, finite), m, n >= 1 and override t (None, an entry of s, or any t >= 0):
    /// tol = t or max(m,n)·eps·max(s ∪ {0}); rank = #{x > tol}; kept_min / dropped_max are the min of the kept and the
    /// max of the dropped entries (None iff that side is empty), so kept_min > tol >= dropped_max; gap_ratio >= 1;
    /// and the report does not depend on the order of s.
    #[test]
    fn rank_report_splits_sigmas_at_the_tolerance(
        s in arb_sigmas(), m in 1usize..10, n in 1usize..10, pick in prop::option::of(any::<prop::sample::Index>()), raw in 0f64..1e3,
    ) {
        let t = pick.map(|i| if s.is_empty() { raw } else { s[i.index(s.len())] });
        let r = rank_report(&s, m, n, t);
        let smax = s.iter().cloned().fold(0.0, f64::max);
        prop_assert_eq!(r.sigma_max, smax);
        prop_assert_eq!(r.tol, t.unwrap_or(m.max(n) as f64 * f64::EPSILON * smax));
        let (kept, dropped): (Vec<f64>, Vec<f64>) = s.iter().partition(|&&x| x > r.tol);
        prop_assert_eq!(r.rank, kept.len());
        prop_assert_eq!(r.sigma_kept_min, kept.iter().cloned().reduce(f64::min));
        prop_assert_eq!(r.sigma_dropped_max, dropped.iter().cloned().reduce(f64::max));
        if let Some(k) = r.sigma_kept_min { prop_assert!(k > r.tol); }
        if let Some(d) = r.sigma_dropped_max { prop_assert!(d <= r.tol); }
        prop_assert!(r.gap_ratio() >= 1.0, "gap_ratio {} < 1", r.gap_ratio());
        let rev: Vec<f64> = s.iter().rev().cloned().collect();
        prop_assert_eq!(rank_report(&rev, m, n, t), r);
    }
}
