use proptest::prelude::*;
use sheaf_cell::num::tol::{H0Policy, ZeroCut, EIG_TOL};

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]
    /// For every rel, t in [0, 1] and bounds / maxima b1, b2, m in [0, 1e12]:
    /// PerComponent{rel}.zero_tol(b) = rel·max(1, b) >= rel and is monotone in b; Absolute(t).zero_tol(b) = t;
    /// both have no global cut. GlobalSigma{rel} / GlobalLambda{rel}.global_cut(m) = Some(rel·max(1, m)), monotone in m.
    #[test]
    fn zero_cut_is_the_stated_scaled_rule_and_monotone(
        rel in 0f64..=1.0, t in 0f64..=1.0, b1 in 0f64..1e12, b2 in 0f64..1e12, m in 0f64..1e12,
    ) {
        let (lo, hi) = (b1.min(b2), b1.max(b2));
        let pc = H0Policy { cut: ZeroCut::PerComponent { rel } };
        prop_assert_eq!(pc.zero_tol(b1), rel * b1.max(1.0));
        prop_assert!(pc.zero_tol(b1) >= rel);
        prop_assert!(pc.zero_tol(lo) <= pc.zero_tol(hi));
        prop_assert_eq!(pc.global_cut(m), None);

        let ab = H0Policy::absolute(t);
        prop_assert_eq!(ab.zero_tol(b1), t);
        prop_assert_eq!(ab.zero_tol(b2), t);
        prop_assert_eq!(ab.global_cut(m), None);

        for g in [H0Policy::m1_svd(rel), H0Policy::m1_eig(rel)] {
            prop_assert_eq!(g.global_cut(m), Some(rel * m.max(1.0)));
            prop_assert!(g.global_cut(lo).unwrap() <= g.global_cut(hi).unwrap());
        }
    }
}

/// The global m1 variants have no per-component cut: `zero_tol` refuses them; the default is ADR-012; names are distinct.
#[test]
fn global_policies_refuse_per_component_cut_and_names_are_distinct() {
    let all = [H0Policy::adr012(), H0Policy::absolute(1e-8), H0Policy::m1_svd(1e-9), H0Policy::m1_eig(1e-10)];
    for p in &all[2..] { assert!(std::panic::catch_unwind(|| p.zero_tol(1.0)).is_err(), "{} gave a per-component cut", p.name()); }
    assert_eq!(H0Policy::default(), H0Policy { cut: ZeroCut::PerComponent { rel: EIG_TOL } });
    for i in 0..all.len() { for j in i + 1..all.len() { assert_ne!(all[i].name(), all[j].name()); } }
}
