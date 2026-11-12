//! `smallest_dense` takes `k` only as a cut: it returns the first `min(k + 1, n)` pairs of the
//! full ascending spectrum. `examples/spectral_errors.rs` picks that `k` as `(n - 1).clamp(1, 6)`.
use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::graph::{coboundary, GraphEdge};
use sheaf_cohom::spectral::{smallest_dense, H0Policy, HodgeOp, KnownKernel};

/// One edge: endpoints and 2 entries in {-1, 0, 1} (the 1x1 maps F_u, F_v).
type RawEdge = (usize, usize, [i8; 2]);

/// 2 to 6 nodes with 1-dimensional stalks, up to 7 edges with u != v.
fn arb_graph() -> impl Strategy<Value = (usize, Vec<RawEdge>)> {
    (2usize..7).prop_flat_map(|n| {
        let edge = (0..n, 0..n, prop::array::uniform2(-1i8..=1))
            .prop_map(move |(u, v, m)| if u == v { (u, (v + 1) % n, m) } else { (u, v, m) });
        (Just(n), prop::collection::vec(edge, 0..8))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    /// ∀ graphs from `arb_graph` (so L0 is n x n, 2 ≤ n ≤ 6), no known kernel, ∀ k in 0..=n+1:
    /// `smallest_dense(.., k, ..)` returns exactly `min(k + 1, n)` pairs and as many vectors, the
    /// values ascend, and pair j and vector j are bit-identical to pair j and vector j of the
    /// full run `smallest_dense(.., n, ..)`.
    #[test]
    fn k_only_truncates_the_ascending_spectrum((n, raw) in arb_graph()) {
        let edges: Vec<GraphEdge> = raw.iter().map(|&(u, v, m)| GraphEdge {
            u, v,
            fu: Mat::from_fn(1, 1, |_, _| m[0] as f64),
            fv: Mat::from_fn(1, 1, |_, _| m[1] as f64),
            label: String::new(),
        }).collect();
        let (d0, _, _) = coboundary(&vec![1; n], &edges);
        let op = HodgeOp::l0(&d0);
        prop_assert_eq!(op.n, n);
        let none = KnownKernel::none(n);
        let full = smallest_dense(&op, &none, n, H0Policy::adr012());
        prop_assert_eq!(full.pairs.len(), n);
        for k in 0..=n + 1 {
            let sp = smallest_dense(&op, &none, k, H0Policy::adr012());
            let m = (k + 1).min(n);
            prop_assert_eq!(sp.pairs.len(), m);
            prop_assert_eq!(sp.vectors.len(), m);
            prop_assert!(sp.pairs.windows(2).all(|w| w[0].value <= w[1].value));
            for j in 0..m {
                prop_assert_eq!(sp.pairs[j].value.to_bits(), full.pairs[j].value.to_bits());
                let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
                prop_assert_eq!(bits(&sp.vectors[j]), bits(&full.vectors[j]));
            }
        }
    }
}

/// ∀ x in 0..=4096 and x = usize::MAX: `x.min(6).max(1) == x.clamp(1, 6)`, the rewrite made in
/// `examples/spectral_errors.rs` (both sides are 1 below 1, 6 above 6, and x in between).
#[test]
fn clamp_1_6_equals_min_6_then_max_1() {
    for x in (0..=4096usize).chain([usize::MAX]) {
        #[allow(clippy::manual_clamp)]
        let old = x.min(6).max(1);
        assert_eq!(old, x.clamp(1, 6));
    }
}
