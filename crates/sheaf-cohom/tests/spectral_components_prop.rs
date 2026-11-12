//! `spectral_components` solves each connected component on its own dofs and scatters every
//! pair counted nonzero back to global dofs. The scattered vector must be an eigenvector of the
//! whole L0 for its reported value, and it must live on one component.
use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::graph::{coboundary, components, GraphEdge};
use sheaf_cell::ops::HodgeOp;
use sheaf_cohom::spectral::{spectral_components, SpectralOpts};

/// One edge before its maps are cut to shape: endpoints, edge stalk dim, 8 entries in {-1, 0, 1}.
type RawEdge = (usize, usize, usize, [i8; 8]);

/// 2 to 6 nodes, node stalk dims in {1, 2}, up to 7 edges with u != v, edge stalk dim in {1, 2}.
fn arb_graph() -> impl Strategy<Value = (Vec<usize>, Vec<RawEdge>)> {
    (2usize..7).prop_flat_map(|n| {
        let edge = (0..n, 0..n, 1usize..3, prop::array::uniform8(-1i8..=1))
            .prop_map(move |(u, v, de, m)| if u == v { (u, (v + 1) % n, de, m) } else { (u, v, de, m) });
        (prop::collection::vec(1usize..3, n), prop::collection::vec(edge, 0..8))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// ∀ graphs from `arb_graph`, with F_u and F_v read row-major from the 8 entries, default opts:
    /// the first `found_zero` vectors are empty; every later vector has length `n`, unit norm,
    /// residual ‖L0 g − λ g‖ ≤ 1e-9 (1 + ‖L0‖_bound) against the GLOBAL L0, and support inside one
    /// connected component; the values counted nonzero ascend and number at most `k + 1`.
    #[test]
    fn nonzero_pairs_are_global_eigenpairs_on_one_component((dims, raw) in arb_graph()) {
        let edges: Vec<GraphEdge> = raw.iter().map(|&(u, v, de, m)| GraphEdge {
            u, v,
            fu: Mat::from_fn(de, dims[u], |r, c| m[r * 2 + c] as f64),
            fv: Mat::from_fn(de, dims[v], |r, c| m[4 + r * 2 + c] as f64),
            label: String::new(),
        }).collect();
        let so = SpectralOpts::default();
        let sp = spectral_components(&dims, &edges, &so).spectral;
        let (d0, on, _) = coboundary(&dims, &edges);
        let l0 = HodgeOp::l0(&d0);
        let (comp, _, _) = components(dims.len(), &edges);
        let n: usize = dims.iter().sum();

        prop_assert_eq!(sp.n, n);
        prop_assert_eq!(sp.pairs.len(), sp.vectors.len());
        prop_assert!(sp.found_zero <= sp.pairs.len());
        prop_assert!(sp.vectors[..sp.found_zero].iter().all(|v| v.is_empty()));
        prop_assert!(sp.pairs.len() - sp.found_zero <= so.k + 1);
        prop_assert!(sp.pairs[sp.found_zero..].windows(2).all(|w| w[0].value <= w[1].value));
        for (p, g) in sp.pairs.iter().zip(&sp.vectors).skip(sp.found_zero) {
            prop_assert_eq!(g.len(), n);
            let norm = g.iter().map(|x| x * x).sum::<f64>().sqrt();
            prop_assert!((norm - 1.0).abs() <= 1e-9, "norm {norm}");
            let res = l0.residual(g, p.value);
            prop_assert!(res <= 1e-9 * (1.0 + sp.lambda_max_bound), "residual {res} at value {}", p.value);
            let mut support = (0..dims.len()).filter(|&i| g[on[i]..on[i + 1]].iter().any(|&x| x != 0.0)).map(|i| comp[i]);
            let first = support.next();
            prop_assert!(first.is_some() && support.all(|c| Some(c) == first));
        }
    }
}
