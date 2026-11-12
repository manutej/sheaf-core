//! Flattening law for `sheaf_cli::m1::to_m1_edges`.
//! Each edge's maps come out row-major, so cutting the flat buffer every `d` entries gives the rows back.
use proptest::prelude::*;
use sheaf_cli::m1::to_m1_edges;
use sheaf_contract::{Edge, EdgeType, Kind, Node, SheafSpec};

const TYPES: [EdgeType; 5] = [EdgeType::About, EdgeType::Related, EdgeType::Logged, EdgeType::ClaimOf, EdgeType::SameSlug];

fn node() -> Node {
    Node { id: String::new(), kind: Kind::Topic, slug: None, title: None, status: None, url: None, parent: None, fetched: None }
}

/// A pair of `de × d` maps with entries in -2..=2, `de` below 4.
fn arb_pair(d: usize) -> impl Strategy<Value = (Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    (0usize..4).prop_flat_map(move |de| {
        let m = || prop::collection::vec(prop::collection::vec((-2i8..=2).prop_map(f64::from), d), de);
        (m(), m())
    })
}

/// `d` in 1..4, below 5 nodes, a map pair for every edge type, below 6 edges with any endpoints.
fn arb_spec() -> impl Strategy<Value = SheafSpec> {
    (1usize..4, 0usize..5).prop_flat_map(|(d, n)| {
        let maps = prop::collection::vec(arb_pair(d), TYPES.len());
        let edges = prop::collection::vec((any::<usize>(), any::<usize>(), 0..TYPES.len()), 0..6);
        (maps, edges).prop_map(move |(maps, edges)| SheafSpec {
            d,
            nodes: (0..n).map(|_| node()).collect(),
            x: vec![vec![0.0; d]; n],
            edges: edges.into_iter().map(|(u, v, t)| Edge { u, v, ty: TYPES[t] }).collect(),
            maps: TYPES.into_iter().zip(maps).collect(),
        })
    })
}

fn rows(flat: &[f64], d: usize) -> Vec<Vec<f64>> { flat.chunks(d).map(<[f64]>::to_vec).collect() }

proptest! {
    /// ∀ specs from `arb_spec` (every edge type has a `de × d` map pair, `d ≥ 1`; endpoints are not
    /// range-checked): `dv` is `d` repeated once per node, there is one tuple per edge in edge order,
    /// tuple j is `(u_j, v_j, de, fu, fv)` with `de` the row count of the edge type's maps, and
    /// `fu`, `fv` cut every `d` entries are exactly that type's two maps.
    #[test]
    fn to_m1_edges_flattens_maps_row_major(s in arb_spec()) {
        let (dv, raw) = to_m1_edges(&s);
        prop_assert_eq!(dv, vec![s.d; s.nodes.len()]);
        prop_assert_eq!(raw.len(), s.edges.len());
        for (e, (u, v, de, fu, fv)) in s.edges.iter().zip(&raw) {
            let (mu, mv) = &s.maps[&e.ty];
            prop_assert_eq!((*u, *v, *de), (e.u, e.v, mu.len()));
            prop_assert_eq!(&rows(fu, s.d), mu);
            prop_assert_eq!(&rows(fv, s.d), mv);
        }
    }
}
