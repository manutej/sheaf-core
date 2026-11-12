use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::graph::{self, GraphEdge};

/// Restriction maps R^d -> R^{d_e} by kind: 0 F_u = F_v = I (injective, equal); 1 F_u = F_v = [I; I] (tall,
/// injective, equal); 2 F_u = F_v = I with column 0 zeroed (equal, NOT injective); 3 F_u = I, F_v = 2I (unequal).
/// The kind is stored in `label`. Entries are 0, 1, 2, so every δ⁰ product is exact in f64.
fn edge(u: usize, v: usize, d: usize, kind: u8) -> GraphEdge {
    let id = |r: usize, s: f64| Mat::from_fn(r, d, move |i, j| if i % d == j { s } else { 0.0 });
    let (fu, fv) = match kind {
        0 => (id(d, 1.0), id(d, 1.0)),
        1 => (id(2 * d, 1.0), id(2 * d, 1.0)),
        2 => { let m = Mat::from_fn(d, d, |i, j| if i == j && j > 0 { 1.0 } else { 0.0 }); (m.clone(), m) }
        _ => (id(d, 1.0), id(d, 2.0)),
    };
    GraphEdge { u, v, fu, fv, label: kind.to_string() }
}

/// n in 2..8 nodes, all of stalk dim d in 1..=3, up to 10 edges with u != v (parallel and u > v edges allowed).
fn arb_graph() -> impl Strategy<Value = (usize, usize, Vec<GraphEdge>)> {
    (2usize..8, 1usize..=3).prop_flat_map(|(n, d)| {
        let e = (0..n, 1..n, 0u8..4).prop_map(move |(u, k, kind)| edge(u, (u + k) % n, d, kind));
        (Just(n), Just(d), prop::collection::vec(e, 0..10))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]
    /// `components(n, E) = (comp, ok, nc)`: comp[i] == comp[j] iff i, j are joined by a path in E (checked against a
    /// BFS model); labels are 0..nc numbered by smallest node; ok[c] iff every edge of c has kind 0 or 1; and for
    /// every ok component, the constant cochain (the same vector on each of its nodes, zero elsewhere) has δ⁰x = 0.
    #[test]
    fn components_partition_and_known_constants_are_sections((n, d, es) in arb_graph(), c in prop::collection::vec(-4i8..=4, 3)) {
        let (comp, ok, nc) = graph::components(n, &es);
        let mut adj = vec![vec![]; n];
        for e in &es { adj[e.u].push(e.v); adj[e.v].push(e.u); }
        let mut model = vec![usize::MAX; n]; let mut m = 0;
        for s in 0..n {
            if model[s] != usize::MAX { continue; }
            let mut stack = vec![s]; model[s] = m;
            while let Some(a) = stack.pop() { for &b in &adj[a] { if model[b] == usize::MAX { model[b] = m; stack.push(b); } } }
            m += 1;
        }
        prop_assert_eq!(&comp, &model);
        prop_assert_eq!((nc, ok.len()), (m, m));
        for (k, &got) in ok.iter().enumerate() {
            let want = es.iter().filter(|e| comp[e.u] == k).all(|e| e.label == "0" || e.label == "1");
            prop_assert_eq!(got, want, "component {}", k);
        }
        let (d0, on, _) = graph::coboundary(&vec![d; n], &es);
        for k in (0..nc).filter(|&k| ok[k]) {
            let mut x = vec![0.0; d0.ncols];
            for i in (0..n).filter(|&i| comp[i] == k) { for j in 0..d { x[on[i] + j] = c[j] as f64; } }
            prop_assert!(d0.mul(&x).iter().all(|&v| v == 0.0), "constants of known-constant component {} are not sections", k);
        }
    }
}
