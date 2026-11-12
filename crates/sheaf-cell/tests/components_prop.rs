//! Partition law for `graph::components`.
//! Labels are the undirected connected components, numbered by minimum vertex.
use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::graph::{self, GraphEdge};

/// BFS labels in the order a scan from vertex 0 assigns them.
/// A new id starts at the least unvisited vertex, so ids follow minimum vertex index.
/// Self-loops are not adjacencies. Parallel edges do not change the adjacency.
fn min_index_bfs(n: usize, edges: &[(usize, usize)]) -> (Vec<usize>, usize) {
    let mut adj = vec![Vec::new(); n];
    for &(u, v) in edges {
        if u != v {
            adj[u].push(v);
            adj[v].push(u);
        }
    }
    let mut comp = vec![usize::MAX; n];
    let mut nc = 0usize;
    for i in 0..n {
        if comp[i] != usize::MAX {
            continue;
        }
        comp[i] = nc;
        let mut q = vec![i];
        let mut h = 0;
        while h < q.len() {
            let x = q[h];
            h += 1;
            for &y in &adj[x] {
                if comp[y] == usize::MAX {
                    comp[y] = nc;
                    q.push(y);
                }
            }
        }
        nc += 1;
    }
    (comp, nc)
}

/// Every n in `0..12`, and every list of endpoints inside `0..n` (self-loops and parallel edges included).
fn arb_endpoints() -> impl Strategy<Value = (usize, Vec<(usize, usize)>)> {
    prop_oneof![
        1 => Just((0usize, Vec::new())),
        11 => (1usize..12).prop_flat_map(|n| {
            prop::collection::vec((0..n, 0..n), 0..24).prop_map(move |ends| (n, ends))
        }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// ∀ n, ∀ edges with endpoints in `0..n`:
    /// `components` labels equal min-index BFS, the count equals that partition's size,
    /// `ok.len() == nc`, and every edge stays inside one label.
    #[test]
    fn components_labels_match_min_index_bfs((n, ends) in arb_endpoints()) {
        let edges: Vec<GraphEdge> = ends.iter().map(|&(u, v)| GraphEdge {
            u,
            v,
            fu: Mat::<f64>::zeros(0, 0),
            fv: Mat::<f64>::zeros(0, 0),
            label: String::new(),
        }).collect();
        let (comp, ok, nc) = graph::components(n, &edges);
        let (want, want_nc) = min_index_bfs(n, &ends);
        prop_assert_eq!(&comp, &want);
        prop_assert_eq!(nc, want_nc);
        prop_assert_eq!(ok.len(), nc);
        for &(u, v) in &ends {
            prop_assert_eq!(comp[u], comp[v]);
        }
    }
}

/// Every directed graph on `n <= 3` vertices, self-loops and parallel edges included
/// (each possible endpoint pair is either absent or present once).
#[test]
fn components_labels_match_bfs_for_every_graph_on_at_most_3_vertices() {
    for n in 0..4 {
        let pairs: Vec<(usize, usize)> = (0..n).flat_map(|u| (0..n).map(move |v| (u, v))).collect();
        let masks = 1usize << pairs.len();
        for mask in 0..masks {
            let ends: Vec<(usize, usize)> = pairs.iter().copied().enumerate().filter(|(i, _)| mask & (1 << i) != 0).map(|(_, e)| e).collect();
            let edges: Vec<GraphEdge> = ends.iter().map(|&(u, v)| GraphEdge {
                u, v, fu: Mat::<f64>::zeros(0, 0), fv: Mat::<f64>::zeros(0, 0), label: String::new(),
            }).collect();
            let (comp, ok, nc) = graph::components(n, &edges);
            let (want, want_nc) = min_index_bfs(n, &ends);
            assert_eq!(comp, want);
            assert_eq!(nc, want_nc);
            assert_eq!(ok.len(), nc);
            for &(u, v) in &ends {
                assert_eq!(comp[u], comp[v]);
            }
        }
    }
}
