//! 1-skeleton component law for `CochainComplex::from_sheaf`.
//! Labels are the undirected connected components, numbered by minimum vertex.
use faer::Mat;
use proptest::prelude::*;
use sheaf_cell::sheaf::Sheaf;
use sheaf_cohom::cohomology::CochainComplex;

/// BFS labels in the order a scan from vertex 0 assigns them.
/// A new id starts at the least unvisited vertex, so ids follow minimum vertex index.
/// Self-loops are not adjacencies. Parallel edges do not change the adjacency.
fn min_index_bfs(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
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
    comp
}

/// Every n in `0..10`, vertex stalk dims in `1..4`, and endpoint lists inside `0..n` (no self-loops).
/// Edge stalk dimension `de` is shared by every edge.
fn arb_graph() -> impl Strategy<Value = (Vec<usize>, Vec<(usize, usize)>, usize)> {
    prop_oneof![
        1 => Just((Vec::new(), Vec::new(), 1usize)),
        1 => (1usize..4).prop_map(|d| (vec![d], Vec::new(), 1usize)),
        8 => (2usize..10).prop_flat_map(|n| {
            (
                prop::collection::vec(1usize..4, n),
                prop::collection::vec(
                    (0..n, 0..n).prop_map(move |(u, v)| if u == v { (u, (v + 1) % n) } else { (u, v) }),
                    0..12,
                ),
                1usize..3,
            )
        }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// ∀ vertex dims, ∀ edges with endpoints in range and u ≠ v, ∀ edge stalk dim `de`:
    /// `n_comp` and `c0_comp` equal min-index BFS of the 1-skeleton (each vertex stalk
    /// repeated `dims[v]` times), every edge stays in one label, and `c1_comp` is that
    /// label repeated `de` times per edge, in edge order. `CellComplex::graph` records
    /// the first endpoint first, and the edge joins the two endpoints, so the label is
    /// the component of either endpoint.
    #[test]
    fn skeleton_components_match_min_index_bfs((dims, ends, de) in arb_graph()) {
        let n = dims.len();
        let edges: Vec<_> = ends.iter().map(|&(u, v)| {
            (u, v, Mat::<f64>::zeros(de, dims[u]), Mat::<f64>::zeros(de, dims[v]))
        }).collect();
        let s = Sheaf::from_graph(&dims, &edges).unwrap();
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        let vcomp = min_index_bfs(n, &ends);
        let nc = vcomp.iter().copied().max().map(|m| m + 1).unwrap_or(0);
        prop_assert_eq!(cc.n_comp, nc);
        let c0: Vec<usize> = (0..n).flat_map(|v| std::iter::repeat(vcomp[v]).take(dims[v])).collect();
        prop_assert_eq!(&cc.c0_comp, &c0);
        for &(u, v) in &ends {
            prop_assert_eq!(vcomp[u], vcomp[v]);
        }
        let c1: Vec<usize> = ends.iter().flat_map(|&(u, _)| std::iter::repeat(vcomp[u]).take(de)).collect();
        prop_assert_eq!(&cc.c1_comp, &c1);
    }
}
