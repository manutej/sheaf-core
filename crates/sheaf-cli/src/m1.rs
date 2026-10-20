//! m1 adapter (review §5: the only place sheaf-contract meets the cell/cohom stack). sheaf-contract stays
//! byte-identical and faer-free; this maps its `SheafSpec` onto a sheaf-cell graph sheaf, keeping m1's edge order,
//! parallel edges and orientation (edge j = 1-cell j). Parity is pinned by `tests/contract_parity.rs`.
use faer::Mat;
use sheaf_cell::graph::GraphEdge;
use sheaf_contract::SheafSpec;

fn mat(r: &[Vec<f64>]) -> Mat<f64> { Mat::from_fn(r.len(), r.first().map_or(0, |x| x.len()), |i, j| r[i][j]) }

/// (node stalk dims, edges) of an m1 spec. Caller should have run `SheafSpec::validate`.
pub fn to_graph(s: &SheafSpec) -> (Vec<usize>, Vec<GraphEdge>) {
    let edges = s.edges.iter().map(|e| {
        let (fu, fv) = &s.maps[&e.ty];
        GraphEdge { u: e.u, v: e.v, fu: mat(fu), fv: mat(fv), label: format!("{:?}", e.ty) }
    }).collect();
    (vec![s.d; s.nodes.len()], edges)
}
