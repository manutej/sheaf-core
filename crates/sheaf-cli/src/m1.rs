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

// ---- m1 math on Rust (ADR-016) ----------------------------------------------------------------------------------
use sheaf_cell::m1 as cell;
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::tol::{H0Policy, M1_RANK_TOL, EIG_TOL};

/// The m1 contract's quantities for one spec, computed by the Rust port (sheaf_cell::m1 + sheaf_cohom::m1).
#[derive(Debug, Clone)]
pub struct M1Math {
    pub delta: Csr,
    pub laplacian: Csr,
    /// per-edge discord, m1 edge order
    pub discord: Vec<f64>,
    /// xᵀ L₀ x
    pub energy: f64,
    /// dim H⁰ by the SVD route (`global_sections`, H0Policy::m1_svd(1e-9)) and by the eigen route (H0Policy::m1_eig(1e-10))
    pub h0_svd: usize,
    pub h0_eig: usize,
    pub policy_svd: String,
    pub policy_eig: String,
}

/// One m1 edge with owned row-major maps: (u, v, de, fu, fv).
pub type RawM1Edge = (usize, usize, usize, Vec<f64>, Vec<f64>);

/// Row-major flat maps per m1 edge: (dv, [(u, v, de, fu, fv)]).
pub fn to_m1_edges(s: &SheafSpec) -> (Vec<usize>, Vec<RawM1Edge>) {
    let e = s.edges.iter().map(|e| { let (fu, fv) = &s.maps[&e.ty];
        (e.u, e.v, fu.len(), fu.iter().flatten().copied().collect(), fv.iter().flatten().copied().collect()) }).collect();
    (vec![s.d; s.nodes.len()], e)
}

/// The m1 path. R-NONID-04 (and every other m1 rule) is enforced FIRST via `SheafSpec::validate`: the m1 contract never
/// defaults a map to the identity (maps are required per edge type, R-TYPE-03) and rejects an identity/identity pair,
/// unlike v1, whose omitted maps default to the identity (ADR-016; DECISIONS open item F8 resolved for the m1 path).
pub fn analyze(s: &SheafSpec) -> Result<M1Math, String> {
    s.validate()?;
    let (dv, raw) = to_m1_edges(s);
    let edges: Vec<cell::M1Edge> = raw.iter().map(|(u, v, de, fu, fv)| cell::M1Edge { u: *u, v: *v, de: *de, fu, fv }).collect();
    cell::check_edges(&dv, &edges).map_err(|e| e.to_string())?;
    let delta = cell::coboundary(&dv, &edges);
    let laplacian = cell::drop_zeros(cell::laplacian(&delta));
    let x: Vec<f64> = s.x.iter().flatten().copied().collect();
    let de: Vec<usize> = edges.iter().map(|e| e.de).collect();
    let discord = cell::edge_discord(&delta, &x, &de);
    let energy = cell::energy(&laplacian, &x);
    let (ps, pe) = (H0Policy::m1_svd(M1_RANK_TOL), H0Policy::m1_eig(EIG_TOL));
    let h0_svd = sheaf_cohom::m1::global_sections(&delta, ps).map_err(|e| e.to_string())?.ncols();
    let h0_eig = sheaf_cohom::m1::h0_dim_from_l(&laplacian, pe).map_err(|e| e.to_string())?;
    Ok(M1Math { delta, laplacian, discord, energy, h0_svd, h0_eig, policy_svd: ps.describe(), policy_eig: pe.describe() })
}
