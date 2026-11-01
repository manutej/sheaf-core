//! Typed v1 report (review F10 / ADR-004: the types are the source of truth). Serialised through
//! `serde_json::Value` so object keys come out sorted (deterministic, as before); floats are pre-rounded to
//! 12 significant digits by the core, and NaN/inf serialise as null.
//! ADR-013: every field here is the exact-cohomology family ("method": "exact"); no ordered/interval
//! approximation type exists in this crate, so none can be written into these fields.
use serde::Serialize;

/// The only family value this report can carry (ADR-013): the `Exact` marker type from sheaf-cell (serialises as
/// "exact"); an `Approximate<T>` cannot be stored in any field here.
pub use sheaf_cell::num::family::Exact;

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Report {
    pub format: &'static str,
    pub method: Exact,
    pub input: Input,
    pub h0: H0Report,
    pub h1: H1Report,
    pub spectral: SpectralReport,
    pub hodge: HodgeReport,
    #[serde(rename = "energy_xLx", skip_serializing_if = "Option::is_none")]
    pub energy_xlx: Option<f64>,
    pub edges_by_harmonic_energy: Vec<EdgeRow>,
    pub edges_by_discord: Vec<EdgeRow>,
    pub nodes_by_abs_potential: Vec<NodePotential>,
    pub nodes_by_incident_discord: Vec<NodeDiscord>,
    /// present only with `--boundary` (harmonic extension, sheaf-solve + uniqueness via sheaf-cohom::relative)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extension: Option<ExtensionReport>,
    pub assumptions: Vec<&'static str>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Input {
    pub nodes: usize,
    pub edges: usize,
    #[serde(rename = "dim_C0")] pub dim_c0: usize,
    #[serde(rename = "dim_C1")] pub dim_c1: usize,
    pub components: usize,
    pub known_constant_components: usize,
    pub has_x: bool,
    pub has_w: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Pair { pub value: f64, pub residual: f64 }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct H0Report {
    pub method: Exact,
    pub dim: usize,
    pub known_kernel_dim: usize,
    pub found_zero_by_solver: usize,
    /// the policy echoed: name + human-readable rule (ADR-015)
    pub zero_policy: &'static str,
    pub zero_rule: String,
    /// largest / smallest per-component cut actually applied
    pub zero_tol: f64,
    pub zero_tol_min: f64,
    pub largest_counted_zero: Option<Pair>,
    pub smallest_counted_nonzero: Option<Pair>,
    pub cut_consistent_with_zero_tol: bool,
    pub near_zero_ambiguous: bool,
    pub cut_margin_ratio: Option<f64>,
    pub lower_bound_only: bool,
    /// solver certainty for h0, separate from near_zero_ambiguous (sheaf-glue-graph-001 :292)
    pub status: H0Status,
    pub wording: &'static str,
}

/// `lower_bound_only` wins over `not_converged`; `converged` only when the solver converged and h0 is not a lower bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum H0Status { Converged, NotConverged, LowerBoundOnly }
impl H0Status {
    pub fn from_flags(converged: bool, lower_bound_only: bool) -> Self {
        if lower_bound_only { H0Status::LowerBoundOnly } else if !converged { H0Status::NotConverged } else { H0Status::Converged }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct H1Report { pub method: Exact, pub dim: usize, pub rule: &'static str }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct InnerCg { pub solves: usize, pub iters: usize, pub worst_rel_residual: f64 }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct L1Graph { pub kernel_dim: usize, pub nonzero_spectrum: &'static str, pub mapped_pairs_checked: usize, pub max_mapped_residual: f64 }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SpectralReport {
    pub method: Exact,
    pub solver: &'static str,
    #[serde(rename = "L0_smallest")] pub l0_smallest: Vec<f64>,
    pub solver_pairs_residuals: Vec<f64>,
    pub gap_k: Option<f64>,
    pub spectral_gap: Option<f64>,
    pub lambda_max_bound: f64,
    pub rounds: usize,
    pub converged: bool,
    pub inner_cg: Option<InnerCg>,
    #[serde(rename = "L1_graph")] pub l1_graph: L1Graph,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct NullRandom {
    pub reps: usize,
    pub mean: Option<f64>,
    pub p05: f64,
    pub p95: f64,
    pub frac_ge_observed: Option<f64>,
    #[serde(rename = "analytic_expectation_h1_over_dimC1")] pub analytic_expectation_h1_over_dim_c1: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct NullShuffle {
    pub degenerate: bool,
    pub note: &'static str,
    pub reps: usize,
    pub mean: Option<f64>,
    pub p05: f64,
    pub p95: f64,
    pub frac_ge_observed: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HodgeReport {
    pub method: Exact,
    pub cochain: &'static str,
    pub norm2: f64,
    pub exact_share: f64,
    pub harmonic_share: f64,
    pub cg_iters: usize,
    pub cert_coclosed: f64,
    pub converged: bool,
    pub null_random_gaussian: NullRandom,
    pub null_shuffled_values: NullShuffle,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct EdgeRow { pub edge: usize, pub u: String, pub v: String, pub label: String, pub harmonic_energy: f64, pub discord: f64, pub harmonic_frac_of_total: f64 }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct NodePotential { pub node: String, pub abs_potential: f64 }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct NodeDiscord { pub node: String, pub incident_discord: f64 }

impl Report {
    /// Deterministic JSON (sorted keys via `Value`, pretty, trailing newline).
    pub fn to_json(&self) -> String {
        let v = serde_json::to_value(self).expect("report serialises");
        serde_json::to_string_pretty(&v).expect("pretty") + "\n"
    }
    pub fn to_value(&self) -> serde_json::Value { serde_json::to_value(self).expect("report serialises") }
}

/// Harmonic extension of boundary data (`--boundary`): x_I minimising ‖δx‖² with x_B fixed (sheaf-solve), plus the
/// uniqueness certificate H⁰(X,B;F) (sheaf-cohom::relative) under the same H0 policy as `h0`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ExtensionReport {
    pub method: Exact,
    pub boundary_nodes: usize,
    pub boundary_dofs: usize,
    /// dim H⁰(X,B;F); None if a block exceeded dense_max (not computed, never guessed)
    pub h0_relative: Option<usize>,
    /// h0_relative == 0
    pub unique: Option<bool>,
    pub uniqueness_rule: String,
    pub note: &'static str,
    pub route: sheaf_solve::Route,
    pub converged: bool,
    pub residual_rel: f64,
    pub iters: usize,
    pub energy: f64,
    pub values: Vec<NodeValue>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct NodeValue { pub id: String, pub boundary: bool, pub value: Vec<f64> }
