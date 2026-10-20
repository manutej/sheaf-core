//! The ONE zero-cut policy for H⁰ (ADR-012 per-component cut; ADR-015 makes it the only H⁰ rule).
//!
//! An eigenvalue λ of L₀ (equivalently a singular value σ of δ⁰ with λ = σ²) counts as zero iff
//!     λ ≤ cut(c) = rel · max(1, ‖L₀|_c‖_bound)        (`ZeroCut::PerComponent`, default rel = 1e-10)
//! where c is the connected component the eigenvector lives on and ‖·‖_bound is `HodgeOp::norm_bound`
//! (‖δ‖₁‖δ‖∞ on that component). L₀ is block diagonal over components, so each block gets its own cut.
//! `ZeroCut::Absolute(t)` is the per-context "permitted delta" (λ ≤ t on every component), chosen by the caller.
//! Every report echoes `H0Policy::describe()`. Superseded and removed from the v1/cohomology paths: the global R-H0-08 count (`dim_global_rule`,
//! re-admitted ONLY as the contract-named m1 variants below, ADR-016),
//! the per-component "strict" count (`dim_strict`), and the numpy-default SVD rank cut for δ⁰ in `cohomology(0)`.
//! Ranks of δ¹ (H¹/H² on 2-complexes) still use `linalg::rank_report` (numpy default); that is not an H⁰ rule.
//! The labelling is a numerical policy; eigenvalues carry certified residuals (Weyl), not a proof of the kernel.
//!
//! ADR-016: the gated m1 contract (SPEC v1.18 R-H0-08) is a NAMED variant of this same type, chosen by contract, so
//! the m1 port keeps the gated behaviour: `ZeroCut::GlobalSigma` (SVD route: σ is nonzero iff σ > rel·max(1, σ_max),
//! σ_max the largest singular value of δ over ALL components, default rel = 1e-9) and `ZeroCut::GlobalLambda` (eigen
//! route: λ is zero iff λ ≤ rel·max(1, λ_max), λ_max = max|λ| over all components, default rel = 1e-10). The values are
//! computed per component; only the cut is global. They are used by `sheaf_cohom::m1` only; the v1 / exact-cohomology
//! entry points keep `PerComponent` / `Absolute`.

/// Default relative cut (ADR-012).
pub const EIG_TOL: f64 = 1e-10;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ZeroCut {
    /// λ ≤ rel · max(1, ‖L₀|_c‖_bound), per connected component c (ADR-012).
    PerComponent { rel: f64 },
    /// λ ≤ t on every component: an explicit, context-chosen permitted delta.
    Absolute(f64),
    /// m1 R-H0-08 SVD route: σ nonzero iff σ > rel · max(1, σ_max over all components).
    GlobalSigma { rel: f64 },
    /// m1 R-H0-08 eigen route: λ zero iff λ ≤ rel · max(1, max|λ| over all components).
    GlobalLambda { rel: f64 },
}

/// R-H0-08 SVD rank tolerance (`RANK_TOL` in the reference).
pub const M1_RANK_TOL: f64 = 1e-9;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct H0Policy { pub cut: ZeroCut }

impl Default for H0Policy { fn default() -> Self { H0Policy { cut: ZeroCut::PerComponent { rel: EIG_TOL } } } }

impl H0Policy {
    pub const fn adr012() -> Self { H0Policy { cut: ZeroCut::PerComponent { rel: EIG_TOL } } }
    pub const fn absolute(t: f64) -> Self { H0Policy { cut: ZeroCut::Absolute(t) } }
    /// m1 contract, SVD route (R-H0-08, `global_sections`).
    pub const fn m1_svd(rel: f64) -> Self { H0Policy { cut: ZeroCut::GlobalSigma { rel } } }
    /// m1 contract, eigen route (R-H0-08, `h0_dim_from_L`).
    pub const fn m1_eig(rel: f64) -> Self { H0Policy { cut: ZeroCut::GlobalLambda { rel } } }
    /// Global variants: the cut for the GLOBAL maximum `global_max` (σ units for GlobalSigma, λ units for GlobalLambda).
    /// Returns None for the per-component variants.
    pub fn global_cut(&self, global_max: f64) -> Option<f64> {
        match self.cut { ZeroCut::GlobalSigma { rel } | ZeroCut::GlobalLambda { rel } => Some(rel * global_max.max(1.0)), _ => None }
    }
    /// Zero cut on λ for ONE component whose L₀ block has norm bound `component_bound`.
    /// # Panics
    /// On the global m1 variants (they have no per-component cut; use `global_cut`).
    pub fn zero_tol(&self, component_bound: f64) -> f64 {
        match self.cut {
            ZeroCut::PerComponent { rel } => rel * component_bound.max(1.0), ZeroCut::Absolute(t) => t,
            // the m1 variants have no per-component cut; callers of the per-component entry points must not pass them
            ZeroCut::GlobalSigma { .. } | ZeroCut::GlobalLambda { .. } => panic!("zero_tol: {} is a global m1 policy (use sheaf_cohom::m1)", self.name()),
        }
    }
    /// Machine-readable name echoed in reports.
    pub fn name(&self) -> &'static str {
        match self.cut { ZeroCut::PerComponent { .. } => "per_component_eig", ZeroCut::Absolute(_) => "absolute",
                         ZeroCut::GlobalSigma { .. } => "m1_global_svd", ZeroCut::GlobalLambda { .. } => "m1_global_eig" }
    }
    pub fn describe(&self) -> String {
        match self.cut {
            ZeroCut::PerComponent { rel } => format!("per component c: lambda <= {rel:e} * max(1, ||L0|_c||_bound) (ADR-012/ADR-015)"),
            ZeroCut::Absolute(t) => format!("absolute: lambda <= {t:e} on every component (caller-chosen permitted delta, ADR-015)"),
            ZeroCut::GlobalSigma { rel } => format!("m1 R-H0-08 SVD route: sigma nonzero iff sigma > {rel:e} * max(1, sigma_max over all components) (ADR-016)"),
            ZeroCut::GlobalLambda { rel } => format!("m1 R-H0-08 eigen route: lambda zero iff lambda <= {rel:e} * max(1, max|lambda| over all components) (ADR-016)"),
        }
    }
}
