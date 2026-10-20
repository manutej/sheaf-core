//! The ONE zero-cut policy for H⁰ (ADR-012 per-component cut; ADR-015 makes it the only H⁰ rule).
//!
//! An eigenvalue λ of L₀ (equivalently a singular value σ of δ⁰ with λ = σ²) counts as zero iff
//!     λ ≤ cut(c) = rel · max(1, ‖L₀|_c‖_bound)        (`ZeroCut::PerComponent`, default rel = 1e-10)
//! where c is the connected component the eigenvector lives on and ‖·‖_bound is `HodgeOp::norm_bound`
//! (‖δ‖₁‖δ‖∞ on that component). L₀ is block diagonal over components, so each block gets its own cut.
//! `ZeroCut::Absolute(t)` is the per-context "permitted delta" (λ ≤ t on every component), chosen by the caller.
//! Every report echoes `H0Policy::describe()`. Superseded and removed: the global R-H0-08 count (`dim_global_rule`),
//! the per-component "strict" count (`dim_strict`), and the numpy-default SVD rank cut for δ⁰ in `cohomology(0)`.
//! Ranks of δ¹ (H¹/H² on 2-complexes) still use `linalg::rank_report` (numpy default); that is not an H⁰ rule.
//! The labelling is a numerical policy; eigenvalues carry certified residuals (Weyl), not a proof of the kernel.

/// Default relative cut (ADR-012).
pub const EIG_TOL: f64 = 1e-10;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ZeroCut {
    /// λ ≤ rel · max(1, ‖L₀|_c‖_bound), per connected component c (ADR-012).
    PerComponent { rel: f64 },
    /// λ ≤ t on every component: an explicit, context-chosen permitted delta.
    Absolute(f64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct H0Policy { pub cut: ZeroCut }

impl Default for H0Policy { fn default() -> Self { H0Policy { cut: ZeroCut::PerComponent { rel: EIG_TOL } } } }

impl H0Policy {
    pub const fn adr012() -> Self { H0Policy { cut: ZeroCut::PerComponent { rel: EIG_TOL } } }
    pub const fn absolute(t: f64) -> Self { H0Policy { cut: ZeroCut::Absolute(t) } }
    /// Zero cut on λ for ONE component whose L₀ block has norm bound `component_bound`.
    pub fn zero_tol(&self, component_bound: f64) -> f64 {
        match self.cut { ZeroCut::PerComponent { rel } => rel * component_bound.max(1.0), ZeroCut::Absolute(t) => t }
    }
    /// Machine-readable name echoed in reports.
    pub fn name(&self) -> &'static str { match self.cut { ZeroCut::PerComponent { .. } => "per_component_eig", ZeroCut::Absolute(_) => "absolute" } }
    pub fn describe(&self) -> String {
        match self.cut {
            ZeroCut::PerComponent { rel } => format!("per component c: lambda <= {rel:e} * max(1, ||L0|_c||_bound) (ADR-012/ADR-015)"),
            ZeroCut::Absolute(t) => format!("absolute: lambda <= {t:e} on every component (caller-chosen permitted delta, ADR-015)"),
        }
    }
}
