//! Policies in types (ADR-013, review F4).
//! * Family: `Exact` = the exact-cohomology family (linear algebra on the exact complex; the numbers are still f64
//!   with an explicit policy cut and certified residuals). `Approximate` = the ordered/interval layer (ADR-013),
//!   which must never be written into an exact field. Report fields that hold exact results are typed with the
//!   marker `Exact`, so an `Approximate<T>` cannot be put there.
//! * Arithmetic: `F64` (floating point + tolerance policy) vs `Fp(p)` (exact arithmetic over a DIFFERENT field; for an
//!   integral δ, rank over F_p ≤ rank over R, so betti_p ≥ betti_R: a rigorous one-sided check, not "the answer").
use serde::{Serialize, Serializer};

/// Marker serialised as "exact".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Exact;
impl Serialize for Exact { fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> { s.serialize_str("exact") } }

/// Marker serialised as "approximate".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ApproximateTag;
impl Serialize for ApproximateTag { fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> { s.serialize_str("approximate") } }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Family { Exact, Approximate }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arithmetic { F64, Fp(u32) }
impl Serialize for Arithmetic {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self { Arithmetic::F64 => s.serialize_str("f64"), Arithmetic::Fp(p) => s.serialize_str(&format!("F_{p}")) }
    }
}

/// An exact-family result computed in f64 under a stated policy.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Numerical<T> { pub value: T, pub method: Exact, pub arithmetic: Arithmetic }
impl<T> Numerical<T> {
    pub fn f64(value: T) -> Self { Numerical { value, method: Exact, arithmetic: Arithmetic::F64 } }
}

/// An exact-family result in exact arithmetic over F_p.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ExactFp<T> { pub value: T, pub method: Exact, pub arithmetic: Arithmetic }
impl<T> ExactFp<T> {
    pub fn new(value: T, p: u32) -> Self { ExactFp { value, method: Exact, arithmetic: Arithmetic::Fp(p) } }
    pub fn p(&self) -> u32 { match self.arithmetic { Arithmetic::Fp(p) => p, Arithmetic::F64 => unreachable!("ExactFp is always F_p") } }
}

/// An approximate (ordered/interval layer) result: carries its assumptions; never convertible into the exact types.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Approximate<T> { pub value: T, pub method: ApproximateTag, pub assumptions: Vec<String>, pub warning: &'static str }
impl<T> Approximate<T> {
    pub fn new(value: T, assumptions: Vec<String>) -> Self {
        Approximate { value, method: ApproximateTag, assumptions, warning: "APPROXIMATION, not exact cohomology (ADR-013): does not work in all cases; conflicts are witnesses, not a cohomology group" }
    }
}
