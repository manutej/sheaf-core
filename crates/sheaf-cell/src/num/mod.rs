//! Shared numerics: the H0 tolerance policy (`tol`), dense SVD / rank reports (`linalg`), CSR (`sparse`), CG
//! (`krylov`) and the family / arithmetic labels every result carries (`family`).
pub mod family;
pub mod krylov;
pub mod linalg;
pub mod sparse;
pub mod tol;

pub fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() }
pub fn norm(a: &[f64]) -> f64 { dot(a, a).sqrt() }
