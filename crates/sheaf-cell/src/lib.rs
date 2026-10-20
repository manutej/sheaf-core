//! sheaf-cell: the objects. Cell complex / face poset (`complex`: `graph()` keeps order, multi-edges and orientation;
//! `simplicial()` sorts and dedups), cellular sheaves with per-incidence restriction maps (`sheaf`), coboundaries
//! δᵏ, matrix-free Laplacians L_k, energy and discord (`ops`), graph sheaves with per-edge maps (`graph`), and the
//! shared numerics `num` (H0Policy tolerance type, dense SVD/rank, CSR, CG) plus the policy-carrying result types
//! (`num::family`: Exact vs Approximate, f64 vs F_p). Bottom of the DAG: depends on faer only (no other sheaf crate).
//! Convention everywhere: (δx)_e = F_v x_v − F_u x_u, L₀ = δ⁰ᵀδ⁰ (SPEC §1, Hansen–Ghrist).
pub mod complex;
pub mod fixtures;
pub mod graph;
pub mod m1;
pub mod num;
pub mod ops;
pub mod sheaf;
