//! sheaf-cli: composes the workspace. `v1::analyze` (pure core) returns `(Report, Vec<Effect>)`; `src/main.rs` is the
//! shell (ADR-005). `m1` adapts the gated m1 contract (sheaf-contract) to the Rust m1 math (sheaf_cell::m1 + sheaf_cohom::m1, ADR-016);
//! `ffi` is its C ABI for the Rust-backed m1 builder shim.
pub mod ffi;
pub mod m1;
pub mod report;
pub mod v1;
