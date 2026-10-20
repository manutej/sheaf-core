//! sheaf-cli: composes the workspace. `v1::analyze` (pure core) returns `(Report, Vec<Effect>)`; `src/main.rs` is the
//! shell (ADR-005). `m1` adapts the gated m1 contract (sheaf-contract) to sheaf-cell graph sheaves.
pub mod m1;
pub mod report;
pub mod v1;
