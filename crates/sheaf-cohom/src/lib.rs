//! sheaf-cohom: real (over R) cellular sheaf cohomology. H^k with harmonic bases (`cohomology`), relative
//! H^0(X,B;F) (`relative`: the uniqueness certificate for harmonic extension), the Hodge split and the nearest global
//! section (`cohomology::hodge`, `CochainComplex::nearest_section`), the sparse graph Hodge split (`hodge_graph`),
//! scalable spectral H^0 on faer (`spectral`, incl. the per-component driver), and exact Betti numbers over F_p (`fp`).
//! Depends on sheaf-cell only. Reported with certified residuals and rank gaps (never "certified rank", ADR-008).
//! Everything here is family `Exact` (ADR-013); no ordered/interval approximation layer lives here.
//! The m1 contract (sheaf-contract) is bridged in sheaf-cli (`tests/contract_parity.rs`), not here.
pub mod cohomology;
pub mod fp;
pub mod hodge_graph;
pub mod relative;
pub mod spectral;
