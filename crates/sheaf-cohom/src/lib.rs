//! sheaf-cohom: real (over R) cellular sheaf cohomology. H^k with harmonic bases (`cohomology`), relative
//! H^0(X,B;F) (`relative`: the uniqueness certificate for harmonic extension), the Hodge split and the nearest global
//! section (`cohomology::hodge`, `CochainComplex::nearest_section`), the sparse graph Hodge split (`hodge_graph`),
//! scalable spectral H^0 on faer (`spectral`, incl. the per-component driver), and exact Betti numbers over F_p (`fp`).
//! Depends on sheaf-cell only. Reported with certified residuals and rank gaps (never "certified rank", ADR-008).
//! Everything here is family `Exact` (ADR-013); no ordered/interval approximation layer lives here.
//! `m1` is the m1 contract's own H⁰ (R-H0-08 global cut, ADR-016); the contract types (sheaf-contract) are bridged in
//! sheaf-cli, not here.
pub mod cohomology;
pub mod fp;
pub mod hodge_graph;
pub mod m1;
pub mod relative;
pub mod spectral;
