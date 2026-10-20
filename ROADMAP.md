# Roadmap

1. **Port base to faer (top TODO).** The m1 Rust crate (`sheaf-contract`) does not use faer; m1's δ⁰/L₀/energy/H⁰ are
   Python-only. Port them to Rust on faer *through the gate* so the gated behaviour stays identical (R-RUST-20 parity,
   mutation table, hidden holdout). Not done in the repo setup on purpose.
2. ~~Layout decision~~ **Done (ADR-014 accepted):** sheaf-cell / sheaf-cohom / sheaf-solve / sheaf-contract / sheaf-cli,
   harmonic extension in sheaf-solve, relative H⁰(X,B;F) in sheaf-cohom, `sheaf analyze --boundary`.
3. **Harmonic extension at scale:** relative H⁰ is per-block dense (refuses blocks > `dense_max`); add a sparse
   (shift-invert on L_II) path, and an x0 option in the CLI.
4. **Formatting/lint gate:** the box toolchain has no rustfmt/clippy, so CI runs `cargo fmt --check` and `clippy`
   as advisory steps. Format sheaf-cell/cohom/solve/cli once with rustfmt, then make both blocking (sheaf-contract stays
   byte-identical to the gate and is excluded from formatting).
5. **Review flags** (F4 done: `num::family` types; F12 partly: `krylov::SolveTol`) F5 (`ValidSheaf` newtype), F8 (two named spec parsers,
   m1 vs v1; R-NONID-04 vs identity default), F9 (one sparse `CochainComplex`; faer `SparseRowMat`), F11 (library
   `expect`s → `Result`), F13 (rayon feature for wasm32).
   sheaf-core uses its own CG (`sheaf_cell::num::krylov`); the graph Hodge CG and the shift-invert inner PCG still have
   their own loops (dedup pending). If faer's matrix-free CG is adopted: use `InitialGuessStatus::MaybeNonZero` (the `Zero` path misreports in 0.24.4);
   pad `lsmr_scratch`; `try_new_from_triplets` sums duplicates (multi-edges must keep distinct rows, as `graph()` does).
6. **Faces in v1** (2-complexes in the CLI), exact eigenvalue count (inertia) to certify Krylov completeness,
   FB15k rerun on the final code, speed tuning.
7. **Ordered/interval layer** stays out of this repo until it has a caller; when it lands it is labelled
   approximate and never writes into exact fields (ADR-013).
