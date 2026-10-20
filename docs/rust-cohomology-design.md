> sheaf-core note (2026-10-19): copied from FRONTIER `research/`. Since ADR-015, H⁰ and rank(δ⁰) use the single `H0Policy` (per-component cut), not the numpy-default rule described under "Field" below; rank(δᵏ), k ≥ 1, still uses it. Wording is "certified residual" (ADR-008).

# Real cellular-sheaf cohomology in Rust: survey, math, design, spike

Date: 2026-10-19 · Status: spike working (`sheaf-experiments/cohomology-spike/`, `cargo test --release` green)
Labels: **[V]** verified on this box (ran it / read the source at path:line) · **[C]** claimed by docs/authors, not run · **[S]** speculation/estimate.

## 0. Six-line summary
1. "Real cohomology" = ranks and kernels of the sheaf coboundaries δ0, δ1 on a 2-complex, with δ1δ0 = 0 checked, H^k = ker δk / im δk-1, harmonic reps = ker L_k. Over R it is numerical linear algebra plus a reported rank gap, never a proof. **[V]**
2. Nobody ships this in Rust. The Rust "sheaf" crates are graph-only or float Gaussian elimination with hard-coded 1e-12 pivots. lophat/oat_rust do real cohomology, but over F2/F_p/Q for constant coefficients only, with no restriction maps. CellularSheaves.jl is graph-only, H0 only. **[V]**
3. The sweet spot is a small typed core: `CellComplex` (signed incidences) → `Sheaf` (codim-1 maps) → `CochainComplex` (d²=0 verified on build, or a typed error) → `cohomology(k)`/`hodge(k, ω)`, each returning a certificate. Dense SVD for certification at ≤ ~5k dims, matrix-free CG/LSMR at scale, F_p rank as an exact cross-check. **[V spike]**
4. The spike passes 7 known-answer tests and 5 property tests (200 cases each). It matches numpy SVD on all 9 fixtures plus 3 numpy-from-scratch complexes, and matches CellularSheaves.jl on all 5 graph fixtures. It shows field dependence (RP² H¹ = 0 over R, 1 over F2) and rejects non-functorial maps with `NotAComplex`. **[V]**
5. Scale (box, 8 cores, one run each): the matrix-free graph Hodge split runs in 0.12 s at 100k V / 300k E. It runs in 0.42 s at ogbn-arxiv *size* (169,343 V / 1.34M E, synthetic random chords, not the real dataset), with certificate ‖δ0ᵀh‖/‖ω‖ ≈ 1e-10. Dense eigen is never needed for the split; it is needed only for explicit H¹ *bases*, which do not scale (b1 ≈ 1.17M there). **[V]**
6. Estimates: v1 graph-H¹ + Hodge split productised in about 3 days, 2-complexes in about 4 more, scale hardening in about 5 more. **[S]** Key finding: on a single cycle a holonomy that kills H0 also kills H1 (χ = dimC0 − dimC1 = 0). "Loop-level" signal therefore lives in the Hodge harmonic part of the *data* cochain, measured against the expected (identity) sheaf, or in graphs with χ < 0. **[V test]**

## 1. What exists (real, runnable code)

| Code | Language | Computes | Coefficients / field | Cells | Exactness / scale | Label |
|---|---|---|---|---|---|---|
| AlgebraicJulia **CellularSheaves.jl** (registry, installed here; `src/network_sheaves/EuclideanSheaves.jl`) | Julia | coboundary (`coboundary_map` L254), Laplacian, **H0 via sparse chordal LDLᵀ** (`nullspace_ldlt` L559), harmonic extension, nearest global section (CG/LDL/pinv L469-532) | R, arbitrary per-edge restriction maps | **graphs only**, no faces, no H1 API | numerical pivot threshold; sparse | [V] ran `crosscheck.jl`: 5/5 match |
| **pysheaf** (kb1dds; `pysheaf/pysheaf.py`) | Python | consistency radius (`ComputeConsistency` L159), filtrations, `FuseAssignment` L618 (optimisation) | arbitrary (possibly nonlinear) maps on a poset | cell poset | no cohomology at all (grep: zero hits for cohomology/coboundary) | [V] source |
| **neural-sheaf-diffusion** (twitter-research; `lib/laplace.py`) | PyTorch | sparse sheaf Laplacian L0 (`build_sheaf_laplacian` L27), normalised variants, diffusion | R, learned maps | graphs | L0 only; no H1 | [V] source |
| **gluing-lc** (`sheafint/core/hypergraph.py:171` `algebraic_cohomology`; `experiments/utils/hodge.py:216` `hodge_2d`) | Python/numpy | H0/H1/H2 of the **constant** sheaf via `matrix_rank` (correct formula H1 = (nE − rk d1) − rk d0, L190); Hodge split by `lstsq` on (nE, k) arrays, Kronecker-factored | R | 2-complexes (constant sheaf only); the general-sheaf `ScalableSheaf` path is rank-capped (see gluing-lc review) | dense; `lstsq`, not pinv (comment at L219) | [V] source + earlier replication |
| crates.io **lau-sheaf-cohomology 0.1.0 / sheaf-cohomology 0.3.0** (near-identical source, 2-line diff) | Rust | rank/kernel by float Gaussian elimination, pivot `> 1e-12` hard-coded (`src/lib.rs:189,224`) | R (f64) | claims cellular | no tolerance policy and no gap report; about 25–63 downloads | [V] source |
| crates.io **sheaf-spectral 0.1.0** | Rust (nalgebra DMatrix) | dense coboundary (`sheaf.rs:88`), synchronization, "persistent H0" | R | graphs | dense | [V] source, not run |
| **lophat 0.14.1** (Chaplin) | Rust | persistent homology of arbitrary filtered chain complexes, R = DV, clearing, anti-transpose for cohomology | **F2 only** (column = set of indices, `columns/vec.rs`) | any chain complex you supply | exact | [V] source, not run |
| **oat_rust 0.2.0** (Open Applied Topology) | Rust | chain complexes, boundary matrices, (persistent) homology, U-match | BooleanField, `PrimeOrderField`, native types incl. `num::Ratio` (Q) (`algebra/rings/types/`) | simplicial (VR etc.) | exact | [V] source, not run |
| GUDHI / PHAT / Ripser | C++ | persistent (co)homology | F_p | simplicial/cubical | exact, very scalable | [C] not checked here |
| **faer 0.24.4** | Rust | dense SVD/EVD/col-pivoted QR (`linalg/solvers.rs:364-457`); sparse Cholesky (LLᵀ/LDLᵀ/LBLᵀ), LU, **QR without pivoting** (`sparse/linalg/qr.rs` uses `qr::no_pivoting`); matrix-free **CG, LSMR, BiCGSTAB, partial self-adjoint eigen, partial SVD** (`src/operator/*`) | R/C | n/a | **no sparse rank-revealing QR**; compiles on rustc 1.85.1 | [V] source + compiled |

Bottom line **[V]**: nothing in Rust or Julia computes H¹ (or harmonic reps, or the Hodge split) of a sheaf with non-trivial restriction maps on a 2-complex with a stated tolerance. The math is standard, the building blocks (faer) exist, and the field-exact tools (oat_rust/lophat) cover only constant coefficients.

## 2. The math (precise)

**Complex.** X is a finite regular CW complex of dimension ≤ 2 (spike: simplicial). Each cell carries a fixed orientation, and incidence numbers [σ:τ] ∈ {−1, 0, +1} for dim τ = dim σ + 1 satisfy Σ_τ [σ:τ][τ:ρ] = 0 (∂² = 0). For a simplicial complex with vertex order, ∂[a,b] = b − a and ∂[a,b,c] = [b,c] − [a,c] + [a,b].

**Sheaf.** F assigns a stalk F(σ) = R^{d_σ} to each cell and a linear map F_{σ⊴τ}: F(σ) → F(τ) to each codim-1 incidence. Functoriality: for every diamond σ ⊴ τ, τ' ⊴ ρ (dim ρ = dim σ + 2), F_{τ⊴ρ}F_{σ⊴τ} = F_{τ'⊴ρ}F_{σ⊴τ'}. On a regular complex it is equivalent to δ1δ0 = 0 for the assembled coboundary, which is what we check.

**Cochains.** C^k = ⊕_{dim σ = k} F(σ). δ_k: C^k → C^{k+1} has block (τ, σ) = [σ:τ] F_{σ⊴τ}. Then δ_{k+1}δ_k = 0.

**Cohomology.** H^k = ker δ_k / im δ_{k−1}, so dim H^k = dim C^k − rank δ_k − rank δ_{k−1}. Euler: Σ(−1)^k dim C^k = Σ(−1)^k dim H^k.

**Hodge.** L_k = δ_kᵀδ_k + δ_{k−1}δ_{k−1}ᵀ (for k = 0 this is the sheaf Laplacian L0 = δ0ᵀδ0). ker L_k = ker δ_k ∩ ker δ_{k−1}ᵀ ≅ H^k (harmonic representatives, w.r.t. the standard inner products). Every ω ∈ C^k splits orthogonally:
ω = δ_{k−1}x (exact, "explained by a vertex assignment") + h (harmonic, "loop-level / global obstruction") + δ_kᵀy (coexact, "explained by face-level curl"), with ‖ω‖² = sum of the three energies.

**Graph-only complex.** C² = 0, so H¹ = coker δ0, dim H¹ = dim C¹ − rank δ0, and ω = exact + harmonic. The harmonic part is the edge disagreement that *no* vertex assignment can explain. This is exactly the "explained vs loop-level" split gluing-lc's `hodge_2d` computes, and what product case 4 needs.

**Field.** Over **R** (default): rank = #{σ_i > tol}, tol = max(m,n)·ε·σ_max (the numpy/LAPACK default). Each rank decision is reported with σ_kept_min / σ_dropped_max (the gap). It gives orthogonal Hodge projections, works with real-valued maps and scales with Krylov methods, but rank is ill-posed near the threshold. Over **Q / F_p** (exact): ranks are exact, but you need integral or rational maps (fitted rotations aren't), there are no orthogonal projections (no Hodge), coefficients grow over Q, and **the answer can differ from R**. Verified: RP² has H¹ = 0 over R and F3 but 1 over F2 (torsion). Use F_p only for certifying integral fixtures, or for constant/sign sheaves (e.g. signed-graph balance).

## 3. Assumptions (stated, each enforced or reported)
- **A1** Finite regular cell complex, dim ≤ 2 in the spike (the code is k-generic in `inc`). *Enforced by builder.*
- **A2** Fixed cell orientations, with incidence signs supplied by the builder. H^k dims are orientation-independent; signs of components are not. *By construction.*
- **A3** ∂² = 0 on the complex. *Checked* (`CellComplex::validate`, error `BoundaryOfBoundaryNonzero`).
- **A4** Finite-dimensional real stalks; stalk dims may vary per cell. *Shape-checked* (`SheafError::MapShape`).
- **A5** Only codim-1 maps are stored. Functoriality is **not assumed**: δ_{k+1}δ_k = 0 is verified with relative Frobenius residual ≤ 1e-10, else `CohomError::NotAComplex`. A twisted map on a filled face is rejected (test). **Consequence:** a non-trivial holonomy around a *filled* face is not a sheaf. Holonomy lives on unfilled loops.
- **A6** Numerical rank with the tolerance above, and every rank comes with its gap. **No result is called a proof.** The dimension is computed two independent ways (rank-nullity vs the nullspace of [δ_k; δ_{k−1}ᵀ]), and disagreement is surfaced (`Cohomology::consistent`).
- **A7** F_p results are only computed when every coboundary entry is an integer (`fp::to_int` returns None otherwise).
- **A8** Inner products are the standard Euclidean ones on each stalk. Harmonic reps and the Hodge split depend on this choice of weighting (dims do not). Weighted/normalised variants (NSD-style D^{-1/2}) are a v2 option.

## 4. Crate API (spike = implemented; v1 = proposed delta)

```rust
// pure core, no IO (ADR-005): everything returns values + certificates; effects (logs, files) are the caller's.
pub struct CellComplex { n_cells: Vec<usize>, inc: Vec<Vec<Incidence{facet, cofacet, sign}>>, verts }   // [V]
impl CellComplex { fn simplicial(n, edges, faces) -> Self; fn validate(&self) -> Result<(), ComplexError> }
pub struct Sheaf { complex, stalk: Vec<Vec<usize>>, maps: Vec<Vec<Mat<f64>>> }                         // [V]
impl Sheaf { fn constant(cx, d); fn validate(); fn coboundary_triplets(k) -> (rows, cols, Vec<(i,j,v)>) }
pub struct CochainComplex { delta: Vec<Mat<f64>>, dims, d2_residual }                                 // [V]
impl CochainComplex {
  fn from_sheaf(&Sheaf) -> Result<Self, CohomError /* Sheaf(..) | NotAComplex{k, rel_residual} */>;
  fn cohomology(&self, k) -> Cohomology { dim, rank_in: RankReport, rank_out: RankReport, harmonic_basis, laplacian_rank };
  fn hodge(&self, k, &[f64]) -> Hodge { exact, harmonic, coexact, potential, cert: HodgeCertificate };
}
RankReport { rank, tol, sigma_max, sigma_kept_min, sigma_dropped_max } .gap_ratio()                  // [V]
HodgeCertificate { reconstruction, ortho_*, harmonic_closed, harmonic_coclosed }  // all relative     // [V]
fp::betti_mod_p(&delta, &dims, p) -> Option<Vec<usize>>                                               // [V]
sparse::graph_hodge_cg(&Csr, &w, rtol, max_iter) -> GraphHodge { exact, harmonic, potential, iters, cert_coclosed } // [V]
// v1 additions [S]: Backend enum {DenseSvd, Krylov{rtol}, ExactFp(p)}; hodge on sparse 2-complexes via two
// LSMR solves (min ||δ0 x − ω||, then min ||δ1ᵀ y − (ω − δ0x)||); per-edge harmonic energy ranking (the audit's
// "loop-level discord"); random-cochain null b1/dimC1 reported next to every harmonic fraction; serde JSON in/out
// matching sheaf-experiments/contract/sheaf.schema.json, extended with `faces` and per-incidence maps.
```

## 5. Numerical strategy at scale
- **Hodge split (the product need): never form a basis.** Exact part = δ0 x*, with x* solving L0 x = δ0ᵀω by CG, which is consistent because the RHS ∈ range(L0). Coexact part = a second solve with δ1 (LSMR on δ1ᵀ y ≈ residual). Certificate: ‖δ0ᵀh‖/‖ω‖ and ‖δ1h‖/‖ω‖. **[V]** for the graph case (spike `sparse.rs` uses its own CSR plus CG; faer's `operator::conjugate_gradient`/`lsmr` exist **[V source]** and should replace it in v1, with an AMD/sparse-Cholesky preconditioner if iterations grow **[S]**).
- **dim H^k at scale:** H0 = nullity(L0). Use LDLᵀ zero-pivots (what CellularSheaves.jl does; reliable for graph Laplacians, heuristic in general **[C]**), or faer `partial_self_adjoint_eigen` for the smallest eigenvalues with the spectral gap reported. Then dim H¹ (graph) = dim C¹ − dim C⁰ + dim H⁰: **one** nullity gives H¹ for free. On 2-complexes you also need rank δ1. Use the Euler characteristic plus H⁰ and H² (the nullity of δ1δ1ᵀ on C², often small).
- **Explicit H¹ bases** do not scale (b1 ≈ 1.17M at ogbn-arxiv size). Return only the top-m harmonic directions that carry the data's energy (project h onto a Krylov subspace) **[S]**.
- **Certification:** dense SVD for dims ≤ ~5k, plus F_p rank for integral fixtures. faer has **no sparse rank-revealing QR [V]**. SuiteSparseQR (C, via FFI) would be the only sparse RRQR option and is out of scope **[S]**.
- Measured **[V]**: `cargo run --release --example scale` on 100k V / 300k E (d=1) gives build 0.14 s, CG 0.12 s, 56 iterations, cert 1.6e-10. At 169,343 V / 1,335,510 E, build 0.83 s, CG 0.42 s, 36 iterations, cert 8.9e-11. Harmonic fraction 0.889 vs random null b1/E = 0.873. This is a synthetic random graph (expander, so CG converges fast). Real citation graphs with long chains will need more iterations **[S]**.

## 6. Tests (all [V], `cargo test --release`)
- Known answers over R: triangle graph (1,1) vs filled (1,0,0); 5-cycle R² (2,2); S² (1,0,1); 7-vertex torus (1,2,1); RP² (1,0,0); triangle with 60° rotation holonomy (0,0); reflection holonomy (1,1); K4 with rotation (0,4). Euler characteristic checked for each.
- F_p: RP² over F2 = (1,1,1), over F3 = (1,0,0); torus over F2 = (1,2,1).
- `NotAComplex` on a rotation twist inside a filled face.
- Translation loop: 3-cycle, +0.1 drift per hop. Exact part = 0, harmonic energy = 0.03 (= 0.3²/3), and pure δ0x data has 0 harmonic; the potential recovers x up to H⁰.
- Holonomy: K4 with rotation has H⁰ = 0, the constant assignment shows discord, and a random cochain has a harmonic part.
- Sparse CG harmonic matches dense SVD (50-cycle, 1e-8).
- proptest (200 cases each): d² = 0 on random graph sheaves with random maps; rank-nullity + Euler + two-way dim agreement; Hodge orthogonality, reconstruction, Pythagoras, and harmonic ∈ span(basis); gauge-twisted sheaves F_{σ⊴τ} = A_τA_σ⁻¹ on random 2-complexes keep d² = 0 and the constant sheaf's Betti numbers.
- proptest (Grok r1 #5): generic random maps on filled-face complexes (triangle, S²) are always rejected as `NotAComplex`.
- `max_iter = 0` CG is reported `converged = false`.
- Cross-checks: `python3 crosscheck.py` gives ALL_OK (9 fixtures re-ranked by numpy from Rust triplets, plus torus/RP²/S² built from scratch in numpy). `julia crosscheck.jl` gives ALL_OK (5 graph fixtures rebuilt natively in CellularSheaves.jl; H0 by `nullspace_ldlt` equals H0 by SVD).

## 7a. Grok adversarial review (round 1 of max 2; `grok` CLI headless, raw text in `cohomology-spike/out/grok_review_r1.md`)
Two attempts that inlined the whole design and spike (~50 KB) timed out with no output (420 s and 760 s). A trimmed prompt (`cohomology.rs`, `sparse.rs`, `props.rs`, §2–3 of this doc) returned 10 findings, **none CRITICAL** ("dense rank-nullity and Hodge formulas are right when δ² = 0"). Severity rule applied: no repro, no critical.
| # | Grok sev | Finding | Disposition |
|---|---|---|---|
| 1 | HIGH | CG has no curvature guard on singular L0 | **Fixed**: breakdown if pᵀLp ≤ 1e-14‖p‖² (`sparse.rs`) |
| 2 | HIGH | Unconverged CG returned as a split; stop rule ≠ certificate | **Fixed**: `converged` is now judged on the certificate ‖δ0ᵀh‖/‖ω‖ ≤ cert_tol, plus a `breakdown` flag; test for max_iter = 0 |
| 3 | HIGH | `n − rank − rank` usize can wrap | **Fixed**: `cohomology()` returns `Result`, with `CohomError::RankConflict` |
| 4 | HIGH | rank-nullity/Euler props are identities | **Accepted, partly**: annotated as sanity checks. The non-trivial checks are the two-way dim agreement and the golden values in `known.rs` (Grok was not shown that file) |
| 5 | HIGH | d² prop never stresses faces with non-constant maps | **Fixed**: new prop, random maps on triangle/S² must give `NotAComplex`. The gauge prop already covers valid random maps on faces |
| 6 | HIGH | gauge-vs-constant can't detect dropped maps | **Already covered** by golden non-isomorphic cases in `known.rs` (rot60 → (0,0), K4 → (0,4), reflection → (1,1)) |
| 7 | MED | d² tolerance (rel. Frobenius 1e-10) not tied to the rank cut | **Open, documented**: the ⟨exact, coexact⟩ certificate catches the consequence. v1: per-degree residuals tied to σ-tolerance |
| 8 | MED | `reconstruction` is zero by construction; 3 ortho checks ≈ 1 number | **Fixed**: replaced with `potential_residual` = ‖δ_{k−1}x − exact‖/‖ω‖ (an independent identity) |
| 9 | MED | CSR doesn't merge duplicate triplets | **Rejected**: CSR matvec sums duplicates, which is the same semantics as the dense `+=` assembly. No defect |
| 10 | MED | potential uses 1/σ down to the cutoff | **Fixed (reported)**: `Hodge::potential_cond` = σ_max/σ_kept_min |
Round 2 was not run (time budget); the fixes above are covered by the re-run test suite and cross-checks.

## 7. Verified vs speculation
Verified: everything in §6, the scale timings in §5, the survey rows marked [V]. Claimed: GUDHI/PHAT capabilities, LDLᵀ-nullity reliability. Speculation: effort estimates, real ogbn-arxiv iteration counts, preconditioner need, the API delta for v1.

## 8. Effort (days, one engineer + agents) [S]
| Milestone | Scope | Days |
|---|---|---|
| v1 graph-H¹ + Hodge split | port to faer `operator::conjugate_gradient`; per-edge harmonic energy ranking; random-cochain null; JSON in/out on the contract schema with per-incidence maps (fixes the `from_spec` map-per-type gap); Python parity ≤1e-9 vs sheafcore; wire into audit-local | 3 |
| 2-complexes | sparse δ1, LSMR coexact solve, face specification in schema (which loops are "filled" = declared-consistent cycles), H² via δ1 nullity, F_p certifier on integral fixtures | 4 |
| Scale (ogbn-arxiv real data) | streaming CSR build, preconditioner (incomplete Cholesky / AMG via FFI if needed), partial eigen for H0, memory budget (≈ (dimC1·nnz-per-row)·12 B), benchmark harness | 5 |

## 9. How it plugs in
- **product/audit-local:** facts become vertices, "same key across files" relations become edges, and unit-conversion maps become restriction maps. Today we rank `edge_discord`. Add `harmonic energy per edge` = the part of each conflict that no single "true value" assignment can reconcile (a cycle of mutually inconsistent conversions), separate from plain one-file-is-wrong conflicts (exact part, whose potential x points at the outlier vertex). The random-cochain null b1/dimC1 calibrates the score.
- **Translation-loop case (product/cases 04):** use the *expected* sheaf (identity round trip) and put observed per-hop drift on edges. Then the harmonic part = the loop circulation spread over the loop. It is zero for any drift a per-language offset explains. Alternatively, put the fitted per-hop maps on the edges and read H⁰: a holonomy without eigenvalue 1 kills H⁰. **But on a single cycle it also kills H¹ (χ = 0)**, so "nonzero loop-level component" must come from the data-cochain view, or from graphs with more edges than vertices (K4 test: H⁰ = 0, H¹ = 4).
- **Work Graph:** Task Board ↔ Obligations duplicates form 2-cycles and longer loops. The harmonic part flags inconsistency that propagates around a loop (CCAF-like chains) rather than one stale row.
