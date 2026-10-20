# sheaf-core decisions

ADR-004, -005, -008, -010, -012 and -013 are copied verbatim from FRONTIER `product/DECISIONS.md` (2026-10-19) for context;
FRONTIER stays the record for those. ADR-014 and ADR-015 are new here.

## ADR-004 — Stack for v0
- **Date:** 2026-10-19 03:55 Bahia
- **Decision:**
  - **Frontend = Leptos (ACCEPTED, Manu 03:35 Bahia: "Let's go with leptos").** Rust, CSR/static build via trunk or cargo-leptos, hydration later. Visuals: three.js WebGPURenderer + TSL and cosmos.gl through thin wasm-bindgen/JS interop. Rust types are shared end to end between core and UI.
  - **Rest (PROPOSED):**
    - Core: Rust `frontier-core` on faer, pure with no I/O, built to WASM (wasm-bindgen) and to a PyO3/maturin wheel.
    - Contract: Rust `serde` + `schemars` types are the single source of truth, and JSON Schema, TS and Python types are generated from them.
    - Deploy: Cloudflare Pages static, mirrored to a Hugging Face Static Space.
    - Temporal: none in v0. Temporal (Rust SDK) is planned for when the builder→verifier→reviewer gate is hosted as a service.
    - Review cap: at most 3 adversarial rounds per artefact, diff-only re-review.
- **Options:** SvelteKit-static or Astro (stack research recommendation, superseded for the frontend by Manu); React (excluded, Manu); Restate/DBOS instead of Temporal (re-evaluate when hosting the gate).
- **Rationale:** `product/research/stack-2026-10.md` §4. No backend in v0, so there are zero ops and the demo runs the real core in the browser. Temporal is overkill until there are long-running hosted workflows; our 50-minute gate cutoffs are its first real use case. The 3-round cap applies the lesson of m1's 18 rounds.
- **Risks:** Leptos has lower adoption and a smaller ecosystem, and coding agents default to React. Mitigations: keep JS interop thin and limited to the visuals, share types end to end, and pin the Leptos docs in AGENTS.md. faer SIMD on wasm32 is unproven.
- **Spikes:** S1 faer → wasm32 size and speed, plus Leptos + core WASM bundle size. S2 three.js WebGPU/TSL and cosmos.gl interop from Leptos.
- **Status:** Frontend ACCEPTED. Everything else PROPOSED.

## ADR-005 — Functional core, effects as data
- **Date:** 2026-10-19 03:55 Bahia
- **Decision:** All code inside scope is pure and never performs side effects. It returns a typed `Effect` enum (EmitReport, WriteFile, HttpPost, GitPush, NotionWrite). A single audited interpreter (the imperative shell) validates each effect against the contract, the caller's `Capability` and policy, executes it, appends a content-hashed, hash-chained `EventRecord`, and feeds the `Outcome` back as input. State is a fold over the log, so it can be replayed. Capabilities are typed, scoped, expiring and MAC'd by the interpreter, so pure code cannot forge them, and they map onto the builder, evaluator and reviewer lanes. Later, Temporal workflows are deterministic folds and effects become activities; Leptos uses pure signals with effects only via explicit actions. Keep it lean: one enum, one interpreter, no monad stacks. Types are in PRODUCT_SPEC §7.6.
- **Options:** (a) ad-hoc I/O in the core; (b) effect systems or monad stacks (rejected as heavy); (c) effects as plain data plus one interpreter (chosen).
- **Rationale:** Manu, 3:49 AM Bahia: all externally modifying events from inside scope should involve functional programming. It also matches the existing lane separation and the approval rule for external actions. The interpreter is the one place that needs auditing.
- **Status:** PROPOSED (reviewed in adversarial round 2 / focused pass).


## ADR-008 — Wording: "certified residual", never "verified" or "proof"
- **Date:** 2026-10-19 03:55 Bahia · **Status:** PROPOSED (default)
- **Decision:** The report's numeric block (h0_dim, energy, kernel residual ≤ τ) is called the **certified residual**. It certifies that the consistency computation is correct, not that a flag is true. It is not a headline badge in the demo.
- **Rationale:** Round 1, product-skeptic-4: the badge read as proof of a contradiction. META_PLAN Q9.


## ADR-010 — Repo stays private
- **Date:** 2026-10-19 03:55 Bahia · **Status:** PROPOSED (default)
- **Decision:** manutej/frontier-sheaf stays private until Manu says otherwise; the hosted demo publishes only the static build and the public corpus. No license is chosen yet.


## ADR-012 — v1 spectral method and tolerance policy (2026-10-19 Bahia) [ACCEPTED for v1, spike-level]

- **Context:** Manu asked whether the spectral eigenvalue decomposition was nailed. Before this, it was dense-only. v1 also needed per-edge maps, which the per-type `from_spec` couldn't express.
- **Decision: spectral method.** The smallest-k eigenpairs of L0 (and of a graph's L1) are computed per connected component, since L0 is block diagonal.
  - A component with ≤ 1500 dofs gets a dense EVD of all eigenvalues.
  - A larger component gets faer 0.24's thick-restart Lanczos (`matrix_free::eigen::partial_self_adjoint_eigen`). It only returns the largest-magnitude eigenpairs, so it runs on the shift-invert transform `P (L + σI)⁻¹ P`, with σ = 1e-6·‖L‖_bound and an inner Jacobi-PCG solve to a relative residual of 1e-12.
  - `P` projects out the exact known kernel (component constants when every edge has `Fu == Fv` injective) and the locked vectors. Lock and restart until a fresh run adds nothing at or below the current (k+1)-th value, then one Rayleigh–Ritz pass on L itself.
  - For a graph, L1's nonzero spectrum is L0's. It's checked by mapping `u = δv/√λ` and measuring `‖L1u − λu‖` directly. L1's kernel dim comes from rank-nullity.
- **Options considered:**
  - (a) Polynomial transform `cI − L` with no inner solves. Kept as `--method polynomial`. Measured too slow on path-like graphs: a 600-dof cycle-plus-chords sheaf gave no result after about 2.5 minutes and was killed. Shift-invert does the same case in 0.17 s.
  - (b) LOBPCG. faer has none, and the round-3 scipy LOBPCG ran out of memory on ogbn-arxiv.
  - (c) Shift-invert with an exact factorization. There's no sparse LDLᵀ in our stack, and it's memory-heavy at 1M edges.
  - (d) Shift-invert with PCG. Chosen.
- **Decision: tolerance policy.**
  - An eigenvalue λ counts as zero iff λ ≤ 1e-10·max(1, ‖L0_c‖_bound), using the bound of its own component c.
  - `dim_global_rule` reports the sheafcore R-H0-08 count, which uses the global λ_max, for parity.
  - `dim_strict` counts λ ≤ max(64ε‖L0_c‖, 10·residual) per component.
  - `near_zero_ambiguous` is set when the three counts disagree.
  - The per-component rule replaced the global one after Grok v1 finding #1. Repro: components at scales 1e5 and 1e-1 give a true H0 of 2, but the global rule says 3. That case is now a test.
  - Each eigenvalue carries its residual r. A true eigenvalue lies in [λ − r, λ + r] (Weyl). The report calls this "certified residual" (ADR-008). The zero/nonzero labelling is a policy, and Krylov completeness (no missed eigenvalue below the cut) is **not** certified. The H0 field is `cut_consistent_with_zero_tol`, never "proof".
- **Evidence (box runs, Oct 19, single run each):**
  - Lanczos vs dense: max relative eigenvalue error 9.9e-16 over 200 seeded random graph sheaves (k+1 = 7). Fixtures: 6.7e-16, H0 18/18.
  - One H0 disagreement out of 200 is a policy split, not a solver error: λ = 1.8e-10 counts as zero by the eigen rule, but its σ = 1.36e-5 is nonzero by the spike's SVD rule. `dim_strict` agrees with the SVD rule on 200/200.
  - Python sheafcore parity: worst relative error 4.7e-16 on 5 fixtures (tolerance 1e-9).
  - Scale timings are in `sheaf-experiments/cohomology-spike/V1.md` and `out/bench/*.json`.
- **Status:** ACCEPTED for v1. Grok review: 1 round, 10 findings, 1 critical with a repro (fixed). Open items: widening k inside `smallest` itself (only the component driver widens), deflating inside the inner PCG, and an inertia (Sylvester) count to certify completeness.

## ADR-013: Ordered/interval constraint layer is labelled an approximation (2026-10-19)
- Status: accepted (Manu, 2026-10-19: "it doesn't work in all cases, flag it for awareness as something that is an approximation, and not the same as exact case").
- Decision: the ordered-sheaf / interval constraint layer (research/ordered-sheaves-note.md) is a heuristic approximation, kept strictly separate from the exact cohomology core. Its outputs carry `"method": "approximate"` and a warning line in every report; they never feed H0/H1/H2 or Hodge figures, and are never described as proof or as exact cohomology.
- Known limits: obstructions are witnesses, not a group; propagation misses unbounded infeasible loops (needs a negative-cycle check); results that depend on assumed rules must name the assumption.

## ADR-014 — Repo layout: one repo, five crates, not microservices (2026-10-19 Bahia) [ACCEPTED — Manu, 2026-10-19]
- **Context:** Manu (2026-10-19): ship the sheaf code as its own private repo `sheaf-core`, Rust-only on faer with Python as backup; asked whether m1 and the cohomology work are "both sheaf-core" (yes) and whether they should be microservices (no). The first three-crate proposal (`sheaf-base`, `sheaf-cohom`, `sheaf-cli`) went to a category-theory architecture review (FRONTIER `research/sheaf-core-architecture-review.md`); Manu approved its layout, with `sheaf-base` renamed `sheaf-contract`.
- **Decision:** one Cargo workspace, five crates, dependency DAG `sheaf-cell → {sheaf-cohom, sheaf-solve} → sheaf-cli` and `sheaf-contract → sheaf-cli`, no other edges, no cycles (CI checks the forbidden edges with `cargo tree`).
  - `sheaf-cell` (bottom): cell complex / face poset (`graph()` order/multi-edge/orientation-preserving, `simplicial()`), sheaf validation, cochains, δᵏ, L_k, energy/discord, graph sheaves (the single δ⁰ assembly), and `num`: `H0Policy` (ADR-015), dense SVD, CSR, CG, and the policy types (`Exact` vs `Approximate` family; `Numerical<T>` f64 vs `ExactFp<T>` F_p).
  - `sheaf-cohom`: Hᵏ, relative H⁰(X,B;F) (uniqueness certificate for harmonic extension, Hansen–Ghrist Prop 4.1), Hodge split, nearest global section, spectral (incl. the per-component H⁰ driver), mod-p. Depends on cell only.
  - `sheaf-solve`: harmonic extension (moved from `pending/sheaf-harmonic`). Depends on cell only, never on cohom: whether an extension is unique is a cohomology question, composed by the caller.
  - `sheaf-contract`: the m1 core, renamed from `sheaf-base`; `src/lib.rs` byte-identical to the gate (sha256 `7128f2b4…6de`, CI-checked), no faer, changes only through the gate.
  - `sheaf-cli`: composes everything: v1 parser/analyze (typed `Report`), the m1 adapter, `--boundary` (extension + `h0_relative`/`unique`), the `sheaf` binary.
- **Not microservices:** the crates share sheaf types and the coboundary; a network boundary would copy large matrices and duplicate the core math. A service wraps one crate later, only when something needs it.
- **Rust-only on faer.** Python (`python/`) is a reference for parity checks only and is never a runtime dependency.
- **Binary name:** `sheaf` (`sheaf analyze <spec.json>`). Wire formats keep the ids `frontier-sheaf/v1` / `frontier-sheaf/report-v1` (plus `frontier-sheaf/boundary-v1`).
- **Fact (recorded):** the m1 Rust crate does **not** use faer. It is a serde mirror + `validate()` + per-type `discord()`. m1's δ⁰/L₀/energy/H⁰ exist only in the Python reference `sheafcore.py`; "18-round validated" applies to that contract, not to a Rust δ⁰/L₀/H⁰.

## ADR-015 — One H⁰ tolerance rule, typed (2026-10-19 Bahia) [ACCEPTED in sheaf-core; amends ADR-012's tolerance policy]
- **Context:** the architecture review (F2) found three H⁰ rules in use: SPEC R-H0-08's global cut (`dim_global_rule`), the spike's numpy-default SVD rank for δ⁰ in `cohomology(0)`, and ADR-012's per-component cut (plus a per-component `dim_strict` count). They could disagree at the margin.
- **Decision:** `sheaf_cell::num::tol::H0Policy` is the only H⁰ rule. Default `ZeroCut::PerComponent { rel: 1e-10 }`: λ counts as zero iff λ ≤ 1e-10·max(1, ‖L₀|_c‖_bound) on its own connected component c (ADR-012). `ZeroCut::Absolute(t)` is the per-context permitted delta Manu asked for (v1 option `h0_absolute_tol`). Used by `cohomology(0)`, rank(δ⁰) inside `cohomology(1)` and `hodge(0|1)`, `spectral::smallest*`, the per-component driver (`sheaf_cohom::spectral::spectral_components`) and relative H⁰ (`sheaf_cohom::relative`). Every report echoes `zero_policy` and `zero_rule`.
- **Removed:** the `dim_global_rule` and `dim_strict` counts and the numpy-default cut for δ⁰. `near_zero_ambiguous` is now a flag on the same rule: the cut is not consistent with the certified residuals, or a counted-zero pair is numerically resolvable from 0 (λ > max(64ε‖L‖, 10r), a "near section"), or the nearest nonzero eigenvalue is within 10× of its cut.
- **Not changed:** rank(δᵏ) for k ≥ 1 (H¹/H² of 2-complexes) still uses the numpy-default SVD cut; that is not an H⁰ rule and stays documented in `tol.rs`. The m1 Python reference keeps R-H0-08 (it is gated).
- **Wording:** "certified residual" (ADR-008), never "certified rank". The Weyl check is `cut_consistent_with_residuals()`.

## ADR-016 — m1 contract math ported to Rust on faer (2026-10-20 Bahia) [ACCEPTED as a port; gated m1 behaviour unchanged]
- **Context:** Manu approved ROADMAP item 1 (2026-10-20): port m1's δ⁰, L₀, energy/discord and H⁰, which existed only in the Python reference `sheafcore.py`, to Rust on faer, through the gate, without changing gated behaviour; `sheaf-contract/src/lib.rs` stays byte-identical.
- **Where:**
  - `sheaf_cell::m1`: δ⁰ from per-edge row-major maps with heterogeneous stalks (`coboundary`; nonzero map entries only, self-loop blocks summed with a cancelled entry kept as a stored 0.0 like scipy `sum_duplicates`, sorted columns), `laplacian` (Gustavson δᵀδ) + `drop_zeros` (scipy SpGEMM never stores an exact 0), `edge_discord`/`edge_discord_dense`, `energy`/`energy_dense`, `pattern_components` (bipartite nonzero pattern) and `adjacency_components` (L's undirected nonzero adjacency).
  - `sheaf_cohom::m1`: `global_sections` (pass 0: every component's R-PERF-18 caps before anything is densified; pass 1: singular values only, global σ_max, ranks; basis-size cap; pass 2: one dense block at a time, V only, kernel = V[:, r..]) and `h0_dim_from_l` (per-component `self_adjoint_eigenvalues(Lower)` = numpy `eigvalsh` UPLO='L'). Component order = the reference's (scipy labels: by smallest row, row-less columns after), so G's column blocks come in the same order. Dense SVD/EVD only: no CG (so the faer `InitialGuessStatus::Zero` bug cannot apply) and no LSMR (no `lsmr_scratch`).
  - `sheaf_cli::m1::analyze(&SheafSpec)`: the Rust m1 path; `sheaf_cli::ffi`: a C ABI (cdylib `libsheaf_cli.so`) for `python/m1/sheafimpl_rs.py`, the Rust-backed builder shim used to run the gate. The gate's builder sandbox forbids fork/exec, so the Rust core is reached in-process via ctypes. The shim is the reference file with only math lines replaced; all 27 `# M:` guard lines, their messages and their order are unchanged.
- **H⁰ policy (amends ADR-015, still ONE typed rule):** `H0Policy` gains the contract-named variants `ZeroCut::GlobalSigma{rel}` (`H0Policy::m1_svd(1e-9)`: σ nonzero iff σ > rel·max(1, σ_max over all components)) and `ZeroCut::GlobalLambda{rel}` (`H0Policy::m1_eig(1e-10)`: λ zero iff λ ≤ rel·max(1, max|λ|)), i.e. R-H0-08 exactly. Only `sheaf_cohom::m1` accepts them (and it refuses any other policy); `zero_tol()` panics on them (documented) so they can never leak into the per-component v1 paths. The tolerances come from the caller (`tol` argument of the reference API).
- **R-NONID-04 (open item F8, resolved for the m1 path):** `analyze` runs `SheafSpec::validate()` first, so an identity/identity map pair is rejected (`R-NONID-04`) before any math, on every edge type; a one-sided identity is accepted. m1 maps are required per edge type (R-TYPE-03), so there is no identity default on the m1 path. v1's identity default for omitted maps is unchanged and v1-only. Tests: `m1_path_follows_r_nonid_04`, and the 5 `identity_pair_*` rows of the mutation table through `analyze`.
- **Parity (box runs 2026-10-20):** see `crates/sheaf-cli/tests/m1_parity.rs` (tolerances stated there) and `python/m1/rust_parity.py`. Results: 67 generated cases + the m1 spec fixture + the 114-row mutation table (0 mismatches, 24 accepted rows compared value by value, 5 R-NONID-04 rejections); worst |δ| 2.2e-16, rel L 3.6e-16, rel discord 7.6e-16, |E| 1.4e-14, projector 8.4e-11; H⁰ dims exact on both routes. Through ctypes: 187 value cases (incl. caller tols −1/0/1e-3, raw duplicate CSR, cancellation, self-loops, G column order) + 25 invalid-input cases (same exception type and message), 0 mismatches. Rust tests: 60 pass (54 before + 6 in `m1_parity.rs`).
- **Gate (FRONTIER `harness/gate.sh m1-rust-core`, started 2026-10-20 09:49 Bahia, shim + `libsheaf_cli.so` installed as `builder/sheafimpl.py`; the original Python builder is restored afterwards):** public suite 297/297 passed, models mutants 19/19, Rust mutants 23/23, gated sheaf-contract tests 10/10. conform reports two FAILs that do not depend on the builder: R-DATA-23 (corpus topics `evals`/`model-merging`/`latent-reasoning`/`compilers` not in the locked validate.py taxonomy) and R-ID-01 (rebuilding spec.json from the current corpus gives 496 nodes vs the frozen 382). Both come from the 2026-10-19 08:31 research batch (FRONTIER c58c599), which also edited `tools/contract/validate.py` without relocking, so at FRONTIER HEAD conform stops at R-LOCK-21 for ANY builder, the Python reference included. For this run validate.py was temporarily set to its locked bytes (sha 220c9763…, as at R18) and put back afterwards; nothing in the gate was relaxed. Holdout and the 8-reviewer consensus: see the final report of this change (they were still running at commit time). Note: the reviewer bundle contains only the shim, not the Rust sources (lanes.sh copies a fixed file list), so the reviewers cannot inspect the Rust math.
- **Review (Grok CLI, read-only strict sandbox, trimmed prompt; round 1 of ≤2):** 1 critical with repro, fixed: a wide block's zero-padded singular values were not counted when the caller's tol < 0 (`global_sections([[1,1]], tol=-1)` gave (2,1) vs the reference's (2,0)); now padded to `cols` entries, pinned by a test. Majors: FFI size arithmetic unchecked (fixed: checked_mul/add); `_csr_args` could tell Rust nnz > buffer length for a tampered `indptr` (fixed: length check → R-DIM-02 ValueError); "laplacian drops exact zeros unlike scipy" refuted (scipy's SpGEMM stores 2 entries on the repro, as we do; pinned by `m1_laplacian_cancellation_pattern`). Minors fixed: parity gaps (caller tols, raw duplicate CSR, cancellation, self-loops, G column-block order) and non-convergence → `np.linalg.LinAlgError`.
- **Known differences (documented, not gated):** a non-float (e.g. integer) δ gives a float64 L (the reference keeps the input dtype); L's CSR has sorted columns (scipy's SpGEMM order is unspecified); a multi-dimensional kernel basis is a different orthonormal basis of the same subspace (singular vectors are unique only up to rotation).

## Open items
- **R-NONID-04 vs v1 identity default (review F8):** resolved for the m1 path by ADR-016 (`sheaf_cli::m1::analyze` enforces R-NONID-04 via `validate()` before any math). Still open: two named parsers (`m1`, `v1`) and the v1 report stating which contract was used.
- Review flags still open after the restructure: F5 (`ValidSheaf` newtype), F9 (one sparse `CochainComplex`; cohomology is still dense), F11 (library `expect`s → `Result`), F13 (rayon as a feature), F15 (inner product). F4 is done (`sheaf_cell::num::family`); F12 partly (`krylov::SolveTol`).
- The graph Hodge CG (`sheaf_cohom::hodge_graph`) and the shift-invert inner PCG keep their own loops; only harmonic extension uses `sheaf_cell::num::krylov::cg` so far.
