# sheaf-products — product specification v0.1

Branch `verify-loop/products-v0` (worktree of `manutej/sheaf-core`). Date 2026-11-02.
Status: DRAFT for build. Everything below the line "core is untouched" is new code; `crates/` is read-only here.

---

## 0. The thesis in one paragraph

sheaf-core answers one question that almost no data tool can answer with attribution: **is this set of
related facts globally consistent, and if not, is the problem one bad source (fixable) or a loop of
claims that cannot all be true (irreducible)?** The Hodge split of the observed edge cochain is the
mechanism: the *exact* part is explained by some node assignment (the potential points at the bad node),
the *harmonic* part is not explainable by any assignment (the loop itself is wrong). Harmonic extension
with the relative-H⁰ certificate adds the second capability: **fill unknowns from trusted anchors, and say
whether the fill is determined.** Every product below is one domain adapter on these two capabilities.
The products sell the *verdict with a certificate*, never a bare score.

## 1. Evidence base (run 2026-11-02, this worktree)

| check | result |
|---|---|
| `cargo build --release --workspace --locked` (rustc stable 1.97; pin 1.85.1 unreachable) | OK |
| `cargo test --release --workspace --locked` | 70/70 pass |
| `products/experiments/run_experiments.py` (black-box tool use, no core changes) | **25/25 checks pass**, 0.54 s |
| E4b: 40 random graph sheaves, random dense restriction maps, dims 1–3, vs numpy oracle | H⁰ exact on 40/40; harmonic share worst abs err 4.5e-13; extension energy worst rel err 2.5e-12 |

Experiment ↔ product map (each is an acceptance test for its product; specs + reports are in `experiments/out/`):

| id | product | what it proves | ground truth |
|---|---|---|---|
| E1a | P1 | one corrupted claim → ranked #1 by incident discord; exact share 1.0 | hand |
| E1b | P1 | odd cycle of "contradicts" maps → H⁰ = 0; pin A=1 → B=C=−1/3, irreducible energy 4/3, unique | hand |
| E1c | P1 | fill unknown claims from anchors; orphan claim → `unique=false`, h0_relative = 1 | hand |
| E2a | P2 | triangular FX arbitrage → top-3 harmonic edges = the loop; arbitrage-free control → share 1e-33 | numpy oracle |
| E2b | P2 | one system's unit slip across USD/EUR/MXN maps → wrong system ranked #1 | hand |
| E3a | P3 | outlier agent among 5 (3-d beliefs) → ranked #1; consensus space H⁰ = 3 | hand |
| E3b/c | P3 | frame holonomy: 120° loop → H⁰ = H¹ = 2; 90° loop → H⁰ = H¹ = 0 (no consensus in any frame) | hand |
| E4a/b | P4 | extension on a path; 40-case randomized oracle sweep | hand + numpy |

**Findings that change the design** (from the runs, not from reading):

1. **Harmonic share is the wrong alert metric.** In E2a the arbitrage is correctly localised but
   `harmonic_share = 1.96e-5`, because large exact components (log prices) dominate ‖w‖². Alerts must use
   absolute per-loop harmonic energy against a domain tolerance (`harmonic_energy ≥ τ_domain`), with share
   reported only as context. This is encoded in the classifier laws (§4.2).
2. **Harmonic attribution is per-loop, not per-edge.** In E2a the three loop edges carry *equal* harmonic
   energy. The UI must say "this loop is inconsistent, no single edge is to blame", never "edge X is wrong".
   Single-edge blame is only valid for the exact part (node potential / incident discord).
3. **Incident discord, not |potential|, is the default node ranking** whenever maps rescale units
   (E2b, A11 in the report): `|potential|` is in node-stalk units and mis-ranks across currencies.
4. **`x` and `w` are different products.** With only node data `x`, `w := δx` is exact by construction (A9),
   so harmonic findings are impossible. Loop conflicts require independent edge observations `w`. Every
   adapter must declare which it supplies; the UI must not show a "loop conflicts: none" badge for x-only
   inputs (it would read as a negative finding when it is structurally impossible).
5. The harness itself had a bug (wrong JSON key) caught on first run. Black-box adapters need a typed
   report binding generated from `report.rs`, not hand-written key strings (→ §3.2).

## 2. Architecture: the core as a tool

```
L3 products       P1 Memory Integrity · P2 Reconciliation Auditor · P3 Consensus Monitor · P4 Studio · P5+ adapters
L2 surfaces       MCP server (agent tool calls) · CLI · HTTP (later) · Studio web app (WASM later)
L1 sheaf-kit      domain compilers (domain → frontier-sheaf/v1 spec)  ·  verdict classifier (report → Finding)
                  ├── verdict logic: Bend (LAWS.bend + PROOF.bend), compiled to .mjs/.c   ← proof gate G1
                  └── adapters/IO: TypeScript (MCP SDK) calling the `sheaf` binary          ← oracle gate G2
L0 sheaf-core     crates/* UNTOUCHED, pinned by commit 825a23b + sha256 of crates/**     ← invariance gate G0
```

**Why this split.** Bend has `Nat/U32/F32` and no linear algebra, so it cannot (and should not) re-prove
the numerics; the Rust core already carries certified residuals for those (ADR-008). What *can* be proven
mechanically, and is where agent-written code actually goes wrong, is the **discrete decision layer**: how
numbers become verdicts, when the system is allowed to claim anything, and the structural invariants of the
compiled spec. That layer is small, pure, total and exactly what an auditor reads. Bend compiles it to an ES
module (`bend verdict.bend -o verdict.mjs`) which the TS surface imports, so the proven code *is* the shipped
code, not a parallel model of it.

Process boundary: L1 calls L0 as a subprocess (`sheaf analyze --out --boundary`). No FFI in v0: the binary is
the stable contract, it keeps the core byte-identical, and it gives free crash isolation (F11 `expect` panics
become a typed `CoreCrash` finding instead of killing the server).

## 3. Shared contracts

### 3.1 Domain graph (adapter input, product-agnostic)
```jsonc
{ "format": "sheaf-kit/domain-v1",
  "entities":  [{ "id": "crm:acct-42", "kind": "system|claim|source|agent|...", "dim": 1, "frame": "USD" }],
  "relations": [{ "u": "...", "v": "...", "kind": "supports|contradicts|converts|agrees|...",
                  "map": { "u": [[..]], "v": [[..]] } | "identity" | { "scale_u": 1.0, "scale_v": 0.92 },
                  "observation": [..] | null,   // independent edge evidence -> w (enables loop findings)
                  "source": "provenance string" }],
  "values":    { "crm:acct-42": [1000.0] },      // node evidence -> x
  "anchors":   { "s1": [1.0] },                   // trusted values -> --boundary
  "tolerance": { "harmonic_abs": 1e-6, "discord_abs": 1e-6, "h0_absolute_tol": null } }
```
Compiler rules (each is a law candidate, §4.2): every relation endpoint exists; map shapes are
`de × dim(u)` / `de × dim(v)`; `observation` present on all relations or none (mixed → reject, P_INPUT);
ids are stable and round-trip; the compiled spec's node order is the entity order.

### 3.2 Finding (classifier output, what every UI and agent consumes)
```jsonc
{ "format": "sheaf-kit/finding-v1",
  "verdict": "consistent | fixable_nodes | loop_conflict | underdetermined | inconclusive",
  "certificate": { "core_commit": "825a23b", "spec_sha256": "...", "report_sha256": "...",
                   "converged": true, "cut_consistent": true, "near_zero_ambiguous": false,
                   "laws_checked": "ALL PROOFS CHECK", "laws_sha256": "..." },
  "fixable":  [{ "entity": "ERP", "incident_discord": 1.2e5, "rank": 1 }],
  "loops":    [{ "edges": ["USD->EUR","EUR->MXN","MXN->USD"], "harmonic_energy": 3.0e-4, "blame": "loop" }],
  "fills":    [{ "entity": "c1", "value": [1.0], "determined": true }],
  "context":  { "h0": 1, "h1": 2, "harmonic_share": 1.96e-5, "null_frac_ge_observed": 0.02 },
  "x_only":   false,
  "wording":  "certified residual; not a proof that any fact is true" }
```
Report bindings: generate TS types from `report.rs` (serde field names) once, commit them, and add a G2
check that a fresh report validates against them. No hand-written key strings (finding 5).

### 3.3 Verdict rules (the classifier; formalised in Bend, §4.2)
Evaluated in order, first match wins:
1. core crashed, or `spectral.converged = false`, or `hodge.converged = false`, or
   `h0.cut_consistent_with_zero_tol = false`, or `h0.near_zero_ambiguous = true` → **inconclusive**
2. `--boundary` given and `extension.unique = false` → **underdetermined** (report the free entities)
3. not `x_only` and any loop with `harmonic_energy ≥ τ_h` → **loop_conflict** (plus fixable nodes if any)
4. any node with `incident_discord ≥ τ_d` → **fixable_nodes**
5. otherwise → **consistent**

Rule 1 is the certificate discipline made executable: the product may never claim a conflict, or its
absence, on an uncertified computation.

## 4. Products

### P1 — Memory Integrity Layer  (wedge; build first)
**User.** Teams and solo builders running agent memory / knowledge bases (llm-wiki, graph-wiki, Notion
research DBs, Obsidian, sheaf-memory). Buyer: anyone whose agents read their own notes back as truth.
**Job.** "Before my agent cites its memory, tell me which claims contradict their sources, which clusters
of claims cannot all be true, and fill missing statuses only where the evidence determines them."
**Model.** Entities = sources, claims, topics. Stalk dim 1 (status in [−1, 1]: refuted … supported), or
dim k for k-way categorical status. Maps: `supports` = identity; `contradicts` = −1 (E1b); `refines` =
projection. `values` = current statuses; `observation` = pairwise judgments when an extractor emits them
(e.g. an NLI score between two claims) → enables loop findings. `anchors` = human-verified facts.
**MCP tools** (agent-facing; the "programmatic tool call" path):
- `kb_check(graph | path, tolerance?) -> Finding`
- `kb_fill(graph, anchors) -> Finding` (fills + `determined` per entity)
- `kb_explain(finding, entity | loop) -> { why, evidence_edges, what_would_resolve }` (pure, from Finding)
- `kb_diff(finding_a, finding_b) -> { new_conflicts, resolved, drifted }` (status drift over time)
**Adapters v0:** llm-wiki / graph-wiki markdown (frontmatter `status`, `[[links]]` typed by link label);
generic JSON. v1: Notion database via MCP connector; Obsidian vault.
**UI (Studio skin).** Hero = verdict card with certificate chip. Left: ranked fixable claims with their
source evidence. Right: loop cards ("these 3 claims cannot all be true", blame = loop). Bottom: fill table
with determined / free badges. No Betti numbers on the hero; they live in a "math" drawer.
**Acceptance.** E1a, E1b, E1c pass through the MCP tool (not just the CLI); a 300-note graph-wiki import
classifies in < 2 s; `kb_explain` never names a single edge as the cause of a loop (law L-WORD-1).

### P2 — Multi-System Reconciliation Auditor
**User.** Finance ops, RevOps, data-quality leads; delivered through enterprise training/consulting.
**Job.** "The same customer/amount/SKU lives in CRM, ERP, billing and three spreadsheets. Tell me which
system is wrong, and which conversion chains can't all be right."
**Model.** Entity = (system, record key); dim = number of reconciled fields. Maps = unit/currency/FX
conversions into a common edge frame (E2b). `observation` = independently sourced rates or transfer
records (E2a: log-rates make arbitrage a harmonic cochain). Tolerances in business units (cents).
**Surfaces.** CLI `sheaf-recon run recon.yaml` (CSV/Parquet in, Finding + CSV of exceptions out); MCP tools
`recon_check`, `recon_explain`; scheduled run producing a dated exception report.
**recon.yaml**: systems, key columns, field→dim mapping, conversion table (static or dated), tolerance.
**Scale target.** 100k records × 5 systems: one sheaf per key is tiny; batch them as disconnected
components in one spec (the core is per-component already). Measure; if > 60 s, shard by key range.
**Acceptance.** E2a, E2b; a synthetic 10k-key fixture with 50 planted slips and 5 planted rate loops:
recall ≥ 0.98 on slips at #1 rank, all 5 loops reported with blame = loop, zero findings on a clean twin.

### P3 — Multi-Agent Consensus Monitor
**User.** People running multi-agent panels (your MoE/multi-agent-eval seats, Fable/Opus fleets).
**Job.** "Do my agents agree? Who is the outlier? Is disagreement a bad agent or a framing problem no
agent can fix?"
**Model.** Entity = agent (per question); dim = belief vector (probabilities over options, or a scored
rubric). Maps = identity when agents share a rubric; frame changes (rotations/permutations/rescalings)
when rubrics differ (E3b/c: holonomy decides whether consensus exists at all).
**Surfaces.** MCP tool `consensus_check(beliefs, frames?)`; a hook for the multi-agent-eval skill that
appends a Finding to each merge step; live panel in Studio.
**Acceptance.** E3a/b/c; on a recorded 5-seat eval: the seat the human merger overruled is ranked #1 in
≥ 80% of rounds (measured, reported, not asserted).

### P4 — Sheaf Knowledge Studio (shared UI shell + teaching tool)
**User.** You, your Silver cohort and enterprise classes first; then every P1–P3 customer as their UI.
**Job.** Upload or draw a graph + maps, see H⁰/H¹, the Hodge split, per-edge harmonic energy, pin
boundary values and watch the extension, read the Finding in plain language.
**Build.** v0: static web app (TS) + a tiny local runner that shells out to `sheaf` (works today). v1:
WASM core in-browser per ADR-004 (needs the core's F13 rayon feature gate: a core change request, §6).
Reuse stalks-and-sections' WebGL explorer and SheafGraph JSON as the visual layer.
**Views.** (1) Graph canvas with edge heat = harmonic energy, node halo = incident discord. (2) Verdict
card (same component as P1). (3) Boundary mode: click to pin, live extension, `determined` badges.
(4) Math drawer: spectrum, cut, residuals, nulls. Steps are numbered vertically; hovering a step number
highlights the matching graph region.
**Acceptance.** E4a/b reproduced by clicking; every number on screen traceable to a report field.

### P5+ — explored use cases (adapter only, same core)
| use case | entities / stalks | maps | x or w | verdict that sells |
|---|---|---|---|---|
| Sensor fusion / IoT drift | sensors, readings | calibration gains | w = pairwise co-measurements | which sensor drifted vs calibration loop broken |
| Ontology / schema alignment | concepts in 2+ ontologies | mapping matrices | w = alignment scores | inconsistent alignment cycles across 3 vocabularies |
| Config drift across services | service × setting | unit/format conversion | x | which service's config disagrees |
| RAG citation integrity | answer spans, citations, sources | supports/contradicts | w = NLI judgments | the answer's claims form an unsatisfiable loop |
| Financial close / intercompany | entities × accounts | FX | w = intercompany transfers | elimination loops that don't net to zero |
| Supply chain units | SKU per system | pack/unit conversions | x | the system with the unit slip |
| LLM eval rubric consistency | graders | rubric frame maps | x | grader outlier vs rubric that can't be reconciled |
| Translation memory / multilingual KB | term per language | translation maps | w | inconsistent translation cycles (en→es→fr→en) |

## 5. Editions (how the core ships)
| edition | form | depends on | when |
|---|---|---|---|
| E-Tool | MCP server + CLI over the `sheaf` binary | nothing new in core | now (v0) |
| E-Lib | PyO3 wheel `sheafkit` | core F11 (`expect` → `Result`) | after core CR-1 |
| E-Studio | WASM in browser | core F13 (rayon feature), wasm32 spike S1 | after CR-3 |
| E-Hosted | HTTP service, per-tenant | ADR-005 effect interpreter, Temporal later | after paying P2 pilot |

## 6. Core change requests (NOT done here; each goes through sheaf-core's own review)
- CR-1 F11: library `expect` on SVD non-convergence → `Result` (blocks E-Lib).
- CR-2 Sparse relative H⁰ (roadmap 3): `h0_relative = null` above `dense_max` today; P1 at scale needs it.
- CR-3 F13 rayon feature for wasm32 (blocks E-Studio v1).
- CR-4 Incremental recompute (edit one edge without a full rerun) for the Studio's live loop.
- CR-5 Report option to emit per-loop (cycle-basis) harmonic energy, so P2 can list loops directly instead
  of re-deriving cycles from edge rankings in L1.

## 7. Milestones (each closes only through the verification loop, VERIFY_LOOP.md)
| m | deliverable | gates |
|---|---|---|
| M0 | this branch: experiments 25/25, SPEC, loop design | G0, G2 ✅ · G1 blocked (Bend checker not yet runnable here) |
| M1 | `sheaf-kit` verdict classifier in Bend with laws L-SOUND-1..3, L-TOTAL, L-UNDET, L-MONO proven; compiled `.mjs` | G1 ALL PROOFS CHECK + `--verdict` |
| M2 | domain compiler + report bindings (TS); E1–E4 pass through it | G0 G2 G3 |
| M3 | MCP server exposing kb_*, recon_*, consensus_* tools | G0–G4 |
| M4 | Studio v0 (local runner) | G0–G4 + design review |
| M5 | P1 graph-wiki adapter on the real wiki; P2 synthetic 10k fixture | G0–G5 |
