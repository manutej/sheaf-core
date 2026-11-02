# Operadic gap pass on the law ontology (meta-operad moves 1–5) + sheaf-core self-check

Run 2026-11-02, before the three independent checker seats report. Author: the ontology's builder, so this is
NOT the independent check — the seats are. Kept separate so the two can be compared.

## 1 TREE — every edge typed
Holds. `root_wiring_is_well_typed_and_every_slot_is_fed_once` passes: 19 wired slots, each typed, each fed once.

## 2 COMPOSE witness — the composite must reproduce the product
| id | severity | finding | evidence | fix |
|---|---|---|---|---|
| MO-01 | P0 | **Extend's output is consumed by nothing.** Classify's `bounded` bit (extension unique) and Attribute's fill/write-back have no wired source, so `A.classify.undetermined`, `A.attribute.fill`, `P.memory.writeback` describe data the pipeline never delivers. | `grep -c "From(Extend)" src/wiring.rs` → 0 (Join also 0, legitimately: it is the root) | Merge Extend into Analyze: `Analyze : (ValidSpec, Boundary) → Report`, Boundary possibly empty. Justified by the core's own law `empty_boundary_equals_absolute_h0_on_the_catalogue` (cite it as the glue law). Add test `no_dangling_outputs` (every op output consumed, except the root). |
| MO-02 | P1 | Test suite cannot catch MO-01: no law or test states "every output is consumed except the root". | tests/ontology.rs has no such test | add `no_dangling_outputs` |

## 3 COLLAPSE / 4 CHECK — direct vs decomposed
- Depth-2 collapse of the pipeline = one core call with `--boundary` on the whole spec (direct) vs stage-by-stage (decomposed). For the core-owned part this IS G6 (batch vs split, differential tier), so the collapse check exists as a law but is not yet executed: open.
- On the spec itself the collapse is "state the global law directly" vs "derive it from its parts". The derivation checks + the weakest-link tier valuation run as tests; they caught 3 overclaims (G1, G2, S.circulate.exact) in the previous round.

## 5 VALUATE — same tree, other algebras
| algebra | status |
|---|---|
| tier (weakest link) | executed as a test |
| evidence support | executed with sheaf-core (below) |
| cost / budget per stage | **missing** — no stage carries a budget; Touchstone's loop control needs one (MO-03, P1) |
| provenance (which artifact discharges which law) | partial: `evidence` names, matched by string (MO-04, P2: the probe is untyped; a renamed Rust test silently drops evidence to 0) |

## sheaf-core run on the ontology itself (products/laws/selfsheaf/, emitted by `cargo run -- selfsheaf`)
**Coverage sheaf** (92 nodes: ROOT + 68 laws + 11 stages + 12 colors; 207 edges; constant R¹; 43 leaf laws pinned to evidence: 8 executed here = 1, 8 Rust-checked with Bend pending = 0.5, 27 planned = 0):
- H⁰ = 1, components = 1 → no island: every law, stage and color is connected to the root.
- extension unique (h0_relative = 0), converged → every node rests on at least one pinned leaf.
- support scores (harmonic extension, 0 = rests on planned leaves only, 1 = fully executed):
  G2 0.045 · Encode 0.055 · Attribute 0.088 · G0 0.107 · G7 0.107 · Quantize 0.117 · Circulate 0.119 · Validate 0.152 · Certify 0.175 · Receipt 0.204 · Join 0.281 · G1 0.337 · Classify 0.344 · G4 0.354 · G5 0.482 · Extend 0.54 · G3 0.544 · G6 0.589 · Analyze 0.70.
  Reading: the build frontier is Encode → Attribute → Quantize/Circulate (the product side); the core-owned stages are carried by sheaf-core's own tests. G0 scores low only because its gate script is not yet registered as evidence (the script ran: G0 PASS).
- The coverage sheaf cannot see MO-01 (edges are undirected), which is why the COMPOSE witness is a separate check.

**Contradiction sheaf** (our 8 globals + the laws they touch + Spec v1's laws L1–L10, its 5 % share threshold and its Python tool layer; agrees = +1 map, contradicts = −1 map; our globals pinned adopted):
- extension unique → every v1 law gets a forced value.
- survive (+1): L1, L2, L3, L4, L8, L9, L10, null-p95-reported-not-gated.
- retract (−1): L5, L6, L7 (share-based verdicts), the 5 % share threshold, the Python tool layer.
- residual energy 1.5e-30 → after those retractions nothing irreducible remains: the surviving laws of both specs are consistent.
- Caveat: the agree/contradict labels are inputs written by the builder; the run checks the label set is consistent, not that each label is right.
