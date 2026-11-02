# Consolidated disposition — ontology v0.1 → v0.2 (2026-11-02)

Inputs: three independent adversarial seats (Fable 5.1, fresh context each) + the builder's own operadic pass.
- seat 1 math soundness (`seat1-math.md`, 19 findings) · seat 2 Touchstone conformance (`seat2-touchstone.md`, 18) ·
  seat 3 spec coverage + operadic structure (`seat3-coverage.md`, 20 + a 40-row capability trace) · builder pass (`meta-operad.md`, 4)
- Overlap: the dangling Extend output was found by seat 3 (S3-01) and the builder (MO-01) independently. The
  missing Consistent-branch law was found by seat 2 (S2-04, by planting a mutant) and seat 1 (S1-01, by the
  clean-twin argument), from different directions.

Dispositions: **FIXED** (the law, test or file named), **HOLE** (now an explicit Hole law; propagates upward),
**TICKET** (stated as a law, build deferred to the named ticket), **REJECTED** (with reason).

| finding | severity | disposition |
|---|---|---|
| S1-01 G5 false in w-mode (discord = ‖w_e‖²) | P0 | FIXED: A.quantize.disc-source; G5 restated and re-derived |
| S1-02 G4 only for scalar maps | P0 | FIXED: S.circulate.domain restricts Circulate; G4 restated; derives from it; single-oracle flag |
| S1-03 two-oracle law wrong by cycle length | P0 | FIXED: I.circulate.quantize restated with c^T (C C^T)^-1 c and its eigenvalue interval |
| S1-04 v1 report bytes vary with thread count; cited test was m1 | P0 | FIXED: G3 restated (pinned threads), S.analyze.thread-pin, S.bind.binding binds the count; m1 citation removed |
| S1-05 extension pins rounded to 12 sig digits | P0 | FIXED: S.analyze.extension-pins compares sig(pin) |
| S1-06 basis size ≠ core h1 in general | P0 | FIXED: A.circulate.basis-size restricted to the constant R¹ identity sheaf |
| S1-07 batching masks small-component loops | P0 | HOLE: S.analyze.hodge-local; G6 is a Hole by propagation; P.recon.per-key-join carries the scale caveat |
| S1-08 extension certificate and policy ignored by Cert | P1 | FIXED: A.certify.converged, A.certify.unique-null, S.certify.policy |
| S1-09 holonomy invisible to Classify | P1 | TICKET T2: A.classify.h0-policy (MapConflict) + P.recon.rate-triangle, P.consensus.frame-holonomy, P.memory.contradiction-parity |
| S1-10 15 of 25 Formal leaves had no Bend law | P0 | FIXED as a check: `formal_ids_match_across_bend_haskell_rust` + LAWS.md backing table (10 of 41 backed now); the rest are tickets T2–T7 |
| S1-11 PROOF.bend did not check as committed | P0 | FIXED: adopted the seat's checked copy verbatim (arg order, qualified constructors, mono_tolerance moved to OPEN.bend) |
| S1-12 quantize laws underdetermined, overflow | P1 | FIXED: A.quantize.compare-once, A.quantize.saturate, S.quantize.floor restated with report rounding and row filter |
| S1-13 shares sum to 0 on w = 0 | P1 | FIXED: S.certify.wf restated (norm2 = 0 branch, 2·cert bound) |
| S1-14 no component identification | P1 | FIXED: I.judge.join (union-find on the spec sent, count checked) |
| S1-15 circ_only weaker than prose | P2 | FIXED: G4 proviso (still certified, floor cleared); derives from S.quantize.floor |
| S1-16 classify_total vacuous as Formal | P2 | FIXED: S.classify.total is Types; join_is_max registered as S.join.is-max (Property) |
| S1-17 cycle basis depends on forest choice | P2 | FIXED: A.circulate.forest; G3 derives from it |
| S1-18 core prefers w when both sent | P2 | FIXED: S.encode.one-cochain |
| S1-19 extension routes and caps differ | P2 | FIXED: folded into A.certify.unique-null and S.analyze.extension-minimal |
| S2-01 laws frozen before they ran | P0 | FIXED: core.sha256 covers the core only; laws.lock written by the gate only after ALL PROOFS CHECK; OPEN.bend holds open claims |
| S2-02 no Contract / Oracle / bound Receipt | P0 | FIXED in the ontology: Judge (TsVerdict) and Bind (Receipt bound to hashes) are separate stages; S.bind.*; **TICKET T6–T8** for the Rust objects; Contract object: TICKET T11 |
| S2-03 builder writes laws, tests, golden, mutants | P0 | HOLE: G10 (V.sep.signed, V.sep.fence, V.sep.family); TICKET T11 |
| S2-04 Consistent→BROKEN survived every layer | P0 | FIXED: consistent_iff + holds_iff_consistent in Bend, Haskell and Rust; 5 constructor mutants now killed by laws (`kit/receipts/constructor-mutants.txt`); Bend M6/M7, Haskell H6/H7 in the gate |
| S2-05 mutation layers 1, 3, 4 missing | P1 | FIXED for L1 (gate 7b, 7 Haskell mutants); HOLE for L3, L4 (V.mut.L3, V.mut.L4 → G9) |
| S2-06 mutant judged by text diff vs a failing baseline | P1 | FIXED: baseline must be ALL PROOFS CHECK; killed iff the mutant no longer checks |
| S2-07 bridge was 10 rows | P1 | FIXED: bridge B over the full domain, 19 234 rows, Haskell == Rust == golden (sha256 f92e5865…); 10 rows kept as smoke |
| S2-08 interface.lock lacked type hash, effect | P1 | FIXED: type_fnv64, effect, bridge field; auto-approve refused while any hash is null |
| S2-09 laws without oracles | P1 | FIXED as reporting: LAWS.md status + backing tables; leaves keep their evidence field; full Oracle objects TICKET T6 |
| S2-10 no monotone HOLDS ledger | P1 | FIXED: every gate step appends a hash-chained ledger line with its log hash; monotone check over the ledger: TICKET T11 |
| S2-11 PARTIAL exited 0 | P1 | FIXED: no PARTIAL; exit 0 / 1 BROKEN / 3 BLOCKED |
| S2-12 "proved" wording before any check | P1 | FIXED: kit wording; S.bind.tier; LAWS.bend header states exactly which laws were checked and by whom |
| S2-13 stale Python / .mjs leaking from older docs | P1 | FIXED: G11; experiments moved to exploration-python with a README; VERIFY_LOOP marked superseded; SPEC status note; QUEUE rewritten; test `no_python_on_the_gate_path` |
| S2-14 join discards a BROKEN under a HOLE | P2 | FIXED in the ontology: S.join.no-drop, I.join.bind; Rust object TICKET T7 |
| S2-15 agree-bit source and cover source | P2 | FIXED: I.circulate.quantize, I.judge.join |
| S2-16 join_is_max not in Bend | P2 | FIXED: registered as Property S.join.is-max |
| S2-17 same model family verifies | P2 | HOLE: V.sep.family (G10) |
| S2-18 gate logs discarded | P2 | FIXED: per-step logs under products/ledger/runs/<commit>/ |
| S3-01 Extend dangling | P0 | FIXED: Analyze absorbs Extend (`Analyze : ValidSpec × CorePolicy × Boundary → Outcome`); test `no_dangling_outputs` |
| S3-02 no Explain stage | P0 | FIXED: Explain stage + S.explain.pure, S.explain.wording, A.explain.loop-blame, A.explain.resolve |
| S3-03 no Diff stage | P0 | FIXED: Diff stage + S.diff.keyed/partition/reflexive/exit |
| S3-04 no crash path | P0 | FIXED: Outcome = Report \| CoreError; S.analyze.outcome; I.analyze.certify-error |
| S3-05 write-back unrepresentable | P0 | FIXED: Fill stage + S.fill.determined-only/ask/idempotent; P.memory.writeback reparented |
| S3-06 no adapter stage | P0 | FIXED: Ingest stage + S.ingest.*, I.ingest.encode, P.memory.status-table |
| S3-07 unique is three-valued | P1 | FIXED: A.certify.unique-null |
| S3-08 Attribute could not name a loop | P1 | FIXED: Attribute takes Circ; I.circulate.attribute |
| S3-09 no surface equivalence | P1 | FIXED: G8 + S.surface.envelope |
| S3-10 11 wired edges without glue | P1 | FIXED: 22 of 22 wired edges have a glue law; test `every_wired_edge_has_a_glue_law` |
| S3-11 `top` truncation | P1 | FIXED: S.analyze.top-complete, I.analyze.attribute |
| S3-12 policy side channel | P1 | FIXED: CorePolicy and DomainPolicy are separate typed inputs; G3 restated over (spec, CorePolicy, Boundary) |
| S3-13 nulls had no law | P1 | FIXED: I.analyze.quantize whitelist; A.attribute.context-only |
| S3-14 schema, translation, dated rates | P1 | FIXED: Prod::Schema, Prod::Translation with leaves; P.recon.dated-rate |
| S3-15 consensus attribution | P1 | FIXED: P.consensus.outlier-rank, P.consensus.round-diff |
| S3-16 Certify and Validate too shallow | P2 | FIXED for Certify (three atoms); Validate keeps its differential law with core-cited evidence |
| S3-17 misnamed glue laws | P2 | FIXED: I.circulate.quantize, I.attribute.judge |
| S3-18 G0 derivation did not imply it | P2 | FIXED: G0 derives from S.bind.manifest (the gate check) |
| S3-19 Studio interaction laws | P2 | FIXED: P.studio.steps, P.studio.pin-is-boundary |
| S3-20 no spec-coverage fixture | P2 | FIXED: src/capabilities.rs (41 rows, 2 deferred with reasons) + test `every_spec_capability_has_stage_and_law` |
| MO-01 Extend dangling | P0 | = S3-01 |
| MO-02 no dangling-output test | P1 | = S3-01 test |
| MO-03 no cost/budget valuation | P1 | TICKET: budgets are a Touchstone Contract field (T11) |
| MO-04 evidence matched by string | P2 | FIXED in part: Formal evidence now must match law, def, Haskell and Rust names at once |

Counts (tallied from this table, not by hand): 61 findings · 54 FIXED (5 of them leave a part as an explicit HOLE or TICKET) · 3 HOLE · 2 TICKET-only · 2 duplicates of another row · 0 REJECTED.
Result of re-verification after the fixes: ontology 17/17 checks, sheaf-kit 21/21, coverage sheaf H⁰ = 1 and
extension unique on 178 nodes, contradiction sheaf unchanged (5 Spec v1 laws retracted, residual 1.5e-30).

What this does NOT establish: the two new Bend laws are unchecked until the Mac gate runs; Haskell has not been run
in this environment; "no gaps" means no gap these four reviewers, the typed checks and the two sheaf runs could find.
