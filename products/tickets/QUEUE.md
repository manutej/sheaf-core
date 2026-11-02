# Ticket queue (v0.2, after the three review seats). Both model arms run the same queue.

Every ticket: goal · acceptance · forbidden. Every ticket forbids editing `crates/**` (G0) and, once signed,
`products/bend/LAWS.bend` and `laws.lock` (V.sep.fence). Law ids refer to products/laws (LAWS.md).

| id | goal | acceptance |
|---|---|---|
| T0 | Run `products/gates/atom_gate.sh` on the Mac; fix only PROOF.bend / main.bend syntax if the two new laws (consistent_iff, holds_iff_consistent) fail | exit 0; ledger lines HOLDS for every step; laws.lock written |
| T1b | Prove mono_tolerance (OPEN.bend): Nat.is_ge lemma, then move the law back into LAWS.bend | ALL PROOFS CHECK with it included; A.classify.mono backed 1 of 1 |
| T2 | Bend + Haskell + Rust for the remaining Formal leaves of Classify/Judge: h0_map_conflict (adds MapConflict) | formal_ids_match passes; MapConflict mutant killed |
| T3 | Circulate atoms: circ_exact_zero, circ_linear (signed values as Bool+Nat), forest_deterministic, rate_loop_closes | backed at all four layers; I.circulate.quantize differential on E2a and S2 fixtures |
| T4 | Quantize: compare-once, saturate, floor, disc-source (Rust types + properties) | G5 clean-twin fixtures HOLDS in x-mode and w-mode |
| T5 | Encode + Ingest: encode_valid, encode_index_bij, encode_mode, ingest_ids_stable, status_table | I.encode.validate: every fixture accepted by the core binary |
| T6 | Certify: report bindings from report.rs, Outcome sum, unique three-valued | a renamed report field -> BLOCKED(report-schema) |
| T7 | Attribute + Explain + Diff + Fill in Rust with their laws | capability rows C02-C07 backed |
| T8 | MCP server (kb_check / kb_fill / kb_explain / kb_diff) over the binary; CLI = MCP (G8) | byte-identical Findings over both surfaces |
| T9 | Typed Rust oracle replacing exploration-python (independent linear algebra, not faer) | differential vs core on the 40-sheaf sweep |
| T10 | V.mut.L3 interface mutation and V.mut.L4 planted-bug corpus with hidden holdout | G9 leaves Hole |
| T11 | V.sep: hook fence, laws.lock signing flow, one non-Claude review seat | G10 leaves Hole |
