# Ticket queue (both arms run the same queue; VERIFY_LOOP §4)

Each ticket: goal · acceptance (gates + named cases) · forbidden. Every ticket forbids editing `crates/**`
and `products/bend/LAWS.bend` (G0). Dependencies run top to bottom.

| id | goal | acceptance |
|---|---|---|
| T0 | Make `bend products/bend/PROOF.bend` run; fix draft syntax in verdict.bend/PROOF.bend only | G1 output is exactly the expected SOME PROOFS FAIL naming only xonly_no_loop, mono_tolerance |
| T1a | Prove `xonly_no_loop` | G1 ALL-but-mono; `--verdict` agrees |
| T1b | Prove `mono_tolerance` (needs a `Nat.is_ge` monotonicity lemma) | G1 ALL PROOFS CHECK + `--verdict` |
| T2 | Domain compiler `domain-v1 → frontier-sheaf/v1` (TS) + propose its laws (human approves text) | E1–E4 rebuilt from domain-v1 inputs; G2 |
| T3 | Report bindings generated from `report.rs`; G2 validates every report | a renamed field fails G2 |
| T4 | MCP server: kb_check, kb_fill, kb_explain, kb_diff over the `sheaf` binary, classifier = compiled verdict.mjs | E1a–c through MCP; CoreCrash path tested |
| T5 | graph-wiki / llm-wiki markdown adapter | 300-note import < 2 s; findings on planted contradictions |
| T6 | recon CLI + synthetic 10k-key fixture (50 slips, 5 rate loops, clean twin) | recall ≥ 0.98 at rank 1; 5/5 loops, blame = loop; 0 on clean twin |
| T7 | consensus_check tool + multi-agent-eval hook | E3a–c; recorded-panel measurement reported |
| T8 | Studio v0: verdict card, graph heat, boundary mode, numbered vertical steps | E4a reproducible by clicking; design review |
