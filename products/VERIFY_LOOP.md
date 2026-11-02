# Verification loop engine — design v0.1

Goal: every change an agent proposes reaches `main` only after it passes gates that the agent cannot
talk its way through, and every pass leaves a hash-chained record you can read later (the "history"
you asked to be able to view). Same loop for both model arms of the Opus-vs-Fable experiment.

## 1. Gates (run in order; any FAIL stops the ticket and returns the failure text to the builder)

| gate | checks | who/what decides | output |
|---|---|---|---|
| **G0 invariance** | `crates/**`, `Cargo.lock`, `rust-toolchain.toml` byte-identical to pinned commit 825a23b (sha256 manifest `gates/core.sha256`) | script | PASS/FAIL + diff list |
| **G1 laws** | `bend PROOF.bend` prints ALL PROOFS CHECK, then `bend PROOF.bend --verdict` (BendTT kernel, proved in Lean) also ALL PROOFS CHECK; `LAWS.bend` unchanged unless the ticket is a law-change ticket signed off by Manu | Bend checker | verdict text + `laws.sha256` |
| **G2 oracle** | `experiments/run_experiments.py` 25/25 + per-ticket acceptance cases; report bindings validate | numpy oracle, hand truths | results.json |
| **G3 build/test** | `cargo test` for the untouched core (70/70) + TS typecheck + unit tests of L1 | toolchains | logs |
| **G4 independent review** | a reviewer agent of the *other* model, with no access to the builder's transcript, reads the diff + LAWS + finding schema and must either APPROVE or file numbered findings with a repro; at most 3 rounds (ADR-004 cap), diff-only re-review | reviewer agent | review.md |
| **G5 publish** | append a ledger entry; tag; (optionally) `bend --publish` the verdict package to BendHub (public and permanent; needs your explicit OK each time) | script + Manu | ledger line |

**Laws are written by the human side, proofs by the agent side** (Bend convention). A builder agent may
not edit `LAWS.bend`; G0 extends to it. A ticket that needs a new law is split: "propose law" (human
approves the text) → "prove law" (agent).

## 2. The ledger (viewable history)
`products/ledger/ledger.jsonl`, one line per gate run:
`{ts, ticket, arm: "opus"|"fable", commit, gate, verdict, artifact_sha256, prev_sha256, wall_s, tokens?}`.
`prev_sha256` chains entries (ADR-005 EventRecord shape), so the history cannot be silently rewritten.
A Studio page renders it as a vertical, numbered timeline per ticket and per arm.

## 3. What Bend proves here (and what it does not)
Provable now, in `products/bend/`:
- **L-SOUND-1** not converged ⇒ inconclusive. **L-SOUND-2** ambiguous cut ⇒ inconclusive.
- **L-UNDET** (converged, unambiguous, boundary given, not unique) ⇒ underdetermined.
- **L-XONLY** x-only input never yields loop_conflict (finding 4).
- **L-MONO** raising the harmonic tolerance never turns `consistent` into `loop_conflict` (open claim: first ticket).
- **L-TOTAL** follows from the types: `classify` is a total function to a closed enum (Bend checks totality and termination).
- Later: compiler laws on the domain graph (endpoint existence, shape agreement, id round-trip) over `List`/`Map`.

Not provable in Bend and not claimed: floating-point eigen/Hodge numerics. Those stay "certified residual"
from the Rust core (ADR-008) and are cross-checked by the G2 oracle. The Finding's `wording` field says so.
Quantities enter the classifier as `Nat` in fixed units (e.g. harmonic energy in 1e-12 units, ceil'd), so
comparisons are exact; the quantisation step is TS code covered by G2 property tests, not by G1.

## 4. Opus-vs-Fable long-horizon measurement (independent test)
- **Setup.** Two worktrees off the same commit: `arm/opus` and `arm/fable`. Same ticket queue
  (`products/tickets/`), same gates, same wall-clock and turn budget. Each arm's builder is one model; its G4
  reviewer is the *other* model (so no model grades itself). Ledger records the arm on every line.
- **Tickets (long horizon, dependent chain):** T1 prove L-MONO; T2 domain compiler + its laws; T3 report
  bindings + G2 binding check; T4 MCP server (kb_*); T5 graph-wiki adapter; T6 P2 recon CLI + 10k fixture;
  T7 consensus hook; T8 Studio v0 verdict card + boundary mode.
- **Metrics (all read from the ledger, none self-reported):** tickets closed; first-try gate pass rate per
  gate; gate-failure count before pass; G4 findings per KLOC and severity; regressions (a later ticket
  re-failing an earlier ticket's acceptance cases); law-edit attempts (G0 blocks); wall time; tokens.
- **Ablation (the Bend question).** Run T2–T6 twice per arm: with G1 on, and with G1 off (laws as prose
  comments only). The delta in G4 findings and regressions is the measured value of proof gates.
- **Blinding.** The final quality review (end of chain) is done by a fresh agent given both repos with arm
  labels stripped.

## 5. Running it continuously
A scheduled task per arm: pull the next open ticket, run builder → gates → reviewer, append the ledger,
push the arm branch, and post a one-line status. It stops on: G0 failure (an agent touched the core or the
laws), three consecutive failures on one ticket, or any G5 publish (needs your OK).

## 6. Status right now
- G0 ✅ (no file under crates/ changed on this branch). G2 ✅ 25/25. G3 ✅ 70/70.
- **G1 blocked:** Bend is vendored from GitHub (bendlang/bend) and bun is installed, but executing the Bend
  checker was refused by this session's safety guard (it is code fetched from outside your repo); the
  installer host bend-lang.com is also outside the network allowlist. `products/bend/*.bend` are therefore
  **UNCHECKED drafts**: they have never been through `bend`, and nothing in this repo may describe them as proven.
