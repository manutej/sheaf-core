# sheaf-contract — m1 contract core (LOCKED; formerly sheaf-base)

`src/lib.rs` is byte-identical to FRONTIER `sheaf-experiments/contract/rust/src/lib.rs` (sha256 `7128f2b4…6de`, the hash in the
gate's `LOCK.public`; CI checks it). It changes **only through the gate** (builder → hidden evaluator → 8 reviewers).

What it is: the serde mirror of `sheaf.schema.json`, `validate()` returning requirement IDs, and per-edge-type `discord()`.
What it is not: it has no faer and no δ⁰/L₀/energy/H⁰ — those are in the Python reference. Porting them to faer is
ROADMAP item 1 and must go through the gate. Tests read `fixtures/m1/` (synthetic) via `.cargo/config.toml`.
Only sheaf-cli depends on it (the `m1` adapter); the crate rename did not touch `src/lib.rs` (Rust lib name `sheaf_contract`).
