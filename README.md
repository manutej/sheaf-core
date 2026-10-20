# sheaf-core

Cellular sheaf consistency in Rust on [faer](https://github.com/sarah-quinones/faer-rs): the m1 contract core, real
sheaf cohomology (H⁰/H¹/H²) with a Hodge split and scalable spectral analysis, and the `sheaf` CLI.
Python under `python/` is a **reference only** (parity checks); nothing at runtime depends on it.

Private. All rights reserved (see `LICENSE`). Layout per ADR-014 — **ACCEPTED** (Manu, 2026-10-19, after the
category-theory architecture review): five crates.

## Crates

```
sheaf-cell ──┬──> sheaf-cohom ──┐
             └──> sheaf-solve ──┼──> sheaf-cli
sheaf-contract ─────────────────┘
```
No other edges: sheaf-solve never depends on sheaf-cohom, sheaf-cell never depends on sheaf-contract.

| crate | what | status |
|---|---|---|
| `crates/sheaf-cell` | The objects. `CellComplex` (face poset; `graph()` keeps order, multi-edges and orientation; `simplicial()` sorts/dedups), `Sheaf` (per-incidence restriction maps, `from_graph`, validation), cochains, δᵏ, matrix-free L_k (`ops::HodgeOp`), energy and per-edge discord, graph sheaves with per-edge maps (`graph`: the one δ⁰ assembly), and `num`: the H⁰ policy type (`num::tol::H0Policy`), dense SVD/rank, CSR, CG (`num::krylov`), and the policy-carrying result types (`num::family`: `Exact` vs `Approximate`, `Numerical<T>` (f64) vs `ExactFp<T>` (F_p)). | Extension, not gated. |
| `crates/sheaf-cohom` | Hᵏ with harmonic bases and rank gaps, relative H⁰(X,B;F) (`relative::h0_relative`: the uniqueness certificate for harmonic extension), Hodge split with certificates, nearest global section, sparse CG graph Hodge split, per-component shift-invert Lanczos spectral solver (faer 0.24.4), exact Betti numbers over F_p. Depends on sheaf-cell only. | Extension, **not** through the 18-round gate. Exact family (ADR-013), certified residuals (ADR-008). |
| `crates/sheaf-solve` | Harmonic extension: x_I minimising ‖δx‖² with x_B fixed (dense Cholesky / CG, certified residual). Depends on sheaf-cell only. | Extension. Uniqueness is decided by the caller via sheaf-cohom. |
| `crates/sheaf-contract` | (formerly `sheaf-base`) m1 contract core, SPEC v1.18: serde mirror of `sheaf.schema.json`, `validate()` (requirement IDs), per-edge-type `discord()`. No faer. | **LOCKED.** `src/lib.rs` is byte-identical to the gated FRONTIER crate (sha256 checked in CI). Changes only through the gate. |
| `crates/sheaf-cli` | Composes everything. Pure v1 core `v1::analyze(&Spec) -> (Report, Vec<Effect>)` with a typed `Report`; `m1` adapter (sheaf-contract → sheaf-cell); `sheaf` binary = the imperative shell (ADR-005), incl. `--boundary`. | v1, graph-only. |

Wording: "certified residual" — each eigenvalue lies within ±residual of a true eigenvalue (Weyl). The zero/nonzero
labelling is a numerical policy; Krylov completeness is not certified. Never "proof", never "certified rank".

### Does m1 use faer?
No. The m1 Rust crate (sheaf-contract) depends only on serde/serde_json; m1's δ⁰/L₀/energy/H⁰ exist in the Python
reference `python/sheafcore.py`. sheaf-cli's `m1` adapter maps an m1 spec onto sheaf-cell, and `tests/contract_parity.rs`
pins m1 `discord()` to sheaf-cell's δ⁰ per edge. See ROADMAP: **port base to faer** is the top TODO, and must keep m1 behaviour identical
to the gated version.

## Build and test

```
cargo build --release --workspace
cargo test --release --workspace
cargo run --release -p sheaf-cli --bin sheaf -- analyze spec.json [--out report.json] [--k K] [--null-reps R] [--seed S] [--top N] [--method shift-invert|polynomial] [--boundary boundary.json]
```

`.cargo/config.toml` points sheaf-contract's gated tests at `fixtures/m1/` (env `SPEC_JSON`, `MUTATIONS_JSON`, `PY_DISCORD`,
`PY_DISCORD_RAND`). Those fixtures are **synthetic** (`python/m1/make_fixtures.py`): the frozen spec the FRONTIER gate uses
is a private Notion research-graph export and is not shipped. The mutation table is regenerated with the unchanged
generator, and `python/m1/check_table.py` checks it against the Python mirrors (114/114 rows).
The FRONTIER gate's fresh-seed parity draws are fixed seeds here.

Reference checks (need numpy + scipy; optional): `python python/parity/parity.py` (Python sheafcore vs Rust v1:
coboundary, Laplacian, energy, per-edge discord to 1e-9, H⁰ exact; last run worst 4.7e-16), `python python/m1/check_table.py`.

Spec format, report fields and measured scale numbers: `docs/V1.md`. Design and assumptions A1–A11:
`docs/rust-cohomology-design.md`. m1 contract: `docs/SPEC.md` (the gate's anti-cheat rules and hidden holdout are not in
this repo).

## H⁰ tolerance (ADR-015)
One rule: λ counts as zero iff λ ≤ 1e-10·max(1, ‖L₀|_c‖_bound) on its own connected component, or a caller-chosen
absolute permitted delta (`"options": {"h0_absolute_tol": t}`). Reports echo `zero_policy` and `zero_rule`.

## Harmonic extension (`--boundary`)
`sheaf analyze spec.json --boundary boundary.json` with `{"boundary": [{"node": "a", "value": [..]}, ...]}` adds an
`extension` section: the minimiser of ‖δx‖² with the boundary fixed (sheaf-solve), its certified residual and energy,
and `h0_relative` = dim H⁰(X,B;F) under the same H⁰ policy (sheaf-cohom); `unique` iff it is 0. Without `--boundary`
the report is unchanged.
