# python/ — REFERENCE ONLY

Not a runtime dependency of any crate. Kept so the Rust code can be checked against the m1 Python oracle.

- `sheafcore.py`, `models.py`, `sheaf.schema.json`: copies of the gated m1 reference (FRONTIER `sheaf-experiments/contract/`, SPEC v1.18).
- `parity/parity.py`: Python sheafcore vs the Rust v1 core (`cargo run -p sheaf-cli --example parity_dump`) on 5 shared fixtures; tolerance 1e-9, H⁰ exact.
- `m1/make_fixtures.py`: deterministic SYNTHETIC m1 fixtures (`fixtures/m1/`), the generator copy `m1/gen_mutations.py` (paths only changed) and `m1/check_table.py` (mutation table vs pydantic models + JSON Schema).
- `crosscheck.py` / `crosscheck.jl`: the spike's independent numpy / CellularSheaves.jl cross-checks (read `out/fixtures.json` from `cargo run -p sheaf-cohom --example dump`; the fixtures live in `sheaf_cell::fixtures`).

Needs numpy, scipy, pydantic, jsonschema.
