# sheaf (Go): the graph-sheaf subset

A pure-Go port of the graph subset of sheaf-core's `sheaf analyze`, for embedding in Go hosts (first consumer:
the meta-suite v2 kernel, gate G3). The Rust workspace is the reference and the test oracle. This package is checked
against golden reports that Rust produced, and nothing here depends on Rust at build or test time. See ADR-018 in
`../DECISIONS.md`.

The package is self-contained. It imports nothing from the rest of this repository, contains no code from the
non-subset parts (m1, FFI, Lanczos, CG, mod-p), and its testdata holds only shared, synthetic, random and adversarial
fixtures. It is intended for release under MIT inside the public meta-suite repo. That release (license file and
module path) happens there, not here.

## Before pushing: run the gate

```
make -C go check          # gofmt, go vet, go test (golden parity + properties), CGO_ENABLED=0 build   (~5 s)
make -C go cross          # the same, plus cross-builds for linux/darwin/windows with sizes
```

This local check is the gate. The GitHub Actions `go` job runs the same steps but is advisory
(`continue-on-error: true`) and is not a required check.

## Use

```
CGO_ENABLED=0 go build -o sheafgo ./cmd/sheafgo
sheafgo analyze spec.json [--out report.json] [--k K] [--null-reps R] [--seed S] [--top N] [--method M] [--boundary boundary.json]
```

The CLI takes the same input JSON (`frontier-sheaf/v1`, `frontier-sheaf/boundary-v1`, documented in `../docs/V1.md`),
flags, exit codes (2 for input errors, with one `P_INPUT_01 ...` line on stderr) and report shape
(`frontier-sheaf/report-v1`: same keys, sorted, 12 significant digits) as the Rust `sheaf analyze`.
Library entry points: `sheaf.Parse`, `sheaf.ParseBoundary`, `sheaf.Analyze`, `(*Report).JSON`, and the pieces
`Coboundary`, `Components`, `SpectralComponents`, `Merge`, `H0Basis`, `Hodge`, `RelativeH0`, `HarmonicExtension`.

## What is ported

δ⁰ with per-edge restriction maps (any stalk dims, including 0; multi-edges; both orientations; identity default),
connected components, a per-component dense symmetric eigendecomposition with the ADR-015 cut
(λ ≤ 1e-10·max(1, ‖L₀|_c‖_bound) on the component's own bound, or `h0_absolute_tol`), H⁰/H¹ by counting
(rank-nullity), energy xᵀL₀x, per-edge and per-node discord, the Hodge exact/harmonic split and the minimum-norm
potential (per component, with the rank set by the same cut: never a global cut), the Gaussian and shuffled nulls
(same xorshift stream and seed mixing, so the numbers match), the report flags (`near_zero_ambiguous`,
`cut_consistent_with_zero_tol`, `lower_bound_only`), relative H⁰(X,B;F), and the dense harmonic extension.

## Known, documented differences from the reference

| field | reference | Go | why |
|---|---|---|---|
| `spectral.solver`, `rounds`, `inner_cg`, residual lists | Lanczos above `dense_max` dofs | always `dense-evd`, `inner_cg: null` | Lanczos is not ported. The reported eigenvalues are the same (the k+1 smallest, widened ×4 up to twice while all are zero). |
| `hodge.cg_iters`, `cert_coclosed` | CG solve | 0, and the certificate of the dense split | CG is not ported. `converged` = certificate ≤ 1e-8, as in the reference. |
| Hodge split when the cut is ambiguous or the reference's CG did not certify | CG without the zero cut | split at the policy's rank | Go follows the dense `CochainComplex::hodge`, which the tests use as the oracle for these cases. |
| `extension.route` / `iters` / `residual_rel` | `cg` when L_II is singular or > 400 dofs | `dense`, `dense-pinv` (minimum norm), `trivial` | No CG. The values and energy match. |
| `--k 0` | panics (exit 101) | exit 2, `P_INPUT_01 flag --k: must be >= 1` | Input errors are never panics. |
| JSON syntax / OS error text | serde / Rust std | Go encoding/json / os | The kind (`P_INPUT_01 json:` / `read`) matches. `spec:` messages match byte for byte. |

## Parity tests

- `testdata/specs/`: 57 valid specs. 5 shared (`shared_*`, copied from `python/parity/fixtures`), 7 toys (`toy_*`,
  synthetic math fixtures), 24 seeded random sheaves (`rand_*`, mixed stalk dims 0–3, de 0–3), 21 adversarial cases
  (`adv_*`: two-scale components, near-cut eigenvalues, multi-edges, zero maps, isolated and dim-0 vertices, the
  empty graph, a component above `dense_max` (the reference runs Lanczos), widening and an exhausted-widening lower bound, absolute tolerance, boundaries, delta0 decoupled inside a
  component). There are also 36 error
  cases (`err_*` invalid specs, `errb_*` invalid boundaries, and flag cases from `flag_cases.txt`).
- `testdata/golden/`: for each valid spec, `<name>.report.json` (the Rust `sheaf analyze` output with the arguments
  in `golden/args.txt`) and `<name>.aux.json` (the `go_oracle` example: the full L₀ spectrum, the dense Hodge split,
  the potential, and bases of H⁰ and relative H⁰). For each invalid one, `<name>.error.json` (exit code and stderr).
- Regenerate them with `make -C go golden` (needs cargo; it uses `RUSTUP_TOOLCHAIN` if the pinned toolchain is
  unavailable). Generated specs come from `go run ./internal/gen`. Shared fixtures and toys are committed as they are.

The tolerances are in `sheaf/parity_test.go` and ADR-018. The test prints the worst observed error for each quantity.
