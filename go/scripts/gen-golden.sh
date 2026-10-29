#!/usr/bin/env bash
# Regenerate go/testdata/golden from the Rust reference (the oracle). Go CI never runs this: it only reads the
# committed goldens. Needs cargo (the workspace toolchain, or RUSTUP_TOOLCHAIN=stable where 1.85.1 is unavailable).
#   go/scripts/gen-golden.sh            # build the oracle, regenerate specs + goldens
# Env: SHEAF_BIN / ORACLE_BIN to use prebuilt binaries.
set -euo pipefail
GO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$(cd "$GO_DIR/.." && pwd)"
cd "$ROOT"
if [ -z "${SHEAF_BIN:-}" ] || [ -z "${ORACLE_BIN:-}" ]; then
  cargo build --release --locked -p sheaf-cli --bin sheaf --example go_oracle
fi
SHEAF="${SHEAF_BIN:-$ROOT/target/release/sheaf}"
ORACLE="${ORACLE_BIN:-$ROOT/target/release/examples/go_oracle}"
cd "$GO_DIR"
# 1. specs: shared fixtures are copied by hand (see README), toys are committed; the rest is generated
go run ./internal/gen testdata/specs >/dev/null
# 2. goldens
rm -f testdata/golden/*.json
read -r -a ARGS < testdata/golden/args.txt
n_rep=0; n_err=0
errjson() { python3 -c 'import json,sys; print(json.dumps({"exit": int(sys.argv[1]), "stderr": open(sys.argv[2]).read()}, indent=1, sort_keys=True))' "$1" "$2"; }
tmp="$(mktemp)"; trap 'rm -f "$tmp"' EXIT
for spec in testdata/specs/*.json; do
  case "$spec" in *.boundary.json) continue;; esac
  name="$(basename "$spec" .json)"
  bnd="testdata/specs/$name.boundary.json"
  extra=(); [ -f "$bnd" ] && extra=(--boundary "$bnd")
  case "$name" in errf_*) continue;; esac
  set +e; "$SHEAF" analyze "$spec" "${ARGS[@]}" "${extra[@]}" --out "testdata/golden/$name.report.json" 2>"$tmp"; rc=$?; set -e
  if [ $rc -ne 0 ]; then
    rm -f "testdata/golden/$name.report.json"; errjson $rc "$tmp" > "testdata/golden/$name.error.json"; n_err=$((n_err+1)); continue
  fi
  if [ -f "$bnd" ]; then "$ORACLE" "$spec" "$bnd" > "testdata/golden/$name.aux.json"; else "$ORACLE" "$spec" > "testdata/golden/$name.aux.json"; fi
  n_rep=$((n_rep+1))
done
while IFS='|' read -r tag args; do
  tag="$(echo "$tag" | xargs)"; [ -z "$tag" ] && continue; case "$tag" in \#*) continue;; esac
  read -r -a fa <<< "$args"
  set +e; "$SHEAF" analyze testdata/specs/errf_flags.json "${fa[@]}" >/dev/null 2>"$tmp"; rc=$?; set -e
  # keep only the first stderr line: a Rust panic adds a backtrace that is not part of the contract
  head -n 1 "$tmp" > "$tmp.1"; mv "$tmp.1" "$tmp"
  errjson $rc "$tmp" > "testdata/golden/errf_$tag.error.json"; n_err=$((n_err+1))
done < testdata/flag_cases.txt
echo "goldens: $n_rep reports (+aux), $n_err error cases -> go/testdata/golden"
