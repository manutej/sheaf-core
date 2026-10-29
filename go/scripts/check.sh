#!/usr/bin/env bash
# The gate for the Go port: run this before pushing. Fast (seconds), needs only Go (no Rust, no network once the
# gonum module is cached). Remote CI runs the same steps but is advisory (continue-on-error).
#   go/scripts/check.sh            # gofmt, go vet, go test (golden parity + properties), CGO_ENABLED=0 build
#   CROSS=1 go/scripts/check.sh    # also cross-build the CLI for linux/darwin/windows and print sizes
set -euo pipefail
cd "$(dirname "$0")/.."
echo "== gofmt";  out="$(gofmt -l .)"; if [ -n "$out" ]; then echo "gofmt needed:"; echo "$out"; exit 1; fi
echo "== go vet"; go vet ./...
echo "== go test (golden parity vs the Rust oracle + property tests)"
log="$(mktemp)"; trap 'rm -f "$log"' EXIT
if ! go test -count=1 -v ./... > "$log" 2>&1; then grep -E -- '--- FAIL|_test.go|^FAIL|panic' "$log" | head -60; echo "check: FAILED"; exit 1; fi
awk '/^parity:/{p=1} p && /^(parity:|  )/' "$log"
awk '/^=== RUN/{next} /^--- PASS/{t++} /^    +--- PASS/{s++} END{printf "tests passed: %d top-level, %d subtests (golden reports, error cases, property seeds)\n", t, s}' "$log"
echo "== build (CGO_ENABLED=0)"
mkdir -p bin
CGO_ENABLED=0 go build -trimpath -ldflags='-s -w' -o bin/sheafgo ./cmd/sheafgo
if [ "${CROSS:-0}" = 1 ]; then
  for t in linux/amd64 linux/arm64 darwin/amd64 darwin/arm64 windows/amd64; do
    os="${t%/*}"; arch="${t#*/}"; ext=""; [ "$os" = windows ] && ext=".exe"
    CGO_ENABLED=0 GOOS="$os" GOARCH="$arch" go build -trimpath -ldflags='-s -w' -o "bin/sheafgo-$os-$arch$ext" ./cmd/sheafgo
    printf '%-24s %8d bytes\n' "$os/$arch" "$(wc -c < "bin/sheafgo-$os-$arch$ext")"
  done
fi
echo "check: OK"
