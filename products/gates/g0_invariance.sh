#!/usr/bin/env bash
# G0: the core (crates/, lockfile, toolchain pin) and the human-owned laws must be byte-identical to the manifest.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
sha256sum -c --quiet products/gates/core.sha256 && echo "G0 PASS" || { echo "G0 FAIL"; exit 1; }
