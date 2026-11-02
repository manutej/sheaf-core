#!/usr/bin/env bash
# G1 on Manu's Mac (Apple GPU). Run from the repo root of the verify-loop/products-v0 checkout:
#   bash products/gates/g1_bend_mac.sh            # check only
#   bash products/gates/g1_bend_mac.sh --install  # first time: installs bend via the official installer
# Appends one hash-chained line to products/ledger/ledger.jsonl and exits non-zero unless both checks say ALL PROOFS CHECK.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
if [[ "${1:-}" == "--install" ]]; then curl -fsSL https://bend-lang.com/install.sh | sh; fi
command -v bend >/dev/null || { echo "bend not on PATH (run with --install, then open a new shell)"; exit 2; }

products/gates/g0_invariance.sh || exit 1
mkdir -p products/ledger/g1
ts=$(date -u +%Y-%m-%dT%H:%M:%SZ); commit=$(git rev-parse --short HEAD)
t0=$(date +%s)
(cd products/bend && bend PROOF.bend)            > products/ledger/g1/$commit.check.txt   2>&1; c1=$?
(cd products/bend && bend PROOF.bend --verdict)  > products/ledger/g1/$commit.verdict.txt 2>&1; c2=$?
wall=$(( $(date +%s) - t0 ))
ok1=$(grep -c "ALL PROOFS CHECK" products/ledger/g1/$commit.check.txt)
ok2=$(grep -c "ALL PROOFS CHECK" products/ledger/g1/$commit.verdict.txt)
verdict=$([[ $ok1 -ge 1 && $ok2 -ge 1 ]] && echo PASS || echo FAIL)
art=$(cat products/ledger/g1/$commit.*.txt products/bend/*.bend | shasum -a 256 | cut -d' ' -f1)
prev=$(tail -n1 products/ledger/ledger.jsonl 2>/dev/null | shasum -a 256 | cut -d' ' -f1)
echo "{\"ts\":\"$ts\",\"ticket\":\"${TICKET:-adhoc}\",\"arm\":\"${ARM:-manual}\",\"commit\":\"$commit\",\"gate\":\"G1\",\"verdict\":\"$verdict\",\"bend_exit\":[$c1,$c2],\"artifact_sha256\":\"$art\",\"prev_sha256\":\"$prev\",\"wall_s\":$wall,\"host\":\"$(uname -sm)\"}" >> products/ledger/ledger.jsonl
cat products/ledger/g1/$commit.check.txt; echo "---- verdict kernel ----"; cat products/ledger/g1/$commit.verdict.txt
echo "G1 $verdict"; [[ $verdict == PASS ]]
