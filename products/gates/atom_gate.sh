#!/usr/bin/env bash
# Touchstone atom gate for the discrete verdict layer (Classify / Judge / Join). Runs on the Mac:
# needs bend (2.0.39+), lean 4.34.0 for --verdict, ghc (runghc), cargo, perl. No Python anywhere.
#   bash products/gates/atom_gate.sh
# Exit 0 = every step HOLDS; 1 = something BROKEN; 3 = BLOCKED (a tool or a precondition is missing).
# There is no PARTIAL: open claims live in products/bend/OPEN.bend, outside the checked file (seat 2 S2-01, S2-11).
# Every step appends a hash-chained line to products/ledger/ledger.jsonl and keeps its full log (S2-10, S2-18).
# On success it writes products/bend/laws.lock (hashes of the checked files); signed_by stays null until a human signs.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
B=products/bend; H=products/haskell; K=products/kit
commit=$(git rev-parse --short HEAD); runs=products/ledger/runs/$commit; mkdir -p "$runs"
ledger=products/ledger/ledger.jsonl; touch "$ledger"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
worst=0  # 0 HOLDS, 3 BLOCKED, 1 BROKEN (BROKEN dominates)
sha() { shasum -a 256 "$@" | cut -d' ' -f1; }
record() {  # verdict step detail logfile
  local v=$1 step=$2 detail=$3 log=${4:-/dev/null}
  printf '%-8s %-22s %s\n' "$v" "$step" "$detail"
  case $v in BROKEN) worst=1;; BLOCKED) [[ $worst == 0 ]] && worst=3;; esac
  local prev; prev=$(tail -n1 "$ledger" | shasum -a 256 | cut -d' ' -f1)
  printf '{"ts":"%s","commit":"%s","step":"%s","verdict":"%s","detail":"%s","log_sha256":"%s","prev_sha256":"%s","host":"%s","builder_family":"claude"}\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$commit" "$step" "$v" "${detail//\"/\'}" "$(sha "$log" 2>/dev/null || echo none)" "$prev" "$(uname -sm)" >> "$ledger"
}
need() { command -v "$1" >/dev/null || { record BLOCKED "tool:$1" "$1 not on PATH"; exit 3; }; }
need bend; need runghc; need cargo; need perl

# 1 core untouched (laws are no longer in this manifest: they are frozen in laws.lock only after they check)
products/gates/g0_invariance.sh > "$runs/g0.log" 2>&1 && record HOLDS G0 "core matches manifest" "$runs/g0.log" || record BROKEN G0 "core manifest mismatch" "$runs/g0.log"

# 2 Bend proofs: must be ALL PROOFS CHECK (OPEN.bend is not imported)
(cd $B && bend PROOF.bend) > "$runs/bend.log" 2>&1
if grep -q "ALL PROOFS CHECK" "$runs/bend.log"; then record HOLDS bend "ALL PROOFS CHECK ($(grep -c '^law ' $B/LAWS.bend) laws)" "$runs/bend.log"
else record BROKEN bend "SOME PROOFS FAIL: $(grep -m3 -E 'Location|fail' "$runs/bend.log" | tr '\n' ' ')" "$runs/bend.log"; fi

# 3 kernel recheck (BendTT, proved in Lean)
(cd $B && bend PROOF.bend --verdict) > "$runs/verdict.log" 2>&1
grep -q "ALL PROOFS CHECK" "$runs/verdict.log" && record HOLDS bend-verdict "kernel agrees" "$runs/verdict.log" || record BROKEN bend-verdict "kernel does not agree" "$runs/verdict.log"

# 4 bridge A (smoke): Bend's 10 rows == golden
bend $B/main.bend > "$runs/bend.table" 2>&1
diff -q "$runs/bend.table" $K/bridge.golden >/dev/null && record HOLDS bridgeA "Bend table == golden (10 rows)" "$runs/bend.table" || record BROKEN bridgeA "Bend table differs" "$runs/bend.table"

# 5 Haskell correlate: same law ids, exhaustive small domain
runghc $H/Laws.hs laws > "$runs/hs.laws" 2>&1 && record HOLDS haskell "$(grep -c '^HOLDS' "$runs/hs.laws") laws HOLD" "$runs/hs.laws" || record BROKEN haskell "$(grep BROKEN "$runs/hs.laws" | tr '\n' ' ')" "$runs/hs.laws"

# 6 bridge B (full domain, 19234 rows): Haskell == Rust == golden
runghc $H/Laws.hs fulltable > "$runs/hs.full" 2>&1
(cd $K && cargo run --release -q --bin bridge-table -- full) > "$runs/rs.full" 2>/dev/null
if diff -q "$runs/hs.full" $K/bridge_full.golden >/dev/null && diff -q "$runs/rs.full" $K/bridge_full.golden >/dev/null; then
  record HOLDS bridgeB "Haskell == Rust == golden, $(wc -l < $K/bridge_full.golden | tr -d ' ') rows" "$runs/hs.full"
else record BROKEN bridgeB "full tables differ (see $runs/hs.full, $runs/rs.full)" "$runs/hs.full"; fi
(cd $K && cargo test --release -q) > "$runs/rust.log" 2>&1 && record HOLDS rust "sheaf-kit law tests pass" "$runs/rust.log" || record BROKEN rust "sheaf-kit tests fail" "$runs/rust.log"

# 7 mutation. Each mutant must be killed by a LAW (not by a table). perl substitutions with unique context (portable).
#   L0 Bend: killed iff the checker no longer prints ALL PROOFS CHECK (valid only because the baseline does).
bend_mutants=(
  'M1 join(Broken,Blocked)->Broken|s/(def join.*?case Broken\{\}:.*?case Blocked\{\}:\n)          Blocked\{\}/${1}          Broken{}/s'
  'M2 circulation compare swapped|s/Nat\.is_ge\(circ, tol_c\)/Nat.is_ge(tol_c, circ)/'
  'M3 disagreement->Blocked|s/(def to_receipt.*?case False\{\}:\n)      Hole\{\}/${1}      Blocked{}/s'
  'M4 x-only->LoopConflict|s/(def loops.*?case True\{\}:\n)      nodes\(disc_over\)/${1}      LoopConflict{}/s'
  'M5 nodes(true)->Consistent|s/(def nodes.*?case True\{\}:\n)      FixableNodes\{\}/${1}      Consistent{}/s'
  'M6 to_receipt(Consistent)->Broken|s/(case Consistent\{\}:\n)          Holds\{\}/${1}          Broken{}/'
  'M7 unbounded->Inconclusive|s/(case False\{\}:\n)          Underdetermined\{\}/${1}          Inconclusive{}/'
)
if grep -q "ALL PROOFS CHECK" "$runs/bend.log"; then
  for m in "${bend_mutants[@]}"; do
    id=${m%%|*}; expr=${m#*|}
    rm -rf "$tmp/m"; cp -r $B "$tmp/m"; perl -0pi -e "$expr" "$tmp/m/main.bend"
    if cmp -s "$tmp/m/main.bend" "$B/main.bend"; then record BROKEN "L0 $id" "mutant did not apply"; continue; fi
    (cd "$tmp/m" && bend PROOF.bend) > "$tmp/m.log" 2>&1
    if grep -q "ALL PROOFS CHECK" "$tmp/m.log"; then record BROKEN "L0 $id" "survived every proof: a law is missing" "$tmp/m.log"
    else record HOLDS "L0 $id" "killed: $(grep -m1 -oE 'Laws\.[a-z_]+' "$tmp/m.log" || echo 'a proof failed')" "$tmp/m.log"; fi
  done
else record BLOCKED "L0" "baseline does not check, so mutant results would be meaningless"; fi

#   L1 Haskell: killed iff `Laws.hs laws` exits non-zero (some law BROKEN).
hs_mutants=(
  'H1 join(Broken,Blocked)->Broken|s/Broken -> Broken; Blocked -> Blocked;/Broken -> Broken; Blocked -> Broken;/'
  'H2 circulation compare swapped|s/loops xonly \(c >= tc\)/loops xonly (tc >= c)/'
  'H3 disagreement->Blocked|s/toReceipt False _ = Hole/toReceipt False _ = Blocked/'
  'H4 x-only->LoopConflict|s/\| xonly = nodes discOver/| xonly = LoopConflict/'
  'H5 nodes(true)->Consistent|s/nodes d = if d then FixableNodes else Consistent/nodes d = if d then Consistent else Consistent/'
  'H6 toReceipt(Consistent)->Broken|s/  Consistent -> Holds/  Consistent -> Broken/'
  'H7 unbounded->Inconclusive|s/\| not bounded = Underdetermined/| not bounded = Inconclusive/'
)
for m in "${hs_mutants[@]}"; do
  id=${m%%|*}; expr=${m#*|}
  cp $H/Laws.hs "$tmp/Laws.hs"; perl -0pi -e "$expr" "$tmp/Laws.hs"
  if cmp -s "$tmp/Laws.hs" $H/Laws.hs; then record BROKEN "L1 $id" "mutant did not apply"; continue; fi
  if runghc "$tmp/Laws.hs" laws > "$tmp/h.log" 2>&1; then record BROKEN "L1 $id" "survived every Haskell law" "$tmp/h.log"
  else record HOLDS "L1 $id" "killed: $(grep -m1 '^BROKEN' "$tmp/h.log" | cut -c8-40)" "$tmp/h.log"; fi
done

# 8 freeze: only now, and only if everything held, write laws.lock (signed_by stays null for the human).
if [[ $worst == 0 ]]; then
  {
    echo "{"
    echo "  \"format\": \"sheaf-products/laws-lock-v1\", \"commit\": \"$commit\", \"bend\": \"$(bend --version 2>&1 | head -1)\","
    for f in main SPEC LAWS PROOF; do echo "  \"$f.bend\": \"$(sha $B/$f.bend)\","; done
    echo "  \"Laws.hs\": \"$(sha $H/Laws.hs)\", \"bridge_full\": \"$(sha $K/bridge_full.golden)\","
    echo "  \"checker_log\": \"$(sha "$runs/bend.log")\", \"verdict_log\": \"$(sha "$runs/verdict.log")\","
    echo "  \"signed_by\": null"
    echo "}"
  } > $B/laws.lock
  record HOLDS freeze "laws.lock written; awaiting a human signed_by" $B/laws.lock
fi
exit $worst
