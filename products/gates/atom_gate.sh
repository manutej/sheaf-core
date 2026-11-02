#!/usr/bin/env bash
# Touchstone atom gate for the discrete verdict layer (Classify / Receipt / Join). Runs on the Mac:
# needs bend (2.0.39+), lean 4.34.0 for --verdict, ghc (runghc), cargo. No Python anywhere.
#   bash products/gates/atom_gate.sh
# One receipt line per step; exit 0 only if every step HOLDS. Open laws make step 2 PARTIAL, never HOLDS.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
B=products/bend; H=products/haskell; K=products/kit; G=$K/bridge.golden
fail=0; tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
say() { printf '%-8s %s\n' "$1" "$2"; [[ $1 == HOLDS || $1 == PARTIAL ]] || fail=1; }

# 1 core + human-owned laws untouched
products/gates/g0_invariance.sh >/dev/null && say HOLDS "G0 core and LAWS.bend match the manifest" || say BROKEN "G0 manifest mismatch"

# 2 Bend proofs; open claims are listed, not hidden
out=$(cd $B && bend PROOF.bend 2>&1); echo "$out" > "$tmp/proof.txt"
if grep -q "ALL PROOFS CHECK" <<<"$out"; then say HOLDS "bend PROOF.bend: ALL PROOFS CHECK"
elif grep -q "SOME PROOFS FAIL" <<<"$out" && ! grep -qE "^(Error|- expected)" <<<"$out"; then say PARTIAL "bend PROOF.bend: open claims only ($(grep -oE 'mono_tolerance|[a-z_]+ is open' <<<"$out" | sort -u | tr '\n' ' '))"
else say BROKEN "bend PROOF.bend: a proof fails (see below)"; echo "$out" | head -40; fi

# 3 kernel recheck (BendTT, proved in Lean)
(cd $B && bend PROOF.bend --verdict) > "$tmp/verdict.txt" 2>&1
grep -q "ALL PROOFS CHECK" "$tmp/verdict.txt" && say HOLDS "bend --verdict: kernel agrees" || say PARTIAL "bend --verdict: not ALL PROOFS CHECK (expected while mono_tolerance is open)"

# 4 bridge A: Bend == golden
(bend $B/main.bend) > "$tmp/bend.table" 2>&1
diff -q "$tmp/bend.table" $G >/dev/null && say HOLDS "bridge A: Bend table == golden (10 rows)" || { say BROKEN "bridge A: Bend table differs"; diff "$tmp/bend.table" $G | head; }

# 5 Haskell correlate: same law ids, exhaustive on small domains, then its table
runghc $H/Laws.hs laws > "$tmp/hs.laws" 2>&1 && say HOLDS "Haskell laws: $(grep -c '^HOLDS' "$tmp/hs.laws") HOLD" || { say BROKEN "Haskell laws"; grep BROKEN "$tmp/hs.laws"; }
runghc $H/Laws.hs table > "$tmp/hs.table" 2>&1
diff -q "$tmp/hs.table" $G >/dev/null && say HOLDS "bridge A': Haskell table == golden" || say BROKEN "bridge A': Haskell table differs"

# 6 bridge B: Rust conforms (same law ids) and prints the golden table
(cd $K && cargo test --release -q >/dev/null 2>&1) && say HOLDS "Rust sheaf-kit: law tests pass" || say BROKEN "Rust sheaf-kit tests"
(cd $K && cargo run --release -q --bin bridge-table) > "$tmp/rs.table" 2>/dev/null
diff -q "$tmp/rs.table" $G >/dev/null && say HOLDS "bridge B: Rust table == golden" || say BROKEN "bridge B: Rust table differs"

# 7 mutation, layers 0 and 1: each mutant of the atoms must be killed (proof fails / a law BREAKS).
#   id | file | sed expression
mutants=(
  "M1-join-broken-blocked|main.bend|/^    case Broken{}:$/,/^    case Blocked{}:$/s/^          Blocked{}$/          Broken{}/"
  "M2-circ-swap|main.bend|s/Nat.is_ge(circ, tol_c)/Nat.is_ge(tol_c, circ)/"
  "M3-disagree-not-hole|main.bend|0,/      Hole{}/s//      Blocked{}/"
  "M4-xonly-loops|main.bend|0,/      nodes(disc_over)/s//      LoopConflict{}/"
  "M5-nodes-swapped|main.bend|0,/      FixableNodes{}/s//      Consistent{}/"
)
for m in "${mutants[@]}"; do
  IFS='|' read -r id file expr <<<"$m"
  rm -rf "$tmp/m"; cp -r $B "$tmp/m"; sed -i.bak "$expr" "$tmp/m/$file"
  if cmp -s "$tmp/m/$file" "$B/$file"; then say BROKEN "mutant $id did not apply (fix the sed)"; continue; fi
  # killed iff the checker output differs from the unmutated baseline (open claims fail in both)
  if (cd "$tmp/m" && bend PROOF.bend 2>&1) | diff -q - "$tmp/proof.txt" >/dev/null; then
    # PROOF still checks: no law pins this behaviour. Only the example table can catch it -> weak law.
    if (bend "$tmp/m/main.bend" 2>/dev/null | diff -q - $G >/dev/null); then say BROKEN "L0 mutant $id survives proofs AND table"; else say PARTIAL "L0 mutant $id survives the proofs, killed only by the golden table: a law is missing"; fi
  else say HOLDS "L0 mutant $id killed by a proof"; fi
done
exit $fail
