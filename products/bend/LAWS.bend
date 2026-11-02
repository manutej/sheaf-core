# sheaf-products: the laws. Human-owned (Touchstone): a builder agent proposes, a human signs (V.sep.signed).
# Law ids mirror products/laws (the typed ontology); the same ids name the Haskell properties and Rust tests.
# Status: the 16 laws up to fixable_iff_disc were checked ALL PROOFS CHECK by an independent review seat
# (seat 1, Bend 2.0.39 from source, 2026-11-02) after three mechanical fixes; consistent_iff and holds_iff_consistent
# were added afterwards and are UNCHECKED until products/gates/atom_gate.sh runs. Open claims live in OPEN.bend.

import Base
import ./main.bend as Main
import ./SPEC.bend as Spec

# A.classify.uncertified -- an uncertified run never yields a claim
law sound_uncertified:
  for b: Bool
  for x: Bool
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Main.classify(False{}, b, x, c, tc, d, td) == Main.Inconclusive{} : Main.Verdict}

# A.classify.undetermined -- certified, boundary given, not unique
law undetermined:
  for x: Bool
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Main.classify(True{}, False{}, x, c, tc, d, td) == Main.Underdetermined{} : Main.Verdict}

# A.classify.xonly -- node data only can never report a loop conflict
law xonly_no_loop:
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Main.classify(True{}, True{}, True{}, c, tc, d, td) != Main.LoopConflict{} : Main.Verdict}

# A.classify.circ-only (the local form of G4, gauge invariance) -- the loop decision reads only
# whether the max circulation clears its tolerance
law circ_only:
  for x: Bool
  for c1: Nat
  for c2: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  for e: {Nat.is_ge(c1, tc) == Nat.is_ge(c2, tc) : Bool}
  {Main.classify(True{}, True{}, x, c1, tc, d, td) == Main.classify(True{}, True{}, x, c2, tc, d, td) : Main.Verdict}

# S.receipt.map -- oracles disagreeing is a HOLE, whatever the verdict
law receipt_map_hole:
  for v: Main.Verdict
  {Main.to_receipt(False{}, v) == Main.Hole{} : Main.Receipt}

# S.receipt.map + A.classify.uncertified (the formal core of G1) -- uncertified is BLOCKED
law receipt_map_uncertified:
  for b: Bool
  for x: Bool
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Main.to_receipt(True{}, Main.classify(False{}, b, x, c, tc, d, td)) == Main.Blocked{} : Main.Receipt}

# S.join.monoid -- receipts form a bounded join-semilattice; HOLDS is the unit, HOLE absorbs
law join_assoc:
  for a: Main.Receipt
  for b: Main.Receipt
  for c: Main.Receipt
  {Main.join(a, Main.join(b, c)) == Main.join(Main.join(a, b), c) : Main.Receipt}

law join_comm:
  for a: Main.Receipt
  for b: Main.Receipt
  {Main.join(a, b) == Main.join(b, a) : Main.Receipt}

law join_idem:
  for a: Main.Receipt
  {Main.join(a, a) == a : Main.Receipt}

law join_unit_l:
  for a: Main.Receipt
  {Main.join(Main.Holds{}, a) == a : Main.Receipt}

law join_unit_r:
  for a: Main.Receipt
  {Main.join(a, Main.Holds{}) == a : Main.Receipt}

law join_hole_absorbs:
  for a: Main.Receipt
  {Main.join(a, Main.Hole{}) == Main.Hole{} : Main.Receipt}

# S.join.cover -- a missing component anywhere in the cover makes the system receipt a HOLE
law join_cover_missing_is_hole:
  for xs: Main.Cover
  for ys: Main.Cover
  {Main.join_all(Main.cat(xs, Main.Missing{ys})) == Main.Hole{} : Main.Receipt}

# A.classify.loop-iff -- certified and bounded: a loop is reported exactly when it fires
law loop_iff_circ:
  for x: Bool
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Spec.is_loop(Main.classify(True{}, True{}, x, c, tc, d, td)) == Spec.fires(x, Nat.is_ge(c, tc)) : Bool}

# A.classify.disc-iff -- certified and bounded: fixable exactly when discord clears its tolerance and no loop fires
# (added after mutation testing: swapping the discord comparison survived every earlier law)
law fixable_iff_disc:
  for x: Bool
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Spec.is_fixable(Main.classify(True{}, True{}, x, c, tc, d, td)) == Bool.and(Nat.is_ge(d, td), Bool.not(Spec.fires(x, Nat.is_ge(c, tc)))) : Bool}

# A.classify.consistent-iff -- certified and bounded: consistent exactly when no loop fires and discord is below
# tolerance (added after seat 2 showed the Consistent branch was pinned by no law)
law consistent_iff:
  for x: Bool
  for c: Nat
  for tc: Nat
  for d: Nat
  for td: Nat
  {Spec.is_consistent(Main.classify(True{}, True{}, x, c, tc, d, td)) == Bool.and(Bool.not(Spec.fires(x, Nat.is_ge(c, tc))), Bool.not(Nat.is_ge(d, td))) : Bool}

# S.judge.map -- with agreeing oracles, HOLDS exactly for a consistent verdict (seat 2 mutant A: Consistent -> BROKEN
# survived every earlier law)
law holds_iff_consistent:
  for v: Main.Verdict
  {Spec.is_holds(Main.to_receipt(True{}, v)) == Spec.is_consistent(v) : Bool}
