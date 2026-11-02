# sheaf-products: open claims. Not imported by PROOF.bend, so `bend PROOF.bend` can reach ALL PROOFS CHECK while
# these stay open (seat 2 S2-01: an open law must not sit in the frozen file). Each is a ticket.
#
# T1b  A.classify.mono -- a theorem, not a conjecture (seat 1 S1-11): it needs a lemma
#      Nat.is_ge(c, tc) == False  ->  Nat.is_ge(c, Nat.add(tc, k)) == False, by induction on tc and c.

import Base
import ./main.bend as Main

law mono_tolerance:
  for c: Nat
  for tc: Nat
  for k: Nat
  for d: Nat
  for td: Nat
  {Main.classify(True{}, True{}, False{}, c, tc, d, td) == Main.Consistent{} : Main.Verdict} ->
  {Main.classify(True{}, True{}, False{}, c, Nat.add(tc, k), d, td) == Main.Consistent{} : Main.Verdict}
