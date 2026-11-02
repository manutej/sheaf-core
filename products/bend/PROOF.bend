# sheaf-products: the proofs. Agent-owned. Finite laws are closed by case analysis ({==} after the
# match fixes every constructor); the cover law by induction with one lemma.
# UNCHECKED until the Mac gate runs.

import Base
import ./main.bend as Main
import ./SPEC.bend as Spec
import ./LAWS.bend as Laws

def Laws.sound_uncertified(b, x, c, tc, d, td):
  {==}

def Laws.undetermined(x, c, tc, d, td):
  {==}

# nodes(b) is never LoopConflict: rewrite the clash through Spec.not_loop, which is Empty only there
def nodes_not_loop(b: Bool) -> {Main.nodes(b) != Main.LoopConflict{} : Main.Verdict}:
  match b:
    case True{}:
      e =>
        %e : Spec.not_loop(_)
        Unit{}
    case False{}:
      e =>
        %e : Spec.not_loop(_)
        Unit{}

# classify(T, T, T, ..) computes to nodes(is_ge(d, td))
def Laws.xonly_no_loop(c, tc, d, td):
  nodes_not_loop(Nat.is_ge(d, td))

# both sides compute to loops(x, is_ge(ci, tc), is_ge(d, td)); rewrite the right one with e
def Laws.circ_only(x, c1, c2, tc, d, td, e):
  %e : {Main.loops(x, Nat.is_ge(c1, tc), Nat.is_ge(d, td)) == Main.loops(x, _, Nat.is_ge(d, td)) : Main.Verdict}
  {==}

def Laws.receipt_map_hole(v):
  {==}

def Laws.receipt_map_uncertified(b, x, c, tc, d, td):
  {==}

def Laws.join_assoc(a, b, c):
  match a b c:
    case Main.Holds{} Main.Holds{} Main.Holds{}:
      {==}
    case Main.Holds{} Main.Holds{} Main.Broken{}:
      {==}
    case Main.Holds{} Main.Holds{} Main.Blocked{}:
      {==}
    case Main.Holds{} Main.Holds{} Main.Hole{}:
      {==}
    case Main.Holds{} Main.Broken{} Main.Holds{}:
      {==}
    case Main.Holds{} Main.Broken{} Main.Broken{}:
      {==}
    case Main.Holds{} Main.Broken{} Main.Blocked{}:
      {==}
    case Main.Holds{} Main.Broken{} Main.Hole{}:
      {==}
    case Main.Holds{} Main.Blocked{} Main.Holds{}:
      {==}
    case Main.Holds{} Main.Blocked{} Main.Broken{}:
      {==}
    case Main.Holds{} Main.Blocked{} Main.Blocked{}:
      {==}
    case Main.Holds{} Main.Blocked{} Main.Hole{}:
      {==}
    case Main.Holds{} Main.Hole{} Main.Holds{}:
      {==}
    case Main.Holds{} Main.Hole{} Main.Broken{}:
      {==}
    case Main.Holds{} Main.Hole{} Main.Blocked{}:
      {==}
    case Main.Holds{} Main.Hole{} Main.Hole{}:
      {==}
    case Main.Broken{} Main.Holds{} Main.Holds{}:
      {==}
    case Main.Broken{} Main.Holds{} Main.Broken{}:
      {==}
    case Main.Broken{} Main.Holds{} Main.Blocked{}:
      {==}
    case Main.Broken{} Main.Holds{} Main.Hole{}:
      {==}
    case Main.Broken{} Main.Broken{} Main.Holds{}:
      {==}
    case Main.Broken{} Main.Broken{} Main.Broken{}:
      {==}
    case Main.Broken{} Main.Broken{} Main.Blocked{}:
      {==}
    case Main.Broken{} Main.Broken{} Main.Hole{}:
      {==}
    case Main.Broken{} Main.Blocked{} Main.Holds{}:
      {==}
    case Main.Broken{} Main.Blocked{} Main.Broken{}:
      {==}
    case Main.Broken{} Main.Blocked{} Main.Blocked{}:
      {==}
    case Main.Broken{} Main.Blocked{} Main.Hole{}:
      {==}
    case Main.Broken{} Main.Hole{} Main.Holds{}:
      {==}
    case Main.Broken{} Main.Hole{} Main.Broken{}:
      {==}
    case Main.Broken{} Main.Hole{} Main.Blocked{}:
      {==}
    case Main.Broken{} Main.Hole{} Main.Hole{}:
      {==}
    case Main.Blocked{} Main.Holds{} Main.Holds{}:
      {==}
    case Main.Blocked{} Main.Holds{} Main.Broken{}:
      {==}
    case Main.Blocked{} Main.Holds{} Main.Blocked{}:
      {==}
    case Main.Blocked{} Main.Holds{} Main.Hole{}:
      {==}
    case Main.Blocked{} Main.Broken{} Main.Holds{}:
      {==}
    case Main.Blocked{} Main.Broken{} Main.Broken{}:
      {==}
    case Main.Blocked{} Main.Broken{} Main.Blocked{}:
      {==}
    case Main.Blocked{} Main.Broken{} Main.Hole{}:
      {==}
    case Main.Blocked{} Main.Blocked{} Main.Holds{}:
      {==}
    case Main.Blocked{} Main.Blocked{} Main.Broken{}:
      {==}
    case Main.Blocked{} Main.Blocked{} Main.Blocked{}:
      {==}
    case Main.Blocked{} Main.Blocked{} Main.Hole{}:
      {==}
    case Main.Blocked{} Main.Hole{} Main.Holds{}:
      {==}
    case Main.Blocked{} Main.Hole{} Main.Broken{}:
      {==}
    case Main.Blocked{} Main.Hole{} Main.Blocked{}:
      {==}
    case Main.Blocked{} Main.Hole{} Main.Hole{}:
      {==}
    case Main.Hole{} Main.Holds{} Main.Holds{}:
      {==}
    case Main.Hole{} Main.Holds{} Main.Broken{}:
      {==}
    case Main.Hole{} Main.Holds{} Main.Blocked{}:
      {==}
    case Main.Hole{} Main.Holds{} Main.Hole{}:
      {==}
    case Main.Hole{} Main.Broken{} Main.Holds{}:
      {==}
    case Main.Hole{} Main.Broken{} Main.Broken{}:
      {==}
    case Main.Hole{} Main.Broken{} Main.Blocked{}:
      {==}
    case Main.Hole{} Main.Broken{} Main.Hole{}:
      {==}
    case Main.Hole{} Main.Blocked{} Main.Holds{}:
      {==}
    case Main.Hole{} Main.Blocked{} Main.Broken{}:
      {==}
    case Main.Hole{} Main.Blocked{} Main.Blocked{}:
      {==}
    case Main.Hole{} Main.Blocked{} Main.Hole{}:
      {==}
    case Main.Hole{} Main.Hole{} Main.Holds{}:
      {==}
    case Main.Hole{} Main.Hole{} Main.Broken{}:
      {==}
    case Main.Hole{} Main.Hole{} Main.Blocked{}:
      {==}
    case Main.Hole{} Main.Hole{} Main.Hole{}:
      {==}

def Laws.join_comm(a, b):
  match a b:
    case Main.Holds{} Main.Holds{}:
      {==}
    case Main.Holds{} Main.Broken{}:
      {==}
    case Main.Holds{} Main.Blocked{}:
      {==}
    case Main.Holds{} Main.Hole{}:
      {==}
    case Main.Broken{} Main.Holds{}:
      {==}
    case Main.Broken{} Main.Broken{}:
      {==}
    case Main.Broken{} Main.Blocked{}:
      {==}
    case Main.Broken{} Main.Hole{}:
      {==}
    case Main.Blocked{} Main.Holds{}:
      {==}
    case Main.Blocked{} Main.Broken{}:
      {==}
    case Main.Blocked{} Main.Blocked{}:
      {==}
    case Main.Blocked{} Main.Hole{}:
      {==}
    case Main.Hole{} Main.Holds{}:
      {==}
    case Main.Hole{} Main.Broken{}:
      {==}
    case Main.Hole{} Main.Blocked{}:
      {==}
    case Main.Hole{} Main.Hole{}:
      {==}

def Laws.join_idem(a):
  match a:
    case Main.Holds{}:
      {==}
    case Main.Broken{}:
      {==}
    case Main.Blocked{}:
      {==}
    case Main.Hole{}:
      {==}

def Laws.join_unit_r(a):
  match a:
    case Main.Holds{}:
      {==}
    case Main.Broken{}:
      {==}
    case Main.Blocked{}:
      {==}
    case Main.Hole{}:
      {==}

def Laws.join_hole_absorbs(a):
  match a:
    case Main.Holds{}:
      {==}
    case Main.Broken{}:
      {==}
    case Main.Blocked{}:
      {==}
    case Main.Hole{}:
      {==}

def Laws.join_unit_l(a):
  {==}

def join_hole_r(r: Main.Receipt) -> {Main.join(r, Main.Hole{}) == Main.Hole{} : Main.Receipt}:
  Laws.join_hole_absorbs(r)

# induction on xs: End and Missing compute to Hole; Got{r, rest} is join(r, J) with J == Hole by the
# hypothesis, turned around by Equal.sym so the rewrite replaces J with Hole
def Laws.join_cover_missing_is_hole(xs, ys):
  match xs:
    case Main.End{}:
      {==}
    case Main.Got{r, rest}:
      %Equal.sym(Main.Receipt, Main.join_all(Main.cat(rest, Main.Missing{ys})), Main.Hole{}, Laws.join_cover_missing_is_hole(rest, ys)) : {Main.join(r, _) == Main.Hole{} : Main.Receipt}
      join_hole_r(r)
    case Main.Missing{rest}:
      {==}

# the partition lemmas over arbitrary Bools (a computed value cannot be matched in place, so the
# law passes the computed comparisons in as parameters)
def loop_lemma(x: Bool, co: Bool, dov: Bool) -> {Spec.is_loop(Main.loops(x, co, dov)) == Spec.fires(x, co) : Bool}:
  match x co dov:
    case True{} True{} True{}:
      {==}
    case True{} True{} False{}:
      {==}
    case True{} False{} True{}:
      {==}
    case True{} False{} False{}:
      {==}
    case False{} True{} True{}:
      {==}
    case False{} True{} False{}:
      {==}
    case False{} False{} True{}:
      {==}
    case False{} False{} False{}:
      {==}

def fixable_lemma(x: Bool, co: Bool, dov: Bool) -> {Spec.is_fixable(Main.loops(x, co, dov)) == Bool.and(dov, Bool.not(Spec.fires(x, co))) : Bool}:
  match x co dov:
    case True{} True{} True{}:
      {==}
    case True{} True{} False{}:
      {==}
    case True{} False{} True{}:
      {==}
    case True{} False{} False{}:
      {==}
    case False{} True{} True{}:
      {==}
    case False{} True{} False{}:
      {==}
    case False{} False{} True{}:
      {==}
    case False{} False{} False{}:
      {==}

def Laws.loop_iff_circ(x, c, tc, d, td):
  loop_lemma(x, Nat.is_ge(c, tc), Nat.is_ge(d, td))

def Laws.fixable_iff_disc(x, c, tc, d, td):
  fixable_lemma(x, Nat.is_ge(c, tc), Nat.is_ge(d, td))

def consistent_lemma(x: Bool, co: Bool, dov: Bool) -> {Spec.is_consistent(Main.loops(x, co, dov)) == Bool.and(Bool.not(Spec.fires(x, co)), Bool.not(dov)) : Bool}:
  match x co dov:
    case True{} True{} True{}:
      {==}
    case True{} True{} False{}:
      {==}
    case True{} False{} True{}:
      {==}
    case True{} False{} False{}:
      {==}
    case False{} True{} True{}:
      {==}
    case False{} True{} False{}:
      {==}
    case False{} False{} True{}:
      {==}
    case False{} False{} False{}:
      {==}

def Laws.consistent_iff(x, c, tc, d, td):
  consistent_lemma(x, Nat.is_ge(c, tc), Nat.is_ge(d, td))

def Laws.holds_iff_consistent(v):
  match v:
    case Main.Consistent{}:
      {==}
    case Main.FixableNodes{}:
      {==}
    case Main.LoopConflict{}:
      {==}
    case Main.Underdetermined{}:
      {==}
    case Main.Inconclusive{}:
      {==}
