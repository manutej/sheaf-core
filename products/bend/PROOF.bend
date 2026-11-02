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
    case Holds{} Holds{} Holds{}:
      {==}
    case Holds{} Holds{} Broken{}:
      {==}
    case Holds{} Holds{} Blocked{}:
      {==}
    case Holds{} Holds{} Hole{}:
      {==}
    case Holds{} Broken{} Holds{}:
      {==}
    case Holds{} Broken{} Broken{}:
      {==}
    case Holds{} Broken{} Blocked{}:
      {==}
    case Holds{} Broken{} Hole{}:
      {==}
    case Holds{} Blocked{} Holds{}:
      {==}
    case Holds{} Blocked{} Broken{}:
      {==}
    case Holds{} Blocked{} Blocked{}:
      {==}
    case Holds{} Blocked{} Hole{}:
      {==}
    case Holds{} Hole{} Holds{}:
      {==}
    case Holds{} Hole{} Broken{}:
      {==}
    case Holds{} Hole{} Blocked{}:
      {==}
    case Holds{} Hole{} Hole{}:
      {==}
    case Broken{} Holds{} Holds{}:
      {==}
    case Broken{} Holds{} Broken{}:
      {==}
    case Broken{} Holds{} Blocked{}:
      {==}
    case Broken{} Holds{} Hole{}:
      {==}
    case Broken{} Broken{} Holds{}:
      {==}
    case Broken{} Broken{} Broken{}:
      {==}
    case Broken{} Broken{} Blocked{}:
      {==}
    case Broken{} Broken{} Hole{}:
      {==}
    case Broken{} Blocked{} Holds{}:
      {==}
    case Broken{} Blocked{} Broken{}:
      {==}
    case Broken{} Blocked{} Blocked{}:
      {==}
    case Broken{} Blocked{} Hole{}:
      {==}
    case Broken{} Hole{} Holds{}:
      {==}
    case Broken{} Hole{} Broken{}:
      {==}
    case Broken{} Hole{} Blocked{}:
      {==}
    case Broken{} Hole{} Hole{}:
      {==}
    case Blocked{} Holds{} Holds{}:
      {==}
    case Blocked{} Holds{} Broken{}:
      {==}
    case Blocked{} Holds{} Blocked{}:
      {==}
    case Blocked{} Holds{} Hole{}:
      {==}
    case Blocked{} Broken{} Holds{}:
      {==}
    case Blocked{} Broken{} Broken{}:
      {==}
    case Blocked{} Broken{} Blocked{}:
      {==}
    case Blocked{} Broken{} Hole{}:
      {==}
    case Blocked{} Blocked{} Holds{}:
      {==}
    case Blocked{} Blocked{} Broken{}:
      {==}
    case Blocked{} Blocked{} Blocked{}:
      {==}
    case Blocked{} Blocked{} Hole{}:
      {==}
    case Blocked{} Hole{} Holds{}:
      {==}
    case Blocked{} Hole{} Broken{}:
      {==}
    case Blocked{} Hole{} Blocked{}:
      {==}
    case Blocked{} Hole{} Hole{}:
      {==}
    case Hole{} Holds{} Holds{}:
      {==}
    case Hole{} Holds{} Broken{}:
      {==}
    case Hole{} Holds{} Blocked{}:
      {==}
    case Hole{} Holds{} Hole{}:
      {==}
    case Hole{} Broken{} Holds{}:
      {==}
    case Hole{} Broken{} Broken{}:
      {==}
    case Hole{} Broken{} Blocked{}:
      {==}
    case Hole{} Broken{} Hole{}:
      {==}
    case Hole{} Blocked{} Holds{}:
      {==}
    case Hole{} Blocked{} Broken{}:
      {==}
    case Hole{} Blocked{} Blocked{}:
      {==}
    case Hole{} Blocked{} Hole{}:
      {==}
    case Hole{} Hole{} Holds{}:
      {==}
    case Hole{} Hole{} Broken{}:
      {==}
    case Hole{} Hole{} Blocked{}:
      {==}
    case Hole{} Hole{} Hole{}:
      {==}

def Laws.join_comm(a, b):
  match a b:
    case Holds{} Holds{}:
      {==}
    case Holds{} Broken{}:
      {==}
    case Holds{} Blocked{}:
      {==}
    case Holds{} Hole{}:
      {==}
    case Broken{} Holds{}:
      {==}
    case Broken{} Broken{}:
      {==}
    case Broken{} Blocked{}:
      {==}
    case Broken{} Hole{}:
      {==}
    case Blocked{} Holds{}:
      {==}
    case Blocked{} Broken{}:
      {==}
    case Blocked{} Blocked{}:
      {==}
    case Blocked{} Hole{}:
      {==}
    case Hole{} Holds{}:
      {==}
    case Hole{} Broken{}:
      {==}
    case Hole{} Blocked{}:
      {==}
    case Hole{} Hole{}:
      {==}

def Laws.join_idem(a):
  match a:
    case Holds{}:
      {==}
    case Broken{}:
      {==}
    case Blocked{}:
      {==}
    case Hole{}:
      {==}

def Laws.join_unit_r(a):
  match a:
    case Holds{}:
      {==}
    case Broken{}:
      {==}
    case Blocked{}:
      {==}
    case Hole{}:
      {==}

def Laws.join_hole_absorbs(a):
  match a:
    case Holds{}:
      {==}
    case Broken{}:
      {==}
    case Blocked{}:
      {==}
    case Hole{}:
      {==}

def Laws.join_unit_l(a):
  {==}

def join_hole_r(r: Main.Receipt) -> {Main.join(r, Main.Hole{}) == Main.Hole{} : Main.Receipt}:
  Laws.join_hole_absorbs(r)

# induction on xs: End and Missing compute to Hole; Got{r, rest} is join(r, J) with J == Hole by the
# hypothesis, turned around by Equal.sym so the rewrite replaces J with Hole
def Laws.join_cover_missing_is_hole(xs, ys):
  match xs:
    case End{}:
      {==}
    case Got{r, rest}:
      %Equal.sym(Main.Receipt, Main.join_all(Main.cat(rest, Main.Missing{ys})), Main.Hole{}, Laws.join_cover_missing_is_hole(rest, ys)) : {Main.join(r, _) == Main.Hole{} : Main.Receipt}
      join_hole_r(r)
    case Missing{rest}:
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
