# sheaf-products: the discrete atoms of the pipeline (ontology: products/laws, LAWS.md).
#   Classify : (Cert, Evidence) -> Verdict     A.classify.*
#   Receipt  : (Verdict, oracles agree) -> Receipt   S.receipt.map
#   Join     : Receipt x Receipt -> Receipt    S.join.monoid, S.join.cover
# Numerics never enter here: Quantize hands over Nat ticks measured relative to the domain tolerance.
# main prints the bridge-A table that the Haskell correlate must reproduce line for line.
# UNCHECKED until the Mac gate runs.
#
#   bend products/bend/main.bend

import Base

# Verdict (domain)
# ----------------

type Verdict is Data:
  Consistent{}
  FixableNodes{}
  LoopConflict{}
  Underdetermined{}
  Inconclusive{}

def nodes(disc_over: Bool) -> Verdict:
  match disc_over:
    case True{}:
      FixableNodes{}
    case False{}:
      Consistent{}

# xonly: no independent edge observations, so circulations are zero by construction
def loops(xonly: Bool, circ_over: Bool, disc_over: Bool) -> Verdict:
  match xonly:
    case True{}:
      nodes(disc_over)
    case False{}:
      match circ_over:
        case True{}:
          LoopConflict{}
        case False{}:
          nodes(disc_over)

# certified: converged AND cut consistent AND NOT ambiguous AND report well-formed (Certify)
# bounded:   False iff a boundary was given and the extension is not unique
# circ:      max fundamental-cycle circulation, in ticks of the loop tolerance (G4: never a share)
def classify(certified: Bool, bounded: Bool, xonly: Bool, circ: Nat, tol_c: Nat, disc: Nat, tol_d: Nat) -> Verdict:
  match certified:
    case False{}:
      Inconclusive{}
    case True{}:
      match bounded:
        case False{}:
          Underdetermined{}
        case True{}:
          loops(xonly, Nat.is_ge(circ, tol_c), Nat.is_ge(disc, tol_d))

# Receipt (Touchstone), ordered HOLDS < BROKEN < BLOCKED < HOLE
# ------------------------------------------------------------

type Receipt is Data:
  Holds{}
  Broken{}
  Blocked{}
  Hole{}

# agree: the two oracles (core Hodge split, product Circulate) agree; disagreement is a HOLE, never averaged
def to_receipt(agree: Bool, v: Verdict) -> Receipt:
  match agree:
    case False{}:
      Hole{}
    case True{}:
      match v:
        case Consistent{}:
          Holds{}
        case FixableNodes{}:
          Broken{}
        case LoopConflict{}:
          Broken{}
        case Underdetermined{}:
          Blocked{}
        case Inconclusive{}:
          Blocked{}

def join(a: Receipt, b: Receipt) -> Receipt:
  match a:
    case Holds{}:
      b
    case Hole{}:
      Hole{}
    case Broken{}:
      match b:
        case Holds{}:
          Broken{}
        case Broken{}:
          Broken{}
        case Blocked{}:
          Blocked{}
        case Hole{}:
          Hole{}
    case Blocked{}:
      match b:
        case Holds{}:
          Blocked{}
        case Broken{}:
          Blocked{}
        case Blocked{}:
          Blocked{}
        case Hole{}:
          Hole{}

# a cover: one entry per declared component; Missing is a component that produced no receipt
type Cover is Data:
  End{}
  Got{r: Receipt, rest: Cover}
  Missing{rest: Cover}

def join_all(xs: Cover) -> Receipt:
  match xs:
    case End{}:
      Holds{}
    case Got{r, rest}:
      join(r, join_all(rest))
    case Missing{rest}:
      Hole{}

def cat(xs: Cover, ys: Cover) -> Cover:
  match xs:
    case End{}:
      ys
    case Got{r, rest}:
      Got{r, cat(rest, ys)}
    case Missing{rest}:
      Missing{cat(rest, ys)}

# Bridge-A table
# --------------

def Verdict.show(v: Verdict) -> String:
  match v:
    case Consistent{}:
      "consistent"
    case FixableNodes{}:
      "fixable"
    case LoopConflict{}:
      "loop"
    case Underdetermined{}:
      "underdetermined"
    case Inconclusive{}:
      "inconclusive"

def Receipt.show(r: Receipt) -> String:
  match r:
    case Holds{}:
      "HOLDS"
    case Broken{}:
      "BROKEN"
    case Blocked{}:
      "BLOCKED"
    case Hole{}:
      "HOLE"

def Bool.bit(b: Bool) -> String:
  match b:
    case True{}:
      "1"
    case False{}:
      "0"

def main() -> IO(Unit):
  do IO<Unit>:
    IO.print(Verdict.show(classify(False{}, True{}, False{}, 9n, 1n, 9n, 1n)))
    IO.print(Verdict.show(classify(True{}, False{}, False{}, 9n, 1n, 9n, 1n)))
    IO.print(Verdict.show(classify(True{}, True{}, True{}, 9n, 1n, 0n, 1n)))
    IO.print(Verdict.show(classify(True{}, True{}, True{}, 9n, 1n, 2n, 1n)))
    IO.print(Verdict.show(classify(True{}, True{}, False{}, 3n, 3n, 0n, 1n)))
    IO.print(Verdict.show(classify(True{}, True{}, False{}, 2n, 3n, 5n, 5n)))
    IO.print(Verdict.show(classify(True{}, True{}, False{}, 2n, 3n, 4n, 5n)))
    IO.print(Receipt.show(to_receipt(False{}, LoopConflict{})))
    IO.print(Receipt.show(join(Broken{}, Blocked{})))
    IO.print(Receipt.show(join_all(Got{Holds{}, Missing{Got{Holds{}, End{}}}})))
