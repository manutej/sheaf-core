# UNCHECKED DRAFT — has never been run through `bend`. See products/VERIFY_LOOP.md §6.
# sheaf-kit verdict classifier: report evidence -> Verdict. Pure, total; compiled to verdict.mjs for the TS layer.
# Quantities are Nat in fixed units (harmonic energy and tolerance in 1e-12 units, ceil'd by the TS quantiser).
import Base

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

# certified: converged AND cut consistent AND not ambiguous (rule 1 of SPEC §3.3, folded by the caller)
# bounded:   false iff a boundary was given and the extension is not unique
# xonly:     no independent edge observations (w := delta x, harmonic part is zero by construction)
def loops(xonly: Bool, harm_over: Bool, disc_over: Bool) -> Verdict:
  match xonly:
    case True{}:
      nodes(disc_over)
    case False{}:
      match harm_over:
        case True{}:
          LoopConflict{}
        case False{}:
          nodes(disc_over)

def classify(certified: Bool, bounded: Bool, xonly: Bool, +harm: Nat, +tol_h: Nat, +disc: Nat, +tol_d: Nat) -> Verdict:
  match certified:
    case False{}:
      Inconclusive{}
    case True{}:
      match bounded:
        case False{}:
          Underdetermined{}
        case True{}:
          loops(xonly, Nat.is_ge(harm, tol_h), Nat.is_ge(disc, tol_d))
