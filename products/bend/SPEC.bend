# sheaf-products: specification helpers for LAWS.bend. Human-owned, like LAWS.bend.
# UNCHECKED until the Mac gate runs `bend PROOF.bend` and `bend PROOF.bend --verdict`.

import Base
import ./main.bend as Main

# Motives for refuting constructor clashes (GUIDE: rewrite e through disc(_), answer Unit{}).
# not_loop(v) is Empty exactly at LoopConflict, so {nodes(b) == LoopConflict} proves Empty.
def not_loop(v: Main.Verdict) -> Type:
  match v:
    case Consistent{}:
      Unit
    case FixableNodes{}:
      Unit
    case LoopConflict{}:
      Empty
    case Underdetermined{}:
      Unit
    case Inconclusive{}:
      Unit

# Observers used to state the exact verdict partition (A.classify.loop-iff, A.classify.disc-iff)
def is_loop(v: Main.Verdict) -> Bool:
  match v:
    case LoopConflict{}:
      True{}
    case Consistent{}:
      False{}
    case FixableNodes{}:
      False{}
    case Underdetermined{}:
      False{}
    case Inconclusive{}:
      False{}

def is_fixable(v: Main.Verdict) -> Bool:
  match v:
    case FixableNodes{}:
      True{}
    case Consistent{}:
      False{}
    case LoopConflict{}:
      False{}
    case Underdetermined{}:
      False{}
    case Inconclusive{}:
      False{}

# the loop fires: observations exist (not x-only) and the max circulation clears its tolerance
def fires(x: Bool, co: Bool) -> Bool:
  Bool.and(Bool.not(x), co)
