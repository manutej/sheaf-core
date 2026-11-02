# sheaf-products: specification helpers for LAWS.bend. Human-owned, like LAWS.bend.
# UNCHECKED until the Mac gate runs `bend PROOF.bend` and `bend PROOF.bend --verdict`.

import Base
import ./main.bend as Main

# Motives for refuting constructor clashes (GUIDE: rewrite e through disc(_), answer Unit{}).
# not_loop(v) is Empty exactly at LoopConflict, so {nodes(b) == LoopConflict} proves Empty.
def not_loop(v: Main.Verdict) -> Type:
  match v:
    case Main.Consistent{}:
      Unit
    case Main.FixableNodes{}:
      Unit
    case Main.LoopConflict{}:
      Empty
    case Main.Underdetermined{}:
      Unit
    case Main.Inconclusive{}:
      Unit

# Observers used to state the exact verdict partition (A.classify.loop-iff, A.classify.disc-iff)
def is_loop(v: Main.Verdict) -> Bool:
  match v:
    case Main.LoopConflict{}:
      True{}
    case Main.Consistent{}:
      False{}
    case Main.FixableNodes{}:
      False{}
    case Main.Underdetermined{}:
      False{}
    case Main.Inconclusive{}:
      False{}

def is_fixable(v: Main.Verdict) -> Bool:
  match v:
    case Main.FixableNodes{}:
      True{}
    case Main.Consistent{}:
      False{}
    case Main.LoopConflict{}:
      False{}
    case Main.Underdetermined{}:
      False{}
    case Main.Inconclusive{}:
      False{}

# the loop fires: observations exist (not x-only) and the max circulation clears its tolerance
def fires(x: Bool, co: Bool) -> Bool:
  Bool.and(Bool.not(x), co)

def is_consistent(v: Main.Verdict) -> Bool:
  match v:
    case Main.Consistent{}:
      True{}
    case Main.FixableNodes{}:
      False{}
    case Main.LoopConflict{}:
      False{}
    case Main.Underdetermined{}:
      False{}
    case Main.Inconclusive{}:
      False{}

def is_holds(r: Main.Receipt) -> Bool:
  match r:
    case Main.Holds{}:
      True{}
    case Main.Broken{}:
      False{}
    case Main.Blocked{}:
      False{}
    case Main.Hole{}:
      False{}
