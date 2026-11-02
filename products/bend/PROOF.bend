# UNCHECKED DRAFT — has never been run through `bend`. Agent-owned.
# Expected gate output once runnable: SOME PROOFS FAIL (xonly_no_loop and mono_tolerance are open by design:
# they are tickets T1a/T1b, the first work items for the Opus and Fable arms).
import Base
import ./LAWS.bend as Laws

# Both reduce by computation: classify matches on the literal first arguments.
def Laws.sound_uncertified(b, x, h, th, d, td):
  {==}

def Laws.undetermined(x, h, th, d, td):
  {==}
