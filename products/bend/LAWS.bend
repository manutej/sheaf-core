# UNCHECKED DRAFT — has never been run through `bend`. Human-owned: builder agents may not edit this file (gate G0).
# The verdict classifier's contract (SPEC §3.3, VERIFY_LOOP §3).
import Base
import ./verdict.bend as V

# L-SOUND-1/2: an uncertified computation (not converged, cut inconsistent, or ambiguous) never yields a claim.
law sound_uncertified:
  for b: Bool
  for x: Bool
  for h: Nat
  for th: Nat
  for d: Nat
  for td: Nat
  {V.classify(False{}, b, x, h, th, d, td) == V.Inconclusive{} : V.Verdict}

# L-UNDET: certified, boundary given, extension not unique -> underdetermined, whatever the energies.
law undetermined:
  for x: Bool
  for h: Nat
  for th: Nat
  for d: Nat
  for td: Nat
  {V.classify(True{}, False{}, x, h, th, d, td) == V.Underdetermined{} : V.Verdict}

# L-XONLY: with node data only (w := delta x), a loop conflict is impossible, so it is never reported.
law xonly_no_loop:
  for h: Nat
  for th: Nat
  for d: Nat
  for td: Nat
  {V.classify(True{}, True{}, True{}, h, th, d, td) != V.LoopConflict{} : V.Verdict}

# L-MONO: raising the harmonic tolerance never turns a consistent verdict into a loop conflict.
law mono_tolerance:
  for h: Nat
  for th: Nat
  for k: Nat
  for d: Nat
  for td: Nat
  {V.classify(True{}, True{}, False{}, h, th, d, td) == V.Consistent{} : V.Verdict} ->
  {V.classify(True{}, True{}, False{}, h, Nat.add(th, k), d, td) == V.Consistent{} : V.Verdict}
