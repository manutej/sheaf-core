-- sheaf-products: the Haskell correlate (Touchstone step 3, "ENFORCE"). Same law ids as LAWS.bend and
-- products/laws. Base only, no packages:
--   runghc products/haskell/Laws.hs table   -- bridge-A table (must equal `bend products/bend/main.bend`)
--   runghc products/haskell/Laws.hs laws    -- every law, exhaustive over finite types and Nat <= bound
-- Exit code 1 if any law fails. Nat bound is small on purpose: the Bend proofs cover all Nat; this layer
-- is the independent second check that the same ids hold, and the mutation target (layer 1).
module Main (main) where

import Control.Monad (forM_, when)
import Data.IORef
import System.Environment (getArgs)
import System.Exit (exitFailure)

-- the atoms (transcribable core: first-order ADTs and total functions only)
data Verdict = Consistent | FixableNodes | LoopConflict | Underdetermined | Inconclusive
  deriving (Eq, Show, Enum, Bounded)
data Receipt = Holds | Broken | Blocked | Hole
  deriving (Eq, Ord, Show, Enum, Bounded)
data Cover = End | Got Receipt Cover | Missing Cover
  deriving (Eq, Show)

type Nat = Integer -- only non-negative values are ever constructed below

nodes :: Bool -> Verdict
nodes d = if d then FixableNodes else Consistent

loops :: Bool -> Bool -> Bool -> Verdict
loops xonly circOver discOver
  | xonly = nodes discOver
  | circOver = LoopConflict
  | otherwise = nodes discOver

classify :: Bool -> Bool -> Bool -> Nat -> Nat -> Nat -> Nat -> Verdict
classify certified bounded xonly c tc d td
  | not certified = Inconclusive
  | not bounded = Underdetermined
  | otherwise = loops xonly (c >= tc) (d >= td)

toReceipt :: Bool -> Verdict -> Receipt
toReceipt False _ = Hole
toReceipt True v = case v of
  Consistent -> Holds
  FixableNodes -> Broken
  LoopConflict -> Broken
  Underdetermined -> Blocked
  Inconclusive -> Blocked

-- written as the Bend match is written (not as `max`), so mutants of either side are comparable
join :: Receipt -> Receipt -> Receipt
join Holds b = b
join Hole _ = Hole
join Broken b = case b of Holds -> Broken; Broken -> Broken; Blocked -> Blocked; Hole -> Hole
join Blocked b = case b of Hole -> Hole; _ -> Blocked

joinAll :: Cover -> Receipt
joinAll End = Holds
joinAll (Got r rest) = join r (joinAll rest)
joinAll (Missing _) = Hole

cat :: Cover -> Cover -> Cover
cat End ys = ys
cat (Got r rest) ys = Got r (cat rest ys)
cat (Missing rest) ys = Missing (cat rest ys)

-- domains
bools :: [Bool]
bools = [False, True]
receipts :: [Receipt]
receipts = [minBound .. maxBound]
verdicts :: [Verdict]
verdicts = [minBound .. maxBound]
nats :: [Nat]
nats = [0 .. 6]
covers :: Int -> [Cover]
covers 0 = [End]
covers n = End : concat [[Got r c | r <- receipts] ++ [Missing c] | c <- covers (n - 1)]

-- law runner: the id, then every counterexample (none expected)
type Check = (String, [String])

laws :: [Check]
laws =
  [ ("sound_uncertified",
      [ show (b, x, c, tc, d, td) | b <- bools, x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats
      , classify False b x c tc d td /= Inconclusive ])
  , ("undetermined",
      [ show (x, c, tc, d, td) | x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats
      , classify True False x c tc d td /= Underdetermined ])
  , ("xonly_no_loop",
      [ show (c, tc, d, td) | c <- nats, tc <- nats, d <- nats, td <- nats
      , classify True True True c tc d td == LoopConflict ])
  , ("circ_only",
      [ show (x, c1, c2, tc, d, td) | x <- bools, c1 <- nats, c2 <- nats, tc <- nats, d <- nats, td <- nats
      , (c1 >= tc) == (c2 >= tc)
      , classify True True x c1 tc d td /= classify True True x c2 tc d td ])
  , ("mono_tolerance",
      [ show (c, tc, k, d, td) | c <- nats, tc <- nats, k <- nats, d <- nats, td <- nats
      , classify True True False c tc d td == Consistent
      , classify True True False c (tc + k) d td /= Consistent ])
  , ("receipt_map_hole", [ show v | v <- verdicts, toReceipt False v /= Hole ])
  , ("receipt_map_uncertified",
      [ show (b, x, c, tc, d, td) | b <- bools, x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats
      , toReceipt True (classify False b x c tc d td) /= Blocked ])
  , ("join_assoc", [ show (a, b, c) | a <- receipts, b <- receipts, c <- receipts, join a (join b c) /= join (join a b) c ])
  , ("join_comm", [ show (a, b) | a <- receipts, b <- receipts, join a b /= join b a ])
  , ("join_idem", [ show a | a <- receipts, join a a /= a ])
  , ("join_unit_l", [ show a | a <- receipts, join Holds a /= a ])
  , ("join_unit_r", [ show a | a <- receipts, join a Holds /= a ])
  , ("join_hole_absorbs", [ show a | a <- receipts, join a Hole /= Hole ])
  , ("join_cover_missing_is_hole",
      [ show (xs, ys) | xs <- covers 3, ys <- covers 2, joinAll (cat xs (Missing ys)) /= Hole ])
  , ("loop_iff_circ",
      [ show (x, c, tc, d, td) | x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats
      , (classify True True x c tc d td == LoopConflict) /= (not x && c >= tc) ])
  , ("fixable_iff_disc",
      [ show (x, c, tc, d, td) | x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats
      , (classify True True x c tc d td == FixableNodes) /= (d >= td && not (not x && c >= tc)) ])
  , ("consistent_iff",
      [ show (x, c, tc, d, td) | x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats
      , (classify True True x c tc d td == Consistent) /= (not (not x && c >= tc) && not (d >= td)) ])
  , ("holds_iff_consistent", [ show v | v <- verdicts, (toReceipt True v == Holds) /= (v == Consistent) ])
    -- the join is the max of the order HOLDS < BROKEN < BLOCKED < HOLE (S.join.monoid, stated directly)
  , ("join_is_max", [ show (a, b) | a <- receipts, b <- receipts, join a b /= max a b ])
  ]

showV :: Verdict -> String
showV v = case v of
  Consistent -> "consistent"; FixableNodes -> "fixable"; LoopConflict -> "loop"
  Underdetermined -> "underdetermined"; Inconclusive -> "inconclusive"

showR :: Receipt -> String
showR r = case r of Holds -> "HOLDS"; Broken -> "BROKEN"; Blocked -> "BLOCKED"; Hole -> "HOLE"

-- bridge-A table: the same ten rows, same order, as main.bend
table :: [String]
table =
  [ showV (classify False True False 9 1 9 1)
  , showV (classify True False False 9 1 9 1)
  , showV (classify True True True 9 1 0 1)
  , showV (classify True True True 9 1 2 1)
  , showV (classify True True False 3 3 0 1)
  , showV (classify True True False 2 3 5 5)
  , showV (classify True True False 2 3 4 5)
  , showR (toReceipt False LoopConflict)
  , showR (join Broken Blocked)
  , showR (joinAll (Got Holds (Missing (Got Holds End))))
  ]

main :: IO ()
main = do
  args <- getArgs
  case args of
    ["table"] -> mapM_ putStrLn table
    ["fulltable"] -> mapM_ putStrLn fullTable
    ["laws"] -> do
      failed <- newIORef False
      forM_ laws $ \(name, cex) ->
        if null cex then putStrLn ("HOLDS  " ++ name)
        else do
          putStrLn ("BROKEN " ++ name ++ "  e.g. " ++ head cex)
          writeIORef failed True
      bad <- readIORef failed
      when bad exitFailure
    _ -> putStrLn "usage: Laws.hs table|laws" >> exitFailure

-- Bridge B, full domain (seat 2 S2-07): every classify input over Bool^3 x {0..6}^4, every toReceipt input, every
-- join input. Must equal `cargo run --bin bridge-table -- full` byte for byte.
fullTable :: [String]
fullTable =
  [ unwords [b c', b bo, b x, show c, show tc, show d, show td, showV (classify c' bo x c tc d td)]
  | c' <- bools, bo <- bools, x <- bools, c <- nats, tc <- nats, d <- nats, td <- nats ]
  ++ [ unwords ["receipt", b a, showV v, showR (toReceipt a v)] | a <- bools, v <- verdicts ]
  ++ [ unwords ["join", showR p, showR q, showR (join p q)] | p <- receipts, q <- receipts ]
  where b t = if t then "1" else "0"
