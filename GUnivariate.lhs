GUnivariate.lhs

> module GUnivariate where

> import Data.Ratio
> import Data.Maybe
> import Data.Either
> import Data.List
> import Control.Monad

> type Q = Ratio Int
> type Graph = [(Q,Q)]

--

For general, non-sequential input.

to do:
Need try-catch-throw type exception-handling. 
Control.Monad.Catch

> graph :: (a -> b) -> [a] -> [(a,b)] 
> graph _ []     = []
> graph f (a:as) = (a, f a): graph f as

> graph' :: (a -> b) -> [a] -> [([a], b)]
> graph' _ []     = []
> graph' f (x:xs) = ([x], f x) : graph' f xs

--

For non-sequential input-output.
Assume
  [(x, f x) | x <- domain]
as the given data.

> divDif :: ([Q], Q) -> ([Q], Q) -> ([Q], Q)
> divDif ([x0], f0) ([x1], f1) = ([x0,x1], (f1 - f0) / (x1 - x0))
> divDif (xs  , f0) (xs' , f1) = ((z:xs'), (f1-f0) / range)
>   where
>     a = last xs'
>     z = head xs
>     range = a - z
>
> map' :: (a -> a -> b) -> [a] -> [b]
> map' _ []            = []
> map' _ [a]           = []
> -- map' f (a:bb@(b:bs)) = (f a b) : map' f bb
> map' f as = zipWith f as (tail as)
> 
> aStepOfDivDif :: [([Q], Q)] -> [([Q], Q)]
> aStepOfDivDif = map' divDif
>
> aStepOfDivDif' []            = []
> aStepOfDivDif' [a,b]         = [divDif a b]
> aStepOfDivDif' (a:bb@(b:bs)) = (divDif a b) : aStepOfDivDif' bb

The example from 1st ref. in page 46:
  *Univariate> let ds = [([0%1],1%1), ([1%1],3%1), ([3%1],2%1)] :: [([Q],Q)]
  *Univariate> aStepOfDivDif ds
  [([0 % 1,1 % 1],2 % 1),([1 % 1,3 % 1],(-1) % 2)]
  *Univariate> aStepOfDivDif it
  [([0 % 1,1 % 1,3 % 1],(-5) % 6)]

> isZeros :: Int -> [([Q], Q)] -> Bool
> isZeros n ds = let ds' = map snd ds in
>   all (==0%1) $ take n ds'
>
> finiteDifferences :: [([Q],Q)] -> [[([Q],Q)]]
> finiteDifferences graph 
>   = if isZeros 3 graph 
>       then [graph]
>       else graph : finiteDifferences g'
>         where
>           g' = aStepOfDivDif graph

An Example; x^2 see page 48.
  *Univariate> let ds = [([0%1],0%1), ([1%1],1%1), ([2%1],4%1), ([3%1],9%1), ([4%1], 16%1), ([5%1], 25%1)] :: [([Q],Q)]
  *Univariate> finiteDifferences ds
  [[([0%1],0%1),([1%1],1%1),([2%1],4%1),([3%1],9%1),([4%1],16%1),([5%1],25%1)]
  ,[([0 % 1,1 % 1],1 % 1),([1 % 1,2 % 1],3 % 1),([2 % 1,3 % 1],5 % 1),([3 % 1,4 % 1],7 % 1),([4 % 1,5 % 1],9 % 1)]
  ,[([0 % 1,1 % 1,2 % 1],1 % 1),([1 % 1,2 % 1,3 % 1],1 % 1),([2 % 1,3 % 1,4 % 1],1 % 1),([3 % 1,4 % 1,5 % 1],1 % 1)]
  ,[([0 % 1,1 % 1,2 % 1,3 % 1],0 % 1),([1 % 1,2 % 1,3 % 1,4 % 1],0 % 1),([2 % 1,3 % 1,4 % 1,5 % 1],0 % 1)]]

  *Univariate> let func x = (1 + x*2 + (1%3)*x^3)
  *Univariate> let ds = graph' func [1,3,4,5,8,9]
  *Univariate> finiteDifferences ds
  [[([1 % 1],10 % 3),([3 % 1],16 % 1),([4 % 1],91 % 3),([5 % 1],158 % 3),([8 % 1],563 % 3),([9 % 1],262 % 1)]
  ,[([1 % 1,3 % 1],19 % 3),([3 % 1,4 % 1],43 % 3),([4 % 1,5 % 1],67 % 3),([5 % 1,8 % 1],45 % 1),([8 % 1,9 % 1],223 % 3)]
  ,[([1 % 1,3 % 1,4 % 1],8 % 3),([3 % 1,4 % 1,5 % 1],4 % 1),([4 % 1,5 % 1,8 % 1],17 % 3),([5 % 1,8 % 1,9 % 1],22 % 3)]
  ,[([1 % 1,3 % 1,4 % 1,5 % 1],1 % 3),([3 % 1,4 % 1,5 % 1,8 % 1],1 % 3),([4 % 1,5 % 1,8 % 1,9 % 1],1 % 3)]
  ,[([1 % 1,3 % 1,4 % 1,5 % 1,8 % 1],0 % 1),([3 % 1,4 % 1,5 % 1,8 % 1,9 % 1],0 % 1)]
  ]

--

> data Diff = Diff Int (Q,Q) Q
>   deriving (Show, Read)

> toDiff :: (Q,Q) -> Diff
> toDiff (x,fx) = Diff 0 (x,x) fx

> diffStep :: Diff -> Diff -> Either String Diff
> diffStep (Diff n (x,x') f) (Diff m (y,y') g)
>   | n == m    = Right $ Diff (n+1) (x,y') ((f - g) / (x - y'))
>   | otherwise = Left "diffStep: mixed degree"

  *GUnivariate> let gs = map (\x -> (x,x^2)) [1,3,4,7,8] :: Graph
  *GUnivariate> gs
  [(1 % 1,1 % 1),(3 % 1,9 % 1),(4 % 1,16 % 1),(7 % 1,49 % 1),(8 % 1,64 % 1)]
  *GUnivariate> map toDiff gs
  [Diff 0 (1 % 1,1 % 1) (1 % 1),Diff 0 (3 % 1,3 % 1) (9 % 1),Diff 0 (4 % 1,4 % 1) (16 % 1),Diff 0 (7 % 1,7 % 1) (49 % 1),Diff 0 (8 % 1,8 % 1) (64 % 1)]
  *GUnivariate> map' diffStep it
  [Right (Diff 1 (1 % 1,3 % 1) (4 % 1)),Right (Diff 1 (3 % 1,4 % 1) (7 % 1)),Right (Diff 1 (4 % 1,7 % 1) (11 % 1)),Right (Diff 1 (7 % 1,8 % 1) (15 % 1))]
  *GUnivariate> sequence it
  Right [Diff 1 (1 % 1,3 % 1) (4 % 1),Diff 1 (3 % 1,4 % 1) (7 % 1),Diff 1 (4 % 1,7 % 1) (11 % 1),Diff 1 (7 % 1,8 % 1) (15 % 1)]
  *GUnivariate> fmap (map' diffStep) it
  Right [Right (Diff 2 (1 % 1,4 % 1) (1 % 1)),Right (Diff 2 (3 % 1,7 % 1) (1 % 1)),Right (Diff 2 (4 % 1,8 % 1) (1 % 1))]
  *GUnivariate> fmap sequence it
  Right (Right [Diff 2 (1 % 1,4 % 1) (1 % 1),Diff 2 (3 % 1,7 % 1) (1 % 1),Diff 2 (4 % 1,8 % 1) (1 % 1)])
  *GUnivariate> join it
  Right [Diff 2 (1 % 1,4 % 1) (1 % 1),Diff 2 (3 % 1,7 % 1) (1 % 1),Diff 2 (4 % 1,8 % 1) (1 % 1)]
  
  *GUnivariate> let gs = map (\x -> (x,x^2)) [1,3,4,7,8] :: Graph
  *GUnivariate> let d0 = map toDiff gs
  *GUnivariate> d0
  [Diff 0 (1 % 1,1 % 1) (1 % 1),Diff 0 (3 % 1,3 % 1) (9 % 1),Diff 0 (4 % 1,4 % 1) (16 % 1),Diff 0 (7 % 1,7 % 1) (49 % 1),Diff 0 (8 % 1,8 % 1) (64 % 1)]
  *GUnivariate> let d1 = sequence . map' diffStep $ d0
  *GUnivariate> d1
  Right [Diff 1 (1 % 1,3 % 1) (4 % 1),Diff 1 (3 % 1,4 % 1) (7 % 1),Diff 1 (4 % 1,7 % 1) (11 % 1),Diff 1 (7 % 1,8 % 1) (15 % 1)]
  *GUnivariate> join . fmap (sequence . map' diffStep) $ d1
  Right [Diff 2 (1 % 1,4 % 1) (1 % 1),Diff 2 (3 % 1,7 % 1) (1 % 1),Diff 2 (4 % 1,8 % 1) (1 % 1)]
  *GUnivariate> join . fmap (sequence . map' diffStep) $ it
  Right [Diff 3 (1 % 1,7 % 1) (0 % 1),Diff 3 (3 % 1,8 % 1) (0 % 1)]
  *GUnivariate> join . fmap (sequence . map' diffStep) $ it
  Right [Diff 4 (1 % 1,8 % 1) (0 % 1)]
  *GUnivariate> join . fmap (sequence . map' diffStep) $ it
  Right []

> firstStep = sequence . map' diffStep
> secondStep = join . fmap firstStep

  *GUnivariate> let gs = map (\x -> (x,x^2)) [1,3,4,7,8] :: Graph
  *GUnivariate> gs
  [(1 % 1,1 % 1),(3 % 1,9 % 1),(4 % 1,16 % 1),(7 % 1,49 % 1),(8 % 1,64 % 1)]
  *GUnivariate> let d0 = map toDiff gs
  *GUnivariate> d0
  [Diff 0 (1 % 1,1 % 1) (1 % 1),Diff 0 (3 % 1,3 % 1) (9 % 1),Diff 0 (4 % 1,4 % 1) (16 % 1),Diff 0 (7 % 1,7 % 1) (49 % 1),Diff 0 (8 % 1,8 % 1) (64 % 1)]
  *GUnivariate> firstStep d0
  Right [Diff 1 (1 % 1,3 % 1) (4 % 1),Diff 1 (3 % 1,4 % 1) (7 % 1),Diff 1 (4 % 1,7 % 1) (11 % 1),Diff 1 (7 % 1,8 % 1) (15 % 1)]
  *GUnivariate> secondStep it
  Right [Diff 2 (1 % 1,4 % 1) (1 % 1),Diff 2 (3 % 1,7 % 1) (1 % 1),Diff 2 (4 % 1,8 % 1) (1 % 1)]
  *GUnivariate> secondStep it
  Right [Diff 3 (1 % 1,7 % 1) (0 % 1),Diff 3 (3 % 1,8 % 1) (0 % 1)]

> {-
> -- newDiffs :: Diff -> [Diff] -> [Either String Diff]
> newDiffs _ []     = []
> newDiffs f (g:gs) = f : f' : newDiffs f' gs
>   where
>     f' = diffStep f g -- diffStep :: Diff -> Diff -> Either [Char] Diff
> -}



