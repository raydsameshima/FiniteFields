GUnivariate.lhs

> module GUnivariate where

> import Data.Ratio
> import Data.Maybe
> import Data.Either
> import Data.List
> import Control.Monad

> import Polynomials
> import Ffield

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
> -- f [a,b,c ..] -> [(f a b), (f b c) ..]
> map' :: (a -> a -> b) -> [a] -> [b]
> -- map' _ []            = []
> -- map' _ [a]           = []
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

> -- using record syntax
> data Diff = Diff { rank   :: Int
>                  , points :: (Q,Q)
>                  , value  :: Q
>                  }
>   deriving (Show, Read)

> toDiff :: (Q,Q) -> Diff
> toDiff (x,fx) = Diff 0 (x,x) fx

> diffStep :: Diff -> Diff -> Either String Diff
> diffStep (Diff n (x,x') f) (Diff m (y,y') g)
>   | x == y'   = Left "diffStep: detected an infinity "
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

> zerothStep :: Graph -> Either e [Diff]
> zerothStep = return . map toDiff 
>
> firstStep :: [Diff] -> Either String [Diff]
> firstStep = sequence . map' diffStep
>
> secondStep :: Either String [Diff] -> Either String [Diff]
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

> initialize :: Graph -> [Diff]
> initialize = map toDiff . take 3

> isConstsDiffs :: Int -> [Diff] -> Either String Bool
> isConstsDiffs n ds
>   | n < 0         = Left "isConsts: argument should be positive" 
>   | length ds < n = Left "isConsts: need more data points"
> isConstsDiffs n ds = let (l:ls) = map value ds in
>                      return $ all (==l) $ take (n-1) ls

  *GUnivariate> let gs = map (\x -> (x,x^2)) [1,3,4,7,8,11,13,17] :: Graph
  *GUnivariate> gs
  [(1 % 1,1 % 1),(3 % 1,9 % 1),(4 % 1,16 % 1),(7 % 1,49 % 1),(8 % 1,64 % 1),(11 % 1,121 % 1),(13 % 1,169 % 1),(17 % 1,289 % 1)]
  *GUnivariate> let d0 = map toDiff gs
  *GUnivariate> let d1 = firstStep d0
  *GUnivariate> d1 >>= isConstsDiffs 3
  Right False
  *GUnivariate> let d2 = secondStep d1
  *GUnivariate> d2 >>= isConstsDiffs 3
  Right True
  *GUnivariate> d2
  Right [Diff {rank = 2, points = (1 % 1,4 % 1), value = 1 % 1},Diff {rank = 2, points = (3 % 1,7 % 1), value = 1 % 1},Diff {rank = 2, points = (4 % 1,8 % 1), value = 1 % 1},Diff {rank = 2, points = (7 % 1,11 % 1), value = 1 % 1},Diff {rank = 2, points = (8 % 1,13 % 1), value = 1 % 1},Diff {rank = 2, points = (11 % 1,17 % 1), value = 1 % 1}]
  *GUnivariate> fmap (map value) d2
  Right [1 % 1,1 % 1,1 % 1,1 % 1,1 % 1,1 % 1]

  *GUnivariate> zerothStep gs
  Right [Diff {rank = 0, points = (1 % 1,1 % 1), value = 1 % 1},Diff {rank = 0, points = (3 % 1,3 % 1), value = 9 % 1},Diff {rank = 0, points = (4 % 1,4 % 1), value = 16 % 1},Diff {rank = 0, points = (7 % 1,7 % 1), value = 49 % 1},Diff {rank = 0, points = (8 % 1,8 % 1), value = 64 % 1},Diff {rank = 0, points = (11 % 1,11 % 1), value = 121 % 1},Diff {rank = 0, points = (13 % 1,13 % 1), value = 169 % 1},Diff {rank = 0, points = (17 % 1,17 % 1), value = 289 % 1}]
  *GUnivariate> secondStep it
  Right [Diff {rank = 1, points = (1 % 1,3 % 1), value = 4 % 1},Diff {rank = 1, points = (3 % 1,4 % 1), value = 7 % 1},Diff {rank = 1, points = (4 % 1,7 % 1), value = 11 % 1},Diff {rank = 1, points = (7 % 1,8 % 1), value = 15 % 1},Diff {rank = 1, points = (8 % 1,11 % 1), value = 19 % 1},Diff {rank = 1, points = (11 % 1,13 % 1), value = 24 % 1},Diff {rank = 1, points = (13 % 1,17 % 1), value = 30 % 1}]
  *GUnivariate> zerothStep gs >>= firstStep 
  Right [Diff {rank = 1, points = (1 % 1,3 % 1), value = 4 % 1},Diff {rank = 1, points = (3 % 1,4 % 1), value = 7 % 1},Diff {rank = 1, points = (4 % 1,7 % 1), value = 11 % 1},Diff {rank = 1, points = (7 % 1,8 % 1), value = 15 % 1},Diff {rank = 1, points = (8 % 1,11 % 1), value = 19 % 1},Diff {rank = 1, points = (11 % 1,13 % 1), value = 24 % 1},Diff {rank = 1, points = (13 % 1,17 % 1), value = 30 % 1}]
  *GUnivariate> it >>= isConstsDiffs 3
  Right False
  *GUnivariate> zerothStep gs >>= firstStep >>= firstStep 
  Right [Diff {rank = 2, points = (1 % 1,4 % 1), value = 1 % 1},Diff {rank = 2, points = (3 % 1,7 % 1), value = 1 % 1},Diff {rank = 2, points = (4 % 1,8 % 1), value = 1 % 1},Diff {rank = 2, points = (7 % 1,11 % 1), value = 1 % 1},Diff {rank = 2, points = (8 % 1,13 % 1), value = 1 % 1},Diff {rank = 2, points = (11 % 1,17 % 1), value = 1 % 1}]
  *GUnivariate> it >>= isConstsDiffs 3
  Right True

--

Instead of treating Either monad, let me make a working small codes.

> diffStep' :: Diff -> Diff -> Diff
> diffStep' (Diff n (x,x') f) (Diff m (y,y') g) 
>   | n == m = Diff (n+1) (x,y') ((f-g)/(x-y'))
>
> isConstsDiffs' :: Int -> [Diff] -> Bool
> isConstsDiffs' n ds
>   | length ds < n = False    
> isConstsDiffs' n ds
>   = all (==l) $ take (n-1) ls
>   where (l:ls) = map value ds

> tsrif3 :: Graph -> [Diff]
> tsrif3 = reverse . initialize

> -- trial :: Graph -> [[Diff]]
> trial gs = let sg3 = reverse . take 3 $ gs in
>   if isConstsDiffs' 3 sg3 
>     then [sg3]
>     else sg3 : (trial $ map' diffStep' $ gs)

  *GUnivariate> let gs = map (\x -> (x,x^2)) [1,3,4,7,8,11,13,17] :: Graph
  *GUnivariate> gs
  [(1 % 1,1 % 1),(3 % 1,9 % 1),(4 % 1,16 % 1),(7 % 1,49 % 1),(8 % 1,64 % 1),(11 % 1,121 % 1),(13 % 1,169 % 1),(17 % 1,289 % 1)]
  *GUnivariate> :t trial 
  trial :: [Diff] -> [[Diff]]
  *GUnivariate> trial . map toDiff $ gs
  [[Diff {rank = 0, points = (4 % 1,4 % 1), value = 16 % 1},Diff {rank = 0, points = (3 % 1,3 % 1), value = 9 % 1},Diff {rank = 0, points = (1 % 1,1 % 1), value = 1 % 1}]
  ,[Diff {rank = 1, points = (4 % 1,7 % 1), value = 11 % 1},Diff {rank = 1, points = (3 % 1,4 % 1), value = 7 % 1},Diff {rank = 1, points = (1 % 1,3 % 1), value = 4 % 1}]
  ,[Diff {rank = 2, points = (4 % 1,8 % 1), value = 1 % 1},Diff {rank = 2, points = (3 % 1,7 % 1), value = 1 % 1},Diff {rank = 2, points = (1 % 1,4 % 1), value = 1 % 1}]
  ]

> plant :: Diff -> [Diff] -> [[Diff]]
> plant d ds
>   | isConstsDiffs' 3 ds = [ds]
>   | otherwise           = dd : [map' diffStep' dd]
>   where dd = (d:ds)

  *GUnivariate> let gs = map (\x -> (x,x^2)) [1,3,4,7,8,11,13,17] :: Graph
  *GUnivariate> map toDiff gs
  [Diff {rank = 0, points = (1 % 1,1 % 1), value = 1 % 1},Diff {rank = 0, points = (3 % 1,3 % 1), value = 9 % 1},Diff {rank = 0, points = (4 % 1,4 % 1), value = 16 % 1},Diff {rank = 0, points = (7 % 1,7 % 1), value = 49 % 1},Diff {rank = 0, points = (8 % 1,8 % 1), value = 64 % 1},Diff {rank = 0, points = (11 % 1,11 % 1), value = 121 % 1},Diff {rank = 0, points = (13 % 1,13 % 1), value = 169 % 1},Diff {rank = 0, points = (17 % 1,17 % 1), value = 289 % 1}]
  *GUnivariate> head it
  Diff {rank = 0, points = (1 % 1,1 % 1), value = 1 % 1}
  *GUnivariate> let d0 = head . map toDiff $ gs
  *GUnivariate> d0
  Diff {rank = 0, points = (1 % 1,1 % 1), value = 1 % 1}
  *GUnivariate> let d0s = tail . map toDiff $ gs
  *GUnivariate> d0s
  [Diff {rank = 0, points = (3 % 1,3 % 1), value = 9 % 1},Diff {rank = 0, points = (4 % 1,4 % 1), value = 16 % 1},Diff {rank = 0, points = (7 % 1,7 % 1), value = 49 % 1},Diff {rank = 0, points = (8 % 1,8 % 1), value = 64 % 1},Diff {rank = 0, points = (11 % 1,11 % 1), value = 121 % 1},Diff {rank = 0, points = (13 % 1,13 % 1), value = 169 % 1},Diff {rank = 0, points = (17 % 1,17 % 1), value = 289 % 1}]
  *GUnivariate> plant d0 d0s
  [[Diff {rank = 0, points = (1 % 1,1 % 1), value = 1 % 1},Diff {rank = 0, points = (3 % 1,3 % 1), value = 9 % 1},Diff {rank = 0, points = (4 % 1,4 % 1), value = 16 % 1},Diff {rank = 0, points = (7 % 1,7 % 1), value = 49 % 1},Diff {rank = 0, points = (8 % 1,8 % 1), value = 64 % 1},Diff {rank = 0, points = (11 % 1,11 % 1), value = 121 % 1},Diff {rank = 0, points = (13 % 1,13 % 1), value = 169 % 1},Diff {rank = 0, points = (17 % 1,17 % 1), value = 289 % 1}]
  ,[Diff {rank = 1, points = (1 % 1,3 % 1), value = 4 % 1},Diff {rank = 1, points = (3 % 1,4 % 1), value = 7 % 1},Diff {rank = 1, points = (4 % 1,7 % 1), value = 11 % 1},Diff {rank = 1, points = (7 % 1,8 % 1), value = 15 % 1},Diff {rank = 1, points = (8 % 1,11 % 1), value = 19 % 1},Diff {rank = 1, points = (11 % 1,13 % 1), value = 24 % 1},Diff {rank = 1, points = (13 % 1,17 % 1), value = 30 % 1}]
  ]


--
Newton Backward interpolation

> -- backward
> bdiffStep :: Diff -> Diff -> Diff
> bdiffStep (Diff m (y,y') g) (Diff n (x,x') f)
>   | n == m = Diff (n+1) (x,y') ((f-g)/(x-y'))
>
> -- backward
> add1 :: Diff -> [[Diff]] -> [[Diff]]
> add1 f [gs] = fgs : [zipWith bdiffStep fgs gs] -- singleton
>   where 
>     fgs = f:gs
> add1 f (gg@(g:gs) : hhs) -- gg is reversed order
>             = (f:gg) : add1 fg hhs
>   where
>     fg = bdiffStep f g
>
> -- backward
> stair :: [[Diff]] -> [[Diff]]
> stair [fs] = fs : [zipWith bdiffStep fs (tail fs)]
>
> newtonTriangle :: [Diff] -> [[Diff]]
> newtonTriangle fs
>   | length fs < 3 = []
>   | otherwise = helper [sf3] (drop 3 fs)
>   where
>     sf3 = reverse . take 3 $ fs -- [[f2,f1,f0]]
>     helper fss [] = error "newtonBT: need more evaluation" 
>     helper fss (f:fs)
>       | isConstsDiffs' 3 . last $ fss = fss
>       | otherwise                     = helper (add1 f fss) fs
>
> bNewtonCoeff :: Graph -> [Diff]
> bNewtonCoeff = map head . newtonTriangle . map toDiff
>
> bnewton2canonical :: [Diff] -> [Q] -- using Polynomial.hs
> bnewton2canonical [d]    = [value d]
> bnewton2canonical (d:ds) = (z * next) - (zd .* next) + [value d]
>   where 
>     zd = fst . points $ d
>     next = bnewton2canonical ds
>
> -- for double check
> fNewtonCoeff :: Graph -> [Diff]
> fNewtonCoeff = map last . newtonTriangle . map toDiff

> fnewton2canonical :: [Diff] -> [Q] -- using Polynomial.hs
> fnewton2canonical [d]    = [value d]
> fnewton2canonical (d:ds) = (z * next) - (zd .* next) + [value d]
>   where 
>     zd = snd . points $ d
>     next = fnewton2canonical ds
  
  *GUnivariate> let gs = map (\x -> (x,(2%3) + (4%11)*x + (13%2)*x^8)) [1,3,4,7,8,11,13,17,19,20,21,25,29,31,33,34,35] :: Graph
  *GUnivariate> bnewton2canonical . bNewtonCoeff $ gs
  [2 % 3,4 % 11,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,13 % 2]

--
From now on, we combine finite fields technique to control overflow.

> -- 
> sample :: Int   -- prime
>        -> Graph -- increasing
>        -> Graph 
> sample p = filter ((< (fromIntegral p)) . fst)
>
> -- eliminate (1%p) type "fake" infinity
> check :: Int 
>       -> Graph 
>       -> Graph -- safe data sets
> check p gs = filter (not . isDanger p) gs
>
> isDanger :: Int -- prime 
>        -> (Q,Q) -> Bool
> isDanger p (_, fx) = (d `rem` p) == 0
>    where d = denominator fx
> 
> project :: Int -> (Q,Q) -> (Int, Maybe Int)
> project p (x, fx)
>   | denominator x == 1 = (numerator x, fx `modp` p)
>   | otherwise          = error "project: integer input?"
>

