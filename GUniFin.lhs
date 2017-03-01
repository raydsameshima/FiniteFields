GUniFin.lhs

Non sequential inputs Newton-interpolation with finite fields.
Our target is a function
  f :: Q -> Q
which is equivalent to determine (canonical) coefficients.

> module GUniFin where
> --
> import Data.Ratio
> import Data.Maybe
> import Data.Either
> import Data.List
> import Control.Monad
> --
> import Polynomials
> import Ffield
> --
> type Q = Ratio Int   -- Rational fields
> type Graph = [(Q,Q)] -- [(x, f x) | x <- some range]
> --
> -- f [a,b,c ..] -> [(f a b), (f b c) ..]
> map' :: (a -> a -> b) -> [a] -> [b]
> map' f as = zipWith f as (tail as)
>
> -- To select Z_p valid inputs.
> sample :: Int   -- prime
>        -> Graph -- increasing
>        -> Graph 
> sample p = filter ((< (fromIntegral p)) . fst)
>
> -- To eliminate (1%p) type "fake" infinity.
> check :: Int 
>       -> Graph 
>       -> Graph -- safe data sets
> check p gs = filter (not . isDanger p) gs
>
> -- To detect (1%p) type infinity.
> isDanger 
>   :: Int -- prime 
>   -> (Q,Q) -> Bool
> isDanger p (_, fx) = (d `rem` p) == 0
>    where d = denominator fx
> 
> project :: Int -> (Q,Q) -> (Int, Int)
> project p (x, fx)
>   | denominator x == 1 = (numerator x, fx `modp'` p)
>   | otherwise          = error "project: integer input?"
 
  *> let gs = map (\x -> (x,(2%3) + (4%11)*x + (13%2)*x^8)) 
                         [1,3,4,7,8,11,13,17,19,20,21,25,29,31,33,34,35] :: Graph
  *> sample 23 gs
  [(1 % 1,497 % 66),(3 % 1,2814785 % 66),(4 % 1,14057542 % 33)
  ,(7 % 1,2473099841 % 66),(8 % 1,3598712950 % 33),(11 % 1,8359996387 % 6)
  ,(13 % 1,349948479665 % 66),(17 % 1,2992599942641 % 66)
  ,(19 % 1,7285948545089 % 66),(20 % 1,5491200000262 % 33)
  ,(21 % 1,16226006666417 % 66)
  ]
  *> check 23 it
  [(1 % 1,497 % 66),(3 % 1,2814785 % 66),(4 % 1,14057542 % 33)
  ,(7 % 1,2473099841 % 66),(8 % 1,3598712950 % 33),(11 % 1,8359996387 % 6)
  ,(13 % 1,349948479665 % 66),(17 % 1,2992599942641 % 66),(19 % 1,7285948545089 % 66)
  ,(20 % 1,5491200000262 % 33),(21 % 1,16226006666417 % 66)
  ]
  *> map (project 23) it
  [(1,Just 3),(3,Just 8),(4,Just 8),(7,Just 15),(8,Just 1),(11,Just 3)
  ,(13,Just 17),(17,Just 20),(19,Just 3),(20,Just 10),(21,Just 17)
  ]

>{- 
> -- From Graph to Zp values.
> onZp 
>   :: Int                -- base prime
>   -> Graph
>   -> [(Int, Maybe Int)] -- in-out on Zp value 
> onZp p = map (project p) . check p . sample p

  *GUniFin> let gs = map (\x -> (x,(2%3) + (4%11)*x + (13%2)*x^8)) [1,3,4,7,8,11,13,17,19,20,21,25,29,31,33,34,35] :: Graph
  *GUniFin> onZp 101 gs
  [(1,Just 6),(3,Just 89),(4,Just 63),(7,Just 98),(8,Just 7),(11,Just 75),(13,Just 51),(17,Just 97),(19,Just 64),(20,Just 59),(21,Just 76),(25,Just 72),(29,Just 66),(31,Just 70),(33,Just 68),(34,Just 73),(35,Just 2)]

> -- using record syntax
> data PDiff 
>   = PDiff { points    :: (Int, Int) -- end points
>           , value     :: Maybe Int
>           , basePrime :: Int
>           }
>   | PReci -- for reciprocal differences
>   deriving (Show, Read)
>
> toPDiff 
>   :: Int 
>   -> (Int, Maybe Int) 
>   -> PDiff
> toPDiff p (x,fx) = PDiff (x,x) fx p

  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) 
                         [1,2,4,5,9,10,11] :: Graph 
  *GUniFin> check 11 . sample 11 $ gs
  [(1 % 1,11 % 6),(2 % 1,16 % 3),(4 % 1,55 % 3),(5 % 1,167 % 6),(9 % 1,515 % 6),(10 % 1,316 % 3)]
  *GUniFin> onZp 11 it
  [(1,Just 0),(2,Just 9),(4,Just 0),(5,Just 4),(9,Just 7),(10,Just 10)]
  *GUniFin> map (toPDiff 11) it
  [PDiff {points = (1,1), value = Just 0, basePrime = 11},PDiff {points = (2,2), value = Just 9, basePrime = 11},PDiff {points = (4,4), value = Just 0, basePrime = 11},PDiff {points = (5,5), value = Just 4, basePrime = 11},PDiff {points = (9,9), value = Just 7, basePrime = 11},PDiff {points = (10,10), value = Just 10, basePrime = 11}]

> -- newtonTriangleZp :: Int -> [Diff] -> [[Diff]]
> newtonTriangleZp :: [PDiff] -> [[PDiff]]
> newtonTriangleZp fs
>   | length fs < 3 = []
>   | otherwise = helper [sf3] (drop 3 fs)
>   where
>     sf3 = reverse . take 3 $ fs -- [[f2,f1,f0]]
>     helper fss [] = error "newtonBT: need more evaluation" 
>     helper fss (f:fs)
>       | isConsts 3 . last $ fss = fss
>       | otherwise               = helper (add1 f fss) fs
>
> isConsts :: Int -> [PDiff] -> Bool
> isConsts n ds
>   | length ds < n = False    
> isConsts n ds
>   = all (==l) $ take (n-1) ls
>   where (l:ls) = map value ds
>
> -- backward
> add1 :: PDiff -> [[PDiff]] -> [[PDiff]]
> add1 f [gs] = fgs : [zipWith bdiffStep fgs gs] -- singleton
>   where 
>     fgs = f:gs
> add1 f (gg@(g:gs) : hhs) -- gg is reversed order
>             = (f:gg) : add1 fg hhs
>   where
>     fg = bdiffStep f g
>
> -- backward
> bdiffStep :: PDiff -> PDiff -> PDiff
> bdiffStep (PDiff (y,y') g q) (PDiff (x,x') f p)
>   | p == q    = PDiff (x,y') finiteDiff p
>   | otherwise = error "bdiffStep: different primes?"
>   where
>     finiteDiff = (%) <$> fg <*> xy' >>= (`modp` p)
>     xy' = return . fromIntegral $ (x - y' `mod` p)
>     fg = (`mod` p) <$> ((-) <$> f <*> g)

  *GUniFin> gs
  [(1 % 1,11 % 6),(2 % 1,16 % 3),(4 % 1,55 % 3),(5 % 1,167 % 6),(9 % 1,515 % 6),(10 % 1,316 % 3),(11 % 1,761 % 6)]
  *GUniFin> onZp 101 . check 101 . sample 101 $ gs
  [(1,Just 86),(2,Just 39),(4,Just 52),(5,Just 11),(9,Just 69),(10,Just 38),(11,Just 9)]
  *GUniFin> :t it
  it :: [(Int, Maybe Int)]
  *GUniFin> map (toPDiff 101) it
  [PDiff {points = (1,1), value = Just 86, basePrime = 101},PDiff {points = (2,2), value = Just 39, basePrime = 101},PDiff {points = (4,4), value = Just 52, basePrime = 101},PDiff {points = (5,5), value = Just 11, basePrime = 101},PDiff {points = (9,9), value = Just 69, basePrime = 101},PDiff {points = (10,10), value = Just 38, basePrime = 101},PDiff {points = (11,11), value = Just 9, basePrime = 101}]
  *GUniFin> :t it
  it :: [PDiff]
  *GUniFin> newtonTriangleZp it
  [[PDiff {points = (9,9), value = Just 69, basePrime = 101},PDiff {points = (5,5), value = Just 11, basePrime = 101},PDiff {points = (4,4), value = Just 52, basePrime = 101},PDiff {points = (2,2), value = Just 39, basePrime = 101},PDiff {points = (1,1), value = Just 86, basePrime = 101}],[PDiff {points = (5,9), value = Just 65, basePrime = 101},PDiff {points = (4,5), value = Just 60, basePrime = 101},PDiff {points = (2,4), value = Just 57, basePrime = 101},PDiff {points = (1,2), value = Just 54, basePrime = 101}],[PDiff {points = (4,9), value = Just 1, basePrime = 101},PDiff {points = (2,5), value = Just 1, basePrime = 101},PDiff {points = (1,4), value = Just 1, basePrime = 101}]]
  *GUniFin> map head it
  [PDiff {points = (9,9), value = Just 69, basePrime = 101},PDiff {points = (5,9), value = Just 65, basePrime = 101},PDiff {points = (4,9), value = Just 1, basePrime = 101}]

> graph2Zp :: Int -> Graph -> [(Int, Maybe Int)]
> graph2Zp p = onZp p . check p . sample p 
>
> graph2PDiff :: Int -> Graph -> [PDiff]
> graph2PDiff p = map (toPDiff p) . graph2Zp p
>
> newtonCoeffZp :: Int -> Graph -> [PDiff]
> newtonCoeffZp p = map head . newtonTriangleZp . graph2PDiff p
>
> -}














> {-
> bnewton2canonical :: [Diff] -> [Q] -- using Polynomial.hs
> bnewton2canonical [d]    = [value d]
> bnewton2canonical (d:ds) = (z * next) - (zd .* next) + [value d]
>   where 
>     zd = fst . points $ d
>     next = bnewton2canonical ds
> -}


> {-

Instead of treating Either monad, let me make a working small codes.

> diffStep' :: Diff -> Diff -> Diff
> diffStep' (Diff n (x,x') f) (Diff m (y,y') g) 
>   | n == m = Diff (n+1) (x,y') ((f-g)/(x-y'))
>

> -- trial :: Graph -> [[Diff]]
> trial gs = let sg3 = reverse . take 3 $ gs in
>   if isConstsDiffs' 3 sg3 
>     then [sg3]
>     else sg3 : (trial $ map' diffStep' $ gs)

> plant :: Diff -> [Diff] -> [[Diff]]
> plant d ds
>   | isConstsDiffs' 3 ds = [ds]
>   | otherwise           = dd : [map' diffStep' dd]
>   where dd = (d:ds)

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

> bNewtonCoeff :: Graph -> [Diff]
> bNewtonCoeff = map head . newtonTriangle . map toDiff
>
> bnewton2canonical :: [Diff] -> [Q] -- using Polynomial.hs
> bnewton2canonical [d]    = [value d]
> bnewton2canonical (d:ds) = (z * next) - (zd .* next) + [value d]
>   where 
>     zd = fst . points $ d
>     next = bnewton2canonical ds

--

> -}
