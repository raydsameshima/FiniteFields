GUniFin.lhs

Non sequential inputs Newton-interpolation with finite fields.
Our target is a function
  f :: Q -> Q
which is equivalent to determine (canonical) coefficients.
Accessible input is pairs of in-out, i.e., a (sub) graph of f.

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
> type Graph = [(Q,Q)] -- [(x, f x) | x <- someFinieRange]
> --
> -- f [a,b,c ..] -> [(f a b), (f b c) ..]
> map' :: (a -> a -> b) -> [a] -> [b]
> map' f as = zipWith f as (tail as)
>
> -- To select Z_p valid inputs.
> sample :: Int   -- prime
>        -> Graph -- increasing input
>        -> Graph 
> sample p = filter ((< (fromIntegral p)) . fst)
>
> -- To eliminate (1%p) type "fake" infinity.
> -- After eliminating these, we can freely use modp', primed version.
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
> project p (x, fx) -- for simplicity
>   | denominator x == 1 = (numerator x, fx `modp'` p)
>   | otherwise          = error "project: integer input?"
>
> -- From Graph to Zp (safe) values.
> onZp 
>   :: Int                -- base prime
>   -> Graph
>   -> [(Int, Int)] -- in-out on Zp value 
> onZp p = map (project p) . check p . sample p
>
> -- using record syntax
> data PDiff 
>   = PDiff { points    :: (Int, Int) -- end points
>           , value     :: Int -- Zp value
>           , basePrime :: Int
>           }
>   deriving (Show, Read)
>
> toPDiff 
>   :: Int 
>   -> (Int, Int) -- in and out mod p 
>   -> PDiff
> toPDiff p (x,fx) = PDiff (x,x) fx p
>
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
>     finiteDiff = ((fg % xy') `modp'` p)
>     xy' = (x - y' `mod` p)
>     fg = ((f-g) `mod` p) 
>
> graph2Zp :: Int -> Graph -> [(Int, Int)]
> graph2Zp p = onZp p . check p . sample p 
>
> graph2PDiff :: Int -> Graph -> [PDiff]
> graph2PDiff p = map (toPDiff p) . graph2Zp p
>
> newtonTriangleZp' :: Int -> Graph -> [[PDiff]]
> newtonTriangleZp' p = newtonTriangleZp . graph2PDiff p
> 
> newtonCoeffZp :: Int -> Graph -> [PDiff]
> newtonCoeffZp p = map head . newtonTriangleZp' p

  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) [1,2,4,5,9,10,11] :: Graph 
  *GUniFin> newtonCoeffZp 101 gs
  [PDiff {points = (9,9), value = 69, basePrime = 101},PDiff {points = (5,9), value = 65, basePrime = 101},PDiff {points = (4,9), value = 1, basePrime = 101}]
  *GUniFin> map (\x -> (Just . value $ x, basePrime x)) it
  [(Just 69,101),(Just 65,101),(Just 1,101)]

> n2cZp :: [PDiff] -> ([Int], Int)
> n2cZp graph = (helper graph, p)
>   where
>     p = basePrime . head $ graph
>     helper [d]    = [value d]
>     helper (d:ds) = map (`mod` p) $ ([value d] + (z * next)) - (map (`mod` p) (zd .* next))
>       where 
>         zd = fst . points $ d
>         next = helper ds
>
> format :: ([Int],Int) -> [(Maybe Int, Int)]
> format (as,p) = [(return a,p) | a <- as]
  
  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) [0,2,3,5,7,8,11] :: Graph 
  *GUniFin> newtonCoeffZp 10007 gs
  [PDiff {points = (7,7), value = 8392, basePrime = 10007},PDiff {points = (5,7), value = 5016, basePrime = 10007},PDiff {points = (3,7), value = 1, basePrime = 10007}]
  *GUniFin> n2cZp it
  ([3336,5004,1],10007)
  *GUniFin> format it
  [(Just 3336,10007),(Just 5004,10007),(Just 1,10007)]
  *GUniFin> map guess it
  [Just (1 % 3,10007),Just (1 % 2,10007),Just (1 % 1,10007)]

  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) [0,2,3,5,7,8,11] :: Graph 
  *GUniFin> map guess . format . n2cZp . newtonCoeffZp 10007 $ gs
  [Just (1 % 3,10007),Just (1 % 2,10007),Just (1 % 1,10007)]
  *GUniFin> let gs = map (\x -> (x,x^5 + x^2 + (1%2)*x + 1%3)) [0,2,3,5,7,8,11,13,17,18,19,21,24,28,31,33,34] :: Graph 
  *GUniFin> map guess . format . n2cZp . newtonCoeffZp 10007 $ gs
  [Just (1 % 3,10007),Just (1 % 2,10007),Just (1 % 1,10007),Just (0 % 1,10007),Just (0 % 1,10007),Just (1 % 1,10007)] 

> preTrial gs p = format . n2cZp . newtonCoeffZp p $ gs

  *GUniFin> let gs = map (\x -> (x,x^5 + x^2 + (1%2)*x + 1%3)) [0,2,3,5,7,8,11,13,17,18,19,21,24,28,31,33,34] :: Graph 
  *GUniFin> map reconstruct . transpose . map (preTrial gs) $ bigPrimes 
  [Just (1 % 3),Just (1 % 2),Just (1 % 1),Just (0 % 1),Just (0 % 1),Just (1 % 1)]

Here is "a" final version, the univariate polynomial reconstruction with finite fields.

> uniPolCoeff :: Graph -> Maybe [(Ratio Integer)]
> uniPolCoeff gs = sequence . map reconstruct . transpose . map (preTrial gs) $ bigPrimes

  *GUniFin> let gs = map (\x -> (x,x^5 + x^2 + (1%2)*x + 1%3)) [0,2,3,5,7,8,11,13,17,18,19,21,24,28,31,33,34] :: Graph
  *GUniFin> gs
  [(0 % 1,1 % 3),(2 % 1,112 % 3),(3 % 1,1523 % 6),(5 % 1,18917 % 6),(7 % 1,101159 % 6),(8 % 1,98509 % 3),(11 % 1,967067 % 6),(13 % 1,2228813 % 6),(17 % 1,8520929 % 6),(18 % 1,5669704 % 3),(19 % 1,14858819 % 6),(21 % 1,24507317 % 6),(24 % 1,23889637 % 3),(28 % 1,51633499 % 3),(31 % 1,171780767 % 6),(33 % 1,234818993 % 6),(34 % 1,136309792 % 3)]
  *GUniFin> uniPolCoeff gs
  Just [1 % 3,1 % 2,1 % 1,0 % 1,0 % 1,1 % 1]

--

Non sequential inputs Thiele-interpolation with finite fields.

to do list
Zipper

> type ListZipper a = ([a], [a])
> goForward, goBack :: ListZipper a -> ListZipper a
> goForward (x:xs, ys) = (xs, x:ys)
> goBack (xs, y:ys) = (y:xs, ys)
>
> list2zipper :: [a] -> ListZipper a 
> list2zipper xs = (xs, [])
> zipper2list :: ListZipper a -> [a]
> zipper2list (xs, [])   = xs
> zipper2list (xs, y:ys) = zipper2list (y:xs, ys)

> removeHeads :: ListZipper [a] -> ListZipper [a]
> removeHeads (xs, ys) = (xs, map tail ys)

> add1'
>   :: (Eq a) =>
>      (a -> a -> a) -> a -> ListZipper [a] -> ListZipper [a]
> add1' bo f zs@(((g:gs):ggs), bs)
>   | f /= g = add1' bo fg (ggs, (f:g:gs):bs) 
>   | f == g = removeHeads zs
>   | otherwise = error "???"
>   where fg = f `bo` g

--
Let me start naive rho:

> rho :: Graph -> Int -> [Q]
> rho gs 0 = map snd gs
> rho gs 1 = zipWith (/) xs' fs'
>   where
>     xs' = zipWith (-) xs (tail xs)
>     xs = map fst gs
>     fs' = zipWith (-) fs (tail fs)
>     fs = map snd gs
> rho gs n = zipWith (+) twoAbove oneAbove
>   where
>     twoAbove = zipWith (/) xs' rs'
>     xs' = zipWith (-) xs (drop n xs)
>     xs = map fst gs
>     rs' = zipWith (-) rs (tail rs)
>     rs = rho gs (n-1)
>     oneAbove = tail $ rho gs (n-2)

  *GUniFin> let func x = (1+x+2*x^2)/(3+2*x +(1%4)*x^2)
  *GUniFin> let fs = map (\x -> (x, func x)) [0,1,3,4,6,7,9,10,11,13,14,15,17,19,20] :: Graph 
  *GUniFin> let r = rho fs
  *GUniFin> r 0
  [1 % 3,16 % 21,88 % 45,37 % 15,79 % 24,424 % 117,688 % 165,211 % 48,1016 % 221,1408 % 285,407 % 80,1864 % 357,2384 % 437,424 % 75,821 % 143]
  *GUniFin> r 1
  [7 % 3,315 % 188,45 % 23,80 % 33,936 % 311,6435 % 1756,880 % 199,10608 % 2137,62985 % 10804,4560 % 671,28560 % 3821,156009 % 18260,32775 % 3244,10725 % 943]
  *GUniFin> r 2
  [(-604) % 159,5116 % 405,9458 % 1065,18962 % 2253,75244 % 9171,117388 % 14439,174700 % 21603,243084 % 30151,329516 % 40955,436876 % 54375,559148 % 69659,26491 % 3303,138404 % 17267]
  *GUniFin> r 3
  [900 % 469,585 % 938,(-5805) % 938,(-19323) % 938,(-23418) % 469,(-165867) % 1876,(-295485) % 1876,(-111560) % 469,(-651015) % 1876,(-977265) % 1876,(-199317) % 268,(-278589) % 268]
  *GUniFin> r 4
  [8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1]
