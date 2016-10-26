> module FROverZp where

Functional Reconstruction over finite field Z_p

> import Data.Ratio
> import Data.Maybe
> import Data.Numbers.Primes
> import Data.List (null)
>
> import Ffield (modp, inversep, bigPrimes, recCRT, recCRT')
> import Univariate ((^-), stirlingC, fall2pol, npol2pol)

Univariate Polynomial case
Our target is a univariate polynomial
  f :: (Integral a) =>
       Ratio a -> Ratio a -- Real?

> -- Function-modular.
> fmodp :: Integral c => (a -> Ratio c) -> c -> a -> c
> f `fmodp` p = (`modp` p) . f

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> let fs = map f [0..]
  *FROverZp> take 5 $ map (f `fmodp` 101) [0..]
  [34,93,87,16,82]
  *FROverZp> take 5 $ map (`modp` 101) fs
  [34,93,87,16,82]

Difference analysis over Z_p
Every arithmetic should be on Z_p, i.e., (`mod` p).

> accessibleData :: (Ratio Int -> Ratio Int) -> Int -> [Int]
> accessibleData f p = take p $ map (f `fmodp` p) [0..]
> 
> accessibleData' :: [Ratio Int] -> Int -> [Int]
> accessibleData' fs p = take p $ map (`modp` p) fs
>
> difsp :: Integral b => b -> [b] -> [b]
> difsp p xs = map (`mod` p) (zipWith (-) (tail xs) xs)

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> take 5 $ accessibleData f 101
  [34,93,87,16,82]
  *FROverZp> difsp 101 it
  [59,95,30,66]
  *FROverZp> difsp 101 it
  [36,36,36]
  *FROverZp> difsp 101 it
  [0,0]

> difListsp :: Integral b => b -> [[b]] -> [[b]]
> difListsp _ [] = []
> difListsp p xx@(xs:xxs) =
>   if isConst xs then xx
>                 else difListsp p $ difsp p xs : xx
>   where
>     isConst (i:jj@(j:js)) = all (==i) jj
>     isConst _ = error "difListsp: "

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> map head $ difListsp 101 [accessibleData f 101]
  [36,59,34]

Degree, eager and lazy versions

> degreep' p xs = length (difListsp p [xs]) -1
> degreep'Lazy p xs = helper xs 0
>   where
>     helper as@(a:b:c:_) n
>       | a==b && b==c = n -- two times matching
>       | otherwise    = helper (difsp p as) (n+1)
>
> degreep :: Integral b => b -> [b] -> Int
> degreep p xs = let l = degreep'Lazy p xs in
>   degreep' p $ take (l+2) xs

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> let myDeg p = degreep p $ accessibleData f p
  *FROverZp> myDeg 101
  2
  *FROverZp> myDeg 103
  2
  *FROverZp> myDeg 107
  2
  *FROverZp> degreep 101 $ accessibleData (\n -> (1%2)+(2%3)*n+(3%4)*n^2+(6%7)*n^7) 101
  7

> firstDifsp :: Integral a => a -> [a] -> [a]
> firstDifsp p xs = reverse $ map head $ difListsp p [xs']
>   where
>     xs' = take n xs
>     n   = 2+ degreep p xs

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> firstDifsp 101 $ accessibleData f 101
  [34,59,36]
  *FROverZp> firstDifsp 101 $ accessibleData (\n -> (1%2)+(2%3)*n+(3%4)*n^2+(6%7)*n^7) 101
  [51,66,59,33,29,58,32,78]

Our target is this diff-list, since once we reconstruct the diflists from several prime fields to rational field, we can fully convert it to canonical form in Q, by applying Univariate.npol2pol.

> wellOrd :: [[a]] -> [[a]]
> wellOrd xss 
>   | null (head xss) = [] 
>   | otherwise       = map head xss : wellOrd (map tail xss)

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> let fps p = accessibleData f p
  *FROverZp> let ourData p = firstDifsp p (fps p)
  *FROverZp> let fivePrimes = take 5 bigPrimes 
  *FROverZp> map (\p -> zip (ourData p) (repeat p)) fivePrimes 
  [[(299158,897473),(867559,897473),(299160,897473)]
  ,[(299166,897497),(329084,897497),(299168,897497)]
  ,[(598333,897499),(388918,897499),(598335,897499)]
  ,[(598345,897517),(29919,897517),(598347,897517)]
  ,[(299176,897527),(329095,897527),(299178,897527)]
  ]
  *FROverZp> wellOrd it
  [[(299158,897473),(299166,897497),(598333,897499)
   ,(598345,897517),(299176,897527)]
  ,[(867559,897473),(329084,897497),(388918,897499)
   ,(29919,897517),(329095,897527)]
  ,[(299160,897473),(299168,897497),(598335,897499)
   ,(598347,897517),(299178,897527)]
  ]
  *FROverZp> :t it
  it :: [[(Int, Int)]]

We need to transform
  Int -> Integer
to use recCRT :: Integral a => [(a, a)] -> Ratio a

  *FROverZp> let impCasted = 
  [[(299158,897473),(299166,897497),(598333,897499)
   ,(598345,897517),(299176,897527)]
  ,[(867559,897473),(329084,897497),(388918,897499)
   ,(29919,897517),(329095,897527)]
  ,[(299160,897473),(299168,897497),(598335,897499)
   ,(598347,897517),(299178,897527)]
  ]
  *FROverZp> :t impCasted 
  impCasted :: (Num t1, Num t) => [[(t, t1)]]
  *FROverZp> map recCRT impCasted 
  [1 % 3,53 % 30,7 % 3]
  *FROverZp> map recCRT' impCastet 
  [(1 % 3,897473),(53 % 30,897473),(7 % 3,897473)]

This result is consistent:
  *Univariate> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *Univariate> firstDifs (map f [0..10])
  [1 % 3,53 % 30,7 % 3] 

Let us define above casting function

> toInteger2 :: (Integral a1, Integral a) => (a, a1) -> (Integer, Integer)
> toInteger2 (a,b) = (toInteger a, toInteger b)

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> let fps p = accessibleData f p
  *FROverZp> let ourData p = firstDifsp p (fps p)
  *FROverZp> let longList' = map (\p -> zip (ourData p) (repeat p)) bigPrimes 
  *FROverZp> let longList = wellOrd longList' 
  *FROverZp> :t longList
  longList :: [[(Int, Int)]]
  *FROverZp> let longList'' = map (map toInteger2) longList
  *FROverZp> :t longList''
  longList'' :: [[(Integer, Integer)]]
  *FROverZp> map recCRT longList''
  [1 % 3,53 % 30,7 % 3]
  *FROverZp> map recCRT' longList''
  [(1 % 3,897473),(53 % 30,897473),(7 % 3,897473)]
  *FROverZp> let f x = (1%3) + (3%5)*x + (7%6)*x^2
  *FROverZp> let fps p = accessibleData f p
  *FROverZp> let ourData p = firstDifsp p (fps p)
  *FROverZp> let longList' = map (\p -> zip (ourData p) (repeat p)) bigPrimes 
  *FROverZp> let longList = wellOrd longList' 
  *FROverZp> :t longList
  longList :: [[(Int, Int)]]
  *FROverZp> let longList'' = map (map toInteger2) longList
  *FROverZp> :t longList''
  longList'' :: [[(Integer, Integer)]]
  *FROverZp> map recCRT longList''
  [1 % 3,53 % 30,7 % 3]
  *FROverZp> map recCRT' longList''
  [(1 % 3,897473),(53 % 30,897473),(7 % 3,897473)]

Let us try another example:

  *FROverZp> let f x = (895 % 922) + (1080 % 6931)*x + (2323 % 1248)*x^2
  *FROverZp> let fps p = accessibleData f p
  *FROverZp> let longList = map (map toInteger2) $ wellOrd $ map (\p -> zip (firstDifsp p (fps p)) (repeat p)) bigPrimes 
  *FROverZp> map recCRT' longList 
  [(895 % 922,805479325081),(17448553 % 8649888,722916888780872419),(2323 % 624,805479325081)]

This result is consistent to that of on Q:

  *FROverZp> :l Univariate
  [1 of 2] Compiling Polynomials      ( Polynomials.hs, interpreted )
  [2 of 2] Compiling Univariate       ( Univariate.lhs, interpreted )
  Ok, modules loaded: Univariate, Polynomials.
  *Univariate> let f x = (895 % 922) + (1080 % 6931)*x + (2323 % 1248)*x^2
  *Univariate> firstDifs (map f [0..20])
  [895 % 922,17448553 % 8649888,2323 % 624]


