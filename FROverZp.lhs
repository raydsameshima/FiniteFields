FROverZp.lhs

> module FROverZp where

Functional Reconstruction over finite field Z_p

> import Data.Ratio
> import Data.Maybe
> import Data.Numbers.Primes
> import Data.List (null)
> import Control.Monad (sequence)
>
> import Ffield (modp)
> -- , inversep, bigPrimes, recCRT, recCRT')
> -- import Univariate (npol2pol, newtonC)

Univariate Polynomial case
Our target is a univariate polynomial
  f :: (Integral a) =>
       Ratio a -> Ratio a -- Real?

> -- Function-modular, now our modp function is wrapped by Maybe.
> fmodp :: (a -> Ratio Int) -> Int -> a -> Maybe Int
> f `fmodp` p = (`modp` p) . f

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%13)*x^2
  *FROverZp> take 10 $ map (f `fmodp` 13)  [0..]
  [Just 9,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing]
  *FROverZp> take 10 $ map (f `fmodp` 19)  [0..]
  [Just 13,Just 8,Just 7,Just 10,Just 17,Just 9,Just 5,Just 5,Just 9,Just 17]

Difference analysis over Z_p
Every arithmetic should be on Z_p, i.e., (`mod` p).

> accessibleData :: (Num a, Enum a) => (a -> Ratio Int) -> Int -> [Maybe Int]
> accessibleData f p = take p $ map (f `fmodp` p) [0..]
> 
> accessibleData' :: [Ratio Int] -> Int -> [Maybe Int]
> accessibleData' fs p = take p $ map (`modp` p) fs

  *FROverZp> let helper x y = (-) <$> x <*> y
  *FROverZp> :type helper 
  helper :: (Applicative f, Num b) => f b -> f b -> f b
  *FROverZp> :t map helper
  map helper :: (Applicative f, Num b) => [f b] -> [f b -> f b]
  *FROverZp> let myDif xs = zipWith helper (tail xs) xs
  *FROverZp> myDif [Just 5,Just 0,Just 2,Just 4,Just 6,Just 1,Just 3]
  [Just (-5),Just 2,Just 2,Just 2,Just (-5),Just 2]
  *FROverZp> map (fmap (`mod` 13)) it
  [Just 8,Just 2,Just 2,Just 2,Just 8,Just 2]
  *FROverZp> map (fmap (`mod` 13)) . myDif $ [Just 5,Just 0,Just 2,Just 4,Just 6,Just 1,Just 3, Nothing, Just 1, Just 2]
  [Just 8,Just 2,Just 2,Just 2,Just 8,Just 2,Nothing,Nothing,Just 1]

> -- difsp :: (Applicative f, Integral b) => b -> [f b] -> [f b]
> difsp :: Applicative f => Int -> [f Int] -> [f Int]
> difsp p = map (fmap (`mod` p)) . difsp' 
>   where
>     difsp' :: (Applicative f, Num b) => [f b] -> [f b]
>     difsp' xs = zipWith helper (tail xs) xs
>     helper :: (Applicative f, Num b) => f b -> f b -> f b
>     helper x y = (-) <$> x <*> y

  *FROverZp> difsp 13 [Just 5,Just 0,Just 2,Just 4,Just 6,Just 1,Just 3, Nothing, Just 1, Just 2]
  [Just 8,Just 2,Just 2,Just 2,Just 8,Just 2,Nothing,Nothing,Just 1]

> difListsp :: (Applicative f, Eq (f Int)) => Int -> [[f Int]] -> [[f Int]]
> difListsp _ [] = []
> difListsp p xx@(xs:xss) =
>   if isConst xs
>     then xx
>     else difListsp p $ difsp p xs : xx
>   where
>     isConst (i:jj@(j:js)) = all (==i) jj
>     isConst _ = error "difListsp: "
  
  *FROverZp> let f x = (1%3) + (3%5)*x + (7%13)*x^2
  *FROverZp> let ds = accessibleData f 101
  *FROverZp> map head $ difListsp 101 [ds]
  [Just 71,Just 26,Just 34]

Degree, eager and lazy versions

> degreep' :: (Applicative f, Eq (f Int)) => Int -> [f Int] -> Int
> degreep' p xs = length (difListsp p [xs]) -1
>
> degreepLazy :: (Applicative f, Num t, Eq (f Int)) => Int -> [f Int] -> t
> degreepLazy p xs = helper xs 0
>   where
>     helper as@(a:b:c:_) n
>       | a==b && b==c = n -- two times matching
>       | otherwise    = helper (difsp p as) (n+1)
>
> degreep :: (Applicative f, Eq (f Int)) => Int -> [f Int] -> Int
> degreep p xs = let l = degreepLazy p xs in
>   degreep' p $ take (l+2) xs

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%13)*x^2
  *FROverZp> degreep 101 $ accessibleData f 101
  2
  *FROverZp> degreep 103 $ accessibleData f 103
  2
  *FROverZp> degreep 107 $ accessibleData f 107
  2
  *FROverZp> degreep 11 $ accessibleData f 11
  2
  *FROverZp> degreep 13 $ accessibleData f 13
  1
  *FROverZp> degreep 17 $ accessibleData f 17
  2

> -- firstDifsp :: Integral a => a -> [a] -> [a]
> firstDifsp :: (Applicative f, Eq (f Int)) => Int -> [f Int] -> [f Int]
> firstDifsp p xs = reverse $ map head $ difListsp p [xs']
>   where
>     xs' = take n xs
>     n   = 2 + degreep p xs

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%13)*x^2
  *FROverZp> let fs p = accessibleData f p
  *FROverZp> firstDifsp 101 (fs 101)
  [Just 34,Just 26,Just 71]
  *FROverZp> firstDifsp 103 (fs 103)
  [Just 69,Just 36,Just 9]
  *FROverZp> firstDifsp 107 (fs 107)
  [Just 36,Just 39,Just 34]

  *FROverZp> map ourData [11,13,17,19,101,103,107]
  [[Just 4,Just 3,Just 7]
  ,[Just 9,Nothing]
  ,[Just 6,Just 15,Just 5]
  ,[Just 13,Just 14,Just 4]
  ,[Just 34,Just 26,Just 71]
  ,[Just 69,Just 36,Just 9]
  ,[Just 36,Just 39,Just 34]
  ]

  *FROverZp> let ourData' = sequence . ourData
  *FROverZp> map ourData' ourPrimes 
  [Just [4,3,7]
  ,Nothing
  ,Just [6,15,5]
  ,Just [13,14,4]
  ,Just [34,26,71]
  ,Just [69,36,9]
  ,Just [36,39,34]
  ]

  *FROverZp> zip (map (sequence . ourData) smallPrimes) smallPrimes 
  [(Just [4,3,7],11)
  ,(Nothing,13)
  ,(Just [6,15,5],17)
  ,(Just [13,14,4],19)
  ,(Just [34,26,71],101)
  ,(Just [69,36,9],103)
  ,(Just [36,39,34],107)
  ]






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
  [(895 % 922,805479325081)
  ,(17448553 % 8649888,722916888780872419)
  ,(2323 % 624,805479325081)
  ]

This result is consistent to that of on Q:

  *FROverZp> :l Univariate
  [1 of 2] Compiling Polynomials      ( Polynomials.hs, interpreted )
  [2 of 2] Compiling Univariate       ( Univariate.lhs, interpreted )
  Ok, modules loaded: Univariate, Polynomials.
  *Univariate> let f x = (895 % 922) + (1080 % 6931)*x + (2323 % 1248)*x^2
  *Univariate> firstDifs (map f [0..20])
  [895 % 922,17448553 % 8649888,2323 % 624]

> {-
> list2firstDifZp' fs = map (recCRT' . map toInteger2) $ wellOrd $ map helper bigPrimes
>   where 
>     helper p = zip (firstDifsp p (accessibleData' fs p)) (repeat p)

  *FROverZp> let f x = (895 % 922) + (1080 % 6931)*x + (2323 % 1248)*x^2
  *FROverZp> let fs = map f [0..]
  *FROverZp> list2firstDifZp' fs
  [(895 % 922,805479325081)
  ,(17448553 % 8649888,722916888780872419)
  ,(2323 % 624,805479325081)
  ]
  *FROverZp> map fst it
  [895 % 922,17448553 % 8649888,2323 % 624]
  *FROverZp> newtonC it
  [895 % 922,17448553 % 8649888,2323 % 1248]
  *FROverZp> npol2pol it
  [895 % 922,1080 % 6931,2323 % 1248]

> list2polZp :: [Ratio Int] -> [Ratio Integer]
> list2polZp = npol2pol . newtonC . (map fst) . list2firstDifZp'

--
Univariate Rational function case
Since thiele2coef uses only (*), (+) and (-) operations, we don't have to do these calculation over prime fields.
So, our target should be rho function (matrix?) calculation.

Reciprocal difference

> rhoZp :: Integral a => [Ratio a] -> a -> Int -> a -> a
> rhoZp fs 0 i p = (fs !! i) `modp` p
> rhoZp fs n i p 
>   | n <= 0     = 0
>   | otherwise = (n*inv + rhoZp fs (n-2) (i+1) p) `mod` p
>   where
>     inv = fromJust inv'
>     inv'= (rhoZp fs (n-1) (i+1) p - rhoZp fs (n-1) i p) `inversep` p
>
> aZp :: Integral a => [Ratio a] -> a -> a -> a
> aZp fs 0 p = head fs `modp` p
> aZp fs n p = (rhoZp fs n 0 p - rhoZp fs (n-2) 0 p) `mod` p
>
> tDegreeZp fs p = helper fs 0 p
>   where
>     helper fs n p
>       | isConst fs' = n
>       | otherwise   = helper fs (n+1) p
>       where
>         fs' = map (rhoZp fs n p) [0..]
>     isConst (i:j:_) = i==j

  *FROverZp> let h t = (3+6*t+18*t^2)%(1+2*t+20*t^2)
  *FROverZp> let hs = map h [0..]
  *FROverZp> take 5 $ map (\n -> rhoZp hs 0 n 101) [0..]
  [3,89,64,8,16]
  *FROverZp> take 5 $ map (\n -> rhoZp hs 1 n 101) [0..]
  [74,4,9,38,65]
  *FROverZp> take 5 $ map (\n -> rhoZp hs 2 n 101) [0..]
  [*** Exception: Maybe.fromJust: Nothing

> -}
