FROverZp.lhs

> module FROverZp where

Functional Reconstruction over finite field Z_p

> import Data.Ratio
> import Data.Maybe
> import Data.Numbers.Primes
> import Data.List (null, transpose)
> import Control.Monad (sequence)
>
> import Ffield (modp, bigPrimes, reconstruct)
> -- , inversep, bigPrimes, recCRT, recCRT')
> import Univariate (newtonC, firstDifs, list2npol)

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

Since we have chosen 
  *Ffield> last bigPrimes 
  10939
we can put a finite number of outputs as our accessible "row" data.

Difference analysis over Z_p
Every arithmetic should be on Z_p, i.e., (`mod` p).

> accessibleData :: (Num a, Enum a) => 
>                   (a -> Ratio Int) -> Int -> [Maybe Int]
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
> degreep' p xs = length (difListsp p [xs]) - 1
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

> firstDifsp :: (Applicative f, Eq (f Int)) => Int -> [f Int] -> [f Int]
> firstDifsp p xs = reverse $ map head $ difListsp p [xs']
>   where
>     xs' = take n xs
>     n   = 2 + degreep p xs
>
> makeAPair :: [Int] -> [Ratio Int] -> [(Int, Maybe [Int])]
> makeAPair ps fs = zip ps . map (sequence . (\p -> firstDifsp p (fsp p))) $ ps 
>   where
>     fsp = accessibleData' fs

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%13)*x^2
  *FROverZp> let fs = map f [0..]
  *FROverZp> let smallerPrimes = filter isPrime [11, 13, 19, 23, 101, 103, 107]
  *FROverZp> makeAPair smallerPrimes  fs
  [(11,Just [4,3,7])
  ,(13,Nothing)
  ,(19,Just [13,14,4])
  ,(23,Just [8,16,17])
  ,(101,Just [34,26,71])
  ,(103,Just [69,36,9])
  ,(107,Just [36,39,34])
  ]

  *FROverZp>  makeAPair smallerPrimes  fs
  [(11,Just [4,3,7])
  ,(13,Nothing)
  ,(19,Just [13,14,4])
  ,(23,Just [8,16,17])
  ,(101,Just [34,26,71])
  ,(103,Just [69,36,9])
  ,(107,Just [36,39,34])
  ]
  *FROverZp> filter (isJust . snd) it
  [(11,Just [4,3,7])
  ,(19,Just [13,14,4])
  ,(23,Just [8,16,17])
  ,(101,Just [34,26,71])
  ,(103,Just [69,36,9])
  ,(107,Just [36,39,34])
  ]
  *FROverZp> map (\(p, xs) -> (zip (sequence xs) (repeat p))) it
  [[(Just 4,11),(Just 3,11),(Just 7,11)]
  ,[(Just 13,19),(Just 14,19),(Just 4,19)]
  ,[(Just 8,23),(Just 16,23),(Just 17,23)]
  ,[(Just 34,101),(Just 26,101),(Just 71,101)]
  ,[(Just 69,103),(Just 36,103),(Just 9,103)]
  ,[(Just 36,107),(Just 39,107),(Just 34,107)]
  ]
  *FROverZp> transpose it
  [[(Just 4,11),(Just 13,19),(Just 8,23),(Just 34,101),(Just 69,103),(Just 36,107)]
  ,[(Just 3,11),(Just 14,19),(Just 16,23),(Just 26,101),(Just 36,103),(Just 39,107)]
  ,[(Just 7,11),(Just 4,19),(Just 17,23),(Just 71,101),(Just 9,103),(Just 34,107)]
  ]
  *FROverZp> :t it
  it :: [[(Maybe Int, Int)]]
  *FROverZp> :t reconstruct 
  reconstruct :: [(Maybe Int, Int)] -> Maybe (Ratio Integer, Integer)
  *FROverZp> map reconstruct it
  [Just (1 % 3,209),Just (74 % 65,485507),Just (14 % 13,209)]

  *FROverZp> let f x = (1%3) + (3%5)*x + (7%13)*x^2
  *FROverZp> let fs = map f [0..]
  *FROverZp> list2npol fs
  [1 % 3,74 % 65,7 % 13]

> firstDifs' :: [Ratio Int] -> [Maybe (Ratio Integer, Integer)]
> firstDifs' = map reconstruct . transpose . map (\(p,xs) -> (zip (sequence xs) (repeat p))) . 
>              filter (isJust .snd) . makeAPair bigPrimes 

  *FROverZp> let g x = 1%153 + x*(133%122) + (x^2)*(1%199) + (x^3)*(922%855)
  *FROverZp> let gs = map g [0..]
  *FROverZp> fmap (newtonC . map fst) . sequence . firstDifs' $ take 100 gs
  Just [1 % 153,45117911 % 20757690,183763 % 56715,922 % 855]
  *FROverZp> list2npol gs
  [1 % 153,45117911 % 20757690,183763 % 56715,922 % 855]

> -- Eager-version, i.e., input should be finite.
> list2npolp :: [Ratio Int] -> Maybe [Ratio Integer]
> list2npolp = fmap (newtonC . map fst) . sequence . firstDifs' 

  *FROverZp> let g x = 1%153 + x*(133%122) + (x^2)*(1%199) + (x^3)*(922%855)
  *FROverZp> let gs = map g [0..100]
  *FROverZp> list2npol gs
  [1 % 153,45117911 % 20757690,183763 % 56715,922 % 855]
  *FROverZp> list2npolp gs
  Just [1 % 153,45117911 % 20757690,183763 % 56715,922 % 855]
 

--
Univariate Rational function case
Since thiele2coef uses only (*), (+) and (-) operations, we don't have to do these calculation over prime fields.
So, our target should be rho function (matrix?) calculation.

Reciprocal difference

> {-
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
