> module FROverZp where

Functional Reconstruction over finite field Z_p

> import Data.Ratio
> import Data.Numbers.Primes
>
> import Ffield (modp, guess, matches3, bigPrimes, reconstruct)
> import Univariate ((^-), stirlingC, fall2pol, npol2pol)

Univariate Polynomial case
Our target is a univariate polynomial
  f :: (Integral a) =>
       Ratio a -> Ratio a -- Real?

Let us consider

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> let fs = map f [0..]

So fs is our accessible data.
First, we should map (`modp` p) over this list, and take p elements from fs.

  *FROverZp> let fsp p = map (`modp` p) $ take p fs
  *FROverZp> take 10 $ fsp 101
  [34,82,43,18,7,10,27,58,2,61]
  *FROverZp> map (`mod` 101) $ difs it
  [48,62,76,90,3,17,31,45,59]
  *FROverZp> map (`mod` 101) $ difs it
  [14,14,14,14,14,14,14,14]
  
So, on Z_101, f is 2nd degree polynomial and is 
  34*(x ^- 0) + 48*(x ^- 1) + 14/(2!) * (x ^- 2)
   == 34 + 48*x + 7*(x ^-2)
   == 34 + 48*x + 7*x*(x-1)
   == 34 + 41*x + 7*x^2 (mod 101)

> -- Function-modular.
> fmodp :: Integral c => (a -> Ratio c) -> c -> a -> c
> f `fmodp` p = (`modp` p) . f

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> let fp = f `fmodp` 101
  *FROverZp> :t fp 
  fp :: Integral c => Ratio c -> c
  *FROverZp> take 10 $ map (f `fmodp` 101) [0..]
  [34,82,43,18,7,10,27,58,2,61]

Difference analysis over Z_p

> accessibleData :: (Ratio Int -> Ratio Int) -> Int -> [Int]
> accessibleData f p = take p $ map (f `fmodp` p) [0..]
> 
> accessibleData' :: [Ratio Int] -> Int -> [Int]
> accessibleData' fs p = take p $ map (`modp` p) fs
>
> difsp :: Integral b => b -> [b] -> [b]
> difsp p xs = map (`mod` p) (zipWith (-) (tail xs) xs)

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> let fs = map f [0..]
  *FROverZp> accessibleData' fs 101 == accessibleData f 101
  True
  *FROverZp> take 5 $ accessibleData f 101
  [34,82,43,18,7]
  *FROverZp> difsp 101 it
  [48,62,76,90]
  *FROverZp> difsp 101 it
  [14,14,14] 

> difListsp :: Integral b => b -> [[b]] -> [[b]]
> difListsp _ [] = []
> difListsp p xx@(xs:xxs) =
>   if isConst xs then xx
>                 else difListsp p $ difsp p xs : xx
>   where
>     isConst (i:jj@(j:js)) = all (==i) jj
>     isConst _ = error "difListsp: "

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> map head $ difListsp 101 [(accessibleData f 101)]
  [14,48,34]

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

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> let myDeg p = degreep p $ accessibleData f p
  *FROverZp> myDeg 101
  2
  *FROverZp> myDeg 103
  2
  *FROverZp> myDeg 107
  2

> firstDifsp :: Integral a => a -> [a] -> [a]
> firstDifsp p xs = reverse $ map head $ difListsp p [xs]

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> firstDifsp 101 $ accessibleData f 101
  [34,48,14]

> newtonCp :: (Integral a, Integral t) => a -> [t] -> [t]
> newtonCp p xs = [x `div` factorial k | (x,k) <- zip xs [0..(p-1)]]
>   where
>     factorial k = product [1.. fromIntegral k]




We guess these differences (at 0) then transform it as canonical form.






We ready to guess these data:

  *FROverZp> map (\p -> zip (npol2pol . fsp $ p) (repeat p)) [101,103,107]
  [[(34,101),(41,101),(7,101)],[(69,103),(83,103),(7,103)],[(36,107),(22,107),(7,107)]]
  *FROverZp> :t reconstruct
  reconstruct :: Integral a => [(a, a)] -> Ratio a
  *FROverZp> map head it
  [(34,101),(69,103),(36,107)]
  *FROverZp> reconstruct it
  1 % 3
  *FROverZp> map (\p -> zip (npol2pol . fsp $ p) (repeat p)) [101,103,107]
  [[(34,101),(41,101),(7,101)],[(69,103),(83,103),(7,103)],[(36,107),(22,107),(7,107)]]

  *FROverZp> map (head . tail) it
  [(41,101),(83,103),(22,107)]
  *FROverZp> reconstruct it
  3 % 5
  *FROverZp> map (\p -> zip (npol2pol . fsp $ p) (repeat p)) [101,103,107]
  [[(34,101),(41,101),(7,101)],[(69,103),(83,103),(7,103)],[(36,107),(22,107),(7,107)]]
  *FROverZp> map last it
  [(7,101),(7,103),(7,107)]
  *FROverZp> reconstruct it
  7 % 1

> wellOrd :: Eq a => [[a]] -> [[a]]
> wellOrd xs 
>   | head xs == [] = [] 
>   | otherwise     = map head xs : wellOrd (map tail xs)

  *FROverZp> wellOrd [[(34,101),(41,101),(7,101)],[(69,103),(83,103),(7,103)],[(36,107),(22,107),(7,107)]]
  [[(34,101),(69,103),(36,107)],[(41,101),(83,103),(22,107)],[(7,101),(7,103),(7,107)]]
  *FROverZp> map reconstruct it
  [1 % 3,3 % 5,7 % 1]

Here is the step-by-step usage of above functions:

  *FROverZp> let g x = (1%3) + (3%5)*x + (5%7)*x^2 + (7%9)*x^3
  *FROverZp> let gs = map g [0..]
  *FROverZp> let gsp p = list2npolp p $ accessibleData g p
  *FROverZp> let gData = map (\p -> zip (npol2pol . gsp $ p) (repeat p)) bigPrimes 
  *FROverZp> let gData' = wellOrd gData 
  *FROverZp> map reconstruct gData'
  [1 % 3,1 % 10,19 % 7,7 % 9]


