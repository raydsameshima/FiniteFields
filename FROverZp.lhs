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
       Ratio a -> Ratio

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
  *FROverZp> let fps = map (f `fmodp` 101) [0..]
  *FROverZp> take 5 fps
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
>       | a==b && b==c = n
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
  
> newtonCp :: (Integral a, Integral t) => a -> [t] -> [t]
> newtonCp p xs = [x `div` factorial k | (x,k) <- zip xs [0..(p-1)]]
>   where
>     factorial k = product [1.. fromIntegral k]

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> difListsp 101 [take 10 $ accessibleData f 101]
  [[14,14,14,14,14,14,14,14]
  ,[48,62,76,90,3,17,31,45,59]
  ,[34,82,43,18,7,10,27,58,2,61]
  ]
  *FROverZp> reverse $ map head it
  [34,48,14]
  *FROverZp> newtonCp 101 it
  [34,48,7]

> firstDifsp :: Integral a => a -> [a] -> [a]
> firstDifsp p xs = reverse $ map head $ difListsp p [xs]

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> firstDifsp 101 $ accessibleData f 101
  [34,48,14]

> list2npolp :: Integral t => t -> [t] -> [t]
> list2npolp p xs = newtonCp p $ firstDifsp p $ take n xs
>   where n = (degreep p xs) + 2

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> let fsp p = list2npolp p $ accessibleData f p
  *FROverZp> fsp 101
  [34,48,7]
  *FROverZp> npol2pol $ fsp 101
  [34,41,7]
  *FROverZp> npol2pol $ fsp 103
  [69,83,7]
  *FROverZp> npol2pol $ fsp 107
  [36,22,7]

> list2polp' :: Integral t => t -> [t] -> [t]
> list2polp' p xs = npol2pol $ list2npolp p xs

  *FROverZp> list2polp' 101 $ map (f `fmodp` 101) [0..]
  [34,41,7]
  *FROverZp> list2polp' 103 $ map (f `fmodp` 103) [0..]
  [69,83,7]
  *FROverZp> list2polp' 107 $ map (f `fmodp` 107) [0..]
  [36,22,7]

We ready to guess these data:
























> guessPol' :: Integral t => t -> [t] -> [Ratio t]
> guessPol' p xs = map (fst . (\a -> guess (a, p))) $ list2polp' p xs  

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> guessPol' 101 $ map (f `fmodp` 101) [0..]
  [1 % 3,3 % 5,7 % 1]
  *FROverZp> guessPol' 103 $ map (f `fmodp` 103) [0..]
  [1 % 3,3 % 5,7 % 1]
  *FROverZp> guessPol' 107 $ map (f `fmodp` 107) [0..]
  [1 % 3,3 % 5,7 % 1]

  *FROverZp> map (\p -> guessPol' p  (map (`fmodp f p) [0..])) [101,103,107]
  [[1 % 3,3 % 5,7 % 1],[1 % 3,3 % 5,7 % 1],[1 % 3,3 % 5,7 % 1]]
  *FROverZp> matches3 $ map (\p -> guessPol' p  (map (f `fmodp` p) [0..])) primes
  [1 % 3,3 % 5,7 % 1]

> reconstructPol fs = 
>   matches3 $ map (\p -> guessPol' p (accessibleData' fs p)) bigPrimes

  *FROverZp> let f x = (1%3) + (3%5)*x + 7*x^2
  *FROverZp> let fs = map f [0..]
  *FROverZp> reconstructPol fs
  [1 % 3,3 % 5,7 % 1]
  







