Univariate.lhs

> module Univariate where

References
J. Stoer, R. Bulirsch
  Introduction to Numerical Analysis (2nd edition)
L. M. Milne-Thomson
  THE CALCULUS OF FINITE DIFFERENCES

> import Control.Applicative
> import Control.Monad
> import Data.Ratio
> import Data.Maybe
> import Data.List
> -- import Control.Monad.Catch
>
> import Polynomials 

From the output list
  map f [0..]
of a polynomial
  f :: Int -> Ratio Int
we reconstrunct the canonical form of f.

> difs :: (Num a) => [a] -> [a]
> difs [] = []
> difs [_] = []
> difs (i:jj@(j:js)) = j-i : difs jj
>
> difLists :: (Eq a, Num a) => [[a]] -> [[a]]
> difLists []          = []
> difLists xx@(xs:_) =
>   if isConst xs 
>     then xx
>     else difLists $ difs xs : xx
>   where
>     isConst (i:jj@(j:_)) = all (==i) jj
>     isConst _ = error "difLists: lack of data, or not a polynomial"
>
> -- This degree function is "strict", so only take finite list.
> degree' :: (Eq a, Num a) => [a] -> Int
> degree' xs = length (difLists [xs]) - 1
>
> -- This degree function can compute the degree of infinite list.
> degreeLazy :: (Eq a, Num a) => [a] -> Int
> degreeLazy xs = helper xs 0
>   where
>     helper as@(a:b:c:_) n
>       | a==b && b==c = n
>       | otherwise    = helper (difs as) (n+1)
>
> -- This is a hyblid version, safe and lazy.
> degree :: (Num a, Eq a) => [a] -> Int
> degree xs = let l = degreeLazy xs in
>   degree' $ take (l+2) xs

> -- m-times match version
> degreeTimes :: (Num a, Eq a) => Int -> [a] -> Int
> degreeTimes m xs = helper xs 0
>   where
>     helper aa@(a:as) n
>       | all (== a) (take (m-1) as) = n
>       | otherwise                  = helper (difs aa) (n+1)

Newton interpolation formula
First we introduce a new infix symbol for the operation 
of taking a falling power.

> infixr 8 ^- -- falling power
> (^-) :: (Eq a, Num a) => a -> a -> a
> x ^- 0 = 1
> x ^- n = (x ^- (n-1)) * (x - n + 1)

Claim (Newton interpolation formula):
  A polynomial f of degree n is expressed as
    f(z) = \sum_{k=0}^n  (diff^n(f)(0)/k!) * (x ^- k)
  where diff^n(f) is the n-th difference of f.

Example
Consider a polynomial f(x) = 2*x^3+3*x.

In general, we have no prior knowledge of this form, 
but we know the sequences as a list of outputs (map f [0..]):

  Univariate> let f x = 2*x^3+3*x
  Univariate> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  Univariate> degree $ take 10 $ map f [0..]
  3

Let us try to get differences:

  Univariate> difs $ take 10 $ map f [0..]
  [5,17,41,77,125,185,257,341,437]
  Univariate> difs it
  [12,24,36,48,60,72,84,96]
  Univariate> difs it
  [12,12,12,12,12,12,12]

Or more simply take difLists:

  Univariate> difLists [take 10 $ map f [0..]]
  [[12,12,12,12,12,12,12]
  ,[12,24,36,48,60,72,84,96]
  ,[5,17,41,77,125,185,257,341,437]
  ,[0,5,22,63,140,265,450,707,1048,1485]
  ]

What we need is the heads of above lists.

  Univariate> map head it
  [12,12,5,0]
  
Newton interpolation formula gives
  f' x = 0*(x ^- 0) `div` (0!) + 5*(x ^- 1) `div` (1!) 
         + 12*(x ^- 2) `div` (2!) + 12*(x ^- 3) `div` (3!)
       = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)
So

  Univariate> let f x = 2*x^3+3*x
  Univariate> let f' x = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)
  Univariate> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  Univariate> take 10 $ map f' [0..]
  [0,5,22,63,140,265,450,707,1048,1485]

Assume the differences are given in a list
  [x_0, x_1 ..]
where x_k = diff^k(f)(0).
Then the implementation of the Newton interpolation formula is as follows:

> newtonC 
>   :: (Fractional t, Enum t) => 
>      [t] -- first differences
>   -> [t] -- Newton coefficients
> newtonC xs = [x / factorial k | (x,k) <- zip xs [0..]]
>   where
>     factorial k = product [1..fromInteger k]

  Univariate> let f x = 2*x^3+3*x
  Univariate> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  Univariate> difLists [it]
  [[12,12,12,12,12,12,12]
  ,[12,24,36,48,60,72,84,96]
  ,[5,17,41,77,125,185,257,341,437]
  ,[0,5,22,63,140,265,450,707,1048,1485]
  ]
  Univariate> reverse $ map head it
  [0,5,12,12]
  Univariate> newtonC it
  [0 % 1,5 % 1,6 % 1,2 % 1]

The list of first differences can be computed as follows:

> firstDifs 
>   :: (Eq a, Num a) => 
>      [a] -- map f [0..]
>   -> [a]
> firstDifs xs = reverse . map head . difLists $ [xs]

Mapping a list of integers to a Newton representation:

> -- This implementation can take infinite list.
> list2npol :: (Integral a) => [Ratio a] -> [Ratio a]
> list2npol xs = newtonC . firstDifs $ take n xs
>   where n = (degree xs) + 2
>
> -- m-times matches version
> list2npolTimes :: (Integral a) => Int -> [Ratio a] -> [Ratio a]
> list2npolTimes m xs = newtonC . firstDifs $ take n xs
>   where n = (degreeTimes m xs) + 2

  *Univariate> let f x = x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)
  *Univariate> let fs = map f [0..]
  *Univariate> list2npolTimes 10 fs
  [0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,1 % 1]
  *Univariate> npol2pol it
  [0 % 1,(-120) % 1,274 % 1,(-225) % 1,85 % 1,(-15) % 1,1 % 1]
  *Univariate> list2npol fs
  [0 % 1]

We need to map Newton falling powers to standard powers.  
This is a matter of applying combinatorics, by means of a convention formula 
that uses the so-called Stirling cyclic numbers (of the first kind.)
Its defining relation is
  (x ^- n) = \sum_{k=1}^n (stirlingC n k) * (-1)^(n-k) * x^k.
The key equation is
  (x ^- n) = (x ^- (n-1)) * (x-n+1)
           = x*(x ^- (n-1)) - (n-1)*(x ^- (n-1))
  
Therefore, an implementation is as follows:

> stirlingC :: (Integral a) => a -> a -> a
> stirlingC 0 0 = 1
> stirlingC 0 _ = 0
> stirlingC n k = stirlingC (n-1) (k-1) + (n-1) * stirlingC (n-1) k  

This definition can be used to convert from falling powers to standard powers.

> fall2pol :: (Integral a) => a -> [a]
> fall2pol 0 = [1]
> fall2pol n = 0   -- No constant term. 
>            : [(-1)^(n-k) * stirlingC n k| k<-[1..n]]

We use this to convert Newton representations to standard polynomials 
in coefficients list representation.
Here we have uses 
  sum 
to collect same order terms in list representation.

> npol2pol :: (Ord t, Num t) => [t] -> [t]
> npol2pol xs = sum [ [x] * map fromInteger (fall2pol k)
>                   | (x,k) <- zip xs [0..]
>                   ]

Finally, here is the function for reconstruction the polynomial
from an output sequence:

> list2pol :: (Integral a) => [Ratio a] -> [Ratio a]
> list2pol = npol2pol . list2npol

Reconstruction as curve fitting

  *Univariate> let f x = 2*x^3 + 3*x + 1%5
  *Univariate> take 10 $ map f [0..]
  [1%5, 26%5, 111%5, 316%5, 701%5, 1326%5, 2251%5, 3536%5, 5241%5, 7426%5]
  *Univariate> list2npol it
  [1 % 5,5 % 1,6 % 1,2 % 1]
  *Univariate> list2npol $ map f [0..]
  [1 % 5,5 % 1,6 % 1,2 % 1]
  *Univariate> list2pol $ map (\n -> 1%3 + (3%5)*n + (5%7)*n^2) [0..]
  [1 % 3,3 % 5,5 % 7]
  *Univariate>  list2pol [0,1,5,14,30,55]
  [0 % 1,1 % 6,1 % 2,1 % 3]
  *Univariate> map (p2fct $ list2pol [0,1,5,14,30,55]) [0..6]
  [0 % 1,1 % 1,5 % 1,14 % 1,30 % 1,55 % 1,91 % 1]

Here is n-times match version:

> list2polTimes :: Integral a => Int -> [Ratio a] -> [Ratio a]
> list2polTimes n = npol2pol . list2npolTimes n

--

Thiele's interpolation formula
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#Haskell
http://mathworld.wolfram.com/ThielesInterpolationFormula.html

reciprocal difference
Using the same notation of 
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#C

> rho :: (Integral a) => 
>        [Ratio a] -- A list of output of f :: a -> Ratio a 
>     -> a -> Int  -- "matrix"
>     -> Maybe (Ratio a) -- Nothing means 1/0 type infinity 
> rho fs 0 i = Just $ fs !! i
> rho fs n i 
>   | n < 0         = Just 0
>   | num == Just 0 = Nothing -- "infinity"
>   | otherwise     = (+) <$> recipro <*> rho fs (n-2) (i+1)
>   where
>     recipro = ((%) . (* n) <$> den) <*> num -- (den*n)%num
> --            (%) <$> (*n) <$> den <*> num -- functor law
>     num  = numerator <$> next
>     den  = denominator <$> next
> --    next = (-) <$> rho fs (n-1) (i+1) <*> rho fs (n-1) i
>     next = x `seq` y `seq` (-) <$> x <*> y
>       where x = rho fs (n-1) (i+1) 
>             y = rho fs (n-1) i

Note that (%) has the following type,
  (%) :: Integral a => a -> a -> Ratio a

  *Univariate> (%) <$> (*2) <$> Just 5 <*> Just 3
  Just (10 % 3)

The follwoing reciprocal differences match the table of 
Milne-Thompson[1951] page 106:

  *Univariate> map (\p -> map (rho (map (\t -> 1%(1+t^2)) [0..]) p) [0..3]) [0..5]
  [[Just (1 % 1)   ,Just (1 % 2)    ,Just (1 % 5)    ,Just (1 % 10)]
  ,[Just ((-2) % 1),Just ((-10) % 3),Just ((-10) % 1),Just ((-170) % 7)]
  ,[Just ((-1) % 1),Just ((-1) % 10),Just ((-1) % 25),Just ((-1) % 46)]
  ,[Just (0 % 1)   ,Just (40 % 1)   ,Just (140 % 1)  ,Just (324 % 1)]
  ,[Just (0 % 1)   ,Just (0 % 1)    ,Just (0 % 1)    ,Just (0 % 1)]
  ,[Nothing,Nothing,Nothing,Nothing]
  ]

> -- Thiele coefficients (continuous fraction)
> a :: (Integral a) => [Ratio a] -> a -> Maybe (Ratio a)
> a fs 0 = Just $ head fs
> a fs n = (-) <$> rho fs n 0 <*> rho fs (n-2) 0
>
> -- shifted Thiele coefficients
> a' :: Integral a => [Ratio a] -> Int -> a -> Maybe (Ratio a)
> a' fs p 0 = Just $ fs !! p
> a' fs p n = (-) <$> rho fs n p <*> rho fs (n-2) p

  *Univariate> map (\p -> map (rho (map (\t -> t%(1+t^2)) [0..]) p) [0..5]) [0..5]
  [[Just (0%1),Just (1%2)    ,Just (2%5)    ,Just (3%10)     ,Just (4%17)     ,Just (5%26)]
  ,[Just (2%1),Just ((-10)%1),Just ((-10)%1),Just ((-170)%11),Just ((-442)%19),Just ((-962)%29)]
  ,[Just (1%3),Nothing       ,Just ((-1)%15),Just ((-1)%48)  ,Just ((-1)%105) ,Just ((-1)%192)]
  ,[Nothing   ,Nothing       ,Just (50%1)   ,Just (242%1)    ,Just (662%1)    ,Just (1430%1)]
  ,[Nothing   ,Nothing       ,Just (0%1)    ,Just (0%1)      ,Just (0%1)      ,Just (0%1)]
  ,[Nothing   ,Nothing       ,Nothing       ,Nothing         ,Nothing         ,Nothing]
  ]

Here, the consecutive Just ((-10) % 1) in second list make "fake" infinity (Nothing).

  *Univariate> let f t = t%(1+t^2)
  *Univariate> let fs = map f [0..]
  *Univariate> let aMat = [map (a' fs i) [0..] | i <- [0..]]
  *Univariate> take 20 $ map (length . takeWhile isJust) $ aMat 
  [3,2,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

> -- Thiele coefficients with shifts.
> aMatrix :: Integral a => [Ratio a] -> [[Maybe (Ratio a)]]
> aMatrix fs = [map (a' fs i) [0..] | i <- [0..]]
> 
> tDegree :: Integral a => [Ratio a] -> Int
> tDegree = isConsts' 3 . map (length . takeWhile isJust) . aMatrix
> 
> -- To find constant sub sequence.
> isConsts' :: Eq t => Int -> [t] -> t
> isConsts' n (l:ls)  
>   | all (==l) $ take (n-1) ls = l 
>   | otherwise                 = isConsts' n ls

we also need the shift, in this case, p=2 to get full Thiele coefficients.

> shiftaMatrix 
>   :: Integral a => 
>      [Ratio a] -> [Maybe [Ratio a]]
> shiftaMatrix gs = map (sequence . (\q -> map (a' gs q) [0..(thieleD-1)])) [0..]
>   where
>     thieleD = fromIntegral $ tDegree gs
>
> shiftAndThieleC 
>   :: Integral a => 
>      [Ratio a] -> (Maybe Int, Maybe [Ratio a])
> shiftAndThieleC fs = (findIndex isJust gs, join $ find isJust gs)
>   where
>     gs = shiftaMatrix fs

  *Univariate> take 10 $ map sequence $ transpose $ take (tDegree fs) m
  [Just [0 % 1,1 % 2,2 % 5,3 % 10,4 % 17]
  ,Just [2 % 1,(-10) % 1,(-10) % 1,(-170) % 11,(-442) % 19]
  ,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing]

Packed version, this scans the given data only once.

> degSftTC
>   :: Integral a =>
>      [Ratio a] -> (Int, Maybe Int, Maybe (Maybe [Ratio a]))
> --                 |    |          + Thiele coefficients
> --                 |    + shift
> --                 + degree
> degSftTC fs = (d,s,ts)
>   where
>     m = [map (a' fs i) [0..] | i <- [0..]]
>     d = isConsts' 3 . map (length . takeWhile isJust) $ m -- 3 times match
>     m' = map (sequence . take d) m
>     s = findIndex isJust m'
>     ts = find isJust m'

  *Univariate Control.Monad> let g t = t%(1+t^2)
  *Univariate Control.Monad> let gs = map g [0..]
  *Univariate Control.Monad> shiftAndThieleC $ shiftaMatrix gs
  (Just 2,Just [2 % 5,(-10) % 1,(-7) % 15,60 % 1,1 % 15])
  *Univariate Control.Monad> let f t = 1%(1+t^2)
  *Univariate Control.Monad> let fs = map f [0..]
  *Univariate Control.Monad> shiftAndThieleC $ shiftaMatrix fs
  (Just 0,Just [1 % 1,(-2) % 1,(-2) % 1,2 % 1,1 % 1])

We need a convertor from this thiele sequence to continuous fractional form of rational function.

> nextStep [a0,a1] (v:_)  = a0 + v/a1
> nextStep (a:as)  (v:vs) = a + (v / nextStep as vs)
>
> -- From thiele sequence to (rational) function.
> thiele2ratf :: Integral a => [Ratio a] -> (Ratio a -> Ratio a)
> thiele2ratf as x
>   | x == 0    = head as -- only constant term
>   | otherwise = nextStep as [x,x-1 ..]

  *Univariate> let h t = (3+6*t+18*t^2)%(1+2*t+20*t^2)
  *Univariate> let hs = map h [0..]
  *Univariate> let as = thieleC hs
  *Univariate> as
  [3 % 1,(-23) % 42,(-28) % 13,767 % 14,7 % 130]
  *Univariate> let th x = thiele2ratf as x
  *Univariate> take 5 hs
  [3 % 1,27 % 23,87 % 85,183 % 187,45 % 47]
  *Univariate> map th [0..5]
  [3 % 1,27 % 23,87 % 85,183 % 187,45 % 47,69 % 73]

We represent a rational function by a tuple of coefficient lists:
  (ns,ds) :: ([Ratio Int],[Ratio Int])
where ns and ds are coef-list-rep of numerator polynomial and denominator polynomial. 
Here is a translator from coefficients lists to rational function.

> -- similar to p2fct
> lists2ratf :: (Integral a) => 
>               ([Ratio a],[Ratio a]) -> (Ratio a -> Ratio a)
> lists2ratf (ns,ds) x = p2fct ns x / p2fct ds x

  *Univariate> let frac x = lists2ratf ([1,1%2,1%3],[2,2%3]) x
  *Univariate> take 10 $ map frac [0..]
  [1 % 2,11 % 16,1 % 1,11 % 8,25 % 14,71 % 32,8 % 3,25 % 8,79 % 22,65 % 16]
  *Univariate> let ffrac x = (1+(1%2)*x+(1%3)*x^2)/(2+(2%3)*x)
  *Univariate> take 10 $ map ffrac [0..]
  [1 % 2,11 % 16,1 % 1,11 % 8,25 % 14,71 % 32,8 % 3,25 % 8,79 % 22,65 % 16]

The following canonicalizer reduces the tuple-rep of rational function in canonical form
That is, the coefficien of the lowest degree term of the denominator to be 1.
However, since our input starts from 0 and this means firstNonzero is the same as head.

> canonicalize :: (Integral a) => ([Ratio a],[Ratio a]) -> ([Ratio a],[Ratio a])
> canonicalize rat@(ns,ds)
>   | dMin == 1 = rat
>   | otherwise = (map (/ dMin) ns, map (/ dMin) ds)
>   where
>     dMin = firstNonzero ds
>     firstNonzero [a] = a -- head
>     firstNonzero (a:as)
>       | a /= 0    = a
>       | otherwise = firstNonzero as

What we need is a translator from Thiele coefficients to this tuple-rep.

> thiele2coef :: (Integral a) => [Ratio a] -> ([Ratio a],[Ratio a])
> thiele2coef as = canonicalize $ t2r as 0
>   where
>     t2r [an,an'] n = ([an*an'-n,1],[an'])
>     t2r (a:as)   n = ((a .* num) + ([-n,1] * den), num)
>       where
>         (num, den) = t2r as (n+1)
>

  *Univariate> let h t = (3+6*t+18*t^2)%(1+2*t+20*t^2)
  *Univariate> let hs = map h [0..]
  *Univariate> take 5 hs
  [3 % 1,27 % 23,87 % 85,183 % 187,45 % 47]
  *Univariate> let th x = thiele2ratf as x
  *Univariate> map th [0..5]
  [3 % 1,27 % 23,87 % 85,183 % 187,45 % 47,69 % 73]
  *Univariate> as
  [3 % 1,(-23) % 42,(-28) % 13,767 % 14,7 % 130]
  *Univariate> thiele2coef as
  ([3 % 1,6 % 1,18 % 1],[1 % 1,2 % 1,20 % 1])

> thiele2coef' -- shifted version (0 -> sft)
>   :: Integral a => 
>      Ratio a -> [Ratio a] -> ([Ratio a], [Ratio a])
> thiele2coef' sft [a] = ([a],1)
> thiele2coef' sft as = canonicalize $ t2r as sft
>   where
>     t2r [an,an'] n = (([an*an'-n] + z),[an'])
>     t2r (a:as)   n = ((a .* num) + ((z - [n]) * den), num)
>       where
>         (num, den) = t2r as (n+1)

  *Univariate> let f x = x^2%(1+x^2)
  *Univariate> shiftAndThieleC $ map f [0..]
  (Just 0,Just [0 % 1,2 % 1,2 % 1,(-2) % 1,(-1) % 1])
  *Univariate> take 3 $ shiftaMatrix (map f [0..])
  [Just [0 % 1,2 % 1,2 % 1,(-2) % 1,(-1) % 1]
  ,Just [1 % 2,10 % 3,3 % 5,(-130) % 3,(-1) % 10]
  ,Just [4 % 5,10 % 1,6 % 25,(-150) % 1,(-1) % 25]
  ]
  *Univariate> thiele2coef' 0 [0 % 1,2 % 1,2 % 1,(-2) % 1,(-1) % 1]
  ([0 % 1,0 % 1,1 % 1],[1 % 1,0 % 1,1 % 1])
  *Univariate> thiele2coef' 1 [1 % 2,10 % 3,3 % 5,(-130) % 3,(-1) % 10]
  ([0 % 1,0 % 1,1 % 1],[1 % 1,0 % 1,1 % 1])

> shiftAndThiele2coef (Just sft, Just ts) = Just $ thiele2coef' (fromIntegral sft) ts
> shiftAndThiele2coef _                   = Nothing
>
> list2rat' :: (Integral a) => [Ratio a] -> Maybe ([Ratio a], [Ratio a])
> list2rat' = shiftAndThiele2coef . shiftAndThieleC
>
> list2rat'' lst = let (_,s,ts) = degSftTC lst in
>   shiftAndThiele2coef (s, join ts)  
 
  *Univariate> let f t = t%(1+t^2)
  *Univariate> let fs = map (\t -> 1%(1+t^2)) [0..]
  *Univariate> list2rat' fs
  Just ([1 % 1,0 % 1,0 % 1],[1 % 1,0 % 1,1 % 1])
  *Univariate> shiftAndThieleC fs
  (Just 0,Just [1 % 1,(-2) % 1,(-2) % 1,2 % 1,1 % 1])
  *Univariate> let fs = map (\t -> t%(1+t^2)) [0..]
  *Univariate> let gs = map (\t -> t%(1+t^2)) [0..]
  *Univariate> shiftAndThieleC gs
  (Just 2,Just [2 % 5,(-10) % 1,(-7) % 15,60 % 1,1 % 15])
  *Univariate> list2rat' gs
  Just ([0 % 1,1 % 1,0 % 1],[1 % 1,0 % 1,1 % 1])
  *Univariate> let hs = map (\t -> t^2%(1+t^2)) [0..]
  *Univariate> shiftAndThieleC hs
  (Just 0,Just [0 % 1,2 % 1,2 % 1,(-2) % 1,(-1) % 1])
  *Univariate> list2rat' hs
  Just ([0 % 1,0 % 1,1 % 1],[1 % 1,0 % 1,1 % 1])

  *Univariate> let f t = t%(1+t^2)
  *Univariate> let fs = map f [0..]
  *Univariate> let aMat = [map (\j -> a' fs j i) [0..] | i <- [0..]]
  *Univariate> take 6 $ map (take 5) aMat
  [[Just (0 % 1),Just (1 % 2),Just (2 % 5),Just (3 % 10),Just (4 % 17)]
  ,[Just (2 % 1),Just ((-10) % 1),Just ((-10) % 1),Just ((-170) % 11),Just ((-442) % 19)]
  ,[Just (1 % 3),Nothing,Just ((-7) % 15),Just ((-77) % 240),Just ((-437) % 1785)]
  ,[Nothing,Nothing,Just (60 % 1),Just (2832 % 11),Just (13020 % 19)]
  ,[Nothing,Nothing,Just (1 % 15),Just (1 % 48),Just (1 % 105)]
  ,[Nothing,Nothing,Nothing,Nothing,Nothing]
  ] 

--

For general, non-sequential input.

to do:
Need try-catch-throw type exception-handling. 
Control.Monad.Catch

> graph :: (a -> b) -> [a] -> [(a,b)] 
> graph _ []     = []
> graph f (a:as) = (a, f a): graph f as

--

For non-sequential input-output.
Assume
  [(x, f x) | x<-domain]
as the given data.

> type Q = Ratio Int
>
> divDif :: ([Q], Q) -> ([Q], Q) -> ([Q], Q)
> divDif ([x0], f0) ([x1], f1) = ([x0,x1], (f1 - f0) / (x1 - x0))
> divDif (xs  , f0) (xs' , f1) = ((z:xs'), (f1-f0) / range)
>   where
>     a = last xs'
>     z = head xs
>     range = a - z
>
> map' :: (a -> a -> a) -> [a] -> [a]
> map' _ []            = []
> map' _ [a]           = []
> map' f (a:bb@(b:bs)) = (f a b) : map' f bb
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
  [[([0 % 1],0 % 1),([1 % 1],1 % 1),([2 % 1],4 % 1),([3 % 1],9 % 1),([4 % 1],16 % 1),([5 % 1],25 % 1)]
  ,[([0 % 1,1 % 1],1 % 1),([1 % 1,2 % 1],3 % 1),([2 % 1,3 % 1],5 % 1),([3 % 1,4 % 1],7 % 1),([4 % 1,5 % 1],9 % 1)]
  ,[([0 % 1,1 % 1,2 % 1],1 % 1),([1 % 1,2 % 1,3 % 1],1 % 1),([2 % 1,3 % 1,4 % 1],1 % 1),([3 % 1,4 % 1,5 % 1],1 % 1)]
  ,[([0 % 1,1 % 1,2 % 1,3 % 1],0 % 1),([1 % 1,2 % 1,3 % 1,4 % 1],0 % 1),([2 % 1,3 % 1,4 % 1,5 % 1],0 % 1)]]

