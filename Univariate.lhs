Univariate.lhs

References
J. Stoer, R. Bulirsch
  Introduction to Numerical Analysis (2nd edition)
L. M. Milne-Thomson
  THE CALCULUS OF FINITE DIFFERENCES

> module Univariate where
> import Control.Applicative
> import Control.Monad
> import Data.Ratio
> import Data.Maybe
> import Data.List
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
>
> degreeTimes :: (Num a, Eq a) => Int -> [a] -> Int
> degreeTimes m xs = helper xs 0
>   where
>     helper aa@(a:as) n
>       | all (== a) (take m as) = n
>       | otherwise              = helper (difs aa) (n+1)

Newton interpolation formula
First we introduce a new infix symbol for the operation of taking a falling power.

> infixr 8 ^- -- falling power
> (^-) :: (Eq a, Num a) => a -> a -> a
> x ^- 0 = 1
> x ^- n = (x ^- (n-1)) * (x - n + 1)

Claim (Newton interpolation formula)
A polynomial f of degree n is expressed as
  f(z) = \sum_{k=0}^n  (diff^n(f)(0)/k!) * (x ^- n)
where diff^n(f) is the n-th difference of f.

Example
Consider a polynomial f = 2*x^3+3*x.

In general, we have no prior knowledge of this form, but we know the sequences as a list of outputs:

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
  f' x = 0*(x ^- 0) `div` (0!) + 5*(x ^- 1) `div` (1!) + 12*(x ^- 2) `div` (2!) + 12*(x ^- 3) `div` (3!)
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
where x_i = diff^k(f)(0).
Then the implementation of the Newton interpolation formula is as follows:

> newtonC :: (Fractional t, Enum t) => [t] -> [t]
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

> firstDifs :: (Eq a, Num a) => [a] -> [a]
> firstDifs xs = reverse $ map head $ difLists [xs]

Mapping a list of integers to a Newton representation:

> -- This implementation can take infinite list.
> list2npol :: (Integral a) => [Ratio a] -> [Ratio a]
> list2npol xs = newtonC . firstDifs $ take n xs
>   where n = (degree xs) + 2
>
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
This is a matter of applying combinatorics, by means of a convention formula that uses the so-called Stirling cyclic numbers (of the first kind.)
Its defining relation is
  (x ^- n) = \sum_{k=1}^n (stirlingC n k) * (-1)^(n-k) * x^k.
The key equation is
  (x ^- n) = (x ^- (n-1)) * (x-n+1)
           = x*(x ^- (n-1)) - (n-1)*(x ^- (n-1))
  
Therefore, an implementation is as follows:

> stirlingC :: (Integral a) => a -> a -> a
> stirlingC 0 0 = 1
> stirlingC 0 _ = 0
> stirlingC n k = stirlingC (n-1) (k-1) + (n-1)*stirlingC (n-1) k  

This definition can be used to convert from falling powers to standard powers.

> fall2pol :: (Integral a) => a -> [a]
> fall2pol 0 = [1]
> fall2pol n = 0   -- No constant term. 
>            : [(-1)^(n-k) * stirlingC n k| k<-[1..n]]

We use this to convert Newton representations to standard polynomials in coefficients list representation.
Here we have uses sum to collect same order terms in list representation.

> -- For later convenience, we relax the type annotation.
> -- npol2pol :: (Integral a) => [Ratio a] -> [Ratio a]
> npol2pol :: (Ord t, Num t) => [t] -> [t]
> npol2pol xs = sum [ [x] * map fromInteger (fall2pol k)
>                   | (x,k) <- zip xs [0..]
>                   ]

Finally, here is the function for computing a polynomial from an output sequence:

> list2pol :: (Integral a) => [Ratio a] -> [Ratio a]
> list2pol = npol2pol . list2npol

Reconstruction as curve fitting
  *Univariate> let f x = 2*x^3 + 3*x + 1%5
  *Univariate> take 10 $ map f [0..]
  [1 % 5,26 % 5,111 % 5,316 % 5,701 % 5,1326 % 5,2251 % 5,3536 % 5,5241 % 5,7426 % 5]
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

--

Thiele's interpolation formula
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#Haskell
http://mathworld.wolfram.com/ThielesInterpolationFormula.html

reciprocal difference
Using the same notation of 
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#C

> rho :: (Integral a) => 
>        [Ratio a] -- A list of output of f :: a -> Ratio a 
>     -> a -> Int -> Maybe (Ratio a)
> rho fs 0 i = Just $ fs !! i
> rho fs n i 
>   | n < 0          = Just 0
>   | num == Just 0  = Nothing
>   | otherwise      = (+) <$> recipro <*> rho fs (n-2) (i+1)
>   where
>     recipro = (%) <$> (*n) <$> den <*> num
>     num  = numerator <$> next
>     den  = denominator <$> next
>     next = (-) <$> rho fs (n-1) (i+1) <*> rho fs (n-1) i

Note that (%) has the following type,
  (%) :: Integral a => a -> a -> Ratio a

  *Univariate> (%) <$> (*2) <$> Just 5 <*> Just 3
  Just (10 % 3)

This reciprocal difference rho matches the table of Milne-Thompson[1951] page 106:

  *Univariate> map (\p -> map (rho (map (\t -> 1%(1+t^2)) [0..]) p) [0..3]) [0..5]
  [[Just (1 % 1),Just (1 % 2),Just (1 % 5),Just (1 % 10)]
  ,[Just ((-2) % 1),Just ((-10) % 3),Just ((-10) % 1),Just ((-170) % 7)]
  ,[Just ((-1) % 1),Just ((-1) % 10),Just ((-1) % 25),Just ((-1) % 46)]
  ,[Just (0 % 1),Just (40 % 1),Just (140 % 1),Just (324 % 1)]
  ,[Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)]
  ,[Nothing,Nothing,Nothing,Nothing]
  ]

> a :: (Integral a) => [Ratio a] -> a -> Maybe (Ratio a)
> a fs 0 = Just $ head fs
> a fs n = (-) <$> rho fs n 0 <*> rho fs (n-2) 0
>
> -- shifted Thiele coefficients
> a' :: Integral a => [Ratio a] -> Int -> a -> Maybe (Ratio a)
> a' fs p 0 = Just $ fs !! p
> a' fs p n = (-) <$> rho fs n p <*> rho fs (n-2) p

  *Univariate> map (\p -> map (rho (map (\t -> t%(1+t^2)) [0..]) p) [0..5]) [0..5]
  [[Just (0 % 1),Just (1 % 2),Just (2 % 5),Just (3 % 10),Just (4 % 17),Just (5 % 26)]
  ,[Just (2 % 1),Just ((-10) % 1),Just ((-10) % 1),Just ((-170) % 11),Just ((-442) % 19),Just ((-962) % 29)]
  ,[Just (1 % 3),Nothing,Just ((-1) % 15),Just ((-1) % 48),Just ((-1) % 105),Just ((-1) % 192)]
  ,[Nothing,Nothing,Just (50 % 1),Just (242 % 1),Just (662 % 1),Just (1430 % 1)]
  ,[Nothing,Nothing,Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)]
  ,[Nothing,Nothing,Nothing,Nothing,Nothing,Nothing]
  ]

> tDegree :: Integral a => [Ratio a] -> a
> tDegree fs = helper fs 0
>   where
>     helper fs n
>       | areNothings fs' = n-1 -- We have taken one more, so need (-1)!
>       | otherwise       = helper fs (n+1)
>       where
>         fs' = map (rho fs n) [0..]
>     areNothings js = all (== Nothing) $ take 10 js -- 10 for practical reason

  *Univariate> map (\p -> map (rho (map (\t -> t^2%(1+t^2)) [0..]) p) [0..5]) [0..5]
  [[Just (0 % 1),Just (1 % 2),Just (4 % 5),Just (9 % 10),Just (16 % 17),Just (25 % 26)]
  ,[Just (2 % 1),Just (10 % 3),Just (10 % 1),Just (170 % 7),Just (442 % 9),Just (962 % 11)]
  ,[Just (2 % 1),Just (11 % 10),Just (26 % 25),Just (47 % 46),Just (74 % 73),Just (107 % 106)]
  ,[Just (0 % 1),Just ((-40) % 1),Just ((-140) % 1),Just ((-324) % 1),Just ((-616) % 1),Just ((-1040) % 1)]
  ,[Just (1 % 1),Just (1 % 1),Just (1 % 1),Just (1 % 1),Just (1 % 1),Just (1 % 1)]
  ,[Nothing,Nothing,Nothing,Nothing,Nothing,Nothing]
  ]

Using primed-a (a'), we can simply shift and reconstruct functions,
  
  *Univariate> let f = \t -> t^2%(1+t^2)
  *Univariate> let fs = map f [0..]
  *Univariate> tDegree fs
  4
  *Univariate> map (a' fs 0) [0..4]
  [Just (0 % 1),Just (2 % 1),Just (2 % 1),Just ((-2) % 1),Just ((-1) % 1)]
  *Univariate> map (a' fs 1) [0..4]
  [Just (1 % 2),Just (10 % 3),Just (3 % 5),Just ((-130) % 3),Just ((-1) % 10)]

Using Maxima, we get same results

  (%i17) f(t) := 0 + (x/(2 + (x-1)/(2 + (x-2)/(-2 + (x-3)/(-1))))), ratsimp;
  (%o17) f(t):=x^2/(1+x^2)
  (%i18) ff(t) := 1/2 + (x-1)/(10/3 + (x-2)/(3/5 + (x-3)/(-130/3 + (x-4)/(-1/10)))), ratsimp;
  (%o18) ff(t):=x^2/(1+x^2)

  *Univariate> let g t = t%(1+t^2)
  *Univariate> let gs = map g [0..]
  *Univariate> tDegree gs
  4
  *Univariate> map (a' gs 0) [0..4]
  [Just (0 % 1),Just (2 % 1),Just (1 % 3),Nothing,Nothing]
  *Univariate> map (a' gs 1) [0..4]
  [Just (1 % 2),Just ((-10) % 1),Nothing,Nothing,Nothing]
  *Univariate> map (a' gs 2) [0..4]
  [Just (2 % 5),Just ((-10) % 1),Just ((-7) % 15),Just (60 % 1),Just (1 % 15)]

  (%i19) g(x) := 2/5 + (x-2)/(-10 + (x-3)/(-7/15 + (x-4)/(60 + (x-5)/(1/15)))), ratsimp;
  (%o19) g(x):=x/(1+x^2)

  *Univariate> map (\q -> map (a' gs q) [0..4]) [0..5]
  [[Just (0 % 1),Just (2 % 1),Just (1 % 3),Nothing,Nothing]
  ,[Just (1 % 2),Just ((-10) % 1),Nothing,Nothing,Nothing]
  ,[Just (2 % 5),Just ((-10) % 1),Just ((-7) % 15),Just (60 % 1),Just (1 % 15)]
  ,[Just (3 % 10),Just ((-170) % 11),Just ((-77) % 240),Just (2832 % 11),Just (1 % 48)]
  ,[Just (4 % 17),Just ((-442) % 19),Just ((-437) % 1785),Just (13020 % 19),Just (1 % 105)]
  ,[Just (5 % 26),Just ((-962) % 29),Just ((-493) % 2496),Just (42432 % 29),Just (1 % 192)]
  ]
  *Univariate> map sequence it
  [Nothing
  ,Nothing
  ,Just [2 % 5,(-10) % 1,(-7) % 15,60 % 1,1 % 15]
  ,Just [3 % 10,(-170) % 11,(-77) % 240,2832 % 11,1 % 48]
  ,Just [4 % 17,(-442) % 19,(-437) % 1785,13020 % 19,1 % 105]
  ,Just [5 % 26,(-962) % 29,(-493) % 2496,42432 % 29,1 % 192]
  ]

We also need the shift, in this case, 2 to get full Thiele coefficients.

  *Univariate> findIndex isJust it
  Just 2

> rhoMatrix :: Integral a => [Ratio a] -> [Maybe [Ratio a]]
> rhoMatrix gs = map (sequence . (\q -> map (a' gs q) [0..thieleD])) [0..]
>   where
>     thieleD = tDegree gs
>
> shiftAndThieleC :: Integral a => 
>                    [Ratio a] -> (Maybe Int, Maybe [Ratio a])
> shiftAndThieleC fs = (findIndex isJust gs, join $ find isJust gs)
>   where
>     gs = rhoMatrix fs

  *Univariate Control.Monad> let g t = t%(1+t^2)
  *Univariate Control.Monad> let gs = map g [0..]
  *Univariate Control.Monad> shiftAndThieleC $ rhoMatrix gs
  (Just 2,Just [2 % 5,(-10) % 1,(-7) % 15,60 % 1,1 % 15])
  *Univariate Control.Monad> let f t = 1%(1+t^2)
  *Univariate Control.Monad> let fs = map f [0..]
  *Univariate Control.Monad> shiftAndThieleC $ rhoMatrix fs
  (Just 0,Just [1 % 1,(-2) % 1,(-2) % 1,2 % 1,1 % 1])

We need a convertor from this thiele sequence to continuous form of rational function.

> nextStep [a0,a1] (v:_)  = a0 + v/a1
> nextStep (a:as)  (v:vs) = a + (v / nextStep as vs)
>
> -- From thiele sequence to (rational) function.
> thiele2ratf :: Integral a => [Ratio a] -> (Ratio a -> Ratio a)
> thiele2ratf as x
>   | x == 0    = head as
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

> thiele2coef' :: Integral a => 
>                 Ratio a -> [Ratio a] -> ([Ratio a], [Ratio a])
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
  *Univariate> take 3 $ rho
  rho        rhoMatrix
  *Univariate> take 3 $ rhoMatrix (map f [0..])
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

> list2rat' :: (Integral a) => [Ratio a] -> Maybe ([Ratio a], [Ratio a])
> list2rat' = shiftAndThiele2coef . shiftAndThieleC
 
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
