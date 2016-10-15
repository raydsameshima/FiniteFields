Univariate.lhs

> module Univariate where
> import Data.Ratio
> import Polynomials 

From the output list
  map f [0..]
of a polynomial
  f :: Int -> Ratio Int
we reconstrunct the canonical form of f.

> -- difference analysis
> difs :: (Num a) => [a] -> [a]
> difs [] = []
> difs [_] = []
> difs (i:jj@(j:js)) = j-i : difs jj
>
> difLists :: (Eq a, Num a) => [[a]] -> [[a]]
> difLists [] = []
> difLists xx@(xs:xss) =
>   if isConst xs then xx
>                 else difLists $ difs xs : xx
>   where
>     isConst (i:jj@(j:js)) = all (==i) jj
>     isConst _ = error "difLists: lack of data, or not a polynomial"
>
> -- This degree function is "strict", so only take finite list.
> degree' :: (Eq a, Num a) => [a] -> Int
> degree' xs = length (difLists [xs]) -1
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

  > let f x = 2*x^3+3*x
  > take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  > degree $ take 10 $ map f [0..]
  3

Let us try to get differences:

  > difs $ take 10 $ map f [0..]
  [5,17,41,77,125,185,257,341,437]
  > difs it
  [12,24,36,48,60,72,84,96]
  > difs it
  [12,12,12,12,12,12,12]

Or more simply take difLists:

  > difLists [take 10 $ map f [0..]]
  [[12,12,12,12,12,12,12]
  ,[12,24,36,48,60,72,84,96]
  ,[5,17,41,77,125,185,257,341,437]
  ,[0,5,22,63,140,265,450,707,1048,1485]
  ]

What we need is the heads of above lists.

  > map head it
  [12,12,5,0]
  
Newton interpolation formula gives
  f' x = 0*(x ^- 0) `div` (0!) + 5*(x ^- 1) `div` (1!) + 12*(x ^- 2) `div` (2!) + 12*(x ^- 3) `div` (3!)
       = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)
So

  > let f x = 2*x^3+3*x
  > let f' x = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)
  > take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  > take 10 $ map f' [0..]
  [0,5,22,63,140,265,450,707,1048,1485]

Assume the differences are given in a list
  [x_0, x_1 ..]
where x_i = diff^k(f)(0).
Then the implementation of the Newton interpolation formula is as follows:

> newtonC :: (Fractional t, Enum t) => [t] -> [t]
> newtonC xs = [x / factorial k | (x,k) <- zip xs [0..]]
>   where
>     factorial k = product [1..fromInteger k]

  > let f x = 2*x^3+3*x
  > take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  > difLists [it]
  [[12,12,12,12,12,12,12],[12,24,36,48,60,72,84,96],[5,17,41,77,125,185,257,341,437],[0,5,22,63,140,265,450,707,1048,1485]]
  > reverse $ map head it
  [0,5,12,12]
  > newtonC it
  [0 % 1,5 % 1,6 % 1,2 % 1]

The list of first differences can be computed as follows:

> firstDifs :: (Eq a, Num a) => [a] -> [a]
> firstDifs xs = reverse $ map head $ difLists [xs]

Mapping a list of integers to a Newton representation:

> list2npol :: (Integral a) => [Ratio a] -> [Ratio a]
> list2npol = newtonC . firstDifs

  > take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  > list2npol it
  [0 % 1,5 % 1,6 % 1,2 % 1]

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

> npol2pol :: (Integral a) => [Ratio a] -> [Ratio a]
> npol2pol xs = sum [ [x] * map fromInteger (fall2pol k)
>                   | (x,k) <- zip xs [0..]
>                   ]

Finally, here is the function for computing a polynomial from an output sequence:

> list2pol :: (Integral a) => [Ratio a] -> [Ratio a]
> list2pol = npol2pol . list2npol

Reconstruction as curve fitting
  > list2pol $ map (\n -> 7*n^2+3*n-4) [0..100]
  [(-4) % 1,3 % 1,7 % 1]

  > list2pol [0,1,5,14,30]
  [0 % 1,1 % 6,1 % 2,1 % 3]
  > map (\n -> n%6 + n^2%2 + n^3%3) [0..4]
  [0 % 1,1 % 1,5 % 1,14 % 1,30 % 1]

  > map (p2fct $ list2pol [0,1,5,14,30]) [0..8]
  [0 % 1,1 % 1,5 % 1,14 % 1,30 % 1,55 % 1,91 % 1,140 % 1,204 % 1]

-----------------------------------------------------------------------

Thiele's interpolation formula
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#Haskell
http://mathworld.wolfram.com/ThielesInterpolationFormula.html

reciprocal difference
Using the same notation of 
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#C

> rho :: (Integral a) => 
>        [Ratio a] -- A list of output of f :: a -> Ratio a 
>     -> a -> Int -> Ratio a
> rho fs 0 i = fs !! i
> rho fs n _ 
>   | n < 0 = 0
> rho fs n i = (n*den)%num + rho fs (n-2) (i+1)
>   where
>     num  = numerator next
>     den  = denominator next
>     next = rho fs (n-1) (i+1) - rho fs (n-1) i

Note that (%) has the following type,
  (%) :: Integral a => a -> a -> Ratio a

> a :: (Integral a) => [Ratio a] -> a -> Ratio a
> a fs 0 = head fs
> a fs n = rho fs n 0 - rho fs (n-2) 0

Consider the following continuous fraction form.
  (%i25) f(x) := 1+(x/(2+(x-1)/(3+(x-2)/4)));
  (%o25) f(x):=x/(2+(x-1)/(3+(x-2)/4))+1
  (%i26) ratsimp(f(x));
  (%o26) (x^2+16*x+16)/(16+6*x)

  *Univariate> map (a fs) [0..]
  [1 % 1,2 % 1,3 % 1,4 % 1,*** Exception: Ratio has zero denominator

  *Univariate> let func x = (x^2 + 16*x + 16)%(6*x + 16)
  *Univariate> let fs = map func [0..]
  *Univariate> take 5 $ map (rho fs 0) [0..]
  [1 % 1,3 % 2,13 % 7,73 % 34,12 % 5]
  *Univariate> take 5 $ map (rho fs 1) [0..]
  [2 % 1,14 % 5,238 % 69,170 % 43,230 % 53]
  *Univariate> take 5 $ map (rho fs 2) [0..]
  [4 % 1,79 % 16,269 % 44,667 % 88,413 % 44]
  *Univariate> take 5 $ map (rho fs 3) [0..]
  [6 % 1,6 % 1,6 % 1,6 % 1,6 % 1]

> tDegree :: Integral a => [Ratio a] -> a
> tDegree fs = helper fs 0
>   where
>     helper fs n
>       | isConstants fs' = n
>       | otherwise       = helper fs (n+1)
>       where
>         fs' = map (rho fs n) [0..]
>     isConstants (i:j:_) = i==j -- 2 times match
> --  isConstants (i:j:k_) = i==j && j==k -- 3 times match

  *Univariate> let h t = (3+6*t+18*t^2)%(1+2*t+20*t^2)
  *Univariate> let hs = map h [0..]
  *Univariate> tDegree hs
  4
  *Univariate> map (a hs) [0..(tDegree hs)]
  [3 % 1,(-23) % 42,(-28) % 13,767 % 14,7 % 130]

With Maxima,
  (%i35) h(t) := 3+t/((-23/42)+(t-1)/((-28/13)+(t-2)/((767/14)+(t-3)/(7/130))));

  (%o35) h(t):=t/((-23)/42+(t-1)/((-28)/13+(t-2)/(767/14+(t-3)/(7/130))))+3
  (%i36) ratsimp(h(t));

  (%o36) (18*t^2+6*t+3)/(1+2*t+20*t^2)

> thieleC :: (Integral a) => [Ratio a] -> [Ratio a]
> thieleC lst = map (a lst) [0..(tDegree lst)]

  *Univariate> thieleC hs
  [3 % 1,(-23) % 42,(-28) % 13,767 % 14,7 % 130]

We need a convertor from this thiele sequence to continuous form of rational function.

> nextStep [a0,a1] (v:_)  = a0 + v/a1
> nextStep (a:as)  (v:vs) = a + (v / nextStep as vs)
>
> -- From thiele sequence to (rational) function.
> thiele2ratf :: Integral a => [Ratio a] -> (Ratio a -> Ratio a)
> thiele2ratf as x
>   | x == 0 = head as
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
>   | otherwise = (map (/dMin) ns, map (/dMin) ds)
>   where
>     dMin = firstNonzero ds
>     firstNonzero [a] = a -- head
>     firstNonzero (a:as)
>       | a /= 0 = a
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
> lists2rat :: (Integral a) => [Ratio a] -> ([Ratio a], [Ratio a])
> lists2rat = thiele2coef . thieleC

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
