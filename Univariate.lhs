Univariate.lhs

> module Univariate where
> import Data.Ratio
>
> -- polynomials, as coefficients lists
> instance (Num a, Ord a) => Num [a] where
>   fromInteger c = [fromInteger c] 
>   -- operator overloading
>   negate []     = []
>   negate (f:fs) = negate f : negate fs
> 
>   signum [] = []
>   signum gs 
>     | signum (last gs) < 0 = negate z
>     | otherwise = z
> 
>   abs [] = []
>   abs gs 
>     | signum gs == z = gs
>     | otherwise      = negate gs
> 
>   fs     + []     = fs
>   []     + gs     = gs
>   (f:fs) + (g:gs) = f+g : fs+gs
> 
>   fs     * []     = []
>   []     * gs     = []
>   (f:fs) * gg@(g:gs) = f*g : (f .* gs + fs * gg)
>
> --
>
> -- scalar multiplication
> infixl 7 .*
> (.*) :: Num a => a -> [a] -> [a]
> c .* []     = []
> c .* (f:fs) = c*f : c .* fs
>
> -- z of f(z), variable
> z :: Num a => [a]
> z = [0,1]
>
> -- difference analysis
> difs :: (Integral n) => [n] -> [n]
> difs [] = []
> difs [_] = []
> difs (i:jj@(j:js)) = j-i : difs jj
>
> difLists :: (Integral n) => [[n]] -> [[n]]
> difLists [] = []
> difLists xx@(xs:xss) =
>   if isConst xs then xx
>                 else difLists $ difs xs : xx
>   where
>     isConst (i:jj@(j:js)) = all (==i) jj
>     isConst _ = error "difLists: lack of data, or not a polynomial"
>
> degree :: (Integral n) => [n] -> Int
> degree xs = length (difLists [xs]) -1
>
> p2fct :: Num a => [a] -> a -> a
> p2fct [] x = 0
> p2fct (a:as) x = a + (x * p2fct as x)
>

Newton interpolation formula
First we introduce a new infix symbol for the operation of taking a falling power.

> infixr 8 ^- -- falling power
> (^-) :: (Integral a) => a -> a -> a
> x ^- 0 = 1
> x ^- n = (x ^- (n-1)) * (x - n + 1)

Claim (Newton interpolation formula)
A polynomial f of degree n is expressed as
  f(z) = \sum_{k=0}^n  (diff^n(f)(0)/k!) * (x ^- n)
where diff^n(f) is the n-th difference of f.

Example
Consider

> f x = 2*x^3+3*x

In general, we have no prior knowledge of the form of this polynomial, but we know the sequences as a list of outputs:

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

> f' x = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)

  > take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  > take 10 $ map f' [0..]
  [0,5,22,63,140,265,450,707,1048,1485]

Assume the differences are given in a list
  [x_0, x_1 ..]
where x_i = diff^k(f)(0).
Then the implementation of the Newton interpolation formula is as follows:

> newton :: Integral a => [a] -> [Ratio a]
> newton xs = [x % factorial k | (x,k) <- zip xs [0..]]
>   where
>     factorial k = product [1..fromInteger k]

  > take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  > difLists [it]
  [[12,12,12,12,12,12,12],[12,24,36,48,60,72,84,96],[5,17,41,77,125,185,257,341,437],[0,5,22,63,140,265,450,707,1048,1485]]
  > reverse $ map head it
  [0,5,12,12]
  > newton it
  [0 % 1,5 % 1,6 % 1,2 % 1]

The list of first differences can be computed as follows:

> firstDifs :: [Integer] -> [Integer]
> firstDifs xs = reverse $ map head $ difLists [xs]

Mapping a list of integers to a Newton representation:

> list2npol :: [Integer] -> [Rational]
> list2npol = newton . map fromInteger . firstDifs

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

> stirlingC :: Integer -> Integer -> Integer
> stirlingC 0 0 = 1
> stirlingC 0 _ = 0
> stirlingC n k = (n-1)*(stirlingC (n-1) k) + stirlingC (n-1) (k-1)

This definition can be used to convert from falling powers to standard powers.

> fall2pol :: Integer -> [Integer]
> fall2pol 0 = [1]
> fall2pol n = 0 : [(stirlingC n k)*(-1)^(n-k) | k<-[1..n]]

We use this to convert Newton representations to standard polynomials in coefficients list representation.
Here we have uses sum to collect same order terms in list representation.

> npol2pol :: (Ord t, Num t) => [t] -> [t]
> npol2pol xs = sum [ [x] * (map fromInteger $ fall2pol k)
>                   | (x,k) <- zip xs [0..]
>                   ]

Finally, here is the function for computing a polynomial from an output sequence:

> list2pol :: [Integer] -> [Rational]
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

Thiele's interpolation formula
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#Haskell
http://mathworld.wolfram.com/ThielesInterpolationFormula.html

reciprocal difference

> rho :: [Ratio Int] -> Int -> Int -> Ratio Int
> rho fs 0 i = fs !! i
> rho fs n i = (n*den)%num + rho fs (n-2) (i+1)
>   where
>     num = numerator next
>     den = denominator next
>     next = (rho fs (n-1) (i+1)) - (rho fs (n-2) (i+1))

Note that (%) has the following type,
  (%) :: Integral a => a -> a -> Ratio a


