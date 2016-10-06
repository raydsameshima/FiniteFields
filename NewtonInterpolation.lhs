NewtonInterpolation.lhs

http://homepages.cwi.nl/~jve/HR/PolAddendum.pdf

> module NewtonInterpolation where
> import POL
> import Polynomials
> import Data.Ratio

First we introduce a new infix symbol for the operation of taking a falling power.

> infixr 8 ^-
> (^-) :: (Integral a) => a -> a -> a
> x ^- 0 = 1
> x ^- n = (x ^- (n-1)) * (x - n + 1)

In a similar way, we can define rising powers.

> infixr 8 ^+
> (^+) :: (Integral a) => a -> a -> a
> x ^+ 0 = 1
> x ^+ n = (x ^+ (n-1)) * (x + n - 1)

Claim (Newton interpolation formula)
A polynomial f of degree n is expressed as
  f(z) = \sum_{k=0}^n  (diff^n(f)(0)/k!) * (x ^- n)
where diff^n(f) is the n-th difference of f.

Example
Consider

> f x = 2*x^3+3*x

In general, we have no prior knowledge of the form of this polynomial, but we know the sequences as a list of outputs:

  *NewtonInterpolation> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  *NewtonInterpolation> degree $ take 10 $ map f [0..]
  3

Let us try to get differences:

  *NewtonInterpolation> difs $ take 10 $ map f [0..]
  [5,17,41,77,125,185,257,341,437]
  *NewtonInterpolation> difs it
  [12,24,36,48,60,72,84,96]
  *NewtonInterpolation> difs it
  [12,12,12,12,12,12,12]

Or more simply take difLists:

  *NewtonInterpolation> difLists [take 10 $ map f [0..]]
  [[12,12,12,12,12,12,12]
  ,[12,24,36,48,60,72,84,96]
  ,[5,17,41,77,125,185,257,341,437]
  ,[0,5,22,63,140,265,450,707,1048,1485]
  ]

What we need is the heads of above lists.

  *NewtonInterpolation> map head it
  [12,12,5,0]
  
Newton interpolation formula gives
  f' x = 0*(x ^- 0) `div` (0!) + 5*(x ^- 1) `div` (1!) + 12*(x ^- 2) `div` (2!) + 12*(x ^- 3) `div` (3!)
       = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)
So

> f' x = 5*(x ^- 1) + 6*(x ^- 2) + 2*(x ^- 3)

  *NewtonInterpolation> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  *NewtonInterpolation> take 10 $ map f' [0..]
  [0,5,22,63,140,265,450,707,1048,1485]

Assume the differences are given in a list
  [x_0, x_1 ..]
where x_i = diff^k(f)(0).
Then the implementation of the Newton interpolation formula is as follows:

> newton :: Integral a => [a] -> [Ratio a]
> newton xs = [x % factorial k | (x,k) <- zip xs [0..]]
>   where
>     factorial k = product [1..fromInteger k]

  *NewtonInterpolation> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  *NewtonInterpolation> difLists [it]
  [[12,12,12,12,12,12,12],[12,24,36,48,60,72,84,96],[5,17,41,77,125,185,257,341,437],[0,5,22,63,140,265,450,707,1048,1485]]
  *NewtonInterpolation> reverse $ map head it
  [0,5,12,12]
  *NewtonInterpolation> newton it
  [0 % 1,5 % 1,6 % 1,2 % 1]

The list of first differences can be computed as follows:

> firstDifs :: [Integer] -> [Integer]
> firstDifs xs = reverse $ map head $ difLists [xs]

Mapping a list of integers to a Newton representation:

> list2npol :: [Integer] -> [Rational]
> list2npol = newton . map fromInteger . firstDifs

  *NewtonInterpolation> take 10 $ map f [0..]
  [0,5,22,63,140,265,450,707,1048,1485]
  *NewtonInterpolation> list2npol it
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
  *NewtonInterpolation> list2pol $ map (\n -> 7*n^2+3*n-4) [0..100]
  [(-4) % 1,3 % 1,7 % 1]

  *NewtonInterpolation> list2pol [0,1,5,14,30]
  [0 % 1,1 % 6,1 % 2,1 % 3]
  *NewtonInterpolation> map (\n -> n%6 + n^2%2 + n^3%3) [0..4]
  [0 % 1,1 % 1,5 % 1,14 % 1,30 % 1]

  *NewtonInterpolation> map (p2fct $ list2pol [0,1,5,14,30]) [0..8]
  [0 % 1,1 % 1,5 % 1,14 % 1,30 % 1,55 % 1,91 % 1,140 % 1,204 % 1]
