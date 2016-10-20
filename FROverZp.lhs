> module FROverZp where
> import Data.Ratio
> import Data.Numbers.Primes
>
> import Ffield (modp, guess)
> import Univariate

Univariate Polynomial case
Our target is a univariate polynomial
  f :: (Integral a) =>
       Ratio a -> Ratio

