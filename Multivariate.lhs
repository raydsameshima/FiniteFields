Multivariate.lhs

> module Multivariate 
>   where

> import Univariate ( degree, list2pol
>                   , thiele2ratf, lists2ratf, thiele2coef, lists2rat
>                   )

Multivariate Polynomial case
Example:
Consider a polynomial of the numeratorin (3.23).

  *Univariate> let f z1 z2 = 3+2*z1+4*z2+7*z1^2+5*z1*z2+6*z2^2
  *Univariate> :t f
  f :: Num a => a -> a -> a
  *Univariate> let f1s = map (`f` 0) [0..]
  *Univariate> degreeLazy f1s
  2
  *Univariate> degree $ take 20 f1s
  2
  *Univariate> list2pol $ take 20 f1s
  [3 % 1,2 % 1,7 % 1]
  *Univariate> let f2s = map (f 0) [0..]
  *Univariate> degreeLazy f2s
  2
  *Univariate> degree $ take 20 f2s 
  2
  *Univariate> list2pol $ take 20 f2s
  [3 % 1,4 % 1,6 % 1]
  *Univariate> let f1s' = map (`f` 1) [0..]
  *Univariate> list2pol $ take 20 f1s'
  [13 % 1,7 % 1,7 % 1]
  *Univariate> list2pol $ take 20 $ map (`f` 2) [0..]
  [35 % 1,12 % 1,7 % 1]
  *Univariate> degree [2,7,11]
  *** Exception: difLists: lack of data, or not a polynomial
  CallStack (from HasCallStack):
    error, called at Univariate.lhs:61:19 in main:Univariate
  *Univariate> degree [2,7,12]
  1
  *Univariate> list2pol [2,7,12]
  [2 % 1,5 % 1]
  *Univariate> degree [3, 13, 35]
  *** Exception: difLists: lack of data, or not a polynomial
  CallStack (from HasCallStack):
    error, called at Univariate.lhs:61:19 in main:Univariate
  *Univariate> list2pol $ take 20 $ map (`f` 3) [0..]
  [69 % 1,17 % 1,7 % 1]
  *Univariate> degree [3, 13, 35, 69]
  2
  *Univariate> list2pol [3,13,35,69]
  [3 % 1,4 % 1,6 % 1]


