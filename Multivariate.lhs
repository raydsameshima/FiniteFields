Multivariate.lhs

> module Multivariate 
>   where

> import Data.Ratio
> import Univariate 
>   ( degree, list2pol
>   , thiele2ratf, lists2ratf, thiele2coef, list2rat
>   )

Let us start 2-variate polynomials.

  *Multivariate> let f z1 z2 = 3+2*z1+4*z2+7*z1^2+5*z1*z2+6*z2^2
  *Multivariate> [[f x y | y <- [0..9]] | x <- [0..9]]
  [[3,13,35,69,115,173,243,325,419,525]
  ,[12,27,54,93,144,207,282,369,468,579]
  ,[35,55,87,131,187,255,335,427,531,647]
  ,[72,97,134,183,244,317,402,499,608,729]
  ,[123,153,195,249,315,393,483,585,699,825]
  ,[188,223,270,329,400,483,578,685,804,935]
  ,[267,307,359,423,499,587,687,799,923,1059]
  ,[360,405,462,531,612,705,810,927,1056,1197]
  ,[467,517,579,653,739,837,947,1069,1203,1349]
  ,[588,643,710,789,880,983,1098,1225,1364,1515]
  ]

Assuming the list of lists is a matrix of 2-variate function's values, (f i j).

> tablize :: (Enum t1, Num t1) => (t1 -> t1 -> t) -> Int -> [[t]]
> tablize f n = [[f x y | y <- range] | x <- range]
>   where
>     range = take n [0..]

  *Multivariate> let fTable = tablize f 10
  *Multivariate> map list2pol fTable 
  [[3 % 1,4 % 1,6 % 1]
  ,[12 % 1,9 % 1,6 % 1]
  ,[35 % 1,14 % 1,6 % 1]
  ,[72 % 1,19 % 1,6 % 1]
  ,[123 % 1,24 % 1,6 % 1]
  ,[188 % 1,29 % 1,6 % 1]
  ,[267 % 1,34 % 1,6 % 1]
  ,[360 % 1,39 % 1,6 % 1]
  ,[467 % 1,44 % 1,6 % 1]
  ,[588 % 1,49 % 1,6 % 1]
  ]

> wellOrd :: [[a]] -> [[a]]
> wellOrd xss 
>   | null (head xss) = [] 
>   | otherwise       = map head xss : wellOrd (map tail xss)

  *Multivariate> let f z1 z2 = 3+2*z1+4*z2+7*z1^2+5*z1*z2+6*z2^2
  *Multivariate> let fTable = tablize f 10
  *Multivariate> map list2pol fTable 
  [[3 % 1,4 % 1,6 % 1]
  ,[12 % 1,9 % 1,6 % 1]
  ,[35 % 1,14 % 1,6 % 1]
  ,[72 % 1,19 % 1,6 % 1]
  ,[123 % 1,24 % 1,6 % 1]
  ,[188 % 1,29 % 1,6 % 1]
  ,[267 % 1,34 % 1,6 % 1]
  ,[360 % 1,39 % 1,6 % 1]
  ,[467 % 1,44 % 1,6 % 1]
  ,[588 % 1,49 % 1,6 % 1]
  ]
  *Multivariate> wellOrd it
  [[3 % 1,12 % 1,35 % 1,72 % 1,123 % 1,188 % 1,267 % 1,360 % 1,467 % 1,588 % 1]
  ,[4 % 1,9 % 1,14 % 1,19 % 1,24 % 1,29 % 1,34 % 1,39 % 1,44 % 1,49 % 1]
  ,[6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1]
  ]
  *Multivariate> map list2pol it
  [[3 % 1,2 % 1,7 % 1]
  ,[4 % 1,5 % 1]
  ,[6 % 1]]

> table2pol :: [[Ratio Integer]] -> [[Ratio Integer]]
> table2pol = map list2pol . wellOrd . map list2pol

  *Multivariate> let g x y = 1+7*x + 8*y + 10*x^2 + x*y+9*y^2
  *Multivariate> table2pol $ tablize g 5
  [[1 % 1,7 % 1,10 % 1],[8 % 1,1 % 1],[9 % 1]]

--

Next, 2-variate rational functions.
  
  *Multivariate> let h x y = (3+2*x+4*y+7*x^2+5*x*y+6*y^2) % (1+7*x+8*y+10*x^2+x*y+9*y^2)
  *Multivariate> let auxh x y t = h (t*x) (t*y)
  *Multivariate> let h x y = (3+2*x+4*y+7*x^2+5*x*y+6*y^2)% (1+7*x+8*y+10*x^2+x*y+9*y^2)
  *Multivariate> let auxh x y t = h (t*x) (t*y)

Using the homogenious property, we just take x=1:

  *Multivariate> let auxhs = [map (auxh 1 y) [0..5] | y <- [0..5]]
  *Multivariate> auxhs
  [[3 % 1,2 % 3,7 % 11,9 % 14,41 % 63,94 % 143]
  ,[3 % 1,3 % 4,29 % 37,183 % 226,105 % 127,161 % 192]
  ,[3 % 1,3 % 4,187 % 239,201 % 251,233 % 287,77 % 94]
  ,[3 % 1,31 % 42,335 % 439,729 % 940,425 % 543,1973 % 2506]
  ,[3 % 1,8 % 11,59 % 79,291 % 385,681 % 895,528 % 691]
  ,[3 % 1,23 % 32,155 % 211,1707 % 2302,1001 % 1343,4663 % 6236]
  ]
  
  *Multivariate> map list2rat auxhs
  [([3 % 1,2 % 1,7 % 1],[1 % 1,7 % 1,10 % 1])
  ,([3 % 1,6 % 1,18 % 1],[1 % 1,15 % 1,20 % 1])
  ,([3 % 1,10 % 1,41 % 1],[1 % 1,23 % 1,48 % 1])
  ,([3 % 1,14 % 1,76 % 1],[1 % 1,31 % 1,94 % 1])
  ,([3 % 1,18 % 1,123 % 1],[1 % 1,39 % 1,158 % 1])
  ,([3 % 1,22 % 1,182 % 1],[1 % 1,47 % 1,240 % 1])
  ]
  *Multivariate> map fst it
  [[3 % 1,2 % 1,7 % 1]
  ,[3 % 1,6 % 1,18 % 1]
  ,[3 % 1,10 % 1,41 % 1]
  ,[3 % 1,14 % 1,76 % 1]
  ,[3 % 1,18 % 1,123 % 1]
  ,[3 % 1,22 % 1,182 % 1]
  ]
  *Multivariate> wellOrd it
  [[3 % 1,3 % 1,3 % 1,3 % 1,3 % 1,3 % 1]
  ,[2 % 1,6 % 1,10 % 1,14 % 1,18 % 1,22 % 1]
  ,[7 % 1,18 % 1,41 % 1,76 % 1,123 % 1,182 % 1]
  ]
  *Multivariate> map list2pol it
  [[3 % 1],[2 % 1,4 % 1],[7 % 1,5 % 1,6 % 1]]

So, the numerator is given by

  *Multivariate> map list2pol . wellOrd . map (fst . list2rat) $ auxhs
  [[3 % 1],[2 % 1,4 % 1],[7 % 1,5 % 1,6 % 1]]
  
and the denominator is

  *Multivariate> map list2pol . wellOrd . map (snd . list2rat) $ auxhs
  [[1 % 1],[7 % 1,8 % 1],[10 % 1,1 % 1,9 % 1]]

> table2ratf table = (t2r fst table, t2r snd table)
>   where
>     t2r third = map list2pol . wellOrd . map (third . list2rat)
  
  *Multivariate> table2ratf auxhs
  ([[3 % 1],[2 % 1,4 % 1],[7 % 1,5 % 1,6 % 1]],[[1 % 1],[7 % 1,8 % 1],[10 % 1,1 % 1,9 % 1]])
