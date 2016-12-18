Multivariate.lhs

> module Multivariate where

> import Data.Ratio
> import Data.List (transpose)
> import Univariate ( list2pol, list2npolTimes
>                   , tDegree, list2rat'
>                   )

Let us start 2-variate polynomials.

  *Multivariate> let f z1 z2 = 3+2*z1+4*z2+7*z1^2+5*z1*z2+6*z2^2
  *Multivariate> [[f x y | y <- [0..9]] | x <- [0..9]]
  [[3,13,5,69,115,173,243,325,419,525]
  ,[12,27,54,93,144,207,282,369,468,579]
  ,[35,55,87,131,187,255,335,427,531,647]
  ,[72,97,134,183,244,317,402,499,608,729]
  ,[123,153,195,249,315,393,483,585,699,825]
  ,[188,223,270,329,400,483,578,685,804,935]
  ,[267,37,359,423,499,587,687,799,923,1059]
  ,[360,405,462,531,612,705,810,927,1056,1197]
  ,[467,517,579,653,739,837,947,1069,1203,1349]
  ,[588,643,710,789,880,983,1098,1225,1364,1515]
  ]

Assuming the list of lists is a matrix of 2-variate function's values,
  f i j

> tablize :: (Enum t1, Num t1) => (t1 -> t1 -> t) -> Int -> [[t]]
> tablize f n = [[f x y | y <- range] | x <- range]
>   where
>     range = take n [0..]
  
  *Multivariate> tablize (\x y -> (x,y)) 4
  [[(0,0),(0,1),(0,2),(0,3)]
  ,[(1,0),(1,1),(1,2),(1,3)]
  ,[(2,0),(2,1),(2,2),(2,3)]
  ,[(3,0),(3,1),(3,2),(3,3)]
  ]

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

Let us take the transpose of this "matrix" to see the behavior of coefficients.

  *Multivariate> let f z1 z2 = 3+2*z1+4*z2+7*z1^2+5*z1*z2+6*z2^2
  *Multivariate> let fTable = tablize f 10
  *Multivariate> map list2pol fTable 
  [[3 % 1,4 % 1,6 % 1]
  ,[2 % 1,9 % 1,6 % 1]
  ,[35 % 1,14 % 1,6 % 1]
  ,[72 % 1,19 % 1,6 % 1]
  ,[123 % 1,24 % 1,6 % 1]
  ,[188 % 1,29 % 1,6 % 1]
  ,[267 % 1,34 % 1,6 % 1]
  ,[360 % 1,39 % 1,6 % 1]
  ,[467 % 1,44 % 1,6 % 1]
  ,[588 % 1,49 % 1,6 % 1]
  ]
  *ultivariate> transpose it
  [[3 % 1,12 % 1,35 % 1,72 % 1,123 % 1,188 % 1,267 % 1,360 % 1,467 % 1,588 % 1]
  [4 % 1,9 % 1,14 % 1,19 % 1,24 % 1,29 % 1,34 % 1,39 % 1,44 % 1,49 % 1]
  ,[6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1,6 % 1]
  ]
  *Multivariate> map list2pol it
  [[3 % 1,2 % 1,7 % 1]
  ,[4 % 1,5 % 1]
  ,[6 % 1]]

> table2pol :: [[Ratio Integer]] -> [[Ratio Integer]]
> table2pol = map list2pol . transpose . map list2pol

  *Multivariate> let g x y = 1+7*x + 8*y + 10*x^2 + x*y+9*y^2
  *Multivariate> table2pol $ tablize g 5
  [[1 % 1,7 % 1,10 % 1],[8 % 1,1 % 1],[9 % 1]]

There are some bad-behavior polynomials;
  *Multivariate> table2pol $ tablize (\x y -> x*y) 20
  [[0 % 1],[1 % 1,1 % 1]]
  *Multivariate> tablize (\x y -> (x,y)) 5
  [[(0,0),(0,1),(0,2),(0,3),(0,4)]
  ,[(1,0),(1,1),(1,2),(1,3),(1,4)]
  ,[(2,0),(2,1),(2,2),(2,3),(2,4)]
  ,[(3,0),(3,1),(3,2),(3,3),(3,4)]
  ,[(4,0),(4,1),(4,2),(4,3),(4,4)]
  ]
  *Multivariate> tablize (\x y -> (x*y)) 5
  [[0,0,0,0,0]
  ,[0,1,2,3,4]
  ,[0,2,4,6,8]
  ,[0,3,6,9,12]
  ,[0,4,8,12,16]
  ]

Here we have assumed that the list of functions has the same length, but

  *Multivariate> map list2pol $ tablize (\x y -> x*y) 5
  [[0 % 1],[0 % 1,1 % 1],[0 % 1,2 % 1],[0 % 1,3 % 1],[0 % 1,4 % 1]]

So, we should repeat 0's if we have zero-function.

> xyDegree f = (dX, dY)
>   where
>     dX = length . list2pol $ map (\t -> f t 1) [0..] 
>     dY = length . list2pol $ map (\t -> f 1 t) [0..]

  *Multivariate> let test x y = x^2*(2*y + y^3)
  *Multivariate> uncurry (*) . xyDegree $ test
  6
  *Multivariate> maximum . map (length . list2pol) . tablize test $ 6
  4
  *Multivariate> map (take 4 . (++ (repeat (0%1))) . list2pol) . tablize test $ 6
  [[0 % 1,0 % 1,0 % 1,0 % 1]
  ,[0 % 1,2 % 1,0 % 1,1 % 1]
  ,[0 % 1,8 % 1,0 % 1,4 % 1]
  ,[0 % 1,18 % 1,0 % 1,9 % 1]
  ,[0 % 1,32 % 1,0 % 1,16 % 1]
  ,[0 % 1,50 % 1,0 % 1,25 % 1]
  ]
  *Multivariate> map list2pol . transpose $ it
  [[0 % 1],[0 % 1,0 % 1,2 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]

  *Multivariate> let test x y = (1%3)*x^2*((2%5)*y + ((3%4)*x*y^3)) 
                           -- = (2%15)*x^2*y + (1%4)*x^3*y^3
  *Multivariate> xyDegree  test
  (3,3)
  *Multivariate> map (take 4 . (++ (repeat (0%1))) . list2pol) . tablize test $ 9
  [[0 % 1,0 % 1,0 % 1,0 % 1]
  ,[0 % 1,2 % 15,0 % 1,1 % 4]
  ,[0 % 1,8 % 15,0 % 1,2 % 1]
  ,[0 % 1,6 % 5,0 % 1,27 % 4]
  ,[0 % 1,32 % 15,0 % 1,16 % 1]
  ,[0 % 1,10 % 3,0 % 1,125 % 4]
  ,[0 % 1,24 % 5,0 % 1,54 % 1]
  ,[0 % 1,98 % 15,0 % 1,343 % 4]
  ,[0 % 1,128 % 15,0 % 1,128 % 1]
  ]
  *Multivariate> map list2pol . transpose $ it
  [[0 % 1],[0 % 1,0 % 1,2 % 15],[0 % 1],[0 % 1,0 % 1,0 % 1,1 % 4]]

> xyPol2Coef :: (Enum t, Integral a, Num t) =>
>               (t -> t -> Ratio a) -> [[Ratio a]]
> xyPol2Coef f = map list2pol . transpose . map (take num . (++ (repeat (0%1))) . list2pol) . tablize f $ rank
>   where
>     rank = uncurry (*) . xyDegree $ f
>     num  = maximum . map (length . list2pol) . tablize f $ rank 

  *Multivariate> let test x y = (1%3)*x^2*((2%5)*y + ((3%4)*x*y^3))
  *Multivariate> xyPol2Coef test
  [[0 % 1],[0 % 1,0 % 1,2 % 15],[0 % 1],[0 % 1,0 % 1,0 % 1,1 % 4]]
  *Multivariate> let test2 x y = x*y
  *Multivariate> xyPol2Coef test2
  [[0 % 1],[0 % 1,1 % 1]]
  *Multivariate> let test3 x y = x^3*y^4
  *Multivariate> xyPol2Coef test3
  [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,0 % 1,1 % 1]]

> table2pol' :: Integral a => [[Ratio a]] -> [[Ratio a]]
> table2pol' tbl = map list2pol . transpose . map (take num . (++ (repeat (0%1))) . list2pol) $ tbl
>   where
>     num = maximum . map (length . list2pol) $ tbl

--
2-variable rational functions

  *Univariate> let h x y = (1+2*x+4*y+7*x^2+5*x*y+6*y^2) % (1+7*x+8*y+10*x^2+x*y+9*y^2)
  *Univariate> let auxh x y t = h (t*x) (t*y)
  *Univariate> let auxhs = [map (auxh 1 y) [0..] | y <- [0..]]
  *Univariate> take 10 $ map list2rat' auxhs 
  [Just ([1 % 1,2 % 1,7 % 1],[1 % 1,7 % 1,10 % 1])
  ,Just ([1 % 1,6 % 1,18 % 1],[1 % 1,15 % 1,20 % 1])
  ,Just ([1 % 1,10 % 1,41 % 1],[1 % 1,23 % 1,48 % 1])
  ,Just ([1 % 1,14 % 1,76 % 1],[1 % 1,31 % 1,94 % 1])
  ,Just ([1 % 1,18 % 1,123 % 1],[1 % 1,39 % 1,158 % 1])
  ,Just ([1 % 1,22 % 1,182 % 1],[1 % 1,47 % 1,240 % 1])
  ,Just ([1 % 1,26 % 1,253 % 1],[1 % 1,55 % 1,340 % 1])
  ,Just ([1 % 1,30 % 1,336 % 1],[1 % 1,63 % 1,458 % 1])
  ,Just ([1 % 1,34 % 1,431 % 1],[1 % 1,71 % 1,594 % 1])
  ,Just ([1 % 1,38 % 1,538 % 1],[1 % 1,79 % 1,748 % 1])
  ]
  *Univariate> sequence it
  Just [([1 % 1,2 % 1,7 % 1],[1 % 1,7 % 1,10 % 1])
       ,([1 % 1,6 % 1,18 % 1],[1 % 1,15 % 1,20 % 1])
       ,([1 % 1,10 % 1,41 % 1],[1 % 1,23 % 1,48 % 1])
       ,([1 % 1,14 % 1,76 % 1],[1 % 1,31 % 1,94 % 1])
       ,([1 % 1,18 % 1,123 % 1],[1 % 1,39 % 1,158 % 1])
       ,([1 % 1,22 % 1,182 % 1],[1 % 1,47 % 1,240 % 1])
       ,([1 % 1,26 % 1,253 % 1],[1 % 1,55 % 1,340 % 1])
       ,([1 % 1,30 % 1,336 % 1],[1 % 1,63 % 1,458 % 1])
       ,([1 % 1,34 % 1,431 % 1],[1 % 1,71 % 1,594 % 1])
       ,([1 % 1,38 % 1,538 % 1],[1 % 1,79 % 1,748 % 1])
       ]
  *Univariate> fmap (map fst) it
  Just [[1 % 1,2 % 1,7 % 1]
       ,[1 % 1,6 % 1,18 % 1]
       ,[1 % 1,10 % 1,41 % 1]
       ,[1 % 1,14 % 1,76 % 1]
       ,[1 % 1,18 % 1,123 % 1]
       , .. 
       ]
  *Univariate> fmap transpose it
  Just [[1 % 1,1 % 1,1 % 1,1 % 1,1 % 1,1 % 1,1 % 1,1 % 1,1 % 1,1 % 1]
       ,[2 % 1,6 % 1,10 % 1,14 % 1,18 % 1,22 % 1,26 % 1,30 % 1,34 % 1,38 % 1]
       ,[7 % 1,18 % 1,41 % 1,76 % 1,123 % 1,182 % 1,253 % 1,336 % 1,431 % 1,538 % 1]
       ]
  *Univariate> fmap (map list2pol) it
  Just [[1 % 1],[2 % 1,4 % 1],[7 % 1,5 % 1,6 % 1]]

  *Multivariate> let h x y = (1+2*x+4*y+7*x^2+5*x*y+(6%13)*y^2) / (1+(7%3)*x+8*y+10*x^2+x*y+9*y^2)
  *Multivariate> let auxh x y t = h (t*x) (t*y)
  *Multivariate> let auxhs = [map (auxh 1 y) [0..100] | y <- [0..100]]
  *Multivariate> fmap (map list2pol . transpose . map fst) . sequence . map list2rat' $ auxhs
  Just [[1 % 1],[2 % 1,4 % 1],[7 % 1,5 % 1,6 % 13]]
  *Multivariate> fmap (map list2pol . transpose . map snd) . sequence . map list2rat' $ auxhs
  Just [[1 % 1],[7 % 3,8 % 1],[10 % 1,1 % 1,9 % 1]]

> -- table2ratf :: Integral a => [[Ratio a]] -> ([[Ratio a]], [[Ratio a]])
> table2ratf :: Integral a =>
>               [[Ratio a]] -> (Maybe [[Ratio a]], Maybe [[Ratio a]])
> table2ratf table = (t2r fst table, t2r snd table)
>   where
> --    t2r' third = fmap (map third) . sequence . map list2rat'
>
>     myMax Nothing   = 0
>     myMax (Just ns) = ns
>
>     t2r third = fmap (map list2pol . transpose . map (take num . (++ (repeat (0%1))) . third)) . 
>                 sequence . map list2rat'   
>       where
>         num = myMax . fmap (maximum . map (length . fst)) . sequence . map list2rat' $ table
> 
> -- fmap (maximum . map (length . fst)) . sequence . map list2rat'
> -- map (take num . (++ (repeat (0%1))) . list2pol)
>
> tablizer :: (Num a, Enum a) => (a -> a -> b) -> a -> [[b]]
> tablizer f n = [map (f_t 1 y) [0..(n-1)] | y <- [0..(n-1)]]
>   where
>     f_t x y t = f (t*x) (t*y)

  *Multivariate> let h x y = (1+2*x+4*y+7*x^2+5*x*y+(6%13)*y^2) / (1+(7%3)*x+8*y+10*x^2+x*y+9*y^2)
  *Multivariate> let hTable = tablizer h 20
  *Multivariate> table2ratf hTable 
  (Just [[1 % 1],[2 % 1,4 % 1],[7 % 1,5 % 1,6 % 13]],Just [[1 % 1],[7 % 3,8 % 1],[10 % 1,1 % 1,9 % 1]])

Note that, the sampling points for n=10 case are
  
  *Multivariate> tablizer (\x y -> (x,y)) 10
  [[(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(6,0),(7,0),(8,0),(9,0)]
  ,[(0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(6,6),(7,7),(8,8),(9,9)]
  ,[(0,0),(1,2),(2,4),(3,6),(4,8),(5,10),(6,12),(7,14),(8,16),(9,18)]
  ,[(0,0),(1,3),(2,6),(3,9),(4,12),(5,15),(6,18),(7,21),(8,24),(9,27)]
  ,[(0,0),(1,4),(2,8),(3,12),(4,16),(5,20),(6,24),(7,28),(8,32),(9,36)]
  ,[(0,0),(1,5),(2,10),(3,15),(4,20),(5,25),(6,30),(7,35),(8,40),(9,45)]
  ,[(0,0),(1,6),(2,12),(3,18),(4,24),(5,30),(6,36),(7,42),(8,48),(9,54)]
  ,[(0,0),(1,7),(2,14),(3,21),(4,28),(5,35),(6,42),(7,49),(8,56),(9,63)]
  ,[(0,0),(1,8),(2,16),(3,24),(4,32),(5,40),(6,48),(7,56),(8,64),(9,72)]
  ,[(0,0),(1,9),(2,18),(3,27),(4,36),(5,45),(6,54),(7,63),(8,72),(9,81)]
  ]

  *Multivariate> let f x y = (1 + 2*x + 3*y + 4*x^2 + (1%5)*x*y + (1%6)*y^2) / (7 + 8*x + (1%9)*y + x^2 + x*y + 10*y^2)
  *Multivariate> let g x y = (11 + 10*x + 9*y) / (8 + 7*x^2 + (1%6)*x*y + 5*y^2)
  *Multivariate> table2ratf $ tablizer f 20
  (Just [[1 % 7],[2 % 7,3 % 7],[4 % 7,1 % 35,1 % 42]],Just [[1 % 1],[8 % 7,1 % 63],[1 % 7,1 % 7,10 % 7]])
  *Multivariate> table2ratf $ tablizer g 20
  (Just [[11 % 8],[5 % 4,9 % 8],[0 % 1]],Just [[1 % 1],[0 % 1],[7 % 8,1 % 48,5 % 8]])
  *Multivariate> let h x y = (f x y) / (g x y)
  *Multivariate> table2ratf $ tablizer h 20
  (Just [[8 % 77],[16 % 77,24 % 77],[39 % 77,53 % 2310,19 % 231],[2 % 11,64 % 231,3 % 22,15 % 77],[4 % 11,31 % 1155,106 % 385,37 % 2772,5 % 462]]
  ,Just [[1 % 1],[158 % 77,578 % 693],[13 % 11,757 % 693,111 % 77],[10 % 77,19 % 77,109 % 77,90 % 77]]
  )
  *Multivariate> table2ratf $ tablizer f 10
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer f 13
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer f 15
  (Just [[1 % 7],[2 % 7,3 % 7],[4 % 7,1 % 35,1 % 42]],Just [[1 % 1],[8 % 7,1 % 63],[1 % 7,1 % 7,10 % 7]])
  *Multivariate> table2ratf $ tablizer g 11
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer g 13
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer g 15
  (Just [[11 % 8],[5 % 4,9 % 8],[0 % 1]],Just [[1 % 1],[0 % 1],[7 % 8,1 % 48,5 % 8]])
  *Multivariate> table2ratf $ tablizer (\x y -> x/(1+x^2)) 10
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer (\x y -> x/(1+x^2)) 13
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer (\x y -> x/(1+x^2)) 20
  (Just [[0 % 1],[1 % 1],[0 % 1]],Just [[1 % 1],[0 % 1],[1 % 1]])
  *Multivariate> table2ratf $ tablizer (\x y -> x/(1+x^2)) 13
  (*** Exception: Prelude.!!: index too large
  *Multivariate> table2ratf $ tablizer (\x y -> x/(1+x^2)) 15
  (Just [[0 % 1],[1 % 1],[0 % 1]],Just [[1 % 1],[0 % 1],[1 % 1]])

> wilFunc x y = (x^2*y^2) % ((1 + y)^3)

  *Multivariate> table2ratf $ tablizer wilFunc 20
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1,3 % 1],[0 % 1,0 % 1,3 % 1],[0 % 1,0 % 1,0 % 1,1 % 1],[0 % 1]]) 
