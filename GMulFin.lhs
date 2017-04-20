GMulFin.lhs

> module GMulFin where

Assume we can access
  f :: Q -> Q -> Q
of two-variable function.

> import Data.Ratio
> import Control.Monad (join)

> import GUniFin (Q, uniPolCoeff, uniRatCoeff, ratFunc2Coeff)
> import Multivariate (transposeWith)

> -- a test function
> wilFunc :: Q -> Q -> Q
> wilFunc x y = (x^2*y^2)/(1+y)^3

To track in-out correspondence, we should generalize the concept of graph:

> homogeneous 
>   :: (Q -> Q -> Q) -- 2var rational function
>   -> Q 
>   -> Q 
>   -> (Q -> Q)      -- 1var rational function
> homogeneous f x y t = f (t*x) (t*y)

  *GMulFin> :t homogeneous wilFunc 1 2
  homogeneous wilFunc 1 2 :: Q -> Q
  *GMulFin> ratFunc2Coeff (homogeneous wilFunc 1 2) 
  (Just [0 % 1,0 % 1,0 % 1,0 % 1,4 % 1],Just [1 % 1,6 % 1,12 % 1,8 % 1])
  *GMulFin> ratFunc2Coeff (homogeneous wilFunc 1 0) 
  (Just [0 % 1],Just [1 % 1])
  *GMulFin> ratFunc2Coeff (homogeneous wilFunc 1 1) 
  (Just [0 % 1,0 % 1,0 % 1,0 % 1,1 % 1],Just [1 % 1,3 % 1,3 % 1,1 % 1])
  *GMulFin> ratFunc2Coeff (homogeneous wilFunc 1 2) 
  (Just [0 % 1,0 % 1,0 % 1,0 % 1,4 % 1],Just [1 % 1,6 % 1,12 % 1,8 % 1])
  *GMulFin> ratFunc2Coeff (homogeneous wilFunc 1 3) 
  (Just [0 % 1,0 % 1,0 % 1,0 % 1,9 % 1],Just [1 % 1,9 % 1,27 % 1,27 % 1])
  *GMulFin> ratFunc2Coeff (homogeneous wilFunc 1 4) 
  (Just [0 % 1,0 % 1,0 % 1,0 % 1,16 % 1],Just [1 % 1,12 % 1,48 % 1,64 % 1])

We introduce homogeneous-function, and apply univariate rational function reconstruction.

  *GMulFin> map (\y -> ratFunc2Coeff (homogeneous wilFunc 1 y)) [0,1,3,5,6,8,9,11,13]
  [(Just [0 % 1],Just [1 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,1 % 1],Just [1 % 1,3 % 1,3 % 1,1 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,9 % 1],Just [1 % 1,9 % 1,27 % 1,27 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,25 % 1],Just [1 % 1,15 % 1,75 % 1,125 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,36 % 1],Just [1 % 1,18 % 1,108 % 1,216 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,64 % 1],Just [1 % 1,24 % 1,192 % 1,512 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,81 % 1],Just [1 % 1,27 % 1,243 % 1,729 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,121 % 1],Just [1 % 1,33 % 1,363 % 1,1331 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,169 % 1],Just [1 % 1,39 % 1,507 % 1,2197 % 1])
  ]

For simplicity, take numerator only.

  *GMulFin> map fst it
  [Just [0 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,1 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,9 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,25 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,36 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,64 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,81 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,121 % 1]
  ,Just [0 % 1,0 % 1,0 % 1,0 % 1,169 % 1]
  ]
  *GMulFin> sequence it
  Just [[0 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,1 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,9 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,25 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,36 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,64 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,81 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,121 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,169 % 1]
       ]

Technically, this transposeWith function is a key.

  *GMulFin> fmap (transposeWith (0%1)) it
  Just [[0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1]
       ,[0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1]
       ,[0 % 1,1 % 1,9 % 1,25 % 1,36 % 1,64 % 1,81 % 1,121 % 1,169 % 1]
       ]
  *GMulFin> fmap (map (zip [0,1,3,5,6,8,9,11,13])) it
  Just [[(0,0 % 1),(1,0 % 1),(3,0 % 1),(5,0 % 1),(6,0 % 1),(8,0 % 1),(9,0 % 1),(11,0 % 1),(13,0 % 1)]
       ,[(0,0 % 1),(1,0 % 1),(3,0 % 1),(5,0 % 1),(6,0 % 1),(8,0 % 1),(9,0 % 1),(11,0 % 1),(13,0 % 1)]
       ,[(0,0 % 1),(1,0 % 1),(3,0 % 1),(5,0 % 1),(6,0 % 1),(8,0 % 1),(9,0 % 1),(11,0 % 1),(13,0 % 1)]
       ,[(0,0 % 1),(1,0 % 1),(3,0 % 1),(5,0 % 1),(6,0 % 1),(8,0 % 1),(9,0 % 1),(11,0 % 1),(13,0 % 1)]
       ,[(0,0 % 1),(1,1 % 1),(3,9 % 1),(5,25 % 1),(6,36 % 1),(8,64 % 1),(9,81 % 1),(11,121 % 1),(13,169 % 1)]
       ]
  *GMulFin> it :: Maybe [Graph]
  Just [[(0 % 1,0 % 1),(1 % 1,0 % 1),(3 % 1,0 % 1),(5 % 1,0 % 1),(6 % 1,0 % 1),(8 % 1,0 % 1),(9 % 1,0 % 1),(11 % 1,0 % 1),(13 % 1,0 % 1)]
       ,[(0 % 1,0 % 1),(1 % 1,0 % 1),(3 % 1,0 % 1),(5 % 1,0 % 1),(6 % 1,0 % 1),(8 % 1,0 % 1),(9 % 1,0 % 1),(11 % 1,0 % 1),(13 % 1,0 % 1)]
       ,[(0 % 1,0 % 1),(1 % 1,0 % 1),(3 % 1,0 % 1),(5 % 1,0 % 1),(6 % 1,0 % 1),(8 % 1,0 % 1),(9 % 1,0 % 1),(11 % 1,0 % 1),(13 % 1,0 % 1)]
       ,[(0 % 1,0 % 1),(1 % 1,0 % 1),(3 % 1,0 % 1),(5 % 1,0 % 1),(6 % 1,0 % 1),(8 % 1,0 % 1),(9 % 1,0 % 1),(11 % 1,0 % 1),(13 % 1,0 % 1)]
       ,[(0 % 1,0 % 1),(1 % 1,1 % 1),(3 % 1,9 % 1),(5 % 1,25 % 1),(6 % 1,36 % 1),(8 % 1,64 % 1),(9 % 1,81 % 1),(11 % 1,121 % 1),(13 % 1,169 % 1)]
       ]

Then we can apply polynomial reconstruction for each "coefficient".

  *GMulFin> fmap (map uniPolCoeff) it
  Just [Just [0 % 1],Just [0 % 1],Just [0 % 1],Just [0 % 1],Just [0 % 1,0 % 1,1 % 1]]
  *GMulFin> fmap sequence it
  Just (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]])
  *GMulFin> Control.Monad.join it
  Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]

This means that the numerator has t^4, and it clealy is x^2*y^2.

  *GMulFin> map (\t -> ratFunc2Coeff (homogeneous wilFunc 1 t)) [0,1,3,5,6,8,9,11,13]
  [(Just [0 % 1],Just [1 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,1 % 1],Just [1 % 1,3 % 1,3 % 1,1 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,9 % 1],Just [1 % 1,9 % 1,27 % 1,27 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,25 % 1],Just [1 % 1,15 % 1,75 % 1,125 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,36 % 1],Just [1 % 1,18 % 1,108 % 1,216 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,64 % 1],Just [1 % 1,24 % 1,192 % 1,512 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,81 % 1],Just [1 % 1,27 % 1,243 % 1,729 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,121 % 1],Just [1 % 1,33 % 1,363 % 1,1331 % 1])
  ,(Just [0 % 1,0 % 1,0 % 1,0 % 1,169 % 1],Just [1 % 1,39 % 1,507 % 1,2197 % 1])
  ]
  *GMulFin> join . fmap (sequence . (map (uniPolCoeff . (zip [0,1,3,5,6,8,9,11,13]))) . (transposeWith (0%1))) . sequence . map fst $ it 
  Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]
  
> twoVariableRational
>   :: (Q -> Q -> Q) -- 2var function
>   -> [Q]           -- safe ys
>   -> (Maybe [[Ratio Int]], Maybe [[Ratio Int]]) 
> twoVariableRational f ys = (num, den)
>   where
>     num = helper fst
>     den = helper snd
>     helper third = join . fmap (mapM (uniPolCoeff . (zip ys)) 
>                    . transposeWith (0 % 1)) . mapM third $ gs
>     gs = map (\y -> ratFunc2Coeff (homogeneous f 1 y)) ys
>                    
> --  helper third = join . fmap (sequence . (map (uniPolCoeff . (zip ys))) 
> --                 . (transposeWith (0%1))) . sequence . map third $ gs
>
 
  GMulFin.lhs:139:35: Warning: Use mapM
  Found:
    sequence . (map (uniPolCoeff . (zip ys))) . (transposeWith (0 % 1))
  Why not:
    mapM (uniPolCoeff . (zip ys)) . transposeWith (0 % 1)

  *GMulFin> twoVariableRational wilFunc [0..10]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1,3 % 1],[0 % 1,0 % 1,3 % 1],[0 % 1,0 % 1,0 % 1,1 % 1]]
  )
  *GMulFin> twoVariableRational wilFunc [10..20]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1,3 % 1],[0 % 1,0 % 1,3 % 1],[0 % 1,0 % 1,0 % 1,1 % 1]]
  )
  *GMulFin> twoVariableRational wilFunc [1,3..21]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1,3 % 1],[0 % 1,0 % 1,3 % 1],[0 % 1,0 % 1,0 % 1,1 % 1]]
  )  
-- wilFunc x y = (x^2*y^2)/(1+y)^3
  
  *GMulFin> twoVariableRational wilFunc [1,2,4,6,9,11,13]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1,3 % 1],[0 % 1,0 % 1,3 % 1],[0 % 1,0 % 1,0 % 1,1 % 1]]
  )
  *GMulFin> twoVariableRational (\x y -> (x^3*y)/(1 + (x-y)^2)) [0..20]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,1 % 1]],*** Exception: newtonTriangleZp: need more evaluation
  CallStack (from HasCallStack):
    error, called at ./GUniFin.lhs:80:23 in main:GUniFin
  *GMulFin> twoVariableRational (\x y -> (x^3*y)/(1 + (x-y)^2)) [10..30]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1],[1 % 1,(-2) % 1,1 % 1],[0 % 1]]
  )
  
  *GMulFin> twoVariableRational (\x y -> (x^3*y)/(1 + (x-y)^2)) [1,3..9]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,1 % 1]],*** Exception: newtonTriangleZp: need more evaluation
  CallStack (from HasCallStack):
    error, called at ./GUniFin.lhs:80:23 in main:GUniFin
  *GMulFin> twoVariableRational (\x y -> (x^3*y)/(1 + (x-y)^2)) [1,3..11]
  (Just [[0 % 1],[0 % 1],[0 % 1],[0 % 1],[0 % 1,1 % 1]]
  ,Just [[1 % 1],[0 % 1],[1 % 1,(-2) % 1,1 % 1],[0 % 1]]
  )

--

> wilFunc2 x y = (x^4*y^2)*(1+y+y^2)^2 / (1+y)^4

