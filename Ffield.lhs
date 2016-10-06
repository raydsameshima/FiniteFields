Ffield.lhs

https://arxiv.org/pdf/1608.01902.pdf

> module Ffield where

> import Data.Ratio 
> import Data.Maybe
> import Data.Numbers.Primes

> coprime :: Integral a => a -> a -> Bool
> coprime a b = gcd a b == 1

Consider a finite ring
  Z_n := [0..(n-1)]

> haveInverse :: Integral a => a -> [Bool]
> haveInverse n = map (coprime n) [0..(n-1)]

  *Ffield> haveInverse 8
  [False,True,False,True,False,True,False,True]
  *Ffield> zip [0..] $ haveInverse 8
  [(0,False),(1,True),(2,False),(3,True),(4,False),(5,True),(6,False),(7,True)]

If any non-zero element has its multiplication inverse, then the ring is a field:

> isField' :: Integral a => a -> Bool
> isField' n = and $ tail $ haveInverse n

Or more efficiently,

> isField :: Integral a => a -> Bool
> isField = isPrime

  zip [2..] $ map isField [2..13]
  [(2,True),(3,True),(4,False),(5,True),(6,False),(7,True),(8,False),(9,False),(10,False),(11,True),(12,False),(13,True)]

Here we would like to implement the extended Euclidean algorithm.
See the algorithm, examples, and pseudo code at:

  https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

I've asked at Qiita and get some solutions:

  http://qiita.com/bra_cat_ket/items/205c19611e21f3d422b7

> exGCD' :: (Integral n) => n -> n -> ([n], [n], [n], [n])
> exGCD' a b = (qs, rs, ss, ts)
>   where
>     qs = zipWith quot rs (tail rs)
>     rs = takeUntil (==0) r'
>     r' = steps a b
>     ss = steps 1 0
>     ts = steps 0 1
>     steps a b = rr
>       where rr@(_:rs) = a:b: zipWith (-) rr (zipWith (*) qs rs)
>
> takeUntil :: (a -> Bool) -> [a] -> [a]
> takeUntil p = foldr func []
>   where
>     func x xs 
>       | p x = []
>       | otherwise = x : xs

This example is from wikipedia:

  *Ffield> exGCD' 240 46
  ([5,4,1,1,2],[240,46,10,6,4,2],[1,0,1,-4,5,-9,23],[0,1,-5,21,-26,47,-120])
  *Ffield> gcd 240 46
  2
  *Ffield> 240*(-9) + 46*(47)
  2

> -- a*x + b*y = gcd a b
> exGcd a b = (g, x, y)
>   where
>     (_,r,s,t) = exGCD' a b
>     g = last r
>     x = last . init $ s
>     y = last . init $ t

  *Ffield> exGcd 46 240
  (2,47,-9)
  *Ffield> 46*47 + 240*(-9)
  2
  *Ffield> gcd 46 240
  2

Example Z_{11}

  *Ffield> isField 11
  True
  *Ffield> map (exGcd 11) [0..10]
  [(11,1,0),(1,0,1),(1,1,-5),(1,-1,4),(1,-1,3),(1,1,-2),(1,-1,2),(1,2,-3),(1,3,-4),(1,-4,5),(1,1,-1)]

  *Ffield> map ((`mod` 11) . (\(_,_,x)->x) . exGcd 11) [1..10] 
  [1,6,4,3,9,2,8,7,5,10]
  *Ffield> zip [1..10] it
  [(1,1),(2,6),(3,4),(4,3),(5,9),(6,2),(7,8),(8,7),(9,5),(10,10)]

> inverses :: Int -> Maybe [(Int, Int)]
> inverses n
>   | isPrime n = Just lst -- isPrime n
>   | otherwise = Nothing
>   where
>     lst' = map ((`mod` n) . (\(_,_,c)->c) . exGcd n) [1..(n-1)]
>     lst = zip [1..] lst'
>     
> inversep :: Int -> Int -> Maybe Int
> inversep p a = do
>   l <- inverses p
>   let a' = a `mod` p
>   return $ snd $ l !! (a'-1)
   
  *Ffield> inverses 11
  Just [(1,1),(2,6),(3,4),(4,3),(5,9),(6,2),(7,8),(8,7),(9,5),(10,10)]
  
A map from Q to Z_p.

> -- modp :: Int -> Ratio Int -> Int
> modp :: Ratio Int -> Int -> Int
> q `modp` p = (a * (bi `mod` p)) `mod` p
>   where
>     (a,b) = (numerator q, denominator q)
>     bi = fromJust $ inversep p b

Example: on Z_{11}
Consider (3 % 7).

  *Ffield Data.Ratio> let q = 3 % 7
  *Ffield Data.Ratio> 3 `mod` 11
  3
  *Ffield Data.Ratio> 7 `mod` 11
  7
  *Ffield Data.Ratio> inverses 11
  Just [(1,1),(2,6),(3,4),(4,3),(5,9),(6,2),(7,8),(8,7),(9,5),(10,10)]
  *Ffield Data.Ratio> 7*8 == 11*5+1
  True

on Z_{11}, (7^{-1} `mod` 11) is equal to (8 `mod` 11) and
  (3%7) |-> (3 * (7^{-1} `mod` 11) `mod` 11)
             == (3*8 `mod` 11) 
             == 2 ` mod 11

  *Ffield Data.Ratio> modp q 11
  2

Example: on Z_{5}
  *Ffield Data.Ratio> 3 `mod` 5
  3
  *Ffield Data.Ratio> 7 `mod` 5
  2
  *Ffield Data.Ratio> inverses 5
  Just [(1,1),(2,3),(3,2),(4,4)]
  *Ffield Data.Ratio> modp q 5 
  4

Functional reconstruction

Here is a nice imprementation for Thiele's interpolation formula:
https://rosettacode.org/wiki/Thiele%27s_interpolation_formula#Haskell
