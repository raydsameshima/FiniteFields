Ffield.lhs

https://arxiv.org/pdf/1608.01902.pdf

> module Ffield where

> import Data.Ratio 
> import Data.Maybe
> import Data.Numbers.Primes
>
> import System.Random
> import Test.QuickCheck

> coprime :: Integral a => a -> a -> Bool
> coprime a b = gcd a b == 1

Consider a finite ring
  Z_n := [0..(n-1)]
If any non-zero element has its multiplication inverse, then the ring is a field:

> isField :: Integral a => a -> Bool
> isField = isPrime

Here we would like to implement the extended Euclidean algorithm.
See the algorithm, examples, and pseudo code at:

  https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
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
>       where 
>         rr@(_:rs) = a:b: zipWith (-) rr (zipWith (*) qs rs)
>
> takeUntil :: (a -> Bool) -> [a] -> [a]
> takeUntil p = foldr func []
>   where
>     func x xs 
>       | p x = []
>       | otherwise = x : xs

This example is from wikipedia:

  *Ffield> exGCD' 240 46
  ([5,4,1,1,2]
  ,[240,46,10,6,4,2]
  ,[1,0,1,-4,5,-9,23]
  ,[0,1,-5,21,-26,47,-120]
  )
  *Ffield> gcd 240 46
  2
  *Ffield> 240*(-9) + 46*(47)
  2

> -- a*x + b*y = gcd a b
> exGCD :: Integral t => t -> t -> (t, t, t)
> exGCD a b = (g, x, y)
>   where
>     (_,r,s,t) = exGCD' a b
>     g = last r
>     x = last . init $ s
>     y = last . init $ t

Example Z_{11}

  *Ffield> isField 11
  True
  *Ffield> map (exGCD 11) [0..10]
  [(11,1,0),(1,0,1),(1,1,-5),(1,-1,4),(1,-1,3)
  ,(1,1,-2),(1,-1,2),(1,2,-3),(1,3,-4),(1,-4,5),(1,1,-1)
  ]

  *Ffield> map ((`mod` 11) . (\(_,_,x)->x) . exGCD 11) [1..10] 
  [1,6,4,3,9,2,8,7,5,10]

We get non-zero elements with its inverse:

  *Ffield> zip [1..10] it
  [(1,1),(2,6),(3,4),(4,3),(5,9),(6,2),(7,8),(8,7),(9,5),(10,10)]
     
> -- a^{-1} (in Z_p) == a `inversep` p
> inversep :: Integral a => a -> a -> Maybe a
> a `inversep` p = let (g,x,_) = exGCD a p in
>   if (g == 1) then Just (x `mod` p) -- g==1 <=> coprime a p
>               else Nothing
>
> inversesp :: Integral a => a -> [Maybe a]
> inversesp p = map (`inversep` p) [1..(p-1)]

A map from Q to Z_p.

> -- p should be prime.
> modp :: Integral a => Ratio a -> a -> a
> q `modp` p = (a * (bi `mod` p)) `mod` p
>   where
>     (a,b) = (numerator q, denominator q)
>     bi    = fromJust (b `inversep` p)

Example: on Z_{11}
Consider (3 % 7).

  *Ffield> let q = 3%7
  *Ffield> 3 `mod` 11
  3
  *Ffield> 7 `inversep` 11
  Just 8
  *Ffield> fromJust it * 3 `mod` 11
  2

On Z_{11}, (7^{-1} `mod` 11) is equal to (8 `mod` 11) and

  (3%7) |-> (3 * (7^{-1} `mod` 11) `mod` 11)
             == (3*8 `mod` 11) 
             == 2 ` mod 11

  *Ffield> (3%7) `modp` 11
  2

Example of reconstruction Z_p -> Q

  *Ffield> let q = (1%3)
  *Ffield> take 3 $ dropWhile (<100) primes
  [101,103,107]
  *Ffield> q `modp` 101
  34
  *Ffield> let rec x = exGCD' (q `modp` x) x
  *Ffield> rec 101
  ([0,2,1,33],[34,101,34,33,1],[1,0,1,-2,3,-101],[0,1,0,1,-1,34])
  *Ffield> rec 103
  ([0,1,2,34],[69,103,69,34,1],[1,0,1,-1,3,-103],[0,1,0,1,-2,69])
  *Ffield> rec 107
  ([0,2,1,35],[36,107,36,35,1],[1,0,1,-2,3,-107],[0,1,0,1,-1,36])  

> -- This is guess function without Chinese Reminder Theorem.
> guess :: Integral t => 
>          (t, t)       -- (q `modp` p, p)
>       -> (Ratio t, t)
> guess (a, p) = let (_,rs,ss,_) = exGCD' a p in
>   (select rs ss p, p)
>     where
>       select :: Integral t => [t] -> [t] -> t -> Ratio t
>       select [] _ _ = 0%1
>       select (r:rs) (s:ss) p
>         | s /= 0 && r*r <= p && s*s <= p = r%s
>         | otherwise = select rs ss p
>
> -- Hard code of big primes
> -- For chinese reminder theorem we declare it as [Integer].
> bigPrimes :: (Integral a) => [a]
> bigPrimes = dropWhile (< 897473) $ takeWhile (< 978948) primes  
>
> matches3 :: Eq a => [a] -> a
> matches3 (a:bb@(b:c:cs))
>   | a == b && b == c = a
>   | otherwise        = matches3 bb

What we know is a list of (q `modp` p) and prime p.

  *Ffield> let q = 10%19
  *Ffield> let knownData = zip (map (modp q) bigPrimes) bigPrimes  
  *Ffield> matches3 $  map (fst . guess) knownData 
  10 % 19

> -- This function does not use CRT, so it can fail (O(10^3)%O(10^3)).
> reconstruct :: Integral a => 
>                [(a, a)]  -- :: [(Z_p, primes)]
>             -> Ratio a
> reconstruct aps = matches3 $ map (fst . guess) aps

Here is a naive test:
  > let qs = [1 % 3, 10 % 19, 41 % 17, 30 % 311, 311 % 32
             ,869 % 232, 778 % 123, 331 % 739
             ]
  > let func q = zip (map (modp q) bigPrimes) bigPrimes 
  > let longList = map func qs
  > map reconstruct longList 
  [1 % 3,10 % 19,41 % 17,30 % 311,311 % 32
  ,869 % 232,778 % 123,331 % 739
  ]
  > it == qs
  True

> matches3' :: Eq a => [(a, t)] -> (a, t)
> matches3' (a0@(a,_):bb@((b,_):(c,_):cs))
>   | a == b && b == c = a0
>   | otherwise        = matches3' bb

  *Ffield> let q = (331%739)
  (0.01 secs, 48,472 bytes)
  *Ffield> let knownData = zip (map (modp q) primes) primes
  (0.02 secs, 39,976 bytes)
  *Ffield> matches3' $ map guess knownData 
  (331 % 739,614693)
  (19.92 secs, 12,290,852,136 bytes)

--
Chinese Remainder Theorem, and its usage
 
> imagesAndPrimes :: (Integral b, Integral a) => Ratio a -> [(a, b)]
> imagesAndPrimes q = zip (map (modp q) bigPrimes) bigPrimes

  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> let [(a1,p1),(a2,p2)] = take 2 knownData 
  *Ffield> take 2 knownData 
  [(882873,897473),(365035,897497)]
  *Ffield> map guess it
  [((-854) % 123,897473),((-656) % 327,897497)]
  
> crtRec' :: Integral t => (t, t) -> (t, t) -> (t, t)
> crtRec' (a1,p1) (a2,p2) = (a,p)
>   where
>     a = (a1*p2*m2 + a2*p1*m1) `mod` p
>     m1 = fromJust (p1 `inversep` p2) 
>     m2 = fromJust (p2 `inversep` p1)
>     p = p1*p2
>
> pile :: (a -> a -> a) -> [a] -> [a]
> pile f [] = []
> pile f dd@(d:ds) = d : zipWith' f (pile f dd) ds
> -- pile f [d0,d1,d2 ..] == [d0, (f d0 d1), (f (f d0 d1) d2) ..]
>
> -- Strict zipWith, from:
> --   http://d.hatena.ne.jp/kazu-yamamoto/touch/20100624/1277348961
> zipWith' :: (a -> b -> c) -> [a] -> [b] -> [c]
> zipWith' f (a:as) (b:bs) = (x `seq` x) : zipWith' f as bs
>   where x = f a b
> zipWith' _ _      _      = []

  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> take 4 knownData 
  [(882873,897473)
  ,(365035,897497)
  ,(705735,897499)
  ,(511060,897517)
  ]
  *Ffield> pile crtRec' it
  [(882873,897473)
  ,(86488560937,805479325081)
  ,(397525881357811624,722916888780872419)
  ,(232931448259966259937614,648830197267942270883623)
  ]
  *Ffield> map guess it
  [((-854) % 123,897473)
  ,(895 % 922,805479325081)
  ,(895 % 922,722916888780872419)
  ,(895 % 922,648830197267942270883623)
  ]

  *Ffield> reconstruct knownData'
  895 % 922
  *Ffield> let knownData' = pile crtRec' knownData 
  *Ffield> matches3' $ map guess knownData'
  (895 % 922,805479325081)

> recCRT :: Integral a => [(a,a)] -> Ratio a
> recCRT = reconstruct . pile crtRec'

> recCRT' = matches3' . map guess . pile crtRec'

  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> recCRT knownData 
  895 % 922
  *Ffield> recCRT' knownData 
  (895 % 922,805479325081)

--
todo: use QuickCheck

> trial = do
>   n <- randomRIO (0,10000) :: IO Integer
>   d <- randomRIO (1,10000) :: IO Integer
>   let q = (n%d)
>   putStrLn $ "input: " ++ show q
>   return $ recCRT' . imagesAndPrimes $ q

  *Ffield> trial
  input: 1080 % 6931
  (1080 % 6931,805479325081)
  *Ffield> trial
  input: 2323 % 1248
  (2323 % 1248,805479325081)
  *Ffield> trial
  input: 6583 % 1528
  (6583 % 1528,805479325081)
  *Ffield> trial
  input: 721 % 423
  (721 % 423,897473)
  *Ffield> trial
  input: 9967 % 7410
  (9967 % 7410,805479325081)
