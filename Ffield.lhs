Ffield.lhs

https://arxiv.org/pdf/1608.01902.pdf

> module Ffield where

> import Data.Ratio 
> import Data.Maybe
> import Data.Numbers.Primes
> import System.Random

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
> a `inversep` p = let (_,x,_) = exGCD a p in
>   if isPrime p then Just (x `mod` p)
>                else Nothing
>
> inversesp :: Integral a => a -> [Maybe a]
> inversesp p = map (`inversep` p) [1..(p-1)]

A map from Q to Z_p.

> -- p should be prime.
> modp :: Integral a => Ratio a -> a -> a
> q `modp` p = (a * (bi `mod` p)) `mod` p
>   where
>     (a,b) = (numerator q, denominator q)
>     bi = fromJust (b `inversep` p)

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
> -- Hard code of big primes.
> -- bigPrimes :: [Int]
> bigPrimes = dropWhile (< (897473 :: Integer)) $ takeWhile (< (978948 :: Integer)) primes  
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
Test, we first tried to use QuickCheck, but it does not work.

> {-
> trial = do
>   n <- randomRIO (0,1000) :: IO Int
>   d <- randomRIO (1,1000) :: IO Int
>   putStrLn $ "input: " ++ show (n%d)
>   let knownData = zip (map (modp (n%d)) bigPrimes) bigPrimes
>   return $ matches3' $ map guess knownData
>   putStrLn $ show $ (n%d) == fst (matches3' $ map guess knownData)

Our choice of bigPrimes are sometimes fail:

  *Ffield> trial
  input: 895 % 922
  ^[[A^?^?*** Exception: Ffield.lhs:(224,3)-(226,37): Non-exhaustive patterns in function matches3'

> trial' = do
>   n <- randomRIO (0,1000) :: IO Int
>   d <- randomRIO (1,1000) :: IO Int
>   putStrLn $ "input: " ++ show (n%d)
>   let knownData = zip (map (modp (n%d)) bigger) bigger
>   return $ matches3' $ map guess knownData
>   putStrLn $ show (matches3' $ map guess knownData) 
>   putStrLn $ show $ (n%d) == fst (matches3' $ map guess knownData)
 
> bigger = dropWhile (<897473) primes
> -}
  
  *Ffield> trial'
  input: 125 % 399
  (125 % 399,897473)
  True
  (0.25 secs, 310,621,352 bytes)
  *Ffield> trial'
  input: 112 % 939
  (112 % 939,909383)
  True
  (0.40 secs, 378,062,424 bytes)
  *Ffield> trial'
  input: 297 % 391
  (297 % 391,897473)
  True
  (0.01 secs, 2,101,240 bytes)
  *Ffield> trial'
  input: 17 % 16
  (17 % 16,897473)
  True
  (0.01 secs, 2,103,728 bytes)
  *Ffield> trial'
  input: 125 % 102
  (125 % 102,897473)
  True
  (0.01 secs, 2,103,848 bytes)

--
Chinese Remeinder Theorem, and its usage

  *Ffield> let q = 895 % 922
  *Ffield> let knownData = zip (map (modp q) bigPrimes) bigPrimes 
  *Ffield> take 10 knownData 
  [(882873,897473),(365035,897497),(705735,897499),(511060,897517),(526641,897527),(30179,897553),(760296,897557),(330016,897563),(567554,897571),(137266,897577)]
  *Ffield> map guess it
  [((-854) % 123,897473),((-656) % 327,897497),((-192) % 805,897499),((-491) % 497,897517),((-856) % 121,897527),((-194) % 803,897553),((-161) % 837,897557),((-559) % 427,897563),((-493) % 495,897571),((-891) % 85,897577)]
  *Ffield> last knownData 
  (491598,978947)
  *Ffield> guess it
  ((-391) % 691,978947)
  *Ffield> uncurry exGCD' (882873,897473)
  ([0,1,60,2,8,20,1,4,1,6]
  ,[882873,897473,882873,14600,6873,854,41,34,7,6,1]
  ,[1,0,1,-1,61,-123,1045,-21023,22068,-109295,131363,-897473]
  ,[0,1,0,1,-60,121,-1028,20681,-21709,107517,-129226,882873]
  )

  *Ffield> take 2 knownData 
  [(882873,897473),(365035,897497)]
  *Ffield> let [(a1,p1),(a2,p2)] = it
  *Ffield> (inversep p2 p1)
  Just 261763
  *Ffield> 261763 `mod` p1
  f261763
  *Ffield> it * p2
  234931507211
  *Ffield> let m1 = it * p2
  *Ffield> m1
  210850322927350867
  *Ffield> let m1 = (261763 `mod` p1) * p2
  *Ffield> m1
  234931507211
  *Ffield> inversep p1 p2
  Just 635727
  *Ffield> let m2 = (635727 `mod` p2) * p1
  *Ffield> m2
  570547817871
  *Ffield> let a = m1*a1 + m2*a2
  *Ffield> a
  415684607262437688
  *Ffield> p1*p2
  805479325081
  *Ffield> a `mod` (p1*p2)
  86488560937
  *Ffield> guess (86488560937, p1*p2)
  (86488560937 % 1,805479325081)
  *Ffield> uncurry exGCD' (86488560937, p1*p2)
  ([0,9,3,5,6,976113,1,1,1,1,2,1,1,1,8,2]
  ,[86488560937,805479325081,86488560937,27082276648,5241730993,873621683,895,548,347,201,146,55,36,19,17,2,1]
  ,[1,0,1,-9,28,-149,922,-899976335,899977257,-1799953592,2699930849,-4499884441,11699699731,-16199584172,27899283903,-44098868075,380690228503,-805479325081]
  ,[0,1,0,1,-3,16,-99,96635203,-96635302,193270505,-289905807,483176312,-1256258431,1739434743,-2995693174,4735127917,-40876716510,86488560937]
  )
  *Ffield> uncurry exGCD' (a, p1*p2)
  ([516071,9,3,5,6,976113,1,1,1,1,2,1,1,1,8,2],[415684607262437688,805479325081,86488560937,27082276648,5241730993,873621683,895,548,347,201,146,55,36,19,17,2,1],[1,0,1,-9,28,-149,922,-899976335,899977257,-1799953592,2699930849,-4499884441,11699699731,-16199584172,27899283903,-44098868075,380690228503,-805479325081],[0,1,-516071,4644640,-14449991,76894595,-475817561,464451783814988,-464452259632549,928904043447537,-1393356303080086,2322260346527623,-6037876996135332,8360137342662955,-14398014338798287,22758151681461242,-196463227790488223,415684607262437688])

  *Ffield> let q = (895%922)
  *Ffield> let knownData = zip (map (modp q) bigPrimes) bigPrimes 
  *Ffield> let [(a1,p1),(a2,p2)] = take 2 knownData 
  *Ffield> it
  (895 % 922,805479325081)
  *Ffield> let m1 = (fromJust $ inversep p2 p1)*p2
  *Ffield> let m2 = (fromJust $ inversep p1 p2)*p1
  *Ffield> let a = (m1*a1 + m2*a2) `mod` (p1*p2)
  *Ffield> a
  86488560937
  *Ffield> let p = p1*p2
  *Ffield> p
  805479325081
  *Ffield> guess (a,p)
  (86488560937 % 1,805479325081)
  *Ffield> p
  805479325081
  *Ffield> a
  86488560937
  *Ffield> guess (86488560937, 805479325081)
  (895 % 922,805479325081)

> {-
> imagesAndPrimes :: Ratio Int -> [(Int, Int)]
> imagesAndPrimes q = zip (map (modp q) bigPrimes) bigPrimes
> -}

  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> let [(a1,p1),(a2,p2)] = take 2 knownData 
  *Ffield> take 2 knownData 
  [(882873,897473),(365035,897497)]
  *Ffield> map guess it
  [((-854) % 123,897473),((-656) % 327,897497)]
  
  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> let [(a1,p1),(a2,p2)] = take 2 knownData 
  *Ffield> let m1 = fromJust (p2 `inversep` p1) * p2
  *Ffield> let m2 = fromJust (p1 `inversep` p2) * p1
  *Ffield> let p = p1*p2
  *Ffield> let a = (m1*a1 + m2*a2) `mod` p
  *Ffield> a
  86488560937
  *Ffield> p
  805479325081
  *Ffield> guess (a,p)
  (86488560937 % 1,805479325081) -- overflow!

  *Ffield> guess (86488560937, 805479325081)
  (895 % 922,805479325081)

> crtRec' (a1,p1) (a2,p2) = (a,p)
>   where
>     a = (m1*a1 + m2*a2) `mod` p
>     m1 = fromJust (p2 `inversep` p1) * p2
>     m2 = fromJust (p1 `inversep` p2) * p1
>     p = p1*p2
>
> pile :: (a -> a -> a) -> [a] -> [a]
> pile f [] = []
> pile f dd@(d:ds) = d : zipWith f (pile f dd) ds

