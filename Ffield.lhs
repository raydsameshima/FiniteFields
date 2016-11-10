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
of some Int number.
If any non-zero element has its multiplication inverse, then the ring is a field:

> -- Our target should be in Int.
> isField :: Int -> Bool
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
>       | p x       = []
>       | otherwise = x : xs
> -- a*x + b*y = gcd a b
> exGCD :: Integral t => t -> t -> (t, t, t)
> exGCD a b = (g, x, y)
>   where
>     (_,r,s,t) = exGCD' a b
>     g = last r
>     x = last . init $ s
>     y = last . init $ t
>
> -- a^{-1} (in Z_p) == a `inversep` p
> inversep :: Integral a => a -> a -> Maybe a -- We also use in CRT.
> a `inversep` p = let (g,x,_) = exGCD a p in
>   if (g == 1) 
>     then Just (x `mod` p) -- g==1 <=> coprime a p
>     else Nothing
>
> inversesp :: Int -> [Maybe Int]
> inversesp p = map (`inversep` p) [1..(p-1)]
>
> -- A map from Q to Z_p.
> -- p should be prime.
> modp :: Ratio Int -> Int -> Maybe Int
> q `modp` p 
>   | coprime b p = Just $ (a * (bi `mod` p)) `mod` p
>   | otherwise   = Nothing
>   where
>     (a,b)   = (numerator q, denominator q)
>     Just bi = b `inversep` p
>
> -- This is guess function without Chinese Reminder Theorem.
> guess :: Integral t => 
>          (Maybe t, t)       -- (q `modp` p, p)
>       -> Maybe (Ratio t, t)
> guess (Nothing, _) = Nothing
> guess (Just a, p) = let (_,rs,ss,_) = exGCD' a p in
>   Just (select rs ss p, p)
>     where
>       select :: Integral t => [t] -> [t] -> t -> Ratio t
>       select [] _ _ = 0%1
>       select (r:rs) (s:ss) p
>         | s /= 0 && r*r <= p && s*s <= p = r%s
>         | otherwise = select rs ss p
>
> -- Hard code of big primes
> bigPrimes :: [Int]
> bigPrimes = take 10 $ dropWhile (<10^6) primes
> -- 10 is just for "practical" reason.

  *Ffield> let knownData q = zip (map (modp q) bigPrimes) bigPrimes
  *Ffield> let ds = knownData (112%113)
  *Ffield> map guess ds
  [Just (112 % 113,1000003)
  ,Just (112 % 113,1000033)
  ,Just (112 % 113,1000037) ..

  *Ffield> let ds' = knownData (895%922)
  *Ffield> map guess ds'
  [Just ((-309) % 799,1000003)
  ,Just ((-509) % 593,1000033)
  ,Just ((-476) % 627,1000037)
  ,Just ((-907) % 183,1000039)
  ,Just (895 % 922,1000081)
  ,Just ((-412) % 693,1000099)
  ,Just ((-711) % 385,1000117)
  ,Just ((-678) % 419,1000121)
  ,Just ((-579) % 521,1000133)
  ,Just ((-878) % 213,1000151)
  ]

--
Chinese Remainder Theorem, and its usage
 
> imagesAndPrimes :: Ratio Int -> [(Maybe Int, Int)]
> imagesAndPrimes q = zip (map (modp q) bigPrimes) bigPrimes

  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> let [(a1,p1),(a2,p2)] = take 2 knownData 
  *Ffield> take 2 knownData 
  [(882873,897473),(365035,897497)]
  *Ffield> map guess it
  [((-854) % 123,897473),((-656) % 327,897497)]
  
Our data is a list of the type
  [(Maybe Int, Int)]
In order to use CRT, we should cast its type.

> toInteger2 :: [(Maybe Int, Int)] -> [(Maybe Integer, Integer)]
> toInteger2 = map helper
>   where 
>     helper (x,y) = (fmap toInteger x, toInteger y)
>
> pile :: (a -> a -> a) -> [a] -> [a]
> pile f [] = []
> pile f dd@(d:ds) = d : zipWith f (pile f dd) ds
>
> crtRec':: Integral a => (Maybe a, a) -> (Maybe a, a) -> (Maybe a, a)
> crtRec' (Nothing,p) (_,q)       = (Nothing, p*q)
> crtRec' (_,p)       (Nothing,q) = (Nothing, p*q)
> crtRec' (Just a1,p1) (Just a2,p2) = (Just a,p)
>   where
>     a = (a1*p2*m2 + a2*p1*m1) `mod` p
>     Just m1 = p1 `inversep` p2 
>     Just m2 = p2 `inversep` p1
>     p = p1*p2

  *Ffield> let knownData q = zip (map (modp q) bigPrimes) bigPrimes
  *Ffield> let ds = knownData (1123%1135)
  *Ffield> let dsI = toInteger2 ds

  *Ffield> pile crtRec' dsI
  [(Just 294275,1000003)
  ,(Just 451998650266,1000036000099)
  ,(Just 386812376764943268,1000073001431003663) ..

  *Ffield> map guess it
  [Just ((-138) % 751,1000003)
  ,Just (1123 % 1135,1000036000099)
  ,Just (1123 % 1135,1000073001431003663)
  ,Just (1123 % 1135,1000112004278059472142857) ..

> -- Here is super auxiliary function.
> matches3 :: Eq a => [Maybe (a,b)] -> Maybe (a,b)
> matches3  (b1@(Just (q1,p1)):bb@((Just (q2,_)):(Just (q3,_)):_))
>   | q1==q2 && q2==q3 = b1
>   | otherwise        = matches3 bb
> matches3 _ = Nothing
 


> {-
> -- Before we apply matches3, (filter isJust)
> matches3 (a0@(a,_):bb@((b,_):(c,_):cs)) 
>   | a == b && b == c = a0
>   | otherwise        = matches3 bb


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

> -}
