Ffield.lhs

https://arxiv.org/pdf/1608.01902.pdf

> module Ffield where

> import Data.Ratio 
> import Data.Maybe
> import Data.Numbers.Primes
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
>
> -- a*x + b*y = gcd a b
> exGCD 
>   :: Integral t => 
>      t -> t -> (t, t, t)
> exGCD a b = (g, x, y)
>   where
>     (_,r,s,t) = exGCD' a b
>     g = last r
>     x = last . init $ s
>     y = last . init $ t
>
> -- a^{-1} (in Z_p) == a `inversep` p
> inversep 
>   :: Integral a => 
>      a -> a -> Maybe a -- We also use in CRT.
> a `inversep` p = let (g,x,_) = exGCD a p in
>   if (g == 1) 
>     then Just (x `mod` p) -- g==1 <=> coprime a p
>     else Nothing
>
> -- If a is "safe" value, we can use this.
> inversep' :: Int -> Int -> Int
> a `inversep'` p = (x `mod` p)
>   where (_,x,_) = exGCD a p
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
> -- When the denominator of q does not factor p, use this.
> modp' :: Ratio Int -> Int -> Int
> q `modp'` p = (a * (bi `mod` p)) `mod` p
>   where
>     (a,b)   = (numerator q, denominator q)
>     bi = b `inversep'` p

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
> -- We have chosen a finite number version.
> bigPrimes :: [Int]
> bigPrimes = take 100 $ dropWhile (<10^4) primes
> -- bigPrimes = dropWhile (<10^4) primes

  *Ffield> let knownData q = zip (map (modp q) bigPrimes) bigPrimes
  *Ffield> let ds = knownData (12%13)
  *Ffield> map guess ds
  [Just (12 % 13,10007)
  ,Just (12 % 13,10009)
  ,Just (12 % 13,10037)
  ,Just (12 % 13,10039) ..

  *Ffield> let ds = knownData (112%113)
  *Ffield> map guess ds
  [Just ((-39) % 50,10007)
  ,Just ((-41) % 48,10009)
  ,Just ((-69) % 20,10037)
  ,Just ((-71) % 18,10039) ..

--

Chinese Remainder Theorem, and its usage
 
> imagesAndPrimes :: Ratio Int -> [(Maybe Int, Int)]
> imagesAndPrimes q = zip (map (modp q) bigPrimes) bigPrimes
 
Our data is a list of the type
  [(Maybe Int, Int)]
In order to use CRT, we should cast its type.

> toInteger2 :: [(Maybe Int, Int)] -> [(Maybe Integer, Integer)]
> toInteger2 = map helper
>   where 
>     helper (x,y) = (fmap toInteger x, toInteger y)
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
>
> matches3 :: Eq a => [Maybe (a,b)] -> Maybe (a,b)
> matches3 (b1@(Just (q1,p1)):bb@((Just (q2,_)):(Just (q3,_)):_))
>   | q1==q2 && q2==q3 = b1
>   | otherwise        = matches3 bb
> matches3 _ = Nothing

  *Ffield> let ds = imagesAndPrimes (1123%1135)
  *Ffield> map guess ds
  [Just (25 % 52,10007)
  ,Just ((-81) % 34,10009)
  ,Just ((-88) % 63,10037) ..

  *Ffield> matches3 it
  Nothing

  *Ffield> scanl1 crtRec' ds

  *Ffield> scanl1 crtRec' . toInteger2 $ ds
  [(Just 3272,10007)
  ,(Just 14913702,100160063)
  ,(Just 298491901442,1005306552331) ..

  *Ffield> map guess it
  [Just (25 % 52,10007)
  ,Just (1123 % 1135,100160063)
  ,Just (1123 % 1135,1005306552331)
  ,Just (1123 % 1135,10092272478850909) ..

  *Ffield> matches3 it
  Just (1123 % 1135,100160063)

The final reconstruction function takes a list of Z_p values and returns the three times matched guess.

> -- reconstruct :: [(Maybe Int, Int)] -> Maybe (Ratio Integer, Integer)
> -- reconstruct = matches3 . map guess . scanl1 crtRec' . toInteger2 
> -- reconstruct = matches3 . map guess . scanl1 crtRec' . filter (isJust . fst) . toInteger2 

  *Ffield> let q = 10937 % 10939
  *Ffield> let ds = imagesAndPrimes q
  *Ffield> map (fmap fst . guess) . scanl1 crtRec' . toInteger2 . filter (isJust . fst) $ ds
  [Just ((-19) % 24)
  ,Just ((-4977) % 4180)
  ,Just (10937 % 10939)
  ,Just (10937 % 10939)
  ,Just (10937 % 10939)
  ,Just (10937 % 10939)
  ..
  ]

> reconstruct :: [(Maybe Int, Int)] -> Maybe (Ratio Integer)
> reconstruct = matches 10 . makeList -- 10 times match
>   where
>     matches n (a:as)
>       | all (a==) $ take (n-1) as = a
>       | otherwise                     = matches n as
>     makeList = map (fmap fst . guess) . scanl1 crtRec' . toInteger2 . filter (isJust . fst)

--

  *Ffield> let q = 513197683989569 % 1047805145658 :: Ratio Int
  *Ffield> let ds = imagesAndPrimes q
  *Ffield> let answer = fmap fromRational . reconstruct $ ds
  *Ffield> answer :: Maybe (Ratio Int)
  Just (513197683989569 % 1047805145658)

> prop_rec :: Ratio Int -> Bool
> prop_rec q = Just q == answer
>   where
>    answer :: Maybe (Ratio Int)
>    answer = fmap fromRational . reconstruct $ ds
>    ds = imagesAndPrimes q

  *Ffield> quickCheckWith stdArgs { maxSuccess = 100000 } prop_rec 
  +++ OK, passed 100000 tests.


