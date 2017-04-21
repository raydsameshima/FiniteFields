Ffield.lhs

https://arxiv.org/pdf/1608.01902.pdf

> module Ffield where

> import Data.Ratio 
> import Data.Maybe
> import Data.Numbers.Primes
> import Test.QuickCheck

> -- Eucledian algorithm.
> myGCD :: Integral a => a -> a -> a
> myGCD a b
>   | b < 0 = myGCD a (-b)
> myGCD a b
>   | a == b = a
>   | b >  a = myGCD b a
>   | b <  a = myGCD (a-b) b

Consider a finite ring
  Z_n := [0..(n-1)]
of some Int number.
If any non-zero element has its multiplication inverse, 
then the ring is a field:

> -- Our target should be in Int.
> isField 
>   :: Int -> Bool
> isField = isPrime

Here we would like to implement the extended Euclidean algorithm.
See the algorithm, examples, and pseudo code at:

  https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
  http://qiita.com/bra_cat_ket/items/205c19611e21f3d422b7

> exGCD' 
>   :: (Integral n) => 
>      n -> n -> ([n], [n], [n], [n])
> exGCD' a b = (qs, rs, ss, ts)
>   where
>     qs = zipWith quot rs (tail rs)
>     rs = takeBefore (==0) r'
>     r' = steps a b
>     ss = steps 1 0
>     ts = steps 0 1
>
>     steps a b = rr
>       where 
>         rr@(_:rs) = a:b: zipWith (-) rr (zipWith (*) qs rs)
>
> takeBefore 
>   :: (a -> Bool) -> [a] -> [a]
> takeBefore p = foldr func []
>   where
>     func x xs 
>       | p x       = []
>       | otherwise = x : xs
>
> -- Bezout's identity a*x + b*y = gcd a b 
> exGCD 
>   :: Integral t => 
>      t -> t -> (t, t, t)
> exGCD a b = (g, x, y) 
>   where
>     (_,r,s,t) = exGCD' a b
>     g = last r
>     x = last . init $ s
>     y = last . init $ t

> -- We use built-in function gcd.
> coprime 
>   :: Integral a => 
>      a -> a -> Bool
> coprime a b = gcd a b == 1

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
> inversep' 
>   :: Int -> Int -> Int
> 0 `inversep'` _ = error "inversep': zero division"
> a `inversep'` p = (x `mod` p)
>   where 
>     (_,x,_) = exGCD a p
>
> -- Returns a list of inveres of given ring Z_p.
> inversesp 
>   :: Int -> [Maybe Int]
> inversesp p = map (`inversep` p) [1..(p-1)]
>
> -- A map from Q to Z_p, where p is a prime.
> modp 
>   :: Ratio Int -> Int -> Maybe Int
> q `modp` p 
>   | coprime b p = Just $ (a * (bi `mod` p)) `mod` p
>   | otherwise   = Nothing
>   where
>     (a,b)   = (numerator q, denominator q)
>     Just bi = b `inversep` p
>
> -- When the denominator of q is not proprtional to p, use this.
> modp' 
>   :: Ratio Int -> Int -> Int
> q `modp'` p = (a * (bi `mod` p)) `mod` p
>   where
>     (a,b)   = (numerator q, denominator q)
>     bi = b `inversep'` p
>
> -- This is guess function without Chinese Reminder Theorem.
> guess 
>   :: Integral t => 
>      (Maybe t, t)       -- (q `modp` p, p)
>   -> Maybe (Ratio t, t)
> guess (Nothing, _) = Nothing
> guess (Just a, p) = let (_,rs,ss,_) = exGCD' a p in
>   Just (select rs ss p, p)
>     where
>       select 
>         :: Integral t => 
>            [t] -> [t] -> t -> Ratio t
>       select [] _ _ = 0%1
>       select (r:rs) (s:ss) p
>         | s /= 0 && r*r <= p && s*s <= p = r%s
>         | otherwise                      = select rs ss p
>
> -- Hard code of big primes
> -- We have chosen a finite number (100) version.
> bigPrimes :: [Int]
> bigPrimes = take 100 $ dropWhile (<10^4) primes
> -- bigPrimes = take 100 $ dropWhile (< 10^6) primes

  *Ffield> bigPrimes 
  [10007,10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,10103
  ,10111,10133,10139,10141,10151,10159,10163,10169,10177,10181,10193,10211
  ,10223,10243,10247,10253,10259,10267,10271,10273,10289,10301,10303,10313
  ,10321,10331,10333,10337,10343,10357,10369,10391,10399,10427,10429,10433
  ,10453,10457,10459,10463,10477,10487,10499,10501,10513,10529,10531,10559
  ,10567,10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,10663
  ,10667,10687,10691,10709,10711,10723,10729,10733,10739,10753,10771,10781
  ,10789,10799,10831,10837,10847,10853,10859,10861,10867,10883,10889,10891
  ,10903,10909,10937,10939
  ]

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
 
> imagesAndPrimes 
>   :: Ratio Int -> [(Maybe Int, Int)]
> imagesAndPrimes q = zip (map (modp q) bigPrimes) bigPrimes

  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> let [(a1,p1),(a2,p2)] = take 2 knownData
  *Ffield> take 2 knownData
  [(Just 6003,10007),(Just 9782,10009)]
  *Ffield> map guess it
  [Just ((-6) % 5,10007),Just (21 % 44,10009)]
 
Our data is a list of the type
  [(Maybe Int, Int)]
In order to use CRT, we should cast its type.

> toInteger2 
>   :: [(Maybe Int, Int)] -> [(Maybe Integer, Integer)]
> toInteger2 = map helper
>   where 
>     helper (x,y) = (fmap toInteger x, toInteger y)
>
> crtRec'
>   :: Integral a => 
>      (Maybe a, a) -> (Maybe a, a) -> (Maybe a, a)
> crtRec' (Nothing,p) (_,q)       = (Nothing, p*q)
> crtRec' (_,p)       (Nothing,q) = (Nothing, p*q)
> crtRec' (Just a1,p1) (Just a2,p2) = (Just a,p)
>   where
>     a = (a1*p2*m2 + a2*p1*m1) `mod` p
>     Just m1 = p1 `inversep` p2 
>     Just m2 = p2 `inversep` p1
>     p = p1*p2
>
> matches3 
>   :: Eq a => 
>      [Maybe (a,b)] -> Maybe (a,b)
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

We should determine the number of matches to cover the range of machine size
Integer, i.e., Int of Haskell.

  *Ffield> let mI = maxBound :: Int
  *Ffield> mI == 2^63-1
  True
  *Ffield> logBase 10 (fromIntegral mI)
  18.964889726830812

Since our choice of bigPrimes are
  O(10^4)
5 times is enough to cover the machine size integers.

> reconstruct 
>   :: [(Maybe Int, Int)] -> Maybe (Ratio Integer)
> -- reconstruct = matches 10 . makeList -- 10 times match
> reconstruct = matches 5 . makeList -- 5 times match
>   where
>     matches n (a:as)
>       | all (a==) $ take (n-1) as = a
>       | otherwise                     = matches n as
>
>     makeList = map (fmap fst . guess) . scanl1 crtRec' . toInteger2 
>                . filter (isJust . fst)
>
> -- cast version
> reconstruct' 
>   :: [(Maybe Int, Int)] -> Maybe (Ratio Int)
> reconstruct' = fmap coersion . reconstruct
>   where
>     coersion :: Ratio Integer -> Ratio Int
>     coersion q = (fromInteger . numerator $ q) 
>                    % (fromInteger . denominator $ q)
  
  *Ffield> let q = 895%922
  *Ffield> let knownData = imagesAndPrimes q
  *Ffield> reconstruct knownData 
  Just (895 % 922)

-- QuickCheck

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
