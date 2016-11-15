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
> bigPrimes = dropWhile (<10^4) primes
> -- bigPrimes = take 10 $ dropWhile (<10^4) primes
> -- 10 is just for "practical" reason.

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
> matches3  (b1@(Just (q1,p1)):bb@((Just (q2,_)):(Just (q3,_)):_))
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

> reconstruct :: [(Maybe Int, Int)] -> Maybe (Ratio Integer, Integer)
> reconstruct = matches3 . map guess . scanl1 crtRec' . toInteger2 

--

> aux :: [(Maybe Int, Int)] -> Ratio Int
> aux = fromRational . fst . fromJust . reconstruct
>
> prop_rec :: Ratio Int -> Bool
> prop_rec q = q == aux ds
>   where
>     ds = imagesAndPrimes q

  *Ffield> quickCheck prop_rec 
  +++ OK, passed 100 tests.

  *Ffield> verboseCheck prop_rec 
  Passed:  
  0 % 1
  Passed: 
  1329564279259 % 689966835691
  Passed:  
  (-4749106441063) % 5354272769957
  Passed:  
  24393617614032 % 3429811214201
  Passed:  
  11317350866338 % 2376141260619
  Passed:  
  (-39223552414434) % 7602900239041
  Passed:  
  11646784360896 % 1038773197715
  Passed:  
  (-49349171845055) % 3752783115927
  Passed:  
  (-2655767154053) % 2023542479048
  Passed:  
  46861895449961 % 9804154238025
  Passed:   
  1432195821037 % 1223607419737
  Passed:   
  36703237585644 % 1206154505383
  Passed:   
  85403617062851 % 2493712836857
  Passed:   
  (-932690014756) % 2008903765881
  Passed:   
  (-1361504657430) % 5532549345619
  Passed:   
  (-29205487754653) % 2287912091229
  Passed:   
  (-41868516551014) % 4135899903201
  Passed:   
  (-165759270559576) % 2472941838409
  Passed:   
  133806538589611 % 4499617533687
  Passed:   
  (-17620738360931) % 6240948007885
  Passed:   
  39459205935 % 497205472096
  Passed:   
  117939084165907 % 6428931743209
  Passed:   
  (-39205881063583) % 914594476839
  Passed:   
  85731534523169 % 177493065396
  Passed:   
  22844351052383 % 625751286364
  Passed:   
  (-41929177046869) % 7717513998384
  Passed:   
  (-239104997136076) % 6797515055957
  Passed:   
  5217952768959 % 58314961636
  Passed:   
  (-16757243373579) % 51622205471
  Passed:   
  (-20666882197275) % 1468620478768
  Passed:   
  (-136479419289527) % 6830615084400
  Passed:   
  240244445293065 % 3788154937217
  Passed:   
  (-33372393134923) % 371591225618
  Passed:   
  (-101391102970683) % 1939379047834
  Passed:   
  (-309270804437319) % 8714166158161
  Passed:   
  51669639177678 % 6275195507707
  Passed:   
  (-55827295694032) % 179632819809
  Passed:   
  (-114812798039573) % 700296073152
  Passed:   
  (-178448558145370) % 7607296029601
  Passed:   
  23572878370179 % 1057593174922
  Passed:   
  20398146164919 % 2939076586582
  Passed:   
  (-291447760874087) % 185244184723
  Passed:   
  329433704076007 % 1794485745740
  Passed:   
  (-336625542209231) % 1599122245520
  Passed:   
  (-19873460891647) % 1422075885668
  Passed:   
  (-152680242051531) % 4981173771512
  Passed:   
  (-424983220468388) % 7464862254053
  Passed:   
  (-7574293347638) % 1657127826349
  Passed:   
  (-203671927794075) % 6166500691591
  Passed:   
  (-152003488487855) % 4130192768847
  Passed:   
  (-63707410228315) % 4680142309544
  Passed:   
  (-10545640427299) % 2918419691242
  Passed:   
  (-53037754892432) % 2711758675081
  Passed:   
  213538441367933 % 1513965455185
  Passed:   
  (-431088172705814) % 1197611173437
  Passed:   
  55694987217273 % 1329086881162
  Passed:   
  (-39053206137251) % 861738484033
  Passed:   
  135704343771202 % 568815043263
  Passed:   
  250530684965279 % 1811055503527
  Passed:   
  43497985509891 % 2566986060653
  Passed:   
  141711472568173 % 98312031422
  Passed:   
  121278120735174 % 7291914184429
  Passed:   
  (-147562462382707) % 1771917424401
  Passed:   
  (-324880460707855) % 3174757048141
  Passed:   
  (-524141830329171) % 4257616728179
  Passed:   
  402302293029904 % 6006135033215
  Passed:   
  109760026030327 % 9235795361984
  Passed:   
  299680927494929 % 2702925540266
  Passed:   
  (-543368241802641) % 2107283788940
  Passed:   
  247271165585597 % 5434442463210
  Passed:   
  (-139088214814555) % 1552244371702
  Passed:   
  599098456669089 % 385667799901
  Passed:   
  (-203684593105028) % 2274842279547
  Passed:   
  42525242468213 % 48119845245
  Passed:   
  (-105766208155669) % 1558441423963
  Passed:   
  (-412407833111131) % 4188726274137
  Passed:   
  225552331901157 % 8818333576793
  Passed:   
  (-129129125253969) % 1828301954564
  Passed:   
  281086187565442 % 6778036064417
  Passed:   
  (-566341489448993) % 7383824017887
  Passed:   
  10836500052331 % 198735103052
  Passed:   
  807851942917132 % 4690669292593
  Passed:   
  59853908586384 % 436894123951
  Passed:   
  (-373353472029451) % 2562548788743
  Passed:   
  28725866844478 % 3261694948765
  Passed:   
  (-277342703389826) % 9273579280753
  Passed:   
  (-538083608317924) % 6023228287859
  Passed:   
  (-705036515362431) % 9703276577398
  Passed:   
  (-95322573966737) % 2904547141235
  Passed:   
  77781855022925 % 2086621157657
  Passed:   
  (-199176972632101) % 474722076748
  Passed:   
  827853802907850 % 4443227503433
  Passed:   
  286584402406559 % 1003462078535
  Passed:   
  (-195595743485494) % 4643981758567
  Passed:   
  244786496629207 % 2885416954375
  Passed:   
  (-830786791620) % 23171578669
  Passed:   
  (-886324448647) % 5672970176
  Passed:   
  1249713170911 % 159111156329
  Passed:   
  (-99490220475093) % 1255980833900
  Passed:   
  (-445170783092013) % 2379282629701
  +++ OK, passed 100 tests.

