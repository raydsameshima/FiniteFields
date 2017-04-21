GUniFin.lhs

Non sequential inputs Newton-interpolation with finite fields.
Our target is a function
  f :: Q -> Q
which means to determine (canonical) coefficients.
Accessible input is pairs of in-out, i.e., a (sub) graph of f.

> module GUniFin where
> --
> import Data.Ratio
> import Data.Maybe
> import Data.Either
> import Data.List
> import Control.Monad
> --
> import Polynomials
> import Ffield hiding (takeUntil)
> --
> type Q = Ratio Int   -- Rational fields
> type Graph = [(Q,Q)] -- [(x, f x) | x <- someFinieRange]
> --
> -- f [a,b,c ..] -> [(f a b), (f b c) ..]
> -- pair wise application
> map' :: (a -> a -> b) -> [a] -> [b]
> map' f as = zipWith f as (tail as)
>
> -- To select Z_p valid inputs.
> sample :: Int   -- prime
>        -> Graph -- increasing input
>        -> Graph 
> sample p = filter ((< (fromIntegral p)) . fst)
>
> -- To eliminate (1%p) type "fake" infinity.
> -- After eliminating these, we can freely use 'modp', primed version.
> check :: Int   -- prime
>       -> Graph 
>       -> Graph -- safe data sets
> check p gs = filter (not . isDanger p) gs
>   where
>     isDanger -- To detect (1%p) type infinity.
>       :: Int -- prime 
>       -> (Q,Q) -> Bool
>     isDanger p (_, fx) = (d `rem` p) == 0
>       where 
>         d = denominator fx
> 
> project :: Int -> (Q,Q) -> (Int, Int)
> project p (x, fx) -- for simplicity
>   | denominator x == 1 = (numerator x, fx `modp'` p)
>   | otherwise          = error "project: integer input?"
>
> -- From Graph to Zp (safe) values.
> onZp 
>   :: Int                -- base prime
>   -> Graph
>   -> [(Int, Int)] -- in-out on Zp value 
> onZp p = map (project p) . check p . sample p
>
> -- using record syntax
> data PDiff 
>   = PDiff { points    :: (Int, Int) -- end points
>           , value     :: Int        -- Zp value
>           , basePrime :: Int
>           }
>   deriving (Show, Read)
>
> toPDiff 
>   :: Int        -- prime
>   -> (Int, Int) -- in and out mod p 
>   -> PDiff
> toPDiff p (x,fx) = PDiff (x,x) fx p
>
> newtonTriangleZp :: [PDiff] -> [[PDiff]]
> newtonTriangleZp fs
>   | length fs < 3 = []
>   | otherwise     = helper [sf3] (drop 3 fs)
>   where
>     sf3 = reverse . take 3 $ fs -- [[f2,f1,f0]]
>     helper fss [] = error "newtonTriangleZp: need more evaluation" 
>     helper fss (f:fs)
>       | isConsts 3 . last $ fss = fss
>       | otherwise               = helper (add1 f fss) fs
>
> isConsts 
>   :: Int -- 3times match
>   -> [PDiff] -> Bool
> isConsts n ds
>   | length ds < n = False    
> -- isConsts n ds     = all (==l) $ take (n-1) ls
>   | otherwise     = all (==l) $ take (n-1) ls
>   where 
>     (l:ls) = map value ds
>
> -- backward, each [PDiff] is decreasing inputs (i.e., reversed)
> add1 :: PDiff -> [[PDiff]] -> [[PDiff]]
> add1 f [gs] = fgs : [zipWith bdiffStep fgs gs] -- singleton
>   where 
>     fgs = f:gs
> add1 f (gg@(g:gs) : hhs) -- gg is reversed order
>             = (f:gg) : add1 fg hhs
>   where
>     fg = bdiffStep f g
>
> -- backward
> bdiffStep :: PDiff -> PDiff -> PDiff
> bdiffStep (PDiff (y,y') g q) (PDiff (x,x') f p)
>   | p == q    = PDiff (x,y') finiteDiff p
>   | otherwise = error "bdiffStep: different primes?"
>   where
>     finiteDiff = ((fg % xy') `modp'` p)
>     xy' = (x - y' `mod` p)
>     fg = ((f-g) `mod` p) 
>
> graph2Zp :: Int -> Graph -> [(Int, Int)]
> graph2Zp p = onZp p . check p . sample p 
>
> graph2PDiff :: Int -> Graph -> [PDiff]
> graph2PDiff p = map (toPDiff p) . graph2Zp p
>
> newtonTriangleZp' :: Int -> Graph -> [[PDiff]]
> newtonTriangleZp' p = newtonTriangleZp . graph2PDiff p
> 
> newtonCoeffZp :: Int -> Graph -> [PDiff]
> newtonCoeffZp p = map head . newtonTriangleZp' p

  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) 
                         [1,2,4,5,9,10,11] :: Graph 
  *GUniFin> newtonCoeffZp 101 gs
  [PDiff {points = (9,9), value = 69, basePrime = 101}
  ,PDiff {points = (5,9), value = 65, basePrime = 101}
  ,PDiff {points = (4,9), value = 1, basePrime = 101}
  ]
  *GUniFin> map (\x -> (Just . value $ x, basePrime x)) it
  [(Just 69,101),(Just 65,101),(Just 1,101)]

We take formally the canonical form on Zp, 
then apply rational "number" reconstruction.

> n2cZp :: [PDiff] -> ([Int], Int)
> n2cZp graph = (helper graph, p)
>   where
>     p = basePrime . head $ graph
>     helper [d]    = [value d]
>     helper (d:ds) = map (`mod` p) $ ([value d] + (z * next)) 
>                                    - (map (`mod` p) (zd .* next))
>       where 
>         zd = fst . points $ d
>         next = helper ds
>
> format :: ([Int],Int) -> [(Maybe Int, Int)]
> format (as,p) = [(return a,p) | a <- as]
  
  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) 
                         [0,2,3,5,7,8,11] :: Graph 
  *GUniFin> newtonCoeffZp 10007 gs
  [PDiff {points = (7,7), value = 8392, basePrime = 10007}
  ,PDiff {points = (5,7), value = 5016, basePrime = 10007}
  ,PDiff {points = (3,7), value = 1, basePrime = 10007}
  ]
  *GUniFin> n2cZp it
  ([3336,5004,1],10007)
  *GUniFin> format it
  [(Just 3336,10007),(Just 5004,10007),(Just 1,10007)]
  *GUniFin> map guess it
  [Just (1 % 3,10007),Just (1 % 2,10007),Just (1 % 1,10007)]

  *GUniFin> let gs = map (\x -> (x,x^2 + (1%2)*x + 1%3)) 
                         [0,2,3,5,7,8,11] :: Graph 
  *GUniFin> map guess . format . n2cZp . newtonCoeffZp 10007 $ gs
  [Just (1 % 3,10007),Just (1 % 2,10007),Just (1 % 1,10007)]
  *GUniFin> let gs = map (\x -> (x,x^5 + x^2 + (1%2)*x + 1%3)) 
                         [0,2,3,5,7,8,11,13,17,18,19,21,24,28,31,33,34] 
                         :: Graph 
  *GUniFin> map guess . format . n2cZp . newtonCoeffZp 10007 $ gs
  [Just (1 % 3,10007),Just (1 % 2,10007),Just (1 % 1,10007)
  ,Just (0 % 1,10007),Just (0 % 1,10007),Just (1 % 1,10007)
  ] 

> preTrial gs p = format . n2cZp . newtonCoeffZp p $ gs

  *GUniFin> let gs = map (\x -> (x,x^5 + x^2 + (1%2)*x + 1%3)) 
                         [0,2,3,5,7,8,11,13,17,18,19,21,24,28,31,33,34] 
                         :: Graph 
  *GUniFin> map reconstruct . transpose . map (preTrial gs) $ bigPrimes 
  [Just (1 % 3),Just (1 % 2),Just (1 % 1)
  ,Just (0 % 1),Just (0 % 1),Just (1 % 1)
  ]

Here is "a" final version, the univariate polynomial reconstruction 
with finite fields.

> uniPolCoeff :: Graph -> Maybe [(Ratio Int)]
> uniPolCoeff gs = sequence . map reconstruct' . transpose . 
>                  map (preTrial gs) $ bigPrimes

  *GUniFin> let gs = map (\x -> (x,x^5 + x^2 + (1%2)*x + 1%3)) 
                         [0,2,3,5,7,8,11,13,17,18,19,21,24,28,31,33,34]
                         :: Graph
  *GUniFin> gs
  [(0 % 1,1 % 3),(2 % 1,112 % 3),(3 % 1,1523 % 6),(5 % 1,18917 % 6)
  ,(7 % 1,101159 % 6),(8 % 1,98509 % 3),(11 % 1,967067 % 6)
  ,(13 % 1,2228813 % 6),(17 % 1,8520929 % 6),(18 % 1,5669704 % 3)
  ,(19 % 1,14858819 % 6),(21 % 1,24507317 % 6),(24 % 1,23889637 % 3)
  ,(28 % 1,51633499 % 3),(31 % 1,171780767 % 6),(33 % 1,234818993 % 6)
  ,(34 % 1,136309792 % 3)
  ]
  *GUniFin> uniPolCoeff gs
  Just [1 % 3,1 % 2,1 % 1,0 % 1,0 % 1,1 % 1]

  *GUniFin> let fs = map (\x -> (x,(3+x+(1%3)*x^9)/(1))) 
                         [1,3..101] :: Graph
  *GUniFin> uniPolCoeff fs
  Just [3 % 1,1 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,0 % 1,1 % 3]
  *GUniFin> let fs = map (\x -> (x,(3+x+(1%3)*x^10)/(1))) 
                         [1,3..101] :: Graph
  *GUniFin> uniPolCoeff fs
  *** Exception: newtonBT: need more evaluation
  CallStack (from HasCallStack):
    error, called at GUniFin.lhs:79:23 in main:GUniFin
  *GUniFin> let fs = map (\x -> (x,(3+x+(1%3)*x^10)/(1))) 
                         [1,3..1001] :: Graph
  *GUniFin> uniPolCoeff fs
  *** Exception: newtonBT: need more evaluation
  CallStack (from HasCallStack):
    error, called at GUniFin.lhs:79:23 in main:GUniFin

Rough estimation says, in 64-bits system with sequential inputs,
the upper limit of degree is about 15.
If we use non sequential inputs, this upper limit will go down.

-- up to here, polinomials
--
-- from now on, rational functions

Non sequential inputs Thiele-interpolation with finite fields.

Let me start naive rho with non-sequential inputs:

> -- over rational (infinite) field
> rho :: Graph -> Int -> [Q]
> rho gs 0 = map snd gs
> rho gs 1 = zipWith (/) xs' fs'
>   where
>     xs' = zipWith (-) xs (tail xs)
>     xs  = map fst gs
>     fs' = zipWith (-) fs (tail fs)
>     fs  = map snd gs
> rho gs n = zipWith (+) twoAbove oneAbove
>   where
>     twoAbove = zipWith (/) xs' rs'
>     xs' = zipWith (-) xs (drop n xs)
>     xs  = map fst gs
>     rs' = zipWith (-) rs (tail rs)
>     rs  = rho gs (n-1)
>     oneAbove = tail $ rho gs (n-2)

This works 

  *GUniFin> let func x = (1+x+2*x^2)/(3+2*x +(1%4)*x^2)
  *GUniFin> let fs = map (\x -> (x, func x)) 
                         [0,1,3,4,6,7,9,10,11,13,14,15,17,19,20] :: Graph 
  *GUniFin> let r = rho fs
  *GUniFin> r 0
  [1 % 3,16 % 21,88 % 45,37 % 15,79 % 24,424 % 117,688 % 165,211 % 48
  ,1016 % 221,1408 % 285,407 % 80,1864 % 357,2384 % 437,424 % 75,821 % 143
  ]
  *GUniFin> r 1
  [7 % 3,315 % 188,45 % 23,80 % 33,936 % 311,6435 % 1756,880 % 199
  ,10608 % 2137,62985 % 10804,4560 % 671,28560 % 3821,156009 % 18260
  ,32775 % 3244,10725 % 943
  ]
  *GUniFin> r 2
  [(-604) % 159,5116 % 405,9458 % 1065,18962 % 2253,75244 % 9171
  ,117388 % 14439,174700 % 21603,243084 % 30151,329516 % 40955
  ,436876 % 54375,559148 % 69659,26491 % 3303,138404 % 17267
  ]
  *GUniFin> r 3
  [900 % 469,585 % 938,(-5805) % 938,(-19323) % 938,(-23418) % 469
  ,(-165867) % 1876,(-295485) % 1876,(-111560) % 469,(-651015) % 1876
  ,(-977265) % 1876,(-199317) % 268,(-278589) % 268
  ]
  *GUniFin> r 4
  [8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1,8 % 1]

But here is a corner case, an accidental match.
We should detect them and handle them safely.

  *GUniFin> let func x = (x%(1+x^2))
  *GUniFin> let fs = map (\x -> (x%1,func x)) [0..10]
  *GUniFin> let r = rho fs
  *GUniFin> r 0
  [0 % 1,1 % 2,2 % 5,3 % 10,4 % 17,5 % 26
  ,6 % 37,7 % 50,8 % 65,9 % 82,10 % 101
  ]
  *GUniFin> r 1
  [2 % 1,(-10) % 1,(-10) % 1,(-170) % 11,(-442) % 19
  ,(-962) % 29,(-1850) % 41,(-650) % 11,(-5330) % 71,(-8282) % 89
  ]
  *GUniFin> r 2
  [1 % 3,*** Exception: Ratio has zero denominator

--

> -- We have assumed our out-puts are safe, i.e. no fake infinity.
> initialThieleZp :: [PDiff] -> [[PDiff]]
> initialThieleZp fs 
>   | isConsts 3 fs = [first]   
>   | otherwise = [first, second]
>   where
>     first  = reverse . take 4 $ fs
>     second = map' reciproDiff first

> {-
> -- To make safe initials.
> initialThieleTriangleZp :: [PDiff] -> [[PDiff]]
> initialThieleTriangleZp ff@(f:fs)
>   | isConsts 3 ff = [reverse $ take 3 ff]
>   | otherwise     = [[g,f],[h]]
>   where
>     g = firstDifferent f fs
>
>     firstDifferent _ [] 
>       = error "initialThieleTriangleZp: need more points"
>     firstDifferent f (g:gs)
>       = if (g' /= f') then g 
>                       else firstDifferent f gs
>       where
>         f' = value f
>         g' = value g
>     
>     h = reciproDiff g f
>
> -- to make safe first two stairs.
> initialThieleZp' :: [PDiff] -> [[PDiff]]
> initialThieleZp' fs 
>   | isConsts 3 fs = [reverse $ take 3 fs]
>   | otherwise     = [firsts, seconds]
>   where
>     firsts = undefined
>     seconds = undefined
> -}

> -- reversed order
> reciproDiff :: PDiff -> PDiff -> PDiff
> reciproDiff (PDiff (_,z') u p) (PDiff (w,_) v q)
>   | p /= q = error "reciproDiff: wrong base prime"
>   | otherwise = PDiff (w,z') r p
>   where
>     r = ((zw) * (uv `inversep'` p)) `mod` p
>     zw = (z' - w) `mod` p
>     uv = (u  - v) `mod` p -- assuming (u-v) is not "zero"
>
> -- reciproAdd1 :: PDiff -> ListZipper [PDiff] -> ListZipper [PDiff]
> -- reciproAdd1 _ ([], sb) = ([], sb) -- reversed order
> -- reciproAdd1 f ((gs@(g:_) : hs : iss), []) 
> --                        = reciproAdd1 k (iss, [(j:hs), (f:gs)])
> --   where
> --     j = reciproDiff f g
> --     k = addZp' j g
>
> addZp' :: PDiff -> PDiff -> PDiff
> addZp' (PDiff (x,y) v p) (PDiff (_,_) w q)
>   | p /= q = error "addZp': wrong primes"
>   | otherwise = PDiff (x,y) vw p   
>   where 
>     vw = (v + w) `mod` p
>
> -- This takes new point and the heads, and returns the new heads.
> thieleHeads 
>   :: PDiff   -- a new element  rho8
>   -> [PDiff] -- oldies         [rho7, rho67, rho567, rho4567 ..]
>   -> [PDiff] --                [rho8, rho78, rho678, rho5678 ..]
> thieleHeads _ []        = []
> thieleHeads f gg@(g:gs) = f : fg : helper fg gg
>   where
>     fg  = reciproDiff f g
>
>     helper _ bs
>       | length bs < 3 = []
>     helper a (b:bs@(c:d:_)) = e : helper e bs
>       where
>         e = addZp' (reciproDiff a c) 
>                    b
>
>     tHs :: PDiff -> [PDiff] -> [PDiff] -- reciprocal diff. part
>     tHs _ [] = []
>     tHs f' hh@(h:hs) = fh : tHs fh hs
>       where 
>         fh = reciproDiff f' h
>
> thieleTriangle' :: [PDiff] -> [[PDiff]]
> thieleTriangle' fs 
>   | length fs < 4 = []
>   | otherwise     = helper fourThree (drop 4 fs)
>   where
>     fourThree = initialThieleZp fs
>     helper fss [] 
>       | isConsts 3 . last $ fss = fss
>       | otherwise               = error "thieleTriangle: need more inputs"
>     helper fss (g:gs) 
>       | isConsts 3 . last $ fss = fss
>       | otherwise               = helper gfss gs
>       where
>         gfss = thieleComp g fss
>
> thieleComp :: PDiff -> [[PDiff]] -> [[PDiff]]
> thieleComp g fss = wholeButLast ++ [three]
>   where
>     wholeButLast = zipWith (:) hs fss
>     hs = thieleHeads g (map head fss)
>     three = fiveFour2three $ last2 wholeButLast
>     -- Finally from two stairs (5 and 4 elements),
>     -- we create the bottom 3 elements.
>
>     last2 :: [a] -> [a]
>     last2 [a,b] = [a,b]
>     last2 (_:bb@(_:_)) = last2 bb
>
> fiveFour2three -- This works!
>   :: [[PDiff]] -- 5 and 4, under last2
>   -> [PDiff]   -- 3
> fiveFour2three [ff@(_:fs), gg] = zipWith addZp' (map' reciproDiff gg) fs
>
> thieleTriangle :: Graph -> Int -> [[PDiff]]
> thieleTriangle fs p = thieleTriangle' $ graph2PDiff p fs
>
> thieleCoeff' :: Graph -> Int -> [PDiff]
> thieleCoeff' fs = map last . thieleTriangle fs 
>
> -- thieleCoeff'' fs p = a:b:(zipWith subZp bs as)
> thieleCoeff'' fs p 
>   | length (thieleCoeff' fs p) == 1 = thieleCoeff' fs p
>   | otherwise = a:b:(zipWith subZp bs as)
>   where
>     as@(a:b:bs) = thieleCoeff' fs p
> --  as@(a:b:bs) = firstReciprocalDifferences fs p
>
>     subZp :: PDiff -> PDiff -> PDiff
>     subZp (PDiff (x,y) v p) (PDiff (_,_) w q)
>       | p /= q    = error "thileCoeff: different primes"
>       | otherwise = PDiff (x,y) ((v-w) `mod` p) p
>

> t2cZp 
>   :: [PDiff]              -- thieleCoeff'' fs p
>   -> (([Int],[Int]), Int) -- (rat-func, basePrime)
> t2cZp gs = (helper gs, p)
>   where
>     p = basePrime . head $ gs
>     helper [n]   = ([(value n) `mod` p], [1])
>     helper [d,e] = ([de',1], [e']) -- base case
>       where 
>         de' = ((d'*e' `mod` p) - xd) `mod` p
>         d'  = value d
>         e'  = value e
>         xd  = snd . points $ d
>     helper (d:ds) = (den', num)
>       where
>         (num, den) = helper ds
>         den'  = map (`mod` p) $ (z * den) + (map (`mod` p) $ num''-den'')
>         num'' = map (`mod` p) ((value d) .* num)
>         den'' = map (`mod` p) ((snd . points $ d) .* den)
>
> -- pre "canonicalizer"
> beforeFormat' :: (([Int], [Int]), Int) -> (([Int], [Int]), Int)
> beforeFormat' ((num,(d:ds)), p) = ((num',den'), p)
>   where
>     num' = map (`mod` p) $ di .* num
>     den' = 1: (map (`mod` p) $ di .* ds)
>     di   = d `inversep'` p
>
> format'
>   :: (([Int], [Int]), Int)
>   -> ([(Maybe Int, Int)], [(Maybe Int, Int)])
> format' ((num,den), p) = (format (num, p), format (den, p))
  
  *GUniFin> let fs = map (\x -> (x,(1+x)/(2+x))) [0,2,3,4,6,8,9] :: Graph 
  *GUniFin> thieleCoeff'' fs 101
  [PDiff {points = (0,0), value = 51, basePrime = 101}
  ,PDiff {points = (0,2), value = 8, basePrime = 101}
  ,PDiff {points = (0,3), value = 51, basePrime = 101}
  ]
  *GUniFin> t2cZp it
  (([1,1],[2,1]),101)
  *GUniFin> format' it
  ([(Just 1,101),(Just 1,101)]
  ,[(Just 2,101),(Just 1,101)]
  )
  *GUniFin> format' . t2cZp . thieleCoeff'' fs $ 101
  ([(Just 1,101),(Just 1,101)],[(Just 2,101),(Just 1,101)])
  *GUniFin> format' . t2cZp . thieleCoeff'' fs $ 103
  ([(Just 1,103),(Just 1,103)],[(Just 2,103),(Just 1,103)])
  *GUniFin> format' . t2cZp . thieleCoeff'' fs $ 107
  ([(Just 1,107),(Just 1,107)],[(Just 2,107),(Just 1,107)])

> ratCanZp
>   :: Graph -> Int -> ([(Maybe Int, Int)], [(Maybe Int, Int)])
> ratCanZp fs = format' . beforeFormat' . t2cZp . thieleCoeff'' fs 

  *GUniFin> let fivePrimes = take 5 bigPrimes 
  *GUniFin> let fs = map (\x -> (x,(1+x)/(2+x))) [0,2,3,4,6,8,9] :: Graph 
  *GUniFin> map (ratCanZp fs) fivePrimes 
  [([(Just 5004,10007),(Just 5004,10007)]
   ,[(Just 1,10007),(Just 5004,10007)]
   )
  ,([(Just 5005,10009),(Just 5005,10009)]
   ,[(Just 1,10009),(Just 5005,10009)]
   )
  ,([(Just 5019,10037),(Just 5019,10037)]
   ,[(Just 1,10037),(Just 5019,10037)]
   )
  ,([(Just 5020,10039),(Just 5020,10039)]
   ,[(Just 1,10039),(Just 5020,10039)]
   )
  ,([(Just 5031,10061),(Just 5031,10061)]
   ,[(Just 1,10061),(Just 5031,10061)]
   )
  ]
  *GUniFin> map fst it
  [[(Just 5004,10007),(Just 5004,10007)]
  ,[(Just 5005,10009),(Just 5005,10009)]
  ,[(Just 5019,10037),(Just 5019,10037)]
  ,[(Just 5020,10039),(Just 5020,10039)]
  ,[(Just 5031,10061),(Just 5031,10061)]
  ]
  *GUniFin> transpose it
  [[(Just 5004,10007),(Just 5005,10009),(Just 5019,10037)
   ,(Just 5020,10039),(Just 5031,10061)
   ]
  ,[(Just 5004,10007),(Just 5005,10009),(Just 5019,10037)
   ,(Just 5020,10039),(Just 5031,10061)
   ]
  ]
  *GUniFin> map reconstruct it
  [Just (1 % 2),Just (1 % 2)]
  *GUniFin> map (ratCanZp fs) fivePrimes 
  [([(Just 5004,10007),(Just 5004,10007)]
   ,[(Just 1,10007),(Just 5004,10007)]
   )
  ,([(Just 5005,10009),(Just 5005,10009)]
   ,[(Just 1,10009),(Just 5005,10009)]
   )
  ,([(Just 5019,10037),(Just 5019,10037)]
   ,[(Just 1,10037),(Just 5019,10037)]
   )
  ,([(Just 5020,10039),(Just 5020,10039)]
   ,[(Just 1,10039),(Just 5020,10039)]
   )
  ,([(Just 5031,10061),(Just 5031,10061)]
   ,[(Just 1,10061),(Just 5031,10061)]
   )
  ]
  *GUniFin> map snd it
  [[(Just 1,10007),(Just 5004,10007)]
  ,[(Just 1,10009),(Just 5005,10009)]
  ,[(Just 1,10037),(Just 5019,10037)]
  ,[(Just 1,10039),(Just 5020,10039)]
  ,[(Just 1,10061),(Just 5031,10061)]
  ]
  *GUniFin> transpose it
  [[(Just 1,10007),(Just 1,10009),(Just 1,10037)
   ,(Just 1,10039),(Just 1,10061)
   ]
  ,[(Just 5004,10007),(Just 5005,10009),(Just 5019,10037)
   ,(Just 5020,10039),(Just 5031,10061)
   ]
  ]
  *GUniFin> map reconstruct it
  [Just (1 % 1),Just (1 % 2)]
  
> -- uniPolCoeff :: Graph -> Maybe [(Ratio Integer)]
> -- uniPolCoeff gs = sequence . map reconstruct . transpose . map (preTrial gs) $ bigPrimes
  
> -- Clearly this is double running implementation.
> uniRatCoeff
>   :: Graph -> ([Maybe (Ratio Int)], [Maybe (Ratio Int)])
> uniRatCoeff gs = (num, den)
>   where
>     (num,den) = (helper fst, helper snd)
>     helper third 
>       = map reconstruct' . transpose 
>         . map (third . ratCanZp gs) $ bigPrimes

  *GUniFin> let fs = map (\x -> (x,(1+2*x+x^10)/(1+(3%2)*x+x^5))) [0..101] :: Graph
  (0.01 secs, 44,232 bytes)
  *GUniFin> uniRatCoeff fs
  ([Just (1 % 1),Just (2 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ,Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ,Just (1 % 1)
   ]
  ,[Just (1 % 1),Just (3 % 2),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ,Just (1 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ]
  )
  (1.72 secs, 1,424,003,616 bytes)

> isJustZero n = Just (0%1) == n
> 
> uniRatCoeffShort gs = (num', den')
>   where
>     (num, den) = uniRatCoeff gs
>     (num', den') = (helper num, helper den)
>     helper nd = filter (not . isJustZero . fst) $ zip nd [0..]
  
  *GUniFin> let fs = map (\x -> (x,(1+2*x+x^10)/(1+(3%2)*x+x^5))) [0..101] :: Graph
  (0.01 secs, 44,320 bytes)
  *GUniFin> uniRatCoeff fs
  ([Just (1 % 1),Just (2 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ,Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ,Just (1 % 1)
   ]
  ,[Just (1 % 1),Just (3 % 2),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ,Just (1 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1),Just (0 % 1)
   ]
  )
  (1.72 secs, 1,424,009,472 bytes)
  *GUniFin> uniRatCoeffShort  fs
  ([(Just (1 % 1),0),(Just (2 % 1),1),(Just (1 % 1),10)]
  ,[(Just (1 % 1),0),(Just (3 % 2),1),(Just (1 % 1),5)]
  )
  (1.74 secs, 1,422,577,184 bytes)

> uniRatCoeff'
>   :: Graph -> (Maybe [Ratio Int], Maybe [Ratio Int])
> uniRatCoeff' gs = (num', den')
>   where
>     (num, den) = uniRatCoeff gs
>     num' = sequence num
>     den' = sequence den

> func2graph :: (Q -> Q) -> [Q] -> Graph
> func2graph f xs = [(x, f x) | x <- xs]

  *GUniFin> func2graph g [0,3..30]
  [(0 % 1,0 % 1),(3 % 1,27 % 64),(6 % 1,216 % 343),(9 % 1,729 % 1000)
  ,(12 % 1,1728 % 2197),(15 % 1,3375 % 4096),(18 % 1,5832 % 6859)
  ,(21 % 1,9261 % 10648),(24 % 1,13824 % 15625),(27 % 1,19683 % 21952)
  ,(30 % 1,27000 % 29791)
  ]
  (0.01 secs, 363,080 bytes)
  *GUniFin> uniRatCoeffShort it
  ([(Just (1 % 1),3)]
  ,[(Just (1 % 1),0),(Just (3 % 1),1),(Just (3 % 1),2),(Just (1 % 1),3)]
  )
  (0.30 secs, 231,980,488 bytes) 

> -- Up to degree~100 version.
> ratFunc2Coeff
>   :: (Q -> Q) -- rational function
>   -> (Maybe [Ratio Int], Maybe [Ratio Int])
> ratFunc2Coeff f = uniRatCoeff' . func2graph f $ [0..100]

--

We want to use safe list, i.e., the given graph as much as possible.
So, the easiest way could be
  pick a prime
  construct Thiele triangle up to consts.
  if we face fake infinity before it matches, 
    then return Nothing and use another prime
Since we have a lot of bigPrimes.
  
  *GUniFin> let g x = x^4 / (1+x)^3
  *GUniFin> func2graph g [0..10]
  [(0 % 1,0 % 1),(1 % 1,1 % 8),(2 % 1,16 % 27),(3 % 1,81 % 64)
  ,(4 % 1,256 % 125),(5 % 1,625 % 216),(6 % 1,1296 % 343),(7 % 1,2401 % 512)
  ,(8 % 1,4096 % 729),(9 % 1,6561 % 1000),(10 % 1,10000 % 1331)
  ]
  *GUniFin> let p = head bigPrimes 
  *GUniFin> p
  10007
  *GUniFin> graph2PDiff p $ func2graph g [0..10]
  [PDiff {points = (0,0), value = 0, basePrime = 10007}
  ,PDiff {points = (1,1), value = 1251, basePrime = 10007}
  ,PDiff {points = (2,2), value = 2595, basePrime = 10007}
  ,PDiff {points = (3,3), value = 6412, basePrime = 10007}
  ,PDiff {points = (4,4), value = 1363, basePrime = 10007}
  ,PDiff {points = (5,5), value = 2273, basePrime = 10007}
  ,PDiff {points = (6,6), value = 1375, basePrime = 10007}
  ,PDiff {points = (7,7), value = 8624, basePrime = 10007}
  ,PDiff {points = (8,8), value = 9038, basePrime = 10007}
  ,PDiff {points = (9,9), value = 7782, basePrime = 10007}
  ,PDiff {points = (10,10), value = 7150, basePrime = 10007}
  ]

We need Maybe-wrapped version of reciprocal (inverse) difference.

> -- normal order for rho-matrix
> inverseDiff 
>   :: PDiff -> PDiff -> Maybe PDiff
> inverseDiff (PDiff (w,_) v p) (PDiff (_,z') u _) -- z in reserved
>   | v == u    = Nothing
>   | otherwise = return $ PDiff (w,z') r p
>   where
>     r = ((zw) * (uv `inversep'` p)) `mod` p
>     zw = (z' - w) `mod` p
>     uv = (u - v) `mod` p  
>
> inverseDiff' :: Maybe PDiff -> Maybe PDiff -> Maybe PDiff
> inverseDiff' Nothing  _        = Nothing
> inverseDiff' _        Nothing  = Nothing
> inverseDiff' (Just a) (Just b) = inverseDiff a b

> -- rho-matrix version
> -- This implementation is quite straightforward, but no error handling.
> inverseDiffs 
>   :: Int     -- a prime
>   -> Int     -- "degree" or the depth of thiele TRAPEZOID
>   -> Graph
>   -> [Maybe PDiff]
> inverseDiffs p 0 fs = map return $ graph2PDiff p fs 
> inverseDiffs p 1 fs = map' inverseDiff $ graph2PDiff p fs
> inverseDiffs p n fs 
>   = zipWith addPDiff (map' inverseDiff' (inverseDiffs p (n-1) fs))
>                      (tail $ inverseDiffs p (n-2) fs)
>
> addPDiff :: Maybe PDiff -> Maybe PDiff -> Maybe PDiff
> addPDiff Nothing _ = Nothing
> addPDiff _ Nothing = Nothing
> addPDiff (Just a) (Just b) = return $ addZp' a b
  
  *GUniFin> let f x = x / (1+x^2)
  *GUniFin> let fs = func2graph f [0..10]
  *GUniFin> sequence $ filter isJust $ inverseDiffs 10007 4 fs
  Just [PDiff {points = (2,6), value = 0, basePrime = 10007}
       ,PDiff {points = (3,7), value = 0, basePrime = 10007}
       ,PDiff {points = (4,8), value = 0, basePrime = 10007}
       ,PDiff {points = (5,9), value = 0, basePrime = 10007}
       ,PDiff {points = (6,10), value = 0, basePrime = 10007}
       ]
  *GUniFin> fmap (isConsts 3) it
  Just True
    
> isConsts' 
>   :: Int -> [Maybe PDiff] -> Bool
> isConsts' n fs
>   | Just True == fmap (isConsts n) fs' = True
>   | otherwise                          = False
>   where
>     fs' = sequence . filter isJust $ fs
>
> -- This is the main function which returns Nothing when we face 
> -- so many fake infinities with really bad prime.
> thieleTrapezoid 
>   :: Graph -> Int -> Maybe [[Maybe PDiff]]
> thieleTrapezoid fs p
>   | any (isConsts' 3) gs = return gs'
> --  | or $ map (isConsts' 3) gs = return gs' 
>   | otherwise                 = Nothing
>   where
>     gs' = aMatrix fs p
>     gs  = map (filter isJust) gs'
>
>     aMatrix 
>       :: Graph -> Int -> [[Maybe PDiff]]
>     aMatrix fs p = takeUntil (isConsts' 3) 
>                      [inverseDiffs p n fs | n <- [0..]]
>
> takeUntil 
>   :: (a -> Bool) -> [a] -> [a]
> takeUntil _ []     = []
> takeUntil f (x:xs) 
>   | not (f x) = x : takeUntil f xs
>   | f x       = [x]

Finally, we need the Thiele coefficients!

  *GUniFin> fmap head . join . fmap (sequence . filter isJust . map sequence . transpose) . thieleTrapezoid  fs $ 10007
  Just [PDiff {points = (2,2), value = 8006, basePrime = 10007},PDiff {points = (2,3), value = 9997, basePrime = 10007},PDiff {points = (2,4), value = 5337, basePrime = 10007},PDiff {points = (2,5), value = 50, basePrime = 10007},PDiff {points = (2,6), value = 0, basePrime = 10007}]

> thieleCoefficients 
>   :: Graph -> Int -> Maybe [PDiff]
> thieleCoefficients fs
>   = fmap head . join 
>     . fmap (sequence . filter isJust . map sequence . transpose) 
>     . thieleTrapezoid fs 

  *GUniFin> let f x = x / (1+x^2)
  *GUniFin> let fs = func2graph f [0..10]
  *GUniFin> :t thieleCoefficients fs 10007
  thieleCoefficients fs 10007 :: Maybe [PDiff]
  *GUniFin> thieleCoefficients fs 10007
  Just [PDiff {points = (2,2), value = 8006, basePrime = 10007}
       ,PDiff {points = (2,3), value = 9997, basePrime = 10007}
       ,PDiff {points = (2,4), value = 5337, basePrime = 10007}
       ,PDiff {points = (2,5), value = 50, basePrime = 10007}
       ,PDiff {points = (2,6), value = 0, basePrime = 10007}
       ]

> thieleCoefficients' Nothing   = Nothing
> thieleCoefficients' (Just cs) = return (a:b:zipWith subZp bs as)
>   where
>     as@(a:b:bs) = cs
>
>     subZp :: PDiff -> PDiff -> PDiff
>     subZp (PDiff (x,y) v p) (PDiff (_,_) w _)
>       = PDiff (x,y) ((v-w) `mod` p) p 

  *GUniFin> let f x = x / (1+x^2)
  *GUniFin> let fs = func2graph f [0..10]
  *GUniFin> thieleCoefficients fs 10007
  Just [PDiff {points = (2,2), value = 8006, basePrime = 10007}
       ,PDiff {points = (2,3), value = 9997, basePrime = 10007}
       ,PDiff {points = (2,4), value = 5337, basePrime = 10007}
       ,PDiff {points = (2,5), value = 50, basePrime = 10007}
       ,PDiff {points = (2,6), value = 0, basePrime = 10007}
       ]
  *GUniFin> thieleCoefficients' it
  Just [PDiff {points = (2,2), value = 8006, basePrime = 10007}
       ,PDiff {points = (2,3), value = 9997, basePrime = 10007}
       ,PDiff {points = (2,4), value = 7338, basePrime = 10007}
       ,PDiff {points = (2,5), value = 60, basePrime = 10007}
       ,PDiff {points = (2,6), value = 4670, basePrime = 10007}
       ]
  *GUniFin> fmap t2cZp it
  Just (([0,1,0],[1,0,1]),10007)
  *GUniFin> fmap format' it
  Just ([(Just 0,10007),(Just 1,10007),(Just 0,10007)]
       ,[(Just 1,10007),(Just 0,10007),(Just 1,10007)]
       )
  
  *GUniFin> fmap (format' . t2cZp) . thieleCoefficients' . thieleCoefficients fs $ 10007
  Just ([(Just 0,10007),(Just 1,10007),(Just 0,10007)],[(Just 1,10007),(Just 0,10007),(Just 1,10007)])

> ratCanZp'
>   :: Graph -> Int -> Maybe ([(Maybe Int, Int)], [(Maybe Int, Int)])
> -- ratCanZp' fs = fmap (format' . t2cZp) . thieleCoefficients' 
> ratCanZp' fs 
>   = fmap (format' . beforeFormat' . t2cZp) . thieleCoefficients' 
>     . thieleCoefficients fs
  
  *GUniFin> let fivePrimes = take 5 bigPrimes 
  *GUniFin> let f x = x / (1+x^2)
  *GUniFin> let fs = func2graph f [0..10]
  *GUniFin> map (ratCanZp' fs) five
  fiveFour2three  fivePrimes
  *GUniFin> map (ratCanZp' fs) fivePrimes 
  [Just ([(Just 0,10007),(Just 1,10007),(Just 0,10007)],[(Just 1,10007),(Just 0,10007),(Just 1,10007)]),Just ([(Just 0,10009),(Just 1,10009),(Just 0,10009)],[(Just 1,10009),(Just 0,10009),(Just 1,10009)]),Just ([(Just 0,10037),(Just 1,10037),(Just 0,10037)],[(Just 1,10037),(Just 0,10037),(Just 1,10037)]),Just ([(Just 0,10039),(Just 1,10039),(Just 0,10039)],[(Just 1,10039),(Just 0,10039),(Just 1,10039)]),Just ([(Just 0,10061),(Just 1,10061),(Just 0,10061)],[(Just 1,10061),(Just 0,10061),(Just 1,10061)])]
  *GUniFin> sequence it
  Just [([(Just 0,10007),(Just 1,10007),(Just 0,10007)],[(Just 1,10007),(Just 0,10007),(Just 1,10007)]),([(Just 0,10009),(Just 1,10009),(Just 0,10009)],[(Just 1,10009),(Just 0,10009),(Just 1,10009)]),([(Just 0,10037),(Just 1,10037),(Just 0,10037)],[(Just 1,10037),(Just 0,10037),(Just 1,10037)]),([(Just 0,10039),(Just 1,10039),(Just 0,10039)],[(Just 1,10039),(Just 0,10039),(Just 1,10039)]),([(Just 0,10061),(Just 1,10061),(Just 0,10061)],[(Just 1,10061),(Just 0,10061),(Just 1,10061)])]
  *GUniFin> fmap (map fst) it
  Just [[(Just 0,10007),(Just 1,10007),(Just 0,10007)],[(Just 0,10009),(Just 1,10009),(Just 0,10009)],[(Just 0,10037),(Just 1,10037),(Just 0,10037)],[(Just 0,10039),(Just 1,10039),(Just 0,10039)],[(Just 0,10061),(Just 1,10061),(Just 0,10061)]]
  *GUniFin> fmap transpose it
  Just [[(Just 0,10007),(Just 0,10009),(Just 0,10037),(Just 0,10039),(Just 0,10061)],[(Just 1,10007),(Just 1,10009),(Just 1,10037),(Just 1,10039),(Just 1,10061)],[(Just 0,10007),(Just 0,10009),(Just 0,10037),(Just 0,10039),(Just 0,10061)]]
  *GUniFin> fmap (map reconstruct) it
  Just [Just (0 % 1),Just (1 % 1),Just (0 % 1)]
  
  *GUniFin> fmap (map reconstruct . transpose . map fst) . sequence . map (ratCanZp' fs) $ fivePrimes 
  Just [Just (0 % 1),Just (1 % 1),Just (0 % 1)]

> -- need "less data pts" error handling
> uniRatCoeffm
>   :: Graph -> (Maybe [Ratio Integer], Maybe [Ratio Integer]) 
> uniRatCoeffm fs = (num, den)
>   where
>     num = helper fst
>     den = helper snd
>     helper third 
>       = join . fmap (sequence . map reconstruct . transpose . map third) 
>         . sequence . map (ratCanZp' fs) $ bigPrimes
  
  *GUniFin> let f x = x^3 / (1+x)^4
  (0.01 secs, 48,440 bytes)
  *GUniFin> let fs = func2graph f [0..20]
  (0.02 secs, 48,472 bytes)
  *GUniFin> uniRatCoeffm fs
  (Just [0 % 1,0 % 1,0 % 1,1 % 1,0 % 1]
  ,Just [1 % 1,4 % 1,6 % 1,4 % 1,1 % 1]
  )
  (10.98 secs, 8,836,586,776 bytes)

