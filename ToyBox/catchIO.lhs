catchIO.lhs

This is IO-wrapping.

> import System.IO
> import Control.Exception 
>          ( ArithException (..)
>          , catch
>          )
> import Data.Ratio

> func x = 1%x
>
> saferEvalIO 
>   :: (Show t1, Show t) =>
>      (t -> t1) -> t -> IO ()
> saferEvalIO f x = do
>   print $ (f x) `seq` (x, f x) 
>   `catch` \RatioZeroDenominator 
>             -> return ()
> --             -> print $ "zero division"

  *Main> mapM_ (saferEvalIO func) [0..5]
  (1,1 % 1)
  (2,1 % 2)
  (3,1 % 3)
  (4,1 % 4)
  (5,1 % 5)

makeGraph_
  :: (Foldable t, Show t2, Show t1) => (t2 -> t1) -> t t2 -> IO ()

> makeGraph_ :: (Ratio Int -> Ratio Int) -> [Ratio Int] -> IO ()
> makeGraph_ f xs = mapM_ (saferEvalIO f) xs

  *Main> let f x = (1%1)/((x-1)*(x-2))
  *Main> makeGraph_ f [0..4]
  (0 % 1,1 % 2)
  (3 % 1,1 % 2)
  (4 % 1,1 % 6)

> gunc x = (1%1)/((x-1)*(x-2))

> main = do
>   arg <- getLine
>   let d = read arg :: Int
>   if d < 0 
>     then return ()
>     else
>       makeGraph_ gunc $ take d [0..]


