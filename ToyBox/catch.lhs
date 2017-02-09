catch.lhs

This is IO-wrapping.

> import Control.Exception 
>          ( ArithException (..)
>          , catch
>          )
> import Data.Ratio

> func x = 1%x
>
> main = do
>   putStrLn "give me an input: "
>   input <- getLine
>   do
>     print $ (func (read input :: Int))
>     `catch` \RatioZeroDenominator -> print $ "zero division"
>   print "the end"
