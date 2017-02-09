catchIO.lhs

This is IO-wrapping.

> import Control.Exception 
>          ( ArithException (..)
>          , catch
>          )
> import Data.Ratio

> func x = 1%x
>
> saferEvalIO f x = do
>   print $ f x
>   `catch` \RatioZeroDenominator -> print $ "zero division"

