> {-# LANGUAGE BangPatterns #-}

catchEither.lhs

> import Control.Monad.Catch
> import Control.Exception 
>          ( ArithException (..)
>          , evaluate
>          )
> import Data.Ratio

> func x = 1%x
>
> saferEvalIO :: Show a => (t -> a) -> t -> IO ()
> saferEvalIO f x = do
>   print $ f x
>   `catch` \RatioZeroDenominator -> print $ "zero division"
>
> saferEvalEither 
>   :: (t -> b) -> t -> Either SomeException b
> saferEvalEither f x = do
>   fx <- try (return $ f x)
>   case fx of
>     Right val -> Right val
>     Left  ex  -> Left  ex
>
> saferEvalEither' f x = do
>   fx <- try (let !val = f x in return val)
>   case fx of 
>     Right val -> Right val
>     Left  ex  -> Left  ex

  *Main> try (evaluate $ func 0) 
           :: IO (Either SomeException (Ratio Int))
  Left Ratio has zero denominator
  *Main> try (evaluate $ func 1) 
           :: IO (Either SomeException (Ratio Int))
  Right (1 % 1)
