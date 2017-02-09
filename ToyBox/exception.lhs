exception.lhs

> import Control.Exception (ArithException (..))
> import Control.Monad.Catch

> main :: IO ()
> main = do
>   denominator <- getLine
> --   let x = 10 `safeDiv` 0
>   let x = 10 `safeDiv` (read denominator)
>   case x of
>     Left  a -> print (a :: SomeException)
>     Right b -> print b
>
> safeDiv 
>   :: MonadThrow m =>
>      Int -> Int -> m Int
> a `safeDiv` b =
>   if b == 0 
>     then throwM RatioZeroDenominator
>     else return $ a `div` b

