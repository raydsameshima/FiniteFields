safeEval.lhs

> import Prelude hiding(catch)
> -- import Control.Exception (ArithException (..))
> import Control.Monad.Catch
>
> import Data.Ratio
>
> catchError = catch (error "error occur")
>                    (\e -> return (e :: SomeException)
>                    >> print "CAUGHT ERROR.")
>

