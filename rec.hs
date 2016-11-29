-- rec.hs

import System.Environment (getArgs)
import System.IO 
import Data.Ratio

import Multivariate (table2ratf)

main = do
  args <- getArgs
  let fileName = head args

  handle <- openFile fileName ReadMode
  contents <- hGetContents handle
  putStrLn $ "The data from " ++ fileName ++ " is: "
  putStrLn contents
  putStrLn "The coefficients: "
  let d = read contents :: [[Ratio Integer]]
  putStrLn . show $ table2ratf d
  hClose handle
  
