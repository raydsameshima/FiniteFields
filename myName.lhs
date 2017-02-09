myName.lhs

> main = do
>   putStrLn "Ray D. Sameshima"
>   putStrLn "New York City College of Technology, CUNY"

This Haskell code is written in so called
  literal
style.
In literal-haskell, the normal lines will be ignored by compiler.
Only the lines which are followed by
  >
symbol will be read.

Haskell
An advanced, purely functional programming language
(standardized, general purpose, statically typed)

One of the striking features is its lazy evaluation.
(Call-by-need strategy cf. call-by-name or "eager")
www.haskell.org

Imperative v.s. Declarative

Imperative (e.g., C/C++, Java, Python etc)
procedure ~ how to do it
loop

Declarative (e.g., Haskell)
definition ~ what is it
recursion

Example:

> mySum []     = 0
> mySum (n:ns) = n + mySum ns

