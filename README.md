# Finite Fields 
Here is an attempt to make my version of functional reconstruction, i.e., interpolation of polynomials and rational functions on discrete samplings.

## Finite fields Z_p
### Ffield.lhs
For rational number reconstruction.
Extended Euclidean algorithm, chinese remainder theorem.
It passed QuickCheck.

# Function reconstruction
## Univariate (1-variable) cases
Newton interpolation and Thiele interpolation.

### Univariate.lhs
Following "The Haskell Road to Logic, Maths and Programming (Kees Doets, Jan van Eijck), arbitrary precision version of Newton algorithm (finite difference analysis).
It determines the degree of the input polynomial (or the outputs), then constructs the list of coefficients. 

For rational function, it uses reciprocal (inverse) difference analysis, so called Thiele interpolation.
Using arbitrary precision rational numbers and Maybe wrapping it can treat 1/0 type infinity in intermediate steps.
Once it hits a pole at the sampling phase, it simply halts.

### GUniFin.lhs
Non-sequential inputs, finite fields.
Time consuming but still works.

## Multivariate cases
More than 2 variables.
Here merely is 2 variables cases.

### GMulFin.lhs
Heavy recursion (adding one layer more).
Up to certain total degree it works, say x^2*y^2/(1+y)^3, but it also requires clever choice of samplings.

