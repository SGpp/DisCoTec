"""module provides Gauss-Laguerre weights and knots for integration"""

from math import exp
from operator import mul
from scipy.special.orthogonal import l_roots

 # This routine is adapted from GENE  routine.
 # n size of weights/knots array
 # x2 upper bound
 # x  Gauss-Laguerre knots
 # w  Gauss-Leguerre weights
 
def gaulag(x2,n):
    """=========================================================
    Purpose : Compute the zeros of Laguerre polynomial Ln(x)
    in the interval [0,oo] and scale to [0,x2], and the corresponding
    weighting coefficients for Gauss-Laguerre quadrature
    integration
         Input :
         x2   --- Upper bound
         n    --- Order of the Laguerre polynomial
         Output :
         x(n) --- Zeros of the Laguerre polynomial
         w(n) --- Corresponding weighting coefficients
    ========================================================="""

    # find Laguerre polynomials roots and weights
    (x, w) = l_roots(n)

    # scale to [0,x2]
    w = map(mul,w,map(exp,x))
    fac = x2/sum(w)
    w = [e*fac for e in w]
    x = [e*fac for e in x]
    
    return (x, w)
