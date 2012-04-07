"This module provides support for calculating AHP weights and consistency ratios (CR)"

from decimal import Decimal
from numpy import array
from numpy.linalg import eig


def calculateWeights(arr, rounding=4):
   """Given pairwise comparisons array, calculate weights using the AHP method

   >>> arr = array([ [1, 1./3, 5], [3,1,7], [1./5,1./7,1] ])
   >>> vector = calculateWeights(arr, rounding=4)
   >>> vector
   [Decimal('0.2790'), Decimal('0.6491'), Decimal('0.0719')]

   """

   PLACES = Decimal(10) ** -(rounding)
   
   # get eigenvalues and vectors
   evas, eves = eig(arr)

   # get primary eigenvalue and vector
   eva = max(evas)
   eva_idx = evas.tolist().index(eva)
   eve = eves.take((eva_idx,), axis=1)

   # priority vector = normalized primary eigenvector

   normalized = eve / sum(eve)

   # turn into list of real part values
   vector = [abs(e[0]) for e in normalized]

   # return nice rounded Decimal values with labels
   return [ Decimal( str(v) ).quantize(PLACES) for v in vector ]


def calculateConsistency(arr):
   "given pairwise comparisons array, calculate consistency ratio (CR) for the comparisons"



if __name__ == "__main__":
   import doctest
   doctest.testmod()