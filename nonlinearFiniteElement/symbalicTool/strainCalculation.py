#! /usr/bin/python3

import numpy as np
from sympy import *

# define symbols
N1x, N2x, N3x, N4x, N1y, N2y, N3y, N4y, d1, d2, d3, d4, d5, d6, d7, d8 = symbols('N1x N2x N3x N4x N1y N2y N3y N4y d1 d2 d3 d4 d5 d6 d7 d8')

# transformation matrix
N = np.array([[N1x, 0, N2x, 0, N3x, 0, N4x, 0],
              [0, N1y, 0, N2y, 0, N3y, 0, N4y],
              [N1y, N1x, N2y, N2x, N3y, N3x, N4y, N4x]])
# mat matrix 
D = np.array([d1,d2,d3,d4,d5,d6,d7,d8])
D = D.transpose()

result = (N).dot(D)


print(simplify(result))
#print(result[0][1] == result[0][1])

