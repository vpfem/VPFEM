import numpy as np
from sympy import *

# define symbols
c, s, D11, D22, D33, D12, D13, D23 = symbols('c s D11 D22 D33 D12 D13 D23')

# transformation matrix
T = np.array([[c**2, s**2, c*s],
              [s**2, c**2, -c*s],
              [-2*c*s, 2*c*s, c**2-s**2]])
# mat matrix 
D = np.array([[D11, D12, D13],
              [D12, D22, D23],
              [D13, D23, D33]])
T2 = T.transpose()

result = (T2.dot(D)).dot(T)

print(T2)
print(D)
print(T)
print(simplify(result))
print(result[0][1] == result[0][1])

