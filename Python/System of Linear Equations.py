
# https://stackabuse.com/solving-systems-of-linear-equations-with-pythons-numpy/

# The ultimate goal of solving a system of linear equations is to find the values of the 
# unknown variables. Here is an example of a system of linear equations with two unknown variables, 
# x and y:

# Equation 1:

# 4x  + 3y = 20
# -5x + 9y = 26


# To solve the above system of linear equations, we need to find the values of the x and y variables. 
# There are multiple ways to solve such a system, such as Elimination of Variables, 
# Cramer's Rule, Row Reduction Technique, 
# and the Matrix Solution. In this article we will cover the matrix solution.

# In the matrix solution, the system of linear equations to be solved is represented in the form 
# of matrix AX = B. For instance, 
# we can represent Equation 1 in the form of a matrix as follows:

'''
A = [[ 4   3]
     [-5   9]]

X = [[x]
     [y]]

B = [[20]
     [26]]
'''

# To find the value of x and y variables in Equation 1, 
# we need to find the values in the matrix X. 
# To do so, we can take the dot product of the inverse of matrix A, and the matrix B as shown below:

# X = inverse(A).B

import numpy as np

m_list = [[4, 3], 
         [-5, 9]]
A = np.array(m_list)

inv_A = np.linalg.inv(A)

print(inv_A) 

B = np.array([20, 26])
X = np.linalg.inv(A).dot(B)

print(X) # [2. 4.]

# 20x + 10y = 350
# 17x + 22y = 500

A = np.array([[20, 10], [17, 22]])
B = np.array([350, 500])
X = np.linalg.solve(A,B)


print(X) # [10. 15.]


# 4x + 3y + 2z = 25
# -2x + 2y + 3z = -10
# 3x -5y + 2z = -4

A = np.array([[4, 3, 2], [-2, 2, 3], [3, -5, 2]])
B = np.array([25, -10, -4])
X = np.linalg.inv(A).dot(B)

print(X) # [ 5.  3. -2.]

A = np.array([[-6, 3], [4, 5]])
B = np.array([6, 6])
X = np.linalg.solve(A,B)
print(X)


import numpy as np
from scipy.optimize import fsolve

def myFunction(z):
   x = z[0]
   y = z[1]

   F = np.empty((2))
   F[0] = -12*x + 3*y 	
   F[1] = 4*x - 1*y 	
   return F

zGuess = np.array([5,2])
z = fsolve(myFunction,zGuess)

print(z)

import sympy as sym

x,y = sym.symbols('x,y')

f = sym.Eq(-6*x+3*y,6*x)
g = sym.Eq(4*x + 5*y,6*y)

z = sym.solve([f,g],(x,y))

print(z)


x,y,z = sym.symbols('x,y,z')

f = sym.Eq(2*x,-x)
g = sym.Eq(4*y + 5*z,-y)
h = sym.Eq(4*y + 3*z,-z)

z = sym.solve([f,g,h],(x,y,z))
print(z)




