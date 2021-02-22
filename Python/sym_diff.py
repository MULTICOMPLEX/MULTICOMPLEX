from sympy import *
x, y, z, p = symbols("x y z p", complex=True)
from sympy import I
import time

import numpy as np  
import matplotlib.pyplot as plt

#from sympy.printing.theanocode import theano_function
#import theano.tensor as T

#x = Symbol('x')
#y = Symbol('y')
#z = Symbol('z')

#print (cos(z/integrate(integrate(integrate(exp(x+1) * sqrt(y) -2, x),y),x)))

exa = diff(diff(x*x * x * y*y, x),x)
exb = diff(diff(x*x + x * y*y, x),y)
exc = diff(diff(x*x + x * y*y, y),x)
exd = diff(diff(x*x + x * y*y, y),y)

print (exa)
print (exb)
print (exc)
print (exd)

#print (sym.diff(sym.sin(x), x))

print ("\n")

f1 = lambdify((x,y), exa); print(f1(.7,.9))
f2 = lambdify((x,y), exb); print(f2(.7,.9))
f3 = lambdify((x,y), exc); print(f3(.7,.9))
f4 = lambdify((x,y), exd); print(f4(.7,.9))

def HermitePoly(x, n):
  # uses the identity H_n(x) == (-1)^n exp(x^2) (d/dx)^n exp(-x^2)
  p = exp( -(x*x));
  return (pow(-1, n) * diff(p, x, n) / p)

print ("\n")

for i in range(7):
    print ("H_" , i , "(z) == " , HermitePoly(z,i))

x, y = symbols('x y')
#p = x * exp(x)-1
p = 1-exp(x)
#p = exp(x) - x
#r = solve(p, x)

bounds = lambda i: (3.14*i, 3.14*(i + 1))
r = nsolve(p, bounds(100), solver='bisect', verify=False)
#r = nsolve(p, I)

#p.subs(x, r[0])
#print ("\n",r)
#print ("\n",p.subs(x, r[0]).evalf())
print ("\n", r)


eq = x*z*x+sqrt(x*z*y*x+y*z)*pow(y,6)*x*z
v = list(ordered(eq.free_symbols)); 

t0 = time.time()
k = hessian (eq,v)
t1 = time.time()

print ("\n",k)

print ("\nDuration ",t1-t0, "Sec")

fl = lambdify((x,y,z), k)
lf = fl(.7+0j,.9+0j,.2+0j)
t1 = time.time()

print("\n",lf)
print ("\nDuration + evaluation",t1-t0, "Sec")

my_array = np.array(lf)
t1 = time.time()
print("\nDeterminant",np.linalg.det(my_array))

print ("\nDuration + evaluation + determinant",t1-t0, "Sec")

exa = integrate(y*x * exp(z) -6*x* z*z, z)
print ("\n", exa)

from sympy.vector import CoordSys3D, scalar_potential
R = CoordSys3D('R')

conservative_field = 4*R.x*R.y*R.z*R.i + 2*R.x**2*R.z*R.j + 2*R.x**2*R.y*R.k

print ("\n",scalar_potential(conservative_field, R))

  

