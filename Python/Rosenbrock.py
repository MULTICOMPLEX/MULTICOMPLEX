import numpy as np
from scipy.optimize import minimize
"""
Example of minimizing the Rosenbrock function with different
methods using SciPy.

"""


def rosenbrock(x):
    """The Rosenbrock function"""
    return sum(100.0 * (x[1:] - x[:-1]**2.0)**2.0 + (1 - x[:-1])**2.0)


# Randomize initial starting points for N-dimensions
# when  1.0 <= x <= 10.0
N = 9
x0 = np.random.uniform(low=1.0, high=10.0, size=N)
print(x0)

# Nelder-Mead method
print('\nResults of minimizing with Nelder-Mead method: ')
res_NM = minimize(rosenbrock, x0, method='Nelder-Mead',
                  options={'disp': True, 'maxiter': 10**4})
print('Solution: ')
print(res_NM.x)

# Sequential Least SQuares Programming (SLSQP) method
print('\nResults of minimizing with SLSQP method: ')
res_SLSQP = minimize(rosenbrock, x0, method='SLSQP',
                     options={'disp': True, 'maxiter': 10**4})
print('Solution: ')
print(res_SLSQP.x)

# Powell method
print('\nResults of minimizing with Powell method: ')
res_Powell = minimize(rosenbrock, x0, method='Powell',
                      options={'disp': True, 'maxiter': 10**4})
print('Solution: ')
print(res_Powell.x)

# Hessianmethod

def rosen(x):

    """The Rosenbrock function"""

    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)


def rosen_der(x):

    xm = x[1:-1]

    xm_m1 = x[:-2]

    xm_p1 = x[2:]

    der = np.zeros_like(x)

    der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)

    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])

    der[-1] = 200*(x[-1]-x[-2]**2)

    return der



def rosen_hess(x):

    x = np.asarray(x)

    H = np.diag(-400*x[:-1],1) - np.diag(400*x[:-1],-1)

    diagonal = np.zeros_like(x)

    diagonal[0] = 1200*x[0]**2-400*x[1]+2

    diagonal[-1] = 200

    diagonal[1:-1] = 202 + 1200*x[1:-1]**2 - 400*x[2:]

    H = H + np.diag(diagonal)

    return H
    
res_Newton = minimize(rosen, x0, method='Newton-CG',

               jac=rosen_der, hess=rosen_hess,

               options={'xtol': 1e-8, 'disp': True})

print('Solution: ')
print(res_Newton.x)             

d = ({'Nelder-Mead': res_NM.nit, 'SLSQP': res_SLSQP.nit,
      'Powell': res_Powell.nit,'Newton': res_Newton.nit})
print('\nMethod with the least amount of iterations was: ' + min(d, key=d.get))
