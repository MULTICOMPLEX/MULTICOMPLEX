import sys

import numpy as np
np.set_printoptions(threshold=np.inf)
np.set_printoptions(suppress=True) #prevent numpy exponential 

from typing import Tuple,Callable
from scipy.integrate import quadrature, odeint, solve_ivp
import matplotlib.pyplot as plt


def func(y, x):
    return  x * np.sqrt(np.abs(y)) + np.sin(x * np.pi/2)**3 - 5 * (x > 2)

def func2(x, y):
    return  x * np.sqrt(np.abs(y)) + np.sin(x * np.pi/2)**3 - 5 * (x > 2)


def mc_int(func: Callable, domain: Tuple, n_samples: int):
    samples = np.random.uniform(low=domain[0], high=domain[1], size=(n_samples,))
    #print(samples)
    volume = abs(domain[1] - domain[0])
    return np.mean(func(samples)) * volume


def mc_ode_solve(func, y0, t, n_samples=1000):
    sols = [y0]
    for lo, hi in zip(t[:-1], t[1:]) :
        part_func = lambda v: func(x=v, y=sols[-1])
        
        sols.append(sols[-1] + mc_int(part_func, (lo, hi), n_samples=n_samples))
    return np.asarray(sols)


base2 = np.linspace(-4, 5, 501)
y0 = 4.
ys_mc = []

ys_mc.append(mc_ode_solve(func, y0, base2))

c = np.vstack((base2, ys_mc)).T

print(c)


y_ode = solve_ivp(func2, (-4,5), [y0], method='LSODA', t_eval = base2, dense_output=True)
y_ode3 = solve_ivp(func2, (-4,5), [y0], method='RK23', t_eval = base2, dense_output=True)
z = y_ode.sol(base2)
z1 = y_ode3.sol(base2)
y_ode2 = odeint(func, y0, base2)

#print(y_ode)

fig, ax = plt.subplots()

ax.plot(base2, y_ode2, 'k', label="odeint", linewidth=1, color = 'green')
ax.plot(base2, z.T, 'k', label="solve_ivp method LSODA", linewidth=1, color = 'red')
ax.plot(base2, z1.T, 'k', label="solve_ivp method RK23", linewidth=1, color = 'orange')
ax.plot(base2, np.asarray(ys_mc).flatten(), 'k', label="monte carlo ode solve", linewidth=1, color = 'blue')

ax.legend()
plt.show()



