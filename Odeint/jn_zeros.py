from scipy.special import jv, jn_zeros
from scipy.optimize import brentq
import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.16f}".format(x)})

n = 1

x = jn_zeros(n, 6)

print(x,"\n")

BESSEL_ROOTS = [jn_zeros(m, 10) for m in range(10)]
print(BESSEL_ROOTS, "\n")

def Wave_function(x, m):
    psi = jv(m, x)
    return psi

def find_all_zeroes(x,y, m):
    all_zeroes = []
    s = np.sign(y)
    for i in range(len(y)-1):
        if s[i]+s[i+1] == 0:
            zero = brentq(Wave_function, x[i], x[i+1], m)
            all_zeroes.append(zero)
    return all_zeroes

en = np.linspace(0, 20, 20)   

for m in range(10):
    psi_b = jv(m, en)   
    E_zeroes = find_all_zeroes(en, psi_b, m)   
    print(E_zeroes)

MODES = (
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 1),
    (2, 2),
    (2, 3)
)

def get_vmin_vmax(m, n):
    vmax = np.max(jv(m, np.linspace(0, BESSEL_ROOTS[m][n], 50)))
    return -vmax, vmax

m, n = MODES[4]
vmin, vmax = get_vmin_vmax(m, n)

print("\n", vmin)


  
    