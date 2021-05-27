# -*- coding: utf-8 -*-
"""
Created on Mon Sep 03 21:02:59 2014
 
@author: Pero
 
1D Schrödinger Equation in a finite square well.
 
Program calculates bound states and energies for a finite potential well. It will find eigenvalues
in a given range of energies and plot wave function for each state.
 
For a given energy vector e, program will calculate 1D wave function using the Schrödinger equation
in a finite square well defined by the potential V(x). If the wave function diverges on x-axis, the
energy e represents an unstable state and will be discarded. If the wave function converges on x-axis,
energy e is taken as an eigenvalue of the Hamiltonian (i.e. it is alowed energy and wave function
represents allowed state). 
 
Program uses differential equation solver &amp;amp;amp;quot;odeint&amp;amp;amp;quot; to calculate Sch. equation and optimization
tool &amp;amp;amp;quot;brentq&amp;amp;amp;quot; to find the root of the function. Both tools are included in the Scipy module.
The following functions are provided:
 
    - V(x) is a potential function of the square well. For a given x it returns the value of the potential
    - SE(psi, x) creates the state vector for a Schrödinger differential equation. Arguments are:
        psi - previous state of the wave function
        x - the x-axis
    - Wave_function(energy) calculates wave function using SE and &amp;amp;amp;quot;odeint&amp;amp;amp;quot;. It returns the wave-function
        at the value b far outside of the square well, so we can estimate convergence of the wave function.
    - find_all_zeroes(x,y) finds the x values where y(x) = 0 using &amp;amp;amp;quot;brentq&amp;amp;amp;quot; tool.    
 
Values of m and L are taken so that h-bar^2/m*L^2 is 1.
 
v2 adds feature of computational solution of analytical model from the usual textbooks. As a result,
energies computed by the program are printed and compared with those gained by the previous program.
"""
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq
import matplotlib.pyplot as plt
 
def V(x):
    """
    Potential function in the finite square well. Width is L and value is global variable Vo
    """
    L = 1
    if abs(x) > L:
        return 0
    else:
        return Vo
 
def SE(psi, x):
    """
    Returns derivatives for the 1D schrodinger eq.
    Requires global value E to be set somewhere. State0 is first derivative of the
    wave function psi, and state1 is its second derivative.
    """
    state0 = psi[1]
    state1 = -2.0*(V(x) - E)*psi[0]
    return [state0, state1]
 
def Wave_function(energy):
    """
    Calculates wave function psi for the given value
    of energy E and returns value at point b
    """
    global psi
    global E
    E = energy
    psi = odeint(SE, psi0, x)
    return psi[-1,0]
 
def find_all_zeroes(x,y):
    """
    Gives all zeroes in y = Psi(x)
    """
    all_zeroes = []
    s = np.sign(y)
    for i in range(len(y)-1):
        if s[i]+s[i+1] == 0:
            zero = brentq(Wave_function, x[i], x[i+1])
            all_zeroes.append(zero)
    return all_zeroes
    
N = 1000                  # number of points to take
psi = np.zeros([N,2])     # Wave function values and its derivative (psi and psi')
psi0 = [0,1]   # Wave function initial states
Vo = 20
E = 0.0                   # global variable Energy  needed for Sch.Eq, changed in function "Wave function"
b = 2                     # point outside of well where we need to check if the function diverges
x = np.linspace(-b, b, N)    # x-axis
 

# main program        


en = np.linspace(0, Vo, 100)   # vector of energies where we look for the stable states
 
psi_b = []      # vector of wave function at x = b for all of the energies in en
for e1 in en:
    psi_b.append(Wave_function(e1))     # for each energy e1 find the the psi(x) at x = b
    E_zeroes = find_all_zeroes(en, psi_b)   # now find the energies where psi(b) = 0 


for E in E_zeroes:
        print("%.2f" %E)
 
plt.grid()


for E in E_zeroes:
        Wave_function(E)
        plt.plot(x, psi[:,0], 'k', label="E = %.6f"%E, linewidth=1, color = 'red')

plt.legend()
plt.show()



