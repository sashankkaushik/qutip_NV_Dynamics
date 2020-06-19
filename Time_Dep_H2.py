#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:27:02 2020

@author: sashankkaushik
"""

import numpy as np
from numpy import pi
from numpy import linalg as LA 
from qutip import *
import matplotlib.pyplot as plt
g    = 2   #g constant
me = 9.1*(10**(-31))  #mass of an electron
e = 1.6*(10**(-19))   #charge of an electron
theta = 0.2 * pi       # qubit angle from sigma_z axis (toward sigma_x axis)
gamma1 = 0.0     # qubit relaxation rate
gamma2 = 0.0     # qubit dephasing rate
h = 1.054*(10**(-34))
Bx = 0.5
By = 0.5
Bz = 0.01
gmr = 28024.951 #gyromagnetic ratio in Mhz/T 

options = Options()
options.nsteps = 5000
# initial state
psi0 = basis(2,1)
#psi0 = (np.sqrt(2)/2)*basis(2,1)+(np.sqrt(2)/2)*basis(2,0)
t = np.linspace(0,500,10000)

# Hamiltonian
sx = sigmax()
sy = sigmay()
sz = sigmaz()
sm = sigmam()
H0 = (Bz*sz)
c_op_list = []  # collapse operators
 # evolve and calculate expectation values
omega = gmr * Bx * (10 ** (-6))/ pi
def H1_Coeff(t,args):
    return np.cos(omega*t)
H1 = (Bx*sx)
def H2_Coeff(t,args):
    return np.sin(omega*t)
H2 = -(By*sy)
H = [H0,[H1,H1_Coeff],[H2, H2_Coeff]]
output = mesolve(H, psi0, t, c_op_list,[sx,sy,sz], options=options)       
sx = output.expect[0]
sy = output.expect[1]
sz = output.expect[2]

plt.plot(t,sx)
plt.plot(t,sy)
plt.plot(t,sz)
plt.xlabel("t")
plt.ylabel(['Expectation values (omega =', str(omega), ')'])
plt.legend(["sx", "sy", "sz"], loc=4)
"""sphere=Bloch()
sphere.add_points([sx,sy,sz])
sphere.show()"""