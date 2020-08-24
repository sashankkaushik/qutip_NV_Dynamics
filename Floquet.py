#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 00:30:36 2020

@author: sashankkaushik
"""

import numpy as np
from numpy import pi
from numpy import linalg as LA 
from qutip import *
from qutip.qip.device import Processor
from qutip.bloch import Bloch
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from qutip.ipynbtools import plot_animation

g = 2   #g constant
me = 9.1*(10**(-31))  #mass of an electron
e = 1.6*(10**(-19))   #charge of an electron
theta = 0.2 * pi       # qubit angle from sigma_z axis (toward sigma_x axis)
gamma1 = 0.0     # qubit relaxation rate
gamma2 = 0.0     # qubit dephasing rate
h = 1.054*(10**(-34))
gmr = 28024.951 #gyromagnetic ratio in Mhz/T 

j = 1
k = 30.0
alpha = pi/1000
epsilon = 0.000001
theta = pi/2
phi = pi/2
sine = np.sin(theta/2)
cosine = np.cos(theta/2)
phase = exp(1j*phi)
# initial state
#psi0 = (basis(3,0)+basis(3,2)).unit() #GHZ
psi0 = (sine**2)*(phase**2)*basis(3,0)+np.sqrt(2)*sine*cosine*phase*basis(3,1)+(cosine**2)*basis(3,2) #Coherent


# Propagator
sx = spin_Jx(j)
sy = spin_Jy(j)
sz = spin_Jz(j)
sm = spin_Jm(j)
Ularm = propagator(alpha*sz, 1)
Unl = propagator(k*sy*sy/(2*j), 1)
U = Unl*Ularm

Ularm1 = propagator((alpha+epsilon)*sz, 1)
Unl1 = propagator(k*sy*sy/(2*j), 1)
U1 = Unl1*Ularm1
Udag = U1.dag()

timelist = arange(1,100,1)
psi = [psi0]
qfi = []

## Kicked evolution QFI calculation
for tmax in range(1,100,1):
    tlist = arange(1, tmax)
    tlist1 = arange(1, 2*tmax)

    for t in tlist:
        psi.append(U*psi[t-1])

    for t in tlist:
        psi.append(Udag*psi[t+len(tlist)-1])

    #expectation = [expect(sx, psi), expect(sy, psi), expect(sz, psi)]
    F = fidelity(psi0, psi[2*len(tlist)])
    Q = 8*(1-np.sqrt(F))/(epsilon**2)
    qfi.append(Q)
    
## Unperturbed evolution QFI calculation
"""for t in timelist:
    Q = 4*variance(t*sz, psi0)
    qfi.append(Q)"""

timelist = np.asarray(timelist, dtype=float)
qfi = np.asarray(qfi, dtype=float)
#logt = np.log10(timelist[1:])
#logQ = np.log10(qfi[1:])
axes = plt.gca()
plt.xscale("log")
plt.yscale("log")
#coefficients = np.polyfit(logt, logQ, 1)
#polynomial = np.poly1d(coefficients)
#ys = polynomial(logt)
plt.plot(timelist, qfi, '.')
#plt.plot(logt, ys)
plt.grid()
plt.xlabel("log(t)")
plt.ylabel("log(QFI)")
#print(coefficients)
#plt.legend(["sx", "sy", "sz"], loc=3)

#print(qfi)