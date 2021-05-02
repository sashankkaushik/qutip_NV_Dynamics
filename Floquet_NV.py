#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 00:30:36 2020

@author: sashankkaushik
"""

import numpy as np
import scipy as sp
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
from qutip.piqs import *

g = 2   #g constant
me = 9.1*(10**(-31))  #mass of an electron
e = 1.6*(10**(-19))   #charge of an electron
theta = 0.2 * pi       # qubit angle from sigma_z axis (toward sigma_x axis)
gamma1 = 0.0     # qubit relaxation rate
gamma2 = 0.0     # qubit dephasing rate
h = 1.054*(10**(-34))
gmr = 28024.951 #gyromagnetic ratio in Mhz/T 
D = 2*pi*28.7

N = 100
J = N/2
k = 30
alpha = 1e-1
epsilon = 1e-7

## initial state


def C(n,k): #Binomial coefficient
    return sp.special.binom(n,k)

def su2coherent(N, theta=pi/2, phi=pi/2): #Coherent
    return sum(np.sqrt(C(N,(N-M)//2))*(np.cos(theta/2)**((N+M)//2))*(np.sin(theta/2)*np.exp(1j*phi))**((N-M)//2)*basis(N+1,(N-M)//2).data for M in range(-N, N+1, 2))

def fidel(psiA, psiB):
    innp = psiA.overlap(psiB)
    fid = np.abs(innp)**2
    fid = np.clip(fid, 0.0, 1.0)
    return fid, innp

#psi0 = (basis(N+1,0)+basis(N+1,N)).unit() #GHZ              
psi0 = Qobj(su2coherent(N))
#print(psi0)
# Propagator
sx = jmat(J, 'x')
sy = jmat(J, 'y')
sz = jmat(J, 'z')
sm = jmat(J, '-')
I = Qobj(np.eye(N+1))

Ularm = (-1j*(alpha*sz+D*(sz*sz-(2/3)*I))).expm()
Unl = (-1j*k*sy*sy/(N+1)).expm()
#Unl = Unl/np.abs(np.linalg.det(Unl))
U = Unl*Ularm

Ularm1 = (-1j*((alpha+epsilon)*sz+D*(sz*sz-(2/3)*I))).expm()
U1 = Unl*Ularm1
Udag = U1.dag()
#prop1 = Qobj(np.linalg.matrix_power(U,1024))
#print(np.abs(np.linalg.det(prop1)))
timelist = np.geomspace(1, 262144, num=19, dtype=int)
#timelist1 = arange(1, 10000, 100)
qfi = []

## Kicked evolution QFI calculation
for tmax in timelist:
    prop1 = Qobj(np.linalg.matrix_power(U,tmax))
    prop2 = Qobj(np.linalg.matrix_power(Udag,tmax))
    psif = Qobj(prop1 * psi0)
    psi = Qobj(prop2 * psif.unit())
    psi = psi.unit()
    F, innp = fidel(psi0,psi)
    Q = 8*(1-F)/(epsilon**2)
    qfi.append(Q)
    print(np.ceil(np.log2(tmax))+1, Q)

"""prop1 = np.linalg.matrix_power(U,1024)
prop2 = np.linalg.matrix_power(Udag,1024)
psi = Qobj(prop2*prop1*psi0)
print(psi0)"""
    
## Unperturbed evolution QFI calculation
"""for t in timelist:
    Q = 4*variance(t*sz, psi0)
    qfi.append(Q)"""

#qfi = np.asarray(qfi, dtype=float)

axes = plt.gca()
plt.xscale("log")
plt.yscale("log")

"""logtlin = np.log10(timelist[6:11])
logtquad = np.log10(timelist[11:26])
logQlin = np.log10(qfi[6:11])
logQquad = np.log10(qfi[11:26])
lincoeff = np.polyfit(logtlin, logQlin, 1)
linpoly = np.poly1d(lincoeff)
yl = linpoly(logtlin)
quadcoeff = np.polyfit(logtquad, logQquad, 1)
quadpoly = np.poly1d(quadcoeff)
yq = quadpoly(logtquad)"""

plt.plot(timelist, qfi, '.')
#plt.plot(logtlin, yl)
#plt.plot(logtquad, yq)
plt.grid()
plt.xlabel("t")
plt.ylabel("QFI for NV [N = "+ str(N)+"]")

#print(lincoeff)
#print(quadcoeff)"""