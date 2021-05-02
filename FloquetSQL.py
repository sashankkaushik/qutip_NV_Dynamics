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

N = 6
J = N/2
k = 20
alpha = 5e-1
epsilon = 5e-6

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

psi0 = (basis(N+1,(N+1)//2)+basis(N+1,N)).unit() #Entangled
#psi0 = (basis(N+1,0)+basis(N+1,N)).unit() #GHZ                  
#psi0 = Qobj(su2coherent(N))
#print(psi0)
# Propagator
sx = jmat(J, 'x')
sy = jmat(J, 'y')
sz = jmat(J, 'z')
sm = jmat(J, '-')


Ularm = (-1j*alpha*sz).expm()
Unl = (-1j*k*sy*sy/(N+1)).expm()
#Unl = Unl/np.abs(np.linalg.det(Unl))
U = Unl*Ularm

#Ularm1 = (-1j*(alpha+epsilon)*sz).expm()
Unl1 = (-1j*(k+epsilon)*sy*sy/(N+1)).expm()
U1 = Unl1*Ularm
Udag = U1.dag()
#prop1 = Qobj(np.linalg.matrix_power(U,1024))
#print(np.abs(np.linalg.det(prop1)))
timelist = np.geomspace(2, 262144, num=18, dtype=int)
Nlist = np.array([2,4,11,16,32,64,128,256,512]) #np.geomspace(2, 256, num=8, dtype=int) #
#timelist1 = arange(1, 10000, 100)
qfi = []
qfiu = []
## Kicked evolution QFI calculation
for tmax in timelist:
    prop1 = Qobj(np.linalg.matrix_power(U,tmax)) #Qobj((-1j*k*sy*sy*tmax).expm())#
    prop2 = Qobj(np.linalg.matrix_power(Udag,tmax)) #Qobj((-1j*(k+epsilon)*sy*sy*tmax).expm()) #
    prop2 = prop2.dag()
    psif = Qobj(prop1 * psi0)
    psi = Qobj(prop2 * psif.unit())
    psi = psi.unit()
    F, innp = fidel(psi0,psi)
    Q = 8*(1-F)/(epsilon**2)
    qfi.append(Q)
    print(np.ceil(np.log2(tmax))+1, Q)

"""tmax = 100
for n in Nlist:
    psi0 = (basis(n+1,n//2)+basis(n+1,n)).unit() #Entangled
    #psi0 = (basis(n+1,0)+basis(n+1,n)).unit() #GHZ                  
    #psi0 = Qobj(su2coherent(n))
    sy = jmat(n/2, 'y')
    sz = jmat(n/2, 'z')
    Unl = (-1j*k*sy*sy/(n+1)).expm()
    Ularm = (-1j*alpha*sz).expm()
    U=Unl*Ularm
    Unl1 = (-1j*(k+epsilon)*sy*sy/(n+1)).expm()
    U1 = Unl1*Ularm
    prop1 = Qobj(np.linalg.matrix_power(U,tmax)) #Qobj((-1j*k*sy*sy*tmax).expm()) #
    prop2 = Qobj(np.linalg.matrix_power(U1,tmax)) #Qobj((-1j*(k+epsilon)*sy*sy*tmax).expm()) #
    prop2 = prop2.dag()
    psif = Qobj(prop1 * psi0)
    psi = Qobj(prop2 * psif.unit())
    psi = psi.unit()
    F, innp = fidel(psi0,psi)
    Q = 8*(1-F)/(epsilon**2)
    qfi.append(Q)
    print(n, Q) #np.ceil(np.log2(n))"""

  
## Unperturbed evolution QFI calculation
"""for n in Nlist:
    psi0 = (basis(n+1,0)+basis(n+1,n)).unit() #GHZ                  
    #psi0 = Qobj(su2coherent(n))
    sy = jmat(n/2, 'y')
    Q = 4*variance(tmax*sy*sy/(n+1), psi0)
    qfiu.append(np.abs(Q))"""

#qfi = np.asarray(qfi, dtype=float)
#print(qfi)
qfi = np.asarray(qfi, dtype=float)
#np.savetxt("/qfi.CSV",b,delimiter=',')

axes = plt.gca()

plt.plot(np.log10(timelist), np.log10(qfi), '.')
"""lognlin = np.log10(Nlist[3:])
logQlin = np.log10(qfi[3:])
lincoeff = np.polyfit(lognlin, logQlin, 1)
linpoly = np.poly1d(lincoeff)
yl = linpoly(lognlin)
plt.plot(lognlin, yl)
print(lincoeff)
logQquad = np.log10(qfi[:4])
lognquad = np.log10(Nlist[:4]) #logtquad = np.log10(timelist[2:])
quadcoeff = np.polyfit(lognquad, logQquad, 1)
quadpoly = np.poly1d(quadcoeff)
yq = quadpoly(lognquad)
plt.plot(lognquad, yq)
print(quadcoeff)"""
plt.grid()
plt.xlabel("N") #plt.xlabel("t") #
plt.ylabel("QFI [t = "+ str(tmax)+"]") #plt.ylabel("QFI [N = "+ str(N)+"]") #"""