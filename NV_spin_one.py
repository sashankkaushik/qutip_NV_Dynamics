#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 23:17:47 2020

@author: sashankkaushik

Spin-1
Current setting: Rabi
"""

import numpy as np
from numpy import pi
from qutip import *
from qutip.qip.device import Processor
from qutip.bloch import Bloch
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from qutip.ipynbtools import plot_animation
from qutip.states import qutrit_basis

s11, s22, s33, s12, s23, s31 = qutrit_ops()
is12 = s12 * -1j
is23 = s23 * -1j
sz = s11 - s33
sx = (np.sqrt(2)/2)*(s12 + s12.dag() + s23 + s23.dag())
sy = (np.sqrt(2)/2)*(is12 + is12.dag() + is23 + is23.dag())
I = s11+s22+s33
#print(sx, sy, sz)

h = 1.054*(10**(-34))
D = 2*pi*2.87*10**9
B1 = 0.00007 #Max simulable ratio of Bz/By (i.e., B0/B1) is 1000. Time steps in tlist must be increased by factor of 10 accordingly 
B0 = 0.007
gmr = 28.024951*10**9 #gyromagnetic ratio in Hz/T 

options = Options(store_states=True, store_final_state=True)
c_op_list = []

omega0 = gmr * B0 * 2 * pi
omega1 = gmr * B1 * 2 * pi
tlist = np.linspace(0,pi/(np.sqrt(2)*omega1),10000)

NV = Processor(N=1, dims=[3], spline_kind="cubic" ) #t1 = 10**-3, t2 = 10**-6 (characteristic of NV)
NV.add_control(D*(sz*sz-(2/3)*I) + omega0*sz, label="Bz field") #D*(sz*sz-(2/3)*I) added for zfs term
NV.add_control(omega1*sx, label="Bx pulse")
NV.add_control(omega1*sy, label="By pulse")
NV.pulses[0].coeff = np.array([-1. for t in tlist])
cf1 = np.array([-np.cos((omega0-D)*t) for t in tlist[:5000]])
cf1 = np.append(cf1,np.zeros(5000))
NV.pulses[1].coeff = cf1 #change to np.array([-np.cos((omega0-D)*t) for t in tlist]) for Rabi
cf2 = np.array([np.sin((omega0-D)*t) for t in tlist[:5000]])
cf2 = np.append(cf2,np.zeros(5000))
NV.pulses[2].coeff = cf2 #change to np.array([np.sin((omega0-D)*t) for t in tlist]) for Rabi
NV.pulses[0].tlist = tlist
NV.pulses[1].tlist = tlist
NV.pulses[2].tlist = tlist

#NV.plot_pulses(title="Pulse amplitudes")

init=basis(3,1)
basis1 = ket2dm(basis(3,0))
basis0 = ket2dm(basis(3,1))
basism1 = ket2dm(basis(3,2))
result = NV.run_state(init, analytical=False,  e_ops = [sx,sy,sz,basis1,basis0,basism1], options=options)
expectation = result.expect
states = result.states
#print(states)
#plot_expectation_values(result) #Plots Pauli matrices expectation values

f1_norm = "{:0.3f}".format(omega1/(2*pi*10**6))
f0_norm = "{:0.3f}".format(omega0/(2*pi*10**6))
plt.figure(dpi=800)
axes = plt.gca()
axes.set_ylim([-1.1,1.1])
plt.plot(tlist,expectation[0])
plt.plot(tlist,expectation[1])
plt.plot(tlist,expectation[2])
plt.grid()
plt.xlabel("t")
plt.ylabel("Expectation values (f0 = "+str(f0_norm)+" MHz)")
plt.legend(["sx", "sy", "sz"], loc=3)

plt.figure(dpi=800)
axes = plt.gca()
axes.set_ylim([-0.1,1.1])
plt.plot(tlist,expectation[3])
plt.plot(tlist,expectation[4])
plt.plot(tlist,expectation[5])
plt.grid()
plt.xlabel("t")
plt.ylabel("Transition probabilities (f0 = "+str(f0_norm)+" MHz)")
plt.legend(["P(m=1)", "P(m=0)", "P(m=-1)"], loc=3)