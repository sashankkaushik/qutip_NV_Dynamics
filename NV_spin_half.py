#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:53:12 2020

@author: sashankkaushik

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

h = 1.054*(10**(-34))
B1 = 0.0007 #Max simulable ratio of Bz/By (i.e., B0/B1) is 1000. Time steps in tlist must be increased by factor of 10 accordingly 
B0 = 0.007
gmr = 28024.951*10**6 #gyromagnetic ratio in Hz/T 

sx = sigmax()
sy = sigmay()
sz = sigmaz()
sm = sigmam()

options = Options()
c_op_list = []

omega0 = gmr * B0 * 2 * pi
omega1 = gmr * B1 * 2 * pi
tlist = np.linspace(0,pi/omega1,1000)

NV = Processor(N=1,spline_kind="cubic" ) #t1 = 10**-3, t2 = 10**-6 (characteristic of NV)
NV.add_control(0.5*omega0*sz, label="Bz field")
NV.add_control(0.5*omega1*sx, label="Bx pulse")
NV.add_control(0.5*omega1*sy, label="By pulse")
NV.pulses[0].coeff = np.array([-1. for t in tlist])
cf1 = np.array([-np.cos(omega0*t) for t in tlist[:500]])
cf1 = np.append(cf1,np.zeros(2500)) 
NV.pulses[1].coeff = np.array([-np.cos(omega0*t) for t in tlist]) #change to np.array([-np.cos(omega0*t) for t in tlist]) for Rabi
cf2 = np.array([np.sin(omega0*t) for t in tlist[:500]])
cf2 = np.append(cf2,np.zeros(2500)) 
NV.pulses[2].coeff = np.array([np.sin(omega0*t) for t in tlist]) #change to np.array([np.sin(omega0*t) for t in tlist]) for Rabi
NV.pulses[0].tlist = tlist
NV.pulses[1].tlist = tlist
NV.pulses[2].tlist = tlist

#NV.plot_pulses(title="Pulse amplitudes")

basis0 = (2*basis(2,0)+basis(2,1)).unit()
result = NV.run_state(basis0, e_ops = [sx, sy, sz], options=options)
expectation = [result.expect[0], result.expect[1], result.expect[2]]
#print(omega0, omega1)
#plot_expectation_values(result) #Plots Pauli matrices expectation values

f1_norm = "{:0.3f}".format(omega1/(2*pi*10**6))
plt.figure(dpi=800)
axes = plt.gca()
axes.set_ylim([-1.0,1.0])
plt.plot(tlist,expectation[0])
plt.plot(tlist,expectation[1])
plt.plot(tlist,expectation[2])
plt.xlabel("t")
plt.ylabel('Expectation values (f1 ='+str(f1_norm)+' MHz)')
plt.legend(["sx", "sy", "sz"], loc=4)

### For plotting Pauli matrices expectation values wiht better resolution by splitting into 10 graphs ###
"""for i in range(0,9001,1000):
    plt.figure(dpi=800)
    axes = plt.gca()
    axes.set_ylim([-1.0,1.0])
    plt.plot(tlist[i:i+999],expectation[0][i:i+999])
    plt.plot(tlist[i:i+999],expectation[1][i:i+999])
    plt.plot(tlist[i:i+999],expectation[2][i:i+999])
    plt.xlabel("t")
    plt.ylabel(['Expectation values (omega1 =', str(omega1), ')'])
    plt.legend(["sx", "sy", "sz"], loc=4)"""

### Plots basic Bloch sphere representation of evolution ###
"""b = Bloch()
b.clear()
b.add_points(expectation, meth = 'l')
b.show()"""


### Generates frames for Bloch sphere animation (takes lot of time for len(tlist)>1000) ###
"""for i in range(len(expectation[0])):
    b = Bloch()
    b.point_color = ['r']
    b.clear()
    b.view = [-60-omega0*tlist[i]*180/pi, 30]
    b.add_vectors([expectation[0][i],expectation[1][i],expectation[2][i]])
    b.add_points([expectation[0][:i+1],expectation[1][:i+1],expectation[2][:i+1]], meth='l')
    b.save(dirc='Ramsey_T1',name='/Users/sashankkaushik/.spyder-py3/Ramsey_T1/bloch_'+str(i))"""
