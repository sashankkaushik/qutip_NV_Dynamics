#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:53:12 2020

@author: sashankkaushik
"""
import numpy as np
from numpy import pi
from qutip import *
from qutip.qip.device import Processor
from qutip.bloch import Bloch
import matplotlib.pyplot as plt

h = 1.054*(10**(-34))
Bx = 0.0001
By = 0.0001
Bz = 0.1
gmr = 28024.951*10**6 #gyromagnetic ratio in Hz/T 

sx = sigmax()
sy = sigmay()
sz = sigmaz()
sm = sigmam()

options = Options()
c_op_list = []

omega0 = gmr * Bz * 2 * pi
omega1 = gmr * Bx * 2 * pi
tlist = np.linspace(0,pi/(omega1),10000)

NV = Processor(N=1, spline_kind="cubic" ) #T1 = 10**-3 T2 = 10**-6
NV.add_control(0.5*omega0*sz, label="Bz field")
NV.add_control(0.5*omega1*sx, label="Bx pulse")
NV.add_control(0.5*omega1*sy, label="By pulse")
NV.pulses[0].coeff = np.array([-1. for t in tlist])
NV.pulses[1].coeff = np.array([-np.cos(omega0*t) for t in tlist])
NV.pulses[2].coeff = np.array([np.sin(omega0*t) for t in tlist])
NV.pulses[0].tlist = tlist
NV.pulses[1].tlist = tlist
NV.pulses[2].tlist = tlist

#NV.plot_pulses(title="Pulse amplitudes")

basis0 = basis(2,0)
result = NV.run_state(basis0, e_ops = [sx, sy, sz], options=options)
expectation = [result.expect[0], result.expect[1], result.expect[2]]
""""print(expectation)"""

for i in range(0,9001,1000):
    plt.figure(dpi=800)
    axes = plt.gca()
    axes.set_ylim([-1.0,1.0])
    plt.plot(tlist[i:i+999],expectation[0][i:i+999])
    plt.plot(tlist[i:i+999],expectation[1][i:i+999])
    plt.plot(tlist[i:i+999],expectation[2][i:i+999])
    plt.xlabel("t")
    plt.ylabel(['Expectation values (omega1 =', str(omega1), ')'])
    plt.legend(["sx", "sy", "sz"], loc=4)

"""b = Bloch()
b.add_points(expectation)
b.show()"""