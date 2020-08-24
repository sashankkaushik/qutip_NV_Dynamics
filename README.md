# qutip_NV_Dynamics
Python code using QuTip to simulate NV dynamics in spin magnetic resonance, along with pulses, non-linear dynamics and decoherence.

Time_Dep_H2.py simulates the Time-dependent Rabi Hamiltonian of a spin-1/2 system with qutip.mesolve().
NV_spin_half.py simulates Magnetic Resonance dynamics for a spin-1/2 system, with corresponding T1 and T2 decay times.
NV_spin_one.py simulates the NV ground state Spin Hamiltonian with zfs and changeable tuning of transitions based on resonant frequencies. T1 and T2 decay functionalities haven't been added yet.
Floquet.py simulates quantum chaotic dynamics via the Kicked Top Hamiltonian implemented with discrete Floquet propagators and generates the QFI for different evolution periods.
