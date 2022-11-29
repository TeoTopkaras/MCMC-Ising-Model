import numpy as np
import matplotlib.pyplot as plt
from functions import *
from numpy import random
'''
MCMC simulation of an Ising model. The algorithm is able to generate and evolve an NxN lattice and measure
thermodynamic properties of the system, such as mean magnetisation and energy per spin, as well as, the
system's correlation time, magnetic susceptibility per spin and specific heat per spin. The algorithm allows
the simulation to generate lattices of different number of particles, while allowing the user to change the
system's coupling constant (J) and external field (H), by changing the relevant input parameters. The user is
able to determine the range of temperatures they wish to explore, as well as other aspects of the simulation,
by changing the relevant simulation constants. For a more detailed analysis of the functionality of each
constant please read the report. The simulation returns a list of values for all the output parameters
relating to each simulated system per temperature. The algorith allows for an independent call of the main
function from an external programm and the output format allows for the user to perform their prefered method
of statistical analysis. This program also performs bootstrap statistics, to account for cases with a small
size, and plots the thermodynamic behaviour of the system as a function of temperature via its plot functions.

'''
#-------------------------------------------------------------
#Input Parameters
#-------------------------------------------------------------
n_particles=50	#Number of particles in each side of the lattice 
coupl_const=1	#Value of the coupling constant of the system
ext_field=0	#Value of the external magnetic field
#-------------------------------------------------------------
#Simulation Constants
#-------------------------------------------------------------
temp_min=1.0		#Lower temperature limit
temp_max=4.0		#Upper temperature limit
temp_step=0.2		#Temperature increment
n_sys=10		#Number of systems simulated per temperature
n_training=100		#Number of systems to be evolved for the training algorithm
n_eq_lim=75		#Number of expected iterations for the system to be reaching equilibrium
n_bootstrap=1000	#Number of bootstrap iterations
#-------------------------------------------------------------
temp,corr_time,magn_spin,sigma_magn_spin,energy_spin,sigma_energy_spin,magn_sus,spec_heat=ising_mc_simu(n_particles,ext_field,coupl_const,temp_min,temp_max,temp_step,n_sys,n_training,n_eq_lim)

plot_corr_time(temp,corr_time,n_bootstrap)
plot_magn_spin(temp,magn_spin,sigma_magn_spin)
plot_energy_spin(temp,energy_spin,sigma_energy_spin)
plot_magn_sus(temp,magn_sus,n_bootstrap)
plot_spec_heat(temp,spec_heat,n_bootstrap)
