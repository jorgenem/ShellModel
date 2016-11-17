"""
Script to batch run shell model calculations
Made November 2016 by J{\o}rgen Eriksson Midtb{\o}
j.e.midtbo//at//fys.uio.no
github.com/jorgenem
"""
from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 
import subprocess
import os
import shutil
import sys
from shell_model_utilities_v1 import *



# === Set up lists of nuclei with corresponding lists of parameters for SM calculation for each:
nucleus_list = [
# A   Z   MJ   par  totalJ	      single-particle energies     					 min_part            max_part      Nstates	  
# TODO: Make it work for only one type of valence particle
# [50, 22, 0,    1,      0,	 [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [6,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2] ]
[52,  22, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [0,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2]      ,#SPE energies should probably be adjusted! 
[54,  22, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [0,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2]      ,
[58,  28, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [6,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2]      ,
[60,  28, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [6,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2]      
]

# Experimental values - must match with nucleus ordering and Nstates above
exp_values = [
# energies  	J      parity (1=+,0=-)
[[0, 1.05006],  [0,2], [1,1]], # 52,22 Ti
[[0, 1.4948],  [0,2], [1,1]],  # 54,22 Ti
[[0, 1.45421],  [0,2], [1,1]], # 58,28 Ni
[[0, 1.332514], [0,2], [1,1]]  # 60,28 Ni
]

# Set run mode: 0 for calculation of levels, 1 for comparison of previous calculations to experimental levels
mode = 1


if mode == 0:
	for nucleus in nucleus_list:
		# === Set parameters for current nucleus ===
		A = nucleus[0]
		Z = nucleus[1]
		MJ = nucleus[2]
		parity = nucleus[3] # 1=+, 0=-
		totalJ = nucleus[4]
		SPE = nucleus[5] # Single particle energies ordered as [f7_p, p3_p, f5_p, p1_p, p3_n, f5_n, p1_n, g9_n] where f7 means f7/2 orbital, etc, and _p and _n signify proton and neutron, respectively.
		min_part = nucleus[6] # Minimum number of particles in each orbital. Same ordering as SPE.
		max_part = nucleus[7] # Maximum number of particles in each orbital. Same ordering as SPE.
		Nstates = nucleus[8] # Desired number of eigenstates
	
		
		# === Start calculation ===
		# Copy current TBMEs and current inputfile to rundir
		if not os.path.exists('rundir'):
		    os.makedirs('rundir')

		# Set path to executable depending on whether it is identical particles or p&n:
		if Z == 20 or A-Z == 28:
			path_to_PAR_lanczo = '/home/jorgenem/gitrepos/CENS-fork/FCI/parallel/IdenticalParticles/src/PAR-lanczo'
		else:
			path_to_PAR_lanczo = '/home/jorgenem/gitrepos/CENS-fork/FCI/parallel/pnCase/src/PAR-lanczo'
		makeinputfile(A,Z,MJ,parity,totalJ,SPE,min_part,max_part,Nstates)

		# sys.exit(0)

		makeTBMEfile(A)
		os.chdir('rundir')
		# TODO: Figure out a way to use MPI in this. Turns out to be complicated. Consider asking Anders Kvellestad for advice. The below idea does not work.
		# Nprocs = 4 # Number of CPU cores to use for mpirun -np Nprocs call.
		# print 'mpirun -np {:d} {:s}'.format(Nprocs,path_to_PAR_lanczo)
		# sys.exit(0)
		subprocess.call([path_to_PAR_lanczo, "input.dat"])
		os.chdir('../')
		if not os.path.exists('resultsdir'):
			os.makedirs('resultsdir')
		shutil.copy('rundir/outputRANK0-result.dat','resultsdir/result-A{:d}-Z{:d}.dat'.format(A,Z))
		shutil.copy('rundir/outputRANK0-eigen-vectors.dat1','resultsdir/eigen_vectors-A{:d}-Z{:d}.dat1'.format(A,Z))


elif mode == 1:
	energies_list = []
	nucleus_counter = 1
	for nucleus in nucleus_list:
		filename = "resultsdir/result-A{:d}-Z{:d}.dat".format(nucleus[0], nucleus[1])
		E_gs, energies = read_energy_levels(filename)
		print "Nucleus #{:d}. A={:d}, Z={:d}. GS energy {:10.6f}. Excitation energies: \nlevel      energy(calc)    energy(exp)       <J^2>(calc)   J^2(exp)".format(nucleus_counter, nucleus[0], nucleus[1], E_gs)
		# print energies
		for i in range(len(energies)):
			print "{:3d}    {:12.6f}     {:12.6f}         {:3.0f}           {:3.0f}".format(i,energies[i][0],exp_values[nucleus_counter-1][0][i],np.abs(energies[i][1]),exp_values[nucleus_counter-1][1][i]*(exp_values[nucleus_counter-1][1][i]+1))

		energies_list.append(energies)	
		nucleus_counter += 1