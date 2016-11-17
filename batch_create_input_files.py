""" 
Script to create a set of input files to the Oslo Shell Model Code
Made Nov 2016 by J{\o}rgen Eriksson Midtb{\o}
j.e.midtbo // at // fys.uio.no
github.com/jorgenem
"""
from __future__ import division
import numpy as np
import shell_model_utilities_v1 as smutil


# === Set up lists of nuclei with corresponding lists of parameters for SM calculation for each:
nucleus_list = [
# A   Z   MJ   par  totalJ	      single-particle energies     					 min_part            max_part      Nstates	 
[50,  22, 0,    1,     1,	[0.0,2.9447,4.487,7.2411,0,0,0,0], 				[0,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2], # NB! Only one type of valence particle, so input file generated is different and should be run with IdenticalParticles program
[52,  22, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [0,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2],#SPE energies should probably be adjusted! 
[54,  22, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [0,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2],
[58,  28, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [6,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2],
[60,  28, 0,    1,    -1,   [0.0,2.9447,4.487,7.2411,0.0,1.5423,2.754,3.7], [6,0,0,0,0,0,0,0], [8,2,0,0,4,6,2,4],    2]      
]


for nucleus in nucleus_list:
	# === Set parameters for current nucleus ===
	A = nucleus[0]
	Z = nucleus[1]
	MJ = nucleus[2]
	parity = nucleus[3] # 1=+, 0=-
	totalJ = nucleus[4] # Binary switch indicating even/odd total J. Only relevant for IdenticalParticles runs.
	SPE = nucleus[5] # Single particle energies ordered as [f7_p, p3_p, f5_p, p1_p, p3_n, f5_n, p1_n, g9_n] where f7 means f7/2 orbital, etc, and _p and _n signify proton and neutron, respectively.
	min_part = nucleus[6] # Minimum number of particles in each orbital. Same ordering as SPE.
	max_part = nucleus[7] # Maximum number of particles in each orbital. Same ordering as SPE.
	Nstates = nucleus[8] # Desired number of eigenstates

	input_filename = "input-A{:d}-Z{:d}.dat".format(A,Z)

	smutil.makeinputfile(A,Z,MJ,parity,totalJ,SPE,min_part,max_part,Nstates,input_filename=input_filename)