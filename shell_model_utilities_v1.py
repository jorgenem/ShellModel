"""
This library contains a set of functions to facilitate running the 
Oslo Shell-Model code of M. Hjorth-Jensen et. al. , which can be found at
https://github.com/ManyBodyPhysics/CENS


Made November 2016 by J{\o}rgen Eriksson Midtb{\o}
University of Oslo
Postboks 1048 Blindern
0316 OSLO
Norway

j.e.midtbo//at//fys.uio.no
github.com/jorgenem


Copyright (C) 2016  Joergen Eriksson Midtboe

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
from __future__ import division
import numpy as np 
import os

def makeTBMEfile(A,nn_orig=0,pn_orig=0,pp_orig=0):
	"""
	Function to create three files of two-body matrix elements for shell-model
	calculations. It is assumed that they correspond to a Ca48 core, and they
	are scaled by (A/48)^-0.3 relative to their original values, as prescribed
	in Honma et. al., Nucl. Phys. A (2002).
	"""
	# Read original (master) TBME files
	if nn_orig==0 or pn_orig==0 or pp_orig==0:
		nn_orig = np.genfromtxt("f5pg9shell_nn.dat", skip_header=1)
		pn_orig = np.genfromtxt("fpg9shell_pn.dat", skip_header=1)
		pp_orig = np.genfromtxt("fpshell_pp.dat", skip_header=1)

	# Copy into new variables and update the TBME values column
	nn_current = nn_orig
	pn_current = pn_orig
	pp_current = pp_orig
	nn_current[:,14] = np.power(A/48, -0.3) * nn_orig[:,14]
	pn_current[:,14] = np.power(A/48, -0.3) * pn_orig[:,14]
	pp_current[:,14] = np.power(A/48, -0.3) * pp_orig[:,14]

	# Write to files
	try:
		os.remove("nn_current.dat")
	except OSError:
		pass
	nn_file = open("nn_current.dat", "w")
	nn_file.write("%d\n"%nn_current.shape[0])
	for i in range(nn_current.shape[0]):
		nn_file.write("%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%14.6f\n"%(nn_current[i,0],nn_current[i,1],nn_current[i,2],nn_current[i,3],nn_current[i,4],nn_current[i,5],nn_current[i,6],nn_current[i,7],nn_current[i,8],nn_current[i,9],nn_current[i,10],nn_current[i,11],nn_current[i,12],nn_current[i,13],nn_current[i,14]))
	nn_file.close()
	try:
		os.remove("pn_current.dat")
	except OSError:
		pass
	pn_file = open("pn_current.dat", "w")
	pn_file.write("%d\n"%pn_current.shape[0])
	for i in range(pn_current.shape[0]):
		pn_file.write("%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%14.6f\n"%(pn_current[i,0],pn_current[i,1],pn_current[i,2],pn_current[i,3],pn_current[i,4],pn_current[i,5],pn_current[i,6],pn_current[i,7],pn_current[i,8],pn_current[i,9],pn_current[i,10],pn_current[i,11],pn_current[i,12],pn_current[i,13],pn_current[i,14]))
	pn_file.close()
	try:
		os.remove("pp_current.dat")
	except OSError:
		pass
	pp_file = open("pp_current.dat", "w")
	pp_file.write("%d\n"%pp_current.shape[0])
	for i in range(pp_current.shape[0]):
		pp_file.write("%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%14.6f\n"%(pp_current[i,0],pp_current[i,1],pp_current[i,2],pp_current[i,3],pp_current[i,4],pp_current[i,5],pp_current[i,6],pp_current[i,7],pp_current[i,8],pp_current[i,9],pp_current[i,10],pp_current[i,11],pp_current[i,12],pp_current[i,13],pp_current[i,14]))
	pp_file.close()






def makeinputfile(A,Z,MJ,parity,totalJ,SPE,min_part,max_part,Nstates,max_iterations=1000,mem_size=1000,file_size=10,input_filename="input.dat"):
	""" 
	Function to make input files to the shell-model code.

	Explanation of variables:
	A = nucleon number
	Z = proton number
	MJ = twice total projection of angular momentum
	totalJ = only applies in cases where only protons OR neutrons are in valence orbitals, is a binary switch 1/0 corresponding to even/odd total J of nucleus
	parity = 1 for +, 0 for -
	 Single particle energies ordered as [f7_p, p3_p, f5_p, p1_p, p3_n, f5_n, p1_n, g9_n] where f7 means f7/2 orbital, etc, and _p and _n signify proton and neutron, respectively.
	SPE = single particle energies, ordered as [f7_p, p3_p, f5_p, p1_p, p3_n, f5_n, p1_n, g9_n] where f7 means f7/2 orbital, etc, and _p and _n signify proton and neutron, respectively.
	min_part = minimum number of particles in each orbital. Same ordering as SPE.
	max_part = maximum number of particles in each orbital. Same ordering as SPE.
	Nstates = desired number of eigenstates to calculate
	max_iterations = maximum number of allowed Lanczo iterations, recommended by code authors to be <1000
	mem_size = maximum memory to store nondiagonal <pn_SD'|OP|pn_SD> in MB
	file_size = maximum file size of lanczo eigenvectors in GB
	"""



	# Run consistency checks:
	if A <= 48:
		raise ValueError("Error in function makeinputfile: A<=48 is not supported, there must be at least one particle on top of the Ca48 core.")
	if not parity in [0,1]:
		raise ValueError("Error in function makeinputfile: Parity value {0} is not allowed.".format(parity))
	if Z == 20 or A-Z == 28:
		if not totalJ in [0,1]:
			raise ValueError("Error in function makeinputfile: Seems to be a case of only proton or neutron excitations, but totalJ is {:d}, should be 1 (even) or 0 (odd)".format(totalJ))


	# Check if this is a case with only proton or neutron excitations, or both:
	if Z == 20 or A-Z == 28:
		# This is a case of only proton or neutron exicitations. Will make the input file and run the code for identical particles.
		if Z == 20:
			case = "n" # identifies this as a case of neutrons only
			particle_number = A-Z-28 # number of valence nucleons
		else:
			case = "p" # case of protons only
			particle_number = Z-20

		# Add header 
		input_string = """
		   ** All input data are identified through a text-line which
           ** starts with the two symbols '/' and '*' and ends with ":" or ";"
           ** Additional informative text may be included in the form
           **        <........ Informative text....... >
           ** 


\* Shell model calculation of identical particle system -- title: output \n"""

			
		# Specify particle number
		input_string += "\* The particle number: {:d}\n\n".format(particle_number)
		# Specify total angular momentum, projection and parity 
		input_string += "\* Total angular momentum J is (even, odd): {:s}\n".format("even" if totalJ else "odd") # The variable MJ for identical particles assumes role of binary switch for even/odd total J
		input_string += "\* Twice total projection of angular momentum: {:d}\n".format(MJ)
		input_string += "\* Total parity (+, -): {}\n\n".format("+" if parity else "-")
		# Set up orbital configuration
		input_string += "\* The number of particle j-orbits: 4\n\n" # Four orbitals for both p and n
		input_string += "	       <n  l  2*j   min_part  max_part energy:\n"
		if case == "p":
			input_string += "\* Orbit:   1  1   1       {:d}         {:d}     {:f}\n".format(min_part[3],max_part[3],SPE[3])
			input_string += "\* Orbit:   0  3   5       {:d}         {:d}     {:f}\n".format(min_part[2],max_part[2],SPE[2])
			input_string += "\* Orbit:   1  1   3       {:d}         {:d}     {:f}\n".format(min_part[1],max_part[1],SPE[1])
			input_string += "\* Orbit:   0  3   7       {:d}         {:d}     {:f}\n\n".format(min_part[0],max_part[0],SPE[0])
		elif case == "n":
			input_string += "\* Orbit:   0  4   9       {:d}         {:d}     {:f}\n".format(min_part[7],max_part[7],SPE[7])	
			input_string += "\* Orbit:   1  1   1       {:d}         {:d}     {:f}\n".format(min_part[6],max_part[6],SPE[6])
			input_string += "\* Orbit:   0  3   5       {:d}         {:d}     {:f}\n".format(min_part[5],max_part[5],SPE[5])
			input_string += "\* Orbit:   1  1   3       {:d}         {:d}     {:f}\n".format(min_part[4],max_part[4],SPE[4])
		else: 
			raise ValueError("Error in function makeinputfile: Seems to be a case of only one type of valence particle, but neither protons nor neutrons. New Physics? (No, it has to be a bug.)")

		# Specify type of calculation (assume this is held fixed)
		input_string += """<***** Type of interaction ******  value ****>
<** Two-part   veffJ interaction:             1 **>
<** Two-part   veffM interaction:             2 **>
<** Converted Two-part pluss three-part elem  3 **>

\* Type of effective interaction: 1\n"""

		# Specify TBME file 
		if case == "p":
			input_string += "\* Input V-effective matrix element: pp_current.dat\n\n\n"
		elif case == "n":
			input_string += "\* Input V-effective matrix element: nn_current.dat\n\n\n"
		else: 
			raise ValueError("Error in function makeinputfile: Seems to be a case of only one type of valence particle, but neither protons nor neutrons. New Physics? (Just kidding, it has to be a bug.)")

		# Add more compulsory info which we will not alter:
		input_string += """\* Center_of_Mass matrix elements included(yes = 1, no = 0): 0
\* If(yes): Filename Center_of_Mass matrix elements: 

<**************** Type of calculation ****************     Symbol   ********>
    <** Calculation to find the shell model dimension    --> dimension	  **>
    <** Lanczos iteration based on random initial vector --> random-start **>

\* Type of calculation process: random-start


<** Temporary files for storage of Lanczos vector **>
"""
		
		# Specify max file size for lanczos vectors, mem size, max iterations and no of eigenstates
		input_string += "\* Max. file size for lanczos vectors in Gb    (default 10): {:d}\n".format(file_size)
		input_string += "\* Temporary memory to store nondiag <SD'|OP|SD> in MB (default 100 MB): {:d}\n\n\n".format(mem_size)
		input_string += "<** Energy eigenvalue parameters **>\n\n"
		input_string += "\* Maximum Lanczos iterations (<1000):  {:d}\n".format(max_iterations)
		input_string += "\* Wanted number of converged eigenstates: {:d}\n\n".format(Nstates)
		input_string += "<** END_OF_INPUT_FILE **>"










	else:
		# This is a case of both proton and neutron excitations. Will make the input file and run the code for proton&neutron case.
		# Start filling input file:
		
		# Add first nine lines, always the same:
		input_string = """			** All input data are identified through a text-line which
	        ** starts with the two symbols '/' and '*' and ends with ":" or ";"
	        ** Additional informative text may be included in the form
	        **        <........ Informative text....... >
	        ** 
	
\* Shell model calculation of proton/neutron systemA -- title: output\n\n\n"""
		# Add proton and neutron number:
		input_string += "\* The proton number:   {0}\n".format(Z-20)
		input_string += "\* The neutron number:  {0}\n\n".format(A-Z-28)
		
		# Add angular momentum and parity info:
		input_string += "\* Twice total projection of angular momentum: {0}\n".format(MJ)
		input_string += "\* Total parity (+, -): {}\n\n".format("+" if parity else "-")
		
		# Add orbital info:
		input_string += "\* The number of proton particle j-orbits: 4\n"
		input_string += "             <N  n  l  2*j   min_part  max_part energy\n"
		input_string += "\* Orbit_Z:   3  1  1   1       {:d}         {:d}     {:f}\n".format(min_part[3],max_part[3],SPE[3])
		input_string += "\* Orbit_Z:   3  0  3   5       {:d}         {:d}     {:f}\n".format(min_part[2],max_part[2],SPE[2])
		input_string += "\* Orbit_Z:   3  1  1   3       {:d}         {:d}     {:f}\n".format(min_part[1],max_part[1],SPE[1])
		input_string += "\* Orbit_Z:   3  0  3   7       {:d}         {:d}     {:f}\n\n".format(min_part[0],max_part[0],SPE[0])
		input_string += "\* The number of neutron particle j-orbits: 4\n"
		input_string += "             <N  n  l  2*j   min_part  max_part energy\n"
		input_string += "\* Orbit_N:   4  0  4   9       {:d}         {:d}     {:f}\n".format(min_part[7],max_part[7],SPE[7])	
		input_string += "\* Orbit_N:   3  1  1   1       {:d}         {:d}     {:f}\n".format(min_part[6],max_part[6],SPE[6])
		input_string += "\* Orbit_N:   3  0  3   5       {:d}         {:d}     {:f}\n".format(min_part[5],max_part[5],SPE[5])
		input_string += "\* Orbit_N:   3  1  1   3       {:d}         {:d}     {:f}\n".format(min_part[4],max_part[4],SPE[4])
	
		# Add stuff about type of program process:
		input_string += """
	
<************* type of program process ***********>

 <** Effective interaction Lanczo process     0  **>
 <** Calculation of shell-model dimension     1	 **>


\* Type of calculation process: 0

      <***** Type of interaction ******>

    <** Two-part veffJ    0 **>
    <** Two-part veffM    1 **>
    <** Three-part veffM  2 **>

\* Type of effective interaction: 0

\* Input proton-proton   V-effective:    pp_current.dat
\* Input neutron-neutron V-effective:    nn_current.dat
\* Input proton-neutron  V-effective:    pn_current.dat


<** Temp. files for storage of nondiag_matr_elem and Lanczos vector **>
	
	"""
	
		# Set memory size and max file size:
		input_string += "\* Memory size to store nondiag <pn_SD'|OP|pn_SD> in MB (default 500 MB):{:d}\n\n".format(mem_size)
		input_string += "\* Max file size for lanczos vectors    in Gb (default 10 GB): {:d}\n\n".format(file_size)
	
		# Set wanted number of states and maximum number of iterations:
		input_string += "<** Energy eigenvalue parameters **>\n\n"
		input_string += "\* Maximum Lanczos iterations (<1000): {:d}\n".format(max_iterations)
		input_string += "\* Wanted number of converged eigenstates: {:d}\n\n".format(Nstates)
	
		# Puh! Finally just add end of file statement:
		input_string += "<*** END_OF_INPUT_FILE ***> "
	
	
	# END if test for only one type of particle or both
	
	# Write to file
	try:
		os.remove(input_filename)
	except OSError:
		pass
	current_file = open(input_filename,'w')
	current_file.write(input_string)
	current_file.close()





def read_energy_levels(filename):
	"""
	Function to read calculated energy levels from shell model output file
	"""

	file_tmp = open(filename, 'r')
	data = file_tmp.readlines()
	file_tmp.close()

	# Loop over all lines
	E_gs = -999999 # Ground state energy (allocate dummy value)
	i_line = 0 # Current line number
	for line in data:
		words = line.split()
		# Get ground state energy
		if len(words) > 3:
			if words[1] == "ground" and words[2] == "state" and words[3] == "energy":
				E_gs = float(words[5])
				i_line_gs = i_line # Store the line number at which the ground state energy is given. This is used to count down to the lines where the other energy levels are.

				# Check what type of file this is: identical particles or proton-neutron
				identicalparticles = True
				if data[i_line_gs+3].split()[0] == "Proton": # This is the signature of a file containing p&n results
					identicalparticles = False

				# Stop here, get the other energy levels manually
				break 

		if i_line == len(data)-1:
			raise Exception("Error in function read_energy_levels: Reached end of results file without finding levels. Did the run succeed?")

		i_line += 1

	# print "E_gs =", E_gs

	# Calculate the number of energy levels present in the file:
	if identicalparticles:
		N_levels = (len(data) - i_line_gs - 4)/3 # This seems to match the formatting of the output file such that N_levels becomes a whole number
		print N_levels
		levels = np.zeros((int(N_levels), 2)) # Matrix to store level info (energy relative to gs, J**2)
		for i_level in range(int(N_levels)):
			words = data[i_line_gs + 1 + 3*i_level + 1].split()
			levels[i_level, :] = (float(words[1]), float(words[4]))
	else:
		N_levels = (len(data) - i_line_gs - 4)/6 # This seems to match the formatting of the output file such that N_levels becomes a whole number
		levels = np.zeros((int(N_levels), 2)) # Matrix to store level info (energy relative to gs, J**2)
		for i_level in range(int(N_levels)):
			words = data[i_line_gs + 1 + 6*i_level + 1].split()
			levels[i_level, :] = (float(words[1]), float(words[4]))
	return E_gs, levels




def plot_energy_levels(filename):
	"""
	Function to visualize calculated energy levels by plotting them as horizontal lines above each other
	"""
	Levels = read_energy_levels(filename)
	# Plot energy levels
	plt.figure(0)
	plt.hold('on')
	plt.xlim([0,1])
	plt.ylim([-0.1*max(Levels[:,0]), 1.1*max(Levels[:,0])])
	for i_level in range(len(Levels[:,0])):
		plt.plot([0.2,0.6], [Levels[i_level,0], Levels[i_level,0]], color='k', lw=1.5)
	plt.text(0.7, Levels[-5,0], '$^{95}\mathrm{Zr}$, \n$N_\mathrm{levels} =%d$' %N_levels, fontsize=15)
	plt.show()