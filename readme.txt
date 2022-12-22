Tourbillon performs simulations of polarisation transfer among protons using
low-order correlations in Liouville space. 

This is a modified version of the code originally distributed as supplementary 
information for "Ab initio simulation of proton spin diffusion" by J.N. Dumez, 
M.C. Butler, E. Salager, B. Elena and L. Emsley. 
Please cite this paper if you use this scheme or any part of it.

This code has been modified to add restricted basis sets, to use unordered maps
to store the density matrices, and to include chemical shift evolution.

--

Compilation:

To compile the program, simply indicate a c++ compiler in the makefile; only
the standard c++ library is used.

--

Usage: 

ssd <name>
The input file <name>.ssdi is read, the simulation is performed, and two files
are created, <name>.ssdo and <name>.dat

--

Files:

.ssdi file : contains all the parameters for the simulation. Detailed
explanations for new parameters can be found in the example input files 
"FUMTEM.ssdi" and "FUMTEM_TP.ssdi", while explanations for parameters
distributed with the original version may be found in the SI for that paper.

.ssdo file : contains a summary of the processed parameters.

.dat file : contains the polarisation curves for all the protons in the
simulation. The first column contains the time (in ms) and the subsequent
columns the quantity < I_{iz} > / || I_{1z} ||^2 for every proton i in the
system. 
The file is named according to the following convention:
<name>_<rho0>.dat
where <rho0> is the indices of the nuclei initially polarised

--

The average number of spin neighbours per spin may be changed by adjusting
#define NMAX
in crystal.cpp
