#!/usr/bin/env python
import numpy as np

# param.py for Al2O3
# parameters to initialize the system
num_clust = 100        # num of Ptn clusters with n > 1
num_single_atom = 0  # num of single atom Pt
largest_cluster = 6  # largest initial cluster size in the system 

# total number of steps, write step, and temperature
MMAX   = 100000
wstep  =  1 
T      = 750

# unit cell parameters in A
unitcell_a = 13.216 
unitcell_b = 12.243
unitcell_c = 26.362

# primitive cell parameters in A (used for PES)
primcell_a = 4.84025 #these are the lengths of the a and b vectors- not the coordinates of then in cartesian
primcell_b = 4.84025 #would it be a good idea to add in the angle between the unit cells as well, and make that  something used to walk around on the mesh?
N_mesh = 11*11
xstep_max  =  primcell_a / 10
ystep_max  =  primcell_b / 10
numprimecellFactor = 5.0

# x and y boundaries of simulation
fixFactorForY = np.sqrt(3)/2 
num_primcellInOneDirection = 25 # this have to be an integer
maxx   = num_primcellInOneDirection * primcell_a
minx   =  0.0 #this needs to be defined differently - I think if i define the vector rel to (0,0) then scale the vector appropriately- probably just define it in the cdes
maxy   =  primcell_b * fixFactorForY  * num_primcellInOneDirection
miny   =  0.0 * fixFactorForY

# constants
hartree_to_ev = 27.2114
k_J           = 1.38064852e-23
k_ev          = 8.6173303e-5
k_hartree     = k_ev/hartree_to_ev
kT            = k_hartree*T
beta          = 1.0/kT

# VDW radius of Pt in A
Ratom  =  1.75

# Control
CounterLimit = 50000
SinteringResultPlot = False  # plot the result of sintering
LimitForOverlap = 5  # limit of total number of overlap process in one step. 
