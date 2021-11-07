#!/usr/bin/env python

# parameters to initialize the system
num_clust = 100        # num of Ptn clusters with n > 1
num_single_atom = 0  # num of single atom Pt
largest_cluster = 3  # largest initial cluster size in the system 

# total number of steps, write step, and temperature
MMAX   = 10000
wstep  =  1 
T      = 700

# unit cell parameters in A
unitcell_a = 13.216 
unitcell_b = 12.243
unitcell_c = 26.362

# primitive cell parameters in A (used for PES)
primcell_a = 6.607835831
primcell_b = 3.060749557
N_mesh = 11*11
xstep_max  =  primcell_a / 10.0
ystep_max  =  primcell_b / 10.0
numprimecellFactor = 2.0

# x and y boundaries of simulation
maxx   =  primcell_a * 17
minx   =  0.0
maxy   =  primcell_b * 17
miny   =  0.0

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
SinteringResultPlot = True # plot the result of sintering
LimitForOverlap = 1  # limit of total number of overlap process in one step. 
