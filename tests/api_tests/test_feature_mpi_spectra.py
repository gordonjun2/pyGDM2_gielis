#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import time

## --- for benchmark: limit openmp to 1 thread (for parallel scipy/numpy routines)
import os
nthreads = 1
os.environ["MKL_NUM_THREADS"] = "{}".format(int(nthreads))
os.environ["NUMEXPR_NUM_THREADS"] = "{}".format(int(nthreads))
os.environ["OMP_NUM_THREADS"] = "{}".format(int(nthreads))


from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core

import numpy as np


## --- It is not necessary to load mpi4py within the simulation script, this 
## --- will be done automatically by pyGDM2 prior to the actual MPI-simulation. 
## --- We do it however at this point to do some output to stdout only from 
## --- within the master process (rank == 0).
from mpi4py import MPI
rank = MPI.COMM_WORLD.rank



#==============================================================================
# Configure simulation
#==============================================================================
## ---------- Setup structure
mesh = 'cube'
step = 20.0
radius = 3.5
geometry = structures.sphere(step, R=radius, mesh=mesh)
material = materials.dummy(2.0)  # dummy material with constant refindex n=2.0

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization(mesh))


## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
wavelengths = np.linspace(400, 800, 20)   # spectrum
kwargs = dict(theta = [0.0])              # one polarizations
#kwargs = dict(theta = np.linspace(0.0, 90, 3), kSign=[-1, 1]) # several configurations
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)


## ---------- Simulation initialization
sim = core.simulation(struct, efield)



#==============================================================================
# Rum the simulation
#==============================================================================
## --- mpi
if rank == 0: 
    print "performing MPI parallel simulation... "
    t0_mpi = time.time()
core.scatter_mpi(sim, verbose=True, method='dyson', nthreads=nthreads)


## --- no mpi, for comparison
if rank ==0:
    t1_mpi = time.time()
    print '\n'*2
    print "performing sequential run for comparison...",
    sim2 = copy.deepcopy(sim)
    t0_seq = time.time()
    core.scatter(sim2, method='dyson', nthreads=nthreads)
    t1_seq = time.time()
    print " Done."






#%%
#==============================================================================
# Print results, compare timings
#==============================================================================
if rank == 0:
    print '\n'*6
    print "                -------- results -------- "
    print 
    print 20*'-'
    print "MPI"
    print 20*'-'
    print np.shape(sim.E)
    for i in sim.E:
        print ' config =', i[0], ',     N_dp=', len(i[1])
    print
    print
    
    print 20*'-'
    print "non-parallel execution"
    print 20*'-'
    print np.shape(sim2.E)
    for i in sim2.E:
        print ' config =', i[0], ',     N_dp=', len(i[1])
        
    print 
    print 
    print 20*'-'
    print "timing"
    print 20*'-'
    print "MPI-run:        ", round(t1_mpi-t0_mpi, 2), 's'
    print "sequantial-run: ", round(t1_seq-t0_seq, 2), 's'
    print "speed-up:        x{}".format(round((t1_seq-t0_seq)/(t1_mpi-t0_mpi), 1))
