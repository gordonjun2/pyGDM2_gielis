# encoding: utf-8
import os, sys
import pickle
 
import time

import numpy as np
import matplotlib.pyplot as plt
#~ from mayavi import mlab


from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core
from pyGDM2 import linear





from pyGDM import visu
from pyGDM import tools
from pyGDM import main as mainOLD



## ---------- Setup structure
mesh = 'cube'
step = 20.0
geom1 = structures.rect_wire(step, L=10,H=3,W=4, mesh=mesh)
geom2 = structures.rect_wire(step, L=10,H=3,W=4, mesh=mesh)
geom2.T[1] += 80
#geom2.T[1] += 180
geometry = np.concatenate([geom1, geom2])


## -- block1: silicon, block2: gold
mat1 = len(geom1)*[materials.silicon()]
mat2 = len(geom2)*[materials.gold()]
material = mat1 + mat2


n1, n2 = 1.0, 1.0  # constant environment
struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))

print "N dipoles:", len(struct.geometry)




## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
kwargs = dict(theta = [0.0, 90.0])              # several polarizations
wavelengths = [500]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)





#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
EFIELD = core.scatter(sim, verbose=True)



    

#%%
#==============================================================================
# Compare Fields old/new
#==============================================================================
xm,ym,zm = geometry.T
print "Using:", EFIELD[0][0]
X,Y,Z = EFIELD[0][1].T
NFnew = np.array([xm,ym,zm, X.real,X.imag, Y.real,Y.imag, Z.real,Z.imag]).T



#==============================================================================
# Plot old/new
#==============================================================================
#visu.plotNF2Dreal(NFnew, tit='NEW, Efield, Re')

a = visu.animateNearfield(NFnew)
#a = visu.animateNearfield3D(NFnew, scale=20)
