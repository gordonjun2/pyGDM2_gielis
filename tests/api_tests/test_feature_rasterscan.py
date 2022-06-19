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
from pyGDM2 import nonlinear
from pyGDM2 import visu
from pyGDM2 import tools





#from pyGDM import visu
#from pyGDM import tools
#from pyGDM import main as mainOLD



## ---------- Setup structure
mesh = 'cube'
#mesh = 'hex'
step = 20.0
geometry = structures.rect_wire(step, L=50,H=3,W=3, mesh=mesh)
geometry = structures.hexagon(step, NSIDE=12, H=1, mesh=mesh)
material = materials.dummy(2.0); mat_string_old='20'    # dummy material with constant and real dielectric function
material = materials.silicon(); mat_string_old='si'       # silicon dielectric function
material = materials.gold(); mat_string_old='au'          # gold dielectric function

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))
struct = structures.center_struct(struct)



## ---------- Setup incident field
field_generator = fields.focused_planewave        # planwave excitation
field_generator = fields.gaussian        # planwave excitation
#kwargs = dict(theta = [0.0, 45.0, 90.0], kSign=[-1,1])              # several configurations
#wavelengths = np.linspace(400,1000,25)                     # spectrum
xSpot = np.linspace(-400, 400, 31)
ySpot = np.linspace(-400, 400, 31)
kwargs = dict(theta = 0.0, spotsize=[200], kSign=-1, xSpot=xSpot, ySpot=ySpot)
wavelengths = [800]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)

visu.structure(sim.struct.geometry, scale=0.5)
print "N dipoles:", len(sim.struct.geometry)



rasterscan_fieldconfigs = tools.get_possible_field_params_rasterscan(sim)

print 'available rasterscan configurations' 
for p in rasterscan_fieldconfigs:
    print p



#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
core.scatter(sim, method='LU', verbose=True)


#%%


TPL = tools.calculate_rasterscan(sim, 0, nonlinear.tpl_ldos, nonlin_order=2)


visu.scalarfield(TPL, tit='TPL signal', cmap='viridis')

