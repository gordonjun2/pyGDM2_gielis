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
from pyGDM2 import visu

from pyGDM2 import tools
from pyGDM2.tools import get_geometry, get_step_from_geometry, get_geometry_2d_projection





#from pyGDM import visu
#from pyGDM import tools
#from pyGDM import main as mainOLD



## ---------- Setup structure
mesh = 'cube'
mesh = 'hex'
step = 10
geometry = structures.rect_wire(step, L=50,H=3,W=3, mesh=mesh)
geometry = structures.hexagon(step, NSIDE=10, H=1, mesh=mesh)
geometry = structures.prism(step, NSIDE=40, H=1, mesh=mesh)



#material = materials.dummy(2.0)
material = materials.silicon()
material = materials.gold()
material = materials.alu()

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))
struct = structures.center_struct(struct)



## ---------- Setup incident field
field_generator = fields.focused_planewave
#field_generator = fields.gaussian
xSpot = np.linspace(-200, 200, 51)
ySpot = np.linspace(-200, 200, 51)
xSpot = [0]
ySpot = [0]
kwargs = dict(theta = 0.0, spotsize=[250], kSign=-1, xSpot=xSpot, ySpot=ySpot)
wavelengths = [700]                     # one single wavelength


#field_generator = fields.planewave
#kwargs = dict(theta = 90.0)
#wavelengths = [600]                     # one single wavelength


efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)


## ---------- Simulation initialization
sim = core.simulation(struct, efield)




#remap = adapt_map_to_structure_mesh(sim)
MAP = tools.generate_NF_map(-250,250,81, -250,250,81, sim.struct.geometry.T[2].min())
remap = tools.adapt_map_to_structure_mesh(MAP, sim, projection='XY', 
                                          min_dist=1.5, verbose=True)
MAP2 = tools.generate_NF_map(-250,250,81, -250,250,81, sim.struct.geometry.T[2].max())
remap2 = tools.adapt_map_to_structure_mesh(MAP2, sim, projection='XY', 
                                          min_dist=1.5, verbose=True)

visu.structure(sim.struct.geometry, projection='xz', scale=0.5)
visu.structure(sim.struct.geometry, scale=0.5)
visu.structure(remap.T, scale=0.25)
print "N dipoles:", len(sim.struct.geometry)


#%%


#rasterscan_fieldconfigs = tools.get_possible_field_params_rasterscan(sim)
#
#print 'available rasterscan configurations' 
#for p in rasterscan_fieldconfigs:
#    print p



#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
core.scatter(sim, method='LU', verbose=True)


##%%

#E_rasterscan = tools.get_rasterscan_fields(sim, 0)

#posRS, sigRS = tools.calculate_rasterscan(sim, 0, linear.nearfield, r_probe=[0,0,200])

Es1, Et1, Bs, Bt = linear.nearfield(sim, 0, remap)
Es2, Et2, Bs, Bt = linear.nearfield(sim, 0, remap2)

plt.subplot(aspect='equal')
visu.vectorfield_color(Et1, show=False)
visu.vectorfield_color(Et2, show=False)
#visu.structure(sim.struct.geometry, color='w', scale=0.3, show=False)
#visu.structure(remap.T, scale=0.3, show=False)


plt.show()


