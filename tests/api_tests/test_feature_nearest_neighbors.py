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
from pyGDM2 import tools
from pyGDM2 import visu








## ---------- Setup structure
n1, n2 = 1.0, 1.0  # constant environment
mesh = 'cube'
mesh = 'hex'
step = 10.0
geometry = structures.rect_wire(step, L=6,H=3,W=7, mesh=mesh, ORIENTATION=2)
geometry = structures.rect_dimer(step, L=13,H=3,W=4,G=2)
#geometry = structures.nanodisc(step, R=10,H=3, mesh=mesh)
#geometry = structures.hexagon(step, NSIDE=8,H=5, mesh=mesh)
#geometry = structures.split_ring(step, R=20,W=6,H=1,G=10, mesh=mesh)
#geometry = structures.prism(step, NSIDE=16,H=1, mesh=mesh)
#
## double splitring
#geometry = structures.split_ring(step, R=10,W=4,H=1,G=4, mesh=mesh)
#geometry.T[1] += 250
#geometry.T[0] += 250
#geometry = np.concatenate([structures.split_ring(step, R=10,W=4,H=1,G=4, mesh=mesh), geometry])

material = materials.dummy(2.0); mat_string_old='20'    # dummy material with constant and real dielectric function
material = materials.silicon(); mat_string_old='si'       # silicon dielectric function
material = materials.gold(); mat_string_old='au'          # gold dielectric function
#material = materials.alu(); mat_string_old='al'          # gold dielectric function



geo1 = structures.rect_wire(step, L=6,H=3,W=7, mesh=mesh, ORIENTATION=1)
geo2 = structures.rect_wire(step, L=6,H=3,W=7, mesh=mesh, ORIENTATION=2)

geo1 = structures.sphere(step, R=6, mesh=mesh, ORIENTATION=1)
geo2 = structures.sphere(step, R=6, mesh=mesh, ORIENTATION=2)


#%%
geo_test = structures.rect_dimer(step, L=5,H=1,W=2,G=0.5)
visu.structure(geo_test, scale=4)
print 'test geometry:'
print geo_test
#%%

#geo1 = structures.prism(step, NSIDE=9, H=1, mesh=mesh, ORIENTATION=1)
#geo2 = structures.prism(step, NSIDE=9, H=1, mesh=mesh, ORIENTATION=2)

struct1 = structures.struct(step, geo1, material, n1,n2, structures.get_normalization(mesh))
struct2 = structures.struct(step, geo2, material, n1,n2, structures.get_normalization(mesh))









## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
kwargs = dict(theta = [0.0])              # single polarizations
#kwargs = dict(theta = np.linspace(0,180, 90))              # several polarizations
wavelengths = [800]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)


## ---------- Simulation initialization
sim1 = core.simulation(struct1, efield)
sim2 = core.simulation(struct2, efield)
print "N dipoles:", len(sim1.struct.geometry)



#td1 = tools.get_geometry_2d_projection(geo1)
#td2 = tools.get_geometry_2d_projection(geo2)
td1 = tools.get_geometry_2d_projection(sim1)
td2 = tools.get_geometry_2d_projection(sim2)

SF1, SFv1 = tools.get_surface_meshpoints(td1, NN_bulk=3, max_bound=np.sqrt(2))
SF2, SFv2 = tools.get_surface_meshpoints(td2, NN_bulk=4, max_bound=np.sqrt(1.3))

visu.structure(sim1, color='g', scale=0.5, show=False)
try:
    visu.structure(SF1, color='r',  scale=1.5, show=False)
    plt.quiver(SF1.T[0], SF1.T[1], SFv1.T[0], SFv1.T[1])
except:
    pass

visu.structure_contour(sim1, color='b', style='lines', input_mesh='hex1',  show=0, borders=10)
#visu.structure_contour(sim1, color='b', style='lines', input_mesh='cube',  show=0, borders=10)
visu.structure_contour(sim1, color='b', style='dots', input_mesh='hex1',  show=0, borders=10)
plt.show()

visu.structure(sim2, color='g', scale=0.5, show=False)
try:
    visu.structure(SF2, color='r',  scale=1.5, show=False)
    plt.quiver(SF2.T[0], SF2.T[1], SFv2.T[0], SFv2.T[1])
except:
    pass


visu.structure_contour(sim2, color='b', style='lines', input_mesh='hex2',  show=0, borders=10)
#visu.structure_contour(sim2, color='b', style='lines', input_mesh='cube',  show=0, borders=10)
visu.structure_contour(sim2, color='b', style='dots', input_mesh='hex2',  show=0, borders=10)
plt.show()




