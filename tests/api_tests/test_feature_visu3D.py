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
from pyGDM2 import visu3d








## ---------- Setup structure
n1, n2 = 1.0, 1.0  # constant environment
mesh = 'cube'
#mesh = 'hex'
step = 5.0
#geometry = structures.rect_wire(step, L=6,H=3,W=7, mesh=mesh)
#geometry = structures.rect_dimer(step, L=5,H=3,W=4,G=1)
#geometry = structures.nanodisc(step, R=10,H=3, mesh=mesh)
#geometry = structures.hexagon(step, NSIDE=8,H=5, mesh=mesh)
#geometry = structures.split_ring(step, R=10,W=6,H=3,G=5, mesh=mesh)
geometry = structures.prism(step, NSIDE=32,H=3, mesh=mesh, ORIENTATION=1)
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


struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))








## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
kwargs = dict(theta = [0.0])              # single polarizations
#kwargs = dict(theta = np.linspace(0,180, 90))              # several polarizations
wavelengths = [800]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)


## ---------- Simulation initialization
sim = core.simulation(struct, efield)

print "N dipoles:", len(sim.struct.geometry)





#%%
core.scatter(sim, verbose=1)


#c = visu3d.structure(sim, scale=0.5)

NF = tools.get_field_as_list_by_fieldindex(sim, 0)
#c = visu3d.vectorfield(NF)
#c = visu3d.vectorfield_by_fieldindex(sim, 0)


#c = visu3d.vectorfield_color(NF)
c = visu3d.vectorfield_color_by_fieldindex(sim, 0)


#%% --- animate
#NF = tools.get_field_as_list_by_fieldindex(sim, 0)
#ani = visu3d.animate_vectorfield(NF, Nframes=50, draw_struct=1, save_anim=1, figsize=(1200, 800))


