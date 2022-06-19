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
step = 10.0
geometry1 = structures.rect_wire(step, L=26,H=1,W=7, mesh=mesh)
#geometry = structures.rectDimer(step, L=13,H=3,W=4,G=2)
geometry2 = structures.nanodisc(step, R=10,H=3, mesh=mesh)
#geometry = structures.hexagon(step, NSIDE=8,H=5, mesh=mesh)
#geometry = structures.split_ring(step, R=20,W=6,H=1,G=10, mesh=mesh)
#geometry = structures.prism(step, NSIDE=16,H=1, mesh=mesh)

geometry2.T[0] -= 50.0  # shift X
geometry2.T[1] += 150.0  # shift Y

geometry = np.concatenate([geometry1, geometry2])

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

## -- structure geometry, 2D-projections
s = visu.structure(sim, scale=0.75, show=0)
#c = visu.structure_contour(sim, color='r',style='dots', dashes=[2,2], show=1)
c = visu.structure_contour(sim, color='r', style='dots', show=1)







#%%
########################################################################
##############            Test Simulations               ###############
########################################################################

print "\n\nSIM ----"
EFIELD = core.scatter(sim, method='LU', verbose=True)



    
    

#%%
#==============================================================================
# test visulalization routines
#==============================================================================





## -- 3 ways to plot the nearfield inside the structure:
plt.figure(figsize=(12,3))
plt.subplot(131)
v = visu.vectorfield(sim.E[0], sim, complex_part='real', show=0)
plt.subplot(132)
v = visu.vectorfield_by_fieldindex(sim, 0, complex_part='real', show=0)

NF = tools.get_field_as_list(sim.E[0], sim)
NF = tools.get_field_as_list_by_fieldindex(sim, 0)
plt.subplot(133)
v = visu.vectorfield(NF, complex_part='real', show=0)
plt.show()







## -- nearfield maps, streamplot
map_def = tools.generate_NF_map(-400,400,51, -400,400,51, Z0=sim.struct.geometry.T[2].max()+3*step)
nf_map_Es, nf_map_E1, nf_map_Bs, nf_map_B1 = linear.nearfield(sim, 0, map_def)

plt.figure(figsize=(13,4))
plt.subplot(131)
plt.title("Real part, E-Field")
c = visu.structure_contour(sim, color='.5', show=0, zorder=5)
visu.vectorfield(nf_map_Es, show=0)

plt.subplot(132)
plt.title("Streamlines, E-Field")
c = visu.structure_contour(sim, color='.5', show=0, zorder=5)
visu.vectorfield_fieldlines(nf_map_Es, show=0)

plt.subplot(133)
plt.title("Intensity of E-Field")
c = visu.structure_contour(sim, color='w', show=0, zorder=5)
im = visu.vectorfield_color(nf_map_E1, fieldComp='I', show=0)
plt.colorbar()

plt.show()
#%%





## -- farfield back-scattering BFP image
tetalist, philist, I_sc, I_tot, I0 = linear.farfield(
                                              sim, field_index=0,
                                              r=10000, tetamin=0, tetamax=np.pi/2., 
                                              Nteta=20, Nphi=36)

c= visu.farfield_pattern_2D(tetalist, philist, I_sc, show=False)
plt.colorbar()
plt.show()
#%%
q = linear.heat(sim, field_index=0, return_value='structure')

im = visu.scalarfield(q, tit='heat distribution', show=False)
plt.colorbar(label=r"heat (nW/nm^3)")
plt.show()


#%% --- animate
import matplotlib

NF = tools.get_field_as_list_by_fieldindex(sim, 0)
ani = visu.animate_vectorfield(NF, show=False)
ani.save('video.mp4', writer=matplotlib.animation.FFMpegWriter(fps=25))   # save video to file
plt.show()
