# encoding: utf-8
import os, sys
 
import time

import numpy as np
import matplotlib.pyplot as plt
#~ from mayavi import mlab


from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core
from pyGDM2 import tools
from pyGDM2 import linear
from pyGDM2 import visu





#==============================================================================
# Testscript
#==============================================================================
## ---------- Setup structure
mesh = 'cube'
step = 20.0
radius = 3.
geometry = structures.sphere(step, R=radius, mesh=mesh)
geometry = structures.nanodisc(step, R=radius, H=9, mesh=mesh)
geometry = structures.split_ring(step, R=8, W=4, H=3,  alphaG=30)
geometry = structures.lshape_round(step, L=15,W=8,H=2, DELTA=3, RAD=3)
material = materials.dummy(2.0+1j)  # dummy material with constant refindex n=2.0
material = materials.gold()  # dummy material with constant refindex n=2.0

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization(mesh))




## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
wavelengths = np.linspace(400, 1000, 5)           # spectrum
kwargs = dict(theta = [0.0,45,90], kSign=[-1,1])              # several polarizations

kwargs = [dict(theta = 0.0, kSign=-1),
          dict(theta = 33.0, kSign=-1),
          dict(theta = 71.0, kSign=1),
          dict(theta = 90.0, kSign=-1)]

efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)


visu.structure(sim)
print "N dipoles:", len(sim.struct.geometry)



img = structures.struct_to_image(sim, projection='XY')
plt.imshow(img, interpolation="none")
plt.show()





#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
EFIELD = core.scatter(sim, method='LU', verbose=True)



## spectrum via tools
spec_all = []
for field_kwargs in tools.get_possible_field_params_spectra(sim):
    wl, spec_ext = tools.calculate_spectrum(sim, field_kwargs, linear.extinct)
    spec_all.append([field_kwargs, spec_ext.T])



tools.save_simulation(sim, "test_save.sim")
sim2 = tools.load_simulation("test_save.sim")


tools.print_sim_info(sim)
print '\n'*3
tools.print_sim_info(sim2)


#%%
#==============================================================================
# Plot old/new
#==============================================================================
R = float(radius*step)

for aext_entry in spec_all:
    plt.plot(wl, aext_entry[1][0]/(np.pi*R**2), label=str(aext_entry[0]))
plt.legend(loc='best', fontsize=10)

plt.xlabel("wavelength (nm)")
    
plt.show()



#%% --- all field-configurations
field_confs = tools.get_field_indices(sim)
for i, fc in enumerate(field_confs):
    print i, fc

#%% --- get closest simulation entry

idx = tools.get_closest_field_index(sim, dict(theta=80, kSign=-1, wavelength=750))
print 'closest match --> index:', idx, ', dict:', sim.E[idx][0]




