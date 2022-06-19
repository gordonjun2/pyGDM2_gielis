# encoding: utf-8
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






#==============================================================================
# Testscript
#==============================================================================
## ---------- Setup structure
mesh = 'hex'
step = 15.0
radius = 5.5
geometry = structures.sphere(step, R=radius, mesh=mesh)
material = materials.dummy(2.0)
material = materials.gold()

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization(mesh))




## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
wavelengths = np.logspace(2, 3, 50)[:40]           # spectrum
#wavelengths = np.linspace(800, 1000, 5)           # spectrum
wavelengths = np.linspace(400, 1000, 30)           # spectrum
kwargs = dict(theta = [0.0])              # several polarizations
#kwargs = dict(theta = np.linspace(0,180,90))              # several polarizations
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)


#visu.plot2Dstruct(sim.struct.geometry, scale=10)
print "N dipoles:", len(sim.struct.geometry)



#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
print "\n\nNEW CODE ----"
EFIELD = core.scatter(sim, method='lu', verbose=True)
#EFIELD = core.scatter(sim, method='LU', verbose=True)
#tools.save_simulation(sim, "testsim.sim")

#t0 = time.time()
#sim = tools.load_simulation("testsim.sim")
#EFIELD = sim.E
#print "time reload: {:.3f}s".format(time.time() - t0)

#%%


## spectrum via tools
wl, spec_ext = tools.calculate_spectrum(sim, 0, linear.extinct)
aext, asca, aabs = spec_ext.T


plt.subplot()
plt.plot(wl, asca, label='sca')
plt.plot(wl, aext, label='ext')
plt.plot(wl, aabs, label='abs')
plt.show()


