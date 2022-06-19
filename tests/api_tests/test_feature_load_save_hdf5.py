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




def saveh5(sim, fname):
    tables_filter_kwargs = dict(complib='blosc', complevel=5)
    import tables
    
    ## write E-field to h5
    FILTERS = tables.Filters(**tables_filter_kwargs)
    f = tables.open_file(fname + '.h5', 'w', filters=FILTERS)
    
    grp_E = f.create_group('/', 'E')
    
    E_dicts = []
    for i, ef in enumerate(sim.E):
        Edict = ef[0]
        E_dicts.append([Edict, None])   # remove E-field-data
        E = ef[1]
        f.create_carray(grp_E, 'E{}'.format(i), obj=E)
    
    ## write rest (without sim.E) to pickle file
    E = sim.E
    sim.E = E_dicts
    pickle.dump( sim, open( fname, "wb" ) )
    sim.E = E
    
    f.close()


def loadh5(fname):
    import tables, pickle
    
    ## load meta-data from pickle
    sim = pickle.load( open( fname, "rb" ) )
    
    ### efields from hdf5
    f = tables.open_file(fname + '.h5', 'r')
    for ef in f.root.E:
        enr = int(ef.name[1:])
        E = ef.read()
        sim.E[enr][1] = E
    f.close()
    
    return sim
    

#==============================================================================
# Testscript
#==============================================================================
## ---------- Setup structure
mesh = 'cube'
step = 20.0
radius = 3.
geometry = structures.sphere(step, R=radius, mesh=mesh)
material = materials.gold()  # dummy material with constant refindex n=2.0

n1, n2 = 1.5, 1.3  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization(mesh))




## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
wavelengths = np.linspace(400, 1000, 10)           # spectrum
kwargs = dict(theta = np.linspace(0.0,90,10), kSign=[-1,1])    # several configs

efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)

## ---------- Simulation initialization
sim = core.simulation(struct, efield)





#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
EF = core.scatter(sim, method='LU', verbose=True)


tools.save_simulation(sim, "test_save.sim", mode='pickle')
tools.save_simulation(sim, "test_save2.sim", mode='h5')
sim2 = tools.load_simulation("test_save.sim")
sim3 = tools.load_simulation("test_save2.sim")

#saveh5(sim, "test_save2.sim")
#sim3 = loadh5("test_save2.sim")


tools.print_sim_info(sim)
print '\n'*3
tools.print_sim_info(sim2)


#%%
#==============================================================================
# Plot old/new
#==============================================================================
R = float(radius*step)
wl, exti1 = tools.calculate_spectrum(sim, 0, linear.extinct)
qex1,qsc1,qab1 = exti1.T
wl, exti2 = tools.calculate_spectrum(sim2, 0, linear.extinct)
qex2,qsc2,qab2 = exti2.T
wl, exti3 = tools.calculate_spectrum(sim3, 0, linear.extinct)
qex3,qsc3,qab3 = exti3.T


plt.plot(wl, qex1)
plt.plot(wl, qex2, dashes=[2,2], label='pcl')
plt.plot(wl, qex3, dashes=[4,4], c='C3', label='h5')
plt.legend()
plt.show()
