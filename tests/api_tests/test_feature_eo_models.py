# encoding: utf-8
"""
Created on May 6, 2017

@author: P. R. Wiecha

Example script demonstrating `pyGDM2.EO`:
    evolutionary optimizaton of plasmonic nanostructure geometry for 
    directional scattering

"""

import numpy as np
import matplotlib.pyplot as plt


from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core

from pyGDM2.EO import models as eo_models
from pyGDM2.EO import problems as eo_problems
from pyGDM2.EO.core import do_eo


from PyGMO import algorithm


#==============================================================================
# Setup pyGDM part
#==============================================================================
## ---------- Setup structure
mesh = 'cube'
step = 40
material = materials.dummy(2.0)    # material: constant and real dielectric function
material = materials.gold()        # material: gold
n1, n2 = 1.0, 1.0         # constant environment

## --- Empty dummy-geometry, will be replaced on run-time by EO trial geometries
geometry = []       

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))



## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
kwargs = dict(theta = [0.0])              # several polarizations
wavelengths = [800]                       # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)




#==============================================================================
# setup evolutionary optimization
#==============================================================================
## --- output folder and file-names
results_folder = 'eo_out'
results_suffix = 'test_eo'


## --- structure model and optimizaiton problem

#model = eo_models.BlockModel(sim, 
#                             30, 1,1,1, [-10,10],
#                             forbidden=[], symmetric=False)


model = eo_models.MultiRectAntenna(sim, 
                             N_antennas=15, limits_W=[2,10], limits_L=[2,10], 
                             limits_pos_x=[-100,100], limits_pos_y=[-100,100], 
                             height=3)




model.print_info()
model.plot_structure(scale=0.4)


















