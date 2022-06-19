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




# =============================================================================
# NEW CODE
# =============================================================================
## ---------- Setup structure
mesh = 'cube'
step = 7.0   # discretization in nm
L = 21       # cube side-length in nm
geometry = structures.rect_wire(step, L=L/step,H=L/step, W=L/step, mesh=mesh)
geometry.T[2] += step/2.

#material = materials.dummy(1.5); MATERIAL = 'TT'
material = materials.dummy(2.0); MATERIAL = '20'
#material = materials.silicon(); MATERIAL = 'si'

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))
struct = structures.center_struct(struct)



## ---------- Setup incident field
field_generator = fields.dipole_electric
field_generator = fields.dipole_magnetic
x0 = np.linspace(-250, 250, 25)
y0 = np.linspace(-250, 250, 25)
z0 = struct.geometry.T[2].max() + step/2. + 15

kwargs = dict(x0=x0, y0=y0, z0=z0,    # positions where to evaluate
              mx=0,my=0,mz=0)         # dipole orientation = (0,0,0) --> placeholder
wavelengths = [500]
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)

visu.structure(sim.struct.geometry, scale=0.5)
print "N dipoles:", len(sim.struct.geometry)



#%%
## --- main simulation
#core.scatter(sim, method='lu', verbose=True)
SBB_all = core.decay_rate(sim, method='lu', verbose=True)



#%%
# =============================================================================
# OLD CODE
# =============================================================================
from pyGDM import main, decay
from pyGDM import tools as tools_old

SPACING = 5000
N1 = n1; N2 = n2; N3 = N2
ELAMBDA = wavelengths
ATHETA = [0]

SIM = main.genSimDict(step, SPACING, geometry, N1,N2,N3, 
                        ELAMBDA, ATHETA, material=MATERIAL, mesh=mesh, DTYPE='f')
#print "polarizability per cell OLD:", main.getChi(SIM,Ilambda=1)


MAP = tools_old.genNFmapXY(x0.min(),x0.max(),len(x0), y0.min(),y0.max(),len(y0), Z0=z0)

## scan electric dipole  
if sim.efield.field_generator == fields.dipole_electric:
    magnetic = 0
    dp_type = 'electric'
elif sim.efield.field_generator == fields.dipole_magnetic:
    magnetic = 1
    dp_type = 'magnetic'


## --- main simulation
SBB = decay.decayRatePropa(SIM, 1, MAP, magnetic=magnetic, method='lu', verbose=1)




#%%
# =============================================================================
# Plot decay for several dipole orientations
# =============================================================================

## list of dipole test orientations
dp_list = [[1, 0, 0], 
           [0, 1, 0], 
           [0, 0, 1]]


plt.figure(figsize=(8,8))
for i, (mx,my,mz) in enumerate(dp_list):

    SBB_map_old = decay.calcDecayMap(SIM, 1, MAP, SBB, mx,my,mz)
    NF_MAP_old, ext = tools_old.mapToGridXY(MAP, SBB_map_old)
    
    SBB_map_new = linear.decay_eval(sim, SBB_all[0], mx,my,mz, verbose=1)
    NF_MAP_new, ext = tools_old.mapToGridXY(MAP, SBB_map_new.T[-1])
    
    
    plt.subplot(3,2,2*i+1, aspect='equal')
    plt.title("old code\n{} dipole || ({},{},{})".format(dp_type,mx,my,mz))
    plt.imshow(NF_MAP_old, extent=ext, cmap='jet')
    plt.colorbar(label='gamnma / gamma_0')
    
    plt.subplot(3,2,2*i+2, aspect='equal')
    plt.title("new code\n{} dipole || ({},{},{})".format(dp_type,mx,my,mz))
    plt.imshow(NF_MAP_new, extent=ext, cmap='jet')
    #visu.scalarfield(SBB_map_new, cmap='jet', show=False)
    plt.colorbar(label='gamnma / gamma_0')
    

plt.tight_layout()
plt.show()











