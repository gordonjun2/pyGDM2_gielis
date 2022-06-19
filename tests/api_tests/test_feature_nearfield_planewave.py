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





from pyGDM import visu
from pyGDM import tools as toolsOld
from pyGDM import main as mainOLD



## ---------- Setup structure
mesh = 'cube'
step = 30.0
geometry = structures.rect_wire(step, L=30,H=3,W=4, mesh=mesh)
#geometry = structures.nanodisc(step, R=4,H=3, mesh=mesh)
#geometry = structures.hexagon(step, NSIDE=8,H=3, mesh=mesh)
material = materials.dummy(2.0); mat_string_old='20'    # dummy material with constant and real dielectric function
material = materials.silicon(); mat_string_old='si'       # silicon dielectric function
#material = materials.gold(); mat_string_old='au'          # gold dielectric function
#material = materials.alu(); mat_string_old='al'          # gold dielectric function

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))




## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
kwargs = dict(theta = [30.0])              # single polarizations
#kwargs = dict(theta = np.linspace(0,180, 90))              # several polarizations
wavelengths = [500]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)

print "N dipoles:", len(sim.struct.geometry)
visu.plot2Dstruct(geometry)





## --- TESTING - old pyGDM code
simOLD = mainOLD.genSimDict(step, 5000., geometry, n1,n2,n2, 
               wavelengths, kwargs['theta'], material=mat_string_old, mesh=mesh)


#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
print "\n\nOLD CODE ----"
EFIELDold = mainOLD.dipoles(simOLD, method='LU', verbose=True)

print "\n\nNEW CODE ----"
EFIELD = core.scatter(sim, method='LU', verbose=True)



    
    

#%%
#==============================================================================
# Compare Fields old/new
#==============================================================================
NFold = toolsOld.getFieldAsList(simOLD, EFIELDold, 1, 1, withIntensity=False)

xm,ym,zm = geometry.T
X,Y,Z = EFIELD[0][1].T
NFnew = np.array([xm,ym,zm, X.real,X.imag, Y.real,Y.imag, Z.real,Z.imag]).T


MAP = toolsOld.genNFmapXY(-500,500,50, -500,500,50, Z0=sim.struct.geometry.T[2].max()+sim.struct.step*1.5)


print "calc nearfield using 'old' pyGDM...", 
t0 = time.time()
NFmapOld = mainOLD.nearfield(simOLD, EFIELDold, MAP, 1,1, )
t1 = time.time()
NFmapOld = NFmapOld.T[:9].T
print "done. (time = {}ms)".format(1000.*(t1-t0))


#%% -- Group by incident field parameter-set
print "calc nearfield using pyGDM2...", 
t0 = time.time()
field_index = tools.get_closest_field_index(sim, dict(wavelength=400, theta=0))
NFmapEs, NFmapE1, NFmapBs, NFmapB1 = linear.nearfield(sim, field_index, MAP)
t1 = time.time()


X,Y,Z, EX,EY,EZ = NFmapE1.T
NFmapNew = np.array([X.real,Y.real,Z.real, 
                     EX.real,EX.imag, EY.real,EY.imag, EZ.real,EZ.imag]).T
print "done. (time = {}ms)".format(1000.*(t1-t0))

print "alpha pyGDM2: ", sim.struct.getPolarizability(wavelengths[0])[0]
print "alpha fortran:", mainOLD.getChi(simOLD, 1, SI=False, dim='3D')


## -- calculate relative deviation between new/old implementation
Inew = np.abs(NFmapNew.T[3]+1j*NFmapNew.T[4])**2 + \
       np.abs(NFmapNew.T[5]+1j*NFmapNew.T[6])**2 + \
       np.abs(NFmapNew.T[7]+1j*NFmapNew.T[8])**2
Iold = np.abs(NFmapOld.T[3]+1j*NFmapOld.T[4])**2 + \
       np.abs(NFmapOld.T[5]+1j*NFmapOld.T[6])**2 + \
       np.abs(NFmapOld.T[7]+1j*NFmapOld.T[8])**2

maxI = max([Inew.max(), Iold.max()])
NormMaxInew = np.array(NFmapNew.T[3:]/maxI**0.5)
NormMaxIold = np.array(NFmapOld.T[3:]/maxI**0.5)
RelDevIntensity = np.concatenate([NFmapNew.T[:3], (NormMaxInew-NormMaxIold)]).T

#%%

#==============================================================================
# Plot old/new
#==============================================================================
plt.figure(figsize=(10,4))
plt.subplot(121)
visu.plotNF2Dreal(NFmapOld, tit='OLD, Efield, Realpart', show=False)
plt.subplot(122)
visu.plotNF2Dreal(NFmapNew, tit='NEW, Efield, Realpart', show=False)
plt.show()

plt.figure(figsize=(13,3.5))
plt.subplot(131)
visu.plotNFcolor2Dintensity(NFmapOld, tit='OLD, Efield, Intensity', show=False)
plt.colorbar()
plt.subplot(132)
visu.plotNFcolor2Dintensity(NFmapNew, tit='NEW, Efield, Intensity', show=False)
plt.colorbar()
plt.subplot(133)
visu.plotNFcolor2Dintensity(RelDevIntensity, tit='Rel. difference, E-Intensity', show=False)
plt.colorbar()
plt.show()


