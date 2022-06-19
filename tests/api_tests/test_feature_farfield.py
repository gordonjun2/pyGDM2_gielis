# encoding: utf-8
import os, sys
import pickle
 
import time
import copy

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




#from pyGDM import visu
from pyGDM import main as mainOLD
from pyGDM import tools as toolsOld
from pyGDM import fields as fieldsOld




## ---------- Setup structure
mesh = 'cube'
step = 30.0
geometry = structures.rect_wire(step, L=10,H=4,W=4, mesh=mesh)
#geometry = structures.nanodisc(step, R=4,H=3, mesh=mesh)
#geometry = structures.hexagon(step, NSIDE=8,H=3, mesh=mesh)
geometry.T[2] += step
geometry = structures.center_struct(geometry)


#geometry = np.loadtxt("structure_opt_antenna_planewave.txt"); step = 40.0

material = materials.dummy(2.0); mat_string_old='20'    # dummy material with constant and real dielectric function
#material = materials.silicon(); mat_string_old='si'     # silicon dielectric function
#material = materials.gold(); mat_string_old='au'        # gold dielectric function
#material = materials.alu(); mat_string_old='al'         # alu dielectric function




n1, n2 = 1.0, 1.0  # constant environment

#geometry = structures.centerStruct(geometry)
struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))
#struct = structures.centerStruct(struct)



## ---------- Setup incident field
field_generator = fields.dipole_electric        # planwave excitation
kwargs = dict(x0=0.,y0=50.,z0=70., mx=1.,my=0.,mz=0.)
kwargs_old = copy.deepcopy(kwargs)
#field_generator = fields.planewave        # planwave excitation
#kwargs = dict(theta = [0.0, 90.0])              # single polarizations

wavelengths = [800.]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)
sim = structures.center_struct(sim)

print "N dipoles:", len(sim.struct.geometry)
visu.structure(sim)



#kwargs_old = {}
Ein = fieldsOld.getE0elecDipole
## --- TESTING - old pyGDM code
simOLD = mainOLD.genSimDict(step, 5000., geometry, n1,n2,n2, 
               wavelengths, [0.0, 90.0], material=mat_string_old, mesh=mesh)


#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
print "\n\nOLD CODE ----"
EFIELDold = mainOLD.dipoles(simOLD, method='LU', Ein=Ein, einKwargs=kwargs_old, verbose=True)

print "\n\nNEW CODE ----"
EFIELD = core.scatter(sim, method='LU', verbose=True)
print "DONE"

#%%
NF = tools.get_field_as_list_by_fieldindex(sim, 0)

    
    

#%%
#==============================================================================
# Compare Fields old/new
#==============================================================================
R=10000.; Nteta=50; Nphi=72
tetamax = np.pi/2.
## --- OLD
print "scatter old"
n_subst = n1; alambda = wavelengths[0]
ilambda,iteta, = 1, 1
theta2d,phi2d, I2d = mainOLD.dipemFF_plus_zero_dipole(
                                simOLD, EFIELDold, ilambda,iteta, EinKwargs=kwargs_old,
                                R=R,NTETASCA=Nteta,NPHI=Nphi+1, NPOL=4,
                                add_zero=0,
                                tetamax=tetamax, nthreads=-1)

## --- NEW
print "scatter new"
tetalist, philist, I_sc, I_tot, I0 = linear.farfield(
                                sim, field_index=0,
                                r=R, tetamin=0, tetamax=tetamax, Nteta=Nteta, Nphi=Nphi, 
                                polarizerangle='none')
print "done"


Intens = I_sc
#Intens = I_tot
#Intens = I0




### --- for plotting only: Add 360degrees (copy of 0 degrees)
#PHI = philist; THETA = tetalist; I = Intens
Nteta, Nphi = I_sc.shape
THETA = np.concatenate([tetalist.T, [tetalist.T[-1]]]).T
PHI = np.concatenate([philist.T, [np.ones(philist.T[-1].shape) * 2*np.pi]]).T - np.pi/float(Nphi)
I = np.concatenate([Intens.T, [Intens.T[-1]]]).T

THETAo = theta2d.reshape((Nteta, Nphi+1))
PHIo = phi2d.reshape((Nteta, Nphi+1)) - np.pi/float(Nphi+1)
Io = I2d.reshape((Nteta, Nphi+1))



## --- test-plot
plt.figure(figsize=(12,4))
plt.subplot(131, polar=True); plt.title('OLD')
plt.pcolormesh(PHIo, THETAo*180./np.pi, Io, edgecolors='face', rasterized=True, snap=True)
plt.colorbar()

plt.subplot(132, polar=True); plt.title('NEW')
visu.farfield_pattern_2D(tetalist, philist, Intens, show=False)
plt.colorbar()

plt.subplot(133, polar=True); plt.title('new-old\nrel. diff. (1=100%)')
maxI = max([I.max(), Io.max()])
plt.pcolormesh(PHI, THETA*180./np.pi, (I/maxI-Io/maxI) , edgecolors='face', rasterized=True, snap=True)
plt.colorbar()



plt.show()






## =============================================================================
## NEW DATA ONLY
## =============================================================================
#I = Intens
#
#
#tools.print_sim_info(sim, verbose=1)
### --- plot
#plt.figure(figsize=(10,5))
##    plt.subplot(121, aspect="equal")
##    visu.structure(sim, scale=0.15, show=False)
#
#plt.figure(figsize=(12,3.5))
#
#
#plt.subplot(231, aspect="equal")
#plt.title("phase X")
#NF_imagX = tools.get_field_as_list_by_fieldindex(sim, 0).T[ [0,1,2,3,4,5] ]
#colorcode1 = np.angle(NF_imagX[3]) * 180./np.pi
#
#plt.scatter(NF_imagX[0], NF_imagX[1], c=colorcode1, cmap='bwr', marker='s', s=10)
#plt.colorbar(label=r'phase of $E_x$ (deg)')
#
#plt.subplot(234, aspect="equal")
#plt.title("phase Y")
#NF_imagX = tools.get_field_as_list_by_fieldindex(sim, 0).T[ [0,1,2,3,4,5] ]
#colorcode1 = np.angle(NF_imagX[4]) * 180./np.pi
#
#plt.scatter(NF_imagX[0], NF_imagX[1], c=colorcode1, cmap='bwr', marker='s', s=10)
#plt.colorbar(label=r'phase of $E_y$ (deg)')
##    plt.clim([0, 180])
#
#plt.subplot(132, aspect="equal")
#plt.title("nearfield energy")
#NF_imagX = tools.get_field_as_list_by_fieldindex(sim, 0).T[ [0,1,2,3,4,5] ]
#colorcode2 = np.abs(NF_imagX[3]**2 + NF_imagX[4]**2 + NF_imagX[5]**2)
#
#plt.scatter(NF_imagX[0], NF_imagX[1], c=colorcode2, cmap='jet', marker='s', s=10)
#plt.colorbar(label=r'|E|^2')
##    plt.clim([0, 180])
#
#
#
#plt.subplot(133, polar=True)
#plt.title("farfield intensity")
#visu.farfield_pattern_2D(tetalist, philist, I, show=False)
#plt.colorbar()
#
#
#
#plt.tight_layout()
#plt.show()
#
#





