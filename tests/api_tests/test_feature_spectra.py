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





from pyGDM import visu
#from pyGDM import tools
from pyGDM import main as mainOLD


#==============================================================================
# Config
#==============================================================================
DoOld = 1    # whether or not compare to old code
#DoOld = 0     # whether or not compare to old code





#==============================================================================
# Testscript
#==============================================================================
## ---------- Setup structure
mesh = 'hex'
step = 15.0
radius = 3.5
geometry = structures.sphere(step, R=radius, mesh=mesh)
material = materials.dummy(2.0); matold='20'  # dummy material with constant refindex n=2.0
#material = materials.gold(); matold='au'  # dummy material with constant refindex n=2.0
#material = materials.gold(3); matold='au'  # dummy material with constant refindex n=2.0
#material = materials.fromFile("refindex_gold.txt")

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








## --- TESTING - compare with "old" pyGDM
if DoOld:
    simOLD = mainOLD.genSimDict(step, 5000., geometry, n1,n2,n2, 
                   wavelengths, kwargs['theta'], material=matold, mesh=mesh)

#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
if DoOld:
    print "\n\nOLD CODE ----"
    EFIELDold = mainOLD.dipoles(simOLD, method='LU', verbose=True)

print "\n\nNEW CODE ----"
#EFIELD = core.scatter(sim, method='CG', verbose=True)
EFIELD = core.scatter(sim, method='lu', verbose=True)
#tools.save_simulation(sim, "testsim.sim")

#t0 = time.time()
#sim = tools.load_simulation("testsim.sim")
#EFIELD = sim.E
#print "time reload: {:.3f}s".format(time.time() - t0)

#%%
## --- spectra
if DoOld:
    print "\n\n\n------new extinct"
    wl, aextOLD, ascaOLD, aabsOLD = mainOLD.extinct(simOLD, EFIELDold)

    print "\n\n\n------new extinct"

## calculate spectrum manually
aext2 = []; asca2 = []; aabs2 = []
t0=time.time()
for i, E in enumerate(sim.E):
    params = E[0]
    if params['theta'] == 0:
        s = linear.extinct(sim, i)
        aext2.append(s[0])
        asca2.append(s[1])
        aabs2.append(s[2])


## spectrum via tools
field_kwargs = tools.get_possible_field_params_spectra(sim)[0]
wl, spec_ext = tools.calculate_spectrum(sim, field_kwargs, linear.extinct)
wl, spec_heat = tools.calculate_spectrum(sim, field_kwargs, linear.heat)


plt.figure(figsize=(8,6))

plt.subplot(221); plt.title("extinction")
plt.plot(sim.efield.wavelengths, aext2, 'g-')
plt.plot(wl, spec_ext.T[0], 'k--')

plt.subplot(222); plt.title("scattering")
plt.plot(sim.efield.wavelengths, asca2, 'g-')
plt.plot(wl, spec_ext.T[1], 'k--')

plt.subplot(223); plt.title("absorption")
plt.plot(sim.efield.wavelengths, aabs2, 'g-')
plt.plot(wl, spec_ext.T[2], 'k--')

plt.subplot(224); plt.title("heat")
plt.plot(wl, spec_heat/1E3, 'r')

plt.tight_layout()
plt.show()
#%%
#==============================================================================
# Plot old/new
#==============================================================================
#R = float((radius)*step)
#
#plt.figure(figsize=(13,3))
#if DoOld:
#    plt.subplot(131); plt.title("old implementation")
#    for aextOLD_entry in aextOLD:
#        plt.plot(wl/R, aextOLD_entry/(np.pi*R**2), label='ext')
#    for aextOLD_entry in aabsOLD:
#        plt.plot(wl/R, aextOLD_entry/(np.pi*R**2), label='abs')
#    for aextOLD_entry in ascaOLD:
#        plt.plot(wl/R, aextOLD_entry/(np.pi*R**2), label='sca')
#    plt.legend(loc='best')
#    #plt.ylabel("x-section (nm^2)")
#    plt.ylabel("ext. efficiency")
#    
#    plt.subplot(133); plt.title("rel. diff")
#    for aext_entry,aext_entryO in zip(aext, aextOLD):
#        plt.plot(aext_entry[1]/R, (aext_entry[2]-aext_entryO)*2./(aext_entry[2].max()+aext_entryO.max()))
#    plt.legend(loc='best')
#
#if DoOld:
#    plt.subplot(132); plt.title("new implementation")
#else:
#    plt.subplot()
#
#for aext_entry in aext:
#    plt.plot(aext_entry[1], aext_entry[2], label='ext'+str(aext_entry[0]))
#for aext_entry in aabs:
#    plt.plot(aext_entry[1], aext_entry[2], label='abs'+str(aext_entry[0]))
#for aext_entry in asca:
#    plt.plot(aext_entry[1], aext_entry[2], label='sca'+str(aext_entry[0]))
#plt.legend(loc='best', fontsize=10)
##plt.xlabel("wavelength (nm)")
#plt.xlabel("wavelength/radius")
#    
#plt.show()
#
#
#np.savetxt("sphere_scata_test.txt", np.array([aext_entry[1]/R, aext_entry[2]/(np.pi * R**2)]), 
#           header='R={}'.format(R), fmt="%.5g")
#
#print EFIELD[0][0]
#print np.shape(EFIELD[0][1])
#print len(EFIELD[0][1])
#
#
#
#tools.save_simulation(sim, "test_save.sim")
#sim2 = tools.load_simulation("test_save.sim")
#
#tools.print_sim_info(sim2)

