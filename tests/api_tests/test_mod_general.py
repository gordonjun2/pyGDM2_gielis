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





from pyGDM import visu
from pyGDM import tools
from pyGDM import main as mainOLD


## --- whether or not compare with old pyGDM implementation
doOld = 1





## ---------- Setup structure
mesh = 'hex'
step = 10.0
geometry = structures.rect_wire(step, L=40,H=2,W=5, mesh=mesh)
geometry = structures.nanodisc(step,   R=5,H=3,ELONGATED=10, mesh=mesh)
material = materials.dummy(2.0); mat_string_old='20'    # dummy material with constant and real dielectric function
#material = materials.silicon(); mat_string_old='si'       # silicon dielectric function

# gold dielectric function
#material = materials.gold(); mat_string_old='au' # linear interpolation (numpy)
#material = materials.gold(interpolate_order=3); mat_string_old='au' # cubic interpolation

n1, n2 = 1.0, 1.0  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))




## ---------- Setup incident field
field_generator = fields.planewave        # planwave excitation
#kwargs = dict(theta = [0.0, 45.0, 90.0], kSign=[-1,1])              # several configurations
#wavelengths = np.linspace(400,1000,25)                     # spectrum
kwargs = dict(theta = [0.0])              # several polarizations
wavelengths = [400]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)

plt.subplot(aspect='equal')
visu.plot2Dstruct(sim.struct.geometry, scale=0.5, show=False)
#plt.savefig("structure.png")
plt.show()
visu.plot2Dstruct(sim.struct.geometry, projection='YZ',scale=0.5)

print "N dipoles:", len(sim.struct.geometry)
print "Approx. size of sim. object before simulations:", sys.getsizeof(pickle.dumps(sim))/1024, "kB"


#np.savetxt("stadium_R150nm_L400nm_S10nm_H3meshpoints_hex.txt", geometry, fmt="%.5g")
#np.savetxt("stadium_R150nm_L400nm_S10nm_H3meshpoints_cube.txt", geometry, fmt="%.5g")
#np.savetxt("stadium_R150nm_L400nm_S10nm_H1meshpoints_cube.txt", geometry, fmt="%.5g")






## --- TESTING - old pyGDM code
if doOld:
    simOLD = mainOLD.genSimDict(step, 5000., geometry, n1,n2,n2, 
                   wavelengths, kwargs['theta'], material=mat_string_old, mesh=mesh)
    print "alpha pyGDM2: ", sim.struct.getPolarizability(wavelengths[0])[0]
    print "alpha fortran:", mainOLD.getChi(simOLD, 1, SI=False, dim='3D')

#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
if doOld:
    print "\n\nOLD CODE ----"
    EFIELDold = mainOLD.dipoles(simOLD, method='LU', verbose=True, nthreads=1)

print "\n\nNEW CODE ----"
EFIELD = core.scatter(sim, method='lu', verbose=True)
#EFIELD = core.scatter(sim, method='numpyinv', verbose=True, nthreads=-1)
#EFIELD = core.scatter(sim, method='scipyinv', verbose=True, nthreads=-1)
#EFIELD = core.scatter(sim, method='pinv2', verbose=True, nthreads=-1)
#EFIELD = core.scatter(sim, method='superlu', verbose=True, nthreads=-1)
#EFIELD = core.scatter(sim, method='dyson', verbose=True, nthreads=-1)
#EFIELD = core.scatter(sim, method='cg', verbose=True, nthreads=-1)


xm,ym,zm = geometry.T
X,Y,Z = EFIELD[0][1].T
NFnew = np.array([xm,ym,zm, X.real,X.imag, Y.real,Y.imag, Z.real,Z.imag]).T


#==============================================================================
# Compare CAM0 Matrices (testing, modified output of dipoles/scatter to "CAM0")
#==============================================================================
if doOld:
    import scipy
    if type(EFIELD) == scipy.sparse.csc.csc_matrix:
        count = counttot = 0
        reltolerance = 1E-6
        x = (EFIELD - EFIELDold)
        rows,cols = x.nonzero()
        for row,col in zip(rows,cols):
            counttot += 1
            if np.abs(x[row,col]) > reltolerance: 
                count += 1
                print "dev >", reltolerance, ":", ((row,col), x[row,col],np.abs(x[row,col]))
        print "Elements with realtive derivation > ", reltolerance, ": ", count, "/" , counttot

    
    
    #%%
    #==============================================================================
    # Compare Fields old/new
    #==============================================================================
    NFold = tools.getFieldAsList(simOLD, EFIELDold, 1, 1, withIntensity=False)
    
    
    
    
    
    #==============================================================================
    # Compare relative errors of fields new/old
    #==============================================================================
    count = counttot = 0
    reltolerance = 1E-5
    x = (NFnew - NFold).T
    totabs = np.abs(x[3] + 1j*x[4]) + np.abs(x[5] + 1j*x[6]) + np.abs(x[7] + 1j*x[8])
    for i, absval in enumerate(totabs):
        counttot += 1
        if absval > reltolerance: 
            count += 1
    #        if absval > reltolerance*10.0: 
    #            np.set_printoptions(precision=4)
    #            print ">x10 rel. diff: idx", i,":", absval
    #            print " -- New:", NFnew[i][3:]
    #            print " -- Old:", NFold[i][3:]
    #            print ""
    print "Elements with relative deviation > ", reltolerance, ": ", count, "/" , counttot
    
    

#%% -- Group by incident field parameter-set
print '\n ----- calculated spectra:'
sortedFields = core._sortByParameterConfig(sim, EFIELD)
for i in sortedFields:
    print 
    for j in i:
        print j[0]
#%%

#==============================================================================
# Plot old/new
#==============================================================================
plt.figure(figsize=(10,4))
if doOld:
    plt.subplot(121, aspect='equal')
    visu.plotNF2Dreal(NFold, tit='OLD, Efield, Re', show=False)
    plt.subplot(122, aspect='equal')
visu.plotNF2Dreal(NFnew, tit='NEW, Efield, Re', show=False)
plt.show()


