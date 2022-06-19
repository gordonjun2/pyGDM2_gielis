#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from mayavi import mlab

from mayavi_objects import *
from pyGDM import structures
from pyGDM2 import visu3d


mesh = 'hex'
STEP=10
NSIDE=20

S1 = structures.nanorod(STEP,L=14,R=4,caps='round', mesh=mesh)
S1 = structures.rotateXY(S1, 0)
S1.T[0] += 14*STEP
S1.T[1] += 10*STEP

S2 = structures.prism(STEP,NSIDE=NSIDE,H=5, mesh=mesh)
S2 = structures.rotateXY(S2, 0)
S2.T[0] -= 0*STEP
S2.T[1] += 0.42*NSIDE*STEP
S2.T[1] -= 18*STEP

S3 = structures.sphere(STEP,R=5.2, mesh=mesh)
S3.T[0] -= 10*STEP
S3.T[1] += 10*STEP

S = np.concatenate([S1,S2,S3])
S.T[2] += STEP

np.savetxt('struct_hex.txt', S, fmt='%.3g')
print( "Ndipoles =", len(S))
visu3d.structure(S, scale=.8, color=ColorGold, opacity=1, show=1)



## --- Higher resolution Version
mesh='cube'
STEP2 = 10
FACTOR = STEP/STEP2
NSIDE = int(NSIDE * FACTOR)

S1 = structures.nanorod(STEP2,L=14*FACTOR,R=4*FACTOR,caps='round', mesh=mesh)
S1 = structures.rotateXY(S1, 0)
S1.T[0] += 14*FACTOR*STEP2
S1.T[1] += 10*FACTOR*STEP2

S2 = structures.prism(STEP2,NSIDE=NSIDE,H=5*FACTOR, mesh=mesh)
S2 = structures.rotateXY(S2, 0)
S2.T[0] -= 0*STEP2
S2.T[1] += 0.42*NSIDE*STEP2
S2.T[1] -= 18*FACTOR*STEP2

S3 = structures.sphere(STEP2,R=5.2*FACTOR, mesh=mesh)
S3.T[0] -= 10*FACTOR*STEP2
S3.T[1] += 10*FACTOR*STEP2

S = np.concatenate([S1,S2,S3])
#S.T[2] += STEP2

np.savetxt('struct_cube.txt', S, fmt='%.3g')
print( "Ndipoles =", len(S))


## test plot
#drawRect(-700.,700., -600.,600., -20., 0., ColorGrey)


#mlab.show()

