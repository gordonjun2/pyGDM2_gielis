#!/usr/bin/python
# -*- coding: utf-8 -*-

## General
import numpy as np
import sys
import os

## Mayavi Imports
from mayavi import mlab
import mayavi.tools.pipeline
from tvtk.tools import visual


## Custom 3D Objects
from mayavi_objects import *




########################################################################
## VISUAL CONFIGURATION
########################################################################
## disable performance intensive stuff
TestDraw = False


## Config
FIGSIZE = (1000,500)    # [width=0.8\textwidth]  -  (long) golden ratio: ~1.3 x 0.8 
#FIGSIZE = (800, 800)    # use in subfigures -  (square) ratio: 1 x 1 


BG_COLOR = (1,1,1)


## camera position
AZIM_SC=45; ELEV_SC=87; DIST_SC=4200
FOCPT_SC=[0, 0, 400]

COLOR1 = (255/255., 170/255., 0/255.)
COLOR2 = (217/255., 0/255., 255/255.)
COLOR3 = (0/255., 205/255., 19/255.)
COLOR4 = (0/255., 137/255., 205/255.)
COLOR = ColorGold



#%%
########################################################################
## Load / Generate data
########################################################################
from pyGDM2 import structures

STEP = 25

S1 = structures.rect_wire(STEP, 14,3,9)


S1 = structures.sphere(STEP, 10)
S1 = structures.center_struct(S1)
S1.T[0] += 950
S1.T[1] += 950

S2 = structures.sphere(STEP,6)
S2 = structures.center_struct(S2)
S2.T[0] += 1050

S3 = structures.sphere(STEP, 4)
S3 = structures.center_struct(S3)
S3.T[0] += 50
S3.T[1] -= 600

S4 = structures.sphere(STEP, 8)
S4 = structures.center_struct(S4)
S4.T[0] -= 550
S4.T[1] += 550

S1 = np.concatenate([S1,S2,S3, S4])
S1.T[2] -= STEP



#%%
########################################################################
## MAIN PLOT
########################################################################

#==============================================================================
# Setup Figure
#==============================================================================
fig=mlab.figure(size=FIGSIZE, bgcolor=BG_COLOR)
visual.set_viewer(fig)


## Adapt Lights
#fig.scene.light_manager.lights[0].move_to(45.0, 10.0)
#fig.scene.light_manager.lights[0].intensity=1.0
#
#fig.scene.light_manager.lights[1].move_to(-30.0, -60.0)
#fig.scene.light_manager.lights[1].intensity=0.7
#
#fig.scene.light_manager.lights[2].move_to(-30.0, 60.0)
#fig.scene.light_manager.lights[2].intensity=0.3






#==============================================================================
# Main Scene
#==============================================================================
## Sample Table
drawRect(-1500.,1500., -1500.,1500., -350.,-15., COLOR=(0.7,0.7,0.7), OPACITY=0.3, rotAng=(45), rotAx=(0,0,1))
drawRect(-1500.,1500., -1500.,1500., -15.,900., COLOR=(0.7,0.7,0.9), OPACITY=0.15, rotAng=(45), rotAx=(0,0,1))
drawRect(-1500.,1500., -1500.,1500., 900.,1400., COLOR=(0.7,0.7,0.7), OPACITY=0.3, rotAng=(45), rotAx=(0,0,1))


### Objective
#drawObjective(Xcenter=0, Ycenter=0, Zoffset=700)


### Coordinate System
#Arrow_From_A_to_B(600,600,700, 1100,600,700, Lcone=0.2,Rcone=0.05,Rshaft=0.015, color=(0.2,0.2,0.2))
#Arrow_From_A_to_B(600,600,700, 600,1100,700, Lcone=0.2,Rcone=0.05,Rshaft=0.015, color=(0.2,0.2,0.2))
#Arrow_From_A_to_B(600,600,700, 600,600,1200, Lcone=0.2,Rcone=0.05,Rshaft=0.015, color=(0.2,0.2,0.2))



## Sample
#for iy in range(4):
#    for ix in range(10):
Xs,Ys,Zs = S1.T
ix = iy = iz = 0
mlab.points3d(Xs+ix*STEP*4,Ys-iy*STEP*8,Zs+STEP+5, scale_factor=STEP*0.8, color=COLOR, mode='cube')
#mlab.points3d([0], [0], [50], scale_factor=5, color=(0,0,0), mode='cube')






## Testing: skip expensive rendering, don't save to file
if TestDraw:
    mlab.view(AZIM_SC, ELEV_SC, DIST_SC, FOCPT_SC)
    mlab.show()

#==============================================================================
# Performance intensive stuff
#==============================================================================
# ...















mlab.view(AZIM_SC, ELEV_SC, DIST_SC, FOCPT_SC)



## Save to file
WD = os.path.basename(os.getcwd())
filename = 'sketch_reference_system_raw.png'.format(WD)
mlab.savefig(filename, figure=fig)#, magnification=2)
mlab.show()














