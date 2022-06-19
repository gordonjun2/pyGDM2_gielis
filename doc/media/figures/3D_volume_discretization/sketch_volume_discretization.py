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
sys.path.append("../../")
from mayavi_objects import *




########################################################################
## VISUAL CONFIGURATION
########################################################################
## disable performance intensive stuff
TestDraw = False


## Config
#FIGSIZE = (1300,800)    # [width=0.8\textwidth]  -  (long) golden ratio: ~1.3 x 0.8 
FIGSIZE = (800, 600)    # use in subfigures -  (square) ratio: 1 x 1 


BG_COLOR = (1,1,1)


BG_COLOR = (1,1,1)


## camera position
AZIM_SC=135; ELEV_SC=45; DIST_SC=750
FOCPT_SC=[25, 25, 30]






#%%
########################################################################
## Load / Generate data
########################################################################
STRUCT = np.loadtxt("struct_cube.txt")
SAVE_SUFFIX='cube'
scaling = 0.8

STRUCT = np.loadtxt("struct_hex.txt")
SAVE_SUFFIX='hex'
scaling = 0.7





#%%
########################################################################
## MAIN PLOT
########################################################################

#==============================================================================
# Setup Figure
#==============================================================================
fig=mlab.figure(size=FIGSIZE, bgcolor=BG_COLOR)
visual.set_viewer(fig)

#
### Adapt Lights
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
drawRect(-350.,350., -350.,450., -20., 0., OPACITY=0.5, COLOR=(0.8, 0.8, 0.9))


### Objective
#drawObjective(Xcenter=0, Ycenter=0, Zoffset=700)


### Coordinate System
#Arrow_From_A_to_B(600,600,700, 1100,600,700, Lcone=0.2,Rcone=0.05,Rshaft=0.015, color=(0.2,0.2,0.2))
#Arrow_From_A_to_B(600,600,700, 600,1100,700, Lcone=0.2,Rcone=0.05,Rshaft=0.015, color=(0.2,0.2,0.2))
#Arrow_From_A_to_B(600,600,700, 600,600,1200, Lcone=0.2,Rcone=0.05,Rshaft=0.015, color=(0.2,0.2,0.2))



## Sample
Xs,Ys,Zs = STRUCT.T
mlab.points3d(Xs,Ys,Zs, scale_factor=10*scaling, color=ColorGold, mode='cube')






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
filename = 'structure_discretized_{}.jpg'.format(SAVE_SUFFIX)
mlab.savefig(filename, figure=fig)#, magnification=2)
mlab.show()














