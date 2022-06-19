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
TestDraw = 0


## Config
FIGSIZE = (1300,800)    # [width=0.8\textwidth]  -  (long) golden ratio: ~1.3 x 0.8 
FIGSIZE = (800, 600)    # use in subfigures -  (square) ratio: 1 x 1 



BG_COLOR = (1,1,1)


## camera position
#AZIM_SC=-45; ELEV_SC=55; DIST_SC=750
AZIM_SC=135; ELEV_SC=45; DIST_SC=750
FOCPT_SC=[25, 25, 30]






########################################################################
## Load data
########################################################################
#Xs,Ys,Zs = np.loadtxt("struct.txt").T





########################################################################
## MAIN PLOT
########################################################################

#==============================================================================
# Setup Figure
#==============================================================================
fig=mlab.figure(size=FIGSIZE, bgcolor=BG_COLOR)
visual.set_viewer(fig)


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
#drawRect(-350.,350., -350.,450., -20., 0.,  COLOR=ColorDarkGrey)
drawRect(-350.,350., -350.,450., -20., 0., OPACITY=0.5, COLOR=(0.8, 0.8, 0.9))

## Short Nano-Rod
drawNanorod(Xcenter=140,Ycenter=100,Zoffset=40, Lzyl=140, Rzylinder=45, COLOR=ColorGold)
drawSphere(Xcenter=-100, Ycenter=100, Zoffset=50, R=60, COLOR=ColorGold)
drawTriangle(Xcenter=0, Ycenter=-190, L=105, H=50, Zoffset=0, COLOR=ColorGold)





## Testing: skip expensive rendering, don't save to file
if TestDraw:
    mlab.view(AZIM_SC, ELEV_SC, DIST_SC, FOCPT_SC)
    mlab.show()
else:

    mlab.view(AZIM_SC, ELEV_SC, DIST_SC, FOCPT_SC)
    
    
    
    ## Save to file
    WD = os.path.basename(os.getcwd())
    filename = 'structure_3Dobjects.jpg'.format(WD)
    mlab.savefig(filename, figure=fig)#, magnification=2)
    mlab.show()














