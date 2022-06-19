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
#sys.path.append("../../")
#from mayavi_objects import *

from pyGDM2 import structures, materials, core, tools, linear, fields
from pyGDM2 import visu3d


########################################################################
## VISUAL CONFIGURATION
########################################################################
## disable performance intensive stuff
TestDraw = 0


## Config
FIGSIZE = (800,800)    # [width=0.8\textwidth]  -  (long) golden ratio: ~1.3 x 0.8 

BG_COLOR = (1,1,1)

## camera position
AZIM_SC=45; ELEV_SC=78; DIST_SC=400
FOCPT_SC=[0, 0, 75]



## suffix for filename

#plottype = 'geo'
#savename_prefix = 'nanosphere_discretized'
#
#plottype = 'geo_zoom'
#savename_prefix = 'nanosphere_discretized_zoom'

plottype = 'vectors'
savename_prefix = 'nanosphere_discretized_vectors'

plottype = 'vectors_large'
savename_prefix = 'nanosphere_discretized_vectors_large'


#plottype = 'single_red_arrow'
#savename_prefix = 'single_red_arrow'

#plottype = 'single_blue_arrow'
#savename_prefix = 'single_blue_arrow'


########################################################################
## Generate data
########################################################################
1
## ---------- structure: nano-sphere
step = 25.
R = 3.0

mesh = 'cube'
geo = structures.sphere(step, R=R, mesh=mesh)
geo2 = structures.sphere(step/2, R=2*R, mesh=mesh)
geo.T[2] += step/2
geo = structures.center_struct(geo)
geo2 = structures.center_struct(geo2)

## --------- simulation for quiverplot
material = materials.silicon()
material2 = materials.silicon()
#material2 = materials.dummy(5)
n1, n2 = 1.0, 1.0

struct = structures.struct(step, geo, material, n1,n2, 
                                   structures.get_normalization(mesh))
struct2 = structures.struct(step/2, geo2, material2, n1,n2, 
                                   structures.get_normalization(mesh))

#wavelengths = np.arange(650,770,20)
wavelengths = [5000]
field_generator = fields.planewave
kwargs = dict(theta = [0])
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)
wavelengths2 = [600]
kwargs2 = dict(theta = [-45])
efield2 = fields.efield(field_generator, wavelengths=wavelengths2, kwargs=kwargs2)

sim = core.simulation(struct, efield)
sim_large = core.simulation(struct2, efield2)
core.scatter(sim)
core.scatter(sim_large)



#==============================================================================
# Sphere
#==============================================================================
def drawSphere(Xcenter=0, Ycenter=0, Zoffset=0, R=50,
                OPACITY=1, COLOR=(.8,.8,.8), Nphi=35, Nteta=50):

    phi    = np.linspace(0, 2.*np.pi, Nphi) + 1
    theta  = np.linspace(0, np.pi, Nteta)

    ## sphere
    zSphcap1 = [[Zoffset + R*np.sin(t)*np.cos(p) for p in phi] for t in theta]
    ySphcap1 = [[Ycenter + R*np.sin(t)*np.sin(p) for p in phi] for t in theta]
    xSphcap1 = [[Xcenter + R*np.cos(t) for p in phi] for t in theta]
    
    mlab.mesh(xSphcap1, ySphcap1, zSphcap1, color=COLOR, opacity=OPACITY)
    
    return 1




########################################################################
## MAIN PLOT
########################################################################

#==============================================================================
# Setup Figure
#==============================================================================
fig=mlab.figure(size=FIGSIZE, bgcolor=BG_COLOR)
visual.set_viewer(fig)



#==============================================================================
# Main Scene
#==============================================================================
ColorGrey = (0.2, 0.2, 0.2)



if plottype == 'geo':
    drawSphere(Xcenter=0, Ycenter=0, Zoffset=geo.T[2].max()/2 + step/4, R=98,
                    OPACITY=.10, COLOR=ColorGrey, Nphi=100, Nteta=100)
    mlab.points3d(*geo.T, scale_factor=22., color=ColorGrey, mode='cube')#, opacity=.25)


if plottype == 'geo_zoom':
    drawSphere(Xcenter=0, Ycenter=0, Zoffset=geo.T[2].max()/2 + step/4, R=98,
                    OPACITY=.06, COLOR=ColorGrey, Nphi=100, Nteta=100)
    mlab.points3d(*geo.T, scale_factor=22., color=ColorGrey, mode='cube', opacity=.05)
    geo_2 = geo[[96,46]]
    mlab.points3d(*geo_2.T, scale_factor=22., color=ColorGrey, mode='cube')#, opacity=.25)


if plottype in ['vectors', 'vectors_large']:
    
    if not plottype=='vectors_large':
        drawSphere(Xcenter=0, Ycenter=0, Zoffset=geo.T[2].max()/2 + step/4, R=98,
                        OPACITY=.05, COLOR=ColorGrey, Nphi=100, Nteta=100)
        mlab.points3d(*geo.T, scale_factor=1., color=ColorGrey, mode='cube')#, opacity=.25)
        mlab.points3d(*geo.T, scale_factor=22., color=ColorGrey, mode='cube', opacity=.05)
        NF = tools.get_field_as_list_by_fieldindex(sim, 0)
        mlab.quiver3d(*NF.T.real, color=(1,0,0), scale_factor=60, 
    #                  scale_mode="none",
                      mode='arrow')
    else:
        drawSphere(Xcenter=0, Ycenter=0, Zoffset=geo.T[2].max()/2 - step/4, R=90,
                        OPACITY=.05, COLOR=ColorGrey, Nphi=100, Nteta=100)
#        mlab.points3d(*geo2.T, scale_factor=1., color=ColorGrey, mode='cube')#, opacity=.25)
        mlab.points3d(*geo2.T, scale_factor=11., color=ColorGrey, mode='cube', opacity=.015)
        NF = tools.get_field_as_list_by_fieldindex(sim_large, 0)
        mlab.quiver3d(*NF.T.real, color=(1,0,0), scale_factor=20, 
    #                  scale_mode="none",
                      mode='arrow')

mlab.view(AZIM_SC, ELEV_SC, DIST_SC, FOCPT_SC)




if plottype == 'single_red_arrow':
    NF = np.array([[-100, 0, geo.T[2].max()/2 + step/4,
                   10, 0, 0]])
    mlab.quiver3d(*NF.T, color=(1,0,0), scale_factor=20,
                  mode='arrow')
    mlab.view(90, 45, DIST_SC, FOCPT_SC)
    
if plottype == 'single_blue_arrow':
    NF = np.array([[-100, 0, geo.T[2].max()/2 + step/4,
                   10, 0, 0]])
    mlab.quiver3d(*NF.T, color=(0,0,1), scale_factor=20,
                  mode='arrow')
    mlab.view(90, 45, DIST_SC, FOCPT_SC)






filename = '{}_cells.png'.format(savename_prefix)
mlab.savefig(filename, figure=fig)#, magnification=2)
mlab.show()














