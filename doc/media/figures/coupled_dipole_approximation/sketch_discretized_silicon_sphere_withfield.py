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
FIGSIZE = (1200,1200)    # [width=0.8\textwidth]  -  (long) golden ratio: ~1.3 x 0.8 

BG_COLOR = (1,1,1)

## camera position
AZIM_SC=45; ELEV_SC=78; DIST_SC=600
FOCPT_SC=[0, 0, 90]



## suffix for filename

#plottype = 'geo'
#savename_prefix = 'nanosphere_discretized'
#
#plottype = 'geo_zoom'
#savename_prefix = 'nanosphere_discretized_zoom'

plottype = 'single_vector_with_emission'
savename_prefix = 'single_vector_with_emission'
double_dp = "single"
#double_dp = "double"
double_dp = "triple"
#double_dp = "all"


#plottype = 'single_red_arrow'
#savename_prefix = 'single_red_arrow'


########################################################################
## Generate data
########################################################################
1
## ---------- structure: nano-sphere
step = 25.
R = 3.0

mesh = 'cube'
geo = structures.sphere(step, R=R, mesh=mesh)
geo.T[2] += step/2
geo = structures.center_struct(geo)

## --------- simulation for quiverplot
material = materials.silicon()
n1, n2 = 1.0, 1.0

struct = structures.struct(step, geo, material, n1,n2, 
                                   structures.get_normalization(mesh))

#wavelengths = np.arange(650,770,20)
wavelengths = [5000]
field_generator = fields.planewave
kwargs = dict(theta = [135])
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)

sim = core.simulation(struct, efield)
core.scatter(sim)




## dipole emission field
r_probe = tools.generate_NF_map(-200,200,61, -200,200,61, 5, projection='XZ')
r_probe = structures.rotate_XY(r_probe.T, -45).T
#r_probe = r_probe.T[np.sqrt(r_probe[0]**2 +r_probe[1]**2 +r_probe[2]**2)>40].T


field_generator = fields.dipole_electric
x0_offset = 20 if double_dp=='all' else 0
kwargs = dict(x0=-0 + x0_offset, y0=0, z0=0, mx=1, my=-1, mz=0)
kwargs2 = dict(x0=65 + x0_offset, y0=-65, z0=30, mx=1, my=-1, mz=0)
kwargs3 = dict(x0=30 + x0_offset, y0=-30, z0=-60, mx=1, my=-1, mz=0)

r_probe = r_probe.T[np.sqrt((r_probe[0]-kwargs['x0']+6)**2 +
                   (r_probe[1]-kwargs['y0']-1)**2 +
                   (r_probe[2]-kwargs['z0'])**2)>25].T
if double_dp=="double" or double_dp=="triple" or double_dp=="all":
    r_probe = r_probe.T[np.sqrt((r_probe[0]-kwargs2['x0']+6)**2 +
                       (r_probe[1]-kwargs2['y0']-1)**2 +
                       (r_probe[2]-kwargs2['z0'])**2)>25].T
if double_dp=="triple" or double_dp=="all":
    r_probe = r_probe.T[np.sqrt((r_probe[0]-kwargs3['x0']+6)**2 +
                       (r_probe[1]-kwargs3['y0']-1)**2 +
                       (r_probe[2]-kwargs3['z0'])**2)>25].T
 
## --- layered environment
n1,n2,n3 = 1,1,1
spacing = 5000      # thickness of interface layer n2 (in nm)
## -- evaluation of the field-generator
wavelength = 70
NF_dp = tools.evaluate_incident_field(field_generator, wavelength, kwargs, r_probe,
                            n1=n1,n2=n2,n3=n3, spacing=spacing)

NF_dp = np.nan_to_num(NF_dp)
NF_dp.T[3:] /= 6.981276e-07
NF_dp.T[3:] *= np.sqrt((r_probe[0]-kwargs['x0'])**2 +
                   (r_probe[1]-kwargs['y0'])**2 +(r_probe[2]-kwargs['z0'])**2)**(.35)
if double_dp=="double" or double_dp=="triple" or double_dp=="all":
    NF_dp2 = tools.evaluate_incident_field(field_generator, wavelength, kwargs2, r_probe,
                                n1=n1,n2=n2,n3=n3, spacing=spacing)
    NF_dp2 = np.nan_to_num(NF_dp2)
    NF_dp2.T[3:] /= 6.981276e-07
    NF_dp2.T[3:] *= np.sqrt((r_probe[0]-kwargs2['x0'])**2 +
                   (r_probe[1]-kwargs2['y0'])**2 +(r_probe[2]-kwargs2['z0'])**2)**(.35)
    NF_dp.T[3:] = NF_dp.T[3:] + NF_dp2.T[3:]
        
if double_dp=="triple" or double_dp=="all":
    NF_dp3 = tools.evaluate_incident_field(field_generator, wavelength, kwargs3, r_probe,
                                n1=n1,n2=n2,n3=n3, spacing=spacing)
    NF_dp3 = np.nan_to_num(NF_dp3)
    NF_dp3.T[3:] /= 6.981276e-07
    NF_dp3.T[3:] *= np.sqrt((r_probe[0]-kwargs3['x0'])**2 +
                   (r_probe[1]-kwargs3['y0'])**2 +(r_probe[2]-kwargs3['z0'])**2)**(.35)
    NF_dp.T[3:] = NF_dp.T[3:] + NF_dp3.T[3:]
        
#%%
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



if plottype == 'single_vector_with_emission':
    drawSphere(Xcenter=0, Ycenter=0, Zoffset=geo.T[2].max()/2 + step/4, R=98,
                    OPACITY=.1, COLOR=ColorGrey, Nphi=100, Nteta=100)
    
    if double_dp=="triple":
        idxlist = [46, 96, 110]
    elif double_dp=="double":
        idxlist = [46, 96]
    elif double_dp=="single":
        idxlist = [46]
    elif double_dp=="all":
        idxlist = np.arange(len(geo))
    geo_2 = geo[idxlist]
    NF_2 = tools.get_field_as_list_by_fieldindex(sim, 0)[idxlist]
    mlab.points3d(*geo_2.T, scale_factor=1., color=ColorGrey, mode='cube')#, opacity=.25)
    mlab.points3d(*geo_2.T, scale_factor=22., color=ColorGrey, mode='cube', opacity=.05)
    sf_polar = 60 if double_dp != 'all' else 20
    mlab.quiver3d(*NF_2.T.real, color=(1,0,0), scale_factor=sf_polar, 
                  scale_mode="none",
                  mode='arrow')
    
    ## emission
    NF_dp.T[0] += geo_2[0][0]
    NF_dp.T[1] += geo_2[0][1]
    NF_dp.T[2] += geo_2[0][2]
#    visu3d.vectorfield_color(NF_dp, show=0)
#    mlab.flow(*np.nan_to_num(NF_dp).T.real)
    mlab.quiver3d(*NF_dp.T.real, scale_factor=.02, 
                  colormap='Blues',
#                  color=(0,0,1),
#                  scale_mode="none",
                  mode='arrow')
    
#    , color=(0,0,1), scale_factor=60, 
#                  scale_mode="none",
#                  mode='arrow')

mlab.view(AZIM_SC, ELEV_SC, DIST_SC, FOCPT_SC)






filename = '{}_{}.png'.format(savename_prefix, double_dp)
mlab.savefig(filename, figure=fig)#, magnification=2)
mlab.show()














