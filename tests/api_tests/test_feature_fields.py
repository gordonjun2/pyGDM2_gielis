# encoding: utf-8
import os, sys
import pickle
 
import time

import numpy as np
import matplotlib.pyplot as plt
#~ from mayavi import mlab


from pyGDM2 import structures
from pyGDM2 import fields

from pyGDM2 import tools
from pyGDM2 import visu


from pyGDM import fields as fieldsOld








## ---------- Setup incident field
field_generator = fields.planewave
kwargs = dict(theta = 0.0, kSign=1)

field_generator = fields.gaussian
kwargs = dict(theta=0, xSpot=0, ySpot=0, zSpot=0, NA=-1, 
              spotsize=200, kSign=1, paraxial=0, returnField='E')

field_generator = fields.focused_planewave
kwargs = dict(theta=0, xSpot=0, ySpot=0, 
              NA=-1, spotsize=200, kSign=-1, returnField='E')

field_generator = fields.evanescent_planewave
kwargs = dict(theta_inc=44, polar='s', returnField='B')
kwargs = dict(theta_inc=44, polar='p', returnField='E')

#field_generator = fields.dipole_electric
#kwargs = dict(x0=0,y0=0,z0=100, mx=1,my=0,mz=0, returnField='E')
#kwargs = dict(x0=0,y0=0,z0=-100, mx=0,my=1,mz=0, returnField='H')

#field_generator = fields.dipole_magnetic
#kwargs = dict(x0=0,y0=0,z0=100, mx=1,my=0,mz=0, returnField='E')
#kwargs = dict(x0=0,y0=0,z0=100, mx=0,my=1,mz=0, returnField='E')
#kwargs = dict(x0=0,y0=0,z0=100, mx=0,my=0,mz=1, returnField='E')





n3 = 1.0
n2 = 1.0
n1 = 1.5 
spacing = 5000.0

n3 = 1.0
n2 = 1.05 + 1.8j 
n1 = 1.5
spacing = 40.0

wavelength = 500                     # one single wavelength
projection = 'XZ'
## --- 2D evaluation volume (plane)
r_probe = tools.generate_NF_map(-1000,1000,30, -1000,1000,30,0, projection=projection)
r_probe = tools.generate_NF_map(-500,500,30, -400,600,30,0, projection=projection)

## --- 3D evaluation volume
#r_probe = structures.rect_wire(70, L=30,H=30,W=20)
#r_probe.T[2] -= r_probe.T[2].max()/2.

NF = tools.evaluate_incident_field(field_generator, wavelength, kwargs, r_probe, 
                            n1=n1,n2=n2,n3=n3, spacing=spacing)











#%%
#==============================================================================
# visualization 
#==============================================================================
def plot_layer():
    fc_color_layer='b'
    if n1!=n2:
        if n2==n3: 
            plt.axhline(0, color='w'); plt.axhline(0, dashes=[2,2], color='k')
        else:
            plt.axhline(0, color='w',lw=1); plt.axhline(0, dashes=[2,2], color='k',lw=1)
            plt.axhspan(0,spacing, color='k', ls='--', fc=fc_color_layer, alpha=0.35)
            plt.axhline(spacing, color='w',lw=1); plt.axhline(spacing, dashes=[2,2], color='k',lw=1)
            


### -- 3 ways to plot the nearfield inside the structure:
plt.figure(figsize=(12,5))
plt.subplot(141)
v = visu.vectorfield(NF, complex_part='real', projection=projection, tit=projection+' real part', show=0)
plot_layer()

plt.subplot(142)
v = visu.vectorfield(NF, complex_part='imag', cmap=plt.cm.Reds, projection=projection, tit=projection+' imag part',  show=0)
plot_layer()




plt.subplot(243); plt.axis('off')
v = visu.vectorfield_color(NF, projection=projection, tit='intensity', show=0)
plt.colorbar(label=r'$|E|^2/|E0|^2$')
plot_layer()

plt.subplot(244); plt.axis('off')
v = visu.vectorfield_color(NF, projection=projection, fieldComp='ex', tit='Ex', show=0)
plt.colorbar(label=r'$E_x$')
plot_layer()

plt.subplot(247); plt.axis('off')
v = visu.vectorfield_color(NF, projection=projection, fieldComp='ey', tit='Ey', show=0)
plt.colorbar(label=r'$E_y$')
plot_layer()

plt.subplot(248); plt.axis('off')
v = visu.vectorfield_color(NF, projection=projection, fieldComp='ez', tit='Ez', show=0)
plt.colorbar(label=r'$E_z$')
plot_layer()




#from pyGDM2 import visu3d
#visu3d.animate_vectorfield(NF)
#
#
#ani = visu.animate_vectorfield(NF, projection='XZ')

#
#
### ---- OLD CODE FOR COMPARISON (case: dipolar emitter)
#if field_generator in [fields.dipole_electric, fields.dipole_magnetic]:
#    ## ---------- Setup dummy-structure for fundamental field-calculation
#    mesh = 'cube'
#    step = 30.0
#    geometry = structures.rect_wire(step, L=40,H=65,W=1, mesh=mesh)
#    geometry.T[2] -= 490
#    
#
#    SIMold = dict(elambda=[wavelength],
#                  atheta=[0],
#                  xm=geometry.T[0], ym=geometry.T[1], zm=geometry.T[2],
#                  n1 = complex(n1), n2 = complex(n2), n3 = complex(n2)
#                  )
#    
#    if field_generator == fields.dipole_electric:
#        Eold = fieldsOld.getE0elecDipole(SIMold, 1, 1, **kwargs)
#    else:
#        Eold = fieldsOld.getE0magDipole(SIMold, 1, 1, **kwargs)
#    
#    Eold = np.reshape(Eold, (len(Eold)/3, 3))
#    NFold = np.concatenate([geometry, Eold], axis=1)
#    
#    
#    print "NEW ", NF[0]
#    print "OLD", NFold[0]
#    plt.figure(figsize=(6,4))
#    plt.subplot(121)
#    v = visu.vectorfield(NFold, complex_part='real', projection='XZ', tit='OLD real', show=0)
#    
#    plt.subplot(122)
#    v = visu.vectorfield(NFold, complex_part='imag', cmap=plt.cm.Reds, projection='XZ', tit='OLD imag',  show=0)
#
#
#plt.show()
#


