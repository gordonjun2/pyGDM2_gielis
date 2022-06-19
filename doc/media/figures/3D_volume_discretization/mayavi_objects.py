#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np


## Mayavi Imports
from mayavi import mlab
from mayavi.sources.builtin_surface import BuiltinSurface
from mayavi.modules.surface import Surface
from mayavi.filters.transform_data import TransformData
from tvtk.util import ctf
from tvtk.tools import visual




########################################################################
## CONSTANTS (eg. colors)
########################################################################


ColorGold = tuple(1.2 * np.array((207, 181, 59))/255.)
ColorRed = (1, 0, 0)
ColorGreen = (0, 1, 0)
ColorBlue = (0, 0, 1)
ColorGrey = (0.3, 0.3, 0.3)
ColorDarkGrey = (0.2, 0.2, 0.2)
ColorLightGrey = (0.6, 0.6, 0.6)



########################################################################
## PREDEFINED STRUCTURES AND VOLUME DISCRETIZATIONS
########################################################################


#==============================================================================
# Substrate
#==============================================================================
def drawRect(X0,X1, Y0,Y1, Z0,Z1, rotAx=(1,0,0), rotAng=(0),
                OPACITY=1, COLOR=(.1,.1,.1)):
    engine = mlab.get_engine()
    def rotMat3D(axis, angle, tol=1e-12):
        """Return the rotation matrix for 3D rotation by angle `angle` degrees about an
        arbitrary axis `axis`.
        """
        t = np.radians(angle)
        x, y, z = axis
        R = (np.cos(t))*np.eye(3) +\
            (1-np.cos(t))*np.matrix(((x**2,x*y,x*z),(x*y,y**2,y*z),(z*x,z*y,z**2))) + \
            np.sin(t)*np.matrix(((0,-z,y),(z,0,-x),(-y,x,0)))
        R[np.abs(R)<tol]=0.0
        return R

    # Main code
    # Add a cubic builtin source
    rect_src = BuiltinSurface()
    engine.add_source(rect_src)
    rect_src.source = 'cube'
    rect_src.data_source.center = np.array([ (X1+X0)/2.,  (Y1+Y0)/2.,  (Z1+Z0)/2.])
    rect_src.data_source.x_length=X1-X0
    rect_src.data_source.y_length=Y1-Y0
    rect_src.data_source.z_length=Z1-Z0
    #~ rect_src.data_source.capping = False
    #~ rect_src.data_source.resolution = 250


    # Add transformation filter to rotate cylinder about an axis
    transform_data_filter = TransformData()
    engine.add_filter(transform_data_filter, rect_src)
    Rt = np.eye(4)
    Rt[0:3,0:3] = rotMat3D(rotAx, rotAng) # in homogeneous coordinates
    Rtl = list(Rt.flatten()) # transform the rotation matrix into a list
    transform_data_filter.transform.matrix.__setstate__({'elements': Rtl})
    transform_data_filter.widget.set_transform(transform_data_filter.transform)
    transform_data_filter.filter.update()
    transform_data_filter.widget.enabled = False   # disable the rotation control further.

    # Add surface module to the cylinder source
    rect_surface = Surface()
    engine.add_filter(rect_surface, transform_data_filter)
    # add color property
    rect_surface.actor.property.color = COLOR
    rect_surface.actor.property.opacity = OPACITY
    return 1




#==============================================================================
# Nanorod
#==============================================================================
def drawNanorod(Xcenter=0, Ycenter=0, Zoffset=50, Lzyl=500, Rzylinder=200, RzylMin=200,
                OPACITY=1, COLOR=(.8,.8,.8), direction='X'):

    Lz     = np.linspace(Xcenter-Lzyl/2., Xcenter+Lzyl/2., 2)
    
    phi    = np.linspace(0, 2.*np.pi, 35) - 2
    theta  = np.linspace(0, 0.5*np.pi, 35)
    


    ## Zylindrical main Objective
    xZyl = [[ix for p in phi] for ix in Lz]
    yZyl = [[Ycenter + Rzylinder*np.sin(p) for p in phi] for iz in Lz]
    zZyl = [[Zoffset + Rzylinder*np.cos(p) for p in phi] for iz in Lz]
    
    ## Semi-sphere
    zSphcap1 = [[Zoffset + Rzylinder*np.sin(t)*np.cos(p) for p in phi] for t in theta]
    ySphcap1 = [[Ycenter + Rzylinder*np.sin(t)*np.sin(p) for p in phi] for t in theta]
    ySphcap2 = ySphcap1
    xSphcap1 = [[Lzyl/2. + Xcenter + Rzylinder*np.cos(t) for p in phi] for t in theta]
    xSphcap2 = [[Xcenter - Lzyl/2. - Rzylinder*np.cos(t) for p in phi] for t in theta]
    
    if direction=='Y':
        _x = yZyl
        yZyl = xZyl
        xZyl = _x
        _x = ySphcap1
        ySphcap1 = xSphcap1
        xSphcap1 = _x
        _x = ySphcap2
        ySphcap2 = xSphcap2
        xSphcap2 = _x
    
    
    mlab.mesh(xZyl, yZyl, zZyl, color=COLOR, opacity=OPACITY)
    mlab.mesh(xSphcap1, ySphcap1, zSphcap1, color=COLOR, opacity=OPACITY)
    mlab.mesh(xSphcap2, ySphcap2, zSphcap1, color=COLOR, opacity=OPACITY)
    
    return 1


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


#==============================================================================
# Triangle
#==============================================================================
def drawTriangle(Xcenter=0, Ycenter=0, Zoffset=0, aspect=1, stop=1, L=100, H=50,
                OPACITY=1, COLOR=(.8,.8,.8)):

    L1 = L*np.sqrt(3.)
    
    X = [[-L+Xcenter,L+Xcenter,(1-stop)*L+Xcenter,-(1-stop)*L+Xcenter,-L+Xcenter],
         [-L+Xcenter,L+Xcenter,(1-stop)*L+Xcenter,-(1-stop)*L+Xcenter,-L+Xcenter]]
    Y = [[Ycenter,Ycenter,L1+Ycenter,L1+Ycenter,Ycenter],
         [Ycenter,Ycenter,L1+Ycenter,L1+Ycenter,Ycenter]]
    Z = [[Zoffset,Zoffset,Zoffset,Zoffset,Zoffset],
         [H+Zoffset,H+Zoffset,H+Zoffset,H+Zoffset,H+Zoffset]]
    Y = np.array(Y)*aspect
    
    Xcap = [[-L+Xcenter,L+Xcenter],[-(1-stop)*L+Xcenter,(1-stop)*L+Xcenter]]
    Ycap = [[Ycenter,Ycenter],[L1+Ycenter,L1+Ycenter]]
    Zcap1 = np.ones(np.shape(Xcap))*0 + Zoffset
    Zcap2 = np.ones(np.shape(Xcap))*H + Zoffset
    Ycap = np.array(Ycap)*aspect
    
    
    mlab.mesh(X,Y,Z, color=COLOR, opacity=OPACITY)
    mlab.mesh(Xcap,Ycap,Zcap1, color=COLOR, opacity=OPACITY)
    mlab.mesh(Xcap,Ycap,Zcap2, color=COLOR, opacity=OPACITY)

    return 1




#==============================================================================
# Arrows
#==============================================================================
def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2, Lcone=0.25,Rcone=0.08, Rshaft=0.03, color=(1,0,0)):
    ar1=visual.arrow(x=x1, y=y1, z=z1, color=color)
    ar1.length_cone=Lcone
    
    ar1.radius_cone = Rcone
    ar1.radius_shaft = Rshaft
    
    arrow_length=np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.actor.scale=[arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos/arrow_length
    ar1.axis = [x2-x1, y2-y1, z2-z1]
    
    
    return ar1