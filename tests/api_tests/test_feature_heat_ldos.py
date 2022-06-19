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
from pyGDM2 import nonlinear
from pyGDM2 import visu
from pyGDM2 import tools





## ---------- Setup structure
mesh = 'cube'
mesh = 'hex'
step = 30.0
geometry = structures.rect_wire(step, L=500/step-1,H=1, W=500/step, mesh=mesh)
geometry = structures.hexagon(step, NSIDE=10, H=2, mesh=mesh)
geometry.T[2] += step

material = materials.dummy(2.0)
material = materials.silicon()
material = materials.gold()
material = materials.alu()
n1, n2 = 1.5, 1.33  # constant environment

struct = structures.struct(step, geometry, material, n1,n2, structures.get_normalization(mesh))
struct = structures.center_struct(struct)



## ---------- Setup incident field
field_generator = fields.focused_planewave        # planwave excitation
#field_generator = fields.gaussian        # planwave excitation
#kwargs = dict(theta = [0.0, 45.0, 90.0], kSign=[-1,1])              # several configurations
#wavelengths = np.linspace(400,1000,25)                     # spectrum
xSpot = np.linspace(-500, 500, 21)
ySpot = np.linspace(-500, 500, 21)
kwargs = dict(theta = [ 0], spotsize=[50], kSign=-1, xSpot=xSpot, ySpot=ySpot)
wavelengths = [750]                     # one single wavelength
efield = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs)



## ---------- Simulation initialization
sim = core.simulation(struct, efield)

visu.structure(sim.struct.geometry, scale=0.5)
print "N dipoles:", len(sim.struct.geometry)



rasterscan_fieldconfigs = tools.get_possible_field_params_rasterscan(sim)

print '\n\navailable rasterscan configurations' 
for p in rasterscan_fieldconfigs:
    print p



#%%
########################################################################
##############            Test Simulations               ###############
########################################################################
core.scatter(sim, method='lu', verbose=True)
#tools.save_simulation(sim, "rasterscan_sim.pcl")
#sim = tools.load_simulation("rasterscan_sim.pcl")

#%%
rasterscan_field_params = rasterscan_fieldconfigs[0]
TPL_xy = []
Q_xy   = []
DT_xy  = []
DT_scan_probe_xy = []




Zprobe = sim.struct.geometry.T[2].max()+sim.struct.step/2.
r_probe = (0, 0, 150 + Zprobe)
#r_probe = (0, 0, 150)

## --- calculate TPL signal via indices at each position
key_x_pos='xSpot'; key_y_pos='ySpot'
E_rasterscan_indices = tools.get_rasterscan_field_indices(sim, rasterscan_field_params)

t0 = time.time()
for i in E_rasterscan_indices:
    params = sim.E[i][0]
    Itpl = nonlinear.tpl_ldos(sim, i)
    TPL_xy.append([params[key_x_pos], params[key_y_pos], Itpl])
print "time TPL        : {:.2}s".format(time.time()-t0)
t0 = time.time()
for i in E_rasterscan_indices:
    params = sim.E[i][0]
    Q = linear.heat(sim, i)
    Q_xy.append([params[key_x_pos], params[key_y_pos], Q])
print "time heat       : {:.2}s".format(time.time()-t0)
t0 = time.time()
for i in E_rasterscan_indices:
    params = sim.E[i][0]
    DT = linear.temperature(sim, i, r_probe, kappa_env=0.6, kappa_subst=0.8)
    DT_xy.append([params[key_x_pos], params[key_y_pos], DT])
print "time temperature: {:.2}s".format(time.time()-t0)
t0 = time.time()
TPL_xy = np.array(TPL_xy)
Q_xy = np.array(Q_xy)
DT_xy = np.array(DT_xy)


## --- do a temperature probe rasterscan with fixed incident beam position
#xFixed, yFixed = 300, 300
#idx_E = tools.get_closest_field_index(sim, dict(xSpot=xFixed, ySpot=yFixed))
#xFixed, yFixed = sim.E[idx_E][0]["xSpot"], sim.E[idx_E][0]["ySpot"]
#scan_pos = tools.generate_NF_map_XY(xSpot.min(), xSpot.max(), len(xSpot),
#                                    ySpot.min(), ySpot.max(), len(ySpot),
#                                    Z0=Zprobe).T
#for r_prb in scan_pos:
#    DT = linear.temperature(sim, idx_E, r_prb)
#    DT_scan_probe_xy.append([r_prb[0],r_prb[1], DT])
#print "time temp. scan : {:.2}s".format(time.time()-t0)
#t0 = time.time()
#DT_scan_probe_xy = np.array(DT_scan_probe_xy)


#%%
plt.figure(figsize=(15,5))

if len(TPL_xy)!=0:
    TPL_map, extent = tools.list_to_grid(TPL_xy)

    plt.subplot(141, aspect='equal')
    plt.title("TPL")
    plt.imshow(TPL_map, extent=extent, interpolation='none')
    plt.xlabel("Y (nm)"); plt.ylabel("Y (nm)")
    plt.colorbar(orientation='horizontal', label='intensity (a.u.)')

if len(Q_xy)!=0:
    Q_map, extent = tools.list_to_grid(Q_xy)

    plt.subplot(142, aspect='equal')
    plt.title("total deposited Q")
    plt.imshow(Q_map/1E3, extent=extent, interpolation='none')
    plt.xlabel("Y (nm)"); plt.ylabel("Y (nm)")
    plt.colorbar(orientation='horizontal', label='heat (micro watt)')

if len(DT_xy)!=0:
    DT_map, extent = tools.list_to_grid(DT_xy)
    
    plt.subplot(143, aspect='equal')
    plt.title(r"$\Delta T$ at {}".format(r_probe))
    plt.imshow(DT_map, extent=extent, interpolation='none')
    plt.xlabel("Y (nm)"); plt.ylabel("Y (nm)")
    plt.colorbar(orientation='horizontal', label=r'$\Delta$T (Kelvin)')

if len(DT_scan_probe_xy)!=0:
    DT_scan_map, extent = tools.list_to_grid(DT_scan_probe_xy)

    plt.subplot(144, aspect='equal')
    plt.title(r"$\Delta T$ scan, focused beam fixed at {}".format( (xFixed, yFixed) ))
    plt.imshow(DT_scan_map, extent=extent, interpolation='none')
    plt.xlabel("Y (nm)"); plt.ylabel("Y (nm)")
    plt.colorbar(orientation='horizontal', label=r'$\Delta$T (Kelvin)')



plt.show()


