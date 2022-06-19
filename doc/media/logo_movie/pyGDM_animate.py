# encoding: utf-8
import matplotlib.pyplot as plt

from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core
from pyGDM2 import tools
from pyGDM2 import visu


#===============================================================
# Setup the simulation
#===============================================================
## --- structure / environment
n1, n2 = 1.0, 1.0  # vacuum env.
step = 8.0
geometry = structures.image_to_struct("pyGDM_logo_static.png", 
                      useDarkPixel=1, threshold=100, H=1,
                      nm_per_pixel=1.*step, stepsize=step)
material = materials.gold()
struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization('cube'))

## --- incident field
field_generator = fields.planewave
kwargs = dict(theta = [45.0])
wavelengths = [700]
efield = fields.efield(field_generator, wavelengths=wavelengths, 
                                                   kwargs=kwargs)

## --- simulation object
sim = core.simulation(struct, efield)
visu.structure(sim)

#===============================================================
# Run the simulation
#===============================================================
E = core.scatter(sim)
NF = tools.get_field_as_list_by_fieldindex(sim, 0)
    

#===============================================================
# create the field-animation
#===============================================================
## setup figure / axes
plt.figure(figsize=(6.0,2.5))
ax = plt.subplot(aspect='equal')
plt.axis('off')
plt.subplots_adjust(left=0, right=1, bottom=0,top=1)

## geometry
s = visu.structure(sim, scale=0.1, color='.75', show=0)

## field-animation
config_vectorfield = dict(cmin=0.5, cmap=plt.cm.Blues, 
                          borders=50, vecwidth=0.8)
ani = visu.animate_vectorfield(NF, Nframes=100, scale=12,
                               kwargs=config_vectorfield, 
                               ax=ax, show=False)
# save video to file
ani.save('pyGDM_logo.mp4', writer="ffmpeg", 
         codec='h264', bitrate=1500)
ani.save('pyGDM_logo.webm', writer="ffmpeg", 
         codec='libvpx', bitrate=1500)
plt.show()