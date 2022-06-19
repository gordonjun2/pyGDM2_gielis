# encoding: utf-8
import matplotlib.pyplot as plt

from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core
from pyGDM2 import tools
from pyGDM2 import visu


scale = 1
img_file_in  = "quickstart_static.png"
mov_file_out = "quickstart"

img_file_in  = "installation_static.png"
mov_file_out = "installation"

#img_file_in  = "overview_static.png"
#mov_file_out = "overview"
#
#img_file_in  = "tutorials_static.png"
#mov_file_out = "tutorials"
#
#img_file_in  = "examples_static.png"
#mov_file_out = "examples"
#
#img_file_in  = "gallery_static.png"
#mov_file_out = "gallery"
#
img_file_in  = "pygdmUI_static.png"
mov_file_out = "pygdmUI_static"
#

#scale = 2
#img_file_in  = "API_static.png"
#mov_file_out = "API"

#===============================================================
# Setup the simulation
#===============================================================
## --- structure / environment
n1, n2 = 1.0, 1.0  # vacuum env.
step = 5.0
geometry = structures.image_to_struct(img_file_in, 
                      useDarkPixel=1, threshold=100, H=1,
                      nm_per_pixel=1.*step, stepsize=step)
#material = materials.gold()
material = materials.alu()
struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization('cube'))

## --- incident field
field_generator = fields.planewave
kwargs = dict(theta = [-45.0])

#field_generator = fields.dipole_electric
#kwargs = dict(x0=-350,y0=10,z0=-500, mx=1,my=0,mz=0)
#
#field_generator = fields.focused_planewave
#kwargs = dict(theta = [90.0], xSpot=-350, ySpot=10, spotsize=500)

wavelengths = [700]
efield = fields.efield(field_generator, wavelengths=wavelengths, 
                                                   kwargs=kwargs)

## --- simulation object
sim = core.simulation(struct, efield)
#visu.structure(sim, scale=0.5)
print("N dipoles {}".format(len(geometry)))

#%%
#===============================================================
# Run the simulation
#===============================================================
E = core.scatter(sim)
NF = tools.get_field_as_list_by_fieldindex(sim, 0)
    
#%%
#===============================================================
# create the field-animation
#===============================================================
# setup figure / axes
plt.figure(figsize=(6.0,2.5))
ax = plt.subplot(aspect='equal')
plt.axis('off')
plt.subplots_adjust(left=0, right=1, bottom=0,top=1)

## geometry
s = visu.structure(sim, scale=0.4*scale, color='.75', show=0)

## field-animation
config_vectorfield = dict(cmin=0.5, cmap=plt.cm.Blues, 
                          borders=50, vecwidth=0.8*scale)
ani = visu.animate_vectorfield(NF, Nframes=100, scale=15/scale,
                               kwargs=config_vectorfield, 
                               ax=ax, show=False)
plt.xlim(-550,550)
plt.ylim(-150,150)
# save video to file
#ani.save(mov_file_out+'.mp4', writer="ffmpeg", 
#         codec='h264', bitrate=1500)
ani.save(mov_file_out+'.webm', writer="ffmpeg", 
         codec='libvpx-vp9', bitrate=1500)

plt.show()

