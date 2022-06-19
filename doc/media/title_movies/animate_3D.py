# encoding: utf-8
import matplotlib.pyplot as plt

from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields

from pyGDM2 import core
from pyGDM2 import tools
from pyGDM2 import visu



scale = 2
img_file_in  = "3D_static.png"
mov_file_out = "3D"

#===============================================================
# Setup the simulation
#===============================================================
## --- structure / environment
n1, n2 = 1.0, 1.0  # vacuum env.
step = 5.0
geometry = structures.image_to_struct(img_file_in, 
                      useDarkPixel=1, threshold=100, H=3,
                      nm_per_pixel=0.75*step, stepsize=step)
#material = materials.gold()
material = materials.alu()
struct = structures.struct(step, geometry, material, n1,n2, 
                           structures.get_normalization('cube'))

## --- incident field
field_generator = fields.planewave
kwargs = dict(theta = [-45.0])

#field_generator = fields.dipole_electric
#kwargs = dict(x0=-0,y0=0,z0=50, mx=1,my=1,mz=0)
#
#field_generator = fields.focused_planewave
#kwargs = dict(theta = [90.0], xSpot=-350, ySpot=10, spotsize=500)

wavelengths = [1000]
efield = fields.efield(field_generator, wavelengths=wavelengths, 
                                                   kwargs=kwargs)

## --- simulation object
sim = core.simulation(struct, efield)
visu.structure(sim, scale=0.5)
print "N dipoles {}".format(len(geometry))

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
from pyGDM2 import visu3d
from mayavi import mlab

fig = mlab.figure( size=(600, 300), bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.,0.,0.) )

## structure
visu3d.structure(sim, axis_labels=False, draw_substrate=False, 
                 opacity=0.1, show=False)


## 3D field-animation
ani2 = visu3d.animate_vectorfield(NF, Nframes=100, scale=8, 
                          draw_struct=False, 
                          draw_substrate=False, substrate_size=1.1, 
                          colormap='Blues', clim=[0.0, 0.5], 
                          fig=fig, view=(85, -45, 350, (0,0,-15)),
#                          ffmpeg_args="-b:v 1.5M -c:v libx264", mov_file="3D.mp4",
                          ffmpeg_args="-b:v 1.5M -c:v libvpx", mov_file="3D.webm",
                          save_anim=True, 
                          opacity=0.5)

