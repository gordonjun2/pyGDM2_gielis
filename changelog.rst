Changelog
******************

[v1.1] (unreleased)
=====================
breaking
--------------
- planned: change in the field generator API to be more consistent within the "efield" class and to facilitate field evaluation. This will affect the usage of field generators outside pyGDM simulations

added
--------------
- full support of pyGDM-UI GUI
- full support of callback functions
- unit tests for all routines

changes
--------------
- full implementation of the API in pure python



[v1.0.11] - 2020-02-25
=====================
breaking
--------------
- internal API changes in `core`: re-structured `scatter` and `get_generalized_propagator`. The latter now takes a `sim` instance as input. The order of kwargs was changed.

added
--------------
- added a geometry consistency check to structure class
- new module `linear_py` with experimental pure python implementations of all linear functions
- optional radiative correction term in "linear_py.extinct" (can improve absorption section with large stepsizes)
- `visu`: 2D plotting functions try to determine the best 2D projection automatically
- started writing unit tests
- callback function support for `core.scatter` 

changes
--------------
- conversion to pure python of some helper functions
- some internal modifications for pygdmUI GUI support

fixes
--------------
- fixed geometry consistency-test routine for multi-material structures
- fixes in autoscaling in `visu.structure` (TODO: adapt to screen dpi)
- fixed bug in "linear.farfield" in cases when n2>n1 ("environment" optically denser than "substrate")
- **potentially breaking!!**: fixed several structures, where "hex" meshing gave a wrong height. *Attention*, following structure generators may now produce different heights: `lshape_rect`, `lshape_rect_nonsym`, `lshape_round`, `split_ring`, `rect_split_ring`
- some code cleaning
- minor fixes in several visualization tools



[v1.0.10.1] - 2019-10-08
=====================
added
--------------
- `materials.hyperdopedFromFile`: materials class, which adds doping to any tabulated dielectric permittivity (contributed by C. Majorel)

fixes
--------------
- fixed bug in "linear.farfield", causing zero scattering at angles teta > 3pi / 2



[v1.0.10] - 2019-10-02
=====================
added
--------------
- `linear.optical_chirality`: chirality of the electromagnetic field (experimental feature)
- new structure generator `polygon`
- `tools.combine_simulations`: tool to combine the structures of several simulations into a single simulation. Combining simulations with calculated fields, this also allows to analyze how structures behave if optical interactions are artificially deactivated.
- added support for "cupy" GPU solver (req. version 7+) as alternative to "pycuda"
- added experimental pure-python implementation of propagators and coupling matrix initialization routines based on "numba" (by default pygdm is still using the former fortran implementation)
- added experimental support for tensorial material permittivity

fixes
--------------
- critical fix in "linear.extinct": Works correctly now for environments n2!=1
- corrected phase of B-field in field generators "planewave" and "gaussian"
- fieldindex-search: works now with strings as fieldgenerator kwargs
- added exception handling to "linear.farfield" for simulation configurations where the underlying approximations don't apply



[v1.0.9] - 2019-08-22
=====================
no more compiled binaries for python 2.7 (compilation from source still possible)

fixes
--------------
- critical fix in linear.farfield: works correctly now also for non-vacuum environment above the substrate (refractive index n2 != 1)



[v1.0.8] - 2019-06-07
=====================
added
--------------
- multipole decomposition (dipole and quadrupole moments)
- elliptic polarization in field generators "planewave", "focused" and "gaussian"
- new materials: *hyperdopedSilicon* and *hyperdopedConstantDielectric* (contributed by C. Majorel)
- extended capabilities for "visu3d.animate_vectorfield" and according documentation
- zero-field generator

fixes
--------------
- linear.farfield: scattering into a substrate now correctly calculated (contributed by C. Majorel)
- python3 compatibility: fixed structure generator problem with hexagonal meshing and some float parameters. Also fixed the python 3.X compatibility of the examplescripts.
- fixed a bug in silver dispersion
- numerous small fixes and docstring improvements



[v1.0.7] - 2018-11-20
=====================
added
--------------
- experimental CUDA support for matrix inversion on GPU (method "cuda")
- structure generators:
    - "prism" now supports truncated edges
    - "spheroid"

fixes
--------------
- MAJOR: fix absolute import error in "visu3d"module, which was broken in former version
- minor fix in struct class, treats lists of wavelengths correctly now (was not affecting pyGDM itself. Failed if a `struct` instance was externally used with a list of wavelengths)



[v1.0.6] - 2018-10-31
=====================
added
--------------
- compatibility with python3 (compatible with python 2 and 3 now)
- default inversion method is now in-place LU decomposition: reduces memory requirement by ~40%
- added some tools to simplify intensity calculation

fixes
--------------
- fix in visu.animate: Works now glitch-less with any user-defined framerate
- minor fix: all classes now initialize with single precision by default. 



[v1.0.5] - 2018-07-9
=====================
fixes
--------------
- critical fix in hdf5 saving / loading. hdf5-data was corrupted during saving/reloading. Works now.

minor
--------------
- by default, multithreading disabled in MPI-version of "scatter". Using SLURM, MPI and pathos seems to conflict which results in major performance drop



[v1.0.4] - 2018-06-07
=====================
added
--------------
- multi-threading support via "thanos" in generalized propagator operations. 
  This drastically increases the speed of raster-scan simulations on multi-core systems.
- hdf5 support for saving/loading simulations
    - doubles the speed for writing, triples speed for reading
    - by default, using "blosc" compression, reduces the filesize by ~ 50%
- hexagonal meshing support in "image_to_struct"
- support for scipy < V0.17 in "decay"



[v1.0.3] - 2018-04-06
=====================
added
--------------
- intallation instructions for MacOSX



[v1.0.2] - 2018-03-29
=====================
added
--------------
- "visu.structure" does automatic multi-structure plots
- compile option for compilation without openmp
- several structure models
- hardcoded silver dielectric function

fixes
--------------
- in "visu.vectorfield_color", fixed an error in the calculation of the field intensity



[v1.0.1] - 2018-02-13
=====================
fixes
--------------
- fixes in "setup.py" script
