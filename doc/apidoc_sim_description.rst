Simulation description
*******************************



structures
=========================

**structures** contains the main *struct* object, structure generators and helper functions related to the generation of structure geometries.

structure definition class *struct*
--------------------------------------------

.. autoclass:: pyGDM2.structures.struct
   :members:
   :special-members: __init__
   

Structure generators
-----------------------

.. autofunction:: pyGDM2.structures.rect_wire

.. autofunction:: pyGDM2.structures.rect_dimer

.. autofunction:: pyGDM2.structures.prism

.. autofunction:: pyGDM2.structures.rhombus

.. autofunction:: pyGDM2.structures.hexagon

.. autofunction:: pyGDM2.structures.double_hexagon

.. autofunction:: pyGDM2.structures.diabolo

.. autofunction:: pyGDM2.structures.polygon

.. autofunction:: pyGDM2.structures.sphere

.. autofunction:: pyGDM2.structures.spheroid

.. autofunction:: pyGDM2.structures.nanodisc

.. autofunction:: pyGDM2.structures.nanorod

.. autofunction:: pyGDM2.structures.lshape_rect

.. autofunction:: pyGDM2.structures.lshape_rect_nonsym

.. autofunction:: pyGDM2.structures.lshape_round

.. autofunction:: pyGDM2.structures.split_ring

.. autofunction:: pyGDM2.structures.rect_split_ring



Mesher
-----------------------

.. autofunction:: pyGDM2.structures._meshCubic

.. autofunction:: pyGDM2.structures._meshHexagonalCompact


Other functions
-----------------------

.. autofunction:: pyGDM2.structures.image_to_struct

.. autofunction:: pyGDM2.structures.struct_to_image

.. autofunction:: pyGDM2.structures.get_normalization

.. autofunction:: pyGDM2.structures.center_struct

.. autofunction:: pyGDM2.structures.rotate_XY


materials
=========================

**materials** is a collection of classes defining the dispersion of materials.


*fromFile* class
-----------------------

.. autoclass:: pyGDM2.materials.fromFile
   :members:
   :special-members: __init__

Predefined materials
-----------------------

.. autoclass:: pyGDM2.materials.dummy

.. autoclass:: pyGDM2.materials.gold

.. autoclass:: pyGDM2.materials.alu

.. autoclass:: pyGDM2.materials.silver

.. autoclass:: pyGDM2.materials.silicon

.. autoclass:: pyGDM2.materials.hyperdopedConstantDielectric

.. autoclass:: pyGDM2.materials.hyperdopedFromFile

.. autoclass:: pyGDM2.materials.hyperdopedSilicon






fields
=========================

**fields** contains the main *efield* object and field-generator functions for frequently used fundamental fields.

Incident electro-magnetic field class *efield*
--------------------------------------------------------------------

.. autoclass:: pyGDM2.fields.efield
   :members:
   :special-members: __init__


field generators
-----------------------

.. autofunction:: pyGDM2.fields.planewave

.. autofunction:: pyGDM2.fields.evanescent_planewave

.. autofunction:: pyGDM2.fields.focused_planewave

.. autofunction:: pyGDM2.fields.gaussian

.. autofunction:: pyGDM2.fields.double_gaussian

.. autofunction:: pyGDM2.fields.dipole_electric

.. autofunction:: pyGDM2.fields.dipole_magnetic



