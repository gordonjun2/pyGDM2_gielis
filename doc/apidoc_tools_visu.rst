Visualization, tools
************************************


.. _label-visu:

visu
=========================

**visu** contains 2D visualization functions.

.. autofunction:: pyGDM2.visu.structure

.. autofunction:: pyGDM2.visu.structure_contour


vector fields
-----------------------

.. autofunction:: pyGDM2.visu.vectorfield

.. autofunction:: pyGDM2.visu.vectorfield_by_fieldindex

.. autofunction:: pyGDM2.visu.vectorfield_fieldlines



scalar fields
-----------------------

.. autofunction:: pyGDM2.visu.vectorfield_color

.. autofunction:: pyGDM2.visu.vectorfield_color_by_fieldindex

.. autofunction:: pyGDM2.visu.scalarfield

.. autofunction:: pyGDM2.visu.farfield_pattern_2D



*animate_vectorfield*
++++++++++++++++++++++++++++++++++++++

.. autofunction:: pyGDM2.visu.animate_vectorfield





visu3d
=========================

**visu3d** contains 3D visualization functions via `mayavi.mlab`. The API is mostly sytnax-compatible to the :ref:`label-visu` module.

.. autofunction:: pyGDM2.visu3d.structure


vector fields
-----------------------

.. autofunction:: pyGDM2.visu3d.vectorfield

.. autofunction:: pyGDM2.visu3d.vectorfield_by_fieldindex


scalar fields
-----------------------

.. autofunction:: pyGDM2.visu3d.vectorfield_color

.. autofunction:: pyGDM2.visu3d.vectorfield_color_by_fieldindex

.. autofunction:: pyGDM2.visu3d.scalarfield


*animate_vectorfield*
++++++++++++++++++++++++++++++++++++++

.. autofunction:: pyGDM2.visu3d.animate_vectorfield






tools
=========================

**tools** contains data-processing tools and helper functions


Save / Load / Print info
-------------------------------------

.. autofunction:: pyGDM2.tools.save_simulation

.. autofunction:: pyGDM2.tools.load_simulation

.. autofunction:: pyGDM2.tools.print_sim_info


Field-index searching
-------------------------------------

.. autofunction:: pyGDM2.tools.get_field_indices

.. autofunction:: pyGDM2.tools.get_closest_field_index


Structure geometry tools
-------------------------------------

.. autofunction:: pyGDM2.tools.test_geometry

.. autofunction:: pyGDM2.tools.get_geometry

.. autofunction:: pyGDM2.tools.get_step_from_geometry

.. autofunction:: pyGDM2.tools.get_geometry_2d_projection

.. autofunction:: pyGDM2.tools.get_geometric_cross_section

.. autofunction:: pyGDM2.tools.get_surface_meshpoints

.. autofunction:: pyGDM2.tools.get_closest_field_index


Spectra tools
-------------------------------------

.. autofunction:: pyGDM2.tools.get_possible_field_params_spectra

.. autofunction:: pyGDM2.tools.calculate_spectrum


Raster-scan tools
-------------------------------------

.. autofunction:: pyGDM2.tools.get_possible_field_params_rasterscan

.. autofunction:: pyGDM2.tools.get_rasterscan_fields

.. autofunction:: pyGDM2.tools.get_rasterscan_field_indices

.. autofunction:: pyGDM2.tools.calculate_rasterscan

.. autofunction:: pyGDM2.tools.get_possible_field_params_rasterscan

.. autofunction:: pyGDM2.tools.get_rasterscan_fields


Coordinate-list transformations
-------------------------------------

.. autofunction:: pyGDM2.tools.unique_rows

.. autofunction:: pyGDM2.tools.unique_rows_3D

.. autofunction:: pyGDM2.tools.adapt_map_to_structure_mesh

.. autofunction:: pyGDM2.tools.map_to_grid_2D

.. autofunction:: pyGDM2.tools.list_to_grid

.. autofunction:: pyGDM2.tools.grid_to_list

.. autofunction:: pyGDM2.tools.get_field_as_list

.. autofunction:: pyGDM2.tools.get_field_as_list_by_fieldindex

.. autofunction:: pyGDM2.tools.get_intensity_from_fieldlist


2D-mapping coordinate generators
-------------------------------------

.. autofunction:: pyGDM2.tools.generate_NF_map

.. autofunction:: pyGDM2.tools.generate_NF_map_XY

.. autofunction:: pyGDM2.tools.generate_NF_map_XZ

.. autofunction:: pyGDM2.tools.generate_NF_map_YZ


Other
-------------------------------------

.. autofunction:: pyGDM2.tools.evaluate_incident_field

.. autofunction:: pyGDM2.tools.get_intensity_from_fieldlist

.. autofunction:: pyGDM2.tools.combine_simulations





