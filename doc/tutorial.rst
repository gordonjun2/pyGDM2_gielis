..  _tutorials:
.. raw:: html

    <video autoplay loop poster="_static/tutorials_static.png">
        <source src="_static/tutorials.mp4" type="video/mp4">
        <source src="_static/tutorials.webm" type="video/webm">
        Sorry, your browser doesn't support HTML5 video.
    </video>



Tutorials
**************************

in contrast to the :ref:`examples<examples>`, which demonstrate pyGDM simulations on 
actual physical problems, this section contains more technical tutorials. 
*Note:* The syntax of the example and tutorial notebooks is python 2 only.

.. note:: You can download also :download:`the source code of the example and tutorial files <_downloads/examples.tar.gz>` (these are **python 2 and 3** compatible).


Use pyGDM 
=========================================

.. toctree:: 
   :maxdepth: 1
   :caption: Tutorials about the main pyGDM toolkit:

   tutorials/02_spectra.ipynb
   tutorials/03_field_maps.ipynb
   tutorials/04_rasterscans_tpl_heat.ipynb
   tutorials/05_rasterscans_ldos.ipynb
   tutorials/06_visualize_incident_fields.ipynb
   tutorials/07_MPI_spectra.ipynb
   tutorials/08_multi_materials_structures.ipynb
   
   
   
Extend pyGDM
===============================================

.. toctree:: 
   :maxdepth: 1
   :caption: Tutorials on how to extend the main pyGDM toolkit:

   tutorials/extend_01_own_structure_and_material.ipynb
   tutorials/extend_03_own_field_generator.ipynb
   


  

Use Evolutionary Optimization Toolkit
===============================================

.. note:: The tutorials on the **EO**-module are not deterministic.
          It is completely normal that the results differ from run to run.


.. toctree:: 
   :maxdepth: 1
   :caption: Tutorials about the evolutionary optimization interface in pyGDM:

   tutorials/EO_01_own_problem.ipynb
   tutorials/EO_02_multi_objective_problem.ipynb


