..  _examples:
.. raw:: html

    <video autoplay loop poster="_static/examples_static.png">
        <source src="_static/examples.mp4" type="video/mp4">
        <source src="_static/examples.webm" type="video/webm">
        Sorry, your browser doesn't support HTML5 video.
    </video>


    
Examples
**************************

This section contains examples in which pyGDM is used on specific physical problems. 
For more technical examples, see the :ref:`tutorials<tutorials>`.
*Note:* The syntax of the example and tutorial notebooks is python 2 only.

.. note:: You can download also :download:`the source code of the example and tutorial files <_downloads/examples.tar.gz>` (these are **python 2 and 3** compatible).
             


Comparison to Mie theory
=========================================

.. note:: The spectra calculated by Mie theory used in the corresponding 
             examples can be :download:`downloaded here <_downloads/mie_example_spectra.tar.gz>`.


.. toctree:: 
   :maxdepth: 1

   examples/example01_mie01.ipynb
   examples/example02_mie02.ipynb
   examples/example03_mie03.ipynb



Comparison to selected publications
=========================================

.. toctree:: 
   :maxdepth: 1

   examples/example04_SiSphere_FWBW.ipynb
   examples/example05_splitring_with_dipole.ipynb
   examples/example06_polarization_conversion.ipynb
   examples/example07_heat_generation.ipynb
   examples/example08_decay_rate.ipynb
   examples/example09_multipole_decomposition.ipynb



Evolutionary optimization
=========================================

.. note:: The examples on the evolutionary optimization module (**EO**) are not deterministic.
             If the results differ slightly from run to run, this is completely normal and actually the 
             consequence of the main drawback of evolutionary optimization: 
             *Convergence to the global optimum can never be guaranteed!*. 
             Nevertheless, usually the optimizations converge very closely to the optimum, which
             can be tested by repeated runs of the same optimization.


.. toctree:: 
   :maxdepth: 1

   examples/exampleEO_01_scattering.ipynb
   examples/exampleEO_02_nearfield.ipynb
