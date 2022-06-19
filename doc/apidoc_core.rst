Core
*******************************

core
=========================

**core** contains the main *simulation* object and the *scatter* routine to calculate the self-consistent fields inside the nano-structure.


simulation description class
--------------------------------------

.. autoclass:: pyGDM2.core.simulation
   :members:
   :special-members: __init__

electric fields inside particle
-----------------------------------------------

.. autofunction:: pyGDM2.core.scatter

.. autofunction:: pyGDM2.core.scatter_mpi


decay rate of dipole transition
-----------------------------------------------

.. autofunction:: pyGDM2.core.decay_rate





Other functions
-----------------------

.. autofunction:: pyGDM2.core.get_side_by_side

.. autofunction:: pyGDM2.core.get_general_propagator

.. autofunction:: pyGDM2.core.get_efield_by_cg

