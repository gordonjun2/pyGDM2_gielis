Evolutionary Optimization
************************************

The evolutionary optimization module **EO** provides an interface between pyGDM and the `pagmo/pygmo <https://esa.github.io/pagmo2/>`_.
**pygmo** is capable to do massively parallelized optimizations via what their authors call the "generalized island model", however, due to the recent transition from pygmo 1.x to version 2.x, this parallelization technique is not yet supported in pyGDM.
A pyGDM-interface to pagmo/pygmo 1.x exists, with support for the generalized island model. 
We can provide this older interface upon request.


EO.core
=========================

`core` contains the main optimization routine.

*run_eo*: Launch an optimization
-------------------------------------

.. autofunction:: pyGDM2.EO.core.run_eo

*continue_eo*: continue a former optimization
----------------------------------------------

.. autofunction:: pyGDM2.EO.core.continue_eo



EO.models
=========================

`models` contains the structure models (geometries to optimize).


*BaseModel*: Base class for model definitions
----------------------------------------------------

.. autoclass:: pyGDM2.EO.models.BaseModel
   :members:


Predefined models
----------------------------------------------------

.. autoclass:: pyGDM2.EO.models.RectangularAntenna

.. autoclass:: pyGDM2.EO.models.CrossAntenna

.. autoclass:: pyGDM2.EO.models.MultiRectAntenna

.. autoclass:: pyGDM2.EO.models.BlockModel





EO.problems
=========================

`problems` contains the problems defining the optimization target(s).

*BaseProblem*: Base class for problem definitions
----------------------------------------------------

.. autoclass:: pyGDM2.EO.problems.BaseProblem
   :members:

Predefined problems
----------------------------------------------------

.. autoclass:: pyGDM2.EO.problems.ProblemScat

.. autoclass:: pyGDM2.EO.problems.ProblemNearfield

.. autoclass:: pyGDM2.EO.problems.ProblemDirectivity





EO.tools
=========================

`tools` contains various EO helper functions.

Get results from stored optimization
-----------------------------------------------------------------------

.. autofunction:: pyGDM2.EO.tools.reload_eo

.. autofunction:: pyGDM2.EO.tools.get_best_candidate

.. autofunction:: pyGDM2.EO.tools.get_best_candidate_f_x

.. autofunction:: pyGDM2.EO.tools.get_model

.. autofunction:: pyGDM2.EO.tools.get_problem

.. autofunction:: pyGDM2.EO.tools.get_population


Results from Pareto multi-objective optimization
-------------------------------------------------------------------------------

.. autofunction:: pyGDM2.EO.tools.get_pareto_fronts

.. autofunction:: pyGDM2.EO.tools.plot_pareto_2d

.. autofunction:: pyGDM2.EO.tools.plot_all_pareto_fronts_2d


Other tools related to the Evolutionary Optimization module
---------------------------------------------------------------------------------------

.. autofunction:: pyGDM2.EO.tools.calculate_solid_angle_by_dir_index


