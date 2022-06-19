.. raw:: html

    <video autoplay loop poster="_static/pygdmUI_static.png">
        <source src="_static/pygdmUI_static.mp4" type="video/mp4">
        <source src="_static/pygdmUI_static.webm" type="video/webm">
        Sorry, your browser doesn't support HTML5 video.
    </video>



Here we give an overview of the `pygdmUI <https://gitlab.com/wiechapeter/pygdm-ui>`_ package, an graphical user interface to pyGDM2.

*Note:* pygdmUI requires minimum pyGDM2 version v1.0.11 (experimental support). Stable support with pyGDM2 v1.1+.


.. contents:: Table of Contents



pygdm-UI
**************************

pyGDM-UI is a pyqt based, pure python graphical user interface to pyGDM2.

Requirements
================================

    - **python** (tested with 3.5+, `link <https://www.python.org/>`_)
    - **pyGDM2** (v1.1+, `link <https://wiechapeter.gitlab.io/pyGDM2-doc/>`_)
    - **pyqt5** (`link <https://pypi.org/project/PyQt5/>`_)
    - **matplotlib** (`link <https://matplotlib.org/>`_)
    - **mayavi** (`link <http://docs.enthought.com/mayavi/mayavi/mlab.html>`_)

Installation
=============================================

Via pip
-------------------------------

Install from pypi repository via

.. code-block:: bash
    
    $ pip install pygdmUI


Manual
-------------------------------

Download the latest code from `gitlab <https://gitlab.com/wiechapeter/pygdm-ui>`_ or clone the git repository:

.. code-block:: bash
    
    $ git clone https://gitlab.com/wiechapeter/pygdm-ui.git

    

Run
=============================================

Run in the terminal

.. code-block:: bash
    
    $ python3 -m pygdmUI.main
    
..     
..     
.. is based on the Green Dyadic Method (GDM) which resolves an optical Lippmann-Schwinger equation and **calculates the total field** :math:`\mathbf{E}(\mathbf{r}, \omega)` **, inside a nanostructure**, embedded in a fixed environment, upon illumination by an incident electromagnetic field :math:`\mathbf{E}_0(\mathbf{r}, \omega)`. The environment is described by the Green's tensor :math:`\mathbf{G}(\mathbf{r}, \mathbf{r'}, \omega)`:
.. 
.. .. math::
..     \mathbf{E}(\mathbf{r}, \omega)  = 
..      \mathbf{E}_0(\mathbf{r}, \omega) + 
..          \int \mathbf{G}(\mathbf{r}, \mathbf{r'}, \omega) \cdot 
..               \chi_{\text{e}} \cdot \mathbf{E}(\mathbf{r'}, \omega) \text{d} \mathbf{r'} 
.. 
.. 
.. 
.. Structure compositor
.. =============================================
.. .. 
.. 
.. In pyGDM, the discretization is done either on a cubic grid (center panel in image below), or on a hexagonal compact grid (right panel in image below).
.. 
.. .. figure:: _static/3D_volume_discretization_discretized.jpg
.. ..    :align: left
.. 
.. 
.. Simulation configurator
.. =============================================
.. 
.. One advantage of GDM compared to domain discretization techniques such as FDTD is that only the nanostructure is 
.. 
.. 
.. Simulation evaluation
.. =============================================
.. 
.. Spectra
.. -------------------------------
.. 
.. Spatial profiles / mappings
.. -------------------------------
.. 
.. Raster-scan simulations
.. -------------------------------
.. 

 

Author
=============================================

   - P\. R. Wiecha
