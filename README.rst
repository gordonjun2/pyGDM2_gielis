***********************************
Requirements / Installation
***********************************

pyGDM2 is available on `pypi <https://pypi.python.org/pypi/pygdm2/>`_ and `gitlab <https://gitlab.com/wiechapeter/pyGDM2>`_. 

Detailed documentation with many examples is avaiable at the `pyGDM2 documentation website <https://wiechapeter.gitlab.io/pyGDM2-doc/>`_. See also the `documentation paper on arXiv (1802.04071) <https://arxiv.org/abs/1802.04071>`_




Requirements
================================

Python
------------------
    - **python** (2.7 or 3.5+, `python <https://www.python.org/>`_)
    - **numpy** (`numpy <http://www.numpy.org/>`_)
    - **python headers** (under ubuntu, install the package *python-dev* or *python-devel*)


Fortran
------------------
    - *fortran* compiler (tested with **gfortran**. `gcc <https://gcc.gnu.org/fortran/>`_)
    - **openmp** (usually comes with fortran. `openmp <http://www.openmp.org/>`_)
    - **f2py** (comes with **numpy**. `link <http://www.numpy.org/>`_)


Optional Python packages
-------------------------------------
    - **scipy** >= v0.17.0, lower versions supported with restrictions (*Strongly recommended*. Used for standard solver LU decomposition and several tools. `scipy <https://www.scipy.org/>`_)
    - **matplotlib** (*Strongly recommended*. For all 2D visualization tools. `matplotlib <https://matplotlib.org/>`_)
    - **mpi4py** (for MPI parallelized calculation of spectra. `mpi4py <http://mpi4py.readthedocs.io/en/stable/>`_)
    - **mayavi** (for all 3D visualization. `mayavi <http://docs.enthought.com/mayavi/mayavi/mlab.html>`_)
    - **PIL** (image processing. `PIL <https://pypi.python.org/pypi/PIL>`_)
    - **pathos** (for multi-threaded generalized propagator operations. `pathos <https://pypi.org/project/pathos/>`_)
    - **pytables** (v3.x recommended. For hdf5 saving/loading of simulations. `pytables <https://www.pytables.org/>`_)
    - **PaGMO / PyGMO** (version 2.4+. *Required* for the **EO** submodule. `pagmo <https://esa.github.io/pagmo2/>`_)
    - **pycuda** (tested with version 2018.1, for GPU-based matrix inversion. `pyCUDA <https://documen.tician.de/pycuda/>`_, for problems during installation with pip, see `solution proposed here <https://codeyarns.com/2015/07/31/pip-install-error-with-pycuda/>`_)
    - **scikit-cuda** (tested with version 0.5, for GPU-based matrix inversion. `scikit-cuda <https://scikit-cuda.readthedocs.io/en/latest/>`_)
    - alternatively: **cupy** (version 7+, tested with version 7.0.0b4, alternative to *pycuda* for GPU-based matrix inversion. `cupy <https://docs-cupy.chainer.org/en/stable/index.html>`_)

(all available via `pip <https://pypi.python.org/pypi/pip>`_)



Installation under linux
=============================================

Via pip
-------------------------------

Install from pypi repository via

.. code-block:: bash
    
    $ pip install pygdm2



Via setup script
-------------------------------

The easiest possibility to compile (and install) pyGDM is via the 
setup-script, which uses the extended *distutils* from *numpy*. 

To install pyGDM, run in the source directory:

.. code-block:: bash
    
    $ python setup.py install

To install to a user-defined location, use the *prefix* option:

.. code-block:: bash
    
    $ python setup.py install --prefix=/some/specific/location


To only compile without installation, use

.. code-block:: bash
    
    $ python setup.py build




Manual compilation
-------------------------------------------------------------

1. clone git:

   .. code-block:: bash
    
        $ git clone https://gitlab.com/wiechapeter/pyGDM2.git

2.a *python 2.7*: compile fortran parts:

   .. code-block:: bash
    
        $ cd fortranBase
        $ make
        
2.b *python 3.5+*:

   .. code-block:: bash
    
        $ cd fortranBase
        $ make python3

3. *optional, for system-wide usage* add to **path** and **pythonpath**, 
   e.g. add following lines to file "/home/USER/.profile", where 
   "path_of_pyGDM_folder" is the pyGDM installation directory:
  
   .. code-block:: bash
    
        PATH="path_of_pyGDM_folder:$PATH"
        export PATH
        
        PYTHONPATH="path_of_pyGDM_folder:$PYTHONPATH"
        export PYTHONPATH

        


Installation under windows
=============================================

For windows, we also recommend `Anaconda <https://www.anaconda.com/download/#windows>`_ in which pyGDM can be installed easily via pip. See also the MacOS X instructions, but you can skip all steps for installing the gcc compilers, since the windows version of pyGDM comes as pre-compiled binary package.

Via pip
-------------------------------

We provide a 64bit windows binary on the pypi repository (tested on Win7 and Win10). Install via

.. code-block:: bash
    
    $ pip install pygdm2

    
Compile using the Anaconda distribution (tested with anaconda3)
------------------------------------------------------------------------------------------
    
1. get the repo (e.g. download from gitlab)

2. install gcc compiler:

   .. code-block:: bash
    
        $ conda install m2w64-toolchain libpython

3. compile fortran parts:

   .. code-block:: bash
    
        $ python setupy.py build

4. install:

   .. code-block:: bash
    
        $ python setupy.py install





Installation under Mac OS X
=============================================

Using the Anaconda distribution
-------------------------------------------------------------

The default compiler on OSX uses a clang which does not support OpenMP. Hence compilation might fail. We therefore suggest using `Anaconda (Mac) <https://www.anaconda.com/download/#macos>`_ and install gcc from the conda repository in a virtualenv (Here the example of python2. python3 was not tested on OSX so far):

   .. code-block:: bash
    
        $ conda create -n python2 python=2.7 anaconda

"anaconda" at the end will copy the whole anaconda distribution to the virtial env. You can omit this option and create a "blank" virtual environment to install only selected packages. 

Next activate the virtualenv and install the required software:

   .. code-block:: bash

        $ source activate python2
        $ xcode-select --install
        $ conda install pip
        $ conda install gcc
        $ pip install pygdm2
        
Also make sure you have the latest versions of numpy and scipy:

   .. code-block:: bash
    
        $ pip install numpy scipy --upgrade
        




Without Anaconda
-------------------------------------------------------------

Alternatively, you can download the latest version and compile it manually without OpenMP support, which should work with the default OSX compiler:

.. code-block:: bash
    
    $ python setup.py install --no-openmp






Authors
=========================

Python code
------------------------
   - P\. R. Wiecha
   - contributions by C\. Majorel


Fortran code
-------------------------
   - C\. Girard
   - A\. Arbouet
   - R\. Marty
   - P\. R. Wiecha



   


