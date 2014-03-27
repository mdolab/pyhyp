.. pyHyp documentation master file, created by
   sphinx-quickstart on Sun Oct 13 13:46:01 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _pyhyp:

=====
pyHyp
=====

Contents:

.. toctree::
   :maxdepth: 2

.. _pyhyp_introduction:

Introduction
============

`pyHyp` is hyperbolic mesh generator that is capable of automatically
generating two or three dimensional 3D meshes around simple geometric
configurations. The basic idea is to start with an initial surface (or
line) corresponding to the geometry of interest and then *grow* or
*extrude* the mesh in successive layers until it reaches a sufficient
distance from the original surface. In the process, the entire space
surrounding the geometry is meshed. In the usage_ section.

.. _pyhyp_installation:

Installation 
=============

Prerequisites
------------- 

`pyHyp` depends heavily on other packages to do much of the underlying
"heavy lifting". The following external components are required for
pyHyp:


1. `CGNS Library <http://cgns.sourceforge.net>`_
 
The CGNS library is used to provide output CGNS functionality for
pyHyp. The latest CGNS version 3.1 is recommended, but the older
version 2.5.x may also be used. The After downloading, untar::

   $ tar -xzf cgnslib_3.1.4-2.tar.gz 

.. NOTE:: cmake and ccmake are required for building the CGNS
   library.

Enter cgnslib_3.1.4 and type::

   $ cmake . 

By default, the CGNS library does not include the Fortran bindings
that are required for `pyHyp`. This needs to be enabled using the
cmake configure utility, `ccmake`.::

   $ ccamke .

A "GUI" appears and toggle ENABLE_FORTRAN by pressing [enter]. Type
'c' to reconfigure and 'g' to generate and exit. Then build the
library using::

   $ make

If you are compiling on machine without root access, you're done.
If you are on a local machine with root access, it is usually
easiest to 'install' into the default directory using::

    $ sudo make install

.. NOTE:: It is also suggested that the BUILD_CGNSTOOLS option be
   used in the configuration, but this is optional. 

2. `PETSc <http://www.mcs.anl.gov/petsc/index.html>`_ 

Download the latest tarball and extract::

   $ tar -xzf petsc-3.4.1.tar.gz

PETSc must be first configured. There are a wide variety of
options. The only ones that is strictly necessary are
--with-shared-libraries and --with-fortran-interfaces. However, it is
recommend the following options are used since some are required for
pyWarp. To compile without debugging use: --with-debugging=no. It is
HIGHLY recommended to use debugging until ready to perform production
runs. In that case, you should compile with a separate architecture
for the optimized build. In this case use '--PETSC_ARCH=real-opt' for
example::

       $ ./configure --with-shared-libraries --download-superlu_dist=yes --download-parmetis=yes --download-metis=yes --with-fortran-interfaces=1 --with-debuggig=yes --with-scalar-type=real --PETSC_ARCH=real-debug 

After the configuration step, PETSc must be built. This is
accomplished with the command provided at the end of the configure
script. It will look something like below::

   $ make PETSC_DIR=/home/gaetan/Downloads/petsc-3.4.1 PETSC_ARCH=real-debug all

The last step is to add PETSC_DIR and PETSC_ARCH entries to your
.bashrc file. This is essential. It should look something like
this: (Make sure the CORRECT directory is used!)::

    export PETSC_ARCH=real-debug
    export PETSC_DIR=/home/user/packages/petsc-3.4.1
   
Make sure the .bashrc file is sourced before pyHyp is compiled using::

   $ source ~/.bashrc

Or simply open a new terminal before compiling pyHyp.

See the documentation of each of these packages for further
information. 

.. NOTE:: A working MPI is not strictly required. However, in most
   cases PETSc should be configured with MPI. `pyHyp` itself, however,
   does not currently run in parallel.

Compilation 
------------ 

`pyHyp` is built using a standard `make` utility on unix-like
systems. To begin compilation, type `make`::

    $ make

A list of previously tested architectures should show::

    LINUX_GFORTRAN_OPENMPI
    LINUX_INTEL_OPENMPI

Before starting compilation, copy a default configuration file from
the `config/defaults` directory into `config`. From the root `pyHyp`
directory run (or example)::

  $ cp config/defaults/config.LINUX_INTEL_OPENMPI config/

Next, edit the configuration file. This will allow you to set which
 compiler to use the compilation flags as well as the location for
 the CGNS include files and library. 

Other platforms/compiler combination are most likely possible but have
 not been tested. On a new platform, simply modify one of the existing
 configuration files to suit your needs. 

After the configuration file has been modified, compilation can be
started with (for example)::

    $ make LINUX_GFORTRAN_OPENMPI

Compilation is successful when the following two lines are displayed
at the end of the compilation::

    Testing if module hyp can be imported...
    Module hyp was successfully imported.

Finally, add the directory containing the pyhyp folder to the
$PYTHONPATH variable in your bashrc file::

    export $PYTHONPATH:/path/to/folder/containing/pyhyp/

.. _pyhyp_theory:

Theory
======
Most the theory for `pyHyp` was taken from Chan.

.. _pyhyp_usage:

Usage
=====

A complete sample script to generate a grid is given below. This
particular example is available under `python/examples/sphere/runSphere.py`::

  from pyhyp import pyHyp
  options= {
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 73, 
    's0':1e-4,
    'rMin':10,
 
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':1e-4,
    'pGridRatio':1.1,
    'cMax': 5,
    
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'volCoef': 0.16,
    'volBlend': 0.0005,
    'volSmoothIter': 15,

    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 1500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

  hyp = pyHyp.pyHyp('3d',fileName='uneven_sphere.fmt', options=options)
  hyp.run()
  hyp.writeCGNS('sphere.cgns')

Each section of the example is now described::

  from pyhyp import pyHyp

is the preferred way of importing the `pyHyp` module into python.

The options dictionary is used to provide all run-time options to
pyHyp to control the generation of the grid. Each option, its
description and the general effect on the grid generation process is
now described:

===================  ========  =======================================================================================
Parameter              Type      Description
===================  ========  =======================================================================================
`N`                   `int`    Number of grid levels to march. This determines the grid dimension 
                               in the off-wall direction. Typically this should be a "multi-grid" friendly number.

`s0`                 `float`   Initial off-wall (normal) spacing of grid. This is taken to
                               be constant across the entire geometry.  The units are consistent 
			       with the rest of the geometry.

`rMin`                `float`  Relative distance in the normal direction to march. It is 
                               specified as a multiple of the radius of the sphere 
			       enclosing the initial surface geometry. If symmetry is specified, 
			       the full mirrored geometry is used to compute the sphere's radius. 
			       Most wing geometries will have `rMin` between 10 and 20, that is the 
			       farfield boundary is 20 spans away from the geometry. 

`cMax`               `float`   The maximum permissible ratio of marching direction length to 
                               the any other in-plane edge.  This parameter effectively operates 
			       as a CFL-type limit. If a step would require a step which would result 
			       in a ratio `c` greater than `cMax`, the step is automatically split internally 
			       to respect this user-supplied limit. Typical values of `cMax` are around 6-8.
			       Increased robustness can be achieved at the expense of computational cost by lowering `cMax`.   

`nonLinear`          `float`   Use the nonlinear formulation. This is experimental and not 
                               currently recommended and may not work at all.

`slExp`              `float`   Exponent for the :math:`S_l` compuitation.
                               The :math:`S_l` value serves the same purpose as found in Chan et al. 
			       but the computation is different. The :math:`S_l` computation in Chan is 
			       given as :math:`\sqrt{\frac{N-1}{l-1}}` for :math:`l > 2`. 
			       
`ps0`                `float`   Initial pseudo offwall spacing. This spacing *must* be less than or equal to `s0`. This
                               is actual spacing the hyperbolic scheme uses. The solver may take many pseudo steps
			       before the first real grid level at `s0`. 

`pGridRatio`         `float`   The ratio between sucessive levels in the pseudo grid. This will be typically somewhere
                               between ~1.05 for large grids to 1.2 for small grids. This number is *not* the actual grid
			       spacing of the final grid; that spacing ratio is computed and displayed at the beginning
			       of a calculation. The pGridRatio *must* be smaller than that number. 

`epsE`                `float`  The explict smoothing parameter. See the :ref:`Theory<pyhyp_theory>` section for more information.
                               Typical values are approximately 1.0. Increasing the explict smoothing may result in  a
			       smoother grid, at the expense of orhtogonality. If the geometry is very sharp corners,
                               too much explict smoothing will cause the solver to rapidly "soften" the corner and the 
			       grid will fold back on itself. 

`epsI`                `float`  Implict smoothing parameter. See the :ref:`Theory<pyhyp_theory>` section for more information.
                               Typical values are from 2.0 to 6.0. Generally increasing the implict coefficient results
			       in a more stable solution procedure.

`theta`               `float`  Kinsley-Barth coefficient See the :ref:`Theory<pyhyp_theory>` section for more information.
                               Only a single theta value is used for both directions. Typical values are ~2.0 to ~4.0.

`volCoef`             `float`  Coefficient used in point-Jacobi local volume smoothing algorithm. Typically this
                               value is 0.16 and need not be modified. Use more `volSmoothIter` for stronger local smoothing.
			    
`volBlend`            `float`  The gloabl volume blending coefficinet. See the :ref:`Theory<pyhyp_theory>` section for more information.
                               This value will typically be very small, especially if you widely varrying cell sizes. 
			       Typically values are from ~0 to 0.001. Default is 0.0001

`volSmoothIter`       `int`    The number of point-Jacobi local volume smoothing iterations to perform. Typical values
                               are ~5 to ~25. Default is 10.

`kspRelTol`           `float`  Tolerance for the solution of the linear system at each iteration. Typically :math:`1\times 10^{-8}` 
                               is sufficient. Ver dificult cases may benefit from a tigher convergence tolerance.

`kspMaxIts`           `int`    Maximum number of iterations to perform for each step. Default is 500 which should be sufficient
                               for most cases. 

`preConLag`           `int`    Lag the update of the preconditioner by this number of iterations. The default value of 10 
                               will typically not need to be changed. 

`kspSubspaceSize`     `int`    Size of the ksp subspace. Default is 50. Very large and difficult problems may befefit
                               from a larger subspace size. 

`writeMetrics`       `bool`    Flag to write the mesh gradients to the solution file. This option should only be used
                               for debugging purposes. 
===================  ========  =======================================================================================

The next line of code::

  hyp = pyHyp.pyHyp('3d',fileName='uneven_sphere.fmt', options=options)

generates the pyHyp object. Here we have a 3D problem and we are using
the surface mesh `uneven_sphere.fmt`. This is a plot3d file which has
been generated in an external meshing program, typically ICEMCFD. 

.. NOTE:: When exporting a surface mesh from ICEMCFD in plot3d format
          use the following options:
	  
	  * Formatted
	  * Whole
	  * Double
	  * No IBLANK Array
	    
.. WARNING:: It is essential that the normals of each of the surface
   patches point in the *OUTWARD* direction, i.e. the marching
   direction. `pyHyp` does not check for normal consistency and the
   mesh generation process will fail with inconsistent normals. 

The next two lines performs the actual generation and writes the
resulting grid to a cgns file::

  hyp.run()
  hyp.writeCGNS('sphere.cgns')


The output of the run should look similar to the following::

  Loading ascii plot3D file: uneven_sphere_large.fmt ...
   -> nSurf = 6
   -> Surface Points: 6534
  Problem sucessfully setup!
  #--------------------#
  Grid Ratio:  1.1429
  #--------------------#
  #----------------------------------------------------------------------------------------------------------------------------
  # Grid  |     CPU    | Sub  | KSP  |     Sl     | Grid       | Grid       |     Min    |   deltaS   |    cMax   |   min R   | 
  # Level |     Time   | Iter | Its  |            | Sensor Max | Sensor Min |  Quality   |            |           |           | 
  #----------------------------------------------------------------------------------------------------------------------------
       2  0.38516E-01      1     13  0.17656E+00  0.10000E+01  0.10000E+01  0.89103E+00  0.11000E-03  0.23373E-01  0.95337E-04
       3  0.76673E-01      2     13  0.21128E+00  0.10008E+01  0.99723E+00  0.89176E+00  0.13310E-03  0.28286E-01  0.31556E-03
       4  0.98730E-01      1     13  0.22227E+00  0.10009E+01  0.99697E+00  0.89177E+00  0.14641E-03  0.31126E-01  0.44246E-03
       5  0.12045E+00      1     13  0.23160E+00  0.10010E+01  0.99668E+00  0.89185E+00  0.16105E-03  0.34252E-01  0.58204E-03
       6  0.14223E+00      1     13  0.23988E+00  0.10012E+01  0.99634E+00  0.89191E+00  0.17716E-03  0.37695E-01  0.73558E-03
       7  0.16397E+00      1     13  0.24743E+00  0.10013E+01  0.99596E+00  0.89196E+00  0.19487E-03  0.41487E-01  0.90448E-03
       8  0.18566E+00      1     13  0.25447E+00  0.10014E+01  0.99554E+00  0.89200E+00  0.21436E-03  0.45663E-01  0.10903E-02
       9  0.20734E+00      1     13  0.26111E+00  0.10015E+01  0.99505E+00  0.89203E+00  0.23579E-03  0.50263E-01  0.12946E-02
      10  0.24493E+00      2     13  0.27357E+00  0.10016E+01  0.99388E+00  0.89203E+00  0.28531E-03  0.60914E-01  0.17667E-02
      < iterations skipped for brevity> 
      70  0.32058E+01      3     34  0.94224E+00  0.93695E+00  0.90972E+00  0.38760E+00  0.33863E+00  0.50000E+01  0.67257E+01
      71  0.33147E+01      3     35  0.96129E+00  0.93688E+00  0.91092E+00  0.38618E+00  0.38543E+00  0.50000E+01  0.76861E+01
      72  0.34279E+01      3     35  0.98067E+00  0.93687E+00  0.91211E+00  0.38490E+00  0.43909E+00  0.50000E+01  0.87798E+01
      73  0.35392E+01      3     35  0.10004E+01  0.93690E+00  0.91329E+00  0.38377E+00  0.50063E+00  0.50000E+01  0.10027E+02


Several important parameters are displayed to inform the user of the
solution progress. The most of important of which is the `Min Quality`
column. This column displays the minimum quality of all the cells in
the most recently computed layer of cells. For a valid mesh, these
must be all greater than zero. The appearance of negative values 

Code Reference
==============

..
   .. automodule:: hypmodule
      :members:





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex` 
* :ref:`search`

 
