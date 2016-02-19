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
Most the theory for `pyHyp` was taken from Chan and Steger. Check the references folder

.. _pyhyp_usage:

Usage
=====

A complete sample script to generate a grid is given below. This
particular example is available under `python/examples/sphere/runSphere.py`::

  from pyhyp import pyHyp
  fileName = 'even_sphere.fmt'
  #fileName = 'uneven_sphere.fmt'
  #fileName = 'uneven_sphere_large.fmt'

  options= {

      # ---------------------------
      #        Input File
      # ---------------------------
      'inputFile':fileName,
      'fileType':'Plot3d',

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
      'epsI': 2.0,
      'theta': 3.0,
      'volCoef': .16,
      'volBlend': 0.0001,
      'volSmoothIter': 20,

      # ---------------------------
      #   Solution Parameters
      # ---------------------------
      'kspRelTol': 1e-10,
      'kspMaxIts': 1500,
      'preConLag': 5,
      'kspSubspaceSize':50,
      }

  hyp = pyHyp(options=options)
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
`inputFile`           `char`   Name of the file that contains the surface mesh. Here we have a 3D problem and we are
                               using the surface mesh `uneven_sphere.fmt`. This is a plot3d file which has
                               been generated in an external meshing program, typically ICEMCFD.

`fileType`            `char`   Type of the input file. Use either 'Plot3d' or 'CGNS'.

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
			       in a more stable solution procedure. Usually this value should be twice the explicit smoothing parameter.

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

  hyp = pyHyp(options=options)

generates the pyHyp object.

.. NOTE:: When exporting a surface mesh from ICEMCFD in plot3d format
          use the following options:
	  
	  * Formatted
	  * Whole
	  * Double
	  * No IBLANK Array
	    
.. WARNING:: It is essential that the normals of each of the surface
   patches point in the *OUTWARD* direction, i.e. the marching
   direction.

The next two lines perform the actual generation and write the
resulting grid to a cgns file::

  hyp.run()
  hyp.writeCGNS('sphere.cgns')


The output of the run should look similar to the following::

 #--------------------#
  Total Nodes:     150 
  Unique Nodes:     98 
  Total Faces:      96 
 #--------------------#
  Normal orientation check ...
  Normals seem correct!
 #--------------------#
 Grid Ratio:  1.1420 
 #--------------------#
 #-------------------------------------------------------------------------------------------------------------------------------------------
 # Grid  |     CPU    | Sub  | KSP  |     Sl     | Grid       | Grid       |     Min    |   deltaS   |    cMax    |    min R   |    max     | 
 # Level |     Time   | Iter | Its  |            | Sensor Max | Sensor Min |  Quality   |            |            |            |  KStretch  | 
 #-------------------------------------------------------------------------------------------------------------------------------------------

       2  0.26410E-02      1     10  0.17783E+00  0.10000E+01  0.10000E+01  0.65528E+00  0.11000E-03  0.32833E-03  0.10000E-03  0.00000E+00 
       3  0.60689E-02      2     10  0.21280E+00  0.99943E+00  0.99923E+00  0.64792E+00  0.13310E-03  0.39724E-03  0.33100E-03  0.10552E+01 
       4  0.81701E-02      1     10  0.22387E+00  0.99941E+00  0.99919E+00  0.64731E+00  0.14641E-03  0.43693E-03  0.46410E-03  0.11350E+01 
       5  0.10073E-01      1     10  0.23327E+00  0.99938E+00  0.99915E+00  0.64681E+00  0.16105E-03  0.48059E-03  0.61051E-03  0.11364E+01 
       6  0.11975E-01      1     10  0.24160E+00  0.99934E+00  0.99910E+00  0.64637E+00  0.17716E-03  0.52861E-03  0.77156E-03  0.11371E+01 
       7  0.13715E-01      1     10  0.24921E+00  0.99930E+00  0.99904E+00  0.64596E+00  0.19487E-03  0.58143E-03  0.94872E-03  0.11376E+01 
       8  0.15457E-01      1     10  0.25630E+00  0.99926E+00  0.99898E+00  0.64558E+00  0.21436E-03  0.63951E-03  0.11436E-02  0.11379E+01 
       9  0.17204E-01      1     10  0.26299E+00  0.99921E+00  0.99891E+00  0.64521E+00  0.23579E-03  0.70340E-03  0.13579E-02  0.11381E+01 
      10  0.20547E-01      2     10  0.27554E+00  0.99909E+00  0.99874E+00  0.64482E+00  0.28531E-03  0.85092E-03  0.18531E-02  0.11379E+01
      < iterations skipped for brevity> 
      70  0.17639E+00      1     13  0.94933E+00  0.91466E+00  0.90983E+00  0.44618E+00  0.70716E+00  0.85351E+00  0.70706E+01  0.10857E+01 
      71  0.17845E+00      1     13  0.96300E+00  0.91433E+00  0.90996E+00  0.44321E+00  0.77788E+00  0.89981E+00  0.77778E+01  0.10933E+01 
      72  0.18323E+00      2     15  0.99094E+00  0.91399E+00  0.91040E+00  0.44094E+00  0.94123E+00  0.99767E+00  0.94113E+01  0.10859E+01 
      73  0.18649E+00      1     15  0.10052E+01  0.91396E+00  0.91071E+00  0.43894E+00  0.10354E+01  0.10493E+01  0.10353E+02  0.10874E+01


Several important parameters are displayed to inform the user of the
solution progress. The most of important of which is the `Min Quality`
column. This column displays the minimum quality of all the cells in
the most recently computed layer of cells. For a valid mesh, these
must be all greater than zero.

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

 
