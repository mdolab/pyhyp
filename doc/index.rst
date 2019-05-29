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

`pyHyp` follows the standard MDO Lab build procedure.
To start, find a configuration file close to your current setup in::

    $ config/defaults

and copy it to ''config/config.mk''. For example::

    $ cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

If you are a beginner user installing the packages on a linux desktop, 
you should use the ``config.LINUX_GFORTRAN.mk`` versions of the configuration 
files. The ``config.LINUX_INTEL.mk`` versions are usually used on clusters.

Once you have copied the config file, compile :ref:`pyHyp` by running::

    $ make

If everything was successful, the following lines will be printed to
the screen (near the end)::

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

Usage with Plot3d Files
=======================

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
			       grid will fold back on itself. In concave corners, additional smoothing will prevent lines
			       from crossing (avooiding negative cells).

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

.. pyHyp boundary conditions example.
   Written by: Ney Seco (February 2016)
   Editted by: 

.. _pyhyp_cgns:

Usage with CGNS Files
=====================

If the initial surface is given in a CGNS file, we can specify boundary conditions
at each open edge of the geometry. The boundary conditions currently supported are:

* Constant X, Y, or Z planes;
* Symmetry X, Y, or Z planes;
* Splay (free edge).

This section will show how we can use ICEM to specify boundary conditions in a CGNS file.

.. NOTE::
   It is still not possible to specify boundary conditions when plot3d files are
   used as inputs. In this case, the surface should be entirely closed or it should
   end at a symmetry plane with the 'mirror' option enabled.

Flat square example
-----------------------------------

.. NOTE::
  If you have a surface geometry, this step is not required, and you
  can proceed to the next section, 'Create Parts with Edges.'

This subsection will show how to create a flat square surface mesh in ICEM. This geometry
is the one used as an example of boundary conditions setup.

1. **Prepare your workspace and open ICEMCFD**

   Create an empty folder anywhere in your computer. Navigate to this folder using the
   terminal and type the following command::

     $ icemcfd

   This will open ICEM, and all the files will be stored in this folder.

2. **Create the surface geometry**

   Under the *Geometry* tab, select *Create/Modify Surface*, as shown below:

      .. image:: images/Figure_CreateSurface.png

   A new menu will show up on the lower-left corner of the screen. Select *Standart Shapes*, then
   *Box*.
   
   Finally, type '1 1 0' in the *X Y Z size* field and click on 'Apply', just as shown below:

      .. image:: images/Figure_CreateBox.png

   This function would usually generate the 6 surfaces of a box, but since we used Z size = 0, it will
   automatically give just a single surface in this case as the upper and lower sides of the box coincide.

3. **Create parts for the edges**

   ICEM groups geometry components into *Parts*. We need to create Parts for the edges so that we can
   easily set up the boundary conditions later on.
   Look at the model tree on the left side of the screen and click with the right-mouse-button on the *Parts*
   brach. Select the *Create Parts* option. Also make sure that the *Surfaces* box in the *Geometry* brach in unchecked
   (otherwise it will be hard to select just the edges):

      .. image:: images/Figure_CreatePart.png

   On the lower-left corner menu change the name of the part to 'EDGE1', then click on the *Create Part by Selection*,
   as shown below:

      .. image:: images/Figure_PartDef.png

   Click on one edge with the left-mouse-button in order to highlight it (see figure below) then click with the
   middle-mouse (scroll) button anywhere to confirm your selection. Now you can check if you have the *EDGE1* component
   under the *Parts* branch in the model tree.

      .. image:: images/Figure_PartSelect.png

   We can repeat the process to create the 'EDGE2' part. Remember to change the part name before selecting another edge.
   This time I chose the upper edge:

      .. image:: images/Figure_Edge2.png

   You can group multiple edges under the same Part if you want to apply the same boundary condition to all of them. For
   instance, we will create 'EDGE3' by grouping the remaining two edges. When selecting the edges, click on both of them
   with the left-mouse-button and only then you click with the middle-mouse-button to create the part.

      .. image:: images/Figure_Edge3.png

   In the end, the *Parts* branch of the model tree should have three components: EDGE1, EDGE2, and EDGE3. Now you can also
   turn the *Surface* box back on. The model tree should look like this:

      .. image:: images/Figure_ModelTree.png

4. **Create blocking**

   Now we need to create the blocks used by ICEM to generate meshes. Under the *Blocking* tab, choose *Create Block*:

      .. image:: images/Figure_CreateBlock.png

   On the lower-left corner menu choose '2D Surface Blocking' as *Type*, and select the 'All Quad' option under *Free Mesh Type*.
   Later, click on the *Select surfaces(s)* button:

      .. image:: images/Figure_BlockOptions.png

   Now click on the flat square in order to highlight it (it should turn white), then click with the middle-mouse-button to
   confirm your selection. Finally, click on the *Apply* button to create the block (the square should turn back to blue). We
   can check if the blocking was done correctly by expanding the *Blocking* branch of the model tree. Click with the
   right-mouse-button on the *Edges* component and turn on the 'Counts' option. This will show how many nodes we have on
   each edge. For now, we should have two nodes per edge, as shown below:

      .. image:: images/Figure_EdgeCount1.png

   Let's increase the number of nodes to make things more interesting. Under the *Blocking* tab, choose *Pre-Mesh Params*:

      .. image:: images/Figure_PreMeshParams.png

   Select the *Scale Sizes* feature on the lower-left menu. Adjust the *Factor* option to 10 and click on *Apply* to globally
   refine the mesh:

      .. image:: images/Figure_PreMeshOptions.png

   Now each edge should show 11 nodes:

      .. image:: images/Figure_EdgeCount2.png

5. **Generate the Pre-mesh**

   It is time to generate the Pre-mesh. Just check the *Pre-mesh* box under the *Blocking* branch of the model tree and choose *Yes*.
   You should see all the surface cells now:

      .. image:: images/Figure_PreMesh.png

   See if everything looks right in your Pre-mesh.

Preparing to export the mesh
-----------------------------------

Just to recap, we have done the following procedures: 

* Created our geometry in ICEM;
* Created Parts grouping edges that will share same boundary conditions;
* Added the surface blocks;
* Generated the Pre-mesh.

.. NOTE::
   The procedures described from now on apply to any geometry. However, make sure you have followed all these steps above
   if you are working with your own geometry.

Now that we confirmed that the Pre-mesh looks right, we can generate the structured mesh. Click with the right-mouse-button on the
*Pre-mesh* component under the *Blocking* branch of the model tree and choose *Convert to MultiBlock Mesh*. Save the project in the
folder you just created. A new branch named *Mesh* should show up in your model tree.

Next, we should select the export format. Under the *Output* folder, click on the *Select Solver* button (a red toolbox):

      .. image:: images/Figure_SelectSolver.png

Choose the following options in the lower-left menu and click *Apply*:

      .. image:: images/Figure_SolverMenu.png

This will apply the CGNS format to our output. We will select the boundary conditions in the next step.

Applying boundary conditions
----------------------------

Under the *Output* folder, click on the *Boundary condition* button:

      .. image:: images/Figure_BCbutton.png

A new window should pop up. As we will apply boundary conditions (BCs) to the edges, expand everything under the *Edges* branch.
You should see the edges Parts we defined previously (EDGE1, EDGE2, and EDGE3) as shown below. In some cases, they may end up under the
*Mixed/Unknown* branch.

      .. image:: images/Figure_BCwindow.png

1. **Symmetry plane**

Let's add a symmetry plane boundary condition to EDGE1. Click on the *Create New* branch under *EDGE1*. Another window will show up,
where you should choose 'BCType' and click on 'Okay'.

      .. image:: images/Figure_BCselection.png

Now go back to the Boundary conditions window. If you click on the green button, you will see several types of Boundary Conditions. The
currently supported types are:

* BCExtrapolate -> Splay;
* BCSymmetryPlane -> Symmetry X, Y, or Z planes;
* BCWall -> Constant X, Y, or Z planes.

For this example, let's choose 'BCSymmetryPlane' for this edge.

      .. image:: images/Figure_BCEdge1.png

Now we need to specify if we want an X, Y, or Z symmetry plane. We will use a 'Velocity' node from the CGNS format in order to specify
the reference plane. For instance, if we add a 'VelocityX' node to this boundary condition, we have a X symmetry plane, and the
same for the other coordinates.

Do the following steps to add a X symmetry plane to the EDGE1 Part:

* Click again on the *Create New* branch under *EDGE1*;
* This time, select 'BCDataset_t' on the new window, and click on 'Okay';
* Now back to the BC window, click on the green button near 'BCTypeSimple' and select 'BCSymmetryPlane';
* Click on the green button that corresponds to the 'Data-Name Identifier (1)' field and choose 'VelocityX';
* Change the 'Data value (1)' field to 1.0.

In the end, the options should look like the figure below:

      .. image:: images/Figure_BCSymOptions.png

In the case of a Y symmetry plane, you should select 'VelocityY' instead, and similarly for a Z symmetry plane. Do not click on
*Accept* yet, otherwise it will close the window. Now let's see the other edges.

2. **Constant Plane**

We will add a constant Y plane to EDGE2. Follow these steps:

* Click on the *Create New* branch under *EDGE2*;
* Select 'BCType' on the new window, and click on 'Okay';
* Click on the green button of the BC window and select 'BCWall'.

We still need to specify which coordinate (X, Y, or Z) should remain constant in this boundary condition. We will do this by adding
a 'Velocity' option to this boundary condition.

* Click again on the *Create New* branch under *EDGE2*;
* This time, select 'BCDataset_t' on the new window, and click on 'Okay';
* Now back to the BC window, click on the green button near 'BCTypeSimple' and select 'BCWall';
* Click on the green button that corresponds to the 'Data-Name Identifier (1)' field and choose 'VelocityY';
* Change the 'Data value (1)' field to 1.0.

In the end, the options should look like the figure below:

      .. image:: images/Figure_BCWallOptions.png

In the case of a constant X plane, you should select 'VelocityX' instead, and similarly for a constant Z plane. Do not click on
*Accept* yet, because we still have one more boundary condition to go!

3. **Splay**

We'll finally add a Splay boundary condition for the two edges included in the EDGE3 Part.

* Click on the *Create New* branch under *EDGE3*;
* Select 'BCType' on the new window, and click on 'Okay';
* Click on the green button of the BC window and select 'BCExtrapolate'.

This one was easier! The same BC will be applied to all edges in this Part. In the end, your boundary conditions tree should
look like this:

      .. image:: images/Figure_BCEdge3.png

Now we can click on *Accept* as we finished adding all the boundary conditions.

Exporting the mesh
------------------

We are ready to export the mesh! Click on the *Write input* button under the *Output* tab:

      .. image:: images/Figure_WriteInput.png

Save the project if asked for. Next we need to select the Multiblock mesh file. The name shown by default should be correct, so just
click on 'Open'. In the next window, click on 'All'.

Another window with CGNS export options will show up. The default options should work fine, but you can compare it with the ones below.
Make sure you are exporting a structured mesh:

      .. image:: images/Figure_CGNSOptions.png

Click on 'Done' to finally conclude export procedure! A CGNS file should appear on your working folder. It will have the same name
of the project file. This CGNS is ready to be used by pyHyp.

Checking the CGNS Structure
---------------------------

You can check if the CGNS file is correct by looking at its structure. If you have cgnslib installed, you
can open the *cgnsview* GUI with the following command::

     $ cgnsview

Open the newly generated CGNS file and expand its tree. For the flat square case, we have the follwing structure:

      .. image:: images/Figure_CGNSView.png

Note the VelocityX node indicating a X symmetry plane boundary condition and a VelocityY node indicating a
constant Y plane boundary condition.

Boundary condition priorities
-----------------------------

The corner nodes share two edges. Each edge may have different boundary conditions. The boundary condition at
the corner node is chosen according to the following priority:

1. Constant X, Y, or Z planes
2. Symmetry X, Y, or Z planes
3. Splay

Therefore, if one edge that has a Splay BC is connected to another edge with a symmetry plane BC, the shared corner node
will be computed with the symmetry plane BC.

Running pyHyp with the generated mesh
-------------------------------------

Create another empty folder and copy the CGNS file exported by ICEM to it. We can add the following Python script to
the same folder (This script is also available in `python/examples/plate/generate_grid.py`, and just the file name was
adjusted for this example)::

  from pyhyp import pyHyp

  fileName = 'plate.cgns'
  fileType = 'CGNS'

  options= {
      # ---------------------------
      #        Input File
      # ---------------------------
      'inputFile': fileName,
      'fileType': fileType,

      # ---------------------------
      #        Grid Parameters
      # ---------------------------
      'N': 65, 
      's0': 1e-6,
      'rMin': 2.5,

      # ---------------------------
      #   Pseudo Grid Parameters
      # ---------------------------
      'ps0': 1e-6,
      'pGridRatio': 1.15,
      'cMax': 5,

      # ---------------------------
      #   Smoothing parameters
      # ---------------------------
      'epsE': 1.0,
      'epsI': 2.0,
      'theta': 0.0,
      'volCoef': 0.3,
      'volBlend': 0.001,
      'volSmoothIter': 10,

      # ---------------------------
      #   BC parameters
      # ---------------------------
      'sigmaSplay': 0.4,
      'nuSplay': 0.95,

      # ---------------------------
      #   Solution Parameters
      # ---------------------------
      'kspRelTol': 1e-15,
      'kspMaxIts': 1500,
      'preConLag': 10,
      'kspSubspaceSize':50,
      'writeMetrics': False,
      }


  hyp = pyHyp(options=options)
  hyp.run()
  hyp.writeCGNS('plate3D.cgns')

Save this script with the name 'generate_grid.py'. Then, navigate to the folder using the terminal and write the following command::

  $ python generate_grid.py

This script will read the plate.cgns file (which contains the surface mesh and the boundary conditions) and will generate the
plate3D.cgns file with the volume mesh. It is important to check the 'MinQuality' column of the screen output. A valid mesh should
have only positive values.

You can also run pyHyp in parallel with the following command::

  $ mpirun -np 4 python generate_grid.py

The option '-np 4' indicates that 4 processors will be used. The results may vary slight due to the parallel solution of the linear system.

Visualizing the mesh in TecPlot 360
-------------------------------------

If you have TecPlot 360 installed in your computer you can visualize the volume mesh. Open a terminal and navigate to the folder
than constains the newly generated CGNS file with the volume mesh. Then type the following command::

  $ tec360

This will open TecPlot 360. On the upper menu select 'File' > 'Load Data Files', then choose your CGNS file. Next, check the 'Mesh' box on the left panel, and click 'Yes'. You will be able to visualize the mesh as shown below:

      .. image:: images/Figure_MeshIso.png

We can see that the boundary conditions where correctly applied. The image below shows a bottom view of the mesh:

      .. image:: images/Figure_MeshBottom.png

Try playing with the different parameters to see their impact on the final mesh. In this case, it is helpful to save a TecPlot
layout file. For instance, place the mesh in a position you want and click on 'File' > 'Save Layout File' and save it with the
name you want (let's say layout_hyp.lay). Then you can open you mesh directly from the command line by typing::

  $ tec360 layout_hyp.lay

Then you don't have to go over the menus all over again!

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

 
