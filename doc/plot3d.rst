.. _pyhyp_plot3d:

Usage with Plot3d Files
=======================

A complete sample script to generate a grid is given below. This
particular example is available under `examples/sphere/runSphere.py`::

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
pyHyp to control the generation of the grid.
A description of each option and its general effect on the grid generation process is
explained in :ref:`pyhyp_options`.

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