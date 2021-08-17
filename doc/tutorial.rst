.. _pyhyp_tutorial:

Tutorial
========

A complete sample script to generate a grid is given below. This
particular example is available under `pyhyp_examples/BWB/runBWB.py`.

.. literalinclude:: ../pyhyp_examples/BWB/runBWB.py

Each section of the example is now described.

.. literalinclude:: ../pyhyp_examples/BWB/runBWB.py
   :start-after: # rst import (start)
   :end-before: # rst import (end)

This is the preferred way of importing the `pyHyp` module into python.

The options dictionary is used to provide all run-time options to
pyHyp to control the generation of the grid.
A description of each option and its general effect on the grid generation process is
explained in :ref:`pyhyp_options`.

The next line of code

.. literalinclude:: ../pyhyp_examples/BWB/runBWB.py
   :start-after: # rst object
   :end-before: # rst run

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
resulting grid to a cgns file:

.. literalinclude:: ../pyhyp_examples/BWB/runBWB.py
   :start-after: # rst run

The output of the run should look similar to the following::

  #--------------------#
   Total Nodes:   14105
   Unique Nodes:  13457
   Total Faces:   13376
  #--------------------#
   Normal orientation check ...
   Normals are consistent!
   Determining topology ...
   Topology complete.
  #--------------------#
  Grid Ratio:  1.2532
  #--------------------#
  #-------------------------------------------------------------------------------------------------------------------#
  # Grid | CPU  | Sub | KSP | nAvg |  Sl  | Sensor | Sensor | Min     | Min     |  deltaS  | March    | cMax  | Ratio |
  # Lvl  | Time | Its | Its |      |      | Max    | Min    | Quality | Volume  |          | Distance |       | kMax  |
  #-------------------------------------------------------------------------------------------------------------------#
        2    0.1     2    11      0  0.055  1.00006  0.98387  0.38387  0.359E-10  0.314E-05  0.451E-05  0.0011  0.0000
        3    0.2     2    11      0  0.064  1.00010  0.96703  0.35947  0.657E-10  0.493E-05  0.116E-04  0.0018  1.9184
        4    0.3     1    11      0  0.067  1.00012  0.96111  0.35650  0.103E-09  0.618E-05  0.165E-04  0.0022  1.5882
        5    0.4     2    11      0  0.074  1.00024  0.90864  0.35575  0.174E-09  0.787E-05  0.305E-04  0.0035  1.7282
        6    0.4     1    11      0  0.076  1.00031  0.89199  0.35610  0.237E-09  0.987E-05  0.383E-04  0.0035  1.3561
        7    0.5     1    11      0  0.079  1.00037  0.87508  0.35613  0.371E-09  0.124E-04  0.482E-04  0.0044  1.5634
        8    0.5     2    11      0  0.084  1.00083  0.78268  0.35729  0.548E-09  0.155E-04  0.761E-04  0.0069  1.4708
        9    0.6     1    11      0  0.087  1.00108  0.75129  0.36201  0.717E-09  0.194E-04  0.916E-04  0.0069  1.3003
       10    0.6     1    11      0  0.089  1.00124  0.74002  0.36599  0.104E-08  0.243E-04  0.111E-03  0.0086  1.4222
       < iterations skipped for brevity >
       78   33.9    22    29      0  0.904  0.98475  0.97354  0.32556  0.347E+04  0.730E+01  0.562E+03  2.5000  1.2603
       79   36.1    22    29      0  0.936  0.98435  0.97427  0.32565  0.688E+04  0.930E+01  0.708E+03  2.5000  1.2570
       80   38.0    21    29      0  0.968  0.98388  0.97469  0.32573  0.138E+05  0.119E+02  0.885E+03  2.5000  1.2551
       81   39.9    20    29      0  1.000  0.98341  0.97476  0.32579  0.276E+05  0.152E+02  0.110E+04  2.5000  1.2540

Several important parameters are displayed to inform the user of the
solution progress. The most of important of which is the `Min Quality`
column. This column displays the minimum quality of all the cells in
the most recently computed layer of cells. For a valid mesh, these
must be all greater than zero.