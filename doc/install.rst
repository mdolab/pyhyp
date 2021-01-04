.. _pyhyp_installation:

Installation 
=============

Prerequisites
------------- 

`pyHyp` depends heavily on other packages to do much of the underlying
"heavy lifting". The following external components are required for
pyHyp:

- CGNS Libarary
- PETSc

See the MDO Lab installation guide `here <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/installInstructions/install3rdPartyPackages.html#installthirdpartypackages>`_ for the supported versions and installation instructions.

.. NOTE:: A working MPI is not strictly required. However, in most
   cases PETSc should be configured with MPI.

Compilation 
------------ 
``pyHyp`` follows the standard MDO Lab build procedure.
To start, first clone the repo. For stability we recommend checking out a tagged release.

Next, find a configuration file close to your current setup in::

    $ config/defaults

and copy it to ''config/config.mk''. For example::

    $ cp config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk config/config.mk

If you are a beginner user installing the packages on a Linux desktop, 
you should use the ``config.LINUX_GFORTRAN_OPENMPI.mk`` versions of the configuration 
files. The ``config.LINUX_INTEL.mk`` versions are usually used on clusters.

Once you have copied the config file, compile :ref:`pyHyp` by running::

    $ make

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module hyp can be imported...
   Module hyp was successfully imported.

Finally, install the Python interface with::

    pip install .

Testing Your Installation
-------------------------

To test your installation, you can run some of the scripts in the `python/examples` folder.
Some of these require you to have `cgnsUtilities <https://github.com/mdolab/cgnsutilities>`_ installed.

After you run some of the files, you will get a message like this::

  *** The MPI_Attr_get() function was called after MPI_FINALIZE was invoked.
  *** This is disallowed by the MPI standard.
  *** Your MPI job will now abort.
  [MDO-John:7977] Local abort after MPI_FINALIZE completed successfully; not able to aggregate error messages, and not able to guarantee that all other processes were killed!
  
Despite its scary look, this is a non-issue and means that the script successfully finished.
