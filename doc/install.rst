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
   cases PETSc should be configured with MPI.

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