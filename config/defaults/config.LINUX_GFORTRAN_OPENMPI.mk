# ----------------------------------------------------------------------
# Config file for Gfortran  with OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpifort
CC   = mpicc


# ------- Define CGNS Inlcude and linker flags -------------------------
# Define the CGNS include directory and linking flags for the CGNS library.
# We are assuming that HDF5 came from PETSc so it is included in ${PETSC_LIB}.
# Otherwise you will have to specify the HDF5 library.
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

# ------- Define Compiler Flags ----------------------------------------
FF90_GEN_FLAGS = -fPIC
CC_GEN_FLAGS   = -fPIC

FF90_OPT_FLAGS   =  -fPIC -fdefault-real-8 -O2
CC_OPT_FLAGS     = -O2

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS =

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/lib/petsc/conf/variables # PETSc 3.6
#include ${PETSC_DIR}/conf/variables # PETSc 3.5
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py
