# ----------------------------------------------------------------------
# Config file for Intel ifort  with OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpif90
CC   = mpicc

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS = -I$(HOME)/../kenway/packages/cgnslib_2.5
CGNS_LINKER_FLAGS  = -Wl,-rpath,$(HOME)/../kenway/packages/cgnslib_2.5/LINUX -L$(HOME)/../kenway/packages/cgnslib_2.5/LINUX -lcgns

# ------- Define Compiler Flags ----------------------------------------
FF90_GEN_FLAGS = 
CC_GEN_FLAGS   =

FF90_OPT_FLAGS   =  -fPIC -r8 -O2
CC_OPT_FLAGS     = -O2 -fPIC

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS = -nofor_main

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}


