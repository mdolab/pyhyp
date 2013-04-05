#      ******************************************************************
#      *                                                                *
#      * File:          Hyp_Common.mk                                   *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 04-05-2013                                      *
#      * Last modified: 04-05-2013                                      *
#      *                                                                *
#      ******************************************************************

HYP_MODDIR = $(HYP_DIR)/mod
HYP_OBJDIR = $(HYP_DIR)/obj
HYP_LIBDIR = $(HYP_DIR)/lib

#      ******************************************************************
#      *                                                                *
#      * Include the file describing the compiler settings.             *
#      *                                                                *
#      ******************************************************************

HYP_COMPILERS = $(HYP_DIR)/config.mk
ifneq ($(MAKECMDGOALS),clean)
include ${HYP_COMPILERS}
endif

#      ******************************************************************
#      *                                                                *
#      * Redefine .SUFFIXES to be sure all the desired ones are         *
#      * included.                                                      *
#      *                                                                *
#      ******************************************************************

.SUFFIXES: .o .f .F .f90 .F90 

#      ******************************************************************
#      *                                                                *
#      * Arguments of make clean.                                       *
#      *                                                                *
#      ******************************************************************

MAKE_CLEAN_ARGUMENTS = *~ *.o *.mod *.il *.stb c_*

#      ******************************************************************
#      *                                                                *
#      * Compiler flags to compile the sources.                         *
#      * The macro's ADDITIONAL_FF90_FLAGS and ADDITIONAL_CC_FLAGS      *
#      * make it possible that every subdirectory adds its own specific *
#      * compiler flags, if necessary.                                  *
#      *                                                                *
#      ******************************************************************

FF90_ALL_FLAGS   = -I$(HYP_MODDIR) $(CGNS_INCLUDE_FLAGS) \
		   $(FF90_GEN_FLAGS) $(FF90_OPT_FLAGS) $(PETSC_INCLUDE_FLAGS) 

CC_ALL_FLAGS = -I$(HYP_MODDIR) $(CGNS_INCLUDE_FLAGS) \
		   $(CC_GEN_FLAGS) $(CC_OPT_FLAGS) $(PETSC_INCLUDE_FLAGS)
