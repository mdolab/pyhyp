#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 03-05-2013                                      *
#      * Last modified: 03-05-2013                                      *
#      *                                                                *
#      ******************************************************************

HYP_DIR = ../../..
HYP_COMMON_FILE = $(HYP_DIR)/Hyp_Common.mk
HYP_RULES_FILE  = $(HYP_DIR)/rulesSources.mk
include ${HYP_COMMON_FILE}
include ${HYP_RULES_FILE}

LINKER_ALL_FLAGS = -L$(HYP_LIBDIR) -lhyp $(CGNS_LINKER_FLAGS) \
		    $(PETSC_LINKER_FLAGS) $(LINKER_FLAGS) $(CGNS_LINKER_FLAGS)

PYTHON_OBJECTS = fortranobject.o \
		 hypmodule.o \
		 hyp-f2pywrappers2.o \
# This file is sometimes for F77 wrapper
#		 hyp-f2pywrappers.o \

default: all

all:

# Generate Python inlude directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Include Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print(numpy.get_include())'))
	 @echo "#------------------------------------------------------#"
	 @echo Numpy Include Directory: $(NUMPY_INCLUDES)
	 @echo "#------------------------------------------------------#"
# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) get_f2py.py))
	 @echo "#------------------------------------------------------#"
	 @echo f2py root directory: $(F2PY_ROOT)
	 @echo "#------------------------------------------------------#"

#       need to have f2py_f2cmap in this folder....
	cp f2py_f2cmap.ref .f2py_f2cmap
	$(F2PY) hyp.pyf
	@echo " "
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c hypmodule.c
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) -c \
	$(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# # This file is for F77 wrappers if they exist
#	$(FF90) $(FF90_ALL_FLAGS) -I$(HYP_MODDIR) -c hyp-f2pywrappers.f

	$(FF90) $(FF90_ALL_FLAGS) -I$(HYP_MODDIR) -c hyp-f2pywrappers2.f90
	$(FF90) -shared $(PYTHON_OBJECTS) $(LINKER_ALL_FLAGS) -o hyp.so
	$(PYTHON) importTest.py

	mv hyp.so ../../../pyhyp
# Clean up temporary files
	-rm *.o
	-rm *.f90 *.c
	-rm .f2py_f2cmap

clean:
	@echo "        Making clean in f2py... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src:   $(FF90_OBJECTS) $(CC_OBJECTS)
