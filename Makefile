#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 04-05-2013                                      *
#      * Last modified: 04-05-2013                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    = src/modules\
		src/2D\

default:
	@echo "Usage: make <arch>"
	@echo "Supported architectures: LINUX_GFORTRAN_OPENMPI"
	@echo "                         LINUX_INTEL_OPENMPI"

all:	 default

clean:
	@echo " Making clean ... "

	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f lib/lib* mod/* obj/*

hyp:
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)

#      ******************************************************************
#      *                                                                *
#      * Platform specific targets.                                     *
#      *                                                                *
#      ******************************************************************

LINUX_GFORTRAN_OPENMPI:
	mkdir -p obj
	mkdir -p mod
	if [ ! -f "config/config.LINUX_GFORTRAN_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_GFORTRAN_OPENMPI.mk config.mk
	make hyp
	(cd src/python/f2py && make)

LINUX_INTEL_OPENMPI:
	mkdir -p obj
	mkdir -p mod
	if [ ! -f "config/config.LINUX_INTEL_OPENMPI.mk" ]; then cp "config/defaults/config.LINUX_INTEL_OPENMPI.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_OPENMPI.mk config.mk
	make hyp
	(cd src/python/f2py && make)
