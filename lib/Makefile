#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 04-05-2013                                      *
#      * Last modified: 04-05-2013                                      *
#      *                                                                *
#      ******************************************************************
#
HYP_DIR = ..
HYP_COMMON_FILE = $(HYP_DIR)/Hyp_Common.mk
include ${HYP_COMMON_FILE}

default: all

all: 	libhyp

clean:
	@echo "        Making clean in lib ... "
	rm -f lib/*

libhyp:
	@echo "        creating library ..."
	@echo
	$(AR) $(AR_FLAGS) libhyp.a $(HYP_OBJDIR)/*.o
	@echo
	@echo "        library libhyp.a created."
	@echo
