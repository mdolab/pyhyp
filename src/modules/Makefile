#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 04-05-2013                                      *
#      * Last modified: 04-05-2013                                      *
#      *                                                                *
#      ******************************************************************
HYP_DIR = ../..
HYP_COMMON_FILE = $(HYP_DIR)/Hyp_Common.mk
HYP_RULES_FILE  = $(HYP_DIR)/rulesSources.mk
include ${HYP_COMMON_FILE}
include ${HYP_RULES_FILE}
vpath %.o $(HYP_OBJDIR)

# Note: precision MUST be compiled first.
FF90_OBJECTS_1 = precision.o\

FF90_OBJECTS_2 =	communication.o\
			hypInput.o\
			hypData.o\
			panel.o

default: all

all:
	make     src_1
	make     src_2
	cp -f *.mod $(HYP_MODDIR)
clean:
	@echo "        Making clean in modules... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src_1:	$(FF90_OBJECTS_1) 
src_2:	$(FF90_OBJECTS_2) 
