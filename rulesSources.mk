#      ******************************************************************
#      *                                                                *
#      * File:          rulesSources.mk                                 *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 04-05-2013                                      *
#      * Last modified: 04-05-2013                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Rules to make the objects. These are the general  *
#      *              rules. If in a subdirectory different rules must  *
#      *              be used, this file should not be included.        *
#      *                                                                *
#      ******************************************************************

.F90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

.f.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

.F.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

.c.o:   Makefile
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
