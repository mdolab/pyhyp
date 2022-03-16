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

%.o: %.F90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

%.o: %.f90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

%.o: %.f
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

%.o: %.F
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

%.o: %.c
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(HYP_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
