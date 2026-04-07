FC               = mpifort # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -cpp -w -O3 -std=legacy -g -fmax-errors=1#-fallow-argument-mismatch -fcheck=all
CPFLAGS         += -E
ifeq ($(TEST),1)
    PPICLF_FFLAGS += -DTEST -g -fcheck=bounds
endif

ifeq ($(PERF),1)
    PPICLF_FFLAGS += -DPERF
endif

TOOLCHAIN=GNU
PREPROCDEFINES+=-DNDIM=3
 	           
INSTALL_LOCATION = .

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY PATHS
SOURCE_ROOT_PPICLF=$(INSTALL_LOCATION)/ppiclF/source
SOURCE_ROOT_USER=$(INSTALL_LOCATION)/user_files
BUILD_ROOT_PPICLF=$(INSTALL_LOCATION)/ppiclF/build
BUILD_ROOT_USER=$(BUILD_ROOT_PPICLF)/user




# COMPILER FLAGS
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)/include -I$(SOURCE_ROOT_USER) -I$(SOURCE_ROOT_USER)/include

PPICLF_FFLAGS+=$(PPICLF_IFLAGS)

ifeq ($(TOOLCHAIN),GNU)
    PPICLF_FFLAGS += -ffree-line-length-none -Wno-missing-include-dirs -fmax-errors=5 -g -gdwarf-4 -O0 -fno-omit-frame-pointer -cpp -Wall
    PPICLF_USER_FFLAGS:=$(PPICLF_FFLAGS) -J$(BUILD_ROOT_USER) -I$(BUILD_ROOT_PPICLF)
    PPICLF_FFLAGS+=-J$(BUILD_ROOT_PPICLF) -I$(BUILD_ROOT_USER)
    PREPROC_FLAG = -E
else ifeq ($(TOOLCHAIN),NVHPC)
    PPICLF_FFLAGS += -g -O0 -Wall
    PPICLF_USER_FFLAGS:=$(PPICLF_FFLAGS) -module $(BUILD_ROOT_USER) -I$(BUILD_ROOT_PPICLF)
    PPICLF_FFLAGS+=-module $(BUILD_ROOT_PPICLF) -I$(BUILD_ROOT_USER)
    PREPROC_FLAG = -cpp
else ifeq ($(TOOLCHAIN),INTEL)
    PPICLF_FFLAGS += -g -O0 -free -fpp  -Wp,-macro_expand=cpp
    PPICLF_USER_FFLAGS:=$(PPICLF_FFLAGS) -module $(BUILD_ROOT_USER) -I$(BUILD_ROOT_PPICLF)
    PPICLF_FFLAGS+=-module $(BUILD_ROOT_PPICLF) -I$(BUILD_ROOT_USER)
    PREPROC_FLAG = -preprocess-only
else
    PPICLF_FFLAGS += UNKNOWN_TOOLCHAIN_IN_PPICLF_MAKE
endif



# GSLIB LIBRARY
SOURCE_ROOT_GSLIB=$(INSTALL_LOCATION)/ppiclF/3rd_party/gslib
SOURCE_ROOT_GSLIB_OBJ=$(SOURCE_ROOT_GSLIB)/gslib/src
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include

SRCUSR =	$(SOURCE_ROOT_USER)/ppiclf_user_particle.F90 \
			$(SOURCE_ROOT_USER)/ppiclf_user_random.F90 \
            $(SOURCE_ROOT_USER)/ppiclf_user_MapProjPart.f90    \
			$(SOURCE_ROOT_USER)/ppiclf_m_user_data.F90 \
			$(SOURCE_ROOT_USER)/ppiclf_m_user_RFLUdata.F90 \
			$(SOURCE_ROOT_USER)/ppiclf_user_SetYdot.F90		\
            $(SOURCE_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_InitZero.F90     \
			$(SOURCE_ROOT_USER)/ppiclf_user_NearestNeighbor.f90 \
     		$(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.F90       \
			$(SOURCE_ROOT_USER)/ppiclf_m_user_ForceModels.F90 \
            $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.F90   \
            $(SOURCE_ROOT_USER)/ppiclf_user_Lift.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_Torque.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_HeatTransfer.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_Fluctuations.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_QuasiSteady.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_ViscousUnsteady.F90     \
            $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.F90   \
            $(SOURCE_ROOT_USER)/ppiclf_user_AddedMass.F90     \
			
UNUSED=     $(SOURCE_ROOT_USER)/ppiclf_user_unitproblems.F90     \
            $(SOURCE_ROOT_USER)/ppiclf_user_debug.F90     \

SRC=		$(SOURCE_ROOT_PPICLF)/ppiclf_m_types.F90	\
			$(SOURCE_ROOT_PPICLF)/ppiclf_user.F90	\
			$(SOURCE_ROOT_PPICLF)/ppiclf_initsolve.F90       \
			$(SOURCE_ROOT_PPICLF)/ppiclf_data.F90 \
            $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.F90	\
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_particledata.F90	\
			$(SOURCE_ROOT_PPICLF)/ppiclf_op.F90       \
			$(SOURCE_ROOT_PPICLF)/ppiclf_op_imp.F90 \
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_particle_ops.F90 \
		    $(SOURCE_ROOT_PPICLF)/ppiclf_m_comm.F90     \
            $(SOURCE_ROOT_PPICLF)/ppiclf_io.F90       \
            $(SOURCE_ROOT_PPICLF)/ppiclf_solve.F90    \
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_wrapped_types.F90 \
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_transfers.F90 \
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_wrapped_types_imp.F90 \

OBJ = 	$(patsubst $(SOURCE_ROOT_PPICLF)/%.F90,$(BUILD_ROOT_PPICLF)/%.o,$(SRC))
OBJUSR =  	$(patsubst $(SOURCE_ROOT_USER)/%.F90,$(BUILD_ROOT_USER)/%.o,$(SRCUSR))

PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h      \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PERFORMANCE.h
POS = $(SOURCE_ROOT_PPICLF)/PPICLF

.PHONY: default, test, testprog

default: libppicl
test: $(BUILD_ROOT_PPICLF)/ppiclf_m_transfers.o


testprog: $(BUILD_ROOT_PPICLF)/testprog
	@$(BUILD_ROOT_PPICLF)/testprog

$(BUILD_ROOT_PPICLF)/testprog: $(OBJ) $(SOURCE_ROOT_PPICLF)/testfile.F90
	$(FC) $(SOURCE_ROOT_PPICLF)/testfile.F90 -o $@ $(PPICLF_FFLAGS) $(BUILD_ROOT_PPICLF)/ppiclf_m_types.o $(BUILD_ROOT_USER)/ppiclf_user_particle.o

# default: $(BUILD_ROOT_USER)/ppiclf_user_imp.o

libppicl: $(BUILD_ROOT_PPICLF)/libppiclF.a

libppiclcore: $(BUILD_ROOT_PPICLF)/libppiclcore.a

# core library build rule
$(BUILD_ROOT_PPICLF)/libppiclcore.a: $(OBJ) makeThird
	ar crv $(BUILD_ROOT_PPICLF)/libppiclcore.a $(OBJ) $(SOURCE_ROOT_GSLIB_OBJ)/*.o

$(BUILD_ROOT_PPICLF)/testfile.f90: $(SOURCE_ROOT_PPICLF)/include/PPICLF_PARTMACROS.fypp

# module dependencies
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_io.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_particle_ops.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_transfers.o



$(BUILD_ROOT_PPICLF)/ppiclf_op_imp.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o

$(BUILD_ROOT_PPICLF)/ppiclf_m_particle_ops.o: $(BUILD_ROOT_USER)/ppiclf_user_particle.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_particle_ops.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o

$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o
$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o

$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o

$(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o
$(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_particledata.o
$(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o

$(BUILD_ROOT_PPICLF)/ppiclf_m_particledata.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_types.o

$(BUILD_ROOT_PPICLF)/ppiclf_data.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_wrapped_types.o

$(BUILD_ROOT_PPICLF)/ppiclf_m_wrapped_types.o: $(BUILD_ROOT_USER)/ppiclf_user_particle.o

$(BUILD_ROOT_PPICLF)/ppiclf_user.o: $(BUILD_ROOT_USER)/ppiclf_user_particle.o

$(BUILD_ROOT_PPICLF)/ppiclf_m_transfers.o: $(SOURCE_ROOT_PPICLF)/include/Generic_Procedures/ppiclf_alltoallv.subroutine.F90

# user file dependencies 
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_m_user_ForceModels.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_user_subbinMap.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o

$(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o: $(BUILD_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.o
$(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_Torque.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_Torque.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_Lift.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_Lift.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_HeatTransfer.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_HeatTransfer.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_Fluctuations.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_Fluctuations.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_ViscousUnsteady.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_ViscousUnsteady.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_AddedMass.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_AddedMass.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_AddedMass.o: $(BUILD_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.o

$(BUILD_ROOT_USER)/ppiclf_user_particle.o:	$(BUILD_ROOT_PPICLF)/ppiclf_m_types.o
$(BUILD_ROOT_USER)/ppiclf_user_particle.o:	$(BUILD_ROOT_PPICLF)/ppiclf_op.o

$(BUILD_ROOT_USER)/ppiclf_user_subbinMap.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o


$(BUILD_ROOT_PPICLF)/%.o : $(SOURCE_ROOT_PPICLF)/%.F90
	@echo called standard build on $<
	@mkdir -vp $(BUILD_ROOT_PPICLF)
	@$(FC) $< $(PREPROC_FLAG) -o $(BUILD_ROOT_PPICLF)/$*.i $(PPICLF_FFLAGS)
	@$(FC) $< -c -o $@ $(PPICLF_FFLAGS)
$(BUILD_ROOT_PPICLF)/%.f90 : $(SOURCE_ROOT_PPICLF)/%.F90
	@echo called standard preprocess on $<
	@mkdir -vp $(BUILD_ROOT_PPICLF)
	@$(FC) $< $(PREPROC_FLAG) -o $@ $(PPICLF_FFLAGS)

$(BUILD_ROOT_PPICLF)/ppiclf_io.f90 $(BUILD_ROOT_PPICLF)/ppiclf_io.o : $(SOURCE_ROOT_PPICLF)/ppiclf_io.F90
	@echo called custom ppiclf_io build
	@mkdir -vp $(BUILD_ROOT_PPICLF)
	@$(FC) $(SOURCE_ROOT_PPICLF)/ppiclf_io.F90 $(PREPROC_FLAG) -o $(BUILD_ROOT_PPICLF)/ppiclf_io.i $(PPICLF_FFLAGS) -P -g1
	@sed 's/RE_NEWLINE/\n/g; s/RE_AMPERSAND/\&/g' $(BUILD_ROOT_PPICLF)/ppiclf_io.i > $(BUILD_ROOT_PPICLF)/ppiclf_io.f90
	@echo $@
	@echo $<
	@$(if $(and $(findstring ./$@,$(BUILD_ROOT_PPICLF)/ppiclf_io.o),$(findstring $(BUILD_ROOT_PPICLF)/ppiclf_io.o,./$@)),@$(FC) $(BUILD_ROOT_PPICLF)/ppiclf_io.f90 -c -o $@ $(PPICLF_FFLAGS))
# mpifort source/ppiclf_data.f90 -cpp -o ppiclf_data.o -Isource/include -I../user_files -ffree-line-length-none





# full library build rule, adds the user files to the core archive generated above
# this separation allows for the core library to be build without the user files (except the .h file)
$(BUILD_ROOT_PPICLF)/libppiclF.a: $(BUILD_ROOT_PPICLF)/libppiclcore.a $(OBJUSR)
	cp $(BUILD_ROOT_PPICLF)/libppiclcore.a $(BUILD_ROOT_PPICLF)/libppiclF.a
	ar rv $(BUILD_ROOT_PPICLF)/libppiclF.a $(OBJUSR)



$(BUILD_ROOT_USER)/%.o : $(SOURCE_ROOT_USER)/%.F90
	@echo called standard build on $<
	@mkdir -vp $(BUILD_ROOT_USER)
	@$(FC) $< -c -o $@ $(PPICLF_USER_FFLAGS)



# $(BUILD_ROOT_PPICLF)/%.o : $(BUILD_ROOT_USER)/%.f90
# 	@mkdir -vp $(BUILD_ROOT_USER)
# 	@$(FC) $< $(PREPROC_FLAG) -o $(BUILD_ROOT_USER)/$*.i $(PPICLF_FFLAGS)
# 	@$(FC) $< -c -o $@ $(PPICLF_FFLAGS)






# Old Make commands
defaultold: makeThird preProcess getObjs libObjs 

libObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.o
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(SOURCE_ROOT_PPICLF)/ppiclf.o $(SOURCE_ROOT_GSLIB_OBJ)/*.o
	@echo "                       "
	@echo "***********************"
	@echo "*** LIBRARY SUCCESS ***"
	@echo "***********************"
	@echo "                       "

getObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.f
	$(FC) $(FFLAGS) -c $(SOURCE_ROOT_PPICLF)/ppiclf.f $(GSLIB_IFLAGS) 
	mv *.o $(SOURCE_ROOT_PPICLF)
	@echo "                              "
	@echo "******************************"
	@echo "*** LIBRARY OBJECT SUCCESS ***"
	@echo "******************************"
	@echo "                              "

preProcess: $(PRE)
	$(CC) $(CPFLAGS) $(PRE) > $(POS)_tmp
	grep "^[^#]" $(POS)_tmp > $(POS)
	rm $(POS)_tmp
	echo '#include "PPICLF_USER.h"' > $(SOURCE_ROOT_PPICLF)/ppiclf.f
	echo '#include "PPICLF_STD.h"' >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	cat $(SRC)  >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	@echo "                       "
	@echo "********************"
	@echo "*** PREPROCESSED ***"
	@echo "********************"
	@echo "                       "
	

makeThird: $(SOURCE_ROOT_GSLIB)/install
	cd $(SOURCE_ROOT_GSLIB); \
	./install $(CC) $(FC)
	@echo "                       "
	@echo "***********************"
	@echo "*** INSTALLED GSLIB ***"
	@echo "***********************"
	@echo "                       "

clean:             
	@rm -rf $(SOURCE_ROOT_GSLIB)/gslib                       
	@rm -rf $(SOURCE_ROOT_GSLIB)/lib                         
	@rm -rf $(SOURCE_ROOT_GSLIB)/include                     
	@rm -rf $(SOURCE_ROOT_GSLIB)/*.tar.gz                    
	@rm -rf $(INSTALL_LOCATION)/short_tests/*.pyc            
	@rm -rf $(INSTALL_LOCATION)/short_tests/lib/*.pyc        
	@rm -rf $(INSTALL_LOCATION)/short_tests/test1/test1.log*
	@rm -rf $(INSTALL_LOCATION)/short_tests/test1/*.vtu
	@rm -rf $(BUILD_ROOT_PPICLF)
	@mkdir $(BUILD_ROOT_PPICLF)
	@mkdir $(BUILD_ROOT_USER)
	@mkdir $(BUILD_ROOT_USER)/Added_Mass
	@cd $(INSTALL_LOCATION)/ppiclF/short_tests/test1/; \
	make clean


.SECONDARY: $(patsubst $(SOURCE_ROOT_PPICLF)/%.fpp,$(BUILD_ROOT_PPICLF)/%.f90,$(SRC)) $(patsubst $(SOURCE_ROOT_USER)/%.fpp,$(BUILD_ROOT_USER)/%.f90,$(SRCUSR))

FORCE:
