FC               = mpifort
CC               = mpicc  
CPFLAGS         += -E

PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)/include -I$(SOURCE_ROOT_USER) -I$(SOURCE_ROOT_USER)/include 
PPICLF_FFLAGS+=$(PPICLF_IFLAGS) -cpp -Wall -ffree-line-length-none -Wno-missing-include-dirs -fmax-errors=5 -O3 -g
PPICLF_FFLAGS+=-I$(BUILD_ROOT_PPICLF) -I$(BUILD_ROOT_USER)

# PPICLF LIBRARY
LIBPPICL_ROOT=../
SOURCE_ROOT_PPICLF= $(LIBPPICL_ROOT)/ppiclF/source
BUILD_ROOT_PPICLF=$(LIBPPICL_ROOT)/ppiclF/build
SOURCE_ROOT_USER=$(LIBPPICL_ROOT)/user_files
BUILD_ROOT_USER=$(BUILD_ROOT_PPICLF)/user
LFLAGS+=-L$(BUILD_ROOT_PPICLF) -lppiclF

      
MODULE_SRC=	unitTestSrc/UT_m_data.F90 \
			unitTestSrc/UT_m_setup.F90 \

MODULE_OBJ= $(patsubst unitTestSrc/%.F90,%.o,$(MODULE_SRC))

EXE_SRC = 	unitTestSrc/UT_CreateBin.F90 \
			unitTestSrc/UT_InterpProj.F90 \
			unitTestSrc/UT_P2PNN.F90 \

UNUSED=			unitTestSrc/UT_CodePerformance.F90 \

EXE=$(patsubst %.F90,%_UnitTest,$(patsubst unitTestSrc/UT_%,%,$(EXE_SRC)))

# Make commands
default: $(EXE)

%_UnitTest: unitTestSrc/UT_%.F90 $(MODULE_OBJ)
	@echo Building $@ from $<
	@$(FC) $< -o $@ $(PPICLF_FFLAGS) $(MODULE_OBJ) $(LFLAGS)

%.o: unitTestSrc/%.F90
	@echo Building module $@ from $<
	@$(FC) $< -c -o $@ $(PPICLF_FFLAGS)

# Cleaning stuff
clean:
	@rm -f $(MODULE_OBJ) $(EXE)
	@rm -f *.mod
	@rm -f *.smod

.SECONDARY: $(MODULE_OBJ)
