.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = librembo/include
bindir = librembo/bin
libdir = libs
testdir = tests

# Command-line options at make call
debug ?= 0

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

FC = gfortran

INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf

COORDROOT = /Users/robinson/models/EURICE/coordinates/libcoordinates
INC_COORD = -I${COORDROOT}/include
LIB_COORD = -L${COORDROOT}/include -lcoordinates

MKLROOT = /opt/intel/mkl
INC_MKL = -I${MKLROOT}/include
# For sequential running (no tbb library)
LIB_MKL =  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

LISROOT = /Users/robinson/apps/lis/lis
INC_LIS = -I${LISROOT}/include 
LIB_LIS = -L${LISROOT}/lib/ -llis

FFLAGS  = -ffree-line-length-none -I$(objdir) -J$(objdir) $(INC_COORD)
LFLAGS  = $(LIB_NC) $(LIB_COORD) $(LIB_LIS)
DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg


# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of rembo rules and source files
##
###############################################

include config/Makefile_rembo.mk

# GLOBAL RULES ####
$(objdir)/nml.o : $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
	
$(objdir)/ncio.o : $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/ncio_transpose.o : $(libdir)/ncio_transpose.f90 $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o : $(libdir)/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
##
## Compilation of complete programs
##
###############################################

# MAIN PROGRAMS ####

rembo: $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/ncio_transpose.o $(objdir)/interp1D.o \
       $(rembo_base) $(srcdir)/standalone_main.f90
	$(FC) $(DFLAGS) $(FFLAGS) -o librembo/bin/rembo.x $^ $(LFLAGS)
	@echo " "
	@echo "    librembo/bin/rembo.x is ready."
	@echo " "

# Static libraries 
yelmo-static: $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ar rc $(objdir)/libyelmo.a $(yelmo_libs) $(yelmo_physics) $(yelmo_base) $(yelmo_tests)
	ranlib $(objdir)/libyelmo.a
	@echo " "
	@echo "    $(objdir)/libyelmo.a is ready."
	@echo " "

yelmo_test : yelmo-static
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/yelmo_test.x yelmo_test.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_test.x is ready."
		@echo " "

yelmo_benchmarks : yelmo-static
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/yelmo_benchmarks.x tests/yelmo_benchmarks.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_benchmarks.x is ready."
		@echo " "

yelmo_halfarmed : yelmo-static
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/yelmo_halfarmed.x tests/yelmo_halfarmed.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_halfarmed.x is ready."
		@echo " "

yelmo_mismip : yelmo-static
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/yelmo_mismip.x tests/yelmo_mismip.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_mismip.x is ready."
		@echo " "

yelmo_stommel : yelmo-static $(objdir)/stommel.o
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/yelmo_stommel.x yelmo_stommel.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o $(objdir)/stommel.o
		@echo " "
		@echo "    yelmo_stommel.x is ready."
		@echo " "

yelmo_optbeta : yelmo-static
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/yelmo_optbeta.x yelmo_optbeta.f90 \
			-L${CURDIR}/libyelmo/include -lyelmo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    yelmo_optbeta.x is ready."
		@echo " "

test_icetemp : $(yelmo_thermo)
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/test_icetemp.x tests/test_icetemp.f90 \
			$(LFLAGS) $(yelmo_thermo) $(objdir)/nml.o
		@echo " "
		@echo "    test_icetemp.x is ready."
		@echo " "

check: $(yelmo_libs) $(yelmo_tests)
	@echo $(yelmo_libs)
	@echo $(yelmo_tests)
	@echo ""

# NOTE: eventually ncio should be extracted from all external libraries (like coord?), and then
# all external libraries should be linked by adding them to the end of the program compilation 
# command, instead of just `$(objdir)/nml.o` as it is now. 

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make yelmo : compiles test_yelmo.x, yelmo on a given domain."
	@echo " make clean     : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
