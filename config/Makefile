.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir   = src
objdir   = librembo/include
bindir   = librembo/bin
libdir   = libs

# Command-line options at make call
debug  ?= 0
openmp ?= 0 

## COMPILER CONFIGURATION ##
# Loaded from a config/<host> fragment by `python config.py config/<host>`,
# or assembled by `configme`. The fragment defines the compiler (FC,
# FFLAGS, FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and
# any host-specific overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC, and LFLAGS_EXTRA set there. Defines the coordinates
# and LIS paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

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
# nml, ncio, ncio_transpose, interp1D and gaussian_filter are provided by
# fesm-utils (libfesmutils); rembo1 no longer vendors its own copies.
rembo_libs =

###############################################
##
## Compilation of complete programs
##
###############################################

# MAIN PROGRAMS ####

# Static library compilation
rembo-static: $(rembo_libs) $(rembo_base)
	ar rc $(objdir)/librembo.a $(rembo_libs) $(rembo_base)
	ranlib $(objdir)/librembo.a
	@echo " "
	@echo "    $(objdir)/librembo.a is ready."
	@echo " "

rembo: rembo-static
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_REMBO) -o librembo/bin/test_rembo.x \
		test_rembo.f90 $(LIB_REMBO) $(LFLAGS)
	@echo " "
	@echo "    librembo/bin/test_rembo.x is ready."
	@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make rembo : compiles rembo.x"
	@echo " make clean : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
