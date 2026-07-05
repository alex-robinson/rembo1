# Shared build configuration for rembo1 (dependency wiring).
#
# Loaded after the compiler and machine fragments. References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# fesm-utils provides all shared modules rembo1 needs: the utility modules
# (ncio, nml, interp1D, gaussian_filter, ncio_transpose) and the coordinate/grid
# modules (formerly the standalone `coordinates` library, now folded into
# libfesmutils). rembo1 links this directly instead of vendoring its own copies.
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

LISROOT = fesm-utils/lis/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib -llis

# rembo1 needs the fesm-utils includes in its compile flags.
FFLAGS += $(INC_FESMUTILS)

ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    LISROOT = fesm-utils/lis/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib -llis

    FFLAGS += $(FFLAGS_OPENMP)
endif

# rembo build outputs, for downstream linking by an orchestrator.
REMBOROOT = ${CURDIR}
INC_REMBO = -I${REMBOROOT}/librembo/include
LIB_REMBO = -L${REMBOROOT}/librembo/include -lrembo

LFLAGS_EXTRA ?= -Wl,-zmuldefs
LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_LIS) $(LFLAGS_EXTRA)
