# Shared build configuration for rembo1 (dependency wiring).
#
# Loaded after the compiler and machine fragments. References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# fesm-utils provides the shared utility modules (ncio, nml, interp1D,
# gaussian_filter, ncio_transpose). rembo1 relies on these directly instead of
# vendoring its own copies. Listed before coordinates so its modules win on the
# include search path where names overlap.
FESMUTILSROOT = fesm-utils/utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

# coordinates is linked into libs/coordinates (see this package's
# [[package.links]] in configme: rembo1 -> coordinates at libs/coordinates).
COORDROOT = libs/coordinates
INC_COORD = -I${COORDROOT}/libcoordinates/include
LIB_COORD = -L${COORDROOT}/libcoordinates/include -lcoordinates

LISROOT = fesm-utils/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib/ -llis

# rembo1 needs the fesm-utils and coordinates includes in its compile flags
# (fesm-utils first so its modules take precedence over coordinates' copies).
FFLAGS += $(INC_FESMUTILS)
FFLAGS += $(INC_COORD)

ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    LISROOT = fesm-utils/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib/ -llis

    FFLAGS += $(FFLAGS_OPENMP)
endif

# rembo build outputs, for downstream linking by an orchestrator.
REMBOROOT = ${CURDIR}
INC_REMBO = -I${REMBOROOT}/librembo/include
LIB_REMBO = -L${REMBOROOT}/librembo/include -lrembo

LFLAGS_EXTRA ?= -Wl,-zmuldefs
LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_COORD) $(LIB_LIS) $(LFLAGS_EXTRA)
