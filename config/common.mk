# Shared build configuration for rembo1 (dependency wiring).
#
# Loaded after the compiler and machine fragments. References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# coordinates is linked into libs/coordinates (see this package's
# [[package.links]] in configme: rembo1 -> coordinates at libs/coordinates).
COORDROOT = libs/coordinates
INC_COORD = -I${COORDROOT}/libcoordinates/include
LIB_COORD = -L${COORDROOT}/libcoordinates/include -lcoordinates

LISROOT = fesm-utils/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib/ -llis

# rembo1 needs the coordinates include in its compile flags.
FFLAGS += $(INC_COORD)

ifeq ($(openmp), 1)
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
LFLAGS = $(LIB_NC) $(LIB_COORD) $(LIB_LIS) $(LFLAGS_EXTRA)
