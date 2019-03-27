.SUFFIXES: .f .F .F90 .f90 .o
.SHELL: /bin/sh

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make rembo      : compiles the main program rembo.x"
	@echo " make clean      : cleans all object files
	@echo ""
	
# Driver related source files and directories
maindir = ./
o1      = $(maindir).obj/
objdir  = $(o1)

netcdf_inc = /opt/local/include
netcdf_lib = /opt/local/lib -lnetcdff -lnetcdf
#netcdf_inc = /home/fispalma25/apps/netcdf/netcdf/include
#netcdf_lib = /home/fispalma25/apps/netcdf/netcdf/lib -lnetcdf -lnetcdff
netcdf_inc_ifort = /home/robinson/apps/netcdf/netcdf/include
netcdf_lib_ifort = /home/robinson/apps/netcdf/netcdf/lib -lnetcdf

eurice_root = /Users/robinson/models/EURICE
#eurice_root = /p/projects/tumble/robinson/EURICE
tracer_inc  = ${eurice_root}/tracer/libtracer/include
tracer_lib  = ${eurice_root}/tracer/libtracer/include

ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I$(netcdf_inc_ifort)
	LFLAGS		 = -L${tracer_lib} -ltracer -L$(netcdf_lib_ifort)

	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all
	    # -w 
	else
	    DFLAGS   = -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I$(netcdf_inc)
	#LFLAGS		 = -L${tracer_lib} -ltracer -L$(netcdf_lib)
	LFLAGS		 = -L$(netcdf_lib)

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
	else
	    DFLAGS   = -O3
	endif
endif

LDFLAGS = $(FLAGS) $(DFLAGS)
LIB     = $(LFLAGS)

# Include emb related files and rules
include src_rembo/Makefile_emb.mk

.PHONY : files
files:
	@echo $(src_sico)

# MAIN PROGRAMS ####

rembo: $(o1)nml.o $(o1)ncio.o $(o1)ncio_transpose.o $(o1)parameters.o $(o1)exchange.o $(o1)interp1D.o \
       $(emb_OBJ) $(embdriver_OBJ) $(embdriver_PROG)
	$(F77) $(LDFLAGS) -o rembo.x $^ $(LIB)
	@echo " "
	@echo "    rembo.x is ready."
	@echo " "

# DONE MAIN PROGRAMS ####

# GLOBAL RULES ####
$(o1)parameters.o : $(main_dir)parameters.f90
	$(F77) $(LDFLAGS) -c -o $@ $<

$(o1)exchange.o : $(main_dir)exchange.f90 $(o1)parameters.o
	$(F77) $(LDFLAGS) -c -o $@ $<

$(o1)nml.o : $(main_dir)nml.f90
	$(F77) $(LDFLAGS) -c -o $@ $<
	
$(o1)ncio.o : $(main_dir)ncio.f90
	$(F77) $(LDFLAGS) -c -o $@ $<

$(o1)ncio_transpose.o : $(main_dir)ncio_transpose.f90 $(o1)ncio.o
	$(F77) $(LDFLAGS) -c -o $@ $<

$(o1)interp1D.o : $(main_dir)interp1D.f90
	$(F77) $(LDFLAGS) -c -o $@ $<
	
# -Dsimenv_run_char="' '"
clean:
	rm -f rembo.x $(emb_OBJ) $(embdriver_OBJ) $(o1)parameters.o $(o1)ncio.o $(o1)exchange.o \
						   $(o1)interp1D.o $(o1)*.mod

