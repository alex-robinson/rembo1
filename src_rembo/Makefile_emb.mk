



# EMB related source files
emb_dir   = src_rembo/
emb_files = emb_global.f90 emb_pdesolver.f90 projector.f90 \
            emb_functions.f90 rembo_functions.f90 \
            smb_itm.f90 \
            rembo_main.f90 climate.f90 

# Add the appropriate subdirectory prefixes on to the files and
# change the file suffix (eg, .o object files)  
emb_SRC = $(patsubst %, $(emb_dir)%, $(emb_files) )
TMP = $(patsubst %.f90, %.o, $(emb_files) )
emb_OBJ = $(patsubst %, $(o1)%, $(TMP) )


# EMB RULES ####
$(o1)emb_global.o : $(emb_dir)emb_global.f90 $(o1)ncio.o
	$(F77) $(LDFLAGS) -c -o $@ $< 
$(o1)emb_pdesolver.o : $(emb_dir)emb_pdesolver.f90 $(o1)emb_global.o 
	$(F77) $(LDFLAGS) -c -o $@ $<
$(o1)projector.o : $(emb_dir)projector.f90 $(o1)emb_global.o 
	$(F77) $(LDFLAGS) -c -o $@ $<
$(o1)emb_functions.o : $(emb_dir)emb_functions.f90 $(o1)emb_global.o $(o1)ncio.o 
	$(F77) $(LDFLAGS) -c -o $@ $<
$(o1)rembo_functions.o : $(emb_dir)rembo_functions.f90 $(o1)emb_global.o
	$(F77) $(LDFLAGS) -c -o $@ $<
$(o1)smb_itm.o : $(emb_dir)smb_itm.f90 $(o1)emb_global.o $(o1)emb_functions.o
	$(F77) $(LDFLAGS) -c -o $@ $<
$(o1)rembo_main.o : $(emb_dir)rembo_main.f90 $(o1)rembo_functions.o $(o1)emb_global.o \
$(o1)emb_functions.o $(o1)emb_pdesolver.o $(o1)projector.o
	$(F77) $(LDFLAGS) -c -o $@ $<
$(o1)climate.o : $(emb_dir)climate.f90 $(o1)emb_global.o $(o1)emb_functions.o $(o1)rembo_main.o 
	$(F77) $(LDFLAGS) -c -o $@ $<
# DONE EMB RULES


# EMB DRIVER related source files
embdriver_OBJ  = $(o1)sinsol_orbit.o
embdriver_PROG = $(emb_dir)standalone_main.f90

$(o1)sinsol_orbit.o : $(emb_dir)sinsol_orbit.f
	$(F77) $(FFLAGS) -c -o $@ $<





