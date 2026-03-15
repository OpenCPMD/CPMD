AR     = /usr/bin/ar ruv
RANLIB = /usr/bin/ranlib
FFLAGS =  -g -traceback -qopenmp -axCORE-AVX2,CORE-AVX512 -mavx -unroll -funroll-loops -O2
CPPFLAGS =

FC = mpiifx

GRIMME_LIB = libgrimmevdw.a
SRC_GRIMME = vdw_calculator.F90 vdw_interface.F90 vdw_param.F90
.SUFFIXES: .F90 .f90

%.mod:

%.o:    %.F90 
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(@F:.o=.F90)


OBJ_GRIMME     = $(SRC_GRIMME:%.F90=%.o)
$(GRIMME_LIB): $(OBJ_GRIMME)
	$(AR) $(GRIMME_LIB) $(OBJ_GRIMME)
	$(RANLIB) $(GRIMME_LIB)


vdw_interface.mod: vdw_interface.o
	@true
vdw_param.mod: vdw_param.o
	@true
vdw_interface.o: vdw_param.mod vdw_calculator.mod vdw_interface.F90
vdw_calculator.o: vdw_param.mod vdw_calculator.F90


clean: 
	rm -f *.o *.a
