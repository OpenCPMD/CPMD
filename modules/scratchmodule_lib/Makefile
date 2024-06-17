AR = ar -r
RANLIB = ranlib
#AR     = /usr/bin/ar ruv
#RANLIB = /usr/bin/ranlib
FFLAGS =  -D_OPENMP_OFFLOAD -qopenmp -fopenmp-targets=spir64 -g -traceback -cpp -axCORE-AVX2,CORE-AVX512 -mavx -unroll -funroll-loops -O2 -IPO  #-O0 -check  bounds,assume,contiguous,format,pointers,shape,stack,udio_iostat,uninit   #
#FFLAGS =  -cpp -D_DEBUG -warn all -O0 -check bounds,assume,contiguous,format,pointers,shape,stack,udio_iostat,uninit   #-g -traceback -cpp -axCORE-AVX2 -mavx -unroll -funroll-loops -O2 -IPO  #-O0 -check  bounds,assume,contiguous,format,pointers,shape,stack,udio_iostat,uninit   #
#FFLAGS =  -g -fbacktrace  -cpp -unroll -funroll-loops -O2   #-O0 -check  bounds,assume,contiguous,format,pointers,shape,stack,udio_iostat,uninit   #
CPPFLAGS =

#FC = gfortran
FC = ifx

SCRATCHMODULE_LIB = libscratch_module.a
SRC_SCRATCHMODULE = data_managment_utils.f90 segment_managment_utils.f90 pool_managment_utils.f90 scratch_interface.f90
.SUFFIXES: .f90

%.mod:

%.o:    %.f90 
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(@F:.o=.f90)


OBJ_SCRATCHMODULE     = $(SRC_SCRATCHMODULE:%.f90=%.o)
$(SCRATCHMODULE_LIB): $(OBJ_SCRATCHMODULE)
	$(AR) $(SCRATCHMODULE_LIB) $(OBJ_SCRATCHMODULE)
	$(RANLIB) $(SCRATCHMODULE_LIB)


data_managment_utils.mod: data_managment_utils.o
	@true
segment_managment_utils.mod: segment_managment_utils.o
	@true
pool_managment_utils.mod: pool_managment_utils.o
	@true
segment_managment_utils.o: data_managment_utils.mod segment_managment_utils.f90
pool_managment_utils.o:  segment_managment_utils.mod pool_managment_utils.f90
scratch_interface.o: pool_managment_utils.mod scratch_interface.f90


clean: 
	rm -f *.o *.mod *.a
