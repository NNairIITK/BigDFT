include ../../Makefile.examples.am

examplesdir=@prefix@/examples/fortran/add_array

if BUILD_EXAMPLES
buildexe = add_array
endif

examples_DATA = \
	GPU.config \
	add_array.F90 \
	type.h \
	add_kernel.cu \
	binding.c \
	Makefile.standalone \
	$(buildexe)

EXTRA_DIST = \
	GPU.config \
	add_array.F90 \
	type.h \
	add_kernel.cu \
	binding.c

MPIFC = mpif90 -ffree-line-length-none -fno-underscoring

#all: add_array

add_array: add_array.o add_kernel.o binding.o
	$(MPIFC) $^ -o $@ $(MYLDFLAGS) $(MYLIBS)

add_array.o: add_array.F90 
	$(MPIFC) $(MYCFLAGS) -c $< -o $@

binding.o: binding.c type.h
	$(CC) $(MYCFLAGS) -c $< -o $@

add_kernel.o: add_kernel.cu type.h
	$(NVCC) $(NVCCFLAGS) -c add_kernel.cu -o $@

clean:
	rm -f *.o *~ add_array
