include ../../Makefile.examples.am

examplesdir=@prefix@/examples/c/add_array

if BUILD_EXAMPLES
buildexe = add_array
endif

examples_DATA = \
	GPU.config \
	add_array.c \
	type.h \
	add_kernel.cu \
	Makefile.standalone \
	$(buildexe)

EXTRA_DIST = \
	GPU.config \
	add_array.c \
	type.h \
	add_kernel.cu

MPICC = mpicc

#all: add_array

add_array: add_array.o add_kernel.o
	$(MPICC) add_array.o add_kernel.o -o $@ $(MYLDFLAGS) $(MYLIBS)

add_array.o: add_array.c type.h
	$(MPICC) $(MYCFLAGS) -c add_array.c -o $@

add_kernel.o: add_kernel.cu type.h
	$(NVCC) $(NVCCFLAGS) -c add_kernel.cu -o $@

clean:
	rm -f *.o *~ add_array
