include ../../Makefile.examples.am

examplesdir=@prefix@/examples/c/simple

examples_DATA = \
Makefile.standalone \
square.cu \
square_not_shared.cu \
GPU.config \
square \
square_not_shared


all: square square_not_shared

square: square.cu
	$(NVCC) $(NVCCFLAGS) $(MYCFLAGS) $^ -o $@ $(MYLDFLAGS) $(MYLIBS)

square_not_shared: square_not_shared.cu
	$(NVCC) $(NVCCFLAGS) $(MYCFLAGS) $^ -o $@ $(MYLDFLAGS) $(MYLIBS)

clean:
	rm -f *.o *~ square square_not_shared

