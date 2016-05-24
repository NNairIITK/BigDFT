#ifndef __compressh__
#define __compressh__

template<typename T>
int uncompressgpu(int n1, int n2, int n3,
		  T *psicf,T *psig, int *keys);

template<typename T>
int compressgpu(int n1, int n2, int n3, 
		T *psig,T *psicf, int *keys);


#endif
