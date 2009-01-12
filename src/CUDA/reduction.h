#ifndef   	REDUCTION_H
#define   	REDUCTION_H

/*float reducearrays(int n,
		   int ndat,
		   float *psi,
		   float *vpsi,
		   float *epot);*/


template<typename T>
float reducearrays(int n,
		   int ndat,
		   T *psi,
		   T *vpsi,
		   T *epot);

float reducearrays_d(int n,
		     int ndat,
		     float *psi,
		     float *vpsi,
		     double *epot);


template<typename T>
float reducearrays(int n,
		   int ndat,
		   T *psi,
		   T *vpsi,
		   T *epot);

#endif
