#ifndef   	REDUCTION_H
#define   	REDUCTION_H

float reducearrays(int n,
		   int ndat,
		   float *psi,
		   float *vpsi,
		   float *epot);

float reducearrays_d(int n,
		     int ndat,
		     float *psi,
		     float *vpsi,
		     double *epot);


#endif
