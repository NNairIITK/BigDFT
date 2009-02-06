#ifndef _locpot_
#define _locpot_

template<typename T>
int magicfilterpot(int n1,int n2, int n3,
		   T *psi,
		   T *work,
		   T *pot,
		   T *epot);


template<typename T>
int mf1d(int ndat, int n3,
	 T *psi,
	 T *out);

//this function is called by fortran programs
extern "C"
void localpotential_(int *n1,
		     int *n2,
		     int *n3,
		     float **psi,
		     float **work,
		     float **pot,
		     float *epot); 


extern "C" 
void localpotentiald_(int *n1,
		     int *n2,
		     int *n3,
		     double **psi,
		     double **work,
		     double **pot,
		     double *epot);



#endif
