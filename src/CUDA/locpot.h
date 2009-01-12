#ifndef _locpot_
#define _locpot_


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
