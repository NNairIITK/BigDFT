#ifndef _kinetich_
#define _kinetich_





extern "C" 
void kineticterm_(int *n1,
		  int *n2,
		  int *n3,
		  float *hx,
		  float *hy,
		  float *hz,
		  float *c,
		  float **x,
		  float **y,
		  float **workx,
		  float **worky,
		  float *ekin);

int kineticfilter(int n1,
		  int n2, 
		  int n3,
		  float h1,
		  float h2,
		  float h3,
		  float c,
		  float *x,
		  float *workx,
		  float *y,
		  float *worky,
		  float *ekin);



#endif
