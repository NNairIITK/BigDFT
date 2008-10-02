/****c* CUDA/convSeria.h
** 
** AUTHOR
**  Matthieu Ospici
** 
** CHANGELOG
**  Started on  Wed Apr  9 15:36:25 2008 Matthieu Ospici
**
**  Last update Wed Apr  9 15:36:25 2008 Matthieu Ospici
**
** SOURCE
*/

#ifndef   	CONVSERIA_H_
# define   	CONVSERIA_H_

//#include "class_convolution.h"
class convSeria
{
public:
  static void convSerial(unsigned int n1,unsigned int n2,double *tab_out, double *tab_in,const double *in_f,int nf);
  static void convSerial_s(unsigned int n1,unsigned int n2,float *tab_out, float *tab_in,const float *in_f,int nf);
  static void convSerial_sd(unsigned int n1,unsigned int n2,float *tab_out, float *tab_in,const double *in_f,int nf);


  /*  static void conv3dCPU_s(multiTab_t* m_dataIn,
			  multiTab_t* m_dataOut,
			  const float *f_data,
			  int fsize);*/

  
  //  static void conv3DCPU(convolution& in, convolution& out);
};




#endif 	    /* !CONVSERIA_H_ */

/****/
