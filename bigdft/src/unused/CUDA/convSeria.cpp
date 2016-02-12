#include <stdio.h>

#include "convSeria.h"

//#include  "structUtil.h"
void convSeria::convSerial(unsigned int n1,unsigned int n2,double *tab_out, double *tab_in,const double *in_f,int nf)
{  //ng*N  N*ng  nf
  
  for(unsigned int i=0;i < n1;++i)
    {
      for(unsigned int j=0;j <  n2;++j)
	{
	  double tmp = 0;
	  for(int k=0,mod = j ;k < nf;++k)
	    {
	      if(mod >= (int)n2)
		mod = 0;
	      
	      //  tmp += tab_in[i*n2 + mod]*in_f[k];
	      tmp += tab_in[i + mod*n1]*in_f[k];


	      // tmp += in_g->tab[i*in_g->n2 + (k+j)%( in_g->n2_real)]*in_f[k];
	   
		 ++mod;
		 //	 printf("i : %i, j : %i, mod %i,  val : %.40f\n",i,j,mod,tmp);
	    }
	  
	  //  tab_out[j*n1 + i] = tmp;
	  tab_out[j + i*n2] = tmp; 

 
	}
    }
}


void convSeria::convSerial_s(unsigned int n1,unsigned int n2,float *tab_out, float *tab_in,const float *in_f,int nf)
{  //ng*N  N*ng  nf
  
  for(unsigned int i=0;i < n1;++i)
    {
      for(unsigned int j=0;j <  n2;++j)
	{
	  float tmp = 0;
	  for(int k=0,mod = j ;k < nf;++k)
	    {
	      if(mod >= (int)n2)
		mod = 0;
	      
	   
	      //  tmp += tab_in[i*n2 + mod]*in_f[k];
	      tmp += tab_in[i + mod*n1]*in_f[k];
	      // printf("%i %i %i\n",mod,i,j,k);
	      ++mod;
		
	    }

	   
	  // tab_out[j*n1 + i] = tmp;
	   tab_out[j + i*n2] = tmp;
	   // printf("---\n");
	}
    }
}



void convSeria::convSerial_sd(unsigned int n1,unsigned int n2,float *tab_out, float *tab_in,const double *in_f,int nf)
{  //ng*N  N*ng  nf
 
  for(unsigned int i=0;i < n1;++i)
    {
      for(unsigned int j=0;j <  n2;++j)
	{
	  double tmp = 0;
	  for(int k=0,mod = j ;k < nf;++k)
	    {
	      if(mod >= (int)n2)
		mod = 0;
	      
	   
	      //  tmp += tab_in[i*n2 + mod]*in_f[k];
	      tmp += tab_in[i + mod*n1]*in_f[k];

	      ++mod;
		
	    }

	   
	  //  tab_out[j*n1 + i] = tmp;
	  tab_out[j + i*n2] = tmp;
 
	}
    }
}

/*void  convSeria::conv3dCPU_s(multiTab_t* m_dataIn,
				   multiTab_t* m_dataOut,
				   const float *f_data,
				   int fsize)
{
  unsigned int n1 = m_dataIn->n1;
  unsigned int n2 = m_dataIn->n2 * m_dataIn->n3; 
  convSerial_s(n1,n2,m_dataOut->tab,m_dataIn->tab,f_data,fsize);

  
  n1 = m_dataIn->n2;
  n2 = m_dataIn->n3 * m_dataIn->n1;
  convSerial_s(n1,n2,m_dataIn->tab,m_dataOut->tab,f_data,fsize);


  n1 = m_dataIn->n3;
  n2 = m_dataIn->n1 * m_dataIn->n2;
  convSerial_s(n1,n2,m_dataOut->tab,m_dataIn->tab,f_data,fsize);

}

void convSeria::conv3DCPU(convolution& in, convolution& out)
{



  multiTab_t m_dataIn;
  in.extractMultiTab(m_dataIn); 
  
  multiTab_t m_dataOut;
  out.extractMultiTab(m_dataOut);


  float *fdata = in.fdata;
  int fsize = in.fdata_size;



  conv3dCPU_s(&m_dataIn,
	      &m_dataOut,
	      fdata,
	      fsize);
  

}*/




