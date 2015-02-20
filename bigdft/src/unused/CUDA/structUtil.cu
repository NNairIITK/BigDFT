/****u* CUDA/structUtil.cu
**
** AUTHOR
**  Matthieu Ospici
**
** CHANGELOG
** Started on  Mon Feb 11 13:39:47 2008 Matthieu Ospici
**
** Last update Sun May 12 01:17:25 2002 Speed Blue
**
** SOURCE
*/

#include <stdio.h>
#include "structUtil.h"

// =================================
/*void printTab(const  multiTab_t *toPrint,int ib ,int jb)
{

  int iStop,jStop;

  if(ib == 0)
    iStop = toPrint->n1;
  else
    iStop = ib;

  if(jb == 0)
    jStop = toPrint->n2;
  else
    jStop = jb;

  for(unsigned int i=0;i < iStop;++i)
    {
      printf("| ");
      for(unsigned int j=0; j < jStop; ++j)
	{
	  // printf(" %f ",toPrint->tab[i*toPrint->n2 +j]);    
	  printf(" %f ",getTab(toPrint,i,j));
	}
      printf(" |\n");
    }
printf("\n");
}

//===================================
bool compareTab(const multiTab_t *a, const multiTab_t *b,int *i_,int *j_,int ib,int jb)
{
  int iStop,jStop;

  if(ib == 0)
    iStop = a->n1;
  else
    iStop = ib;

  if(jb == 0)
    jStop = a->n2;
  else
    jStop = jb;


    for(unsigned int i=0;i < iStop;++i)
    {
      for(unsigned int j=0; j < jStop; ++j)
	{
	 
	  if (getTab(a,i,j) != getTab(b,i,j))
	    {
	      *i_ = i;
	      *j_ = j;
	      // printf(" i : %i, j: %i\n",i,j);    
	      return false;
	    }
	}
     
    }
    return true;
}*/
/****/
