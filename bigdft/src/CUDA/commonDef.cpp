//uniformise the tabular of the number of elements treated by each
//thread (counterpart of uniformiseTab)
void correctSequence(int thds,int elem,int * tab)
{
  int j,i;
  //put to zero all the values;
  for(j=0;j< elem; ++j)
    {
      tab[j]=0;
    }

  //then start to fill consecutively until reaching of the end
  //if elem > thds no element will be zero
  //this is the most balanced choice
  for(i=0;i< elem; ++i)
    {
      tab[i % thds]+=1; 
    }
}
