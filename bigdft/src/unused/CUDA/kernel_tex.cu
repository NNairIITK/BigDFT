__global__ void conv1d_stride_tex_m(int n,int ndat, float *psi_in, float *psi_out)
{

  //line treated by the given block
  unsigned int lineOffset = min(blockIdx.y*NUM_LINES,ndat-NUM_LINES);
  //starting element treated by the block
  unsigned int elemOffset = min(blockIdx.x*par.ElementsPerBlock,n-par.ElementsPerBlock);

  //line treated by the given block
  //unsigned int lineOffset = min(blockIdx.y*par.LinesPerBlock,ndat-par.LinesPerBlock);
  //starting element treated by the block
  //unsigned int elemOffset = min(blockIdx.x*par.ElementsPerBlock,n-par.ElementsPerBlock);

  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  //shared memory array
  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread in ndat axis
  //which is the input base element
  unsigned int BaseElem = par.thline[tid_hw] + lineOffset;
  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = par.thelem[tid_hw] + par.hwoffset_copy[hwid];

  unsigned int ShBaseElem = tid_hw + NUM_LINES*par.hwoffset_copy[hwid];
  //unsigned int ShBaseElem = tid_hw + par.LinesPerBlock*par.hwoffset_copy[hwid];

  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-par.lowfil+thelem;i < par.hwelem_copy[hwid] ; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;

      float x=((float) (BaseElem) + 0.5f )/((float) (ndat));
      float y=((float) (npos) + 0.5f )/((float) (n));


      psi_sh[ShBaseElem]= tex2D(psi_tex, x,y); // -+-+ psi_in[BaseElem+ndat*npos];
      //psi_sh[ShBaseElem]=tex1Dfetch(psi_tex,BaseElem+ndat*npos);

      ShBaseElem += HALF_WARP_SIZE;
      ipos += HW_ELEM;
      //ipos += par.ElementsPerHalfWarp;
      
    }

  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = par.thelem[tid_hw] + par.hwoffset_calc[hwid];
  //base element for the given thread in shared memory
  ShBaseElem = tid_hw + NUM_LINES*par.hwoffset_calc[hwid];
  //ShBaseElem = tid_hw + par.LinesPerBlock*par.hwoffset_calc[hwid];

  //output base element, from the input one
  BaseElem =  n*BaseElem+ thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)

  for(int i=0;i < par.hwelem_calc[hwid]; ++i)
    {
      //values of the convolution
      register float conv = 
	//hand-unrolled loop (16 elements for this filter)
	//order changed for increasing the precision

	MFIL0 *psi_sh[ShBaseElem               ] +
	MFIL15*psi_sh[ShBaseElem + 15*NUM_LINES] +
	MFIL1 *psi_sh[ShBaseElem +   NUM_LINES ] +
	MFIL14*psi_sh[ShBaseElem + 14*NUM_LINES] +
	MFIL2 *psi_sh[ShBaseElem + 2*NUM_LINES ] +
	MFIL13*psi_sh[ShBaseElem + 13*NUM_LINES] +
	MFIL3 *psi_sh[ShBaseElem + 3*NUM_LINES ] +
	MFIL12*psi_sh[ShBaseElem + 12*NUM_LINES] +
	MFIL4 *psi_sh[ShBaseElem + 4*NUM_LINES ] +
	MFIL11*psi_sh[ShBaseElem + 11*NUM_LINES] +
	MFIL5 *psi_sh[ShBaseElem + 5*NUM_LINES ] +
	MFIL10*psi_sh[ShBaseElem + 10*NUM_LINES] +
	MFIL6 *psi_sh[ShBaseElem + 6*NUM_LINES ] +
	MFIL9 *psi_sh[ShBaseElem + 9*NUM_LINES ] +
	MFIL7 *psi_sh[ShBaseElem + 7*NUM_LINES ] +
	MFIL8 *psi_sh[ShBaseElem + 8*NUM_LINES ] ;

	/*
	par.fil[0]*psi_sh[ShBaseElem               ] +
	par.fil[15]*psi_sh[ShBaseElem + 15*NUM_LINES] +
	par.fil[1]*psi_sh[ShBaseElem +   NUM_LINES ] +
	par.fil[14]*psi_sh[ShBaseElem + 14*NUM_LINES] +
	par.fil[2]*psi_sh[ShBaseElem + 2*NUM_LINES ] +
	par.fil[13]*psi_sh[ShBaseElem + 13*NUM_LINES] +
	par.fil[3]*psi_sh[ShBaseElem + 3*NUM_LINES ] +
	par.fil[12]*psi_sh[ShBaseElem + 12*NUM_LINES] +
	par.fil[4]*psi_sh[ShBaseElem + 4*NUM_LINES ] +
	par.fil[11]*psi_sh[ShBaseElem + 11*NUM_LINES] +
	par.fil[5]*psi_sh[ShBaseElem + 5*NUM_LINES ] +
	par.fil[10]*psi_sh[ShBaseElem + 10*NUM_LINES] +
	par.fil[6]*psi_sh[ShBaseElem + 6*NUM_LINES ] +
	par.fil[9]*psi_sh[ShBaseElem + 9*NUM_LINES ] +
	par.fil[7]*psi_sh[ShBaseElem + 7*NUM_LINES ] +
	par.fil[8]*psi_sh[ShBaseElem + 8*NUM_LINES ] ;

      
	par.fil[0]*psi_sh[ShBaseElem               ] +
	par.fil[15]*psi_sh[ShBaseElem + 15*par.LinesPerBlock] +
	par.fil[1]*psi_sh[ShBaseElem +   par.LinesPerBlock ] +
	par.fil[14]*psi_sh[ShBaseElem + 14*par.LinesPerBlock] +
	par.fil[2]*psi_sh[ShBaseElem + 2*par.LinesPerBlock ] +
	par.fil[13]*psi_sh[ShBaseElem + 13*par.LinesPerBlock] +
	par.fil[3]*psi_sh[ShBaseElem + 3*par.LinesPerBlock ] +
	par.fil[12]*psi_sh[ShBaseElem + 12*par.LinesPerBlock] +
	par.fil[4]*psi_sh[ShBaseElem + 4*par.LinesPerBlock ] +
	par.fil[11]*psi_sh[ShBaseElem + 11*par.LinesPerBlock] +
	par.fil[5]*psi_sh[ShBaseElem + 5*par.LinesPerBlock ] +
	par.fil[10]*psi_sh[ShBaseElem + 10*par.LinesPerBlock] +
	par.fil[6]*psi_sh[ShBaseElem + 6*par.LinesPerBlock ] +
	par.fil[9]*psi_sh[ShBaseElem + 9*par.LinesPerBlock ] +
	par.fil[7]*psi_sh[ShBaseElem + 7*par.LinesPerBlock ] +
	par.fil[8]*psi_sh[ShBaseElem + 8*par.LinesPerBlock ] ;


	par.fil[0]*psi_sh[ShBaseElem               ] +
	par.fil[1]*psi_sh[ShBaseElem +   par.LinesPerBlock ] +
	par.fil[2]*psi_sh[ShBaseElem + 2*par.LinesPerBlock ] +
	par.fil[3]*psi_sh[ShBaseElem + 3*par.LinesPerBlock ] +
	par.fil[4]*psi_sh[ShBaseElem + 4*par.LinesPerBlock ] +
	par.fil[5]*psi_sh[ShBaseElem + 5*par.LinesPerBlock ] +
	par.fil[6]*psi_sh[ShBaseElem + 6*par.LinesPerBlock ] +
	par.fil[7]*psi_sh[ShBaseElem + 7*par.LinesPerBlock ] +
	par.fil[8]*psi_sh[ShBaseElem + 8*par.LinesPerBlock ] +
	par.fil[9]*psi_sh[ShBaseElem + 9*par.LinesPerBlock ] +
	par.fil[10]*psi_sh[ShBaseElem + 10*par.LinesPerBlock] +
	par.fil[11]*psi_sh[ShBaseElem + 11*par.LinesPerBlock] +
	par.fil[12]*psi_sh[ShBaseElem + 12*par.LinesPerBlock] +
	par.fil[13]*psi_sh[ShBaseElem + 13*par.LinesPerBlock] +
	par.fil[14]*psi_sh[ShBaseElem + 14*par.LinesPerBlock] +
	par.fil[15]*psi_sh[ShBaseElem + 15*par.LinesPerBlock];
      */

      psi_out[BaseElem]=conv;
      //psi_sh[ShBaseElem+par.lowfil*par.LinesPerBlock]; //for testing only

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      //BaseElem += par.ElementsPerHalfWarp;

      
    }

 
}
