//! @file
//! uncompress kernel in OpenCL
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"
#include "Uncompress.h"
#include "Initialize.h"
#include <math.h>

char * uncompress_program="\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
__kernel void uncompress_coarseKernel_d(uint n1, uint n2, uint n3, uint nseg_c, uint nvctr_c, __global const uint * restrict keyg_c, __global const uint *  restrict keyv_c, __global const double *  restrict psi_c, __global double *  restrict psi_g) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_c ) : ig;\n\
size_t length = nseg_c;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_c;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    first = middle+1;\n\
    length = length - half - 1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while (length > 0);\n\
first--;\n\
size_t iseg = first - keyv_c;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_c[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
//psi_g(i1,1,i2,1,i3,1)=psi_c[ig]\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = psi_c[ig];\n\
};\n\
\n\
__kernel void uncompress_fineKernel_d(uint n1, uint n2, uint n3, uint nseg_f, uint nvctr_f, __global const uint *  restrict keyg_f, __global const uint const *  restrict keyv_f, __global const double *  restrict psi_f, __global double *  restrict psi_g, __local double * tmp) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_f ) : ig;\n\
size_t length = nseg_f;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_f;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    first = middle+1;\n\
    length = length - half - 1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while (length > 0);\n\
first--;\n\
size_t iseg = first - keyv_f;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_f[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
/*\n\
 psig(i1,2,i2,1,i3,1)=psi_f(1,i-i0+jj)\n\
 psig(i1,1,i2,2,i3,1)=psi_f(2,i-i0+jj)\n\
 psig(i1,2,i2,2,i3,1)=psi_f(3,i-i0+jj)\n\
 psig(i1,1,i2,1,i3,2)=psi_f(4,i-i0+jj)\n\
 psig(i1,2,i2,1,i3,2)=psi_f(5,i-i0+jj)\n\
 psig(i1,1,i2,2,i3,2)=psi_f(6,i-i0+jj)\n\
 psig(i1,2,i2,2,i3,2)=psi_f(7,i-i0+jj)*/\n\
size_t i = get_local_id(0);\n\
size_t igr = ig - i;\n\
__local double * tmp_o = tmp + i;\n\
tmp_o[0*64] = psi_f[igr * 7 + i];\n\
tmp_o[1*64] = psi_f[igr * 7 + i + 1*64];\n\
tmp_o[2*64] = psi_f[igr * 7 + i + 2*64];\n\
tmp_o[3*64] = psi_f[igr * 7 + i + 3*64];\n\
tmp_o[4*64] = psi_f[igr * 7 + i + 4*64];\n\
tmp_o[5*64] = psi_f[igr * 7 + i + 5*64];\n\
tmp_o[6*64] = psi_f[igr * 7 + i + 6*64];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
tmp_o = tmp + i*7;\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[0];\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = tmp_o[1];\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[2];\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = tmp_o[3];\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[4];\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = tmp_o[5];\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[6];\n\
};\n\
__kernel void scale_psi_coarseKernel_d(uint nvctr_c, double GPUscal0, __global double * psi_c) {\n\
  size_t ig = get_global_id(0);\n\
  if(ig<nvctr_c)\n\
    psi_c[ig] *= GPUscal0;\n\
}\n\
__kernel void uncompress_scale_coarseKernel_d(uint n1, uint n2, uint n3, uint nseg_c, uint nvctr_c, double GPUscal0, __global const uint *  restrict keyg_c, __global const uint const *  restrict keyv_c, __global const double *  restrict psi_c, __global double *  restrict psi_g) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_c ) : ig;\n\
size_t length = nseg_c;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_c;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    first = middle+1;\n\
    length = length - half - 1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while ( length > 0);\n\
first--;\n\
size_t iseg = first - keyv_c;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_c[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
//psi_g(i1,1,i2,1,i3,1)=psi_c[ig]\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = psi_c[ig]*GPUscal0;\n\
};\n\
__kernel void scale_psi_fineKernel_d(uint nvctr_f, double GPUscal1, double GPUscal2, double GPUscal3, double GPUscal4, double GPUscal5, double GPUscal6, double GPUscal7, __global double *  restrict psi_f) {\n\
  size_t ig = get_global_id(0);\n\
  if(ig<nvctr_f){\n\
    psi_f[ig*7] *= GPUscal1;\n\
    psi_f[ig*7+1] *= GPUscal2;\n\
    psi_f[ig*7+2] *= GPUscal3;\n\
    psi_f[ig*7+3] *= GPUscal4;\n\
    psi_f[ig*7+4] *= GPUscal5;\n\
    psi_f[ig*7+5] *= GPUscal6;\n\
    psi_f[ig*7+6] *= GPUscal7;\n\
  }\n\
}\n\
__kernel void uncompress_scale_fineKernel_d(uint n1, uint n2, uint n3, uint nseg_f, uint nvctr_f, double GPUscal1, double GPUscal2, double GPUscal3, double GPUscal4, double GPUscal5, double GPUscal6, double GPUscal7, __global const uint *  restrict keyg_f, __global const uint const *  restrict keyv_f, __global const double *  restrict psi_f, __global double *  restrict psi_g, __local double * tmp) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_f ) : ig;\n\
size_t length = nseg_f;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_f;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    first = middle+1;\n\
    length = length - half - 1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while (length > 0);\n\
first--;\n\
size_t iseg = first - keyv_f;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_f[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
/*\n\
 psig(i1,2,i2,1,i3,1)=psi_f(1,i-i0+jj)\n\
 psig(i1,1,i2,2,i3,1)=psi_f(2,i-i0+jj)\n\
 psig(i1,2,i2,2,i3,1)=psi_f(3,i-i0+jj)\n\
 psig(i1,1,i2,1,i3,2)=psi_f(4,i-i0+jj)\n\
 psig(i1,2,i2,1,i3,2)=psi_f(5,i-i0+jj)\n\
 psig(i1,1,i2,2,i3,2)=psi_f(6,i-i0+jj)\n\
 psig(i1,2,i2,2,i3,2)=psi_f(7,i-i0+jj)*/\n\
size_t i = get_local_id(0);\n\
size_t igr = ig - i;\n\
__local double * tmp_o = tmp + i;\n\
tmp_o[0*64] = psi_f[igr * 7 + i];\n\
tmp_o[1*64] = psi_f[igr * 7 + i + 1*64];\n\
tmp_o[2*64] = psi_f[igr * 7 + i + 2*64];\n\
tmp_o[3*64] = psi_f[igr * 7 + i + 3*64];\n\
tmp_o[4*64] = psi_f[igr * 7 + i + 4*64];\n\
tmp_o[5*64] = psi_f[igr * 7 + i + 5*64];\n\
tmp_o[6*64] = psi_f[igr * 7 + i + 6*64];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
tmp_o = tmp + i*7;\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[0]*GPUscal1;\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = tmp_o[1]*GPUscal2;\n\
psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[2]*GPUscal3;\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = tmp_o[3]*GPUscal4;\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[4]*GPUscal5;\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] = tmp_o[5]*GPUscal6;\n\
psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] = tmp_o[6]*GPUscal7;\n\
};\n\
";

char * compress_program="\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
__kernel void compress_coarseKernel_d(uint n1, uint n2, uint n3, uint nseg_c, uint nvctr_c, __global const uint *  restrict keyg_c, __global const uint * restrict keyv_c, __global double *  restrict psi_c, __global const double *  restrict psi_g) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_c ) : ig;\n\
size_t length = nseg_c;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_c;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    first = middle+1;\n\
    length = length - half - 1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while ( length > 0);\n\
first--;\n\
size_t iseg = first - keyv_c;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_c[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
//psi_g(i1,1,i2,1,i3,1)=psi_c[ig]\n\
psi_c[ig] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1];\n\
};\n\
__kernel void compress_scale_coarseKernel_d(uint n1, uint n2, uint n3, uint nseg_c, uint nvctr_c, double GPUscal0, __global const uint *  restrict keyg_c, __global const uint *  restrict keyv_c, __global double *  restrict psi_c, __global const double *  restrict psi_g) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_c ) : ig;\n\
size_t length = nseg_c;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_c;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    first = middle+1;\n\
    length = length - half - 1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while ( length > 0);\n\
first--;\n\
size_t iseg = first - keyv_c;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_c[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
//psi_g(i1,1,i2,1,i3,1)=psi_c[ig]\n\
psi_c[ig] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1]*GPUscal0;\n\
};\n\
\n\
//hypothesis : nseg_f > 0\n\
__kernel void compress_fineKernel_d(uint n1, uint n2, uint n3, uint nseg_f, uint nvctr_f, __global const uint *  restrict keyg_f, __global const uint *  restrict keyv_f, __global double *  restrict psi_f, __global const double *  restrict psi_g, __local double * tmp) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_f ) : ig;\n\
size_t length = nseg_f;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_f;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    length = length - half - 1;\n\
    first = middle+1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while ( length > 0);\n\
first--;\n\
size_t iseg = first - keyv_f;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_f[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
size_t i = get_local_id(0);\n\
__local double * tmp_o = tmp + i*7;\n\
tmp_o[0] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1];\n\
tmp_o[1] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1];\n\
tmp_o[2] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1];\n\
tmp_o[3] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1];\n\
tmp_o[4] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1];\n\
tmp_o[5] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1];\n\
tmp_o[6] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
size_t igr = ig - i;\n\
tmp_o = tmp + i;\n\
psi_f[igr * 7 + i] = tmp_o[0*64];\n\
psi_f[igr * 7 + i + 1*64] = tmp_o[1*64];\n\
psi_f[igr * 7 + i + 2*64] = tmp_o[2*64];\n\
psi_f[igr * 7 + i + 3*64] = tmp_o[3*64];\n\
psi_f[igr * 7 + i + 4*64] = tmp_o[4*64];\n\
psi_f[igr * 7 + i + 5*64] = tmp_o[5*64];\n\
psi_f[igr * 7 + i + 6*64] = tmp_o[6*64];\n\
};\n\
__kernel void compress_scale_fineKernel_d(uint n1, uint n2, uint n3, uint nseg_f, uint nvctr_f, double GPUscal1, double GPUscal2, double GPUscal3, double GPUscal4, double GPUscal5, double GPUscal6, double GPUscal7, __global const uint *  restrict keyg_f, __global const uint *  restrict keyv_f, __global double *  restrict psi_f, __global const double *  restrict psi_g, __local double * tmp) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - nvctr_f ) : ig;\n\
size_t length = nseg_f;\n\
__global const uint * first;\n\
__global const uint * middle;\n\
size_t half;\n\
first = keyv_f;\n\
do {\n\
  half = length / 2;\n\
  middle = first + half;\n\
  if( *middle-1 <= ig ) {\n\
    length = length - half - 1;\n\
    first = middle+1;\n\
  } else {\n\
    length = half;\n\
  }\n\
} while ( length > 0);\n\
first--;\n\
size_t iseg = first - keyv_f;\n\
size_t j0,ii,i1,i2,i3;\n\
j0=keyg_f[iseg*2];\n\
ii=j0-1;\n\
i3=ii/((n1)*(n2));\n\
ii=ii-i3*(n1)*(n2);\n\
i2=ii/(n1);\n\
i1=ii-i2*(n1)+ig-(*first-1);\n\
size_t i = get_local_id(0);\n\
__local double * tmp_o = tmp + i*7;\n\
tmp_o[0] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] * GPUscal1;\n\
tmp_o[1] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] * GPUscal2;\n\
tmp_o[2] = psi_g[ ( ( ( (0 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] * GPUscal3;\n\
tmp_o[3] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] * GPUscal4;\n\
tmp_o[4] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 0 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] * GPUscal5;\n\
tmp_o[5] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 0 ) * n1  + i1] * GPUscal6;\n\
tmp_o[6] = psi_g[ ( ( ( (1 * n3 + i3 ) * 2 + 1 ) * n2 + i2 ) * 2 + 1 ) * n1  + i1] * GPUscal7;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
size_t igr = ig - i;\n\
tmp_o = tmp + i;\n\
psi_f[igr * 7 + i] = tmp_o[0*64];\n\
psi_f[igr * 7 + i + 1*64] = tmp_o[1*64];\n\
psi_f[igr * 7 + i + 2*64] = tmp_o[2*64];\n\
psi_f[igr * 7 + i + 3*64] = tmp_o[3*64];\n\
psi_f[igr * 7 + i + 4*64] = tmp_o[4*64];\n\
psi_f[igr * 7 + i + 5*64] = tmp_o[5*64];\n\
psi_f[igr * 7 + i + 6*64] = tmp_o[6*64];\n\
};\n\
";

#define A2 3.55369228991319019
#define B2 24.8758460293923314
inline void scale_psi_coarse_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *nvctr_c, double *h, double *c, cl_mem *psi_c) {
  double hh1 = 0.125/(h[0]*h[0]);
  double hh2 = 0.125/(h[1]*h[1]);
  double hh3 = 0.125/(h[2]*h[2]);
  double GPUscal0 = 1/sqrt(A2*hh1 + A2*hh2 + A2*hh3 + *c);
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i = 0;
  clSetKernelArg(kernel, i++, sizeof(*nvctr_c), (void*)nvctr_c);
  clSetKernelArg(kernel, i++, sizeof(GPUscal0), (void*)&GPUscal0);
  clSetKernelArg(kernel, i++, sizeof(*psi_c), (void*)psi_c);
  size_t localWorkSize[]= { block_size_i };
  size_t globalWorkSize[] = { shrRoundUp(block_size_i, *nvctr_c) } ;
  ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue scale_psi_coarse kernel!");
}

inline void un_compress_scale_coarse_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions, double *h, double *c,
                                    cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                    cl_mem *psi_c, cl_mem * psi) {
  double hh1 = 0.125/(h[0]*h[0]);
  double hh2 = 0.125/(h[1]*h[1]);
  double hh3 = 0.125/(h[2]*h[2]);
  double GPUscal0 = 1.0/sqrt(A2*hh1 + A2*hh2 + A2*hh3 + *c);
    cl_int ciErrNum;
    size_t block_size_i=64;
    assert( *nvctr_c >= block_size_i );
    cl_uint i = 0;
    clSetKernelArg(kernel, i++, sizeof(dimensions[0]), (void*)&(dimensions[0]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[1]), (void*)&(dimensions[1]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[2]), (void*)&(dimensions[2]));
    clSetKernelArg(kernel, i++, sizeof(*nseg_c), (void*)nseg_c);
    clSetKernelArg(kernel, i++, sizeof(*nvctr_c), (void*)nvctr_c);
    clSetKernelArg(kernel, i++, sizeof(GPUscal0), (void*)&GPUscal0);
    clSetKernelArg(kernel, i++, sizeof(*keyg_c), (void*)keyg_c);
    clSetKernelArg(kernel, i++, sizeof(*keyv_c), (void*)keyv_c);
    clSetKernelArg(kernel, i++, sizeof(*psi_c), (void*)psi_c);
    clSetKernelArg(kernel, i++, sizeof(*psi), (void*)psi);
    size_t localWorkSize[]= { block_size_i };
    size_t globalWorkSize[] = { shrRoundUp(block_size_i, *nvctr_c) } ;
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue uncompress_scale_coarse kernel!");

}

inline void scale_psi_fine_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *nvctr_f, double *h, double *c, cl_mem *psi_f) {
  double hh1 = 0.125/(h[0]*h[0]);
  double hh2 = 0.125/(h[1]*h[1]);
  double hh3 = 0.125/(h[2]*h[2]);
  double GPUscal1 = 1.0/sqrt(B2*hh1 + A2*hh2 + A2*hh3 + *c);
  double GPUscal2 = 1.0/sqrt(A2*hh1 + B2*hh2 + A2*hh3 + *c);
  double GPUscal3 = 1.0/sqrt(B2*hh1 + B2*hh2 + A2*hh3 + *c);
  double GPUscal4 = 1.0/sqrt(A2*hh1 + A2*hh2 + B2*hh3 + *c);
  double GPUscal5 = 1.0/sqrt(B2*hh1 + A2*hh2 + B2*hh3 + *c);
  double GPUscal6 = 1.0/sqrt(A2*hh1 + B2*hh2 + B2*hh3 + *c);
  double GPUscal7 = 1.0/sqrt(B2*hh1 + B2*hh2 + B2*hh3 + *c);
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i = 0;
  clSetKernelArg(kernel, i++, sizeof(*nvctr_f), (void*)nvctr_f);
  clSetKernelArg(kernel, i++, sizeof(GPUscal1), (void*)&GPUscal1);
  clSetKernelArg(kernel, i++, sizeof(GPUscal2), (void*)&GPUscal2);
  clSetKernelArg(kernel, i++, sizeof(GPUscal3), (void*)&GPUscal3);
  clSetKernelArg(kernel, i++, sizeof(GPUscal4), (void*)&GPUscal4);
  clSetKernelArg(kernel, i++, sizeof(GPUscal5), (void*)&GPUscal5);
  clSetKernelArg(kernel, i++, sizeof(GPUscal6), (void*)&GPUscal6);
  clSetKernelArg(kernel, i++, sizeof(GPUscal7), (void*)&GPUscal7);
  clSetKernelArg(kernel, i++, sizeof(*psi_f), (void*)psi_f);
  size_t localWorkSize[]= { block_size_i };
  size_t globalWorkSize[] = { shrRoundUp(block_size_i, *nvctr_f) } ;
  ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue scale_psi_coarse kernel!");
}

inline void un_compress_scale_fine_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions, double *h, double *c,
                                    cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                    cl_mem *psi_f, cl_mem * psi) {
  double hh1 = 0.125/(h[0]*h[0]);
  double hh2 = 0.125/(h[1]*h[1]);
  double hh3 = 0.125/(h[2]*h[2]);
  double GPUscal1 = 1.0/sqrt(B2*hh1 + A2*hh2 + A2*hh3 + *c);
  double GPUscal2 = 1.0/sqrt(A2*hh1 + B2*hh2 + A2*hh3 + *c);
  double GPUscal3 = 1.0/sqrt(B2*hh1 + B2*hh2 + A2*hh3 + *c);
  double GPUscal4 = 1.0/sqrt(A2*hh1 + A2*hh2 + B2*hh3 + *c);
  double GPUscal5 = 1.0/sqrt(B2*hh1 + A2*hh2 + B2*hh3 + *c);
  double GPUscal6 = 1.0/sqrt(A2*hh1 + B2*hh2 + B2*hh3 + *c);
  double GPUscal7 = 1.0/sqrt(B2*hh1 + B2*hh2 + B2*hh3 + *c);
    cl_int ciErrNum;
    size_t block_size_i=64;
    assert( *nvctr_f >= block_size_i );
    cl_uint i = 0;
    clSetKernelArg(kernel, i++, sizeof(dimensions[0]), (void*)&(dimensions[0]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[1]), (void*)&(dimensions[1]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[2]), (void*)&(dimensions[2]));
    clSetKernelArg(kernel, i++, sizeof(*nseg_f), (void*)nseg_f);
    clSetKernelArg(kernel, i++, sizeof(*nvctr_f), (void*)nvctr_f);
    clSetKernelArg(kernel, i++, sizeof(GPUscal1), (void*)&GPUscal1);
    clSetKernelArg(kernel, i++, sizeof(GPUscal2), (void*)&GPUscal2);
    clSetKernelArg(kernel, i++, sizeof(GPUscal3), (void*)&GPUscal3);
    clSetKernelArg(kernel, i++, sizeof(GPUscal4), (void*)&GPUscal4);
    clSetKernelArg(kernel, i++, sizeof(GPUscal5), (void*)&GPUscal5);
    clSetKernelArg(kernel, i++, sizeof(GPUscal6), (void*)&GPUscal6);
    clSetKernelArg(kernel, i++, sizeof(GPUscal7), (void*)&GPUscal7);
    clSetKernelArg(kernel, i++, sizeof(*keyg_f), (void*)keyg_f);
    clSetKernelArg(kernel, i++, sizeof(*keyv_f), (void*)keyv_f);
    clSetKernelArg(kernel, i++, sizeof(*psi_f), (void*)psi_f);
    clSetKernelArg(kernel, i++, sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++, sizeof(double)*block_size_i*7, NULL);
    size_t localWorkSize[]= { block_size_i };
    size_t globalWorkSize[] = { shrRoundUp(block_size_i, *nvctr_f) } ;
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue uncompress_scale_fine kernel!");

}

inline void un_compress_coarse_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions,
                                    cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                    cl_mem *psi_c, cl_mem * psi) {
    cl_int ciErrNum;
    size_t block_size_i=64;
    assert( *nvctr_c >= block_size_i );
    cl_uint i = 0;
    clSetKernelArg(kernel, i++, sizeof(dimensions[0]), (void*)&(dimensions[0]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[1]), (void*)&(dimensions[1]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[2]), (void*)&(dimensions[2]));
    clSetKernelArg(kernel, i++, sizeof(*nseg_c), (void*)nseg_c);
    clSetKernelArg(kernel, i++, sizeof(*nvctr_c), (void*)nvctr_c);
    clSetKernelArg(kernel, i++, sizeof(*keyg_c), (void*)keyg_c);
    clSetKernelArg(kernel, i++, sizeof(*keyv_c), (void*)keyv_c);
    clSetKernelArg(kernel, i++, sizeof(*psi_c), (void*)psi_c);
    clSetKernelArg(kernel, i++, sizeof(*psi), (void*)psi);
    size_t localWorkSize[]= { block_size_i };
    size_t globalWorkSize[] = { shrRoundUp(block_size_i, *nvctr_c) } ;
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue compress_coarse kernel!");

}

inline void un_compress_fine_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions,
                                    cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                    cl_mem *psi_f, cl_mem * psi) {
    cl_int ciErrNum;
    size_t block_size_i=64;
    assert( *nvctr_f >= block_size_i );
    cl_uint i = 0;
    clSetKernelArg(kernel, i++, sizeof(dimensions[0]), (void*)&(dimensions[0]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[1]), (void*)&(dimensions[1]));
    clSetKernelArg(kernel, i++, sizeof(dimensions[2]), (void*)&(dimensions[2]));
    clSetKernelArg(kernel, i++, sizeof(*nseg_f), (void*)nseg_f);
    clSetKernelArg(kernel, i++, sizeof(*nvctr_f), (void*)nvctr_f);
    clSetKernelArg(kernel, i++, sizeof(*keyg_f), (void*)keyg_f);
    clSetKernelArg(kernel, i++, sizeof(*keyv_f), (void*)keyv_f);
    clSetKernelArg(kernel, i++, sizeof(*psi_f), (void*)psi_f);
    clSetKernelArg(kernel, i++, sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++, sizeof(double)*block_size_i*7, NULL);
    size_t localWorkSize[]= { block_size_i };
    size_t globalWorkSize[] = { shrRoundUp(block_size_i, *nvctr_f) } ;
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue compress_fine kernel!");

}

cl_program uncompressProgram;
cl_program compressProgram;

void create_uncompress_kernels(struct bigdft_kernels * kernels){
    cl_int ciErrNum = CL_SUCCESS;
    kernels->uncompress_coarse_kernel_d=clCreateKernel(uncompressProgram,"uncompress_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 1!");
    kernels->uncompress_fine_kernel_d=clCreateKernel(uncompressProgram,"uncompress_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 2!");
    kernels->uncompress_scale_coarse_kernel_d=clCreateKernel(uncompressProgram,"uncompress_scale_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 3!");
    kernels->uncompress_scale_fine_kernel_d=clCreateKernel(uncompressProgram,"uncompress_scale_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 4!");
    kernels->scale_psi_fine_kernel_d=clCreateKernel(uncompressProgram,"scale_psi_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 5!");
    kernels->scale_psi_coarse_kernel_d=clCreateKernel(uncompressProgram,"scale_psi_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 6!");
    kernels->compress_coarse_kernel_d=clCreateKernel(compressProgram,"compress_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 7!");
    kernels->compress_fine_kernel_d=clCreateKernel(compressProgram,"compress_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 8!");
    kernels->compress_scale_coarse_kernel_d=clCreateKernel(compressProgram,"compress_scale_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 9!");
    kernels->compress_scale_fine_kernel_d=clCreateKernel(compressProgram,"compress_scale_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel 10!");
}

void build_uncompress_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    uncompressProgram = clCreateProgramWithSource(*context,1,(const char**) &uncompress_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(uncompressProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build uncompress program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(uncompressProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }

    compressProgram = clCreateProgramWithSource(*context,1,(const char**) &compress_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(compressProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build compress program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(compressProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void FC_FUNC_(scale_psi_d,SCALE_PSI_D)(bigdft_command_queue *command_queue, cl_uint *nvctr_c, cl_uint *nvctr_f, double *h, double *c, cl_mem *psi_c,  cl_mem *psi_f) {
  scale_psi_coarse_generic((*command_queue)->kernels.scale_psi_coarse_kernel_d, (*command_queue)->command_queue, nvctr_c, h, c, psi_c);
  if(nvctr_f == 0) return;
  scale_psi_fine_generic((*command_queue)->kernels.scale_psi_fine_kernel_d, (*command_queue)->command_queue, nvctr_f, h, c, psi_f);
}

void FC_FUNC_(uncompress_d,UNCOMPRESS_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
				       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out) {
  cl_uint full_size = dimensions[0] * dimensions[1] * dimensions[2] * 8;
  double init = 0.0;
  v_initialize_generic((*command_queue)->kernels.v_initialize_kernel_d, (*command_queue)->command_queue, &full_size, psi_out, &init);
  un_compress_coarse_generic((*command_queue)->kernels.uncompress_coarse_kernel_d, (*command_queue)->command_queue, dimensions, nseg_c, nvctr_c, keyg_c, keyv_c, psi_c, psi_out);
  if(nvctr_f == 0) return;
  un_compress_fine_generic((*command_queue)->kernels.uncompress_fine_kernel_d, (*command_queue)->command_queue, dimensions, nseg_f, nvctr_f, keyg_f, keyv_f, psi_f, psi_out);
}

void FC_FUNC_(compress_d,COMPRESS_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi) {
  un_compress_coarse_generic((*command_queue)->kernels.compress_coarse_kernel_d, (*command_queue)->command_queue, dimensions, nseg_c, nvctr_c, keyg_c, keyv_c, psi_c, psi);
  if(nvctr_f == 0) return;
  un_compress_fine_generic((*command_queue)->kernels.compress_fine_kernel_d, (*command_queue)->command_queue, dimensions, nseg_f, nvctr_f, keyg_f, keyv_f, psi_f, psi);
}

void FC_FUNC_(uncompress_scale_d,UNCOMPRESS_SCALE_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, double *h, double *c,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out) {
  cl_uint full_size = dimensions[0] * dimensions[1] * dimensions[2] * 8;
  double init = 0.0;
  v_initialize_generic((*command_queue)->kernels.v_initialize_kernel_d, (*command_queue)->command_queue, &full_size, psi_out, &init);
  un_compress_scale_coarse_generic((*command_queue)->kernels.uncompress_scale_coarse_kernel_d, (*command_queue)->command_queue, dimensions, h, c, nseg_c, nvctr_c, keyg_c, keyv_c, psi_c, psi_out);
  if(nvctr_f == 0) return;
  un_compress_scale_fine_generic((*command_queue)->kernels.uncompress_scale_fine_kernel_d, (*command_queue)->command_queue, dimensions, h, c, nseg_f, nvctr_f, keyg_f, keyv_f, psi_f, psi_out);
}
void FC_FUNC_(compress_scale_d,COMPRESS_SCALE_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, double *h, double *c,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi) {
  un_compress_scale_coarse_generic((*command_queue)->kernels.compress_scale_coarse_kernel_d, (*command_queue)->command_queue, dimensions, h, c, nseg_c, nvctr_c, keyg_c, keyv_c, psi_c, psi);
  if(nvctr_f == 0) return;
  un_compress_scale_fine_generic((*command_queue)->kernels.compress_scale_fine_kernel_d, (*command_queue)->command_queue, dimensions, h, c, nseg_f, nvctr_f, keyg_f, keyv_f, psi_f, psi);
}
void clean_uncompress_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->uncompress_coarse_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->uncompress_fine_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->uncompress_scale_coarse_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->uncompress_scale_fine_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->compress_coarse_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->compress_fine_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->compress_scale_coarse_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->compress_scale_fine_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->scale_psi_fine_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->scale_psi_coarse_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}
void clean_uncompress_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(uncompressProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram(compressProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
