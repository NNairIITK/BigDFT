/*
** cudafct.h
** 
** Made by Matthieu Ospici
** Login   <mo219174@badiane>
** 
** Started on  Wed Apr 16 16:51:28 2008 Matthieu Ospici
** Last update Wed Apr 16 16:51:28 2008 Matthieu Ospici
*/

#ifndef   	CUDAFCT_H_
# define   	CUDAFCT_H_

int c_cudaGetDeviceCount(int *s_gpuCount);
int c_cudaSetDevice(int device);


int c_cudaMalloc(void**,size_t memsize);

int c_cudaMallocHost(void**,size_t memsize);
int c_cudaMallocHost(float**,size_t memsize);

int c_cuda_gpu_send_pi(void *dest, const void *src,  size_t memByte);

int c_cuda_gpu_recv_pi(void *dest, const void *src,  size_t memByte);

int c_cuda_setdevice(int device);
int c_cuda_get_device(int *dev);

int c_cuda_setdevice_ctx(int device);
int c_cuda_ctxpopcur(void *contex);
int c_cuda_ctxpushcur(void *contex);

#endif 	    /* !CUDAFCT_H_ */
