!!****p* OpenCL/conv_check
!! FUNCTION
!!    Program test for the convolution in GPU
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! CREATION DATE
!!    Septembre 2008
!!
!! SOURCE
!!

program conv_check_fft
  use module_base
  use Poisson_Solver
  implicit none
  integer  :: n1,n2,n3 
  real(gp) :: hx,r2,sigma2,x,y,z,maxdiff,arg
  real(wp), dimension(:,:,:,:), allocatable :: psi_in,psi_out 
  !local variables
  character(len=*), parameter :: subname='conv_check_fft'
  character(len=50) :: chain
  integer :: i,i_stat,j,i1,i2,i3,ntimes,i_all
  integer :: l,ierror
  real(wp) :: tt,scale
  real(gp) :: CPUtime,GPUtime,ekin,ehartree
  integer, dimension(:), allocatable :: modarr
  real(kind=8), dimension(:), allocatable :: psi 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda,v_cuda,v_cuda_str,v_cuda_str1
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_str 
  real(kind=4), dimension(:,:,:,:), allocatable :: psi_cuda_l,v_cuda_l  
  real(kind=8) :: ekinGPUd,ehartreeGPUd
  real(kind=8) :: psi_GPU,work_GPU,work2_GPU,k_GPU
  real(kind=8) :: context,queue
  integer, parameter :: lowfilK=-14,lupfilK=14 
  real(kind=8), dimension(lowfilK:lupfilK) :: fil
  integer(kind=8) :: tsc0, tsc1
  !objects for the 3D Poisson solver
  real(kind=8), dimension(:), allocatable :: rhopot,rhopot2,rhopot1
  real(kind=8), dimension(:), pointer :: pkernel,pkernel2
  integer :: plan,plan1,plan1_,plan2,plan3,plan3_,size1,size2,sizek,switch_alg,geo1,geo2,geo3
  real(kind=8) :: scal
  integer :: count_0

  write(chain,'(a6)') 'fort.1'
  open(unit=1,file=trim(chain),status='old')

  read(unit=1,fmt=*,iostat=ierror) n1,n2,n3,ntimes
  if (ierror /= 0) then
     write(*,*) "In a file 'fort.1', put a line with:"
     write(*,*) "n1 n2 n3 ntimes"
     write(*,*) "where:"
     write(*,*) "- n1 n2 n3 are the dimensions of the real space"
     write(*,*) "- ntimes is the number of convolutions"
     stop
  end if

 !n1=150
 !n2=150
 !n3=150

 !ntimes=1

  hx=0.1e0_gp
  !n(c) hy=0.1e0_gp
  !n(c) hz=0.1e0_gp

  
   !allocate arrays
   allocate(psi_in(2,n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in,'psi_in',subname)
   allocate(psi_out(2,n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_out,'psi_out',subname)

   allocate(psi_cuda(2,n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda,'psi_cuda',subname)
   allocate(psi_cuda_str(2,n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_str,'psi_cuda_str',subname)
   allocate(v_cuda(2,n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda,'v_cuda',subname)
   allocate(v_cuda_str(2,n1,n2*n3,2+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_str,'v_cuda_str',subname)
   allocate(v_cuda_str1(2,n2*n3,n1,2+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_str1,'v_cuda_str1',subname)
   allocate(psi_cuda_l(2,n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_l,'psi_cuda_l',subname)
   allocate(v_cuda_l(2,n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_l,'v_cuda_l',subname)
   allocate(rhopot(n1*n2*n3+ndebug),stat=i_stat)
   call memocc(i_stat,rhopot,'rhopot',subname)
   allocate(rhopot1(n1*n2*n3+ndebug),stat=i_stat)
   call memocc(i_stat,rhopot1,'rhopot1',subname)
   allocate(rhopot2(n1*n2*n3+ndebug),stat=i_stat)
   call memocc(i_stat,rhopot2,'rhopot2',subname)

   !initialise array
   sigma2=0.25d0*((n1*hx)**2)
   do i=1,n2*n3
      do i1=1,n1
        do i2=1,2
          x=hx*real(i1-n1/2-1,kind=8)
          !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
          r2=x**2
          arg=0.5d0*r2/sigma2
          tt=dexp(-arg)
          call random_number(tt)
          psi_in(i2,i1,i,1)=tt
        end do
      end do
   end do

   !initialize rhopots
   call vcopy(n1*n2*n3,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3,psi_in(1,1,1,1),2,rhopot2(1),1)

   !the input and output arrays must be reverted in this implementation
   do i=1,n2*n3
      do i1=1,n1
        do i2=1,2
          v_cuda_str1(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_str1(i2,i,i1,2)=0.0
          v_cuda(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_l(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
        end do
      end do
   end do


  write(*,'(a,i6,i6)')'CPU FFT, dimensions:',n1,n2*n3

   call nanosec(tsc0);
   i3=1
   call fft_1d_ctoc(-1,n2*n3,n1,v_cuda_str1,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   write(*,'(a,i6,i6)')'GPU FFT, dimensions:',n1,n2*n3

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   call cuda_1d_plan(n1,n2*n3,plan)
   call cudamalloc(size2,work_GPU)
   call cudamalloc(size2,psi_GPU)
   call reset_gpu_data(size2,psi_in,work_GPU)   

   call nanosec(tsc0);
   do i=1,ntimes
    call cuda_1d_forward(plan,work_GPU,psi_GPU)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size2,psi_cuda,psi_GPU)
   call cudafree(work_GPU)
   call cudafree(psi_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_2D_cplx_results_t( n1, n2*n3,  psi_cuda, v_cuda_str1(1,1,1,i3), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   do i=1,n2*n3
      do i1=1,n1
        do i2=1,2
          v_cuda_str(i2,i1,i,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_str(i2,i1,i,2)=0.0
          v_cuda(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_l(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
        end do
      end do
   end do


  write(*,'(a,i6,i6,i6)')'CPU 3D FFT, dimensions:',n1,n2,n3

   call nanosec(tsc0);
   i3=1
   call FFT(n1,n2,n3,n1,n2,n3,v_cuda_str,-1,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)


   write(*,'(a,i6,i6,i6)')'GPU 3D FFT, dimensions:',n1,n2,n3

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   call cuda_3d_plan(n1,n2,n3,plan)
   call cudamalloc(size2,work_GPU)
   call cudamalloc(size2,psi_GPU)
   call reset_gpu_data(size2,psi_in,work_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_forward(plan,work_GPU,psi_GPU)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size2,psi_cuda,psi_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   call compare_3D_cplx_results(n1, n2, n3, v_cuda_str(1,1,1,i3), psi_cuda, maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * (log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)
 
  write(*,'(a,i6,i6,i6)')'CPU 3D Reverse FFT, dimensions:',n1,n2,n3

   call nanosec(tsc0);
   i3=1
   call FFT(n1,n2,n3,n1,n2,n3,v_cuda_str,1,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)

   write(*,'(a,i6,i6,i6)')'GPU 3D Reverse FFT, dimensions:',n1,n2,n3

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   call cuda_3d_plan(n1,n2,n3,plan)
   call cudamalloc(size2,work_GPU)
   call cudamalloc(size2,psi_GPU)
   call reset_gpu_data(size2,psi_cuda,work_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_inverse(n1,n2,n3,plan,work_GPU,psi_GPU)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size2,psi_cuda,psi_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   call compare_3D_cplx_results(n1, n2, n3, psi_in, psi_cuda, maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * (log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   
   !Poisson Solver - periodic boundary
    write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Periodic), dimensions:',n1,n2,n3
   !calculate the kernel in parallel for each processor
   call createKernel(0,1,'P',n1,n2,n3,0.2d0,0.2d0,0.2d0,16,pkernel,.false.,0) 

   call nanosec(tsc0);
   call H_potential('P','D',0,1,n1,n2,n3,0.2d0,0.2d0,0.2d0,&
        rhopot,pkernel,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)


   !here the GPU part
   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Periodic), dimensions:',n1,n2,n3

   call transpose_kernel_forGPU('P',n1,n2,n3,pkernel,pkernel2)

   geo1=1
   geo2=1
   geo3=1

   scal=1.d0/dble(n1*n2*n3)

   size1=n1*n2*n3
   size2=2*n1*n2*n3
   sizek=(n1/2+1)*n2*n3
   call cuda_3d_psolver_general_plan(n1,n2,n3,plan1,plan1_,plan2,plan3,plan3_,switch_alg,geo1,geo2,geo3)
   call cudamalloc(size2,work_GPU)
   call cudamalloc(size2,psi_GPU)
   call cudamalloc(sizek,k_GPU)
   call reset_gpu_data(size1,rhopot2,work_GPU)
   call reset_gpu_data(sizek,pkernel2,k_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_psolver_general(n1,n2,n3,plan1,plan1_,plan2,plan3,plan3_,work_GPU,psi_GPU,k_GPU,switch_alg,geo1,geo2,geo3,scal)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size1,rhopot2,work_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)
   call cudafree(k_GPU)
   i_all=-product(shape(pkernel))*kind(pkernel)
   deallocate(pkernel,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)
   i_all=-product(shape(pkernel2))*kind(pkernel2)
   deallocate(pkernel2,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)
 
   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1, n2, n3, rhopot(1), rhopot2(1), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   !Poisson Solver - Free boundary
    write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Free), dimensions:',n1,n2,n3

   !initialize rhopots
   call vcopy(n1*n2*n3/8,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3/8,psi_in(1,1,1,1),2,rhopot2(1),1)

   !calculate the kernel in parallel for each processor
   call createKernel(0,1,'F',n1/2,n2/2,n3/2,1.0d0,1.0d0,1.0d0,16,pkernel,.false.,0)
 
   call nanosec(tsc0);
   call H_potential('F','D',0,1,n1/2,n2/2,n3/2,1.0d0,1.0d0,1.0d0,&
        rhopot,pkernel,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)


   !here the GPU part
   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Free), dimensions:',n1,n2,n3

   call createKernel(0,1,'F',n1/2,n2/2,n3/2,1.0d0,1.0d0,1.0d0,16,pkernel,.false.,1)

   call transpose_kernel_forGPU_free('F',n1/2,n2/2,n3/2,pkernel,pkernel2)

   geo1=0
   geo2=0
   geo3=0

   scal=(1.0d0*1.0d0*1.0d0)/dble(n1*n2*n3)

   size1=n1*n2*n3/8
   size2=2*n1*n2*n3
   sizek=(n1/2+1)*n2*n3
   call cuda_3d_psolver_general_plan(n1,n2,n3,plan1,plan1_,plan2,plan3,plan3_,switch_alg,geo1,geo2,geo3)
   call cudamalloc(size2,work_GPU)
   call cudamalloc(size2,psi_GPU)
   call cudamalloc(sizek,k_GPU)
   call reset_gpu_data(size1,rhopot2,work_GPU)
   call reset_gpu_data(sizek,pkernel2,k_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_psolver_general(n1,n2,n3,plan1,plan1_,plan2,plan3,plan3_,work_GPU,psi_GPU,k_GPU,switch_alg,geo1,geo2,geo3,scal)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size1,rhopot2,work_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)
   call cudafree(k_GPU)
   i_all=-product(shape(pkernel))*kind(pkernel)
   deallocate(pkernel,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)
   i_all=-product(shape(pkernel2))*kind(pkernel2)
   deallocate(pkernel2,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)
 
   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1/2, n2/2, n3/2, rhopot(1), rhopot2(1), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)


   !Poisson Solver - Surface boundary
    write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver (Surface), dimensions:',n1,n2,n3

   !initialize rhopots
   call vcopy(n1*n2*n3/2,psi_in(1,1,1,1),2,rhopot(1),1)
   call vcopy(n1*n2*n3/2,psi_in(1,1,1,1),2,rhopot2(1),1)

   !calculate the kernel in parallel for each processor
   call createKernel(0,1,'S',n1,n2/2,n3,0.2d0,0.2d0,0.2d0,16,pkernel,.false.,0)
   !pkernel(1:size(pkernel))=1.0_dp

   call nanosec(tsc0);
   call H_potential('S','D',0,1,n1,n2/2,n3,0.2d0,0.2d0,0.2d0,&
        rhopot,pkernel,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)

   !here the GPU part
   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver (Surface), dimensions:',n1,n2,n3

   call createKernel(0,1,'S',n1,n2/2,n3,0.2d0,0.2d0,0.2d0,16,pkernel,.false.,1)

   call transpose_kernel_forGPU_free('S',n1,n2/2,n3,pkernel,pkernel2)

   geo1=1
   geo2=0
   geo3=1

   scal=(-16.0_dp*atan(1.0_dp)*dble(0.2_dp))/dble(n1*n2*n3)

   size1=n1*n2*n3/2
   size2=2*n1*n2*n3
   sizek=(n1/2+1)*n2*n3
   call cuda_3d_psolver_general_plan(n1,n2,n3,plan1,plan1_,plan2,plan3,plan3_,switch_alg,geo1,geo2,geo3)
   call cudamalloc(size2,work_GPU)
   call cudamalloc(size2,psi_GPU)
   call cudamalloc(sizek,k_GPU)
   call reset_gpu_data(size1,rhopot2,work_GPU)
   call reset_gpu_data(sizek,pkernel2,k_GPU)

   call nanosec(tsc0);
   do i=1,ntimes
     call cuda_3d_psolver_general(n1,n2,n3,plan1,plan1_,plan2,plan3,plan3_,work_GPU,psi_GPU,k_GPU,switch_alg,geo1,geo2,geo3,scal)
   end do
   call synchronize()
   call nanosec(tsc1)

   call get_gpu_data(size1,rhopot2,work_GPU)

   call cudafree(work_GPU)
   call cudafree(psi_GPU)
   call cudafree(k_GPU)
   i_all=-product(shape(pkernel))*kind(pkernel)
   deallocate(pkernel,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)
   i_all=-product(shape(pkernel2))*kind(pkernel2)
   deallocate(pkernel2,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)


   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1, n2/2, n3, rhopot(1), rhopot2(1), maxdiff, 3.0d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

  i_all=-product(shape(rhopot))*kind(rhopot)
  deallocate(rhopot,stat=i_stat)
  call memocc(i_stat,i_all,'rhopot',subname)
  i_all=-product(shape(rhopot2))*kind(rhopot2)
  deallocate(rhopot2,stat=i_stat)
  call memocc(i_stat,i_all,'rhopot2',subname)

contains

  subroutine print_time(time,nbelem,nop,ntimes)
    implicit none
    real(gp),intent(in)::time
    integer,intent(in)::nbelem,ntimes
    real(gp),intent(in)::nop

    write(*,'(a,f9.4,1pe12.5)')'Finished. Time(ms), GFlops',&
      time*1.d3/real(ntimes,kind=8),&
      real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9)

  end subroutine print_time

  subroutine compare_time(REFtime,TESTtime,nbelem,nop,ntimes,maxdiff,threshold)
    implicit none
    real(gp),intent(in)::REFtime,TESTtime,maxdiff,threshold
    integer,intent(in)::nbelem,ntimes
    real(gp),intent(in)::nop

    write(*,'(a,i10,f9.5,1pe12.5,2(0pf12.4,0pf12.4))',advance='no')&
      'nbelem,REF/TEST ratio,Time,Gflops: REF,TEST',&
       nbelem,REFtime/TESTtime,maxdiff,&
       REFtime*1.d3/real(ntimes,kind=8),&
       real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(REFtime*1.d9),&
       TESTtime*1.d3/real(ntimes,kind=8),&
       real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(TESTtime*1.d9)
    if (maxdiff <= threshold) then
      write(*,'(a)')''
    else
      write(*,'(a)')'<<<< WARNING' 
    end if
  end subroutine compare_time

  subroutine compare_3D_results(dim1, dim2, dim3, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2, dim3
    real(gp),intent(in):: psi_ref(dim1,dim2,dim3), psi(dim1,dim2,dim3)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,i3

    maxdiff=0.d0
    do i3=1,dim3
      do i2=1,dim2
        do i1=1,dim1
          comp=abs(psi_ref(i1,i2,i3)-psi(i1,i2,i3))
          if(comp > printdiff) then
            write(*,*)i3,i2,i1,psi_ref(i1,i2,i3),psi(i1,i2,i3)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
        end do
      end do
    end do
  end subroutine compare_3D_results

  subroutine compare_3D_cplx_results(dim1, dim2, dim3, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2, dim3
    real(gp),intent(in):: psi_ref(2,dim1,dim2,dim3), psi(2,dim1,dim2,dim3)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,i3

    maxdiff=0.d0
    do i3=1,dim3
      do i2=1,dim2
        do i1=1,dim1
          comp=abs(psi_ref(1,i1,i2,i3)-psi(1,i1,i2,i3))
          if(comp > printdiff) then
            write(*,*)i3,i2,i1,1,psi_ref(1,i1,i2,i3),psi(1,i1,i2,i3)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
          comp=abs(psi_ref(2,i1,i2,i3)-psi(2,i1,i2,i3))
          if(comp > printdiff) then
            write(*,*)i3,i2,i1,2,psi_ref(2,i1,i2,i3),psi(2,i1,i2,i3)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if

        end do
      end do
    end do
  end subroutine compare_3D_cplx_results


  subroutine compare_2D_results(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(dim1,dim2), psi(dim1,dim2)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        comp=abs(psi_ref(i1,i2)-psi(i1,i2))
        if(comp > printdiff) then
          write(*,*)i2,i1,psi_ref(i1,i2),psi(i1,i2)
        endif
        if (comp > maxdiff) then
          maxdiff=comp
        end if
      end do
    end do
  end subroutine compare_2D_results

  subroutine compare_2D_cplx_results(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(2,dim1,dim2), psi(2,dim1,dim2)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,c

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        do c=1,2
          comp=abs(psi_ref(c,i1,i2)-psi(c,i1,i2))
          if(comp > printdiff) then
            write(*,*)c,i2,i1,psi_ref(c,i1,i2),psi(c,i1,i2)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
        end do
      end do
    end do
  end subroutine compare_2D_cplx_results


  subroutine compare_2D_results_t(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(dim1,dim2), psi(dim2,dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        comp=abs(psi_ref(i1,i2)-psi(i2,i1))
        if(comp > printdiff) then
          write(*,*)i2,i1,psi_ref(i1,i2),psi(i2,i1)
        endif
        if (comp > maxdiff) then
          maxdiff=comp
        end if
      end do
    end do
  end subroutine compare_2D_results_t

  subroutine compare_2D_cplx_results_t(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(2,dim1,dim2), psi(2,dim2,dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2,c

    maxdiff=0.d0
    do i2=1,dim2
      do i1=1,dim1
        do c=1,2
          comp=abs(psi_ref(c,i1,i2)-psi(c,i2,i1))
          if(comp > printdiff) then
            write(*,*)c,i2,i1,psi_ref(c,i1,i2),psi(c,i2,i1)
          endif
          if (comp > maxdiff) then
            maxdiff=comp
          end if
        enddo
      end do
    end do
  end subroutine compare_2D_cplx_results_t


  subroutine compare_1D_results(dim1, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1
    real(gp),intent(in):: psi_ref(dim1), psi(dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1

    maxdiff=0.d0
    do i1=1,dim1
      comp=abs(psi_ref(i1)-psi(i1))
      if(comp > printdiff) then
        write(*,*)i1,psi_ref(i1),psi(i1)
      endif
      if (comp > maxdiff) then
        maxdiff=comp
      end if
    end do
  end subroutine compare_1D_results

  subroutine transpose_kernel_forGPU(geocode,n01,n02,n03,pkernel,pkernel2)
    use module_base
    use Poisson_Solver
    implicit none
    integer, intent(in) :: n01,n02,n03
    character(len=*), intent(in) :: geocode
    real(dp), dimension(*), intent(in) :: pkernel
    real(dp), dimension(:), pointer :: pkernel2
    !local variables
    character(len=*), parameter :: subname='transpose_kernel_forGPU'
    integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,nproc,i1,i2,i3,ind,indt
    integer :: j1,j2,j3,i_stat

    !only makes sense in serial for the moment
    nproc=1
    !calculate the dimensions wrt the geocode
    if (geocode == 'P') then
       call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
    else if (geocode == 'S') then
       call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0)
    else if (geocode == 'F') then
       call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0)
    else
       stop 'ERROR(transpose_kernel_forGPU): geometry code not admitted'
    end if

    allocate(pkernel2(n1*n2*n3+ndebug),stat=i_stat)
    call memocc(i_stat,pkernel2,'pkernel2',subname)

    do i3=1,n3
       j3=i3+(i3/(n3/2+2))*(n3+2-2*i3)!injective dimension
       do i2=1,n2
          j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)!injective dimension
          do i1=1,n1
             j1=i1+(i1/(n1/2+2))*(n1+2-2*i1)!injective dimension
             !injective index
             ind=j1+(j2-1)*nd1+(j3-1)*nd1*nd2
             !unfolded index
             !indt=j1+(i3-1)*nd1+(i2-1)*nd1*n3             
             !indt=i2+(i1-1)*n2+(i3-1)*n1*n2             
             indt=i2+(j1-1)*n2+(i3-1)*nd1*n2
             pkernel2(indt)=pkernel(ind)
          end do
       end do
    end do
    !offset to zero
    pkernel2(1)=0.0_dp

  end subroutine transpose_kernel_forGPU
 
  subroutine transpose_kernel_forGPU_free(geocode,n01,n02,n03,pkernel,pkernel2)
    use module_base
    use Poisson_Solver
    implicit none
    integer, intent(in) :: n01,n02,n03
    character(len=*), intent(in) :: geocode
    real(dp), dimension(*), intent(in) :: pkernel
    real(dp), dimension(:), pointer :: pkernel2
    !local variables
    character(len=*), parameter :: subname='transpose_kernel_forGPU'
    integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3,nproc,i1,i2,i3,ind,indt
    integer :: j1,j2,j3,i_stat

    !only makes sense in serial for the moment
    nproc=1
    !calculate the dimensions wrt the geocode
    if (geocode == 'P') then
       call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
    else if (geocode == 'S') then
       call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1)
    else if (geocode == 'F') then
       call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,1)
    else
       stop 'ERROR(transpose_kernel_forGPU): geometry code not admitted'
    end if

    allocate(pkernel2(n1*n2*n3+ndebug),stat=i_stat)
    call memocc(i_stat,pkernel2,'pkernel2',subname)

    do i3=1,n3
       j3=i3+(i3/(n3/2+2))*(n3+2-2*i3)!injective dimension
       do i2=1,n2
          j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)!injective dimension
          do i1=1,n1
             j1=i1+(i1/(n1/2+2))*(n1+2-2*i1)!injective dimension
             !injective index
             ind=j1+(j2-1)*nd1+(j3-1)*nd1*nd2
             !unfolded index
             !indt=j1+(i3-1)*nd1+(i2-1)*nd1*n3
             !indt=i2+(i1-1)*n2+(i3-1)*n1*n2
             indt=i2+(j1-1)*n2+(i3-1)*nd1*n2
             pkernel2(indt)=pkernel(ind)
          end do
       end do
    end do
    !offset to zero
    !pkernel2(1)=0.0_dp

  end subroutine transpose_kernel_forGPU_free

end program conv_check_fft


!!***
