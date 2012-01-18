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
  integer  :: n1,n2,n3 !n(c) n1bis,n2bis,n3bis
  real(gp) :: hx,r2,sigma2,x,y,z,maxdiff,arg !n(c) epot,hy,hz
  real(wp), dimension(:,:,:,:), allocatable :: psi_in,psi_out !n(c) pot,psir,psi_out_s
  !n(c) real(wp), dimension(:,:,:,:,:), allocatable :: psi_k_in, psi_k_out
  !n(c) real(wp), dimension(:,:,:,:), allocatable :: psi_k_in_a, psi_k_out_a,pot_a
  !n(c) real(wp), dimension(:,:,:,:), allocatable :: psi_in_s,psi_out_t,psi_in_t,psi_gemm,psi_gemmsy,psi_cuda_gemm
  !local variables
  character(len=*), parameter :: subname='conv_check_fft'
  !n(c) character(len=50) :: chain
  integer :: i,i_stat,j,i1,i2,i3,ntimes,i_all !n(c) i_all,i1_max,i_max,it0,it1,ndim,itimes
  integer :: l,ierror,device_number !n(c) i1s,i1e,count_rate,count_max
  real(wp) :: tt,scale
  real(gp) :: CPUtime,GPUtime,ekin,ehartree !n(c) v,p,comp
  !n(c) real(gp), dimension(3) :: hgridh
  !n(c) real(gp), dimension(8) :: scal
  integer, dimension(:), allocatable :: modarr !n(c) keyv
  !n(c) integer, dimension(:,:), allocatable :: keyg
  real(kind=8), dimension(:), allocatable :: psi !temporary in view of wp !n(c) psi_d
  !n(c) real(kind=4), dimension(:), allocatable :: psi_l !temporary in view of wp 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda,v_cuda,v_cuda_str !temporary in view of wp !n(c) psi_cuda_s,psi_cuda_t,v_cuda_s,v_cuda_t
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_str 
  !n(c) real(kind=8), dimension(:,:,:,:,:), allocatable :: psi_cuda_k_in,psi_cuda_k_out,psi_cuda_k_in_bis 
  !n(c) real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_k_in_a,psi_cuda_k_out_a 
  real(kind=4), dimension(:,:,:,:), allocatable :: psi_cuda_l,v_cuda_l !temporary in view of wp 
  real(kind=8) :: ekinGPUd,ehartreeGPUd
  !n(c) real(kind=4) :: t0,t1,epotGPU,ekinGPU
  real(kind=8) :: psi_GPU,work_GPU,work2_GPU,k_GPU !pointer to the GPU  memory addresses (with norb=1) !n(c) keys_GPU,v_GPU
  !n(c) real(kind=8) :: psi_c_GPU, psi_f_GPU,keyg_GPU, keyv_GPU
  real(kind=8) :: context,queue
  !n(c) integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
  !n(c) integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
  integer, parameter :: lowfilK=-14,lupfilK=14 ! kinetic term
  real(kind=8), dimension(lowfilK:lupfilK) :: fil
  integer(kind=8) :: tsc0, tsc1
  !objects for the 3D Poisson solver
  real(kind=8), dimension(:), allocatable :: rhopot,rhopot2
  real(kind=8), dimension(:), pointer :: pkernel,pkernel2
 
!!!  !Use arguments
!!!  call getarg(1,chain)
!!!  read(unit=chain,fmt=*) n1
!!!  call getarg(2,chain)
!!!  read(unit=chain,fmt=*) n2*n3

  read(unit=2,fmt=*,iostat=ierror) n1,n2,n3,ntimes
  if (ierror /= 0) then
     write(*,*) "In a file 'fort.1', put a line with:"
     write(*,*) "n1 n2 n3 ntimes"
     write(*,*) "where:"
     write(*,*) "- n1 n2 n3 are the dimensions of the real space"
     write(*,*) "- ntimes is the number of convolutions"
     stop
  end if

  call ocl_create_gpu_context(context,device_number)
  call customize_fft((/n1,n2,n3/))
  call ocl_build_programs(context)
  call ocl_create_command_queue(queue,context)
  call init_event_list

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
   allocate(psi_cuda_l(2,n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_l,'psi_cuda_l',subname)
   allocate(v_cuda_l(2,n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_l,'v_cuda_l',subname)
   allocate(rhopot(n1*n2*n3+ndebug),stat=i_stat)
   call memocc(i_stat,rhopot,'rhopot',subname)
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
          v_cuda_str(i2,i1,i,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_str(i2,i1,i,2)=0.0
          v_cuda(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
          v_cuda_l(i2,i,i1,1)=real(psi_in(i2,i1,i,1),kind=8)
        end do
      end do
   end do


  write(*,'(a,i6,i6)')'CPU FFT, dimensions:',n1,n2*n3

   call nanosec(tsc0);
   i3=1
   call fft_1d_ctoc(-1,n2*n3,n1,v_cuda_str,i3)
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   write(*,'(a,i6,i6)')'GPU FFT, dimensions:',n1,n2*n3

   call ocl_create_write_buffer(context, 2*n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, 2*n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, 2*n1*n2*n3*8, psi_in)!v_cuda)!

   call nanosec(tsc0);
   do i=1,ntimes
      call fft1d_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, 2*n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

!   call compare_2D_cplx_results( n1, n2*n3, v_cuda_str(1,1,1,i3), psi_cuda, maxdiff, 3.d-7)
   call compare_2D_cplx_results_t( n1, n2*n3,  psi_cuda, v_cuda_str(1,1,1,i3), maxdiff, 3.d-7)
!           call compare_2D_cplx_results(n1, n2*n3, psi_in, psi_cuda, maxdiff, 3.d-7)
!           call compare_2D_cplx_results_t(n2*n3, n1, v_cuda, psi_cuda, maxdiff, 3.d-7)

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

   call ocl_create_read_write_buffer(context, 2*n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, 2*n1*n2*n3*8, work_GPU)
   call ocl_create_read_write_buffer(context, 2*n1*n2*n3*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, 2*n1*n2*n3*8, psi_in)!v_cuda)!

   call nanosec(tsc0);
   do i=1,ntimes
      call fft3d_d(queue,(/n1,n2,n3/),work_GPU,psi_GPU,work2_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, 2*n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

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

   call ocl_create_read_write_buffer(context, 2*n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, 2*n1*n2*n3*8, work_GPU)
   call ocl_create_read_write_buffer(context, 2*n1*n2*n3*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, 2*n1*n2*n3*8, psi_cuda)!v_cuda)!

   call nanosec(tsc0);
   do i=1,ntimes
      call fft3d_r_d(queue,(/n1,n2,n3/),work_GPU,psi_GPU,work2_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, 2*n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)

   call compare_3D_cplx_results(n1, n2, n3, psi_in, psi_cuda, maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,5 * (log(real(n1,kind=8))+&
     log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)

   
   !Poisson Solver
    write(*,'(a,i6,i6,i6)')'CPU 3D Poisson Solver, dimensions:',n1,n2,n3
   !calculate the kernel in parallel for each processor
   call createKernel(0,1,'P',n1,n2,n3,0.2d0,0.2d0,0.2d0,16,pkernel,.false.) 

   !call to_zero(size(pkernel),pkernel(1))
   !pkernel(1:size(pkernel))=1.0_dp

   call nanosec(tsc0);
   call H_potential('P','D',0,1,n1,n2,n3,0.2d0,0.2d0,0.2d0,&
        rhopot,pkernel,rhopot,ehartree,0.0d0,.false.,quiet='yes') !optional argument
   call nanosec(tsc1);
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9*ntimes
   call print_time(CPUtime,n1*n2*n3,5 *( log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes)

   !here the GPU part
   write(*,'(a,i6,i6,i6)')'GPU 3D Poisson Solver, dimensions:',n1,n2,n3

   !transpose the kernel before copying
   call transpose_kernel_forGPU('P',n1,n2,n3,pkernel,pkernel2)
   !pkernel2(1:size(pkernel2))=1.0_dp
   call ocl_create_read_write_buffer(context, 2*n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_create_read_write_buffer(context, 2*n1*n2*n3*8, work2_GPU)
   call ocl_create_read_write_buffer(context, n1*n2*n3*8, k_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, rhopot2)!v_cuda)!
   call ocl_enqueue_write_buffer(queue, k_GPU, n1*n2*n3*8, pkernel2)!v_cuda)!

   call nanosec(tsc0);
   do i=1,ntimes
      call fft3d_k_r2c_d(queue,(/n1,n2,n3/),work_GPU,psi_GPU,work2_GPU,k_GPU)
      call fft3d_r_c2r_d(queue,(/n1,n2,n3/),psi_GPU,work2_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, work2_GPU, n1*n2*n3*8, rhopot2(1))
   call ocl_release_mem_object(k_GPU)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time(GPUtime,n1*n2*n3*3,5 * log(real(n1,kind=8))/log(real(2,kind=8)),ntimes)
   call compare_3D_results(n1, n2, n3, rhopot(1), rhopot2(1), maxdiff, 3.d-7)
   call compare_time(CPUtime,GPUtime,n1*n2*n3,2*5 * (log(real(n1,kind=8))+&
        log(real(n2,kind=8))+log(real(n3,kind=8)))/log(real(2,kind=8)),ntimes,maxdiff,3.d-7)



   i_all=-product(shape(pkernel))*kind(pkernel)
   deallocate(pkernel,stat=i_stat)
   call memocc(i_stat,i_all,'pkernel',subname)
  i_all=-product(shape(rhopot))*kind(rhopot)
  deallocate(rhopot,stat=i_stat)
  call memocc(i_stat,i_all,'rhopot',subname)
  i_all=-product(shape(rhopot2))*kind(rhopot2)
  deallocate(rhopot2,stat=i_stat)
  call memocc(i_stat,i_all,'rhopot2',subname)
  i_all=-product(shape(pkernel2))*kind(pkernel2)
  deallocate(pkernel2,stat=i_stat)
  call memocc(i_stat,i_all,'pkernel2',subname)


  call print_event_list
  call ocl_clean_command_queue(queue)
  call ocl_clean(context)

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
       call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
    else if (geocode == 'F') then
       call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
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
             indt=i1+(i3-1)*n1+(i2-1)*n1*n3             
             pkernel2(indt)=pkernel(ind)
          end do
       end do
    end do
    !offset to zero
    pkernel2(1)=0.0_dp

  end subroutine transpose_kernel_forGPU
 
end program conv_check_fft


!!***
