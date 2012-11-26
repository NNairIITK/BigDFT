!> @file
!!   Test program for convolutions with OpenCL
!! @author
!!    Copyright (C) 2008-2012 BigDFT group (LG) (BV)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!


!> Program test for the convolution in GPU (OpenCL version)
program conv_check_ocl
  use yaml_output
   use module_base
   implicit none
   integer  :: n1,n2,n3,n1bis,n2bis,n3bis
   real(gp) :: hx,r2,sigma2,x,maxdiff,arg !n(c) hy,hz
   real(wp), dimension(:,:,:), allocatable :: psi_in,psi_out,psi_out_s
   real(wp), dimension(:,:,:,:,:), allocatable :: psi_k_in, psi_k_out
   real(wp), dimension(:,:,:,:), allocatable :: psi_k_in_a, psi_k_out_a, pot_a
   real(wp), dimension(:,:,:), allocatable :: psi_in_s,psi_out_t,psi_in_t,psi_gemm,psi_gemmsy,psi_cuda_gemm
   !local variables
   character(len=*), parameter :: subname='conv_check'
   integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,itimes !n(c) ndim
   integer :: l,ierror,i1s,i1e,device_number
   integer :: nvctr_cf,nseg,iseg !n(c) n1s,n1e,ndats,ndate
   real(wp) :: tt,scale
   real(gp) :: CPUtime,GPUtime,ekin
   real(gp), dimension(8) :: scal
   integer, dimension(:), allocatable :: keyv,modarr
   integer, dimension(:,:), allocatable :: keyg
   real(kind=8), dimension(:), allocatable :: psi, psi_d !temporary in view of wp 
   real(kind=4), dimension(:), allocatable :: psi_l !temporary in view of wp 
   real(kind=8), dimension(:,:,:), allocatable :: psi_cuda,v_cuda,psi_cuda_s,v_cuda_s,psi_cuda_t,v_cuda_t,v_cuda_str !temporary in view of wp 
   real(kind=8), dimension(:,:,:), allocatable :: psi_cuda_str 
   real(kind=8), dimension(:,:,:,:,:), allocatable :: psi_cuda_k_in,psi_cuda_k_out,psi_cuda_k_in_bis 
   real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_k_in_a,psi_cuda_k_out_a 
   real(kind=4), dimension(:,:,:), allocatable :: psi_cuda_l,v_cuda_l !temporary in view of wp 
   real(kind=8) :: ekinGPUd
   real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU !pointer to the GPU  memory addresses (with norb=1)
   real(kind=8) :: psi_c_GPU, psi_f_GPU, keyg_GPU, keyv_GPU
   real(kind=8) :: context,queue
   !n(c) integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
   !n(c) integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
   integer, parameter :: lowfilK=-14,lupfilK=14 ! kinetic term
   real(kind=8), dimension(lowfilK:lupfilK) :: fil
   integer(kind=8) :: tsc0, tsc1
   character(len=500) :: field

   field=repeat(' ',len(field))

   !!!  !Use arguments
   !!!  call getarg(1,chain)
   !!!  read(unit=chain,fmt=*) n1
   !!!  call getarg(2,chain)
   !!!  read(unit=chain,fmt=*) ndat

   read(unit=1,fmt=*,iostat=ierror) n1,n2,n3,ntimes
   if (ierror /= 0) then
      write(*,*) "In a file 'fort.1', put a line with:"
      write(*,*) "n1 n2 n3 ntimes"
      write(*,*) "where:"
      write(*,*) "- n1 n2 n3 are the dimensions of the real space"
      write(*,*) "- ntimes is the number of convolutions"
      stop
   end if

   call yaml_new_document()

   call ocl_create_gpu_context(context,device_number)
   call ocl_build_programs(context)
   call ocl_create_command_queue(queue,context)
   call init_event_list

   hx=0.1e0_gp
   !n(c) hy=0.1e0_gp
   !n(c) hz=0.1e0_gp

   scale=real(-.5_gp/hx**2,wp)

   ! second derivative filters for Daubechies 16
   fil(0)=   -3.5536922899131901941296809374e0_wp*scale
   fil(1)=    2.2191465938911163898794546405e0_wp*scale
   fil(2)=   -0.6156141465570069496314853949e0_wp*scale
   fil(3)=    0.2371780582153805636239247476e0_wp*scale
   fil(4)=   -0.0822663999742123340987663521e0_wp*scale
   fil(5)=    0.02207029188482255523789911295638968409e0_wp*scale
   fil(6)=   -0.409765689342633823899327051188315485e-2_wp*scale
   fil(7)=    0.45167920287502235349480037639758496e-3_wp*scale
   fil(8)=   -0.2398228524507599670405555359023135e-4_wp*scale
   fil(9)=    2.0904234952920365957922889447361e-6_wp*scale
   fil(10)=  -3.7230763047369275848791496973044e-7_wp*scale
   fil(11)=  -1.05857055496741470373494132287e-8_wp*scale
   fil(12)=  -5.813879830282540547959250667e-11_wp*scale
   fil(13)=   2.70800493626319438269856689037647576e-13_wp*scale
   fil(14)=  -6.924474940639200152025730585882e-18_wp*scale

   !!!  ! second derivative filters for Daubechies 16
   !!!  fil(0)=    0.e-3_wp*scale
   !!!  fil(1)=    1.e-3_wp*scale
   !!!  fil(2)=    2.e-3_wp*scale
   !!!  fil(3)=    3.e-3_wp*scale
   !!!  fil(4)=    4.e-3_wp*scale
   !!!  fil(5)=    5.e-3_wp*scale
   !!!  fil(6)=    6.e-3_wp*scale
   !!!  fil(7)=    7.e-3_wp*scale
   !!!  fil(8)=    8.e-3_wp*scale
   !!!  fil(9)=    9.e-3_wp*scale
   !!!  fil(10)=  10.e-3_wp*scale
   !!!  fil(11)=  11.e-3_wp*scale
   !!!  fil(12)=  12.e-3_wp*scale
   !!!  fil(13)=  13.e-3_wp*scale
   !!!  fil(14)=  14.e-3_wp*scale


   do i=1,14
      fil(-i)=fil(i)
   enddo

   ekin=0.0_wp

   !allocate arrays
   allocate(psi_in(n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in,'psi_in',subname)
   allocate(psi_out(n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_out,'psi_out',subname)
   allocate(psi_gemm(n1,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in,'psi_gemm',subname)
   allocate(psi_gemmsy(n1,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in,'psi_gemmsy',subname)
   allocate(psi_cuda_gemm(n1,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_out,'psi_cuda_gemm',subname)

   allocate(psi_cuda(n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda,'psi_cuda',subname)
   allocate(psi_cuda_str(n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_str,'psi_cuda_str',subname)
   allocate(v_cuda(n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda,'v_cuda',subname)
   allocate(v_cuda_str(n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_str,'v_cuda_str',subname)
   allocate(psi_cuda_l(n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_l,'psi_cuda_l',subname)
   allocate(v_cuda_l(n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_l,'v_cuda_l',subname)

   !initialise array
   sigma2=0.25d0*((n1*hx)**2)
   do i=1,n2*n3
      do i1=1,n1
         x=hx*real(i1-n1/2-1,kind=8)
         !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
         r2=x**2
         arg=0.5d0*r2/sigma2
         tt=dexp(-arg)
         call random_number(tt)
         psi_in(i1,i,1)=tt
      end do
   end do

   !the input and output arrays must be reverted in this implementation
   do i=1,n2*n3
      do i1=1,n1
         v_cuda_str(i1,i,1)=real(psi_in(i1,i,1),kind=8)
         v_cuda(i,i1,1)=real(psi_in(i1,i,1),kind=8)
         v_cuda_l(i,i1,1)=real(psi_in(i1,i,1),kind=4)
      end do
   end do

!   write(*,'(a,i10)')'GPU FLOPS double Benchmark, dimension:',n1*n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call benchmark_flops_d(queue,n1*n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time('GPU FLOPS double Benchmark',(/n1*n2*n3/),&
        GPUtime,n1*n2*n3,16*128*2,ntimes)

!   write(*,'(a,i10)')'GPU MOPS double Benchmark, dimension:',n1*n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call benchmark_mops_d(queue,n1*n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time('GPU MOPS double Benchmark',(/n1*n2*n3/),&
        GPUtime,n1*n2*n3,8*2,ntimes)

!   write(*,'(a,i10,i10)')'GPU transpose FLOPS double Benchmark, dimension:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call transpose_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time('GPU transpose FLOPS double Benchmark',(/n1,n2*n3/),&
        GPUtime,n1*n2*n3,32,ntimes)

!   write(*,'(a,i10,i10)')'GPU notranspose FLOPS double Benchmark, dimension:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call notranspose_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time('GPU notranspose FLOPS double Benchmark',(/n1,n2*n3/),&
        GPUtime,n1*n2*n3,32,ntimes)

!   write(*,'(a,i10,i10)')'GPU copy FLOPS double Benchmark, dimension:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call copy_d(queue,n1*n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   call print_time('GPU copy FLOPS double Benchmark',(/n1,n2*n3/),&
        GPUtime,n1*n2*n3,32,ntimes)



   !write(*,'(a,i6,i6)')'CPU Convolutions, dimensions:',n1,n2*n3

   call nanosec(tsc0);
   do i=1,ntimes
      call convrot_n_per(n1-1,n2*n3,psi_in,psi_out)
   end do
   call nanosec(tsc1);

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


   !write(field,'(a,i6,i6)')'Convolutions, dimensions:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call magicfilter1d_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results_t(n2*n3, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

   call compare_time('Convolutions',(/n1,n2*n3/),&
        CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,3.d-7)
   field=repeat(' ',len(field))

   !write(field,'(a,i6,i6)')'Convolutions (straight), dimensions:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda_str)

   call nanosec(tsc0)
   do i=1,ntimes
      call magicfilter1d_straight_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue)
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda_str)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results(n2*n3, n1, psi_out, psi_cuda_str, maxdiff, 3.d-7)

   call compare_time('Convolutions (str)',(/n1,n2*n3/),&
        CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,3.d-7)

!   field=repeat(' ',len(field))
!   write(field,'(a,i6)')'Convolutions (block)'

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0);
   do i=1,ntimes
      call magicfilter1d_block_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1);

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results_t(n2*n3, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

   call compare_time('Convolutions (block)',(/n1,n2*n3/),&
        CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,3.d-7)


!   write(*,'(a,i6,i6)')'CPU Convolutions T, dimensions:',n1,n2*n3

   call nanosec(tsc0)
   do i=1,ntimes
      call convrot_t_per(n1-1,n2*n3,psi_in,psi_out)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!   write(*,'(a,i6,i6)')'GPU Convolutions T, dimensions:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)


   call nanosec(tsc0)
   do i=1,ntimes
      call magicfilter1d_t_d(queue,n1,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results_t(n2*n3, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

   call compare_time('Convolutions T',(/n1,n2*n3/),&
        CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,3.d-7)

!   write(*,'(a,i6,i6,i6)')'CPU GEMM, dimensions:',n1,n1,n2*n3

   call nanosec(tsc0)
   do i=1,ntimes
      call DGEMM('n','n',n1,n1,n2*n3,1.2d0, psi_in(:,:,1), n1, psi_in(:,:,1), n2*n3, 0.0d0, psi_gemm(:,:,1), n1)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9

   do j=1,n1
      do i=1,j
         psi_gemmsy(i,j,1) = psi_gemm(j,i,1)
         psi_gemmsy(j,i,1) = psi_gemm(j,i,1)
      end do
   end do

!   write(*,'(a,i6,i6,i6)')'GPU GEMM, dimensions:',n1,n1,n2*n3

   call ocl_create_write_buffer(context, n1*n1*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda_str)

   call nanosec(tsc0)
   do i=1,ntimes
      call gemm_d(queue,'n','n',n1,n1,n2*n3,1.2d0,work_GPU,n1,work_GPU,n2*n3,0.0d0,psi_GPU, n1)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n1*8, psi_cuda_gemm)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results(n1, n1, psi_gemm, psi_cuda_gemm, maxdiff, 3.d-7)

   call compare_time('GEMM',(/n1,n1,n2*n3/),&
        CPUtime,GPUtime,n1*n1,n2*n3*2,ntimes,maxdiff,3.d-7)

   !           write(*,'(a,i6,i6,i6)')'GPU GEMM (volkov), dimensions:',n1,n1,n2*n3
   !
   !           call ocl_create_write_buffer(context, n1*n1*8, psi_GPU)
   !           call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   !           call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda_str)
   !
   !           call nanosec(tsc0)
   !           do i=1,ntimes
   !              call gemm_volkov_d(queue,'n','n',n1,n1,n2*n3,1.2d0,work_GPU,n1,work_GPU,n2*n3,0.0d0,psi_GPU, n1)
   !           end do
   !           call ocl_finish(queue);
   !           call nanosec(tsc1)
   !
   !           call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n1*8, psi_cuda_gemm)
   !           call ocl_release_mem_object(psi_GPU)
   !           call ocl_release_mem_object(work_GPU)
   !
   !           GPUtime=real(tsc1-tsc0,kind=8)*1d-9
   !           call print_time(GPUtime,n1*n1,n2*n3*2,ntimes)
   !
   !           call compare_2D_results(n1, n1, psi_gemm, psi_cuda_gemm, maxdiff, 3.d-7)
   !
   !           call compare_time(field,(/n1/),CPUtime,GPUtime,n1*n1,n2*n3*2,ntimes,maxdiff,3.d-7)

!   write(*,'(a,i6,i6,i6)')'GPU GEMM (block), dimensions:',n1,n1,n2*n3

   call ocl_create_write_buffer(context, n1*n1*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda_str)

   call nanosec(tsc0)
   do i=1,ntimes
      call gemm_block_d(queue,'n','n',n1,n1,n2*n3,1.2d0,work_GPU,n1,work_GPU,n2*n3,0.0d0,psi_GPU, n1)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n1*8, psi_cuda_gemm)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_2D_results(n1, n1, psi_gemm, psi_cuda_gemm, maxdiff, 3.d-7)

   call compare_time('GEMM (block)',(/n1,n1,n2*n3/),&
        CPUtime,GPUtime,n1*n1,n2*n3*2,ntimes,maxdiff,3.d-7)

!   write(*,'(a,i6,i6,i6)')'GPU GEMMSY, dimensions:',n1,n1,n2*n3

   call ocl_create_write_buffer(context, n1*n1*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda_str)

   call nanosec(tsc0)
   do i=1,ntimes
      call gemmsy_d(queue,'n','n',n1,n1,n2*n3,1.2d0,work_GPU,n1,work_GPU,n2*n3,0.0d0,psi_GPU, n1)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n1*8, psi_cuda_gemm)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   !call compare_2D_results(n1, n1, psi_gemmsy, psi_cuda_gemm, maxdiff, 3.d-7)

   call compare_time('GEMMSY',(/n1,n1,n2*n3/),&
        CPUtime,GPUtime,n1*n1,n2*n3*2,ntimes,maxdiff,3.d-7)

!   write(*,'(a,i6,i6,i6)')'CPU ZGEMM, dimensions:',n1/2,n1/2,n2*n3

   call nanosec(tsc0)
   do i=1,ntimes
      call ZGEMM('n','n',n1/2,n1/2,n2*n3,(/1.2d0,1.3d0/),psi_in(:,:,1),n1/2,&
         &   psi_in(:,:,1),n2*n3,&
      (/0.0d0,0.0d0/),psi_gemm(:,:,1),n1/2)
   end do
   call nanosec(tsc1)
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!   write(*,'(a,i6,i6,i6)')'GPU ZGEMM, dimensions:',n1/2,n1/2,n2*n3

   call ocl_create_read_write_buffer(context, (n1/2)*(n1/2)*8*2, psi_GPU)
   call ocl_create_read_buffer(context, (n1/2)*n2*n3*8*2, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, (n1/2)*n2*n3*8*2, v_cuda_str)

   call nanosec(tsc0)
   do i=1,ntimes
      call gemm_z(queue,'n','n',n1/2,n1/2,n2*n3,(/1.2d0,1.3d0/),work_GPU,n1/2,&
         &   work_GPU,n2*n3,&
      (/0.0d0,0.0d0/),psi_GPU, n1/2)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, (n1/2)*(n1/2)*8*2, psi_cuda_gemm)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_2D_results(n1, n1/2, psi_gemm, psi_cuda_gemm, maxdiff, 3.d-7)

   call compare_time('ZGEMM',(/n1/2,n1/2,n2*n3/),&
        CPUtime,GPUtime,n1*n1/4,n2*n3*8,ntimes,maxdiff,3.d-7)

!   write(*,'(a,i8)')'CPU Reduction, dimensions:',n1*n2*n3

   call nanosec(tsc0)
   do i=1,ntimes
      ekin = sum(psi_in*psi_in)
   end do
   call nanosec(tsc1)
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,i8)')'GPU Reduction, dimensions:',n1*n2*n3

   call ocl_create_read_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_write_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_create_read_write_buffer(context, n1*n2*n3*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, psi_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0)
   do i=1,ntimes
      call nrm2sq_d(queue, n1*n2*n3, psi_GPU, work_GPU, work2_GPU, ekinGPUd )
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   maxdiff=abs(ekin/real(n1*n2*n3,kind=8) - ekinGPUd/real(n1*n2*n3,kind=8))

   call compare_time('Reduction',(/n1*n2*n3/),CPUtime,GPUtime,n1*n2*n3,2,ntimes,maxdiff,3.d-7)

!   write(*,'(a,i8)')'CPU Reduction Dot, dimensions:',n1*n2*n3

   call nanosec(tsc0)
   do i=1,ntimes
      ekin = sum(psi_in*psi_in)
   end do
   call nanosec(tsc1)
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,i8)')'GPU Reduction Dot, dimensions:',n1*n2*n3

   call ocl_create_read_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_write_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_create_read_write_buffer(context, n1*n2*n3*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, psi_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0)
   do i=1,ntimes
      call dot_d(queue, n1*n2*n3, psi_GPU, psi_GPU, work_GPU, work2_GPU, ekinGPUd )
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   maxdiff=abs(ekin/real(n1*n2*n3,kind=8) - ekinGPUd/real(n1*n2*n3,kind=8))

   call compare_time('Dot Product',(/n1*n2*n3/),CPUtime,GPUtime,n1*n2*n3,2,ntimes,maxdiff,3.d-7)

   n1bis = n1
   n2bis = n2
   n3bis = n3
!   write(*,'(a,i6,i6,i6)')'CPU Convolutions 3D, dimensions:',n1bis,n2bis,n3bis
   allocate(psi_k_in_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_k_in_a,'psi_k_in_a',subname)
   allocate(psi_cuda_k_in_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_k_in_a,'psi_cuda_k_in_a',subname)
   allocate(pot_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,pot_a,'pot_a',subname)

   allocate(psi_k_out_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_k_out_a,'psi_k_out_a',subname)
   allocate(psi_cuda_k_out_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_k_out_a,'psi_cuda_k_out_a',subname)

   sigma2=0.25d0*((n1bis*hx)**2)
   do i3=1,n3bis
      do i2=1,n2bis
         do i1=1,n1bis
            pot_a(i1,i2,i3,1) = sin(i1+i2+i3+.7d0)
            x=hx*real(i1-n1bis/2-1,kind=8)
            !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
            r2=x**2
            arg=0.5d0*r2/sigma2
            tt=dexp(-arg)
            psi_k_in_a(i1,i2,i3,1)=tt
            psi_cuda_k_in_a(i1,i2,i3,1) = tt;
         end do
      end do
   end do

   call nanosec(tsc0)
   do itimes=1,ntimes
      call convolut_magic_n_per(n1bis-1,n2bis-1,n3bis-1,psi_k_in_a,psi_k_out_a,psi_cuda_k_out_a)
   end do
   call nanosec(tsc1)
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!   write(*,'(a,i6,i6,i6)')'GPU Convolutions 3D, dimensions:',n1bis,n2bis,n3bis


   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

   call nanosec(tsc0)
   do itimes=1,ntimes
      call magicfilter_n_d(queue,(/n1bis,n2bis,n3bis/),work2_GPU, work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

   call compare_time('Convolution 3D',(/n1bis,n2bis,n3bis/),CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,3.d-9)

!   write(*,'(a,i6,i6,i6)')'GPU Convolutions 3D (block), dimensions:',n1bis,n2bis,n3bis


   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

   call nanosec(tsc0)
   do itimes=1,ntimes
      call magicfilter_n_block_d(queue,(/n1bis,n2bis,n3bis/),work2_GPU, work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

   call compare_time('Convolution 3D (B)',(/n1bis,n2bis,n3bis/),CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,3.d-9)


!   write(*,'(a,i6,i6,i6)')'GPU Convolutions 3D (straight), dimensions:',n1bis,n2bis,n3bis


   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

   call nanosec(tsc0)
   do itimes=1,ntimes
      call magicfilter_n_straight_d(queue,(/n1bis,n2bis,n3bis/),work2_GPU, work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

   call compare_time('Convolution 3D (S)',(/n1bis,n2bis,n3bis/),&
        CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,3.d-9)

!   write(*,'(a,i6,i6,i6)')'CPU Convolutions T 3D, dimensions:',n1bis,n2bis,n3bis

   call nanosec(tsc0)
   do itimes=1,ntimes
      call convolut_magic_t_per(n1bis-1,n2bis-1,n3bis-1,psi_k_in_a,psi_k_out_a)
   end do
   call nanosec(tsc1)
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,i6,i6,i6)')'GPU Convolutions T 3D, dimensions:',n1bis,n2bis,n3bis


   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
   call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

   call nanosec(tsc0)
   do itimes=1,ntimes
      call magicfilter_t_d(queue,(/n1bis,n2bis,n3bis/),work2_GPU, work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

   call compare_time('Convolution T 3D',(/n1bis,n2bis,n3bis/),&
        CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,3.d-9)

   if(modulo(n1bis,2)==0 .AND. modulo(n2bis,2)==0 .AND. modulo(n3bis,2)==0) then
      !write(*,'(a,i6,i6,i6)')'CPU Potential, dimensions:',n1bis,n2bis,n3bis

      call nanosec(tsc0)
      do itimes=1,ntimes
         call convolut_magic_n_per(n1bis-1,n2bis-1,n3bis-1,psi_k_in_a,psi_k_out_a,psi_cuda_k_out_a)             
         psi_cuda_k_out_a = psi_k_out_a * pot_a
         call convolut_magic_t_per_self(n1bis-1,n2bis-1,n3bis-1,psi_cuda_k_out_a,psi_k_out_a)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9


      !write(*,'(a,i6,i6,i6)')'GPU Potential, dimensions:',n1bis,n2bis,n3bis

      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
      call ocl_create_read_buffer(context, n1bis*n2bis*n3bis*8, v_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)
      call ocl_enqueue_write_buffer(queue, v_GPU, n1bis*n2bis*n3bis*8, pot_a)

      call nanosec(tsc0)
      do itimes=1,ntimes
         call potential_application_d(queue,(/n1bis/2,n2bis/2,n3bis/2/),work2_GPU, work_GPU,psi_GPU,v_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)
      call ocl_release_mem_object(work2_GPU)
      call ocl_release_mem_object(v_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9


      call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

      call compare_time('Potential 3D',(/n1bis,n2bis,n3bis/),&
           CPUtime,GPUtime,n1bis*n2bis*n3bis,6*32+1,ntimes,maxdiff,3.d-9)

!      write(*,'(a,i6,i6,i6)')'CPU Kinetic 3D, dimensions:',n1bis,n2bis,n3bis
      psi_k_out_a = 0.d0

      call nanosec(tsc0)
      do itimes=1,ntimes
         call convolut_kinetic_per_c(n1bis-1,n2bis-1,n3bis-1,(/0.1d0,0.1d0,0.1d0/),&
              psi_k_in_a,psi_k_out_a,1.d0)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!      write(*,'(a,i6,i6,i6)')'GPU Kinetic 3D, dimensions:',n1bis,n2bis,n3bis

      !pot_a = 0.d0

      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, v_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_c_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_f_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)
      call ocl_enqueue_write_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

      call nanosec(tsc0)
      do itimes=1,ntimes
         call kinetic_stable_d(queue,(/n1bis/2,n2bis/2,n3bis/2/),(/0.1d0,.1d0,.1d0/),&
            &   work_GPU,psi_GPU,work2_GPU,v_GPU,psi_c_GPU,psi_f_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, v_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
      call ocl_enqueue_read_buffer(queue, work2_GPU, n1bis*n2bis*n3bis*8, pot_a)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)
      call ocl_release_mem_object(work2_GPU)
      call ocl_release_mem_object(v_GPU)
      call ocl_release_mem_object(psi_c_GPU)
      call ocl_release_mem_object(psi_f_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9

      call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

      call compare_time('Kinetic 3D',(/n1bis,n2bis,n3bis/),&
           CPUtime,GPUtime,n1bis*n2bis*n3bis,3*45,ntimes,maxdiff,3.d-9)

   endif

   i_all=-product(shape(psi_k_in_a))
   deallocate(psi_k_in_a,stat=i_stat)
   call memocc(i_stat,i_all,'psi_k_in_a',subname)
   i_all=-product(shape(psi_cuda_k_in_a))
   deallocate(psi_cuda_k_in_a,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_k_in_a',subname)
   i_all=-product(shape(pot_a))
   deallocate(pot_a,stat=i_stat)
   call memocc(i_stat,i_all,'pot_a',subname)

   i_all=-product(shape(psi_k_out_a))
   deallocate(psi_k_out_a,stat=i_stat)
   call memocc(i_stat,i_all,'psi_k_out_a',subname)
   i_all=-product(shape(psi_cuda_k_out_a))
   deallocate(psi_cuda_k_out_a,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_k_out_a',subname)



!   write(*,'(a,i6,i6)')'CPU Convolutions shrink, dimensions:',n1-15,n2*n3

   allocate(psi_in_s(n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in_s,'psi_in_s',subname)
   allocate(v_cuda_s(n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_s,'v_cuda_s',subname)

   allocate(psi_in_t(n2*n3,n1,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in_t,'psi_in_t',subname)
   allocate(v_cuda_t(n1,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,v_cuda_t,'v_cuda_t',subname)

   allocate(psi_out_s(n2*n3,n1-15,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_out_s,'psi_out_s',subname)
   allocate(psi_cuda_s(n1-15,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_s,'psi_cuda_s',subname)

   allocate(psi_out_t(n1-15,n2*n3,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_out_t,'psi_out_t',subname)
   allocate(psi_cuda_t(n2*n3,n1-15,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_t,'psi_cuda_t',subname)

   do i=1,n2*n3
      do i1=1,n1
         v_cuda_s(i,i1,1)=psi_in(i1,i,1)
         psi_in_s(i1,i,1)=psi_in(i1,i,1)
      end do
   end do

   call nanosec(tsc0)
   do i=1,ntimes
      call convrot_shrink(n1-16,n2*n3,psi_in_s,psi_out_s)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!  write(*,'(a,i6,i6)')'GPU Convolutions shrink, dimensions:',n1-15,n2*n3


   call ocl_create_write_buffer(context, (n1-15)*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda_s)

   call nanosec(tsc0)
   do i=1,ntimes
      call magicfiltershrink1d_d(queue,n1-15,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, (n1-15)*n2*n3*8, psi_cuda_s)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results_t(n2*n3, n1-15, psi_out_s, psi_cuda_s, maxdiff, 1d-9)

   call compare_time('MF Shrink',(/n1-15,n2*n3/),CPUtime,GPUtime,(n1-15)*n2*n3,32,ntimes,maxdiff,1.d-9)

   do i=1,n2*n3
      do i1=1,n1-15
         psi_out_t(i1,i,1) = psi_out_s(i,i1,1)
         psi_cuda_t(i,i1,1) = psi_cuda_s(i1,i,1)
      end do
   end do

!   write(*,'(a,i6,i6)')'CPU Convolutions grow, dimensions:',n1-15,n2*n3

   call nanosec(tsc0)
   do i=1,ntimes
      call convrot_grow(n1-16,n2*n3,psi_out_t,psi_in_t)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!   write(*,'(a,i6,i6)')'GPU Convolutions grow, dimensions:',n1-15,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, (n1-15)*n2*n3*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, (n1-15)*n2*n3*8, psi_cuda_t)

   call nanosec(tsc0)
   do i=1,ntimes
      call magicfiltergrow1d_d(queue,n1-15,n2*n3,work_GPU,psi_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, v_cuda_t)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   call compare_2D_results_t(n2*n3, n1, psi_in_t, v_cuda_t, maxdiff, 1d-9)

   call compare_time('MF Grow',(/n1-15,n2*n3/),CPUtime,GPUtime,(n1-15)*n2*n3,32,ntimes,maxdiff,1.d-9)

   i_all=-product(shape(psi_out_s))
   deallocate(psi_out_s,stat=i_stat)
   call memocc(i_stat,i_all,'psi_out_s',subname)
   i_all=-product(shape(psi_cuda_s))
   deallocate(psi_cuda_s,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_s',subname)

   i_all=-product(shape(psi_in_s))
   deallocate(psi_in_s,stat=i_stat)
   call memocc(i_stat,i_all,'psi_in_s',subname)
   i_all=-product(shape(v_cuda_s))
   deallocate(v_cuda_s,stat=i_stat)
   call memocc(i_stat,i_all,'v_cuda_s',subname)

   i_all=-product(shape(psi_in_t))
   deallocate(psi_in_t,stat=i_stat)
   call memocc(i_stat,i_all,'psi_in_t',subname)
   i_all=-product(shape(v_cuda_t))
   deallocate(v_cuda_t,stat=i_stat)
   call memocc(i_stat,i_all,'v_cuda_t',subname)

   i_all=-product(shape(psi_out_t))
   deallocate(psi_out_t,stat=i_stat)
   call memocc(i_stat,i_all,'psi_out_t',subname)
   i_all=-product(shape(psi_cuda_t))
   deallocate(psi_cuda_t,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_t',subname)



   !**************************************************kinetic term
   n1bis = n1
   n2bis = n2
   n3bis = n3
!   write(*,'(a,i6,i6,i6)')'CPU Kinetic k 3D, dimensions:',n1bis,n2bis,n3bis

   allocate(psi_k_in(2,n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_k_in,'psi_k_in',subname)
   allocate(psi_cuda_k_in(2,n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_k_in,'psi_cuda_k_in',subname)
   allocate(psi_cuda_k_in_bis(2,n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda_k_in_bis,'psi_cuda_k_in_bis',subname)

   allocate(psi_k_out(2,n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_k_out,'psi_k_out',subname)
   allocate(psi_cuda_k_out(2,n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_k_out,'psi_cuda_k_out',subname)

   sigma2=0.25d0*((n1bis*hx)**2)
   do i3=1,n3bis
      do i2=1,n2bis
         do i1=1,n1bis
            x=hx*real(i1-n1bis/2-1,kind=8)
            !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
            r2=x**2
            arg=0.5d0*r2/sigma2
            tt=dexp(-arg)
            psi_k_in(1,i1,i2,i3,1)=tt
            psi_cuda_k_in(1,i1,i2,i3,1) = tt;
            arg=0.55d0*r2/sigma2
            tt=dexp(-arg)
            psi_k_in(2,i1,i2,i3,1)=tt
            psi_cuda_k_in(2,i1,i2,i3,1) = tt;
         end do
      end do
   end do

   call nanosec(tsc0)
   do itimes=1,ntimes
      call convolut_kinetic_per_c_k(n1bis-1,n2bis-1,n3bis-1,(/0.1d0,.1d0,.1d0/),&
         &   psi_k_in,psi_k_out,0.2d0,0.1d0,0.2d0,0.3d0)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,i6,i6,i6)')'GPU Kinetic k 3D, dimensions:',n1bis,n2bis,n3bis


   call ocl_create_read_write_buffer(context, 2*n1bis*n2bis*n3bis*8, psi_GPU)
   call ocl_create_read_write_buffer(context, 2*n1bis*n2bis*n3bis*8, work_GPU)
   call ocl_create_read_write_buffer(context, 2*n1bis*n2bis*n3bis*8, work2_GPU)
   call ocl_create_read_write_buffer(context, 2*n1bis*n2bis*n3bis*8, v_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, 2*n1bis*n2bis*n3bis*8, psi_cuda_k_in)

   call nanosec(tsc0)
   do itimes=1,ntimes
      call kinetic_k_d(queue,(/n1bis,n2bis,n3bis/),(/0.1d0,.1d0,.1d0/),&
         &   work_GPU,psi_GPU,work2_GPU,v_GPU,0.2d0,(/0.1d0,0.2d0,0.3d0/))
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, 2*n1bis*n2bis*n3bis*8, psi_cuda_k_out)
   call ocl_enqueue_read_buffer(queue, work_GPU, 2*n1bis*n2bis*n3bis*8, psi_cuda_k_in_bis)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)
   call ocl_release_mem_object(v_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9

   
   call compare_3D_results(n1bis*2, n2bis, n3bis, psi_k_out, psi_cuda_k_out, maxdiff, 1.d-9)

   call compare_time('Kinetic k 3D',(/n1bis,n2bis,n3bis/),&
        CPUtime,GPUtime,n1bis*n2bis*n3bis,4*3*45,ntimes,maxdiff,1.d-9)

   i_all=-product(shape(psi_k_in))
   deallocate(psi_k_in,stat=i_stat)
   call memocc(i_stat,i_all,'psi_k_in',subname)
   i_all=-product(shape(psi_cuda_k_in))
   deallocate(psi_cuda_k_in,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_k_in',subname)
   i_all=-product(shape(psi_cuda_k_in_bis))
   deallocate(psi_cuda_k_in_bis,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_k_in_bis',subname)

   i_all=-product(shape(psi_k_out))
   deallocate(psi_k_out,stat=i_stat)
   call memocc(i_stat,i_all,'psi_k_out',subname)
   i_all=-product(shape(psi_cuda_k_out))
   deallocate(psi_cuda_k_out,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_k_out',subname)




!   write(*,'(a,i6,i6)')'CPU Kinetic, dimensions:',n1,n2*n3

   allocate(modarr(lowfilK:n1-1+lupfilK+ndebug),stat=i_stat)
   call memocc(i_stat,modarr,'modarr',subname)
   call fill_mod_arr(modarr,lowfilK,n1-1+lupfilK,n1)


   psi_out=0.d0

   call nanosec(tsc0)
   do itimes=1,ntimes
      ekin=0.0_gp
      call conv_kin_x(psi_in,psi_out,n2*n3,ekin)   

   end do
   call nanosec(tsc1)

   i_all=-product(shape(modarr))
   deallocate(modarr,stat=i_stat)
   call memocc(i_stat,i_all,'modarr',subname)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,i6,i6)')'GPU Kinetic, dimensions:',n1,n2*n3

   call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
   call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
   call ocl_create_write_buffer(context, n1*n2*n3*8, work2_GPU)
   call ocl_create_read_write_buffer(context, n1*n2*n3*8, v_GPU)
   call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

   call nanosec(tsc0)
   do i=1,ntimes
      call kinetic1d_d(queue,n1,n2*n3,real(hx,kind=8),real(0.d0,kind=8),&
         &   work_GPU,psi_GPU,work2_GPU,v_GPU,ekinGPUd)
   end do
   call ocl_finish(queue)
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
   call ocl_release_mem_object(psi_GPU)
   call ocl_release_mem_object(work_GPU)
   call ocl_release_mem_object(work2_GPU)
   call ocl_release_mem_object(v_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_2D_results_t(n2*n3,n1,psi_out,psi_cuda, maxdiff, 1d-9)

   call compare_time('Kinetic',(/n1,n2*n3/),&
        CPUtime,GPUtime,n1*n2*n3,45,ntimes,maxdiff,1.d-9)

   !**************************************************wavelet transformations
   if (modulo(n1,2) == 0) then
      n1bis = n1
      n2bis = n2
      n3bis = n3
!      write(*,'(a,i6,i6,i6)')'CPU Analysis 3D, dimensions:',n1bis,n2bis,n3bis

      allocate(psi_k_in_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_k_in_a,'psi_k_in_a',subname)
      allocate(psi_cuda_k_in_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_cuda_k_in_a,'psi_cuda_k_in_a',subname)

      allocate(psi_k_out_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_k_out_a,'psi_k_out_a',subname)
      allocate(psi_cuda_k_out_a(n1bis,n2bis,n3bis,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_k_out_a,'psi_cuda_k_out_a',subname)

      sigma2=0.25d0*((n1bis*hx)**2)
      do i3=1,n3bis
         do i2=1,n2bis
            do i1=1,n1bis
               x=hx*real(i1-n1bis/2-1,kind=8)
               !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
               r2=x**2
               arg=0.5d0*r2/sigma2
               tt=dexp(-arg)
               psi_k_in_a(i1,i2,i3,1)=tt
               psi_cuda_k_in_a(i1,i2,i3,1) = tt;
            end do
         end do
      end do

      call nanosec(tsc0)
      do itimes=1,ntimes
         call analyse_per(n1bis/2-1,n2bis/2-1,n3bis/2-1,psi_cuda_k_out_a,psi_k_in_a,psi_k_out_a)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!      write(*,'(a,i6,i6,i6)')'GPU Analysis 3D, dimensions:',n1bis,n2bis,n3bis

      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

      call nanosec(tsc0)
      do itimes=1,ntimes
         call ana_d(queue,(/n1bis/2,n2bis/2,n3bis/2/), work2_GPU, work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)
      call ocl_release_mem_object(work2_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9


      call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

      call compare_time('Analysis 3D',(/n1bis,n2bis,n3bis/),&
           CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,1.d-9)

!      write(*,'(a,i6,i6,i6)')'GPU Analysis 3D (block), dimensions:',n1bis,n2bis,n3bis

      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

      call nanosec(tsc0)
      do itimes=1,ntimes
         call ana_block_d(queue,(/n1bis/2,n2bis/2,n3bis/2/), work2_GPU, work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)
      call ocl_release_mem_object(work2_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9


      call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

      call compare_time('Analysis 3D (B)',(/n1bis,n2bis,n3bis/),&
           CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,1.d-9)

      !write(*,'(a,i6,i6,i6)')'CPU Synthesis 3D, dimensions:',n1bis,n2bis,n3bis

      sigma2=0.25d0*((n1bis*hx)**2)
      do i3=1,n3bis
         do i2=1,n2bis
            do i1=1,n1bis
               x=hx*real(i1-n1bis/2-1,kind=8)
               !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
               r2=x**2
               arg=0.5d0*r2/sigma2
               tt=dexp(-arg)
               psi_k_in_a(i1,i2,i3,1)=tt
               psi_cuda_k_in_a(i1,i2,i3,1) = tt;
            end do
         end do
      end do

      call nanosec(tsc0)
      do itimes=1,ntimes
         call synthese_per(n1bis/2-1,n2bis/2-1,n3bis/2-1,psi_cuda_k_out_a,psi_k_in_a,psi_k_out_a)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!      write(*,'(a,i6,i6,i6)')'GPU Synthesis 3D, dimensions:',n1bis,n2bis,n3bis

      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, psi_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work_GPU)
      call ocl_create_read_write_buffer(context, n1bis*n2bis*n3bis*8, work2_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_in_a)

      call nanosec(tsc0)
      do itimes=1,ntimes
         call syn_d(queue,(/n1bis/2,n2bis/2,n3bis/2/), work2_GPU, work_GPU, psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1bis*n2bis*n3bis*8, psi_cuda_k_out_a)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)
      call ocl_release_mem_object(work2_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9


      call compare_3D_results(n1bis, n2bis, n3bis, psi_k_out_a, psi_cuda_k_out_a, maxdiff, 1d-9)

      call compare_time('Synthesis 3D',(/n1bis,n2bis,n3bis/),&
           CPUtime,GPUtime,n1bis*n2bis*n3bis,3*32,ntimes,maxdiff,1.d-9)

      i_all=-product(shape(psi_k_in_a))
      deallocate(psi_k_in_a,stat=i_stat)
      call memocc(i_stat,i_all,'psi_k_in_a',subname)
      i_all=-product(shape(psi_cuda_k_in_a))
      deallocate(psi_cuda_k_in_a,stat=i_stat)
      call memocc(i_stat,i_all,'psi_cuda_k_in_a',subname)

      i_all=-product(shape(psi_k_out_a))
      deallocate(psi_k_out_a,stat=i_stat)
      call memocc(i_stat,i_all,'psi_k_out_a',subname)
      i_all=-product(shape(psi_cuda_k_out_a))
      deallocate(psi_cuda_k_out_a,stat=i_stat)
      call memocc(i_stat,i_all,'psi_cuda_k_out_a',subname)

!      write(*,'(a,i6,i6)')'CPU Analysis, dimensions:',n1,n2*n3

      call nanosec(tsc0)
      do i=1,ntimes
         call ana_rot_per(n1/2-1,n2*n3,psi_in,psi_out)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9


 !     write(*,'(a,i6,i6)')'GPU Analysis, dimensions:',n1,n2*n3

      call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
      call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

      call nanosec(tsc0)
      do i=1,ntimes
         call ana1d_d(queue,n1/2,n2*n3,work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9


      call compare_2D_results_t(n2*n3, n1, psi_out, psi_cuda, maxdiff, 1d-9)

      call compare_time('Analysis',(/n1,n2*n3/),&
           CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,1.d-9)

!      write(*,'(a,i6,i6)')'GPU Analysis (block), dimensions:',n1,n2*n3

      call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
      call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

      call nanosec(tsc0)
      do i=1,ntimes
         call ana1d_block_d(queue,n1/2,n2*n3,work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9

      call compare_2D_results_t(n2*n3, n1, psi_out, psi_cuda, maxdiff, 1d-9)

      call compare_time('Analysis (B)',(/n1,n2*n3/),&
           CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,1.d-9)

      !write(*,'(a,i6,i6)')'CPU Synthesis, dimensions:',n1,n2*n3
 
      call nanosec(tsc0)
      do i=1,ntimes*100
         !call syn_rot_per(n1/2-1,n2*n3,psi_in,psi_out)
         !call syn_rot_per_temp(n1/2,n2*n3,psi_in,psi_out)
         !call syn_rot_per_simple(n1/2-1,n2*n3,psi_in,psi_out)
         call synthesis_per_u5(n1/2,n2*n3,psi_in,psi_out)
         !call synthesis_per_u2(n1/2,n2*n3,psi_in,psi_out)
         !call synthesis_per_u12(n1/2,n2*n3,psi_in,psi_out)
         !call synthesis_per_u24(n1/2,n2*n3,psi_in,psi_out)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9

      !write(*,'(a,i6,i6)')'GPU Synthesis, dimensions:',n1,n2*n3

      call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
      call ocl_create_read_buffer(context, n1*n2*n3*8, work_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, n1*n2*n3*8, v_cuda)

      call nanosec(tsc0)
      do i=1,ntimes
         call syn1d_d(queue,n1/2,n2*n3,work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, psi_cuda)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9


      call compare_2D_results_t(n2*n3, n1, psi_out, psi_cuda, maxdiff, 1d-9)

      call compare_time('Synthesis',(/n1,n2*n3/),&
           CPUtime,GPUtime,n1*n2*n3,32,ntimes,maxdiff,1.d-9)

!      write(*,'(a,i6,i6)')'CPU Analysis shrink, dimensions:',n1-14,n2*n3

      allocate(psi_in_s(n1,n2*n3,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_in_s,'psi_in_s',subname)
      allocate(v_cuda_s(n2*n3,n1,1+ndebug),stat=i_stat)
      call memocc(i_stat,v_cuda_s,'v_cuda_s',subname)

      allocate(psi_in_t(n2*n3,n1,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_in_t,'psi_in_t',subname)
      allocate(v_cuda_t(n1,n2*n3,1+ndebug),stat=i_stat)
      call memocc(i_stat,v_cuda_t,'v_cuda_t',subname)

      allocate(psi_out_s(n2*n3,n1-14,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_out_s,'psi_out_s',subname)
      allocate(psi_cuda_s(n1-14,n2*n3,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_cuda_s,'psi_cuda_s',subname)

      allocate(psi_out_t(n1-14,n2*n3,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_out_t,'psi_out_t',subname)
      allocate(psi_cuda_t(n2*n3,n1-14,1+ndebug),stat=i_stat)
      call memocc(i_stat,psi_cuda_t,'psi_cuda_t',subname)

      do i=1,n2*n3
         do i1=1,n1
            v_cuda_s(i,i1,1)=psi_in(i1,i,1)
            psi_in_s(i1,i,1)=psi_in(i1,i,1)
         end do
      end do

      call nanosec(tsc0)
      do i=1,ntimes
         call ana_rot_shrink((n1-14)/2-1,n2*n3,psi_in_s,psi_out_s)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9

!      write(*,'(a,i6,i6)')'GPU Analysis shrink, dimensions:',n1-14,n2*n3

      call ocl_create_write_buffer(context, (n1-14)*n2*n3*8, psi_GPU)
      call ocl_create_read_buffer(context, (n1)*n2*n3*8, work_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, (n1)*n2*n3*8, v_cuda_s)

      call nanosec(tsc0)
      do i=1,ntimes
         call anashrink1d_d(queue,(n1-14)/2,n2*n3,work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, (n1-14)*n2*n3*8, psi_cuda_s)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9

      call compare_2D_results_t(n2*n3, n1-14, psi_out_s, psi_cuda_s, maxdiff, 1d-9)

      call compare_time('Analysis shrink',(/n1-14,n2*n3/),&
           CPUtime,GPUtime,(n1-14)*n2*n3,32,ntimes,maxdiff,1.d-9)

      do i=1,n2*n3
         do i1=1,(n1-14)
            psi_out_t(i1,i,1) = psi_out_s(i,i1,1)
            psi_cuda_t(i,i1,1) = psi_cuda_s(i1,i,1)
         end do
      end do

!      write(*,'(a,i6,i6)')'CPU Synthesis grow, dimensions:',n1-14,n2*n3
      call nanosec(tsc0)
      do i=1,ntimes*100
         !call syn_rot_grow((n1-14)/2-1,n2*n3,psi_out_t,psi_in_t)
         !call synthesis_free((n1-14)/2,n2*n3,psi_out_t,psi_in_t)
         call synthesis_free_u4((n1-14)/2,n2*n3,psi_out_t,psi_in_t)
      end do
      call nanosec(tsc1)

      CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!      write(*,'(a,i6,i6)')'GPU Synthesis grow, dimensions:',n1-14,n2*n3

      call ocl_create_write_buffer(context, n1*n2*n3*8, psi_GPU)
      call ocl_create_read_buffer(context, (n1-14)*n2*n3*8, work_GPU)
      call ocl_enqueue_write_buffer(queue, work_GPU, (n1-14)*n2*n3*8, psi_cuda_t)

      call nanosec(tsc0)
      do i=1,ntimes
         call syngrow1d_d(queue,(n1-14)/2,n2*n3,work_GPU,psi_GPU)
      end do
      call ocl_finish(queue);
      call nanosec(tsc1)

      call ocl_enqueue_read_buffer(queue, psi_GPU, n1*n2*n3*8, v_cuda_t)
      call ocl_release_mem_object(psi_GPU)
      call ocl_release_mem_object(work_GPU)

      GPUtime=real(tsc1-tsc0,kind=8)*1d-9

      call compare_2D_results_t(n2*n3, n1, psi_in_t, v_cuda_t, maxdiff, 1d-9)

      call compare_time('Synthesis grow',(/n1-14,n2*n3/),&
           CPUtime,GPUtime,(n1-14)*n2*n3,32,ntimes,maxdiff,1.d-9)

      i_all=-product(shape(psi_out_s))
      deallocate(psi_out_s,stat=i_stat)
      call memocc(i_stat,i_all,'psi_out_s',subname)
      i_all=-product(shape(psi_cuda_s))
      deallocate(psi_cuda_s,stat=i_stat)
      call memocc(i_stat,i_all,'psi_cuda_s',subname)

      i_all=-product(shape(psi_in_s))
      deallocate(psi_in_s,stat=i_stat)
      call memocc(i_stat,i_all,'psi_in_s',subname)
      i_all=-product(shape(v_cuda_s))
      deallocate(v_cuda_s,stat=i_stat)
      call memocc(i_stat,i_all,'v_cuda_s',subname)

      i_all=-product(shape(psi_in_t))
      deallocate(psi_in_t,stat=i_stat)
      call memocc(i_stat,i_all,'psi_in_t',subname)
      i_all=-product(shape(v_cuda_t))
      deallocate(v_cuda_t,stat=i_stat)
      call memocc(i_stat,i_all,'v_cuda_t',subname)

      i_all=-product(shape(psi_out_t))
      deallocate(psi_out_t,stat=i_stat)
      call memocc(i_stat,i_all,'psi_out_t',subname)
      i_all=-product(shape(psi_cuda_t))
      deallocate(psi_cuda_t,stat=i_stat)
      call memocc(i_stat,i_all,'psi_cuda_t',subname)



   end if

   i_all=-product(shape(psi_out))
   deallocate(psi_out,stat=i_stat)
   call memocc(i_stat,i_all,'psi_out',subname)
   i_all=-product(shape(psi_in))
   deallocate(psi_in,stat=i_stat)
   call memocc(i_stat,i_all,'psi_in',subname)
   i_all=-product(shape(psi_cuda))
   deallocate(psi_cuda,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda',subname)
   i_all=-product(shape(v_cuda))
   deallocate(v_cuda,stat=i_stat)
   call memocc(i_stat,i_all,'v_cuda',subname)
   i_all=-product(shape(psi_cuda_l))
   deallocate(psi_cuda_l,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda_l',subname)
   i_all=-product(shape(v_cuda_l))
   deallocate(v_cuda_l,stat=i_stat)
   call memocc(i_stat,i_all,'v_cuda_l',subname)

   !**************************************************compression-decompression
   !create keys arrays
   !cubic domain
   nseg=(n1+1)*(n1+1)

!   print *,'nseg=',nseg

   allocate(keyg(2,nseg+ndebug),stat=i_stat)
   call memocc(i_stat,keyg,'keyg',subname)
   allocate(keyv(nseg+ndebug),stat=i_stat)
   call memocc(i_stat,keyv,'keyv',subname)

   !take a rectangle of a cubic region
   i1s=5
   i1e=n1-5
   keyv(1)=1
   do i3=0,n1
      do i2=0,n1
         iseg=i2+1+(n1+1)*i3
         keyg(1,iseg)=i3*((n1+1)*(n1+1)) + i2*(n1+1) + i1s+1
         keyg(2,iseg)=i3*((n1+1)*(n1+1)) + i2*(n1+1) + i1e+1
         if (iseg >1) keyv(iseg)=keyv(iseg-1)+i1e-i1s+1
      end do
   end do
   nvctr_cf=keyv(nseg)+i1e-i1s

   allocate(psi(8*nvctr_cf+ndebug),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)
   allocate(psi_l(8*nvctr_cf+ndebug),stat=i_stat)
   call memocc(i_stat,psi_l,'psi_l',subname)
   allocate(psi_d(8*nvctr_cf+ndebug),stat=i_stat)
   call memocc(i_stat,psi_d,'psi_d',subname)
   !determine the values for psi function
   do i=1,nvctr_cf
      psi(i)=real(i,kind=8)
      psi_l(i)=real(i,kind=4)
   end do
   do i=1,nvctr_cf
      do j=1,7
         psi(nvctr_cf+7*(i-1)+j)=real(i,kind=8)+0.1d0*real(j,kind=8)
         psi_l(nvctr_cf+7*(i-1)+j)=real(i,kind=4)+0.10*real(j,kind=4)
      end do
   end do

   allocate(psi_in((2*n1+2),(2*n1+2),(2*n1+2)+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in,'psi_in',subname)
   allocate(psi_out((2*n1+2),(2*n1+2),(2*n1+2)+ndebug),stat=i_stat)
   call memocc(i_stat,psi_in,'psi_out',subname)

   allocate(psi_cuda((2*n1+2),(2*n1+2),(2*n1+2)+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda,'psi_cuda',subname)

   allocate(psi_cuda_l((2*n1+2),(2*n1+2),(2*n1+2)+ndebug),stat=i_stat)
   call memocc(i_stat,psi_cuda,'psi_cuda_l',subname)


!  write(*,'(a,3(i6))')'CPU Uncompress, dimensions:',n1,n1,n1

   !take timings
   call nanosec(tsc0)
   do i=1,ntimes
      call uncompress(n1,n1,n1,nseg,nvctr_cf,keyg,keyv,  & 
      nseg,nvctr_cf,keyg,keyv,psi(1),psi(nvctr_cf+1),psi_in)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!  write(*,'(a,3(i6))')'GPU Uncompress, dimensions:',n1,n1,n1

   call ocl_create_read_buffer(context, nvctr_cf*8, psi_c_GPU)
   call ocl_create_read_buffer(context, 7*nvctr_cf*8, psi_f_GPU)
   call ocl_create_read_buffer(context, nseg*4*2, keyg_GPU)
   call ocl_create_read_buffer(context, nseg*4, keyv_GPU)
   call ocl_create_write_buffer(context, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, psi_c_GPU, nvctr_cf*8, psi)
   call ocl_enqueue_write_buffer(queue, psi_f_GPU, 7*nvctr_cf*8, psi(nvctr_cf+1))
   call ocl_enqueue_write_buffer(queue, keyg_GPU, nseg*2*4, keyg)
   call ocl_enqueue_write_buffer(queue, keyv_GPU, nseg*4, keyv)

   call nanosec(tsc0)
   do i=1,ntimes
      call uncompress_d(queue , (/n1+1,n1+1,n1+1/),&
         &   nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
      nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
         &   psi_c_GPU, psi_f_GPU, work_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, work_GPU, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, psi_cuda)
   call ocl_release_mem_object(psi_c_GPU)
   call ocl_release_mem_object(psi_f_GPU)
   call ocl_release_mem_object(keyg_GPU)
   call ocl_release_mem_object(keyv_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_3D_results(n1*2+2, n1*2+2, n1*2+2, psi_in, psi_cuda, maxdiff,1d-9)

   call compare_time('Uncompress',(/n1,n1,n1/),CPUtime,GPUtime,8*nvctr_cf,1,ntimes,maxdiff,1.d-9)

!   write(*,'(a,3(i6))')'CPU Compress, dimensions:',n1,n1,n1

   call nanosec(tsc0)
   do i=1,ntimes
      call compress(n1,n1,0,n1,0,n1,0,n1,nseg,nvctr_cf,keyg,keyv,  & 
      nseg,nvctr_cf,keyg,keyv,psi_in,psi(1),psi(nvctr_cf+1))
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,3(i6))')'GPU Compress, dimensions:',n1,n1,n1

   call ocl_create_write_buffer(context, nvctr_cf*8, psi_c_GPU)
   call ocl_create_write_buffer(context, 7*nvctr_cf*8, psi_f_GPU)
   call ocl_create_read_buffer(context, nseg*8*2, keyg_GPU)
   call ocl_create_read_buffer(context, nseg*8, keyv_GPU)
   call ocl_create_read_buffer(context, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, keyg_GPU, nseg*2*4, keyg)
   call ocl_enqueue_write_buffer(queue, keyv_GPU, nseg*4, keyv)
   call ocl_enqueue_write_buffer(queue, work_GPU, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, psi_cuda)

   call nanosec(tsc0)
   do i=1,ntimes
      call compress_d(queue , (/n1+1, n1+1, n1+1/),&
         &   nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
      nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
         &   psi_c_GPU, psi_f_GPU, work_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_c_GPU, nvctr_cf*8, psi_d)
   call ocl_enqueue_read_buffer(queue, psi_f_GPU, 7*nvctr_cf*8, psi_d(nvctr_cf+1))
   call ocl_release_mem_object(psi_c_GPU)
   call ocl_release_mem_object(psi_f_GPU)
   call ocl_release_mem_object(keyg_GPU)
   call ocl_release_mem_object(keyv_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_1D_results(8*nvctr_cf, psi, psi_d, maxdiff,1d-9)

   call compare_time('Compress',(/n1,n1,n1/),CPUtime,GPUtime,8*nvctr_cf,1,ntimes,maxdiff,1.d-9)

!   write(*,'(a,3(i6))')'CPU Uncompress Scal, dimensions:',n1,n1,n1

   call wscal_init_per(scal,0.1d0,0.2d0,0.3d0,0.4d0)
   call nanosec(tsc0)
   do i=1,ntimes
      call uncompress_scal(n1,n1,n1,nseg,nvctr_cf,keyg,keyv,  & 
      nseg,nvctr_cf,keyg,keyv,psi(1),psi(nvctr_cf+1),psi_in,scal)
   end do
   call nanosec(tsc1)
   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,3(i6))')'GPU Uncompress Scal, dimensions:',n1,n1,n1

   call ocl_create_read_buffer(context, nvctr_cf*8, psi_c_GPU)
   call ocl_create_read_buffer(context, 7*nvctr_cf*8, psi_f_GPU)
   call ocl_create_read_buffer(context, nseg*4*2, keyg_GPU)
   call ocl_create_read_buffer(context, nseg*4, keyv_GPU)
   call ocl_create_write_buffer(context, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, psi_c_GPU, nvctr_cf*8, psi)
   call ocl_enqueue_write_buffer(queue, psi_f_GPU, 7*nvctr_cf*8, psi(nvctr_cf+1))
   call ocl_enqueue_write_buffer(queue, keyg_GPU, nseg*2*4, keyg)
   call ocl_enqueue_write_buffer(queue, keyv_GPU, nseg*4, keyv)

   call nanosec(tsc0)
   do i=1,ntimes
      call uncompress_scale_d(queue , (/n1+1,n1+1,n1+1/),(/0.1d0/2.0d0,0.2d0/2.0d0,0.3d0/2.0d0/),0.4d0,&
         &   nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
      nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
         &   psi_c_GPU, psi_f_GPU, work_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, work_GPU, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, psi_cuda)
   call ocl_release_mem_object(psi_c_GPU)
   call ocl_release_mem_object(psi_f_GPU)
   call ocl_release_mem_object(keyg_GPU)
   call ocl_release_mem_object(keyv_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_3D_results(n1*2+2, n1*2+2, n1*2+2, psi_in, psi_cuda, maxdiff,1d-9)

   call compare_time('Uncompress Scal',(/n1,n1,n1/),CPUtime,GPUtime,8*nvctr_cf,1,ntimes,maxdiff,1.d-9)

!   write(*,'(a,3(i6))')'CPU Compress Scal, dimensions:',n1,n1,n1

   call nanosec(tsc0)
   do i=1,ntimes
      call compress_scal(n1,n1,n1,nseg,nvctr_cf,keyg,keyv,  & 
      nseg,nvctr_cf,keyg,keyv,psi_in,psi(1),psi(nvctr_cf+1),scal)
   end do
   call nanosec(tsc1)

   CPUtime=real(tsc1-tsc0,kind=8)*1d-9


!   write(*,'(a,3(i6))')'GPU Compress Scal, dimensions:',n1,n1,n1

   call ocl_create_write_buffer(context, nvctr_cf*8, psi_c_GPU)
   call ocl_create_write_buffer(context, 7*nvctr_cf*8, psi_f_GPU)
   call ocl_create_read_buffer(context, nseg*8*2, keyg_GPU)
   call ocl_create_read_buffer(context, nseg*8, keyv_GPU)
   call ocl_create_read_buffer(context, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, work_GPU)
   call ocl_enqueue_write_buffer(queue, keyg_GPU, nseg*2*4, keyg)
   call ocl_enqueue_write_buffer(queue, keyv_GPU, nseg*4, keyv)
   call ocl_enqueue_write_buffer(queue, work_GPU, (2*n1+2)*(2*n1+2)*(2*n1+2)*8, psi_cuda)

   call nanosec(tsc0)
   do i=1,ntimes
      call compress_scale_d(queue , (/n1+1, n1+1, n1+1/),(/0.1d0/2.0d0,0.2d0/2.0d0,0.3d0/2.0d0/),0.4d0,&
         &   nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
      nseg, nvctr_cf, keyg_GPU, keyv_GPU,&
         &   psi_c_GPU, psi_f_GPU, work_GPU)
   end do
   call ocl_finish(queue);
   call nanosec(tsc1)

   call ocl_enqueue_read_buffer(queue, psi_c_GPU, nvctr_cf*8, psi_d)
   call ocl_enqueue_read_buffer(queue, psi_f_GPU, 7*nvctr_cf*8, psi_d(nvctr_cf+1))
   call ocl_release_mem_object(psi_c_GPU)
   call ocl_release_mem_object(psi_f_GPU)
   call ocl_release_mem_object(keyg_GPU)
   call ocl_release_mem_object(keyv_GPU)
   call ocl_release_mem_object(work_GPU)

   GPUtime=real(tsc1-tsc0,kind=8)*1d-9


   call compare_1D_results(8*nvctr_cf, psi, psi_d, maxdiff,1d-9)

   call compare_time('Compress Scal',(/n1,n1,n1/),CPUtime,GPUtime,8*nvctr_cf,1,ntimes,maxdiff,1.d-9)

   i_all=-product(shape(psi))
   deallocate(psi,stat=i_stat)
   call memocc(i_stat,i_all,'psi',subname)
   i_all=-product(shape(psi_in))
   deallocate(psi_in,stat=i_stat)
   call memocc(i_stat,i_all,'psi_in',subname)
   i_all=-product(shape(psi_cuda))
   deallocate(psi_cuda,stat=i_stat)
   call memocc(i_stat,i_all,'psi_cuda',subname)
   i_all=-product(shape(keyg))
   deallocate(keyg,stat=i_stat)
   call memocc(i_stat,i_all,'keyg',subname)
   i_all=-product(shape(keyv))
   deallocate(keyv,stat=i_stat)
   call memocc(i_stat,i_all,'keyv',subname)

   !Need in order to have the last lines (TD,2011-11-10)
   call flush(6)

   call print_event_list
   call ocl_clean_command_queue(queue)
   call ocl_clean(context)

    !release the yaml document
    call yaml_release_document()


   contains

   subroutine print_time(field,dims,time,nbelem,nop,ntimes)
     use yaml_output
      implicit none
      character(len=*), intent(in) :: field
      integer, dimension(:), intent(in) :: dims
      real(gp),intent(in)::time
      integer,intent(in)::nbelem,nop,ntimes
      logical, parameter :: debug=.true.

      if (debug) then
!         write(*,'(a,f9.4,f10.2)')'Finished. Time(ms), GFlops',&
!            &   time*1.d3/real(ntimes,kind=8),&
!         real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9)
         call yaml_open_map(trim(field))
         call yaml_map('ms',time*1.d3/real(ntimes,kind=8),fmt='(f10.2)')
         call yaml_map('Gflop/s',real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9),fmt='(f10.2)')
         call yaml_map('Dims',dims,fmt='(i8)')
         call yaml_close_map()

      end if

   END SUBROUTINE print_time

   subroutine compare_time(field,dims,REFtime,TESTtime,nbelem,nop,ntimes,maxdiff,threshold)
     use yaml_output
      implicit none
      character(len=*), intent(in) :: field
      integer, dimension(:), intent(in) :: dims
      real(gp),intent(in)::REFtime,TESTtime,maxdiff,threshold
      integer,intent(in)::nbelem,nop,ntimes

!!$      write(*,'(1x,a)')'| CPU: ms  |  Gflops  || GPU:  ms |  GFlops  || Ratio  | No. Elements | Max. Diff. |'
!!$
!!$      write(*,'(1x,2(2(a,f10.2),a),a,f8.3,a,i14,a,1pe12.4,a)',advance='no')&
!!$         &   '|',REFtime*1.d3/real(ntimes,kind=8),'|',&
!!$         & real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(REFtime*1.d9),'|',&
!!$         &   '|',TESTtime*1.d3/real(ntimes,kind=8),'|',&
!!$         & real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(TESTtime*1.d9),'|',&
!!$         &   '|',REFtime/TESTtime,'|',&
!!$         & nbelem,'|',&
!!$         &   maxdiff,'|'
!!$      if (maxdiff <= threshold) then
!!$         write(*,'(a)')''
!!$      else
!!$         write(*,'(a)')'<<<< WARNING' 
!!$      end if

      !yaml output
      call yaml_comment(trim(field),hfill='-')
      call yaml_open_map(trim(field),flow=.true.)
!      call yaml_newline()
      call yaml_open_map('CPU')
       call yaml_map('ms',REFtime*1.d3/real(ntimes,kind=8),fmt='(f10.2)')
       call yaml_map('Gflop/s',real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(REFtime*1.d9),fmt='(f10.2)')
      call yaml_close_map()
      call yaml_open_map('GPU')
       call yaml_map('ms',TESTtime*1.d3/real(ntimes,kind=8),fmt='(f10.2)')
       call yaml_map('Gflop/s',real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(TESTtime*1.d9),fmt='(f10.2)')
      call yaml_close_map()
      call yaml_newline()
      call yaml_map('Ratio',REFtime/TESTtime,fmt='(f8.3)')
      call yaml_map('Dims',dims,fmt='(i8)')
      call yaml_map('No. Elems',nbelem,fmt='(i14)')
      call yaml_map('Max. Diff.',maxdiff,fmt='(1pe12.4)')
      call yaml_close_map(advance='no')
      if (maxdiff <= threshold) then
         call yaml_newline()   
      else
         call yaml_comment('<<<< WARNING')
      end if


   END SUBROUTINE compare_time

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
   END SUBROUTINE compare_3D_results

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
   END SUBROUTINE compare_2D_results

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
   END SUBROUTINE compare_2D_results_t

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
   END SUBROUTINE compare_1D_results

   subroutine conv_kin_x(x,y,ndat,ekin)
      implicit none
      integer,intent(in)::ndat
      real(wp),intent(in):: x(0:n1-1,ndat)
      real(wp),intent(out)::y(ndat,0:n1-1)
      real(wp),intent(inout)::ekin
      real(wp) tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12

      !$omp do
      do i=0,ndat/12-1
         do i1=0,n1-1
            tt1=0.e0_wp
            tt2=0.e0_wp
            tt3=0.e0_wp
            tt4=0.e0_wp
            tt5=0.e0_wp
            tt6=0.e0_wp
            tt7=0.e0_wp
            tt8=0.e0_wp
            tt9 =0.e0_wp
            tt10=0.e0_wp
            tt11=0.e0_wp
            tt12=0.e0_wp

            do l=lowfilK,lupfilK
               j=modarr(i1+l)

               tt1=tt1+x(j,i*12+1)*fil(l)
               tt2=tt2+x(j,i*12+2)*fil(l)
               tt3=tt3+x(j,i*12+3)*fil(l)
               tt4=tt4+x(j,i*12+4)*fil(l)
               tt5=tt5+x(j,i*12+5)*fil(l)
               tt6=tt6+x(j,i*12+6)*fil(l)
               tt7=tt7+x(j,i*12+7)*fil(l)
               tt8=tt8+x(j,i*12+8)*fil(l)
               tt9 =tt9 +x(j,i*12+9 )*fil(l)
               tt10=tt10+x(j,i*12+10)*fil(l)
               tt11=tt11+x(j,i*12+11)*fil(l)
               tt12=tt12+x(j,i*12+12)*fil(l)
            enddo
            y(i*12+1 ,i1)=tt1;     ekin=ekin+tt1*x(i1,i*12+1)
            y(i*12+2 ,i1)=tt2;     ekin=ekin+tt2*x(i1,i*12+2)
            y(i*12+3 ,i1)=tt3;     ekin=ekin+tt3*x(i1,i*12+3)
            y(i*12+4 ,i1)=tt4;     ekin=ekin+tt4*x(i1,i*12+4)
            y(i*12+5 ,i1)=tt5;     ekin=ekin+tt5*x(i1,i*12+5)
            y(i*12+6 ,i1)=tt6;     ekin=ekin+tt6*x(i1,i*12+6)
            y(i*12+7 ,i1)=tt7;     ekin=ekin+tt7*x(i1,i*12+7)
            y(i*12+8 ,i1)=tt8;     ekin=ekin+tt8*x(i1,i*12+8)
            y(i*12+9 ,i1)=tt9 ;    ekin=ekin+tt9 *x(i1,i*12+9 )
            y(i*12+10,i1)=tt10;    ekin=ekin+tt10*x(i1,i*12+10)
            y(i*12+11,i1)=tt11;    ekin=ekin+tt11*x(i1,i*12+11)
            y(i*12+12,i1)=tt12;    ekin=ekin+tt12*x(i1,i*12+12)
         enddo
      enddo
      !$omp end do

      !$omp do
      do i=(ndat/12)*12+1,ndat
         do i1=0,n1-1
            tt=0.e0_wp
            do l=lowfilK,lupfilK
               j=modarr(i1+l)
               tt=tt+x(j   ,i)*fil(l)
            enddo
            y(i,i1)=tt ; ekin=ekin+tt*x(i1,i)
         enddo
      enddo
      !$omp end do
   END SUBROUTINE conv_kin_x

include 'syn.txt'

END PROGRAM conv_check_ocl
