!> @file
!!    Test of convolution in GPU (conv_check)
!! @author
!!    Copyright (C) 2008-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Program test for the convolution in GPU
program conv_check
  use module_base
  implicit none
  integer  :: n1,n2,n3,n1bis,n2bis,n3bis
  integer :: iproc,nproc,ierr
  real(gp) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(wp), dimension(:,:,:), allocatable :: pot,psir,psi_in,psi_out,psi_out_s
  real(wp), dimension(:,:,:,:), allocatable :: psi_3d_in,psi_3d_out
  real(wp), dimension(:,:,:,:,:), allocatable :: psi_k_in, psi_k_out
  real(wp), dimension(:,:,:,:), allocatable :: psi_k_in_a, psi_k_out_a, pot_a
  real(wp), dimension(:,:,:), allocatable :: psi_in_s,psi_out_t,psi_in_t,psi_gemm,psi_cuda_gemm
  !local variables
  character(len=*), parameter :: subname='conv_check'
  character(len=50) :: chain
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat,i1_max,i_max,it0,it1,ndim,itimes
  integer :: count_rate,count_max,l,ierror,i1s,i1e
  integer :: n1s,n1e,ndats,ndate,nvctr_cf,nseg,iseg
  real(wp) :: tt,scale
  real(gp) :: v,p,CPUtime,GPUtime,comp,ekin
  real(gp), dimension(3) :: hgridh
  real(gp), dimension(8) :: scal
  integer, dimension(:), allocatable :: keyv,modarr
  integer, dimension(:,:), allocatable :: keyg
  real(kind=8), dimension(:), allocatable :: psi, psi_d !<temporary in view of wp 
  real(kind=4), dimension(:), allocatable :: psi_l      !<temporary in view of wp 
  real(kind=8), dimension(:,:,:), allocatable :: psi_cuda,v_cuda,psi_cuda_s,v_cuda_s,psi_cuda_t,v_cuda_t,v_cuda_str !<temporary in view of wp 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_3d_cuda, psi_3d_tmp 
  real(kind=8), dimension(:,:,:,:,:), allocatable :: psi_cuda_k_in,psi_cuda_k_out,psi_cuda_k_in_bis 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_k_in_a,psi_cuda_k_out_a 
  real(kind=4), dimension(:,:,:), allocatable :: psi_cuda_l,v_cuda_l !<temporary in view of wp 
  real(kind=8) :: ekinGPUd
  real(kind=4) :: t0,t1,epotGPU,ekinGPU
  real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU,keys_GPU !<pointer to the GPU  memory addresses (with norb=1)
  real(kind=8) :: psi_c_GPU, psi_f_GPU, keyg_GPU, keyv_GPU
  real(kind=8) :: context,queue
  integer, parameter :: lowfil1=-8,lupfil1=7   !<for GPU computation
  integer, parameter :: lowfil2=-7,lupfil2=8   !<for GPU computation
  integer, parameter :: lowfilK=-14,lupfilK=14 !< kinetic term
  real(kind=8), dimension(lowfilK:lupfilK) :: fil
  integer(kind=8) :: tsc0, tsc1

      INTEGER, PARAMETER :: PAPI_L1_DCM        = ((-2147483647) - 1)
      INTEGER, PARAMETER :: PAPI_L1_ICM        = -2147483647
      INTEGER, PARAMETER :: PAPI_L2_DCM        = -2147483646
      INTEGER, PARAMETER :: PAPI_L2_ICM        = -2147483645
      INTEGER, PARAMETER :: PAPI_L3_DCM        = -2147483644
      INTEGER, PARAMETER :: PAPI_L3_ICM        = -2147483643
      INTEGER, PARAMETER :: PAPI_L1_TCM        = -2147483642
      INTEGER, PARAMETER :: PAPI_L2_TCM        = -2147483641
      INTEGER, PARAMETER :: PAPI_L3_TCM        = -2147483640
      INTEGER, PARAMETER :: PAPI_CA_SNP        = -2147483639
      INTEGER, PARAMETER :: PAPI_CA_SHR        = -2147483638
      INTEGER, PARAMETER :: PAPI_CA_CLN        = -2147483637
      INTEGER, PARAMETER :: PAPI_CA_INV        = -2147483636
      INTEGER, PARAMETER :: PAPI_CA_ITV        = -2147483635
      INTEGER, PARAMETER :: PAPI_L3_LDM        = -2147483634
      INTEGER, PARAMETER :: PAPI_L3_STM        = -2147483633
      INTEGER, PARAMETER :: PAPI_BRU_IDL       = -2147483632
      INTEGER, PARAMETER :: PAPI_FXU_IDL       = -2147483631
      INTEGER, PARAMETER :: PAPI_FPU_IDL       = -2147483630
      INTEGER, PARAMETER :: PAPI_LSU_IDL       = -2147483629
      INTEGER, PARAMETER :: PAPI_TLB_DM        = -2147483628
      INTEGER, PARAMETER :: PAPI_TLB_IM        = -2147483627
      INTEGER, PARAMETER :: PAPI_TLB_TL        = -2147483626
      INTEGER, PARAMETER :: PAPI_L1_LDM        = -2147483625
      INTEGER, PARAMETER :: PAPI_L1_STM        = -2147483624
      INTEGER, PARAMETER :: PAPI_L2_LDM        = -2147483623
      INTEGER, PARAMETER :: PAPI_L2_STM        = -2147483622
      INTEGER, PARAMETER :: PAPI_BTAC_M        = -2147483621
      INTEGER, PARAMETER :: PAPI_PRF_DM        = -2147483620
      INTEGER, PARAMETER :: PAPI_L3_DCH        = -2147483619
      INTEGER, PARAMETER :: PAPI_TLB_SD        = -2147483618
      INTEGER, PARAMETER :: PAPI_CSR_FAL       = -2147483617
      INTEGER, PARAMETER :: PAPI_CSR_SUC       = -2147483616
      INTEGER, PARAMETER :: PAPI_CSR_TOT       = -2147483615
      INTEGER, PARAMETER :: PAPI_MEM_SCY       = -2147483614
      INTEGER, PARAMETER :: PAPI_MEM_RCY       = -2147483613
      INTEGER, PARAMETER :: PAPI_MEM_WCY       = -2147483612
      INTEGER, PARAMETER :: PAPI_STL_ICY       = -2147483611
      INTEGER, PARAMETER :: PAPI_FUL_ICY       = -2147483610
      INTEGER, PARAMETER :: PAPI_STL_CCY       = -2147483609
      INTEGER, PARAMETER :: PAPI_FUL_CCY       = -2147483608
      INTEGER, PARAMETER :: PAPI_HW_INT        = -2147483607
      INTEGER, PARAMETER :: PAPI_BR_UCN        = -2147483606
      INTEGER, PARAMETER :: PAPI_BR_CN         = -2147483605
      INTEGER, PARAMETER :: PAPI_BR_TKN        = -2147483604
      INTEGER, PARAMETER :: PAPI_BR_NTK        = -2147483603
      INTEGER, PARAMETER :: PAPI_BR_MSP        = -2147483602
      INTEGER, PARAMETER :: PAPI_BR_PRC        = -2147483601
      INTEGER, PARAMETER :: PAPI_FMA_INS       = -2147483600
      INTEGER, PARAMETER :: PAPI_TOT_IIS       = -2147483599
      INTEGER, PARAMETER :: PAPI_TOT_INS       = -2147483598
      INTEGER, PARAMETER :: PAPI_INT_INS       = -2147483597
      INTEGER, PARAMETER :: PAPI_FP_INS        = -2147483596
      INTEGER, PARAMETER :: PAPI_LD_INS        = -2147483595
      INTEGER, PARAMETER :: PAPI_SR_INS        = -2147483594
      INTEGER, PARAMETER :: PAPI_BR_INS        = -2147483593
      INTEGER, PARAMETER :: PAPI_VEC_INS       = -2147483592
      INTEGER, PARAMETER :: PAPI_RES_STL       = -2147483591
      INTEGER, PARAMETER :: PAPI_FP_STAL       = -2147483590
      INTEGER, PARAMETER :: PAPI_TOT_CYC       = -2147483589
      INTEGER, PARAMETER :: PAPI_LST_INS       = -2147483588
      INTEGER, PARAMETER :: PAPI_SYC_INS       = -2147483587
      INTEGER, PARAMETER :: PAPI_L1_DCH        = -2147483586
      INTEGER, PARAMETER :: PAPI_L2_DCH        = -2147483585
      INTEGER, PARAMETER :: PAPI_L1_DCA        = -2147483584
      INTEGER, PARAMETER :: PAPI_L2_DCA        = -2147483583
      INTEGER, PARAMETER :: PAPI_L3_DCA        = -2147483582
      INTEGER, PARAMETER :: PAPI_L1_DCR        = -2147483581
      INTEGER, PARAMETER :: PAPI_L2_DCR        = -2147483580
      INTEGER, PARAMETER :: PAPI_L3_DCR        = -2147483579
      INTEGER, PARAMETER :: PAPI_L1_DCW        = -2147483578
      INTEGER, PARAMETER :: PAPI_L2_DCW        = -2147483577
      INTEGER, PARAMETER :: PAPI_L3_DCW        = -2147483576
      INTEGER, PARAMETER :: PAPI_L1_ICH        = -2147483575
      INTEGER, PARAMETER :: PAPI_L2_ICH        = -2147483574
      INTEGER, PARAMETER :: PAPI_L3_ICH        = -2147483573
      INTEGER, PARAMETER :: PAPI_L1_ICA        = -2147483572
      INTEGER, PARAMETER :: PAPI_L2_ICA        = -2147483571
      INTEGER, PARAMETER :: PAPI_L3_ICA        = -2147483570
      INTEGER, PARAMETER :: PAPI_L1_ICR        = -2147483569
      INTEGER, PARAMETER :: PAPI_L2_ICR        = -2147483568
      INTEGER, PARAMETER :: PAPI_L3_ICR        = -2147483567
      INTEGER, PARAMETER :: PAPI_L1_ICW        = -2147483566
      INTEGER, PARAMETER :: PAPI_L2_ICW        = -2147483565
      INTEGER, PARAMETER :: PAPI_L3_ICW        = -2147483564
      INTEGER, PARAMETER :: PAPI_L1_TCH        = -2147483563
      INTEGER, PARAMETER :: PAPI_L2_TCH        = -2147483562
      INTEGER, PARAMETER :: PAPI_L3_TCH        = -2147483561
      INTEGER, PARAMETER :: PAPI_L1_TCA        = -2147483560
      INTEGER, PARAMETER :: PAPI_L2_TCA        = -2147483559
      INTEGER, PARAMETER :: PAPI_L3_TCA        = -2147483558
      INTEGER, PARAMETER :: PAPI_L1_TCR        = -2147483557
      INTEGER, PARAMETER :: PAPI_L2_TCR        = -2147483556
      INTEGER, PARAMETER :: PAPI_L3_TCR        = -2147483555
      INTEGER, PARAMETER :: PAPI_L1_TCW        = -2147483554
      INTEGER, PARAMETER :: PAPI_L2_TCW        = -2147483553
      INTEGER, PARAMETER :: PAPI_L3_TCW        = -2147483552
      INTEGER, PARAMETER :: PAPI_FML_INS       = -2147483551
      INTEGER, PARAMETER :: PAPI_FAD_INS       = -2147483550
      INTEGER, PARAMETER :: PAPI_FDV_INS       = -2147483549
      INTEGER, PARAMETER :: PAPI_FSQ_INS       = -2147483548
      INTEGER, PARAMETER :: PAPI_FNV_INS       = -2147483547
      INTEGER, PARAMETER :: PAPI_FP_OPS        = -2147483546
      INTEGER, PARAMETER :: PAPI_SP_OPS        = -2147483545
      INTEGER, PARAMETER :: PAPI_DP_OPS        = -2147483544
      INTEGER, PARAMETER :: PAPI_VEC_SP        = -2147483543
      INTEGER, PARAMETER :: PAPI_VEC_DP        = -2147483542

  integer, parameter :: event_number = 5
  integer, dimension(event_number) :: events
  data events /&
               PAPI_L3_TCM,&
!               PAPI_L1_DCM,&
!               PAPI_L2_TCM,&
!               PAPI_L3_DCM,&
               PAPI_L2_DCM,&
!               PAPI_L2_DCM,&
               PAPI_L3_TCA,&
               PAPI_TOT_CYC,&
               PAPI_TOT_INS&
                /
  integer(kind=8), dimension(event_number) :: counters
  character(*), parameter, dimension(event_number) :: event_name = (/&
                'PAPI_L3_TCM ',&
!                'PAPI_L1_DCM ',&
!                'PAPI_L2_TCM ',&
!                'PAPI_L3_TCM ',&
                'PAPI_L2_DCM ',&
!                'PAPI_L2_DCM ',&
                'PAPI_L3_TCA ',&
                'PAPI_TOT_CYC',&
                'PAPI_TOT_INS'&
                 /)
 
!!!  !Use arguments
!!!  call getarg(1,chain)
!!!  read(unit=chain,fmt=*) n1
!!!  call getarg(2,chain)
!!!  read(unit=chain,fmt=*) ndat

  read(unit=1,fmt=*,iostat=ierror) ndim,n1s,n1e,ndats,ndate,ntimes
  if (ierror /= 0) then
     write(*,*) "In a file 'fort.1', put a line with:"
     write(*,*) "ndim n1s n1e ndats ndate ntimes"
     write(*,*) "where:"
     write(*,*) "- ndim (1 or 3) is the dimension of the real space"
     write(*,*) "       1 do convolution from n1s to n1e (ntimes * (ndate-ndats+1))"
     write(*,*) "       3 do convolution n1=(from n1s to n1e), n2=ndats, n3=ndate"
     write(*,*) "- ntimes is the number of convolutions"
     stop
  end if

  hx=0.1e0_gp
  hy=0.1e0_gp
  hz=0.1e0_gp

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

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ekin=0.0_wp
  
  !call set_gpu_double() !after this call, all memory operations are in double precision, call set_gpu_simple() in order to have simple memory operations
!  call init_thread_engine();

  !one dimensional case
  if (ndim == 1) then
     do ndat=ndats,ndate
        do n1=n1s,n1e
           n2 = n1-2
           n3 = n1+2
           !set of one-dimensional convolutions
           !allocate arrays
           psi_in = f_malloc((/ n1, ndat, 1 /),id='psi_in')
           psi_out = f_malloc((/ ndat, n1, 1 /),id='psi_out')
           psi_cuda = f_malloc((/ ndat, n1, 1 /),id='psi_cuda')
           psi_3d_in = f_malloc((/ n1, n2, n3, ntimes /),id='psi_3d_in')
           psi_3d_out = f_malloc((/ n1, n2, n3, ntimes /),id='psi_3d_out')
           psi_3d_cuda = f_malloc((/ n1, n2, n3, ntimes /),id='psi_3d_cuda')
           psi_3d_tmp = f_malloc((/ n1, n2, n3, ntimes /),id='psi_3d_tmp')


           !initialise array
           sigma2=0.25d0*((n1*hx)**2)
           do i=1,ndat
              do i1=1,n1
!                 x=hx*real(i1-n1/2-1,kind=8)
                 !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
!                 r2=x**2
!                 arg=0.5d0*r2/sigma2
!                 tt=dexp(-arg)
                 call random_number(tt)
                 psi_in(i1,i,1)=tt
              end do
           end do

           do itimes=1,ntimes
!           call convrot_n_per_3d_simple(n1,n2,n3,psi_3d_in(1,1,1,itimes),psi_3d_out(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
!           psi_3d_in(:,:,:,:)=0.0e0_wp
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 call random_number(tt)
                 psi_3d_in(i1,i2,i3,itimes)=tt
               end do
             end do
           end do
           end do
            
           call start_counters(events,event_number,ierror)

           if (iproc == 0) then
             write(*,'(a,i7,i7,i7)')'CPU 3D Simple Convolutions, dimensions:',n1,n2,n3
           end if

           call nanosec(tsc0);
           do itimes=1,ntimes
              call convrot_n_per_3d_simple(n1,n2,n3,psi_3d_in(1,1,1,itimes),&
                   psi_3d_out(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
           end do
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(CPUtime,n1*n2*n3,32*3,ntimes)

           if (iproc == 0) then
             write(*,'(a,i7,i7,i7)')'CPU 3D Simple Transposed Convolutions, dimensions:',n1,n2,n3
           end if

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           do itimes=1,ntimes
              call convrot_n_per_3d_simple_transpose(n1,n2,n3,psi_3d_in(1,1,1,itimes),&
                   psi_3d_cuda(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           call read_print_stop_counters(counters,event_name,event_number,ierror)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)


           if (iproc == 0) then
             write(*,'(a,i7,i7,i7)')'CPU 3D Convolutions, dimensions:',n1,n2,n3
           end if

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           do itimes=1,ntimes
              call convrot_n_per_3d(n1,n2,n3,psi_3d_in(1,1,1,itimes),psi_3d_cuda(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           call read_print_stop_counters(counters,event_name,event_number,ierror)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)

           if (iproc == 0) then
             write(*,'(a,i7,i7,i7)')'CPU 3D Transposed Convolutions, dimensions:',n1,n2,n3
           end if

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           do itimes=1,ntimes
              call convrot_n_per_3d_transpose(n1,n2,n3,psi_3d_in(1,1,1,itimes),&
                   psi_3d_cuda(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
           end do
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)

           if (iproc == 0) then
             write(*,'(a,i7,i7,i7)')'CPU 3D Transposed sse Convolutions, dimensions:',n1,n2,n3
           end if

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           do itimes=1,ntimes
              call convrot_n_per_3d_sse(n1,n2,n3,psi_3d_in(1,1,1,itimes),psi_3d_cuda(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           call read_print_stop_counters(counters,event_name,event_number,ierror)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)

           if (iproc == 0) then
             write(*,'(a,i7,i7,i7)')'CPU 3D Transposed sse Convolutions T, dimensions:',n1,n2,n3
           end if

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           do itimes=1,ntimes
              call convrot_t_per_3d_sse(n1,n2,n3,psi_3d_in(1,1,1,itimes),psi_3d_cuda(1,1,1,itimes),psi_3d_tmp(1,1,1,itimes))
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           call read_print_stop_counters(counters,event_name,event_number,ierror)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)



 
           if (iproc == 0) then
             write(*,'(a,i7,i7)')'CPU Convolutions, dimensions:',n1,ndat
           end if

           !take timings
           !call system_clock(it0,count_rate,count_max)
           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_n_per(n1-1,ndat,psi_in,psi_out)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(CPUtime,n1*ndat,32,ntimes)

           !the input and output arrays must be reverted in this implementation
           !take timings

           if (iproc == 0) then
             write(*,'(a,i7,i7)')'CPU naive Convolutions, dimensions:',n1,ndat
           end if

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_naive(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
          
           if (iproc == 0) then
             write(*,'(a,i7,i7)')'CPU sse Convolutions, dimensions:',n1,ndat
           end if

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_sse(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
          
           if (iproc == 0) then
             write(*,'(a,i7,i7)')'CPU Convolutions T, dimensions:',n1,ndat
           end if

           !take timings
           !call system_clock(it0,count_rate,count_max)
           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_t_per(n1-1,ndat,psi_in,psi_out)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(CPUtime,n1*ndat,32,ntimes)

           !take timings

           if (iproc == 0) then
             write(*,'(a,i7,i7)')'CPU naive Convolutions T, dimensions:',n1,ndat
           end if

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_t_naive(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)


           if (iproc == 0) then
             write(*,'(a,i7,i7)')'CPU sse Convolutions T, dimensions:',n1,ndat
           end if

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_t_sse(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)

 
        end do
     end do
    
  else 
     print *,'wrong ndim',ndim
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call MPI_FINALIZE(ierr)

contains

  subroutine start_counters(events,event_number,ierror)
    implicit none
    integer, intent(in) :: event_number
    integer, dimension(event_number), intent(in) :: events
    integer, intent(inout) :: ierror

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call PAPIF_start_counters(events,event_number,ierror)

  end subroutine start_counters

  subroutine read_print_stop_counters(counters,event_name,event_number,ierror)
    implicit none
    integer, intent(in) :: event_number
    integer, intent(inout) :: ierror
    integer(kind=8), dimension(event_number), intent(inout) :: counters
    character(*), dimension(event_number), intent(in) :: event_name
    integer :: i

    call PAPIF_stop_counters(counters,event_number,ierror)
    if (iproc == 0) then
      do i=1,event_number
        print *,event_name(i),':',counters(i)
      enddo
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine read_print_stop_counters

  subroutine print_time(time,nbelem,nop,ntimes)
    implicit none
    real(gp),intent(in)::time
    integer,intent(in)::nbelem,nop,ntimes
  
    if (iproc == 0) then
      write(*,'(a,f9.4,1pe12.5)')'Finished. Time(ms), GFlops',&
        time*1.d3/real(ntimes,kind=8),&
        real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9)
    end if
  END SUBROUTINE print_time

  subroutine compare_time(REFtime,TESTtime,nbelem,nop,ntimes,maxdiff,threshold)
    implicit none
    real(gp),intent(in)::REFtime,TESTtime,maxdiff,threshold
    integer,intent(in)::nbelem,nop,ntimes

    if (iproc == 0) then
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

    if (iproc == 0) then
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
    end if
  END SUBROUTINE compare_3D_results

  subroutine compare_2D_results(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(dim1,dim2), psi(dim1,dim2)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2

    if (iproc == 0) then
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
    end if
  END SUBROUTINE compare_2D_results

  subroutine compare_2D_results_t(dim1, dim2, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1, dim2
    real(gp),intent(in):: psi_ref(dim1,dim2), psi(dim2,dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1,i2

    if (iproc == 0) then
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
    end if
  END SUBROUTINE compare_2D_results_t

  subroutine compare_1D_results(dim1, psi_ref, psi, maxdiff, printdiff)
    implicit none
    integer,intent(in):: dim1
    real(gp),intent(in):: psi_ref(dim1), psi(dim1)
    real(gp),intent(out):: maxdiff
    real(gp),intent(in):: printdiff
    real(gp)::comp
    integer::i1

    if (iproc == 0) then
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
    end if
  END SUBROUTINE compare_1D_results


 
subroutine convrot_n_per_simple(n1,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:(n1-1),ndat), intent(in) :: x
  real(wp), dimension(ndat,0:(n1-1)), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l
  real(wp) :: tt

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /
  do j=1,ndat
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)  
           tt=tt+x(  k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE convrot_n_per_simple

subroutine convrot_n_per_3d_simple_transpose(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(n1*n2*n3), intent(in) :: x
  real(wp), dimension(n1*n2*n3), intent(out) :: y,tmp

  call convrot_n_per_simple(n1,n2*n3,x,y)
  call convrot_n_per_simple(n2,n1*n3,y,tmp)
  call convrot_n_per_simple(n3,n1*n2,tmp,y)
  
END SUBROUTINE convrot_n_per_3d_simple_transpose

subroutine convrot_n_per_3d_transpose(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(n1*n2*n3), intent(in) :: x
  real(wp), dimension(n1*n2*n3), intent(out) :: y,tmp

  call convrot_n_per(n1-1,n2*n3,x,y)
  call convrot_n_per(n2-1,n1*n3,y,tmp)
  call convrot_n_per(n3-1,n1*n2,tmp,y)
  
END SUBROUTINE convrot_n_per_3d_transpose

subroutine convrot_n_per_3d_simple(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(in) :: x
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(out) :: y,tmp
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l,m
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  do j=0,n3-1
    do m=0,n2-1
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)   
           tt=tt+x(  k,m,j)*fil(l)
        enddo
        y(i,m,j)=tt
     enddo
    enddo
  enddo

  do j=0,n3-1
    do i=0,n1-1
     do m=0,n2-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(m+l,n2)   
           tt=tt+y(  i,k,j)*fil(l)
        enddo
        tmp(i,m,j)=tt
     enddo
    enddo
  enddo
  do m=0,n2-1
    do i=0,n1-1
     do j=0,n3-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(j+l,n3)   
           tt=tt+tmp(  i,m,k)*fil(l)
        enddo
        y(i,m,j)=tt
     enddo
    enddo
  enddo

END SUBROUTINE convrot_n_per_3d_simple

subroutine convrot_n_per_3d_sse(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(n1*n2*n3), intent(in) :: x
  real(wp), dimension(n1*n2*n3), intent(out) :: y,tmp

  call magicfilter1d_sse(n1,n2*n3,x,y)
  call magicfilter1d_sse(n2,n1*n3,y,tmp)
  call magicfilter1d_sse(n3,n1*n2,tmp,y)

END SUBROUTINE convrot_n_per_3d_sse

subroutine convrot_t_per_3d_sse(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(n1*n2*n3), intent(in) :: x
  real(wp), dimension(n1*n2*n3), intent(out) :: y,tmp

  call magicfilter1d_t_sse(n1,n2*n3,x,y)
  call magicfilter1d_t_sse(n2,n1*n3,y,tmp)
  call magicfilter1d_t_sse(n3,n1*n2,tmp,y)

end subroutine convrot_t_per_3d_sse


subroutine convrot_n_per_3d(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(in) :: x
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(out) :: y,tmp
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l,m
  real(wp) :: fill,tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  integer mod_arr1(lowfil:n1-1+lupfil)
  integer mod_arr2(lowfil:n2-1+lupfil)
  integer mod_arr3(lowfil:n3-1+lupfil)
  call fill_mod_arr(mod_arr1,lowfil,n1-1+lupfil,n1)
  call fill_mod_arr(mod_arr2,lowfil,n2-1+lupfil,n2)
  call fill_mod_arr(mod_arr3,lowfil,n3-1+lupfil,n3)

  do j=0,n3-1
    do m=0,n2/8-1
     do i=0,n1-1
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        do l=lowfil,lupfil
           k=mod_arr1(i+l)
           fill=fil(l)
           tt1=tt1+x(  k,m*8+0,j)*fill
           tt2=tt2+x(  k,m*8+1,j)*fill
           tt3=tt3+x(  k,m*8+2,j)*fill
           tt4=tt4+x(  k,m*8+3,j)*fill
           tt5=tt5+x(  k,m*8+4,j)*fill
           tt6=tt6+x(  k,m*8+5,j)*fill
           tt7=tt7+x(  k,m*8+6,j)*fill
           tt8=tt8+x(  k,m*8+7,j)*fill
        enddo
        y(i,m*8+0,j)=tt1
        y(i,m*8+1,j)=tt2
        y(i,m*8+2,j)=tt3
        y(i,m*8+3,j)=tt4
        y(i,m*8+4,j)=tt5
        y(i,m*8+5,j)=tt6
        y(i,m*8+6,j)=tt7
        y(i,m*8+7,j)=tt8
     enddo
    enddo
    do m=(n2/8)*8, n2-1
     do i=0,n1-1
       tt=0.e0_wp
       do l=lowfil,lupfil
           k=mod_arr1(i+l)
           fill=fil(l)
           tt=tt+x(  k,m,j)*fill
       enddo
       y(i,m,j)=tt
     enddo
    enddo
  enddo
  do j=0,n3-1
    do i=0,n1/8-1
     do m=0,n2-1
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        do l=lowfil,lupfil
           k=mod_arr2(m+l)
           fill=fil(l)
           tt1=tt1+y(  i*8+0,k,j)*fill
           tt2=tt2+y(  i*8+1,k,j)*fill
           tt3=tt3+y(  i*8+2,k,j)*fill
           tt4=tt4+y(  i*8+3,k,j)*fill
           tt5=tt5+y(  i*8+4,k,j)*fill
           tt6=tt6+y(  i*8+5,k,j)*fill
           tt7=tt7+y(  i*8+6,k,j)*fill
           tt8=tt8+y(  i*8+7,k,j)*fill
        enddo
        tmp(i*8+0,m,j)=tt1
        tmp(i*8+1,m,j)=tt2
        tmp(i*8+2,m,j)=tt3
        tmp(i*8+3,m,j)=tt4
        tmp(i*8+4,m,j)=tt5
        tmp(i*8+5,m,j)=tt6
        tmp(i*8+6,m,j)=tt7
        tmp(i*8+7,m,j)=tt8
     enddo
    enddo
    do i=(n1/8)*8,n1-1
     do m=0, n2-1
       tt=0.e0_wp
       do l=lowfil,lupfil
           k=mod_arr2(m+l)
           fill=fil(l)
           tt=tt+y(  i,k,j)*fill
       enddo
       tmp(i,m,j)=tt
     enddo
    enddo

  enddo


  do m=0,n2-1
    do i=0,n1/8-1
     do j=0,n3-1
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        do l=lowfil,lupfil
           k=mod_arr3(j+l)
           fill=fil(l)
           tt1=tt1+tmp(  i*8+0,m,k)*fill
           tt2=tt2+tmp(  i*8+1,m,k)*fill
           tt3=tt3+tmp(  i*8+2,m,k)*fill
           tt4=tt4+tmp(  i*8+3,m,k)*fill
           tt5=tt5+tmp(  i*8+4,m,k)*fill
           tt6=tt6+tmp(  i*8+5,m,k)*fill
           tt7=tt7+tmp(  i*8+6,m,k)*fill
           tt8=tt8+tmp(  i*8+7,m,k)*fill
        enddo
        y(i*8+0,m,j)=tt1
        y(i*8+1,m,j)=tt2
        y(i*8+2,m,j)=tt3
        y(i*8+3,m,j)=tt4
        y(i*8+4,m,j)=tt5
        y(i*8+5,m,j)=tt6
        y(i*8+6,m,j)=tt7
        y(i*8+7,m,j)=tt8
     enddo
    enddo
    do i=(n1/8)*8,n1-1
     do j=0, n3-1
       tt=0.e0_wp
       do l=lowfil,lupfil
           k=mod_arr3(j+l)
           fill=fil(l)
           tt=tt+tmp(  i,m,k)*fill
       enddo
       y(i,m,j)=tt
     enddo
    enddo


  enddo

END SUBROUTINE convrot_n_per_3d


end program conv_check


