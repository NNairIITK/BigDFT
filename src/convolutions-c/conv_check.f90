!>    Program test for the convolution in GPU
!!
!!
!! Author:
!!
!!    Luigi Genovese
!!
!!
!! Copyright:
!!
!!    Copyright (C) 2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! CREATION DATE
!!    Septembre 2008
!!
!!
!!

program conv_check
  use module_base
  implicit none
  integer  :: n1,n2,n3,n1bis,n2bis,n3bis
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
  real(kind=8), dimension(:), allocatable :: psi, psi_d !temporary in view of wp 
  real(kind=4), dimension(:), allocatable :: psi_l !temporary in view of wp 
  real(kind=8), dimension(:,:,:), allocatable :: psi_cuda,v_cuda,psi_cuda_s,v_cuda_s,psi_cuda_t,v_cuda_t,v_cuda_str !temporary in view of wp 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_3d_cuda, psi_3d_tmp 
  real(kind=8), dimension(:,:,:,:,:), allocatable :: psi_cuda_k_in,psi_cuda_k_out,psi_cuda_k_in_bis 
  real(kind=8), dimension(:,:,:,:), allocatable :: psi_cuda_k_in_a,psi_cuda_k_out_a 
  real(kind=4), dimension(:,:,:), allocatable :: psi_cuda_l,v_cuda_l !temporary in view of wp 
  real(kind=8) :: ekinGPUd
  real(kind=4) :: t0,t1,epotGPU,ekinGPU
  real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU,keys_GPU !pointer to the GPU  memory addresses (with norb=1)
  real(kind=8) :: psi_c_GPU, psi_f_GPU, keyg_GPU, keyv_GPU
  real(kind=8) :: context,queue
  integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
  integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
  integer, parameter :: lowfilK=-14,lupfilK=14 ! kinetic term
  real(kind=8), dimension(lowfilK:lupfilK) :: fil
  integer(kind=8) :: tsc0, tsc1

 
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
           allocate(psi_in(n1,ndat,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_in,'psi_in',subname)
           allocate(psi_out(ndat,n1,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_out,'psi_out',subname)
           allocate(psi_cuda(ndat,n1,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_cuda,'psi_cuda',subname)
           allocate(psi_3d_in(n1,n2,n3,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_3d_in,'psi_3d_in',subname)
           allocate(psi_3d_out(n1,n2,n3,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_3d_out,'psi_3d_out',subname)
           allocate(psi_3d_cuda(n1,n2,n3,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_3d_cuda,'psi_3d_cuda',subname)
           allocate(psi_3d_tmp(n1,n2,n3,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_3d_tmp,'psi_3d_tmp',subname)


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

           call convrot_n_per_3d_simple(n1,n2,n3,psi_3d_in,psi_3d_out,psi_3d_tmp)
!           psi_3d_in(:,:,:,:)=0.0e0_wp
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 call random_number(tt)
                 psi_3d_in(i1,i2,i3,1)=tt
               end do
             end do
           end do
            
           write(*,'(a,i7,i7,i7)')'CPU 3D Simple Convolutions, dimensions:',n1,n2,n3

           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_n_per_3d_simple(n1,n2,n3,psi_3d_in,psi_3d_out,psi_3d_tmp)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(CPUtime,n1*n2*n3,32*3,ntimes)

           write(*,'(a,i7,i7,i7)')'CPU 3D Simple Transposed Convolutions, dimensions:',n1,n2,n3

           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_n_per_3d_simple_transpose(n1,n2,n3,psi_3d_in,psi_3d_cuda,psi_3d_tmp)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)

           write(*,'(a,i7,i7,i7)')'CPU 3D Transposed Convolutions, dimensions:',n1,n2,n3

           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_n_per_3d_transpose(n1,n2,n3,psi_3d_in,psi_3d_cuda,psi_3d_tmp)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)


           write(*,'(a,i7,i7,i7)')'CPU 3D Convolutions, dimensions:',n1,n2,n3

           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_n_per_3d(n1,n2,n3,psi_3d_in,psi_3d_cuda,psi_3d_tmp)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)

           write(*,'(a,i7,i7,i7)')'CPU 3D Transposed sse Convolutions, dimensions:',n1,n2,n3

           call nanosec(tsc0);
           do i=1,ntimes
              call convrot_n_per_3d_sse(n1,n2,n3,psi_3d_in,psi_3d_cuda,psi_3d_tmp)
           end do
           call nanosec(tsc1);
           !call system_clock(it1,count_rate,count_max)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*n2*n3,32*3,ntimes)

           call compare_3D_results(n1,n2,n3,psi_3d_out, psi_3d_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*n2*n3,32*3,ntimes,maxdiff,3.d-7)



 
           write(*,'(a,i7,i7)')'CPU Convolutions, dimensions:',n1,ndat

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

           write(*,'(a,i7,i7)')'CPU naive Convolutions, dimensions:',n1,ndat

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_naive(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
          
           write(*,'(a,i7,i7)')'CPU sse Convolutions, dimensions:',n1,ndat

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_sse(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
          
           write(*,'(a,i7,i7)')'CPU Convolutions T, dimensions:',n1,ndat

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

           write(*,'(a,i7,i7)')'CPU naive Convolutions T, dimensions:',n1,ndat

           call nanosec(tsc0);
           do i=1,ntimes
              call magicfilter1d_t_naive(n1,ndat,psi_in,psi_cuda)
           end do
           call nanosec(tsc1);

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)

           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)


           write(*,'(a,i7,i7)')'CPU sse Convolutions T, dimensions:',n1,ndat

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

contains

  subroutine print_time(time,nbelem,nop,ntimes)
    implicit none
    real(gp),intent(in)::time
    integer,intent(in)::nbelem,nop,ntimes

    write(*,'(a,f9.4,1pe12.5)')'Finished. Time(ms), GFlops',&
      time*1.d3/real(ntimes,kind=8),&
      real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9)

  END SUBROUTINE print_time

  subroutine compare_time(REFtime,TESTtime,nbelem,nop,ntimes,maxdiff,threshold)
    implicit none
    real(gp),intent(in)::REFtime,TESTtime,maxdiff,threshold
    integer,intent(in)::nbelem,nop,ntimes

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

subroutine convrot_n_per_3d(n1,n2,n3,x,y,tmp)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(in) :: x
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(out) :: y,tmp
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l,m
  real(wp) :: fill,tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8
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


