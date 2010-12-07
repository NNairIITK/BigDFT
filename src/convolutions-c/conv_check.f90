!!****p* convolutions-c/conv_check
!! FUNCTION
!!    Program test for the convolution in GPU
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2008-2010 BigDFT group 
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
program conv_check
  use module_base
  implicit none
  integer  :: n1,n2,n3,n1bis,n2bis,n3bis
  real(gp) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(wp), dimension(:,:,:), allocatable :: pot,psir,psi_in,psi_out,psi_out_s
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
  call init_thread_engine();

  !one dimensional case
  if (ndim == 1) then
     do ndat=ndats,ndate
        do n1=n1s,n1e
           !set of one-dimensional convolutions
           !allocate arrays
           allocate(psi_in(n1,ndat,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_in,'psi_in',subname)
           allocate(psi_out(ndat,n1,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_out,'psi_out',subname)
           allocate(psi_cuda(ndat,n1,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_cuda,'psi_cuda',subname)

           !initialise array
           sigma2=0.25d0*((n1*hx)**2)
           do i=1,ndat
              do i1=1,n1
                 x=hx*real(i1-n1/2-1,kind=8)
                 !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
                 r2=x**2
                 arg=0.5d0*r2/sigma2
                 tt=dexp(-arg)

                 psi_in(i1,i,1)=tt
              end do
           end do

           write(*,'(a,i6,i6)')'CPU Convolutions, dimensions:',n1,ndat

           !take timings
           !call system_clock(it0,count_rate,count_max)
           call cpu_time(t0)
           call rdtsc(tsc0)
           do i=1,ntimes
              call convrot_n_per(n1-1,ndat,psi_in,psi_out)
           end do
           call rdtsc(tsc1)
           call cpu_time(t1)
           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                CPUtime*1.d3/real(ntimes,kind=8),&
                real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)
           print *,"tsc : ",tsc1-tsc0

           !the input and output arrays must be reverted in this implementation
           !take timings

           write(*,'(a,i6,i6)')'CPU C Convolutions, dimensions:',n1,ndat

           call cpu_time(t0)
           call rdtsc(tsc0)
           do i=1,ntimes
              call magicfilter1d_d_par(n1,ndat,psi_in,psi_cuda)
           end do
           call rdtsc(tsc1)
           call cpu_time(t1)
           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           print *,"tsc : ",tsc1-tsc0
           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                GPUtime*1.d3/real(ntimes,kind=8),&
                real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

           !check the differences between the results
           maxdiff=0.d0
           i1_max=1
           i_max=1
           do i=1,ndat
              do i1=1,n1
                 comp=abs(psi_out(i,i1,1)-real(psi_cuda(i,i1,1),kind=8))
!                 write(*,*)psi_out(i,i1,1),psi_cuda(i,i1,1)
                 if (comp > maxdiff) then
                    maxdiff=comp
                    i1_max=i1
                    i_max=i
                 end if
              end do
           end do

           if (maxdiff <= 3.d-7) then
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
           else
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9),&
                   '<<<< WARNING' 
           end if

        end do
     end do
    
  else 
     print *,'wrong ndim',ndim
  end if

contains

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
  end subroutine conv_kin_x

END PROGRAM conv_check
!!***
