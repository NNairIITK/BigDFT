!> @file
!!  Test of convolutions
!! @author
!!    Copyright (C) 2008-2014 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Program test for the convolution in GPU with CUDA
program conv_check_cuda
  use module_base
  implicit none
  integer  :: n1 !n(c) n2,n3
  real(gp) :: hx,r2,sigma2,x,y,z,maxdiff,arg !n(c) epot,hy,hz
  real(wp), dimension(:,:,:), allocatable :: psi_in,psi_out !n(c) pot,psir
  !local variables
  character(len=*), parameter :: subname='conv_check'
  !n(c) character(len=50) :: chain
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat,it0,it1,ndim,itimes !n(c) i1_max,i_max
  integer :: l,ierror,i1s,i1e !n(c) count_rate,count_max
  integer :: n1s,n1e,ndats,ndate,nvctr_cf,nseg,iseg
  real(wp) :: tt,scale
  real(gp) :: CPUtime,GPUtime,comp,ekin !n(c) v,p
  !n(c) real(gp), dimension(3) :: hgridh
  integer, dimension(:), allocatable :: keyv,modarr
  integer, dimension(:,:), allocatable :: keyg
  real(kind=8), dimension(:), allocatable :: psi !temporary in view of wp 
  real(kind=8), dimension(:,:,:), allocatable :: psi_cuda,v_cuda !temporary in view of wp 
  real(kind=8) :: ekinGPUd
  !n(c) real(kind=4) :: t0,t1,epotGPU,ekinGPU
  real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU,keys_GPU !pointer to the GPU  memory addresses (with norb=1)
  !n(c) integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
  !n(c) integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
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
  
  !call set_gpu_double() !after this call, all memory operations are in double precision, call set_gpu_simple() in order to have simple memory operations


  !one dimensional case
  if (ndim == 1) then
     do ndat=ndats,ndate
        do n1=n1s,n1e
           !set of one-dimensional convolutions
           !allocate arrays
           psi_in = f_malloc((/ n1, ndat, 1 /),id='psi_in')
           psi_out = f_malloc((/ ndat, n1, 1 /),id='psi_out')

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

           call nanosec_cuda(tsc0)

           do i=1,ntimes
              call convrot_n_per(n1-1,ndat,psi_in,psi_out)
           end do

           call nanosec_cuda(tsc1)

           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(CPUtime,n1*ndat,32,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                CPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

           psi_cuda = f_malloc((/ n1, ndat, 1 /),id='psi_cuda')
           v_cuda = f_malloc((/ ndat, n1, 1 /),id='v_cuda')

           !the input and output arrays must be reverted in this implementation
           do i=1,ndat
              do i1=1,n1
                 v_cuda(i,i1,1)=real(psi_in(i1,i,1),kind=8)
              end do
           end do

!!!           call GPU_allocate(n1*ndat,psi_GPU,i_stat)
!!!           call GPU_allocate(n1*ndat,work_GPU,i_stat)
           call sg_gpu_malloc(psi_GPU,n1*ndat,8,i_stat)
           call sg_gpu_malloc(work_GPU,n1*ndat,8,i_stat)


           call sg_send_mem_instantaneously(work_GPU,v_cuda,n1*ndat,8,i_stat)
           !call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

           !now the CUDA part
           !take timings

           write(*,'(a,i6,i6)')'GPU Convolutions, dimensions:',n1,ndat

           call nanosec_cuda(tsc0)

           do i=1,ntimes
              call magicfilter1d(n1-1,ndat,work_GPU,psi_GPU)
           end do
           call nanosec_cuda(tsc1)
           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,32,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                GPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

           call sg_recv_mem_instantaneously(psi_cuda,psi_GPU,n1*ndat,8,i_stat)

           !call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)

           call sg_gpu_free(psi_GPU,i_stat)
           call sg_gpu_free(work_GPU,i_stat)

!!!           call GPU_deallocate(psi_GPU,i_stat)
!!!           call GPU_deallocate(work_GPU,i_stat)

           !check the differences between the results
           maxdiff=0.d0
           !n(c) i1_max=1
           !n(c) i_max=1
           do i=1,ndat
              do i1=1,n1
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                 comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    !n(c) i1_max=i1
                    !n(c) i_max=i
                 end if
              end do
           end do
           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
!!$           if (maxdiff <= 3.d-7) then
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
!!$           else
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9),&
!!$                   '<<<< WARNING' 
!!$           end if

           !**************************************************kinetic term

           write(*,'(a,i6,i6)')'CPU Kinetic, dimensions:',n1,ndat

           modarr = f_malloc(lowfilK.to.n1-1+lupfilK,id='modarr')
           call fill_mod_arr(modarr,lowfilK,n1-1+lupfilK,n1)


           psi_out=0.d0
           !take timings
           call nanosec_cuda(tsc0)
           do itimes=1,ntimes
              ekin=0.0_gp
!!!              do i2=1,ndat
!!!                 do i1=1,n1
!!!                    tt=0.0_wp
!!!                    do l=lowfilK,lupfilK
!!!                       j=modulo(i1-1+l,n1)+1
!!!                       tt=tt+psi_in(j   ,i2,1)*fil(l)
!!!                    enddo
!!!                    psi_out(i2,i1,1)=tt
!!!                    ekin=ekin+psi_in(i1,i2,1)*tt
!!!                 enddo
!!!              end do
              call conv_kin_x(psi_in,psi_out,ndat,ekin)   

           end do
           call nanosec_cuda(tsc1)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call f_free(modarr)

           call print_time(CPUtime,n1*ndat,45,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                CPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)


           print *,'ekin',ekin

           call sg_gpu_malloc(psi_GPU,n1*ndat,8,i_stat)
           call sg_gpu_malloc(work_GPU,n1*ndat,8,i_stat)
           call sg_gpu_malloc(work2_GPU,n1*ndat,8,i_stat)
           call sg_gpu_malloc(v_GPU,n1*ndat,8,i_stat)


!!!           call GPU_allocate(n1*ndat,psi_GPU,i_stat)
!!!           call GPU_allocate(n1*ndat,work_GPU,i_stat)
!!!           call GPU_allocate(n1*ndat,work2_GPU,i_stat)
!!!           call GPU_allocate(n1*ndat,v_GPU,i_stat)

           call sg_send_mem_instantaneously(work_GPU,v_cuda,n1*ndat,8,i_stat)

           !now the CUDA part
           !take timings
           write(*,'(a,i6,i6)')'GPU Kinetic, dimensions:',n1,ndat

           call nanosec_cuda(tsc0)
           do i=1,ntimes
              call kinetic1d(n1-1,ndat,hx,0.d0,&
                   work_GPU,psi_GPU,work2_GPU,v_GPU,ekinGPUd)
           end do
           call nanosec_cuda(tsc1)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,n1*ndat,45,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                GPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

           !call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)
           call sg_recv_mem_instantaneously(psi_cuda,psi_GPU,n1*ndat,8,i_stat)

           call sg_gpu_free(v_GPU,i_stat)
           call sg_gpu_free(psi_GPU,i_stat)
           call sg_gpu_free(work_GPU,i_stat)
           call sg_gpu_free(work2_GPU,i_stat)

           print *,'ekinGPU',ekinGPUd

           !check the differences between the results
           maxdiff=0.d0
           !n(c) i1_max=1
           !n(c) i_max=1
           do i=1,ndat
              do i1=1,n1
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                 comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    !n(c) i1_max=i1
                    !n(c) i_max=i
                 end if
              end do
           end do
           call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-4)
!!$           if (maxdiff <= 3.d-4) then
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
!!$           else
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9),&
!!$                   '<<<< WARNING' 
!!$           end if
           

           !**************************************************wavelet transformations
           if (modulo(n1,2) == 0) then

              write(*,'(a,i6,i6)')'CPU Analysis, dimensions:',n1,ndat

              !take timings
              !call system_clock(it0,count_rate,count_max)
              call nanosec_cuda(tsc0)
              do i=1,ntimes
                 call ana_rot_per(n1/2-1,ndat,psi_in,psi_out)
              end do
              call nanosec_cuda(tsc1)

              CPUtime=real(tsc1-tsc0,kind=8)*1d-9

              call print_time(CPUtime,n1*ndat,32,ntimes)
!!$              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

!!!              call GPU_allocate(n1*ndat,psi_GPU,i_stat)
!!!              call GPU_allocate(n1*ndat,work_GPU,i_stat)
              call sg_gpu_malloc(psi_GPU,n1*ndat,8,i_stat)
              call sg_gpu_malloc(work_GPU,n1*ndat,8,i_stat)

              !call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)
              call sg_send_mem_instantaneously(work_GPU,v_cuda,n1*ndat,8,i_stat)

              !now the CUDA part
              !take timings

              write(*,'(a,i6,i6)')'GPU Analysis, dimensions:',n1,ndat

              call nanosec_cuda(tsc0)
              do i=1,ntimes
                 call ana1d(n1/2-1,ndat,work_GPU,psi_GPU)
              end do
              call nanosec_cuda(tsc1)
              GPUtime=real(tsc1-tsc0,kind=8)*1d-9

              call print_time(GPUtime,n1*ndat,32,ntimes)
!!$              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

              call sg_recv_mem_instantaneously(psi_cuda,psi_GPU,n1*ndat,8,i_stat)

              call sg_gpu_free(psi_GPU,i_stat)
              call sg_gpu_free(work_GPU,i_stat)

              !check the differences between the results
              maxdiff=0.d0
              !n(c) i1_max=1
              !n(c) i_max=1
              do i=1,ndat
                 do i1=1,n1
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                    comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                    !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                    if (comp > maxdiff) then
                       maxdiff=comp
                       !n(c) i1_max=i1
                       !n(c) i_max=i
                    end if
                 end do
              end do
              call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
!!$              if (maxdiff <= 3.d-7) then
!!$                 write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
!!$                      'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                      n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                      CPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                      GPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
!!$              else
!!$                 write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
!!$                      'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                      n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                      CPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                      GPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9),&
!!$                      '<<<< WARNING' 
!!$              end if

              write(*,'(a,i6,i6)')'CPU Synthesis, dimensions:',n1,ndat

              !take timings
              !call system_clock(it0,count_rate,count_max)
              call nanosec_cuda(tsc0)
              do i=1,ntimes
                 call syn_rot_per(n1/2-1,ndat,psi_in,psi_out)
              end do
              call nanosec_cuda(tsc1)

              CPUtime=real(tsc1-tsc0,kind=8)*1d-9

              call print_time(CPUtime,n1*ndat,32,ntimes)
!!$              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

!!!              do i=1,ndat
!!!                 do i1=1,n1
!!!                    v_cuda(i,i1,1)=real(i1,kind=8)+1.d-4*real(i,kind=8)
!!!                 end do
!!!              end do

              call sg_gpu_malloc(psi_GPU,n1*ndat,8,i_stat)
              call sg_gpu_malloc(work_GPU,n1*ndat,8,i_stat)

              call sg_send_mem_instantaneously(work_GPU,v_cuda,n1*ndat,8,i_stat)

              !now the CUDA part
              !take timings

              write(*,'(a,i6,i6)')'GPU Synthesis, dimensions:',n1,ndat

              call nanosec_cuda(tsc0)
              do i=1,ntimes
                 call syn1d(n1/2-1,ndat,work_GPU,psi_GPU)
              end do
              call nanosec_cuda(tsc1)
              GPUtime=real(tsc1-tsc0,kind=8)*1d-9

              call print_time(GPUtime,n1*ndat,32,ntimes)
!!$              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

              call sg_recv_mem_instantaneously(psi_cuda,psi_GPU,n1*ndat,8,i_stat)

              call sg_gpu_free(psi_GPU,i_stat)
              call sg_gpu_free(work_GPU,i_stat)

              !check the differences between the results
              maxdiff=0.d0
              !n(c) i1_max=1
              !n(c) i_max=1
              do i=1,ndat
                 do i1=1,n1
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                    comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                    !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                    if (comp > maxdiff) then
                       maxdiff=comp
                       !n(c) i1_max=i1
                       !n(c) i_max=i
                    end if
                 end do
              end do
              call compare_time(CPUtime,GPUtime,n1*ndat,32,ntimes,maxdiff,3.d-7)
!!$              if (maxdiff <= 3.d-7) then
!!$                 write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
!!$                      'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                      n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                      CPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                      GPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
!!$              else
!!$                 write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
!!$                      'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                      n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                      CPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
!!$                      GPUtime*1.d3/real(ntimes,kind=8),&
!!$                      real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9),&
!!$                      '<<<< WARNING' 
!!$              end if
           end if

           call f_free(psi_out)
           call f_free(psi_in)
           call f_free(psi_cuda)
           call f_free(v_cuda)

           !**************************************************compression-decompression
           !create keys arrays
           !cubic domain
           nseg=(n1+1)*(n1+1)

           print *,'nseg=',nseg

           keyg = f_malloc((/ 2, nseg /),id='keyg')
           keyv = f_malloc(nseg,id='keyv')

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
           nvctr_cf=keyv(nseg)+i1e-i1s+1

           psi = f_malloc(8*nvctr_cf,id='psi')
           !determine the values for psi function
           do i=1,nvctr_cf
              psi(i)=real(i,kind=8)
           end do
           do i=1,nvctr_cf
              do j=1,7
                 psi(nvctr_cf+7*(i-1)+j)=real(i,kind=8)+0.1d0*real(j,kind=8)
              end do
           end do

           psi_in = f_malloc((/ (2*n1+2), (2*n1+2), (2*n1+2) /),id='psi_in')
           psi_cuda = f_malloc((/ (2*n1+2), (2*n1+2), (2*n1+2) /),id='psi_cuda')

           
           write(*,'(a,3(i6))')'CPU Uncompress, dimensions:',n1,n1,n1

           !take timings
           call nanosec_cuda(tsc0)
           do i=1,ntimes
              call uncompress(n1,n1,n1,nseg,nvctr_cf,keyg,keyv,  & 
                   nseg,nvctr_cf,keyg,keyv,psi(1),psi(nvctr_cf+1),psi_in)

              !call compress(n1,n1,0,n1,0,n1,0,n1,nseg,mvctr_cf,keyg,keyv,  & 
              !     nseg,mvctr_cf,keyg,keyv,psi_in,psi(1),psi(nvctr_cf+1))
           end do
           call nanosec_cuda(tsc1)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call print_time(CPUtime,8*nvctr_cf,1,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
!!$                CPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9)


           !now the CUDA part
           call sg_gpu_malloc(psi_GPU,8*nvctr_cf,8,i_stat)
           call sg_gpu_malloc(work_GPU,(2*n1+2)*(2*n1+2)*(2*n1+2),8,i_stat)

           call sg_send_mem_instantaneously(psi_GPU,psi,8*nvctr_cf,8,i_stat)

           !assign the keys values
           call adjust_keys_for_gpu(nseg,nseg,keyv,keyg,keyv,keyg,nvctr_cf,keys_GPU)

           !now the CUDA part
           !take timings

           write(*,'(a,3(i6))')'GPU Uncompress, dimensions:',n1,n1,n1

           call nanosec_cuda(tsc0)
           do i=1,ntimes
              call uncompressgpu(n1,n1,n1,psi_GPU,work_GPU,keys_GPU)
           end do
           call nanosec_cuda(tsc1)
           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,8*nvctr_cf,1,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
!!$                GPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)

           call sg_recv_mem_instantaneously(psi_cuda,work_GPU,(2*n1+2)*(2*n1+2)*(2*n1+2),8,i_stat)

           call sg_gpu_free(psi_GPU,i_stat)
           call sg_gpu_free(work_GPU,i_stat)
           call sg_gpu_free(keys_GPU,i_stat)


           !check the differences between the results
           maxdiff=0.d0
           !n(c) i1_max=1
           !n(c) i_max=1
           do i3=1,2*n1+2
              do i2=1,2*n1+2
                 do i1=1,2*n1+2
                    !write(17,'(3(i6),2(1pe24.17))')i1,i2,i3,&
                    !     psi_in(i1,i2,i3),psi_cuda(i1,i2,i3)
                 comp=abs(psi_in(i1,i2,i3)-real(psi_cuda(i1,i2,i3),kind=8))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    !n(c) i1_max=i1
                    !n(c) i_max=i
                 end if
                 end do
              end do
           end do
           call compare_time(CPUtime,GPUtime,n1*ndat,1,ntimes,maxdiff,1.d-12)
!!$           if (maxdiff <= 1.d-12) then
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)
!!$           else
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9),&
!!$                   '<<<< WARNING' 
!!$           end if

           call f_free(psi_cuda)
           psi_cuda = f_malloc((/ 8*nvctr_cf, 1, 1 /),id='psi_cuda')


           write(*,'(a,3(i6))')'CPU Compress, dimensions:',n1,n1,n1

           !take timings
           call nanosec_cuda(tsc0)
           do i=1,ntimes
              call compress_plain(n1,n1,0,n1,0,n1,0,n1,nseg,nvctr_cf,keyg,keyv,  & 
                   nseg,nvctr_cf,keyg,keyv,psi_in,psi(1),psi(nvctr_cf+1))
           end do
           call nanosec_cuda(tsc1)

           CPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(CPUtime,8*nvctr_cf,1,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
!!$                CPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9)


           !now the CUDA part
           call sg_gpu_malloc(psi_GPU,8*nvctr_cf,8,i_stat)
           call sg_gpu_malloc(work_GPU,(2*n1+2)*(2*n1+2)*(2*n1+2),8,i_stat)

           call sg_send_mem_instantaneously(work_GPU,psi_in,(2*n1+2)*(2*n1+2)*(2*n1+2),8,i_stat)

           call adjust_keys_for_gpu(nseg,nseg,keyv,keyg,keyv,keyg,nvctr_cf,keys_GPU)

           !now the CUDA part
           !take timings

           write(*,'(a,3(i6))')'GPU Compress, dimensions:',n1,n1,n1

           call nanosec_cuda(tsc0)
           do i=1,ntimes
              call compressgpu(n1,n1,n1,work_GPU,psi_GPU,keys_GPU)
           end do
           call nanosec_cuda(tsc1)

           GPUtime=real(tsc1-tsc0,kind=8)*1d-9

           call print_time(GPUtime,8*nvctr_cf,1,ntimes)
!!$           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
!!$                GPUtime*1.d3/real(ntimes,kind=8),&
!!$                real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)

           call sg_recv_mem_instantaneously(psi_cuda,psi_GPU,8*nvctr_cf,8,i_stat)

           call sg_gpu_free(psi_GPU,i_stat)
           call sg_gpu_free(work_GPU,i_stat)
           call sg_gpu_free(keys_GPU,i_stat)


           !check the differences between the results
           maxdiff=0.d0
           !n(c) i1_max=1
           !n(c) i_max=1
           do i1=1,8*nvctr_cf
              !write(17,'(i6,2(1pe24.17))')i1,psi(i1),psi_cuda(i1,1,1)
              comp=abs(psi(i1)-real(psi_cuda(i1,1,1),kind=8))
              if (comp > maxdiff) then
                 maxdiff=comp
                 !n(c) i1_max=i1
                 !n(c) i_max=i
              end if
           end do
           call compare_time(CPUtime,GPUtime,n1*ndat,1,ntimes,maxdiff,1.d-12)
!!$           if (maxdiff <= 1.d-12) then
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)
!!$           else
!!$              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
!!$                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
!!$                   n1,ndat,CPUtime/GPUtime,maxdiff,&
!!$                   CPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
!!$                   GPUtime*1.d3/real(ntimes,kind=8),&
!!$                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9),&
!!$                   '<<<< WARNING' 
!!$           end if

           call f_free(psi)
           call f_free(psi_in)
           call f_free(psi_cuda)
           call f_free(keyg)
           call f_free(keyv)


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
    logical, parameter :: debug=.false.

    if (debug) then
       write(*,'(a,f9.4,f10.2)')'Finished. Time(ms), GFlops',&
      time*1.d3/real(ntimes,kind=8),&
      real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(time*1.d9)
    end if

  END SUBROUTINE print_time

  subroutine compare_time(REFtime,TESTtime,nbelem,nop,ntimes,maxdiff,threshold)
    implicit none
    real(gp),intent(in)::REFtime,TESTtime,maxdiff,threshold
    integer,intent(in)::nbelem,nop,ntimes

    write(*,'(1x,a)')'| CPU: ms  |  Gflops  || GPU:  ms |  GFlops  || Ratio  | No. Elements | Max. Diff. |'

    write(*,'(1x,2(2(a,f10.2),a),a,f8.3,a,i14,a,1pe12.4,a)',advance='no')&
                                !'nbelem,REF/TEST ratio,Time,Gflops: REF,TEST',&
         '|',REFtime*1.d3/real(ntimes,kind=8),'|',& ! Time CPU (ms)
         real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(REFtime*1.d9),'|',& !Gflops CPU (ms)
         '|',TESTtime*1.d3/real(ntimes,kind=8),'|',& ! Time GPU (ms)
         real(ntimes,kind=8)*real(nbelem,kind=8)*real(nop,kind=8)/(TESTtime*1.d9),'|',&!Gflops GPU (ms)
         '|',REFtime/TESTtime,'|',& ! ratio
         nbelem,'|',& !No. elements
         maxdiff,'|' ! maxdiff
    if (maxdiff <= threshold) then
       write(*,'(a)')''
    else
       write(*,'(a)')'<<<< WARNING' 
    end if
  END SUBROUTINE compare_time


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

end program conv_check_cuda
