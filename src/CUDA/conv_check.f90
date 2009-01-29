!!****p* CUDA/conv_check
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
program conv_check
  use module_base
  implicit none
  integer  :: n1,n2,n3
  real(gp) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(wp), dimension(:,:,:), allocatable :: pot,psir,psi_in,psi_out
  !local variables
  character(len=*), parameter :: subname='conv_check'
  character(len=50) :: chain
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat,i1_max,i_max,it0,it1,ndim
  integer :: count_rate,count_max,l,i1s,i1e
  integer :: n1s,n1e,ndats,ndate,nvctr_cf,nseg,iseg
  real(wp) :: tt,scale
  real(gp) :: v,p,CPUtime,GPUtime,comp,ekin
  real(gp), dimension(3) :: hgridh
  integer, dimension(:), allocatable :: keyv
  integer, dimension(:,:), allocatable :: keyg
  real(kind=8), dimension(:), allocatable :: psi !temporary in view of wp 
  real(kind=8), dimension(:,:,:), allocatable :: psi_cuda,v_cuda !temporary in view of wp 
  real(kind=8) :: ekinGPUd
  real(kind=4) :: t0,t1,epotGPU,ekinGPU
  real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU,keys_GPU !pointer to the GPU  memory addresses (with norb=1)
  integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
  integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
  integer, parameter :: lowfilK=-14,lupfilK=14 ! kinetic term
  real(kind=8), dimension(lowfilK:lupfilK) :: fil

 
!!$  !Use arguments
!!$  call getarg(1,chain)
!!$  read(unit=chain,fmt=*) n1
!!$  call getarg(2,chain)
!!$  read(unit=chain,fmt=*) ndat

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

!!$  ! second derivative filters for Daubechies 16
!!$  fil(0)=    0.e-3_wp*scale
!!$  fil(1)=    1.e-3_wp*scale
!!$  fil(2)=    2.e-3_wp*scale
!!$  fil(3)=    3.e-3_wp*scale
!!$  fil(4)=    4.e-3_wp*scale
!!$  fil(5)=    5.e-3_wp*scale
!!$  fil(6)=    6.e-3_wp*scale
!!$  fil(7)=    7.e-3_wp*scale
!!$  fil(8)=    8.e-3_wp*scale
!!$  fil(9)=    9.e-3_wp*scale
!!$  fil(10)=  10.e-3_wp*scale
!!$  fil(11)=  11.e-3_wp*scale
!!$  fil(12)=  12.e-3_wp*scale
!!$  fil(13)=  13.e-3_wp*scale
!!$  fil(14)=  14.e-3_wp*scale


  do i=1,14
     fil(-i)=fil(i)
  enddo

  ekin=0.0_wp
  
  call set_gpu_double() !after this call, all memory operations are in double precision, call set_gpu_simple() in order to have simple memory operations


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
           do i=1,ntimes
              call convrot_n_per(n1-1,ndat,psi_in,psi_out)
           end do
           call cpu_time(t1)
           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                CPUtime*1.d3/real(ntimes,kind=8),&
                real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

           allocate(psi_cuda(n1,ndat,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_cuda,'psi_cuda',subname)
           allocate(v_cuda(ndat,n1,1+ndebug),stat=i_stat)
           call memocc(i_stat,v_cuda,'v_cuda',subname)

           !the input and output arrays must be reverted in this implementation
           do i=1,ndat
              do i1=1,n1
                 v_cuda(i,i1,1)=real(psi_in(i1,i,1),kind=8)
              end do
           end do

           call GPU_allocate(n1*ndat,psi_GPU,i_stat)
           call GPU_allocate(n1*ndat,work_GPU,i_stat)

           call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

           !now the CUDA part
           !take timings

           write(*,'(a,i6,i6)')'GPU Convolutions, dimensions:',n1,ndat

           call cpu_time(t0)
           do i=1,ntimes
              call magicfilter1d(n1-1,ndat,work_GPU,psi_GPU)
           end do
           call cpu_time(t1)
           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                GPUtime*1.d3/real(ntimes,kind=8),&
                real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

           call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)

           call GPU_deallocate(psi_GPU,i_stat)
           call GPU_deallocate(work_GPU,i_stat)

           !check the differences between the results
           maxdiff=0.d0
           i1_max=1
           i_max=1
           do i=1,ndat
              do i1=1,n1
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                 comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
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

           !**************************************************kinetic term

           write(*,'(a,i6,i6)')'CPU Kinetic, dimensions:',n1,ndat
           psi_out=0.d0
           !take timings
           call cpu_time(t0)
           do i=1,ntimes
              ekin=0.0_gp
              do i2=1,ndat
                 do i1=1,n1
                    tt=0.0_wp
                    do l=lowfilK,lupfilK
                       j=modulo(i1-1+l,n1)+1
                       tt=tt+psi_in(j   ,i2,1)*fil(l)
                    enddo
                    psi_out(i2,i1,1)=tt
                    ekin=ekin+psi_in(i1,i2,1)*tt
                 enddo
              end do
           end do
           call cpu_time(t1)

           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                CPUtime*1.d3/real(ntimes,kind=8),&
                real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)


           print *,'ekin',ekin

           call GPU_allocate(n1*ndat,psi_GPU,i_stat)
           call GPU_allocate(n1*ndat,work_GPU,i_stat)
           call GPU_allocate(n1*ndat,work2_GPU,i_stat)
           call GPU_allocate(n1*ndat,v_GPU,i_stat)

           call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

           !now the CUDA part
           !take timings
           write(*,'(a,i6,i6)')'GPU Kinetic, dimensions:',n1,ndat

           call cpu_time(t0)
           do i=1,ntimes
              call kinetic1d(n1-1,ndat,hx,0.d0,&
                   work_GPU,psi_GPU,work2_GPU,v_GPU,ekinGPUd)
           end do
           call cpu_time(t1)

           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                GPUtime*1.d3/real(ntimes,kind=8),&
                real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

           call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)

           call GPU_deallocate(v_GPU,i_stat)
           call GPU_deallocate(psi_GPU,i_stat)
           call GPU_deallocate(work_GPU,i_stat)
           call GPU_deallocate(work2_GPU,i_stat)

           print *,'ekinGPU',ekinGPUd

           !check the differences between the results
           maxdiff=0.d0
           i1_max=1
           i_max=1
           do i=1,ndat
              do i1=1,n1
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                 comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    i1_max=i1
                    i_max=i
                 end if
              end do
           end do

           if (maxdiff <= 3.d-4) then
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
           

           !**************************************************wavelet transformations
           if (modulo(n1,2) == 0) then

              write(*,'(a,i6,i6)')'CPU Analisys, dimensions:',n1,ndat

              !take timings
              !call system_clock(it0,count_rate,count_max)
              call cpu_time(t0)
              do i=1,ntimes
                 call ana_rot_per(n1/2-1,ndat,psi_in,psi_out)
              end do
              call cpu_time(t1)

              CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

              call GPU_allocate(n1*ndat,psi_GPU,i_stat)
              call GPU_allocate(n1*ndat,work_GPU,i_stat)

              call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

              !now the CUDA part
              !take timings

              write(*,'(a,i6,i6)')'GPU Analisys, dimensions:',n1,ndat

              call cpu_time(t0)
              do i=1,ntimes
                 call ana1d(n1/2-1,ndat,work_GPU,psi_GPU)
              end do
              call cpu_time(t1)
              GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

              call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)

              call GPU_deallocate(psi_GPU,i_stat)
              call GPU_deallocate(work_GPU,i_stat)

              !check the differences between the results
              maxdiff=0.d0
              i1_max=1
              i_max=1
              do i=1,ndat
                 do i1=1,n1
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                    comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                    !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
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

              write(*,'(a,i6,i6)')'CPU Synthesis, dimensions:',n1,ndat

              !take timings
              !call system_clock(it0,count_rate,count_max)
              call cpu_time(t0)
              do i=1,ntimes
                 call syn_rot_per(n1/2-1,ndat,psi_in,psi_out)
              end do
              call cpu_time(t1)

              CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

!!$              do i=1,ndat
!!$                 do i1=1,n1
!!$                    v_cuda(i,i1,1)=real(i1,kind=8)+1.d-4*real(i,kind=8)
!!$                 end do
!!$              end do

              call GPU_allocate(n1*ndat,psi_GPU,i_stat)
              call GPU_allocate(n1*ndat,work_GPU,i_stat)

              call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

              !now the CUDA part
              !take timings

              write(*,'(a,i6,i6)')'GPU Synthesis, dimensions:',n1,ndat

              call cpu_time(t0)
              do i=1,ntimes
                 call syn1d(n1/2-1,ndat,work_GPU,psi_GPU)
              end do
              call cpu_time(t1)
              GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

              write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

              call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)

              call GPU_deallocate(psi_GPU,i_stat)
              call GPU_deallocate(work_GPU,i_stat)

              !check the differences between the results
              maxdiff=0.d0
              i1_max=1
              i_max=1
              do i=1,ndat
                 do i1=1,n1
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                    !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                    comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                    !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
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


           !**************************************************compression-decompression
           !create keys arrays
           !cubic domain
           nseg=(n1+1)*(n1+1)

           print *,'nseg=',nseg

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
           nvctr_cf=keyv(nseg)+i1e-i1s+1

           allocate(psi(8*nvctr_cf),stat=i_stat)
           call memocc(i_stat,psi,'psi',subname)
           !determine the values for psi function
           do i=1,nvctr_cf
              psi(i)=real(i,kind=8)
           end do
           do i=1,nvctr_cf
              do j=1,7
                 psi(nvctr_cf+7*(i-1)+j)=real(i,kind=8)+0.1d0*real(j,kind=8)
              end do
           end do

           allocate(psi_in((2*n1+2),(2*n1+2),(2*n1+2)+ndebug),stat=i_stat)
           call memocc(i_stat,psi_in,'psi_in',subname)

           allocate(psi_cuda((2*n1+2),(2*n1+2),(2*n1+2)+ndebug),stat=i_stat)
           call memocc(i_stat,psi_cuda,'psi_cuda',subname)


           
           write(*,'(a,3(i6))')'CPU Uncompress, dimensions:',n1,n1,n1

           !take timings
           call cpu_time(t0)
           do i=1,ntimes
              call uncompress(n1,n1,n1,nseg,nvctr_cf,keyg,keyv,  & 
                   nseg,nvctr_cf,keyg,keyv,psi(1),psi(nvctr_cf+1),psi_in)

              !call compress(n1,n1,n1,0,n1,0,n1,0,n1,nseg,mvctr_cf,keyg,keyv,  & 
              !     nseg,mvctr_cf,keyg,keyv,psi_in,psi(1),psi(nvctr_cf+1))
           end do
           call cpu_time(t1)

           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
                CPUtime*1.d3/real(ntimes,kind=8),&
                real(8*nvctr_cf*ntimes,kind=8)*32.d0/(CPUtime*1.d9)


           !now the CUDA part
           call GPU_allocate(8*nvctr_cf,psi_GPU,i_stat)
           call GPU_allocate((2*n1+2)*(2*n1+2)*(2*n1+2),work_GPU,i_stat)

           call GPU_send(8*nvctr_cf,psi,psi_GPU,i_stat)

           !assign the keys values
           call adjust_keys_for_gpu(nseg,nseg,keyv,keyg,keyv,keyg,nvctr_cf,keys_GPU)

           !now the CUDA part
           !take timings

           write(*,'(a,3(i6))')'GPU Uncompress, dimensions:',n1,n1,n1

           call cpu_time(t0)
           do i=1,ntimes
              call uncompressgpu(n1,n1,n1,psi_GPU,work_GPU,keys_GPU)
           end do
           call cpu_time(t1)
           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
                GPUtime*1.d3/real(ntimes,kind=8),&
                real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)

           call GPU_receive((2*n1+2)*(2*n1+2)*(2*n1+2),psi_cuda,work_GPU,i_stat)

           call GPU_deallocate(psi_GPU,i_stat)
           call GPU_deallocate(work_GPU,i_stat)
           call GPU_deallocate(keys_GPU,i_stat)


           !check the differences between the results
           maxdiff=0.d0
           i1_max=1
           i_max=1
           do i3=1,2*n1+2
              do i2=1,2*n1+2
                 do i1=1,2*n1+2
                    !write(17,'(3(i6),2(1pe24.17))')i1,i2,i3,&
                    !     psi_in(i1,i2,i3),psi_cuda(i1,i2,i3)
                 comp=abs(psi_in(i1,i2,i3)-real(psi_cuda(i1,i2,i3),kind=8))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    i1_max=i1
                    i_max=i
                 end if
                 end do
              end do
           end do
           if (maxdiff <= 1.d-12) then
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)
           else
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9),&
                   '<<<< WARNING' 
           end if

           i_all=-product(shape(psi_cuda))
           deallocate(psi_cuda,stat=i_stat)
           call memocc(i_stat,i_all,'psi_cuda',subname)
           allocate(psi_cuda(8*nvctr_cf,1,1),stat=i_stat)
           call memocc(i_stat,psi_cuda,'psi_cuda',subname)


           write(*,'(a,3(i6))')'CPU Compress, dimensions:',n1,n1,n1

           !take timings
           call cpu_time(t0)
           do i=1,ntimes
              call compress(n1,n1,n1,0,n1,0,n1,0,n1,nseg,nvctr_cf,keyg,keyv,  & 
                   nseg,nvctr_cf,keyg,keyv,psi_in,psi(1),psi(nvctr_cf+1))
           end do
           call cpu_time(t1)

           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
                CPUtime*1.d3/real(ntimes,kind=8),&
                real(8*nvctr_cf*ntimes,kind=8)*32.d0/(CPUtime*1.d9)


           !now the CUDA part
           call GPU_allocate(8*nvctr_cf,psi_GPU,i_stat)
           call GPU_allocate((2*n1+2)*(2*n1+2)*(2*n1+2),work_GPU,i_stat)

           call GPU_send((2*n1+2)*(2*n1+2)*(2*n1+2),psi_in,work_GPU,i_stat)

           call adjust_keys_for_gpu(nseg,nseg,keyv,keyg,keyv,keyg,nvctr_cf,keys_GPU)

           !now the CUDA part
           !take timings

           write(*,'(a,3(i6))')'GPU Compress, dimensions:',n1,n1,n1

           call cpu_time(t0)
           do i=1,ntimes
              call compressgpu(n1,n1,n1,work_GPU,psi_GPU,keys_GPU)
           end do
           call cpu_time(t1)
           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GCopy',&
                GPUtime*1.d3/real(ntimes,kind=8),&
                real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)

           call GPU_receive(8*nvctr_cf,psi_cuda,psi_GPU,i_stat)

           call GPU_deallocate(psi_GPU,i_stat)
           call GPU_deallocate(work_GPU,i_stat)
           call GPU_deallocate(keys_GPU,i_stat)


           !check the differences between the results
           maxdiff=0.d0
           i1_max=1
           i_max=1
           do i1=1,8*nvctr_cf
              !write(17,'(i6,2(1pe24.17))')i1,psi(i1),psi_cuda(i1,1,1)
              comp=abs(psi(i1)-real(psi_cuda(i1,1,1),kind=8))
              if (comp > maxdiff) then
                 maxdiff=comp
                 i1_max=i1
                 i_max=i
              end if
           end do
           if (maxdiff <= 1.d-12) then
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9)
           else
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),&
                   real(8*nvctr_cf*ntimes,kind=8)/(GPUtime*1.d9),&
                   '<<<< WARNING' 
           end if

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



        end do
     end do
    
     !three-dimensional case
  else if (ndim == 3) then

     n1=n1s
     n2=n1e
     n3=ndats

     !-------------------------3d case---------------------------
     !allocate arrays
     allocate(psi_in((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,psi_in,'psi_in',subname)
     allocate(psi_out((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,psi_out,'psi_out',subname)
     allocate(psir((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     allocate(pot((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)

     ! Wavefunction expressed everywhere in fine scaling functions 

     !fake initialisation, random numbers
     !here the grid spacings are the small ones
     sigma2=0.25d0*(((n1+1)*hx)**2+((n2+1)*hy)**2+((n3+1)*hz)**2)
     do i3=1,n3+1
        z=hz*real(i3-n3/2-1,kind=8)
        do i2=1,n2+1
           y=hy*real(i2-n2/2-1,kind=8)
           do i1=1,n1+1
              x=hx*real(i1-n1/2-1,kind=8)
              !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
              r2=x**2+y**2+z**2
              arg=0.5d0*r2/sigma2
              tt=dexp(-arg)
              !same initialisation for psi and pot
              psi_in(i1,i2,i3)=tt
              pot(i1,i2,i3)=tt
           end do
        end do
     end do

     write(*,'(a,i6,i6,i6)')'CPU Convolutions, dimensions:',n1,n2,n3

     !take timings
     call cpu_time(t0)

     do j=1,ntimes
        !traditional convolution, CPU
        ! psir serves as a work array	   
        call convolut_magic_n_per(n1,n2,n3,psi_in,psir,psi_out) 

        ! Initialisation of potential energy  
        epot=0.0_gp
        do i3=1,n3+1
           do i2=1,n2+1
              do i1=1,n1+1
                 v=real(pot(i1,i2,i3),gp)
                 p=real(psir(i1,i2,i3),gp)
                 tt=pot(i1,i2,i3)*psir(i1,i2,i3)
                 epot=epot+p*v*p
                 psir(i1,i2,i3)=tt
              end do
           end do
        end do

        call convolut_magic_t_per_self(n1,n2,n3,psir,psi_out)

     end do

     call cpu_time(t1)

     CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

     write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
          CPUtime*1.d3/real(ntimes,kind=8),&
          real(n1*n2*n3*ntimes,kind=8)*192.d0/(CPUtime*1.d9)

     print *,'epot=',epot

     !save the values of psi_out and of epot for GPU comparison

     !create the CUDA values
     allocate(psi_cuda((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,psi_cuda,'psi_cuda',subname)
     allocate(v_cuda((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,v_cuda,'v_cuda',subname)

     psi_cuda=real(psi_in,kind=4)
     v_cuda=real(pot,kind=4)

     ! Initialisation of potential energy  
     epot=0.0_gp

     !allocate the GPU memory
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),psi_GPU,i_stat)
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),v_GPU,i_stat)
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),work_GPU,i_stat)

     call GPU_send((n1+1)*(n2+1)*(n3+1),psi_cuda,psi_GPU,i_stat)
     call GPU_send((n1+1)*(n2+1)*(n3+1),v_cuda,v_GPU,i_stat)

     !copy the data on GPU(must be separated)
     !call CUDA_ALLOC_MEM(1,2*n1+1,2*n2+1,2*n3+1,psi_cuda,v_cuda,psi_GPU,v_GPU,work_GPU)

     write(*,'(a,i6,i6,i6)')'GPU Convolutions, dimensions:',n1,n2,n3

     call cpu_time(t0)
     do i=1,ntimes

        !calculate the potential application on GPU
        !call cuda_psi_to_vpsi(1,n1,n2,n3,psi_GPU,v_GPU,work_GPU,&
        !filCUDA1,filCUDA2,lowfil1,lupfil1,lowfil2,lupfil2)

        call localpotential(n1,n2,n3,psi_GPU,work_GPU,v_GPU,epotGPU)

     end do
     call cpu_time(t1)

     !copy vpsi on the CPU
     call GPU_receive((n1+1)*(n2+1)*(n3+1),psi_cuda,psi_GPU,i_stat)
     !call GPU_receive((n1+1)*(n2+1)*(n3+1),psi_cuda,work_GPU,i_stat)
     !     call cuda_fetch_vpsi(1,2*n1+1,2*n2+1,2*n3+1,psi_GPU,psi_cuda)

     GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

     write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
          GPUtime*1.d3/real(ntimes,kind=8),&
          real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9)

     !deallocate GPU memory
     call GPU_deallocate(psi_GPU,i_stat)
     call GPU_deallocate(v_GPU,i_stat)
     call GPU_deallocate(work_GPU,i_stat)

     !call CUDA_DEALLOCATE_MEM(1,psi_GPU,v_GPU,work_GPU)

     print *,'epot=',epotGPU
     !print *,'epot=',epot
     !check the differences between the results
     maxdiff=0.d0
     do i3=1,n3+1
        do i2=1,n2+1
           do i1=1,n1+1
              !write(17,*),i1,i2,i3,psi_out(i1,i2,i3),psi_cuda(i1,i2,i3)
              maxdiff=max(abs(psi_out(i1,i2,i3)-real(psi_cuda(i1,i2,i3),kind=8)),maxdiff)
           end do
        end do
     end do

     if (maxdiff <= 3.d-6) then
        write(*,'(a,i6,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
             'n1,n2,n3,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
             n1,n2,n3,CPUtime/GPUtime,maxdiff,&
             CPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(CPUtime*1.d9),&
             GPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9)
     else
        write(*,'(a,i6,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
             'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
             n1,n2,n3,CPUtime/GPUtime,maxdiff,&
             CPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(CPUtime*1.d9),&
             GPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9),&
             '<<<< WARNING' 
     end if



     i_all=-product(shape(psi_cuda))
     deallocate(psi_cuda,stat=i_stat)
     call memocc(i_stat,i_all,'psi_cuda',subname)
     i_all=-product(shape(v_cuda))
     deallocate(v_cuda,stat=i_stat)
     call memocc(i_stat,i_all,'v_cuda',subname)
     
     !********************************************kinetic

     ! Wavefunction expressed everywhere in fine scaling functions 

     !fake initialisation, random numbers
     !here the grid spacings are the small ones
     sigma2=0.25d0*(((n1+1)*hx)**2+((n2+1)*hy)**2+((n3+1)*hz)**2)
     do i3=1,n3+1
        z=hz*real(i3-n3/2-1,kind=8)
        do i2=1,n2+1
           y=hy*real(i2-n2/2-1,kind=8)
           do i1=1,n1+1
              x=hx*real(i1-n1/2-1,kind=8)
              !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
              r2=x**2+y**2+z**2
              arg=0.5d0*r2/sigma2
              tt=dexp(-arg)
              !same initialisation for psi and pot
              psi_in(i1,i2,i3)=tt
           end do
        end do
     end do

     hgridh(1)=hx
     hgridh(2)=hy
     hgridh(3)=hz

     psi_out=0.0_wp

     write(*,'(a,i6,i6,i6)')'CPU Kinetic, dimensions:',n1,n2,n3

     !take timings
     call cpu_time(t0)

     do j=1,ntimes

        ! compute the kinetic part and add  it to psi_out
        ! the kinetic energy is calculated at the same time
        call convolut_kinetic_per_T(n1,n2,n3,hgridh,psi_in,psi_out,ekin)

     end do

     call cpu_time(t1)

     print *,'ekin=',ekin
     ekin=0.d0
     !check the differences between the results
     maxdiff=0.d0
     do i3=1,n3+1
        do i2=1,n2+1
           do i1=1,n1+1
              ekin=ekin+psi_in(i1,i2,i3)*psi_out(i1,i2,i3)
           end do
        end do
     end do


     CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

     write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
          CPUtime*1.d3/real(ntimes,kind=8),&
          real(n1*n2*n3*ntimes,kind=8)*192.d0/(CPUtime*1.d9)

     print *,'ekin=',ekin

     !save the values of psi_out and of epot for GPU comparison

     !create the CUDA values
     allocate(psi_cuda((n1+1),(n2+1),(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,psi_cuda,'psi_cuda',subname)

     psi_cuda=real(psi_in,kind=4)

     ! Initialisation of potential energy  
     epot=0.0_gp

     !allocate the GPU memory
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),psi_GPU,i_stat)
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),v_GPU,i_stat)
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),work_GPU,i_stat)
     call GPU_allocate((n1+1)*(n2+1)*(n3+1),work2_GPU,i_stat)

     call GPU_send((n1+1)*(n2+1)*(n3+1),psi_cuda,work_GPU,i_stat)

     write(*,'(a,i6,i6,i6)')'GPU Kinetic, dimensions:',n1,n2,n3

     call cpu_time(t0)
     do i=1,ntimes

        call kineticterm(n1,n2,n3,&
             real(hx,kind=4),real(hy,kind=4),real(hz,kind=4),0.e0,&
             work_GPU,psi_GPU,work2_GPU,v_GPU,ekinGPU)

     end do
     call cpu_time(t1)

     !copy vpsi on the CPU
     call GPU_receive((n1+1)*(n2+1)*(n3+1),psi_cuda,psi_GPU,i_stat)
     !call GPU_receive((n1+1)*(n2+1)*(n3+1),psi_cuda,work_GPU,i_stat)

     GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

     write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
          GPUtime*1.d3/real(ntimes,kind=8),&
          real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9)
 
     !deallocate GPU memory
     call GPU_deallocate(psi_GPU,i_stat)
     call GPU_deallocate(v_GPU,i_stat)
     call GPU_deallocate(work_GPU,i_stat)
     call GPU_deallocate(work2_GPU,i_stat)


     print *,'ekin=',ekinGPU
     !print *,'ekin=',ekin 
     ekin=0.d0
     !check the differences between the results
     maxdiff=0.d0
     do i3=1,n3+1
        do i2=1,n2+1
           do i1=1,n1+1
              !write(17,*),i1,i2,i3,psi_out(i1,i2,i3),psi_cuda(i1,i2,i3)
              maxdiff=max(abs(psi_out(i1,i2,i3)-real(psi_cuda(i1,i2,i3),kind=8)),maxdiff)
              ekin=ekin+psi_in(i1,i2,i3)*real(psi_cuda(i1,i2,i3),kind=8)
           end do
        end do
     end do

     print *,'ekin=',ekin
     if (maxdiff <= 3.d-4) then
        write(*,'(a,i6,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
             'n1,n2,n3,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
             n1,n2,n3,CPUtime/GPUtime,maxdiff,&
             CPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(CPUtime*1.d9),&
             GPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9)
     else
        write(*,'(a,i6,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
             'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
             n1,n2,n3,CPUtime/GPUtime,maxdiff,&
             CPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(CPUtime*1.d9),&
             GPUtime*1.d3/real(ntimes,kind=8),&
             real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9),&
             '<<<< WARNING' 
     end if



     i_all=-product(shape(psi_cuda))
     deallocate(psi_cuda,stat=i_stat)
     call memocc(i_stat,i_all,'psi_cuda',subname)


     i_all=-product(shape(psi_in))
     deallocate(psi_in,stat=i_stat)
     call memocc(i_stat,i_all,'psi_in',subname)
     i_all=-product(shape(psi_out))
     deallocate(psi_out,stat=i_stat)
     call memocc(i_stat,i_all,'psi_out',subname)
     i_all=-product(shape(psir))
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)
     i_all=-product(shape(pot))
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)

  else 
     print *,'wrong ndim',ndim
  end if
 
end program conv_check

!!***
