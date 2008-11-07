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
  integer :: count_rate,count_max
  integer :: n1s,n1e,ndats,ndate
  real(wp) :: tt
  real(gp) :: v,p,CPUtime,GPUtime,comp
  real(kind=4), dimension(:,:,:), allocatable :: psi_cuda,v_cuda !temporary in view of wp 
  real(kind=4) :: t0,t1,epotGPU
  real(kind=8) :: psi_GPU,v_GPU,work_GPU !pointer to the GPU  memory addresses (with norb=1)
  integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
  integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
  real(kind=4) filCUDA1(lowfil1:lupfil1) !array of filters to be passed to CUDA interface
  data filCUDA1 / &
       8.4334247333529341094733325815816e-7_4,&
       -0.1290557201342060969516786758559028e-4_4,&
       0.8762984476210559564689161894116397e-4_4,&
       -0.30158038132690463167163703826169879e-3_4,&
       0.174723713672993903449447812749852942e-2_4,&
       -0.942047030201080385922711540948195075e-2_4,&
       0.2373821463724942397566389712597274535e-1_4,&
       0.612625895831207982195380597e-1_4,&
       0.9940415697834003993178616713_4,&
       -0.604895289196983516002834636e-1_4, &
       -0.2103025160930381434955489412839065067e-1_4,&
       0.1337263414854794752733423467013220997e-1_4,&
       -0.344128144493493857280881509686821861e-2_4,&
       0.49443227688689919192282259476750972e-3_4,&
       -0.5185986881173432922848639136911487e-4_4,&
       2.72734492911979659657715313017228e-6_4 /
  real(kind=4) filCUDA2(lowfil2:lupfil2) !array of filters to be passed to CUDA interface
  data filCUDA2 / &
       2.72734492911979659657715313017228e-6_4,&
       -0.5185986881173432922848639136911487e-4_4,&
       0.49443227688689919192282259476750972e-3_4,&
       -0.344128144493493857280881509686821861e-2_4,&
       0.1337263414854794752733423467013220997e-1_4,&
       -0.2103025160930381434955489412839065067e-1_4,&
       -0.604895289196983516002834636e-1_4,&
       0.9940415697834003993178616713_4,&
       0.612625895831207982195380597e-1_4,&
       0.2373821463724942397566389712597274535e-1_4,&
       -0.942047030201080385922711540948195075e-2_4,&
       0.174723713672993903449447812749852942e-2_4,&
       -0.30158038132690463167163703826169879e-3_4,&
       0.8762984476210559564689161894116397e-4_4,&
       -0.1290557201342060969516786758559028e-4_4,&
       8.4334247333529341094733325815816e-7_4 /

!!$  !Use arguments
!!$  call getarg(1,chain)
!!$  read(unit=chain,fmt=*) n1
!!$  call getarg(2,chain)
!!$  read(unit=chain,fmt=*) ndat

  read(1,*)ndim,n1s,n1e,ndats,ndate,ntimes

  hx=0.1e0_gp
  hy=0.1e0_gp
  hz=0.1e0_gp


  !one dimensional case
  if (ndim == 1) then
     do ndat=ndats,ndate
        do n1=n1s,n1e
           !values for the grid
!!$  n1=127
           n2=50
           n3=50
           !ntimes=1


           !set of one-dimensional convolutions
           !allocate arrays
           !ndat=256!(n2+1)*(n3+1)
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

           !CPUtime=real(it1-it0,kind=8)/real(count_rate*ntimes,kind=8)
           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                CPUtime*1.d3/real(ntimes,kind=8),real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

           !print *,'starting CUDA with n,ndat',n1,ndat

           allocate(psi_cuda(n1,ndat,1+ndebug),stat=i_stat)
           call memocc(i_stat,psi_cuda,'psi_cuda',subname)
           allocate(v_cuda(ndat,n1,1+ndebug),stat=i_stat)
           call memocc(i_stat,v_cuda,'v_cuda',subname)


           !print *,'one'
           !the input and output arrays must be reverted in this implementation
           do i=1,ndat
              do i1=1,n1
                 v_cuda(i,i1,1)=real(psi_in(i1,i,1),kind=4)
                 !write(16,'(2(i6),2(1pe24.17)')i,i1,v_cuda(i,i1,1),psi_in(i1,i,1)
              end do
           end do
           !print *,'two'

           call GPU_allocate(n1*ndat,psi_GPU,i_stat)
           !print *,'three'
           call GPU_allocate(n1*ndat,work_GPU,i_stat)

           !call GPU_send(n1*ndat,v_cuda,psi_GPU,i_stat)
           call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

           !print *,'starting CUDA'

           !now the CUDA part
           !take timings

           write(*,'(a,i6,i6)')'GPU Convolutions, dimensions:',n1,ndat

           !call system_clock(it0,count_rate,count_max)
           call cpu_time(t0)
           do i=1,ntimes
              !call cuda C interface
              !call m1dconv(n1-1,ndat,work_GPU,psi_GPU,filCUDA1,lowfil1,lupfil1)
              !call n1dconv(n1-1,ndat,work_GPU,psi_GPU,filCUDA1,-lowfil1,lupfil1)
              call g1dconv(n1-1,ndat,work_GPU,psi_GPU,filCUDA1,-lowfil1,lupfil1)

              !call localpotential(n1-1,0,ndat-1,work_GPU,psi_GPU,psi_GPU,epotGPU)
           end do
           call cpu_time(t1)
           !call system_clock(it1,count_rate,count_max)

           !GPUtime=real(it1-it0,kind=8)/real(count_rate*ntimes,kind=8)
           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

           write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
                GPUtime*1.d3/real(ntimes,kind=8),real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
           !print *,'ending CUDA'

           call GPU_receive(n1*ndat,psi_cuda,psi_GPU,i_stat)

           call GPU_deallocate(psi_GPU,i_stat)
           call GPU_deallocate(work_GPU,i_stat)



           !check the differences between the results
           maxdiff=0.d0
           i1_max=1
           i_max=1
           do i=1,ndat
              do i1=1,n1
                 write(17,'(2(i6),2(1pe24.17))')i,i1,v_cuda(i,i1,1),psi_cuda(i1,i,1)
                 !write(17,'(2(i6),2(1pe24.17))')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
                 !comp=abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8))
                 comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    i1_max=i1
                    i_max=i
                 end if
              end do
           end do

           !!print *,'timings,difference',CPUtime,GPUtime,maxdiff

           !print *,'i1,maxdiff',i_max,i1_max,v_cuda(i_max,i1_max,1),psi_cuda(i1_max,i_max,1)
           !!print *,'i1,maxdiff',i_max,i1_max,psi_out(i_max,i1_max,1),psi_cuda(i1_max,i_max,1)

           if (maxdiff <= 3.d-7) then
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)
           else
              write(*,'(a,i6,i6,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
                   'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
                   n1,ndat,CPUtime/GPUtime,maxdiff,&
                   CPUtime*1.d3/real(ntimes,kind=8),real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9),&
                   GPUtime*1.d3/real(ntimes,kind=8),real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9),&
                   '<<<< WARNING' 
           end if

           i_all=-product(shape(psi_in))
           deallocate(psi_in,stat=i_stat)
           call memocc(i_stat,i_all,'psi_in',subname)
           i_all=-product(shape(psi_out))
           deallocate(psi_out,stat=i_stat)
           call memocc(i_stat,i_all,'psi_out',subname)

           i_all=-product(shape(psi_cuda))
           deallocate(psi_cuda,stat=i_stat)
           call memocc(i_stat,i_all,'psi_cuda',subname)
           i_all=-product(shape(v_cuda))
           deallocate(v_cuda,stat=i_stat)
           call memocc(i_stat,i_all,'v_cuda',subname)
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

     CPUtime=real(t1-t0,kind=8)/real(ntimes,kind=8)


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

     GPUtime=real(t1-t0,kind=8)/real(ntimes,kind=8)

     write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
          GPUtime*1.d3/real(ntimes,kind=8),&
          real(n1*n2*n3*ntimes,kind=8)*192.d0/(GPUtime*1.d9)

     !deallocate GPU memory
     call GPU_deallocate(psi_GPU,i_stat)
     call GPU_deallocate(v_GPU,i_stat)
     call GPU_deallocate(work_GPU,i_stat)

     !call CUDA_DEALLOCATE_MEM(1,psi_GPU,v_GPU,work_GPU)

     print *,'epot=',epotGPU
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
