!!****p* CUDA/GPUham
!! FUNCTION
!!    Test the hamiltonian operation for GPU and compare with CPU
!!
!! AUTHOR
!!    Luigi Genovese, Matthieu Ospici
!!
!! COPYRIGHT
!!    Copyright (C) 2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! CREATION DATE
!!    December 2008
!!
!! SOURCE
!!
program GPUham
  use module_base
  implicit none
  integer  :: n1,n2,n3
  real(gp) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(wp), dimension(:,:,:), allocatable :: pot,psir,psi_in,psi_out
  !local variables
  character(len=*), parameter :: subname='conv_check'
  character(len=50) :: chain
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat,i1_max,i_max,it0,it1,ndim
  integer :: count_rate,count_max,l
  integer :: n1s,n1e,ndats,ndate
  real(wp) :: tt,scale
  real(gp) :: v,p,CPUtime,GPUtime,comp,ekin,CPUGflops,GPUGflops
  real(gp), dimension(3) :: hgridh
!  real(kind=4), dimension(:,:,:), allocatable :: psi_cuda,v_cuda !temporary in view of wp 
  real(kind=4) :: t0,t1,epotGPU,ekinGPU
real(kind=8) :: epotGPU2,ekinGPU2
  !pointer to the GPU  memory addresses (with norb=1)
  real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU,work3_GPU

  read(1,*)n1,n2,n3,ntimes

  hx=0.1e0_gp
  hy=0.1e0_gp
  hz=0.1e0_gp

  scale=real(-.5_gp/hx**2,wp)


  hgridh(1)=hx
  hgridh(2)=hy
  hgridh(3)=hz


  ekin=0.0_wp
  
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

  write(*,'(a,i6,i6,i6)')'CPU local hamiltonian, dimensions:',n1,n2,n3

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

     call convolut_kinetic_per_T(n1,n2,n3,hgridh,psi_in,psi_out,ekin)

  end do

  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)
  CPUGflops=real(n1*n2*n3*ntimes,kind=8)*366.d0/(CPUtime*1.d9)

  write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
       CPUtime*1.d3/real(ntimes,kind=8),CPUGflops

  print *,'ekin,epot=',ekin,epot





  ! Initialisation of potential energy  
  epot=0.0_gp

  call set_gpu_double() !after this call, all memory operations are in double precision, call set_gpu_simple() in order to have simple memory operations

  !allocate the GPU memory
  call GPU_allocate((n1+1)*(n2+1)*(n3+1),psi_GPU,i_stat)
  call GPU_allocate((n1+1)*(n2+1)*(n3+1),v_GPU,i_stat)
  call GPU_allocate((n1+1)*(n2+1)*(n3+1),work_GPU,i_stat)
  call GPU_allocate((n1+1)*(n2+1)*(n3+1),work2_GPU,i_stat)
  call GPU_allocate((n1+1)*(n2+1)*(n3+1),work3_GPU,i_stat)


  call GPU_send((n1+1)*(n2+1)*(n3+1),psi_in,psi_GPU,i_stat)
  call GPU_send((n1+1)*(n2+1)*(n3+1),pot,v_GPU,i_stat)

  write(*,'(a,i6,i6,i6)')'GPU Hamiltonian, dimensions:',n1,n2,n3

  call cpu_time(t0)
  do i=1,ntimes


     call kinetictermd(n1,n2,n3,&
          hx,hy,hz,0.e0,&
          psi_GPU,work2_GPU,work_GPU,work3_GPU,ekinGPU2)
     

 
     call localpotentiald(n1,n2,n3,work_GPU,psi_GPU,v_GPU,epotGPU2)

  end do
  call cpu_time(t1)

  !copy vpsi on the CPU

  call GPU_receive((n1+1)*(n2+1)*(n3+1),psi_in,work_GPU,i_stat)
  call GPU_receive((n1+1)*(n2+1)*(n3+1),pot,work2_GPU,i_stat)

  GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)
  GPUGflops=real(n1*n2*n3*ntimes,kind=8)*366.d0/(GPUtime*1.d9)

  write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
       GPUtime*1.d3/real(ntimes,kind=8),GPUGflops

  !deallocate GPU memory
  call GPU_deallocate(psi_GPU,i_stat)
  call GPU_deallocate(v_GPU,i_stat)
  call GPU_deallocate(work_GPU,i_stat)
  call GPU_deallocate(work2_GPU,i_stat)
  call GPU_deallocate(work3_GPU,i_stat)

  print *,'ekin,epot=',ekinGPU2,epotGPU2


  !add Vpsi to Kpsi
  psi_in=psi_in+pot

  !check the differences between the results
  maxdiff=0.d0
  do i3=1,n3+1
     do i2=1,n2+1
        do i1=1,n1+1
           write(17,*),i1,i2,i3,psi_out(i1,i2,i3),psi_in(i1,i2,i3)
           maxdiff=max(abs(psi_out(i1,i2,i3)-psi_in(i1,i2,i3)),maxdiff)
        end do
     end do
  end do

  if (maxdiff <= 3.d-6) then
     write(*,'(a,i6,i6,i6,3x,f9.5,1pe12.5,2(0pf9.2,0pf12.4))')&
          'n1,n2,n3,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
          n1,n2,n3,CPUtime/GPUtime,maxdiff,&
          CPUtime*1.d3/real(ntimes,kind=8),CPUGflops,&
          GPUtime*1.d3/real(ntimes,kind=8),GPUGflops
  else
     write(*,'(a,i6,i6,i6,3x,f9.5,1pe12.5,2(0pf9.2,0pf12.4),a)')&
          'n,ndat,GPU/CPU ratio,Time,Gflops: CPU,GPU',&
          n1,n2,n3,CPUtime/GPUtime,maxdiff,&
          CPUtime*1.d3/real(ntimes,kind=8),CPUGflops,&
          GPUtime*1.d3/real(ntimes,kind=8),GPUGflops,&
          '<<<< WARNING' 
  end if



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

end program GPUham

!!***
