!> @file
!! Program of comparing of the application of the Hamiltonian for GPU and CPU
!! @author
!!    Luigi Genovese, Matthieu Ospici
!!    Copyright (C) 2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>    Test the hamiltonian operation for GPU and compare with CPU
program GPUham
  use module_base
  use iso_c_binding 

  implicit none
  integer  :: n1,n2,n3
  real(gp) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(wp), dimension(:,:,:), pointer :: psir,psifscf,psi_out,pot,psi_in


  !local variables
  character(len=*), parameter :: subname='conv_check'
  character(len=50) :: chain
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat,i1_max,i_max,it0,it1,ndim
  integer :: count_rate,count_max,l
  integer :: n1s,n1e,ndats,ndate
  real(wp) :: tt,scale
  real(gp) :: v,p,CPUtime,GPUtime,comp,ekin,CPUGflops,GPUGflops,epot_p,ALLOCtime, SENDtime,RECVtime
  real(gp), dimension(3) :: hgridh
!  real(kind=4), dimension(:,:,:), allocatable :: psi_cuda,v_cuda !temporary in view of wp 
  real(kind=4) :: t0,t1,epotGPU,ekinGPU
  real(kind=8) :: epotGPU2,ekinGPU2
  !pointer to the GPU  memory addresses (with norb=1)
  real(kind=8) :: psi_GPU,v_GPU,work_GPU,work2_GPU,out_GPU
  integer::   res

  type(C_PTR) :: cptr_psi_in, cptr_pot, cptr_psifscf


  read(1,*)n1,n2,n3,ntimes
  hx=0.21376778958023715_gp
  !hx=0.8e0_gp
  hy=hx
  hz=hx

  scale=real(-.5_gp/hx**2,wp)


  hgridh(1)=hx
  hgridh(2)=hy
  hgridh(3)=hz


  ekin=0.0_wp
  

  !init cublas
  call cublas_init()


  call set_gpu_double() 


  !allocate arrays
  psi_in = f_malloc_ptr((/ 2*(n1+1), 2*(n2+1), 2*(n3+1) /),id='psi_in')


 ! call cpu_pinned_allocation(2*(n1+1)*2*(n2+1)*2*(n3+1)+ndebug,cptr_psi_in,i_stat)
 ! call c_f_pointer(cptr_psi_in, psi_in, (/ 2*(n1+1) ,2*(n2+1),2*(n3+1)+ndebug /) )

  psi_out = f_malloc_ptr((/ 2*(n1+1), 2*(n2+1), 2*(n3+1) /),id='psi_out')
  psir = f_malloc_ptr((/ 2*(n1+1), 2*(n2+1), 2*(n3+1) /),id='psir')
  psifscf = f_malloc_ptr((/ 2*(n1+1), 2*(n2+1), 2*(n3+1) /),id='psifscf')

! call cpu_pinned_allocation(2*(n1+1)*2*(n2+1)*2*(n3+1)+ndebug,cptr_psifscf,i_stat)

!  call c_f_pointer(cptr_psifscf, psifscf, (/ 2*(n1+1),2*(n2+1),2*(n3+1)+ndebug /) )


  pot = f_malloc_ptr((/ 2*(n1+1), 2*(n2+1), 2*(n3+1) /),id='pot')

 ! call cpu_pinned_allocation(2*(n1+1)*2*(n2+1)*2*(n3+1)+ndebug,cptr_pot,i_stat)

 ! call c_f_pointer(cptr_pot, pot, (/ 2*(n1+1),2*(n2+1),2*(n3+1)+ndebug /) )



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
           psi_in(2*i1,2*i2,2*i3)=tt
           psi_in(2*i1-1,2*i2-1,2*i3-1)=tt
           pot(2*i1,2*i2,2*i3)=1.d0!tt
           pot(2*i1-1,2*i2-1,2*i3-1)=1.d0!tt

        end do
     end do
  end do

  pot=1.d0

  psi_in=1.d0/sqrt(real(8*(n1+1)*(n2+1)*(n3+1),wp))

  psifscf=psi_in

  write(*,'(a,i6,i6,i6)')'CPU local hamiltonian, dimensions:',2*n1+1,2*n2+1,2*n3+1

  !take timings
  call cpu_time(t0)

  do j=1,ntimes
     !traditional convolution, CPU
     
     ! calculate fine scaling functions
     call synthese_per_self(n1,n2,n3,psifscf,psi_out)

     ! psir serves as a work array	   
     call convolut_magic_n_per(2*n1+1,2*n2+1,2*n3+1,psi_out,psir,psifscf) 

     !$omp parallel default(private)&
     !$omp shared(pot,psir,n1,n2,n3,epot)

     epot_p=0._gp
     !$omp do
     do i3=1,2*n3+2
        do i2=1,2*n2+2
           do i1=1,2*n1+2
              v=real(pot(i1,i2,i3),gp)
              p=real(psir(i1,i2,i3),gp)
              tt=pot(i1,i2,i3)*psir(i1,i2,i3)
              epot_p=epot_p+p*v*p
              psir(i1,i2,i3)=tt
           end do
        end do
     end do
     !$omp end do

     !$omp critical
     epot=epot+epot_p 
     !$omp end critical

     !$omp end parallel

     call convolut_magic_t_per_self(2*n1+1,2*n2+1,2*n3+1,psir,psifscf)

     ! compute the kinetic part and add  it to psifscf
     ! the kinetic energy is calculated at the same time

     call convolut_kinetic_per_t(2*n1+1,2*n2+1,2*n3+1,hgridh,psi_out,psifscf,ekin)

     !the output is psi_out
     call analyse_per_self(n1,n2,n3,psifscf,psi_out)

  end do

  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)
  CPUGflops=8.d0*real(n1*n2*n3*ntimes,kind=8)*366.d0/(CPUtime*1.d9)

  write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
       CPUtime*1.d3/real(ntimes,kind=8),CPUGflops

  print *,'ekin,epot=',ekin,epot


  ! Initialisation of potential energy  
  epot=0.0_gp

  !after this call, all memory operations are in double precision, 
  !call set_gpu_simple() in order to have simple memory operations

  call cpu_time(t0)


  !allocate the GPU memory
  call GPU_allocate(8*(n1+1)*(n2+1)*(n3+1),psi_GPU,i_stat)
  call GPU_allocate(8*(n1+1)*(n2+1)*(n3+1),v_GPU,i_stat)
  call GPU_allocate(8*(n1+1)*(n2+1)*(n3+1),work_GPU,i_stat)
  call GPU_allocate(8*(n1+1)*(n2+1)*(n3+1),work2_GPU,i_stat)
  call GPU_allocate(8*(n1+1)*(n2+1)*(n3+1),out_GPU,i_stat)

  call cpu_time(t1)
  ALLOCtime=real(t1-t0,kind=8)

  call cpu_time(t0)

  call GPU_send(8*(n1+1)*(n2+1)*(n3+1),psi_in,psi_GPU,i_stat)
  call GPU_send(8*(n1+1)*(n2+1)*(n3+1),pot,v_GPU,i_stat)

  call cpu_time(t1)

  SENDtime=real(t1-t0,kind=8)
  write(*,'(a,i6,i6,i6)')'GPU Hamiltonian, dimensions:',2*n1+1,2*n2+1,2*n3+1

  call cpu_time(t0)
  do i=1,ntimes

     call gpuhamilt(n1,n2,n3,&
          hx,hy,hz,&
          psi_GPU,out_GPU,v_GPU,work_GPU,work2_GPU,epotGPU2,ekinGPU2)
 
  end do
  call cpu_time(t1)

 
  GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)
  GPUGflops=real(n1*n2*n3*ntimes,kind=8)*366.d0/(GPUtime*1.d9)

  psifscf=15.d0

  call cpu_time(t0)
  call GPU_receive(8*(n1+1)*(n2+1)*(n3+1),psifscf,out_GPU,i_stat)
  call cpu_time(t1)

  RECVtime=real(t1-t0,kind=8)


  write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
       GPUtime*1.d3/real(ntimes,kind=8),GPUGflops

  !deallocate GPU memory
  call GPU_deallocate(psi_GPU,i_stat)
  call GPU_deallocate(v_GPU,i_stat)
  call GPU_deallocate(work_GPU,i_stat)
  call GPU_deallocate(work2_GPU,i_stat)
  call GPU_deallocate(out_GPU,i_stat)

  print *,'ekin,epot=',ekinGPU2,epotGPU2

  print *,'ALOCtime, SENDtime, GPUtime, RECVtime : ',ALLOCtime*1d3,SENDtime*1d3,GPUtime*1d3,RECVtime*1d3
  !check the differences between the results
  maxdiff=0.d0
  do i3=1,2*(n3+1)
     do i2=1,2*(n2+1)
        do i1=1,2*(n1+1)
           write(17,*),i1,i2,i3,psi_out(i1,i2,i3),psifscf(i1,i2,i3)
           maxdiff=max(abs(psi_out(i1,i2,i3)-psifscf(i1,i2,i3)),maxdiff)
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



 ! call cpu_pinned_deallocation(cptr_psi_in,i_stat)

  call f_free_ptr(psi_in)


  call f_free_ptr(psi_out)
  call f_free_ptr(psir)
  call f_free_ptr(psifscf)
  call f_free_ptr(pot)

 ! call cpu_pinned_deallocation(cptr_pot,i_stat)



 call cublas_shutdown()

end program GPUham
