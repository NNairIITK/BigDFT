!!****f* CUDA/conv_check
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


!if errorFound becomes 1, error during calculation
subroutine check_init(errorFound)
 ! use module_base
  implicit none
  integer  :: n1,n2,n3
  real(kind=8) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(kind=8), dimension(:,:,:), allocatable :: pot,psir,psi_in,psi_out
  !local variables
  character(len=*), parameter :: subname='conv_check'
  character(len=50) :: chain
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat,i1_max,i_max,it0,it1,ndim,itimes
  integer :: count_rate,count_max,l,ierror,i1s,i1e
  integer :: n1s,n1e,ndats,ndate,nvctr_cf,nseg,iseg
  real(kind=8) :: tt,scale
  real(kind=8) :: v,p,CPUtime,GPUtime,comp,ekin
  real(kind=8), dimension(3) :: hgridh
  integer, dimension(:), allocatable :: keyv,modarr
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

  integer, intent(out) :: errorFound
 
!!$  !Use arguments
!!$  call getarg(1,chain)
!!$  read(unit=chain,fmt=*) n1
!!$  call getarg(2,chain)
!!$  read(unit=chain,fmt=*) ndat

 ! read(unit=1,fmt=*,iostat=ierror) ndim,n1s,n1e,ndats,ndate,ntimes


  ndim = 1
  n1s = 64
  n1e = 64
  ndats = 2048
  ndate =  2048
  ntimes = 1


  hx=0.1e0
  hy=0.1e0
  hz=0.1e0

  scale=real(-.5/hx**2,kind=8)

  ! second derivative filters for Daubechies 16
  fil(0)=   -3.5536922899131901941296809374e0*scale
  fil(1)=    2.2191465938911163898794546405e0*scale
  fil(2)=   -0.6156141465570069496314853949e0*scale
  fil(3)=    0.2371780582153805636239247476e0*scale
  fil(4)=   -0.0822663999742123340987663521e0*scale
  fil(5)=    0.02207029188482255523789911295638968409e0*scale
  fil(6)=   -0.409765689342633823899327051188315485e-2*scale
  fil(7)=    0.45167920287502235349480037639758496e-3*scale
  fil(8)=   -0.2398228524507599670405555359023135e-4*scale
  fil(9)=    2.0904234952920365957922889447361e-6*scale
  fil(10)=  -3.7230763047369275848791496973044e-7*scale
  fil(11)=  -1.05857055496741470373494132287e-8*scale
  fil(12)=  -5.813879830282540547959250667e-11*scale
  fil(13)=   2.70800493626319438269856689037647576e-13*scale
  fil(14)=  -6.924474940639200152025730585882e-18*scale

!!$  ! second derivative filters for Daubechies 16
!!$  fil(0)=    0.e-3*scale
!!$  fil(1)=    1.e-3*scale
!!$  fil(2)=    2.e-3*scale
!!$  fil(3)=    3.e-3*scale
!!$  fil(4)=    4.e-3*scale
!!$  fil(5)=    5.e-3*scale
!!$  fil(6)=    6.e-3*scale
!!$  fil(7)=    7.e-3*scale
!!$  fil(8)=    8.e-3*scale
!!$  fil(9)=    9.e-3*scale
!!$  fil(10)=  10.e-3*scale
!!$  fil(11)=  11.e-3*scale
!!$  fil(12)=  12.e-3*scale
!!$  fil(13)=  13.e-3*scale
!!$  fil(14)=  14.e-3*scale


  do i=1,14
     fil(-i)=fil(i)
  enddo

  ekin=0.0
  
  call set_gpu_double() !after this call, all memory operations are in double precision, call set_gpu_simple() in order to have simple memory operations


  !one dimensional case
  if (ndim == 1) then
     do ndat=ndats,ndate
        do n1=n1s,n1e
           !set of one-dimensional convolutions
           !allocate arrays
           allocate(psi_in(n1,ndat,1),stat=i_stat)
         
           allocate(psi_out(ndat,n1,1),stat=i_stat)


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


         !  write(*,'(a,i6,i6)')'CPU Convolutions, dimensions:',n1,ndat
           !take timings
           !call system_clock(it0,count_rate,count_max)
           call cpu_time(t0)
           do i=1,ntimes
              call convrot_n_per(n1-1,ndat,psi_in,psi_out)
           end do
           call cpu_time(t1)
           !call system_clock(it1,count_rate,count_max)

           CPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

       !    write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
       !         CPUtime*1.d3/real(ntimes,kind=8),&
      !          real(n1*ndat*ntimes,kind=8)*32.d0/(CPUtime*1.d9)

           allocate(psi_cuda(n1,ndat,1),stat=i_stat)
    
           allocate(v_cuda(ndat,n1,1),stat=i_stat)
     

           !the input and output arrays must be reverted in this implementation
           do i=1,ndat
              do i1=1,n1
                 v_cuda(i,i1,1)=real(psi_in(i1,i,1),kind=8)
              end do
           end do

          ! write(*,'(a,i6,i6)')'ALLOCATE GPU',0,0
           call GPU_allocate(n1*ndat,psi_GPU,i_stat)
           call GPU_allocate(n1*ndat,work_GPU,i_stat)
if(i_stat == 1) then
	errorFound = 1
	write(*,*) 'Alloc error'
	return
	end if
           call GPU_send(n1*ndat,v_cuda,work_GPU,i_stat)

           !now the CUDA part
           !take timings

          ! write(*,'(a,i6,i6)')'exec kernel',0,0
         !  write(*,'(a,i6,i6)')'GPU Convolutions, dimensions:',n1,ndat

           call cpu_time(t0)
           do i=1,ntimes
              call magicfilter1d(n1-1,ndat,work_GPU,psi_GPU)
           end do
           call cpu_time(t1)
           GPUtime=real(t1-t0,kind=8)!/real(ntimes,kind=8)

         ! write(*,'(a,i6,i6)')'FINISH KERNEL',0,0
      !     write(*,'(a,f9.2,1pe12.5)')'Finished. Time(ms), GFlops',&
     !           GPUtime*1.d3/real(ntimes,kind=8),&
     !           real(n1*ndat*ntimes,kind=8)*32.d0/(GPUtime*1.d9)

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
                 comp=abs(psi_out(i,i1,1)-psi_cuda(i1,i,1))
                 !comp=abs(v_cuda(i,i1,1)-psi_cuda(i1,i,1))
                 if (comp > maxdiff) then
                    maxdiff=comp
                    i1_max=i1
                    i_max=i
                 end if
              end do
           end do

           if (maxdiff <= 3.d-7) then
              errorFound = 0 ! no error
              
         
           else
              errorFound = 1 !  error found
              
       
           end if

        
           i_all=-product(shape(psi))
           deallocate(psi,stat=i_stat)
       
           i_all=-product(shape(psi_in))
           deallocate(psi_in,stat=i_stat)
    
           i_all=-product(shape(psi_cuda))
           deallocate(psi_cuda,stat=i_stat)
     



           i_all=-product(shape(keyg))
           deallocate(keyg,stat=i_stat)
     
           i_all=-product(shape(keyv))
           deallocate(keyv,stat=i_stat)
      



        end do
     end do
    
  else 
     print *,'wrong ndim',ndim
  end if

 
end subroutine check_init

!!***
