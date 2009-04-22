!prepare the keys to be passed in GPU global or constant memory (presumably texture)
subroutine adjust_keys_for_gpu(nseg_c,nseg_f,keyv_c,keyg_c,keyv_f,keyg_f,nvctr_c,keys_GPU)
  use module_base
  implicit none
  integer, intent(in) :: nseg_c,nseg_f,nvctr_c
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f  
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(kind=8), intent(out) :: keys_GPU
  !local variables
  character(len=*), parameter :: subname='adjust_keys_for_gpu'
  !maximum number of elements to be passed into a block
  integer :: iseg,i_stat,i_all,segment_elements,nseg_tot,j,nseggpu,nblocks,nseg_block
  integer :: elems_block,nblocks_max,ncuts
  integer, dimension(:,:), allocatable :: keys
  
  call readkeysinfo(elems_block,nblocks_max)
  
  nseg_tot=0
  !assign the keys values, coarse part
  do iseg=1,nseg_c
     segment_elements=keyg_c(2,iseg)-keyg_c(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts
           nseg_tot=nseg_tot+1
        end do
     end if
  end do
  !assign the keys values, fine part
  do iseg=1,nseg_f
     segment_elements=keyg_f(2,iseg)-keyg_f(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts
           nseg_tot=nseg_tot+1
        end do
     end if
  end do

  !number of segments treated by each block;
  nseg_block =(nseg_tot+1)/nblocks_max +1

  nblocks=min(nblocks_max,nseg_tot+1)

  !print *,'nseg,nblocks',nseg_block,nblocks
  call copynseginfo(nblocks,nseg_block)

  nseggpu=nseg_block*nblocks


  !allocate the array to be copied
  allocate(keys(4,nseggpu+ndebug),stat=i_stat)
  call memocc(i_stat,keys,'keys',subname)
  
  nseg_tot=0
  !assign the keys values, coarse part
  do iseg=1,nseg_c
     segment_elements=keyg_c(2,iseg)-keyg_c(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1;
        keys(1,nseg_tot)=1 !which means coarse segment
        keys(2,nseg_tot)=segment_elements
        keys(3,nseg_tot)=keyv_c(iseg)
        keys(4,nseg_tot)=keyg_c(1,iseg)
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts-1
           nseg_tot=nseg_tot+1
           keys(1,nseg_tot)=1
           keys(2,nseg_tot)=elems_block
           keys(3,nseg_tot)=keyv_c(iseg)+(j-1)*elems_block
           keys(4,nseg_tot)=keyg_c(1,iseg)+(j-1)*elems_block
        end do
        j=ncuts
        nseg_tot=nseg_tot+1
        keys(1,nseg_tot)=1
        keys(2,nseg_tot)=segment_elements-(j-1)*elems_block
        keys(3,nseg_tot)=keyv_c(iseg)+(j-1)*elems_block
        keys(4,nseg_tot)=keyg_c(1,iseg)+(j-1)*elems_block
     end if
  end do
  !assign the keys values, fine part
  do iseg=1,nseg_f
     segment_elements=keyg_f(2,iseg)-keyg_f(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1;
        keys(1,nseg_tot)=nvctr_c !which means fine segment
        keys(2,nseg_tot)=segment_elements
        keys(3,nseg_tot)=keyv_f(iseg)
        keys(4,nseg_tot)=keyg_f(1,iseg)
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts-1
           nseg_tot=nseg_tot+1
           keys(1,nseg_tot)=nvctr_c
           keys(2,nseg_tot)=elems_block
           keys(3,nseg_tot)=keyv_f(iseg)+(j-1)*elems_block
           keys(4,nseg_tot)=keyg_f(1,iseg)+(j-1)*elems_block
        end do
        j=ncuts
        nseg_tot=nseg_tot+1
        keys(1,nseg_tot)=nvctr_c
        keys(2,nseg_tot)=segment_elements-(j-1)*elems_block
        keys(3,nseg_tot)=keyv_f(iseg)+(j-1)*elems_block
        keys(4,nseg_tot)=keyg_f(1,iseg)+(j-1)*elems_block
     end if
  end do

  !fill the last part of the segments
  do iseg=nseg_tot+1,nblocks*nseg_block
     keys(1,iseg)=0
     keys(2,iseg)=0
     keys(3,iseg)=0
     keys(4,iseg)=0
  end do

  !allocate the gpu pointer and copy the values
  call GPU_int_allocate(4*nseggpu,keys_GPU,i_stat)

  call GPU_int_send(4*nseggpu,keys,keys_GPU,i_stat)

  i_all=-product(shape(keys))*kind(keys)
  deallocate(keys,stat=i_stat)
  call memocc(i_stat,i_all,'keys',subname)


end subroutine adjust_keys_for_gpu

subroutine prepare_gpu_for_locham(n1,n2,n3,nspin,hx,hy,hz,wfd,orbs,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='prepare_gpu_for_locham'
  integer :: i_stat,iorb

  !after this call, all memory operations are in double precision, 
  !call set_gpu_simple() in order to have simple memory operations
  call set_gpu_double() 

  call adjust_keys_for_gpu(wfd%nseg_c,wfd%nseg_f,wfd%keyv(1),wfd%keyg(1,1),&
       wfd%keyv(wfd%nseg_c+1),wfd%keyg(1,wfd%nseg_c+1),wfd%nvctr_c,GPU%keys)

  !create parameters
  call creategpuparameters(n1,n2,n3,hx,hy,hz)

  !allocate the number of GPU pointers for the wavefunctions
  allocate(GPU%psi(orbs%norbp),stat=i_stat)
  call memocc(i_stat,GPU%psi,'GPU%psi',subname)

  !allocate space on the card
  !allocate the compressed wavefunctions such as to be used as workspace
  do iorb=1,orbs%norbp
     !print *,iorb
     call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2)*orbs%nspinor,GPU%psi(iorb),i_stat)
  end do
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work1,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work2,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work3,i_stat)
  !here spin value should be taken into account
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2)*nspin,GPU%rhopot,i_stat)
  !needed for the preconditioning
  call GPU_allocate((wfd%nvctr_c+7*wfd%nvctr_f),GPU%r,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%d,i_stat)

  
end subroutine prepare_gpu_for_locham



subroutine free_gpu(GPU,norbp)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='free_GPU'
  integer :: i_stat,iorb,norbp,i_all
  

  call GPU_deallocate(GPU%r,i_stat)
  call GPU_deallocate(GPU%d,i_stat)
  call GPU_deallocate(GPU%work1,i_stat)
  call GPU_deallocate(GPU%work2,i_stat)
  call GPU_deallocate(GPU%work3,i_stat)
  call GPU_deallocate(GPU%keys,i_stat)
  call GPU_deallocate(GPU%rhopot,i_stat)

  do iorb=1,norbp
     call GPU_deallocate(GPU%psi(iorb),i_stat)
  end do

  i_all=-product(shape(GPU%psi))*kind(GPU%psi)
  deallocate(GPU%psi,stat=i_stat)
  call memocc(i_stat,i_all,'GPU%psi',subname)

end subroutine free_gpu

subroutine gpu_locden(lr,nspin,hxh,hyh,hzh,orbs,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables

  call gpulocden(lr%d%n1,lr%d%n2,lr%d%n3,orbs%norbp,nspin,&
       hxh,hyh,hzh,&
       orbs%occup(min(orbs%isorb+1,orbs%norb)),&
       orbs%spinsgn(min(orbs%isorb+1,orbs%norb)),&
       GPU%psi,GPU%keys,&
       GPU%work1,GPU%work2,GPU%rhopot)

end subroutine gpu_locden


subroutine gpu_locham(n1,n2,n3,hx,hy,hz,orbs,GPU,ekin_sum,epot_sum)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  real(gp), intent(out) :: ekin_sum,epot_sum
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: i_stat,iorb
  real(gp) :: ekin,epot
  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  do iorb=1,orbs%norbp

     call gpulocham(n1,n2,n3,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
          GPU%psi(iorb),GPU%rhopot,GPU%keys,&
          GPU%work1,GPU%work2,GPU%work3,epot,ekin)

     ekin_sum=ekin_sum+orbs%occup(iorb+orbs%isorb)*ekin
     epot_sum=epot_sum+orbs%occup(iorb+orbs%isorb)*epot

  end do

end subroutine gpu_locham


subroutine gpu_precond(lr,hx,hy,hz,GPU,norbp,ncong,eval,gnrm)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: norbp,ncong
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(norbp), intent(in) :: eval
  real(wp), intent(out) :: gnrm
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: i_stat,iorb
  real(wp) :: gnrm_gpu

  gnrm=0.0_wp
  do iorb=1,norbp

     !use rhopot as a work array here
     call gpuprecond(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,&
          0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
          GPU%psi(iorb),&
          GPU%keys,GPU%r,GPU%rhopot,GPU%d,GPU%work1,GPU%work2,GPU%work3,&
          0.5_wp,ncong,gnrm_gpu)

     gnrm=gnrm+gnrm_gpu

  end do

end subroutine gpu_precond


