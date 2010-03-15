subroutine allocate_data_OCL(n1,n2,n3,periodic,nspin,hx,hy,hz,wfd,orbs,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  integer, dimension(3), intent(in) :: periodic
  !local variables
  character(len=*), parameter :: subname='prepare_gpu_for_locham'
  integer :: i_stat,iorb
  integer :: n1b, n2b, n3b

  n1b = (n1+1) * 2
  n2b = (n2+1) * 2
  n3b = (n3+1) * 2
  if (periodic(1)==0) then
    n1b = n1b + 2*7 + 15
  endif
  if (periodic(2)==0) then
    n2b = n2b + 2*7 + 15
  endif
  if (periodic(3)==0) then
    n3b = n3b + 2*7 + 15
  endif


  !allocate the number of GPU pointers for the wavefunctions
  !allocate(GPU%psi(orbs%norbp+ndebug),stat=i_stat)
  !call memocc(i_stat,GPU%psi,'GPU%psi',subname)

  !allocate space on the card
  !allocate the compressed wavefunctions such as to be used as workspace
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f);
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work1)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work2)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work3)
  !here spin value should be taken into account
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*nspin*8,GPU%rhopot)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%d)

  !allocate and copy the compression-decompression keys
  call ocl_create_read_buffer(GPU%context,wfd%nseg_c*4*2,GPU%keyg_c)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_c*4,GPU%keyv_c)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_f*4*2,GPU%keyg_f)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_f*4,GPU%keyv_f)
  call ocl_enqueue_write_buffer(GPU%queue,GPU%keyg_c,wfd%nseg_c*2*4,wfd%keyg)
  call ocl_enqueue_write_buffer(GPU%queue,GPU%keyv_c,wfd%nseg_c*4,wfd%keyv)
  if (wfd%nseg_f > 0) then
     call ocl_enqueue_write_buffer(GPU%queue,GPU%keyg_f,wfd%nseg_f*2*4,wfd%keyg(1,wfd%nseg_c+1))
     call ocl_enqueue_write_buffer(GPU%queue,GPU%keyv_f,wfd%nseg_f*4,wfd%keyv(wfd%nseg_c+1))
  end if

end subroutine allocate_data_OCL


subroutine free_gpu_OCL(GPU,norbp)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='free_GPU'
  integer :: i_stat,iorb,norbp,i_all
  

  !call ocl_release_mem_object(GPU%r)
  call ocl_release_mem_object(GPU%d)
  call ocl_release_mem_object(GPU%work1)
  call ocl_release_mem_object(GPU%work2)
  call ocl_release_mem_object(GPU%work3)
  call ocl_release_mem_object(GPU%rhopot)
  call ocl_release_mem_object(GPU%keyg_c)
  call ocl_release_mem_object(GPU%keyv_c)
  call ocl_release_mem_object(GPU%keyg_f)
  call ocl_release_mem_object(GPU%keyv_f)
  call ocl_release_mem_object(GPU%psi_c)
  call ocl_release_mem_object(GPU%psi_f)

end subroutine free_gpu_OCL


subroutine local_hamiltonian_OCL(iproc,orbs,periodic,lr,hx,hy,hz,&
     nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nspin
  integer, dimension(3), intent(in) :: periodic
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,orbs%norbp), intent(inout) :: psi
  real(wp), dimension(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,orbs%norbp), intent(out) :: hpsi

  type(GPU_pointers), intent(inout) :: GPU
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_GPU'
  integer :: i_stat,iorb,isf,i
  real(gp), dimension(3) :: hgrids
  !stream ptr array
  real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr
  real(kind=8) :: stream_ptr_first_trsf
  integer :: n1, n2, n3
  
  n1 = (lr%d%n1+1) * 2
  n2 = (lr%d%n2+1) * 2
  n3 = (lr%d%n3+1) * 2
  if (periodic(1)==0) then
    n1 = n1 + 2*7 + 15
  endif
  if (periodic(2)==0) then
    n2 = n2 + 2*7 + 15
  endif
  if (periodic(3)==0) then
    n3 = n3 + 2*7 + 15
  endif

  call ocl_enqueue_write_buffer(GPU%queue,GPU%rhopot,n1*n2*n3*8,pot) 
 
  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

  hgrids(1)=0.5_gp*hx
  hgrids(2)=0.5_gp*hy
  hgrids(3)=0.5_gp*hz


  do iorb=1,orbs%norbp

     !if orbs%nspinor /= 1 this implementation should be rediscussed
     if (GPU%full_locham) then
        call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*orbs%nspinor*8,&
             psi(1,iorb)) 
        call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*orbs%nspinor*8,&
             psi(isf,iorb))
     else
        stop 'ONLY FULL LOCHAM IS IMPLEMENTED!'
     end if
     
     !calculate the local hamiltonian
     !WARNING: the difference between full_locham and normal locham is inside
     call ocl_fulllocham_generic(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          (/periodic(1),periodic(2),periodic(3)/),&
          hgrids,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,& 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,& 
          GPU%psi_c,GPU%psi_f,GPU%rhopot, &
          GPU%work1,GPU%work2,GPU%work3,GPU%d,&
          epot_sum,ekin_sum)

     call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*orbs%nspinor*8,hpsi(1,iorb))
     call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*orbs%nspinor*8,hpsi(isf,iorb))
     
  end do
  
end subroutine local_hamiltonian_OCL
