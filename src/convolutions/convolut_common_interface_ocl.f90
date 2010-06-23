subroutine release_acceleration_OCL(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  call ocl_clean_command_queue(GPU%queue)
  call ocl_clean(GPU%context)
end subroutine release_acceleration_OCL

subroutine init_acceleration_OCL(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU

  call ocl_create_gpu_context(GPU%context)
  !call ocl_create_command_queue(GPU%queue,GPU%context)
  call ocl_build_programs(GPU%context)
  call ocl_create_command_queue_id(GPU%queue,GPU%context,GPU%id_proc)
  call init_event_list
end subroutine init_acceleration_OCL

subroutine allocate_data_OCL(n1,n2,n3,geocode,nspin,hx,hy,hz,wfd,orbs,GPU)
  use module_base
  use module_types
  implicit none
  character(len=1), intent (in) :: geocode
  integer, intent(in) :: n1,n2,n3,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='allocate_data_OCL'
  integer :: i_stat,iorb
  integer :: n1b, n2b, n3b
  integer, dimension(3) :: periodic

  if (geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif

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

  !for preconditioner
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_r);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_r);
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_b);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_b);
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_d);
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_d);

  !full_locham stategy (always true for the moment)
  GPU%full_locham=.true.

end subroutine allocate_data_OCL


subroutine free_gpu_OCL(GPU,norbp)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='free_gpu_OCL'
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
  !for preconditioner
  call ocl_release_mem_object(GPU%psi_c_r)
  call ocl_release_mem_object(GPU%psi_f_r)
  call ocl_release_mem_object(GPU%psi_c_b)
  call ocl_release_mem_object(GPU%psi_f_b)
  call ocl_release_mem_object(GPU%psi_c_d)
  call ocl_release_mem_object(GPU%psi_f_d)

end subroutine free_gpu_OCL


subroutine local_hamiltonian_OCL(iproc,orbs,geocode,lr,hx,hy,hz,&
     nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,orbs%norbp), intent(inout) :: psi
  real(wp), dimension(lr%d%n1i,lr%d%n2i,lr%d%n3i,nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor,orbs%norbp), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_OCL'
  integer :: i_stat,iorb,isf,i
  real(gp), dimension(3) :: hgrids
  integer, dimension(3) :: periodic
  !stream ptr array
  real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr
  real(kind=8) :: stream_ptr_first_trsf
  integer :: n1, n2, n3
  real(gp) :: epot, ekin

  if (geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif

  
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

  epot_sum=0.0_gp
  ekin_sum=0.0_gp

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
          epot,ekin)
     epot_sum = epot_sum + orbs%occup(orbs%isorb+iorb)*epot
     ekin_sum = ekin_sum + orbs%occup(orbs%isorb+iorb)*ekin

     call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*orbs%nspinor*8,hpsi(1,iorb))
     call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*orbs%nspinor*8,hpsi(isf,iorb))
     
  end do
  
end subroutine local_hamiltonian_OCL

subroutine preconditionall_OCL(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,GPU)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: iproc,nproc,ncong
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  real(dp), intent(out) :: gnrm
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='preconditionall_OCL'
  integer ::  ierr,iorb,i_stat,ncplx,i_all,inds,isf
  real(wp) :: scpr
  type(GPU_pointers), intent(inout) :: GPU
  type(workarr_precond) :: w
  real(wp), dimension(:,:), allocatable :: b
  real(gp), dimension(0:7) :: scal
  !stream ptr array
  real(kind=8), dimension(orbs%norbp) :: tab_stream_ptr

  ncplx=1
  
  call allocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,lr%d,w)

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if
 
     !arrays for the CG procedure
     allocate(b(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),orbs%norbp+ndebug),stat=i_stat)
     call memocc(i_stat,b,'b',subname)

     gnrm=0.0_dp

     do iorb=1,orbs%norbp
        do inds=1,orbs%nspinor,ncplx !the streams should be more if nspinor>1
           !the nrm2 function can be replaced here by ddot
           scpr=nrm2(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),hpsi(1,inds,iorb),1)
           gnrm=gnrm+orbs%kwgts(orbs%iokpt(iorb))*scpr**2

        call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,0.5_wp,w,&
                hpsi(1,inds,iorb),b(1,iorb))

        call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*orbs%nspinor*8,&
             hpsi(1,inds,iorb))
        call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*orbs%nspinor*8,&
             hpsi(isf,inds,iorb))

        call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_b,lr%wfd%nvctr_c*orbs%nspinor*8,&
             b(1,iorb))
        call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_b,7*lr%wfd%nvctr_f*orbs%nspinor*8,&
             b(isf,iorb))

        call ocl_preconditioner(GPU%queue,&
                                (/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
                                (/0.5_gp*hx,0.5_gp*hy,0.5_gp*hz/),&
                                0.5_wp,&
                                ncong,&
                                lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
                                lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
                                GPU%psi_c,GPU%psi_f,&
                                GPU%psi_c_r,GPU%psi_f_r,&
                                GPU%psi_c_b,GPU%psi_f_b,&
                                GPU%psi_c_d,GPU%psi_f_d,&
                                GPU%d,GPU%work1,GPU%work2,GPU%work3)

        call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*orbs%nspinor*8,hpsi(1,inds,iorb))
        call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*orbs%nspinor*8,hpsi(isf,inds,iorb))

        end do
     end do


     !end of dynamic repartition
 

  i_all=-product(shape(b))*kind(b)
  deallocate(b,stat=i_stat)
  call memocc(i_stat,i_all,'b',subname)

  call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)


end subroutine preconditionall_OCL

subroutine local_partial_density_OCL(iproc,nproc,orbs,&
     nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nrhotot
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
 
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  real(dp), dimension(lr%d%n1i,lr%d%n2i,nrhotot,max(nspin,orbs%nspinor)), intent(inout) :: rho_p
  type(GPU_pointers), intent(inout) :: GPU
  
  integer:: iorb,i_stat,isf,iaddjmp
  real(kind=8) :: stream_ptr
  real(gp) :: hfac


  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

  call set_d(GPU%queue, lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin , 0.0d0,  GPU%rhopot)
  !copy the wavefunctions on GPU
  do iorb=1,orbs%norbp
     
    call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*orbs%nspinor*8,&
             psi(1,(iorb-1)*orbs%nspinor+1)) 
    call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*orbs%nspinor*8,&
             psi(isf,(iorb-1)*orbs%nspinor+1))
 
    hfac=orbs%occup(min(orbs%isorb+1,orbs%norb)+iorb-1)/(hxh*hyh*hzh);
    if (orbs%spinsgn(min(orbs%isorb+1,orbs%norb)+iorb-1) > 0.0) then
      iaddjmp = 0
    else
      iaddjmp = 8*(lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1)
    endif
  !calculate the density
   call ocl_locden(GPU%queue, (/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
                          hfac, iaddjmp,&
                          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
                          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
                          GPU%psi_c,GPU%psi_f,&
                          GPU%work1,GPU%work2,GPU%work3,&
                          GPU%rhopot)

  
  end do
  !copy back the results and leave the uncompressed wavefunctions on the card
  

  call ocl_enqueue_read_buffer(GPU%queue,GPU%rhopot,lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin*8,rho_p)

end subroutine local_partial_density_OCL


