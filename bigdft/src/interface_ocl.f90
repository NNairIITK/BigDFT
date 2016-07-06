!> @file
!!  Interface routines to do GPU convolution with OpenCL
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine release_acceleration_OCL(GPU)
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  call ocl_clean_command_queue(GPU%queue)
  call ocl_clean(GPU%context)
END SUBROUTINE release_acceleration_OCL


subroutine init_acceleration_OCL(matacc,GPU)
  use module_input_keys, only: material_acceleration
  use module_types
  implicit none
  type(material_acceleration), intent(in) :: matacc
  type(GPU_pointers), intent(out) :: GPU
  integer(kind=8) :: context_address

  call ocl_create_context(GPU%context, matacc%OCL_platform, matacc%OCL_devices, matacc%iacceleration,&
                     GPU%ndevices)
  !call ocl_create_gpu_context(GPU%context,GPU%ndevices)
  !call ocl_create_command_queue(GPU%queue,GPU%context)
  !to avoid a representation of the address which is lower than tiny(1.d0)
  context_address=transfer(GPU%context,context_address)
  !if (GPU%context /= 0.) then
  if (context_address /= int(0,kind=8)) then
     call ocl_build_programs(GPU%context)
     call ocl_create_command_queue_id(GPU%queue,GPU%context,GPU%id_proc)
     call init_event_list(GPU%context)
  end if
END SUBROUTINE init_acceleration_OCL


subroutine allocate_data_OCL(n1,n2,n3,geocode,nspin,wfd,orbs,GPU)
  use module_base
  use module_types
  implicit none
  character(len=1), intent (in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: n1,n2,n3,nspin
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='allocate_data_OCL'
  logical, parameter :: pin=.false.
  integer :: n1b, n2b, n3b, iorb,ispinor
  integer, dimension(3) :: periodic

  call f_routine(id=subname)

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
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c)
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work1)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work2)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work3)
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%d)

  if ( orbs%nspinor == 2) then
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_i)
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work1_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work2_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%work3_i)
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%d_i)
  end if
  !here spin value should be taken into account
  call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%rhopot_up)
  if( nspin == 2 ) then
    call ocl_create_read_write_buffer(GPU%context, n1b*n2b*n3b*8,GPU%rhopot_down)
  end if

  !allocate and copy the compression-decompression keys
  call ocl_create_read_buffer(GPU%context,wfd%nseg_c*4*2,GPU%keyg_c)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_c*4,GPU%keyv_c)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_f*4*2,GPU%keyg_f)
  call ocl_create_read_buffer(GPU%context,wfd%nseg_f*4,GPU%keyv_f)
  
  if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,wfd%nseg_c*4*2,wfd%keygloc,GPU%keyg_c_host)
  call ocl_enqueue_write_buffer(GPU%queue,GPU%keyg_c,wfd%nseg_c*2*4,wfd%keygloc)
  if (pin) call ocl_release_mem_object(GPU%keyg_c_host)

  if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,wfd%nseg_c*4,wfd%keyvloc,GPU%keyv_c_host)
  call ocl_enqueue_write_buffer(GPU%queue,GPU%keyv_c,wfd%nseg_c*4,wfd%keyvloc)
  if (pin) call ocl_release_mem_object(GPU%keyv_c_host)

  if (wfd%nseg_f > 0) then
     if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,wfd%nseg_f*4*2,wfd%keygloc(1,wfd%nseg_c+1),GPU%keyg_f_host)
     call ocl_enqueue_write_buffer(GPU%queue,GPU%keyg_f,wfd%nseg_f*2*4,wfd%keygloc(1,wfd%nseg_c+1))
     if (pin) call ocl_release_mem_object(GPU%keyg_f_host)
  
     if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,wfd%nseg_f*4,wfd%keyvloc(wfd%nseg_c+1),GPU%keyv_f_host)
     call ocl_enqueue_write_buffer(GPU%queue,GPU%keyv_f,wfd%nseg_f*4,wfd%keyvloc(wfd%nseg_c+1))
     if (pin) call ocl_release_mem_object(GPU%keyv_f_host)
  end if

  !for preconditioner
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_r)
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_r)
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_b)
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_b)
  call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_d)
  call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_d)
  if ( orbs%nspinor == 2) then
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_r_i)
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_r_i)
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_b_i)
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_b_i)
    call ocl_create_read_write_buffer(GPU%context,wfd%nvctr_c*8,GPU%psi_c_d_i)
    call ocl_create_read_write_buffer(GPU%context,7*wfd%nvctr_f*8,GPU%psi_f_d_i)
  end if
  !full_locham stategy (always true for the moment)
  GPU%full_locham=.true.

  GPU%ekin=f_malloc_ptr((/2,orbs%norbp/),id='ekin')
  GPU%epot=f_malloc_ptr((/2,orbs%norbp/),id='epot')
  if (pin) then
     GPU%ekinpot_host=f_malloc_ptr((/orbs%nspinor,orbs%norbp,2/),id='ekinpot_host')
     GPU%psicf_host=f_malloc_ptr((/2*orbs%nspinor,orbs%norbp/),id='psicf_host')
     GPU%hpsicf_host=f_malloc_ptr((/2*orbs%nspinor,orbs%norbp/),id='hpsicf_host')
     GPU%bprecond_host=f_malloc_ptr(2*orbs%nspinor,id='bprecond_host')
  end if

  !pin the memory of the orbitals energies
  if (pin) then
     do iorb=1,orbs%norbp
        do ispinor=1,orbs%nspinor
           call ocl_pin_write_buffer_async(GPU%context,GPU%queue,8,GPU%ekin(ispinor,iorb),GPU%ekinpot_host(ispinor,iorb,1))
           call ocl_pin_write_buffer_async(GPU%context,GPU%queue,8,GPU%epot(ispinor,iorb),GPU%ekinpot_host(ispinor,iorb,2))
        end do
     end do
  end if

  nullify(GPU%hpsi_ASYNC)

  call f_release_routine()

END SUBROUTINE allocate_data_OCL


subroutine free_gpu_OCL(GPU,orbs,nspin)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  character(len=*), parameter :: subname='free_gpu_OCL'
  logical, parameter :: pin=.false.
  integer :: iorb,ispinor

  call f_free_ptr(GPU%ekin)
  call f_free_ptr(GPU%epot)


  call ocl_release_mem_object(GPU%d)
  call ocl_release_mem_object(GPU%work1)
  call ocl_release_mem_object(GPU%work2)
  call ocl_release_mem_object(GPU%work3)
  call ocl_release_mem_object(GPU%rhopot_up)
  if ( nspin == 2 ) then
    call ocl_release_mem_object(GPU%rhopot_down)
  endif
  call ocl_release_mem_object(GPU%keyg_c)
  call ocl_release_mem_object(GPU%keyv_c)
  call ocl_release_mem_object(GPU%keyg_f)
  call ocl_release_mem_object(GPU%keyv_f)
  call ocl_release_mem_object(GPU%psi_c)
  call ocl_release_mem_object(GPU%psi_f)
  if ( orbs%nspinor == 2) then
    call ocl_release_mem_object(GPU%psi_c_i)
    call ocl_release_mem_object(GPU%psi_f_i)
    call ocl_release_mem_object(GPU%work1_i)
    call ocl_release_mem_object(GPU%work2_i)
    call ocl_release_mem_object(GPU%work3_i)
    call ocl_release_mem_object(GPU%d_i)
  endif
  !for preconditioner
  call ocl_release_mem_object(GPU%psi_c_r)
  call ocl_release_mem_object(GPU%psi_f_r)
  call ocl_release_mem_object(GPU%psi_c_b)
  call ocl_release_mem_object(GPU%psi_f_b)
  call ocl_release_mem_object(GPU%psi_c_d)
  call ocl_release_mem_object(GPU%psi_f_d)
  if ( orbs%nspinor == 2) then
    call ocl_release_mem_object(GPU%psi_c_r_i)
    call ocl_release_mem_object(GPU%psi_f_r_i)
    call ocl_release_mem_object(GPU%psi_c_b_i)
    call ocl_release_mem_object(GPU%psi_f_b_i)
    call ocl_release_mem_object(GPU%psi_c_d_i)
    call ocl_release_mem_object(GPU%psi_f_d_i)
  endif

  if(associated(GPU%hpsi_ASYNC)) nullify(GPU%hpsi_ASYNC)

  if (pin) then
     !for pinning tracing
     do iorb=1,orbs%norbp
        do ispinor=1,orbs%nspinor
           call ocl_release_mem_object(GPU%ekinpot_host(ispinor,iorb,1))
           call ocl_release_mem_object(GPU%ekinpot_host(ispinor,iorb,2))
        end do
     end do
     call f_free_ptr(GPU%ekinpot_host)
     call f_free_ptr(GPU%psicf_host)
     call f_free_ptr(GPU%hpsicf_host)
     call f_free_ptr(GPU%bprecond_host)

  end if


END SUBROUTINE free_gpu_OCL


subroutine daub_to_isf_OCL(lr,psi,psi_r,GPU)
  use module_base
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  type(GPU_pointers), intent(inout) :: GPU
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psi_r
  
  integer, dimension(3) :: periodic
  integer :: isf

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif 

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if


  call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1))
  call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf))

  call ocl_daub_to_isf(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          periodic,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
          GPU%psi_c,GPU%psi_f,&
          GPU%work1,GPU%work2,GPU%work3,GPU%d)

  call ocl_enqueue_read_buffer(GPU%queue,GPU%work2,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,&
          psi_r(1))

END SUBROUTINE daub_to_isf_OCL


subroutine isf_to_daub_OCL(lr,psi_r,psi,GPU)
  use module_base
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  type(GPU_pointers), intent(inout) :: GPU
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)), intent(out) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(in) :: psi_r
  
  integer, dimension(3) :: periodic
  integer :: isf

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif 

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if


  call ocl_enqueue_write_buffer(GPU%queue,GPU%work1,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,&
          psi_r(1))

  call ocl_isf_to_daub(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          periodic,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
          GPU%psi_c,GPU%psi_f,&
          GPU%work1,GPU%work2,GPU%work3,GPU%d)

  call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1))
  call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf))

END SUBROUTINE isf_to_daub_OCL


subroutine local_hamiltonian_OCL(orbs,lr,hx,hy,hz,&
     nspin,pot,psi,hpsi,ekin_sum,epot_sum,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
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
  logical, parameter :: pin=.false.
  integer :: iorb,isf
  real(gp), dimension(3) :: hgrids
  integer, dimension(3) :: periodic
  !stream ptr array
  real(kind=8) :: rhopot
  integer :: n1, n2, n3

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (lr%geocode /= 'F') then
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

  !define the pinned adresses for the pinning of the interesting objects

  if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,n1*n2*n3*8,pot,GPU%rhopot_up_host)
  call ocl_enqueue_write_buffer_async(GPU%queue,GPU%rhopot_up,n1*n2*n3*8,pot) 
  if (pin) call ocl_release_mem_object(GPU%rhopot_up_host)
  if( nspin == 2 ) then
     if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,n1*n2*n3*8,pot(1,1,1,2),GPU%rhopot_down_host)
     call ocl_enqueue_write_buffer_async(GPU%queue,GPU%rhopot_down,n1*n2*n3*8,pot(1,1,1,2)) 
     if (pin) call ocl_release_mem_object(GPU%rhopot_down_host)
  end if
 
  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

  hgrids(1)=0.5_gp*hx
  hgrids(2)=0.5_gp*hy
  hgrids(3)=0.5_gp*hz

!!$  epot_sum=0.0_gp
!!$  ekin_sum=0.0_gp
  do iorb=1,orbs%norbp
     if (orbs%nspinor == 2) then
     end if
  enddo

  if (pin .and. orbs%norbp > 0) call ocl_create_write_buffer_host( GPU%context, &
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8, hpsi, GPU%hpsicf_host(1,1) )

!  call ocl_create_read_buffer_host( GPU%context, &
!       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8, psi, GPU%psicf_host(1,1) )

  do iorb=1,orbs%norbp

     if (orbs%spinsgn(orbs%isorb+iorb) > 0.0) then
       rhopot = GPU%rhopot_up
     else
       rhopot = GPU%rhopot_down
     endif
     !if orbs%nspinor /= 1 this implementation should be rediscussed
     if (.not. GPU%full_locham) then
        stop 'ONLY FULL LOCHAM IS IMPLEMENTED!'
     end if

     !pin the adresses of the wavefucntions
     if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,psi(1,iorb),GPU%psicf_host(1,iorb))
!     call ocl_map_write_buffer_async(GPU%queue, GPU%psicf_host(1,1), &
!          (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*(iorb-1)*orbs%nspinor, lr%wfd%nvctr_c*8)
     call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1,iorb))
!     call ocl_unmap_mem_object(GPU%queue, GPU%psicf_host(1,1), psi(1,iorb))
     if (pin) call ocl_release_mem_object(GPU%psicf_host(1,iorb))

     if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,psi(isf,iorb),GPU%psicf_host(2,iorb))
!     call ocl_map_write_buffer_async(GPU%queue, GPU%psicf_host(1,1), &
!          (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*(iorb-1)*orbs%nspinor+lr%wfd%nvctr_c*8, 7*lr%wfd%nvctr_f*8)
     call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf,iorb))
!     call ocl_unmap_mem_object(GPU%queue, GPU%psicf_host(1,1), psi(isf,iorb))
     if (pin) call ocl_release_mem_object(GPU%psicf_host(2,iorb))
     if (orbs%nspinor == 2) then
        if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
             psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb),GPU%psicf_host(3,iorb))
!       call ocl_map_write_buffer_async(GPU%queue, GPU%psicf_host(1,1), &
!            (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*((iorb-1)*orbs%nspinor+1), lr%wfd%nvctr_c*8)
       call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,&
            psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb))
!       call ocl_unmap_mem_object(GPU%queue, GPU%psicf_host(1,1), psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb))
       if (pin) call ocl_release_mem_object(GPU%psicf_host(3,iorb))

       if (pin) call ocl_pin_read_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
            psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb),GPU%psicf_host(4,iorb))
!       call ocl_map_write_buffer_async(GPU%queue, GPU%psicf_host(1,1), &
!            (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*((iorb-1)*orbs%nspinor+1)+lr%wfd%nvctr_c*8, 7*lr%wfd%nvctr_f*8)
       call ocl_enqueue_write_buffer_async(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,&
            psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb))
!       call ocl_unmap_mem_object(GPU%queue, GPU%psicf_host(1,1), psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb))
       if (pin) call ocl_release_mem_object(GPU%psicf_host(4,iorb))
     end if
     !calculate the local hamiltonian
     !WARNING: the difference between full_locham and normal locham is inside
     call ocl_fulllocham_generic_k(GPU%queue,(/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          periodic,&
          hgrids,&
          (/orbs%kpts(1,orbs%iokpt(iorb)),orbs%kpts(2,orbs%iokpt(iorb)),orbs%kpts(3,orbs%iokpt(iorb))/),&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
          GPU%psi_c,GPU%psi_f,&
          GPU%psi_c_i,GPU%psi_f_i,&
          rhopot,&
          GPU%work1,GPU%work2,GPU%work3,&
          GPU%work1_i,GPU%work2_i,GPU%work3_i,&
          GPU%d,&
          GPU%d_i,&
          orbs%nspinor,&
          GPU%epot(1,iorb),GPU%ekin(1,iorb))

!     call ocl_pin_write_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,hpsi(1,iorb),GPU%hpsicf_host(1,iorb))
     if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1), &
          (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*(iorb-1)*orbs%nspinor, lr%wfd%nvctr_c*8)
     call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,hpsi(1,iorb))
     if (pin) call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(1,iorb))
!     call ocl_release_mem_object(GPU%hpsicf_host(1,iorb))
!     call ocl_pin_write_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,hpsi(isf,iorb),GPU%hpsicf_host(2,iorb))
     if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1), &
          (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*(iorb-1)*orbs%nspinor+lr%wfd%nvctr_c*8, 7*lr%wfd%nvctr_f*8)
     call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,hpsi(isf,iorb))
     if (pin) call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(isf,iorb))
!     call ocl_release_mem_object(GPU%hpsicf_host(2,iorb))

     if (orbs%nspinor == 2) then
!       call ocl_pin_write_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
!             hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb),GPU%hpsicf_host(3,iorb))
        if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1), &
            (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*((iorb-1)*orbs%nspinor+1), lr%wfd%nvctr_c*8)
       call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb))
       if (pin) call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+1,iorb))
!       call ocl_release_mem_object(GPU%hpsicf_host(3,iorb))
!       call ocl_pin_write_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
!            hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb),GPU%hpsicf_host(4,iorb))
       if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1), &
            (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*((iorb-1)*orbs%nspinor+1)+lr%wfd%nvctr_c*8, 7*lr%wfd%nvctr_f*8)
       call ocl_enqueue_read_buffer_async(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb))
       if (pin) call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+isf,iorb))
!       call ocl_release_mem_object(GPU%hpsicf_host(4,iorb))
     end if
  end do
!  call ocl_release_mem_object(GPU%psicf_host(1,1))
  if (pin .and. orbs%norbp > 0) call ocl_release_mem_object(GPU%hpsicf_host(1,1))
  if (.not. ASYNCconv) then
     call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU)
  endif

END SUBROUTINE local_hamiltonian_OCL


subroutine finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  real(gp), intent(out) :: ekin_sum,epot_sum
  type(GPU_pointers), intent(inout) :: GPU

  integer :: iorb

  call ocl_finish(GPU%queue)
  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  do iorb=1,orbs%norbp
    ekin_sum = ekin_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*((GPU%ekin(1,iorb)+GPU%ekin(2,iorb))&
                 - (GPU%epot(1,iorb)+GPU%epot(2,iorb)))
    epot_sum = epot_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*(GPU%epot(1,iorb)+GPU%epot(2,iorb))
  end do

  !free pinning information for wavefunctions
!  do iorb=1,orbs%norbp
!     do ispinor=1,orbs%nspinor
!        call ocl_release_mem_object(GPU%psicf_host(1+(ispinor-1)*2,iorb))
!        call ocl_release_mem_object(GPU%psicf_host(2+(ispinor-1)*2,iorb))
!        call ocl_release_mem_object(GPU%hpsicf_host(1+(ispinor-1)*2,iorb))
!        call ocl_release_mem_object(GPU%hpsicf_host(2+(ispinor-1)*2,iorb))
!     end do
!  end do
  !free pinning information for potential

END SUBROUTINE finish_hamiltonian_OCL


subroutine preconditionall_OCL(orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero,GPU)
  use module_base
  use module_types
  use locreg_operations, only: workarr_precond,allocate_work_arrays,deallocate_work_arrays
  implicit none
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: ncong
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  real(dp), intent(out) :: gnrm,gnrm_zero
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='preconditionall_OCL'
  logical, parameter :: pin=.false.
  integer ::  iorb,jorb,ncplx,inds,isf,ikpt,ispinor,ioff_c,ioff_f,ioff_ci,ioff_fi
  real(wp) :: scpr
  real(gp) :: cprecr,eval_zero,evalmax
  type(GPU_pointers), intent(inout) :: GPU
  type(workarr_precond) :: w
  integer, dimension(3) :: periodic
  real(wp), dimension(:,:), allocatable :: b
  real(gp), dimension(0:7) :: scal
  !stream ptr array

  !the eval array contains all the values
  !take the max for all k-points
  !one may think to take the max per k-point

  call f_routine(id=subname)

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif


  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

     !arrays for the CG procedure
!!$     allocate(b(orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),orbs%norbp+ndebug),stat=i_stat)
!!$     call memocc(i_stat,b,'b',subname)
  b=f_malloc((/orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1/),id='b')

!!$     allocate(b(orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1+ndebug),stat=i_stat)
!!$     call memocc(i_stat,b,'b',subname)

  if (pin) then
     call ocl_pin_read_buffer(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,b(1,1),GPU%bprecond_host(1))
     call ocl_pin_read_buffer(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,b(isf,1),GPU%bprecond_host(2))
     if(orbs%nspinor == 2) then
        call ocl_pin_read_buffer(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
             b(1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1),GPU%bprecond_host(3))
        call ocl_pin_read_buffer(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
             b(isf+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1),GPU%bprecond_host(4))
     end if

     if (orbs%norbp >0) then
        call ocl_create_write_buffer_host(GPU%context, &
             orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8, hpsi, GPU%hpsicf_host(1,1))
        !call ocl_create_read_buffer_host(GPU%context, &
        !     orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8, hpsi, GPU%hpsicf_host(1,1))
     end if

  end if

     gnrm=0.0_dp
     gnrm_zero=0.0_dp
  call allocate_work_arrays(lr%geocode,lr%hybrid_on,orbs%nspinor,lr%d,w)
  if (orbs%norbp >0) ikpt=orbs%iokpt(1)
  do iorb=1,orbs%norbp
     !if it is the first orbital or the k-point has changed calculate the max
     if (orbs%iokpt(iorb) /= ikpt .or. iorb == 1) then
        !the eval array contains all the values
        !take the max for all k-points
        !one may think to take the max per k-point
        evalmax=orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+1)
        do jorb=1,orbs%norb
           evalmax=max(orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+jorb),evalmax)
        enddo
        eval_zero=evalmax
        ikpt=orbs%iokpt(iorb)
     end if


       if (orbs%kpts(1,orbs%iokpt(iorb))**2+orbs%kpts(2,orbs%iokpt(iorb))**2+&
           orbs%kpts(3,orbs%iokpt(iorb))**2 > 0.0_gp .or. orbs%nspinor==2 ) then
         ncplx=2
       else
         ncplx=1
       end if

        call cprecr_from_eval(lr%geocode,eval_zero,orbs%eval(orbs%isorb+iorb),cprecr)

        do inds=1,orbs%nspinor,ncplx !the streams should be more if nspinor>1

           !offsets for pinning
           ioff_c=8*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*(orbs%nspinor*(iorb-1)+(inds-1))
           ioff_f=ioff_c+8*lr%wfd%nvctr_c
           !case for ncplx=2
           ioff_ci=ioff_f+7*lr%wfd%nvctr_f*8
           ioff_fi=ioff_ci+8*lr%wfd%nvctr_c

           !the nrm2 function can be replaced here by ddot
           scpr=nrm2(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),hpsi(1,inds,iorb),1)
           if (orbs%occup(orbs%isorb+iorb) == 0.0_gp) then
              gnrm_zero=gnrm_zero+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
           else
              !write(17,*)'iorb,gnrm',orbs%isorb+iorb,scpr**2
              gnrm=gnrm+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
           end if
           call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,&
                hpsi(1,inds,iorb),b(1,1))!iorb))

           if (pin) then
              call ocl_pin_read_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
                   hpsi(1,inds,iorb),GPU%psicf_host(1,iorb))
              !call ocl_pin_read_write_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
              !  hpsi(1,inds,iorb),GPU%hpsicf_host(1,iorb))

              !call ocl_pin_read_write_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
              !  hpsi(1,inds,iorb),GPU%hpsicf_host(1,iorb))
              !call ocl_map_read_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_c,lr%wfd%nvctr_c*8)
              !call ocl_map_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_c,lr%wfd%nvctr_c*8)
           end if
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
                hpsi(1,inds,iorb))
           if (pin) then
              !call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(1,inds,iorb))
              call ocl_release_mem_object(GPU%psicf_host(1,iorb))
           end if

           if (pin) then
              call ocl_pin_read_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
                hpsi(isf,inds,iorb),GPU%psicf_host(2,iorb))
              !call ocl_pin_read_write_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
              !  hpsi(isf,inds,iorb),GPU%hpsicf_host(2,iorb))
              !call ocl_map_read_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_f,7*lr%wfd%nvctr_f*8)
              !call ocl_map_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_f,7*lr%wfd%nvctr_f*8)
           end if
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
                hpsi(isf,inds,iorb))
           if (pin) then
              !call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(isf,inds,iorb))
              call ocl_release_mem_object(GPU%psicf_host(2,iorb))
           end if

           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_b,lr%wfd%nvctr_c*8,&
                b(1,1))!iorb))
           call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_b,7*lr%wfd%nvctr_f*8,&
                b(isf,1))!iorb))

           if(ncplx == 2) then
              if (pin) then
                 call ocl_pin_read_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
                      hpsi(1,inds+1,iorb),GPU%psicf_host(3,iorb))
                 !call ocl_pin_read_write_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,&
                 !  hpsi(1,inds+1,iorb),GPU%hpsicf_host(3,iorb))
                 !call ocl_map_read_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_ci,lr%wfd%nvctr_c*8)
                 !call ocl_map_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_ci,lr%wfd%nvctr_c*8)
              end if
              call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,&
                   hpsi(1,inds+1,iorb))
              if (pin) then
                 !call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(1,inds+1,iorb))
                 call ocl_release_mem_object(GPU%psicf_host(3,iorb))
              end if

              if (pin) then
                 call ocl_pin_read_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
                      hpsi(isf,inds+1,iorb),GPU%psicf_host(4,iorb))
                 !call ocl_pin_read_write_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,&
                 !  hpsi(isf,inds+1,iorb),GPU%hpsicf_host(4,iorb))
                 !call ocl_map_read_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_fi, 7*lr%wfd%nvctr_f*8)
                 !call ocl_map_write_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_fi, 7*lr%wfd%nvctr_f*8)
              end if
              call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,&
                   hpsi(isf,inds+1,iorb))
              if (pin) then
                 !call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(isf,inds+1,iorb))
                 call ocl_release_mem_object(GPU%psicf_host(4,iorb))
              end if

             call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c_b_i,lr%wfd%nvctr_c*8,&
                  b(1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1))!iorb))
             call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f_b_i,7*lr%wfd%nvctr_f*8,&
                  b(isf+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),1))!iorb))
           endif
           call ocl_preconditioner_generic_k(GPU%queue,&
                (/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
                periodic,&
                (/0.5_gp*hx,0.5_gp*hy,0.5_gp*hz/),&
                (/orbs%kpts(1,orbs%iokpt(iorb)),orbs%kpts(2,orbs%iokpt(iorb)),orbs%kpts(3,orbs%iokpt(iorb))/),&
                cprecr,&
                ncong,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
                lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
                GPU%psi_c,GPU%psi_f,&
                GPU%psi_c_i,GPU%psi_f_i,&
                GPU%psi_c_r,GPU%psi_f_r,&
                GPU%psi_c_r_i,GPU%psi_f_r_i,&
                GPU%psi_c_b,GPU%psi_f_b,&
                GPU%psi_c_b_i,GPU%psi_f_b_i,&
                GPU%psi_c_d,GPU%psi_f_d,&
                GPU%psi_c_d_i,GPU%psi_f_d_i,&
                GPU%d,GPU%work1,GPU%work2,GPU%work3,&
                GPU%d_i,GPU%work1_i,GPU%work2_i,GPU%work3_i,&
                ncplx,GPU%ekin) !buffer for scalars resulting from reductions

           if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_c,lr%wfd%nvctr_c*8)
           call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,hpsi(1,inds,iorb))
           if (pin) then
              !call ocl_release_mem_object(GPU%hpsicf_host(1,iorb))
              call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(1,inds,iorb))
           end if

           if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_f,7*lr%wfd%nvctr_f*8)
           call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,hpsi(isf,inds,iorb))
           if (pin) then
              !call ocl_release_mem_object(GPU%hpsicf_host(2,iorb))
              call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(isf,inds,iorb))
           end if

           if ( ncplx == 2 ) then
              if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_ci,lr%wfd%nvctr_c*8)
              call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_c_i,lr%wfd%nvctr_c*8,hpsi(1,inds+1,iorb))
              if (pin) then
                 !call ocl_release_mem_object(GPU%hpsicf_host(3,iorb))
                 call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(1,inds+1,iorb))
             end if

             if (pin) call ocl_map_read_buffer_async(GPU%queue, GPU%hpsicf_host(1,1),ioff_fi, 7*lr%wfd%nvctr_f*8)
             call ocl_enqueue_read_buffer(GPU%queue,GPU%psi_f_i,7*lr%wfd%nvctr_f*8,hpsi(isf,inds+1,iorb))
             if (pin) then
                !call ocl_release_mem_object(GPU%hpsicf_host(4,iorb))
                call ocl_unmap_mem_object(GPU%queue, GPU%hpsicf_host(1,1), hpsi(isf,inds+1,iorb))
             end if
           endif

        end do
     end do

     if (pin .and. orbs%norbp > 0) then
        do ispinor=1,orbs%nspinor
           call ocl_release_mem_object(GPU%bprecond_host(1+(ispinor-1)*2))
           call ocl_release_mem_object(GPU%bprecond_host(2+(ispinor-1)*2))
        end do
        call ocl_release_mem_object(GPU%hpsicf_host(1,1))
     end if

     call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)

!!$  i_all=-product(shape(b))*kind(b)
!!$  deallocate(b,stat=i_stat)
!!$  call memocc(i_stat,i_all,'b',subname)

  call f_free(b)

  call f_release_routine()

END SUBROUTINE preconditionall_OCL


subroutine local_partial_density_OCL(orbs,&
     nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nrhotot
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  real(dp), dimension(lr%d%n1i,lr%d%n2i,nrhotot,nspin), intent(inout) :: rho_p
  type(GPU_pointers), intent(inout) :: GPU
  !local variables
  logical, parameter :: pin=.false.
  integer:: iorb,iorb_r,isf,ispinor
  real(gp) :: hfac
  integer, dimension(3) :: periodic
  real(kind=8) :: rhopot

  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif

  if (lr%wfd%nvctr_f > 0) then
     isf=lr%wfd%nvctr_c+1
  else
     isf=lr%wfd%nvctr_c
  end if

  call set_d(GPU%queue, lr%d%n1i*lr%d%n2i*lr%d%n3i , 1.d-20,  GPU%rhopot_up)
  if ( nspin == 2 ) then
    call set_d(GPU%queue, lr%d%n1i*lr%d%n2i*lr%d%n3i , 1.d-20,  GPU%rhopot_down)
  end if
  if (pin .and. orbs%norbp > 0) call ocl_create_read_buffer_host( GPU%context, &
       orbs%norbp*orbs%nspinor*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8, psi, GPU%psicf_host(1,1) )
  !copy the wavefunctions on GPU
  do iorb=1,orbs%norbp*orbs%nspinor
     iorb_r = (iorb-1)/orbs%nspinor + 1
     ispinor=iorb-orbs%nspinor*(iorb_r-1)
     !print *,'here',iorb,iorb_r,ispinor
!     call ocl_pin_read_buffer_async(GPU%context,GPU%queue,lr%wfd%nvctr_c*8,psi(1,iorb),&
!          GPU%psicf_host(1+(ispinor-1)*2,iorb_r))
     if (pin) call ocl_map_write_buffer_async(GPU%queue, GPU%psicf_host(1,1), &
          (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*(iorb-1), lr%wfd%nvctr_c*8)
     call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_c,lr%wfd%nvctr_c*8,&
          psi(1,iorb))
     if (pin) call ocl_unmap_mem_object(GPU%queue, GPU%psicf_host(1,1), psi(1,iorb))
!     call ocl_release_mem_object(GPU%psicf_host(1+(ispinor-1)*2,iorb_r))

!     call ocl_pin_read_buffer_async(GPU%context,GPU%queue,7*lr%wfd%nvctr_f*8,psi(isf,iorb),&
!          GPU%psicf_host(2+(ispinor-1)*2,iorb_r))
     if (pin) call ocl_map_write_buffer_async(GPU%queue, GPU%psicf_host(1,1), &
          (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8*(iorb-1)+lr%wfd%nvctr_c*8, 7*lr%wfd%nvctr_f*8)
     call ocl_enqueue_write_buffer(GPU%queue,GPU%psi_f,7*lr%wfd%nvctr_f*8,&
          psi(isf,iorb))
     if (pin) call ocl_unmap_mem_object(GPU%queue, GPU%psicf_host(1,1), psi(isf,iorb))
!     call ocl_release_mem_object(GPU%psicf_host(2+(ispinor-1)*2,iorb_r))

     hfac=orbs%kwgts(orbs%iokpt(iorb_r))*orbs%occup(orbs%isorb+iorb_r)/(hxh*hyh*hzh)
     if (orbs%spinsgn(orbs%isorb+iorb_r) > 0.0) then
        rhopot = GPU%rhopot_up
     else
        rhopot = GPU%rhopot_down
     endif
     !calculate the density
     call ocl_locden_generic(GPU%queue, (/lr%d%n1+1,lr%d%n2+1,lr%d%n3+1/),&
          periodic,&
          hfac,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,GPU%keyg_c,GPU%keyv_c,&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,GPU%keyg_f,GPU%keyv_f,&
          GPU%psi_c,GPU%psi_f,&
          GPU%work1,GPU%work2,GPU%work3,&
          rhopot)

  end do

  !copy back the results and leave the uncompressed wavefunctions on the card
  
  if (pin) call ocl_pin_write_buffer_async(GPU%context,GPU%queue,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,rho_p,GPU%rhopot_up_host)
  call ocl_enqueue_read_buffer(GPU%queue,GPU%rhopot_up,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,rho_p)
  if (pin) call ocl_release_mem_object(GPU%rhopot_up_host)
  if( nspin == 2 ) then
     if (pin) call ocl_pin_write_buffer_async(GPU%context,GPU%queue,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,&
          rho_p(1,1,1,2),GPU%rhopot_down_host)
    call ocl_enqueue_read_buffer(GPU%queue,GPU%rhopot_down,lr%d%n1i*lr%d%n2i*lr%d%n3i*8,rho_p(1,1,1,2))
    if (pin) call ocl_release_mem_object(GPU%rhopot_down_host)
  endif

  !free pinning information for potential
  if (pin .and. orbs%norbp>0) call ocl_release_mem_object(GPU%psicf_host(1,1))

END SUBROUTINE local_partial_density_OCL
