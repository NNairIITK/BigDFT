!> @file
!!  Daubechies to Interpolation scaling functions routines
!! @author
!!    Copyright (C) 2010-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Initialize work arrays for local hamiltonian
subroutine initialize_work_arrays_locham(nlr,lr,nspinor,allocate_arrays,w)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nlr, nspinor
  type(locreg_descriptors), dimension(nlr), intent(in) :: lr
  logical,intent(in) :: allocate_arrays
  type(workarr_locham), intent(out) :: w
  !local variables
  character(len=*), parameter :: subname='initialize_work_arrays_locham'
  integer :: ilr
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i,nw,nww,nf
  character(len=1) :: geo
  logical :: hyb

  ! Determine the maximum array sizes for all locregs 1,..,nlr
  ! If the sizes for a specific locreg are needed, simply call the routine with nlr=1
  ! For the moment the geocode of all locregs must be the same
  n1=0
  n2=0
  n3=0
  n1i=0
  n2i=0
  n3i=0
  nfl1=1000000000
  nfl2=1000000000
  nfl3=1000000000
  nfu1=0
  nfu2=0
  nfu3=0
  geo=lr(1)%geocode
  hyb=lr(1)%hybrid_on
  do ilr=1,nlr
      n1=max(n1,lr(ilr)%d%n1)
      n2=max(n2,lr(ilr)%d%n2)
      n3=max(n3,lr(ilr)%d%n3)
      n1i=max(n1i,lr(ilr)%d%n1i)
      n2i=max(n2i,lr(ilr)%d%n2i)
      n3i=max(n3i,lr(ilr)%d%n3i)
      nfl1=min(nfl1,lr(ilr)%d%nfl1)
      nfl2=min(nfl2,lr(ilr)%d%nfl2)
      nfl3=min(nfl3,lr(ilr)%d%nfl3)
      nfu1=max(nfu1,lr(ilr)%d%nfu1)
      nfu2=max(nfu2,lr(ilr)%d%nfu2)
      nfu3=max(nfu3,lr(ilr)%d%nfu3)
      if (lr(ilr)%geocode /= geo) stop 'lr(ilr)%geocode/=geo'
      if (lr(ilr)%hybrid_on .neqv. hyb) stop 'lr(ilr)%hybrid_on .neqv. hyb'
  end do


  if (allocate_arrays) then !this might create memory leaks if there is no check performed
     !if (associated(w%xc)) &
     !     call f_err_throw('Error in initialize_work_arrays_locham: arrays already allocated',&
     !     err_name='BIGDFT_RUNTIME_ERROR')
     nullify(w%w1)
     nullify(w%w2)
     nullify(w%x_c)
     nullify(w%y_c)
     nullify(w%x_f1)
     nullify(w%x_f2)
     nullify(w%x_f3)
     nullify(w%x_f)
     nullify(w%y_f)
  end if


  select case(geo)
  case('F')
     !dimensions of work arrays
     ! shrink convention: nw1>nw2
     w%nw1=max((n3+1)*(2*n1+31)*(2*n2+31),&
           (n1+1)*(2*n2+31)*(2*n3+31),&
           2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
           2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     w%nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
           4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
           (n1+1)*(n2+1)*(2*n3+31),&
           (2*n1+31)*(n2+1)*(n3+1))
     w%nyc=(n1+1)*(n2+1)*(n3+1)
     w%nyf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxc=(n1+1)*(n2+1)*(n3+1)
     w%nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxf1=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxf2=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxf3=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

     !allocation of work arrays
     if (allocate_arrays) then
         w%y_c = f_malloc_ptr((/ w%nyc, nspinor /),id='w%y_c')
         w%y_f = f_malloc_ptr((/ w%nyf, nspinor /),id='w%y_f')
         w%x_c = f_malloc_ptr((/ w%nxc, nspinor /),id='w%x_c')
         w%x_f = f_malloc_ptr((/ w%nxf, nspinor /),id='w%x_f')
         w%w1 = f_malloc_ptr(w%nw1,id='w%w1')
         w%w2 = f_malloc_ptr(w%nw2,id='w%w2')
         w%x_f1 = f_malloc_ptr((/ w%nxf1, nspinor /),id='w%x_f1')
         w%x_f2 = f_malloc_ptr((/ w%nxf2, nspinor /),id='w%x_f2')
         w%x_f3 = f_malloc_ptr((/ w%nxf3, nspinor /),id='w%x_f3')
     end if

     !initialisation of the work arrays
     call f_zero(w%x_f1)
     call f_zero(w%x_f2)
     call f_zero(w%x_f3)
     call f_zero(w%x_c)
     call f_zero(w%x_f)
     call f_zero(w%y_c)
     call f_zero(w%y_f)

  case('S')
     w%nw1=0
     w%nw2=0
     w%nyc=n1i*n2i*n3i
     w%nyf=0
     w%nxc=n1i*n2i*n3i
     w%nxf=0
     w%nxf1=0
     w%nxf2=0
     w%nxf3=0

     !allocation of work arrays
     if (allocate_arrays) then
         w%x_c = f_malloc_ptr((/ w%nxc, nspinor /),id='w%x_c')
         w%y_c = f_malloc_ptr((/ w%nyc, nspinor /),id='w%y_c')
     end if

  case('P')
     if (hyb) then
        ! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        nf=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

        nw=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*n1+2),(2*n1+2)*(n2+2)*(n3+2))
        nw=max(nw,2*(n3+1)*(n1+1)*(n2+1))      ! for the comb_shrink_hyb_c
        nw=max(nw,4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1)) ! for the _f

        nww=max(2*(nfu3-nfl3+1)*(2*n1+2)*(2*n2+2),(n3+1)*(2*n1+2)*(2*n2+2))
        nww=max(nww,4*(n2+1)*(n3+1)*(n1+1))   ! for the comb_shrink_hyb_c   
        nww=max(nww,2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1)) ! for the _f

        w%nw1=nw
        w%nw2=nww
        w%nxc=(n1+1)*(n2+1)*(n3+1)
        w%nyc=(n1+1)*(n2+1)*(n3+1)
        w%nxf=7*nf
        w%nyf=7*nf
        w%nxf1=nf
        w%nxf2=nf
        w%nxf3=nf

        w%y_c = f_malloc_ptr((/ w%nyc, nspinor /),id='w%y_c')
        w%y_f = f_malloc_ptr((/ w%nyf, nspinor /),id='w%y_f')
        w%x_c = f_malloc_ptr((/ w%nxc, nspinor /),id='w%x_c')
        w%x_f = f_malloc_ptr((/ w%nxf, nspinor /),id='w%x_f')
        w%w1 = f_malloc_ptr(w%nw1,id='w%w1')
        w%w2 = f_malloc_ptr(w%nw2,id='w%w2')
        w%x_f1 = f_malloc_ptr((/ w%nxf1, nspinor /),id='w%x_f1')
        w%x_f2 = f_malloc_ptr((/ w%nxf2, nspinor /),id='w%x_f2')
        w%x_f3 = f_malloc_ptr((/ w%nxf3, nspinor /),id='w%x_f3')

     else

        w%nw1=0
        w%nw2=0
        w%nyc=n1i*n2i*n3i
        w%nyf=0
        w%nxc=n1i*n2i*n3i
        w%nxf=0
        w%nxf1=0
        w%nxf2=0
        w%nxf3=0

        if (allocate_arrays) then
            w%x_c = f_malloc_ptr((/ w%nxc, nspinor /),id='w%x_c')
            w%y_c = f_malloc_ptr((/ w%nyc, nspinor /),id='w%y_c')
        end if
     endif
  end select

END SUBROUTINE initialize_work_arrays_locham


subroutine memspace_work_arrays_locham(lr,memwork) !n(c) nspinor (arg:2)
  !n(c) use module_base
  use module_types
  implicit none
  !n(c) integer, intent(in) :: nspinor
  type(locreg_descriptors), intent(in) :: lr
  integer(kind=8), intent(out) :: memwork
  !local variables
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i,nw,nww,nf
  integer :: nw1,nw2,nxc,nxf,nyc,nyf,nxf1,nxf2,nxf3

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3

  select case(lr%geocode) 
  case('F')
     !dimensions of work arrays
     ! shrink convention: nw1>nw2
     nw1=max((n3+1)*(2*n1+31)*(2*n2+31),&
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

     nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))

     nyc=(n1+1)*(n2+1)*(n3+1)
     nyf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     nxc=(n1+1)*(n2+1)*(n3+1)
     nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     nxf1=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     nxf2=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     nxf3=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

  case('S')
     nw1=0
     nw2=0
     nyc=n1i*n2i*n3i
     nyf=0
     nxc=n1i*n2i*n3i
     nxf=0
     nxf1=0
     nxf2=0
     nxf3=0

  case('P')
     if (lr%hybrid_on) then
        ! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        nf=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

        nw=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*n1+2),(2*n1+2)*(n2+2)*(n3+2))
        nw=max(nw,2*(n3+1)*(n1+1)*(n2+1))      ! for the comb_shrink_hyb_c
        nw=max(nw,4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1)) ! for the _f

        nww=max(2*(nfu3-nfl3+1)*(2*n1+2)*(2*n2+2),(n3+1)*(2*n1+2)*(2*n2+2))
        nww=max(nww,4*(n2+1)*(n3+1)*(n1+1))   ! for the comb_shrink_hyb_c   
        nww=max(nww,2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1)) ! for the _f

        nw1=nw
        nw2=nww
        nxc=(n1+1)*(n2+1)*(n3+1)
        nyc=(n1+1)*(n2+1)*(n3+1)
        nxf=7*nf
        nyf=7*nf
        nxf1=nf
        nxf2=nf
        nxf3=nf

     else

        nw1=0
        nw2=0
        nyc=n1i*n2i*n3i
        nyf=0
        nxc=n1i*n2i*n3i
        nxf=0
        nxf1=0
        nxf2=0
        nxf3=0

     endif
  end select

  memwork=nw1+nw2+nxc+nxf+nyc+nyf+nxf1+nxf2+nxf3

END SUBROUTINE memspace_work_arrays_locham


!> Set to zero the work arrays for local hamiltonian
subroutine zero_work_arrays_locham(lr,nspinor,w)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspinor
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(out) :: w
  !local variables
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3

  select case(lr%geocode)

  case('F')

     w%nyc=(n1+1)*(n2+1)*(n3+1)
     w%nyf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxc=(n1+1)*(n2+1)*(n3+1)
     w%nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxf1=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxf2=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
     w%nxf3=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

     !initialisation of the work arrays
     call f_zero(w%x_f1)
     call f_zero(w%x_f2)
     call f_zero(w%x_f3)
     call f_zero(w%x_c)
     call f_zero(w%x_f)
     call f_zero(w%y_c)
     call f_zero(w%y_f)

  case('S')

  case('P')

  end select

END SUBROUTINE zero_work_arrays_locham


subroutine deallocate_work_arrays_locham(w)
  use module_base
  use module_types
  implicit none
  type(workarr_locham), intent(inout) :: w
  !local variables
  character(len=*), parameter :: subname='deallocate_work_arrays_locham'
  
  call f_free_ptr(w%y_c)
  call f_free_ptr(w%x_c)
  call f_free_ptr(w%x_f1)
  call f_free_ptr(w%x_f2)
  call f_free_ptr(w%x_f3)
  call f_free_ptr(w%y_f)
  call f_free_ptr(w%x_f)
  call f_free_ptr(w%w1)
  call f_free_ptr(w%w2)
END SUBROUTINE deallocate_work_arrays_locham


subroutine psi_to_tpsi(hgrids,kptv,nspinor,lr,psi,w,hpsi,ekin,k_strten)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspinor
  real(gp), dimension(3), intent(in) :: hgrids,kptv
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(inout) :: w
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(in) :: psi
  real(gp), intent(out) :: ekin
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(inout) :: hpsi
  real(wp), dimension(6), optional :: k_strten
  !Local variables
  logical, parameter :: transpose=.false.
  logical :: usekpts
  integer :: idx,i,i_f,iseg_f,ipsif,isegf
  real(gp) :: ekino
  real(wp), dimension(0:3) :: scal
  real(gp), dimension(3) :: hgridh
  real(wp), dimension(6) :: kstrten,kstrteno


  !control whether the k points are to be used
  !real k-point different from Gamma still not implemented
  usekpts = nrm2(3,kptv(1),1) > 0.0_gp .or. nspinor == 2

  hgridh=.5_gp*hgrids

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+i_f
  isegf=lr%wfd%nseg_c+iseg_f

    !call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
  ekin=0.0_gp
  
  kstrten=0.0_wp
  select case(lr%geocode)
  case('F')

     !here kpoints cannot be used (for the moment, to be activated for the 
     !localisation region scheme
     if (usekpts) stop 'K points not allowed for Free BC locham'

     do idx=1,nspinor
        call uncompress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  & 
             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),  & 
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
             scal,psi(1,idx),psi(ipsif,idx),  &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx))

        call f_zero(w%nyc,w%y_c(1,idx))
        call f_zero(w%nyf,w%y_f(1,idx))

        call ConvolkineticT(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
             hgrids(1),&        !here the grid spacings are supposed to be equal
             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f, &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%y_c(1,idx),w%y_f(1,idx),ekino, &
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx))
        ekin=ekin+ekino

        !new compression routine in standard form
        call compress_and_accumulate_standard(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             w%y_c(1,idx),w%y_f(1,idx),&
             hpsi(1,idx),hpsi(ipsif,idx))

     end do

  case('S')

     if (usekpts) then
        !first calculate the proper arrays then transpose them before passing to the
        !proper routine
        do idx=1,nspinor
           call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))
        end do

        !Transposition of the work arrays (use y_c as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%x_c,w%y_c,.true.)
        call f_zero(nspinor*w%nyc,w%y_c(1,1))

        ! compute the kinetic part and add  it to psi_out
        ! the kinetic energy is calculated at the same time
        ! do this thing for both components of the spinors
        do idx=1,nspinor,2
           call convolut_kinetic_slab_T_k(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino,kptv(1),kptv(2),kptv(3))
        ekin=ekin+ekino        
        end do

        !re-Transposition of the work arrays (use x_c as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%y_c,w%x_c,.false.)

        do idx=1,nspinor
           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),w%x_c(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

        end do
        
     else
        do idx=1,nspinor
           call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))

           call f_zero(w%nyc,w%y_c(1,idx))
           ! compute the kinetic part and add  it to psi_out
           ! the kinetic energy is calculated at the same time
           call convolut_kinetic_slab_T(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino)
           ekin=ekin+ekino

           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),w%x_c(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))
        end do
     end if

  case('P')
     
     if (lr%hybrid_on) then

        !here kpoints cannot be used, such BC are used in general to mimic the Free BC
        if (usekpts) stop 'K points not allowed for hybrid BC locham'

        !here the grid spacing is not halved
        hgridh=hgrids
        do idx=1,nspinor
           call uncompress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%x_f(1,idx),&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3)

           call f_zero(w%nyc,w%y_c(1,idx))
           call f_zero(w%nyf,w%y_f(1,idx))

           call convolut_kinetic_hyb_T(lr%d%n1,lr%d%n2,lr%d%n3, &
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
                hgridh,w%x_c(1,idx),w%x_f(1,idx),w%y_c(1,idx),w%y_f(1,idx),kstrteno,&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),lr%bounds%kb%ibyz_f,&
                lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
            kstrten=kstrten+kstrteno
            !ekin=ekin+ekino

            call compress_and_accumulate_standard(lr%d,lr%wfd,&
                 lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                 lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                 w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

        end do
     else

        if (usekpts) then

           do idx=1,nspinor
              call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                   lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                   lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                   psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))
           end do

           if (transpose) then
              !Transposition of the work arrays (use psir as workspace)
              call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                   w%x_c,w%y_c,.true.)

              call f_zero(w%y_c)
              ! compute the kinetic part and add  it to psi_out
              ! the kinetic energy is calculated at the same time
              do idx=1,nspinor,2
                 !print *,'AAA',2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,hgridh

                 call convolut_kinetic_per_T_k(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                      hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kptv(1),kptv(2),kptv(3))
                 kstrten=kstrten+kstrteno
                 !ekin=ekin+ekino
              end do

              !Transposition of the work arrays (use psir as workspace)
              call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                   w%y_c,w%x_c,.false.)

           else
              call f_zero(w%y_c)
              do idx=1,nspinor,2
                 call convolut_kinetic_per_T_k_notranspose(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                      hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kptv(1),kptv(2),kptv(3))
                 kstrten=kstrten+kstrteno
              end do
           end if

           do idx=1,nspinor

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),w%x_c(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                 lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                 lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                 w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

           end do
        else
           !first calculate the proper arrays then transpose them before passing to the
           !proper routine
           do idx=1,nspinor
              call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                   lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                   lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                   psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))

              call f_zero(w%nyc,w%y_c(1,idx))
              ! compute the kinetic part and add  it to psi_out
              ! the kinetic energy is calculated at the same time
              call convolut_kinetic_per_t(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno)
              kstrten=kstrten+kstrteno

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),w%x_c(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                   lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                   lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                   w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

           end do
        end if

     end if
     ekin=ekin+kstrten(1)+kstrten(2)+kstrten(3)
     if (present(k_strten)) k_strten=kstrten 

  end select

END SUBROUTINE psi_to_tpsi


subroutine daub_to_isf_locham(nspinor,lr,w,psi,psir)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspinor
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(out) :: psir
  !local variables
  integer :: idx,i,i_f,iseg_f
  real(wp), dimension(0:3) :: scal

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)

  !call f_zero((2*n1+31)*(2*n2+31)*(2*n3+31)*nspinor,psir)
  !call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
  select case(lr%geocode)
  case('F')
     call f_zero(psir)
     !call timing(iproc,'CrtDescriptors','ON') !temporary
     do idx=1,nspinor  
        call uncompress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  & 
             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),  & 
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
             scal,psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),  &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx))
        !call timing(iproc,'CrtDescriptors','OF') !temporary
        !call timing(iproc,'CrtLocPot     ','ON') !temporary
        call comb_grow_all(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%w1,w%w2,w%x_c(1,idx),w%x_f(1,idx), & 
             psir(1,idx),lr%bounds%kb%ibyz_c,lr%bounds%gb%ibzxx_c,&
             lr%bounds%gb%ibxxyy_c,lr%bounds%gb%ibyz_ff,&
             lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)
        !call timing(iproc,'CrtLocPot     ','OF') !temporary
     end do
     
  case('S')

     do idx=1,nspinor
        call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
             psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),w%x_c(1,idx),psir(1,idx))

        call convolut_magic_n_slab(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,w%x_c(1,idx),psir(1,idx),w%y_c(1,idx)) 

     end do

  case('P')
     
     if (lr%hybrid_on) then

        do idx=1,nspinor
           call uncompress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
                psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),w%x_c(1,idx),w%x_f(1,idx),&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3)
           
           call comb_grow_all_hybrid(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
                w%nw1,w%nw2,&
                w%w1,w%w2,w%x_c(1,idx),w%x_f(1,idx),psir(1,idx),lr%bounds%gb)
        end do
        
     else

        do idx=1,nspinor
           call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
                psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),w%x_c(1,idx),psir(1,idx))

           call convolut_magic_n_per(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,w%x_c(1,idx),psir(1,idx),w%y_c(1,idx)) 
        end do

     end if
  end select
  
END SUBROUTINE daub_to_isf_locham


subroutine isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,w,psir,hpsi,ekin,k_strten)
  use module_defs, only: gp,wp
  use module_types
  implicit none
  integer, intent(in) :: nspinor
  real(gp), intent(in) :: hx,hy,hz,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(inout) :: w
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(in) :: psir
  real(gp), intent(out) :: ekin
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(inout) :: hpsi
  real(wp), dimension(6), optional :: k_strten
  !Local variables
  logical :: usekpts
  integer :: idx,i,i_f,iseg_f,ipsif,isegf
  real(gp) :: ekino
  real(wp), dimension(0:3) :: scal
  real(gp), dimension(3) :: hgridh
  real(wp), dimension(6) :: kstrten,kstrteno


  !control whether the k points are to be used
  !real k-point different from Gamma still not implemented
  usekpts = kx**2+ky**2+kz**2 > 0.0_gp .or. nspinor == 2

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+i_f
  isegf=lr%wfd%nseg_c+iseg_f

    !call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
  ekin=0.0_gp
  
  kstrten=0.0_wp
  select case(lr%geocode)
  case('F')

     !here kpoints cannot be used (for the moment, to be activated for the 
     !localisation region scheme
     if (usekpts) stop 'K points not allowed for Free BC locham'

     do idx=1,nspinor

        call comb_shrink(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%w1,w%w2,psir(1,idx),&
             lr%bounds%kb%ibxy_c,lr%bounds%sb%ibzzx_c,lr%bounds%sb%ibyyzz_c,&
             lr%bounds%sb%ibxy_ff,lr%bounds%sb%ibzzx_f,lr%bounds%sb%ibyyzz_f,&
             w%y_c(1,idx),w%y_f(1,idx))

        call ConvolkineticT(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
             hx,&        !here the grid spacings are supposed to be equal
             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f, &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%y_c(1,idx),w%y_f(1,idx),ekino, &
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx))
        ekin=ekin+ekino

        !new compression routine in standard form
        call compress_and_accumulate_standard(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             w%y_c(1,idx),w%y_f(1,idx),&
             hpsi(1,idx),hpsi(ipsif,idx))
!!$        call compress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
!!$             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$             lr%wfd%keygloc(1,1),lr%wfd%keyv(1),&
!!$             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),   &
!!$             scal,w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx))

     end do

  case('S')

     if (usekpts) then
        !first calculate the proper arrays then transpose them before passing to the
        !proper routine
        do idx=1,nspinor
           call convolut_magic_t_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                psir(1,idx),w%y_c(1,idx))
        end do

        !Transposition of the work arrays (use psir as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%x_c,psir,.true.)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%y_c,psir,.true.)

        ! compute the kinetic part and add  it to psi_out
        ! the kinetic energy is calculated at the same time
        ! do this thing for both components of the spinors
        do idx=1,nspinor,2
           call convolut_kinetic_slab_T_k(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino,kx,ky,kz)
        ekin=ekin+ekino        
        end do

        !re-Transposition of the work arrays (use psir as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%y_c,psir,.false.)

        do idx=1,nspinor
           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),psir(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$           call compress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                lr%wfd%keygloc(1,1),lr%wfd%keyv(1),   & 
!!$                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),   & 
!!$                w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
        end do
        
     else
        do idx=1,nspinor
           call convolut_magic_t_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                psir(1,idx),w%y_c(1,idx))

           ! compute the kinetic part and add  it to psi_out
           ! the kinetic energy is calculated at the same time
           call convolut_kinetic_slab_T(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino)
           ekin=ekin+ekino

           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),psir(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$           call compress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                lr%wfd%keygloc(1,1),lr%wfd%keyv(1),   & 
!!$                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),   & 
!!$                w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
        end do
     end if

  case('P')
     
     if (lr%hybrid_on) then

        !here kpoints cannot be used, such BC are used in general to mimic the Free BC
        if (usekpts) stop 'K points not allowed for hybrid BC locham'

        !here the grid spacing is not halved
        hgridh(1)=hx
        hgridh(2)=hy
        hgridh(3)=hz
        do idx=1,nspinor
           call comb_shrink_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
                w%w2,w%w1,psir(1,idx),w%y_c(1,idx),w%y_f(1,idx),lr%bounds%sb)

           call convolut_kinetic_hyb_T(lr%d%n1,lr%d%n2,lr%d%n3, &
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
                hgridh,w%x_c(1,idx),w%x_f(1,idx),w%y_c(1,idx),w%y_f(1,idx),kstrteno,&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),lr%bounds%kb%ibyz_f,&
                lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
            kstrten=kstrten+kstrteno
            !ekin=ekin+ekino

            call compress_and_accumulate_standard(lr%d,lr%wfd,&
                 lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                 lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                 w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$           call compress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                lr%wfd%keygloc(1,1),lr%wfd%keyv(1),& 
!!$                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f), & 
!!$                w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),&
!!$                lr%d%nfl1,lr%d%nfl2,lr%d%nfl3,lr%d%nfu1,lr%d%nfu2,lr%d%nfu3)
        end do
     else

        if (usekpts) then
           !first calculate the proper arrays then transpose them before passing to the
           !proper routine
           do idx=1,nspinor
              call convolut_magic_t_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   psir(1,idx),w%y_c(1,idx))
           end do

           !Transposition of the work arrays (use psir as workspace)
           call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                w%x_c,psir,.true.)
           call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                w%y_c,psir,.true.)


           ! compute the kinetic part and add  it to psi_out
           ! the kinetic energy is calculated at the same time
           do idx=1,nspinor,2
              !print *,'AAA',2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,hgridh
              
              call convolut_kinetic_per_T_k(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kx,ky,kz)
              kstrten=kstrten+kstrteno
              !ekin=ekin+ekino
           end do

           !Transposition of the work arrays (use psir as workspace)
           call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                w%y_c,psir,.false.)

           do idx=1,nspinor

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),psir(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                 lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                 lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                 psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$              call compress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                   lr%wfd%keygloc(1,1),lr%wfd%keyv(1),& 
!!$                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                   lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),&
!!$                   w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
           end do
        else
           !first calculate the proper arrays then transpose them before passing to the
           !proper routine
           do idx=1,nspinor
              call convolut_magic_t_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   psir(1,idx),w%y_c(1,idx))
              ! compute the kinetic part and add  it to psi_out
              ! the kinetic energy is calculated at the same time
              call convolut_kinetic_per_t(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno)
              kstrten=kstrten+kstrteno

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),psir(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                   lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                   lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                   psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$              call compress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                   lr%wfd%keygloc(1,1),lr%wfd%keyv(1),& 
!!$                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                   lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),& 
!!$                   w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
           end do
        end if

     end if
     ekin=ekin+kstrten(1)+kstrten(2)+kstrten(3)
     if (present(k_strten)) k_strten=kstrten 

  end select

END SUBROUTINE isf_to_daub_kinetic


subroutine initialize_work_arrays_sumrho(nlr,lr,allocate_arrays,w)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nlr
  type(locreg_descriptors), dimension(nlr), intent(in) :: lr
  logical, intent(in) :: allocate_arrays
  type(workarr_sumrho), intent(out) :: w
  !local variables
  character(len=*), parameter :: subname='initialize_work_arrays_sumrho'
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3!n(c) n1i,n2i,n3i
  integer :: ilr
  character(len=1) :: geo
  logical :: hyb

  ! Determine the maximum array sizes for all locregs 1,..,nlr
  ! If the sizes for a specific locreg are needed, simply call the routine with nlr=1
  ! For the moment the geocode of all locregs must be the same
  
  n1=0
  n2=0
  n3=0
  nfl1=1000000000
  nfl2=1000000000
  nfl3=1000000000
  nfu1=0
  nfu2=0
  nfu3=0
  geo=lr(1)%geocode
  hyb=lr(1)%hybrid_on
  do ilr=1,nlr
      n1=max(n1,lr(ilr)%d%n1)
      n2=max(n2,lr(ilr)%d%n2)
      n3=max(n3,lr(ilr)%d%n3)
      nfl1=min(nfl1,lr(ilr)%d%nfl1)
      nfl2=min(nfl2,lr(ilr)%d%nfl2)
      nfl3=min(nfl3,lr(ilr)%d%nfl3)
      nfu1=max(nfu1,lr(ilr)%d%nfu1)
      nfu2=max(nfu2,lr(ilr)%d%nfu2)
      nfu3=max(nfu3,lr(ilr)%d%nfu3)
      if (lr(ilr)%geocode /= geo) then
          write(*,*) 'lr(ilr)%geocode, geo', lr(ilr)%geocode, geo
          stop 'lr(ilr)%geocode/=geo'
      end if
      if (lr(ilr)%hybrid_on .neqv. hyb) stop 'lr(ilr)%hybrid_on .neqv. hyb'
  end do

  if (allocate_arrays) then
     nullify(w%x_c)
     nullify(w%x_f)
     nullify(w%w1)
     nullify(w%w2)
  end if

  select case(geo)
  case('F')
     !dimension of the work arrays
     ! shrink convention: nw1>nw2
     w%nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     w%nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))
     w%nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
     w%nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
  case('S')
     !dimension of the work arrays
     w%nw1=1
     w%nw2=1
     w%nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
     w%nxf=1
  case('P')
     if (hyb) then
        ! hybrid case:
        w%nxc=(n1+1)*(n2+1)*(n3+1)
        w%nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

        w%nw1=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*n1+2),(2*n1+2)*(n2+2)*(n3+2))
        w%nw1=max(w%nw1,2*(n3+1)*(n1+1)*(n2+1))      ! for the comb_shrink_hyb_c
        w%nw1=max(w%nw1,4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1)) ! for the _f

        w%nw2=max(2*(nfu3-nfl3+1)*(2*n1+2)*(2*n2+2),(n3+1)*(2*n1+2)*(2*n2+2))
        w%nw2=max(w%nw2,4*(n2+1)*(n3+1)*(n1+1))   ! for the comb_shrink_hyb_c   
        w%nw2=max(w%nw2,2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1)) ! for the _f
     else
        !dimension of the work arrays, fully periodic case
        w%nw1=1
        w%nw2=1
        w%nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
        w%nxf=1
     endif

  end select
  !work arrays
  if (allocate_arrays) then
      w%x_c = f_malloc_ptr(w%nxc,id='w%x_c')
      w%x_f = f_malloc_ptr(w%nxf,id='w%x_f')
      w%w1 = f_malloc_ptr(w%nw1,id='w%w1')
      w%w2 = f_malloc_ptr(w%nw2,id='w%w2')
  end if
  

  if (geo == 'F') then
     call f_zero(w%x_c)
     call f_zero(w%x_f)
  end if


END SUBROUTINE initialize_work_arrays_sumrho


subroutine memspace_work_arrays_sumrho(lr,memwork)
  !n(c) use module_base
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  integer(kind=8), intent(out) :: memwork
  !local variables
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer :: nw1,nw2,nxc,nxf

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3

  select case(lr%geocode)
  case('F')
     !dimension of the work arrays
     ! shrink convention: nw1>nw2
     nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))
     nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
  case('S')
     !dimension of the work arrays
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
     nxf=1
  case('P')
     if (lr%hybrid_on) then
        ! hybrid case:
        nxc=(n1+1)*(n2+1)*(n3+1)
        nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

        nw1=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*n1+2),(2*n1+2)*(n2+2)*(n3+2))
        nw1=max(nw1,2*(n3+1)*(n1+1)*(n2+1))      ! for the comb_shrink_hyb_c
        nw1=max(nw1,4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1)) ! for the _f

        nw2=max(2*(nfu3-nfl3+1)*(2*n1+2)*(2*n2+2),(n3+1)*(2*n1+2)*(2*n2+2))
        nw2=max(nw2,4*(n2+1)*(n3+1)*(n1+1))   ! for the comb_shrink_hyb_c   
        nw2=max(nw2,2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1)) ! for the _f
     else
        !dimension of the work arrays, fully periodic case
        nw1=1
        nw2=1
        nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
        nxf=1
     endif

  end select
  memwork=nxc+nxf+nw1+nw2

END SUBROUTINE memspace_work_arrays_sumrho


subroutine deallocate_work_arrays_sumrho(w)
  use module_base
  use module_types
  implicit none
  type(workarr_sumrho), intent(inout) :: w
  !local variables
  character(len=*), parameter :: subname='deallocate_work_arrays_sumrho'

  call f_free_ptr(w%x_c)
  call f_free_ptr(w%x_f)
  call f_free_ptr(w%w1)
  call f_free_ptr(w%w2)
  
END SUBROUTINE deallocate_work_arrays_sumrho


subroutine psig_to_psir_free(n1,n2,n3,work,psig_psir)
 implicit none
 integer, intent(in) :: n1,n2,n3 !< dimensions in the daubechies grid
 real(kind=8),dimension((2*n1+31)*(2*n2+31)*(2*n3+31)), intent(inout):: work !< enlarged buffer 
 real(kind=8),dimension((2*n1+31)*(2*n2+31)*(2*n3+31)), intent(inout):: psig_psir  !< final result, containing psig data but big enough to contain psir (used as work array)

 call synthese_free_self(n1,n2,n3,psig_psir,work)
 call convolut_magic_n_free_self(2*n1+15,2*n2+15,2*n3+15,work,psig_psir)

end subroutine psig_to_psir_free


!> Transform a daubechies function in compressed form to a function in real space via
!! the Magic Filter operation
!! do this for a single component (spinorial and/or complex)
subroutine daub_to_isf(lr,w,psi,psir)
  use module_base
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_sumrho), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psir
  !local variables
  integer :: i,i_f,iseg_f
  real(wp), dimension(0:3) :: scal

  i_f=min(lr%wfd%nvctr_f,1)
  iseg_f=min(lr%wfd%nseg_f,1)

  do i=0,3
     scal(i)=1.0_wp
  enddo

  select case(lr%geocode)
  case('F')
     call f_zero(psir)

     call uncompress_forstandard_short(lr%d%n1,lr%d%n2,lr%d%n3,&
          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),  & 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), &
          scal,psi(1),psi(lr%wfd%nvctr_c+i_f),&
          w%x_c,w%x_f)

     call comb_grow_all(lr%d%n1,lr%d%n2,lr%d%n3,&
          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
          w%w1,w%w2,w%x_c,w%x_f,  & 
          psir,lr%bounds%kb%ibyz_c,&
          lr%bounds%gb%ibzxx_c,lr%bounds%gb%ibxxyy_c,&
          lr%bounds%gb%ibyz_ff,lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)

  case('P')
     if (lr%hybrid_on) then
        ! hybrid case
        call uncompress_per_f_short(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,&
             lr%wfd%nvctr_c,lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),& 
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), &
             psi(1),psi(lr%wfd%nvctr_c+i_f),&
             w%x_c,w%x_f,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3)

        call comb_grow_all_hybrid(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%nw1,w%nw2,w%w1,w%w2,w%x_c,w%x_f,psir,lr%bounds%gb)
     else
        call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,&
             lr%wfd%nvctr_c,lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),&
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),&
             psi,psi(lr%wfd%nvctr_c+i_f),&
             w%x_c,psir)

        call convolut_magic_n_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
             w%x_c,psir) 
     endif

  case('S')
     call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
          lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),&
          lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), &
          psi(1),psi(lr%wfd%nvctr_c+i_f),w%x_c,&
          psir)

     call convolut_magic_n_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,w%x_c,&
          psir) 

  end select

END SUBROUTINE daub_to_isf


!> Transforms a wavefunction written in real space basis into a 
!! wavefunction in Daubechies form
!! does the job for all supported BC
!! @warning: 
!!  - the psir is destroyed for some BCs (slab and periodic)
!!  - psi must already be initialized (to zero) before entering this routine
subroutine isf_to_daub(lr,w,psir,psi)
  use module_defs, only: wp
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_sumrho), intent(inout) :: w
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(inout) :: psir
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(inout) :: psi
  !local variables
  integer :: i,i_f,iseg_f,isegf,ipsif
  real(wp), dimension(0:3) :: scal

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+i_f
  isegf=lr%wfd%nseg_c+iseg_f

  select case(lr%geocode)
  case('F')
     call comb_shrink(lr%d%n1,lr%d%n2,lr%d%n3,&
          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
          w%w1,w%w2,psir,&
          lr%bounds%kb%ibxy_c,lr%bounds%sb%ibzzx_c,lr%bounds%sb%ibyyzz_c,&
          lr%bounds%sb%ibxy_ff,lr%bounds%sb%ibzzx_f,lr%bounds%sb%ibyyzz_f,&
          w%x_c,w%x_f)
     
     call compress_and_accumulate_standard(lr%d,lr%wfd,&
          lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
          lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
          w%x_c,w%x_f,psi(1),psi(ipsif))

!!$     call compress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
!!$          lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$          lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),&
!!$          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$          lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
!!$          scal,w%x_c,w%x_f,psi(1),psi(lr%wfd%nvctr_c+i_f))
  case('S')

     call convolut_magic_t_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
          psir(1),w%x_c(1))
     
     call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
          w%x_c,psir(1))
     call compress_and_accumulate_mixed(lr%d,lr%wfd,&
          lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
          lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
          psir(1),psi(1),psi(ipsif))

!!$     call compress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$          lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$          lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   & 
!!$          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$          lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   & 
!!$          w%x_c(1),psi(1),psi(lr%wfd%nvctr_c+i_f),psir(1))

  case('P')
     
     if (lr%hybrid_on) then

        call comb_shrink_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%w2,w%w1,psir(1),w%x_c(1),w%x_f(1),lr%bounds%sb)
        call compress_and_accumulate_standard(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             w%x_c(1),w%x_f(1),psi(1),psi(ipsif))

!!$        call compress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),& 
!!$             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), & 
!!$             w%x_c(1),w%x_f(1),psi(1),psi(lr%wfd%nvctr_c+i_f),&
!!$             lr%d%nfl1,lr%d%nfl2,lr%d%nfl3,lr%d%nfu1,lr%d%nfu2,lr%d%nfu3)
     else

        call convolut_magic_t_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
             psir(1),w%x_c(1))

        call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
             w%x_c(1),psir(1))
        call compress_and_accumulate_mixed(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             psir(1),psi(1),psi(ipsif))

!!$        call compress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),& 
!!$             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),& 
!!$ w%x_c(1),psi(1),psi(lr%wfd%nvctr_c+i_f),psir(1)) 
     end if
        
  end select

END SUBROUTINE isf_to_daub
