module locreg_operations
  use module_base
  use locregs
  implicit none

  private

  !> Information for the confining potential to be used in TMB scheme
  !! The potential is supposed to be defined as prefac*(r-rC)**potorder
  type, public :: confpot_data
     integer :: potorder                !< Order of the confining potential
     integer, dimension(3) :: ioffset   !< Offset for the coordinates of potential lr in global region
     real(gp) :: prefac                 !< Prefactor
     real(gp), dimension(3) :: hh       !< Grid spacings in ISF grid
     real(gp), dimension(3) :: rxyzConf !< Confining potential center in global coordinates
     real(gp) :: damping                !< Damping factor to be used after the restart
  end type confpot_data

  !> Contains the work arrays needed for expressing wavefunction in real space
  !! with all the BC
  type, public :: workarr_sumrho
     integer :: nw1,nw2,nxc,nxf
     real(wp), dimension(:), pointer :: x_c,x_f,w1,w2
  end type workarr_sumrho

  !> Contains the work arrays needed for hamiltonian application with all the BC
  type, public :: workarr_locham
     integer :: nw1,nw2,nxc,nyc,nxf1,nxf2,nxf3,nxf,nyf
     real(wp), dimension(:), pointer :: w1,w2
     !for the periodic BC case, these arrays substitute 
     !psifscf,psifscfk,psig,ww respectively
     real(wp), dimension(:,:), pointer :: x_c,y_c,x_f1,x_f2,x_f3,x_f,y_f
  end type workarr_locham

  !> Contains the work arrays needed for th preconditioner with all the BC
  !! Take different pointers depending on the boundary conditions
  type, public :: workarr_precond
     integer, dimension(:), pointer :: modul1,modul2,modul3
     real(wp), dimension(:), pointer :: psifscf,ww,x_f1,x_f2,x_f3,kern_k1,kern_k2,kern_k3
     real(wp), dimension(:,:), pointer :: af,bf,cf,ef
     real(wp), dimension(:,:,:), pointer :: xpsig_c,ypsig_c,x_c
     real(wp), dimension(:,:,:,:), pointer :: xpsig_f,ypsig_f,x_f,y_f
     real(wp), dimension(:,:,:,:,:), pointer :: z1,z3 ! work array for FFT
  end type workarr_precond

  type, public :: workarrays_quartic_convolutions
     real(wp), dimension(:,:,:), pointer :: xx_c, xy_c, xz_c
     real(wp), dimension(:,:,:), pointer :: xx_f1
     real(wp), dimension(:,:,:), pointer :: xy_f2
     real(wp), dimension(:,:,:), pointer :: xz_f4
     real(wp), dimension(:,:,:,:), pointer :: xx_f, xy_f, xz_f
     real(wp), dimension(:,:,:), pointer :: y_c
     real(wp), dimension(:,:,:,:), pointer :: y_f
     ! The following arrays are work arrays within the subroutine
     real(wp), dimension(:,:), pointer :: aeff0array, beff0array, ceff0array, eeff0array
     real(wp), dimension(:,:), pointer :: aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array
     real(wp), dimension(:,:), pointer :: aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray
     real(wp), dimension(:,:,:), pointer :: xya_c, xyc_c
     real(wp), dimension(:,:,:), pointer :: xza_c, xzc_c
     real(wp), dimension(:,:,:), pointer :: yza_c, yzb_c, yzc_c, yze_c
     real(wp), dimension(:,:,:,:), pointer :: xya_f, xyb_f, xyc_f, xye_f
     real(wp), dimension(:,:,:,:), pointer :: xza_f, xzb_f, xzc_f, xze_f
     real(wp), dimension(:,:,:,:), pointer :: yza_f, yzb_f, yzc_f, yze_f
  end type workarrays_quartic_convolutions


  public :: psir_to_vpsi,isf_to_daub_kinetic
  public :: nullify_confpot_data 
  public :: Lpsi_to_global2
  public :: global_to_local_parallel
  public :: boundary_weight
  public :: psi_to_locreg2
  public :: global_to_local
  public :: initialize_work_arrays_sumrho,deallocate_work_arrays_sumrho
  public :: initialize_work_arrays_locham,deallocate_work_arrays_locham
  public :: memspace_work_arrays_sumrho,memspace_work_arrays_locham
  public :: allocate_work_arrays,init_local_work_arrays,deallocate_work_arrays
  public :: deallocate_workarrays_quartic_convolutions,zero_local_work_arrays

  contains

    !> Initialize work arrays for local hamiltonian
    subroutine initialize_work_arrays_locham(nlr,lr,nspinor,allocate_arrays,w)
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
      implicit none
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
      implicit none
      integer, intent(in) :: nspinor
      type(locreg_descriptors), intent(in) :: lr
      type(workarr_locham), intent(inout) :: w
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


    subroutine initialize_work_arrays_sumrho(nlr,lr,allocate_arrays,w)
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
      implicit none
      type(workarr_sumrho), intent(inout) :: w
      !local variables
      character(len=*), parameter :: subname='deallocate_work_arrays_sumrho'

      call f_free_ptr(w%x_c)
      call f_free_ptr(w%x_f)
      call f_free_ptr(w%w1)
      call f_free_ptr(w%w2)

    END SUBROUTINE deallocate_work_arrays_sumrho

    subroutine allocate_work_arrays(geocode,hybrid_on,ncplx,d,w)
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      logical, intent(in) :: hybrid_on
      integer, intent(in) :: ncplx
      type(grid_dimensions), intent(in) :: d
      type(workarr_precond), intent(out) :: w
      !local variables
      character(len=*), parameter :: subname='allocate_work_arrays'
      integer, parameter :: lowfil=-14,lupfil=14
      integer :: nd1,nd2,nd3
      integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
      integer :: nf

      if (geocode == 'F') then

         nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)
         !allocate work arrays
         w%xpsig_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%xpsig_c')
         w%xpsig_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%xpsig_f')
         w%ypsig_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%ypsig_c')
         w%ypsig_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%ypsig_f')

         w%x_f1 = f_malloc_ptr(nf,id='w%x_f1')
         w%x_f2 = f_malloc_ptr(nf,id='w%x_f2')
         w%x_f3 = f_malloc_ptr(nf,id='w%x_f3')

      else if (geocode == 'P') then

         if (hybrid_on) then

            call dimensions_fft(d%n1,d%n2,d%n3,&
                 nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

            nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)

            w%kern_k1 = f_malloc_ptr(0.to.d%n1,id='w%kern_k1')
            w%kern_k2 = f_malloc_ptr(0.to.d%n2,id='w%kern_k2')
            w%kern_k3 = f_malloc_ptr(0.to.d%n3,id='w%kern_k3')
            w%z1 = f_malloc_ptr((/ 2, nd1b, nd2, nd3, 2 /),id='w%z1')
            w%z3 = f_malloc_ptr((/ 2, nd1, nd2, nd3f, 2 /),id='w%z3')
            w%x_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%x_c')

            w%x_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%x_f')
            w%x_f1 = f_malloc_ptr(nf,id='w%x_f1')
            w%x_f2 = f_malloc_ptr(nf,id='w%x_f2')
            w%x_f3 = f_malloc_ptr(nf,id='w%x_f3')
            w%y_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%y_f')
            w%ypsig_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%ypsig_c')


         else 

            if (ncplx == 1) then
               !periodic, not k-points
               w%modul1 = f_malloc_ptr(lowfil.to.d%n1+lupfil,id='w%modul1')
               w%modul2 = f_malloc_ptr(lowfil.to.d%n2+lupfil,id='w%modul2')
               w%modul3 = f_malloc_ptr(lowfil.to.d%n3+lupfil,id='w%modul3')
               w%af = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%af')
               w%bf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%bf')
               w%cf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%cf')
               w%ef = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%ef')
            end if

            w%psifscf = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+2)*(2*d%n3+2),id='w%psifscf')
            w%ww = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+2)*(2*d%n3+2),id='w%ww')

         end if

      else if (geocode == 'S') then

         if (ncplx == 1) then
            w%modul1 = f_malloc_ptr(lowfil.to.d%n1+lupfil,id='w%modul1')
            w%modul3 = f_malloc_ptr(lowfil.to.d%n3+lupfil,id='w%modul3')
            w%af = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%af')
            w%bf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%bf')
            w%cf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%cf')
            w%ef = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%ef')
         end if

         w%psifscf = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+16)*(2*d%n3+2),id='w%psifscf')
         w%ww = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+16)*(2*d%n3+2),id='w%ww')

      end if

    END SUBROUTINE allocate_work_arrays


    subroutine memspace_work_arrays_precond(geocode,hybrid_on,ncplx,d,memwork)
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      logical, intent(in) :: hybrid_on
      integer, intent(in) :: ncplx
      type(grid_dimensions), intent(in) :: d
      integer(kind=8), intent(out) :: memwork
      !local variables
      integer, parameter :: lowfil=-14,lupfil=14
      integer :: nd1,nd2,nd3
      integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
      integer :: nf


      if (geocode == 'F') then

         nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)

         memwork=2*(d%n1+1)*(d%n2+1)*(d%n3+1)+2*7*(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)+3*nf


      else if (geocode == 'P') then

         if (hybrid_on) then

            call dimensions_fft(d%n1,d%n2,d%n3,&
                 nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

            nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)

            memwork=(d%n1+1)+(d%n2+1)+(d%n3+1)+2*nd1b*nd2*nd3*2+2*nd1*nd2*nd3f*2+&
                 (d%n1+1)*(d%n2+1)*(d%n3+1)+2*7*(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)+3*nf

         else 

            memwork=0
            if (ncplx == 1) then
               memwork=d%n1+d%n2+d%n3+15*(lupfil-lowfil+1)
            end if
            memwork=memwork+2*ncplx*(2*d%n1+2)*(2*d%n2+2)*(2*d%n3+2)

         end if

      else if (geocode == 'S') then

         memwork=0
         if (ncplx == 1) then
            memwork=d%n1+d%n3+14*(lupfil-lowfil+1)
         end if
         memwork=memwork+2*ncplx*(2*d%n1+2)*(2*d%n2+16)*(2*d%n3+2)
      end if

    END SUBROUTINE memspace_work_arrays_precond


    subroutine deallocate_work_arrays(geocode,hybrid_on,ncplx,w)
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      logical, intent(in) :: hybrid_on
      integer, intent(in) :: ncplx
      type(workarr_precond), intent(inout) :: w
      !local variables
      character(len=*), parameter :: subname='deallocate_work_arrays'

      if (geocode == 'F') then

         call f_free_ptr(w%xpsig_c)
         call f_free_ptr(w%ypsig_c)
         call f_free_ptr(w%xpsig_f)
         call f_free_ptr(w%ypsig_f)
         call f_free_ptr(w%x_f1)
         call f_free_ptr(w%x_f2)
         call f_free_ptr(w%x_f3)

      else if ((geocode == 'P' .and. .not. hybrid_on) .or. geocode == 'S') then

         if (ncplx == 1) then
            call f_free_ptr(w%modul1)
            if (geocode /= 'S') then
               call f_free_ptr(w%modul2)
            end if
            call f_free_ptr(w%modul3)
            call f_free_ptr(w%af)
            call f_free_ptr(w%bf)
            call f_free_ptr(w%cf)
            call f_free_ptr(w%ef)
         end if

         call f_free_ptr(w%psifscf)
         call f_free_ptr(w%ww)

      else if (geocode == 'P' .and. hybrid_on) then

         call f_free_ptr(w%z1)
         call f_free_ptr(w%z3)
         call f_free_ptr(w%kern_k1)
         call f_free_ptr(w%kern_k2)
         call f_free_ptr(w%kern_k3)
         call f_free_ptr(w%x_c)
         call f_free_ptr(w%x_f)
         call f_free_ptr(w%x_f1)
         call f_free_ptr(w%x_f2)
         call f_free_ptr(w%x_f3)
         call f_free_ptr(w%y_f)
         call f_free_ptr(w%ypsig_c)


      end if

    END SUBROUTINE deallocate_work_arrays

    subroutine deallocate_workarrays_quartic_convolutions(work)
      implicit none

      ! Calling arguments
      type(workarrays_quartic_convolutions),intent(inout):: work

      ! Local variables
      integer:: iall, istat


      call f_free_ptr(work%xx_c)

      call f_free_ptr(work%xy_c)

      call f_free_ptr(work%xz_c)

      call f_free_ptr(work%xx_f1)

      call f_free_ptr(work%xx_f)

      call f_free_ptr(work%xy_f2)

      call f_free_ptr(work%xy_f)

      call f_free_ptr(work%xz_f4)

      call f_free_ptr(work%xz_f)

      call f_free_ptr(work%y_c)

      call f_free_ptr(work%y_f)

      call f_free_ptr(work%aeff0array)

      call f_free_ptr(work%beff0array)

      call f_free_ptr(work%ceff0array)

      call f_free_ptr(work%eeff0array)

      call f_free_ptr(work%aeff0_2array)

      call f_free_ptr(work%beff0_2array)

      call f_free_ptr(work%ceff0_2array)

      call f_free_ptr(work%eeff0_2array)

      call f_free_ptr(work%aeff0_2auxarray)

      call f_free_ptr(work%beff0_2auxarray)

      call f_free_ptr(work%ceff0_2auxarray)

      call f_free_ptr(work%eeff0_2auxarray)

      call f_free_ptr(work%xya_c)

      call f_free_ptr(work%xyc_c)

      call f_free_ptr(work%xza_c)

      call f_free_ptr(work%xzc_c)

      call f_free_ptr(work%yza_c)

      call f_free_ptr(work%yzb_c)

      call f_free_ptr(work%yzc_c)

      call f_free_ptr(work%yze_c)

      call f_free_ptr(work%xya_f)

      call f_free_ptr(work%xyb_f)

      call f_free_ptr(work%xyc_f)

      call f_free_ptr(work%xye_f)

      call f_free_ptr(work%xza_f)

      call f_free_ptr(work%xzb_f)

      call f_free_ptr(work%xzc_f)

      call f_free_ptr(work%xze_f)

      call f_free_ptr(work%yza_f)

      call f_free_ptr(work%yzb_f)

      call f_free_ptr(work%yzc_f)

      call f_free_ptr(work%yze_f)

    end subroutine deallocate_workarrays_quartic_convolutions


    subroutine init_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work)
      implicit none

      ! Calling arguments
      integer,intent(in)::n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
      logical,intent(in):: with_confpot
      type(workarrays_quartic_convolutions),intent(inout):: work

      ! Local variables
      integer:: i, istat
      integer,parameter :: lowfil=-14,lupfil=14

      work%xx_c = f_malloc0_ptr((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='work%xx_c')
      work%xy_c = f_malloc0_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xy_c')
      work%xz_c = f_malloc0_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xz_c')

      work%xx_f1 = f_malloc0_ptr((/ nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%xx_f1')
      work%xx_f = f_malloc0_ptr((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%xx_f')


      work%xy_f2 = f_malloc0_ptr((/ nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xy_f2')
      work%xy_f = f_malloc0_ptr((/ 1.to.7, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xy_f')


      work%xz_f4 = f_malloc0_ptr((/ nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xz_f4')
      work%xz_f = f_malloc0_ptr((/ 1.to.7, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xz_f')


      work%y_c = f_malloc0_ptr((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='work%y_c')

      work%y_f = f_malloc0_ptr((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%y_f')

      i=max(n1,n2,n3)
      work%aeff0array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0array')
      work%beff0array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0array')
      work%ceff0array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0array')
      work%eeff0array = f_malloc0_ptr((/ lowfil.to.lupfil, 0.to.i /),id='work%eeff0array')

      work%aeff0_2array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0_2array')
      work%beff0_2array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0_2array')
      work%ceff0_2array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0_2array')
      work%eeff0_2array = f_malloc0_ptr((/ lowfil.to.lupfil, 0.to.i /),id='work%eeff0_2array')

      work%aeff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0_2auxarray')
      work%beff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0_2auxarray')
      work%ceff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0_2auxarray')
      work%eeff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%eeff0_2auxarray')

      work%xya_c = f_malloc_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xya_c')
      work%xyc_c = f_malloc_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xyc_c')
      if(with_confpot) then
         call f_zero(work%xya_c)
         call f_zero(work%xyc_c)
      end if

      work%xza_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xza_c')
      work%xzc_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xzc_c')
      if(with_confpot) then
         call f_zero(work%xza_c)
         call f_zero(work%xzc_c)
      end if

      work%yza_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yza_c')
      work%yzb_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yzb_c')
      work%yzc_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yzc_c')
      work%yze_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yze_c')
      if(with_confpot) then
         call f_zero(work%yza_c)
         call f_zero(work%yzb_c)
         call f_zero(work%yzc_c)
         call f_zero(work%yze_c)
      end if

      work%xya_f = f_malloc_ptr((/ 1.to.3, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xya_f')
      work%xyb_f = f_malloc_ptr((/ 1.to.4, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xyb_f')
      work%xyc_f = f_malloc_ptr((/ 1.to.3, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xyc_f')
      work%xye_f = f_malloc_ptr((/ 1.to.4, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xye_f')
      if(with_confpot) then
         call f_zero(work%xya_f)
         call f_zero(work%xyb_f)
         call f_zero(work%xyc_f)
         call f_zero(work%xye_f)
      end if

      work%xza_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xza_f')
      work%xzb_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xzb_f')
      work%xzc_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xzc_f')
      work%xze_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xze_f')
      if(with_confpot) then
         call f_zero(work%xza_f)
         call f_zero(work%xzb_f)
         call f_zero(work%xzc_f)
         call f_zero(work%xze_f)
      end if

      work%yza_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yza_f')
      work%yzb_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yzb_f')
      work%yzc_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yzc_f')
      work%yze_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yze_f')
      if(with_confpot) then
         call f_zero(work%yza_f)
         call f_zero(work%yzb_f)
         call f_zero(work%yzc_f)
         call f_zero(work%yze_f)
      end if


    END SUBROUTINE init_local_work_arrays


    subroutine zero_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work, subname)
      implicit none

      ! Calling arguments
      integer,intent(in)::n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
      logical,intent(in):: with_confpot
      type(workarrays_quartic_convolutions),intent(inout):: work
      character(len=*),intent(in):: subname

      ! Local variables
      integer:: i, istat
      integer,parameter :: lowfil=-14,lupfil=14

      call f_routine(id='zero_local_work_arrays')

      call f_zero(work%xx_c)
      call f_zero(work%xy_c)
      call f_zero(work%xz_c)

      call f_zero(work%xx_f1)
      call f_zero(work%xx_f)


      call f_zero(work%xy_f2)
      call f_zero(work%xy_f)


      call f_zero(work%xz_f4)
      call f_zero(work%xz_f)


      call f_zero(work%y_c)

      call f_zero(work%y_f)

      i=max(n1,n2,n3)
      call f_zero(work%aeff0array)
      call f_zero(work%beff0array)
      call f_zero(work%ceff0array)
      call f_zero(work%eeff0array)

      call f_zero(work%aeff0_2array)
      call f_zero(work%beff0_2array)
      call f_zero(work%ceff0_2array)
      call f_zero(work%eeff0_2array)

      call f_zero(work%aeff0_2auxarray)
      call f_zero(work%beff0_2auxarray)
      call f_zero(work%ceff0_2auxarray)
      call f_zero(work%eeff0_2auxarray)

      if(with_confpot) then
         call f_zero(work%xya_c)
         call f_zero(work%xyc_c)
      end if

      if(with_confpot) then
         call f_zero(work%xza_c)
         call f_zero(work%xzc_c)
      end if

      if(with_confpot) then
         call f_zero(work%yza_c)
         call f_zero(work%yzb_c)
         call f_zero(work%yzc_c)
         call f_zero(work%yze_c)
      end if

      if(with_confpot) then
         call f_zero(work%xya_f)
         call f_zero(work%xyb_f)
         call f_zero(work%xyc_f)
         call f_zero(work%xye_f)
      end if

      if(with_confpot) then
         call f_zero(work%xza_f)
         call f_zero(work%xzb_f)
         call f_zero(work%xzc_f)
         call f_zero(work%xze_f)
      end if

      if(with_confpot) then
         call f_zero(work%yza_f)
         call f_zero(work%yzb_f)
         call f_zero(work%yzc_f)
         call f_zero(work%yze_f)
      end if

      call f_release_routine()

    END SUBROUTINE zero_local_work_arrays

    !> Tranform wavefunction between localisation region and the global region
    !!!!!#######!> This routine only works if both locregs have free boundary conditions.
    !! @warning 
    !! WARNING: Make sure psi is set to zero where Glr does not collide with Llr (or everywhere)
    subroutine Lpsi_to_global2(iproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, lpsi, psi)    
     implicit none
    
      ! Subroutine Scalar Arguments
      integer,intent(in):: iproc
      integer :: Gdim          ! dimension of psi 
      integer :: Ldim          ! dimension of lpsi
      integer :: norb          ! number of orbitals
      integer :: nspinor       ! number of spinors
      integer :: nspin         ! number of spins 
      type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
      type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
      
      !Subroutine Array Arguments
      real(wp),dimension(Gdim),intent(inout) :: psi       !Wavefunction (compressed format)
      real(wp),dimension(Ldim),intent(in) :: lpsi         !Wavefunction in localization region
      
      !local variables
      integer :: igrid,isegloc,isegG,ix!,iorbs
      integer :: lmin,lmax,Gmin,Gmax
      integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
      integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
      integer :: length      ! Length of the overlap between Lseg and Gseg
      integer :: lincrement  ! Increment for writing orbitals in loc_psi
      integer :: Gincrement  ! Increment for reading orbitals in psi
      integer :: nseg        ! total number of segments in Llr
      !integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
      character(len=*), parameter :: subname='Lpsi_to_global'
      integer :: i_all
      integer :: start,Gstart,Lindex
      integer :: lfinc,Gfinc,spinshift,ispin,Gindex,isegstart
      integer :: istart
      !integer :: i_stat
    
      call f_routine(id=subname)
    
      !!! This routine is only intended for conversions between locregs with the same boundary conditions.
      !!if (glr%geocode/= 'F' .or. llr%geocode/='F') then
      !!    call f_err_throw('Lpsi_to_global2 can only be used for locregs with free boundary conditions', &
      !!         err_name='BIGDFT_RUNTIME_ERROR')
      !!end if
    
      if(nspin/=1) stop 'not fully implemented for nspin/=1!'
    
    ! Define integers
      nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
      lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
      Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      icheck = 0
      spinshift = Gdim / nspin
     
    ! Get the keymask: shift for every segment of Llr (with respect to Glr)
    ! allocate(keymask(2,nseg),stat=i_stat)
      !keymask = f_malloc((/2,nseg/),id='keymask')
    
      !call shift_locreg_indexes(Glr,Llr,keymask,nseg)
      !call shift_locreg_indexes_global(Glr,Llr,keymask,nseg)
      !!keymask = llr%wfd%keyglob
    
    !####################################################
    ! Do coarse region
    !####################################################
      isegstart=1
    
     
      !$omp parallel default(private) &
      !$omp shared(Glr,Llr, lpsi,icheck,psi,norb) &
      !$omp firstprivate(isegstart,nseg,lincrement,Gincrement,spinshift,nspin) 
    
      !$omp do reduction(+:icheck)
      local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
         lmin = llr%wfd%keyglob(1,isegloc)
         lmax = llr%wfd%keyglob(2,isegloc)
         istart = llr%wfd%keyvglob(isegloc)-1
    
         
         global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
            Gmin = Glr%wfd%keyglob(1,isegG)
            Gmax = Glr%wfd%keyglob(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            !if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax) exit global_loop_c
    
            !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_c
            if(lmin > Gmax)  cycle global_loop_c
    
            ! Define the offset between the two segments
            offset = lmin - Gmin
            if(offset < 0) then
               offset = 0
            end if
    
            ! Define the length of the two segments
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            !Find the common elements and write them to the new global wavefunction
            icheck = icheck + (length + 1)
    
            ! WARNING: index goes from 0 to length because it is the offset of the element
    
            do ix = 0,length     
               istart = istart + 1
               do ispin=1,nspin
                  Gindex = Glr%wfd%keyvglob(isegG)+offset+ix+spinshift*(ispin-1)
                  Lindex = istart+lincrement*norb*(ispin-1)
                  psi(Gindex) = lpsi(Lindex) 
               end do
            end do
         end do global_loop_c
      end do local_loop_c
      !$omp end do
    
    
    !##############################################################
    ! Now do fine region
    !##############################################################
    
      start = Llr%wfd%nvctr_c
      Gstart = Glr%wfd%nvctr_c
      lfinc  = Llr%wfd%nvctr_f
      Gfinc = Glr%wfd%nvctr_f
    
      isegstart=Glr%wfd%nseg_c+1
    
      !$omp do reduction(+:icheck)
      local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
         lmin = llr%wfd%keyglob(1,isegloc)
         lmax = llr%wfd%keyglob(2,isegloc)
         istart = llr%wfd%keyvglob(isegloc)-1
    
         global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
    
            Gmin = Glr%wfd%keyglob(1,isegG)
            Gmax = Glr%wfd%keyglob(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            ! if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax)  exit global_loop_f
            !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f
            if(lmin > Gmax)  cycle global_loop_f
    
            offset = lmin - Gmin
            if(offset < 0) offset = 0
    
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            !Find the common elements and write them to the new global wavefunction
            ! First set to zero those elements which are not copied. WARNING: will not work for npsin>1!!
     
            icheck = icheck + (length + 1)
    
            ! WARNING: index goes from 0 to length because it is the offset of the element
            do ix = 0,length
            istart = istart + 1
               do igrid=1,7
                  do ispin = 1, nspin
                     Gindex = Gstart + (Glr%wfd%keyvglob(isegG)+offset+ix-1)*7+igrid + spinshift*(ispin-1)
                     Lindex = start+(istart-1)*7+igrid + lincrement*norb*(ispin-1) 
                     psi(Gindex) = lpsi(Lindex) 
                  end do
               end do
            end do
         end do global_loop_f
      end do local_loop_f
      !$omp end do
    
      !$omp end parallel
    
      !Check if the number of elements in loc_psi is valid
      if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
        write(*,*)'There is an error in Lpsi_to_global2: sum of fine and coarse points used',icheck
        write(*,*)'is not equal to the sum of fine and coarse points in the region',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
        stop
      end if
    
      !!call f_free(keymask)
    
      call f_release_routine()
    
    END SUBROUTINE Lpsi_to_global2


    !> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
    !! @warning       
    !!    The quantity must not be stored in a compressed form.
    subroutine global_to_local_parallel(Glr,Llr,size_rho,size_Lrho,rho,Lrho,i1s,i1e,i2s,i2e,i3s,i3e,ni1,ni2, &
               i1shift, i2shift, i3shift, ise)
     implicit none
    
     ! Arguments
     type(locreg_descriptors),intent(in) :: Llr   !< Local localization region
     type(locreg_descriptors),intent(in) :: Glr   !< Global localization region
     integer, intent(in) :: size_rho  ! size of rho array
     integer, intent(in) :: size_Lrho ! size of Lrho array
     real(wp),dimension(size_rho),intent(in) :: rho  !< quantity in global region
     real(wp),dimension(size_Lrho),intent(out) :: Lrho !< piece of quantity in local region
     integer,intent(in) :: i1s, i1e, i2s, i2e
     integer,intent(in) :: i3s, i3e ! starting and ending indices on z direction (related to distribution of rho when parallel)
     integer,intent(in) :: ni1, ni2 ! x and y extent of rho
     integer,intent(in) :: i1shift, i2shift, i3shift
     integer,dimension(6) :: ise
    
    ! Local variable
     integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
     integer :: indSmall, indSpin, indLarge ! indexes for the arrays
     integer :: ist2S,ist3S, ist2L, ist3L, istsa, ists, istl
     integer :: ii1shift, ii2shift, ii3shift, i1glob, i2glob, i3glob
     integer :: iii1, iii2, iii3
    
     !THIS ROUTINE NEEDS OPTIMIZING
    
     !!write(*,'(a,8i8)') 'in global_to_local_parallel: i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2', i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2
     
     ! Cut out a piece of the quantity (rho) from the global region (rho) and
     ! store it in a local region (Lrho).
     indSmall=0
     indSpin=0
     ! Deactivate the spin for the moment
     do ispin=1,1!nspin
         !$omp parallel default(none) &
         !$omp shared(Glr, Llr, Lrho, rho, indSpin, i1s, i1e, i2s, i2e, i3s, i3e) &
         !$omp shared(i1shift, i2shift, i3shift, ni1, ni2, ise) &
         !$omp private(ii1, ii2, ii3, i1glob, i2glob, i3glob, ii1shift, ii2shift, ii3shift) &
         !$omp private(ist3S, ist3L, istsa, ist2S, ist2L, ists, istl, indSmall, indLarge) &
         !$omp private(iii1, iii2, iii3)
         !$omp do
         do ii3=i3s,i3e
             i3glob = ii3+ise(5)-1
             !i3=modulo(i3glob-1,glr%d%n3i)+1
             if (modulo(ii3-1,glr%d%n3i)+1>modulo(i3e-1,glr%d%n3i)+1) then
                 !This is a line before the wrap around, i.e. one needs a shift since 
                 ii3shift = i3shift
             else
                 ii3shift = 0
             end if
             if (i3glob<=glr%d%n3i) then
                 iii3=ii3+i3shift
             else
                 iii3=modulo(i3glob-1,glr%d%n3i)+1
             end if
             ist3S = (ii3-i3s)*Llr%d%n2i*Llr%d%n1i
             ist3L = (iii3-1)*ni2*ni1
             istsa=ist3S-i1s+1
             do ii2=i2s,i2e
                 i2glob = ii2+ise(3)-1
                 !i2=modulo(i2glob-1,glr%d%n2i)+1
                 if (modulo(ii2-1,glr%d%n2i)+1>modulo(i2e-1,glr%d%n2i)+1) then
                     !This is a line before the wrap around, i.e. one needs a shift since 
                     !the potential in the global region starts with the wrapped around part
                     ii2shift = i2shift
                 else
                     ii2shift = 0
                 end if
                 if (i2glob<=glr%d%n2i) then
                     iii2=ii2+i2shift
                 else
                     iii2=modulo(i2glob-1,glr%d%n2i)+1
                 end if
                 ist2S = (ii2-i2s)*Llr%d%n1i 
                 ist2L = (iii2-1)*ni1
                 ists=istsa+ist2S
                 istl=ist3L+ist2L
                 do ii1=i1s,i1e
                     i1glob = ii1+ise(1)-1
                     !i1=modulo(i1glob-1,glr%d%n1i)+1
                     if (modulo(ii1-1,glr%d%n1i)+1>modulo(i1e-1,glr%d%n1i)+1) then
                         !This is a line before the wrap around, i.e. one needs a shift since 
                         !the potential in the global region starts with the wrapped around part
                         ii1shift = i1shift
                     else
                         ii1shift = 0
                     end if
                     if (i1glob<=glr%d%n1i) then
                         iii1=ii1+i1shift
                     else
                         iii1=modulo(i1glob-1,glr%d%n1i)+1
                     end if
                     ! indSmall is the index in the local localization region
                     indSmall=ists+ii1
                     ! indLarge is the index in the global localization region. 
                     indLarge= iii1+istl
                     Lrho(indSmall)=rho(indLarge+indSpin)
                     !write(600+bigdft_mpi%iproc,'(a,14i7,2es16.8)') 'i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, i1shift, i2shift, i3shift, indsmall, indlarge, val, testval', &
                     !    i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, i1shift, i2shift, i3shift, indsmall, indlarge, Lrho(indSmall), real((i1+(i2-1)*glr%d%n1i+(i3-1)*glr%d%n1i*glr%d%n2i),kind=8)
                     !if (abs(Lrho(indSmall)-real((i1+(i2-1)*glr%d%n1i+(i3-1)*glr%d%n1i*glr%d%n2i),kind=8))>1.d-3) then
                     !    write(700+bigdft_mpi%iproc,'(a,11i7,2es16.8)') 'i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, indsmall, indlarge, val, testval', &
                     !        i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, indsmall, indlarge, Lrho(indSmall), real((i1+(i2-1)*glr%d%n1i+(i3-1)*glr%d%n1i*glr%d%n2i),kind=8)
                     !end if
                 end do
             end do
         end do
         !$omp end do
         !$omp end parallel
         indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
     end do
    
    END SUBROUTINE global_to_local_parallel

    function boundary_weight(hgrids,glr,lr,rad,psi) result(weight_normalized)
      implicit none
      real(gp), intent(in) :: rad
      real(gp), dimension(3) :: hgrids
      type(locreg_descriptors), intent(in) :: glr,lr
      real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
      real(gp) :: weight_normalized
      !local variables
      integer :: iorb, iiorb, ilr, iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ind 
      integer :: ij3, ij2, ij1, jj3, jj2, jj1, ijs3, ijs2, ijs1, ije3, ije2, ije1
      real(kind=8) :: h, x, y, z, d, weight_inside, weight_boundary, points_inside, points_boundary, ratio
      real(kind=8) :: boundary
      logical :: perx, pery, perz, on_boundary


      ! mean value of the grid spacing
      h = sqrt(hgrids(1)**2+hgrids(2)**2+hgrids(3)**2)

      ! periodicity in the three directions
      perx=(glr%geocode /= 'F')
      pery=(glr%geocode == 'P')
      perz=(glr%geocode /= 'F')

      ! For perdiodic boundary conditions, one has to check also in the neighboring
      ! cells (see in the loop below)
      if (perx) then
         ijs1 = -1
         ije1 = 1
      else
         ijs1 = 0
         ije1 = 0
      end if
      if (pery) then
         ijs2 = -1
         ije2 = 1
      else
         ijs2 = 0
         ije2 = 0
      end if
      if (perz) then
         ijs3 = -1
         ije3 = 1
      else
         ijs3 = 0
         ije3 = 0
      end if


      boundary = min(rad,lr%locrad)

      weight_boundary = 0.d0
      weight_inside = 0.d0
      points_inside = 0.d0
      points_boundary = 0.d0
      ind = 0
      do iseg=1,lr%wfd%nseg_c
         jj=lr%wfd%keyvglob(iseg)
         j0=lr%wfd%keyglob(1,iseg)
         j1=lr%wfd%keyglob(2,iseg)
         ii=j0-1
         i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
         ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
         i2=ii/(glr%d%n1+1)
         i0=ii-i2*(glr%d%n1+1)
         i1=i0+j1-j0
         do i=i0,i1
            ind = ind + 1
            on_boundary = .false.
            do ij3=ijs3,ije3!-1,1
               jj3=i3+ij3*(glr%d%n3+1)
               z = real(jj3,kind=8)*hgrids(3)
               do ij2=ijs2,ije2!-1,1
                  jj2=i2+ij2*(glr%d%n2+1)
                  y = real(jj2,kind=8)*hgrids(2)
                  do ij1=ijs1,ije1!-1,1
                     jj1=i+ij1*(glr%d%n1+1)
                     x = real(i,kind=8)*hgrids(1)
                     d = sqrt((x-lr%locregcenter(1))**2 + &
                          (y-lr%locregcenter(2))**2 + &
                          (z-lr%locregcenter(3))**2)
                     if (abs(d-boundary)<h) then
                        on_boundary=.true.
                     end if
                  end do
               end do
            end do
            if (on_boundary) then
               ! This value is on the boundary
               !write(*,'(a,2f9.2,3i8,3es16.8)') 'on boundary: boundary, d, i1, i2, i3, x, y, z', &
               !    boundary, d, i, i2, i3, x, y, z
               weight_boundary = weight_boundary + psi(ind)**2
               points_boundary = points_boundary + 1.d0
            else
               weight_inside = weight_inside + psi(ind)**2
               points_inside = points_inside + 1.d0
            end if
         end do
      end do
      ! fine part, to be done only if nseg_f is nonzero
      do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
         jj=lr%wfd%keyvglob(iseg)
         j0=lr%wfd%keyglob(1,iseg)
         j1=lr%wfd%keyglob(2,iseg)
         ii=j0-1
         i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
         ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
         i2=ii/(glr%d%n1+1)
         i0=ii-i2*(glr%d%n1+1)
         i1=i0+j1-j0
         do i=i0,i1
            ind = ind + 7
            on_boundary = .false.
            do ij3=ijs3,ije3!-1,1
               jj3=i3+ij3*(glr%d%n3+1)
               z = real(jj3,kind=8)*hgrids(3)
               do ij2=ijs2,ije2!-1,1
                  jj2=i2+ij2*(glr%d%n2+1)
                  y = real(jj2,kind=8)*hgrids(2)
                  do ij1=ijs1,ije1!-1,1
                     jj1=i+ij1*(glr%d%n1+1)
                     x = real(i,kind=8)*hgrids(1)
                     d = sqrt((x-lr%locregcenter(1))**2 + &
                          (y-lr%locregcenter(2))**2 + &
                          (z-lr%locregcenter(3))**2)
                     if (abs(d-boundary)<h) then
                        on_boundary=.true.
                     end if
                  end do
               end do
            end do
            if (on_boundary) then
               ! This value is on the boundary
               !write(*,'(a,f9.2,3i8,3es16.8)') 'on boundary: d, i1, i2, i3, x, y, z', d, i, i2, i3, x, y, z
               weight_boundary = weight_boundary + psi(ind-6)**2
               weight_boundary = weight_boundary + psi(ind-5)**2
               weight_boundary = weight_boundary + psi(ind-4)**2
               weight_boundary = weight_boundary + psi(ind-3)**2
               weight_boundary = weight_boundary + psi(ind-2)**2
               weight_boundary = weight_boundary + psi(ind-1)**2
               weight_boundary = weight_boundary + psi(ind-0)**2
               points_boundary = points_boundary + 7.d0
            else
               weight_inside = weight_inside + psi(ind-6)**2
               weight_inside = weight_inside + psi(ind-5)**2
               weight_inside = weight_inside + psi(ind-4)**2
               weight_inside = weight_inside + psi(ind-3)**2
               weight_inside = weight_inside + psi(ind-2)**2
               weight_inside = weight_inside + psi(ind-1)**2
               weight_inside = weight_inside + psi(ind-0)**2
               points_inside = points_inside + 7.d0
            end if
         end do
      end do
      ! Ratio of the points on the boundary with resepct to the total number of points
      ratio = points_boundary/(points_boundary+points_inside)
      weight_normalized = weight_boundary/ratio
      !write(*,'(a,2f9.1,4es16.6)') 'iiorb, pi, pb, weight_inside, weight_boundary, ratio, xi', &
      !    points_inside, points_boundary, weight_inside, weight_boundary, &
      !    points_boundary/(points_boundary+points_inside), &
      !    weight_boundary/ratio
      
    end function boundary_weight

!!$    !> Check the relative weight which the support functions have at the
!!$    !! boundaries of the localization regions.
!!$    subroutine get_boundary_weight(iproc, nproc, orbs, lzd, atoms, crmult, nsize_psi, psi, crit)
!!$      use module_types, only: orbitals_data, local_zone_descriptors
!!$      use module_atoms, only: atoms_data
!!$      use yaml_output
!!$      implicit none
!!$
!!$      ! Calling arguments
!!$      integer,intent(in) :: iproc, nproc
!!$      type(orbitals_data),intent(in) :: orbs
!!$      type(local_zone_descriptors),intent(in) :: lzd
!!$      type(atoms_data),intent(in) :: atoms
!!$      real(kind=8),intent(in) :: crmult
!!$      integer,intent(in) :: nsize_psi
!!$      real(kind=8),dimension(nsize_psi),intent(in) :: psi
!!$      real(kind=8),intent(in) :: crit
!!$
!!$      ! Local variables
!!$      integer :: iorb, iiorb, ilr, iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ind, iat, iatype
!!$      integer :: ij3, ij2, ij1, jj3, jj2, jj1, ijs3, ijs2, ijs1, ije3, ije2, ije1, nwarnings
!!$      real(kind=8) :: h, x, y, z, d, weight_inside, weight_boundary, points_inside, points_boundary, ratio
!!$      real(kind=8) :: atomrad, rad, boundary, weight_normalized, maxweight, meanweight
!!$      real(kind=8),dimension(:),allocatable :: maxweight_types, meanweight_types
!!$      integer,dimension(:),allocatable :: nwarnings_types, nsf_per_type
!!$      logical :: perx, pery, perz, on_boundary
!!$
!!$      call f_routine(id='get_boundary_weight')
!!$
!!$      maxweight_types = f_malloc0(atoms%astruct%ntypes,id='maxweight_types')
!!$      meanweight_types = f_malloc0(atoms%astruct%ntypes,id='maxweight_types')
!!$      nwarnings_types = f_malloc0(atoms%astruct%ntypes,id='nwarnings_types')
!!$      nsf_per_type = f_malloc0(atoms%astruct%ntypes,id='nsf_per_type')
!!$
!!$      if (iproc==0) then
!!$         call yaml_sequence(advance='no')
!!$      end if
!!$
!!$      ! mean value of the grid spacing
!!$      h = sqrt(lzd%hgrids(1)**2+lzd%hgrids(2)**2+lzd%hgrids(3)**2)
!!$
!!$      ! periodicity in the three directions
!!$      perx=(lzd%glr%geocode /= 'F')
!!$      pery=(lzd%glr%geocode == 'P')
!!$      perz=(lzd%glr%geocode /= 'F')
!!$
!!$      ! For perdiodic boundary conditions, one has to check also in the neighboring
!!$      ! cells (see in the loop below)
!!$      if (perx) then
!!$         ijs1 = -1
!!$         ije1 = 1
!!$      else
!!$         ijs1 = 0
!!$         ije1 = 0
!!$      end if
!!$      if (pery) then
!!$         ijs2 = -1
!!$         ije2 = 1
!!$      else
!!$         ijs2 = 0
!!$         ije2 = 0
!!$      end if
!!$      if (perz) then
!!$         ijs3 = -1
!!$         ije3 = 1
!!$      else
!!$         ijs3 = 0
!!$         ije3 = 0
!!$      end if
!!$
!!$      nwarnings = 0
!!$      maxweight = 0.d0
!!$      meanweight = 0.d0
!!$      if (orbs%norbp>0) then
!!$         ind = 0
!!$         do iorb=1,orbs%norbp
!!$            iiorb = orbs%isorb + iorb
!!$            ilr = orbs%inwhichlocreg(iiorb)
!!$
!!$            iat = orbs%onwhichatom(iiorb)
!!$            iatype = atoms%astruct%iatype(iat)
!!$            atomrad = atoms%radii_cf(iatype,1)*crmult
!!$            rad = atoms%radii_cf(atoms%astruct%iatype(iat),1)*crmult
!!$
!!$              boundary = min(rad,lzd%llr(ilr)%locrad)
!!$              !write(*,*) 'rad, locrad, boundary', rad, lzd%llr(ilr)%locrad, boundary
!!$
!!$              nsf_per_type(iatype) = nsf_per_type(iatype ) + 1
!!$
!!$              weight_boundary = 0.d0
!!$              weight_inside = 0.d0
!!$              points_inside = 0.d0
!!$              points_boundary = 0.d0
!!$              do iseg=1,lzd%llr(ilr)%wfd%nseg_c
!!$                 jj=lzd%llr(ilr)%wfd%keyvglob(iseg)
!!$                 j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
!!$                 j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
!!$                 ii=j0-1
!!$                 i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
!!$                 ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
!!$                 i2=ii/(lzd%glr%d%n1+1)
!!$                 i0=ii-i2*(lzd%glr%d%n1+1)
!!$                 i1=i0+j1-j0
!!$                 do i=i0,i1
!!$                    ind = ind + 1
!!$                    on_boundary = .false.
!!$                    do ij3=ijs3,ije3!-1,1
!!$                       jj3=i3+ij3*(lzd%glr%d%n3+1)
!!$                       z = real(jj3,kind=8)*lzd%hgrids(3)
!!$                       do ij2=ijs2,ije2!-1,1
!!$                          jj2=i2+ij2*(lzd%glr%d%n2+1)
!!$                          y = real(jj2,kind=8)*lzd%hgrids(2)
!!$                          do ij1=ijs1,ije1!-1,1
!!$                             jj1=i+ij1*(lzd%glr%d%n1+1)
!!$                             x = real(i,kind=8)*lzd%hgrids(1)
!!$                             d = sqrt((x-lzd%llr(ilr)%locregcenter(1))**2 + &
!!$                                  (y-lzd%llr(ilr)%locregcenter(2))**2 + &
!!$                                  (z-lzd%llr(ilr)%locregcenter(3))**2)
!!$                             if (abs(d-boundary)<h) then
!!$                                on_boundary=.true.
!!$                             end if
!!$                          end do
!!$                       end do
!!$                    end do
!!$                    if (on_boundary) then
!!$                       ! This value is on the boundary
!!$                       !write(*,'(a,2f9.2,3i8,3es16.8)') 'on boundary: boundary, d, i1, i2, i3, x, y, z', &
!!$                       !    boundary, d, i, i2, i3, x, y, z
!!$                       weight_boundary = weight_boundary + psi(ind)**2
!!$                       points_boundary = points_boundary + 1.d0
!!$                    else
!!$                       weight_inside = weight_inside + psi(ind)**2
!!$                       points_inside = points_inside + 1.d0
!!$                    end if
!!$                 end do
!!$              end do
!!$              ! fine part, to be done only if nseg_f is nonzero
!!$              do iseg=lzd%llr(ilr)%wfd%nseg_c+1,lzd%llr(ilr)%wfd%nseg_c+lzd%llr(ilr)%wfd%nseg_f
!!$                 jj=lzd%llr(ilr)%wfd%keyvglob(iseg)
!!$                 j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
!!$                 j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
!!$                 ii=j0-1
!!$                 i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
!!$                 ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
!!$                 i2=ii/(lzd%glr%d%n1+1)
!!$                 i0=ii-i2*(lzd%glr%d%n1+1)
!!$                 i1=i0+j1-j0
!!$                 do i=i0,i1
!!$                    ind = ind + 7
!!$                    on_boundary = .false.
!!$                    do ij3=ijs3,ije3!-1,1
!!$                       jj3=i3+ij3*(lzd%glr%d%n3+1)
!!$                       z = real(jj3,kind=8)*lzd%hgrids(3)
!!$                       do ij2=ijs2,ije2!-1,1
!!$                          jj2=i2+ij2*(lzd%glr%d%n2+1)
!!$                          y = real(jj2,kind=8)*lzd%hgrids(2)
!!$                          do ij1=ijs1,ije1!-1,1
!!$                             jj1=i+ij1*(lzd%glr%d%n1+1)
!!$                             x = real(i,kind=8)*lzd%hgrids(1)
!!$                             d = sqrt((x-lzd%llr(ilr)%locregcenter(1))**2 + &
!!$                                  (y-lzd%llr(ilr)%locregcenter(2))**2 + &
!!$                                  (z-lzd%llr(ilr)%locregcenter(3))**2)
!!$                             if (abs(d-boundary)<h) then
!!$                                on_boundary=.true.
!!$                             end if
!!$                          end do
!!$                       end do
!!$                    end do
!!$                    if (on_boundary) then
!!$                       ! This value is on the boundary
!!$                       !write(*,'(a,f9.2,3i8,3es16.8)') 'on boundary: d, i1, i2, i3, x, y, z', d, i, i2, i3, x, y, z
!!$                       weight_boundary = weight_boundary + psi(ind-6)**2
!!$                       weight_boundary = weight_boundary + psi(ind-5)**2
!!$                       weight_boundary = weight_boundary + psi(ind-4)**2
!!$                       weight_boundary = weight_boundary + psi(ind-3)**2
!!$                       weight_boundary = weight_boundary + psi(ind-2)**2
!!$                       weight_boundary = weight_boundary + psi(ind-1)**2
!!$                       weight_boundary = weight_boundary + psi(ind-0)**2
!!$                       points_boundary = points_boundary + 7.d0
!!$                    else
!!$                       weight_inside = weight_inside + psi(ind-6)**2
!!$                       weight_inside = weight_inside + psi(ind-5)**2
!!$                       weight_inside = weight_inside + psi(ind-4)**2
!!$                       weight_inside = weight_inside + psi(ind-3)**2
!!$                       weight_inside = weight_inside + psi(ind-2)**2
!!$                       weight_inside = weight_inside + psi(ind-1)**2
!!$                       weight_inside = weight_inside + psi(ind-0)**2
!!$                       points_inside = points_inside + 7.d0
!!$                    end if
!!$                 end do
!!$              end do
!!$              ! Ratio of the points on the boundary with resepct to the total number of points
!!$              ratio = points_boundary/(points_boundary+points_inside)
!!$              weight_normalized = weight_boundary/ratio
!!$              meanweight = meanweight + weight_normalized
!!$              maxweight = max(maxweight,weight_normalized)
!!$              meanweight_types(iatype) = meanweight_types(iatype) + weight_normalized
!!$              maxweight_types(iatype) = max(maxweight_types(iatype),weight_normalized)
!!$              if (weight_normalized>crit) then
!!$                 nwarnings = nwarnings + 1
!!$                 nwarnings_types(iatype) = nwarnings_types(iatype) + 1
!!$              end if
!!$              !write(*,'(a,i7,2f9.1,4es16.6)') 'iiorb, pi, pb, weight_inside, weight_boundary, ratio, xi', &
!!$              !    iiorb, points_inside, points_boundary, weight_inside, weight_boundary, &
!!$              !    points_boundary/(points_boundary+points_inside), &
!!$              !    weight_boundary/ratio
!!$           end do
!!$           if (ind/=nsize_psi) then
!!$              call f_err_throw('ind/=nsize_psi ('//trim(yaml_toa(ind))//'/='//trim(yaml_toa(nsize_psi))//')', &
!!$                   err_name='BIGDFT_RUNTIME_ERROR')
!!$           end if
!!$        end if
!!$
!!$        ! Sum up among all tasks... could use workarrays
!!$        if (nproc>1) then
!!$           call mpiallred(nwarnings, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$           call mpiallred(meanweight, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$           call mpiallred(maxweight, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
!!$           call mpiallred(nwarnings_types, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$           call mpiallred(meanweight_types, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$           call mpiallred(maxweight_types, mpi_max, comm=bigdft_mpi%mpi_comm)
!!$           call mpiallred(nsf_per_type, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$        end if
!!$        meanweight = meanweight/real(orbs%norb,kind=8)
!!$        do iatype=1,atoms%astruct%ntypes
!!$           meanweight_types(iatype) = meanweight_types(iatype)/real(nsf_per_type(iatype),kind=8)
!!$        end do
!!$        if (iproc==0) then
!!$           call yaml_sequence_open('Check boundary values')
!!$           call yaml_sequence(advance='no')
!!$           call yaml_mapping_open(flow=.true.)
!!$           call yaml_map('type','overall')
!!$           call yaml_map('mean / max value',(/meanweight,maxweight/),fmt='(2es9.2)')
!!$           call yaml_map('warnings',nwarnings)
!!$           call yaml_mapping_close()
!!$           do iatype=1,atoms%astruct%ntypes
!!$              call yaml_sequence(advance='no')
!!$              call yaml_mapping_open(flow=.true.)
!!$              call yaml_map('type',trim(atoms%astruct%atomnames(iatype)))
!!$              call yaml_map('mean / max value',(/meanweight_types(iatype),maxweight_types(iatype)/),fmt='(2es9.2)')
!!$              call yaml_map('warnings',nwarnings_types(iatype))
!!$              call yaml_mapping_close()
!!$           end do
!!$           call yaml_sequence_close()
!!$        end if
!!$
!!$        ! Print the warnings
!!$        if (nwarnings>0) then
!!$           if (iproc==0) then
!!$              call yaml_warning('The support function localization radii might be too small, got'&
!!$                   &//trim(yaml_toa(nwarnings))//' warnings')
!!$           end if
!!$        end if
!!$
!!$        call f_free(maxweight_types)
!!$        call f_free(meanweight_types)
!!$        call f_free(nwarnings_types)
!!$        call f_free(nsf_per_type)
!!$
!!$        call f_release_routine()
!!$
!!$      end subroutine get_boundary_weight

    !> Tranform one wavefunction between Global region and localisation region
    subroutine psi_to_locreg2(iproc, ldim, gdim, Llr, Glr, gpsi, lpsi)
    
     implicit none
    
      ! Subroutine Scalar Arguments
      integer,intent(in) :: iproc                  ! process ID
      integer,intent(in) :: ldim          ! dimension of lpsi 
      integer,intent(in) :: gdim          ! dimension of gpsi 
      type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
      type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
      
      !Subroutine Array Arguments
      real(wp),dimension(gdim),intent(in) :: gpsi       !Wavefunction (compressed format)
      real(wp),dimension(ldim),intent(out) :: lpsi   !Wavefunction in localization region
      
      !local variables
      integer :: igrid,isegloc,isegG,ix!,iorbs
      integer :: lmin,lmax,Gmin,Gmax
      integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
      integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
      integer :: length      ! Length of the overlap between Lseg and Gseg
      integer :: lincrement  ! Increment for writing orbitals in loc_psi
      integer :: Gincrement  ! Increment for reading orbitals in psi
      integer :: nseg        ! total number of segments in Llr
      integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
      character(len=*), parameter :: subname='psi_to_locreg'
    !  integer :: i_stat,i_all
      integer :: start,Gstart
      integer :: isegstart,istart
    
      call f_routine(id=subname)
    
    ! Define integers
      nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
      lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
      Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      icheck = 0
    
    ! Initialize loc_psi
      call f_zero(lpsi)
    
    ! Get the keymask: shift for every segment of Llr (with respect to Glr)
    ! allocate(keymask(2,nseg),stat=i_stat)
      keymask = f_malloc((/ 2, nseg /),id='keymask')
    
      call shift_locreg_indexes(Glr,Llr,keymask,nseg)
    
    
    !####################################################
    ! Do coarse region
    !####################################################
      isegstart=1
      icheck = 0
    
    
    !$omp parallel default(private) &
    !$omp shared(icheck,lpsi,gpsi,Glr,Llr,keymask,lincrement,Gincrement,Gstart) &
    !$omp firstprivate(isegstart,nseg)
    
      !$omp do reduction(+:icheck)
      local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
         lmin = keymask(1,isegloc)
         lmax = keymask(2,isegloc)
         istart = llr%wfd%keyvloc(isegloc)-1
     
         global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
            Gmin = Glr%wfd%keygloc(1,isegG)
            Gmax = Glr%wfd%keygloc(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            ! if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax) exit global_loop_c
            if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c
            
            ! Define the offset between the two segments
            offset = lmin - Gmin
            if(offset < 0) then
               offset = 0
            end if
        
            ! Define the length of the two segments
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            icheck = icheck + (length + 1)
     
            !Find the common elements and write them to the new localized wavefunction
            ! WARNING: index goes from 0 to length because it is the offset of the element
    
            do ix = 0,length
               istart = istart + 1
               lpsi(istart) = gpsi(Glr%wfd%keyvloc(isegG)+offset+ix)
            end do
         end do global_loop_c
      end do local_loop_c
      !$omp end do
    
    ! Check if the number of elements in loc_psi is valid
     ! if(icheck .ne. Llr%wfd%nvctr_c) then
       ! write(*,*)'There is an error in psi_to_locreg2: number of coarse points used',icheck
       ! write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
     ! end if
    
    !##############################################################
    ! Now do fine region
    !##############################################################
    
      !icheck = 0
      start = Llr%wfd%nvctr_c
      Gstart = Glr%wfd%nvctr_c
    
      isegstart=Glr%wfd%nseg_c+1
    
      !$omp do reduction(+:icheck)
      local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
         lmin = keymask(1,isegloc)
         lmax = keymask(2,isegloc)
         istart = llr%wfd%keyvloc(isegloc)-1
     
         global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
    
            Gmin = Glr%wfd%keygloc(1,isegG)
            Gmax = Glr%wfd%keygloc(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            ! if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax)  exit global_loop_f
            if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f
    
            offset = lmin - Gmin
            if(offset < 0) offset = 0
    
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            icheck = icheck + (length + 1)
    
            !Find the common elements and write them to the new localized wavefunction
            ! WARNING: index goes from 0 to length because it is the offset of the element
            do ix = 0,length
               istart = istart+1
               do igrid=1,7
                  lpsi(start+(istart-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid)
               end do
            end do
         end do global_loop_f
      end do local_loop_f
      !$omp end do
    
      !$omp end parallel
    
     !! Check if the number of elements in loc_psi is valid
      if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
        write(*,'(a,i0,a,i0)')'process ',iproc,': There is an error in psi_to_locreg: number of fine points used ',icheck
        write(*,'(a,i0)')'is not equal to the number of fine points in the region ',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
      end if
    
    
    
    !  i_all=-product(shape(keymask))*kind(keymask)
    ! deallocate(keymask,stat=i_stat)
      call f_free(keymask)
      call f_release_routine()
    
    END SUBROUTINE psi_to_locreg2



    !> Find the shift necessary for the indexes of every segment of Blr
    !!   to make them compatible with the indexes of Alr. These shifts are
    !!   returned in the array keymask(nseg), where nseg should be the number
    !!   of segments in Blr.
    !! @warning 
    !!   This routine supposes that the region Blr is contained in the region Alr.
    !!   This should always be the case, if we concentrate on the overlap between two regions.
    subroutine shift_locreg_indexes(Alr,Blr,keymask,nseg)
     implicit none
    
    ! Arguments
     type(locreg_descriptors),intent(in) :: Alr,Blr   ! The two localization regions
     integer,intent(in) :: nseg
     integer,intent(out) :: keymask(2,nseg)
    
    ! Local variable
     integer :: iseg      !integer for the loop
     integer :: Bindex    !starting index of segments in Blr
     integer :: x,y,z     !coordinates of start of segments in Blr 
     integer :: shift(3)  !shift between the beginning of the segment in Blr and the origin of Alr
     integer ::  tmp
    
    
     ! This routine is only intended for conversions between locregs with the same boundary conditions.
     if (blr%geocode/='F') then
         call f_err_throw('shift_locreg_indexes can only be used for locregs with free boundary conditions', &
              err_name='BIGDFT_RUNTIME_ERROR')
     end if
    
    !Big loop on all segments
    !$omp parallel do default(private) shared(Blr,nseg,Alr,keymask)
     do iseg=1,nseg
    
    !##########################################
    ! For the Starting index
        Bindex = Blr%wfd%keygloc(1,iseg)
        tmp = Bindex -1
        z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
        tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
        y   = tmp / (Blr%d%n1+1)
        x   = tmp - y * (Blr%d%n1+1)
     
    ! Shift between the beginning of the segment and the start of the Alr region
        shift(1) = x + Blr%ns1 - Alr%ns1
        shift(2) = y + Blr%ns2 - Alr%ns2
        shift(3) = z + Blr%ns3 - Alr%ns3
    
    ! Write the shift in index form
        keymask(1,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1
    
    !######################################
    ! For the ending index
    
        Bindex = Blr%wfd%keygloc(2,iseg)
        tmp = Bindex -1
        z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
        tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
        y   = tmp / (Blr%d%n1+1)
        x   = tmp - y * (Blr%d%n1+1)
    
    ! Shift between the beginning of the segment and the start of the Alr region
        shift(1) = x + Blr%ns1 - Alr%ns1
        shift(2) = y + Blr%ns2 - Alr%ns2
        shift(3) = z + Blr%ns3 - Alr%ns3
    
    ! Write the shift in index form
        keymask(2,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1
     end do
    !$omp end parallel do
    
    END SUBROUTINE shift_locreg_indexes


    !> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
    !! @warning: The quantity must not be stored in a compressed form.
    subroutine global_to_local(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho)
     
     implicit none
    
    ! Arguments
     type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
     type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
     integer, intent(in) :: size_rho  ! size of rho array
     integer, intent(in) :: size_Lrho ! size of Lrho array
     integer, intent(in) :: nspin  !number of spins
     real(wp),dimension(size_rho),intent(in) :: rho  ! quantity in global region
     real(wp),dimension(size_Lrho),intent(out) :: Lrho ! piece of quantity in local region
    
    ! Local variable
     integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
     integer :: indSmall, indSpin, indLarge ! indexes for the arrays
     logical:: z_inside, y_inside, x_inside
     integer:: iz, iy, m
     
    ! Cut out a piece of the quantity (rho) from the global region (rho) and
    ! store it in a local region (Lrho).
    
     if(Glr%geocode == 'F') then
         ! Use loop unrolling here
         indSmall=0
         indSpin=0
         do ispin=1,nspin
             ! WARNING: I added the factors 2.
             do i3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
                 iz=(i3-1)*Glr%d%n2i*Glr%d%n1i
                 do i2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
                     iy=(i2-1)*Glr%d%n1i
                     m=mod(Llr%d%n1i+Llr%nsi1-Llr%nsi1,4)
                     if(m/=0) then
                         do i1=Llr%nsi1+1,Llr%nsi1+m
                            indSmall=indSmall+1
                            indLarge=iz+iy+i1
                            Lrho(indSmall)=rho(indLarge+indSpin)
                         end do
                      end if
                      do i1=Llr%nsi1+1+m,Llr%d%n1i+Llr%nsi1,4
                         Lrho(indSmall+1)=rho(iz+iy+i1+0+indSpin)
                         Lrho(indSmall+2)=rho(iz+iy+i1+1+indSpin)
                         Lrho(indSmall+3)=rho(iz+iy+i1+2+indSpin)
                         Lrho(indSmall+4)=rho(iz+iy+i1+3+indSpin)
                         indSmall=indSmall+4
                      end do
                 end do
             end do
             indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
         end do
     else
         ! General case
         indSmall=0
         indSpin=0
         do ispin=1,nspin
             ! WARNING: I added the factors 2.
             do ii3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
                 i3 = mod(ii3-1,Glr%d%n3i)+1
                 z_inside = (i3>0 .and. i3<=Glr%d%n3i+1)
                 iz=(i3-1)*Glr%d%n2i*Glr%d%n1i
                 do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
                     i2 = mod(ii2-1,Glr%d%n2i)+1
                     y_inside = (i2>0 .and. i2<=Glr%d%n2i+1)
                     iy=(i2-1)*Glr%d%n1i
                     do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                         i1 = mod(ii1-1,Glr%d%n1i)+1 
                         x_inside = (i1 > 0 .and. i1 <= Glr%d%n1i+1)
                         ! indSmall is the index in the local localization region
                         indSmall=indSmall+1
                         !!if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                                       !This initializes the buffers of locreg to zeros if outside the simulation box.
                         !!    i3 <= Glr%d%n3i+1 .and. i2 <= Glr%d%n2i+1 .and. i1 <= Glr%d%n1i+1) then       !Should use periodic image instead... MUST FIX THIS.
                         !!   ! indLarge is the index in the global localization region. 
                         !!   indLarge=(i3-1)*Glr%d%n2i*Glr%d%n1i + (i2-1)*Glr%d%n1i + i1
                         if(z_inside .and. y_inside .and. x_inside) then
                            indLarge=iz+iy+i1
                            Lrho(indSmall)=rho(indLarge+indSpin)
                         else
                            Lrho(indSmall)= 0.0_wp
                         end if
                     end do
                 end do
             end do
             indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
         end do
     end if
    
    END SUBROUTINE global_to_local

    pure subroutine nullify_confpot_data(c)
      use module_defs, only: UNINITIALIZED
      implicit none
      type(confpot_data), intent(out) :: c
      c%potorder=0
      !the rest is not useful
      c%prefac     =UNINITIALIZED(c%prefac)     
      c%hh(1)      =UNINITIALIZED(c%hh(1))      
      c%hh(2)      =UNINITIALIZED(c%hh(2))      
      c%hh(3)      =UNINITIALIZED(c%hh(3))      
      c%rxyzConf(1)=UNINITIALIZED(c%rxyzConf(1))
      c%rxyzConf(2)=UNINITIALIZED(c%rxyzConf(2))
      c%rxyzConf(3)=UNINITIALIZED(c%rxyzConf(3))
      c%ioffset(1) =UNINITIALIZED(c%ioffset(1)) 
      c%ioffset(2) =UNINITIALIZED(c%ioffset(2)) 
      c%ioffset(3) =UNINITIALIZED(c%ioffset(3)) 
      c%damping    =UNINITIALIZED(c%damping)

    end subroutine nullify_confpot_data


    !> apply the potential to the psir wavefunction and calculate potential energy
    subroutine psir_to_vpsi(npot,nspinor,lr,pot,vpsir,epot,confdata,vpsir_noconf,econf)
      use dynamic_memory
      implicit none
      integer, intent(in) :: npot,nspinor
      type(locreg_descriptors), intent(in) :: lr !< localization region of the wavefunction
      real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,npot), intent(in) :: pot
      real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout) :: vpsir
      real(gp), intent(out) :: epot
      type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
      real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout), optional :: vpsir_noconf !< wavefunction with  the potential without confinement applied
      real(gp), intent(out),optional :: econf !< confinement energy
      !local variables
      logical :: confining
      integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation

      call f_routine(id='psir_to_vpsi')

      !write(*,'(a,a4,2l5)') 'in psir_to_vpsi: lr%geocode, present(vpsir_noconf), present(econf)', lr%geocode, present(vpsir_noconf), present(econf)

      epot=0.0_gp
      ishift=(/0,0,0/)
      confining=present(confdata)
      if (confining) confining= (confdata%potorder /=0)

      if (confining) then
         if (lr%geocode == 'F') then
            if (present(vpsir_noconf)) then
               if (.not.present(econf)) stop 'ERROR: econf must be present when vpsir_noconf is present!'
               !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
               call apply_potential_lr_conf_noconf(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                    lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                    ishift,lr%d%n2,lr%d%n3,&
                    nspinor,npot,vpsir,pot,epot,&
                    confdata,lr%bounds%ibyyzz_r,vpsir_noconf,econf)
               !confdata=confdata,ibyyzz_r=lr%bounds%ibyyzz_r,psir_noconf=vpsir_noconf,econf=econf)
            else
               !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
               call apply_potential_lr_conf(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                    lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                    ishift,lr%d%n2,lr%d%n3,&
                    nspinor,npot,vpsir,pot,epot,&
                    confdata,lr%bounds%ibyyzz_r)
               !confdata=confdata,ibyyzz_r=lr%bounds%ibyyzz_r)
            end if
         else
!!!if (present(vpsir_noconf)) then
!!!if (.not.present(econf)) stop 'ERROR: econf must be present when vpsir_noconf is present!'
!!!!call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!!    call apply_potential_lr_conf_noconf_nobounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!!         lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!!         ishift,lr%d%n2,lr%d%n3,&
!!!         nspinor,npot,vpsir,pot,epot,&
!!!         confdata,vpsir_noconf,econf)
!!!         !confdata=confdata)
!!!else
            call apply_potential_lr_conf_nobounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 ishift,lr%d%n2,lr%d%n3,&
                 nspinor,npot,vpsir,pot,epot,&
                 confdata)
            !confdata=confdata)
!!! end if
         end if

      else

         if (lr%geocode == 'F') then
            !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
            call apply_potential_lr_bounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 ishift,lr%d%n2,lr%d%n3,&
                 nspinor,npot,vpsir,pot,epot,&
                 lr%bounds%ibyyzz_r)
            !     ibyyzz_r=lr%bounds%ibyyzz_r)
         else
            !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
            call apply_potential_lr_nobounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 ishift,lr%d%n2,lr%d%n3,&
                 nspinor,npot,vpsir,pot,epot)
         end if
      end if

      call f_release_routine()

    end subroutine psir_to_vpsi
   
    subroutine isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,w,psir,hpsi,ekin,k_strten)
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
                 hx,hy,hz,&      !here the grid spacings are supposed to be equal.  SM: not any more
                 lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
                 lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f, &
                 w%x_c(1,idx),w%x_f(1,idx),&
                 w%y_c(1,idx),w%y_f(1,idx),ekino, &
                 w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),111)
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


end module locreg_operations
