!!****m* art/bigdft_forces
!! FUNCTION
!!   Module which contains information for bigdft run inside art
!!
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
module bigdft_forces

  use module_base!, only : gp,wp,dp,bohr2ang
  use module_types
  use module_interfaces
  use defs, only : iproc
  implicit none

  private

  ! Storage of required variables for a SCF loop calculation.
  logical :: initialised = .false.
  logical :: first_time  = .True.
  integer :: nproc, me
  type(atoms_data) :: at
  type(input_variables) :: in
  type(restart_objects) :: rst
  real(gp), parameter :: ht2ev = 27.2113834_gp

  public :: bigdft_init
  public :: calcforce
  public :: mingeo
  public :: bigdft_finalise

contains
!!***


!!****f* bigdft_forces/bigdft_init
!! FUNCTION
!!   Routine to initialize all BigDFT stuff
!! SOURCE
!!
subroutine bigdft_init( nat, typa, posa, const_, boxl, boxt, nproc_, me_ )

  implicit none

  !Arguments
  integer,      intent(out) :: nat
  integer,      pointer     :: typa(:)
  real(kind=8), pointer     :: posa(:)
  integer,      pointer     :: const_(:)
  real(kind=8), dimension(3), intent(out) :: boxl
  character(len=1), intent(out) :: boxt
  integer,      intent(in)  :: nproc_
  integer,      intent(in)  :: me_

  !Local variables
  character(len=*), parameter :: subname='bigdft_init'
  real(gp), dimension(:,:), pointer     :: rxyz
  real(gp), dimension(:,:), allocatable :: fcart
  integer  :: i, infocode, ierror
    real(gp) :: energy,fnoise

  nproc = nproc_
  me = me_
                                      ! Initialize memory counting.
                                      ! Obsolete !!
!  call memocc( 0, me, 'count', 'start' )

                                      ! Read inputs.
  call read_input_variables( me_, "posinp", "input.dft",  &
       & "input.kpt", "input.geopt", "input.perf", in, at, rxyz )

                                      ! Transfer at data to ART variables.
  nat = at%nat
  boxt = at%geocode

  allocate(posa(3 * nat))
  allocate(typa(nat))
  allocate(const_(nat))
 
  do i = 1, nat, 1
     posa(i)           = rxyz(1, i) * bohr2ang
     posa(i + nat)     = rxyz(2, i) * bohr2ang
     posa(i + 2 * nat) = rxyz(3, i) * bohr2ang
     typa(i) = at%iatype(i)
  end do

  boxl(1) = at%alat1 * bohr2ang
  boxl(2) = at%alat2 * bohr2ang
  boxl(3) = at%alat3 * bohr2ang
                                      ! Blocked atoms 
  const_ = 0                          ! Initialization, everyone is free.

  do i = 1, at%nat, 1
     const_(i) = at%ifrztyp(i)  
     !if (at%ifrztyp(i) == 1) then     ! all the atoms' components are frozen.
     !   f_mask(i) = 1 
     !   f_mask(nat + i) = 1
     !   f_mask(2 * nat + i) = 1 
     !else if (at%ifrztyp(i) == 2) then ! atom can only move in the XZ plane 
     !   f_mask(nat + i) = 1
     !else if (at%ifrztyp(i) == 3) then ! atom can only move along the Y axis.
     !   f_mask(i) = 1 
     !   f_mask(2 * nat + i) = 1 
     !end if
  end do
                                      ! The BigDFT restart structure.
  call init_restart_objects(me, in%iacceleration, at, rst, subname)

  END SUBROUTINE bigdft_init
!!***


!!****f* bigdft_forces/calcforce
!! FUNCTION
!!   Calculation of forces
!! SOURCE
!!
subroutine calcforce( nat, posa, boxl, forca, energy, evalf_number )

  implicit none

  !Arguments
  integer,      intent(in)                            :: nat
  real(kind=8), intent(in),  dimension(3*nat), target :: posa
  real(kind=8), dimension(3), intent(in)              :: boxl
  real(kind=8), intent(out), dimension(3*nat), target :: forca
  real(kind=8), intent(out)                           :: energy
  integer,      intent(inout)                         :: evalf_number

  !Local variables
  integer :: infocode, i, ierror 
  real(gp) :: fnoise
  real(gp), allocatable :: xcart(:,:), fcart(:,:)
  real(gp) :: tsumx, tsumy, tsumz

                                      ! We transfer acell into 'at'
  at%nat   = nat
  at%alat1 = boxl(1)/bohr2ang
  at%alat2 = boxl(2)/bohr2ang
  at%alat3 = boxl(3)/bohr2ang
                                      ! Need to transform posa into xcart
                                      ! 1D -> 2D array
  allocate(xcart(3, at%nat))
  do i = 1, at%nat, 1
     xcart(:, i) = (/ posa(i), posa(nat + i), posa(2 * nat + i) /) / bohr2ang
  end do

  allocate(fcart(3, at%nat))

  if ( first_time ) then

     in%inputPsiId = 0
     call MPI_Barrier(MPI_COMM_WORLD,ierror)
     call call_bigdft( nproc, me, at, xcart, in, energy, fcart, fnoise, rst, infocode )
     evalf_number = evalf_number + 1
     
     in%inputPsiId = 1
     initialised   = .true.
     first_time    = .False.

  else 

     if ( .not. initialised ) then
      write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
      write(*,*) "No previous call to bigdft_init(). On strike, refuse to work."
      stop
     end if

                                      ! Get into BigDFT
     call MPI_Barrier(MPI_COMM_WORLD,ierror)
     call call_bigdft( nproc, me, at, xcart, in, energy, fcart, fnoise, rst, infocode )
     evalf_number = evalf_number + 1

  end if
                                      ! Energy in eV 
  energy = energy * ht2ev
                                      ! clean the forces for blocked atoms 
  do i = 1, nat, 1
     if ( at%ifrztyp(i) .ne. 0 ) fcart(:,i) = 0.0d0 
  end do 

  call center_f( fcart, nat )         ! We remove the net force.

  do i = 1, nat, 1                    ! Forces into ev/ang and in 1D array.
     forca( i )           = fcart(1, i) * ht2ev / bohr2ang
     forca( nat + i )     = fcart(2, i) * ht2ev / bohr2ang
     forca( 2 * nat + i ) = fcart(3, i) * ht2ev / bohr2ang
  end do

  deallocate(xcart)
  deallocate(fcart)

END SUBROUTINE calcforce
!!***


!!****f* bigdft_forces/mingeo
!! FUNCTION
!!   Minimise geometry
!! SOURCE
!!
subroutine mingeo( nat, boxl, posa, evalf_number, total_energy, success )

  implicit none

  !Arguments
  integer,      intent(in)                      :: nat
  real(kind=8), intent(in), dimension(3)        :: boxl
  real(kind=8), intent(inout), dimension(3*nat) :: posa
  integer,      intent(inout)                   :: evalf_number
  real(kind=8), intent(out)                     :: total_energy
  logical,      intent(out)                     :: success

  !Local variables
  integer :: i, ierror, ncount_bigdft
  integer :: infocode
  real(gp), allocatable :: xcart(:,:), fcart(:,:)

  success = .True.                    ! success will be .False. if:
                                      !ncount_bigdft > in%ncount_cluster_x

  at%nat   = nat                      ! We transfer acell into at
  at%alat1 = boxl(1)/bohr2ang
  at%alat2 = boxl(2)/bohr2ang
  at%alat3 = boxl(3)/bohr2ang
                                      ! Need to transform posa into xcart
                                      ! 1D -> 2D array
  allocate(xcart(3, at%nat))
  do i = 1, at%nat, 1
     xcart(:, i) = (/ posa(i), posa(nat + i), posa(2 * nat + i) /) / bohr2ang
  end do

  allocate(fcart(3, at%nat))

  if ( first_time ) then

     in%inputPsiId = 0
     call MPI_Barrier(MPI_COMM_WORLD,ierror)
     call call_bigdft( nproc, me, at, xcart, in, total_energy, fcart, fnoise, rst, infocode )
     evalf_number = evalf_number + 1
     call geopt( nproc, me, xcart, at, fcart, total_energy, rst, in, ncount_bigdft )
     evalf_number = evalf_number + ncount_bigdft 

     in%inputPsiId = 1
     initialised   = .true.
     first_time    = .False.

     if (ncount_bigdft > in%ncount_cluster_x) success = .False.

  else

     if ( .not. initialised ) then
        write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
        write(*,*) "No previous call to bigdft_init(). On strike, refuse to work."
        stop
     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierror)
     call geopt( nproc, me, xcart, at, fcart, total_energy, rst, in, ncount_bigdft )
     evalf_number = evalf_number + ncount_bigdft 
     if (ncount_bigdft > in%ncount_cluster_x) success = .False.

  end if

  total_energy = total_energy * ht2ev
                                      ! Positions into ang.
  do i = 1, nat, 1
     posa(i)           = xcart(1, i) * bohr2ang
     posa(nat + i)     = xcart(2, i) * bohr2ang
     posa(2 * nat + i) = xcart(3, i) * bohr2ang
  end do
  
  deallocate(xcart)
  deallocate(fcart)

END SUBROUTINE mingeo
!!***


!!****f* bigdft_forces/bigdft_finalise
!! FUNCTION
!!   Routine to finalise all BigDFT stuff
!! SOURCE
!!
subroutine bigdft_finalise ( )

  implicit none

  !Local variable
  character(len=*), parameter :: subname='bigdft_finalise'
                                      ! Warning: for what this ??
  call free_restart_objects ( rst, subname )
                                      ! Warning, there is this note of Damian :
                                      ! To be completed
                                      ! but what ??? 

  call memocc( 0, 0, 'count', 'stop' )  ! finalize memory counting.

  !deallocate(f_mask)

END SUBROUTINE bigdft_finalise
!!***

!!****f* bigdft_forces/center_f
!! FUNCTION
!!   Removes the net force taking into account the blocked atoms
!!
!! SOURCE
!!
subroutine center_f( vector, natoms )

  implicit none

  !Arguments
  integer, intent(in) :: natoms 
  real(kind=8), dimension(3,natoms), intent(inout), target :: vector

  !Local variables
  integer :: i
  integer :: natoms_f                              ! degrees of freedom
!  integer :: natoms_x, natoms_y, natoms_z         ! degrees of freedom 
  real(kind=8) :: xtotal, ytotal, ztotal
!  real(kind=8), dimension(:), pointer :: x, y, z  ! Pointers for coordinates
!  integer, dimension(:), pointer :: mx, my, mz    ! Pointer for mask 
  logical, dimension(natoms) :: mask

!  ! We first set-up pointers for the x, y, z components 
!  x => vector(1:natoms)
!  y => vector(natoms+1:2*natoms)
!  z => vector(2*natoms+1:3*natoms)

!  ! pointers for the mask components
!  mx => f_mask(1:natoms)
!  my => f_mask(natoms+1:2*natoms)
!  mz => f_mask(2*natoms+1:3*natoms)

  ! degrees of freedom 
  mask = at%ifrztyp .eq. 0
  natoms_f = count(mask)

  ! degrees of freedom in each direction
  ! natoms_x = count(mask)
  ! mask = my .eq. 0
  ! natoms_y = count(mask)
  ! mask = mz .eq. 0
  ! natoms_z = count(mask)

  xtotal = 0.0d0
  ytotal = 0.0d0
  ztotal = 0.0d0

  do i = 1, natoms
     xtotal = xtotal + vector(1,i)
     ytotal = ytotal + vector(2,i)
     ztotal = ztotal + vector(3,i)
  enddo 

  ! The average is only over the degrees of freedom in each direction 
  xtotal = xtotal / natoms_f 
  ytotal = ytotal / natoms_f
  ztotal = ztotal / natoms_f

  do i = 1, natoms, 1
     if ( at%ifrztyp(i) == 0 ) then
        vector(1,i) = vector(1,i) - xtotal
        vector(2,i) = vector(2,i) - ytotal
        vector(3,i) = vector(3,i) - ztotal
     end if 
  end do 

  !where( mx == 0 ) 
  !  x = x - xtotal
  !end where
  !where( my == 0 ) 
  !  y = y - ytotal
  !end where
  !where( mz == 0 ) 
  !  z = z - ztotal
  !end where

END SUBROUTINE center_f
!!***
END MODULE bigdft_forces
!!***
