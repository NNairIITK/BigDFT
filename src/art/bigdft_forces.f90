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
subroutine bigdft_init( nat, typa, posa, boxl, nproc_, me_ )

  implicit none

  !Arguments
  integer,      intent(out) :: nat
  integer,      pointer     :: typa(:)
  real(kind=8), pointer     :: posa(:)
  real(kind=8), intent(out) :: boxl
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

  allocate(posa(3 * nat))
  allocate(typa(nat))
  do i = 1, nat, 1
     posa(i)           = rxyz(1, i) * bohr2ang
     posa(i + nat)     = rxyz(2, i) * bohr2ang
     posa(i + 2 * nat) = rxyz(3, i) * bohr2ang
     typa(i) = at%iatype(i)
  end do
  boxl = at%alat1 * bohr2ang

                                      ! The BigDFT restart structure.
    call init_restart_objects(me, in%iacceleration, at, rst, subname)

  deallocate(rxyz)

  END SUBROUTINE bigdft_init
!!***


!!****f* bigdft_forces/calcforce
!! FUNCTION
!!   Calculation of forces
!! SOURCE
!!
subroutine calcforce( nat, posa, boxl, forca, energy )

    implicit none

  !Arguments
  integer,      intent(in)                            :: nat
  real(kind=8), intent(in),  dimension(3*nat), target :: posa
  real(kind=8), intent(in)                            :: boxl
  real(kind=8), intent(out), dimension(3*nat), target :: forca
  real(kind=8), intent(out)                           :: energy

  !Local variables
  integer :: infocode, i, ierror 
  real(gp) :: fnoise
  real(gp), allocatable :: xcart(:,:), fcart(:,:)
  real(gp) :: tsumx, tsumy, tsumz

                                      ! We transfer acell into 'at'
  at%nat   = nat
  at%alat1 = boxl
  at%alat2 = boxl
  at%alat3 = boxl
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
     call call_bigdft( nproc, me, at, xcart, in, energy, fcart, rst, infocode )

  end if

                                      ! Energy in eV 
  energy = energy * ht2ev

                                      ! We remove the net force.
  tsumx = 0.0_dp 
  tsumy = 0.0_dp 
  tsumz = 0.0_dp

  do i= 1, nat 
     tsumx = tsumx + fcart(1,i)
     tsumy = tsumy + fcart(2,i)
     tsumz = tsumz + fcart(3,i)
  end do

  tsumx = tsumx/Real(nat,dp)
  tsumy = tsumy/Real(nat,dp)
  tsumz = tsumz/Real(nat,dp)

  do i = 1, nat
     fcart(1,i) = fcart(1,i) - tsumx
     fcart(2,i) = fcart(2,i) - tsumy
     fcart(3,i) = fcart(3,i) - tsumz
  end do

                                      ! Forces into ev/ang and in 1D array.
  do i = 1, nat, 1
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
subroutine mingeo( nat, boxl, posa, evalf_number, forca, total_energy )

  implicit none

  !Arguments
  integer,      intent(in)                      :: nat
  real(kind=8), intent(in)                      :: boxl
  real(kind=8), intent(out)                     :: total_energy
  real(kind=8), intent(inout), dimension(3*nat) :: posa
  integer,      intent(out)                     :: evalf_number
  real(kind=8), intent(out),   dimension(3*nat) :: forca

  !Local variables
  integer :: i, ierror
  integer :: infocode
  real(gp), allocatable :: xcart(:,:), fcart(:,:)

                                      ! We transfer acell into at
  at%nat = nat
  at%alat1 = boxl
  at%alat2 = boxl
  at%alat3 = boxl
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
     call call_bigdft( nproc, me, at, xcart, in, total_energy, fcart, rst, infocode )
     call geopt( nproc, me, xcart, at, fcart, total_energy, rst, in, evalf_number )

     in%inputPsiId = 1
     initialised   = .true.
     first_time    = .False.

  else

     if ( .not. initialised ) then
      write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
      write(*,*) "No previous call to bigdft_init(). On strike, refuse to work."
      stop
     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierror)
     call geopt( nproc, me, xcart, at, fcart, total_energy, rst, in, evalf_number )

  end if

  total_energy = total_energy * ht2ev

                                      ! Positions into ang.
  do i = 1, nat, 1
     posa(i)           = xcart(1, i) * bohr2ang
     posa(nat + i)     = xcart(2, i) * bohr2ang
     posa(2 * nat + i) = xcart(3, i) * bohr2ang
  end do

                                      ! Forces into ev/ang.
  do i = 1, nat, 1
     forca(i)           = fcart(1, i) * ht2ev / bohr2ang
     forca(nat + i)     = fcart(2, i) * ht2ev / bohr2ang
     forca(2 * nat + i) = fcart(3, i) * ht2ev / bohr2ang
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

                                      ! finalize memory counting.
  call memocc( 0, 0, 'count', 'stop' )

END SUBROUTINE bigdft_finalise
!!***
end module bigdft_forces
!!***
