!!****m* bart/bigdft_forces
!! FUNCTION
!!   Module which contains information for bigdft run inside art
!! COPYRIGHT
!!   Copyright (C) 2010 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
module bigdft_forces

  use module_base!, only : gp,wp,dp,bohr2ang
  use module_types
  use module_interfaces
  !use BigDFT_API

  implicit none

  private

  ! Storage of required variables for a SCF loop calculation.
  logical :: initialised = .false.
  integer :: nproc, me
  type(atoms_data) :: at
  type(input_variables) :: in
  type(restart_objects) :: rst
  real(gp), parameter :: ht2ev=27.2113834_gp

  public :: bigdft_init
  public :: calcforce
  public :: mingeo
!!$  public :: bigdft_output
  public :: bigdft_finalise

contains
!!***


!!****f* bigdft_forces/bigdft_init
!! FUNCTION
!!   Routine to initialize all BigDFT stuff
!! SOURCE
!!
  subroutine bigdft_init(nat,typa,posa,boxl, nproc_, me_)

    implicit none

    integer, intent(in) :: nproc_, me_
    integer, intent(out) :: nat
    integer, pointer :: typa(:)
    real*8, pointer :: posa(:)
    real*8, intent(out) :: boxl

    character(len=*), parameter :: subname='bigdft_init'
    real(gp), dimension(:,:), pointer :: rxyz
    real(gp), dimension(:,:), allocatable :: fcart
    integer :: i, infocode, ierror
    real(gp) :: energy

    nproc = nproc_
    me = me_

    ! Initialize memory counting
    call memocc(0,me,'count','start')

    call read_input_variables(me_, "posinp", "input.dft", "input.kpt", &
         & "input.geopt", in, at, rxyz)

    ! Transfer at data to ART variables
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
    
    ! The restart structure.
    call init_restart_objects(me,in%iacceleration,at,rst,subname)

    ! First calculation of forces.
    allocate(fcart(3, at%nat))
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call call_bigdft(nproc,me,at,rxyz,in,energy,fcart,rst,infocode)
    deallocate(rxyz)
    deallocate(fcart)

    in%inputPsiId=1

    initialised = .true.
  end subroutine bigdft_init
!!***

!!****f* bigdft_forces/calcforce
!! FUNCTION
!!   Calculation of forces
!! SOURCE
!!
  subroutine calcforce(nat,typa,posa,boxl,forca,energy)

    implicit none

    integer, intent(in) :: nat
    integer, intent(in), dimension(nat) :: typa
    real*8, intent(in), dimension(3*nat),target :: posa
    real*8, intent(in) :: boxl
    real*8, intent(out), dimension(3*nat),target :: forca
    real*8, intent(out) :: energy

    integer :: infocode, i, ierror
    real(dp), allocatable :: xcart(:,:), fcart(:,:)

    if (.not. initialised) then
       write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
       stop
    end if

    ! We transfer acell into at
    at%nat = nat
    at%alat1 = boxl
    at%alat2 = boxl
    at%alat3 = boxl

    ! need to transform xred into xcart
    allocate(xcart(3, at%nat))
    do i = 1, at%nat, 1
       xcart(:, i) = (/ posa(i), posa(nat + i), posa(2 * nat + i) /) / bohr2ang
    end do

    allocate(fcart(3, at%nat))

    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call call_bigdft(nproc,me,at,xcart,in,energy,fcart,rst,infocode)

    ! need to transform the forces into ev/ang.
    do i = 1, nat, 1
       forca(i)           = fcart(1, i) * ht2ev / bohr2ang
       forca(nat + i)     = fcart(2, i) * ht2ev / bohr2ang
       forca(2 * nat + i) = fcart(3, i) * ht2ev / bohr2ang
    end do
    
    deallocate(xcart)
    deallocate(fcart)
  end subroutine calcforce
!!***


!!****f* bigdft_forces/mingeo
!! FUNCTION
!!   Minimise geometry
!! SOURCE
!!
  subroutine mingeo(nat, boxl, posa, evalf_number, forca, total_energy)

    implicit none

    integer, intent(in) :: nat
    real*8, intent(in) :: boxl
    real*8, intent(out) :: total_energy
    real*8, intent(inout), dimension(3*nat) :: posa
    integer, intent(out) :: evalf_number
    real*8, intent(out), dimension(3*nat) :: forca

    integer :: i, ierror
    real(dp), allocatable :: xcart(:,:), fcart(:,:)

    if (.not. initialised) then
       write(0,*) "No previous call to bigdft_init(). On strike, refuse to work."
       stop
    end if

    ! We transfer acell into at
    at%nat = nat
    at%alat1 = boxl
    at%alat2 = boxl
    at%alat3 = boxl

    ! need to transform xred into xcart
    allocate(xcart(3, at%nat))
    do i = 1, at%nat, 1
       xcart(:, i) = (/ posa(i), posa(nat + i), posa(2 * nat + i) /) / bohr2ang
    end do

    allocate(fcart(3, at%nat))

    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call geopt(nproc,me,xcart,at,fcart,total_energy,rst,in,evalf_number)

    ! need to transform the forces into ev/ang.
    do i = 1, nat, 1
       forca(i)           = fcart(1, i) * ht2ev / bohr2ang
       forca(nat + i)     = fcart(2, i) * ht2ev / bohr2ang
       forca(2 * nat + i) = fcart(3, i) * ht2ev / bohr2ang
    end do
    
    deallocate(xcart)
    deallocate(fcart)
  end subroutine mingeo
!!***


!!****f* bigdft_forces/bigdft_finalise
!! FUNCTION
!!   Routine to finalise all BigDFT stuff
!! SOURCE
!!
  subroutine bigdft_finalise()

    implicit none
    !Local variable
    character(len=*), parameter :: subname='bigdft_finalise'

    call free_restart_objects(rst,subname)
    ! To be completed

    !Finalize memory counting
    call memocc(0,0,'count','stop')

  end subroutine bigdft_finalise


end module bigdft_forces
!!***
