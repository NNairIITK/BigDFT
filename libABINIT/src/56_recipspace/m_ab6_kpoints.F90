!* * Fortran90 source file *
!*
!* Copyright (C) 2008-2011 ABINIT Group (Damien Caliste)
!* All rights reserved.
!*
!* This file is part of the ABINIT software package. For license information,
!* please see the COPYING file in the top-level directory of the ABINIT source
!* distribution.
!*
!*

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

module m_ab6_kpoints

  use defs_basis
  use m_ab6_symmetry

  implicit none

  private

  logical, private, parameter :: AB_DBG = .false.

  public :: kpoints_get_irreductible_zone

  public :: kpoints_get_mp_k_grid
  public :: kpoints_get_auto_k_grid

  public :: kpoints_binding_mp_k_1
  public :: kpoints_binding_mp_k_2
  public :: kpoints_binding_auto_k_1
  public :: kpoints_binding_auto_k_2

contains

  subroutine kpoints_get_irreductible_zone(irrzon, phnons, &
       & n1, n2, n3, nsppol, nspden, symid, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_56_recipspace
!End of the abilint section

    integer, intent(in)   :: symid
    integer, intent(in)   :: n1, n2, n3, nsppol, nspden
    integer, intent(out)  :: irrzon(n1*n2*n3,2,(nspden/nsppol)-3*(nspden/4))
    real(dp), intent(out) :: phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))
    integer, intent(out)  :: errno

    type(symmetry_type), pointer  :: sym

    if (AB_DBG) write(std_err,*) "AB kpoints: call get irreductible zone."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    if (sym%withSpin /= nspden) then
       errno = AB7_ERROR_ARG
       return
    end if

    call irrzg(irrzon, nspden, nsppol, sym%nSym, n1, n2, n3, phnons, &
         & sym%symAfm, sym%sym, sym%transNon)
  end subroutine kpoints_get_irreductible_zone



  subroutine kpoints_binding_mp_k_1(symid, nkpt, ngkpt, &
       & kptrlatt, kptrlen, nshiftk, shiftk, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_56_recipspace
!End of the abilint section

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(in) :: ngkpt(3)
    integer, intent(inout) :: nshiftk
    real(dp), intent(inout) :: shiftk(3, 8)
    real(dp), intent(out) :: kptrlen
    integer, intent(out) :: kptrlatt(3,3)
    integer, intent(out) :: nkpt

    type(symmetry_type), pointer  :: sym
    real(dp) :: kpt(3,1), wkpt(1)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get k grid1."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    ! First, compute the number of kpoints
    kptrlatt(:,:) = 0
    kptrlatt(1,1) = ngkpt(1)
    kptrlatt(2,2) = ngkpt(2)
    kptrlatt(3,3) = ngkpt(3)
    kptrlen = 20.

    call getkgrid(6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, 0, nkpt, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
  end subroutine kpoints_binding_mp_k_1

  subroutine kpoints_binding_mp_k_2(symid, nkpt, kpt, wkpt, &
       & kptrlatt, kptrlen, nshiftk, shiftk, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_56_recipspace
!End of the abilint section

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(inout) :: nshiftk
    real(dp), intent(inout) :: shiftk(3, 8)
    integer, intent(in) :: nkpt
    real(dp), intent(out) :: kpt(3,nkpt), wkpt(nkpt)
    real(dp), intent(inout) :: kptrlen
    integer, intent(inout) :: kptrlatt(3,3)

    type(symmetry_type), pointer  :: sym
    integer :: nkpt_

    if (AB_DBG) write(std_err,*) "AB symmetry: call get k grid2."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    ! Then, we call it again to get the actual values for the k points.
    call getkgrid(6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
  end subroutine kpoints_binding_mp_k_2


  subroutine kpoints_get_mp_k_grid(symid, nkpt, kpt, wkpt, &
       & ngkpt, nshiftk, shiftk, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(in) :: ngkpt(3)
    integer, intent(in) :: nshiftk
    real(dp), intent(in) :: shiftk(3, nshiftk)
    integer, intent(out) :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)

    real(dp) :: kptrlen
    integer :: kptrlatt(3,3)
    integer :: nshiftk_
    real(dp) :: shiftk_(3, 8)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get k grid."

    nshiftk_ = nshiftk
    shiftk_(:,1:nshiftk_) = shiftk(:,:)

    call kpoints_binding_mp_k_1(symid, nkpt, ngkpt, kptrlatt, kptrlen, &
         & nshiftk_, shiftk_, errno)
    if (errno /= AB7_NO_ERROR) return
    allocate(kpt(3, nkpt))
    allocate(wkpt(nkpt))
    call kpoints_binding_mp_k_2(symid, nkpt, kpt, wkpt, &
       & kptrlatt, kptrlen, nshiftk_, shiftk_, errno)
  end subroutine kpoints_get_mp_k_grid



  subroutine kpoints_binding_auto_k_1(symid, nkpt, kptrlatt, kptrlen, &
       & nshiftk, shiftk, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_56_recipspace
!End of the abilint section

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(out) :: nkpt
    real(dp), intent(inout) :: kptrlen
    integer, intent(out) :: nshiftk
    real(dp), intent(out) :: shiftk(3, 8)
    integer, intent(out) :: kptrlatt(3,3)

    type(symmetry_type), pointer  :: sym
    real(dp), allocatable :: kpt(:,:), wkpt(:)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get auto k grid1."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    !  The parameters of the k lattice are not known, compute
    !  kptrlatt, nshiftk, shiftk.
    call testkgrid(sym%bravais,6,kptrlatt,kptrlen,&
         & AB6_MAX_SYMMETRIES,nshiftk,sym%nSym,0,sym%rprimd,&
         & shiftk,sym%symAfm,sym%sym,sym%vacuum)
    if (AB_DBG) write(std_err,*) "AB symmetry: testkgrid -> kptrlatt=", kptrlatt
    
    !LG: the array kpt seems not allocated here.
    !I would classify it as a bug
    !therefore:
    nkpt=0
    allocate(kpt(3, nkpt))
    allocate(wkpt(nkpt))

    call getkgrid(6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, 0, nkpt, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
    deallocate(kpt,wkpt)
    if (AB_DBG) write(std_err,*) "AB symmetry: getkgrid -> nkpt=", nkpt
  end subroutine kpoints_binding_auto_k_1


  subroutine kpoints_binding_auto_k_2(symid, nkpt, kpt, wkpt, kptrlatt, kptrlen, &
       & nshiftk, shiftk, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_56_recipspace
!End of the abilint section

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(in) :: nkpt
    real(dp), intent(out) :: kpt(3,nkpt), wkpt(nkpt)
    real(dp), intent(inout) :: kptrlen
    integer, intent(inout) :: nshiftk
    real(dp), intent(inout) :: shiftk(3, 8)
    integer, intent(inout) :: kptrlatt(3,3)

    type(symmetry_type), pointer  :: sym
    integer :: nkpt_

    if (AB_DBG) write(std_err,*) "AB symmetry: call get auto k grid2."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    ! Then, we call it again to get the actual values for the k points.
    call getkgrid(6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
  end subroutine kpoints_binding_auto_k_2

  subroutine kpoints_get_auto_k_grid(symid, nkpt, kpt, wkpt, &
       & kptrlen, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(out) :: nkpt
    real(dp), intent(in) :: kptrlen
    real(dp), pointer :: kpt(:,:), wkpt(:)

    real(dp) :: kptrlen_
    integer :: kptrlatt(3,3)
    integer :: nshiftk
    real(dp) :: shiftk(3, 8)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get auto k grid."

    kptrlen_ = kptrlen
    call kpoints_binding_auto_k_1(symid, nkpt, kptrlatt, kptrlen_, &
       & nshiftk, shiftk, errno)
    if (errno /= AB7_NO_ERROR) return
    allocate(kpt(3, nkpt))
    allocate(wkpt(nkpt))
    call kpoints_binding_auto_k_2(symid, nkpt, kpt, wkpt, kptrlatt, kptrlen_, &
       & nshiftk, shiftk, errno)
  end subroutine kpoints_get_auto_k_grid

end module m_ab6_kpoints
