!> @file
!! Check the module DIIS
!! @author
!!    Copyright (C) 2015-2015 BigDFT group (TD)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test DIIS to find k of resonant states for a well potential
program DIIS_test
  use diis_sd_optimization
  implicit none
  type(DIIS_obj) :: diis
  real(kind=8) :: V0      !< Depth of the well potential
  real(kind=8) :: L       !< Length of the well potential
  real(kind=8) :: alphaSD !< SD step size
  integer :: idsx         !< Length of DIIS history
  integer :: ndim_psi

  V0 = 1.d0
  L = 1.d0
  idsx = 1
  alphaSD = 0.1d0
  ndim_psi = 1
  !Allocate the diis objects
 ! call DIIS_set(idsx,alphaSD,ndim_psi,ngrpp,diis)
 
 ! call diis_opt(iproc,nproc,ngrp,isgrp,ngrpp,igrpproc,ncomp_grp,ndim_psi,psi,hpsi,diis)

 ! call diis_free(diis)

contains

  function k2(k,V0)
    implicit none
    complex(kind=8), intent(in) :: k
    real(kind=8), intent(in) :: V0
    complex(kind=8) :: k2
    k2 = sqrt(k**2 + 2.0d0*V0)
  end function k2

  function diff_k2(k,V0)
    implicit none
    complex(kind=8), intent(in) :: k
    real(kind=8), intent(in) :: V0
    complex(kind=8) :: diff_k2
    diff_k2 = -2.0*k/sqrt(k**2 + 2.0d0*V0)
  end function diff_k2

  function e_l1(k,V0,L)
    implicit none
    complex(kind=8), intent(in) :: k
    real(kind=8), intent(in) :: V0,L
    complex(kind=8) :: e_l1
    e_l1 = k + cmplx(0.d0,1.d0,kind=8)*k2(k,V0)*tan(0.5d0*L*k2(k,V0))
  end function e_l1

  function diff_e_l1(k,V0,L)
    implicit none
    complex(kind=8), intent(in) :: k
    real(kind=8), intent(in) :: V0,L
    complex(kind=8) :: diff_e_l1
    diff_e_l1 = 1.d0 + cmplx(0.d0,1.d0,kind=8)*diff_k2(k,V0)*tan(0.5d0*L*k2(k,V0)) &
            & + cmplx(0.d0,1.d0,kind=8)*k2(k,V0)*0.5d0*L*diff_k2(k,V0)/cos(0.5d0*L*k2(k,V0))
  end function diff_e_l1

end program DIIS_test
