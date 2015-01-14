!! @section LICENCE                                                    
!!    Copyright (C) 2015 BigDFT group                                  
!!    This file is distributed under the terms of the                  
!!    GNU General Public License, see ~/COPYING file                   
!!    or http://www.gnu.org/copyleft/gpl.txt .                         
!!    For the list of contributors, see ~/AUTHORS

!   
!   Copyright (C) 1999-2006 David J. Wales
!   This file was part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!*************************************************************************
!
!  Subroutine MORSE_bulk calculates the energy and gradient analytically for 
!  the Morse potential with periodic boundary conditions. The potential has
!  a cutoff and is shifted to make the potential continuous.  It is not smooth
!  though.
!
!  subroutine morse_bulk_wrapper makes calling it from potentials.f a bit
!  simpler.
!
!  use this potential by assigning the atoms the label "M" as with normal morse
!  potential, but additionally pass the keyword BULK
!
!  the options for the potential are passed using the keyword 
!
!  PARAMS rho boxlx boxly boxlz rcut
!  
!  where rho is the inverse width of the well.  boxlx, boxly, boxlz are the box
!  lengths.  And rcut is the cutoff.
!
!*************************************************************************
!
module module_morse_bulk
    use module_base
    implicit none

    !default private
    private

    !parameter variables
    real(gp),save :: rho
    real(gp),save :: rcut
    real(gp),save :: R0
    real(gp),save :: A
    logical, save :: periodic
    logical, save :: use_cutoff
    logical,save :: initialized=.false.

    public morse_bulk_wrapper
    contains 
subroutine init_morse_bulk(rhoIn,rcutIn)
    use module_base
    use yaml_output
    implicit none
    !parameters
    real(gp) :: rhoIn
    real(gp) :: rcutIn
    !local

    initialized=.true.
    rho = rhoIn
    rcut = rcutIn
    R0 = 1.0_gp
    A = 1.0_gp
    periodic = .true.
    use_cutoff = .true.
    if (rcut / R0 < 1.e-1_gp .and. bigdft_mpi%iproc==0) then
       call yaml_warning('morse_bulk: warning the cutoff is very small')
    endif
end subroutine init_morse_bulk
subroutine morse_bulk_wrapper(nat,alat,rxyz, fxyz, epot)
    use module_base
    implicit none 
    !parameter
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    real(gp), intent(in) :: rxyz(3*nat)
    real(gp), intent(out) :: fxyz(3*nat), epot
    !local
    real(gp) :: alatint(3)
    if(.not.initialized)then
        call f_err_throw('Potential "morse_bulk" not initialized',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
    
!    rho = param1
!    boxvec(1) = param2
!    boxvec(2) = param3
!    boxvec(3) = param4
!    rcut = param5
    
!    R0 = 1.0_gp
!    A = 1.0_gp
!    periodic = .true.
!    use_cutoff = .true.
    
    call morse_bulk(rxyz(1),fxyz(1),epot, nat, rho, R0, A, periodic, & 
       alatint, use_cutoff, rcut)
end subroutine morse_bulk_wrapper

      SUBROUTINE MORSE_BULK(X,V,EMORSE, natoms, rho, R0, A, periodic, &
         boxvec, use_cutoff, rcut)
      ! R0 is the position of the bottom of the well
      ! rho is the width of the well and has units of inverse length
      ! A is the energy scale
!      USE commons
      implicit none 
      logical, intent(in) :: periodic, use_cutoff
      integer, intent(in) :: NATOMS
      real(gp), intent(in) :: X(3*NATOMS), rho, R0, A, boxvec(3), rcut
      real(gp), intent(out) :: V(3*NATOMS), EMORSE
      integer ::  J1, J2, J3, J4
      real(gp) :: DIST, R, DUMMY, &
                       RR(NATOMS,NATOMS), &
                       XMUL2, iboxvec(3), dx(3), eshift
!     logical EVAP, evapreject
!     COMMON /EV/ EVAP, evapreject
      if (periodic) iboxvec(:) = 1.0_gp / boxvec(:)

      if (use_cutoff) then
         Eshift = (1.0_gp - exp(rho * (r0 - rcut)))**2 - 1.0_gp
         !write(*,*) "Eshift", eshift, rcut
      endif

!     EVAP=.FALSE.
      V(:) = 0.0_gp
      EMORSE=0.0_gp
      DO J1=1,NATOMS
         J3=3*J1
         RR(J1,J1)=0.0_gp
         DO J2=J1+1,NATOMS
            J4=3*J2
            dx(:) = X(J3-2:j3)-X(J4-2:j4)
            if (periodic) then
               dx = dx - boxvec * nint(dx * iboxvec)
            endif
            dist = max(sqrt(sum(dx**2)), 1.0e-5_gp)

            if (use_cutoff .and. dist.ge.rcut) cycle

            R=exp(RHO*R0-RHO*DIST)
            DUMMY=R*(R-2.0_gp)
            EMORSE=EMORSE+DUMMY - Eshift

!            if (gtest) then
               xmul2 = 2.0_gp*R*(R-1.0_gp)/DIST * A
               V(J3-2:j3) = V(j3-2:j3) - xmul2 * dx
               V(J4-2:j4) = V(j4-2:j4) + xmul2 * dx
!            endif
         ENDDO
      ENDDO
      EMORSE = EMORSE * A

      RETURN
      END subroutine morse_bulk
end module module_morse_bulk
