!> @file
!!    Define module and routines for minimization scheme for ART method
!! @author
!!    June 2001  Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!

!> This module defines a number of parameters used during the minimization
!! for ART method
MODULE minimization
  implicit none
  save

  integer, parameter :: MAX_ITER = 1000
  real(8), parameter :: FTHRESHOLD = 1.0D-1
  real(8), parameter :: STEPSIZE = 0.0001  ! Size in angstroems

  real(8), parameter :: FTHRESH2 = FTHRESHOLD * FTHRESHOLD
END MODULE minimization


!> Minimizes the energy at constant volume. It uses a steepest descent
!! algorithm which is fast and precise enough for our needs.
!!
!! This minimization is done with only a minimal knowledge of the physics
!! of the problem so that it is portable
!!
!! The minimization uses a simple steepest descent with variable step size.
subroutine min_converge
  use defs
  use minimization
  use bigdft_forces
  implicit none
 
!!$  integer :: iter, i, npart
!!$  real(8) :: current_energy, ftot,ftot2, step, boxl, delr
!!$  real(8), dimension(:),allocatable :: posb, forceb
!!$
!!$  allocate(posb(vecsize))
!!$  allocate(forceb(vecsize))
!!$
!!$  evalf_number = 0   ! Set the number of force evaluations to zero
!!$
!!$  ! We compute at constant volume
!!$  boxl = box * scala
!!$
!!$  call calcforce(NATOMS,pos,boxl,force,total_energy)
!!$  evalf_number = evalf_number + 1
!!$  current_energy = total_energy
!!$  if (iproc .eq. 0 ) write(*,*)  'current energy :', total_energy
!!$
!!$  step = STEPSIZE
!!$  do iter = 1, MAX_ITER
!!$  
!!$    posb = pos + step * force
!!$    call calcforce(NATOMS,posb,boxl,forceb,total_energy)
!!$    evalf_number = evalf_number + 1
!!$    
!!$    ftot2 = 0.0d0
!!$    do i=1, VECSIZE
!!$      ftot2 = ftot2 + force(i) * force(i)
!!$    end do
!!$
!!$    if(total_energy < current_energy ) then
!!$      pos = posb
!!$      force = forceb
!!$      step = 1.2 * step
!!$      current_energy = total_energy
!!$
!!$      call displacement(posref, pos, delr,npart)
!!$      if ( mod(iter,1)  == 0 ) then
!!$         if (iproc .eq. 0 ) write(*, "(' ','it: ',i5,' ener: ', f12.4,' ftot: ', f12.6, ' step: ', &
!!$          & f12.6,' evalf: ',i4,' delr: ', f12.6,' npart: ', i4)") iter, &
!!$          & total_energy, sqrt(ftot2),step, evalf_number, delr, npart 
!!$      endif
!!$    else 
!!$      step = 0.6 * step
!!$    endif
!!$
!!$    if(ftot2 < FTHRESH2) exit
!!$  end do
!!$
!!$  ftot = sqrt(ftot2)  
!!$  if (iproc .eq. 0 ) then 
!!$     if (ftot < FTHRESHOLD ) then
!!$        write(*,*) 'Minimization successful   ftot : ', ftot
!!$     else
!!$        write(*,*) 'Minimization failed   ftot : ', ftot
!!$     endif
!!$  endif
!!$
!!$  deallocate(posb)
!!$  deallocate(forceb)
  call mingeo(NATOMS, box * scala, pos, evalf_number, force, total_energy)
END SUBROUTINE min_converge
