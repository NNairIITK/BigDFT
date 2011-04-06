!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> ART initialize_potential
subroutine initialize_potential( )

  use defs
  use bigdft_forces
  implicit None

  !Local variables
  integer                   :: nat_test
  integer, pointer          :: typa(:)    ! Atomic type
  integer, pointer          :: const_(:)  ! Constraints
  real(kind=8), pointer     :: posa(:)    ! Working positions of the atoms
  real(kind=8),dimension(3) :: boxref_    ! Reference box from posinp file

  call bigdft_init( nat_test, typa, posa, const_, boxref_, boundary, nproc, iproc, my_gnrm )

                                      ! test nat_test and nat
  if ( nat_test /= NATOMS ) then
     write(*,*) "Different number of atoms"
     call end_art()
  end if

  if ( .not. restart ) then                
                                      ! transfer.
     typat(:)   = typa(:)
     pos(:)     = posa(:)
     constr(:)  = const_(:)
  else
                                      ! our box corresponds to that one in restart file
     boxref_(:) = boxref(:)
  end if
                                      ! First force calculation.
  call calcforce( NATOMS, pos, boxref_, force, total_energy, evalf_number, .false. )
 
                                      ! The real box is given by call_bigdft 
  boxref(:) = boxref_(:) 

  deallocate(posa)
  deallocate(typa)
  deallocate(const_)

END SUBROUTINE initialize_potential


!> ART finalise_potential
subroutine finalise_potential( )

  use bigdft_forces

  implicit none  
  
  call bigdft_finalise( )

END SUBROUTINE finalise_potential
