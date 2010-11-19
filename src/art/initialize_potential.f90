!!****f* art/initialize_potential
!! FUNCTION
!!   Initialize the potential
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
subroutine initialize_potential( )

  use defs
  use bigdft_forces
  implicit None

  !Local variables
  integer               :: nat_test
  integer, pointer      :: typa(:)    ! Atomic type
  integer, pointer      :: const_(:)  ! constrains
  real(kind=8), pointer :: posa(:)    ! Working positions of the atoms

  call bigdft_init( nat_test, typa, posa, const_, boxref, boundary, nproc, iproc )

  call geopt_set_verbosity(0)
                                      ! test nat_test and nat
  if ( nat_test /= NATOMS ) stop "Different number of atoms"

  if ( .not. restart ) then                
                    ! transfer.
     typat(:)   = typa(:)
     pos(:)     = posa(:)
     constr(:)  = const_(:)
     
  end if
  deallocate(posa)
  deallocate(typa)
  deallocate(const_)

END SUBROUTINE initialize_potential
!!***


!!****f* bart/finalise_potential
!! FUNCTION
!!   Finalize the potential
!! SOURCE
!!
subroutine finalise_potential( )

  use bigdft_forces

  implicit none  
  
  call bigdft_finalise( )

END SUBROUTINE finalise_potential
!!***
