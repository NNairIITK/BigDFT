!!****f* bart/initialize_potential
!! FUNCTION
!!   Initialize the potential
!! COPYRIGHT
!!   Copyright (C) 2010 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine initialize_potential()

  use defs
  use bigdft_forces

  implicit none

  integer :: ierr, nat_test
  integer, pointer  :: typa(:) ! Atomic type
  real(8), pointer  :: posa(:)  ! Working positions of the atoms

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  call bigdft_init(nat_test,typa,posa,box, nproc, iproc)
  if (nat_test /= NATOMS) stop "NAT"
  ! todo: test nat_test and nat
  type(:) = typa(:)
  pos(:) = posa(:)

end subroutine initialize_potential
!!***


!!****f* bart/finalise_potential
!! FUNCTION
!!   Finalize the potential
!! SOURCE
!!
subroutine finalise_potential()

  use bigdft_forces

  implicit none  
  
  call bigdft_finalise()
end subroutine finalise_potential
!!***
