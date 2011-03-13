!> @file
!!  Routines to initialize the potential
!! @author
!!   Copyright (C) 2010 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!>   Initialize the potential
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
  call geopt_set_verbosity(0)
  if (nat_test /= NATOMS) stop "NAT"
  ! todo: test nat_test and nat
  type(:) = typa(:)
  pos(:) = posa(:)

END SUBROUTINE initialize_potential


!>   Finalize the potential
subroutine finalise_potential()

  use bigdft_forces

  implicit none  
  
  call bigdft_finalise()
END SUBROUTINE finalise_potential
