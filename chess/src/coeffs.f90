!> @file
!!   File containing the routines to calculate the coefficients (eigenvectors) of the Hamiltonian matrix
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


!> Module coefficients
module coeffs

  implicit none

  private

  !> Public routines
  public :: get_coeffs_diagonalization

  contains

    subroutine get_coeffs_diagonalization(iproc, nproc, comm, nfvctr, norbu, norbd, norb, blocksize_pdsyev, &
               smats, smatm, ovrlp, ham, coeff, eval_all, eval_occup, info_coeff)
      use futile
      use wrapper_mpi
      use wrapper_linalg, only: vcopy
      use sparsematrix_base
      use sparsematrix, only: uncompress_matrix2, diagonalizehamiltonian2
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, nfvctr, norbu, norbd, norb, blocksize_pdsyev
      type(sparse_matrix),intent(in) :: smats, smatm
      type(matrices),intent(in) :: ovrlp, ham
      real(mp),dimension(nfvctr,norb),intent(out) :: coeff
      real(mp),dimension(nfvctr*smats%nspin),intent(out) :: eval_all
      real(mp),dimension(norb),intent(out) :: eval_occup
      integer,intent(out) :: info_coeff
    
      ! Local variables
      integer :: ispin, iorb
      real(mp) :: maxdiff
      real(mp),dimension(:),allocatable :: eval
      real(mp),dimension(:,:,:),allocatable :: ovrlp_full, ham_full
    
      ovrlp_full = f_malloc((/nfvctr,nfvctr,smats%nspin/),id='ovrlp_full')
      ham_full = f_malloc((/nfvctr,nfvctr,smatm%nspin/),id='ham_full')
    
      call uncompress_matrix2(iproc, nproc, comm, smats, ovrlp%matrix_compr, ovrlp_full)
      call uncompress_matrix2(iproc, nproc, comm, smatm, ham%matrix_compr, ham_full)
    
      ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
      eval = f_malloc(nfvctr,id='eval')
    
      do ispin=1,smats%nspin
          if (iproc==0) call yaml_map('method','diagonalization')
          call diagonalizeHamiltonian2(iproc, nproc, comm, &
               blocksize_pdsyev, nfvctr, &
               ham_full(1,1,ispin), ovrlp_full(1,1,ispin), eval)
    
          ! Broadcast the results (eigenvectors and eigenvalues) from task 0 to
          ! all other tasks (in this way avoiding that different MPI tasks have different values)
          if (nproc>1) then
              if (iproc==0) call yaml_mapping_open('Cross-check among MPI tasks')
              call mpibcast(ham_full(:,:,ispin), comm=comm, maxdiff=maxdiff)
              if (iproc==0) call yaml_map('max diff of eigenvectors',maxdiff,fmt='(es8.2)')
              call mpibcast(eval, comm=comm, maxdiff=maxdiff)
              if (iproc==0) call yaml_map('max diff of eigenvalues',maxdiff,fmt='(es8.2)')
              if (iproc==0) call yaml_mapping_close()
          end if
    
          ! copy all the eigenvalues
          call vcopy(nfvctr, eval(1), 1, eval_all((ispin-1)*nfvctr+1), 1)
          ! copy the eigenvalues of the occupied states
          if (ispin==1) then
              call vcopy(norbu, eval(1), 1, eval_occup(1), 1)
          else
              call vcopy(norbd, eval(1), 1, eval_occup(norbu+1), 1)
          end if
    
          ! Make sure that the eigenvectors have the same sign on all MPI tasks.
          ! To do so, ensure that the first entry is always positive.
          do iorb=1,nfvctr
              if (ham_full(1,iorb,1)<0.d0) then
                  call dscal(nfvctr, -1.d0, ham_full(1,iorb,1), 1)
              end if
          end do
    
          ! Copy the diagonalized matrix to the coeff array.
          if (smats%nspin/=1) then
              ! Only copy the occupied states
              if (ispin==1) then
                  call vcopy(norbu*nfvctr, ham_full(1,1,ispin), 1, coeff(1,1), 1)
              else if (ispin==2) then
                  call vcopy(norbd*nfvctr, ham_full(1,1,ispin), 1, coeff(1,norbu+1), 1)
              end if
          else
              ! Copy all states
              call vcopy(norb*nfvctr, ham_full(1,1,1), 1, coeff(1,1), 1)
          end if
          info_coeff = 0
      end do
    
      call f_free(eval)
      call f_free(ovrlp_full)
      call f_free(ham_full)
    
    end subroutine get_coeffs_diagonalization

end module coeffs
