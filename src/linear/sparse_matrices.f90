!> @file
!!  Linear version: Handle Sparse Matrices
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine check_matrix_compression(iproc,sparsemat)
  use module_base
  use module_types
  use module_interfaces
  use yaml_output
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: compress_matrix_for_allreduce, uncompressMatrix
  implicit none
  integer,intent(in) :: iproc
  type(sparse_matrix),intent(inout) :: sparsemat
  !Local variables
  integer :: i_stat, i_all, jorb, irow, icol, iseg, ii
  character(len=*),parameter :: subname='check_matrix_compression'
  real(kind=8) :: maxdiff
  real(kind=8), parameter :: tol=1.e-10


  !!allocate(sparsemat%matrix(sparsemat%nfvctr,sparsemat%nfvctr),stat=i_stat)
  !!call memocc(i_stat,sparsemat%matrix,'sparsemat%matrix',subname)

  !!allocate(sparsemat%matrix_compr(sparsemat%nvctr),stat=i_stat)
  !!call memocc(i_stat,sparsemat%matrix_compr,'sparsemat%matrix_compr',subname)

  sparsemat%matrix=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctr/),id='sparsemat%matrix')
  sparsemat%matrix_compr=f_malloc_ptr(sparsemat%nvctr,id='sparsemat%matrix_compr')

  call to_zero(sparsemat%nfvctr**2,sparsemat%matrix(1,1))
  do iseg = 1, sparsemat%nseg
     do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
        call get_indecies(jorb,irow,icol)
        !print *,'irow,icol',irow, icol,test_value_matrix(sparsemat%nfvctr, irow, icol)
        sparsemat%matrix(irow,icol) = test_value_matrix(sparsemat%nfvctr, irow, icol)
     end do
  end do
  
  call compress_matrix_for_allreduce(iproc,sparsemat)

  maxdiff = 0.d0
  do iseg = 1, sparsemat%nseg
     ii=0
     do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
        call get_indecies(jorb,irow,icol)
        maxdiff = max(abs(sparsemat%matrix_compr(sparsemat%keyv(iseg)+ii)&
             -test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff)
        ii=ii+1
     end do
  end do

  if (iproc==0) call yaml_map('Tolerances for this check',tol,fmt='(1pe25.17)')

  if(iproc==0) then
    if (maxdiff > tol) then
       call yaml_warning('COMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
    else
       call yaml_map('Maxdiff for compress', maxdiff,fmt='(1pe25.17)')
    end if
  end if

  call uncompressMatrix(iproc,sparsemat)

  maxdiff = 0.d0
  do iseg = 1, sparsemat%nseg
     do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
        call get_indecies(jorb,irow,icol)
        maxdiff = max(abs(sparsemat%matrix(irow,icol)-test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff) 
     end do
  end do

  if(iproc==0) then
    if (maxdiff > tol) then
       call yaml_warning('UNCOMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
    else
       call yaml_map('Maxdiff for uncompress', maxdiff,fmt='(1pe25.17)')
    end if
  end if

  !!i_all = -product(shape(sparsemat%matrix))*kind(sparsemat%matrix)
  !!deallocate(sparsemat%matrix,stat=i_stat)
  !!call memocc(i_stat,i_all,'sparsemat%matrix',subname)
  !!i_all = -product(shape(sparsemat%matrix_compr))*kind(sparsemat%matrix_compr)
  !!deallocate(sparsemat%matrix_compr,stat=i_stat)
  !!call memocc(i_stat,i_all,'sparsemat%matrix_compr',subname)
  call f_free_ptr(sparsemat%matrix)
  call f_free_ptr(sparsemat%matrix_compr)

contains
   !> define a value for the wavefunction which is dependent of the indices
   function test_value_matrix(norb,iorb,jorb)
      use module_base
      implicit none
      integer, intent(in) :: norb,iorb,jorb
      real(kind=8) :: test_value_matrix

      test_value_matrix = norb*(iorb-1)+jorb
      !print *,iorb,jorb,test_value_matrix
   END FUNCTION test_value_matrix

   subroutine get_indecies(ind,irow,icol)
     implicit none
     integer, intent(in) :: ind
     integer, intent(out) :: irow, icol

     icol = (ind - 1) / sparsemat%nfvctr + 1
     irow = ind - (icol-1)*sparsemat%nfvctr
     !print *,'irow,icol',irow,icol
   END SUBROUTINE get_indecies 
end subroutine check_matrix_compression
