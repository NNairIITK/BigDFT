!!subroutine copy_orthon_data(orthpar_in, orthpar_out, subname)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  type(orthon_data),intent(in):: orthpar_in
!!  type(orthon_data),intent(out):: orthpar_out
!!  character(len=*),intent(in):: subname
!!
!!
!!  orthpar_out%directDiag=orthpar_in%directDiag
!!  orthpar_out%norbpInguess=orthpar_in%norbpInguess
!!  orthpar_out%bsLow=orthpar_in%bsLow
!!  orthpar_out%bsUp=orthpar_in%bsUp
!!  orthpar_out%methOrtho=orthpar_in%methOrtho
!!  orthpar_out%iguessTol=orthpar_in%iguessTol
!!  orthpar_out%methTransformOverlap=orthpar_in%methTransformOverlap
!!  orthpar_out%nItOrtho=orthpar_in%nItOrtho
!!  orthpar_out%blocksize_pdsyev=orthpar_in%blocksize_pdsyev
!!  orthpar_out%blocksize_pdgemm=orthpar_in%blocksize_pdgemm
!!
!!end subroutine copy_orthon_data


!only copying sparsity pattern here, not copying whole matrix
subroutine sparse_copy_pattern_new(sparseMat_in, sparseMat_out, iproc, subname)
  use module_base
  use module_types
  use module_interfaces, except_this_one => sparse_copy_pattern
  implicit none

  ! Calling arguments
  type(sparseMatrix),intent(in):: sparseMat_in
  type(sparseMatrix),intent(out):: sparseMat_out
  integer, intent(in) :: iproc
  character(len=*),intent(in):: subname

  ! Local variables
  integer :: i1
  !integer:: iis1, iie1, iis2, iie2, i2, istat, iall

  call timing(iproc,'sparse_copy','ON')

  !assume sparsemat_out is nullified, take advantage of pointers

sparsemat_out=sparsemat_in

  nullify(sparsemat_out%matrix)
  nullify(sparsemat_out%matrix_compr)

return

  sparsemat_out%nvctr = sparsemat_in%nvctr
  sparsemat_out%nseg = sparsemat_in%nseg
  sparsemat_out%full_dim1 = sparsemat_in%full_dim1
  sparsemat_out%full_dim2 = sparsemat_in%full_dim2

  nullify(sparsemat_out%matrix)
  nullify(sparsemat_out%matrix_compr)

  if(associated(sparsemat_in%keyv)) then
     sparsemat_out%keyv = sparsemat_in%keyv
  end if

  if(associated(sparsemat_in%nsegline)) then
     sparsemat_out%nsegline(i1) = sparsemat_in%nsegline(i1)
  end if

  if(associated(sparsemat_in%istsegline)) then
     sparsemat_out%istsegline = sparsemat_in%istsegline
  end if

  if(associated(sparsemat_in%keyg)) then
     sparsemat_out%keyg = sparsemat_in%keyg
  end if

  !!if(associated(sparsemat_in%matrixindex_in_compressed)) then
  !!   sparsemat_out%matrixindex_in_compressed = sparsemat_in%matrixindex_in_compressed
  !!end if
  if(associated(sparsemat_in%matrixindex_in_compressed_arr)) then
     sparsemat_out%matrixindex_in_compressed_arr = sparsemat_in%matrixindex_in_compressed_arr
  end if

  if(associated(sparsemat_in%orb_from_index)) then
     sparsemat_out%orb_from_index = sparsemat_in%orb_from_index
  end if

  call timing(iproc,'sparse_copy','OF')

end subroutine sparse_copy_pattern_new
