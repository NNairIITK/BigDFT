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
