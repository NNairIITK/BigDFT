!!!!subroutine copy_linearInputParameters_to_linearParameters(ntypes, nlr, input, lin)
!!!!  use module_base
!!!!  use module_types
!!!!  implicit none
!!!!
!!!!  ! Calling arguments
!!!!  integer,intent(in):: ntypes, nlr
!!!!  type(input_variables),intent(in):: input
!!!!  type(linearParameters),intent(out):: lin
!!!!
!!!!  ! Local variables
!!!!  integer:: itype, ilr
!!!!
!!!!  lin%nit_lowaccuracy = input%lin%nit_lowaccuracy
!!!!  lin%nit_highaccuracy = input%lin%nit_highaccuracy
!!!!  lin%nItBasis_lowaccuracy = input%lin%nItBasis_lowaccuracy
!!!!  lin%nItBasis_highaccuracy = input%lin%nItBasis_highaccuracy
!!!!  lin%nItInnerLoop = input%lin%nItInnerLoop
!!!!  lin%convCrit = input%lin%convCrit
!!!!  lin%DIISHistMin = input%lin%DIISHistMin
!!!!  lin%DIISHistMax = input%lin%DIISHistMax
!!!!  lin%alphaDIIS = input%lin%alphaDIIS
!!!!  lin%alphaSD = input%lin%alphaSD
!!!!  lin%nItPrecond = input%lin%nItPrecond
!!!!  lin%locregShape = input%lin%locregShape
!!!!  lin%blocksize_pdsyev = input%lin%blocksize_pdsyev
!!!!  lin%blocksize_pdgemm = input%lin%blocksize_pdgemm
!!!!  lin%nproc_pdsyev = input%lin%nproc_pdsyev
!!!!  lin%nproc_pdgemm = input%lin%nproc_pdgemm
!!!!  lin%methTransformOverlap = input%lin%methTransformOverlap
!!!!  lin%nItOrtho = input%lin%nItOrtho
!!!!  lin%correctionOrthoconstraint = input%lin%correctionOrthoconstraint
!!!!  lin%mixingMethod = input%lin%mixingMethod
!!!!  lin%mixHist_lowaccuracy = input%lin%mixHist_lowaccuracy
!!!!  lin%nItSCCWhenOptimizing_lowaccuracy = input%lin%nItSCCWhenOptimizing_lowaccuracy
!!!!  lin%nItSCCWhenFixed_lowaccuracy = input%lin%nItSCCWhenFixed_lowaccuracy
!!!!  lin%mixHist_highaccuracy = input%lin%mixHist_highaccuracy
!!!!  lin%nItSCCWhenOptimizing_highaccuracy = input%lin%nItSCCWhenOptimizing_highaccuracy
!!!!  lin%nItSCCWhenFixed_highaccuracy = input%lin%nItSCCWhenFixed_highaccuracy
!!!!  lin%alphaMixWhenOptimizing_lowaccuracy = input%lin%alphaMixWhenOptimizing_lowaccuracy
!!!!  lin%alphaMixWhenFixed_lowaccuracy = input%lin%alphaMixWhenFixed_lowaccuracy
!!!!  lin%convCritMix = input%lin%convCritMix
!!!!  lin%alphaMixWhenOptimizing_highaccuracy = input%lin%alphaMixWhenOptimizing_highaccuracy
!!!!  lin%alphaMixWhenFixed_highaccuracy = input%lin%alphaMixWhenFixed_highaccuracy
!!!!  lin%lowaccuray_converged = input%lin%lowaccuray_converged
!!!!  lin%useDerivativeBasisFunctions = input%lin%useDerivativeBasisFunctions
!!!!  lin%ConfPotOrder = input%lin%ConfPotOrder
!!!!  lin%nItInguess = input%lin%nItInguess
!!!!  lin%memoryForCommunOverlapIG = input%lin%memoryForCommunOverlapIG
!!!!  lin%plotBasisFunctions = input%lin%plotBasisFunctions
!!!!  lin%transformToGlobal = input%lin%transformToGlobal
!!!!  lin%norbsPerProcIG = input%lin%norbsPerProcIG
!!!!  lin%mixedmode = input%lin%mixedmode
!!!!  do itype=1,ntypes
!!!!      lin%norbsPerType(itype) = input%lin%norbsPerType(itype)
!!!!      lin%potentialPrefac_lowaccuracy(itype) = input%lin%potentialPrefac_lowaccuracy(itype)
!!!!      lin%potentialPrefac_highaccuracy(itype) = input%lin%potentialPrefac_highaccuracy(itype)
!!!!  end do
!!!!
!!!!  ! Initialize lin%potentialPrefac to some value (will be adjusted later)
!!!!  lin%potentialPrefac=-1.d0
!!!!
!!!!  ! Assign the localization radius to each atom.
!!!!  do ilr=1,nlr
!!!!      lin%locrad(ilr) = input%lin%locrad(ilr)
!!!!  end do
!!!!  
!!!!end subroutine copy_linearInputParameters_to_linearParameters


!!subroutine compressMatrix(norb, mad, mat, lmat)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: norb
!!  type(matrixDescriptors),intent(in):: mad
!!  real(8),dimension(norb**2),intent(in):: mat
!!  real(8),dimension(mad%nvctr),intent(out):: lmat
!!  
!!  ! Local variables
!!  integer:: iseg, jj, jorb, iiorb, jjorb
!!  
!!  
!!  jj=0
!!  do iseg=1,mad%nseg
!!      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
!!          jj=jj+1
!!          lmat(jj)=mat(jorb)
!!      end do
!!  end do
!!  if(jj/=mad%nvctr) then
!!      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
!!      stop
!!  end if
!!  
!!end subroutine compressMatrix


