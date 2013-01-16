!!!subroutine checkLinearParameters(iproc, nproc, lin)
!!!!
!!!! Purpose:
!!!! ========
!!!!  Checks some values contained in the variable lin on errors.
!!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(linearParameters),intent(inout):: lin
!!!
!!!! Local variables
!!!integer:: norbTarget, nprocIG, ierr
!!!
!!!
!!!  if(lin%DIISHistMin>lin%DIISHistMax) then
!!!      if(iproc==0) write(*,'(1x,a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
!!!      & DIISHistMax, but you chose ', lin%DIISHistMin, ' and ', lin%DIISHistMax, '!'
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!  !!if(trim(lin%getCoeff)/='min' .and. trim(lin%getCoeff)/='diag') then
!!!  !!    if(iproc==0) write(*,'(1x,a,a,a)') "ERROR: lin%getCoeff can have the values 'diag' or 'min', &
!!!  !!        & but we found '", trim(lin%getCoeff), "'!"
!!!  !!    call mpi_barrier(mpi_comm_world, ierr)
!!!  !!    stop
!!!  !!end if
!!!
!!!  if(lin%methTransformOverlap<0 .or. lin%methTransformOverlap>2) then
!!!      if(iproc==0) write(*,'(1x,a,i0,a)') 'ERROR: lin%methTransformOverlap must be 0,1 or 2, but you specified ', &
!!!                               lin%methTransformOverlap,'.'
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!  !!if(trim(lin%getCoeff)=='diag') then
!!!  !!    if(trim(lin%diagMethod)/='seq' .and. trim(lin%diagMethod)/='par') then
!!!  !!        if(iproc==0) write(*,'(1x,a,a,a)') "ERROR: lin%diagMethod can have the values 'seq' or 'par', &
!!!  !!            & but we found '", trim(lin%diagMethod), "'!"
!!!  !!        call mpi_barrier(mpi_comm_world, ierr)
!!!  !!        stop
!!!  !!    end if
!!!  !!end if
!!!
!!!  if(lin%confPotOrder/=4 .and. lin%confPotOrder/=6) then
!!!      if(iproc==0) write(*,'(1x,a,i0,a)') 'ERROR: lin%confPotOrder can have the values 4 or 6, &
!!!          & but we found ', lin%confPotOrder, '!'
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!
!!!  ! Determine the number of processes we need for the minimization of the trace in the input guess.
!!!  if(lin%norbsPerProcIG>lin%orbs%norb) then
!!!      norbTarget=lin%orbs%norb
!!!  else
!!!      norbTarget=lin%norbsperProcIG
!!!  end if
!!!  nprocIG=ceiling(dble(lin%orbs%norb)/dble(norbTarget))
!!!  nprocIG=min(nprocIG,nproc)
!!!
!!!  if( nprocIG/=nproc .and. ((lin%methTransformOverlap==0 .and. (lin%blocksize_pdsyev>0 .or. lin%blocksize_pdgemm>0)) .or. &
!!!      (lin%methTransformOverlap==1 .and. lin%blocksize_pdgemm>0)) ) then
!!!      if(iproc==0) then
!!!          write(*,'(1x,a)') 'ERROR: You want to use some routines from scalapack. This is only possible if all processes are &
!!!                     &involved in these calls, which is not the case here.'
!!!          write(*,'(1x,a)') 'To avoid this problem you have several possibilities:'
!!!          write(*,'(3x,a,i0,a)') "-set 'lin%norbsperProcIG' to a value not greater than ",floor(dble(lin%orbs%norb)/dble(nproc)), &
!!!              ' (recommended; probably only little influence on performance)'
!!!          write(*,'(3x,a)') "-if you use 'lin%methTransformOverlap==1': set 'lin%blocksize_pdgemm' to a negative value &
!!!              &(may heavily affect performance)"
!!!          write(*,'(3x,a)') "-if you use 'lin%methTransformOverlap==0': set 'lin%blocksize_pdsyev' and 'lin%blocksize_pdsyev' &
!!!              &to negative values (may very heavily affect performance)"
!!!      end if
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!  if(lin%nproc_pdsyev>nproc) then
!!!      if(iproc==0) write(*,'(1x,a)') 'ERROR: lin%nproc_pdsyev can not be larger than nproc'
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!  if(lin%nproc_pdgemm>nproc) then
!!!      if(iproc==0) write(*,'(1x,a)') 'ERROR: lin%nproc_pdgemm can not be larger than nproc'
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!  if(lin%locregShape/='c' .and. lin%locregShape/='s') then
!!!      if(iproc==0) write(*,*) "ERROR: lin%locregShape must be 's' or 'c'!"
!!!      call mpi_barrier(mpi_comm_world, ierr)
!!!      stop
!!!  end if
!!!
!!!  if(lin%mixedmode .and. .not.lin%useDerivativeBasisFunctions) then
!!!      if(iproc==0) write(*,*) 'WARNING: will set lin%useDerivativeBasisFunctions to true, &
!!!                               &since this is required if lin%mixedmode is true!'
!!!      lin%useDerivativeBasisFunctions=.true.
!!!  end if
!!!
!!!
!!!
!!!end subroutine checkLinearParameters

!!!subroutine deallocateLinear(iproc, lin, lphi, coeff)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Deallocates all array related to the linear scaling version which have not been 
!!!!   deallocated so far.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!
!!!use module_base
!!!use module_types
!!!use module_interfaces, exceptThisOne => deallocateLinear
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc
!!!type(linearParameters),intent(inout):: lin
!!!real(8),dimension(:),pointer,intent(inout):: lphi
!!!real(8),dimension(:,:),pointer,intent(inout):: coeff
!!!
!!!! Local variables
!!!integer:: istat, iall
!!!character(len=*),parameter:: subname='deallocateLinear'
!!!
!!!
!!!  iall=-product(shape(lphi))*kind(lphi)
!!!  deallocate(lphi, stat=istat)
!!!  call memocc(istat, iall, 'lphi', subname)
!!!
!!!  iall=-product(shape(coeff))*kind(coeff)
!!!  deallocate(coeff, stat=istat)
!!!  call memocc(istat, iall, 'coeff', subname)
!!!
!!!  call deallocate_linearParameters(lin, subname)
!!!
!!!end subroutine deallocateLinear

!!subroutine allocateBasicArrays(lin, ntypes)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(linearParameters),intent(inout):: lin
!!  integer, intent(in) :: ntypes
!!  
!!  ! Local variables
!!  integer:: istat
!!  character(len=*),parameter:: subname='allocateBasicArrays'
!!  
!!  allocate(lin%norbsPerType(ntypes), stat=istat)
!!  call memocc(istat, lin%norbsPerType, 'lin%norbsPerType', subname)
!!  
!!  allocate(lin%potentialPrefac(ntypes), stat=istat)
!!  call memocc(istat, lin%potentialPrefac, 'lin%potentialPrefac', subname)
!!
!!  allocate(lin%potentialPrefac_lowaccuracy(ntypes), stat=istat)
!!  call memocc(istat, lin%potentialPrefac_lowaccuracy, 'lin%potentialPrefac_lowaccuracy', subname)
!!
!!  allocate(lin%potentialPrefac_highaccuracy(ntypes), stat=istat)
!!  call memocc(istat, lin%potentialPrefac_highaccuracy, 'lin%potentialPrefac_highaccuracy', subname)
!!
!!  allocate(lin%locrad(lin%nlr),stat=istat)
!!  call memocc(istat,lin%locrad,'lin%locrad',subname)
!!
!!end subroutine allocateBasicArrays

!!!subroutine deallocateBasicArrays(lin)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  
!!!  ! Calling arguments
!!!  type(linearParameters),intent(inout):: lin
!!!  
!!!  ! Local variables
!!!  integer:: i_stat,i_all
!!!  character(len=*),parameter:: subname='deallocateBasicArrays'
!!! 
!!!  if(associated(lin%potentialPrefac)) then
!!!    !print *,'lin%potentialPrefac',associated(lin%potentialPrefac)
!!!    i_all = -product(shape(lin%potentialPrefac))*kind(lin%potentialPrefac)
!!!    !print *,'i_all',i_all
!!!    deallocate(lin%potentialPrefac,stat=i_stat)
!!!    call memocc(i_stat,i_all,'lin%potentialPrefac',subname)
!!!    nullify(lin%potentialPrefac)
!!!  end if 
!!!  if(associated(lin%norbsPerType)) then
!!!    !print *,'lin%norbsPerType',associated(lin%norbsPerType)
!!!    i_all = -product(shape(lin%norbsPerType))*kind(lin%norbsPerType)
!!!    deallocate(lin%norbsPerType,stat=i_stat)
!!!    call memocc(i_stat,i_all,'lin%norbsPerType',subname)
!!!    nullify(lin%norbsPerType)
!!!  end if 
!!!  if(associated(lin%locrad)) then
!!!    !print *,'lin%locrad',associated(lin%locrad)
!!!    i_all = -product(shape(lin%locrad))*kind(lin%locrad)
!!!    deallocate(lin%locrad,stat=i_stat)
!!!    call memocc(i_stat,i_all,'lin%locrad',subname)
!!!    nullify(lin%locrad)
!!!  end if 
!!!
!!!end subroutine deallocateBasicArrays

!!!!> Allocate the coefficients for the linear combinations of the  orbitals and initialize
!!!!! them at random.
!!!!! Do this only on the root, since the calculations to determine coeff are not yet parallelized.
!!!subroutine initCoefficients(iproc, orbs, lin, coeff)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  
!!!  ! Calling arguments
!!!  integer,intent(in):: iproc
!!!  type(orbitals_data),intent(in):: orbs
!!!  type(linearParameters),intent(in):: lin
!!!  real(8),dimension(:,:),pointer,intent(out):: coeff
!!!  
!!!  ! Local variables
!!!  integer:: iorb, jorb, istat
!!!  real:: ttreal
!!!  character(len=*),parameter:: subname='initCoefficients'
!!!  
!!!  
!!!  allocate(coeff(lin%lb%orbs%norb,orbs%norb), stat=istat)
!!!  call memocc(istat, coeff, 'coeff', subname)
!!!  
!!!  call initRandomSeed(0, 1)
!!!  if(iproc==0) then
!!!      do iorb=1,orbs%norb
!!!         do jorb=1,lin%lb%orbs%norb
!!!            call random_number(ttreal)
!!!            coeff(jorb,iorb)=real(ttreal,kind=8)
!!!         end do
!!!      end do
!!!  end if
!!!
!!!end subroutine initCoefficients

!!!!subroutine create_new_locregs(iproc, nproc, nlr, hx, hy, hz, lorbs, glr, locregCenter, locrad, nscatterarr, withder, &
!!!!           inwhichlocreg_reference, ldiis, &
!!!!           lphilarge, lhphilarge, lhphilargeold, lphilargeold,tmb)
!!!!use module_base
!!!!use module_types
!!!!use module_interfaces, except_this_one => create_new_locregs
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, nlr
!!!!real(8),intent(in):: hx, hy, hz
!!!!type(orbitals_data),intent(in):: lorbs
!!!!type(locreg_descriptors),intent(in):: glr
!!!!real(8),dimension(3,nlr),intent(in):: locregCenter
!!!!real(8),dimension(nlr):: locrad
!!!!integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!!logical,intent(in):: withder
!!!!integer,dimension(lorbs%norb),intent(in):: inwhichlocreg_reference
!!!!type(localizedDIISParameters),intent(inout):: ldiis
!!!!real(8),dimension(:),pointer,intent(out):: lphilarge, lhphilarge, lhphilargeold, lphilargeold
!!!!type(DFT_wavefunction),intent(out):: tmb
!!!!
!!!!! Local variables
!!!!integer:: tag, norbu, norbd, nspin, iorb, iiorb, ilr, npsidim, ii, istat, iall, ierr, norb
!!!!integer,dimension(:),allocatable:: orbsperlocreg
!!!!character(len=*),parameter:: subname='create_new_locregs'
!!!!logical:: reallocate
!!!!
!!!!
!!!!   if(iproc==0) write(*,'(x,a)') 'creating new locregs...'
!!!!
!!!!   call update_locreg(iproc, nproc, nlr, locrad, inwhichlocreg_reference, locregCenter, glr, &
!!!!        withder, nscatterarr, hx, hy, hz, &
!!!!        lorbs, tmb%lzd, tmb%orbs, tmb%op, tmb%comon, tmb%comgp, tmb%comsr, tmb%mad, tmb%collcom)
!!!!   if(withder) stop 'withder is true'
!!!!
!!!!
!!!!end subroutine create_new_locregs

!!!subroutine enlarge_locreg(iproc, nproc, hx, hy, hz, withder, lzd, locrad, &
!!!           ldiis, denspot, tmb)
!!!  use module_base
!!!  use module_types
!!!  use module_interfaces, except_this_one => enlarge_locreg
!!!  implicit none
!!!  
!!!  ! Calling arguments
!!!  integer,intent(in):: iproc, nproc
!!!  real(8),intent(in):: hx, hy, hz
!!!  logical,intent(in):: withder
!!!  type(local_zone_descriptors),intent(inout):: lzd
!!!  real(8),dimension(lzd%nlr),intent(in):: locrad
!!!  type(localizedDIISParameters),intent(inout):: ldiis
!!!  type(DFT_local_fields),intent(inout):: denspot
!!!  type(DFT_wavefunction),intent(inout):: tmb
!!!  
!!!  ! Local variables
!!!  type(orbitals_data):: orbs_tmp
!!!  type(locreg_descriptors):: glr_tmp
!!!  type(local_zone_descriptors):: lzd_tmp
!!!  real(8),dimension(:),pointer:: lphilarge, lhphilarge, lhphilargeold, lphilargeold, lhphi, lhphiold, lphiold
!!!  real(8),dimension(:,:),allocatable:: locregCenter
!!!  integer,dimension(:),allocatable:: inwhichlocreg_reference, onwhichatom_reference
!!!  integer:: istat, iall, iorb, ilr
!!!  character(len=*),parameter:: subname='enlarge_locreg'
!!!  
!!!  
!!!  allocate(locregCenter(3,lzd%nlr), stat=istat)
!!!  call memocc(istat, locregCenter, 'locregCenter', subname)
!!!  
!!!  call nullify_orbitals_data(orbs_tmp)
!!!  call nullify_locreg_descriptors(glr_tmp)
!!!  call nullify_local_zone_descriptors(lzd_tmp)
!!!  call copy_orbitals_data(tmb%orbs, orbs_tmp, subname)
!!!  call copy_locreg_descriptors(tmb%lzd%glr, glr_tmp, subname)
!!!  call copy_local_zone_descriptors(tmb%lzd, lzd_tmp, subname)
!!!  
!!!  do iorb=1,tmb%orbs%norb
!!!      ilr=tmb%orbs%inwhichlocreg(iorb)
!!!      locregCenter(:,ilr)=lzd%llr(ilr)%locregCenter
!!!  end do
!!!  
!!!  call destroy_new_locregs(iproc, nproc, tmb)
!!!  call update_locreg(iproc, nproc, lzd_tmp%nlr, locrad, orbs_tmp%inwhichlocreg, locregCenter, lzd_tmp%glr, &
!!!       withder, denspot%dpcom%nscatterarr, hx, hy, hz, &
!!!       orbs_tmp, tmb%lzd, tmb%orbs, tmb%op, tmb%comon, &
!!!       tmb%comgp, tmb%comsr, tmb%mad, tmb%collcom)
!!!  
!!!  allocate(lphilarge(tmb%orbs%npsidim_orbs), stat=istat)
!!!  call memocc(istat, lphilarge, 'lphilarge', subname)
!!!  call small_to_large_locreg(iproc, nproc, lzd_tmp, tmb%lzd, orbs_tmp, tmb%orbs, tmb%psi, lphilarge)
!!!  iall=-product(shape(tmb%psi))*kind(tmb%psi)
!!!  deallocate(tmb%psi, stat=istat)
!!!  call memocc(istat, iall, 'tmb%psi', subname)
!!!  allocate(tmb%psi(tmb%orbs%npsidim_orbs), stat=istat)
!!!  call memocc(istat, tmb%psi, 'tmb%psi', subname)
!!!  call dcopy(tmb%orbs%npsidim_orbs, lphilarge(1), 1, tmb%psi(1), 1)
!!!  !nphi=tmb%orbs%npsidim_orbs
!!!  
!!!  call update_ldiis_arrays(tmb, subname, ldiis)
!!!  call vcopy(tmb%orbs%norb, orbs_tmp%onwhichatom(1), 1, tmb%orbs%onwhichatom(1), 1)
!!!  iall=-product(shape(lphilarge))*kind(lphilarge)
!!!  deallocate(lphilarge, stat=istat)
!!!  call memocc(istat, iall, 'lphilarge', subname)
!!!  
!!!  iall=-product(shape(locregCenter))*kind(locregCenter)
!!!  deallocate(locregCenter, stat=istat)
!!!  call memocc(istat, iall, 'locregCenter', subname)
!!!  
!!!  call deallocate_orbitals_data(orbs_tmp, subname)
!!!  call deallocate_locreg_descriptors(glr_tmp, subname)
!!!  call deallocate_local_zone_descriptors(lzd_tmp, subname)
!!!
!!!end subroutine enlarge_locreg



!!subroutine initCompressedMatmul(iproc, nproc, norb, mad)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, norb
!!  type(matrixDescriptors),intent(inout):: mad
!!  
!!  ! Local variables
!!  integer:: iorb, jorb, ii, j, istat, iall, ij, iseg
!!  logical:: segment
!!  integer,dimension(:),allocatable:: row, column
!!  character(len=*),parameter:: subname='initCompressedMatmul'
!!  
!!  
!!  allocate(row(norb), stat=istat)
!!  call memocc(istat, row, 'row', subname)
!!  allocate(column(norb), stat=istat)
!!  call memocc(istat, column, 'column', subname)
!!  
!!  
!!  segment=.false.
!!  mad%nsegmatmul=0
!!  mad%nvctrmatmul=0
!!  do iorb=1,norb
!!      do jorb=1,norb
!!          ! Get an array of this line and column indicating whether
!!          ! there are nonzero numbers at these positions. Since the localization
!!          ! within the matrix is symmetric, we can use both time the same subroutine.
!!          call getRow(norb, mad, iorb, row) 
!!          call getRow(norb, mad, jorb, column) 
!!          !!if(iproc==0) write(*,'(a,i4,4x,100i4)') 'iorb, row', iorb, row
!!          !!if(iproc==0) write(*,'(a,i4,4x,100i4)') 'jorb, row', jorb, column
!!          ii=0
!!          do j=1,norb
!!              ii=ii+row(j)*column(j)
!!          end do
!!          if(ii>0) then
!!              ! This entry of the matrix will be different from zero.
!!              mad%nvctrmatmul=mad%nvctrmatmul+1
!!              if(.not. segment) then
!!                  ! This is the start of a new segment
!!                  segment=.true.
!!                  mad%nsegmatmul=mad%nsegmatmul+1
!!              end if
!!          else
!!              if(segment) then
!!                  ! We reached the end of a segment
!!                  segment=.false.
!!              end if
!!          end if
!!      end do
!!  end do
!!  
!!  allocate(mad%keygmatmul(2,mad%nsegmatmul), stat=istat)
!!  allocate(mad%keyvmatmul(mad%nsegmatmul), stat=istat)
!!  
!!  ! Now fill the descriptors.
!!  segment=.false.
!!  ij=0
!!  iseg=0
!!  do iorb=1,norb
!!      do jorb=1,norb
!!          ij=ij+1
!!          ! Get an array of this line and column indicating whether
!!          ! there are nonzero numbers at these positions. Since the localization
!!          ! within the matrix is symmetric, we can use both time the same subroutine.
!!          call getRow(norb, mad, iorb, row) 
!!          call getRow(norb, mad, jorb, column) 
!!          ii=0
!!          do j=1,norb
!!              ii=ii+row(j)*column(j)
!!          end do
!!          if(ii>0) then
!!              ! This entry of the matrix will be different from zero.
!!              if(.not. segment) then
!!                  ! This is the start of a new segment
!!                  segment=.true.
!!                  iseg=iseg+1
!!                  mad%keygmatmul(1,iseg)=ij
!!              end if
!!              mad%keyvmatmul(iseg)=mad%keyvmatmul(iseg)+1
!!          else
!!              if(segment) then
!!                  ! We reached the end of a segment
!!                  segment=.false.
!!                  mad%keygmatmul(2,iseg)=ij-1
!!              end if
!!          end if
!!      end do
!!  end do
!!
!!  ! Close the last segment if required.
!!  if(segment) then
!!      mad%keygmatmul(2,iseg)=ij
!!  end if
!!  
!!  
!!  iall=-product(shape(row))*kind(row)
!!  deallocate(row, stat=istat)
!!  call memocc(istat, iall, 'row', subname)
!!  iall=-product(shape(column))*kind(column)
!!  deallocate(column, stat=istat)
!!  call memocc(istat, iall, 'column', subname)
!!
!!end subroutine initCompressedMatmul



!!!subroutine initCompressedMatmul2(norb, nseg, keyg, nsegmatmul, keygmatmul, keyvmatmul)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!
!!!  ! Calling arguments
!!!  integer,intent(in):: norb, nseg
!!!  integer,dimension(2,nseg),intent(in):: keyg
!!!  integer,intent(out):: nsegmatmul
!!!  integer,dimension(:,:),pointer,intent(out):: keygmatmul
!!!  integer,dimension(:),pointer,intent(out):: keyvmatmul
!!!
!!!  ! Local variables
!!!  integer:: iorb, jorb, ii, j, istat, iall, ij, iseg, i
!!!  logical:: segment
!!!  character(len=*),parameter:: subname='initCompressedMatmul2'
!!!  real(8),dimension(:),allocatable:: mat1, mat2, mat3
!!!
!!!
!!!
!!!  allocate(mat1(norb**2), stat=istat)
!!!  call memocc(istat, mat1, 'mat1', subname)
!!!  allocate(mat2(norb**2), stat=istat)
!!!  call memocc(istat, mat2, 'mat2', subname)
!!!  allocate(mat3(norb**2), stat=istat)
!!!  call memocc(istat, mat2, 'mat2', subname)
!!!
!!!  mat1=0.d0
!!!  mat2=0.d0
!!!  do iseg=1,nseg
!!!      do i=keyg(1,iseg),keyg(2,iseg)
!!!          ! the localization region is "symmetric"
!!!          mat1(i)=1.d0
!!!          mat2(i)=1.d0
!!!      end do
!!!  end do
!!!
!!!  call dgemm('n', 'n', norb, norb, norb, 1.d0, mat1, norb, mat2, norb, 0.d0, mat3, norb)
!!!
!!!  segment=.false.
!!!  nsegmatmul=0
!!!  do iorb=1,norb**2
!!!      if(mat3(iorb)>0.d0) then
!!!          ! This entry of the matrix will be different from zero.
!!!          if(.not. segment) then
!!!              ! This is the start of a new segment
!!!              segment=.true.
!!!              nsegmatmul=nsegmatmul+1
!!!          end if
!!!      else
!!!          if(segment) then
!!!              ! We reached the end of a segment
!!!              segment=.false.
!!!          end if
!!!      end if
!!!  end do
!!!
!!!
!!!  allocate(keygmatmul(2,nsegmatmul), stat=istat)
!!!  call memocc(istat, keygmatmul, 'keygmatmul', subname)
!!!  allocate(keyvmatmul(nsegmatmul), stat=istat)
!!!  call memocc(istat, keyvmatmul, 'keyvmatmul', subname)
!!!  keyvmatmul=0
!!!  ! Now fill the descriptors.
!!!  segment=.false.
!!!  ij=0
!!!  iseg=0
!!!  do iorb=1,norb**2
!!!      ij=iorb
!!!      if(mat3(iorb)>0.d0) then
!!!          ! This entry of the matrix will be different from zero.
!!!          if(.not. segment) then
!!!              ! This is the start of a new segment
!!!              segment=.true.
!!!              iseg=iseg+1
!!!              keygmatmul(1,iseg)=ij
!!!          end if
!!!          keyvmatmul(iseg)=keyvmatmul(iseg)+1
!!!      else
!!!          if(segment) then
!!!              ! We reached the end of a segment
!!!              segment=.false.
!!!              keygmatmul(2,iseg)=ij-1
!!!          end if
!!!      end if
!!!  end do
!!!  ! Close the last segment if required.
!!!  if(segment) then
!!!      keygmatmul(2,iseg)=ij
!!!  end if
!!!
!!!
!!!iall=-product(shape(mat1))*kind(mat1)
!!!deallocate(mat1, stat=istat)
!!!call memocc(istat, iall, 'mat1', subname)
!!!iall=-product(shape(mat2))*kind(mat2)
!!!deallocate(mat2, stat=istat)
!!!call memocc(istat, iall, 'mat2', subname)
!!!iall=-product(shape(mat3))*kind(mat3)
!!!deallocate(mat3, stat=istat)
!!!call memocc(istat, iall, 'mat3', subname)
!!!
!!!
!!!end subroutine initCompressedMatmul2



!!subroutine initCompressedMatmul3(iproc, norb, mad)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, norb
!!  type(matrixDescriptors),intent(inout):: mad
!!
!!  ! Local variables
!!  integer:: iorb, jorb, ii, j, istat, iall, ij, iseg, i
!!  logical:: segment
!!  character(len=*),parameter:: subname='initCompressedMatmul3'
!!  real(8),dimension(:),allocatable:: mat1, mat2, mat3
!!
!!  call timing(iproc,'initMatmulComp','ON')
!!
!!  allocate(mat1(norb**2), stat=istat)
!!  call memocc(istat, mat1, 'mat1', subname)
!!  allocate(mat2(norb**2), stat=istat)
!!  call memocc(istat, mat2, 'mat2', subname)
!!  allocate(mat3(norb**2), stat=istat)
!!  call memocc(istat, mat2, 'mat2', subname)
!!
!!  call to_zero(norb**2, mat1(1))
!!  call to_zero(norb**2, mat2(1))
!!  do iseg=1,mad%nseg
!!      do i=mad%keyg(1,iseg),mad%keyg(2,iseg)
!!          ! the localization region is "symmetric"
!!          mat1(i)=1.d0
!!          mat2(i)=1.d0
!!      end do
!!  end do
!!
!!  call dgemm('n', 'n', norb, norb, norb, 1.d0, mat1, norb, mat2, norb, 0.d0, mat3, norb)
!!
!!  segment=.false.
!!  mad%nsegmatmul=0
!!  do iorb=1,norb**2
!!      if(mat3(iorb)>0.d0) then
!!          ! This entry of the matrix will be different from zero.
!!          if(.not. segment) then
!!              ! This is the start of a new segment
!!              segment=.true.
!!              mad%nsegmatmul=mad%nsegmatmul+1
!!          end if
!!      else
!!          if(segment) then
!!              ! We reached the end of a segment
!!              segment=.false.
!!          end if
!!      end if
!!  end do
!!
!!
!!  allocate(mad%keygmatmul(2,mad%nsegmatmul), stat=istat)
!!  call memocc(istat, mad%keygmatmul, 'mad%keygmatmul', subname)
!!  allocate(mad%keyvmatmul(mad%nsegmatmul), stat=istat)
!!  call memocc(istat, mad%keyvmatmul, 'mad%keyvmatmul', subname)
!!  mad%keyvmatmul=0
!!  ! Now fill the descriptors.
!!  segment=.false.
!!  ij=0
!!  iseg=0
!!  do iorb=1,norb**2
!!      ij=iorb
!!      if(mat3(iorb)>0.d0) then
!!          ! This entry of the matrix will be different from zero.
!!          if(.not. segment) then
!!              ! This is the start of a new segment
!!              segment=.true.
!!              iseg=iseg+1
!!              mad%keygmatmul(1,iseg)=ij
!!          end if
!!          mad%keyvmatmul(iseg)=mad%keyvmatmul(iseg)+1
!!      else
!!          if(segment) then
!!              ! We reached the end of a segment
!!              segment=.false.
!!              mad%keygmatmul(2,iseg)=ij-1
!!          end if
!!      end if
!!  end do
!!  ! Close the last segment if required.
!!  if(segment) then
!!      mad%keygmatmul(2,iseg)=ij
!!  end if
!!
!!
!!iall=-product(shape(mat1))*kind(mat1)
!!deallocate(mat1, stat=istat)
!!call memocc(istat, iall, 'mat1', subname)
!!iall=-product(shape(mat2))*kind(mat2)
!!deallocate(mat2, stat=istat)
!!call memocc(istat, iall, 'mat2', subname)
!!iall=-product(shape(mat3))*kind(mat3)
!!deallocate(mat3, stat=istat)
!!call memocc(istat, iall, 'mat3', subname)
!!
!!call timing(iproc,'initMatmulComp','OF')
!!
!!end subroutine initCompressedMatmul3



!!subroutine initCollectiveComms(iproc, nproc, lzd, input, orbs, collcomms)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(local_zone_descriptors),intent(in):: lzd
!!type(input_variables),intent(in):: input
!!type(orbitals_data),intent(inout):: orbs
!!type(collectiveComms),intent(out):: collcomms
!!
!!! Local variables
!!integer:: iorb, ilr, kproc, jproc, ii, ncount, iiorb, istat, gdim, ldim, ist
!!integer:: n1l, n2l, n3l, n1g, n2g, n3g, nshift1, nshift2, nshift3, ind, i, is, ie
!!integer:: transform_index, iseg, offset, iall
!!integer,dimension(:),allocatable:: work_int
!!character(len=*),parameter:: subname='initCollectiveComms'
!!integer:: ii1s, ii1e, ii5s, ii5e, i1, i5
!!logical:: stop1, stop5
!!
!!! Allocate all arrays
!!allocate(collComms%nvctr_par(orbs%norb,0:nproc-1), stat=istat)
!!call memocc(istat, collComms%nvctr_par, 'collComms%nvctr_par', subname)
!!
!!allocate(collComms%sendcnts(0:nproc-1), stat=istat)
!!call memocc(istat, collComms%sendcnts, 'collComms%sendcnts', subname)
!!
!!allocate(collComms%senddspls(0:nproc-1), stat=istat)
!!call memocc(istat, collComms%senddspls, 'collComms%senddspls', subname)
!!
!!allocate(collComms%recvcnts(0:nproc-1), stat=istat)
!!call memocc(istat, collComms%recvcnts, 'collComms%recvcnts', subname)
!!
!!allocate(collComms%recvdspls(0:nproc-1), stat=istat)
!!call memocc(istat, collComms%recvdspls, 'collComms%recvdspls', subname)
!!
!!
!!! Distribute the orbitals among the processes.
!!do iorb=1,orbs%norb
!!    ilr=orbs%inwhichlocreg(iorb)
!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    ! All processes get ii elements
!!    ii=ncount/nproc
!!    do jproc=0,nproc-1
!!        collComms%nvctr_par(iorb,jproc)=ii
!!    end do
!!    ! Process from 0 to kproc get one additional element
!!    kproc=mod(ncount,nproc)-1
!!    do jproc=0,kproc
!!        collComms%nvctr_par(iorb,jproc)=collComms%nvctr_par(iorb,jproc)+1
!!    end do
!!    !write(*,'(a,3i6,i12)') 'iorb, iproc, ncount, collComms%nvctr_par(iorb,iproc)', iorb, iproc, ncount, collComms%nvctr_par(iorb,iproc)
!!end do
!!
!!! Determine the amount of data that has has to be sent to each process
!!collComms%sendcnts=0
!!do jproc=0,nproc-1
!!    do iorb=1,orbs%norbp
!!        iiorb=orbs%isorb+iorb
!!        collComms%sendcnts(jproc) = collComms%sendcnts(jproc) + collComms%nvctr_par(iiorb,jproc)
!!    end do
!!    !write(*,'(a,2i6,i12)') 'jproc, iproc, collComms%sendcnts(jproc)', jproc, iproc, collComms%sendcnts(jproc)
!!end do
!!
!!! Determine the displacements for the send operation
!!collComms%senddspls(0)=0
!!do jproc=1,nproc-1
!!    collComms%senddspls(jproc) = collComms%senddspls(jproc-1) + collComms%sendcnts(jproc-1)
!!    !write(*,'(a,2i6,i12)') 'jproc, iproc, collComms%senddspls(jproc)', jproc, iproc, collComms%senddspls(jproc)
!!end do
!!
!!! Determine the amount of data that each process receives
!!collComms%recvcnts=0
!!do jproc=0,nproc-1
!!    do iorb=1,orbs%norb_par(jproc,0)
!!        iiorb=orbs%isorb_par(jproc)+iorb
!!        collComms%recvcnts(jproc) = collComms%recvcnts(jproc) + collComms%nvctr_par(iiorb,iproc)
!!    end do
!!    !write(*,'(a,2i6,i12)') 'jproc, iproc, collComms%recvcnts(jproc)', jproc, iproc, collComms%recvcnts(jproc)
!!end do
!!
!!! Determine the displacements for the receive operation
!!collComms%recvdspls(0)=0
!!do jproc=1,nproc-1
!!   collComms%recvdspls(jproc) = collComms%recvdspls(jproc-1) + collComms%recvcnts(jproc-1)
!!    !write(*,'(a,2i6,i12)') 'jproc, iproc, collComms%recvdspls(jproc)', jproc, iproc, collComms%recvdspls(jproc)
!!end do
!!
!!! Modify orbs%npsidim, if required
!!ii=0
!!do jproc=0,nproc-1
!!    ii=ii+collComms%recvcnts(jproc)
!!end do
!!!orbs%npsidim=max(orbs%npsidim,ii)
!!orbs%npsidim_orbs = max(orbs%npsidim_orbs,ii) 
!!orbs%npsidim_comp = max(orbs%npsidim_comp,ii)
!!
!!
!!ii1s=0
!!ii5s=0
!!ii1e=0
!!ii5e=0
!!
!!! Get the global indices of all elements
!!allocate(collComms%indexarray(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!call memocc(istat, collComms%indexarray, 'collComms%indexarray', subname)
!!ist=1
!!ind=1
!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    ilr=orbs%inwhichlocreg(iiorb)
!!    ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    !!!gdim=lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
!!    !!!call index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, &
!!    !!!     lzd%glr, lzd%llr(ilr), collComms%indexarray(ist))
!!    n1l=lzd%llr(ilr)%d%n1
!!    n2l=lzd%llr(ilr)%d%n2
!!    n3l=lzd%llr(ilr)%d%n3
!!    n1g=lzd%glr%d%n1
!!    n2g=lzd%glr%d%n2
!!    n3g=lzd%glr%d%n3
!!    !write(*,'(a,i8,6i9)') 'ilr, n1l, n2l, n3l, n1g, n2g, n3g', ilr, n1l, n2l, n3l, n1g, n2g, n3g
!!    nshift1=lzd%llr(ilr)%ns1-lzd%glr%ns1
!!    nshift2=lzd%llr(ilr)%ns2-lzd%glr%ns2
!!    nshift3=lzd%llr(ilr)%ns3-lzd%glr%ns3
!!
!!    if(iiorb==1) then
!!        ii1s=ind
!!    else if(iiorb==5) then
!!        ii5s=ind
!!    end if
!!    do iseg=1,lzd%llr(ilr)%wfd%nseg_c
!!        is=lzd%llr(ilr)%wfd%keygloc(1,iseg)
!!        ie=lzd%llr(ilr)%wfd%keygloc(2,iseg)
!!        !write(800+iiorb,'(a,i9,3i12,6i7)') 'ilr, iseg, is, ie, n1l, n2l, n3l, nshift1, nshift2, nshift3', &
!!        !      ilr, iseg, is, ie, n1l, n2l, n3l, nshift1, nshift2, nshift3
!!        do i=is,ie
!!            collComms%indexarray(ind)=transform_index(i, n1l, n2l, n3l, n1g, n2g, n3g, nshift1, nshift2, nshift3)
!!            !!!! DEBUG !!
!!            !!collComms%indexarray(ind)=iiorb
!!            !!!! DEBUG !!
!!            !!write(900+iiorb,'(a,i9,3i12,6i7,i10)') 'ilr, iseg, is, ie, n1l, n2l, n3l, nshift1, &
!!            !!    &nshift2, nshift3, collComms%indexarray(ind)', &
!!            !!    ilr, iseg, is, ie, n1l, n2l, n3l, nshift1, nshift2, nshift3, collComms%indexarray(ind)
!!            ind=ind+1
!!        end do
!!    end do
!!    !if(iiorb==1) then
!!    !    ii1e=ind-1
!!    !else if(iiorb==5) then
!!    !    ii5e=ind-1
!!    !end if
!!
!!
!!    offset=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
!!    do iseg=1,lzd%llr(ilr)%wfd%nseg_f
!!        is=lzd%llr(ilr)%wfd%keygloc(1,iseg+lzd%llr(ilr)%wfd%nseg_c)
!!        ie=lzd%llr(ilr)%wfd%keygloc(2,iseg+lzd%llr(ilr)%wfd%nseg_c)
!!        do i=is,ie
!!            ii=transform_index(i, n1l, n2l, n3l, n1g, n2g, n3g, nshift1, nshift2, nshift3)
!!
!!            collComms%indexarray(ind  ) = offset + 7*(ii-1)+1
!!            collComms%indexarray(ind+1) = offset + 7*(ii-1)+2
!!            collComms%indexarray(ind+2) = offset + 7*(ii-1)+3
!!            collComms%indexarray(ind+3) = offset + 7*(ii-1)+4
!!            collComms%indexarray(ind+4) = offset + 7*(ii-1)+5
!!            collComms%indexarray(ind+5) = offset + 7*(ii-1)+6
!!            collComms%indexarray(ind+6) = offset + 7*(ii-1)+7
!!            ind=ind+7
!!        end do
!!    end do
!!    if(iiorb==1) then
!!        ii1e=ind-1
!!    else if(iiorb==5) then
!!        ii5e=ind-1
!!    end if
!!
!!    !do istat=0,ldim-1
!!    !    write(200+iproc,*) ist+istat, collComms%indexarray(ist+istat)
!!    !end do
!!
!!    ist=ist+ldim
!!end do
!!
!!
!!!! ATTENTION: This will not work for nproc=1, so comment it.
!!!! As a consequence, the transposition will not work correctly.
!!
!!!!! Transpose the index array
!!!!allocate(work_int(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, work_int, 'work_int', subname)
!!!!call transpose_linear_int(iproc, 0, nproc-1, orbs, collComms, collComms%indexarray, mpi_comm_world, work_int)
!!!!iall=-product(shape(work_int))*kind(work_int)
!!!!deallocate(work_int, stat=istat)
!!!!call memocc(istat, iall, 'work_int', subname)
!!
!!
!!
!!end subroutine initCollectiveComms



! This subroutine is VERY similar to compressMatrix2...
subroutine getCommunArraysMatrixCompression(iproc, nproc, orbs, mad, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(matrixDescriptors),intent(in) :: mad
  integer,dimension(0:nproc-1),intent(out) :: sendcounts, displs
  
  ! Local variables
  integer :: iseg, jj, jorb, jjorb, jjproc, jjprocold, ncount
  
  sendcounts=0
  displs=0
  
  jj=0
  ncount=0
  jjprocold=0
  displs(0)=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          ncount=ncount+1
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(jjproc>jjprocold) then
              ! This part of the matrix is calculated by a new MPI process.
              sendcounts(jjproc-1)=ncount-1
              displs(jjproc)=displs(jjproc-1)+sendcounts(jjproc-1)
              ncount=1
              jjprocold=jjproc
          end if
      end do
  end do
  !sendcounts(nproc-1)=ncount
  sendcounts(jjproc)=ncount !last process
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if

  if(sum(sendcounts)/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix2: sum(sendcounts)/=mad%nvctr',sum(sendcounts),mad%nvctr
      stop
  end if

end subroutine getCommunArraysMatrixCompression
  
!!!subroutine initialize_comms_sumrho(iproc,nproc,nscatterarr,lzd,orbs,comsr)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  
!!!  ! Calling arguments
!!!  integer,intent(in) :: iproc,nproc
!!!  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!  type(local_zone_descriptors),intent(in) :: lzd
!!!  type(orbitals_data),intent(in) :: orbs
!!!  type(p2pComms),intent(out) :: comsr
!!!  
!!!  ! Local variables
!!!  integer :: istat,jproc,is,ie,ioverlap,i3s,i3e,ilr,iorb,is3ovrlp,n3ovrlp,iiproc,isend
!!!  integer :: jlr, jorb, istr, tag
!!!  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,p2p_tag,istsource,ncount
!!!  character(len=*),parameter :: subname='initialize_comms_sumrho'
!!!
!!!
!!!  call timing(iproc,'init_commSumro','ON')
!!!  
!!!  ! Buffer sizes 
!!!  call ext_buffers(lzd%Glr%geocode /= 'F',nbl1,nbr1)
!!!  call ext_buffers(lzd%Glr%geocode == 'P',nbl2,nbr2)
!!!  call ext_buffers(lzd%Glr%geocode /= 'F',nbl3,nbr3)
!!!  
!!!  ! First count the number of overlapping orbitals for each slice.
!!!  allocate(comsr%noverlaps(0:nproc-1),stat=istat)
!!!  call memocc(istat,comsr%noverlaps,'comsr%noverlaps',subname)
!!!  isend=0
!!!  do jproc=0,nproc-1
!!!      is=nscatterarr(jproc,3) 
!!!      ie=is+nscatterarr(jproc,1)-1
!!!      ioverlap=0
!!!      do iorb=1,orbs%norb
!!!          ilr=orbs%inWhichLocreg(iorb)
!!!          iiproc=orbs%onwhichmpi(iorb)
!!!          i3s=lzd%Llr(ilr)%nsi3 
!!!          i3e=i3s+lzd%Llr(ilr)%d%n3i-1
!!!          if(i3s<=ie .and. i3e>=is) then
!!!              ioverlap=ioverlap+1        
!!!              if(iproc==iiproc) isend=isend+1
!!!          end if
!!!          !For periodicity
!!!          if(i3e > Lzd%Glr%nsi3 + Lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F') then
!!!            i3s = Lzd%Glr%nsi3
!!!            i3e = mod(i3e-1,Lzd%Glr%d%n3i) + 1 + Lzd%Glr%nsi3
!!!            if(i3s<=ie .and. i3e>=is) then
!!!                ioverlap=ioverlap+1
!!!                if(iproc==iiproc) isend=isend+1
!!!            end if
!!!          end if
!!!      end do
!!!      comsr%noverlaps(jproc)=ioverlap
!!!  end do
!!!  
!!!  ! Do the initialization concerning the calculation of the charge density.
!!!  allocate(comsr%overlaps(comsr%noverlaps(iproc)),stat=istat)
!!!  call memocc(istat,comsr%overlaps,'comsr%overlaps',subname)
!!!  
!!!  allocate(comsr%comarr(6,maxval(comsr%noverlaps),0:nproc-1),stat=istat)
!!!  call memocc(istat,comsr%comarr,'comsr%comarr',subname)
!!!  allocate(comsr%ise3(comsr%noverlaps(iproc),2), stat=istat)
!!!  call memocc(istat, comsr%ise3, 'comsr%ise3', subname)
!!!  allocate(comsr%requests(max(comsr%noverlaps(iproc),isend),2),stat=istat)
!!!  call memocc(istat,comsr%requests,'comsr%requests',subname)
!!!  
!!!  comsr%nrecvBuf=0
!!!  do jproc=0,nproc-1
!!!     is=nscatterarr(jproc,3)
!!!     ie=is+nscatterarr(jproc,1)-1
!!!     ioverlap=0
!!!     istr=1
!!!     do iorb=1,orbs%norb
!!!        ilr=orbs%inWhichLocreg(iorb)
!!!        i3s=lzd%Llr(ilr)%nsi3
!!!        i3e=i3s+lzd%Llr(ilr)%d%n3i-1
!!!        if(i3s<=ie .and. i3e>=is) then
!!!           ioverlap=ioverlap+1
!!!           !tag=tag+1
!!!           tag=p2p_tag(jproc)
!!!           is3ovrlp=max(is,i3s) !start of overlapping zone in z direction
!!!           n3ovrlp=min(ie,i3e)-max(is,i3s)+1  !extent of overlapping zone in z direction
!!!           is3ovrlp=is3ovrlp-lzd%Llr(ilr)%nsi3+1
!!!           if(jproc == iproc) then
!!!              comsr%ise3(ioverlap,1) = max(is,i3s) 
!!!              comsr%ise3(ioverlap,2) = min(ie,i3e)
!!!           end if
!!!           istsource=1
!!!           do jorb=orbs%isorb_par(orbs%onwhichmpi(iorb))+1,iorb-1
!!!               jlr=orbs%inwhichlocreg(jorb)
!!!               istsource = istsource + lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n3i
!!!           end do
!!!           jlr=orbs%inwhichlocreg(iorb)
!!!           istsource = istsource + lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*(is3ovrlp-1)
!!!           ncount=lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*n3ovrlp
!!!           call setCommsParameters(orbs%onwhichmpi(iorb), jproc, istsource, istr, ncount, tag, comsr%comarr(1,ioverlap,jproc))
!!!           if(iproc==jproc) then
!!!              comsr%nrecvBuf = comsr%nrecvBuf + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*n3ovrlp
!!!              comsr%overlaps(ioverlap)=iorb
!!!           end if
!!!           istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*n3ovrlp
!!!        end if
!!!        !For periodicity
!!!        if(i3e > Lzd%Glr%nsi3 + Lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F') then
!!!           i3s = Lzd%Glr%nsi3
!!!           i3e = mod(i3e-1,Lzd%Glr%d%n3i) + 1 + Lzd%Glr%nsi3
!!!           if(i3s<=ie .and. i3e>=is) then
!!!              ioverlap=ioverlap+1
!!!              !tag=tag+1
!!!              tag=p2p_tag(jproc)
!!!              is3ovrlp=max(is,i3s) !start of overlapping zone in z direction
!!!              n3ovrlp=min(ie,i3e)-max(is,i3s)+1  !extent of overlapping zone in z direction
!!!              is3ovrlp=is3ovrlp + lzd%Glr%d%n3i-lzd%Llr(ilr)%nsi3+1 
!!!              if(jproc == iproc) then
!!!                 comsr%ise3(ioverlap,1) = max(is,i3s) 
!!!                 comsr%ise3(ioverlap,2) = min(ie,i3e)
!!!              end if
!!!              istsource=1
!!!              do jorb=orbs%isorb_par(orbs%onwhichmpi(iorb))+1,iorb-1
!!!                  jlr=orbs%inwhichlocreg(jorb)
!!!                  istsource = istsource + lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n3i
!!!              end do
!!!              jlr=orbs%inwhichlocreg(iorb)
!!!              istsource = istsource + lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*(is3ovrlp-1)
!!!              ncount=lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*n3ovrlp
!!!              call setCommsParameters(orbs%onwhichmpi(iorb), jproc, istsource, istr, ncount, tag, comsr%comarr(1,ioverlap,jproc))
!!!              if(iproc==jproc) then
!!!                 comsr%nrecvBuf = comsr%nrecvBuf + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*n3ovrlp
!!!                 comsr%overlaps(ioverlap)=iorb
!!!              end if
!!!              istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*n3ovrlp
!!!           end if
!!!           !For periodicity
!!!           if(i3e > Lzd%Glr%nsi3 + Lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F') then
!!!              i3s = Lzd%Glr%nsi3
!!!              i3e = mod(i3e-1,Lzd%Glr%d%n3i) + 1 + Lzd%Glr%nsi3
!!!              if(i3s<=ie .and. i3e>=is) then
!!!                 ioverlap=ioverlap+1
!!!                 !tag=tag+1
!!!                 tag=p2p_tag(jproc)
!!!                 is3ovrlp=max(is,i3s) !start of overlapping zone in z direction
!!!                 n3ovrlp=min(ie,i3e)-max(is,i3s)+1  !extent of overlapping zone in z direction
!!!                 is3ovrlp=is3ovrlp + lzd%Glr%d%n3i-lzd%Llr(ilr)%nsi3+1 
!!!                 if(jproc == iproc) then
!!!                    comsr%ise3(ioverlap,1) = max(is,i3s) 
!!!                    comsr%ise3(ioverlap,2) = min(ie,i3e)
!!!                 end if
!!!                 istsource=1
!!!                 do jorb=orbs%isorb_par(orbs%onwhichmpi(iorb))+1,iorb-1
!!!                     jlr=orbs%inwhichlocreg(jorb)
!!!                     istsource = istsource + lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n3i
!!!                 end do
!!!                 jlr=orbs%inwhichlocreg(iorb)
!!!                 istsource = istsource + lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*(is3ovrlp-1)
!!!                 ncount=lzd%llr(jlr)%d%n1i*lzd%llr(jlr)%d%n2i*n3ovrlp
!!!                 call setCommsParameters(orbs%onwhichmpi(iorb), jproc, istsource, istr, ncount, tag, comsr%comarr(1,ioverlap,jproc))
!!!                 if(iproc==jproc) then
!!!                    comsr%nrecvBuf = comsr%nrecvBuf + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*n3ovrlp
!!!                    comsr%overlaps(ioverlap)=iorb
!!!                 end if
!!!                 istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*n3ovrlp
!!!              end if
!!!           end if
!!!        end if
!!!     end do
!!!  end do
!!!  
!!!  ! To avoid allocations with size 0.
!!!  comsr%nrecvbuf=max(comsr%nrecvbuf,1)
!!!  
!!!  
!!!  
!!!  ! Calculate the dimension of the wave function for each process.
!!!  ! Do it for both the compressed ('npsidim') and for the uncompressed real space
!!!  ! ('npsidimr') case.
!!!  comsr%nsendBuf=0
!!!  do iorb=1,orbs%norbp
!!!      ilr=orbs%inWhichLocreg(orbs%isorb+iorb)
!!!      comsr%nsendBuf=comsr%nsendBuf+lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i*orbs%nspinor
!!!  end do
!!!
!!!  ! To indicate that no communication is going on.
!!!  comsr%communication_complete=.true.
!!!  comsr%messages_posted=.false.
!!!
!!!
!!!  call timing(iproc,'init_commSumro','OF')
!!!  
!!!end subroutine initialize_comms_sumrho
