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



