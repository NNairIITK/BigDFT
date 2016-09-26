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


subroutine dgemm_compressed_parallel(iproc, nproc, norb, nsegline, nseglinemax, keygline, &
           nsegmatmul, keygmatmul, norb_par, isorb_par, norbp, a, b, c)
!! ATTENTION: A MUST BE SYMMETRIC
use module_base
use module_types
implicit none

! Calling arguments
integer, intent(in) :: iproc, nproc, norb, norbp, nseglinemax, nsegmatmul
integer, dimension(2,nsegmatmul), intent(in) :: keygmatmul
integer, dimension(norb) :: nsegline
integer, dimension(2,nseglinemax,norb) :: keygline
integer, dimension(0:nproc-1), intent(in) :: norb_par, isorb_par
real(kind=8), dimension(norb,norb), intent(in) :: a, b
real(kind=8), dimension(norb,norb), intent(out) :: c

! Local variables
integer :: iseg, i, irow, icolumn, ii
integer :: ierr, istart, iend, iiseg, jjseg, ncount, jproc, istat, iall, iicolumn
real(kind=8) :: tt, ddot
logical :: iistop, jjstop
integer, dimension(:), allocatable :: sendcounts, displs
real(kind=8), dimension(:,:), allocatable :: c_loc
character(len=*),parameter :: subname='dgemm_compressed_parallel'


allocate(c_loc(norb,norbp), stat=istat)
call memocc(istat, c_loc, 'c_loc', subname)

!c=0.d0
c_loc=0.d0
if(norbp>0) call to_zero(norb*norbp, c_loc(1,1))
ii=0
do iseg=1,nsegmatmul
    do i=keygmatmul(1,iseg),keygmatmul(2,iseg)
        ii=ii+1
        ! Get the row and column index
        !irow=(i-1)/norb+1
        !icolumn=i-(irow-1)*norb
        icolumn=(i-1)/norb+1
        irow=i-(icolumn-1)*norb
        !if(irow>isorb_par(iproc) .and. irow<=isorb_par(min(iproc+1,nproc-1))) then
        if((icolumn>isorb_par(iproc) .and. icolumn<=isorb_par(min(iproc+1,nproc-1)))&
             .or. (iproc==nproc-1 .and. icolumn>isorb_par(iproc))) then
            !iirow=irow-isorb_par(iproc)
            iicolumn=icolumn-isorb_par(iproc)
            ! This process handles this entry of the matrix
            !c(irow,icolumn)=ddot(norb, a(1,irow), 1, b(1,icolumn), 1)
            iiseg=1
            jjseg=1
            iistop=.false.
            jjstop=.false.
            !write(*,'(a,3(i0,a))') 'process ',iproc,' calculates entry (',irow,',',iicolumn,')'
            do
                istart=max(keygline(1,iiseg,irow),keygline(1,jjseg,icolumn))
                iend=min(keygline(2,iiseg,irow),keygline(2,jjseg,icolumn))
                ncount=iend-istart+1

                if(ncount>0) then
                    tt=ddot(ncount, a(istart,irow), 1, b(istart,icolumn), 1)
                    !tt=ddot(ncount, a(istart,icolumn), 1, b(istart,irow), 1)
                else
                    tt=0.d0
                end if
                !c(irow,icolumn) = c(irow,icolumn) + tt
                !c_loc(icolumn,iirow) = c_loc(icolumn,iirow) + tt
                c_loc(irow,iicolumn) = c_loc(irow,iicolumn) + tt
                if(iiseg==nsegline(irow)) iistop=.true.
                if(jjseg==nsegline(icolumn)) jjstop=.true.
                if(iistop .and. jjstop) exit
                if((keygline(1,iiseg,irow)<=keygline(1,jjseg,icolumn) .or. jjstop) .and. .not.iistop) then
                    iiseg=iiseg+1
                else
                    jjseg=jjseg+1
                end if
            end do
            !write(*,'(5(a,i0),a,es15.6)') 'process ',iproc,': c_loc(',irow,',',iicolumn,')=c(',irow,',',icolumn,')=',c_loc(irow,iicolumn)
        end if
    end do
end do
!write(*,*) 'ii, norb**2', ii, norb**2

! Communicate the matrix.
allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)

displs(0)=0
do jproc=0,nproc-1
    sendcounts(jproc)=norb*norb_par(jproc)
    if(jproc>0) displs(jproc)=displs(jproc-1)+sendcounts(jproc-1)
end do
if (nproc > 1) then
   call mpi_allgatherv(c_loc(1,1), sendcounts(iproc), mpi_double_precision, c(1,1), sendcounts, displs, &
        mpi_double_precision, mpi_comm_world, ierr)
else
   call vcopy(sendcounts(iproc),c_loc(1,1),1,c(1,1),1)
end if

iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)
iall=-product(shape(c_loc))*kind(c_loc)
deallocate(c_loc, stat=istat)
call memocc(istat, iall, 'c_loc', subname)

!!do icolumn=1,norb
!!    do irow=1,norb
!!        if(iproc==0) write(201,*) icolumn, irow, c(irow,icolumn)
!!    end do
!!end do


end subroutine dgemm_compressed_parallel




subroutine compressMatrixPerProcess(iproc, nproc, orbs, mad, mat, size_lmat, lmat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc, nproc, size_lmat
  type(orbitals_data), intent(in) :: orbs
  type(matrixDescriptors), intent(in) :: mad
  real(kind=8), dimension(orbs%norb**2), intent(in) :: mat
  real(kind=8), dimension(size_lmat), intent(out) :: lmat
  
  ! Local variables
  integer :: iseg, jj, jorb, jjorb, jjproc
  
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(iproc==jjproc) then
              jj=jj+1
              lmat(jj)=mat(jorb)
          end if
      end do
  end do
  
end subroutine compressMatrixPerProcess



subroutine dgemm_compressed2(iproc, nproc, norb, nsegline, nseglinemax, keygline, nsegmatmul, keygmatmul, a, b, c)
!! ATTENTION: A MUST BE SYMMETRIC
use module_base
use module_types
implicit none

! Calling arguments
integer, intent(in) :: iproc, nproc, norb, nseglinemax, nsegmatmul
integer, dimension(2,nsegmatmul), intent(in) :: keygmatmul
integer, dimension(norb) :: nsegline
!integer, dimension(2,maxval(nsegline),norb) :: keygline
integer, dimension(2,nseglinemax,norb) :: keygline
real(kind=8), dimension(norb,norb), intent(in) :: a, b
real(kind=8), dimension(norb,norb), intent(out) :: c

! Local variables
integer :: iseg, i, irow, icolumn,  ii
integer :: istart, iend, iiseg, jjseg, ncount
real(kind=8) :: tt, ddot
logical :: iistop, jjstop


!!c=0.d0
call to_zero(norb**2, c(1,1))
ii=0
do iseg=1,nsegmatmul
    do i=keygmatmul(1,iseg),keygmatmul(2,iseg)
        ii=ii+1
        ! Get the row and column index
        irow=(i-1)/norb+1
        icolumn=i-(irow-1)*norb
        !c(irow,icolumn)=ddot(norb, a(1,irow), 1, b(1,icolumn), 1)
        iiseg=1
        jjseg=1
        iistop=.false.
        jjstop=.false.
        do
            istart=max(keygline(1,iiseg,irow),keygline(1,jjseg,icolumn))
            iend=min(keygline(2,iiseg,irow),keygline(2,jjseg,icolumn))
            ncount=iend-istart+1

            if(ncount>0) then
                tt=ddot(ncount, a(istart,irow), 1, b(istart,icolumn), 1)
            else
                tt=0.d0
            end if
            c(irow,icolumn) = c(irow,icolumn) + tt
            if(iiseg==nsegline(irow)) iistop=.true.
            if(jjseg==nsegline(icolumn)) jjstop=.true.
            if(iistop .and. jjstop) exit
            if((keygline(1,iiseg,irow)<=keygline(1,jjseg,icolumn) .or. jjstop) .and. .not.iistop) then
                iiseg=iiseg+1
            else
                jjseg=jjseg+1
            end if
        end do
        !if(iproc==0) write(*,'(3(a,i0),a,es15.6)') 'process ',iproc,': c(',irow,',',icolumn,')=',c(irow,icolumn)
    end do
end do
!write(*,*) 'ii, norb**2', ii, norb**2
!!do icolumn=1,norb
!!    do irow=1,norb
!!        if(iproc==0) write(200,*) icolumn, irow, c(irow,icolumn)
!!    end do
!!end do

end subroutine dgemm_compressed2






subroutine compressMatrix2(iproc, nproc, orbs, mad, mat, lmat, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc, nproc
  type(orbitals_data), intent(in) :: orbs
  type(matrixDescriptors), intent(in) :: mad
  real(kind=8), dimension(orbs%norb**2), intent(in) :: mat
  real(kind=8), dimension(mad%nvctr), intent(out) :: lmat
  integer, dimension(0:nproc-1), intent(out) :: sendcounts, displs
  
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
          lmat(jj)=mat(jorb)
          
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

  !if(iproc==0) then
  !    do jjproc=0,nproc-1
  !        write(*,'(a,3i8)') 'jjproc, displs(jjproc), sendcounts(jjproc)', jjproc, displs(jjproc), sendcounts(jjproc)
  !    end do
  !end if
  
end subroutine compressMatrix2

