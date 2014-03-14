!> @file
!!  Linear version: Handle Sparse Matrices
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 

subroutine initCommsOrtho(iproc, nproc, lzd, orbs, locregShape, noverlaps, overlaps) 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initCommsOrtho
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  character(len=1),intent(in) :: locregShape
  integer,dimension(:),pointer,intent(out):: noverlaps
  integer,dimension(:,:),pointer,intent(out):: overlaps

  ! Local variables
  integer ::  istat
  character(len=*),parameter :: subname='initCommsOrtho'


  call timing(iproc,'init_commOrtho','ON')

  ! Allocate the arrays that count the number of overlaps per process (comon%noverlaps)
  ! and per orbital (noverlaps)
  allocate(noverlaps(orbs%norb), stat=istat)
  call memocc(istat, noverlaps, 'noverlaps',subname)

  ! Count how many overlaping regions each orbital / process has.
  if(locregShape=='c') then
     stop "ERROR: locregShape=='c' is deprecated!"
  else if(locregShape=='s') then
     call determine_overlap_from_descriptors(iproc, nproc, orbs, orbs, lzd, lzd, noverlaps, overlaps)
  end if

  call timing(iproc,'init_commOrtho','OF')

end subroutine initCommsOrtho


!> Count for each orbital and each process the number of overlapping orbitals.
subroutine determine_overlap_from_descriptors(iproc, nproc, orbs, orbsig, lzd, lzdig, op_noverlaps, op_overlaps)
  use module_base
  use module_types
  use module_interfaces, except_this_one => determine_overlap_from_descriptors
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs, orbsig
  type(local_zone_descriptors),intent(in) :: lzd, lzdig
  integer,dimension(orbs%norb),intent(out):: op_noverlaps
  integer,dimension(:,:),pointer,intent(out):: op_overlaps

  ! Local variables
  integer :: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold
  !integer :: is1, ie1, is2, ie2, is3, ie3
  !integer :: js1, je1, js2, je2, js3, je3
  integer :: iiorb, istat, iall, noverlaps, ierr, ii
  !integer :: noverlapsmaxp
  !logical :: ovrlpx, ovrlpy, ovrlpz
  logical :: isoverlap
  !integer :: i1, i2
  integer :: onseg
  logical,dimension(:,:),allocatable :: overlapMatrix
  integer,dimension(:),allocatable :: noverlapsarr, displs, recvcnts, comon_noverlaps
  integer,dimension(:,:),allocatable :: overlaps_op
  integer,dimension(:,:,:),allocatable :: overlaps_nseg
  !integer,dimension(:,:,:),allocatable :: iseglist, jseglist
  character(len=*),parameter :: subname='determine_overlap_from_descriptors'

  allocate(overlapMatrix(orbsig%norb,maxval(orbs%norb_par(:,0))), stat=istat)
  call memocc(istat, overlapMatrix, 'overlapMatrix', subname)
  allocate(noverlapsarr(orbs%norbp), stat=istat)
  call memocc(istat, noverlapsarr, 'noverlapsarr', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  allocate(recvcnts(0:nproc-1), stat=istat)
  call memocc(istat, recvcnts, 'recvcnts', subname)
  allocate(overlaps_nseg(orbsig%norb,orbs%norbp,2), stat=istat)
  call memocc(istat, overlaps_nseg, 'overlaps_nseg', subname)

  allocate(comon_noverlaps(0:nproc-1), stat=istat)
  call memocc(istat, comon_noverlaps, 'comon_noverlaps',subname)

  overlapMatrix=.false.
  overlaps_nseg = 0
  ioverlapMPI=0 ! counts the overlaps for the given MPI process.
  ilrold=-1
  do iorb=1,orbs%norbp
     ioverlaporb=0 ! counts the overlaps for the given orbital.
     iiorb=orbs%isorb+iorb
     ilr=orbs%inWhichLocreg(iiorb)
!     call get_start_and_end_indices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
     do jorb=1,orbsig%norb
        jlr=orbsig%inWhichLocreg(jorb)
        call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzdig%llr(jlr),isoverlap)
!        call get_start_and_end_indices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
        !write(*,*) 'second: iproc, ilr, jlr, isoverlap', iproc, ilr, jlr, isoverlap
        if(isoverlap) then
           ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
           ! Now explicitely check whether there is an overlap by using the descriptors.
           !!call overlapbox_from_descriptors(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           !!     lzd%llr(jlr)%d%n1, lzd%llr(jlr)%d%n2, lzd%llr(jlr)%d%n3, &
           !!     lzd%glr%d%n1, lzd%glr%d%n2, lzd%glr%d%n3, &
           !!     lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
           !!     lzd%llr(jlr)%ns1, lzd%llr(jlr)%ns2, lzd%llr(jlr)%ns3, &
           !!     lzd%glr%ns1, lzd%glr%ns2, lzd%glr%ns3, &
           !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c, &
           !!     lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, lzd%llr(jlr)%wfd%keygloc, lzd%llr(jlr)%wfd%keyvloc, &
           !!     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)
           call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_c,&
                lzd%llr(ilr)%wfd%keyglob, lzdig%llr(jlr)%wfd%keyglob, &
                isoverlap, onseg)
           if(isoverlap) then
              ! There is really an overlap
              overlapMatrix(jorb,iorb)=.true.
              ioverlaporb=ioverlaporb+1
              overlaps_nseg(ioverlaporb,iorb,1)=onseg
              if(ilr/=ilrold) then
                 ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                 ! would get the same orbitals again. Therefore the counter is not increased
                 ! in that case.
                 ioverlapMPI=ioverlapMPI+1
              end if
           else
              overlapMatrix(jorb,iorb)=.false.
           end if
        else
           overlapMatrix(jorb,iorb)=.false.
        end if
     end do
     noverlapsarr(iorb)=ioverlaporb
     ilrold=ilr
  end do
  noverlaps=ioverlapMPI

  !call mpi_allreduce(overlapMatrix, orbs%norb*maxval(orbs%norb_par(:,0))*nproc, mpi_sum bigdft_mpi%mpi_comm, ierr)

  ! Communicate op_noverlaps and comon_noverlaps
  if (nproc > 1) then
     call mpi_allgatherv(noverlapsarr, orbs%norbp, mpi_integer, op_noverlaps, orbs%norb_par, &
          orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(orbs%norb,noverlapsarr(1),1,op_noverlaps(1),1)
  end if
  do jproc=0,nproc-1
     recvcnts(jproc)=1
     displs(jproc)=jproc
  end do
  if (nproc > 1) then
     call mpi_allgatherv(noverlaps, 1, mpi_integer, comon_noverlaps, recvcnts, &
          displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
  else
     comon_noverlaps=noverlaps
  end if

  allocate(op_overlaps(maxval(op_noverlaps),orbs%norb), stat=istat)
  call memocc(istat, op_overlaps, 'op_overlaps', subname)
  allocate(overlaps_op(maxval(op_noverlaps),orbs%norbp), stat=istat)
  call memocc(istat, overlaps_op, 'overlaps_op', subname)
  if(orbs%norbp>0) call to_zero(maxval(op_noverlaps)*orbs%norbp,overlaps_op(1,1))

  ! Now we know how many overlaps have to be calculated, so determine which orbital overlaps
  ! with which one. This is essentially the same loop as above, but we use the array 'overlapMatrix'
  ! which indicates the overlaps.

  ! Initialize to some value which will never be used.
  op_overlaps=-1

  iiorb=0
  ilrold=-1
  do iorb=1,orbs%norbp
     ioverlaporb=0 ! counts the overlaps for the given orbital.
     iiorb=orbs%isorb+iorb
     ilr=orbs%inWhichLocreg(iiorb)
     do jorb=1,orbsig%norb
        jlr=orbsig%inWhichLocreg(jorb)
        if(overlapMatrix(jorb,iorb)) then
           ioverlaporb=ioverlaporb+1
           ! Determine the number of segments of the fine grid in the overlap
           call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_f,&
                lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1+lzdig%llr(jlr)%wfd%nseg_c), &
                isoverlap, overlaps_nseg(ioverlaporb,iorb,2))
           !op_overlaps(ioverlaporb,iiorb)=jorb
           overlaps_op(ioverlaporb,iorb)=jorb
        end if
     end do 
     ilrold=ilr
  end do


  displs(0)=0
  recvcnts(0)=comon_noverlaps(0)
  do jproc=1,nproc-1
      recvcnts(jproc)=comon_noverlaps(jproc)
      displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
  end do

  ii=maxval(op_noverlaps)
  displs(0)=0
  recvcnts(0)=ii*orbs%norb_par(0,0)
  do jproc=1,nproc-1
      recvcnts(jproc)=ii*orbs%norb_par(jproc,0)
      displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
  end do
  if (nproc > 1) then
     call mpi_allgatherv(overlaps_op, ii*orbs%norbp, mpi_integer, op_overlaps, recvcnts, &
          displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(ii*orbs%norbp,overlaps_op(1,1),1,op_overlaps(1,1),1)
  end if


  iall=-product(shape(overlapMatrix))*kind(overlapMatrix)
  deallocate(overlapMatrix, stat=istat)
  call memocc(istat, iall, 'overlapMatrix', subname)

  iall=-product(shape(noverlapsarr))*kind(noverlapsarr)
  deallocate(noverlapsarr, stat=istat)
  call memocc(istat, iall, 'noverlapsarr', subname)

  iall=-product(shape(overlaps_op))*kind(overlaps_op)
  deallocate(overlaps_op, stat=istat)
  call memocc(istat, iall, 'overlaps_op', subname)

  iall=-product(shape(displs))*kind(displs)
  deallocate(displs, stat=istat)
  call memocc(istat, iall, 'displs', subname)

  iall=-product(shape(recvcnts))*kind(recvcnts)
  deallocate(recvcnts, stat=istat)
  call memocc(istat, iall, 'recvcnts', subname)

  iall=-product(shape(overlaps_nseg))*kind(overlaps_nseg)
  deallocate(overlaps_nseg, stat=istat)
  call memocc(istat, iall, 'overlaps_nseg', subname)

  iall=-product(shape(comon_noverlaps))*kind(comon_noverlaps)
  deallocate(comon_noverlaps, stat=istat)
  call memocc(istat, iall, 'comon_noverlaps', subname)


end subroutine determine_overlap_from_descriptors



!> Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
function compressed_index(irow, jcol, norb, sparsemat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: irow, jcol, norb
  type(sparseMatrix),intent(in) :: sparsemat
  integer :: compressed_index

  ! Local variables
  integer :: ii, iseg

  ii=(jcol-1)*norb+irow

  iseg=sparsemat%istsegline(jcol)
  do
      if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
          ! The matrix element is in this segment
           compressed_index = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
          return
      end if
      iseg=iseg+1
      if (iseg>sparsemat%nseg) exit
      if (ii<sparsemat%keyg(1,iseg)) then
          compressed_index=0
          return
      end if
  end do

  ! Not found
  compressed_index=0

end function compressed_index


subroutine compress_matrix_for_allreduce(iproc,sparsemat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc
  type(sparseMatrix),intent(inout) :: sparsemat

  ! Local variables
  integer :: jj, irow, jcol, jjj, ierr

  call timing(iproc,'compress_uncom','ON')

  if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
     !$omp parallel do default(private) shared(sparsemat)
     do jj=1,sparsemat%nvctr
        irow = sparsemat%orb_from_index(1,jj)
        jcol = sparsemat%orb_from_index(2,jj)
        sparsemat%matrix_compr(jj)=sparsemat%matrix(irow,jcol)
     end do
     !$omp end parallel do
  else if (sparsemat%parallel_compression==1) then
     call to_zero(sparsemat%nvctr, sparsemat%matrix_compr(1))
     !$omp parallel do default(private) shared(sparsemat)
     do jj=1,sparsemat%nvctrp
        jjj=jj+sparsemat%isvctr
        irow = sparsemat%orb_from_index(1,jjj)
        jcol = sparsemat%orb_from_index(2,jjj)
        sparsemat%matrix_compr(jjj)=sparsemat%matrix(irow,jcol)
     end do
     !$omp end parallel do
     call mpiallred(sparsemat%matrix_compr(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
     sparsemat%matrix_comprp=f_malloc_ptr((sparsemat%nvctrp),id='sparsemat%matrix_comprp')
     !$omp parallel do default(private) shared(sparsemat)
     do jj=1,sparsemat%nvctrp
        jjj=jj+sparsemat%isvctr
        irow = sparsemat%orb_from_index(1,jjj)
        jcol = sparsemat%orb_from_index(2,jjj)
        sparsemat%matrix_comprp(jj)=sparsemat%matrix(irow,jcol)
     end do
     !$omp end parallel do
     call mpi_allgatherv(sparsemat%matrix_comprp, sparsemat%nvctrp, mpi_double_precision, sparsemat%matrix_compr, &
          sparsemat%nvctr_par(:), sparsemat%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     call f_free_ptr(sparsemat%matrix_comprp)
  end if

  call timing(iproc,'compress_uncom','OF')

end subroutine compress_matrix_for_allreduce



subroutine uncompressMatrix(iproc,sparsemat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc
  type(sparseMatrix), intent(inout) :: sparsemat
  
  ! Local variables
  integer :: ii, irow, jcol, iii, ierr

  call timing(iproc,'compress_uncom','ON')

  if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
     call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
     !$omp parallel do default(private) shared(sparsemat)
     do ii=1,sparsemat%nvctr
        irow = sparsemat%orb_from_index(1,ii)
        jcol = sparsemat%orb_from_index(2,ii)
        sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(ii)
     end do
     !$omp end parallel do
  else if (sparsemat%parallel_compression==1) then
     call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
     !$omp parallel do default(private) shared(sparsemat)
     do ii=1,sparsemat%nvctrp
        iii=ii+sparsemat%isvctr
        irow = sparsemat%orb_from_index(1,iii)
        jcol = sparsemat%orb_from_index(2,iii)
        sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(iii)
     end do
     !$omp end parallel do
     call mpiallred(sparsemat%matrix(1,1), sparsemat%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
     sparsemat%matrixp=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctrp/),id='sparsemat%matrixp')
     call to_zero(sparsemat%nfvctr*sparsemat%nfvctrp, sparsemat%matrixp(1,1))
     !$omp parallel do default(private) shared(sparsemat)
     do ii=1,sparsemat%nvctrp
        iii=ii+sparsemat%isvctr
        irow = sparsemat%orb_from_index(1,iii)
        jcol = sparsemat%orb_from_index(2,iii) - sparsemat%isfvctr
        sparsemat%matrixp(irow,jcol)=sparsemat%matrix_compr(iii)
     end do
     !$omp end parallel do
     call mpi_allgatherv(sparsemat%matrixp, sparsemat%nfvctr*sparsemat%nfvctrp, mpi_double_precision, sparsemat%matrix, &
          sparsemat%nfvctr*sparsemat%nfvctr_par(:), sparsemat%nfvctr*sparsemat%isfvctr_par, mpi_double_precision, &
          bigdft_mpi%mpi_comm, ierr)
     call f_free_ptr(sparsemat%matrixp)
  end if
  sparsemat%can_use_dense=.true.  

  call timing(iproc,'compress_uncom','OF')

end subroutine uncompressMatrix



subroutine check_matrix_compression(iproc,sparsemat)
  use module_base
  use module_types
  use module_interfaces
  use yaml_output
  use sparsematrix_base, only: sparseMatrix
  implicit none
  integer,intent(in) :: iproc
  type(sparseMatrix),intent(inout) :: sparsemat
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



integer function matrixindex_in_compressed(sparsemat, iorb, jorb)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none

  ! Calling arguments
  type(sparseMatrix),intent(in) :: sparsemat
  integer,intent(in) :: iorb, jorb

  ! Local variables

  if (sparsemat%store_index) then
      ! Take the value from the array
      matrixindex_in_compressed = sparsemat%matrixindex_in_compressed_arr(iorb,jorb)
  else
      ! Recalculate the value
      matrixindex_in_compressed = compressed_index_fn(iorb, jorb, sparsemat%nfvctr, sparsemat)
  end if

  contains
    ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
    integer function compressed_index_fn(irow, jcol, norb, sparsemat)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: irow, jcol, norb
      type(sparseMatrix),intent(in) :: sparsemat
    
      ! Local variables
      integer :: ii, iseg
    
      ii=(jcol-1)*norb+irow
    
      iseg=sparsemat%istsegline(jcol)
      do
          if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
              ! The matrix element is in sparsemat segment
               compressed_index_fn = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
              return
          end if
          iseg=iseg+1
          if (iseg>sparsemat%nseg) exit
          if (ii<sparsemat%keyg(1,iseg)) then
              compressed_index_fn=0
              return
          end if
      end do
    
      ! Not found
      compressed_index_fn=0
    
    end function compressed_index_fn
end function matrixindex_in_compressed






subroutine init_indices_in_compressed(store_index, norb, sparsemat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none

  ! Calling arguments
  logical,intent(in) :: store_index
  integer,intent(in) :: norb
  type(sparseMatrix),intent(inout) :: sparsemat

  ! Local variables
  integer :: iorb, jorb, compressed_index, istat
  character(len=*),parameter :: subname='init_indices_in_compressed'

  if (store_index) then
      ! store the indices of the matrices in the sparse format
      sparsemat%store_index=.true.

      ! initialize sparsemat%matrixindex_in_compressed
      !allocate(sparsemat%matrixindex_in_compressed_arr(norb,norb),stat=istat)
      !call memocc(istat, sparsemat%matrixindex_in_compressed_arr, 'sparsemat%matrixindex_in_compressed_arr', subname)
      sparsemat%matrixindex_in_compressed_arr=f_malloc_ptr((/norb,norb/),id='sparsemat%matrixindex_in_compressed_arr')

      do iorb=1,norb
         do jorb=1,norb
            sparsemat%matrixindex_in_compressed_arr(iorb,jorb)=compressed_index(iorb,jorb,norb,sparsemat)
         end do
      end do

  else
      ! otherwise alwyas calculate them on-the-fly
      sparsemat%store_index=.false.
      nullify(sparsemat%matrixindex_in_compressed_arr)
  end if

end subroutine init_indices_in_compressed



subroutine init_orbs_from_index(sparsemat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none

  ! Calling arguments
  type(sparseMatrix),intent(inout) :: sparsemat

  ! local variables
  integer :: ind, iseg, segn, iorb, jorb, istat
  character(len=*),parameter :: subname='init_orbs_from_index'

  !allocate(sparsemat%orb_from_index(2,sparsemat%nvctr),stat=istat)
  !call memocc(istat, sparsemat%orb_from_index, 'sparsemat%orb_from_index', subname)
  sparsemat%orb_from_index=f_malloc_ptr((/2,sparsemat%nvctr/),id='sparsemat%orb_from_index')

  ind = 0
  do iseg = 1, sparsemat%nseg
     do segn = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
        ind=ind+1
        iorb = (segn - 1) / sparsemat%nfvctr + 1
        jorb = segn - (iorb-1)*sparsemat%nfvctr
        sparsemat%orb_from_index(1,ind) = jorb
        sparsemat%orb_from_index(2,ind) = iorb
     end do
  end do

end subroutine init_orbs_from_index




subroutine check_kernel_cutoff(iproc, orbs, atoms, lzd)
  use module_base
  use module_types
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: atoms
  type(local_zone_descriptors),intent(inout) :: lzd

  ! Local variables
  integer :: iorb, ilr, iat, iatype
  real(kind=8) :: cutoff_sf, cutoff_kernel
  character(len=20) :: atomname
  logical :: write_data
  logical,dimension(atoms%astruct%ntypes) :: write_atomtype

  write_atomtype=.true.

  if (iproc==0) then
      call yaml_open_sequence('check of kernel cutoff radius')
  end if

  do iorb=1,orbs%norb
      ilr=orbs%inwhichlocreg(iorb)

      ! cutoff radius of the support function, including shamop region
      cutoff_sf=lzd%llr(ilr)%locrad+8.d0*lzd%hgrids(1)

      ! cutoff of the density kernel
      cutoff_kernel=lzd%llr(ilr)%locrad_kernel

      ! check whether the date for this atomtype has already shoudl been written
      iat=orbs%onwhichatom(iorb)
      iatype=atoms%astruct%iatype(iat)
      if (write_atomtype(iatype)) then
          if (iproc==0) then
              write_data=.true.
          else
              write_data=.false.
          end if
          write_atomtype(iatype)=.false.
      else
          write_data=.false.
      end if

      ! Adjust if necessary
      if (write_data) then
          call yaml_sequence(advance='no')
          call yaml_open_map(flow=.true.)
          atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
          call yaml_map('atom type',atomname)
      end if
      if (cutoff_sf>cutoff_kernel) then
          if (write_data) then
              call yaml_map('adjustment required',.true.)
              call yaml_map('new value',cutoff_sf,fmt='(f6.2)')
          end if
          lzd%llr(ilr)%locrad_kernel=cutoff_sf
      else
          if (write_data) then
              call yaml_map('adjustment required',.false.)
          end if
      end if
      if (write_data) then
          call yaml_close_map()
      end if
  end do

  if (iproc==0) then
      call yaml_close_sequence
  end if


end subroutine check_kernel_cutoff
