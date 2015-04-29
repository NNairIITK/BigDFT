!> @file
!! Intialization of the collective communications for the linear version
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_pulay_overlap(iproc, nproc, orbs1, orbs2, collcom1, collcom2, psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs1, orbs2
  type(comms_linear),intent(in) :: collcom1, collcom2
  real(kind=8),dimension(collcom1%ndimind_c),intent(in) :: psit_c1
  real(kind=8),dimension(collcom2%ndimind_c),intent(in) :: psit_c2
  real(kind=8),dimension(7*collcom1%ndimind_f),intent(in) :: psit_f1
  real(kind=8),dimension(7*collcom2%ndimind_f),intent(in) :: psit_f2
  real(kind=8),dimension(orbs1%norb,orbs2%norb),intent(out) :: ovrlp
  
  ! Local variables
  integer :: i0, j0, ipt, ii, iiorb, j, jj, jjorb, i, ierr  

  call timing(iproc,'ovrlptransComp','ON') !lr408t
  call f_zero(ovrlp)
  if(collcom1%nptsp_c/=collcom2%nptsp_c) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_c/=collcom2%nptsp_c'
      stop
  end if
  if(collcom1%nptsp_f/=collcom2%nptsp_f) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_f/=collcom2%nptsp_f'
      stop
  end if

  i0=0
  j0=0
  do ipt=1,collcom1%nptsp_c 
      ii=collcom1%norb_per_gridpoint_c(ipt)
      jj=collcom2%norb_per_gridpoint_c(ipt)
      do i=1,ii
          iiorb=collcom1%indexrecvorbital_c(i0+i)
          do j=1,jj
              jjorb=collcom2%indexrecvorbital_c(j0+j)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_c1(i0+i)*psit_c2(j0+j)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  i0=0
  j0=0
  do ipt=1,collcom1%nptsp_f 
      ii=collcom1%norb_per_gridpoint_f(ipt)
      jj=collcom2%norb_per_gridpoint_f(ipt)
      do i=1,ii
          iiorb=collcom1%indexrecvorbital_f(i0+i)
          do j=1,jj
              jjorb=collcom2%indexrecvorbital_f(j0+j)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-6)*psit_f2(7*(j0+j)-6)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-5)*psit_f2(7*(j0+j)-5)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-4)*psit_f2(7*(j0+j)-4)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-3)*psit_f2(7*(j0+j)-3)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-2)*psit_f2(7*(j0+j)-2)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-1)*psit_f2(7*(j0+j)-1)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-0)*psit_f2(7*(j0+j)-0)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc>1) then
      call mpiallred(ovrlp, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if
  call timing(iproc,'ovrlptransComm','OF') !lr408t
end subroutine calculate_pulay_overlap





subroutine check_grid_point_from_boxes(i1, i2, i3, lr, overlap_possible)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: i1, i2, i3
  type(locreg_descriptors),intent(in) :: lr  
  logical,intent(out) :: overlap_possible

  ! Local variables
  logical :: ovrlpx, ovrlpy, ovrlpz
  
  ovrlpx = (i1>=lr%ns1 .and. i1<=lr%ns1+lr%d%n1)
  ovrlpy = (i2>=lr%ns2 .and. i2<=lr%ns2+lr%d%n2)
  ovrlpz = (i3>=lr%ns3 .and. i3<=lr%ns3+lr%d%n3)
  if(ovrlpx .and. ovrlpy .and. ovrlpz) then
      overlap_possible=.true.
  else
      overlap_possible=.true.
  end if

end subroutine check_grid_point_from_boxes

!!subroutine get_reverse_indices(n, indices, reverse_indices)
!!  use module_base
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in) :: n
!!  integer,dimension(n),intent(in) :: indices
!!  integer,dimension(n),intent(out) :: reverse_indices
!!
!!  ! Local variables
!!  integer :: i, j, m, j0, j1, j2, j3
!!
!!  !$omp parallel default(private) &
!!  !$omp shared(n, m, indices, reverse_indices)
!!
!!  m=mod(n,4)
!!  if (m/=0) then
!!      do i=1,m
!!          j=indices(i)
!!          reverse_indices(j)=i
!!      end do
!!  end if
!!
!!  !$omp do
!!  do i=m+1,n,4
!!      j0=indices(i+0)
!!      reverse_indices(j0)=i+0
!!      j1=indices(i+1)
!!      reverse_indices(j1)=i+1
!!      j2=indices(i+2)
!!      reverse_indices(j2)=i+2
!!      j3=indices(i+3)
!!      reverse_indices(j3)=i+3
!!  end do
!!  !$omp end do
!!
!!  !$omp end parallel
!!
!!  !!do i=1,n
!!  !!    j=indices(i)
!!  !!    reverse_indices(j)=i
!!  !!end do
!!
!!end subroutine get_reverse_indices





subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, orbs, collcom, collcom_shamop, &
           collcom_sr, sparsemat)
  use module_base
  use module_types
  !use module_interfaces, except_this_one => init_matrixindex_in_compressed_fortransposed
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed!compressed_index
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(comms_linear),intent(in) :: collcom, collcom_shamop, collcom_sr
  type(sparse_matrix), intent(inout) :: sparsemat
  
  ! Local variables
  integer :: iorb, jorb, istat, imin, imax, nmiddle, imax_old, imin_old, iiorb, jjorb
  integer :: ii, imin_new, imax_new, i, nlen, j
  !integer :: kproc,jproc,jjorbold,jjorb,isend,irecv,ilr,ijorb,iiorb,ind,ierr, irow, irowold, iseg
  !integer :: compressed_index
!  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
  character(len=*),parameter :: subname='init_sparse_matrix'

  call f_routine(id='init_matrixindex_in_compressed_fortransposed')


  ! for the calculation of overlaps and the charge density
  !imin=minval(collcom%indexrecvorbital_c)
  !imin=min(imin,minval(collcom%indexrecvorbital_f))
  !imin=min(imin,minval(collcom_shamop%indexrecvorbital_c))
  !imin=min(imin,minval(collcom_shamop%indexrecvorbital_f))
  !imin=min(imin,minval(collcom_sr%indexrecvorbital_c))
  !imax=maxval(collcom%indexrecvorbital_c)
  !imax=max(imax,maxval(collcom%indexrecvorbital_f))
  !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_c))
  !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_f))
  !imax=max(imax,maxval(collcom_sr%indexrecvorbital_c))

  nmiddle = sparsemat%nfvctr/2 + 1

  imin_old = huge(1)
  imax_old = 0
  imin_new = huge(1)
  imax_new = 0
  do i=1,size(collcom%indexrecvorbital_c)
      ii = mod(collcom%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
      imin_old = min(imin_old,ii)
      imax_old = max(imax_old,ii)
      if (ii>nmiddle) then
          imin_new = min(imin_new,ii)
      else
          imax_new = max(imax_new,ii+sparsemat%nfvctr)
      end if
  end do
  do i=1,size(collcom%indexrecvorbital_f)
      ii = mod(collcom%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
      imin_old = min(imin_old,ii)
      imax_old = max(imax_old,ii)
      if (ii>nmiddle) then
          imin_new = min(imin_new,ii)
      else
          imax_new = max(imax_new,ii+sparsemat%nfvctr)
      end if
  end do
  do i=1,size(collcom_shamop%indexrecvorbital_c)
      ii = mod(collcom_shamop%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
      imin_old = min(imin_old,ii)
      imax_old = max(imax_old,ii)
      if (ii>nmiddle) then
          imin_new = min(imin_new,ii)
      else
          imax_new = max(imax_new,ii+sparsemat%nfvctr)
      end if
  end do
  do i=1,size(collcom_shamop%indexrecvorbital_f)
      ii = mod(collcom_shamop%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
      imin_old = min(imin_old,ii)
      imax_old = max(imax_old,ii)
      if (ii>nmiddle) then
          imin_new = min(imin_new,ii)
      else
          imax_new = max(imax_new,ii+sparsemat%nfvctr)
      end if
  end do
  do i=1,size(collcom_sr%indexrecvorbital_c)
      ii = mod(collcom_sr%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
      imin_old = min(imin_old,ii)
      imax_old = max(imax_old,ii)
      if (ii>nmiddle) then
          imin_new = min(imin_new,ii)
      else
          imax_new = max(imax_new,ii+sparsemat%nfvctr)
      end if
  end do


  write(*,*) 'iproc, imin_old, imax_old', iproc, imin_old, imax_old
  write(*,*) 'iproc, imin_new, imax_new', iproc, imin_new, imax_new

  !! values regardless of the spin
  !imin=mod(imin-1,sparsemat%nfvctr)+1
  !imax=mod(imax-1,sparsemat%nfvctr)+1


  ! Determine with which size the array should be allocated
  if (imax_new-imin_new<0) then
      ! everything in either first or second half
      imin = imin_old
      imax = imax_old
      !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
  else
      ! in both half
      if (imax_old-imin_old>imax_new-imin_new) then
          ! wrap around
          imin = imin_new
          imax = imax_new
          !sparsemat%offset_matrixindex_in_compressed_fortransposed = imin_new
      else
          ! no wrap around
          imin = imin_old
          imax = imax_old
          !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
      end if
  end if

  !!! Check
  !!if (sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1) then
  !!    stop 'sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1'
  !!end if

  nlen = imax - imin + 1
  sparsemat%offset_matrixindex_in_compressed_fortransposed = imin
  write(*,*) 'iproc, imin, imax, nlen', iproc, imin, imax, nlen

  !!! This is a temporary solution for spin polarized systems
  !!imax=min(imax,orbs%norbu)



  !!allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
  !!call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
  !sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/imin.to.imax,imin.to.imax/),&
  !    id='sparsemat%matrixindex_in_compressed_fortransposed')
  sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/nlen,nlen/),&
      id='sparsemat%matrixindex_in_compressed_fortransposed')

  !$omp parallel do default(private) shared(sparsemat,orbs,imin,imax)
  do iorb=imin,imax
      i = iorb - imin + 1
      do jorb=imin,imax
          j = jorb - imin + 1
          !@ii=(jorb-1)*sparsemat%nfvctr+iorb
          !@ispin=(ii-1)/sparsemat%nfvctr+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)
          !@iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
          !@jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
          !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iiorb,jjorb,orbs%norbu,sparsemat)
          iiorb = mod(iorb-1,sparsemat%nfvctr)+1
          jjorb = mod(jorb-1,sparsemat%nfvctr)+1
          !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
          sparsemat%matrixindex_in_compressed_fortransposed(i,j)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
          !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
          !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
      end do
  end do
  !$omp end parallel do

  !@! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
  !@if (ispin==2) then
  !@    matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
  !@end if

  call f_release_routine()

end subroutine init_matrixindex_in_compressed_fortransposed


subroutine synchronize_matrix_taskgroups(iproc, nproc, smat, mat)
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(sparse_matrix),intent(in) :: smat
  type(matrices),intent(in) :: mat

  ! Local variables
  integer :: ncount, itg, iitg, ispin, ishift, ist_send, ist_recv
  integer,dimension(:),allocatable :: request
  real(kind=8),dimension(:),allocatable :: recvbuf

  if (nproc>1) then
      request = f_malloc(smat%ntaskgroupp,id='request')
      ncount = 0
      do itg=1,smat%ntaskgroupp
          iitg = smat%taskgroupid(itg)
          ncount = ncount + smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
      end do
      recvbuf = f_malloc(ncount,id='recvbuf')
      do ispin=1,smat%nspin
          ishift = (ispin-1)*smat%nvctrp_tg

          ncount = 0
          do itg=1,smat%ntaskgroupp
              iitg = smat%taskgroupid(itg)
              ist_send = smat%taskgroup_startend(1,1,iitg) - smat%isvctrp_tg
              ist_recv = ncount + 1
              ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
              !!call mpi_iallreduce(mat%matrix_compr(ist_send), recvbuf(ist_recv), ncount, &
              !!     mpi_double_precision, mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg), ierr)
              if (nproc>1) then
                  call mpiiallred(mat%matrix_compr(ishift+ist_send), recvbuf(ist_recv), ncount, &
                       mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg))
              else
                  call vcopy(ncount, mat%matrix_compr(ishift+ist_send), 1, recvbuf(ist_recv), 1)
              end if
          end do
          if (nproc>1) then
              call mpiwaitall(smat%ntaskgroupp, request)
          end if
          ncount = 0
          do itg=1,smat%ntaskgroupp
              iitg = smat%taskgroupid(itg)
              ist_send = smat%taskgroup_startend(1,1,iitg) - smat%isvctrp_tg
              ist_recv = ncount + 1
              ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
              !call vcopy(ncount, recvbuf(ist_recv), 1, mat%matrix_compr(ishift+ist_send), 1)
              call dcopy(ncount, recvbuf(ist_recv), 1, mat%matrix_compr(ishift+ist_send), 1)
          end do
      end do
      call f_free(request)
      call f_free(recvbuf)
  end if
end subroutine synchronize_matrix_taskgroups
