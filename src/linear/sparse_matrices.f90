subroutine initSparseMatrix(iproc, nproc, lzd, at, input, orbs, collcom, sparsemat)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(atoms_data),intent(in) :: at
  type(input_variables),intent(in) :: input
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  type(sparseMatrix), intent(out) :: sparsemat
  
  ! Local variables
  integer :: jproc, iorb, jorb, iiorb, jjorb, ijorb, jjorbold, istat, iseg, nseg, irow, irowold, isegline, ilr, jlr
  integer :: iwa, jwa, itype, jtype, ierr, nseglinemax, iall
  integer,dimension(:,:,:),pointer:: keygline
  logical :: seg_started
  real(kind=8) :: tt, cut
  character(len=*),parameter :: subname='initSparseMatrix'
  
  call timing(iproc,'init_matrCompr','ON')

  call nullify_sparsematrix(sparsemat)
  call initCommsOrtho(iproc, nproc, input%nspin, lzd, orbs, 's', sparsemat%noverlaps, sparsemat%overlaps)

  sparsemat%nseg=0
  sparsemat%nvctr=0
  jjorbold=-1
  irowold=0
  allocate(sparsemat%nsegline(orbs%norb), stat=istat)
  call memocc(istat, sparsemat%nsegline, 'sparsemat%nsegline', subname)
  allocate(sparsemat%istsegline(orbs%norb), stat=istat)
  call memocc(istat, sparsemat%istsegline, 'sparsemat%istsegline', subname)
  sparsemat%nsegline=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc,0)
          iiorb=orbs%isorb_par(jproc)+iorb
          ilr=orbs%inWhichLocreg(iiorb)
          ijorb=(iiorb-1)*orbs%norb
          do jorb=1,sparsemat%noverlaps(iiorb)
              jjorb=sparsemat%overlaps(jorb,iiorb)+ijorb
              if(jjorb==jjorbold+1 .and. jorb/=1) then
                  ! There was no zero element in between, i.e. we are in the same segment.
                  jjorbold=jjorb
                  sparsemat%nvctr=sparsemat%nvctr+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      sparsemat%nsegline(irow)=sparsemat%nsegline(irow)+1
                      irowold=irow
                  end if

              else
                  ! There was a zero segment in between, i.e. we are in a new segment
                  sparsemat%nseg=sparsemat%nseg+1
                  sparsemat%nvctr=sparsemat%nvctr+1
                  jjorbold=jjorb
                  
                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  sparsemat%nsegline(irow)=sparsemat%nsegline(irow)+1
                  irowold=irow
                  if (jorb==1) then
                      ! Starting segment for this line
                      sparsemat%istsegline(iiorb)=sparsemat%nseg
                  end if
              end if
          end do
      end do
  end do

  if (iproc==0) then
      write(*,'(a,i0)') 'total elements: ',orbs%norb**2
      write(*,'(a,i0)') 'non-zero elements: ',sparsemat%nvctr
      write(*,'(a,f5.2,a)') 'sparsity: ',1.d2*dble(orbs%norb**2-sparsemat%nvctr)/dble(orbs%norb**2),'%'
  end if

  nseglinemax=0
  do iorb=1,orbs%norb
      if(sparsemat%nsegline(iorb)>nseglinemax) then
          nseglinemax=sparsemat%nsegline(iorb)
      end if
  end do

  allocate(sparsemat%keyv(sparsemat%nseg), stat=istat)
  call memocc(istat, sparsemat%keyv, 'sparsemat%keyv', subname)
  allocate(sparsemat%keyg(2,sparsemat%nseg), stat=istat)
  call memocc(istat, sparsemat%keyg, 'sparsemat%keyg', subname)
  allocate(keygline(2,nseglinemax,orbs%norb), stat=istat)
  call memocc(istat, keygline, 'keygline', subname)

  nseg=0
  sparsemat%keyv(1)=1
  jjorbold=-1
  irow=0
  isegline=0
  irowold=0
  keygline=0
  sparsemat%keyg=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc,0)
          iiorb=orbs%isorb_par(jproc)+iorb
          ilr=orbs%inWhichLocreg(iiorb)
          ijorb=(iiorb-1)*orbs%norb
          do jorb=1,sparsemat%noverlaps(iiorb)
              jjorb=sparsemat%overlaps(jorb,iiorb)+ijorb
              if(jjorb==jjorbold+1 .and. jorb/=1) then
                  ! There was no zero element in between, i.e. we are in the same segment.

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line, so close the last segment and start the new one
                      keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                      isegline=1
                      keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  end if
                  jjorbold=jjorb
              else
                  ! There was a zero segment in between, i.e. we are in a new segment.
                  ! First determine the end of the previous segment.
                  if(jjorbold>0) then
                      sparsemat%keyg(2,nseg)=jjorbold
                      keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                  end if
                  ! Now add the new segment.
                  nseg=nseg+1
                  sparsemat%keyg(1,nseg)=jjorb
                  jjorbold=jjorb
                  if(nseg>1) then
                      sparsemat%keyv(nseg) = sparsemat%keyv(nseg-1) + sparsemat%keyg(2,nseg-1) - sparsemat%keyg(1,nseg-1) + 1
                  end if

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      isegline=1
                      keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  else
                      ! We are in the same line
                      isegline=isegline+1
                      keygline(1,isegline,irow)=mod(jjorb-1,orbs%norb)+1
                      irowold=irow
                  end if
              end if
          end do
      end do
  end do
  ! Close the last segment
  sparsemat%keyg(2,nseg)=jjorb
  keygline(2,isegline,orbs%norb)=mod(jjorb-1,orbs%norb)+1

  iall=-product(shape(keygline))*kind(keygline)
  deallocate(keygline, stat=istat)
  call memocc(istat, iall, 'keygline', subname)

  call init_matrixindex_in_compressed(iproc, nproc, orbs, collcom, sparsemat)

  call timing(iproc,'init_matrCompr','OF')


end subroutine initSparseMatrix






subroutine initCommsOrtho(iproc, nproc, nspin, lzd, orbs, locregShape, noverlaps, overlaps) 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initCommsOrtho
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspin
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


subroutine init_matrixindex_in_compressed(iproc, nproc, orbs, collcom, sparsemat)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(inout) :: sparsemat
  type(collective_comms),intent(in) :: collcom
  
  ! Local variables
  integer :: ii, istat, iorb, iiorb, ilr, iall, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, ierr
  integer :: ipt, jproc, nvalp_c, nvalp_f, imin, imax, jorb, kproc, jjorb, isend, irecv
  integer :: compressed_index
  real(kind=8) :: tt, t1, t2
  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
  character(len=*),parameter :: subname='init_matrixindex_in_compressed'
  
  call timing(iproc,'init_collcomm ','ON')

  ! matrix index in the compressed format
  imin=minval(collcom%indexrecvorbital_c)
  imin=min(imin,minval(collcom%indexrecvorbital_f))
  imax=maxval(collcom%indexrecvorbital_c)
  imax=max(imax,maxval(collcom%indexrecvorbital_f))

  allocate(sparsemat%matrixindex_in_compressed(orbs%norb,imin:imax), stat=istat)
  call memocc(istat, sparsemat%matrixindex_in_compressed, 'sparsemat%matrixindex_in_compressed', subname)

  allocate(sendbuf(orbs%norb,orbs%norbp), stat=istat)
  call memocc(istat, sendbuf, 'sendbuf', subname)

  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jorb=1,orbs%norb
          sendbuf(jorb,iorb)=compressed_index(iiorb,jorb,orbs%norb, sparsemat)
      end do
  end do

  allocate(iminmaxarr(2,0:nproc-1), stat=istat)
  call memocc(istat, iminmaxarr, 'iminmaxarr', subname)
  call to_zero(2*nproc, iminmaxarr(1,0))
  iminmaxarr(1,iproc)=imin
  iminmaxarr(2,iproc)=imax
  call mpiallred(iminmaxarr(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  allocate(requests(maxval(orbs%norb_par(:,0))*nproc,2), stat=istat)
  call memocc(istat, requests, 'requests', subname)

  if (nproc>1) then

      isend=0
      irecv=0
      do jproc=0,nproc-1
          do jorb=1,orbs%norb_par(jproc,0)
              jjorb=jorb+orbs%isorb_par(jproc)
              do kproc=0,nproc-1
                  if (jjorb>=iminmaxarr(1,kproc) .and. jjorb<=iminmaxarr(2,kproc)) then
                      ! send from jproc to kproc
                      if (iproc==jproc) then
                          isend=isend+1
                          call mpi_isend(sendbuf(1,jorb), orbs%norb, &
                               mpi_integer, kproc, jjorb, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                      end if
                      if (iproc==kproc) then
                          irecv=irecv+1
                          call mpi_irecv(sparsemat%matrixindex_in_compressed(1,jjorb), orbs%norb, &
                               mpi_integer, jproc, jjorb, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                      end if
                  end if
              end do
          end do
      end do

      call mpi_waitall(isend, requests(1,1), mpi_statuses_ignore, ierr)
      call mpi_waitall(irecv, requests(1,2), mpi_statuses_ignore, ierr)

  else
      call vcopy(orbs%norb*orbs%norbp, sendbuf(1,1), 1, sparsemat%matrixindex_in_compressed(1,1), 1)
  end if

  iall=-product(shape(iminmaxarr))*kind(iminmaxarr)
  deallocate(iminmaxarr, stat=istat)
  call memocc(istat, iall, 'iminmaxarr', subname)

  iall=-product(shape(requests))*kind(requests)
  deallocate(requests, stat=istat)
  call memocc(istat, iall, 'requests', subname)

  iall=-product(shape(sendbuf))*kind(sendbuf)
  deallocate(sendbuf, stat=istat)
  call memocc(istat, iall, 'sendbuf', subname)
  
call timing(iproc,'init_collcomm ','OF')

  
end subroutine init_matrixindex_in_compressed




! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
function compressed_index(iiorb, jjorb, norb, sparsemat)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iiorb, jjorb, norb
  type(sparseMatrix),intent(in) :: sparsemat
  integer :: compressed_index

  ! Local variables
  integer :: ii, iseg

  ii=(iiorb-1)*norb+jjorb

  iseg=sparsemat%istsegline(iiorb)
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


subroutine compress_matrix_for_allreduce(n, sparsemat, mat, mat_compr)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: n
  type(sparseMatrix),intent(in) :: sparsemat
  real(kind=8),dimension(n**2),intent(in) :: mat
  real(kind=8),dimension(sparsemat%nvctr),intent(out) :: mat_compr

  ! Local variables
  integer :: jj, iseg, jorb

  !$omp parallel do default(private) shared(sparsemat,mat_compr,mat)
  do iseg=1,sparsemat%nseg
      jj=1
      do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
          mat_compr(sparsemat%keyv(iseg)+jj-1)=mat(jorb)
          jj=jj+1
      end do
  end do
  !$omp end parallel do

end subroutine compress_matrix_for_allreduce





subroutine uncompressMatrix(norb, sparsemat, lmat, mat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: norb
  type(sparseMatrix), intent(in) :: sparsemat
  real(kind=8), dimension(sparsemat%nvctr), intent(in) :: lmat
  real(kind=8), dimension(norb**2), intent(out) :: mat
  
  ! Local variables
  integer :: iseg, ii, jorb

  
  call to_zero(norb**2, mat(1))
  
  !$omp parallel do default(private) shared(sparsemat,lmat,mat)

  do iseg=1,sparsemat%nseg
      ii=0
      do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
          mat(jorb)=lmat(sparsemat%keyv(iseg)+ii)
          ii=ii+1
      end do
  end do

  !$omp end parallel do
  
end subroutine uncompressMatrix
