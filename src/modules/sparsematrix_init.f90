module sparsematrix_init
  use module_base
  use sparsematrix_base
  implicit none

  private


  !> Public routines
  public :: initSparseMatrix
  public :: init_sparsity_from_distance
  public :: compressed_index
  public :: matrixindex_in_compressed
  public :: check_kernel_cutoff

  contains


    !> Currently assuming square matrices
    subroutine initSparseMatrix(iproc, nproc, lzd, orbs, input, sparsemat)
      use module_base
      use module_types
      use module_interfaces
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      type(input_variables),intent(in) :: input
      type(sparse_matrix), intent(out) :: sparsemat
      
      ! Local variables
      character(len=*), parameter :: subname='initSparseMatrix'
      integer :: jproc, iorb, jorb, iiorb, jjorb, ijorb, jjorbold, istat, nseg, irow, irowold, isegline, ilr, segn, ind, iseg
      integer :: nseglinemax, iall, ierr, jorbe, jorbs, jorbold, ii
      !integer :: compressed_index, matrixindex_in_compressed
      integer :: matrixindex_in_compressed
      integer, dimension(:,:,:), pointer :: keygline
      integer, dimension(:), pointer :: noverlaps
      integer, dimension(:,:), pointer :: overlaps
      !integer, dimension(:,:), allocatable :: sendbuf, requests, iminmaxarr
      !integer :: ierr, imin, imax, irecv, isend, kproc
      
      call timing(iproc,'init_matrCompr','ON')
    
      sparsemat=sparsematrix_null()
      call determine_overlaps(iproc, nproc, lzd, orbs, 's', noverlaps, overlaps)
    
      sparsemat%nfvctr=orbs%norb
      sparsemat%nfvctrp=orbs%norbp
      sparsemat%isfvctr=orbs%isorb
      sparsemat%nfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nfvctr_par')
      sparsemat%isfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isfvctr_par')
      call vcopy(nproc,orbs%isorb_par(0),1,sparsemat%isfvctr_par(0),1)
      call vcopy(nproc,orbs%norb_par(0,0),1,sparsemat%nfvctr_par(0),1)
    
      sparsemat%nseg=0
      sparsemat%nvctr=0
      jjorbold=-1
      irowold=0
      call allocate_sparsematrix_basic(input%store_index, orbs%norb, nproc, sparsemat)
      sparsemat%nsegline=0
      do jproc=0,nproc-1
          do iorb=1,orbs%norb_par(jproc,0)
              iiorb=orbs%isorb_par(jproc)+iorb
              ilr=orbs%inWhichLocreg(iiorb)
              ijorb=(iiorb-1)*orbs%norb
              do jorb=1,noverlaps(iiorb)
                  jjorb=overlaps(jorb,iiorb)+ijorb
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
          !!write(*,'(a,i0)') 'total elements: ',orbs%norb**2
          !!write(*,'(a,i0)') 'non-zero elements: ',sparsemat%nvctr
          !!write(*,'(a,f5.2,a)') 'sparsity: ',1.d2*dble(orbs%norb**2-sparsemat%nvctr)/dble(orbs%norb**2),'%'
          call yaml_map('total elements',orbs%norb**2)
          call yaml_map('non-zero elements',sparsemat%nvctr)
          call yaml_map('sparsity in %',1.d2*dble(orbs%norb**2-sparsemat%nvctr)/dble(orbs%norb**2),fmt='(f5.2)')
      end if
    
      nseglinemax=0
      do iorb=1,orbs%norb
          if(sparsemat%nsegline(iorb)>nseglinemax) then
              nseglinemax=sparsemat%nsegline(iorb)
          end if
      end do
    
      call allocate_sparsematrix_keys(sparsemat)
    
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
              do jorb=1,noverlaps(iiorb)
                  jjorb=overlaps(jorb,iiorb)+ijorb
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
    
      iall=-product(shape(noverlaps))*kind(noverlaps)
      deallocate(noverlaps, stat=istat)
      call memocc(istat, iall, 'noverlaps', subname)
    
      iall=-product(shape(overlaps))*kind(overlaps)
      deallocate(overlaps, stat=istat)
      call memocc(istat, iall, 'overlaps', subname)
    
      if (input%store_index) then
          ! store the indices of the matrices in the sparse format
          sparsemat%store_index=.true.
    
          !!! initialize sparsemat%matrixindex_in_compressed
          do iorb=1,orbs%norb
             do jorb=1,orbs%norb
                sparsemat%matrixindex_in_compressed_arr(iorb,jorb)=compressed_index(iorb,jorb,orbs%norb,sparsemat)
             end do
          end do
    
      else
          ! Otherwise alwyas calculate them on-the-fly
          sparsemat%store_index=.false.
          !!nullify(sparsemat%matrixindex_in_compressed_arr)
      end if
    
      !!allocate(sparsemat%orb_from_index(2,sparsemat%nvctr), stat=istat)
      !!call memocc(istat, sparsemat%orb_from_index, 'sparsemat%orb_from_index', subname)
    
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
    
      ! parallelization of matrices, following same idea as norb/norbp/isorb
    
      !most equal distribution, but want corresponding to norbp for second column
      !call kpts_to_procs_via_obj(nproc,1,sparsemat%nvctr,sparsemat%nvctr_par)
      !sparsemat%nvctrp=sparsemat%nvctr_par(iproc)
      !sparsemat%isvctr=0
      !do jproc=0,iproc-1
      !    sparsemat%isvctr=sparsemat%isvctr+sparsemat%nvctr_par(jproc)
      !end do
      !sparsemat%isvctr_par=0
      !do jproc=0,nproc-1
      !    if(iproc==jproc) then
      !       sparsemat%isvctr_par(jproc)=sparsemat%isvctr
      !    end if
      !end do
      !if(nproc >1) call mpiallred(sparsemat%isvctr_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      do jproc=0,nproc-1
         jorbs=sparsemat%isfvctr_par(jproc)+1
         jorbold=0
         do ii=1,sparsemat%nvctr
            jorb = sparsemat%orb_from_index(2,ii)
            if (jorb/=jorbold .and. jorb==jorbs) then
               sparsemat%isvctr_par(jproc)=ii-1
               exit
            end if
            jorbold=jorb
         end do
      end do
      do jproc=0,nproc-1
         if (jproc==nproc-1) then
            sparsemat%nvctr_par(jproc)=sparsemat%nvctr-sparsemat%isvctr_par(jproc)
         else
            sparsemat%nvctr_par(jproc)=sparsemat%isvctr_par(jproc+1)-sparsemat%isvctr_par(jproc)
         end if
         if (iproc==jproc) sparsemat%isvctr=sparsemat%isvctr_par(jproc)
         if (iproc==jproc) sparsemat%nvctrp=sparsemat%nvctr_par(jproc)
      end do
      !do jproc=0,nproc-1
      !   write(*,*) iproc,jproc,sparsemat%isvctr_par(jproc),sparsemat%isvctr,sparsemat%nvctr_par(jproc),sparsemat%nvctrp,sparsemat%nvctr
      !end do
    
      ! 0 - none, 1 - mpiallred, 2 - allgather
      sparsemat%parallel_compression=0
      sparsemat%can_use_dense=.false.
    
    
      call timing(iproc,'init_matrCompr','OF')
    
    
    end subroutine initSparseMatrix



    subroutine determine_overlaps(iproc, nproc, lzd, orbs, locregShape, noverlaps, overlaps) 
      use module_base
      use module_types
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
      character(len=*),parameter :: subname='determine_overlaps'
    
    
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
    
    end subroutine determine_overlaps


    !> Count for each orbital and each process the number of overlapping orbitals.
    subroutine determine_overlap_from_descriptors(iproc, nproc, orbs, orbsig, lzd, lzdig, op_noverlaps, op_overlaps)
      use module_base
      use module_types
      use module_interfaces
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



    !> Routine that initiakizes the sparsity pattern for a matrix, bases solely on
    !  the distance between the locreg centers
    subroutine init_sparsity_from_distance(iproc, nproc, orbs, lzd, input, sparsemat)
      use module_base
      use module_types
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(input_variables),intent(in) :: input
      type(sparse_matrix),intent(out) :: sparsemat
    
      ! Local variables
      integer :: iorb, jorb, isegline, iseg, ivctr, ijorb, istat
      logical :: segment_started, overlap, newline
      character(len=*),parameter :: subname='init_sparsity_from_distance'
    
    
      !call nullify_sparsematrix(sparsemat)
      sparsemat=sparsematrix_null()
    
      sparsemat%nfvctr=orbs%norb
      sparsemat%nfvctrp=orbs%norbp
      sparsemat%isfvctr=orbs%isorb
    
      !allocate(sparsemat%nfvctr_par(0:nproc-1),stat=istat)
      !call memocc(istat, sparsemat%nfvctr_par, 'sparsemat%nfvctr_par', subname)
      sparsemat%nfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nfvctr_par')
      !allocate(sparsemat%isfvctr_par(0:nproc-1),stat=istat)
      !call memocc(istat, sparsemat%isfvctr_par, 'sparsemat%isfvctr_par', subname)
      sparsemat%isfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isfvctr_par')
    
      call vcopy(nproc,orbs%isorb_par(0),1,sparsemat%isfvctr_par(0),1)
      call vcopy(nproc,orbs%norb_par(0,0),1,sparsemat%nfvctr_par(0),1)
    
      !allocate(sparsemat%nsegline(orbs%norb),stat=istat)
      !call memocc(istat, sparsemat%nsegline, 'sparsemat%nsegline', subname)
      sparsemat%nsegline=f_malloc_ptr(orbs%norb,id='sparsemat%nsegline')
      !allocate(sparsemat%istsegline(orbs%norb),stat=istat)
      !call memocc(istat, sparsemat%istsegline, 'sparsemat%istsegline', subname)
      sparsemat%istsegline=f_malloc_ptr(orbs%norb,id='sparsemat%istsegline')
    
      sparsemat%nseg=0
      sparsemat%nvctr=0
    
      do iorb=1,orbs%norb
          ! Always start a new segment for each line
          segment_started=.false.
          isegline=0
          newline=.true.
          do jorb=1,orbs%norb
              overlap=check_overlap(iorb,jorb)
              if (overlap) then
                  if (segment_started) then
                      ! there is no "hole" in between, i.e. we are in the same segment
                      sparsemat%nvctr=sparsemat%nvctr+1
                  else
                      ! there was a "hole" in between, i.e. we are in a new segment
                      sparsemat%nseg=sparsemat%nseg+1
                      isegline=isegline+1
                      sparsemat%nvctr=sparsemat%nvctr+1
                      if (newline) sparsemat%istsegline(iorb)=sparsemat%nseg
                      newline=.false.
                  end if
                  segment_started=.true.
              else
                  segment_started=.false.
              end if
          end do
          sparsemat%nsegline(iorb)=isegline
      end do
    
      if (iproc==0) then
          call yaml_map('total elements',orbs%norb**2)
          call yaml_map('non-zero elements',sparsemat%nvctr)
          call yaml_map('sparsity in %',1.d2*dble(orbs%norb**2-sparsemat%nvctr)/dble(orbs%norb**2),fmt='(f5.2)')
      end if
    
    
      !allocate(sparsemat%keyv(sparsemat%nseg),stat=istat)
      !call memocc(istat, sparsemat%keyv, 'sparsemat%keyv', subname)
      sparsemat%keyv=f_malloc_ptr(sparsemat%nseg,id='sparsemat%keyv')
      !allocate(sparsemat%keyg(2,sparsemat%nseg),stat=istat)
      !call memocc(istat, sparsemat%keyg, 'sparsemat%keyg', subname)
      sparsemat%keyg=f_malloc_ptr((/2,sparsemat%nseg/),id='sparsemat%keyg')
    
    
      iseg=0
      ivctr=0
      do iorb=1,orbs%norb
          ! Always start a new segment for each line
          segment_started=.false.
          do jorb=1,orbs%norb
              overlap=check_overlap(iorb,jorb)
              ijorb=(iorb-1)*orbs%norb+jorb
              if (overlap) then
                  if (segment_started) then
                      ! there is no "hole" in between, i.e. we are in the same segment
                      ivctr=ivctr+1
                  else
                      ! there was a "hole" in between, i.e. we are in a new segment.
                      iseg=iseg+1
                      ivctr=ivctr+1
                      ! open the current segment
                      sparsemat%keyg(1,iseg)=ijorb
                      ! start of  the current segment
                      sparsemat%keyv(iseg)=ivctr
                  end if
                  segment_started=.true.
              else
                  if (segment_started) then
                      ! close the previous segment
                      sparsemat%keyg(2,iseg)=ijorb-1
                  end if
                  segment_started=.false.
              end if
          end do
          ! close the last segment on the line if necessary
          if (segment_started) then
              sparsemat%keyg(2,iseg)=iorb*orbs%norb
          end if
      end do
    
      ! check whether the number of segments and elements agrees
      if (iseg/=sparsemat%nseg) then
          write(*,'(a,2i8)') 'ERROR: iseg/=sparsemat%nseg', iseg, sparsemat%nseg
          stop
      end if
      if (ivctr/=sparsemat%nvctr) then
          write(*,'(a,2i8)') 'ERROR: ivctr/=sparsemat%nvctr', ivctr, sparsemat%nvctr
          stop
      end if
    
    
      call init_indices_in_compressed(input%store_index, orbs%norb, sparsemat)
      call init_orbs_from_index(sparsemat)
      call init_matrix_parallelization(iproc, nproc, sparsemat)
    
    
      
    
      contains
    
        function check_overlap(iiorb,jjorb)
          ! Calling arguments
          integer,intent(in) :: iiorb, jjorb
          logical :: check_overlap
    
          ! Local variables
          integer :: ilr, jlr
          real(kind=8) :: maxdist, dist
          
          ilr=orbs%inwhichlocreg(iiorb)
          jlr=orbs%inwhichlocreg(jjorb)
          maxdist = (lzd%llr(ilr)%locrad_kernel+lzd%llr(jlr)%locrad_kernel)**2
          dist = (lzd%llr(ilr)%locregcenter(1) - lzd%llr(jlr)%locregcenter(1))**2 &
                +(lzd%llr(ilr)%locregcenter(2) - lzd%llr(jlr)%locregcenter(2))**2 &
                +(lzd%llr(ilr)%locregcenter(3) - lzd%llr(jlr)%locregcenter(3))**2
          if (dist<=maxdist) then
              check_overlap=.true.
          else
              check_overlap=.false.
          end if
        end function check_overlap
    
    
    end subroutine init_sparsity_from_distance


    subroutine init_matrix_parallelization(iproc, nproc, sparsemat)
      use module_base
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(inout) :: sparsemat
    
      ! Local variables
      integer :: jproc, jorbs, jorbold, ii, jorb, istat
      character(len=*),parameter :: subname='init_matrix_parallelization'
    
      ! parallelization of matrices, following same idea as norb/norbp/isorb
      !allocate(sparsemat%nvctr_par(0:nproc-1),stat=istat)
      !call memocc(istat, sparsemat%nvctr_par, 'sparsemat%nvctr_par', subname)
      sparsemat%nvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nvctr_par')
      !allocate(sparsemat%isvctr_par(0:nproc-1),stat=istat)
      !call memocc(istat, sparsemat%isvctr_par, 'sparsemat%isvctr_par', subname)
      sparsemat%isvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isvctr_par')
    
      !most equal distribution, but want corresponding to norbp for second column
      !call kpts_to_procs_via_obj(nproc,1,sparsemat%nvctr,sparsemat%nvctr_par)
      !sparsemat%nvctrp=sparsemat%nvctr_par(iproc)
      !sparsemat%isvctr=0
      !do jproc=0,iproc-1
      !    sparsemat%isvctr=sparsemat%isvctr+sparsemat%nvctr_par(jproc)
      !end do
      !sparsemat%isvctr_par=0
      !do jproc=0,nproc-1
      !    if(iproc==jproc) then
      !       sparsemat%isvctr_par(jproc)=sparsemat%isvctr
      !    end if
      !end do
      !if(nproc >1) call mpiallred(sparsemat%isvctr_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      do jproc=0,nproc-1
         jorbs=sparsemat%isfvctr_par(jproc)+1
         jorbold=0
         do ii=1,sparsemat%nvctr
            jorb = sparsemat%orb_from_index(2,ii)
            if (jorb/=jorbold .and. jorb==jorbs) then
               sparsemat%isvctr_par(jproc)=ii-1
               exit
            end if
            jorbold=jorb
         end do
      end do
      do jproc=0,nproc-1
         if (jproc==nproc-1) then
            sparsemat%nvctr_par(jproc)=sparsemat%nvctr-sparsemat%isvctr_par(jproc)
         else
            sparsemat%nvctr_par(jproc)=sparsemat%isvctr_par(jproc+1)-sparsemat%isvctr_par(jproc)
         end if
         if (iproc==jproc) sparsemat%isvctr=sparsemat%isvctr_par(jproc)
         if (iproc==jproc) sparsemat%nvctrp=sparsemat%nvctr_par(jproc)
      end do
      !do jproc=0,nproc-1
      !   write(*,*) iproc,jproc,sparsemat%isvctr_par(jproc),sparsemat%isvctr,sparsemat%nvctr_par(jproc),sparsemat%nvctrp,sparsemat%nvctr
      !end do
    
      ! 0 - none, 1 - mpiallred, 2 - allgather
      sparsemat%parallel_compression=0
      sparsemat%can_use_dense=.false.
    
    end subroutine init_matrix_parallelization


    subroutine init_indices_in_compressed(store_index, norb, sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparseMatrix
      implicit none
    
      ! Calling arguments
      logical,intent(in) :: store_index
      integer,intent(in) :: norb
      type(sparse_matrix),intent(inout) :: sparsemat
    
      ! Local variables
      !integer :: iorb, jorb, compressed_index, istat
      integer :: iorb, jorb, istat
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


    !> Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
    function compressed_index(irow, jcol, norb, sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparseMatrix
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: irow, jcol, norb
      type(sparse_matrix),intent(in) :: sparsemat
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


    integer function matrixindex_in_compressed(sparsemat, iorb, jorb)
      use module_base
      use module_types
      use sparsematrix_base, only: sparseMatrix
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: sparsemat
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
          type(sparse_matrix),intent(in) :: sparsemat
        
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


    subroutine init_orbs_from_index(sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparseMatrix
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(inout) :: sparsemat
    
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


end module sparsematrix_init
