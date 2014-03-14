module sparsematrix_init
  use module_base
  use sparsematrix_base
  implicit none

  private


  !> Public routines
  public :: initSparseMatrix

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
      type(sparseMatrix), intent(out) :: sparsemat
      
      ! Local variables
      character(len=*), parameter :: subname='initSparseMatrix'
      integer :: jproc, iorb, jorb, iiorb, jjorb, ijorb, jjorbold, istat, nseg, irow, irowold, isegline, ilr, segn, ind, iseg
      integer :: nseglinemax, iall, ierr, jorbe, jorbs, jorbold, ii
      integer :: compressed_index, matrixindex_in_compressed
      integer, dimension(:,:,:), pointer :: keygline
      integer, dimension(:), pointer :: noverlaps
      integer, dimension(:,:), pointer :: overlaps
      !integer, dimension(:,:), allocatable :: sendbuf, requests, iminmaxarr
      !integer :: ierr, imin, imax, irecv, isend, kproc
      
      call timing(iproc,'init_matrCompr','ON')
    
      !call nullify_sparsematrix(sparsemat)
      sparsemat=sparsematrix_null()
      call initCommsOrtho(iproc, nproc, lzd, orbs, 's', noverlaps, overlaps)
    
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
    
    !!  ! for the calculation of overlaps and the charge density
    !!  imin=minval(collcom%indexrecvorbital_c)
    !!  imin=min(imin,minval(collcom%indexrecvorbital_f))
    !!  imax=maxval(collcom%indexrecvorbital_c)
    !!  imax=max(imax,maxval(collcom%indexrecvorbital_f))
    !!
    !!  allocate(collcom%matrixindex_in_compressed(orbs%norb,imin:imax), stat=istat)
    !!  call memocc(istat, collcom%matrixindex_in_compressed, 'collcom%matrixindex_in_compressed', subname)
    !!
    !!  allocate(sendbuf(orbs%norb,orbs%norbp), stat=istat)
    !!  call memocc(istat, sendbuf, 'sendbuf', subname)
    !!
    !!  do iorb=1,orbs%norbp
    !!      iiorb=orbs%isorb+iorb
    !!      do jorb=1,orbs%norb
    !!          sendbuf(jorb,iorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
    !!      end do
    !!  end do
    !!
    !!  allocate(iminmaxarr(2,0:nproc-1), stat=istat)
    !!  call memocc(istat, iminmaxarr, 'iminmaxarr', subname)
    !!  call to_zero(2*nproc, iminmaxarr(1,0))
    !!  iminmaxarr(1,iproc)=imin
    !!  iminmaxarr(2,iproc)=imax
    !!  call mpiallred(iminmaxarr(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    !!
    !!  allocate(requests(maxval(orbs%norb_par(:,0))*nproc,2), stat=istat)
    !!  call memocc(istat, requests, 'requests', subname)
    !!
    !!  if (nproc>1) then
    !!
    !!      isend=0
    !!      irecv=0
    !!      do jproc=0,nproc-1
    !!          do jorb=1,orbs%norb_par(jproc,0)
    !!              jjorb=jorb+orbs%isorb_par(jproc)
    !!              do kproc=0,nproc-1
    !!                  if (jjorb>=iminmaxarr(1,kproc) .and. jjorb<=iminmaxarr(2,kproc)) then
    !!                      ! send from jproc to kproc
    !!                      if (iproc==jproc) then
    !!                          isend=isend+1
    !!                          call mpi_isend(sendbuf(1,jorb), orbs%norb, &
    !!                               mpi_integer, kproc, jjorb, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
    !!                      end if
    !!                      if (iproc==kproc) then
    !!                          irecv=irecv+1
    !!                          call mpi_irecv(collcom%matrixindex_in_compressed(1,jjorb), orbs%norb, &
    !!                               mpi_integer, jproc, jjorb, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
    !!                      end if
    !!                  end if
    !!              end do
    !!          end do
    !!      end do
    !!
    !!      call mpi_waitall(isend, requests(1,1), mpi_statuses_ignore, ierr)
    !!      call mpi_waitall(irecv, requests(1,2), mpi_statuses_ignore, ierr)
    !!
    !!  else
    !!      call vcopy(orbs%norb*orbs%norbp, sendbuf(1,1), 1, collcom%matrixindex_in_compressed(1,1), 1)
    !!  end if
    !!
    !!  !!do iorb=imin,imax
    !!  !!    do jorb=1,orbs%norb
    !!  !!        write(200+iproc,*) iorb,jorb,collcom%matrixindex_in_compressed(jorb,iorb)
    !!  !!    end do
    !!  !!end do
    !!
    !!
    !!  iall=-product(shape(iminmaxarr))*kind(iminmaxarr)
    !!  deallocate(iminmaxarr, stat=istat)
    !!  call memocc(istat, iall, 'iminmaxarr', subname)
    !!
    !!  iall=-product(shape(requests))*kind(requests)
    !!  deallocate(requests, stat=istat)
    !!  call memocc(istat, iall, 'requests', subname)
    !!
    !!  iall=-product(shape(sendbuf))*kind(sendbuf)
    !!  deallocate(sendbuf, stat=istat)
    !!  call memocc(istat, iall, 'sendbuf', subname)
    
    
    
      call timing(iproc,'init_matrCompr','OF')
    
    
    end subroutine initSparseMatrix


end module sparsematrix_init
