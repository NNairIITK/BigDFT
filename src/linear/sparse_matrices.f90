subroutine initMatrixCompression2(iproc, nproc, ndim, lzd, at, input, orbs, noverlaps, overlaps, mad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ndim
  type(local_zone_descriptors),intent(in) :: lzd
  type(atoms_data),intent(in) :: at
  type(input_variables),intent(in) :: input
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(orbs%norb),intent(in) :: noverlaps
  integer,dimension(ndim,orbs%norb),intent(in) :: overlaps
  type(matrixDescriptors),intent(out) :: mad
  
  ! Local variables
  integer :: jproc, iorb, jorb, iiorb, jjorb, ijorb, jjorbold, istat, iseg, nseg, irow, irowold, isegline, ilr, jlr
  integer :: iwa, jwa, itype, jtype, ierr, nseglinemax, iall
  integer,dimension(:,:,:),pointer:: keygline
  logical :: seg_started
  real(kind=8) :: tt, cut
  logical,dimension(:,:),allocatable :: kernel_locreg
  character(len=*),parameter :: subname='initMatrixCompression'
!  integer :: ii, iseg
  
  call timing(iproc,'init_matrCompr','ON')

  call nullify_matrixDescriptors(mad)

  mad%nseg=0
  mad%nvctr=0
  jjorbold=-1
  irowold=0
  allocate(mad%nsegline(orbs%norb), stat=istat)
  call memocc(istat, mad%nsegline, 'mad%nsegline', subname)
  allocate(mad%istsegline(orbs%norb), stat=istat)
  call memocc(istat, mad%istsegline, 'mad%istsegline', subname)
  mad%nsegline=0
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
                  mad%nvctr=mad%nvctr+1

                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  if(irow/=irowold) then
                      ! We are in a new line
                      mad%nsegline(irow)=mad%nsegline(irow)+1
                      irowold=irow
                  end if

              else
                  ! There was a zero segment in between, i.e. we are in a new segment
                  mad%nseg=mad%nseg+1
                  mad%nvctr=mad%nvctr+1
                  jjorbold=jjorb
                  
                  ! Segments for each row
                  irow=(jjorb-1)/orbs%norb+1
                  mad%nsegline(irow)=mad%nsegline(irow)+1
                  irowold=irow
                  if (jorb==1) then
                      ! Starting segment for this line
                      mad%istsegline(iiorb)=mad%nseg
                  end if
              end if
          end do
      end do
  end do

  if (iproc==0) then
      write(*,'(a,i0)') 'total elements: ',orbs%norb**2
      write(*,'(a,i0)') 'non-zero elements: ',mad%nvctr
      write(*,'(a,f5.2,a)') 'sparsity: ',1.d2*dble(orbs%norb**2-mad%nvctr)/dble(orbs%norb**2),'%'
  end if

  nseglinemax=0
  do iorb=1,orbs%norb
      if(mad%nsegline(iorb)>nseglinemax) then
          nseglinemax=mad%nsegline(iorb)
      end if
  end do

  allocate(mad%keyv(mad%nseg), stat=istat)
  call memocc(istat, mad%keyv, 'mad%keyv', subname)
  allocate(mad%keyg(2,mad%nseg), stat=istat)
  call memocc(istat, mad%keyg, 'mad%keyg', subname)
  allocate(keygline(2,nseglinemax,orbs%norb), stat=istat)
  call memocc(istat, keygline, 'keygline', subname)

  nseg=0
  mad%keyv(1)=1
  jjorbold=-1
  irow=0
  isegline=0
  irowold=0
  keygline=0
  mad%keyg=0
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
                      mad%keyg(2,nseg)=jjorbold
                      keygline(2,isegline,irowold)=mod(jjorbold-1,orbs%norb)+1
                  end if
                  ! Now add the new segment.
                  nseg=nseg+1
                  mad%keyg(1,nseg)=jjorb
                  jjorbold=jjorb
                  if(nseg>1) then
                      mad%keyv(nseg) = mad%keyv(nseg-1) + mad%keyg(2,nseg-1) - mad%keyg(1,nseg-1) + 1
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
  mad%keyg(2,nseg)=jjorb
  keygline(2,isegline,orbs%norb)=mod(jjorb-1,orbs%norb)+1

  iall=-product(shape(keygline))*kind(keygline)
  deallocate(keygline, stat=istat)
  call memocc(istat, iall, 'keygline', subname)


  ! Initialize kernel_locreg
  allocate(kernel_locreg(orbs%norbp,orbs%norb), stat=istat)
  call memocc(istat, kernel_locreg, 'kernel_locreg', subname)
  allocate(mad%kernel_nseg(orbs%norb), stat=istat)
  call memocc(istat, mad%kernel_nseg, 'mad%kernel_nseg', subname)
  call to_zero(orbs%norb, mad%kernel_nseg(1))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      iwa=orbs%onwhichatom(iiorb)
      itype=at%iatype(iwa)
      mad%kernel_nseg(iiorb)=0
      seg_started=.false.
      do jjorb=1,orbs%norb
          jlr=orbs%inwhichlocreg(jjorb)
          jwa=orbs%onwhichatom(jjorb)
          jtype=at%iatype(jwa)
          tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 + &
               (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 + &
               (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
          cut = input%lin%kernel_cutoff(itype)+input%lin%kernel_cutoff(jtype)
          tt=sqrt(tt)
          if (tt<=cut) then
              kernel_locreg(iorb,jjorb)=.true.
              if (.not.seg_started) then
                  mad%kernel_nseg(iiorb)=mad%kernel_nseg(iiorb)+1
              end if
              seg_started=.true.
          else
              kernel_locreg(iorb,jjorb)=.false.
              seg_started=.false.
          end if
      end do
  end do
  call mpiallred(mad%kernel_nseg(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  allocate(mad%kernel_segkeyg(2,maxval(mad%kernel_nseg),orbs%norb), stat=istat)
  call memocc(istat, mad%kernel_segkeyg, 'mad%kernel_segkeyg', subname)
  call to_zero(2*maxval(mad%kernel_nseg)*orbs%norb, mad%kernel_segkeyg(1,1,1))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      iseg=0
      seg_started=.false.
      do jjorb=1,orbs%norb
          if(kernel_locreg(iorb,jjorb)) then
              if (.not.seg_started) then
                  iseg=iseg+1
                  mad%kernel_segkeyg(1,iseg,iiorb)=jjorb
              end if
              seg_started=.true.
          else
              if (seg_started) then
                  mad%kernel_segkeyg(2,iseg,iiorb)=jjorb-1
              end if
              seg_started=.false.
          end if
      end do
      if (seg_started) then
          mad%kernel_segkeyg(2,iseg,iiorb)=orbs%norb
      end if
  end do
  call mpiallred(mad%kernel_segkeyg(1,1,1), 2*maxval(mad%kernel_nseg)*orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  iall = -product(shape(kernel_locreg))*kind(kernel_locreg) 
  deallocate(kernel_locreg,stat=istat)
  call memocc(istat,iall,'kernel_locreg',subname)

  call timing(iproc,'init_matrCompr','OF')


end subroutine initMatrixCompression2






subroutine initCommsOrtho(iproc, nproc, nspin, lzd, orbs, locregShape, op) 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initCommsOrtho
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspin
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  character(len=1),intent(in) :: locregShape
  type(overlapParameters),intent(out) :: op

  ! Local variables
  !!integer :: iorb, jorb, iiorb, ierr, jjorb, nsub
  integer ::  istat
  character(len=*),parameter :: subname='initCommsOrtho'


  call timing(iproc,'init_commOrtho','ON')

  call nullify_overlapParameters(op)

  ! Allocate the arrays that count the number of overlaps per process (comon%noverlaps)
  ! and per orbital (op%noverlaps)
  allocate(op%noverlaps(orbs%norb), stat=istat)
  call memocc(istat, op%noverlaps, 'op%noverlaps',subname)

  ! Count how many overlaping regions each orbital / process has.
  if(locregShape=='c') then
     stop "ERROR: locregShape=='c' is deprecated!"
  else if(locregShape=='s') then
     call determine_overlap_from_descriptors(iproc, nproc, orbs, orbs, lzd, lzd, op)
  end if

  call timing(iproc,'init_commOrtho','OF')

end subroutine initCommsOrtho

