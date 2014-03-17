subroutine transform_sparse_matrix(smat, lmat, cmode)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(sparse_matrix),intent(inout) :: smat, lmat
  character(len=14),intent(in) :: cmode

  ! Local variables
  integer :: imode, icheck, isseg, isstart, isend, ilseg, ilstart, ilend
  integer :: iostart, ioend, ilength, isoffset, iloffset, iscostart, ilcostart, i
  integer,parameter :: SMALL_TO_LARGE=1
  integer,parameter :: LARGE_TO_SMALL=2


  ! determine the case:
  ! SMALL_TO_LARGE -> transform from large sparsity pattern to small one
  ! LARGE_TO_SMALL -> transform from small sparsity pattern to large one
  if (cmode=='small_to_large' .or. cmode=='SMALL_TO_LARGE') then
      imode=SMALL_TO_LARGE
  else if (cmode=='large_to_small' .or. cmode=='LARGE_TO_SMALL') then
      imode=LARGE_TO_SMALL
  else
      stop 'wrong cmode'
  end if

  select case (imode)
  case (SMALL_TO_LARGE)
      call to_zero(lmat%nvctr,lmat%matrix_compr(1))
  case (LARGE_TO_SMALL)
      call to_zero(smat%nvctr,smat%matrix_compr(1))
  case default
      stop 'wrong imode'
  end select


  icheck=0
  sloop: do isseg=1,smat%nseg
      isstart=smat%keyg(1,isseg)
      isend=smat%keyg(2,isseg)
      lloop: do ilseg=1,lmat%nseg
          ilstart=lmat%keyg(1,ilseg)
          ilend=lmat%keyg(2,ilseg)

          !write(*,*) 'isstart, isend, ilstart, ilend', isstart, isend, ilstart, ilend
          ! check whether there is an overlap:
          ! if not, increase loop counters
          if (ilstart>isend) exit lloop
          if (isstart>ilend) cycle lloop
          ! if yes, determine start end end of overlapping segment (in uncompressed form)
          iostart=max(isstart,ilstart)
          ioend=min(isend,ilend)
          !write(*,*) 'iostart, ioend', iostart, ioend
          ilength=ioend-iostart+1

          ! offset with respect to the starting point of the segment
          isoffset=iostart-smat%keyg(1,isseg)
          iloffset=iostart-lmat%keyg(1,ilseg)

          ! determine start end and of the overlapping segment in compressed form
          iscostart=smat%keyv(isseg)+isoffset
          ilcostart=lmat%keyv(ilseg)+iloffset

          ! copy the elements
          select case (imode)
          case (SMALL_TO_LARGE) 
              do i=0,ilength-1
                  lmat%matrix_compr(ilcostart+i)=smat%matrix_compr(iscostart+i)
              end do
          case (LARGE_TO_SMALL) 
              do i=0,ilength-1
                  smat%matrix_compr(iscostart+i)=lmat%matrix_compr(ilcostart+i)
              end do
          case default
              stop 'wrong imode'
          end select
          icheck=icheck+ilength
      end do lloop
  end do sloop

  ! all elements of the small matrix must have been processed, no matter in
  ! which direction the transformation has been executed
  if (icheck/=smat%nvctr) then
      write(*,'(a,2i8)') 'ERROR: icheck/=smat%nvctr', icheck, smat%nvctr
      stop
  end if

end subroutine transform_sparse_matrix
