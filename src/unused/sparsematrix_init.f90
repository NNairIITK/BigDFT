    !> Routine that initiakizes the sparsity pattern for a matrix, bases solely on
    !  the distance between the locreg centers
    subroutine init_sparsity_from_distance(iproc, nproc, orbs, lzd, input, &
               nseg, nsegline, istsegline, keyg, sparsemat)
      use module_base
      use module_types
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nseg
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(input_variables),intent(in) :: input
      integer,dimension(orbs%norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(out) :: sparsemat
    
      ! Local variables
      integer :: iorb, jorb, isegline, iseg, ivctr, ijorb, istat
      logical :: segment_started, overlap, newline
      character(len=*),parameter :: subname='init_sparsity_from_distance'
    
    
      !call nullify_sparse_matrix(sparsemat)
      sparsemat=sparse_matrix_null()
    
      sparsemat%nfvctr=orbs%norb
      sparsemat%nfvctrp=orbs%norbp
      sparsemat%isfvctr=orbs%isorb
    
      sparsemat%nfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nfvctr_par')
      sparsemat%isfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isfvctr_par')
    
      call vcopy(nproc,orbs%isorb_par(0),1,sparsemat%isfvctr_par(0),1)
      call vcopy(nproc,orbs%norb_par(0,0),1,sparsemat%nfvctr_par(0),1)
    
      sparsemat%nsegline=f_malloc_ptr(orbs%norb,id='sparsemat%nsegline')
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
    
    
      sparsemat%keyv=f_malloc_ptr(sparsemat%nseg,id='sparsemat%keyv')
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
    
    
      ! Initialize the parameters for the spare matrix matrix multiplication
      call init_sparse_matrix_matrix_multiplication(orbs%norb, orbs%norbp, orbs%isorb, nseg, &
               nsegline, istsegline, keyg, sparsemat)
      
    
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

