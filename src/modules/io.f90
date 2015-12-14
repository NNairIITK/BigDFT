module io
  use public_enums
  implicit none

  private

  !> Public routines
  public :: read_linear_matrix_dense
  public :: writeonewave_linear
  public :: writemywaves_linear
  public :: writemywaves_linear_fragments
  public :: read_coeff_minbasis
  public :: io_read_descr_linear
  public :: write_sparse_matrix
  public :: write_dense_matrix
  public :: write_linear_matrices
  public :: writeLinearCoefficients
  public :: write_linear_coefficients
  public :: read_linear_coefficients
  public :: write_partial_charges

  public :: io_error, io_warning, io_open
  public :: io_read_descr, read_psi_compress
  public :: io_gcoordToLocreg
  public :: read_psig
  public :: read_sparse_matrix
  public :: read_dense_matrix
  public :: dist_and_shift
  public :: find_neighbours
  public :: write_ccs_matrix


  contains


    !> Write all my wavefunctions in files by calling writeonewave
    subroutine writemywaves_linear(iproc,filename,iformat,npsidim,Lzd,orbs,nelec,at,rxyz,psi,nfvctr,coeff)
      use module_types
      use module_base
      use yaml_output
      use module_interfaces, only: open_filename_of_iorb
      implicit none
      integer, intent(in) :: iproc,iformat,npsidim,nelec,nfvctr
      !integer, intent(in) :: norb   !< number of orbitals, not basis functions
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(in) :: orbs         !< orbs describing the basis functions
      type(local_zone_descriptors), intent(in) :: Lzd
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      real(wp), dimension(npsidim), intent(in) :: psi  ! Should be the real linear dimension and not the global
      real(wp), dimension(nfvctr,orbs%norb), intent(in) :: coeff
      character(len=*), intent(in) :: filename
      !Local variables
      logical :: binary, is_etsf
      integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr,shift,ii,iat,unitwf
      integer :: jorb,jlr
      real(kind=4) :: tr0,tr1
      real(kind=8) :: tel
    
      unitwf=99
      binary=(iformat/=WF_FORMAT_PLAIN)
      is_etsf=(iformat==WF_FORMAT_ETSF)
    
      if (iproc == 0) call yaml_map('Write wavefunctions to file', trim(filename)//'.*')
      !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
    
      !if (binary) then
      if (is_etsf) then
         call f_err_throw('Linear scaling with ETSF writing not implemented yet')
    !     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
      else
         call cpu_time(tr0)
         call system_clock(ncount1,ncount_rate,ncount_max)
    
         ! Write the TMBs in the Plain BigDFT files.
         ! Use same ordering as posinp and llr generation
         ii = 0
         do iat = 1, at%astruct%nat
            do iorb=1,orbs%norbp
               if(iat == orbs%onwhichatom(iorb+orbs%isorb)) then
                  shift = 1
                  do jorb = 1, iorb-1 
                     jlr = orbs%inwhichlocreg(jorb+orbs%isorb)
                     shift = shift + Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f
                  end do
                  ii = ii + 1
                  ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
                  do ispinor=1,orbs%nspinor
                     call open_filename_of_iorb(unitwf,binary,filename, &
                        & orbs,iorb,ispinor,iorb_out)
                     call writeonewave_linear(unitwf,.not. binary,iorb_out,&
                        & Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,&
                        & Lzd%Llr(ilr)%ns1,Lzd%Llr(ilr)%ns2,Lzd%Llr(ilr)%ns3,& 
                        & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                        & Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
                        & at%astruct%nat,rxyz,Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                        & Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                        & Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                        & Lzd%Llr(ilr)%wfd%keygloc(1:,Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & Lzd%Llr(ilr)%wfd%keyvloc(Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & psi(shift),psi(Lzd%Llr(ilr)%wfd%nvctr_c+shift),orbs%eval(iorb+orbs%isorb),&
                        & orbs%onwhichatom(iorb+orbs%isorb))
                     call f_close(unitwf)
                  end do
               end if
            enddo
         end do
    
        ! Now write the coefficients to file
        ! Must be careful, the orbs%norb is the number of basis functions
        ! while the norb is the number of orbitals.
        if(iproc == 0) then
           call f_open_file(unitwf,file=filename//'_coeff.bin',&
                binary=binary)
           !if(iformat == WF_FORMAT_PLAIN) then
           !  open(99, file=filename//'_coeff.bin', status='unknown',form='formatted')
           !else
           !open(99, file=filename//'_coeff.bin', status='unknown',form='unformatted')
           !end if
          call writeLinearCoefficients(unitwf,.not. binary,at%astruct%nat,rxyz,orbs%norb,&
               nelec,nfvctr,coeff,orbs%eval)
          call f_close(unitwf)
       end if
         call cpu_time(tr1)
         call system_clock(ncount2,ncount_rate,ncount_max)
         tel=dble(ncount2-ncount1)/dble(ncount_rate)
         if (iproc == 0) then
            call yaml_sequence_open('Write Waves Time')
            call yaml_sequence(advance='no')
            call yaml_mapping_open(flow=.true.)
            call yaml_map('Process',iproc)
            call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
            call yaml_mapping_close()
            call yaml_sequence_close()
         end if
         !write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
         !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
      end if
    
    END SUBROUTINE writemywaves_linear


    !> Write all my wavefunctions for fragments
    ! NB spin needs fixing!
    subroutine writemywaves_linear_fragments(iproc,filename,iformat,npsidim,Lzd,orbs,nelec,at,rxyz,psi,coeff,&
         dir_output,input_frag,ref_frags,linmat,methTransformOverlap,max_inversion_error,orthpar,num_neighbours)
      use module_types
      use module_atoms, only: astruct_dump_to_file, deallocate_atomic_structure, nullify_atomic_structure
      use module_base
      use module_fragments
      use yaml_output
      use module_interfaces, only: open_filename_of_iorb
      use sparsematrix_base, only: sparsematrix_malloc_ptr, assignment(=), DENSE_FULL, &
           matrices_null, allocate_matrices, deallocate_matrices
      use sparsematrix, only: uncompress_matrix2
      !use matrix_operations, only: overlapPowerGeneral, check_taylor_order
      implicit none
      integer, intent(in) :: iproc,iformat,npsidim,nelec
      !integer, intent(in) :: norb   !< number of orbitals, not basis functions
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(in) :: orbs         !< orbs describing the basis functions
      type(local_zone_descriptors), intent(in) :: Lzd
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      real(wp), dimension(npsidim), intent(in) :: psi  ! Should be the real linear dimension and not the global
      real(wp), dimension(orbs%norb,orbs%norb), intent(in) :: coeff !SM: IS this correcy even with spin?
      character(len=*), intent(in) :: dir_output, filename
      type(fragmentInputParameters), intent(in) :: input_frag
      type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
      type(linear_matrices), intent(inout) :: linmat
      integer,intent(inout) :: methTransformOverlap
      real(kind=8),intent(in) :: max_inversion_error
      type(orthon_data),intent(in) :: orthpar
      integer, intent(in) :: num_neighbours
      !Local variables
      integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr,shift,ii,iat
      integer :: jorb,jlr,isforb,isfat,ifrag,ifrag_ref,iforb,iiorb,iorbp,iiat,unitwf,ityp,iatt
      integer :: onwhichatom_frag,ifr,iatf,ntmb_frag_and_env,nelec_frag,ifr_ref,ia,ja,io,jo,unitm,iatnf
      integer :: num_neighbours_tot, num_after, num_before, jat
      integer, allocatable, dimension(:,:) :: map_frag_and_env, frag_map   
      real(kind=8), allocatable, dimension(:,:) :: coeff_frag, kernel_frag, rxyz_frag
      real(kind=8), allocatable, dimension(:) :: eval_frag
      real(kind=4) :: tr0,tr1
      real(kind=8) :: tel
      character(len=256) :: full_filename
      logical, allocatable, dimension(:) :: fragment_written
      logical :: binary
      !logical, parameter :: write_overlap=.true. ! want this to include rotation taking into account environment?  also do some check after restart? - not sure this is possible

      type(matrices),dimension(1) :: ovrlp_half_
      real(kind=8) :: mean_error, max_error
      logical :: calc_sks=.false.
      real(kind=8), dimension(:,:,:), pointer :: kernel_orthog
    
      if (iproc == 0) call yaml_map('Write wavefunctions to file', trim(filename)//'.*')
      !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
    
      if (iformat == WF_FORMAT_ETSF) then
          stop 'Linear scaling with ETSF writing not implemented yet'
      !    call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
      else
         call cpu_time(tr0)
         call system_clock(ncount1,ncount_rate,ncount_max)
    
         ! Write the TMBs in the Plain BigDFT files.
         ! For now only output one (first) set per reference fragment
        
         ! array to check if we already outputted tmbs for this fragment type
         fragment_written=f_malloc((/input_frag%nfrag_ref/),id='fragment_written')
         fragment_written=.false.
    
         unitwf=99
         isforb=0
         isfat=0
         loop_ifrag: do ifrag=1,input_frag%nfrag
            ! find reference fragment this corresponds to and check if we already outputted tmbs for this reference fragment
            ifrag_ref=input_frag%frag_index(ifrag)

            if (fragment_written(ifrag_ref)) then
               isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
               isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat   
               cycle
            end if
            fragment_written(ifrag_ref)=.true.

            ! loop over orbitals of this fragment
            !loop_iforb: do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb

            ! always have at least 1 tmb/atom so this is safe
            frag_map=f_malloc0((/ref_frags(ifrag_ref)%fbasis%forbs%norb,3/),id='frag_map')

               loop_iorb: do iorbp=1,orbs%norbp
                  iiorb=iorbp+orbs%isorb
    
                  ! check if this ref frag orbital corresponds to the orbital we want
                  !if (iiorb/=iforb+isforb) cycle
                  !NO LONGER TRUE - can reactivate the above once input tmb order is guaranteed to be in fragment order
                  !for now do more complicated atom based testing
                  !This also means coeffs and evals are printed wrong, but we're neglecting these anyway so ignore until we tidy
                  ilr=orbs%inwhichlocreg(iiorb)
                  iiat=orbs%onwhichatom(iiorb)

                  if (iiat<=isfat .or. iiat>isfat+ref_frags(ifrag_ref)%astruct_frg%nat) cycle
                  ! there might be an easier way to figure this out, but in calculations from scratch there's no guarantee the tmbs are arranged by fragment order
                  ! alternatively modify so they are in this order, but not sure there's any benefit to that approach?
                  iat=0
                  onwhichatom_frag=-1
                  do ifr=1,input_frag%nfrag
                     ifr_ref=input_frag%frag_index(ifr)
                     do iatf=1,ref_frags(ifr_ref)%astruct_frg%nat
                        iat=iat+1
                        if (iat==iiat) then
                           onwhichatom_frag=iatf
                           ! double check
                           if (ifr/=ifrag) stop 'Error error abort warning' 
                           exit
                        end if
                     end do
                     if (onwhichatom_frag/=-1) exit
                  end do

                  iforb=0
                  do io=1,iiorb
                     ia=orbs%onwhichatom(io)
                     if (ia>isfat .and. ia<=isfat+ref_frags(ifrag_ref)%astruct_frg%nat) iforb=iforb+1
                  end do
                  !totally arbitrary order here but at least we know we're on the right fragment!
                  !iforb = iforb+1
                  !tmb frag -> tmb full
                  frag_map(iforb,1)=iiorb
                  !tmb frag -> atom frag
                  frag_map(iforb,2)=onwhichatom_frag
                  !atom frag -> atom full
                  !fill this in after allreduce otherwise might be modifying on more than one mpi
                  !frag_map(onwhichatom_frag,3)=iiat

                  !debug
                  !print*,'iiorb,ifrag,ifrag_ref,iiat,onwhichatom_frag',iiorb,ifrag,ifrag_ref,iiat,onwhichatom_frag,iforb
                  !write(*,'(A,7(1x,I4))') 'iproc,iiorb,ifrag,ifrag_ref,iiat,onwhichatom_frag',&
                  !     iproc,iiorb,ifrag,ifrag_ref,iiat,onwhichatom_frag,iforb

                  shift = 1
                  do jorb = 1, iorbp-1 
                     jlr = orbs%inwhichlocreg(jorb+orbs%isorb)
                     shift = shift + Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f
                  end do
    
                  loop_ispinor: do ispinor=1,orbs%nspinor
                     ! as this is a fragment calculation frag%dirname should contain fragment directory (otherwise it would be empty - should add check)
                     ! bit of a hack to use orbs here not forbs, but different structures so this is necessary - to clean somehow
                     full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)
    
                     call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),full_filename, &
                          & orbs,iorbp,ispinor,iorb_out,iforb)
    
                     !also what to do with eval? - at the moment completely arbitrary
                     call writeonewave_linear(unitwf,(iformat == WF_FORMAT_PLAIN),iorb_out,&
                        & Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,&
                        & Lzd%Llr(ilr)%ns1,Lzd%Llr(ilr)%ns2,Lzd%Llr(ilr)%ns3,& 
                        & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                        & Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
                        & ref_frags(ifrag_ref)%astruct_frg%nat,rxyz(:,isfat+1:isfat+ref_frags(ifrag_ref)%astruct_frg%nat),&
                        & Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                        & Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                        & Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                        & Lzd%Llr(ilr)%wfd%keygloc(1:,Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & Lzd%Llr(ilr)%wfd%keyvloc(Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & psi(shift),psi(Lzd%Llr(ilr)%wfd%nvctr_c+shift),-0.5d0, & !orbs%eval(iiorb),&
                        & onwhichatom_frag)
                        !& orbs%onwhichatom(iiorb)-isfat) ! only works if reading the rewriting fragment tmbs
    
                     close(unitwf)
    
                  end do loop_ispinor
               end do loop_iorb
            !end do loop_iforb

            if (bigdft_mpi%nproc > 1) then
               call mpiallred(frag_map, mpi_sum, comm=bigdft_mpi%mpi_comm)
            end if

            ! reconstruct atom->atom mapping part
            do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
               !frag_map(iforb,1)=iiorb
               !frag_map(iforb,2)=onwhichatom_frag
               iiorb=frag_map(iforb,1)
               iiat=orbs%onwhichatom(iiorb)
               onwhichatom_frag=frag_map(iforb,2)
               frag_map(onwhichatom_frag,3)=iiat
            end do

            !debug
            !if (iproc==0) then
            !   do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
            !      if (iorb<=ref_frags(ifrag_ref)%astruct_frg%nat) then
            !         write(*,'(A,6(1x,I4),3(1x,F12.6))') 'iproc,ifrag,ifrag_ref',iproc,ifrag,ifrag_ref,frag_map(iorb,:),rxyz(:,frag_map(iorb,3))
            !      else
            !         write(*,'(A,6(1x,I4))') 'iproc,ifrag,ifrag_ref',iproc,ifrag,ifrag_ref,frag_map(iorb,:)
            !      end if
            !   end do
            !end if

            ! write environment coordinates 
            ! want to avoid several MPI ranks writing to the same file
            ! and since the number of reference frags is expected to be relatively small as well as the file sizes,
            ! write all files from iproc=0
            ! make use of frag_map(onwhichatom_frag,3)=iiat (atom frag -> atom full)
            if (num_neighbours/=0) then
               ! column 1 tmbs->tmbs, column 2 tmbs->atoms, column 3 atoms->atoms (as for frag_map)
               ! size orbs%norb is overkill but we don't know norb for env yet
               map_frag_and_env = f_malloc((/orbs%norb,3/),id='map_frag_and_env')

               call find_neighbours(num_neighbours,at,rxyz,orbs,ref_frags(ifrag_ref),frag_map,ntmb_frag_and_env,map_frag_and_env)
               !full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)//'_env'
               full_filename=trim(dir_output)//trim(input_frag%label(ifrag_ref))//'_env'
               if (iproc==0) then
                  !open file to make sure we overwrite rather than append
                  open(99,file=trim(full_filename)//'.xyz',status='replace')
                  call astruct_dump_to_file(ref_frags(ifrag_ref)%astruct_env,full_filename,'# fragment environment',unit=99)
                  close(99)
               end if
            end if
   
            ! NEED to think about this - just make it diagonal for now? or random?  or truncate so they're not normalized?  or normalize after truncating?
            ! Or maybe don't write coeffs at all but assume we're always doing frag to frag and can use isolated frag coeffs?
    
            ! Now write the coefficients to file
            ! Must be careful, the orbs%norb is the number of basis functions
            ! while the norb is the number of orbitals.

            ! in from scratch case ref_frags(ifrag_ref)%nelec will be empty (zero)
            ! need to recalculate in this case, so may as well recalculate anyway
            ! probably also an easier way to do this...
            nelec_frag=0
            do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
               !ityp=ref_frags(ifrag_ref)%astruct_frg%iatype(iat)
               ityp=at%astruct%iatype(iat+isfat)
               nelec_frag=nelec_frag+at%nelpsp(ityp)
            end do

            if (nelec_frag/=ref_frags(ifrag_ref)%nelec .and. ref_frags(ifrag_ref)%nelec/=0) &
                 call f_err_throw('Problem with nelec in fragment output')

            ! order of tmbs not guaranteed so can't do this the simple way
            ! in fact it makes no sense to be cropping coeffs at all in this case as ks states not divided between fragments
            ! but don't want to keep all of them
            ! so either do cropped kernel (to be added) or support function 'charges' as diagonal part
            ! for the latter assume orthogonal and just take K_aa

            ! should eventually extract this from sparse?
            linmat%kernel_%matrix = sparsematrix_malloc_ptr(linmat%l,iaction=DENSE_FULL,id='linmat%kernel_%matrix')
    
            call uncompress_matrix2(iproc, bigdft_mpi%nproc, linmat%l, &
                 linmat%kernel_%matrix_compr, linmat%kernel_%matrix)

            !might be better to move this outside of the routine?
            !if (calc_sks) then
            !   ovrlp_half_(1) = matrices_null()
            !   call allocate_matrices(linmat%l, allocate_full=.true., matname='ovrlp_half_(1)', mat=ovrlp_half_(1))

            !   !for case where tmbs aren't properly orthogonal, calculate S^1/2 K S^1/2
            !   call overlapPowerGeneral(iproc, bigdft_mpi%nproc, methTransformOverlap, 1, (/2/), &
            !        orthpar%blocksize_pdgemm, &
            !        imode=1, ovrlp_smat=linmat%s, inv_ovrlp_smat=linmat%l, &
            !        ovrlp_mat=linmat%ovrlp_, inv_ovrlp_mat=ovrlp_half_, &
            !        check_accur=.true., mean_error=mean_error, max_error=max_error)!!, &
            !   !if (iproc==0) call yaml_map('max error',max_error)
            !   !if (iproc==0) call yaml_map('mean error',mean_error)
            !   call check_taylor_order(mean_error, max_inversion_error, methTransformOverlap)

            !   call uncompress_matrix2(iproc, bigdft_mpi%nproc, linmat%l, &
            !        ovrlp_half_(1)%matrix_compr, ovrlp_half_(1)%matrix)

            !   kernel_orthog=f_malloc_ptr(src_ptr=linmat%kernel_%matrix,id='kernel_orthog')

            !   ! do this with sparse instead... fix after testing
            !   call calculate_coeffMatcoeff(bigdft_mpi%nproc,linmat%kernel_%matrix,orbs,orbs,ovrlp_half_(1)%matrix,kernel_orthog)
            !else
               kernel_orthog => linmat%kernel_%matrix
            !end if

            coeff_frag=f_malloc0((/ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%fbasis%forbs%norb/),id='coeff_frag')
            do io=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
               coeff_frag(io,io)=kernel_orthog(frag_map(io,1),frag_map(io,1),1)!ispin)
            end do
   
            ! want this to make sense also
            rxyz_frag=f_malloc0((/3,ref_frags(ifrag_ref)%astruct_frg%nat/),id='rxyz_frag')
            do ia=1,ref_frags(ifrag_ref)%astruct_frg%nat
               rxyz_frag(:,ia)=rxyz(:,frag_map(ia,3))
            end do

            ! evals irrelevant so fill with fake value
            eval_frag=f_malloc0(ref_frags(ifrag_ref)%fbasis%forbs%norb,id='eval_frag')
            do io=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
               eval_frag(io)=-0.5d0
            end do

            !kernel_frag=f_malloc0((/ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%fbasis%forbs%norb/),id='kernel_frag')
            !do io=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
            !   do jo=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
            !      kernel_frag(io,jo)=kernel_orthog(frag_map(io,1),frag_map(jo,1),1)!ispin)
            !   end do
            !end do

            if(iproc == 0) then
               full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)
     
               call f_open_file(unitwf,file=trim(full_filename)//'_coeff.bin',&
                    binary=(iformat /= WF_FORMAT_PLAIN))
               !if(iformat == WF_FORMAT_PLAIN) then
               !   open(unitwf, file=trim(full_filename)//'_coeff.bin', status='unknown',form='formatted')
               !else
               !   open(unitwf, file=trim(full_filename)//'_coeff.bin', status='unknown',form='unformatted')
               !end if
               
               ! Not sure whether this is correct for nspin=2...

               call writeLinearCoefficients(unitwf, (iformat == WF_FORMAT_PLAIN), ref_frags(ifrag_ref)%astruct_frg%nat, &
                    rxyz_frag, ref_frags(ifrag_ref)%fbasis%forbs%norb, nelec_frag, &
                    ref_frags(ifrag_ref)%fbasis%forbs%norb, coeff_frag, eval_frag)

               !call writeLinearCoefficients(unitwf,(iformat == WF_FORMAT_PLAIN),ref_frags(ifrag_ref)%astruct_frg%nat,&
               !     rxyz(:,isfat+1:isfat+ref_frags(ifrag_ref)%astruct_frg%nat),ref_frags(ifrag_ref)%fbasis%forbs%norb,&
               !     nelec_frag,&
               !     ref_frags(ifrag_ref)%fbasis%forbs%norb, &
               !     coeff(isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb,&
               !     isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb),&
               !     orbs%eval(isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb)) !-0.5d0
               call f_close(unitwf)
            end if
            call f_free(coeff_frag)
            !call f_free(kernel_frag)
            call f_free(eval_frag)


            ! put this in a routine and reuse in write_linear_matrices
            ! FIX SPIN
            unitm=99
            binary=(iformat /= WF_FORMAT_PLAIN)
            if(iproc == 0) then
               full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))
               call f_open_file(unitm, file=trim(full_filename)//'density_kernel.bin',&
                    binary=binary)
    
               if (.not. binary) then
                   write(unitm,'(a,3i10,a)') '#  ', ref_frags(ifrag_ref)%fbasis%forbs%norb, &
                       ref_frags(ifrag_ref)%astruct_frg%nat, linmat%m%nspin, &
                       '    number of basis functions, number of atoms, number of spins'
               else
                   write(unitm) '#  ', ref_frags(ifrag_ref)%fbasis%forbs%norb, &
                       ref_frags(ifrag_ref)%astruct_frg%nat, linmat%m%nspin, &
                       '    number of basis functions, number of atoms, number of spins'
               end if
               do ia=1,ref_frags(ifrag_ref)%astruct_frg%nat
                   if (.not. binary) then
                       write(unitm,'(a,3es24.16)') '#  ',rxyz_frag(1:3,ia)
                   else
                       write(unitm) '#  ',rxyz_frag(1:3,ia)
                   end if
               end do
 
               ! need to fix spin
               !do ispin=1,linmat%m%nspin
                  do io=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                     ia=frag_map(io,2)
                     do jo=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                        ja=frag_map(jo,2)
                        if (.not. binary) then
                           write(unitm,'(2(i6,1x),e19.12,2(1x,i6))') io,jo,&
                                kernel_orthog(frag_map(io,1),frag_map(jo,1),1),ia,ja
                        else
                           write(unitm) io,jo,kernel_orthog(frag_map(io,1),frag_map(jo,1),1),ia,ja
                        end if
                     end do
                  end do
               !end do
               call f_close(unitm)
            end if

            call f_free(rxyz_frag)

            ! also output 'environment' kernel
            if (ref_frags(ifrag_ref)%astruct_env%nat/=ref_frags(ifrag_ref)%astruct_frg%nat) then
               ! FIX SPIN
               unitm=99
               binary=(iformat /= WF_FORMAT_PLAIN)
               if(iproc == 0) then
                  full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))
                  call f_open_file(unitm, file=trim(full_filename)//'density_kernel_env.bin',&
                       binary=binary)
    
                  if (.not. binary) then
                      write(unitm,'(a,3i10,a)') '#  ', ntmb_frag_and_env, ref_frags(ifrag_ref)%astruct_env%nat, &
                          linmat%m%nspin, '    number of basis functions, number of atoms, number of spins'
                  else
                      write(unitm) '#  ', ntmb_frag_and_env, ref_frags(ifrag_ref)%astruct_env%nat, &
                          linmat%m%nspin, '    number of basis functions, number of atoms, number of spins'
                  end if
                  do ia=1,ref_frags(ifrag_ref)%astruct_env%nat
                      if (.not. binary) then
                          write(unitm,'(a,3es24.16)') '#  ',ref_frags(ifrag_ref)%astruct_env%rxyz(1:3,ia)
                      else
                          write(unitm) '#  ',ref_frags(ifrag_ref)%astruct_env%rxyz(1:3,ia)
                      end if
                  end do
 
                  ! need to fix spin
                  !do ispin=1,linmat%m%nspin
                     do io=1,ntmb_frag_and_env
                        ia=map_frag_and_env(io,2)
                        do jo=1,ntmb_frag_and_env
                           ja=map_frag_and_env(jo,2)
                           if (.not. binary) then
                              write(unitm,'(2(i6,1x),e19.12,2(1x,i6))') io,jo,&
                                   kernel_orthog(map_frag_and_env(io,1),map_frag_and_env(jo,1),1),ia,ja
                           else
                              write(unitm) io,jo,kernel_orthog(map_frag_and_env(io,1),map_frag_and_env(jo,1),1),ia,ja
                           end if
                        end do
                     end do
                  !end do
                  call f_close(unitm)
               end if
               call f_free(map_frag_and_env)

               ! deallocate/nullify here to be safe
               nullify(ref_frags(ifrag_ref)%astruct_env%atomnames)
               call deallocate_atomic_structure(ref_frags(ifrag_ref)%astruct_env)
               call nullify_atomic_structure(ref_frags(ifrag_ref)%astruct_env)
            end if


            !pretty sure these should be moved out of fragment loop...
            call f_free_ptr(linmat%kernel_%matrix)
            if (calc_sks) then
               call deallocate_matrices(ovrlp_half_(1))
               call f_free_ptr(kernel_orthog)
            end if

            nullify(kernel_orthog)

            call cpu_time(tr1)
            call system_clock(ncount2,ncount_rate,ncount_max)
            tel=dble(ncount2-ncount1)/dble(ncount_rate)
            if (iproc == 0) then
               call yaml_sequence_open('Write Waves Time')
               call yaml_sequence(advance='no')
               call yaml_mapping_open(flow=.true.)
               call yaml_map('Process',iproc)
               call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
               call yaml_mapping_close()
               call yaml_sequence_close()
            end if
    
            call f_free(frag_map)

            isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
            isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
         end do loop_ifrag
      end if
    
      call f_free(fragment_written)
    
    
    END SUBROUTINE writemywaves_linear_fragments


    !might not be the best place for it - used above and in restart
    function dist_and_shift(periodic,alat,A,B)
      implicit none
      real(kind=8), intent(in) :: A
      real(kind=8), intent(inout) :: B
      real(kind=8) :: dist_and_shift
      real(kind=8), intent(in) :: alat
      logical, intent(in) :: periodic
      !local variables
      integer :: i
      real(kind=8) :: shift

      !shift the B vector on its periodic image
      dist_and_shift = A - B
      if (periodic) then
         !periodic image, if distance is bigger than half of the box
         shift = 0.0d0
         if (dist_and_shift > 0.5d0*alat) then
           shift = alat
         else if (dist_and_shift < -0.5d0*alat) then
           shift = -alat
         end if
         !shift = real(floor(dist_and_shift+0.5d0),kind=8)
         dist_and_shift = dist_and_shift - shift
         B = B + shift
      end if

    end function dist_and_shift



   !> finds environment atoms and fills ref_frag%astruct_env accordingly
   !! also returns mapping array for fragment and environment
   subroutine find_neighbours(num_neighbours,at,rxyz,orbs,ref_frag,frag_map,ntmb_frag_and_env,map_frag_and_env)
      use module_fragments
      use module_types
      use module_atoms, only: deallocate_atomic_structure, nullify_atomic_structure
      use module_base
      implicit none
      integer, intent(in) :: num_neighbours
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(orbitals_data), intent(in) :: orbs
      type(system_fragment), intent(inout) :: ref_frag
      integer, dimension(ref_frag%fbasis%forbs%norb,3), intent(in) :: frag_map
      integer, intent(out) :: ntmb_frag_and_env
      integer, dimension(orbs%norb,3), intent(out) :: map_frag_and_env

      !local variables
      integer :: iatnf, iat, ityp, iatf, iatt, num_before, num_after, jat, num_neighbours_tot, iorb
      integer, allocatable, dimension(:) :: map_not_frag, ipiv
      integer, allocatable, dimension(:) :: num_neighbours_type, atype_not_frag

      real(kind=8), parameter :: tol=1.d-3
      real(kind=8), dimension(3) :: frag_centre
      real(kind=8), allocatable, dimension(:) :: dist
      real(kind=8), allocatable, dimension(:,:) :: rxyz_frag, rxyz_not_frag, rxyz_frag_and_env

      logical :: on_frag, perx, pery, perz


      rxyz_frag = f_malloc((/ 3,ref_frag%astruct_frg%nat /),id='rxyz_frag') 
      rxyz_not_frag = f_malloc((/ 3,at%astruct%nat-ref_frag%astruct_frg%nat /),id='rxyz_not_frag')
      map_not_frag = f_malloc(at%astruct%nat-ref_frag%astruct_frg%nat,id='map_not_frag')

      ! we actually want num_neighbours of each atom type (if possible)
      ! figure out how many atoms that actually is by counting how many of each type are not in fragment
      atype_not_frag=f_malloc0(at%astruct%ntypes,id='atype_not_frag')
      num_neighbours_type=f_malloc(at%astruct%ntypes,id='num_neighbours_type')
             
      iatnf=0
      do iat=1,at%astruct%nat
         on_frag=.false.
         do iatf=1,ref_frag%astruct_frg%nat
            if (frag_map(iatf,3)==iat) then
               on_frag=.true.
               rxyz_frag(:,iatf)=rxyz(:,iat)
               exit
            end if
         end do
         if (.not. on_frag) then
            iatnf=iatnf+1
            rxyz_not_frag(:,iatnf) = rxyz(:,iat)
            map_not_frag(iatnf)=iat
            atype_not_frag(at%astruct%iatype(iat)) = atype_not_frag(at%astruct%iatype(iat)) + 1
         end if
      end do

      num_neighbours_tot=0
      do ityp=1,at%astruct%ntypes
         num_neighbours_type(ityp) = min(atype_not_frag(ityp),num_neighbours)
         num_neighbours_tot = num_neighbours_tot + num_neighbours_type(ityp)
         !write(*,'(a,6(i3,2x))') 'type',ifrag,ifrag_ref,ityp,at%astruct%ntypes,num_neighbours_type(ityp),num_neighbours_tot
      end do
      call f_free(atype_not_frag)

      ! find distances from centre of fragment - not necessarily the best approach but fine if fragment has 1 atom (most likely case)
      frag_centre=frag_center(ref_frag%astruct_frg%nat,rxyz_frag)

      dist = f_malloc(at%astruct%nat-ref_frag%astruct_frg%nat,id='dist')
      ipiv = f_malloc(at%astruct%nat-ref_frag%astruct_frg%nat,id='ipiv')
      ! find distances from this atom BEFORE shifting
      perx=(at%astruct%geocode /= 'F')
      pery=(at%astruct%geocode == 'P')
      perz=(at%astruct%geocode /= 'F')

      ! if coordinates wrap around (in periodic), take this into account
      ! assume that the fragment and environment are written in free bc
      ! think about other periodic cases that might need fixing...
      ! no need to shift wrt centre as this is done when reading in
      do iat=1,at%astruct%nat-ref_frag%astruct_frg%nat
         dist(iat) = dist_and_shift(perx,at%astruct%cell_dim(1),frag_centre(1),rxyz_not_frag(1,iat))**2
         dist(iat) = dist(iat) + dist_and_shift(pery,at%astruct%cell_dim(2),frag_centre(2),rxyz_not_frag(2,iat))**2
         dist(iat) = dist(iat) + dist_and_shift(perz,at%astruct%cell_dim(3),frag_centre(3),rxyz_not_frag(3,iat))**2

         dist(iat) = -dsqrt(dist(iat))
!!$         write(*,'(A,2(I3,2x),F12.6,3x,2(3(F12.6,1x),2x))') 'ifrag,iat,dist',ifrag,iat,dist(iat),&
!!$              at%astruct%cell_dim(:),rxyz_new_all(:,iat)

         !rxyz_not_frag(:,iat) = rxyz_new_all(:,iat)-frag_trans_frag(ifrag)%rot_center_new
      end do     

      ! sort atoms into neighbour order
      call sort_positions(at%astruct%nat-ref_frag%astruct_frg%nat,dist,ipiv)

      rxyz_frag_and_env = f_malloc((/ 3,ref_frag%astruct_frg%nat+num_neighbours_tot /),id='rxyz_frag_and_env')

      ! take fragment and closest neighbours (assume that environment atoms were originally the closest)
      do iat=1,ref_frag%astruct_frg%nat
         rxyz_frag_and_env(:,iat) = rxyz_frag(:,iat)
         map_frag_and_env(iat,3) = frag_map(iat,3)
      end do
      call f_free(rxyz_frag)

      iatf=0
      do ityp=1,at%astruct%ntypes
         iatt=0
         if (num_neighbours_type(ityp)==0) cycle
         do iat=1,at%astruct%nat-ref_frag%astruct_frg%nat
            if (at%astruct%iatype(map_not_frag(ipiv(iat)))/=ityp) cycle
            iatf=iatf+1
            iatt=iatt+1
            rxyz_frag_and_env(:,iatf+ref_frag%astruct_frg%nat) = rxyz_not_frag(:,ipiv(iat))
            map_frag_and_env(iatf+ref_frag%astruct_frg%nat,3) = map_not_frag(ipiv(iat))
            !print*,'iatt',ityp,iatt,num_neighbours_type(ityp),iatf,iat
            if (iatt==num_neighbours_type(ityp)) exit
         end do

         !never cut off a shell - either include all at that distance or none (depending on if at least at minimum?)
         ! - check distances of next point to see...
         num_after=0
         do jat=iat+1,at%astruct%nat-ref_frag%astruct_frg%nat
            if (at%astruct%iatype(map_not_frag(ipiv(jat)))/=ityp) cycle
            if (abs(dist(ipiv(jat))-dist(ipiv(iat))) < tol) then
               num_after=num_after+1
            else
               exit
            end if
         end do
         if (num_after==0) cycle
         !also check before
         num_before=0
         do jat=iat-1,1,-1
            if (at%astruct%iatype(map_not_frag(ipiv(jat)))/=ityp) cycle
            if (abs(dist(ipiv(jat))-dist(ipiv(iat))) < tol) then
               num_before=num_before+1
            else
               exit
            end if
         end do
         !get rid of them (assuming we will still have at least one neighbour of this type left)
         if (num_neighbours_type(ityp)>num_before+1) then
            num_neighbours_type(ityp)=num_neighbours_type(ityp)-num_before-1
            num_neighbours_tot=num_neighbours_tot-num_before-1
            iatf=iatf-num_before-1
         !add neighbours until shell is complete
         else
            num_neighbours_type(ityp)=num_neighbours_type(ityp)+num_after
            num_neighbours_tot=num_neighbours_tot+num_after
            do jat=iat+1,at%astruct%nat-ref_frag%astruct_frg%nat
               if (at%astruct%iatype(map_not_frag(ipiv(jat)))/=ityp) cycle
               iatf=iatf+1
               iatt=iatt+1
               rxyz_frag_and_env(:,iatf+ref_frag%astruct_frg%nat) = rxyz_not_frag(:,ipiv(jat))
               map_frag_and_env(iatf+ref_frag%astruct_frg%nat,3) = map_not_frag(ipiv(jat))
               if (iatt==num_neighbours_type(ityp)) exit 
            end do
         end if
      end do
      if (iatf/=num_neighbours_tot) stop 'Error num_neighbours_tot/=iatf in write_linear_fragments'
      call f_free(num_neighbours_type)
      call f_free(rxyz_not_frag)
      call f_free(map_not_frag)
      call f_free(dist)
      call f_free(ipiv)

      ! fill in map for tmbs
      ntmb_frag_and_env = 0
      do iat=1,ref_frag%astruct_frg%nat+num_neighbours_tot
         do iorb=1,orbs%norb
            if (orbs%onwhichatom(iorb)==map_frag_and_env(iat,3)) then
               ntmb_frag_and_env = ntmb_frag_and_env+1
               map_frag_and_env(ntmb_frag_and_env,1)=iorb
               map_frag_and_env(ntmb_frag_and_env,2)=iat
            end if
         end do
      end do

      ! put rxyz_env into astruct_env - first delete pre-existing structure
      if (ref_frag%astruct_env%nat/=0) then
         call deallocate_atomic_structure(ref_frag%astruct_env)
         call nullify_atomic_structure(ref_frag%astruct_env)
      end if

      ! copy some stuff from astruct_frg
      ref_frag%astruct_env%inputfile_format = ref_frag%astruct_frg%inputfile_format
      ref_frag%astruct_env%units = ref_frag%astruct_frg%units
      ref_frag%astruct_env%geocode = ref_frag%astruct_frg%geocode

      ref_frag%astruct_env%nat = ref_frag%astruct_frg%nat+num_neighbours_tot
      ref_frag%astruct_env%rxyz = f_malloc_ptr((/3,ref_frag%astruct_env%nat/),id='ref_frag%astruct_env%rxyz')
      call vcopy(3*ref_frag%astruct_env%nat,rxyz_frag_and_env(1,1),1,ref_frag%astruct_env%rxyz(1,1),1)
      call f_free(rxyz_frag_and_env)

      ! now deal with atom types
      ref_frag%astruct_env%ntypes = at%astruct%ntypes
      ref_frag%astruct_env%atomnames => at%astruct%atomnames

      ! can't just point to full version due to atom reordering
      ref_frag%astruct_env%iatype = f_malloc_ptr(ref_frag%astruct_env%nat,id='ref_frag%astruct_env%iatype')

      ! polarization etc is irrelevant
      ref_frag%astruct_env%input_polarization = f_malloc0_ptr(ref_frag%astruct_env%nat,&
           id='ref_frag%astruct_env%input_polarization')
      ref_frag%astruct_env%ifrztyp = f_malloc0_ptr(ref_frag%astruct_env%nat,id='ref_frag%astruct_env%ifrztyp')

      do iat=1,ref_frag%astruct_env%nat
         ref_frag%astruct_env%iatype(iat) = at%astruct%iatype(map_frag_and_env(iat,3))
         ref_frag%astruct_env%input_polarization(iat) = 100
         !if (iproc==0) print*,iat,trim(ref_frag%astruct_env%atomnames(ref_frag%astruct_env%iatype(iat)))
      end do

   end subroutine find_neighbours


    subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,ns1,ns2,ns3,hx,hy,hz,locregCenter,&
         locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  &
         nseg_f,nvctr_f,keyg_f,keyv_f, &
         psi_c,psi_f,eval,onwhichatom)
      use module_base
      use yaml_output
      implicit none
      logical, intent(in) :: useFormattedOutput
      integer, intent(in) :: unitwf,iorb,n1,n2,n3,ns1,ns2,ns3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
      real(gp), intent(in) :: hx,hy,hz,locrad,confPotprefac
      real(wp), intent(in) :: eval
      integer, dimension(nseg_c), intent(in) :: keyv_c
      integer, dimension(nseg_f), intent(in) :: keyv_f
      integer, dimension(2,nseg_c), intent(in) :: keyg_c
      integer, dimension(2,nseg_f), intent(in) :: keyg_f
      real(wp), dimension(nvctr_c), intent(in) :: psi_c
      real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
      real(gp), dimension(3,nat), intent(in) :: rxyz
      real(gp), dimension(3), intent(in) :: locregCenter
      integer, intent(in) :: onwhichatom
      !local variables
      integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j,np,n1p1
      real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
    
      if (useFormattedOutput) then
         write(unitwf,*) iorb,eval
         write(unitwf,*) hx,hy,hz
         write(unitwf,*) n1,n2,n3
         write(unitwf,*) ns1,ns2,ns3
         write(unitwf,*) locregCenter(1),locregCenter(2),locregCenter(3),onwhichatom,locrad,&
              confPotOrder,confPotprefac
         write(unitwf,*) nat
         do iat=1,nat
            write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
         enddo
         write(unitwf,*) nvctr_c, nvctr_f
      else
         write(unitwf) iorb,eval
         write(unitwf) hx,hy,hz
         write(unitwf) n1,n2,n3
         write(unitwf) ns1,ns2,ns3
         write(unitwf) locregCenter(1),locregCenter(2),locregCenter(3),onwhichatom,locrad,&
              confPotOrder,confPotprefac
         write(unitwf) nat
         do iat=1,nat
            write(unitwf) (rxyz(j,iat),j=1,3)
         enddo
         write(unitwf) nvctr_c, nvctr_f
      end if
    
      n1p1=n1+1
      np=n1p1*(n2+1)
    
      ! coarse part
      do iseg=1,nseg_c
         jj=keyv_c(iseg)
         j0=keyg_c(1,iseg)
         j1=keyg_c(2,iseg)
         ii=j0-1
         i3=ii/np
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         do i=i0,i1
            tt=psi_c(i-i0+jj)
            if (useFormattedOutput) then
               write(unitwf,'(3(i4),1x,e19.12)') i,i2,i3,tt
            else
               write(unitwf) i,i2,i3,tt
            end if
         enddo
      enddo
    
      ! fine part
      do iseg=1,nseg_f
         jj=keyv_f(iseg)
         j0=keyg_f(1,iseg)
         j1=keyg_f(2,iseg)
         ii=j0-1
         i3=ii/np
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         do i=i0,i1
            t1=psi_f(1,i-i0+jj)
            t2=psi_f(2,i-i0+jj)
            t3=psi_f(3,i-i0+jj)
            t4=psi_f(4,i-i0+jj)
            t5=psi_f(5,i-i0+jj)
            t6=psi_f(6,i-i0+jj)
            t7=psi_f(7,i-i0+jj)
            if (useFormattedOutput) then
               write(unitwf,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
               write(unitwf) i,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
         enddo
      enddo
    
      if (verbose >= 2 .and. bigdft_mpi%iproc==0) call yaml_map('Wavefunction written No.',iorb)
      !if (verbose >= 2) write(*,'(1x,i0,a)') iorb,'th wavefunction written'
    
    END SUBROUTINE writeonewave_linear


    subroutine read_linear_matrix_dense(iunit, ntmb, nat, matrix, rxyz, on_which_atom)
      use module_base
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iunit, ntmb, nat
      real(kind=8),dimension(ntmb,ntmb),intent(out) :: matrix
      real(kind=8),dimension(3,nat),intent(out),optional :: rxyz
      integer,dimension(ntmb),intent(out),optional :: on_which_atom
    
      ! Local variables
      integer :: itmb, jtmb, ii, jj, iat, ntmb_check, nat_check
      logical :: read_rxyz, read_on_which_atom
      real(kind=8),dimension(3) :: dummy
      character(len=128) :: dummy_char
    
      read_on_which_atom = present(on_which_atom)
      read_rxyz = present(rxyz)

      read(iunit,*) dummy_char, ntmb_check, nat_check
      if (ntmb/=ntmb_check) then
          call f_err_throw('number of basis function specified ('//trim(yaml_toa(ntmb,fmt='(i0)'))//&
              &') does not agree with the number indicated in the file ('//trim(yaml_toa(ntmb_check,fmt='(i0)'))//')',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (nat/=nat_check) then
          call f_err_throw('number of atoms specified ('//trim(yaml_toa(nat,fmt='(i0)'))//&
              &') does not agree with the number indicated in the file ('//trim(yaml_toa(nat_check))//')',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end if
      do iat=1,nat
          if (read_rxyz) then
              read(iunit,*) dummy_char, rxyz(1:3,iat)
          else
              read(iunit,*) dummy_char, dummy(1:3)
          end if
      end do
    
      do itmb=1,ntmb
          do jtmb=1,ntmb
              if(read_on_which_atom .and. jtmb==1) then
                  read(iunit,*) ii, jj, matrix(ii,jj), on_which_atom(itmb)
              else
                  read(iunit,*) ii, jj, matrix(ii,jj)
              end if
              if (ii/=itmb) call f_err_throw('ii/=itmb',err_name='BIGDFT_RUNTIME_ERROR')
              if (jj/=jtmb) call f_err_throw('jj/=jtmb',err_name='BIGDFT_RUNTIME_ERROR')
          end do
      end do
    
    end subroutine read_linear_matrix_dense


    subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, &
           & lstat, error, nat, rxyz_old)
        use module_base
        use module_types
        !use internal_io
        implicit none
        integer, intent(in) :: unitwf
        logical, intent(in) :: formatted
        integer, intent(out) :: norb_old, ntmb_old
        logical, intent(out) :: lstat
        character(len =256), intent(out) :: error
        ! Optional arguments
        integer, intent(in), optional :: nat
        real(gp), dimension(:,:), intent(out), optional :: rxyz_old
    
        integer :: i, iat, i_stat, nat_
        real(gp) :: rxyz(3)
    
        lstat = .false.
        write(error, "(A)") "cannot read coeff description."
        if (formatted) then
           read(unitwf,*,iostat=i_stat) ntmb_old, norb_old
           if (i_stat /= 0) return
           !write(*,*) 'reading ',nat,' atomic positions'
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
              enddo
           else
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,*,iostat=i_stat)
                 if (i_stat /= 0) return
              enddo
           end if
           !read(unitwf,*,iostat=i_stat) i, iat
           !if (i_stat /= 0) return
        else
           read(unitwf,iostat=i_stat) ntmb_old, norb_old
           if (i_stat /= 0) return
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
              enddo
           else
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,iostat=i_stat) rxyz
                 if (i_stat /= 0) return
              enddo
           end if
           !read(unitwf,iostat=i_stat) i, iat
           !if (i_stat /= 0) return
        end if
        lstat = .true.
    END SUBROUTINE io_read_descr_coeff


    subroutine read_coeff_minbasis(unitwf,useFormattedInput,iproc,ntmb,norb_old,nfvctr,coeff,eval,nat,rxyz_old)
      use module_base
      use module_types
      !use internal_io
      use yaml_output
      implicit none
      logical, intent(in) :: useFormattedInput
      integer, intent(in) :: unitwf,iproc,ntmb,nfvctr
      integer, intent(out) :: norb_old
      real(wp), dimension(nfvctr,ntmb), intent(out) :: coeff
      real(wp), dimension(ntmb), intent(out) :: eval
      integer, optional, intent(in) :: nat
      real(gp), dimension(:,:), optional, intent(out) :: rxyz_old
    
      !local variables
      character(len = 256) :: error
      logical :: lstat
      integer :: i_stat
      integer :: ntmb_old, i1, i2,i,j,iorb,iorb_old
      real(wp) :: tt
    
      call io_read_descr_coeff(unitwf, useFormattedInput, norb_old, ntmb_old, &
           & lstat, error, nat, rxyz_old)
      if (.not. lstat) call io_error(trim(error))
    
      if (ntmb_old /= ntmb) then
         if (iproc == 0) write(error,"(A)") 'error in read coeffs, ntmb_old/=ntmb'
         call io_error(trim(error))
      end if
    
      ! read the eigenvalues
      if (useFormattedInput) then
         do iorb=1,ntmb
            read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
            if (iorb_old /= iorb) stop 'read_coeff_minbasis'
         enddo
      else 
         do iorb=1,ntmb
            read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
            if (iorb_old /= iorb) stop 'read_coeff_minbasis'
         enddo
         if (i_stat /= 0) stop 'Problem reading the eigenvalues'
      end if
    
      !if (iproc == 0) write(*,*) 'coefficients need NO reformatting'
    
      ! Now read the coefficients
      do i = 1, ntmb
         do j = 1,nfvctr
            if (useFormattedInput) then
               read(unitwf,*,iostat=i_stat) i1,i2,tt
            else
               read(unitwf,iostat=i_stat) i1,i2,tt
            end if
            if (i_stat /= 0) stop 'Problem reading the coefficients'
            coeff(j,i) = tt  
         end do
      end do
    
      ! rescale so first significant element is +ve
      do i = 1, ntmb
         do j = 1,nfvctr
            if (abs(coeff(j,i))>1.0e-1) then
               if (coeff(j,i)<0.0_gp) call dscal(ntmb,-1.0_gp,coeff(1,i),1)
               exit
            end if
         end do
         !if (j==ntmb+1) print*,'Error finding significant coefficient, coefficients not scaled to have +ve first element'
      end do
    
    END SUBROUTINE read_coeff_minbasis

    !1subroutine read_kernel(unitwf,useFormattedInput,iproc,ntmb,norb_old,kernel,nat,rxyz_old)
    !1  use module_base
    !1  use module_types
    !1  !use internal_io
    !1  use yaml_output
    !1  implicit none
    !1  logical, intent(in) :: useFormattedInput
    !1  integer, intent(in) :: unitwf,iproc,ntmb
    !1  integer, intent(out) :: norb_old
    !1  real(wp), dimension(ntmb,ntmb), intent(out) :: kernel
    !1  !real(wp), dimension(ntmb), intent(out) :: eval !eliminate use of eval?
    !1  integer, optional, intent(in) :: nat
    !1  real(gp), dimension(:,:), optional, intent(out) :: rxyz_old
    !1
    !1  !local variables
    !1  character(len = 256) :: error
    !1  logical :: lstat
    !1  integer :: i_stat
    !1  integer :: ntmb_old, i1, i2,i,j,iorb,iorb_old
    !1  real(wp) :: tt
    !1
    !1  call io_read_descr_kernel(unitwf, useFormattedInput, norb_old, ntmb_old, &
    !1       & lstat, error, nat, rxyz_old)
    !1  if (.not. lstat) call io_error(trim(error))
    !1
    !1  if (ntmb_old /= ntmb) then
    !1     if (iproc == 0) write(error,"(A)") 'error in read coeffs, ntmb_old/=ntmb'
    !1     call io_error(trim(error))
    !1  end if
    !1
    !1  ! read the eigenvalues
    !1  if (useFormattedInput) then
    !1     do iorb=1,ntmb
    !1        read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
    !1        if (iorb_old /= iorb) stop 'read_coeff_minbasis'
    !1     enddo
    !1  else 
    !1     do iorb=1,ntmb
    !1        read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
    !1        if (iorb_old /= iorb) stop 'read_coeff_minbasis'
    !1     enddo
    !1     if (i_stat /= 0) stop 'Problem reading the eigenvalues'
    !1  end if
    !1
    !1  !if (iproc == 0) write(*,*) 'coefficients need NO reformatting'
    !1
    !1  ! Now read the coefficients
    !1  do i = 1, ntmb
    !1     do j = 1,nfvctr
    !1        if (useFormattedInput) then
    !1           read(unitwf,*,iostat=i_stat) i1,i2,tt
    !1        else
    !1           read(unitwf,iostat=i_stat) i1,i2,tt
    !1        end if
    !1        if (i_stat /= 0) stop 'Problem reading the coefficients'
    !1        coeff(j,i) = tt  
    !1     end do
    !1  end do
    !1
    !1  ! rescale so first significant element is +ve
    !1  do i = 1, ntmb
    !1     do j = 1,nfvctr
    !1        if (abs(coeff(j,i))>1.0e-1) then
    !1           if (coeff(j,i)<0.0_gp) call dscal(ntmb,-1.0_gp,coeff(1,i),1)
    !1           exit
    !1        end if
    !1     end do
    !1     !if (j==ntmb+1) print*,'Error finding significant coefficient, coefficients not scaled to have +ve first element'
    !1  end do
    !1
    !1END SUBROUTINE read_kernel


    subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n_old1, n_old2, n_old3, &
           & ns_old1, ns_old2, ns_old3, hgrids_old, lstat, error, onwhichatom, locrad, locregCenter, &
           & confPotOrder, confPotprefac, nvctr_c_old, nvctr_f_old, nat, rxyz_old)
        use module_base
        use module_types
        !use internal_io
        use yaml_output
        implicit none
    
        integer, intent(in) :: unitwf
        logical, intent(in) :: formatted
        integer, intent(out) :: iorb_old
        integer, intent(out) :: n_old1, n_old2, n_old3, ns_old1, ns_old2, ns_old3
        real(gp), dimension(3), intent(out) :: hgrids_old
        logical, intent(out) :: lstat
        real(wp), intent(out) :: eval
        real(gp), intent(out) :: locrad
        real(gp), dimension(3), intent(out) :: locregCenter
        character(len =256), intent(out) :: error
        integer, intent(out) :: onwhichatom
        integer, intent(out) :: confPotOrder
        real(gp), intent(out) :: confPotprefac
        ! Optional arguments
        integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
        integer, intent(in), optional :: nat
        real(gp), dimension(:,:), intent(out), optional :: rxyz_old
    
        integer :: i, iat, i_stat, nat_
        real(gp) :: rxyz(3)
    
        lstat = .false.
        write(error, "(A)") "cannot read psi description."
        if (formatted) then
           read(unitwf,*,iostat=i_stat) iorb_old,eval
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) hgrids_old(1),hgrids_old(2),hgrids_old(3)
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) n_old1,n_old2,n_old3
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) ns_old1,ns_old2,ns_old3
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
                locrad,confPotOrder, confPotprefac
           if (i_stat /= 0) return
           !call yaml_map('Reading atomic positions',nat)
           !write(*,*) 'reading ',nat,' atomic positions' !*
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
    
              enddo
           else
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,*,iostat=i_stat)
                 if (i_stat /= 0) return
              enddo
           end if
           if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
              read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
              if (i_stat /= 0) return
           else
              read(unitwf,*,iostat=i_stat) i, iat
              if (i_stat /= 0) return
           end if
        else
           read(unitwf,iostat=i_stat) iorb_old,eval
           if (i_stat /= 0) return
    
           read(unitwf,iostat=i_stat) hgrids_old(1),hgrids_old(2),hgrids_old(3)
           if (i_stat /= 0) return
           read(unitwf,iostat=i_stat) n_old1,n_old2,n_old3
           if (i_stat /= 0) return
           read(unitwf,iostat=i_stat) ns_old1,ns_old2,ns_old3
           if (i_stat /= 0) return
           read(unitwf,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
                locrad,confPotOrder, confPotprefac
           if (i_stat /= 0) return
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
              enddo
           else
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,iostat=i_stat) rxyz
                 if (i_stat /= 0) return
              enddo
           end if
           if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
              read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
              if (i_stat /= 0) return
           else
              read(unitwf,iostat=i_stat) i, iat
              if (i_stat /= 0) return
           end if
        end if
        lstat = .true.
    
    END SUBROUTINE io_read_descr_linear


    subroutine io_error(error)
      use module_base, only: bigdft_mpi
  
      implicit none
  
      character(len = *), intent(in) :: error
      integer :: ierr
  
      call io_warning(error)
      call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
    END SUBROUTINE io_error
  
  
    subroutine io_warning(error)
      use module_defs
  
      implicit none
  
      character(len = *), intent(in) :: error
  
      write(0,"(2A)") "WARNING! ", trim(error)
    END SUBROUTINE io_warning
  
  
    !> Read the input/output descriptors (for a wavefunction for instance)
    subroutine io_read_descr(unitwf, formatted, iorb_old, eval, n1_old, n2_old, n3_old, &
         & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
      use module_base
      use module_types
  
      implicit none
  
      integer, intent(in) :: unitwf
      logical, intent(in) :: formatted
      integer, intent(out) :: iorb_old
      integer, intent(out) :: n1_old, n2_old, n3_old
      real(gp), intent(out) :: hx_old, hy_old, hz_old
      logical, intent(out) :: lstat
      real(wp), intent(out) :: eval
      character(len =256), intent(out) :: error
      ! Optional arguments
      integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
      integer, intent(in), optional :: nat
      real(gp), dimension(:,:), intent(out), optional :: rxyz_old
  
      integer :: i, iat, i_stat, nat_
      real(gp) :: rxyz(3)
  
      lstat = .false.
      write(error, "(A)") "cannot read psi description."
      if (formatted) then
         read(unitwf,*,iostat=i_stat) iorb_old,eval
         if (i_stat /= 0) return
         read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
         if (i_stat /= 0) return
         read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
         if (i_stat /= 0) return
         !write(*,*) 'reading ',nat,' atomic positions'
         if (present(nat) .And. present(rxyz_old)) then
            read(unitwf,*,iostat=i_stat) nat_
            if (i_stat /= 0) return
            ! Sanity check
            if (size(rxyz_old, 2) /= nat) call io_error("Mismatch in coordinate array size.")
            if (nat_ /= nat) call io_error("Mismatch in coordinate array size.")
            do iat=1,nat
               read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
               if (i_stat /= 0) return
            enddo
         else
            read(unitwf,*,iostat=i_stat) nat_
            if (i_stat /= 0) return
            do iat=1,nat_
               read(unitwf,*,iostat=i_stat)
               if (i_stat /= 0) return
            enddo
         end if
         if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
            read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
            if (i_stat /= 0) return
         else
            read(unitwf,*,iostat=i_stat) i, iat
            if (i_stat /= 0) return
         end if
      else
         read(unitwf,iostat=i_stat) iorb_old,eval
         if (i_stat /= 0) return
         read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
         if (i_stat /= 0) return
         read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
         if (i_stat /= 0) return
         if (present(nat) .And. present(rxyz_old)) then
            read(unitwf,iostat=i_stat) nat_
            if (i_stat /= 0) return
            ! Sanity check
            if (size(rxyz_old, 2) /= nat) call io_error("Mismatch in coordinate array size.")
            if (nat_ /= nat) call io_error("Mismatch in coordinate array size.")
            do iat=1,nat
               read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
               if (i_stat /= 0) return
            enddo
         else
            read(unitwf,iostat=i_stat) nat_
            if (i_stat /= 0) return
            do iat=1,nat_
               read(unitwf,iostat=i_stat) rxyz
               if (i_stat /= 0) return
            enddo
         end if
         if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
            read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
            if (i_stat /= 0) return
         else
            read(unitwf,iostat=i_stat) i, iat
            if (i_stat /= 0) return
         end if
      end if
      lstat = .true.
    END SUBROUTINE io_read_descr
  
  
    subroutine io_gcoordToLocreg(n1, n2, n3, nvctr_c, nvctr_f, gcoord_c, gcoord_f, lr)
      use module_base
      use locregs
  
      implicit none
      !Arguments
      integer, intent(in) :: n1, n2, n3, nvctr_c, nvctr_f
      integer, dimension(3, nvctr_c), intent(in) :: gcoord_c
      integer, dimension(3, nvctr_f), intent(in) :: gcoord_f
      type(locreg_descriptors), intent(out) :: lr
      !Local variables
      character(len = *), parameter :: subname = "io_gcoordToLocreg"
      integer :: i
      logical, dimension(:,:,:), allocatable :: logrid_c, logrid_f
  
      call f_routine(id=subname)
  
      call nullify_locreg_descriptors(lr)
  
      lr%geocode = "P"
      lr%hybrid_on = .false.
  
      lr%ns1 = 0
      lr%ns2 = 0
      lr%ns3 = 0
  
      lr%d%n1 = n1
      lr%d%n2 = n2
      lr%d%n3 = n3
  
      lr%d%n1i = 2 * n1 + 2
      lr%d%n2i = 2 * n2 + 2
      lr%d%n3i = 2 * n3 + 2
  
      logrid_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_c')
      logrid_f = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_f')
  
      lr%d%nfl1 = n1
      lr%d%nfl2 = n2
      lr%d%nfl3 = n3
      lr%d%nfu1 = 0
      lr%d%nfu2 = 0
      lr%d%nfu3 = 0
  
      logrid_c(:,:,:) = .false.
      do i = 1, nvctr_c, 1
         logrid_c(gcoord_c(1, i), gcoord_c(2, i), gcoord_c(3, i)) = .true.
      end do
      logrid_f(:,:,:) = .false.
      do i = 1, nvctr_f, 1
         logrid_f(gcoord_f(1, i), gcoord_f(2, i), gcoord_f(3, i)) = .true.
         lr%d%nfl1 = min(lr%d%nfl1, gcoord_f(1, i))
         lr%d%nfl2 = min(lr%d%nfl2, gcoord_f(2, i))
         lr%d%nfl3 = min(lr%d%nfl3, gcoord_f(3, i))
         lr%d%nfu1 = max(lr%d%nfu1, gcoord_f(1, i))
         lr%d%nfu2 = max(lr%d%nfu2, gcoord_f(2, i))
         lr%d%nfu3 = max(lr%d%nfu3, gcoord_f(3, i))
      end do
  
      !correct the values of the delimiter if there are no wavelets
      if (lr%d%nfl1 == n1 .and. lr%d%nfu1 == 0) then
         lr%d%nfl1 = n1 / 2
         lr%d%nfu1 = n1 / 2
      end if
      if (lr%d%nfl2 == n2 .and. lr%d%nfu2 == 0) then
         lr%d%nfl2 = n2 / 2
         lr%d%nfu2 = n2 / 2
      end if
      if (lr%d%nfl3 == n3 .and. lr%d%nfu3 == 0) then
         lr%d%nfl3 = n3 / 2
         lr%d%nfu3 = n3 / 2
      end if
  
      call wfd_from_grids(logrid_c, logrid_f, .true., lr)
  
      call f_free(logrid_c)
      call f_free(logrid_f)
  
      call f_release_routine()
  
    END SUBROUTINE io_gcoordToLocreg
  
    subroutine read_psi_compress(unitwf, formatted, nvctr_c, nvctr_f, psi, lstat, error, gcoord_c, gcoord_f)
      use module_base
      use module_types
  
      implicit none
  
      integer, intent(in) :: unitwf, nvctr_c, nvctr_f
      logical, intent(in) :: formatted
      real(wp), dimension(nvctr_c+7*nvctr_f), intent(out) :: psi
      logical, intent(out) :: lstat
      character(len =256), intent(out) :: error
      integer, dimension(3, nvctr_c), optional, intent(out) :: gcoord_c
      integer, dimension(3, nvctr_f), optional, intent(out) :: gcoord_f
  
      integer :: j, i1, i2, i3, i_stat
      real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
  
      lstat = .false.
      write(error, "(A)") "cannot read psi values."
      if (present(gcoord_c)) then
         do j=1,nvctr_c
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,tt
            end if
            if (i_stat /= 0) return
            psi(j)=tt
            gcoord_c(:, j) = (/ i1, i2, i3 /)
         enddo
      else
         do j=1,nvctr_c
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,tt
            end if
            if (i_stat /= 0) return
            psi(j)=tt
         enddo
      end if
      if (present(gcoord_f)) then
         do j=1,7*nvctr_f-6,7
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
            if (i_stat /= 0) return
            psi(nvctr_c+j+0)=t1
            psi(nvctr_c+j+1)=t2
            psi(nvctr_c+j+2)=t3
            psi(nvctr_c+j+3)=t4
            psi(nvctr_c+j+4)=t5
            psi(nvctr_c+j+5)=t6
            psi(nvctr_c+j+6)=t7
            gcoord_f(:, (j - 1) / 7 + 1) = (/ i1, i2, i3 /)
         enddo
      else
         do j=1,7*nvctr_f-6,7
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
            if (i_stat /= 0) return
            psi(nvctr_c+j+0)=t1
            psi(nvctr_c+j+1)=t2
            psi(nvctr_c+j+2)=t3
            psi(nvctr_c+j+3)=t4
            psi(nvctr_c+j+4)=t5
            psi(nvctr_c+j+5)=t6
            psi(nvctr_c+j+6)=t7
         enddo
      end if
      lstat = .true.
    END SUBROUTINE read_psi_compress
  
  
    subroutine read_psig(unitwf, formatted, nvctr_c, nvctr_f, n1, n2, n3, psig, lstat, error)
      use module_base
      use module_types
  
      implicit none
  
      integer, intent(in) :: unitwf, nvctr_c, nvctr_f, n1, n2, n3
      logical, intent(in) :: formatted
      real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psig
      logical, intent(out) :: lstat
      character(len =256), intent(out) :: error
  
      integer :: i1, i2, i3, i_stat, iel
      real(wp) :: tt, t1, t2, t3, t4, t5, t6, t7
  
      lstat = .false.
      write(error, "(A)") "cannot read psig values."
  
      call f_zero(psig)
      do iel=1,nvctr_c
         if (formatted) then
            read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
         else
            read(unitwf,iostat=i_stat) i1,i2,i3,tt
         end if
         if (i_stat /= 0) return
         psig(i1,1,i2,1,i3,1)=tt
      enddo
      do iel=1,nvctr_f
         if (formatted) then
            read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         else
            read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         end if
         if (i_stat /= 0) return
         psig(i1,2,i2,1,i3,1)=t1
         psig(i1,1,i2,2,i3,1)=t2
         psig(i1,2,i2,2,i3,1)=t3
         psig(i1,1,i2,1,i3,2)=t4
         psig(i1,2,i2,1,i3,2)=t5
         psig(i1,1,i2,2,i3,2)=t6
         psig(i1,2,i2,2,i3,2)=t7
      enddo
      lstat = .true.
    END SUBROUTINE read_psig
  
    subroutine io_open(unitwf, filename, formatted)
      use f_utils, only: f_open_file
      implicit none
      character(len = *), intent(in) :: filename
      logical, intent(in) :: formatted
      integer, intent(out) :: unitwf
  
      integer :: i_stat
  
      ! We open the Fortran file
      unitwf = 99
      call f_open_file(unitwf,file=trim(filename),binary=.not. formatted)
  !!$    if (.not. formatted) then
  !!$       open(unit=unitwf,file=trim(filename),status='unknown',form="unformatted", iostat=i_stat)
  !!$    else
  !!$       open(unit=unitwf,file=trim(filename),status='unknown', iostat=i_stat)
  !!$    end if
  !!$    if (i_stat /= 0) then
  !!$       call io_warning("Cannot open file '" // trim(filename) // "'.")
  !!$       unitwf = -1
  !!$       return
  !!$    end if
    END SUBROUTINE io_open


    !> Write a sparse matrix to disk.
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! in distributed among the matrix taksgroups
    subroutine write_sparse_matrix(nat, ntypes, iatype, rxyz, nzatom, nelpsp, atomnames, smat, mat, filename)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, &
                                   assignment(=), sparsematrix_malloc
      use sparsematrix, only: gather_matrix_from_taskgroups
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: nat, ntypes
      integer,dimension(nat),intent(in) :: iatype
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      integer,dimension(ntypes),intent(in) :: nzatom, nelpsp
      character(len=*),dimension(ntypes),intent(in) :: atomnames
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, iat, jat, ind, ispin, itype
      real(kind=8),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_sparse_matrix')

      matrix_compr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='matrix_compr')
      call gather_matrix_from_taskgroups(bigdft_mpi%iproc, bigdft_mpi%nproc, &
           smat, mat%matrix_compr, matrix_compr)

      if (bigdft_mpi%iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)

          write(iunit,'(i10,2i6,3x,a,a)') nat, ntypes, smat%nspin, smat%geocode, &
              '   # number of atoms, number of atom types, nspin, geocode'
          do itype=1,ntypes
              write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
                  '   # nz, nelpsp, name'
          end do
          do iat=1,nat
              write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          end do
          write(iunit,'(3i12,a)') smat%nfvctr, smat%nseg, smat%nvctr, '   # nfvctr, nseg, nvctr'
          do iseg=1,smat%nseg
              write(iunit,'(5i12,a)') smat%keyv(iseg), smat%keyg(1,1,iseg), smat%keyg(2,1,iseg), &
                  smat%keyg(1,2,iseg), smat%keyg(2,2,iseg), '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
          end do
          ind = 0
          do ispin=1,smat%nspin
              do iseg=1,smat%nseg
                  icol = smat%keyg(1,2,iseg)
                  iat = smat%on_which_atom(icol)
                  do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                      irow = jorb
                      jat = smat%on_which_atom(irow)
                      ind = ind + 1
                      write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
                  end do
              end do
          end do

          call f_close(iunit)

          call f_free(matrix_compr)

      end if

      call f_release_routine()

    end subroutine write_sparse_matrix


    !> Write a dense matrix to disk.
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! in distributed among the matrix taksgroups
    subroutine write_dense_matrix(nat, ntypes, iatype, rxyz, nzatom, nelpsp, atomnames, smat, mat, filename, binary, orbs)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, DENSE_FULL, &
                                   assignment(=), sparsematrix_malloc_ptr
      use sparsematrix, only: uncompress_matrix2
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: nat, ntypes
      integer,dimension(nat),intent(in) :: iatype
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      integer,dimension(ntypes),intent(in) :: nzatom, nelpsp
      character(len=*),dimension(ntypes),intent(in) :: atomnames
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(inout) :: mat
      character(len=*),intent(in) :: filename
      logical, intent(in) :: binary
      type(orbitals_data), intent(in) :: orbs !to be eventually eliminated? only needed for onwhichatom

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, iat, jat, ind, ispin, itype, iorb
      real(kind=8),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_dense_matrix')


      mat%matrix = sparsematrix_malloc_ptr(smat, iaction=DENSE_FULL, id='mat%matrix')
      call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, smat, &
           mat%matrix_compr, mat%matrix)

      if (bigdft_mpi%iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=binary)

          !write(iunit,'(i10,2i6,a)') nat, ntypes, smat%nspin, &
          !    '   # number of atoms, number of atom types, nspin'
          !do itype=1,ntypes
          !    write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
          !        '   # nz, nelpsp, name'
          !end do
          !do iat=1,nat
          !    write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          !end do
          !write(iunit,'(3i12,a)') smat%nfvctr, smat%nseg, smat%nvctr, '   # nfvctr, nseg, nvctr'
          !do iseg=1,smat%nseg
          !    write(iunit,'(5i12,a)') smat%keyv(iseg), smat%keyg(1,1,iseg), smat%keyg(2,1,iseg), &
          !        smat%keyg(1,2,iseg), smat%keyg(2,2,iseg), '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
          !end do
          !ind = 0
          !do ispin=1,smat%nspin
          !    do iseg=1,smat%nseg
          !        icol = smat%keyg(1,2,iseg)
          !        iat = smat%on_which_atom(icol)
          !        do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
          !            irow = jorb
          !            jat = smat%on_which_atom(irow)
          !            ind = ind + 1
          !            write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
          !        end do
          !    end do
          !end do

          !unify the structure with write_sparse?
          if (.not. binary) then
              write(iunit,'(a,3i10,a)') '#  ',smat%nfvctr, nat, smat%nspin, &
                  '    number of basis functions, number of atoms, number of spins'
          else
              write(iunit) '#  ',smat%nfvctr, nat, smat%nspin, &
                  '    number of basis functions, number of atoms, number of spins'
          end if
          do iat=1,nat
              if (.not. binary) then
                  write(iunit,'(a,3es24.16,a,i4.4)') '#  ',rxyz(1:3,iat), '   # position of atom no. ',iat
              else
                  write(iunit) '#  ',rxyz(1:3,iat)
              end if
          end do
    
          do ispin=1,smat%nspin
             do iorb=1,smat%nfvctr
                iat=orbs%onwhichatom(iorb)
                do jorb=1,smat%nfvctr
                   jat=orbs%onwhichatom(jorb)
                   if (.not. binary) then
                      write(iunit,'(2(i6,1x),es19.12,2(1x,i6),a)') iorb,jorb,mat%matrix(iorb,jorb,ispin),iat,jat, &
                          '   # i, j, mat(i,j), iat, jat'
                   else
                      write(iunit) iorb,jorb,mat%matrix(iorb,jorb,ispin),iat,jat
                   end if
                end do
             end do
          end do

          call f_close(iunit)

          call f_free_ptr(mat%matrix)
      end if

      call f_release_routine()

    end subroutine write_dense_matrix


    subroutine read_sparse_matrix(filename, nspin, geocode, nfvctr, nseg, nvctr, keyv, keyg, mat_compr, &
               nat, ntypes, nzatom, nelpsp, atomnames, iatype, rxyz, on_which_atom)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, nseg, nvctr
      character(len=1),intent(out) :: geocode
      integer,dimension(:),pointer,intent(out) :: keyv
      integer,dimension(:,:,:),pointer,intent(out) :: keyg
      real(kind=8),dimension(:),pointer,intent(out) :: mat_compr
      integer,intent(out),optional :: nat, ntypes
      integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      integer,dimension(:),pointer,intent(inout),optional :: on_which_atom

      ! Local variables
      integer :: iunit, dummy_int, iseg, icol, irow, jorb, ind, ispin, iat, ntypes_, nat_, itype
      real(kind=8) :: dummy_double
      character(len=20) :: dummy_char
      logical :: read_rxyz, read_on_which_atom

      call f_routine(id='read_sparse_matrix')

      if (present(nat) .and. present(ntypes) .and. present(nzatom) .and.  &
          present(nelpsp) .and. present(atomnames) .and. present(iatype) .and. present(rxyz)) then
          read_rxyz = .true.
      else if (present(nat) .or. present(ntypes) .or. present(nzatom) .or.  &
          present(nelpsp) .or. present(atomnames) .or. present(iatype) .or. present(rxyz)) then
          call f_err_throw("not all optional arguments were given", &
               err_name='BIGDFT_RUNTIME_ERROR')
      else
          read_rxyz = .false.
      end if
      
      if (present(on_which_atom)) then
          read_on_which_atom = .true.
      else
          read_on_which_atom = .false.
      end if

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      if (read_rxyz) then
          read(iunit,*) nat, ntypes, nspin, geocode
          nzatom = f_malloc_ptr(ntypes,id='nzatom')
          nelpsp = f_malloc_ptr(ntypes,id='nelpsp')
          atomnames = f_malloc0_str_ptr(len(atomnames),ntypes,id='atomnames')

          do itype=1,ntypes
              read(iunit,*) nzatom(itype), nelpsp(itype), atomnames(itype)
          end do
          rxyz = f_malloc_ptr((/3,nat/),id='rxyz')
          iatype = f_malloc_ptr(nat,id='iatype')
          do iat=1,nat
              read(iunit,*) iatype(iat), rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
          end do
      else
          read(iunit,*) nat_, ntypes_, nspin, geocode
          do itype=1,ntypes_
              read(iunit,*) dummy_int, dummy_int, dummy_char
          end do
          do iat=1,nat_
              read(iunit,*) dummy_int, dummy_double, dummy_double, dummy_double
          end do
      end if
      read(iunit,*) nfvctr, nseg, nvctr
      keyv = f_malloc_ptr(nseg,id='keyv')
      keyg = f_malloc_ptr((/2,2,nseg/),id='keyg')
      do iseg=1,nseg
          read(iunit,*) keyv(iseg), keyg(1,1,iseg), keyg(2,1,iseg), keyg(1,2,iseg), keyg(2,2,iseg)
      end do
      mat_compr = f_malloc_ptr(nvctr*nspin,id='mat_compr')
      if (read_on_which_atom) then
          nullify(on_which_atom)
          on_which_atom = f_malloc_ptr(nfvctr,id='on_which_atom')
          ind = 0
          do ispin=1,nspin
              do iseg=1,nseg
                  icol = keyg(1,2,iseg)
                  do jorb=keyg(1,1,iseg),keyg(2,1,iseg)
                      irow = jorb
                      ind = ind + 1
                      read(iunit,*) mat_compr(ind), on_which_atom(irow), on_which_atom(icol)
                  end do
              end do
          end do
      else
          ind = 0
          do ispin=1,nspin
              do iseg=1,nseg
                  icol = keyg(1,2,iseg)
                  do jorb=keyg(1,1,iseg),keyg(2,1,iseg)
                      irow = jorb
                      ind = ind + 1
                      read(iunit,*) mat_compr(ind), dummy_int, dummy_int
                  end do
              end do
          end do
      end if

      call f_close(iunit)

      call f_release_routine()

    end subroutine read_sparse_matrix

    !subroutine read_dense_matrix(nat, ntypes, iatype, rxyz, nzatom, nelpsp, atomnames, smat, mat, filename, binary, orbs)
    !eventually have the same header information as write_sparse?
    !differs in behaviour though as assumes all matrices are alrady allocated to correct size - think more about whether this is the best approach...
    subroutine read_dense_matrix(filename, binary, nspin, ntmb, mat, nat, rxyz, on_which_atom) 
               !ntypes, nzatom, nelpsp, atomnames, iatype)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      logical, intent(in) :: binary
      integer, intent(in) :: ntmb, nat, nspin
      real(kind=8),dimension(ntmb,ntmb,nspin),intent(inout) :: mat
      !integer,intent(out),optional :: ntypes
      !integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      !character(len=*),dimension(:),pointer,intent(inout),optional :: atomnames
      real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      integer,dimension(:),pointer,intent(inout),optional :: on_which_atom

      ! Local variables
      integer :: iunit, dummy_int, ispin, iat, iorb, jorb, ntmb_old, nat_old, nspin_old
      real(kind=8) :: dummy_double
      character(len=20) :: dummy_char
      logical :: read_rxyz, read_on_which_atom

      call f_routine(id='read_dense_matrix')

      if (present(rxyz)) then
          read_rxyz = .true.
      else
          read_rxyz = .false.
      end if
      
      if (present(on_which_atom)) then
          read_on_which_atom = .true.
      else
          read_on_which_atom = .false.
      end if

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=binary)

      !check heading information is consistent
      if (.not. binary) then
          read(iunit,*) dummy_char, ntmb_old, nat_old, nspin_old
      else
          read(iunit) dummy_char, ntmb_old, nat_old, nspin_old
      end if
      if (ntmb_old/=ntmb) call f_err_throw("Number of tmbs incorrect in read_dense_matrix", &
           err_name='BIGDFT_RUNTIME_ERROR')
      if (nat_old/=nat) call f_err_throw("Number of atoms incorrect in read_dense_matrix", &
           err_name='BIGDFT_RUNTIME_ERROR')
      if (nspin_old/=nspin) call f_err_throw("Number of spins incorrect in read_dense_matrix", &
           err_name='BIGDFT_RUNTIME_ERROR')

      if (read_rxyz) then
          do iat=1,nat
              if (.not. binary) then
                  read(iunit,*) dummy_char,rxyz(1:3,iat)
              else
                  read(iunit) dummy_char,rxyz(1:3,iat)
              end if
          end do  
      else
          do iat=1,nat
              read(iunit,*) dummy_char, dummy_double, dummy_double, dummy_double
          end do
      end if

      if (read_on_which_atom) then
          do ispin=1,nspin
             do iorb=1,ntmb
                do jorb=1,ntmb
                   if (.not. binary) then
                      read(iunit,*) dummy_int,dummy_int,mat(iorb,jorb,ispin),on_which_atom(iorb),on_which_atom(jorb)
                   else
                      read(iunit) dummy_int,dummy_int,mat(iorb,jorb,ispin),on_which_atom(iorb),on_which_atom(jorb)
                   end if
                end do
             end do
          end do
      else
          do ispin=1,nspin
             do iorb=1,ntmb
                do jorb=1,ntmb
                   if (.not. binary) then
                      read(iunit,*) dummy_int,dummy_int,mat(iorb,jorb,ispin),dummy_int,dummy_int
                   else
                      read(iunit) dummy_int,dummy_int,mat(iorb,jorb,ispin),dummy_int,dummy_int
                   end if
                end do
             end do
          end do
      end if

      call f_close(iunit)

      call f_release_routine()

    end subroutine read_dense_matrix


    !> Write Hamiltonian, overlap and kernel matrices in tmb basis
    subroutine write_linear_matrices(iproc,nproc,imethod_overlap,filename,iformat,tmb,at,rxyz,norder_taylor, &
               calculate_onsite_overlap, write_SminusonehalfH)
      use module_types
      use module_base
      use yaml_output
      use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc0_ptr, matrices_null, &
                                   DENSE_FULL, SPARSE_TASKGROUP, assignment(=), &
                                   deallocate_matrices
      use sparsematrix, only: uncompress_matrix2, transform_sparse_matrix_local, matrix_matrix_mult_wrapper
      use matrix_operations, only: overlapPowerGeneral
      implicit none
      integer, intent(in) :: iproc,nproc,imethod_overlap,norder_taylor
      integer,intent(in) :: iformat !< 1: plain sparse, 11: plain dense, 21: plain both (later extend to other than plain formats...)
      character(len=*), intent(in) :: filename 
      type(DFT_wavefunction), intent(inout) :: tmb
      type(atoms_data), intent(in) :: at
      real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz
      logical,intent(in) :: calculate_onsite_overlap, write_SminusonehalfH
      !local variables
      logical :: binary
      integer :: ispin, iorb, jorb, iat, jat,unitm
      !!integer :: i_stat, i_all
      character(len=*),parameter :: subname='write_linear_matrices'
      logical :: write_sparse, write_dense
      type(matrices),dimension(1) :: SminusonehalfH
      real(kind=8),dimension(:),pointer :: ham_large, tmp_large
      real(kind=8) :: max_error, mean_error
    
      call f_routine(id='write_linear_matrices')

      write_sparse = (1<=iformat .and. iformat<=9) .or. (21<=iformat .and.  iformat<=29)
      write_dense = (11<=iformat .and. iformat<=19) .or. (21<=iformat .and.  iformat<=29)
    
      unitm=99
      binary=(mod(iformat,10) /= WF_FORMAT_PLAIN)
    
      if (write_dense) then
          call write_dense_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
               at%nzatom, at%nelpsp, at%astruct%atomnames, &
               tmb%linmat%m, tmb%linmat%ham_, trim(filename//'hamiltonian.bin'), binary, tmb%orbs)
      end if

      if (write_sparse) then
          call write_sparse_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
               at%nzatom, at%nelpsp, at%astruct%atomnames, &
               tmb%linmat%m, tmb%linmat%ham_, trim(filename//'hamiltonian_sparse.bin'))
      end if
    
    
      if (write_dense) then
          call write_dense_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
               at%nzatom, at%nelpsp, at%astruct%atomnames, &
               tmb%linmat%s, tmb%linmat%ovrlp_, trim(filename//'overlap.bin'), binary, tmb%orbs)
      end if

      if (write_sparse) then
          call write_sparse_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
               at%nzatom, at%nelpsp, at%astruct%atomnames, &
               tmb%linmat%s, tmb%linmat%ovrlp_, filename//'overlap_sparse.bin')
      end if
    
    
      if (write_dense) then
          call write_dense_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
               at%nzatom, at%nelpsp, at%astruct%atomnames, &
               tmb%linmat%l, tmb%linmat%kernel_, trim(filename//'density_kernel.bin'), binary, tmb%orbs)
      end if

      if (write_sparse) then
          call write_sparse_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
               at%nzatom, at%nelpsp, at%astruct%atomnames, &
               tmb%linmat%l, tmb%linmat%kernel_, filename//'density_kernel_sparse.bin')
      end if
    
    
      if (calculate_onsite_overlap) then
          ! calculate 'onsite' overlap matrix as well - needs double checking
    
          if (write_dense) then
              tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, &
                                         id='tmb%linmat%ovrlp_%matrix')
    
              call tmb_overlap_onsite(iproc, nproc, imethod_overlap, at, tmb, rxyz)
              !call tmb_overlap_onsite_rotate(iproc, nproc, at, tmb, rxyz)
    
              if (iproc==0) then
                 !if(iformat == WF_FORMAT_PLAIN) then
                 call f_open_file(unitm,file=filename//'overlap_onsite.bin',&
                      binary=binary)
                 !else
                 !open(99, file=filename//'overlap_onsite.bin', status='unknown',form='unformatted')
                 !end if
    
                 if (.not. binary) then
                     write(unitm,'(a,2i10,a)') '#  ',tmb%linmat%m%nfvctr, at%astruct%nat, &
                         '    number of basis functions, number of atoms'
                 else
                     write(unitm) '#  ',tmb%linmat%m%nfvctr, at%astruct%nat, &
                         '    number of basis functions, number of atoms'
                 end if
                 do iat=1,at%astruct%nat
                     if (.not. binary) then
                         write(unitm,'(a,3es24.16)') '#  ',rxyz(1:3,iat)
                     else
                         write(unitm) '#  ',rxyz(1:3,iat)
                     end if
                 end do
    
                 do ispin=1,tmb%linmat%l%nspin
                    do iorb=1,tmb%linmat%l%nfvctr
                       iat=tmb%orbs%onwhichatom(iorb)
                       do jorb=1,tmb%linmat%l%nfvctr
                          jat=tmb%orbs%onwhichatom(jorb)
                          if (.not. binary) then
                             write(unitm,'(2(i6,1x),e19.12,2(1x,i6))') iorb,jorb,tmb%linmat%ovrlp_%matrix(iorb,jorb,ispin),iat,jat
                          else
                             write(unitm) iorb,jorb,tmb%linmat%ovrlp_%matrix(iorb,jorb,ispin),iat,jat
                          end if
                       end do
                    end do
                 end do
    
                 call f_close(unitm)
    
              end if

              call f_free_ptr(tmb%linmat%ovrlp_%matrix)

          end if

      end if


      if (write_SminusonehalfH) then
          SminusonehalfH(1) = matrices_null()
          SminusonehalfH(1)%matrix_compr = &
              sparsematrix_malloc0_ptr(tmb%linmat%l,iaction=SPARSE_TASKGROUP,id='SminusonehalfH%matrix_compr')
          ham_large = sparsematrix_malloc0_ptr(tmb%linmat%l,iaction=SPARSE_TASKGROUP,id='ham_large')
          tmp_large = sparsematrix_malloc0_ptr(tmb%linmat%l,iaction=SPARSE_TASKGROUP,id='tmp_large')
          call transform_sparse_matrix_local(tmb%linmat%m, tmb%linmat%l, &
               tmb%linmat%ham_%matrix_compr, ham_large, 'small_to_large')
          ! calculate S^-1/2
          call overlapPowerGeneral(iproc, nproc, norder_taylor, 1, (/-2/), -1, &
               imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
               ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=SminusonehalfH(1), &
               check_accur=.true., max_error=max_error, mean_error=mean_error)
          ! Calculate S^-1/2 * H
          call f_memcpy(src=SminusonehalfH(1)%matrix_compr,dest=tmp_large)
          call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, tmp_large, ham_large, SminusonehalfH(1)%matrix_compr)
          call f_free_ptr(ham_large)
          call f_free_ptr(tmp_large)
          if (write_dense) then
              call write_dense_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
                   at%nzatom, at%nelpsp, at%astruct%atomnames, &
                   tmb%linmat%l, SminusonehalfH(1), trim(filename//'SminusonehalfH.bin'), binary, tmb%orbs)
          end if

          if (write_sparse) then
              call write_sparse_matrix(at%astruct%nat, at%astruct%ntypes, at%astruct%iatype, at%astruct%rxyz, &
                   at%nzatom, at%nelpsp, at%astruct%atomnames, &
                   tmb%linmat%s, SminusonehalfH(1), filename//'SminusonehalfH_sparse.bin')
          end if
          call deallocate_matrices(SminusonehalfH(1))
      end if

      call f_release_routine()
    
    end subroutine write_linear_matrices


    subroutine writeLinearCoefficients(unitwf,useFormattedOutput,nat,rxyz,&
               ntmb,norb,nfvctr,coeff,eval)
      use module_base
      use yaml_output
      implicit none
      logical, intent(in) :: useFormattedOutput
      integer, intent(in) :: unitwf,nat,ntmb,norb,nfvctr
      real(wp), dimension(nfvctr,ntmb), intent(in) :: coeff
      real(wp), dimension(ntmb), intent(in) :: eval
      real(gp), dimension(3,nat), intent(in) :: rxyz
      !local variables
      integer :: iat,i,j,iorb
      real(wp) :: tt
    
      ! Write the Header
      if (useFormattedOutput) then
         write(unitwf,*) ntmb,norb
         write(unitwf,*) nat
         do iat=1,nat
         write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
         enddo
         do iorb=1,ntmb
         write(unitwf,*) iorb,eval(iorb)
         enddo
      else
         write(unitwf) ntmb, norb
         write(unitwf) nat
         do iat=1,nat
         write(unitwf) (rxyz(j,iat),j=1,3)
         enddo
         do iorb=1,ntmb
         write(unitwf) iorb,eval(iorb)
         enddo
      end if
    
      ! Now write the coefficients
      do i = 1, ntmb
         ! first element always positive, for consistency when using for transfer integrals
         ! unless 1st element below some threshold, in which case first significant element
         do j=1,nfvctr
            if (abs(coeff(j,i))>1.0e-1) then
               if (coeff(j,i)<0.0_gp) call dscal(ntmb,-1.0_gp,coeff(1,i),1)
               exit
            end if
         end do
         !if (j==ntmb+1)print*,'Error finding significant coefficient, coefficients not scaled to have +ve first element'
    
         do j = 1,nfvctr
              tt = coeff(j,i)
              if (useFormattedOutput) then
                 write(unitwf,'(2(i6,1x),e19.12)') i,j,tt
              else
                 write(unitwf) i,j,tt
              end if
         end do
      end do  
      if (verbose >= 2 .and. bigdft_mpi%iproc==0) call yaml_map('Wavefunction coefficients written',.true.)
    
    END SUBROUTINE writeLinearCoefficients


    !> Basically the same as writeLinearCoefficients, but with a slightly different format
    subroutine write_linear_coefficients(iroot, filename, at, rxyz, nfvctr, ntmb, nspin, coeff, eval)
      use module_base
      use module_types
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(atoms_data),intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      integer,intent(in) :: iroot, nfvctr, ntmb, nspin
      real(wp), dimension(nfvctr,ntmb), intent(in) :: coeff
      real(wp), dimension(ntmb), intent(in) :: eval
      ! Local variables
      integer :: iunit, itype, iat, i, j
      logical :: scaled

      call f_routine(id='write_linear_coefficients')


      if (bigdft_mpi%iproc==iroot) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)
    
          ! Write the Header
          write(iunit,'(i10,2i6,a)') at%astruct%nat, at%astruct%ntypes, nspin, &
              '   # number of atoms, number of atom types, nspin'
          do itype=1,at%astruct%ntypes
              write(iunit,'(2i8,3x,a,a)') at%nzatom(itype), at%nelpsp(itype), trim(at%astruct%atomnames(itype)), &
                  '   # nz, nelpsp, name'
          end do
          do iat=1,at%astruct%nat
              write(iunit,'(i5, 3es24.16,a,i0)') at%astruct%iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          end do
          write(iunit,'(2i12,a)') nfvctr, ntmb, '   # nfvctr, ntmb'
          do i=1,ntmb
              write(iunit,'(es24.16,a,i0)') eval(i), '   # eval no. ', i
          enddo
    
          ! Now write the coefficients
          do i=1,ntmb
             ! First element always positive, for consistency when using for transfer integrals;
             ! unless 1st element below some threshold, in which case first significant element.
             scaled = .false.
             do j=1,nfvctr
                if (abs(coeff(j,i))>1.0d-3) then
                   if (coeff(j,i)<0.0_gp) call dscal(ntmb,-1.0_gp,coeff(1,i),1)
                   scaled = .true.
                   exit
                end if
             end do
             if (.not.scaled) then
                 call yaml_warning('Consistency between the written coefficients not guaranteed')
             end if
    
             do j = 1,nfvctr
                 write(iunit,'(es24.16,2i9,a)') coeff(j,i), j, i, '   # coeff, j, i'
             end do
          end do  
          if (verbose >= 2 .and. bigdft_mpi%iproc==0) call yaml_map('Wavefunction coefficients written',.true.)

          call f_close(iunit)

      end if

      call f_release_routine()
    
    end subroutine write_linear_coefficients


    subroutine read_linear_coefficients(filename, nspin, nfvctr, ntmb, coeff, &
               nat, ntypes, nzatom, nelpsp, iatype, atomnames, rxyz, eval)
      use module_base
      use module_types
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, ntmb
      real(kind=8),dimension(:,:),pointer,intent(inout) :: coeff
      integer,intent(out),optional :: nat, ntypes
      integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      real(kind=8),dimension(:),pointer,intent(inout),optional :: eval
      ! Local variables
      real(kind=8) :: dummy_double
      character(len=20) :: dummy_char
      integer :: iunit, itype, iat, i, j, dummy_int, ntypes_, nat_
      logical :: scaled, read_rxyz, read_eval

      call f_routine(id='write_linear_coefficients')

      if (present(nat) .and. present(ntypes) .and. present(nzatom) .and.  &
          present(nelpsp) .and. present(atomnames) .and. present(iatype) .and. present(rxyz)) then
          read_rxyz = .true.
      else if (present(nat) .or. present(ntypes) .or. present(nzatom) .or.  &
          present(nelpsp) .or. present(atomnames) .or. present(iatype) .or. present(rxyz)) then
          call f_err_throw("not all optional arguments were given", &
               err_name='BIGDFT_RUNTIME_ERROR')
      else
          read_rxyz = .false.
      end if

      if (present(eval)) then
          read_eval = .true.
      else
          read_eval = .false.
      end if

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)
    
      ! Read the Header
      if (read_rxyz) then
          read(iunit,*) nat, ntypes, nspin
          nzatom = f_malloc_ptr(ntypes,id='nzatom')
          nelpsp = f_malloc_ptr(ntypes,id='nelpsp')
          atomnames = f_malloc0_str_ptr(len(atomnames),ntypes,id='atomnames')

          do itype=1,ntypes
              read(iunit,*) nzatom(itype), nelpsp(itype), atomnames(itype)
          end do
          rxyz = f_malloc_ptr((/3,nat/),id='rxyz')
          iatype = f_malloc_ptr(nat,id='iatype')
          do iat=1,nat
              read(iunit,*) iatype(iat), rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
          end do
      else
          read(iunit,*) nat_, ntypes_, nspin
          do itype=1,ntypes_
              read(iunit,*) dummy_int, dummy_int, dummy_char
          end do
          do iat=1,nat_
              read(iunit,*) dummy_int, dummy_double, dummy_double, dummy_double
          end do
      end if

      read(iunit,*) nfvctr, ntmb

      if (read_eval) then
          eval = f_malloc_ptr(ntmb,id='eval')
          do i=1,ntmb
              read(iunit,*) eval(i)
          end do
      else
          do i=1,ntmb
              read(iunit,*) dummy_double
          end do
      end if
    
      ! Now read the coefficients
      coeff = f_malloc_ptr((/nfvctr,ntmb/),id='coeff')

      do i=1,ntmb

         do j = 1,nfvctr
             read(iunit,*) coeff(j,i)
         end do

         ! First element always positive, for consistency when using for transfer integrals;
         ! unless 1st element below some threshold, in which case first significant element.
         scaled = .false.
         do j=1,nfvctr
            if (abs(coeff(j,i))>1.0d-3) then
               if (coeff(j,i)<0.0_gp) call dscal(ntmb,-1.0_gp,coeff(1,i),1)
               scaled = .true.
               exit
            end if
         end do
         if (.not.scaled) then
             call yaml_warning('Consistency between the written coefficients not guaranteed')
         end if
    
      end do  

      call f_close(iunit)


      call f_release_routine()
    
    end subroutine read_linear_coefficients



    subroutine write_ccs_matrix(filename, nfvctr, nvctr, row_ind, col_ptr, mat_compr)
      use module_base
      implicit none
      !Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: nfvctr !number of rows/columns
      integer,intent(in) :: nvctr !number of non-zero elements
      integer,dimension(nvctr),intent(in) :: row_ind
      integer,dimension(nfvctr),intent(in) :: col_ptr
      real(kind=8),dimension(nvctr),intent(in) :: mat_compr
      ! Local variables
      integer :: i, iunit

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      write(iunit,'(4(i0,1x))') nfvctr, nfvctr, nvctr, 0
      write(iunit,'(100000(i0,1x))') (col_ptr(i),i=1,nfvctr),nvctr+1
      write(iunit,'(100000(i0,1x))') (row_ind(i),i=1,nvctr)
      do i=1,nvctr
          write(iunit,'(es24.15)') mat_compr(i)
      end do

      call f_close(iunit)

    end subroutine write_ccs_matrix


    subroutine write_partial_charges(atoms, charge_per_atom, write_gnuplot)
      use module_base
      use module_types
      use yaml_output
      ! Calling arguments
      type(atoms_data),intent(in) :: atoms
      real(kind=8),dimension(atoms%astruct%nat),intent(in) :: charge_per_atom
      logical,intent(in) :: write_gnuplot
      ! Local variables
      integer :: iat, itypes, iitype, nntype, intype, iunit
      real(kind=8) :: total_charge, total_net_charge, frac_charge, range_min, range_max
      character(len=20) :: atomname, colorname
      real(kind=8),dimension(2) :: charges
      character(len=128) :: output
      character(len=2) :: backslash
      integer,parameter :: ncolors = 12 !7
      !character(len=20),dimension(ncolors),parameter :: colors=(/'violet', &
      !                                                           'blue  ', &
      !                                                           'cyan  ', &
      !                                                           'green ', &
      !                                                           'yellow', &
      !                                                           'orange', &
      !                                                           'red   '/)
      ! Presumably well suited colorschemes from colorbrewer2.org
      character(len=20),dimension(ncolors),parameter :: colors=(/'#a6cee3', &
                                                                 '#1f78b4', &
                                                                 '#b2df8a', &
                                                                 '#33a02c', &
                                                                 '#fb9a99', &
                                                                 '#e31a1c', &
                                                                 '#fdbf6f', &
                                                                 '#ff7f00', &
                                                                 '#cab2d6', &
                                                                 '#6a3d9a', &
                                                                 '#ffff99', &
                                                                 '#b15928'/)

      call yaml_sequence_open('Charge analysis (charge / net charge)')
      total_charge=0.d0
      total_net_charge=0.d0
      do iat=1,atoms%astruct%nat
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          atomname=atoms%astruct%atomnames(atoms%astruct%iatype(iat))
          charges(1)=-charge_per_atom(iat)
          charges(2)=-(charge_per_atom(iat)-real(atoms%nelpsp(atoms%astruct%iatype(iat)),kind=8))
          total_charge = total_charge + charges(1)
          total_net_charge = total_net_charge + charges(2)
          call yaml_map(trim(atomname),charges,fmt='(1es20.12)')
          call yaml_mapping_close(advance='no')
          call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
      end do
      call yaml_sequence(advance='no')
      call yaml_map('total charge',total_charge,fmt='(es16.8)')
      call yaml_sequence(advance='no')
      call yaml_map('total net charge',total_net_charge,fmt='(es16.8)')
      call yaml_sequence_close()

      if (write_gnuplot) then
          output='chargeanalysis.gp'
          call yaml_map('output file',trim(output))
          iunit=100
          call f_open_file(iunit, file=trim(output), binary=.false.)
          write(iunit,'(a)') '# plot the fractional charge as a normalized sum of Gaussians'
          write(iunit,'(a)') 'set samples 1000'
          range_min = minval(-(charge_per_atom(:)-real(atoms%nelpsp(atoms%astruct%iatype(:)),kind=8))) - 0.1d0
          range_max = maxval(-(charge_per_atom(:)-real(atoms%nelpsp(atoms%astruct%iatype(:)),kind=8))) + 0.1d0
          write(iunit,'(a,2(es12.5,a))') 'set xrange[',range_min,':',range_max,']'
          write(iunit,'(a)') 'sigma=0.005'
          write(backslash,'(a)') '\ '
          do itypes=1,atoms%astruct%ntypes
              nntype = 0
              do iat=1,atoms%astruct%nat
                  iitype = (atoms%astruct%iatype(iat))
                  if (iitype==itypes) then
                      nntype = nntype + 1
                  end if
              end do
              write(iunit,'(a,i0,a,i0,2a)') 'f',itypes,'(x) = 1/',nntype,'.0*( '//trim(backslash)
              intype = 0
              do iat=1,atoms%astruct%nat
                  iitype = (atoms%astruct%iatype(iat))
                  if (iitype==itypes) then
                      intype = intype + 1
                      frac_charge = -(charge_per_atom(iat)-real(atoms%nelpsp(atoms%astruct%iatype(iat)),kind=8))
                      if (intype<nntype) then
                          write(iunit,'(a,es16.9,a)') '  1.0*exp(-(x-',frac_charge,')**2/(2*sigma**2)) + '//trim(backslash)
                      else
                          write(iunit,'(a,es16.9,a)') '  1.0*exp(-(x-',frac_charge,')**2/(2*sigma**2)))'
                      end if
                  end if
              end do
              atomname=atoms%astruct%atomnames(itypes)
              if (itypes<ncolors) then
                  colorname = colors(itypes)
              else
                  colorname = 'color'
              end if
              if (itypes==1) then
                  write(iunit,'(a,i0,5a)') "plot f",itypes,"(x) lc rgb '",trim(colorname), &
                      "' lt 1 lw 2 w l title '",trim(atomname),"'"
              else
                  write(iunit,'(a,i0,5a)') "replot f",itypes,"(x) lc rgb '",trim(colorname), &
                      "' lt 1 lw 2 w l title '",trim(atomname),"'"
              end if
          end do
      end if
    end subroutine write_partial_charges


end module io
