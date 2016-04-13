module sparsematrix_wrappers
  use module_base
  use sparsematrix_base
  implicit none

  private

  !> Public routines
  public :: init_sparse_matrix_wrapper
  public :: init_sparse_matrix_for_KSorbs
  public :: check_kernel_cutoff

  contains

    subroutine init_sparse_matrix_wrapper(iproc, nproc, nspin, orbs, lzd, astruct, &
               store_index, init_matmul, imode, smat, smat_ref)
      use module_types, only: orbitals_data, local_zone_descriptors, atomic_structure
      use sparsematrix_init, only: init_sparse_matrix
      implicit none
  
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin, imode
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(atomic_structure),intent(in) :: astruct
      logical,intent(in) :: store_index, init_matmul
      type(sparse_matrix),intent(out) :: smat
      type(sparse_matrix),intent(in),optional :: smat_ref !< reference sparsity pattern, in case smat must be at least as large as smat_ref
  
      ! Local variables
      integer :: nnonzero, nnonzero_mult, ilr
      integer,dimension(:,:),pointer :: nonzero, nonzero_mult
      real(kind=8),dimension(:),allocatable :: locrad
      logical :: present_smat_ref
      integer,parameter :: KEYS=1
      integer,parameter :: DISTANCE=2
  
      call f_routine(id='init_sparse_matrix_wrapper')
  
      present_smat_ref = present(smat_ref)
  
      locrad = f_malloc(lzd%nlr,id='lzd%nlr')
  
      if (imode==KEYS) then
         call determine_sparsity_pattern(iproc, nproc, orbs, lzd, nnonzero, nonzero)
      else if (imode==DISTANCE) then
          do ilr=1,lzd%nlr
              locrad(ilr) = lzd%llr(ilr)%locrad_kernel
          end do
         if (present_smat_ref) then
            call determine_sparsity_pattern_distance(orbs, lzd, astruct, locrad, nnonzero, nonzero, smat_ref)
         else
            call determine_sparsity_pattern_distance(orbs, lzd, astruct, locrad, nnonzero, nonzero)
         end if
      else
         stop 'wrong imode'
      end if
  
      ! Make sure that the cutoff for the multiplications is larger than the kernel cutoff
      do ilr=1,lzd%nlr
         !write(*,*) 'lzd%llr(ilr)%locrad_mult, lzd%llr(ilr)%locrad_kernel', lzd%llr(ilr)%locrad_mult, lzd%llr(ilr)%locrad_kernel
         if (lzd%llr(ilr)%locrad_mult<lzd%llr(ilr)%locrad_kernel) then
            call f_err_throw('rloc_kernel_foe ('//trim(yaml_toa(lzd%llr(ilr)%locrad_mult,fmt='(f5.2)'))//&
                 &') too small, must be at least as big as rloc_kernel('&
                 &//trim(yaml_toa(lzd%llr(ilr)%locrad_kernel,fmt='(f5.2)'))//')', err_id=BIGDFT_RUNTIME_ERROR)
         end if
      end do
  
      do ilr=1,lzd%nlr
          locrad(ilr) = lzd%llr(ilr)%locrad_mult
      end do
      if (present_smat_ref) then
         call determine_sparsity_pattern_distance(orbs, lzd, astruct, locrad, &
              nnonzero_mult, nonzero_mult, smat_ref)
      else
         call determine_sparsity_pattern_distance(orbs, lzd, astruct, locrad, &
              nnonzero_mult, nonzero_mult)
      end if
      call init_sparse_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
           orbs%norbu, nnonzero, nonzero, nnonzero_mult, nonzero_mult, smat, &
           init_matmul=init_matmul, nspin=nspin, geocode=astruct%geocode, &
           cell_dim=astruct%cell_dim, norbup=orbs%norbup, &
           isorbu=orbs%isorbu, store_index=store_index, on_which_atom=orbs%onwhichatom)
      call f_free_ptr(nonzero)
      call f_free_ptr(nonzero_mult)
      call f_free(locrad)
  
      call f_release_routine()
  
    end subroutine init_sparse_matrix_wrapper
  
  
    subroutine check_kernel_cutoff(iproc, orbs, atoms, hamapp_radius_incr, lzd)
      use module_types
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, hamapp_radius_incr
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
          call yaml_sequence_open('check of kernel cutoff radius')
      end if
    
      do iorb=1,orbs%norb
          ilr=orbs%inwhichlocreg(iorb)
    
          ! cutoff radius of the support function, including shamop region
          cutoff_sf=lzd%llr(ilr)%locrad+real(hamapp_radius_incr,kind=8)*lzd%hgrids(1)
    
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
              call yaml_mapping_open(flow=.true.)
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
              call yaml_mapping_close()
          end if
      end do
    
      if (iproc==0) then
          call yaml_sequence_close
      end if
    
    
    end subroutine check_kernel_cutoff
  
  
  
    subroutine determine_sparsity_pattern(iproc, nproc, orbs, lzd, nnonzero, nonzero)
          use module_types
          use locregs, only: check_overlap_cubic_periodic,check_overlap_from_descriptors_periodic
          implicit none
          ! Calling arguments
          integer, intent(in) :: iproc, nproc
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: lzd
          integer, intent(out) :: nnonzero
          integer, dimension(:,:), pointer,intent(out) :: nonzero
          ! Local variables
          integer :: iorb, jorb, ioverlaporb, ilr, jlr, ilrold
          integer :: iiorb, ii
          !!integer :: istat
          logical :: isoverlap
          integer :: onseg
          logical, dimension(:,:), allocatable :: overlapMatrix
          integer, dimension(:), allocatable :: noverlapsarr
          integer, dimension(:,:), allocatable :: overlaps_op
          !character(len=*), parameter :: subname='determine_overlap_from_descriptors'
    
          call f_routine('determine_sparsity_pattern')
          call timing(iproc,'determinespars','ON')
        
          overlapMatrix = f_malloc((/orbs%norbu,maxval(orbs%norbu_par(:,0))/),id='overlapMatrix')
          noverlapsarr = f_malloc(orbs%norbup,id='noverlapsarr')
        
          overlapMatrix=.false.
  
          do iorb=1,orbs%norbup
             ioverlaporb=0 ! counts the overlaps for the given orbital.
             iiorb=orbs%isorbu+iorb
             ilr=orbs%inWhichLocreg(iiorb)
  
             !$omp parallel default(none) &
             !$omp private(jorb,jlr,isoverlap,onseg) &
             !$omp shared(orbs,lzd,iorb,ilr,overlapMatrix,ioverlaporb)
             !$omp do schedule(guided) reduction(+:ioverlaporb)
             do jorb=1,orbs%norbu
                jlr=orbs%inWhichLocreg(jorb)
                call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzd%llr(jlr),isoverlap)
                !write(*,'(a,3(6i6,4x),l4)') 'is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ns1, ne1, ns2, ne2, ns3, ne3, isoverlap', &
                !    lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns1+lzd%llr(ilr)%d%n1, &
                !    lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns2+lzd%llr(ilr)%d%n2, &
                !    lzd%llr(ilr)%ns3, lzd%llr(ilr)%ns3+lzd%llr(ilr)%d%n3, &
                !    lzd%llr(jlr)%ns1, lzd%llr(jlr)%ns1+lzd%llr(jlr)%d%n1, &
                !    lzd%llr(jlr)%ns2, lzd%llr(jlr)%ns2+lzd%llr(jlr)%d%n2, &
                !    lzd%llr(jlr)%ns3, lzd%llr(jlr)%ns3+lzd%llr(jlr)%d%n3, &
                !    lzd%glr%ns1, lzd%glr%ns1+lzd%glr%d%n1, &
                !    lzd%glr%ns2, lzd%glr%ns2+lzd%glr%d%n2, &
                !    lzd%glr%ns3, lzd%glr%ns3+lzd%glr%d%n3, &
                !    isoverlap
                if(isoverlap) then
                   ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
                   ! Now explicitly check whether there is an overlap by using the descriptors.
                   call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c,&
                        lzd%llr(ilr)%wfd%keyglob, lzd%llr(jlr)%wfd%keyglob, &
                        isoverlap, onseg)
                   if(isoverlap) then
                      ! There is really an overlap
                      overlapMatrix(jorb,iorb)=.true.
                      ioverlaporb=ioverlaporb+1
                   else
                      overlapMatrix(jorb,iorb)=.false.
                   end if
                else
                   overlapMatrix(jorb,iorb)=.false.
                end if
                !!write(*,'(a,2i8,l4)') 'iiorb, jorb, isoverlap', iiorb, jorb, isoverlap
             end do
             !$omp end do
             !$omp end parallel     
             noverlapsarr(iorb)=ioverlaporb
          end do
   
  
  
          overlaps_op = f_malloc((/maxval(noverlapsarr),orbs%norbup/),id='overlaps_op')
        
          ! Now we know how many overlaps have to be calculated, so determine which orbital overlaps
          ! with which one. This is essentially the same loop as above, but we use the array 'overlapMatrix'
          ! which indicates the overlaps.
          !$omp parallel if (orbs%norbu>1000) &
          !$omp default(none) &
          !$omp shared(orbs,overlapMatrix,overlaps_op) &
          !$omp private(iorb,ioverlaporb,iiorb)
          !$omp do schedule(guided)
          do iorb=1,orbs%norbup
             ioverlaporb=0 ! counts the overlaps for the given orbital.
             iiorb=orbs%isorbu+iorb
             do jorb=1,orbs%norbu
                if(overlapMatrix(jorb,iorb)) then
                   ioverlaporb=ioverlaporb+1
                   overlaps_op(ioverlaporb,iorb)=jorb
                end if
             end do 
          end do
          !$omp end do
          !$omp end parallel
    
    
          nnonzero=0
          do iorb=1,orbs%norbup
              nnonzero=nnonzero+noverlapsarr(iorb)
          end do
          nonzero = f_malloc_ptr((/2,nnonzero/),id='nonzero')
          ii=0
          do iorb=1,orbs%norbup
              iiorb=orbs%isorbu+iorb
              do jorb=1,noverlapsarr(iorb)
                  ii=ii+1
                  nonzero(1,ii)=overlaps_op(jorb,iorb)
                  nonzero(2,ii)=iiorb
              end do
          end do
    
          call f_free(overlapMatrix)
          call f_free(noverlapsarr)
          call f_free(overlaps_op)
        
          call timing(iproc,'determinespars','OF')
          call f_release_routine()
    
    end subroutine determine_sparsity_pattern
  
  
    !> Initializes a sparse matrix type compatible with the ditribution of the KS orbitals
    subroutine init_sparse_matrix_for_KSorbs(iproc, nproc, orbs, input, geocode, cell_dim, nextra, smat, smat_extra)
      use module_types
      use module_interfaces, only: orbitals_descriptors
      use public_enums
      use sparsematrix_init, only: init_sparse_matrix
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: iproc, nproc, nextra
      type(orbitals_data), intent(in) :: orbs
      type(input_variables), intent(in) :: input
      character(len=1),intent(in) :: geocode
      real(kind=8),dimension(3),intent(in) :: cell_dim
      type(sparse_matrix),dimension(:),pointer,intent(out) :: smat, smat_extra
    
      ! Local variables
      integer :: i, iorb, iiorb, jorb, ind, norb, norbp, isorb, ispin
      integer,dimension(:,:),allocatable :: nonzero
      type(orbitals_data) :: orbs_aux
      character(len=*), parameter :: subname='init_sparse_matrix_for_KSorbs'
    
      call f_routine('init_sparse_matrix_for_KSorbs')
    
    
      allocate(smat(input%nspin))
      allocate(smat_extra(input%nspin))
    
    
      ! First the type for the normal KS orbitals distribution
      do ispin=1,input%nspin
    
          smat(ispin) = sparse_matrix_null()
          smat_extra(ispin) = sparse_matrix_null()
    
          if (ispin==1) then
              norb=orbs%norbu
              norbp=orbs%norbup
              isorb=orbs%isorbu
          else
              norb=orbs%norbd
              norbp=orbs%norbdp
              isorb=orbs%isorbd
          end if
    
          nonzero = f_malloc((/2,norb*norbp/), id='nonzero')
          i=0
          do iorb=1,norbp
              iiorb=isorb+iorb
              do jorb=1,norb
                  i=i+1
                  ind=(iiorb-1)*norb+jorb
                  nonzero(1,i)=jorb
                  nonzero(2,i)=iiorb
              end do
          end do
          call init_sparse_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
               norb, norb*norbp, nonzero, norb*norbp, nonzero, smat(ispin), &
               nspin=input%nspin, geocode=geocode, cell_dim=cell_dim, norbup=norbp, isorbu=isorb, &
               store_index=input%store_index, on_which_atom=orbs%onwhichatom, print_info=.false.)
          call f_free(nonzero)
    
    
          !SM: WARNING: not tested whether the spin works here! Mainly just to create a
          !spin down part and make the compiler happy at another location.
          ! Now the distribution for the KS orbitals including the extr states. Requires
          ! first to calculate a corresponding orbs type.
          call nullify_orbitals_data(orbs_aux)
          call orbitals_descriptors(iproc, nproc, norb+nextra, norb+nextra, 0, input%nspin, orbs%nspinor,&
               input%gen_nkpt, input%gen_kpt, input%gen_wkpt, orbs_aux, LINEAR_PARTITION_NONE)
          nonzero = f_malloc((/2,orbs_aux%norbu*orbs_aux%norbup/), id='nonzero')
          !write(*,*) 'iproc, norb, norbp, norbu, norbup', iproc, orbs_aux%norb, orbs_aux%norbp, orbs_aux%norbu, orbs_aux%norbup
          i=0
          do iorb=1,orbs_aux%norbup
              iiorb=orbs_aux%isorbu+iorb
              do jorb=1,orbs_aux%norbu
                  i=i+1
                  ind=(iiorb-1)*orbs_aux%norbu+jorb
                  nonzero(1,i)=jorb
                  nonzero(2,i)=iiorb
              end do
          end do
          !!call init_sparse_matrix(iproc, nproc, input%nspin, orbs_aux%norb, orbs_aux%norbp, orbs_aux%isorb, &
          !!     orbs%norbu, orbs%norbup, orbs%isorbu, input%store_index, &
          !!     orbs_aux%norbu*orbs_aux%norbup, nonzero, orbs_aux%norbu, nonzero, smat_extra, print_info_=.false.)
          !!call init_sparse_matrix(iproc, nproc, input%nspin, orbs_aux%norb, orbs_aux%norbp, orbs_aux%isorb, &
          !!     norb, norbp, isorb, input%store_index, &
          !!     orbs_aux%norbu*orbs_aux%norbup, nonzero, orbs_aux%norbu, nonzero, smat_extra(ispin), print_info_=.false.)
          call init_sparse_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
               orbs_aux%norb, orbs_aux%norbu*orbs_aux%norbup, nonzero, &
               orbs_aux%norbu*orbs_aux%norbup, nonzero, smat_extra(ispin), &
               nspin=input%nspin, geocode=geocode, cell_dim=cell_dim, norbup=orbs_aux%norbp, isorbu=orbs_aux%isorb, &
               store_index=input%store_index, on_which_atom=orbs_aux%onwhichatom, print_info=.false.)
          call f_free(nonzero)
          call deallocate_orbitals_data(orbs_aux)
    
      end do
    
      call f_release_routine()
    
    end subroutine init_sparse_matrix_for_KSorbs


    subroutine determine_sparsity_pattern_distance(orbs, lzd, astruct, cutoff, nnonzero, nonzero, smat_ref)
      use module_types
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none
    
      ! Calling arguments
      type(orbitals_data), intent(in) :: orbs
      type(local_zone_descriptors), intent(in) :: lzd
      type(atomic_structure), intent(in) :: astruct
      real(kind=8),dimension(lzd%nlr), intent(in) :: cutoff
      integer, intent(out) :: nnonzero
      integer, dimension(:,:), pointer,intent(out) :: nonzero
      type(sparse_matrix),intent(in),optional :: smat_ref !< reference sparsity pattern, in case the sparisty pattern to be calculated must be at least be as large as smat_ref
    
      ! Local variables
      logical :: overlap
      integer :: i1, i2, i3
      integer :: iorb, iiorb, ilr, iwa, itype, jjorb, jlr, jwa, jtype, ii
      integer :: ijs1, ije1, ijs2, ije2, ijs3, ije3, ind
      real(kind=8) :: tt, cut, xi, yi, zi, xj, yj, zj, x0, y0, z0
      logical :: perx, pery, perz, present_smat_ref

      call f_routine('determine_sparsity_pattern_distance')
      call timing(bigdft_mpi%iproc,'determinespars','ON')    

      present_smat_ref = present(smat_ref)
    
      ! periodicity in the three directions
      perx=(lzd%glr%geocode /= 'F')
      pery=(lzd%glr%geocode == 'P')
      perz=(lzd%glr%geocode /= 'F')
      ! For periodic boundary conditions, one has to check also in the neighboring
      ! cells (see in the loop below)
      if (perx) then
          ijs1 = -1
          ije1 = 1
      else
          ijs1 = 0
          ije1 = 0
      end if
      if (pery) then
          ijs2 = -1
          ije2 = 1
      else
          ijs2 = 0
          ije2 = 0
      end if
      if (perz) then
          ijs3 = -1
          ije3 = 1
      else
          ijs3 = 0
          ije3 = 0
      end if

      nnonzero=0
      do iorb=1,orbs%norbup
         iiorb=orbs%isorbu+iorb
         ilr=orbs%inwhichlocreg(iiorb)
         iwa=orbs%onwhichatom(iiorb)
         itype=astruct%iatype(iwa)
         xi=lzd%llr(ilr)%locregcenter(1)
         yi=lzd%llr(ilr)%locregcenter(2)
         zi=lzd%llr(ilr)%locregcenter(3)

         !$omp parallel default(none) &
         !$omp private(jjorb,ind,overlap,jlr,jwa,jtype,x0,y0,z0,cut,i3,i2,i1,zj,yj,xj,tt) &
         !$omp shared(xi,yi,zi,iiorb,ilr,orbs,nnonzero,lzd,cutoff,astruct) &
         !$omp shared(ijs3,ije3,ijs2,ije2,ijs1,ije1,present_smat_ref,smat_ref) 
         !$omp do schedule(guided) reduction(+:nnonzero)
         do jjorb=1,orbs%norbu
            if (present_smat_ref) then
                ind = matrixindex_in_compressed(smat_ref,jjorb,iiorb)
            else
                ind = 0
            end if
            if (ind>0) then
                ! There is an overlap in the reference sparsity pattern
                overlap = .true.
            else
                ! Check explicitly whether there is an overlap
                jlr=orbs%inwhichlocreg(jjorb)
                jwa=orbs%onwhichatom(jjorb)
                jtype=astruct%iatype(jwa)
                x0=lzd%llr(jlr)%locregcenter(1)
                y0=lzd%llr(jlr)%locregcenter(2)
                z0=lzd%llr(jlr)%locregcenter(3)
                cut = (cutoff(ilr)+cutoff(jlr))**2
                overlap = .false.
                do i3=ijs3,ije3!-1,1
                    zj=z0+i3*(lzd%glr%d%n3+1)*lzd%hgrids(3)
                    tt = (zi-zj)**2
                    if (tt<cut) then
                        ! For the z coordinate an overlap is possible, check the remaining directions
                        do i2=ijs2,ije2!-1,1
                            yj=y0+i2*(lzd%glr%d%n2+1)*lzd%hgrids(2)
                            tt = (zi-zj)**2 + (yi-yj)**2
                            if (tt<cut) then
                                ! For the y and z coordinate an overlap is possible, check the remaining direction
                                do i1=ijs1,ije1!-1,1
                                    xj=x0+i1*(lzd%glr%d%n1+1)*lzd%hgrids(1)
                                    tt = (zi-zj)**2 + (yi-yj)**2 + (xi-xj)**2
                                    if (tt<cut) then
                                        ! There is an overlap taking into account all three direction
                                        overlap=.true.
                                    end if
                                end do
                            end if
                        end do
                    end if
                end do
            end if
            if (overlap) then
               nnonzero=nnonzero+1
            end if
         end do
         !$omp end do
         !$omp end parallel
      end do

      !call mpiallred(nnonzero, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      nonzero = f_malloc_ptr((/2,nnonzero/),id='nonzero')
    
      ii=0
      !!do iorb=1,orbs%norbup
      !!   iiorb=orbs%isorbu+iorb
      !!   ilr=orbs%inwhichlocreg(iiorb)
      !!   iwa=orbs%onwhichatom(iiorb)
      !!   itype=astruct%iatype(iwa)
      !!   do jjorb=1,orbs%norbu
      !!      jlr=orbs%inwhichlocreg(jjorb)
      !!      jwa=orbs%onwhichatom(jjorb)
      !!      jtype=astruct%iatype(jwa)
      !!      tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 + &
      !!           (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 + &
      !!           (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
      !!      cut = cutoff(ilr)+cutoff(jlr)!+2.d0*incr
      !!      tt=sqrt(tt)
      !!      if (tt<=cut) then
      !!         ii=ii+1
      !!         nonzero(1,ii)=jjorb
      !!         nonzero(2,ii)=iiorb
      !!      end if
      !!   end do
      !!end do
      do iorb=1,orbs%norbup
         iiorb=orbs%isorbu+iorb
         ilr=orbs%inwhichlocreg(iiorb)
         iwa=orbs%onwhichatom(iiorb)
         itype=astruct%iatype(iwa)
         xi=lzd%llr(ilr)%locregcenter(1)
         yi=lzd%llr(ilr)%locregcenter(2)
         zi=lzd%llr(ilr)%locregcenter(3)

         do jjorb=1,orbs%norbu
            if (present_smat_ref) then
                ind = matrixindex_in_compressed(smat_ref,jjorb,iiorb)
            else
                ind = 0
            end if
            if (ind>0) then
                ! There is an overlap in the reference sparsity pattern
                overlap = .true.
            else
                ! Check explicitly whether there is an overlap
                jlr=orbs%inwhichlocreg(jjorb)
                jwa=orbs%onwhichatom(jjorb)
                jtype=astruct%iatype(jwa)
                x0=lzd%llr(jlr)%locregcenter(1)
                y0=lzd%llr(jlr)%locregcenter(2)
                z0=lzd%llr(jlr)%locregcenter(3)
                cut = (cutoff(ilr)+cutoff(jlr))**2
                overlap = .false.
                do i3=ijs3,ije3!-1,1
                    zj=z0+i3*(lzd%glr%d%n3+1)*lzd%hgrids(3)
                    tt = (zi-zj)**2
                    if (tt<cut) then
                        ! For the z coordinate an overlap is possible, check the remaining directions
                        do i2=ijs2,ije2!-1,1
                            yj=y0+i2*(lzd%glr%d%n2+1)*lzd%hgrids(2)
                            tt = (zi-zj)**2 + (yi-yj)**2
                            if (tt<cut) then
                                ! For the y and z coordinate an overlap is possible, check the remaining direction
                                do i1=ijs1,ije1!-1,1
                                    xj=x0+i1*(lzd%glr%d%n1+1)*lzd%hgrids(1)
                                    tt = (zi-zj)**2 + (yi-yj)**2 + (xi-xj)**2
                                    if (tt<cut) then
                                        ! There is an overlap taking into account all three direction
                                        overlap=.true.
                                    end if
                                end do
                            end if
                        end do
                    end if
                end do
            end if
            if (overlap) then
               ii=ii+1
               nonzero(1,ii)=jjorb
               nonzero(2,ii)=iiorb
            end if
         end do
      end do

      if (ii/=nnonzero) stop 'ii/=nnonzero'

      call timing(bigdft_mpi%iproc,'determinespars','OF')
      call f_release_routine()
    
    end subroutine determine_sparsity_pattern_distance


end module sparsematrix_wrappers
