module postprocessing_linear
  use public_enums
  implicit none

  private

  !> Public routines
  public :: loewdin_charge_analysis
  public :: loewdin_charge_analysis_core
  public :: support_function_multipoles
  public :: build_ks_orbitals
  public :: calculate_theta
  !public :: supportfunction_centers
  public :: projector_for_charge_analysis

  contains

    subroutine loewdin_charge_analysis(iproc,tmb,atoms,denspot, &
               calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap,blocksize,&
               ntheta, istheta, theta)
      use module_base
      use module_types
      use module_interfaces
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   DENSE_FULL, assignment(=), &
                                   matrices_null, allocate_matrices, deallocate_matrices
      use sparsematrix, only: compress_matrix, uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                              uncompress_matrix2
      use transposed_operations, only: calculate_overlap_transposed
      use matrix_operations, only: overlapPowerGeneral
      use yaml_output
      implicit none
      integer,intent(in) :: iproc, blocksize
      type(dft_wavefunction),intent(inout) :: tmb
      type(atoms_data),intent(in) :: atoms
      type(DFT_local_fields), intent(inout) :: denspot
      logical,intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
      integer,intent(in) :: meth_overlap
      integer,intent(in),optional :: ntheta, istheta
      real(kind=8),dimension(:,:),intent(in),optional :: theta ! must have dimension (atoms%astruct%nat,ndim_theta)
    
      !local variables
      !integer :: ifrag,ifrag_ref,isforb,jorb
      integer :: iorb,ierr
      real(kind=gp), allocatable, dimension(:,:,:) :: proj_mat
      real(kind=gp), allocatable, dimension(:,:) :: proj_ovrlp_half, weight_matrixp
      character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
      real(kind=gp) :: max_error, mean_error
      type(matrices),dimension(1) :: inv_ovrlp
    
      ! new variables
      integer :: iat
      real(kind=8),dimension(:,:),allocatable :: weight_matrix
      !real(kind=gp),dimension(:,:),pointer :: ovrlp
      real(kind=8) :: total_charge, total_net_charge
      real(kind=8),dimension(:),allocatable :: charge_per_atom
      !logical :: psit_c_associated, psit_f_associated
      logical :: optionals_present
    
    
      ! needs parallelizing/converting to sparse
      ! re-use overlap matrix if possible either before or after
    
      call f_routine(id='loewdin_charge_analysis')


      if (present(theta)) then
          if (.not.present(ntheta)) then
              call f_err_throw('ntheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else if (.not.present(istheta)) then
              call f_err_throw('istheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else
              optionals_present = .true.
          end if
          if (size(theta,1)/=atoms%astruct%nat) then
              call f_err_throw('wrong first dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (size(theta,2)/=ntheta) then
              call f_err_throw('wrong second dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      else
          optionals_present = .false.
      end if

    
      !inv_ovrlp(1) = matrices_null()
      !call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))
    
    
    
      if (calculate_overlap_matrix) then
         if(.not.tmb%can_use_transposed) then
             !!if(.not.associated(tmb%psit_c)) then
             !!    tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
             !!    psit_c_associated=.false.
             !!else
             !!    psit_c_associated=.true.
             !!end if
             !!if(.not.associated(tmb%psit_f)) then
             !!    tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
             !!    psit_f_associated=.false.
             !!else
             !!    psit_f_associated=.true.
             !!end if
             call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                  TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
             tmb%can_use_transposed=.true.
         end if
    
         call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
              tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
         !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
         ! This can then be deleted if the transition to the new type has been completed.
         !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr
    
    
         !!if (.not.psit_c_associated) then
         !!   call f_free_ptr(tmb%psit_c)
         !!   tmb%can_use_transposed=.false.
         !!end if
         !!if (.not.psit_f_associated) then
         !!   call f_free_ptr(tmb%psit_f)
         !!   tmb%can_use_transposed=.false.
         !!end if
      end if

      !call loewdin_charge_analysis_core(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%isorb, &
      !         tmb%orbs%norb_par, tmb%orbs%isorb_par, meth_overlap, tmb%linmat%s, tmb%linmat%l, atoms, &
      !         tmb%linmat%kernel_, tmb%linmat%ovrlp_)
      if (optionals_present) then
          call loewdin_charge_analysis_core(bigdft_mpi%iproc, bigdft_mpi%nproc, &
               tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, tmb%linmat%s%isfvctr, &
               tmb%linmat%s%nfvctr_par, tmb%linmat%s%isfvctr_par, &
               meth_overlap, blocksize, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_, &
               ntheta=ntheta, istheta=istheta, theta=theta)
      else
          call loewdin_charge_analysis_core(bigdft_mpi%iproc, bigdft_mpi%nproc, &
               tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, tmb%linmat%s%isfvctr, &
               tmb%linmat%s%nfvctr_par, tmb%linmat%s%isfvctr_par, &
               meth_overlap, blocksize, tmb%linmat%s, tmb%linmat%l, atoms, &
               tmb%linmat%kernel_, tmb%linmat%ovrlp_)
      end if
    
      call f_release_routine()
    
    
    end subroutine loewdin_charge_analysis




    subroutine loewdin_charge_analysis_core(iproc, nproc, norb, norbp, isorb, &
            norb_par, isorb_par, meth_overlap, blocksize, smats, smatl, atoms, kernel, ovrlp, &
               ntheta, istheta, theta)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   assignment(=), sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                                   DENSE_FULL, SPARSE_TASKGROUP, SPARSE_FULL, &
                                   deallocate_matrices, matrices_null, allocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: uncompress_matrix2, matrix_matrix_mult_wrapper, gather_matrix_from_taskgroups
      use matrix_operations, only: overlapPowerGeneral
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norb, norbp, isorb, meth_overlap, blocksize
      integer,dimension(0:nproc-1),intent(in) :: norb_par, isorb_par
      type(sparse_matrix),intent(inout) :: smats, smatl
      type(atoms_data),intent(in) :: atoms
      type(matrices),intent(inout) :: kernel
      type(matrices),intent(inout) :: ovrlp
      integer,intent(in),optional :: ntheta, istheta
      real(kind=8),dimension(:,:),intent(in),optional :: theta !dimension (atoms%astruct%nat,ntheta)

      ! Local variables
      integer :: ierr, iorb, iat, ind, ist, ishift, ispin, iiorb
      type(matrices),dimension(1) :: inv_ovrlp
      real(kind=8),dimension(:,:,:),allocatable :: proj_mat
      real(kind=8),dimension(:,:),allocatable :: weight_matrix, weight_matrixp, proj_ovrlp_half
      real(kind=8),dimension(:),allocatable :: charge_per_atom, proj_ovrlp_half_compr
      real(kind=8),dimension(:),allocatable :: weight_matrix_compr_tg, weight_matrix_compr
      real(kind=8) :: mean_error, max_error
      integer,parameter :: DENSE=101, SPARSE=102
      integer :: imode=SPARSE
      logical :: optionals_present

      call f_routine(id='loewdin_charge_analysis_core')

      if (present(theta)) then
          if (.not.present(ntheta)) then
              call f_err_throw('ntheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else if (.not.present(istheta)) then
              call f_err_throw('istheta not present',err_name='BIGDFT_RUNTIME_ERROR')
          else
              optionals_present = .true.
          end if
          if (size(theta,1)/=atoms%astruct%nat) then
              call f_err_throw('wrong first dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (size(theta,2)/=ntheta) then
              call f_err_throw('wrong second dimension of theta',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      else
          optionals_present = .false.
      end if


      if (imode==DENSE) then

          call f_err_throw('Dense mode is deprecated',err_name='BIGDT_RUNTIME_ERROR')

          inv_ovrlp(1) = matrices_null()
          call allocate_matrices(smatl, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))

          ovrlp%matrix = sparsematrix_malloc_ptr(smats, iaction=DENSE_FULL, id='ovrlp%matrix')
          call uncompress_matrix2(iproc, nproc, smats, &
               ovrlp%matrix_compr, ovrlp%matrix)
          call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), -1, &
               imode=2, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
               ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
               max_error=max_error, mean_error=mean_error)
          call f_free_ptr(ovrlp%matrix)
    
          ! optimize this to just change the matrix multiplication?
          proj_mat = sparsematrix_malloc0(smatl,iaction=DENSE_FULL,id='proj_mat')
    
          call uncompress_matrix2(iproc, nproc, smatl, kernel%matrix_compr, proj_mat)

          proj_ovrlp_half=f_malloc((/norb,norbp/),id='proj_ovrlp_half')
          if (norbp>0) then
             call dgemm('n', 'n', norb, norbp, &
                    norb, 1.d0, &
                    proj_mat(1,1,1), norb, &
                    inv_ovrlp(1)%matrix(1,isorb+1,1), norb, 0.d0, &
                    proj_ovrlp_half(1,1), norb)
          end if
          call f_free(proj_mat)
          weight_matrixp=f_malloc((/norb,norbp/), id='weight_matrixp')
          if (norbp>0) then
             call dgemm('n', 'n', norb, norbp, &
                  norb, 1.d0, &
                  inv_ovrlp(1)%matrix(1,1,1), norb, &
                  proj_ovrlp_half(1,1), norb, 0.d0, &
                  weight_matrixp(1,1), norb)
          end if
          call f_free(proj_ovrlp_half)
          weight_matrix=f_malloc((/norb,norb/), id='weight_matrix')
          if (nproc>1) then
             call mpi_allgatherv(weight_matrixp, norb*norbp, mpi_double_precision, weight_matrix, &
                  norb*norb_par(:), norb*isorb_par, &
                  mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          else
             call vcopy(norb*norb,weight_matrixp(1,1),1,weight_matrix(1,1),1)
          end if
          call f_free(weight_matrixp)
    
          charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')
    
          do iorb=1,norb
              iat=smats%on_which_atom(iorb)
              charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(iorb,iorb)
          end do
          if (iproc==0) then
              call write_partial_charges(atoms, charge_per_atom, .true.)
              call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
              call calculate_dipole(iproc, atoms, charge_per_atom)
              call calculate_quadropole(iproc, atoms, charge_per_atom)
              call yaml_sequence_close()
          end if
          !!call support_function_multipoles()
    
          call deallocate_matrices(inv_ovrlp(1))
          call f_free(charge_per_atom)
          call f_free(weight_matrix)

      else if (imode==SPARSE) then


          inv_ovrlp(1) = matrices_null()
          inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='inv_ovrlp(1)%matrix_compr')

          call overlapPowerGeneral(iproc, nproc, meth_overlap, 1, (/2/), blocksize, &
               imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
               ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
               max_error=max_error, mean_error=mean_error)
          call f_free_ptr(ovrlp%matrix)

          proj_ovrlp_half_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='proj_ovrlp_half_compr')
          weight_matrix_compr_tg = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='weight_matrix_compr_tg')
          do ispin=1,smatl%nspin
              ist = (ispin-1)*smatl%nvctrp_tg + 1
              if (norbp>0) then
                 call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                      kernel%matrix_compr(ist:), inv_ovrlp(1)%matrix_compr(ist:), proj_ovrlp_half_compr(ist:))
              end if
              if (norbp>0) then
                 call matrix_matrix_mult_wrapper(iproc, nproc, smatl, &
                      inv_ovrlp(1)%matrix_compr(ist:), proj_ovrlp_half_compr(ist:), weight_matrix_compr_tg(ist:))
              end if
          end do
          call f_free(proj_ovrlp_half_compr)
    
          call deallocate_matrices(inv_ovrlp(1))
          charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')

          ! Maybe this can be improved... not really necessary to gather the entire matrix
          weight_matrix_compr = sparsematrix_malloc0(smatl,iaction=SPARSE_FULL,id='weight_matrix_compr')
          call gather_matrix_from_taskgroups(iproc, nproc, smatl, weight_matrix_compr_tg, weight_matrix_compr)

          if (optionals_present) then
              do ispin=1,smatl%nspin
                  ishift = (ispin-1)*smatl%nvctr
                  do iorb=1,ntheta
                      iiorb = iorb + istheta
                      iiorb = modulo(iiorb-1,smatl%nfvctr)+1
                      ind = matrixindex_in_compressed(smatl, iorb, iorb)
                      !write(*,*) 'iorb, trace charge', iorb, weight_matrix_compr(ind)
                      do iat=1,atoms%astruct%nat
                          ind = ind + ishift
                          charge_per_atom(iat) = charge_per_atom(iat) + theta(iat,iorb)*weight_matrix_compr(ind)
                       end do
                  end do
              end do
          else 
              do ispin=1,smatl%nspin
                  ishift = (ispin-1)*smatl%nvctr
                  do iorb=1,norb
                      iiorb = modulo(iorb-1,smatl%nfvctr)+1
                      iat=smats%on_which_atom(iiorb)
                      ind = matrixindex_in_compressed(smatl, iorb, iorb)
                      ind = ind + ishift
                      !write(*,*) 'iorb, trace charge', iorb, weight_matrix_compr(ind)
                      charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix_compr(ind)
                  end do
              end do
          end if

          if (iproc==0) then
              call write_partial_charges(atoms, charge_per_atom, .true.)
              call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
              call calculate_dipole(iproc, atoms, charge_per_atom)
              call calculate_quadropole(iproc, atoms, charge_per_atom)
              call yaml_sequence_close()
          end if
    
          call f_free(charge_per_atom)
          call f_free(weight_matrix_compr_tg)
          call f_free(weight_matrix_compr)

      else
          call f_err_throw('wrong value for imode',err_name='BIGDFT_RUNTIME_ERROR')
      end if

      call f_release_routine()

    end subroutine loewdin_charge_analysis_core


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


    subroutine calculate_dipole(iproc, atoms, charge_per_atom)
      use module_base
      use module_types
      use yaml_output
      ! Calling arguments
      integer,intent(in) :: iproc
      type(atoms_data),intent(in) :: atoms
      real(kind=8),dimension(atoms%astruct%nat),intent(in) :: charge_per_atom
      ! Local variables
      integer :: iat
      real(kind=8),dimension(3) :: dipole_elec, dipole_cores, dipole_net
    
      dipole_cores(1:3)=0._gp
      do iat=1,atoms%astruct%nat
         dipole_cores(1:3)=dipole_cores(1:3)+atoms%nelpsp(atoms%astruct%iatype(iat))*atoms%astruct%rxyz(1:3,iat)
      end do
    
      dipole_elec=0.d0
      do iat=1,atoms%astruct%nat
          dipole_elec(1:3) = dipole_elec(1:3) -charge_per_atom(iat)*atoms%astruct%rxyz(1:3,iat)
      end do
    
      dipole_net=dipole_cores+dipole_elec
    
      if (iproc==0) then
          !!call yaml_map('core dipole', dipole_cores)
          !!call yaml_map('electronic dipole', dipole_elec)
          call yaml_map('net dipole', dipole_net,fmt='(es12.5)')
      end if
    
    end subroutine calculate_dipole


    subroutine calculate_quadropole(iproc, atoms, charge_per_atom)
      use module_base
      use module_types
      use yaml_output
      ! Calling arguments
      integer,intent(in) :: iproc
      type(atoms_data),intent(in) :: atoms
      real(kind=8),dimension(atoms%astruct%nat),intent(in) :: charge_per_atom
      ! Local variables
      real(kind=8),dimension(3,3) :: quadropole_elec, quadropole_cores, quadropole_net
      real(kind=8),dimension(3) :: charge_center_cores, charge_center_charge
      integer :: iat, i, j
      real(kind=8) :: delta_term, rj, ri, q, qtot
    
    
      ! charge center of the cores
      charge_center_cores(1:3)=0.d0
      qtot=0.d0
      do iat=1,atoms%astruct%nat
          q=atoms%nelpsp(atoms%astruct%iatype(iat))
          charge_center_cores(1:3) = charge_center_cores(1:3) + q*atoms%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_cores=charge_center_cores/qtot
    
    
      ! charge center of the charge
      charge_center_charge(1:3)=0.d0
      qtot=0.d0
      do iat=1,atoms%astruct%nat
          q=-charge_per_atom(iat)
          charge_center_charge(1:3) = charge_center_charge(1:3) + q*atoms%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_charge=charge_center_charge/qtot
    
    
      quadropole_cores(1:3,1:3)=0._gp
      do iat=1,atoms%astruct%nat
         q=atoms%nelpsp(atoms%astruct%iatype(iat))
         do i=1,3
             do j=1,3
                 if (i==j) then
                     delta_term = atoms%astruct%rxyz(1,iat)**2 + atoms%astruct%rxyz(2,iat)**2 + atoms%astruct%rxyz(3,iat)**2
                 else
                     delta_term=0.d0
                 end if
                 rj=atoms%astruct%rxyz(j,iat)
                 ri=atoms%astruct%rxyz(i,iat)
                 quadropole_cores(j,i) = quadropole_cores(j,i) + q*(3.d0*rj*ri-delta_term)
                 !!quadropole_cores(j,i) = quadropole_cores(j,i) + &
                 !!                        atoms%nelpsp(atoms%astruct%iatype(iat))* &
                 !!                          (3.d0*atoms%astruct%rxyz(j,iat)*atoms%astruct%rxyz(i,iat)-delta_term)
             end do
         end do
      end do
    
    
      quadropole_elec(1:3,1:3)=0._gp
      do iat=1,atoms%astruct%nat
         q=-charge_per_atom(iat)
         do i=1,3
             do j=1,3
                 if (i==j) then
                     delta_term = (atoms%astruct%rxyz(1,iat)+(charge_center_cores(1)-charge_center_charge(1)))**2 + &
                                  (atoms%astruct%rxyz(2,iat)+(charge_center_cores(2)-charge_center_charge(2)))**2 + &
                                  (atoms%astruct%rxyz(3,iat)+(charge_center_cores(3)-charge_center_charge(3)))**2
                 else
                     delta_term=0.d0
                 end if
                 rj=atoms%astruct%rxyz(j,iat)+(charge_center_cores(j)-charge_center_charge(j))
                 ri=atoms%astruct%rxyz(i,iat)+(charge_center_cores(i)-charge_center_charge(i))
                 quadropole_elec(j,i) = quadropole_elec(j,i) + q*(3.d0*rj*ri-delta_term)
                 !!quadropole_elec(j,i) = quadropole_elec(j,i) + &
                 !!                       -charge_per_atom(iat)* &
                 !!                         (3.d0*atoms%astruct%rxyz(j,iat)*atoms%astruct%rxyz(i,iat)-delta_term)
             end do
         end do
      end do
    
      quadropole_net=quadropole_cores+quadropole_elec
    
      if (iproc==0) then
          !!call yaml_sequence_open('core quadropole')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(quadropole_cores(i,1:3),fmt='(es12.5)')))
          !!end do
          !!call yaml_sequence_close()
    
          !!call yaml_sequence_open('electronic quadropole')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(quadropole_elec(i,1:3),fmt='(es12.5)')))
          !!end do
          !!call yaml_sequence_close()
    
          call yaml_sequence_open('net quadropole')
          do i=1,3
             call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3),fmt='(es12.5)')))
          end do
          call yaml_sequence(advance='no')
          call yaml_map('trace of quadropole matrix',&
               quadropole_net(1,1)+quadropole_net(2,2)+quadropole_net(3,3),fmt='(es12.2)')
          call yaml_sequence_close()
      end if
    
    end subroutine calculate_quadropole


    subroutine support_function_multipoles(iproc, tmb, atoms, denspot)
      use module_base
      use module_types
      use yaml_output
      
      ! Calling arguments
      integer,intent(in) :: iproc
      type(DFT_wavefunction),intent(in) :: tmb
      type(atoms_data),intent(in) :: atoms
      type(DFT_local_fields), intent(inout) :: denspot
    
      integer :: ist, istr, iorb, iiorb, ilr, i, iat
      real(kind=8),dimension(3) :: charge_center_elec
      real(kind=8),dimension(:),allocatable :: phir
      type(workarr_sumrho) :: w
      character(len=20) :: atomname
      real(kind=8),dimension(:,:),allocatable :: dipole_net
      real(kind=8),dimension(:,:,:),allocatable :: quadropole_net
    
      call f_routine(id='support_function_multipoles')
    
      phir = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir')
      dipole_net = f_malloc0((/3,tmb%orbs%norb/),id='dipole_net')
      quadropole_net = f_malloc0((/3,3,tmb%orbs%norb/),id='quadropole_net')
    
      !call to_zero(3*tmb%orbs%norb, dipole_net(1,1))
      !call to_zero(9*tmb%orbs%norb, quadropole_net(1,1,1))
    
      ist=1
      istr=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          iat=tmb%orbs%onwhichatom(iiorb)
          call initialize_work_arrays_sumrho(1,tmb%lzd%Llr(ilr),.true.,w)
          ! Transform the support function to real space
          call daub_to_isf(tmb%lzd%llr(ilr), w, tmb%psi(ist), phir(istr))
          call deallocate_work_arrays_sumrho(w)
          ! Calculate the charge center
          call charge_center(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
               denspot%dpbox%hgrids, phir(istr), charge_center_elec)
          !write(*,*) 'ilr, tmb%lzd%llr(ilr)%locregcenter', iat, tmb%lzd%llr(ilr)%locregcenter
          atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
          call calculate_multipoles(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
               denspot%dpbox%hgrids, phir(istr), charge_center_elec, tmb%lzd%llr(ilr)%locregcenter, &
               dipole_net(:,iiorb), quadropole_net(:,:,iiorb))
          !write(*,*) 'charge_center', charge_center_elec
          ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
          stop
      end if
    
    
      if (bigdft_mpi%nproc>1) then
          call mpiallred(dipole_net, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(quadropole_net, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      if (iproc==0) then
          call yaml_sequence_open('Support functions moments')
          do iorb=1,tmb%orbs%norb
              iat=tmb%orbs%onwhichatom(iorb)
              atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
              call yaml_sequence_open('number'//trim(yaml_toa(iorb))// &
                   ' (atom number ='//trim(yaml_toa(iat))//', type = '//trim(atomname)//')')
              call yaml_sequence(advance='no')
              call yaml_map('net dipole',dipole_net(:,iorb),fmt='(es18.10)')
              call yaml_sequence(advance='no')
              call yaml_sequence_open('net quadropole')
              do i=1,3
                 call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3,iorb),fmt='(es15.8)')))
              end do
              call yaml_sequence_close()
              call yaml_sequence_close()
          end do
          call yaml_sequence_close()
      end if
    
      call f_free(phir)
      call f_free(dipole_net)
      call f_free(quadropole_net)
      call f_release_routine()
    
    end subroutine support_function_multipoles


    subroutine calculate_multipoles(n1i, n2i, n3i, hgrids, phir, charge_center_elec, rxyz_center, &
               dipole_net, quadropole_net)
      use yaml_output
      implicit none
      ! Calling arguments
      integer,intent(in) :: n1i, n2i, n3i
      real(kind=8),dimension(3),intent(in) :: hgrids
      real(kind=8),dimension(n1i*n2i*n3i),intent(in) :: phir
      real(kind=8),dimension(3),intent(in) :: charge_center_elec, rxyz_center
      real(kind=8),dimension(3),intent(out) :: dipole_net
      real(kind=8),dimension(3,3),intent(out) :: quadropole_net
    
      integer :: i1, i2, i3, jj, iz, iy, ix, ii, i, j
      real(kind=8) :: q, x, y, z, qtot, ri, rj, delta_term
      real(kind=8),dimension(3) :: dipole_center, dipole_el
      real(kind=8),dimension(3,3) :: quadropole_center, quadropole_el
    
    
      !!call yaml_map('rxyz_center',rxyz_center,fmt='(es16.6)')
      !!call yaml_map('charge_center_elec',charge_center_elec,fmt='(es16.6)')
      !!call yaml_map('sum phir',sum(phir),fmt='(es16.6)')
    
      ! Dipole and quadropole of the support function
      dipole_el=0.d0
      quadropole_el=0.d0
      qtot=0.d0
      jj=0
      do i3=1,n3i
          do i2=1,n2i
              do i1=1,n1i
                  jj=jj+1
                  ! z component of point jj
                  iz=jj/(n2i*n1i)
                  ! Subtract the 'lower' xy layers
                  ii=jj-iz*(n2i*n1i)
                  ! y component of point jj
                  iy=ii/n1i
                  ! Subtract the 'lower' y rows
                  ii=ii-iy*n1i
                  ! x component
                  ix=ii
    
                  ! Shift the values due to the convolutions bounds
                  ix=ix-14
                  iy=iy-14
                  iz=iz-14
    
                  q = phir(jj)**2 * product(hgrids)
                  x = ix*hgrids(1) + (rxyz_center(1)-charge_center_elec(1))
                  y = iy*hgrids(2) + (rxyz_center(2)-charge_center_elec(2))
                  z = iz*hgrids(3) + (rxyz_center(3)-charge_center_elec(3))
    
                  ! Dipole part
                  dipole_el(1) = dipole_el(1) + q*x
                  dipole_el(2) = dipole_el(2) + q*y
                  dipole_el(3) = dipole_el(3) + q*z
                  qtot=qtot+q
    
                  ! Quadrupole part
                  do i=1,3
                      ri=get_r(i, x, y, z)
                      do j=1,3
                          rj=get_r(j, x, y, z)
                          if (i==j) then
                              delta_term = x**2 + y**2 + z**2
                          else
                              delta_term=0.d0
                          end if
                          quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
                      end do
                  end do
              end do
          end do
      end do

      !call yaml_map('rxyz_center',rxyz_center)
      !call yaml_map('charge_center_elec',charge_center_elec)
      !call yaml_map('qtot',qtot)
    
      ! Dipole of the center
      dipole_center(1) = -qtot*rxyz_center(1)
      dipole_center(2) = -qtot*rxyz_center(2)
      dipole_center(3) = -qtot*rxyz_center(3)
    
      ! Quadropole of the center
      quadropole_center=0.d0
      do i=1,3
          ri=rxyz_center(i)
          do j=1,3
              rj=rxyz_center(j)
              if (i==j) then
                  delta_term = rxyz_center(1)**2 + rxyz_center(2)**2 + rxyz_center(3)**2
              else
                  delta_term=0.d0
              end if
              quadropole_center(j,i) = quadropole_center(j,i) -qtot*(3.d0*rj*ri-delta_term)
          end do
      end do
    
      ! Net dipole and quadropole
      !call yaml_map('dipole_el',dipole_el)
      !call yaml_map('dipole_center',dipole_center)
      dipole_net = dipole_el + dipole_center
      quadropole_net = quadropole_el + quadropole_center
    
    
    !  call yaml_sequence_open(trim(yaml_toa(it))//'('//trim(atomname)//')')
    !  !call yaml_map('qtot',qtot)
    !  call yaml_sequence(advance='no')
    !  !call yaml_map('center dipole',dipole_center,fmt='(es16.6)')
    !  !call yaml_map('electronic dipole',dipole_el,fmt='(es18.10)')
    !  call yaml_map('net dipole',dipole_net,fmt='(es18.10)')
    !  call yaml_sequence(advance='no')
    !  !call yaml_sequence_open('center quadropole')
    !  !do i=1,3
    !  !   call yaml_sequence(trim(yaml_toa(quadropole_center(i,1:3),fmt='(es15.8)')))
    !  !end do
    !  !call yaml_sequence_close()
    !  !call yaml_sequence_open('electronic quadropole')
    !  !do i=1,3
    !  !   call yaml_sequence(trim(yaml_toa(quadropole_el(i,1:3),fmt='(es15.8)')))
    !  !end do
    !  !call yaml_sequence_close()
    !  call yaml_sequence_open('net quadropole')
    !  do i=1,3
    !     call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3),fmt='(es15.8)')))
    !  end do
    !  call yaml_sequence_close()
    !  call yaml_sequence_close()
    
      contains
    
        function get_r(i, x, y, z)
          integer,intent(in) :: i
          real(kind=8),intent(in) :: x, y, z
          real(kind=8) :: get_r
    
          select case (i)
          case (1)
              get_r=x
          case (2)
              get_r=y
          case (3)
              get_r=z
          case default
              stop 'wrong value of i'
          end select
        end function get_r
    
    end subroutine calculate_multipoles


    subroutine build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, order_taylor, &
               energy, energyDiff, energyold)
      use module_base
      use module_types
      use module_interfaces
      use communications_base, only: comms_cubic
      use communications_init, only: orbitals_communicators
      use communications, only: transpose_v, untranspose_v
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), SPARSE_FULL
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use yaml_output
      use rhopotential, only: updatePotential, sumrho_for_TMBs, corrections_for_negative_charge
      use locreg_operations, only: small_to_large_locreg
      implicit none
      
      ! Calling arguments
      integer:: iproc, nproc
      type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers), intent(inout) :: GPU
      type(energy_terms),intent(inout) :: energs
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      type(input_variables),intent(in) :: input
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(out) :: energy, energyDiff
      real(kind=8), intent(inout) :: energyold
    
      ! Local variables
      type(orbitals_data) :: orbs
      type(comms_cubic) :: comms
      real(gp) :: fnrm
      logical :: rho_negative
      integer :: infoCoeff, nvctrp, npsidim_global
      real(kind=8),dimension(:),pointer :: phi_global, phiwork_global
      real(kind=8),dimension(:),allocatable :: tmparr
      character(len=*),parameter :: subname='build_ks_orbitals'
      real(wp), dimension(:,:,:), pointer :: mom_vec_fake
      type(work_mpiaccumulate) :: energs_work
      integer,dimension(:,:),allocatable :: ioffset_isf
    
    
      nullify(mom_vec_fake)
    
      energs_work = work_mpiaccumulate_null()
      energs_work%ncount = 4
      call allocate_work_mpiaccumulate(energs_work)
    
    
      !debug
      !integer :: iorb, jorb, ist, jst, ierr, i
      !real(kind=8) :: ddot, tt
    
    
      ! Get the expansion coefficients
      ! Start with a "clean" density, i.e. without legacy from previous mixing steps
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
    
      tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
      call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
      call f_free(tmparr)
    
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
    
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
    
      tmb%can_use_transposed=.false.
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor,input%lin%max_inversion_error,input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap, energs_work, .false.)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
    
      if (bigdft_mpi%iproc ==0) then
         call write_eigenvalues_data(0.1d0,KSwfn%orbs,mom_vec_fake)
      end if
    
    
      !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
      !     max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      !call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
      !tmb%can_use_transposed=.false.
      !call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
      !energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      !energyDiff=energy-energyold
      !energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    
      ! Create communication arrays for support functions in the global box
      
      call nullify_orbitals_data(orbs)
      call copy_orbitals_data(tmb%orbs, orbs, subname)
      call orbitals_communicators(iproc, nproc, tmb%lzd%glr, orbs, comms)
    
    
      ! Transform the support functions to the global box
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
                         tmb%orbs%norb*comms%nvctr_par(iproc,0)*orbs%nspinor)
      phi_global = f_malloc_ptr(npsidim_global,id='phi_global')
      phiwork_global = f_malloc_ptr(npsidim_global,id='phiwork_global')
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
           tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), tmb%lzd, &
           KSwfn%lzd, tmb%orbs, tmb%psi, phi_global, to_global=.true.)
      call transpose_v(iproc, nproc, orbs, tmb%lzd%glr%wfd, comms, phi_global(1), phiwork_global(1))
    
    
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
      call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
                 tmb%orbs%norb, 0.d0, phiwork_global, nvctrp)
      
      call untranspose_v(iproc, nproc, KSwfn%orbs, tmb%lzd%glr%wfd, KSwfn%comms, phiwork_global(1), phi_global(1))  
    
      call f_free_ptr(phi_global)
    
      !!ist=1
      !!do iorb=1,KSwfn%orbs%norbp
      !!    do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!        write(800+iproc,*) iorb, i, phiwork_global(ist)
      !!        ist=ist+1
      !!    end do
      !!end do
    
    
      !!ierr=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!do i=1,KSwfn%orbs%norb*ierr
      !!    write(401,*) i, phiwork_global(i)
      !!end do
      !!write(*,*) 'GLOBAL DDOT',ddot(KSwfn%orbs%norb*ierr, phi_global, 1, phi_global, 1)
    
      !!do i=1,KSwfn%orbs%norb*ierr
      !!     tmb%psi(i)=phi_global(i)
      !!end do
      !!call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, tmb%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !!     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
    
      !!do i=1,KSwfn%orbs%norb*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)
      !!    write(600,'(i10,es16.7)') i, tmb%psi(i)
      !!end do
    
    
      !!write(*,*) 'iproc, input%output_wf_format',iproc, WF_FORMAT_PLAIN
      call writemywaves(iproc,trim(input%dir_output)//"wavefunction", WF_FORMAT_PLAIN, &
           KSwfn%orbs, KSwfn%Lzd%Glr%d%n1, KSwfn%Lzd%Glr%d%n2, KSwfn%Lzd%Glr%d%n3, &
           KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           at, rxyz, KSwfn%Lzd%Glr%wfd, phiwork_global)
    
      if (input%write_orbitals==2) then
          call write_orbital_density(iproc, .false., mod(input%lin%plotBasisFunctions,10), 'KSDens', &
               KSwfn%orbs%npsidim_orbs, phiwork_global, input, KSwfn%orbs, KSwfn%lzd, at, rxyz)
      end if
    
      if (input%wf_extent_analysis) then
          ioffset_isf = f_malloc0((/3,KSwfn%orbs%norbp/),id='ioffset_isf')
          call analyze_wavefunctions('Kohn Sham orbitals extent analysis', 'global', &
               KSwfn%lzd, KSwfn%orbs, KSwfn%orbs%npsidim_orbs, phiwork_global, ioffset_isf)
          call f_free(ioffset_isf)
      end if
    
    
    
       call f_free_ptr(phiwork_global)
       call deallocate_orbitals_data(orbs)
       call deallocate_comms_cubic(comms)
    
      ! To get consistent values of the energy and the Kohn-Sham residue with those
      ! which will be calculated by the cubic restart.
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
      call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
      call f_free(tmparr)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
      tmb%can_use_transposed=.false.
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor, input%lin%max_inversion_error, input%purification_quickreturn, &
           input%calculate_KS_residue, input%calculate_gap, energs_work, .false., updatekernel=.false.)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      energyDiff=energy-energyold
      energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    
      call deallocate_work_mpiaccumulate(energs_work)
    
    end subroutine build_ks_orbitals
    

    !TEMPORARY, to be cleaned/removed
    subroutine build_ks_orbitals_laura_tmp(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, order_taylor, &
               energy, energyDiff, energyold, npsidim_global, phiwork_global)
      use module_base
      use module_types
      use module_interfaces
      use communications_base, only: comms_cubic
      use communications_init, only: orbitals_communicators
      use communications, only: transpose_v, untranspose_v
      use sparsematrix_base, only: sparse_matrix
      use yaml_output
      use rhopotential, only: updatepotential, sumrho_for_TMBs, corrections_for_negative_charge
      use locreg_operations, only: small_to_large_locreg
      implicit none
      
      ! Calling arguments
      integer:: iproc, nproc
      type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers), intent(inout) :: GPU
      type(energy_terms),intent(inout) :: energs
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      type(input_variables),intent(in) :: input
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(out) :: energy, energyDiff
      real(kind=8), intent(inout) :: energyold
    integer, intent(in) :: npsidim_global
    real(kind=8),dimension(:),pointer :: phiwork_global
    
      ! Local variables
      type(orbitals_data) :: orbs
      type(comms_cubic) :: comms
      real(gp) :: fnrm
      logical :: rho_negative
      integer :: infoCoeff, nvctrp
      real(kind=8),dimension(:),pointer :: phi_global
      real(kind=8),dimension(:,:),allocatable :: coeffs_occs
      character(len=*),parameter :: subname='build_ks_orbitals'
      real(wp), dimension(:,:,:), pointer :: mom_vec_fake
      integer :: iorb, itmb
      type(work_mpiaccumulate) :: energs_work
    
      nullify(mom_vec_fake)
    
      energs_work = work_mpiaccumulate_null()
      energs_work%ncount = 4
      call allocate_work_mpiaccumulate(energs_work)
    
      !debug
      !integer :: iorb, jorb, ist, jst, ierr, i
      !real(kind=8) :: ddot, tt
    
    
      ! Get the expansion coefficients
      ! Start with a "clean" density, i.e. without legacy from previous mixing steps
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
    
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
    
      tmb%can_use_transposed=.false.
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor,input%lin%max_inversion_error,input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap,energs_work, .false.)
    
      if (bigdft_mpi%iproc ==0) then
         call write_eigenvalues_data(0.1d0,KSwfn%orbs,mom_vec_fake)
      end if
    
    
      !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
      !     max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      !call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
      !tmb%can_use_transposed=.false.
      !call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
      !energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      !energyDiff=energy-energyold
      !energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    
      ! Create communication arrays for support functions in the global box
      
      call nullify_orbitals_data(orbs)
      call copy_orbitals_data(tmb%orbs, orbs, subname)
      call orbitals_communicators(iproc, nproc, tmb%lzd%glr, orbs, comms)
    
    
      ! Transform the support functions to the global box
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      !npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
      !                   tmb%orbs%norb*comms%nvctr_par(iproc,0)*orbs%nspinor)
      phi_global = f_malloc_ptr(npsidim_global,id='phi_global')
      !phiwork_global = f_malloc_ptr(npsidim_global,id='phiwork_global')
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
           tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), tmb%lzd, &
           KSwfn%lzd, tmb%orbs, tmb%psi, phi_global, to_global=.true.)
      call transpose_v(iproc, nproc, orbs, tmb%lzd%glr%wfd, comms, phi_global(1), phiwork_global(1))
    
    !NOT PRINTING, 2xcorrect charge?!
    !print*,'ntmb,ntmbp,norb,norbp',tmb%orbs%norb,tmb%orbs%norbp,KSwfn%orbs%norb,KSwfn%orbs%norbp
      ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
      coeffs_occs=f_malloc0((/tmb%orbs%norb,KSwfn%orbs%norb/),id='coeffs_occs')
      !call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, KSwfn%orbs%occup(1), nvctrp, tmb%coeff(1,1), &
      !     tmb%orbs%norb, 0.d0, coeffs_occs, nvctrp)
                 do iorb=1,KSwfn%orbs%norbp
                    do itmb=1,tmb%orbs%norb
                        coeffs_occs(itmb,iorb) = sqrt(KSwfn%orbs%occup(KSwfn%orbs%isorb+iorb))*tmb%coeff(itmb,KSwfn%orbs%isorb+iorb)
    !print*,KSwfn%orbs%occup(KSwfn%orbs%isorb+iorb)
                    end do
                 end do
              call mpiallred(coeffs_occs, mpi_sum, comm=bigdft_mpi%mpi_comm)
    
      nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
      !call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
      call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, coeffs_occs(1,1), &
                 tmb%orbs%norb, 0.d0, phiwork_global, nvctrp)
      call f_free(coeffs_occs)  
    
      call untranspose_v(iproc, nproc, KSwfn%orbs, tmb%lzd%glr%wfd, KSwfn%comms, phiwork_global(1), phi_global(1))  
    
      call f_free_ptr(phi_global)
    
      !!ist=1
      !!do iorb=1,KSwfn%orbs%norbp
      !!    do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!        write(800+iproc,*) iorb, i, phiwork_global(ist)
      !!        ist=ist+1
      !!    end do
      !!end do
    
    
      !!ierr=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!do i=1,KSwfn%orbs%norb*ierr
      !!    write(401,*) i, phiwork_global(i)
      !!end do
      !!write(*,*) 'GLOBAL DDOT',ddot(KSwfn%orbs%norb*ierr, phi_global, 1, phi_global, 1)
    
      !!do i=1,KSwfn%orbs%norb*ierr
      !!     tmb%psi(i)=phi_global(i)
      !!end do
      !!call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, tmb%orbs, at, rxyz, denspot, GPU, infoCoeff, &
      !!     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
    
      !!do i=1,KSwfn%orbs%norb*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)
      !!    write(600,'(i10,es16.7)') i, tmb%psi(i)
      !!end do
    
    
      !!write(*,*) 'iproc, input%output_wf_format',iproc, WF_FORMAT_PLAIN
      !call writemywaves(iproc,trim(input%dir_output)//"wavefunction", WF_FORMAT_PLAIN, &
      !     KSwfn%orbs, KSwfn%Lzd%Glr%d%n1, KSwfn%Lzd%Glr%d%n2, KSwfn%Lzd%Glr%d%n3, &
      !     KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     at, rxyz, KSwfn%Lzd%Glr%wfd, phiwork_global)
    
       !call f_free_ptr(phiwork_global)
       call deallocate_orbitals_data(orbs)
       call deallocate_comms_cubic(comms)
    
    if (.false.) then
      ! To get consistent values of the energy and the Kohn-Sham residue with those
      ! which will be calculated by the cubic restart.
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
           max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
      tmb%can_use_transposed=.false.
      call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
           energs, nlpsp, input%SIC, tmb, fnrm, .true., .true., .false., .true., 0, 0, 0, 0, &
           order_taylor, input%lin%max_inversion_error, input%purification_quickreturn, &
           input%calculate_KS_residue, input%calculate_gap, energs_work, .false., updatekernel=.false.)
      energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      energyDiff=energy-energyold
      energyold=energy
    
      !!if(tmb%can_use_transposed) then
      !!    call f_free_ptr(tmb%psit_c)
      !!    call f_free_ptr(tmb%psit_f)
      !!end if
    end if
      call allocate_work_mpiaccumulate(energs_work)
    end subroutine build_ks_orbitals_laura_tmp



    subroutine calculate_theta(nat, rxyz, nphidim, phi, nphirdim, orbs, lzd, theta)
      use module_base
      use module_types, only: orbitals_data, local_zone_descriptors, workarr_sumrho
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: nat, nphidim, nphirdim
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      real(kind=8),dimension(nphidim),intent(in) :: phi
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(nat,orbs%norbp),intent(out) :: theta

      ! Local variables
      real(kind=8),dimension(:),allocatable :: psir
      type(workarr_sumrho) :: w
      integer :: ist, istr, iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, iat, iiat, l, m
      real(kind=8),dimension(3) :: com
      real(kind=8),dimension(-1:1) :: dipole
      real(kind=8) :: weight, tt, x, y, z, r2, hxh, hyh, hzh, q, qtot, monopole, r
      real(kind=8),parameter :: sigma2=0.1d0

      call f_zero(theta)

      ! Transform the support functions to real space
      psir = f_malloc(max(nphirdim,1),id='psir')
      ist=1
      istr=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call initialize_work_arrays_sumrho(1,lzd%Llr(ilr),.true.,w)
          call daub_to_isf(lzd%Llr(ilr), w, phi(ist), psir(istr))
          call deallocate_work_arrays_sumrho(w)
          !write(*,'(a,4i8,es16.6)') 'INITIAL: iproc, iiorb, n, istr, ddot', &
          !    iproc, iiorb, lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, &
          !    istr, ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1)
          !testarr(1,iiorb) = ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1) 
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=nphirdim+1) then
          call f_err_throw('ERROR on process '//adjustl(trim(yaml_toa(bigdft_mpi%iproc)))//': istr/=nphirdim+1', &
               err_name='BIGDFT_RUNTIME_ERROR')
          stop
      end if

      hxh = 0.5d0*lzd%hgrids(1)
      hyh = 0.5d0*lzd%hgrids(2)
      hzh = 0.5d0*lzd%hgrids(3)


      !! METHOD 1 ####################################################
      !istr = 1
      !do iorb=1,orbs%norbp
      !    iiorb=orbs%isorb+iorb
      !    ilr=orbs%inwhichlocreg(iiorb)
      !    write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
      !    com(1:3) = 0.d0
      !    weight = 0.d0
      !    do i3=1,lzd%llr(ilr)%d%n3i
      !        ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
      !        z = ii3*hzh
      !        do i2=1,lzd%llr(ilr)%d%n2i
      !            ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
      !            y = ii2*hyh
      !            do i1=1,lzd%llr(ilr)%d%n1i
      !                ii1 = lzd%llr(ilr)%nsi1 + i1 - 14 - 1
      !                x = ii1*hxh
      !                tt = psir(istr)**2
      !                com(1) = com(1) + x*tt
      !                com(2) = com(2) + y*tt
      !                com(3) = com(3) + z*tt
      !                weight = weight + tt
      !                istr = istr + 1
      !            end do
      !        end do
      !    end do
      !    call yaml_map('weight',weight)

      !    weight = 0.d0
      !    do iat=1,nat
      !        write(*,*) 'com, rxzy', com, rxyz(:,iat)
      !        r2 = (rxyz(1,iat)-com(1))**2 + (rxyz(2,iat)-com(2))**2 + (rxyz(3,iat)-com(3))**2
      !        tt = exp(-r2/(2*sigma2))
      !        theta(iat,iorb) = tt
      !        weight = weight + tt
      !    end do
      !    call yaml_map('weight',weight)
      !    tt = 1.d0/weight
      !    call dscal(nat, tt, theta(1,iorb), 1)

      !    call yaml_map('theta '//adjustl(trim(yaml_toa(iorb))),theta(:,iorb))

      !end do
      !! END METHOD 1 ################################################


      ! METHOD 2 ####################################################
      istr = 1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          iiat = orbs%onwhichatom(iiorb)
          !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
          com(1:3) = 0.d0
          dipole(:) = 0.d0
          weight = 0.d0
          qtot = 0.d0
          do i3=1,lzd%llr(ilr)%d%n3i
              ii3 = lzd%llr(ilr)%nsi3 + i3 - 14 - 1
              z = ii3*hzh
              !write(*,*) 'i3, ii3, z, d', i3, ii3, z, z-rxyz(3,iiat)
              do i2=1,lzd%llr(ilr)%d%n2i
                  ii2 = lzd%llr(ilr)%nsi2 + i2 - 14 - 1
                  y = ii2*hyh
                  do i1=1,lzd%llr(ilr)%d%n1i
                      ii1 = lzd%llr(ilr)%nsi1 + i1 - 14 - 1
                      x = ii1*hxh
                      q = psir(istr)**2!*hxh*hyh*hzh
                      com(1) = com(1) + x*q
                      com(2) = com(2) + y*q
                      com(3) = com(3) + z*q
                      qtot = qtot + q
                      !r2 = (rxyz(1,iiat)-x)**2 + (rxyz(2,iiat)-y)**2 + (rxyz(3,iiat)-z)**2
                      !tt = 1.d0/r2**4
                      !!r2 = (rxyz(1,iiat)-x)**2 + (rxyz(2,iiat)-y)**2 + (rxyz(3,iiat)-z)**2
                      !!if (r2/(2.d0*0.2d0)>300) then
                      !!    tt = 0.d0
                      !!else
                      !!    tt = safe_exp(-r2/(2.d0*0.2d0))
                      !!end if
                      !!if (tt*psir(istr)**2<1.d-300) then
                      !!    q =0.d0
                      !!else
                      !!    q = tt*psir(istr)**2
                      !!end if
                      do iat=1,nat
                          !write(*,*) 'com, rxzy', com, rxyz(:,iat)
                          !r2 = (rxyz(1,iat)-com(1))**2 + (rxyz(2,iat)-com(2))**2 + (rxyz(3,iat)-com(3))**2
                          r2 = (rxyz(1,iat)-x)**2 + (rxyz(2,iat)-y)**2 + (rxyz(3,iat)-z)**2
                          if (r2/(2*sigma2)>300) then
                              tt = 0.d0
                          else
                              theta(iat,iorb) = theta(iat,iorb) + q*safe_exp(-r2/(2*sigma2))
                          end if
                      end do
                      istr = istr + 1
                  end do
              end do
          end do
          dipole(-1) = com(2)
          dipole( 0) = com(3)
          dipole( 1) = com(1)
          com(1:3) = com(1:3)/qtot
          monopole = qtot - 1.d0
          dipole(:) = dipole(:) - qtot*rxyz(1:3,iiat)
          call yaml_map('monopole',monopole)
          call yaml_map('dipole',dipole)
          call yaml_map('qtot', qtot)
          call yaml_map('rxyz(..,iiat)', rxyz(1:3,iiat))
          call yaml_map('com', com)
          !if (iiat==1 .and. iiorb==5) then
          !    theta(:,iorb) = 0.d0
          !    theta(iiat,iorb) = 1.d0
          !end if

          do iat=1,nat
              theta(iat,iorb) = 0.d0
              !x = rxyz(1,iat) - rxyz(1,iiat)
              !y = rxyz(2,iat) - rxyz(2,iiat)
              !z = rxyz(3,iat) - rxyz(3,iiat)
              x = rxyz(1,iat) - com(1)
              y = rxyz(2,iat) - com(2)
              z = rxyz(3,iat) - com(3)
              r = sqrt( x**2 + y**2 + z**2 )
              do l=0,0!1
                  do m=-l,l
                      tt = spherical_harmonic(l, m, x, y, z)*gaussian(1.d0, r)
                      if(l==1) tt = tt * dipole(m)
                      write(*,*) 'iorb, iat, m, r, rxyz, tt', iorb, iat, m, r, rxyz(:,iat), tt
                      theta(iat,iorb) = theta(iat,iorb) + tt**2
                  end do
              end do
          end do

          tt = sum(theta(1:nat,iorb))
          !write(*,*) 'tt',tt
          tt = 1.d0/tt
          call dscal(nat, tt, theta(1,iorb), 1)

          call yaml_map('theta '//adjustl(trim(yaml_toa(iorb))),theta(:,iorb))

      end do
      ! END METHOD 2 ################################################

    end subroutine calculate_theta


    function gaussian(sigma, r) result(g)
      use module_base, only: pi => pi_param
      implicit none
      ! Calling arguments
      real(kind=8),intent(in) :: sigma, r
      real(kind=8) :: g
      
      g = exp(-r**2/(2.d0*sigma**2))
      g = g/sqrt(2.d0*pi*sigma**2)**3
      !g = g/(sigma**3*sqrt(2.d0*pi)**3)
    end function gaussian

    !> Calculates the real spherical harmonic for given values of l, m, x, y, z.
    function spherical_harmonic(l, m, x, y, z) result(sh)
      !use module_base, only: pi => pi_param
      use module_base, pi => pi_param
      implicit none
      ! Calling arguments
      integer,intent(in) :: l, m
      real(kind=8),intent(in) :: x, y, z
      real(kind=8) :: sh

      ! Local variables
      integer,parameter :: l_max=2
      real(kind=8) :: r, r2, rnorm

      if (l<0) call f_err_throw('l must be non-negative',err_name='BIGDFT_RUNTIME_ERROR')
      if (l>l_max) call f_err_throw('spherical harmonics only implemented up to l='//trim(yaml_toa(l_max)),&
          err_name='BIGDFT_RUNTIME_ERROR')
      if (abs(m)>l) call f_err_throw('abs(m) must not be larger than l',err_name='BIGDFT_RUNTIME_ERROR')


      ! Normalization for a sphere of radius rmax
      select case (l)
      case (0)
          sh = 0.5d0*sqrt(1/pi)
      case (1)
          r = sqrt(x**2+y**2+z**2)
          ! fix for small r (needs proper handling later...)
          if (r==0.d0) r=1.d-20
          select case (m)
          case (-1)
              sh = sqrt(3.d0/(4.d0*pi))*y/r
          case (0)
              sh = sqrt(3.d0/(4.d0*pi))*z/r
          case (1)
              sh = sqrt(3.d0/(4.d0*pi))*x/r
          end select
      case (2)
          r2 = x**2+y**2+z**2
          ! fix for small r2 (needs proper handling later...)
          if (r2==0.d0) r2=1.d-20
          select case (m)
          case (-2)
              sh = 0.5d0*sqrt(15.d0/pi)*x*y/r2
          case (-1)
              sh = 0.5d0*sqrt(15.d0/pi)*y*z/r2
          case (0)
              sh = 0.25d0*sqrt(5.d0/pi)*(-x**2-y**2+2*z**2)/r2
          case (1)
              sh = 0.5d0*sqrt(15.d0/pi)*z*x/r2
          case (2)
              sh = 0.25d0*sqrt(15.d0/pi)*(x**2-y**2)/r2
          end select
      end select

    end function spherical_harmonic



    subroutine supportfunction_centers(nat, rxyz, nphidim, phi, nphirdim, &
               norb, norbp, isorb, in_which_locreg, lzd, com)
      use module_base
      use module_types, only: local_zone_descriptors, workarr_sumrho
      use bounds, only: geocode_buffers
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: nat, nphidim, nphirdim, norb, norbp, isorb
      integer,dimension(norb),intent(in) :: in_which_locreg
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      real(kind=8),dimension(nphidim),intent(in) :: phi
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(3,norbp),intent(out) :: com

      ! Local variables
      real(kind=8),dimension(:),allocatable :: psir
      type(workarr_sumrho) :: w
      integer :: ist, istr, iorb, iiorb, ilr, i1, i2, i3, ii1, ii2, ii3, iat, iiat, l, m, nl1, nl2, nl3
      real(kind=8),dimension(-1:1) :: dipole
      real(kind=8) :: weight, tt, x, y, z, r2, hxh, hyh, hzh, q, qtot, monopole, r
      real(kind=8),parameter :: sigma2=0.1d0

      call f_routine(id='supportfunction_centers')

      ! Transform the support functions to real space
      psir = f_malloc(max(nphirdim,1),id='psir')
      ist=1
      istr=1
      do iorb=1,norbp
          iiorb=isorb+iorb
          ilr=in_which_locreg(iiorb)
          call initialize_work_arrays_sumrho(1,lzd%Llr(ilr),.true.,w)
          call daub_to_isf(lzd%Llr(ilr), w, phi(ist), psir(istr))
          call deallocate_work_arrays_sumrho(w)
          !write(*,'(a,4i8,es16.6)') 'INITIAL: iproc, iiorb, n, istr, ddot', &
          !    iproc, iiorb, lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, &
          !    istr, ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1)
          !testarr(1,iiorb) = ddot(lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i, psir(istr), 1, psir(istr), 1) 
          ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=nphirdim+1) then
          call f_err_throw('ERROR on process '//adjustl(trim(yaml_toa(bigdft_mpi%iproc)))//': istr/=nphirdim+1', &
               err_name='BIGDFT_RUNTIME_ERROR')
          stop
      end if

      hxh = 0.5d0*lzd%hgrids(1)
      hyh = 0.5d0*lzd%hgrids(2)
      hzh = 0.5d0*lzd%hgrids(3)

      istr = 1
      do iorb=1,norbp
          iiorb=isorb+iorb
          ilr=in_which_locreg(iiorb)
          call geocode_buffers(lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)
          !write(*,*) 'iorb, iiorb, ilr', iorb, iiorb, ilr
          com(1:3,iorb) = 0.d0
          weight = 0.d0
          do i3=1,lzd%llr(ilr)%d%n3i
              ii3 = lzd%llr(ilr)%nsi3 + i3 - nl3 - 1
              z = ii3*hzh
              do i2=1,lzd%llr(ilr)%d%n2i
                  ii2 = lzd%llr(ilr)%nsi2 + i2 - nl2 - 1
                  y = ii2*hyh
                  do i1=1,lzd%llr(ilr)%d%n1i
                      ii1 = lzd%llr(ilr)%nsi1 + i1 - nl1 - 1
                      x = ii1*hxh
                      tt = psir(istr)**2
                      com(1,iorb) = com(1,iorb) + x*tt
                      com(2,iorb) = com(2,iorb) + y*tt
                      com(3,iorb) = com(3,iorb) + z*tt
                      weight = weight + tt
                      istr = istr + 1
                  end do
              end do
          end do
          !call yaml_map('weight',weight)
          com(1:3,iorb) = com(1:3,iorb)/weight

      end do

      call f_free(psir)

      call f_release_routine()

    end subroutine supportfunction_centers






    subroutine projector_for_charge_analysis(at, smats, smatm, smatl, ovrlp_, ham_, kernel_, &
               lzd,  rxyz, norbs_per_type, centers_provided, &
               nphirdim, psi, orbs, com_)
      use module_base
      use module_types, only: local_zone_descriptors, orbitals_data
      use module_atoms, only: atoms_data
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   sparsematrix_malloc, sparsematrix_malloc0, &
                                   sparsematrix_malloc_ptr, sparsematrix_malloc0_ptr, &
                                   SPARSE_TASKGROUP, assignment(=), &
                                   matrices_null, deallocate_matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix
      use matrix_operations, only: overlapPowerGeneral, overlap_plus_minus_one_half_exact
      use yaml_output
      implicit none

      ! Calling arguments
      type(atoms_data),intent(in) :: at
      type(sparse_matrix),intent(inout) :: smats, smatl !< should be intent(in)...
      type(sparse_matrix),intent(in) :: smatm
      type(local_zone_descriptors),intent(in) :: lzd
      type(matrices),intent(inout) :: ovrlp_ !< should be intent(in)...
      type(matrices),intent(in) :: ham_, kernel_
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      integer,dimension(at%astruct%ntypes),intent(in) :: norbs_per_type
      logical,intent(in) :: centers_provided
      integer,intent(in),optional :: nphirdim
      real(kind=8),dimension(:),intent(in),optional :: psi
      type(orbitals_data),intent(in),optional :: orbs
      real(kind=8),dimension(:,:),pointer,intent(in),optional :: com_

      ! Local variables
      integer :: kat, iat, jat, i, j, ii, jj, icheck, n, indm, inds, ntot, ist, ind, iq, itype, ieval, ij, nmax, indl
      integer :: k, l, iatold, isat, natp, kkat, istot, ntotp, i1, i2, i3, is1, ie1, is2, ie2, is3, ie3, j1, j2, j3
      real(kind=8) :: r2, cutoff2, rr2, tt, ef, q, occ, max_error, mean_error
      real(kind=8) :: xi, xj, yi, yj, zi, zj, ttx, tty, ttz, xx, yy, zz, x, y, z
      real(kind=8),dimension(:),allocatable :: projector_compr
      real(kind=8),dimension(:,:),pointer :: com
      real(kind=8),dimension(:,:),allocatable :: ham, ovrlp, proj, ovrlp_tmp
      real(kind=8),dimension(:,:,:),allocatable :: coeff_all, ovrlp_onehalf_all
      integer,dimension(:,:,:,:),allocatable :: ilup
      real(kind=8),dimension(:),allocatable :: eval, eval_all, ovrlp_large, tmpmat1, tmpmat2, kerneltilde, charge_per_atom
      real(kind=8),dimension(:,:,:),allocatable :: tmpmat2d
      integer,dimension(:),allocatable :: id_all, n_all, itmparr
      real(kind=8),dimension(3) :: rr
      logical,dimension(:,:),allocatable :: neighbor
      type(matrices),dimension(1) :: ovrlp_onehalf_
      logical :: perx, pery, perz
      real(kind=8),parameter :: kT = 1.d-2


      call f_routine(id='projector_for_charge_analysis')

      ! Check the arguments
      if (centers_provided) then
          ! The centers of the support functions are already given
          if (.not.present(com_)) then
              call f_err_throw('com_ not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (size(com_,1)/=3) then
              call f_err_throw('wrong first dimension of com_',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (size(com_,2)/=smats%nfvctr) then
              call f_err_throw('wrong second dimension of com_',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          com => com_
      else
          ! Must calculate the centers of the support functions
          if (.not.present(nphirdim)) then
              call f_err_throw('nphirdim not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(psi)) then
              call f_err_throw('psi not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(orbs)) then
              call f_err_throw('orbs not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          com = f_malloc0_ptr((/3,orbs%norb/),id='com')
      end if

      ! Determine the periodicity...
      perx=(lzd%glr%geocode /= 'F')
      pery=(lzd%glr%geocode == 'P')
      perz=(lzd%glr%geocode /= 'F')
      if (perx) then
          is1 = -1
          ie1 = 1
      else
          is1 = 0
          ie1 = 0
      end if
      if (pery) then
          is2 = -1
          ie2 = 1
      else
          is2 = 0
          ie2 = 0
      end if
      if (perz) then
          is3 = -1
          ie3 = 1
      else
          is3 = 0
          ie3 = 0
      end if



      ! Parallelization over the number of atoms
      ii = at%astruct%nat/bigdft_mpi%nproc
      natp = ii
      jj = at%astruct%nat - bigdft_mpi%nproc*natp
      if (bigdft_mpi%iproc<jj) then
          natp = natp + 1
      end if
      isat = (bigdft_mpi%iproc)*ii + min(bigdft_mpi%iproc,jj)


      ! Determine the sum of the size of all submatrices (i.e. the total number of eigenvalues we will have)
      ! and the maximal value for one atom.
      neighbor = f_malloc((/smats%nfvctr,natp/),id='neighbor')
      neighbor(:,:) = .false.
      ntot = 0
      nmax = 0
      iatold = 0
      do kat=1,natp
          kkat = kat + isat
          n = 0
          do i=1,smats%nfvctr
               iat = smats%on_which_atom(i)
               ! Only do the following for the first TMB per atom
               if (iat==iatold) cycle
               iatold = iat
               if (iat==kkat) then
                   do j=1,smats%nfvctr
                       inds =  matrixindex_in_compressed(smats, i, j)
                       if (inds/=0) then
                          neighbor(j,kat) = .true.
                          n = n + 1
                       end if
                   end do
               end if
          end do
          ntot = ntot + n
          nmax = max(nmax, n)
      end do
      itmparr = f_malloc0(0.to.bigdft_mpi%nproc-1,id='itmparr')
      itmparr(bigdft_mpi%iproc) = ntot
      ntotp = ntot
      if (bigdft_mpi%nproc>1) then
          call mpiallred(itmparr, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(ntot, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(nmax, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
      end if
      istot = 0
      do i=0,bigdft_mpi%iproc-1
          istot = istot + itmparr(i)
      end do
      call f_free(itmparr)
      if (bigdft_mpi%iproc==0) then
          call yaml_map('maximal size of a submatrix',nmax)
      end if
      
      eval_all = f_malloc0(ntot,id='eval_all')
      id_all = f_malloc0(ntot,id='id_all')
      coeff_all = f_malloc((/nmax,nmax,natp/),id='coeff_all')
      ovrlp_onehalf_all = f_malloc((/nmax,nmax,natp/),id='ovrlp_onehalf_all')
      ilup = f_malloc((/2,nmax,nmax,natp/),id='ilup')
      n_all = f_malloc(natp,id='n_all')


      ! Centers of the support functions
      if (.not.centers_provided) then
          if (orbs%norb>0) then
              !call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, tmb%collcom_sr%ndimpsi_c, &
              !     orbs%norb, orbs%norbp, orbs%isorb, orbs%in_which_locreg, lzd, com(1:,orbs%isorb+1:))
              call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, nphirdim, &
                   orbs%norb, orbs%norbp, orbs%isorb, orbs%inwhichlocreg, lzd, com(1:,orbs%isorb+1:))
              if (bigdft_mpi%nproc>1) then
                  call mpiallred(com, mpi_sum, comm=bigdft_mpi%mpi_comm)
              end if
          end if
      end if


      ist = 0
      do kat=1,natp
          kkat = kat + isat

          ! Determine the size of the submatrix
          n = 0
          do j=1,smats%nfvctr
              if (neighbor(j,kat)) then
                  n = n + 1
              end if
          end do
          n_all(kat) = n


          ! Extract the submatrices
          ham = f_malloc0((/n,n/),id='ham')
          ovrlp = f_malloc0((/n,n/),id='ovrlp')
          proj = f_malloc0((/n,n/),id='proj')
          eval = f_malloc0((/n/),id='eval')
          icheck = 0
          ii = 0
          do i=1,smats%nfvctr
              if (neighbor(i,kat)) then
                  jj = 0
                  do j=1,smats%nfvctr
                      if (neighbor(j,kat)) then
                          icheck = icheck + 1
                          jj = jj + 1
                          if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                          inds =  matrixindex_in_compressed(smats, i, j)
                          if (inds>0) then
                              ovrlp(jj,ii) = ovrlp_%matrix_compr(inds)
                          else
                              ovrlp(jj,ii) = 0.d0
                          end if
                          indm =  matrixindex_in_compressed(smatm, i, j)
                          if (indm>0) then
                              ham(jj,ii) = ham_%matrix_compr(indm)
                          else
                              ham(jj,ii) = 0.d0
                          end if
                          !!! Search the midpoint between the two TMBs, taking into account the periodicity
                          !!r2 = huge(r2)
                          !!do i3=is3,ie3
                          !!    zi = com(3,i) + i3*at%astruct%cell_dim(3)
                          !!    zj = com(3,j) + i3*at%astruct%cell_dim(3)
                          !!    zz = 0.5d0*(zi+zj)
                          !!    ttz = (zi-zj)**2
                          !!    do i2=is2,ie2
                          !!        yi = com(2,i) + i2*at%astruct%cell_dim(2)
                          !!        yj = com(2,j) + i2*at%astruct%cell_dim(2)
                          !!        yy = 0.5d0*(yi+yj)
                          !!        tty = (yi-yj)**2
                          !!        do i1=is1,ie1
                          !!            xi = com(1,i) + i1*at%astruct%cell_dim(1)
                          !!            xj = com(1,j) + i1*at%astruct%cell_dim(1)
                          !!            xx = 0.5d0*(xi+xj)
                          !!            ttx = (xi-xj)**2
                          !!            tt = ttx + tty + ttz
                          !!            if (tt<r2) then
                          !!                rr(3) = zz
                          !!                rr(2) = yy
                          !!                rr(1) = xx
                          !!            end if
                          !!        end do
                          !!    end do
                          !!end do
                          !!! Determine the distance between the midpoint and the atom, taking into account the periodicity
                          !!rr2 = huge(rr2)
                          !!do i3=is3,ie3
                          !!    zz = rr(3) + i3*at%astruct%cell_dim(3)
                          !!    ttz = (zz-rxyz(3,kkat))**2
                          !!    do i2=is2,ie2
                          !!        yy = rr(2) + i2*at%astruct%cell_dim(2)
                          !!        tty = (yy-rxyz(2,kkat))**2
                          !!        do i1=is1,ie1
                          !!            xx = rr(1) + i1*at%astruct%cell_dim(1)
                          !!            ttx = (xx-rxyz(1,kkat))**2
                          !!            tt = ttx + tty + ttz
                          !!            if (tt<rr2) then
                          !!                rr2 = tt
                          !!            end if
                          !!        end do
                          !!    end do
                          !!end do
                          rr2 = huge(rr2)
                          do i3=is3,ie3
                              zi = com(3,i) + i3*at%astruct%cell_dim(3)
                              zj = com(3,j) !+ i3*at%astruct%cell_dim(3)
                              !zz = 0.5d0*(zi+zj)
                              zz = modulo(0.5d0*(zi+zj),at%astruct%cell_dim(3))
                              do i2=is2,ie2
                                  yi = com(2,i) + i2*at%astruct%cell_dim(2)
                                  yj = com(2,j) !+ i2*at%astruct%cell_dim(2)
                                  !yy = 0.5d0*(yi+yj)
                                  yy = modulo(0.5d0*(yi+yj),at%astruct%cell_dim(2))
                                  do i1=is1,ie1
                                      xi = com(1,i) + i1*at%astruct%cell_dim(1)
                                      xj = com(1,j) !+ i1*at%astruct%cell_dim(1)
                                      !xx = 0.5d0*(xi+xj)
                                      xx = modulo(0.5d0*(xi+xj),at%astruct%cell_dim(1))
                                      do j3=is3,ie3
                                          z = rxyz(3,kkat) + j3*at%astruct%cell_dim(3)
                                          ttz = (zz-z)**2
                                          do j2=is2,ie2
                                              y = rxyz(2,kkat) + j2*at%astruct%cell_dim(2)
                                              tty = (yy-y)**2
                                              do j1=is1,ie1
                                                  x = rxyz(1,kkat) + j1*at%astruct%cell_dim(1)
                                                  ttx = (xx-x)**2
                                                  tt = ttx + tty + ttz
                                                  if (tt<rr2) then
                                                      rr2 = tt
                                                  end if
                                              end do
                                          end do
                                      end do
                                  end do
                              end do
                          end do
                          !!rr2 = huge(rr2)
                          !!do i3=is3,ie3
                          !!    zz = rr(3) + i3*at%astruct%cell_dim(3)
                          !!    ttz = (zz-rxyz(3,kkat))**2
                          !!    do i2=is2,ie2
                          !!        yy = rr(2) + i2*at%astruct%cell_dim(2)
                          !!        tty = (yy-rxyz(2,kkat))**2
                          !!        do i1=is1,ie1
                          !!            xx = rr(1) + i1*at%astruct%cell_dim(1)
                          !!            ttx = (xx-rxyz(1,kkat))**2
                          !!            tt = ttx + tty + ttz
                          !!            if (tt<rr2) then
                          !!                rr2 = tt
                          !!            end if
                          !!        end do
                          !!    end do
                          !!end do
                          !rr(1) = 0.5d0*(com(1,i)+com(1,j))
                          !rr(2) = 0.5d0*(com(2,i)+com(2,j))
                          !rr(3) = 0.5d0*(com(3,i)+com(3,j))
                          !rr2 = (rr(1)-rxyz(1,kkat))**2 + (rr(2)-rxyz(2,kkat))**2 + (rr(3)-rxyz(3,kkat))**2
                          !write(*,*) 'kat, i, j, ii, jj, iat, jat, rr2', kat, i, j, ii, jj, rr2
                          ham(jj,ii) = ham(jj,ii) + 1.d0*(0.5d0)*rr2**3*ovrlp(jj,ii)
                          ilup(1,jj,ii,kat) = j
                          ilup(2,jj,ii,kat) = i
                      end if
                  end do
              end if
          end do
          if (icheck>n**2) then
              call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
                  &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
          end if

          ! Calculate ovrlp^1/2. The last argument is wrong, clean this.
          ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
          call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
          call overlap_plus_minus_one_half_exact(1, n, -1, .true., ovrlp_tmp, smats)
          do i=1,n
              call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
          end do
          call f_free(ovrlp_tmp)

          call diagonalizeHamiltonian2(bigdft_mpi%iproc, n, ham, ovrlp, eval)
          do i=1,n
              ii = ist + i
              eval_all(istot+ii) = eval(i)
              id_all(istot+ii) = kkat
              call vcopy(n, ham(1,i), 1, coeff_all(1,i,kat), 1)
          end do

          ist = ist + n

          call f_free(ham)
          call f_free(ovrlp)
          call f_free(proj)
          call f_free(eval)

      end do

      if (ist/=ntotp) call f_err_throw('ist/=ntotp',err_name='BIGDFT_RUNTIME_ERROR')

      if (bigdft_mpi%nproc>1) then
          call mpiallred(eval_all, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(id_all, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if


      ! Order the eigenvalues and IDs
      do i=1,ntot
          ! add i-1 since we are only searching in the subarray
          ind = minloc(eval_all(i:ntot),1) + (i-1)
          tt = eval_all(i)
          eval_all(i) = eval_all(ind)
          eval_all(ind) = tt
          ii = id_all(i)
          id_all(i) = id_all(ind)
          id_all(ind) = ii
      end do
    
    
      ! Calculate how many states should be included
      q = 0.d0
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          q = q + ceiling(0.5d0*real(at%nelpsp(itype),kind=8))
      end do
      iq = nint(q)
      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_open('Calculating projector for charge analysis')
          call yaml_map('number of states to be occupied (without smearing)',iq)
      end if
    
    
      ! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
      ef = eval_all(1)
      do
          ef = ef + 1.d-3
          occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
          if (abs(occ-1.d0)<1.d-8) exit
      end do
      if (bigdft_mpi%iproc==0) then
          call yaml_map('Pseudo Fermi level for occupations',ef)
      end if
    
    
      if (bigdft_mpi%iproc==0) then
          call yaml_sequence_open('ordered eigenvalues and occupations')
          ii = 0
          do i=1,ntot
              occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
              if (occ>1.d-100) then
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open(flow=.true.)
                  call yaml_map('eval',eval_all(i),fmt='(es13.4)')
                  call yaml_map('atom',id_all(i),fmt='(i5.5)')
                  call yaml_map('occ',occ,fmt='(1pg13.5e3)')
                  call yaml_mapping_close(advance='no')
                  call yaml_comment(trim(yaml_toa(i,fmt='(i5.5)')))
              else
                  ii = ii + 1
              end if
          end do
          if (ii>0) then
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('remaining states',ii)
              call yaml_map('occ','<1.d-100')
              call yaml_mapping_close()
          end if
          call yaml_sequence_close()
      end if
    
    
      projector_compr = sparsematrix_malloc0(iaction=SPARSE_TASKGROUP, smat=smatl, id='projector_compr')
      ! Calculate the projector. First for each single atom, then insert it into the big one.
      do kat=1,natp
          kkat = kat + isat
          n = n_all(kat)
          proj = f_malloc0((/n,n/),id='proj')
          ij = 0
          do ieval=1,ntot
              if (id_all(ieval)/=kkat) cycle
              ij = ij + 1
              occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
              do i=1,n
                  do j=1,n
                      proj(j,i) = proj(j,i) + occ*coeff_all(j,ij,kat)*coeff_all(i,ij,kat)
                  end do
             end do
          end do
          tmpmat2d = f_malloc((/n,n,1/),id='tmppmat2d')
          call gemm('n', 'n', n, n, n, 1.d0, proj(1,1), n, ovrlp_onehalf_all(1,1,kat), nmax, 0.d0, tmpmat2d(1,1,1), n)
          call gemm('n', 'n', n, n, n, 1.d0, ovrlp_onehalf_all(1,1,kat), nmax, tmpmat2d(1,1,1), n, 0.d0, proj(1,1), n)
          call f_free(tmpmat2d)
          do i=1,n
              do j=1,n
                   jj = ilup(1,j,i,kat)
                   ii = ilup(2,j,i,kat)
                   indl=matrixindex_in_compressed(smatl, ii, jj)
                   if (indl>0) then
                       ! Within the sparsity pattern
                       projector_compr(indl) = projector_compr(indl) + proj(j,i)
                   end if
              end do
          end do
          call f_free(proj)
      end do

      if (bigdft_mpi%nproc>1) then
          call mpiallred(projector_compr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      call f_free(coeff_all)
      call f_free(ilup)
      call f_free(n_all)
    
      tmpmat1 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat1')
      tmpmat2 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat2')
      kerneltilde = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='kerneltilde')


      ! Calculate S^1/2
      ovrlp_onehalf_(1) = matrices_null()
      ovrlp_onehalf_(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_onehalf_(1)%matrix_compr')
      call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, 1020, 1, (/2/), -1, &
            imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
            ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_(1), &
            check_accur=.true., max_error=max_error, mean_error=mean_error)
    
      ! Calculate S^1/2 * K * S^1/2  * P'
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           ovrlp_onehalf_(1)%matrix_compr, projector_compr, tmpmat1)
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           kernel_%matrix_compr, tmpmat1, tmpmat2)
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           ovrlp_onehalf_(1)%matrix_compr, tmpmat2, kerneltilde)

      call deallocate_matrices(ovrlp_onehalf_(1))

      call f_free(projector_compr)
      call f_free(tmpmat1)
      call f_free(tmpmat2)
    
      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_close()
      end if
    
    
      ! Calculate the partial traces
      charge_per_atom = f_malloc0(at%astruct%nat,id='charge_per_atom')
      call determine_atomic_charges(smatl, at%astruct%nat, kerneltilde, charge_per_atom)
    
      if (bigdft_mpi%iproc==0) then
          call write_partial_charges(at, charge_per_atom, write_gnuplot=.false.)
      end if

      call f_free(charge_per_atom)
      call f_free(kerneltilde)
      call f_free(neighbor)
      call f_free(eval_all)
      call f_free(id_all)
      call f_free(ovrlp_onehalf_all)
      if (.not.centers_provided) then
          call f_free_ptr(com)
      end if

      call f_release_routine()

  end subroutine projector_for_charge_analysis



  subroutine determine_atomic_charges(smat, nat, effective_kernel, charge_per_atom)
    use module_base
    use sparsematrix_base, only: sparse_matrix
    use sparsematrix_init, only: matrixindex_in_compressed
    implicit none
  
    ! Calling arguments
    type(sparse_matrix),intent(in) :: smat
    integer,intent(in) :: nat
    real(kind=8),dimension(smat%nvctr*smat%nspin),intent(in) :: effective_kernel
    real(kind=8),dimension(nat),intent(out) :: charge_per_atom
  
    ! Local variables
    integer :: ispin, ishift, iorb, iiorb, ind, iat
  
    call f_zero(charge_per_atom)
    
    do ispin=1,smat%nspin
        ishift = (ispin-1)*smat%nvctr
        do iorb=1,smat%nfvctr
            iiorb = modulo(iorb-1,smat%nfvctr)+1
            iat = smat%on_which_atom(iiorb)
            ind = matrixindex_in_compressed(smat, iorb, iorb)
            ind = ind + ishift
            !write(*,*) 'iorb, ind, val', iorb, ind, effective_kernel(ind)
            charge_per_atom(iat) = charge_per_atom(iat) + effective_kernel(ind)
        end do
    end do
  
  end subroutine determine_atomic_charges

end module postprocessing_linear
