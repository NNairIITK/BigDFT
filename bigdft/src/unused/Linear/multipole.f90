    subroutine projector_for_charge_analysis(smmd, smats, smatm, smatl, &
               ovrlp_, ham_, kernel_, rxyz, calculate_centers, write_output, ortho, mode, &
               lzd, nphirdim, psi, orbs, &
               multipoles, &
               natpx, isatx, nmaxx, nx, projx, neighborx, &
               rpower_matrix, only_sizes, psppar)
      use module_base
      use module_types, only: local_zone_descriptors, orbitals_data
      use module_atoms, only: atoms_data
      use sparsematrix_base, only: sparse_matrix, matrices, &
                                   sparsematrix_malloc, sparsematrix_malloc0, &
                                   sparsematrix_malloc_ptr, sparsematrix_malloc0_ptr, &
                                   SPARSE_TASKGROUP, assignment(=), &
                                   matrices_null, deallocate_matrices, &
                                   sparse_matrix_metadata
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix
      use sparsematrix_highlevel, only: trace_AB
      use matrix_operations, only: overlapPowerGeneral, overlap_plus_minus_one_half_exact
      use yaml_output
      use multipole_base, only: lmax
      use io, only: write_partial_charges
      implicit none

      ! Calling arguments
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smats, smatl
      type(sparse_matrix),intent(in) :: smatm
      type(matrices),intent(in) :: ovrlp_
      type(matrices),intent(in) :: ham_, kernel_
      real(kind=8),dimension(3,smmd%nat),intent(in) :: rxyz
      logical,intent(in) :: calculate_centers, write_output
      character(len=*),intent(in) :: ortho, mode
      type(local_zone_descriptors),intent(in),optional :: lzd
      integer,intent(in),optional :: nphirdim
      real(kind=8),dimension(:),intent(in),optional :: psi
      type(orbitals_data),intent(in),optional :: orbs
      real(kind=8),dimension(-lmax:lmax,0:lmax,1:smats%nfvctr),intent(in),optional :: multipoles
      integer,intent(out),optional :: natpx, isatx, nmaxx
      integer,dimension(:),pointer,intent(out),optional :: nx
      real(kind=8),dimension(:,:),allocatable,intent(out),optional :: projx
      logical,dimension(:,:),allocatable,intent(out),optional :: neighborx
      type(matrices),dimension(24),intent(in),optional :: rpower_matrix
      logical,intent(in),optional :: only_sizes
      real(kind=8),dimension(0:4,0:6,1:smmd%ntypes),intent(in),optional :: psppar

      ! Local variables
      integer :: kat, iat, jat, i, j, ii, jj, icheck, n, indm, inds, ntot, ist, ind, iq, itype, ieval, ij, nmax, indl, lwork
      integer :: k, l, iatold, isat, natp, kkat, istot, ntotp, i1, i2, i3, is1, ie1, is2, ie2, is3, ie3, j1, j2, j3, ikT, info
      integer :: ialpha, ilr, isshift, ilshift, ispin
      integer, dimension(1) :: power
      real(kind=8) :: r2, cutoff2, rr2, tt, ef, q, occ, max_error, mean_error, rr2i, rr2j, ttxi, ttyi, ttzi, ttxj, ttyj, ttzj
      real(kind=8) :: tti, ttj, charge_net, charge_total, rloc, charge, sigma2
      real(kind=8) :: xi, xj, yi, yj, zi, zj, ttx, tty, ttz, xx, yy, zz, x, y, z
      real(kind=8),dimension(:),allocatable :: work, occ_all, ef_atom
      real(kind=8),dimension(:,:),allocatable :: com
      real(kind=8),dimension(:,:),allocatable :: ham, ovrlp, proj, ovrlp_tmp, ovrlp_minusonehalf, kp, ktilde
      real(kind=8),dimension(:,:,:),allocatable :: coeff_all, ovrlp_onehalf_all, penaltymat
      integer,dimension(:,:,:,:),allocatable :: ilup
      real(kind=8),dimension(:),allocatable :: eval, eval_all, ovrlp_large, tmpmat1, tmpmat2, kerneltilde, charge_per_atom
      real(kind=8),dimension(:,:,:),allocatable :: tmpmat2d
      integer,dimension(:),allocatable :: id_all, n_all, itmparr
      real(kind=8),dimension(3) :: rr
      logical,dimension(:,:),allocatable :: neighbor
      integer,dimension(:),allocatable :: locregs_ID
      type(matrices),dimension(1) :: ovrlp_onehalf_
      logical :: perx, pery, perz, final, bound_low_ok, bound_up_ok, only_sizes_
      !real(kind=8),parameter :: kT = 5.d-2
      real(kind=8) :: kT, ttt, tr_KS
      !real(kind=8),parameter :: alpha = 5.d-1
      real(kind=8) :: alpha, alpha_up, alpha_low, convergence_criterion
      real(kind=8),dimension(:,:,:),allocatable :: multipoles_fake
      real(kind=8),dimension(:,:),allocatable :: penalty_matrices
      real(kind=8),dimension(:),allocatable :: alpha_calc
      !character(len=*),parameter :: mode='old'
      real(kind=8),dimension(3) :: target_charges
      character(len=*),parameter :: determine_ef='new'

      call f_routine(id='projector_for_charge_analysis')


      select case (trim(mode))
      case ('simple')
          ! everything ok
      case ('full')
          ! check the optional arguments
          if (.not.present(rpower_matrix)) then
              call f_err_throw('rpower_matrix not present')
          end if
          if (.not.present(orbs)) then
              call f_err_throw('orbs not present')
          end if
          if (.not.present(lzd)) then
              call f_err_throw('lzd not present')
          end if
          if (.not.present(psppar)) then
              call f_err_throw('psppar not present')
          end if
      case default
          call f_err_throw('wrong value of mode')
      end select

      only_sizes_ = .false.
      if (present(only_sizes)) only_sizes_ = only_sizes

      if (bigdft_mpi%iproc==0 .and. .not.only_sizes_) then
          call yaml_comment('Calculate projector for multipole analysis',hfill='~')
      end if

      if (present(natpx) .or. present(isatx) .or. present(nmaxx) .or. &
          present(projx) .or. present(neighborx) .or. present(nx)) then
          if (.not. (present(natpx) .and. present(isatx) .and. present(nmaxx) .and. &
              present(projx) .and. present(neighborx) .and.  present(nx))) then
              call f_err_throw('not all optional arguments present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      end if

      ! Determine the overall target charge, by calculating tr(KS). Only for ispin=1, the rest might not work
      ispin=1
      isshift=(ispin-1)*smats%nvctrp_tg
      ilshift=(ispin-1)*smatl%nvctrp_tg
      !tr_KS = trace_sparse(bigdft_mpi%iproc, bigdft_mpi%nproc, smats, smatl, &
      !       ovrlp_%matrix_compr(isshift+1:), &
      !       kernel_%matrix_compr(ilshift+1:), ispin)
      tr_KS = trace_AB(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, smats, smatl, ovrlp_, kernel_, ispin)



      kT = 1.d-2

      ! Convergence criterion: one million-th of the total charge
      tt = 0.d0
      do iat=1,smmd%nat
          tt = tt + real(smmd%nelpsp(smmd%iatype(iat)),kind=8)
      end do
      convergence_criterion = max(1.d-6*abs(tt),1.d-4)

      ! Check the arguments
      if (calculate_centers) then
      !!    ! The centers of the support functions are already given
      !!    if (.not.present(com_)) then
      !!        call f_err_throw('com_ not present',err_name='BIGDFT_RUNTIME_ERROR')
      !!    end if
      !!    if (size(com_,1)/=3) then
      !!        call f_err_throw('wrong first dimension of com_',err_name='BIGDFT_RUNTIME_ERROR')
      !!    end if
      !!    if (size(com_,2)/=smats%nfvctr) then
      !!        call f_err_throw('wrong second dimension of com_',err_name='BIGDFT_RUNTIME_ERROR')
      !!    end if
      !!    com => com_
      !!else
          ! Must calculate the centers of the support functions
          if (.not.present(lzd)) then
              call f_err_throw('lzd not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(nphirdim)) then
              call f_err_throw('nphirdim not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(psi)) then
              call f_err_throw('psi not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(orbs)) then
              call f_err_throw('orbs not present',err_name='BIGDFT_RUNTIME_ERROR')
          end if
      end if


      kerneltilde = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='kerneltilde')
      if (.not.only_sizes_) then
          if (ortho=='yes') then
              ! Calculate S^1/2
              ovrlp_onehalf_(1) = matrices_null()
              ovrlp_onehalf_(1)%matrix_compr = &
                  sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_onehalf_(1)%matrix_compr')
              power=2
              call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
                   1020, 1, power, -1, &
                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
                   ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_(1), &
                   check_accur=.false.)

              ! Calculate S^1/2 * K * S^1/2 = Ktilde
              tmpmat1 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat1')
              !tmpmat2 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat2')
              call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
                   kernel_%matrix_compr, ovrlp_onehalf_(1)%matrix_compr, tmpmat1)
              call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
                   ovrlp_onehalf_(1)%matrix_compr, tmpmat1, kerneltilde)
          else if (ortho=='no') then
              ovrlp_large = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_large')
              call transform_sparse_matrix(bigdft_mpi%iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
                   smat_in=ovrlp_%matrix_compr, lmat_out=ovrlp_large)
              call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
                   kernel_%matrix_compr, ovrlp_large, kerneltilde)
              call f_free(ovrlp_large)
          end if
      end if


      ! Parallelization over the number of atoms
      ii = smmd%nat/bigdft_mpi%nproc
      natp = ii
      jj = smmd%nat - bigdft_mpi%nproc*natp
      if (bigdft_mpi%iproc<jj) then
          natp = natp + 1
      end if
      isat = (bigdft_mpi%iproc)*ii + min(bigdft_mpi%iproc,jj)

      if (present(natpx)) natpx = natp
      if (present(isatx)) isatx = isat


      ! Determine the sum of the size of all submatrices (i.e. the total number of eigenvalues we will have)
      ! and the maximal value for one atom.
      neighbor = f_malloc((/smats%nfvctr,natp/),id='neighbor')
      neighbor(:,:) = .false.
      ntot = 0
      nmax = 0
      do kat=1,natp
          iatold = 0
          kkat = kat + isat
          n = 0
          do i=1,smats%nfvctr
               iat = smmd%on_which_atom(i)
               ! Only do the following for the first TMB per atom
               if (iat==iatold) cycle
               iatold = iat
               if (iat==kkat) then
                   do j=1,smats%nfvctr
                       inds =  matrixindex_in_compressed(smats, j, i)
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

      if (present(neighborx)) then
          neighborx = f_malloc((/smats%nfvctr,natpx/),id='neighborx')
          !call f_memcpy(src=neighbor,dest=neighborx)
          !neighborx = neighbor
          do iat=1,natpx
              do i=1,smats%nfvctr
                  neighborx(i,iat) = neighbor(i,iat)
              end do
          end do
      end if

      if (present(nx)) then
          nx = f_malloc_ptr(natpx,id='nx')
      end if


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
      
      eval_all = f_malloc0(ntot,id='eval_all')
      occ_all = f_malloc0(ntot,id='occ_all')
      id_all = f_malloc0(ntot,id='id_all')
      coeff_all = f_malloc((/nmax,nmax,natp/),id='coeff_all')
      ovrlp_onehalf_all = f_malloc((/nmax,nmax,natp/),id='ovrlp_onehalf_all')
      ovrlp_minusonehalf = f_malloc((/nmax,nmax/),id='ovrlp_minusonehalf')
      !ilup = f_malloc((/2,nmax,nmax,natp/),id='ilup')
      n_all = f_malloc(natp,id='n_all')

      penalty_matrices = f_malloc((/nmax**2,natp/),id='penalty_matrices')
      alpha_calc = f_malloc(natp,id='alpha_calc')


      ! Centers of the support functions
      com = f_malloc0((/3,smats%nfvctr/),id='com')
      if (calculate_centers) then
          if (orbs%norb>0) then
              !call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, tmb%collcom_sr%ndimpsi_c, &
              !     orbs%norb, orbs%norbp, orbs%isorb, orbs%in_which_locreg, lzd, com(1:,orbs%isorb+1:))
              call supportfunction_centers(smmd%nat, rxyz, size(psi), psi, nphirdim, &
                   orbs%norb, orbs%norbp, orbs%isorb, orbs%inwhichlocreg, lzd, com(1:,orbs%isorb+1:))
              if (bigdft_mpi%nproc>1) then
                  call mpiallred(com, mpi_sum, comm=bigdft_mpi%mpi_comm)
              end if
          end if
      else
          do i=1,smats%nfvctr
              iat = smmd%on_which_atom(i)
              com(1:3,i) = rxyz(1:3,iat)
          end do
      end if


      charge_per_atom = f_malloc0(smmd%nat,id='charge_per_atom')
      locregs_ID = f_malloc(smats%nfvctr,id='locregs_ID')


      ! Calculate how many states should be included
      q = 0.d0
      tt = 0.d0
      do iat=1,smmd%nat
          itype = smmd%iatype(iat)
          if (determine_ef=='old') then
              q = q + ceiling(0.5d0*real(smmd%nelpsp(itype),kind=8))
              tt = tt + real(smmd%nelpsp(itype),kind=8)
          else
              q = q + 0.5d0*real(smmd%nelpsp(itype),kind=8)
              tt = tt + real(smmd%nelpsp(itype),kind=8)
          end if
          !!if (at%nelpsp(itype)<=2) then
          !!    q = q + 1.d0
          !!else if (at%nelpsp(itype)<=8) then
          !!    q = q + 4.d0
          !!else if (at%nelpsp(itype)<=18) then
          !!    q = q + 9.d0
          !!else
          !!    call f_err_throw('strange electronic configuration')
          !!end if
      end do
      q = q + 0.5d0*(tr_KS - tt)
      iq = nint(q)
      if (bigdft_mpi%iproc==0 .and. .not.only_sizes_) then
          call yaml_mapping_open('Calculating projector for charge analysis')
          call yaml_map('convergence criterion',convergence_criterion)
          call yaml_map('maximal size of a submatrix',nmax)
          call yaml_sequence_open('Searching alpha for charge neutrality')
      end if

      ! Initial guess for the bisection bounds
      alpha_low = 1.d-3
      alpha_up = 1.d1
      bound_low_ok = .false.
      bound_up_ok = .false.

      eF = -1.d0

      ! Calculate the matrices <phi|r**x|phi>

      ef_atom = f_malloc(smmd%nat,id='ef_atom')

      alpha_loop: do ialpha=1,10000

          if (bigdft_mpi%iproc==0 .and. .not.only_sizes_) then
              call yaml_sequence(advance='no')
          end if

          if (.not.bound_low_ok) then
              alpha = alpha_low
          else if (.not.bound_up_ok) then
              alpha = alpha_up
          else
              alpha = 0.5d0*(alpha_low+alpha_up)
          end if

          !if (ialpha==0) alpha = 2.d-1
          if (ialpha==0) alpha = 2.d-1 !1.d-1

          charge_net = 0.d0
          call f_zero(eval_all)
          call f_zero(id_all)

          ist = 0
          atoms_loop: do kat=1,natp
              kkat = kat + isat
    
              ! Determine the size of the submatrix
              n = 0
              do j=1,smats%nfvctr
                  if (neighbor(j,kat)) then
                      n = n + 1
                      if (mode=='full') then
                          locregs_ID(n) = orbs%inwhichlocreg(j)
                          !write(*,*) 'j, n, ID', j, n, locregs_ID(n)
                      end if
                  end if
              end do
              n_all(kat) = n
              if (present(nx)) nx(kat) = n
    
    
              ! Extract the submatrices
              ham = f_malloc0((/n,n/),id='ham')
              ovrlp = f_malloc0((/n,n/),id='ovrlp')
              proj = f_malloc0((/n,n/),id='proj')
              penaltymat = f_malloc0((/n,n,24/),id='penaltymat')
              eval = f_malloc0((/n/),id='eval')
              call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1,kat), n_all(kat), ovrlp)!, ilup)
              call extract_matrix(smatm, ham_%matrix_compr, neighbor(1,kat), n_all(kat), ham)



              if (ortho=='yes' .and. .not.only_sizes_) then
                  ! Calculate ovrlp^1/2 and ovrlp^-1/2. The last argument is wrong, clean this.
                  ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                  call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                  ! Passing 0 as comm... not best practice
                  call overlap_plus_minus_one_half_exact(0, 1, 0, n, -1, .true., ovrlp_tmp, smats)
                  do i=1,n
                      call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
                  end do
                  call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                  ! Passing 0 as comm... not best practice
                  call overlap_plus_minus_one_half_exact(0, 1, 0, n, -1, .false., ovrlp_tmp, smats)
                  do i=1,n
                      call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_minusonehalf(1,i), 1)
                  end do
                  call f_free(ovrlp_tmp)
    
                  ! Calculate S^-1/2 * H * S^-1/2
                  tmpmat2d = f_malloc((/n,n,1/),id='tmppmat2d')
                  call gemm('n', 'n', n, n, n, 1.d0, ham(1,1), n, ovrlp_minusonehalf(1,1), nmax, 0.d0, tmpmat2d(1,1,1), n)
                  call gemm('n', 'n', n, n, n, 1.d0, ovrlp_minusonehalf(1,1), nmax, tmpmat2d(1,1,1), n, 0.d0, ham(1,1), n)
                  call f_free(tmpmat2d)
              end if



              ! Add the penalty term
              if (mode=='full') then
                  ! @ NEW #####################################################################
                  !write(*,*) 'call extract_matrix with penaltymat'
                  !write(*,*) 'orbs%onwhichatom',orbs%onwhichatom
                  !write(*,*) 'rxyz(:,kkat)',rxyz(:,kkat)
                  !write(*,*) 'HACK: SET ALPHA TO 0.5d0'
                  !alpha = 0.02d0
                  do i=1,24
                      call extract_matrix(smats, rpower_matrix(i)%matrix_compr, &
                          neighbor(1,kat), n_all(kat), penaltymat(1,1,i))
                  end do
                  !tt = sqrt(rxyz(1,kkat)**2+rxyz(2,kkat)**2+rxyz(3,kkat)**2)
                  tt = rxyz(1,kkat)**2 + rxyz(2,kkat)**2 + rxyz(3,kkat)**2
                  tmpmat2d = f_malloc0((/n,n,2/),id='tmppmat2d')
                  do i=1,n
                      do j=1,n
                          !if (i==j .and. orbs%onwhichatom(i)/=kkat) then
                          if (i==j) then
                              !!ttt = penaltymat(j,i,4) &
                              !!    - 2.d0*(rxyz(1,kkat)*penaltymat(j,i,1) &
                              !!          + rxyz(2,kkat)*penaltymat(j,i,2) &
                              !!          + rxyz(3,kkat)*penaltymat(j,i,3)) &
                              !!    + tt*ovrlp(j,i)
                              !!ttt = penaltymat(j,i,2) &
                              !!      + penaltymat(j,i,6) &
                              !!      + penaltymat(j,i,10) &
                              !!    - 2.d0*(rxyz(1,kkat)*penaltymat(j,i,1) &
                              !!          + rxyz(2,kkat)*penaltymat(j,i,5) &
                              !!          + rxyz(3,kkat)*penaltymat(j,i,9)) &
                              !!    + tt*ovrlp(j,i)
                              ttt = penaltymat(j,i,4) - 4.d0*rxyz(1,kkat)*penaltymat(j,i,3) &
                                    + 6.d0*rxyz(1,kkat)**2*penaltymat(j,i,2) - 4.d0*rxyz(1,kkat)**3*penaltymat(j,i,1) &
                                    + rxyz(1,kkat)**4*ovrlp(j,i) &
                                    + penaltymat(j,i,8) - 4.d0*rxyz(2,kkat)*penaltymat(j,i,7) &
                                    + 6.d0*rxyz(2,kkat)**2*penaltymat(j,i,6) - 4.d0*rxyz(2,kkat)**3*penaltymat(j,i,5) &
                                    + rxyz(2,kkat)**4*ovrlp(j,i) &
                                    + penaltymat(j,i,12) - 4.d0*rxyz(3,kkat)*penaltymat(j,i,11) &
                                    + 6.d0*rxyz(3,kkat)**2*penaltymat(j,i,10) - 4.d0*rxyz(3,kkat)**3*penaltymat(j,i,9) &
                                    + rxyz(3,kkat)**4*ovrlp(j,i) &
                                    + 2.d0*(penaltymat(j,i,16) &
                                            - 2.d0*rxyz(2,kkat)*penaltymat(j,i,14) &
                                            + rxyz(2,kkat)**2*penaltymat(j,i,2) &
                                            - 2.d0*rxyz(1,kkat)*penaltymat(j,i,15) &
                                            + 4.d0*rxyz(1,kkat)*rxyz(2,kkat)*penaltymat(j,i,13) &
                                            - 2.d0*rxyz(1,kkat)*rxyz(2,kkat)**2*penaltymat(j,i,1) &
                                            + rxyz(1,kkat)**2*penaltymat(j,i,6) &
                                            - 2.d0*rxyz(1,kkat)**2*rxyz(2,kkat)*penaltymat(j,i,5) &
                                            + rxyz(1,kkat)**2*rxyz(2,kkat)**2*ovrlp(j,i) &
                                            + penaltymat(j,i,20) &
                                            - 2.d0*rxyz(3,kkat)*penaltymat(j,i,18) &
                                            + rxyz(3,kkat)**2*penaltymat(j,i,2) &
                                            - 2.d0*rxyz(1,kkat)*penaltymat(j,i,19) &
                                            + 4.d0*rxyz(1,kkat)*rxyz(3,kkat)*penaltymat(j,i,17) &
                                            - 2.d0*rxyz(1,kkat)*rxyz(3,kkat)**2*penaltymat(j,i,1) &
                                            + rxyz(1,kkat)**2*penaltymat(j,i,10) &
                                            - 2.d0*rxyz(1,kkat)**2*rxyz(3,kkat)*penaltymat(j,i,9) &
                                            + rxyz(1,kkat)**2*rxyz(3,kkat)**2*ovrlp(j,i) &
                                            + penaltymat(j,i,24) &
                                            - 2.d0*rxyz(3,kkat)*penaltymat(j,i,22) &
                                            + rxyz(3,kkat)**2*penaltymat(j,i,6) &
                                            - 2.d0*rxyz(2,kkat)*penaltymat(j,i,23) &
                                            + 4.d0*rxyz(2,kkat)*rxyz(3,kkat)*penaltymat(j,i,21) &
                                            - 2.d0*rxyz(2,kkat)*rxyz(3,kkat)**2*penaltymat(j,i,5) &
                                            + rxyz(2,kkat)**2*penaltymat(j,i,10) &
                                            - 2.d0*rxyz(2,kkat)**2*rxyz(3,kkat)*penaltymat(j,i,9) &
                                            + rxyz(2,kkat)**2*rxyz(3,kkat)**2*ovrlp(j,i) )
                              itype = smmd%iatype(kkat)
                              rloc = psppar(0,0,itype)
                              ilr = locregs_ID(i)
                              !rloc = lzd%llr(ilr)%locrad
                              !write(*,*) 'kkat, itype, at%psppar(0,0,itype), rloc', kkat, itype, at%psppar(0,0,itype), rloc
                              !if (kkat==1) then
                                  ttt = ttt*alpha/rloc**4
                              !else
                              !    !if (i==1 .or. i==6) ttt=3.d0*ttt
                              !    ttt = 1.4641d0*alpha*ttt
                              !    !ttt = 5.0d0*alpha*ttt
                              !end if
                              !ttt = alpha*penaltymat(j,i,2) - alpha*2.d0*tt*penaltymat(j,i,1) + alpha*tt**2*ovrlp(j,i)
                              !ham(j,i) = ham(j,i) + ttt
                              tmpmat2d(j,i,1) = ttt
                              !write(*,*) 'kkat, j, i, owa(j), owa(i), alpha, tt, pm1, pm2, pm3, pm4, ovrlp, ttt', &
                              !            kkat, j, i, orbs%onwhichatom(j), orbs%onwhichatom(i), &
                              !            alpha, tt, penaltymat(j,i,1), penaltymat(j,i,2), &
                              !            penaltymat(j,i,3), penaltymat(j,i,4), &
                              !            ovrlp(j,i), ttt
                          end if
                      end do
                  end do
                  ! Additinal term
                  !call add_penalty_term(at%astruct%geocode, smats%nfvctr, neighbor(1:,kat), rxyz(1:,kkat), &
                  !     at%astruct%cell_dim, com, 10.d0*alpha, n, ovrlp, tmpmat2d)
                  if (ortho=='no') then
                      ! Calculate ovrlp^1/2. The last argument is wrong, clean this.
                      ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                      call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                      ! Passing 0 as comm... not best practice
                      call overlap_plus_minus_one_half_exact(0, 1, 0, n, -1, .true., ovrlp_tmp, smats)
                      do i=1,n
                          call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
                      end do
                      call f_free(ovrlp_tmp)

                      ! Calculate S^1/2 * penalty * S^1/2
                      call gemm('n', 'n', n, n, n, 1.d0, tmpmat2d(1,1,1), n, &
                           ovrlp_onehalf_all(1,1,kat), nmax, 0.d0, tmpmat2d(1,1,2), n)
                      call gemm('n', 'n', n, n, n, 1.d0, ovrlp_onehalf_all(1,1,kat), nmax, &
                           tmpmat2d(1,1,2), n, 0.d0, tmpmat2d(1,1,1), n)
                   end if
                   call axpy(n**2, 1.d0, tmpmat2d(1,1,1), 1, ham(1,1), 1)
                   call f_free(tmpmat2d)
                  ! @ END NEW #################################################################
              else if (mode=='simple') then
                  if (ortho=='yes') then
                      ! directly add the penalty terms to ham
                      call add_penalty_term(smmd%geocode, smats%nfvctr, neighbor(1,kat), rxyz(1,kkat), &
                           smmd%cell_dim, com, alpha, n, ovrlp, ham)
                   else if (ortho=='no') then
                          ! Calculate ovrlp^1/2. The last argument is wrong, clean this.
                          ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                          call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
                          ! Passing 0 as comm... not best practice
                          call overlap_plus_minus_one_half_exact(0, 1, 0, n, -1, .true., ovrlp_tmp, smats)
                          do i=1,n
                              call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
                          end do
                          call f_free(ovrlp_tmp)
                          ! Calculate the penaly term separately and then calculate S^1/2*penalty*S^1/2
                          tmpmat2d = f_malloc0((/n,n,2/),id='tmppmat2d')
                          call add_penalty_term(smmd%geocode, smats%nfvctr, neighbor(1,kat), rxyz(1,kkat), &
                               smmd%cell_dim, com, alpha, n, ovrlp, tmpmat2d(1,1,1))

                          ! Calculate S^1/2 * penalty * S^1/2
                          call gemm('n', 'n', n, n, n, 1.d0, tmpmat2d(1,1,1), n, &
                               ovrlp_onehalf_all(1,1,kat), nmax, 0.d0, tmpmat2d(1,1,2), n)
                          call gemm('n', 'n', n, n, n, 1.d0, ovrlp_onehalf_all(1,1,kat), nmax, &
                               tmpmat2d(1,1,2), n, 0.d0, tmpmat2d(1,1,1), n)
                          call axpy(n**2, 1.d0, tmpmat2d(1,1,1), 1, ham(1,1), 1)
                          call f_free(tmpmat2d)
                      end if
              else if (mode=='new') then
                  multipoles_fake = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.smats%nfvctr/),id='multipoles_fake')
                  multipoles_fake = 0.d0
                  multipoles_fake(0,0,:) = 1.d0
                  if (ialpha==1) then
                      if (present(multipoles)) then
                          write(*,*) 'call with multipoles'
                          call add_penalty_term_new(smmd%geocode, smmd%nat, smats%nfvctr, &
                               neighbor(1,kat), rxyz(1,kkat), smmd%on_which_atom, &
                               multipoles, smmd%cell_dim, com, alpha, n, ham, &
                               nmax, penalty_matrices(1,kat))
                      else
                          write(*,*) 'call with multipoles_fake'
                          call add_penalty_term_new(smmd%geocode, smmd%nat, smats%nfvctr, &
                               neighbor(1,kat), rxyz(1,kkat), smmd%on_which_atom, &
                               multipoles_fake, smmd%cell_dim, com, alpha, n, ham, &
                               nmax, penalty_matrices(1,kat))
                      end if
                      alpha_calc(kat) = alpha
                  else
                      tt = alpha/alpha_calc(kat)
                      !write(*,*) 'tt',tt
                      !do i=1,n
                      !    do j=1,n
                      !        write(*,*) 'i, j, penmat', i, j, penalty_matrices(j,i,kat)
                      !    end do
                      !end do
                      !ham(1:n,1:n) = ham(1:n,1:n) + tt*penalty_matrices(1:n,1:n,kat)
                      call axpy(n**2, tt, penalty_matrices(1,kat), 1, ham(1,1), 1)
                  end if
                  call f_free(multipoles_fake)
              end if

    
    
              !!call diagonalizeHamiltonian2(bigdft_mpi%iproc, n, ham, ovrlp, eval)
              lwork = 100*n
              work = f_malloc(lwork,id='work')
              if (ortho=='yes') then
                  call syev('v', 'l', n, ham(1,1), n, eval(1), work(1), lwork, info)
              else if (ortho=='no') then
                  ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
                  call f_memcpy(src=ovrlp,dest=ovrlp_tmp)
                  call sygv(1, 'v', 'l', n, ham(1,1), n, ovrlp_tmp(1,1), n, eval(1), work(1), lwork, info)
                  call f_free(ovrlp_tmp)
              end if
              call f_free(work)
              do i=1,n
                  ii = ist + i
                  eval_all(istot+ii) = eval(i)
                  id_all(istot+ii) = kkat
                  call vcopy(n, ham(1,i), 1, coeff_all(1,i,kat), 1)
              end do
    
              ist = ist + n
    
              call f_free(ham)
              call f_free(ovrlp)
              call f_free(penaltymat)
              call f_free(proj)
              call f_free(eval)

    
          end do atoms_loop

    
          if (ist/=ntotp) call f_err_throw('ist/=ntotp',err_name='BIGDFT_RUNTIME_ERROR')


          if (ialpha==1) then
              if (present(nmaxx)) nmaxx = maxval(n_all)
              if (present(projx)) projx = f_malloc((/nmaxx**2,natpx/),id='projx')
          end if

          ! This is maybe not very elegant in this way...
          if (only_sizes_) then
              call f_free(neighbor)
              call f_free(eval_all)
              call f_free(occ_all)
              call f_free(id_all)
              call f_free(coeff_all)
              call f_free(ovrlp_onehalf_all)
              call f_free(ovrlp_minusonehalf)
              !call f_free(ilup)
              call f_free(n_all)
              call f_free(penalty_matrices)
              call f_free(alpha_calc)
              call f_free(com)
              call f_free(charge_per_atom)
              call f_free(locregs_ID)
              call f_free(ef_atom)
              call f_free(kerneltilde)
              call f_release_routine()
              return
          end if
    
          if (bigdft_mpi%nproc>1) then
              call mpiallred(eval_all, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(id_all, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
    
    
          ! Order the eigenvalues and IDs
          call order_eigenvalues(ntot, eval_all, id_all)
        
        


    
    
              !ikT = ikT + 1

              !kT = 1.d-1*abs(eval_all(1))
              if (mode=='full') then
                  kT = max(1.d-1*abs(eval_all(iq)),1.d-1)
                  !kT = 1.d-2
                  !write(*,*) 'adjust kT to',kT
              end if
    
              call f_zero(charge_per_atom)
    
              ! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
              if (determine_ef=='old') then
                  if (ialpha>=0) then
                      ef = eval_all(1)
                      tt = ef
                      do
                          ef = ef + max(1.d-3,1.d-3*tt)
                          occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
                          if (abs(occ-1.d0)<1.d-8) exit
                      end do
                      !write(*,*) 'HACK: SET EF TO 0.0'
                      !ef = 0.0
                      !!eF = eF + 0.001d0
                      do ieval=1,ntot
                          occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
                          occ_all(ieval) = occ
                      end do
    !!$$$!target_charges(1) = 6.d0
    !!$$$!target_charges(2) = 1.d0
    !!$$$!target_charges(3) = 1.d0
    !!$$$kT = 1.d-2
    !!$$$!do itype=1,at%astruct%ntypes
    !!$$$!do iat=1,at%astruct%nat
    !!$$$call f_zero(ef_atom)
    !!$$$do kat=1,natp
    !!$$$    kkat = kat + isat
    !!$$$    itype=at%astruct%iatype(kkat)
    !!$$$    tt = 0.d0
    !!$$$    do ieval=1,ntot
    !!$$$        if(id_all(ieval)==itype) then
    !!$$$            tt = tt + 2.d0
    !!$$$        end if
    !!$$$        !if (tt>=target_charges(itype)-1.d-1) then
    !!$$$        if (tt>=at%nelpsp(itype)-1.d-2) then
    !!$$$            iq = ieval
    !!$$$            exit
    !!$$$        end if
    !!$$$    end do
    !!$$$
    !!$$$    ef = eval_all(1)
    !!$$$    do
    !!$$$        ef = ef + 1.d-3
    !!$$$        occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
    !!$$$        if (abs(occ-1.d0)<1.d-8) exit
    !!$$$    end do
    !!$$$    !if (bigdft_mpi%iproc==0) write(*,*) 'itype, ef', itype, ef
    !!$$$    ef_atom(kkat) = ef
    !!$$$    !write(*,*) 'HACK: SET EF TO 0.0'
    !!$$$    !ef = 0.0
    !!$$$    !!eF = eF + 0.001d0
    !!$$$    do ieval=1,ntot
    !!$$$        if (id_all(ieval)==itype) then
    !!$$$            occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
    !!$$$            !if (bigdft_mpi%iproc==0) write(*,*) 'itype, eval_all(ieval), ef, kT, occ', itype, itype, eval_all(ieval), ef, kT, occ
    !!$$$            occ_all(ieval) = occ
    !!$$$        end if
    !!$$$    end do
    !!$$$end do
    !!$$$if (bigdft_mpi%nproc>0) then
    !!$$$    call mpiallred(ef_atom, mpi_sum, comm=bigdft_mpi%mpi_comm)
    !!$$$end if
    !!$$$do ieval=1,ntot
    !!$$$ !!!$   if (target_charges(id_all(ieval))>1.d-10) then
    !!$$$ !!!$       occ = 1.d0
    !!$$$ !!!$       target_charges(id_all(ieval)) = target_charges(id_all(ieval)) - 2.d0
    !!$$$ !!!$   else
    !!$$$ !!!$       occ = 0.d0
    !!$$$ !!!$   end if
    !!$$$ !!!$   occ_all(ieval) = occ
    !!$$$    !occ = max(1.d0,0.5d0*target_charges(id_all(ieval)))
    !!$$$    !occ = real(ceiling(occ),kind=8)
    !!$$$    !target_charges(id_all(ieval)) = target_charges(id_all(ieval)) - 2.d0*occ
    !!$$$    !if (bigdft_mpi%iproc==0) write(*,*) 'ieval, eval, ID, occ', ieval, eval_all(ieval), id_all(ieval), occ_all(ieval)
    !!$$$end do
                      !!!if (bigdft_mpi%iproc==0) then
                      !!!    call yaml_sequence_close()
                      !!!    call yaml_map('number of states to be occupied (without smearing)',iq)
                      !!!    call yaml_map('Pseudo Fermi level for occupations',ef)
                      !!!    call yaml_sequence_open('ordered eigenvalues and occupations')
                      !!!    ii = 0
                      !!!    do i=1,ntot
                      !!!        !occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
                      !!!        occ = occ_all(i)
                      !!!        if (.true. .or. occ>1.d-100) then
                      !!!            call yaml_sequence(advance='no')
                      !!!            call yaml_mapping_open(flow=.true.)
                      !!!            call yaml_map('eval',eval_all(i),fmt='(es13.4)')
                      !!!            call yaml_map('atom',id_all(i),fmt='(i5.5)')
                      !!!            call yaml_map('occ',occ,fmt='(1pg13.5e3)')
                      !!!            call yaml_mapping_close(advance='no')
                      !!!            call yaml_comment(trim(yaml_toa(i,fmt='(i5.5)')))
                      !!!        else
                      !!!            ii = ii + 1
                      !!!        end if
                      !!!    end do
                      !!!    if (ii>0) then
                      !!!        call yaml_sequence(advance='no')
                      !!!        call yaml_mapping_open(flow=.true.)
                      !!!        call yaml_map('remaining states',ii)
                      !!!        call yaml_map('occ','<1.d-100')
                      !!!        call yaml_mapping_close()
                      !!!    end if
                      !!!    call yaml_sequence_close()
                      !!!end if
                  end if
              end if

!!ef = 0.d0
!!kT = 1.d-3
!!ef_loop: do        
!!    do ieval=1,ntot
!!        occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
!!        occ_all(ieval) = occ
!!    end do

if (determine_ef=='new') then
    ef = eval_all(iq)
    ii = 0
    do ieval=iq,ntot
        if (eval_all(ieval)<=eval_all(iq)+1.d0) then
            ii = ii + 1
        else
            exit
        end if
    end do
    sigma2 = min(1.d-1/real(ii,kind=8),2.d-3)
    do
        ef = ef + 1.d-4
        tt = 0.d0
        do ieval=iq,ntot
            if (eval_all(ieval)-eval_all(iq)>1.d0) exit
            tt = tt + safe_exp(-0.5d0*(ef-eval_all(ieval))**2/sigma2)
        end do
        if (tt<1.d-6) exit
    end do
    !if (bigdft_mpi%iproc==0) write(*,*) 'sigma2, ef',sigma2, ef
    do ieval=1,ntot
        if (eval_all(ieval)<=ef) then
            occ_all(ieval) = 1.d0
        else
            occ_all(ieval) = 0.d0
        end if
        !occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
        !occ_all(ieval) = occ
    end do
end if
        
              ! Calculate the projector. First for each single atom, then insert it into the big one.
              charge_total = 0.d0
              call f_zero(charge_per_atom)
              do kat=1,natp
                  kkat = kat + isat
                  n = n_all(kat)
                  proj = f_malloc0((/n,n/),id='proj')
                  call calculate_projector(n, ntot, nmax, kkat, id_all, eval_all, &
                       coeff_all(1,1,kat), occ_all, proj)
                  if (present(projx)) then
                      call vcopy(n**2, proj(1,1), 1, projx(1,kat), 1)
                  end if
    
    
                  !@ TEMPORARY ############################################
                  ! Extract ktilde
                  ktilde = f_malloc0((/n,n/),id='ktilde')
                  !if (ortho=='yes') then
                      call extract_matrix(smatl, kerneltilde, neighbor(1,kat), n_all(kat), ktilde)
                  !else if (ortho=='no') then
                  !    call extract_matrix(smatl, kernel_%matrix_compr, neighbor(1:,kat), n, nmax, ktilde)
                  !end if
                  kp = f_malloc((/n,n/),id='kp')
                  call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, proj(1,1), n, 0.d0, kp(1,1), n)
                  if (ortho=='no') then
                      call f_memcpy(src=kp,dest=ktilde)
                      ovrlp = f_malloc0((/n,n/),id='ovrlp')
                      call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1,kat), n_all(kat), ovrlp)
                      call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, ovrlp(1,1), n, 0.d0, kp(1,1), n)
                      call f_free(ovrlp)
                  end if
                  tt = 0
                  do i=1,n
                      tt = tt + kp(i,i)
                  end do
                  charge_per_atom(kkat) = tt
                  charge_total = charge_total + tt
                  call f_free(proj)
                  call f_free(ktilde)
                  call f_free(kp)
    
    
              end do
    
    

          if (bigdft_mpi%nproc>1) then
              call mpiallred(charge_per_atom, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(charge_total, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          if (bigdft_mpi%iproc==0) then
              !do iat=1,at%astruct%nat
              !    write(*,*) 'iat, cpa',iat,charge_per_atom(iat)
              !end do
              !write(*,*) 'charge_total',charge_total
          end if
          !charge_net = 0.d0
          charge = 0.d0
          do iat=1,smmd%nat
              !charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
              charge = charge - charge_per_atom(iat)
          end do
          charge_net = charge + tr_KS
          if (bigdft_mpi%iproc==0) then
              !write(*,*) 'net charge', charge_net
              call yaml_mapping_open(flow=.true.)
              call yaml_map('alpha',alpha,fmt='(es12.4)')
              call yaml_map('charge',charge,fmt='(es12.4)')
              call yaml_map('net charge',charge_net,fmt='(es12.4)')
              call yaml_map('bisection bounds ok',(/bound_low_ok,bound_up_ok/))
              if (determine_ef=='old') then
                  call yaml_map('kT',kT,fmt='(es12.4)')
              else
                  call yaml_map('sigma2',sigma2,fmt='(es12.4)')
              end if
              call yaml_map('eF',ef,fmt='(es12.4)')
              call yaml_mapping_close()
          end if

!!if (bigdft_mpi%iproc==0) then
!!call yaml_sequence_open('ordered eigenvalues and occupations')
!!ii = 0
!!do i=1,ntot
!!    !occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
!!    occ = occ_all(i)
!!    if (occ>1.d-100) then
!!        call yaml_sequence(advance='no')
!!        call yaml_mapping_open(flow=.true.)
!!        call yaml_map('eval',eval_all(i),fmt='(es13.4)')
!!        call yaml_map('atom',id_all(i),fmt='(i5.5)')
!!        call yaml_map('eF',ef_atom(id_all(i)),fmt='(es13.4)')
!!        call yaml_map('occ',occ,fmt='(1pg13.5e3)')
!!        call yaml_mapping_close(advance='no')
!!        call yaml_comment(trim(yaml_toa(i,fmt='(i5.5)')))
!!    else
!!        ii = ii + 1
!!    end if
!!end do
!!if (ii>0) then
!!    call yaml_sequence(advance='no')
!!    call yaml_mapping_open(flow=.true.)
!!    call yaml_map('remaining states',ii)
!!    call yaml_map('occ','<1.d-100')
!!    call yaml_mapping_close()
!!end if
!!call yaml_sequence_close()
!!end if


          if (abs(charge_net)<convergence_criterion .or. ialpha==10000) then
          !if (.true.) then
              if (bigdft_mpi%iproc==0) then
                  call yaml_sequence_close()
                  call yaml_map('number of states to be occupied (without smearing)',iq)
                  call yaml_map('Pseudo Fermi level for occupations',ef)
                  call yaml_sequence_open('ordered eigenvalues and occupations')
                  ii = 0
                  do i=1,ntot
                      !occ = 1.d0/(1.d0+safe_exp( (eval_all(i)-ef)*(1.d0/kT) ) )
                      occ = occ_all(i)
                      if (occ<1.d-100) then
                          ii = ii + 1
                      end if
                      !if (occ>1.d-100 .or. .true.) then
                      if (ii<=10) then
                          call yaml_sequence(advance='no')
                          call yaml_mapping_open(flow=.true.)
                          call yaml_map('eval',eval_all(i),fmt='(es13.4)')
                          call yaml_map('atom',id_all(i),fmt='(i5.5)')
                          !call yaml_map('eF',ef_atom(id_all(i)),fmt='(es13.4)')
                          call yaml_map('occ',occ,fmt='(1pg13.5e3)')
                          call yaml_mapping_close(advance='no')
                          call yaml_comment(trim(yaml_toa(i,fmt='(i5.5)')))
                      !!else
                      !!    ii = ii + 1
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
              exit alpha_loop
          end if

!!    if (charge_net<0.d0) then
!!        ef = ef - 1.d-4
!!    else
!!        ef = ef + 1.d-4
!!    end if
!!
!!end do ef_loop

          ! If we are still searching the boundaries for the bisection...
          if (.not.bound_low_ok) then
              if (charge_net<0.d0) then
                  ! this is a lower bound
                  alpha_low = alpha
                  bound_low_ok = .true.
              else
                  alpha_low = 0.5d0*alpha
              end if
              cycle alpha_loop
          else if (.not.bound_up_ok) then
              if (charge_net>0.d0) then
                  ! this is an upper bound
                  alpha_up = alpha
                  bound_up_ok = .true.
              else
                  alpha_up = 2.0d0*alpha
              end if
              cycle alpha_loop
          end if

          if (charge_net>0.d0) then
              ! Too few electrons, i.e. confinement should be smaller
              !alpha = alpha*0.80
              alpha_up = alpha
          else if (charge_net<0.d0) then
              ! Too many electrons, i.e. confinement should be larger
              !alpha = alpha*1.2
              alpha_low = alpha
          end if



      end do alpha_loop

      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_close()
          call yaml_comment('Projector calculation finished',hfill='~')
      end if


      !call f_free(tmpmat2)
      if (ortho=='yes') then
          call f_free(tmpmat1)
          call deallocate_matrices(ovrlp_onehalf_(1))
      end if
          call f_free(kerneltilde)
      call f_free(coeff_all)
      !call f_free(ilup)
      call f_free(n_all)
      call f_free(ovrlp_minusonehalf)
      call f_free(penalty_matrices)
      call f_free(alpha_calc)
    
      !if (bigdft_mpi%iproc==0) then
      !    call yaml_mapping_close()
      !end if
    
    
    
      if (write_output .and. bigdft_mpi%iproc==0) then
          call f_err_throw('writing the data here is deprecated')
          !call write_partial_charges(at, charge_per_atom, write_gnuplot=.true.)
      end if

      call f_free(ef_atom)
      call f_free(locregs_ID)
      call f_free(charge_per_atom)
      call f_free(neighbor)
      call f_free(eval_all)
      call f_free(occ_all)
      call f_free(id_all)
      call f_free(ovrlp_onehalf_all)
      call f_free(com)

      call f_release_routine()

  end subroutine projector_for_charge_analysis

  !!subroutine multipole_analysis_driver(iproc, nproc, lmax, ixc, smmd, smats, smatm, smatl, &
  !!           ovrlp, ham, kernel, rxyz, method, do_ortho, projectormode, &
  !!           calculate_multipole_matrices, do_check, &
  !!           nphi, lphi, nphir, hgrids, orbs, collcom, collcom_sr, &
  !!           lzd, at, denspot, orthpar, shift, multipole_matrix_in, ice_obj)
  !!  use module_base
  !!  use module_types, only: orbitals_data, comms_linear, local_zone_descriptors, orthon_data, DFT_local_fields, comms_linear
  !!  use sparsematrix_base, only: sparse_matrix, matrices, sparsematrix_malloc0, assignment(=), &
  !!                               sparsematrix_malloc, matrices_null, sparsematrix_malloc_ptr, deallocate_matrices, &
  !!                               SPARSE_TASKGROUP, sparse_matrix_metadata
  !!  use sparsematrix, only: matrix_matrix_mult_wrapper, transform_sparse_matrix
  !!  use communications, only: transpose_localized
  !!  use orthonormalization, only: orthonormalizelocalized,overlap_matrix
  !!  use module_atoms, only: atoms_data
  !!  use yaml_output
  !!  use multipole_base, only: external_potential_descriptors_null, multipole_set_null, multipole_null, &
  !!       deallocate_external_potential_descriptors
  !!  use orbitalbasis
  !!  use matrix_operations, only: overlapPowerGeneral
  !!  !use Poisson_Solver, only: H_potential
  !!  use Poisson_Solver, except_dp => dp, except_gp => gp
  !!  use foe_base, only: foe_data
  !!  use box
  !!  implicit none
  !!  ! Calling arguments
  !!  integer,intent(in) :: iproc, nproc, lmax, ixc
  !!  type(sparse_matrix_metadata),intent(in) :: smmd
  !!  type(sparse_matrix),intent(in) :: smats
  !!  type(sparse_matrix),intent(in) :: smatm
  !!  type(sparse_matrix),intent(in) :: smatl
  !!  type(matrices),intent(in) :: ovrlp
  !!  type(matrices),intent(in) :: ham
  !!  type(matrices),intent(in) :: kernel
  !!  real(kind=8),dimension(3,smmd%nat),intent(in) :: rxyz
  !!  character(len=*),intent(in) :: method
  !!  character(len=*),intent(in) :: do_ortho
  !!  character(len=*),intent(in) :: projectormode
  !!  logical,intent(in) :: calculate_multipole_matrices, do_check
  !!  integer,intent(in),optional :: nphi, nphir
  !!  real(kind=8),dimension(:),intent(in),optional :: lphi
  !!  real(kind=8),dimension(3),intent(in),optional :: hgrids
  !!  type(orbitals_data),intent(in),optional :: orbs
  !!  type(comms_linear),intent(in),optional :: collcom, collcom_sr
  !!  type(local_zone_descriptors),intent(in),optional :: lzd
  !!  type(atoms_data),intent(in),optional :: at
  !!  type(DFT_local_fields),intent(inout),optional :: denspot
  !!  type(orthon_data),intent(in),optional :: orthpar
  !!  real(kind=8),dimension(3),intent(in),optional :: shift
  !!  type(matrices),dimension(-lmax:lmax,0:lmax),intent(in),target,optional :: multipole_matrix_in
  !!  type(foe_data),intent(inout),optional :: ice_obj

  !!  ! Local variables
  !!  integer :: methTransformOverlap, iat, ind, ispin, ishift, iorb, jorb, iiorb, l, m, itype, natpx, isatx, nmaxx, kat, n, i, kkat
  !!  integer :: ilr, impl, mm, lcheck, nelpsp, psp_source, j, lwork, ii
  !!  integer, dimension(1) :: power
  !!  logical :: can_use_transposed, all_norms_ok
  !!  real(kind=8),dimension(:),pointer :: phit_c, phit_f
  !!  real(kind=8),dimension(:),allocatable :: phi_ortho, Qmat, kernel_ortho, Qmat_tmp,Slmphi
  !!  real(kind=8),dimension(:),allocatable :: eval, work, newoverlap, newmultipole_matrix_large, newoverlap_large
  !!  real(kind=8),dimension(:,:),allocatable :: Qmat_tilde, kp, locregcenter, overlap_small, tmpmat, tempmat
  !!  real(kind=8),dimension(:,:,:),pointer :: atomic_multipoles
  !!  real(kind=8),dimension(:),pointer :: atomic_monopoles_analytic
  !!  real(kind=8),dimension(:,:,:),allocatable :: test_pot
  !!  real(kind=8),dimension(:,:,:,:),allocatable :: lmp_extracted
  !!  real(kind=8),dimension(:,:),allocatable :: projx
  !!  real(kind=8),dimension(:,:),allocatable :: kernel_extracted, multipole_extracted
  !!  real(kind=8) :: q, tt, rloc, max_error, mean_error
  !!  type(matrices) :: multipole_matrix
  !!  !type(matrices),target :: multipole_matrix_
  !!  type(matrices) :: newovrlp, ovrlp_large, multipole_matrix_large
  !!  type(matrices),dimension(-1:1,0:1) :: lower_multipole_matrices
  !!  type(matrices),dimension(1) :: inv_ovrlp
  !!  logical :: perx, pery, perz
  !!  logical,dimension(:,:),allocatable :: neighborx
  !!  integer,dimension(:),pointer :: nx
  !!  character(len=20),dimension(:),allocatable :: names
  !!  real(kind=8) :: rr1, rr2, rr3
  !!  real(kind=8),dimension(3) :: dipole_check
  !!  real(kind=8),dimension(3,3) :: quadrupole_check
  !!  type(external_potential_descriptors) :: ep_check
  !!  type(matrices),dimension(24) :: rpower_matrix
  !!  type(orbital_basis) :: psi_ob
  !!  real(gp), dimension(3) :: acell, center
  !!  character(len=*),parameter :: no='no', yes='yes'
  !!  type(external_potential_descriptors) :: ep
  !!  !character(len=*),parameter :: projectormode='verynew'!'old'
  !!  !character(len=*),parameter :: do_ortho = no!yes
  !!  integer :: is1, ie1, is2, ie2, is3, ie3, ioffset, icheck
  !!  real(kind=8),dimension(:,:,:),allocatable :: rho_exact, rho_mp, pot_exact, pot_mp
  !!  integer,parameter :: ncheck = 5
  !!  real(kind=8),dimension(ncheck),parameter :: check_threshold = [ 1.d-12 , &
  !!                                                                  1.d-10 , &
  !!                                                                  1.d-8 , &
  !!                                                                  1.d-6 , &
  !!                                                                  1.d-4]
  !!  real(kind=8),dimension(ncheck) :: charge_error, charge_total, potential_error, potential_total
  !!  type(cell) :: mesh


  !!  call f_routine(id='multipole_analysis_driver')

  !!  perx=(smmd%geocode /= 'F')
  !!  pery=(smmd%geocode == 'P')
  !!  perz=(smmd%geocode /= 'F')

  !!  ! Check that the proper optional arguments are present
  !!  if (trim(do_ortho)==yes .and. calculate_multipole_matrices) then
  !!      if (.not.present(nphi) .or. &
  !!          .not.present(orbs) .or. &
  !!          .not.present(lzd) .or. &
  !!          .not.present(collcom) .or. &
  !!          .not.present(lphi) .or. &
  !!          .not.present(orthpar)) then
  !!          call f_err_throw('do_ortho: not all required optional arguments are present')
  !!      end if
  !!  end if
  !!  if (trim(projectormode)=='full') then
  !!      if (.not.present(orbs) .or. &
  !!          .not.present(lzd) .or. &
  !!          .not.present(nphi) .or. &
  !!          .not.present(lphi) .or. &
  !!          .not.present(collcom) .or. &
  !!          .not.present(collcom_sr)) then
  !!          call f_err_throw('projectormode==full: not all required optional arguments are present')
  !!      end if
  !!  end if
  !!  if (lmax>0) then
  !!      if (.not.present(orbs) .or. &
  !!          .not.present(lzd)) then
  !!          call f_err_throw('lmax>0: not all required optional arguments are present')
  !!      end if
  !!      if (.not.calculate_multipole_matrices) then
  !!          call f_err_throw('The multipole matrices must be calculated in-situ for lmax>0')
  !!      end if
  !!  end if

  !!  if (calculate_multipole_matrices) then
  !!      if (.not.present(orbs) .or. &
  !!          .not.present(lzd) .or. &
  !!          .not.present(nphi) .or. &
  !!          .not.present(nphir) .or. &
  !!          .not.present(lphi) .or. &
  !!          .not.present(hgrids) .or. &
  !!          .not.present(collcom)) then
  !!          call f_err_throw('calculate_multipole_matrices .true.: not all required optional arguments are present')
  !!      end if
  !!  else
  !!      if (.not.present(multipole_matrix_in)) then
  !!          call f_err_throw('multipole_matrix_in .false.: not all required optional arguments are present')
  !!      end if
  !!  end if

  !!  if (do_check) then
  !!      if (.not.present(denspot) .or. &
  !!          .not.present(shift) .or. &
  !!          .not.present(lzd) .or. &
  !!          .not.present(at)) then
  !!          call f_err_throw('calculate_multipole_matrices .true.: not all required optional arguments are present')
  !!      end if
  !!  end if

  !!  if (present(lphi) .and. present(nphi)) then
  !!      if (size(lphi)<nphi) then
  !!          call f_err_throw('wrong size of lphi')
  !!      end if
  !!  end if


  !!  if (iproc==0) then
  !!      call yaml_comment('Atomic multipole analysis, new approach',hfill='=')
  !!      call yaml_map('Method',trim(method))
  !!      call yaml_map('Projector mode',trim(projectormode))
  !!      call yaml_map('Orthogonalized support functions',trim(do_ortho))
  !!  end if

  !!  if (calculate_multipole_matrices) then
  !!      call unitary_test_multipoles(iproc, nproc, nphi, nphir, orbs, lzd, smmd, smats, collcom, hgrids)
  !!  end if

  !!  ! Check the method
  !!  if (trim(method)/='projector' .and. trim(method)/='loewdin') then
  !!      call f_err_throw('wrong method',err_name='BIGDFT_RUNTIME_ERROR')
  !!  end if
  !!  if (trim(do_ortho)/='no' .and. trim(do_ortho)/='yes') then
  !!      call f_err_throw('wrong do_ortho',err_name='BIGDFT_RUNTIME_ERROR')
  !!  end if
  !!  select case (trim(projectormode))
  !!  case ('none','simple','full')
  !!      ! everything ok
  !!  case default
  !!      call f_err_throw('wrong projectormode',err_name='BIGDFT_RUNTIME_ERROR')
  !!  end select



  !!  !!multipole_matrix_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='multipole_matrix_large')
  !!  kernel_ortho = sparsematrix_malloc0(smatl,iaction=SPARSE_TASKGROUP,id='kernel_ortho')

  !!  if (do_ortho==yes) then
  !!      methTransformOverlap = 1020
  !!      call matrix_for_orthonormal_basis(iproc, nproc, methTransformOverlap, smats, smatl, &
  !!           ovrlp, kernel, 'plus', kernel_ortho)
  !!   end if

  !!  if (do_ortho==yes .and. calculate_multipole_matrices) then
  !!      ! Orthogonalize the support functions
  !!      can_use_transposed = .false.
  !!      methTransformOverlap = 1020
  !!      phit_c = f_malloc_ptr(collcom%ndimind_c,id='phit_c')
  !!      phit_f = f_malloc_ptr(7*collcom%ndimind_f,id='phit_f')
  !!      phi_ortho = f_malloc(nphi,id='phi_ortho')
  !!      call vcopy(nphi, lphi(1), 1, phi_ortho(1), 1)
  !!      if (iproc==0) then
  !!          call yaml_comment('Orthonormalizing the support functions',hfill='~')
  !!      end if
  !!      call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, &
  !!           1.d-8, nphi, orbs, lzd, &
  !!           smats, smatl, collcom, orthpar, &
  !!           phi_ortho, phit_c, phit_f, &
  !!           can_use_transposed)
  !!      call f_free_ptr(phit_c)
  !!      call f_free_ptr(phit_f)
  !!  end if


  !!  if (trim(method)=='projector') then
  !!      if (smatl%nspin/=1) then
  !!          call f_err_throw('projector not tested for spin polarized calculations, better to stop here')
  !!      end if

  !!      if (projectormode=='full') then
  !!          do i=1,24
  !!              rpower_matrix(i) = matrices_null()
  !!              rpower_matrix(i)%matrix_compr = &
  !!                  sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='rpower_matrix(i)%matrix_compr')
  !!          end do

  !!          !call yaml_warning('NOT SURE WHAT HAPPENS IN THE CASE OF ORTHOGONALIZED SUPPORT FUNCTIONS')
  !!          if (do_ortho==yes) then
  !!              call calculate_rpowerx_matrices(iproc, nproc, nphi, collcom_sr%ndimpsi_c, lzd, &
  !!                   orbs, collcom, phi_ortho, smats, rpower_matrix)
  !!          else if (do_ortho==no) then
  !!              call calculate_rpowerx_matrices(iproc, nproc, nphi, collcom_sr%ndimpsi_c, lzd, &
  !!                   orbs, collcom, lphi, smats, rpower_matrix)
  !!          end if
  !!          ! Calculate the projector using the penalty term
  !!          call projector_for_charge_analysis(smmd, smats, smatm, smatl, &
  !!               ovrlp, ham, kernel, rxyz, calculate_centers=.false., write_output=.false., ortho=do_ortho, mode=projectormode, &
  !!               lzd=lzd, orbs=orbs, natpx=natpx, isatx=isatx, nmaxx=nmaxx, nx=nx, projx=projx, neighborx=neighborx, &
  !!               rpower_matrix=rpower_matrix)
  !!      else
  !!          ! Calculate the projector using the penalty term
  !!          call projector_for_charge_analysis(smmd, smats, smatm, smatl, &
  !!               ovrlp, ham, kernel, rxyz, calculate_centers=.false., write_output=.false., ortho=do_ortho, mode=projectormode, &
  !!               natpx=natpx, isatx=isatx, nmaxx=nmaxx, nx=nx, projx=projx, neighborx=neighborx)
  !!      end if
  !!      if (projectormode=='full') then
  !!          do i=1,24
  !!              call deallocate_matrices(rpower_matrix(i))
  !!          end do
  !!      end if
  !!  else
  !!      ! Just to get the sizes...
  !!      call projector_for_charge_analysis(smmd, smats, smatm, smatl, &
  !!           ovrlp, ham, kernel, rxyz, calculate_centers=.false., write_output=.false., ortho=do_ortho, mode='simple', &
  !!           natpx=natpx, isatx=isatx, nmaxx=nmaxx, nx=nx, projx=projx, neighborx=neighborx, &
  !!           only_sizes=.true.)
  !!      inv_ovrlp(1) = matrices_null()
  !!      inv_ovrlp(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, id='inv_ovrlp%matrix_compr')
  !!      if (do_ortho==yes) then
  !!          ! Calculate the new overlap matrix (which should be more or less the idendity) for the orthogonalized case.
  !!          ! This one is given by S^-1/2*S*S^-1/2
  !!          newoverlap = sparsematrix_malloc(smats, SPARSE_TASKGROUP, id='newoverlap')
  !!          methTransformOverlap = 1020
  !!          newoverlap_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='newoverlap_large')
  !!          ovrlp_large = matrices_null()
  !!          ovrlp_large%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, id='ovrlp_large%matrix_compr')
  !!          call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
  !!               smat_in=ovrlp%matrix_compr, lmat_out=ovrlp_large%matrix_compr)
  !!          call matrix_for_orthonormal_basis(iproc, nproc, methTransformOverlap, smats, smatl, &
  !!               ovrlp, ovrlp_large, 'minus', newoverlap_large)
  !!          call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'large_to_small', &
  !!               lmat_in=newoverlap_large, smat_out=newoverlap)
  !!          call deallocate_matrices(ovrlp_large)
  !!          call f_free(newoverlap_large)
  !!          newovrlp = matrices_null()
  !!          newovrlp%matrix_compr = sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='newovrlp%matrix_compr')
  !!          call f_memcpy(src=newoverlap, dest=newovrlp%matrix_compr)
  !!          power=1
  !!          if (present(ice_obj)) then
  !!              call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
  !!                   1020, 1, power, -1, &
  !!                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
  !!                   ovrlp_mat=newovrlp, inv_ovrlp_mat=inv_ovrlp, &
  !!                   check_accur=.false., ice_obj=ice_obj)
  !!          else
  !!              call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
  !!                   1020, 1, power, -1, &
  !!                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
  !!                   ovrlp_mat=newovrlp, inv_ovrlp_mat=inv_ovrlp, &
  !!                   check_accur=.false.)
  !!          end if
  !!          call deallocate_matrices(newovrlp)
  !!          call f_free(newoverlap)
  !!      else
  !!          power(1)=1
  !!          if (present(ice_obj)) then
  !!              call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
  !!                   1020, 1, power, -1, &
  !!                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
  !!                   ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, &
  !!                   check_accur=.false., ice_obj=ice_obj)
  !!          else
  !!              call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
  !!                   1020, 1, power, -1, &
  !!                   imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
  !!                   ovrlp_mat=ovrlp, inv_ovrlp_mat=inv_ovrlp, &
  !!                   check_accur=.false.)
  !!          end if
  !!      end if
  !!  end if

  !!  !write(*,*) 'before allocate Qmat'

  !!  Qmat = sparsematrix_malloc(smatl,iaction=SPARSE_TASKGROUP,id='Qmat')
  !!  atomic_multipoles = f_malloc0_ptr((/-lmax.to.lmax,0.to.lmax,1.to.smmd%nat/),id='atomic_multipoles')
  !!  !atomic_monopoles_analytic = f_malloc0_ptr(1.to.at%astruct%nat,id='atomic_monopoles_analytic')


  !!  multipole_matrix = matrices_null()
  !!  !multipole_matrix => multipole_matrix_
  !!  multipole_matrix%matrix_compr = sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='multipole_matrix%matrix_compr')



  !!  ! Choose as reference point the midpoint of the simulation cell, in order to avoid
  !!  ! problems with periodic BC (in this way the reference point is always the same and never a periodic image)
  !!  center(1:3) = 0.5d0*smmd%cell_dim(1:3)
  !!  if (calculate_multipole_matrices) then
  !!      locregcenter = f_malloc((/3,lzd%nlr/),id='locregcenter')
  !!      do ilr=1,lzd%nlr
  !!          locregcenter(1:3,ilr) = lzd%llr(ilr)%locregcenter(1:3) !+ (/1.d0,2.d0,3.d0/)
  !!      end do
  !!  end if

  !!  acell(1)=smmd%cell_dim(1)
  !!  acell(2)=smmd%cell_dim(2)
  !!  acell(3)=smmd%cell_dim(3)

  !!  do l=0,lmax
  !!      do m=-l,l

  !!          call f_zero(multipole_matrix%matrix_compr)

  !!          ! Calculate the multipole matrix
  !!          if (calculate_multipole_matrices) then
  !!              if (do_ortho==yes) then
  !!                  call calculate_multipole_matrix(iproc, nproc, l, m, nphi, phi_ortho, phi_ortho, nphir, hgrids, &
  !!                       orbs, collcom, lzd, smmd, smats, locregcenter, 'box', multipole_matrix)
  !!              else if (do_ortho==no) then
  !!                  call calculate_multipole_matrix(iproc, nproc, l, m, nphi, lphi, lphi, nphir, hgrids, &
  !!                       orbs, collcom, lzd, smmd, smats, locregcenter, 'box', multipole_matrix) 
  !!              end if
  !!          else
  !!              !multipole_matrix => multipole_matrix_in(m,l)
  !!              if (do_ortho==yes)  then
  !!                  methTransformOverlap = 1020
  !!                  newmultipole_matrix_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='newmultipole_matrix_large')
  !!                  !multipole_matrix_large = sparsematrix_malloc(smatl, SPARSE_TASKGROUP, id='multipole_matrix_large')
  !!                  multipole_matrix_large = matrices_null()
  !!                  multipole_matrix_large%matrix_compr = sparsematrix_malloc_ptr(smatl, SPARSE_TASKGROUP, &
  !!                      id='multipole_matrix_large%matrix_compr')
  !!                  call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
  !!                       smat_in=multipole_matrix_in(m,l)%matrix_compr, lmat_out=multipole_matrix_large%matrix_compr)
  !!                  call matrix_for_orthonormal_basis(iproc, nproc, methTransformOverlap, smats, smatl, &
  !!                       ovrlp, multipole_matrix_large, 'minus', newmultipole_matrix_large)
  !!                  call transform_sparse_matrix(iproc, smats, smatl, SPARSE_TASKGROUP, 'small_to_large', &
  !!                       lmat_in=newmultipole_matrix_large, smat_out=multipole_matrix%matrix_compr)
  !!                  !call f_free(multipole_matrix_large)
  !!                  call deallocate_matrices(multipole_matrix_large)
  !!                  call f_free(newmultipole_matrix_large)
  !!              else
  !!                  call f_memcpy(src=multipole_matrix_in(m,l)%matrix_compr, dest=multipole_matrix%matrix_compr)
  !!              end if
  !!          end if

  !!          if (l<=1) then
  !!              lower_multipole_matrices(m,l) = matrices_null()
  !!              lower_multipole_matrices(m,l)%matrix_compr = &
  !!                  sparsematrix_malloc_ptr(smats, SPARSE_TASKGROUP, id='lower_multipole_matrix%matrix_compr')
  !!              call f_memcpy(src=multipole_matrix%matrix_compr,dest=lower_multipole_matrices(m,l)%matrix_compr)
  !!          end if


  !!          !write(*,*) 'before f_zero'
  !!          call f_zero(Qmat)

  !!          if (trim(method)=='projector' .or. trim(method)=='loewdin') then
  !!              do kat=1,natpx
  !!                  kkat = kat + isatx
  !!                  n = nx(kat)
  !!                  qmat_tilde = f_malloc((/n,n/),id='qmat_tilde')
  !!                  kp = f_malloc((/n,n/),id='kp')
  !!                   kernel_extracted = f_malloc((/n,n/),id='kernel_extracted')
  !!                   multipole_extracted = f_malloc((/n,n/),id='multipole_extracted')
  !!                   !write(*,*) 'before extract'
  !!                   call extract_matrix(smats, multipole_matrix%matrix_compr, &
  !!                       neighborx(1,kat), n, nmaxx, multipole_extracted)
  !!                   ! The minus sign is required since the phi*S_lm*phi represent the electronic charge which is a negative quantity
  !!                   call dscal(n**2, -1.d0, multipole_extracted(1,1), 1)
  !!                   if (do_ortho==no) then
  !!                       call extract_matrix(smatl, kernel%matrix_compr, neighborx(1,kat), n, nmaxx, kernel_extracted)
  !!                   else
  !!                       call extract_matrix(smatl, kernel_ortho, neighborx(1,kat), n, nmaxx, kernel_extracted)
  !!                   end if
  !!                   if (l>0) then
  !!                       call correct_multipole_origin(smmd%nat, l, m, n, lzd%nlr, natpx, nmaxx, kat, kkat, &
  !!                            smats, orbs, rxyz, neighborx, perx, pery, perz, acell, &
  !!                            lower_multipole_matrices, locregcenter, multipole_extracted)
  !!                   end if
  !!                   !!do i=1,n
  !!                   !!    do j=1,n
  !!                   !!        write(*,*) 'i, j, multipole_extracted(j,i)', i, j, multipole_extracted(j,i)
  !!                   !!    end do
  !!                   !!end do
  !!                   call gemm('n', 'n', n, n, n, 1.d0, kernel_extracted(1,1), n, &
  !!                        multipole_extracted(1,1), n, 0.d0, qmat_tilde(1,1), n)
  !!                   call f_free(kernel_extracted)
  !!                   call f_free(multipole_extracted)
  !!                   if (trim(method)=='loewdin') then
  !!                           call extract_matrix(smatl, inv_ovrlp(1)%matrix_compr, neighborx(1,kat), n, nmaxx, projx(1,kat))
  !!                           iiorb = 0
  !!                           do iorb=1,smats%nfvctr
  !!                               if (neighborx(iorb,kat)) then
  !!                                   iiorb = iiorb + 1
  !!                                   if (smmd%on_which_atom(iorb)/=kkat) then
  !!                                       do jorb=1,n
  !!                                           projx((iiorb-1)*n+jorb,kat) = 0.d0
  !!                                       end do
  !!                                   end if
  !!                               end if
  !!                           end do
  !!                   end if
  !!                   !!do i=1,n**2
  !!                   !!    write(*,*) 'i, j, projx(i,kat)', i, j, projx(i,kat)
  !!                   !!end do
  !!                   !write(*,*) 'sum(projx)', sum(projx)
  !!                  call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, projx(1,kat), n, 0.d0, kp(1,1), n)
  !!                  !!!! LUIGI'S NEW IDEA #############################################################
  !!                  !!!do i=1,n
  !!                  !!!    kp(i,i) = 0.d0
  !!                  !!!    do j=1,n
  !!                  !!!        if (orbs%onwhichatom(j)==kkat .and. orbs%onwhichatom(i)==kkat) then
  !!                  !!!            kp(i,i) = kp(i,i) - 0.25*qmat_tilde(j,i)**2
  !!                  !!!        end if
  !!                  !!!    end do
  !!                  !!!end do
  !!                  !!!! END OF LUIGI'S NEW IDEA ######################################################
  !!                  if (do_ortho==no) then
  !!                      overlap_small = f_malloc((/n,n/),id='overlap_small')
  !!                      call extract_matrix(smats, ovrlp%matrix_compr, neighborx(1,kat), n, nmaxx, overlap_small)
  !!                      call f_memcpy(src=kp,dest=qmat_tilde)
  !!                      call gemm('n', 'n', n, n, n, 1.d0, qmat_tilde(1,1), n, overlap_small(1,1), n, 0.d0, kp(1,1), n)
  !!                      call f_free(overlap_small)
  !!                  end if
  !!                  tt = 0.d0
  !!                  do i=1,n
  !!                      tt = tt + kp(i,i)
  !!                      !write(*,*) 'kat, i, owa(i), qmat_tilde(i,i), kp(i,i)', &
  !!                      !    kat, i, smmd%on_which_atom(i), qmat_tilde(i,i), kp(i,i)
  !!                  end do
  !!                  atomic_multipoles(m,l,kkat) = tt
  !!                  !write(*,*) 'm, l, kkat, atomic_multipoles(m,l,kkat)', m, l, kkat, atomic_multipoles(m,l,kkat)
  !!                  call f_free(qmat_tilde)
  !!                  call f_free(kp)
  !!              end do
  !!          end if

  !!      end do
  !!  end do


  !!  if (calculate_multipole_matrices) then
  !!      call f_free(locregcenter)
  !!  end if

  !!  call mpiallred(atomic_multipoles, mpi_sum, comm=bigdft_mpi%mpi_comm)


  !!  ! The monopole term should be the net charge, i.e. add the positive atomic charges
  !!  do iat=1,smmd%nat
  !!      itype = smmd%iatype(iat)
  !!      q = real(smmd%nelpsp(itype),kind=8)
  !!      atomic_multipoles(0,0,iat) = atomic_multipoles(0,0,iat) + q
  !!  end do


  !!  names = f_malloc_str(len(names),smmd%nat,id='names')
  !!  do iat=1,smmd%nat
  !!      itype = smmd%iatype(iat)
  !!      names(iat) = smmd%atomnames(itype)
  !!  end do


  !!  ep = external_potential_descriptors_null()
  !!  ep%nmpl = smmd%nat
  !!  allocate(ep%mpl(ep%nmpl))
  !!  do impl=1,ep%nmpl
  !!      ep%mpl(impl) = multipole_set_null()
  !!      allocate(ep%mpl(impl)%qlm(0:lmax))
  !!      ep%mpl(impl)%rxyz = smmd%rxyz(1:3,impl)
  !!      ep%mpl(impl)%sym = trim(names(impl))
  !!      if (present(at)) then
  !!          call get_psp_info(ep%mpl(impl)%sym, ixc, smmd, nelpsp, psp_source, rloc, at%psppar)
  !!      else
  !!          call get_psp_info(ep%mpl(impl)%sym, ixc, smmd, nelpsp, psp_source, rloc)
  !!      end if
  !!      if (psp_source/=0 .and. iproc==0) then
  !!          call yaml_warning('Taking internal PSP information for multipole '//trim(yaml_toa(impl)))
  !!      end if
  !!      ep%mpl(impl)%nzion = nelpsp
  !!      ep%mpl(impl)%sigma(0:lmax) = rloc
  !!      do l=0,lmax
  !!          ep%mpl(impl)%qlm(l) = multipole_null()
  !!          !if (l>=3) cycle
  !!          ep%mpl(impl)%qlm(l)%q = f_malloc0_ptr(2*l+1,id='q')
  !!          mm = 0
  !!          !if (l>0) cycle
  !!          do m=-l,l
  !!              mm = mm + 1
  !!              ep%mpl(impl)%qlm(l)%q(mm) = atomic_multipoles(m,l,impl)
  !!          end do
  !!      end do
  !!  end do

  !!  if (iproc==0) then
  !!      call yaml_comment('Final result of the multipole analysis',hfill='~')
  !!      call write_multipoles_new(ep, lmax, smmd%units)
  !!  end if


  !!  if (do_check) then
  !!      ! Calculate the total dipole moment resulting from the previously calculated multipoles.
  !!      ! This is done by calling the following routine (which actually calculates the potential, but also
  !!      ! has the option to calculate the dipole on the fly).
  !!      test_pot = f_malloc((/size(denspot%V_ext,1),size(denspot%V_ext,2),size(denspot%V_ext,3)/),id='test_pot')
  !!      if (iproc==0) call yaml_sequence_open('Checking the total multipoles based on the atomic multipoles')
  !!      is1 = 1
  !!      ie1 = denspot%dpbox%mesh%ndims(1)
  !!      is2 = 1
  !!      ie2 = denspot%dpbox%mesh%ndims(2)
  !!      is3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1
  !!      ie3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
  !!            denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2)
  !!      rho_exact = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='rho_exact')
  !!      rho_mp = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='rho_mp')
  !!      pot_exact = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='pot_exact')
  !!      pot_mp = f_malloc((/is1.to.ie1,is2.to.ie2,is3.to.ie3/),id='pot_mp')
  !!      do lcheck=0,lmax
  !!          ep_check = external_potential_descriptors_null()
  !!          ep_check%nmpl = ep%nmpl
  !!          allocate(ep_check%mpl(ep_check%nmpl))
  !!          do impl=1,ep_check%nmpl
  !!              ep_check%mpl(impl) = multipole_set_null()
  !!              allocate(ep_check%mpl(impl)%qlm(0:lmax))
  !!              ep_check%mpl(impl)%rxyz = ep%mpl(impl)%rxyz
  !!              ep_check%mpl(impl)%sym = ep%mpl(impl)%sym
  !!              ep_check%mpl(impl)%nzion = ep%mpl(impl)%nzion
  !!              do l=0,lmax
  !!                  ep_check%mpl(impl)%sigma(l) = ep%mpl(impl)%sigma(l)
  !!                  ep_check%mpl(impl)%qlm(l) = multipole_null()
  !!                  if (l>lcheck) cycle
  !!                  ep_check%mpl(impl)%qlm(l)%q = f_malloc0_ptr(2*l+1,id='q')
  !!                  mm = 0
  !!                  do m=-l,l
  !!                      mm = mm + 1
  !!                      ep_check%mpl(impl)%qlm(l)%q(mm) = ep%mpl(impl)%qlm(l)%q(mm)
  !!                  end do
  !!              end do
  !!          end do
  !!          call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
  !!               denspot%V_ext(1,1,1,1), 1, test_pot(1,1,1), 1)
  !!          call potential_from_charge_multipoles(iproc, nproc, at, denspot, ep_check, 1, &
  !!               denspot%dpbox%mesh%ndims(1), 1, denspot%dpbox%mesh%ndims(2), &
  !!               denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
  !!               denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+&
  !!               denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
  !!               denspot%dpbox%mesh%hgrids(1),denspot%dpbox%mesh%hgrids(2),denspot%dpbox%mesh%hgrids(3), &
  !!               shift, verbosity=0, ixc=ixc, lzd=lzd, pot=test_pot, &
  !!               rxyz=rxyz, dipole_total=dipole_check, quadrupole_total=quadrupole_check, &
  !!               all_norms_ok=all_norms_ok, &
  !!               rho_mp=rho_mp, pot_mp=pot_mp)
  !!          if (.not. all_norms_ok) then
  !!              call f_err_throw('When checking the previously calculated multipoles, all norms should be ok')
  !!          end if
  !!          dipole_check=dipole_check/Debye_AU  ! au2debye              

  !!          !# NEW: compare the density and potential ##########################
  !!          if (smatl%nspin/=1) then
  !!              call f_err_throw('Multipole analysis check not yet ready for nspin>1')
  !!          end if
  !!          ! Get the exact charge density
  !!          ioffset = denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*&
  !!                    denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,4)
  !!          !write(*,*) 'MP: ioffset', ioffset
  !!          call f_memcpy(n=(ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), &
  !!               src=denspot%rhov(ioffset+1), dest=rho_exact(is1,is2,is3))
  !!          call f_memcpy(src=rho_exact, dest=pot_exact)
!!!!             call H_potential('D',denspot%pkernel,pot_exact,denspot%V_ext,tt,0.0_dp,.true.,&
!!!!                  quiet='yes')
  !!          call Electrostatic_Solver(denspot%pkernel,pot_exact,pot_ion=denspot%V_ext)
  !!          !mesh=cell_new(smmd%geocode,denspot%pkernel%ndims,denspot%pkernel%hgrids)
  !!          call compare_charge_and_potential(denspot%dpbox%bitp,&!iproc, is1, ie1, is2, ie2, is3, ie3, &
  !!               smmd%nat, &
  !!               rho_exact, rho_mp, pot_exact, pot_mp, denspot%pkernel, rxyz, &
  !!               ncheck, check_threshold, charge_error, charge_total, potential_error, potential_total)
  !!          !# NEW: compare the density and potential ##########################
  !!          if (iproc==0) then
  !!              call yaml_sequence(advance='no')
  !!              call yaml_mapping_open('Up to multipole l='//trim(yaml_toa(lcheck)))
  !!              call yaml_mapping_open('Electric Dipole Moment (Debye)')
  !!              call yaml_map('P vector',dipole_check(1:3),fmt='(1es13.4)')
  !!              call yaml_map('norm(P)',sqrt(sum(dipole_check**2)),fmt='(1es14.6)')
  !!              call yaml_mapping_close()
  !!              call yaml_mapping_open('Quadrupole Moment (AU)')
  !!              call yaml_map('Q matrix',quadrupole_check,fmt='(1es13.4)')
  !!              call yaml_map('trace',quadrupole_check(1,1)+quadrupole_check(2,2)+quadrupole_check(3,3),fmt='(es12.2)')
  !!              call yaml_mapping_close()
  !!              call yaml_sequence_open('Average relative error of resulting potential in the Exterior region')
  !!              !call yaml_sequence_open('density threshold for check')
  !!              do icheck=1,ncheck
  !!                  !call yaml_mapping_open('density threshold for check',check_threshold(icheck))
  !!                  call yaml_sequence(advance='no')
  !!                  call yaml_mapping_open(flow=.true.)
  !!                  call yaml_map('Thr',check_threshold(icheck),fmt='(es9.1)')
  !!                  call yaml_map('Ext. Vol. %',charge_total(icheck)/&
  !!                       (denspot%dpbox%mesh%volume_element*product(real(denspot%dpbox%mesh%ndims,gp))),&
  !!                       fmt='(2pf5.1)')
  !!                  call yaml_map('int(V)',potential_total(icheck),fmt='(es10.3)')
  !!                  call yaml_map('Err %',potential_error(icheck)/potential_total(icheck),fmt='(2pf5.1)')
  !!                  call yaml_map('int(rho)',charge_error(icheck),fmt='(es10.3)')
!!!!                     !call yaml_mapping_open('density threshold for check is'//&
!!!!                     !     &trim(yaml_toa(check_threshold(icheck),fmt='(es9.2)')))
!!!!                     call yaml_mapping_open('rho',flow=.true.)
!!!!                     call yaml_map('int(q-q_exact))',charge_error(icheck),fmt='(es10.3)')
!!!!                     call yaml_map('int(q_exact)',charge_total(icheck),fmt='(es10.3)')
!!!!                     call yaml_map('ratio',charge_error(icheck)/charge_total(icheck),fmt='(es10.3)')
!!!!                     call yaml_mapping_close()
!!!!                     call yaml_mapping_open('pot',flow=.true.)
!!!!                     call yaml_map('int(V-V_exact))',potential_error(icheck),fmt='(es10.3)')
!!!!                     call yaml_map('int(V_exact)',potential_total(icheck),fmt='(es10.3)')
!!!!                     call yaml_map('ratio',potential_error(icheck)/potential_total(icheck),fmt='(es10.3)')
!!!!                     call yaml_mapping_close()
  !!                  call yaml_mapping_close()
  !!              end do
  !!              call yaml_sequence_close()
  !!              !call yaml_mapping_close()
  !!              call yaml_mapping_close()
  !!           end if
  !!           call deallocate_external_potential_descriptors(ep_check)
  !!        end do
  !!      if (iproc==0) call yaml_sequence_close()
  !!      call f_free(test_pot)
  !!      call f_free(rho_exact)
  !!      call f_free(rho_mp)
  !!      call f_free(pot_exact)
  !!      call f_free(pot_mp)
  !!  end if

  !!  do l=0,min(1,lmax)
  !!      do m=-l,l
  !!          call deallocate_matrices(lower_multipole_matrices(m,l))
  !!      end do
  !!  end do

  !!  if (trim(method)=='loewdin') then
  !!      call deallocate_matrices(inv_ovrlp(1))
  !!  end if

  !!  call f_free_str(len(names),names)
  !!  call deallocate_matrices(multipole_matrix)
  !!  call f_free(kernel_ortho)
  !!  call f_free(Qmat)
  !!  if (do_ortho==yes .and. calculate_multipole_matrices) then
  !!      call f_free(phi_ortho)
  !!  end if
  !!  call f_free(projx)
  !!  call f_free_ptr(nx)
  !!  call f_free(neighborx)
  !!  call f_free_ptr(atomic_multipoles)
  !!  !!call f_free(multipole_matrix_large)
  !!  call deallocate_external_potential_descriptors(ep)

  !!  if (iproc==0) then
  !!      call yaml_comment('Atomic multipole analysis done',hfill='=')
  !!  end if

  !!  call f_release_routine()

  !!end subroutine multipole_analysis_driver
