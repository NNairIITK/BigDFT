
    subroutine projector_for_charge_analysis(at, smats, smatm, smatl, &
               ovrlp_, ham_, kernel_, rxyz, calculate_centers, &
               lzd, nphirdim, psi, orbs)
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
      type(matrices),intent(inout) :: ovrlp_ !< should be intent(in)...
      type(matrices),intent(in) :: ham_, kernel_
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      logical,intent(in) :: calculate_centers
      type(local_zone_descriptors),intent(in),optional :: lzd
      integer,intent(in),optional :: nphirdim
      real(kind=8),dimension(:),intent(in),optional :: psi
      type(orbitals_data),intent(in),optional :: orbs

      ! Local variables
      integer :: kat, iat, jat, i, j, ii, jj, icheck, n, indm, inds, ntot, ist, ind, iq, itype, ieval, ij, nmax, indl, lwork
      integer :: k, l, iatold, isat, natp, kkat, istot, ntotp, i1, i2, i3, is1, ie1, is2, ie2, is3, ie3, j1, j2, j3, ikT, info
      integer :: ialpha
      real(kind=8) :: r2, cutoff2, rr2, tt, ef, q, occ, max_error, mean_error, rr2i, rr2j, ttxi, ttyi, ttzi, ttxj, ttyj, ttzj
      real(kind=8) :: tti, ttj, charge_net, charge_total
      real(kind=8) :: xi, xj, yi, yj, zi, zj, ttx, tty, ttz, xx, yy, zz, x, y, z
      real(kind=8),dimension(:),allocatable :: work, alpha_arr
      real(kind=8),dimension(:,:),allocatable :: com, alpha_arr_bounds
      real(kind=8),dimension(:,:),allocatable :: ham, ovrlp, proj, ovrlp_tmp, ovrlp_minusonehalf, kp, ktilde
      real(kind=8),dimension(:,:,:),allocatable :: coeff_all, ovrlp_onehalf_all
      integer,dimension(:,:,:,:),allocatable :: ilup
      real(kind=8),dimension(:),allocatable :: eval, eval_all, ovrlp_large, tmpmat1, tmpmat2, kerneltilde, charge_per_atom
      real(kind=8),dimension(:,:,:),allocatable :: tmpmat2d
      integer,dimension(:),allocatable :: id_all, n_all, itmparr
      real(kind=8),dimension(3) :: rr
      logical,dimension(:,:),allocatable :: neighbor
      type(matrices),dimension(1) :: ovrlp_onehalf_
      logical :: perx, pery, perz, final, bound_low_ok, bound_up_ok
      !real(kind=8),parameter :: kT = 5.d-2
      real(kind=8) :: kT, eval_target
      !real(kind=8),parameter :: alpha = 5.d-1
      real(kind=8) :: alpha, alpha_up, alpha_low, convergence_criterion
      real(kind=8) :: ef_up, ef_low


      call f_routine(id='projector_for_charge_analysis')

      alpha_arr = f_malloc(at%astruct%ntypes,id='alpha_arr')
      alpha_arr_bounds = f_malloc((/at%astruct%ntypes,2/),id='alpha_arr_bounds')

      alpha_arr(1) = 17.777777778d0
      alpha_arr(2) = 5.25099769d0
      !alpha_arr(:) = 17.777777778d0
      !alpha_arr(:) = 5.25099769d0
      !alpha_arr(1) = 63.2d0*10.d-1
      !alpha_arr(2) = 5.51459534d0*10.d-1
      !alpha_arr(:) = 63.2d0*10.d-1
      !alpha_arr(:) = 5.51459534d0*10.d-1
      alpha_arr_bounds(:,1) = 0.1d0*alpha_arr(:)
      alpha_arr_bounds(:,2) = 10.d0*alpha_arr(:)

      kT = 1.d-2

      ! Convergence criterion: one million-th of the total charge
      tt = 0.d0
      do iat=1,at%astruct%nat
          tt = tt + real(at%nelpsp(at%astruct%iatype(iat)),kind=8)
      end do
      convergence_criterion = 1.d-6*abs(tt)

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


      ! Calculate S^1/2
      ovrlp_onehalf_(1) = matrices_null()
      ovrlp_onehalf_(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_onehalf_(1)%matrix_compr')
      call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, 1020, 1, (/2/), -1, &
            imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
            ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_(1), &
            check_accur=.true., max_error=max_error, mean_error=mean_error)

      ! Calculate S^1/2 * K * S^1/2 = Ktilde
      tmpmat1 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat1')
      !tmpmat2 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat2')
      kerneltilde = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='kerneltilde')
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           kernel_%matrix_compr, ovrlp_onehalf_(1)%matrix_compr, tmpmat1)
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           ovrlp_onehalf_(1)%matrix_compr, tmpmat1, kerneltilde)


      ! Determine the periodicity...
      perx=(smats%geocode /= 'F')
      pery=(smats%geocode == 'P')
      perz=(smats%geocode /= 'F')
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
      do kat=1,natp
          iatold = 0
          kkat = kat + isat
          n = 0
          do i=1,smats%nfvctr
               iat = smats%on_which_atom(i)
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
      id_all = f_malloc0(ntot,id='id_all')
      coeff_all = f_malloc((/nmax,nmax,natp/),id='coeff_all')
      ovrlp_onehalf_all = f_malloc((/nmax,nmax,natp/),id='ovrlp_onehalf_all')
      ovrlp_minusonehalf = f_malloc((/nmax,nmax/),id='ovrlp_minusonehalf')
      ilup = f_malloc((/2,nmax,nmax,natp/),id='ilup')
      n_all = f_malloc(natp,id='n_all')


      ! Centers of the support functions
      com = f_malloc0((/3,smats%nfvctr/),id='com')
      if (calculate_centers) then
          if (orbs%norb>0) then
              !call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, tmb%collcom_sr%ndimpsi_c, &
              !     orbs%norb, orbs%norbp, orbs%isorb, orbs%in_which_locreg, lzd, com(1:,orbs%isorb+1:))
              call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, nphirdim, &
                   orbs%norb, orbs%norbp, orbs%isorb, orbs%inwhichlocreg, lzd, com(1:,orbs%isorb+1:))
              if (bigdft_mpi%nproc>1) then
                  call mpiallred(com, mpi_sum, comm=bigdft_mpi%mpi_comm)
              end if
          end if
      else
          do i=1,smats%nfvctr
              iat = smats%on_which_atom(i)
              com(1:3,i) = rxyz(1:3,iat)
          end do
      end if


      charge_per_atom = f_malloc0(at%astruct%nat,id='charge_per_atom')


      ! Calculate how many states should be included
      q = 0.d0
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          q = q + ceiling(0.5d0*real(at%nelpsp(itype),kind=8))
      end do
      iq = nint(q)
      if (bigdft_mpi%iproc==0) then
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

      alpha_low = 0.d-0
      alpha_up = 0.d-0

      alpha_loop: do ialpha=1,100

          if (bigdft_mpi%iproc==0) then
              call yaml_sequence(advance='no')
          end if
          
          if (bigdft_mpi%iproc==0) then
              write(*,*) 'bound_low_ok, bound_up_ok', bound_low_ok, bound_up_ok
              write(*,*) 'alpha_arr_bounds(:,1)',alpha_arr_bounds(:,1)
              write(*,*) 'alpha_arr_bounds(:,2)',alpha_arr_bounds(:,2)
          end if

          if (.not.bound_low_ok) then
              alpha = alpha_low
              alpha_arr(:) = alpha_arr_bounds(:,1)
          else if (.not.bound_up_ok) then
              alpha = alpha_up
              alpha_arr(:) = alpha_arr_bounds(:,2)
          else
              alpha = 0.5d0*(alpha_low+alpha_up)
              alpha_arr(:) = 0.5d0*(alpha_arr_bounds(:,1)+alpha_arr_bounds(:,2))
          end if

          charge_net = 0.d0
          call f_zero(eval_all)
          call f_zero(id_all)

          if (bigdft_mpi%iproc==0) write(*,*) 'ialpha, alpha_arr', ialpha, alpha_arr

          ist = 0
          do kat=1,natp
              kkat = kat + isat

              alpha = alpha_arr(at%astruct%iatype(kkat))
              !write(*,*) 'kkat, name, alpha', kkat, at%astruct%atomnames(at%astruct%iatype(kkat)), alpha
              !if (kkat==2) then
              !    !O
              !    alpha = 5.51459534d0*10.d-1
              !    !alpha = 63.2d0*5.d-1
              !else
              !    !H
              !    alpha = 63.2d0*10.d-1
              !end if
    
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
              call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1:,kat), n, nmax, ovrlp, ilup)
              call extract_matrix(smatm, ham_%matrix_compr, neighbor(1:,kat), n, nmax, ham)
              !!icheck = 0
              !!ii = 0
              !!do i=1,smats%nfvctr
              !!    if (neighbor(i,kat)) then
              !!        jj = 0
              !!        do j=1,smats%nfvctr
              !!            if (neighbor(j,kat)) then
              !!                icheck = icheck + 1
              !!                jj = jj + 1
              !!                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
              !!                inds =  matrixindex_in_compressed(smats, i, j)
              !!                if (inds>0) then
              !!                    ovrlp(jj,ii) = ovrlp_%matrix_compr(inds)
              !!                else
              !!                    ovrlp(jj,ii) = 0.d0
              !!                end if
              !!                indm =  matrixindex_in_compressed(smatm, i, j)
              !!                if (indm>0) then
              !!                    ham(jj,ii) = ham_%matrix_compr(indm)
              !!                else
              !!                    ham(jj,ii) = 0.d0
              !!                end if
              !!            end if
              !!        end do
              !!    end if
              !!end do
              !!if (icheck>n**2) then
              !!    call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
              !!        &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
              !!end if


              ! Calculate ovrlp^1/2 and ovrlp^-1/2. The last argument is wrong, clean this.
              ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
              call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
              call overlap_plus_minus_one_half_exact(1, n, -1, .true., ovrlp_tmp, smats)
              do i=1,n
                  call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
              end do
              call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
              call overlap_plus_minus_one_half_exact(1, n, -1, .false., ovrlp_tmp, smats)
              do i=1,n
                  call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_minusonehalf(1,i), 1)
              end do
              call f_free(ovrlp_tmp)
    
              ! Calculate S^-1/2 * H * S^-1/2
              tmpmat2d = f_malloc((/n,n,1/),id='tmppmat2d')
              call gemm('n', 'n', n, n, n, 1.d0, ham(1,1), n, ovrlp_minusonehalf(1,1), nmax, 0.d0, tmpmat2d(1,1,1), n)
              call gemm('n', 'n', n, n, n, 1.d0, ovrlp_minusonehalf(1,1), nmax, tmpmat2d(1,1,1), n, 0.d0, ham(1,1), n)
              call f_free(tmpmat2d)

              ! Add the penalty term
              call add_penalty_term(smats%geocode, smats%nfvctr, neighbor(1:,kat), rxyz(1:,kkat), &
                   at%astruct%cell_dim, com, alpha, n, ovrlp, ham)
              !!icheck = 0
              !!ii = 0
              !!do i=1,smats%nfvctr
              !!    if (neighbor(i,kat)) then
              !!        jj = 0
              !!        do j=1,smats%nfvctr
              !!            if (neighbor(j,kat)) then
              !!                icheck = icheck + 1
              !!                jj = jj + 1
              !!                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
              !!                if (i==j) then
              !!                    rr2 = huge(rr2)
              !!                    do i3=is3,ie3
              !!                        z = rxyz(3,kkat) + i3*at%astruct%cell_dim(3)
              !!                        ttz = (com(3,i)-z)**2
              !!                        do i2=is2,ie2
              !!                            y = rxyz(2,kkat) + i2*at%astruct%cell_dim(2)
              !!                            tty = (com(2,i)-y)**2
              !!                            do i1=is1,ie1
              !!                                x = rxyz(1,kkat) + i1*at%astruct%cell_dim(1)
              !!                                ttx = (com(1,i)-x)**2
              !!                                tt = ttx + tty + ttz
              !!                                if (tt<rr2) then
              !!                                    rr2 = tt
              !!                                end if
              !!                            end do
              !!                        end do
              !!                    end do
              !!                    ham(jj,ii) = ham(jj,ii) + alpha*rr2**3*ovrlp(jj,ii)
              !!                end if
              !!                ilup(1,jj,ii,kat) = j
              !!                ilup(2,jj,ii,kat) = i
              !!            end if
              !!        end do
              !!    end if
              !!end do
              !!if (icheck>n**2) then
              !!    call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
              !!        &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
              !!end if
    
    
              !!call diagonalizeHamiltonian2(bigdft_mpi%iproc, n, ham, ovrlp, eval)
              lwork = 10*n
              work = f_malloc(lwork,id='work')
              call syev('v', 'l', n, ham(1,1), n, eval(1), work(1), lwork, info)
              call f_free(work)
              do i=1,n
                  ii = ist + i
                  eval_all(istot+ii) = eval(i)
                  !write(*,*) 'kkat, alpha, i, eval', kkat, alpha, i, eval(i)
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
          call order_eigenvalues(ntot, eval_all, id_all)
          !do i=1,ntot
          !    ! add i-1 since we are only searching in the subarray
          !    ind = minloc(eval_all(i:ntot),1) + (i-1)
          !    tt = eval_all(i)
          !    eval_all(i) = eval_all(ind)
          !    eval_all(ind) = tt
          !    ii = id_all(i)
          !    id_all(i) = id_all(ind)
          !    id_all(ind) = ii
          !end do

          ii = 0
          do iat=1,at%astruct%nat
              itype = at%astruct%iatype(iat)
              ii = ii + at%nelpsp(itype)
          end do
          ii = ii/2 !closed shell. quick and dirty
          !eval_target = -2.63d-1+0.05d0
          !eval_target = -2.627921676582d-1+.1d0 !standard
          eval_target = -2.637375544845d-1+.1d0 !large
          if (bigdft_mpi%iproc==0) write(*,*) 'eval_target, ii, eval_all(ii)', eval_target, ii, eval_all(ii)
          if (abs(eval_all(ii)-eval_target)<1.d-4) then
              exit alpha_loop
          end if

          !if (eval_all(ii)>eval_target) then
          !    alpha_arr(:) = 0.8d0*alpha_arr(:)
          !else
          !    alpha_arr(:) = 1.2d0*alpha_arr(:)
          !end if


          ! If we are still searching the boundaries for the bisection...
          if (.not.bound_low_ok) then
              if (eval_all(ii)<eval_target) then
                  ! this is an lower bound
                  alpha_arr_bounds(:,1) = alpha_arr(:)
                  bound_low_ok = .true.
              else
                  alpha_arr_bounds(:,1) = 0.5d0*alpha_arr_bounds(:,1)
              end if
              cycle alpha_loop
          else if (.not.bound_up_ok) then
              if (eval_all(ii)>eval_target) then
                  ! this is an upper bound
                  alpha_arr_bounds(:,2) = alpha_arr(:)
                  bound_up_ok = .true.
              else
                  alpha_arr_bounds(:,2) = 2.0d0*alpha_arr_bounds(:,2)
              end if
              cycle alpha_loop
          end if

          if (eval_all(ii)<eval_target) then
              ! new lower bound
              if (bigdft_mpi%iproc==0) write(*,*) 'new lower bound'
              alpha_arr_bounds(:,1) = alpha_arr(:)
          else if (eval_all(ii)>eval_target) then
              ! new upper bound
              if (bigdft_mpi%iproc==0) write(*,*) 'new upper bound'
              alpha_arr_bounds(:,2) = alpha_arr(:)
          end if
        
        
      end do alpha_loop
    
    
          !!! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
          !!ef = eval_all(1)
          !!do
          !!    ef = ef + 1.d-3
          !!    occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
          !!    if (abs(occ-1.d0)<1.d-8) exit
          !!end do
          !!if (bigdft_mpi%iproc==0) then
          !!    call yaml_map('Pseudo Fermi level for occupations',ef)
          !!end if
        
          !!final = .true.
          !!ikT = 0
          !!kT_loop: do
    
              !ikT = ikT + 1
    
              call f_zero(charge_per_atom)
    
              ! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
              ef = eval_all(1)
              do
                  ef = ef + 1.d-3
                  occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
                  if (abs(occ-1.d0)<1.d-8) exit
              end do
              ef = -6.0e-4

              bound_low_ok = .false.
              bound_up_ok = .false.

              ef_low = -10.d0
              ef_up = 10.d0



              ef_loop: do

              call f_zero(charge_per_atom)

                  if (.not.bound_low_ok) then
                      ef = ef_low
                  else if (.not.bound_up_ok) then
                      ef = ef_up
                  else
                      ef = 0.5d0*(ef_low+ef_up)
                  end if
        
                  ! Calculate the projector. First for each single atom, then insert it into the big one.
                  charge_total = 0.d0
                  do kat=1,natp
                      kkat = kat + isat
                      n = n_all(kat)
                      proj = f_malloc0((/n,n/),id='proj')
                      call calculate_projector(n, ntot, nmax, kkat, id_all, eval_all, &
                           coeff_all(1:,1:,kat), ef, kT, proj)
                      !ij = 0
                      !do ieval=1,ntot
                      !    if (id_all(ieval)/=kkat) cycle
                      !    ij = ij + 1
                      !    occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
                      !    do i=1,n
                      !        do j=1,n
                      !            proj(j,i) = proj(j,i) + occ*coeff_all(j,ij,kat)*coeff_all(i,ij,kat)
                      !        end do
                      !   end do
                      !end do
                      !tt = 0.d0
                      !do i=1,n
                      !    tt = tt + proj(i,i)
                      !end do
    
    
                      !@ TEMPORARY ############################################
                      ! Extract ktilde
                      ktilde = f_malloc0((/n,n/),id='ktilde')
                      call extract_matrix(smatl, kerneltilde, neighbor(1:,kat), n, nmax, ktilde)
                      kp = f_malloc((/n,n/),id='kp')
                      !ii = 0
                      !do i=1,smats%nfvctr
                      !    if (neighbor(i,kat)) then
                      !        jj = 0
                      !        do j=1,smats%nfvctr
                      !            if (neighbor(j,kat)) then
                      !                icheck = icheck + 1
                      !                jj = jj + 1
                      !                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                      !                indl =  matrixindex_in_compressed(smatl, i, j)
                      !                if (indl>0) then
                      !                    ktilde(jj,ii) = kerneltilde(indl)
                      !                else
                      !                    ktilde(jj,ii) = 0.d0
                      !                end if
                      !            end if
                      !        end do
                      !    end if
                      !end do
    
                      ! Calculate ktilde * proj
                      call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, proj(1,1), n, 0.d0, kp(1,1), n)
                      tt = 0
                      do i=1,n
                          tt = tt + kp(i,i)
                      end do
                      !!if (bigdft_mpi%iproc==0) then
                      !!    do i=1,n
                      !!        do j=1,n
                      !!            write(*,'(a,2i5,3es13.3)') 'i, j, kt, proj, kp', i, j, ktilde(i,j), proj(j,i), kp(j,i)
                      !!        end do
                      !!    end do
                      !!    write(*,*) 'kkat, trace, sum(proj)', kkat, tt, sum(proj)
                      !!end if
                      charge_per_atom(kkat) = tt
                      !write(*,*) 'alpha, kkat, tt', alpha, kkat, tt
                      charge_total = charge_total + tt
                      call f_free(proj)
                      call f_free(ktilde)
                      call f_free(kp)
    
    
                  end do
    
    
                  !!if (final) exit kT_loop
    
                  !!charge_net = 0.d0
                  !!do iat=1,at%astruct%nat
                  !!    charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
                  !!end do
                  !!!!if (bigdft_mpi%iproc==0) then
                  !!!!    call yaml_map('kT, ef, net_charge',(/kT,ef,charge_net/))
                  !!!!end if
                  !!if (abs(charge_net)<1.d0 .or. ikT==100) then
                  !!    final = .true.
                  !!else if (charge_net<0.d0) then
                  !!    kT = kT*0.95d0
                  !!    !ef = ef + 1.d-3
                  !!else if (charge_net>0.d0) then
                  !!    kT = kT*1.05d0
                  !!    !ef = ef - 1.d-3
                  !!end if
    
              !!end do kT_loop

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
              charge_net = 0.d0
              do iat=1,at%astruct%nat
                  charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
              end do
              if (bigdft_mpi%iproc==0) then
                  !write(*,*) 'net charge', charge_net
                  call yaml_mapping_open(flow=.true.)
                  call yaml_map('alpha',alpha,fmt='(es12.4)')
                  call yaml_map('ef',ef,fmt='(es12.4)')
                  call yaml_map('net charge',charge_net,fmt='(es12.4)')
                  call yaml_map('bisection bounds ok',(/bound_low_ok,bound_up_ok/))
                  call yaml_mapping_close()
              end if


              if (abs(charge_net)<convergence_criterion) then
                  if (bigdft_mpi%iproc==0) then
                      call yaml_sequence_close()
                      call yaml_map('number of states to be occupied (without smearing)',iq)
                      call yaml_map('Pseudo Fermi level for occupations',ef)
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
                  !exit alpha_loop
                  exit ef_loop
              end if


              ! If we are still searching the boundaries for the bisection...
              if (.not.bound_up_ok) then
                  if (charge_net<0.d0) then
                      ! this is an upper bound
                      ef_up = ef
                      bound_up_ok = .true.
                  else
                      !ef_low = 0.5d0*ef
                      ef_low = ef + 0.5d0
                  end if
                  cycle ef_loop
              else if (.not.bound_low_ok) then
                  if (charge_net>0.d0) then
                      ! this is a lower bound
                      ef_low = ef
                      bound_low_ok = .true.
                  else
                      !ef_low = 2.0d0*ef
                      ef_low = ef - 0.5d0
                  end if
                  cycle ef_loop
              end if

              if (charge_net>0.d0) then
                  ! Too few electrons, i.e. confinement should be smaller
                  !alpha = alpha*0.80
                  ef_low = ef
              else if (charge_net<0.d0) then
                  ! Too many electrons, i.e. confinement should be larger
                  !alpha = alpha*1.2
                  ef_up = ef
              end if

          end do ef_loop



      !end do alpha_loop

      if (bigdft_mpi%iproc==0) then
          call yaml_mapping_close()
      end if

      call deallocate_matrices(ovrlp_onehalf_(1))
      call f_free(tmpmat1)
      !call f_free(tmpmat2)
      call f_free(kerneltilde)
      call f_free(coeff_all)
      call f_free(ilup)
      call f_free(n_all)
      call f_free(ovrlp_minusonehalf)
      call f_free(alpha_arr)
      call f_free(alpha_arr_bounds)
    
      !if (bigdft_mpi%iproc==0) then
      !    call yaml_mapping_close()
      !end if
    
    
    
      if (bigdft_mpi%iproc==0) then
          call write_partial_charges(at, charge_per_atom, write_gnuplot=.false.)
      end if

      call f_free(charge_per_atom)
      call f_free(neighbor)
      call f_free(eval_all)
      call f_free(id_all)
      call f_free(ovrlp_onehalf_all)
      call f_free(com)

      call f_release_routine()

  end subroutine projector_for_charge_analysis


    subroutine projector_for_charge_analysis_old(at, smats, smatm, smatl, &
               ovrlp_, ham_, kernel_, rxyz, calculate_centers, &
               lzd, nphirdim, psi, orbs)
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
      type(matrices),intent(inout) :: ovrlp_ !< should be intent(in)...
      type(matrices),intent(in) :: ham_, kernel_
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      logical,intent(in) :: calculate_centers
      type(local_zone_descriptors),intent(in),optional :: lzd
      integer,intent(in),optional :: nphirdim
      real(kind=8),dimension(:),intent(in),optional :: psi
      type(orbitals_data),intent(in),optional :: orbs

      ! Local variables
      integer :: kat, iat, jat, i, j, ii, jj, icheck, n, indm, inds, ntot, ist, ind, iq, itype, ieval, ij, nmax, indl, lwork
      integer :: k, l, iatold, isat, natp, kkat, istot, ntotp, i1, i2, i3, is1, ie1, is2, ie2, is3, ie3, j1, j2, j3, ikT, info
      integer :: ialpha
      real(kind=8) :: r2, cutoff2, rr2, tt, ef, q, occ, max_error, mean_error, rr2i, rr2j, ttxi, ttyi, ttzi, ttxj, ttyj, ttzj
      real(kind=8) :: tti, ttj, charge_net, charge_total
      real(kind=8) :: xi, xj, yi, yj, zi, zj, ttx, tty, ttz, xx, yy, zz, x, y, z
      real(kind=8),dimension(:),allocatable :: work
      real(kind=8),dimension(:,:),allocatable :: com
      real(kind=8),dimension(:,:),allocatable :: ham, ovrlp, proj, ovrlp_tmp, ovrlp_minusonehalf, kp, ktilde
      real(kind=8),dimension(:,:,:),allocatable :: coeff_all, ovrlp_onehalf_all
      integer,dimension(:,:,:,:),allocatable :: ilup
      real(kind=8),dimension(:),allocatable :: eval, eval_all, ovrlp_large, tmpmat1, tmpmat2, kerneltilde, charge_per_atom
      real(kind=8),dimension(:,:,:),allocatable :: tmpmat2d
      integer,dimension(:),allocatable :: id_all, n_all, itmparr
      real(kind=8),dimension(3) :: rr
      logical,dimension(:,:),allocatable :: neighbor
      type(matrices),dimension(1) :: ovrlp_onehalf_
      logical :: perx, pery, perz, final, bound_low_ok, bound_up_ok
      !real(kind=8),parameter :: kT = 5.d-2
      real(kind=8) :: kT
      !real(kind=8),parameter :: alpha = 5.d-1
      real(kind=8) :: alpha, alpha_up, alpha_low, convergence_criterion


      call f_routine(id='projector_for_charge_analysis')

      kT = 1.d-2

      ! Convergence criterion: one million-th of the total charge
      tt = 0.d0
      do iat=1,at%astruct%nat
          tt = tt + real(at%nelpsp(at%astruct%iatype(iat)),kind=8)
      end do
      convergence_criterion = 1.d-6*abs(tt)

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


      ! Calculate S^1/2
      ovrlp_onehalf_(1) = matrices_null()
      ovrlp_onehalf_(1)%matrix_compr = sparsematrix_malloc_ptr(smatl, iaction=SPARSE_TASKGROUP, id='ovrlp_onehalf_(1)%matrix_compr')
      call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, 1020, 1, (/2/), -1, &
            imode=1, ovrlp_smat=smats, inv_ovrlp_smat=smatl, &
            ovrlp_mat=ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_(1), &
            check_accur=.true., max_error=max_error, mean_error=mean_error)

      ! Calculate S^1/2 * K * S^1/2 = Ktilde
      tmpmat1 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat1')
      !tmpmat2 = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='tmpmat2')
      kerneltilde = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=smatl, id='kerneltilde')
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           kernel_%matrix_compr, ovrlp_onehalf_(1)%matrix_compr, tmpmat1)
      call matrix_matrix_mult_wrapper(bigdft_mpi%iproc, bigdft_mpi%nproc, smatl, &
           ovrlp_onehalf_(1)%matrix_compr, tmpmat1, kerneltilde)


      ! Determine the periodicity...
      perx=(smats%geocode /= 'F')
      pery=(smats%geocode == 'P')
      perz=(smats%geocode /= 'F')
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
      do kat=1,natp
          iatold = 0
          kkat = kat + isat
          n = 0
          do i=1,smats%nfvctr
               iat = smats%on_which_atom(i)
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
      id_all = f_malloc0(ntot,id='id_all')
      coeff_all = f_malloc((/nmax,nmax,natp/),id='coeff_all')
      ovrlp_onehalf_all = f_malloc((/nmax,nmax,natp/),id='ovrlp_onehalf_all')
      ovrlp_minusonehalf = f_malloc((/nmax,nmax/),id='ovrlp_minusonehalf')
      ilup = f_malloc((/2,nmax,nmax,natp/),id='ilup')
      n_all = f_malloc(natp,id='n_all')


      ! Centers of the support functions
      com = f_malloc0((/3,smats%nfvctr/),id='com')
      if (calculate_centers) then
          if (orbs%norb>0) then
              !call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, tmb%collcom_sr%ndimpsi_c, &
              !     orbs%norb, orbs%norbp, orbs%isorb, orbs%in_which_locreg, lzd, com(1:,orbs%isorb+1:))
              call supportfunction_centers(at%astruct%nat, rxyz, size(psi), psi, nphirdim, &
                   orbs%norb, orbs%norbp, orbs%isorb, orbs%inwhichlocreg, lzd, com(1:,orbs%isorb+1:))
              if (bigdft_mpi%nproc>1) then
                  call mpiallred(com, mpi_sum, comm=bigdft_mpi%mpi_comm)
              end if
          end if
      else
          do i=1,smats%nfvctr
              iat = smats%on_which_atom(i)
              com(1:3,i) = rxyz(1:3,iat)
          end do
      end if


      charge_per_atom = f_malloc0(at%astruct%nat,id='charge_per_atom')


      ! Calculate how many states should be included
      q = 0.d0
      do iat=1,at%astruct%nat
          itype = at%astruct%iatype(iat)
          q = q + ceiling(0.5d0*real(at%nelpsp(itype),kind=8))
      end do
      iq = nint(q)
      if (bigdft_mpi%iproc==0) then
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

      ialpha = 0
      alpha_loop: do! ialpha=1,100

          ialpha = ialpha + 1

          if (bigdft_mpi%iproc==0) then
              call yaml_sequence(advance='no')
          end if

          if (.not.bound_low_ok) then
              alpha = alpha_low
          else if (.not.bound_up_ok) then
              alpha = alpha_up
          else
              alpha = 0.5d0*(alpha_low+alpha_up)
          end if

          if (kkat==2) then
              !O
              alpha = 5.51459534d0*2.d-2
          else
              !H
              alpha = 63.2d0*2.d-2
          end if

          charge_net = 0.d0
          call f_zero(eval_all)
          call f_zero(id_all)

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
              call extract_matrix(smats, ovrlp_%matrix_compr, neighbor(1:,kat), n, nmax, ovrlp, ilup)
              call extract_matrix(smatm, ham_%matrix_compr, neighbor(1:,kat), n, nmax, ham)
              !!icheck = 0
              !!ii = 0
              !!do i=1,smats%nfvctr
              !!    if (neighbor(i,kat)) then
              !!        jj = 0
              !!        do j=1,smats%nfvctr
              !!            if (neighbor(j,kat)) then
              !!                icheck = icheck + 1
              !!                jj = jj + 1
              !!                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
              !!                inds =  matrixindex_in_compressed(smats, i, j)
              !!                if (inds>0) then
              !!                    ovrlp(jj,ii) = ovrlp_%matrix_compr(inds)
              !!                else
              !!                    ovrlp(jj,ii) = 0.d0
              !!                end if
              !!                indm =  matrixindex_in_compressed(smatm, i, j)
              !!                if (indm>0) then
              !!                    ham(jj,ii) = ham_%matrix_compr(indm)
              !!                else
              !!                    ham(jj,ii) = 0.d0
              !!                end if
              !!            end if
              !!        end do
              !!    end if
              !!end do
              !!if (icheck>n**2) then
              !!    call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
              !!        &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
              !!end if


              ! Calculate ovrlp^1/2 and ovrlp^-1/2. The last argument is wrong, clean this.
              ovrlp_tmp = f_malloc((/n,n/),id='ovrlp_tmp')
              call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
              call overlap_plus_minus_one_half_exact(1, n, -1, .true., ovrlp_tmp, smats)
              do i=1,n
                  call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_onehalf_all(1,i,kat), 1)
              end do
              call f_memcpy(src=ovrlp, dest=ovrlp_tmp)
              call overlap_plus_minus_one_half_exact(1, n, -1, .false., ovrlp_tmp, smats)
              do i=1,n
                  call vcopy(n, ovrlp_tmp(1,i), 1, ovrlp_minusonehalf(1,i), 1)
              end do
              call f_free(ovrlp_tmp)
    
              ! Calculate S^-1/2 * H * S^-1/2
              tmpmat2d = f_malloc((/n,n,1/),id='tmppmat2d')
              call gemm('n', 'n', n, n, n, 1.d0, ham(1,1), n, ovrlp_minusonehalf(1,1), nmax, 0.d0, tmpmat2d(1,1,1), n)
              call gemm('n', 'n', n, n, n, 1.d0, ovrlp_minusonehalf(1,1), nmax, tmpmat2d(1,1,1), n, 0.d0, ham(1,1), n)
              call f_free(tmpmat2d)

              ! Add the penalty term
              call add_penalty_term(smats%geocode, smats%nfvctr, neighbor(1:,kat), rxyz(1:,kkat), &
                   at%astruct%cell_dim, com, alpha, n, ovrlp, ham)
              !!icheck = 0
              !!ii = 0
              !!do i=1,smats%nfvctr
              !!    if (neighbor(i,kat)) then
              !!        jj = 0
              !!        do j=1,smats%nfvctr
              !!            if (neighbor(j,kat)) then
              !!                icheck = icheck + 1
              !!                jj = jj + 1
              !!                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
              !!                if (i==j) then
              !!                    rr2 = huge(rr2)
              !!                    do i3=is3,ie3
              !!                        z = rxyz(3,kkat) + i3*at%astruct%cell_dim(3)
              !!                        ttz = (com(3,i)-z)**2
              !!                        do i2=is2,ie2
              !!                            y = rxyz(2,kkat) + i2*at%astruct%cell_dim(2)
              !!                            tty = (com(2,i)-y)**2
              !!                            do i1=is1,ie1
              !!                                x = rxyz(1,kkat) + i1*at%astruct%cell_dim(1)
              !!                                ttx = (com(1,i)-x)**2
              !!                                tt = ttx + tty + ttz
              !!                                if (tt<rr2) then
              !!                                    rr2 = tt
              !!                                end if
              !!                            end do
              !!                        end do
              !!                    end do
              !!                    ham(jj,ii) = ham(jj,ii) + alpha*rr2**3*ovrlp(jj,ii)
              !!                end if
              !!                ilup(1,jj,ii,kat) = j
              !!                ilup(2,jj,ii,kat) = i
              !!            end if
              !!        end do
              !!    end if
              !!end do
              !!if (icheck>n**2) then
              !!    call f_err_throw('icheck('//adjustl(trim(yaml_toa(icheck)))//') > n**2('//&
              !!        &adjustl(trim(yaml_toa(n**2)))//')',err_name='BIGDFT_RUNTIME_ERROR')
              !!end if
    
    
              !!call diagonalizeHamiltonian2(bigdft_mpi%iproc, n, ham, ovrlp, eval)
              lwork = 10*n
              work = f_malloc(lwork,id='work')
              call syev('v', 'l', n, ham(1,1), n, eval(1), work(1), lwork, info)
              call f_free(work)
              do i=1,n
                  ii = ist + i
                  eval_all(istot+ii) = eval(i)
                  write(*,*) 'kkat, alpha, i, eval', kkat, alpha, i, eval(i)
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
          call order_eigenvalues(ntot, eval_all, id_all)
          !do i=1,ntot
          !    ! add i-1 since we are only searching in the subarray
          !    ind = minloc(eval_all(i:ntot),1) + (i-1)
          !    tt = eval_all(i)
          !    eval_all(i) = eval_all(ind)
          !    eval_all(ind) = tt
          !    ii = id_all(i)
          !    id_all(i) = id_all(ind)
          !    id_all(ind) = ii
          !end do
        
        
    
    
          !!! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
          !!ef = eval_all(1)
          !!do
          !!    ef = ef + 1.d-3
          !!    occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
          !!    if (abs(occ-1.d0)<1.d-8) exit
          !!end do
          !!if (bigdft_mpi%iproc==0) then
          !!    call yaml_map('Pseudo Fermi level for occupations',ef)
          !!end if
        
          !!final = .true.
          !!ikT = 0
          !!kT_loop: do
    
              !ikT = ikT + 1
    
              call f_zero(charge_per_atom)
    
              ! Determine the "Fermi level" such that the iq-th state is still fully occupied even with a smearing
              ef = eval_all(1)
              do
                  ef = ef + 1.d-3
                  !occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
                  occ = 1.d0/(1.d0+safe_exp( (eval_all(iq)-ef)*(1.d0/kT) ) )
                  !occ = 1.d0/(1.d0+safe_exp( (-0.2627921562586862d0-ef)*(1.d0/kT) ) ) ! $ 4-1
                  !occ = 1.d0/(1.d0+safe_exp( (-0.2637375288617387d0-ef)*(1.d0/kT) ) ) ! $ 9-4
                  if (abs(occ-1.d0)<1.d-8) exit
              end do

              !ef = -0.00d0
        
              ! Calculate the projector. First for each single atom, then insert it into the big one.
              charge_total = 0.d0
              do kat=1,natp
                  kkat = kat + isat
                  n = n_all(kat)
                  proj = f_malloc0((/n,n/),id='proj')
                  call calculate_projector(n, ntot, nmax, kkat, id_all, eval_all, &
                       coeff_all(1:,1:,kat), ef, kT, proj)
                  !ij = 0
                  !do ieval=1,ntot
                  !    if (id_all(ieval)/=kkat) cycle
                  !    ij = ij + 1
                  !    occ = 1.d0/(1.d0+safe_exp( (eval_all(ieval)-ef)*(1.d0/kT) ) )
                  !    do i=1,n
                  !        do j=1,n
                  !            proj(j,i) = proj(j,i) + occ*coeff_all(j,ij,kat)*coeff_all(i,ij,kat)
                  !        end do
                  !   end do
                  !end do
                  !tt = 0.d0
                  !do i=1,n
                  !    tt = tt + proj(i,i)
                  !end do
    
    
                  !@ TEMPORARY ############################################
                  ! Extract ktilde
                  ktilde = f_malloc0((/n,n/),id='ktilde')
                  call extract_matrix(smatl, kerneltilde, neighbor(1:,kat), n, nmax, ktilde)
                  kp = f_malloc((/n,n/),id='kp')
                  !ii = 0
                  !do i=1,smats%nfvctr
                  !    if (neighbor(i,kat)) then
                  !        jj = 0
                  !        do j=1,smats%nfvctr
                  !            if (neighbor(j,kat)) then
                  !                icheck = icheck + 1
                  !                jj = jj + 1
                  !                if (jj==1) ii = ii + 1 !new column if we are at the first line element of a a column
                  !                indl =  matrixindex_in_compressed(smatl, i, j)
                  !                if (indl>0) then
                  !                    ktilde(jj,ii) = kerneltilde(indl)
                  !                else
                  !                    ktilde(jj,ii) = 0.d0
                  !                end if
                  !            end if
                  !        end do
                  !    end if
                  !end do
    
                  ! Calculate ktilde * proj
                  call gemm('n', 'n', n, n, n, 1.d0, ktilde(1,1), n, proj(1,1), n, 0.d0, kp(1,1), n)
                  tt = 0
                  do i=1,n
                      tt = tt + kp(i,i)
                  end do
                  !!if (bigdft_mpi%iproc==0) then
                  !!    do i=1,n
                  !!        do j=1,n
                  !!            write(*,'(a,2i5,3es13.3)') 'i, j, kt, proj, kp', i, j, ktilde(i,j), proj(j,i), kp(j,i)
                  !!        end do
                  !!    end do
                  !!    write(*,*) 'kkat, trace, sum(proj)', kkat, tt, sum(proj)
                  !!end if
                  charge_per_atom(kkat) = tt
                  !write(*,*) 'alpha, kkat, tt', alpha, kkat, tt
                  charge_total = charge_total + tt
                  call f_free(proj)
                  call f_free(ktilde)
                  call f_free(kp)
    
    
              end do
    
    
              !!if (final) exit kT_loop
    
              !!charge_net = 0.d0
              !!do iat=1,at%astruct%nat
              !!    charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
              !!end do
              !!!!if (bigdft_mpi%iproc==0) then
              !!!!    call yaml_map('kT, ef, net_charge',(/kT,ef,charge_net/))
              !!!!end if
              !!if (abs(charge_net)<1.d0 .or. ikT==100) then
              !!    final = .true.
              !!else if (charge_net<0.d0) then
              !!    kT = kT*0.95d0
              !!    !ef = ef + 1.d-3
              !!else if (charge_net>0.d0) then
              !!    kT = kT*1.05d0
              !!    !ef = ef - 1.d-3
              !!end if
    
          !!end do kT_loop

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
          charge_net = 0.d0
          do iat=1,at%astruct%nat
              charge_net = charge_net -(charge_per_atom(iat)-real(at%nelpsp(at%astruct%iatype(iat)),kind=8))
          end do
          if (bigdft_mpi%iproc==0) then
              !write(*,*) 'net charge', charge_net
              call yaml_mapping_open(flow=.true.)
              call yaml_map('alpha',alpha,fmt='(es12.4)')
              call yaml_map('net charge',charge_net,fmt='(es12.4)')
              call yaml_map('bisection bounds ok',(/bound_low_ok,bound_up_ok/))
              call yaml_mapping_close()
          end if


          if (abs(charge_net)<convergence_criterion .or. ialpha==50 .or. .true.) then
              if (bigdft_mpi%iproc==0) then
                  call yaml_sequence_close()
                  call yaml_map('number of states to be occupied (without smearing)',iq)
                  call yaml_map('Pseudo Fermi level for occupations',ef)
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
              exit alpha_loop
          end if

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
      end if

      call deallocate_matrices(ovrlp_onehalf_(1))
      call f_free(tmpmat1)
      !call f_free(tmpmat2)
      call f_free(kerneltilde)
      call f_free(coeff_all)
      call f_free(ilup)
      call f_free(n_all)
      call f_free(ovrlp_minusonehalf)
    
      !if (bigdft_mpi%iproc==0) then
      !    call yaml_mapping_close()
      !end if
    
    
    
      if (bigdft_mpi%iproc==0) then
          call write_partial_charges(at, charge_per_atom, write_gnuplot=.false.)
      end if

      call f_free(charge_per_atom)
      call f_free(neighbor)
      call f_free(eval_all)
      call f_free(id_all)
      call f_free(ovrlp_onehalf_all)
      call f_free(com)

      call f_release_routine()

  end subroutine projector_for_charge_analysis_old
