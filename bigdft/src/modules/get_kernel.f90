module get_kernel
  implicit none
  private

  !> Public routines
  public :: get_coeff
  public :: reconstruct_kernel
  public :: renormalize_kernel
  public :: reorthonormalize_coeff
  public :: calculate_gap_FOE
  public :: coeff_weight_analysis
  public :: order_coeffs_by_energy

  contains


    subroutine get_coeff(iproc,nproc,scf_mode,orbs,at,rxyz,denspot,GPU,infoCoeff,&
        energs,nlpsp,SIC,tmb,fnrm,calculate_overlap_matrix,invert_overlap_matrix,&
        calculate_pspandkin,communicate_phi_for_lsumrho,&
        calculate_ham,extra_states,itout,it_scc,it_cdft,order_taylor,max_inversion_error,&
        calculate_KS_residue,calculate_gap,energs_work,remove_coupling_terms,factor,tel,occopt,&
        pexsi_npoles,pexsi_nproc_per_pole,pexsi_mumin,pexsi_mumax,pexsi_mu,&
        pexsi_DeltaE,pexsi_temperature, pexsi_tol_charge, pexsi_np_sym_fact, &
        pexsi_do_inertia_count, pexsi_max_iter, pexsi_verbosity, &
        convcrit_dmin,nitdmin,curvefit_dmin,ldiis_coeff,reorder,cdft,updatekernel,hphi_pspandkin,eproj,ekin)
      use module_base
      use module_types
      use module_interfaces, only: LocalHamiltonianApplication, SynchronizeHamiltonianApplication
      use Poisson_Solver, except_dp => dp, except_gp => gp
      use diis_sd_optimization
      use yaml_output
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized, start_onesided_communication, communicate_basis_for_density_collective
      use constrained_dft
      use rhopotential, only: full_local_potential
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, sparsematrix_malloc, &
                                   DENSE_FULL, DENSE_PARALLEL, DENSE_MATMUL, assignment(=), SPARSE_FULL
      use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                              extract_taskgroup_inplace, uncompress_matrix_distributed2, gather_matrix_from_taskgroups, &
                              extract_taskgroup, uncompress_matrix2, &
                              write_sparsematrix, delete_coupling_terms, transform_sparse_matrix, &
                              max_asymmetry_of_matrix, diagonalizeHamiltonian2
      use sparsematrix_init, only: sparsebigdft_to_ccs
      use transposed_operations, only: calculate_overlap_transposed
      use parallel_linalg, only: dsygv_parallel
      use matrix_operations, only: deviation_from_unity_parallel
      !use foe, only: fermi_operator_expansion_new
      use sparsematrix_highlevel, only: matrix_fermi_operator_expansion
      use public_enums
      use locregs_init, only: small_to_large_locreg
      use locreg_operations, only: confpot_data
      use foe_base, only: foe_data_get_real
      use pexsi, only: pexsi_wrapper !pexsi_driver
      use coeffs, only: get_coeffs_diagonalization, calculate_kernel_and_energy
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, scf_mode, itout, it_scc, it_cdft, occopt
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(orbitals_data),intent(inout) :: orbs
      type(atoms_data),intent(in) :: at
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers),intent(inout) :: GPU
      integer,intent(out) :: infoCoeff
      type(energy_terms),intent(inout) :: energs
      real(kind=8),intent(inout) :: fnrm
      type(DFT_PSP_projectors),intent(inout) :: nlpsp
      type(SIC_data),intent(in) :: SIC
      type(DFT_wavefunction),intent(inout) :: tmb
      logical,intent(in):: calculate_overlap_matrix, invert_overlap_matrix, calculate_pspandkin
      logical,intent(in):: communicate_phi_for_lsumrho
      logical,intent(in) :: calculate_ham, calculate_KS_residue, calculate_gap
      type(work_mpiaccumulate),intent(inout) :: energs_work
      logical,intent(in) :: remove_coupling_terms
      real(kind=8), intent(in) :: factor, tel
      integer,intent(in) :: pexsi_npoles, pexsi_nproc_per_pole, pexsi_np_sym_fact, pexsi_max_iter, pexsi_verbosity
      logical,intent(in) :: pexsi_do_inertia_count
      real(kind=8),intent(in) :: pexsi_mumin,pexsi_mumax,pexsi_mu,pexsi_DeltaE,pexsi_temperature, pexsi_tol_charge
      type(DIIS_obj),intent(inout),optional :: ldiis_coeff ! for dmin only
      integer, intent(in), optional :: nitdmin ! for dmin only
      real(kind=gp), intent(in), optional :: convcrit_dmin ! for dmin only
      logical, intent(in), optional :: curvefit_dmin ! for dmin only
      type(cdft_data),intent(inout),optional :: cdft
      integer, intent(in) :: extra_states
      logical, optional, intent(in) :: reorder
      logical, optional, intent(in) :: updatekernel
      ! The foloowing array contains the psp and kinetic part of the Hamiltonian appllication.. Can be used to spped up the code in case phi does not change between calls, but only the potential
      real(kind=8),dimension(tmb%ham_descr%npsidim_orbs),intent(inout),optional :: hphi_pspandkin
      real(kind=8),intent(inout),optional :: eproj, ekin
    
      ! Local variables 
      integer :: iorb, info, ishift, ispin, ii, jorb, i, ishifts, ishiftm, jproc, j
      real(kind=8),dimension(:),allocatable :: hpsit_c, hpsit_f, eval, tmparr1, tmparr2, tmparr, ovrlp_large, ham_large, eval2
      real(kind=8),dimension(:,:),allocatable :: ovrlp_fullp, tempmat
      real(kind=8),dimension(:,:,:),allocatable :: matrixElements
      type(confpot_data),dimension(:),allocatable :: confdatarrtmp
      logical :: update_kernel, auxiliary_arguments_present
      character(len=*),parameter :: subname='get_coeff'
      real(kind=gp) :: tmprtr
      real(kind=8) :: max_deviation, mean_deviation, KSres, max_deviation_p,  mean_deviation_p, maxdiff, tt
      real(kind=8) :: asymm_S, asymm_H, asymm_K
      integer,dimension(:),allocatable :: row_ind, col_ptr, n3p
    
      call f_routine(id='get_coeff')
    
      if(calculate_ham) then
          !!energs_work = work_mpiaccumulate_null()
          !!energs_work%ncount = 4
          !!call allocate_work_mpiaccumulate(energs_work)
    
          call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
          n3p = f_malloc(0.to.nproc-1,id='n3p')
          do jproc=0,nproc-1
              n3p(jproc) = max(denspot%dpbox%nscatterarr(jproc,2),1)
          end do
          call start_onesided_communication(iproc, nproc, denspot%dpbox%mesh%ndims(1), denspot%dpbox%mesh%ndims(2), &
               n3p, denspot%rhov, &
               tmb%ham_descr%comgp%nrecvbuf*tmb%ham_descr%comgp%nspin, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, &
               tmb%ham_descr%lzd)
          call f_free(n3p)
      end if
    
      ! Option to only calculate the energy without updating the kernel
      if (present(updatekernel)) then
          update_kernel=updatekernel
      else
          update_kernel=.true.
      end if
    
      ! This is simply not yet implemented
      if (.not.update_kernel .and. scf_mode==LINEAR_FOE) then
          stop 'ERROR: for the moment, update_kernel must be true for FOE'
      end if
    
       if (iproc==0) call yaml_mapping_open('Kernel update')
    
    
      ! Calculate the Hamiltonian matrix if it is not already present.
      if(calculate_ham) then
    
          !!call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
          !!call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
          !!     tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)
    
          if (iproc==0) then
              call yaml_map('Hamiltonian application required',.true.)
          end if
    
          allocate(confdatarrtmp(tmb%orbs%norbp))
          call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)
    
          call small_to_large_locreg(iproc, tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%isorb, tmb%orbs%inwhichlocreg, &
               tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
               tmb%psi, tmb%ham_descr%psi)
    
          ! Check the optional arguments
          if (any((/present(hphi_pspandkin),present(eproj),present(ekin)/))) then
              if (all((/present(hphi_pspandkin),present(eproj),present(ekin)/))) then
                  auxiliary_arguments_present = .true.
              else
                  call f_err_throw('The arguments hphi_pspandkin, eproj and ekin miust be present at the same time',&
                       err_name='BIGDFT_RUNTIME_ERROR')
              end if
          else
              auxiliary_arguments_present = .false.
          end if
    
          if (.not.calculate_pspandkin) then
              if (.not.auxiliary_arguments_present) then
                  call f_err_throw('The optionals arguments hphi_pspandkin, eproj and ekin must be present &
                       &when calculate_pspandkin is wrong')
              end if
          end if
    
    
          if(calculate_pspandkin) then
    
              if (tmb%ham_descr%npsidim_orbs > 0) call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
    
              call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
                   tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,tmb%hpsi,energs%eproj,tmb%paw)
              ! only kinetic as waiting for communications
              call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
                   tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,&
                   & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,3,denspot%xc,&
                   & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
                   & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
              if (auxiliary_arguments_present) then
                  if (tmb%ham_descr%npsidim_orbs > 0) then
                      call f_memcpy(src=tmb%hpsi, dest=hphi_pspandkin)
                      eproj = energs%eproj
                      ekin = energs%ekin
                  end if
              end if
              if (iproc==0) then
                  call yaml_map('PSP and kinetic Hamiltonian application','recalculated')
              end if
          else
              call f_memcpy(src=hphi_pspandkin, dest=tmb%hpsi)
              energs%eproj = eproj
              energs%ekin = ekin
              if (iproc==0) then
                  call yaml_map('PSP and kinetic Hamiltonian application','from memory')
              end if
          end if
          !!do i=1,tmb%ham_descr%comgp%nrecvbuf
          !!    write(8000+iproc,'(a,i8,es16.6)') 'i, recvbuf(i)', i, tmb%ham_descr%comgp%recvbuf(i)
          !!end do
          !call wait_p2p_communication(iproc, nproc, tmb%ham_descr%comgp)
          ! only potential
          call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
               & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,&
               & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
               & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
               & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
          call timing(iproc,'glsynchham1','ON')
          call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,&
               & tmb%orbs,tmb%ham_descr%lzd,GPU,denspot%xc,tmb%hpsi,&
               energs, energs_work)
          call timing(iproc,'glsynchham1','OF')
          deallocate(confdatarrtmp)
    
    
    
          !!if (iproc==0) write(*,'(a,5es20.12)') 'ekin, eh, epot, eproj, eex', &
          !!              energs%ekin, energs%eh, energs%epot, energs%eproj, energs%exc
    
          !DEBUG
          !if(iproc==0) then
          !  print *,'Ekin,Epot,Eproj,Eh,Exc,Evxc',energs%ekin,energs%epot,energs%eproj,energs%eh,energs%exc,energs%evxc
          !end if
          !END DEBUG
    
          call f_free_ptr(denspot%pot_work)
    
          ! Calculate the matrix elements <phi|H|phi>.
          if(.not.tmb%ham_descr%can_use_transposed) then
              call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   TRANSPOSE_FULL, tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
              tmb%ham_descr%can_use_transposed=.true.
          end if
    
          hpsit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
          hpsit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_FULL, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
    
    
          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%smat(2), tmb%linmat%auxm, tmb%linmat%ham_)
          !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
    
    
    !!  call diagonalize_subset(iproc, nproc, tmb%orbs, tmb%linmat%smat(1), tmb%linmat%ovrlp_, tmb%linmat%smat(2), tmb%linmat%ham_)
    !!  if (iproc==0) then
    !!      do iorb=1,tmb%orbs%norb
    !!          write(*,*) 'iorb, tmb%orbs%eval(iorb)',iorb,tmb%orbs%eval(iorb)
    !!      end do
    !!  end if
    
      else
          !!if(iproc==0) write(*,*) 'No Hamiltonian application required.'
          if (iproc==0) then
              call yaml_map('Hamiltonian application required',.false.)
          end if
      end if
    
    
      ! Calculate the overlap matrix if required.
      if(calculate_overlap_matrix) then
          if(.not.tmb%can_use_transposed) then
              call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                   TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
              tmb%can_use_transposed=.true.
          end if
    
          if (iproc==0) call yaml_map('calculate overlap matrix',.true.)
          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
               tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%ovrlp_)
          !call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
      end if
    
      ovrlp_fullp = sparsematrix_malloc(tmb%linmat%smat(3),iaction=DENSE_PARALLEL,id='ovrlp_fullp')
      max_deviation=0.d0
      mean_deviation=0.d0
      do ispin=1,tmb%linmat%smat(1)%nspin
          ishift=(ispin-1)*tmb%linmat%smat(1)%nvctrp_tg
          call uncompress_matrix_distributed2(iproc, tmb%linmat%smat(1), DENSE_PARALLEL, &
               tmb%linmat%ovrlp_%matrix_compr(ishift+1:), ovrlp_fullp)
          call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%smat(1)%nfvctr, tmb%linmat%smat(1)%nfvctrp, &
               tmb%linmat%smat(1)%isfvctr, ovrlp_fullp, &
               tmb%linmat%smat(1), max_deviation_p, mean_deviation_p)
          !!call deviation_from_unity_parallel_new(iproc, nproc, bigdft_mpi%mpi_comm, &
          !!     tmb%linmat%ovrlp_%matrix_compr(ishift+1:), tmb%linmat%smat(1), &
          !!     max_deviation_p, mean_deviation_p)
          max_deviation = max_deviation + max_deviation_p/real(tmb%linmat%smat(1)%nspin,kind=8)
          mean_deviation = mean_deviation + mean_deviation_p/real(tmb%linmat%smat(1)%nspin,kind=8)
      end do
      call f_free(ovrlp_fullp)
      !if (iproc==0) then
      !    call yaml_map('max dev from unity',max_deviation,fmt='(es9.2)')
      !    call yaml_map('mean dev from unity',mean_deviation,fmt='(es9.2)')
      !end if
    
      ! Post the p2p communications for the density. (must not be done in inputguess)
      if(communicate_phi_for_lsumrho) then
          call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
               tmb%orbs, tmb%psi, tmb%collcom_sr)
      end if
    
      ! CDFT: add V*w_ab to Hamiltonian here - assuming ham and weight matrix have the same sparsity...
      if (present(cdft)) then
         call timing(iproc,'constraineddft','ON')
         call daxpy(tmb%linmat%smat(2)%nvctr,cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmb%linmat%ham_%matrix_compr,1)
         call timing(iproc,'constraineddft','OF') 
      end if
    
      if (remove_coupling_terms) then
          call delete_coupling_terms(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%smmd, tmb%linmat%smat(1), tmb%linmat%ovrlp_%matrix_compr)
          call delete_coupling_terms(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%smmd, tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr)
          call delete_coupling_terms(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%smmd, tmb%linmat%smat(3), tmb%linmat%kernel_%matrix_compr)
      end if
      ! Calculate the asymmetry of S and H
      call max_asymmetry_of_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(1), tmb%linmat%ovrlp_%matrix_compr, asymm_S)
      call max_asymmetry_of_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr, asymm_H)
    
      !!if (scf_mode/=LINEAR_FOE .and. scf_mode/=LINEAR_PEXSI) then
      !!    tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(2), iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
      !!    !call yaml_map('Ham1com',tmb%linmat%ham_%matrix_compr)
      !!    !call yaml_map('ovrlpcom',tmb%linmat%ovrlp_%matrix_compr)
    
      !!    tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(1), iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      !!    !call uncompress_matrix(iproc, tmb%linmat%smat(2), &
      !!    !     inmat=tmb%linmat%ham_%matrix_compr, outmat=tmb%linmat%ham_%matrix)
      !!    !call uncompress_matrix(iproc, tmb%linmat%smat(1), &
      !!    !     inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
      !!    do ispin=1,tmb%linmat%smat(2)%nspin
      !!        ishifts = (ispin-1)*tmb%linmat%smat(1)%nvctrp_tg
      !!        ishiftm = (ispin-1)*tmb%linmat%smat(2)%nvctrp_tg
      !!        call f_zero(tmb%linmat%smat(2)%nfvctr**2, tmb%linmat%ham_%matrix(1,1,ispin))
      !!        tempmat = sparsematrix_malloc(tmb%linmat%smat(2), iaction=DENSE_PARALLEL, id='tempmat')
      !!        call uncompress_matrix_distributed2(iproc, tmb%linmat%smat(2), DENSE_PARALLEL, &
      !!             tmb%linmat%ham_%matrix_compr(ishiftm+1:ishiftm+tmb%linmat%smat(2)%nvctrp_tg), tempmat)
      !!        if (tmb%linmat%smat(2)%nfvctrp>0) then
      !!            call vcopy(tmb%linmat%smat(2)%nfvctr*tmb%linmat%smat(2)%nfvctrp, tempmat(1,1), 1, &
      !!                 tmb%linmat%ham_%matrix(1,tmb%linmat%smat(2)%isfvctr+1,ispin), 1)
      !!        end if
      !!        call f_free(tempmat)
      !!        if (nproc>1) then
      !!            call mpiallred(tmb%linmat%ham_%matrix(1,1,ispin), tmb%linmat%smat(2)%nfvctr**2, &
      !!                 mpi_sum, comm=bigdft_mpi%mpi_comm)
      !!        end if
    
      !!        call f_zero(tmb%linmat%smat(1)%nfvctr**2, tmb%linmat%ovrlp_%matrix(1,1,ispin))
      !!        tempmat = sparsematrix_malloc(tmb%linmat%smat(1), iaction=DENSE_PARALLEL, id='tempmat')
      !!        call uncompress_matrix_distributed2(iproc, tmb%linmat%smat(1), DENSE_PARALLEL, &
      !!             tmb%linmat%ovrlp_%matrix_compr(ishifts+1:), tempmat)
      !!        if (tmb%linmat%smat(2)%nfvctrp>0) then
      !!            call vcopy(tmb%linmat%smat(1)%nfvctr*tmb%linmat%smat(1)%nfvctrp, tempmat(1,1), 1, &
      !!                 tmb%linmat%ovrlp_%matrix(1,tmb%linmat%smat(1)%isfvctr+1,ispin), 1)
      !!        end if
      !!        call f_free(tempmat)
      !!        if (nproc>1) then
      !!            call mpiallred(tmb%linmat%ovrlp_%matrix(1,1,ispin), tmb%linmat%smat(1)%nfvctr**2, &
      !!                 mpi_sum, comm=bigdft_mpi%mpi_comm)
      !!        end if
      !!    end do
      !!end if
    
      if(scf_mode==LINEAR_MIXPOT_SIMPLE .or. scf_mode==LINEAR_MIXDENS_SIMPLE) then
          call get_coeffs_diagonalization(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%smat(2)%nfvctr, orbs%norbu, orbs%norbd, orbs%norb, tmb%orthpar%blocksize_pdsyev, &
               tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%coeff, &
               tmb%orbs%eval, orbs%eval, infoCoeff)
      else if (scf_mode==LINEAR_DIRECT_MINIMIZATION) then
         if(.not.present(ldiis_coeff)) &
              call f_err_throw('ldiis_coeff must be present for scf_mode==LINEAR_DIRECT_MINIMIZATION',&
              err_name='BIGDFT_RUNTIME_ERROR')
         ! call routine which updates coeffs for tmb%orbs%norb or orbs%norb depending on whether or not extra states are required
         if (iproc==0) call yaml_map('method','directmin')
         tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(1), iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
         tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(2), iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
         call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, tmb%linmat%smat(1), &
              tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
         call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, tmb%linmat%smat(2), &
              tmb%linmat%ham_%matrix_compr, tmb%linmat%ham_%matrix)
         if (extra_states>0) then
            call optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, convcrit_dmin, nitdmin, energs%ebs, &
                 curvefit_dmin, factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder, extra_states)
         else
            call optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, convcrit_dmin, nitdmin, energs%ebs, &
                 curvefit_dmin, factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder)
         end if
         call f_free_ptr(tmb%linmat%ovrlp_%matrix)
         call f_free_ptr(tmb%linmat%ham_%matrix)
      end if
    
      ! CDFT: subtract V*w_ab from Hamiltonian so that we are calculating the correct energy
      if (present(cdft)) then
         call timing(iproc,'constraineddft','ON')
         tmparr = sparsematrix_malloc(tmb%linmat%smat(2),iaction=SPARSE_FULL,id='tmparr')
         call gather_matrix_from_taskgroups(iproc, nproc, bigdft_mpi%mpi_comm, &
              tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr, tmparr)
         call daxpy(tmb%linmat%smat(2)%nvctr*tmb%linmat%smat(2)%nspin,-cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmparr,1)
         call extract_taskgroup(tmb%linmat%smat(2), tmparr, tmb%linmat%ham_%matrix_compr)
         call f_free(tmparr)
         call timing(iproc,'constraineddft','OF') 
      end if
    
      if (scf_mode/=LINEAR_FOE .and. scf_mode/=LINEAR_PEXSI) then
          call evaltoocc(iproc,nproc,.false.,tel,orbs,occopt)
          if (scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
             call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),tmb%linmat%smat(2), &
                  tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
                  tmb%coeff, orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup, update_kernel)
          else if (present(cdft)) then
             ! for directmin we have the kernel already, but only the CDFT function not actual energy for CDFT
             call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),tmb%linmat%smat(2), &
                  tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
                  tmb%coeff,orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup, .false.)
          end if
    
      else ! foe or pexsi
    
          !same as for directmin/diag
          ! CDFT: add V*w_ab to Hamiltonian here - assuming ham and weight matrix have the same sparsity...
          if (present(cdft)) then
             call timing(iproc,'constraineddft','ON')
             call daxpy(tmb%linmat%smat(2)%nvctr,cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,&
                  tmb%linmat%ham_%matrix_compr,1)
             call timing(iproc,'constraineddft','OF') 
          end if
    
          if (calculate_gap) then
              tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(2), &
                  iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
              !!tmparr = sparsematrix_malloc(tmb%linmat%smat(2),iaction=SPARSE_FULL,id='tmparr')
              !!call gather_matrix_from_taskgroups(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr, tmparr)
              call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
                   tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr, tmb%linmat%ham_%matrix)
              !!call f_free(tmparr)
              tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(1), &
                  iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
              call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
                   tmb%linmat%smat(1), tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
              ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
              matrixElements=f_malloc((/tmb%orbs%norb,tmb%orbs%norb,2/),id='matrixElements')
              !SM: need to fix the spin here
              call vcopy(tmb%orbs%norb**2, tmb%linmat%ham_%matrix(1,1,1), 1, matrixElements(1,1,1), 1)
              call vcopy(tmb%orbs%norb**2, tmb%linmat%ovrlp_%matrix(1,1,1), 1, matrixElements(1,1,2), 1)
              call diagonalizeHamiltonian2(iproc, nproc, bigdft_mpi%mpi_comm, &
                   tmb%orthpar%blocksize_pdsyev, tmb%orbs%norb, &
                   matrixElements(1,1,1), matrixElements(1,1,2), tmb%orbs%eval)
              if (iproc==0) call yaml_map('gap',tmb%orbs%eval(orbs%norb+1)-tmb%orbs%eval(orbs%norb))
              if (iproc==0) call yaml_map('lowest eigenvalue',tmb%orbs%eval(1))
              if (iproc==0) call yaml_map('highest eigenvalue',tmb%orbs%eval(tmb%orbs%norb))
              call f_free(matrixElements)
              call f_free_ptr(tmb%linmat%ham_%matrix)
              call f_free_ptr(tmb%linmat%ovrlp_%matrix)
          end if
    
          tmprtr=0.d0
    
          if (scf_mode==LINEAR_PEXSI) then
              if (iproc==0) call yaml_map('method','PEXSI')
              call pexsi_wrapper(iproc, nproc, bigdft_mpi%mpi_comm, &
                   tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
                   tmb%linmat%ovrlp_, tmb%linmat%ham_, &
                   foe_data_get_real(tmb%foe_obj,"charge",1), pexsi_npoles, pexsi_nproc_per_pole, &
                   pexsi_mumin, pexsi_mumax, &
                   pexsi_mu, pexsi_DeltaE, pexsi_temperature, pexsi_tol_charge, pexsi_np_sym_fact, &
                   pexsi_do_inertia_count, pexsi_max_iter, pexsi_verbosity, &
                   tmb%linmat%kernel_, energs%ebs)
          else if (scf_mode==LINEAR_FOE) then
              if (iproc==0) call yaml_map('method','FOE')
    
              call matrix_fermi_operator_expansion(iproc, nproc, bigdft_mpi%mpi_comm, &
                   tmb%foe_obj, tmb%ice_obj, tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
                   tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%ovrlppowers_(2), tmb%linmat%kernel_, &
                   energs%ebs, &
                   calculate_minusonehalf=invert_overlap_matrix, foe_verbosity=1, symmetrize_kernel=.true.)
          end if
    
          ! Eigenvalues not available, therefore take -.5d0
          tmb%orbs%eval=-.5d0
    
          !same as for directmin/diag
          ! CDFT: subtract V*w_ab from Hamiltonian so that we are calculating the correct energy
          if (present(cdft)) then
             call timing(iproc,'constraineddft','ON')
             tmparr = sparsematrix_malloc(tmb%linmat%smat(2),iaction=SPARSE_FULL,id='tmparr')
             call gather_matrix_from_taskgroups(iproc, nproc, bigdft_mpi%mpi_comm, &
                  tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr, tmparr)
             call daxpy(tmb%linmat%smat(2)%nvctr*tmb%linmat%smat(2)%nspin,-cdft%lag_mult,&
                  cdft%weight_matrix_%matrix_compr,1,tmparr,1)
             call extract_taskgroup(tmb%linmat%smat(2), tmparr, tmb%linmat%ham_%matrix_compr)
             call f_free(tmparr)
             call timing(iproc,'constraineddft','OF') 
    
             ! we have ebs and the kernel already, but only the CDFT function not actual energy for CDFT
             ! maybe pass both Hamiltonians to FOE to save calculation ebs twice?
             call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),tmb%linmat%smat(2), &
                  tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
                  tmb%coeff,orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup, .false.)
    
          end if
    
      end if
    
      call max_asymmetry_of_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(3), tmb%linmat%kernel_%matrix_compr, asymm_K)
    
      if (iproc==0) then
          call yaml_mapping_open('Asymmetry of the matrices')
          call yaml_map('Overlap',asymm_S,fmt='(es8.2)')
          call yaml_map('Hamiltonian',asymm_H,fmt='(es8.2)')
          call yaml_map('Kernel (possibly symmetrized)',asymm_K,fmt='(es8.2)')
          call yaml_mapping_close()
      end if
    
      if (calculate_ham) then
          if (calculate_KS_residue) then
              call get_KS_residue(iproc, nproc, tmb, orbs, hpsit_c, hpsit_f, KSres)
              if (iproc==0) call yaml_map('Kohn-Sham residue',KSres,fmt='(es10.3)')
          else
              if (iproc==0) call yaml_map('Kohn-Sham residue','not calculated')
          end if
          call f_free(hpsit_c)
          call f_free(hpsit_f)
      end if
      if (iproc==0) call yaml_map('Coefficients available',(scf_mode /= LINEAR_FOE .and. scf_mode /= LINEAR_PEXSI))
    
      if (calculate_ham) then
          if (nproc>1) then
              ! Wait for the communication of energs_work on root
              call mpi_fenceandfree(energs_work%window)
          end if
    
          ! Copy the value, only necessary on root
          if (iproc==0) then
              energs%ekin = energs_work%receivebuf(1)
              energs%epot = energs_work%receivebuf(2)
              energs%eproj = energs_work%receivebuf(3)
              energs%evsic = energs_work%receivebuf(4)
          end if
    
          !!call deallocate_work_mpiaccumulate(energs_work)
      end if
    
      if (iproc==0) call yaml_mapping_close() !close kernel update
    
    
      call f_release_routine()
    
    end subroutine get_coeff


    subroutine reconstruct_kernel(iproc, nproc, inversion_method, blocksize_dsyev, blocksize_pdgemm, &
               orbs, tmb, overlap_calculated)
      use module_base
      use module_types
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized
      use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
      use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, uncompress_matrix2
      use transposed_operations, only: calculate_overlap_transposed
      use coeffs, only: calculate_density_kernel
      implicit none
    
      ! Calling arguments
      integer,intent(in):: iproc, nproc, inversion_method, blocksize_dsyev, blocksize_pdgemm
      type(orbitals_data),intent(in):: orbs
      type(DFT_wavefunction),intent(inout):: tmb
      logical,intent(inout):: overlap_calculated
    
      ! Local variables
      !integer:: istat, iall
      character(len=*),parameter:: subname='reconstruct_kernel'
      integer :: i, j, ispin
    
      !call timing(iproc,'renormCoefComp','ON')
    
      ! Calculate the overlap matrix between the TMBs.
      if(.not. overlap_calculated) then
         if(.not.tmb%can_use_transposed) then
             !!if(associated(tmb%psit_c)) then
             !!    call f_free_ptr(tmb%psit_c)
             !!end if
             !!if(associated(tmb%psit_f)) then
             !!    call f_free_ptr(tmb%psit_f)
             !!end if
             !!tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
             !!tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
             call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                  TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
             tmb%can_use_transposed=.true.
         end if
         !call timing(iproc,'renormCoefComp','OF')
    
         call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
              tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, &
              tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%ovrlp_)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
         ! This can then be deleted if the transition to the new type has been completed.
         !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr
    
         !call timing(iproc,'renormCoefComp','ON')
         overlap_calculated=.true.
      end if
    
      tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(1), iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      !!do i=1,tmb%linmat%smat(1)%nspin*tmb%linmat%smat(1)%nvctr
      !!    write(2500,'(a,i8,es16.5)') 'i, tmb%linmat%ovrlp_%matrix_compr(i)', i, tmb%linmat%ovrlp_%matrix_compr(i)
      !!end do
      call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(1), tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
      !!do ispin=1,tmb%linmat%smat(1)%nspin
      !!    do i=1,tmb%linmat%smat(1)%nfvctr
      !!        do j=1,tmb%linmat%smat(1)%nfvctr
      !!            write(2600,'(a,3i8,es16.5)') 'ispin, i, j, tmb%linmat%ovrlp_%matrix(j,i,ispin)', ispin, i, j, tmb%linmat%ovrlp_%matrix(j,i,ispin)
      !!        end do
      !!    end do
      !!end do
      call reorthonormalize_coeff(iproc, nproc, orbs%norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, &
           tmb%orbs, tmb%linmat%smat(1), tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, orbs)
    
      call f_free_ptr(tmb%linmat%ovrlp_%matrix)
    
    
      ! Recalculate the kernel
      call calculate_density_kernel(iproc, nproc, bigdft_mpi%mpi_comm, .true., &
           orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup, &
           tmb%coeff, tmb%linmat%smat(3), tmb%linmat%kernel_)
      !call transform_sparse_matrix(tmb%linmat%denskern, tmb%linmat%denskern_large, 'large_to_small')
    
    end subroutine reconstruct_kernel
    
    !> Passing sparse ovrlp, but for now assuming ovrlp%matrix will be allocated and filled if using dense
    subroutine reorthonormalize_coeff(iproc, nproc, norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, basis_orbs, &
               basis_overlap, KS_overlap, basis_overlap_mat, coeff, orbs)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, &
           deallocate_matrices
      use yaml_output, only: yaml_newline, yaml_map
      use matrix_operations, only: overlapPowerGeneral, overlap_minus_one_half_serial, deviation_from_unity_parallel
      use orthonormalization, only: gramschmidt_coeff_trans
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: iproc, nproc, norb
      integer, intent(in) :: blocksize_dsyev, blocksize_pdgemm, inversion_method
      type(orbitals_data), intent(in) :: basis_orbs   !number of basis functions
      type(sparse_matrix),intent(in) :: basis_overlap
      type(sparse_matrix),dimension(basis_overlap%nspin),intent(in) :: KS_overlap
      type(matrices),intent(inout) :: basis_overlap_mat
      real(kind=8),dimension(basis_overlap%nfvctr,norb),intent(inout) :: coeff
      type(orbitals_data), intent(in) :: orbs   !Kohn-Sham orbitals that will be orthonormalized and their parallel distribution
      ! Local variables
      integer :: ierr, ind, iorb, korb, llorb, jorb, ist
      integer :: npts_per_proc, ind_start, ind_end, indc, ispin, norbx, iseg, i
      real(kind=8), dimension(:,:), allocatable :: coeff_tmp, coefftrans
      real(kind=8), dimension(:,:), allocatable :: ovrlp_coeff
      real(kind=8),dimension(:,:),pointer :: ovrlp_matrix, inv_ovrlp_matrix
      character(len=*),parameter:: subname='reorthonormalize_coeff'
      type(matrices) :: KS_ovrlp_
      type(matrices),dimension(1) :: inv_ovrlp_
      integer,dimension(2) :: irowcol
      integer,dimension(1) :: power
      !integer :: iorb, jorb !DEBUG
      real(kind=8) :: tt, max_error, mean_error!, tt2, tt3, ddot   !DEBUG
      !logical :: dense
      integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
      integer, parameter :: communication_strategy=ALLGATHERV
      logical,parameter :: dense=.true.
      logical,parameter :: check_accuracy=.false.
    
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr) ! to check timings
      call timing(iproc,'renormCoefCom1','ON')
    
      ! For a spin polarized calculation, the specified value of norb must be
      ! consistent with orbs%norb (not sure whether everything will work otherwise)
      if (basis_overlap%nspin==2 .and. norb/=orbs%norb) then
          stop 'ERROR: for spin polarized systems, norb must be the same as orbs%norb'
      end if
    
    
      !if (present(orbs)) then
      !   communication_strategy=ALLREDUCE
      !else
      !   communication_strategy=ALLGATHERV
      !end if
    
      spin_loop: do ispin=1,basis_overlap%nspin
    
          ! choose the correct number of KS orbitals. A bit ugly, maybe this routine
          ! should not be used with norb not conform with orbs%norb...
          if (orbs%norb/=norb) then
              norbx=norb
              ist=1
          else
              if (ispin==1) then
                  norbx=orbs%norbu
                  ist=1
              else
                  norbx=orbs%norbd
                  ist=orbs%norbu+1
              end if
          end if
    
          ovrlp_coeff=f_malloc((/norbx,norbx/), id='ovrlp_coeff')
    
          !!if(iproc==0) then
          !!    write(*,'(a)',advance='no') 'coeff renormalization...'
          !!end if
    
          !dense=.true.
    
          KS_ovrlp_ = matrices_null()
          ! can not use the wrapper since it cannot distinguish between up and down spin
          !call allocate_matrices(KS_overlap, allocate_full=.true., matname='KS_ovrlp_', mat=KS_ovrlp_)
          KS_ovrlp_%matrix = f_malloc_ptr((/norbx,norbx,1/))
    
          if (dense) then
             coeff_tmp=f_malloc((/basis_overlap%nfvctrp,max(norbx,1)/), id='coeff_tmp')
    
             ! Calculate the overlap matrix among the coefficients with respect to basis_overlap.
             if (basis_overlap%nfvctrp>0) then
                 !coeff_tmp=0.d0
                 !!do iorb=1,basis_overlap%nfvctr
                 !!    do jorb=1,basis_overlap%nfvctr
                 !!        write(2300+iproc,'(a,2i9,es13.5)') 'iorb, jorb, basis_overlap_mat%matrix(jorb,iorb,ispin)', iorb, jorb, basis_overlap_mat%matrix(jorb,iorb,ispin)
                 !!    end do
                 !!end do
                 !!do iorb=1,norb
                 !!    do jorb=1,basis_overlap%nfvctr
                 !!        write(2400+iproc,'(a,2i9,es13.5)') 'iorb, jorb, coeff(jorb,iorb)', iorb, jorb, coeff(jorb,iorb)
                 !!    end do
                 !!end do
                 !!call dgemm('n', 'n', basis_orbs%norbp, norb, basis_orbs%norb, 1.d0, basis_overlap_mat%matrix(basis_orbs%isorb+1,1,1), &
                 !!     basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
                 call dgemm('n', 'n', basis_overlap%nfvctrp, norbx, basis_overlap%nfvctr, &
                      1.d0, basis_overlap_mat%matrix(basis_overlap%isfvctr+1,1,ispin), &
                      basis_overlap%nfvctr, coeff(1,ist), basis_overlap%nfvctr, 0.d0, coeff_tmp, basis_overlap%nfvctrp)
                 !!do iorb=1,norbx
                 !!    do jorb=1,basis_overlap%nfvctrp
                 !!        write(2100+iproc,'(a,2i9,es13.5)') 'iorb, jorb, coeff_tmp(jorb,iorb)', iorb, jorb, coeff_tmp(jorb,iorb)
                 !!    end do
                 !!end do
                 !!call dgemm('t', 'n', norb, norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
                 !!     basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ovrlp_coeff, norb)
                 call dgemm('t', 'n', norbx, norbx, basis_overlap%nfvctrp, 1.d0, coeff(basis_overlap%isfvctr+1,ist), &
                      basis_overlap%nfvctr, coeff_tmp, basis_overlap%nfvctrp, 0.d0, ovrlp_coeff, norbx)
                 !!do iorb=1,norbx
                 !!    do jorb=1,norbx
                 !!        write(2200+iproc,'(a,i3,3x,2i9,es13.5)') 'ispin, iorb, jorb, ovrlp_coeff(jorb,iorb)', &
                 !!            ispin, iorb, jorb, ovrlp_coeff(jorb,iorb)
                 !!    end do
                 !!end do
              else
                 call f_zero(ovrlp_coeff)
              end if
    
              call f_free(coeff_tmp)
          else ! sparse - still less efficient than dense, also needs moving to a subroutine
    
             stop 'reorthonormalize_coeff: sparse version needs reworking'
             !also a problem with sparse at the moment - result not stored in correct arrays/allreduce etc
    
             !SM: need to fix the spin here
             call f_zero(KS_ovrlp_%matrix)
             npts_per_proc = nint(real(basis_overlap%nvctr + basis_overlap%nfvctr,dp) / real(nproc*2,dp))
             ind_start = 1+iproc*npts_per_proc
             ind_end = (iproc+1)*npts_per_proc
             if (iproc==nproc-1) ind_end = basis_overlap%nvctr!ceiling(0.5d0*real(basis_overlap%nvctr + basis_overlap%nfvctr,dp))
    
             indc=0
             do iseg=1,basis_overlap%nseg
                 ind=basis_overlap%keyv(iseg)
                 ! A segment is always on one line, therefore no double loop
                 do i = basis_overlap%keyg(1,1,iseg),basis_overlap%keyg(2,1,iseg)
                    !korb = basis_overlap%orb_from_index(1,ind)
                    !llorb = basis_overlap%orb_from_index(2,ind)
                    if (i<basis_overlap%keyg(1,2,iseg)) cycle ! so still only doing half
                    indc = indc + 1
                    if (indc < ind_start .or. indc > ind_end) cycle
    
                    do iorb=1,norb
                         if (basis_overlap%keyg(1,2,iseg)==i) then
                            tt=basis_overlap_mat%matrix_compr(ind)*coeff(i,iorb)
                            do jorb=iorb,norb
                                !SM: need to fix the spin here
                                KS_ovrlp_%matrix(jorb,iorb,1)=KS_ovrlp_%matrix(jorb,iorb,1) &
                                     +coeff(basis_overlap%keyg(1,2,iseg),jorb)*tt
                            end do
                         else
                            do jorb=iorb,norb
                                !SM: need to fix the spin here
                                KS_ovrlp_%matrix(jorb,iorb,1)=KS_ovrlp_%matrix(jorb,iorb,1) &
                                     +(coeff(basis_overlap%keyg(1,2,iseg),iorb)*coeff(i,jorb) &
                                     + coeff(basis_overlap%keyg(1,2,iseg),jorb)*coeff(i,iorb))&
                                     *basis_overlap_mat%matrix_compr(ind)
                            end do
                         end if
                     end do
                     ind=ind+1
                 end do
             end do
    
             ! use symmetry to calculate other half
             do iorb=1,norb
                do jorb=iorb+1,norb
                   KS_ovrlp_%matrix(iorb,jorb,1) = KS_ovrlp_%matrix(jorb,iorb,1)
                end do
             end do
    
          end if !sparse/dense
    
    
          if (nproc > 1) then
              call timing(iproc,'renormCoefCom1','OF')
              call timing(iproc,'renormCoefComm','ON')
              call mpiallred(ovrlp_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call timing(iproc,'renormCoefComm','OF')
              call timing(iproc,'renormCoefCom1','ON')
          end if
    
          !!if (iproc==0) call yaml_map('ovrlp_coeff',ovrlp_coeff)
          !!do iorb=1,norbx
          !!    do jorb=1,norbx
          !!        write(8000+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, ovrlp_coeff(jorb,iorb)',iorb, jorb, ovrlp_coeff(jorb,iorb)
          !!    end do
          !!end do
    
          ! Recalculate the coefficients
          call timing(iproc,'renormCoefCom1','OF')
    
          ! check whether this routine will be stable. Parallelization for nspin/=1 not done
          ! SM: I think one should rather pass KS_overlap instead of basis_overlap...
          if (norb==orbs%norb .and. basis_overlap%nspin==1) then
              if (orbs%norbp>0) then
                 call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
                      norbx, orbs%norbp, orbs%isorb, &
                      ovrlp_coeff(1:orbs%norb,orbs%isorb+1:orbs%isorb+orbs%norbp), &
                      basis_overlap, max_error, mean_error)
                 !!call deviation_from_unity_parallel(iproc, nproc, norbx, orbs%norbp, orbs%isorb, &
                 !!     ovrlp_coeff, &
                 !!     basis_overlap, max_error, mean_error)
              else
                 ! It is necessary to call the routine since it has a built-in mpiallred.
                 ! Use the first element of ovrlp_coeff; thanks to orbs%norbp==0 this should be safe
                 call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
                      orbs%norb, orbs%norbp, orbs%isorb, &
                      ovrlp_coeff(1:orbs%norb,1:orbs%norb), &
                      basis_overlap, max_error, mean_error)
              end if
          else
             ! SM: I think it's not good to call with iproc but 1 instead of nproc
             call deviation_from_unity_parallel(iproc, 1, bigdft_mpi%mpi_comm, &
                  norbx, norbx, 0, ovrlp_coeff(1,1), &
                  basis_overlap, max_error, mean_error)    
          end if
    
          if (iproc==0) then
              call yaml_newline()
              if (basis_overlap%nspin==1) then
                  call yaml_map('Max deviation from unity in reorthonormalize_coeff',max_error,fmt='(es8.2)')
                  call yaml_map('Mean deviation from unity in reorthonormalize_coeff',mean_error,fmt='(es8.2)')
              else
                  if (ispin==1) then
                      call yaml_map('spin up, Max deviation from unity in reorthonormalize_coeff',max_error,fmt='(es8.2)')
                      call yaml_map('spin up, Mean deviation from unity in reorthonormalize_coeff',mean_error,fmt='(es8.2)')
                  else if (ispin==2) then
                      call yaml_map('spin down, Max deviation from unity in reorthonormalize_coeff',max_error,fmt='(es8.2)')
                      call yaml_map('spin down, Mean deviation from unity in reorthonormalize_coeff',mean_error,fmt='(es8.2)')
                  end if
              end if
          end if
    
          if (max_error>5.0d0.and.orbs%norb==norb) then
             if (iproc==0) print*,'Error in reorthonormalize_coeff too large, reverting to gram-schmidt orthonormalization'
             ! gram-schmidt as too far from orthonormality to use iterative schemes for S^-1/2
             call f_free(ovrlp_coeff)
             call timing(iproc,'renormCoefCom2','ON')
             call gramschmidt_coeff_trans(iproc,nproc,orbs%norbu,orbs%norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
             call timing(iproc,'renormCoefCom2','OF')
          else
             ! standard lowdin
             ! Not clean to use twice basis_overlap, but it should not matter as everything
             ! is done using the dense version
    
             inv_ovrlp_(1) = matrices_null()
             ! can not use the wrapper since it cannot distinguish between up and down spin
             !call allocate_matrices(KS_overlap, allocate_full=.true., matname='inv_ovrlp_', mat=inv_ovrlp_)
             inv_ovrlp_(1)%matrix = f_malloc_ptr((/norbx,norbx,1/),id='inv_ovrlp_%matrix')
             
    
             if (norb==orbs%norb) then
                 !SM: need to fix the spin here
                 if (dense) call vcopy(norbx**2, ovrlp_coeff(1,1), 1, KS_ovrlp_%matrix(1,1,1), 1)
                 !!do iorb=1,norbx
                 !!    do jorb=1,norbx
                 !!        write(2000+iproc,'(a,2i9,es13.5)') 'iorb, jorb, KS_ovrlp_%matrix(jorb,iorb,1)', iorb, jorb, KS_ovrlp_%matrix(jorb,iorb,1)
                 !!    end do
                 !!end do
                 power(1)=-2
                 call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
                      inversion_method, 1, power, &
                      blocksize_dsyev, imode=2, ovrlp_smat=KS_overlap(ispin), inv_ovrlp_smat=KS_overlap(ispin), &
                      ovrlp_mat=KS_ovrlp_, inv_ovrlp_mat=inv_ovrlp_, &
                      check_accur=.false., nspinx=1)
                 !!do iorb=1,norbx
                 !!    do jorb=1,norbx
                 !!        write(8100+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, inv_ovrlp_%matrix(jorb,iorb,1)',iorb, jorb, inv_ovrlp_%matrix(jorb,iorb,1)
                 !!    end do
                 !!end do
    
             else
                 ! It is not possible to use the standard parallelization scheme, so do serial
                 ovrlp_matrix = f_malloc_ptr((/norbx,norbx/), id='ovrlp_matrix')
                 inv_ovrlp_matrix = f_malloc_ptr((/norbx,norbx/), id='inv_ovrlp_matrix')
                 call vcopy(norbx**2, ovrlp_coeff(1,1), 1, ovrlp_matrix(1,1), 1)
                 ! SM: I think it's not good to call with iproc but 1 instead of nproc
                 call overlap_minus_one_half_serial(iproc, 1, bigdft_mpi%mpi_comm, inversion_method, -2, blocksize_dsyev, &       
                      norbx, ovrlp_matrix, inv_ovrlp_matrix, check_accur=.false., smat=basis_overlap)
                 call f_free_ptr(ovrlp_matrix)
             !    call overlapPowerGeneral(iproc, 1, inversion_method, -2, &
             !         blocksize_dsyev, norb, orbs, imode=2, ovrlp_smat=basis_overlap, inv_ovrlp_smat=basis_overlap, &
             !         ovrlp_mat=basis_overlap_mat, inv_ovrlp_mat=inv_ovrlp, &
             !         check_accur=.false., ovrlp=ovrlp_coeff, inv_ovrlp=ovrlp_coeff2)
             end if
    
             call timing(iproc,'renormCoefCom2','ON')
    
             call f_free(ovrlp_coeff)
    
             ! Build the new linear combinations
             if (communication_strategy==ALLREDUCE) then
                if (basis_overlap%nspin/=1) then
                    stop 'reorthonormalize_coeff: for nspin/=1, the ALLREDUCE option is not implemented!'
                end if
                coeff_tmp=f_malloc0((/basis_overlap%nfvctr,norbx/), id='coeff_tmp')
    
                if (orbs%norbp>0) then
                    if (norb==orbs%norb) then
                        !SM: need to fix the spin here
                        call dgemm('n', 't', basis_orbs%norb, orbs%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), basis_orbs%norb, &
                             inv_ovrlp_(1)%matrix(1,orbs%isorb+1,1), orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
                        !@WARNING: THE FOLLOWING CALL IS NOT TESTED AND MIGHT BE WRONG!!
                        !!call dgemm('n', 't', basis_overlap%nfvctrp, norbx, norbx, 1.d0, coeff(basis_overlap%isfvctr,1), basis_overlap%nfvctr, &
                        !!     inv_ovrlp_%matrix(1,1,1), norbx, 0.d0, coeff_tmp(basis_overlap%isfvctr,1), basis_overlap%nfvctr)
                    else !surely this isn't correct??
                        call dgemm('n', 't', basis_orbs%norb, orbs%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), basis_orbs%norb, &
                             inv_ovrlp_matrix(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
                    end if
                else
                   call f_zero(coeff_tmp)
                end if
    
                if (nproc > 1) then
                   call mpiallred(coeff_tmp, mpi_sum, comm=bigdft_mpi%mpi_comm)
                end if
                call vcopy(basis_overlap%nfvctr*norbx,coeff_tmp(1,1),1,coeff(1,1),1)
             else
                coeff_tmp=f_malloc((/norbx,max(1,basis_overlap%nfvctrp)/), id='coeff_tmp')
                ! need to transpose so we can allgather - NOT VERY ELEGANT
                if (basis_orbs%norbp>0) then
                    if (norb==orbs%norb) then
                        !SM: need to fix the spin here
                        !!call dgemm('n', 't', norb, basis_orbs%norbp, norb, 1.d0, inv_ovrlp_%matrix(1,1,1), norb, &
                        !!    coeff(1+basis_orbs%isorb,1), basis_orbs%norb, 0.d0, coeff_tmp(1,1), norb)
                        call dgemm('n', 't', norbx, basis_overlap%nfvctrp, norbx, 1.d0, inv_ovrlp_(1)%matrix(1,1,1), norbx, &
                            coeff(1+basis_overlap%isfvctr,ist), basis_overlap%nfvctr, 0.d0, coeff_tmp(1,1), norbx)
                    else
                        call dgemm('n', 't', norb, basis_orbs%norbp, norb, 1.d0, inv_ovrlp_matrix(1,1), norb, &
                            coeff(1+basis_orbs%isorb,1), basis_orbs%norb, 0.d0, coeff_tmp(1,1), norb)
                    end if
                end if
                 !!do iorb=1,basis_overlap%nfvctrp
                 !!    do jorb=1,norbx
                 !!        write(8200+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, coeff_tmp(jorb,iorb)',iorb, jorb, coeff_tmp(jorb,iorb)
                 !!    end do
                 !!end do
    
                coefftrans=f_malloc((/norbx,basis_overlap%nfvctr/), id='coefftrans')
    
                ! gather together
                if(nproc > 1) then
                   call mpi_allgatherv(coeff_tmp(1,1), basis_overlap%nfvctrp*norbx, &
                        mpi_double_precision, coefftrans(1,1), &
                        norbx*basis_overlap%nfvctr_par(:), norbx*basis_overlap%isfvctr_par, &
                        mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                else
                   call vcopy(basis_overlap%nfvctrp*norbx,coeff_tmp(1,1),1,coefftrans(1,1),1)
                end if
                 !!do iorb=1,basis_overlap%nfvctr
                 !!    do jorb=1,norbx
                 !!        write(8300+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, coefftrans(jorb,iorb)',iorb, jorb, coefftrans(jorb,iorb)
                 !!    end do
                 !!end do
    
                ! untranspose coeff
                !$omp parallel do default(private) shared(coeff,coefftrans,norbx,basis_overlap,ist)
                do iorb=1,norbx
                   do jorb=1,basis_overlap%nfvctr
                      coeff(jorb,ist+iorb-1) = coefftrans(iorb,jorb)
                   end do
                end do
                !$omp end parallel do
    
                 !!do iorb=1,norbx
                 !!    do jorb=1,basis_overlap%nfvctr
                 !!        write(8400+10*iproc+ispin,'(a,3i8,es16.6)') 'iorb, jorb, ist, coeff(jorb,ist+iorb-1)',iorb, jorb, ist, coeff(jorb,ist+iorb-1)
                 !!    end do
                 !!end do
    
                call f_free(coefftrans)
             end if
    
             call timing(iproc,'renormCoefCom2','OF')
    
             call deallocate_matrices(inv_ovrlp_(1))
             if (norb/=orbs%norb) then
                 call f_free_ptr(inv_ovrlp_matrix)
             end if
    
             call f_free(coeff_tmp)
          end if
    
          if (check_accuracy) then
             ovrlp_coeff=f_malloc((/norbx,norbx/), id='ovrlp_coeff')
             coeff_tmp=f_malloc((/basis_overlap%nfvctrp,max(norbx,1)/), id='coeff_tmp')
             ! Calculate the overlap matrix among the coefficients with respect to basis_overlap.
             if (basis_overlap%nfvctrp>0) then
                coeff_tmp=0.d0
                !!call dgemm('n', 'n', basis_orbs%norbp, norb, basis_orbs%norb, 1.d0, basis_overlap_mat%matrix(basis_orbs%isorb+1,1,1), &
                !!     basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
                call dgemm('n', 'n', basis_overlap%nfvctrp, norbx, basis_overlap%nfvctr, &
                     1.d0, basis_overlap_mat%matrix(basis_overlap%isfvctr+1,1,1), &
                     basis_overlap%nfvctr, coeff(1,1), basis_overlap%nfvctr, 0.d0, coeff_tmp, basis_overlap%nfvctrp)
                !!call dgemm('t', 'n', norb, norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
                !!     basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ovrlp_coeff, norb)
                call dgemm('t', 'n', norbx, norbx, basis_overlap%nfvctrp, 1.d0, coeff(basis_overlap%isfvctr+1,1), &
                     basis_overlap%nfvctr, coeff_tmp, basis_overlap%nfvctrp, 0.d0, ovrlp_coeff, norbx)
             else
                call f_zero(ovrlp_coeff)
             end if
    
             call f_free(coeff_tmp)
    
             if (nproc>1) then
                call mpiallred(ovrlp_coeff, MPI_SUM, comm=bigdft_mpi%mpi_comm)
             end if
    
             if (norb==orbs%norb .and. basis_overlap%nspin==1) then
                ! Parallelization for nspin/=1 not done
                if (orbs%norbp>0) then
                   call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
                        orbs%norb, orbs%norbp, orbs%isorb, &
                        ovrlp_coeff(1:orbs%norb,orbs%isorb+1:orbs%isorb+orbs%norbp), &
                        basis_overlap, max_error, mean_error)
                else
                   ! It is necessary to call the routine since it has a built-in mpiallred.
                   ! Use the first element of ovrlp_coeff; thanks to orbs%norbp==0 this should be safe
                   call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
                        orbs%norb, orbs%norbp, &
                        orbs%isorb, ovrlp_coeff(1:orbs%norb,1:orbs%norb), &
                        basis_overlap, max_error, mean_error)
                end if
             else
                ! SM: I think it's not good to call with iproc but 1 instead of nproc
                call deviation_from_unity_parallel(iproc, 1, bigdft_mpi%mpi_comm, &
                     norbx, norbx, 0, &
                     ovrlp_coeff(1,1), basis_overlap, max_error, mean_error)
             end if
    
             if (iproc==0) print*,'Max deviation from unity following reorthonormalize_coeff',max_error
             if (iproc==0) print*,'Mean deviation from unity following reorthonormalize_coeff',mean_error
    
             !do iorb=1,norb
             !   do jorb=1,norb
             !      if (iproc==0) print*,jorb,iorb,ovrlp_coeff(jorb,iorb)
             !   end do
             !end do
    
             call f_free(ovrlp_coeff)
          end if
    
          call deallocate_matrices(KS_ovrlp_)
    
      end do spin_loop
    
      !!do iorb=1,norb
      !!    do jorb=1,basis_overlap%nfvctr
      !!        write(8500+10*iproc,'(a,2i8,es16.6)') 'iorb, jorb, coeff(jorb,iorb)',iorb, jorb, coeff(jorb,iorb)
      !!    end do
      !!end do
    
    
    end subroutine reorthonormalize_coeff


    subroutine renormalize_kernel(iproc, nproc, order_taylor, max_inversion_error, tmb, ovrlp, ovrlp_old)
      use module_base
      use module_types
      use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, assignment(=), &
                                   SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, &
                                   matrices, ONESIDED_FULL, ONESIDED_POST, ONESIDED_GATHER
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: uncompress_matrix
      use matrix_operations, only: overlapPowerGeneral, check_taylor_order
      use foe_common, only: retransform_ext
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(DFT_wavefunction),intent(inout) :: tmb
      type(matrices),intent(inout) :: ovrlp, ovrlp_old
    
      ! Local variables
      real(kind=8) :: max_error, mean_error, tt, tt1, tt2
      type(matrices) :: inv_ovrlp
      real(kind=8),dimension(:,:),pointer :: inv_ovrlpp, tempp
      real(kind=8),dimension(:),allocatable :: inv_ovrlp_compr_seq, kernel_compr_seq
      integer,dimension(3) :: power
      integer :: ispin, ilshift, ilshiftpp, isl
      integer,dimension(:,:),allocatable :: windowsx1
      real(kind=8),dimension(:),allocatable :: kernelpp_work1, kernelpp_work2
      real(kind=8),dimension(:),pointer :: matrix_local1
    
      call f_routine(id='renormalize_kernel')
    
      windowsx1 = f_malloc((/tmb%linmat%smat(3)%ntaskgroup,tmb%linmat%smat(3)%nspin/),id='windowsx1')
      !windowsx2 = f_malloc(tmb%linmat%smat(3)%ntaskgroup,id='windowsx2')
      kernelpp_work1 = f_malloc(tmb%linmat%smat(3)%smmm%nvctrp*tmb%linmat%smat(3)%nspin,id='kernelpp_work1')
      !kernelpp_work2 = f_malloc(tmb%linmat%smat(3)%smmm%nvctrp*tmb%linmat%smat(3)%nspin,id='kernelpp_work2')
      matrix_local1 = f_malloc_ptr(tmb%linmat%smat(3)%smmm%nvctrp_mm*tmb%linmat%smat(3)%nspin,id='matrix_local1')
    
      ! Calculate S^1/2 * K * S^1/2. Take the value of S^1/2 from memory (was
      ! calculated in the last call to this routine or (it it is the first call)
      ! just before the call.
      do ispin=1,tmb%linmat%smat(3)%nspin
          ilshift = (ispin-1)*tmb%linmat%smat(3)%nvctrp_tg
          ilshiftpp = (ispin-1)*tmb%linmat%smat(3)%smmm%nvctrp
          isl = (ispin-1)*tmb%linmat%smat(3)%smmm%nvctrp_mm
          call retransform_ext(iproc, nproc, tmb%linmat%smat(3), ONESIDED_POST, kernelpp_work1(ilshiftpp+1:), &
               tmb%linmat%ovrlppowers_(1)%matrix_compr(ilshift+1:), tmb%linmat%kernel_%matrix_compr(ilshift+1:), &
               matrix_localx=matrix_local1(isl+1:isl+tmb%linmat%smat(3)%smmm%nvctrp_mm), windowsx=windowsx1(:,ispin))
      end do
    
      ! Calculate S^1/2 for the overlap matrix
      power=(/2,-2,1/)
      call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
           order_taylor, 3, power, -1, &
           imode=1, ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
           ovrlp_mat=ovrlp, inv_ovrlp_mat=tmb%linmat%ovrlppowers_, &
           verbosity=0, &
           check_accur=order_taylor<1000, max_error=max_error, mean_error=mean_error, &
           ice_obj=tmb%ice_obj)
      call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
      
      !!if (iproc==0) write(*,*) 'in sub: sum(s-1)',sum(tmb%linmat%ovrlppowers_(1)%matrix_compr)
    
      do ispin=1,tmb%linmat%smat(3)%nspin
          ilshift = (ispin-1)*tmb%linmat%smat(3)%nvctrp_tg
          ilshiftpp = (ispin-1)*tmb%linmat%smat(3)%smmm%nvctrp
          isl = (ispin-1)*tmb%linmat%smat(3)%smmm%nvctrp_mm
          call retransform_ext(iproc, nproc, tmb%linmat%smat(3), ONESIDED_GATHER, kernelpp_work1(ilshiftpp+1:), &
               tmb%linmat%ovrlppowers_(1)%matrix_compr(ilshift+1:), tmb%linmat%kernel_%matrix_compr(ilshift+1:), &
               matrix_localx=matrix_local1(isl+1:isl+tmb%linmat%smat(3)%smmm%nvctrp_mm), windowsx=windowsx1(:,ispin))
      end do
    
      ! Calculate S^-1/2 * K * S^-1/2
      do ispin=1,tmb%linmat%smat(3)%nspin
          ilshift = (ispin-1)*tmb%linmat%smat(3)%nvctrp_tg
          ilshiftpp = (ispin-1)*tmb%linmat%smat(3)%smmm%nvctrp
          call retransform_ext(iproc, nproc, tmb%linmat%smat(3), ONESIDED_FULL, kernelpp_work1(ilshiftpp+1:), &
               tmb%linmat%ovrlppowers_(2)%matrix_compr(ilshift+1:), tmb%linmat%kernel_%matrix_compr(ilshift+1:), &
               windowsx=windowsx1(:,ispin))
      end do
    
      call f_free(windowsx1)
      !call f_free(windowsx2)
      call f_free(kernelpp_work1)
      !call f_free(kernelpp_work2)
      call f_free_ptr(matrix_local1)
    
      call f_release_routine()
    
    end subroutine renormalize_kernel



    subroutine calculate_gap_FOE(iproc, nproc, input, orbs_KS, tmb)
      use module_base
      use module_types    
      use foe_base, only: foe_data, foe_data_null, foe_data_get_real, foe_data_set_real, foe_data_deallocate
      !use foe, only:  fermi_operator_expansion_new
      use sparsematrix_highlevel, only: matrix_fermi_operator_expansion
      use sparsematrix_base, only: matrices_null, sparsematrix_malloc_ptr, deallocate_matrices, &
                                   SPARSE_TASKGROUP, assignment(=)
      use yaml_output
      use public_enums
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(input_variables),intent(in) :: input
      type(orbitals_data),intent(in) :: orbs_KS
      type(DFT_wavefunction),intent(inout) :: tmb
      
      ! Local variables
      type(foe_data) :: foe_obj
      real(kind=8) :: dq, ebs, tt, qq
      real(kind=8) :: e_homo, e_lumo
      integer :: ispin, norder_taylor, ind, iorb, norb, ntmb, iqq, iispin, jspin
      type(matrices),dimension(2) :: kernel
      logical :: calculation_possible, calculation_possible_all
      integer,dimension(input%nspin) :: ntmb_spin
      real(kind=8),dimension(input%nspin) :: e_homo_spin, e_lumo_spin
      logical,dimension(input%nspin) :: calculate_spin_channels
    
      if (iproc==0) call yaml_comment('FOE calculation for HOMO-LUMO analysis',hfill='=')
    
      ! Check whether the gap calculation is possible.
      ! This is the case if there are more support functions than occupied orbitals.
      !if (iproc==0) call yaml_mapping_open('Check possibility to calculate the gap')
      calculation_possible = .true.
      calculation_possible_all = .true.
      !do ispin=1,input%nspin
    
          if (iproc==0) call yaml_mapping_open('Check possibility to calculate the gap')
    
          !!if (ispin==1) then
          !!    norb = orbs_KS%norbu
          !!    ntmb = tmb%orbs%norbu
          !!else if (ispin==2) then
          !!    norb = orbs_KS%norbd
          !!    ntmb = tmb%orbs%norbd
          !!end if
          !!do ispin=1,input%nspin
              !!qq_spin(ispin) = 0.d0
              qq = 0.d0
              !if (ispin==1) then
                  !!ntmb_spin(ispin) = tmb%orbs%norbu
                  !ntmb = tmb%orbs%norbu
                  ntmb = tmb%orbs%norb
                  !do iorb=1,orbs_KS%norbu
                  do iorb=1,orbs_KS%norb
                      !qq_spin(ispin) = qq_spin(ispin) + orbs_KS%occup(iorb)
                      qq = qq + orbs_KS%occup(iorb)
                  end do
              !else if (ispin==2) then
              !    !ntmb_spin(ispin) = tmb%orbs%norbd
              !    ntmb = tmb%orbs%norbd
              !    do iorb=orbs_KS%norbu+1,orbs_KS%norbu+orbs_KS%norbd
              !        !qq_spin(ispin) = qq_spin(ispin) + orbs_KS%occup(iorb)
              !        qq = qq + orbs_KS%occup(iorb)
              !    end do
              !end if
          !!end do
          !!iispin = maxloc(qq_spin,1)
          !!ntmb = ntmb_spin(iispin)
          !!qq = qq_spin(iispin)
          if (input%nspin==1) then
              if (2*ntmb<=qq) then
                  calculation_possible = .false.
              end if
          else if (input%nspin==2) then
              if (ntmb<=qq) then
                  calculation_possible = .false.
              end if
          end if
          if (iproc==0) then
              !!call yaml_mapping_open('Checking individual spin component')
              !!call yaml_map('ispin',ispin)
              !call yaml_map('spin charges',qq_spin)
              !call yaml_map('spin channel',ispin)
              call yaml_map('charge for calculation',qq)
              call yaml_map('ntmb',ntmb)
              call yaml_map('Calculation possible',calculation_possible)
              !!call yaml_mapping_close()
          end if
          !end do
          if (iproc==0) then
              !!call yaml_map('Calculation possible',all(calculation_possible))
              call yaml_mapping_close()
          end if
                              
          if (calculation_possible) then
    
              !do jspin=1,input%nspin
              !    calculate_spin_channels(jspin) = (jspin==ispin)
              !end do
    
              ! To determine the HOMO/LUMO, subtract/add one electrom for closed shell
              ! systems of one half electron for open shell systems.
              !if (input%nspin==1) then
                  dq = 1.d0
              !else if (input%nspin==2) then
              !    dq = 0.5d0
              !end if
    
              ! determine the HOMO
              kernel(1) = matrices_null()
              kernel(1)%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(3), &
                  iaction=SPARSE_TASKGROUP, id='kernel%matrix_compr')
              foe_obj = foe_data_null()
              call init_foe_wrapper(iproc, nproc, input, orbs_KS, 0.d0, foe_obj)
    
              ! Round up the target charge (required for systems with non-integer charge)
              call foe_data_set_real(foe_obj,"charge",qq,1)
              !!qq = foe_data_get_real(foe_obj,"charge",1)
              iqq = nint(qq)
              if (real(iqq,kind=8)<qq) then
                  qq = real(iqq+1,kind=8)
              else
                  qq = real(iqq,kind=8)
              end if
              call foe_data_set_real(foe_obj,"charge",qq,1)
    
              call foe_data_set_real(foe_obj,"charge",foe_data_get_real(foe_obj,"charge",1)-dq,1)
              call foe_data_set_real(foe_obj,"fscale",1.d-2)
              norder_taylor = input%lin%order_taylor
              !call fermi_operator_expansion_new(iproc, nproc, &
              !     ebs, &
              !     .true., 2, &
              !     tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
              !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(1), foe_obj)
              if (iproc==0) then
                  call yaml_mapping_open('calculate HOMO kernel')
                  call yaml_map('target charge',foe_data_get_real(foe_obj,"charge",1))
              end if
              call matrix_fermi_operator_expansion(iproc, nproc, bigdft_mpi%mpi_comm, &
                   foe_obj, tmb%ice_obj, tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
                   tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%ovrlppowers_(2), kernel(1), &
                   ebs, calculate_minusonehalf=.true., foe_verbosity=1, symmetrize_kernel=.true.)!, &
                   !calculate_spin_channels=calculate_spin_channels)
              !call fermi_operator_expansion(iproc, nproc, &
              !     ebs, norder_taylor, input%lin%max_inversion_error, &
              !     .true., 2, &
              !     'HOMO', tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
              !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(1), foe_obj)
              !do ispin=1,input%nspin
                  !e_homo(ispin) = foe_data_get_real(foe_obj,"ef",ispin)
                  e_homo = foe_data_get_real(foe_obj,"ef",1)
              !end do
              call foe_data_deallocate(foe_obj)
              if (iproc==0) call yaml_mapping_close()
    
              ! determine the LUMO
              kernel(2) = matrices_null()
              kernel(2)%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(3), &
                  iaction=SPARSE_TASKGROUP, id='kernel%matrix_compr')
              foe_obj = foe_data_null()
              call init_foe_wrapper(iproc, nproc, input, orbs_KS, 0.d0, foe_obj)
    
              ! Round up the target charge (required for systems with non-integer charge)
              call foe_data_set_real(foe_obj,"charge",qq,1)
              !!qq = foe_data_get_real(foe_obj,"charge",1)
              iqq = nint(qq)
              if (real(iqq,kind=8)<qq) then
                  qq = real(iqq+1,kind=8)
              else
                  qq = real(iqq,kind=8)
              end if
              call foe_data_set_real(foe_obj,"charge",qq,1)
    
              call foe_data_set_real(foe_obj,"charge",foe_data_get_real(foe_obj,"charge",1)+dq,1)
              call foe_data_set_real(foe_obj,"fscale",1.d-2)
              norder_taylor = input%lin%order_taylor
              !call fermi_operator_expansion_new(iproc, nproc, &
              !     ebs, &
              !     .true., 2, &
              !     tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
              !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(2), foe_obj)
              if (iproc==0) then
                  call yaml_mapping_open('calculate LUMO kernel')
                  call yaml_map('target charge',foe_data_get_real(foe_obj,"charge",1))
              end if
              call matrix_fermi_operator_expansion(iproc, nproc, bigdft_mpi%mpi_comm, &
                   foe_obj, tmb%ice_obj, tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
                   tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%ovrlppowers_(2), kernel(1), &
                   ebs, calculate_minusonehalf=.true., foe_verbosity=1, symmetrize_kernel=.true.)!, &
                   !calculate_spin_channels=calculate_spin_channels)
              !call fermi_operator_expansion(iproc, nproc, &
              !     ebs, norder_taylor, input%lin%max_inversion_error, &
              !     .true., 2, &
              !     'LUMO', tmb%linmat%smat(1), tmb%linmat%smat(2), tmb%linmat%smat(3), &
              !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(2), foe_obj)
              !do ispin=1,input%nspin
                  !e_lumo(ispin) = foe_data_get_real(foe_obj,"ef",ispin)
                  e_lumo = foe_data_get_real(foe_obj,"ef",1)
              !end do
              call foe_data_deallocate(foe_obj)
              if (iproc==0) call yaml_mapping_close()
    
              !!e_homo_spin(ispin) = e_homo
              !!e_lumo_spin(ispin) = e_lumo
    
    
    
          !!if (iproc==0) then
          !!    tt = sqrt(kernel(2)%matrix_compr(1)-kernel(1)%matrix_compr(1))
          !!    do iorb=1,tmb%orbs%norb
          !!        write(*,*) 'iorb, val', iorb, (kernel(2)%matrix_compr(iorb)-kernel(1)%matrix_compr(iorb))/tt
          !!    end do
          !!end if
    
          call deallocate_matrices(kernel(1))
          call deallocate_matrices(kernel(2))
    
          else
              calculation_possible_all = .false.
          end if
      !!end do
    
      if (calculation_possible_all) then
          if (iproc==0) then
              call yaml_mapping_open('HOMO-LUMO analysis')
              !do ispin=1,input%nspin
                  !if (ispin==1) then
                  !    call yaml_mapping_open('spin channel up')
                  !else
                  !    call yaml_mapping_open('spin channel down')
                  !end if
                  !call yaml_map('spin_channel',ispin)
                  !!call yaml_map('HOMO energy',e_homo_spin(ispin))
                  !!call yaml_map('LUMO energy',e_lumo_spin(ispin))
                  !!call yaml_map('HOMO-LUMO gap (Ha)',e_lumo_spin(ispin)-e_homo_spin(ispin))
                  !!call yaml_map('HOMO-LUMO gap (eV)',(e_lumo_spin(ispin)-e_homo_spin(ispin))*Ha_eV)
                  call yaml_map('HOMO energy',e_homo)
                  call yaml_map('LUMO energy',e_lumo)
                  call yaml_map('HOMO-LUMO gap (Ha)',e_lumo-e_homo)
                  call yaml_map('HOMO-LUMO gap (eV)',(e_lumo-e_homo)*Ha_eV)
                  !call yaml_mapping_close()
              !end do
              call yaml_mapping_close()
          end if
      end if
    
    end subroutine calculate_gap_FOE


    subroutine get_KS_residue(iproc, nproc, tmb, KSorbs, hpsit_c, hpsit_f, KSres)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                                   matrices_null, allocate_matrices, deallocate_matrices, &
                                   sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
      use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, uncompress_matrix2
      use transposed_operations, only: calculate_overlap_transposed
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(DFT_wavefunction) :: tmb
      type(orbitals_data),intent(in) :: KSorbs
      real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(in) :: hpsit_c
      real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(in) :: hpsit_f
      real(kind=8),intent(out) :: KSres
    
      ! Local variables
      integer :: iorb, iiorb, ii, ispin!, ierr,  jorb
      real(kind=8) :: norbtot, scale_factor
      type(matrices) :: gradmat 
      real(kind=8),dimension(:,:,:),allocatable ::KH, KHKH, Kgrad
      character(len=*),parameter :: subname='get_KS_residue'
    
      call f_routine(id='get_KS_residue')
     
      !call nullify_sparse_matrix(gradmat)
      gradmat=matrices_null()
      call allocate_matrices(tmb%linmat%smat(2), allocate_full=.true., &
           matname='gradmat', mat=gradmat)
    
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
           hpsit_c, hpsit_c, hpsit_f, hpsit_f, tmb%linmat%smat(2), tmb%linmat%auxm, gradmat)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), gradmat)
    
    
      !gradmat%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='gradmat%matrix')
      tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(2), iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
      tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(3), iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
      call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(2), gradmat%matrix_compr, gradmat%matrix)
      call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(2), tmb%linmat%ham_%matrix_compr, tmb%linmat%ham_%matrix)
      call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(3), tmb%linmat%kernel_%matrix_compr, tmb%linmat%kernel_%matrix)
      KH=f_malloc0((/tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nspin/),id='KH')
      KHKH=f_malloc0((/tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nspin/),id='KHKH')
    
      ! scale_factor takes into account the occupancies which are present in the kernel
      if (KSorbs%nspin==1) then
          ! closed shell, i.e. factor 2 is included in the kernel
          scale_factor=0.5d0
      else
          scale_factor=1.0d0
      end if
    
      !!KHKH=0.d0
      !!call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, 1.0d0, tmb%linmat%denskern%matrix, &
      !!     tmb%orbs%norb, tmb%linmat%ham%matrix, tmb%orbs%norb, 0.d0, KH, tmb%orbs%norb)
      !!call dgemm('n', 't', tmb%orbs%norbp, tmb%orbs%norbp, tmb%orbs%norb, scale_factor, KH, &
      !!     tmb%orbs%norb, KH, tmb%orbs%norb, 0.d0, KHKH, tmb%orbs%norb)
      !!call mpiallred(KHKH(1,1), tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.0d0, tmb%linmat%denskern%matrix, &
      !!     tmb%orbs%norb, gradmat%matrix, tmb%orbs%norb, 0.d0, Kgrad, tmb%orbs%norb)
    
      call timing(iproc,'ks_residue','ON')
      ! Parallelized version
      if (tmb%orbs%norbp>0) then
          !call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.0d0, tmb%linmat%kernel_%matrix, &
          !     tmb%orbs%norb, tmb%linmat%ham_%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, &
          !     0.d0, KH(1,tmb%orbs%isorb+1), tmb%orbs%norb)
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              if (tmb%orbs%spinsgn(iiorb)>0.d0) then
                  ispin=1
              else
                  ispin=2
              end if
              ii=mod(iiorb-1,tmb%linmat%smat(3)%nfvctr)+1
              call dgemm('n', 'n', tmb%linmat%smat(3)%nfvctr, 1, &
                  tmb%linmat%smat(3)%nfvctr, 1.0d0, tmb%linmat%kernel_%matrix(1,1,ispin), &
                   tmb%linmat%smat(3)%nfvctr, tmb%linmat%ham_%matrix(1,ii,ispin), tmb%linmat%smat(3)%nfvctr, &
                   0.d0, KH(1,ii,ispin), tmb%linmat%smat(3)%nfvctr)
          end do
      end if
    
      if (nproc > 1) then
          call mpiallred(KH, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      if (tmb%orbs%norbp>0) then
          !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, scale_factor, KH, &
          !!     tmb%orbs%norb, KH(1,tmb%orbs%isorb+1), tmb%orbs%norb, &
          !!     0.d0, KHKH(1,tmb%orbs%isorb+1), tmb%orbs%norb)
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              if (tmb%orbs%spinsgn(iiorb)>0.d0) then
                  ispin=1
              else
                  ispin=2
              end if
              ii=mod(iiorb-1,tmb%linmat%smat(3)%nfvctr)+1
              call dgemm('n', 'n', tmb%linmat%smat(3)%nfvctr, 1, tmb%linmat%smat(3)%nfvctr, scale_factor, KH(1,1,ispin), &
                   tmb%linmat%smat(3)%nfvctr, KH(1,ii,ispin), tmb%linmat%smat(3)%nfvctr, &
                   0.d0, KHKH(1,ii,ispin), tmb%linmat%smat(3)%nfvctr)
          end do
      end if
    
      if (nproc > 1) then
          call mpiallred(KHKH, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      call f_free(KH)
      Kgrad=f_malloc0((/tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nspin/),id='Kgrad')
      if (tmb%orbs%norbp>0) then
          !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.0d0, tmb%linmat%kernel_%matrix, &
          !!     tmb%orbs%norb, gradmat%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, &
          !!     0.d0, Kgrad(1,tmb%orbs%isorb+1), tmb%orbs%norb)
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              if (tmb%orbs%spinsgn(iiorb)>0.d0) then
                  ispin=1
              else
                  ispin=2
              end if
              ii=mod(iiorb-1,tmb%linmat%smat(3)%nfvctr)+1
              call dgemm('n', 'n', tmb%linmat%smat(3)%nfvctr, 1, tmb%linmat%smat(3)%nfvctr, &
                   1.0d0, tmb%linmat%kernel_%matrix(1,1,ispin), &
                   tmb%linmat%smat(3)%nfvctr, gradmat%matrix(1,ii,ispin), tmb%linmat%smat(3)%nfvctr, &
                   0.d0, Kgrad(1,ii,ispin), tmb%linmat%smat(3)%nfvctr)
          end do
      end if
    
      if (nproc > 1) then
          call mpiallred(Kgrad, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      !!if (iproc==0) then
      !!  do iorb=1,tmb%orbs%norb
      !!    do jorb=1,tmb%orbs%norb
      !!      write(200,*) iorb, jorb, KHKH(jorb,iorb), Kgrad(jorb,iorb)
      !!    end do
      !!  end do
      !!end if
    
      !!! Sequential version
      !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.0d0, tmb%linmat%denskern%matrix, &
      !!     tmb%orbs%norb, tmb%linmat%ham%matrix(1,1), tmb%orbs%norb, 0.d0, KH, tmb%orbs%norb)
      !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, scale_factor, KH, &
      !!     tmb%orbs%norb, KH(1,1), tmb%orbs%norb, 0.d0, KHKH, tmb%orbs%norb)
      !!call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.0d0, tmb%linmat%denskern%matrix, &
      !!     tmb%orbs%norb, gradmat%matrix(1,1), tmb%orbs%norb, 0.d0, Kgrad, tmb%orbs%norb)
      !!if (iproc==0) then
      !!  do iorb=1,tmb%orbs%norb
      !!    do jorb=1,tmb%orbs%norb
      !!      write(201,*) iorb, jorb, KHKH(jorb,iorb), Kgrad(jorb,iorb)
      !!    end do
      !!  end do
      !!end if
    
    
      norbtot=0.d0
      do iorb=1,KSorbs%norb
          norbtot=norbtot+KSorbs%occup(iorb)
      end do
    
      KSres=0.d0
      do ispin=1,tmb%linmat%smat(3)%nspin
          do iorb=1,tmb%linmat%smat(3)%nfvctr
              KSres=KSres+Kgrad(iorb,iorb,ispin)-KHKH(iorb,iorb,ispin)
          end do
      end do
      KSres=sqrt(KSres/norbtot)
      !!if (iproc==0) write(*,*) 'KSgrad',sqrt(KSgrad/norbtot)
      call timing(iproc,'ks_residue','OF')
    
      !call f_free_ptr(gradmat%matrix)
      call f_free_ptr(tmb%linmat%ham_%matrix)
      call f_free_ptr(tmb%linmat%kernel_%matrix)
      !call f_free_ptr(gradmat%matrix_compr)
      call deallocate_matrices(gradmat)
    
    
      call f_free(KHKH)
      call f_free(Kgrad)
    
      call f_release_routine()
    
    end subroutine get_KS_residue


    !>  Coefficients are defined for Ntmb KS orbitals so as to maximize the number
    !!  of orthonormality constraints. This should speedup the convergence by
    !!  reducing the effective number of degrees of freedom.
    subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, fnrm_crit, itmax, energy, sd_fit_curve, &
        factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder, num_extra)
      use module_base
      use module_types
      !use module_interfaces, fake_name => optimize_coeffs
      use diis_sd_optimization
      use yaml_output
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use coeffs, only: calculate_kernel_and_energy
      implicit none
    
      ! Calling arguments
      integer,intent(in):: iproc, nproc, itmax, itout, it_scc, it_cdft
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(orbitals_data),intent(inout):: orbs
      type(DFT_wavefunction),intent(inout):: tmb
      type(DIIS_obj), intent(inout) :: ldiis_coeff
      real(kind=gp),intent(in):: fnrm_crit
      real(kind=gp),intent(out):: fnrm
      real(kind=gp), intent(inout) :: energy
      logical, intent(in) :: sd_fit_curve
      real(kind=gp), intent(in) :: factor
      integer, optional, intent(in) :: num_extra
      logical, optional, intent(in) :: reorder
    
      ! Local variables
      integer:: iorb, jorb, iiorb, ierr, it, itlast, ispin
      integer, dimension(1) :: iproc_arr, ncomp
      real(kind=gp),dimension(:,:),allocatable:: grad, grad_cov_or_coeffp !coeffp, grad_cov
      real(kind=gp),dimension(:),allocatable:: mat_coeff_diag, occup_tmp
      real(kind=gp) :: tt, ddot, energy0, pred_e
    
    
      !!if (present(num_extra) .and. tmb%linmat%smat(2)%nspin==2) then
      !!    stop 'ERROR: optimize_coeffs not yet implemented for nspin=2 and num_extra present!'
      !!end if
    
      call f_routine(id='optimize_coeffs')
    
    
      if (present(num_extra)) then
         if (num_extra > 0) then
            occup_tmp=f_malloc(orbs%norb,id='occup_tmp')
            ! make a copy of 'correct' occup
            call vcopy(orbs%norb, orbs%occup(1), 1, occup_tmp(1), 1)
         end if
      end if
    
    
      if (ldiis_coeff%idsx == 0 .and. sd_fit_curve) then
         ! calculate initial energy for SD line fitting and printing (maybe don't need to (re)calculate kernel here?)
         !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
         !!call extract_taskgroup_inplace(tmb%linmat%smat(2), tmb%linmat%ham_)
         call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),tmb%linmat%smat(2), &
              tmb%linmat%kernel_, tmb%linmat%ham_, energy0,&
              tmb%coeff, orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup,.true.)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
         !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
      else
         energy0=energy
      end if
    
      grad=f_malloc((/tmb%linmat%smat(2)%nfvctr,orbs%norbp/), id='grad')
      grad_cov_or_coeffp=f_malloc((/tmb%linmat%smat(2)%nfvctr,orbs%norbp/), id='grad_cov_or_coeffp')
    
      if (iproc==0) then
          call yaml_newline()
          call yaml_sequence_open('expansion coefficients optimization',label=&
               'it_coeff'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'_'//&
               trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))//&
               '_'//trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))))
      end if
    
    
      do it=1,itmax
    
          if (iproc==0) then
              call yaml_newline()
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_comment('it coeff:'//yaml_toa(it,fmt='(i6)'),hfill='-')
          end if
    
         !!do iorb=1,orbs%norbp
         !!    iiorb=orbs%isorb+iorb
         !!    if (orbs%spinsgn(iiorb)>0.d0) then
         !!        ispin=1
         !!    else
         !!        ispin=2
         !!    end if
         !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
         !!        write(5300+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
         !!    end do
         !!end do
    
         if (present(num_extra)) then
            if (num_extra > 0) then
               if (tmb%linmat%smat(2)%nspin==1) then
                   tt = 2.0_gp
               else if (tmb%linmat%smat(2)%nspin==2) then
                   tt = 1.0_gp
               end if
               do iorb=1,orbs%norb
                  orbs%occup(iorb)=tt
               end do
            end if
         end if
    
         call calculate_coeff_gradient(iproc,nproc,tmb,order_taylor,max_inversion_error,orbs,grad_cov_or_coeffp,grad)
    
         if (present(num_extra)) then
            if (num_extra > 0) then
               call vcopy(orbs%norb, occup_tmp(1), 1, orbs%occup(1), 1)
            end if
         end if
    
         ! scale the gradient by a factor of 2 in case of spin polarized calculation
         ! in order to make the gradient analogous to the case of no polarization
         ! (in that case the gradient is included in the kernel).
         if (tmb%linmat%smat(3)%nspin==2) then
             call dscal(orbs%norbp*tmb%linmat%smat(2)%nfvctr, 2.d0, grad, 1)
         end if
    
         !!do iorb=1,orbs%norbp
         !!    iiorb=orbs%isorb+iorb
         !!    if (orbs%spinsgn(iiorb)>0.d0) then
         !!        ispin=1
         !!    else
         !!        ispin=2
         !!    end if
         !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
         !!        write(5400+10*iproc+ispin,'(a,2i8,2es18.7)') 'iiorb, jorb, coeff(jorb,iiorb),grad(jorb,iorb)', iiorb, jorb, tmb%coeff(jorb,iiorb), grad(jorb,iorb)
         !!    end do
         !!end do
    
         ! Precondition the gradient (only making things worse...)
         !call precondition_gradient_coeff(tmb%orbs%norb, orbs%norbp, tmb%linmat%ham%matrix, tmb%linmat%ovrlp%matrix, grad)
    
         call timing(iproc,'dirmin_sddiis','ON')
    
         !For fnrm, we only sum on the occupied KS orbitals
         tt=0.d0
         do iorb=1,orbs%norbp
             tt=tt+ddot(tmb%linmat%smat(2)%nfvctr, grad_cov_or_coeffp(1,iorb), 1, grad(1,iorb), 1)
         end do
    
         if (nproc > 1) then
            call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
         end if
         fnrm=2.0_gp*tt
    
         !scale the gradient - useful for fragments/constrained
         call dscal(tmb%linmat%smat(2)%nfvctr*orbs%norbp,factor,grad,1)
    
         !write(*,*) 'present(num_extra), ldiis_coeff%idsx', present(num_extra), ldiis_coeff%idsx
         if (ldiis_coeff%idsx > 0) then !do DIIS
            !TO DO: make sure DIIS works
            ldiis_coeff%mids=mod(ldiis_coeff%ids,ldiis_coeff%idsx)+1
            ldiis_coeff%ids=ldiis_coeff%ids+1
    
            if (orbs%norbp>0) call vcopy(tmb%linmat%smat(2)%nfvctr*orbs%norbp,tmb%coeff(1,orbs%isorb+1), &
                1,grad_cov_or_coeffp(1,1),1)
    
            iproc_arr(1)=iproc
            ncomp(1)=tmb%linmat%smat(2)%nfvctr*orbs%norbp
            call diis_opt(iproc,nproc,1,0,1,iproc_arr,ncomp,tmb%linmat%smat(2)%nfvctr*orbs%norbp, &
                 grad_cov_or_coeffp,grad,ldiis_coeff) 
         else  !steepest descent with curve fitting for line minimization
            call timing(iproc,'dirmin_sddiis','OF')
            if (sd_fit_curve) call find_alpha_sd(iproc,nproc,ldiis_coeff%alpha_coeff,tmb,orbs,&
                 grad_cov_or_coeffp,grad,energy0,fnrm,pred_e,order_taylor)
            call timing(iproc,'dirmin_sddiis','ON')
            do iorb=1,orbs%norbp
               iiorb = orbs%isorb + iorb
               do jorb=1,tmb%linmat%smat(2)%nfvctr
                  grad_cov_or_coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*grad(jorb,iorb)
               end do
            end do
         end if
    
         !!do iorb=1,orbs%norbp
         !!    iiorb=orbs%isorb+iorb
         !!    if (orbs%spinsgn(iiorb)>0.d0) then
         !!        ispin=1
         !!    else
         !!        ispin=2
         !!    end if
         !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
         !!        write(4500+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, grad_cov_or_coeffp(jorb,iorb)', iiorb, jorb, grad_cov_or_coeffp(jorb,iorb)
         !!    end do
         !!end do
    
         call timing(iproc,'dirmin_sddiis','OF')
         
    
         call timing(iproc,'dirmin_allgat','ON')
         if(nproc > 1) then 
            call mpi_allgatherv(grad_cov_or_coeffp, tmb%linmat%smat(2)%nfvctr*orbs%norbp, mpi_double_precision, tmb%coeff, &
               tmb%linmat%smat(2)%nfvctr*orbs%norb_par(:,0), tmb%linmat%smat(2)%nfvctr*orbs%isorb_par, mpi_double_precision, &
               bigdft_mpi%mpi_comm, ierr)
         else
            call vcopy(tmb%linmat%smat(2)%nfvctr*orbs%norb,grad_cov_or_coeffp(1,1),1,tmb%coeff(1,1),1)
         end if
    
         call timing(iproc,'dirmin_allgat','OF')
    
         fnrm=sqrt(fnrm/dble(orbs%norb))
         ! Multiply the gradient by sqrt(2) to make it homologous to the case of
         ! no spin polarization (since with polarization norb is doubled).
         if (tmb%linmat%smat(3)%nspin==2) then
             fnrm=fnrm*sqrt(2.d0)
         end if
    
    
         !!do iorb=1,orbs%norbp
         !!    iiorb=orbs%isorb+iorb
         !!    if (orbs%spinsgn(iiorb)>0.d0) then
         !!        ispin=1
         !!    else
         !!        ispin=2
         !!    end if
         !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
         !!        write(4600+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
         !!    end do
         !!end do
    
         !! experimenting with calculating cHc and diagonalizing
         !if (present(num_extra).and.present(reorder)) then
         !   call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
         !   !call reorthonormalize_coeff(iproc, nproc, orbs%norb+num_extra, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
         !   call reordering_coeffs(iproc, nproc, num_extra, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
         !   if (reorder) call reordering_coeffs(iproc, nproc, num_extra, orbs, tmb%orbs, &
         !        tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
         !else if (present(reorder)) then
         !   call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
         !   !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
         !   call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
         !   if (reorder) call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, reorder)
         !else
         !   call reordering_coeffs(iproc, nproc, 0, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, tmb%coeff, .false.)
         !end if
    
         !!do iorb=1,orbs%norbp
         !!    iiorb=orbs%isorb+iorb
         !!    if (orbs%spinsgn(iiorb)>0.d0) then
         !!        ispin=1
         !!    else
         !!        ispin=2
         !!    end if
         !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
         !!        write(5200+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
         !!    end do
         !!end do
    
    
         ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
         ! instead of twice could add some criterion to check accuracy?
         ! first need to check if we're at least close to normalized, otherwise reorthonormalize_coeff may explode (should maybe incorporate this into reorthonormalize_coeff)
    
         ! should really check this first, for now just normalize anyway
         !mat_coeff_diag=f_malloc(orbs%norb,id='mat_coeff_diag')
         !call calculate_coeffMatcoeff_diag(tmb%linmat%ovrlp_%matrix,tmb%orbs,orbs,tmb%coeff,mat_coeff_diag)
         !do iorb=1,orbs%norb
         !   call dscal(tmb%orbs%norb,1.0d0/sqrt(mat_coeff_diag(iorb)),tmb%coeff(1,iorb),1)
         !end do
         !call f_free(mat_coeff_diag)
    
         call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, order_taylor, &
              tmb%orbs, tmb%linmat%smat(1), tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, orbs)
         !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
         !!!!!!!!!!!!!!!!!!!!!!!!
         !can't put coeffs directly in ksorbs%eval as intent in, so change after - problem with orthonormality of coeffs so adding extra
         !call find_eval_from_coeffs(iproc, nproc, tmb%orthpar%methTransformOverlap, orbs, tmb%orbs, tmb%linmat%ham, tmb%linmat%ovrlp, &
         !     tmb%coeff, tmb%orbs%eval, .true., .true.)
         !ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
         !call order_coeffs_by_energy(orbs%norb,tmb%orbs%norb,tmb%coeff,tmb%orbs%eval,ipiv)
         !call f_free(ipiv)
         !!!!!!!!!!!!!!!!!!!!!!!!
    
         !!do iorb=1,orbs%norbp
         !!    iiorb=orbs%isorb+iorb
         !!    if (orbs%spinsgn(iiorb)>0.d0) then
         !!        ispin=1
         !!    else
         !!        ispin=2
         !!    end if
         !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
         !!        write(5100+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
         !!    end do
         !!end do
    
         !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
         !!call extract_taskgroup_inplace(tmb%linmat%smat(2), tmb%linmat%ham_)
         call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),tmb%linmat%smat(2),&
              tmb%linmat%kernel_, tmb%linmat%ham_, energy,&
              tmb%coeff,orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup,.true.)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
    
         !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
         !write(127,*) ldiis_coeff%alpha_coeff,energy
         !close(127)
    
         ! can't check ebs only, need to check Etot in linearScaling, but if we do it>1 without curve fitting will still need to update here
         !if (ldiis_coeff%idsx == 0 .and. (.not. sd_fit_curve) .and. energy0/=0.0_gp) then ! only update alpha after first iteration
         !   ! apply a cap so that alpha_coeff never goes below around 1.d-2 or above 2
         !   if ((energy-energy0)<0.d0 .and. ldiis_coeff%alpha_coeff < 1.8d0) then
         !      ldiis_coeff%alpha_coeff=1.1d0*ldiis_coeff%alpha_coeff
         !   else if (ldiis_coeff%alpha_coeff > 1.7d-3) then
         !      ldiis_coeff%alpha_coeff=0.5d0*ldiis_coeff%alpha_coeff
         !   end if
         !   if (iproc==0) print*,'EBSdiff,alpha',energy-energy0,ldiis_coeff%alpha_coeff,energy,energy0
         !end if
    
         !if (iproc==0) write(*,*) ''
         if (sd_fit_curve .and. ldiis_coeff%idsx == 0) then
            !!if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha, pred E, diff',&
            !!     it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff,pred_e,pred_e-energy
            if (iproc==0) then
                call yaml_map('method','DminSD')
                call yaml_newline()
                call yaml_map('iter',it)
                call yaml_map('fnrm',fnrm,fmt='(es9.2)')
                call yaml_map('eBS',energy0,fmt='(es24.17)')
                call yaml_map('D',energy-energy0,fmt='(es10.3)')
                call yaml_map('alpha',ldiis_coeff%alpha_coeff,fmt='(es10.3)')
                call yaml_map('predicted energy',pred_e,fmt='(es24.17)')
                call yaml_map('D',pred_e-energy,fmt='(es10.3)')
            end if
         else if (ldiis_coeff%idsx == 0) then
            !!if (iproc==0) write(*,'(a,I4,2x,4(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha',&
            !!     it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff
            if (iproc==0) then
                call yaml_map('method','DminSD')
                call yaml_newline()
                call yaml_map('iter',it)
                call yaml_map('fnrm',fnrm,fmt='(es9.2)')
                call yaml_map('eBS',energy0,fmt='(es24.17)')
                call yaml_map('D',energy-energy0,fmt='(es10.3)')
                call yaml_map('alpha',ldiis_coeff%alpha_coeff,fmt='(es10.3)')
            end if
         else
            !!if (iproc==0) write(*,'(a,I4,2x,3(ES16.6e3,2x))')'DminDIIS: it, fnrm, ebs, ebsdiff',&
            !!     it,fnrm,energy0,energy-energy0
            if (iproc==0) then
                call yaml_map('method','DminDIIS')
                call yaml_newline()
                call yaml_map('iter',it)
                call yaml_map('fnrm',fnrm,fmt='(es9.2)')
                call yaml_map('eBS',energy0,fmt='(es24.17)')
                call yaml_map('D',energy-energy0,fmt='(es10.3)')
            end if
         end if
    
         energy0=energy
    
         if (iproc==0) then
             call yaml_mapping_close()
             call yaml_flush_document()
             !call bigdft_utils_flush(unit=6)
         end if
    
    
         itlast=it
         if (fnrm<fnrm_crit) exit
    
      end do
    
      !!if (iproc==0) then
      !!    call yaml_newline()
      !!    call yaml_sequence(label='final_coeff'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'_'//&
      !!         trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))//'_'//&
      !!         trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))),advance='no')
      !!    call yaml_mapping_open(flow=.true.)
      !!    call yaml_comment('iter:'//yaml_toa(itlast,fmt='(i6)'),hfill='-')
      !!    call yaml_map('iter',itlast,fmt='(i6)')
      !!    call yaml_map('fnrm',fnrm,fmt='(es9.2)')
      !!    call yaml_map('eBS',energy0,fmt='(es24.17)')
      !!    call yaml_map('D',energy-energy0,fmt='(es10.3)')
      !!    call yaml_mapping_close()
      !!    call bigdft_utils_flush(unit=6)
      !!end if
    
    
      if (iproc==0) then
          call yaml_sequence_close()
          call yaml_newline()
      end if
    
    
      call f_free(grad_cov_or_coeffp)
      call f_free(grad)
    
      if (present(num_extra)) then
         if (num_extra > 0) then
            call f_free(occup_tmp)
         end if
      end if
    
    
      call f_release_routine()
    
    end subroutine optimize_coeffs
    
    
    !> subset of reordering coeffs - need to arrange this routines better but taking the lazy route for now
    !! (also assuming we have no extra - or rather number of extra bands come from input.mix not input.lin)
    subroutine coeff_weight_analysis(iproc, nproc, input, ksorbs, tmb, ref_frags)
      use module_base
      use module_types
      use module_fragments
      use constrained_dft
      use yaml_output
      use sparsematrix_base, only: sparse_matrix, matrices, sparse_matrix_null, deallocate_sparse_matrix, &
                                   sparsematrix_malloc_ptr, DENSE_FULL, SPARSE_FULL, assignment(=), &
                                   matrices_null, allocate_matrices, deallocate_matrices,copy_sparse_matrix
      use sparsematrix, only: uncompress_matrix, uncompress_matrix2
      use matrix_operations, only: overlapPowerGeneral
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: iproc, nproc
      type(orbitals_data), intent(in) :: ksorbs
      type(dft_wavefunction), intent(inout) :: tmb
      type(input_variables),intent(in) :: input
      type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags
    
      integer :: iorb, istat, iall, ifrag
      integer, dimension(2) :: ifrag_charged
      integer, dimension(1) :: power
      !real(kind=8), dimension(:,:,:), allocatable :: weight_coeff
      real(kind=8), dimension(:,:), allocatable :: weight_coeff_diag
      real(kind=8), dimension(:,:), pointer :: ovrlp_half
      real(kind=8), dimension(:), allocatable :: weight_sum
      real(kind=8) :: max_error, mean_error, occsum
      type(sparse_matrix) :: weight_matrix
      type(matrices) :: weight_matrix_
      type(matrices),dimension(1) :: inv_ovrlp
      character(len=256) :: subname='coeff_weight_analysis'
    
      call timing(iproc,'weightanalysis','ON')
    
      call f_routine('coeff_weight_analysis')
    
      !call nullify_sparse_matrix(weight_matrix)
      weight_matrix=sparse_matrix_null()
      weight_matrix_=matrices_null()
      !call sparse_copy_pattern(tmb%linmat%smat(2), weight_matrix, iproc, subname)
      call copy_sparse_matrix(tmb%linmat%smat(2), weight_matrix)
      !weight_matrix_%matrix_compr=f_malloc_ptr(weight_matrix%nvctr,id='weight_matrix%matrix_compr')
      weight_matrix_%matrix_compr=sparsematrix_malloc_ptr(weight_matrix,iaction=SPARSE_FULL,id='weight_matrix%matrix_compr')
    
      inv_ovrlp(1) = matrices_null()
      call allocate_matrices(tmb%linmat%smat(3), allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))
    
      !weight_coeff=f_malloc((/ksorbs%norb,ksorbs%norb,input%frag%nfrag/), id='weight_coeff')
      weight_coeff_diag=f_malloc((/ksorbs%norb,input%frag%nfrag/), id='weight_coeff')
      !ovrlp_half=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='ovrlp_half')
      tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(1), DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      call uncompress_matrix2(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smat(1), tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
      power(1)=2
      call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
           input%lin%order_taylor, 1, power, &
           tmb%orthpar%blocksize_pdsyev, imode=2, &
           ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
           ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
           max_error=max_error, mean_error=mean_error)
      call f_free_ptr(tmb%linmat%ovrlp_%matrix)
    
      do ifrag=1,input%frag%nfrag
         ifrag_charged(1)=ifrag
         call calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,1,ifrag_charged,tmb,input%frag,ref_frags,&
              .false.,.true.,input%lin%order_taylor)
         weight_matrix_%matrix = sparsematrix_malloc_ptr(weight_matrix,iaction=DENSE_FULL,id='weight_matrix%matrix')
         call uncompress_matrix(iproc,nproc,weight_matrix,weight_matrix_%matrix_compr,weight_matrix_%matrix)
         !call calculate_coeffMatcoeff(nproc,weight_matrix%matrix,tmb%orbs,ksorbs,tmb%coeff,weight_coeff(1,1,ifrag))
         call calculate_coeffMatcoeff_diag(weight_matrix_%matrix,tmb%orbs,ksorbs,tmb%coeff,weight_coeff_diag(1,ifrag))
         call f_free_ptr(weight_matrix_%matrix)
      end do
      !call f_free_ptr(ovrlp_half)
    
      if (iproc==0) call yaml_sequence_open('Weight analysis',flow=.true.)
      if (iproc==0) call yaml_newline()
      if (iproc==0) call yaml_comment ('coeff, occ, eval, frac for each frag')
      ! only care about diagonal elements
      weight_sum=f_malloc(input%frag%nfrag,id='weight_sum')
      weight_sum=0.0d0
      do iorb=1,ksorbs%norb
         if (iproc==0) then
             call yaml_mapping_open(flow=.true.)
             call yaml_map('iorb',iorb,fmt='(i4)')
             call yaml_map('occ',KSorbs%occup(iorb),fmt='(f6.4)')
             call yaml_map('eval',tmb%orbs%eval(iorb),fmt='(f10.6)')
         end if
         do ifrag=1,input%frag%nfrag
            if (iproc==0) call yaml_map('frac',weight_coeff_diag(iorb,ifrag),fmt='(f6.4)')
            weight_sum(ifrag)=weight_sum(ifrag)+KSorbs%occup(iorb)*weight_coeff_diag(iorb,ifrag)
         end do
         if (iproc==0) call yaml_mapping_close()
         if (iproc==0) call yaml_newline()
      end do
         
      !sum for each fragment
      occsum=sum(KSorbs%occup(1:ksorbs%norb))
      if (iproc==0) then
         call yaml_mapping_open(flow=.true.)
         call yaml_map('sum',occsum,fmt='(f8.4)')
      end if
      do ifrag=1,input%frag%nfrag
         if (iproc==0) call yaml_map('frac',weight_sum(ifrag),fmt='(f6.4)')
      end do
      if (iproc==0) call yaml_mapping_close()
      if (iproc==0) call yaml_newline()
      call f_free(weight_sum)
    
      if (iproc==0) call yaml_sequence_close()
    
      call deallocate_matrices(inv_ovrlp(1))
    
      call deallocate_sparse_matrix(weight_matrix)
      call deallocate_matrices(weight_matrix_)
      call f_free(weight_coeff_diag)
      !call f_free(weight_coeff)
    
      call f_release_routine()
    
      call timing(iproc,'weightanalysis','OF')
    
    
    end subroutine coeff_weight_analysis
    
    
    !!! subset of reordering coeffs - need to arrange this routines better but taking the lazy route for now
    !!! (also assuming we have no extra - or rather number of extra bands come from input.mix not input.lin)
    !!subroutine find_eval_from_coeffs(iproc, nproc, meth_overlap, ksorbs, basis_orbs, ham, ovrlp, coeff, eval, calc_overlap, diag)
    !!  use module_base
    !!  use module_types
    !!  use module_interfaces
    !!  use sparsematrix_base, only: sparse_matrix
    !!  implicit none
    !!
    !!  ! Calling arguments
    !!  integer, intent(in) :: iproc, nproc, meth_overlap
    !!  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
    !!  type(sparse_matrix),intent(in) :: ham
    !!  type(sparse_matrix),intent(inout) :: ovrlp
    !!  real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
    !!  real(kind=8),dimension(ksorbs%norb),intent(inout) :: eval
    !!  logical, intent(in) :: diag, calc_overlap
    !!
    !!  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb
    !!  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ham_coeff,  ovrlp_coeff
    !!  real(kind=8) :: offdiagsum, offdiagsum2, coeff_orthog_threshold
    !!  character(len=256) :: subname='reordering_coeffs'
    !!
    !!  coeff_orthog_threshold=1.0d-3
    !!
    !!  allocate(ham_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
    !!  call memocc(istat, ham_coeff, 'ham_coeff', subname)
    !!
    !!  call calculate_coeffMatcoeff(ham%matrix,basis_orbs,ksorbs,coeff,ham_coeff)
    !!
    !!  if (calc_overlap) then
    !!     allocate(ovrlp_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
    !!     call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
    !!     call calculate_coeffMatcoeff(ovrlp%matrix,basis_orbs,ksorbs,coeff,ovrlp_coeff)
    !!  end if
    !!
    !!  ! above is overkill, actually just want diagonal elements but print off as a test out of curiosity
    !!  offdiagsum=0.0d0
    !!  offdiagsum2=0.0d0
    !!  do iorb=1,ksorbs%norb
    !!     do jorb=1,ksorbs%norb
    !!        if (iorb==jorb) then
    !!           eval(iorb)=ham_coeff(iorb,iorb)
    !!        else
    !!           offdiagsum=offdiagsum+abs(ham_coeff(iorb,jorb))
    !!           if (calc_overlap) offdiagsum2=offdiagsum2+abs(ovrlp_coeff(iorb,jorb))
    !!        end if
    !!     end do
    !!  end do
    !!  offdiagsum=offdiagsum/(ksorbs%norb**2-ksorbs%norb)
    !!  if (calc_overlap) offdiagsum2=offdiagsum2/(ksorbs%norb**2-ksorbs%norb)
    !!  if (calc_overlap.and.iproc==0) print*,''
    !!  if (calc_overlap) then
    !!     if (iproc==0) print*,'offdiagsum (ham,ovrlp):',offdiagsum,offdiagsum2
    !!  else
    !!     if (iproc==0) print*,'offdiagsum (ham):',offdiagsum
    !!  end if
    !!
    !!  ! if coeffs are too far from orthogonality
    !!  if (calc_overlap .and. offdiagsum2>coeff_orthog_threshold) then
    !!     call reorthonormalize_coeff(iproc, nproc, ksorbs%norb, -8, -8, meth_overlap, basis_orbs, ovrlp, coeff, ksorbs)
    !!  end if
    !!
    !!  if (diag.or.offdiagsum>1.0d-2) then
    !!     ! diagonalize within the space of occ+extra
    !!     if (.not.calc_overlap) then
    !!        ! assume ovrlp_coeff is orthgonal for now
    !!        allocate(ovrlp_coeff(ksorbs%norb,ksorbs%norb), stat=istat)
    !!        call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
    !!        call calculate_coeffMatcoeff(ovrlp%matrix,basis_orbs,ksorbs,coeff,ovrlp_coeff)
    !!        do iorb=1,ksorbs%norb
    !!           do jorb=1,ksorbs%norb
    !!              if (iorb==jorb) then
    !!                 ovrlp_coeff(iorb,jorb)=1.0d0
    !!              else
    !!                 ovrlp_coeff(iorb,jorb)=0.0d0
    !!              end if
    !!           end do
    !!        end do
    !!     end if
    !!     ! diagonalize within the space of occ+extra
    !!     call diagonalizeHamiltonian2(iproc, ksorbs%norb, ham_coeff, ovrlp_coeff, eval)
    !!     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norb/),id='coeff_tmp')
    !!
    !!     ! multiply new eigenvectors by coeffs
    !!     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norb, ksorbs%norb, 1.d0, coeff(1,1), &
    !!          basis_orbs%norb, ham_coeff, ksorbs%norb, 0.d0, coeff_tmp, basis_orbs%norb)
    !! 
    !!     call vcopy(basis_orbs%norb*(ksorbs%norb),coeff_tmp(1,1),1,coeff(1,1),1)
    !!     call f_free(coeff_tmp)
    !!  end if
    !!
    !!  if (calc_overlap.or.diag) then
    !!     iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
    !!     deallocate(ovrlp_coeff,stat=istat)
    !!     call memocc(istat,iall,'ovrlp_coeff',subname)
    !!  end if
    !!
    !!  iall=-product(shape(ham_coeff))*kind(ham_coeff)
    !!  deallocate(ham_coeff,stat=istat)
    !!  call memocc(istat,iall,'ham_coeff',subname)
    !!
    !!end subroutine find_eval_from_coeffs
    
    
    subroutine calculate_coeffMatcoeff(nproc,matrix,basis_orbs,ksorbs,coeff,mat_coeff)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: nproc
      type(orbitals_data), intent(in) :: basis_orbs, ksorbs
      real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(in) :: matrix
      real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
      real(kind=8),dimension(ksorbs%norb,ksorbs%norb),intent(inout) :: mat_coeff
    
      integer :: iall, istat, ierr
      real(kind=8), dimension(:,:), allocatable :: coeff_tmp
      character(len=256) :: subname='calculate_coeffMatcoeff'
    
      coeff_tmp = f_malloc((/ basis_orbs%norbp, ksorbs%norb /),id='coeff_tmp')
    
      if (basis_orbs%norbp>0) then
         call dgemm('n', 'n', basis_orbs%norbp, ksorbs%norb, basis_orbs%norb, 1.d0, matrix(basis_orbs%isorb+1,1), &
              basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
         call dgemm('t', 'n', ksorbs%norb, ksorbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
              basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, mat_coeff, ksorbs%norb)
      else
         call f_zero(mat_coeff)
      end if
    
      call f_free(coeff_tmp)
    
      if (nproc>1) then
          call mpiallred(mat_coeff, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
    end subroutine calculate_coeffMatcoeff
    
    
    !> Same as above but only calculating diagonal elements
    subroutine calculate_coeffMatcoeff_diag(matrix,basis_orbs,ksorbs,coeff,mat_coeff_diag)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(orbitals_data), intent(in) :: basis_orbs, ksorbs
      real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(in) :: matrix
      real(kind=8),dimension(basis_orbs%norb,ksorbs%norb),intent(inout) :: coeff
      real(kind=8),dimension(ksorbs%norb),intent(inout) :: mat_coeff_diag
    
      integer :: iall, istat, ierr, iorb
      real(kind=8), dimension(:,:), allocatable :: coeff_tmp
      real(kind=8), dimension(:), allocatable :: mat_coeff_diagp
      logical, parameter :: allgather=.true.
    
    
      if (allgather) then
         mat_coeff_diagp=f_malloc((/ksorbs%norbp/), id='mat_coeff_diagp')
      else
         call f_zero(mat_coeff_diag)
      end if
      if (ksorbs%norbp>0) then
         coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norbp/), id='coeff_tmp')
         call dgemm('n', 'n', basis_orbs%norb, ksorbs%norbp, basis_orbs%norb, 1.d0, matrix(1,1), &
              basis_orbs%norb, coeff(1,ksorbs%isorb+1), basis_orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
         if (allgather) then
            do iorb=1,ksorbs%norbp
               call dgemm('t', 'n', 1, 1, basis_orbs%norb, 1.d0, coeff(1,ksorbs%isorb+iorb), basis_orbs%norb, &
                    coeff_tmp(1,iorb), basis_orbs%norb, 0.d0, mat_coeff_diagp(iorb), ksorbs%norb)
            end do
         else
            do iorb=1,ksorbs%norbp
               call dgemm('t', 'n', 1, 1, basis_orbs%norb, 1.d0, coeff(1,ksorbs%isorb+iorb), basis_orbs%norb, &
                    coeff_tmp(1,iorb), basis_orbs%norb, 0.d0, mat_coeff_diag(ksorbs%isorb+iorb), ksorbs%norb)
            end do
         end if
         call f_free(coeff_tmp)
      end if
    
      if (bigdft_mpi%nproc>1) then
         if (allgather) then
            call mpi_allgatherv(mat_coeff_diagp, ksorbs%norbp, mpi_double_precision, mat_coeff_diag, &
                 ksorbs%norb_par(:,0), ksorbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         else
            call mpiallred(mat_coeff_diag, mpi_sum, comm=bigdft_mpi%mpi_comm)
         end if
      else
         if (allgather) then
            call vcopy(ksorbs%norb, mat_coeff_diagp(1), 1, mat_coeff_diag(1), 1)
         end if
      end if
    
      if (allgather) then
         call f_free(mat_coeff_diagp)
     end if
    
    end subroutine calculate_coeffMatcoeff_diag
     
    !> not really fragment related so prob should be moved - reorders coeffs by eval
      ! output ipiv in case want to use it for something else
      subroutine order_coeffs_by_energy(nstate,ntmb,coeff,eval,ipiv)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: nstate, ntmb
      real(kind=gp), dimension(ntmb,nstate), intent(inout) :: coeff
      real(kind=gp), dimension(nstate), intent(inout) :: eval
        integer, dimension(nstate), intent(out) :: ipiv
    
      integer :: itmb, jorb
        !integer, allocatable, dimension(:) :: ipiv
      real(gp), dimension(:), allocatable :: tmp_array
      real(gp), dimension(:,:), allocatable :: tmp_array2
    
        !ipiv=f_malloc(nstate,id='coeff_final')
      tmp_array=f_malloc(nstate,id='tmp_array')
    
      do itmb=1,nstate
         tmp_array(itmb)=-eval(itmb)
      end do
    
      call sort_positions(nstate,tmp_array,ipiv)
    
      do itmb=1,nstate
         eval(itmb)=-tmp_array(ipiv(itmb))
      end do
    
      call f_free(tmp_array)
    
      tmp_array2=f_malloc((/ntmb,nstate/),id='tmp_array2')
    
      do jorb=1,nstate
         do itmb=1,ntmb
            tmp_array2(itmb,jorb)=coeff(itmb,jorb)
         end do
      end do
    
      do jorb=1,nstate
         do itmb=1,ntmb
            coeff(itmb,jorb)=tmp_array2(itmb,ipiv(jorb))
         end do
      end do
    
      call f_free(tmp_array2)
        !call f_free(ipiv)
    
    end subroutine order_coeffs_by_energy
    
    !!! experimental - for num_extra see if extra states are actually lower in energy and reorder (either by sorting or diagonalization)
    !!! could improve efficiency by only calculating cSc or cHc up to norb (i.e. ksorbs%norb+extra)
    !!subroutine reordering_coeffs(iproc, nproc, num_extra, ksorbs, basis_orbs, ham, ovrlp, coeff, reorder)
    !!  use module_base
    !!  use module_types
    !!  use module_interfaces
    !!  use sparsematrix_base, only: sparse_matrix
    !!  implicit none
    !!
    !!  ! Calling arguments
    !!  integer, intent(in) :: iproc, nproc, num_extra
    !!  type(orbitals_data), intent(in) :: basis_orbs, ksorbs
    !!  type(sparse_matrix),intent(in) :: ham, ovrlp
    !!  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff
    !!  logical, intent(in) :: reorder
    !!
    !!  integer :: iorb, jorb, istat, iall, ierr, itmb, jtmb
    !!  integer, allocatable, dimension(:) :: ipiv
    !!  real(gp), dimension(:), allocatable :: tmp_array, eval
    !!  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ham_coeff, tmp_array2, ham_coeff_small, ovrlp_coeff_small, ovrlp_coeff
    !!  real(kind=8) :: offdiagsum, offdiagsum2
    !!  character(len=256) :: subname='reordering_coeffs'
    !!
    !!  allocate(ham_coeff(basis_orbs%norb,basis_orbs%norb), stat=istat)
    !!  call memocc(istat, ham_coeff, 'ham_coeff', subname)
    !!
    !!  allocate(coeff_tmp(basis_orbs%norbp,max(basis_orbs%norb,1)), stat=istat)
    !!  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
    !!
    !!  if (basis_orbs%norbp>0) then
    !!     call dgemm('n', 'n', basis_orbs%norbp, basis_orbs%norb, basis_orbs%norb, 1.d0, ham%matrix(basis_orbs%isorb+1,1), &
    !!          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
    !!     call dgemm('t', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
    !!          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ham_coeff, basis_orbs%norb)
    !!  else
    !!     call to_zero(basis_orbs%norb**2, ham_coeff(1,1))
    !!  end if
    !!
    !!  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
    !!  deallocate(coeff_tmp,stat=istat)
    !!  call memocc(istat,iall,'coeff_tmp',subname)
    !!
    !!  if (nproc>1) then
    !!      call mpiallred(ham_coeff(1,1), basis_orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    !!  end if
    !!
    !!  allocate(ovrlp_coeff(basis_orbs%norb,basis_orbs%norb), stat=istat)
    !!  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
    !!
    !!  allocate(coeff_tmp(basis_orbs%norbp,max(basis_orbs%norb,1)), stat=istat)
    !!  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
    !!
    !!  if (basis_orbs%norbp>0) then
    !!     call dgemm('n', 'n', basis_orbs%norbp, basis_orbs%norb, basis_orbs%norb, 1.d0, ovrlp%matrix(basis_orbs%isorb+1,1), &
    !!          basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
    !!     call dgemm('t', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
    !!          basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ovrlp_coeff, basis_orbs%norb)
    !!  else
    !!     call to_zero(basis_orbs%norb**2, ovrlp_coeff(1,1))
    !!  end if
    !!
    !!  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
    !!  deallocate(coeff_tmp,stat=istat)
    !!  call memocc(istat,iall,'coeff_tmp',subname)
    !!
    !!  if (nproc>1) then
    !!      call mpiallred(ovrlp_coeff(1,1), basis_orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    !!  end if
    !!
    !!  ! above is overkill, actually just want diagonal elements but print off as a test out of curiosity
    !!  offdiagsum=0.0d0
    !!  offdiagsum2=0.0d0
    !!  do iorb=1,ksorbs%norb+num_extra!basis_orbs%norb
    !!     do jorb=1,ksorbs%norb+num_extra!basis_orbs%norb
    !!        if (iorb==jorb) cycle
    !!        offdiagsum=offdiagsum+abs(ham_coeff(iorb,jorb))
    !!        offdiagsum2=offdiagsum2+abs(ovrlp_coeff(iorb,jorb))
    !!     end do
    !!  end do
    !!  offdiagsum=offdiagsum/(basis_orbs%norb**2-basis_orbs%norb)
    !!  offdiagsum2=offdiagsum2/(basis_orbs%norb**2-basis_orbs%norb)
    !!  if (iproc==0) print*,'offdiagsum (ham,ovrlp):',offdiagsum,offdiagsum2
    !!
    !!  ! sort the states - really need just ks+extra not all, otherwise sloshing!
    !!  ipiv=f_malloc(ksorbs%norb+num_extra,id='coeff_final')
    !!  tmp_array=f_malloc(ksorbs%norb+num_extra,id='tmp_array')
    !!  do itmb=1,ksorbs%norb+num_extra
    !!     tmp_array(itmb)=-ham_coeff(itmb,itmb)
    !!  end do
    !!
    !!  do itmb=1,ksorbs%norb+num_extra
    !!     ipiv(itmb)=itmb
    !!  end do
    !!  !call sort_positions(ksorbs%norb+num_extra,tmp_array,ipiv)
    !!  do jtmb=1,ksorbs%norb+num_extra
    !!     ham_coeff(jtmb,jtmb)=-tmp_array(ipiv(jtmb))
    !!  end do
    !!
    !!  tmp_array2=f_malloc((/basis_orbs%norb,ksorbs%norb+num_extra/),id='tmp_array2')
    !!  do jtmb=1,ksorbs%norb+num_extra
    !!     do itmb=1,basis_orbs%norb
    !!        tmp_array2(itmb,jtmb)=coeff(itmb,jtmb)
    !!     end do
    !!  end do
    !!
    !!  do jtmb=1,ksorbs%norb+num_extra
    !!     do itmb=1,basis_orbs%norb
    !!        coeff(itmb,jtmb)=tmp_array2(itmb,ipiv(jtmb))
    !!     end do
    !!  end do
    !!  call f_free(tmp_array2)
    !!
    !!  if (reorder) then
    !!     ! diagonalize within the space of occ+extra
    !!     eval=f_malloc((/ksorbs%norb+num_extra/),id='eval')
    !!     ham_coeff_small=f_malloc((/ksorbs%norb+num_extra,ksorbs%norb+num_extra/),id='ham_coeff_small')
    !!     ovrlp_coeff_small=f_malloc((/ksorbs%norb+num_extra,ksorbs%norb+num_extra/),id='ovrlp_coeff_small')
    !!     ! assume ovrlp_coeff is orthgonal for now
    !!     do iorb=1,ksorbs%norb+num_extra
    !!        do jorb=1,ksorbs%norb+num_extra
    !!           ham_coeff_small(iorb,jorb)=ham_coeff(iorb,jorb)
    !!           if (iorb==jorb) then
    !!              ovrlp_coeff_small(iorb,jorb)=1.0d0
    !!           else
    !!              ovrlp_coeff_small(iorb,jorb)=0.0d0
    !!           end if
    !!        end do
    !!     end do
    !!     ! diagonalize within the space of occ+extra
    !!     call diagonalizeHamiltonian2(iproc, ksorbs%norb+num_extra, ham_coeff_small, ovrlp_coeff_small, eval)
    !!     call f_free(ovrlp_coeff_small)
    !!     coeff_tmp=f_malloc((/basis_orbs%norb,ksorbs%norb+num_extra/),id='coeff_tmp')
    !!
    !!     ! multiply new eigenvectors by coeffs
    !!     call dgemm('n', 'n', basis_orbs%norb, ksorbs%norb+num_extra, ksorbs%norb+num_extra, 1.d0, coeff(1,1), &
    !!          basis_orbs%norb, ham_coeff_small, ksorbs%norb+num_extra, 0.d0, coeff_tmp, basis_orbs%norb)
    !! 
    !!     call f_free(ham_coeff_small)
    !!     call vcopy(basis_orbs%norb*(ksorbs%norb+num_extra),coeff_tmp(1,1),1,coeff(1,1),1)
    !!     call f_free(coeff_tmp)
    !!
    !!     do iorb=1,basis_orbs%norb
    !!        if (iorb<=ksorbs%norb) then
    !!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
    !!                ksorbs%occup(iorb),ham_coeff(iorb,iorb),eval(iorb)
    !!        else if (iorb<=ksorbs%norb+num_extra) then
    !!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
    !!                0.0d0,ham_coeff(iorb,iorb),eval(iorb)
    !!        !else
    !!        !   if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,ham_coeff(iorb,iorb),basis_orbs%occup(iorb),&
    !!        !        0.0d0,ham_coeff(iorb,iorb)
    !!        end if
    !!     end do
    !!     call f_free(eval)
    !!   else
    !!      do iorb=1,basis_orbs%norb
    !!        if (iorb<=ksorbs%norb) then
    !!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
    !!                ksorbs%occup(iorb),ham_coeff(iorb,iorb)
    !!        else if (iorb<=ksorbs%norb+num_extra) then
    !!           if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,-tmp_array(ipiv(iorb)),basis_orbs%occup(iorb),&
    !!                0.0d0,ham_coeff(iorb,iorb)
    !!        !else
    !!        !   if (iproc==0) write(*,*) 'optimize coeffs eval',iorb,ham_coeff(iorb,iorb),basis_orbs%occup(iorb),&
    !!        !        0.0d0,ham_coeff(iorb,iorb)
    !!        end if
    !!     end do
    !!  end if
    !!
    !!  call f_free(ipiv)
    !!  call f_free(tmp_array)
    !!  iall=-product(shape(ham_coeff))*kind(ham_coeff)
    !!  deallocate(ham_coeff,stat=istat)
    !!  call memocc(istat,iall,'ham_coeff',subname)
    !!
    !!end subroutine reordering_coeffs
    
    
    subroutine find_alpha_sd(iproc,nproc,alpha,tmb,orbs,coeffp,grad,energy0,fnrm,pred_e,taylor_order)
      use module_base
      use module_types
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use coeffs, only: calculate_kernel_and_energy
      implicit none
      integer, intent(in) :: iproc, nproc, taylor_order
      real(kind=gp), intent(inout) :: alpha
      type(DFT_wavefunction) :: tmb
      type(orbitals_data), intent(in) :: orbs
      real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(inout) :: coeffp
      real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(in) :: grad
      real(kind=gp), intent(in) :: energy0, fnrm
      real(kind=gp), intent(out) :: pred_e
      integer :: iorb, iiorb, jorb, ierr
      real(kind=gp) :: energy1, a, b, c, alpha_old
      real(kind=gp),dimension(:,:),allocatable :: coeff_tmp
    
      call timing(iproc,'dirmin_sdfit','ON')
    
      ! take an initial step to get 2nd point
      coeff_tmp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_tmp')
      do iorb=1,orbs%norbp
         iiorb = orbs%isorb + iorb
         do jorb=1,tmb%linmat%smat(2)%nfvctr
            coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-alpha*grad(jorb,iorb)
         end do
      end do
    
      if(nproc > 1) then 
         call mpi_allgatherv(coeffp, tmb%linmat%smat(2)%nfvctr*orbs%norbp, mpi_double_precision, coeff_tmp, &
              tmb%linmat%smat(2)%nfvctr*orbs%norb_par(:,0), tmb%linmat%smat(2)%nfvctr*orbs%isorb_par, &
              mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call vcopy(tmb%linmat%smat(2)%nfvctr*orbs%norb,coeffp(1,1),1,coeff_tmp(1,1),1)
      end if
    
      ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
      ! instead of twice could add some criterion to check accuracy?
      call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, taylor_order, tmb%orbs, &
           tmb%linmat%smat(1), tmb%linmat%ks, tmb%linmat%ovrlp_, coeff_tmp, orbs)
      !call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, &
      !     tmb%linmat%smat(1), tmb%linmat%ks, tmb%linmat%ovrlp_, coeff_tmp, orbs)
      !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
      !!call extract_taskgroup_inplace(tmb%linmat%smat(2), tmb%linmat%ham_)
      call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),tmb%linmat%smat(2),&
           tmb%linmat%kernel_, tmb%linmat%ham_, energy1,&
           coeff_tmp,orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup,.true.)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
      !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
      call f_free(coeff_tmp)
    
      ! find ideal alpha using both points
      alpha_old=alpha
      a=fnrm/alpha_old+(energy1-energy0)/alpha_old**2
      b=-fnrm
      c=energy0
      alpha=-0.5_gp*b/a
      ! don't update if we have found a maximum, or negative alpha is predicted
      ! do something better here - don't just want to do the same thing twice, so at least check if energy has decreased
      if (alpha<0.0_gp .or. a<0.0_gp) alpha=alpha_old
      pred_e=a*alpha**2+b*alpha+c
    
      !open(127)
      !write(127,*) '#',a,b,c,(energy1-energy0)/alpha_old,b-(energy1-energy0)/alpha_old
      !write(127,*) 0.0_gp,energy0
      !write(127,*) alpha_old,energy1
    
      call timing(iproc,'dirmin_sdfit','OF')
    
    end subroutine find_alpha_sd
    
    
    
    
    !> calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
    !! then grad=S^-1grad_cov
    subroutine calculate_coeff_gradient(iproc,nproc,tmb,order_taylor,max_inversion_error,KSorbs,grad_cov,grad)
      use module_base
      use module_types
      use sparsematrix_base, only: matrices, sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), &
                                   matrices_null, allocate_matrices, deallocate_matrices
      use sparsematrix, only: extract_taskgroup_inplace, gather_matrix_from_taskgroups_inplace
      use matrix_operations, only: overlapPowerGeneral, check_taylor_order
      use parallel_linalg, only: dgesv_parallel
      use coeffs, only: calculate_density_kernel
      implicit none
    
      integer, intent(in) :: iproc, nproc
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(DFT_wavefunction), intent(inout) :: tmb
      type(orbitals_data), intent(in) :: KSorbs
      real(gp), dimension(tmb%linmat%smat(2)%nfvctr,KSorbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp
    
      integer :: iorb, iiorb, info, ierr, ispin
      integer, dimension(1) :: power
      real(gp),dimension(:,:,:),allocatable :: sk, skh, skhp
      real(gp),dimension(:,:),pointer :: inv_ovrlp
      real(kind=gp), dimension(:,:), allocatable:: grad_full
      character(len=*),parameter:: subname='calculate_coeff_gradient'
      type(matrices),dimension(1) :: inv_ovrlp_
    
      integer :: itmp, itrials, jorb
      integer :: ncount1, ncount_rate, ncount_max, ncount2
      real(kind=4) :: tr0, tr1
      real(kind=8) :: deviation, time, max_error, mean_error, maxerror
    
    
      call f_routine(id='calculate_coeff_gradient')
      call timing(iproc,'dirmin_lagmat1','ON')
    
      inv_ovrlp_(1) = matrices_null()
      call allocate_matrices(tmb%linmat%smat(3), allocate_full=.true., matname='inv_ovrlp_', mat=inv_ovrlp_(1))
    
    
      ! we have the kernel already, but need it to not contain occupations so recalculate here
      ! don't want to lose information in the compress/uncompress process - ideally need to change sparsity pattern of kernel
      !call calculate_density_kernel(iproc, nproc, .false., KSorbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
      tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(3), iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
    
      !!do iorb=1,KSorbs%norbp
      !!    iiorb=KSorbs%isorb+iorb
      !!    if (KSorbs%spinsgn(iiorb)>0.d0) then
      !!        ispin=1
      !!    else
      !!        ispin=2
      !!    end if
      !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
      !!        write(4900+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
      !!    end do
      !!end do
    
      !!do iiorb=1,KSorbs%norb
      !!    if (KSorbs%spinsgn(iiorb)>0.d0) then
      !!        ispin=1
      !!    else
      !!        ispin=2
      !!    end if
      !!    do jorb=1,tmb%linmat%smat(2)%nfvctr
      !!        write(4800+10*iproc+ispin,'(a,2i8,es18.7)') 'iiorb, jorb, coeff(jorb,iiorb)', iiorb, jorb, tmb%coeff(jorb,iiorb)
      !!    end do
      !!end do
    
     !call uncompress_matrix(iproc,tmb%linmat%denskern)
      !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
      call calculate_density_kernel(iproc, nproc, bigdft_mpi%mpi_comm, .false., &
           KSorbs%norbp, KSorbs%isorb, KSorbs%norbu, KSorbs%norb, KSorbs%occup, &
           tmb%coeff, tmb%linmat%smat(3), tmb%linmat%kernel_, .true.)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
    
      sk=f_malloc0((/tmb%linmat%smat(2)%nfvctrp,tmb%linmat%smat(2)%nfvctr,tmb%linmat%smat(2)%nspin/), id='sk')
    
      ! calculate I-S*K - first set sk to identity
      do ispin=1,tmb%linmat%smat(2)%nspin
          do iorb=1,tmb%linmat%smat(2)%nfvctrp
              iiorb=mod(tmb%linmat%smat(2)%isfvctr+iorb-1,tmb%linmat%smat(2)%nfvctr)+1
              sk(iorb,iiorb,ispin) = 1.d0
          end do 
      end do
    
      if (tmb%linmat%smat(2)%nfvctrp>0) then
          do ispin=1,tmb%linmat%smat(2)%nspin
              call dgemm('t', 'n', tmb%linmat%smat(2)%nfvctrp, tmb%linmat%smat(2)%nfvctr, tmb%linmat%smat(2)%nfvctr, -1.d0, &
                   tmb%linmat%ovrlp_%matrix(1,tmb%linmat%smat(2)%isfvctr+1,ispin), tmb%linmat%smat(2)%nfvctr, &
                   tmb%linmat%kernel_%matrix(1,1,ispin), tmb%linmat%smat(2)%nfvctr, 1.d0, sk(1,1,ispin), tmb%linmat%smat(2)%nfvctrp)
          end do
      end if
      !!do ispin=1,tmb%linmat%smat(2)%nspin
      !!    do iorb=1,tmb%linmat%smat(2)%nfvctr
      !!        do jorb=1,tmb%linmat%smat(2)%nfvctr
      !!            write(6300+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, tmb%linmat%ovrlp_%matrix(jorb,iorb,ispin)', iorb, jorb, tmb%linmat%ovrlp_%matrix(jorb,iorb,ispin)
      !!            write(6400+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, tmb%linmat%kernel_%matrix(jorb,iorb,ispin)', iorb, jorb, tmb%linmat%kernel_%matrix(jorb,iorb,ispin)
      !!        end do
      !!    end do
      !!end do
      !!do ispin=1,tmb%linmat%smat(2)%nspin
      !!    do iorb=1,tmb%linmat%smat(2)%nfvctr
      !!        do jorb=1,tmb%linmat%smat(2)%nfvctrp
      !!            write(6100+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, sk(jorb,iorb,ispin)', iorb, jorb, sk(jorb,iorb,ispin)
      !!        end do
      !!    end do
      !!end do
    
      ! coeffs and therefore kernel will change, so no need to keep it
      call f_free_ptr(tmb%linmat%kernel_%matrix)
    
      skhp=f_malloc((/tmb%linmat%smat(2)%nfvctr,max(tmb%linmat%smat(2)%nfvctrp,1),tmb%linmat%smat(2)%nspin/), id='skhp')
    
      ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
      if (tmb%linmat%smat(2)%nfvctrp>0) then
          do ispin=1,tmb%linmat%smat(2)%nspin
              call dgemm('t', 't', tmb%linmat%smat(2)%nfvctr, tmb%linmat%smat(2)%nfvctrp, tmb%linmat%smat(2)%nfvctr, &
                   1.d0, tmb%linmat%ham_%matrix(1,1,ispin), tmb%linmat%smat(2)%nfvctr, &
                   sk(1,1,ispin), tmb%linmat%smat(2)%nfvctrp, 0.d0, skhp(1,1,ispin), tmb%linmat%smat(2)%nfvctr)
          end do
      end if
      !!do ispin=1,tmb%linmat%smat(2)%nspin
      !!    do iorb=1,tmb%linmat%smat(2)%nfvctrp
      !!        do jorb=1,tmb%linmat%smat(2)%nfvctr
      !!            write(6200+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, skhp(jorb,iorb,ispin)', iorb, jorb, skhp(jorb,iorb,ispin)
      !!        end do
      !!    end do
      !!end do
    
      call f_free(sk)
    
      skh=f_malloc((/tmb%linmat%smat(2)%nfvctr,tmb%linmat%smat(2)%nfvctr,tmb%linmat%smat(2)%nspin/), id='skh')
    
      call timing(iproc,'dirmin_lagmat1','OF')
      call timing(iproc,'dirmin_lagmat2','ON')
    
      ! gather together
      if(nproc > 1) then
          do ispin=1,tmb%linmat%smat(2)%nspin
             call mpi_allgatherv(skhp(1,1,ispin), tmb%linmat%smat(2)%nfvctr*tmb%linmat%smat(2)%nfvctrp, &
                mpi_double_precision, skh(1,1,ispin), tmb%linmat%smat(2)%nfvctr*tmb%linmat%smat(2)%nfvctr_par, &
                tmb%linmat%smat(2)%nfvctr*tmb%linmat%smat(2)%isfvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
          end do
      else
         call vcopy(tmb%linmat%smat(2)%nfvctr*tmb%linmat%smat(2)%nfvctrp*tmb%linmat%smat(2)%nspin,skhp(1,1,1),1,skh(1,1,1),1)
      end if
    
      !!do ispin=1,tmb%linmat%smat(2)%nspin
      !!    do iorb=1,tmb%linmat%smat(2)%nfvctr
      !!        do jorb=1,tmb%linmat%smat(2)%nfvctr
      !!            write(6000+10*iproc+ispin,'(a,2i8,es16.6)') 'iorb, jorb, skh(jorb,iorb,ispin)', iorb, jorb, skh(jorb,iorb,ispin)
      !!        end do
      !!    end do
      !!end do
    
      call timing(iproc,'dirmin_lagmat2','OF')
      call timing(iproc,'dirmin_lagmat1','ON')
    
      call f_free(skhp)
    
      ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
      if (KSorbs%norbp>0) then
         !call dgemm('t', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
         !     tmb%orbs%norb, tmb%coeff(1,KSorbs%isorb+1), tmb%orbs%norb, 0.d0, grad_cov(1,1), tmb%orbs%norb)
         do iorb=1,KSorbs%norbp
             iiorb=KSorbs%isorb+iorb
             if (KSorbs%spinsgn(iiorb)>0.d0) then
                 ispin=1
             else
                 ispin=2
             end if
             call dgemm('t', 'n', tmb%linmat%smat(2)%nfvctr, 1, tmb%linmat%smat(2)%nfvctr, 1.d0, skh(1,1,ispin), &
                  tmb%linmat%smat(2)%nfvctr, tmb%coeff(1,KSorbs%isorb+iorb), tmb%linmat%smat(2)%nfvctr, &
                  0.d0, grad_cov(1,iorb), tmb%linmat%smat(2)%nfvctr)
         end do
      end if
    
      call f_free(skh)
    
      ! multiply by f_i to get grad_i^a
      do iorb=1,KSorbs%norbp
         iiorb=KSorbs%isorb+iorb
         grad_cov(:,iorb)=grad_cov(:,iorb)*KSorbs%occup(iiorb)
      end do
    
      call timing(iproc,'dirmin_lagmat1','OF')
      call timing(iproc,'dirmin_dgesv','ON')
    
      ! Solve the linear system ovrlp*grad=grad_cov
      if(tmb%orthpar%blocksize_pdsyev<0) then
         call timing(iproc,'dirmin_dgesv','OF')
         inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
         power(1)=1
         call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
              order_taylor, 1, power, -8, &
              imode=2, &
              ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
              ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp_, check_accur=.true., &
              max_error=max_error, mean_error=mean_error)
         call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
    
         !!!DEBUG checking S^-1 etc.
         !!!test dense version of S^-1
         !!inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
         !!if (iproc==0)write(*,*) ''
         !!do itmp=0,20
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr0)
         !!   call system_clock(ncount1,ncount_rate,ncount_max)
         !!   maxerror=0.0d0
         !!   do itrials=1,50
         !!      call overlapPowerGeneral(iproc, nproc, itmp, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, &
         !!           inv_ovrlp, error, tmb%orbs)
         !!      maxerror=max(maxerror,error)
         !!   end do
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr1)
         !!   call system_clock(ncount2,ncount_rate,ncount_max)
         !!   if (iproc==0) then
         !!      write(*,*) 'order,maxerror,time',itmp,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
         !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
         !!   end if
         !!end do
         !!call f_free_ptr(inv_ovrlp)
    
         !!! test S^+/-1/2
         !!inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
         !!if (iproc==0)write(*,*) ''
         !!do itmp=0,20
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr0)
         !!   call system_clock(ncount1,ncount_rate,ncount_max)
         !!   maxerror=0.0d0
         !!   do itrials=1,50
         !!      call overlapPowerGeneral(iproc, nproc, itmp, 2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, &
         !!           inv_ovrlp, error, tmb%orbs)
         !!      maxerror=max(maxerror,error)
         !!   end do
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr1)
         !!   call system_clock(ncount2,ncount_rate,ncount_max)
         !!   if (iproc==0) then
         !!      write(*,*) 'order,maxerror,time',itmp,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
         !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
         !!   end if
         !!end do
    
         !!if (iproc==0)write(*,*) ''
         !!do itmp=0,20
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr0)
         !!   call system_clock(ncount1,ncount_rate,ncount_max)
         !!   maxerror=0.0d0
         !!   do itrials=1,50
         !!      call overlapPowerGeneral(iproc, nproc, itmp, -2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, &
         !!           inv_ovrlp, error, tmb%orbs)
         !!      maxerror=max(maxerror,error)
         !!   end do
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr1)
         !!   call system_clock(ncount2,ncount_rate,ncount_max)
         !!   if (iproc==0) then
         !!      write(*,*) 'order,maxerror,time',itmp,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
         !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
         !!   end if
         !!end do
    
         !!!test sparse version of S^-1
         !!deallocate(tmb%linmat%ovrlp%matrix)
    
         !!if (iproc==0)write(*,*) ''
         !!call f_free_ptr(inv_ovrlp)
         !!tmb%linmat%inv_ovrlp%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp%nvctr,id='tmb%linmat%inv_ovrlp%matrix_compr')
         !!do itmp=0,2
         !!   tmb%linmat%inv_ovrlp%parallel_compression=itmp
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr0)
         !!   call system_clock(ncount1,ncount_rate,ncount_max)
         !!   maxerror=0.0d0
         !!   do itrials=1,50
         !!      call overlapPowerGeneral(iproc, nproc, 1, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp, &
         !!           error, tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp)
         !!      maxerror=max(maxerror,error)
         !!   end do
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr1)
         !!   call system_clock(ncount2,ncount_rate,ncount_max)
         !!   if (iproc==0) then
         !!      write(*,*) 'order,maxerror,time',1,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
         !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
         !!   end if
         !!end do
         !!call f_free_ptr(tmb%linmat%inv_ovrlp%matrix_compr)
    
         !!!test sparse version of S^1/2
         !!if (iproc==0)write(*,*) ''
         !!nullify(inv_ovrlp)
         !!tmb%linmat%inv_ovrlp%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp%nvctr,id='tmb%linmat%inv_ovrlp%matrix_compr')
         !!do itmp=0,2
         !!   tmb%linmat%inv_ovrlp%parallel_compression=itmp
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr0)
         !!   call system_clock(ncount1,ncount_rate,ncount_max)
         !!   maxerror=0.0d0
         !!   do itrials=1,50
         !!      call overlapPowerGeneral(iproc, nproc, 1, 2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp, &
         !!           error, tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp)
         !!      maxerror=max(maxerror,error)
         !!   end do
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr1)
         !!   call system_clock(ncount2,ncount_rate,ncount_max)
         !!   if (iproc==0) then
         !!      write(*,*) 'order,maxerror,time',1,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
         !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
         !!   end if
         !!end do
         !!call f_free_ptr(tmb%linmat%inv_ovrlp%matrix_compr)
    
         !!!test sparse version of S^-1/2
         !!if (iproc==0)write(*,*) ''
         !!nullify(inv_ovrlp)
         !!tmb%linmat%inv_ovrlp%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp%nvctr,id='tmb%linmat%inv_ovrlp%matrix_compr')
         !!do itmp=0,2
         !!   tmb%linmat%inv_ovrlp%parallel_compression=itmp
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr0)
         !!   call system_clock(ncount1,ncount_rate,ncount_max)
         !!   maxerror=0.0d0
         !!   do itrials=1,50
         !!      call overlapPowerGeneral(iproc, nproc, 1, -2, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp, &
         !!           error, tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp)
         !!      maxerror=max(maxerror,error)
         !!   end do
         !!   call mpi_barrier(mpi_comm_world,ierr)
         !!   call cpu_time(tr1)
         !!   call system_clock(ncount2,ncount_rate,ncount_max)
         !!   if (iproc==0) then
         !!      write(*,*) 'order,maxerror,time',1,maxerror,real(tr1-tr0,kind=8)/real(itrials-1,kind=8),&
         !!           (dble(ncount2-ncount1)/dble(ncount_rate))/real(itrials-1,kind=8)
         !!   end if
         !!end do
         !!call f_free_ptr(tmb%linmat%inv_ovrlp%matrix_compr)
    
         !!call mpi_finalize(bigdft_mpi%mpi_comm)
         !!stop
         !!END DEBUG
    
         call timing(iproc,'dirmin_dgesv','ON')
         if (KSorbs%norbp>0) then
             !call dgemm('n', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, inv_ovrlp_%matrix(1,1,1), &
             !     tmb%orbs%norb, grad_cov(1,1), tmb%orbs%norb, 0.d0, grad(1,1), tmb%orbs%norb)
             do iorb=1,KSorbs%norbp
                 iiorb=KSorbs%isorb+iorb
                 if (KSorbs%spinsgn(iiorb)>0.d0) then
                     ispin=1
                 else
                     ispin=2
                 end if
                 call dgemm('n', 'n', tmb%linmat%smat(2)%nfvctr, 1, tmb%linmat%smat(2)%nfvctr, &
                      1.d0, inv_ovrlp_(1)%matrix(1,1,ispin), &
                      tmb%linmat%smat(2)%nfvctr, grad_cov(1,iorb), tmb%linmat%smat(2)%nfvctr, &
                      0.d0, grad(1,iorb), tmb%linmat%smat(2)%nfvctr)
             end do
         end if
         call f_free_ptr(inv_ovrlp)
      else
         info = 0 ! needed for when some processors have orbs%norbp=0
         grad_full=f_malloc((/tmb%linmat%smat(2)%nfvctr,KSorbs%norb/),id='grad_full')
         ! do allgather instead of allred so we can keep grad as per proc
         if(nproc > 1) then 
            call mpi_allgatherv(grad_cov, tmb%linmat%smat(2)%nfvctr*KSorbs%norbp, mpi_double_precision, grad_full, &
                 tmb%linmat%smat(2)%nfvctr*KSorbs%norb_par(:,0), tmb%linmat%smat(2)%nfvctr*KSorbs%isorb_par, &
                 mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         else
            call vcopy(tmb%linmat%smat(2)%nfvctr*KSorbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
         end if
    
         call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
              tmb%linmat%smat(2)%nfvctr, KSorbs%norb, tmb%linmat%ovrlp_%matrix, tmb%linmat%smat(2)%nfvctr, &
              grad_full, tmb%linmat%smat(2)%nfvctr, info)
         call vcopy(tmb%linmat%smat(2)%nfvctr*KSorbs%norbp,grad_full(1,KSorbs%isorb+1),1,grad(1,1),1)
    
         call f_free(grad_full)
         if(info/=0) then
            write(*,'(a,i0)') 'ERROR in dgesv: info=',info
            stop
         end if
      end if
    
      call deallocate_matrices(inv_ovrlp_(1))
    
      call timing(iproc,'dirmin_dgesv','OF')
      call f_release_routine()
    
    end subroutine calculate_coeff_gradient
    
    
    !!!> calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
    !!!! then grad=S^-1grad_cov
    !!subroutine calculate_coeff_gradient_extra(iproc,nproc,num_extra,tmb,order_taylor,max_inversion_error,KSorbs,grad_cov,grad)
    !!  use module_base
    !!  use module_types
    !!  use module_interfaces, only: calculate_density_kernel
    !!  use sparsematrix_base, only: matrices, sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), &
    !!                               matrices_null, allocate_matrices, deallocate_matrices
    !!  use sparsematrix, only: extract_taskgroup_inplace, gather_matrix_from_taskgroups_inplace
    !!  use matrix_operations, only: overlapPowerGeneral, check_taylor_order
    !!  use parallel_linalg, only: dgesv_parallel
    !!  implicit none
    !!
    !!  integer, intent(in) :: iproc, nproc, num_extra
    !!  integer,intent(inout) :: order_taylor
    !!  real(kind=8),intent(in) :: max_inversion_error
    !!  type(DFT_wavefunction), intent(inout) :: tmb
    !!  type(orbitals_data), intent(in) :: KSorbs
    !!  real(gp), dimension(tmb%linmat%smat(2)%nfvctr,tmb%orbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp
    !!
    !!  integer :: iorb, iiorb, info, ierr, ispin
    !!  integer, dimension(1) :: power
    !!  real(gp),dimension(:,:,:),allocatable :: sk, skhp, skh
    !!  real(gp),dimension(:,:),pointer ::  inv_ovrlp
    !!  real(kind=gp), dimension(:), allocatable:: occup_tmp
    !!  real(kind=gp), dimension(:,:), allocatable:: grad_full
    !!  real(kind=gp) :: max_error, mean_error
    !!  type(matrices),dimension(1) :: inv_ovrlp_
    !!
    !!  !!if (tmb%linmat%smat(2)%nspin==2) stop 'ERROR: calculate_coeff_gradient_extra not yet implemented for npsin==2'
    !!
    !!  call f_routine(id='calculate_coeff_gradient_extra')
    !!  call timing(iproc,'dirmin_lagmat1','ON')
    !!
    !!  inv_ovrlp_(1) = matrices_null()
    !!  call allocate_matrices(tmb%linmat%smat(3), allocate_full=.true., matname='inv_ovrlp_', mat=inv_ovrlp_(1))
    !!
    !!
    !!  occup_tmp=f_malloc(tmb%orbs%norb,id='occup_tmp')
    !!  call vcopy(tmb%orbs%norb,tmb%orbs%occup(1),1,occup_tmp(1),1)
    !!
    !!  call f_zero(tmb%orbs%norb,tmb%orbs%occup(1))
    !!  do iorb=1,KSorbs%norb+num_extra
    !!     tmb%orbs%occup(iorb)=1.0d0
    !!  end do
    !!
    !!  ! we have the kernel already, but need it to not contain occupations so recalculate here
    !!  !call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
    !!  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%smat(3), iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
    !!  !call uncompress_matrix(iproc,tmb%linmat%denskern)
    !!  !!call calculate_density_kernel_uncompressed (iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%kernel_%matrix)
    !!  !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
    !!  call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, &
    !!       tmb%coeff, tmb%linmat%smat(3), tmb%linmat%kernel_, .true.)
    !!  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
    !!
    !!
    !!  call vcopy(tmb%orbs%norb,occup_tmp(1),1,tmb%orbs%occup(1),1)
    !!  call f_free(occup_tmp)
    !!
    !!  sk=f_malloc0((/tmb%linmat%smat(3)%nfvctrp,tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nspin/), id='sk')
    !!
    !!  ! calculate I-S*K - first set sk to identity
    !!  do ispin=1,tmb%linmat%smat(3)%nspin
    !!     do iorb=1,tmb%linmat%smat(3)%nfvctrp
    !!        iiorb=mod(tmb%linmat%smat(3)%isfvctr+iorb-1,tmb%linmat%smat(3)%nfvctr)+1
    !!        sk(iorb,iiorb,ispin) = 1.d0
    !!     end do 
    !!  end do
    !!
    !!
    !!  if (tmb%linmat%smat(3)%nfvctrp>0) then
    !!     do ispin=1,tmb%linmat%smat(3)%nspin
    !!        call dgemm('t', 'n', tmb%linmat%smat(3)%nfvctrp, tmb%linmat%smat(3)%nfvctr, tmb%linmat%smat(3)%nfvctr, -1.d0, &
    !!             tmb%linmat%ovrlp_%matrix(1,tmb%linmat%smat(3)%isfvctr+1,ispin), tmb%linmat%smat(3)%nfvctr, &
    !!             tmb%linmat%kernel_%matrix(1,1,ispin), tmb%linmat%smat(3)%nfvctr, 1.d0, sk(1,1,ispin), tmb%linmat%smat(3)%nfvctrp)
    !!      end do
    !!  end if
    !!
    !!
    !!  ! coeffs and therefore kernel will change, so no need to keep it
    !!  call f_free_ptr(tmb%linmat%kernel_%matrix)
    !!
    !!  skhp=f_malloc((/tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nfvctrp,tmb%linmat%smat(3)%nspin/), id='skhp')
    !!
    !!  ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
    !!  if (tmb%linmat%smat(3)%nfvctrp>0) then
    !!     do ispin=1,tmb%linmat%smat(3)%nspin
    !!        call dgemm('t', 't', tmb%linmat%smat(3)%nfvctr, tmb%linmat%smat(3)%nfvctrp, tmb%linmat%smat(3)%nfvctr, &
    !!             1.d0, tmb%linmat%ham_%matrix(1,1,ispin), &
    !!             tmb%linmat%smat(3)%nfvctr, sk(1,1,ispin), tmb%linmat%smat(3)%nfvctrp, 0.d0, skhp(1,1,ispin), tmb%linmat%smat(3)%nfvctr)
    !!      end do
    !!  end if
    !!
    !!
    !!  call f_free(sk)
    !!
    !!  skh=f_malloc((/tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nfvctr,tmb%linmat%smat(3)%nspin/), id='skh')
    !!
    !!  call timing(iproc,'dirmin_lagmat1','OF')
    !!  call timing(iproc,'dirmin_lagmat2','ON')
    !!
    !!  ! gather together
    !!  if(nproc > 1) then
    !!     do ispin=1,tmb%linmat%smat(3)%nspin
    !!        call mpi_allgatherv(skhp(1,1,ispin), tmb%linmat%smat(3)%nfvctr*tmb%linmat%smat(3)%nfvctrp, &
    !!            mpi_double_precision, skh(1,1,ispin), &
    !!            tmb%linmat%smat(3)%nfvctr*tmb%linmat%smat(3)%nfvctr_par(:), &
    !!            tmb%linmat%smat(3)%nfvctr*tmb%linmat%smat(3)%isfvctr_par, &
    !!            mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    !!     end do
    !!  else
    !!     call vcopy(tmb%linmat%smat(3)%nfvctrp*tmb%linmat%smat(3)%nfvctr*tmb%linmat%smat(3)%nspin,skhp(1,1,1),1,skh(1,1,1),1)
    !!  end if
    !!
    !!
    !!  call timing(iproc,'dirmin_lagmat2','OF')
    !!  call timing(iproc,'dirmin_lagmat1','ON')
    !!
    !!  call f_free(skhp)
    !!
    !!  ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
    !!  if (tmb%linmat%smat(3)%nfvctrp>0) then
    !!     do iorb=1,tmb%linmat%smat(3)%nfvctrp
    !!        iiorb=tmb%linmat%smat(3)%isfvctr+iorb
    !!        if (tmb%orbs%spinsgn(iiorb)>0.d0) then
    !!            ispin=1
    !!        else
    !!            ispin=2
    !!        end if
    !!        call dgemm('t', 'n', tmb%linmat%smat(3)%nfvctr, 1, tmb%linmat%smat(3)%nfvctr, 1.d0, skh(1,1,ispin), &
    !!             tmb%linmat%smat(3)%nfvctr, tmb%coeff(1,tmb%linmat%smat(3)%isfvctr+iorb), tmb%linmat%smat(3)%nfvctr, &
    !!             0.d0, grad_cov(1,iorb), tmb%linmat%smat(3)%nfvctr)
    !!      end do
    !!  end if
    !!
    !!  call f_free(skh)
    !!
    !!  ! multiply by f_i to get grad_i^a
    !!  do iorb=1,tmb%orbs%norbp
    !!     iiorb=tmb%orbs%isorb+iorb
    !!     grad_cov(:,iorb)=grad_cov(:,iorb)*tmb%orbs%occup(iiorb)
    !!  end do
    !!
    !!  call timing(iproc,'dirmin_lagmat1','OF')
    !!  call timing(iproc,'dirmin_dgesv','ON')
    !!
    !!
    !!  info = 0 ! needed for when some processors have orbs%norbp=0
    !!  ! Solve the linear system ovrlp*grad=grad_cov
    !!  if(tmb%orthpar%blocksize_pdsyev<0) then
    !!     !! keep the covariant gradient to calculate fnrm correctly
    !!     !call vcopy(tmb%orbs%norb*tmb%orbs%norbp,grad_cov,1,grad,1)
    !!     !if (tmb%orbs%norbp>0) then
    !!     !   ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
    !!     !   call dgesv(tmb%orbs%norb, tmb%orbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
    !!     !        grad(1,1), tmb%orbs%norb, info)
    !!     !   call f_free(ipiv)
    !!     !end if
    !!     !!inv_ovrlp=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
    !!     power(1)=1
    !!     call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
    !!          order_taylor, 1, power, -8, &
    !!          imode=2, ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
    !!          ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp_, check_accur=.true., &
    !!          max_error=max_error, mean_error=mean_error)
    !!     call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
    !!
    !!     if (tmb%orbs%norbp>0) then
    !!        call dgemm('n', 'n', tmb%linmat%smat(3)%nfvctr, tmb%orbs%norbp, tmb%linmat%smat(3)%nfvctr, 1.d0, inv_ovrlp_(1)%matrix(1,1,1), &
    !!             tmb%linmat%smat(3)%nfvctr, grad_cov(1,1), tmb%linmat%smat(3)%nfvctr, 0.d0, grad(1,1), tmb%linmat%smat(3)%nfvctr)
    !!     end if
    !!     !!call f_free_ptr(inv_ovrlp)
    !!  else
    !!      grad_full=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='grad_full')
    !!      ! do allgather instead of allred so we can keep grad as per proc
    !!      if(nproc > 1) then 
    !!         call mpi_allgatherv(grad_cov, tmb%linmat%smat(3)%nfvctr*tmb%orbs%norbp, &
    !!              mpi_double_precision, grad_full, &
    !!              tmb%linmat%smat(3)%nfvctr*tmb%orbs%norb_par(:,0), &
    !!              tmb%linmat%smat(3)%nfvctr*tmb%orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    !!      else
    !!         call vcopy(tmb%linmat%smat(3)%nfvctr*tmb%orbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
    !!      end if
    !!      !call mpiallred(grad(1,1), tmb%orbs%norb*tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    !!
    !!      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
    !!           tmb%orbs%norb, tmb%orbs%norb, tmb%linmat%ovrlp_%matrix, tmb%orbs%norb, grad_full, tmb%orbs%norb, info)
    !!
    !!      call vcopy(tmb%linmat%smat(3)%nfvctr*tmb%orbs%norbp,grad_full(1,tmb%orbs%isorb+1),1,grad(1,1),1)
    !!
    !!      call f_free(grad_full)
    !!
    !!
    !!  end if
    !!
    !!  if(info/=0) then
    !!      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
    !!      stop
    !!  end if
    !!
    !!  call deallocate_matrices(inv_ovrlp_(1))
    !!
    !!  call timing(iproc,'dirmin_dgesv','OF')
    !!  call f_release_routine()
    !!
    !!end subroutine calculate_coeff_gradient_extra
    
    
    subroutine precondition_gradient_coeff(ntmb, norb, ham, ovrlp, grad)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in):: ntmb, norb
      real(8),dimension(ntmb,ntmb),intent(in):: ham, ovrlp
      real(8),dimension(ntmb,norb),intent(inout):: grad
      
      ! Local variables
      integer:: iorb, itmb, jtmb, info, istat, iall
      complex(8),dimension(:,:),allocatable:: mat
      complex(8),dimension(:,:),allocatable:: rhs
      integer,dimension(:),allocatable:: ipiv
      character(len=*),parameter:: subname='precondition_gradient_coeff'
      
      mat = f_malloc((/ ntmb, ntmb /),id='mat')
      rhs = f_malloc((/ ntmb, norb /),id='rhs')
      
      ! Build the matrix to be inverted
      do itmb=1,ntmb
          do jtmb=1,ntmb
              mat(jtmb,itmb) = cmplx(ham(jtmb,itmb)+.5d0*ovrlp(jtmb,itmb),0.d0,kind=8)
          end do
          mat(itmb,itmb)=mat(itmb,itmb)+cmplx(0.d0,-1.d-1,kind=8)
          !mat(itmb,itmb)=mat(itmb,itmb)-cprec
      end do
      do iorb=1,norb
          do itmb=1,ntmb
              rhs(itmb,iorb)=cmplx(grad(itmb,iorb),0.d0,kind=8)
          end do
      end do
      
      ipiv = f_malloc(ntmb,id='ipiv')
      
      call zgesv(ntmb, norb, mat(1,1), ntmb, ipiv, rhs(1,1), ntmb, info)
      if(info/=0) then
          stop 'ERROR in zgesv'
      end if
      !call vcopy(nel, rhs(1), 1, grad(1), 1)
      do iorb=1,norb
          do itmb=1,ntmb
              grad(itmb,iorb)=real(rhs(itmb,iorb))
          end do
      end do
      
      call f_free(ipiv)
      call f_free(mat)
      call f_free(rhs)
    
    end subroutine precondition_gradient_coeff
    
    
    
    subroutine initialize_DIIS_coeff(isx, ldiis)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in):: isx
      type(localizedDIISParameters),intent(inout):: ldiis
      
      ! Local variables
      character(len=*),parameter:: subname='initialize_DIIS_coeff'
        
      ldiis%isx=isx
      ldiis%is=0
      ldiis%switchSD=.false.
      ldiis%trmin=1.d100
      ldiis%trold=1.d100
      ldiis%alpha_coeff=0.1d0
    
    end subroutine initialize_DIIS_coeff
    
    
    subroutine allocate_DIIS_coeff(tmb, ldiis)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      type(DFT_wavefunction),intent(in):: tmb
      type(localizedDIISParameters),intent(inout):: ldiis
      
      ! Local variables
      integer:: ii, istat
      character(len=*),parameter:: subname='allocate_DIIS_coeff'
    
      ldiis%mat = f_malloc_ptr((/ldiis%isx,ldiis%isx,tmb%orbs%norbp/),id='ldiis%mat')
    
      ii=ldiis%isx*tmb%orbs%norb*tmb%orbs%norbp
      ldiis%phiHist = f_malloc_ptr(ii,id='ldiis%phiHist')
      ldiis%hphiHist = f_malloc_ptr(ii,id='ldiis%hphiHist')
    
    end subroutine allocate_DIIS_coeff


end module get_kernel
