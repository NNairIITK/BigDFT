!> @file
!! Linear version: Handle local basis set
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine get_coeff(iproc,nproc,scf_mode,orbs,at,rxyz,denspot,GPU,infoCoeff,&
    energs,nlpsp,SIC,tmb,fnrm,calculate_overlap_matrix,invert_overlap_matrix,communicate_phi_for_lsumrho,&
    calculate_ham,extra_states,itout,it_scc,it_cdft,order_taylor,max_inversion_error,&
    calculate_KS_residue,calculate_gap,energs_work,remove_coupling_terms,factor,&
    pexsi_npoles,pexsi_mumin,pexsi_mumax,pexsi_mu,pexsi_temperature, pexsi_tol_charge,&
    convcrit_dmin,nitdmin,curvefit_dmin,ldiis_coeff,reorder,cdft,updatekernel)
  use module_base
  use module_types
  use module_interfaces, only: LocalHamiltonianApplication, SynchronizeHamiltonianApplication, optimize_coeffs
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use constrained_dft
  use diis_sd_optimization
  use yaml_output
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, start_onesided_communication
  use rhopotential, only: full_local_potential
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, sparsematrix_malloc, &
                               DENSE_FULL, DENSE_PARALLEL, DENSE_MATMUL, assignment(=), SPARSE_FULL
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          extract_taskgroup_inplace, uncompress_matrix_distributed2, gather_matrix_from_taskgroups, &
                          extract_taskgroup, uncompress_matrix2, &
                          write_sparsematrix, delete_coupling_terms, transform_sparse_matrix, &
                          max_asymmetry_of_matrix
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
  use pexsi, only: pexsi_driver
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, scf_mode, itout, it_scc, it_cdft
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
  logical,intent(in):: calculate_overlap_matrix, invert_overlap_matrix
  logical,intent(in):: communicate_phi_for_lsumrho
  logical,intent(in) :: calculate_ham, calculate_KS_residue, calculate_gap
  type(work_mpiaccumulate),intent(inout) :: energs_work
  logical,intent(in) :: remove_coupling_terms
  real(kind=8), intent(in) :: factor
  integer,intent(in) :: pexsi_npoles
  real(kind=8),intent(in) :: pexsi_mumin,pexsi_mumax,pexsi_mu,pexsi_temperature, pexsi_tol_charge
  type(DIIS_obj),intent(inout),optional :: ldiis_coeff ! for dmin only
  integer, intent(in), optional :: nitdmin ! for dmin only
  real(kind=gp), intent(in), optional :: convcrit_dmin ! for dmin only
  logical, intent(in), optional :: curvefit_dmin ! for dmin only
  type(cdft_data),intent(inout),optional :: cdft
  integer, intent(in) :: extra_states
  logical, optional, intent(in) :: reorder
  logical, optional, intent(in) :: updatekernel

  ! Local variables 
  integer :: iorb, info, ishift, ispin, ii, jorb, i, ishifts, ishiftm, jproc
  real(kind=8),dimension(:),allocatable :: hpsit_c, hpsit_f, eval, tmparr1, tmparr2, tmparr, ovrlp_large, ham_large
  real(kind=8),dimension(:,:),allocatable :: ovrlp_fullp, tempmat
  real(kind=8),dimension(:,:,:),allocatable :: matrixElements
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  logical :: update_kernel
  character(len=*),parameter :: subname='get_coeff'
  real(kind=gp) :: tmprtr
  real(kind=8) :: max_deviation, mean_deviation, KSres, max_deviation_p,  mean_deviation_p, maxdiff, tt
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
      call start_onesided_communication(iproc, nproc, denspot%dpbox%ndims(1), denspot%dpbox%ndims(2), &
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

      call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%psi, tmb%ham_descr%psi)

      if (tmb%ham_descr%npsidim_orbs > 0) call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))

      call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,tmb%hpsi,energs%eproj,tmb%paw)
      ! only kinetic as waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,&
           & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,3,denspot%xc,&
           & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
           & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
           & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
      !!do i=1,tmb%ham_descr%comgp%nrecvbuf
      !!    write(8000+iproc,'(a,i8,es16.6)') 'i, recvbuf(i)', i, tmb%ham_descr%comgp%recvbuf(i)
      !!end do
      !call wait_p2p_communication(iproc, nproc, tmb%ham_descr%comgp)
      ! only potential
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
           tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%m, tmb%linmat%ham_)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)


!!  call diagonalize_subset(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%ovrlp_, tmb%linmat%m, tmb%linmat%ham_)
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
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
      !call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
  end if

  ovrlp_fullp = sparsematrix_malloc(tmb%linmat%l,iaction=DENSE_PARALLEL,id='ovrlp_fullp')
  max_deviation=0.d0
  mean_deviation=0.d0
  do ispin=1,tmb%linmat%s%nspin
      ishift=(ispin-1)*tmb%linmat%s%nvctrp_tg
      call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
           tmb%linmat%ovrlp_%matrix_compr(ishift+1:), ovrlp_fullp)
      call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, &
           tmb%linmat%s%isfvctr, ovrlp_fullp, &
           tmb%linmat%s, max_deviation_p, mean_deviation_p)
      max_deviation = max_deviation + max_deviation_p/real(tmb%linmat%s%nspin,kind=8)
      mean_deviation = mean_deviation + mean_deviation_p/real(tmb%linmat%s%nspin,kind=8)
  end do
  call f_free(ovrlp_fullp)
  if (iproc==0) then
      call yaml_map('max dev from unity',max_deviation,fmt='(es9.2)')
      call yaml_map('mean dev from unity',mean_deviation,fmt='(es9.2)')
  end if

  ! Post the p2p communications for the density. (must not be done in inputguess)
  if(communicate_phi_for_lsumrho) then
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
           tmb%orbs, tmb%psi, tmb%collcom_sr)
  end if

  ! CDFT: add V*w_ab to Hamiltonian here - assuming ham and weight matrix have the same sparsity...
  if (present(cdft)) then
     call timing(iproc,'constraineddft','ON')
     call daxpy(tmb%linmat%m%nvctr,cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmb%linmat%ham_%matrix_compr,1)
     call timing(iproc,'constraineddft','OF') 
  end if

  if (remove_coupling_terms) then
      call delete_coupling_terms(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smmd, tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr)
      call delete_coupling_terms(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smmd, tmb%linmat%m, tmb%linmat%ham_%matrix_compr)
      call delete_coupling_terms(iproc, nproc, bigdft_mpi%mpi_comm, &
           tmb%linmat%smmd, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr)
  end if

  if (scf_mode/=LINEAR_FOE .and. scf_mode/=LINEAR_PEXSI) then
      tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
      !call yaml_map('Ham1com',tmb%linmat%ham_%matrix_compr)
      !call yaml_map('ovrlpcom',tmb%linmat%ovrlp_%matrix_compr)

      tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      !call uncompress_matrix(iproc, tmb%linmat%m, &
      !     inmat=tmb%linmat%ham_%matrix_compr, outmat=tmb%linmat%ham_%matrix)
      !call uncompress_matrix(iproc, tmb%linmat%s, &
      !     inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
      do ispin=1,tmb%linmat%m%nspin
          ishifts = (ispin-1)*tmb%linmat%s%nvctrp_tg
          ishiftm = (ispin-1)*tmb%linmat%m%nvctrp_tg
          call f_zero(tmb%linmat%m%nfvctr**2, tmb%linmat%ham_%matrix(1,1,ispin))
          tempmat = sparsematrix_malloc(tmb%linmat%m, iaction=DENSE_PARALLEL, id='tempmat')
          call uncompress_matrix_distributed2(iproc, tmb%linmat%m, DENSE_PARALLEL, &
               tmb%linmat%ham_%matrix_compr(ishiftm+1:ishiftm+tmb%linmat%m%nvctrp_tg), tempmat)
          if (tmb%linmat%m%nfvctrp>0) then
              call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctrp, tempmat(1,1), 1, &
                   tmb%linmat%ham_%matrix(1,tmb%linmat%m%isfvctr+1,ispin), 1)
          end if
          call f_free(tempmat)
          if (nproc>1) then
              call mpiallred(tmb%linmat%ham_%matrix(1,1,ispin), tmb%linmat%m%nfvctr**2, &
                   mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if

          call f_zero(tmb%linmat%s%nfvctr**2, tmb%linmat%ovrlp_%matrix(1,1,ispin))
          tempmat = sparsematrix_malloc(tmb%linmat%s, iaction=DENSE_PARALLEL, id='tempmat')
          call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
               tmb%linmat%ovrlp_%matrix_compr(ishifts+1:), tempmat)
          if (tmb%linmat%m%nfvctrp>0) then
              call vcopy(tmb%linmat%s%nfvctr*tmb%linmat%s%nfvctrp, tempmat(1,1), 1, &
                   tmb%linmat%ovrlp_%matrix(1,tmb%linmat%s%isfvctr+1,ispin), 1)
          end if
          call f_free(tempmat)
          if (nproc>1) then
              call mpiallred(tmb%linmat%ovrlp_%matrix(1,1,ispin), tmb%linmat%s%nfvctr**2, &
                   mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
      end do
  end if

  ! Diagonalize the Hamiltonian.
!  if (iproc==0) call yaml_sequence_open('kernel method')
  if(scf_mode==LINEAR_MIXPOT_SIMPLE .or. scf_mode==LINEAR_MIXDENS_SIMPLE) then
      ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
      matrixElements = f_malloc((/ tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr,2 /),id='matrixElements')
      eval = f_malloc(tmb%linmat%l%nfvctr,id='eval')

      do ispin=1,tmb%linmat%s%nspin
          !if (ispin==1) then
          !    ishift=0
          !    !ii=orbs%norbu
          !    tmb%linmat%
          !else
          !    ishift=orbs%norbu
          !    ii=orbs%norbd
          !end if
          !ishift=(ispin-1)*tmb%linmat%s%nfvctr
          call vcopy(tmb%linmat%m%nfvctr**2, tmb%linmat%ham_%matrix(1,1,ispin), 1, matrixElements(1,1,1), 1)
          call vcopy(tmb%linmat%m%nfvctr**2, tmb%linmat%ovrlp_%matrix(1,1,ispin), 1, matrixElements(1,1,2), 1)
          if (iproc==0) call yaml_map('method','diagonalization')
          if(tmb%orthpar%blocksize_pdsyev<0) then
              if (iproc==0) call yaml_map('mode','sequential')
              !!if (iproc==0) then
              !!    do iorb=1,tmb%linmat%m%nfvctr
              !!        do jorb=1,tmb%linmat%m%nfvctr
              !!            write(690+ispin,'(a,2i8,2f15.8)') 'iorb, jorb, vals', iorb, jorb, matrixElements(jorb,iorb,:)
              !!        end do
              !!    end do
              !!end if
              call diagonalizeHamiltonian2(iproc, tmb%linmat%m%nfvctr, &
                   matrixElements(1,1,1), matrixElements(1,1,2), eval)
          else
              if (iproc==0) call yaml_map('mode','parallel')
              call dsygv_parallel(iproc, nproc, bigdft_mpi%mpi_comm, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%nproc_pdsyev, &
                   1, 'v', 'l', tmb%linmat%m%nfvctr, &
                   matrixElements(1,1,1), tmb%linmat%m%nfvctr, matrixElements(1,1,2), tmb%linmat%m%nfvctr, &
                   eval, info)
          end if

          ! Broadcast the results (eigenvectors and eigenvalues) from task 0 to
          ! all other tasks (in this way avoiding that different MPI tasks have different values)
          if (nproc>1) then
              if (iproc==0) call yaml_mapping_open('Cross-check among MPI tasks')
              call mpibcast(matrixElements(:,:,1), comm=bigdft_mpi%mpi_comm, maxdiff=maxdiff)
              if (iproc==0) call yaml_map('max diff of eigenvectors',maxdiff,fmt='(es8.2)')
              call mpibcast(eval, comm=bigdft_mpi%mpi_comm, maxdiff=maxdiff)
              if (iproc==0) call yaml_map('max diff of eigenvalues',maxdiff,fmt='(es8.2)')
              if (iproc==0) call yaml_mapping_close()
          end if
          !!if (iproc==0) then
          !!    do iorb=1,tmb%orbs%norb
          !!        write(*,*) 'ind, val', iorb, matrixElements(iorb,orbs%norb+1,1)
          !!    end do
          !!end if
          !if (iproc==0) write(*,'(a,3i6,100f9.2)') 'ispin, ishift+1, ishift+ii, evals', ispin, ishift+1, ishift+ii, tmb%orbs%eval(ishift+1:ishift+ii)

          ! copy all the eigenvalues
          !tmb%orbs%eval(ishift+1:ishift+ii) = eval(1:ii-ishift)
          call vcopy(tmb%linmat%m%nfvctr, eval(1), 1, tmb%orbs%eval((ispin-1)*tmb%linmat%m%nfvctr+1), 1)
          ! copy the eigenvalues of the occupied states
          if (ispin==1) then
              call vcopy(orbs%norbu, eval(1), 1, orbs%eval(1), 1)
          else
              call vcopy(orbs%norbd, eval(1), 1, orbs%eval(orbs%norbu+1), 1)
          end if

          ! Make sure that the eigenvectors have the same sign on all MPI tasks.
          ! To do so, ensure that the first entry is always positive.
          do iorb=1,tmb%linmat%m%nfvctr
              if (matrixElements(1,iorb,1)<0.d0) then
                  call dscal(tmb%linmat%m%nfvctr, -1.d0, matrixElements(1,iorb,1), 1)
              end if
          end do

          ! Copy the diagonalized matrix to the coeff array.
          ! In principle I would prefer to copy orbs%norbu/orbs%norbd states.
          ! However this is not possible since the extra states are not included in there (WHY?!)
          ! Therefore as a workaround I use the following dirty solution with different cases.
          if (tmb%linmat%l%nspin/=1) then
              if (extra_states>0) stop 'extra states and spin polarization not possible at the moment'
              ! Only copy the occupied states
              if (ispin==1) then
                  call vcopy(orbs%norbu*tmb%linmat%m%nfvctr, matrixElements(1,1,1), 1, tmb%coeff(1,1), 1)
              else if (ispin==2) then
                  call vcopy(orbs%norbd*tmb%linmat%m%nfvctr, matrixElements(1,1,1), 1, tmb%coeff(1,orbs%norbu+1), 1)
              end if
          else
              ! Copy all states
              call vcopy(tmb%orbs%norb*tmb%linmat%m%nfvctr, matrixElements(1,1,1), 1, tmb%coeff(1,1), 1)
          end if
          infoCoeff=0


          ! keep the eigenvalues for the preconditioning - instead should take h_alpha,alpha for both cases
          ! instead just use -0.5 everywhere
          !tmb%orbs%eval(:) = -0.5_dp
      end do

      call f_free(eval)

      !!if (iproc==0) then
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,tmb%linmat%m%nfvctr
      !!            if (orbs%spinsgn(iorb)>0.d0) then
      !!                write(620,*) 'iorb, jorb, val', iorb, jorb, tmb%coeff(jorb,iorb)
      !!            else
      !!                write(621,*) 'iorb, jorb, val', iorb, jorb, tmb%coeff(jorb,iorb)
      !!            end if
      !!        end do
      !!    end do
      !!end if

      call f_free(matrixElements)
  else if (scf_mode==LINEAR_DIRECT_MINIMIZATION) then
     if(.not.present(ldiis_coeff)) &
          call f_err_throw('ldiis_coeff must be present for scf_mode==LINEAR_DIRECT_MINIMIZATION',&
          err_name='BIGDFT_RUNTIME_ERROR')
     ! call routine which updates coeffs for tmb%orbs%norb or orbs%norb depending on whether or not extra states are required
     if (iproc==0) call yaml_map('method','directmin')
     if (extra_states>0) then
        call optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, convcrit_dmin, nitdmin, energs%ebs, &
             curvefit_dmin, factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder, extra_states)
     else
        call optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, convcrit_dmin, nitdmin, energs%ebs, &
             curvefit_dmin, factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder)
     end if
  end if

  ! CDFT: subtract V*w_ab from Hamiltonian so that we are calculating the correct energy
  if (present(cdft)) then
     call timing(iproc,'constraineddft','ON')
     tmparr = sparsematrix_malloc(tmb%linmat%m,iaction=SPARSE_FULL,id='tmparr')
     call gather_matrix_from_taskgroups(iproc, nproc, bigdft_mpi%mpi_comm, &
          tmb%linmat%m, tmb%linmat%ham_%matrix_compr, tmparr)
     call daxpy(tmb%linmat%m%nvctr*tmb%linmat%m%nspin,-cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmparr,1)
     call extract_taskgroup(tmb%linmat%m, tmparr, tmb%linmat%ham_%matrix_compr)
     call f_free(tmparr)
     call timing(iproc,'constraineddft','OF') 
  end if

  if (scf_mode/=LINEAR_FOE .and. scf_mode/=LINEAR_PEXSI) then
      ! Calculate the band structure energy and update kernel
      if (scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
         !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
         !!call extract_taskgroup_inplace(tmb%linmat%m, tmb%linmat%ham_)
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
              tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
              tmb%coeff,orbs,tmb%orbs,update_kernel)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
      else if (present(cdft)) then
         ! for directmin we have the kernel already, but only the CDFT function not actual energy for CDFT
         !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
         !!call extract_taskgroup_inplace(tmb%linmat%m, tmb%linmat%ham_)
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
              tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
              tmb%coeff,orbs,tmb%orbs,.false.)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
      end if
      call f_free_ptr(tmb%linmat%ham_%matrix)
      call f_free_ptr(tmb%linmat%ovrlp_%matrix)

      !!if (iproc==0) then
      !!    ii=0
      !!    do ispin=1,tmb%linmat%l%nspin
      !!        do iorb=1,tmb%linmat%l%nvctr
      !!            ii=ii+1
      !!            if (ispin==1) then
      !!                write(630,'(a,2i8,f14.7)') 'ispin, iorb, val', ispin, iorb, tmb%linmat%kernel_%matrix_compr(ii)
      !!            else
      !!                write(631,'(a,2i8,f14.7)') 'ispin, iorb, val', ispin, iorb, tmb%linmat%kernel_%matrix_compr(ii)
      !!            end if
      !!        end do
      !!    end do
      !!end if

  else ! foe or pexsi

      !same as for directmin/diag
      ! CDFT: add V*w_ab to Hamiltonian here - assuming ham and weight matrix have the same sparsity...
      if (present(cdft)) then
         call timing(iproc,'constraineddft','ON')
         call daxpy(tmb%linmat%m%nvctr,cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmb%linmat%ham_%matrix_compr,1)
         call timing(iproc,'constraineddft','OF') 
      end if

      ! TEMPORARY #################################################
      if (calculate_gap) then
          tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
          !!tmparr = sparsematrix_malloc(tmb%linmat%m,iaction=SPARSE_FULL,id='tmparr')
          !!call gather_matrix_from_taskgroups(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_%matrix_compr, tmparr)
          call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%m, tmb%linmat%ham_%matrix_compr, tmb%linmat%ham_%matrix)
          !!call f_free(tmparr)
          tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
          call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
          ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
          matrixElements=f_malloc((/tmb%orbs%norb,tmb%orbs%norb,2/),id='matrixElements')
          !SM: need to fix the spin here
          call vcopy(tmb%orbs%norb**2, tmb%linmat%ham_%matrix(1,1,1), 1, matrixElements(1,1,1), 1)
          call vcopy(tmb%orbs%norb**2, tmb%linmat%ovrlp_%matrix(1,1,1), 1, matrixElements(1,1,2), 1)
          call diagonalizeHamiltonian2(iproc, tmb%orbs%norb, matrixElements(1,1,1), matrixElements(1,1,2), tmb%orbs%eval)
          if (iproc==0) call yaml_map('gap',tmb%orbs%eval(orbs%norb+1)-tmb%orbs%eval(orbs%norb))
          if (iproc==0) call yaml_map('lowest eigenvalue',tmb%orbs%eval(1))
          if (iproc==0) call yaml_map('highest eigenvalue',tmb%orbs%eval(tmb%orbs%norb))
          call f_free(matrixElements)
          call f_free_ptr(tmb%linmat%ham_%matrix)
          call f_free_ptr(tmb%linmat%ovrlp_%matrix)
      end if
      ! END TEMPORARY #############################################

      tmprtr=0.d0

      !!if (iproc==0) then
      !!    call yaml_map('write overlap matrix',.true.)
      !!    call write_sparsematrix('overlap.dat', tmb%linmat%s, tmb%linmat%ovrlp_)
      !!end if
      if (scf_mode==LINEAR_PEXSI) then
          if (iproc==0) call yaml_map('method','PEXSI')
          !call write_pexsi_matrices(iproc, nproc, tmb%linmat%m, tmb%linmat%s, tmb%linmat%ham_%matrix_compr, tmb%linmat%ovrlp_%matrix_compr)
          row_ind = f_malloc(tmb%linmat%l%nvctr,id='row_ind')
          col_ptr = f_malloc(tmb%linmat%l%nfvctr,id='col_ptr')
          call sparsebigdft_to_ccs(tmb%linmat%l%nfvctr, tmb%linmat%l%nvctr, tmb%linmat%l%nseg, tmb%linmat%l%keyg, row_ind, col_ptr)
          ! AT the moment not working for nspin>1
          ovrlp_large = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='ovrlp_large')
          ham_large = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='ham_large')
          !!call transform_sparsity_pattern(tmb%linmat%l%nfvctr, tmb%linmat%s%nvctr, 0, &
          !!     tmb%linmat%s%nseg, tmb%linmat%s%keyv, tmb%linmat%s%keyg, tmb%linmat%s%smmm%line_and_column_mm, &
          !!     tmb%linmat%l%nvctr, 0, tmb%linmat%l%smmm%nseg, &
          !!     tmb%linmat%l%smmm%keyv, tmb%linmat%l%smmm%keyg, &
          !!     tmb%linmat%l%smmm%istsegline, 'small_to_large', tmb%linmat%ovrlp_%matrix_compr, ovrlp_large)
          !!call transform_sparsity_pattern(tmb%linmat%l%nfvctr, tmb%linmat%m%smmm%nvctrp_mm, tmb%linmat%m%smmm%isvctr_mm, &
          !!     tmb%linmat%m%nseg, tmb%linmat%m%keyv, tmb%linmat%m%keyg, tmb%linmat%m%smmm%line_and_column_mm, &
          !!     tmb%linmat%l%smmm%nvctrp, tmb%linmat%l%smmm%isvctr, tmb%linmat%l%smmm%nseg, &
          !!     tmb%linmat%l%smmm%keyv, tmb%linmat%l%smmm%keyg, &
          !!     tmb%linmat%l%smmm%istsegline, 'small_to_large', tmb%linmat%ham_%matrix_compr, ham_large)
          if (tmb%linmat%s%ntaskgroup/=1 .or. tmb%linmat%m%ntaskgroup/=1 .or. tmb%linmat%l%ntaskgroup/=1) then
              call f_err_throw('PEXSI is not yet tested with matrix taskgroups', err_name='BIGDFT_RUNTIME_ERROR')
          end if
          call transform_sparse_matrix(iproc, tmb%linmat%s, tmb%linmat%l, SPARSE_FULL, 'small_to_large', &
               smat_in=tmb%linmat%ovrlp_%matrix_compr, lmat_out=ovrlp_large)
          call transform_sparse_matrix(iproc, tmb%linmat%m, tmb%linmat%l, SPARSE_FULL, 'small_to_large', &
               smat_in=tmb%linmat%ham_%matrix_compr, lmat_out=ham_large)
          !write(*,*) 'iproc, ham_large', iproc, ham_large
          call pexsi_driver(iproc, nproc, tmb%linmat%l%nfvctr, tmb%linmat%l%nvctr, row_ind, col_ptr, &
               ham_large, ovrlp_large, foe_data_get_real(tmb%foe_obj,"charge",1), pexsi_npoles, &
               pexsi_mumin, pexsi_mumax, pexsi_mu, pexsi_temperature, pexsi_tol_charge, &
               tmb%linmat%kernel_%matrix_compr, energs%ebs)

          call f_free(ovrlp_large)
          call f_free(ham_large)
          call f_free(row_ind)
          call f_free(col_ptr)
      else if (scf_mode==LINEAR_FOE) then
          if (iproc==0) call yaml_map('method','FOE')
          !call fermi_operator_expansion_new(iproc, nproc, &
          !     energs%ebs, &
          !     invert_overlap_matrix, 2, &
          !     tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, tmb%linmat%ham_, &
          !     tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), tmb%linmat%kernel_, tmb%foe_obj)
          call max_asymmetry_of_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%m, tmb%linmat%ham_%matrix_compr, tt)
          if (iproc==0) call yaml_map('max assymetry of H',tt)
          call max_asymmetry_of_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr, tt)
          if (iproc==0) call yaml_map('max assymetry of S',tt)
          call matrix_fermi_operator_expansion(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%foe_obj, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
               tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%ovrlppowers_(2), tmb%linmat%kernel_, &
               energs%ebs, invert_overlap_matrix, 2)
          call max_asymmetry_of_matrix(iproc, nproc, bigdft_mpi%mpi_comm, &
               tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, tt)
          if (iproc==0) call yaml_map('max assymetry of K',tt)

          !!call fermi_operator_expansion(iproc, nproc, &
          !!     energs%ebs, order_taylor, max_inversion_error, &
          !!     invert_overlap_matrix, 2, &
          !!     trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'-'//trim(adjustl(yaml_toa(it_cdft,fmt='(i3.3)')))&
          !!     //'-'//trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))), &
          !!     tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, tmb%linmat%ham_, &
          !!     tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), tmb%linmat%kernel_, tmb%foe_obj)
      end if

      ! Eigenvalues not available, therefore take -.5d0
      tmb%orbs%eval=-.5d0

      !same as for directmin/diag
      ! CDFT: subtract V*w_ab from Hamiltonian so that we are calculating the correct energy
      if (present(cdft)) then
         call timing(iproc,'constraineddft','ON')
         tmparr = sparsematrix_malloc(tmb%linmat%m,iaction=SPARSE_FULL,id='tmparr')
         call gather_matrix_from_taskgroups(iproc, nproc, bigdft_mpi%mpi_comm, &
              tmb%linmat%m, tmb%linmat%ham_%matrix_compr, tmparr)
         call daxpy(tmb%linmat%m%nvctr*tmb%linmat%m%nspin,-cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmparr,1)
         call extract_taskgroup(tmb%linmat%m, tmparr, tmb%linmat%ham_%matrix_compr)
         call f_free(tmparr)
         call timing(iproc,'constraineddft','OF') 

         ! we have ebs and the kernel already, but only the CDFT function not actual energy for CDFT
         ! maybe pass both Hamiltonians to FOE to save calculation ebs twice?
         !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
         !!call extract_taskgroup_inplace(tmb%linmat%m, tmb%linmat%ham_)
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
              tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
              tmb%coeff,orbs,tmb%orbs,.false.)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
         !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)

      end if

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



subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
        fnrm_tmb,infoBasisFunctions,nlpsp,scf_mode,ldiis,SIC,tmb,energs_base,do_iterative_orthogonalization,sf_per_type,&
    nit_precond,target_function,&
    correction_orthoconstraint,nit_basis,&
    ratio_deltas,ortho_on,extra_states,itout,conv_crit,experimental_mode,early_stop,&
    gnrm_dynamic, min_gnrm_for_dynamic, can_use_ham, order_taylor, max_inversion_error, kappa_conv, &
    correction_co_contra, &
    precond_convol_workarrays, precond_workarrays, &
    wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi, fnrm, energs_work, frag_calc, &
    cdft, input_frag, ref_frags)
  !
  ! Purpose:
  ! ========
  !   Calculates the localized basis functions phi. These basis functions are obtained by adding a
  !   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
  !   determined by minimizing the trace until the gradient norm is below the convergence criterion.
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, only: LocalHamiltonianApplication, SynchronizeHamiltonianApplication, &
       & calculate_density_kernel, hpsitopsi, write_energies
  use communications_base, only: work_transpose, TRANSPOSE_FULL, TRANSPOSE_POST, TRANSPOSE_GATHER
  use communications, only: transpose_localized, untranspose_localized, start_onesided_communication, &
                            synchronize_onesided_communication
  use rhopotential, only: full_local_potential
  use sparsematrix_base, only: assignment(=), sparsematrix_malloc, sparsematrix_malloc_ptr, SPARSE_FULL, &
                               SPARSE_TASKGROUP
  use constrained_dft, only: cdft_data
  use fragment_base, only: fragmentInputParameters
  use module_fragments, only: system_fragment
  use sparsematrix,only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
  use transposed_operations, only: calculate_overlap_transposed
  use matrix_operations, only: overlapPowerGeneral, check_taylor_order
  use public_enums
  use locreg_operations
  use locregs_init, only: small_to_large_locreg
  !  use Poisson_Solver
  !use allocModule
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  integer,intent(out) :: infoBasisFunctions
  type(atoms_data), intent(in) :: at
  type(orbitals_data) :: orbs
  real(kind=8),dimension(3,at%astruct%nat) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers), intent(inout) :: GPU
  real(kind=8),intent(out) :: trH, fnrm_tmb
  real(kind=8),intent(inout) :: trH_old
  type(DFT_PSP_projectors),intent(inout) :: nlpsp
  integer,intent(in) :: scf_mode
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb
  type(SIC_data) :: SIC !<parameters for the SIC methods
  type(energy_terms),intent(in) :: energs_base
  logical,intent(in) :: do_iterative_orthogonalization
  integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
  real(kind=8),intent(out) :: ratio_deltas
  logical, intent(inout) :: ortho_on
  integer, intent(in) :: extra_states
  integer,intent(in) :: itout
  real(kind=8),intent(in) :: conv_crit, early_stop, gnrm_dynamic, min_gnrm_for_dynamic, kappa_conv
  logical,intent(in) :: experimental_mode
  logical,intent(out) :: can_use_ham
  logical,intent(in) :: correction_co_contra
  type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
  type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
  type(work_transpose),intent(inout) :: wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi
  type(work_mpiaccumulate),intent(inout) :: fnrm, energs_work
  logical, intent(in) :: frag_calc
  !these must all be present together
  type(cdft_data),intent(inout),optional :: cdft
  type(fragmentInputParameters),optional,intent(in) :: input_frag
  type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
 
  ! Local variables
  integer :: iorb, it, it_tot, ncount, ncharge, ii, kappa_satur, nit_exit, ispin, jproc
  !integer :: jorb, nspin
  !real(kind=8),dimension(:),allocatable :: occup_tmp
  real(kind=8) :: meanAlpha, ediff_best, alpha_max, delta_energy, delta_energy_prev, ediff
  real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS
  real(kind=8),dimension(:),allocatable :: hpsit_c_tmp, hpsit_f_tmp, hpsi_tmp, psidiff, tmparr1, tmparr2
  real(kind=8),dimension(:),allocatable :: delta_energy_arr, hpsi_noprecond, kernel_compr_tmp, kernel_best, hphi_nococontra
  logical :: energy_increased, overlap_calculated, energy_diff, energy_increased_previous, complete_reset, even
  logical :: calculate_inverse, allow_increase
  real(kind=8),dimension(:),pointer :: lhphiold, lphiold, hpsit_c, hpsit_f, hpsi_small
  type(energy_terms) :: energs
  real(kind=8), dimension(2):: reducearr
  real(gp) :: econf, dynamic_convcrit, kappa_mean
  real(kind=8) :: energy_first, trH_ref, charge, fnrm_old
  real(kind=8),dimension(3),save :: kappa_history
  integer,save :: nkappa_history
  logical,save :: has_already_converged
  logical,dimension(7) :: exit_loop
  type(matrices) :: ovrlp_old
  integer :: iiorb, ilr, i, ist
  real(kind=8) :: max_error, mean_error
  integer,dimension(1) :: power
  integer,dimension(:),allocatable :: n3p
  interface
     subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
          ldiis, fnrmOldArr, fnrm_old, alpha, trH, trHold, fnrm, alpha_mean, alpha_max, &
          energy_increased, tmb, lhphiold, overlap_calculated, &
          energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
          hpsi_small, experimental_mode, calculate_inverse, correction_co_contra, hpsi_noprecond, &
          norder_taylor, max_inversion_error, precond_convol_workarrays, precond_workarrays,&
          wt_hphi, wt_philarge, wt_hpsinoprecond, &
          cdft, input_frag, ref_frags)
       use module_defs, only: gp,dp,wp
       use module_types
       use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
       use communications_base, only: work_transpose
       use sparsematrix_base, only: matrices
       use constrained_dft, only: cdft_data
       use module_fragments, only: system_fragment,fragmentInputParameters
       implicit none
       integer, intent(in) :: iproc, nproc, it
       integer,intent(inout) :: norder_taylor
       real(kind=8),intent(in) :: max_inversion_error
       type(DFT_wavefunction), target, intent(inout):: tmb
       type(localizedDIISParameters), intent(inout) :: ldiis
       real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: fnrmOldArr
       real(kind=8),intent(inout) :: fnrm_old
       real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
       real(kind=8), intent(out):: trH, alpha_mean, alpha_max
       type(work_mpiaccumulate), intent(inout):: fnrm
       real(kind=8), intent(in):: trHold
       logical,intent(out) :: energy_increased
       real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
       logical, intent(inout):: overlap_calculated
       type(energy_terms), intent(in) :: energs
       real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c) :: hpsit_c
       real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f) :: hpsit_f
       integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
       logical, intent(in) :: experimental_mode, calculate_inverse, correction_co_contra
       real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
       real(kind=8), dimension(tmb%npsidim_orbs),intent(out) :: hpsi_noprecond
       type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
       type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
       type(work_transpose),intent(inout) :: wt_hphi
       type(work_transpose),intent(inout) :: wt_philarge
       type(work_transpose),intent(out) :: wt_hpsinoprecond
       type(cdft_data),intent(inout),optional :: cdft
       type(fragmentInputParameters), optional, intent(in) :: input_frag
       type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
     END SUBROUTINE calculate_energy_and_gradient_linear
  end interface
  interface
     subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, at, do_iterative_orthonormalization, sf_per_type, &
          lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
          experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer,intent(in) :: iproc, nproc, it
       integer,intent(inout) :: order_taylor
       real(kind=8),intent(in) :: max_inversion_error
       type(localizedDIISParameters),intent(inout):: ldiis
       type(DFT_wavefunction),target,intent(inout):: tmb
       type(atoms_data),intent(in) :: at
       logical,intent(in) :: do_iterative_orthonormalization
       integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type 
       real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lphiold
       real(8),intent(in):: trH, meanAlpha, alpha_max
       real(8),dimension(tmb%orbs%norbp),intent(inout):: alpha, alphaDIIS
       real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(inout) :: hpsi_small
       real(kind=8),dimension(tmb%orbs%npsidim_orbs),optional,intent(out) :: psidiff
       logical, intent(in) :: ortho, experimental_mode
       real(kind=8),intent(out) :: trH_ref
       real(kind=8),dimension(tmb%linmat%l%nvctrp_tg*tmb%linmat%l%nspin),intent(inout) :: kernel_best
       logical,intent(out) :: complete_reset
     END SUBROUTINE hpsitopsi_linear
  end interface



  call f_routine(id='getLocalizedBasis')

  !!fnrm = work_mpiaccumulate_null()
  !!fnrm%ncount = 1
  !!call allocate_work_mpiaccumulate(fnrm)

  !!energs_work = work_mpiaccumulate_null()
  !!energs_work%ncount = 4
  !!call allocate_work_mpiaccumulate(energs_work)

  energs = energy_terms_null()
  delta_energy_arr=f_malloc(nit_basis+6,id='delta_energy_arr')
  kernel_best=sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_TASKGROUP,id='kernel_best')
  energy_diff=.false.

  ovrlp_old%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%s, &
                           iaction=SPARSE_TASKGROUP, id='ovrlp_old%matrix_compr')

  ! Allocate all local arrays.
  call allocateLocalArrays()


  call timing(iproc,'getlocbasinit','ON')
  tmb%can_use_transposed=.false.

  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS
  ldiis%resetDIIS=.false.
  ldiis%immediateSwitchToSD=.false.
  allow_increase=.false.
 
  call timing(iproc,'getlocbasinit','OF')

  overlap_calculated=.false.
  it=0
  it_tot=0
  !ortho=.true.
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
  n3p = f_malloc(0.to.nproc-1,id='n3p')
  do jproc=0,nproc-1
      n3p(jproc) = max(denspot%dpbox%nscatterarr(jproc,2),1)
  end do
  call start_onesided_communication(iproc, nproc, denspot%dpbox%ndims(1), denspot%dpbox%ndims(2), &
       n3p, denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf*tmb%ham_descr%comgp%nspin, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, &
       tmb%ham_descr%lzd)
  call f_free(n3p)

  delta_energy_prev=1.d100

  energy_increased_previous=.false.
  ratio_deltas=1.d0
  ediff_best=1.d0
  ediff=1.d0
  delta_energy_prev=1.d0
  delta_energy_arr=1.d0
  trH_ref=trH_old
  dynamic_convcrit=1.d-100
  kappa_satur=0
  fnrm_old = 0.d0


  ! Count whether there is an even or an odd number of electrons
  charge=0.d0
  do iorb=1,orbs%norb
      charge=charge+orbs%occup(iorb)
  end do
  ncharge=nint(charge)
  even=(mod(ncharge,2)==0)

  !!! Purify the initial kernel (only when necessary and if there is an even number of electrons)
  !!if (target_function/=TARGET_FUNCTION_IS_TRACE .and. even .and. &
  !!    (scf_mode==LINEAR_FOE .or. scf_mode==LINEAR_PEXSI)) then
  !!    if (iproc==0) then
  !!        call yaml_sequence(advance='no')
  !!        call yaml_mapping_open(flow=.true.)
  !!        call yaml_map('Initial kernel purification',.true.)
  !!    end if
  !!    overlap_calculated=.true.
  !!    do ispin=1,tmb%linmat%l%nspin
  !!        call purify_kernel(iproc, nproc, tmb, overlap_calculated, 1, 30, order_taylor, &
  !!             max_inversion_error, purification_quickreturn, ispin)
  !!    end do
  !!    if (iproc==0) call yaml_mapping_close()
  !!end if

  if (itout==0) then
      nkappa_history=0
      kappa_history=0.d0
      has_already_converged=.false.
  end if

  iterLoop: do


      it=it+1
      it=max(it,1) !since it could become negative (2 is subtracted if the loop cycles)
      it_tot=it_tot+1

      fnrm%sendbuf(1)=0.d0
      fnrm%receivebuf(1)=0.d0
  
      if (iproc==0) then
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_comment('iter:'//yaml_toa(it,fmt='(i6)'),hfill='-')
          if (target_function==TARGET_FUNCTION_IS_TRACE) then
              call yaml_map('target function','TRACE')
          else if (target_function==TARGET_FUNCTION_IS_ENERGY) then
              call yaml_map('target function','ENERGY')
          else if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              call yaml_map('target function','HYBRID')
          end if
      end if

      ! Synchronize the mpi_get before starting a new communication
      call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)

      ! Start the communication
      call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
           TRANSPOSE_POST, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd, wt_phi)

      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if (tmb%ham_descr%npsidim_orbs > 0)  call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%psi, tmb%ham_descr%psi)

      ! Start the nonblocking transposition (the results will be gathered in
      ! orthoconstraintNonorthogonal)
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           TRANSPOSE_POST, tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd, &
           wt_philarge)

      call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,tmb%hpsi,energs%eproj,tmb%paw)
      ! only kinetic because waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
           & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,3,denspot%xc,&
           & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
           & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
           & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
      ! only potential
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call vcopy(tmb%ham_descr%npsidim_orbs, tmb%hpsi(1), 1, hpsi_tmp(1), 1)
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
               & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
               & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
               & potential=denspot%rhov,comgp=tmb%ham_descr%comgp,&
               hpsi_noconf=hpsi_tmp,econf=econf)

          !!if (nproc>1) then
          !!    call mpiallred(econf, 1, mpi_sum, bigdft_mpi%mpi_comm)
          !!end if

      else
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,&
               & denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
               & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
               & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      end if


      !!if (target_function==TARGET_FUNCTION_IS_HYBRID .and. iproc==0) then
      !!    write(*,*) 'econf, econf/tmb%orbs%norb',econf, econf/tmb%orbs%norb
      !!end if

      call timing(iproc,'glsynchham2','ON')
      call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,denspot%xc,tmb%hpsi,&
           energs,energs_work)
      call timing(iproc,'glsynchham2','OF')

      if (iproc==0) then
          call yaml_map('Hamiltonian Applied',.true.)
      end if

      !if (iproc==0) write(*,'(a,5es16.6)') 'ekin, eh, epot, eproj, eex', &
      !              energs%ekin, energs%eh, energs%epot, energs%eproj, energs%exc


      ! Start the communication
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_POST, hpsi_tmp, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
      else
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_POST, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
      end if

      ! Gather the data
      call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
           TRANSPOSE_GATHER, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd, wt_phi)

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

      ! Use this subroutine to write the energies, with some fake number
      ! to prevent it from writing too much
      if (iproc==0) then
          call write_energies(0,energs,0.d0,0.d0,'',only_energies=.true.)
      end if
      if (iproc==0) then
          call yaml_map('Orthoconstraint',.true.)
      end if

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          tmb%ham_descr%can_use_transposed=.false.
          do ispin=1,tmb%linmat%s%nspin
              call vcopy(tmb%linmat%s%nvctrp_tg, &
                   tmb%linmat%ovrlp_%matrix_compr((ispin-1)*tmb%linmat%s%isvctrp_tg+1), 1, &
                   ovrlp_old%matrix_compr((ispin-1)*tmb%linmat%s%nvctrp_tg+1), 1)
          end do
          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
               tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
          !if (iproc==0) call yaml_newline()
          !if (iproc==0) call yaml_sequence_open('kernel update by renormalization')
          if (it==1 .or. energy_increased .or. .not.experimental_mode) then
              ! Calculate S^1/2, as it can not be taken from memory
              power(1)=2
              call overlapPowerGeneral(iproc, nproc,bigdft_mpi%mpi_comm,&
                   order_taylor, 1, power, -1, &
                   imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
                   ovrlp_mat=ovrlp_old, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(1), &
                   verbosity=0, &
                   check_accur=.true., max_error=max_error, mean_error=mean_error)
              call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
          end if
          call renormalize_kernel(iproc, nproc, order_taylor, max_inversion_error, tmb, tmb%linmat%ovrlp_, ovrlp_old)
          !if (iproc==0) call yaml_sequence_close()
      else
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_GATHER, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
      end if

      ! Gather the data in case it has not been done before
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_GATHER, hpsi_tmp, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
      end if
      ncount=tmb%ham_descr%collcom%ndimind_c
      if(ncount>0) call vcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*tmb%ham_descr%collcom%ndimind_f
      if(ncount>0) call vcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)

      ! optimize the tmbs for a few extra states
      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          kernel_compr_tmp = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_TASKGROUP, id='kernel_compr_tmp')
          !call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
          call vcopy(tmb%linmat%l%nvctrp_tg*tmb%linmat%l%nspin, &
               tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !allocate(occup_tmp(tmb%orbs%norb), stat=istat)
          !call memocc(istat, occup_tmp, 'occup_tmp', subname)
          !call vcopy(tmb%orbs%norb, tmb%orbs%occup(1), 1, occup_tmp(1), 1)
          !call f_zero(tmb%orbs%norb,tmb%orbs%occup(1))
          !call vcopy(orbs%norb, orbs%occup(1), 1, tmb%orbs%occup(1), 1)
          !! occupy the next few states - don't need to preserve the charge as only using for support function optimization
          !do iorb=1,tmb%orbs%norb
          !   if (tmb%orbs%occup(iorb)==1.0_gp) then
          !      tmb%orbs%occup(iorb)=2.0_gp
          !   else if (tmb%orbs%occup(iorb)==0.0_gp) then
          !      do jorb=iorb,min(iorb+extra_states-1,tmb%orbs%norb)
          !         tmb%orbs%occup(jorb)=2.0_gp
          !      end do
          !      exit
          !   end if
          !end do
          call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, &
               tmb%linmat%l, tmb%linmat%kernel_)
          !call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
          !call transform_sparse_matrix(tmb%linmat%denskern, tmb%linmat%denskern_large, 'large_to_small')
      end if

      ! use hpsi_tmp as temporary array for hpsi_noprecond, even if it is allocated with a larger size
      !write(*,*) 'calling calc_energy_and.., correction_co_contra',correction_co_contra
      calculate_inverse = (target_function/=TARGET_FUNCTION_IS_HYBRID)
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call calculate_energy_and_gradient_linear(iproc, nproc, it, ldiis, fnrmOldArr, &
           fnrm_old, alpha, trH, trH_old, fnrm, &
           meanAlpha, alpha_max, energy_increased, tmb, lhphiold, overlap_calculated, energs_base, &
           hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, hpsi_small, &
           experimental_mode, calculate_inverse, &
           correction_co_contra, hpsi_noprecond=hpsi_tmp, norder_taylor=order_taylor, &
           max_inversion_error=max_inversion_error, &
           precond_convol_workarrays=precond_convol_workarrays, precond_workarrays=precond_workarrays, &
           wt_hphi=wt_hphi, wt_philarge=wt_philarge, wt_hpsinoprecond=wt_hpsinoprecond,&
           cdft=cdft, input_frag=input_frag, ref_frags=ref_frags)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      !fnrm_old=fnrm


      if (experimental_mode) then
          if (it_tot==1) then
              energy_first=trH
          end if
          if (iproc==0) call yaml_map('rel D',(trH-energy_first)/energy_first,fmt='(es9.2)')
          if ((trH-energy_first)<0.d0 .and. abs((trH-energy_first)/energy_first)>early_stop .and. itout>0) then
              energy_diff=.true.
          end if
      end if

      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          if (tmb%linmat%l%nspin>1) stop 'THIS IS NOT TESTED FOR SPIN POLARIZED SYSTEMS!'
          call vcopy(tmb%linmat%l%nvctrp_tg*tmb%linmat%l%nspin, kernel_compr_tmp(1), 1, &
               tmb%linmat%kernel_%matrix_compr(1), 1)
          call f_free(kernel_compr_tmp)
      end if

      ediff=trH-trH_old
      ediff_best=trH-trH_ref

      if (it>1 .and. (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode)) then
          if (.not.energy_increased .and. .not.energy_increased_previous) then
              if (.not.ldiis%switchSD) then
                  ratio_deltas=ediff_best/delta_energy_prev
              else
                  ratio_deltas=ediff_best/delta_energy_arr(ldiis%itBest)
              end if
          else
              ! use a default value
              if (iproc==0) then
                  call yaml_warning('use a fake value for kappa')
                  call yaml_newline()
              end if
              ratio_deltas=0.5d0
          end if
          if (ldiis%switchSD) then
              !!ratio_deltas=0.5d0
              !!if (iproc==0) write(*,*) 'WARNING: TEMPORARY FIX for ratio_deltas!'
          end if
          if (iproc==0) call yaml_map('kappa',ratio_deltas,fmt='(es10.3)')
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              !if (ratio_deltas>0.d0) then
              !if (ratio_deltas>1.d-12) then
              if (.not.energy_increased .and. .not.energy_increased_previous) then
                  if (iproc==0) call yaml_map('kappa to history',.true.)
                  nkappa_history=nkappa_history+1
                  ii=mod(nkappa_history-1,3)+1
                  kappa_history(ii)=ratio_deltas
              else
                  if (iproc==0) call yaml_map('kappa to history',.false.)
              end if
              !!if (nkappa_history>=3) then
              !!    kappa_mean=sum(kappa_history)/3.d0
              !!    if (iproc==0) call yaml_map('mean kappa',kappa_mean,fmt='(es10.3)')
              !!    dynamic_convcrit=conv_crit/kappa_mean
              !!    if (iproc==0) call yaml_map('dynamic conv crit',dynamic_convcrit,fmt='(es9.2)')
              !!end if
          end if
      end if
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          if (nkappa_history>=3) then
              kappa_mean=sum(kappa_history)/3.d0
              if (iproc==0) call yaml_map('mean kappa',kappa_mean,fmt='(es10.3)')
              !dynamic_convcrit=conv_crit/kappa_mean
              dynamic_convcrit=gnrm_dynamic/kappa_mean
              if (iproc==0) call yaml_map('dynamic conv crit',dynamic_convcrit,fmt='(es9.2)')
          end if
      end if

      if (energy_increased) then
          energy_increased_previous=.true.
      else
          energy_increased_previous=.false.
      end if



      !!delta_energy_prev=delta_energy

      ! Wait for the communication of fnrm on root
      if (nproc>1) then
          call mpi_fenceandfree(fnrm%window)
      end if
      fnrm%receivebuf(1)=sqrt(fnrm%receivebuf(1)/dble(tmb%orbs%norb))

      ! The other processes need to get fnrm as well. The fence will be later as only iproc=0 has to write.
      if (nproc>1) then
          if (iproc==0) fnrm%sendbuf(1) = fnrm%receivebuf(1)
          fnrm%window = mpiwindow(1, fnrm%sendbuf(1), bigdft_mpi%mpi_comm)
          if (iproc/=0) then
              call mpiget(fnrm%receivebuf(1), 1, 0, int(0,kind=mpi_address_kind), fnrm%window)
          end if
      end if

      if (energy_increased .and. ldiis%isx==0 .and. (.not. allow_increase)) then
          !if (iproc==0) write(*,*) 'WARNING: ENERGY INCREASED'
          !if (iproc==0) call yaml_warning('The target function increased, D='&
          !              //trim(adjustl(yaml_toa(trH-ldiis%trmin,fmt='(es10.3)'))))
          if (nproc>1) then
              call mpi_fenceandfree(fnrm%window)
          end if
          fnrm_old=fnrm%receivebuf(1)
          if (iproc==0) then
              call yaml_newline()
              call yaml_map('iter',it,fmt='(i5)')
              call yaml_map('fnrm',fnrm%receivebuf(1),fmt='(es9.2)')
              call yaml_map('Omega',trH,fmt='(es22.15)')
              call yaml_map('D',ediff,fmt='(es9.2)')
              call yaml_map('D best',ediff_best,fmt='(es9.2)')
          end if
          tmb%ham_descr%can_use_transposed=.false.
          call vcopy(tmb%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
          can_use_ham=.false.
          call vcopy(tmb%linmat%l%nvctrp_tg*tmb%linmat%l%nspin, kernel_best(1), 1, &
               tmb%linmat%kernel_%matrix_compr(1), 1)
          trH_old=0.d0
          it=it-2 !go back one iteration (minus 2 since the counter was increased)
          overlap_calculated=.false.
          ! print info here anyway for debugging
          if (it_tot<2*nit_basis) then ! just in case the step size is the problem
              call yaml_mapping_close()
              call yaml_flush_document()
              !call bigdft_utils_flush(unit=6)
              ! This is to avoid memory leaks
              call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   TRANSPOSE_GATHER, hpsit_c, hpsit_f, hpsi_tmp, tmb%ham_descr%lzd, wt_hpsinoprecond)

              ! for fragment calculations tmbs may be far from orthonormality so allow an increase in energy for a few iterations
              if (frag_calc .and. itout<4) then
                 allow_increase=.true.
              end if
              cycle
          else if(it_tot<3*nit_basis) then ! stop orthonormalizing the tmbs
             if (iproc==0) call yaml_newline()
             if (iproc==0) call yaml_warning('Energy increasing, switching off orthonormalization of tmbs')
             ortho_on=.false.
             alpha=alpha*5.0d0/3.0d0 ! increase alpha to make up for decrease from previous iteration
          end if
      else
          can_use_ham=.true.
      end if 


      ! information on the progress of the optimization
      if (iproc==0) then
          call yaml_newline()
          call yaml_map('iter',it,fmt='(i5)')
          call yaml_map('fnrm',fnrm%receivebuf(1),fmt='(es9.2)')
          call yaml_map('Omega',trH,fmt='(es22.15)')
          call yaml_map('D',ediff,fmt='(es9.2)')
          call yaml_map('D best',ediff_best,fmt='(es9.2)')
      end if

      ! Add some extra iterations if DIIS failed (max 6 failures are allowed before switching to SD)
      nit_exit=min(nit_basis+ldiis%icountDIISFailureTot,nit_basis+6)

      ! Normal case
      if (.not.energy_increased .or. ldiis%isx/=0) then
          if (nproc>1) then
              call mpi_fenceandfree(fnrm%window)
          end if
          fnrm_old=fnrm%receivebuf(1)
      end if

      ! Determine whether the loop should be exited
      exit_loop(1) = (it>=nit_exit)
      exit_loop(2) = (it_tot>=3*nit_basis)
      exit_loop(3) = energy_diff
      exit_loop(4) = (fnrm%receivebuf(1)<conv_crit .and. experimental_mode)
      exit_loop(5) = (experimental_mode .and. fnrm%receivebuf(1)<dynamic_convcrit .and. fnrm%receivebuf(1)<min_gnrm_for_dynamic &
                     .and. (it>1 .or. has_already_converged)) ! first overall convergence not allowed in a first iteration
      exit_loop(6) = (itout==0 .and. it>1 .and. ratio_deltas<kappa_conv .and.  ratio_deltas>0.d0)
      if (ratio_deltas>0.d0 .and. ratio_deltas<1.d-1) then
          kappa_satur=kappa_satur+1
      else
          kappa_satur=0
      end if
      exit_loop(7) = (.false. .and. itout>0 .and. kappa_satur>=2)

      if(any(exit_loop)) then
          if(exit_loop(1)) then
              infoBasisFunctions=-1
              if(iproc==0) call yaml_map('exit criterion','net number of iterations')
          end if
          if (exit_loop(2)) then
              infoBasisFunctions=-2
              if (iproc==0) call yaml_map('exit criterion','total number of iterations')
          end if
          if (exit_loop(3)) then
              infoBasisFunctions=it
              if (iproc==0) call yaml_map('exit criterion','energy difference')
          end if
          if (exit_loop(4)) then
              if (iproc==0) call yaml_map('exit criterion','gradient')
              infoBasisFunctions=it
          end if
          if (exit_loop(5)) then
              if (iproc==0) call yaml_map('exit criterion','dynamic gradient')
              infoBasisFunctions=it
              has_already_converged=.true.
          end if
          if (exit_loop(6)) then
              infoBasisFunctions=it
              if (iproc==0) call yaml_map('exit criterion','extended input guess')
          end if
          if (exit_loop(7)) then
              infoBasisFunctions=it
              if (iproc==0) call yaml_map('exit criterion','kappa')
          end if
          if (can_use_ham) then
              ! Calculate the Hamiltonian matrix, since we have all quantities ready. This matrix can then be used in the first
              ! iteration of get_coeff.
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
                   tmb%ham_descr%psit_c, hpsit_c_tmp, tmb%ham_descr%psit_f, hpsit_f_tmp, tmb%linmat%m, tmb%linmat%ham_)
              !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
          end if

          if (iproc==0) then
              !yaml output
              call yaml_mapping_close() !iteration
              call yaml_flush_document()
              !call bigdft_utils_flush(unit=6)
          end if

          ! This is to avoid memory leaks
          ! Gather together the data (was posted in orthoconstraintNonorthogonal)
          ! Give hpsit_c and hpsit_f, this should not matter if GATHER is specified.
          ! To be modified later
          call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_GATHER, hpsit_c, hpsit_f, hpsi_tmp, tmb%ham_descr%lzd, wt_hpsinoprecond)

          exit iterLoop
      end if
      trH_old=trH

      if (ldiis%isx>0) then
          ldiis%mis=mod(ldiis%is,ldiis%isx)+1 !to store the energy at the correct location in the history
      end if
      call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, at, do_iterative_orthogonalization, sf_per_type, &
           lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho_on, psidiff, &
           experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)


      overlap_calculated=.false.
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmb%ham_descr%can_use_transposed) then
          tmb%ham_descr%can_use_transposed=.false.
      end if

      ! Gather together the data (was posted in orthoconstraintNonorthogonal)
      ! Give hpsit_c and hpsit_f, this should not matter if GATHER is specified.
      ! To be modified later
      hphi_nococontra = f_malloc(tmb%ham_descr%npsidim_orbs,id='hphi_nococontra')
      call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           TRANSPOSE_GATHER, hpsit_c, hpsit_f, hphi_nococontra, tmb%ham_descr%lzd, wt_hpsinoprecond)
      call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, hphi_nococontra, hpsi_tmp)
      call f_free(hphi_nococontra)


      ! Estimate the energy change, that is to be expected in the next optimization
      ! step, given by the product of the force and the "displacement" .
      if (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode) then
          call estimate_energy_change(tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%l%nspin, psidiff, &
               hpsi_tmp,delta_energy)
          ! This is a hack...
          if (energy_increased) then
              delta_energy=1.d100
              !ratio_deltas=1.d100
          end if
          !if (iproc==0) write(*,*) 'delta_energy', delta_energy
          delta_energy_prev=delta_energy
          delta_energy_arr(max(it,1))=delta_energy !max since the counter was decreased if there are problems, might lead to wrong results otherwise
      end if


      ! Only need to reconstruct the kernel if it is actually used.
      if ((target_function/=TARGET_FUNCTION_IS_TRACE .or. scf_mode==LINEAR_DIRECT_MINIMIZATION) &
           .and. .not.complete_reset ) then
          if(scf_mode/=LINEAR_FOE .and. scf_mode/=LINEAR_PEXSI) then
              call reconstruct_kernel(iproc, nproc, order_taylor, tmb%orthpar%blocksize_pdsyev, &
                   tmb%orthpar%blocksize_pdgemm, orbs, tmb, overlap_calculated)
              if (iproc==0) call yaml_map('reconstruct kernel',.true.)
          else if (experimental_mode .and. .not.complete_reset) then
          end if
      end if

      if (iproc==0) then
          call yaml_mapping_close() !iteration
          call yaml_flush_document()
          !call bigdft_utils_flush(unit=6)
      end if


  end do iterLoop

  ! Write the final results
  if (iproc==0) then
      call yaml_sequence(label='final_supfun'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))),advance='no')
      call yaml_mapping_open(flow=.true.)
      call yaml_comment('iter:'//yaml_toa(it,fmt='(i6)'),hfill='-')
      if (target_function==TARGET_FUNCTION_IS_TRACE) then
          call yaml_map('target function','TRACE')
      else if (target_function==TARGET_FUNCTION_IS_ENERGY) then
          call yaml_map('target function','ENERGY')
      else if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call yaml_map('target function','HYBRID')
      end if
      call write_energies(0,energs,0.d0,0.d0,'',only_energies=.true.)
      call yaml_newline()
      call yaml_map('iter',it,fmt='(i5)')
      call yaml_map('fnrm',fnrm%receivebuf(1),fmt='(es9.2)')
      call yaml_map('Omega',trH,fmt='(es22.15)')
      call yaml_map('D',ediff,fmt='(es9.2)')
      call yaml_map('D best',ediff_best,fmt='(es9.2)')
      call yaml_mapping_close() !iteration
      call yaml_flush_document()
      !call bigdft_utils_flush(unit=6)
  end if


  if (iproc==0) then
      call yaml_comment('Support functions created')
  end if


  ! Deallocate potential
  call f_free_ptr(denspot%pot_work)


  ! Keep the values for the next iteration
  reducearr(1)=0.d0
  reducearr(2)=0.d0
  do iorb=1,tmb%orbs%norbp
      reducearr(1)=reducearr(1)+alpha(iorb)
      reducearr(2)=reducearr(2)+alphaDIIS(iorb)
  end do

  if (nproc > 1) then
      call mpiallred(reducearr, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

  reducearr(1)=reducearr(1)/dble(tmb%orbs%norb)
  reducearr(2)=reducearr(2)/dble(tmb%orbs%norb)

  ldiis%alphaSD=reducearr(1)
  ldiis%alphaDIIS=reducearr(2)


  ! Deallocate all local arrays.
  call deallocateLocalArrays()
  call f_free(delta_energy_arr)
  call f_free(kernel_best)
  call f_free_ptr(ovrlp_old%matrix_compr)

  fnrm_tmb = fnrm%receivebuf(1)
  !!call deallocate_work_mpiaccumulate(fnrm)
  !!call deallocate_work_mpiaccumulate(energs_work)

  call f_release_routine()

contains


    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
    logical :: with_confpot
    integer :: iiorb, ilr, ncplx
    real(gp) :: kx, ky, kz

      alpha = f_malloc(tmb%orbs%norbp,id='alpha')
      alphaDIIS = f_malloc(tmb%orbs%norbp,id='alphaDIIS')
      fnrmOldArr = f_malloc(tmb%orbs%norbp,id='fnrmOldArr')
      hpsi_small = f_malloc_ptr(max(tmb%npsidim_orbs, tmb%npsidim_comp),id='hpsi_small')
      lhphiold = f_malloc_ptr(max(tmb%npsidim_orbs, tmb%npsidim_comp),id='lhphiold')
      lphiold = f_malloc_ptr(size(tmb%psi),id='lphiold')
      hpsit_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
      hpsit_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
      hpsit_c_tmp = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c_tmp')
      hpsit_f_tmp = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f_tmp')
      hpsi_tmp = f_malloc(tmb%ham_descr%npsidim_orbs,id='hpsi_tmp')
      psidiff = f_malloc(tmb%npsidim_orbs,id='psidiff')
      !hpsi_noprecond = f_malloc(tmb%npsidim_orbs,id='hpsi_noprecond')


      !!allocate(precond_convol_workarrays(tmb%orbs%norbp))
      !!allocate(precond_workarrays(tmb%orbs%norbp))
      !!do iorb=1,tmb%orbs%norbp
      !!    iiorb=tmb%orbs%isorb+iorb
      !!    ilr=tmb%orbs%inwhichlocreg(iiorb)
      !!    with_confpot = (tmb%confdatarr(iorb)%prefac/=0.d0)
      !!    call init_local_work_arrays(tmb%lzd%llr(ilr)%d%n1, tmb%lzd%llr(ilr)%d%n2, tmb%lzd%llr(ilr)%d%n3, &
      !!         tmb%lzd%llr(ilr)%d%nfl1, tmb%lzd%llr(ilr)%d%nfu1, &
      !!         tmb%lzd%llr(ilr)%d%nfl2, tmb%lzd%llr(ilr)%d%nfu2, &
      !!         tmb%lzd%llr(ilr)%d%nfl3, tmb%lzd%llr(ilr)%d%nfu3, &
      !!         with_confpot, precond_convol_workarrays(iorb))
      !!    kx=tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb))
      !!    ky=tmb%orbs%kpts(2,tmb%orbs%iokpt(iorb))
      !!    kz=tmb%orbs%kpts(3,tmb%orbs%iokpt(iorb))
      !!    if (kx**2+ky**2+kz**2 > 0.0_gp .or. tmb%orbs%nspinor==2 ) then
      !!       ncplx=2
      !!    else
      !!       ncplx=1
      !!    end if
      !!    call allocate_work_arrays(tmb%lzd%llr(ilr)%geocode, tmb%lzd%llr(ilr)%hybrid_on, &
      !!         ncplx, tmb%lzd%llr(ilr)%d, precond_workarrays(iorb))
      !!end do


    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !
    integer :: iiorb, ilr, ncplx
    real(gp) :: kx, ky, kz

    call f_free(alpha)
    call f_free(alphaDIIS)
    call f_free(fnrmOldArr)
    call f_free_ptr(hpsi_small)
    call f_free_ptr(lhphiold)
    call f_free_ptr(lphiold)
    call f_free_ptr(hpsit_c)
    call f_free_ptr(hpsit_f)
    call f_free(hpsit_c_tmp)
    call f_free(hpsit_f_tmp)
    call f_free(hpsi_tmp)
    call f_free(psidiff)
    !call f_free(hpsi_noprecond)
    !!do iorb=1,tmb%orbs%norbp
    !!    iiorb=tmb%orbs%isorb+iorb
    !!    ilr=tmb%orbs%inwhichlocreg(iiorb)
    !!    call deallocate_workarrays_quartic_convolutions(precond_convol_workarrays(iorb))
    !!    kx=tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb))
    !!    ky=tmb%orbs%kpts(2,tmb%orbs%iokpt(iorb))
    !!    kz=tmb%orbs%kpts(3,tmb%orbs%iokpt(iorb))
    !!    if (kx**2+ky**2+kz**2 > 0.0_gp .or. tmb%orbs%nspinor==2 ) then
    !!       ncplx=2
    !!    else
    !!       ncplx=1
    !!    end if
    !!    call deallocate_work_arrays(tmb%lzd%llr(ilr)%geocode, tmb%lzd%llr(ilr)%hybrid_on, &
    !!         ncplx, precond_workarrays(iorb))
    !!end do
    !!deallocate(precond_convol_workarrays)
    !!deallocate(precond_workarrays)

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis



subroutine improveOrbitals(iproc, nproc, tmb, nspin, ldiis, alpha, gradient, experimental_mode)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspin
  type(DFT_wavefunction),intent(inout) :: tmb
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmb%orbs%norbp),intent(in) :: alpha
  real(kind=wp),dimension(tmb%npsidim_orbs),intent(inout) :: gradient
  logical,intent(in) :: experimental_mode
  
  ! Local variables
  integer :: istart, iorb, iiorb, ilr, ncount

  call f_routine(id='improveOrbitals')

  if(ldiis%isx==0) then ! steepest descents
      call timing(iproc,'optimize_SD   ','ON')
      istart=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          call daxpy(ncount, -alpha(iorb), gradient(istart), 1, tmb%psi(istart), 1)
          istart=istart+ncount
      end do
      call timing(iproc,'optimize_SD   ','OF')
  else! DIIS
      ldiis%mis=mod(ldiis%is,ldiis%isx)+1
      ldiis%is=ldiis%is+1
      if(ldiis%alphaDIIS/=1.d0) then
          if (tmb%orbs%norbp>0) call dscal(tmb%npsidim_orbs, ldiis%alphaDIIS, gradient, 1)
      end if
      call optimizeDIIS(iproc, nproc, max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, nspin, tmb%lzd, gradient, tmb%psi, ldiis, &
           experimental_mode)
  end if

  call f_release_routine()

end subroutine improveOrbitals



subroutine diagonalizeHamiltonian2(iproc, norb, HamSmall, ovrlp, eval)
  !
  ! Purpose:
  ! ========
  !   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
  !   the same result. This is done by requiring that the first entry of each vector
  !   is positive.
  !
  ! Calling arguments:
  ! ==================
  !   Input arguments:
  !   ----------------
  !     iproc     process ID
  !     nproc     number of MPI processes
  !     orbs      type describing the physical orbitals psi
  !   Input / Putput arguments
  !     HamSmall  on input: the Hamiltonian
  !               on exit: the eigenvectors
  !   Output arguments
  !     eval      the associated eigenvalues 
  !
  use module_base
  use module_types
  use yaml_output, only: yaml_map
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, norb
  real(kind=8),dimension(norb, norb),intent(inout) :: HamSmall
  real(kind=8),dimension(norb, norb),intent(inout) :: ovrlp
  real(kind=8),dimension(norb),intent(out) :: eval

  ! Local variables
  integer :: lwork, info
  real(kind=8),dimension(:),allocatable :: work
  character(len=*),parameter :: subname='diagonalizeHamiltonian'
  !!real(8),dimension(:,:),pointer :: hamtmp, ovrlptmp, invovrlp, tmpmat, tmpmat2
  !!real(8) :: tt, tt2
  !!integer :: nproc
  !!real(8),dimension(norb,norb) :: kernel

  !!allocate(hamtmp(norb,norb))
  !!allocate(ovrlptmp(norb,norb))
  !!allocate(invovrlp(norb,norb))
  !!allocate(tmpmat(norb,norb))
  !!allocate(tmpmat2(norb,norb))

  !!call mpi_comm_size(mpi_comm_world,nproc,istat)

  !!hamtmp=HamSmall
  !!ovrlptmp=ovrlp
  !!call overlapPowerGeneral(iproc, nproc, 100, -2, -1, norb, ovrlptmp, invovrlp, tt)

  !!call dgemm('n', 'n', norb, norb, norb, 1.d0, invovrlp, norb, hamtmp, norb, 0.d0, tmpmat, norb)
  !!call dgemm('n', 'n', norb, norb, norb, 1.d0, tmpmat, norb, invovrlp, norb, 0.d0, tmpmat2, norb)

  !!lwork=10000
  !!allocate(work(lwork))
  !!call dsyev('v', 'l', norb, tmpmat2, norb, eval, work, lwork, info)
  !!deallocate(work)

  !!ovrlptmp=ovrlp
  !!tmpmat=tmpmat2
  !!call overlapPowerGeneral(iproc, nproc, 100, -2, -1, norb, ovrlptmp, invovrlp, tt)
  !!!call dgemm('n', 'n', norb, norb, norb, 1.d0, invovrlp, norb, tmpmat, norb, 0.d0, tmpmat2, norb)
  !!!if (iproc==0) then
  !!!    do istat=1,norb
  !!!        do iall=1,norb
  !!!            write(200,*) tmpmat2(iall,istat)
  !!!        end do
  !!!    end do
  !!!end if

  !!call dgemm('n', 't', norb, norb, 28, 1.d0, tmpmat2, norb, tmpmat2, norb, 0.d0, kernel, norb)
  !!if (iproc==0) then
  !!    tt=0.d0
  !!    tt2=0.d0
  !!    do istat=1,norb
  !!        do iall=1,norb
  !!            write(300,*) kernel(iall,istat)
  !!            if (istat==iall) tt=tt+kernel(iall,istat)
  !!            tt2=tt2+kernel(iall,istat)*ovrlp(iall,istat)
  !!        end do
  !!    end do
  !!    write(*,*) 'Before: trace(K)',tt
  !!    write(*,*) 'Before: trace(KS)',tt2
  !!end if

  !!call dgemm('n', 'n', norb, norb, norb, 1.d0, invovrlp, norb, kernel, norb, 0.d0, tmpmat, norb)
  !!call dgemm('n', 'n', norb, norb, norb, 1.d0, tmpmat, norb, invovrlp, norb, 0.d0, kernel, norb)
  !!if (iproc==0) then
  !!    tt=0.d0
  !!    tt2=0.d0
  !!    do istat=1,norb
  !!        do iall=1,norb
  !!            write(305,*) kernel(iall,istat)
  !!            if (istat==iall) tt=tt+kernel(iall,istat)
  !!            tt2=tt2+kernel(iall,istat)*ovrlp(iall,istat)
  !!        end do
  !!    end do
  !!    write(*,*) 'After: trace(K)',tt
  !!    write(*,*) 'After: trace(KS)',tt2
  !!end if


  call timing(iproc,'diagonal_seq  ','ON')
  call f_routine(id='diagonalizeHamiltonian2')

  ! DEBUG: print hamiltonian and overlap matrices
  !if (iproc==0) then
  !   open(10)
  !   open(11)
  !   do iorb=1,orbs%norb
  !      do jorb=1,orbs%norb
  !         write(10,*) iorb,jorb,HamSmall(iorb,jorb)
  !         write(11,*) iorb,jorb,ovrlp(iorb,jorb)
  !      end do
  !      write(10,*) ''
  !      write(11,*) ''
  !   end do
  !   close(10)
  !   close(11)
  !end if
  ! DEBUG: print hamiltonian and overlap matrices

  !call yaml_map('Hamiltonian before',HamSmall)
  ! Get the optimal work array size
  lwork=-1 
  work = f_malloc(100,id='work')
  call dsygv(1, 'v', 'l', norb, HamSmall(1,1), norb, ovrlp(1,1), norb, eval(1), work(1), lwork, info) 
  lwork=int(work(1))

  ! Deallocate the work array and reallocate it with the optimal size
  call f_free(work)
  work = f_malloc(lwork,id='work')

  ! Diagonalize the Hamiltonian
  call dsygv(1, 'v', 'l', norb, HamSmall(1,1), norb, ovrlp(1,1), norb, eval(1), work(1), lwork, info) 
  if(info/=0)then
    write(*,*) 'ERROR: dsygv in diagonalizeHamiltonian2, info=',info,'N=',norb
  end if
  !!if (iproc==0) then
  !!    do istat=1,norb
  !!        do iall=1,norb
  !!            write(201,*) hamsmall(iall,istat)
  !!        end do
  !!    end do
  !!end if

  call f_free(work)

  call f_release_routine()
  call timing(iproc,'diagonal_seq  ','OF')

end subroutine diagonalizeHamiltonian2



subroutine large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
       orbs, philarge, phismall)
  use module_base
  use module_types
  use locreg_operations, only: psi_to_locreg2
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, npsidim_orbs_small, npsidim_orbs_large
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(npsidim_orbs_large),intent(in) :: philarge
  real(kind=8),dimension(npsidim_orbs_small),intent(out) :: phismall
  
  ! Local variables
  integer :: istl, ists, ilr, ldim, gdim, iorb
       call timing(iproc,'large2small','ON') ! lr408t   
  ! Transform back to small locreg
  ! No need to this array to zero, since all values will be filled with a value during the copy.
  !!call f_zero(npsidim_orbs_small, phismall(1))
  ists=1
  istl=1
  do iorb=1,orbs%norbp
      ilr = orbs%inWhichLocreg(orbs%isorb+iorb)
      ldim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      gdim=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
      call psi_to_locreg2(iproc, ldim, gdim, lzdsmall%llr(ilr), lzdlarge%llr(ilr), &
           philarge(istl:istl+gdim-1), phismall(ists:ists+ldim-1))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
  end do

  if(orbs%norbp>0 .and. ists/=npsidim_orbs_small+1) stop 'ists/=npsidim_orbs_small+1'
  if(orbs%norbp>0 .and. istl/=npsidim_orbs_large+1) stop 'istl/=npsidim_orbs_large+1'
       call timing(iproc,'large2small','OF') ! lr408t 
end subroutine large_to_small_locreg







subroutine communicate_basis_for_density_collective(iproc, nproc, lzd, npsidim, orbs, lphi, collcom_sr)
  use module_base
  use module_types
  !use module_interfaces, except_this_one => communicate_basis_for_density_collective
  use communications, only: transpose_switch_psir, transpose_communicate_psir, transpose_unswitch_psirt
  use locreg_operations
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(npsidim),intent(in) :: lphi
  type(comms_linear),intent(inout) :: collcom_sr
  
  ! Local variables
  integer :: ist, istr, iorb, iiorb, ilr
  real(kind=8),dimension(:),allocatable :: psir, psirwork, psirtwork
  type(workarr_sumrho) :: w
  character(len=*),parameter :: subname='comm_basis_for_dens_coll'

  call timing(iproc,'commbasis4dens','ON')
  call f_routine(id='communicate_basis_for_density_collective')

  psir = f_malloc(collcom_sr%ndimpsi_c,id='psir')

  ! Allocate the communication buffers for the calculation of the charge density.
  !call allocateCommunicationbufferSumrho(iproc, comsr, subname)
  ! Transform all orbitals to real space.
  ist=1
  istr=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inWhichLocreg(iiorb)
      call initialize_work_arrays_sumrho(lzd%Llr(ilr),.true.,w)
      call daub_to_isf(lzd%Llr(ilr), w, lphi(ist), psir(istr))
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
  end do
  if(istr/=collcom_sr%ndimpsi_c+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
      stop
  end if

  psirwork = f_malloc(collcom_sr%ndimpsi_c,id='psirwork')

  call transpose_switch_psir(collcom_sr, psir, psirwork)

  call f_free(psir)

  psirtwork = f_malloc(collcom_sr%ndimind_c,id='psirtwork')

  call transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)

  call f_free(psirwork)

  call transpose_unswitch_psirt(collcom_sr, psirtwork, collcom_sr%psit_c)

  call f_free(psirtwork)

  call f_release_routine()
  call timing(iproc,'commbasis4dens','OF')

end subroutine communicate_basis_for_density_collective




subroutine DIISorSD(iproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt, trH_ref, kernel_best, complete_reset)
  use module_base
  use module_types
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, it
  real(kind=8),intent(in) :: trH
  type(DFT_wavefunction),intent(inout) :: tmbopt
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmbopt%orbs%norbp),intent(inout) :: alpha, alphaDIIS
  real(kind=8),dimension(max(tmbopt%npsidim_orbs,tmbopt%npsidim_comp)),intent(out):: lphioldopt
  real(kind=8),intent(out) :: trH_ref
  real(kind=8),dimension(tmbopt%linmat%l%nvctrp_tg*tmbopt%linmat%l%nspin),intent(inout) :: kernel_best
  logical,intent(out) :: complete_reset
  
  ! Local variables
  integer :: idsx, ii, offset, istdest, iorb, iiorb, ilr, ncount, istsource
  character(len=2) :: numfail_char
  character(len=10) :: stepsize_char
  

  ! Purpose:
  ! ========
  !   This subroutine decides whether one should use DIIS or variable step size
  !   steepest descent to improve the orbitals. In the beginning we start with DIIS
  !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
  !   steepest descent. If the steepest descent iterations are successful, we switch
  !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
  !   history length is limited to be larger or equal than lin%DIISHistMin.

  ! indicates whether both the support functions and the kernel have been reset
  complete_reset=.false.

  ! history of the energy
  if (ldiis%isx>0) then
      ldiis%energy_hist(ldiis%mis)=trH
  end if
  !!write(*,'(a,10es14.6)') 'ldiis%energy_hist', ldiis%energy_hist

  ! If we swicthed to SD in the previous iteration, reset this flag.
  if(ldiis%switchSD) ldiis%switchSD=.false.
  !if(iproc==0) write(*,'(a,2es15.6,l5)') 'trH, ldiis%trmin, ldiis%resetDIIS', trH, ldiis%trmin, ldiis%resetDIIS

  ! Now come some checks whether the trace is descreasing or not. This further decides
  ! whether we should use DIIS or SD.

  ! Determine wheter the trace is decreasing (as it should) or increasing.
  ! This is done by comparing the current value with diisLIN%energy_min, which is
  ! the minimal value of the trace so far.
  !if(iproc==0) write(*,*) 'trH, ldiis%trmin', trH, ldiis%trmin
  if(trH<=ldiis%trmin+1.d-12*abs(ldiis%trmin) .and. .not.ldiis%resetDIIS) then !1.d-12 is here to tolerate some noise...
      ! Everything ok
      ldiis%trmin=trH
      ldiis%switchSD=.false.
      ldiis%itBest=it
      ldiis%icountSDSatur=ldiis%icountSDSatur+1
      ldiis%icountDIISFailureCons=0
      trH_ref=trH
      call vcopy(tmbopt%linmat%l%nvctrp_tg*tmbopt%linmat%l%nspin, &
           tmbopt%linmat%kernel_%matrix_compr(1), 1, kernel_best(1), 1)
      !if(iproc==0) write(*,*) 'everything ok, copy last psi...'
      call vcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)

      ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
      ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
      if(ldiis%icountSDSatur>=10 .and. ldiis%isx==0 .or. ldiis%immediateSwitchToSD) then
          ldiis%icountSwitch=ldiis%icountSwitch+1
          idsx=max(ldiis%DIISHistMin,ldiis%DIISHistMax-ldiis%icountSwitch)
          if(idsx>0) then
              if(iproc==0) call yaml_map('Switch to DIIS with new history length',idsx)
              !write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
              ldiis%icountSDSatur=0
              ldiis%icountSwitch=0
              ldiis%icountDIISFailureTot=0
              ldiis%icountDIISFailureCons=0
              ldiis%is=0
              ldiis%switchSD=.false.
              ldiis%trmin=1.d100
              ldiis%trold=1.d100
              alpha=ldiis%alphaSD
              alphaDIIS=ldiis%alphaDIIS
              ldiis%icountDIISFailureTot=0
              ldiis%icountDIISFailureCons=0
              ldiis%immediateSwitchToSD=.false.
          end if
      end if
  else
      ! The trace is growing.
      ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
      ! total failures or after 2 consecutive failures.
      if (ldiis%isx>0) then
          ldiis%icountDIISFailureCons=ldiis%icountDIISFailureCons+1
          ldiis%icountDIISFailureTot=ldiis%icountDIISFailureTot+1
      end if
      ldiis%icountSDSatur=0
      if((ldiis%icountDIISFailureCons>=4 .or. ldiis%icountDIISFailureTot>=6 .or. ldiis%resetDIIS) .and. ldiis%isx>0) then
          ! Switch back to SD.
          alpha=ldiis%alphaSD
          if(iproc==0) then
              !if(ldiis%icountDIISFailureCons>=4) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
              !    ldiis%icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
              write(numfail_char,'(i2.2)') ldiis%icountDIISFailureCons
              write(stepsize_char,'(es10.3)') alpha(1)
              if(ldiis%icountDIISFailureCons>=4) then
                  call yaml_warning('DIIS failed '//numfail_char//' times consecutively. &
                       &Switch to SD with stepsize'//stepsize_char//'.')
                  call yaml_newline()
                  !!write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  !!ldiis%icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
              end if
              !!if(ldiis%icountDIISFailureTot>=6) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
              !!    ldiis%icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
              if(ldiis%icountDIISFailureTot>=6) then
                  call yaml_warning('DIIS failed '//numfail_char//' times in total. &
                       &Switch to SD with stepsize'//stepsize_char//'.' )
                  call yaml_newline()
              end if
              if(ldiis%resetDIIS) then
                  call yaml_warning('reset DIIS due to flag')
                  call yaml_newline()
                  !write(*,'(1x,a)') 'reset DIIS due to flag'
              end if
              
          end if
          if(ldiis%resetDIIS) then
              ldiis%resetDIIS=.false.
              ldiis%immediateSwitchToSD=.true.
              ldiis%trmin=1.d100
          end if
          ! Otherwise there could be problems due to the orthonormalization (which sligtly increases 
          ! value of the target function)
          ldiis%trmin=1.d100
          ! Try to get back the orbitals of the best iteration. This is possible if
          ! these orbitals are still present in the DIIS history.
          if(it-ldiis%itBest<ldiis%isx) then
              if(iproc==0) then
                  !!if(iproc==0) write(*,'(1x,a,i0,a)')  'Recover the orbitals from iteration ', &
                  !!    ldiis%itBest, ' which are the best so far.'
                  if (iproc==0) then
                      call yaml_map('Take best TMBs from history',ldiis%itBest)
                  end if
              end if
              ii=modulo(ldiis%mis-(it-ldiis%itBest)-1,ldiis%isx)+1
              !if (iproc==0) write(*,*) 'ii',ii
              offset=0
              istdest=1
              !if(iproc==0) write(*,*) 'copy DIIS history psi...'
              do iorb=1,tmbopt%orbs%norbp
                  iiorb=tmbopt%orbs%isorb+iorb
                  ilr=tmbopt%orbs%inWhichLocreg(iiorb)
                  ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
                  istsource=offset+(ii-1)*ncount+1
                  call vcopy(ncount, ldiis%phiHist(istsource), 1, tmbopt%psi(istdest), 1)
                  call vcopy(ncount, ldiis%phiHist(istsource), 1, lphioldopt(istdest), 1)
                  !if (iproc==0 .and. iorb==1) write(*,*) 'istsource, istdest, val', istsource, istdest, tmbopt%psi(istdest)
                  offset=offset+ldiis%isx*ncount
                  istdest=istdest+ncount
              end do
              trH_ref=ldiis%energy_hist(ii)
              !!if (iproc==0) write(*,*) 'take energy from entry',ii
              call vcopy(tmbopt%linmat%l%nvctrp_tg*tmbopt%linmat%l%nspin, &
                   kernel_best(1), 1, tmbopt%linmat%kernel_%matrix_compr(1), 1)
              !!call vcopy(tmbopt%linmat%l%nvctr, kernel_best(1), 1, tmbopt%linmat%denskern_large%matrix_compr(1), 1)
              complete_reset=.true.
          else
              !if(iproc==0) write(*,*) 'copy last psi...'
              call vcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)
              trH_ref=trH
          end if
          ldiis%isx=0
          ldiis%switchSD=.true.
      end if
      ! to indicate that no orthonormalization is required... (CHECK THIS!)
      if(ldiis%isx==0) ldiis%switchSD=.true. 
  end if

end subroutine DIISorSD

subroutine reconstruct_kernel(iproc, nproc, inversion_method, blocksize_dsyev, blocksize_pdgemm, &
           orbs, tmb, overlap_calculated)
  use module_base
  use module_types
  use module_interfaces, only: calculate_density_kernel
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, uncompress_matrix2
  use transposed_operations, only: calculate_overlap_transposed
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
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
     ! This can then be deleted if the transition to the new type has been completed.
     !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr

     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
  end if

  tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
  !!do i=1,tmb%linmat%s%nspin*tmb%linmat%s%nvctr
  !!    write(2500,'(a,i8,es16.5)') 'i, tmb%linmat%ovrlp_%matrix_compr(i)', i, tmb%linmat%ovrlp_%matrix_compr(i)
  !!end do
  call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
       tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
  !!do ispin=1,tmb%linmat%s%nspin
  !!    do i=1,tmb%linmat%s%nfvctr
  !!        do j=1,tmb%linmat%s%nfvctr
  !!            write(2600,'(a,3i8,es16.5)') 'ispin, i, j, tmb%linmat%ovrlp_%matrix(j,i,ispin)', ispin, i, j, tmb%linmat%ovrlp_%matrix(j,i,ispin)
  !!        end do
  !!    end do
  !!end do
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, &
       tmb%orbs, tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, orbs)

  call f_free_ptr(tmb%linmat%ovrlp_%matrix)


  ! Recalculate the kernel
  call calculate_density_kernel(iproc, nproc, .true., orbs, tmb%orbs, tmb%coeff, tmb%linmat%l, tmb%linmat%kernel_)
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
             !!        write(2200+iproc,'(a,2i9,es13.5)') 'iorb, jorb, ovrlp_coeff(jorb,iorb)', iorb, jorb, ovrlp_coeff(jorb,iorb)
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

      ! should convert this to yaml (LG: easily done)
      if (iproc==0) call yaml_newline()
      if (basis_overlap%nspin==1) then
          if (iproc==0) call yaml_map('Max deviation from unity in reorthonormalize_coeff',max_error,fmt='(es8.2)')
          if (iproc==0) call yaml_map('Mean deviation from unity in reorthonormalize_coeff',mean_error,fmt='(es8.2)')
      else
          if (ispin==1) then
              if (iproc==0) call yaml_map('spin up, Max deviation from unity in reorthonormalize_coeff',max_error,fmt='(es8.2)')
              if (iproc==0) call yaml_map('spin up, Mean deviation from unity in reorthonormalize_coeff',mean_error,fmt='(es8.2)')
          else if (ispin==2) then
              if (iproc==0) call yaml_map('spin down, Max deviation from unity in reorthonormalize_coeff',max_error,fmt='(es8.2)')
              if (iproc==0) call yaml_map('spin down, Mean deviation from unity in reorthonormalize_coeff',mean_error,fmt='(es8.2)')
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


!> Estimate the energy change, given by the product of the force and the "displacement" .
subroutine estimate_energy_change(npsidim_orbs, orbs, lzd, nspin, psidiff, hpsi_noprecond, delta_energy)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer, intent(in) :: npsidim_orbs, nspin
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  real(kind=8),dimension(npsidim_orbs),intent(in) :: psidiff, hpsi_noprecond
  real(kind=8),intent(out) :: delta_energy

  ! Local variables
  integer :: ist, iorb, iiorb, ilr, ncount
  real(kind=8) :: tt, ddot

  call f_routine(id='estimate_energy_change')

  ist=1
  delta_energy=0.d0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      tt=ddot(ncount, psidiff(ist), 1, hpsi_noprecond(ist), 1)
      delta_energy=delta_energy+2.0d0*tt
      ist=ist+ncount
  end do

  if (nspin==1) then
      delta_energy = 2.d0*delta_energy
  end if

  if (bigdft_mpi%nproc > 1) then
      call mpiallred(delta_energy, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if

  call f_release_routine()

end subroutine estimate_energy_change




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
  call allocate_matrices(tmb%linmat%m, allocate_full=.true., &
       matname='gradmat', mat=gradmat)

  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
       hpsit_c, hpsit_c, hpsit_f, hpsit_f, tmb%linmat%m, gradmat)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, gradmat)


  !gradmat%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='gradmat%matrix')
  tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
  call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
       tmb%linmat%m, gradmat%matrix_compr, gradmat%matrix)
  call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
       tmb%linmat%m, tmb%linmat%ham_%matrix_compr, tmb%linmat%ham_%matrix)
  call uncompress_matrix2(iproc, nproc, bigdft_mpi%mpi_comm, &
       tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, tmb%linmat%kernel_%matrix)
  KH=f_malloc0((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr,tmb%linmat%l%nspin/),id='KH')
  KHKH=f_malloc0((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr,tmb%linmat%l%nspin/),id='KHKH')

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
          ii=mod(iiorb-1,tmb%linmat%l%nfvctr)+1
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, 1, tmb%linmat%l%nfvctr, 1.0d0, tmb%linmat%kernel_%matrix(1,1,ispin), &
               tmb%linmat%l%nfvctr, tmb%linmat%ham_%matrix(1,ii,ispin), tmb%linmat%l%nfvctr, &
               0.d0, KH(1,ii,ispin), tmb%linmat%l%nfvctr)
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
          ii=mod(iiorb-1,tmb%linmat%l%nfvctr)+1
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, 1, tmb%linmat%l%nfvctr, scale_factor, KH(1,1,ispin), &
               tmb%linmat%l%nfvctr, KH(1,ii,ispin), tmb%linmat%l%nfvctr, &
               0.d0, KHKH(1,ii,ispin), tmb%linmat%l%nfvctr)
      end do
  end if

  if (nproc > 1) then
      call mpiallred(KHKH, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if
  call f_free(KH)
  Kgrad=f_malloc0((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr,tmb%linmat%l%nspin/),id='Kgrad')
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
          ii=mod(iiorb-1,tmb%linmat%l%nfvctr)+1
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, 1, tmb%linmat%l%nfvctr, 1.0d0, tmb%linmat%kernel_%matrix(1,1,ispin), &
               tmb%linmat%l%nfvctr, gradmat%matrix(1,ii,ispin), tmb%linmat%l%nfvctr, &
               0.d0, Kgrad(1,ii,ispin), tmb%linmat%l%nfvctr)
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
  do ispin=1,tmb%linmat%l%nspin
      do iorb=1,tmb%linmat%l%nfvctr
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



subroutine renormalize_kernel(iproc, nproc, order_taylor, max_inversion_error, tmb, ovrlp, ovrlp_old)
  use module_base
  use module_types
  use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, assignment(=), &
                               SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, &
                               matrices
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
  real(kind=8) :: max_error, mean_error
  type(matrices) :: inv_ovrlp
  real(kind=8),dimension(:,:),pointer :: inv_ovrlpp, tempp
  real(kind=8),dimension(:),allocatable :: inv_ovrlp_compr_seq, kernel_compr_seq
  integer,dimension(3) :: power
  !!real(8) :: tr
  !!integer :: ind, iorb


  call f_routine(id='renormalize_kernel')

  !!inv_ovrlp%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%l, &
  !!                         iaction=SPARSE_FULL, id='inv_ovrlp%matrix_compr')
  inv_ovrlpp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
  tempp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
  inv_ovrlp_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
  kernel_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')




  ! Calculate S^1/2 * K * S^1/2. Take the value of S^1/2 from memory (was
  ! calculated in the last call to this routine or (it it is the first call)
  ! just before the call.
  !!call retransform_local(tmb%linmat%ovrlppowers_(1))
  !!call retransform_ext(iproc, nproc, tmb%linmat%l, &
  !!     tmb%linmat%kernel_%matrix_compr, tmb%linmat%ovrlppowers_(1)%matrix_compr)
  call retransform_ext(iproc, nproc, tmb%linmat%l, &
       tmb%linmat%ovrlppowers_(1)%matrix_compr, tmb%linmat%kernel_%matrix_compr)


  ! Calculate S^1/2 for the overlap matrix
  power=(/2,-2,1/)
  call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
       order_taylor, 3, power, -1, &
       imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
       ovrlp_mat=ovrlp, inv_ovrlp_mat=tmb%linmat%ovrlppowers_, &
       verbosity=0, &
       check_accur=.true., max_error=max_error, mean_error=mean_error)
  call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)


  !!tr=0.d0
  !!do iorb=1,tmb%orbs%norb
  !!    ind=tmb%linmat%l%matrixindex_in_compressed_fortransposed(iorb,iorb)
  !!    tr = tr + tmb%linmat%kernel_%matrix_compr(ind)
  !!end do
  !!write(*,*) 'trace',tr

  !!! Calculate S^-1/2 for the new overlap matrix
  !!call overlapPowerGeneral(iproc, nproc, order_taylor, 1, (/-2/), -1, &
  !!     imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
  !!     ovrlp_mat=ovrlp, inv_ovrlp_mat=tmb%linmat%ovrlppowers_, &
  !!     check_accur=.true., max_error=max_error, mean_error=mean_error)
  !!call check_taylor_order(mean_error, max_inversion_error, order_taylor)

  ! Calculate S^-1/2 * K * S^-1/2
  !!call retransform_local(tmb%linmat%ovrlppowers_(2))
  !!call retransform_ext(iproc, nproc, tmb%linmat%l, &
  !!     tmb%linmat%kernel_%matrix_compr, tmb%linmat%ovrlppowers_(2)%matrix_compr)
  call retransform_ext(iproc, nproc, tmb%linmat%l, &
       tmb%linmat%ovrlppowers_(2)%matrix_compr, tmb%linmat%kernel_%matrix_compr)


  call f_free_ptr(inv_ovrlpp)
  call f_free_ptr(tempp)
  call f_free(inv_ovrlp_compr_seq)
  call f_free(kernel_compr_seq)
  !!call f_free_ptr(inv_ovrlp%matrix_compr)

  call f_release_routine()

  !!contains

  !!    subroutine retransform_local(mat)
  !!        use sparsematrix, only: sequential_acces_matrix_fast2, sparsemm, &
  !!             uncompress_matrix_distributed2, compress_matrix_distributed, &
  !!             sequential_acces_matrix_fast
  !!        type(matrices),intent(in) :: mat
  !!        integer :: ncount

  !!        call f_routine(id='retransform_local')

  !!        call sequential_acces_matrix_fast2(tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, kernel_compr_seq)
  !!        call sequential_acces_matrix_fast2(tmb%linmat%l, &
  !!             mat%matrix_compr, inv_ovrlp_compr_seq)
  !!        call uncompress_matrix_distributed2(iproc, tmb%linmat%l, DENSE_MATMUL, &
  !!             mat%matrix_compr, inv_ovrlpp)

  !!        ncount=tmb%linmat%l%nfvctr*tmb%linmat%l%smmm%nfvctrp
  !!        if (ncount>0) then
  !!            call f_zero(ncount, tempp(1,1))
  !!        end if
  !!        call sparsemm(tmb%linmat%l, kernel_compr_seq, inv_ovrlpp, tempp)
  !!        if (ncount>0) then
  !!            call f_zero(ncount, inv_ovrlpp(1,1))
  !!        end if
  !!        call sparsemm(tmb%linmat%l, inv_ovrlp_compr_seq, tempp, inv_ovrlpp)

  !!        !call f_zero(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1))
  !!        call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
  !!             inv_ovrlpp, tmb%linmat%kernel_%matrix_compr)

  !!        call f_release_routine()

  !!    end subroutine retransform_local

end subroutine renormalize_kernel


subroutine allocate_precond_arrays(orbs, lzd, confdatarr, precond_convol_workarrays, precond_workarrays)
  use module_base, only: gp
  use module_types
  use locreg_operations
  implicit none
  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(confpot_data),dimension(orbs%norbp),intent(in) ::  confdatarr
  type(workarrays_quartic_convolutions),dimension(:),pointer,intent(inout) :: precond_convol_workarrays
  type(workarr_precond),dimension(:),pointer,intent(inout) :: precond_workarrays

  ! Local variables
  integer :: iorb, iiorb, ilr, ncplx
  real(kind=8) :: kx, ky, kz
  logical :: with_confpot

  allocate(precond_convol_workarrays(orbs%norbp))
  allocate(precond_workarrays(orbs%norbp))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      with_confpot = (confdatarr(iorb)%prefac/=0.d0)
      call init_local_work_arrays(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           with_confpot, precond_convol_workarrays(iorb))
      kx=orbs%kpts(1,orbs%iokpt(iorb))
      ky=orbs%kpts(2,orbs%iokpt(iorb))
      kz=orbs%kpts(3,orbs%iokpt(iorb))
      if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2 ) then
         ncplx=2
      else
         ncplx=1
      end if
      call allocate_work_arrays(lzd%llr(ilr)%geocode, lzd%llr(ilr)%hybrid_on, &
           ncplx, lzd%llr(ilr)%d, precond_workarrays(iorb))
  end do

end subroutine allocate_precond_arrays


subroutine deallocate_precond_arrays(orbs, lzd, precond_convol_workarrays, precond_workarrays)
  use module_base, only: gp
  use module_types
  use locreg_operations
  implicit none
  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(workarrays_quartic_convolutions),dimension(:),pointer,intent(inout) :: precond_convol_workarrays
  type(workarr_precond),dimension(:),pointer,intent(inout) :: precond_workarrays

  ! Local variables
  integer :: iorb, iiorb, ilr, ncplx
  real(kind=8) :: kx, ky, kz

  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      call deallocate_workarrays_quartic_convolutions(precond_convol_workarrays(iorb))
      kx=orbs%kpts(1,orbs%iokpt(iorb))
      ky=orbs%kpts(2,orbs%iokpt(iorb))
      kz=orbs%kpts(3,orbs%iokpt(iorb))
      if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2 ) then
         ncplx=2
      else
         ncplx=1
      end if
      call deallocate_work_arrays(lzd%llr(ilr)%geocode, lzd%llr(ilr)%hybrid_on, &
           ncplx, precond_workarrays(iorb))
  end do
  deallocate(precond_convol_workarrays)
  deallocate(precond_workarrays)

end subroutine deallocate_precond_arrays



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
  real(kind=8) :: dq, ebs, tt
  real(kind=8),dimension(input%nspin) :: e_homo, e_lumo
  integer :: ispin, norder_taylor, ind, iorb, norb, ntmb
  type(matrices),dimension(2) :: kernel
  logical,dimension(input%nspin) :: calculation_possible

  if (iproc==0) call yaml_comment('FOE calculation for HOMO-LUMO analysis',hfill='=')

  ! Check whether the gap calculation is possible.
  ! This is the case if there are more support functions than occupied orbitals.
  if (iproc==0) call yaml_mapping_open('Check possibility to calculate the gap')
  calculation_possible(1:input%nspin) = .true.
  do ispin=1,input%nspin
      if (ispin==1) then
          norb = orbs_KS%norbu
          ntmb = tmb%orbs%norbu
      else if (ispin==2) then
          norb = orbs_KS%norbd
          ntmb = tmb%orbs%norbd
      end if
      if (ntmb<=norb) then
          calculation_possible(ispin) = .false.
      end if
      if (iproc==0) then
          call yaml_mapping_open('Checking individual spin component')
          call yaml_map('ispin',ispin)
          call yaml_map('norb',norb)
          call yaml_map('ntmb',ntmb)
          call yaml_map('Calculation possible',calculation_possible(ispin))
          call yaml_mapping_close()
      end if
  end do
  if (iproc==0) then
      call yaml_map('Calculation possible',all(calculation_possible))
      call yaml_mapping_close()
  end if
                      
  if (all(calculation_possible)) then

      ! To determine the HOMO/LUMO, subtract/add one electrom for closed shell
      ! systems of one half electron for open shell systems.
      if (input%nspin==1) then
          dq = 1.d0
      else if (input%nspin==2) then
          dq = 0.5d0
      end if

      ! determine the HOMO
      if (iproc==0) call yaml_mapping_open('calculate HOMO kernel')
      kernel(1) = matrices_null()
      kernel(1)%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=SPARSE_TASKGROUP, id='kernel%matrix_compr')
      foe_obj = foe_data_null()
      call init_foe_wrapper(iproc, nproc, input, orbs_KS, 0.d0, foe_obj)
      do ispin=1,input%nspin
          call foe_data_set_real(foe_obj,"charge",foe_data_get_real(foe_obj,"charge",ispin)-dq,ispin)
      end do
      call foe_data_set_real(foe_obj,"fscale",1.d-2)
      norder_taylor = input%lin%order_taylor
      !call fermi_operator_expansion_new(iproc, nproc, &
      !     ebs, &
      !     .true., 2, &
      !     tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
      !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(1), foe_obj)
      call matrix_fermi_operator_expansion(iproc, nproc, bigdft_mpi%mpi_comm, &
           foe_obj, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
           tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%ovrlppowers_(2), kernel(1), &
           ebs, .true., 2)
      !call fermi_operator_expansion(iproc, nproc, &
      !     ebs, norder_taylor, input%lin%max_inversion_error, &
      !     .true., 2, &
      !     'HOMO', tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
      !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(1), foe_obj)
      do ispin=1,input%nspin
          e_homo(ispin) = foe_data_get_real(foe_obj,"ef",ispin)
      end do
      call foe_data_deallocate(foe_obj)
      if (iproc==0) call yaml_mapping_close()

      ! determine the LUMO
      if (iproc==0) call yaml_mapping_open('calculate LUMO kernel')
      kernel(2) = matrices_null()
      kernel(2)%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=SPARSE_TASKGROUP, id='kernel%matrix_compr')
      foe_obj = foe_data_null()
      call init_foe_wrapper(iproc, nproc, input, orbs_KS, 0.d0, foe_obj)
      do ispin=1,input%nspin
          call foe_data_set_real(foe_obj,"charge",foe_data_get_real(foe_obj,"charge",ispin)+dq,ispin)
      end do
      call foe_data_set_real(foe_obj,"fscale",1.d-2)
      norder_taylor = input%lin%order_taylor
      !call fermi_operator_expansion_new(iproc, nproc, &
      !     ebs, &
      !     .true., 2, &
      !     tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
      !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(2), foe_obj)
      call matrix_fermi_operator_expansion(iproc, nproc, bigdft_mpi%mpi_comm, &
           foe_obj, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
           tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%ovrlppowers_(2), kernel(1), &
           ebs, .true., 2)
      !call fermi_operator_expansion(iproc, nproc, &
      !     ebs, norder_taylor, input%lin%max_inversion_error, &
      !     .true., 2, &
      !     'LUMO', tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
      !     tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%ovrlppowers_(2), kernel(2), foe_obj)
      do ispin=1,input%nspin
          e_lumo(ispin) = foe_data_get_real(foe_obj,"ef",ispin)
      end do
      call foe_data_deallocate(foe_obj)
      if (iproc==0) call yaml_mapping_close()

      if (iproc==0) then
          call yaml_mapping_open('HOMO-LUMO analysis')
          call yaml_map('HOMO energy',e_homo)
          call yaml_map('LUMO energy',e_lumo)
          call yaml_map('HOMO-LUMO gap (Ha)',e_lumo-e_homo)
          call yaml_map('HOMO-LUMO gap (eV)',(e_lumo-e_homo)*Ha_eV)
          call yaml_mapping_close()
      end if

      !!if (iproc==0) then
      !!    tt = sqrt(kernel(2)%matrix_compr(1)-kernel(1)%matrix_compr(1))
      !!    do iorb=1,tmb%orbs%norb
      !!        write(*,*) 'iorb, val', iorb, (kernel(2)%matrix_compr(iorb)-kernel(1)%matrix_compr(iorb))/tt
      !!    end do
      !!end if

      call deallocate_matrices(kernel(1))
      call deallocate_matrices(kernel(2))

  end if

end subroutine calculate_gap_FOE


! WARNING: WILL MOST PROBABLY NOT WORK IN PARALLEL!!!!!!!!!!!
subroutine write_pexsi_matrices(iproc, nproc, smat_h, smat_s, matrix_compr_h, matrix_compr_s)
  use module_base
  use sparsematrix_init, only: sparsebigdft_to_ccs
  use sparsematrix_io, only:  write_ccs_matrix
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, &
                               assignment(=), SPARSE_FULL
  use sparsematrix, only: transform_sparsity_pattern
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(sparse_matrix),intent(in) :: smat_h, smat_s
  real(kind=8),dimension(smat_h%smmm%nvctrp),intent(inout) :: matrix_compr_h
  real(kind=8),dimension(smat_s%smmm%nvctrp),intent(inout) :: matrix_compr_s

  ! Local variables
  real(kind=8),dimension(:),pointer :: matrix_compr_sl
  integer,dimension(:),pointer :: row_ind, col_ptr

  call f_routine(id='write_pexsi_matrices')

  stop 'not correct'

  if (nproc/=1) then
      !call f_err_throw('not yet tested in parallel')
      call yaml_warning('not yet tested in parallel')
  end if

  matrix_compr_sl = sparsematrix_malloc_ptr(smat_h, iaction=SPARSE_FULL, id='matrix_compr_sl')
  call transform_sparsity_pattern(iproc, smat_h%nfvctr, smat_s%smmm%nvctrp_mm, smat_s%smmm%isvctr_mm, &
       smat_s%nseg, smat_s%keyv, smat_s%keyg, smat_s%smmm%line_and_column_mm, &
       smat_h%smmm%nvctrp, smat_h%smmm%isvctr, smat_h%smmm%nseg, smat_h%smmm%keyv, smat_h%smmm%keyg, &
       smat_h%smmm%istsegline, 'small_to_large', matrix_compr_s, matrix_compr_sl)

  row_ind = f_malloc_ptr(smat_h%nvctr,id='row_ind')
  col_ptr = f_malloc_ptr(smat_h%nfvctr,id='col_ptr')

  call sparsebigdft_to_ccs(smat_h%nfvctr, smat_h%nvctr, smat_h%nseg, smat_h%keyg, row_ind, col_ptr)

  call write_ccs_matrix('overlap_sparse_PEXSI.bin', smat_h%nfvctr, smat_h%nvctr, row_ind, col_ptr, matrix_compr_sl)
  call write_ccs_matrix('hamiltonian_sparse_PEXSI.bin', smat_h%nfvctr, smat_h%nvctr, row_ind, col_ptr, matrix_compr_h)

  call f_free_ptr(matrix_compr_sl)
  call f_free_ptr(row_ind)
  call f_free_ptr(col_ptr)

  call f_release_routine()

end subroutine write_pexsi_matrices
