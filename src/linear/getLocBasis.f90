!> @file
!! Linear version: Handle local basis set
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine get_coeff(iproc,nproc,scf_mode,orbs,at,rxyz,denspot,GPU,infoCoeff,&
    energs,nlpsp,SIC,tmb,fnrm,calculate_overlap_matrix,communicate_phi_for_lsumrho,&
    calculate_ham,extra_states,itout,it_scc,it_cdft,order_taylor,max_inversion_error,purification_quickreturn, &
    calculate_KS_residue,calculate_gap,&
    convcrit_dmin,nitdmin,curvefit_dmin,ldiis_coeff,reorder,cdft, updatekernel)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => get_coeff
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use constrained_dft
  use diis_sd_optimization
  use yaml_output
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, sparsematrix_malloc, &
                               DENSE_FULL, DENSE_PARALLEL, DENSE_MATMUL, assignment(=)
  use sparsematrix, only: uncompress_matrix, uncompress_matrix_distributed
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
  logical,intent(in):: calculate_overlap_matrix, communicate_phi_for_lsumrho, purification_quickreturn
  logical,intent(in) :: calculate_ham, calculate_KS_residue, calculate_gap
  type(DIIS_obj),intent(inout),optional :: ldiis_coeff ! for dmin only
  integer, intent(in), optional :: nitdmin ! for dmin only
  real(kind=gp), intent(in), optional :: convcrit_dmin ! for dmin only
  logical, intent(in), optional :: curvefit_dmin ! for dmin only
  type(cdft_data),intent(inout),optional :: cdft
  integer, intent(in) :: extra_states
  logical, optional, intent(in) :: reorder
  logical, optional, intent(in) :: updatekernel

  ! Local variables 
  integer :: iorb, info, ishift, ispin, ii, jorb, i, ishifts, ishiftm
  real(kind=8),dimension(:),allocatable :: hpsit_c, hpsit_f, eval
  real(kind=8),dimension(:,:),allocatable :: ovrlp_fullp, tempmat
  real(kind=8),dimension(:,:,:),allocatable :: matrixElements
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  logical :: update_kernel
  character(len=*),parameter :: subname='get_coeff'
  real(kind=gp) :: tmprtr, factor
  real(kind=8) :: max_deviation, mean_deviation, KSres, max_deviation_p,  mean_deviation_p

  call f_routine(id='get_coeff')

  if(calculate_ham) then
      call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
      call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,1), denspot%rhov, &
           tmb%ham_descr%comgp%nrecvbuf*tmb%ham_descr%comgp%nspin, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, &
           tmb%ham_descr%lzd)
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
  ! should eventually make this an input variable
  if (scf_mode==LINEAR_DIRECT_MINIMIZATION) then
      ! maybe need this for fragment calculations also, or make it an input?
     if (present(cdft)) then
        ! factor for scaling gradient
        factor=0.1d0
     else
        factor=1.0d0
     end if
  end if


  ! Calculate the overlap matrix if required.
  if(calculate_overlap_matrix) then
      if(.not.tmb%can_use_transposed) then
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
          tmb%can_use_transposed=.true.
      end if

      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
  end if

  ovrlp_fullp = sparsematrix_malloc(tmb%linmat%l,iaction=DENSE_PARALLEL,id='ovrlp_fullp')
  max_deviation=0.d0
  mean_deviation=0.d0
  do ispin=1,tmb%linmat%s%nspin
      ishift=(ispin-1)*tmb%linmat%s%nvctr
      call uncompress_matrix_distributed(iproc, tmb%linmat%s, DENSE_PARALLEL, &
           tmb%linmat%ovrlp_%matrix_compr(ishift+1:ishift+tmb%linmat%s%nvctr), ovrlp_fullp)
      call deviation_from_unity_parallel(iproc, nproc, tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, &
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

      if (tmb%ham_descr%npsidim_orbs > 0) call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))

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
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
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
               tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
          tmb%ham_descr%can_use_transposed=.true.
      end if

      hpsit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
      hpsit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)


      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
           tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%m, tmb%linmat%ham_)


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

  if (scf_mode/=LINEAR_FOE) then
      tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
      tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      !call uncompress_matrix(iproc, tmb%linmat%m, &
      !     inmat=tmb%linmat%ham_%matrix_compr, outmat=tmb%linmat%ham_%matrix)
      !call uncompress_matrix(iproc, tmb%linmat%s, &
      !     inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
      do ispin=1,tmb%linmat%m%nspin
          ishifts = (ispin-1)*tmb%linmat%s%nvctr
          ishiftm = (ispin-1)*tmb%linmat%m%nvctr
          call to_zero(tmb%linmat%m%nfvctr**2, tmb%linmat%ham_%matrix(1,1,ispin))
          tempmat = sparsematrix_malloc(tmb%linmat%m, iaction=DENSE_MATMUL, id='tempmat')
          call uncompress_matrix_distributed(iproc, tmb%linmat%m, DENSE_MATMUL, &
               tmb%linmat%ham_%matrix_compr(ishiftm+1:ishiftm+tmb%linmat%m%nvctr), tempmat)
          if (tmb%linmat%m%smmm%nfvctrp>0) then
              call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%smmm%nfvctrp, tempmat(1,1), 1, &
                   tmb%linmat%ham_%matrix(1,tmb%linmat%m%smmm%isfvctr+1,ispin), 1)
          end if
          call f_free(tempmat)
          if (nproc>1) then
              call mpiallred(tmb%linmat%ham_%matrix(1,1,ispin), tmb%linmat%m%nfvctr**2, &
                   mpi_sum, bigdft_mpi%mpi_comm)
          end if

          call to_zero(tmb%linmat%s%nfvctr**2, tmb%linmat%ovrlp_%matrix(1,1,ispin))
          tempmat = sparsematrix_malloc(tmb%linmat%s, iaction=DENSE_MATMUL, id='tempmat')
          call uncompress_matrix_distributed(iproc, tmb%linmat%s, DENSE_MATMUL, &
               tmb%linmat%ovrlp_%matrix_compr(ishifts+1:ishifts+tmb%linmat%s%nvctr), tempmat)
          if (tmb%linmat%m%smmm%nfvctrp>0) then
              call vcopy(tmb%linmat%s%nfvctr*tmb%linmat%s%smmm%nfvctrp, tempmat(1,1), 1, &
                   tmb%linmat%ovrlp_%matrix(1,tmb%linmat%s%smmm%isfvctr+1,ispin), 1)
          end if
          call f_free(tempmat)
          if (nproc>1) then
              call mpiallred(tmb%linmat%ovrlp_%matrix(1,1,ispin), tmb%linmat%s%nfvctr**2, &
                   mpi_sum, bigdft_mpi%mpi_comm)
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
              call dsygv_parallel(iproc, nproc, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%nproc_pdsyev, &
                   bigdft_mpi%mpi_comm, 1, 'v', 'l', tmb%linmat%m%nfvctr, &
                   matrixElements(1,1,1), tmb%linmat%m%nfvctr, matrixElements(1,1,2), tmb%linmat%m%nfvctr, &
                   eval, info)
          end if
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
                  call dscal(tmb%orbs%norb, -1.d0, matrixElements(1,iorb,1), 1)
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
     call daxpy(tmb%linmat%m%nvctr,-cdft%lag_mult,cdft%weight_matrix_%matrix_compr,1,tmb%linmat%ham_%matrix_compr,1)
     call timing(iproc,'constraineddft','OF') 
  end if

  if (scf_mode/=LINEAR_FOE) then
      ! Calculate the band structure energy and update kernel
      if (scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
              tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
              tmb%coeff,orbs,tmb%orbs,update_kernel)
      else if (present(cdft)) then
         ! for directmin we have the kernel already, but only the CDFT function not actual energy for CDFT
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
              tmb%linmat%kernel_, tmb%linmat%ham_, energs%ebs,&
              tmb%coeff,orbs,tmb%orbs,.false.)
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

  else ! foe


      ! TEMPORARY #################################################
      if (calculate_gap) then
          tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
          call uncompress_matrix(iproc, tmb%linmat%m, &
               inmat=tmb%linmat%ham_%matrix_compr, outmat=tmb%linmat%ham_%matrix)
          tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
          call uncompress_matrix(iproc, tmb%linmat%s, &
               inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
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

      if (iproc==0) call yaml_map('method','FOE')
      tmprtr=0.d0
      call foe(iproc, nproc, tmprtr, &
           energs%ebs, itout,it_scc, order_taylor, max_inversion_error, purification_quickreturn, &
           calculate_overlap_matrix, 1, FOE_ACCURATE, tmb, tmb%foe_obj)
      ! Eigenvalues not available, therefore take -.5d0
      tmb%orbs%eval=-.5d0

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
  if (iproc==0) call yaml_map('Coefficients available',scf_mode /= LINEAR_FOE)


  if (iproc==0) call yaml_mapping_close() !close kernel update

  call f_release_routine()

end subroutine get_coeff



subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
    fnrm,infoBasisFunctions,nlpsp,scf_mode,ldiis,SIC,tmb,energs_base,&
    nit_precond,target_function,&
    correction_orthoconstraint,nit_basis,&
    ratio_deltas,ortho_on,extra_states,itout,conv_crit,experimental_mode,early_stop,&
    gnrm_dynamic, min_gnrm_for_dynamic, can_use_ham, order_taylor, max_inversion_error, kappa_conv, method_updatekernel,&
    purification_quickreturn, correction_co_contra, cdft, input_frag, ref_frags)
  !
  ! Purpose:
  ! ========
  !   Calculates the localized basis functions phi. These basis functions are obtained by adding a
  !   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
  !   determined by minimizing the trace until the gradient norm is below the convergence criterion.
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => getLocalizedBasis, except_this_one_A => writeonewave
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_base, only: assignment(=), sparsematrix_malloc, sparsematrix_malloc_ptr, SPARSE_FULL
  use constrained_dft, only: cdft_data
  use module_fragments, only: system_fragment
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
  real(kind=8),intent(out) :: trH, fnrm
  real(kind=8),intent(inout) :: trH_old
  type(DFT_PSP_projectors),intent(inout) :: nlpsp
  integer,intent(in) :: scf_mode
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb
  type(SIC_data) :: SIC !<parameters for the SIC methods
  type(energy_terms),intent(in) :: energs_base
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
  real(kind=8),intent(out) :: ratio_deltas
  logical, intent(inout) :: ortho_on
  integer, intent(in) :: extra_states
  integer,intent(in) :: itout
  real(kind=8),intent(in) :: conv_crit, early_stop, gnrm_dynamic, min_gnrm_for_dynamic, kappa_conv
  logical,intent(in) :: experimental_mode, purification_quickreturn
  logical,intent(out) :: can_use_ham
  integer,intent(in) :: method_updatekernel
  logical,intent(in) :: correction_co_contra
  !these must all be present together
  type(cdft_data),intent(inout),optional :: cdft
  type(fragmentInputParameters),optional,intent(in) :: input_frag
  type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
 
  ! Local variables
  integer :: iorb, it, it_tot, ncount, ncharge, ii, kappa_satur, nit_exit, ispin
  !integer :: jorb, nspin
  !real(kind=8),dimension(:),allocatable :: occup_tmp
  real(kind=8) :: fnrmMax, meanAlpha, ediff_best, alpha_max, delta_energy, delta_energy_prev, ediff
  real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS, hpsit_c_tmp, hpsit_f_tmp, hpsi_tmp, psidiff
  real(kind=8),dimension(:),allocatable :: delta_energy_arr, hpsi_noprecond, kernel_compr_tmp, kernel_best
  logical :: energy_increased, overlap_calculated, energy_diff, energy_increased_previous, complete_reset, even
  logical :: calculate_inverse
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
  type(workarrays_quartic_convolutions),dimension(:),allocatable :: precond_convol_workarrays
  type(workarr_precond),dimension(:),allocatable :: precond_workarrays
  integer :: iiorb, ilr, i, ist

  call f_routine(id='getLocalizedBasis')

  energs = energy_terms_null()
  delta_energy_arr=f_malloc(nit_basis+6,id='delta_energy_arr')
  kernel_best=f_malloc(tmb%linmat%l%nvctr,id='kernel_best')
  energy_diff=.false.

  ovrlp_old%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%l, &
                           iaction=SPARSE_FULL, id='ovrlp_old%matrix_compr')

  ! Allocate all local arrays.
  call allocateLocalArrays()


  call timing(iproc,'getlocbasinit','ON')
  tmb%can_use_transposed=.false.

  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS
  ldiis%resetDIIS=.false.
  ldiis%immediateSwitchToSD=.false.
 
  call timing(iproc,'getlocbasinit','OF')

  overlap_calculated=.false.
  it=0
  it_tot=0
  !ortho=.true.
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
  call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,1), denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf*tmb%ham_descr%comgp%nspin, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, &
       tmb%ham_descr%lzd)

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

  ! Purify the initial kernel (only when necessary and if there is an even number of electrons)
  if (target_function/=TARGET_FUNCTION_IS_TRACE .and. even .and. scf_mode==LINEAR_FOE) then
      if (iproc==0) then
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_map('Initial kernel purification',.true.)
      end if
      overlap_calculated=.true.
      do ispin=1,tmb%linmat%l%nspin
          call purify_kernel(iproc, nproc, tmb, overlap_calculated, 1, 30, order_taylor, &
               max_inversion_error, purification_quickreturn, ispin)
      end do
      if (iproc==0) call yaml_mapping_close()
  end if

  if (itout==0) then
      nkappa_history=0
      kappa_history=0.d0
      has_already_converged=.false.
  end if

  iterLoop: do


      it=it+1
      it=max(it,1) !since it could become negative (2 is subtracted if the loop cycles)
      it_tot=it_tot+1

      fnrmMax=0.d0
      fnrm=0.d0
  
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


      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if (tmb%ham_descr%npsidim_orbs > 0)  call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%psi, tmb%ham_descr%psi)

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
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham2','OF')

      if (iproc==0) then
          call yaml_map('Hamiltonian Applied',.true.)
      end if

      ! Use this subroutine to write the energies, with some fake number
      ! to prevent it from writing too much
      if (iproc==0) then
          call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
      end if

      !if (iproc==0) write(*,'(a,5es16.6)') 'ekin, eh, epot, eproj, eex', &
      !              energs%ekin, energs%eh, energs%epot, energs%eproj, energs%exc

      if (iproc==0) then
          call yaml_map('Orthoconstraint',.true.)
      end if


      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               hpsi_tmp, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
          if (method_updatekernel==UPDATE_BY_FOE .or. method_updatekernel==UPDATE_BY_RENORMALIZATION) then
              if (method_updatekernel==UPDATE_BY_FOE) then
                  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                       tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
                  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
                       tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%m, tmb%linmat%ham_)
                  tmb%ham_descr%can_use_transposed=.true.
              else
                  tmb%ham_descr%can_use_transposed=.false.
              end if
              call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                   tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
              call vcopy(tmb%linmat%s%nspin*tmb%linmat%s%nvctr, tmb%linmat%ovrlp_%matrix_compr(1), 1, ovrlp_old%matrix_compr(1), 1)
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
                   tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
              !!write(*,*) 'calling FOE: sums(ovrlp)', &
              !!    sum(tmb%linmat%ovrlp_%matrix_compr(1:tmb%linmat%s%nvctr)), sum(tmb%linmat%ovrlp_%matrix_compr(tmb%linmat%s%nvctr+1:2*tmb%linmat%s%nvctr))
              if (iproc==0) call yaml_newline()
              if (iproc==0) call yaml_sequence_open('kernel update by FOE')
              if (method_updatekernel==UPDATE_BY_RENORMALIZATION) then
                  call renormalize_kernel(iproc, nproc, order_taylor, max_inversion_error, tmb, tmb%linmat%ovrlp_, ovrlp_old)
              else if (method_updatekernel==UPDATE_BY_FOE) then
                  call foe(iproc, nproc, 0.d0, &
                       energs%ebs, -1, -10, order_taylor, max_inversion_error, purification_quickreturn, &
                       .true., 0, &
                       FOE_FAST, tmb, tmb%foe_obj)
              end if
              if (iproc==0) call yaml_sequence_close()
          end if
      else
          !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
          !!    ist=0
          !!    do iorb=1,tmb%orbs%norbp
          !!        iiorb=tmb%orbs%isorb+iorb
          !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
          !!        ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
          !!        do i=1,ncount
          !!            ist=ist+1
          !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
          !!                write(4101,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%hpsi(ist)
          !!            else
          !!                write(4102,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%hpsi(ist)
          !!            end if
          !!        end do
          !!    end do
          !!end if
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
      end if

      ncount=tmb%ham_descr%collcom%ndimind_c
      if(ncount>0) call vcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*tmb%ham_descr%collcom%ndimind_f
      if(ncount>0) call vcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)

      ! optimize the tmbs for a few extra states
      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          kernel_compr_tmp = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='kernel_compr_tmp')
          call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !allocate(occup_tmp(tmb%orbs%norb), stat=istat)
          !call memocc(istat, occup_tmp, 'occup_tmp', subname)
          !call vcopy(tmb%orbs%norb, tmb%orbs%occup(1), 1, occup_tmp(1), 1)
          !call to_zero(tmb%orbs%norb,tmb%orbs%occup(1))
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
          !call transform_sparse_matrix(tmb%linmat%denskern, tmb%linmat%denskern_large, 'large_to_small')
      end if



      ! use hpsi_tmp as temporary array for hpsi_noprecond, even if it is allocated with a larger size
      !write(*,*) 'calling calc_energy_and.., correction_co_contra',correction_co_contra
      calculate_inverse = (target_function/=TARGET_FUNCTION_IS_HYBRID .or. method_updatekernel/=UPDATE_BY_RENORMALIZATION)
      call calculate_energy_and_gradient_linear(iproc, nproc, it, ldiis, fnrmOldArr, &
           fnrm_old, alpha, trH, trH_old, fnrm, fnrmMax, &
           meanAlpha, alpha_max, energy_increased, tmb, lhphiold, overlap_calculated, energs_base, &
           hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, hpsi_small, &
           experimental_mode, calculate_inverse, &
           correction_co_contra, hpsi_noprecond=hpsi_tmp, norder_taylor=order_taylor, &
           max_inversion_error=max_inversion_error, method_updatekernel=method_updatekernel, &
           precond_convol_workarrays=precond_convol_workarrays, precond_workarrays=precond_workarrays, &
           cdft=cdft, input_frag=input_frag, ref_frags=ref_frags)
      !fnrm_old=fnrm


      if (experimental_mode) then
          if (it_tot==1) then
              energy_first=trH
          end if
          if (iproc==0) call yaml_map('rel D',(trH-energy_first)/energy_first,fmt='(es9.2)')
          if ((trH-energy_first)/energy_first>early_stop .and. itout>0) then
              energy_diff=.true.
          end if
      end if

      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          call vcopy(tmb%linmat%l%nvctr, kernel_compr_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
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

      if (energy_increased .and. ldiis%isx==0) then
          !if (iproc==0) write(*,*) 'WARNING: ENERGY INCREASED'
          !if (iproc==0) call yaml_warning('The target function increased, D='&
          !              //trim(adjustl(yaml_toa(trH-ldiis%trmin,fmt='(es10.3)'))))
          if (iproc==0) then
              call yaml_newline()
              call yaml_map('iter',it,fmt='(i5)')
              call yaml_map('fnrm',fnrm,fmt='(es9.2)')
              call yaml_map('Omega',trH,fmt='(es22.15)')
              call yaml_map('D',ediff,fmt='(es9.2)')
              call yaml_map('D best',ediff_best,fmt='(es9.2)')
          end if
          tmb%ham_descr%can_use_transposed=.false.
          call vcopy(tmb%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
          can_use_ham=.false.
          call vcopy(tmb%linmat%l%nvctr, kernel_best(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
          trH_old=0.d0
          it=it-2 !go back one iteration (minus 2 since the counter was increased)
          overlap_calculated=.false.
          ! print info here anyway for debugging
          if (it_tot<2*nit_basis) then ! just in case the step size is the problem
              call yaml_mapping_close()
              call yaml_flush_document()
              !call bigdft_utils_flush(unit=6)
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
          call yaml_map('fnrm',fnrm,fmt='(es9.2)')
          call yaml_map('Omega',trH,fmt='(es22.15)')
          call yaml_map('D',ediff,fmt='(es9.2)')
          call yaml_map('D best',ediff_best,fmt='(es9.2)')
      end if

      ! Add some extra iterations if DIIS failed (max 6 failures are allowed before switching to SD)
      nit_exit=min(nit_basis+ldiis%icountDIISFailureTot,nit_basis+6)

      ! Determine whether the loop should be exited
      exit_loop(1) = (it>=nit_exit)
      exit_loop(2) = (it_tot>=3*nit_basis)
      exit_loop(3) = energy_diff
      exit_loop(4) = (fnrm<conv_crit .and. experimental_mode)
      exit_loop(5) = (experimental_mode .and. fnrm<dynamic_convcrit .and. fnrm<min_gnrm_for_dynamic &
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
          end if

          if (iproc==0) then
              !yaml output
              call yaml_mapping_close() !iteration
              call yaml_flush_document()
              !call bigdft_utils_flush(unit=6)
          end if

          exit iterLoop
      end if
      trH_old=trH

      if (ldiis%isx>0) then
          ldiis%mis=mod(ldiis%is,ldiis%isx)+1 !to store the energy at the correct location in the history
      end if
      call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
           lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho_on, psidiff, &
           experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)


      overlap_calculated=.false.
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmb%ham_descr%can_use_transposed) then
          tmb%ham_descr%can_use_transposed=.false.
      end if

      ! Estimate the energy change, that is to be expected in the next optimization
      ! step, given by the product of the force and the "displacement" .
      if (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode) then
          call estimate_energy_change(tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%l%nspin, psidiff, &
               hpsi_noprecond=hpsi_tmp, delta_energy=delta_energy)
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
          if(scf_mode/=LINEAR_FOE) then
              call reconstruct_kernel(iproc, nproc, order_taylor, tmb%orthpar%blocksize_pdsyev, &
                   tmb%orthpar%blocksize_pdgemm, orbs, tmb, overlap_calculated)
              if (iproc==0) call yaml_map('reconstruct kernel',.true.)
          else if (experimental_mode .and. .not.complete_reset) then
              if (method_updatekernel==UPDATE_BY_PURIFICATION) then
                  if (iproc==0) then
                      call yaml_map('purify kernel',.true.)
                      call yaml_newline()
                  end if
                  do ispin=1,tmb%linmat%l%nspin
                      call purify_kernel(iproc, nproc, tmb, overlap_calculated, 1, 30, &
                           order_taylor, max_inversion_error, purification_quickreturn, ispin)
                  end do
              else if (method_updatekernel==UPDATE_BY_FOE) then
                  if (iproc==0) then
                      call yaml_map('purify kernel',.false.)
                  end if
              end if
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
      call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
      call yaml_newline()
      call yaml_map('iter',it,fmt='(i5)')
      call yaml_map('fnrm',fnrm,fmt='(es9.2)')
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
      call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
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


      allocate(precond_convol_workarrays(tmb%orbs%norbp))
      allocate(precond_workarrays(tmb%orbs%norbp))
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          with_confpot = (tmb%confdatarr(iorb)%prefac/=0.d0)
          call init_local_work_arrays(tmb%lzd%llr(ilr)%d%n1, tmb%lzd%llr(ilr)%d%n2, tmb%lzd%llr(ilr)%d%n3, &
               tmb%lzd%llr(ilr)%d%nfl1, tmb%lzd%llr(ilr)%d%nfu1, &
               tmb%lzd%llr(ilr)%d%nfl2, tmb%lzd%llr(ilr)%d%nfu2, &
               tmb%lzd%llr(ilr)%d%nfl3, tmb%lzd%llr(ilr)%d%nfu3, &
               with_confpot, precond_convol_workarrays(iorb))
          kx=tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb))
          ky=tmb%orbs%kpts(2,tmb%orbs%iokpt(iorb))
          kz=tmb%orbs%kpts(3,tmb%orbs%iokpt(iorb))
          if (kx**2+ky**2+kz**2 > 0.0_gp .or. tmb%orbs%nspinor==2 ) then
             ncplx=2
          else
             ncplx=1
          end if
          call allocate_work_arrays(tmb%lzd%llr(ilr)%geocode, tmb%lzd%llr(ilr)%hybrid_on, &
               ncplx, tmb%lzd%llr(ilr)%d, precond_workarrays(iorb))
      end do


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
    do iorb=1,tmb%orbs%norbp
        iiorb=tmb%orbs%isorb+iorb
        ilr=tmb%orbs%inwhichlocreg(iiorb)
        call deallocate_workarrays_quartic_convolutions(precond_convol_workarrays(iorb))
        kx=tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb))
        ky=tmb%orbs%kpts(2,tmb%orbs%iokpt(iorb))
        kz=tmb%orbs%kpts(3,tmb%orbs%iokpt(iorb))
        if (kx**2+ky**2+kz**2 > 0.0_gp .or. tmb%orbs%nspinor==2 ) then
           ncplx=2
        else
           ncplx=1
        end if
        call deallocate_work_arrays(tmb%lzd%llr(ilr)%geocode, tmb%lzd%llr(ilr)%hybrid_on, &
             ncplx, precond_workarrays(iorb))
    end do
    deallocate(precond_convol_workarrays)
    deallocate(precond_workarrays)

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis



subroutine improveOrbitals(iproc, nproc, tmb, nspin, ldiis, alpha, gradient, experimental_mode)
  use module_base
  use module_types
  use module_interfaces, except_this_one => improveOrbitals
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspin
  type(DFT_wavefunction),intent(inout) :: tmb
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmb%orbs%norbp),intent(in) :: alpha
  real(kind=wp),dimension(max(tmb%npsidim_orbs,tmb%npsidim_comp)),intent(inout) :: gradient
  logical,intent(in) :: experimental_mode
  
  ! Local variables
  integer :: istart, iorb, iiorb, ilr, ncount

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
          call dscal(max(tmb%npsidim_orbs,tmb%npsidim_comp), ldiis%alphaDIIS, gradient, 1)
      end if
      call optimizeDIIS(iproc, nproc, max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, nspin, tmb%lzd, gradient, tmb%psi, ldiis, &
           experimental_mode)
  end if

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
  use module_interfaces
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, norb
  real(kind=8),dimension(norb, norb),intent(inout) :: HamSmall
  real(kind=8),dimension(norb, norb),intent(in) :: ovrlp
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

subroutine small_to_large_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
       orbs, phismall, philarge, to_global)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, npsidim_orbs_small, npsidim_orbs_large
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(npsidim_orbs_small),intent(in) :: phismall
  real(kind=8),dimension(npsidim_orbs_large),intent(out) :: philarge
  logical,intent(in),optional :: to_global
  
  ! Local variables
  integer :: ists, istl, iorb, ilr, sdim, ldim, nspin
  logical :: global

  call f_routine(id='small_to_large_locreg')

  if (present(to_global)) then
      global=to_global
  else
      global=.false.
  end if

  call timing(iproc,'small2large','ON') ! lr408t 
  ! No need to put arrays to zero, Lpsi_to_global2 will handle this.
  call to_zero(npsidim_orbs_large, philarge(1))
  ists=1
  istl=1
  do iorb=1,orbs%norbp
      ilr = orbs%inWhichLocreg(orbs%isorb+iorb)
      sdim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      if (global) then
          ldim=lzdsmall%glr%wfd%nvctr_c+7*lzdsmall%glr%wfd%nvctr_f
      else
          ldim=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
      end if
      nspin=1 !this must be modified later
      if (global) then
          call Lpsi_to_global2(iproc, sdim, ldim, orbs%norb, orbs%nspinor, nspin, lzdsmall%glr, &
               lzdsmall%llr(ilr), phismall(ists), philarge(istl))
      else
          call Lpsi_to_global2(iproc, sdim, ldim, orbs%norb, orbs%nspinor, nspin, lzdlarge%llr(ilr), &
               lzdsmall%llr(ilr), phismall(ists), philarge(istl))
      end if
      ists=ists+sdim
      istl=istl+ldim
  end do
  if(orbs%norbp>0 .and. ists/=npsidim_orbs_small+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',ists,'=ists /= npsidim_orbs_small+1=',npsidim_orbs_small+1
      stop
  end if
  if(orbs%norbp>0 .and. istl/=npsidim_orbs_large+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',istl,'=istl /= npsidim_orbs_large+1=',npsidim_orbs_large+1
      stop
  end if
       call timing(iproc,'small2large','OF') ! lr408t 
  call f_release_routine()
end subroutine small_to_large_locreg


subroutine large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
       orbs, philarge, phismall)
  use module_base
  use module_types
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
  !!call to_zero(npsidim_orbs_small, phismall(1))
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
  use module_interfaces, except_this_one => communicate_basis_for_density_collective
  use communications, only: transpose_switch_psir, transpose_communicate_psir, transpose_unswitch_psirt
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
      call initialize_work_arrays_sumrho(1,lzd%Llr(ilr),.true.,w)
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
  real(kind=8),dimension(tmbopt%linmat%l%nvctr),intent(out) :: kernel_best
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
      call vcopy(tmbopt%linmat%l%nvctr, tmbopt%linmat%kernel_%matrix_compr(1), 1, kernel_best(1), 1)
      !if(iproc==0) write(*,*) 'everything ok, copy last psi...'
      call vcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)

      ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
      ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
      if(ldiis%icountSDSatur>=10 .and. ldiis%isx==0 .or. ldiis%immediateSwitchToSD) then
          ldiis%icountSwitch=ldiis%icountSwitch+1
          idsx=max(ldiis%DIISHistMin,ldiis%DIISHistMax-ldiis%icountSwitch)
          if(idsx>0) then
              if(iproc==0) write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
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
              call vcopy(tmbopt%linmat%l%nvctr, kernel_best(1), 1, tmbopt%linmat%kernel_%matrix_compr(1), 1)
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
  use module_interfaces, except_this_one => reconstruct_kernel
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix
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
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if
     !call timing(iproc,'renormCoefComp','OF')

     call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     ! This can then be deleted if the transition to the new type has been completed.
     !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr

     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
  end if

  tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
  !!do i=1,tmb%linmat%s%nspin*tmb%linmat%s%nvctr
  !!    write(2500,'(a,i8,es16.5)') 'i, tmb%linmat%ovrlp_%matrix_compr(i)', i, tmb%linmat%ovrlp_%matrix_compr(i)
  !!end do
  call uncompress_matrix(iproc, tmb%linmat%s, &
       inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
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
  use module_interfaces, except_this_one => reorthonormalize_coeff
  use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, &
       allocate_matrices, deallocate_matrices
  use sparsematrix, only: orb_from_index
  use yaml_output, only: yaml_newline, yaml_map
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, norb
  integer, intent(in) :: blocksize_dsyev, blocksize_pdgemm, inversion_method
  type(orbitals_data), intent(in) :: basis_orbs   !number of basis functions
  type(sparse_matrix),intent(inout) :: basis_overlap
  type(sparse_matrix),dimension(basis_overlap%nspin),intent(inout) :: KS_overlap
  type(matrices),intent(inout) :: basis_overlap_mat
  real(kind=8),dimension(basis_overlap%nfvctr,norb),intent(inout) :: coeff
  type(orbitals_data), intent(in) :: orbs   !Kohn-Sham orbitals that will be orthonormalized and their parallel distribution
  ! Local variables
  integer :: ierr, ind, iorb, korb, llorb, jorb, ist
  integer :: npts_per_proc, ind_start, ind_end, indc, ispin, norbx, iseg, i
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, coefftrans
  real(kind=8), dimension(:,:), pointer :: ovrlp_coeff
  real(kind=8),dimension(:,:),pointer :: ovrlp_matrix, inv_ovrlp_matrix
  character(len=*),parameter:: subname='reorthonormalize_coeff'
  type(matrices) :: KS_ovrlp_, inv_ovrlp_
  integer,dimension(2) :: irowcol
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

      ovrlp_coeff=f_malloc_ptr((/norbx,norbx/), id='ovrlp_coeff')

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
             call to_zero(norbx**2,ovrlp_coeff(1,1))
          end if

          call f_free(coeff_tmp)
      else ! sparse - still less efficient than dense, also needs moving to a subroutine

         stop 'reorthonormalize_coeff: sparse version needs reworking'
         !also a problem with sparse at the moment - result not stored in correct arrays/allreduce etc

         !SM: need to fix the spin here
         call to_zero(norb**2, KS_ovrlp_%matrix(1,1,1))
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
          call mpiallred(ovrlp_coeff(1,1), norbx**2, mpi_sum, bigdft_mpi%mpi_comm)
          call timing(iproc,'renormCoefComm','OF')
          call timing(iproc,'renormCoefCom1','ON')
      end if

      !if (iproc==0) call yaml_map('ovrlp_coeff',ovrlp_coeff)
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
             call deviation_from_unity_parallel(iproc, nproc, norbx, orbs%norbp, orbs%isorb, &
                  ovrlp_coeff(1:orbs%norb,orbs%isorb+1:orbs%isorb+orbs%norbp), &
                  basis_overlap, max_error, mean_error)
             !!call deviation_from_unity_parallel(iproc, nproc, norbx, orbs%norbp, orbs%isorb, &
             !!     ovrlp_coeff, &
             !!     basis_overlap, max_error, mean_error)
          else
             ! It is necessary to call the routine since it has a built-in mpiallred.
             ! Use the first element of ovrlp_coeff; thanks to orbs%norbp==0 this should be safe
             call deviation_from_unity_parallel(iproc, nproc, orbs%norb, orbs%norbp, orbs%isorb, &
                  ovrlp_coeff(1:orbs%norb,1:orbs%norb), &
                  basis_overlap, max_error, mean_error)
          end if
      else
         call deviation_from_unity_parallel(iproc, 1, norbx, norbx, 0, ovrlp_coeff(1:norbx,1:norbx), &
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
         call f_free_ptr(ovrlp_coeff)
         call timing(iproc,'renormCoefCom2','ON')
         call gramschmidt_coeff_trans(iproc,nproc,orbs%norb,basis_orbs,basis_overlap,basis_overlap_mat,coeff)
         call timing(iproc,'renormCoefCom2','OF')
      else
         ! standard lowdin
         ! Not clean to use twice basis_overlap, but it should not matter as everything
         ! is done using the dense version

         inv_ovrlp_ = matrices_null()
         ! can not use the wrapper since it cannot distinguish between up and down spin
         !call allocate_matrices(KS_overlap, allocate_full=.true., matname='inv_ovrlp_', mat=inv_ovrlp_)
         inv_ovrlp_%matrix = f_malloc_ptr((/norbx,norbx,1/),id='inv_ovrlp_%matrix')
         

         if (norb==orbs%norb) then
             !SM: need to fix the spin here
             if (dense) call vcopy(norbx**2, ovrlp_coeff(1,1), 1, KS_ovrlp_%matrix(1,1,1), 1)
             !!do iorb=1,norbx
             !!    do jorb=1,norbx
             !!        write(2000+iproc,'(a,2i9,es13.5)') 'iorb, jorb, KS_ovrlp_%matrix(jorb,iorb,1)', iorb, jorb, KS_ovrlp_%matrix(jorb,iorb,1)
             !!    end do
             !!end do
             call overlapPowerGeneral(iproc, nproc, inversion_method, -2, &
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
             call overlap_minus_one_half_serial(iproc, 1, inversion_method, -2, blocksize_dsyev, &       
                  norbx, ovrlp_matrix, inv_ovrlp_matrix, check_accur=.false., smat=basis_overlap)
             call f_free_ptr(ovrlp_matrix)
         !    call overlapPowerGeneral(iproc, 1, inversion_method, -2, &
         !         blocksize_dsyev, norb, orbs, imode=2, ovrlp_smat=basis_overlap, inv_ovrlp_smat=basis_overlap, &
         !         ovrlp_mat=basis_overlap_mat, inv_ovrlp_mat=inv_ovrlp, &
         !         check_accur=.false., ovrlp=ovrlp_coeff, inv_ovrlp=ovrlp_coeff2)
         end if

         call timing(iproc,'renormCoefCom2','ON')

         call f_free_ptr(ovrlp_coeff)

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
                         inv_ovrlp_%matrix(1,orbs%isorb+1,1), orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
                    !@WARNING: THE FOLLOWING CALL IS NOT TESTED AND MIGHT BE WRONG!!
                    !!call dgemm('n', 't', basis_overlap%nfvctrp, norbx, norbx, 1.d0, coeff(basis_overlap%isfvctr,1), basis_overlap%nfvctr, &
                    !!     inv_ovrlp_%matrix(1,1,1), norbx, 0.d0, coeff_tmp(basis_overlap%isfvctr,1), basis_overlap%nfvctr)
                else !surely this isn't correct??
                    call dgemm('n', 't', basis_orbs%norb, orbs%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), basis_orbs%norb, &
                         inv_ovrlp_matrix(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
                end if
            else
               call to_zero(basis_overlap%nfvctr*norbx, coeff_tmp(1,1))
            end if

            if (nproc > 1) then
               call mpiallred(coeff_tmp(1,1), basis_overlap%nfvctr*norbx, mpi_sum, bigdft_mpi%mpi_comm)
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
                    call dgemm('n', 't', norbx, basis_overlap%nfvctrp, norbx, 1.d0, inv_ovrlp_%matrix(1,1,1), norbx, &
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

         call deallocate_matrices(inv_ovrlp_)
         if (norb/=orbs%norb) then
             call f_free_ptr(inv_ovrlp_matrix)
         end if

         call f_free(coeff_tmp)
      end if

      if (check_accuracy) then
         ovrlp_coeff=f_malloc_ptr((/norbx,norbx/), id='ovrlp_coeff')
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
            call to_zero(norbx**2,ovrlp_coeff(1,1))
         end if

         call f_free(coeff_tmp)

         if (nproc>1) then
            call mpiallred(ovrlp_coeff(1,1), norbx**2, mpi_sum, bigdft_mpi%mpi_comm)
         end if

         if (norb==orbs%norb) then
            ! Parallelization for nspin/=1 not done
            if (orbs%norbp>0) then
               call deviation_from_unity_parallel(iproc, nproc, orbs%norb, orbs%norbp, orbs%isorb, &
                    ovrlp_coeff(1:orbs%norb,orbs%isorb+1:orbs%isorb+orbs%norbp), &
                    basis_overlap, max_error, mean_error)
            else
               ! It is necessary to call the routine since it has a built-in mpiallred.
               ! Use the first element of ovrlp_coeff; thanks to orbs%norbp==0 this should be safe
               call deviation_from_unity_parallel(iproc, nproc, orbs%norb, orbs%norbp, &
                    orbs%isorb, ovrlp_coeff(1:orbs%norb,1:orbs%norb), &
                    basis_overlap, max_error, mean_error)
            end if
         else
            call deviation_from_unity_parallel(iproc, 1, norbx, norbx, 0, &
                 ovrlp_coeff(1:norbx,1:norbx), basis_overlap, max_error, mean_error)    
         end if

         if (iproc==0) print*,'Max deviation from unity following reorthonormalize_coeff',max_error
         if (iproc==0) print*,'Mean deviation from unity following reorthonormalize_coeff',mean_error

         !do iorb=1,norb
         !   do jorb=1,norb
         !      if (iproc==0) print*,jorb,iorb,ovrlp_coeff(jorb,iorb)
         !   end do
         !end do

         call f_free_ptr(ovrlp_coeff)
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
      call mpiallred(delta_energy, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  call f_release_routine()

end subroutine estimate_energy_change



subroutine purify_kernel(iproc, nproc, tmb, overlap_calculated, it_shift, it_opt, order_taylor, &
           max_inversion_error, purification_quickreturn, ispin)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => purify_kernel
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), matrices, &
                               matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix
  use foe_base, only: foe_data_get_real
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  type(DFT_wavefunction),intent(inout):: tmb
  logical,intent(inout):: overlap_calculated
  integer,intent(in) :: it_shift, it_opt
  logical,intent(in) :: purification_quickreturn
  integer,intent(in) :: ispin

  ! Local variables
  integer :: it, iorb, jorb, jsegstart, jsegend, jseg, jjorb, iiorb !info, lwork, 
  integer :: ishift, isshift, ilshift
  real(kind=8) :: trace_sparse, alpha, shift
  real(kind=8),dimension(:,:),allocatable :: ks, ksk, ksksk, kernel_prime
  !real(kind=8),dimension(:),allocatable :: eval, work
  character(len=*),parameter :: subname='purify_kernel'
  real(kind=8) :: diff, tr_KS, chargediff, max_error, mean_error
  !logical :: overlap_associated, inv_ovrlp_associated
  real(kind=8),dimension(2) :: bisec_bounds
  logical,dimension(2) :: bisec_bounds_ok
  !real(kind=8),dimension(:,:),pointer :: ovrlp_onehalf, ovrlp_minusonehalf
  type(matrices) :: ovrlp_onehalf_, ovrlp_minusonehalf_



  if (purification_quickreturn) then
      if (iproc==0) call yaml_warning('quick return in purification')
      if (iproc==0) call yaml_newline()
      return
  end if

  call f_routine(id='purify_kernel')

  isshift=(ispin-1)*tmb%linmat%s%nvctr
  ilshift=(ispin-1)*tmb%linmat%l%nvctr

  ovrlp_onehalf_ = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='ovrlp_onehalf_', mat=ovrlp_onehalf_)
  ovrlp_minusonehalf_ = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='ovrlp_minusonehalf_', mat=ovrlp_minusonehalf_)


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
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if
     !call timing(iproc,'renormCoefComp','OF')

     call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     ! This can then be deleted if the transition to the new type has been completed.
     !!tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr

     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
  end if


  tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l,iaction=DENSE_FULL,id='tmb%linmat%kernel_%matrix')
  call uncompress_matrix(iproc,tmb%linmat%s, &
       inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
  call uncompress_matrix(iproc, tmb%linmat%l, &
       inmat=tmb%linmat%kernel_%matrix_compr, outmat=tmb%linmat%kernel_%matrix)

  ks=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ks')
  ksk=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctrp/),id='ksk')
  ksksk=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ksksk')
  kernel_prime=f_malloc([tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr],id='kernel_prime')

  !ovrlp_onehalf=f_malloc_ptr((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ovrlp_onehalf')
  !ovrlp_minusonehalf=f_malloc_ptr((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ovrlp_minusonehalf')




  call timing(iproc,'purify_kernel ','ON') 

  if (tmb%linmat%l%nspin==1) then
      call dscal(tmb%linmat%l%nfvctr**2, 0.5d0, tmb%linmat%kernel_%matrix, 1)
  end if



  !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
  tr_KS=trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
        tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
        tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), ispin)
  if (iproc==0) then
      call yaml_map('tr(KS) before purification',tr_KS)
      call yaml_newline
  end if


  alpha=1.d-4
  chargediff=0.d0
  
  !!if (.not.associated(tmb%linmat%inv_ovrlp_large%matrix_compr)) then
  !!    inv_ovrlp_associated=.false.
  !!    !!allocate(tmb%linmat%inv_ovrlp_large%matrix_compr(tmb%linmat%inv_ovrlp_large%nvctr),stat=istat)
  !!    !!call memocc(istat,tmb%linmat%inv_ovrlp_large%matrix_compr,'tmb%linmat%inv_ovrlp_large%matrix_compr',subname)
  !!    tmb%linmat%inv_ovrlp_large%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp_large%nvctr,&
  !!        id='tmb%linmat%inv_ovrlp_large%matrix_compr')
  !!else
  !!    inv_ovrlp_associated=.true.
  !!end if

  if (it_shift>1) then
      call calculate_overlap_onehalf()
      call to_zero(tmb%linmat%l%nfvctr**2, kernel_prime(1,1))
      if (tmb%linmat%l%nfvctrp>0) then
          !SM: need to fix the spin here
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                     1.d0, tmb%linmat%kernel_%matrix, tmb%linmat%l%nfvctr, &
                     ovrlp_onehalf_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                     0.d0, ksksk, tmb%linmat%l%nfvctr) 
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                     1.d0, ovrlp_onehalf_%matrix, tmb%linmat%l%nfvctr, &
                     ksksk, tmb%linmat%l%nfvctr, &
                     0.d0, kernel_prime(1,tmb%linmat%l%isfvctr+1), tmb%linmat%l%nfvctr) 
      end if

      if (nproc > 1) then
          call mpiallred(kernel_prime(1,1), tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
      end if
  end if


  shift=0.d0
  bisec_bounds=0.d0
  bisec_bounds_ok=.false.

  shift_loop: do ishift=1,it_shift

  if (iproc==0) call yaml_newline()
  if (iproc==0) call yaml_map('shift of eigenvalues',shift,fmt='(es10.3)')

  if (iproc==0) call yaml_sequence_open('purification process')

      ! shift the eigenvalues of the density kernel, using ks as temporary variable
      if (shift/=0.d0) then
          if (ishift==1) stop 'eigenvalue shift not allowed for first iteration'
          do iorb=1,tmb%linmat%l%nfvctr
              do jorb=1,tmb%linmat%l%nfvctr
                  if (jorb==iorb) then
                      ks(jorb,iorb)=kernel_prime(jorb,iorb)+shift
                  else
                      ks(jorb,iorb)=kernel_prime(jorb,iorb)
                  end if
              end do
          end do
          !SM: need to fix the spin here
          call to_zero(tmb%linmat%l%nfvctr**2, tmb%linmat%kernel_%matrix(1,1,1))
          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ks, tmb%linmat%l%nfvctr, &
                         ovrlp_minusonehalf_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                         0.d0, ksksk, tmb%linmat%l%nfvctr) 
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ovrlp_minusonehalf_%matrix, tmb%linmat%l%nfvctr, &
                         ksksk, tmb%linmat%l%nfvctr, &
                         0.d0, tmb%linmat%kernel_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr) 
          end if
    

          if (nproc > 1) then
             !SM: need to fix the spin here
             call mpiallred(tmb%linmat%kernel_%matrix(1,1,1), tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
          end if
      end if


      do it=1,it_opt

          call to_zero(tmb%linmat%l%nfvctr**2, ks(1,1))
          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, tmb%linmat%kernel_%matrix(1,1,1), tmb%linmat%l%nfvctr, &
                         tmb%linmat%ovrlp_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                         0.d0, ks(1,tmb%linmat%l%isfvctr+1), tmb%linmat%l%nfvctr) 
          end if

          if (nproc > 1) then
              call mpiallred(ks(1,1), tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
          end if

          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ks(1,1), tmb%linmat%l%nfvctr, &
                         tmb%linmat%kernel_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                         0.d0, ksk(1,1), tmb%linmat%l%nfvctr)
          end if
          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ks(1,1), tmb%linmat%l%nfvctr, &
                         ksk(1,1), tmb%linmat%l%nfvctr, 0.d0, ksksk(1,1), tmb%linmat%l%nfvctr)
          end if


          diff=0.d0
          do iorb=tmb%linmat%l%isfvctr+1,tmb%linmat%l%isfvctr+tmb%linmat%l%nfvctrp
              iiorb=iorb-tmb%linmat%l%isfvctr
              jsegstart=tmb%linmat%l%istsegline(iorb)
              if (iorb<tmb%linmat%l%nfvctr) then
                  jsegend=tmb%linmat%l%istsegline(iorb+1)-1
              else
                  jsegend=tmb%linmat%l%nseg
              end if
              do jseg=jsegstart,jsegend
                  ! A segment is always on one line, therefore no double loop
                  do jorb=tmb%linmat%l%keyg(1,1,jseg),tmb%linmat%l%keyg(2,1,jseg)
                      jjorb=jorb
                      diff = diff + (ksk(jjorb,iiorb)-tmb%linmat%kernel_%matrix(jjorb,iorb,1))**2
                  end do
              end do
          end do

          call compress_matrix(iproc,tmb%linmat%l, &
               inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)
          !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
          tr_KS=trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
                tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
                tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), ispin)
          if (tmb%linmat%l%nspin==2) then
              chargediff=tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
          else if (tmb%linmat%l%nspin==1) then
              chargediff=2.d0*tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
          end if

          if (nproc > 1) then
              call mpiallred(diff, 1, mpi_sum, bigdft_mpi%mpi_comm)
          end if

          diff=sqrt(diff)
          if (iproc==0) then
              call yaml_newline()
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('iter',it)
              call yaml_map('diff from idempotency',diff,fmt='(es9.3)')
              call yaml_map('charge diff',chargediff,fmt='(es10.3)')
              !call yaml_map('alpha',alpha,fmt='(es8.2)')
              call yaml_mapping_close()
          end if

          call to_zero(tmb%linmat%l%nfvctr**2, tmb%linmat%kernel_%matrix(1,1,1))
          do iorb=1,tmb%linmat%l%nfvctrp
              iiorb=iorb+tmb%linmat%l%isfvctr
              do jorb=1,tmb%linmat%l%nfvctr
                  tmb%linmat%kernel_%matrix(jorb,iiorb,1) = 3.d0*ksk(jorb,iorb) - 2.d0*ksksk(jorb,iorb)
              end do
          end do

          if (nproc > 1) then
              call mpiallred(tmb%linmat%kernel_%matrix(1,1,1), tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
          end if

          if (diff<1.d-10) exit

      end do

      call compress_matrix(iproc,tmb%linmat%l, &
           inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)
      !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
      tr_KS=trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
            tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
            tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), ispin)
      if (tmb%linmat%l%nspin==2) then
          chargediff=tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
      else if (tmb%linmat%l%nspin==1) then
          chargediff=2.d0*tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
      end if

      if (iproc==0) call yaml_sequence_close

      if (abs(chargediff)<1.d-6) exit shift_loop

      if (chargediff>0) then
          ! make this the new upper bound for the bisection
          bisec_bounds(2)=shift
          ! choose new shift, based on whether the lower bound is known or not
          if (bisec_bounds_ok(1)) then
              shift=0.5d0*(bisec_bounds(1)+bisec_bounds(2))
          else
              shift=bisec_bounds(2)-0.01d0
          end if
          bisec_bounds_ok(2)=.true.
      end if
      if (chargediff<0) then
          ! make this the new lower bound for the bisection
          bisec_bounds(1)=shift
          ! choose new shift, based on whether the upper bound is known or not
          if (bisec_bounds_ok(2)) then
              shift=0.5d0*(bisec_bounds(1)+bisec_bounds(2))
          else
              shift=bisec_bounds(1)+0.01d0
          end if
          bisec_bounds_ok(1)=.true.
      end if



  end do shift_loop

  !if (iproc==0) call yaml_sequence_close

  if (tmb%linmat%l%nspin==1) then
      call dscal(tmb%linmat%l%nfvctr**2, 2.0d0, tmb%linmat%kernel_%matrix, 1)
  end if

  call timing(iproc,'purify_kernel ','OF') 

  call f_free(ks)
  call f_free(ksk)
  call f_free(ksksk)
  call f_free(kernel_prime)

  !call f_free_ptr(ovrlp_onehalf)
  !call f_free_ptr(ovrlp_minusonehalf)

  !if (.not.inv_ovrlp_associated) then
  !    call f_free_ptr(tmb%linmat%inv_ovrlp_large%matrix_compr)
  !end if


  call compress_matrix(iproc, tmb%linmat%l, inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)

  !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
  tr_KS=trace_sparse(iproc, nproc, tmb%orbs, tmb%linmat%s, tmb%linmat%l, &
        tmb%linmat%ovrlp_%matrix_compr(isshift+1:isshift+tmb%linmat%s%nvctr), &
        tmb%linmat%kernel_%matrix_compr(ilshift+1:ilshift+tmb%linmat%l%nvctr), ispin)
  if (iproc==0) then
      call yaml_newline()
      call yaml_map('tr(KS) after purification',tr_KS)
  end if


  call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  call f_free_ptr(tmb%linmat%kernel_%matrix)

  call deallocate_matrices(ovrlp_onehalf_)
  call deallocate_matrices(ovrlp_minusonehalf_)

  call f_release_routine()



      contains

        subroutine calculate_overlap_onehalf()
          ! Taylor approximation of S^1/2 and S^-1/2 up to higher order

          call overlapPowerGeneral(iproc, nproc, order_taylor, 2, -1, &
               imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
               ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_, check_accur=.true., &
               max_error=max_error, mean_error=mean_error)
          call check_taylor_order(mean_error, max_inversion_error, order_taylor)
          call overlapPowerGeneral(iproc, nproc, order_taylor, -2, -1, &
               imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
               ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=ovrlp_minusonehalf_, check_accur=.true., &
               max_error=max_error, mean_error=mean_error)
          call check_taylor_order(mean_error, max_inversion_error, order_taylor)
          !if (iproc==0) then
          !    call yaml_map('max error of S^-1/2',max_error,fmt='(es9.2)')
          !    call yaml_map('mean error of S^-1/2',mean_error,fmt='(es9.2)')
          !end if
      end subroutine calculate_overlap_onehalf

end subroutine purify_kernel



subroutine get_KS_residue(iproc, nproc, tmb, KSorbs, hpsit_c, hpsit_f, KSres)
  use module_base
  use module_types
  use module_interfaces, except_this_one => get_KS_residue
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices, &
                               sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix
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


  !gradmat%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='gradmat%matrix')
  tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
  call uncompress_matrix(iproc, tmb%linmat%m, inmat=gradmat%matrix_compr, outmat=gradmat%matrix)
  call uncompress_matrix(iproc, tmb%linmat%m, &
       inmat=tmb%linmat%ham_%matrix_compr, outmat=tmb%linmat%ham_%matrix)
  call uncompress_matrix(iproc, tmb%linmat%l, inmat=tmb%linmat%kernel_%matrix_compr, outmat=tmb%linmat%kernel_%matrix)
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
      call mpiallred(KH(1,1,1), tmb%linmat%l%nspin*tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
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
      call mpiallred(KHKH(1,1,1), tmb%linmat%l%nspin*tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
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
      call mpiallred(Kgrad(1,1,1), tmb%linmat%l%nspin*tmb%linmat%l%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
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
  use module_interfaces, only: overlapPowerGeneral
  use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, assignment(=), &
                               SPARSE_FULL, DENSE_FULL, DENSE_MATMUL, SPARSEMM_SEQ, &
                               matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed, &
                          uncompress_matrix_distributed

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
  !!real(8) :: tr
  !!integer :: ind, iorb

  call f_routine(id='renormalize_kernel')

  inv_ovrlp%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%l, &
                           iaction=SPARSE_FULL, id='inv_ovrlp%matrix_compr')
  inv_ovrlpp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
  tempp = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_MATMUL, id='inv_ovrlpp')
  inv_ovrlp_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')
  kernel_compr_seq = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSEMM_SEQ, id='inv_ovrlp_compr_seq')


  ! Calculate S^1/2 for the old overlap matrix
  call overlapPowerGeneral(iproc, nproc, order_taylor, 2, -1, &
       imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
       ovrlp_mat=ovrlp_old, inv_ovrlp_mat=inv_ovrlp, &
       check_accur=.true., max_error=max_error, mean_error=mean_error)
  call check_taylor_order(mean_error, max_inversion_error, order_taylor)

  ! Calculate S^1/2 * K * S^1/2
  call retransform_local()
  !!tr=0.d0
  !!do iorb=1,tmb%orbs%norb
  !!    ind=tmb%linmat%l%matrixindex_in_compressed_fortransposed(iorb,iorb)
  !!    tr = tr + tmb%linmat%kernel_%matrix_compr(ind)
  !!end do
  !!write(*,*) 'trace',tr

  ! Calculate S^-1/2 for the new overlap matrix
  if (iproc==0) write(*,*) 'in renormalize_kernel'
  call overlapPowerGeneral(iproc, nproc, order_taylor, -2, -1, &
       imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
       ovrlp_mat=ovrlp, inv_ovrlp_mat=tmb%linmat%ovrlp_minusonehalf_, &
       check_accur=.true., max_error=max_error, mean_error=mean_error)
  call check_taylor_order(mean_error, max_inversion_error, order_taylor)

  ! Calculate S^-1/2 * K * S^-1/2
  call retransform_local()


  call f_free_ptr(inv_ovrlpp)
  call f_free_ptr(tempp)
  call f_free(inv_ovrlp_compr_seq)
  call f_free(kernel_compr_seq)
  call f_free_ptr(inv_ovrlp%matrix_compr)

  call f_release_routine()

  contains

      subroutine retransform_local()
          use sparsematrix, only: sequential_acces_matrix_fast, sparsemm, &
               & uncompress_matrix_distributed, compress_matrix_distributed
          integer :: ncount

          call f_routine(id='retransform_local')

          call sequential_acces_matrix_fast(tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, kernel_compr_seq)
          call sequential_acces_matrix_fast(tmb%linmat%l, &
               inv_ovrlp%matrix_compr, inv_ovrlp_compr_seq)
          call uncompress_matrix_distributed(iproc, tmb%linmat%l, DENSE_MATMUL, &
               inv_ovrlp%matrix_compr, inv_ovrlpp)

          ncount=tmb%linmat%l%nfvctr*tmb%linmat%l%smmm%nfvctrp
          if (ncount>0) then
              call to_zero(ncount, tempp(1,1))
          end if
          call sparsemm(tmb%linmat%l, kernel_compr_seq, inv_ovrlpp, tempp)
          if (ncount>0) then
              call to_zero(ncount, inv_ovrlpp(1,1))
          end if
          call sparsemm(tmb%linmat%l, inv_ovrlp_compr_seq, tempp, inv_ovrlpp)

          !call to_zero(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1))
          call compress_matrix_distributed(iproc, nproc, tmb%linmat%l, DENSE_MATMUL, &
               inv_ovrlpp, tmb%linmat%kernel_%matrix_compr(tmb%linmat%l%isvctrp_tg+1:))

          call f_release_routine()

      end subroutine retransform_local

end subroutine renormalize_kernel
