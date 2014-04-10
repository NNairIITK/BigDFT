!> @file
!! Constrained DFT (based on linear version)
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module to perform constrained DFT calculations
module constrained_dft
  use module_base, only: gp,wp,f_err_throw
  use module_types
  use dynamic_memory
  implicit none

  private

  !> Need to avoid copying and pasting interfaces here - maybe shouldn't be a module?
  interface
       subroutine LocalHamiltonianApplication(iproc,nproc,at,npsidim_orbs,orbs,&
            Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
            energs,SIC,GPU,PotOrKin,pkernel,orbsocc,psirocc,dpbox,potential,comgp,hpsi_noconf,econf)
         use module_base
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: PotOrKin !< if true, only the potential operator is applied
         integer, intent(in) :: iproc,nproc,npsidim_orbs
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd 
         type(SIC_data), intent(in) :: SIC
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
         type(confpot_data), dimension(orbs%norbp) :: confdatarr
         real(wp), dimension(:), pointer :: pot
         !real(wp), dimension(*) :: pot
         type(energy_terms), intent(inout) :: energs
         real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout) :: hpsi
         type(GPU_pointers), intent(inout) :: GPU
         type(coulomb_operator), intent(in), optional :: pkernel
         type(orbitals_data), intent(in), optional :: orbsocc
         real(wp), dimension(:), pointer, optional :: psirocc
         type(denspot_distribution),intent(in),optional :: dpbox
         real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
         type(p2pComms),intent(inout), optional:: comgp
         real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout),optional :: hpsi_noconf
         real(gp),intent(out),optional :: econf
       end subroutine LocalHamiltonianApplication

       subroutine overlapPowerGeneral(iproc, nproc, iorder, power, blocksize, norb, orbs, imode, &
                  check_accur, ovrlp, inv_ovrlp, error, &
          ovrlp_smat, inv_ovrlp_smat, foe_obj)
         use module_base
         use module_types
         use sparsematrix_base, only: sparse_matrix
         use sparsematrix, only: compress_matrix, uncompress_matrix, transform_sparse_matrix
         implicit none
         integer,intent(in) :: iproc, nproc, iorder, blocksize, norb, power
         type(orbitals_data),intent(in) :: orbs
         integer,intent(in) :: imode
         logical,intent(in) :: check_accur
         real(kind=8),dimension(:,:),pointer,optional :: ovrlp
         real(kind=8),dimension(:,:),pointer,optional :: inv_ovrlp
         type(sparse_matrix), optional, intent(inout) :: ovrlp_smat, inv_ovrlp_smat
         type(foe_data),intent(in),optional :: foe_obj
         real(kind=8),intent(out),optional :: error
       end subroutine overlapPowerGeneral

  end interface


  type, public :: cdft_data
     real(wp), dimension(:), pointer :: weight_function ! the weight function defining the constraint
     type(sparse_matrix) :: weight_matrix ! matrix elements of the weight function between tmbs
     integer :: ndim_dens ! the dimension of the weight function
     real(gp) :: charge ! defines the value of the charge which is to be constrained
     real(gp) :: lag_mult ! the Lagrange multiplier used to enforce the constraint
     character(len=100) :: method
     integer, dimension(2) :: ifrag_charged ! make it allocatable eventually to allow for more charged fragments
     integer :: nfrag_charged
  end type cdft_data

  public :: nullify_cdft_data, cdft_data_allocate, cdft_data_free, cdft_data_init
  public :: calculate_weight_function, calculate_weight_matrix_using_density
  public :: calculate_weight_matrix_lowdin,  calculate_weight_matrix_lowdin_wrapper

contains

  !> CDFT: calculates the weight matrix w_ab via the expression S^1/2PS^1/2, where S is the overlap of the whole system
  !! CDFT: and P is a projector matrix onto the tmbs of the desired fragment
  !! CDFT: for standalone CDFT calculations, assuming one charged fragment, for transfer integrals assuming two fragments
  !! CDFT: where we constrain the difference - should later generalize this
  subroutine calculate_weight_matrix_lowdin_wrapper(cdft,tmb,input,ref_frags,calculate_overlap_matrix,meth_overlap)
    use module_fragments
    implicit none
    integer, intent(in) :: meth_overlap
    type(cdft_data), intent(inout) :: cdft
    type(input_variables),intent(in) :: input
    type(dft_wavefunction), intent(inout) :: tmb
    logical, intent(in) :: calculate_overlap_matrix
    type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

    integer :: nfrag_charged
    real(kind=gp), dimension(:,:), pointer :: ovrlp_half
    character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'

    ! wrapper here so we can modify charge and avoid restructuring the code as much whilst still using the routine elsewhere
    if (.not. input%lin%calc_transfer_integrals) then
       cdft%charge=ref_frags(input%frag%frag_index(cdft%ifrag_charged(1)))%nelec-input%frag%charge(cdft%ifrag_charged(1))
    end if

    if (input%lin%calc_transfer_integrals) then
       nfrag_charged=2
    else
       nfrag_charged=1
    end if

    ovrlp_half=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='ovrlp_half')
    call calculate_weight_matrix_lowdin(cdft%weight_matrix,nfrag_charged,cdft%ifrag_charged,tmb,input,&
         ref_frags,calculate_overlap_matrix,.true.,meth_overlap,ovrlp_half)
    call f_free_ptr(ovrlp_half)

  end subroutine calculate_weight_matrix_lowdin_wrapper


  subroutine calculate_weight_matrix_lowdin(weight_matrix,nfrag_charged,ifrag_charged,tmb,input,ref_frags,&
       calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap,ovrlp_half)
    use module_fragments
    use communications, only: transpose_localized
    use sparsematrix_base, only: sparse_matrix
    use sparsematrix, only: compress_matrix, uncompress_matrix
    implicit none
    type(sparse_matrix), intent(inout) :: weight_matrix
    type(input_variables),intent(in) :: input
    type(dft_wavefunction), intent(inout) :: tmb
    logical, intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
    type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags
    integer, intent(in) :: nfrag_charged, meth_overlap
    integer, dimension(2), intent(in) :: ifrag_charged
    real(kind=gp), dimension(:,:), pointer :: ovrlp_half
    !local variables
    integer :: ifrag,iorb,ifrag_ref,isforb,istat,ierr
    real(kind=gp), allocatable, dimension(:,:) :: proj_mat, proj_ovrlp_half, weight_matrixp
    character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
    real(kind=gp) :: error

    ! needs parallelizing/converting to sparse
    ! re-use overlap matrix if possible either before or after

    call f_routine(id='calculate_weight_matrix_lowdin')

    if (calculate_overlap_matrix) then
       if(.not.tmb%can_use_transposed) then
           if(.not.associated(tmb%psit_c)) then
               allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
               call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
           end if
           if(.not.associated(tmb%psit_f)) then
               allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
               call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
           end if
           call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
           tmb%can_use_transposed=.true.
       end if

       call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
            tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)
    end if   

    if (calculate_ovrlp_half) then
       tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='tmb%linmat%ovrlp%matrix')
       call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%ovrlp)
       call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 2, &
            tmb%orthpar%blocksize_pdsyev, tmb%orbs%norb, tmb%orbs, &
            imode=2, check_accur=.true., ovrlp=tmb%linmat%ovrlp%matrix, inv_ovrlp=ovrlp_half, error=error)
       call f_free_ptr(tmb%linmat%ovrlp%matrix)
    end if

    ! optimize this to just change the matrix multiplication?
    proj_mat=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='proj_mat')

    call to_zero(tmb%orbs%norb**2,proj_mat(1,1))
    isforb=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)
       if (ifrag==ifrag_charged(1)) then
          do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
             proj_mat(iorb+isforb,iorb+isforb)=1.0_gp
          end do
       end if
       if (nfrag_charged==2) then
          if (ifrag==ifrag_charged(2)) then
             do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                proj_mat(iorb+isforb,iorb+isforb)=-1.0_gp
             end do
          end if
       end if
       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
    end do

    proj_ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='proj_ovrlp_half')
    if (tmb%orbs%norbp>0) then
       call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
              tmb%orbs%norb, 1.d0, &
              proj_mat(1,1), tmb%orbs%norb, &
              ovrlp_half(1,tmb%orbs%isorb+1), tmb%orbs%norb, 0.d0, &
              proj_ovrlp_half(1,1), tmb%orbs%norb)
    end if
    call f_free(proj_mat)
    weight_matrixp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='weight_matrixp')
    if (tmb%orbs%norbp>0) then
       call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, & 
            tmb%orbs%norb, 1.d0, &
            ovrlp_half(1,1), tmb%orbs%norb, &
            proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
            weight_matrixp(1,1), tmb%orbs%norb)
    end if
    call f_free(proj_ovrlp_half)
    weight_matrix%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/), id='weight_matrix%matrix')
    if (bigdft_mpi%nproc>1) then
       call mpi_allgatherv(weight_matrixp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, weight_matrix%matrix, &
            tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
            mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    else
       call vcopy(tmb%orbs%norb*tmb%orbs%norb,weight_matrixp(1,1),1,weight_matrix%matrix(1,1),1)
    end if
    call f_free(weight_matrixp)
    call compress_matrix(bigdft_mpi%iproc,weight_matrix)
    call f_free_ptr(weight_matrix%matrix)
    call f_release_routine()

  end subroutine calculate_weight_matrix_lowdin


  !> CDFT: calculates the weight matrix w_ab given w(r)
  !! for the moment putting densities in global box and ignoring parallelization
  subroutine calculate_weight_matrix_using_density(iproc,cdft,tmb,at,input,GPU,denspot)
    use module_fragments
    use communications, only: transpose_localized, start_onesided_communication
    implicit none
    integer,intent(in) :: iproc
    type(cdft_data), intent(inout) :: cdft
    type(atoms_data), intent(in) :: at
    type(input_variables),intent(in) :: input
    type(dft_wavefunction), intent(inout) :: tmb
    type(DFT_local_fields), intent(inout) :: denspot
    type(GPU_pointers),intent(inout) :: GPU

    integer :: iall, istat
    !integer :: iorb, jorb
    real(kind=gp),dimension(:),allocatable :: hpsit_c, hpsit_f
    type(confpot_data),dimension(:),allocatable :: confdatarrtmp
    type(energy_terms) :: energs
    character(len=*),parameter :: subname='calculate_weight_matrix_using_density'

    call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))
    call start_onesided_communication(bigdft_mpi%iproc,bigdft_mpi%nproc,max(denspot%dpbox%ndimpot,1),cdft%weight_function, &
         tmb%ham_descr%comgp%nrecvbuf,tmb%ham_descr%comgp%recvbuf,tmb%ham_descr%comgp,tmb%ham_descr%lzd)

    allocate(confdatarrtmp(tmb%orbs%norbp))
    call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)

    call small_to_large_locreg(bigdft_mpi%iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
         tmb%orbs, tmb%psi, tmb%ham_descr%psi)

    if (tmb%ham_descr%npsidim_orbs > 0) call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))

    call full_local_potential(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%lzd,2, &
         denspot%dpbox,cdft%weight_function,denspot%pot_work,tmb%ham_descr%comgp)
    call LocalHamiltonianApplication(bigdft_mpi%iproc,bigdft_mpi%nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
         tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
         energs,input%SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
         potential=cdft%weight_function,comgp=tmb%ham_descr%comgp)

    deallocate(confdatarrtmp)

    print *,'CDFT Ekin,Epot,Eproj,Eh,Exc,Evxc',energs%ekin,energs%epot,energs%eproj,energs%eh,energs%exc,energs%evxc

    iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
    deallocate(denspot%pot_work, stat=istat)
    call memocc(istat, iall, 'denspot%pot_work', subname)

    ! calculate w_ab
    ! Calculate the matrix elements <phi|H|phi>.
    if(.not.tmb%ham_descr%can_use_transposed) then
        if(associated(tmb%ham_descr%psit_c)) then
            iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
            deallocate(tmb%ham_descr%psit_c, stat=istat)
            call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
        end if
        if(associated(tmb%ham_descr%psit_f)) then
            iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
            deallocate(tmb%ham_descr%psit_f, stat=istat)
            call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
        end if

        allocate(tmb%ham_descr%psit_c(tmb%ham_descr%collcom%ndimind_c), stat=istat)
        call memocc(istat, tmb%ham_descr%psit_c, 'tmb%ham_descr%psit_c', subname)
        allocate(tmb%ham_descr%psit_f(7*tmb%ham_descr%collcom%ndimind_f), stat=istat)
        call memocc(istat, tmb%ham_descr%psit_f, 'tmb%ham_descr%psit_f', subname)
        call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs, &
             tmb%ham_descr%collcom,tmb%ham_descr%psi,tmb%ham_descr%psit_c,tmb%ham_descr%psit_f,tmb%ham_descr%lzd)
        tmb%ham_descr%can_use_transposed=.true.
    end if

    allocate(hpsit_c(tmb%ham_descr%collcom%ndimind_c))
    call memocc(istat, hpsit_c, 'hpsit_c', subname)
    allocate(hpsit_f(7*tmb%ham_descr%collcom%ndimind_f))
    call memocc(istat, hpsit_f, 'hpsit_f', subname)
    call transpose_localized(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,  &
         tmb%ham_descr%collcom,tmb%hpsi,hpsit_c,hpsit_f,tmb%ham_descr%lzd)
    call calculate_overlap_transposed(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%collcom, &
         tmb%ham_descr%psit_c,hpsit_c,tmb%ham_descr%psit_f, hpsit_f,cdft%weight_matrix)
    iall=-product(shape(hpsit_c))*kind(hpsit_c)
    deallocate(hpsit_c, stat=istat)
    call memocc(istat, iall, 'hpsit_c', subname)
    iall=-product(shape(hpsit_f))*kind(hpsit_f)
    deallocate(hpsit_f, stat=istat)
    call memocc(istat, iall, 'hpsit_f', subname)

    ! debug
    !allocate(cdft%weight_matrix%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
    !call memocc(istat, cdft%weight_matrix%matrix, 'cdft%weight_matrix%matrix', subname)
    !call uncompress_matrix(bigdft_mpi%iproc,cdft%weight_matrix)
    !do iorb=1,tmb%orbs%norb
    !   do jorb=1,tmb%orbs%norb
    !      write(87,*) iorb,jorb,cdft%weight_matrix%matrix(iorb,jorb)
    !   end do
    !end do
    !deallocate(cdft%weight_matrix%matrix, stat=istat)
    !call memocc(istat, iall, 'cdft%weight_matrix%matrix', subname)
    ! end debug

  end subroutine calculate_weight_matrix_using_density


  !> CDFT: calculates the weight function w(r)
  !! for the moment putting densities in global box and ignoring parallelization
  subroutine calculate_weight_function(in,ref_frags,cdft,ndimrho_all_fragments,rho_all_fragments,tmb,atoms,rxyz,denspot)
    use module_fragments
    implicit none
    type(input_variables), intent(in) :: in
    type(system_fragment), dimension(in%frag%nfrag_ref), intent(inout) :: ref_frags
    type(cdft_data), intent(inout) :: cdft
    integer, intent(in) :: ndimrho_all_fragments
    real(kind=wp),dimension(ndimrho_all_fragments),intent(in) :: rho_all_fragments
    type(dft_wavefunction), intent(in) :: tmb
    type(atoms_data), intent(in) :: atoms ! just for plotting
    real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz ! just for plotting
    type(DFT_local_fields), intent(in) :: denspot ! just for plotting

    integer :: ifrag, ifrag_charged, ref_frag_charged, iorb_start, ipt

    ! unecessary recalculation here, but needs generalizing anyway
    iorb_start=1
    do ifrag=1,in%frag%nfrag
       if (in%frag%charge(ifrag)/=0) then
          ifrag_charged=ifrag
          exit
       end if
       iorb_start=iorb_start+ref_frags(in%frag%frag_index(ifrag))%fbasis%forbs%norb
    end do
    ref_frag_charged=in%frag%frag_index(ifrag_charged)

    ! check both densities are the same size (for now calculating fragment density in global box)
    if (ndimrho_all_fragments /= cdft%ndim_dens) stop 'Fragment density dimension does not match global density dimension'

    ! only need to calculate the fragment density for the fragment with a charge (stored in frag%fbasis%density)
    ! ideally would use already calculated kernel here, but charges could vary so would need an actual fragment
    ! rather than a reference fragment - could point to original fragment but modify nks...
    ! also converting input charge to integer which might not be the case
    call calculate_fragment_density(ref_frags(ref_frag_charged),ndimrho_all_fragments,tmb,iorb_start,&
         int(in%frag%charge(ifrag_charged)),atoms,rxyz,denspot)

    ! the weight function becomes the fragment density of ifrag_charged divided by the sum of all fragment densities
    ! using a cutoff of 1.e-6 for fragment density and 1.e-12 for denominator
    do ipt=1,ndimrho_all_fragments
        if (denspot%rhov(ipt) >= 1.0e-12) then
           cdft%weight_function(ipt)=ref_frags(ref_frag_charged)%fbasis%density(ipt)/denspot%rhov(ipt)
        else
           cdft%weight_function(ipt)=0.0_gp
        end if
    end do

    ! COME BACK TO THIS
    cdft%charge=ref_frags(ref_frag_charged)%nelec-in%frag%charge(ref_frag_charged)

    ! plot the weight function, fragment density and initial total density (the denominator) to check
    call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'fragment_density.cube', &
         atoms,rxyz,denspot%dpbox,1,ref_frags(ref_frag_charged)%fbasis%density)

    call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'weight_function.cube', &
         atoms,rxyz,denspot%dpbox,1,cdft%weight_function)

    call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'initial_density.cube', &
         atoms,rxyz,denspot%dpbox,1,denspot%rhov)

    ! deallocate fragment density here as we don't need it any more
    if (associated(ref_frags(ref_frag_charged)%fbasis%density)) call f_free_ptr(ref_frags(ref_frag_charged)%fbasis%density)

  end subroutine calculate_weight_function

  subroutine nullify_cdft_data(cdft)
    use sparsematrix_base, only: sparse_matrix_null
    implicit none
    type(cdft_data), intent(out) :: cdft
    cdft%charge=0
    cdft%lag_mult=0.0_gp
    cdft%ndim_dens=0
    nullify(cdft%weight_function)
    !call nullify_sparse_matrix(cdft%weight_matrix)
    cdft%weight_matrix=sparse_matrix_null()
  end subroutine nullify_cdft_data

  subroutine cdft_data_free(cdft)
    use sparsematrix_base, only: deallocate_sparse_matrix
    implicit none
    type(cdft_data), intent(inout) :: cdft

    character(len=200), parameter :: subname='cdft_data_free'

    !if (associated(cdft%weight_function)) call f_free_ptr(cdft%weight_function)
    call deallocate_sparse_matrix(cdft%weight_matrix, subname)
    call nullify_cdft_data(cdft)
  end subroutine cdft_data_free

  subroutine cdft_data_allocate(cdft,ham)
    use sparsematrix_base, only: sparse_matrix
    implicit none
    type(cdft_data), intent(inout) :: cdft
    type(sparse_matrix), intent(in) :: ham

    character(len=200), parameter :: subname='cdft_data_allocate'
    !!integer :: istat

    call f_routine(id='cdft_data_allocate')
    call sparse_copy_pattern(ham, cdft%weight_matrix, bigdft_mpi%iproc, subname)
    !cdft%weight_matrix%matrix_compr=f_malloc_ptr(cdft%weight_matrix%nvctr,id='cdft%weight_matrix%matrix_compr')
    !!allocate(cdft%weight_matrix%matrix_compr(cdft%weight_matrix%nvctr), stat=istat)
    !!call memocc(istat, cdft%weight_matrix%matrix_compr, 'cdft%weight_matrix%matrix_compr', subname)
    cdft%weight_matrix%matrix_compr=f_malloc_ptr(cdft%weight_matrix%nvctr,id='cdft%weight_matrix%matrix_compr')
    call f_release_routine()

  end subroutine cdft_data_allocate

  subroutine cdft_data_init(cdft,input_frag,ndimrho,transfer_int)
    implicit none
    type(cdft_data), intent(inout) :: cdft
    type(fragmentInputParameters), intent(in) :: input_frag
    integer, intent(in) :: ndimrho
    logical, intent(in) :: transfer_int

    integer :: ifrag, icharged

    ! For non-transfer integral calculation only one fragment should be charged
    ! For transfer integral calculation two should have charge
    ! the value is interpreted as the charge difference and so both should have the same charge
    ! we therefore do a calculation with a +ve difference followed by a -ve difference
    cdft%nfrag_charged=0
    do ifrag=1,input_frag%nfrag
       if (input_frag%charge(ifrag)/=0) cdft%nfrag_charged=cdft%nfrag_charged+1
    end do

    if (transfer_int) then
       if (cdft%nfrag_charged/=2) then
          call f_err_throw(&
               'Error in constrained DFT, two fragments must be charged for transfer integral calculation')
          return
       end if
    else ! could generalize this later (by summing charges and fragments), but for now keep as simplest scenario
       if (cdft%nfrag_charged/=1) then
          call f_err_throw(&
               'Error in constrained DFT, exactly one fragment must have a non-zero charge value'//&
               ' unless this is a transfer integral calculation')
          return
       end if
    end if

    icharged=1
    do ifrag=1,input_frag%nfrag
       if (input_frag%charge(ifrag)/=0) then
           cdft%ifrag_charged(icharged)=ifrag
           icharged=icharged+1
       end if
    end do

    if (cdft%nfrag_charged==2) then
       if (input_frag%charge(cdft%ifrag_charged(1))/=input_frag%charge(cdft%ifrag_charged(2))) then
          call f_err_throw('Error in constrained DFT, both fragments should have the same charge, '//& 
               'which is interpreted as the charge difference between the two')
          return
       end if
    end if

    cdft%charge=input_frag%charge(cdft%ifrag_charged(1))

    cdft%ndim_dens=ndimrho ! either size of fragment psi (add to fragment structure?), or size of entire simulation cell

    if (cdft%charge<0) then
       cdft%lag_mult=-0.05_gp ! pick some sensible initial value here
    else
       cdft%lag_mult=0.05_gp ! pick some sensible initial value here
    end if

    cdft%method='lowdin'
    !cdft%method='fragment_density'

  end subroutine cdft_data_init




end module constrained_dft
