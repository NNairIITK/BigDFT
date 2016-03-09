!> @file
!! Datatypes and associated methods relative to the localization regions
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Datatypes for localization regions descriptors
module orbitalbasis
  use module_defs, only: gp,wp
  use locregs
  use f_enums
  use module_types, only: orbitals_data,local_zone_descriptors
  use communications_base, only: comms_linear, comms_cubic
  use dictionaries, only: f_err_throw
  use locreg_operations, only: confpot_data
  use f_utils
  implicit none
  private

  type, public :: transposed_descriptor
     integer :: nspin !< spin multiplicity
     integer :: nup,ndown !<orbitals up and down
     integer :: nkpts !< number of k points
     !> lookup array indicating the dimensions
     !! for the overlap matrices in each quantum numbers
     integer, dimension(:,:), pointer :: ndim_ovrlp
     type(comms_cubic), pointer :: comms           !< communication objects for the cubic approach (can be checked to be
     type(comms_linear), pointer :: collcom        !< describes collective communication
     type(comms_linear), pointer :: collcom_sr     !< describes collective communication for the calculation of th
  end type transposed_descriptor

  !> iterator in the subspace of the wavefunctions
  type, public :: subspace
     integer :: ncplx !<  real(1) and complex (2) arrays
     integer :: nvctr !< number of components
     integer :: norb !< number of orbitals
     integer :: ispin !< spin index
     integer :: ikpt !< k-point id
     real(gp) :: kwgt !< k-point weight
     !>pointer on the occupation numbers of the subspace
     real(gp), dimension(:), pointer :: occup_ptr
     !>metadata
     integer :: ispsi
     integer :: ispsi_prev
     integer :: ikptp
     integer :: ise !<for occupation numbers
     integer :: ise_prev
     type(orbital_basis), pointer :: ob
  end type subspace

!  type, public :: support_function_descriptor
  type, public :: direct_descriptor
  !> wavelet localisation region (@todo gaussian part) for sfd)
     type(locreg_descriptors), pointer :: lr
!!$     integer :: ilr !inwhichlocreg or ilr
  end type direct_descriptor

  type, public :: ket ! support_function
     integer :: nphidim  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: ispsi !< Shift in the global array to store phi_wvl
     type(locreg_descriptors), pointer :: lr
     !id
     integer :: iorb
     !> spin
     integer:: nspin,nspinor
     real(gp), dimension(3) :: kpoint
     real(gp) :: kwgt,occup,spinval
     type(confpot_data) :: confdata
     real(wp), dimension(:), pointer :: phi_wvl !<coefficient
     !>starting address of potential
     integer :: ispot
     !> original orbital basis
     type(orbital_basis), pointer :: ob => null() !to be checked if it implies explicit save attribute

     !> number of localisation regions and current ilr
     integer :: nlrp,ilr,ilr_max,ilr_min,iorbp,ikpt,ikpt_max
  end type ket

  type, public :: orbital_basis
!!$     integer :: nbasis !< number of basis elements
!!$     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     !> descripto of each support function, of size nbasis
     type(direct_descriptor), dimension(:), pointer :: dd
     type(transposed_descriptor) :: td
     type(orbitals_data), pointer :: orbs !<metadata for the application of the hamiltonian
     type(confpot_data), dimension(:), pointer :: confdatarr !< data for the confinement potential
     real(wp), dimension(:), pointer :: phis_wvl !<coefficients in compact form for all the local sf
  end type orbital_basis

  public :: ob_ket_map,orbital_basis_iterator,ket_next_locreg,ket_next,local_hamiltonian_ket
  public :: orbital_basis_associate,orbital_basis_release,test_iterator,ket_next_kpt
  public :: ob_subket_ptr,precondition_ket

contains

  pure subroutine nullify_transposed_descriptor(td)
    implicit none
    type(transposed_descriptor), intent(out) :: td

    td%nspin     =f_none()
    td%nup       =f_none()
    td%ndown     =f_none()
    td%nkpts     =f_none()
    td%ndim_ovrlp=f_none()
    nullify(td%comms)
    nullify(td%collcom)
    nullify(td%collcom_sr)
  end subroutine nullify_transposed_descriptor

  pure subroutine nullify_orbital_basis(ob)
    implicit none
    type(orbital_basis), intent(out) :: ob
    nullify(ob%dd)
    nullify(ob%orbs)
    call nullify_transposed_descriptor(ob%td)
    nullify(ob%confdatarr)
    nullify(ob%phis_wvl)
  end subroutine nullify_orbital_basis

  !>this subroutine is not reinitializing each component of the
  !! iterator as some of them has to be set by the 'next' functions
  pure subroutine nullify_ket(k)
    use module_defs, only: UNINITIALIZED
    use locreg_operations, only: nullify_confpot_data
    implicit none
    type(ket), intent(inout) :: k

    !the orbital id
    k%iorb=-1
    k%nspin=-1
    k%nspinor=-1
    k%kpoint(1)=UNINITIALIZED(k%kpoint(1))
    k%kpoint(2)=UNINITIALIZED(k%kpoint(2))
    k%kpoint(3)=UNINITIALIZED(k%kpoint(3))
    k%kwgt=UNINITIALIZED(k%kwgt)
    k%occup=UNINITIALIZED(k%occup)
    k%spinval=UNINITIALIZED(k%spinval)
    call nullify_confpot_data(k%confdata)
    k%ispot=-1
    nullify(k%phi_wvl)
    nullify(k%ob)
  end subroutine nullify_ket

  function orbital_basis_iterator(ob) result(it)
    implicit none
    type(orbital_basis), intent(in), target :: ob
    type(ket) :: it

    call nullify_ket(it)
    it%ob => ob
    !number of parallel localization regions
    if (associated(ob%dd)) then
       it%nlrp = size(ob%dd)
    else
       it%nlrp=0
    end if
    !minimum value of locreg
    it%ilr_min=minval(ob%orbs%inwhichlocreg(ob%orbs%isorb+1:ob%orbs%isorb+ob%orbs%norbp))
    !zero value
    it%ilr=it%ilr_min-1
    !last value
    it%ilr_max=maxval(ob%orbs%inwhichlocreg(ob%orbs%isorb+1:ob%orbs%isorb+ob%orbs%norbp))
    !start orbital
    it%iorbp=0
    !start kpoint
    if (ob%orbs%norbp>0) then
       it%ikpt=ob%orbs%iokpt(1)-1
    else
       it%ikpt=-1
    end if
    !end kpoint
    it%ikpt_max=maxval(ob%orbs%iokpt)

  end function orbital_basis_iterator

  !>change the status of the iterator in the begininng of the next locreg (in the same kpoint)
  function ket_next_locreg(it,ikpt) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    !>desired kpoint. When absent, the function is legal only
    !! if the number of k-point is equal to one
    integer, intent(in), optional :: ikpt
    logical :: ok
    !local variablesmul
    integer :: ikpt_tmp

    ok=ket_is_valid(it)
!print *,'valid',ok,it%ilr,it%ilr_max,it%ikpt,it%ikpt_max,associated(it%ob)
    if (.not. ok) return

    if (present(ikpt)) then
       ikpt_tmp=ikpt
    else
       if (it%ob%orbs%nkpts > 1 .and. it%ilr_max > 1) &
            call f_err_throw('When more thank one kpt is present function'//&
            ' ket_next_locreg must have kpt specified',&
            err_name='BIGDFT_RUNTIME_ERROR')
       if (it%ob%orbs%norbp > 0) then
          ikpt_tmp=max(it%ikpt,minval(it%ob%orbs%iokpt))
       else
          ikpt_tmp=1
       end if
       it%ikpt=ikpt_tmp
    end if
!print *,'there',ok,it%ilr,ikpt_tmp
    find_next_lr: do while(it%ilr <= it%ilr_max)
       it%ilr=it%ilr+1
       if (dosome(it, it%ilr,ikpt_tmp)) exit find_next_lr
    end do find_next_lr
    ok= it%ilr <= it%ilr_max
!print *,'here',ok,it%ilr
    if (.not. ok ) then
       if (present(ikpt)) then
          it%ilr=it%ilr_min
       else
          call nullify_ket(it)
          return
       end if
    end if

    if (ok) call update_ket(it)

    it%iorbp=0

  end function ket_next_locreg

  !>change the status of the iterator in the begininng of the next kpoint
  function ket_next_kpt(it) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    logical :: ok

    ok=ket_is_valid(it)
!print *,ok,associated(it%ob),it%ilr,it%ilr_max,it%ikpt,it%ikpt_max
    if (.not. ok) return

    find_next_kpt: do while(it%ikpt <= it%ikpt_max)
       it%ikpt=it%ikpt+1
       it%ilr=it%ilr_min
       if (dosome_lr(it,it%ikpt)) exit find_next_kpt
    end do find_next_kpt
!print *,'here',it%ikpt,it%ikpt_max
    it%ilr=it%ilr-1 !for future usage with next locreg
    ok= it%ikpt <= it%ikpt_max
    if (.not. ok) then
       call nullify_ket(it)
       return
    end if

    if (ok) call update_ket(it)

    it%iorbp=0

  end function ket_next_kpt

  function ket_is_valid(it) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    logical :: ok
    !local variables
    logical :: noorb
    ok=.not. associated(it%ob) .or. it%ilr > it%ilr_max .or. &
         it%ikpt > it%ikpt_max
    ok=.not. ok
    noorb=.false.
    if (associated(it%ob)) noorb= it%ob%orbs%norbp==0
    !if we are at the end, nullify the iterator
    if ((it%ilr > it%ilr_max .and. it%ikpt > it%ikpt_max) .or. noorb)&
         call nullify_ket(it)

  end function ket_is_valid


  function ket_next(it,ilr,ikpt) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    integer, intent(in), optional :: ilr,ikpt
    logical :: ok
    !local variables
    integer :: ikpt_tmp

    ok=ket_is_valid(it)
!print *,'ket next',ok
    if (.not. ok) return

    !this is the most particular case, everything should match
    if (present(ilr) .and. present(ikpt)) then
       ok=dosome(it,ilr,ikpt)
    !the ilr is provided, therefore we know the kpt
    !(at least it is the sole one)
    else if (present(ilr)) then
!print *,'here',ilr,it%ilr,it%ikpt
       ok=dosome(it,ilr,it%ikpt)
!if the locregs are over but there are more than one kpt then increase
       if (.not. ok .and. it%ikpt < it%ikpt_max) then
          !this means that next_locreg has not been called
          !therefore re-increase the locreg index
          it%ikpt=it%ikpt+1
          !this should be one as otherwise both the iterators should be called
          if (it%ilr /= it%ilr_min) &
               call f_err_throw('Internal error in ket_next',err_name='BIGDFT_RUNTIME_ERROR')
          ok=dosome(it,ilr,it%ikpt)
       end if
!print *,'there',ilr,it%ilr,it%ikpt,it%iorbp
    !if only the kpt is provided, we might have to find the first valid
    !locreg
    else if (present(ikpt)) then
!print *,'again',it%iorbp,it%ikpt,it%ilr
       ok = dosome_lr(it,ikpt)
       if (.not. ok .and. it%ilr < it%ilr_max) then
          !this means that next_locreg has not been called
          !therefore re-increase the locreg index
          it%ilr=it%ilr+1
          ok = dosome_lr(it,ikpt)
       end if
       if (.not. ok) it%ilr=it%ilr_min
!print *,'again',it%iorbp,it%ikpt,it%ilr,ok
    else
    !increment the localization region
    !first increment the orbital, otherwise the locreg, otherwise k-point
       ikpt_tmp=max(it%ikpt,minval(it%ob%orbs%iokpt)) !to avoid the first useless passage
!print *,'here',it%ilr,ikpt_tmp,it%iorbp
       find_next_ikpt: do while( .not. dosome_lr(it,ikpt_tmp) .and.&
            ikpt_tmp <= it%ikpt_max)
          ikpt_tmp=ikpt_tmp+1
          it%ilr=it%ilr_min
       end do find_next_ikpt
!print *,'there',it%ilr,ikpt_tmp,it%iorbp
       it%ikpt=ikpt_tmp
       ok= it%ikpt <= it%ikpt_max
       !the nullification of the ket here makes sense only for the unspecified case
       if (.not. ok) then
          call nullify_ket(it)
          return
       end if
    end if

    !at this point the iorbp and ilr are determined, the iterator can be updated
    if (ok) call update_ket(it)

  end function ket_next


  !> Find next valid iorb in this lr
  !! put iorbp to zero if ilr_tmp is not found
  function dosome(it,ilr_tmp,ikpt_tmp)
    implicit none
    type(ket), intent(inout) :: it
    integer, intent(in) :: ilr_tmp,ikpt_tmp
    logical :: dosome
    dosome=.false.
    find_iorb: do while(.not. dosome .and. it%iorbp < it%ob%orbs%norbp)
       it%iorbp=it%iorbp+1
       !check if this localisation region is used by one of the orbitals
       dosome = (it%ob%orbs%inwhichlocreg(it%iorbp+it%ob%orbs%isorb) == ilr_tmp) .and. &
            (it%ob%orbs%iokpt(it%iorbp) == ikpt_tmp)
!print *,'dosome',dosome,it%iorbp,ikpt_tmp,ilr_tmp,it%ob%orbs%iokpt(it%iorbp)
    end do find_iorb
    if (.not. dosome) it%iorbp=0
  end function dosome


  function dosome_lr(it,ikpt) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    integer, intent(in) :: ikpt
    logical :: ok
    !local variables
    integer :: ilr_tmp

    ilr_tmp=it%ilr
    find_next_ilr_k: do while(.not. dosome(it,ilr_tmp,ikpt) .and. &
         ilr_tmp <= it%ilr_max)
       ilr_tmp=ilr_tmp+1
    end do find_next_ilr_k
    it%ilr=ilr_tmp
    ok= it%ilr <= it%ilr_max
    if (.not. ok) it%iorbp=0
  end function dosome_lr


  subroutine update_ket(k)
    implicit none
    type(ket), intent(inout) :: k
    !local variables
    integer :: ikpt,iorbq,nvctr
    !the orbital id
    k%iorb=k%ob%orbs%isorb+k%iorbp
    k%nspin=k%ob%orbs%nspin
    k%nspinor=k%ob%orbs%nspinor
    !k-point, spin and confinement
    ikpt=k%ob%orbs%iokpt(k%iorbp)
    if (ikpt /= k%ikpt) call f_err_throw('Internal error in update ket',err_name='BIGDFT_RUNTIME_ERROR')
    k%kpoint=k%ob%orbs%kpts(:,ikpt)
    k%kwgt=k%ob%orbs%kwgts(ikpt)
    k%occup=k%ob%orbs%occup(k%iorb)
    k%spinval=k%ob%orbs%spinsgn(k%iorb)
    if (associated(k%ob%confdatarr)) k%confdata=k%ob%confdatarr(k%iorbp)
    !shifts metadata
    k%ispot=k%ob%orbs%ispot(k%iorbp)
    !find the psi shift for the association
    k%ispsi=1
    if (associated(k%ob%dd)) then
       do iorbq=1,k%iorbp-1
          nvctr=k%ob%dd(iorbq)%lr%wfd%nvctr_c+7*k%ob%dd(iorbq)%lr%wfd%nvctr_f
          k%ispsi=k%ispsi+nvctr*k%nspinor
       end do
       k%lr=>k%ob%dd(k%iorbp)%lr
       k%nphidim=(k%lr%wfd%nvctr_c+7*k%lr%wfd%nvctr_f)*k%nspinor
    end if
    if (associated(k%ob%phis_wvl)) k%phi_wvl=>ob_ket_map(k%ob%phis_wvl,k)
  end subroutine update_ket

  function ob_ket_map(ob_ptr,it)
    use f_precisions, only: f_address,f_loc
    implicit none
    real(wp), dimension(:), target :: ob_ptr !<coming from orbital_basis
    type(ket), intent(in) :: it
    real(wp), dimension(:), pointer :: ob_ket_map

    !might add the calculation of ispsi here

    !here we can add the check of the pertinence of value of ispsi
    ob_ket_map => ob_ptr(it%ispsi:it%ispsi+it%nphidim-1)

    !also, the f_loc function might be used to determine if the association
    !corresponds to a shallow or a deep copy

  end function ob_ket_map
  !the iterator must go in order of localization regions

  function ob_subket_ptr(it,ispinor) result(k)
    use f_precisions, only: f_address,f_loc
    use dynamic_memory
    implicit none
    type(ket), intent(in) :: it
    integer, intent(in) :: ispinor
    real(wp), dimension(:), pointer :: k
    !local variables
    integer :: istart,nvctr
    
    nvctr=it%lr%wfd%nvctr_c+7*it%lr%wfd%nvctr_f
    istart=(ispinor-1)*nvctr

    k => f_subptr(it%phi_wvl,istart+1 .to. istart+nvctr)

  end function ob_subket_ptr



  !> This function gives the number of components
  !! if the ket is in non-colinear spin description, this value is four
  !! otherwise it counts the number of complex components of the key
  !! result is one for real functions
  pure function nspinor(spin_enum)
    implicit none
    type(f_enumerator), intent(in) :: spin_enum
    integer :: nspinor
!!$    if (spin_enum.hasattr. 'NONCOLL') then
!!$       nspinor=4
!!$    else
!!$    end if
    nspinor=1
  end function nspinor

  function subspace_iterator(ob) result(ss)
    implicit none
    type(orbital_basis), intent(in), target :: ob
    type(subspace) :: ss

    call nullify_subspace(ss)
    ss%ob => ob
    ss%ispin=0
    ss%ikptp=1
    ss%ispsi_prev=1
    ss%ise_prev=0
  end function subspace_iterator

  !case of subspace iterators
  pure subroutine nullify_subspace(ss)
    implicit none
    type(subspace), intent(out) :: ss
    ss%ncplx      =f_none()
    ss%nvctr      =f_none()
    ss%norb       =f_none()
    ss%ispin      =f_none()
    ss%ikpt       =f_none()
    ss%kwgt       =f_none()
    ss%occup_ptr  =f_none()
    ss%ispsi  =f_none()
    ss%ispsi_prev =f_none()
    ss%ise        =f_none()
    ss%ise_prev   =f_none()
    nullify(ss%ob)
  end subroutine nullify_subspace

  pure function subspace_is_valid(it) result(ok)
    implicit none
    type(subspace), intent(in) :: it
    logical :: ok

    ok=associated(it%ob)
  end function subspace_is_valid

  function subspace_next(it) result(ok)
    implicit none
    type(subspace), intent(inout) :: it
    logical  :: ok

    ok=subspace_is_valid(it)
    if (ok) call increment_subspace(it)
    ok=subspace_is_valid(it)
  end function subspace_next

  !pure
  subroutine increment_subspace(it)
    use module_base, only: bigdft_mpi
    implicit none
    type(subspace), intent(inout) :: it
    !local variables
    integer :: ncomp,nvctrp,nspinor,ist,norbs

    do
       if (it%ispin < it%ob%td%nspin) then
          it%ispin=it%ispin+1
       else if (it%ikptp < it%ob%orbs%nkptsp) then
          it%ikptp=it%ikptp+1
          it%ispin=1
       else
          call nullify_subspace(it)
          exit
       end if

       it%ispsi=it%ispsi_prev
       it%ikpt=it%ob%orbs%iskpts+it%ikptp

       call orbitals_and_components(bigdft_mpi%iproc,it%ikpt,it%ispin,&
            it%ob%orbs,it%ob%td%comms,&
            nvctrp,it%norb,norbs,ncomp,nspinor)
       if (nvctrp == 0) cycle

       if (it%ispin==1) then
          it%ise=0
       else
          it%ise=it%ise_prev
       end if

       it%ise_prev=it%norb

       it%ncplx=1
       it%nvctr=ncomp*nvctrp
       if (nspinor/=1) it%ncplx=2

       it%kwgt=it%ob%orbs%kwgts(it%ikpt)
       ist=(it%ikpt-1)*it%ob%orbs%norb+1+it%ise
       it%occup_ptr=>it%ob%orbs%occup(ist:ist+it%norb-1)

       it%ispsi_prev=it%ispsi_prev+nvctrp*it%norb*nspinor
       exit
    end do

  end subroutine increment_subspace

  subroutine local_hamiltonian_ket(psi,hgrids,ipotmethod,xc,pkernel,wrk_lh,psir,vsicpsir,hpsi,pot,eSIC_DCi,alphaSIC,epot,ekin)
    use module_xc, only: xc_info, xc_exctXfac
    use locreg_operations, only: workarr_locham,psir_to_vpsi, isf_to_daub_kinetic
    use Poisson_Solver, only: coulomb_operator
    use dynamic_memory, only : f_memcpy
    use wrapper_linalg, only: axpy
    implicit none
    type(ket), intent(in) :: psi
    integer, intent(in) :: ipotmethod
    type(xc_info), intent(in) :: xc
    type(workarr_locham), intent(inout) :: wrk_lh
    real(gp), dimension(3), intent(in) :: hgrids
    real(gp), intent(in) :: alphaSIC
    real(gp), intent(out) :: eSIC_DCi, epot, ekin
    real(wp), dimension(psi%lr%d%n1i*psi%lr%d%n2i*psi%lr%d%n3i,psi%nspinor), intent(inout) :: psir !to be unified with vsicpsir
    real(wp), dimension(*), intent(in) :: pot
    !> the PSolver kernel which should be associated for the SIC schemes
    type(coulomb_operator), intent(in) :: pkernel
    real(wp), dimension(:,:), allocatable, intent(inout) :: vsicpsir
    real(wp), dimension(psi%nphidim), intent(inout) :: hpsi
    !local variables
    integer :: npoints, ispot, npot
    real(gp) :: eSICi, exctXcoeff,fi,hfac
    real(gp), dimension(3) :: hh

    if (psi%iorbp==0) call f_err_throw('Illegal iorbp for local hamiltonian ket',err_name='BIGDFT_RUNTIME_ERROR')

    epot=0.d0
    ekin=0.d0

    !aliases
    hh=0.5_gp*hgrids
    npoints=psi%lr%d%n1i*psi%lr%d%n2i*psi%lr%d%n3i

    !components of the potential (four or one, depending on the spin)
    npot=psi%nspinor
    if (psi%nspinor == 2) npot=1

    exctXcoeff=xc_exctXfac(xc)
    call daub_to_isf_locham(psi%nspinor,psi%lr,wrk_lh,psi%phi_wvl,psir)

    !calculate the ODP, to be added to VPsi array

    !Perdew-Zunger SIC scheme
    eSIC_DCi=0.0_gp
    if (ipotmethod == 2) then
       !in this scheme the application of the potential is already done
!!$       call PZ_SIC_potential(psi%iorb,psi%lr,psi%ob%orbs,xc,&
!!$            hh,pkernel,psir,vsicpsir,eSICi,eSIC_DCi)
       fi=psi%kwgt*psi%occup
       hfac=fi/product(hh)

       call PZ_SIC_potential(psi%nspin,psi%nspinor,hfac,psi%spinval,psi%lr,xc,&
            hh,pkernel,psir,vsicpsir,eSICi,eSIC_DCi)

       !NonKoopmans' correction scheme
    else if (ipotmethod == 3) then
       !in this scheme first we have calculated the potential then we apply it
!!$       call vcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,&
!!$            psir(1,1),1,vsicpsir(1,1),1)
       call f_memcpy(src=psir,dest=vsicpsir)
       !for the moment the ODP is supposed to be valid only with one lr
       call psir_to_vpsi(npot,psi%nspinor,psi%lr,&
            pot(npoints*psi%nspin+&
            (psi%iorbp-1)*npoints*psi%nspinor+1),&
            vsicpsir,eSICi)
    end if

    call psir_to_vpsi(npot,psi%nspinor,psi%lr,&
         pot(psi%ispot),psir,epot,confdata=psi%confdata)

    !ODP treatment (valid only for the nlr=1 case)
    if (ipotmethod==1) then !Exact Exchange
       ispot=1+npoints*(psi%nspin+psi%iorbp-1)
       !add to the psir function the part of the potential coming from the exact exchange
       call axpy(npoints,exctXcoeff,pot(ispot),1,psir(1,1),1)
    else if (ipotmethod == 2) then !PZ scheme
       !subtract the sic potential from the vpsi function
       call axpy(npoints*psi%nspinor,-alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
       !add the SIC correction to the potential energy
       epot=epot-alphaSIC*eSICi
       !accumulate the Double-Counted SIC energy
       !>>>done outside eSIC_DC=eSIC_DC+alphaSIC*eSIC_DCi
    else if (ipotmethod == 3) then !NK scheme
       !add the sic potential from the vpsi function
       call axpy(npoints*psi%nspinor,alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
       epot=epot+alphaSIC*eSICi
       !accumulate the Double-Counted SIC energy
!!!!done eSIC_DC=eSIC_DC+alphaSIC*psi%kwgt*psi%occup*eSICi
       !eSICi=psi%kwgt*psi%occup*eSICi
       eSIC_DCi=psi%kwgt*psi%occup*eSICi
    end if
    call isf_to_daub_kinetic(hgrids(1),hgrids(2),hgrids(3),&
         psi%kpoint(1),psi%kpoint(2),psi%kpoint(3),psi%nspinor,psi%lr,wrk_lh,&
         psir(1,1),hpsi,ekin)

  end subroutine local_hamiltonian_ket

  subroutine precondition_ket(ncong,confdata,ncplx,hgrids,kpoint,lr,eval,eval_zero,hpsi,scpr,&
       lin_prec_conv_work,lin_prec_work)
    use module_base,only: bigdft_mpi,nrm2
    use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
    implicit none
    integer, intent(in) :: ncong,ncplx
    real(gp), intent(in) :: eval,eval_zero
    type(confpot_data), intent(in) :: confdata
    type(locreg_descriptors), intent(in) :: lr
    real(gp), dimension(3), intent(in) :: hgrids,kpoint
    real(wp), intent(out) :: scpr
    real(wp), dimension(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)), intent(in) :: hpsi
    !> convolution workarrays for the linear case
    type(workarrays_quartic_convolutions), intent(inout), optional :: lin_prec_conv_work
    type(workarr_precond) ,intent(inout), optional :: lin_prec_work
    !local variables
    real(wp) :: cprecr

    !the nrm2 function can be replaced here by ddot
    scpr=nrm2(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),hpsi(1),1)
    if (scpr == 0.0_wp) return

    call cprecr_from_eval(lr%geocode,eval_zero,eval,cprecr)
    !cases with no CG iterations, diagonal preconditioning
    !for Free BC it is incorporated in the standard procedure
    if (ncong == 0 .and. lr%geocode /= 'F') then
       select case(lr%geocode)
       case('F')
       case('S')
          call prec_fft_slab(lr%d%n1,lr%d%n2,lr%d%n3, &
               lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
               lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
               cprecr,hgrids(1),hgrids(2),hgrids(3),hpsi)
       case('P')
          call prec_fft(lr%d%n1,lr%d%n2,lr%d%n3, &
               lr%wfd%nseg_c,lr%wfd%nvctr_c,&
               lr%wfd%nseg_f,lr%wfd%nvctr_f,&
               lr%wfd%keygloc,lr%wfd%keyvloc, &
               cprecr,hgrids(1),hgrids(2),hgrids(3),hpsi)
       end select
    else !normal preconditioner
       !case active only in the linear scaling case
       if(confdata%prefac > 0.0_gp .or. confdata%potorder > 0)then
          if (.not.present(lin_prec_conv_work)) then
             call f_err_throw("linear_precond_convol_workarrays must be present when calling the linear preconditioner", &
                  err_name='BIGDFT_RUNTIME_ERROR')
          end if
          if (.not.present(lin_prec_work)) then
             call f_err_throw("linear_precond_workarrays must be present when calling the linear preconditioner", &
                  err_name='BIGDFT_RUNTIME_ERROR')
          end if
          call solvePrecondEquation(bigdft_mpi%iproc,bigdft_mpi%nproc,lr,ncplx,ncong,&
               cprecr,&
               hgrids(1),hgrids(2),hgrids(3),&
               kpoint(1),kpoint(2),kpoint(3),&
               hpsi,&
               lr%locregCenter,&
               confdata%prefac,&
               confdata%potorder,&
               lin_prec_conv_work,lin_prec_work)
       else
          call precondition_residue(lr,ncplx,ncong,cprecr,&
               hgrids(1),hgrids(2),hgrids(3),&
               kpoint(1),kpoint(2),kpoint(3),hpsi)
       end if
    end if

  end subroutine precondition_ket
  

  subroutine orbital_basis_associate(ob,orbs,Lzd,Glr,comms,confdatarr,&
       nspin,phis_wvl)
    use dynamic_memory
    implicit none
    type(orbital_basis), intent(inout) :: ob
    integer, intent(in), optional :: nspin
    type(comms_cubic), intent(in), optional, target :: comms
    type(orbitals_data), intent(in), optional, target :: orbs
    type(local_zone_descriptors), intent(in), optional, target :: Lzd
    type(locreg_descriptors), intent(in), optional, target :: Glr !< in the case where only one Lrr is needed
    type(confpot_data), dimension(:), optional, intent(in), target :: confdatarr
    real(wp), dimension(:), target, optional :: phis_wvl
    !other elements have to be added (comms etc)
    !local variables
    integer :: ilr,iorb

    !nullification
    call nullify_orbital_basis(ob)

    if (present(orbs)) then
       ob%orbs => orbs
       if (present(comms)) then
          ob%td%comms => comms
          ob%td%nspin=orbs%nspin
          if (present(nspin)) ob%td%nspin=nspin
          ob%td%ndim_ovrlp = f_malloc_ptr([1.to.ob%td%nspin, 0.to.orbs%nkpts],id='ndim_ovrlp')
          call dimension_ovrlp(ob%td%nspin,ob%orbs,ob%td%ndim_ovrlp)
       end if
    end if

    if (.not. present(orbs) .and. (present(Lzd) .or. present(Glr))) &
         call f_err_throw('orbs should be present with lzd or glr',err_name='BIGDFT_RUNTIME_ERROR')

    if (present(Lzd)) then
       allocate(ob%dd(orbs%norbp))
       do iorb=1,orbs%norbp
          ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
          ob%dd(iorb)%lr => Lzd%Llr(ilr)
       end do
    else if (present(Glr)) then
       allocate(ob%dd(orbs%norbp))
       do iorb=1,orbs%norbp
          ob%dd(iorb)%lr => Glr
       end do
    end if

    if (present(confdatarr))  ob%confdatarr=>confdatarr

    if (present(phis_wvl)) ob%phis_wvl => phis_wvl

    !before using any iterator whatsoever, let us probe
    !its behaviour
    call probe_iterator(ob)

  end subroutine orbital_basis_associate

  subroutine orbital_basis_release(ob)
    use dynamic_memory
    implicit none
    type(orbital_basis), intent(inout) :: ob

    !nullification and reference counting (when available)
    if (associated(ob%dd)) then
       deallocate(ob%dd)
       nullify(ob%dd)
    end if
    call f_free_ptr(ob%td%ndim_ovrlp)
    call nullify_orbital_basis(ob)
  end subroutine orbital_basis_release


!!>  !here we should prepare the API to iterate on the transposed orbitals
!!>  !the orbitals should be arranged in terms of the quantum number
!!>  alag = f_malloc0(ob%nmatrix,id='alag')
!!>  qn=subspace_iterator(ob)
!!>  do while(subspace_next(qn))
!!>     psi_ptr=>ob_qn_map(psi,qn)
!!>     hpsi_ptr=>ob_qn_map(hpsi,qn)
!!>     lambda_ptr=>ob_qn_matrix_map(alag,qn)
!!>
!!>     call subspace_matrix(symm,psi_ptr,hpsi_ptr,&
!!>          qn%ncmplx,qn%nvctr,qn%norb,lambda_ptr)
!!>
!!>  end do
!!>
!!>  !allreduce
!!>  qn=subspace_iterator(ob)
!!>  do while(subspace_next(qn))
!!>     psi_ptr=>ob_qn_map(psi,qn)
!!>     hpsi_ptr=>ob_qn_map(hpsi,qn)
!!>
!!>     lambda_ptr=>ob_qn_matrix_map(alag,qn)
!!>
!!>     call lagrange_multiplier(symm,qn%occup_ptr,&
!!>          qn%ncmplx,qn%norb,lambda_ptr,trace)
!!>
!!>      if (qn%accumulate) then
!!>         occ=qn%kwgt*real(3-orbs%nspin,gp)
!!>         if (nspinor == 4) occ=qn%kwgt
!!>         scprsum=scprsum+occ*trace
!!>      end if
!!>
!!>      call subspace_update(qn%ncmplx,qn%nvctr,qn%norb,&
!!>           hpsi_ptr,lambda_ptr,psi_ptr)
!!>
!!>
!!>  end do
!!>

!!$
!!$  it=orbital_basis_iter(orb_basis,onatom='Pb')
!!$  do while(iter_next(it))
!!$
!!$  end do
!!$
!!$
!!$
!!$  call f_pointer_alias(start,end,src=phis_wvl,dest=phi_wvl,fallback=phi_add)
!!$
!!$  phi_wvl => phis_wvl(istart:iend)
!!$
!!$
!!$  subroutine overlap_matrix(phi,chi)
!!$ type(orbitals_basis), intent(inout) :: phi,chi
!!$
!!$ associated(phi%td%comms,target=chi%td%comms)
!!$
!!$  subroutine orthogonalize(orb_basis)

!!$  subroutine find_fragment_electrons(astruct_frag,nelpsp,astruct_glob)
!!$    use module_atoms
!!$    !iterate over the atoms of the fragment
!!$    !iterate above atoms
!!$    it=atoms_iter(astruct_frg)
!!$    !python metod
!!$    nelec=0
!!$    do while(atoms_iter_next(it))
!!$       it_glob=atoms_iter(astruct_glob,ontypes=.true.)
!!$       do while(atoms_iter_next(it_glob))
!!$          if (it_glob%name == it%name) then
!!$             nelec=nelec+nelpsp(it_glob%ityp)
!!$             exit
!!$          end if
!!$       end do
!!$    end do
!!$
!!$  end subroutine find_fragment_electrons
!!$
!!$  subroutine find_nn(astruct_frg,rxyz0,rxyz_nn)
!!$    dist = f_malloc(astruct_frg%nat,id='dist')
!!$    ipiv = f_malloc(astruct_frg%nat,id='ipiv')
!!$
!!$
!!$
!!$    !iterate over the atoms of the fragment
!!$    !iterate above atoms
!!$    it=atoms_iter(astruct_frg)
!!$    !python metod
!!$    do while(atoms_iter_next(it))
!!$       dist(it%iat)=-sqrt((it%rxyz-rxyz0)**2)
!!$    end do
!!$    ! sort atoms into neighbour order
!!$    call sort_positions(ref_frags(ifrag_ref)%astruct_frg%nat,dist,ipiv)
!!$
!!$    do iat=1,min(4,ref_frags(ifrag_ref)%astruct_frg%nat)
!!$       rxyz4_ref(:,iat)=rxyz_ref(:,ipiv(iat))
!!$       rxyz4_new(:,iat)=rxyz_new(:,ipiv(iat))
!!$    end do
!!$
!!$    call f_free(ipiv)
!!$    call f_free(dist)
!!$
!!$  end subroutine find_nn
!!$
!!$  !this loop fills the data of the support functions
!!$
!!$  if (cubic) then
!!$     call find_cubic_trans(cubic_tr)
!!$  end if
!!$
!!$  it=orbital_basis_iterator(orbs_basis)
!!$  do while(iter_is_valid(it))
!!$
!!$     !find transformation to apply to SF
!!$     if (linear) then
!!$        if (infrag) then
!!$           rxyz_ref
!!$           nat_frag
!!$        else
!!$           call find_nn(...,it%locregCenter,rxyz_ref)
!!$           nat_frag=min(4,astruct%nat)
!!$        end if
!!$        call find_frag_trans(nat_frag,rxyz_ref,rxyz_new,frag_trans)
!!$     else
!!$        frag_trans= cubic_tr
!!$     end if
!!$

!!$     if (disk) then
!!$        call readold()
!!$     end if
!!$
!!$     if (restart) then
!!$        call reformat_support_function(reformat_strategy,frag_trans,lr_old(it%ialpha),psi_old(it%ialpha),&
!!$             it%lr,it%sf)
!!$     end if
!!$
!!$  end do

  subroutine test_iterator(ob)
    use yaml_output
    implicit none
    type(orbital_basis), intent(in) :: ob
    type(ket) :: it

    !iterate over all the orbitals
    !iterate over the orbital_basis
    it=orbital_basis_iterator(ob)
    do while(ket_next(it))
       call yaml_map('Locreg, orbs',[it%ilr,it%iorb,it%iorbp])
       call yaml_map('associated lr',associated(it%lr))
    end do

    !iterate over the orbital_basis
    it=orbital_basis_iterator(ob)
    loop_lr: do while(ket_next_locreg(it))
       call yaml_map('Locreg2',it%ilr)
       call yaml_map('associated lr',associated(it%lr))
       loop_iorb: do while(ket_next(it,ilr=it%ilr))
          call yaml_map('Locreg2, orbs',[it%ilr,it%iorb,it%iorbp])
       end do loop_iorb
    end do loop_lr
  end subroutine test_iterator

  !>confirm the iterator coincides with the provided information
  function verify_subspace(ss,ispin,ikpt,ikptp,nvctr,&
       ncplx,norb,ispsi,ise,kwgt) result(ok)
    implicit none
    type(subspace), intent(in) :: ss
    integer, intent(in) :: ispin,ikpt,ikptp,nvctr,ncplx,norb,ispsi,ise
    real(gp), intent(in) :: kwgt
    logical :: ok
    !local variables
    logical, parameter :: debug_flag=.false.

    ok=.true.

    ok=ss%ispin==ispin
    if (debug_flag) print *,'ispin',ss%ispin,ispin,ok
    if (.not. ok) return
    ok=ss%ikpt==ikpt
    if (debug_flag) print *,'ikpt',ss%ikpt,ikpt,ok
    if (.not. ok) return
    ok=ss%ikptp==ikptp
    if (debug_flag) print *,'ikptp',ss%ikptp,ikptp,ok
    if (.not. ok) return
    ok=ss%nvctr==nvctr
    if (debug_flag) print *,'nvctr',ss%nvctr,nvctr,ok
    if (.not. ok) return
    ok=ss%ncplx==ncplx
    if (debug_flag) print *,'ncplx',ss%ncplx,ncplx,ok
    if (.not. ok) return
    ok=ss%norb==norb
    if (debug_flag) print *,'norb',ss%norb,norb,ok
    if (.not. ok) return
    ok=ss%ispsi==ispsi
    if (debug_flag) print *,'ispsi',ss%ispsi,ispsi,ok
    if (.not. ok) return
    ok=ss%ise==ise
    if (debug_flag) print *,'ise',ss%ise,ise,ok
    if (.not. ok) return
    ok=ss%kwgt==kwgt
    if (debug_flag) print *,'kwgt',ss%kwgt,kwgt,ok
    if (.not. ok) return
  end function verify_subspace

  !>confirm the iterator coincides with the provided information
  function verify_iterator(it,iorb,ilr,ikpt,ispsi) result(ok)
    implicit none
    type(ket), intent(in) :: it
    integer, intent(in), optional :: iorb,ilr,ispsi,ikpt
    logical :: ok
    !local variables
    logical, parameter :: debug_flag=.false.

    ok=.true.

    if (present(iorb)) then
       ok=it%iorbp==iorb
       if (debug_flag) print *,'iorb',it%iorbp,iorb,ok
    end if
    if (.not. ok) return
    if (present(ilr)) then
       ok= it%ilr==ilr
       if (debug_flag) print *,'ilr',it%ilr,ilr,ok
    end if
    if (.not. ok) return
    if (present(ikpt)) then
       ok= it%ikpt==ikpt
       if (debug_flag) print *,'ikpt',it%ikpt,ikpt,ok
    end if
    if (.not. ok) return
    if (present(ispsi)) then
       ok= it%ispsi == ispsi
       if (debug_flag) print *,'ispsi',it%ispsi,ispsi,ok
    end if
  end function verify_iterator


  !> verify that the iterator performs the same operations of the (nested) loops
  subroutine probe_iterator(ob)
    use yaml_strings
    use dynamic_memory
    use module_base, only: bigdft_mpi
    implicit none
    type(orbital_basis), intent(in) :: ob
    type(ket) :: it
    type(subspace) :: ss
    !local variables
    !logical, parameter :: debug_flag=.true.
    logical :: doso,increase
    integer :: iorb,ikpt,ilr,isorb,ieorb,ispsi,ispsi_k,nsp,nkptp,ikptp
    integer :: ispin,ise,nvctrp,norb,norbs,ncomp,nspinor,ncomponents,ncomplex
    logical, dimension(:), allocatable :: totreat

    totreat=f_malloc(ob%orbs%norbp,id='totreat')

    !first test the different API of the iterators

    !loop over all the orbitals regardless of the order
    !this should always work
    totreat=.true.
    it=orbital_basis_iterator(ob)
    do while(ket_next(it))
       if (.not. totreat(it%iorbp)) &
            call f_err_throw('Error for iterator (1), iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
       totreat(it%iorbp)=.false.
    end do
!print *,'totreat',totreat
    if (any(totreat)) &
         call f_err_throw('Error for iterator (1), not all orbitals treated',&
         err_name='BIGDFT_RUNTIME_ERROR')

    !other version of the loop, by locreg
    !this would only work if there is only one k-point
    if (ob%orbs%nkpts ==1 .or. maxval(ob%orbs%inwhichlocreg)==1) then
       totreat=.true.
       it=orbital_basis_iterator(ob)
       do while(ket_next_locreg(it))
          !this would cycle on the locreg of this kpoint
          do while(ket_next(it,ilr=it%ilr))
             if (.not. totreat(it%iorbp)) &
                  call f_err_throw('Error for iterator (2), iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
             totreat(it%iorbp)=.false.
          end do
       end do
       if (any(totreat)) &
            call f_err_throw('Error for iterator (2), not all orbitals treated',&
            err_name='BIGDFT_RUNTIME_ERROR')
    end if

    !this also should always work
    totreat=.true.
    it=orbital_basis_iterator(ob)
    do while(ket_next_kpt(it))
       !this loop would cycle into the orbitals of this kpoint
       !regardless of the locreg
       do while(ket_next(it,ikpt=it%ikpt))
!print *,'iorb',it%iorbp,it%ilr,it%ikpt
          if (.not. totreat(it%iorbp)) &
               call f_err_throw('Error for iterator (3), iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
          totreat(it%iorbp)=.false.
       end do
    end do
!print *,'totreat',totreat
    if (any(totreat)) &
         call f_err_throw('Error for iterator (3), not all orbitals treated',&
         err_name='BIGDFT_RUNTIME_ERROR')


    !this is the most explicit case and should always be OK
    totreat=.true.
    it=orbital_basis_iterator(ob)
    do while(ket_next_kpt(it))
       !this loop would cycle into the locreg of this kpoint
       do while (ket_next_locreg(it,ikpt=it%ikpt))
          !this loop would cycle into the orbitals of this kpoint
          !and of this locreg
          do while(ket_next(it,ikpt=it%ikpt,ilr=it%ilr))
             if (.not. totreat(it%iorbp)) &
                  call f_err_throw('Error for iterator (4), iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
             totreat(it%iorbp)=.false.
          end do
       end do
    end do
!print *,'totreat',totreat
    if (any(totreat)) &
         call f_err_throw('Error for iterator (4), not all orbitals treated',&
         err_name='BIGDFT_RUNTIME_ERROR')

    nkptp=0
    if (ob%orbs%norbp > 0) nkptp=ob%orbs%iokpt(ob%orbs%norbp)-ob%orbs%iokpt(1)+1
    totreat=.true.
    it=orbital_basis_iterator(ob)
    do while(ket_next_kpt(it))
       ikptp=it%ikpt-ob%orbs%iokpt(1)+1
       if (.not. totreat(ikptp)) &
            call f_err_throw('Error for iterator (5), iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
       totreat(ikptp)=.false.
    end do
    !print *,'totreat',totreat
    if (any(totreat(1:nkptp))) &
         call f_err_throw('Error for iterator (5), not all orbitals treated',&
         err_name='BIGDFT_RUNTIME_ERROR')

    !then there is also the case without the orbitals, iterate only on k_points

    !reverse testing

    totreat=.true.
    !iterate over all the orbitals
    it=orbital_basis_iterator(ob)
    !first test the behaviour of the iterator in the normal loop, unordered
    do iorb=1,ob%orbs%norbp
       increase=ket_next(it)
       !check that all the
       if (.not. totreat(it%iorbp) .and. increase) then
          call f_err_throw('Error for iterator, iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
          return
       end if
       totreat(it%iorbp)=.false.
    end do
    if (ket_next(it) .or. any(totreat)) then
       call f_err_throw('Error for iterator, still valid at the end of the number of iterations',&
            err_name='BIGDFT_RUNTIME_ERROR')
       return
    end if


    !then test the iterator by ordering the wavefunctions in terms of localisation regions, for only one kpoint
    if ((ob%orbs%nkpts == 1 .or. maxval(ob%orbs%inwhichlocreg)==1) .and. associated(ob%dd)) then
       isorb=1
       ieorb=ob%orbs%norbp
       ispsi_k=1
       totreat=.true.
       it=orbital_basis_iterator(ob)
       loop_lr: do ilr=1,it%ilr_max !here the traditional loop was up to nlr, but now we do not know it anymore
          !do something only if at least one of the orbitals lives in the ilr
          doso=.false.
          do iorb=isorb,ieorb
             doso = (ob%orbs%inwhichlocreg(iorb+ob%orbs%isorb) == ilr)
             if (doso) exit
          end do
          if (.not. doso) cycle loop_lr
          increase=ket_next_locreg(it)
          if (.not. (verify_iterator(it,ilr=ilr) .and. increase)) then
             call f_err_throw('Error for iterator, ilr='+ilr,err_name='BIGDFT_RUNTIME_ERROR')
             return
          end if
          ispsi=ispsi_k
          !loop for the wavefunctions
          do iorb=isorb,ieorb
             if (ob%orbs%inwhichlocreg(iorb+ob%orbs%isorb) /= ilr) then
                !increase ispsi to meet orbital index
                ispsi=ispsi+(ob%dd(iorb)%lr%wfd%nvctr_c+&
                     7*ob%dd(iorb)%lr%wfd%nvctr_f)*ob%orbs%nspinor
                cycle
             end if
             increase=ket_next(it,ilr=it%ilr)
             if (.not. (verify_iterator(it,iorb=iorb,ispsi=ispsi) .and. increase ) .and. .not. totreat(it%iorbp)) then
                call f_err_throw('Error for iterator, ilr='+ilr//' iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
                return
             end if
             totreat(it%iorbp)=.false.
             ispsi=ispsi+(ob%dd(iorb)%lr%wfd%nvctr_c+&
                  7*ob%dd(iorb)%lr%wfd%nvctr_f)*ob%orbs%nspinor
          end do
          if (ket_next(it,ilr=it%ilr)) then
             call f_err_throw('Error for iterator, still valid at the end of the number of iterations, ilr='+ilr,&
                  err_name='BIGDFT_RUNTIME_ERROR')
             return
          end if
          ispsi_k=ispsi
       end do loop_lr
       if (ket_next_locreg(it) .or. any(totreat)) then
          call f_err_throw('Error for iterator, still valid at the end of the number of ilr iterations, ilr='+ilr,&
               err_name='BIGDFT_RUNTIME_ERROR')
          return
       end if
    end if

    !let us now probe the behaviour of the iterator with respect to the
    !kpoints only
    if (associated(ob%dd)) then
       if (ob%orbs%norbp >0)ikpt=ob%orbs%iokpt(1)
       ispsi_k=1
       totreat=.true.
       it=orbital_basis_iterator(ob)
       loop_kpt: do
          !specialized treatment
          if (ob%orbs%norbp ==0) then
             increase=ket_next_kpt(it)
             increase=ket_next(it,ikpt=it%ikpt)
             if (increase) then
                call f_err_throw(&
                     'Error for iterator, still valid at the end of the zero number of kpt iterations, ikpt='+ikpt,&
                     err_name='BIGDFT_RUNTIME_ERROR')
                return
             end if
             exit loop_kpt
          end if
          call orbs_in_kpt(ikpt,ob%orbs,isorb,ieorb,nsp)

          increase=ket_next_kpt(it)
          if (.not. (verify_iterator(it,ikpt=ikpt) .and. increase)) then
             call f_err_throw('Error for iterator, ikpt='+ikpt,err_name='BIGDFT_RUNTIME_ERROR')
             return
          end if

          !localisation regions loop
          loop_lrk: do ilr=1,it%ilr_max
             !do something only if at least one of the orbitals lives in the ilr
             doso=.false.
             do iorb=isorb,ieorb
                doso = (ob%orbs%inwhichlocreg(iorb+ob%orbs%isorb) == ilr)
                if (doso) exit
             end do
             if (.not. doso) cycle loop_lrk
             ispsi=ispsi_k
             do iorb=isorb,ieorb
                if (ob%orbs%inwhichlocreg(iorb+ob%orbs%isorb) /= ilr) then
                   !increase ispsi to meet orbital index
                   ispsi=ispsi+(ob%dd(iorb)%lr%wfd%nvctr_c+&
                        7*ob%dd(iorb)%lr%wfd%nvctr_f)*nsp
                   cycle
                end if
                increase=ket_next(it,ikpt=it%ikpt)
                if (.not. (verify_iterator(it,iorb=iorb,ispsi=ispsi) .and. increase ) .and. .not. totreat(it%iorbp)) then
                   call f_err_throw('Error for iterator, ikpt='+ikpt//' iorb='+iorb,err_name='BIGDFT_RUNTIME_ERROR')
                   return
                end if
                totreat(it%iorbp)=.false.

                ispsi=ispsi+(ob%dd(iorb)%lr%wfd%nvctr_c+&
                     7*ob%dd(iorb)%lr%wfd%nvctr_f)*nsp
             end do

          end do loop_lrk

          !last k-point has been treated
          if (ieorb == ob%orbs%norbp) exit loop_kpt

          ikpt=ikpt+1
          ispsi_k=ispsi

       end do loop_kpt
       if (ket_next_kpt(it) .or. any(totreat)) then
          call f_err_throw('Error for iterator, still valid at the end of the number of kpt iterations, ikpt='+ikpt,&
               err_name='BIGDFT_RUNTIME_ERROR')
          return
       end if
    end if
    call f_free(totreat)

    !When it is possible, perform also the check for the transposed operations
    if (.not. associated(ob%td%comms)) return

    !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case

    totreat=f_malloc(ob%orbs%nkptsp*ob%td%nspin,id='totreat')
    totreat=.true.

    ss=subspace_iterator(ob)
    !old format, check the values
    ispsi=1
    do ikptp=1,ob%orbs%nkptsp
       ikpt=ob%orbs%iskpts+ikptp
       do ispin=1,ob%td%nspin
          if (ispin==1) ise=0
          call orbitals_and_components(bigdft_mpi%iproc,ikpt,ispin,ob%orbs,ob%td%comms,&
               nvctrp,norb,norbs,ncomp,nspinor)
          totreat(ispin+(ikptp-1)*ob%td%nspin)=.false.
          if (nvctrp == 0) cycle

          ncomplex=1
          ncomponents=ncomp*nvctrp
          if (nspinor/=1) ncomplex=2

          increase=subspace_next(ss)
          if (.not. (verify_subspace(ss,ispin,ikpt,ikptp,nvctrp,&
               ncomplex,norb,ispsi,ise,ob%orbs%kwgts(ikpt))) &
               .and. increase) then
               call f_err_throw('Error for direct subspace iterator, ikpt='+&
                ikpt+', ispin='+ispin,err_name='BIGDFT_RUNTIME_ERROR')
               return
          end if

          ise=norb
          ispsi=ispsi+nvctrp*norb*nspinor
       end do
    end do
    increase=subspace_next(ss)
    if (increase .or. any(totreat)) then
      call f_err_throw('Error for direct subspace iterator,'//&
      ' still valid at the end of the number of kpt iterations',&
           err_name='BIGDFT_RUNTIME_ERROR')
      return
    end if
    call f_free(totreat)

  end subroutine probe_iterator


end module orbitalbasis
