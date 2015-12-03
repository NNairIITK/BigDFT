!> @file
!! Datatypes and associated methods relative to the localization regions
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
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
  implicit none
  private

  type, public :: transposed_descriptor
     integer :: nspin
     integer :: nup,ndown
!     type(k_point_iterator) :: kit_p !< iterator on k-points belonging to the processor
     type(comms_cubic), pointer :: comms           !< communication objects for the cubic approach (can be checked to be
     type(comms_linear), pointer :: collcom        !< describes collective communication
     type(comms_linear), pointer :: collcom_sr     !< describes collective communication for the calculation of th
  end type transposed_descriptor

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
     integer :: nlrp,ilr,ilr_max,iorbp
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
  public :: orbital_basis_associate,orbital_basis_release,test_iterator

contains

  pure subroutine nullify_orbital_basis(ob)
    implicit none
    type(orbital_basis), intent(out) :: ob
    nullify(ob%dd)
    nullify(ob%orbs)
    !type(transposed_descriptor) :: td
    nullify(ob%confdatarr)
    nullify(ob%phis_wvl)
  end subroutine nullify_orbital_basis

  

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
    it%nlrp = size(ob%dd)
    !zero value
    it%ilr=minval(ob%orbs%inwhichlocreg(ob%orbs%isorb+1:ob%orbs%isorb+ob%orbs%norbp))-1
    !last value
    it%ilr_max=maxval(ob%orbs%inwhichlocreg(ob%orbs%isorb+1:ob%orbs%isorb+ob%orbs%norbp))
    !start orbital
    it%iorbp=0

  end function orbital_basis_iterator

  function ket_next_locreg(it) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    logical :: ok

    ok=ket_is_valid(it)
    if (.not. ok) return

    find_next_lr: do while(it%ilr <= it%ilr_max)
       it%ilr=it%ilr+1
       if (dosome(it, it%ilr)) exit find_next_lr
    end do find_next_lr
    ok= it%ilr <= it%ilr_max
    if (.not. ok) then
       call nullify_ket(it)
       return
    end if

    if (ok) call update_ket(it)

    it%iorbp=0

  end function ket_next_locreg

  function ket_is_valid(it) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    logical :: ok

    ok=.true.
    if (.not. associated(it%ob)) then
       call nullify_ket(it)
       ok=.false.
       return
    end if
    if (it%ilr > it%ilr_max) then
       ok =.false.
       call nullify_ket(it)
       return
    end if

  end function ket_is_valid

  function ket_next(it,ilr) result(ok)
    implicit none
    type(ket), intent(inout) :: it
    integer, intent(in), optional :: ilr
    logical :: ok
    !local variables
    integer :: ilr_tmp

    ok=ket_is_valid(it)
    if (.not. ok) return

    if (present(ilr)) then
       ok=dosome(it,ilr)
    else
    !increment the localization region
    !first increment the orbital, otherwise the locreg
       ilr_tmp=it%ilr
       find_next_ilr: do while(.not. dosome(it,ilr_tmp) .and. ilr_tmp <= it%ilr_max)
          ilr_tmp=ilr_tmp+1
       end do find_next_ilr
       it%ilr=ilr_tmp
       ok= it%ilr <= it%ilr_max 
       if (.not. ok) then
          call nullify_ket(it)
          return
       end if
    end if

    !at this point the iorbp and ilr are determined, the iterator can be updated
    if (ok) call update_ket(it)
    
  end function ket_next

  !> find next valid iorb in this lr
  !! put iorbp to zero if ilr_tmp is not found
  function dosome(it,ilr_tmp)
    implicit none
    type(ket), intent(inout) :: it
    integer, intent(in) :: ilr_tmp
    logical :: dosome
    dosome=.false.
    find_iorb: do while(.not. dosome .and. it%iorbp < it%ob%orbs%norbp)
       it%iorbp=it%iorbp+1
       !check if this localisation region is used by one of the orbitals
       dosome = (it%ob%orbs%inwhichlocreg(it%iorbp+it%ob%orbs%isorb) == ilr_tmp)
    end do find_iorb
    if (.not. dosome) it%iorbp=0
  end function dosome


  subroutine update_ket(k)
    implicit none
    type(ket), intent(inout) :: k
    !local variables
    integer :: ikpt,ispsi,iorbq,ilr_orb,nvctr
    !the orbital id
    k%iorb=k%ob%orbs%isorb+k%iorbp
    k%nspin=k%ob%orbs%nspin
    k%nspinor=k%ob%orbs%nspinor
    !k-point, spin and confinement
    ikpt=k%ob%orbs%iokpt(k%iorbp)
    k%kpoint=k%ob%orbs%kpts(:,ikpt)
    k%kwgt=k%ob%orbs%kwgts(ikpt)
    k%occup=k%ob%orbs%occup(k%iorb)
    k%spinval=k%ob%orbs%spinsgn(k%iorb)
    if (associated(k%ob%confdatarr)) k%confdata=k%ob%confdatarr(k%iorbp)
    !shifts metadata
    k%ispot=k%ob%orbs%ispot(k%iorbp)
    !find the psi shift for the association
    k%ispsi=1
    do iorbq=1,k%iorbp-1
       nvctr=k%ob%dd(iorbq)%lr%wfd%nvctr_c+7*k%ob%dd(iorbq)%lr%wfd%nvctr_f
       k%ispsi=k%ispsi+nvctr*k%nspinor
    end do
    k%lr=>k%ob%dd(k%iorbp)%lr
    k%nphidim=(k%lr%wfd%nvctr_c+7*k%lr%wfd%nvctr_f)*k%nspinor

    k%phi_wvl=>ob_ket_map(k%ob%phis_wvl,k)
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


  !> this function gives the number of componenets
  !! if the ket is in non-colinear spin description, this value is four
  !! otherwise it counts the number of complex componenet of the key
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
  
  subroutine orbital_basis_associate(ob,orbs,Lzd,confdatarr,phis_wvl)
    implicit none
    type(orbital_basis), intent(inout) :: ob
    type(orbitals_data), intent(in), optional, target :: orbs
    type(local_zone_descriptors), intent(in), optional, target :: Lzd
    type(confpot_data), dimension(:), optional, intent(in), target :: confdatarr
    real(wp), dimension(:), target, optional :: phis_wvl
    !other elements have to be added (comms etc)
    !local variables
    integer :: ilr,iorb

    !nullification
    call nullify_orbital_basis(ob)

    if (present(orbs)) ob%orbs => orbs
    
    if (present(Lzd)) then
       if (.not. present(orbs)) call f_err_throw('orbs should be present with lzd',err_name='BIGDFT_RUNTIME_ERROR')
       allocate(ob%dd(orbs%norbp))
       do iorb=1,orbs%norbp
          ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
          ob%dd(iorb)%lr => Lzd%Llr(ilr)
       end do
    end if

    if (present(confdatarr))  ob%confdatarr=>confdatarr

    if (present(phis_wvl)) ob%phis_wvl => phis_wvl
  end subroutine orbital_basis_associate

  subroutine orbital_basis_release(ob)
    implicit none
    type(orbital_basis), intent(inout) :: ob

    !nullification and reference counting (when available)
    if (associated(ob%dd)) then
       deallocate(ob%dd)
       nullify(ob%dd)
    end if
    call nullify_orbital_basis(ob)
  end subroutine orbital_basis_release
  

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




end module orbitalbasis
