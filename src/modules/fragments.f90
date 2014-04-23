!> @file
!!  Module to handle the fragments of a system
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Module which defines some important structures and methods to manipulate embedding systems
module module_fragments
  use module_base, only: gp,wp
  use module_types
  use dynamic_memory
  use module_atoms
  implicit none

  private

   interface

      subroutine reorthonormalize_coeff(iproc, nproc, norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, basis_orbs, &
                 basis_overlap, basis_overlap_mat, coeff, orbs)
        use module_base
        use module_types
        use sparsematrix_base, only: sparse_matrix
        implicit none
        integer, intent(in) :: iproc, nproc, norb
        integer, intent(in) :: blocksize_dsyev, blocksize_pdgemm, inversion_method
        type(orbitals_data), intent(in) :: basis_orbs   !number of basis functions
        type(orbitals_data), intent(in) :: orbs   !Kohn-Sham orbitals that will be orthonormalized and their parallel distribution
        type(sparse_matrix),intent(inout) :: basis_overlap
        type(matrices),intent(inout) :: basis_overlap_mat
        real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff
      end subroutine reorthonormalize_coeff

      subroutine write_eigenvalues_data(etol,orbs,mom_vec)
        use module_base
        use module_types
        implicit none
        real(gp), intent(in) :: etol
        type(orbitals_data), intent(in) :: orbs
        real(gp), dimension(:,:,:), intent(in), pointer :: mom_vec
      end subroutine write_eigenvalues_data

   end interface


  !> information about the basis set metadata shared between cubic and ASF approaches
  type, public :: minimal_orbitals_data
     integer :: norb          !< Total number of orbitals per k point
     integer :: norbp         !< Total number of orbitals for the given processors
     integer :: isorb         !< Total number of orbitals for the given processors
     integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
     integer, dimension(:), pointer :: isorb_par,ispot
     integer, dimension(:,:), pointer :: norb_par
  end type minimal_orbitals_data

  type, public :: phi_array
     real(wp), dimension(:,:,:,:,:,:), pointer :: psig
  end type phi_array

  type, public :: fragment_basis
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(local_zone_descriptors) :: Lzd
     type(minimal_orbitals_data) :: forbs
     real(wp), dimension(:), pointer :: psi_full
     type(phi_array), dimension(:), pointer :: phi
     real(wp), dimension(:), pointer :: density
  end type fragment_basis

  type, public :: fragment_orbitals
     integer :: norb, nphi

  end type fragment_orbitals

  !> Defines the minimal information to identify a system building block
  type, public :: system_fragment
     integer :: nat_env !< environment atoms which complete fragment specifications
     real(gp), dimension(:,:), pointer :: rxyz_env !< position of atoms in environment (AU), external reference frame
     type(atomic_structure) :: astruct_frg !< Number of atoms, positions, atom type etc for fragment
     type(fragment_basis) :: fbasis !< fragment basis, associated only if coherent with positions, pointer - do we really want this to be a pointer?
     ! add coeffs and or kernel
     integer :: nelec
     real(gp), dimension(:,:), pointer :: coeff
     real(gp), dimension(:,:), pointer :: kernel
     real(gp), dimension(:), pointer :: eval
  end type system_fragment

  !> Contains the rotation and translation (possibly deformation) which have to be applied to a given fragment
  type, public :: fragment_transformation
     !real(gp), dimension(3) ::  !< translation of fragment center
     real(gp), dimension(3) :: rot_center_new !< positions of the centers
     real(gp), dimension(3) :: rot_center !< center of rotation in original coordinates (might be fragment center or not)
     real(gp), dimension(3) :: rot_axis !< unit rotation axis (should have modulus one)
     real(gp) :: theta !< angle of rotation 
  end type fragment_transformation

  !public operator(*)

  public :: fragment_null, fragment_free, init_fragments, minimal_orbitals_data_null, rotate_vector
  public :: frag_center, find_frag_trans, calculate_fragment_density, fragment_coeffs_to_kernel

contains

  ! nned to check somewhere if fragment linking is consistent - initialize linear from file?!
  ! initializes reference fragments (already nullified), if it isn't a fragment calculation sets to appropriate dummy values
  ! ignoring environment for now
  subroutine init_fragments(in,orbs,astruct,ref_frags)
    use module_types
    implicit none
    type(input_variables), intent(in) :: in
    type(orbitals_data), intent(in) :: orbs ! orbitals of full system, needed to set 'dummy' values
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system, needed to set 'dummy' values
    type(system_fragment), dimension(in%frag%nfrag_ref), intent(inout) :: ref_frags

    ! local variables
    integer :: ifrag

    if (in%lin%fragment_calculation) then
        ! read fragment posinps and initialize fragment, except for psi and lzds
        do ifrag=1,in%frag%nfrag_ref
           call init_fragment_from_file(ref_frags(ifrag),trim(in%dir_output)//trim(in%frag%label(ifrag)),&
                in,astruct)
        end do

        ! check that fragments are sensible, i.e. correct number of atoms, atom types etc.
        call check_fragments(in,ref_frags,astruct)

     ! set appropriate default values as this is not a fragment calculation
     else
        ! nullify fragment
        ref_frags(1)=fragment_null()

        ! want all components of ref_frags(1)%fbasis%forbs to point towards appropriate component of orbs of full system
        call orbs_to_min_orbs_point(orbs,ref_frags(1)%fbasis%forbs)

        ! environment and coeffs
        call fragment_allocate(ref_frags(1))

        ! rest of fbasis isn't needed (I think!) so can remain nullified
  
        ! astruct - fill in other bits later
        ref_frags(1)%astruct_frg%nat=astruct%nat

     end if

  end subroutine init_fragments


  !> Initializes all of fragment except lzd using the fragment posinp and tmb files
  subroutine init_fragment_from_file(frag,frag_name,input,astruct) ! switch this to pure if possible
    use module_types
    !use module_interfaces
    implicit none
    type(system_fragment), intent(inout) :: frag
    character(len=*), intent(in) :: frag_name
    type(input_variables), intent(in) :: input
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system


    ! nullify fragment
    frag=fragment_null()

    ! read fragment positions
    call set_astruct_from_file(frag_name(1:len(frag_name)),bigdft_mpi%iproc,frag%astruct_frg)

    ! iproc, nproc, nspinor not needed yet, add in later
    call init_minimal_orbitals_data(bigdft_mpi%iproc, bigdft_mpi%nproc, 1, input, frag%astruct_frg, &
         frag%fbasis%forbs,astruct)
    !call init_minimal_orbitals_data(iproc, nproc, nspinor, input, frag%astruct_frg, frag%fbasis%forbs)

    ! environment and coeffs
    call fragment_allocate(frag)

    ! allocate/initialize fragment basis...
    !currently orbitals are initialized via initAndUtils/init_orbitals_data_for_linear
    !which calls wavefunctions/orbitals_descriptors to assign parallel bits and superfluous stuff
    !and locreg_orbitals/assignToLocreg2 to give inwhichlocreg and onwhichatom - but we should really take onwhichatom from file?!

    !type, public :: fragment_basis
    !   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
    !   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
    !   type(local_zone_descriptors) :: Lzd
    !   type(minimal_orbitals_data) :: forbs

  end subroutine init_fragment_from_file


  ! sanity check on fragment definitions
  subroutine check_fragments(input,ref_frags,astruct)
    use module_types
    implicit none
    type(input_variables), intent(in) :: input
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system, needed to check fragments are sensible
    type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

    integer :: ifrag, ifrag_ref, tot_frag_ats, isfat, iat
    logical :: fragments_ok

    fragments_ok=.true.

    tot_frag_ats=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref = input%frag%frag_index(ifrag)
       tot_frag_ats = tot_frag_ats + ref_frags(ifrag_ref)%astruct_frg%nat
    end do

    if (tot_frag_ats /= astruct%nat) then
       fragments_ok=.false.
       write(*,*) 'Sum of atoms in fragments not equal to total atoms in system',tot_frag_ats,astruct%nat
    end if

    isfat=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)

       do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
          if (trim(ref_frags(ifrag_ref)%astruct_frg%atomnames(ref_frags(ifrag_ref)%astruct_frg%iatype(iat))) &
               /= trim(astruct%atomnames(astruct%iatype(iat+isfat)))) then
             fragments_ok=.false.
             write(*,*) 'Atom type for fragment ',ifrag,', reference fragment ',ifrag_ref,' atom number ',iat,&
                  ' does not match ',&
                  trim(ref_frags(ifrag_ref)%astruct_frg%atomnames(ref_frags(ifrag_ref)%astruct_frg%iatype(iat))),&
                  trim(astruct%atomnames(astruct%iatype(iat+isfat)))
          end if
       end do

       isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     
    end do

    if (.not. fragments_ok) stop 'Problem in fragment definitions'

  end subroutine check_fragments


    !type, public :: minimal_orbitals_data
    !   integer :: norb          !< Total number of orbitals per k point
    !   integer :: norbp         !< Total number of orbitals for the given processors
    !   integer :: isorb         !< Total number of orbitals for the given processors
    !   integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
    !   integer, dimension(:), pointer :: isorb_par,ispot
    !   integer, dimension(:,:), pointer :: norb_par
  !> just initializing norb for now, come back and do the rest later
  subroutine init_minimal_orbitals_data(iproc, nproc, nspinor, input, astruct, forbs, astruct_full)
    use module_base
    use module_types
    implicit none
  
    ! Calling arguments
    integer,intent(in) :: iproc, nproc, nspinor
    type(input_variables),intent(in) :: input
    type(atomic_structure),intent(in) :: astruct
    type(minimal_orbitals_data),intent(out) :: forbs
    type(atomic_structure), intent(in) :: astruct_full ! atomic structure of full system
  
    ! Local variables
    integer :: norb, ityp, iat, jtyp
    !integer :: norbu, norbd, nlr, ilr, iall, iorb, istat
    !integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
    !real(kind=8),dimension(:,:),allocatable :: locregCenter
    !character(len=*),parameter :: subname='init_minimal_orbitals_data'

    call timing(iproc,'init_orbs_lin ','ON')
 
    ! Count the number of basis functions.
    !allocate(norbsPerAtom(astruct%nat), stat=istat)
    !call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
    norb=0
    !nlr=0
    do iat=1,astruct%nat
       ityp=astruct%iatype(iat)
       !norbsPerAtom(iat)=input%lin%norbsPerType(ityp)
       do jtyp=1,astruct_full%ntypes
          if (astruct_full%atomnames(jtyp)==astruct%atomnames(ityp)) exit
       end do
       if (jtyp==astruct_full%ntypes+1) then
          print*, 'Error in fragment_init_orbitals, atom type ',astruct%atomnames(ityp),' does not exist in full structure'
          stop
       end if
       norb=norb+input%lin%norbsPerType(jtyp)
       !nlr=nlr+input%lin%norbsPerType(ityp)
    end do

    ! Distribute the basis functions among the processors.
    !norbu=norb
    !norbd=0
 
    forbs%norb=norb

    !call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
    !     input%nkpt, input%kpt, input%wkpt, lorbs,.true.) !simple repartition
 

    !allocate(locregCenter(3,nlr), stat=istat)
    !call memocc(istat, locregCenter, 'locregCenter', subname)
  
    !ilr=0
    !do iat=1,astruct%nat
    !    ityp=astruct%iatype(iat)
    !    do iorb=1,input%lin%norbsPerType(ityp)
    !        ilr=ilr+1
    !        locregCenter(:,ilr)=astruct%rxyz(:,iat)
    !        ! DEBUGLR write(10,*) iorb,locregCenter(:,ilr)
    !    end do
    !end do
 
    !allocate(norbsPerLocreg(nlr), stat=istat)
    !call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)
    !norbsPerLocreg=1 !should be norbsPerLocreg
    
    !iall=-product(shape(forbs%inWhichLocreg))*kind(forbs%inWhichLocreg)
    !deallocate(forbs%inWhichLocreg, stat=istat)
    !call memocc(istat, iall, 'forbs%inWhichLocreg', subname)
    !call assignToLocreg2(iproc, nproc, forbs%norb, forbs%norb_par, astruct%nat, nlr, &
    !     input%nspin, norbsPerLocreg, locregCenter, forbs%inwhichlocreg)

    !iall=-product(shape(forbs%onwhichatom))*kind(forbs%onwhichatom)
    !deallocate(forbs%onwhichatom, stat=istat)
    !call memocc(istat, iall, 'forbs%onwhichatom', subname)
    !call assignToLocreg2(iproc, nproc, forbs%norb, forbs%norb_par, astruct%nat, astruct%nat, &
    !     input%nspin, norbsPerAtom, astruct%rxyz, forbs%onwhichatom)
  

    !iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
    !deallocate(norbsPerLocreg, stat=istat)
    !call memocc(istat, iall, 'norbsPerLocreg', subname)
  
    !iall=-product(shape(locregCenter))*kind(locregCenter)
    !deallocate(locregCenter, stat=istat)
    !call memocc(istat, iall, 'locregCenter', subname)

    !iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
    !deallocate(norbsPerAtom, stat=istat)
    !call memocc(istat, iall, 'norbsPerAtom', subname)


    call timing(iproc,'init_orbs_lin ','OF')

  end subroutine init_minimal_orbitals_data

  ! copy a given full orbs structure to minimal orbs structure
  subroutine orbs_to_min_orbs_copy(orbs,forbs)
    use module_types
    implicit none
    ! Calling arguments
    type(orbitals_data),intent(in):: orbs
    type(minimal_orbitals_data),intent(inout):: forbs

    ! Local variables
    integer:: iis1, iie1, iis2, iie2, i1, i2, istat, iall
    character(len=200) :: subname

    subname='orbs_to_min_orbs_copy'

    forbs%norb = orbs%norb
    forbs%norbp = orbs%norbp
    forbs%isorb = orbs%isorb

    if(associated(forbs%norb_par)) then
        iall=-product(shape(forbs%norb_par))*kind(forbs%norb_par)
        deallocate(forbs%norb_par, stat=istat)
        call memocc(istat, iall, 'forbs%norb_par', subname)
    end if
    if(associated(orbs%norb_par)) then
        iis1=lbound(orbs%norb_par,1)
        iie1=ubound(orbs%norb_par,1)
        iis2=lbound(orbs%norb_par,2)
        iie2=ubound(orbs%norb_par,2)
        allocate(forbs%norb_par(iis1:iie1,iis2:iie2), stat=istat)
        call memocc(istat, forbs%norb_par, 'forbs%norb_par', subname)
        do i1=iis1,iie1
           do i2 = iis2,iie2
            forbs%norb_par(i1,i2) = orbs%norb_par(i1,i2)
           end do
        end do
    end if

    if(associated(forbs%inwhichlocreg)) then
        iall=-product(shape(forbs%inwhichlocreg))*kind(forbs%inwhichlocreg)
        deallocate(forbs%inwhichlocreg, stat=istat)
        call memocc(istat, iall, 'forbs%inwhichlocreg', subname)
    end if
    if(associated(orbs%inwhichlocreg)) then
        iis1=lbound(orbs%inwhichlocreg,1)
        iie1=ubound(orbs%inwhichlocreg,1)
        allocate(forbs%inwhichlocreg(iis1:iie1), stat=istat)
        call memocc(istat, forbs%inwhichlocreg, 'forbs%inwhichlocreg', subname)
        do i1=iis1,iie1
            forbs%inwhichlocreg(i1) = orbs%inwhichlocreg(i1)
        end do
    end if

    if(associated(forbs%onwhichatom)) then
        iall=-product(shape(forbs%onwhichatom))*kind(forbs%onwhichatom)
        deallocate(forbs%onwhichatom, stat=istat)
        call memocc(istat, iall, 'forbs%onwhichatom', subname)
    end if
    if(associated(orbs%onwhichatom)) then
        iis1=lbound(orbs%onwhichatom,1)
        iie1=ubound(orbs%onwhichatom,1)
        allocate(forbs%onwhichatom(iis1:iie1), stat=istat)
        call memocc(istat, forbs%onwhichatom, 'forbs%onwhichatom', subname)
        do i1=iis1,iie1
            forbs%onwhichatom(i1) = orbs%onwhichatom(i1)
        end do
    end if

    if(associated(forbs%isorb_par)) then
        iall=-product(shape(forbs%isorb_par))*kind(forbs%isorb_par)
        deallocate(forbs%isorb_par, stat=istat)
        call memocc(istat, iall, 'forbs%isorb_par', subname)
    end if
    if(associated(orbs%isorb_par)) then
        iis1=lbound(orbs%isorb_par,1)
        iie1=ubound(orbs%isorb_par,1)
        allocate(forbs%isorb_par(iis1:iie1), stat=istat)
        call memocc(istat, forbs%isorb_par, 'forbs%isorb_par', subname)
        do i1=iis1,iie1
            forbs%isorb_par(i1) = orbs%isorb_par(i1)
        end do
    end if

    if(associated(forbs%ispot)) then
        iall=-product(shape(forbs%ispot))*kind(forbs%ispot)
        deallocate(forbs%ispot, stat=istat)
        call memocc(istat, iall, 'forbs%ispot', subname)
    end if
    if(associated(orbs%ispot)) then
        iis1=lbound(orbs%ispot,1)
        iie1=ubound(orbs%ispot,1)
        allocate(forbs%ispot(iis1:iie1), stat=istat)
        call memocc(istat, forbs%ispot, 'forbs%ispot', subname)
        do i1=iis1,iie1
            forbs%ispot(i1) = orbs%ispot(i1)
        end do
    end if

  end subroutine orbs_to_min_orbs_copy

  ! point minimal orbs structure to a given full orbs structure
  subroutine orbs_to_min_orbs_point(orbs,forbs)
    use module_types
    implicit none
    ! Calling arguments
    type(orbitals_data),intent(in):: orbs
    type(minimal_orbitals_data),intent(inout):: forbs

    forbs%norb = orbs%norb
    forbs%norbp = orbs%norbp
    forbs%isorb = orbs%isorb

    forbs%norb_par => orbs%norb_par
    forbs%inwhichlocreg => orbs%inwhichlocreg
    forbs%onwhichatom => orbs%onwhichatom
    forbs%isorb_par => orbs%isorb_par
    forbs%ispot => orbs%ispot
  end subroutine orbs_to_min_orbs_point



  subroutine calculate_fragment_density(frag,ndimrho,tmb,iorb_start,charge,atoms,rxyz,denspot)
    use module_types
    implicit none
    type(system_fragment), intent(inout) :: frag
    integer, intent(in) :: ndimrho ! add to fragment structure?
    integer, intent(in) :: iorb_start ! the first tmb for this fragment
    integer, intent(in) :: charge ! charge on this fragment for calculating correct kernel?!
    type(dft_wavefunction), intent(in) :: tmb
    type(atoms_data), intent(in) :: atoms ! just for plotting
    real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz ! just for plotting
    type(DFT_local_fields), intent(in) :: denspot ! just for plotting

    ! local variables
    integer :: iiorb, ilr, ind, indg, indr, ipt, jjorb, iorb, itmb
    real(kind=gp),dimension(:,:), allocatable :: kernel, fcoeff
    real(kind=wp),dimension(:), allocatable :: gpsi
    real(kind=wp),dimension(:), allocatable :: psir
    real(kind=gp) :: total_charge, tt, tt1, factor, hxh, hyh, hzh
    type(workarr_sumrho) :: w

    ! calculate fragment density and store in frag%fbasis%density
    ! allocate to cover whole simulation cell or just the fragment region?
    call f_routine(id='calculate_fragment_density')
    frag%fbasis%density=f_malloc_ptr(ndimrho,id='frag%fbasis%density')
    kernel=f_malloc((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='kernel')
    fcoeff=f_malloc((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='fcoeff')

    ! calculate density kernel for this fragment (ideally would use routine but don't have access to correct orbs and doing monoproc)
    ! add occup
    do iorb=1,frag%fbasis%forbs%norb
       do itmb=1,frag%fbasis%forbs%norb
          if (iorb<=floor((frag%nelec-charge)/2.0_gp)) then
             fcoeff(itmb,iorb) = 2.0_gp*frag%coeff(itmb,iorb)
          else if (iorb<=ceiling((frag%nelec-charge)/2.0_gp)) then
             fcoeff(itmb,iorb) = 1.0_gp*frag%coeff(itmb,iorb)
          else
             fcoeff(itmb,iorb) = 0.0_gp
          end if
       end do
    end do

    call dgemm('n', 't', frag%fbasis%forbs%norb, frag%fbasis%forbs%norb, &
         nint((frag%nelec-charge)/2.0_gp), 1.d0, &
         frag%coeff(1,1), frag%fbasis%forbs%norb, &
         fcoeff(1,1), frag%fbasis%forbs%norb, 0.d0, &
         kernel(1,1), frag%fbasis%forbs%norb)

    call f_free(fcoeff)

    ! expand all tmbs to global box to simplify density calculation and convert to real space
    ! as we don't have psi_full allocated, just using tmb%psi for now, but need to know where to start
    ind=1
    indg=1
    indr=1

    gpsi=f_malloc(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,id='gpsi')
    call to_zero(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,gpsi)
    call initialize_work_arrays_sumrho(tmb%lzd%glr, w)
    psir=f_malloc(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i*frag%fbasis%forbs%norb,id='psir')

    do iiorb=1,tmb%orbs%norb
       ilr = tmb%orbs%inwhichlocreg(iiorb)
       if (iiorb < iorb_start) then
          ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
          cycle
       else if (iiorb >= iorb_start + frag%fbasis%forbs%norb) then
          exit
       end if

       call Lpsi_to_global2(bigdft_mpi%iproc, tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f, &
            tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f, &
            1, 1, 1, tmb%Lzd%glr, tmb%Lzd%Llr(ilr), tmb%psi(ind), gpsi(indg))

       call daub_to_isf(tmb%lzd%glr, w, gpsi(indg), psir(indr))

       !indg = indg + tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f
       ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
       indr = indr + tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i
    end do

    call deallocate_work_arrays_sumrho(w)
    call f_free(gpsi)


    ! calculate density using kernel and psi
    hxh=.5d0*tmb%lzd%hgrids(1)
    hyh=.5d0*tmb%lzd%hgrids(2)
    hzh=.5d0*tmb%lzd%hgrids(3)
    factor=1.d0/(hxh*hyh*hzh)
    total_charge=0.d0
    !call to_zero(ndimrho,frag%fbasis%density)
    do ipt=1,ndimrho
       tt=1.e-20_dp
       do iiorb=1,frag%fbasis%forbs%norb
          tt1=psir(ipt+(iiorb-1)*ndimrho)
          tt=tt+kernel(iiorb,iiorb)*tt1**2
          do jjorb=iiorb+1,frag%fbasis%forbs%norb
             tt=tt+2.0_dp*kernel(jjorb,iiorb)*tt1*psir(ipt+(jjorb-1)*ndimrho)
          end do
       end do  
        tt=factor*tt
        total_charge=total_charge+tt
        ! apply a cut-off
        if (tt >= 1.e-6) then
           frag%fbasis%density(ipt)=tt
        else
           frag%fbasis%density(ipt)=0.0_gp
        end if
    end do

    print*,'Fragment density total charge',total_charge*hxh*hyh*hzh

    !call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'fragment_density.cube', &
    !     atoms,rxyz,denspot%dpbox,1,frag%fbasis%density)

    ! Define some constant factors.
    !hxh=.5d0*hx
    !hyh=.5d0*hy
    !hzh=.5d0*hz
    !factor=1.d0/(hxh*hyh*hzh)

    !do ipt=1,collcom_sr%nptsp_c
    !    ii=collcom_sr%norb_per_gridpoint_c(ipt)
    !    i0=collcom_sr%isptsp_c(ipt)
    !    tt=1.e-20_dp
    !    do i=1,ii
    !        iiorb=collcom_sr%indexrecvorbital_c(i0+i)
    !        tt1=collcom_sr%psit_c(i0+i)
    !        ind=denskern%matrixindex_in_compressed(iiorb,iiorb)
    !        tt=tt+denskern%matrix_compr(ind)*tt1*tt1
    !        do j=i+1,ii
    !            jjorb=collcom_sr%indexrecvorbital_c(i0+j)
    !            ind=denskern%matrixindex_in_compressed(jjorb,iiorb)
    !            if (ind==0) cycle
    !            tt=tt+2.0_dp*denskern%matrix_compr(ind)*tt1*collcom_sr%psit_c(i0+j)
    !        end do
    !    end do
    !    tt=factor*tt
    !    total_charge=total_charge+tt
    !    rho_local(ipt)=tt
    !end do

    call f_free(psir)
    call f_free(kernel)
    call f_release_routine()


  end subroutine calculate_fragment_density


  pure function minimal_orbitals_data_null() result(forbs)
    implicit none
    type(minimal_orbitals_data) :: forbs
    call nullify_minimal_orbitals_data(forbs)
  end function minimal_orbitals_data_null

  pure subroutine nullify_minimal_orbitals_data(forbs)
    implicit none
    type(minimal_orbitals_data), intent(out) :: forbs

    forbs%norb=0
    forbs%norbp=0
    forbs%isorb=0
    nullify(forbs%inwhichlocreg)
    nullify(forbs%onwhichatom)
    nullify(forbs%isorb_par)
    nullify(forbs%ispot)
    nullify(forbs%norb_par)
  end subroutine nullify_minimal_orbitals_data


  pure function fragment_null() result(frag)
    implicit none
    type(system_fragment) :: frag

    frag%nat_env=0
    nullify(frag%rxyz_env)
    frag%nelec=0
    nullify(frag%coeff)
    nullify(frag%kernel)
    nullify(frag%eval)
    call nullify_atomic_structure(frag%astruct_frg)
    ! nullify fragment basis
    call nullify_fragment_basis(frag%fbasis)

  end function fragment_null

  pure function fragment_basis_null() result(basis)
    implicit none
    type(fragment_basis) :: basis
    call nullify_fragment_basis(basis)
  end function fragment_basis_null

  pure subroutine nullify_fragment_basis(basis)
    implicit none
    type(fragment_basis), intent(out) :: basis

    basis%npsidim_orbs=0
    basis%npsidim_comp=0
    call nullify_local_zone_descriptors(basis%lzd)
    call nullify_minimal_orbitals_data(basis%forbs)
    !basis%forbs=minimal_orbitals_data_null()
    nullify(basis%psi_full)
    nullify(basis%phi)
    nullify(basis%density)
  end subroutine nullify_fragment_basis

  subroutine minimal_orbitals_data_free(forbs)
    implicit none
    type(minimal_orbitals_data), intent(inout) :: forbs

    if (associated(forbs%inwhichlocreg)) call f_free_ptr(forbs%inwhichlocreg)
    if (associated(forbs%onwhichatom)) call f_free_ptr(forbs%onwhichatom)
    if (associated(forbs%isorb_par)) call f_free_ptr(forbs%isorb_par)
    if (associated(forbs%ispot)) call f_free_ptr(forbs%ispot)
    if (associated(forbs%norb_par)) call f_free_ptr(forbs%norb_par)
    forbs=minimal_orbitals_data_null()
  end subroutine minimal_orbitals_data_free

  subroutine fragment_basis_free(basis)
    implicit none
    type(fragment_basis), intent(inout) :: basis
    character(len=200) :: subname

    subname='fragment_basis_free'
    call deallocate_local_zone_descriptors(basis%lzd,subname)
    call minimal_orbitals_data_free(basis%forbs)
    if (associated(basis%psi_full)) call f_free_ptr(basis%psi_full)
    if (associated(basis%density)) call f_free_ptr(basis%density)
    basis=fragment_basis_null()
  end subroutine fragment_basis_free

  subroutine fragment_free(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag

    call deallocate_atomic_structure(frag%astruct_frg)
    frag%astruct_frg=atomic_structure_null()
    call f_free_ptr(frag%rxyz_env)
    call f_free_ptr(frag%coeff)
    call f_free_ptr(frag%kernel)
    call f_free_ptr(frag%eval)
    call fragment_basis_free(frag%fbasis)
    frag=fragment_null()

  end subroutine fragment_free

  subroutine fragment_allocate(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag

    call f_routine(id='fragment_allocate')

    frag%rxyz_env=f_malloc_ptr((/3,frag%nat_env/),id='frag%rxyz_env')
    frag%coeff=f_malloc_ptr((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='frag%coeff')
    frag%kernel=f_malloc_ptr((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='frag%kernel')
    frag%eval=f_malloc_ptr(frag%fbasis%forbs%norb,id='frag%eval')

    call f_release_routine()

  end subroutine fragment_allocate


  pure function frg_center(frag)
    implicit none
    type(system_fragment), intent(in) :: frag
    real(gp), dimension(3) :: frg_center
    !local variables
    integer :: iat

    frg_center=0.0_gp
    do iat=1,frag%astruct_frg%nat
       frg_center=frg_center+frag%astruct_frg%rxyz(:,iat)
    end do
    frg_center=frg_center/real(frag%astruct_frg%nat,gp)

  end function frg_center

  !function transform_fragment(trans,frag) result(frag_new)
  !  implicit none
  !  type(fragment_transformation), intent(in) :: trans
  !  type(system_fragment), intent(in) :: frag
  !  type(system_fragment) :: frag_new
  !
  !  ! local variables
  !  integer :: iat
  !
  !  frag_new=fragment_null()
  !  frag_new%astruct_frg%nat=frag%astruct_frg%nat
  !  frag_new%nat_env=frag%nat_env
  !  frag_new%astruct_frg%ntypes=frag%astruct_frg%ntypes
  !
  !  ! allocate arrays here, leave fragment_basis nullified
  !  call fragment_allocate(frag_new)
  !
  !  ! do fragment first, then environment
  !  do iat=1,frag%astruct_frg%nat
  !     frag_new%astruct_frg%rxyz(:,iat)=rotate_vector(trans%rot_axis,&
  !          trans%theta,frag%astruct_frg%rxyz(:,iat)-trans%rot_center(:))
  !     frag_new%astruct_frg%rxyz(:,iat)=frag_new%astruct_frg%rxyz(:,iat)+trans%rot_center(:)+trans%dr(:)
  !  end do
  !
  !  do iat=1,frag%nat_env
  !     frag_new%rxyz_env(:,iat)=rotate_vector(trans%rot_axis,trans%theta,frag%rxyz_env(:,iat)-trans%rot_center(:))
  !     frag_new%rxyz_env(:,iat)=frag_new%rxyz_env(:,iat)+trans%rot_center(:)+trans%dr(:)
  !  end do
  !
  !  ! to complete should copy across iatype and atomnames but only using as a check to see how effective trans is
  !
  !end function transform_fragment



  !> Express the coordinates of a vector into a rotated reference frame
  pure function rotate_vector(newz,theta,vec) result(vecn)
     use module_base
     implicit none
     real(gp), intent(in) :: theta
     real(gp), dimension(3), intent(in) :: newz,vec
     real(gp), dimension(3) :: vecn
     !local variables
     real(gp) :: sint,cost,onemc,x,y,z

     !save recalculation
     sint=sin(theta)
     cost=cos(theta)
     onemc=1.0_gp-cost
     x=vec(1)
     y=vec(2)
     z=vec(3)

     vecn(1)=x*(cost + onemc*newz(1)**2) + y*(onemc*newz(1)*newz(2) - sint*newz(3)) &
          + z*(sint*newz(2) + onemc*newz(1)*newz(3))
     vecn(2)=y*(cost + onemc*newz(2)**2) + x*(onemc*newz(1)*newz(2) + sint*newz(3)) &
          + z*(-(sint*newz(1)) + onemc*newz(2)*newz(3))
     vecn(3)=z*(cost + onemc*newz(3)**2) + x*(onemc*newz(1)*newz(3) - sint*newz(2)) &
          + y*(sint*newz(1) + onemc*newz(2)*newz(3))

  end function rotate_vector

  !pure function transform_fragment_basis(trans,basis) result(basis_new)
  !  implicit none
  !  type(fragment_transformation), intent(in) :: trans
  !  type(fragment_basis), intent(in) :: basis
  !  type(fragment_basis) :: basis_new

  !  basis_new=fragment_basis_null()

  !  ! minimal_orbitals_data should remain the same
  !! want to have defined new lzd already?  in which case not a pure function...

  !!   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
  !!   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
  !!   type(local_zone_descriptors) :: Lzd
  !!   type(minimal_orbitals_data) :: forbs
  !!   real(wp), dimension(:), pointer :: psi

  !end function transform_fragment_basis

  subroutine find_frag_trans(nat,rxyz_ref,rxyz_new,frag_trans)
    use module_base
    use yaml_output
    implicit none
    integer, intent(in) :: nat !< fragment size
    real(gp), dimension(3,nat), intent(in) :: rxyz_ref,rxyz_new !<coordinates measured wrt rot_center
    type(fragment_transformation), intent(inout) :: frag_trans
    !local variables
    integer, parameter :: lwork=7*3
    integer :: info,iat,i_stat!,i
    real(gp) :: dets,J
    real(gp), dimension(3) :: SM_arr !< array of SVD and M array
    real(gp), dimension(lwork) :: work !< array of SVD and M array
    real(gp), dimension(3,nat) :: J_arr !< matrix for calculating Wahba's cost function
    real(gp), dimension(3,3) :: B_mat,R_mat,U_mat,VT_mat !<matrices of Wahba's problem
    character(len=100) :: subname

    subname='find_frag_trans'

    B_mat=0.0_gp
    R_mat=0.0_gp

    !all positions are of weight one for the moment
    call dgemm('N','T',3,3,nat,1.0_gp,rxyz_new,3,rxyz_ref,3,0.0_gp,B_mat,3)

    !find matrix of svd
    call dgesvd('A','A',3,3,B_mat,3,SM_arr,U_mat,3,VT_mat,3,work,lwork,info)
    if (f_err_raise(info/=0,'Problem in DGESVD')) return
    
    !multiply last line of VT_mat by det(U)*det(V)
    dets=det_33(U_mat)*det_33(VT_mat)
    VT_mat(3,:)=VT_mat(3,:)*dets

    !find rotation matrix
    call dgemm('N','N',3,3,3,1.0_gp,U_mat,3,VT_mat,3,0.0_gp,R_mat,3)

    !find the angle from R matrix
    frag_trans%theta=theta_from_r(R_mat)

    !find rot_axis
    frag_trans%rot_axis=axis_from_r(R_mat)

    !print*,'Rmat:',frag_trans%theta
    !do i=1,3
    !   write(*,'(3(F12.6,2x))') R_mat(i,:)
    !end do

    !print*,'Rcalc:',frag_trans%theta
    !write(*,'(3(F12.6,2x))') cos(frag_trans%theta)+frag_trans%rot_axis(1)**2*(1.0_gp-cos(frag_trans%theta)),&
    !     frag_trans%rot_axis(1)*frag_trans%rot_axis(2)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(3)*sin(frag_trans%theta),&
    !     frag_trans%rot_axis(1)*frag_trans%rot_axis(3)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(2)*sin(frag_trans%theta)
    !write(*,'(3(F12.6,2x))') frag_trans%rot_axis(2)*frag_trans%rot_axis(1)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(3)*sin(frag_trans%theta),&
    !     cos(frag_trans%theta)+frag_trans%rot_axis(2)**2*(1.0_gp-cos(frag_trans%theta)),&
    !     frag_trans%rot_axis(2)*frag_trans%rot_axis(3)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(1)*sin(frag_trans%theta)
    !write(*,'(3(F12.6,2x))') frag_trans%rot_axis(3)*frag_trans%rot_axis(1)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(2)*sin(frag_trans%theta),&
    !     frag_trans%rot_axis(3)*frag_trans%rot_axis(2)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(1)*sin(frag_trans%theta),&
    !     cos(frag_trans%theta)+frag_trans%rot_axis(3)**2*(1.0_gp-cos(frag_trans%theta))
    
    !to be verified that the cost function of Wahba's problem is little
    J_arr=rxyz_new
    call dgemm('N','N',3,nat,3,-1.0_gp,R_mat,3,rxyz_ref,3,1.0_gp,J_arr,3)
    J=0.0_gp
    do iat=1,nat
       J=J+J_arr(1,iat)**2+J_arr(2,iat)**2+J_arr(3,iat)**2
    end do

    if (J>1.0e-3) then
       write(*,'(a,2es18.8)') "Error, Wahba's cost function is too big",J,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
    end if

    !check the pertinence of the suggested rotation
    !if (abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp)) print*,'frag_trans%theta=',frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
    !if  (f_err_raise(abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp),'Angle frag_trans%theta not optimal (frag_trans%theta= '//&
    !      yaml_toa(frag_trans%theta)//' )')) return

  end subroutine find_frag_trans

  pure function theta_from_r(R_mat) result(theta)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R_mat

    real(gp) :: theta

    !local variables
    real(gp) :: tracem1

    tracem1=R_mat(1,1)+R_mat(2,2)+R_mat(3,3)-1.0_gp

    if (abs(tracem1) - 2.0_gp > -1.e-14_gp) then
      if (tracem1 > 0.0_gp) theta = 0.0_gp
      if (tracem1 <= 0.0_gp) theta = 180.0_gp*(4.0_gp*atan(1.d0)/180.0_gp)
    else
       theta=acos(0.5_gp*tracem1)
    end if

  end function theta_from_r

  function axis_from_r(R_mat) result(rot_axis)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R_mat
    real(gp), dimension(3) :: rot_axis

    !local variables
    real(gp) :: dnrm2, norm
    integer :: i

    rot_axis(1)=R_mat(3,2)-R_mat(2,3)    
    rot_axis(2)=R_mat(1,3)-R_mat(3,1)
    rot_axis(3)=R_mat(2,1)-R_mat(1,2)    

    !normalize it

    norm=dnrm2(3,rot_axis,1)
    !print*,'axis_from_r',norm,rot_axis
    if (norm>=1.e-5_gp) then
       call dscal(3,1.0_gp/norm,rot_axis,1)
       !print*,'axis_from_r2',norm,rot_axis
    else
       ! squares of rot_axis are diag 0.5(R+I), signs as before
       !this is not good as it gives a array of modulus bigger than one
       !rot_axis(1:2)=0.0_gp
       !rot_axis(3)=1.0_gp
       do i=1,3
          if (R_mat(i,i)<-1.0_gp.and.R_mat(i,i)>-1.0_gp-1.0e-5_gp) then
                rot_axis(i)=0.0_gp
          else if (R_mat(i,i)<=-1.0_gp-1.0e-5_gp) then
             stop 'Problem in assigning axis from Rmat'
          else
             rot_axis(i)=sign(dsqrt(0.5_gp*(R_mat(i,i)+1.0_gp)),rot_axis(i))
          end if
       end do

       !print*,'axis_from_r3a',rot_axis,R_mat(1,1),R_mat(2,2),R_mat(3,3)
       norm=dnrm2(3,rot_axis,1)
       call dscal(3,1.0_gp/norm,rot_axis,1)
       !print*,'axis_from_r3',norm,rot_axis
    end if

  end function axis_from_r



  pure function frag_center(nat,rxyz) result(cen)
    implicit none
    integer, intent(in) :: nat 
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(gp), dimension(3) :: cen
    !local variables
    integer :: iat,i
    
    cen=0.0_gp
    if (nat > 0) then
       do iat=1,nat
          do i=1,3
             cen(i)=cen(i)+rxyz(i,iat)
          end do
       end do
       cen=cen/real(nat,gp)
    end if
    
  end function frag_center

  !>determinant of a 3x3 matrix
  pure function det_33(a) result(det)
    implicit none
    real(gp), dimension(3,3), intent(in) :: a
    real(gp) :: det

    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
         + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3))  &
         + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function det_33

  subroutine fragment_coeffs_to_kernel(iproc,input,input_frag_charge,ref_frags,tmb,ksorbs,overlap_calculated,&
    nstates_max,cdft)
    use yaml_output
    use module_base
    use communications, only: transpose_localized
    use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
    use sparsematrix, only: uncompress_matrix
    implicit none
    type(DFT_wavefunction), intent(inout) :: tmb
    type(input_variables), intent(in) :: input
    type(system_fragment), dimension(input%frag%nfrag_ref), intent(inout) :: ref_frags
    type(orbitals_data), intent(inout) :: ksorbs
    logical, intent(inout) :: overlap_calculated
    real(kind=gp), dimension(input%frag%nfrag), intent(in) :: input_frag_charge
    integer, intent(in) :: iproc
    integer, intent(out) :: nstates_max ! number of states in total if we consider all partially occupied fragment states to be fully occupied
    logical, intent(in) :: cdft

    integer :: iorb, isforb, jsforb, ifrag, ifrag_ref, itmb, jtmb, iall, istat, num_extra_per_frag, ierr
    integer, allocatable, dimension(:) :: ipiv
    real(gp), dimension(:,:), allocatable :: coeff_final, ks, ksk
    !*real(gp), dimension(:), allocatable :: kernel_final
    real(gp) :: nonidem, nelecorbs, nelecfrag_tot, jstate_max, homo_diff, lag_mult
    real(gp), dimension(:), allocatable :: eval_tmp, eval_tmp2
    character(len=*), parameter :: subname='fragment_coeffs_to_kernel'

    integer :: rand_size
    integer, allocatable, dimension(:) :: rand_seed
    real(kind=dp) :: rtime, random_noise, rmax
    character(len=10) :: sys_time
    logical :: random, completely_random

    real(wp), dimension(:,:,:), pointer :: mom_vec_fake

    call timing(iproc,'kernel_init','ON')
    call f_routine(id='fragment_coeffs_to_kernel')

    ! need to do this properly/rearrange routines
    if (cdft) then
       ! otherwise doesn't make sense
       if (input%frag%nfrag==2) homo_diff=(ref_frags(1)%eval(ceiling(ref_frags(1)%nelec/2.0_gp))&
            -ref_frags(2)%eval(ceiling(ref_frags(2)%nelec/2.0_gp)))/2.0d0
       !if (cdft%charge<0) lag_mult=-0.5, otherwise +0.5
       lag_mult=-0.05d0
    else
       homo_diff=0.0d0
       lag_mult=0.0d0
    end if

    ! adding random noise to starting to help with local minima problem
    random=.false. ! add a bit of noise
    completely_random=.false. ! completely random start for coeffs

    rmax=0.2d0
    random_noise=0.0d0
    rtime=0.0d0
    if (random .or. completely_random) then
       call random_seed(size=rand_size)
       allocate(rand_seed(1:rand_size))
       call date_and_time(time=sys_time)
       ! coeffs need to be the same across processors
       if (iproc==0) read(sys_time,*) rtime
       if (bigdft_mpi%nproc>1) call mpiallred(rtime, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
       rand_seed=int(rtime*1000.0_dp)
       call random_seed(put=rand_seed)
       deallocate(rand_seed) 
    end if
    nstates_max=0
    nelecfrag_tot=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)
       nelecfrag_tot=nelecfrag_tot+ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag)
    end do

    if (completely_random) then
       if (bigdft_mpi%iproc==0) print*,'Starting coeffs are replaced with a random guess'
    else if (random) then
       if (bigdft_mpi%iproc==0) print*,'Random noise added to starting coeffs'
    end if

    ! in theory we could add/remove states depending on their energies, but for now we force the user to specify
    ! need to include occupations as we actually want to compare number of electrons here?
    nelecorbs=0
    do iorb=1,ksorbs%norb
       nelecorbs=nelecorbs+ksorbs%occup(iorb)
    end do

    if (nint(nelecorbs)/=nelecfrag_tot) then
       print*,'User must specify which fragments charges are added to/removed from in charged fragment calculation',&
            nelecfrag_tot,nelecorbs,ksorbs%norb
       stop
    end if

    if (mod(input%norbsempty,input%frag%nfrag)/=0) then
       if (bigdft_mpi%iproc==0) print*,'Warning, number of extra bands does not divide evenly among fragments'
       num_extra_per_frag=(input%norbsempty-mod(input%norbsempty,input%frag%nfrag))/input%frag%nfrag
    else
       num_extra_per_frag=input%norbsempty/input%frag%nfrag
    end if


    eval_tmp=f_malloc(tmb%orbs%norb,id='eval_tmp')
    coeff_final=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_final')
    !*kernel_final=f_malloc(tmb%linmat%denskern%nvctr,id='kernel_final')
    !ref_frags(ifrag_ref)%kernel=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ref_frags(ifrag_ref)%kernel')

    ! Calculate the overlap matrix between the TMBs.
    if(.not. overlap_calculated) then
       if(.not.tmb%can_use_transposed) then
           if(associated(tmb%psit_c)) then
               iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
               deallocate(tmb%psit_c, stat=istat)
               call memocc(istat, iall, 'tmb%psit_c', subname)
           end if
           if(associated(tmb%psit_f)) then
               iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
               deallocate(tmb%psit_f, stat=istat)
               call memocc(istat, iall, 'tmb%psit_f', subname)
           end if
           allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
           call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
           allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
           call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
           call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
           tmb%can_use_transposed=.true.
       end if
       !call timing(iproc,'renormCoefComp','OF')

       call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, &
            tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
       ! This can then be deleted if the transition to the new type has been completed.
       !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr


       !call timing(iproc,'renormCoefComp','ON')
       overlap_calculated=.true.
    end if

    ! copy from coeff fragment to global coeffs - occupied states only
    isforb=0
    jsforb=0
    call to_zero(tmb%orbs%norb*tmb%orbs%norb,coeff_final(1,1))
    !*call to_zero(tmb%linmat%denskern%nvctr,kernel_final(1))
    !!tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%ovrlp%matrix')
    !!call uncompress_matrix(iproc,tmb%linmat%ovrlp)
    do ifrag=1,input%frag%nfrag
       ! find reference fragment this corresponds to
       ifrag_ref=input%frag%frag_index(ifrag)
       call to_zero(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1))

       jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
       !jstate_max=ref_frags(ifrag_ref)%nelec/2.0_gp+num_extra_per_frag
       do jtmb=1,ceiling(jstate_max)

          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=1,isforb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jtmb)=random_noise
             end do
          end if

          do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
             if (random.or.completely_random) then
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
             end if
              if (.not. completely_random) then
                tmb%coeff(isforb+itmb,jtmb)=ref_frags(ifrag_ref)%coeff(itmb,jtmb)+random_noise
             else
                tmb%coeff(isforb+itmb,jtmb)=random_noise
             end if
             tmb%orbs%eval(jsforb+jtmb)=ref_frags(ifrag_ref)%eval(jtmb)-((-1)**(ifrag))*lag_mult-homo_diff
          end do

          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb+1,tmb%orbs%norb
                 call random_number(random_noise)
                 random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jtmb)=random_noise
             end do
          end if

          if (ceiling(jstate_max)/=jstate_max.and.jtmb==ceiling(jstate_max)) then
             tmb%orbs%occup(jtmb+jsforb)=(jstate_max*2.0d0)-2*(ceiling(jstate_max)-1)
             ! want partly occupied states to be last of unoccupied
             eval_tmp(jsforb+jtmb)=tmb%orbs%eval(jsforb+jtmb)+10.0d0
          else
             tmb%orbs%occup(jtmb+jsforb)=2.0d0
             eval_tmp(jsforb+jtmb)=tmb%orbs%eval(jsforb+jtmb)
          end if
          !if (bigdft_mpi%iproc==0) print*,'ifrag,jtmb,occ,iorb',ifrag,jtmb,tmb%orbs%occup(jtmb+jsforb),jtmb+jsforb
       end do
       nstates_max=nstates_max+ceiling(jstate_max)

       ! debug
       !do itmb=1,tmb%orbs%norb
       !   do jtmb=1,tmb%orbs%norb
       !      write(40+ifrag,*) itmb,jtmb,tmb%coeff(itmb,jtmb)
       !  end do
       !end do
       ! end debug

       !call to_zero(tmb%linmat%denskern%nvctr,tmb%linmat%denskern%matrix_compr(1))

       ! should correct the occupation for kernel here, but as we replace the smaller kernel with the correct bigger kernel
       ! don't worry about this for now

       ! reorthonormalize the coeffs for each fragment - don't need unoccupied states here
       tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, &
                                  iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
       call reorthonormalize_coeff(bigdft_mpi%iproc, bigdft_mpi%nproc, &
            ceiling((ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp), &
            tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, tmb%orthpar%methTransformOverlap,&
            tmb%orbs, tmb%linmat%ovrlp, tmb%linmat%ovrlp_, tmb%coeff, ksorbs)
       call f_free_ptr(tmb%linmat%ovrlp_%matrix)

       !! debug
       !!output final kernel
       !! 20 - if just calculate, 21 if reconstruct total, 22 if reconstruct then sum
       !tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
       !call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%denskern)
       !!do itmb=1,tmb%orbs%norb
       !!   do jtmb=1,tmb%orbs%norb
       !!      write(30+ifrag,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),tmb%coeff(itmb,jtmb)
       !!   end do
       !!end do
       !do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
       !   do jtmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
       !      ref_frags(ifrag_ref)%kernel(itmb,jtmb)=tmb%linmat%denskern%matrix(itmb,jtmb).....
       !   end do
       !end do
       !call f_free_ptr(tmb%linmat%denskern%matrix)    
       !! end debug

       ! assemble complete kernel from separate fragment kernels
       !call daxpy(tmb%linmat%denskern%nvctr,1.0d0,tmb%linmat%denskern%matrix_compr(1),1,kernel_final(1),1)

       ! update coeff_final matrix following coeff reorthonormalization
       do jtmb=1,ceiling(jstate_max)
          do itmb=1,tmb%orbs%norb
             coeff_final(itmb,jsforb+jtmb)=tmb%coeff(itmb,jtmb)
          end do
       end do

       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
       jsforb=jsforb+ceiling(jstate_max)
    end do
    !!call f_free_ptr(tmb%linmat%ovrlp%matrix)

    !*call vcopy(tmb%linmat%denskern%nvctr,kernel_final(1),1,tmb%linmat%denskern%matrix_compr(1),1)
    call vcopy(tmb%orbs%norb*tmb%orbs%norb,coeff_final(1,1),1,tmb%coeff(1,1),1)

    !*call f_free(kernel_final)
    call f_free(coeff_final)

    ! debug
    !output final kernel
    ! 20 - if just calculate, 21 if reconstruct total, 22 if reconstruct then sum
    !tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
    !call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%denskern)
    !do itmb=1,tmb%orbs%norb
    !   do jtmb=1,tmb%orbs%norb
    !      write(22,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),tmb%coeff(itmb,jtmb)
    !   end do
    !end do

    ! check final kernel is idempotent
    !tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%ovrlp%matrix')
    !ks=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ks')
    !ksk=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ksk')
    !call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%ovrlp)
    !call dgemm('n', 't', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, &
    !           tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, 0.d0, ks(1,1), tmb%orbs%norb) 
    !call dgemm('n', 't', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, ks(1,1), tmb%orbs%norb, &
    !           tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 0.d0, ksk(1,1), tmb%orbs%norb)

    !nonidem=0
    !do itmb=1,tmb%orbs%norb
    !   do jtmb=1,tmb%orbs%norb
    !      write(60,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),ksk(itmb,jtmb),&
    !           tmb%linmat%denskern%matrix(itmb,jtmb)-ksk(itmb,jtmb)
    !      nonidem=nonidem+tmb%linmat%denskern%matrix(itmb,jtmb)-ksk(itmb,jtmb)
    !   end do
    !end do
    !print*,'non idempotency',nonidem/tmb%orbs%norb**2

    !call f_free(ks) 
    !call f_free(ksk) 
    !call f_free_ptr(tmb%linmat%ovrlp%matrix)   
    !call f_free_ptr(tmb%linmat%denskern%matrix)    
    ! end debug

    ! add unoccupied states to complete coeffs
    isforb=0
    do ifrag=1,input%frag%nfrag
       ! find reference fragment this corresponds to
       ifrag_ref=input%frag%frag_index(ifrag)
       jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
       !jstate_max=ref_frags(ifrag_ref)%nelec/2.0_gp+num_extra_per_frag
       do jtmb=ceiling(jstate_max)+1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=1,isforb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max))=random_noise
             end do
          end if

          do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
             if (random.or.completely_random) then
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
             end if
             if (.not. completely_random) then
                tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max))=ref_frags(ifrag_ref)%coeff(itmb,jtmb)+random_noise
             else
                tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max))=random_noise
             end if
             tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max))=ref_frags(ifrag_ref)%eval(jtmb)-((-1)**(ifrag))*lag_mult-homo_diff
             eval_tmp(jsforb+jtmb-ceiling(jstate_max))=tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max))+20.0d0
          end do

          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb+1,tmb%orbs%norb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max))=random_noise
             end do
          end if

          tmb%orbs%occup(jsforb+jtmb-ceiling(jstate_max))=0.0d0
          !if (bigdft_mpi%iproc==0) print*,'ifrag,jtmb,occ,iorb',ifrag,jtmb,0.0,jsforb+jtmb-ceiling(jstate_max)
       end do

       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
       jsforb=jsforb+ref_frags(ifrag_ref)%fbasis%forbs%norb-ceiling(jstate_max)
    end do
    if (.not. completely_random) then
       !!print*,'nstates_max:',nstates_max,ksorbs%norb,tmb%orbs%norb
       if (bigdft_mpi%iproc==0) then
           call yaml_map('nstates_max',nstates_max)
           call yaml_map('ksorbs%norb',ksorbs%norb)
           call yaml_map('tmb%orbs%norb',tmb%orbs%norb)
       end if

       ! reorder unoccupied states so that extra states functions correctly
       ! still needed just in case number of empty bands doesn't divide by number of fragments
       ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
       call order_coeffs_by_energy(tmb%orbs%norb-nstates_max,tmb%orbs%norb,tmb%coeff(1,nstates_max+1),&
            eval_tmp(nstates_max+1),ipiv(1))!,tmb%orbs%eval(nstates_max+1))
       eval_tmp2=f_malloc(tmb%orbs%norb-nstates_max,id='eval_tmp2')
       call vcopy(tmb%orbs%norb-nstates_max,tmb%orbs%eval(nstates_max+1),1,eval_tmp2(1),1)
       do itmb=nstates_max+1,tmb%orbs%norb
          tmb%orbs%eval(itmb)=eval_tmp2(ipiv(itmb-nstates_max))
       end do
       call vcopy(tmb%orbs%norb-nstates_max,tmb%orbs%occup(nstates_max+1),1,eval_tmp2(1),1)
       do itmb=nstates_max+1,tmb%orbs%norb
          tmb%orbs%occup(itmb)=eval_tmp2(ipiv(itmb-nstates_max))
       end do
       call f_free(eval_tmp2)
       ! reorder ksorbs%norb states by energy - no longer taking charge as input
       call order_coeffs_by_energy(ksorbs%norb,tmb%orbs%norb,tmb%coeff(1,1),eval_tmp(1),ipiv(1))!,tmb%orbs%eval(1))
               !eval_tmp2=f_malloc(tmb%orbs%norb,id='eval_tmp2')
               !call vcopy(tmb%orbs%norb,tmb%orbs%occup(1),1,eval_tmp2(1),1)
               !do itmb=1,ksorbs%norb
               !   tmb%orbs%occup(itmb)=eval_tmp2(ipiv(itmb))
               !end do
               !call vcopy(tmb%orbs%norb,tmb%orbs%eval(1),1,eval_tmp2(1),1)
               !call vcopy(tmb%orbs%norb,eval_tmp(1),1,tmb%orbs%eval(1),1)
               !nullify(mom_vec_fake)
               !if (bigdft_mpi%iproc==0) then 
               !   call write_eigenvalues_data(0.1d0,tmb%orbs,mom_vec_fake)
               !end if
               !call vcopy(tmb%orbs%norb,eval_tmp2(1),1,tmb%orbs%eval(1),1)
               !call f_free(eval_tmp2)
       call vcopy(ksorbs%norb,tmb%orbs%eval(1),1,eval_tmp(1),1)
       do itmb=1,ksorbs%norb
          tmb%orbs%eval(itmb)=eval_tmp(ipiv(itmb))
       end do
       call f_free(eval_tmp)
       call f_free(ipiv)
       ! debug
       !do itmb=1,tmb%orbs%norb
       !   do jtmb=1,tmb%orbs%norb
       !      write(50,*) itmb,jtmb,tmb%coeff(itmb,jtmb)
       !   end do
       !end do
       ! end debug
       ! print starting eigenvalues
       if(bigdft_mpi%iproc==0) then
          !!write(*,'(1x,a)') '-------------------------------------------------'
          !!write(*,'(1x,a)') 'some selected eigenvalues:'
          !!do iorb=1,tmb%orbs%norb!max(ksorbs%norb-8,1),min(ksorbs%norb+8,tmb%orbs%norb)
          !!    if(iorb==ksorbs%norb) then
          !!        write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- last occupied orbital'
          !!    else if(iorb==ksorbs%norb+1) then
          !!        write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- first virtual orbital'
          !!    else
          !!        write(*,'(3x,a,i0,a,es20.12)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb)
          !!    end if
          !!end do
          !!write(*,'(1x,a)') '-------------------------------------------------'
          !!write(*,'(1x,a,2es24.16)') 'lowest, highest ev:',tmb%orbs%eval(1),tmb%orbs%eval(tmb%orbs%norb)

          call yaml_open_sequence('TMB eigenvalues',flow=.true.)
          call yaml_newline()
          do iorb=1,tmb%orbs%norb
              call yaml_open_map(flow=.true.)
              call yaml_map('index',iorb)
              call yaml_map('value',tmb%orbs%eval(iorb),fmt='(es20.12)')
              call yaml_close_map()
              if(iorb==ksorbs%norb) then
                  !!write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- last occupied orbital'
                  call yaml_comment('  <-- last occupied orbital')
              else if(iorb==ksorbs%norb+1) then
                  !!write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- first virtual orbital'
                  call yaml_comment('  <-- first virtual orbital')
              else
                  !!write(*,'(3x,a,i0,a,es20.12)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb)
              end if
              call yaml_newline()
          end do
          call yaml_close_sequence()
          !!write(*,'(1x,a)') '-------------------------------------------------'
          !!write(*,'(1x,a,2es24.16)') 'lowest, highest ev:',tmb%orbs%eval(1),tmb%orbs%eval(tmb%orbs%norb)

       end if

       if (nstates_max/=ksorbs%norb) then
          if (bigdft_mpi%iproc==0) print*,'Warning, number of states with non-zero occupation in fragments (',nstates_max,&
               ') differs from number of KS states (',ksorbs%norb,') - might have convergence problems'
       end if

       !!!!!!!!!!!!!!!
       ! need the eigenvalues to be in ksorbs%eval
       call vcopy(ksorbs%norb,tmb%orbs%eval(1),1,ksorbs%eval(1),1)
       call evaltoocc(bigdft_mpi%iproc,bigdft_mpi%nproc,.false.,input%tel,ksorbs,input%occopt)

       nullify(mom_vec_fake)
       if (bigdft_mpi%iproc ==0) then 
          call write_eigenvalues_data(0.1d0,ksorbs,mom_vec_fake)
       end if
       !!!!!!!!!!!!!!!
    end if ! completely random
    call f_release_routine()
    call timing(iproc,'kernel_init','OF')

  end subroutine fragment_coeffs_to_kernel


end module module_fragments












