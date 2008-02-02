module module_types

  implicit none

  !-Input variable structure
  !structure of the variables read by input.dat file
  type, public :: input_variables
     integer :: ncount_cluster_x
     real(kind=8) :: frac_fluct,randdis,betax
     integer :: ixc,ncharge,itermax,ncong,idsx,ncongt,inputPsiId,nspin,mpol
     real(kind=8) :: hgrid,crmult,frmult,cpmult,fpmult,elecfield,gnrm_cv,rbuf
     logical :: output_grid,output_wf,calc_tail
  end type input_variables

  type, public :: kinetic_bounds
     integer, dimension(:,:,:), pointer :: ibyz_c,ibxz_c,ibxy_c
     integer, dimension(:,:,:), pointer :: ibyz_f,ibxz_f,ibxy_f
  end type kinetic_bounds

  type, public :: shrink_bounds
     integer, dimension(:,:,:), pointer :: ibzzx_c,ibyyzz_c
     integer, dimension(:,:,:), pointer :: ibxy_ff,ibzzx_f,ibyyzz_f
  end type shrink_bounds

  type, public :: grow_bounds
     integer, dimension(:,:,:), pointer :: ibzxx_c,ibxxyy_c
     integer, dimension(:,:,:), pointer :: ibyz_ff,ibzxx_f,ibxxyy_f
  end type grow_bounds

  type, public :: convolutions_bounds
     type(kinetic_bounds) :: kb
     type(shrink_bounds) :: sb
     type(grow_bounds) :: gb
     integer, dimension(:,:,:), pointer ::ibyyzz_r ! real space border
  end type convolutions_bounds

  type, public :: wavefunctions_descriptors
     integer :: nvctr_c,nvctr_f,nseg_c,nseg_f
     integer, dimension(:,:), pointer :: keyg
     integer, dimension(:), pointer :: keyv
  end type wavefunctions_descriptors

  type, public :: nonlocal_psp_descriptors
     !number of projectors and number of elements
     integer :: nproj,nprojel
     ! projector segments on real space grid
     integer, dimension(:), pointer :: nvctr_p,nseg_p,keyv_p
     integer, dimension(:,:), pointer :: keyg_p 
     ! Parameters for the boxes containing the projectors
     integer, dimension(:,:,:), pointer :: nboxp_c,nboxp_f
  end type nonlocal_psp_descriptors

  type, public :: atoms_data
     integer :: nat,ntypes,natsc
     character(len=20), dimension(:), pointer :: atomnames
     logical, dimension(:), pointer :: lfrztyp
     integer, dimension(:), pointer :: iatype,iasctype,nspinat,nelpsp,npspcode,nzatom
     real(kind=8), dimension(:,:,:), pointer :: psppar
  end type atoms_data

contains


  subroutine allocate_wfd(wfd,routine)
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    allocate(wfd%keyg(2,wfd%nseg_c+wfd%nseg_f),stat=i_stat)
    call memocc(i_stat,product(shape(wfd%keyg))*kind(wfd%keyg),'keyg',routine)
    allocate(wfd%keyv(wfd%nseg_c+wfd%nseg_f),stat=i_stat)
    call memocc(i_stat,product(shape(wfd%keyv))*kind(wfd%keyv),'keyv',routine)

  end subroutine allocate_wfd

  subroutine deallocate_wfd(wfd,routine)
    implicit none
    type(wavefunctions_descriptors) :: wfd
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(wfd%keyg))*kind(wfd%keyg)
    deallocate(wfd%keyg,stat=i_stat)
    call memocc(i_stat,i_all,'keyg',routine)
    i_all=-product(shape(wfd%keyv))*kind(wfd%keyv)
    deallocate(wfd%keyv,stat=i_stat)
    call memocc(i_stat,i_all,'keyv',routine)

  end subroutine deallocate_wfd

  subroutine deallocate_bounds(bounds,routine)
    implicit none
    type(convolutions_bounds) :: bounds
    character(len=*), intent(in) :: routine
    !local variables
    integer :: i_all,i_stat

    i_all=-product(shape(bounds%kb%ibyz_c))*kind(bounds%kb%ibyz_c)
    deallocate(bounds%kb%ibyz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_c',routine)
    i_all=-product(shape(bounds%kb%ibxz_c))*kind(bounds%kb%ibxz_c)
    deallocate(bounds%kb%ibxz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxz_c',routine)
    i_all=-product(shape(bounds%kb%ibxy_c))*kind(bounds%kb%ibxy_c)
    deallocate(bounds%kb%ibxy_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_c',routine)
    i_all=-product(shape(bounds%kb%ibyz_f))*kind(bounds%kb%ibyz_f)
    deallocate(bounds%kb%ibyz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_f',routine)
    i_all=-product(shape(bounds%kb%ibxz_f))*kind(bounds%kb%ibxz_f)
    deallocate(bounds%kb%ibxz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxz_f',routine)
    i_all=-product(shape(bounds%kb%ibxy_f))*kind(bounds%kb%ibxy_f)
    deallocate(bounds%kb%ibxy_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_f',routine)

    i_all=-product(shape(bounds%sb%ibzzx_c))*kind(bounds%sb%ibzzx_c)
    deallocate(bounds%sb%ibzzx_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibzzx_c',routine)
    i_all=-product(shape(bounds%sb%ibyyzz_c))*kind(bounds%sb%ibyyzz_c)
    deallocate(bounds%sb%ibyyzz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_c',routine)
    i_all=-product(shape(bounds%sb%ibxy_ff))*kind(bounds%sb%ibxy_ff)
    deallocate(bounds%sb%ibxy_ff,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_ff',routine)
    i_all=-product(shape(bounds%sb%ibzzx_f))*kind(bounds%sb%ibzzx_f)
    deallocate(bounds%sb%ibzzx_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibzzx_f',routine)
    i_all=-product(shape(bounds%sb%ibyyzz_f))*kind(bounds%sb%ibyyzz_f)
    deallocate(bounds%sb%ibyyzz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_f',routine)

    i_all=-product(shape(bounds%gb%ibzxx_c))*kind(bounds%gb%ibzxx_c)
    deallocate(bounds%gb%ibzxx_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibzxx_c',routine)
    i_all=-product(shape(bounds%gb%ibxxyy_c))*kind(bounds%gb%ibxxyy_c)
    deallocate(bounds%gb%ibxxyy_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxxyy_c',routine)
    i_all=-product(shape(bounds%gb%ibyz_ff))*kind(bounds%gb%ibyz_ff)
    deallocate(bounds%gb%ibyz_ff,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_ff',routine)
    i_all=-product(shape(bounds%gb%ibzxx_f))*kind(bounds%gb%ibzxx_f)
    deallocate(bounds%gb%ibzxx_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibzxx_f',routine)
    i_all=-product(shape(bounds%gb%ibxxyy_f))*kind(bounds%gb%ibxxyy_f)
    deallocate(bounds%gb%ibxxyy_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxxyy_f',routine)

    i_all=-product(shape(bounds%ibyyzz_r))*kind(bounds%ibyyzz_r)
    deallocate(bounds%ibyyzz_r,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_r',routine)

  end subroutine deallocate_bounds

end module module_types
