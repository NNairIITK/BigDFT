! TODO: Write me please !
module xc_f90_types_m
  implicit none

  type xc_f90_pointer_t
     integer :: i
  end type xc_f90_pointer_t

end module xc_f90_types_m

module libxc_funcs_m
  implicit none

  integer, parameter :: XC_FAMILY_LDA = 0, XC_FAMILY_GGA = 1, XC_FAMILY_HYB_GGA = 2
  integer, parameter :: XC_EXCHANGE = 0, XC_CORRELATION = 1, XC_EXCHANGE_CORRELATION = 2

  integer, parameter :: XC_HYB_GGA_XC_PBEH = 0

  integer, parameter :: XC_FLAGS_HAVE_EXC = 0, XC_FLAGS_HAVE_FXC = 1
end module libxc_funcs_m

module xc_f90_lib_m
  use xc_f90_types_m

  implicit none

  interface
     integer function xc_f90_info_flags(xc_info)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: xc_info
     end function xc_f90_info_flags

     integer function xc_f90_info_kind(xc_info)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: xc_info
     end function xc_f90_info_kind

     subroutine xc_f90_info_name(xc_info, name)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: xc_info
       character(len = *), intent(out) :: name
     end subroutine xc_f90_info_name

     subroutine xc_f90_info_refs(xc_info, ii, str, name)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: xc_info, str
       integer, intent(inout) :: ii
       character(len = *), intent(out) :: name
     end subroutine xc_f90_info_refs

     integer function xc_f90_family_from_id(xc_id)
       integer, intent(in) :: xc_id
     end function xc_f90_family_from_id

     subroutine xc_f90_func_init(conf, info, ixc, nspden)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(out) :: conf, info
       integer, intent(in) :: ixc, nspden
     end subroutine xc_f90_func_init

     subroutine xc_f90_func_end(conf)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(inout) :: conf
     end subroutine xc_f90_func_end

     subroutine xc_f90_lda_vxc(conf,npts,rho,vxc)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: conf
       integer, intent(in) :: npts
       double precision, intent(in) :: rho
       double precision, intent(out) :: vxc
     end subroutine xc_f90_lda_vxc

     subroutine xc_f90_lda_fxc(conf,npts,rho,fxc)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: conf
       integer, intent(in) :: npts
       double precision, intent(in) :: rho
       double precision, intent(out) :: fxc
     end subroutine xc_f90_lda_fxc

     subroutine xc_f90_lda_exc_vxc(conf,npts,rho,exc,vxc)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: conf
       integer, intent(in) :: npts
       double precision, intent(in) :: rho
       double precision, intent(out) :: exc,vxc
     end subroutine xc_f90_lda_exc_vxc

     subroutine xc_f90_gga_vxc(conf,npts,rho,sigma,vxc,vsigma)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: conf
       integer, intent(in) :: npts
       double precision, intent(in) :: rho, sigma
       double precision, intent(out) :: vxc, vsigma
     end subroutine xc_f90_gga_vxc

     subroutine xc_f90_gga_fxc(conf,npts,rho,sigma,fxc,fxc2,fxc3)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: conf
       integer, intent(in) :: npts
       double precision, intent(in) :: rho,sigma
       double precision, intent(out) :: fxc,fxc2,fxc3
     end subroutine xc_f90_gga_fxc

     subroutine xc_f90_gga_exc_vxc(conf,npts,rho,sigma,exc,vxc,vsigma)
       use xc_f90_types_m
       type(xc_f90_pointer_t), intent(in) :: conf
       integer, intent(in) :: npts
       double precision, intent(in) :: rho,sigma
       double precision, intent(out) :: exc,vxc,vsigma
     end subroutine xc_f90_gga_exc_vxc
  end interface

end module xc_f90_lib_m

integer function xc_f90_info_flags(xc_info)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: xc_info

  xc_f90_info_flags = 0
  write(0, *) "No LibXC support at compile time, abort."
  stop
end function xc_f90_info_flags

integer function xc_f90_info_kind(xc_info)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: xc_info

  xc_f90_info_kind = 0
  write(0, *) "No LibXC support at compile time, abort."
  stop
end function xc_f90_info_kind

subroutine xc_f90_info_name(xc_info, name)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: xc_info
  character(len = *), intent(out) :: name

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_info_name

subroutine xc_f90_info_refs(xc_info, ii, str, name)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: xc_info, str
  integer, intent(inout) :: ii
  character(len = *), intent(out) :: name

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_info_refs

integer function xc_f90_family_from_id(xc_id)
  integer, intent(in) :: xc_id

  xc_f90_family_from_id = 0
  write(0, *) "No LibXC support at compile time, abort."
  stop
end function xc_f90_family_from_id

subroutine xc_f90_func_init(conf, info, ixc, nspden)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(out) :: conf, info
  integer, intent(in) :: ixc, nspden

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_func_init

subroutine xc_f90_func_end(conf)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(inout) :: conf

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_func_end

subroutine xc_f90_lda_vxc(conf,npts,rho,vxc)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: conf
  integer, intent(in) :: npts
  double precision, intent(in) :: rho
  double precision, intent(out) :: vxc

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_lda_vxc

subroutine xc_f90_lda_fxc(conf,npts,rho,fxc)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: conf
  integer, intent(in) :: npts
  double precision, intent(in) :: rho
  double precision, intent(out) :: fxc

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_lda_fxc

subroutine xc_f90_lda_exc_vxc(conf,npts,rho,exc,vxc)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: conf
  integer, intent(in) :: npts
  double precision, intent(in) :: rho
  double precision, intent(out) :: exc,vxc

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_lda_exc_vxc

subroutine xc_f90_gga_vxc(conf,npts,rho,sigma,vxc,vsigma)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: conf
  integer, intent(in) :: npts
  double precision, intent(in) :: rho, sigma
  double precision, intent(out) :: vxc, vsigma

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_gga_vxc

subroutine xc_f90_gga_fxc(conf,npts,rho,sigma,fxc,fxc2,fxc3)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: conf
  integer, intent(in) :: npts
  double precision, intent(in) :: rho,sigma
  double precision, intent(out) :: fxc,fxc2,fxc3

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_gga_fxc

subroutine xc_f90_gga_exc_vxc(conf,npts,rho,sigma,exc,vxc,vsigma)
  use xc_f90_types_m
  type(xc_f90_pointer_t), intent(in) :: conf
  integer, intent(in) :: npts
  double precision, intent(in) :: rho,sigma
  double precision, intent(out) :: exc,vxc,vsigma

  write(0, *) "No LibXC support at compile time, abort."
  stop
end subroutine xc_f90_gga_exc_vxc
