!!****f* BigDFT/write_etsf_density
!! FUNCTION
!!   Write a field in the ISF basis in the ETSF format
!! SOURCE
!!
subroutine write_etsf_density(filename,message,at,rxyz,n1,n2,n3,n1i,n2i,n3i,hxh,hyh,hzh,&
     x)
  use module_base
  use module_types

  implicit none
  character(len=*), intent(in) :: filename,message
  integer, intent(in) :: n1,n2,n3,n1i,n2i,n3i
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  real(wp), dimension(n1i,n2i,n3i), intent(in) :: x
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables

  write(0, "(A)") "Illegal call to write_etsf_density(), not compiled with ETSF_IO support."
  stop
END SUBROUTINE write_etsf_density
!!***

!!****f* BigDFT/read_etsf
!! FUNCTION
!!   Read a field in the ISF basis in the ETSF format
!!
!! SOURCE
!!
subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz)
  use module_base
  use module_types
  use etsf_io
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:), pointer :: rho
  real(gp), dimension(:,:), pointer, optional :: rxyz
  integer, intent(out), optional ::  nat

  write(0, "(A)") "Illegal call to read_etsf(), not compiled with ETSF_IO support."
  stop
END SUBROUTINE read_etsf
!!***
