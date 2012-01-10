subroutine memocc_report()
  use m_profiling, only: mreport => memocc_report

  call mreport()
end subroutine memocc_report

subroutine deallocate_double(array)
  use module_base
  implicit none

  double precision, dimension(:,:), pointer :: array
  integer :: i_all, i_stat

  if (associated(array)) then
     i_all=-product(shape(array))*kind(array)
     deallocate(array,stat=i_stat)
     call memocc(i_stat,i_all,'array',"deallocate_double")
  end if
end subroutine deallocate_double

subroutine glr_new(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  allocate(glr)
end subroutine glr_new
subroutine glr_free(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  deallocate(glr)
end subroutine glr_free
subroutine glr_get_n(glr, n)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  integer, dimension(3), intent(out) :: n

  n(1) = glr%d%n1
  n(2) = glr%d%n2
  n(3) = glr%d%n3
end subroutine glr_get_n

subroutine inputs_new(in)
  use module_types
  implicit none
  type(input_variables), pointer :: in

  allocate(in)
  call default_input_variables(in)
end subroutine inputs_new
subroutine inputs_free(in)
  use module_types
  implicit none
  type(input_variables), pointer :: in

  deallocate(in)
end subroutine inputs_free
subroutine inputs_set_radical(in, rad, ln)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  integer, intent(in) :: ln
  character(len = ln), intent(in) :: rad

  call standard_inputfile_names(in, rad)
end subroutine inputs_set_radical
subroutine inputs_parse_params(in, iproc)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  integer, intent(in) :: iproc

  ! Parse all values independant from atoms.
  call perf_input_variables(iproc,.false.,trim(in%file_perf),in)
  call dft_input_variables_new(iproc,.false.,trim(in%file_dft),in)
  call mix_input_variables_new(iproc,.false.,trim(in%file_mix),in)
  call geopt_input_variables_new(iproc,.false.,trim(in%file_geopt),in)
  call tddft_input_variables_new(iproc,.false.,trim(in%file_tddft),in)
  call sic_input_variables_new(iproc,.false.,trim(in%file_sic),in)
end subroutine inputs_parse_params
subroutine inputs_get_dft(in, hx, hy, hz, crmult, frmult, ixc, chg, efield, nspin, mpol, &
     & gnrm, itermax, nrepmax, ncong, idsx, dispcorr, inpsi, outpsi, outgrid, &
     & rbuf, ncongt, davidson, nvirt, nplottedvirt, sym)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  real(gp), intent(out) :: hx, hy, hz, crmult, frmult, efield(3), gnrm, rbuf
  integer, intent(out) :: ixc, chg, nspin, mpol, itermax, nrepmax, ncong, idsx, &
       & dispcorr, inpsi, outpsi, outgrid, ncongt, davidson, nvirt, nplottedvirt, sym
  
  hx = in%hx
  hy = in%hy
  hz = in%hz
  crmult = in%crmult
  frmult = in%frmult
  ixc = in%ixc
  chg = in%ncharge
  efield = in%elecfield
  nspin = in%nspin
  mpol = in%mpol
  gnrm = in%gnrm_cv
  itermax = in%itrpmax
  nrepmax = in%nrepmax
  ncong = in%ncong
  idsx = in%idsx
  dispcorr = in%dispersion
  inpsi = in%inputPsiId
  outpsi = in%output_wf_format
  outgrid = in%output_denspot
  rbuf = in%rbuf
  ncongt = in%ncongt
  davidson = in%norbv
  nvirt = in%nvirt
  nplottedvirt = in%nplot
  if (in%disableSym) then
     sym = 1
  else
     sym = 0
  end if
end subroutine inputs_get_dft
subroutine inputs_get_files(in, files)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(out) :: files

  files = in%files
end subroutine inputs_get_files
