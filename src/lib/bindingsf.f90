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
  nullify(glr%wfd%keyg)
  nullify(glr%wfd%keyv)

  nullify(glr%bounds%kb%ibyz_f)
  nullify(glr%bounds%kb%ibyz_c)
end subroutine glr_new
subroutine glr_free(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  call deallocate_lr(glr, "glr_free")
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
subroutine glr_set_wave_descriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
      &   crmult,frmult,Glr)
   use module_base
   use module_types
   use module_interfaces
   implicit none
   !Arguments
   type(atoms_data), intent(in) :: atoms
   integer, intent(in) :: iproc
   real(gp), intent(in) :: hx,hy,hz,crmult,frmult
   real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
   real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
   type(locreg_descriptors), intent(inout) :: Glr

   call createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
      &   crmult,frmult,Glr)
end subroutine glr_set_wave_descriptors

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

  call free_input_variables(in)
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
subroutine inputs_parse_params(in, iproc, dump)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  integer, intent(in) :: iproc
  logical, intent(in) :: dump

  ! Parse all values independant from atoms.
  call perf_input_variables(iproc,dump,trim(in%file_perf),in)
  call dft_input_variables_new(iproc,dump,trim(in%file_dft),in)
  call mix_input_variables_new(iproc,dump,trim(in%file_mix),in)
  call geopt_input_variables_new(iproc,dump,trim(in%file_geopt),in)
  call tddft_input_variables_new(iproc,dump,trim(in%file_tddft),in)
  call sic_input_variables_new(iproc,dump,trim(in%file_sic),in)
end subroutine inputs_parse_params
subroutine inputs_parse_add(in, atoms, iproc, dump)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: iproc
  logical, intent(in) :: dump

  ! Read k-points input variables (if given)
  call kpt_input_variables_new(iproc,dump,trim(in%file_kpt),in,atoms)
end subroutine inputs_parse_add
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
END SUBROUTINE inputs_get_dft
subroutine inputs_get_mix(in, iscf, itrpmax, norbsempty, occopt, alphamix, rpnrm_cv, &
     & gnrm_startmix, Tel, alphadiis)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(out) :: iscf, itrpmax, norbsempty, occopt
  real(gp), intent(out) :: alphamix, rpnrm_cv, gnrm_startmix, Tel, alphadiis
  
  iscf = in%iscf
  itrpmax = in%itrpmax
  norbsempty = in%norbsempty
  occopt = in%occopt

  alphamix = in%alphamix
  rpnrm_cv = in%rpnrm_cv
  gnrm_startmix = in%gnrm_startmix
  Tel = in%Tel
  alphadiis = in%alphadiis
END SUBROUTINE inputs_get_mix
subroutine inputs_get_files(in, files)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(out) :: files

  files = in%files
END SUBROUTINE inputs_get_files

subroutine orbs_new(orbs)
  use module_types
  implicit none
  type(orbitals_data), pointer :: orbs

  allocate(orbs)
END SUBROUTINE orbs_new
subroutine orbs_free(orbs)
  use module_types
  implicit none
  type(orbitals_data), pointer :: orbs

  call deallocate_orbs(orbs,"orbs_free")
  deallocate(orbs)
END SUBROUTINE orbs_free
subroutine orbs_get_dimensions(orbs, norb, norbp, norbu, norbd, nspin, nspinor, npsidim, &
     & nkpts, nkptsp, isorb, iskpts)
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: norb, norbp, norbu, norbd, nspin, nspinor, npsidim, &
     & nkpts, nkptsp, isorb, iskpts
  
  norb = orbs%norb
  norbp = orbs%norbp
  norbu = orbs%norbu
  norbd = orbs%norbd
  nspin = orbs%nspin
  nspinor = orbs%nspinor
  npsidim = orbs%npsidim
  nkpts = orbs%nkpts
  nkptsp = orbs%nkptsp
  isorb = orbs%isorb
  iskpts = orbs%iskpts
END SUBROUTINE orbs_get_dimensions
