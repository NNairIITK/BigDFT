subroutine memocc_report()
  use m_profiling, only: mreport => memocc_report

  call mreport()
end subroutine memocc_report

subroutine f90_pointer_1D_init(pt_c, size_c)
  implicit none
  double precision, intent(in) :: pt_c
  integer, intent(in) :: size_c

  double precision, dimension(:), pointer :: pt_f
  interface
     subroutine inquire_pointer1(pt_c, pt_f, size_c)
       double precision, dimension(:), pointer :: pt_f
       double precision, intent(in) :: pt_c
       integer, intent(in) :: size_c
     end subroutine inquire_pointer1
  end interface

  nullify(pt_f)
  call inquire_pointer1(pt_c, pt_f, size_c)
end subroutine f90_pointer_1D_init

subroutine f90_pointer_2D_init(pt_c, size_c)
  implicit none
  double precision, intent(in) :: pt_c
  integer, intent(in) :: size_c

  double precision, dimension(:,:), pointer :: pt_f
  interface
     subroutine inquire_pointer2(pt_c, pt_f, size_c)
       double precision, dimension(:,:), pointer :: pt_f
       double precision, intent(in) :: pt_c
       integer, intent(in) :: size_c
     end subroutine inquire_pointer2
  end interface

  nullify(pt_f)
  call inquire_pointer2(pt_c, pt_f, size_c)
end subroutine f90_pointer_2D_init

subroutine f90_pointer_3D_init(pt_c, size_c)
  implicit none
  double precision, intent(in) :: pt_c
  integer, intent(in) :: size_c

  double precision, dimension(:,:,:), pointer :: pt_f
  interface
     subroutine inquire_pointer3(pt_c, pt_f, size_c)
       double precision, dimension(:,:,:), pointer :: pt_f
       double precision, intent(in) :: pt_c
       integer, intent(in) :: size_c
     end subroutine inquire_pointer3
  end interface

  nullify(pt_f)
  call inquire_pointer3(pt_c, pt_f, size_c)
end subroutine f90_pointer_3D_init

subroutine f90_pointer_4D_init(pt_c, size_c)
  implicit none
  double precision, intent(in) :: pt_c
  integer, intent(in) :: size_c

  double precision, dimension(:,:,:,:), pointer :: pt_f
  interface
     subroutine inquire_pointer4(pt_c, pt_f, size_c)
       double precision, dimension(:,:,:,:), pointer :: pt_f
       double precision, intent(in) :: pt_c
       integer, intent(in) :: size_c
     end subroutine inquire_pointer4
  end interface

  nullify(pt_f)
  call inquire_pointer4(pt_c, pt_f, size_c)
end subroutine f90_pointer_4D_init

subroutine f90_pointer_5D_init(pt_c, size_c)
  implicit none
  double precision, intent(in) :: pt_c
  integer, intent(in) :: size_c

  double precision, dimension(:,:,:,:,:), pointer :: pt_f
  interface
     subroutine inquire_pointer5(pt_c, pt_f, size_c)
       double precision, dimension(:,:,:,:,:), pointer :: pt_f
       double precision, intent(in) :: pt_c
       integer, intent(in) :: size_c
     end subroutine inquire_pointer5
  end interface

  nullify(pt_f)
  call inquire_pointer5(pt_c, pt_f, size_c)
end subroutine f90_pointer_5D_init

subroutine createKernel(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,kernel,wrtmsg)
  use Poisson_Solver, only: ck => createKernel
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,itype_scf,iproc,nproc
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), pointer :: kernel(:)
  logical, intent(in) :: wrtmsg

  call ck(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,kernel,wrtmsg)
end subroutine createKernel

subroutine deallocate_double_1D(array)
  use module_base
  implicit none

  double precision, dimension(:), pointer :: array
  integer :: i_all, i_stat

  if (associated(array)) then
     i_all=-product(shape(array))*kind(array)
     deallocate(array,stat=i_stat)
     call memocc(i_stat,i_all,'array',"deallocate_double")
  end if
end subroutine deallocate_double_1D
subroutine deallocate_double_2D(array)
  use module_base
  implicit none

  double precision, dimension(:,:), pointer :: array
  integer :: i_all, i_stat

  if (associated(array)) then
     i_all=-product(shape(array))*kind(array)
     deallocate(array,stat=i_stat)
     call memocc(i_stat,i_all,'array',"deallocate_double")
  end if
end subroutine deallocate_double_2D

subroutine glr_new(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  allocate(glr)
  call nullify_locreg_descriptors(glr)
end subroutine glr_new
subroutine glr_free(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  call deallocate_lr(glr, "glr_free")
  deallocate(glr)
end subroutine glr_free
subroutine glr_get_dimensions(glr, geocode, n, ni)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  character, intent(out) :: geocode(1)
  integer, dimension(3), intent(out) :: n, ni

  write(geocode, "(A1)") glr%geocode
  n(1) = glr%d%n1
  n(2) = glr%d%n2
  n(3) = glr%d%n3
  ni(1) = glr%d%n1i
  ni(2) = glr%d%n2i
  ni(3) = glr%d%n3i
end subroutine glr_get_dimensions
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
  character, intent(in) :: rad(ln)

  character(len = 1024) :: rad_
  integer :: i

  write(rad_, "(A)") " "
  do i = 1, ln
     write(rad_(i:i), "(A1)") rad(i)
  end do
  call standard_inputfile_names(in, rad_)
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
subroutine inputs_get_geopt(in, geopt_approach, ncount_cluster_x, frac_fluct, forcemax, &
     & randdis, betax, history, ionmov, dtion, strtarget, qmass)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  character(len = 10), intent(out) :: geopt_approach
  integer, intent(out) :: ncount_cluster_x, history, ionmov
  real(gp), intent(out) :: frac_fluct, forcemax, randdis, betax, dtion, strtarget(6)
  real(gp), pointer :: qmass(:)
  
  geopt_approach = in%geopt_approach
  ncount_cluster_x = in%ncount_cluster_x
  frac_fluct = in%frac_fluct
  forcemax = in%forcemax
  randdis = in%randdis
  betax = in%betax
  history = in%history
  ionmov = in%ionmov
  dtion = in%dtion
  strtarget(:) = in%strtarget(:)
  if (associated(in%qmass)) then
     qmass => in%qmass
  else
     nullify(qmass)
  end if
END SUBROUTINE inputs_get_geopt
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
  npsidim = max(orbs%npsidim_orbs,orbs%npsidim_comp)
  nkpts = orbs%nkpts
  nkptsp = orbs%nkptsp
  isorb = orbs%isorb
  iskpts = orbs%iskpts
END SUBROUTINE orbs_get_dimensions

subroutine proj_new(nlpspd)
  use module_types
  implicit none
  type(nonlocal_psp_descriptors), pointer :: nlpspd

  allocate(nlpspd)
END SUBROUTINE proj_new
subroutine proj_free(nlpspd, proj)
  use module_types
  use m_profiling
  implicit none
  type(nonlocal_psp_descriptors), pointer :: nlpspd
  real(kind=8), dimension(:), pointer :: proj

  integer :: i_stat, i_all

  call deallocate_proj_descr(nlpspd,"proj_free")
  i_all=-product(shape(proj))*kind(proj)
  deallocate(proj,stat=i_stat)
  call memocc(i_stat,i_all,'proj',"proj_free")
END SUBROUTINE proj_free
subroutine proj_get_dimensions(nlpspd, nproj, nprojel)
  use module_types
  implicit none
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  integer, intent(out) :: nproj, nprojel
  
  nproj = nlpspd%nproj
  nprojel = nlpspd%nprojel
END SUBROUTINE proj_get_dimensions

subroutine denspot_new(denspotd, rhodsc)
  use module_types
  implicit none
  type(denspot_distribution), pointer :: denspotd
  type(rho_descriptors), pointer :: rhodsc

  allocate(denspotd)
  allocate(rhodsc)
END SUBROUTINE denspot_new
subroutine denspot_free(denspotd, rhodsc, pot_ion, rhopot, potxc, rhocore)
  use module_types
  use m_profiling
  implicit none
  type(denspot_distribution), pointer :: denspotd
  type(rho_descriptors), pointer :: rhodsc
  real(dp), dimension(:), pointer :: pot_ion, rhopot
  real(kind = 8), dimension(:), pointer :: rhocore
  real(kind = 8), dimension(:,:,:,:), pointer :: potxc

  integer :: i_stat, i_all

  call deallocate_rho_descriptors(rhodsc, "denspot_free")
  deallocate(rhodsc)

  call deallocate_denspot_distribution(denspotd, "denspot_free")
  deallocate(denspotd)

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',"denspot_free")

  i_all=-product(shape(rhopot))*kind(rhopot)
  deallocate(rhopot,stat=i_stat)
  call memocc(i_stat,i_all,'rhopot',"denspot_free")

  i_all=-product(shape(rhocore))*kind(rhocore)
  deallocate(rhocore,stat=i_stat)
  call memocc(i_stat,i_all,'rhocore',"denspot_free")

  i_all=-product(shape(potxc))*kind(potxc)
  deallocate(potxc,stat=i_stat)
  call memocc(i_stat,i_all,'potxc',"denspot_free")
END SUBROUTINE denspot_free
