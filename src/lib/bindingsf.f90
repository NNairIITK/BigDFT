subroutine memocc_report()
  use m_profiling, only: mreport => memocc_report
  implicit none
  call mreport()
end subroutine memocc_report

subroutine memocc_verbose()
  use m_profiling, only: mstate => memocc_set_state
  implicit none
  call mstate(2)
end subroutine memocc_verbose

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
end subroutine glr_new
subroutine glr_init(glr, d)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout), target :: glr
  type(grid_dimensions), pointer :: d

  call nullify_locreg_descriptors(glr)
  d => glr%d
end subroutine glr_init
subroutine glr_get_data(glr, d)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout), target :: glr
  type(grid_dimensions), pointer :: d

  d => glr%d
end subroutine glr_get_data
subroutine glr_free(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  deallocate(glr)
end subroutine glr_free
subroutine glr_empty(glr)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout) :: glr

  call deallocate_locreg_descriptors(glr, "glr_empty")
end subroutine glr_empty
subroutine glr_get_dimensions(glr , n, ni)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  integer, dimension(3), intent(out) :: n, ni

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
subroutine lzd_new(lzd)
  use module_types
  implicit none
  type(local_zone_descriptors), pointer :: lzd

  allocate(lzd)
end subroutine lzd_new
subroutine lzd_init(lzd, glr)
  use module_types
  implicit none
  type(local_zone_descriptors), target, intent(inout) :: lzd
  type(locreg_descriptors), pointer :: glr

  call nullify_local_zone_descriptors(lzd)
  glr => lzd%glr
end subroutine lzd_init
subroutine lzd_get_data(lzd, glr)
  use module_types
  implicit none
  type(local_zone_descriptors), target, intent(inout) :: lzd
  type(locreg_descriptors), pointer :: glr

  glr => lzd%glr
end subroutine lzd_get_data
subroutine lzd_free(lzd)
  use module_types
  implicit none
  type(local_zone_descriptors), pointer :: lzd

  call deallocate_local_zone_descriptors(lzd, "lzd_free")
  deallocate(lzd)
end subroutine lzd_free
subroutine lzd_empty(lzd)
  use module_types
  implicit none
  type(local_zone_descriptors), intent(inout) :: lzd

  call deallocate_Lzd_except_Glr(lzd, "lzd_empty")
END SUBROUTINE lzd_empty
subroutine lzd_get_hgrids(Lzd, hgrids)
  use module_base
  use module_types
  implicit none
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), intent(out) :: hgrids(3)
  !initial values
  hgrids = Lzd%hgrids
END SUBROUTINE lzd_get_hgrids

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
subroutine inputs_set_radical(in, nproc, rad, ln)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  integer, intent(in) :: ln, nproc
  character, intent(in) :: rad(ln)

  character(len = 1024) :: rad_
  integer :: i

  write(rad_, "(A)") " "
  do i = 1, ln
     write(rad_(i:i), "(A1)") rad(i)
  end do
  call standard_inputfile_names(in, rad_, nproc)
end subroutine inputs_set_radical

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
subroutine orbs_init(orbs)
  use module_types
  implicit none
  type(orbitals_data), intent(inout) :: orbs

  call nullify_orbitals_data(orbs)
END SUBROUTINE orbs_init
subroutine orbs_free(orbs)
  use module_types
  use m_profiling
  implicit none
  type(orbitals_data), pointer :: orbs

  deallocate(orbs)
END SUBROUTINE orbs_free
subroutine orbs_empty(orbs)
  use module_types
  use m_profiling
  implicit none
  type(orbitals_data), intent(inout) :: orbs

  integer :: i_all, i_stat

  if (associated(orbs%norb_par)) then
     call deallocate_orbs(orbs,"orbs_empty")
  end if
  if (associated(orbs%eval)) then
     i_all=-product(shape(orbs%eval))*kind(orbs%eval)
     deallocate(orbs%eval,stat=i_stat)
     call memocc(i_stat,i_all,'orbs%eval',"orbs_empty")
  end if
END SUBROUTINE orbs_empty
subroutine orbs_comm_new(comms)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  type(communications_arrays), pointer :: comms

  allocate(comms)
  nullify(comms%nvctr_par)
end subroutine orbs_comm_new
subroutine orbs_comm_init(comms, orbs, lr, iproc, nproc)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(inout) :: comms

  call orbitals_communicators(iproc,nproc,lr,orbs,comms)
end subroutine orbs_comm_init
subroutine orbs_comm_free(comms)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  type(communications_arrays), pointer :: comms

  deallocate(comms)
end subroutine orbs_comm_free
subroutine orbs_comm_empty(comms)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  type(communications_arrays), intent(inout) :: comms

  if (associated(comms%nvctr_par)) then
     call deallocate_comms(comms,"orbs_comm_empty")
  end if
end subroutine orbs_comm_empty
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
subroutine orbs_get_eval(orbs, eval)
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(wp), dimension(:), pointer :: eval
  
  eval => orbs%eval
END SUBROUTINE orbs_get_eval
subroutine orbs_get_occup(orbs, occup)
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(gp), dimension(:), pointer :: occup
  
  occup => orbs%occup
END SUBROUTINE orbs_get_occup
subroutine orbs_get_kpts(orbs, kpts)
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(gp), dimension(:,:), pointer :: kpts
  
  kpts => orbs%kpts
END SUBROUTINE orbs_get_kpts
subroutine orbs_get_kwgts(orbs, kwgts)
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(gp), dimension(:), pointer :: kwgts
  
  kwgts => orbs%kwgts
END SUBROUTINE orbs_get_kwgts

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

subroutine localfields_new(self, denspotd, rhod, dpbox)
  use module_types
  implicit none
  double precision, intent(in) :: self
  type(DFT_local_fields), pointer :: denspotd
  type(denspot_distribution), pointer :: dpbox
  type(rho_descriptors), pointer :: rhod

  allocate(denspotd)
  rhod => denspotd%rhod
  dpbox => denspotd%dpbox
  denspotd%c_obj = self
END SUBROUTINE localfields_new
subroutine localfields_get_data(denspotd, rhod, dpcom)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in), target :: denspotd
  type(denspot_distribution), pointer :: dpcom
  type(rho_descriptors), pointer :: rhod

  rhod => denspotd%rhod
  dpcom => denspotd%dpcom
END SUBROUTINE localfields_get_data
subroutine localfields_free(denspotd)
  use module_types
  use m_profiling
  implicit none
  type(DFT_local_fields), pointer :: denspotd
  
  character(len = *), parameter :: subname = "localfields_free"
  integer :: i_stat, i_all

  call deallocate_rho_descriptors(denspotd%rhod, subname)
  call deallocate_denspot_distribution(denspotd%dpbox, subname)
  
  if (associated(denspotd%V_ext)) then
     i_all=-product(shape(denspotd%V_ext))*kind(denspotd%V_ext)
     deallocate(denspotd%V_ext,stat=i_stat)
     call memocc(i_stat,i_all,'denspotd%V_ext',subname)
  end if

!!$  if (associated(denspotd%pkernelseq)) then
!!$     i_all=-product(shape(denspotd%pkernelseq))*kind(denspotd%pkernelseq)
!!$     deallocate(denspotd%pkernelseq,stat=i_stat)
!!$     call memocc(i_stat,i_all,'kernelseq',subname)
!!$  end if

  if (associated(denspotd%pkernel)) then
     i_all=-product(shape(denspotd%pkernel))*kind(denspotd%pkernel)
     deallocate(denspotd%pkernel,stat=i_stat)
     call memocc(i_stat,i_all,'kernel',subname)
  end if

  if (associated(denspotd%rhov)) then
     i_all=-product(shape(denspotd%rhov))*kind(denspotd%rhov)
     deallocate(denspotd%rhov,stat=i_stat)
     call memocc(i_stat,i_all,'denspotd%rhov',subname)
  end if

  if (associated(denspotd%V_XC)) then
     i_all=-product(shape(denspotd%V_XC))*kind(denspotd%V_XC)
     deallocate(denspotd%V_XC,stat=i_stat)
     call memocc(i_stat,i_all,'denspotd%V_XC',subname)
  end if

  if(associated(denspotd%rho_C)) then
     i_all=-product(shape(denspotd%rho_C))*kind(denspotd%rho_C)
     deallocate(denspotd%rho_C,stat=i_stat)
     call memocc(i_stat,i_all,'denspotd%rho_C',subname)
  end if

  deallocate(denspotd)
END SUBROUTINE localfields_free
subroutine localfields_copy_metadata(denspot, rhov_is, hgrid, psoffset)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(out) :: rhov_is
  real(gp), intent(out) :: hgrid(3)
  real(dp), intent(out) :: psoffset

  rhov_is = denspot%rhov_is
  hgrid = denspot%dpbox%hgrids
  psoffset = denspot%psoffset
END SUBROUTINE localfields_copy_metadata
subroutine localfields_get_rhov(denspot, rhov)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: rhov

  rhov => denspot%rhov
END SUBROUTINE localfields_get_rhov
subroutine localfields_get_v_ext(denspot, v_ext)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(wp), dimension(:,:,:,:), pointer :: v_ext

  v_ext => denspot%v_ext
END SUBROUTINE localfields_get_v_ext
subroutine localfields_get_v_xc(denspot, v_xc)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(wp), dimension(:,:,:,:), pointer :: v_xc

  v_xc => denspot%v_xc
END SUBROUTINE localfields_get_v_xc
subroutine localfields_get_pkernel(denspot, pkernel)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: pkernel

  pkernel => denspot%pkernel
END SUBROUTINE localfields_get_pkernel
subroutine localfields_get_pkernelseq(denspot, pkernelseq)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: pkernelseq

  pkernelseq => denspot%pkernelseq
END SUBROUTINE localfields_get_pkernelseq
subroutine localfields_get_rho_work(denspot, rho)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: rho

  rho => denspot%rho_work
END SUBROUTINE localfields_get_rho_work
subroutine localfields_get_pot_work(denspot, pot)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: pot

  pot => denspot%pot_work
END SUBROUTINE localfields_get_pot_work

subroutine localfields_full_density(denspot, rho_full, iproc, new)
  use module_base
  use module_types
  use m_profiling
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc
  integer, intent(out) :: new
  real(gp), dimension(:), pointer :: rho_full

  character(len = *), parameter :: subname = "localfields_full_density"
  integer :: i_stat, nslice, ierr, irhodim, irhoxcsh

  new = 0
  nslice = max(denspot%dpcom%ndimpot, 1)
  if (nslice < denspot%dpcom%ndimgrid) then
     if (iproc == 0) then
        !allocate full density in pot_ion array
        allocate(rho_full(denspot%dpcom%ndimgrid*denspot%dpcom%nrhodim+ndebug),stat=i_stat)
        call memocc(i_stat,rho_full,'rho_full',subname)
        new = 1
        
        ! Ask to gather density to other procs.
        call MPI_BCAST(0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if

     if (denspot%dpcom%ndimrhopot > 0) then
        irhoxcsh = nslice / denspot%dpcom%n3p * denspot%dpcom%i3xcsh
     else
        irhoxcsh = 0
     end if     
     do irhodim = 1, denspot%dpcom%nrhodim, 1
        call MPI_GATHERV(denspot%rhov(nslice * (irhodim - 1) + irhoxcsh + 1),&
             nslice,mpidtypd,rho_full(denspot%dpcom%ndimgrid * (irhodim - 1) + 1),&
             denspot%dpcom%ngatherarr(0,1),denspot%dpcom%ngatherarr(0,2),&
             mpidtypd,0,MPI_COMM_WORLD,ierr)
     end do
  else
     rho_full => denspot%rhov
  end if
END SUBROUTINE localfields_full_density
subroutine localfields_full_v_ext(denspot, pot_full, iproc, new)
  use module_base
  use module_types
  use m_profiling
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc
  integer, intent(out) :: new
  real(gp), pointer :: pot_full(:)

  character(len = *), parameter :: subname = "localfields_full_potential"
  integer :: i_stat, ierr

  new = 0
  if (denspot%dpcom%ndimpot < denspot%dpcom%ndimgrid) then
     if (iproc == 0) then
        !allocate full density in pot_ion array
        allocate(pot_full(denspot%dpcom%ndimgrid+ndebug),stat=i_stat)
        call memocc(i_stat,pot_full,'pot_full',subname)
        new = 1
      
        ! Ask to gather density to other procs.
        call MPI_BCAST(1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if

     call MPI_GATHERV(denspot%v_ext(1,1,1,1),max(denspot%dpcom%ndimpot, 1),&
          mpidtypd,pot_full(1),denspot%dpcom%ngatherarr(0,1),&
          denspot%dpcom%ngatherarr(0,2),mpidtypd,0,MPI_COMM_WORLD,ierr)
  else
     pot_full => denspot%rhov
  end if
END SUBROUTINE localfields_full_v_ext

subroutine gpu_new(GPU)
  use module_types
  implicit none
  type(GPU_pointers), pointer :: GPU

  allocate(GPU)
END SUBROUTINE gpu_new
subroutine gpu_free(GPU)
  use module_types
  implicit none
  type(GPU_pointers), pointer :: GPU

  deallocate(GPU)
END SUBROUTINE gpu_free

subroutine wf_new(self, wf, orbs, comm, lzd)
  use module_types
  implicit none
  double precision, intent(in) :: self
  type(DFT_wavefunction), pointer :: wf
  type(orbitals_data), pointer :: orbs
  type(communications_arrays), pointer :: comm
  type(local_zone_descriptors), pointer :: lzd

  allocate(wf)
  wf%c_obj = self
  call wf_init(wf)
  orbs => wf%orbs
  comm => wf%comms
  lzd => wf%Lzd
end subroutine wf_new
subroutine wf_init(wf)
  use module_types
  implicit none
  type(DFT_wavefunction), intent(inout) :: wf

  nullify(wf%psi)
  nullify(wf%hpsi)
  nullify(wf%psit)
  nullify(wf%spsi)
  nullify(wf%comms%nvctr_par)
end subroutine wf_init
subroutine wf_get_data(wf, orbs, comm, lzd)
  use module_types
  implicit none
  type(DFT_wavefunction), target, intent(in) :: wf
  type(orbitals_data), pointer :: orbs
  type(communications_arrays), pointer :: comm
  type(local_zone_descriptors), pointer :: lzd

  orbs => wf%orbs
  comm => wf%comms
  lzd => wf%Lzd
end subroutine wf_get_data
subroutine wf_empty(wf)
  use module_types
  use m_profiling
  implicit none
  type(DFT_wavefunction), intent(inout) :: wf

  integer :: i_all, i_stat

  if (associated(wf%psi)) then
     i_all=-product(shape(wf%psi))*kind(wf%psi)
     deallocate(wf%psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi', "wf_empty")
  end if
  if (associated(wf%psit)) then
     i_all=-product(shape(wf%psit))*kind(wf%psit)
     deallocate(wf%psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit', "wf_empty")
  end if
  if (associated(wf%hpsi)) then
     i_all=-product(shape(wf%hpsi))*kind(wf%hpsi)
     deallocate(wf%hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi', "wf_empty")
  end if
END SUBROUTINE wf_empty
subroutine wf_free(wf)
  use module_types
  use m_profiling
  implicit none
  type(DFT_wavefunction), pointer :: wf

  call orbs_comm_empty(wf%comms)
  call orbs_empty(wf%orbs)
  call deallocate_local_zone_descriptors(wf%lzd, "wf%lzd")
  deallocate(wf)
end subroutine wf_free
subroutine wf_get_psi(wf, psi)
  use module_types
  implicit none
  type(DFT_wavefunction), intent(in) :: wf
  double precision, intent(out) :: psi

  interface
     subroutine inquire_address1(add, pt_f)
       double precision, dimension(:), pointer :: pt_f
       double precision, intent(out) :: add
     end subroutine inquire_address1
  end interface
  call inquire_address1(psi, wf%psi)
end subroutine wf_get_psi
subroutine wf_iorbp_to_psi(psir, psi, lr)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psir
  
  character(len=*), parameter :: subname='wf_orb_to_psi'
  type(workarr_sumrho) :: w

  call initialize_work_arrays_sumrho(lr,w)

  !initialisation
  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if

  call daub_to_isf(lr,w,psi,psir)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE wf_iorbp_to_psi

subroutine orbs_get_iorbp(orbs, iorbp, iproc, ikpt, iorb, ispin, ispinor)
  use module_types
  implicit none

  integer, intent(out) :: iorbp, iproc
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: ikpt, iorb, ispin, ispinor

  integer :: iorbtot

  iorbp = (ikpt - 1) * (orbs%nspinor * orbs%norb)
  if (ispin == 1) iorbp = iorbp + (iorb - 1) * orbs%nspinor
  if (ispin == 2) iorbp = iorbp + orbs%norbu * orbs%nspinor + (iorb - 1) * orbs%nspinor
  iorbp = iorbp + ispinor - 1

  iorbtot = 0
  do iproc = 0, size(orbs%norb_par, 1) - 1, 1
     if (iorbp >= iorbtot .and. iorbp < iorbtot + orbs%norb_par(iproc, 0)) then
        iorbp = iorbp - iorbtot
        return
     end if
     iorbtot = iorbtot + orbs%norb_par(iproc, 0)
  end do

  iorbp = -1;
  iproc = -1;
END SUBROUTINE orbs_get_iorbp

subroutine energs_new(self, energs)
  use module_types
  implicit none
  double precision, intent(in) :: self
  type(energy_terms), pointer :: energs

  allocate(energs)
  energs%c_obj = self
END SUBROUTINE energs_new
subroutine energs_free(energs)
  use module_types
  implicit none
  type(energy_terms), pointer :: energs

  deallocate(energs)
END SUBROUTINE energs_free
subroutine energs_copy_data(energs, eh, exc, evxc, eion, edisp, ekin, epot, &
     & eproj, eexctX, ebs, eKS, trH, evsum, evsic)
  use module_types
  implicit none
  type(energy_terms), intent(in) :: energs
  real(gp), intent(out) :: eh, exc, evxc, eion, edisp, ekin, epot, eproj, &
       & eexctX, ebs, eKS, trH, evsum, evsic

  eh     = energs%eh
  exc    = energs%exc
  evxc   = energs%evxc
  eion   = energs%eion
  edisp  = energs%edisp
  ekin   = energs%ekin
  epot   = energs%epot
  eproj  = energs%eproj
  eexctX = energs%eexctX
  ebs    = energs%ebs
  eKS    = energs%eKS
  trH    = energs%trH
  evsum  = energs%evsum
  evsic  = energs%evsic
END SUBROUTINE energs_copy_data
