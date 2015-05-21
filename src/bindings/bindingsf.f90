!> @file
!! Bindings for BigDFT
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!!$subroutine memocc_report()
!!$  use memory_profiling, only: mreport => memocc_report
!!$  implicit none
!!$  call mreport()
!!$end subroutine memocc_report


!!$subroutine memocc_verbose()
!!$  use memory_profiling, only: mstate => memocc_set_state
!!$  implicit none
!!$  call mstate(2)
!!$end subroutine memocc_verbose


!!$subroutine memocc_set_output(file, ln)
!!$  use memory_profiling, only: mstate => memocc_set_filename
!!$  implicit none
!!$  integer, intent(in) :: ln
!!$  character(len = ln), intent(in) :: file
!!$  call mstate(file)
!!$end subroutine memocc_set_output


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


subroutine open_write_file(unitwf, name, ln, bin)
  implicit none
  integer, intent(in) :: unitwf, ln, bin
  character(len = 1), dimension(ln), intent(in) :: name

  character(len = ln) :: filename
  integer :: i
  
  write(filename, *)
  do i = 1, ln, 1
     filename(i:i) = name(i)(1:1)
  end do
  if (bin == 1) then
     open(unit = unitwf, file = trim(filename), status = 'unknown', form = "unformatted")
  else
     open(unit = unitwf, file = trim(filename), status = 'unknown')
  end if
END SUBROUTINE open_write_file


subroutine close_file(unitwf)
  implicit none
  integer, intent(in) :: unitwf

  close(unit = unitwf)
END SUBROUTINE close_file


!subroutine createKernel(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,kernel,wrtmsg)
!  use Poisson_Solver, only: ck => createKernel
!  implicit none
!  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!  integer, intent(in) :: n01,n02,n03,itype_scf,iproc,nproc
!  real(kind=8), intent(in) :: hx,hy,hz
!  real(kind=8), pointer :: kernel(:)
!  logical, intent(in) :: wrtmsg
!
!  call ck(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,kernel,wrtmsg)
!end subroutine createKernel


subroutine deallocate_double_1D(array)
  use BigDFT_API
  implicit none

  double precision, dimension(:), pointer :: array

  if (associated(array)) then
     call f_free_ptr(array)
  end if
end subroutine deallocate_double_1D


subroutine deallocate_double_2D(array)
  use BigDFT_API
  implicit none

  double precision, dimension(:,:), pointer :: array

  if (associated(array)) then
     call f_free_ptr(array)
  end if
end subroutine deallocate_double_2D


subroutine glr_new(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  allocate(glr)
end subroutine glr_new


subroutine glr_copy(glr, d, wfd, from)
  use module_types
  use locregs, only: copy_locreg_descriptors
  implicit none
  type(locreg_descriptors), pointer :: glr
  type(grid_dimensions), pointer :: d
  type(wavefunctions_descriptors), pointer :: wfd
  type(locreg_descriptors), intent(in) :: from

  allocate(glr)
  call nullify_locreg_descriptors(glr)
  d => glr%d
  wfd => glr%wfd
  call copy_locreg_descriptors(from, glr)
end subroutine glr_copy


subroutine glr_init(glr, d, wfd)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout), target :: glr
  type(grid_dimensions), pointer :: d
  type(wavefunctions_descriptors), pointer :: wfd

  call nullify_locreg_descriptors(glr)
  d => glr%d
  wfd => glr%wfd
end subroutine glr_init


subroutine glr_get_data(glr, d, wfd)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout), target :: glr
  type(grid_dimensions), pointer :: d
  type(wavefunctions_descriptors), pointer :: wfd

  d => glr%d
  wfd => glr%wfd
end subroutine glr_get_data


subroutine glr_free(glr)
  use module_types
  implicit none
  type(locreg_descriptors), pointer :: glr

  deallocate(glr)
end subroutine glr_free


subroutine glr_empty(glr)
  use locregs
  implicit none
  type(locreg_descriptors), intent(inout) :: glr

  call deallocate_locreg_descriptors(glr)
end subroutine glr_empty


subroutine glr_get_dimensions(glr , n, ni, ns, nsi, nfl, nfu, norb)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  integer, dimension(3), intent(out) :: n, ni, ns, nsi, nfl, nfu
  integer, intent(out) :: norb

  n(1) = glr%d%n1
  n(2) = glr%d%n2
  n(3) = glr%d%n3
  ni(1) = glr%d%n1i
  ni(2) = glr%d%n2i
  ni(3) = glr%d%n3i

  nfl(1) = glr%d%nfl1
  nfl(2) = glr%d%nfl2
  nfl(3) = glr%d%nfl3
  nfu(1) = glr%d%nfu1
  nfu(2) = glr%d%nfu2
  nfu(3) = glr%d%nfu3

  ns(1) = glr%ns1
  ns(2) = glr%ns2
  ns(3) = glr%ns3
  nsi(1) = glr%nsi1
  nsi(2) = glr%nsi2
  nsi(3) = glr%nsi3
  
  norb = glr%Localnorb
end subroutine glr_get_dimensions


subroutine glr_set_dimensions(glr, n, ni, ns, nsi, nfl, nfu)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout) :: glr
  integer, dimension(3), intent(in) :: n, ni, ns, nsi, nfl, nfu

  glr%d%n1 = n(1)
  glr%d%n2 = n(2)
  glr%d%n3 = n(3)
  glr%d%n1i = ni(1)
  glr%d%n2i = ni(2)
  glr%d%n3i = ni(3)

  glr%d%nfl1 = nfl(1)
  glr%d%nfl2 = nfl(2)
  glr%d%nfl3 = nfl(3)
  glr%d%nfu1 = nfu(1)
  glr%d%nfu2 = nfu(2)
  glr%d%nfu3 = nfu(3)

  glr%ns1 = ns(1)
  glr%ns2 = ns(2)
  glr%ns3 = ns(3)
  glr%nsi1 = nsi(1)
  glr%nsi2 = nsi(2)
  glr%nsi3 = nsi(3)
end subroutine glr_set_dimensions


subroutine glr_get_locreg_data(glr, locrad, locregCenter)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  double precision, dimension(3), intent(out) :: locregCenter
  double precision, intent(out) :: locrad

  locrad = glr%locrad
  locregCenter = glr%locregCenter
end subroutine glr_get_locreg_data


subroutine glr_set_wfd_dims(glr, nseg_c, nseg_f, nvctr_c, nvctr_f)
  use module_types
  use locregs, only: allocate_wfd
  implicit none
  type(locreg_descriptors), intent(inout) :: glr
  integer, intent(in) :: nseg_c, nseg_f, nvctr_c, nvctr_f

  glr%wfd%nseg_c = nseg_c
  glr%wfd%nseg_f = nseg_f
  glr%wfd%nvctr_c = nvctr_c
  glr%wfd%nvctr_f = nvctr_f
  call allocate_wfd(glr%wfd)
END SUBROUTINE glr_set_wfd_dims


subroutine glr_set_wave_descriptors(iproc,hx,hy,hz,atoms,rxyz,&
      &   crmult,frmult,Glr)
   use module_base, only: gp
   use module_types
   use module_interfaces, only:createWavefunctionsDescriptors
   implicit none
   !Arguments
   type(atoms_data), intent(in) :: atoms
   integer, intent(in) :: iproc
   real(gp), intent(in) :: hx,hy,hz,crmult,frmult
   real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
   !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
   type(locreg_descriptors), intent(inout) :: Glr

   call createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,&
        crmult,frmult,.true.,Glr)
end subroutine glr_set_wave_descriptors


subroutine glr_set_bounds(lr)
  use module_types
  implicit none
  type(locreg_descriptors), intent(inout) :: lr
  
  call locreg_bounds(lr%d%n1,lr%d%n2,lr%d%n3, &
       & lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3, &
       & lr%wfd,lr%bounds)
END SUBROUTINE glr_set_bounds


subroutine glr_wfd_get_data(wfd, nvctr_c, nvctr_f, nseg_c, nseg_f, &
     & keyglob, keygloc, keyvglob, keyvloc)
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, intent(out) :: nvctr_c, nvctr_f, nseg_c, nseg_f
  integer, dimension(:,:), pointer :: keyglob, keygloc
  integer, dimension(:), pointer :: keyvglob, keyvloc

  nvctr_c = wfd%nvctr_c
  nvctr_f = wfd%nvctr_f
  nseg_c  = wfd%nseg_c
  nseg_f  = wfd%nseg_f
  keyglob => wfd%keyglob
  keygloc => wfd%keygloc
  keyvglob => wfd%keyvglob
  keyvloc => wfd%keyvloc
END SUBROUTINE glr_wfd_get_data


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
  lzd%nlr = 0
  glr => lzd%glr
end subroutine lzd_init


subroutine lzd_get_data(lzd, glr)
  use module_types
  implicit none
  type(local_zone_descriptors), target, intent(inout) :: lzd
  type(locreg_descriptors), pointer :: glr

  glr => lzd%glr
end subroutine lzd_get_data


subroutine lzd_copy_data(lzd, nlr)
  use module_types
  implicit none
  type(local_zone_descriptors), intent(in) :: lzd
  integer, intent(out) :: nlr

  nlr = lzd%nlr
end subroutine lzd_copy_data


subroutine lzd_free(lzd)
  use module_types
  implicit none
  type(local_zone_descriptors), pointer :: lzd

  call deallocate_local_zone_descriptors(lzd)
  deallocate(lzd)
end subroutine lzd_free


subroutine lzd_empty(lzd)
  use module_types
  implicit none
  type(local_zone_descriptors), intent(inout) :: lzd

  call deallocate_Lzd_except_Glr(lzd)
END SUBROUTINE lzd_empty


subroutine lzd_set_nlr(lzd, nlr, geocode)
  use locregs
  use module_types, only: local_zone_descriptors
  implicit none
  type(local_zone_descriptors), intent(inout) :: lzd
  integer, intent(in) :: nlr
  character, intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode

  integer :: i
  
  if (lzd%nlr > 0) then
     do i = 1, lzd%nlr, 1
        call deallocate_locreg_descriptors(lzd%Llr(i))
     end do
     deallocate(lzd%llr)
  end if

  lzd%nlr = nlr
  allocate(lzd%Llr(Lzd%nlr))
  do i = 1, nlr, 1
     call nullify_locreg_descriptors(lzd%Llr(i))
     lzd%Llr(i)%geocode = geocode
  end do
END SUBROUTINE lzd_set_nlr


subroutine lzd_get_hgrids(Lzd, hgrids)
  use module_base
  use module_types
  implicit none
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), intent(out) :: hgrids(3)
  !initial values
  hgrids = Lzd%hgrids
END SUBROUTINE lzd_get_hgrids


subroutine lzd_get_llr(Lzd, i, llr)
  use module_base
  use module_types
  implicit none
  type(local_zone_descriptors), intent(in) :: Lzd
  integer, intent(in) :: i
  type(locreg_descriptors), pointer :: llr

  llr => Lzd%Llr(i)
END SUBROUTINE lzd_get_llr


!> Allocation input_variables structure and initialize it with default values
subroutine inputs_new(in)
  use module_types
  use dictionaries
  use f_refcnts, only: nullify_f_ref
  implicit none
  type(input_variables), pointer :: in
  allocate(in)
  call nullify_f_ref(in%refcnt)
  !call default_input_variables(in)
end subroutine inputs_new


subroutine inputs_free(in)
  use module_types
  implicit none
  type(input_variables), pointer :: in

  call free_input_variables(in)
  deallocate(in)
end subroutine inputs_free


subroutine inputs_set_dict(in, level, val)

  use dictionaries
  use module_types
  use yaml_output
  implicit none
  type(input_variables), intent(inout) :: in
  character(len = *), intent(in) :: level
  type(dictionary), pointer :: val

  ! This is a patch for Intel, to be corrected properly later.
  call input_set(in, level(1:len(level)), val%child)
END SUBROUTINE inputs_set_dict


subroutine inputs_get_output(in, dir_output)

  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  character(len = 100), intent(out) :: dir_output

  dir_output = trim(in%dir_output)
END SUBROUTINE inputs_get_output


subroutine inputs_get_dft(in, hx, hy, hz, crmult, frmult, ixc, qcharge, efield, nspin, mpol, &
     & gnrm, itermax, nrepmax, ncong, idsx, dispcorr, inpsi, outpsi, outgrid, &
     & rbuf, ncongt, davidson, nvirt, nplottedvirt, sym, last_run)
  use module_defs, only: gp
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  real(gp), intent(out) :: hx, hy, hz, crmult, frmult, efield(3), gnrm, rbuf, qcharge
  integer, intent(out) :: ixc, nspin, mpol, itermax, nrepmax, ncong, idsx, &
       & dispcorr, inpsi, outpsi, outgrid, ncongt, davidson, nvirt, nplottedvirt, &
       & sym, last_run
  
  hx = in%hx
  hy = in%hy
  hz = in%hz
  crmult = in%crmult
  frmult = in%frmult
  ixc = in%ixc
  qcharge= in%qcharge !only works for integer values
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
  last_run = in%last_run
END SUBROUTINE inputs_get_dft


subroutine inputs_get_mix(in, iscf, itrpmax, norbsempty, occopt, alphamix, rpnrm_cv, &
     & gnrm_startmix, Tel, alphadiis)
  use module_defs, only: gp
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
  use module_defs, only: gp
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


subroutine inputs_get_perf(in, linear)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(out) :: linear
  
  linear = in%linear
END SUBROUTINE inputs_get_perf


!seems not used anymore as the filling of files variables has disappeared in most of the cases
!!$subroutine inputs_get_files(in, files)
!!$  use module_types
!!$  implicit none
!!$  type(input_variables), intent(in) :: in
!!$  integer, intent(out) :: files
!!$
!!$  files = in%files
!!$END SUBROUTINE inputs_get_files


subroutine inputs_get_linear(linear, inputPsiId)
  use module_types
  implicit none
  integer, intent(out) :: linear
  integer, intent(in) :: inputPsiId

  linear = 0
  if (inputPsiId == INPUT_PSI_LINEAR_AO .or. inputPsiId == INPUT_PSI_DISK_LINEAR) linear = 1
END SUBROUTINE inputs_get_linear


subroutine inputs_check_psi_id(inputpsi, input_wf_format, dir_output, ln, orbs, lorbs, iproc, nproc)
  use module_types
  use module_fragments
  use module_interfaces, only: input_check_psi_id
  implicit none
  integer, intent(out) :: input_wf_format
  integer, intent(inout) :: inputpsi
  integer, intent(in) :: iproc, ln, nproc
  character(len = ln), intent(in) :: dir_output
  type(orbitals_data), intent(in) :: orbs, lorbs
  !local variables
  type(system_fragment), dimension(:), pointer :: ref_frags  !< number of orbitals for each fragment
  character(len=100) :: frag_dir !< label for fragment subdirectories (blank if not a fragment calculation)

  nullify(ref_frags)
  frag_dir(1:len(frag_dir))=' '
  
  call input_check_psi_id(inputpsi, input_wf_format, trim(dir_output), orbs, lorbs, iproc, nproc,0, frag_dir, ref_frags)
END SUBROUTINE inputs_check_psi_id


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
!  use memory_profiling
  implicit none
  type(orbitals_data), pointer :: orbs

  deallocate(orbs)
END SUBROUTINE orbs_free


subroutine orbs_empty(orbs)
  use module_types
  implicit none
  type(orbitals_data), intent(inout) :: orbs

  call deallocate_orbitals_data(orbs)
END SUBROUTINE orbs_empty


subroutine orbs_comm_new(comms)
  use module_base
  use module_types
  use communications_base, only: comms_cubic
!  use module_interfaces
  implicit none
  type(comms_cubic), pointer :: comms

  allocate(comms)
  nullify(comms%nvctr_par)
end subroutine orbs_comm_new


subroutine orbs_comm_init(comms, orbs, lr, iproc, nproc)
  use module_base
  use module_types
  use communications_base, only: comms_cubic
  use communications_init, only: orbitals_communicators
  implicit none
  integer, intent(in) :: iproc,nproc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(comms_cubic), intent(inout) :: comms

  call orbitals_communicators(iproc,nproc,lr,orbs,comms)
end subroutine orbs_comm_init


subroutine orbs_comm_free(comms)
  use module_base
  use module_types
  use communications_base, only: comms_cubic
!  use module_interfaces
  implicit none
  type(comms_cubic), pointer :: comms

  deallocate(comms)
end subroutine orbs_comm_free


subroutine orbs_comm_empty(comms)
  use module_base
  use module_types
  use communications_base, only: comms_cubic, deallocate_comms
!  use module_interfaces
  implicit none
  type(comms_cubic), intent(inout) :: comms

  if (associated(comms%nvctr_par)) then
     call deallocate_comms(comms)
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
  use module_defs, only: wp
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(wp), dimension(:), pointer :: eval
  
  eval => orbs%eval
END SUBROUTINE orbs_get_eval


subroutine orbs_get_occup(orbs, occup)
  use module_defs, only: gp
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(gp), dimension(:), pointer :: occup
  
  occup => orbs%occup
END SUBROUTINE orbs_get_occup


subroutine orbs_get_kpts(orbs, kpts)
  use module_defs, only: gp
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(gp), dimension(:,:), pointer :: kpts
  
  kpts => orbs%kpts
END SUBROUTINE orbs_get_kpts


subroutine orbs_get_kwgts(orbs, kwgts)
  use module_defs, only: gp
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  real(gp), dimension(:), pointer :: kwgts
  
  kwgts => orbs%kwgts
END SUBROUTINE orbs_get_kwgts


subroutine orbs_get_inwhichlocreg(orbs, locreg)
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  integer, dimension(:), pointer :: locreg
  
  locreg => orbs%inwhichlocreg
END SUBROUTINE orbs_get_inwhichlocreg


subroutine orbs_get_onwhichatom(orbs, atom)
  use module_types
  implicit none
  type(orbitals_data) :: orbs
  integer, dimension(:), pointer :: atom
  
  atom => orbs%onwhichatom
END SUBROUTINE orbs_get_onwhichatom


subroutine orbs_open_file(orbs, unitwf, name, ln, iformat, iorbp, ispinor)
  use module_types
  use module_interfaces, only: open_filename_of_iorb
  implicit none
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: ln, iformat, iorbp, ispinor
  character(len = 1), dimension(ln), intent(in) :: name
  integer, intent(inout) :: unitwf

  character(len = ln) :: filename
  integer :: i, iorb_out
  
  write(filename, *)
  do i = 1, ln, 1
     filename(i:i) = name(i)(1:1)
  end do
  call open_filename_of_iorb(unitwf, (iformat == WF_FORMAT_BINARY), trim(filename), &
       & orbs, iorbp, ispinor, iorb_out)
END SUBROUTINE orbs_open_file

subroutine proj_new(nlpspd)
  use module_types
  implicit none
  type(nonlocal_psp_descriptors), pointer :: nlpspd

  allocate(nlpspd)
END SUBROUTINE proj_new


subroutine proj_free(nlpspd)
  use psp_projectors
  use module_types
  use memory_profiling
  implicit none
  type(DFT_PSP_projectors), pointer :: nlpspd

  call free_DFT_PSP_projectors(nlpspd)
END SUBROUTINE proj_free


subroutine proj_get_dimensions(nlpspd, nproj, nprojel)
  use module_types
  implicit none
  type(DFT_PSP_projectors), intent(in) :: nlpspd
  integer, intent(out) :: nproj, nprojel
  
  nproj = nlpspd%nproj
  nprojel = nlpspd%nprojel
END SUBROUTINE proj_get_dimensions


subroutine kernel_get_comm(pkernel, igroup, ngroup, iproc_grp, &
     & nproc_grp, mpi_comm)
  use module_types
  implicit none
  type(coulomb_operator), intent(inout) :: pkernel
  integer, intent(out) :: igroup, ngroup, iproc_grp, nproc_grp, mpi_comm
  igroup = pkernel%mpi_env%igroup
  ngroup = pkernel%mpi_env%ngroup
  iproc_grp = pkernel%mpi_env%iproc
  nproc_grp = pkernel%mpi_env%nproc
  mpi_comm = pkernel%mpi_env%mpi_comm
end subroutine kernel_get_comm


subroutine localfields_new(self, denspotd, rhod, dpbox)
  use module_dpbox
  use module_types
  implicit none
  integer(kind = 8), intent(in) :: self
  type(DFT_local_fields), pointer :: denspotd
  type(denspot_distribution), pointer :: dpbox
  type(rho_descriptors), pointer :: rhod

  allocate(denspotd)
  rhod => denspotd%rhod
  dpbox => denspotd%dpbox
  denspotd%c_obj = self
END SUBROUTINE localfields_new


subroutine localfields_get_data(denspotd, rhod, dpbox)
  use module_dpbox
  use module_types
  implicit none
  type(DFT_local_fields), intent(in), target :: denspotd
  type(denspot_distribution), pointer :: dpbox
  type(rho_descriptors), pointer :: rhod

  rhod => denspotd%rhod
  dpbox => denspotd%dpbox
END SUBROUTINE localfields_get_data


subroutine localfields_free(denspotd, fion, fdisp)
  use module_base
  use module_dpbox, only: dpbox_free
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use memory_profiling
  implicit none
  type(DFT_local_fields), pointer :: denspotd
  real(gp), dimension(:,:), pointer :: fion, fdisp
  
  character(len = *), parameter :: subname = "localfields_free"

  call deallocate_rho_descriptors(denspotd%rhod)
  call dpbox_free(denspotd%dpbox)
  
  call f_free_ptr(denspotd%V_ext)
  
  if (associated(denspotd%pkernelseq%kernel,target=denspotd%pkernel%kernel)) then
     nullify(denspotd%pkernelseq%kernel)
  else if (associated(denspotd%pkernelseq%kernel)) then
     call pkernel_free(denspotd%pkernelseq)
  end if
  call pkernel_free(denspotd%pkernel)

  call f_free_ptr(denspotd%rhov)
  call f_free_ptr(denspotd%V_XC)
  call f_free_ptr(denspotd%rho_C)
  
  deallocate(denspotd)

  call f_free_ptr(fion)
  call f_free_ptr(fdisp)
  
END SUBROUTINE localfields_free


subroutine localfields_copy_metadata(denspot, rhov_is, hgrid, ni, psoffset)
  use module_defs, only: gp,dp
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(out) :: rhov_is, ni(3)
  real(gp), intent(out) :: hgrid(3)
  real(dp), intent(out) :: psoffset

  rhov_is = denspot%rhov_is
  hgrid = denspot%dpbox%hgrids
  ni = denspot%dpbox%ndims
  psoffset = denspot%psoffset
END SUBROUTINE localfields_copy_metadata


subroutine localfields_get_rhov(denspot, rhov)
  use module_defs, only: dp
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: rhov

  rhov => denspot%rhov
END SUBROUTINE localfields_get_rhov


subroutine localfields_get_v_ext(denspot, v_ext)
  use module_defs, only: wp
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(wp), dimension(:,:,:,:), pointer :: v_ext

  v_ext => denspot%v_ext
END SUBROUTINE localfields_get_v_ext


subroutine localfields_get_v_xc(denspot, v_xc)
  use module_defs, only: wp
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(wp), dimension(:,:,:,:), pointer :: v_xc

  v_xc => denspot%v_xc
END SUBROUTINE localfields_get_v_xc


subroutine localfields_get_pkernel(denspot, pkernel)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in), target :: denspot
  type(coulomb_operator), pointer :: pkernel

  pkernel => denspot%pkernel
END SUBROUTINE localfields_get_pkernel


subroutine localfields_get_pkernelseq(denspot, pkernelseq)
  use module_types
  implicit none
  type(DFT_local_fields), intent(in), target :: denspot
  type(coulomb_operator), pointer :: pkernelseq

  pkernelseq => denspot%pkernelseq
END SUBROUTINE localfields_get_pkernelseq


subroutine localfields_get_rho_work(denspot, rho)
  use module_defs, only: dp
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: rho

  rho => denspot%rho_work
END SUBROUTINE localfields_get_rho_work


subroutine localfields_get_pot_work(denspot, pot)
  use module_defs, only: dp
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  real(dp), dimension(:), pointer :: pot

  pot => denspot%pot_work
END SUBROUTINE localfields_get_pot_work


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
  use communications_base, only: comms_cubic
  implicit none
  integer(kind = 8), intent(in) :: self
  type(DFT_wavefunction), pointer :: wf
  type(orbitals_data), pointer :: orbs
  type(comms_cubic), pointer :: comm
  type(local_zone_descriptors), pointer :: lzd

  allocate(wf)
  wf%c_obj = self
  call wf_init(wf)
  call orbs_init(wf%orbs)
  call nullify_local_zone_descriptors(wf%lzd)
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
  nullify(wf%comms%nvctr_par)
end subroutine wf_init


subroutine wf_get_data(wf, orbs, comm, lzd)
  use module_types
  use communications_base, only: comms_cubic
  implicit none
  type(DFT_wavefunction), target, intent(in) :: wf
  type(orbitals_data), pointer :: orbs
  type(comms_cubic), pointer :: comm
  type(local_zone_descriptors), pointer :: lzd

  orbs => wf%orbs
  comm => wf%comms
  lzd => wf%Lzd
end subroutine wf_get_data


subroutine wf_empty(wf)
  use module_base, only: f_free_ptr
  use module_types
  use memory_profiling
  implicit none
  type(DFT_wavefunction), intent(inout) :: wf

  call f_free_ptr(wf%psi)
  call f_free_ptr(wf%psit)
  call f_free_ptr(wf%hpsi)
END SUBROUTINE wf_empty


subroutine wf_free(wf)
  use module_types
  use memory_profiling
  implicit none
  type(DFT_wavefunction), pointer :: wf

  call orbs_comm_empty(wf%comms)
  call orbs_empty(wf%orbs)
  call deallocate_local_zone_descriptors(wf%lzd)
  deallocate(wf)
end subroutine wf_free


subroutine wf_get_psi(wf, psi, hpsi)
  use module_types
  implicit none
  type(DFT_wavefunction), intent(in) :: wf
  integer(kind = 8), intent(out) :: psi
  integer(kind = 8), intent(out) :: hpsi

  interface
     subroutine inquire_address1(add, pt_f)
       double precision, dimension(:), pointer :: pt_f
       integer(kind = 8), intent(out) :: add
     end subroutine inquire_address1
  end interface
  call inquire_address1(psi, wf%psi)
  call inquire_address1(hpsi, wf%hpsi)
end subroutine wf_get_psi


subroutine wf_get_psi_size(psi, psiSize)
  use module_defs, only: wp
  use module_types
  implicit none
  real(wp), dimension(:), pointer :: psi
  integer(kind = 8), intent(out) :: psiSize

  psiSize = product(shape(psi))
end subroutine wf_get_psi_size


subroutine wf_iorbp_to_psi(psir, psi, lr)
  use module_base, only: wp,f_zero
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psir
  
  character(len=*), parameter :: subname='wf_orb_to_psi'
  type(workarr_sumrho) :: w

  call initialize_work_arrays_sumrho(1,lr,.true.,w)

  !initialisation
  if (lr%geocode == 'F') then
     call f_zero(psir)
  end if

  call daub_to_isf(lr,w,psi,psir)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE wf_iorbp_to_psi


subroutine orbs_get_iorbp(orbs, iorbp, isorb, iproc, ikpt, iorb, ispin, ispinor)
  use module_types
  implicit none

  integer, intent(out) :: iorbp, isorb, iproc
  type(orbitals_data), intent(in) :: orbs
  integer, intent(in) :: ikpt, iorb, ispin, ispinor

  iorbp = (ikpt - 1) * (orbs%nspinor * orbs%norb)
  if (ispin == 1) iorbp = iorbp + (iorb - 1) * orbs%nspinor
  if (ispin == 2) iorbp = iorbp + orbs%norbu * orbs%nspinor + (iorb - 1) * orbs%nspinor
  iorbp = iorbp + ispinor - 1

  isorb = 0
  do iproc = 0, size(orbs%norb_par, 1) - 1, 1
     if (iorbp >= isorb .and. iorbp < isorb + orbs%norb_par(iproc, 0)) then
        iorbp = iorbp - isorb
        return
     end if
     isorb = isorb + orbs%norb_par(iproc, 0)
  end do

  iorbp = -1
  isorb = -1
  iproc = -1
END SUBROUTINE orbs_get_iorbp


subroutine state_properties_new(self, outs, energs, fxyz, nat)
  use module_defs, only: gp
  use module_types,only: energy_terms
  use bigdft_run
  implicit none
  integer(kind = 8), intent(in) :: self
  type(state_properties), pointer :: outs
  type(energy_terms), pointer :: energs
  real(gp), dimension(:,:), pointer :: fxyz
  integer, intent(in) :: nat

  type(state_properties), pointer :: intern

  allocate(intern)
  call init_state_properties(intern, nat)
  energs => intern%energs
  fxyz => intern%fxyz
  intern%energs%c_obj = self
  outs => intern
END SUBROUTINE state_properties_new


subroutine state_properties_free(outs)
  use module_types
  use bigdft_run
  implicit none
  type(state_properties), pointer :: outs

  call deallocate_state_properties(outs)
  deallocate(outs)
END SUBROUTINE state_properties_free


subroutine state_properties_get(outs, energs, fxyz, fdim, fnoise, pressure, strten, etot)
  use module_defs, only: gp
  use module_types, only: energy_terms
  use bigdft_run
  implicit none
  type(state_properties), intent(in), target :: outs
  type(energy_terms), pointer :: energs
  real(gp), dimension(:,:), pointer :: fxyz
  integer, intent(out) :: fdim
  real(gp), intent(out) :: fnoise, pressure
  real(gp), dimension(6), intent(out) :: strten
  real(gp), intent(out) :: etot

  energs => outs%energs
  fxyz => outs%fxyz
  fdim = outs%fdim

  fnoise = outs%fnoise
  pressure = outs%pressure
  strten = outs%strten

  etot = outs%energy
END SUBROUTINE state_properties_get


subroutine energs_copy_data(energs, eh, exc, evxc, eion, edisp, ekin, epot, &
     & eproj, eexctX, ebs, eKS, trH, evsum, evsic)
  use module_defs, only: gp
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

subroutine optloop_new(self, optloop)
  use module_types
  implicit none
  integer(kind = 8), intent(in) :: self
  type(DFT_optimization_loop), pointer :: optloop

  allocate(optloop)
  optloop%c_obj = self
END SUBROUTINE optloop_new


subroutine optloop_free(optloop)
  use module_types
  implicit none
  type(DFT_optimization_loop), pointer :: optloop

  deallocate(optloop)
END SUBROUTINE optloop_free


subroutine optloop_copy_data(optloop, gnrm_cv, rpnrm_cv, gnrm_startmix, gnrm, rpnrm, &
     &  itrpmax, nrepmax, itermax, itrp, itrep, iter, iscf, infocode)
  use module_defs, only: gp
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(in) :: optloop
  integer, intent(out) :: iscf, itrpmax, nrepmax, itermax, itrp, itrep, iter, infocode
  real(gp), intent(out) :: gnrm, rpnrm, gnrm_cv, rpnrm_cv, gnrm_startmix

  gnrm_cv = optloop%gnrm_cv 
  rpnrm_cv = optloop%rpnrm_cv 
  gnrm_startmix = optloop%gnrm_startmix 
  gnrm = optloop%gnrm 
  rpnrm = optloop%rpnrm 

  itrpmax = optloop%itrpmax 
  nrepmax = optloop%nrepmax 
  itermax = optloop%itermax 
  itrp = optloop%itrp 
  itrep = optloop%itrep 
  iter = optloop%iter 
  iscf = optloop%iscf 
  infocode = optloop%infocode
END SUBROUTINE optloop_copy_data


subroutine optloop_sync_data(optloop, gnrm_cv, rpnrm_cv, gnrm_startmix, gnrm, rpnrm, &
     &  itrpmax, nrepmax, itermax, itrp, itrep, iter, iscf, infocode)
  use module_defs, only: gp
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(inout) :: optloop
  integer, intent(in) :: iscf, itrpmax, nrepmax, itermax, itrp, itrep, iter, infocode
  real(gp), intent(in) :: gnrm, rpnrm, gnrm_cv, rpnrm_cv, gnrm_startmix

  optloop%gnrm_cv = gnrm_cv 
  optloop%rpnrm_cv = rpnrm_cv 
  optloop%gnrm_startmix = gnrm_startmix 
  optloop%gnrm = gnrm 
  optloop%rpnrm = rpnrm 

  optloop%itrpmax = itrpmax 
  optloop%nrepmax = nrepmax 
  optloop%itermax = itermax 
  optloop%itrp = itrp 
  optloop%itrep = itrep 
  optloop%iter = iter 
  optloop%iscf = iscf 
  optloop%infocode = infocode
END SUBROUTINE optloop_sync_data


subroutine optloop_emit_done(optloop, id, energs, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(inout) :: optloop
  type(energy_terms), intent(in) :: energs
  integer, intent(in) :: id, iproc, nproc

  call optloop_emit_iter(optloop, id + OPTLOOP_N_LOOPS, energs, iproc, nproc)
END SUBROUTINE optloop_emit_done


subroutine optloop_emit_iter(optloop, id, energs, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(inout) :: optloop
  type(energy_terms), intent(in) :: energs
  integer, intent(in) :: id, iproc, nproc

  integer, parameter :: SIGNAL_DONE = -1
  integer, parameter :: SIGNAL_WAIT = -2
  integer :: message

  call timing(iproc,'energs_signals','ON')
  if (iproc == 0) then
     ! Only iproc 0 emit the signal. This call is blocking.
     ! All other procs are blocked by the bcast to wait for
     ! possible transfer to proc 0.
     call optloop_emit(optloop%c_obj, id, energs%c_obj)
     if (nproc > 1) then
        ! After handling the signal, iproc 0 broadcasts to other
        ! proc to continue (jproc == -1).
        message = SIGNAL_DONE
        call mpibcast(message,1,comm= bigdft_mpi%mpi_comm)
     end if
  else
     message = SIGNAL_WAIT
     do
        if (message == SIGNAL_DONE) then
           exit
        end if
        call mpibcast(message, 1,comm=bigdft_mpi%mpi_comm)
        
        if (message >= 0) then
           ! sync values from proc 0.
           call optloop_bcast(optloop, iproc)
        end if
     end do
  end if
  call timing(iproc,'energs_signals','OF')
END SUBROUTINE optloop_emit_iter


subroutine optloop_bcast(optloop, iproc)
  use module_base
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(inout) :: optloop
  integer, intent(in) :: iproc

  integer, dimension(4) :: iData
  real(gp), dimension(3) :: rData

  if (iproc == 0) then
     iData(1) = optloop%iscf
     iData(2) = optloop%itrpmax
     iData(3) = optloop%nrepmax
     iData(4) = optloop%itermax

     rData(1) = optloop%gnrm_cv
     rData(2) = optloop%rpnrm_cv
     rData(3) = optloop%gnrm_startmix

     !what is this?
     !zero=0
     !call MPI_BCAST(zero, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)
  end if
  call mpibcast(iData,comm=bigdft_mpi%mpi_comm)
  call mpibcast(rData,comm=bigdft_mpi%mpi_comm)
  if (iproc /= 0) then
     optloop%iscf = iData(1)
     optloop%itrpmax = iData(2)
     optloop%nrepmax = iData(3)
     optloop%itermax = iData(4)

     optloop%gnrm_cv = rData(1)
     optloop%rpnrm_cv = rData(2)
     optloop%gnrm_startmix = rData(3)
  end if
END SUBROUTINE optloop_bcast


subroutine run_objects_new(runObj)
  use bigdft_run
  implicit none
  type(run_objects), pointer :: runObj

  type(run_objects), pointer :: intern

  allocate(intern)
  call nullify_run_objects(intern)
  runObj => intern

  ! Allocate persistent structures.
  allocate(runObj%rst)
  call nullify_QM_restart_objects(runObj%rst)

  ! Allocate persistent structures.
  allocate(runObj%mm_rst)
  call nullify_MM_restart_objects(runObj%mm_rst)
END SUBROUTINE run_objects_new


subroutine run_objects_destroy(runObj)
  use bigdft_run
  implicit none
  type(run_objects), pointer :: runObj

  ! The caller is responsible to nullify attributes he wants to keep.
  call free_run_objects(runObj)
  deallocate(runObj)
end subroutine run_objects_destroy


subroutine run_objects_get(runObj, dict, inputs, atoms)
  use bigdft_run
  use module_types, only: input_variables
  use module_atoms, only: atoms_data
  use dictionaries
  implicit none
  type(run_objects), intent(in) :: runObj
  type(dictionary), pointer :: dict
  type(input_variables), pointer :: inputs
  type(atoms_data), pointer :: atoms

  dict => runObj%user_inputs
  inputs => runObj%inputs
  atoms => runObj%atoms
END SUBROUTINE run_objects_get


subroutine run_objects_dump_to_file(iostat, dict, fname, userOnly,ln)
  use dictionaries, only: dictionary
  use module_input_keys, only: input_keys_dump
  use module_defs, only: UNINITIALIZED, gp
  use yaml_output
  use f_utils, only: f_get_free_unit
  use yaml_strings, only: f_strcpy
  implicit none
  integer, intent(in) :: ln
  integer, intent(out) :: iostat
  type(dictionary), pointer :: dict
  character, dimension(ln), intent(in) :: fname
  logical, intent(in) :: userOnly

  integer, parameter :: iunit_true = 145214 !< Hopefully being unique...
  integer :: iunit_def,iunit,iln
  real(gp), dimension(3), parameter :: dummy = (/ 0._gp, 0._gp, 0._gp /)
  character(len=256) :: filetmp

  !check free unit
  iunit=f_get_free_unit(iunit_true)

  call yaml_get_default_stream(iunit_def)
  if (iunit_def == iunit) then
     iostat = 1
     return
  end if
  !call f_strcpy(src=fname(1:ln),dest=filetmp)
  filetmp=''
  do iln=1,ln
     filetmp(iln:iln)=fname(iln)
  end do

  open(unit = iunit, file =trim(filetmp), iostat = iostat)
  if (iostat /= 0) return
  call yaml_set_stream(unit = iunit, tabbing = 40, record_length = 100, istat = iostat)
  if (iostat /= 0) return
  call yaml_new_document(unit = iunit)
  call input_keys_dump(dict, userOnly)
  call yaml_close_stream(iunit, iostat)
  if (iostat /= 0) return
  close(unit = iunit)
  call yaml_set_default_stream(iunit_def, iostat)
END SUBROUTINE run_objects_dump_to_file

!wrapper to bigdft_state in bigdft run
subroutine bigdft_exec(runObj,outs,infocode)
  use bigdft_run, only: run_objects,state_properties,bigdft_state
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(state_properties), intent(inout) :: outs
  integer, intent(inout) :: infocode
  call bigdft_state(runObj,outs,infocode)

end subroutine bigdft_exec

subroutine run_objects_set_dict(runObj, dict)
  use bigdft_run, only: run_objects
  use dictionaries, only: dictionary
  implicit none
  type(run_objects), intent(inout) :: runObj
  type(dictionary), pointer :: dict

  ! Warning, taking ownership here. Use run_objects_nullify_dict()
  ! to release this ownership without freeing dict.
  runObj%user_inputs => dict
END SUBROUTINE run_objects_set_dict


subroutine run_objects_nullify_dict(runObj)
  use bigdft_run, only: run_objects
  implicit none
  type(run_objects), intent(inout) :: runObj

  nullify(runObj%user_inputs)
END SUBROUTINE run_objects_nullify_dict


subroutine run_objects_nullify_volatile(runObj)
  use f_enums
  use public_enums
  use bigdft_run, only: run_objects
  use module_defs, only: bigdft_mpi
  use yaml_output, only: yaml_sequence_close
  implicit none
  type(run_objects), intent(inout) :: runObj

  if (associated(runObj%run_mode)) then
    if (bigdft_mpi%iproc == 0 .and. runObj%run_mode /= QM_RUN_MODE)&
         call yaml_sequence_close()
  end if

  nullify(runObj%inputs)
  nullify(runObj%run_mode)
  nullify(runObj%atoms)
END SUBROUTINE run_objects_nullify_volatile


subroutine mem_new(mem)
  use module_types, only: memory_estimation
  implicit none
  type(memory_estimation), pointer :: mem

  allocate(mem)
end subroutine mem_new


subroutine mem_destroy(mem)
  use module_types, only: memory_estimation
  implicit none
  type(memory_estimation), pointer :: mem

  deallocate(mem)
  nullify(mem)
end subroutine mem_destroy


subroutine mem_to_c(mem, submat, ncomponents, norb, norbp, oneorb, allpsi_mpi, &
     & psistorage, projarr, grid, workarr, kernel, density, psolver, ham, peak)
  use module_types, only: memory_estimation
  implicit none
  type(memory_estimation), intent(in) :: mem
  double precision, intent(out) :: submat, oneorb, allpsi_mpi, &
     & psistorage, projarr, grid, workarr, kernel, density, psolver, ham, peak
  integer, intent(out) :: ncomponents, norb, norbp

  submat = mem%submat
  ncomponents = mem%ncomponents
  norb = mem%norb
  norbp = mem%norbp
  oneorb = mem%oneorb
  allpsi_mpi = mem%allpsi_mpi
  psistorage = mem%psistorage
  projarr = mem%projarr
  grid = mem%grid
  workarr = mem%workarr
  kernel = mem%kernel
  density = mem%density
  psolver = mem%psolver
  ham = mem%ham
  peak = mem%peak
END SUBROUTINE mem_to_c


subroutine dict_move_to_key(dict, exists, key)
  use dictionaries, only: dictionary, operator(//), has_key
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists
  character(len = *), intent(in) :: key

  exists = has_key(dict, key(1:len(key)))
  if (exists) dict => dict // key(1:len(key))
END SUBROUTINE dict_move_to_key


subroutine dict_move_to_item(dict, exists, id)
  use dictionaries, only: dictionary, operator(//), dict_len
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists
  integer, intent(in) :: id

  exists = (id < dict_len(dict) .and. id >= 0)
  if (exists) dict => dict // id
END SUBROUTINE dict_move_to_item


subroutine dict_insert(dict, key)
  use dictionaries, only: dictionary, operator(//)
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: key

  ! This is a patch for Intel, to be corrected properly later.
  dict => dict // key(1:len(key))
END SUBROUTINE dict_insert


subroutine dict_append(dict)
  use dictionaries, only: dictionary, operator(//), dict_len
  implicit none
  type(dictionary), pointer :: dict

  dict => dict // dict_len(dict)
END SUBROUTINE dict_append


subroutine dict_put(dict, val)
  use dictionaries, only: dictionary, set
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: val

  ! This is a patch for Intel, to be corrected properly later.
  call set(dict, val(1:len(val)))
END SUBROUTINE dict_put


subroutine dict_dump(dict, unit)
  use dictionaries, only: dictionary
  use yaml_output, only: yaml_dict_dump
  implicit none
  type(dictionary), pointer :: dict
  integer, intent(in) :: unit

  if (unit < 0) then
     call yaml_dict_dump(dict)
  else
     call yaml_dict_dump(dict, unit = unit)
  end if
END SUBROUTINE dict_dump


subroutine dict_dump_to_file(dict, path)
  use dictionaries, only: dictionary
  use yaml_output, only: yaml_dict_dump, yaml_set_stream, yaml_close_stream, yaml_get_default_stream
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: path

  integer :: unit

  call yaml_set_stream(filename = trim(path), &
       & record_length = 92, setdefault = .true., tabbing = 0)
  call yaml_get_default_stream(unit = unit)
  call yaml_dict_dump(dict, unit = unit)
  call yaml_close_stream(unit = unit)
END SUBROUTINE dict_dump_to_file


subroutine dict_parse(dict, buf)
  use dictionaries, only: dictionary, operator(//), dict_len,operator(.pop.),dict_free
  use yaml_parse, only: yaml_parse_from_string
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: buf
  type(dictionary), pointer :: dict_load

  nullify(dict_load)
  call yaml_parse_from_string(dict_load, buf)
  if (dict_len(dict_load) == 1) then
     dict => dict_load .pop. 0
  end if
  call dict_free(dict_load)
END SUBROUTINE dict_parse


subroutine dict_pop(dict, exists, key)
  use dictionaries, only: dictionary, has_key, dict_remove, dict_init
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists
  character(len = *), intent(in) :: key

  exists = (has_key(dict, key(1:len(key))))
  if (exists) then
     call dict_remove(dict, key(1:len(key)))
     if (.not. associated(dict)) call dict_init(dict)
  end if
END SUBROUTINE dict_pop


subroutine dict_value_binding(dict, buf)
  use dictionaries, only: dictionary, max_field_length, wrapper => dict_value
  implicit none
  type(dictionary), pointer :: dict
  character(len = max_field_length), intent(out) :: buf
  
  buf = wrapper(dict)
END SUBROUTINE dict_value_binding


subroutine dict_key(dict, buf)
  use dictionaries, only: dictionary, max_field_length, wrapper => dict_key
  implicit none
  type(dictionary), pointer :: dict
  character(len = max_field_length), intent(out) :: buf
  
  buf = wrapper(dict)
END SUBROUTINE dict_key


subroutine dict_iter(dict, exists)
  use dictionaries, only: dictionary, wrapper => dict_iter
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists

  type(dictionary), pointer :: start

  start => wrapper(dict)
  exists = associated(start)
  if (exists) dict => start
END SUBROUTINE dict_iter


subroutine dict_next(dict, exists)
  use dictionaries, only: dictionary, wrapper => dict_next
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists

  type(dictionary), pointer :: next

  next => wrapper(dict)
  exists = associated(next)
  if (exists) dict => next
END SUBROUTINE dict_next


subroutine dict_len(dict, ln)
  use dictionaries, only: dictionary, wrapper => dict_len
  implicit none
  type(dictionary), pointer :: dict
  integer, intent(out) :: ln

  ln = wrapper(dict)
END SUBROUTINE dict_len


subroutine dict_size(dict, ln)
  use dictionaries, only: dictionary, wrapper => dict_size
  implicit none
  type(dictionary), pointer :: dict
  integer, intent(out) :: ln

  ln = wrapper(dict)
END SUBROUTINE dict_size


subroutine dict_copy_binding(dict, ref)
  use dictionaries, only: dictionary, wrapper => dict_copy
  implicit none
  type(dictionary), pointer :: dict, ref

  call wrapper(dict, ref)
END SUBROUTINE dict_copy_binding


subroutine dict_update_binding(dict, ref)
  use dictionaries, only: dictionary, wrapper => dict_update
  implicit none
  type(dictionary), pointer :: dict, ref

  call wrapper(dict, ref)
END SUBROUTINE dict_update_binding


subroutine dict_init_binding(dict)
  use dictionaries, only: dictionary, wrapper => dict_init
  implicit none
  type(dictionary), pointer :: dict

  call wrapper(dict)
END SUBROUTINE dict_init_binding


subroutine err_severe_override(callback)
  use dictionaries, only: f_err_severe_override, f_loc
  implicit none
  external :: callback
  
  write(*,*) f_loc(callback)
  call f_err_severe_override(callback)
end subroutine err_severe_override

subroutine astruct_set_from_dict_binding(astruct, dict)
  use module_input_dicts, only: astruct_set_from_dict
  use dictionaries, only: dictionary
  use module_atoms
  implicit none
  type(dictionary), pointer :: dict !< dictionary of the input variables
  !! the keys have to be declared like input_dicts module
  type(atomic_structure), intent(out) :: astruct          !< Structure created from the file

  call astruct_set_from_dict(dict, astruct)
end subroutine astruct_set_from_dict_binding
