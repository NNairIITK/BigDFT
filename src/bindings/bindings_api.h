#ifndef BINDINGS_API_H
#define BINDINGS_API_H

#undef hz

/* atoms_get_amu src/init/atoms.f90:1824 */
/* Fortran header:
subroutine atoms_get_amu(atoms, amu)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:), pointer :: amu
*/
void FC_FUNC_(atoms_get_amu, ATOMS_GET_AMU)(const _atoms_data *atoms, 
                                            f90_pointer_double *amu);
/* atoms_get_aocc src/init/atoms.f90:1834 */
/* Fortran header:
subroutine atoms_get_aocc(atoms, aocc)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: aocc
*/
void FC_FUNC_(atoms_get_aocc, ATOMS_GET_AOCC)(const _atoms_data *atoms, 
                                              f90_pointer_double_2D *aocc);
/* atoms_get_iasctype src/init/atoms.f90:1752 */
/* Fortran header:
subroutine atoms_get_iasctype(atoms, iasctype)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: iasctype
*/
void FC_FUNC_(atoms_get_iasctype, ATOMS_GET_IASCTYPE)(const _atoms_data *atoms, 
                                                      f90_pointer_int *iasctype);
/* atoms_get_iatype src/init/atoms.f90:1744 */
/* Fortran header:
subroutine atoms_get_iatype(atoms, iatype)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: iatype
*/
void FC_FUNC_(atoms_get_iatype, ATOMS_GET_IATYPE)(const _atoms_data *atoms, 
                                                  f90_pointer_int *iatype);
/* atoms_get_ifrztyp src/init/atoms.f90:1768 */
/* Fortran header:
subroutine atoms_get_ifrztyp(atoms, ifrztyp)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: ifrztyp
*/
void FC_FUNC_(atoms_get_ifrztyp, ATOMS_GET_IFRZTYP)(const _atoms_data *atoms, 
                                                    f90_pointer_int *ifrztyp);
/* atoms_get_ig_nlccpar src/init/atoms.f90:1875 */
/* Fortran header:
subroutine atoms_get_ig_nlccpar(atoms, ig_nlccpar)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: ig_nlccpar
*/
void FC_FUNC_(atoms_get_ig_nlccpar, ATOMS_GET_IG_NLCCPAR)(const _atoms_data *atoms, 
                                                          f90_pointer_double_2D *ig_nlccpar);
/* atoms_get_ixcpsp src/init/atoms.f90:1816 */
/* Fortran header:
subroutine atoms_get_ixcpsp(atoms, ixcpsp)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: ixcpsp
*/
void FC_FUNC_(atoms_get_ixcpsp, ATOMS_GET_IXCPSP)(const _atoms_data *atoms, 
                                                  f90_pointer_int *ixcpsp);
/* atoms_get_natpol src/init/atoms.f90:1760 */
/* Fortran header:
subroutine atoms_get_natpol(atoms, natpol)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: natpol
*/
void FC_FUNC_(atoms_get_natpol, ATOMS_GET_NATPOL)(const _atoms_data *atoms, 
                                                  f90_pointer_int *natpol);
/* atoms_get_nelpsp src/init/atoms.f90:1776 */
/* Fortran header:
subroutine atoms_get_nelpsp(atoms, nelpsp)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nelpsp
*/
void FC_FUNC_(atoms_get_nelpsp, ATOMS_GET_NELPSP)(const _atoms_data *atoms, 
                                                  f90_pointer_int *nelpsp);
/* atoms_get_nlcc_ngc src/init/atoms.f90:1808 */
/* Fortran header:
subroutine atoms_get_nlcc_ngc(atoms, nlcc_ngc)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nlcc_ngc
*/
void FC_FUNC_(atoms_get_nlcc_ngc, ATOMS_GET_NLCC_NGC)(const _atoms_data *atoms, 
                                                      f90_pointer_int *nlcc_ngc);
/* atoms_get_nlcc_ngv src/init/atoms.f90:1800 */
/* Fortran header:
subroutine atoms_get_nlcc_ngv(atoms, nlcc_ngv)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nlcc_ngv
*/
void FC_FUNC_(atoms_get_nlcc_ngv, ATOMS_GET_NLCC_NGV)(const _atoms_data *atoms, 
                                                      f90_pointer_int *nlcc_ngv);
/* atoms_get_nlccpar src/init/atoms.f90:1865 */
/* Fortran header:
subroutine atoms_get_nlccpar(atoms, nlccpar)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: nlccpar
*/
void FC_FUNC_(atoms_get_nlccpar, ATOMS_GET_NLCCPAR)(const _atoms_data *atoms, 
                                                    f90_pointer_double_2D *nlccpar);
/* atoms_get_npspcode src/init/atoms.f90:1784 */
/* Fortran header:
subroutine atoms_get_npspcode(atoms, npspcode)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: npspcode
*/
void FC_FUNC_(atoms_get_npspcode, ATOMS_GET_NPSPCODE)(const _atoms_data *atoms, 
                                                      f90_pointer_int *npspcode);
/* atoms_get_nzatom src/init/atoms.f90:1792 */
/* Fortran header:
subroutine atoms_get_nzatom(atoms, nzatom)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nzatom
*/
void FC_FUNC_(atoms_get_nzatom, ATOMS_GET_NZATOM)(const _atoms_data *atoms, 
                                                  f90_pointer_int *nzatom);
/* atoms_get_psppar src/init/atoms.f90:1855 */
/* Fortran header:
subroutine atoms_get_psppar(atoms, psppar)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:,:), pointer :: psppar
*/
void FC_FUNC_(atoms_get_psppar, ATOMS_GET_PSPPAR)(const _atoms_data *atoms, 
                                                  f90_pointer_double_3D *psppar);
/* atoms_get_radii_cf src/init/atoms.f90:1845 */
/* Fortran header:
subroutine atoms_get_radii_cf(atoms, radii_cf)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: radii_cf
*/
void FC_FUNC_(atoms_get_radii_cf, ATOMS_GET_RADII_CF)(const _atoms_data *atoms, 
                                                      f90_pointer_double_2D *radii_cf);
/* localfields_get_rhov src/bindings/bindingsf.f90:937 */
/* Fortran header:
subroutine localfields_get_rhov(denspot, rhov)
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
real(dp), dimension(:), pointer :: rhov
*/
void FC_FUNC_(localfields_get_rhov, LOCALFIELDS_GET_RHOV)(const _DFT_local_fields *denspot, 
                                                          f90_pointer_double *rhov);
/* localfields_get_v_ext src/bindings/bindingsf.f90:945 */
/* Fortran header:
subroutine localfields_get_v_ext(denspot, v_ext)
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
real(wp), dimension(:,:,:,:), pointer :: v_ext
*/
void FC_FUNC_(localfields_get_v_ext, LOCALFIELDS_GET_V_EXT)(const _DFT_local_fields *denspot, 
                                                            f90_pointer_double_4D *v_ext);
/* localfields_get_v_xc src/bindings/bindingsf.f90:953 */
/* Fortran header:
subroutine localfields_get_v_xc(denspot, v_xc)
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
real(wp), dimension(:,:,:,:), pointer :: v_xc
*/
void FC_FUNC_(localfields_get_v_xc, LOCALFIELDS_GET_V_XC)(const _DFT_local_fields *denspot, 
                                                          f90_pointer_double_4D *v_xc);
/* orbs_get_eval src/bindings/bindingsf.f90:725 */
/* Fortran header:
subroutine orbs_get_eval(orbs, eval)
use module_types
implicit none
type(orbitals_data) :: orbs
real(wp), dimension(:), pointer :: eval
*/
void FC_FUNC_(orbs_get_eval, ORBS_GET_EVAL)(_orbitals_data *orbs, 
                                            f90_pointer_double *eval);
/* orbs_get_inwhichlocreg src/bindings/bindingsf.f90:757 */
/* Fortran header:
subroutine orbs_get_inwhichlocreg(orbs, locreg)
use module_types
implicit none
type(orbitals_data) :: orbs
integer, dimension(:), pointer :: locreg
*/
void FC_FUNC_(orbs_get_inwhichlocreg, ORBS_GET_INWHICHLOCREG)(_orbitals_data *orbs, 
                                                              f90_pointer_int *locreg);
/* orbs_get_kpts src/bindings/bindingsf.f90:741 */
/* Fortran header:
subroutine orbs_get_kpts(orbs, kpts)
use module_types
implicit none
type(orbitals_data) :: orbs
real(gp), dimension(:,:), pointer :: kpts
*/
void FC_FUNC_(orbs_get_kpts, ORBS_GET_KPTS)(_orbitals_data *orbs, 
                                            f90_pointer_double_2D *kpts);
/* orbs_get_kwgts src/bindings/bindingsf.f90:749 */
/* Fortran header:
subroutine orbs_get_kwgts(orbs, kwgts)
use module_types
implicit none
type(orbitals_data) :: orbs
real(gp), dimension(:), pointer :: kwgts
*/
void FC_FUNC_(orbs_get_kwgts, ORBS_GET_KWGTS)(_orbitals_data *orbs, 
                                              f90_pointer_double *kwgts);
/* orbs_get_occup src/bindings/bindingsf.f90:733 */
/* Fortran header:
subroutine orbs_get_occup(orbs, occup)
use module_types
implicit none
type(orbitals_data) :: orbs
real(gp), dimension(:), pointer :: occup
*/
void FC_FUNC_(orbs_get_occup, ORBS_GET_OCCUP)(_orbitals_data *orbs, 
                                              f90_pointer_double *occup);
/* orbs_get_onwhichatom src/bindings/bindingsf.f90:773 */
/* Fortran header:
subroutine orbs_get_onwhichatom(orbs, atom)
use module_types
implicit none
type(orbitals_data) :: orbs
integer, dimension(:), pointer :: atom
*/
void FC_FUNC_(orbs_get_onwhichatom, ORBS_GET_ONWHICHATOM)(_orbitals_data *orbs, 
                                                          f90_pointer_int *atom);
/* orbs_get_onwhichmpi src/bindings/bindingsf.f90:765 */
/* Fortran header:
subroutine orbs_get_onwhichmpi(orbs, mpi)
use module_types
implicit none
type(orbitals_data) :: orbs
integer, dimension(:), pointer :: mpi
*/
void FC_FUNC_(orbs_get_onwhichmpi, ORBS_GET_ONWHICHMPI)(_orbitals_data *orbs, 
                                                        f90_pointer_int *mpi);
/* allocaterhopot src/init/denspotd.f90:424 */
/* Fortran header:
subroutine allocateRhoPot(iproc,Glr,nspin,atoms,rxyz,denspot)
use module_base
use module_types
use module_interfaces, fake_name => allocateRhoPot
use m_ab6_mixing
implicit none
integer, intent(in) :: iproc,nspin
type(locreg_descriptors), intent(in) :: Glr
type(atoms_data), intent(in) :: atoms
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot

character(len = *), parameter :: subname = "allocateRhoPot"
integer :: i_stat
*/
void FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(const int *iproc, 
                                             const _locreg_descriptors *Glr, 
                                             const int *nspin, 
                                             const _atoms_data *atoms, 
                                             const double *rxyz, 
                                             _DFT_local_fields *denspot);
/* atoms_copy_alat src/init/atoms.f90:1937 */
/* Fortran header:
subroutine atoms_copy_alat(atoms, alat1, alat2, alat3)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), intent(out) :: alat1, alat2, alat3
*/
void FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(const _atoms_data *atoms, 
                                                double *alat1, 
                                                double *alat2, 
                                                double *alat3);
/* atoms_copy_geometry_data src/init/atoms.f90:1885 */
/* Fortran header:
subroutine atoms_copy_geometry_data(atoms, geocode, format, units)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
character(len = 1), intent(out) :: geocode
character(len = 5), intent(out) :: format
character(len = 20), intent(out) :: units
*/
void FC_FUNC_(atoms_copy_geometry_data, ATOMS_COPY_GEOMETRY_DATA)(const _atoms_data *atoms, 
                                                                  char *geocode, 
                                                                  char *format, 
                                                                  char *units, 
                                                                  int str_ln_1, 
                                                                  int str_ln_2, 
                                                                  int str_ln_3);
/* atoms_copy_name src/init/atoms.f90:1911 */
/* Fortran header:
subroutine atoms_copy_name(atoms, ityp, name, ln)
use module_types
implicit none

type(atoms_data), intent(in) :: atoms
integer, intent(in) :: ityp
character(len=1), dimension(20), intent(out) :: name

integer, intent(out) :: ln

integer :: i,lname
*/
void FC_FUNC_(atoms_copy_name, ATOMS_COPY_NAME)(const _atoms_data *atoms, 
                                                const int *ityp, 
                                                char *name, 
                                                int *ln, 
                                                int str_ln_1);
/* atoms_copy_nat src/init/atoms.f90:1728 */
/* Fortran header:
subroutine atoms_copy_nat(atoms, nat)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, intent(out) :: nat
*/
void FC_FUNC_(atoms_copy_nat, ATOMS_COPY_NAT)(const _atoms_data *atoms, 
                                              int *nat);
/* atoms_copy_ntypes src/init/atoms.f90:1736 */
/* Fortran header:
subroutine atoms_copy_ntypes(atoms, ntypes)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, intent(out) :: ntypes
*/
void FC_FUNC_(atoms_copy_ntypes, ATOMS_COPY_NTYPES)(const _atoms_data *atoms, 
                                                    int *ntypes);
/* atoms_copy_psp_data src/init/atoms.f90:1899 */
/* Fortran header:
subroutine atoms_copy_psp_data(atoms, natsc, donlcc)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, intent(out) :: natsc
logical, intent(out) :: donlcc
*/
void FC_FUNC_(atoms_copy_psp_data, ATOMS_COPY_PSP_DATA)(const _atoms_data *atoms, 
                                                        int *natsc, 
                                                        int *donlcc);
/* atoms_empty src/init/atoms.f90:1666 */
/* Fortran header:
subroutine atoms_empty(atoms)
use module_types
implicit none
type(atoms_data), intent(inout) :: atoms
*/
void FC_FUNC_(atoms_empty, ATOMS_EMPTY)(_atoms_data *atoms);
/* atoms_free src/init/atoms.f90:23 */
/* Fortran header:
subroutine atoms_free(atoms)
use module_types
implicit none
type(atoms_data), pointer :: atoms
*/
void FC_FUNC_(atoms_free, ATOMS_FREE)(_atoms_data **atoms);
/* atoms_get src/init/atoms.f90:1720 */
/* Fortran header:
subroutine atoms_get(atoms, symObj)
use module_types
implicit none
type(atoms_data), intent(in), target :: atoms
type(symmetry_data), pointer :: symObj
*/
void FC_FUNC_(atoms_get, ATOMS_GET)(const _atoms_data *atoms, 
                                    _symmetry_data **symObj);
/* atoms_new src/init/atoms.f90:12 */
/* Fortran header:
subroutine atoms_new(atoms, sym)
use module_types
implicit none
type(atoms_data), pointer :: atoms
type(symmetry_data), pointer :: sym
*/
void FC_FUNC_(atoms_new, ATOMS_NEW)(_atoms_data **atoms, 
                                    _symmetry_data **sym);
/* atoms_read_variables src/init/atoms.f90:1675 */
/* Fortran header:
subroutine atoms_read_variables(atoms, nspin, occup, ln)
use module_types
use m_profiling
implicit none
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: nspin, ln
character, intent(in) :: occup(ln)

integer :: i
character(len = 1024) :: filename_
*/
void FC_FUNC_(atoms_read_variables, ATOMS_READ_VARIABLES)(_atoms_data *atoms, 
                                                          const int *nspin, 
                                                          const char *occup, 
                                                          const int *ln, 
                                                          int str_ln_1);
/* atoms_set_displacement src/init/atoms.f90:350 */
/* Fortran header:
subroutine atoms_set_displacement(atoms, rxyz, randdis)
use module_types
implicit none
type(atoms_data), intent(inout) :: atoms
real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
real(gp), intent(in) :: randdis

integer :: iat
real(gp) :: tt
*/
void FC_FUNC_(atoms_set_displacement, ATOMS_SET_DISPLACEMENT)(_atoms_data *atoms, 
                                                              double *rxyz, 
                                                              const double *randdis);
/* atoms_set_from_file src/init/atoms.f90:1624 */
/* Fortran header:
subroutine atoms_set_from_file(lstat, atoms, rxyz, filename, ln)
use module_base
use module_types
use module_interfaces
implicit none
logical, intent(out) :: lstat
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: ln
character, intent(in) :: filename(ln)
real(gp), dimension(:,:), pointer :: rxyz

integer :: status, i
character(len = 1024) :: filename_
*/
void FC_FUNC_(atoms_set_from_file, ATOMS_SET_FROM_FILE)(int *lstat, 
                                                        _atoms_data *atoms, 
                                                        f90_pointer_double_2D *rxyz, 
                                                        const char *filename, 
                                                        const int *ln, 
                                                        int str_ln_1);
/* atoms_set_name src/init/atoms.f90:1693 */
/* Fortran header:
subroutine atoms_set_name(atoms, ityp, name)
use module_types
implicit none
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: ityp
character, intent(in) :: name(20)
*/
void FC_FUNC_(atoms_set_name, ATOMS_SET_NAME)(_atoms_data *atoms, 
                                              const int *ityp, 
                                              const char *name, 
                                              int str_ln_1);
/* atoms_set_n_atoms src/init/atoms.f90:191 */
/* Fortran header:
subroutine atoms_set_n_atoms(atoms, rxyz, nat)
use module_types
use m_profiling
implicit none
type(atoms_data), intent(inout) :: atoms
real(gp), dimension(:,:), pointer :: rxyz
integer, intent(in) :: nat

character(len = *), parameter :: subname = "atoms_set_n_atoms"
integer :: i_stat
integer, parameter :: nelecmax=32
*/
void FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(_atoms_data *atoms, 
                                                    f90_pointer_double_2D *rxyz, 
                                                    const int *nat);
/* atoms_set_n_types src/init/atoms.f90:235 */
/* Fortran header:
subroutine atoms_set_n_types(atoms, ntypes)
use module_base
use module_types
implicit none
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: ntypes

character(len = *), parameter :: subname = "atoms_set_n_types"
integer :: i, i_stat
*/
void FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(_atoms_data *atoms, 
                                                    const int *ntypes);
/* atoms_set_symmetries src/init/atoms.f90:280 */
/* Fortran header:
subroutine atoms_set_symmetries(atoms, rxyz, disableSym, tol, elecfield)
use module_base
use module_types
use defs_basis
use m_ab6_symmetry
implicit none
type(atoms_data), intent(inout) :: atoms
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
logical, intent(in) :: disableSym
real(gp), intent(in) :: tol
real(gp), intent(in) :: elecfield(3)

character(len=*), parameter :: subname='atoms_set_symmetries'
integer :: i_stat, ierr, i_all
real(gp) :: rprimd(3, 3)
real(gp), dimension(:,:), allocatable :: xRed
*/
void FC_FUNC_(atoms_set_symmetries, ATOMS_SET_SYMMETRIES)(_atoms_data *atoms, 
                                                          const double *rxyz, 
                                                          const int *disableSym, 
                                                          const double *tol, 
                                                          const double *elecfield);
/* atoms_sync src/init/atoms.f90:1702 */
/* Fortran header:
subroutine atoms_sync(atoms, alat1, alat2, alat3, geocode, format, units)
use module_types
implicit none
type(atoms_data), intent(inout) :: atoms
real(gp), intent(in) :: alat1, alat2, alat3
character, intent(in) :: geocode(1)
character, intent(in) :: format(5)
character, intent(in) :: units(20)
*/
void FC_FUNC_(atoms_sync, ATOMS_SYNC)(_atoms_data *atoms, 
                                      const double *alat1, 
                                      const double *alat2, 
                                      const double *alat3, 
                                      const char *geocode, 
                                      const char *format, 
                                      const char *units, 
                                      int str_ln_1, 
                                      int str_ln_2, 
                                      int str_ln_3);
/* atoms_write src/init/atoms.f90:1647 */
/* Fortran header:
subroutine atoms_write(atoms, filename, filelen, rxyz, forces, energy, comment, ln)
use module_types
use module_interfaces, only: write_atomic_file
implicit none
integer, intent(in) :: ln, filelen
character(len = ln), intent(in) :: comment
character(len = filelen), intent(in) :: filename
type(atoms_data), intent(in) :: atoms
real(gp), intent(in) :: energy
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
real(gp), dimension(:,:), pointer :: forces
*/
void FC_FUNC_(atoms_write, ATOMS_WRITE)(const _atoms_data *atoms, 
                                        const char *filename, 
                                        const int *filelen, 
                                        const double *rxyz, 
                                        f90_pointer_double_2D *forces, 
                                        const double *energy, 
                                        const char *comment, 
                                        const int *ln, 
                                        int str_ln_1, 
                                        int str_ln_2);
/* bigdft_finalize src/external.f90:96 */
/* Fortran header:
subroutine bigdft_finalize(ierr)
use BigDFT_API
implicit none
integer, intent(out) :: ierr
*/
void FC_FUNC_(bigdft_finalize, BIGDFT_FINALIZE)(int *ierr);
/* bigdft_init_mpi_env src/external.f90:56 */
/* Fortran header:
subroutine bigdft_init_mpi_env(mpi_info,mpi_groupsize, ierr)
use BigDFT_API
implicit none

integer, dimension(4), intent(out) :: mpi_info
integer, intent(in) :: mpi_groupsize
integer, intent(out) :: ierr

integer :: iproc,nproc
*/
void FC_FUNC_(bigdft_init_mpi_env, BIGDFT_INIT_MPI_ENV)(int *mpi_info, 
                                                        const int *mpi_groupsize, 
                                                        int *ierr);
/* bigdft_init_mpi_force src/external.f90:87 */
/* Fortran header:
subroutine bigdft_init_mpi_force(igroup, ngroup)
use BigDFT_API
implicit none
integer, intent(in) :: igroup, ngroup
*/
void FC_FUNC_(bigdft_init_mpi_force, BIGDFT_INIT_MPI_FORCE)(const int *igroup, 
                                                            const int *ngroup);
/* bigdft_mpi_init  */
/* bigdft_set_input src/input_variables.f90:16 */
/* Fortran header:
subroutine bigdft_set_input(radical,posinp,rxyz,in,atoms)
use module_base
use module_types
use module_interfaces, except_this_one => bigdft_set_input
use yaml_output
implicit none


character(len=*),intent(in) :: posinp
character(len=*),intent(in) :: radical
type(input_variables), intent(inout) :: in
type(atoms_data), intent(out) :: atoms
real(gp), dimension(:,:), pointer :: rxyz 


integer :: ierr
*/
void FC_FUNC_(bigdft_set_input, BIGDFT_SET_INPUT)(const char *radical, 
                                                  const char *posinp, 
                                                  f90_pointer_double_2D *rxyz, 
                                                  _input_variables *in, 
                                                  _atoms_data *atoms, 
                                                  int str_ln_1, 
                                                  int str_ln_2);
/* call_bigdft src/cluster.f90:12 */
/* Fortran header:
subroutine call_bigdft(runObj,outs,nproc,iproc,infocode)
use module_base
use module_types
use module_interfaces, except_this_one => call_bigdft
use yaml_output
implicit none
integer, intent(in) :: iproc,nproc
type(run_objects), intent(inout) :: runObj
type(DFT_global_output), intent(inout) :: outs
integer, intent(inout) :: infocode


character(len=*), parameter :: subname='call_bigdft'
character(len=40) :: comment
logical :: exists
integer :: i_stat,i_all,ierr,inputPsiId_orig,iat,iorb,istep,i,jproc
real(gp) :: maxdiff
real(gp), dimension(:,:,:), allocatable :: rxyz_glob
*/
void FC_FUNC_(call_bigdft, CALL_BIGDFT)(_run_objects *runObj, 
                                        _DFT_global_output *outs, 
                                        const int *nproc, 
                                        const int *iproc, 
                                        int *infocode);
/* check_linear_and_create_lzd src/linear/initAndUtils.f90:332 */
/* Fortran header:
subroutine check_linear_and_create_Lzd(iproc,nproc,linType,Lzd,atoms,orbs,nspin,rxyz)
use module_base
use module_types
use module_xc
implicit none

integer, intent(in) :: iproc,nproc,nspin
type(local_zone_descriptors), intent(inout) :: Lzd
type(atoms_data), intent(in) :: atoms
type(orbitals_data),intent(inout) :: orbs
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
integer, intent(in) :: linType


character(len=*), parameter :: subname='check_linear_and_create_Lzd'
logical :: linear
integer :: iat,ityp,nspin_ig,i_all,i_stat,ilr
real(gp), dimension(:), allocatable :: locrad
logical,dimension(:),allocatable :: calculateBounds
*/
void FC_FUNC_(check_linear_and_create_lzd, CHECK_LINEAR_AND_CREATE_LZD)(const int *iproc, 
                                                                        const int *nproc, 
                                                                        const int *linType, 
                                                                        _local_zone_descriptors *Lzd, 
                                                                        const _atoms_data *atoms, 
                                                                        _orbitals_data *orbs, 
                                                                        const int *nspin, 
                                                                        const double *rxyz);
/* close_file src/bindings/bindingsf.f90:130 */
/* Fortran header:
subroutine close_file(unitwf)
implicit none
integer, intent(in) :: unitwf
*/
void FC_FUNC_(close_file, CLOSE_FILE)(const int *unitwf);
/* create_dir_output src/input_variables.f90:297 */
/* Fortran header:
subroutine create_dir_output(iproc, in)
use yaml_output
use module_types
use module_base
implicit none
integer, intent(in) :: iproc
type(input_variables), intent(inout) :: in

character(len=100) :: dirname
integer :: i_stat,ierror,ierr
*/
void FC_FUNC_(create_dir_output, CREATE_DIR_OUTPUT)(const int *iproc, 
                                                    _input_variables *in);
/* createeffectiveionicpotential src/init/ionicpot.f90:467 */
/* Fortran header:
subroutine createEffectiveIonicPotential(iproc, nproc, verb, in, atoms, rxyz, shift,  Glr, hxh, hyh, hzh, rhopotd, pkernel, pot_ion, elecfield, psoffset,rholoc)
use module_base
use module_types

implicit none

integer, intent(in) :: iproc,nproc
logical, intent(in) :: verb
real(gp), intent(in) :: hxh,hyh,hzh,psoffset
type(atoms_data), intent(in) :: atoms
type(locreg_descriptors), intent(in) :: Glr
type(input_variables), intent(in) :: in
type(denspot_distribution), intent(in) :: rhopotd
real(gp), intent(in) :: elecfield(3)
real(gp), dimension(3), intent(in) :: shift
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
type(coulomb_operator), intent(in) :: pkernel
real(wp), dimension(*), intent(inout) :: pot_ion
type(rholoc_objects),intent(in)::rholoc  


character(len = *), parameter :: subname = "createEffectiveIonicPotential"
logical :: counterions
integer :: i_stat, i_all
real(dp), dimension(:), allocatable :: counter_ions
*/
void FC_FUNC(createeffectiveionicpotential, CREATEEFFECTIVEIONICPOTENTIAL)(const int *iproc, 
                                                                           const int *nproc, 
                                                                           const int *verb, 
                                                                           const _input_variables *in, 
                                                                           const _atoms_data *atoms, 
                                                                           const double *rxyz, 
                                                                           const double *shift, 
                                                                           const _locreg_descriptors *Glr, 
                                                                           const double *hxh, 
                                                                           const double *hyh, 
                                                                           const double *hzh, 
                                                                           const _denspot_distribution *rhopotd, 
                                                                           const _coulomb_operator *pkernel, 
                                                                           double *pot_ion, 
                                                                           const double *elecfield, 
                                                                           const double *psoffset, 
                                                                           const _rholoc_objects *rholoc);
/* createprojectorsarrays src/init.f90:295 */
/* Fortran header:
subroutine createProjectorsArrays(iproc,lr,rxyz,at,orbs,   radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj_G,proj)
use module_base
use module_types
use gaussians, only: gaussian_basis
implicit none
integer, intent(in) :: iproc
real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
type(locreg_descriptors),intent(in) :: lr
type(atoms_data), intent(in) :: at
type(orbitals_data), intent(in) :: orbs
real(gp), dimension(3,at%nat), intent(in) :: rxyz
real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
type(nonlocal_psp_descriptors), intent(out) :: nlpspd
real(wp), dimension(:), pointer :: proj
type(gaussian_basis),dimension(at%ntypes),intent(in) :: proj_G

character(len=*), parameter :: subname='createProjectorsArrays'
integer :: n1,n2,n3,nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj
integer :: iat,i_stat,i_all,iseg
logical, dimension(:,:,:), allocatable :: logrid
*/
void FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)(const int *iproc, 
                                                             const _locreg_descriptors *lr, 
                                                             const double *rxyz, 
                                                             const _atoms_data *at, 
                                                             const _orbitals_data *orbs, 
                                                             const double *radii_cf, 
                                                             const double *cpmult, 
                                                             const double *fpmult, 
                                                             const double *hx, 
                                                             const double *hy, 
                                                             const double *hz, 
                                                             _nonlocal_psp_descriptors *nlpspd, 
                                                             const _gaussian_basis *proj_G, 
                                                             f90_pointer_double *proj);
/* deallocate_double_1d src/bindings/bindingsf.f90:150 */
/* Fortran header:
subroutine deallocate_double_1D(array)
use module_base
implicit none

double precision, dimension(:), pointer :: array
integer :: i_all, i_stat
*/
void FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(f90_pointer_double *array);
/* deallocate_double_2d src/bindings/bindingsf.f90:163 */
/* Fortran header:
subroutine deallocate_double_2D(array)
use module_base
implicit none

double precision, dimension(:,:), pointer :: array
integer :: i_all, i_stat
*/
void FC_FUNC_(deallocate_double_2d, DEALLOCATE_DOUBLE_2D)(f90_pointer_double_2D *array);
/* density_descriptors src/init/denspotd.f90:547 */
/* Fortran header:
subroutine density_descriptors(iproc,nproc,nspin,crmult,frmult,atoms,dpbox,rho_commun,rxyz,radii_cf,rhodsc)
use module_base
use module_types
use module_xc
use module_interfaces, except_this_one_A => density_descriptors
implicit none
integer, intent(in) :: iproc,nproc,nspin
real(gp), intent(in) :: crmult,frmult
type(atoms_data), intent(in) :: atoms
type(denspot_distribution), intent(in) :: dpbox
character(len=3), intent(in) :: rho_commun
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
type(rho_descriptors), intent(out) :: rhodsc
*/
void FC_FUNC_(density_descriptors, DENSITY_DESCRIPTORS)(const int *iproc, 
                                                        const int *nproc, 
                                                        const int *nspin, 
                                                        const double *crmult, 
                                                        const double *frmult, 
                                                        const _atoms_data *atoms, 
                                                        const _denspot_distribution *dpbox, 
                                                        const char *rho_commun, 
                                                        const double *rxyz, 
                                                        const double *radii_cf, 
                                                        _rho_descriptors *rhodsc, 
                                                        int str_ln_1);
/* denspot_communications src/init/denspotd.f90:181 */
/* Fortran header:
subroutine denspot_communications(iproc,nproc,ixc,nspin,geocode,SICapproach,dpbox)
use module_base
use module_types
use module_interfaces, except_this_one => denspot_communications
implicit none
integer, intent(in) :: ixc,nspin,iproc,nproc
character(len=1), intent(in) :: geocode
character(len=4), intent(in) :: SICapproach
type(denspot_distribution), intent(inout) :: dpbox

character(len = *), parameter :: subname = 'denspot_communications' 
integer :: i_stat
*/
void FC_FUNC_(denspot_communications, DENSPOT_COMMUNICATIONS)(const int *iproc, 
                                                              const int *nproc, 
                                                              const int *ixc, 
                                                              const int *nspin, 
                                                              const char *geocode, 
                                                              const char *SICapproach, 
                                                              _denspot_distribution *dpbox, 
                                                              int str_ln_1, 
                                                              int str_ln_2);
/* denspot_full_density src/init/denspotd.f90:238 src/init/denspotd.f90:331 */
/* Fortran header:
subroutine denspot_full_density(denspot, rho_full, iproc, new)
use module_base
use module_types
use m_profiling
implicit none
type(DFT_local_fields), intent(in) :: denspot
integer, intent(in) :: iproc
integer, intent(out) :: new
real(gp), dimension(:), pointer :: rho_full

character(len = *), parameter :: subname = "denspot_full_density"
integer :: i_stat, nslice, ierr, irhodim, irhoxcsh
*/
void FC_FUNC_(denspot_full_density, DENSPOT_FULL_DENSITY)(const _DFT_local_fields *denspot, 
                                                          f90_pointer_double *rho_full, 
                                                          const int *iproc, 
                                                          int *new);
/* denspot_full_v_ext src/init/denspotd.f90:286 src/init/denspotd.f90:383 */
/* Fortran header:
subroutine denspot_full_v_ext(denspot, pot_full, iproc, new)
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
*/
void FC_FUNC_(denspot_full_v_ext, DENSPOT_FULL_V_EXT)(const _DFT_local_fields *denspot, 
                                                      f90_pointer_double *pot_full, 
                                                      const int *iproc, 
                                                      int *new);
/* dpbox_set_box src/init/denspotd.f90:119 */
/* Fortran header:
subroutine dpbox_set_box(dpbox,Lzd)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: Lzd
type(denspot_distribution), intent(inout) :: dpbox
*/
void FC_FUNC_(dpbox_set_box, DPBOX_SET_BOX)(_denspot_distribution *dpbox, 
                                            const _local_zone_descriptors *Lzd);
/* energs_copy_data src/bindings/bindingsf.f90:1203 */
/* Fortran header:
subroutine energs_copy_data(energs, eh, exc, evxc, eion, edisp, ekin, epot,  eproj, eexctX, ebs, eKS, trH, evsum, evsic)
use module_types
implicit none
type(energy_terms), intent(in) :: energs
real(gp), intent(out) :: eh, exc, evxc, eion, edisp, ekin, epot, eproj,  eexctX, ebs, eKS, trH, evsum, evsic
*/
void FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)(const _energy_terms *energs, 
                                                  double *eh, 
                                                  double *exc, 
                                                  double *evxc, 
                                                  double *eion, 
                                                  double *edisp, 
                                                  double *ekin, 
                                                  double *epot, 
                                                  double *eproj, 
                                                  double *eexctX, 
                                                  double *ebs, 
                                                  double *eKS, 
                                                  double *trH, 
                                                  double *evsum, 
                                                  double *evsic);
/* fill_logrid src/init/gridmanipulation.f90:455 */
/* Fortran header:
subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
use module_base
implicit none

character, intent(in) :: geocode(1)
integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
real(gp), intent(in) :: rmult,hx,hy,hz
integer, dimension(nat), intent(in) :: iatype
real(gp), dimension(ntypes), intent(in) :: radii
real(gp), dimension(3,nat), intent(in) :: rxyz
logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid

real(kind=8), parameter :: eps_mach=1.d-12
integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
real(gp) :: dx,dy2,dz2,rad
*/
void FC_FUNC_(fill_logrid, FILL_LOGRID)(const char *geocode, 
                                        const int *n1, 
                                        const int *n2, 
                                        const int *n3, 
                                        const int *nl1, 
                                        const int *nu1, 
                                        const int *nl2, 
                                        const int *nu2, 
                                        const int *nl3, 
                                        const int *nu3, 
                                        const int *nbuf, 
                                        const int *nat, 
                                        const int *ntypes, 
                                        const int *iatype, 
                                        const double *rxyz, 
                                        const double *radii, 
                                        const double *rmult, 
                                        const double *hx, 
                                        const double *hy, 
                                        const double *hz, 
                                        int *logrid, 
                                        int str_ln_1);
/* free_wave_to_isf src/restart.f90:687 */
/* Fortran header:
subroutine free_wave_to_isf(psiscf)
use module_base
implicit none
real(wp), dimension(:,:,:,:), pointer :: psiscf

integer :: i_all, i_stat
*/
void FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(f90_pointer_double_4D *psiscf);
/* global_output_free src/bindings/bindingsf.f90:1174 */
/* Fortran header:
subroutine global_output_free(outs)
use module_types
implicit none
type(DFT_global_output), pointer :: outs
*/
void FC_FUNC_(global_output_free, GLOBAL_OUTPUT_FREE)(_DFT_global_output **outs);
/* global_output_get src/bindings/bindingsf.f90:1182 */
/* Fortran header:
subroutine global_output_get(outs, energs, fxyz, fdim, fnoise, pressure, strten, etot)
use module_types
implicit none
type(DFT_global_output), intent(in), target :: outs
type(energy_terms), pointer :: energs
real(gp), dimension(:,:), pointer :: fxyz
integer, intent(out) :: fdim
real(gp), intent(out) :: fnoise, pressure
real(gp), dimension(6), intent(out) :: strten
real(gp), intent(out) :: etot
*/
void FC_FUNC_(global_output_get, GLOBAL_OUTPUT_GET)(const _DFT_global_output *outs, 
                                                    _energy_terms **energs, 
                                                    f90_pointer_double_2D *fxyz, 
                                                    int *fdim, 
                                                    double *fnoise, 
                                                    double *pressure, 
                                                    double *strten, 
                                                    double *etot);
/* global_output_new src/bindings/bindingsf.f90:1159 */
/* Fortran header:
subroutine global_output_new(self, outs, energs, fxyz, nat)
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_global_output), pointer :: outs
type(energy_terms), pointer :: energs
real(gp), dimension(:,:), pointer :: fxyz
integer, intent(in) :: nat
*/
void FC_FUNC_(global_output_new, GLOBAL_OUTPUT_NEW)(const long *self, 
                                                    _DFT_global_output **outs, 
                                                    _energy_terms **energs, 
                                                    f90_pointer_double_2D *fxyz, 
                                                    const int *nat);
/* glr_copy src/bindings/bindingsf.f90:184 */
/* Fortran header:
subroutine glr_copy(glr, d, wfd, from)
use module_types
implicit none
type(locreg_descriptors), pointer :: glr
type(grid_dimensions), pointer :: d
type(wavefunctions_descriptors), pointer :: wfd
type(locreg_descriptors), intent(in) :: from
*/
void FC_FUNC_(glr_copy, GLR_COPY)(_locreg_descriptors **glr, 
                                  _grid_dimensions **d, 
                                  _wavefunctions_descriptors **wfd, 
                                  const _locreg_descriptors *from);
/* glr_empty src/bindings/bindingsf.f90:226 */
/* Fortran header:
subroutine glr_empty(glr)
use module_types
implicit none
type(locreg_descriptors), intent(inout) :: glr
*/
void FC_FUNC_(glr_empty, GLR_EMPTY)(_locreg_descriptors *glr);
/* glr_free src/bindings/bindingsf.f90:219 */
/* Fortran header:
subroutine glr_free(glr)
use module_types
implicit none
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(glr_free, GLR_FREE)(_locreg_descriptors **glr);
/* glr_get_data src/bindings/bindingsf.f90:209 */
/* Fortran header:
subroutine glr_get_data(glr, d, wfd)
use module_types
implicit none
type(locreg_descriptors), intent(inout), target :: glr
type(grid_dimensions), pointer :: d
type(wavefunctions_descriptors), pointer :: wfd
*/
void FC_FUNC_(glr_get_data, GLR_GET_DATA)(_locreg_descriptors *glr, 
                                          _grid_dimensions **d, 
                                          _wavefunctions_descriptors **wfd);
/* glr_get_dimensions src/bindings/bindingsf.f90:233 */
/* Fortran header:
subroutine glr_get_dimensions(glr , n, ni, ns, nsi, nfl, nfu, norb)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: glr
integer, dimension(3), intent(out) :: n, ni, ns, nsi, nfl, nfu
integer, intent(out) :: norb
*/
void FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(const _locreg_descriptors *glr, 
                                                      int *n, 
                                                      int *ni, 
                                                      int *ns, 
                                                      int *nsi, 
                                                      int *nfl, 
                                                      int *nfu, 
                                                      int *norb);
/* glr_get_locreg_data src/bindings/bindingsf.f90:290 */
/* Fortran header:
subroutine glr_get_locreg_data(glr, locrad, locregCenter)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: glr
double precision, dimension(3), intent(out) :: locregCenter
double precision, intent(out) :: locrad
*/
void FC_FUNC_(glr_get_locreg_data, GLR_GET_LOCREG_DATA)(const _locreg_descriptors *glr, 
                                                        double *locrad, 
                                                        double *locregCenter);
/* glr_get_psi_size src/init/kswfn.f90:1 */
/* Fortran header:
subroutine glr_get_psi_size(glr, psisize)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: glr
integer, intent(out) :: psisize
*/
void FC_FUNC_(glr_get_psi_size, GLR_GET_PSI_SIZE)(const _locreg_descriptors *glr, 
                                                  int *psisize);
/* glr_init src/bindings/bindingsf.f90:198 */
/* Fortran header:
subroutine glr_init(glr, d, wfd)
use module_types
implicit none
type(locreg_descriptors), intent(inout), target :: glr
type(grid_dimensions), pointer :: d
type(wavefunctions_descriptors), pointer :: wfd
*/
void FC_FUNC_(glr_init, GLR_INIT)(_locreg_descriptors *glr, 
                                  _grid_dimensions **d, 
                                  _wavefunctions_descriptors **wfd);
/* glr_new src/bindings/bindingsf.f90:177 */
/* Fortran header:
subroutine glr_new(glr)
use module_types
implicit none
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(glr_new, GLR_NEW)(_locreg_descriptors **glr);
/* glr_set_bounds src/bindings/bindingsf.f90:329 */
/* Fortran header:
subroutine glr_set_bounds(lr)
use module_types
implicit none
type(locreg_descriptors), intent(inout) :: lr
*/
void FC_FUNC_(glr_set_bounds, GLR_SET_BOUNDS)(_locreg_descriptors *lr);
/* glr_set_dimensions src/bindings/bindingsf.f90:263 */
/* Fortran header:
subroutine glr_set_dimensions(glr, n, ni, ns, nsi, nfl, nfu)
use module_types
implicit none
type(locreg_descriptors), intent(inout) :: glr
integer, dimension(3), intent(in) :: n, ni, ns, nsi, nfl, nfu
*/
void FC_FUNC_(glr_set_dimensions, GLR_SET_DIMENSIONS)(_locreg_descriptors *glr, 
                                                      const int *n, 
                                                      const int *ni, 
                                                      const int *ns, 
                                                      const int *nsi, 
                                                      const int *nfl, 
                                                      const int *nfu);
/* glr_set_wave_descriptors src/bindings/bindingsf.f90:312 */
/* Fortran header:
subroutine glr_set_wave_descriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,   crmult,frmult,Glr)
use module_base
use module_types
use module_interfaces
implicit none

type(atoms_data), intent(in) :: atoms
integer, intent(in) :: iproc
real(gp), intent(in) :: hx,hy,hz,crmult,frmult
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
type(locreg_descriptors), intent(inout) :: Glr
*/
void FC_FUNC_(glr_set_wave_descriptors, GLR_SET_WAVE_DESCRIPTORS)(const int *iproc, 
                                                                  const double *hx, 
                                                                  const double *hy, 
                                                                  const double *hz, 
                                                                  const _atoms_data *atoms, 
                                                                  const double *rxyz, 
                                                                  const double *radii_cf, 
                                                                  const double *crmult, 
                                                                  const double *frmult, 
                                                                  _locreg_descriptors *Glr);
/* glr_set_wfd_dims src/bindings/bindingsf.f90:300 */
/* Fortran header:
subroutine glr_set_wfd_dims(glr, nseg_c, nseg_f, nvctr_c, nvctr_f)
use module_types
implicit none
type(locreg_descriptors), intent(inout) :: glr
integer, intent(in) :: nseg_c, nseg_f, nvctr_c, nvctr_f
*/
void FC_FUNC_(glr_set_wfd_dims, GLR_SET_WFD_DIMS)(_locreg_descriptors *glr, 
                                                  const int *nseg_c, 
                                                  const int *nseg_f, 
                                                  const int *nvctr_c, 
                                                  const int *nvctr_f);
/* glr_wfd_get_data src/bindings/bindingsf.f90:338 */
/* Fortran header:
subroutine glr_wfd_get_data(wfd, nvctr_c, nvctr_f, nseg_c, nseg_f,  keyglob, keygloc, keyvglob, keyvloc)
use module_types
implicit none
type(wavefunctions_descriptors), intent(in) :: wfd
integer, intent(out) :: nvctr_c, nvctr_f, nseg_c, nseg_f
integer, dimension(:,:), pointer :: keyglob, keygloc
integer, dimension(:), pointer :: keyvglob, keyvloc
*/
void FC_FUNC_(glr_wfd_get_data, GLR_WFD_GET_DATA)(const _wavefunctions_descriptors *wfd, 
                                                  int *nvctr_c, 
                                                  int *nvctr_f, 
                                                  int *nseg_c, 
                                                  int *nseg_f, 
                                                  f90_pointer_int_2D *keyglob, 
                                                  f90_pointer_int_2D *keygloc, 
                                                  f90_pointer_int *keyvglob, 
                                                  f90_pointer_int *keyvloc);
/* gpu_free src/bindings/bindingsf.f90:1001 */
/* Fortran header:
subroutine gpu_free(GPU)
use module_types
implicit none
type(GPU_pointers), pointer :: GPU
*/
void FC_FUNC_(gpu_free, GPU_FREE)(_GPU_pointers **GPU);
/* gpu_new src/bindings/bindingsf.f90:994 */
/* Fortran header:
subroutine gpu_new(GPU)
use module_types
implicit none
type(GPU_pointers), pointer :: GPU
*/
void FC_FUNC_(gpu_new, GPU_NEW)(_GPU_pointers **GPU);
/* image_calculate src/images.f90:937 */
/* Fortran header:
subroutine image_calculate(img, iteration, id)
use yaml_output
use module_types
use module_images
use module_interfaces, only: write_atomic_file
implicit none
type(run_image), intent(inout) :: img
integer :: iteration
integer, intent(in) :: id

integer :: ierr, infocode
character(len = 4) :: fn4
*/
void FC_FUNC_(image_calculate, IMAGE_CALCULATE)(_run_image *img, 
                                                int *iteration, 
                                                const int *id);
/* image_free src/images.f90:1044 */
/* Fortran header:
subroutine image_free(img, run, outs)
use module_types
use module_images
implicit none

type(run_image), pointer :: img
type(run_objects), pointer :: run
type(DFT_global_output), pointer :: outs
*/
void FC_FUNC_(image_free, IMAGE_FREE)(_run_image **img, 
                                      _run_objects **run, 
                                      _DFT_global_output **outs);
/* image_get_attributes src/images.f90:1063 */
/* Fortran header:
subroutine image_get_attributes(img, error, F, id)
use module_images
use module_types
implicit none

type(run_image), intent(in) :: img
real(gp), intent(out) :: error, F
integer, intent(out) :: id
*/
void FC_FUNC_(image_get_attributes, IMAGE_GET_ATTRIBUTES)(const _run_image *img, 
                                                          double *error, 
                                                          double *F, 
                                                          int *id);
/* image_new src/images.f90:1018 */
/* Fortran header:
subroutine image_new(img, run, outs, atoms, inputs, rst, rxyz, algorithm)
use module_types
use module_images
implicit none

type(run_image), pointer :: img
type(run_objects), pointer :: run
type(DFT_global_output), pointer :: outs
type(input_variables), intent(in) :: inputs
type(atoms_data), intent(in) :: atoms
type(restart_objects), intent(in) :: rst
real(gp), dimension(:,:), pointer :: rxyz
integer, intent(in) :: algorithm

integer :: i_all, i_stat
*/
void FC_FUNC_(image_new, IMAGE_NEW)(_run_image **img, 
                                    _run_objects **run, 
                                    _DFT_global_output **outs, 
                                    const _atoms_data *atoms, 
                                    const _input_variables *inputs, 
                                    const _restart_objects *rst, 
                                    f90_pointer_double_2D *rxyz, 
                                    const int *algorithm);
/* images_distribute_tasks src/images.f90:967 */
/* Fortran header:
subroutine images_distribute_tasks(igroup, update, nimages, ngroup)
implicit none
integer, intent(in) :: nimages, ngroup
logical, dimension(nimages), intent(in) :: update
integer, dimension(nimages), intent(out) :: igroup

integer :: alpha, beta
integer :: i, l, m, n

n = 0
*/
void FC_FUNC_(images_distribute_tasks, IMAGES_DISTRIBUTE_TASKS)(int *igroup, 
                                                                const int *update, 
                                                                const int *nimages, 
                                                                const int *ngroup);
/* image_update_pos src/images.f90:775 */
/* Fortran header:
subroutine image_update_pos(img, iteration, posm1, posp1, Vm1, Vp1,  km1, kp1, optimization, climbing, neb)
use Minimization_routines
use module_images
implicit none
type(run_image), intent(inout) :: img
integer, intent(in) :: iteration
real(gp), intent(in) :: km1, kp1
real(gp), intent(in) :: Vm1, Vp1
real(gp), dimension(3*img%run%atoms%nat), intent(in) :: posm1, posp1
logical, intent(in) :: optimization, climbing
type(NEB_data), intent(in) :: neb

integer :: ndim
real(gp) :: Lx, Ly, Lz
real(gp), dimension(:), allocatable :: tangent, grad
*/
void FC_FUNC_(image_update_pos, IMAGE_UPDATE_POS)(_run_image *img, 
                                                  const int *iteration, 
                                                  const double *posm1, 
                                                  const double *posp1, 
                                                  const double *Vm1, 
                                                  const double *Vp1, 
                                                  const double *km1, 
                                                  const double *kp1, 
                                                  const int *optimization, 
                                                  const int *climbing, 
                                                  const _NEB_data *neb);
/* image_update_pos_from_file src/images.f90:868 */
/* Fortran header:
subroutine image_update_pos_from_file(img, iteration, filem1, filep1, km1, kp1, climbing, neb)
use Minimization_routines
use module_types
use module_images
use module_interfaces, only: read_atomic_file
implicit none
type(run_image), intent(inout) :: img
integer, intent(in) :: iteration
character(len = *), intent(in) :: filem1, filep1
real(gp), intent(in) :: km1, kp1
logical, intent(in) :: climbing
type(NEB_data), intent(in) :: neb

character(len = *), parameter :: subname = "image_update_pos_from_file"
real(gp), dimension(:,:), pointer :: rxyzm1, rxyzp1
type(atoms_data) :: atoms
real(gp) :: Vm1, Vp1
integer :: stat
*/
void FC_FUNC_(image_update_pos_from_file, IMAGE_UPDATE_POS_FROM_FILE)(_run_image *img, 
                                                                      const int *iteration, 
                                                                      const char *filem1, 
                                                                      const char *filep1, 
                                                                      const double *km1, 
                                                                      const double *kp1, 
                                                                      const int *climbing, 
                                                                      const _NEB_data *neb, 
                                                                      int str_ln_1, 
                                                                      int str_ln_2);
/* init_atomic_values src/init/sysprop.f90:379 */
/* Fortran header:
subroutine init_atomic_values(verb, atoms, ixc)
use module_base
use module_types
implicit none

integer, intent(in) :: ixc
logical, intent(in) :: verb
type(atoms_data), intent(inout) :: atoms


character(len=*), parameter :: subname='init_atomic_values'
real(gp), dimension(3) :: radii_cf
integer :: nlcc_dim, ityp, ig, j, ngv, ngc, i_stat,i_all,ierr
integer :: paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices
logical :: exists, read_radii,exist_all
character(len=27) :: filename
*/
void FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(const int *verb, 
                                                      _atoms_data *atoms, 
                                                      const int *ixc);
/* initialize_dft_local_fields src/init/denspotd.f90:11 */
/* Fortran header:
subroutine initialize_DFT_local_fields(denspot)
use module_base
use module_types
implicit none
type(DFT_local_fields), intent(inout) :: denspot
*/
void FC_FUNC_(initialize_dft_local_fields, INITIALIZE_DFT_LOCAL_FIELDS)(_DFT_local_fields *denspot);
/* init_orbitals_data_for_linear src/linear/initAndUtils.f90:627 */
/* Fortran header:
subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, at, rxyz, lorbs)
use module_base
use module_types
use module_interfaces, except_this_one => init_orbitals_data_for_linear
implicit none


integer,intent(in) :: iproc, nproc, nspinor
type(input_variables),intent(in) :: input
type(atoms_data),intent(in) :: at
real(kind=8),dimension(3,at%nat),intent(in) :: rxyz
type(orbitals_data),intent(out) :: lorbs


integer :: norb, norbu, norbd, ityp, iat, ilr, istat, iall, iorb, nlr
integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
real(kind=8),dimension(:,:),allocatable :: locregCenter
character(len=*),parameter :: subname='init_orbitals_data_for_linear'
*/
void FC_FUNC_(init_orbitals_data_for_linear, INIT_ORBITALS_DATA_FOR_LINEAR)(const int *iproc, 
                                                                            const int *nproc, 
                                                                            const int *nspinor, 
                                                                            const _input_variables *input, 
                                                                            const _atoms_data *at, 
                                                                            const double *rxyz, 
                                                                            _orbitals_data *lorbs);
/* inputs_check_psi_id src/bindings/bindingsf.f90:614 */
/* Fortran header:
subroutine inputs_check_psi_id(inputpsi, input_wf_format, dir_output, ln, orbs, lorbs, iproc, nproc)
use module_types
implicit none
integer, intent(out) :: input_wf_format
integer, intent(inout) :: inputpsi
integer, intent(in) :: iproc, ln, nproc
character(len = ln), intent(in) :: dir_output
type(orbitals_data), intent(in) :: orbs, lorbs
*/
void FC_FUNC_(inputs_check_psi_id, INPUTS_CHECK_PSI_ID)(int *inputpsi, 
                                                        int *input_wf_format, 
                                                        const char *dir_output, 
                                                        const int *ln, 
                                                        const _orbitals_data *orbs, 
                                                        const _orbitals_data *lorbs, 
                                                        const int *iproc, 
                                                        const int *nproc, 
                                                        int str_ln_1);
/* inputs_free src/bindings/bindingsf.f90:456 */
/* Fortran header:
subroutine inputs_free(in)
use module_types
implicit none
type(input_variables), pointer :: in
*/
void FC_FUNC_(inputs_free, INPUTS_FREE)(_input_variables **in);
/* inputs_get_dft src/bindings/bindingsf.f90:504 */
/* Fortran header:
subroutine inputs_get_dft(in, hx, hy, hz, crmult, frmult, ixc, chg, efield, nspin, mpol,  gnrm, itermax, nrepmax, ncong, idsx, dispcorr, inpsi, outpsi, outgrid,  rbuf, ncongt, davidson, nvirt, nplottedvirt, sym)
use module_types
implicit none
type(input_variables), intent(in) :: in
real(gp), intent(out) :: hx, hy, hz, crmult, frmult, efield(3), gnrm, rbuf
integer, intent(out) :: ixc, chg, nspin, mpol, itermax, nrepmax, ncong, idsx,  dispcorr, inpsi, outpsi, outgrid, ncongt, davidson, nvirt, nplottedvirt, sym
*/
void FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(const _input_variables *in, 
                                              double *hx, 
                                              double *hy, 
                                              double *hz, 
                                              double *crmult, 
                                              double *frmult, 
                                              int *ixc, 
                                              int *chg, 
                                              double *efield, 
                                              int *nspin, 
                                              int *mpol, 
                                              double *gnrm, 
                                              int *itermax, 
                                              int *nrepmax, 
                                              int *ncong, 
                                              int *idsx, 
                                              int *dispcorr, 
                                              int *inpsi, 
                                              int *outpsi, 
                                              int *outgrid, 
                                              double *rbuf, 
                                              int *ncongt, 
                                              int *davidson, 
                                              int *nvirt, 
                                              int *nplottedvirt, 
                                              int *sym);
/* inputs_get_files src/bindings/bindingsf.f90:597 */
/* Fortran header:
subroutine inputs_get_files(in, files)
use module_types
implicit none
type(input_variables), intent(in) :: in
integer, intent(out) :: files
*/
void FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(const _input_variables *in, 
                                                  int *files);
/* inputs_get_geopt src/bindings/bindingsf.f90:563 */
/* Fortran header:
subroutine inputs_get_geopt(in, geopt_approach, ncount_cluster_x, frac_fluct, forcemax,  randdis, betax, history, ionmov, dtion, strtarget, qmass)
use module_types
implicit none
type(input_variables), intent(in) :: in
character(len = 10), intent(out) :: geopt_approach
integer, intent(out) :: ncount_cluster_x, history, ionmov
real(gp), intent(out) :: frac_fluct, forcemax, randdis, betax, dtion, strtarget(6)
real(gp), pointer :: qmass(:)
*/
void FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(const _input_variables *in, 
                                                  char *geopt_approach, 
                                                  int *ncount_cluster_x, 
                                                  double *frac_fluct, 
                                                  double *forcemax, 
                                                  double *randdis, 
                                                  double *betax, 
                                                  int *history, 
                                                  int *ionmov, 
                                                  double *dtion, 
                                                  double *strtarget, 
                                                  f90_pointer_double *qmass, 
                                                  int str_ln_1);
/* inputs_get_linear src/bindings/bindingsf.f90:605 */
/* Fortran header:
subroutine inputs_get_linear(linear, inputPsiId)
use module_types
implicit none
integer, intent(out) :: linear
integer, intent(in) :: inputPsiId
*/
void FC_FUNC_(inputs_get_linear, INPUTS_GET_LINEAR)(int *linear, 
                                                    const int *inputPsiId);
/* inputs_get_mix src/bindings/bindingsf.f90:544 */
/* Fortran header:
subroutine inputs_get_mix(in, iscf, itrpmax, norbsempty, occopt, alphamix, rpnrm_cv,  gnrm_startmix, Tel, alphadiis)
use module_types
implicit none
type(input_variables), intent(in) :: in
integer, intent(out) :: iscf, itrpmax, norbsempty, occopt
real(gp), intent(out) :: alphamix, rpnrm_cv, gnrm_startmix, Tel, alphadiis
*/
void FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(const _input_variables *in, 
                                              int *iscf, 
                                              int *itrpmax, 
                                              int *norbsempty, 
                                              int *occopt, 
                                              double *alphamix, 
                                              double *rpnrm_cv, 
                                              double *gnrm_startmix, 
                                              double *Tel, 
                                              double *alphadiis);
/* inputs_get_naming src/bindings/bindingsf.f90:474 */
/* Fortran header:
subroutine inputs_get_naming(in, run_name, file_dft, file_geopt, file_kpt,  file_perf, file_tddft, file_mix, file_sic, file_occnum, file_igpop, file_lin)
use module_types
implicit none
type(input_variables), intent(in) :: in
character(len = 100), intent(out) :: run_name, file_dft, file_geopt, file_kpt,  file_perf, file_tddft, file_mix, file_sic, file_occnum, file_igpop, file_lin
*/
void FC_FUNC_(inputs_get_naming, INPUTS_GET_NAMING)(const _input_variables *in, 
                                                    char *run_name, 
                                                    char *file_dft, 
                                                    char *file_geopt, 
                                                    char *file_kpt, 
                                                    char *file_perf, 
                                                    char *file_tddft, 
                                                    char *file_mix, 
                                                    char *file_sic, 
                                                    char *file_occnum, 
                                                    char *file_igpop, 
                                                    char *file_lin, 
                                                    int str_ln_1, 
                                                    int str_ln_2, 
                                                    int str_ln_3, 
                                                    int str_ln_4, 
                                                    int str_ln_5, 
                                                    int str_ln_6, 
                                                    int str_ln_7, 
                                                    int str_ln_8, 
                                                    int str_ln_9, 
                                                    int str_ln_10, 
                                                    int str_ln_11);
/* inputs_get_output src/bindings/bindingsf.f90:494 */
/* Fortran header:
subroutine inputs_get_output(in, dir_output, writing_directory)
use module_types
implicit none
type(input_variables), intent(in) :: in
character(len = 100), intent(out) :: dir_output
character(len = 500), intent(out) :: writing_directory
*/
void FC_FUNC_(inputs_get_output, INPUTS_GET_OUTPUT)(const _input_variables *in, 
                                                    char *dir_output, 
                                                    char *writing_directory, 
                                                    int str_ln_1, 
                                                    int str_ln_2);
/* inputs_get_perf src/bindings/bindingsf.f90:589 */
/* Fortran header:
subroutine inputs_get_perf(in, linear)
use module_types
implicit none
type(input_variables), intent(in) :: in
integer, intent(out) :: linear
*/
void FC_FUNC_(inputs_get_perf, INPUTS_GET_PERF)(const _input_variables *in, 
                                                int *linear);
/* inputs_new src/bindings/bindingsf.f90:448 */
/* Fortran header:
subroutine inputs_new(in)
use module_types
implicit none
type(input_variables), pointer :: in
*/
void FC_FUNC_(inputs_new, INPUTS_NEW)(_input_variables **in);
/* inputs_parse_add src/init/wavefunctions.f90:741 */
/* Fortran header:
subroutine inputs_parse_add(in, atoms, iproc, dump)
use module_types
use yaml_output
implicit none
type(input_variables), intent(inout) :: in
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: iproc
logical, intent(in) :: dump

integer :: ierr
*/
void FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(_input_variables *in, 
                                                  _atoms_data *atoms, 
                                                  const int *iproc, 
                                                  const int *dump);
/* inputs_parse_params src/init/wavefunctions.f90:697 */
/* Fortran header:
subroutine inputs_parse_params(in, iproc, dump)
use module_types
use module_defs
use module_xc
use yaml_output
implicit none
type(input_variables), intent(inout) :: in
integer, intent(in) :: iproc
logical, intent(in) :: dump

integer :: ierr
*/
void FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(_input_variables *in, 
                                                        const int *iproc, 
                                                        const int *dump);
/* inputs_set_restart src/bindings/bindingsf.f90:625 */
/* Fortran header:
subroutine inputs_set_restart(in, id)
use module_types
implicit none
type(input_variables), intent(inout) :: in
integer, intent(in) :: id
*/
void FC_FUNC_(inputs_set_restart, INPUTS_SET_RESTART)(_input_variables *in, 
                                                      const int *id);
/* input_wf src/init.f90:2346 */
/* Fortran header:
subroutine input_wf(iproc,nproc,in,GPU,atoms,rxyz,denspot,denspot0,nlpspd,proj,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,wfd_old,psi_old,d_old,hx_old,hy_old,hz_old,rxyz_old,tmb_old)
use module_defs
use module_types
use module_interfaces, except_this_one => input_wf
use yaml_output
use gaussians, only:gaussian_basis
implicit none

integer, intent(in) :: iproc, nproc, inputpsi, input_wf_format
type(input_variables), intent(in) :: in
type(GPU_pointers), intent(inout) :: GPU
real(gp), intent(in) :: hx_old,hy_old,hz_old
type(atoms_data), intent(inout) :: atoms
real(gp), dimension(3, atoms%nat), target, intent(in) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(DFT_wavefunction), intent(inout) :: KSwfn,tmb,tmb_old 
real(gp), dimension(*), intent(out) :: denspot0 
type(energy_terms), intent(inout) :: energs 

real(wp), dimension(:), pointer :: psi_old
integer, intent(out) :: norbv
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
real(kind=8), dimension(:), pointer :: proj


type(grid_dimensions), intent(in) :: d_old
real(gp), dimension(3, atoms%nat), intent(inout) :: rxyz_old
type(wavefunctions_descriptors), intent(inout) :: wfd_old

character(len = *), parameter :: subname = "input_wf"
integer :: i_stat, nspin, i_all
type(gaussian_basis) :: Gvirt
real(wp), allocatable, dimension(:) :: norm

integer :: iatyp
type(gaussian_basis),dimension(atoms%ntypes)::proj_G
type(paw_objects)::paw
logical :: overlap_calculated
*/
void FC_FUNC_(input_wf, INPUT_WF)(const int *iproc, 
                                  const int *nproc, 
                                  const _input_variables *in, 
                                  _GPU_pointers *GPU, 
                                  _atoms_data *atoms, 
                                  const double *rxyz, 
                                  _DFT_local_fields *denspot, 
                                  double *denspot0, 
                                  const _nonlocal_psp_descriptors *nlpspd, 
                                  f90_pointer_double *proj, 
                                  _DFT_wavefunction *KSwfn, 
                                  _DFT_wavefunction *tmb, 
                                  _energy_terms *energs, 
                                  const int *inputpsi, 
                                  const int *input_wf_format, 
                                  int *norbv, 
                                  _wavefunctions_descriptors *wfd_old, 
                                  f90_pointer_double *psi_old, 
                                  const _grid_dimensions *d_old, 
                                  const double *hx_old, 
                                  const double *hy_old, 
                                  const double *hz_old, 
                                  double *rxyz_old, 
                                  _DFT_wavefunction *tmb_old);
/* ionicenergyandforces src/init/ionicpot.f90:12 */
/* Fortran header:
subroutine IonicEnergyandForces(iproc,nproc,dpbox,at,elecfield, rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr,n1,n2,n3, pot_ion,pkernel,psoffset)
use module_base
use module_types
use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
use vdwcorrection
use yaml_output
implicit none
type(denspot_distribution), intent(in) :: dpbox
type(atoms_data), intent(in) :: at
integer, intent(in) :: iproc,nproc,n1,n2,n3,dispersion
real(gp), dimension(3), intent(in) :: elecfield
real(gp), dimension(3,at%nat), intent(in) :: rxyz
type(coulomb_operator), intent(in) :: pkernel
real(gp), intent(out) :: eion,edisp,psoffset
real(dp), dimension(6),intent(out) :: ewaldstr
real(gp), dimension(:,:), pointer :: fion,fdisp
real(dp), dimension(*), intent(out) :: pot_ion

character(len=*), parameter :: subname='IonicEnergyandForces'
logical :: slowion=.false.
logical :: perx,pery,perz,gox,goy,goz
integer :: n1i,n2i,n3i,i3s,n3pi
integer :: i,iat,ii,i_all,i_stat,ityp,jat,jtyp,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
integer :: isx,iex,isy,iey,isz,iez,i1,i2,i3,j1,j2,j3,ind,ierr
real(gp) :: ucvol,rloc,twopitothreehalf,pi,atint,shortlength,charge,eself,rx,ry,rz
real(gp) :: fxion,fyion,fzion,dist,fxerf,fyerf,fzerf,cutoff
real(gp) :: hxh,hyh,hzh
real(gp) :: hxx,hxy,hxz,hyy,hyz,hzz,chgprod
real(gp) :: x,y,z,xp,Vel,prefactor,r2,arg,ehart,de

real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd

real(gp), dimension(:,:), allocatable :: fewald,xred
real(gp), dimension(3) :: cc
*/
void FC_FUNC(ionicenergyandforces, IONICENERGYANDFORCES)(const int *iproc, 
                                                         const int *nproc, 
                                                         const _denspot_distribution *dpbox, 
                                                         const _atoms_data *at, 
                                                         const double *elecfield, 
                                                         const double *rxyz, 
                                                         double *eion, 
                                                         f90_pointer_double_2D *fion, 
                                                         const int *dispersion, 
                                                         double *edisp, 
                                                         f90_pointer_double_2D *fdisp, 
                                                         double *ewaldstr, 
                                                         const int *n1, 
                                                         const int *n2, 
                                                         const int *n3, 
                                                         double *pot_ion, 
                                                         const _coulomb_operator *pkernel, 
                                                         double *psoffset);
/* kernel_get_comm src/bindings/bindingsf.f90:831 */
/* Fortran header:
subroutine kernel_get_comm(pkernel, igroup, ngroup, iproc_grp,  nproc_grp, mpi_comm)
use module_types
implicit none
type(coulomb_operator), intent(in) :: pkernel
integer, intent(out) :: igroup, ngroup, iproc_grp, nproc_grp, mpi_comm
*/
void FC_FUNC_(kernel_get_comm, KERNEL_GET_COMM)(const _coulomb_operator *pkernel, 
                                                int *igroup, 
                                                int *ngroup, 
                                                int *iproc_grp, 
                                                int *nproc_grp, 
                                                int *mpi_comm);
/* kswfn_init_comm src/init/kswfn.f90:105 */
/* Fortran header:
subroutine kswfn_init_comm(wfn, in, atoms, dpbox, iproc, nproc)
use module_types
use module_interfaces, except_this_one => kswfn_init_comm
implicit none
integer, intent(in) :: iproc, nproc
type(DFT_wavefunction), intent(inout) :: wfn
type(input_variables), intent(in) :: in
type(atoms_data),intent(in) :: atoms
type(denspot_distribution), intent(in) :: dpbox
*/
void FC_FUNC_(kswfn_init_comm, KSWFN_INIT_COMM)(_DFT_wavefunction *wfn, 
                                                const _input_variables *in, 
                                                const _atoms_data *atoms, 
                                                const _denspot_distribution *dpbox, 
                                                const int *iproc, 
                                                const int *nproc);
/* kswfn_mpi_copy src/init/kswfn.f90:87 */
/* Fortran header:
subroutine kswfn_mpi_copy(psic, jproc, psiStart, psiSize)
use module_base
use module_types
implicit none
integer, intent(in) :: psiSize, jproc, psiStart
real(wp), intent(inout) :: psic(psiSize)

integer :: ierr
integer :: status(MPI_STATUS_SIZE)
*/
void FC_FUNC_(kswfn_mpi_copy, KSWFN_MPI_COPY)(double *psic, 
                                              const int *jproc, 
                                              const int *psiStart, 
                                              const int *psiSize);
/* kswfn_optimization_loop src/cluster.f90:1375 */
/* Fortran header:
subroutine kswfn_optimization_loop(iproc, nproc, opt,  alphamix, idsx, inputpsi, KSwfn, denspot, nlpspd, proj, energs, atoms, rxyz, GPU, xcstr,  in)
use module_base
use module_types
use module_interfaces, except_this_one => kswfn_optimization_loop
use yaml_output
use m_ab6_mixing
implicit none
real(dp), dimension(6), intent(out) :: xcstr
integer, intent(in) :: iproc, nproc, idsx, inputpsi
real(gp), intent(in) :: alphamix
type(DFT_optimization_loop), intent(inout) :: opt
type(DFT_wavefunction), intent(inout) :: KSwfn
type(DFT_local_fields), intent(inout) :: denspot
type(energy_terms), intent(inout) :: energs
type(atoms_data), intent(in) :: atoms
type(GPU_pointers), intent(inout) :: GPU
type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
real(kind=8), dimension(:), pointer :: proj
real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
type(input_variables), intent(in) :: in 

character(len = *), parameter :: subname = "kswfn_optimization_loop"
logical :: endloop, scpot, endlooprp, lcs
integer :: ndiis_sd_sw, idsx_actual_before, linflag, ierr,iter_for_diis
real(gp) :: gnrm_zero
character(len=5) :: final_out
*/
void FC_FUNC_(kswfn_optimization_loop, KSWFN_OPTIMIZATION_LOOP)(const int *iproc, 
                                                                const int *nproc, 
                                                                _DFT_optimization_loop *opt, 
                                                                const double *alphamix, 
                                                                const int *idsx, 
                                                                const int *inputpsi, 
                                                                _DFT_wavefunction *KSwfn, 
                                                                _DFT_local_fields *denspot, 
                                                                _nonlocal_psp_descriptors *nlpspd, 
                                                                f90_pointer_double *proj, 
                                                                _energy_terms *energs, 
                                                                const _atoms_data *atoms, 
                                                                const double *rxyz, 
                                                                _GPU_pointers *GPU, 
                                                                double *xcstr, 
                                                                const _input_variables *in);
/* kswfn_post_treatments src/cluster.f90:1759 */
/* Fortran header:
subroutine kswfn_post_treatments(iproc, nproc, KSwfn, tmb, linear,  fxyz, fnoise, fion, fdisp, fpulay,  strten, pressure, ewaldstr, xcstr,  GPU, energs, denspot, atoms, rxyz, nlpspd, proj,  output_denspot, dir_output, gridformat, refill_proj, calculate_dipole)
use module_base
use module_types
use module_interfaces, except_this_one => kswfn_post_treatments
use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
use yaml_output

implicit none

type(DFT_wavefunction), intent(in) :: KSwfn
type(DFT_wavefunction), intent(inout) :: tmb
type(GPU_pointers), intent(inout) :: GPU
type(energy_terms), intent(in) :: energs
type(DFT_local_fields), intent(inout) :: denspot
type(atoms_data), intent(in) :: atoms
type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
real(kind=8), dimension(:), pointer :: proj
logical, intent(in) :: linear, refill_proj, calculate_dipole
integer, intent(in) :: output_denspot, iproc, nproc
character(len = *), intent(in) :: dir_output
character(len = *), intent(in) :: gridformat
real(gp), dimension(3, atoms%nat), intent(in) :: rxyz
real(gp), dimension(3, atoms%nat), intent(in) :: fdisp, fion, fpulay
real(dp), dimension(6), intent(in) :: ewaldstr, xcstr
real(gp), intent(out) :: fnoise, pressure
real(gp), dimension(6), intent(out) :: strten
real(gp), dimension(3, atoms%nat), intent(out) :: fxyz

character(len = *), parameter :: subname = "kswfn_post_treatments"
integer :: i_all, i_stat, jproc, nsize_psi, imode
real(dp), dimension(6) :: hstrten
real(gp) :: ehart_fake
*/
void FC_FUNC_(kswfn_post_treatments, KSWFN_POST_TREATMENTS)(const int *iproc, 
                                                            const int *nproc, 
                                                            const _DFT_wavefunction *KSwfn, 
                                                            _DFT_wavefunction *tmb, 
                                                            const int *linear, 
                                                            double *fxyz, 
                                                            double *fnoise, 
                                                            const double *fion, 
                                                            const double *fdisp, 
                                                            const double *fpulay, 
                                                            double *strten, 
                                                            double *pressure, 
                                                            const double *ewaldstr, 
                                                            const double *xcstr, 
                                                            _GPU_pointers *GPU, 
                                                            const _energy_terms *energs, 
                                                            _DFT_local_fields *denspot, 
                                                            const _atoms_data *atoms, 
                                                            const double *rxyz, 
                                                            _nonlocal_psp_descriptors *nlpspd, 
                                                            f90_pointer_double *proj, 
                                                            const int *output_denspot, 
                                                            const char *dir_output, 
                                                            const char *gridformat, 
                                                            const int *refill_proj, 
                                                            const int *calculate_dipole, 
                                                            int str_ln_1, 
                                                            int str_ln_2);
/* localfields_copy_metadata src/bindings/bindingsf.f90:924 */
/* Fortran header:
subroutine localfields_copy_metadata(denspot, rhov_is, hgrid, ni, psoffset)
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
integer, intent(out) :: rhov_is, ni(3)
real(gp), intent(out) :: hgrid(3)
real(dp), intent(out) :: psoffset
*/
void FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)(const _DFT_local_fields *denspot, 
                                                                    int *rhov_is, 
                                                                    double *hgrid, 
                                                                    int *ni, 
                                                                    double *psoffset);
/* localfields_free src/bindings/bindingsf.f90:866 */
/* Fortran header:
subroutine localfields_free(denspotd, fion, fdisp)
use module_types
use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
use m_profiling
implicit none
type(DFT_local_fields), pointer :: denspotd
real(gp), dimension(:,:), pointer :: fion, fdisp

character(len = *), parameter :: subname = "localfields_free"
integer :: i_stat, i_all
*/
void FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(_DFT_local_fields **denspotd, 
                                                  f90_pointer_double_2D *fion, 
                                                  f90_pointer_double_2D *fdisp);
/* localfields_get_data src/bindings/bindingsf.f90:856 */
/* Fortran header:
subroutine localfields_get_data(denspotd, rhod, dpbox)
use module_types
implicit none
type(DFT_local_fields), intent(in), target :: denspotd
type(denspot_distribution), pointer :: dpbox
type(rho_descriptors), pointer :: rhod
*/
void FC_FUNC_(localfields_get_data, LOCALFIELDS_GET_DATA)(const _DFT_local_fields *denspotd, 
                                                          _rho_descriptors **rhod, 
                                                          _denspot_distribution **dpbox);
/* localfields_get_pkernel src/bindings/bindingsf.f90:961 */
/* Fortran header:
subroutine localfields_get_pkernel(denspot, pkernel)
use module_types
implicit none
type(DFT_local_fields), intent(in), target :: denspot
type(coulomb_operator), pointer :: pkernel
*/
void FC_FUNC_(localfields_get_pkernel, LOCALFIELDS_GET_PKERNEL)(const _DFT_local_fields *denspot, 
                                                                _coulomb_operator **pkernel);
/* localfields_get_pkernelseq src/bindings/bindingsf.f90:969 */
/* Fortran header:
subroutine localfields_get_pkernelseq(denspot, pkernelseq)
use module_types
implicit none
type(DFT_local_fields), intent(in), target :: denspot
type(coulomb_operator), pointer :: pkernelseq
*/
void FC_FUNC_(localfields_get_pkernelseq, LOCALFIELDS_GET_PKERNELSEQ)(const _DFT_local_fields *denspot, 
                                                                      _coulomb_operator **pkernelseq);
/* localfields_new src/bindings/bindingsf.f90:843 */
/* Fortran header:
subroutine localfields_new(self, denspotd, rhod, dpbox)
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_local_fields), pointer :: denspotd
type(denspot_distribution), pointer :: dpbox
type(rho_descriptors), pointer :: rhod
*/
void FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(const long *self, 
                                                _DFT_local_fields **denspotd, 
                                                _rho_descriptors **rhod, 
                                                _denspot_distribution **dpbox);
/* lzd_copy_data src/bindings/bindingsf.f90:382 */
/* Fortran header:
subroutine lzd_copy_data(lzd, nlr)
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: lzd
integer, intent(out) :: nlr
*/
void FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(const _local_zone_descriptors *lzd, 
                                            int *nlr);
/* lzd_empty src/bindings/bindingsf.f90:398 */
/* Fortran header:
subroutine lzd_empty(lzd)
use module_types
implicit none
type(local_zone_descriptors), intent(inout) :: lzd
*/
void FC_FUNC_(lzd_empty, LZD_EMPTY)(_local_zone_descriptors *lzd);
/* lzd_free src/bindings/bindingsf.f90:390 */
/* Fortran header:
subroutine lzd_free(lzd)
use module_types
implicit none
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(lzd_free, LZD_FREE)(_local_zone_descriptors **lzd);
/* lzd_get_data src/bindings/bindingsf.f90:374 */
/* Fortran header:
subroutine lzd_get_data(lzd, glr)
use module_types
implicit none
type(local_zone_descriptors), target, intent(inout) :: lzd
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(lzd_get_data, LZD_GET_DATA)(_local_zone_descriptors *lzd, 
                                          _locreg_descriptors **glr);
/* lzd_get_hgrids src/bindings/bindingsf.f90:428 */
/* Fortran header:
subroutine lzd_get_hgrids(Lzd, hgrids)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: Lzd
real(gp), intent(out) :: hgrids(3)
*/
void FC_FUNC_(lzd_get_hgrids, LZD_GET_HGRIDS)(const _local_zone_descriptors *Lzd, 
                                              double *hgrids);
/* lzd_get_llr src/bindings/bindingsf.f90:437 */
/* Fortran header:
subroutine lzd_get_llr(Lzd, i, llr)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: Lzd
integer, intent(in) :: i
type(locreg_descriptors), pointer :: llr
*/
void FC_FUNC_(lzd_get_llr, LZD_GET_LLR)(const _local_zone_descriptors *Lzd, 
                                        const int *i, 
                                        _locreg_descriptors **llr);
/* lzd_init src/bindings/bindingsf.f90:364 */
/* Fortran header:
subroutine lzd_init(lzd, glr)
use module_types
implicit none
type(local_zone_descriptors), target, intent(inout) :: lzd
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(lzd_init, LZD_INIT)(_local_zone_descriptors *lzd, 
                                  _locreg_descriptors **glr);
/* lzd_init_llr src/linear/initAndUtils.f90:728 */
/* Fortran header:
subroutine lzd_init_llr(iproc, nproc, input, at, rxyz, orbs, lzd)
use module_base
use module_types
use module_interfaces
implicit none


integer,intent(in) :: iproc, nproc
type(input_variables),intent(in) :: input
type(atoms_data),intent(in) :: at
real(kind=8),dimension(3,at%nat),intent(in) :: rxyz
type(orbitals_data),intent(in) :: orbs
type(local_zone_descriptors),intent(inout) :: lzd


integer :: iat, ityp, ilr, istat, iorb, iall
real(kind=8),dimension(:,:),allocatable :: locregCenter
character(len=*),parameter :: subname='lzd_init_llr'
real(8):: t1, t2
*/
void FC_FUNC_(lzd_init_llr, LZD_INIT_LLR)(const int *iproc, 
                                          const int *nproc, 
                                          const _input_variables *input, 
                                          const _atoms_data *at, 
                                          const double *rxyz, 
                                          const _orbitals_data *orbs, 
                                          _local_zone_descriptors *lzd);
/* lzd_new src/bindings/bindingsf.f90:357 */
/* Fortran header:
subroutine lzd_new(lzd)
use module_types
implicit none
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(lzd_new, LZD_NEW)(_local_zone_descriptors **lzd);
/* lzd_set_hgrids src/init/wavefunctions.f90:686 */
/* Fortran header:
subroutine lzd_set_hgrids(Lzd, hgrids)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(inout) :: Lzd
real(gp), intent(in) :: hgrids(3)
*/
void FC_FUNC_(lzd_set_hgrids, LZD_SET_HGRIDS)(_local_zone_descriptors *Lzd, 
                                              const double *hgrids);
/* lzd_set_nlr src/bindings/bindingsf.f90:405 */
/* Fortran header:
subroutine lzd_set_nlr(lzd, nlr, geocode)
use module_types
implicit none
type(local_zone_descriptors), intent(inout) :: lzd
integer, intent(in) :: nlr
character, intent(in) :: geocode

integer :: i
*/
void FC_FUNC_(lzd_set_nlr, LZD_SET_NLR)(_local_zone_descriptors *lzd, 
                                        const int *nlr, 
                                        const char *geocode, 
                                        int str_ln_1);
/* memoryestimator src/profiling/memoryestimator.f90:12 */
/* Fortran header:
subroutine MemoryEstimator(nproc,idsx,lr,nat,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,peakmem)

use module_base
use module_types
use Poisson_Solver
use yaml_output

implicit none


integer, intent(in) :: nproc,idsx,nat,norb,nspin,nprojel
integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
type(locreg_descriptors), intent(in) :: lr
real(kind=8), intent(out) :: peakmem

character(len=*), parameter :: subname='MemoryEstimator'
real(kind=8), parameter :: eps_mach=1.d-12
integer :: norbp,nvctrp,n1,n2,n3
integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3
integer(kind=8) :: mworkham, mworkrho
real(kind=8) :: omemwf,omemker,omemden,omempot,omemproj,nden,npotden,npotham,narr
real(kind=8) :: tt,tmemker,tmemden,tmemps,tmemha
*/
void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(const int *nproc, 
                                               const int *idsx, 
                                               const _locreg_descriptors *lr, 
                                               const int *nat, 
                                               const int *norb, 
                                               const int *nspinor, 
                                               const int *nkpt, 
                                               const int *nprojel, 
                                               const int *nspin, 
                                               const int *itrpmax, 
                                               const int *iscf, 
                                               double *peakmem);
/* neb_free src/images.f90:1091 */
/* Fortran header:
subroutine neb_free(neb)
use module_images
implicit none
type(NEB_data), pointer :: neb
*/
void FC_FUNC_(neb_free, NEB_FREE)(_NEB_data **neb);
/* neb_new src/images.f90:1077 */
/* Fortran header:
subroutine neb_new(neb)
use module_images
implicit none
type(NEB_data), pointer :: neb
*/
void FC_FUNC_(neb_new, NEB_NEW)(_NEB_data **neb);
/* optloop_bcast src/bindings/bindingsf.f90:1339 */
/* Fortran header:
subroutine optloop_bcast(optloop, iproc)
use module_base
use module_types
implicit none
type(DFT_optimization_loop), intent(inout) :: optloop
integer, intent(in) :: iproc

integer :: iData(4), ierr
real(gp) :: rData(3)
*/
void FC_FUNC_(optloop_bcast, OPTLOOP_BCAST)(_DFT_optimization_loop *optloop, 
                                            const int *iproc);
/* optloop_copy_data src/bindings/bindingsf.f90:1243 */
/* Fortran header:
subroutine optloop_copy_data(optloop, gnrm_cv, rpnrm_cv, gnrm_startmix, gnrm, rpnrm,   itrpmax, nrepmax, itermax, itrp, itrep, iter, iscf, infocode)
use module_types
implicit none
type(DFT_optimization_loop), intent(in) :: optloop
integer, intent(out) :: iscf, itrpmax, nrepmax, itermax, itrp, itrep, iter, infocode
real(gp), intent(out) :: gnrm, rpnrm, gnrm_cv, rpnrm_cv, gnrm_startmix
*/
void FC_FUNC_(optloop_copy_data, OPTLOOP_COPY_DATA)(const _DFT_optimization_loop *optloop, 
                                                    double *gnrm_cv, 
                                                    double *rpnrm_cv, 
                                                    double *gnrm_startmix, 
                                                    double *gnrm, 
                                                    double *rpnrm, 
                                                    int *itrpmax, 
                                                    int *nrepmax, 
                                                    int *itermax, 
                                                    int *itrp, 
                                                    int *itrep, 
                                                    int *iter, 
                                                    int *iscf, 
                                                    int *infocode);
/* optloop_free src/bindings/bindingsf.f90:1236 */
/* Fortran header:
subroutine optloop_free(optloop)
use module_types
implicit none
type(DFT_optimization_loop), pointer :: optloop
*/
void FC_FUNC_(optloop_free, OPTLOOP_FREE)(_DFT_optimization_loop **optloop);
/* optloop_new src/bindings/bindingsf.f90:1227 */
/* Fortran header:
subroutine optloop_new(self, optloop)
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_optimization_loop), pointer :: optloop
*/
void FC_FUNC_(optloop_new, OPTLOOP_NEW)(const long *self, 
                                        _DFT_optimization_loop **optloop);
/* optloop_sync_data src/bindings/bindingsf.f90:1266 */
/* Fortran header:
subroutine optloop_sync_data(optloop, gnrm_cv, rpnrm_cv, gnrm_startmix, gnrm, rpnrm,   itrpmax, nrepmax, itermax, itrp, itrep, iter, iscf, infocode)
use module_types
implicit none
type(DFT_optimization_loop), intent(inout) :: optloop
integer, intent(in) :: iscf, itrpmax, nrepmax, itermax, itrp, itrep, iter, infocode
real(gp), intent(in) :: gnrm, rpnrm, gnrm_cv, rpnrm_cv, gnrm_startmix
*/
void FC_FUNC_(optloop_sync_data, OPTLOOP_SYNC_DATA)(_DFT_optimization_loop *optloop, 
                                                    const double *gnrm_cv, 
                                                    const double *rpnrm_cv, 
                                                    const double *gnrm_startmix, 
                                                    const double *gnrm, 
                                                    const double *rpnrm, 
                                                    const int *itrpmax, 
                                                    const int *nrepmax, 
                                                    const int *itermax, 
                                                    const int *itrp, 
                                                    const int *itrep, 
                                                    const int *iter, 
                                                    const int *iscf, 
                                                    const int *infocode);
/* orbs_comm_empty src/bindings/bindingsf.f90:694 */
/* Fortran header:
subroutine orbs_comm_empty(comms)
use module_base
use module_types
use module_interfaces
implicit none
type(communications_arrays), intent(inout) :: comms
*/
void FC_FUNC_(orbs_comm_empty, ORBS_COMM_EMPTY)(_communications_arrays *comms);
/* orbs_comm_free src/bindings/bindingsf.f90:685 */
/* Fortran header:
subroutine orbs_comm_free(comms)
use module_base
use module_types
use module_interfaces
implicit none
type(communications_arrays), pointer :: comms
*/
void FC_FUNC_(orbs_comm_free, ORBS_COMM_FREE)(_communications_arrays **comms);
/* orbs_comm_init src/bindings/bindingsf.f90:673 */
/* Fortran header:
subroutine orbs_comm_init(comms, orbs, lr, iproc, nproc)
use module_base
use module_types
use module_interfaces
implicit none
integer, intent(in) :: iproc,nproc
type(locreg_descriptors), intent(in) :: lr
type(orbitals_data), intent(inout) :: orbs
type(communications_arrays), intent(inout) :: comms
*/
void FC_FUNC_(orbs_comm_init, ORBS_COMM_INIT)(_communications_arrays *comms, 
                                              _orbitals_data *orbs, 
                                              const _locreg_descriptors *lr, 
                                              const int *iproc, 
                                              const int *nproc);
/* orbs_comm_new src/bindings/bindingsf.f90:663 */
/* Fortran header:
subroutine orbs_comm_new(comms)
use module_base
use module_types
use module_interfaces
implicit none
type(communications_arrays), pointer :: comms
*/
void FC_FUNC_(orbs_comm_new, ORBS_COMM_NEW)(_communications_arrays **comms);
/* orbs_empty src/bindings/bindingsf.f90:656 */
/* Fortran header:
subroutine orbs_empty(orbs)
use module_types
implicit none
type(orbitals_data), intent(inout) :: orbs
*/
void FC_FUNC_(orbs_empty, ORBS_EMPTY)(_orbitals_data *orbs);
/* orbs_free src/bindings/bindingsf.f90:648 */
/* Fortran header:
subroutine orbs_free(orbs)
use module_types
use m_profiling
implicit none
type(orbitals_data), pointer :: orbs
*/
void FC_FUNC_(orbs_free, ORBS_FREE)(_orbitals_data **orbs);
/* orbs_get_dimensions src/bindings/bindingsf.f90:705 */
/* Fortran header:
subroutine orbs_get_dimensions(orbs, norb, norbp, norbu, norbd, nspin, nspinor, npsidim,  nkpts, nkptsp, isorb, iskpts)
use module_types
implicit none
type(orbitals_data), intent(in) :: orbs
integer, intent(out) :: norb, norbp, norbu, norbd, nspin, nspinor, npsidim,  nkpts, nkptsp, isorb, iskpts
*/
void FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(const _orbitals_data *orbs, 
                                                        int *norb, 
                                                        int *norbp, 
                                                        int *norbu, 
                                                        int *norbd, 
                                                        int *nspin, 
                                                        int *nspinor, 
                                                        int *npsidim, 
                                                        int *nkpts, 
                                                        int *nkptsp, 
                                                        int *isorb, 
                                                        int *iskpts);
/* orbs_get_iorbp src/bindings/bindingsf.f90:1132 */
/* Fortran header:
subroutine orbs_get_iorbp(orbs, iorbp, isorb, iproc, ikpt, iorb, ispin, ispinor)
use module_types
implicit none

integer, intent(out) :: iorbp, isorb, iproc
type(orbitals_data), intent(in) :: orbs
integer, intent(in) :: ikpt, iorb, ispin, ispinor
*/
void FC_FUNC_(orbs_get_iorbp, ORBS_GET_IORBP)(const _orbitals_data *orbs, 
                                              int *iorbp, 
                                              int *isorb, 
                                              int *iproc, 
                                              const int *ikpt, 
                                              const int *iorb, 
                                              const int *ispin, 
                                              const int *ispinor);
/* orbs_init src/bindings/bindingsf.f90:641 */
/* Fortran header:
subroutine orbs_init(orbs)
use module_types
implicit none
type(orbitals_data), intent(inout) :: orbs
*/
void FC_FUNC_(orbs_init, ORBS_INIT)(_orbitals_data *orbs);
/* orbs_new src/bindings/bindingsf.f90:634 */
/* Fortran header:
subroutine orbs_new(orbs)
use module_types
implicit none
type(orbitals_data), pointer :: orbs
*/
void FC_FUNC_(orbs_new, ORBS_NEW)(_orbitals_data **orbs);
/* orbs_open_file src/bindings/bindingsf.f90:781 */
/* Fortran header:
subroutine orbs_open_file(orbs, unitwf, name, ln, iformat, iorbp, ispinor)
use module_types
use module_interfaces
implicit none
type(orbitals_data), intent(in) :: orbs
integer, intent(in) :: unitwf, ln, iformat, iorbp, ispinor
character(len = 1), dimension(ln), intent(in) :: name

character(len = ln) :: filename
integer :: i, iorb_out
*/
void FC_FUNC_(orbs_open_file, ORBS_OPEN_FILE)(const _orbitals_data *orbs, 
                                              const int *unitwf, 
                                              const char *name, 
                                              const int *ln, 
                                              const int *iformat, 
                                              const int *iorbp, 
                                              const int *ispinor, 
                                              int str_ln_1);
/* proj_free src/bindings/bindingsf.f90:807 */
/* Fortran header:
subroutine proj_free(nlpspd, proj)
use module_types
use m_profiling
implicit none
type(nonlocal_psp_descriptors), pointer :: nlpspd
real(kind=8), dimension(:), pointer :: proj

integer :: i_stat, i_all
*/
void FC_FUNC_(proj_free, PROJ_FREE)(_nonlocal_psp_descriptors **nlpspd, 
                                    f90_pointer_double *proj);
/* proj_get_dimensions src/bindings/bindingsf.f90:821 */
/* Fortran header:
subroutine proj_get_dimensions(nlpspd, nproj, nprojel)
use module_types
implicit none
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
integer, intent(out) :: nproj, nprojel
*/
void FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(const _nonlocal_psp_descriptors *nlpspd, 
                                                        int *nproj, 
                                                        int *nprojel);
/* proj_new src/bindings/bindingsf.f90:800 */
/* Fortran header:
subroutine proj_new(nlpspd)
use module_types
implicit none
type(nonlocal_psp_descriptors), pointer :: nlpspd
*/
void FC_FUNC_(proj_new, PROJ_NEW)(_nonlocal_psp_descriptors **nlpspd);
/* read_orbital_variables src/init/sysprop.f90:673 */
/* Fortran header:
subroutine read_orbital_variables(iproc,nproc,verb,in,atoms,orbs)
use module_base
use module_types
use module_interfaces
use yaml_output
implicit none
type(input_variables), intent(in) :: in
integer, intent(in) :: iproc,nproc
logical, intent(in) :: verb
type(atoms_data), intent(in) :: atoms
type(orbitals_data), intent(inout) :: orbs

character(len=*), parameter :: subname='read_orbital_variables'
integer, parameter :: nelecmax=32,nmax=6,lmax=4,noccmax=2
logical :: exists
integer :: iat,iunit,norb,norbu,norbd,nspinor,jpst,norbme,norbyou,jproc,ikpts
integer :: norbuempty,norbdempty,nelec
integer :: nt,ntu,ntd,ityp,ierror,ispinsum
integer :: ispol,ichg,ichgsum,norbe,norbat,nspin
integer, dimension(lmax) :: nl
real(gp), dimension(noccmax,lmax) :: occup
character(len=100) :: radical
*/
void FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(const int *iproc, 
                                                              const int *nproc, 
                                                              const int *verb, 
                                                              const _input_variables *in, 
                                                              const _atoms_data *atoms, 
                                                              _orbitals_data *orbs);
/* read_radii_variables src/init/sysprop.f90:622 */
/* Fortran header:
subroutine read_radii_variables(atoms, radii_cf, crmult, frmult, projrad)
use module_base
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), intent(in) :: crmult, frmult, projrad
real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf

integer, parameter :: nelecmax=32,nmax=6,lmax=4
character(len=2) :: symbol
integer :: i,ityp,mxpl,mxchg,nsccode
real(gp) :: rcov,rprb,ehomo,radfine,amu,maxrad
real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
*/
void FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(const _atoms_data *atoms, 
                                                          double *radii_cf, 
                                                          const double *crmult, 
                                                          const double *frmult, 
                                                          const double *projrad);
/* read_wave_descr src/restart.f90:699 */
/* Fortran header:
subroutine read_wave_descr(lstat, filename, ln,  norbu, norbd, iorbs, ispins, nkpt, ikpts, nspinor, ispinor)
use module_types
implicit none
integer, intent(in) :: ln
character, intent(in) :: filename(ln)
integer, intent(out) :: norbu, norbd, nkpt, nspinor
integer, intent(out) :: iorbs, ispins, ikpts, ispinor
logical, intent(out) :: lstat

character(len = 1024) :: filename_
integer :: wave_format_from_filename, iformat, i
character(len = 1024) :: testf
*/
void FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)(int *lstat, 
                                                const char *filename, 
                                                const int *ln, 
                                                int *norbu, 
                                                int *norbd, 
                                                int *iorbs, 
                                                int *ispins, 
                                                int *nkpt, 
                                                int *ikpts, 
                                                int *nspinor, 
                                                int *ispinor, 
                                                int str_ln_1);
/* read_wave_to_isf src/restart.f90:650 */
/* Fortran header:
subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz,  n1, n2, n3, nspinor, psiscf)
use module_base
use module_types
use module_interfaces, except_this_one => read_wave_to_isf

implicit none

integer, intent(in) :: ln
character, intent(in) :: filename(ln)
integer, intent(in) :: iorbp
integer, intent(out) :: n1, n2, n3, nspinor
real(gp), intent(out) :: hx, hy, hz
real(wp), dimension(:,:,:,:), pointer :: psiscf
logical, intent(out) :: lstat

character(len = 1024) :: filename_
integer :: wave_format_from_filename, iformat, i
*/
void FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)(int *lstat, 
                                                  const char *filename, 
                                                  const int *ln, 
                                                  const int *iorbp, 
                                                  double *hx, 
                                                  double *hy, 
                                                  double *hz, 
                                                  int *n1, 
                                                  int *n2, 
                                                  int *n3, 
                                                  int *nspinor, 
                                                  f90_pointer_double_4D *psiscf, 
                                                  int str_ln_1);
/* rst_free src/bindings/bindingsf.f90:1384 */
/* Fortran header:
subroutine rst_free(rst)
use module_types
implicit none
type(restart_objects), pointer :: rst
*/
void FC_FUNC_(rst_free, RST_FREE)(_restart_objects **rst);
/* rst_init src/bindings/bindingsf.f90:1392 */
/* Fortran header:
subroutine rst_init(rst, iproc, atoms, inputs)
use module_types
implicit none
type(restart_objects), intent(out) :: rst
integer, intent(in) :: iproc
type(atoms_data), intent(in) :: atoms
type(input_variables), intent(in) :: inputs
*/
void FC_FUNC_(rst_init, RST_INIT)(_restart_objects *rst, 
                                  const int *iproc, 
                                  const _atoms_data *atoms, 
                                  const _input_variables *inputs);
/* rst_new src/bindings/bindingsf.f90:1375 */
/* Fortran header:
subroutine rst_new(self, rst)
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(restart_objects), pointer :: rst
*/
void FC_FUNC_(rst_new, RST_NEW)(const long *self, 
                                _restart_objects **rst);
/* run_objects_associate src/bindings/bindingsf.f90:1439 */
/* Fortran header:
subroutine run_objects_associate(runObj, inputs, atoms, rst, rxyz)
use module_types
implicit none
type(run_objects), intent(out) :: runObj
type(input_variables), intent(in), target :: inputs
type(atoms_data), intent(in), target :: atoms
type(restart_objects), intent(in), target :: rst
real(gp), dimension(:,:), pointer :: rxyz
*/
void FC_FUNC_(run_objects_associate, RUN_OBJECTS_ASSOCIATE)(_run_objects *runObj, 
                                                            const _input_variables *inputs, 
                                                            const _atoms_data *atoms, 
                                                            const _restart_objects *rst, 
                                                            f90_pointer_double_2D *rxyz);
/* run_objects_destroy src/bindings/bindingsf.f90:1411 */
/* Fortran header:
subroutine run_objects_destroy(runObj)
use module_types
use module_base
use yaml_output
implicit none
type(run_objects), pointer :: runObj
*/
void FC_FUNC_(run_objects_destroy, RUN_OBJECTS_DESTROY)(_run_objects **runObj);
/* run_objects_get src/bindings/bindingsf.f90:1425 */
/* Fortran header:
subroutine run_objects_get(runObj, inputs, atoms, rst, rxyz)
use module_types
implicit none
type(run_objects), intent(in) :: runObj
type(input_variables), pointer :: inputs
type(atoms_data), pointer :: atoms
type(restart_objects), pointer :: rst
real(gp), dimension(:,:), pointer :: rxyz
*/
void FC_FUNC_(run_objects_get, RUN_OBJECTS_GET)(const _run_objects *runObj, 
                                                _input_variables **inputs, 
                                                _atoms_data **atoms, 
                                                _restart_objects **rst, 
                                                f90_pointer_double_2D *rxyz);
/* run_objects_init_from_files src/cluster.f90:275 */
/* Fortran header:
subroutine run_objects_init_from_files(runObj, radical, posinp)
use module_interfaces
use module_types
implicit none
type(run_objects), intent(out) :: runObj
character(len = *), intent(in) :: radical, posinp
*/
void FC_FUNC_(run_objects_init_from_files, RUN_OBJECTS_INIT_FROM_FILES)(_run_objects *runObj, 
                                                                        const char *radical, 
                                                                        const char *posinp, 
                                                                        int str_ln_1, 
                                                                        int str_ln_2);
/* run_objects_new src/bindings/bindingsf.f90:1403 */
/* Fortran header:
subroutine run_objects_new(runObj)
use module_types
implicit none
type(run_objects), pointer :: runObj
*/
void FC_FUNC_(run_objects_new, RUN_OBJECTS_NEW)(_run_objects **runObj);
/* standard_inputfile_names src/input_variables.f90:220 */
/* Fortran header:
subroutine standard_inputfile_names(in, radical, nproc)
use module_types
use module_base
use yaml_output
implicit none
type(input_variables), intent(inout) :: in
character(len = *), intent(in) :: radical
integer, intent(in) :: nproc
integer :: ierr
*/
void FC_FUNC_(standard_inputfile_names, STANDARD_INPUTFILE_NAMES)(_input_variables *in, 
                                                                  const char *radical, 
                                                                  const int *nproc, 
                                                                  int str_ln_1);
/* symmetry_set_irreductible_zone src/init/atoms.f90:1948 */
/* Fortran header:
subroutine symmetry_set_irreductible_zone(sym, geocode, n1i, n2i, n3i, nspin)
use module_base
use module_types
use m_ab6_kpoints
use m_ab6_symmetry
implicit none
type(symmetry_data), intent(inout) :: sym
integer, intent(in) :: n1i, n2i, n3i, nspin
character, intent(in) :: geocode

character(len = *), parameter :: subname = "symmetry_set_irreductible_zone"
integer :: i_stat, nsym, i_all, i_third
integer, dimension(:,:,:), allocatable :: irrzon
real(dp), dimension(:,:,:), allocatable :: phnons
*/
void FC_FUNC_(symmetry_set_irreductible_zone, SYMMETRY_SET_IRREDUCTIBLE_ZONE)(_symmetry_data *sym, 
                                                                              const char *geocode, 
                                                                              const int *n1i, 
                                                                              const int *n2i, 
                                                                              const int *n3i, 
                                                                              const int *nspin, 
                                                                              int str_ln_1);
/* system_createkernels src/init/sysprop.f90:238 */
/* Fortran header:
subroutine system_createKernels(denspot, verb)
use module_types
use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
implicit none
logical, intent(in) :: verb
type(DFT_local_fields), intent(inout) :: denspot
*/
void FC_FUNC_(system_createkernels, SYSTEM_CREATEKERNELS)(_DFT_local_fields *denspot, 
                                                          const int *verb);
/* system_initkernels src/init/sysprop.f90:211 */
/* Fortran header:
subroutine system_initKernels(verb, iproc, nproc, geocode, in, denspot)
use module_types
use module_xc
use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
implicit none
logical, intent(in) :: verb
integer, intent(in) :: iproc, nproc
character, intent(in) :: geocode
type(input_variables), intent(in) :: in
type(DFT_local_fields), intent(inout) :: denspot

integer, parameter :: ndegree_ip = 16
*/
void FC_FUNC_(system_initkernels, SYSTEM_INITKERNELS)(const int *verb, 
                                                      const int *iproc, 
                                                      const int *nproc, 
                                                      const char *geocode, 
                                                      const _input_variables *in, 
                                                      _DFT_local_fields *denspot, 
                                                      int str_ln_1);
/* system_size src/init/gridmanipulation.f90:14 */
/* Fortran header:
subroutine system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
use module_base
use module_types
use yaml_output
implicit none
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: iproc
real(gp), intent(in) :: crmult,frmult
real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
real(gp), intent(inout) :: hx,hy,hz
type(locreg_descriptors), intent(out) :: Glr
real(gp), dimension(3), intent(out) :: shift

integer, parameter :: lupfil=14
real(gp), parameter ::eps_mach=1.e-12_gp
integer :: iat,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
real(gp) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3
character(len=*), parameter :: subname='system_size'
*/
void FC_FUNC_(system_size, SYSTEM_SIZE)(const int *iproc, 
                                        _atoms_data *atoms, 
                                        double *rxyz, 
                                        const double *radii_cf, 
                                        const double *crmult, 
                                        const double *frmult, 
                                        double *hx, 
                                        double *hy, 
                                        double *hz, 
                                        _locreg_descriptors *Glr, 
                                        double *shift);
/* update_wavefunctions_size src/linear/initAndUtils.f90:1027 */
/* Fortran header:
subroutine update_wavefunctions_size(lzd,npsidim_orbs,npsidim_comp,orbs,iproc,nproc)
use module_base
use module_types
implicit none


type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(inout) :: orbs
integer, intent(in) :: iproc, nproc
integer, intent(out) :: npsidim_orbs, npsidim_comp


integer :: npsidim, ilr, iorb
integer :: nvctr_tot,jproc,istat,ierr,iall
integer, allocatable, dimension(:) :: ncntt 
integer, allocatable, dimension(:,:) :: nvctr_par
character(len = *), parameter :: subname = "update_wavefunctions_size"
*/
void FC_FUNC_(update_wavefunctions_size, UPDATE_WAVEFUNCTIONS_SIZE)(const _local_zone_descriptors *lzd, 
                                                                    int *npsidim_orbs, 
                                                                    int *npsidim_comp, 
                                                                    _orbitals_data *orbs, 
                                                                    const int *iproc, 
                                                                    const int *nproc);
/* wf_empty src/bindings/bindingsf.f90:1050 */
/* Fortran header:
subroutine wf_empty(wf)
use module_types
use m_profiling
implicit none
type(DFT_wavefunction), intent(inout) :: wf

integer :: i_all, i_stat
*/
void FC_FUNC_(wf_empty, WF_EMPTY)(_DFT_wavefunction *wf);
/* wf_free src/bindings/bindingsf.f90:1074 */
/* Fortran header:
subroutine wf_free(wf)
use module_types
use m_profiling
implicit none
type(DFT_wavefunction), pointer :: wf
*/
void FC_FUNC_(wf_free, WF_FREE)(_DFT_wavefunction **wf);
/* wf_get_data src/bindings/bindingsf.f90:1038 */
/* Fortran header:
subroutine wf_get_data(wf, orbs, comm, lzd)
use module_types
implicit none
type(DFT_wavefunction), target, intent(in) :: wf
type(orbitals_data), pointer :: orbs
type(communications_arrays), pointer :: comm
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(wf_get_data, WF_GET_DATA)(const _DFT_wavefunction *wf, 
                                        _orbitals_data **orbs, 
                                        _communications_arrays **comm, 
                                        _local_zone_descriptors **lzd);
/* wf_get_psi src/bindings/bindingsf.f90:1085 */
/* Fortran header:
subroutine wf_get_psi(wf, psi, hpsi)
use module_types
implicit none
type(DFT_wavefunction), intent(in) :: wf
integer(kind = 8), intent(out) :: psi
integer(kind = 8), intent(out) :: hpsi
*/
void FC_FUNC_(wf_get_psi, WF_GET_PSI)(const _DFT_wavefunction *wf, 
                                      long *psi, 
                                      long *hpsi);
/* wf_get_psi_size src/bindings/bindingsf.f90:1101 */
/* Fortran header:
subroutine wf_get_psi_size(psi, psiSize)
use module_types
implicit none
real(wp), dimension(:), pointer :: psi
integer(kind = 8), intent(out) :: psiSize
*/
void FC_FUNC_(wf_get_psi_size, WF_GET_PSI_SIZE)(f90_pointer_double *psi, 
                                                long *psiSize);
/* wf_iorbp_to_psi src/bindings/bindingsf.f90:1109 */
/* Fortran header:
subroutine wf_iorbp_to_psi(psir, psi, lr)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: lr
real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psir

character(len=*), parameter :: subname='wf_orb_to_psi'
type(workarr_sumrho) :: w
*/
void FC_FUNC_(wf_iorbp_to_psi, WF_IORBP_TO_PSI)(double *psir, 
                                                const double *psi, 
                                                const _locreg_descriptors *lr);
/* wf_new src/bindings/bindingsf.f90:1009 */
/* Fortran header:
subroutine wf_new(self, wf, orbs, comm, lzd)
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_wavefunction), pointer :: wf
type(orbitals_data), pointer :: orbs
type(communications_arrays), pointer :: comm
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(wf_new, WF_NEW)(const long *self, 
                              _DFT_wavefunction **wf, 
                              _orbitals_data **orbs, 
                              _communications_arrays **comm, 
                              _local_zone_descriptors **lzd);
/* write_extra_info src/init/atoms.f90:1413 */
/* Fortran header:
subroutine write_extra_info(extra,natpol,ifrztyp)
use module_base
implicit none 
integer, intent(in) :: natpol,ifrztyp
character(len=50), intent(out) :: extra

character(len=4) :: frzchain
integer :: ispol,ichg
*/
void FC_FUNC_(write_extra_info, WRITE_EXTRA_INFO)(char *extra, 
                                                  const int *natpol, 
                                                  const int *ifrztyp, 
                                                  int str_ln_1);
/* writeonewave src/wavelib/i-o.f90:901 */
/* Fortran header:
subroutine writeonewave(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,nat,rxyz,  nseg_c,nvctr_c,keyg_c,keyv_c,  nseg_f,nvctr_f,keyg_f,keyv_f, psi_c,psi_f,eval)
use module_base
use yaml_output
implicit none
logical, intent(in) :: useFormattedOutput
integer, intent(inout) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f
real(gp), intent(in) :: hx,hy,hz
real(wp), intent(in) :: eval
integer, dimension(nseg_c), intent(in) :: keyv_c
integer, dimension(nseg_f), intent(in) :: keyv_f
integer, dimension(2,nseg_c), intent(in) :: keyg_c
integer, dimension(2,nseg_f), intent(in) :: keyg_f
real(wp), dimension(nvctr_c), intent(in) :: psi_c
real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
real(gp), dimension(3,nat), intent(in) :: rxyz

integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
*/
void FC_FUNC(writeonewave, WRITEONEWAVE)(int *unitwf, 
                                         const int *useFormattedOutput, 
                                         int *iorb, 
                                         int *n1, 
                                         int *n2, 
                                         int *n3, 
                                         const double *hx, 
                                         const double *hy, 
                                         const double *hz, 
                                         int *nat, 
                                         const double *rxyz, 
                                         int *nseg_c, 
                                         int *nvctr_c, 
                                         const int *keyg_c, 
                                         const int *keyv_c, 
                                         int *nseg_f, 
                                         int *nvctr_f, 
                                         const int *keyg_f, 
                                         const int *keyv_f, 
                                         const double *psi_c, 
                                         const double *psi_f, 
                                         const double *eval);
/* writeonewave_linear src/restart.f90:737 */
/* Fortran header:
subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,locregCenter,locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  nseg_f,nvctr_f,keyg_f,keyv_f, psi_c,psi_f,eval,onwhichatom)
use module_base
use yaml_output
implicit none
logical, intent(in) :: useFormattedOutput
integer, intent(in) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
real(gp), intent(in) :: hx,hy,hz,locrad,confPotprefac
real(wp), intent(in) :: eval
integer, dimension(nseg_c), intent(in) :: keyv_c
integer, dimension(nseg_f), intent(in) :: keyv_f
integer, dimension(2,nseg_c), intent(in) :: keyg_c
integer, dimension(2,nseg_f), intent(in) :: keyg_f
real(wp), dimension(nvctr_c), intent(in) :: psi_c
real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
real(gp), dimension(3,nat), intent(in) :: rxyz
real(gp), dimension(3), intent(in) :: locregCenter
integer, intent(in) :: onwhichatom

integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
*/
void FC_FUNC_(writeonewave_linear, WRITEONEWAVE_LINEAR)(const int *unitwf, 
                                                        const int *useFormattedOutput, 
                                                        const int *iorb, 
                                                        const int *n1, 
                                                        const int *n2, 
                                                        const int *n3, 
                                                        const double *hx, 
                                                        const double *hy, 
                                                        const double *hz, 
                                                        const double *locregCenter, 
                                                        const double *locrad, 
                                                        const int *confPotOrder, 
                                                        const double *confPotprefac, 
                                                        const int *nat, 
                                                        const double *rxyz, 
                                                        const int *nseg_c, 
                                                        const int *nvctr_c, 
                                                        const int *keyg_c, 
                                                        const int *keyv_c, 
                                                        const int *nseg_f, 
                                                        const int *nvctr_f, 
                                                        const int *keyg_f, 
                                                        const int *keyv_f, 
                                                        const double *psi_c, 
                                                        const double *psi_f, 
                                                        const double *eval, 
                                                        const int *onwhichatom);
#endif
