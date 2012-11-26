#ifndef BINDINGS_API_H
#define BINDINGS_API_H

void FC_FUNC_(atoms_new, ATOMS_NEW)(void *atoms, void *sym);
void FC_FUNC_(atoms_free, ATOMS_FREE)(void *atoms);
void FC_FUNC_(atoms_empty, ATOMS_EMPTY)(void *atoms);
void FC_FUNC_(atoms_set_from_file, ATOMS_SET_FROM_FILE)(int *lstat, void *atoms,
                                                        f90_pointer_double_2D *rxyz,
                                                        const gchar *filename, int *ln);
void FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(void *atoms,
                                                    f90_pointer_double_2D *rxyz, int *nat);
void FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(void *atoms, int *ntypes);
void FC_FUNC_(atoms_set_symmetries, ATOMS_SET_SYMMETRIES)(void *atoms, double *rxyz,
                                                          int *disable, double *tol,
                                                          double *elecfield);
void FC_FUNC_(atoms_set_displacement, ATOMS_SET_DISPLACEMENT)(void *atoms, double *rxyz,
                                                              double *randdis);
void FC_FUNC_(atoms_set_name, ATOMS_SET_NAME)(void *atoms, int *ityp, gchar *name);
void FC_FUNC_(atoms_sync, ATOMS_SYNC)(void *atoms, double *alat1, double *alat2,
                                      double *alat3, gchar *geocode, gchar *format,
                                      gchar *units);
void FC_FUNC_(atoms_get_geocode, ATOMS_GET_GEOCODE)(void *atoms, char *geocode);
void FC_FUNC_(atoms_get_iatype, ATOMS_GET_IATYPE)(void *atoms, f90_pointer_int *iatype);
void FC_FUNC_(atoms_get_iasctype, ATOMS_GET_IASCTYPE)(void *atoms, f90_pointer_int *iasctype);
void FC_FUNC_(atoms_get_natpol, ATOMS_GET_NATPOL)(void *atoms, f90_pointer_int *natpol);
void FC_FUNC_(atoms_get_ifrztyp, ATOMS_GET_IFRZTYP)(void *atoms, f90_pointer_int *ifrztyp);
void FC_FUNC_(atoms_get_amu, ATOMS_GET_AMU)(void *atoms, f90_pointer_double *amu);
void FC_FUNC_(atoms_get_aocc, ATOMS_GET_AOCC)(void *atoms, f90_pointer_double_2D *aocc);
void FC_FUNC_(atoms_get_nelpsp, ATOMS_GET_NELPSP)(void *atoms, f90_pointer_int *nelpsp);
void FC_FUNC_(atoms_get_npspcode, ATOMS_GET_NPSPCODE)(void *atoms, f90_pointer_int *npspcode);
void FC_FUNC_(atoms_get_nzatom, ATOMS_GET_NZATOM)(void *atoms, f90_pointer_int *nzatom);
void FC_FUNC_(atoms_get_nlcc_ngv, ATOMS_GET_NLCC_NGV)(void *atoms, f90_pointer_int *nlcc_ngv);
void FC_FUNC_(atoms_get_nlcc_ngc, ATOMS_GET_NLCC_NGC)(void *atoms, f90_pointer_int *nlcc_ngc);
void FC_FUNC_(atoms_get_ixcpsp, ATOMS_GET_IXCPSP)(void *atoms, f90_pointer_int *ixcpsp);
void FC_FUNC_(atoms_get_radii_cf, ATOMS_GET_RADII_CF)(void *atoms, f90_pointer_double_2D *radii_cf);
void FC_FUNC_(atoms_get_psppar, ATOMS_GET_PSPPAR)(void *atoms, f90_pointer_double_3D *psppar);
void FC_FUNC_(atoms_get_nlccpar, ATOMS_GET_NLCCPAR)(void *atoms, f90_pointer_double_2D *nlccpar);
void FC_FUNC_(atoms_get_ig_nlccpar, ATOMS_GET_IG_NLCCPAR)(void *atoms, f90_pointer_double_2D *ig_nlccpar);
void FC_FUNC_(atoms_copy_nat, ATOMS_COPY_NAT)(void *atoms, int *nat);
void FC_FUNC_(atoms_copy_ntypes, ATOMS_COPY_NTYPES)(void *atoms, int *ntypes);
void FC_FUNC_(atoms_copy_geometry_data, ATOMS_COPY_GEOMETRY_DATA)
     (void *atoms, gchar *geocode, gchar *format, gchar *units, guint len_geocode,
      guint len_format, guint len_units);
void FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(void *atoms, double *alat1,
                                                double *alat2, double *alat3);
void FC_FUNC_(atoms_copy_psp_data, ATOMS_COPY_PSP_DATA)(void *atoms, int *natsc, int *donlcc);
void FC_FUNC_(atoms_copy_name, ATOMS_COPY_NAME)(void *atoms, int *ityp, gchar *name, int *ln);
void FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(int *verb, void *atoms, int *ixc);
void FC_FUNC_(atoms_read_variables, ATOMS_READ_VARIABLES)(void *atoms, const guint *nspin,
                                                          const gchar *occup, const guint *ln);
void FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(void *atoms, double *radii_cf,
                                                          double *crmult, double *frmult,
                                                          double *projrad);
void FC_FUNC_(atoms_write, ATOMS_WRITE)(void *atoms, const gchar *filename, int *ln2,
                                        double *rxyz, f90_pointer_double *forces,
                                        const double *energy, const gchar *comment, int *ln);
void FC_FUNC_(fill_logrid, FILL_LOGRID)(const char *geocode, const guint *n1, const guint *n2,
                                        const guint *n3, const guint *nl1, const guint *nu1,
                                        const guint *nl2, const guint *nu2, const guint *nl3,
                                        const guint *nu3, const guint *orig,
                                        const guint *nat, const guint *ntypes,
                                        const guint *iatype, const double *rxyz,
                                        const double *radii, const double *mult,
                                        const double *hx, const double *hy, const double *hz,
                                        int *grid);
void FC_FUNC_(symmetry_set_irreductible_zone, SYMMETRY_SET_IRREDUCTIBLE_ZONE)
     (void *sym, const gchar *geocode, const guint *n1i, const guint *n2i,
      const guint *n3i, const guint *nspin, int ln);

void FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(double *self, void *denspotd,
                                                void *rhod, void *dpcom);
void FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(void *denspotd);
void FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(const guint *iproc,
                                             const void *glr, const int *nspin,
                                             const void *atoms, const double *rxyz,
                                             void *denspotd);
void FC_FUNC(system_createkernels, SYSTEM_CREATEKERNELS)(void *denspot, const guint *verb);


void FC_FUNC_(glr_new, GLR_NEW)(void *glr);
void FC_FUNC_(glr_init, GLR_INIT)(void *glr, void *d, void *wfd);
void FC_FUNC_(glr_get_data, GLR_GET_DATA)(void *glr, void *d, void *wfd);
void FC_FUNC_(system_size, SYSTEM_SIZE)(int *iproc, void *atoms, double *rxyz,
                                        double *radii_cf, double *crmult, double *frmult,
                                        double *hx, double *hy, double *hz,
                                        void *glr, double *shift);
void FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(const void *glr, guint *n, guint *ni,
                                                      guint *ns, guint *nsi, guint *nfl,
                                                      guint *nfu, guint *norb);
void FC_FUNC_(glr_set_dimensions, GLR_SET_DIMENSIONS)(void *glr, const guint *n, const guint *ni,
                                                      const guint *ns, const guint *nsi,
                                                      const guint *nfl, const guint *nfu);
void FC_FUNC_(glr_empty, GLR_EMPTY)(void *glr);
void FC_FUNC_(glr_free, GLR_FREE)(void *glr);
void FC_FUNC_(glr_set_wave_descriptors,
             GLR_SET_WAVE_DESCRIPTORS)(int *iproc, double *hx, double *hy,
                                       double *hz, void *atoms, double *rxyz, double *radii,
                                       double *crmult, double *frmult, void *glr);
void FC_FUNC_(glr_wfd_get_data, GLR_WFD_GET_DATA)(void *wfd, guint *nvctr_c, guint *nvctr_f,
                                                  guint *nseg_c, guint *nseg_f, void *keyglob,
                                                  void *keygloc, void *keyvglob, void * keyvloc);
void FC_FUNC_(lzd_new, LZD_NEW)(void *lzd);
void FC_FUNC_(lzd_free, LZD_FREE)(void *lzd);
void FC_FUNC_(lzd_empty, LZD_EMPTY)(void *lzd);
void FC_FUNC_(lzd_init, LZD_INIT)(void *lzd, void *glr);
void FC_FUNC_(lzd_init_llr, LZD_INIT_LLR)(const guint *iproc, const guint *nproc,
                                          const void *in, const void *at, const double *rxyz,
                                          const void *orbs, void *lzd);
void FC_FUNC_(lzd_set_hgrids, LZD_SET_HGRIDS)(void *lzd, const double *hgrids);
void FC_FUNC_(lzd_get_hgrids, LZD_GET_HGRIDS)(void *lzd, double *hgrids);
void FC_FUNC_(lzd_get_data, LZD_GET_DATA)(void *lzd, void *glr);
void FC_FUNC_(lzd_get_llr, LZD_GET_LLR)(void *lzd, const guint *i, void *llr);
void FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(void *lzd, guint *nlr);
void FC_FUNC_(check_linear_and_create_lzd, CHECK_LINEAR_AND_CREATE_LZD)
     (const guint *iproc, const guint *nproc, const guint *type, void *lzd,
      const void *atoms, void *orbs, const guint *npsin, const double *rxyz);


void FC_FUNC_(orbs_new, ORBS_NEW)(void *orbs);
void FC_FUNC_(orbs_init, ORBS_INIT)(void *orbs);
void FC_FUNC_(orbs_free, ORBS_FREE)(void *orbs);
void FC_FUNC_(orbs_empty, ORBS_EMPTY)(void *orbs);
void FC_FUNC_(orbs_comm_new, ORBS_COMM_NEW)(void *comm);
void FC_FUNC_(orbs_comm_init, ORBS_COMM_INIT)(void *comm, void *orbs, void *lzd,
                                              const guint *iproc, const guint *nproc);
void FC_FUNC_(orbs_comm_free, ORBS_COMM_FREE)(void *comm);
void FC_FUNC_(orbs_comm_empty, ORBS_COMM_EMPTY)(void *comm);
void FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(guint *iproc, guint *nproc,
                                                              int *verb, void *in, void *atoms,
                                                              void *orbs, int *nelec);
void FC_FUNC_(orbs_comm, ORBS_COMM)(void *comm, void *orbs, const void *glr,
                                    const int *iproc, const int *nproc);
void FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(const void *orbs, guint *norb,
                                                        guint *norbp, guint *norbu,
                                                        guint *norbd, guint *nspin,
                                                        guint *nspinor, guint *npsidim,
                                                        guint *nkpts, guint *nkptsp,
                                                        guint *isorb, guint *iskpts);
void FC_FUNC_(orbs_get_eval, ORBS_GET_EVAL)(void *orbs, void *eval);
void FC_FUNC_(orbs_get_occup, ORBS_GET_OCCUP)(void *orbs, void *occup);
void FC_FUNC_(orbs_get_kwgts, ORBS_GET_KWGTS)(void *orbs, void *kwgts);
void FC_FUNC_(orbs_get_kpts, ORBS_GET_KPTS)(void *orbs, void *kpts);
void FC_FUNC_(orbs_get_inwhichlocreg, ORBS_GET_INWHICHLOCREG)(void *orbs, void *locreg);
void FC_FUNC_(orbs_get_onwhichmpi, ORBS_GET_ONWHICHMPI)(void *orbs, void *mpi);
void FC_FUNC_(orbs_get_onwhichatom, ORBS_GET_ONWHICHATOM)(void *orbs, void *atom);


void FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
     (int *lstat, const char* filename, int *ln, int *iorbp,
      double *hx, double *hy, double *hz,
      int *n1, int *n2, int *n3, int *nspinor, f90_pointer_double_4D *psiscf);
void FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(f90_pointer_double_4D *psiscf);

void FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)
     (int *lstat, const char* filename, int *ln, int *norbu,
      int *norbd, int *iorb, int *ispin, int *nkpt, int *ikpt, int *nspinor, int *ispinor);

void FC_FUNC_(wf_iorbp_to_psi, WF_IORBP_TO_PSI)(double *psir, const double *psic, void *glr);

void FC_FUNC_(wf_new, WF_NEW)(double *self, void *wf, void *orbs, void *comm, void *lzd);
void FC_FUNC_(wf_free, WF_FREE)(void *wf);
void FC_FUNC_(wf_empty, WF_EMPTY)(void *wf);
void FC_FUNC_(wf_get_psi, WF_GET_PSI)(void *wf, void *psi, void *hpsi);
void FC_FUNC_(wf_get_data, WF_GET_DATA)(void *wf, void *orbs, void *comm, void *lzd);
void FC_FUNC_(input_wf, INPUT_WF)(const guint *iproc, const guint *nproc,
                                  const void *in, const void *GPU,
                                  const void *atoms, const double *rxyz,
                                  void *denspot, const double *denspot0, const void *nlpspd,
                                  const f90_pointer_double *proj, void *wf,
                                  void *tmb, void *tmbl, void *energs,
                                  const int *inputpsi, const guint *input_wf_format, guint *norbv,
                                  void *lzd_old, void *wfd_old, void *phi_old, void *coeff_old,
                                  void *psi_old, void *d_old, const double *hx_old,
                                  const double *hy_old, const double *hz_old, double *rxyz_old);

void FC_FUNC_(energs_new, ENERGS_NEW)(double *self, void *energs);
void FC_FUNC_(energs_free, ENERGS_FREE)(void *energs);

void FC_FUNC_(optloop_new, OPTLOOP_NEW)(double *self, void *optloop);
void FC_FUNC_(optloop_free, OPTLOOP_FREE)(void *optloop);
void FC_FUNC_(optloop_copy_data, OPTLOOP_COPY_DATA)(void *optloop, double *gnrm_cv,
                                                    double *rpnrm_cv, double *gnrm_startmix,
                                                    double *gnrm, double *rpnrm,
                                                    guint *itrpmax, guint *nrepmax, guint *itermax,
                                                    guint *itrp, guint *itrep, guint *iter,
                                                    int *iscf, int *infocode);
void FC_FUNC_(optloop_sync_data, OPTLOOP_SYNC_DATA)(void *optloop, double *gnrm_cv,
                                                    double *rpnrm_cv, double *gnrm_startmix,
                                                    double *gnrm, double *rpnrm,
                                                    guint *itrpmax, guint *nrepmax, guint *itermax,
                                                    guint *itrp, guint *itrep, guint *iter,
                                                    int *iscf, int *infocode);

void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(const guint *nproc, const guint *idsx,
                                               const void *lr,
                                               const guint *nat, const guint *norb,
                                               const guint *nspinor, const guint *nkpt,
                                               const guint *nprojel, const guint *nspin,
                                               const guint *itrpmax, const guint *iscf,
                                               double *peak);

void FC_FUNC_(inputs_new, INPUTS_NEW)(void *in);
void FC_FUNC_(inputs_free, INPUTS_FREE)(void *in);
void FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(void *in, int *nproc,
                                                      const gchar *rad, int *len);
void FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(const void *in, double *hx, double *hy, double *hz, double *crmult, double *frmult, int *ixc, int *ncharge, double *elecfield, int *nspin, int *mpol, double *gnrm_cv, guint *itermax, guint *nrepmax, int *ncong, guint *idsx, int *dispersion, int *inputPsiId, int *output_wf_format, int *output_grid, double *rbuf, int *ncongt, int *norbv, int *nvirt, int *nplot, int *disableSym);
void FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(void *in, guint *iscf, guint *itrpmax,
                                              int *norbsempty, int *occopt, double *alphamix,
                                              double *rpnrm_cv, double *gnrm_startmix,
                                              double *Tel, double *alphadiis);
void FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(void *in, char *geopt_approach,
                                                  int *ncount_cluster_x, double *frac_fluct,
                                                  double *forcemax, double *randdis,
                                                  double *betax, int *history, int *ionmov,
                                                  double *dtion, double *strtarget,
                                                  f90_pointer_double *qmass);
void FC_FUNC_(inputs_get_perf, INPUTS_GET_PERF)(void *in, guint *linear);
void FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(void *in,
                                                        int *iproc, int *dump);
void FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(const void *in, int *files);
void FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(void *in, const void *atoms,
                                                  int *iproc, int *dump);

#endif
