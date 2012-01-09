#include "bigdft.h"

#include <config.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_DEBUG
#define DBG_MEM(A, T) {int i__; for (i__ = 0; i__ < sizeof(T) / sizeof(void*); i__++) \
                                  fprintf(stderr, "DBG (%2d) -> %p\n", i__, ((void**)A)[i__]); }
#else
#define DBG_MEM(A, T)
#endif

void FC_FUNC_(deallocate_double, DEALLOCATE_DOUBLE)(f90_pointer_double *array);

void FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
     (int *lstat, const char* filename, int *ln, int *iorbp,
      double *hx, double *hy, double *hz,
      int *n1, int *n2, int *n3, int *nspinor, f90_pointer_double *psiscf);
void FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(f90_pointer_double *psiscf);

void FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)
     (int *lstat, const char* filename, int *ln, int *norbu,
      int *norbd, int *iorb, int *ispin, int *nkpt, int *ikpt, int *nspinor, int *ispinor);

f90_pointer_double* bigdft_read_wave_to_isf(const char *filename, int iorbp,
                                            double h[3], int n[3], int *nspinor)
{
  int ln, lstat;
  f90_pointer_double *psiscf;

  psiscf = g_malloc(sizeof(f90_pointer_double));
  memset(psiscf, 0, sizeof(f90_pointer_double));
  
  ln = strlen(filename);
  FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
    (&lstat, filename, &ln, &iorbp, h, h + 1, h + 2, n, n + 1, n + 2, nspinor, psiscf);
  if (!lstat)
    {
      g_free(psiscf);
      psiscf = (f90_pointer_double*)0;
    }

  DBG_MEM(psiscf, f90_pointer_double);

  return psiscf;
}
void bigdft_free_wave_to_isf(f90_pointer_double *psiscf)
{
  FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(psiscf);
  g_free(psiscf);
}

gboolean bigdft_read_wave_descr(const char *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor)
{
  int ln, lstat, norbu_, norbd_, nkpt_, nspinor_;
  int iorb_, ispin_, ikpt_, ispinor_;
  
  ln = strlen(filename);
  FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)
    (&lstat, filename, &ln, &norbu_, &norbd_, &iorb_, &ispin_,
     &nkpt_, &ikpt_, &nspinor_, &ispinor_);
  if (!lstat)
    return FALSE;

  if (norbu)   *norbu   = norbu_;
  if (norbd)   *norbd   = norbd_;
  if (nkpt)    *nkpt    = nkpt_;
  if (nspinor) *nspinor = nspinor_;

  if (iorb)    *iorb    = iorb_;
  if (ispin)   *ispin   = ispin_;
  if (ikpt)    *ikpt    = ikpt_;
  if (ispinor) *ispinor = ispinor_;

  return TRUE;
}

/* The atoms_data structure binding. */
struct f90_pointer_atoms_
{
  void *atoms;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(atoms_new, ATOMS_NEW)(f90_pointer_atoms *atoms);
void FC_FUNC_(atoms_free, ATOMS_FREE)(f90_pointer_atoms *atoms);
void FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(f90_pointer_atoms *atoms, int *nat);
void FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(f90_pointer_atoms *atoms, int *ntypes);
void FC_FUNC_(atoms_get_geocode, ATOMS_GET_GEOCODE)(f90_pointer_atoms *atoms, char *geocode);
void FC_FUNC_(atoms_get_nat, ATOMS_GET_NAT)(f90_pointer_atoms *atoms, int *nat);
void FC_FUNC_(atoms_get_ntypes, ATOMS_GET_NTYPES)(f90_pointer_atoms *atoms, int *ntypes);
void FC_FUNC_(atoms_get_iatype, ATOMS_GET_IATYPE)(f90_pointer_atoms *atoms,
                                                  f90_pointer_int *iatype);
void FC_FUNC_(atoms_get_name, ATOMS_GET_NAME)(f90_pointer_atoms *atoms, int *ityp,
                                              gchar *name, int *ln);
void FC_FUNC_(atoms_get_alat, ATOMS_GET_ALAT)(f90_pointer_atoms *atoms, double *alat1,
                                              double *alat2, double *alat3);
void FC_FUNC_(atoms_new_from_file, ATOMS_NEW_FROM_FILE)(int *lstat, f90_pointer_atoms *atoms,
                                                        f90_pointer_double *rxyz,
                                                        const gchar *filename, int *ln);
void FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(int *iproc, f90_pointer_atoms *atoms,
                                                      int *ixc);
void FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(f90_pointer_atoms *atoms,
                                                          double *radii_cf);
void FC_FUNC_(atoms_set_symmetries, ATOMS_SET_SYMMETRIES)(void *atoms, double *rxyz,
                                                          int *disable, double *elecfield);
void FC_FUNC_(atoms_set_displacement, ATOMS_SET_DISPLACEMENT)(void *atoms, double *rxyz,
                                                              double *randdis);


#define GET_ATTR_INT(name,NAME) { \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int)); \
  FC_FUNC_(atoms_get_ ## name, ATOMS_GET_ ## NAME)(atoms->data->atoms, &tmp); \
  atoms->name = (guint*)tmp.data; \
  }

BigDFT_Atoms* bigdft_atoms_init()
{
  BigDFT_Atoms *atoms;

  atoms = g_malloc(sizeof(BigDFT_Atoms));
  memset(atoms, 0, sizeof(BigDFT_Atoms));
  atoms->data = g_malloc(sizeof(f90_pointer_atoms));
  memset(atoms->data, 0, sizeof(f90_pointer_atoms));
  memset(&atoms->rxyz, 0, sizeof(f90_pointer_double));

  return atoms;
}
void bigdft_atoms_dispose(BigDFT_Atoms *atoms)
{
  guint i;

  g_free(atoms->data);
  FC_FUNC_(deallocate_double, DEALLOCATE_DOUBLE)(&atoms->rxyz);
  if (atoms->atomnames)
    {
      for (i = 0; i < atoms->ntypes; i++)
        if (atoms->atomnames[i])
          g_free(atoms->atomnames[i]);
      g_free(atoms->atomnames);
    }
  g_free(atoms);
}

BigDFT_Atoms* bigdft_atoms_new()
{
  BigDFT_Atoms *atoms;

  atoms = bigdft_atoms_init();
  FC_FUNC_(atoms_new, ATOMS_NEW)(atoms->data);

  return atoms;
}
void bigdft_atoms_free(BigDFT_Atoms *atoms)
{
  FC_FUNC_(atoms_free, ATOMS_FREE)(atoms->data);
  bigdft_atoms_dispose(atoms);
}
BigDFT_Atoms* bigdft_atoms_new_from_file(const gchar *filename)
{
  BigDFT_Atoms *atoms;
  guint lstat, ln, i, j;
  gchar str[20];

  atoms = bigdft_atoms_init();
  ln = strlen(filename);
  FC_FUNC_(atoms_new_from_file, ATOMS_NEW_FROM_FILE)((int*)(&lstat), atoms->data,
                                                     &atoms->rxyz, filename, (int*)(&ln));

  if (!lstat)
    {
      bigdft_atoms_dispose(atoms);
      atoms = (BigDFT_Atoms*)0;
    }
  else
    {
      FC_FUNC_(atoms_get_geocode, ATOMS_GET_GEOCODE)(atoms->data->atoms, (gchar*)(&atoms->geocode));
      FC_FUNC_(atoms_get_nat, ATOMS_GET_NAT)(atoms->data->atoms, (int*)(&atoms->nat));
      FC_FUNC_(atoms_get_ntypes, ATOMS_GET_NTYPES)(atoms->data->atoms, (int*)(&atoms->ntypes));
      atoms->atomnames = g_malloc(sizeof(gchar*) * atoms->ntypes);
      for (i = 0; i < atoms->ntypes; i++)
        {
          j = i + 1;
          FC_FUNC_(atoms_get_name, ATOMS_GET_NAME)(atoms->data->atoms, (int*)(&j), str, (int*)(&ln));
          atoms->atomnames[i] = g_malloc(sizeof(gchar) * (ln + 1));
          memcpy(atoms->atomnames[i], str, sizeof(gchar) * ln);
          atoms->atomnames[i][ln] = '\0';
        }
      FC_FUNC_(atoms_get_alat, ATOMS_GET_ALAT)(atoms->data->atoms, atoms->alat, atoms->alat + 1, atoms->alat + 2);
      GET_ATTR_INT(iatype,IATYPE);
    }

  return atoms;
}
void bigdft_atoms_set_n_atoms(BigDFT_Atoms *atoms, guint nat)
{
  FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(atoms->data->atoms, (int*)(&nat));
  atoms->nat = nat;

  GET_ATTR_INT(iatype,IATYPE);
}
void bigdft_atoms_set_n_types(BigDFT_Atoms *atoms, guint ntypes)
{
  FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(atoms->data, (int*)(&ntypes));
  atoms->ntypes = ntypes;
  atoms->atomnames = g_malloc(sizeof(gchar*) * ntypes);
  memset(atoms->atomnames, 0, sizeof(gchar*) * ntypes);
}

void bigdft_atoms_set_psp(BigDFT_Atoms *atoms, int ixc)
{
  int iproc = 0;
  FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(&iproc, atoms->data->atoms, &ixc);
}

void bigdft_atoms_set_symmetries(BigDFT_Atoms *atoms, gboolean active, double elecfield[3])
{
  int disable;

  disable = (!active);
  FC_FUNC_(atoms_set_symmetries, ATOMS_SET_SYMMETRIES)(atoms->data->atoms, atoms->rxyz.data,
                                                       &disable, elecfield);
}

void bigdft_atoms_set_displacement(BigDFT_Atoms *atoms, double randdis)
{
  FC_FUNC_(atoms_set_displacement, ATOMS_SET_DISPLACEMENT)(atoms->data->atoms, atoms->rxyz.data, &randdis);
}

double* bigdft_atoms_get_radii(const BigDFT_Atoms *atoms)
{
  double *radii_cf;

  radii_cf = g_malloc(sizeof(double) * 3 * atoms->ntypes);
  FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(atoms->data->atoms, radii_cf);
  return radii_cf;
}


/* Wavefunction descriptor part. */
struct f90_pointer_glr_
{
  void *glr;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(glr_new, GLR_NEW)(f90_pointer_glr *glr);
void FC_FUNC_(system_size, SYSTEM_SIZE)(int *iproc, f90_pointer_atoms *atoms, double *rxyz,
                                        double *radii_cf, double *crmult, double *frmult,
                                        double *hx, double *hy, double *hz,
                                        f90_pointer_glr *glr, double *shift);
void FC_FUNC_(glr_get_n, GLR_GET_N)(f90_pointer_glr *glr, int *n);
void FC_FUNC_(glr_free, GLR_FREE)(f90_pointer_glr *glr);
void FC_FUNC(glr_set_wave_descriptors,
             GLR_SET_WAVE_DESCRIPTORS)(int *iproc, double *hx, double *hy,
                                       double *hz, void *atoms, double *rxyz, double *radii,
                                       double *crmult, double *frmult, void *glr);

BigDFT_Glr* bigdft_glr_init(BigDFT_Atoms *atoms, double *radii, double h[3],
                            double crmult, double frmult)
{
  BigDFT_Glr *glr;
  int iproc = 1;

  glr = g_malloc(sizeof(BigDFT_Glr));
  memset(glr, 0, sizeof(BigDFT_Glr));
  glr->data = g_malloc(sizeof(f90_pointer_glr));
  memset(glr->data, 0, sizeof(f90_pointer_glr));

  FC_FUNC_(glr_new, GLR_NEW)(glr->data);
  FC_FUNC_(system_size, SYSTEM_SIZE)(&iproc, atoms->data->atoms, atoms->rxyz.data, radii,
                                     &crmult, &frmult, h, h + 1, h + 2, glr->data->glr,
                                     atoms->shift);
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  FC_FUNC_(glr_get_n, GLR_GET_N)(glr->data->glr, (int*)glr->n);
  FC_FUNC_(atoms_get_alat, ATOMS_GET_ALAT)(atoms->data->atoms, atoms->alat, atoms->alat + 1, atoms->alat + 2);

  return glr;
}
void bigdft_glr_dispose(BigDFT_Glr *glr)
{
  g_free(glr->data);
  g_free(glr);
}

BigDFT_Glr* bigdft_glr_new(BigDFT_Atoms *atoms, double *radii, double h[3],
                           double crmult, double frmult)
{
  BigDFT_Glr *glr;

  glr = bigdft_glr_init(atoms, radii, h, crmult, frmult);
  
  return glr;
}
BigDFT_Glr* bigdft_glr_new_with_wave_descriptors(BigDFT_Atoms *atoms, double *radii,
                                                 double h[3], double crmult, double frmult)
{
  BigDFT_Glr *glr;

  glr = bigdft_glr_init(atoms, radii, h, crmult, frmult);
  bigdft_glr_set_wave_descriptors(glr, atoms, radii, crmult, frmult);
  
  return glr;
}
void bigdft_glr_free(BigDFT_Glr *glr)
{
  FC_FUNC_(glr_free, GLR_FREE)(glr->data);
  bigdft_glr_dispose(glr);
}
void bigdft_glr_set_wave_descriptors(BigDFT_Glr *glr, BigDFT_Atoms *atoms, double *radii,
                                     double crmult, double frmult)
{
  int iproc = 1;

  FC_FUNC(glr_set_wave_descriptors,
          GLR_SET_WAVE_DESCRIPTORS)(&iproc, glr->h, glr->h + 1, glr->h + 2,
                                    atoms->data->atoms, atoms->rxyz.data, radii,
                                    &crmult, &frmult, glr->data->glr);
}

void FC_FUNC_(fill_logrid, FILL_LOGRID)(char *geocode, int *n1, int *n2, int *n3,
                                        int *nl1, int *nu1, int *nl2, int *nu2, int *nl3, int *nu3,
                                        int *orig, int *nat, int *ntypes, int *iatype,
                                        double *rxyz, double *radii, double *mult,
                                        double *hx, double *hy, double *hz, int *grid);
guint* bigdft_fill_logrid(BigDFT_Atoms *atoms, guint n[3], double *radii,
                          double mult, double h[3])
{
  guint *grid;
  int orig = 0;

  grid = g_malloc(sizeof(guint) * (n[0] + 1) * (n[1] + 1) * (n[2] + 1));

  FC_FUNC_(fill_logrid, FILL_LOGRID)(&atoms->geocode, (int*)n, (int*)(n + 1), (int*)(n + 2),
                                     &orig, (int*)n, &orig, (int*)(n + 1), &orig, (int*)(n + 2),
                                     &orig, (int*)(&atoms->nat), (int*)(&atoms->ntypes), (int*)atoms->iatype,
                                     atoms->rxyz.data, radii, &mult, h, h + 1, h + 2, (int*)grid);

  return grid;
}


/* Wavefunction descriptor part. */
struct f90_pointer_inputs_
{
  void *in;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(inputs_new, INPUTS_NEW)(f90_pointer_inputs *in);
void FC_FUNC_(inputs_free, INPUTS_FREE)(f90_pointer_inputs *in);
void FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(f90_pointer_inputs *in,
                                                      const gchar *rad, int *len);
void FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(const f90_pointer_inputs *in, double *hx, double *hy, double *hz, double *crmult, double *frmult, int *ixc, int *ncharge, double *elecfield, int *nspin, int *mpol, double *gnrm_cv, int *itermax, int *nrepmax, int *ncong, int *idsx, int *dispersion, int *inputPsiId, int *output_wf_format, int *output_grid, double *rbuf, int *ncongt, int *norbv, int *nvirt, int *nplot, int *disableSym);
void FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(void *in, int *iscf, int *itrpmax,
                                              int *norbsempty, int *occopt, double *alphamix,
                                              double *rpnrm_cv, double *gnrm_startmix,
                                              double *Tel, double *alphadiis);
void FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(void *in, char *geopt_approach,
                                                  int *ncount_cluster_x, double *frac_fluct,
                                                  double *forcemax, double *randdis,
                                                  double *betax, int *history, int *ionmov,
                                                  double *dtion, double *strtarget,
                                                  f90_pointer_double *qmass);
void FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(f90_pointer_inputs *in,
                                                        int *iproc, int *dump);
void FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(const f90_pointer_inputs *in, int *files);
void FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(void *in, void *atoms,
                                                  int *iproc, int *dump);

BigDFT_Inputs* bigdft_inputs_init()
{
  BigDFT_Inputs *in;

  in = g_malloc(sizeof(BigDFT_Inputs));
  memset(in, 0, sizeof(BigDFT_Inputs));
  in->data = g_malloc(sizeof(f90_pointer_inputs));
  memset(in->data, 0, sizeof(f90_pointer_inputs));
  memset(&in->qmass, 0, sizeof(f90_pointer_double));

  return in;
}
void bigdft_inputs_dispose(BigDFT_Inputs *in)
{
  g_free(in->data);
  g_free(in);
}

BigDFT_Inputs* bigdft_inputs_new(const gchar *naming)
{
  BigDFT_Inputs *in;
  int iproc = 0, len, dump = 0;

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data->in, naming, &len);
    }
  else
    {
      len = 1;
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data->in, " ", &len);
    }
  FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(in->data->in, &iproc, &dump);

  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data->in, &in->files);
  FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(in->data->in, in->h, in->h + 1, in->h + 2,
                                           &in->crmult, &in->frmult, &in->ixc,
                                           &in->ncharge, in->elecfield, &in->nspin,
                                           &in->mpol, &in->gnrm_cv, &in->itermax,
                                           &in->nrepmax, &in->ncong, &in->idsx,
                                           &in->dispersion, &in->inputPsiId,
                                           &in->output_wf_format, &in->output_grid,
                                           &in->rbuf, &in->ncongt, &in->norbv, &in->nvirt,
                                           &in->nplot, &in->disableSym);
  FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(in->data->in, &in->iscf, &in->itrpmax,
                                           &in->norbsempty, (int*)(&in->occopt), &in->alphamix,
                                           &in->rpnrm_cv, &in->gnrm_startmix, &in->Tel,
                                           &in->alphadiis);
  FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(in->data->in, in->geopt_approach,
                                               &in->ncount_cluster_x, &in->frac_fluct,
                                               &in->forcemax, &in->randdis, &in->betax,
                                               &in->history, &in->ionmov, &in->dtion,
                                               in->strtarget, &in->qmass);
  /* FC_FUNC_(inputs_get_sic, INPUTS_GET_SIC)(); */
  /* FC_FUNC_(inputs_get_tddft, INPUTS_GET_TDDFT)(); */
  
  return in;
}
void bigdft_inputs_free(BigDFT_Inputs *in)
{
  FC_FUNC_(inputs_free, INPUTS_FREE)(in->data);
  bigdft_inputs_dispose(in);
}
void bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms)
{
  int iproc = 0, dump = 0;

  FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(in->data->in, atoms->data->atoms, &iproc, &dump);
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data->in, &in->files);

  /* FC_FUNC_(inputs_get_kpt, INPUTS_GET_KPT)(); */
}

/* Orbital descriptor part. */
struct f90_pointer_orbs_
{
  void *orbs;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(orbs_new, ORBS_NEW)(f90_pointer_orbs *orbs);
void FC_FUNC_(orbs_free, ORBS_FREE)(f90_pointer_orbs *orbs);
void FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(int *iproc, int *nproc,
                                                              void *in, void *atoms,
                                                              void *orbs, int *nelec);
void FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(void *orbs, int *norb,
                                                        int *norbp, int *norbu,
                                                        int *norbd, int *nspin,
                                                        int *nspinor, int *npsidim,
                                                        int *nkpts, int *nkptsp,
                                                        int *isorb, int *iskpts);

BigDFT_Orbs* bigdft_orbs_init()
{
  BigDFT_Orbs *orbs;

  orbs = g_malloc(sizeof(BigDFT_Orbs));
  memset(orbs, 0, sizeof(BigDFT_Orbs));
  orbs->data = g_malloc(sizeof(f90_pointer_orbs));
  memset(orbs->data, 0, sizeof(f90_pointer_orbs));

  return orbs;
}
void bigdft_orbs_dispose(BigDFT_Orbs *orbs)
{
  g_free(orbs->data);
  g_free(orbs);
}

BigDFT_Orbs* bigdft_orbs_new(const BigDFT_Atoms *atoms, const BigDFT_Inputs *in,
                             int iproc, int nproc, guint *nelec)
{
  BigDFT_Orbs *orbs;
  int nelec_;

  orbs = bigdft_orbs_init();
  FC_FUNC_(orbs_new, ORBS_NEW)(orbs->data);
  FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(&iproc, &nproc, in->data->in,
                                                           atoms->data->atoms,
                                                           orbs->data->orbs, &nelec_);
  if (nelec)
    *nelec = nelec_;
  
  FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(orbs->data->orbs, &orbs->norb,
                                                     &orbs->norbp, &orbs->norbu,
                                                     &orbs->norbd, &orbs->nspin,
                                                     &orbs->nspinor, &orbs->npsidim,
                                                     &orbs->nkpts, &orbs->nkptsp,
                                                     &orbs->isorb, &orbs->iskpts);
  
  return orbs;
}
void bigdft_orbs_free(BigDFT_Orbs *orbs)
{
  FC_FUNC_(orbs_free, ORBS_FREE)(orbs->data);
  bigdft_orbs_dispose(orbs);
}

void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(int *nproc, int *idsx, void *lr, int *nat,
                                               int *norb, int *nspinor, int *nkpt,
                                               int *nprojel, int *nspin, int *itrpmax,
                                               int *iscf, double *peak);
double bigdft_memory_peak(int nproc, BigDFT_Glr *lr, BigDFT_Inputs *in, BigDFT_Orbs *orbs)
{
  double peak;
  int nat = -1;

  FC_FUNC(memoryestimator, MEMORYESTIMATOR)(&nproc, &in->idsx, lr->data->glr, &nat,
                                            &orbs->norb, &orbs->nspinor, &orbs->nkpts,
                                            &nat, &orbs->nspin, &in->itrpmax, &in->iscf,
                                            &peak);
  return peak;
}
