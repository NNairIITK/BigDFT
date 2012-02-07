/*> @file 
     Header for the public BigDFT API.
    @author
     Copyright (C) 2011-2012 BigDFT group 
     This file is distributed under the terms of the
     GNU General Public License, see ~/COPYING file
     or http://www.gnu.org/copyleft/gpl.txt .
     For the list of contributors, see ~/AUTHORS 
*/

#ifndef BIGDFT_H
#define BIGDFT_H

#include <bigdft_cst.h>

/********************************/
/* BigDFT_Atoms data structure. */
/********************************/
#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_ATOMS_TYPE    (bigdft_atoms_get_type())
#define BIGDFT_ATOMS(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_ATOMS_TYPE, BigDFT_Atoms))
typedef struct BigDFT_AtomsClass_
{
  GObjectClass parent;
} BigDFT_AtomsClass;
#else
#define BIGDFT_ATOMS_TYPE    (999)
#define BIGDFT_ATOMS(obj)    ((BigDFT_Atoms*)obj)
#endif
typedef struct BigDFT_Atoms_
{
#ifdef GLIB_MAJOR_VERSION
  GObject parent;
#endif
  /* Bindings to values, obtained by copy. Update them with
     bigdft_atoms_sync(). */
  gchar geocode, format[6], units[21];
  guint  ntypes, nat, natsc;
  double alat[3];
  gchar **atomnames;
  gboolean donlcc;
  /* Bindings of pointer arrays, access is direct. */
  guint *iatype, *iasctype, *natpol, *nelpsp, *npspcode, *nzatom;
  int *nlcc_ngv, *nlcc_ngc, *ixcpsp, *ifrztyp;
  double *radii_cf, *amu, *aocc, *psppar, *nlccpar, *ig_nlccpar;

  /* Coordinates. */
  f90_pointer_double_2D rxyz;
  double shift[3];

  /* Additional fields. */
  gchar *comment;
  double energy;

  /* Private. */
  void *data;

  /* GObject machinery. */
  gboolean dispose_has_run;
} BigDFT_Atoms;


BigDFT_Atoms* bigdft_atoms_new();
BigDFT_Atoms* bigdft_atoms_new_from_file   (const gchar *filename);
void          bigdft_atoms_free            (BigDFT_Atoms *atoms);
void          bigdft_atoms_set_n_atoms     (BigDFT_Atoms *atoms, guint nat);
void          bigdft_atoms_set_n_types     (BigDFT_Atoms *atoms, guint ntypes);
void          bigdft_atoms_set_psp         (BigDFT_Atoms *atoms, int ixc);
void          bigdft_atoms_set_symmetries  (BigDFT_Atoms *atoms, gboolean active,
                                            double elecfield[3]);
void          bigdft_atoms_set_displacement(BigDFT_Atoms *atoms, double randdis);
void          bigdft_atoms_sync            (BigDFT_Atoms *atoms);
double*       bigdft_atoms_get_radii       (const BigDFT_Atoms *atoms);

/*********************************/
/* BigDFT_Inputs data structure. */
/*********************************/
typedef enum
  {
    SMEARING_DIST_ERF   = 1,
    SMEARING_DIST_FERMI = 2,
    SMEARING_DIST_COLD1 = 3,
    SMEARING_DIST_COLD2 = 4,
    SMEARING_DIST_METPX = 5
  } BigDFT_Smearing;

typedef struct BigDFT_Inputs_
{
  /* TODO: bindings to values... */
  int files;
  
  /* DFT file variables. */
  int ixc, ncharge, nspin, mpol, itermax, nrepmax, ncong, idsx,
    dispersion, inputPsiId, output_wf_format, output_grid, ncongt, norbv, nvirt,
    nplot, disableSym;
  double crmult, frmult, gnrm_cv, rbuf;
  double h[3], elecfield[3];

  /* MIX file variables. */
  int iscf, itrpmax, norbsempty;
  BigDFT_Smearing occopt;
  double alphamix, rpnrm_cv, gnrm_startmix, Tel, alphadiis;

  /* GEOPT file variables. */
  char geopt_approach[10];
  int ncount_cluster_x, history, ionmov;
  double frac_fluct, forcemax, randdis, betax, dtion;
  double strtarget[6];
  f90_pointer_double qmass;

  /* Private. */
  void *data;
} BigDFT_Inputs;

#define BIGDFT_INPUTS_NONE    0
#define BIGDFT_INPUTS_DFT     1
#define BIGDFT_INPUTS_GEOPT   2
#define BIGDFT_INPUTS_PERF    4
#define BIGDFT_INPUTS_KPT     8
#define BIGDFT_INPUTS_MIX    16
#define BIGDFT_INPUTS_TDDFT  32
#define BIGDFT_INPUTS_SIC    64
#define BIGDFT_INPUTS_FREQ  128

BigDFT_Inputs* bigdft_inputs_new             (const gchar *radical);
void           bigdft_inputs_free            (BigDFT_Inputs *in);
void           bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms);

/******************************/
/* BigDFT_Glr data structure. */
/******************************/
typedef struct BigDFT_glr_
{
  gchar geocode;
  double h[3];
  guint n[3], ni[3];
  
  /* TODO: bindings to values... */

  /* Private. */
  void *d;
  void *data;
} BigDFT_Glr;

BigDFT_Glr* bigdft_glr_new                 (BigDFT_Atoms *atoms, double *radii, double h[3],
                                            double crmult, double frmult);
BigDFT_Glr* bigdft_glr_new_with_wave_descriptors(BigDFT_Atoms *atoms, double *radii,
                                                 double h[3], double crmult, double frmult);
void        bigdft_glr_free                (BigDFT_Glr *glr);
void        bigdft_glr_set_wave_descriptors(BigDFT_Glr *glr, BigDFT_Atoms *atoms,
                                            double *radii, double crmult, double frmult);

/*******************************/
/* BigDFT_Orbs data structure. */
/*******************************/
typedef struct f90_pointer_orbs_ f90_pointer_orbs;
typedef struct BigDFT_orbs_
{
  /* TODO: bindings to values... */
  int norb, norbp, norbu, norbd;
  int nspin, nspinor, npsidim;
  int nkpts, nkptsp;
  int isorb, iskpts;

  double efermi, HLgap, eTS;

     /* integer, dimension(:), pointer :: iokpt,ikptproc */
     /* integer, dimension(:,:), pointer :: norb_par */
     /* real(wp), dimension(:), pointer :: eval */
     /* real(gp), dimension(:), pointer :: occup,spinsgn,kwgts */
     /* real(gp), dimension(:,:), pointer :: kpts */

  /* Private. */
  f90_pointer_orbs *data;
} BigDFT_Orbs;

BigDFT_Orbs* bigdft_orbs_new (const BigDFT_Atoms *atoms, const BigDFT_Inputs *in,
                              const BigDFT_Glr *glr, int iproc, int nproc, guint *nelec);
void         bigdft_orbs_free(BigDFT_Orbs *orbs);

/*******************************/
/* BigDFT_Proj data structure. */
/*******************************/
typedef struct f90_pointer_nlpspd_ f90_pointer_nlpspd;
typedef struct BigDFT_proj_
{
  /* TODO: bindings to values... */
  guint nproj, nprojel;

  /* Additional pointers. */
  f90_pointer_double proj;

  /* Private. */
  f90_pointer_nlpspd *nlpspd;
} BigDFT_Proj;

BigDFT_Proj* bigdft_proj_new (const BigDFT_Atoms *atoms, const BigDFT_Glr *glr,
                              const BigDFT_Orbs *orbs, double *radii, double frmult);
void         bigdft_proj_free(BigDFT_Proj *proj);

/**********************************/
/* BigDFT_DensPot data structure. */
/**********************************/
typedef enum
  {
    BIGDFT_RHO_IS_EMPTY              = -1980,
    BIGDFT_RHO_IS_ELECTRONIC_DENSITY = -1979,
    BIGDFT_RHO_IS_CHARGE_DENSITY     = -1978,
    BIGDFT_RHO_IS_KS_POTENTIAL       = -1977,
    BIGDFT_RHO_IS_HARTREE_POTENTIAL  = -1976
  } BigDFT_RhoIs;

#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_LOCALFIELDS_TYPE    (bigdft_localfields_get_type())
#define BIGDFT_LOCALFIELDS(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_LOCALFIELDS_TYPE, BigDFT_LocalFields))
typedef struct BigDFT_LocalFieldsClass_
{
  GObjectClass parent;
} BigDFT_LocalFieldsClass;
#else
#define BIGDFT_LOCALFIELDS_TYPE    (999)
#define BIGDFT_LOCALFIELDS(obj)    ((BigDFT_LocalFields*)obj)
#endif
typedef struct BigDFT_LocalFields_
{
#ifdef GLIB_MAJOR_VERSION
  GObject parent;
#endif
  /* bindings to values... */
  BigDFT_RhoIs rhov_is;
  double psoffset;
  double h[3];

  /* Additional pointers. */
  f90_pointer_double rhov, rho_full, pot_full;
  f90_pointer_double_2D rho_psi;
  f90_pointer_double_4D rho_c, v_ext, v_xc, vloc_ks, f_xc;
  f90_pointer_double pkernel, pkernelseq;

  /* Private. */
  void *rhod;
  void *dpcom;
  void *data;

  /* GObject machinery. */
  gboolean dispose_has_run;
} BigDFT_LocalFields;

BigDFT_LocalFields* bigdft_localfields_new (const BigDFT_Atoms *atoms,
                                            const BigDFT_Glr *glr,
                                            const BigDFT_Inputs *in,
                                            const double *radii,
                                            guint iproc, guint nproc);
void                bigdft_localfields_free(BigDFT_LocalFields *denspotd);

/******************/
/* Miscellaneous. */
/******************/
double bigdft_memory_peak(int nproc, BigDFT_Glr *lr, BigDFT_Inputs *in,
                          BigDFT_Orbs *orbs, BigDFT_Proj *proj);
guint* bigdft_fill_logrid(BigDFT_Atoms *atoms, guint n[3], double *radii,
                          double mult, double h[3]);
f90_pointer_double_4D* bigdft_read_wave_to_isf(const gchar *filename, int iorbp,
                                            double h[3], int n[3], int *nspinor);
void bigdft_free_wave_to_isf(f90_pointer_double_4D *psiscf);
gboolean bigdft_read_wave_descr(const gchar *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor);


#endif
