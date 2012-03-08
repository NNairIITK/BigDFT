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

typedef struct BigDFT_lzd_ BigDFT_Lzd;
typedef struct BigDFT_orbs_ BigDFT_Orbs;
typedef struct BigDFT_proj_ BigDFT_Proj;
typedef struct BigDFT_LocalFields_ BigDFT_LocalFields;

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
  gboolean dispose_has_run;
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
  void *sym;
} BigDFT_Atoms;


BigDFT_Atoms* bigdft_atoms_new();
BigDFT_Atoms* bigdft_atoms_new_from_file   (const gchar *filename);
void          bigdft_atoms_free            (BigDFT_Atoms *atoms);
void          bigdft_atoms_set_n_atoms     (BigDFT_Atoms *atoms, guint nat);
void          bigdft_atoms_set_n_types     (BigDFT_Atoms *atoms, guint ntypes);
gboolean      bigdft_atoms_set_structure_from_file(BigDFT_Atoms *atoms, const gchar *filename);
void          bigdft_atoms_set_psp         (BigDFT_Atoms *atoms, int ixc,
                                            guint nspin, const gchar *occup);
void          bigdft_atoms_set_symmetries  (BigDFT_Atoms *atoms, gboolean active,
                                            double tol, double elecfield[3]);
void          bigdft_atoms_set_displacement(BigDFT_Atoms *atoms, double randdis);
void          bigdft_atoms_sync            (BigDFT_Atoms *atoms);
double*       bigdft_atoms_get_radii       (const BigDFT_Atoms *atoms, double crmult,
                                            double frmult, double projrad);

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

/*********************************/
/* BigDFT_LocReg data structure. */
/*********************************/
#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_LOCREG_TYPE    (bigdft_locreg_get_type())
#define BIGDFT_LOCREG(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_LOCREG_TYPE, BigDFT_LocReg))
typedef struct BigDFT_LocRegClass_
{
  GObjectClass parent;
} BigDFT_LocRegClass;
#else
#define BIGDFT_LOCREG_TYPE    (999)
#define BIGDFT_LOCREG(obj)    ((BigDFT_LocReg*)obj)
#endif
typedef struct BigDFT_locReg_
{
#ifdef GLIB_MAJOR_VERSION
  BigDFT_Atoms parent;
  gboolean dispose_has_run;
#endif

  double h[3];
  guint n[3], ni[3];

  /* Values that have been used to built this localisation region. */
  double *radii;
  double crmult, frmult;
  
  /* TODO: bindings to values... */

  /* Private. */
  void *d;
  void *data;
} BigDFT_LocReg;
typedef enum
  {
    GRID_COARSE,
    GRID_FINE
  } BigDFT_Grid;

BigDFT_LocReg* bigdft_locreg_new                 ();
void           bigdft_locreg_free                (BigDFT_LocReg *glr);
void           bigdft_locreg_set_radii           (BigDFT_LocReg *glr, const double *radii);
void           bigdft_locreg_set_size            (BigDFT_LocReg *glr, double h[3],
                                                  double crmult, double frmult);
void           bigdft_locreg_set_wave_descriptors(BigDFT_LocReg *glr);
gboolean*      bigdft_locreg_get_grid            (const BigDFT_LocReg *glr,
                                                  BigDFT_Grid gridType);
double*        bigdft_locreg_convert_to_isf      (const BigDFT_LocReg *glr,
                                                  const double *psic);

/*********************************/
/* BigDFT_Lzd data structure. */
/*********************************/
#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_LZD_TYPE    (bigdft_lzd_get_type())
#define BIGDFT_LZD(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_LZD_TYPE, BigDFT_Lzd))
typedef struct BigDFT_LzdClass_
{
  GObjectClass parent;
} BigDFT_LzdClass;
#else
#define BIGDFT_LZD_TYPE    (999)
#define BIGDFT_LZD(obj)    ((BigDFT_Lzd*)obj)
#endif
struct BigDFT_lzd_
{
  BigDFT_LocReg parent;
#ifdef GLIB_MAJOR_VERSION
  gboolean dispose_has_run;
#endif

  /* Private. */
  void *data;
};
BigDFT_Lzd* bigdft_lzd_new ();
void        bigdft_lzd_setup_linear(BigDFT_Lzd *lzd, BigDFT_Orbs *orbs,
                                    const BigDFT_Inputs *in, guint iproc, guint nproc);
void        bigdft_lzd_free(BigDFT_Lzd *lzd);


/*******************************/
/* BigDFT_Orbs data structure. */
/*******************************/
#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_ORBS_TYPE    (bigdft_orbs_get_type())
#define BIGDFT_ORBS(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_ORBS_TYPE, BigDFT_Orbs))
typedef struct BigDFT_OrbsClass_
{
  GObjectClass parent;
} BigDFT_OrbsClass;
#else
#define BIGDFT_ORBS_TYPE    (999)
#define BIGDFT_ORBS(obj)    ((BigDFT_Orbs*)obj)
#endif
struct BigDFT_orbs_
{
#ifdef GLIB_MAJOR_VERSION
  GObject parent;
  gboolean dispose_has_run;
#endif

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

  /* Pointers on building objects. */
  const BigDFT_Inputs *in;
  const BigDFT_LocReg *glr;

  /* Private. */
  void *data;
  void *comm;
};

BigDFT_Orbs* bigdft_orbs_new (const BigDFT_LocReg *glr, const BigDFT_Inputs *in,
                              guint iproc, guint nproc, guint *nelec);
void         bigdft_orbs_free(BigDFT_Orbs *orbs);

/*****************************/
/* BigDFT_Wf data structure. */
/*****************************/
#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_WF_TYPE    (bigdft_wf_get_type())
#define BIGDFT_WF(obj)                                          \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_WF_TYPE, BigDFT_Wf))
typedef struct BigDFT_WfClass_
{
  GObjectClass parent;
} BigDFT_WfClass;
#else
#define BIGDFT_WF_TYPE    (999)
#define BIGDFT_WF(obj)    ((BigDFT_Wf*)obj)
#endif
typedef struct BigDFT_wf_
{
  BigDFT_Orbs parent;
#ifdef GLIB_MAJOR_VERSION
  gboolean dispose_has_run;
#endif

  /* Pointers on building objects. */
  const BigDFT_Lzd *lzd;

  /* Private. */
  f90_pointer_double psi, hpsi, psit;
} BigDFT_Wf;
typedef enum
  {
    BIGDFT_SPIN_UP,
    BIGDFT_SPIN_DOWN
  } BigDFT_Spin;
typedef enum
  {
    BIGDFT_REAL,
    BIGDFT_IMAG,
    BIGDFT_PARTIAL_DENSITY
  } BigDFT_Spinor;

BigDFT_Wf* bigdft_wf_new (BigDFT_Lzd *lzd, BigDFT_Inputs *in,
                          guint iproc, guint nproc, guint *nelec);
void       bigdft_wf_free(BigDFT_Wf *wf);
void       bigdft_wf_calculate_psi0(BigDFT_Wf *wf, BigDFT_LocalFields *denspot,
                                    BigDFT_Proj *proj, guint iproc, guint nproc);
const double* bigdft_wf_get_psi_compress(const BigDFT_Wf *wf, guint ikpt, guint iorb,
                                         BigDFT_Spin ispin, BigDFT_Spinor ispinor,
                                         guint *psiSize, guint iproc);
double*    bigdft_wf_convert_to_isf(const BigDFT_Wf *wf, guint ikpt, guint iorb,
                                    BigDFT_Spin ispin, BigDFT_Spinor ispinor, guint iproc);

/*******************************/
/* BigDFT_Proj data structure. */
/*******************************/
struct BigDFT_proj_
{
  /* TODO: bindings to values... */
  guint nproj, nprojel;

  /* Additional pointers. */
  f90_pointer_double proj;

  /* Private. */
  void *nlpspd;
};

BigDFT_Proj* bigdft_proj_new (const BigDFT_LocReg *glr, const BigDFT_Orbs *orbs,
                              double frmult);
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
struct BigDFT_LocalFields_
{
#ifdef GLIB_MAJOR_VERSION
  GObject parent;
  gboolean dispose_has_run;
#endif
  /* bindings to values... */
  BigDFT_RhoIs rhov_is;
  double psoffset;
  double h[3];

  /* Additional pointers. */
  double *rhov, *v_ext, *v_xc;
  /* TODO, see when these are associated. */
  /* double *rho_full, *pot_full, *rho_psi, *rho_c, *vloc_ks, *f_xc; */

  /* Pointers on building objects. */
  const BigDFT_LocReg *glr;

  /* Private. */
  double *pkernel, *pkernelseq;
  void *rhod;
  void *dpcom;
  void *data;
};

BigDFT_LocalFields* bigdft_localfields_new (const BigDFT_LocReg *glr,
                                            const BigDFT_Inputs *in,
                                            guint iproc, guint nproc);
void                bigdft_localfields_free(BigDFT_LocalFields *denspotd);
void bigdft_localfields_create_effective_ionic_pot(BigDFT_LocalFields *denspot,
                                                   const BigDFT_Inputs *in,
                                                   guint iproc, guint nproc);

/******************/
/* Miscellaneous. */
/******************/
double bigdft_memory_get_peak(int nproc, const BigDFT_LocReg *lr, const BigDFT_Inputs *in,
                              const BigDFT_Orbs *orbs, const BigDFT_Proj *proj);

f90_pointer_double_4D* bigdft_read_wave_to_isf(const gchar *filename, int iorbp,
                                               double h[3], int n[3], int *nspinor);
void bigdft_free_wave_to_isf(f90_pointer_double_4D *psiscf);
gboolean bigdft_read_wave_descr(const gchar *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor);


#endif
