#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#endif

#include <bigdft.h>

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

/* #undef G_THREADS_ENABLED */

#define MAX_FORTRAN_OUTPUT 4096

typedef struct bigdft_data
{
  BigDFT_Inputs      *in;
  BigDFT_Proj        *proj;
  BigDFT_LocalFields *denspot;
  BigDFT_Wf          *wf;
  BigDFT_Energs      *energs;

#ifdef HAVE_GLIB
  GMainLoop *loop;
#endif
} BigDFT_Data;


#ifdef HAVE_GLIB
static gboolean exit_loop(gpointer data);
static void onVExtReady(BigDFT_LocalFields *denspot, gpointer data);
static void onDensReady(BigDFT_LocalFields *denspot, guint istep, gpointer data);
static void onEksReady(BigDFT_Energs *energs, guint iter, gpointer data);
#endif
static int redirect_init(int out_pipe[2]);
static void redirect_dump(int out_pipe[2], int stdout_fileno_old);
static BigDFT_Data* run_bigdft(BigDFT_Inputs *in, BigDFT_Proj *proj,
                               BigDFT_LocalFields *denspot, BigDFT_Wf *wf, gpointer data);

static void output_inputs(const BigDFT_Inputs *in)
{
  fprintf(stdout, " Read 'test.dft' file: %d\n", (in->files & BIGDFT_INPUTS_DFT));
  fprintf(stdout, " Read 'test.lin' file: %d\n", (in->files & BIGDFT_INPUTS_LIN));
  fprintf(stdout, " Input variables are %f %f %f  -  %f %f  -  %d\n",
          in->h[0], in->h[1], in->h[2], in->crmult, in->frmult, in->ixc);
}
static void output_atoms(const BigDFT_Atoms *atoms)
{
  guint i;

  for (i = 0; i  < atoms->nat; i++)
    fprintf(stdout, " Atom %d, coord. %f %f %f, type %d\n", i,
            atoms->rxyz.data[3 * i], atoms->rxyz.data[3 * i + 1],
            atoms->rxyz.data[3 * i + 2], atoms->iatype[i]);
  fprintf(stdout, " Box [%c], is in %f %f %f\n", atoms->geocode,
          atoms->alat[0], atoms->alat[1], atoms->alat[2]);
}
static void output_locreg(const BigDFT_LocReg *glr)
{
  guint i, n;
  gboolean *cgrid, *fgrid, valid;
  BigDFT_LocRegIter iter;

  for (i = 0; i < BIGDFT_ATOMS(glr)->ntypes; i++)
    fprintf(stdout, " Type %d, radii %f %f %f\n", i,
            glr->radii[i], glr->radii[BIGDFT_ATOMS(glr)->ntypes + i],
            glr->radii[BIGDFT_ATOMS(glr)->ntypes * 2 + i]);
  for (i = 0; i  < BIGDFT_ATOMS(glr)->nat; i++)
    fprintf(stdout, " Atoms %d, coord. %10.6f %10.6f %10.6f '%2s', type %d\n",
            i, BIGDFT_ATOMS(glr)->rxyz.data[3 * i],
            BIGDFT_ATOMS(glr)->rxyz.data[3 * i + 1],
            BIGDFT_ATOMS(glr)->rxyz.data[3 * i + 2],
            BIGDFT_ATOMS(glr)->atomnames[BIGDFT_ATOMS(glr)->iatype[i] - 1],
            BIGDFT_ATOMS(glr)->iatype[i]);
  fprintf(stdout, " Box is in   %f %f %f\n",
          BIGDFT_ATOMS(glr)->alat[0], BIGDFT_ATOMS(glr)->alat[1], BIGDFT_ATOMS(glr)->alat[2]);
  fprintf(stdout, " Shift is     %f %f  %f\n",
          BIGDFT_ATOMS(glr)->shift[0], BIGDFT_ATOMS(glr)->shift[1],
          BIGDFT_ATOMS(glr)->shift[2]);
  fprintf(stdout, " Geocode is  %c\n", BIGDFT_ATOMS(glr)->geocode);
  fprintf(stdout, " Grid is     %9d %9d %9d\n", glr->n[0],
          glr->n[1], glr->n[2]);
  fprintf(stdout, " Int grid is %9d %9d %9d\n", glr->ni[0],
          glr->ni[1], glr->ni[2]);
  fprintf(stdout, " H grids are %9.9g %9.9g %9.9g\n", glr->h[0],
          glr->h[1], glr->h[2]);
  fprintf(stdout, " Coarse grid has %7d points %5d segments.\n", glr->nvctr_c, glr->nseg_c);
  fprintf(stdout, " Fine grid has   %7d points %5d segments.\n", glr->nvctr_f, glr->nseg_f);

  for (valid = bigdft_locreg_iter_new(glr, &iter, GRID_FINE); valid && iter.iseg - glr->nseg_c < 10;
       valid = bigdft_locreg_iter_next(&iter))
    fprintf(stdout, " fine seg %3d has bounds (%2d - %2d;%2d;%2d)\n",
            iter.iseg - glr->nseg_c + 1, iter.i0, iter.i1, iter.i2, iter.i3);

  fprintf(stdout, "Test calculation of grid points.\n");
  cgrid = bigdft_locreg_get_grid(glr, GRID_COARSE);
  for (i = 0, n = 0; i < (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1); i++)
    if (cgrid[i])
      n += 1;
  fprintf(stdout, " Coarse grid has %7d points.\n", n);
  fgrid = bigdft_locreg_get_grid(glr, GRID_FINE);
  for (i = 0, n = 0; i < (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1); i++)
    if (fgrid[i])
      n += 1;
  fprintf(stdout, " Fine grid has   %7d points.\n", n);

  g_free(cgrid);
  g_free(fgrid);
}
static void output_llr(const BigDFT_Lzd *lzd)
{
  guint i;

  fprintf(stdout, " Lzd has %d local regions.\n", lzd->nlr);
  for (i = 0; i < lzd->nlr; i++)
    {
      fprintf(stdout, " region %d has (%4d/%4d) coarse/fine segments.\n",
              i, lzd->Llr[i]->nseg_c, lzd->Llr[i]->nseg_f);
      fprintf(stdout, " region %d has (%4d/%4d) coarse/fine elements.\n",
              i, lzd->Llr[i]->nvctr_c, lzd->Llr[i]->nvctr_f);
    }
}

static void init_atoms(BigDFT_Atoms *atoms)
{
  guint i;

  atoms->comment = strdup("Test from memory generation.");
  fprintf(stdout, "Test BigDFT_Atoms structure set n types.\n");
  bigdft_atoms_set_n_types(atoms, 2);
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stdout, " Atom %d -> type %d\n", i, atoms->iatype[i]);
  fprintf(stdout, "Test BigDFT_Atoms structure syncing.\n");
  atoms->geocode = 'F';
  atoms->atomnames[0] = strdup("O");
  atoms->atomnames[1] = strdup("H");
  atoms->alat[0] = 5.;
  atoms->alat[1] = 5.;
  atoms->alat[2] = 5.;
  sprintf(atoms->format, "ascii");
  sprintf(atoms->units, "bohr");
  bigdft_atoms_sync(atoms);
  fprintf(stdout, "Test BigDFT_Atoms structure set n atoms.\n");
  bigdft_atoms_set_n_atoms(atoms, 3);
  fprintf(stdout, "Test BigDFT_Atoms structure direct access.\n");
  atoms->iatype[0] = 1;
  atoms->iatype[1] = 2;
  atoms->iatype[2] = 2;
  atoms->rxyz.data[3 * 0 + 0] = 0.;
  atoms->rxyz.data[3 * 0 + 1] = 0.;
  atoms->rxyz.data[3 * 0 + 2] = 0.;
  atoms->rxyz.data[3 * 1 + 0] = 0.5;
  atoms->rxyz.data[3 * 1 + 1] = 0.;
  atoms->rxyz.data[3 * 1 + 2] = 0.5;
  atoms->rxyz.data[3 * 2 + 0] = -0.5;
  atoms->rxyz.data[3 * 2 + 1] = 0.;
  atoms->rxyz.data[3 * 2 + 2] = 0.5;
  /* bigdft_atoms_set_displacement(atoms, 0.001); */
  fprintf(stdout, "Test BigDFT_Atoms pseudo-potential on the fly.\n");
  bigdft_atoms_set_psp(atoms, 1, 1, (const gchar*)0);
}

int main(int argc, char **argv)
{
  BigDFT_Atoms *atoms;
  guint i, j, nelec, n[3 * 6];
  double *radii, peak;
  double h[3] = {0.45, 0.45, 0.45};
#define CRMULT 5.
#define FRMULT 8.
  BigDFT_Inputs *in;
  BigDFT_Wf *wf;
  BigDFT_LocReg *lr;
  BigDFT_Proj *proj;
  BigDFT_LocalFields *denspot;
  BigDFT_Data *data;
  BigDFT_Energs *energs;
#ifdef HAVE_GLIB
  GMainLoop *loop;
#endif

  int out_pipe[2], stdout_fileno_old;

#ifdef HAVE_GLIB
  /* g_mem_set_vtable (glib_mem_profiler_table); */
  g_type_init();
  g_thread_init(NULL);
  loop = g_main_loop_new(NULL, FALSE);
#endif

  FC_FUNC_(memocc_verbose, MEMOCC_VERBOSE)();

  fprintf(stdout, "Test BigDFT_Wf structure creation.\n");
  wf = bigdft_wf_new(100);
  atoms = BIGDFT_ATOMS(wf->lzd);

  atoms->comment = strdup("Test BigDFT_Atoms reading a wrong file.");
  if (!bigdft_atoms_set_structure_from_file(BIGDFT_ATOMS(wf->lzd), "truc"))
    {
      bigdft_wf_free(wf);

      wf = bigdft_wf_new(100);
      atoms = BIGDFT_ATOMS(wf->lzd);
    }
  fprintf(stdout, " Ok\n");

  init_atoms(atoms);
  bigdft_atoms_write(atoms, "output");
  radii = bigdft_atoms_get_radii(atoms, 0., 0., 0.);
  bigdft_locreg_set_radii(BIGDFT_LOCREG(wf->lzd), radii);
  g_free(radii);
  bigdft_locreg_set_size(BIGDFT_LOCREG(wf->lzd), h, CRMULT, FRMULT);
  bigdft_lzd_init_d(wf->lzd);
  bigdft_locreg_init_wfd(BIGDFT_LOCREG(wf->lzd));
  output_locreg(BIGDFT_LOCREG(wf->lzd));

  if (argc > 1)
    chdir(argv[1]);

  /* Typical cluster run. */
  fprintf(stdout, "Test BigDFT_Inputs structure creation.\n");
  in = bigdft_inputs_new("test");
  fprintf(stdout, " base Ok\n");
  bigdft_atoms_set_symmetries(BIGDFT_ATOMS(wf->lzd), !in->disableSym, -1., in->elecfield);
  bigdft_inputs_parse_additional(in, BIGDFT_ATOMS(wf->lzd));
  fprintf(stdout, " additional Ok\n");
  output_inputs(in);

  fprintf(stdout, "Test BigDFT_Wf define.\n");
  nelec = bigdft_wf_define(wf, in, 0, 1);
  output_llr(wf->lzd);
  fprintf(stdout, "Test BigDFT_Wf free.\n");
  bigdft_wf_free(wf);
  fprintf(stdout, " Ok\n");
  fprintf(stdout, "Test BigDFT_Inputs free.\n");
  bigdft_inputs_free(in);
  fprintf(stdout, " Ok\n");

  /* Test some set methods. */
  wf = bigdft_wf_new(0);
  init_atoms(BIGDFT_ATOMS(wf->lzd));
  radii = bigdft_atoms_get_radii(BIGDFT_ATOMS(wf->lzd), 0., 0., 0.);
  bigdft_locreg_set_radii(BIGDFT_LOCREG(wf->lzd), radii);
  g_free(radii);
  bigdft_locreg_set_size(BIGDFT_LOCREG(wf->lzd), h, CRMULT, FRMULT);
  fprintf(stdout, "Test BigDFT_Lzd set nlr.\n");
  bigdft_lzd_set_n_locreg(wf->lzd, 3);
  for (i = 0; i < wf->lzd->nlr; i++)
    {
      fprintf(stdout, "Test BigDFT_Locreg set dims %d.\n", i);
      for (j = 0; j < 3 * 6; j++)
        n[j] = i * (1 + j % 3);
      bigdft_locreg_set_d_dims(wf->lzd->Llr[i], n, n + 3, n + 6, n + 9, n + 12, n + 15);
      bigdft_locreg_set_wfd_dims(wf->lzd->Llr[i], 3 * i, i, 10 * 3 * i, 10 * i);
    }
  fprintf(stdout, " -> check gives %d\n", bigdft_lzd_check(wf->lzd));
  output_llr(wf->lzd);
  /* Keep a reference on the second locreg. */
  lr = wf->lzd->Llr[1];
  g_object_ref(lr);
  /* Modify lzd to suppress the second and third locreg. */
  fprintf(stdout, "Test BigDFT_Lzd reallocation.\n");
  bigdft_lzd_set_n_locreg(wf->lzd, 1);
  bigdft_locreg_set_d_dims(wf->lzd->Llr[0], n, n + 3, n + 6, n + 9, n + 12, n + 15);
  bigdft_locreg_set_wfd_dims(wf->lzd->Llr[0], 3 * i, i, 10 * 3 * i, 10 * i);
  fprintf(stdout, " -> check gives %d\n", bigdft_lzd_check(wf->lzd));
  fprintf(stdout, "Test BigDFT_Locreg separation.\n");
  fprintf(stdout, " -> check gives %d\n", bigdft_locreg_check(lr));
  bigdft_wf_free(wf);
  g_object_unref(lr);

  if (argc > 2)
    {
      stdout_fileno_old = redirect_init(out_pipe);
      FC_FUNC_(memocc_report, MEMOCC_REPORT)();
      redirect_dump(out_pipe, stdout_fileno_old);
    }

  /* Test restarting the wavefunction definition. */
  in = bigdft_inputs_new("test");
  wf = bigdft_wf_new(in->inputPsiId);
  fprintf(stdout, "Test BigDFT_Atoms structure creation from file.\n");
  if (!bigdft_atoms_set_structure_from_file(BIGDFT_ATOMS(wf->lzd), "posinp.ascii"))
    {
      fprintf(stdout, "Problem with your file.\n");
      return 1;
    }
  output_atoms(BIGDFT_ATOMS(wf->lzd));

  bigdft_atoms_set_symmetries(BIGDFT_ATOMS(wf->lzd), !in->disableSym, -1., in->elecfield);
  bigdft_inputs_parse_additional(in, BIGDFT_ATOMS(wf->lzd));

  fprintf(stdout, "Test BigDFT_Atoms pseudo-potential evaluation.\n");
  bigdft_atoms_set_psp(BIGDFT_ATOMS(wf->lzd), in->ixc, in->nspin, (const gchar*)0);
  radii = bigdft_atoms_get_radii(BIGDFT_ATOMS(wf->lzd), in->crmult, in->frmult, 0.);
  bigdft_locreg_set_radii(BIGDFT_LOCREG(wf->lzd), radii);
  g_free(radii);
  bigdft_locreg_set_size(BIGDFT_LOCREG(wf->lzd), in->h, in->crmult, in->frmult);
  bigdft_lzd_init_d(wf->lzd);
  bigdft_locreg_init_wfd(BIGDFT_LOCREG(wf->lzd));
  output_locreg(BIGDFT_LOCREG(wf->lzd));

  nelec = bigdft_wf_define(wf, in, 0, 1);
  fprintf(stdout, " System has %d electrons.\n", nelec);

  fprintf(stdout, "Test BigDFT_Proj structure creation.\n");
  proj = bigdft_proj_new(BIGDFT_LOCREG(wf->lzd), BIGDFT_ORBS(wf), in->frmult);
  fprintf(stdout, " System has %d projectors, and %d elements.\n",
          proj->nproj, proj->nprojel);

  if (argc > 2)
    {
      fprintf(stdout, "Test memory estimation.\n");
      stdout_fileno_old = redirect_init(out_pipe);
      peak = bigdft_memory_get_peak(4, BIGDFT_LOCREG(wf->lzd), in, BIGDFT_ORBS(wf), proj);
      redirect_dump(out_pipe, stdout_fileno_old);
      fprintf(stdout, " Memory peak will reach %f octets.\n", peak);
    }

  fprintf(stdout, "Test BigDFT_LocalFields creation.\n");
  denspot = bigdft_localfields_new(wf->lzd, in, 0, 1);
  fprintf(stdout, " Meta data are %f %f %f  -  %d  -  %f\n",
          denspot->h[0], denspot->h[1], denspot->h[2],
          denspot->rhov_is, denspot->psoffset);
  bigdft_localfields_create_poisson_kernels(denspot, wf->lzd, in, 0, 1);

  /* Block here in a main loop. */
#ifdef HAVE_GLIB
  data = run_bigdft(in, proj, denspot, wf, loop);
  g_timeout_add(100000, exit_loop, (gpointer)loop);
  g_main_loop_run(loop);
#else
  data = run_bigdft(in, proj, denspot, wf, (gpointer)0);
#endif
  energs = data->energs;
  g_free(data);

  fprintf(stdout, " Total energy after relaxation is %gHt.\n", energs->eKS);

  fprintf(stdout, "Test BigDFT_LocalFields free.\n");
  bigdft_localfields_free(denspot);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Proj free.\n");
  bigdft_proj_free(proj);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Wf free.\n");
  bigdft_wf_free(wf);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Inputs free.\n");
  bigdft_inputs_free(in);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Energs free.\n");
  bigdft_energs_free(energs);
  fprintf(stdout, " Ok\n");

  if (argc > 2)
    {
      stdout_fileno_old = redirect_init(out_pipe);
      FC_FUNC_(memocc_report, MEMOCC_REPORT)();
      redirect_dump(out_pipe, stdout_fileno_old);

#ifdef HAVE_GLIB
      /* g_mem_profile(); */
#endif
    }

  return 0;
}

static int redirect_init(int out_pipe[2])
{
  int stdout_fileno_old;

  /* Flush before redirecting. */
  fflush(stdout);

  /* Make a pipe to redirect stdout. */
  stdout_fileno_old = dup(STDOUT_FILENO);
  pipe(out_pipe);
  dup2(out_pipe[1], STDOUT_FILENO);

  return stdout_fileno_old;
}

static void redirect_dump(int out_pipe[2], int stdout_fileno_old)
{
  gchar foutput[MAX_FORTRAN_OUTPUT];
  ssize_t ncount;
  guint i;
  long flags;

  /* Flush before reconnecting. */
  fflush(stdout);
  /* Reconnect stdout. */
  dup2(stdout_fileno_old, STDOUT_FILENO);

  /* Make the reading pipe non blocking. */
  flags = fcntl(out_pipe[0], F_GETFL);
  flags |= O_NONBLOCK;
  fcntl(out_pipe[0], F_SETFL, flags);

  /* Write Fortran output with prefix... */
  foutput[0] = '\0';
  ncount = read(out_pipe[0], foutput, MAX_FORTRAN_OUTPUT);
  foutput[ncount] = '\0';
  fprintf(stdout, "Fortran: ");
  for (i = 0; foutput[i]; i++)
    {
      putc(foutput[i], stdout);
      if (foutput[i] == '\n' && foutput[i + 1] != '\0')
        fprintf(stdout, "Fortran: ");
    }

  /* Close the pipes. */
  close(out_pipe[0]);
  close(out_pipe[1]);
}

#ifdef HAVE_GLIB
static gboolean exit_loop(gpointer data)
{
  g_main_loop_quit((GMainLoop*)data);
  fprintf(stdout, "Error, signals timeout.\n");
  return FALSE;
}
static void onIterDone(BigDFT_OptLoop *optloop, BigDFT_Energs *energs, gpointer data)
{
  fprintf(stdout, "Callback for 'iter-done-wavefunctions' signal at iter %d -> %gHt (%g).\n",
          optloop->iter, energs->eKS, optloop->gnrm);
}
#endif

static gpointer optimize_psi_thread(gpointer data)
{
  BigDFT_Data *container = (BigDFT_Data*)data;
  BigDFT_OptLoop *p;
  
  fprintf(stdout, " Calculation of optimization started.\n");
  p = bigdft_optloop_new();
  p->itermax = 2;
  g_signal_connect(G_OBJECT(p), "iter-wavefunctions",
                   G_CALLBACK(onIterDone), (gpointer)0);
  g_signal_connect(G_OBJECT(p), "done-wavefunctions",
                   G_CALLBACK(onIterDone), (gpointer)0);
  bigdft_optloop_sync_to_fortran(p);
  bigdft_wf_optimization_loop(container->wf, container->denspot, container->proj,
                              container->energs, p, 0, 1);
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(container->wf));
  g_object_unref(G_OBJECT(container->denspot));
  /* g_object_unref(G_OBJECT(container->proj)); */
  g_idle_add((GSourceFunc)g_main_loop_quit, container->loop);
#endif
  fprintf(stdout, " Calculation of optimization finished.\n");

  bigdft_optloop_free(p);
  return (gpointer)0;
}
static gboolean optimize_psi(gpointer data)
{
  BigDFT_Data *ct = (BigDFT_Data*)data;
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
  GError *error = (GError*)0;
#endif

#ifdef HAVE_GLIB
  g_object_ref(G_OBJECT(ct->denspot));
  g_object_ref(G_OBJECT(ct->wf));
  /* g_object_ref(G_OBJECT(ct->proj)); */
#endif
#ifdef G_THREADS_ENABLED
  ld_thread = g_thread_create(optimize_psi_thread, ct, FALSE, &error);
#else
  optimize_psi_thread(ct);
#endif

  return FALSE;
}

static gpointer calculate_psi_0_thread(gpointer data)
{
  BigDFT_Data *container = (BigDFT_Data*)data;
  
  fprintf(stdout, " Calculation of input guess started.\n");
  bigdft_wf_calculate_psi0(container->wf, container->denspot,
                           container->proj, container->energs, 0, 1);
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(container->wf));
  g_object_unref(G_OBJECT(container->denspot));
  g_object_unref(G_OBJECT(container->energs));
  /* g_object_unref(G_OBJECT(container->proj)); */
  /* Chain up with the SCF loop. */
  g_idle_add(optimize_psi, data);
#endif
  fprintf(stdout, " Calculation of input guess finished.\n");

  return (gpointer)0;
}
static gboolean calculate_psi_0(gpointer data)
{
  BigDFT_Data *ct = (BigDFT_Data*)data;
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
  GError *error = (GError*)0;
#endif

#ifdef HAVE_GLIB
  g_object_ref(G_OBJECT(ct->denspot));
  g_object_ref(G_OBJECT(ct->wf));
  g_object_ref(G_OBJECT(ct->energs));
  /* g_object_ref(G_OBJECT(ct->proj)); */
#endif
#ifdef G_THREADS_ENABLED
  ld_thread = g_thread_create(calculate_psi_0_thread, ct, FALSE, &error);
#else
  calculate_psi_0_thread(ct);
#endif

  return FALSE;
}

#ifdef HAVE_GLIB
static void onVExtReady(BigDFT_LocalFields *denspot, gpointer data)
{
  /* Copy the data of V_Ext to main process memory for later use. */

  /* Chain up with the input guess. */
  g_idle_add(calculate_psi_0, data);
}
static void onDensReady(BigDFT_LocalFields *denspot, guint istep, gpointer data)
{
  guint i;
  double dens;

  /* Do something with the density. */
  fprintf(stdout, "Callback for \"density-ready\" signal at iter %d.\n", istep);
  dens = 0.;
  for (i = 0; i < denspot->ni[0] * denspot->ni[1] * denspot->ni[2]; i++)
    dens += denspot->rhov[i];
  dens *= denspot->h[0] * denspot->h[1] * denspot->h[2];
  fprintf(stdout, " Density calculated by C is %16.16f.\n", dens);
}
static void onPsiReady(BigDFT_Wf *wf, guint iter, gpointer data)
{
  const double *psic;
  double *psir, *psii;
  guint size, i, n;
  double minDens, maxDens;

  fprintf(stdout, "Callback for 'psi-ready' signal at iter %d.\n", iter);

  psic = bigdft_wf_get_psi_compress(wf, 1, 4, BIGDFT_SPIN_UP, BIGDFT_REAL, &size, 0);
  fprintf(stdout, " Band 4 has %ld bytes.\n", size * sizeof(double));
  
  minDens = G_MAXDOUBLE;
  maxDens = 0.;
  n = BIGDFT_LOCREG(wf->lzd)->ni[0] * 
    BIGDFT_LOCREG(wf->lzd)->ni[1] * 
    BIGDFT_LOCREG(wf->lzd)->ni[2];
  psir = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf->lzd), psic);
  if (BIGDFT_ORBS(wf)->nspinor == 2)
    psii = bigdft_locreg_convert_to_isf(BIGDFT_LOCREG(wf->lzd), psic + size);
  for (i = 0; i < n; i++)
    {
      psir[i] *= psir[i];
      if (BIGDFT_ORBS(wf)->nspinor == 2)
        psir[i] += psii[i] * psii[i];
      minDens = MIN(minDens, psir[i]);
      maxDens = MAX(maxDens, psir[i]);
    }
  fprintf(stdout, " Band 4 has min partial density %g and max %g.\n", minDens, maxDens);

  if (iter == 0)
    bigdft_wf_write_psi_compress(wf, "wave", BIGDFT_WF_FORMAT_PLAIN, psic, 1, 4, BIGDFT_SPIN_UP, size);

  g_free(psir);
  if (BIGDFT_ORBS(wf)->nspinor == 2)
    g_free(psii);
}
static void onEksReady(BigDFT_Energs *energs, guint iter, gpointer data)
{
  fprintf(stdout, "Callback for 'eks-ready' signal at iter %d -> %gHt.\n", iter, energs->eKS);
}
#endif

static gpointer calculate_ionic_pot_thread(gpointer data)
{
  BigDFT_Data *container = (BigDFT_Data*)data;
  
  fprintf(stdout, " Calculation of ionic potential started.\n");
  bigdft_localfields_create_effective_ionic_pot(container->denspot, container->wf->lzd,
                                                container->in, 0, 1);
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(container->denspot));
#endif
  fprintf(stdout, " Calculation of ionic potential finished.\n");
  return (gpointer)0;
}
static gboolean calculate_ionic_pot(gpointer data)
{
  BigDFT_Data *ct = (BigDFT_Data*)data;
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
  GError *error = (GError*)0;
#endif

#ifdef HAVE_GLIB
  g_object_ref(G_OBJECT(ct->denspot));
#endif
#ifdef G_THREADS_ENABLED
  ld_thread = g_thread_create(calculate_ionic_pot_thread, ct, FALSE, &error);
#else
  calculate_ionic_pot_thread(ct);
#endif

  return FALSE;
}

static BigDFT_Data* run_bigdft(BigDFT_Inputs *in, BigDFT_Proj *proj,
                               BigDFT_LocalFields *denspot, BigDFT_Wf *wf, gpointer data)
{
  BigDFT_Data *ct;

  ct = g_malloc(sizeof(BigDFT_Data));
  ct->denspot = denspot;
  ct->in      = in;
  ct->proj    = proj;
  ct->wf      = wf;
  ct->energs  = bigdft_energs_new();
#ifdef HAVE_GLIB
  ct->loop    = (GMainLoop*)data;
  g_signal_connect(G_OBJECT(ct->denspot), "v-ext-ready",
                   G_CALLBACK(onVExtReady), (gpointer)ct);
  g_signal_connect(G_OBJECT(ct->denspot), "density-ready",
                   G_CALLBACK(onDensReady), (gpointer)ct);
  g_signal_connect(G_OBJECT(ct->wf), "psi-ready",
                   G_CALLBACK(onPsiReady), (gpointer)ct);
  g_signal_connect(G_OBJECT(ct->energs), "eks-ready",
                   G_CALLBACK(onEksReady), (gpointer)ct);
  g_idle_add(calculate_ionic_pot, (gpointer)ct);
#else
  calculate_ionic_pot((gpointer)ct);
  calculate_psi_0((gpointer)ct);
#endif

  return ct;
}
