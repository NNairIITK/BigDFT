#include <config.h>

#include <glib-object.h>

#include <bigdft.h>

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#define MAX_FORTRAN_OUTPUT 4096

#ifdef HAVE_GLIB
static gboolean exit_loop(gpointer data);
static void onVExtReady(BigDFT_LocalFields *denspot, gpointer data);
#endif
static int redirect_init(int out_pipe[2]);
static void redirect_dump(int out_pipe[2], int stdout_fileno_old);
static void calculate_ionic_pot(BigDFT_LocalFields *denspot, BigDFT_Inputs *in);

int main(guint argc, char **argv)
{
  BigDFT_Atoms *atoms;
  guint i, n, nelec;
  double *radii, peak;
  BigDFT_Lzd *lzd;
  double h[3] = {0.45, 0.45, 0.45};
  gboolean *cgrid, *fgrid;
#define CRMULT 5.
#define FRMULT 8.
  BigDFT_Inputs *in;
  BigDFT_Wf *wf;
  BigDFT_Proj *proj;
  BigDFT_LocalFields *denspot;
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

  fprintf(stdout, "Test BigDFT_Atoms structure creation.\n");
  atoms = bigdft_atoms_new();
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
  bigdft_atoms_write(atoms, "output");
  fprintf(stdout, "Test BigDFT_Atoms pseudo-potential on the fly.\n");
  bigdft_atoms_set_psp(atoms, 1);
  radii = bigdft_atoms_get_radii(atoms);
  for (i = 0; i < atoms->ntypes; i++)
    fprintf(stdout, " Type %d, radii %f %f %f\n", i,
            radii[i], radii[atoms->ntypes + i], radii[atoms->ntypes * 2 + i]);
  g_free(radii);
  fprintf(stdout, "Test BigDFT_Atoms free.\n");
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(atoms));
#else
  bigdft_atoms_free(atoms);
#endif
  fprintf(stdout, " Ok\n");

  if (argc > 1)
    chdir(argv[1]);

  fprintf(stdout, "Test BigDFT_Atoms structure creation from file.\n");
  atoms = bigdft_atoms_new_from_file("posinp.ascii");
  if (!atoms)
    {
      fprintf(stdout, "Problem with your file.\n");
      return 1;
    }
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stdout, " Atom %d, coord. %f %f %f, type %d\n", i,
            atoms->rxyz.data[3 * i], atoms->rxyz.data[3 * i + 1],
            atoms->rxyz.data[3 * i + 2], atoms->iatype[i]);
  fprintf(stdout, " Box [%c], is in %f %f %f\n", atoms->geocode,
          atoms->alat[0], atoms->alat[1], atoms->alat[2]);

  fprintf(stdout, "Test BigDFT_Atoms pseudo-potential evaluation.\n");
  bigdft_atoms_set_psp(atoms, 11);
  radii = bigdft_atoms_get_radii(atoms);
  for (i = 0; i < atoms->ntypes; i++)
    fprintf(stdout, " Type %d, radii %f %f %f\n", i,
            radii[i], radii[atoms->ntypes + i], radii[atoms->ntypes * 2 + i]);
  
  fprintf(stdout, "Test BigDFT_Lzd structure creation.\n");
  lzd = bigdft_lzd_new(atoms, radii, h, CRMULT, FRMULT);
  bigdft_locreg_set_wave_descriptors(BIGDFT_LOCREG(lzd));
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stdout, " Atoms %d, coord. %10.6f %10.6f %10.6f '%2s', type %d\n",
            i, atoms->rxyz.data[3 * i], atoms->rxyz.data[3 * i + 1],
            atoms->rxyz.data[3 * i + 2], atoms->atomnames[atoms->iatype[i] - 1],
            atoms->iatype[i]);
  fprintf(stdout, " Box is in   %f %f %f\n", atoms->alat[0], atoms->alat[1], atoms->alat[2]);
  fprintf(stdout, " Shift is     %f %f  %f\n", atoms->shift[0], atoms->shift[1], atoms->shift[2]);
  fprintf(stdout, " Geocode is  %c\n", BIGDFT_LOCREG(lzd)->geocode);
  fprintf(stdout, " Grid is     %9d %9d %9d\n", BIGDFT_LOCREG(lzd)->n[0], BIGDFT_LOCREG(lzd)->n[1], BIGDFT_LOCREG(lzd)->n[2]);
  fprintf(stdout, " Int grid is %9d %9d %9d\n", BIGDFT_LOCREG(lzd)->ni[0], BIGDFT_LOCREG(lzd)->ni[1], BIGDFT_LOCREG(lzd)->ni[2]);
  fprintf(stdout, " H grids are %9.9g %9.9g %9.9g\n", BIGDFT_LOCREG(lzd)->h[0], BIGDFT_LOCREG(lzd)->h[1], BIGDFT_LOCREG(lzd)->h[2]);

  fprintf(stdout, "Test calculation of grid points.\n");
  cgrid = bigdft_atoms_get_grid(atoms, radii, CRMULT, BIGDFT_LOCREG(lzd)->n);
  for (i = 0, n = 0; i < (BIGDFT_LOCREG(lzd)->n[0] + 1) * (BIGDFT_LOCREG(lzd)->n[1] + 1) * (BIGDFT_LOCREG(lzd)->n[2] + 1); i++)
    if (cgrid[i])
      n += 1;
  fprintf(stdout, " Coarse grid has %7d points.\n", n);
  fgrid = bigdft_atoms_get_grid(atoms, radii + atoms->ntypes, FRMULT, BIGDFT_LOCREG(lzd)->n);
  for (i = 0, n = 0; i < (BIGDFT_LOCREG(lzd)->n[0] + 1) * (BIGDFT_LOCREG(lzd)->n[1] + 1) * (BIGDFT_LOCREG(lzd)->n[2] + 1); i++)
    if (fgrid[i])
      n += 1;
  fprintf(stdout, " Fine grid has   %7d points.\n", n);

  g_free(cgrid);
  g_free(fgrid);

  fprintf(stdout, "Test BigDFT_Inputs structure creation.\n");
  in = bigdft_inputs_new("test");
  fprintf(stdout, " Read 'test.dft' file: %d\n", (in->files & BIGDFT_INPUTS_DFT));
  fprintf(stdout, " Input variables are %f %f %f  -  %f %f  -  %d\n",
          in->h[0], in->h[1], in->h[2], in->crmult, in->frmult, in->ixc);

  bigdft_atoms_set_symmetries(atoms, !in->disableSym, -1., in->elecfield);
  bigdft_inputs_parse_additional(in, atoms);

  fprintf(stdout, "Test BigDFT_Wf structure creation.\n");
  wf = bigdft_wf_new(lzd, in, 0, 1, &nelec);
  fprintf(stdout, " System has %d electrons.\n", nelec);

  fprintf(stdout, "Test BigDFT_Proj structure creation.\n");
  proj = bigdft_proj_new(BIGDFT_LOCREG(lzd), BIGDFT_ORBS(wf), in->frmult);
  fprintf(stdout, " System has %d projectors, and %d elements.\n",
          proj->nproj, proj->nprojel);

  if (argc > 2)
    {
      fprintf(stdout, "Test memory estimation.\n");
      stdout_fileno_old = redirect_init(out_pipe);
      peak = bigdft_memory_get_peak(4, BIGDFT_LOCREG(lzd), in, BIGDFT_ORBS(wf), proj);
      redirect_dump(out_pipe, stdout_fileno_old);
      fprintf(stdout, " Memory peak will reach %f octets.\n", peak);
    }

  fprintf(stdout, "Test BigDFT_LocalFields creation.\n");
  denspot = bigdft_localfields_new(BIGDFT_LOCREG(lzd), in, 0, 1);
  fprintf(stdout, " Meta data are %f %f %f  -  %d  -  %f\n",
          denspot->h[0], denspot->h[1], denspot->h[2],
          denspot->rhov_is, denspot->psoffset);
  fprintf(stdout, " Add linear zone description.\n");
  bigdft_lzd_setup_linear(lzd, BIGDFT_ORBS(wf), in, atoms,0, 1);

  /* Use a thread to generate the ionic potential... */
  fprintf(stdout, " Calculate ionic potential.\n");
  calculate_ionic_pot(denspot, in);

  /* Block here in a main loop. */
#ifdef HAVE_GLIB
  g_signal_connect(G_OBJECT(denspot), "v-ext-ready",
                   G_CALLBACK(onVExtReady), (gpointer)loop);
  g_timeout_add_seconds(5, exit_loop, (gpointer)loop);
  g_main_loop_run(loop);
#endif

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

  fprintf(stdout, "Test BigDFT_Lzd free.\n");
  bigdft_lzd_free(lzd);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Atoms free.\n");
  bigdft_atoms_free(atoms);
  fprintf(stdout, " Ok\n");

  g_free(radii);

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

static void onVExtReady(BigDFT_LocalFields *denspot, gpointer data)
{
  /* Copy the data of V_Ext to main process memory for later use. */
  g_idle_add((GSourceFunc)g_main_loop_quit, data);
}
#endif

struct ionicpot_
{
  BigDFT_LocalFields *denspot;
  BigDFT_Inputs *in;
};
static gpointer calculate_ionic_pot_thread(gpointer data)
{
  struct ionicpot_ *container = (struct ionicpot_*)data;
  
  fprintf(stdout, " Calculation of ionic potential started.\n");
  bigdft_localfields_create_effective_ionic_pot(container->denspot, container->in, 0, 1);
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(container->denspot));
#endif
  fprintf(stdout, " Calculation of ionic potential finished.\n");
  g_free(container);

  return (gpointer)0;
}

static void calculate_ionic_pot(BigDFT_LocalFields *denspot, BigDFT_Inputs *in)
{
#ifdef G_THREADS_ENABLED
  GThread *ld_thread;
  GError *error = (GError*)0;
#endif
  struct ionicpot_ *ct;

  ct = g_malloc(sizeof(struct ionicpot_));
  ct->denspot = denspot;
  ct->in = in;
#ifdef HAVE_GLIB
  g_object_ref(G_OBJECT(denspot));
#endif
#ifdef G_THREADS_ENABLED
  ld_thread = g_thread_create(calculate_ionic_pot_thread, ct, FALSE, &error);
#else
  calculate_ionic_pot_thread(ct);
#endif
}
