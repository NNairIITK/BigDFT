#include <config.h>
#include <bigdft.h>

#include <stdio.h>
#include <unistd.h>

#define MAX_FORTRAN_OUTPUT 4096

int main(int argc, char **argv)
{
  BigDFT_Atoms *atoms;
  guint i, n, nelec;
  double *radii, peak;
  BigDFT_Glr *glr;
  double h[3] = {0.45, 0.45, 0.45};
  int *cgrid, *fgrid;
#define CRMULT 5.
#define FRMULT 8.
  BigDFT_Inputs *in;
  BigDFT_Orbs *orbs;

  int out_pipe[2], stdout_fileno_old;
  gchar foutput[MAX_FORTRAN_OUTPUT];
  ssize_t ncount;

  if (argc > 1)
    chdir(argv[1]);

  fprintf(stdout, "Test BigDFT_Atoms structure creation.\n");
  atoms = bigdft_atoms_new();
  bigdft_atoms_set_n_atoms(atoms, 3);
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stdout, " Atom %d -> type %d\n", i, atoms->iatype[i]);
  bigdft_atoms_free(atoms);

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
  
  fprintf(stdout, "Test BigDFT_Glr structure creation.\n");
  glr = bigdft_glr_new_with_wave_descriptors(atoms, radii, h, CRMULT, FRMULT);
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stdout, " Atoms %d, coord. %10.6f %10.6f %10.6f '%2s', type %d\n",
            i, atoms->rxyz.data[3 * i], atoms->rxyz.data[3 * i + 1],
            atoms->rxyz.data[3 * i + 2], atoms->atomnames[atoms->iatype[i] - 1],
            atoms->iatype[i]);
  fprintf(stdout, " Box is in %f %f %f\n", atoms->alat[0], atoms->alat[1], atoms->alat[2]);
  fprintf(stdout, " Shift is  %f %f %f\n", atoms->shift[0], atoms->shift[1], atoms->shift[2]);
  fprintf(stdout, " Grid is   %9d %9d %9d\n", glr->n[0], glr->n[1], glr->n[2]);

  fprintf(stdout, "Test calculation of grid points.\n");
  cgrid = bigdft_fill_logrid(atoms, glr->n, radii, CRMULT, h);
  for (i = 0, n = 0; i < (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1); i++)
    if (cgrid[i] != 0)
      n += 1;
  fprintf(stdout, " Coarse grid has %7d points.\n", n);
  fgrid = bigdft_fill_logrid(atoms, glr->n, radii + atoms->ntypes, FRMULT, h);
  for (i = 0, n = 0; i < (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1); i++)
    if (fgrid[i] != 0)
      n += 1;
  fprintf(stdout, " Fine grid has   %7d points.\n", n);

  g_free(cgrid);
  g_free(fgrid);

  fprintf(stdout, "Test BigDFT_Inputs structure creation.\n");
  in = bigdft_inputs_new("test");
  fprintf(stdout, " Read 'test.dft' file: %d\n", (in->files & BIGDFT_INPUTS_DFT));
  fprintf(stdout, " Input variables are %f %f %f  -  %f %f  -  %d\n",
          in->h[0], in->h[1], in->h[2], in->crmult, in->frmult, in->ixc);

  bigdft_atoms_set_symmetries(atoms, !in->disableSym, in->elecfield);
  bigdft_inputs_parse_additional(in, atoms);

  fprintf(stdout, "Test BigDFT_Orbs structure creation.\n");
  orbs = bigdft_orbs_new(atoms, in, 0, 1, &nelec);
  fprintf(stderr, " System has %d electrons.\n", nelec);

  fprintf(stdout, "Test memory estimation.\n");
  peak = bigdft_memory_peak(4, glr, in, orbs);
  fprintf(stderr, " Memory peak will reach %f octets.\n", peak);

  fprintf(stdout, "Test BigDFT_Orbs free.\n");
  bigdft_orbs_free(orbs);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Inputs free.\n");
  bigdft_inputs_free(in);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Glr free.\n");
  bigdft_glr_free(glr);
  fprintf(stdout, " Ok\n");

  fprintf(stdout, "Test BigDFT_Atoms free.\n");
  bigdft_atoms_free(atoms);
  fprintf(stdout, " Ok\n");

  fflush(stdout);

  /* Make a pipe to redirect stdout. */
  stdout_fileno_old = dup(STDOUT_FILENO);
  pipe(out_pipe);
  dup2(out_pipe[1], STDOUT_FILENO);

  FC_FUNC_(memocc_report, MEMOCC_REPORT)();
  fflush(stdout);

  /* Reconnect stdout. */
  dup2(stdout_fileno_old, STDOUT_FILENO);

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

  return 0;
}
