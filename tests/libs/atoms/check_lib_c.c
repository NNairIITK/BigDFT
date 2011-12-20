#include <config.h>
#include <bigdft.h>

#include <stdio.h>

int main(int argc, char **argv)
{
  BigDFT_Atoms *atoms;
  guint i, n;
  double *radii;
  BigDFT_Glr *glr;
  double h[3] = {0.45, 0.45, 0.45};
  int *cgrid, *fgrid;
#define CRMULT 5.
#define FRMULT 8.

  atoms = bigdft_atoms_new();
  bigdft_atoms_set_n_atoms(atoms, 3);
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stderr, "%d -> %d\n", i, atoms->iatype[i]);
  bigdft_atoms_free(atoms);

  atoms = bigdft_atoms_new_from_file(argv[1]);
  if (!atoms)
    {
      fprintf(stderr, "Problem with your file.\n");
      return 1;
    }
  bigdft_atoms_set_psp(atoms, 11);
  radii = bigdft_atoms_get_radii(atoms);

  for (i = 0; i  < atoms->nat; i++)
    fprintf(stderr, "%f %f %f %d\n", atoms->rxyz.data[3 * i], atoms->rxyz.data[3 * i + 1],
            atoms->rxyz.data[3 * i + 2], atoms->iatype[i]);
  for (i = 0; i < atoms->ntypes; i++)
    fprintf(stderr, "%d %f %f %f\n", i, radii[i], radii[atoms->ntypes + i], radii[atoms->ntypes * 2 + i]);
  
  glr = bigdft_glr_new(atoms, radii, h, CRMULT, FRMULT);
  for (i = 0; i  < atoms->nat; i++)
    fprintf(stderr, "%f %f %f %d\n", atoms->rxyz.data[3 * i], atoms->rxyz.data[3 * i + 1],
            atoms->rxyz.data[3 * i + 2], atoms->iatype[i]);
  fprintf(stderr, "Grid is in %d %d %d\n", glr->n[0], glr->n[1], glr->n[2]);

  cgrid = bigdft_fill_logrid(atoms, glr->n, radii, CRMULT, h);
  for (i = 0, n = 0; i < (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1); i++)
    if (cgrid[i] != 0)
      n += 1;
  fprintf(stderr, "Coarse grid has %d points.\n", n);
  fgrid = bigdft_fill_logrid(atoms, glr->n, radii + atoms->ntypes, FRMULT, h);
  for (i = 0, n = 0; i < (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1); i++)
    if (fgrid[i] != 0)
      n += 1;
  fprintf(stderr, "Coarse grid has %d points.\n", n);

  g_free(cgrid);
  g_free(fgrid);

  bigdft_glr_free(glr);

  bigdft_atoms_free(atoms);

  FC_FUNC_(memocc_report, MEMOCC_REPORT)();

  return 0;
}
