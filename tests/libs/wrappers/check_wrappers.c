#include <bigdft.h>

int main(int argc, const char **argv)
{
  guint iproc, nproc, igroup, ngroup;
  int ierr;

  const gchar *types[] = {"C", "O", NULL};
  double alat[3] = {10., 0., 10.};
  BigDFT_Atoms *atoms;

  BigDFT_Inputs *ins;

  BigDFT_Run *run;

  BigDFT_Goutput *outs;

  ierr = bigdft_init(&iproc, &nproc, &igroup, &ngroup, 0);

  atoms = bigdft_atoms_new();
  bigdft_atoms_set_types(atoms, types);
  bigdft_atoms_set_n_atoms(atoms, 2);
  atoms->iatype[0] = 1;
  atoms->iatype[1] = 2;
  atoms->rxyz.data[3] = 1.23;
  bigdft_atoms_set_geometry(atoms, 'S', alat, "bohr");
  if (iproc == 0) bigdft_atoms_write(atoms, "posinp", "xyz");

  ins = bigdft_inputs_new("truc", iproc);

  run = bigdft_run_new_from_objects(atoms, ins, NULL, iproc);
  bigdft_inputs_unref(ins);
  bigdft_atoms_unref(atoms);

  outs = bigdft_run_calculate(run, iproc, nproc);

  bigdft_run_unref(run);

  ierr = bigdft_finalize();

  return 0;
}
