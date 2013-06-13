#include <bigdft.h>

int main(int argc, const char **argv)
{
  guint iproc, nproc, igroup, ngroup;
  int ierr;

  const gchar *types[] = {"C", "O", NULL};
  double alat[3] = {10., 0., 10.};
  BigDFT_Atoms *atoms;

  BigDFT_Inputs *ins;
  const gchar *hgrids[] = {"0.55", "0.55", "0.55", NULL};

  BigDFT_Run *run;

  BigDFT_Goutput *outs;

  ierr = bigdft_init(&iproc, &nproc, &igroup, &ngroup, 0);

  atoms = bigdft_atoms_new();
  bigdft_atoms_set_types(atoms, types);
  bigdft_atoms_set_n_atoms(atoms, 2);
  atoms->iatype[0] = 1;
  atoms->iatype[1] = 2;
  atoms->rxyz[3] = 1.23;
  atoms->ifrztyp[1] = 2;
  atoms->natpol[0] = 99;
  atoms->natpol[1] = 101;
  bigdft_atoms_set_geometry(atoms, 'S', alat, "bohr");
  if (iproc == 0) bigdft_atoms_write(atoms, "posinp", "yaml");

  ins = bigdft_inputs_new("truc", nproc);
  bigdft_inputs_set(ins, INPUTS_IXC, "PBE");
  bigdft_inputs_set_array(ins, INPUTS_HGRIDS, hgrids);
  bigdft_inputs_set(ins, INPUTS_NSPIN, "2");
  /* bigdft_inputs_set(ins, INPUTS_INPUTPSIID, "-500"); */

  run = bigdft_run_new_from_objects(atoms, ins, NULL, iproc);
  bigdft_inputs_unref(ins);
  bigdft_atoms_unref(atoms);

  outs = bigdft_run_calculate(run, iproc, nproc);

  bigdft_run_unref(run);

  ierr = bigdft_finalize();

  return 0;
}
