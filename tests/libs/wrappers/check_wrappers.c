#include <bigdft.h>

int main(int argc, const char **argv)
{
  guint iproc, nproc, igroup, ngroup;
  int ierr;

  const gchar *types[] = {"C", "O", NULL};
  double alat[3] = {10., 0., 10.};
  BigDFT_Atoms *atoms;

  BigDFT_Inputs *ins;
  const gchar *hgrids[] = {"2/5", "0.55", "0.55", NULL};
  /* const gchar *ngkpt[] = {"2", "3", "4", NULL}; */
  /* const gchar *shiftk1[] = {"0.", "0.", "0.", NULL}; */
  /* const gchar *shiftk2[] = {"0.1", "-0.2", "0.2", NULL}; */

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

  ins = bigdft_inputs_new("");
  bigdft_inputs_set(ins, INPUTS_IXC, "PBE");
  bigdft_inputs_set_array(ins, INPUTS_HGRIDS, hgrids);
  bigdft_inputs_set(ins, INPUTS_NSPIN, "2");
  bigdft_inputs_set(ins, INPUTS_ITERMAX, "4");
  /* bigdft_inputs_set(ins, INPUTS_DEBUG, "Yes"); */
  /* bigdft_inputs_set(ins, INPUTS_INPUTPSIID, "-500"); */
  /* bigdft_inputs_set(ins, INPUTS_GEOPT_METHOD, "DIIS"); */
  /* bigdft_inputs_set(ins, INPUTS_BETAX, "1."); */
  /* bigdft_inputs_set(ins, INPUTS_KPT_METHOD, "MPGrid"); */
  /* bigdft_inputs_set_array(ins, INPUTS_NGKPT, ngkpt); */
  /* bigdft_inputs_set_array_at(ins, INPUTS_SHIFTK, 0, shiftk1); */
  /* bigdft_inputs_set_array_at(ins, INPUTS_SHIFTK, 1, shiftk2); */

  run = bigdft_run_new_from_objects(atoms, ins, NULL, iproc, TRUE);
  bigdft_inputs_dump(ins, "input.yaml", TRUE);
  bigdft_inputs_unref(ins);
  bigdft_atoms_unref(atoms);

  outs = bigdft_run_calculate(run, iproc, nproc);

  bigdft_run_unref(run);
  bigdft_goutput_unref(outs);

  ierr = bigdft_finalize();

  return 0;
}
