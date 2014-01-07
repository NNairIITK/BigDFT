#include <bigdft.h>

int main(int argc, const char **argv)
{
  guint iproc, nproc, igroup, ngroup;
  int ierr, i;

#define NAT 2
  const gchar *types[] = {"C", "O", NULL};
  guint iatype[NAT] = {0, 1};
  const gchar *xyz[] = {"0.", "0.", "0.", NULL, "1.23", "0.", "0.", NULL};
  const gchar *frz[NAT] = {NULL, "fy"};
  const gchar *sp[NAT] = {"-1", "1"};
  const gchar *alat[] = {"10.", ".inf", "10.", NULL};
  const gchar *hgrids[] = {"2/5", "0.55", "0.55", NULL};
  /* const gchar *ngkpt[] = {"2", "3", "4", NULL}; */
  /* const gchar *shiftk1[] = {"0.", "0.", "0.", NULL}; */
  /* const gchar *shiftk2[] = {"0.1", "-0.2", "0.2", NULL}; */

  BigDFT_Dict *dict;
  BigDFT_DictIter root, coords;

  BigDFT_Atoms *atoms;

  BigDFT_Run *run;

  BigDFT_Goutput *outs;

  ierr = bigdft_init(&iproc, &nproc, &igroup, &ngroup, 0);

  dict = bigdft_dict_new(&root);

  bigdft_dict_insert(dict, "posinp", NULL);
  bigdft_dict_set_array(dict, "Cell", alat);
  bigdft_dict_insert(dict, "Positions", &coords);
  for (i = 0; i < NAT; i++)
    {
      bigdft_dict_append(dict, NULL);
      bigdft_dict_set_array(dict, types[iatype[i]], xyz + (4 * i));
      if (frz[i]) bigdft_dict_set(dict, "Frozen", frz[i]);
      if (sp[i]) bigdft_dict_set(dict, "IGSpin", sp[i]);
      bigdft_dict_move_to(dict, &coords);
    }
  bigdft_dict_move_to(dict, &root);

  bigdft_dict_insert(dict, "dft", NULL);
  bigdft_dict_set(dict, "ixc", "PBE (ABINIT)");
  bigdft_dict_set_array(dict, "hgrids", hgrids);
  bigdft_dict_set(dict, "nspin", "2");
  bigdft_dict_set(dict, "itermax", "4");

  bigdft_dict_dump(dict);

  run = bigdft_run_new_from_dict(dict, TRUE);
  bigdft_dict_unref(dict);

  atoms = bigdft_run_get_atoms(run);
  if (iproc == 0) bigdft_atoms_write(atoms, "posinp", "yaml");
  bigdft_atoms_unref(atoms);

  if (iproc == 0) bigdft_run_dump(run, "input.yaml", TRUE);

  /* bigdft_run_memoryEstimation(run, iproc, nproc); */
  /* print_memory_estimation_(run->mem.data); */

  outs = bigdft_run_calculate(run, iproc, nproc);

  bigdft_run_unref(run);
  bigdft_goutput_unref(outs);

  ierr = bigdft_finalize();

  return 0;
}
