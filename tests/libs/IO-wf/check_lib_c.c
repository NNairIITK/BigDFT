#include <config.h>
#include <bigdft.h>

#include <stdio.h>

#define WAVEFILE "data/wavefunction-k001-UR.bin.b0001"
#define IORBP 23

int main(int argc, char **argv)
{
  char *filename;
  f90_pointer_double_4D *psiscf;
  int n[3], nspinor;
  int i, j, k, ind, ntot, norbu, norbd, nkpt;
  int iorb, ispin, ikpt, ispinor;
  double h[3], nrm;

  if (argc > 1)
    filename = argv[1];
  else
    filename = WAVEFILE;

  fprintf(stdout, " --- Test read_wave_to_isf_etsf() from %s ---\n", filename);

  if (bigdft_read_wave_descr(filename, &norbu, &norbd, &nkpt, &nspinor,
                             &iorb, &ispin, &ikpt, &ispinor))
    fprintf(stdout, " ETSF wavefunction file (nou, nod, nk, sp):   %6d %6d %6d %6d\n",
            norbu, norbd, nkpt, nspinor);
  else
    return 1;

  psiscf = bigdft_read_wave_to_isf(filename, IORBP, h, n, &nspinor);
  if (!psiscf)
    return 1;
  /* for (i = 0; i < F90_POINTER_SIZE; i++) */
  /*   fprintf(stdout, "%d %p\n", i, psiscf->info[i]); */

  fprintf(stdout, " hgrid values for iscf representation:     %9.6f %9.6f %9.6f\n",
          h[0], h[1], h[2]);
  fprintf(stdout, " number of points in iscf representation:  %9d %9d %9d\n",
          n[0], n[1], n[2]);

  nrm  = 0.;
  ind  = 0;
  ntot = n[0] * n[1] * n[2];
  for (k = 0; k < n[2]; k++)
    for (j = 0; j < n[1]; j++)
      for (i = 0; i < n[0]; i++)
        {
          nrm += psiscf->data[ind] * psiscf->data[ind];
          if (nspinor == 2)
            nrm += psiscf->data[ind + ntot] * psiscf->data[ind + ntot];
          ind += 1;
        }
  fprintf(stdout, " norm of orbital %d:                      %12.8f\n", IORBP, nrm);

  bigdft_free_wave_to_isf(psiscf);

  fflush(stdout);

  FC_FUNC_(memocc_report, MEMOCC_REPORT)();

  return 0;
}
