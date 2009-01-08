#ifndef _structDeflocpot_
#define _structDeflocpot_



//maximum size of the shared memory array
//conceived for maximize occupancy on a hardware of compute
//capability 1.2 and higher (1024 threads at same time on a given multiprocessor)
#define MAX_SHARED_SIZE 3072 //16*256 4 kB (should be =~ 3.9 kB, try also 3072)
#define HALF_WARP_SIZE 16 // for all architectures
#define NUM_LINES 16 
#define HW_ELEM 1 //this is HALF_WARP_SIZE/NUM_LINES

//parameter related to the Magic Filter convolution
//lowfil + lupfil + 1  must be a multiple of 16
#define LOWFIL 8
#define LUPFIL 7

//convolution filters
#define MFIL0   8.4334247333529341094733325815816e-7f
#define MFIL1  -0.1290557201342060969516786758559028e-4f
#define MFIL2   0.8762984476210559564689161894116397e-4f
#define MFIL3  -0.30158038132690463167163703826169879e-3f
#define MFIL4   0.174723713672993903449447812749852942e-2f
#define MFIL5  -0.942047030201080385922711540948195075e-2f
#define MFIL6   0.2373821463724942397566389712597274535e-1f
#define MFIL7   0.612625895831207982195380597e-1f
#define MFIL8   0.9940415697834003993178616713f
#define MFIL9  -0.604895289196983516002834636e-1f
#define MFIL10 -0.2103025160930381434955489412839065067e-1f
#define MFIL11  0.1337263414854794752733423467013220997e-1f
#define MFIL12 -0.344128144493493857280881509686821861e-2f
#define MFIL13  0.49443227688689919192282259476750972e-3f
#define MFIL14 -0.5185986881173432922848639136911487e-4f
#define MFIL15  2.72734492911979659657715313017228e-6f


typedef struct  _parMF
{
  unsigned int ElementsPerBlock;

  int thline[HALF_WARP_SIZE]; //line considered by a thread within the half-warp
  int thelem[HALF_WARP_SIZE]; //elements considered by a thread within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps

} parMF_t;

#endif
