
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****f* BigDFT/fourier_dim
!! FUNCTION
!!   Give a number n_next > n compatible for the FFT
!! SOURCE
!!
subroutine fourier_dim(n,n_next)
  implicit none
  !Arguments
  integer, intent(in) :: n
  integer, intent(out) :: n_next

  !Local variables
  integer, parameter :: ndata = 292
  !Radix 2, 3, 4, 5, 6, 7, 8
  integer, dimension(ndata), parameter :: idata = (/   &
       3,      4,      5,      6,      7,      8,      9,     12,  &
      14,     15,     16,     18,     20,     21,     24,     25,  &
      27,     28,     30,     32,     35,     36,     40,     42,  &
      45,     48,     50,     54,     56,     60,     63,     64,  &
      70,     72,     75,     80,     81,     84,     90,     96,  &
     100,    105,    108,    112,    120,    125,    126,    128,  &
     135,    140,    144,    150,    160,    162,    168,    175,  &
     180,    189,    192,    200,    210,    216,    224,    225,  &
     240,    243,    252,    256,    270,    280,    288,    300,  &
     315,    320,    324,    336,    360,    375,    378,    384,  &
     400,    405,    420,    432,    448,    450,    480,    486,  &
     500,    504,    512,    525,    540,    560,    567,    576,  &
     600,    625,    630,    640,    648,    672,    675,    700,  &
     720,    729,    750,    756,    768,    800,    810,    840,  &
     864,    875,    896,    900,    945,    960,    972,   1000,  &
    1008,   1024,   1050,   1080,   1120,   1125,   1134,   1152,  &
    1200,   1215,   1260,   1280,   1296,   1344,   1350,   1400,  &
    1440,   1458,   1500,   1512,   1536,   1575,   1600,   1620,  &
    1680,   1701,   1728,   1792,   1800,   1875,   1890,   1920,  &
    1944,   2000,   2016,   2025,   2048,   2100,   2160,   2240,  &
    2250,   2268,   2304,   2400,   2430,   2500,   2520,   2560,  &
    2592,   2625,   2688,   2700,   2800,   2835,   2880,   3000,  &
    3024,   3072,   3125,   3150,   3200,   3240,   3360,   3375,  &
    3402,   3456,   3500,   3584,   3600,   3750,   3780,   3840,  &
    3888,   4000,   4032,   4050,   4096,   4200,   4320,   4375,  &
    4480,   4500,   4536,   4608,   4725,   4800,   5000,   5040,  &
    5120,   5184,   5250,   5376,   5400,   5600,   5625,   5670,  &
    5760,   6000,   6048,   6144,   6300,   6400,   6480,   6720,  &
    6750,   6912,   7000,   7168,   7200,   7500,   7560,   7680,  &
    7875,   8000,   8064,   8192,   8400,   8640,   8960,   9000,  &
    9072,   9216,   9375,   9450,   9600,  10000,  10080,  10240,  &
   10368,  10500,  10752,  10800,  11200,  11250,  11520,  12000,  &
   12096,  12288,  12500,  12600,  12800,  13125,  13440,  13824,  &
   14000,  14336,  14400,  15000,  15120,  15360,  15625,  15750,  &
   16000,  16128,  16384,  16800,  17280,  17500,  17920,  18000,  &
   18432,  18750,  19200,  20000,  20160,  20480,  21000,  21504,  &
   21875,  22400,  23040,  24000 /)
  integer :: i

  loop_data: do i=1,ndata
     if (n <= idata(i)) then
        n_next = idata(i)
        return
     end if
  end do loop_data
  write(unit=*,fmt=*) "fourier_dim: ",n," is bigger than ",idata(ndata)
  stop
end subroutine fourier_dim
!!***


!  Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .


! --------------------------------------------------------------
!   3-dimensional complex-complex FFT routine: 
!   When compared to the best vendor implementations on RISC architectures 
!   it gives close to optimal performance (perhaps loosing 20 percent in speed)
!   and it is significanly faster than many not so good vendor implementations 
!   as well as other portable FFT's. 
!   On all vector machines tested so far (Cray, NEC, Fujitsu) is 
!   was significantly faster than the vendor routines
! The theoretical background is described in :
! 1) S. Goedecker: Rotating a three-dimensional array in optimal
! positions for vector processing: Case study for a three-dimensional Fast
! Fourier Transform, Comp. Phys. Commun. \underline{76}, 294 (1993)
! Citing of this reference is greatly appreciated if the routines are used 
! for scientific work.


! Presumably good compiler flags:
! IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
! with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ; 
!                   a.out
! DEC: f90 -O3 -arch ev67 -pipeline
! with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ; 
!                   prun -N1 -c4 a.out


!-----------------------------------------------------------

! FFT PART -----------------------------------------------------------------


subroutine ctrig(n,trig,after,before,now,isign,ic)
!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

!     Different factorizations affect the performance
!     Factoring 64 as 4*4*4 might for example be faster on some machines than 8*8.

  implicit real(kind=8) (a-h,o-z)
! Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
! Arguments
  integer :: n,isign,ic
  integer :: after(7),before(7),now(7)
  real(kind=8) :: trig(2,nfft_max)
! Local variables
  integer, parameter :: ndata = 292
  integer, dimension(7,ndata) :: idata 
! The factor 6 is only allowed in the first place!
        data ((idata(i,j),i=1,7),j=1,ndata) /                        &
       3,   3,   1,   1,   1,   1,   1,  &
       4,   4,   1,   1,   1,   1,   1,  &
       5,   5,   1,   1,   1,   1,   1,  &
       6,   6,   1,   1,   1,   1,   1,  &
       7,   7,   1,   1,   1,   1,   1,  &
       8,   8,   1,   1,   1,   1,   1,  &
       9,   3,   3,   1,   1,   1,   1,  &
      12,   4,   3,   1,   1,   1,   1,  &
      14,   7,   2,   1,   1,   1,   1,  &
      15,   5,   3,   1,   1,   1,   1,  &
      16,   4,   4,   1,   1,   1,   1,  &
      18,   6,   3,   1,   1,   1,   1,  &
      20,   5,   4,   1,   1,   1,   1,  &
      21,   7,   3,   1,   1,   1,   1,  &
      24,   8,   3,   1,   1,   1,   1,  &
      25,   5,   5,   1,   1,   1,   1,  &
      27,   3,   3,   3,   1,   1,   1,  &
      28,   7,   4,   1,   1,   1,   1,  &
      30,   6,   5,   1,   1,   1,   1,  &
      32,   8,   4,   1,   1,   1,   1,  &
      35,   7,   5,   1,   1,   1,   1,  &
      36,   4,   3,   3,   1,   1,   1,  &
      40,   8,   5,   1,   1,   1,   1,  &
      42,   7,   6,   1,   1,   1,   1,  &
      45,   5,   3,   3,   1,   1,   1,  &
      48,   4,   4,   3,   1,   1,   1,  &
      50,   5,   5,   2,   1,   1,   1,  &
      54,   6,   3,   3,   1,   1,   1,  &
      56,   8,   7,   1,   1,   1,   1,  &
      60,   5,   4,   3,   1,   1,   1,  &
      63,   7,   3,   3,   1,   1,   1,  &
      64,   8,   8,   1,   1,   1,   1,  &
      70,   7,   5,   2,   1,   1,   1,  &
      72,   8,   3,   3,   1,   1,   1,  &
      75,   5,   5,   3,   1,   1,   1,  &
      80,   5,   4,   4,   1,   1,   1,  &
      81,   3,   3,   3,   3,   1,   1,  &
      84,   7,   4,   3,   1,   1,   1,  &
      90,   6,   5,   3,   1,   1,   1,  &
      96,   8,   4,   3,   1,   1,   1,  &
     100,   5,   5,   4,   1,   1,   1,  &
     105,   7,   5,   3,   1,   1,   1,  &
     108,   4,   3,   3,   3,   1,   1,  &
     112,   7,   4,   4,   1,   1,   1,  &
     120,   8,   5,   3,   1,   1,   1,  &
     125,   5,   5,   5,   1,   1,   1,  &
     126,   7,   6,   3,   1,   1,   1,  &
     128,   8,   4,   4,   1,   1,   1,  &
     135,   5,   3,   3,   3,   1,   1,  &
     140,   7,   5,   4,   1,   1,   1,  &
     144,   6,   8,   3,   1,   1,   1,  &
     150,   6,   5,   5,   1,   1,   1,  &
     160,   8,   5,   4,   1,   1,   1,  &
     162,   6,   3,   3,   3,   1,   1,  &
     168,   8,   7,   3,   1,   1,   1,  &
     175,   7,   5,   5,   1,   1,   1,  &
     180,   5,   4,   3,   3,   1,   1,  &
     189,   7,   3,   3,   3,   1,   1,  &
     192,   6,   8,   4,   1,   1,   1,  &
     200,   8,   5,   5,   1,   1,   1,  &
     210,   7,   6,   5,   1,   1,   1,  &
     216,   8,   3,   3,   3,   1,   1,  &
     224,   8,   7,   4,   1,   1,   1,  &
     225,   5,   5,   3,   3,   1,   1,  &
     240,   6,   8,   5,   1,   1,   1,  &
     243,   3,   3,   3,   3,   3,   1,  &
     252,   7,   4,   3,   3,   1,   1,  &
     256,   8,   8,   4,   1,   1,   1,  &
     270,   6,   5,   3,   3,   1,   1,  &
     280,   8,   7,   5,   1,   1,   1,  &
     288,   8,   4,   3,   3,   1,   1,  &
     300,   5,   5,   4,   3,   1,   1,  &
     315,   7,   5,   3,   3,   1,   1,  &
     320,   5,   4,   4,   4,   1,   1,  &
     324,   4,   3,   3,   3,   3,   1,  &
     336,   7,   4,   4,   3,   1,   1,  &
     360,   8,   5,   3,   3,   1,   1,  &
     375,   5,   5,   5,   3,   1,   1,  &
     378,   7,   6,   3,   3,   1,   1,  &
     384,   8,   4,   4,   3,   1,   1,  &
     400,   5,   5,   4,   4,   1,   1,  &
     405,   5,   3,   3,   3,   3,   1,  &
     420,   7,   5,   4,   3,   1,   1,  &
     432,   4,   4,   3,   3,   3,   1,  &
     448,   8,   8,   7,   1,   1,   1,  &
     450,   6,   5,   5,   3,   1,   1,  &
     480,   8,   5,   4,   3,   1,   1,  &
     486,   6,   3,   3,   3,   3,   1,  &
     500,   5,   5,   5,   4,   1,   1,  &
     504,   8,   7,   3,   3,   1,   1,  &
     512,   8,   8,   8,   1,   1,   1,  &
     525,   7,   5,   5,   3,   1,   1,  &
     540,   5,   4,   3,   3,   3,   1,  &
     560,   7,   5,   4,   4,   1,   1,  &
     567,   7,   3,   3,   3,   3,   1,  &
     576,   4,   4,   4,   3,   3,   1,  &
     600,   8,   5,   5,   3,   1,   1,  &
     625,   5,   5,   5,   5,   1,   1,  &
     630,   7,   6,   5,   3,   1,   1,  &
     640,   8,   5,   4,   4,   1,   1,  &
     648,   8,   3,   3,   3,   3,   1,  &
     672,   8,   7,   4,   3,   1,   1,  &
     675,   5,   5,   3,   3,   3,   1,  &
     700,   7,   5,   5,   4,   1,   1,  &
     720,   5,   4,   4,   3,   3,   1,  &
     729,   3,   3,   3,   3,   3,   3,  &
     750,   6,   5,   5,   5,   1,   1,  &
     756,   7,   4,   3,   3,   3,   1,  &
     768,   4,   4,   4,   4,   3,   1,  &
     800,   8,   5,   5,   4,   1,   1,  &
     810,   6,   5,   3,   3,   3,   1,  &
     840,   8,   7,   5,   3,   1,   1,  &
     864,   8,   4,   3,   3,   3,   1,  &
     875,   7,   5,   5,   5,   1,   1,  &
     896,   8,   7,   4,   4,   1,   1,  &
     900,   5,   5,   4,   3,   3,   1,  &
     945,   7,   5,   3,   3,   3,   1,  &
     960,   5,   4,   4,   4,   3,   1,  &
     972,   4,   3,   3,   3,   3,   3,  &
    1000,   8,   5,   5,   5,   1,   1,  &
    1008,   8,   7,   6,   3,   1,   1,  &
    1024,   4,   4,   4,   4,   4,   1,  &
    1050,   7,   6,   5,   5,   1,   1,  &
    1080,   6,   5,   4,   3,   3,   1,  &
    1120,   8,   7,   5,   4,   1,   1,  &
    1125,   5,   5,   5,   3,   3,   1,  &
    1134,   7,   6,   3,   3,   3,   1,  &
    1152,   6,   4,   4,   4,   3,   1,  &
    1200,   6,   8,   5,   5,   1,   1,  &
    1215,   5,   3,   3,   3,   3,   3,  &
    1260,   7,   5,   4,   3,   3,   1,  &
    1280,   8,   8,   5,   4,   1,   1,  &
    1296,   6,   8,   3,   3,   3,   1,  &
    1344,   8,   7,   6,   4,   1,   1,  &
    1350,   6,   5,   5,   3,   3,   1,  &
    1400,   8,   7,   5,   5,   1,   1,  &
    1440,   6,   5,   4,   4,   3,   1,  &
    1458,   6,   3,   3,   3,   3,   3,  &
    1500,   5,   5,   5,   4,   3,   1,  &
    1512,   8,   7,   3,   3,   3,   1,  &
    1536,   6,   8,   8,   4,   1,   1,  &
    1575,   7,   5,   5,   3,   3,   1,  &
    1600,   8,   8,   5,   5,   1,   1,  &
    1620,   5,   4,   3,   3,   3,   3,  &
    1680,   8,   7,   6,   5,   1,   1,  &
    1701,   7,   3,   3,   3,   3,   3,  &
    1728,   6,   8,   4,   3,   3,   1,  &
    1792,   8,   8,   7,   4,   1,   1,  &
    1800,   6,   5,   5,   4,   3,   1,  &
    1875,   5,   5,   5,   5,   3,   1,  &
    1890,   7,   6,   5,   3,   3,   1,  &
    1920,   6,   5,   4,   4,   4,   1,  &
    1944,   6,   4,   3,   3,   3,   3,  &
    2000,   5,   5,   5,   4,   4,   1,  &
    2016,   8,   7,   4,   3,   3,   1,  &
    2025,   5,   5,   3,   3,   3,   3,  &
    2048,   8,   4,   4,   4,   4,   1,  &
    2100,   7,   5,   5,   4,   3,   1,  &
    2160,   6,   8,   5,   3,   3,   1,  &
    2240,   7,   5,   4,   4,   4,   1,  &
    2250,   6,   5,   5,   5,   3,   1,  &
    2268,   7,   4,   3,   3,   3,   3,  &
    2304,   6,   8,   4,   4,   3,   1,  &
    2400,   6,   5,   5,   4,   4,   1,  &
    2430,   6,   5,   3,   3,   3,   3,  &
    2500,   5,   5,   5,   5,   4,   1,  &
    2520,   8,   7,   5,   3,   3,   1,  &
    2560,   8,   5,   4,   4,   4,   1,  &
    2592,   6,   4,   4,   3,   3,   3,  &
    2625,   7,   5,   5,   5,   3,   1,  &
    2688,   8,   7,   4,   4,   3,   1,  &
    2700,   5,   5,   4,   3,   3,   3,  &
    2800,   7,   5,   5,   4,   4,   1,  &
    2835,   7,   5,   3,   3,   3,   3,  &
    2880,   6,   8,   5,   4,   3,   1,  &
    3000,   6,   5,   5,   5,   4,   1,  &
    3024,   7,   4,   4,   3,   3,   3,  &
    3072,   6,   8,   4,   4,   4,   1,  &
    3125,   5,   5,   5,   5,   5,   1,  &
    3150,   7,   6,   5,   5,   3,   1,  &
    3200,   8,   5,   5,   4,   4,   1,  &
    3240,   6,   5,   4,   3,   3,   3,  &
    3360,   8,   7,   5,   4,   3,   1,  &
    3375,   5,   5,   5,   3,   3,   3,  &
    3402,   7,   6,   3,   3,   3,   3,  &
    3456,   6,   4,   4,   4,   3,   3,  &
    3500,   7,   5,   5,   5,   4,   1,  &
    3584,   8,   8,   8,   7,   1,   1,  &
    3600,   6,   8,   5,   5,   3,   1,  &
    3750,   6,   5,   5,   5,   5,   1,  &
    3780,   7,   5,   4,   3,   3,   3,  &
    3840,   6,   8,   5,   4,   4,   1,  &
    3888,   6,   8,   3,   3,   3,   3,  &
    4000,   8,   5,   5,   5,   4,   1,  &
    4032,   7,   4,   4,   4,   3,   3,  &
    4050,   6,   5,   5,   3,   3,   3,  &
    4096,   8,   8,   4,   4,   4,   1,  &
    4200,   8,   7,   5,   5,   3,   1,  &
    4320,   6,   5,   4,   4,   3,   3,  &
    4375,   7,   5,   5,   5,   5,   1,  &
    4480,   8,   7,   5,   4,   4,   1,  &
    4500,   5,   5,   5,   4,   3,   3,  &
    4536,   8,   7,   3,   3,   3,   3,  &
    4608,   6,   8,   8,   4,   3,   1,  &
    4725,   7,   5,   5,   3,   3,   3,  &
    4800,   6,   8,   5,   5,   4,   1,  &
    5000,   8,   5,   5,   5,   5,   1,  &
    5040,   7,   5,   4,   4,   3,   3,  &
    5120,   8,   8,   5,   4,   4,   1,  &
    5184,   6,   8,   4,   3,   3,   3,  &
    5250,   7,   6,   5,   5,   5,   1,  &
    5376,   7,   4,   4,   4,   4,   3,  &
    5400,   6,   5,   5,   4,   3,   3,  &
    5600,   8,   7,   5,   5,   4,   1,  &
    5625,   5,   5,   5,   5,   3,   3,  &
    5670,   7,   6,   5,   3,   3,   3,  &
    5760,   6,   8,   8,   5,   3,   1,  &
    6000,   6,   8,   5,   5,   5,   1,  &
    6048,   8,   7,   4,   3,   3,   3,  &
    6144,   6,   8,   8,   4,   4,   1,  &
    6300,   7,   5,   5,   4,   3,   3,  &
    6400,   8,   8,   5,   5,   4,   1,  &
    6480,   6,   8,   5,   3,   3,   3,  &
    6720,   7,   5,   4,   4,   4,   3,  &
    6750,   6,   5,   5,   5,   3,   3,  &
    6912,   6,   8,   4,   4,   3,   3,  &
    7000,   8,   7,   5,   5,   5,   1,  &
    7168,   7,   4,   4,   4,   4,   4,  &
    7200,   6,   5,   5,   4,   4,   3,  &
    7500,   5,   5,   5,   5,   4,   3,  &
    7560,   7,   6,   5,   4,   3,   3,  &
    7680,   6,   8,   8,   5,   4,   1,  &
    7875,   7,   5,   5,   5,   3,   3,  &
    8000,   8,   8,   5,   5,   5,   1,  &
    8064,   7,   6,   4,   4,   4,   3,  &
    8192,   8,   8,   8,   4,   4,   1,  &
    8400,   8,   7,   6,   5,   5,   1,  &
    8640,   8,   8,   5,   3,   3,   3,  &
    8960,   8,   8,   7,   5,   4,   1,  &
    9000,   8,   5,   5,   5,   3,   3,  &
    9072,   8,   7,   6,   3,   3,   3,  &
    9216,   6,   8,   8,   8,   3,   1,  &
    9375,   5,   5,   5,   5,   5,   3,  &
    9450,   7,   6,   5,   5,   3,   3,  &
    9600,   8,   5,   5,   4,   4,   3,  &
   10000,   5,   5,   5,   5,   4,   4,  &
   10080,   7,   6,   5,   4,   4,   3,  &
   10240,   8,   8,   8,   5,   4,   1,  &
   10368,   6,   8,   8,   3,   3,   3,  &
   10500,   7,   5,   5,   5,   4,   3,  &
   10752,   8,   8,   7,   6,   4,   1,  &
   10800,   6,   8,   5,   5,   3,   3,  &
   11200,   8,   8,   7,   5,   5,   1,  &
   11250,   6,   5,   5,   5,   5,   3,  &
   11520,   8,   8,   5,   4,   3,   3,  &
   12000,   8,   5,   5,   5,   4,   3,  &
   12096,   8,   7,   6,   4,   3,   3,  &
   12288,   8,   8,   8,   8,   3,   1,  &
   12500,   5,   5,   5,   5,   5,   4,  &
   12600,   7,   6,   5,   5,   4,   3,  &
   12800,   8,   8,   8,   5,   5,   1,  &
   13125,   7,   5,   5,   5,   5,   3,  &
   13440,   7,   6,   5,   4,   4,   4,  &
   13824,   8,   8,   8,   3,   3,   3,  &
   14000,   7,   5,   5,   5,   4,   4,  &
   14336,   8,   7,   4,   4,   4,   4,  &
   14400,   8,   8,   5,   5,   3,   3,  &
   15000,   8,   5,   5,   5,   5,   3,  &
   15120,   8,   7,   6,   5,   3,   3,  &
   15360,   6,   8,   8,   8,   5,   1,  &
   15625,   5,   5,   5,   5,   5,   5,  &
   15750,   7,   6,   5,   5,   5,   3,  &
   16000,   8,   5,   5,   5,   4,   4,  &
   16128,   8,   7,   6,   4,   4,   3,  &
   16384,   8,   8,   8,   8,   4,   1,  &
   16800,   7,   6,   5,   5,   4,   4,  &
   17280,   6,   8,   8,   5,   3,   3,  &
   17500,   7,   5,   5,   5,   5,   4,  &
   17920,   8,   7,   5,   4,   4,   4,  &
   18000,   6,   8,   5,   5,   5,   3,  &
   18432,   8,   8,   8,   4,   3,   3,  &
   18750,   6,   5,   5,   5,   5,   5,  &
   19200,   8,   8,   5,   5,   4,   3,  &
   20000,   8,   5,   5,   5,   5,   4,  &
   20160,   8,   7,   6,   5,   4,   3,  &
   20480,   8,   8,   8,   8,   5,   1,  &
   21000,   7,   6,   5,   5,   5,   4,  &
   21504,   8,   7,   6,   4,   4,   4,  &
   21875,   7,   5,   5,   5,   5,   5,  &
   22400,   8,   7,   5,   5,   4,   4,  &
   23040,   8,   8,   8,   5,   3,   3,  &
   24000,   8,   8,   5,   5,   5,   3   /

  do i=1,ndata
     if (n.eq.idata(1,i)) then
        ic=0
        do j=1,6
           itt=idata(1+j,i)
           if (itt.gt.1) then
              ic=ic+1
              now(j)=idata(1+j,i)
           else
              goto 1000
           endif
        end do
        goto 1000
     endif
  end do
  print *,'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
  write(6,"(15(i5))") (idata(1,i),i=1,ndata)
  stop
1000 continue

  after(1)=1
  before(ic)=1
  do i=2,ic
     after(i)=after(i-1)*now(i-1)
     before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
  end do

!  write(6,"(6(i3))") (after(i),i=1,ic)
!  write(6,"(6(i3))") (now(i),i=1,ic)
!  write(6,"(6(i3))") (before(i),i=1,ic)

  twopi=6.283185307179586d0
  angle=real(isign,kind=8)*twopi/real(n,kind=8)
  if (mod(n,2).eq.0) then
     nh=n/2
     trig(1,1)=1.d0
     trig(2,1)=0.d0
     trig(1,nh+1)=-1.d0
     trig(2,nh+1)=0.d0
     do i=1,nh-1
        trigc=cos(real(i,kind=8)*angle)
        trigs=sin(real(i,kind=8)*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
     end do
  else
     nh=(n-1)/2
     trig(1,1)=1.d0
     trig(2,1)=0.d0
     do i=1,nh
        trigc=cos(real(i,kind=8)*angle)
        trigs=sin(real(i,kind=8)*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
      end do
  endif

END SUBROUTINE ctrig


!ccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fftstp(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

  implicit real(kind=8) (a-h,o-z)
  integer, parameter :: nfft_max=24000
! Arguments
  integer :: mm,nfft,m,nn,n,after,before,isign
  real(kind=8) :: trig(2,nfft_max),zin(2,mm,m),zout(2,nn,n)
! Local variables
  integer :: atn,atb
  atn=after*now
  atb=after*before

!         sqrt(.5d0)
  rt2i=0.7071067811865475d0

! First now == 2
  if (now.eq.2) then
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        do j=1,nfft
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           r2=zin(1,j,nin2)
           s2=zin(2,j,nin2)
           zout(1,j,nout1)= r2 + r1
           zout(2,j,nout1)= s2 + s1
           zout(1,j,nout2)= r1 - r2
           zout(2,j,nout2)= s1 - s2
        end do
     end do

   ! Big loop
     Big_Loop: do ia=2,after
        ias=ia-1
        if (2*ias.eq.after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    zout(1,j,nout1)= r1 - r2
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r2 + r1
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s1 - s2
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s2 + s1
                 end do
              end do
           endif
        else if (4*ias.eq.after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r - s)*rt2i
                    s2=(r + s)*rt2i
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           endif
        else if (4*ias.eq.3*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(r - s)*rt2i
                    zout(1,j,nout1)= r1 - r2
                    zout(2,j,nout1)= s2 + s1
                    zout(1,j,nout2)= r2 + r1
                    zout(2,j,nout2)= s1 - s2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(s - r)*rt2i
                    s2=(r + s)*rt2i
                    zout(1,j,nout1)= r2 + r1
                    zout(2,j,nout1)= s1 - s2
                    zout(1,j,nout2)= r1 - r2
                    zout(2,j,nout2)= s2 + s1
                 end do
              end do
           endif
        else
           itrig=ias*before+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nout1=nout1+atn
              nout2=nout1+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 zout(1,j,nout1)= r2 + r1
                 zout(2,j,nout1)= s2 + s1
                 zout(1,j,nout2)= r1 - r2
                 zout(2,j,nout2)= s1 - s2
              end do
           end do
        endif

     end do Big_Loop
!    End of Big_loop
! End of if (now.eq.2)

  else if (now.eq.4) then
     if (isign.eq.1) then 
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r=r1 + r3
              s=r2 + r4
              zout(1,j,nout1) = r + s
              zout(1,j,nout3) = r - s
              r=r1 - r3
              s=s2 - s4
              zout(1,j,nout2) = r - s 
              zout(1,j,nout4) = r + s
              r=s1 + s3
              s=s2 + s4
              zout(2,j,nout1) = r + s 
              zout(2,j,nout3) = r - s
              r=s1 - s3
              s=r2 - r4
              zout(2,j,nout2) = r + s 
              zout(2,j,nout4) = r - s
           end do
        end do
        do ia=2,after
           ias=ia-1
           if (2*ias.eq.after) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r-s)*rt2i
                    s2=(r+s)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(r + s)*rt2i
                    s4=(r - s)*rt2i
                    r=r1 - r3
                    s=r2 - r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 + r3
                    s=s2 - s4
                    zout(1,j,nout2) = r - s 
                    zout(1,j,nout4) = r + s
                    r=s1 + s3
                    s=s2 + s4
                    zout(2,j,nout1) = r + s 
                    zout(2,j,nout3) = r - s
                    r=s1 - s3
                    s=r2 + r4
                    zout(2,j,nout2) = r + s 
                    zout(2,j,nout4) = r - s
                 end do
              end do
           else
              itt=ias*before
              itrig=itt+1
              cr2=trig(1,itrig)
              ci2=trig(2,itrig)
              itrig=itrig+itt
              cr3=trig(1,itrig)
              ci3=trig(2,itrig)
              itrig=itrig+itt
              cr4=trig(1,itrig)
              ci4=trig(2,itrig)
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=r*cr2 - s*ci2
                    s2=r*ci2 + s*cr2
                    r=zin(1,j,nin3)
                    s=zin(2,j,nin3)
                    r3=r*cr3 - s*ci3
                    s3=r*ci3 + s*cr3
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=r*cr4 - s*ci4
                    s4=r*ci4 + s*cr4
                    r=r1 + r3
                    s=r2 + r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 - r3
                    s=s2 - s4
                    zout(1,j,nout2) = r - s 
                    zout(1,j,nout4) = r + s
                    r=s1 + s3
                    s=s2 + s4
                    zout(2,j,nout1) = r + s 
                    zout(2,j,nout3) = r - s
                    r=s1 - s3
                    s=r2 - r4
                    zout(2,j,nout2) = r + s 
                    zout(2,j,nout4) = r - s
                 end do
              end do
           endif
        end do
     else
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r=r1 + r3
              s=r2 + r4
              zout(1,j,nout1) = r + s
              zout(1,j,nout3) = r - s
              r=r1 - r3
              s=s2 - s4
              zout(1,j,nout2) = r + s
              zout(1,j,nout4) = r - s
              r=s1 + s3
              s=s2 + s4
              zout(2,j,nout1) = r + s
              zout(2,j,nout3) = r - s
              r=s1 - s3
              s=r2 - r4
              zout(2,j,nout2) = r - s
              zout(2,j,nout4) = r + s
           end do
        end do
        do ia=2,after
           ias=ia-1
           if (2*ias.eq.after) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(s - r)*rt2i
                    s4=(r + s)*rt2i
                    r=r1 + r3
                    s=r2 + r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 - r3
                    s=s2 + s4
                    zout(1,j,nout2) = r + s
                    zout(1,j,nout4) = r - s
                    r=s1 - s3
                    s=s2 - s4
                    zout(2,j,nout1) = r + s
                    zout(2,j,nout3) = r - s
                    r=s1 + s3
                    s=r2 - r4
                    zout(2,j,nout2) = r - s
                    zout(2,j,nout4) = r + s
                 end do
              end do
           else
              itt=ias*before
              itrig=itt+1
              cr2=trig(1,itrig)
              ci2=trig(2,itrig)
              itrig=itrig+itt
              cr3=trig(1,itrig)
              ci3=trig(2,itrig)
              itrig=itrig+itt
              cr4=trig(1,itrig)
              ci4=trig(2,itrig)
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=r*cr2 - s*ci2
                    s2=r*ci2 + s*cr2
                    r=zin(1,j,nin3)
                    s=zin(2,j,nin3)
                    r3=r*cr3 - s*ci3
                    s3=r*ci3 + s*cr3
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=r*cr4 - s*ci4
                    s4=r*ci4 + s*cr4
                    r=r1 + r3
                    s=r2 + r4
                    zout(1,j,nout1) = r + s
                    zout(1,j,nout3) = r - s
                    r=r1 - r3
                    s=s2 - s4
                    zout(1,j,nout2) = r + s
                    zout(1,j,nout4) = r - s
                    r=s1 + s3
                    s=s2 + s4
                    zout(2,j,nout1) = r + s
                    zout(2,j,nout3) = r - s
                    r=s1 - s3
                    s=r2 - r4
                    zout(2,j,nout2) = r - s
                    zout(2,j,nout4) = r + s
                 end do
              end do
           endif
        end do
     endif
! End of else if (now.eq.4)
  else if (now.eq.8) then
     if (isign.eq.-1) then 
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nin5=nin4+atb
           nin6=nin5+atb
           nin7=nin6+atb
           nin8=nin7+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           nout5=nout4+after
           nout6=nout5+after
           nout7=nout6+after
           nout8=nout7+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r5=zin(1,j,nin5)
              s5=zin(2,j,nin5)
              r6=zin(1,j,nin6)
              s6=zin(2,j,nin6)
              r7=zin(1,j,nin7)
              s7=zin(2,j,nin7)
              r8=zin(1,j,nin8)
              s8=zin(2,j,nin8)
              r=r1 + r5
              s=r3 + r7
              ap=r + s
              am=r - s
              r=r2 + r6
              s=r4 + r8
              bp=r + s
              bm=r - s
              r=s1 + s5
              s=s3 + s7
              cp=r + s
              cm=r - s
              r=s2 + s6
              s=s4 + s8
              dp1=r + s
              dm=r - s
              zout(1,j,nout1) = ap + bp
              zout(2,j,nout1) = cp + dp1
              zout(1,j,nout5) = ap - bp
              zout(2,j,nout5) = cp - dp1
              zout(1,j,nout3) = am + dm
              zout(2,j,nout3) = cm - bm
              zout(1,j,nout7) = am - dm
              zout(2,j,nout7) = cm + bm
              r=r1 - r5
              s=s3 - s7
              ap=r + s
              am=r - s
              r=s1 - s5
              s=r3 - r7
              bp=r + s
              bm=r - s
              r=s4 - s8
              s=r2 - r6
              cp=r + s
              cm=r - s
              r=s2 - s6
              s=r4 - r8
              dp1=r + s
              dm=r - s
              r = ( cp + dm)*rt2i
              s = ( dm - cp)*rt2i
              cp= ( cm + dp1)*rt2i
              dp1 = ( cm - dp1)*rt2i
              zout(1,j,nout2) = ap + r
              zout(2,j,nout2) = bm + s
              zout(1,j,nout6) = ap - r
              zout(2,j,nout6) = bm - s
              zout(1,j,nout4) = am + cp
              zout(2,j,nout4) = bp + dp1
              zout(1,j,nout8) = am - cp
              zout(2,j,nout8) = bp - dp1
           end do
        end do
        do ia=2,after
           ias=ia-1
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           itrig=itrig+itt
           cr4=trig(1,itrig)
           ci4=trig(2,itrig)
           itrig=itrig+itt
           cr5=trig(1,itrig)
           ci5=trig(2,itrig)
           itrig=itrig+itt
           cr6=trig(1,itrig)
           ci6=trig(2,itrig)
           itrig=itrig+itt
           cr7=trig(1,itrig)
           ci7=trig(2,itrig)
           itrig=itrig+itt
           cr8=trig(1,itrig)
           ci8=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb
              nin6=nin5+atb
              nin7=nin6+atb
              nin8=nin7+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              nout6=nout5+after
              nout7=nout6+after
              nout8=nout7+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=zin(1,j,nin4)
                 s=zin(2,j,nin4)
                 r4=r*cr4 - s*ci4
                 s4=r*ci4 + s*cr4
                 r=zin(1,j,nin5)
                 s=zin(2,j,nin5)
                 r5=r*cr5 - s*ci5
                 s5=r*ci5 + s*cr5
                 r=zin(1,j,nin6)
                 s=zin(2,j,nin6)
                 r6=r*cr6 - s*ci6
                 s6=r*ci6 + s*cr6
                 r=zin(1,j,nin7)
                 s=zin(2,j,nin7)
                 r7=r*cr7 - s*ci7
                 s7=r*ci7 + s*cr7
                 r=zin(1,j,nin8)
                 s=zin(2,j,nin8)
                 r8=r*cr8 - s*ci8
                 s8=r*ci8 + s*cr8
                 r=r1 + r5
                 s=r3 + r7
                 ap=r + s
                 am=r - s
                 r=r2 + r6
                 s=r4 + r8
                 bp=r + s
                 bm=r - s
                 r=s1 + s5
                 s=s3 + s7
                 cp=r + s
                 cm=r - s
                 r=s2 + s6
                 s=s4 + s8
                 dp1=r + s
                 dm=r - s
                 zout(1,j,nout1) = ap + bp
                 zout(2,j,nout1) = cp + dp1
                 zout(1,j,nout5) = ap - bp
                 zout(2,j,nout5) = cp - dp1
                 zout(1,j,nout3) = am + dm
                 zout(2,j,nout3) = cm - bm
                 zout(1,j,nout7) = am - dm
                 zout(2,j,nout7) = cm + bm
                 r=r1 - r5
                 s=s3 - s7
                 ap=r + s
                 am=r - s
                 r=s1 - s5
                 s=r3 - r7
                 bp=r + s
                 bm=r - s
                 r=s4 - s8
                 s=r2 - r6
                 cp=r + s
                 cm=r - s
                 r=s2 - s6
                 s=r4 - r8
                 dp1=r + s
                 dm=r - s
                 r = ( cp + dm)*rt2i
                 s = ( dm - cp)*rt2i
                 cp= ( cm + dp1)*rt2i
                 dp1 = ( cm - dp1)*rt2i
                 zout(1,j,nout2) = ap + r
                 zout(2,j,nout2) = bm + s
                 zout(1,j,nout6) = ap - r
                 zout(2,j,nout6) = bm - s
                 zout(1,j,nout4) = am + cp
                 zout(2,j,nout4) = bp + dp1
                 zout(1,j,nout8) = am - cp
                 zout(2,j,nout8) = bp - dp1
              end do
           end do
        end do
     else ! else for isign.eq.1
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
           nin1=nin1+after
           nin2=nin1+atb
           nin3=nin2+atb
           nin4=nin3+atb
           nin5=nin4+atb
           nin6=nin5+atb
           nin7=nin6+atb
           nin8=nin7+atb
           nout1=nout1+atn
           nout2=nout1+after
           nout3=nout2+after
           nout4=nout3+after
           nout5=nout4+after
           nout6=nout5+after
           nout7=nout6+after
           nout8=nout7+after
           do j=1,nfft
              r1=zin(1,j,nin1)
              s1=zin(2,j,nin1)
              r2=zin(1,j,nin2)
              s2=zin(2,j,nin2)
              r3=zin(1,j,nin3)
              s3=zin(2,j,nin3)
              r4=zin(1,j,nin4)
              s4=zin(2,j,nin4)
              r5=zin(1,j,nin5)
              s5=zin(2,j,nin5)
              r6=zin(1,j,nin6)
              s6=zin(2,j,nin6)
              r7=zin(1,j,nin7)
              s7=zin(2,j,nin7)
              r8=zin(1,j,nin8)
              s8=zin(2,j,nin8)
              r=r1 + r5
              s=r3 + r7
              ap=r + s
              am=r - s
              r=r2 + r6
              s=r4 + r8
              bp=r + s
              bm=r - s
              r=s1 + s5
              s=s3 + s7
              cp=r + s
              cm=r - s
              r=s2 + s6
              s=s4 + s8
              dp1=r + s
              dm=r - s
              zout(1,j,nout1) = ap + bp
              zout(2,j,nout1) = cp + dp1
              zout(1,j,nout5) = ap - bp
              zout(2,j,nout5) = cp - dp1
              zout(1,j,nout3) = am - dm
              zout(2,j,nout3) = cm + bm
              zout(1,j,nout7) = am + dm
              zout(2,j,nout7) = cm - bm
              r= r1 - r5
              s=-s3 + s7
              ap=r + s
              am=r - s
              r=s1 - s5
              s=r7 - r3
              bp=r + s
              bm=r - s
              r=-s4 + s8
              s= r2 - r6
              cp=r + s
              cm=r - s
              r=-s2 + s6
              s= r4 - r8
              dp1=r + s
              dm=r - s
              r = ( cp + dm)*rt2i
              s = ( cp - dm)*rt2i
              cp= ( cm + dp1)*rt2i
              dp1= ( dp1 - cm)*rt2i
              zout(1,j,nout2) = ap + r
              zout(2,j,nout2) = bm + s
              zout(1,j,nout6) = ap - r
              zout(2,j,nout6) = bm - s
              zout(1,j,nout4) = am + cp
              zout(2,j,nout4) = bp + dp1
              zout(1,j,nout8) = am - cp
              zout(2,j,nout8) = bp - dp1
           end do
        end do
        do ia=2,after
           ias=ia-1
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           itrig=itrig+itt
           cr4=trig(1,itrig)
           ci4=trig(2,itrig)
           itrig=itrig+itt
           cr5=trig(1,itrig)
           ci5=trig(2,itrig)
           itrig=itrig+itt
           cr6=trig(1,itrig)
           ci6=trig(2,itrig)
           itrig=itrig+itt
           cr7=trig(1,itrig)
           ci7=trig(2,itrig)
           itrig=itrig+itt
           cr8=trig(1,itrig)
           ci8=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb
              nin6=nin5+atb
              nin7=nin6+atb
              nin8=nin7+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              nout6=nout5+after
              nout7=nout6+after
              nout8=nout7+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=zin(1,j,nin4)
                 s=zin(2,j,nin4)
                 r4=r*cr4 - s*ci4
                 s4=r*ci4 + s*cr4
                 r=zin(1,j,nin5)
                 s=zin(2,j,nin5)
                 r5=r*cr5 - s*ci5
                 s5=r*ci5 + s*cr5
                 r=zin(1,j,nin6)
                 s=zin(2,j,nin6)
                 r6=r*cr6 - s*ci6
                 s6=r*ci6 + s*cr6
                 r=zin(1,j,nin7)
                 s=zin(2,j,nin7)
                 r7=r*cr7 - s*ci7
                 s7=r*ci7 + s*cr7
                 r=zin(1,j,nin8)
                 s=zin(2,j,nin8)
                 r8=r*cr8 - s*ci8
                 s8=r*ci8 + s*cr8
                 r=r1 + r5
                 s=r3 + r7
                 ap=r + s
                 am=r - s
                 r=r2 + r6
                 s=r4 + r8
                 bp=r + s
                 bm=r - s
                 r=s1 + s5
                 s=s3 + s7
                 cp=r + s
                 cm=r - s
                 r=s2 + s6
                 s=s4 + s8
                 dp1=r + s
                 dm=r - s
                 zout(1,j,nout1) = ap + bp
                 zout(2,j,nout1) = cp + dp1
                 zout(1,j,nout5) = ap - bp
                 zout(2,j,nout5) = cp - dp1
                 zout(1,j,nout3) = am - dm
                 zout(2,j,nout3) = cm + bm
                 zout(1,j,nout7) = am + dm
                 zout(2,j,nout7) = cm - bm
                 r= r1 - r5
                 s=-s3 + s7
                 ap=r + s
                 am=r - s
                 r=s1 - s5
                 s=r7 - r3
                 bp=r + s
                 bm=r - s
                 r=-s4 + s8
                 s= r2 - r6
                 cp=r + s
                 cm=r - s
                 r=-s2 + s6
                 s= r4 - r8
                 dp1=r + s
                 dm=r - s
                 r = ( cp + dm)*rt2i
                 s = ( cp - dm)*rt2i
                 cp= ( cm + dp1)*rt2i
                 dp1= ( dp1 - cm)*rt2i
                 zout(1,j,nout2) = ap + r
                 zout(2,j,nout2) = bm + s
                 zout(1,j,nout6) = ap - r
                 zout(2,j,nout6) = bm - s
                 zout(1,j,nout4) = am + cp
                 zout(2,j,nout4) = bp + dp1
                 zout(1,j,nout8) = am - cp
                 zout(2,j,nout8) = bp - dp1
              end do
           end do
        end do
     endif  !end if of isign

  else if (now.eq.3) then
!         .5d0*sqrt(3.d0)
     bb=real(isign,kind=8)*0.8660254037844387d0
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do j=1,nfft
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           r2=zin(1,j,nin2)
           s2=zin(2,j,nin2)
           r3=zin(1,j,nin3)
           s3=zin(2,j,nin3)
           r=r2 + r3
           s=s2 + s3
           zout(1,j,nout1) = r + r1
           zout(2,j,nout1) = s + s1
           r1=r1 - .5d0*r
           s1=s1 - .5d0*s
           r2=bb*(r2-r3)
           s2=bb*(s2-s3)
           zout(1,j,nout2) = r1 - s2 
           zout(2,j,nout2) = s1 + r2
           zout(1,j,nout3) = r1 + s2 
           zout(2,j,nout3) = s1 - r2
        end do
     end do
     loop_3000: do ia=2,after
        ias=ia-1
        if (4*ias.eq.3*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    r3=zin(1,j,nin3)
                    s3=zin(2,j,nin3)
                    r=r3 + r2
                    s=s2 - s3
                    zout(1,j,nout1) = r1 - r
                    zout(2,j,nout1) = s + s1
                    r1=r1 + .5d0*r
                    s1=s1 - .5d0*s        
                    r2=bb*(r2-r3)        
                    s2=bb*(s2+s3)
                    zout(1,j,nout2) = r1 - s2 
                    zout(2,j,nout2) = s1 - r2
                    zout(1,j,nout3) = r1 + s2 
                    zout(2,j,nout3) = s1 + r2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r2=zin(2,j,nin2)
                    s2=zin(1,j,nin2)
                    r3=zin(1,j,nin3)
                    s3=zin(2,j,nin3)
                    r=r2 - r3
                    s=s2 + s3
                    zout(1,j,nout1) = r + r1
                    zout(2,j,nout1) = s1 - s
                    r1=r1 - .5d0*r
                    s1=s1 + .5d0*s        
                    r2=bb*(r2+r3)        
                    s2=bb*(s2-s3)
                    zout(1,j,nout2) = r1 + s2 
                    zout(2,j,nout2) = s1 + r2
                    zout(1,j,nout3) = r1 - s2 
                    zout(2,j,nout3) = s1 - r2
                 end do
              end do
           endif
        else if (8*ias.eq.3*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r - s)*rt2i
                    s2=(r + s)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3) 
                    r=r2 - r3
                    s=s2 + s3
                    zout(1,j,nout1) = r + r1
                    zout(2,j,nout1) = s + s1
                    r1=r1 - .5d0*r
                    s1=s1 - .5d0*s        
                    r2=bb*(r2+r3)        
                    s2=bb*(s2-s3)
                    zout(1,j,nout2) = r1 - s2 
                    zout(2,j,nout2) = s1 + r2
                    zout(1,j,nout3) = r1 + s2 
                    zout(2,j,nout3) = s1 - r2
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=r2 + r3
                    s=s2 - s3
                    zout(1,j,nout1) = r + r1
                    zout(2,j,nout1) = s + s1
                    r1=r1 - .5d0*r
                    s1=s1 - .5d0*s        
                    r2=bb*(r2-r3)        
                    s2=bb*(s2+s3)
                    zout(1,j,nout2) = r1 - s2 
                    zout(2,j,nout2) = s1 + r2
                    zout(1,j,nout3) = r1 + s2 
                    zout(2,j,nout3) = s1 - r2
                 end do
              end do
           endif
        else
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=r2 + r3
                 s=s2 + s3
                 zout(1,j,nout1) = r + r1
                 zout(2,j,nout1) = s + s1
                 r1=r1 - .5d0*r
                 s1=s1 - .5d0*s
                 r2=bb*(r2-r3)
                 s2=bb*(s2-s3)
                 zout(1,j,nout2) = r1 - s2 
                 zout(2,j,nout2) = s1 + r2
                 zout(1,j,nout3) = r1 + s2 
                 zout(2,j,nout3) = s1 - r2
              end do
           end do
        endif
     end do loop_3000
! End of if (now.eq.3)

  else if (now.eq.5) then
!     cos(2.d0*pi/5.d0)
     cos2=0.3090169943749474d0
!     cos(4.d0*pi/5.d0)
     cos4=-0.8090169943749474d0
!     sin(2.d0*pi/5.d0)
     sin2=real(isign,kind=8)*0.9510565162951536d0
!      sin(4.d0*pi/5.d0)
     sin4=real(isign,kind=8)*0.5877852522924731d0
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        do j=1,nfft
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           r2=zin(1,j,nin2)
           s2=zin(2,j,nin2)
           r3=zin(1,j,nin3)
           s3=zin(2,j,nin3)
           r4=zin(1,j,nin4)
           s4=zin(2,j,nin4)
           r5=zin(1,j,nin5)
           s5=zin(2,j,nin5)
           r25 = r2 + r5
           r34 = r3 + r4
           s25 = s2 - s5
           s34 = s3 - s4
           zout(1,j,nout1) = r1 + r25 + r34
           r = r1 + cos2*r25 + cos4*r34
           s = sin2*s25 + sin4*s34
           zout(1,j,nout2) = r - s
           zout(1,j,nout5) = r + s
           r = r1 + cos4*r25 + cos2*r34
           s = sin4*s25 - sin2*s34
           zout(1,j,nout3) = r - s
           zout(1,j,nout4) = r + s
           r25 = r2 - r5
           r34 = r3 - r4
           s25 = s2 + s5
           s34 = s3 + s4
           zout(2,j,nout1) = s1 + s25 + s34
           r = s1 + cos2*s25 + cos4*s34
           s = sin2*r25 + sin4*r34
           zout(2,j,nout2) = r + s
           zout(2,j,nout5) = r - s
           r = s1 + cos4*s25 + cos2*s34
           s = sin4*r25 - sin2*r34
           zout(2,j,nout3) = r + s
           zout(2,j,nout4) = r - s
        end do
     end do
     loop_5000: do ia=2,after
        ias=ia-1
        if (8*ias.eq.5*after) then
           if (isign.eq.1) then
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nin5=nin4+atb        
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 nout5=nout4+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r - s)*rt2i
                    s2=(r + s)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3) 
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(r + s)*rt2i
                    s4=(r - s)*rt2i
                    r5=zin(1,j,nin5)
                    s5=zin(2,j,nin5)
                    r25 = r2 - r5
                    r34 = r3 + r4
                    s25 = s2 + s5
                    s34 = s3 - s4
                    zout(1,j,nout1) = r1 + r25 - r34
                    r = r1 + cos2*r25 - cos4*r34 
                    s = sin2*s25 + sin4*s34
                    zout(1,j,nout2) = r - s
                    zout(1,j,nout5) = r + s
                    r = r1 + cos4*r25 - cos2*r34 
                    s = sin4*s25 - sin2*s34
                    zout(1,j,nout3) = r - s
                    zout(1,j,nout4) = r + s
                    r25 = r2 + r5
                    r34 = r4 - r3
                    s25 = s2 - s5
                    s34 = s3 + s4
                    zout(2,j,nout1) = s1 + s25 + s34
                    r = s1 + cos2*s25 + cos4*s34
                    s = sin2*r25 + sin4*r34
                    zout(2,j,nout2) = r + s
                    zout(2,j,nout5) = r - s
                    r = s1 + cos4*s25 + cos2*s34
                    s = sin4*r25 - sin2*r34
                    zout(2,j,nout3) = r + s
                    zout(2,j,nout4) = r - s
                 end do
              end do
           else
              nin1=ia-after
              nout1=ia-atn
              do ib=1,before
                 nin1=nin1+after
                 nin2=nin1+atb
                 nin3=nin2+atb
                 nin4=nin3+atb
                 nin5=nin4+atb        
                 nout1=nout1+atn
                 nout2=nout1+after
                 nout3=nout2+after
                 nout4=nout3+after
                 nout5=nout4+after
                 do j=1,nfft
                    r1=zin(1,j,nin1)
                    s1=zin(2,j,nin1)
                    r=zin(1,j,nin2)
                    s=zin(2,j,nin2)
                    r2=(r + s)*rt2i
                    s2=(s - r)*rt2i
                    r3=zin(2,j,nin3)
                    s3=zin(1,j,nin3)
                    r=zin(1,j,nin4)
                    s=zin(2,j,nin4)
                    r4=(s - r)*rt2i
                    s4=(r + s)*rt2i
                    r5=zin(1,j,nin5)
                    s5=zin(2,j,nin5)
                    r25 = r2 - r5
                    r34 = r3 + r4
                    s25 = s2 + s5
                    s34 = s4 - s3
                    zout(1,j,nout1) = r1 + r25 + r34
                    r = r1 + cos2*r25 + cos4*r34
                    s = sin2*s25 + sin4*s34
                    zout(1,j,nout2) = r - s
                    zout(1,j,nout5) = r + s
                    r = r1 + cos4*r25 + cos2*r34
                    s = sin4*s25 - sin2*s34
                    zout(1,j,nout3) = r - s
                    zout(1,j,nout4) = r + s
                    r25 = r2 + r5
                    r34 = r3 - r4
                    s25 = s2 - s5
                    s34 = s3 + s4
                    zout(2,j,nout1) = s1 + s25 - s34
                    r = s1 + cos2*s25 - cos4*s34
                    s = sin2*r25 + sin4*r34
                    zout(2,j,nout2) = r + s
                    zout(2,j,nout5) = r - s
                    r = s1 + cos4*s25 - cos2*s34
                    s = sin4*r25 - sin2*r34
                    zout(2,j,nout3) = r + s
                    zout(2,j,nout4) = r - s
                 end do
              end do
           endif
        else !if of (ias...
           ias=ia-1
           itt=ias*before
           itrig=itt+1
           cr2=trig(1,itrig)
           ci2=trig(2,itrig)
           itrig=itrig+itt
           cr3=trig(1,itrig)
           ci3=trig(2,itrig)
           itrig=itrig+itt
           cr4=trig(1,itrig)
           ci4=trig(2,itrig)
           itrig=itrig+itt
           cr5=trig(1,itrig)
           ci5=trig(2,itrig)
           nin1=ia-after
           nout1=ia-atn
           do ib=1,before
              nin1=nin1+after
              nin2=nin1+atb
              nin3=nin2+atb
              nin4=nin3+atb
              nin5=nin4+atb
              nout1=nout1+atn
              nout2=nout1+after
              nout3=nout2+after
              nout4=nout3+after
              nout5=nout4+after
              do j=1,nfft
                 r1=zin(1,j,nin1)
                 s1=zin(2,j,nin1)
                 r=zin(1,j,nin2)
                 s=zin(2,j,nin2)
                 r2=r*cr2 - s*ci2
                 s2=r*ci2 + s*cr2
                 r=zin(1,j,nin3)
                 s=zin(2,j,nin3)
                 r3=r*cr3 - s*ci3
                 s3=r*ci3 + s*cr3
                 r=zin(1,j,nin4)
                 s=zin(2,j,nin4)
                 r4=r*cr4 - s*ci4
                 s4=r*ci4 + s*cr4
                 r=zin(1,j,nin5)
                 s=zin(2,j,nin5)
                 r5=r*cr5 - s*ci5
                 s5=r*ci5 + s*cr5
                 r25 = r2 + r5
                 r34 = r3 + r4
                 s25 = s2 - s5
                 s34 = s3 - s4
                 zout(1,j,nout1) = r1 + r25 + r34
                 r = r1 + cos2*r25 + cos4*r34
                 s = sin2*s25 + sin4*s34
                 zout(1,j,nout2) = r - s
                 zout(1,j,nout5) = r + s
                 r = r1 + cos4*r25 + cos2*r34
                 s = sin4*s25 - sin2*s34
                 zout(1,j,nout3) = r - s
                 zout(1,j,nout4) = r + s
                 r25 = r2 - r5
                 r34 = r3 - r4
                 s25 = s2 + s5
                 s34 = s3 + s4
                 zout(2,j,nout1) = s1 + s25 + s34
                 r = s1 + cos2*s25 + cos4*s34
                 s = sin2*r25 + sin4*r34
                 zout(2,j,nout2) = r + s
                 zout(2,j,nout5) = r - s
                 r = s1 + cos4*s25 + cos2*s34
                 s = sin4*r25 - sin2*r34
                 zout(2,j,nout3) = r + s
                 zout(2,j,nout4) = r - s
              end do
           end do
        endif
     end do loop_5000
! end of if now.eq.5

  else if (now.eq.6) then
!     .5d0*sqrt(3.d0)
     bb=real(isign,kind=8)*0.8660254037844387d0
     ia=1
     nin1=ia-after
     nout1=ia-atn
     do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        do j=1,nfft
           r2=zin(1,j,nin3)
           s2=zin(2,j,nin3)
           r3=zin(1,j,nin5)
           s3=zin(2,j,nin5)
           r=r2 + r3
           s=s2 + s3
           r1=zin(1,j,nin1)
           s1=zin(2,j,nin1)
           ur1 = r + r1
           ui1 = s + s1
           r1=r1 - .5d0*r
           s1=s1 - .5d0*s
           r=r2-r3
           s=s2-s3
           ur2 = r1 - s*bb
           ui2 = s1 + r*bb
           ur3 = r1 + s*bb
           ui3 = s1 - r*bb

           r2=zin(1,j,nin6)
           s2=zin(2,j,nin6)
           r3=zin(1,j,nin2)
           s3=zin(2,j,nin2)
           r=r2 + r3
           s=s2 + s3
           r1=zin(1,j,nin4)
           s1=zin(2,j,nin4)
           vr1 = r + r1
           vi1 = s + s1
           r1=r1 - .5d0*r
           s1=s1 - .5d0*s
           r=r2-r3
           s=s2-s3
           vr2 = r1 - s*bb
           vi2 = s1 + r*bb
           vr3 = r1 + s*bb
           vi3 = s1 - r*bb

           zout(1,j,nout1)=ur1+vr1
           zout(2,j,nout1)=ui1+vi1
           zout(1,j,nout5)=ur2+vr2
           zout(2,j,nout5)=ui2+vi2
           zout(1,j,nout3)=ur3+vr3
           zout(2,j,nout3)=ui3+vi3
           zout(1,j,nout4)=ur1-vr1
           zout(2,j,nout4)=ui1-vi1
           zout(1,j,nout2)=ur2-vr2
           zout(2,j,nout2)=ui2-vi2
           zout(1,j,nout6)=ur3-vr3
           zout(2,j,nout6)=ui3-vi3
        end do
     end do

  else if (now.eq.7) then

     ia=1
     nin1=ia-after
     nout1=ia-atn

     do ib=1,before
      	nin1=nin1+after
      	nin2=nin1+atb
      	nin3=nin2+atb
      	nin4=nin3+atb
      	nin5=nin4+atb
      	nin6=nin5+atb
      	nin7=nin6+atb
      	nout1=nout1+atn
      	nout2=nout1+after
      	nout3=nout2+after
      	nout4=nout3+after
      	nout5=nout4+after
      	nout6=nout5+after
      	nout7=nout6+after
      	do j=1,nfft

           cos2= 0.62348980185873353053d0
           cos4=-0.22252093395631440429d0
           cos6=-0.90096886790241912624d0
           
           sin2= 0.78183148246802980871d0*isign
           sin4= 0.97492791218182360702d0*isign
           sin6= 0.43388373911755812048d0*isign
           
           r1=zin(1,j,nin1);   s1=zin(2,j,nin1)
           r2=zin(1,j,nin2);   s2=zin(2,j,nin2)
           r3=zin(1,j,nin3);   s3=zin(2,j,nin3)
           r4=zin(1,j,nin4);   s4=zin(2,j,nin4)
           r5=zin(1,j,nin5);   s5=zin(2,j,nin5)
           r6=zin(1,j,nin6);   s6=zin(2,j,nin6)
           r7=zin(1,j,nin7);   s7=zin(2,j,nin7)
           
           rp2=r2+r7; rp3=r3+r6; rp4=r4+r5
           rm2=r2-r7; rm3=r3-r6; rm4=r4-r5
           sp2=s2+s7; sp3=s3+s6; sp4=s4+s5
           sm2=s2-s7; sm3=s3-s6; sm4=s4-s5 !3*4=12 flops

           ur2=r1+cos2*rp2+cos4*rp3+cos6*rp4
           vr2=   sin2*sm2+sin4*sm3+sin6*sm4
           
           ui2=s1+cos2*sp2+cos4*sp3+cos6*sp4
           vi2=   sin2*rm2+sin4*rm3+sin6*rm4
                     
           ur3=r1+cos4*rp2+cos6*rp3+cos2*rp4
           vr3=   sin4*sm2-sin6*sm3-sin2*sm4
           
           ui3=s1+cos4*sp2+cos6*sp3+cos2*sp4
           vi3=   sin4*rm2-sin6*rm3-sin2*rm4
           
           ur4=r1+cos6*rp2+cos2*rp3+cos4*rp4
           vr4=   sin6*sm2-sin2*sm3+sin4*sm4
           
           ui4=s1+cos6*sp2+cos2*sp3+cos4*sp4
           vi4=   sin6*rm2-sin2*rm3+sin4*rm4   ! 11*6=66 flops, of them
           !6*6=36 multiplications
           
           zout(1,j,nout1)=r1+rp2+rp3+rp4; zout(2,j,nout1)=s1+sp2+sp3+sp4
           zout(1,j,nout2)=ur2-vr2;        zout(2,j,nout2)=ui2+vi2
           zout(1,j,nout3)=ur3-vr3;        zout(2,j,nout3)=ui3+vi3
           zout(1,j,nout4)=ur4-vr4;        zout(2,j,nout4)=ui4+vi4
           zout(1,j,nout5)=ur4+vr4;        zout(2,j,nout5)=ui4-vi4
           zout(1,j,nout6)=ur3+vr3;        zout(2,j,nout6)=ui3-vi3
           zout(1,j,nout7)=ur2+vr2;        zout(2,j,nout7)=ui2-vi2 
           !(3+6)*2 =18 flops

           ! In total: 12+66+18=96flops

        end do
     end do
  else 
     stop 'error fftstp'
  endif !end of now

END SUBROUTINE fftstp

