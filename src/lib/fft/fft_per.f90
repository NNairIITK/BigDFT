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
!
!
! Presumably good compiler flags:
! IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
! with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ; 
!                   a.out
! DEC: f90 -O3 -arch ev67 -pipeline
! with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ; 
!                   prun -N1 -c4 a.out
!-----------------------------------------------------------


!!****f* BigDFT/fourier_dim
!! FUNCTION
!!   Give a number n_next > n compatible for the FFT
!! COPYRIGHT
!!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!  Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
!!  Copyright (C) 2002-2009 BigDFT group 
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~/COPYING file
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine fourier_dim(n,n_next)
  implicit none
  !Arguments
  integer, intent(in) :: n
  integer, intent(out) :: n_next

  !Local variables
  integer, parameter :: ndata = 180
  !Multiple of 2,3,5
  integer, dimension(ndata), parameter :: idata = (/   &
          3,     4,     5,     6,     8,     9,    12,    15,    16,    18, &
         20,    24,    25,    27,    30,    32,    36,    40,    45,    48, &
         54,    60,    64,    72,    75,    80,    81,    90,    96,   100, &
        108,   120,   125,   128,   135,   144,   150,   160,   162,   180, &
        192,   200,   216,   225,   240,   243,   256,   270,   288,   300, &
        320,   324,   360,   375,   384,   400,   405,   432,   450,   480, &
        486,   500,   512,   540,   576,   600,   625,   640,   648,   675, &
        720,   729,   750,   768,   800,   810,   864,   900,   960,   972, &
       1000,  1024,  1080,  1125,  1152,  1200,  1215,  1280,  1296,  1350, &
       1440,  1458,  1500,  1536,  1600,  1620,  1728,  1800,  1875,  1920, &
       1944,  2000,  2025,  2048,  2160,  2250,  2304,  2400,  2430,  2500, &
       2560,  2592,  2700,  2880,  3000,  3072,  3125,  3200,  3240,  3375, &
       3456,  3600,  3750,  3840,  3888,  4000,  4050,  4096,  4320,  4500, &
       4608,  4800,  5000,  5120,  5184,  5400,  5625,  5760,  6000,  6144, &
       6400,  6480,  6750,  6912,  7200,  7500,  7680,  8000,  8192,  8640, &
       9000,  9216,  9375,  9600, 10000, 10240, 10368, 10800, 11250, 11520, &
      12000, 12288, 12500, 12800, 13824, 14400, 15000, 15360, 15625, 16000, &
      16384, 17280, 18000, 18432, 18750, 19200, 20000, 20480, 23040, 24000   /)                                              
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


!!****f* BigDFT/dimensions_fft
!! FUNCTION
!!   Give the dimensions of fft arrays
!! SOURCE
!!
subroutine dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)
    implicit none
    integer, intent(in)::n1,n2,n3
    integer,intent(out)::nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
! and the same for n2,n3.
    nd1=n1+2;        nd2=n2+2;    nd3=n3+2
    ! n1b>=n1f;   n3f>=n3b
    n1f=(n1+2)/2 
    n3f=(n3+1)/2+1
    
    n1b=(n1+1)/2+1
    n3b=(n3+2)/2
    
    nd1f=n1f+1
    nd3f=n3f+1
    
    nd1b=n1b+1
    nd3b=n3b+1
end subroutine dimensions_fft
!!***


!!****f* BigDFT/FFT
!! FUNCTION
!!  Calculates the discrete fourier transform F(I1,I2,I3)=
!!  S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
!! DESCRIPTION
!!  with optimal performance on vector computer, workstations and 
!!  multiprocessor shared memory computers using OpenMP compiler directives
!!   INPUT:
!!       n1,n2,n3:physical dimension of the transform. It must be a 
!!                product of the prime factors 2,3,5, but greater than 3. 
!!               If two ni's are equal it is recommended to place them 
!!               behind each other.
!!       nd1,nd2,nd3:memory dimension of Z. ndi must always be greater or 
!!                   equal than ni. On a vector machine, it is recomended 
!!                  to chose ndi=ni if ni is odd and ndi=ni+1 if ni is 
!!                  even to obtain optimal execution speed. On RISC 
!!                  machines ndi=ni is usually fine for odd ni, for even 
!!                  ni one should try ndi=ni+1, ni+2, ni+4 to find the 
!!                  optimal performance. 
!!      inzee=1: first part of Z is data (input) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           Z(1,i1,i2,i3,inzee)=real(R(i1,i2,i3))
!!           Z(2,i1,i2,i3,inzee)=imag(R(i1,i2,i3))
!!   OUTPUT:
!!      inzee=1: first part of Z is data (output) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!           imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!      inzee on output is in general different from inzee on input
!!   The input data are always overwritten independently of the 
!!  value of inzee.
!! 
!! PERFORMANCE AND THE NCACHE
!!  The most important feature for performance is the right choice of 
!!  the parameter ncache. On a vector machine ncache has to be put to 0.
!!  On a RISC machine with cache, it is very important to find the optimal 
!!  value of NCACHE. NCACHE determines the size of the work array zw, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly 
!!  half the size of the physical cache in units of real*8 numbers.
!!  If the machine has 2 cache levels it can not be predicted which 
!!  cache level will be the most relevant one for choosing ncache. 
!!  The optimal value of ncache can easily be determined by numerical 
!!  experimentation. A too large value of ncache leads to a dramatic 
!!  and sudden decrease of performance, a too small value to a to a 
!!  slow and less dramatic decrease of performance. If NCACHE is set 
!!  to a value so small, that not even a single one dimensional transform 
!!  can be done in the workarray zw, the program stops with an error 
!!  message.
!! SOURCE
!!
subroutine FFT(n1,n2,n3,nd1,nd2,nd3,z,isign,inzee)
        implicit real*8 (a-h,o-z)
!!!$      interface
!!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
!!!!$        end function omp_get_num_threads
!!!$      end interface
!!!!$      interface
!!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
!!!!$        end function omp_get_thread_num
!!!!$      end interface


        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: zw  
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: trig
        INTEGER, ALLOCATABLE, DIMENSION(:) :: after,now,before
        dimension z(2,nd1*nd2*nd3,2)
        if (max(n1,n2,n3).gt.1024) stop '1024'

! some reasonable values of ncache: 
!   IBM/RS6000/590: 16*1024 ; IBM/RS6000/390: 3*1024 ; 
!   IBM/PwPC: 1*1024 ; SGI/MIPS/R8000: 16*1024 ; DEC/Alpha/EV5 and EV6 6*1024
!   But if you care about performance find the optimal value of ncache yourself!
!       On all vector machines: ncache=0

        ncache=6*1024

! check whether input values are reasonable
    if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
    if (isign.ne.1 .and. isign.ne.-1) stop 'wrong isign'
    if (n1.gt.nd1) stop 'n1>nd1'
    if (n2.gt.nd2) stop 'n2>nd2'
    if (n3.gt.nd3) stop 'n3>nd3'
    

! vector computer with memory banks:
      if (ncache.eq.0) then
        allocate(trig(2,1024),after(20),now(20),before(20))

        call ctrig_(n3,trig,after,before,now,isign,ic)
        nfft=nd1*n2
        mm=nd1*nd2
        do 51093,i=1,ic-1
        call fftstp_(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
                          trig,after(i),now(i),before(i),isign)
51093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
                          trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

        if (n2.ne.n3) call ctrig_(n2,trig,after,before,now,isign,ic)
        nfft=nd3*n1
        mm=nd3*nd1
        do 52093,i=1,ic-1
        call fftstp_(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee), &
                           trig,after(i),now(i),before(i),isign)
52093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee), &
                       trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

        if (n1.ne.n2) call ctrig_(n1,trig,after,before,now,isign,ic)
        nfft=nd2*n3
        mm=nd2*nd3
        do 53093,i=1,ic-1
        call fftstp_(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee), &
                         trig,after(i),now(i),before(i),isign)
53093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee), &
                         trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

! RISC machine with cache:
      else
! INtel IFC does not understand default(private)
!!!!!$omp parallel  default(private) &
!!!!$omp parallel & 
!!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
!!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,isign,inzee,ncache) 
        npr=1
!!!!$       npr=omp_get_num_threads()
        iam=0
!!!!$       iam=omp_get_thread_num()
!        write(6,*) 'npr,iam',npr,iam
! Critical section only necessary on Intel
!!!!$omp critical
        allocate(zw(2,ncache/4,2),trig(2,1024),after(20),now(20),before(20))
!!!!$omp end critical

        inzet=inzee
! TRANSFORM ALONG Z AXIS

        mm=nd1*nd2
        m=nd3
        lot=max(1,ncache/(4*n3))
        nn=lot
        n=n3
        if (2*n*lot*2.gt.ncache) stop 'ncache1'

        call ctrig_(n3,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd1*n2)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd1*n2)
        nfft=mb-ma+1
        j=ma
        jj=j*nd3-nd3+1
        call fftrot(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
                          trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd1*n2)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd1*n2)
        do 1000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd3-nd3+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 1093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
1093        inzeep=3-inzeep
        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
1000        continue
      endif

        inzet=3-inzet

!!!!!!!!!$omp barrier

! TRANSFORM ALONG Y AXIS
        mm=nd3*nd1
        m=nd2
        lot=max(1,ncache/(4*n2))
        nn=lot
        n=n2
        if (2*n*lot*2.gt.ncache) stop 'ncache2'

        if (n2.ne.n3) call ctrig_(n2,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd3*n1)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd3*n1)
        nfft=mb-ma+1
        j=ma
        jj=j*nd2-nd2+1
        call fftrot(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd3*n1)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd3*n1)
        do 2000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd2-nd2+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 2093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
2093        inzeep=3-inzeep

        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
2000        continue
      endif
        inzet=3-inzet

!!!!!!!!$omp barrier

! TRANSFORM ALONG X AXIS
        mm=nd2*nd3
        m=nd1
        lot=max(1,ncache/(4*n1))
        nn=lot
        n=n1
        if (2*n*lot*2.gt.ncache) stop 'ncache3'

        if (n1.ne.n2) call ctrig_(n1,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd2*n3)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd2*n3)
        nfft=mb-ma+1
        j=ma
        jj=j*nd1-nd1+1
        call fftrot(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd2*n3)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd2*n3)
        do 3000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd1-nd1+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 3093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
3093        inzeep=3-inzeep
        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
3000        continue
      endif
        inzet=3-inzet
        
        deallocate(zw,trig,after,now,before)
        if (iam.eq.0) inzee=inzet
!!!!!!!!!!$omp end parallel  

      endif
      return
end subroutine FFT
!!***


!!***f* BigDFT/FFT_for
!! DESCRIPTION
!!  Calculates the discrete Fourier transform F(I1,I2,I3)=
!!  S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
!!  with optimal performance on vector computer, workstations and 
!!  multiprocessor shared memory computers using OpenMP compiler directives
!!   INPUT:
!!       n1,n2,n3:physical dimension of the transform. It must be a 
!!                product of the prime factors 2,3,5, but greater than 3. 
!!               If two ni's are equal it is recommended to place them 
!!               behind each other.
!!       nd1,nd2,nd3:memory dimension of Z. ndi must always be greater or 
!!                   equal than ni. On a vector machine, it is recomended 
!!                  to chose ndi=ni if ni is odd and ndi=ni+1 if ni is 
!!                  even to obtain optimal execution speed. On RISC 
!!                  machines ndi=ni is usually fine for odd ni, for even 
!!                  ni one should try ndi=ni+1, ni+2, ni+4 to find the 
!!                  optimal performance. 
!!      inzee=1: first part of Z is data (input) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           Z(1,i1,i2,i3,inzee)=real(R(i1,i2,i3))
!!           Z(2,i1,i2,i3,inzee)=imag(R(i1,i2,i3))
!!   OUTPUT:
!!      inzee=1: first part of Z is data (output) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!           imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!      inzee on output is in general different from inzee on input
!!   The input data are always overwritten independently of the 
!!  value of inzee.
!! SOURCE
!!
subroutine FFT_for(n1,n2,n3,n1f,n3f,nd1,nd2,nd3,nd1f,nd3f,x0,z1,z3,inzee)
        implicit real*8 (a-h,o-z)
!!!!!!!$      interface
!!!!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
!!!!!!$        end function omp_get_num_threads
!!!!!!$      end interface
!!!!!!!$      interface
!!!!!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
!!!!!!!!$        end function omp_get_thread_num
!!!!!!!!$      end interface

        integer,intent(in)::nd1,nd2,nd3,nd1f,nd3f
        integer,intent(in)::n1,n2,n3,n1f,n3f
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: zw  
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: trig
        INTEGER, ALLOCATABLE, DIMENSION(:) :: after,now,before

        real*8,intent(in):: x0(n1,n2,n3)
        real*8,intent(out)::z3(2,nd1*nd2*nd3f,2)
        real*8              ::z1(2,nd1f*nd2*nd3,2) ! work array

        integer::mm,nffta,isign=1
        if (max(n1,n2,n3).gt.1024) stop '1024'

! some reasonable values of ncache: 
!   IBM/RS6000/590: 16*1024 ; IBM/RS6000/390: 3*1024 ; 
!   IBM/PwPC: 1*1024 ; SGI/MIPS/R8000: 16*1024 ; DEC/Alpha/EV5 and EV6 6*1024
!   But if you care about performance find the optimal value of ncache yourself!
!       On all vector machines: ncache=0

        inzee=1
        ncache=6*1024
!*******************Alexey*********************************************************************
!         ncache=0
!**********************************************************************************************

! check whether input values are reasonable
    if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
    if (isign.ne.1 .and. isign.ne.-1) stop 'wrong isign'
    if (n1.gt.nd1) stop 'n1>nd1'
    if (n2.gt.nd2) stop 'n2>nd2'
    if (n3.gt.nd3) stop 'n3>nd3'
    
        
!        call x0_to_z1_simple(x0,z1,inzee)
        call x0_to_z1(x0,z1,inzee)

      if (ncache.eq.0) then
! vector computer with memory banks:
        allocate(trig(2,1024),after(20),now(20),before(20))

        call ctrig_(n3,trig,after,before,now,isign,ic)

        mm=nd1f*nd2 
        nffta=nd1f*n2

        do 51093,i=1,ic-1
        call fftstp_(mm,nffta,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
                          trig,after(i),now(i),before(i),isign)
51093            inzee=3-inzee
        i=ic
        call fftrot(mm,nffta,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
                          trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

        call z1_to_z3(z1,z3,inzee)
!===============================================================================================
        if (n2.ne.n3) call ctrig_(n2,trig,after,before,now,isign,ic)
        nfft=nd3f*n1
        mm=nd3f*nd1
        do 52093,i=1,ic-1
        call fftstp_(mm,nfft,nd2,mm,nd2,z3(1,1,inzee),z3(1,1,3-inzee), &
                           trig,after(i),now(i),before(i),isign)
52093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd2,mm,nd2,z3(1,1,inzee),z3(1,1,3-inzee), &
                       trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

        if (n1.ne.n2) call ctrig_(n1,trig,after,before,now,isign,ic)
        nfft=nd2*n3f
        mm=nd2*nd3f
        do 53093,i=1,ic-1
        call fftstp_(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
                         trig,after(i),now(i),before(i),isign)
53093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
                         trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

! RISC machine with cache:
      else
! INtel IFC does not understand default(private)
!!!!$omp parallel  default(private) &
!!!!$omp parallel & 
!!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
!!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,isign,inzee,ncache) 
        npr=1
!!!!!!!$       npr=omp_get_num_threads()
        iam=0
!!!!!!!$       iam=omp_get_thread_num()
!!!!!        write(6,*) 'npr,iam',npr,iam
! Critical section only necessary on Intel
!!!!!$omp critical
        allocate(zw(2,ncache/4,2),trig(2,1024),after(20),now(20),before(20))
!!!!!$omp end critical

        inzet=inzee
! TRANSFORM ALONG Z AXIS

        mm=nd1f*nd2
        m=nd3
        lot=max(1,ncache/(4*n3))
        nn=lot
        n=n3
        if (2*n*lot*2.gt.ncache) stop 'ncache1'

        call ctrig_(n3,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd1f*n2)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd1f*n2)
        nfft=mb-ma+1
        j=ma
        jj=j*nd3-nd3+1
        call fftrot(mm,nfft,m,mm,m,z1(1,j,inzet),z1(1,jj,3-inzet), &
                          trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd1f*n2)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd1f*n2)
        do 1000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd3-nd3+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z1(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 1093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
1093        inzeep=3-inzeep
        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z1(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
1000        continue
      endif

        inzet=3-inzet


        call z1_to_z3(z1,z3,inzet)
!!!!!!!!!$omp barrier

! TRANSFORM ALONG Y AXIS
        mm=nd3f*nd1
        m=nd2
        lot=max(1,ncache/(4*n2))
        nn=lot
        n=n2
        if (2*n*lot*2.gt.ncache) stop 'ncache2'

        if (n2.ne.n3) call ctrig_(n2,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd3f*n1)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd3f*n1)
        nfft=mb-ma+1
        j=ma
        jj=j*nd2-nd2+1
        call fftrot(mm,nfft,m,mm,m,z3(1,j,inzet),z3(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd3f*n1)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd3f*n1)
        do 2000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd2-nd2+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z3(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 2093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
2093        inzeep=3-inzeep

        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z3(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
2000        continue
      endif
        inzet=3-inzet

!!!!!!!!!$omp barrier

! TRANSFORM ALONG X AXIS
        mm=nd2*nd3f
        m=nd1
        lot=max(1,ncache/(4*n1))
        nn=lot
        n=n1
        if (2*n*lot*2.gt.ncache) stop 'ncache3'

        if (n1.ne.n2) call ctrig_(n1,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd2*n3f)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd2*n3f)
        nfft=mb-ma+1
        j=ma
        jj=j*nd1-nd1+1
        call fftrot(mm,nfft,m,mm,m,z3(1,j,inzet),z3(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd2*n3f)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd2*n3f)
        do 3000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd1-nd1+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z3(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 3093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
3093        inzeep=3-inzeep
        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z3(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
3000        continue
      endif
        inzet=3-inzet
        
        deallocate(zw,trig,after,now,before)
        if (iam.eq.0) inzee=inzet
!!!!!!!!!$omp end parallel  


      endif

contains

   subroutine x0_to_z1(x0,z1,inzee)
            ! Transform the real array x0 into a complex z1
            ! real      part of z1: elements of x0 with odd  i1
            ! imaginary part of z1: elements of x0 with even i1
            implicit none
            integer,intent(in)::inzee
            real*8,intent(in):: x0(n1,n2,n3)
            real*8,intent(out)::z1(2,nd1f,nd2,nd3,2)
            integer i2,i3

            do i3=1,n3
                do i2=1,n2
                    ! 2*n1f=n1 for even n1
                    ! 2*n1f=n1+1 for odd n1. Then, we copy one more element than
                    ! necessary, but that's no problem.
                    call my_copy(z1(1,1,i2,i3,inzee),x0(1,i2,i3))
                enddo
            enddo

            end subroutine x0_to_z1

            subroutine my_copy(x,y)
                ! copies complex array y into complex array x
                implicit none
                real*8 x(2,n1f),y(2,n1f)
                x=y
            end subroutine my_copy

            subroutine x0_to_z1_simple(x0,z1,inzee)
            ! Transform the real array x0 into a complex z1
            ! real      part of z1: elements of x0 with odd  i1
            ! imaginary part of z1: elements of x0 with even i1
            implicit none
            integer,intent(in)::inzee
            real*8,intent(in):: x0(n1,n2,n3)
            real*8,intent(out)::z1(2,nd1f,nd2,nd3,2)
            integer i1,i2,i3
            
            if (n1f*2.eq.n1) then
                do i3=1,n3
                    do i2=1,n2
                        do i1=1,n1f
                            z1(1,i1,i2,i3,inzee)=x0(2*i1-1,i2,i3)
                            z1(2,i1,i2,i3,inzee)=x0(2*i1  ,i2,i3)
                        enddo
                    enddo
                enddo
            else ! n1=2*n1f-1
                do i3=1,n3
                    do i2=1,n2
                        do i1=1,n1f-1
                            z1(1,i1,i2,i3,inzee)=x0(2*i1-1,i2,i3)
                            z1(2,i1,i2,i3,inzee)=x0(2*i1  ,i2,i3)
                        enddo
                        z1(1,n1f,i2,i3,inzee)=x0(n1,i2,i3)
                    enddo
                enddo
            endif

    end subroutine x0_to_z1_simple


    subroutine z1_to_z3(z1,z3,inzee)
    ! transforms the array z1 that stores elements of z corresponding to even 
    ! and odd values of i1, as symmetric and antisymmetric combinations w.r.t.
    ! flip of i3,
    ! into the array z3 that stores only elements of z with i3=<nd3f
    implicit none
    integer,intent(in)::inzee
    integer i1,i2,i3
    real*8,intent(in):: z1(2,nd3,nd1f,nd2,2)
    real*8,intent(out)::z3( 2,nd3f,nd1 ,nd2,2)
    
    if (n1f*2.eq.n1) then
        ! i3=1
        do i2=1,n2
            do i1=1,n1f
                z3(1,1,2*i1-1,i2,inzee)= 2.d0*z1(1,1,i1,i2,inzee)
                z3(2,1,2*i1-1,i2,inzee)= 0.d0
                z3(1,1,2*i1  ,i2,inzee)= 2.d0*z1(2,1,i1,i2,inzee)
                z3(2,1,2*i1  ,i2,inzee)= 0.d0
            enddo
        enddo
        
        do i2=1,n2
            do i1=1,n1f
                do i3=2,n3f
                    z3(1,i3,2*i1-1,i2,inzee)= z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                    z3(2,i3,2*i1-1,i2,inzee)= z1(2,i3,i1,i2,inzee)-z1(2,n3+2-i3,i1,i2,inzee)
                    z3(1,i3,2*i1  ,i2,inzee)= z1(2,i3,i1,i2,inzee)+z1(2,n3+2-i3,i1,i2,inzee)
                    z3(2,i3,2*i1  ,i2,inzee)=-z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                enddo
            enddo
        enddo
    else ! n1=2*n1f-1
        ! i3=1
        do i2=1,n2
            do i1=1,n1f-1
                z3(1,1,2*i1-1,i2,inzee)= 2.d0*z1(1,1,i1,i2,inzee)
                z3(2,1,2*i1-1,i2,inzee)= 0.d0
                z3(1,1,2*i1  ,i2,inzee)= 2.d0*z1(2,1,i1,i2,inzee)
                z3(2,1,2*i1  ,i2,inzee)= 0.d0
            enddo
        enddo
        
        do i2=1,n2
            do i1=1,n1f-1
                do i3=2,n3f
                    z3(1,i3,2*i1-1,i2,inzee)= z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                    z3(2,i3,2*i1-1,i2,inzee)= z1(2,i3,i1,i2,inzee)-z1(2,n3+2-i3,i1,i2,inzee)
                    z3(1,i3,2*i1  ,i2,inzee)= z1(2,i3,i1,i2,inzee)+z1(2,n3+2-i3,i1,i2,inzee)
                    z3(2,i3,2*i1  ,i2,inzee)=-z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                enddo
            enddo
        enddo

        ! i1=n1f is treated separately: 2*n1f-1=n1, but terms with 2*n1f are
        ! omitted

        do i2=1,n2
            z3(1,1,n1,i2,inzee)= 2.d0*z1(1,1,n1f,i2,inzee)
            z3(2,1,n1,i2,inzee)= 0.d0
        enddo
        
        do i2=1,n2
            do i3=2,n3f
                z3(1,i3,n1,i2,inzee)= z1(1,i3,n1f,i2,inzee)+z1(1,n3+2-i3,n1f,i2,inzee)
                z3(2,i3,n1,i2,inzee)= z1(2,i3,n1f,i2,inzee)-z1(2,n3+2-i3,n1f,i2,inzee)
            enddo
        enddo
    endif

    end subroutine z1_to_z3


end subroutine fft_for
!!***


!!****f* BigDFT/FFT_back
!! FUNCTION
!!  Calculates the discrete fourier transform F(I1,I2,I3)=
!!  S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
subroutine FFT_back(n1,n2,n3,n1f,n1b,n3f,n3b,nd1,nd2,nd3,nd1f,nd1b,nd3f,nd3b,y,z1,z3,inzee)
        implicit real*8 (a-h,o-z)
!!!!!!!$      interface
!!!!!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
!!!!!!$        end function omp_get_num_threads
!!!!!!$      end interface
!!!!!!!$      interface
!!!!!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
!!!!!!!$        end function omp_get_thread_num
!!!!!!!!$      end interface


        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: zw  
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: trig
        INTEGER, ALLOCATABLE, DIMENSION(:) :: after,now,before
        integer,intent(in):: nd3f,n3f,n1b,nd1b,n3b,nd3b
        integer,intent(in):: n1,n2,n3,nd1,nd2,nd3
        real*8,intent(inout)::z3(2,nd1*nd2*nd3b,2)
        real*8                  ::z1(2,nd1b*nd2*nd3,2) ! work array
        real*8,intent(out):: y(n1,n2,n3)

        isign=-1
        if (max(n1,n2,n3).gt.1024) stop '1024'

! some reasonable values of ncache: 
!   IBM/RS6000/590: 16*1024 ; IBM/RS6000/390: 3*1024 ; 
!   IBM/PwPC: 1*1024 ; SGI/MIPS/R8000: 16*1024 ; DEC/Alpha/EV5 and EV6 6*1024
!   But if you care about performance find the optimal value of ncache yourself!
!       On all vector machines: ncache=0

        ncache=6*1024
!        ncache=0

! check whether input values are reasonable
    if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
    if (isign.ne.1 .and. isign.ne.-1) stop 'wrong isign'
    if (n1.gt.nd1) stop 'n1>nd1'
    if (n2.gt.nd2) stop 'n2>nd2'
    if (n3.gt.nd3) stop 'n3>nd3'
    

!        call z3_to_z1(z3,z1,inzee)

      if (ncache.eq.0) then
! vector computer with memory banks:
        allocate(trig(2,1024),after(20),now(20),before(20))

        call ctrig_(n3,trig,after,before,now,isign,ic)
        nfft=nd1b*n2
        mm=nd1b*nd2
        do 51093,i=1,ic-1
        call fftstp_(mm,nfft,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
                          trig,after(i),now(i),before(i),isign)
51093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
                          trig,after(i),now(i),before(i),isign)

        inzee=3-inzee

        if (n2.ne.n3) call ctrig_(n2,trig,after,before,now,isign,ic)
        nfft=nd3*n1b
        mm=nd3*nd1b
        do 52093,i=1,ic-1
        call fftstp_(mm,nfft,nd2,mm,nd2,z1(1,1,inzee),z1(1,1,3-inzee), &
                           trig,after(i),now(i),before(i),isign)
52093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd2,mm,nd2,z1(1,1,inzee),z1(1,1,3-inzee), &
                       trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

        ! here we transform back from z1 to z3
        call z1_to_z3(z1,z3,inzee)

        if (n1.ne.n2) call ctrig_(n1,trig,after,before,now,isign,ic)
        nfft=nd2*n3b
        mm=nd2*nd3b
        do 53093,i=1,ic-1
        call fftstp_(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
                         trig,after(i),now(i),before(i),isign)
53093        inzee=3-inzee
        i=ic
        call fftrot(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
                         trig,after(i),now(i),before(i),isign)
        inzee=3-inzee

        call z3_to_y(z3,y,inzee)

! RISC machine with cache:
      else
! INtel IFC does not understand default(private)
!!!!!!!!!$omp parallel  default(private) &
!!!!!$omp parallel & 
!!!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
!!!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,isign,inzee,ncache) 
        npr=1
!!!!!$       npr=omp_get_num_threads()
        iam=0
!!!$       iam=omp_get_thread_num()
!        write(6,*) 'npr,iam',npr,iam
! Critical section only necessary on Intel
!!!!!!$omp critical
        allocate(zw(2,ncache/4,2),trig(2,1024),after(20),now(20),before(20))
!!!!!!!$omp end critical

        inzet=inzee
! TRANSFORM ALONG Z AXIS

        mm=nd1b*nd2
        m=nd3
        lot=max(1,ncache/(4*n3))
        nn=lot
        n=n3
        if (2*n*lot*2.gt.ncache) stop 'ncache1'

        call ctrig_(n3,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd1b*n2)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd1b*n2)
        nfft=mb-ma+1
        j=ma
        jj=j*nd3-nd3+1
        call fftrot(mm,nfft,m,mm,m,z1(1,j,inzet),z1(1,jj,3-inzet), &
                          trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd1b*n2)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd1b*n2)
        do 1000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd3-nd3+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z1(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 1093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
1093        inzeep=3-inzeep
        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z1(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
1000        continue
      endif

        inzet=3-inzet

!!!!!!!!$omp barrier

! TRANSFORM ALONG Y AXIS
        mm=nd3*nd1b
        m=nd2
        lot=max(1,ncache/(4*n2))
        nn=lot
        n=n2
        if (2*n*lot*2.gt.ncache) stop 'ncache2'

        if (n2.ne.n3) call ctrig_(n2,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd3*n1b)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd3*n1b)
        nfft=mb-ma+1
        j=ma
        jj=j*nd2-nd2+1
        call fftrot(mm,nfft,m,mm,m,z1(1,j,inzet),z1(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd3*n1b)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd3*n1b)
        do 2000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd2-nd2+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z1(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 2093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
2093        inzeep=3-inzeep

        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z1(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
2000        continue
      endif
        inzet=3-inzet

        call z1_to_z3(z1,z3,inzet)
!!!!!!$omp barrier

! TRANSFORM ALONG X AXIS
        mm=nd2*nd3b
        m=nd1
        lot=max(1,ncache/(4*n1))
        nn=lot
        n=n1
        if (2*n*lot*2.gt.ncache) stop 'ncache3'

        if (n1.ne.n2) call ctrig_(n1,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
        i=ic
        lotomp=(nd2*n3b)/npr+1
        ma=iam*lotomp+1
        mb=min((iam+1)*lotomp,nd2*n3b)
        nfft=mb-ma+1
        j=ma
        jj=j*nd1-nd1+1
        call fftrot(mm,nfft,m,mm,m,z3(1,j,inzet),z3(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)

      else

        lotomp=(nd2*n3b)/npr+1
        jompa=iam*lotomp+1
        jompb=min((iam+1)*lotomp,nd2*n3b)
        do 3000,j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft=mb-ma+1
        jj=j*nd1-nd1+1

        i=1
        inzeep=2
        call fftstp_(mm,nfft,m,nn,n,z3(1,j,inzet),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
        inzeep=1

        do 3093,i=2,ic-1
        call fftstp_(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                         trig,after(i),now(i),before(i),isign)
3093        inzeep=3-inzeep
        i=ic
        call fftrot(nn,nfft,n,mm,m,zw(1,1,inzeep),z3(1,jj,3-inzet), &
                         trig,after(i),now(i),before(i),isign)
3000        continue
      endif
        inzet=3-inzet
        
        call z3_to_y(z3,y,inzet)
        deallocate(zw,trig,after,now,before)
        if (iam.eq.0) inzee=inzet
!!!!!!!!!!!$omp end parallel  


      endif

contains

   subroutine z3_to_z1(z3,z1,inzee)
            ! transforms the data from the format z3:
            ! output of fft_for: stores only elements of z with i3=<nd3f
            ! 
            ! to the format z1:
            ! input of fft_back: stores only elements of z with i1=<nd1b
            implicit none
            integer,intent(in)::inzee
            real*8,intent(in):: z3(2,nd1 ,nd2,nd3f,2)
            real*8,intent(out)::z1(2,nd1b,nd2,nd3,2)
            integer i1,i2,i3

            ! i3=1: then z1 is contained in z3 
            do i2=1,n2
                do i1=1,n1b
                    z1(1,i1,i2,1,inzee)=z3(1,i1,i2,1,inzee)
                    z1(2,i1,i2,1,inzee)=z3(2,i1,i2,1,inzee)
                enddo
            enddo    

            do i3=2,n3f
                ! i2=1
                ! i1=1
                z1(1,1,1,i3,inzee)=z3(1,1,1,i3,inzee)
                z1(2,1,1,i3,inzee)=z3(2,1,1,i3,inzee)

                z1(1,1,1,n3+2-i3,inzee)=z3(1,1,1,i3,inzee)
                z1(2,1,1,n3+2-i3,inzee)=-z3(2,1,1,i3,inzee)

                ! i2=1
                do i1=2,n1b
                    z1(1,i1,1,i3,inzee)=z3(1,i1,1,i3,inzee)
                    z1(2,i1,1,i3,inzee)=z3(2,i1,1,i3,inzee)

                    z1(1,i1,1,n3+2-i3,inzee)= z3(1,n1+2-i1,1,i3,inzee)
                    z1(2,i1,1,n3+2-i3,inzee)=-z3(2,n1+2-i1,1,i3,inzee)
                enddo

                do i2=2,n2
                    ! i1=1
                    z1(1,1,i2,i3,inzee)=z3(1,1,i2,i3,inzee)
                    z1(2,1,i2,i3,inzee)=z3(2,1,i2,i3,inzee)

                    z1(1,1,i2,n3+2-i3,inzee)= z3(1,1,n2+2-i2,i3,inzee)
                    z1(2,1,i2,n3+2-i3,inzee)=-z3(2,1,n2+2-i2,i3,inzee)

                    do i1=2,n1b
                        z1(1,i1,i2,i3,inzee)=z3(1,i1,i2,i3,inzee)
                        z1(2,i1,i2,i3,inzee)=z3(2,i1,i2,i3,inzee)

                        z1(1,i1,i2,n3+2-i3,inzee)= z3(1,n1+2-i1,n2+2-i2,i3,inzee)
                        z1(2,i1,i2,n3+2-i3,inzee)=-z3(2,n1+2-i1,n2+2-i2,i3,inzee)
                    enddo
                enddo
            enddo
                                
    end subroutine z3_to_z1

    subroutine z1_to_z3(z1,z3,inzee)
    ! transforms the data from the format z1:
    ! stores the part of z with i1=<nd1b 
    ! to the format z3:
    ! stores the elements of z with even and odd values of i3
    ! as its even and odd parts w.r.t. flip of n1
    implicit none
    integer,intent(in)::inzee
    real*8,intent(in):: z1(2,nd2,nd3,nd1b,2)
    real*8,intent(out)::z3(2,nd2,nd3b,nd1,2)
    integer i1,i2,i3

    if (2*n3b.eq.n3) then 
        ! i1=1
        do i3=1,n3b
            do i2=1,n2
                z3(1,i2,i3,1,inzee)= z1(1,i2,2*i3-1,1,inzee)-z1(2,i2,2*i3,1,inzee)
                z3(2,i2,i3,1,inzee)= z1(2,i2,2*i3-1,1,inzee)+z1(1,i2,2*i3,1,inzee)
            enddo
        enddo

        do i1=2,n1b
            do i3=1,n3b
                do i2=1,n2
                    z3(1,i2,i3,i1     ,inzee)= z1(1,i2,2*i3-1,i1,inzee)-z1(2,i2,2*i3,i1,inzee)
                    z3(2,i2,i3,i1     ,inzee)= z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
                    z3(1,i2,i3,n1+2-i1,inzee)= z1(1,i2,2*i3-1,i1,inzee)+z1(2,i2,2*i3,i1,inzee)
                    z3(2,i2,i3,n1+2-i1,inzee)=-z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
                enddo
            enddo
        enddo
    else  ! 2*n3b=n3+1
        ! i1=1
        do i3=1,n3b-1
            do i2=1,n2
                z3(1,i2,i3,1,inzee)= z1(1,i2,2*i3-1,1,inzee)-z1(2,i2,2*i3,1,inzee)
                z3(2,i2,i3,1,inzee)= z1(2,i2,2*i3-1,1,inzee)+z1(1,i2,2*i3,1,inzee)
            enddo
        enddo

        do i1=2,n1b
            do i3=1,n3b-1
                do i2=1,n2
                    z3(1,i2,i3,i1     ,inzee)= z1(1,i2,2*i3-1,i1,inzee)-z1(2,i2,2*i3,i1,inzee)
                    z3(2,i2,i3,i1     ,inzee)= z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
                    z3(1,i2,i3,n1+2-i1,inzee)= z1(1,i2,2*i3-1,i1,inzee)+z1(2,i2,2*i3,i1,inzee)
                    z3(2,i2,i3,n1+2-i1,inzee)=-z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
                enddo
            enddo
        enddo

        ! i3=n3b is treated separately: 2*n3b-1=n3, but the terms with 2*n3b are
        ! omitted

        ! i1=1
        do i2=1,n2
            z3(1,i2,n3b,1,inzee)= z1(1,i2,n3,1,inzee)
            z3(2,i2,n3b,1,inzee)= z1(2,i2,n3,1,inzee)
        enddo

        do i1=2,n1b
            do i2=1,n2
                z3(1,i2,n3b,i1     ,inzee)= z1(1,i2,n3,i1,inzee)
                z3(2,i2,n3b,i1     ,inzee)= z1(2,i2,n3,i1,inzee)
                z3(1,i2,n3b,n1+2-i1,inzee)= z1(1,i2,n3,i1,inzee)
                z3(2,i2,n3b,n1+2-i1,inzee)=-z1(2,i2,n3,i1,inzee)
            enddo
        enddo

    endif

    end subroutine z1_to_z3

    subroutine z3_to_y(z3,y,inzee)
    ! transforms the output of FFT: z3, for which:
    ! real part of      z3 contains elements of y with odd  i3
    ! imaginary part of z3 contains elements of y with even i3

    ! into the final output: real array y.
    implicit none
    integer,intent(in)::inzee
    real*8,intent(in)::z3(2,nd1,nd2,nd3b,2)
    real*8,intent(out)::y(n1,n2,n3)
    integer i1,i2,i3
    real*8 fac

    fac=.5d0/(n1*n2*n3)

    if (2*n3b.eq.n3) then 
        do i3=1,n3b
            do i2=1,n2
                do i1=1,n1
                    y(i1,i2,2*i3-1)=z3(1,i1,i2,i3,inzee)*fac
                    y(i1,i2,2*i3  )=z3(2,i1,i2,i3,inzee)*fac
                enddo
            enddo
        enddo
    else ! 2*n3b=n3+1
        do i3=1,n3b-1
            do i2=1,n2
                do i1=1,n1
                    y(i1,i2,2*i3-1)=z3(1,i1,i2,i3,inzee)*fac
                    y(i1,i2,2*i3  )=z3(2,i1,i2,i3,inzee)*fac
                enddo
            enddo
        enddo

        ! i3=n3b is treated separately: 2*n3b-1=n3, but the terms with 2*n3b are
        ! omitted

        do i2=1,n2
            do i1=1,n1
                y(i1,i2,n3)=z3(1,i1,i2,n3b,inzee)*fac
            enddo
        enddo
    endif

    end subroutine z3_to_y


end subroutine fft_back
!!***


!!***f* BigDFT/ctrig_
!! SIDE EFFECTS
!!  Different factorizations affect the performance
!!  Factoring 64 as 4*4*4 might for example be faster on some machines than 8*8.
!! SOURCE
!!
subroutine ctrig_(n,trig,after,before,now,isign,ic)

  implicit real(kind=8) (a-h,o-z)
! Maximum number of points for FFT (should be same number in fft3d routine)
  integer, parameter :: nfft_max=24000
! Arguments
  integer :: n,isign,ic
  integer :: after(7),before(7),now(7)
  real(kind=8) :: trig(2,nfft_max)
! Local variables
  integer, parameter :: ndata = 180
  integer, dimension(7,ndata) :: idata 
! The factor 6 is only allowed in the first place!
        data ((idata(i,j),i=1,7),j=1,ndata) /                     &
            3,   3, 1, 1, 1, 1, 1,       4,   4, 1, 1, 1, 1, 1,   &
            5,   5, 1, 1, 1, 1, 1,       6,   6, 1, 1, 1, 1, 1,   &
            8,   8, 1, 1, 1, 1, 1,       9,   3, 3, 1, 1, 1, 1,   &
           12,   4, 3, 1, 1, 1, 1,      15,   5, 3, 1, 1, 1, 1,   &
           16,   4, 4, 1, 1, 1, 1,      18,   6, 3, 1, 1, 1, 1,   &
           20,   5, 4, 1, 1, 1, 1,      24,   8, 3, 1, 1, 1, 1,   &
           25,   5, 5, 1, 1, 1, 1,      27,   3, 3, 3, 1, 1, 1,   &
           30,   6, 5, 1, 1, 1, 1,      32,   8, 4, 1, 1, 1, 1,   &
           36,   4, 3, 3, 1, 1, 1,      40,   8, 5, 1, 1, 1, 1,   &
           45,   5, 3, 3, 1, 1, 1,      48,   4, 4, 3, 1, 1, 1,   &
           54,   6, 3, 3, 1, 1, 1,      60,   5, 4, 3, 1, 1, 1,   &
           64,   8, 8, 1, 1, 1, 1,      72,   8, 3, 3, 1, 1, 1,   &
           75,   5, 5, 3, 1, 1, 1,      80,   5, 4, 4, 1, 1, 1,   &
           81,   3, 3, 3, 3, 1, 1,      90,   6, 5, 3, 1, 1, 1,   &
           96,   8, 4, 3, 1, 1, 1,     100,   5, 5, 4, 1, 1, 1,   &
          108,   4, 3, 3, 3, 1, 1,     120,   8, 5, 3, 1, 1, 1,   &
          125,   5, 5, 5, 1, 1, 1,     128,   8, 4, 4, 1, 1, 1,   &
          135,   5, 3, 3, 3, 1, 1,     144,   6, 8, 3, 1, 1, 1,   &
          150,   6, 5, 5, 1, 1, 1,     160,   8, 5, 4, 1, 1, 1,   &
          162,   6, 3, 3, 3, 1, 1,     180,   5, 4, 3, 3, 1, 1,   &
          192,   6, 8, 4, 1, 1, 1,     200,   8, 5, 5, 1, 1, 1,   &
          216,   8, 3, 3, 3, 1, 1,     225,   5, 5, 3, 3, 1, 1,   &
          240,   6, 8, 5, 1, 1, 1,     243,   3, 3, 3, 3, 3, 1,   &
          256,   8, 8, 4, 1, 1, 1,     270,   6, 5, 3, 3, 1, 1,   &
          288,   8, 4, 3, 3, 1, 1,     300,   5, 5, 4, 3, 1, 1,   &
          320,   5, 4, 4, 4, 1, 1,     324,   4, 3, 3, 3, 3, 1,   &
          360,   8, 5, 3, 3, 1, 1,     375,   5, 5, 5, 3, 1, 1,   &
          384,   8, 4, 4, 3, 1, 1,     400,   5, 5, 4, 4, 1, 1,   &
          405,   5, 3, 3, 3, 3, 1,     432,   4, 4, 3, 3, 3, 1,   &
          450,   6, 5, 5, 3, 1, 1,     480,   8, 5, 4, 3, 1, 1,   &
          486,   6, 3, 3, 3, 3, 1,     500,   5, 5, 5, 4, 1, 1,   &
          512,   8, 8, 8, 1, 1, 1,     540,   5, 4, 3, 3, 3, 1,   &
          576,   4, 4, 4, 3, 3, 1,     600,   8, 5, 5, 3, 1, 1,   &
          625,   5, 5, 5, 5, 1, 1,     640,   8, 5, 4, 4, 1, 1,   &
          648,   8, 3, 3, 3, 3, 1,     675,   5, 5, 3, 3, 3, 1,   &
          720,   5, 4, 4, 3, 3, 1,     729,   3, 3, 3, 3, 3, 3,   &
          750,   6, 5, 5, 5, 1, 1,     768,   4, 4, 4, 4, 3, 1,   &
          800,   8, 5, 5, 4, 1, 1,     810,   6, 5, 3, 3, 3, 1,   &
          864,   8, 4, 3, 3, 3, 1,     900,   5, 5, 4, 3, 3, 1,   &
          960,   5, 4, 4, 4, 3, 1,     972,   4, 3, 3, 3, 3, 3,   &
         1000,   8, 5, 5, 5, 1, 1,    1024,   4, 4, 4, 4, 4, 1,   &
         1080,   6, 5, 4, 3, 3, 1,    1125,   5, 5, 5, 3, 3, 1,   &
         1152,   6, 4, 4, 4, 3, 1,    1200,   6, 8, 5, 5, 1, 1,   &
         1215,   5, 3, 3, 3, 3, 3,    1280,   8, 8, 5, 4, 1, 1,   &
         1296,   6, 8, 3, 3, 3, 1,    1350,   6, 5, 5, 3, 3, 1,   &
         1440,   6, 5, 4, 4, 3, 1,    1458,   6, 3, 3, 3, 3, 3,   &
         1500,   5, 5, 5, 4, 3, 1,    1536,   6, 8, 8, 4, 1, 1,   &
         1600,   8, 8, 5, 5, 1, 1,    1620,   5, 4, 3, 3, 3, 3,   &
         1728,   6, 8, 4, 3, 3, 1,    1800,   6, 5, 5, 4, 3, 1,   &
         1875,   5, 5, 5, 5, 3, 1,    1920,   6, 5, 4, 4, 4, 1,   &
         1944,   6, 4, 3, 3, 3, 3,    2000,   5, 5, 5, 4, 4, 1,   &
         2025,   5, 5, 3, 3, 3, 3,    2048,   8, 4, 4, 4, 4, 1,   &
         2160,   6, 8, 5, 3, 3, 1,    2250,   6, 5, 5, 5, 3, 1,   &
         2304,   6, 8, 4, 4, 3, 1,    2400,   6, 5, 5, 4, 4, 1,   &
         2430,   6, 5, 3, 3, 3, 3,    2500,   5, 5, 5, 5, 4, 1,   &
         2560,   8, 5, 4, 4, 4, 1,    2592,   6, 4, 4, 3, 3, 3,   &
         2700,   5, 5, 4, 3, 3, 3,    2880,   6, 8, 5, 4, 3, 1,   &
         3000,   6, 5, 5, 5, 4, 1,    3072,   6, 8, 4, 4, 4, 1,   &
         3125,   5, 5, 5, 5, 5, 1,    3200,   8, 5, 5, 4, 4, 1,   &
         3240,   6, 5, 4, 3, 3, 3,    3375,   5, 5, 5, 3, 3, 3,   &
         3456,   6, 4, 4, 4, 3, 3,    3600,   6, 8, 5, 5, 3, 1,   &
         3750,   6, 5, 5, 5, 5, 1,    3840,   6, 8, 5, 4, 4, 1,   &
         3888,   6, 8, 3, 3, 3, 3,    4000,   8, 5, 5, 5, 4, 1,   &
         4050,   6, 5, 5, 3, 3, 3,    4096,   8, 8, 4, 4, 4, 1,   &
         4320,   6, 5, 4, 4, 3, 3,    4500,   5, 5, 5, 4, 3, 3,   &
         4608,   6, 8, 8, 4, 3, 1,    4800,   6, 8, 5, 5, 4, 1,   &
         5000,   8, 5, 5, 5, 5, 1,    5120,   8, 8, 5, 4, 4, 1,   &
         5184,   6, 8, 4, 3, 3, 3,    5400,   6, 5, 5, 4, 3, 3,   &
         5625,   5, 5, 5, 5, 3, 3,    5760,   6, 8, 8, 5, 3, 1,   &
         6000,   6, 8, 5, 5, 5, 1,    6144,   6, 8, 8, 4, 4, 1,   &
         6400,   8, 8, 5, 5, 4, 1,    6480,   6, 8, 5, 3, 3, 3,   &
         6750,   6, 5, 5, 5, 3, 3,    6912,   6, 8, 4, 4, 3, 3,   &
         7200,   6, 5, 5, 4, 4, 3,    7500,   5, 5, 5, 5, 4, 3,   &
         7680,   6, 8, 8, 5, 4, 1,    8000,   8, 8, 5, 5, 5, 1,   &
         8192,   8, 8, 8, 4, 4, 1,    8640,   8, 8, 5, 3, 3, 3,   &
         9000,   8, 5, 5, 5, 3, 3,    9216,   6, 8, 8, 8, 3, 1,   &
         9375,   5, 5, 5, 5, 5, 3,    9600,   8, 5, 5, 4, 4, 3,   &
        10000,   5, 5, 5, 5, 4, 4,   10240,   8, 8, 8, 5, 4, 1,   &
        10368,   6, 8, 8, 3, 3, 3,   10800,   6, 8, 5, 5, 3, 3,   &
        11250,   6, 5, 5, 5, 5, 3,   11520,   8, 8, 5, 4, 3, 3,   &
        12000,   8, 5, 5, 5, 4, 3,   12288,   8, 8, 8, 8, 3, 1,   &
        12500,   5, 5, 5, 5, 5, 4,   12800,   8, 8, 8, 5, 5, 1,   &
        13824,   8, 8, 8, 3, 3, 3,   14400,   8, 8, 5, 5, 3, 3,   &
        15000,   8, 5, 5, 5, 5, 3,   15360,   6, 8, 8, 8, 5, 1,   &
        15625,   5, 5, 5, 5, 5, 5,   16000,   8, 5, 5, 5, 4, 4,   &
        16384,   8, 8, 8, 8, 4, 1,   17280,   6, 8, 8, 5, 3, 3,   &
        18000,   6, 8, 5, 5, 5, 3,   18432,   8, 8, 8, 4, 3, 3,   &
        18750,   6, 5, 5, 5, 5, 5,   19200,   8, 8, 5, 5, 4, 3,   &
        20000,   8, 5, 5, 5, 5, 4,   20480,   8, 8, 8, 8, 5, 1,   &
        23040,   8, 8, 8, 5, 3, 3,   24000,   8, 8, 5, 5, 5, 3    /
        
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

end subroutine ctrig_
!!***


!ccccccccccccccccccccccccccccccccccccccccccccccc

         subroutine fftstp_(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

        implicit real*8 (a-h,o-z)
        integer after,before,atn,atb
        dimension trig(2,1024),zin(2,mm,m),zout(2,nn,n)
        atn=after*now
        atb=after*before

!         sqrt(.5d0)
        rt2i=0.7071067811865475d0
        if (now.eq.2) then
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 2001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        do 2001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        zout(1,j,nout1)= r2 + r1
        zout(2,j,nout1)= s2 + s1
        zout(1,j,nout2)= r1 - r2
        zout(2,j,nout2)= s1 - s2
2001        continue
        do 2000,ia=2,after
        ias=ia-1
        if (2*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,j,nout1)= r1 - r2
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r2 + r1
                        zout(2,j,nout2)= s1 - s2
2010                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s1 - s2
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s2 + s1
2020                        continue
                endif
        else if (4*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2030,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2030,j=1,nfft
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
2030                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2040,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2040,j=1,nfft
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
2040                        continue
                endif
        else if (4*ias.eq.3*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2050,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2050,j=1,nfft
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
2050                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2060,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2060,j=1,nfft
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
2060                        continue
                endif
        else
                itrig=ias*before+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                do 2090,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nout1=nout1+atn
                nout2=nout1+after
                do 2090,j=1,nfft
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
2090                continue
        endif
2000        continue
        else if (now.eq.4) then
        if (isign.eq.1) then 
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do 4001,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do 4001,j=1,nfft
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
4001                continue
                do 4000,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 4010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4010,j=1,nfft
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
4010                        continue
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
                        do 4020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4020,j=1,nfft
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
4020                        continue
                endif
4000                continue
        else
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do 4101,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do 4101,j=1,nfft
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
4101                continue
                do 4100,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 4110,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4110,j=1,nfft
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
4110                        continue
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
                        do 4120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4120,j=1,nfft
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
4120                        continue
                endif
4100                continue
        endif
        else if (now.eq.8) then
        if (isign.eq.-1) then 
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do 8120,ib=1,before
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
                        do 8120,j=1,nfft
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
8120                        continue
                do 8000,ia=2,after
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
                        do 8020,ib=1,before
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
                        do 8020,j=1,nfft
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

8020                        continue
8000                continue

        else
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do 8121,ib=1,before
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
                        do 8121,j=1,nfft
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
8121                        continue

                do 8001,ia=2,after
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
                        do 8021,ib=1,before
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
                        do 8021,j=1,nfft
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
8021                        continue
8001                continue

        endif
        else if (now.eq.3) then 
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 3001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do 3001,j=1,nfft
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
3001        continue
        do 3000,ia=2,after
        ias=ia-1
        if (4*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do 3010,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3010,j=1,nfft
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
3010                continue
        else
                nin1=ia-after
                nout1=ia-atn
                do 3020,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3020,j=1,nfft
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
3020                continue
        endif
        else if (8*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do 3030,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3030,j=1,nfft
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
3030                continue
        else
                nin1=ia-after
                nout1=ia-atn
                do 3040,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3040,j=1,nfft
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
3040                continue
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
        do 3090,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do 3090,j=1,nfft
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
3090        continue
        endif
3000        continue
        else if (now.eq.5) then
!         cos(2.d0*pi/5.d0)
        cos2=0.3090169943749474d0
!         cos(4.d0*pi/5.d0)
        cos4=-0.8090169943749474d0
!        sin(2.d0*pi/5.d0)
        sin2=isign*0.9510565162951536d0
!         sin(4.d0*pi/5.d0)
        sin4=isign*0.5877852522924731d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 5001,ib=1,before
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
        do 5001,j=1,nfft
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
5001        continue
        do 5000,ia=2,after
        ias=ia-1
        if (8*ias.eq.5*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 5010,ib=1,before
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
                        do 5010,j=1,nfft
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
5010                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 5020,ib=1,before
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
                        do 5020,j=1,nfft
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
5020                        continue
                endif
        else
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
                do 5100,ib=1,before
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
                do 5100,j=1,nfft
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
5100                continue
        endif
5000        continue
       else if (now.eq.6) then
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0

        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 6120,ib=1,before
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
        do 6120,j=1,nfft
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
6120        continue

        else 
        stop 'error fftstp_'
        endif

        return
        end



         subroutine fftrot(mm,nfft,m,nn,n,zin,zout,trig,after,now,before,isign)
!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1995, 1999
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

        implicit real*8 (a-h,o-z)
        integer after,before,atn,atb
        dimension trig(2,1024),zin(2,mm,m),zout(2,n,nn)
        atn=after*now
        atb=after*before

!         sqrt(.5d0)
        rt2i=0.7071067811865475d0
        if (now.eq.2) then
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 2001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nout1=nout1+atn
        nout2=nout1+after
        do 2001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        zout(1,nout1,j)= r2 + r1
        zout(2,nout1,j)= s2 + s1
        zout(1,nout2,j)= r1 - r2
        zout(2,nout2,j)= s1 - s2
2001        continue
        do 2000,ia=2,after
        ias=ia-1
        if (2*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2010,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,nout1,j)= r1 - r2
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r2 + r1
                        zout(2,nout2,j)= s1 - s2
2010                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2020,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(2,j,nin2)
                        s2=zin(1,j,nin2)
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s1 - s2
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s2 + s1
2020                        continue
                endif
        else if (4*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2030,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2030,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s1 - s2
2030                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2040,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2040,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s1 - s2
2040                        continue
                endif
        else if (4*ias.eq.3*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 2050,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2050,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(r - s)*rt2i
                        zout(1,nout1,j)= r1 - r2
                        zout(2,nout1,j)= s2 + s1
                        zout(1,nout2,j)= r2 + r1
                        zout(2,nout2,j)= s1 - s2
2050                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 2060,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do 2060,j=1,nfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(s - r)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,nout1,j)= r2 + r1
                        zout(2,nout1,j)= s1 - s2
                        zout(1,nout2,j)= r1 - r2
                        zout(2,nout2,j)= s2 + s1
2060                        continue
                endif
        else
                itrig=ias*before+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                do 2090,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nout1=nout1+atn
                nout2=nout1+after
                do 2090,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                zout(1,nout1,j)= r2 + r1
                zout(2,nout1,j)= s2 + s1
                zout(1,nout2,j)= r1 - r2
                zout(2,nout2,j)= s1 - s2
2090                continue
        endif
2000        continue
        else if (now.eq.4) then
        if (isign.eq.1) then 
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do 4001,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do 4001,j=1,nfft
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
                zout(1,nout1,j) = r + s
                zout(1,nout3,j) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,nout2,j) = r - s 
                zout(1,nout4,j) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,nout1,j) = r + s 
                zout(2,nout3,j) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,nout2,j) = r + s 
                zout(2,nout4,j) = r - s
4001                continue
                do 4000,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 4010,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4010,j=1,nfft
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
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 + r3
                        s=s2 - s4
                        zout(1,nout2,j) = r - s 
                        zout(1,nout4,j) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s 
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 + r4
                        zout(2,nout2,j) = r + s 
                        zout(2,nout4,j) = r - s
4010                        continue
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
                        do 4020,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4020,j=1,nfft
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
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r - s 
                        zout(1,nout4,j) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s 
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r + s 
                        zout(2,nout4,j) = r - s
4020                        continue
                endif
4000                continue
        else
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do 4101,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do 4101,j=1,nfft
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
                zout(1,nout1,j) = r + s
                zout(1,nout3,j) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,nout2,j) = r + s
                zout(1,nout4,j) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,nout1,j) = r + s
                zout(2,nout3,j) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,nout2,j) = r - s
                zout(2,nout4,j) = r + s
4101                continue
                do 4100,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 4110,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4110,j=1,nfft
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
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 + s4
                        zout(1,nout2,j) = r + s
                        zout(1,nout4,j) = r - s
                        r=s1 - s3
                        s=s2 - s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 + s3
                        s=r2 - r4
                        zout(2,nout2,j) = r - s
                        zout(2,nout4,j) = r + s
4110                        continue
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
                        do 4120,ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do 4120,j=1,nfft
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
                        zout(1,nout1,j) = r + s
                        zout(1,nout3,j) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,nout2,j) = r + s
                        zout(1,nout4,j) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,nout1,j) = r + s
                        zout(2,nout3,j) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,nout2,j) = r - s
                        zout(2,nout4,j) = r + s
4120                        continue
                endif
4100                continue
        endif
        else if (now.eq.8) then
        if (isign.eq.-1) then 
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do 8120,ib=1,before
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
                        do 8120,j=1,nfft
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
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp1
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp1
                        zout(1,nout3,j) = am + dm
                        zout(2,nout3,j) = cm - bm
                        zout(1,nout7,j) = am - dm
                        zout(2,nout7,j) = cm + bm
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
                        dp1= ( cm - dp1)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp1
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp1
8120                        continue
                do 8000,ia=2,after
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
                        do 8020,ib=1,before
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
                        do 8020,j=1,nfft
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
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp1
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp1
                        zout(1,nout3,j) = am + dm
                        zout(2,nout3,j) = cm - bm
                        zout(1,nout7,j) = am - dm
                        zout(2,nout7,j) = cm + bm
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
                        dp1= ( cm - dp1)*rt2i
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp1
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp1

8020                        continue
8000                continue

        else
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do 8121,ib=1,before
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
                        do 8121,j=1,nfft
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
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp1
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp1
                        zout(1,nout3,j) = am - dm
                        zout(2,nout3,j) = cm + bm
                        zout(1,nout7,j) = am + dm
                        zout(2,nout7,j) = cm - bm
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
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp1
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp1
8121                        continue

                do 8001,ia=2,after
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
                        do 8021,ib=1,before
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
                        do 8021,j=1,nfft
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
                        zout(1,nout1,j) = ap + bp
                        zout(2,nout1,j) = cp + dp1
                        zout(1,nout5,j) = ap - bp
                        zout(2,nout5,j) = cp - dp1
                        zout(1,nout3,j) = am - dm
                        zout(2,nout3,j) = cm + bm
                        zout(1,nout7,j) = am + dm
                        zout(2,nout7,j) = cm - bm
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
                        zout(1,nout2,j) = ap + r
                        zout(2,nout2,j) = bm + s
                        zout(1,nout6,j) = ap - r
                        zout(2,nout6,j) = bm - s
                        zout(1,nout4,j) = am + cp
                        zout(2,nout4,j) = bp + dp1
                        zout(1,nout8,j) = am - cp
                        zout(2,nout8,j) = bp - dp1
8021                        continue
8001                continue

        endif
        else if (now.eq.3) then 
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 3001,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do 3001,j=1,nfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r=r2 + r3
        s=s2 + s3
        zout(1,nout1,j) = r + r1
        zout(2,nout1,j) = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,nout2,j) = r1 - s2 
        zout(2,nout2,j) = s1 + r2
        zout(1,nout3,j) = r1 + s2 
        zout(2,nout3,j) = s1 - r2
3001        continue
        do 3000,ia=2,after
        ias=ia-1
        if (4*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do 3010,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3010,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 + r3
                s=s2 - s3
                zout(1,nout1,j) = r1 - r
                zout(2,nout1,j) = s + s1
                r1=r1 + .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2-r3)        
                s2=bb*(s2+s3)
                zout(1,nout2,j) = r1 - s2 
                zout(2,nout2,j) = s1 - r2
                zout(1,nout3,j) = r1 + s2 
                zout(2,nout3,j) = s1 + r2
3010                continue
        else
                nin1=ia-after
                nout1=ia-atn
                do 3020,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3020,j=1,nfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s1 - s
                r1=r1 - .5d0*r
                s1=s1 + .5d0*s        
                r2=bb*(r2+r3)        
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 + s2 
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 - s2 
                zout(2,nout3,j) = s1 - r2
3020                continue
        endif
        else if (8*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do 3030,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3030,j=1,nfft
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
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2+r3)        
                s2=bb*(s2-s3)
                zout(1,nout2,j) = r1 - s2 
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2 
                zout(2,nout3,j) = s1 - r2
3030                continue
        else
                nin1=ia-after
                nout1=ia-atn
                do 3040,ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do 3040,j=1,nfft
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
                zout(1,nout1,j) = r + r1
                zout(2,nout1,j) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s        
                r2=bb*(r2-r3)        
                s2=bb*(s2+s3)
                zout(1,nout2,j) = r1 - s2 
                zout(2,nout2,j) = s1 + r2
                zout(1,nout3,j) = r1 + s2 
                zout(2,nout3,j) = s1 - r2
3040                continue
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
        do 3090,ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do 3090,j=1,nfft
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
        zout(1,nout1,j) = r + r1
        zout(2,nout1,j) = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,nout2,j) = r1 - s2 
        zout(2,nout2,j) = s1 + r2
        zout(1,nout3,j) = r1 + s2 
        zout(2,nout3,j) = s1 - r2
3090        continue
        endif
3000        continue
        else if (now.eq.5) then
!         cos(2.d0*pi/5.d0)
        cos2=0.3090169943749474d0
!         cos(4.d0*pi/5.d0)
        cos4=-0.8090169943749474d0
!        sin(2.d0*pi/5.d0)
        sin2=isign*0.9510565162951536d0
!         sin(4.d0*pi/5.d0)
        sin4=isign*0.5877852522924731d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 5001,ib=1,before
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
        do 5001,j=1,nfft
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
        zout(1,nout1,j) = r1 + r25 + r34
        r = r1 + cos2*r25 + cos4*r34
        s = sin2*s25 + sin4*s34
        zout(1,nout2,j) = r - s
        zout(1,nout5,j) = r + s
        r = r1 + cos4*r25 + cos2*r34
        s = sin4*s25 - sin2*s34
        zout(1,nout3,j) = r - s
        zout(1,nout4,j) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,nout1,j) = s1 + s25 + s34
        r = s1 + cos2*s25 + cos4*s34
        s = sin2*r25 + sin4*r34
        zout(2,nout2,j) = r + s
        zout(2,nout5,j) = r - s
        r = s1 + cos4*s25 + cos2*s34
        s = sin4*r25 - sin2*r34
        zout(2,nout3,j) = r + s
        zout(2,nout4,j) = r - s
5001        continue
        do 5000,ia=2,after
        ias=ia-1
        if (8*ias.eq.5*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do 5010,ib=1,before
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
                        do 5010,j=1,nfft
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
                        zout(1,nout1,j) = r1 + r25 - r34
                        r = r1 + cos2*r25 - cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,nout2,j) = r - s
                        zout(1,nout5,j) = r + s
                        r = r1 + cos4*r25 - cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,nout3,j) = r - s
                        zout(1,nout4,j) = r + s
                        r25 = r2 + r5
                        r34 = r4 - r3
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,nout1,j) = s1 + s25 + s34
                        r = s1 + cos2*s25 + cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,nout2,j) = r + s
                        zout(2,nout5,j) = r - s
                        r = s1 + cos4*s25 + cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,nout3,j) = r + s
                        zout(2,nout4,j) = r - s
5010                        continue
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do 5020,ib=1,before
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
                        do 5020,j=1,nfft
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
                        zout(1,nout1,j) = r1 + r25 + r34
                        r = r1 + cos2*r25 + cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,nout2,j) = r - s
                        zout(1,nout5,j) = r + s
                        r = r1 + cos4*r25 + cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,nout3,j) = r - s
                        zout(1,nout4,j) = r + s
                        r25 = r2 + r5
                        r34 = r3 - r4
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,nout1,j) = s1 + s25 - s34
                        r = s1 + cos2*s25 - cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,nout2,j) = r + s
                        zout(2,nout5,j) = r - s
                        r = s1 + cos4*s25 - cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,nout3,j) = r + s
                        zout(2,nout4,j) = r - s
5020                        continue
                endif
        else
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
                do 5100,ib=1,before
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
                do 5100,j=1,nfft
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
                zout(1,nout1,j) = r1 + r25 + r34
                r = r1 + cos2*r25 + cos4*r34
                s = sin2*s25 + sin4*s34
                zout(1,nout2,j) = r - s
                zout(1,nout5,j) = r + s
                r = r1 + cos4*r25 + cos2*r34
                s = sin4*s25 - sin2*s34
                zout(1,nout3,j) = r - s
                zout(1,nout4,j) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,nout1,j) = s1 + s25 + s34
                r = s1 + cos2*s25 + cos4*s34
                s = sin2*r25 + sin4*r34
                zout(2,nout2,j) = r + s
                zout(2,nout5,j) = r - s
                r = s1 + cos4*s25 + cos2*s34
                s = sin4*r25 - sin2*r34
                zout(2,nout3,j) = r + s
                zout(2,nout4,j) = r - s
5100                continue
        endif
5000        continue
       else if (now.eq.6) then
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0

        ia=1
        nin1=ia-after
        nout1=ia-atn
        do 6120,ib=1,before
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
        do 6120,j=1,nfft
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

        zout(1,nout1,j)=ur1+vr1
        zout(2,nout1,j)=ui1+vi1
        zout(1,nout5,j)=ur2+vr2
        zout(2,nout5,j)=ui2+vi2
        zout(1,nout3,j)=ur3+vr3
        zout(2,nout3,j)=ui3+vi3
        zout(1,nout4,j)=ur1-vr1
        zout(2,nout4,j)=ui1-vi1
        zout(1,nout2,j)=ur2-vr2
        zout(2,nout2,j)=ui2-vi2
        zout(1,nout6,j)=ur3-vr3
        zout(2,nout6,j)=ui3-vi3
6120        continue

       else
        stop 'error fftrot'
       endif

        return
        end

