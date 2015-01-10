!{\src2tex{textfont=tt}}
!!****f* ABINIT/sqnormm_v
!! NAME
!! sqnormm_v
!!
!! FUNCTION
!! For a series of potentials,
!! compute square of the norm (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden), and sum over them.
!! Need the index of the first potential to be treated, in the provided array
!! of potentials, and the number of potentials to be treated.
!! Might be used to compute just one square of norm, in a big array, such as to avoid 
!! copying a potential from a big array to a temporary place.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2011 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space function on FFT grid is REAL, if 2, COMPLEX
!!  index=index of the first potential to be treated
!!  mpi_comm=the mpi communicator used for the summation
!!  mpi_sumarize=set it to .true. if parallelisation is done over FFT
!!  mult=number of potentials to be treated
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot= third dimension of the potarr array
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potential is stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potential is stored as V, B_x, B_y, Bz  (B=magn. field)
!!  potarr(cplex*nfft,nspden,npot)=array of real space potentials on FFT grid
!!
!! OUTPUT
!!  norm2(mult)= value of the square of the norm of the different potentials
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sqnormm_v(cplex,index,mpi_comm, mpi_summarize,mult,nfft,norm2,npot,nspden,opt_storage,potarr)

 use defs_basis
 use abi_m_xmpi

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,index,mult,nfft,npot,nspden,opt_storage,mpi_comm
 logical, intent(in) :: mpi_summarize
!arrays
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: norm2(mult)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ii,ispden
 real(dp) :: ar
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

 do ii=1,mult
   ar=zero
   do ispden=1,min(nspden,2)
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
       ar=ar + potarr(ifft,ispden,index+ii-1)**2
     end do
!    $OMP END PARALLEL DO
   end do
   norm2(ii)=ar
   if (nspden==4) then
     ar=zero
     do ispden=3,4
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
       do ifft=1,cplex*nfft
         ar=ar + potarr(ifft,ispden,index+ii-1)**2
       end do
!      $OMP END PARALLEL DO
     end do
     if (opt_storage==0) then
       if (cplex==1) then
         norm2(ii)=norm2(ii)+two*ar
       else
         norm2(ii)=norm2(ii)+ar
       end if
     else
       norm2(ii)=half*(norm2(ii)+ar)
     end if
   end if
 end do

 if (mpi_summarize) then
   call abi_xmpi_sum(norm2,mpi_comm ,ierr)
 end if

end subroutine sqnormm_v
!!***
