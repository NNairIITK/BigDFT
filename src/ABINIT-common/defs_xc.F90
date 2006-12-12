!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_xc
!! NAME
!! defs_xc
!!
!! FUNCTION
!! This module contains definitions for the XC routines in the Src_3xc
!! directory and for the routines who call the XC routines, especially when
!! optional arguments are defined.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2005 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_xc

 interface
  subroutine rhohxc(enxc,gsqcut,intxc,ixc,izero,kxc,mpi_enreg,nfft,ngfft,&
&   nkxc,nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,vhartr,&
&   vxc,vxcavg,xccc3d,k3xc)
   use defs_basis
   use defs_datatypes

   integer,  intent(in)  :: intxc,ixc,izero,nfft,nkxc,nspden,n3xccc,option
   integer,  intent(in)  :: ngfft(18)
   real(dp), intent(in)  :: gsqcut
   real(dp), intent(out) :: enxc,vxcavg
   real(dp), intent(in)  :: rhog(2,nfft),rhor(nfft,nspden),&
&   rprimd(3,3),xccc3d(n3xccc)
   real(dp), intent(out) :: kxc(nfft,nkxc),strsxc(6),&
&   vhartr(nfft),vxc(nfft,nspden)
   real(dp), intent(out), optional:: k3xc(1:nfft)
   type(MPI_type),intent(inout) :: mpi_enreg
  end subroutine rhohxc
 end interface

 interface
    subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,nvxcdgr,     & !Mandatory arguments
         &                  dvxc,d2vxc,grho2_updn,vxcgr)    !Optional arguments 
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: ixc,npts,nspden,order
      integer,intent(in) :: ndvxc,nvxcdgr
      !arrays
      real(kind=8),intent(in) :: rho_updn(npts,nspden)
      real(kind=8),intent(in), optional :: grho2_updn(npts,2*nspden-1)
      real(kind=8),intent(out) :: exc(npts),vxc(npts,nspden)
      real(kind=8),intent(out), optional :: d2vxc(npts),dvxc(npts,ndvxc),vxcgr(npts,nvxcdgr)
    end subroutine drivexc
!!$  subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,vxcgr,   & !Mandatory arguments
!!$&                   dvxc,d2vxc,grho2_updn)                             !Optional arguments
!!$   use defs_basis
!!$   integer,  intent(in) :: ixc,npts,nspden,order
!!$   real(dp), intent(in) :: rho_updn(npts,nspden)
!!$   real(dp), intent(in), optional :: grho2_updn(npts,2*nspden-1)
!!$   real(dp), intent(out):: exc(npts),vxc(npts,nspden),vxcgr(npts,3)
!!$   real(dp), intent(out), optional :: dvxc(:,:),d2vxc(npts)
!!$  end subroutine drivexc
 end interface
 interface
    subroutine xctetr(exc,npt,order,rhor,rspts,vxc,& !Mandatory arguments
         &                 d2vxc,dvxc)                    !Optional arguments
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npt,order
      !arrays
      real(kind=8),intent(in) :: rhor(npt),rspts(npt)
      real(kind=8),intent(out) :: exc(npt),vxc(npt)
      real(kind=8),intent(out), optional :: d2vxc(npt),dvxc(npt)
    end subroutine xctetr
    subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,& !Mandatory Arguments
         &                d2vxci,dvxcdgr,dvxci,grho2_updn)       !Optional Arguments
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npts,nspden,option,order
      integer,intent(in) :: ndvxci
      !arrays
      real(kind=8),intent(in) :: rho_updn(npts,nspden)
      real(kind=8),intent(in), optional :: grho2_updn(npts,2*nspden-1)
      real(kind=8),intent(out) :: exci(npts)
      real(kind=8),intent(out), optional :: d2vxci(npts),dvxcdgr(npts,3),dvxci(npts,ndvxci)
      real(kind=8),intent(out) :: vxci(npts,nspden)
    end subroutine xcpbe
    subroutine xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,& !Mandatory arguments
         &                 dvxc)                            !Optional arguments
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npts,nspden,order
      integer,intent(in) :: ndvxc
      !arrays
      real(kind=8),intent(in) :: rspts(npts),zeta(npts)
      real(kind=8),intent(out) :: exc(npts),vxc(npts,nspden)
      real(kind=8),intent(out), optional :: dvxc(npts,ndvxc)
    end subroutine xcspol
    subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
         &                dvxc)                            !Optional arguments
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npt,order
      !arrays
      real(kind=8),intent(in) :: rhor(npt),rspts(npt)
      real(kind=8),intent(out) :: exc(npt),vxc(npt)
      real(kind=8),intent(out), optional :: dvxc(npt)
    end subroutine xcpzca
    subroutine xcwign(exc,npt,order,rhor,rspts,vxc,& !Mandatory arguments
         &                dvxc)                           !Optional arguments
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npt,order
      !arrays
      real(kind=8),intent(in) :: rhor(npt),rspts(npt)
      real(kind=8),intent(out) :: exc(npt),vxc(npt)
      real(kind=8),intent(out), optional :: dvxc(npt)
    end subroutine xcwign
    subroutine xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npt,order
      !arrays
      real(kind=8),intent(in) :: rspts(npt)
      real(kind=8),intent(out) :: exc(npt),vxc(npt)
      real(kind=8),intent(out), optional :: dvxc(npt)
    end subroutine xchelu
    subroutine xcxalp(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: npt,order
      !arrays
      real(kind=8),intent(in) :: rspts(npt)
      real(kind=8),intent(out) :: exc(npt),vxc(npt)
      real(kind=8),intent(out),optional :: dvxc(npt)
    end subroutine xcxalp
 end interface

end module defs_xc
!!***
