!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_pawinit
!! NAME
!! abi_pawinit
!!
!! FUNCTION
!! Initialize some starting values of several arrays used in
!! PAW calculations
!!
!! 1-Initialize data related to angular mesh
!! 2-Tabulate normalized shape function g(r)
!! 3-Compute indklmn indexes giving some l,m,n,lp,mp,np info
!!                           from klmn=[(l,m,n),(lp,mp,np)]
!! 4-Compute various factors/sizes (depending on (l,m,n))
!! 5-Compute $q_ijL=\displaystyle
!!                  \int_{0}^{r_c}{(\phi_i\phi_j-\widetilde{\phi_i}\widetilde{\phi_j}) r^l\,dr}
!!                   Gaunt(l_i m_i,l_j m_j,l m))$
!!           $S_ij=\displaystyle \sqrt{4 \pi} q_ij0$
!! 6-Compute $e_ijkl= vh1_ijkl - Vhatijkl - Bijkl - Cijkl$
!!     With:
!!       $vh1_ijkl =\sum_{L,m} {vh1*Gaunt(i,j,Lm)*Gaunt(k,l,Lm)}$
!!       $Vhat_ijkl=\sum_{L,m} {vhatijL*Gaunt(i,j,Lm)*q_klL}$
!!       $B_ijkl   =\sum_{L,m} {vhatijL*Gaunt(k,l,Lm)*q_ijL}$
!!       $C_ijkl   =\sum_{L,m} {intvhatL*q_ijL*q_klL}$
!!     and:
!!       vh1 according to eq. (A17) in Holzwarth et al., PRB 55, 2005 (1997)
!! 7-Compute Ex-correlation energy for the core density
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gnt_option=flag activated if pawang%gntselect and pawang%realgnt have to be allocated
!!             also determine the size of these pointers
!!  gsqcut_shp=effective cut-off to determine shape functions in reciprocal space
!!  lcutdens=max. l for densities/potentials moments computations
!!  lmix=max. l for which spherical terms will be mixed durinf SCF cycle
!!  mpsang=1+maximum angular momentum
!!  nphi="phi" dimension of paw angular mesh
!!  nsym=Number of symmetry elements in space group
!!  ntheta="theta" dimension of paw angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)=
!!                     selection rules for Gaunt coefficients
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!     %radfact(mesh_size)=Factor used to compute radial integrals
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!     %basis_size=Number of elements for the PAW nl basis
!!     %l_size=Maximum value of l+1 leading to non zero Gaunt coeffs
!!     %lmn_size=Number of (l,m,n) elements for the PAW basis
!!     %lmn2_size=lmn_size*(lmn_size+1)/2
!!     %dltij(lmn2_size)=factors used to compute sums over (ilmn,jlmn)
!!     %phi(mesh_size,basis_size)=PAW all electron wavefunctions
!!     %rshp=shape function radius (radius for compensation charge)
!!     %shape_type=Radial shape function type
!!     %shape_alpha=Alpha parameters in Bessel shape function
!!     %shape_lambda=Lambda parameter in gaussian shape function
!!     %shape_q=Q parameters in Bessel shape function
!!     %shape_sigma=Sigma parameter in gaussian shape function
!!     %tphi(mesh_size,basis_size)=PAW atomic pseudowavefunctions
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1=dev. on moments)
!!  xclevel=XC functional level (1=LDA, 2=GGA)
!!
!! OUTPUT
!!  pawang
!!     %gntselect(l_size_max**2,l_max**2*(l_max**2+1)/2)=selection rules for Gaunt coefficients
!!     %l_max=maximum value of angular momentum l+1
!!     %l_size_max=maximum value of angular momentum l_size=2*l_max-1
!!     %nsym=number of symmetry elements in space group
!!     %ngnt=number of non-zero Gaunt coefficients
!!     %realgnt(pawang%ngnt)=non-zero real Gaunt coefficients
!!     === only if pawxcdev==1 ==
!!       %anginit(3,angl_size)=for each point of the angular mesh, gives the coordinates
!!                             of the corresponding point on an unitary sphere
!!       %angl_size=dimension of paw angular mesh (angl_size=ntheta*nphi)
!!       %angwgth(angl_size)=for each point of the angular mesh, gives the weight
!!                           of the corresponding point on an unitary sphere
!!       %ntheta, nphi=dimensions of paw angular mesh
!!       %ylmr(l_size_max**2,angl_size)=real Ylm calculated in real space
!!       %ylmrgr(1:3,l_size_max**2,angl_size)=first gradients of real Ylm calculated in real space
!!     === only if pawspnorb==1 ==
!!     %ls_ylm(2,l_max**2,l_max**2,4)=LS operator in the real spherical harmonics basis
!!     %use_ls_ylm=flag activated if ls_ylm is allocated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated data read at start:
!!     %lcut_size_=max. value of l+1 leading to non zero Gaunt coeffs modified by lcutdens
!!     %lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
!!     %mqgrid_shp=number of points in reciprocal space for shape function
!!     %indklmn(8,lmn2_size)=array giving klm, kln, abs(il-jl) and (il+jl), ilmn and jlmn for each klmn=(ilmn,jlmn)
!!     %dshpfunc(mesh_size,l_size,4)=derivatives of shape function (used only for numerical shape functions)
!!     %eijkl(lmn2_size,lmn2_size)=part of the Dij that depends only from the projected occupation coeffs
!!     %exccore=Exchange-correlation energy for the core density
!!     %gnorm(l_size)=normalization factor of radial shape function
!!     %phiphj(:,:)=useful product Phi(:,i)*Phi(:,j)
!!     %qgrid_shp(mqgrid_shp)=points in reciprocal space for shape function
!!     %qijl(l_size**2,lmn2_size)=moments of the difference charge density between AE and PS partial wave
!!     %rad_for_spline(mesh_size)=radial grid used for spline (copy of pawrad%rad)
!!     %shapefunc(mesh_size,l_size)=normalized radial shape function
!!     %shapefncg(mqgrid_shp,l_size)=normalized radial shape function in reciprocal space
!!     %sij(lmn2_size)=nonlocal part of the overlap operator
!!     %tphitphj(:,:)=useful product tPhi(:,i)*tPhi(:,j)
!!
!! PARENTS
!!      bethe_salpeter,gstate,respfn,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

#include "../libpaw/libpaw.h"

subroutine abi_pawinit(gnt_option,gsqcut_eff,lcutdens,lmix,mpsang,nphi,nsym,ntheta,&
&                      pawang,pawrad,pawspnorb,pawtab,pawxcdev,xclevel,usepotzero)

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

 use m_paw_numeric, only: paw_spline
 use m_pawpsp, only: pawpsp_nl
 use m_paw_atom, only: atompaw_shpfun
 use m_pawang, only : pawang_type, pawang_init, pawang_free
 use m_pawrad, only : pawrad_type, simp_gen, nderiv_gen, poisson
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_pawinit'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: gnt_option,lcutdens,lmix,mpsang,nphi,nsym,ntheta
 integer,intent(in) :: pawspnorb,pawxcdev,xclevel,usepotzero
 real(dp),intent(in) :: gsqcut_eff
 type(pawang_type),intent(inout) :: pawang
!arrays
 type(pawrad_type),intent(in) :: pawrad(:)
 type(pawtab_type),target,intent(inout) :: pawtab(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: mqgrid_shp_default=300
 integer :: basis_size,i0lm,i0ln,ij_size,il,ilm,ilmn,iln,iq,isel,isel1
 integer :: itypat,j0lm,j0lmn,j0ln,jl,jlm,jlmn,jln,klm,klm1
 integer :: klmn,klmn1,kln,kln1,l_size,ll,lm0,lmax,lmax1,lmin,lmin1,lmn2_size
 integer :: lmn_size,lmnmix,mesh_size,meshsz,mm,ntypat,usexcnhat,use_ls_ylm,use_ylm
 real(dp) :: dq,gnrm,intg,ql,ql1,rg,rg1,vh1,yp1,ypn
 character(len=500) :: message
!arrays
 integer,allocatable :: indl(:,:),klm_diag(:),kmix_tmp(:)
 integer, pointer :: indlmn(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: ff(:),gg(:),hh(:),intvhatl(:)
 integer,allocatable :: indklmn_(:,:)
 real(dp),allocatable :: rad(:),rgl(:,:),vhatijl(:,:),vhatl(:),work(:)

!************************************************************************

 ntypat=size(pawtab)
 if (size(pawrad)/=ntypat) then
    write(message, "(A)") 'pawrad and pawtab should have the same size!'
   MSG_BUG(message)
 end if

 ! Immediately set the value of usepotzero
 ! it will be used later on in this subroutine
 pawtab%usepotzero=usepotzero

!==================================================
!1- INITIALIZE DATA RELATED TO ANGULAR MESH
!* ANGULAR GRID
!* REAL SPHERICAL HARMONICS
!* REAL GAUNT COEFFICIENTS

 use_ylm=0;if (pawxcdev==0) use_ylm=1
 use_ls_ylm=0;if (pawspnorb>0) use_ls_ylm=1
 call pawang_free(pawang)
 call pawang_init(pawang,gnt_option,mpsang-1,nphi,nsym,ntheta,pawxcdev,use_ls_ylm,use_ylm,xclevel)

 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)

!*******************
!Loop on atom types
!*******************
 do itypat=1,ntypat
   mesh_size=pawrad(itypat)%mesh_size
   l_size=pawtab(itypat)%l_size
   lmn_size=pawtab(itypat)%lmn_size
   lmn2_size=pawtab(itypat)%lmn2_size
   basis_size=pawtab(itypat)%basis_size
   ij_size=pawtab(itypat)%ij_size
   indlmn => pawtab(itypat)%indlmn(:,:)
   LIBPAW_ALLOCATE(indklmn_,(8,lmn2_size))
   LIBPAW_ALLOCATE(klm_diag,(lmn2_size))
   LIBPAW_ALLOCATE(ff,(mesh_size))
   LIBPAW_ALLOCATE(gg,(mesh_size))
   LIBPAW_ALLOCATE(hh,(mesh_size))
   LIBPAW_ALLOCATE(rad,(mesh_size))
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

   if (pawtab(itypat)%usexcnhat/=usexcnhat) then
     write(message, '(7a)' )&
&     'You cannot simultaneously use atomic data with different',ch10,&
&     'formulation of XC [using compensation charge in XC or not] !',ch10,&
&     'Action: change at least one of your atomic data (psp) file',ch10,&
&     '        or use usexcnhat keyword in input file.'
     MSG_ERROR(message)
   end if

!  ==================================================
!  2- TABULATE SHAPE FUNCTION

!  Allocated shape function
   if (pawtab(itypat)%shape_type/=-1) then
     if (allocated(pawtab(itypat)%shapefunc))  then
       LIBPAW_DEALLOCATE(pawtab(itypat)%shapefunc)
     end if
     LIBPAW_ALLOCATE(pawtab(itypat)%shapefunc,(mesh_size,l_size))
   else if (.not.allocated(pawtab(itypat)%shapefunc))  then
     message='shapefunc should be allocated with shape_type=-1'
     MSG_ERROR(message)
   end if
   if (allocated(pawtab(itypat)%gnorm))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%gnorm)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%gnorm,(l_size))

!  Compute shape function
   do il=1,l_size
     ll=il-1
     call atompaw_shpfun(ll,pawrad(itypat),gnrm,pawtab(itypat),ff)
     pawtab(itypat)%shapefunc(1:mesh_size,il)=ff(1:mesh_size)
     pawtab(itypat)%gnorm(il)=gnrm
   end do
!  In case of numerical shape function, compute some derivatives
   if (pawtab(itypat)%shape_type==-1) then
     if (allocated(pawtab(itypat)%dshpfunc))  then
       LIBPAW_DEALLOCATE(pawtab(itypat)%dshpfunc)
     end if
     LIBPAW_ALLOCATE(pawtab(itypat)%dshpfunc,(mesh_size,l_size,4))
     LIBPAW_ALLOCATE(work,(mesh_size))
     do il=1,l_size
       call nderiv_gen(pawtab(itypat)%dshpfunc(:,il,1),pawtab(itypat)%shapefunc(:,il),1,pawrad(itypat))
       yp1=pawtab(itypat)%dshpfunc(1,il,1);ypn=pawtab(itypat)%dshpfunc(mesh_size,il,1)
       call paw_spline(rad,pawtab(itypat)%shapefunc(:,il),mesh_size,yp1,ypn,pawtab(itypat)%dshpfunc(:,il,2))
       yp1=pawtab(itypat)%dshpfunc(1,il,2);ypn=pawtab(itypat)%dshpfunc(mesh_size,il,2)
       call paw_spline(rad,pawtab(itypat)%dshpfunc(:,il,1),mesh_size,yp1,ypn,pawtab(itypat)%dshpfunc(:,il,3))
       yp1=pawtab(itypat)%dshpfunc(1,il,3);ypn=pawtab(itypat)%dshpfunc(mesh_size,il,3)
       call paw_spline(rad,pawtab(itypat)%dshpfunc(:,il,2),mesh_size,yp1,ypn,pawtab(itypat)%dshpfunc(:,il,4))
     end do
     LIBPAW_DEALLOCATE(work)
   end if

!  In some cases, has to store radial mesh for shape function in pawtab variable
   if (pawtab(itypat)%shape_type==-1) then
     if (allocated(pawtab(itypat)%rad_for_spline))  then
       LIBPAW_DEALLOCATE(pawtab(itypat)%rad_for_spline)
     end if
     LIBPAW_ALLOCATE(pawtab(itypat)%rad_for_spline,(mesh_size))
     pawtab(itypat)%rad_for_spline(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)
   end if

!  In some cases, has to store shape function in reciprocal space
   if (pawtab(itypat)%has_shapefncg>0) then
     if (gsqcut_eff<tol8) then
       message='Computation of shapefncg only possible when gsqcut>0!'
       MSG_BUG(message)
     end if
     pawtab(itypat)%mqgrid_shp=mqgrid_shp_default
     if (allocated(pawtab(itypat)%shapefncg))  then
       LIBPAW_DEALLOCATE(pawtab(itypat)%shapefncg)
     end if
     if (allocated(pawtab(itypat)%qgrid_shp))  then
       LIBPAW_DEALLOCATE(pawtab(itypat)%qgrid_shp)
     end if
     LIBPAW_ALLOCATE(pawtab(itypat)%shapefncg,(pawtab(itypat)%mqgrid_shp,2,l_size))
     LIBPAW_ALLOCATE(pawtab(itypat)%qgrid_shp,(pawtab(itypat)%mqgrid_shp))
     dq=1.1_dp*sqrt(gsqcut_eff)/dble(pawtab(itypat)%mqgrid_shp-1)
     do iq=1,pawtab(itypat)%mqgrid_shp
       pawtab(itypat)%qgrid_shp(iq)=dble(iq-1)*dq
     end do
     LIBPAW_ALLOCATE(indl,(6,l_size))
     LIBPAW_ALLOCATE(rgl,(mesh_size,il))
     do il=1,l_size
       indl(:,il)=0;indl(1,il)=il-1;indl(5,il)=il
       rgl(1:mesh_size,il)=rad(1:mesh_size)*pawtab(itypat)%shapefunc(1:mesh_size,il)
     end do
     call pawpsp_nl(pawtab(itypat)%shapefncg,indl,l_size,l_size,&
&     pawtab(itypat)%mqgrid_shp,pawtab(itypat)%qgrid_shp,pawrad(itypat),rgl)
     pawtab(itypat)%shapefncg=four_pi*pawtab(itypat)%shapefncg
     LIBPAW_DEALLOCATE(indl)
     LIBPAW_DEALLOCATE(rgl)
   else
     pawtab(itypat)%mqgrid_shp=0
   end if

!  ==================================================
!  3- COMPUTE indklmn INDEXES GIVING klm, kln, abs(il-jl) and (il+jl), ilmn and jlmn
!  for each klmn=(ilmn,jlmn)

   if (allocated(pawtab(itypat)%indklmn))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%indklmn)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%indklmn,(8,lmn2_size))

   klm_diag=0
   do jlmn=1,lmn_size
     jl= indlmn(1,jlmn);jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     j0lm =jlm *(jlm -1)/2
     j0ln =jln *(jln -1)/2
     do ilmn=1,jlmn
       il= indlmn(1,ilmn);ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
       klmn=j0lmn+ilmn
       if (ilm<=jlm) then
         indklmn_(1,klmn)=j0lm+ilm
       else
         i0lm=ilm*(ilm-1)/2
         indklmn_(1,klmn)=i0lm+jlm
       end if
       if (iln<=jln) then
         indklmn_(2,klmn)=j0ln+iln
       else
         i0ln=iln*(iln-1)/2
         indklmn_(2,klmn)=i0ln+jln
       end if
       indklmn_(3,klmn)=min(abs(il-jl),lcutdens)
       indklmn_(4,klmn)=min(il+jl,lcutdens)
       indklmn_(5,klmn)=ilm
       indklmn_(6,klmn)=jlm
       indklmn_(7,klmn)=ilmn
       indklmn_(8,klmn)=jlmn
       pawtab(itypat)%indklmn(:,klmn)=indklmn_(:,klmn)
       if (ilm==jlm) klm_diag(klmn)=1
     end do
   end do

!  ==================================================
!  4- COMPUTE various FACTORS/SIZES (depending on (l,m,n))

   pawtab(itypat)%lcut_size=min(l_size,lcutdens+1)

   if (allocated(pawtab(itypat)%dltij))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%dltij)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%dltij,(lmn2_size))
   pawtab(itypat)%dltij(:)=two
   do ilmn=1,lmn_size
     pawtab(itypat)%dltij(ilmn*(ilmn+1)/2)=one
   end do

   lmnmix=int(zero)
   LIBPAW_ALLOCATE(kmix_tmp,(lmn2_size))
   do jlmn=1,lmn_size
     jl=indlmn(1,jlmn)
     if (jl<=lmix) then
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         il=indlmn(1,ilmn)
         if (il<=lmix) then
           lmnmix=lmnmix+1
           kmix_tmp(lmnmix)=j0lmn+ilmn
         end if
       end do
     end if
   end do
   if (allocated(pawtab(itypat)%kmix))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%kmix)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%kmix,(lmnmix))
   pawtab(itypat)%lmnmix_sz=lmnmix
   pawtab(itypat)%kmix(1:lmnmix)=kmix_tmp(1:lmnmix)
   LIBPAW_DEALLOCATE(kmix_tmp)

!  ==================================================
!  5- COMPUTE Qijl TERMS AND Sij MATRIX

!  Store some usefull quantities
   if (allocated(pawtab(itypat)%phiphj))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%phiphj)
   end if
   if (allocated(pawtab(itypat)%tphitphj))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%tphitphj)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%phiphj,(mesh_size,ij_size))
   LIBPAW_ALLOCATE(pawtab(itypat)%tphitphj,(mesh_size,ij_size))
   do jln=1,basis_size
     j0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=j0ln+iln
       pawtab(itypat)%phiphj  (1:mesh_size,kln)=pawtab(itypat)%phi (1:mesh_size,iln)&
&       *pawtab(itypat)%phi (1:mesh_size,jln)
       pawtab(itypat)%tphitphj(1:mesh_size,kln)=pawtab(itypat)%tphi(1:mesh_size,iln)&
&       *pawtab(itypat)%tphi(1:mesh_size,jln)
     end do
   end do

!  Compute q_ijL and S_ij=q_ij0
   if (allocated(pawtab(itypat)%qijl))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%qijl)
   end if
   if (allocated(pawtab(itypat)%sij))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%sij)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%qijl,(l_size*l_size,lmn2_size))
   LIBPAW_ALLOCATE(pawtab(itypat)%sij,(lmn2_size))
   pawtab(itypat)%qijl=zero
   pawtab(itypat)%sij=zero
   do klmn=1,lmn2_size
     klm=indklmn_(1,klmn);kln=indklmn_(2,klmn)
     lmin=indklmn_(3,klmn);lmax=indklmn_(4,klmn)
     do ll=lmin,lmax,2
       lm0=ll*ll+ll+1;ff(1)=zero
       ff(2:mesh_size)=(pawtab(itypat)%phiphj  (2:mesh_size,kln)&
&       -pawtab(itypat)%tphitphj(2:mesh_size,kln))&
&       *rad(2:mesh_size)**ll
       call simp_gen(intg,ff,pawrad(itypat))
       do mm=-ll,ll
         isel=pawang%gntselect(lm0+mm,klm)
         if (isel>0) pawtab(itypat)%qijl(lm0+mm,klmn)=intg*pawang%realgnt(isel)
       end do
     end do
     if (klm_diag(klmn)==1) pawtab(itypat)%sij(klmn)= &
&     pawtab(itypat)%qijl(1,klmn)*sqrt(four_pi)
   end do

!  ==================================================
!  6- COMPUTE Eijkl TERMS

   if (allocated(pawtab(itypat)%eijkl))  then
     LIBPAW_DEALLOCATE(pawtab(itypat)%eijkl)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%eijkl,(lmn2_size,lmn2_size))

!  Compute:
!  vhatL(r) according to eq. (A14) in Holzwarth et al., PRB 55, 2005 (1997)
!  intvhatL=$\int_{0}^{r_c}{vhatL(r) shapefunc_L(r) r^2\,dr}$
!  vhatijL =$\int_{0}^{r_c}{vhatL(r) \tilde{\phi}_i \tilde{\phi}_j \,dr}$
!  -----------------------------------------------------------------
   LIBPAW_ALLOCATE(vhatl,(mesh_size))
   LIBPAW_ALLOCATE(vhatijl,(lmn2_size,l_size))
   LIBPAW_ALLOCATE(intvhatl,(l_size))
   intvhatl(:)=zero;vhatl(:)=zero;vhatijl(:,:)=zero
   do il=1,l_size
     vhatl(1)=zero;ff(1)=zero
     ff(2:mesh_size)=pawtab(itypat)%shapefunc(2:mesh_size,il)*rad(2:mesh_size)**2
     call poisson(ff,il-1,intg,pawrad(itypat),vhatl)
     vhatl(2:mesh_size)=two*vhatl(2:mesh_size)/rad(2:mesh_size)
     gg(1:mesh_size)=vhatl(1:mesh_size)*ff(1:mesh_size)
     call simp_gen(intvhatl(il),gg,pawrad(itypat))
     do klmn=1,lmn2_size
       kln=indklmn_(2,klmn)
       hh(1:mesh_size)=vhatl(1:mesh_size)*pawtab(itypat)%tphitphj(1:mesh_size,kln)
       call simp_gen(vhatijl(klmn,il),hh,pawrad(itypat))
     end do
   end do
   LIBPAW_DEALLOCATE(vhatl)

!  Compute:
!  eijkl=$ vh1_ijkl - Vhatijkl - Bijkl - Cijkl$
!  With:
!$vh1_ijkl =\sum_{L,m} {vh1*Gaunt(i,j,Lm)*Gaunt(k,l,Lm)}$
!$Vhat_ijkl=\sum_{L,m} {vhatijL*Gaunt(i,j,Lm)*q_klL}$
!$B_ijkl   =\sum_{L,m} {vhatijL*Gaunt(k,l,Lm)*q_ijL}$
!$C_ijkl   =\sum_{L,m} {intvhatL*q_ijL*q_klL}$
!  and:
!  vh1 according to eq. (A17) in Holzwarth et al., PRB 55, 2005 (1997)
!  Warning: compute only eijkl for (i,j)<=(k,l)
!  -----------------------------------------------------------------
   pawtab(itypat)%eijkl(:,:)=zero
   meshsz=pawrad(itypat)%int_meshsz;if (meshsz>mesh_size) ff(meshsz+1:mesh_size)=zero
   do klmn=1,lmn2_size
     klm=indklmn_(1,klmn);kln=indklmn_(2,klmn)
     lmin=indklmn_(3,klmn);lmax=indklmn_(4,klmn)
     do ll=lmin,lmax,2
       lm0=ll*ll+ll+1
       ff(1:meshsz)=pawtab(itypat)%phiphj  (1:meshsz,kln)
       call poisson(ff,ll,intg,pawrad(itypat),gg)
       ff(1:meshsz)=pawtab(itypat)%tphitphj(1:meshsz,kln)
       call poisson(ff,ll,intg,pawrad(itypat),hh)
       do klmn1=klmn,lmn2_size
         klm1=indklmn_(1,klmn1);kln1=indklmn_(2,klmn1)
         lmin1=indklmn_(3,klmn1);lmax1=indklmn_(4,klmn1)
         vh1=zero
         if ((ll.ge.lmin1).and.(ll.le.lmax1)) then
           ff(1)=zero
           ff(2:meshsz)=(pawtab(itypat)%phiphj  (2:meshsz,kln1)*gg(2:meshsz)&
&           -pawtab(itypat)%tphitphj(2:meshsz,kln1)*hh(2:meshsz))&
&           *two/rad(2:meshsz)
           call simp_gen(vh1,ff,pawrad(itypat))
         end if
         do mm=-ll,ll
           isel =pawang%gntselect(lm0+mm,klm)
           isel1=pawang%gntselect(lm0+mm,klm1)
           if (isel>0.and.isel1>0) then
             rg =pawang%realgnt(isel)
             rg1=pawang%realgnt(isel1)
             ql =pawtab(itypat)%qijl(lm0+mm,klmn)
             ql1=pawtab(itypat)%qijl(lm0+mm,klmn1)
             pawtab(itypat)%eijkl(klmn,klmn1)=pawtab(itypat)%eijkl(klmn,klmn1)&
&             +(   vh1                *rg *rg1&      ! vh1_ijkl
&            -    vhatijl(klmn ,ll+1)*rg *ql1&      ! Vhat_ijkl
&            -    vhatijl(klmn1,ll+1)*rg1*ql &      ! B_ijkl
&            -    intvhatl(ll+1)     *ql *ql1&      ! C_ijkl
&            )*two_pi
           end if
         end do
       end do
     end do
   end do
   LIBPAW_DEALLOCATE(vhatijl)
   LIBPAW_DEALLOCATE(intvhatl)
!  ==================================================
!  7- COMPUTE gamma_ij TERMS
!  corrections to get the background right

   if (pawtab(itypat)%usepotzero==1) then 
     if (allocated(pawtab(itypat)%gammaij))  then
       LIBPAW_DEALLOCATE(pawtab(itypat)%gammaij)
     end if
     LIBPAW_ALLOCATE(pawtab(itypat)%gammaij,(lmn2_size))

     LIBPAW_ALLOCATE(work,(mesh_size))
     do klmn=1,lmn2_size
       if (klm_diag(klmn)==1) then
         kln=indklmn_(2,klmn)

         ff(1)=zero
         ff(2:mesh_size)=pawtab(itypat)%phiphj(2:mesh_size,kln)-pawtab(itypat)%tphitphj(2:mesh_size,kln)
         ! First, compute q_ij^00
         call simp_gen(intg,ff,pawrad(itypat))
         ! Second, compute phi^2 - tphi^2 - 4pi*r^2 q_ij^00 g_0(r) 
         ff(2:mesh_size)= ff(2:mesh_size) - intg*pawtab(itypat)%shapefunc(2:mesh_size,1)*rad(2:mesh_size)**2 

         call poisson(ff,0,intg,pawrad(itypat),work)

         ! work is r*vh
         ! it should be then multiplied by r to prepare the integration in the
         ! sphere
         work(1)=zero
         work(2:mesh_size)=work(2:mesh_size)*rad(2:mesh_size)
         ! Third, average it over the sphere
         call simp_gen(intg,work,pawrad(itypat))
         ! Finally, store it in pawtab%gammaij
         pawtab(itypat)%gammaij(klmn)=intg*four_pi
       else
         pawtab(itypat)%gammaij(klmn)=zero
       end if
     end do
     LIBPAW_DEALLOCATE(work)
   end if


!  ***********************
!  End Loop on atom types
!  ***********************
   LIBPAW_DEALLOCATE(ff)
   LIBPAW_DEALLOCATE(gg)
   LIBPAW_DEALLOCATE(hh)
   LIBPAW_DEALLOCATE(indklmn_)
   LIBPAW_DEALLOCATE(klm_diag)
   LIBPAW_DEALLOCATE(rad)
 end do

end subroutine abi_pawinit
!!***

