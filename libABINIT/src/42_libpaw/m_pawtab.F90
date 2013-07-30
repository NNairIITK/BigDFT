!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawtab
!! NAME
!!  m_pawtab
!!
!! FUNCTION
!!  This module contains the definition of the pawtab_type structured datatype,
!!  as well as related functions and methods.
!!  pawtab_type variables define TABulated data for PAW (from pseudopotential)
!!
!! COPYRIGHT
!! Copyright (C) 2013-2013 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"

MODULE m_pawtab

 use defs_basis
! use m_errors
use interfaces_12_hide_mpi
use interfaces_14_hidewrite
use interfaces_16_hideleave
 use m_profiling
 use m_xmpi

 implicit none

 private

!public procedures.
 public :: pawtab_nullify
 public :: pawtab_nullify_array
 public :: pawtab_destroy
 public :: pawtab_destroy_array
 public :: pawtab_print
 public :: pawtab_bcast
 public :: wvlpaw_destroy
 public :: wvlpaw_nullify
!!***

!----------------------------------------------------------------------

!!****t* m_pawtab/rholoc_type
!! NAME
!! rholoc_type
!!
!! FUNCTION
!! Objects for WVL+PAW
!!
!! SOURCE

 type,public :: wvlpaw_rholoc_type

  integer :: msz           ! mesh size 
  real(dp),pointer :: d(:,:) ! local rho and derivatives
  real(dp),pointer :: rad(:) ! radial mesh 

 end type wvlpaw_rholoc_type
!!***

!----------------------------------------------------------------------

!!****t* m_pawtab/wvlpaw_type
!! NAME
!! wvlpaw_type
!!
!! FUNCTION
!! Objects for WVL+PAW 
!!
!! SOURCE

 type,public :: wvlpaw_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: npspcode_init_guess
   ! This is for the PAW-WVL case, only for the initial guess

  integer :: ptotgau
   ! total number of complex gaussians 
   ! for tproj 

  integer,pointer :: pngau(:)
   ! number of complex gaussians per basis element
   ! for tproj 

!Real pointers

  real(dp),pointer :: parg(:,:)
   !argument of Gaussians

  real(dp),pointer :: pfac(:,:)
   !factors of Gaussians

!Other pointers

  type(wvlpaw_rholoc_type) :: rholoc
   ! local density

 end type wvlpaw_type
!!***

!----------------------------------------------------------------------

!!****t* m_pawtab/pawtab_type
!! NAME
!! pawtab_type
!!
!! FUNCTION
!! This structured datatype contains TABulated data for PAW (from pseudopotential)
!! used in PAW calculations.
!!
!! SOURCE

 type,public :: pawtab_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: basis_size
   ! Number of elements for the paw nl basis on the considered atom type

  integer :: has_kij
   ! if 1, onsite matrix elements of the kinetic operator are allocated
   ! if 2, onsite matrix elements of the kinetic operator are computed and stored

  integer :: has_nabla
   ! if 1, onsite matrix elements of the nabla operator are allocated.
   ! if 2, onsite matrix elements of the nabla operator are computed and stored.

  integer :: has_vhntzc
   ! if 1, space for vhntzc is allocated
   ! if 2, vhntzc has been read from PAW file and stored

  integer :: has_vhnzc
   ! if 1, space for vhnzc is allocated
   ! if 2, vhnzc has been computed and stored

  integer :: ij_proj
   ! Number of (i,j) elements for the orbitals on which U acts (PAW+U only)
   ! on the considered atom type (ij_proj=1 (1 projector), 3 (2 projectors)...)
   ! Also used for local exact-exchange

  integer :: ij_size
   ! Number of (i,j) elements for the symetric paw basis
   ! on the considered atom type (ij_size=basis_size*(basis_size+1)/2)

  integer :: lcut_size
   ! Maximum value of l+1 leading to non zero Gaunt coeffs
   ! modified by dtset%pawlcutd
   ! lcut_size=min(2*l_max,dtset%pawlcutd)+1

  integer :: l_size
   ! Maximum value of l+1 leading to non zero Gaunt coeffs
   ! l_size=2*l_max-1

  integer :: lexexch
   ! lexexch gives l on which local exact-exchange is applied for a given type of atom.

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz
   ! lmnmix_sz=number of klmn=(lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix

  integer :: lpawu
   ! lpawu gives l on which U is applied for a given type of atom.

  integer :: nproju
   ! nproju is the number of projectors for orbitals on which paw+u acts.
   ! Also used for local exact-exchange

  integer :: mesh_size
   ! Dimension of radial mesh

  integer :: core_mesh_size
   ! Dimension of radial mesh for core density

  integer :: tnvale_mesh_size
   ! Dimension of radial mesh for tnvale

  integer :: mqgrid
   ! Number of points in the reciprocal space grid on which
   ! the radial functions (tcorespl, tvalespl...) are specified
   ! Same as psps%mqgrid_vl

  integer :: mqgrid_shp
   ! Number of points in the reciprocal space grid on which
   ! the radial shape functions (shapefncg) are given

  integer :: shape_lambda
   ! Lambda parameter in gaussian shapefunction (shape_type=2)

  integer :: shape_type
   ! Radial shape function type
   ! shape_type=-1 ; g(r)=numeric (read from psp file)
   ! shape_type= 1 ; g(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2 if r<=rshp, zero if r>rshp
   ! shape_type= 2 ; g(r)=exp[-(r/sigma)**lambda]
   ! shape_type= 3 ; gl(r)=Alpha(1,l)*jl(q(1,l)*r)+Alpha(2,l)*jl(q(2,l)*r) for each l

  integer :: useexexch
   ! useexexch=0 ; do not use local exact-exchange
   ! useexexch=1 ; use local exact-exchange

  integer :: usepawu
   ! usepawu=0 ; do not use PAW+U formalism
   ! usepawu=1 ; use PAW+U formalism (Full localized limit)
   ! usepawu=2 ; use PAW+U formalism (Around Mean Field)

  integer :: usetcore
   ! Flag controling use of pseudized core density (0 if tncore=zero)

  integer :: usetvale
   ! Flag controling use of pseudized valence density (0 if tnval is unknown)

  integer :: usexcnhat
   ! 0 if compensation charge density is not included in XC terms
   ! 1 if compensation charge density is included in XC terms

!Real (real(dp)) scalars

  real(dp) :: dncdq0
   ! Gives 1/q d(tNcore(q))/dq for q=0
   ! (tNcore(q) = FT of pseudo core density)

  real(dp) :: d2ncdq0
   ! Gives contribution of d2(tNcore(q))/d2q for q=0
   ! \int{(16/15)*pi^5*n(r)*r^6* dr}
   ! (tNcore(q) = FT of pseudo core density)

  real(dp) :: dnvdq0
   ! Gives 1/q d(tNvale(q))/dq for q=0
   ! (tNvale(q) = FT of pseudo valence density)

  real(dp) :: exccore
   ! Exchange-correlation energy for the core density

  real(dp) :: exchmix
   ! mixing of exact exchange; default is 0.25 (PBE0)

  real(dp) :: f4of2_sla
   ! Ratio of Slater Integrals F4 and F2

  real(dp) :: f6of2_sla
   ! Ratio of Slater Integrals F6 and F4

  real(dp) :: jpawu
   ! jpawu
   ! Value of J parameter for paw+u for a given type.

  real(dp) :: rpaw
   ! Radius of PAW sphere

  real(dp) :: rshp
   ! Compensation charge radius (if r>rshp, g(r)=zero)

  real(dp) :: rcore
   ! Radius of core corrections (rcore >= rpaw)

  real(dp) :: shape_sigma
   ! Sigma parameter in gaussian shapefunction (shape_type=2)

  real(dp) :: upawu
   ! upawu
   ! Value of U parameter for paw+u for a given type.

!Objects
  type(wvlpaw_type) :: wvl
   !variable containing objects needed
   !for wvl+paw implementation


!Integer arrays

  integer, pointer :: indklmn(:,:)
   ! indklmn(6,lmn2_size)
   ! Array giving klm, kln, abs(il-jl), (il+jl), ilm and jlm for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn

  integer, pointer :: klmntomn(:,:)
   ! klmntomn(4,lmn2_size)
   ! Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn
   ! NB: klmntomn is an application and not a bijection

  integer, pointer :: kmix(:)
   ! kmix(lmnmix_sz)
   ! Indirect array selecting the klmn=(lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix

  integer, pointer :: lnproju(:)
   ! lnproju(nproju) gives ln (index for phi) for each projectors on which U acts (PAW+U only)
   ! nproju is 1 or 2 and  is the number of projectors for correlated orbitals
   ! Also used for local exact-exchange

  integer, pointer :: orbitals(:)
   ! gives the l quantum number per basis element

!Real (real(dp)) arrays

  real(dp), pointer :: coredens(:)
   ! coredens(mesh_size)
   ! Gives the core density of the atom

  real(dp), pointer :: dij0(:)
   ! dij0(lmn2_size)
   ! Part of the Dij term (non-local operator) completely
   ! calculated in the atomic data part

  real(dp), pointer :: dltij(:)
   ! dltij(lmn2_size)
   ! Factor used to compute sums over klmn=(ilmn,jlmn)
   ! ((ilmn,ilmn) term has to be added once)
   ! dltij(klmn)=1 if ilmn=jlmn, else dltij(klmn)=2

  real(dp), pointer :: dshpfunc(:,:,:)
   ! shapefunc(mesh_size,l_size,4)
   ! Gives the 4 first derivatives of  radial shape function
   ! for each l component; used only if shape_type=-1

  real(dp), pointer :: eijkl(:,:)
   ! eijkl(lmn2_size,lmn2_size)
   ! Part of the Dij term (non-local operator) that depends only from
   ! the projected occupation coeffs in the self-consistent loop

  real(dp), pointer :: fk(:,:)
   ! fk(6,4)
   ! Slater integrals used for local exact exchange

  real(dp), pointer :: gnorm(:)
   ! gnorm(l_size)
   ! Give the the normalization factor of each radial shape function

  real(dp),pointer :: kij(:)
   ! Onsite matrix elements <phi|\kinetic|phj>-<tphi|\kinetic|tphj>

  real(dp),pointer :: nabla_ij(:,:,:)
   ! nabla_ij(3,lmn_size,lmn_size))
   ! Onsite matrix elements <phi|\nabla|phj>-<tphi|\nabla|tphj>

  real(dp), pointer :: phi(:,:)
   ! phi(mesh_size, basis_size)
   ! Gives the paw electron wavefunctions on the radial grid

  real(dp), pointer :: phiphj(:,:)
   ! phiphj(mesh_size,ij_size)
   ! Useful product Phi(:,i)*Phi(:,j)

  real(dp), pointer :: phiphjint(:)
   ! phiphjint(ij_proj)
   ! Integration of Phi(:,i)*Phi(:,j) for LDA+U/local exact-exchange occupation matrix

  real(dp), pointer :: ph0phiint(:)
   ! ph0phjint(ij_proj)
   ! Integration of Phi(:,1)*Phi(:,j) for LDA+DMFT projections

  real(dp), pointer :: qgrid_shp(:)
   ! qgrid_shp(mqgrid_shp)
   ! Grid of points in reciprocal space on which the shape functions are given

  real(dp), pointer :: qijl(:,:)
   ! qijl(l_size**2,lmn2_size)
   ! The qijl are the moments of the charge density difference between
   ! the AE and PS partial wave for each channel (i,j). They take part
   ! to the building of the compensation charge

  real(dp), pointer :: rad_for_spline(:)
   ! rad_for_spline(mesh_size)
   ! Radial mesh used to spline quantities on radial mesh;
   ! Allocated and used only when
   !     shape_type=-1 (numerical shape function)
   !  or usedvloc=1 (use of vloc derivative)

  real(dp), pointer :: rhoij0(:)
   ! rhoij0(lmn2_size)
   ! Initial guess for rhoij

  real(dp), pointer :: shape_alpha(:,:)
   ! shape_alpha(2,l_size)
   ! Alpha_i parameters in Bessel shapefunctions (shape_type=3)

  real(dp), pointer :: shape_q(:,:)
   ! shape_q(2,l_size)
   ! Q_i parameters in Bessel shapefunctions (shape_type=3)

  real(dp), pointer :: shapefunc(:,:)
   ! shapefunc(mesh_size,l_size)
   ! Gives the normalized radial shape function for each l component

  real(dp), pointer :: shapefncg(:,:,:)
   ! shapefncg(mqgrid_shp,2,l_size)
   ! Gives the spherical Fourier transform of the radial shape function
   ! for each l component (for each qgrid_shp(i)) + second derivative

  real(dp), pointer :: sij(:)
   ! sij(lmn2_size)
   ! Nonlocal part of the overlap operator

  real(dp), pointer :: tcoredens(:,:)
   ! tcoredens(mesh_size,1)
   ! Gives the pseudo core density of the atom
   ! In PAW+WVL:
   !  tcoredens(mesh_size,2:6) 
   !  are the first to the fifth derivatives of the pseudo core density.

  real(dp), pointer :: tcorespl(:,:)
   ! tcorespl(mqgrid,2)
   ! Gives the pseudo core density in reciprocal space on a regular grid

  real(dp), pointer :: tphi(:,:)
   ! tphi(mesh_size,basis_size)
   ! Gives, on the radial grid, the paw atomic pseudowavefunctions

  real(dp), pointer :: tphitphj(:,:)
   ! tphitphj(mesh_size,ij_size)
   ! Useful product tPhi(:,i)*tPhi(:,j)

  real(dp), pointer :: tproj(:,:)
   ! non-local projectors

  real(dp), pointer :: tvalespl(:,:)
   ! tvalespl(mqgrid,2)
   ! Gives the pseudo valence density in reciprocal space on a regular grid

  real(dp), pointer :: Vee(:,:,:,:)
   ! PAW+U:
   ! Screened interaction matrix deduced from U and J parameters
   ! computed on the basis of orbitals on which U acts.

  real(dp), pointer :: Vex(:,:,:,:,:)
   ! Local exact-exchange:
   ! Screened interaction matrix deduced from calculation of Slater integrals
   ! computed on the basis of orbitals on which local exact exchange acts.

  real(dp), pointer :: VHntZC(:)
   ! VHntZC(mesh_size)
   ! Hartree potential for pseudized Zc density, v_H[\tilde{n}_{Zc}]
   ! read in from PAW file

  real(dp), pointer :: VHnZC(:)
   ! VHnZC(mesh_size)
   ! Hartree potential for Zc density, v_H[n_{Zc}]
   ! constructed from core density in PAW file (see psp7in.F90)

  real(dp), pointer :: zioneff(:)
   ! zioneff(ij_proj)
   ! "Effective charge"*n "seen" at r_paw, deduced from Phi at r_paw, n:
   ! pricipal quantum number
   ! good approximation to model wave function outside PAW-sphere through

 end type pawtab_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_nullify
!! NAME
!!  pawtab_nullify
!!
!! FUNCTION
!!  Nullify pointers and flags in a pawtab structure
!!
!! SIDE EFFECTS
!!  Pawtab<type(pawtab_type)>=PAW arrays tabulated.
!!                            Nullified in output
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_nullify(Pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtab_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawtab_type),intent(inout) :: Pawtab

!Local variables-------------------------------

! *************************************************************************

 !@Pawtab_type

 nullify(Pawtab%indklmn)
 nullify(Pawtab%klmntomn)
 nullify(Pawtab%kmix)
 nullify(Pawtab%lnproju)
 nullify(Pawtab%coredens)
 nullify(Pawtab%dij0)
 nullify(Pawtab%dltij)
 nullify(Pawtab%dshpfunc)
 nullify(Pawtab%eijkl)
 nullify(Pawtab%fk)
 nullify(Pawtab%gnorm)
 nullify(Pawtab%kij)
 nullify(Pawtab%nabla_ij)
 nullify(Pawtab%orbitals)
 nullify(Pawtab%phi)
 nullify(Pawtab%phiphj)
 nullify(Pawtab%phiphjint)
 nullify(Pawtab%ph0phiint)
 nullify(Pawtab%qgrid_shp)
 nullify(Pawtab%qijl)
 nullify(Pawtab%rad_for_spline)
 nullify(Pawtab%rhoij0)
 nullify(Pawtab%shape_alpha)
 nullify(Pawtab%shape_q)
 nullify(Pawtab%shapefunc)
 nullify(Pawtab%shapefncg)
 nullify(Pawtab%sij)
 nullify(Pawtab%tcoredens)
 nullify(Pawtab%tcorespl)
 nullify(Pawtab%tphi)
 nullify(Pawtab%tphitphj)
 nullify(Pawtab%tproj)
 nullify(Pawtab%tvalespl)
 nullify(Pawtab%Vee)
 nullify(Pawtab%Vex)
 nullify(Pawtab%VHntZC)
 nullify(Pawtab%VHnZC)
 nullify(Pawtab%wvl%parg)
 nullify(Pawtab%wvl%pfac)
 nullify(Pawtab%wvl%pngau)
 nullify(Pawtab%wvl%rholoc%d)
 nullify(Pawtab%wvl%rholoc%rad)
 nullify(Pawtab%zioneff)

 call wvlpaw_nullify(Pawtab%wvl)

 ! === Reset all flags and sizes ===
 Pawtab%has_kij=0
 Pawtab%has_nabla=0
 Pawtab%has_vhntzc=0
 Pawtab%has_vhnzc=0
 Pawtab%usetcore=0
 Pawtab%usetvale=0
 Pawtab%usexcnhat=0
 Pawtab%useexexch=0
 Pawtab%usepawu=0
 Pawtab%mqgrid=0
 Pawtab%mqgrid_shp=0

 Pawtab%basis_size=0
 Pawtab%ij_proj=0
 Pawtab%ij_size=0
 Pawtab%lcut_size=0
 Pawtab%l_size=0
 Pawtab%lexexch=-1
 Pawtab%lmn_size=0
 Pawtab%lmn2_size=0
 Pawtab%lmnmix_sz=0
 Pawtab%lpawu=-1
 Pawtab%nproju=0
 Pawtab%mesh_size=0
 Pawtab%core_mesh_size=0
 Pawtab%tnvale_mesh_size=0
 Pawtab%shape_type=-10

end subroutine pawtab_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_pawtap/pawtab_nullify_array
!! NAME
!!  pawtab_nullify_array
!!
!! FUNCTION
!!  Nullify all pointers in an array of pawtab data structures
!!
!! PARENTS
!!      driver,mblktyp1,mblktyp5,rdddb9,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_nullify_array(Pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtab_nullify_array'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ii,nn

! *************************************************************************

 !@pawtab_type

 nn=size(Pawtab)
 if (nn==0) return

 do ii=1,nn
   call pawtab_nullify(Pawtab(ii))
 end do

end subroutine pawtab_nullify_array
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_destroy
!! NAME
!!  pawtab_destroy
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a pawtab structure
!!
!! SIDE EFFECTS
!!  Pawtab<type(pawtab_type)>=PAW arrays tabulated.
!!  All associated pointers in Pawtab are deallocated
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_destroy(Pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtab_destroy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawtab_type),intent(inout) :: Pawtab

!Local variables-------------------------------

! *************************************************************************

 !@Pawtab_type

 if (associated(Pawtab%indklmn))  then
  ABI_DEALLOCATE(Pawtab%indklmn)
 end if
 if (associated(Pawtab%klmntomn))  then
   ABI_DEALLOCATE(Pawtab%klmntomn)
 end if
 if (associated(Pawtab%kmix))  then
   ABI_DEALLOCATE(Pawtab%kmix)
 end if
 if (associated(Pawtab%lnproju))  then
   ABI_DEALLOCATE(Pawtab%lnproju)
 end if
 if (associated(Pawtab%coredens))  then
   ABI_DEALLOCATE(Pawtab%coredens)
 end if
 if (associated(Pawtab%dij0))  then
   ABI_DEALLOCATE(Pawtab%dij0)
 end if
 if (associated(Pawtab%dltij))  then
   ABI_DEALLOCATE(Pawtab%dltij)
 end if
 if (associated(Pawtab%dshpfunc))  then
   ABI_DEALLOCATE(Pawtab%dshpfunc)
 end if
 if (associated(Pawtab%eijkl))  then
   ABI_DEALLOCATE(Pawtab%eijkl)
 end if
 if (associated(Pawtab%fk))  then
   ABI_DEALLOCATE(Pawtab%fk)
 end if
 if (associated(Pawtab%gnorm))  then
   ABI_DEALLOCATE(Pawtab%gnorm)
 end if
 if (associated(Pawtab%kij))  then
   ABI_DEALLOCATE(Pawtab%kij)
 end if
 if (associated(Pawtab%nabla_ij))  then
   ABI_DEALLOCATE(Pawtab%nabla_ij)
 end if
 if (associated(Pawtab%orbitals)) then
   ABI_DEALLOCATE(Pawtab%orbitals)
 end if
 if (associated(Pawtab%phi))  then
   ABI_DEALLOCATE(Pawtab%phi)
 end if
 if (associated(Pawtab%phiphj))  then
   ABI_DEALLOCATE(Pawtab%phiphj)
 end if
 if (associated(Pawtab%phiphjint))  then
   ABI_DEALLOCATE(Pawtab%phiphjint)
 end if
 if (associated(Pawtab%ph0phiint))  then
   ABI_DEALLOCATE(Pawtab%ph0phiint)
 end if
 if (associated(Pawtab%qgrid_shp))  then
   ABI_DEALLOCATE(Pawtab%qgrid_shp)
 end if
 if (associated(Pawtab%qijl))  then
   ABI_DEALLOCATE(Pawtab%qijl)
 end if
 if (associated(Pawtab%rad_for_spline))  then
   ABI_DEALLOCATE(Pawtab%rad_for_spline)
 end if
 if (associated(Pawtab%rhoij0))  then
   ABI_DEALLOCATE(Pawtab%rhoij0)
 end if
 if (associated(Pawtab%shape_alpha))  then
   ABI_DEALLOCATE(Pawtab%shape_alpha)
 end if
 if (associated(Pawtab%shape_q))  then
   ABI_DEALLOCATE(Pawtab%shape_q)
 end if
 if (associated(Pawtab%shapefunc))  then
   ABI_DEALLOCATE(Pawtab%shapefunc)
 end if
 if (associated(Pawtab%shapefncg))  then
   ABI_DEALLOCATE(Pawtab%shapefncg)
 end if
 if (associated(Pawtab%sij))  then
   ABI_DEALLOCATE(Pawtab%sij)
 end if
 if (associated(Pawtab%tcoredens))  then
   ABI_DEALLOCATE(Pawtab%tcoredens)
 end if
 if (associated(Pawtab%tcorespl))  then
   ABI_DEALLOCATE(Pawtab%tcorespl)
 end if
 if (associated(Pawtab%tphi))  then
   ABI_DEALLOCATE(Pawtab%tphi)
 end if
 if (associated(Pawtab%tphitphj))  then
   ABI_DEALLOCATE(Pawtab%tphitphj)
 end if
 if (associated(Pawtab%tproj)) then
   ABI_DEALLOCATE(Pawtab%tproj)
 end if
 if (associated(Pawtab%tvalespl))  then
   ABI_DEALLOCATE(Pawtab%tvalespl)
 end if
 if (associated(Pawtab%Vee))  then
   ABI_DEALLOCATE(Pawtab%Vee)
 end if
 if (associated(Pawtab%Vex))  then
   ABI_DEALLOCATE(Pawtab%Vex)
 end if
 if (associated(Pawtab%VHntZC))  then
   ABI_DEALLOCATE(Pawtab%VHntZC)
 end if
 if (associated(Pawtab%VHnZC))  then
   ABI_DEALLOCATE(Pawtab%VHnZC)
 end if
 if (associated(Pawtab%zioneff))  then
   ABI_DEALLOCATE(Pawtab%zioneff)
 end if

 call wvlpaw_destroy(Pawtab%wvl)

 ! === Reset all flags and sizes ===
 Pawtab%has_kij=0
 Pawtab%has_nabla=0
 Pawtab%has_vhntzc=0
 Pawtab%has_vhnzc=0
 Pawtab%usetcore=0
 Pawtab%usetvale=0
 Pawtab%usexcnhat=0
 Pawtab%useexexch=0
 Pawtab%usepawu=0
 Pawtab%mqgrid=0
 Pawtab%mqgrid_shp=0

 Pawtab%basis_size=0
 Pawtab%ij_proj=0
 Pawtab%ij_size=0
 Pawtab%lcut_size=0
 Pawtab%l_size=0
 Pawtab%lexexch=-1
 Pawtab%lmn_size=0
 Pawtab%lmn2_size=0
 Pawtab%lmnmix_sz=0
 Pawtab%lpawu=-1
 Pawtab%nproju=0
 Pawtab%mesh_size=0
 Pawtab%core_mesh_size=0
 Pawtab%tnvale_mesh_size=0
 Pawtab%shape_type=-10
 
end subroutine pawtab_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_pawtap/pawtab_destroy_array
!! NAME
!!  pawtab_destroy_array
!!
!! FUNCTION
!!  Destroy (deallocate) all pointers in an array of pawtab data structures
!!
!! PARENTS
!!      driver,mblktyp1,mblktyp5,rdddb9,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_destroy_array(Pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtab_destroy_array'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(pawtab_type),intent(inout) :: Pawtab(:)

!Local variables-------------------------------
 integer :: ii,nn

! *************************************************************************

 !@pawtab_type

 nn=size(Pawtab)
 if (nn==0) return

 do ii=1,nn
   call pawtab_destroy(Pawtab(ii))
 end do

end subroutine pawtab_destroy_array
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_print
!! NAME
!! pawtab_print
!!
!! FUNCTION
!!  Print out the content of a pawtab datastructure
!!
!! INPUTS
!!  Pawtab<pawtab_type> Only for PAW, TABulated data initialized at start
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawtab_print(Pawtab,header,unit,prtvol,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtab_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
!arrays
 type(Pawtab_type) :: Pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: ityp,ntypat,my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(msg,'(6a)')&
&  ' ==================================== ',ch10,&
&  ' ==== Info on PAW TABulated data ==== ',ch10,&
&  ' ==================================== ',ch10
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ntypat=SIZE(Pawtab(:))

 do ityp=1,ntypat

  ! Print out integer values (dimensions)
  write(msg,'(a)')'                                 '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a)')'  ****************************** '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4,a)')'  **** Atom type ',ityp,' ****   '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a)')'  ****************************** '
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Number of (n,l) elements ....................... ',Pawtab(ityp)%basis_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Number of (l,m,n) elements ..................... ',Pawtab(ityp)%lmn_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Number of (i,j) elements (packed form) ......... ',Pawtab(ityp)%ij_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Max L+1 leading to non-zero Gaunt .............. ',Pawtab(ityp)%l_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Max L+1 leading to non-zero Gaunt (pawlcutd) ... ',Pawtab(ityp)%lcut_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  lmn2_size ...................................... ',Pawtab(ityp)%lmn2_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  lmnmix_sz ...................................... ',Pawtab(ityp)%lmnmix_sz
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Size of radial mesh ............................ ',Pawtab(ityp)%mesh_size
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  No of Q-points for tcorespl and tvalespl ....... ',Pawtab(ityp)%mqgrid
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  No of Q-points for the radial shape functions .. ',Pawtab(ityp)%mqgrid_shp
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Radial shape function type ..................... ',Pawtab(ityp)%shape_type
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  shape_lambda ................................... ',Pawtab(ityp)%shape_lambda
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Use pseudized core density ..................... ',Pawtab(ityp)%usetcore
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Use pseudized valence density .................. ',Pawtab(ityp)%usetvale
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Option for the use of hat density in XC terms .. ',Pawtab(ityp)%usexcnhat
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Use LDA+U ...................................... ',Pawtab(ityp)%usepawu
  call wrtout(ab_out,msg,'COLL')
  if (Pawtab(ityp)%usepawu/=0) then
    write(msg,'(a,i4)')'  L on which U is applied ........................ ',Pawtab(ityp)%lpawu
    call wrtout(ab_out,msg,'COLL')
  end if
  write(msg,'(a,i4)')'  Use Local Exact exchange ....................... ',Pawtab(ityp)%useexexch
  call wrtout(ab_out,msg,'COLL')
  if (Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,i4)')'  L on which local exact-exchange is applied ..... ',Pawtab(ityp)%lexexch
    call wrtout(ab_out,msg,'COLL')
  end if
  if (Pawtab(ityp)%usepawu/=0.or.Pawtab(ityp)%useexexch/=0) then
    write(msg,'(a,i4)')'  Number of (i,j) elements for PAW+U or EXX ..... ',Pawtab(ityp)%ij_proj
    call wrtout(ab_out,msg,'COLL')
    write(msg,'(a,i4)')'  Number of projectors on which U or EXX acts .... ',Pawtab(ityp)%nproju
    call wrtout(ab_out,msg,'COLL')
  end if

  ! "Has" flags
  write(msg,'(a,i4)')'  Has kij   ...................................... ',Pawtab(ityp)%has_kij
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Has nabla ...................................... ',Pawtab(ityp)%has_nabla
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Has vhntzc ..................................... ',Pawtab(ityp)%has_vhntzc
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,i4)')'  Has vhnzc ...................................... ',Pawtab(ityp)%has_vhnzc
  call wrtout(ab_out,msg,'COLL')
  !
  ! Real scalars
  write(msg,'(a,es16.8)')'  1/q d(tNcore(q))/dq for q=0 .....................',Pawtab(ityp)%dncdq0
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  1/q d(tNvale(q))/dq for q=0 .....................',Pawtab(ityp)%dnvdq0
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  XC energy for the core density ..................',Pawtab(ityp)%exccore
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  Mixing of exact exchange (PBE0) .................',Pawtab(ityp)%exchmix
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  Radius of the PAW sphere ........................',Pawtab(ityp)%rpaw
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,es16.8)')'  Compensation charge radius (if >rshp, g(r)=0) ...',Pawtab(ityp)%rshp !(if r>rshp, g(r)=zero)
  call wrtout(ab_out,msg,'COLL')
  if (Pawtab(ityp)%shape_type==2) then
   write(msg,'(a,es16.8)')'  Sigma parameter in gaussian shape function ......',Pawtab(ityp)%shape_sigma !(shape_type=2)
   call wrtout(ab_out,msg,'COLL')
  end if
  if (Pawtab(ityp)%usepawu/=0) then
   write(msg,'(a,es16.8)')'  Value of the U parameter [eV] ...................',Pawtab(ityp)%upawu*Ha_eV
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,es16.8)')'  Value of the J parameter [eV] ...................',Pawtab(ityp)%jpawu*Ha_eV
   call wrtout(ab_out,msg,'COLL')
  end if

 if (associated(Pawtab(ityp)%wvl%pngau)) then
   write(msg,'(a,es16.8)')'  WARNING: This Pawtab structure contains WVL data.'
   call wrtout(ab_out,msg,'COLL')
 end if

 end do ! ityp
 !
 ! The other (huge) arrays are not reported..

end subroutine pawtab_print
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/pawtab_bcast
!! NAME
!! pawtab_bcast
!!
!! FUNCTION
!! Communicate pawtab data to all processors
!!
!! INPUTS
!! comm_mpi= communicator used to broadcast data
!!
!! SIDE EFFECTS
!!  pawtab=<type pawtab_type>=a pawtab datastructure
!!
!! PARENTS
!!      pawbcast
!!
!! CHILDREN
!!
!! NOTES
!!      IMPORTANT: NOT ALL THER DATA aRE BROADCASTED !!!!
!!      Only data from pseudopotential file are broadcasted
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"

subroutine pawtab_bcast(pawtab,comm_mpi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtab_bcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm_mpi
 type(pawtab_type),intent(inout) :: pawtab

!Local variables-------------------------------
!scalars
 integer :: ierr,ii,indx,me,nn,nn_int,nn_scalars,isz1
 integer :: isz1_tvalespl,isz1_tproj
 integer :: isz2,isz2_tcoredens,isz2_rholoc_d,isz2_tvalespl
 character (len=500) :: message
!arrays
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)


!*************************************************************************

 me=xcomm_rank(comm_mpi)
 nn_scalars=7 !number of scalar variables passed in list_dpr: see below
 nn_int=0

!Calculate the size of the integers
 if(me==0) then
  if(pawtab%wvl%ptotgau>0) then
    if(associated(pawtab%wvl%pngau)) then
      nn_int=nn_int+size(pawtab%wvl%pngau)
    else
      message='wvl%pngau is not associated '
      call wrtout(std_out,message,'COLL')
      call leave_new('COLL')
    end if
  end if
  if(associated(pawtab%orbitals)) then
    nn_int=nn_int+size(pawtab%orbitals)
  end if
 end if

!Calculate the size of the reals
 if(me==0) then
   nn=nn_scalars 
   if (associated(pawtab%coredens)) then
     nn=nn+size(pawtab%coredens)
   end if
   isz2_tcoredens=0
   if (associated(pawtab%tcoredens)) then
     nn=nn+size(pawtab%tcoredens)
     isz2_tcoredens=size(pawtab%tcoredens,2)
   end if
   if (associated(pawtab%rhoij0)) then
     nn=nn+size(pawtab%rhoij0)
   end if
   if (associated(pawtab%dij0)) then
     nn=nn+size(pawtab%dij0)
   end if
   if (associated(pawtab%kij)) then
     nn=nn+size(pawtab%kij)
   end if
   if (associated(pawtab%vhntzc)) then
     nn=nn+size(pawtab%vhntzc)
   end if
   if (associated(pawtab%vhnzc)) then
     nn=nn+size(pawtab%vhnzc)
   end if
   if (associated(pawtab%shape_alpha)) then
     nn=nn+size(pawtab%shape_alpha)
   end if
   if (associated(pawtab%shape_q)) then
     nn=nn+size(pawtab%shape_q)
   end if
   if (associated(pawtab%shapefunc)) then
     nn=nn+size(pawtab%shapefunc)
   end if
   if (associated(pawtab%phi)) then
     nn=nn+size(pawtab%phi)
   end if
   if (associated(pawtab%tphi)) then
     nn=nn+size(pawtab%tphi)
   end if
   if (associated(pawtab%tcorespl)) then
     nn=nn+size(pawtab%tcorespl)
   end if
   isz1_tvalespl=0 ; isz2_tvalespl=0
   if (associated(pawtab%tvalespl)) then
     nn=nn+size(pawtab%tvalespl)
     isz1_tvalespl=size(pawtab%tvalespl,1)
     isz2_tvalespl=size(pawtab%tvalespl,2)
   end if
   isz2_rholoc_d=0
   if (pawtab%wvl%rholoc%msz>0) then
     if (associated(pawtab%wvl%rholoc%d)) then
       nn=nn+size(pawtab%wvl%rholoc%d)
       isz2_rholoc_d=size(pawtab%wvl%rholoc%d,2)
     end if
     if (associated(pawtab%wvl%rholoc%rad)) then
       nn=nn+size(pawtab%wvl%rholoc%rad)
     end if
   end if
   if(pawtab%wvl%ptotgau>0) then
     if (associated(pawtab%wvl%parg)) then
       nn=nn+size(pawtab%wvl%parg)
     end if
     if (associated(pawtab%wvl%pfac)) then
       nn=nn+size(pawtab%wvl%pfac)
     end if
   end if
   isz1_tproj=0
   if(associated(pawtab%tproj)) then
     nn=nn+size(pawtab%tproj)
     isz1_tproj=size(pawtab%tproj,1)
   end if
 end if

!Broadcast the integers
 ABI_ALLOCATE(list_int,(25))
 if(me==0) then
   list_int(1)=pawtab%basis_size
   list_int(2)=pawtab%lmn_size
   list_int(3)=pawtab%l_size
   list_int(4)=pawtab%lmn2_size
   list_int(5)=pawtab%shape_type
   list_int(6)=pawtab%shape_lambda
   list_int(7)=pawtab%has_vhnzc
   list_int(8)=pawtab%usetcore
   list_int(9)=pawtab%usexcnhat
   list_int(10)=pawtab%usetvale
   list_int(11)=pawtab%has_vhntzc
   list_int(12)=pawtab%mqgrid
   list_int(13)=pawtab%has_kij
   list_int(14)=pawtab%core_mesh_size
   list_int(15)=pawtab%tnvale_mesh_size
   list_int(16)=pawtab%wvl%ptotgau
   list_int(17)=pawtab%wvl%rholoc%msz
   list_int(18)=pawtab%wvl%npspcode_init_guess
   list_int(19)=isz2_tcoredens
   list_int(20)=isz1_tvalespl
   list_int(21)=isz2_tvalespl
   list_int(22)=isz2_rholoc_d
   list_int(23)=isz1_tproj
   list_int(24)=nn
   list_int(25)=nn_int
 end if

 call xcast_mpi(list_int,0,comm_mpi,ierr)

 if(me/=0) then
   pawtab%basis_size=list_int(1)
   pawtab%lmn_size=list_int(2)
   pawtab%l_size=list_int(3)
   pawtab%lmn2_size=list_int(4)
   pawtab%shape_type=list_int(5)
   pawtab%shape_lambda=list_int(6)
   pawtab%has_vhnzc=list_int(7)
   pawtab%usetcore=list_int(8)
   pawtab%usexcnhat=list_int(9)
   pawtab%usetvale=list_int(10)
   pawtab%has_vhntzc=list_int(11)
   pawtab%mqgrid=list_int(12)
   pawtab%has_kij=list_int(13)
   pawtab%core_mesh_size=list_int(14)
   pawtab%tnvale_mesh_size=list_int(15)
   pawtab%wvl%ptotgau=list_int(16)
   pawtab%wvl%rholoc%msz=list_int(17)
   pawtab%wvl%npspcode_init_guess=list_int(18)
   isz2_tcoredens=list_int(19)
   isz1_tvalespl=list_int(20)
   isz2_tvalespl=list_int(21)
   isz2_rholoc_d=list_int(22)
   isz1_tproj=list_int(23)
   nn=list_int(24)
   nn_int=list_int(25)
 end if
 ABI_DEALLOCATE(list_int)

!Broadcast integer arrays:
 if(nn_int>0) then
   ABI_ALLOCATE(list_int,(nn_int))
   if(me==0) then
     indx=1
     if(pawtab%wvl%ptotgau>0) then
       isz1=pawtab%basis_size
       if(.not. associated(pawtab%wvl%pngau)) then
         write(message, '(a,a,a)' )&
&         '  pawtab_bcast: wvl%pngau is not associated ',ch10,&
&         '  BUG: contact ABINIT group.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
       if(size(pawtab%wvl%pngau) /= isz1) then
         write(message, '(a,a,a)' )&
&         '  pawtab_bcast: wvl%pngau invalid size ',ch10,&
&         '  BUG: contact ABINIT group.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
       list_int(1:isz1)=pawtab%wvl%pngau(1:isz1)
       indx=indx+isz1
     end if
     isz1=pawtab%basis_size
     if(.not. associated(pawtab%orbitals)) then
       write(message, '(a,a,a)' )&
&       '  pawtab_bcast: pawtab%orbitals is not associated ',ch10,&
&       '  BUG: contact ABINIT group.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(size(pawtab%orbitals)/=isz1) then
       write(message, '(a,a,a)' )&
&       '  pawtab_bcast: pawtab%orbitals: wrong size ',ch10,&
&       '  BUG: contact ABINIT group.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_int(indx:indx+isz1-1)=pawtab%orbitals(:)
     indx=indx+isz1
   end if

   call xcast_mpi(list_int,0,comm_mpi,ierr)

   if(me/=0) then
     indx=1
     if(pawtab%wvl%ptotgau>0) then
       if(associated(pawtab%wvl%pngau)) then
         ABI_DEALLOCATE(pawtab%wvl%pngau)
       end if
       isz1=pawtab%basis_size
       ABI_ALLOCATE(pawtab%wvl%pngau,(isz1))
       pawtab%wvl%pngau(:)=list_int(indx:indx+isz1-1)
       indx=indx+isz1
     end if
     if(associated(pawtab%orbitals)) then
       ABI_DEALLOCATE(pawtab%orbitals)
     end if
     isz1=pawtab%basis_size
     ABI_ALLOCATE(pawtab%orbitals,(isz1))
     pawtab%orbitals(:)=list_int(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   ABI_DEALLOCATE(list_int)
 end if

!Broadcast the reals
 ABI_ALLOCATE(list_dpr,(nn))
 if(me==0) then
   list_dpr(1)=pawtab%shape_sigma
   list_dpr(2)=pawtab%rshp
   list_dpr(3)=pawtab%rpaw
   list_dpr(4)=pawtab%dncdq0
   list_dpr(5)=pawtab%dnvdq0
   list_dpr(6)=pawtab%exccore
   list_dpr(7)=pawtab%rcore
   indx=nn_scalars+1
   if (associated(pawtab%coredens)) then
     isz1=size(pawtab%coredens)
     if(isz1/=pawtab%core_mesh_size) then
       message='coredens: sz1 /= pawtab%core_mesh_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=pawtab%coredens(:)
     indx=indx+isz1
   end if
   if (associated(pawtab%tcoredens)) then
     isz1=size(pawtab%tcoredens,1)
     if(isz1/=pawtab%core_mesh_size) then
       message='tcoredens: sz1 /= pawtab%core_mesh_size '
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(isz2_tcoredens /=1 .and. isz2_tcoredens /=6) then
       message='tcoredens: sz2 /= 1 and 6 '
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     do ii=1,isz2_tcoredens
       list_dpr(indx:indx+isz1-1)=pawtab%tcoredens(1:isz1,ii)
       indx=indx+isz1
     end do
   end if
   if (associated(pawtab%rhoij0)) then
     isz1=size(pawtab%rhoij0)
     if(isz1/=pawtab%lmn2_size) then
       message='rhoij0: sz1 /= pawtab%lmn2_size '
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=pawtab%rhoij0(:)
     indx=indx+isz1
   end if
   if (associated(pawtab%dij0)) then
     isz1=size(pawtab%dij0)
     if(isz1/=pawtab%lmn2_size) then
       message='dij0: sz1 /= pawtab%lmn2_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=pawtab%dij0(:)
     indx=indx+isz1
   end if
   if (associated(pawtab%kij)) then
     isz1=size(pawtab%kij)
     if(isz1/=pawtab%lmn2_size) then
       message='kij: sz1 /= pawtab%lmn2_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=pawtab%kij(:)
     indx=indx+isz1
   end if
   if (associated(pawtab%vhntzc)) then
     isz1=size(pawtab%vhntzc)
     if(isz1/=pawtab%mesh_size) then
       message='vhntzc: sz1 /= pawtab%mesh_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=pawtab%vhntzc(:)
     indx=indx+isz1
   end if
   if (associated(pawtab%vhnzc)) then
     isz1=size(pawtab%vhnzc)
     if(isz1/=pawtab%mesh_size) then
       message='vhnzc: sz1 /= pawtab%mesh_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=pawtab%vhnzc(:)
     indx=indx+isz1
   end if
   if (associated(pawtab%shape_alpha)) then
     isz1=size(pawtab%shape_alpha)
     if(isz1/=2*pawtab%l_size) then
       message='shape_alpha: sz1 /= 2*pawtab%l_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=reshape(pawtab%shape_alpha(:,:),(/isz1/))
     indx=indx+isz1
   end if
   if (associated(pawtab%shape_q)) then
     isz1=size(pawtab%shape_q)
     if(isz1/=2*pawtab%l_size) then
       message='shape_q: sz1 /= 2*pawtab%l_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=reshape(pawtab%shape_q(:,:),(/isz1/))
     indx=indx+isz1
   end if
   if (associated(pawtab%shapefunc)) then
     isz1=size(pawtab%shapefunc)
     if(isz1/=pawtab%mesh_size*pawtab%l_size) then
       message='shapefunc: sz1 /= pawtab%mesh_size*pawtab%l_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=reshape(pawtab%shapefunc(:,:),(/isz1/))
     indx=indx+isz1
   end if
   if (associated(pawtab%phi)) then
     isz1=size(pawtab%phi)
     if(isz1/=pawtab%mesh_size*pawtab%basis_size) then
       message='phi: sz1 /= pawtab%mesh_size*pawtab%basis_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=reshape(pawtab%phi(:,:),(/isz1/))
     indx=indx+isz1
   end if
   if (associated(pawtab%tphi)) then
     isz1=size(pawtab%tphi)
     if(isz1/=pawtab%mesh_size*pawtab%basis_size) then
       message='tphi: sz1 /= pawtab%mesh_size*pawtab%basis_size'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=reshape(pawtab%tphi(:,:),(/isz1/))
     indx=indx+isz1
   end if
   if (associated(pawtab%tcorespl)) then
     isz1=size(pawtab%tcorespl)
     if(isz1/=2*pawtab%mqgrid) then
       message='tcorespl: sz1 /= 2*pawtab%mqgrid'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     list_dpr(indx:indx+isz1-1)=reshape(pawtab%tcorespl(:,:),(/isz1/))
     indx=indx+isz1
   end if
!  tvalespl:
   if (associated(pawtab%tvalespl)) then
     if(isz2_tvalespl /=1 .and. isz2_tvalespl /=2) then
       message='tvalespl: sz2 /= 1 and 2'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     do ii=1,isz2_tvalespl
       list_dpr(indx:indx+isz1_tvalespl-1)=pawtab%tvalespl(1:isz1_tvalespl,ii)
       indx=indx+isz1_tvalespl
     end do
   end if
!  rholoc:
   if(pawtab%wvl%rholoc%msz<0) then
     message='pawtab%wvl%rholoc%msz < 0'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(pawtab%wvl%rholoc%msz>0) then
     if(.not. associated(pawtab%wvl%rholoc%d) .or. &
&       .not. associated(pawtab%wvl%rholoc%rad)) then
       message='rholoc%msz>0 and pawtab%wvl%rholoc not associated'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if( size(pawtab%wvl%rholoc%rad,1) .ne. pawtab%wvl%rholoc%msz .or.&
&        size(pawtab%wvl%rholoc%d,1) .ne. pawtab%wvl%rholoc%msz ) then
       message='wrong size in pawtab%wvl%rholoc'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     isz1=pawtab%wvl%rholoc%msz
     do ii=1,isz2_rholoc_d
       list_dpr(indx:indx+isz1-1)=pawtab%wvl%rholoc%d(1:isz1,ii)
       indx=indx+isz1
     end do
     list_dpr(indx:indx+isz1-1)=pawtab%wvl%rholoc%rad(1:isz1)
     indx=indx+isz1
   end if
   if(pawtab%wvl%ptotgau>0) then
     if(.not. associated(pawtab%wvl%parg) .or. &
&       .not. associated(pawtab%wvl%pfac)) then
       message='wvl%ptotgau>0 and related arrays not associated'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(size(pawtab%wvl%parg,1) .ne. 2 .or. &
&       size(pawtab%wvl%parg,2) .ne. pawtab%wvl%ptotgau .or. &
&       size(pawtab%wvl%pfac,1) .ne. 2 .or. &
&       size(pawtab%wvl%pfac,2) .ne. pawtab%wvl%ptotgau) then
       message='wrong size in pawtab%wvl objects'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     isz2=pawtab%wvl%ptotgau
     do ii=1,isz2
       list_dpr(indx:indx+1)=pawtab%wvl%parg(1:2,ii)
       indx=indx+2
     end do
     do ii=1,isz2
       list_dpr(indx:indx+1)=pawtab%wvl%pfac(1:2,ii)
       indx=indx+2
     end do
   end if
!  tproj
   if(isz1_tproj>0) then
     isz2=pawtab%basis_size; isz1=isz1_tproj
     do ii=1,isz2
       list_dpr(indx:indx+isz1-1)=pawtab%tproj(:,ii)
       indx=indx+isz1
     end do
   end if
 end if

 call xcast_mpi(list_dpr,0,comm_mpi,ierr)

 if(me/=0) then
   pawtab%shape_sigma=list_dpr(1)
   pawtab%rshp       =list_dpr(2)
   pawtab%rpaw       =list_dpr(3)
   pawtab%dncdq0     =list_dpr(4)
   pawtab%dnvdq0     =list_dpr(5)
   pawtab%exccore    =list_dpr(6)
   pawtab%rcore      =list_dpr(7)
   indx=nn_scalars+1
   if (associated(pawtab%coredens)) then
     ABI_DEALLOCATE(pawtab%coredens)
     isz1=pawtab%core_mesh_size
     ABI_ALLOCATE(pawtab%coredens,(isz1))
     pawtab%coredens(:)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
!  tcoredens:
   if (associated(pawtab%tcoredens)) then
     ABI_DEALLOCATE(pawtab%tcoredens)
   end if
   if(isz2_tcoredens>0) then
     isz1=pawtab%core_mesh_size
     ABI_ALLOCATE(pawtab%tcoredens,(isz1,isz2_tcoredens))
     do ii=1,isz2_tcoredens
       pawtab%tcoredens(:,ii)=list_dpr(indx:indx+isz1-1)
       indx=indx+isz1
     end do
   end if
   if (associated(pawtab%rhoij0)) then
     ABI_DEALLOCATE(pawtab%rhoij0)
     isz1=pawtab%lmn2_size
     ABI_ALLOCATE(pawtab%rhoij0,(isz1))
     pawtab%rhoij0(:)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if (associated(pawtab%dij0)) then
     ABI_DEALLOCATE(pawtab%dij0)
     isz1=pawtab%lmn2_size
     ABI_ALLOCATE(pawtab%dij0,(isz1))
     pawtab%dij0(:)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if (associated(pawtab%kij)) then
     ABI_DEALLOCATE(pawtab%kij)
     isz1=pawtab%lmn2_size
     ABI_ALLOCATE(pawtab%kij,(isz1))
     pawtab%kij(:)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if (associated(pawtab%vhntzc)) then
     ABI_DEALLOCATE(pawtab%vhntzc)
     isz1=pawtab%mesh_size
     ABI_ALLOCATE(pawtab%vhntzc,(isz1))
     pawtab%vhntzc(:)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if (associated(pawtab%vhnzc)) then
     ABI_DEALLOCATE(pawtab%vhnzc)
     isz1=pawtab%mesh_size
     ABI_ALLOCATE(pawtab%vhnzc,(isz1))
     pawtab%vhnzc(:)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if (associated(pawtab%shape_alpha)) then
     ABI_DEALLOCATE(pawtab%shape_alpha)
     isz1=pawtab%l_size
     ABI_ALLOCATE(pawtab%shape_alpha,(2,isz1))
     pawtab%shape_alpha(:,:)=reshape(list_dpr(indx:indx+2*isz1-1),(/2,isz1/))
     indx=indx+2*isz1
   end if
   if (associated(pawtab%shape_q)) then
     ABI_DEALLOCATE(pawtab%shape_q)
     isz1=pawtab%l_size
     ABI_ALLOCATE(pawtab%shape_q,(2,isz1))
     pawtab%shape_q(:,:)=reshape(list_dpr(indx:indx+2*isz1-1),(/2,isz1/))
     indx=indx+2*isz1
   end if
   if (associated(pawtab%shapefunc)) then
     ABI_DEALLOCATE(pawtab%shapefunc)
     isz1=pawtab%mesh_size
     isz2=pawtab%l_size
     ABI_ALLOCATE(pawtab%shapefunc,(isz1,isz2))
     pawtab%shapefunc(:,:)=reshape(list_dpr(indx:indx+isz1*isz2-1),(/isz1,isz2/))
     indx=indx+isz1*isz2
   end if
   if (associated(pawtab%phi)) then
     ABI_DEALLOCATE(pawtab%phi)
     isz1=pawtab%mesh_size
     isz2=pawtab%basis_size
     ABI_ALLOCATE(pawtab%phi,(isz1,isz2))
     pawtab%phi(:,:)=reshape(list_dpr(indx:indx+isz1*isz2-1),(/isz1,isz2/))
     indx=indx+isz1*isz2
   end if
   if (associated(pawtab%tphi)) then
     ABI_DEALLOCATE(pawtab%tphi)
     isz1=pawtab%mesh_size
     isz2=pawtab%basis_size
     ABI_ALLOCATE(pawtab%tphi,(isz1,isz2))
     pawtab%tphi(:,:)=reshape(list_dpr(indx:indx+isz1*isz2-1),(/isz1,isz2/))
     indx=indx+isz1*isz2
   end if
   if (associated(pawtab%tcorespl)) then
     ABI_DEALLOCATE(pawtab%tcorespl)
     isz1=pawtab%mqgrid
     ABI_ALLOCATE(pawtab%tcorespl,(isz1,2))
     pawtab%tcorespl(:,:)=reshape(list_dpr(indx:indx+2*isz1-1),(/isz1,2/))
     indx=indx+2*isz1
   end if
!  tvalespl:
   if (associated(pawtab%tvalespl)) then
     ABI_DEALLOCATE(pawtab%tvalespl)
   end if
   if(isz2_tvalespl>0) then
     ABI_ALLOCATE(pawtab%tvalespl,(isz1_tvalespl,isz2_tvalespl))
     do ii=1,isz2_tvalespl
       pawtab%tvalespl(1:isz1_tvalespl,ii)=list_dpr(indx:indx+isz1_tvalespl-1)
       indx=indx+isz1_tvalespl
     end do
   end if
!  rholoc:
   if(pawtab%wvl%rholoc%msz>0) then
     if(associated(pawtab%wvl%rholoc%d)) then
       ABI_DEALLOCATE(pawtab%wvl%rholoc%d)
     end if
     if(associated(pawtab%wvl%rholoc%rad)) then
       ABI_DEALLOCATE(pawtab%wvl%rholoc%rad)
     end if
     isz1=pawtab%wvl%rholoc%msz
     ABI_ALLOCATE(pawtab%wvl%rholoc%d,(isz1,isz2_rholoc_d))
     ABI_ALLOCATE(pawtab%wvl%rholoc%rad,(isz1))
     do ii=1,isz2_rholoc_d
       pawtab%wvl%rholoc%d(1:isz1,ii)=&
&       list_dpr(indx:indx+isz1-1)
       indx=indx+isz1
     end do
     pawtab%wvl%rholoc%rad(1:isz1)=&
&     list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if(pawtab%wvl%ptotgau>0) then
     if(associated(pawtab%wvl%parg)) then
       ABI_DEALLOCATE(pawtab%wvl%parg)
     end if
     if(associated(pawtab%wvl%pfac)) then
       ABI_DEALLOCATE(pawtab%wvl%pfac)
     end if
     isz2=pawtab%wvl%ptotgau
     ABI_ALLOCATE(pawtab%wvl%parg,(2,isz2))
     ABI_ALLOCATE(pawtab%wvl%pfac,(2,isz2))
     do ii=1,isz2
       pawtab%wvl%parg(1:2,ii)=list_dpr(indx:indx+1)
       indx=indx+2
     end do
     do ii=1,isz2
       pawtab%wvl%pfac(1:2,ii)=list_dpr(indx:indx+1)
       indx=indx+2
     end do
   end if
!  tproj
   if(isz1_tproj>0) then
     if(associated(pawtab%tproj)) then
       ABI_DEALLOCATE(pawtab%tproj)
     end if
     isz1=isz1_tproj; isz2=pawtab%basis_size
     ABI_ALLOCATE(pawtab%tproj,(isz1,isz2))
     do ii=1,isz2
       pawtab%tproj(:,ii)=list_dpr(indx:indx+isz1-1)
       indx=indx+isz1
     end do
   end if
 end if
 ABI_DEALLOCATE(list_dpr)

end subroutine pawtab_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_destroy
!! NAME
!!  wvlpaw_destroy
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a wvlpaw structure
!!
!! SIDE EFFECTS
!!  wvlpaw<type(wvlpaw_type)>=datastructure to be nullified.
!!  All associated pointers are deallocated.
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_destroy(wvlpaw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvlpaw_destroy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wvlpaw_type),intent(inout) :: wvlpaw

! *************************************************************************

 !@wvlpaw_type

 if(associated(wvlpaw%pngau)) then
   ABI_DEALLOCATE(wvlpaw%pngau)
 end if
 if(associated(wvlpaw%parg)) then
   ABI_DEALLOCATE(wvlpaw%parg)
 end if
 if(associated(wvlpaw%pfac)) then
   ABI_DEALLOCATE(wvlpaw%pfac)
 end if
 if(associated(wvlpaw%rholoc%d)) then
   ABI_DEALLOCATE(wvlpaw%rholoc%d)
 end if
 if(associated(wvlpaw%rholoc%rad)) then
   ABI_DEALLOCATE(wvlpaw%rholoc%rad)
 end if

 wvlpaw%ptotgau=0
 wvlpaw%rholoc%msz=0

end subroutine wvlpaw_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_pawtab/wvlpaw_nullify
!! NAME
!!  wvlpaw_nullify
!!
!! FUNCTION
!!  Nullify pointers and flags in a wvlpaw structure
!!
!! SIDE EFFECTS
!!  wvlpaw=datastructure to be nullified
!!
!! PARENTS
!!      m_pawtab
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvlpaw_nullify(wvlpaw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvlpaw_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wvlpaw_type),intent(inout) :: wvlpaw

! *************************************************************************

 !@wvlpaw_type

 nullify(wvlpaw%pngau)
 nullify(wvlpaw%parg)
 nullify(wvlpaw%pfac)
 nullify(wvlpaw%rholoc%d)
 nullify(wvlpaw%rholoc%rad)

 wvlpaw%npspcode_init_guess=0
 wvlpaw%ptotgau=0
 wvlpaw%rholoc%msz=0

end subroutine wvlpaw_nullify
!!***

!----------------------------------------------------------------------

END MODULE m_pawtab
!!***
