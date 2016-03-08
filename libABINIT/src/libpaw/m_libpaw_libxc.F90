!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_libpaw_libxc_funcs
!! NAME
!!  m_libpaw_libxc_funcs
!!
!! FUNCTION
!!  Module containing interfaces to the LibXC library, for exchange
!!  correlation potentials and energies.
!!  A global data-structure (pawxcfuncs) is used to store the XC parameters.
!!  This structure has to be initialized with a negative value of ixc.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MOliveira, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  This file comes directly from m_libxc_functionals.F90 module delivered with ABINIT.
!!
!!  How to initialize (finalize) libXC datastructure:
!!   call libpaw_libxc_init(abinit_ixc,nspin)
!!   call libpaw_libxc_end()
!!
!! SOURCE

#include "libpaw.h"

module m_libpaw_libxc_funcs

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

#if defined LIBPAW_HAVE_LIBXC
 use xc_f90_types_m
 use libxc_funcs_m
 use xc_f90_lib_m
#endif

 implicit none

#if defined LIBPAW_HAVE_LIBXC
 type libpaw_libxc_functional
   private
   integer                :: family  ! LDA, GGA, etc.
   integer                :: id      ! identifier
   integer                :: nspin   ! # of spin components
   logical                :: has_fxc ! TRUE is fxc is available for the functional
   type(xc_f90_pointer_t) :: conf    ! the pointer used to call the library
   type(xc_f90_pointer_t) :: info    ! information about the functional
 end type libpaw_libxc_functional

 type(libpaw_libxc_functional) :: pawxcfuncs(2)

 private

 integer,save :: LIBPAW_IXC = HUGE(0)

 public :: libpaw_libxc_init        ! Initialize the desired XC functional, from LibXC.
 public :: libpaw_libxc_fullname    ! Return full name of the XC functional
 public :: libpaw_libxc_getid       ! Return identifer of a XC functional from its name
 public :: libpaw_libxc_global_ixc  ! The value of ixc used to initialize the global structure pawxcfuncs
 public :: libpaw_libxc_getvxc      ! Return XC potential and energy, from input density (event gradient etc...)
 public :: libpaw_libxc_isgga       ! Flag: is it a XC GGA functional?
 public :: libpaw_libxc_ismgga      ! Flag: is it a XC meta-GGA functional?
 public :: libpaw_libxc_has_kxc     ! Flag: does the functional provide Kxc?
 public :: libpaw_libxc_nspin       ! The number of spin components for the XC functionals
 public :: libpaw_libxc_end         ! End usage of LibXC functional.

contains
!!***

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_init
!! NAME
!!  libpaw_libxc_init
!!
!! FUNCTION
!!  Initialize the desired XC functional, from LibXC.
!!  * Call the LibXC initializer
!!  * Fill preliminary fields in module structures.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE

 subroutine libpaw_libxc_init(ixc,nspden)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_init'
 !use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nspden
 integer, intent(in) :: ixc

!Local variables-------------------------------
!scalars
 integer :: ii,jj,nspden_eff
 character(len=500) :: msg
 type(xc_f90_pointer_t) :: str

! *************************************************************************

!Save abinit value for reference 
 LIBPAW_IXC = ixc

 pawxcfuncs(1)%id = -ixc/1000
 pawxcfuncs(2)%id = -ixc - pawxcfuncs(1)%id*1000

 nspden_eff=min(nspden,2)
 pawxcfuncs(1)%nspin = nspden_eff
 pawxcfuncs(2)%nspin = nspden_eff

 do ii = 1, 2
   if (pawxcfuncs(ii)%id == 0) then
     pawxcfuncs(ii)%family = 0
     pawxcfuncs(ii)%has_fxc= .false.
     cycle
   end if

   ! Get XC functional family
   pawxcfuncs(ii)%family = xc_f90_family_from_id(pawxcfuncs(ii)%id)
   select case (pawxcfuncs(ii)%family)
   case (XC_FAMILY_LDA,XC_FAMILY_GGA,XC_FAMILY_HYB_GGA,XC_FAMILY_MGGA)
     call xc_f90_func_init(pawxcfuncs(ii)%conf,pawxcfuncs(ii)%info,pawxcfuncs(ii)%id,nspden_eff)
     pawxcfuncs(ii)%has_fxc=(iand(xc_f90_info_flags(pawxcfuncs(ii)%info),XC_FLAGS_HAVE_FXC)>0)
   case default
     write(msg, '(a,i8,2a,i8,4a)' )&
&      'Invalid IXC = ',ixc,ch10,&
&      'The LibXC functional family ',pawxcfuncs(ii)%family,&
&      'is currently unsupported!',ch10,&
&      'Please consult the LibXC documentation.',ch10
     MSG_ERROR(msg)
   end select

   if (pawxcfuncs(ii)%id == XC_LDA_C_XALPHA) then
     call xc_f90_lda_c_xalpha_set_par(pawxcfuncs(ii)%conf,zero)
   end if

   !Dump functional information
   call xc_f90_info_name(pawxcfuncs(ii)%info,msg)
   call wrtout(std_out,msg,'COLL')
   jj = 0
   call xc_f90_info_refs(pawxcfuncs(ii)%info,jj,str,msg)
   do while (jj >= 0)
     call wrtout(std_out,msg,'COLL')
     call xc_f90_info_refs(pawxcfuncs(ii)%info,jj,str,msg)
   end do
 end do

end subroutine libpaw_libxc_init
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_end
!! NAME
!!  libpaw_libxc_end
!!
!! FUNCTION
!!  End usage of LibXC functional. Call LibXC end function,
!!  and deallocate module contents.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE
 subroutine libpaw_libxc_end()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_end'
!End of the abilint section

 implicit none

!Local variables-------------------------------
 integer :: ii

! *************************************************************************

 LIBPAW_IXC = HUGE(0)

 do ii = 1, 2
   if (pawxcfuncs(ii)%id == 0) cycle
   call xc_f90_func_end(pawxcfuncs(ii)%conf)
 end do

 end subroutine libpaw_libxc_end
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_fullname
!! NAME
!!  libpaw_libxc_fullname
!!
!! FUNCTION
!!  Return full name of the XC functional
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

 function libpaw_libxc_fullname()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_fullname'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=100) :: libpaw_libxc_fullname

!Local variables-------------------------------
 character(len=100) :: xcname

! *************************************************************************

 if (pawxcfuncs(1)%id == 0) then
   if (pawxcfuncs(2)%id /= 0) then
     call xc_f90_info_name(pawxcfuncs(2)%info,libpaw_libxc_fullname)
   else
     libpaw_libxc_fullname='No XC functional'
   end if
 else if (pawxcfuncs(2)%id == 0) then
   if (pawxcfuncs(1)%id /= 0) then
     call xc_f90_info_name(pawxcfuncs(1)%info,libpaw_libxc_fullname)
   else
     libpaw_libxc_fullname='No XC functional'
   end if
 else
   call xc_f90_info_name(pawxcfuncs(1)%info,libpaw_libxc_fullname)
   call xc_f90_info_name(pawxcfuncs(2)%info,xcname)
   libpaw_libxc_fullname=trim(libpaw_libxc_fullname)//'+'//trim(xcname)
 end if

end function libpaw_libxc_fullname
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_getid
!! NAME
!!  libpaw_libxc_getid
!!
!! FUNCTION
!!  Return identifier of a XC functional from its name
!!  Return -1 if undefined 
!!
!! INPUTS
!!  xcname= string containing the name of a XC functional
!!
!! OUTPUT
!!
!! SOURCE
 function libpaw_libxc_getid(xcname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_getid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: libpaw_libxc_getid
 character(len=*),intent(in) :: xcname

!Local variables-------------------------------
 integer,external :: xc_f90_functional_get_number
 character(len=256) :: str

! *************************************************************************

 str=xcname
 if (xcname(1:3)=="XC_".or.xcname(1:3)=="xc_") str=xcname(4:)

 libpaw_libxc_getid=xc_f90_functional_get_number(trim(str))

end function libpaw_libxc_getid
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_global_ixc
!! NAME
!!  libpaw_libxc_global_ixc
!!
!! FUNCTION
!!  Return the value of ixc used to initialize the global structure pawxcfuncs
!!  Return HUGE(0) if undefined 
!!
!! OUTPUT
!!
!! SOURCE

 function libpaw_libxc_global_ixc()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_global_ixc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: libpaw_libxc_global_ixc

! *************************************************************************

 libpaw_libxc_global_ixc = LIBPAW_IXC

end function libpaw_libxc_global_ixc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_isgga
!! NAME
!!  libpaw_libxc_isgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

 function libpaw_libxc_isgga()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_isgga'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: libpaw_libxc_isgga

! *************************************************************************

 if (any(pawxcfuncs%family == XC_FAMILY_GGA) .or. any(pawxcfuncs%family == XC_FAMILY_HYB_GGA)) then
   libpaw_libxc_isgga = .true.
 else
   libpaw_libxc_isgga = .false.
 end if

end function libpaw_libxc_isgga
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_ismgga
!! NAME
!!  libpaw_libxc_ismgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a Meta-GGA or not
!!
!! SOURCE

function libpaw_libxc_ismgga()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_ismgga'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: libpaw_libxc_ismgga

! *************************************************************************

 if (any(pawxcfuncs%family == XC_FAMILY_MGGA)) then
   libpaw_libxc_ismgga = .true.
 else
   libpaw_libxc_ismgga = .false.
 end if

end function libpaw_libxc_ismgga
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_has_kxc
!! NAME
!!  libpaw_libxc_has_kxc
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  provides Kxc or not fxc in the libXC convention)
!!
!! SOURCE

function libpaw_libxc_has_kxc()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_has_kxc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: libpaw_libxc_has_kxc

!Local variables-------------------------------
 integer :: ii

! *************************************************************************

 libpaw_libxc_has_kxc = .true.

 do ii=1,2
   if (pawxcfuncs(ii)%id/=0) then
     if (.not.pawxcfuncs(ii)%has_fxc) libpaw_libxc_has_kxc = .false.
   end if
 end do

end function libpaw_libxc_has_kxc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_nspin
!! NAME
!!  libpaw_libxc_nspin
!!
!! FUNCTION
!!  Returns the number of spin components for the XC functionals
!!
!! SOURCE

function libpaw_libxc_nspin()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_nspin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: libpaw_libxc_nspin

! *************************************************************************

 if (any(pawxcfuncs%nspin == XC_POLARIZED)) then
   libpaw_libxc_nspin = 2
 else
   libpaw_libxc_nspin = 1
 end if

end function libpaw_libxc_nspin
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_getvxc
!! NAME
!!  libpaw_libxc_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (event gradient etc...)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE

 subroutine libpaw_libxc_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxc,&
&                 grho2,vxcgr,lrho,vxclrho,tau,vxctau,dvxc,d2vxc,xc_tb09_c) ! Optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libpaw_libxc_getvxc'
 !use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ndvxc,nd2vxc,npts,nspden,order
 real(dp),intent(in)  :: rho(npts,nspden)
 real(dp),intent(out) :: vxc(npts,nspden), exc(npts)
 real(dp),intent(in),optional :: grho2(npts,2*min(nspden,2)-1)
 real(dp),intent(out),optional :: vxcgr(npts,3)
 real(dp),intent(in),optional :: lrho(npts,nspden)
 real(dp),intent(out),optional :: vxclrho(npts,nspden)
 real(dp),intent(in),optional :: tau(npts,nspden)
 real(dp),intent(out),optional :: vxctau(npts,nspden)
 real(dp),intent(out),optional :: dvxc(npts,ndvxc)
 real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)
 real(dp),intent(in),optional :: xc_tb09_c

!Local variables -------------------------------
 integer  :: i, ipts
 real(dp) :: c, rhotmp(nspden), exctmp, sigma(3), vsigma(3), vxctmp(nspden)
 real(dp) :: v2rho2(3),v2rhosigma(6),v2sigma2(6),v3rho3(4)
 real(dp) :: lrhotmp(nspden), tautmp(nspden), vxclrhotmp(nspden), vxctautmp(nspden)
 real(dp), allocatable :: gnon(:)
 character(len=500) :: msg

! *************************************************************************

 ! Inititalize all relevant arrays to zero
 vxc=zero
 exc=zero
 vxctmp=zero
 exctmp=zero
 v2rho2=zero
 v2rhosigma=zero
 v2sigma2=zero
 v3rho3=zero
 if (order**2 >1) dvxc=zero
 if (order**2 >4) d2vxc=zero
 if (any(pawxcfuncs%family == XC_FAMILY_GGA) .or. any(pawxcfuncs%family == XC_FAMILY_HYB_GGA)) vxcgr=zero
 if (any(pawxcfuncs%family == XC_FAMILY_MGGA)) then
   vxcgr=zero
   vxclrho=zero
   vxctau=zero
 end if

 !The TB09 MGGA functional requires an extra quantity
 if (any(pawxcfuncs%id == XC_MGGA_X_TB09)) then
   if (PRESENT(xc_tb09_c)) then
     c = xc_tb09_c
   else
     LIBPAW_ALLOCATE(gnon,(npts))
     do ipts = 1, npts
       if (sum(rho(ipts, :)) <= 1e-7_dp) then
         gnon(ipts) = zero
       else
         if (nspden == 1) then
           gnon(ipts) = sqrt(grho2(ipts,1))/rho(ipts, 1)
         else
           gnon(ipts) = sqrt(grho2(ipts,3))/sum(rho(ipts, :))
         end if
       end if
     end do
     c = -0.012_dp + 1.023_dp*sqrt(sum(gnon)/npts)
     LIBPAW_DEALLOCATE(gnon)
   end if
   do i = 1, 2
     if (pawxcfuncs(i)%id == XC_MGGA_X_TB09) then
       call xc_f90_mgga_x_tb09_set_par(pawxcfuncs(i)%conf, c)
       if (PRESENT(xc_tb09_c)) then
         write(msg, '(2a,f9.6)' ) ch10,&
              &     ' In the functional TB09 c is fixed by the user (xc_tb09_c set in input file) and is equal to ', c
       else
         write(msg, '(2a,f9.6)' ) ch10,&
              &     ' In the functional TB09 c = ', c
       end if
       call wrtout(std_out,msg,'COLL')
     end if
   end do
 end if

 !Loop over points
 do ipts = 1, npts

   ! Convert the quantities provided by ABINIT to the ones needed by libxc
   if (nspden == 1) then
     ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc
     ! expects the total density
     rhotmp(1:nspden) = two*rho(ipts,1:nspden)
   else
     rhotmp(1:nspden) = rho(ipts,1:nspden)
   end if
   if (any(pawxcfuncs%family == XC_FAMILY_GGA) .or. any(pawxcfuncs%family == XC_FAMILY_HYB_GGA)&
&      .or. any(pawxcfuncs%family == XC_FAMILY_MGGA)) then
     sigma=zero
     if (nspden==1) then
       ! ABINIT passes |grho_up|^2 while Libxc needs |grho_tot|^2
       sigma(1) = four*grho2(ipts,1)
     else
       ! ABINIT passes |grho_up|^2, |grho_dn|^2, and |grho_tot|^2
       ! while Libxc needs |grho_up|^2, grho_up.grho_dn, and |grho_dn|^2
       sigma(1) = grho2(ipts,1)
       sigma(2) = (grho2(ipts,3) - grho2(ipts,1) - grho2(ipts,2))/two
       sigma(3) = grho2(ipts,2)
     end if
   end if
   if (any(pawxcfuncs%family == XC_FAMILY_MGGA)) then
     if (nspden==1) then
       lrhotmp(1:nspden) = two*lrho(ipts,1:nspden)
       tautmp(1:nspden) = two*tau(ipts,1:nspden)
     else
       lrhotmp(1:nspden) = lrho(ipts,1:nspden)
       tautmp(1:nspden) = tau(ipts,1:nspden)
     end if
   end if

   !Loop over functionals
   do i = 1,2
     if (pawxcfuncs(i)%id == 0) cycle

     !Get the potential (and possibly the energy)
     if (iand(xc_f90_info_flags(pawxcfuncs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
       select case (pawxcfuncs(i)%family)
       case (XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(pawxcfuncs(i)%conf,1,rhotmp(1),exctmp,vxctmp(1))
         if (order**2 > 1) then
           call xc_f90_lda_fxc(pawxcfuncs(i)%conf,1,rhotmp(1),v2rho2(1))
         endif
         if (order**2 > 4) then
           call xc_f90_lda_kxc(pawxcfuncs(i)%conf,1,rhotmp(1),v3rho3(1))
         endif
       case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         call xc_f90_gga_exc_vxc(pawxcfuncs(i)%conf,1,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
         if (order**2 > 1) then
           call xc_f90_gga_fxc(pawxcfuncs(i)%conf,1,rhotmp(1),sigma(1),v2rho2(1),v2rhosigma(1),v2sigma2(1))
         endif

       case (XC_FAMILY_MGGA)
         call xc_f90_mgga_exc_vxc(pawxcfuncs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                     tautmp(1),exctmp,vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
       end select

     else
       exctmp=zero
       select case (pawxcfuncs(i)%family)
       case (XC_FAMILY_LDA)
         call xc_f90_lda_vxc(pawxcfuncs(i)%conf,1,rhotmp(1),vxctmp(1))
       case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         call xc_f90_gga_vxc(pawxcfuncs(i)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))


       case (XC_FAMILY_MGGA)
         call xc_f90_mgga_vxc(pawxcfuncs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                     tautmp(1),vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
       end select
     end if

     exc(ipts) = exc(ipts) + exctmp
     vxc(ipts,1:nspden) = vxc(ipts,1:nspden) + vxctmp(1:nspden)

!    Deal with fxc and kxc
     if (order**2>1) then
       select case (pawxcfuncs(i)%family)
       case (XC_FAMILY_LDA)
         if (nspden==1) then
           if(order>=2) then
             dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
             if(order==3) then
               d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
             endif
           else
             dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
             dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(1)
           endif
         else
           dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
           dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(2)
           dvxc(ipts,3)=dvxc(ipts,3)+v2rho2(3)
           if(order==3) then
             d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
             d2vxc(ipts,2)=d2vxc(ipts,2)+v3rho3(2)
             d2vxc(ipts,3)=d2vxc(ipts,3)+v3rho3(3)
             d2vxc(ipts,4)=d2vxc(ipts,4)+v3rho3(4)
           endif
         endif
       case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         select case(xc_f90_info_kind(pawxcfuncs(i)%info))
         case(XC_EXCHANGE)
           if (nspden==1) then
             dvxc(ipts,1)=v2rho2(1)*two
             dvxc(ipts,2)=dvxc(ipts,1)
             dvxc(ipts,3)=two*two*vsigma(1)
             dvxc(ipts,4)=dvxc(ipts,3)
             dvxc(ipts,5)=four*two*v2rhosigma(1)
             dvxc(ipts,6)=dvxc(ipts,5)
             dvxc(ipts,7)=two*four*four*v2sigma2(1)
             dvxc(ipts,8)=dvxc(ipts,7)
           else
             dvxc(ipts,1)=v2rho2(1)
             dvxc(ipts,2)=v2rho2(3)
             dvxc(ipts,3)=two*vsigma(1)
             dvxc(ipts,4)=two*vsigma(3)
             dvxc(ipts,5)=two*v2rhosigma(1)
             dvxc(ipts,6)=two*v2rhosigma(6)
             dvxc(ipts,7)=four*v2sigma2(1)
             dvxc(ipts,8)=four*v2sigma2(6)
           end if
         case(XC_CORRELATION)
           if (nspden==1) then
             dvxc(ipts,9)=v2rho2(1)
             dvxc(ipts,10)=dvxc(ipts,9)
             dvxc(ipts,11)=dvxc(ipts,9)
             dvxc(ipts,12)=two*vsigma(1)
             dvxc(ipts,13)=two*v2rhosigma(1)
             dvxc(ipts,14)=dvxc(ipts,13)
             dvxc(ipts,15)=four*v2sigma2(1)
           else
             dvxc(ipts,9)=v2rho2(1)
             dvxc(ipts,10)=v2rho2(2)
             dvxc(ipts,11)=v2rho2(3)
             dvxc(ipts,12)=two*vsigma(1)
             dvxc(ipts,13)=two*v2rhosigma(1)
             dvxc(ipts,14)=two*v2rhosigma(6)
             dvxc(ipts,15)=four*v2sigma2(1)
           end if
         case(XC_EXCHANGE_CORRELATION)
           msg=' KXC is not available for GGA XC_EXCHANGE_CORRELATION functionals from LibCXC'
           MSG_ERROR(msg)
         end select
       end select
     end if

     if (any(pawxcfuncs%family == XC_FAMILY_GGA) .or. any(pawxcfuncs%family == XC_FAMILY_HYB_GGA)&
&        .or. any(pawxcfuncs%family == XC_FAMILY_MGGA)) then
       !Convert the quantities returned by Libxc to the ones needed by ABINIT
       if (nspden == 1) then
         vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(1)*two
       else
         vxcgr(ipts,1) = vxcgr(ipts,1) + two*vsigma(1) - vsigma(2)
         vxcgr(ipts,2) = vxcgr(ipts,2) + two*vsigma(3) - vsigma(2)
         vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(2)
       end if
     end if
     if (any(pawxcfuncs%family == XC_FAMILY_MGGA)) then
       vxclrho(ipts,1:nspden) = vxclrho(ipts,1:nspden) + vxclrhotmp(1:nspden)
       vxctau(ipts,1:nspden) = vxctau(ipts,1:nspden) + vxctautmp(1:nspden)
     end if

   end do

 end do

end subroutine libpaw_libxc_getvxc

!----------------------------------------------------------------------

#endif

end module m_libpaw_libxc_funcs
!!***

!======================================================================
!======================================================================

!!****m* ABINIT/m_libpaw_libxc
!! NAME
!!  m_libpaw_libxc
!!
!! FUNCTION
!!  Module used to interface libPAW with host code.
!!  At present, two cases are implemented:
!!   - Use of ABINIT m_libxc_functional module
!!   - Use of embedded m_libpaw_libxc_funcs module
!!
!! COPYRIGHT
!! Copyright (C) 2014-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

module m_libpaw_libxc

#if defined LIBPAW_HAVE_LIBXC

#if defined HAVE_LIBPAW_ABINIT
 use libxc_functionals

#else
 use m_libpaw_libxc_funcs, only : &
& libxc_functionals_init       => libpaw_libxc_init, &
& libxc_functionals_fullname   => libpaw_libxc_fullname, &
& libxc_functionals_getid      => libpaw_libxc_getid, &
& libxc_functionals_global_ixc => libpaw_libxc_global_ixc, &
& libxc_functionals_getvxc     => libpaw_libxc_getvxc, &
& libxc_functionals_isgga      => libpaw_libxc_isgga, &
& libxc_functionals_ismgga     => libpaw_libxc_ismgga, &
& libxc_functionals_has_kxc    => libpaw_libxc_has_kxc, &
& libxc_functionals_nspin      => libpaw_libxc_nspin, &
& libxc_functionals_end        => libpaw_libxc_end

#endif

 implicit none

#endif

end module m_libpaw_libxc
!!***
