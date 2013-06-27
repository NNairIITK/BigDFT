!> @file
!! libXC interfaces for the pseudo program
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining the inetrfaces with libxc library for the pseudo program
!! Test version derived from BigDFTs libABINIT/src/56_xc/m_libxc_functionals.F90 
module libxcModule

  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m

  implicit none

  type libxc_functional
    private
    integer         :: family !< LDA, GGA, etc.
    integer         :: id     !< identifier

    type(xc_f90_pointer_t) :: conf !< the pointer used to call the library
    type(xc_f90_pointer_t) :: info !< information about the functional
  end type libxc_functional

  type(libxc_functional) :: funcs(2)

  private
  public :: libxc_functionals_init, &
&      libxc_functionals_isgga, &
&      libxc_functionals_ismgga, &
&      libxc_functionals_exctXfac, &
&      libxc_functionals_end, &
!      libxc_functionals_getvxc, &
!      is replaced by a simple routine
&      xcfunction

contains


  !> Initialize the desired XC functional, from LibXC.
  subroutine libxc_functionals_init(ixc,nspden)

    implicit none

    !Arguments ------------------------------------
    !scalars
    integer, intent(in) :: nspden
    integer, intent(in) :: ixc

    !Local variables-------------------------------
    !scalars
    integer :: i


    if (ixc < 0) then
       funcs(1)%id = -ixc/1000
       funcs(2)%id = -ixc - funcs(1)%id*1000
    else
       ! DIFFERENCE to ABINIT: Only use LibXC, ignore the sign!
       funcs(1)%id = ixc/1000
       funcs(2)%id = ixc - funcs(1)%id*1000
    end if

    do i = 1, 2
      if (funcs(i)%id == 0) then
        funcs(i)%family = 0
        cycle
      end if

      ! Get XC functional family
      funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)
      select case (funcs(i)%family)
      case (XC_FAMILY_LDA, XC_FAMILY_GGA)!!!, XC_FAMILY_MGGA, XC_FAMILY_HYB_GGA)
        call xc_f90_func_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
      case default

        write(6,*)' libxc_functionals_init : ERROR -'
        write(6,*)'  Invalid IXC = ',ixc
        write(6,*)'  The LibXC functional family ',funcs(i)%family,&
             &    '  is currently not supported by pseudo.'
        write(6,*)'  (-1 means the family is unknown to the LibXC itself)'
        write(6,*)'  Please consult the LibXC documentation'
      end select

!!!!      ! Dump functional information
!!!!      call xc_f90_info_name(funcs(i)%info,message)
!!!!      call wrtout(std_out,message,'COLL')
!!!!      ii = 0
!!!!      call xc_f90_info_refs(funcs(i)%info,ii,str,message)
!!!!      do while (ii >= 0)
!!!!        call wrtout(std_out,message,'COLL')
!!!!        call xc_f90_info_refs(funcs(i)%info,ii,str,message)
!!!!      end do
    end do
  end subroutine libxc_functionals_init


  !> End usage of LibXC functional. Call LibXC end function,
  !! and deallocate module contents.
  subroutine libxc_functionals_end()

    implicit none

    integer :: i

    do i = 1, 2
      if (funcs(i)%id == 0) cycle
      call xc_f90_func_end(funcs(i)%conf)
    end do
  end subroutine libxc_functionals_end


  !> Test function to identify whether the presently used functional
  !! is a GGA or not
  function libxc_functionals_isgga()

    implicit none

    !Arguments ------------------------------------

    !Local variables-------------------------------
    logical :: libxc_functionals_isgga


    if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)) then
      libxc_functionals_isgga = .true.
    else
      libxc_functionals_isgga = .false.
    end if
  end function libxc_functionals_isgga


  real(kind=8) function libxc_functionals_exctXfac()


    implicit none

    !Arguments ------------------------------------

    !Local variables-------------------------------

    if (any(funcs%family == XC_FAMILY_HYB_GGA)) then
       !factors for the exact exchange contribution of different hybrid functionals
       if (any(funcs%id == XC_HYB_GGA_XC_PBEH)) then
          libxc_functionals_exctXfac = 0.25d0 
       end if
    else
      libxc_functionals_exctXfac = 0.d0
    end if

  end function libxc_functionals_exctXfac


  !> Test function to identify whether the presently used functional
  !! is a Meta-GGA or not
  function libxc_functionals_ismgga()

    implicit none

    !Arguments ------------------------------------

    !Local variables-------------------------------

    logical :: libxc_functionals_ismgga

    if (any(funcs%family == XC_FAMILY_MGGA)) then
      libxc_functionals_ismgga = .true.
    else
      libxc_functionals_ismgga = .false.
    end if

  end function libxc_functionals_ismgga

  
  !> Return XC potential and energy, from input density (event gradient etc...)
  !! this version calls BigDFTs XC drivers, which access LibXC as part of the ABINIT XC functions.
  !! I (AW) should really try to put this apart from BigDFT and ABINIT, and directly call libXC. 
  SUBROUTINE XCFUNCTION(nspol,rho,grad,EXC,VXC,dEdg)

!     use libxcModule

      implicit none
      integer :: nspol,i
      real(8) :: EXC,rho(nspol),VXC(nspol),dEdg(nspol),grad(nspol)  ! dummy ARGUMENTS
      real(8) :: EXCi,VXCi(nspol), sigma(3),vsigma(3)!  ! summands and libxc arg

!     These output quantities may be summed over two functionals (X and C)
      EXC  =0.0d0
      VXC  =0.0d0
      dEdg =0.0d0


       
!     convert the gradient to sigma if needed 
      if (libxc_functionals_isgga()) then
         sigma(1)=grad(1)*grad(1)
         if(nspol==2)then
            sigma(2)=grad(1)*grad(2)
            sigma(3)=grad(2)*grad(2)
         end if
      end if

!     libXC can use rather independent parts for exchange and correlation
!     the outer loop goes over both functionals from the two 3 digit codes
      do i = 1,2
          if (funcs(i)%id == 0) cycle
          select case (funcs(i)%family)
!-------------------------------------------------------------------------------
!         LDA XC 
          case (XC_FAMILY_LDA)
             call xc_f90_lda_exc_vxc(funcs(i)%conf,1,rho(1),EXCi,VXCi(1))
             EXC=EXC+EXCi
             VXC=VXC+VXCi

!-------------------------------------------------------------------------------
!         GGA XC
          case (XC_FAMILY_GGA)
             call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rho(1),sigma(1),&
                                     EXCi,VXCi(1),vsigma(1))
             EXC=EXC+EXCi
             VXC=VXC+VXCi

!            vsigma  are derivatives with respect to some products of gradients,
!            we want the derivatives with respect to the up and down gradients.
             if(nspol==1)then
               dEdg(1) = dEdg(1) +  vsigma(1)*grad(1)*2d0

             elseif(nspol==2)then
!              dE/dgup  =          dE/d(gup**2) *2*gup  + dE/d(gup*gdn) * gdn 
               dEdg(1)=dEdg(1) + vsigma(1) *2d0*grad(1) + vsigma(2)*grad(2)
!              dE/dgdn  =          dE/d(gdn**2) *2*gd   + dE/d(gup*gdn) * gup 
               dEdg(2)=dEdg(2) + vsigma(3) *2d0*grad(2) + vsigma(2)*grad(1)
             end if
          end select
      end do



!     this part is to plot the Exc(rho) function for debugging purposes
!      do j=1,nspol
!       write(17,'(i3,4f20.12)')j,rho(j),grad(j),EXC,VXC(j)
!      end do

END SUBROUTINE XCFUNCTION

end module libxcModule
