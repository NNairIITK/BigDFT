! Standalone test routines to compile without LibXC.
! This is a modified version of the module that usually
! interfaces the libXC functionals. This test version 
! simply uses the Teter93 closed shell LDA XC
! from the previous versions of pseudo.
! This should equal "pade" LDA as implemented in CP2K. 


module libxcModule
  implicit none

  type libxc_functional
    private
    integer         :: family ! LDA, GGA, etc.
    integer         :: id     ! identifier

  end type libxc_functional

  type(libxc_functional) :: funcs(2)

  private
  public :: libxc_functionals_init, &
&      libxc_functionals_isgga, &
&      libxc_functionals_ismgga, &
&      libxc_functionals_exctXfac, &
&      libxc_functionals_end,&
&      xcfunction

contains

  subroutine libxc_functionals_init(ixc,nspden)



    implicit none

!Arguments ------------------------------------
!scalars

    integer, intent(in) :: nspden
    integer, intent(in) :: ixc

!Local variables-------------------------------
!scalars

    integer :: i, ii
    character(len=500) :: message

! *************************************************************************

    write(6,*)'NOTE: The program was compiled without libXC.'
    if (ixc == -20 .and. nspden ==1 ) then
       write(6,*)'The Teter93 LDA from libXCfake will be used.'
    else
       write(6,*)'Only the closed shell Teter93 LDA functional is supported.'
       write(6,*)'Cannot proceed. Choose iXC=-20 and nspden=1'
       write(6,*)'or build with libXC to use another functional.'
       stop
    end if
  end subroutine libxc_functionals_init

  subroutine libxc_functionals_end()
  ! the fake does not need any initialization
  end subroutine libxc_functionals_end
!!*** 

  function libxc_functionals_isgga()
    ! the fake is never GGA
    implicit none
    logical :: libxc_functionals_isgga
    libxc_functionals_isgga = .false.
  end function libxc_functionals_isgga
!!*** 


  function libxc_functionals_ismgga()
    ! the fake is never meta-GGA
    implicit none
    logical :: libxc_functionals_ismgga
    libxc_functionals_ismgga = .false.
  end function libxc_functionals_ismgga
!!*** 



SUBROUTINE XCFUNCTION(nspol,rho,grad,EXC,VXC,dEdg)
! simply uses the Teter93 closed shell LDA XC
! from the previous versions of pseudo.
! This equals "pade" LDA as implemented in CP2K. 

!     use libxcModule

      IMPLICIT REAL*8 (A-H,O-Z)
      integer :: nspol,i,j
      real(8) :: abs
      real(8) :: EXC,rho(nspol),VXC(nspol),dEdg(nspol),grad(nspol)  ! dummy ARGUMENTS
      real(8) :: EXCi,VXCi(nspol), sigma(3),vsigma(3)!  ! summands and libxc arg

!-------------------------------------------------------------------------------
!    oLDA XC 
      PARAMETER (A0=0.4581652932831429D0,A1=2.217058676663745D0,&
                 A2=0.7405551735357053D0,A3=0.01968227878617998D0)
      PARAMETER (B1=1.0000000000000000D0,B2=4.504130959426697D0, &
                 B3=1.110667363742916D0,B4=0.02359291751427506D0)
      PARAMETER (O3=1.D0/3.D0)
!     ==--------------------------------------------------------------==

      
      RS=  0.75D0 / 3.141592653589793D+00 / rho(1)
      RS= RS**(1d0/3d0)

      TOP=A0+RS*(A1+RS*(A2+RS*A3))
      DTOP=A1+RS*(2.D0*A2+3.D0*A3*RS)
      BOT=RS*(B1+RS*(B2+RS*(B3+RS*B4)))
      DBOT=B1+RS*(2.D0*B2+RS*(3.D0*B3+RS*4.D0*B4))
      EXC=-TOP/BOT
      VXC=EXC+RS*O3*(DTOP/BOT-TOP*DBOT/(BOT*BOT))
      dEdg=0d0
      write(99,*)rho,Exc,Vxc
      return
      end subroutine

end module 
