!> @file
!! libXC interfaces for the pseudo program
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module (pseudo) defining the inetrfaces with libxc library for the pseudo program
!! Test version derived from BigDFTs libABINIT/src/56_xc/m_libxc_functionals.F90
module libxcmodule
   
   use xc_f90_types_m
   use libxc_funcs_m
   use xc_f90_lib_m
   
   implicit none
   
   type libxc_functional
      private
      integer         :: family ! lda, gga, etc.
      integer         :: id     ! identifier
      
      type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
      type(xc_f90_pointer_t) :: info ! information about the functional
   end type libxc_functional
   
   type(libxc_functional) :: funcs(2)
   
   private
   public :: libxc_functionals_init, &
        &      libxc_functionals_isgga, &
        &      libxc_functionals_ismgga, &
        &      libxc_functionals_exctxfac, &
        &      libxc_functionals_end, &
        !      libxc_functionals_getvxc, &
        !      is replaced by a simple routine
   &      xcfunction
   
   contains
  

   !>  Initialize the desired xc functional, from libxc.
   !!  * call the libxc initializer
   !!  * fill preliminary fields in module structures.
   subroutine libxc_functionals_init(ixc,nspden)
      
      implicit none
      
      !arguments ------------------------------------
      !scalars
      
      integer, intent(in) :: nspden
      integer, intent(in) :: ixc
      
      !local variables-------------------------------
      !scalars
      
      integer :: i
      
      if (ixc < 0) then
         funcs(1)%id = -ixc/1000
         funcs(2)%id = -ixc - funcs(1)%id*1000
      else
         ! difference to abinit: only use libxc, ignore the sign!
         funcs(1)%id = ixc/1000
         funcs(2)%id = ixc - funcs(1)%id*1000
      end if
      
      do i = 1, 2
         if (funcs(i)%id == 0) then
            funcs(i)%family = 0
            cycle
         end if
         
         ! get xc functional family
         funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)
         select case (funcs(i)%family)
         case (xc_family_lda, xc_family_gga)!!!, xc_family_mgga, xc_family_hyb_gga)
            call xc_f90_func_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
         case default
            
            write(6,*)' libxc_functionals_init : error -'
            write(6,*)'  invalid ixc = ',ixc
            write(6,*)'  the libxc functional family ',funcs(i)%family,&
                 &    '  is currently not supported by pseudo.'
            write(6,*)'  (-1 means the family is unknown to the libxc itself)'
            write(6,*)'  please consult the libxc documentation'
         end select
         
         !!!!      ! dump functional information
         !!!!      call xc_f90_info_name(funcs(i)%info,message)
         !!!!      call wrtout(std_out,message,'coll')
         !!!!      ii = 0
         !!!!      call xc_f90_info_refs(funcs(i)%info,ii,str,message)
         !!!!      do while (ii >= 0)
         !!!!        call wrtout(std_out,message,'coll')
         !!!!        call xc_f90_info_refs(funcs(i)%info,ii,str,message)
         !!!!      end do
      end do
   end subroutine libxc_functionals_init
  

   !>  end usage of libxc functional. call libxc end function,
   !!  and deallocate module contents.
   subroutine libxc_functionals_end()
      
      implicit none
      
      integer :: i
      
      do i = 1, 2
         if (funcs(i)%id == 0) cycle
         call xc_f90_func_end(funcs(i)%conf)
      end do
   end subroutine libxc_functionals_end
  

   !> Test function to identify whether the presently used functional
   !! is a gga or not
   function libxc_functionals_isgga()
      
      implicit none
      
      !local variables
      logical :: libxc_functionals_isgga
      
      if (any(funcs%family == xc_family_gga) .or. any(funcs%family == xc_family_hyb_gga)) then
         libxc_functionals_isgga = .true.
      else
         libxc_functionals_isgga = .false.
      end if
   end function libxc_functionals_isgga
   
   
   real(kind=8) function libxc_functionals_exctxfac()
      
      
      implicit none
      
      if (any(funcs%family == xc_family_hyb_gga)) then
         !factors for the exact exchange contribution of different hybrid functionals
         if (any(funcs%id == xc_hyb_gga_xc_pbeh)) then
            libxc_functionals_exctxfac = 0.25d0 
         end if
      else
         libxc_functionals_exctxfac = 0.d0
      end if
      
   end function libxc_functionals_exctxfac
   

   !> Test function to identify whether the presently used functional
   !! is a meta-gga or not
   function libxc_functionals_ismgga()
      
      
      implicit none
      
      logical :: libxc_functionals_ismgga
      
      if (any(funcs%family == xc_family_mgga)) then
         libxc_functionals_ismgga = .true.
      else
         libxc_functionals_ismgga = .false.
      end if
      
   end function libxc_functionals_ismgga
   
   
   !> Return xc potential and energy, from input density (event gradient etc...)
   !! this version calls bigdfts xc drivers, which access libxc as part of the abinit xc functions.
   !! should really try to put this apart from bigdft and abinit, and directly call libxc. 
   subroutine xcfunction(nspol,rho,grad,exc,vxc,dedg)
      
      implicit none
      !Arguments
      integer, intent(in) :: nspol      !< Number of spin components (1 or 2)
      real(kind=8) :: rho(nspol),vxc(nspol),dedg(nspol),grad(nspol)  ! dummy arguments
      real(kind=8), intent(out) :: exc  !< Exchange-Correlation Energy
      !Local variables
      integer :: i
      real(kind=8) :: exci
      real(kind=8), dimension(nspol) :: vxci
      real(kind=8), dimension(3) :: sigma(3),vsigma(3)  ! summands and libxc arg
      ! These output quantities may be summed over two functionals (x and c)
      exc  =0.0d0
      vxc  =0.0d0
      dedg =0.0d0
      
      !     convert the gradient to sigma if needed 
      if (libxc_functionals_isgga()) then
         sigma(1)=grad(1)*grad(1)
         if(nspol==2)then
            sigma(2)=grad(1)*grad(2)
            sigma(3)=grad(2)*grad(2)
         end if
      end if
      
      !     libxc can use rather independent parts for exchange and correlation
      !     the outer loop goes over both functionals from the two 3 digit codes
      do i = 1,2
         if (funcs(i)%id == 0) cycle
         select case (funcs(i)%family)
            !-------------------------------------------------------------------------------
            !         lda xc 
         case (xc_family_lda)
            call xc_f90_lda_exc_vxc(funcs(i)%conf,1,rho(1),exci,vxci(1))
            exc=exc+exci
            vxc=vxc+vxci
            
            !-------------------------------------------------------------------------------
            !         gga xc
         case (xc_family_gga)
            call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rho(1),sigma(1),&
                 exci,vxci(1),vsigma(1))
            exc=exc+exci
            vxc=vxc+vxci
            
            !            vsigma  are derivatives with respect to some products of gradients,
            !            we want the derivatives with respect to the up and down gradients.
            if(nspol==1)then
               dedg(1) = dedg(1) +  vsigma(1)*grad(1)*2d0
               
            elseif(nspol==2)then
               !              de/dgup  =          de/d(gup**2) *2*gup  + de/d(gup*gdn) * gdn 
               dedg(1)=dedg(1) + vsigma(1) *2d0*grad(1) + vsigma(2)*grad(2)
               !              de/dgdn  =          de/d(gdn**2) *2*gd   + de/d(gup*gdn) * gup 
               dedg(2)=dedg(2) + vsigma(3) *2d0*grad(2) + vsigma(2)*grad(1)
            end if
         end select
      end do
      
      !     this part is to plot the exc(rho) function for debugging purposes
      !      do j=1,nspol
      !       write(17,'(i3,4f20.12)')j,rho(j),grad(j),exc,vxc(j)
      !      end do
      
   end subroutine xcfunction
   
end module libxcmodule
