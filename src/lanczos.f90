!> @file
!!  Lanczos diagonalisation used by XANES calculation
!! @author
!!    Copyright (C) 2009-2011 BigDFT group (AM, LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine dirac_hara (rho, E , V)
   use module_base

   implicit none

   !Arguments
   real(gp), intent(in   ) :: rho, E
   real(gp), intent(inout) :: V
   !Local variables
   real(gp), parameter :: f= 1.919158d0  !( 4/(9*pi)**(1/3)  )
   real(gp), parameter :: pi=3.141592653589793_gp

   real(gp) :: Vcorr, rs, xk, EV,x
   integer :: i

   if(rho>1.0e-4) then
      rs = (3.0_gp / (4.0_gp*pi*rho)) ** (1.0_gp/3.0_gp)
   else
      rs=1000.0_gp
   endif

   Vcorr=V

   EV=E-Vcorr
   if(EV<=0) then
      return
   endif

   do i=1,10

      EV=E-Vcorr

      if(EV<=0) then
         return
      endif

      xk =sqrt(2*EV  )

      x = xk *rs/ f

      if ((x-1)  < 1.0D-6) return
      Vcorr =V - (f/pi/rs) * ( log(abs((1+x) / (1-x))) * (1-x**2) / (2*x))

   end do
   V=Vcorr
   return
END SUBROUTINE dirac_hara


function GetBottom(atoms,nspin)

   use module_base
   use ao_inguess, only: &
        iguess_generator,ao_nspin_ig,nmax_occ_ao
   use module_types
   use module_interfaces

   implicit none
   !Arguments
   real(gp) :: GetBottom
   type(atoms_data), intent(in) :: atoms
   integer, intent(in) :: nspin
   !Local variables
   character(len=*), parameter :: subname='GetBottom'
   integer, parameter :: ng=21
   !integer, parameter :: nmax=6

   integer :: ity
   real(gp) , pointer :: expo(:)
   real(gp), dimension(ng,10) :: psi

   real(gp), dimension(nmax_occ_ao) :: gaenes_aux
   integer :: nspinor, iat, nspin_ig

   ! if (in_iat_absorber.ne.0) then

   expo = f_malloc_ptr(ng ,id='expo')

   GetBottom=1.0e4_gp

   nspin_ig=ao_nspin_ig(nspin,nspinor=nspinor)

   do ity=1, atoms%astruct%ntypes
      do iat=1, atoms%astruct%nat
         if (ity.eq.atoms%astruct%iatype(iat)) exit
      end do
      call iguess_generator(atoms%nzatom(ity),atoms%nelpsp(ity),&
           real(atoms%nelpsp(ity),gp),nspin_ig,atoms%aoig(iat)%aocc,atoms%psppar(0:,0:,ity),&
           atoms%npspcode(ity),  &
           atoms%nlcc_ngv(ity),atoms%nlcc_ngc(ity),atoms%nlccpar(0:,ity),&
           ng-1,expo,psi,.false., gaenes_aux=gaenes_aux  )

      if( minval(gaenes_aux ) < GetBottom) GetBottom=minval(gaenes_aux )
   enddo

   call f_free_ptr(expo)

END FUNCTION GetBottom
