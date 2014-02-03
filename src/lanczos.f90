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
   use ao_inguess, only: iguess_generator,ao_nspin_ig,count_atomic_shells
   use module_types
   use module_interfaces

   implicit none
   !Arguments
   real(gp) :: GetBottom
   type(atoms_data), intent(in) :: atoms
   integer, intent(in) :: nspin
   !Local variables
   character(len=*), parameter :: subname='GetBottom'
   integer, parameter :: noccmax=2,lmax=4,nelecmax=32, ng=21
   !integer, parameter :: nmax=6

   integer :: ity,  i_all
   real(gp) , pointer :: expo(:),  occup(:,:)
   real(gp)   psi(ng,5)

   integer :: i_stat
   real(gp) :: gaenes_aux(5)
   integer, dimension(lmax) :: nl
   integer nspinor, iat, noncoll

   ! if (in_iat_absorber.ne.0) then

   allocate(expo(ng +ndebug  ), stat=i_stat)
   call memocc(i_stat,expo,'expo',subname)

   allocate(occup ( noccmax  ,lmax+1+ndebug ), stat=i_stat)
   call memocc(i_stat,occup,'occup',subname)

   GetBottom=1.0e4_gp

   !for the moment, only collinear
   nspinor=1
   !if non-collinear it is like nspin=1 but with the double of orbitals
   if (nspinor == 4) then
      noncoll=2
   else
      noncoll=1
   end if

   do ity=1, atoms%astruct%ntypes
      do iat=1, atoms%astruct%nat
         if (ity.eq.atoms%astruct%iatype(iat)) exit
      end do
      call count_atomic_shells(ao_nspin_ig(nspin,nspinor=nspinor),atoms%aocc(1:,iat),occup,nl)

      call iguess_generator(atoms%nzatom(ity),atoms%nelpsp(ity),& !_modified
         &   real(atoms%nelpsp(ity),gp),atoms%psppar(0:,0:,ity),&
         &   atoms%npspcode(ity),  &
         &   atoms%nlcc_ngv(ity),atoms%nlcc_ngc(ity),atoms%nlccpar(0:,ity),&
         &   ng-1,nl,5,noccmax,lmax,occup,expo,&
         &   psi,.false., gaenes_aux=gaenes_aux  )

      if( minval(gaenes_aux ) < GetBottom) GetBottom=minval(gaenes_aux )
   enddo

   i_all=-product(shape(occup))*kind(occup)
   deallocate(occup,stat=i_stat)
   call memocc(i_stat,i_all,'occup',subname)

   i_all=-product(shape(expo))*kind(expo)
   deallocate(expo,stat=i_stat)
   call memocc(i_stat,i_all,'expo',subname)

END FUNCTION GetBottom
