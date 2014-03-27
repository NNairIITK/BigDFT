!> @file 
!!   energy and forces in linear
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
 

!> Calculates the potential and energy and writes them. This is subroutine is copied
!! from cluster.
subroutine updatePotential(nspin,denspot,ehart,eexcu,vexcu)

use module_base
use module_types
use module_interfaces, exceptThisOne => updatePotential
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
implicit none

! Calling arguments
integer, intent(in) :: nspin
type(DFT_local_fields), intent(inout) :: denspot
real(kind=8),intent(out) :: ehart, eexcu, vexcu

! Local variables
character(len=*), parameter :: subname='updatePotential'
logical :: nullifyVXC
integer :: istat, iall
real(dp), dimension(6) :: xcstr

nullifyVXC=.false.

if(nspin==4) then
   !this wrapper can be inserted inside the poisson solver 
   call PSolverNC(denspot%pkernel%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
        denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
        denspot%dpbox%n3d,denspot%xc,&
        denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
        denspot%rhov,denspot%pkernel%kernel,denspot%V_ext,ehart,eexcu,vexcu,0.d0,.true.,4)

else
   if (.not. associated(denspot%V_XC)) then   
      !Allocate XC potential
      if (denspot%dpbox%n3p >0) then
         allocate(denspot%V_XC(denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%n3p,nspin+ndebug),stat=istat)
         call memocc(istat,denspot%V_XC,'denspot%V_XC',subname)
      else
         allocate(denspot%V_XC(1,1,1,1+ndebug),stat=istat)
         call memocc(istat,denspot%V_XC,'denspot%V_XC',subname)
      end if
      nullifyVXC=.true.
   end if

   call XC_potential(denspot%pkernel%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
        denspot%pkernel%mpi_env%mpi_comm,&
        denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),denspot%xc,&
        denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
        denspot%rhov,eexcu,vexcu,nspin,denspot%rho_C,denspot%V_XC,xcstr)
    
   call H_potential('D',denspot%pkernel,denspot%rhov,denspot%V_ext,ehart,0.0_dp,.true.,&
        quiet=denspot%PSquiet) !optional argument
   
   !sum the two potentials in rhopot array
   !fill the other part, for spin, polarised
   if (nspin == 2) then
      call vcopy(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,denspot%rhov(1),1,&
           denspot%rhov(1+denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p),1)
   end if
   !spin up and down together with the XC part
   call axpy(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p*nspin,1.0_dp,denspot%V_XC(1,1,1,1),1,&
        denspot%rhov(1),1)
   
   if (nullifyVXC) then
      iall=-product(shape(denspot%V_XC))*kind(denspot%V_XC)
      deallocate(denspot%V_XC,stat=istat)
      call memocc(istat,iall,'denspot%V_XC',subname)
   end if

end if

END SUBROUTINE updatePotential
