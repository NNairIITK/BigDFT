
subroutine updatePotential(ixc,nspin,denspot,ehart,eexcu,vexcu)
!
! Purpose:
! ========
!   Calculates the potential and energy and writes them. This is subroutine is copied
!   from cluster.
!
! Calling arguments:
! ==================
!   Input arguments:
!   -----------------
!     iproc       process ID
!     nproc       total number of processes
!     n3d         ??
!     n3p         ??
!     Glr         type describing the localization region
!     orbs        type describing the physical orbitals psi
!     atoms       type containing the parameters for the atoms
!     in          type  containing some very general parameters
!     lin         type containing parameters for the linear version
!     psi         the physical orbitals
!     rxyz        atomic positions
!     rhopot      the charge density
!     nscatterarr ??
!     nlpspd      ??
!     proj        ??
!     pkernelseq  ??
!     radii_cf    coarse and fine radii around the atoms
!     irrzon      ??
!     phnons      ??
!     pkernel     ??
!     pot_ion     the ionic potential
!     rhocore     ??
!     potxc       ??
!     PSquiet     flag to control the output from the Poisson solver
!     eion        ionic energy
!     edisp       dispersion energy
!     fion        ionic forces
!     fdisp       dispersion forces
!   Input / Output arguments
!   ------------------------
!     rhopot      the charge density
!   Output arguments:
!   -----------------


use module_base
use module_types
use module_interfaces, exceptThisOne => updatePotential
use Poisson_Solver
implicit none

! Calling arguments
integer, intent(in) :: ixc,nspin
type(DFT_local_fields), intent(inout) :: denspot
real(8),intent(out):: ehart, eexcu, vexcu

! Local variables
character(len=*), parameter :: subname='updatePotential'
logical :: nullifyVXC
real(8):: ekin_sum, epot_sum, eproj_sum, energybs, energyMod
real(8):: energyMod2, ehartMod, t1, t2, time
integer:: istat, iall, infoCoeff, ilr, ierr, sizeLphir, sizePhibuffr
real(dp), dimension(6) :: xcstr

nullifyVXC=.false.

if(nspin==4) then
   !this wrapper can be inserted inside the poisson solver 
   call PSolverNC(denspot%pkernel%geocode,'D',denspot%pkernel%iproc,denspot%pkernel%nproc,&
        denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
        denspot%dpbox%n3d,ixc,&
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

   call XC_potential(denspot%pkernel%geocode,'D',denspot%pkernel%iproc,denspot%pkernel%nproc,&
        denspot%pkernel%mpi_comm,&
        denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),ixc,&
        denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
        denspot%rhov,eexcu,vexcu,nspin,denspot%rho_C,denspot%V_XC,xcstr)
   
   call H_potential('D',denspot%pkernel,denspot%rhov,denspot%V_ext,ehart,0.0_dp,.true.,&
        quiet=denspot%PSquiet) !optional argument
   
   !sum the two potentials in rhopot array
   !fill the other part, for spin, polarised
   if (nspin == 2) then
      call dcopy(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,denspot%rhov(1),1,&
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

end subroutine updatePotential


