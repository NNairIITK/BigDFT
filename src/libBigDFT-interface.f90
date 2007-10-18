interface
   subroutine CalculateTailCorrection(iproc,nproc,n1,n2,n3,rbuf,norb,norbp,nat,ntypes,&
        nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nproj,nprojel,ncongt,&
        keyv,keyg,nseg_p,keyv_p,keyg_p,nvctr_p,psppar,npspcode,eval,&
        pot,hgrid,rxyz,radii_cf,crmult,frmult,iatype,atomnames,nspin,spinar,&
        proj,psi,occup,output_grid,parallel,ekin_sum,epot_sum,eproj_sum)
     implicit none
     logical, intent(in) :: output_grid,parallel
     character(len=20), dimension(100), intent(in) :: atomnames
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,ncongt,nspin
     integer, intent(in) :: nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nproj,nprojel
     real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
     integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
     integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
     integer, dimension(2,nseg_p(2*nat)), intent(inout) :: keyg_p
     integer, dimension(ntypes), intent(in) :: npspcode
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinar
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(2*n1+31,2*n2+31,2*n3+31,nspin), intent(in) :: pot
     real(kind=8), dimension(nprojel), intent(in) :: proj
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
   end subroutine CalculateTailCorrection
end interface
public :: CalculateTailCorrection

!!$interface
!!$   subroutine orthon_p(iproc,nproc,norb,norbp,nvctrp,psit)
!!$     ! Gram-Schmidt orthogonalisation
!!$     implicit real(kind=8) (a-h,o-z)
!!$     dimension psit(nvctrp,norbp*nproc)
!!$   end subroutine orthon_p
!!$end interface

interface
   subroutine gautowav(iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        nvctr_c,nvctr_f,nseg_c,nseg_f,keyg,keyv,iatype,occup,rxyz,hgrid,psi,eks)

     implicit none
     integer, intent(in) :: norb,norbp,iproc,nproc,nat,ntypes
     integer, intent(in) :: nvctr_c,nvctr_f,n1,n2,n3,nseg_c,nseg_f
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), intent(in) :: hgrid
     real(kind=8), intent(in) :: rxyz(3,nat)
     real(kind=8), dimension(norb), intent(in) :: occup
     real(kind=8), intent(out) :: eks
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
   end subroutine gautowav
end interface
public :: gautowav

interface
   subroutine system_size(nat,rxyz,radii,rmult,iatype,ntypes, &
        cxmin,cxmax,cymin,cymax,czmin,czmax)
     ! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
     implicit real(kind=8) (a-h,o-z)
     dimension rxyz(3,nat),radii(ntypes),iatype(nat)
   end subroutine system_size
end interface
public :: system_size

interface
   subroutine make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
        ibyz_c,ibzxx_c,ibxxyy_c,ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
     ! creates complicated ib arrays
     implicit none
     integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     integer i1,i2,i3,nt,m1,m2,m3,i_stat,i_all

     integer,intent(in):: ibyz_c(2,0:n2,0:n3),ibxy_c(2,0:n1,0:n2)
     integer,intent(in):: ibyz_f(2,0:n2,0:n3),ibxy_f(2,0:n1,0:n2)


     ! for shrink:
     integer,intent(out):: ibzzx_c(2,-14:2*n3+16,0:n1) 
     integer,intent(out):: ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

     integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
     integer,intent(out):: ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
     integer,intent(out):: ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

     ! for grow:
     integer,intent(out):: ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
     integer,intent(out):: ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

     integer,intent(out):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
     integer,intent(out):: ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
     integer,intent(out):: ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

     integer,allocatable,dimension(:,:,:)::ibyx_c,ibxzz_c,ibzzyy_c
     integer,allocatable,dimension(:,:,:)::ibyx_f,ibxzz_f,ibzzyy_f

     logical,allocatable:: logrid_big(:)

     ! for real space:
     integer,intent(out):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

   end subroutine make_all_ib
end interface
public :: make_all_ib

interface
   subroutine pregion_size(rxyz,radii,rmult,iatype,ntypes, &
        hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
     ! finds the size of the smallest subbox that contains a localization region made 
     ! out of atom centered spheres
     implicit real(kind=8) (a-h,o-z)
     dimension rxyz(3),radii(ntypes)
   end subroutine pregion_size
end interface
public :: pregion_size

interface
   subroutine bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
     implicit real(kind=8) (a-h,o-z)
     logical logrid
     dimension logrid(0:n1,0:n2,0:n3)
     dimension ibyz(2,0:n2,0:n3),ibxz(2,0:n1,0:n3),ibxy(2,0:n1,0:n2)
   end subroutine bounds
end interface
public :: bounds

interface
   subroutine eleconf(nzatom,nvalelec,symbol,rcov,rprb,ehomo,neleconf,nsccode)
     implicit none
     ! Arguments
     integer, intent(in) :: nzatom,nvalelec
     character(len=2), intent(out) :: symbol
     real(kind=8), intent(out) :: rcov,rprb,ehomo
     integer, parameter :: nmax=6,lmax=3
     integer, intent(out) :: neleconf(nmax,0:lmax)
     integer, intent(out) :: nsccode
   end subroutine eleconf
end interface
public :: eleconf

interface
   subroutine readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,  & 
        nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
     ! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
     implicit real(kind=8) (a-h,o-z)
     dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
     dimension psi(nvctr_c+7*nvctr_f,norbp)
     dimension rxyz(3,nat),eval(norb),center(3)
   end subroutine readmywaves
end interface
public :: readmywaves

interface
   subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,occup,spinar)
     implicit none
     ! Arguments
     integer, intent(in) :: nelec,nspin,iproc,norb,norbu,norbd,iunit
     real(kind=8), intent(out) :: occup(norb),spinar(norb)
   end subroutine input_occup
end interface
public :: input_occup

interface
   subroutine preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,&
        ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,&
        ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
     ! Calls the preconditioner for each orbital treated by the processor
     implicit real(kind=8) (a-h,o-z)
     dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
     dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
     dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
     dimension hpsi(nvctr_c+7*nvctr_f,norbp),eval(norb)
   end subroutine preconditionall
end interface
public :: preconditionall

interface
   subroutine sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
        nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho,nrho,nscatterarr,nspin,spinar,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f)
     ! Calculates the charge density by summing the square of all orbitals
     ! Input: psi
     ! Output: rho
     implicit real(kind=8) (a-h,o-z)
     logical parallel
     dimension rho(nrho,nspin),occup(norb),spinar(norb)
     dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
     dimension psi(nvctr_c+7*nvctr_f,norbp)
     dimension nscatterarr(0:nproc-1,4)!n3d,n3p,i3s+i3xcsh-1,i3xcsh
     !***************Alexey**************************************************************************
     ! for grow:
     integer ibyz_c(2,0:n2,0:n3)
     integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
     integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

     integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
     integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
     integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
     !***********************************************************************************************
   end subroutine sumrho
end interface
public :: sumrho

!!$interface
!!$   subroutine KStrans_p(iproc,nproc,norb,ndim,nvctrp,occup,  & 
!!$        hpsit,psit,evsum,eval)
!!$     ! at the start each processor has all the Psi's but only its part of the HPsi's
!!$     ! at the end each processor has only its part of the Psi's
!!$     implicit real(kind=8) (a-h,o-z)
!!$     dimension occup(norb),eval(norb)
!!$     dimension psit(nvctrp,ndim),hpsit(nvctrp,ndim)
!!$   end subroutine KStrans_p
!!$end interface

!!$interface
!!$   subroutine orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
!!$     !Effect of orthogonality constraints on gradient 
!!$     implicit real(kind=8) (a-h,o-z)
!!$     dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc),occup(norb)
!!$   end subroutine orthoconstraint_p
!!$end interface

!!$interface
!!$   subroutine orthon(norb,norbp,nvctrp,psi)
!!$     ! Gram-Schmidt orthogonalisation
!!$     implicit real(kind=8) (a-h,o-z)
!!$     dimension psi(nvctrp,norbp)
!!$   end subroutine orthon
!!$end interface

interface
   subroutine readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,nzatom,nelpsp,&
        & psppar,npspcode,norbe,norbsc,atomnames,ntypes,iatype,iasctype,nat,natsc,&
        & scorb,norbsc_arr)
     implicit none
     ! character(len = *), intent(in) :: filename
     integer, intent(in) :: ngx, iproc, ntypes
     integer, intent(in) :: nzatom(ntypes), nelpsp(ntypes)
     real(kind=8), intent(in) :: psppar(0:4,0:6,ntypes)
     integer, intent(in) :: npspcode(ntypes),iasctype(ntypes)
     real(kind=8), intent(out) :: xp(ngx, ntypes), psiat(ngx, 5, ntypes), occupat(5, ntypes)
     integer, intent(out) :: ng(ntypes), nl(4,ntypes)
     character(len = 20), intent(in) :: atomnames(100)
     integer, intent(out) :: norbe,norbsc
     integer, intent(in) :: nat,natsc
     integer, intent(in) :: iatype(nat)
     logical, dimension(4,natsc), intent(out) :: scorb
     integer, dimension(natsc+1), intent(out) :: norbsc_arr
   end subroutine readAtomicOrbitals
end interface
public :: readAtomicOrbitals

!!$interface
!!$   subroutine checkortho(norb,norbp,nvctrp,psi)
!!$     implicit real(kind=8) (a-h,o-z)
!!$     dimension psi(nvctrp,norbp)
!!$   end subroutine checkortho
!!$end interface

interface
   subroutine HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
        psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
        nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,ngatherarr,n3p,&
        potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,spinar,&
        ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
        ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

     implicit none
     logical, intent(in) :: parallel
     character(len=1), intent(in) :: datacode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,nproj,nprojel,n3p
     integer, intent(in) :: nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nspin
     real(kind=8), intent(in) :: hgrid
     integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
     integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
     integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
     integer, dimension(ntypes), intent(in) :: npspcode
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
     integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
     integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
     real(kind=8), dimension(norb), intent(in) :: occup,spinar
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(*), intent(in) :: potential
     real(kind=8), dimension(nprojel), intent(in) :: proj
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: hpsi
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
     !********************Alexey***************************************************************
     !for shrink:
     integer ibzzx_c(2,-14:2*n3+16,0:n1) 
     integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

     integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
     integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
     integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

     !for grow:
     integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
     integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

     integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
     integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
     integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

     !for real space:
     integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
   end subroutine HamiltonianApplication
end interface
public :: HamiltonianApplication

interface
   subroutine numb_proj(ityp,ntypes,psppar,npspcode,mproj)
     ! Determines the number of projectors (valid for GTH and HGH pseudopotentials)
     implicit real(kind=8) (a-h,o-z)
     dimension psppar(0:4,0:6,ntypes),npspcode(ntypes)
   end subroutine numb_proj
end interface
public :: numb_proj

interface
   subroutine MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
        rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid,nspin)
     implicit none
     !Arguments
     logical, intent(in) :: output_grid
     integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin
     integer, dimension(nat), intent(in) :: iatype
     character(len=20), dimension(100), intent(in) :: atomnames
     real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
   end subroutine MemoryEstimator
end interface
public :: MemoryEstimator

interface
   subroutine local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
        n1,n2,n3,n3pi,i3s,rho,pot,floc)
     ! Calculates the local forces acting on the atoms belonging to iproc
     implicit none
     !Arguments---------
     integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3,n3pi,i3s
     real(kind=8), intent(in) :: hgrid
     character(len=20), dimension(100), intent(in) :: atomnames
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: rho,pot
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: nelpsp
     real(kind=8), dimension(3,nat), intent(out) :: floc
   end subroutine local_forces
end interface
public :: local_forces

interface
   subroutine projectors_derivatives(iproc,n1,n2,n3,nboxp_c,nboxp_f, & 
        ntypes,nat,norb,nprojel,nproj,&
        iatype,psppar,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,  &
        keyg,keyv,keyg_p,keyv_p,rxyz,radii_cf,cpmult,fpmult,hgrid,derproj)
     !Calculates the nonlocal forces on all atoms arising from the wavefunctions belonging to iproc and ads them to the force array

     implicit none
     !Arguments-------------
     integer, intent(in) :: iproc,ntypes,nat,norb,nprojel,nproj
     integer, intent(in) :: n1,n2,n3,nseg_c,nseg_f,nvctr_c,nvctr_f
     real(kind=8),intent(in) :: cpmult,fpmult,hgrid 
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
     integer, dimension(2,3,nat), intent(in) :: nboxp_c,nboxp_f
     integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
     integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(nprojel), intent(in) :: proj
     real(kind=8), dimension(nprojel,3), intent(out) :: derproj
   end subroutine projectors_derivatives
end interface
public :: projectors_derivatives

!!$interface
!!$   subroutine KStrans(norb,norbp,nvctrp,occup,hpsi,psi,evsum,eval)
!!$     ! at the start each processor has all the Psi's but only its part of the HPsi's
!!$     ! at the end each processor has only its part of the Psi's
!!$     implicit real(kind=8) (a-h,o-z)
!!$     dimension occup(norb),eval(norb)
!!$     dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp)
!!$   end subroutine KStrans
!!$end interface

interface
   subroutine createAtomicOrbitals(iproc, nproc, atomnames,&
        & nat, rxyz, norbe, norbep, norbsc, occupe, occupat, ngx, xp, psiat, ng, nl, &
        & nvctr_c, nvctr_f, n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, nseg_c, nseg_f, &
        & keyg, keyv, iatype, ntypes, iasctype, natsc, psi, eks, scorb)

     implicit none
     integer, intent(in) :: nat, norbe, norbep, ngx, iproc, nproc
     integer, intent(in) :: nvctr_c, nvctr_f, n1, n2, n3, nseg_c, nseg_f
     integer, intent(in) :: nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, ntypes
     integer, intent(in) :: norbsc,natsc
     logical, dimension(4,natsc), intent(in) :: scorb
     integer, intent(in) :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
     integer, intent(in) :: iatype(nat),iasctype(ntypes)
     real(kind=8), intent(in) :: hgrid
     real(kind=8), intent(out) :: eks
     !character(len = 20), intent(in) :: pspatomnames(npsp)
     character(len = 20), intent(in) :: atomnames(100)
     integer, intent(inout) :: ng(ntypes), nl(4,ntypes)
     real(kind=8), intent(in) :: rxyz(3, nat)
     real(kind=8), intent(inout) :: xp(ngx, ntypes), psiat(ngx, 5, ntypes)
     real(kind=8), intent(inout) :: occupat(5, ntypes)
     real(kind=8), intent(out) :: psi(nvctr_c + 7 * nvctr_f, norbep), occupe(norbe)
   end subroutine createAtomicOrbitals
end interface
public :: createAtomicOrbitals

interface
   subroutine switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psiw)
     implicit none
     integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(nvctrp,norbp,nproc), intent(out) :: psiw
   end subroutine switch_waves
end interface
public :: switch_waves

interface
   subroutine createIonicPotential(iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,hgrid,&
        elecfield,n1,n2,n3,n3pi,i3s,pkernel,pot_ion,eion)
     implicit none
     integer, intent(in) :: iproc,nproc,nat,ntypes,n1,n2,n3,n3pi,i3s
     real(kind=8), intent(in) :: hgrid,elecfield
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: nelpsp
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), intent(out) :: eion
     real(kind=8), dimension(*), intent(out) :: pot_ion
   end subroutine createIonicPotential
end interface
public :: createIonicPotential

interface
   subroutine nonlocal_forces(iproc,ntypes,nat,norb,norbp,nprojel,nproj,&
        iatype,psppar,npspcode,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,derproj,  &
        keyg,keyv,keyg_p,keyv_p,psi,fsep)
     !Calculates the nonlocal forces on all atoms arising from the wavefunctions belonging to iproc and adds them to the force array

     implicit none
     !Arguments-------------
     integer, intent(in) :: iproc,ntypes,nat,norb,norbp,nprojel,nproj
     integer, intent(in) :: nseg_c,nseg_f,nvctr_c,nvctr_f
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
     integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
     integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     integer, dimension(ntypes), intent(in) :: npspcode
     real(kind=8), dimension(norb), intent(in) :: occup
     real(kind=8), dimension(nprojel), intent(in) :: proj
     real(kind=8), dimension(nprojel,3), intent(in) :: derproj
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(3,nat), intent(inout) :: fsep
   end subroutine nonlocal_forces
end interface
public nonlocal_forces

!!$interface
!!$   subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
!!$     ! Calculates the keys describing a wavefunction data structure
!!$     implicit real(kind=8) (a-h,o-z)
!!$     logical logrid,plogrid
!!$     dimension logrid(0:n1,0:n2,0:n3),keyg(2,mseg),keyv(mseg)
!!$   end subroutine segkeys
!!$end interface

interface
   subroutine unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,psi)
     implicit none
     integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp
     real(kind=8), dimension(nvctrp,norbp,nproc), intent(in) :: psiw
     real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
   end subroutine unswitch_waves
end interface
public unswitch_waves

!!$interface
!!$   subroutine orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
!!$     !Effect of orthogonality constraints on gradient 
!!$     implicit real(kind=8) (a-h,o-z)
!!$     logical, parameter :: parallel=.false.
!!$     dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp),occup(norb)
!!$   end subroutine orthoconstraint
!!$end interface

interface
   subroutine reformatmywaves(iproc, norb, norbp, nat, &
        & hgrid_old, nvctr_c_old, nvctr_f_old, n1_old, n2_old, n3_old, rxyz_old, &
        & nseg_c_old, nseg_f_old, keyg_old, keyv_old, psi_old, &
        & hgrid, nvctr_c, nvctr_f, n1, n2, n3, rxyz, &
        & nseg_c, nseg_f, keyg, keyv, psi)
     implicit real(kind=8) (a-h,o-z)
     dimension :: rxyz(3,nat), rxyz_old(3,nat), center(3), center_old(3)
     dimension :: keyg_old(2, nseg_c_old + nseg_f_old), keyv_old(nseg_c_old + nseg_f_old)
     dimension :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
     dimension :: psi_old(nvctr_c_old + 7 * nvctr_f_old, norbp), psi(nvctr_c + 7 * nvctr_f, norbp)
   end subroutine reformatmywaves
end interface
public :: reformatmywaves

interface
   subroutine crtproj(iproc,nterm,n1,n2,n3, & 
        nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
        radius_f,cpmult,fpmult,hgrid,gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
        mvctr_c,mvctr_f,proj_c,proj_f)
     ! returns the compressed form of a Gaussian projector 
     ! x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
     ! in the arrays proj_c, proj_f
     implicit real(kind=8) (a-h,o-z)
     parameter(ntermx=3,nw=16000)
     dimension lx(nterm),ly(nterm),lz(nterm)
     dimension fac_arr(nterm)
     dimension proj_c(mvctr_c),proj_f(7,mvctr_f)
   end subroutine crtproj
end interface
public :: crtproj

interface
   subroutine writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
        nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
     ! write all my wavefunctions in files by calling writeonewave
     implicit real(kind=8) (a-h,o-z)
     character(len=4) f4
     character(len=50) filename
     dimension rxyz(3,nat),eval(norb),center(3)
     dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
     dimension psi(nvctr_c+7*nvctr_f,norbp)
   end subroutine writemywaves
end interface
public :: writemywaves

interface
   subroutine calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)
     implicit none
     integer, intent(in) :: l,i,m,nterm_max
     integer, intent(out) :: nterm
     integer, dimension(nterm_max), intent(out) :: lx,ly,lz
     real(kind=8), dimension(nterm_max), intent(out) :: fac_arr
   end subroutine calc_coeff_proj
end interface
public :: calc_coeff_proj

interface
   subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
     ! calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
     implicit real(kind=8) (a-h,o-z)
     dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
   end subroutine wnrm
end interface
public :: wnrm

interface
   subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
     ! Calculates the length of the keys describing a wavefunction data structure
     implicit real(kind=8) (a-h,o-z)
     logical logrid,plogrid
     dimension logrid(0:n1,0:n2,0:n3)
   end subroutine num_segkeys
end interface
public :: num_segkeys

interface
   subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
        ntypes,iatype,rxyz,radii,rmult,hgrid,logrid)
     ! set up an array logrid(i1,i2,i3) that specifies whether the grid point
     ! i1,i2,i3 is the center of a scaling function/wavelet
     implicit real(kind=8) (a-h,o-z)
     logical logrid
     parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
     dimension rxyz(3,nat),iatype(nat),radii(ntypes)
     dimension logrid(0:n1,0:n2,0:n3)
   end subroutine fill_logrid
end interface
public fill_logrid
