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
subroutine sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
     nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho,nrho,nscatterarr,nspin,spinar,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
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
  ! for real space:
!  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
  integer,intent(in):: ibyyzz_r(2,2*n2+31,2*n3+31)
  !***********************************************************************************************
   end subroutine sumrho
end interface
public :: sumrho

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
