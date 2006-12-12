!!****f* BigDFT/Parxc_energy
!! NAME
!! xc_energy
!!
!! FUNCTION
!! Calculates the exchange-correlation energies, parallel case
!! $ E_{xc} = \int d^3 x \rho(x) e_{xc}(x) $
!! $ V_{xc} = \int d^3 x \rho(x) v_{xc}(x) $
!! where e_xc and v_xc are calculated from ABINIT routines 
!! with the exchance-correlation method chosen
!! Then adds the potential $V_{xc}$ to the distributed ionic potential
!! Substitutes the routine enterdensity of the non-xc case
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme (see ABINIT -> drivexc documentation)
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  nspden= value of the spin polarization
!!  rhopot= charge density
!!  n01,n02,n03= dimension of the real space
!!  hgrid=grid spacing
!!
!! OUTPUT
!!  eht=Hartree energy $E_{h}=1/2\int \rho* V$
!!  exc=$E_{xc}$ defined before
!!  vxc=$V_{xc}$ defined before
!!
!! SOURCE
!!
subroutine Parxc_energy(m1,m2,m3,md1,md2,md3,xcdim,wbdim,mx2,deltaleft,deltaright,wbl,wbr,&
     ixc,hgrid,rhopot,pot_ion,zf,zfpot_ion,exc,vxc,iproc,nproc)
  
  use defs_xc
   
  implicit none
  
  !Arguments----------------------
  integer, intent(in) :: m1,m2,m3,xcdim,wbdim,wbl,wbr,mx2,md1,md2,md3,ixc,iproc,nproc
  integer, intent(in) :: deltaleft,deltaright
  real(kind=8), intent(in) :: hgrid
  real(kind=8), dimension(m1,m3,mx2), intent(inout) :: rhopot
  real(kind=8), dimension(m1,m3,m2), intent(in) :: pot_ion
  real(kind=8), dimension(md1,md3,md2/nproc), intent(out) :: zf,zfpot_ion
  real(kind=8), intent(out) :: exc,vxc
  
  !Local variables----------------
  real(kind=8), dimension(:,:,:), allocatable :: exci,d2vxci
  real(kind=8), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(kind=8), dimension(:,:,:,:,:), allocatable :: gradient
  real(kind=8) :: elocal,vlocal,rho,pot,potion
  integer :: npts,i_all,nspden,order,offset
  integer :: i1,i2,i3,j1,j2,j3,jp2,jpp2,jppp2
  integer :: ndvxc,nvxcdgr
  
  !Body

  !these are always the same
  nspden=1
  order=1

  !Allocations
  allocate(exci(m1,m3,wbdim),stat=i_all)
  if (i_all /=0) stop 'allocation error (exci)'
  allocate(vxci(m1,m3,wbdim,nspden),stat=i_all)
  if (i_all /=0) stop 'allocation error (vxci)'
  !Allocations of the exchange-correlation terms, depending on the ixc value
  if (ixc >= 11 .and. ixc <= 16) allocate(gradient(m1,m3,wbdim,2*nspden-1,0:3),stat=i_all)   
  if (i_all /=0) stop 'allocation error, (gradient)'

  call size_dvxc(ixc,ndvxc,nspden,nvxcdgr,order)

  if (ndvxc/=0) allocate(dvxci(m1,m3,wbdim,ndvxc),stat=i_all)
  if (i_all /=0) stop 'allocation error (dvxci)'
  if (nvxcdgr/=0) allocate(dvxcdgr(m1,m3,wbdim,nvxcdgr),stat=i_all)
  if (i_all /=0) stop 'allocation error (dvxcdgr)'
  if ((ixc==3 .or. (ixc>=7 .and. ixc<=15)) .and. order==3) &
       allocate(d2vxci(m1,m3,wbdim),stat=i_all)
  if (i_all /=0) stop 'allocation error, (d2vxci)'

  if (.not.allocated(gradient) .and. xcdim/=mx2 ) then
     print *,'Parxc_energy: if nx2/=xcdim the gradient must be allocated'
     stop
  end if

 !divide by two the density to applicate it in the ABINIT xc routines
 
  do i3=1,mx2
     do i2=1,m3
        do i1=1,m1
           rhopot(i1,i2,i3)=.5d0*rhopot(i1,i2,i3)
        end do
     end do
  end do

  !computation of the gradient

  if (allocated(gradient)) then

     call calc_gradient(m1,m3,mx2,wbdim,deltaleft,deltaright,rhopot,nspden,&
          hgrid,hgrid,hgrid,gradient)

  end if


  offset=deltaleft+1
  npts=m1*m3*wbdim
  !let us apply ABINIT routines
  !case with gradient
  if (ixc >= 11 .and. ixc <= 16) then
     if (order**2 <= 1 .or. ixc == 16) then
        if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,     & 
                &grho2_updn=gradient,vxcgr=dvxcdgr) 
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,     & 
                &grho2_updn=gradient) 
        end if
     else if (order /= 3) then
        if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,     & 
                &dvxc=dvxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,     & 
                &dvxc=dvxci,grho2_updn=gradient) 
        end if
     else if (order == 3) then
        if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,     & 
                &dvxci,d2vxci,gradient,dvxcdgr) 
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,     & 
                &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient) 
        end if
     end if
 
     !do not calculate the White-Bird term in the Leeuwen Baerends XC case
     if (ixc/=13) call vxcpostprocessing(m1,m3,wbdim,xcdim,wbl,wbr,nspden,nvxcdgr,gradient,&
          hgrid,hgrid,hgrid,dvxcdgr,vxci)

     !cases without gradient
  else
     if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
        call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr)
     else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
        call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,      &
             &dvxc=dvxci,d2vxc=d2vxci)
     else
        call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,nvxcdgr,      &
             &dvxc=dvxci)
     end if
  end if
  
  !distributing the density in the zf array
  !calculating the xc integrated quantities
  !and summing the xc potential into the pot_ion array
 
  exc=0.d0
  vxc=0.d0
  do jp2=1,xcdim
     j2=offset+jp2+wbl-2
     jpp2=iproc*(md2/nproc)+jp2
     jppp2=jp2+wbl-1
     do j3=1,m3
        do j1=1,m1
           rho=rhopot(j1,j3,j2)
           potion=pot_ion(j1,j3,jpp2)
           if (rho < 5.d-21) then
              elocal=0.d0
              vlocal=0.d0
           else
              elocal=exci(j1,j3,jppp2)
              vlocal=vxci(j1,j3,jppp2,1)
           end if
           exc=exc+elocal*rho
           vxc=vxc+vlocal*rho
           zf(j1,j3,jp2)=2.d0*rhopot(j1,j3,j2)
           zfpot_ion(j1,j3,jp2)=potion+vlocal
        end do
        do j1=m1+1,md1
           zf(j1,j3,jp2)=0.d0
           zfpot_ion(j1,j3,jp2)=0.d0
        end do
     end do
     do j3=m3+1,md3
        do j1=1,md1
           zf(j1,j3,jp2)=0.d0
           zfpot_ion(j1,j3,jp2)=0.d0
        end do
     end do
  end do
  do jp2=xcdim+1,md2/nproc
     do j3=1,md3
        do j1=1,md1
           zf(j1,j3,jp2)=0.d0
           zfpot_ion(j1,j3,jp2)=0.d0
        end do
     end do
  end do

  !the two factor is due to the 
  !need of using the density of states in abinit routines
  exc=2.d0*hgrid**3*exc
  vxc=2.d0*hgrid**3*vxc


  !De-allocations
  deallocate(exci,vxci)
  if (allocated(dvxci)) deallocate(dvxci)
  if (allocated(dvxcdgr)) deallocate(dvxcdgr)
  if (allocated(d2vxci)) deallocate(d2vxci)
  if (allocated(gradient)) deallocate(gradient)

end subroutine Parxc_energy



!!****f* BigDFT/xc_energy
!! NAME
!! xc_energy
!!
!! FUNCTION
!! Calculates the exchange-correlation energies
!! $ E_{xc} = \int d^3 x \rho(x) e_{xc}(x) $
!! $ V_{xc} = \int d^3 x \rho(x) v_{xc}(x) $
!! where e_xc and v_xc are calculated from ABINIT routines 
!! with the exchance-correlation method chosen
!! Then adds the potential $V_{xc}$ and the ionic
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme (see ABINIT -> drivexc documentation)
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  nspden= value of the spin polarization
!!  rhopot= charge density
!!  n01,n02,n03= dimension of the real space
!!  hgrid=grid spacing
!!
!! OUTPUT
!!  eht=Hartree energy $E_{h}=1/2\int \rho* V$
!!  exc=$E_{xc}$ defined before
!!  vxc=$V_{xc}$ defined before
!!
!! SOURCE
!!
subroutine xc_energy(n01,n02,n03,nd1,nd2,nd3,ixc,factor,hgrid,rhopot,pot_ion,zarray,&
     eht,exc,vxc)
  
  use defs_xc
   
  implicit none
  
  !Arguments----------------------
  integer, intent(in) :: n01,n02,n03,nd1,nd2,nd3,ixc
  real(kind=8), intent(in) :: hgrid,factor
  real(kind=8), dimension(n01,n02,n03), intent(inout) :: rhopot
  real(kind=8), dimension(n01,n02,n03), intent(in) :: pot_ion
  real(kind=8), dimension(nd1,nd2,nd3), intent(in) :: zarray
  real(kind=8), intent(out) :: eht,exc,vxc
  
  !Local variables----------------
  real(kind=8), dimension(:,:,:), allocatable :: exci,d2vxci
  real(kind=8), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(kind=8), dimension(:,:,:,:,:), allocatable :: gradient
  real(kind=8) :: elocal,vlocal,rho,pot,potion
  integer :: npts,i_all,nspden,order
  integer :: i1,i2,i3
  integer :: ndvxc,nvxcdgr,order_grad
  
  !Body

  !these are always the same
  nspden=1
  order=1

  !order of the finite difference first derivative defining the gradient
  order_grad=4
  
  !Allocations
  allocate(exci(n01,n02,n03),stat=i_all)
  if (i_all /=0) stop 'allocation error (exci)'
  allocate(vxci(n01,n02,n03,nspden),stat=i_all)
  if (i_all /=0) stop 'allocation error (vxci)'
  !Allocations of the exchange-correlation terms, depending on the ixc value
  if (ixc >= 11 .and. ixc <= 16) then
     allocate(gradient(n01,n02,n03,2*nspden-1,0:3),stat=i_all)  
  end if
  if (i_all /=0) stop 'allocation error, (gradient)'
  call size_dvxc(ixc,ndvxc,nspden,nvxcdgr,order)
  if (ndvxc/=0) allocate(dvxci(n01,n02,n03,ndvxc),stat=i_all)
  if (i_all /=0) stop 'allocation error (dvxci)'
  if (nvxcdgr/=0) allocate(dvxcdgr(n01,n02,n03,nvxcdgr),stat=i_all)
  if (i_all /=0) stop 'allocation error (dvxcdgr)'
  if ((ixc==3 .or. (ixc>=7 .and. ixc<=15)) .and. order==3) &
       allocate(d2vxci(n01,n02,n03),stat=i_all)
  if (i_all /=0) stop 'allocation error, (d2vxci)'

  !divide by two the density to applicate it in the ABINIT xc routines
 
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           rhopot(i1,i2,i3)=.5d0*rhopot(i1,i2,i3)
        end do
     end do
  end do


  !computation of the gradient

  if (allocated(gradient)) then
     
     call calc_gradient(n01,n02,n03,n03,0,0,rhopot,nspden,&
          hgrid,hgrid,hgrid,gradient)
     
  end if

  npts=n01*n02*n03
!let us apply ABINIT routines
!case with gradient
  if (ixc >= 11 .and. ixc <= 16) then
     if (order**2 <= 1 .or. ixc == 16) then
        if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,     & 
                &grho2_updn=gradient,vxcgr=dvxcdgr) 
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,     & 
                &grho2_updn=gradient) 
        end if
     else if (order /= 3) then
        if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,     & 
                &dvxc=dvxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,     & 
                &dvxc=dvxci,grho2_updn=gradient) 
        end if
     else if (order == 3) then
        if (ixc /= 13) then             
           call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,     & 
                &dvxci,d2vxci,gradient,dvxcdgr) 
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,     & 
                &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient) 
        end if
     end if

     !do not calculate the White-Bird term in the Leeuwen Baerends XC case
     if (ixc/=13) call vxcpostprocessing(n01,n02,n03,n03,1,1,nspden,nvxcdgr,gradient,&
          hgrid,hgrid,hgrid,dvxcdgr,vxci)

  !cases without gradient
  else
     if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
        call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr)
     else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
        call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,      &
             &dvxc=dvxci,d2vxc=d2vxci)
     else
        call drivexc(exci,ixc,npts,nspden,order,rhopot,vxci,ndvxc,nvxcdgr,      &
             &dvxc=dvxci)
     end if
  end if

  eht=0.d0
  exc=0.d0
  vxc=0.d0
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           rho=rhopot(i1,i2,i3)
           pot=factor*zarray(i1,i2,i3)
           potion=pot_ion(i1,i2,i3)
           eht=eht+rho*pot
           if (rho < 5.d-21) then
              elocal=0.d0
              vlocal=0.d0
           else
              elocal=exci(i1,i2,i3)
              vlocal=vxci(i1,i2,i3,1)
           end if
           exc=exc+elocal*rho
           vxc=vxc+vlocal*rho
           rhopot(i1,i2,i3)=pot+vlocal+potion
        end do
     end do
  end do

  !the two factor is due to the 
  !need of using the density of states in abinit routines
  eht=hgrid**3*eht
  exc=2.d0*hgrid**3*exc
  vxc=2.d0*hgrid**3*vxc

  !De-allocations
  deallocate(exci,vxci)
  if (allocated(dvxci)) deallocate(dvxci)
  if (allocated(dvxcdgr)) deallocate(dvxcdgr)
  if (allocated(d2vxci)) deallocate(d2vxci)
  if (allocated(gradient)) deallocate(gradient)

end subroutine xc_energy


!!****f* BigDFT/vxcpostprocessing
!! NAME
!! vxcpostprocessing
!!
!! FUNCTION
!! Correct the XC potential with the White-Bird formula, to be used for the 
!! GGA case. Works either in parallel of in serial, by proper change of the 
!! arguments.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2006 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!! This routine has been written from rhohxc_coll(DCA, XG, GMR, MF, GZ)
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!
!! OUTPUT
!! 
!!
!! SOURCE
subroutine vxcpostprocessing(n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr,gradient,hx,hy,hz,dvxcdgr,wb_vxc)
  implicit none
  integer, intent(in) :: n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(n01,n02,n03,2*nspden-1,0:3), intent(in) :: gradient
  real(kind=8), dimension(n01,n02,n03,nvxcdgr), intent(in) :: dvxcdgr
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: wb_vxc
  !Local variables
  integer :: i1,i2,i3,dir_i
  real(kind=8) :: dnexcdgog,grad_i
  real(kind=8), dimension(:,:,:,:), allocatable :: f_i

  !Body

  allocate(f_i(n01,n02,n03,3))
  
  !let us first treat the case nspden=1
  if (nspden == 1) then
     !Let us construct the object we have to manipulate with another gradient
     if (nvxcdgr == 3) then
        do dir_i=1,3
           !Let us construct the object we have to manipulate with another gradient
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    dnexcdgog=0.5d0*dvxcdgr(i1,i2,i3,1) + dvxcdgr(i1,i2,i3,3)
                    grad_i=2.d0*gradient(i1,i2,i3,1,dir_i)
                    f_i(i1,i2,i3,dir_i)=dnexcdgog*grad_i
                 end do
              end do
           end do
        end do
     else
        do dir_i=1,3
           !Let us construct the object we have to manipulate with another gradient
           do i3=1,n03
              do i2=1,n02
                 do i1=1,n01
                    dnexcdgog=0.5d0*dvxcdgr(i1,i2,i3,1)
                    grad_i=2.d0*gradient(i1,i2,i3,1,dir_i)
                    f_i(i1,i2,i3,dir_i)=dnexcdgog*grad_i
                 end do
              end do
           end do
        end do
     end if
     !let us now calculate the gradient and correct the result
     call wb_correction(n01,n02,n03,n3eff,wbl,wbr,f_i,&
          hx,hy,hz,nspden,wb_vxc)

  !then the spin-polarized case
  else

     !!to be inserted later, when the non spin-pol case is verified

  !end of spin-polarized if statement
  end if

end subroutine vxcpostprocessing


!!$!!****f* ABINIT/size_dvxc
!!$!! NAME
!!$!! size_dvxc
!!$!!
!!$!! FUNCTION
!!$!! Give the size of the array dvxc(npts,ndvxc) 
!!$!! needed for the allocations depending on the routine which is called from the drivexc routine
!!$!!
!!$!! COPYRIGHT
!!$!! Copyright (C) 1998-2006 ABINIT group (TD)
!!$!! This file is distributed under the terms of the
!!$!! GNU General Public License, see ~abinit/COPYING
!!$!! or http://www.gnu.org/copyleft/gpl.txt .
!!$!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!$!! This routine has been written from rhohxc_coll(DCA, XG, GMR, MF, GZ)
!!$!!
!!$!! INPUTS
!!$!!  ixc= choice of exchange-correlation scheme
!!$!!  order=gives the maximal derivative of Exc computed.
!!$!!    1=usual value (return exc and vxc)
!!$!!    2=also computes the kernel (return exc,vxc,kxc)
!!$!!   -2=like 2, except (to be described)
!!$!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!$!!
!!$!! OUTPUT
!!$!!  ndvxc size of the array dvxc(npts,ndvxc) for allocation
!!$!!
!!$!! PARENTS
!!$!!      pawxc,pawxcm,rhohxc_coll
!!$!!
!!$!! CHILDREN
!!$!!      
!!$!!
!!$!! SOURCE
!!$
!!$subroutine size_dvxc(ixc,ndvxc,nspden,nvxcdgr,order)
!!$
!!$ implicit none
!!$
!!$!Arguments----------------------
!!$ integer, intent(in) :: ixc,nspden,order
!!$ integer, intent(out) :: ndvxc,nvxcdgr
!!$
!!$!Local variables----------------
!!$ nvxcdgr=0
!!$ ndvxc=0
!!$ if (order**2 <= 1) then
!!$    ndvxc=0
!!$    nvxcdgr=0
!!$    if (ixc>=11 .and. ixc<=15 .and. ixc/=13) nvxcdgr=3 
!!$    if (ixc==16) nvxcdgr=2
!!$ else
!!$    if (ixc==1 .or. ixc==21 .or. ixc==22 .or. (ixc>=7 .and. ixc<=10) .or. ixc==13) then
!!$       ! new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option    !routine xcspol
!!$       !routine xcpbe, with different options (optpbe) and orders (order)
!!$       ndvxc=nspden+1
!!$       !if (ixc>=7 .and. ixc<=10) nvxcdgr=3
!!$    else if (ixc>=2 .and. ixc<=6) then
!!$       ! Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)                !routine xcpzca
!!$       ! Teter fit (4/91) to Ceperley-Alder values (no spin-pol)              !routine xctetr
!!$       ! Wigner xc (no spin-pol)           !routine xcwign
!!$       ! Hedin-Lundqvist xc (no spin-pol)           !routine xchelu
!!$       ! X-alpha (no spin-pol)           !routine xcxalp
!!$       ndvxc=1
!!$      
!!$    else if (ixc==12) then
!!$       !routine xcpbe, with optpbe=-2 and different orders (order)
!!$       ndvxc=8
!!$       nvxcdgr=3
!!$    else if (ixc>=11 .and. ixc<=15 .and. ixc/=13) then
!!$       !routine xcpbe, with different options (optpbe) and orders (order)
!!$       ndvxc=15
!!$       nvxcdgr=3
!!$    else if(ixc==16 .or. ((ixc>=30).and.(ixc<=34)) ) then
!!$       !Should be 0
!!$       ndvxc=0
!!$       if (ixc==16) nvxcdgr=2
!!$    end if
!!$ end if
!!$
!!$end subroutine size_dvxc


!comparison subroutine
subroutine compare(npts,vect1,vect2,sum)
  implicit none
  !Arguments ------------------------------------
  integer, intent(in) :: npts
  real(kind=8), dimension(npts), intent(in) :: vect1,vect2
  real(kind=8), intent(out) :: sum
  !Local variables
  integer :: i
  real(kind=8) :: check

  sum=0.d0
  do i=1,npts
     check=abs(vect1(i)-vect2(i))
     sum=max(sum,check)
  end do

end subroutine compare

!fake subroutines
subroutine wrtout(unit,message,mode_paral)
  implicit none

  !Arguments ------------------------------------
  integer,intent(in) :: unit
  character(len=4),intent(in) :: mode_paral
  character(len=500),intent(inout) :: message

  print *,message
end subroutine wrtout

subroutine leave_new(mode_paral)

  implicit none

  !Arguments ------------------------------------
  character(len=4),intent(in) :: mode_paral

  print *,'exiting...'
  stop
end subroutine leave_new

