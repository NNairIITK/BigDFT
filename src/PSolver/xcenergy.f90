subroutine xc_energy(geocode,m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
     nxcl,nxcr,ixc,hx,hy,hz,rhopot,pot_ion,zf,zfionxc,exc,vxc,iproc,nproc)

  implicit none

  !Arguments----------------------
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: m1,m2,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,ixc,iproc,nproc
  integer, intent(in) :: nwbl,nwbr
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(m1,m3,nxt), intent(inout) :: rhopot
  real(kind=8), dimension(m1,m3,nxc), intent(in) :: pot_ion
  real(kind=8), dimension(md1,md3,md2/nproc), intent(out) :: zf,zfionxc
  real(kind=8), intent(out) :: exc,vxc

  !Local variables----------------
  real(kind=8), dimension(:,:,:), allocatable :: exci,d2vxci
  real(kind=8), dimension(:,:,:,:), allocatable :: vxci,dvxci,dvxcdgr
  real(kind=8), dimension(:,:,:,:,:), allocatable :: gradient
  real(kind=8) :: elocal,vlocal,rho,pot,potion,factor,hgrid,facpotion
  integer :: npts,i_all,nspden,order,offset,i_stat
  integer :: i1,i2,i3,j1,j2,j3,jp2,jpp2,jppp2
  integer :: ndvxc,nvxcdgr,ngr2

  !interface with drivexc
  interface
     subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr, & !Mandatory 
          & dvxc,d2vxc,grho2_updn,vxcgr) !Optional
       use defs_basis
       integer,intent(in) :: ixc,npts,nspden,order
       integer,intent(in) :: ndvxc,ngr2,nvxcdgr
       real(dp),intent(in) :: rho_updn(npts,nspden)
       real(dp),intent(in), optional :: grho2_updn(npts,ngr2)
       real(dp),intent(out) :: exc(npts),vxc(npts,nspden)
       real(dp),intent(out), optional :: d2vxc(npts),dvxc(npts,ndvxc),vxcgr(npts,nvxcdgr)
     end subroutine drivexc
  end interface

  !Body

  !check for the dimensions
  if (  nwb/=nxcl+nxc+nxcr-2 .or. nxt/=nwbr+nwb+nwbl) then
     print *,'the XC dimensions are not correct'
     print *,'nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr',nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
     stop
  end if

  !these are always the same
  nspden=1
  order=1
  
  !useful for the freeBC case
  hgrid=max(hx,hy,hz)

  !Allocations
  allocate(exci(m1,m3,nwb),stat=i_stat)
  call memocc(i_stat,product(shape(exci))*kind(exci),'exci','xc_energy')
  allocate(vxci(m1,m3,nwb,nspden),stat=i_stat)
  call memocc(i_stat,product(shape(vxci))*kind(vxci),'vxci','xc_energy')
  !Allocations of the exchange-correlation terms, depending on the ixc value
  if (ixc >= 11 .and. ixc <= 16) then
     allocate(gradient(m1,m3,nwb,2*nspden-1,0:3),stat=i_stat)
     call memocc(i_stat,product(shape(gradient))*kind(gradient),'gradient','xc_energy')
  end if

  call size_dvxc(ixc,ndvxc,ngr2,nspden,nvxcdgr,order)

  if (ndvxc/=0) then
     allocate(dvxci(m1,m3,nwb,ndvxc),stat=i_stat)
     call memocc(i_stat,product(shape(dvxci))*kind(dvxci),'dvxci','xc_energy')
  end if
  if (nvxcdgr/=0) then
     allocate(dvxcdgr(m1,m3,nwb,nvxcdgr),stat=i_stat)
     call memocc(i_stat,product(shape(dvxcdgr))*kind(dvxcdgr),'dvxcdgr','xc_energy')
  end if
  if ((ixc==3 .or. (ixc>=7 .and. ixc<=15)) .and. order==3) then
     allocate(d2vxci(m1,m3,nwb),stat=i_stat)
     call memocc(i_stat,product(shape(d2vxci))*kind(d2vxci),'d2vxci','xc_energy')
  end if

  if (.not.allocated(gradient) .and. nxc/=nxt ) then
     print *,'Parxc_energy: if nx2/=nxc the gradient must be allocated'
     stop
  end if

 !divide by two the density to applicate it in the ABINIT xc routines
  offset=nwbl+1
  if (ixc/=0) then
     !divide by two the density to applicate it in the ABINIT xc routines
     do i3=1,nxt
        do i2=1,m3
           do i1=1,m1
              rhopot(i1,i2,i3)=.5d0*rhopot(i1,i2,i3)
           end do
        end do
     end do

     !computation of the gradient

     if (allocated(gradient)) then
        
        !the calculation of the gradient will depend on the geometry code
        if (geocode=='F') then
           call calc_gradient(m1,m3,nxt,nwb,nwbl,nwbr,rhopot,nspden,&
                hgrid,hgrid,hgrid,gradient)
        else
        print *,'geocode=',geocode,&
             ':the calculation of the gradient is still to be performed in this case'
        stop
        end if

     end if

     !this part can be commented out if you don't want to use ABINIT modules
     npts=m1*m3*nwb
     !let us apply ABINIT routines
     !case with gradient
     if (ixc >= 11 .and. ixc <= 16) then
        if (order**2 <= 1 .or. ixc == 16) then
           if (ixc /= 13) then             
             call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                   &grho2_updn=gradient,vxcgr=dvxcdgr) 
           else
             call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                   &grho2_updn=gradient) 
           end if
        else if (order /= 3) then
           if (ixc /= 13) then             
             call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                   &dvxc=dvxci,grho2_updn=gradient,vxcgr=dvxcdgr) 
           else
             call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                   &dvxc=dvxci,grho2_updn=gradient) 
           end if
        else if (order == 3) then
           if (ixc /= 13) then             
             call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                   &dvxci,d2vxci,gradient,dvxcdgr) 
           else
             call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                   &dvxc=dvxci,d2vxc=d2vxci,grho2_updn=gradient) 
           end if
        end if

     !do not calculate the White-Bird term in the Leeuwen Baerends XC case
        if (ixc/=13 .and. geocode == 'F') then
           call vxcpostprocessing(m1,m3,nwb,nxc,nxcl,nxcr,nspden,nvxcdgr,gradient,&
                hgrid,hgrid,hgrid,dvxcdgr,vxci)
        end if

        !cases without gradient
     else
        if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr)
        else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                &dvxc=dvxci,d2vxc=d2vxci)
        else
           call drivexc(exci,ixc,npts,nspden,order,rhopot(1,1,offset),vxci,ndvxc,ngr2,nvxcdgr,&
                &dvxc=dvxci)
        end if
     end if
     !end of the part that can be commented out

  end if

  if (allocated(dvxci)) then
     i_all=-product(shape(dvxci))*kind(dvxci)
     deallocate(dvxci,stat=i_stat)
     call memocc(i_stat,i_all,'dvxci','xc_energy')
  end if
  if (allocated(dvxcdgr)) then
     i_all=-product(shape(dvxcdgr))*kind(dvxcdgr)
     deallocate(dvxcdgr,stat=i_stat)
     call memocc(i_stat,i_all,'dvxcdgr','xc_energy')
  end if
  if (allocated(d2vxci)) then
     i_all=-product(shape(d2vxci))*kind(d2vxci)
     deallocate(d2vxci,stat=i_stat)
     call memocc(i_stat,i_all,'d2vxci','xc_energy')
  end if
  if (allocated(gradient)) then
     i_all=-product(shape(gradient))*kind(gradient)
     deallocate(gradient,stat=i_stat)
     call memocc(i_stat,i_all,'gradient','xc_energy')
  end if

  !distributing the density in the zf array
  !calculating the xc integrated quantities
  !and summing the xc potential into the zfionxc array

  if(ixc==0) then
     factor=1.d0
     facpotion=0.d0
  else
     factor=2.d0
     facpotion=1.d0
  end if


  exc=0.d0
  vxc=0.d0
  do jp2=1,nxc
     j2=offset+jp2+nxcl-2
     jpp2=jp2     
     jppp2=jp2+nxcl-1
     do j3=1,m3
        do j1=1,m1
           rho=rhopot(j1,j3,j2)
           potion=pot_ion(j1,j3,jpp2)
           if (rho < 5.d-21 .or. ixc == 0) then
              elocal=0.d0
              vlocal=0.d0
           else
              elocal=exci(j1,j3,jppp2)
              vlocal=vxci(j1,j3,jppp2,1)
           end if
           exc=exc+elocal*rho
           vxc=vxc+vlocal*rho
           zf(j1,j3,jp2)=factor*rhopot(j1,j3,j2)
           zfionxc(j1,j3,jp2)=facpotion*potion+vlocal
        end do
        do j1=m1+1,md1
           zf(j1,j3,jp2)=0.d0
           zfionxc(j1,j3,jp2)=0.d0
        end do
     end do
     do j3=m3+1,md3
        do j1=1,md1
           zf(j1,j3,jp2)=0.d0
           zfionxc(j1,j3,jp2)=0.d0
        end do
     end do
  end do
  do jp2=nxc+1,md2/nproc
     do j3=1,md3
        do j1=1,md1
           zf(j1,j3,jp2)=0.d0
           zfionxc(j1,j3,jp2)=0.d0
        end do
     end do
  end do

  !the two factor is due to the 
  !need of using the density of states in abinit routines
  if (geocode== 'F') then
     exc=2.d0*hgrid**3*exc
     vxc=2.d0*hgrid**3*vxc
  else
     exc=2.d0*hx*hy*hz*exc
     vxc=2.d0*hx*hy*hz*vxc
  end if

  !De-allocations
  i_all=-product(shape(exci))*kind(exci)
  deallocate(exci,stat=i_stat)
  call memocc(i_stat,i_all,'exci','xc_energy')
  i_all=-product(shape(vxci))*kind(vxci)
  deallocate(vxci,stat=i_stat)
  call memocc(i_stat,i_all,'vxci','xc_energy')

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
!! SOURCE
subroutine vxcpostprocessing(n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr,gradient,hx,hy,hz,dvxcdgr,wb_vxc)
  implicit none
  integer, intent(in) :: n01,n02,n03,n3eff,wbl,wbr,nspden,nvxcdgr
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), dimension(n01,n02,n03,2*nspden-1,0:3), intent(in) :: gradient
  real(kind=8), dimension(n01,n02,n03,nvxcdgr), intent(in) :: dvxcdgr
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: wb_vxc
  !Local variables
  integer :: i1,i2,i3,dir_i,i_all,i_stat
  real(kind=8) :: dnexcdgog,grad_i
  real(kind=8), dimension(:,:,:,:), allocatable :: f_i

  !Body

  allocate(f_i(n01,n02,n03,3),stat=i_stat)
  call memocc(i_stat,product(shape(f_i))*kind(f_i),'f_i','vxcpostprocessing')
  
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

  i_all=-product(shape(f_i))*kind(f_i)
  deallocate(f_i,stat=i_stat)
  call memocc(i_stat,i_all,'f_i','vxcpostprocessing')

end subroutine vxcpostprocessing
