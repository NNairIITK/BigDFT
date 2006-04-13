

        subroutine input_rho_ion(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,n1,n2,n3,hgrid,rho,eion)
! Creates charge density arising from the ionoc PSP cores
        implicit real*8 (a-h,o-z)
        character*20 :: atomnames(100)
        dimension psppar(0:2,0:4,ntypes),rxyz(3,nat),iatype(nat),nelpsp(ntypes)
        dimension rho(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 
        pi=4.d0*atan(1.d0)
	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

! Ionic charge 
   rholeaked=0.d0
   eion=0.d0
   do iat=1,nat
   ityp=iatype(iat)
   rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
   ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
!    ion-ion interaction
     do jat=1,iat-1
     dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
     jtyp=iatype(jat)
     eion=eion+nelpsp(jtyp)*nelpsp(ityp)/dist
     enddo

     rloc=psppar(0,0,ityp)
     if (iproc.eq.0) write(*,'(a,i3,a,a,a,e10.3)') 'atom ',iat,' of type ',atomnames(ityp),' has an ionic charge with rloc',rloc
     charge=nelpsp(ityp)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

      do i3=iz-ii,iz+ii
      do i2=iy-ii,iy+ii
      do i1=ix-ii,ix+ii
         x=i1*hgridh-rx
         y=i2*hgridh-ry
         z=i3*hgridh-rz
         r2=x**2+y**2+z**2
         arg=r2/rloc**2
        xp=exp(-.5d0*arg)
        if (i3.ge.-14 .and. i3.le.2*n3+16  .and.  & 
            i2.ge.-14 .and. i2.le.2*n2+16  .and.  & 
            i1.ge.-14 .and. i1.le.2*n1+16 ) then
        rho(i1,i2,i3)=rho(i1,i2,i3)-xp*charge
        else
        rholeaked=rholeaked+xp*charge
        endif
      enddo
      enddo
      enddo

    enddo

! Check
	tt=0.d0
        do i3= -14,2*n3+16
        do i2= -14,2*n2+16
        do i1= -14,2*n1+16
        tt=tt+rho(i1,i2,i3)
        enddo
        enddo
        enddo
        tt=tt*hgridh**3
        rholeaked=rholeaked*hgridh**3
	if (iproc.eq.0) write(*,'(a,e21.14,1x,e10.3)') 'total ionic charge,leaked charge: ',tt,rholeaked

        return
	end


        subroutine addlocgauspsp(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,n1,n2,n3,hgrid,pot)
! Add local Gaussian terms of the PSP to pot 
        implicit real*8 (a-h,o-z)
        character*20 :: atomnames(100)
        dimension psppar(0:2,0:4,ntypes),rxyz(3,nat),iatype(nat)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 

   do iat=1,nat
   ityp=iatype(iat)
   rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
   ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
! determine number of local terms
     nloc=0
     do iloc=1,4
     if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo

     rloc=psppar(0,0,ityp)
     if (iproc.eq.0) write(*,'(a,i4,a,a,a,i3,a,1x,e9.2)')  & 
     'atom ',iat,' is of type ',atomnames(ityp),' and has ',nloc,' local terms with rloc',rloc
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

      do i3=max(-14,iz-ii),min(2*n3+16,iz+ii)
      do i2=max(-14,iy-ii),min(2*n2+16,iy+ii)
      do i1=max(-14,ix-ii),min(2*n1+16,ix+ii)
         x=i1*hgridh-rx
         y=i2*hgridh-ry
         z=i3*hgridh-rz
         r2=x**2+y**2+z**2
         arg=r2/rloc**2
        xp=exp(-.5d0*arg)
        tt=psppar(0,nloc,ityp)
	do iloc=nloc-1,1,-1
        tt=arg*tt+psppar(0,iloc,ityp)
        enddo
        pot(i1,i2,i3)=pot(i1,i2,i3)+xp*tt
      enddo
      enddo
      enddo

!! For testing only: Add erf part (in the final version that should be part of Hartree pot)
!      do i3=-14,2*n3+16
!      do i2=-14,2*n2+16
!      do i1=-14,2*n1+16
!         x=(i1-ix)*hgridh
!         y=(i2-iy)*hgridh
!         z=(i3-iz)*hgridh
!         r2=x**2+y**2+z**2
!         r=sqrt(r2)
!         arg=r*(sqrt(.5d0)/rloc)
!         if (arg.lt.1.d-7) then 
!! Taylor expansion
!         x=arg**2
!         tt=   -0.37612638903183752463d0*x + 1.1283791670955125739d0
!         tt=tt*(sqrt(.5d0)/rloc)
!         else
!          tt=derf(arg)/r
!         endif
!        pot(i1,i2,i3)=pot(i1,i2,i3)+nelpsp(ityp)*tt
!      enddo
!      enddo
!      enddo


    enddo

        return
	end
