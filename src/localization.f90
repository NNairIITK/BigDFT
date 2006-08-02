

        subroutine wannier_par(iproc,nat,norb,rxyz,ntypes,iatype,rad_cov,wan_par)
        implicit real*8 (a-h,o-z)
        dimension rxyz(3,nat),iatype(nat),rad_cov(ntypes),wan_par(0:3,norb)
! Wannier centers are in between any pair of atoms whose distance is less than the sum of the covalent radii

        iorb=0
        do iat=1,nat
        ityp=iatype(iat)
        do jat=1,iat-1
        jtyp=iatype(jat)
        dist2=(rxyz(1,iat)-rxyz(1,jat))**2 + (rxyz(2,iat)-rxyz(2,jat))**2+ &
              (rxyz(3,iat)-rxyz(3,jat))**2
! add a safety margin of 20 percent to the covalent radii
        rcut2=(1.2d0*(rad_cov(ityp)+rad_cov(jtyp)))**2

        if (dist2.lt.rcut2) then
           iorb=iorb+1
           if (iorb.gt.norb) stop 'wannier_par, iorb > norb'
           wan_par(1,iorb)=.5d0*(rxyz(1,iat)+rxyz(1,jat))
           wan_par(2,iorb)=.5d0*(rxyz(2,iat)+rxyz(2,jat))
           wan_par(3,iorb)=.5d0*(rxyz(3,iat)+rxyz(3,jat))
           wan_par(0,iorb)=max(rad_cov(ityp),rad_cov(jtyp))
!           if (iorb.eq.norb) goto 1232
        endif
        enddo
        enddo

!1232    continue
        if (iproc.eq.0) then
        write(*,*) ' Wannier centers and their radii'
        do jorb=1,iorb
        write(*,'(i4,4(1x,e10.3))')  &
             jorb,wan_par(1,jorb),wan_par(2,jorb),wan_par(3,jorb),wan_par(0,jorb)
        enddo
        endif
        if (iorb.ne.norb) stop 'wannier_par, iorb <> norb'

        return
        end subroutine wannier_par



        subroutine localizationregion(iproc,nat,norb,rxyz,wan_par,radlocmult,loregion)
! determines the localization region for each orbital in terms of the array loregion
! if loregion(iat,iorb).eq. .true. the atomic centered sphere of atom iat is part 
! of the localization region of orbital iorb
        implicit real*8 (a-h,o-z)
        logical loregion
        dimension rxyz(3,nat),wan_par(0:3,norb),loregion(nat,norb)

        do iorb=1,norb
          ic=0
          do iat=1,nat
            dd=(rxyz(1,iat)-wan_par(1,iorb))**2+(rxyz(2,iat)-wan_par(2,iorb))**2  &
                                               +(rxyz(3,iat)-wan_par(3,iorb))**2 
            if (dd.le.(radlocmult*wan_par(0,iorb))**2) then
              ic=ic+1
              loregion(iat,iorb)=.true.
            else
              loregion(iat,iorb)=.false.
            endif
          enddo
          if (iproc.eq.0) write(*,'(a,i4,a,i3,a)')   & 
             'localization region of orbital ',iorb,' consists of ',ic,' atom centered spheres'
          if (ic.lt.1) stop 'vanishing localization region'
        enddo

        return
        end subroutine localizationregion


        subroutine loregion_size(nat,rxyz,radii,rmult,iatype,ntypes,loregion, &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
! finds the size of the smallest subbox that contains a localization region made 
! out of atom centered spheres
        implicit real*8 (a-h,o-z)
        logical loregion
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),radii(ntypes),iatype(nat),loregion(nat)

        cxmax=-1.d100 ; cxmin=1.d100
        cymax=-1.d100 ; cymin=1.d100
        czmax=-1.d100 ; czmin=1.d100
        do iat=1,nat
        if (loregion(iat)) then
            rad=radii(iatype(iat))*rmult
            cxmax=max(cxmax,rxyz(1,iat)+rad) ; cxmin=min(cxmin,rxyz(1,iat)-rad)
            cymax=max(cymax,rxyz(2,iat)+rad) ; cymin=min(cymin,rxyz(2,iat)-rad)
            czmax=max(czmax,rxyz(3,iat)+rad) ; czmin=min(czmin,rxyz(3,iat)-rad)
!        write(*,*) radii(iatype(iat)),rmult
!        write(*,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
!        write(*,*) 'loregion_size',cxmin,cxmax
!        write(*,*) '             ',cymin,cymax
!        write(*,*) '             ',czmin,czmax
        endif
        enddo
  
      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
      czmax=czmax-eps_mach ; czmin=czmin+eps_mach
      onem=1.d0-eps_mach
      nl1=int(onem+cxmin/hgrid)   
      nl2=int(onem+cymin/hgrid)   
      nl3=int(onem+czmin/hgrid)   
      nu1=int(cxmax/hgrid)  
      nu2=int(cymax/hgrid)  
      nu3=int(czmax/hgrid)  
!        write(*,'(a,6(i4))') 'loregion_size',nl1,nu1,nl2,nu2,nl3,nu3
!        write(*,*) 'loregion_size',cxmin,cxmax
!        write(*,*) '             ',cymin,cymax
!        write(*,*) '             ',czmin,czmax
      if (nl1.lt.0)   stop 'nl1: localization region outside cell'
      if (nl2.lt.0)   stop 'nl2: localization region outside cell'
      if (nl3.lt.0)   stop 'nl3: localization region outside cell'
      if (nu1.gt.n1)   stop 'nu1: localization region outside cell'
      if (nu2.gt.n2)   stop 'nu2: localization region outside cell'
      if (nu3.gt.n3)   stop 'nu3: localization region outside cell'

        return
        end subroutine loregion_size
