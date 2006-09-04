
        subroutine system_size(nat,rxyz,radii,rmult,iatype,ntypes, &
                   cxmin,cxmax,cymin,cymax,czmin,czmax)
! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),radii(ntypes),iatype(nat)

        cxmax=-1.d100 ; cxmin=1.d100
        cymax=-1.d100 ; cymin=1.d100
        czmax=-1.d100 ; czmin=1.d100
        do iat=1,nat
            rad=radii(iatype(iat))*rmult
            cxmax=max(cxmax,rxyz(1,iat)+rad) ; cxmin=min(cxmin,rxyz(1,iat)-rad)
            cymax=max(cymax,rxyz(2,iat)+rad) ; cymin=min(cymin,rxyz(2,iat)-rad)
            czmax=max(czmax,rxyz(3,iat)+rad) ; czmin=min(czmin,rxyz(3,iat)-rad)
        enddo
  
      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
      czmax=czmax-eps_mach ; czmin=czmin+eps_mach

        return
        end subroutine system_size


	subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
! Calculates the length of the keys describing a wavefunction data structure
        implicit real*8 (a-h,o-z)
        logical logrid,plogrid
        dimension logrid(0:n1,0:n2,0:n3)

        mvctr=0
        nsrt=0
        nend=0
        do i3=nl3,nu3 ; do i2=nl2,nu2

        plogrid=.false.
        do i1=nl1,nu1
         if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
             nsrt=nsrt+1
           endif
         else
           if (plogrid .eqv. .true.) then
             nend=nend+1
           endif
         endif
         plogrid=logrid(i1,i2,i3)
        enddo 
           if (plogrid .eqv. .true.) then
             nend=nend+1
           endif
        enddo ; enddo
        if (nend.ne.nsrt) then 
           write(*,*) 'nend , nsrt',nend,nsrt
           stop 'nend <> nsrt'
        endif
        mseg=nend

	return
        end subroutine num_segkeys


	subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
! Calculates the keys describing a wavefunction data structure
        implicit real*8 (a-h,o-z)
        logical logrid,plogrid
        dimension logrid(0:n1,0:n2,0:n3),keyg(2,mseg),keyv(mseg)

        mvctr=0
        nsrt=0
        nend=0
        do i3=nl3,nu3 ; do i2=nl2,nu2

        plogrid=.false.
        do i1=nl1,nu1
         ngridp=i3*((n1+1)*(n2+1)) + i2*(n1+1) + i1+1
         if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
             nsrt=nsrt+1
             keyg(1,nsrt)=ngridp
             keyv(nsrt)=mvctr
           endif
         else
           if (plogrid .eqv. .true.) then
             nend=nend+1
             keyg(2,nend)=ngridp-1
           endif
         endif
         plogrid=logrid(i1,i2,i3)
        enddo 
           if (plogrid .eqv. .true.) then
             nend=nend+1
             keyg(2,nend)=ngridp
           endif
        enddo ; enddo
        if (nend.ne.nsrt) then 
           write(*,*) 'nend , nsrt',nend,nsrt
           stop 'nend <> nsrt'
        endif
        mseg=nend

	return
        end subroutine segkeys




       subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                               ntypes,iatype,rxyz,radii,rmult,hgrid,loregion,logrid)
! set up an array logrid(i1,i2,i3) that specifies whether the grid point
! i1,i2,i3 is the center of a scaling function/wavelet
        implicit real*8 (a-h,o-z)
        logical logrid, loregion
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),iatype(nat),radii(ntypes)
        dimension logrid(0:n1,0:n2,0:n3),loregion(nat)

        do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
         logrid(i1,i2,i3)=.false.
        enddo ; enddo ; enddo

      do iat=1,nat
      if (loregion(iat)) then
        rad=radii(iatype(iat))*rmult
!        write(*,*) 'iat,nat,rad',iat,nat,rad
        onem=1.d0-eps_mach
        ml1=int(onem+(rxyz(1,iat)-rad)/hgrid)  ; mu1=int((rxyz(1,iat)+rad)/hgrid)
        ml2=int(onem+(rxyz(2,iat)-rad)/hgrid)  ; mu2=int((rxyz(2,iat)+rad)/hgrid)
        ml3=int(onem+(rxyz(3,iat)-rad)/hgrid)  ; mu3=int((rxyz(3,iat)+rad)/hgrid)
!        write(*,'(a,6(i4))') 'fill grid',ml1,mu1,ml2,mu2,ml3,mu3
        if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
        if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
        if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
        do i3=ml3,mu3
        dz2=(i3*hgrid-rxyz(3,iat))**2
        do i2=ml2,mu2
        dy2=(i2*hgrid-rxyz(2,iat))**2
        do i1=ml1,mu1
        dx=i1*hgrid-rxyz(1,iat)
        if (dx**2+(dy2+dz2).lt.rad**2) logrid(i1,i2,i3)=.true.
        enddo ; enddo ; enddo
      endif
      enddo

        return
        end subroutine fill_logrid
