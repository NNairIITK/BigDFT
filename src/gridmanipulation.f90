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
END SUBROUTINE system_size




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
  enddo; enddo
  if (nend.ne.nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif
  mseg=nend

  return
END SUBROUTINE num_segkeys


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
  enddo; enddo
  if (nend.ne.nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif
  mseg=nend

  return
END SUBROUTINE segkeys

subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hgrid,logrid)
  ! set up an array logrid(i1,i2,i3) that specifies whether the grid point
  ! i1,i2,i3 is the center of a scaling function/wavelet
  implicit real*8 (a-h,o-z)
  logical logrid
  parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
  dimension rxyz(3,nat),iatype(nat),radii(ntypes)
  dimension logrid(0:n1,0:n2,0:n3)

  do i3=nl3,nu3 
     do i2=nl2,nu2 
        do i1=nl1,nu1
           logrid(i1,i2,i3)=.false.
        enddo
     enddo
  enddo

  do iat=1,nat
     rad=radii(iatype(iat))*rmult+nbuf*hgrid
     !        write(*,*) 'iat,nat,rad',iat,nat,rad
     ml1=int(onem+(rxyz(1,iat)-rad)/hgrid)  ; mu1=int((rxyz(1,iat)+rad)/hgrid)
     ml2=int(onem+(rxyz(2,iat)-rad)/hgrid)  ; mu2=int((rxyz(2,iat)+rad)/hgrid)
     ml3=int(onem+(rxyz(3,iat)-rad)/hgrid)  ; mu3=int((rxyz(3,iat)+rad)/hgrid)

!!$       print *,'values of the mls',ml1,nl1,ml2,nl2,ml3,nl3,nbuf,rmult,rad
!!$       print *,'values of the mus',mu1,nu1,mu2,nu2,mu3,nu3,nbuf,hgrid,rxyz(1,iat)

     if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
     if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
     if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
     do i3=ml3,mu3
        dz2=(i3*hgrid-rxyz(3,iat))**2
        do i2=ml2,mu2
           dy2=(i2*hgrid-rxyz(2,iat))**2
           do i1=ml1,mu1
              dx=i1*hgrid-rxyz(1,iat)
              if (dx**2+(dy2+dz2).le.rad**2) then 
                 logrid(i1,i2,i3)=.true.
              endif
           enddo
        enddo
     enddo
  enddo

  return
END SUBROUTINE fill_logrid

subroutine bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
  implicit real*8 (a-h,o-z)
  logical logrid
  dimension logrid(0:n1,0:n2,0:n3)
  dimension ibyz(2,0:n2,0:n3),ibxz(2,0:n1,0:n3),ibxy(2,0:n1,0:n2)


  do i3=0,n3 
     do i2=0,n2 
        ibyz(1,i2,i3)= 1000
        ibyz(2,i2,i3)=-1000

        do i1=0,n1
           if (logrid(i1,i2,i3)) then 
              ibyz(1,i2,i3)=i1
              goto 10
           endif
        enddo
10      continue
        do i1=n1,0,-1
           if (logrid(i1,i2,i3)) then 
              ibyz(2,i2,i3)=i1
              goto 11
           endif
        enddo
11      continue

     end do
  end do


  do i3=0,n3 
     do i1=0,n1
        ibxz(1,i1,i3)= 1000
        ibxz(2,i1,i3)=-1000

        do i2=0,n2 
           if (logrid(i1,i2,i3)) then 
              ibxz(1,i1,i3)=i2
              goto 20 
           endif
        enddo
20      continue
        do i2=n2,0,-1
           if (logrid(i1,i2,i3)) then 
              ibxz(2,i1,i3)=i2
              goto 21 
           endif
        enddo
21      continue

     end do
  end do


  do i2=0,n2 
     do i1=0,n1 
        ibxy(1,i1,i2)= 1000
        ibxy(2,i1,i2)=-1000

        do i3=0,n3
           if (logrid(i1,i2,i3)) then 
              ibxy(1,i1,i2)=i3
              goto 30 
           endif
        enddo
30      continue
        do i3=n3,0,-1
           if (logrid(i1,i2,i3)) then 
              ibxy(2,i1,i2)=i3
              goto 31 
           endif
        enddo
31      continue

     end do
  end do

  return
END SUBROUTINE bounds
