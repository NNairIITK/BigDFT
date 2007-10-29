subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
  ! Calculates the length of the keys describing a wavefunction data structure
  implicit real(kind=8) (a-h,o-z)
  logical logrid,plogrid
  dimension logrid(0:n1,0:n2,0:n3)

  mvctr=0
  nsrt=0
  nend=0
  do i3=nl3,nu3 
     do i2=nl2,nu2
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
     enddo
  enddo
  if (nend.ne.nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif
  mseg=nend
  
end subroutine num_segkeys

subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
  ! Calculates the keys describing a wavefunction data structure
  implicit real(kind=8) (a-h,o-z)
  logical logrid,plogrid
  dimension logrid(0:n1,0:n2,0:n3),keyg(2,mseg),keyv(mseg)

  mvctr=0
  nsrt=0
  nend=0
  do i3=nl3,nu3 
     do i2=nl2,nu2
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
  implicit real(kind=8) (a-h,o-z)
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
     rad=radii(iatype(iat))*rmult+real(nbuf,kind=8)*hgrid
     !        write(*,*) 'iat,nat,rad',iat,nat,rad
     ml1=ceiling((rxyz(1,iat)-rad)/hgrid - eps_mach)  
     ml2=ceiling((rxyz(2,iat)-rad)/hgrid - eps_mach)   
     ml3=ceiling((rxyz(3,iat)-rad)/hgrid - eps_mach)   
     mu1=floor((rxyz(1,iat)+rad)/hgrid + eps_mach)
     mu2=floor((rxyz(2,iat)+rad)/hgrid + eps_mach)
     mu3=floor((rxyz(3,iat)+rad)/hgrid + eps_mach)
     if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
     if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
     if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
     do i3=ml3,mu3
        dz2=(real(i3,kind=8)*hgrid-rxyz(3,iat))**2
        do i2=ml2,mu2
           dy2=(real(i2,kind=8)*hgrid-rxyz(2,iat))**2
           do i1=ml1,mu1
              dx=real(i1,kind=8)*hgrid-rxyz(1,iat)
              if (dx**2+(dy2+dz2).le.rad**2) then 
                 logrid(i1,i2,i3)=.true.
              endif
           enddo
        enddo
     enddo
  enddo

  return
END SUBROUTINE fill_logrid

subroutine make_bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
  implicit real(kind=8) (a-h,o-z)
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
END SUBROUTINE make_bounds
