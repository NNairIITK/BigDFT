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
     write(*,*)' ERROR: nend <> nsrt',nend,nsrt
     stop 
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

subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
  ! set up an array logrid(i1,i2,i3) that specifies whether the grid point
  ! i1,i2,i3 is the center of a scaling function/wavelet
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
  real(kind=8), intent(in) :: rmult,hx,hy,hz
  integer, dimension(nat), intent(in) :: iatype
  real(kind=8), dimension(ntypes), intent(in) :: radii
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12,onem=1.d0-eps_mach
  integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
  real(kind=8) :: dx,dy2,dz2,rad

  !some checks
  if (geocode /='F') then
     !the nbuf value makes sense only in the case of free BC
     if (nbuf /=0) then
        write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
        stop
     end if
     !the grid spacings must be the same
     if (hx/= hy .or. hy /=hz .or. hx/=hz) then
        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
     end if
  end if

  if (geocode == 'F') then
     do i3=nl3,nu3 
        do i2=nl2,nu2 
           do i1=nl1,nu1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  else !
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  end if

  do iat=1,nat
     rad=radii(iatype(iat))*rmult+real(nbuf,kind=8)*hx
     !        write(*,*) 'iat,nat,rad',iat,nat,rad
     ml1=ceiling((rxyz(1,iat)-rad)/hx - eps_mach)  
     ml2=ceiling((rxyz(2,iat)-rad)/hy - eps_mach)   
     ml3=ceiling((rxyz(3,iat)-rad)/hz - eps_mach)   
     mu1=floor((rxyz(1,iat)+rad)/hx + eps_mach)
     mu2=floor((rxyz(2,iat)+rad)/hy + eps_mach)
     mu3=floor((rxyz(3,iat)+rad)/hz + eps_mach)
     !for Free BC, there must be no incoherences with the previously calculated delimiters
     if (geocode == 'F') then
        if (ml1.lt.nl1) stop 'ml1 < nl1'
        if (ml2.lt.nl2) stop 'ml2 < nl2'
        if (ml3.lt.nl3) stop 'ml3 < nl3'

        if (mu1.gt.nu1) stop 'mu1 > nu1'
        if (mu2.gt.nu2) stop 'mu2 > nu2'
        if (mu3.gt.nu3) stop 'mu3 > nu3'
     end if
     do i3=ml3,mu3
        dz2=(real(i3,kind=8)*hz-rxyz(3,iat))**2
        j3=modulo(i3,n3+1)
        do i2=ml2,mu2
           dy2=(real(i2,kind=8)*hy-rxyz(2,iat))**2
           j2=modulo(i2,n2+1)
           do i1=ml1,mu1
              j1=modulo(i1,n1+1)
              dx=real(i1,kind=8)*hx-rxyz(1,iat)
              if (dx**2+(dy2+dz2).le.rad**2) then 
                 logrid(j1,j2,j3)=.true.
              endif
           enddo
        enddo
     enddo
  enddo

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
