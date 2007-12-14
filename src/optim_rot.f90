
  implicit real*8 (a-h,o-z)
  character*20 atomnames
  parameter(rad=4.d0) ! this rad is only for visualization
  dimension rxyz(3,10000),atomnames(10000),urot(3,3)
  dimension txyz(3,10000)

  open(unit=1,file='posinp')
  read(1,*) nat
  if (nat.gt.10000) stop 'too many atoms'
  do iat=1,nat
     read(1,*) rxyz(:,iat),atomnames(iat)
  enddo
  close(1)

  !! translate molecule such that the origin of the coordinate system is in its center
  !         xmin=1.d100 ; xmax=-1.d100
  !         ymin=1.d100 ; ymax=-1.d100
  !         zmin=1.d100 ; zmax=-1.d100
  !         do iat=1,nat
  !         xmin=min(xmin,rxyz(1,iat)) ; xmax=max(xmax,rxyz(1,iat)) 
  !         ymin=min(ymin,rxyz(2,iat)) ; ymax=max(ymax,rxyz(2,iat)) 
  !         zmin=min(zmin,rxyz(3,iat)) ; zmax=max(zmax,rxyz(3,iat)) 
  !         enddo
  !
  !         do iat=1,nat
  !         rxyz(1,iat)=rxyz(1,iat)-.5d0*(xmax+xmin)
  !         rxyz(2,iat)=rxyz(2,iat)-.5d0*(ymax+ymin)
  !         rxyz(3,iat)=rxyz(3,iat)-.5d0*(zmax+zmin)
  !         enddo

  call volume(nat,rxyz,vol)
  write(*,*) 'initial volume',vol

  it=0
  diag=1.d-2 ! initial small diagonal element allows for search over all angles
100 continue  ! loop over all trial rotations
  diag=diag*1.0001d0 ! increase diag to search over smaller angles
  it=it+1
  if (diag.gt.100.d0) stop ! smaller angle rotations do not make sense

  ! create a random orthogonal (rotation) matrix
  call random_number(urot)
  urot(:,:)=urot(:,:)-.5d0
  do i=1,3
     urot(i,i)=urot(i,i)+diag
  enddo

  s=urot(1,1)**2+urot(2,1)**2+urot(3,1)**2
  s=1.d0/sqrt(s)
  urot(:,1)=s*urot(:,1) 

  s=urot(1,1)*urot(1,2)+urot(2,1)*urot(2,2)+urot(3,1)*urot(3,2)
  urot(:,2)=urot(:,2)-s*urot(:,1)
  s=urot(1,2)**2+urot(2,2)**2+urot(3,2)**2
  s=1.d0/sqrt(s)
  urot(:,2)=s*urot(:,2) 

  s=urot(1,1)*urot(1,3)+urot(2,1)*urot(2,3)+urot(3,1)*urot(3,3)
  urot(:,3)=urot(:,3)-s*urot(:,1)
  s=urot(1,2)*urot(1,3)+urot(2,2)*urot(2,3)+urot(3,2)*urot(3,3)
  urot(:,3)=urot(:,3)-s*urot(:,2)
  s=urot(1,3)**2+urot(2,3)**2+urot(3,3)**2
  s=1.d0/sqrt(s)
  urot(:,3)=s*urot(:,3) 

  ! eliminate reflections
  if (urot(1,1).le.0.d0) urot(:,1)=-urot(:,1)
  if (urot(2,2).le.0.d0) urot(:,2)=-urot(:,2)
  if (urot(3,3).le.0.d0) urot(:,3)=-urot(:,3)

  ! apply the rotation to all atomic positions! 
  do iat=1,nat
     x=rxyz(1,iat) ; y=rxyz(2,iat) ; z=rxyz(3,iat)
     txyz(:,iat)=x*urot(:,1)+y*urot(:,2)+z*urot(:,3)
  enddo
  call volume(nat,txyz,tvol)
  if (tvol.lt.vol) then
     write(*,*) 'new best volume',tvol,it,diag
     rxyz(:,:)=txyz(:,:)
     vol=tvol
     ! write positions into output file
     xmin=1.d100 ; xmax=-1.d100
     ymin=1.d100 ; ymax=-1.d100
     zmin=1.d100 ; zmax=-1.d100
     do iat=1,nat
        xmin=min(xmin,rxyz(1,iat)) ; xmax=max(xmax,rxyz(1,iat)) 
        ymin=min(ymin,rxyz(2,iat)) ; ymax=max(ymax,rxyz(2,iat)) 
        zmin=min(zmin,rxyz(3,iat)) ; zmax=max(zmax,rxyz(3,iat)) 
     enddo
     ! Switch coordinates such that box is longest along z
     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     ! if box longest along x switch x and z
     if (xmax-xmin == dmax)  then
        do  iat=1,nat
           tx=rxyz(1,iat)
           tz=rxyz(3,iat)

           rxyz(1,iat)=tz
           rxyz(3,iat)=tx
        enddo
        tmin=xmin   
        xmin=zmin   
        zmin=tmin   
  
        tmax=xmax
        xmax=zmax
        zmax=tmax
     endif
     ! if box longest along y switch y and z
     if (ymax-ymin == dmax)  then
        do  iat=1,nat
           ty=rxyz(1,iat) ; tz=rxyz(3,iat)
           rxyz(1,iat)=tz ; rxyz(3,iat)=ty
        enddo
        tmin=ymin   
        ymin=zmin   
        zmin=tmin   

        tmax=ymax
        ymax=zmax
        zmax=tmax
     endif

     open(unit=2,file='rotpos.ascii')
     write(2,*) nat
     write(2,*) xmax-xmin+2*rad,' 0. ',ymax-ymin+2.d0*rad
     write(2,*) ' 0.  ',' 0. ',zmax-zmin+2.d0*rad
     do iat=1,nat
        write(2,'(3(1x,e17.10),1x,a20)') rxyz(1,iat)-xmin+rad,rxyz(2,iat)-ymin+rad,rxyz(3,iat)-zmin+rad,atomnames(iat)
     enddo
     close(2)
  endif
  goto 100
end program



subroutine volume(nat,rxyz,vol)
  implicit real*8 (a-h,o-z)
  parameter(rad=10.d0) ! this is the rad for the simulation box size
  dimension rxyz(3,nat)
  ! In the present version the radii around all spheres are taken to be 10 
  ! In a optimized version they shoud be replaced by the real radius, i.e. crmult*atomradius

  xmin=1.d100 ; xmax=-1.d100
  ymin=1.d100 ; ymax=-1.d100
  zmin=1.d100 ; zmax=-1.d100
  do iat=1,nat
     xmin=min(xmin,rxyz(1,iat)) ; xmax=max(xmax,rxyz(1,iat)) 
     ymin=min(ymin,rxyz(2,iat)) ; ymax=max(ymax,rxyz(2,iat)) 
     zmin=min(zmin,rxyz(3,iat)) ; zmax=max(zmax,rxyz(3,iat)) 
  enddo

  vol=(xmax-xmin+2.d0*rad)*(ymax-ymin+2.d0*rad)*(zmax-zmin+2.d0*rad)

end subroutine volume
