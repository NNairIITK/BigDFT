program optim_rot
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
end program optim_rot

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

  !rotate the molecule via an orthogonal matrix in order to minimise the
  !volume of the cubic cell
  subroutine optimise_volume(nat,ntypes,iatype,atomnames,crmult,frmult,hgrid,rxyz,radii_cf)
    use module_base
    implicit none
    integer, intent(in) :: nat,ntypes
    real(kind=8), intent(in) :: crmult,frmult
    real(kind=8), intent(inout) :: hgrid
    character(len=20), dimension(ntypes), intent(in) :: atomnames
    integer(kind=8), dimension(nat), intent(in) :: iatype
    real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
    real(kind=8), dimension(3,nat), intent(inout) :: rxyz
    !local variables
    character(len=*), parameter :: subname='optimise_volume'
    integer :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
    real(kind=8) :: x,y,z,vol,tx,ty,tz
    real(kind=8), dimension(3,3) :: urot
    real(kind=8), dimension(:,:), allocatable :: txyz

    allocate(txyz(3,nat+ndebug),stat=i_stat)
    call memocc(i_stat,txyz,'txyz',subname)
    
    call system_size(1,nat,ntypes,rxyz,radii_cf,crmult,frmult,hgrid,iatype,atomnames, &
         alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)
    !call volume(nat,rxyz,vol)
    vol=alat1*alat2*alat3
    write(*,'(1x,a,1pe16.8)')'Initial volume (Bohr^3)',vol

    it=0
    diag=1.d-2 ! initial small diagonal element allows for search over all angles
    loop_rotations: do  ! loop over all trial rotations
       diag=diag*1.0001d0 ! increase diag to search over smaller angles
       it=it+1
       if (diag.gt.100.d0) exit loop_rotations ! smaller angle rotations do not make sense

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

       call system_size(1,nat,ntypes,txyz,radii_cf,crmult,frmult,hgrid,iatype,atomnames, &
            alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)
       tvol=alat1*alat2*alat3
       !call volume(nat,txyz,tvol)
       if (tvol.lt.vol) then
          write(*,'(1x,a,1pe16.8,1x,i0,1x,f15.5)')'Found new best volume: ',tvol,it,diag
          rxyz(:,:)=txyz(:,:)
          vol=tvol
          dmax=max(alat1,alat2,alat3)
          ! if box longest along x switch x and z
          if (alat1 == dmax)  then
             do  iat=1,nat
                tx=rxyz(1,iat)
                tz=rxyz(3,iat)

                rxyz(1,iat)=tz
                rxyz(3,iat)=tx
             enddo
             ! if box longest along y switch y and z
          else if (alat2 == dmax)  then
             do  iat=1,nat
                ty=rxyz(1,iat) ; tz=rxyz(3,iat)
                rxyz(1,iat)=tz ; rxyz(3,iat)=ty
             enddo
          endif
       endif
    end do loop_rotations

    i_all=-product(shape(txyz))*kind(txyz)
    deallocate(txyz,stat=i_stat)
    call memocc(i_stat,i_all,'txyz',subname)
  end subroutine optimise_volume



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
