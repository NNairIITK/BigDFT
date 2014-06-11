module vector_operations
  contains
    function cross_product(a, b)
      implicit none
      real(kind=8),dimension(3) :: cross_product
      real(kind=8),dimension(3),intent(in) :: a, b
      cross_product(1) = a(2) * b(3) - a(3) * b(2)
      cross_product(2) = a(3) * b(1) - a(1) * b(3)
      cross_product(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product
end module vector_operations


module internal_coordinates

  contains

    !>calculates the dihedral angle between atoms i, j, k,
    !!            and l.  the cartesian coordinates of these atoms
    !!            are in array xyz.
    !!
    !!     dihed is a modified version of a subroutine of the same name
    !!           which was written by dr. w. theil in 1973.
    !!
    subroutine dihed(numat,xyz,i,j,k,l,angle)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: numat, i, j, k, l
      real(kind=8),dimension(3,numat),intent(in) :: xyz
      real(kind=8),intent(out) :: angle
    
      ! Local variables
      real(kind=8) :: xi1, xj1, xl1, yi1, yj1, yl1, zi1, zj1, zl1, dist, cosa, ddd, yxdist
      real(kind=8) :: xi2, xl2, yi2, yl2, xj2, yj2, costh, sinth, cosph, sinph, yi3, yl3
    
      xi1=xyz(1,i)-xyz(1,k)
      xj1=xyz(1,j)-xyz(1,k)
      xl1=xyz(1,l)-xyz(1,k)
      yi1=xyz(2,i)-xyz(2,k)
      yj1=xyz(2,j)-xyz(2,k)
      yl1=xyz(2,l)-xyz(2,k)
      zi1=xyz(3,i)-xyz(3,k)
      zj1=xyz(3,j)-xyz(3,k)
      zl1=xyz(3,l)-xyz(3,k)
      !      rotate around z axis to put kj along y axis
      dist= sqrt(xj1**2+yj1**2+zj1**2)
      cosa=zj1/dist
      if(cosa.gt.1.0d0) cosa=1.0d0
      if(cosa.lt.-1.0d0) cosa=-1.0d0
      ddd=1.0d0-cosa**2
      !if(ddd.le.0.0) go to 10
      if(ddd>0.0) then
          yxdist=dist* sqrt(ddd)
      else
          yxdist=0.d0
      end if
      !if(yxdist.gt.1.0d-9) go to 20
      if(yxdist<=1.0d-9) then
        !10 continue
          xi2=xi1
          xl2=xl1
          yi2=yi1
          yl2=yl1
          costh=cosa
          sinth=0.d0
          !go to 30
      else
        !20 cosph=yj1/yxdist
          cosph=yj1/yxdist
          sinph=xj1/yxdist
          xi2=xi1*cosph-yi1*sinph
          xj2=xj1*cosph-yj1*sinph
          xl2=xl1*cosph-yl1*sinph
          yi2=xi1*sinph+yi1*cosph
          yj2=xj1*sinph+yj1*cosph
          yl2=xl1*sinph+yl1*cosph
          !      rotate kj around the x axis so kj lies along the z axis
          costh=cosa
          sinth=yj2/dist
      end if
    !30 continue
      yi3=yi2*costh-zi1*sinth
      yl3=yl2*costh-zl1*sinth
      call dang(xl2,yl3,xi2,yi3,angle)
      if (angle .lt. 0.) angle=6.2831853d0+angle
      if (angle .ge. 6.2831853d0 ) angle=0.d0
      return
    end subroutine dihed
    
    !> bangle calculates the angle between atoms i,j, and k. the
    !! cartesian coordinates are in xyz.
    subroutine bangle(numat,xyz,i,j,k,angle)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: numat, i, j, k
      real(kind=8),dimension(3,numat),intent(in) :: xyz
      real(kind=8),intent(out) :: angle
    
      ! Local variables
      real(kind=8) :: d2ij, d2jk, d2ik, xy, temp
      
      d2ij = (xyz(1,i)-xyz(1,j))**2+&
           (xyz(2,i)-xyz(2,j))**2+&
           (xyz(3,i)-xyz(3,j))**2
      d2jk = (xyz(1,j)-xyz(1,k))**2+&
           (xyz(2,j)-xyz(2,k))**2+&
           (xyz(3,j)-xyz(3,k))**2
      d2ik = (xyz(1,i)-xyz(1,k))**2+&
           (xyz(2,i)-xyz(2,k))**2+&
           (xyz(3,i)-xyz(3,k))**2
      xy = sqrt(d2ij*d2jk)
      temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
      if (temp .gt. 1.0d0) temp=1.0d0
      if (temp .lt. -1.0d0) temp=-1.0d0
      angle = acos( temp )
    end subroutine bangle
    
    !>    dang  determines the angle between the points (a1,a2), (0,0),
    !!          and (b1,b2).  the result is put in rcos.
    subroutine dang(a1,a2,b1,b2,rcos)
      implicit none
    
      ! Calling arguments
      real(kind=8),intent(inout) :: a1, a2, b1, b2
      real(kind=8),intent(out) :: rcos
    
      ! Local variables
      real(kind=8) :: pi, zero, anorm, bnorm, sinth, costh
    
      pi=2.0d0* asin(1.0d00)
      zero=1.0d-6
      if(( abs(a1).lt.zero.and. abs(a2).lt.zero) .or. &
        ( abs(b1).lt.zero.and. abs(b2).lt.zero) ) then
          rcos=0.0d0
      else
          anorm=1.0d0/ sqrt(a1**2+a2**2)
          bnorm=1.0d0/ sqrt(b1**2+b2**2)
          a1=a1*anorm
          a2=a2*anorm
          b1=b1*bnorm
          b2=b2*bnorm
          sinth=(a1*b2)-(a2*b1)
          costh=a1*b1+a2*b2
          if(costh.gt.1.0d0) costh=1.0d0
          if(costh.lt.-1.0d0) costh=-1.0d0
          rcos= acos(costh)
          if( abs(rcos)>=4.0d-4) then
              if(sinth.gt.0.d0) rcos=6.2831853d0-rcos
              rcos=-rcos
              return
          else
        !10 rcos=0.0d0
              rcos=0.0d0
          end if
      end if
    end subroutine dang
    
    !>   xyzgeo converts coordinates from cartesian to internal.
    !!
    !!     on input xyz  = array of cartesian coordinates
    !!              numat= number of atoms
    !!              na   = numbers of atom to which atoms are related
    !!                     by distance
    !!              nb   = numbers of atom to which atoms are related
    !!                     by angle
    !!              nc   = numbers of atom to which atoms are related
    !!                     by dihedral
    !!
    !!    on output geo  = internal coordinates in angstroms, radians,
    !!                     and radians
    subroutine xyzgeo(xyz,numat,na,nb,nc,degree,geo)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: numat
      real(kind=8),intent(in) :: degree
      real(kind=8),dimension(3,numat),intent(in) :: xyz
      real(kind=8),dimension(3,numat),intent(out) :: geo
      integer,dimension(numat),intent(in) :: na, nb, nc
    
      ! Local variables
      integer :: i, j, k, l, ii
    
      do i=2,numat
         j=na(i)
         k=nb(i)
         l=nc(i)
         geo(1,i)= sqrt((xyz(1,i)-xyz(1,j))**2+&
              (xyz(2,i)-xyz(2,j))**2+&
              (xyz(3,i)-xyz(3,j))**2)
         if(i.lt.3) cycle
         ii=i
         call bangle(numat,xyz,ii,j,k,geo(2,i))
         geo(2,i)=geo(2,i)*degree
         if(i.lt.4) cycle
         call dihed(numat,xyz,ii,j,k,l,geo(3,i))
         geo(3,i)=geo(3,i)*degree
      end do
    
      geo(1,1)=0.d0
      geo(2,1)=0.d0
      geo(3,1)=0.d0
      geo(2,2)=0.d0
      geo(3,2)=0.d0
      geo(3,3)=0.d0
    
    end subroutine xyzgeo
    
    !> xyzint works out the internal coordinates of a molecule.
    !!        the "rules" for the connectivity are as follows:
    !!        atom i is defined as being at a distance from the nearest
    !!        atom j, atom j already having been defined.
    !!        atom i makes an angle with atom j and the atom k, which has
    !!        already been defined, and is the nearest atom to j
    !!        atom i makes a dihedral angle with atoms j, k, and l. l having
    !!        been defined and is the nearest atom to k
    !!
    !!        note that geo and xyz must not be the same in the call.
    !!
    !!   on input xyz    = cartesian array of numat atoms
    !!            degree = 1 if angles are to be in radians
    !!            degree = 57.29578 if angles are to be in radians
    subroutine xyzint(xyz,numat,na,nb,nc,degree,geo)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: numat
      real(kind=8),intent(in) :: degree
      real(kind=8),dimension(3,numat),intent(in) :: xyz
      real(kind=8),dimension(3,numat),intent(out) :: geo
      integer,dimension(numat),intent(out) :: na, nb, nc
    
      ! Local variables
      integer :: nai1, nai2, i, j, im1, k
      real(kind=8) :: sum, r
    
      nai1=0
      nai2=0
      do  i=1,numat
         na(i)=2
         nb(i)=3
         nc(i)=4
         im1=i-1
         if(im1.eq.0) cycle
         sum=100.d0
         do  j=1,im1
            r=(xyz(1,i)-xyz(1,j))**2+&
                 (xyz(2,i)-xyz(2,j))**2+&
                 (xyz(3,i)-xyz(3,j))**2
            if(r.lt.sum.and.na(j).ne.j.and.nb(j).ne.j) then
               sum=r
               k=j
            endif
         end do
         !
         !   atom i is nearest to atom k
         !
         na(i)=k
         if(i.gt.2)nb(i)=na(k)
         if(i.gt.3)nc(i)=nb(k)
         !
         !   find any atom to relate to na(i)
         !
      end do
      na(1)=0
      nb(1)=0
      nc(1)=0
      nb(2)=0
      nc(2)=0
      nc(3)=0
      call xyzgeo(xyz,numat,na,nb,nc,degree,geo)
    
    end subroutine xyzint
    
    
    
    
    subroutine internal_to_cartesian(nat, xyz_int, xyz_cart)
      use vector_operations
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: nat
      real(kind=8),dimension(3,nat),intent(in) :: xyz_int
      real(kind=8),dimension(3,nat),intent(out) :: xyz_cart
    
      ! Local variables
      integer :: iat
      real(kind=8) :: tt
      real(kind=8),dimension(3) :: vector, vector1, vector2
    
      do iat=1,nat
    
          if (iat==1) then
              ! first atom, place it at the origin
              xyz_cart(1,1)=0.d0
              xyz_cart(2,1)=0.d0
              xyz_cart(3,1)=0.d0
          else if (iat==2) then
              ! second atom, put it along the z axis
              xyz_cart(1,2)=0.d0
              xyz_cart(2,2)=0.d0
              xyz_cart(3,2)=xyz_int(1,2)
          else if (iat==3) then
              ! third atom, put at the given distance and angle in the yz-plane
              xyz_cart(1,3) = 0.d0
              !!write(*,*) 'cos(xyz_int(2,3))',cos(xyz_int(2,3))
              !!write(*,*) 'sin(xyz_int(2,3))',sin(xyz_int(2,3))
              !!xyz_cart(2,3) = xyz_cart(3,2) - xyz_int(1,3)*cos(xyz_int(2,3))
              !!xyz_cart(3,3) = xyz_cart(3,2) + xyz_int(1,3)*sin(xyz_int(2,3))
              xyz_cart(2,3) = 0.d0
              xyz_cart(3,3) = xyz_cart(3,2) + xyz_int(1,3)
              !xyz_cart(2,3) = -cos(xyz_int(2,3))*xyz_cart(2,3) + sin(xyz_int(2,3))*xyz_cart(3,3)
              !xyz_cart(3,3) = -cos(xyz_int(2,3))*xyz_cart(3,3) - sin(xyz_int(2,3))*xyz_cart(2,3)
              ! First shift the rotation center (point C) to the origin
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat) - xyz_cart(1:3,iat-1)
              !!write(*,'(a,3f9.2)') 'before rot: xyz_cart(1:3,iat)',xyz_cart(1:3,iat)
              vector1(1)=0.d0 ; vector1(2)=xyz_cart(2,iat) ; vector1(3)=xyz_cart(3,iat)
              vector2(1)=0.d0 ; vector2(2)=-xyz_cart(3,iat) ; vector2(3)=xyz_cart(2,iat)
              xyz_cart(1:3,3) = cos(xyz_int(2,3))*vector1(1:3) + sin(xyz_int(2,3))*vector2(1:3)
              !xyz_cart(2,3) = cos(xyz_int(2,3))*xyz_cart(2,3) - sin(xyz_int(2,3))*xyz_cart(3,3)
              !xyz_cart(3,3) = cos(xyz_int(2,3))*xyz_cart(3,3) + sin(xyz_int(2,3))*xyz_cart(2,3)
              !!write(*,'(a,3f9.2)') 'after rot: xyz_cart(1:3,iat)',xyz_cart(1:3,iat)
              ! Undo the shift
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat) + xyz_cart(1:3,iat-1)
          else
              ! General case
    
              ! First put the new atom D at the distance "bond" away from atom C, extending
              ! along the axis BC.
              vector(1:3) = xyz_cart(1:3,iat-1)-xyz_cart(1:3,iat-2)
              tt = sqrt( vector(1)**2 + vector(2)**2 + vector(3)**2 )
              vector(1:3) = vector(1:3) / tt
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat-1) + xyz_int(1,iat)*vector(1:3)
              !!write(*,'(a,3f9.2)') 'xyz_cart(1:3,iat-1)', xyz_cart(1:3,iat-1)
              !!write(*,'(a,3f9.2)') 'vector', vector(1:3)
              !!write(*,'(a,3f9.2)') 'xyz_cart(1:3,iat)', xyz_cart(1:3,iat)
    
              ! Then rotate around C in the ABC plane
    
              ! Determine the rotation axis
              ! n = AB x bcn / |AB x bcn|, with bcn = BC/|BC|
              vector1(1:3) = xyz_cart(1:3,iat-2) - xyz_cart(1:3,iat-3)
              vector2(1:3) = xyz_cart(1:3,iat-1) - xyz_cart(1:3,iat-2)
              tt = sqrt( vector2(1)**2 + vector2(2)**2 + vector2(3)**2 )
              vector2(1:3) = vector2(1:3) / tt
              vector(1:3) = cross_product(vector1, vector2)
              tt = sqrt( vector(1)**2 + vector(2)**2 + vector(3)**2 )
              vector(1:3) = vector(1:3) / tt
    
              ! Apply the rotation
              ! First shift the rotation center (point C) to the origin
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat) - xyz_cart(1:3,iat-1)
              vector1(1:3) = cross_product(vector, xyz_cart(1:3,iat))
              vector2(1:3) = cross_product(vector1, vector)
              tt = vector(1)*xyz_cart(1,iat) + vector(2)*xyz_cart(2,iat) + vector(3)*xyz_cart(3,iat)
              xyz_cart(1:3,iat) = tt*vector(1:3) + cos(xyz_int(2,iat))*vector2(1:3) + sin(xyz_int(2,iat))*vector1(1:3)
              ! Undo the shift
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat) + xyz_cart(1:3,iat-1)
    
    
              ! Then rotate around BC
              ! Determine the rotation axis
              vector(1:3) = xyz_cart(1:3,iat-1) - xyz_cart(1:3,iat-2)
              tt = sqrt( vector(1)**2 + vector(2)**2 + vector(3)**2 )
              vector(1:3) = vector(1:3) / tt
    
              ! Apply the rotation
              ! First shift the rotation center (point C) to the origin
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat) - xyz_cart(1:3,iat-1)
              vector1(1:3) = cross_product(vector, xyz_cart(1:3,iat))
              vector2(1:3) = cross_product(vector1, vector)
              tt = vector(1)*xyz_cart(1,iat) + vector(2)*xyz_cart(2,iat) + vector(3)*xyz_cart(3,iat)
              xyz_cart(1:3,iat) = tt*vector(1:3) + cos(xyz_int(3,iat))*vector2(1:3) + sin(xyz_int(3,iat))*vector1(1:3)
              ! Undo the shift
              xyz_cart(1:3,iat) = xyz_cart(1:3,iat) + xyz_cart(1:3,iat-1)
    
          end if
    
    
      end do
    
    end subroutine internal_to_cartesian

  end module internal_coordinates
