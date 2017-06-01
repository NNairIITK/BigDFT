!!subroutine getEffectiveFilterSextic(it,parabPrefac,hgrid, x0, eff, filterCode)
!!!
!!! Purpose:
!!! ========
!!!   Calculates the effective filter for the operator [kineticEnergy + (x-x0)^4].
!!!   
!!! Calling arguments:
!!! ==================
!!!   Input arguments:
!!!     hgrid  grid spacing
!!!     x0     the center of the parabolic potential (x-x0)^2
!!!   Output arguments:
!!!     aeff   the effective filter for <phi|Op|phi>
!!!     beff   the effective filter for <psi|Op|phi>
!!!     ceff   the effective filter for <phi|Op|psi>
!!!     eeff   the effective filter for <psi|Op|psi>
!!!
!!!use filterModule2
!!use filterModule
!!implicit none
!!
!!! Calling arguments
!!integer, intent(in):: it
!!real(8),intent(in):: parabPrefac, hgrid, x0
!!real(8),dimension(lb:ub),intent(out):: eff
!!character(len=*):: filterCode
!!
!!! Local variables
!!integer:: i
!!real(8):: fac, fac2, prefac1, prefac2a, hgrid2, hgrid3, hgrid4, hgrid5, x02, x03, x04, x05
!!real(8):: scale
!!
!!scale=1.d0
!!!scale=5.d-2
!!
!!prefac1=-.5d0/hgrid**2
!!!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!!!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
!!fac=parabPrefac*scale
!!fac2=parabPrefac*hgrid*scale
!!hgrid2=hgrid**2
!!hgrid3=hgrid**3
!!hgrid4=hgrid**4
!!hgrid5=hgrid**5
!!x02=x0**2
!!x03=x0**3
!!x04=x0**4
!!x05=x0**5
!!
!!! Determine which filter we have to calculate
!!select case(trim(filterCode))
!!case('a')
!!    do i=lb,ub
!!        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
!!        eff(i)=prefac1*a(i) + fac2*( hgrid5*a6(i) + 6.d0*hgrid4*x0*a5(i) + 15.d0*hgrid3*x02*a4(i) &
!!               + 20.d0*hgrid2*x03*a3(i) + 15.d0*hgrid*x04*a2(i) + 6.d0*x05*a1(i))
!!    end do
!!    !eff(0)=eff(0)+fac*x0**2
!!    eff(0)=eff(0)+fac*x0**6
!!case('b')
!!    do i=lb,ub
!!        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
!!        eff(i)=prefac1*b(i) + fac2*( hgrid5*b6(i) + 6.d0*hgrid4*x0*b5(i) + 15.d0*hgrid3*x02*b4(i) &
!!               + 20.d0*hgrid2*x03*b3(i) + 15.d0*hgrid*x04*b2(i) + 6.d0*x05*b1(i))
!!        !eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
!!    end do
!!case('c')
!!    do i=lb,ub
!!        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
!!        eff(i)=prefac1*c(i) + fac2*( hgrid5*c6(i) + 6.d0*hgrid4*x0*c5(i) + 15.d0*hgrid3*x02*c4(i) &
!!               + 20.d0*hgrid2*x03*c3(i) + 15.d0*hgrid*x04*c2(i) + 6.d0*x05*c1(i))
!!        !eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
!!    end do
!!case('e')
!!    do i=lb,ub
!!        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
!!        eff(i)=prefac1*e(i) + fac2*( hgrid5*e6(i) + 6.d0*hgrid4*x0*e5(i) + 15.d0*hgrid3*x02*e4(i) &
!!               + 20.d0*hgrid2*x03*e3(i) + 15.d0*hgrid*x04*e2(i) + 6.d0*x05*e1(i))
!!        !eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
!!    end do
!!    !eff(0)=eff(0)+fac*x0**2
!!    eff(0)=eff(0)+fac*x0**6
!!case default
!!    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
!!    stop
!!end select
!!
!!
!!end subroutine getEffectiveFilterSextic



!!subroutine getFilterSextic(it,parabPrefac,hgrid, x0, eff, filterCode)
!!!
!!! Purpose:
!!! ========
!!!   Calculates the effective filter for the operator [kineticEnergy + (x-x0)^4].
!!!   
!!! Calling arguments:
!!! ==================
!!!   Input arguments:
!!!     hgrid  grid spacing
!!!     x0     the center of the parabolic potential (x-x0)^2
!!!   Output arguments:
!!!     aeff   the effective filter for <phi|Op|phi>
!!!     beff   the effective filter for <psi|Op|phi>
!!!     ceff   the effective filter for <phi|Op|psi>
!!!     eeff   the effective filter for <psi|Op|psi>
!!!
!!!use filterModule2
!!use filterModule
!!implicit none
!!
!!! Calling arguments
!!integer, intent(in):: it
!!real(8),intent(in):: parabPrefac, hgrid, x0
!!real(8),dimension(lb:ub),intent(out):: eff
!!character(len=*):: filterCode
!!
!!! Local variables
!!integer:: i
!!real(8):: fac, fac2, prefac1, prefac2a, hgrid2, hgrid3, hgrid4, hgrid5, x02, x03, x04, x05
!!real(8):: scale
!!
!!scale=1.d0
!!!scale=5.d-2
!!
!!prefac1=-.5d0/hgrid**2
!!!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!!!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
!!fac=parabPrefac*scale
!!fac2=parabPrefac*hgrid*scale
!!hgrid2=hgrid**2
!!hgrid3=hgrid**3
!!hgrid4=hgrid**4
!!hgrid5=hgrid**5
!!x02=x0**2
!!x03=x0**3
!!x04=x0**4
!!x05=x0**5
!!
!!! Determine which filter we have to calculate
!!select case(trim(filterCode))
!!case('a')
!!    do i=lb,ub
!!        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
!!        eff(i)=fac2*( hgrid5*a6(i) + 6.d0*hgrid4*x0*a5(i) + 15.d0*hgrid3*x02*a4(i) &
!!               + 20.d0*hgrid2*x03*a3(i) + 15.d0*hgrid*x04*a2(i) + 6.d0*x05*a1(i))
!!    end do
!!    !eff(0)=eff(0)+fac*x0**2
!!    eff(0)=eff(0)+fac*x0**6
!!case('b')
!!    do i=lb,ub
!!        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
!!        eff(i)=fac2*( hgrid5*b6(i) + 6.d0*hgrid4*x0*b5(i) + 15.d0*hgrid3*x02*b4(i) &
!!               + 20.d0*hgrid2*x03*b3(i) + 15.d0*hgrid*x04*b2(i) + 6.d0*x05*b1(i))
!!        !eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
!!    end do
!!case('c')
!!    do i=lb,ub
!!        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
!!        eff(i)=fac2*( hgrid5*c6(i) + 6.d0*hgrid4*x0*c5(i) + 15.d0*hgrid3*x02*c4(i) &
!!               + 20.d0*hgrid2*x03*c3(i) + 15.d0*hgrid*x04*c2(i) + 6.d0*x05*c1(i))
!!        !eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
!!    end do
!!case('e')
!!    do i=lb,ub
!!        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
!!        eff(i)=fac2*( hgrid5*e6(i) + 6.d0*hgrid4*x0*e5(i) + 15.d0*hgrid3*x02*e4(i) &
!!               + 20.d0*hgrid2*x03*e3(i) + 15.d0*hgrid*x04*e2(i) + 6.d0*x05*e1(i))
!!        !eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
!!    end do
!!    !eff(0)=eff(0)+fac*x0**2
!!    eff(0)=eff(0)+fac*x0**6
!!case default
!!    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
!!    stop
!!end select
!!
!!
!!end subroutine getFilterSextic
