subroutine md_isokinetic_init(amass, mditemp, natom, vel)

  use defs_basis
  use abi_interfaces_lowlevel

  implicit none

  integer, intent(in) :: natom
  real(dp),intent(in) :: mditemp
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(inout) :: vel(3,natom)

  interface
     function uniformrandom(seed) 
       implicit none
       integer :: seed
       double precision :: uniformrandom
     end function uniformrandom
  end interface

  real(dp),parameter :: v2tol=tol8
  integer :: idim, iatom, idum=-5
  real(dp) :: v2gauss, s1, s2, vtest, sigma2, rescale_vel
  character(len=500) :: message

  !  v2gauss is twice the kinetic energy
  v2gauss=0.0_dp
  do iatom=1,natom
     do idim=1,3
        v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
  end do

  !  If there is no kinetic energy
  if (v2gauss<=v2tol) then
     !   Maxwell-Boltzman distribution
     v2gauss=zero
     vtest=zero
     do iatom=1,natom
        do idim=1,3
           vel(idim,iatom)=sqrt(kb_HaK*mditemp/amass(iatom))*cos(two_pi*uniformrandom(idum))
           vel(idim,iatom)=vel(idim,iatom)*sqrt(-2._dp*log(uniformrandom(idum)))
        end do
     end do

     !   Get rid of center-of-mass velocity
     s1=sum(amass(:))
     do idim=1,3
        s2=sum(amass(:)*vel(idim,:))
        vel(idim,:)=vel(idim,:)-s2/s1
     end do

     !   Recompute v2gauss
     do iatom=1,natom
        do idim=1,3
           v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
           vtest=vtest+vel(idim,iatom)/(3._dp*natom)
        end do
     end do

     !   Now rescale the velocities to give the exact temperature
     rescale_vel=sqrt(3._dp*natom*kb_HaK*mditemp/v2gauss)
     vel(:,:)=vel(:,:)*rescale_vel

     !   Recompute v2gauss with the rescaled velocities
     v2gauss=zero
     do iatom=1,natom
        do idim=1,3
           v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
        end do
     end do

     !   Compute the variance and print
     sigma2=(v2gauss/(3._dp*natom)-amass(1)*vtest**2)/kb_HaK

     write(message, '(a)' )&
          &    ' Rescaling or initializing velocities to initial temperature'
     call abi_wrtout(ab_out,message,'COLL')
     call abi_wrtout(std_out,message,'COLL')
     write(message, '(a,d12.5,a,D12.5)' )&
          &    ' --- Scaling factor :',rescale_vel,' Asked T (K) ',mditemp
     call abi_wrtout(ab_out,message,'COLL')
     call abi_wrtout(std_out,message,'COLL')
     write(message, '(a,d12.5,a,D12.5)' )&
          &    ' --- Effective temperature',v2gauss/(3*natom*kb_HaK),' From variance', sigma2
     call abi_wrtout(ab_out,message,'COLL')
     call abi_wrtout(std_out,message,'COLL')
  end if
end subroutine md_isokinetic_init

subroutine md_isokinetic(acell, amass, dtion, epot, fcart, itime, natom, &
     & mditemp, me, rprimd, vel, vel_nexthalf, xcart, xcart_next, xred_next)

  use defs_basis

  implicit none

  integer, intent(in) :: natom, itime, me
  real(dp), intent(in) :: mditemp, dtion
  real(dp), intent(out) :: epot
  real(dp), intent(in) :: rprimd(3,3), acell(3)
  real(dp), intent(in) :: amass(natom)
  real(dp), intent(inout) :: vel(3,natom), fcart(3, natom)
  real(dp), intent(out) :: vel_nexthalf(3, natom)
  real(dp), intent(in) :: xcart(3, natom)
  real(dp), intent(out) :: xred_next(3,natom), xcart_next(3,natom)

  real(dp), parameter :: v2tol=tol8
  integer :: idim, iatom
  real(dp) :: v2gauss, s1, s2, a, b, sqb, as, s, scdot
  real(dp), allocatable :: fcart_m(:,:)

  !  Application of Gauss' principle of least constraint according to Fei Zhang's algorithm (J. Chem. Phys. 106, 1997, p.6102)
  if (itime == 0) then
     call md_isokinetic_init(amass, mditemp, natom, vel)
     vel_nexthalf(:,:)=vel(:,:)
     xcart_next(:,:)=xcart(:,:)
     call xredxcart(natom,-1,rprimd,xcart_next,xred_next)

     return
  end if

  allocate(fcart_m(3,natom))
  v2gauss = zero
  do iatom=1,natom
     do idim=1,3
        fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
        v2gauss=v2gauss+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
  end do

  !   Computation of vel_nexthalf (4.16 de Ref.1)
  !   Computation of a and b (4.13 de Ref.1)
  a=0.0_dp
  b=0.0_dp
  do iatom=1,natom
     do idim=1,3
        a=a+fcart_m(idim,iatom)*vel(idim,iatom)*amass(iatom)
        b=b+fcart_m(idim,iatom)*fcart_m(idim,iatom)*amass(iatom)
     end do
  end do
  a=a/v2gauss
  b=b/v2gauss
  !   Computation of s and scdot
  sqb=sqrt(b)
  as=sqb*dtion/2.
  s1=cosh(as)
  s2=sinh(as)
  s=a*(s1-1.)/b+s2/sqb
  scdot=a*s2/sqb+s1
  vel_nexthalf(:,:)=(vel(:,:)+fcart_m(:,:)*s)/scdot

  !   Computation of the next positions
  xcart_next(:,:)=xcart(:,:)+vel_nexthalf(:,:)*dtion

  !   Convert back to xred (reduced coordinates)
  call xredxcart(natom,-1,rprimd,xcart_next,xred_next)

  !   Computation of the forces for the new positions
  !   Compute LDA forces (big loop), fcart_m is used as dummy argument for fred
  call scfloop_main(acell, epot, fcart, fcart_m, itime, me, natom, rprimd, xred_next)
  do iatom=1,natom
     do idim=1,3
        fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
  end do

  !   Computation of vel(:,:) at the next positions
  !   Computation of v2gauss
  v2gauss=0.0_dp
  do iatom=1,natom
     do idim=1,3
        v2gauss=v2gauss+vel_nexthalf(idim,iatom)*vel_nexthalf(idim,iatom)*amass(iatom)
     end do
  end do
  !   Calcul de a et b (4.13 de Ref.1)
  a=0.0_dp
  b=0.0_dp
  do iatom=1,natom
     do idim=1,3
        a=a+fcart_m(idim,iatom)*vel_nexthalf(idim,iatom)*amass(iatom)
        b=b+fcart_m(idim,iatom)*fcart_m(idim,iatom)*amass(iatom)
     end do
  end do
  a=a/v2gauss
  b=b/v2gauss
  !   Calcul de s et scdot
  sqb=sqrt(b)
  as=sqb*dtion/2.
  s1=cosh(as)
  s2=sinh(as)
  s=a*(s1-1.)/b+s2/sqb
  scdot=a*s2/sqb+s1
  vel(:,:)=(vel_nexthalf(:,:)+fcart_m(:,:)*s)/scdot
  !FB20090429    
  !  Convert input xred (reduced coordinates) to xcart (cartesian)
  !   call xredxcart(natom,1,rprimd,xcart,xred)
  !FB20090429    
  deallocate(fcart_m)

  !  End of case ionmov = 12

 end subroutine md_isokinetic
