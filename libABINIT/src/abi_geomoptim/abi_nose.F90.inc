subroutine abi_md_nose_init(amass, natom, temp, vel)
  
  use abi_defs_basis
  use abi_interfaces_lowlevel

  implicit none

  integer, intent(in) :: natom
  real(dp),intent(in) :: temp
  real(dp),intent(in) :: amass(natom)
  real(dp),intent(inout) :: vel(3,natom)

  interface
     function abi_uniformrandom(seed) 
       implicit none
       integer :: seed
       double precision :: abi_uniformrandom
     end function abi_uniformrandom
  end interface


  integer :: iatom, idim, idum=-5
  real(dp),parameter :: v2tol=tol8
  real(dp) :: v2nose, rescale_vel
  character(len=500) :: message

  !   Compute twice the kinetic energy of the system, called v2nose
  v2nose=0.0_dp
  do iatom=1,natom
     do idim=1,3
        v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
  end do

  !   If there is no kinetic energy, use a random initial velocity
  if (v2nose<=v2tol) then
     v2nose=0.0_dp
     do iatom=1,natom
        do idim=1,3
           !      abi_uniformrandom returns a uniform random deviate between 0.0 and 1.0
           !      if it were always 0 or 1, then the following expression
           !      would give the requested temperature
           vel(idim,iatom)=(1.0_dp-2.0_dp*abi_uniformrandom(idum))*&
                &       sqrt( temp*kb_HaK / amass(iatom) )
           !      Recompute v2nose
           v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
        end do
     end do
  end if

  !   Now, rescale the velocities to give the proper temperature
  rescale_vel=sqrt(3.0_dp*natom*temp*kb_HaK/v2nose)
  vel(:,:)=vel(:,:)*rescale_vel
  !   Recompute v2nose with the rescaled velocities
  v2nose=0.0_dp
  do iatom=1,natom
     do idim=1,3
        v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
  end do
  write(message, '(a)' )&
       &    ' Rescaling or initializing velocities to initial temperature'
  call abi_wrtout(ab_out,message,'COLL')
  call abi_wrtout(std_out,message,'COLL')
  write(message, '(2(a,es22.14))' )&
       &    ' ---  Scaling factor : ',rescale_vel,' Asked T (K) ', temp
  call abi_wrtout(ab_out,message,'COLL')
  call abi_wrtout(std_out,message,'COLL')
  write(message, '(a,es22.14)' )&
       &    ' ---  Effective temperature',v2nose/(3.0_dp*natom*kb_HaK)
  call abi_wrtout(ab_out,message,'COLL')
  call abi_wrtout(std_out,message,'COLL')

end subroutine abi_md_nose_init

subroutine abi_md_nose(amass, dtion, fcart, fcart_mold, gnose, itime, ktemp, mditemp, &
     & natom, noseinert, rprimd, snose, v2nose, vel, xcart, xcart_next, &
     & xi_nose, xred_next)
  
  use abi_defs_basis
  use abi_interfaces_lowlevel
  use abi_interfaces_geometry

  implicit none

  integer, intent(in) :: natom, itime
  real(dp),intent(in) :: mditemp, dtion, ktemp, noseinert
  real(dp), intent(inout) :: gnose, snose, xi_nose, v2nose
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: amass(natom), fcart(3, natom)
  real(dp),intent(inout) :: vel(3,natom)
  real(dp),intent(in) :: xcart(3, natom)
  real(dp),intent(out) :: xred_next(3,natom), xcart_next(3,natom)
  real(dp), intent(inout) :: fcart_mold(3, natom)

  real(dp),parameter :: nosetol=tol10
  logical :: ready
  integer :: idim, iatom
  real(dp) :: dnose, xin_nose, fsnose, xio, delxi
  real(dp), allocatable :: fcart_m(:,:), finose(:,:), vel_temp(:,:), vonose(:,:)
  real(dp), allocatable :: hnose(:,:), binose(:,:)
  character(len=500) :: message

  if(itime==0)then
     snose=0.0_dp
     xi_nose=0.0_dp

     call abi_md_nose_init(amass, natom, mditemp, vel)
  end if

  allocate(fcart_m(3,natom))
  do iatom=1,natom
     do idim=1,3
        fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
  end do

  !  First step of velocity verlet algorithm
  gnose=3*natom

  !  Calculate nose-hoover force on atoms
  !  If first iteration, no old force are available, so use present forces
!  if (itime==0) fcart_mold(:,:)=fcart_m(:,:)
  fcart_mold(:,:)=fcart_m(:,:)

  allocate(finose(3, natom))
  finose(:,:)=fcart_mold(:,:)-xi_nose*vel(:,:)
  xcart_next(:,:)=xcart(:,:)+dtion*(vel(:,:)+dtion*finose(:,:)/2.0_dp)

  !  Calculate v2nose
  v2nose=0.0_dp
  do iatom=1,natom
     do idim=1,3
        v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
  end do
  vel(:,:)=vel(:,:)+dtion*finose(:,:)/2.0_dp
  deallocate(finose)

  !  Update thermostat
  fsnose=(v2nose-gnose*ktemp)/noseinert
  snose=snose+dtion*(xi_nose+dtion*fsnose/2.0_dp)
  xi_nose=xi_nose+dtion*fsnose/2.0_dp

  !  Second step of the velocity Verlet algorithm, uses the 'new forces'
  !  Calculate v2nose
  v2nose=0.0_dp
  do iatom=1,natom
     do idim=1,3
        v2nose=v2nose+vel(idim,iatom)*vel(idim,iatom)*amass(iatom)
     end do
  end do
  allocate(vel_temp(3, natom))
  vel_temp(:,:)=vel(:,:)

  xin_nose=xi_nose

  !  Start Newton-Raphson loop
  allocate(vonose(3, natom))
  allocate(hnose(3, natom))
  allocate(binose(3, natom))
  ready=.false.
  do while (.not.ready)
     xio=xin_nose
     delxi=0.0D0
     vonose(:,:)=vel_temp(:,:)
     hnose(:,:)=-dtion/2.0_dp*(fcart_m(:,:)-xio*vonose(:,:))-(vel(:,:)-vonose(:,:))
     do iatom=1,natom
        do idim=1,3
           binose(idim,iatom)=vonose(idim,iatom)*dtion/noseinert*amass(iatom) ! a verifier
           delxi=delxi+hnose(idim,iatom)*binose(idim,iatom)
        end do
     end do
     dnose=-(xio*dtion/2.0D0+1.0D0)
     delxi=delxi-dnose*((-v2nose+gnose*ktemp)*dtion/2.0_dp/ &
          &    noseinert-(xi_nose-xio))
     delxi=delxi/(-dtion*dtion/2.0_dp*v2nose/noseinert+dnose)
     v2nose=0.0_dp

     vel_temp(:,:)=vel_temp(:,:)+(hnose+dtion/2.0_dp*vonose(:,:)*delxi)/dnose
     do iatom=1,natom
        do idim=1,3
           v2nose=v2nose+vel_temp(idim,iatom)*vel_temp(idim,iatom)*amass(iatom)
        end do
     end do
     !   New guess for xi
     xin_nose=xio+delxi

     ready=.true.
     !   Test for convergence
     iatom=0
     idim=1
     do while((iatom<=natom).and.(idim<=3).and.ready)
        iatom=iatom+1
        if (iatom>natom) then
           iatom=1
           idim=idim+1
        end if
        if ((iatom<=natom) .and.(idim<=3)) then
           if (abs(vel_temp(idim,iatom))<1.0d-50) vel_temp(idim,iatom)=1.0d-50
           if (abs((vel_temp(idim,iatom)-vonose(idim,iatom))/vel_temp(idim,iatom))&
                &      >nosetol) ready=.false.
        else
           if (xin_nose<1.0d-50) xin_nose=1.0d-50
           if (abs((xin_nose-xio)/xin_nose)>nosetol) ready=.false.
        end if
     end do   ! end of while
     !   Enddo ready
  end do
  deallocate(vonose)
  deallocate(hnose)
  deallocate(binose)

  !  Update velocities to converged value
  vel(:,:)=vel_temp(:,:)
  deallocate(vel_temp)

  write(message, '(a,f12.3,a,f12.3,a)' ) &
       & ' Nose-Hoover updated velocities (T_ask = ', &
       & ktemp/kb_HaK, ' and T_eff = ', v2nose/(3.0_dp*natom*kb_HaK), ')'
  call abi_wrtout(ab_out,message,'COLL')
  call abi_wrtout(std_out,message,'COLL')

  !  Update thermostat
  xi_nose=xin_nose
  !  Convert back to xred_next (reduced coordinates)
  call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)
  !  Store 'new force' as 'old force'
  fcart_mold(:,:)=fcart_m(:,:)
  deallocate(fcart_m)

  !  End of case ionmov =8

end subroutine abi_md_nose

subroutine abi_md_nose_finalise(etotal, gnose, itime, ktemp, noseinert, snose, &
     & v2nose, xi_nose)
  
  use abi_defs_basis
  use abi_interfaces_lowlevel

  implicit none

  integer, intent(in) :: itime
  real(dp), intent(in) :: etotal, gnose, snose, v2nose, xi_nose, ktemp, noseinert

  real(dp) :: hamnose
  character(len=500) :: message

  !  compute conserved quantity
  hamnose=etotal+v2nose/2.0_dp+&
       &   (xi_nose**2*noseinert)/2.0_dp+gnose*ktemp*snose
  write(message, '(a,i6,a,es22.14,a)' )&
       &   ' At the end of Moldyn step',itime,', ham=',hamnose,' Ha.'
  call abi_wrtout(ab_out,message,'COLL')
  call abi_wrtout(std_out,message,'COLL')
end subroutine abi_md_nose_finalise
