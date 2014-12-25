subroutine md_quenched_stop_atoms(amass, dtion, ekin, fcart, iatfix, itime, &
     & natom, nstopped, rprimd, vel, vel_prevhalf, vel_nexthalf, &
     & xcart, xcart_next, xred_next)

  use defs_basis
  use interfaces_14_hidewrite

  implicit none

  integer, intent(in) :: natom, itime
  integer, intent(out) :: nstopped
  integer, intent(in) :: iatfix(3, natom)
  real(dp), intent(in) :: dtion
  real(dp), intent(out) :: ekin
  real(dp), intent(in) :: amass(natom), rprimd(3,3)
  real(dp), intent(in) :: fcart(3, natom)
  real(dp), intent(inout) :: vel(3,natom)
  real(dp), intent(out) :: vel_prevhalf(3,natom), vel_nexthalf(3, natom)
  real(dp), intent(in) :: xcart(3, natom)
  real(dp), intent(out) :: xcart_next(3,natom), xred_next(3, natom)

  integer :: iatom, istopped, ii, idir
  real(dp) :: scprod
  integer, allocatable :: stopped(:)
  character(len=500) :: message
  
  allocate(stopped(natom))
  stopped(:)=0
  do iatom=1,natom
     scprod=fcart(1,iatom)*vel(1,iatom)+&
          &    fcart(2,iatom)*vel(2,iatom)+&
          &    fcart(3,iatom)*vel(3,iatom)
     if(scprod<0.0_dp .and. itime/=0)then
        stopped(iatom)=1
        !    Shift the velocities of the previous half-step and current half-step,
        !    so that the acceleration is correct but the present velocity vanishes.
        vel_prevhalf(:,iatom)=vel_prevhalf(:,iatom)-vel(:,iatom)
        vel_nexthalf(:,iatom)=vel_nexthalf(:,iatom)-vel(:,iatom)
        vel(:,iatom)=0.0_dp
        xcart_next(:,iatom)=xcart(:,iatom)+dtion*vel_nexthalf(:,iatom)
     end if
  end do

  !  Establish a list of stopped atoms
  nstopped=sum(stopped(:))

  if(nstopped/=0)then
     write(message,'(a)') ' List of stopped atoms (ionmov=7) :'
     call abi_wrtout(ab_out,message,'COLL')
     call abi_wrtout(std_out,message,'COLL')
     istopped=1
     do iatom=1,natom
        if(stopped(iatom)==1)then
           stopped(istopped)=iatom
           istopped=istopped+1
        end if
     end do
     do ii=1,nstopped,16
        write(message, '(16i4)' )stopped(ii:min(ii+15,nstopped))
        call abi_wrtout(ab_out,message,'COLL')
        call abi_wrtout(std_out,message,'COLL')
     end do
     !   End of test nstopped/=0
  end if

  deallocate(stopped)

  !   Now, compute the corrected kinetic energy
  ekin=0.0_dp
  do iatom=1,natom
     do idir=1,3
        !     Warning : the fixing of atomis is implemented in reduced
        !     coordinates, so that this expression is wrong
        if (iatfix(idir,iatom) == 0) then
           ekin=ekin+0.5_dp*amass(iatom)*vel(idir,iatom)**2
        end if
     end do
  end do

  !   Generate xred_next from xcart_next
  call xredxcart(natom,-1,rprimd,xcart_next,xred_next)

  !  End of test ionmov==7
end subroutine md_quenched_stop_atoms
