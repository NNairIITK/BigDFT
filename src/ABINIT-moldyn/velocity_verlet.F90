subroutine md_velocity_verlet(acell, acell_next, amass, dtion, fred, &
     & hessin, iatfix, itime, natom, optcell, rprim, &
     & rprim_next, rprimd, rprimd_next, &
     & ucvol, ucvol_next, vel, vel_nexthalf, vel_prevhalf, &
     & xcart, xcart_next, xcart_prev, xred_next)
  
  use defs_basis

  implicit none

  integer, intent(in) :: natom, itime, optcell
  integer, intent(in) :: iatfix(3, natom)
  real(dp), intent(in) :: dtion, ucvol
  real(dp), intent(out) :: ucvol_next
  real(dp), intent(in) :: vel_prevhalf(3, natom), amass(natom)
  real(dp), intent(in) :: rprim(3,3), rprimd(3,3), acell(3)
  real(dp), intent(in) :: hessin(3 * natom, 3 * natom), fred(3, natom)
  real(dp), intent(in) :: xcart_prev(3, natom), xcart(3, natom)
  real(dp), intent(inout) :: vel(3,natom)
  real(dp), intent(out) :: vel_nexthalf(3, natom)
  real(dp), intent(out) :: acell_next(3), rprim_next(3,3), rprimd_next(3,3)
  real(dp), intent(out) :: xred_next(3,natom), xcart_next(3,natom)

  integer :: iatom, idir, ndim, jatom
  real(dp) :: taylor, acc(3)
  real(dp) :: gmet(3,3), gprimd(3,3), rmet(3,3)

  ndim = 3 * natom

  !  Compute next atomic coordinates and cell parameters, using Verlet algorithm
  !  First propagate the position, without acceleration
  if(itime/=0)then
     xcart_next(:,:) = 2 * xcart(:,:) - xcart_prev(:,:)
     taylor=one
  else
     !   Initialisation : no vin_prev is available, but the ionic velocity
     !   is available, in cartesian coordinates
     !   Uses the velocity
     xcart_next(:,:) = xcart(:,:) + dtion * vel(:,:)
     !   Impose no change of acell, ucvol, rprim, and rprimd
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     rprim_next(:,:)=rprim(:,:)
     rprimd_next(:,:)=rprimd(:,:)
     taylor=half
  end if

  !  Now, take into account the acceleration
  do iatom = 1, natom, 1
!!$     acc(:) = zero
!!$     do jatom = 1, natom
!!$        acc(1) = acc(1) + hessin(1 + 3 * (iatom - 1), 1 + 3 * (jatom - 1)) * &
!!$             & fred(1, jatom)
!!$        acc(2) = acc(2) + hessin(2 + 3 * (iatom - 1), 2 + 3 * (jatom - 1)) * &
!!$             & fred(2, jatom)
!!$        acc(3) = acc(3) + hessin(3 + 3 * (iatom - 1), 3 + 3 * (jatom - 1)) * &
!!$             & fred(3, jatom)
!!$     end do
     acc(:) = (rprimd(:,1) * fred(1,iatom) + &
          & rprimd(:,2) * fred(2,iatom) + &
          & rprimd(:,3) * fred(3,iatom)) / amass(iatom)
     !   Note the minus sign: the forces are minus the gradients, contained in vout.
     xcart_next(:, iatom) = xcart_next(:, iatom) - dtion**2 * acc(:) * taylor
     !  Implement fixing of atoms : put back old values for fixed components
     do idir=1,3
        if (iatfix(idir,iatom) == 1) then
           xcart_next(idir, iatom) = xcart(idir, iatom)
        end if
     end do
  end do
  ! Transfert cart into red
  call xredxcart(natom,-1,rprimd,xcart_next,xred_next)

  !  Now, compute the velocity at the next half-step
  !  Get xred_next, and eventually acell_next, ucvol_next, rprim_next and
  !  rprimd_next, from vin_next
  if(optcell/=0)then
     call mkrdim(acell_next,rprim_next,rprimd_next)
     call metric(gmet,gprimd,-1,rmet,rprimd_next,ucvol_next)
  else
     !   Impose no change of acell, ucvol, rprim, and rprimd
     acell_next(:)=acell(:)
     ucvol_next=ucvol
     rprim_next(:,:)=rprim(:,:)
     rprimd_next(:,:)=rprimd(:,:)
  end if

  !  Compute the velocity at half of the new step
  vel_nexthalf(:,:)=(xcart_next(:,:)-xcart(:,:))/dtion

  !  If needed, compute the velocity at present position
  if(itime/=0)then
     vel(:,:)=(vel_nexthalf(:,:)+vel_prevhalf(:,:))*0.5_dp
  end if

  !  End of case ionmov /=8, /=9, /=12, /=13, /=14
end subroutine md_velocity_verlet
