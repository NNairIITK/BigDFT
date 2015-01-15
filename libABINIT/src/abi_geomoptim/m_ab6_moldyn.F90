module m_ab6_moldyn

  use abi_defs_basis

  implicit none

  interface
     subroutine scfloop_main(acell, epot, fcart, grad, itime, me, natom, rprimd, xred)
       use abi_defs_basis
       integer, intent(in) :: natom, itime, me
       real(dp), intent(out) :: epot
       real(dp), intent(in) :: acell(3)
       real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
       real(dp), intent(out) :: fcart(3, natom), grad(3, natom)
     end subroutine scfloop_main
  end interface

  interface
     subroutine scfloop_output(acell, epot, ekin, fred, itime, me, natom, rprimd, vel, xred)
       use abi_defs_basis
       integer, intent(in) :: natom, itime, me
       real(dp), intent(in) :: epot, ekin
       real(dp), intent(in) :: acell(3)
       real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
       real(dp), intent(in) :: fred(3, natom), vel(3, natom)
     end subroutine scfloop_output
  end interface

! abi_mttk_type: dataype used in Martyna et al. (TTK) reversible MD integration scheme
  type abi_mttk_type
   !Real (double precision) scalars
    real(dp) :: glogv              !Logarithm of the volume
    real(dp) :: vlogv              !Derivative of logv
   !Real (double precision) arrays
    real(dp) :: gboxg(3,3)         !Imbalance in pressure (see paper)
    real(dp) :: vboxg(3,3)         !Velocity of log rprimd (see paper)
    real(dp), pointer :: glogs(:)  ! Imbalance of kinetic energy
    real(dp), pointer :: vlogs(:)  ! Velocities of thermostat variables
    real(dp), pointer :: xlogs(:)  ! Positions of thermostat variables
  end type abi_mttk_type

contains

  include "abi_xfpack.F90.inc"
  include "abi_velocity_verlet.F90.inc"
  include "abi_quenched.F90.inc"
  include "abi_langevin.F90.inc"
  include "abi_nose.F90.inc"
  include "abi_isokinetic.F90.inc"
  include "abi_isotemp.F90.inc"
  include "abi_isothermal.F90.inc"
  include "abi_moldyn.F90.inc"

end module m_ab6_moldyn
