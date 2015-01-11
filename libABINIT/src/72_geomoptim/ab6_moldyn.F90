module ab6_moldyn

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

! mttk_type: dataype used in Martyna et al. (TTK) reversible MD integration scheme
  type mttk_type
   !Real (double precision) scalars
    real(dp) :: glogv              !Logarithm of the volume
    real(dp) :: vlogv              !Derivative of logv
   !Real (double precision) arrays
    real(dp) :: gboxg(3,3)         !Imbalance in pressure (see paper)
    real(dp) :: vboxg(3,3)         !Velocity of log rprimd (see paper)
    real(dp), pointer :: glogs(:)  ! Imbalance of kinetic energy
    real(dp), pointer :: vlogs(:)  ! Velocities of thermostat variables
    real(dp), pointer :: xlogs(:)  ! Positions of thermostat variables
  end type mttk_type

contains

  include "xfpack.F90"
  include "velocity_verlet.F90"
  include "quenched.F90"
  include "langevin.F90"
  include "nose.F90"
  include "isokinetic.F90"
  include "isotemp.F90"
  include "isothermal.F90"
  include "moldyn.F90"

end module ab6_moldyn
