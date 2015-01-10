module ab6_moldyn

  use defs_basis

  implicit none

  interface
     subroutine scfloop_main(acell, epot, fcart, grad, itime, me, natom, rprimd, xred)
       use defs_basis
       integer, intent(in) :: natom, itime, me
       real(dp), intent(out) :: epot
       real(dp), intent(in) :: acell(3)
       real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
       real(dp), intent(out) :: fcart(3, natom), grad(3, natom)
     end subroutine scfloop_main
  end interface

  interface
     subroutine scfloop_output(acell, epot, ekin, fred, itime, me, natom, rprimd, vel, xred)
       use defs_basis
       integer, intent(in) :: natom, itime, me
       real(dp), intent(in) :: epot, ekin
       real(dp), intent(in) :: acell(3)
       real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
       real(dp), intent(in) :: fred(3, natom), vel(3, natom)
     end subroutine scfloop_output
  end interface

contains

  include "velocity_verlet.F90"
  include "quenched.F90"
  include "langevin.F90"
  include "nose.F90"
  include "isokinetic.F90"
  include "isotemp.F90"
  include "isothermal.F90"
  include "moldyn.F90"

end module ab6_moldyn
