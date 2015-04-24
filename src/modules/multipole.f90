module multipole
  use module_base
  use multipole_base
  implicit none

  private

  !> Public routines
  public :: potential_from_multipoles

  contains

    !> Calculate the external potential arising from the multipoles provided
    subroutine potential_from_multipoles(ep, n1i, n2i, n3i, hx, hy, hz, shiftx, shifty, shiftz)
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: n1i, n2i, n3i
      real(dp),intent(in) :: hx, hy, hz, shiftx, shifty, shiftz

      ! Local variables
      integer :: i1, i2, i3, impl
      real(dp) :: x, xx, y, yy, z, zz, rnrm1, rnrm2, rnrm3, rnrm5
      real(dp),dimension(3) :: r


      do i3=1,n3i
          z = real(i3,kind=8)*hz
          zz = z + shiftz
          do i2=1,n2i
              y = real(i2,kind=8)*hy
              yy = y + shifty
              do i1=1,n1i
                  x = real(i1,kind=8)*hx
                  xx = x + shiftx
                  do impl=1,ep%nmpl
                      r(1) = ep%mpl(impl)%rxyz(1) - xx
                      r(2) = ep%mpl(impl)%rxyz(2) - yy
                      r(3) = ep%mpl(impl)%rxyz(3) - zz 
                      rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                      rnrm1 = sqrt(rnrm2)
                      rnrm3 = rnrm1*rnrm2
                      rnrm5 = rnrm3*rnrm2
                  end do
              end do
          end do
      end do

    end subroutine potential_from_multipoles

end module multipole
