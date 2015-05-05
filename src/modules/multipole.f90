module multipole
  use module_base
  use multipole_base
  implicit none

  private

  !> Public routines
  public :: potential_from_multipoles

  contains

    !> Calculate the external potential arising from the multipoles provided
    subroutine potential_from_multipoles(ep, is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, pot)
      implicit none
      
      ! Calling arguments
      type(external_potential_descriptors),intent(in) :: ep
      integer,intent(in) :: is1, ie1, is2, ie2, is3, ie3
      real(dp),intent(in) :: hx, hy, hz
      real(dp),dimension(is1:ie1,is2:ie2,is3:ie3),intent(inout) :: pot

      ! Local variables
      integer :: i1, i2, i3, impl, l
      real(dp) :: x, y, z, rnrm1, rnrm2, rnrm3, rnrm5, mp
      real(dp),dimension(3) :: r

      !$omp parallel &
      !$omp default(none) &
      !$omp shared(is1, ie1, is2, ie2, is3, ie3, hx, hy, hz, ep, pot) &
      !$omp private(i1, i2, i3, x, y, z, impl, r, rnrm1, rnrm2, rnrm3, rnrm5, l, mp)
      !$omp do
      do i3=is3,ie3
          z = real(i3,kind=8)*hz
          do i2=is2,ie2
              y = real(i2,kind=8)*hy
              do i1=is1,ie1
                  x = real(i1,kind=8)*hx
                  do impl=1,ep%nmpl
                      r(1) = ep%mpl(impl)%rxyz(1) - x
                      r(2) = ep%mpl(impl)%rxyz(2) - y
                      r(3) = ep%mpl(impl)%rxyz(3) - z 
                      rnrm2 = r(1)**2 + r(2)**2 + r(3)**2
                      rnrm1 = sqrt(rnrm2)
                      rnrm3 = rnrm1*rnrm2
                      rnrm5 = rnrm3*rnrm2
                      mp = 0.0_dp
                      do l=0,lmax-1
                          if (associated(ep%mpl(impl)%qlm(l)%q)) then
                              select case(l)
                              case (0)
                                  mp = mp + calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                                  !write(*,'(a,3es12.4,es16.8)') 'x, y, z, monopole', x, y, z, calc_monopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                              case (1)
                                  mp = mp + calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                                  !write(*,*) 'dipole', calc_dipole(ep%mpl(impl)%qlm(l)%q, r, rnrm3)
                              case (2)
                                  mp = mp + calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                                  !write(*,*) 'quadrupole', calc_quadropole(ep%mpl(impl)%qlm(l)%q, r, rnrm5)
                              case (3)
                                  call f_err_throw('octupole not yet implemented', err_name='BIGDFT_RUNTIME_ERROR')
                                  !multipole_terms(l) = calc_octopole(ep%mpl(impl)%qlm(l)%q, rnrm1)
                              case default
                                  call f_err_throw('Wrong value of l', err_name='BIGDFT_RUNTIME_ERROR')
                              end select
                          end if
                      end do
                      pot(i1,i2,i3) = pot(i1,i2,i3) + mp
                  end do
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel


      contains


        function calc_monopole(q, rnrm1) result(mpm)
          implicit none
          ! Calling arguments
          real(dp),dimension(1),intent(in) :: q
          real(dp),intent(in) :: rnrm1
          real(dp) :: mpm

          mpm = -q(1)/rnrm1

        end function calc_monopole


        function calc_dipole(q, r, rnrm3) result(dpm)
          implicit none
          ! Calling arguments
          real(dp),dimension(3),intent(in) :: q
          real(dp),intent(in) :: rnrm3
          real(dp),dimension(3),intent(in) :: r
          real(dp) :: dpm

          dpm = q(1)*r(1) + q(2)*r(2) + q(3)*r(3)
          dpm = -dpm/rnrm3

        end function calc_dipole


        function calc_quadropole(q, r, rnrm5) result(qpm)
          implicit none
          ! Calling arguments
          real(dp),dimension(5),intent(in) :: q
          real(dp),intent(in) :: rnrm5
          real(dp),dimension(3),intent(in) :: r
          real(dp) :: qpm
          ! Local variables
          real(dp),dimension(3,3) :: qq

          qq(1,1) = q(1)
          qq(2,1) = q(2)
          qq(3,1) = q(3)
          qq(1,2) = qq(2,1)
          qq(2,2) = q(4)
          qq(3,2) = q(5)
          qq(1,3) = qq(3,1)
          qq(2,3) = qq(3,2)
          qq(3,3) = 1.0_dp-qq(1,1)-qq(2,2)

          qpm = qq(1,1)*r(1)*r(1) + &
               qq(2,1)*r(2)*r(1) + &
               qq(3,1)*r(3)*r(1) + &
               qq(1,2)*r(1)*r(2) + &
               qq(2,2)*r(2)*r(2) + &
               qq(3,2)*r(3)*r(2) + &
               qq(1,3)*r(1)*r(3) + &
               qq(2,3)*r(2)*r(3) + &
               qq(3,3)*r(3)*r(3)
          qpm = -0.5_dp*qpm/rnrm5

        end function calc_quadropole


        !function calc_octopole()
        !end function calc_octopole

    end subroutine potential_from_multipoles

end module multipole
