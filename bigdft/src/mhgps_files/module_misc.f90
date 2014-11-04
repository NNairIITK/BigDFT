module module_misc
    implicit none

    private

    public :: almostequal

    contains

!=====================================================================
logical function almostequal( x, y, ulp )
    use module_base
    real(gp), intent(in) :: x
    real(gp), intent(in) :: y
    integer, intent(in) :: ulp
    almostequal = abs(x-y)<( real(ulp,gp)*&
                  spacing(max(abs(x),abs(y))))
end function
!=====================================================================
end module
