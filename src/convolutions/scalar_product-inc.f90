
!> Search the segments which intersect each other
!! @todo Modify this routine to have also the end as result.
pure subroutine hunt_inline(xx,n,x,jlo)
  implicit none
  integer, intent(in) :: x                !< Starting point in grid coordinates
  integer, intent(in) :: n                !< Number of segments
  integer, dimension(n), intent(in) :: xx !< Array of segment starting points
  integer, intent(inout) :: jlo           !< Input: starting segment, 
  !! Output: closest segment corresponding to x
  !! @warning if jlo is outside range, routine is disabled
  !local variables
  integer :: inc,jhi,jm

  !check array extremes
  if (jlo > n) return

  !start searching
  if(x == xx(1))then
     jlo=1
     return
  end if
  if(x == xx(n)) then 
     jlo=n
     return
  end if

  !increment of the segment
  inc=1
  !target is above starting point
  if (x >= xx(jlo)) then
     guess_end: do
        jhi=jlo+inc
        !number of segments is over
        if(jhi > n)then
           jhi=n+1
           exit guess_end
           !increase until the target is below
        else if(x >= xx(jhi))then
           jlo=jhi
           inc=inc+inc
        else
           exit guess_end
        endif
     end do guess_end
  else
     !target is below, invert start and end
     jhi=jlo
     guess_start: do
        jlo=jhi-inc
        !segment are over (from below)
        if (jlo < 1) then
           jlo=0
           exit guess_start
           !decrease until the target is above
        else if(x < xx(jlo))then
           jhi=jlo
           inc=inc+inc
        else
           exit guess_start
        endif
     end do guess_start
  endif

  binary_search: do
     !the end and the beginning are contiguous: segment number has been found
     if (jhi-jlo == 1) exit binary_search
     !mean point (integer division, rounded towards jhi)
     jm=(jhi+jlo)/2
     !restrict search from the bottom of from the top
     if (x >= xx(jm)) then
        jlo=jm
     else
        jhi=jm
     endif
  end do binary_search

END SUBROUTINE hunt_inline
