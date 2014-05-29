module foe_base
  use module_defs, only: uninitialized
  use module_base
  implicit none

  private

  type,public :: foe_data
    real(kind=8) :: ef                     !< Fermi energy for FOE
    real(kind=8) :: evlow, evhigh          !< Eigenvalue bounds for FOE 
    real(kind=8) :: bisection_shift        !< Bisection shift to find Fermi energy (FOE)
    real(kind=8) :: fscale                 !< Length scale for complementary error function (FOE)
    real(kind=8) :: ef_interpol_det        !< FOE: max determinant of cubic interpolation matrix
    real(kind=8) :: ef_interpol_chargediff !< FOE: max charge difference for interpolation
    real(kind=8) :: charge                 !< Total charge of the system
    real(kind=8) :: fscale_lowerbound      !< lower bound for the error function decay length
    real(kind=8) :: fscale_upperbound      !< upper bound for the error function decay length
    integer :: evbounds_isatur, evboundsshrink_isatur, evbounds_nsatur, evboundsshrink_nsatur !< variables to check whether the eigenvalue bounds might be too big
  end type foe_data


  public :: foe_data_null
  public :: foe_data_set_int
  public :: foe_data_get_int
  public :: foe_data_set_real
  public :: foe_data_get_real


  contains
 

    function foe_data_null() result(foe_obj)
      implicit none
      type(foe_data) :: foe_obj
      foe_obj%ef                     = uninitialized(foe_obj%ef)
      foe_obj%evlow                  = uninitialized(foe_obj%evlow)
      foe_obj%bisection_shift        = uninitialized(foe_obj%bisection_shift)
      foe_obj%fscale                 = uninitialized(foe_obj%fscale)
      foe_obj%ef_interpol_det        = uninitialized(foe_obj%ef_interpol_det)
      foe_obj%ef_interpol_chargediff = uninitialized(foe_obj%ef_interpol_chargediff)
      foe_obj%charge                 = uninitialized(foe_obj%charge)
      foe_obj%fscale_lowerbound      = uninitialized(foe_obj%fscale_lowerbound)
      foe_obj%fscale_upperbound      = uninitialized(foe_obj%fscale_upperbound)
      foe_obj%evbounds_isatur        = uninitialized(foe_obj%evbounds_isatur)
      foe_obj%evboundsshrink_isatur  = uninitialized(foe_obj%evboundsshrink_isatur)
      foe_obj%evbounds_nsatur        = uninitialized(foe_obj%evbounds_nsatur)
      foe_obj%evboundsshrink_nsatur  = uninitialized(foe_obj%evboundsshrink_nsatur)
    end function foe_data_null


    subroutine foe_data_set_int(foe_obj, fieldname, val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      integer,intent(in) :: val

      select case (fieldname)
      case ("nseg")
          !!foe_obj%nseg = val
      case ("evbounds_isatur")
          foe_obj%evbounds_isatur = val
      case ("evboundsshrink_isatur")
          foe_obj%evboundsshrink_isatur = val
      case ("evbounds_nsatur")
          foe_obj%evbounds_nsatur = val
      case ("evboundsshrink_nsatur")
          foe_obj%evboundsshrink_nsatur = val
      case default
          stop 'wrong arguments'
      end select

    end subroutine foe_data_set_int


    integer function foe_data_get_int(foe_obj, fieldname) result(val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      !integer,intent(in) :: val

      select case (fieldname)
      case ("nseg")
          !!val = foe_obj%nseg
      case ("evbounds_isatur")
          val = foe_obj%evbounds_isatur
      case ("evboundsshrink_isatur")
          val = foe_obj%evboundsshrink_isatur
      case ("evbounds_nsatur")
          val = foe_obj%evbounds_nsatur
      case ("evboundsshrink_nsatur")
          val = foe_obj%evboundsshrink_nsatur
      case default
          stop 'wrong arguments'
      end select

    end function foe_data_get_int


    subroutine foe_data_set_real(foe_obj, fieldname, val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      real(kind=8),intent(in) :: val

      select case (fieldname)
      case ("ef")
          foe_obj%ef = val
      case ("evlow")
          foe_obj%evlow = val
      case ("evhigh")
          foe_obj%evhigh = val
      case ("bisection_shift")
          foe_obj%bisection_shift = val
      case ("fscale")
          foe_obj%fscale = val
      case ("ef_interpol_det")
          foe_obj%ef_interpol_det = val
      case ("ef_interpol_chargediff")
          foe_obj%ef_interpol_chargediff = val
      case ("charge")
          foe_obj%charge = val
      case ("fscale_lowerbound")
          foe_obj%fscale_lowerbound = val
      case ("fscale_upperbound")
          foe_obj%fscale_upperbound = val
      case default
          stop 'wrong arguments'
      end select

    end subroutine foe_data_set_real


    real(kind=8) function foe_data_get_real(foe_obj, fieldname) result(val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname

      select case (fieldname)
      case ("ef")
          val = foe_obj%ef
      case ("evlow")
          val = foe_obj%evlow
      case ("evhigh")
          val = foe_obj%evhigh
      case ("bisection_shift")
          val = foe_obj%bisection_shift
      case ("fscale")
          val = foe_obj%fscale
      case ("ef_interpol_det")
          val = foe_obj%ef_interpol_det
      case ("ef_interpol_chargediff")
          val = foe_obj%ef_interpol_chargediff
      case ("charge")
          val = foe_obj%charge
      case ("fscale_lowerbound")
          val = foe_obj%fscale_lowerbound
      case ("fscale_upperbound")
          val = foe_obj%fscale_upperbound
      case default
          stop 'wrong arguments'
      end select

    end function foe_data_get_real


end module foe_base
