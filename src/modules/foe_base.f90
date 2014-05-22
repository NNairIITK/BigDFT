module foe_base
  use module_defs, only: uninitialized
  use module_base
  implicit none

  private

  type,public :: foe_data
    real(kind=8),private :: ef                     !< Fermi energy for FOE
    real(kind=8),private :: evlow, evhigh          !< Eigenvalue bounds for FOE 
    real(kind=8),private :: bisection_shift        !< Bisection shift to find Fermi energy (FOE)
    real(kind=8),private :: fscale                 !< Length scale for complementary error function (FOE)
    real(kind=8),private :: ef_interpol_det        !< FOE: max determinant of cubic interpolation matrix
    real(kind=8),private :: ef_interpol_chargediff !< FOE: max charge difference for interpolation
    real(kind=8),private :: charge                 !< Total charge of the system
    real(kind=8),private :: fscale_lowerbound      !< lower bound for the error function decay length
    real(kind=8),private :: fscale_upperbound      !< upper bound for the error function decay length
    integer,private :: evbounds_isatur, evboundsshrink_isatur, evbounds_nsatur, evboundsshrink_nsatur !< variables to check whether the eigenvalue bounds might be too big

    contains
      procedure,pass :: set_int
      procedure,pass :: get_int
      procedure,pass :: set_real
      procedure,pass :: get_real
  end type foe_data


  public :: foe_data_null
  public :: set_int
  public :: get_int
  public :: set_real
  public :: get_real


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


    subroutine set_int(this, fieldname, val)
      class(foe_data) :: this
      character(len=*),intent(in) :: fieldname
      integer,intent(in) :: val

      select case (fieldname)
      case ("nseg")
          !!this%nseg = val
      case ("evbounds_isatur")
          this%evbounds_isatur = val
      case ("evboundsshrink_isatur")
          this%evboundsshrink_isatur = val
      case ("evbounds_nsatur")
          this%evbounds_nsatur = val
      case ("evboundsshrink_nsatur")
          this%evboundsshrink_nsatur = val
      case default
          stop 'wrong arguments'
      end select

    end subroutine set_int


    integer function get_int(this, fieldname) result(val)
      class(foe_data) :: this
      character(len=*),intent(in) :: fieldname
      !integer,intent(in) :: val

      select case (fieldname)
      case ("nseg")
          !!val = this%nseg
      case ("evbounds_isatur")
          val = this%evbounds_isatur
      case ("evboundsshrink_isatur")
          val = this%evboundsshrink_isatur
      case ("evbounds_nsatur")
          val = this%evbounds_nsatur
      case ("evboundsshrink_nsatur")
          val = this%evboundsshrink_nsatur
      case default
          stop 'wrong arguments'
      end select

    end function get_int


    subroutine set_real(this, fieldname, val)
      class(foe_data) :: this
      character(len=*),intent(in) :: fieldname
      real(kind=8),intent(in) :: val

      select case (fieldname)
      case ("ef")
          this%ef = val
      case ("evlow")
          this%evlow = val
      case ("evhigh")
          this%evhigh = val
      case ("bisection_shift")
          this%bisection_shift = val
      case ("fscale")
          this%fscale = val
      case ("ef_interpol_det")
          this%ef_interpol_det = val
      case ("ef_interpol_chargediff")
          this%ef_interpol_chargediff = val
      case ("charge")
          this%charge = val
      case ("fscale_lowerbound")
          this%fscale_lowerbound = val
      case ("fscale_upperbound")
          this%fscale_upperbound = val
      case default
          stop 'wrong arguments'
      end select

    end subroutine set_real


    real(kind=8) function get_real(this, fieldname) result(val)
      class(foe_data) :: this
      character(len=*),intent(in) :: fieldname

      select case (fieldname)
      case ("ef")
          val = this%ef
      case ("evlow")
          val = this%evlow
      case ("evhigh")
          val = this%evhigh
      case ("bisection_shift")
          val = this%bisection_shift
      case ("fscale")
          val = this%fscale
      case ("ef_interpol_det")
          val = this%ef_interpol_det
      case ("ef_interpol_chargediff")
          val = this%ef_interpol_chargediff
      case ("charge")
          val = this%charge
      case ("fscale_lowerbound")
          val = this%fscale_lowerbound
      case ("fscale_upperbound")
          val = this%fscale_upperbound
      case default
          stop 'wrong arguments'
      end select

    end function get_real


end module foe_base
