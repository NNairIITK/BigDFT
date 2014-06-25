module fermi_level
  use module_base
  implicit none

  private

  ! Public routines
  public :: init_fermi_level
  public :: determine_fermi_level
  public :: fermilevel_get_real
  public :: fermilevel_get_logical

  ! Variables local to the module
  logical :: adjust_lower_bound_, adjust_upper_bound_
  real(kind=8) :: target_charge_, bisection_shift_, ef_old_, sumn_old_
  real(kind=8) :: ef_interpol_chargediff_, ef_interpol_det_
  real(kind=8),dimension(2) :: sumnarr_, efarr_
  logical,dimension(2) :: bisection_bounds_ok_
  real(kind=8),dimension(4,4) :: interpol_matrix_
  real(kind=8),dimension(4) :: interpol_vector_
  integer :: it_, it_solver_, verbosity_


  contains
    
    !> Initialize the internal variables
    subroutine init_fermi_level(target_charge, ef, bisection_shift, ef_interpol_chargediff, ef_interpol_det, verbosity)
      implicit none

      ! Calling arguments
      real(kind=8),intent(in) :: target_charge                   !< total charge of the system
      real(kind=8),intent(in) :: ef                              !< initial guess for the fermi level
      real(kind=8),intent(in),optional :: bisection_shift        !< shift to be used for the determination of the bisection bounds
      real(kind=8),intent(in),optional :: ef_interpol_chargediff !< charge difference below which the cubic interpolation is allowed
      real(kind=8),intent(in),optional :: ef_interpol_det        !< determinant of the interpolation matrix above which the cubic interpolation is allowed
      integer,intent(in),optional :: verbosity                   !< verbosity of the output: 0 for no output, 1 for more detailed output

      adjust_lower_bound_ = .true.
      adjust_upper_bound_ = .true.
      target_charge_ = target_charge
      sumnarr_(1) = 0.d0
      sumnarr_(2) = 1.d100
      if (present(bisection_shift)) then
          bisection_shift_ = bisection_shift
      else
          bisection_shift_ = 0.1d0
      end if
      ef_old_ = 0.d0
      sumn_old_ = 0.d0
      efarr_(1) = ef - bisection_shift_
      efarr_(2) = ef + bisection_shift_
      if (present(ef_interpol_chargediff)) then
          ef_interpol_chargediff_ = ef_interpol_chargediff
      else
          ef_interpol_chargediff_ = 1.d0
      end if
      if (present(ef_interpol_det)) then
          ef_interpol_det_ = ef_interpol_det
      else
          ef_interpol_det_ = 1.d-20
      end if
      bisection_bounds_ok_(1) = .false.
      bisection_bounds_ok_(2) = .false.
      interpol_matrix_(:,:) = 0.d0
      interpol_vector_(:) = 0.d0
      it_ = 0
      it_solver_ = 0
      if (present(verbosity)) then
          verbosity_ = verbosity
      else
          verbosity_ = 0
      end if
    end subroutine init_fermi_level



    subroutine determine_fermi_level(sumn, ef, info)
      use yaml_output
      implicit none

      ! Calling arguments
      real(kind=8),intent(in) :: sumn      !< charge of the system (which should be equal to the target charge once the correct Fermi level is found),
                                           !    obtained with the current value of ef
      real(kind=8),intent(inout) :: ef     !< on input: current value of the Fermi level
                                           !  on output: new guess for the Fermi level, depending on the value of info
      integer,intent(out),optional :: info !< info parameter: * -1: adjusting the lower bisection bound, ef not meaningful
                                           !                  * -2: adjusting the upper bisection bound, ef not meaningful
                                           !                  *  0: searching the correct fermi level, ef meaningful
      ! Local variables
      real(kind=8) :: charge_diff
      logical :: interpolation_possible
      integer :: internal_info


      ! Make sure that the bounds for the bisection are negative and positive
      charge_diff = sumn-target_charge_
      if (adjust_lower_bound_) then
          if (charge_diff <= 0.d0) then
              ! Lower bound okay
              adjust_lower_bound_ = .false.
              bisection_shift_ = bisection_shift_*0.9d0
              sumnarr_(1) = sumn
              bisection_bounds_ok_(1) = .true.
          else
              efarr_(1) = efarr_(1)-bisection_shift_
              bisection_shift_ = bisection_shift_*1.1d0
              bisection_bounds_ok_(1) = .false.
          end if
      else if (adjust_upper_bound_) then
          if (charge_diff >= 0.d0) then
              ! Upper bound okay
              adjust_upper_bound_ = .false.
              bisection_shift_ = bisection_shift_*0.9d0
              sumnarr_(2) = sumn
              bisection_bounds_ok_(2) = .true.
          else
              efarr_(2) = efarr_(2)+bisection_shift_
              bisection_shift_ = bisection_shift_*1.1d0
              bisection_bounds_ok_(2)=.false.
          end if
      end if


      internal_info = 0
      if (adjust_lower_bound_) then
          ef = efarr_(1)
          internal_info = -1
      else if (adjust_upper_bound_) then
          ef = efarr_(2)
          internal_info = -2
      end if
      if (present(info)) then
          info = internal_info
      end if

      if (internal_info < 0) return ! no need to proceed further

      ! If we made it here, the bounds are ok (i.e. the lower bound gives a negative charge difference and the upper one a positive one).


      ! Adjust the bounds for the bisection, i.e. make the search interval more narrow
      if (charge_diff < 0.d0) then
          efarr_(1) = ef
          sumnarr_(1) = sumn
      else if (charge_diff >= 0.d0) then
          efarr_(2) = ef
          sumnarr_(2) = sumn
      end if


      it_solver_ = it_solver_+1

      ! Check whether the system behaves reasonably.
      interpolation_possible=.true.
      if (it_solver_ > 1) then
          if (verbosity_ >= 1 .and. bigdft_mpi%iproc==0) then
              call yaml_newline()
              call yaml_open_map('interpol check',flow=.true.)
              call yaml_map('D eF',ef-ef_old_,fmt='(es13.6)')
              call yaml_map('D Tr',sumn-sumn_old_,fmt='(es13.6)')
          end if
          if (ef > ef_old_ .and. sumn < sumn_old_) then
              interpolation_possible = .false.
          end if
          if (ef < ef_old_ .and. sumn > sumn_old_) then
              interpolation_possible = .false.
          end if
          if (interpolation_possible) then
              if (verbosity_>=1 .and. bigdft_mpi%iproc==0) call yaml_map('interpol possible',.true.)
          else
              if (verbosity_>=1 .and. bigdft_mpi%iproc==0) call yaml_map('interpol possible',.false.)
          end if
          if (verbosity_>=1 .and. bigdft_mpi%iproc==0) call yaml_close_map()
          if (verbosity_>=1 .and. bigdft_mpi%iproc==0) call yaml_newline()
      end if
      if (.not.interpolation_possible) then
          ! Set the history for the interpolation to zero.
          it_solver_=0
      end if

      ef_old_ = ef
      sumn_old_ = sumn

      call determine_new_fermi_level()


      contains

        subroutine determine_new_fermi_level()
          implicit none
          integer :: info, i, ii
          real(kind=8) :: determinant, m, b, ef_interpol, det
          real(kind=8),dimension(4,4) :: tmp_matrix
          real(kind=8),dimension(4) :: interpol_solution
          integer,dimension(4) :: ipiv

          ! Shift up the old results.
          if (it_solver_>4) then
              do i=1,4
                  interpol_matrix_(1,i)=interpol_matrix_(2,i)
                  interpol_matrix_(2,i)=interpol_matrix_(3,i)
                  interpol_matrix_(3,i)=interpol_matrix_(4,i)
              end do
              interpol_vector_(1)=interpol_vector_(2)
              interpol_vector_(2)=interpol_vector_(3)
              interpol_vector_(3)=interpol_vector_(4)
          end if
          !LG: if it_solver_==0 this index comes out of bounds!
          ii=max(min(it_solver_,4),1)
          interpol_matrix_(ii,1)=ef**3
          interpol_matrix_(ii,2)=ef**2
          interpol_matrix_(ii,3)=ef
          interpol_matrix_(ii,4)=1
          interpol_vector_(ii)=sumn-target_charge_
        
          ! Solve the linear system interpol_matrix_*interpol_solution=interpol_vector_
          if (it_solver_>=4) then
              do i=1,ii
                  interpol_solution(i)=interpol_vector_(i)
                  tmp_matrix(i,1)=interpol_matrix_(i,1)
                  tmp_matrix(i,2)=interpol_matrix_(i,2)
                  tmp_matrix(i,3)=interpol_matrix_(i,3)
                  tmp_matrix(i,4)=interpol_matrix_(i,4)
              end do
        
              call dgesv(ii, 1, tmp_matrix, 4, ipiv, interpol_solution, 4, info)
              if (info/=0) then
                 if (bigdft_mpi%iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
              end if
        
        
              call get_roots_of_cubic_polynomial(interpol_solution(1), interpol_solution(2), &
                   interpol_solution(3), interpol_solution(4), ef, ef_interpol)
          end if
        
        
        
        
          ! Calculate the new Fermi energy.
          if (verbosity_>=1 .and. bigdft_mpi%iproc==0) then
              call yaml_newline()
              call yaml_open_map('Search new eF',flow=.true.)
          end if
          if (it_solver_>=4 .and.  &
              abs(sumn-target_charge_) < ef_interpol_chargediff_) then
              det=determinant(bigdft_mpi%iproc,4,interpol_matrix_)
              if (verbosity_ >= 1 .and. bigdft_mpi%iproc==0) then
                  call yaml_map('det',det,fmt='(es10.3)')
                  call yaml_map('limit',ef_interpol_det_,fmt='(es10.3)')
              end if
              if(abs(det) > ef_interpol_det_) then
                  ef = ef_interpol
                  if (verbosity_>=1 .and. bigdft_mpi%iproc==0) call yaml_map('method','cubic interpolation')
              else
                  ! linear interpolation
                  m = (interpol_vector_(4)-interpol_vector_(3))/(interpol_matrix_(4,3)-interpol_matrix_(3,3))
                  b = interpol_vector_(4)-m*interpol_matrix_(4,3)
                  ef = -b/m
              end if
          else
              ! Use mean value of bisection and secant method
              ! Secant method solution
              ef = efarr_(2)-(sumnarr_(2)-target_charge_)*(efarr_(2)-efarr_(1))/(sumnarr_(2)-sumnarr_(1))
              ! Add bisection solution
              ef = ef + 0.5d0*(efarr_(1)+efarr_(2))
              ! Take the mean value
              ef = 0.5d0*ef
              if (verbosity_>=1 .and. bigdft_mpi%iproc==0) call yaml_map('method','bisection / secant method')
          end if
          if (verbosity_>=1 .and. bigdft_mpi%iproc==0) then
              call yaml_close_map()
          end if

        end subroutine determine_new_fermi_level


    end subroutine determine_fermi_level


    function fermilevel_get_real(fieldname) result(val)
        ! Calling arguments
        character(len=*),intent(in) :: fieldname
        real(kind=8) :: val

        select case (trim(fieldname))
        case ("efarr(1)")
            val = efarr_(1)
        case ("efarr(2)")
            val = efarr_(2)
        case ("bisection_shift")
            val = bisection_shift_
        case default
            stop 'ERROR: wrong argument'
        end select
    end function fermilevel_get_real

    function fermilevel_get_logical(fieldname) result(val)
        ! Calling arguments
        character(len=*),intent(in) :: fieldname
        logical :: val

        select case (trim(fieldname))
        case ("bisection_bounds_ok(1)")
            val = bisection_bounds_ok_(1)
        case ("bisection_bounds_ok(2)")
            val = bisection_bounds_ok_(2)
        case default
            stop 'ERROR: wrong argument'
        end select
    end function fermilevel_get_logical


end module fermi_level
