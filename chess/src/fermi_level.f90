!> @file
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2014-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Determination of the Fermi level for the density matrix
module fermi_level
  use sparsematrix_base
  implicit none

  private

  ! Public routines
  public :: init_fermi_level
  public :: determine_fermi_level
  public :: fermilevel_get_real
  public :: fermilevel_get_logical
  !!public :: get_roots_of_cubic_polynomial
  !!public :: determinant

  ! Auxiliary structure that holds the required data
  type,public :: fermi_aux
    logical :: adjust_lower_bound, adjust_upper_bound
    real(kind=mp) :: target_charge, bisection_shift, ef_old, sumn_old
    real(kind=mp) :: ef_interpol_chargediff, ef_interpol_det
    real(kind=mp),dimension(2) :: sumnarr, efarr
    logical,dimension(2) :: bisection_bounds_ok
    real(kind=mp),dimension(4,4) :: interpol_matrix
    real(kind=mp),dimension(4) :: interpol_vector
    integer :: it, it_solver, verbosity
  end type fermi_aux


  contains
    
    !> Initialize the internal variables
    subroutine init_fermi_level(target_charge, ef, f, bisection_shift, ef_interpol_chargediff, ef_interpol_det, verbosity)
      implicit none

      ! Calling arguments
      real(kind=mp),intent(in) :: target_charge                   !< total charge of the system
      real(kind=mp),intent(in) :: ef                              !< initial guess for the fermi level
      type(fermi_aux),intent(out) :: f                           !< type that holds the internal data
      real(kind=mp),intent(in),optional :: bisection_shift        !< shift to be used for the determination of the bisection bounds
      real(kind=mp),intent(in),optional :: ef_interpol_chargediff !< charge difference below which the cubic interpolation is allowed
      real(kind=mp),intent(in),optional :: ef_interpol_det        !< determinant of the interpolation matrix above which the cubic interpolation is allowed
      integer,intent(in),optional :: verbosity                   !< verbosity of the output: 0 for no output, 1 for more detailed output

      call f_routine(id='init_fermi_level')

      f%adjust_lower_bound = .true.
      f%adjust_upper_bound = .true.
      f%target_charge = target_charge
      f%sumnarr(1) = 0.d0
      f%sumnarr(2) = 1.d100
      if (present(bisection_shift)) then
          f%bisection_shift = bisection_shift
      else
          f%bisection_shift = 0.1d0
      end if
      f%ef_old = 0.d0
      f%sumn_old = 0.d0
      f%efarr(1) = ef - f%bisection_shift
      f%efarr(2) = ef + f%bisection_shift
      if (present(ef_interpol_chargediff)) then
          f%ef_interpol_chargediff = ef_interpol_chargediff
      else
          f%ef_interpol_chargediff = 1.d0
      end if
      if (present(ef_interpol_det)) then
          f%ef_interpol_det = ef_interpol_det
      else
          f%ef_interpol_det = 1.d-20
      end if
      f%bisection_bounds_ok(1) = .false.
      f%bisection_bounds_ok(2) = .false.
      f%interpol_matrix(:,:) = 0.d0
      f%interpol_vector(:) = 0.d0
      f%it = 0
      f%it_solver = 0
      if (present(verbosity)) then
          f%verbosity = verbosity
      else
          f%verbosity = 0
      end if

      call f_release_routine()

    end subroutine init_fermi_level



    subroutine determine_fermi_level(iproc, f, sumn, ef, info)
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc          !< task ID
      type(fermi_aux),intent(inout) :: f   !< type that holds the internal data
      real(kind=mp),intent(in) :: sumn      !< charge of the system (which should be equal to the target charge once the correct Fermi level is found),
                                           !    obtained with the current value of ef
      real(kind=mp),intent(inout) :: ef     !< on input: current value of the Fermi level
                                           !  on output: new guess for the Fermi level, depending on the value of info
      integer,intent(out),optional :: info !< info parameter: * -1: adjusting the lower bisection bound, ef not meaningful
                                           !                  * -2: adjusting the upper bisection bound, ef not meaningful
                                           !                  *  0: searching the correct fermi level, ef meaningful
      ! Local variables
      real(kind=mp) :: charge_diff
      logical :: interpolation_possible
      integer :: internal_info


!      integer :: iproc,ierr

      call f_routine(id='determine_fermi_level')

      ! Make sure that the bounds for the bisection are negative and positive
      charge_diff = sumn-f%target_charge
!iproc =mpirank(bigdft_mpi%mpi_comm)
!      call mpi_comm_rank(mpi_comm_world, iproc, ierr)
      if (f%adjust_lower_bound) then
          if (charge_diff <= 0.d0) then
              ! Lower bound okay
              f%adjust_lower_bound = .false.
              f%bisection_shift = f%bisection_shift*0.9d0
              f%sumnarr(1) = sumn
              f%bisection_bounds_ok(1) = .true.
          else
              f%efarr(1) = f%efarr(1)-f%bisection_shift
              f%bisection_shift = f%bisection_shift*1.1d0
              f%bisection_bounds_ok(1) = .false.
          end if
      else if (f%adjust_upper_bound) then
          if (charge_diff >= 0.d0) then
              ! Upper bound okay
              f%adjust_upper_bound = .false.
              f%bisection_shift = f%bisection_shift*0.9d0
              f%sumnarr(2) = sumn
              f%bisection_bounds_ok(2) = .true.
          else
              f%efarr(2) = f%efarr(2)+f%bisection_shift
              f%bisection_shift = f%bisection_shift*1.1d0
              f%bisection_bounds_ok(2)=.false.
          end if
      end if


      internal_info = 0
      if (f%adjust_lower_bound) then
          ef = f%efarr(1)
          internal_info = -1
      else if (f%adjust_upper_bound) then
          ef = f%efarr(2)
          internal_info = -2
      end if
      if (present(info)) then
          info = internal_info
      end if

      if (internal_info < 0) then
          call f_release_routine()
          if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','bisec bounds')
          return ! no need to proceed further
      end if

      ! If we made it here, the bounds are ok (i.e. the lower bound gives a negative charge difference and the upper one a positive one).


      ! Adjust the bounds for the bisection, i.e. make the search interval more narrow
      if (charge_diff < 0.d0) then
          f%efarr(1) = ef
          f%sumnarr(1) = sumn
      else if (charge_diff >= 0.d0) then
          f%efarr(2) = ef
          f%sumnarr(2) = sumn
      end if


      f%it_solver = f%it_solver+1

      ! Check whether the system behaves reasonably.
      interpolation_possible=.true.
      if (f%it_solver > 1) then
          if (ef > f%ef_old .and. sumn < f%sumn_old) then
              interpolation_possible = .false.
          else if (ef < f%ef_old .and. sumn > f%sumn_old) then
              interpolation_possible = .false.
          end if
          if (abs(sumn-f%sumn_old)<1.d-10) then
              interpolation_possible = .false.
          end if
          if (f%verbosity >= 2 .and. iproc==0) then
              call yaml_newline()
              call yaml_mapping_open('interpol check',flow=.true.)
                 call yaml_map('D eF',ef-f%ef_old,fmt='(es13.6)')
                 call yaml_map('D Tr',sumn-f%sumn_old,fmt='(es13.6)')
                 call yaml_map('interpol possible',interpolation_possible)
              call yaml_mapping_close()
              call yaml_newline()
           end if
      end if
      if (.not.interpolation_possible) then
          ! Set the history for the interpolation to zero.
          f%it_solver=0
      end if

      f%ef_old = ef
      f%sumn_old = sumn



      call determine_new_fermi_level()

      call f_release_routine()

      contains

        subroutine determine_new_fermi_level()
          implicit none
          integer :: info, i, ii
          real(kind=mp) :: m, b, ef_interpol, det
          real(kind=mp),dimension(4,4) :: tmp_matrix
          real(kind=mp),dimension(4) :: interpol_solution
          integer,dimension(4) :: ipiv
          logical :: interpolation_nonsense, cubicinterpol_possible

          !call yaml_map('sumn',sumn)

          ! Shift up the old results.
          if (f%it_solver>4) then
              do i=1,4
                  f%interpol_matrix(1,i)=f%interpol_matrix(2,i)
                  f%interpol_matrix(2,i)=f%interpol_matrix(3,i)
                  f%interpol_matrix(3,i)=f%interpol_matrix(4,i)
              end do
              f%interpol_vector(1)=f%interpol_vector(2)
              f%interpol_vector(2)=f%interpol_vector(3)
              f%interpol_vector(3)=f%interpol_vector(4)
          end if
          !LG: if f%it_solver==0 this index comes out of bounds!
          ii=max(min(f%it_solver,4),1)
          f%interpol_matrix(ii,1)=ef**3
          f%interpol_matrix(ii,2)=ef**2
          f%interpol_matrix(ii,3)=ef
          f%interpol_matrix(ii,4)=1.d0
          f%interpol_vector(ii)=sumn-f%target_charge
        
          ! Solve the linear system f%interpol_matrix*interpol_solution=f%interpol_vector
          if (f%it_solver>=4) then
              ! Calculate the determinant of the matrix used for the interpolation
              det=determinant(iproc,4,f%interpol_matrix)
              if (abs(det) > f%ef_interpol_det) then
                  cubicinterpol_possible = .true.
                  do i=1,ii
                      interpol_solution(i)=f%interpol_vector(i)
                      tmp_matrix(i,1)=f%interpol_matrix(i,1)
                      tmp_matrix(i,2)=f%interpol_matrix(i,2)
                      tmp_matrix(i,3)=f%interpol_matrix(i,3)
                      tmp_matrix(i,4)=f%interpol_matrix(i,4)
                  end do
                  !if (bigdft_mpi%iproc==0) then
                     !call yaml_map('matrix',tmp_matrix,fmt='(es10.3)')
                     !call yaml_map('interpol_vector',f%interpol_vector,fmt='(es12.5)')
                     !call yaml_newline()
                     !call yaml_map('solution',interpol_solution,fmt='(es10.3)')
                     !call yaml_map('determinant',determinant(bigdft_mpi%iproc,4,f%interpol_matrix),fmt='(es10.3)')
                  call dgesv(ii, 1, tmp_matrix, 4, ipiv, interpol_solution, 4, info)
                  if (info/=0) then
                     if (iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
                  end if
        
                  !if (bigdft_mpi%iproc==0) call yaml_map('a x^3+b x^2 + c x + d',interpol_solution,fmt='(es10.3)')
                  call get_roots_of_cubic_polynomial(interpol_solution(1), interpol_solution(2), &
                       interpol_solution(3), interpol_solution(4), ef, ef_interpol)
                  !if (bigdft_mpi%iproc==0) then
                  !    call yaml_newline()
                  !    call yaml_map('zero of cubic polynomial',ef_interpol,fmt='(es10.3)')
                  !end if
                  ! Sanity check: If the charge was too small, then new new guess
                  ! for the Fermi energy must be larger than the actual value, and
                  ! analogously if the charge was too large.
                  interpolation_nonsense = .false.
                  if (f%interpol_vector(ii)<0) then
                      ! Charge too small, the new guess must be larger
                      if (ef_interpol<ef) interpolation_nonsense = .true.
                  else if (f%interpol_vector(ii)>0) then
                      ! Charge too large, the new guess must be smaller
                      if (ef_interpol>ef) interpolation_nonsense = .true.
                  end if
              else
                  cubicinterpol_possible = .false.
              end if
          end if
        
          ! Calculate the new Fermi energy.
          if (f%verbosity>=2 .and. iproc==0) then
              call yaml_newline()
              call yaml_mapping_open('Search new eF',flow=.true.)
          end if
          if (f%it_solver>=4 .and.  &
              abs(sumn-f%target_charge) < f%ef_interpol_chargediff) then! .and. &
              !.not.interpolation_nonsense) then
              !det=determinant(bigdft_mpi%iproc,4,f%interpol_matrix)
              if (f%verbosity >= 2 .and. iproc==0) then
                  call yaml_map('det',det,fmt='(es10.3)')
                  call yaml_map('limit',f%ef_interpol_det,fmt='(es10.3)')
              end if
              !if(abs(det) > f%ef_interpol_det) then
              if(cubicinterpol_possible .and. .not.interpolation_nonsense) then
                  ef = ef_interpol
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','cubic interpol')
              else
                  ! linear interpolation
                  m = (f%interpol_vector(4)-f%interpol_vector(3))/(f%interpol_matrix(4,3)-f%interpol_matrix(3,3))
                  b = f%interpol_vector(4)-m*f%interpol_matrix(4,3)
                  ef = -b/m
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','linear interpol')
              end if
          else
              ! Use mean value of bisection and secant method if possible,
              ! otherwise only the bisection.
              ! Bisection solution
              ef = 0.5d0*(f%efarr(1)+f%efarr(2))
              !write(*,'(a,i6,5es16.7)') 'iproc, f%efarr, f%sumnarr, f%target_charge', &
              !     bigdft_mpi%iproc, f%efarr, f%sumnarr, f%target_charge
              if (abs(f%sumnarr(2)-f%sumnarr(1))>1.d-6) then !otherwise secant method numerically unstable
                  ! Add Secant method solution
                  ef = ef + f%efarr(2)-(f%sumnarr(2)-f%target_charge)*(f%efarr(2)-f%efarr(1))/(f%sumnarr(2)-f%sumnarr(1))
                  ! Take the mean value
                  ef = 0.5d0*ef
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','bisection/secant')
              else
                  if (f%verbosity>=1 .and. iproc==0) call yaml_map('new eF','bisection')
              end if
          end if
          if (f%verbosity>=2 .and. iproc==0) then
              !call yaml_map('guess for new ef',ef,fmt='(es15.8)')
              call yaml_mapping_close()
          end if


        end subroutine determine_new_fermi_level


    end subroutine determine_fermi_level


    function fermilevel_get_real(f, fieldname) result(val)
        ! Calling arguments
        type(fermi_aux),intent(in) :: f      !< type that holds the internal data
        character(len=*),intent(in) :: fieldname
        real(kind=mp) :: val

        select case (trim(fieldname))
        case ("efarr(1)")
            val = f%efarr(1)
        case ("efarr(2)")
            val = f%efarr(2)
        case ("bisection_shift")
            val = f%bisection_shift
        case default
            call f_err_throw("wrong argument for "//trim(fieldname))
        end select
    end function fermilevel_get_real

    function fermilevel_get_logical(f, fieldname) result(val)
        ! Calling arguments
        type(fermi_aux),intent(in) :: f      !< type that holds the internal data
        character(len=*),intent(in) :: fieldname
        logical :: val

        select case (trim(fieldname))
        case ("bisection_bounds_ok(1)")
            val = f%bisection_bounds_ok(1)
        case ("bisection_bounds_ok(2)")
            val = f%bisection_bounds_ok(2)
        case default
            call f_err_throw("wrong argument for "//trim(fieldname))
        end select
    end function fermilevel_get_logical


    ! Finds the real root of the equation ax**3 + bx**2 + cx + d which is closest to target_solution
    subroutine get_roots_of_cubic_polynomial(a, b, c, d, target_solution, solution)
      implicit none
    
      ! Calling arguments
      real(kind=mp),intent(in) :: a, b, c, d
      real(kind=mp),intent(in) :: target_solution
      real(kind=mp),intent(out) :: solution
    
      ! Local variables
      complex(kind=mp) :: a_c, b_c, c_c, d_c, Q_c, S_c, ttp_c, ttm_c
      complex(kind=mp),dimension(3) :: sol_c
      double complex :: test
      real(kind=mp) :: ttmin, tt
      integer :: i
    
      a_c=cmplx(a,0.d0,kind=mp)
      b_c=cmplx(b,0.d0,kind=mp)
      c_c=cmplx(c,0.d0,kind=mp)
      d_c=cmplx(d,0.d0,kind=mp)
    
      Q_c = sqrt( (2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c)**2 - 4*(b_c**2-3*a_c*c_c)**3 )
      S_c = ( .5d0*(Q_c+2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c) )**(1.d0/3.d0)
      ttp_c = cmplx(1.d0,sqrt(3.d0),kind=mp)
      ttm_c = cmplx(1.d0,-sqrt(3.d0),kind=mp)
    
      sol_c(1) = -b_c/(3*a_c) &
           - S_c/(3*a_c) &
           - (b_c**2-3*a_c*c_c)/(3*a_c*S_c)
      sol_c(2) = -b_c/(3*a_c) + (S_c*ttp_c)/(6*a_c) + ttm_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
      sol_c(3) = -b_c/(3*a_c) + (S_c*ttm_c)/(6*a_c) + ttp_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
      !!if (iproc==0) then
      !!    write(*,*) 'sol 1', sol_c(1)
      !!    write(*,*) 'sol 2', sol_c(2)
      !!    write(*,*) 'sol 3', sol_c(3)
      !!end if
    
      ! Select the real solution that is closest to target_solution
      ttmin=1.d100
      do i=1,3
          if (abs(aimag(sol_c(i)))>1.d-14) cycle !complex solution
          tt=abs(real(sol_c(i),kind=mp)-target_solution)
          if (tt<ttmin) then
              ttmin=tt
              solution=real(sol_c(i),kind=mp)
          end if
      end do
    
    end subroutine get_roots_of_cubic_polynomial


    real(kind=mp) function determinant(iproc, n, mat)
        implicit none
    
        ! Calling arguments
        integer,intent(in) :: iproc, n
        real(kind=mp),dimension(n,n),intent(in) :: mat
    
        ! Local variables
        integer :: i, info
        integer,dimension(n) :: ipiv
        real(kind=mp),dimension(n,n) :: mat_tmp
        real(kind=mp) :: sgn
    
        call vcopy(n**2, mat(1,1), 1, mat_tmp(1,1), 1)
    
        call dgetrf(n, n, mat_tmp, n, ipiv, info)
        if (info/=0) then
            if (iproc==0) write(*,'(a,i0,a)') 'ERROR in dgetrf, info=',info,'. Set determinant to zero.'
            determinant=0
            return
        end if
    
        determinant=1.d0
        do i=1,n
            determinant=determinant*mat_tmp(i,i)
        end do
    
        sgn=1.d0
        do i=1,n
            if(ipiv(i)/=i) then
                sgn=-sgn
            end if
        end do
    
        determinant=sgn*determinant   
    
    end function determinant


end module fermi_level
