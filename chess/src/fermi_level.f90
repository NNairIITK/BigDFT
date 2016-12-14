!> @file
!!   Routines used to find the Fermi level during FOE
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


!> Determination of the Fermi level for the density matrix
module fermi_level
  use dictionaries, only: f_err_throw
  use sparsematrix_base
  use wrapper_linalg
  implicit none

  private

  ! Public parameters
  !> Function to determine the occupation numbers
  integer, parameter, public :: SMEARING_DIST_ERF   = 1  !< Tends to 0 and 1 faster \f$1/2\left[1-erf\left(\frac{E-\mu}{\delta E}\right)\right]\f$
  integer, parameter, public :: SMEARING_DIST_FERMI = 2  !< Normal Fermi distribution i.e.\f$\frac{1}{1+e^{E-\mu}/k_BT}\f$
  integer, parameter, public :: SMEARING_DIST_COLD1 = 3  !< Marzari's cold smearing with a=-.5634 (bumb minimization)
  integer, parameter, public :: SMEARING_DIST_COLD2 = 4  !< Marzari's cold smearing with a=-.8165 (monotonic tail)
  integer, parameter, public :: SMEARING_DIST_METPX = 5  !< Methfessel and Paxton (same as COLD with a=0)

  ! Public routines
  public :: init_fermi_level
  public :: determine_fermi_level
  public :: fermilevel_get_real
  public :: fermilevel_get_logical
  public :: eval_to_occ
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
      use dynamic_memory
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
      use dynamic_memory
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



    !> Finds the fermi level ef for an error function distribution with a width wf
    !! eval are the Kohn Sham eigenvalues and melec is the total number of electrons
    subroutine eval_to_occ(iproc, nproc, norbu, norbd, norb, nkpts, kwgts, &
               eval, occup, filewrite, not_initialized, wf0, occopt, efermi, eTS, &
               norbu_res, norbd_res)
      !use module_base
      use futile
      use yaml_output
      use numerics
      !use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level
      !use public_enums
      !use abi_interfaces_numeric, only: abi_derf_ab
      implicit none
      integer,intent(in) :: nkpts, norbu, norbd, norb
      real(mp),dimension(nkpts),intent(in) :: kwgts
      real(mp),dimension(nkpts*norb),intent(in) :: eval
      real(mp),dimension(nkpts*norb),intent(inout) :: occup
      logical, intent(in) :: filewrite, not_initialized
      integer, intent(in) :: iproc, nproc
      integer, intent(in) :: occopt
      real(mp), intent(in) :: wf0   ! width of Fermi function, i.e. k*T
      !type(orbitals_data), intent(inout) :: orbs
      real(mp),intent(inout) :: efermi, eTS
      integer, intent(in) :: norbu_res, norbd_res !<restricted values of norbu and norbd where the fermi level has to be found
      !local variables
      logical :: exitfermi
      !   real(gp), parameter :: pi=3.1415926535897932d0
      real(mp), parameter :: sqrtpi=sqrt(pi)
      real(mp), dimension(1,1,1) :: fakepsi
      integer :: ikpt,iorb,ii,newnorbu,newnorbd !,info_fermi
      real(mp) :: charge, chargef,wf,deltac
      real(mp) :: ef,electrons,dlectrons,factor,arg,argu,argd,corr,cutoffu,cutoffd,diff,full,res,resu,resd
      real(mp) :: a, x, xu, xd, f, df, tt
      !integer :: ierr
      !type(fermi_aux) :: ft


      exitfermi=.false.
      !if (iproc.lt.1)  write(1000+iproc,*)  'ENTER Fermilevel',norbu,norbd,occopt

      eTS=0.0_mp

      a = 0.d0
      select case (occopt)
      case  (SMEARING_DIST_ERF  )
      case  (SMEARING_DIST_FERMI)
      case  (SMEARING_DIST_COLD1) !Marzari's cold smearing  with a=-.5634 (bumb minimization)
         a=-.5634d0
      case  (SMEARING_DIST_COLD2) !Marzari's cold smearing  with a=-.8165 (monotonic tail)
         a=-.8165d0
      case  (SMEARING_DIST_METPX) !Methfessel and Paxton (same as COLD with a=0)
         a=0.d0
      case default
         call f_err_throw('Unrecognized smearing scheme',err_name='BIGDFT_RUNTIME_ERROR')
         !if(iproc==0) print *, 'unrecognized occopt=', occopt
         !stop
         return
      end select

      if (norbd==0) then
         full=2.d0   ! maximum occupation for closed shell  orbital
      else
         full=1.d0   ! maximum occupation for spin polarized orbital
      endif

      if (nkpts.ne.1 .and. filewrite) then
         call f_err_throw('Fermilevel: CANNOT write input.occ with more than one k-point',&
              err_name='BIGDFT_RUNTIME_ERROR')
         return
         !if (iproc == 0) print *,'Fermilevel: CANNOT write input.occ with more than one k-point'
         !stop
      end if
     
      !newnorbu=norbu
      newnorbu=min(norbu_res,norbu)
      !newnorbd=norbd
      newnorbd=min(norbd_res,norbd)


      charge=0.0_mp
      do ikpt=1,nkpts
         !number of zero orbitals for the given k-point
         !overall charge of the system
         do iorb=1,norb
            charge=charge+occup(iorb+(ikpt-1)*norb) * kwgts(ikpt)
         end do
      end do
      !melec=nint(charge)
      !if (iproc == 0) write(1000+iproc,*) 'charge,wf',charge,melec,wf0
      !call init_fermi_level(charge/full, 0.d0, ft, ef_interpol_det=1.d-12, verbosity=1)

      !!! Send all eigenvalues to all procs (presumably not necessary)
      !!call broadcast_kpt_objects(nproc, nkpts, norb, &
      !!     &   eval, orbs%ikptproc)

      if (wf0 > 0.0_mp) then
         ii=0
         !if (efermi == UNINITIALIZED(efermi)) then
         if (not_initialized) then
            !last value as a guess
            efermi = eval(norbu)
            ! Take initial value at gamma point.
            do iorb = 1, norbu
               if (occup(iorb) < 1.0_mp) then
                  efermi = eval(iorb)
                  exit
               end if
            end do
         end if
         ef=efermi

         ! electrons is N_electons = sum f_i * Wieght_i
         ! dlectrons is dN_electrons/dEf =dN_electrons/darg * darg/dEf= sum df_i/darg /(-wf) , darg/dEf=-1/wf
         !  f:= occupation # for band i ,  df:=df/darg
         wf=wf0
         loop_fermi: do ii=1,100
            !write(1000+iproc,*) 'iteration',ii,' -------------------------------- '
            factor=1.d0/(sqrt(pi)*wf)
            if (ii == 100 .and. iproc == 0) call yaml_warning('Fermilevel could not have been adjusted in the available iterations')
            electrons=0.d0
            dlectrons=0.d0
            do ikpt=1,nkpts
               do iorb=1,norbd+norbu
                  arg=(eval((ikpt-1)*norb+iorb)-ef)/wf
                  if (occopt == SMEARING_DIST_ERF) then
                     call abi_derf_ab(res,arg)
                     f =.5d0*(1.d0-res)
                     df=-safe_exp(-arg**2)/sqrtpi
                  else if (occopt == SMEARING_DIST_FERMI) then
                     f =1.d0/(1.d0+safe_exp(arg))
                     df=-1.d0/(2.d0+safe_exp(arg)+safe_exp(-arg))
                  else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &
                       &  occopt == SMEARING_DIST_METPX ) then
                     x= -arg
                     call abi_derf_ab(res,x)
                     f =.5d0*(1.d0+res +safe_exp(-x**2)*(-a*x**2 + .5d0*a+x)/sqrtpi)
                     df=-safe_exp(-x**2) * (a*x**3 -x**2 -1.5d0*a*x +1.5d0) /sqrtpi   ! df:=df/darg=-df/dx
                  else
                     f  = 0.d0
                     df = 0.d0
                  end if 
                  if (iorb > norbu+newnorbd .or. (iorb <= norbu .and. iorb > newnorbu)) then
                     f  = 0.d0
                     df = 0.d0
                  end if
                  !call yaml_map('arg,f,kwgts(ikpt)',(/arg,f,kwgts(ikpt)/))
                  electrons=electrons+ f  * kwgts(ikpt)  ! electrons := N_e(Ef+corr.)
                  dlectrons=dlectrons+ df * kwgts(ikpt)  ! delectrons:= dN_e/darg ( Well! later we need dN_e/dEf=-1/wf*dN_e/darg
                  !if(iproc==0) write(1000,*) iorb,arg,   f , df,dlectrons
               enddo
            enddo
            !call yaml_map('ef',ef)
            !call yaml_map('electrons',electrons)

            dlectrons=dlectrons/(-wf)  ! df/dEf=df/darg * -1/wf
            diff=-charge/full+electrons
            !if (iproc.lt.1) write(1000+iproc,*) diff,full,melec,real(melec,gp)
            !         if (iproc.lt.1) flush(1000+iproc)
            !if (iproc.lt.1) write(1000+iproc,*) diff,1.d-11*sqrt(electrons),wf
            !if (iproc.lt.1) flush(1000+iproc)
            !Exit criterion satiesfied, Nevertheles do one mor update of fermi level
            if (abs(diff) < 1.d-11*sqrt(electrons) .and. wf == wf0 ) exitfermi=.true.     ! Assume noise grows as sqrt(electrons)

            !alternative solution to avoid division by so high value
            !if (dlectrons == 0.d0) dlectrons=1.d-100  !line to be added
            if (dlectrons == 0.d0) then
               !always enter into first case below
               corr=0.d0
               if (diff > 0.d0) corr=1.d0*wf
               if (diff < 0.d0) corr=-1.d0*wf
               if (ii <= 50 .and. wf < 0.1d0) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature
            else
               corr=diff/abs(dlectrons) ! for case of no-monotonic func. abs is needed
               if (abs(corr).gt.wf) then   !for such a large correction the linear approximation is not any more valid
                  if (corr > 0.d0) corr=1.d0*wf
                  if (corr < 0.d0*wf) corr=-1.d0*wf
                  if (ii <= 50 .and. wf < 0.1d0) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature
               else
                  wf=max(wf0,.5d0*wf)
               endif
            end if
            ef=ef-corr  ! Ef=Ef_guess+corr.
            !if (iproc.lt.1) write(1000+iproc,'(i5,5(1pe17.8))') ii,electrons,ef,dlectrons,abs(dlectrons),corr
            !         if (iproc.lt.1) flush(1000+iproc)
            !call determine_fermi_level(ft, electrons, ef,info_fermi)
            !if (info_fermi /= 0) then
            !   call f_err_throw('Difficulties in guessing the new Fermi energy, info='//trim(yaml_toa(info_fermi)),&
            !        err_name='BIGDFT_RUNTIME_ERROR')
            !end if
            !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr) !debug
            if (exitfermi) exit loop_fermi
         end do loop_fermi

         do ikpt=1,nkpts
            argu=(eval((ikpt-1)*norb+norbu)-ef)/wf0
            argd=(eval((ikpt-1)*norb+norbu+norbd)-ef)/wf0
            if (occopt == SMEARING_DIST_ERF) then
               !error function
               call abi_derf_ab(resu,argu)
               call abi_derf_ab(resd,argd)
               cutoffu=.5d0*(1.d0-resu)
               cutoffd=.5d0*(1.d0-resd)
            else if (occopt == SMEARING_DIST_FERMI) then
               !Fermi function
               cutoffu=1.d0/(1.d0+safe_exp(argu))
               cutoffd=1.d0/(1.d0+safe_exp(argd))
            else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &
                 &  occopt == SMEARING_DIST_METPX ) then
               !Marzari's relation with different a
               xu=-argu
               xd=-argd
               call abi_derf_ab(resu,xu)
               call abi_derf_ab(resd,xd)
               cutoffu=.5d0*(1.d0+resu +safe_exp(-xu**2)*(-a*xu**2 + .5d0*a+xu)/sqrtpi)
               cutoffd=.5d0*(1.d0+resd +safe_exp(-xd**2)*(-a*xd**2 + .5d0*a+xd)/sqrtpi)
            end if
         enddo

         if ((cutoffu > 1.d-12 .or. cutoffd > 1.d-12) .and. iproc == 0) then
            call yaml_warning('Occupation numbers do not fill all available levels' // &
                 ' lastu=' // trim(yaml_toa(cutoffu,fmt='(1pe8.1)')) // &
                 ' lastd=' // trim(yaml_toa(cutoffd,fmt='(1pe8.1)')))
         end if
         !if (iproc.lt.1) write(1000+iproc,'(1x,a,1pe21.14,2(1x,e8.1))') 'Fermi level, Fermi distribution cut off at:  ',ef,cutoffu,cutoffd
         !      if (iproc.lt.1) flush(1000+iproc)
         efermi=ef

         !update the occupation number
         do ikpt=1,nkpts
            do iorb=1,norbu + norbd
               arg=(eval((ikpt-1)*norb+iorb)-ef)/wf0
               if (occopt == SMEARING_DIST_ERF) then
                  call abi_derf_ab(res,arg)
                  f=.5d0*(1.d0-res)
               else if (occopt == SMEARING_DIST_FERMI) then
                  f=1.d0/(1.d0+exp(arg))
               else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &
                    &  occopt == SMEARING_DIST_METPX ) then
                  x=-arg
                  call abi_derf_ab(res,x)
                  f =.5d0*(1.d0+res +exp(-x**2)*(-a*x**2 + .5d0*a+x)/sqrtpi)
               end if
               occup((ikpt-1)*norb+iorb)=full* f
               !if(iproc==0) print*,  eval((ikpt-1)*norb+iorb), occup((ikpt-1)*norb+iorb)
            end do
         end do

         !update electronic entropy S; eTS=T_ele*S is the electronic entropy term the negative of which is added to energy: Free energy = energy-T*S
         eTS=0.0_mp
         do ikpt=1,nkpts
            do iorb=1,norbu + norbd
               if (occopt == SMEARING_DIST_ERF) then
                  !error function
                  eTS=eTS+full*wf0/(2._mp*sqrt(pi))*&
                       safe_exp(-((eval((ikpt-1)*norb+iorb)-ef)/wf0)**2)
               else if (occopt == SMEARING_DIST_FERMI) then
                  !Fermi function
                  tt=occup((ikpt-1)*norb+iorb)
                  eTS=eTS-full*wf0*(tt*log(tt) + (1._mp-tt)*log(1._mp-tt))
               else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &
                    &  occopt == SMEARING_DIST_METPX ) then
                  !cold
                  eTS=eTS+0._mp  ! to be completed if needed
               end if
            end do
         end do
         !!! Sanity check on sum of occup.
         !!chargef=0.0_gp
         !!do ikpt=1,nkpts
         !!   do iorb=1,norb
         !!      chargef=chargef+kwgts(ikpt) * occup(iorb+(ikpt-1)*norb)
         !!   end do
         !!end do
         !!deltac=abs(charge - chargef)
         !!if (deltac > 1.e-9_gp .and. deltac < 1.e-6_gp) then
         !!   if (orbs%nspinor /= 4) call eigensystem_info(iproc,nproc,1.e-8_gp,0,orbs,fakepsi)
         !!   if (iproc==0) call yaml_warning('Failed to determine correctly the occupation number, expected='//yaml_toa(charge)// &
         !!        ', found='//yaml_toa(chargef))
         !!else if (deltac >= 1.e-6_gp) then
         !!   !if (abs(real(melec,gp)- chargef) > 1e-6)  then
         !!   if (orbs%nspinor /= 4) call eigensystem_info(iproc,nproc,1.e-8_gp,0,orbs,fakepsi)
         !!   call f_err_throw('Failed to determine correctly the occupation number, expected='//yaml_toa(charge)// &
         !!        ', found='//yaml_toa(chargef),err_name='BIGDFT_RUNTIME_ERROR')
         !!end if
      else if(full==1.0_mp) then
         !call eFermi_nosmearing(iproc,orbs)
         call eFermi_nosmearing(iproc, nkpts, norbu, norbd, norb, eval, occup, efermi)
         ! no entropic term when electronc temprature is zero
      end if

      !write on file the results if needed
      if (filewrite) then
         open(unit=11,file='input.occ',status='unknown')
         write(11,*)norbu,norbd
         do iorb=1,norb
            write(11,'(i5,e19.12,f10.6)')iorb,occup((ikpt-1)*norb+iorb) &
                 &   ,eval ((ikpt-1)*norb+iorb)
         end do
         close(unit=11)
      end if

    END SUBROUTINE eval_to_occ


    subroutine eFermi_nosmearing(iproc, nkpts, norbu, norbd, norb, eval, occup, efermi)
       use yaml_output
       implicit none
       integer,intent(in) :: iproc, nkpts, norbu, norbd, norb
       real(mp),dimension(nkpts*norb),intent(in) :: eval
       real(mp),dimension(nkpts*norb),intent(inout) :: occup
       !type(orbitals_data), intent(inout) :: orbs
       real(mp),intent(out) :: efermi
       !local variables
       integer :: iu,id,n,nzeroorbs,ikpt,iorb
       real(mp) :: charge
       real(mp) :: eF
    
       !SM: I think iu and id should be initialized to these values, in case the
       ! large if will not be executed.
       iu=norbu
       id=norbd
       eF = 0._mp
       do ikpt=1,nkpts
          !number of zero orbitals for the given k-point
          nzeroorbs=0
          !overall charge of the system
          charge=0.0_mp
          do iorb=1,norb
             if (occup(iorb+(ikpt-1)*norb) == 0.0_mp) then
                nzeroorbs=nzeroorbs+1
             else
                charge=charge+occup(iorb+(ikpt-1)*norb)
             end if
          end do
          if (nzeroorbs /= 0 .and. norbd .gt.0) then
             do iorb=1,norbu-1
                if (eval((ikpt-1)*norb+iorb) > eval((ikpt-1)*norb+iorb+1)) &
                   &   write(*,*) 'wrong ordering of up EVs',iorb,iorb+1
             end do
             do iorb=1,norbd-1
                if (eval((ikpt-1)*norb+iorb+norbu) > eval((ikpt-1)*norb+iorb+1+norbu))&
                   &   write(*,*) 'wrong ordering of dw EVs',iorb+norbu,iorb+1+norbu
             enddo
    
             iu=0
             id=0
             n=0
             do while (real(n,mp) < charge)
                if (eval((ikpt-1)*norb+iu+1) <= eval((ikpt-1)*norb+id+1+norbu)) then
                   iu=iu+1
                   eF=eval((ikpt-1)*norb+iu+1)
                else
                   id=id+1
                   eF=eval((ikpt-1)*norb+id+1+norbu)
                endif
                n=n+1
             enddo
             if (iproc==0) then
                !write(*,'(1x,a,1pe21.14,a,i4)') 'Suggested Homo energy level',eF,', Spin polarization',iu-id
                call yaml_map('Suggested Fermi Level',ef,fmt='(1pe21.14)')
                call yaml_map('Suggested Spin pol.',iu-id,fmt='(i4)')
             end if
             !write(*,*) 'up,down, up-down',iu,id,iu-id
          end if
       end do
       efermi=eF
       !assign the values for the occupation numbers
       do iorb=1,iu
          occup(iorb)=1.0_mp
       end do
       do iorb=iu+1,norbu
          occup(iorb)=0.0_mp
       end do
       do iorb=1,id
          occup(iorb+norbu)=1.0_mp
       end do
       do iorb=id+1,norbd
          occup(iorb+norbu)=0.0_mp
       end do
    
    END SUBROUTINE eFermi_nosmearing


    !!****f* BigDFT/abi_derf_ab
    !! FUNCTION
    !!   Error function in double precision, taken from libABINIT
    !!
    !! SOURCE
    !!
    subroutine abi_derf_ab(derf_yy,yy)
    
     !use abi_defs_basis
     implicit none
     real(mp),intent(in) :: yy
     real(mp),intent(out) :: derf_yy
     integer          ::  done,ii,isw
     real(mp), parameter :: &
           ! coefficients for 0.0 <= yy < .477
           &  pp(5)=(/ 113.8641541510502e0_mp, 377.4852376853020e0_mp,  &
           &           3209.377589138469e0_mp, .1857777061846032e0_mp,  &
           &           3.161123743870566e0_mp /)
      real(mp), parameter :: &
           &  qq(4)=(/ 244.0246379344442e0_mp, 1282.616526077372e0_mp,  &
           &           2844.236833439171e0_mp, 23.60129095234412e0_mp/)
      ! coefficients for .477 <= yy <= 4.0
      real(mp), parameter :: &
           &  p1(9)=(/ 8.883149794388376e0_mp, 66.11919063714163e0_mp,  &
           &           298.6351381974001e0_mp, 881.9522212417691e0_mp,  &
           &           1712.047612634071e0_mp, 2051.078377826071e0_mp,  &
           &           1230.339354797997e0_mp, 2.153115354744038e-8_mp, &
           &           .5641884969886701e0_mp /)
      real(mp), parameter :: &
           &  q1(8)=(/ 117.6939508913125e0_mp, 537.1811018620099e0_mp,  &
           &           1621.389574566690e0_mp, 3290.799235733460e0_mp,  &
           &           4362.619090143247e0_mp, 3439.367674143722e0_mp,  &
           &           1230.339354803749e0_mp, 15.74492611070983e0_mp/)
      ! coefficients for 4.0 < y,
      real(mp), parameter :: &
           &  p2(6)=(/ -3.603448999498044e-01_mp, -1.257817261112292e-01_mp,   &
           &           -1.608378514874228e-02_mp, -6.587491615298378e-04_mp,   &
           &           -1.631538713730210e-02_mp, -3.053266349612323e-01_mp/)
      real(mp), parameter :: &
           &  q2(5)=(/ 1.872952849923460e0_mp   , 5.279051029514284e-01_mp,    &
           &           6.051834131244132e-02_mp , 2.335204976268692e-03_mp,    &
           &           2.568520192289822e0_mp /)
      real(mp), parameter :: &
           &  sqrpi=.5641895835477563e0_mp, xbig=13.3e0_mp, xlarge=6.375e0_mp, xmin=1.0e-10_mp
      real(mp) ::  res,xden,xi,xnum,xsq,xx
    
     xx = yy
     isw = 1
    !Here change the sign of xx, and keep track of it thanks to isw
     if (xx<0.0e0_mp) then
      isw = -1
      xx = -xx
     end if
    
     done=0
    
    !Residual value, if yy < -6.375e0_mp
     res=-1.0e0_mp
    
    !abs(yy) < .477, evaluate approximation for erfc
     if (xx<0.477e0_mp) then
    ! xmin is a very small number
      if (xx<xmin) then
       res = xx*pp(3)/qq(3)
      else
       xsq = xx*xx
       xnum = pp(4)*xsq+pp(5)
       xden = xsq+qq(4)
       do ii = 1,3
        xnum = xnum*xsq+pp(ii)
        xden = xden*xsq+qq(ii)
       end do
       res = xx*xnum/xden
      end if
      if (isw==-1) res = -res
      done=1
     end if
    
    !.477 < abs(yy) < 4.0 , evaluate approximation for erfc
     if (xx<=4.0e0_mp .and. done==0 ) then
      xsq = xx*xx
      xnum = p1(8)*xx+p1(9)
      xden = xx+q1(8)
      do ii=1,7
       xnum = xnum*xx+p1(ii)
       xden = xden*xx+q1(ii)
      end do
      res = xnum/xden
      res = res* exp(-xsq)
      if (isw.eq.-1) then
         res = res-1.0e0_mp
      else
         res=1.0e0_mp-res
      end if
      done=1
     end if
    
    !y > 13.3e0_mp
     if (isw > 0 .and. xx > xbig .and. done==0 ) then
      res = 1.0e0_mp
      done=1
     end if
    
    !4.0 < yy < 13.3e0_mp  .or. -6.375e0_mp < yy < -4.0
    !evaluate minimax approximation for erfc
     if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
      xsq = xx*xx
      xi = 1.0e0_mp/xsq
      xnum= p2(5)*xi+p2(6)
      xden = xi+q2(5)
      do ii = 1,4
       xnum = xnum*xi+p2(ii)
       xden = xden*xi+q2(ii)
      end do
      res = (sqrpi+xi*xnum/xden)/xx
      res = res* exp(-xsq)
      if (isw.eq.-1) then
         res = res-1.0e0_mp
      else
         res=1.0e0_mp-res
      end if
     end if
    
    !All cases have been investigated
     derf_yy = res
    
    end subroutine abi_derf_ab
    !!***

end module fermi_level
