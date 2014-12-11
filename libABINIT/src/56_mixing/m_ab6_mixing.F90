#if defined HAVE_CONFIG_H
#include "config.h"
#endif
  
  module m_ab7_mixing

    use m_profiling
    use defs_basis

    implicit none

    private

    integer, parameter, public :: AB7_MIXING_NONE        = 0
    integer, parameter, public :: AB7_MIXING_EIG         = 1
    integer, parameter, public :: AB7_MIXING_SIMPLE      = 2
    integer, parameter, public :: AB7_MIXING_ANDERSON    = 3
    integer, parameter, public :: AB7_MIXING_ANDERSON_2  = 4
    integer, parameter, public :: AB7_MIXING_CG_ENERGY   = 5
    integer, parameter, public :: AB7_MIXING_CG_ENERGY_2 = 6
    integer, parameter, public :: AB7_MIXING_PULAY       = 7

    integer, parameter, public :: AB7_MIXING_POTENTIAL  = 0
    integer, parameter, public :: AB7_MIXING_DENSITY    = 1

    integer, parameter, public :: AB7_MIXING_REAL_SPACE     = 1
    integer, parameter, public :: AB7_MIXING_FOURRIER_SPACE = 2

    type, public :: ab7_mixing_object
       integer :: iscf
       integer :: nfft, nspden, kind, space

       logical :: useprec
       integer :: mffmem
       character(len = fnlen) :: diskCache
       integer :: n_index, n_fftgr, n_pulayit, n_pawmix
       integer, dimension(:), pointer :: i_rhor, i_vtrial, i_vresid, i_vrespc
       real(dp), dimension(:,:,:), pointer :: f_fftgr, f_atm
       real(dp), dimension(:,:), pointer :: f_paw

       ! Private
       integer :: n_atom
       real(dp), pointer :: xred(:,:), dtn_pc(:,:)
    end type ab7_mixing_object

    public :: ab7_mixing_new
    public :: ab7_mixing_deallocate

    public :: ab7_mixing_use_disk_cache
    public :: ab7_mixing_use_moving_atoms
    public :: ab7_mixing_copy_current_step

    public :: ab7_mixing_eval_allocate
    public :: ab7_mixing_eval
    public :: ab7_mixing_eval_deallocate

  contains

    subroutine init_(mix)
      implicit none

      type(ab7_mixing_object), intent(out) :: mix

      ! Default values.
      mix%iscf      = AB7_MIXING_NONE
      mix%mffmem    = 1
      mix%n_index   = 0
      mix%n_fftgr   = 0
      mix%n_pulayit = 7
      mix%n_pawmix  = 0
      mix%n_atom    = 0
      mix%useprec   = .true.

      call nullify_(mix)
    end subroutine init_

    subroutine nullify_(mix)


      implicit none

      type(ab7_mixing_object), intent(inout) :: mix

      ! Nullify internal pointers.
      nullify(mix%i_rhor)
      nullify(mix%i_vtrial)
      nullify(mix%i_vresid)
      nullify(mix%i_vrespc)
      nullify(mix%f_fftgr)
      nullify(mix%f_atm)
      nullify(mix%f_paw)
      nullify(mix%dtn_pc)
      nullify(mix%xred)
    end subroutine nullify_

    subroutine ab7_mixing_new(mix, iscf, kind, space, nfft, nspden, &
         & npawmix, errid, errmess, npulayit, useprec)
      implicit none

      type(ab7_mixing_object), intent(out) :: mix
      integer, intent(in) :: iscf, kind, space, nfft, nspden, npawmix
      integer, intent(out) :: errid
      character(len = 500), intent(out) :: errmess
      integer, intent(in), optional :: npulayit
      logical, intent(in), optional :: useprec
      
      integer :: ii, i_stat
      character(len = *), parameter :: subname = "ab7_mixing_new"

      ! Set default values.
      call init_(mix)

      ! Argument checkings.
      if (kind /= AB7_MIXING_POTENTIAL .and. kind /= AB7_MIXING_DENSITY) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_set_arrays: ERROR -',ch10,&
              & '  Mixing must be done on density or potential only.'
         return
      end if
      if (space /= AB7_MIXING_REAL_SPACE .and. &
           & space /= AB7_MIXING_FOURRIER_SPACE) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_set_arrays: ERROR -',ch10,&
              & '  Mixing must be done in real or Fourrier space only.'
         return
      end if
      if (iscf /= AB7_MIXING_EIG .and. iscf /= AB7_MIXING_SIMPLE .and. &
           & iscf /= AB7_MIXING_ANDERSON .and. &
           & iscf /= AB7_MIXING_ANDERSON_2 .and. &
           & iscf /= AB7_MIXING_CG_ENERGY .and. &
           & iscf /= AB7_MIXING_PULAY .and. &
           & iscf /= AB7_MIXING_CG_ENERGY_2) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, "(A,I0,A)") "Unknown mixing scheme (", iscf, ")."
         return
      end if
      errid = AB7_NO_ERROR

      ! Mandatory arguments.
      mix%iscf     = iscf
      mix%kind     = kind
      mix%space    = space
      mix%nfft     = nfft
      mix%nspden   = nspden
      mix%n_pawmix = npawmix

      ! Optional arguments.
      if (present(useprec)) mix%useprec = useprec

      ! Set-up internal dimensions.
      !These arrays are needed only in the self-consistent case
      if (iscf == AB7_MIXING_EIG) then
         !    For iscf==1, five additional vectors are needed
         !    The index 1 is attributed to the old trial potential,
         !    The new residual potential, and the new
         !    preconditioned residual potential receive now a temporary index
         !    The indices number 4 and 5 are attributed to work vectors.
         mix%n_fftgr=5 ; mix%n_index=1
      else if(iscf == AB7_MIXING_SIMPLE) then
         !    For iscf==2, three additional vectors are needed.
         !    The index number 1 is attributed to the old trial vector
         !    The new residual potential, and the new preconditioned
         !    residual potential, receive now a temporary index.
         mix%n_fftgr=3 ; mix%n_index=1
         if (.not. mix%useprec) mix%n_fftgr = 2
      else if(iscf == AB7_MIXING_ANDERSON) then
         !    For iscf==3 , four additional vectors are needed.
         !    The index number 1 is attributed to the old trial vector
         !    The new residual potential, and the new and old preconditioned
         !    residual potential, receive now a temporary index.
         mix%n_fftgr=4 ; mix%n_index=2
         if (.not. mix%useprec) mix%n_fftgr = 3
      else if (iscf == AB7_MIXING_ANDERSON_2) then
         !    For iscf==4 , six additional vectors are needed.
         !    The indices number 1 and 2 are attributed to two old trial vectors
         !    The new residual potential, and the new and two old preconditioned
         !    residual potentials, receive now a temporary index.
         mix%n_fftgr=6 ; mix%n_index=3
         if (.not. mix%useprec) mix%n_fftgr = 5
      else if(iscf == AB7_MIXING_CG_ENERGY .or. iscf == AB7_MIXING_CG_ENERGY_2) then
         !    For iscf==5 or 6, ten additional vectors are needed
         !    The index number 1 is attributed to the old trial vector
         !    The index number 6 is attributed to the search vector
         !    Other indices are attributed now. Altogether ten vectors
         mix%n_fftgr=10 ; mix%n_index=3
      else if(iscf == AB7_MIXING_PULAY) then
         !    For iscf==7, lot of additional vectors are needed
         !    The index number 1 is attributed to the old trial vector
         !    The index number 2 is attributed to the old residual
         !    The indices number 2 and 3 are attributed to two old precond. residuals
         !    Other indices are attributed now.
         if (present(npulayit)) mix%n_pulayit = npulayit
         mix%n_fftgr=2+2*mix%n_pulayit ; mix%n_index=1+mix%n_pulayit
         if (.not. mix%useprec) mix%n_fftgr = 1+2*mix%n_pulayit
      end if ! iscf cases

      ! Allocate new arrays.
      allocate(mix%i_rhor(mix%n_index), stat = i_stat)
      call memocc(i_stat, mix%i_rhor, 'mix%i_rhor', subname)
      allocate(mix%i_vtrial(mix%n_index), stat = i_stat)
      call memocc(i_stat, mix%i_vtrial, 'mix%i_vtrial', subname)
      allocate(mix%i_vresid(mix%n_index), stat = i_stat)
      call memocc(i_stat, mix%i_vresid, 'mix%i_vresid', subname)
      allocate(mix%i_vrespc(mix%n_index), stat = i_stat)
      call memocc(i_stat, mix%i_vrespc, 'mix%i_vrespc', subname)

      ! Setup initial values.
      if (iscf == AB7_MIXING_EIG) then
         mix%i_vtrial(1)=1 ; mix%i_vresid(1)=2 ; mix%i_vrespc(1)=3
      else if(iscf == AB7_MIXING_SIMPLE) then
         mix%i_vtrial(1)=1 ; mix%i_vresid(1)=2 ; mix%i_vrespc(1)=3
         if (.not. mix%useprec) mix%i_vrespc(1)=2
      else if(iscf == AB7_MIXING_ANDERSON) then
         mix%i_vtrial(1)=1 ; mix%i_vresid(1)=2
         if (mix%useprec) then
            mix%i_vrespc(1)=3 ; mix%i_vrespc(2)=4
         else
            mix%i_vrespc(1)=2 ; mix%i_vrespc(2)=3
         end if
      else if (iscf == AB7_MIXING_ANDERSON_2) then
         mix%i_vtrial(1)=1 ; mix%i_vtrial(2)=2
         mix%i_vresid(1)=3
         if (mix%useprec) then
            mix%i_vrespc(1)=4 ; mix%i_vrespc(2)=5 ; mix%i_vrespc(3)=6
         else
            mix%i_vrespc(1)=3 ; mix%i_vrespc(2)=4 ; mix%i_vrespc(3)=5
         end if
      else if(iscf == AB7_MIXING_CG_ENERGY .or. &
           & iscf == AB7_MIXING_CG_ENERGY_2) then
         mix%n_fftgr=10 ; mix%n_index=3
         mix%i_vtrial(1)=1
         mix%i_vresid(1)=2 ; mix%i_vresid(2)=4 ; mix%i_vresid(3)=7
         mix%i_vrespc(1)=3 ; mix%i_vrespc(2)=5 ; mix%i_vrespc(3)=8
         mix%i_rhor(2)=9 ; mix%i_rhor(3)=10
      else if(iscf == AB7_MIXING_PULAY) then
         do ii=1,mix%n_pulayit
            mix%i_vtrial(ii)=2*ii-1 ; mix%i_vrespc(ii)=2*ii
         end do
         mix%i_vrespc(mix%n_pulayit+1)=2*mix%n_pulayit+1
         mix%i_vresid(1)=2*mix%n_pulayit+2
         if (.not. mix%useprec) mix%i_vresid(1)=2
      end if ! iscf cases
    end subroutine ab7_mixing_new

    subroutine ab7_mixing_use_disk_cache(mix, fnametmp_fft)


      implicit none


      type(ab7_mixing_object), intent(inout) :: mix
      character(len = *), intent(in) :: fnametmp_fft
      
      if (len(trim(fnametmp_fft)) > 0) then
         mix%mffmem = 0
         write(mix%diskCache, "(A)") fnametmp_fft
      else
         mix%mffmem = 1
      end if
    end subroutine ab7_mixing_use_disk_cache

    subroutine ab7_mixing_use_moving_atoms(mix, natom, xred, dtn_pc)


      type(ab7_mixing_object), intent(inout) :: mix
      integer, intent(in) :: natom
      real(dp), intent(in), target :: dtn_pc(3, natom)
      real(dp), intent(in), target :: xred(3, natom)

      mix%n_atom = natom
      mix%dtn_pc => dtn_pc
      mix%xred => xred
    end subroutine ab7_mixing_use_moving_atoms

    subroutine ab7_mixing_copy_current_step(mix, arr_resid, errid, errmess, &
         & arr_respc, arr_paw_resid, arr_paw_respc, arr_atm)


      type(ab7_mixing_object), intent(inout) :: mix
      real(dp), intent(in) :: arr_resid(mix%space * mix%nfft, mix%nspden)
      integer, intent(out) :: errid
      character(len = 500), intent(out) :: errmess
      real(dp), intent(in), optional :: arr_respc(mix%space * mix%nfft, mix%nspden)
      real(dp), intent(in), optional :: arr_paw_resid(mix%n_pawmix), &
           & arr_paw_respc(mix%n_pawmix)
      real(dp), intent(in), optional :: arr_atm(3, mix%n_atom)

      if (.not. associated(mix%f_fftgr)) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_set_arr_current_step: ERROR -',ch10,&
              & '  Working arrays not yet allocated.'
         return
      end if
      errid = AB7_NO_ERROR

      mix%f_fftgr(:,:,mix%i_vresid(1)) = arr_resid(:,:)
      if (present(arr_respc)) mix%f_fftgr(:,:,mix%i_vrespc(1)) = arr_respc(:,:)
      if (present(arr_paw_resid)) mix%f_paw(:, mix%i_vresid(1)) = arr_paw_resid(:)
      if (present(arr_paw_respc)) mix%f_paw(:, mix%i_vrespc(1)) = arr_paw_respc(:)
      if (present(arr_atm)) mix%f_atm(:,:, mix%i_vresid(1)) = arr_atm(:,:)
    end subroutine ab7_mixing_copy_current_step

    subroutine ab7_mixing_eval_allocate(mix, istep)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_18_timing
!End of the abilint section

      implicit none

      type(ab7_mixing_object), intent(inout) :: mix
      integer, intent(in), optional :: istep

      integer :: istep_, i_stat, usepaw
      real(dp) :: tsec(2)
      character(len = *), parameter :: subname = "ab7_mixing_eval_allocate"

      istep_ = 1
      if (present(istep)) istep_ = istep

      ! Allocate work array.
      if (.not. associated(mix%f_fftgr)) then
         allocate(mix%f_fftgr(mix%space * mix%nfft,mix%nspden,mix%n_fftgr), stat = i_stat)
         call memocc(i_stat, mix%f_fftgr, 'mix%f_fftgr', subname)
         mix%f_fftgr(:,:,:)=zero
         if (mix%mffmem == 0 .and. istep_ > 1) then
            call timab(83,1,tsec)
            open(unit=tmp_unit,file=mix%diskCache,form='unformatted',status='old')
            rewind(tmp_unit)
            read(tmp_unit) mix%f_fftgr
            if (mix%n_pawmix == 0) close(unit=tmp_unit)
            call timab(83,2,tsec)
         end if
      end if
      ! Allocate PAW work array.
      if (.not. associated(mix%f_paw)) then
         usepaw = 0
         if (mix%n_pawmix > 0) usepaw = 1
         allocate(mix%f_paw(max(1,mix%n_pawmix),max(1,mix%n_fftgr * usepaw)), &
              & stat = i_stat)
         call memocc(i_stat, mix%f_paw, 'mix%f_paw', subname)
         if (mix%n_pawmix > 0) then
            mix%f_paw(:,:)=zero
            if (mix%mffmem == 0 .and. istep_ > 1) then
               read(tmp_unit) mix%f_paw
               close(unit=tmp_unit)
               call timab(83,2,tsec)
            end if
         end if
      end if
      ! Allocate atom work array.
      if (.not. associated(mix%f_atm)) then
         allocate(mix%f_atm(3,mix%n_atom,mix%n_fftgr), stat = i_stat)
         call memocc(i_stat, mix%f_atm, 'mix%f_atm', subname)
      end if
    end subroutine ab7_mixing_eval_allocate

    subroutine ab7_mixing_eval_deallocate(mix)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_18_timing
!End of the abilint section

      implicit none

      type(ab7_mixing_object), intent(inout) :: mix

      integer :: i_all, i_stat
      real(dp) :: tsec(2)
      character(len = *), parameter :: subname = "ab7_mixing_eval_deallocate"

      ! Save on disk and deallocate work array in case on disk cache only.
      if (mix%mffmem == 0) then
         call timab(83,1,tsec)
         open(unit=tmp_unit,file=mix%diskCache,form='unformatted',status='unknown')
         rewind(tmp_unit)
         ! VALGRIND complains not all of f_fftgr_disk is initialized
         write(tmp_unit) mix%f_fftgr
         if (mix%n_pawmix > 0) then
            write(tmp_unit) mix%f_paw
         end if
         close(unit=tmp_unit)
         call timab(83,2,tsec)
         i_all = -product(shape(mix%f_fftgr))*kind(mix%f_fftgr)
         deallocate(mix%f_fftgr, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%f_atm', subname)
         nullify(mix%f_fftgr)
         if (associated(mix%f_paw)) then
            i_all = -product(shape(mix%f_paw))*kind(mix%f_paw)
            deallocate(mix%f_paw, stat = i_stat)
            call memocc(i_stat, i_all, 'mix%f_paw', subname)
            nullify(mix%f_paw)
         end if
      end if
    end subroutine ab7_mixing_eval_deallocate

    subroutine ab7_mixing_eval(mix, arr, istep, nfftot, ucvol, &
         & mpi_comm, mpi_summarize, errid, errmess, &
         & reset, isecur, pawarr, pawopt, response, etotal, potden, &
         & resnrm, fnrm, fdot, user_data)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_56_mixing
!End of the abilint section

      implicit none

      type(ab7_mixing_object), intent(inout) :: mix
      integer, intent(in) :: istep, nfftot, mpi_comm
      logical, intent(in) :: mpi_summarize
      real(dp), intent(in) :: ucvol
      real(dp), intent(inout) :: arr(mix%space * mix%nfft,mix%nspden)
      integer, intent(out) :: errid
      character(len = 500), intent(out) :: errmess

      logical, intent(in), optional :: reset
      integer, intent(in), optional :: isecur, pawopt, response
      real(dp), intent(inout), optional, target :: pawarr(mix%n_pawmix)
      real(dp), intent(in), optional :: etotal
      real(dp), intent(in), optional :: potden(mix%space * mix%nfft,mix%nspden)
      real(dp), intent(out), optional :: resnrm
      optional :: fnrm, fdot
      integer, intent(in), optional :: user_data(:)

      interface
         function fdot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
           integer, intent(in) :: cplex,nfft,nspden,opt_denpot
           double precision, intent(in) :: x(*), y(*)
           integer, intent(in) :: user_data(:)

           double precision :: fdot
         end function fdot

         function fnrm(x,cplex,nfft,nspden,opt_denpot,user_data)
           integer, intent(in) :: cplex,nfft,nspden,opt_denpot
           double precision, intent(in) :: x(*)
           integer, intent(in) :: user_data(:)

           double precision :: fnrm
         end function fnrm
      end interface

      character(len = *), parameter :: subname = "ab7_mixing_eval"
      integer :: moveAtm, dbl_nnsclo, initialized, isecur_
      integer :: usepaw, pawoptmix_, response_, i_stat, i_all
      integer :: user_data_(2)
      real(dp) :: resnrm_
      real(dp), pointer :: pawarr_(:)

      ! Argument checkings.
      if (mix%iscf == AB7_MIXING_NONE) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_eval: ERROR -',ch10,&
              & '  No method has been chosen.'
         return
      end if
      if (mix%n_pawmix > 0 .and. .not. present(pawarr)) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_eval: ERROR -',ch10,&
              & '  PAW is used, but no pawarr argument provided.'
         return
      end if
      if (mix%n_atom > 0 .and. (.not. associated(mix%dtn_pc) .or. &
           & .not. associated(mix%xred))) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_eval: ERROR -',ch10,&
              & '  Moving atoms is used, but no xred or dtn_pc attributes provided.'
         return
      end if
      if ((present(fnrm) .or. present(fdot) .or. present(user_data)) .and. &
           & .not. (present(fnrm) .and. present(fdot) .and. present(user_data))) then
         errid = AB7_ERROR_MIXING_ARG
         write(errmess, '(a,a,a,a)' )ch10,&
              & ' ab7_mixing_eval: ERROR -',ch10,&
              & '  Passing optional norm and dot product routines without user_data argument.'
         return
      end if
      errid = AB7_NO_ERROR

      ! Miscellaneous
      moveAtm = 0
      if (mix%n_atom > 0) moveAtm = 1
      initialized = 1
      if (present(reset)) then
         if (reset) initialized = 0
      end if
      isecur_ = 0
      if (present(isecur)) isecur_ = isecur
      usepaw = 0
      if (mix%n_pawmix > 0) usepaw = 1
      pawoptmix_ = 0
      if (present(pawopt)) pawoptmix_ = pawopt
      response_ = 0
      if (present(response)) response_ = response
      if (present(pawarr)) then
         pawarr_ => pawarr
      else
         allocate(pawarr_(1), stat = i_stat)
         call memocc(i_stat, pawarr_, 'pawarr_', subname)
      end if

      ! Norm and dot products.
      if (.not. present(user_data)) then
         user_data_(1) = 0
         if (mpi_summarize) user_data_(1) = 1
         user_data_(2) = mpi_comm
      end if

      ! Do the mixing.
      resnrm_ = 0.d0
      if (mix%iscf == AB7_MIXING_EIG) then
         !  This routine compute the eigenvalues of the SCF operator
         call scfeig(istep, mix%space * mix%nfft, mix%nspden, &
              & mix%f_fftgr(:,:,mix%i_vrespc(1)), arr, &
              & mix%f_fftgr(:,:,1), mix%f_fftgr(:,:,4:5), errid, errmess)
      else if (mix%iscf == AB7_MIXING_SIMPLE .or. &
           & mix%iscf == AB7_MIXING_ANDERSON .or. &
           & mix%iscf == AB7_MIXING_ANDERSON_2 .or. &
           & mix%iscf == AB7_MIXING_PULAY) then
         if (present(user_data)) then
            call scfopt(mix%space, mix%f_fftgr,mix%f_paw,mix%iscf,istep,&
                 & mix%i_vrespc,mix%i_vtrial, mix%nfft,mix%n_pawmix,mix%nspden, &
                 & mix%n_fftgr,mix%n_index,mix%kind,pawoptmix_,usepaw,pawarr_, &
                 & resnrm_, arr, fnrm, fdot, user_data, errid, errmess)
         else
            call scfopt(mix%space, mix%f_fftgr,mix%f_paw,mix%iscf,istep,&
                 & mix%i_vrespc,mix%i_vtrial, mix%nfft,mix%n_pawmix,mix%nspden, &
                 & mix%n_fftgr,mix%n_index,mix%kind,pawoptmix_,usepaw,pawarr_, &
                 & resnrm_, arr, fnrm_default, fdot_default, user_data_, errid, errmess)
         end if
         !  Change atomic positions
         if((istep==1 .or. mix%iscf==AB7_MIXING_SIMPLE) .and. mix%n_atom > 0)then
            !    GAF: 2009-06-03
            !    Apparently there are not reason
            !    to restrict iscf=2 for ionmov=5
            mix%xred(:,:) = mix%xred(:,:) + mix%dtn_pc(:,:)
         end if
      else if (mix%iscf == AB7_MIXING_CG_ENERGY .or. &
           & mix%iscf == AB7_MIXING_CG_ENERGY_2) then
         !  Optimize next vtrial using an algorithm based
         !  on the conjugate gradient minimization of etotal
         if (.not. present(etotal) .or. .not. present(potden)) then
            errid = AB7_ERROR_MIXING_ARG
            write(errmess, '(a,a,a,a)' )ch10,&
                 & ' ab7_mixing_eval: ERROR -',ch10,&
                 & '  Arguments etotal or potden are missing for CG on energy methods.'
            return
         end if
         if (mix%n_atom == 0) then
            allocate(mix%xred(3,0), stat = i_stat)
            call memocc(i_stat, mix%xred, 'mix%xred', subname)
            allocate(mix%dtn_pc(3,0), stat = i_stat)
            call memocc(i_stat, mix%dtn_pc, 'mix%dtn_pc', subname)
         end if
         if (present(user_data)) then
            call scfcge(mix%space,dbl_nnsclo,mix%dtn_pc,etotal,mix%f_atm,&
                 & mix%f_fftgr,initialized,mix%iscf,isecur_,istep,&
                 & mix%i_rhor,mix%i_vresid,mix%i_vrespc,moveAtm,&
                 & mix%n_atom,mix%nfft,nfftot,&
                 & mix%nspden,mix%n_fftgr,mix%n_index,mix%kind,&
                 & response_,potden,ucvol,arr,mix%xred, &
                 & fnrm, fdot, user_data, errid, errmess)
         else
            call scfcge(mix%space,dbl_nnsclo,mix%dtn_pc,etotal,mix%f_atm,&
                 & mix%f_fftgr,initialized,mix%iscf,isecur_,istep,&
                 & mix%i_rhor,mix%i_vresid,mix%i_vrespc,moveAtm,&
                 & mix%n_atom,mix%nfft,nfftot,&
                 & mix%nspden,mix%n_fftgr,mix%n_index,mix%kind,&
                 & response_,potden,ucvol,arr,mix%xred, fnrm_default, &
                 & fdotn_default, user_data_, errid, errmess)
         end if
         if (mix%n_atom == 0) then
            i_all = -product(shape(mix%xred))*kind(mix%xred)
            deallocate(mix%xred, stat = i_stat)
            call memocc(i_stat, i_all, 'mix%xred', subname)
            i_all = -product(shape(mix%dtn_pc))*kind(mix%dtn_pc)
            deallocate(mix%dtn_pc, stat = i_stat)
            call memocc(i_stat, i_all, 'mix%dtn_pc', subname)
         end if
         if (dbl_nnsclo == 1) errid = AB7_ERROR_MIXING_INC_NNSLOOP
      end if
      
      if (present(resnrm)) resnrm = resnrm_
      if (.not. present(pawarr)) then
         i_all = -product(shape(pawarr_))*kind(pawarr_)
         deallocate(pawarr_, stat = i_stat)
         call memocc(i_stat, i_all, 'pawarr_', subname)
      end if
    end subroutine ab7_mixing_eval

    subroutine ab7_mixing_deallocate(mix)
      implicit none

      type(ab7_mixing_object), intent(inout) :: mix

      integer :: i_all, i_stat
      character(len = *), parameter :: subname = "ab7_mixing_deallocate"

      if (associated(mix%i_rhor)) then
         i_all = -product(shape(mix%i_rhor))*kind(mix%i_rhor)
         deallocate(mix%i_rhor, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%i_rhor', subname)
      end if
      if (associated(mix%i_vtrial)) then
         i_all = -product(shape(mix%i_vtrial))*kind(mix%i_vtrial)
         deallocate(mix%i_vtrial, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%i_vtrial', subname)
      end if
      if (associated(mix%i_vresid)) then
         i_all = -product(shape(mix%i_vresid))*kind(mix%i_vresid)
         deallocate(mix%i_vresid, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%i_vresid', subname)
      end if
      if (associated(mix%i_vrespc)) then
         i_all = -product(shape(mix%i_vrespc))*kind(mix%i_vrespc)
         deallocate(mix%i_vrespc, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%i_vrespc', subname)
      end if
      if (associated(mix%f_fftgr)) then
         i_all = -product(shape(mix%f_fftgr))*kind(mix%f_fftgr)
         deallocate(mix%f_fftgr, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%f_fftgr', subname)
      end if
      if (associated(mix%f_paw)) then
         i_all = -product(shape(mix%f_paw))*kind(mix%f_paw)
         deallocate(mix%f_paw, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%f_paw', subname)
      end if
      if (associated(mix%f_atm)) then
         i_all = -product(shape(mix%f_atm))*kind(mix%f_atm)
         deallocate(mix%f_atm, stat = i_stat)
         call memocc(i_stat, i_all, 'mix%f_atm', subname)
      end if

      call nullify_(mix)
    end subroutine ab7_mixing_deallocate

    function fnrm_default(x,cplex,nfft,nspden,opt_denpot,user_data)
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*)
      integer, intent(in) :: user_data(:)

      double precision :: fnrm_default
      real(dp) :: resid_new(1)

      call sqnormm_v(cplex,1,user_data(2),(user_data(1) /= 0),1,&
           & nfft,resid_new,1,nspden,opt_denpot,x)
      fnrm_default = resid_new(1)
    end function fnrm_default

    function fdot_default(x,y,cplex,nfft,nspden,opt_denpot,user_data)
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*), y(*)
      integer, intent(in) :: user_data(:)

      double precision :: fdot_default
      real(dp) :: prod_resid(1)

      call dotprodm_v(cplex,1,prod_resid,1,1,user_data(2),(user_data(1) /= 0),1,1,&
           & nfft,1,1,nspden,opt_denpot,x,y)
      fdot_default = prod_resid(1)
    end function fdot_default

    function fdotn_default(x,y,cplex,nfft,nspden,opt_denpot,user_data)
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*), y(*)
      integer, intent(in) :: user_data(:)

      double precision :: fdotn_default
      real(dp) :: prod_resid(1,1,1)

      call dotprodm_vn(cplex,1,x,prod_resid,1,1,user_data(2),(user_data(1) /= 0),1,1,&
           & 1,nfft,1,nspden,y)
      fdotn_default = prod_resid(1,1,1)
    end function fdotn_default
  end module m_ab7_mixing
