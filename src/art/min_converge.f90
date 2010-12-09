!!****f* art/min_converge
!! FUNCTION
!!   Minimizes the energy at constant volume. 
!!   This minimization is done with only a minimal knowledge of the 
!!   physics of the problem so that it is portable
!!   
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine min_converge ( success )

  use defs
  use bigdft_forces
  implicit none

  !Arguments
  logical, intent(out) :: success

  !Local variables
  integer :: i, ierror
  real(kind=8),dimension(3) :: boxl
  !_______________________

  if ( iproc == 0 ) then              ! Report
     open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
     write(FLOG,'(1X,A)') ' RELAXATION'
     close(FLOG) 
  end if
                                      ! Relaxation subroutine in bigdft_forces.f90                     
  boxl = box * scala                  ! We compute at constant volume.
  call mingeo( NATOMS, boxl, pos, evalf_number, total_energy, success )
                                      ! Report.
  if ( iproc == 0 ) then 
     !open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
     !    & action = 'write', position = 'append', iostat = ierror )
     !write(FLOG,'(1X,A34,(1p,e17.10,0p))')&
     !    & ' - Relaxed Energy (eV)          : ',total_energy
     !write(FLOG,'(1X,A34,I17)') ' - Relax evalf_number           : ', evalf_number  
     !close(FLOG) 
     if ( .not. success ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A)') "Minimization exited before the geometry optimization converged,"
      write(FLOG,'(1X,A)') "this minimum will be rejected."
      close(FLOG) 
     end if
     write(*,"('',' BART: Relaxed energy : ',(1p,e17.10,0p))") total_energy
  end if 

END SUBROUTINE min_converge
!!***


!!****f* art/check_min( )
!! FUNCTION
!! SOURCE
!!
subroutine check_min( )

  use defs
  use lanczos_defs
  implicit none

  !Local variables
  integer :: i, ierror, repetition
  logical :: new_projection
  real(kind=8) :: a1
  !_______________________
                                      ! Report 
   if ( iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Starting Lanczos for minimum'
      write(FLOG,'(1X,A38)') '    Iter     Energy (eV)    Eigenvalue  a1' 
      close(FLOG) 
      write(*,*) 'BART: Starting Lanczos for minimum'
   end if

   new_projection = .true.          ! We do not use any previously computed direction. 

   if ( .not. setup_initial ) then
                                    ! We call lanczos twice.
      repetition = 2
   else 
                                    ! if not, four times.
      repetition = 4
   end if

   do i = 1, repetition
      call lanczos( NVECTOR_LANCZOS, new_projection , a1 )
                                    ! Report
      if ( iproc == 0 ) then 
         open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
             & action = 'write', position = 'append', iostat = ierror )
         write(FLOG,'(I6,3X,(1p,e17.10,0p),F12.6,x,F6.4)') i, total_energy, eigenvalue, a1
         close( FLOG ) 
         write(*,*) 'BART: Iter ', i, ' : ', total_energy,  eigenvalue, a1  
      end if
                                    ! Now we start from the previous direction. 
      new_projection= .false.   
                                    ! let's see the projection
      if ( setup_initial ) call print_proj ( i, 'M', projection, eigenvalue )
   end do
                                    ! Report 
   if ( iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Done Lanczos'
      close(FLOG) 
      write(*,*) 'BART: Done Lanczos'
   end if
 
END SUBROUTINE check_min
!!***
