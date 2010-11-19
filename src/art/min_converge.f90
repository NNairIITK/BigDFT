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
  use lanczos_defs
  use bigdft_forces
  implicit none

  !Arguments
  logical, intent(out) :: success

  !Local variables
  integer :: i, ierror
  logical :: new_projection
  real(kind=8) :: a1
  !_______________________
  if ( iproc == 0 ) then              ! Report
     open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
     write(FLOG,'(1X,A)') ' RELAXATION'
     close(FLOG) 
  end if
                                      ! Relaxation subroutine in bigdft_forces.f90                     
  call mingeo( NATOMS, box * scala, pos, evalf_number, total_energy, success )
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
                                      ! Are we in a real minimum?.
  If_lanc: if ( LANCZOS_MIN .and. success ) then  
                                      ! Report 
     if ( iproc == 0 ) then 
        open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
            & action = 'write', position = 'append', iostat = ierror )
        write(FLOG,*) ' Starting Lanczos for minimum'
        write(FLOG,'(1X,A38)') '    Iter     Energy (eV)    Eigenvalue'
        close(FLOG) 
        write(*,*) 'BART: Starting Lanczos for minimum'
     end if
     new_projection = .true.          ! We do not use any previously computed direction. 
                                      ! We call lanczos twice.
     do i = 1, 2
        call lanczos( NVECTOR_LANCZOS, new_projection , a1 )
                                      ! Report
        if ( iproc == 0 ) then 
           open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
               & action = 'write', position = 'append', iostat = ierror )
           write(FLOG,'(I6,3X,(1p,e17.10,0p),F12.6)') i, total_energy, eigenvalue  
           close( FLOG ) 
           write(*,*) 'BART: Iter ', i, ' : ', total_energy,  eigenvalue  
        end if
                                      ! Now we start from the previous direction. 
        new_projection= .false.   
     end do
                                      ! Write. 
     if ( iproc == 0 ) then 
        open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
            & action = 'write', position = 'append', iostat = ierror )
        write(FLOG,*) ' Done Lanczos'
        close(FLOG) 
        write(*,*) 'BART: Done Lanczos'
     end if

  end if If_lanc
 
END SUBROUTINE min_converge
!!***
