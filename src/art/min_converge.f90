!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> ART min_converge
!!   Minimizes the energy at constant volume. 
!!   This minimization is done with only a minimal knowledge of the 
!!   physics of the problem so that it is portable
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


!> ART check_min
subroutine check_min( stage )

  use defs
  use lanczos_defs
  implicit none

  !Arguments
  character(len=1), intent(in) :: stage

  !Local variables
  integer :: i, ierror, repetition
  logical :: new_projection
  real(kind=8) :: a1
  real(kind=8) :: min_energy ! First reference energy in lanczos
  !_______________________

  ! We check how it changes the energy of the system by applying the projection.
  IN_MINIMUN = .True.
                                      ! Report 

   if ( iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Starting Lanczos'
      close(FLOG) 
   end if

   new_projection = .true.          ! We do not use any previously computed direction. 

   if ( .not. setup_initial ) then
                                    ! We call lanczos twice.
      repetition = 2
   else 
                                    ! if not, four times.
      repetition = 4
   end if

   if ( iproc==0 ) write(*,*) "BART: INIT LANCZOS"  !debug
   do i = 1, repetition
      call lanczos( NVECTOR_LANCZOS, new_projection , a1 )
                                    ! Report
      if ( iproc == 0 ) then 
         open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
             & action = 'write', position = 'append', iostat = ierror )
         if ( i == 1 ) then         ! Our reference energy for the report.
             min_energy = lanc_energy
             write(FLOG,'(1X,A8,(1p,e17.10,0p),A12,1pe8.1,A3)') ' Em= ', min_energy, ' ( gnrm = ', my_gnrm, ' )'
             write(FLOG,'(A39)') '   Iter     Ep-Em (eV)   Eigenvalue  a1' 
         end if 
         write(FLOG,'(I6,3X,F10.3,4X,F12.6,X,F6.4)') i, min_energy-proj_energy, eigenvalue, a1
         close(FLOG) 
         write(*,*) 'BART: Iter ', i, ' : ', lanc_energy, proj_energy,  eigenvalue, a1  
      end if
                                    ! Now we start from the previous direction. 
      new_projection= .false.   
                                    ! let's see the projection
      if ( setup_initial ) call print_proj ( i, stage, projection, eigenvalue, DEL_LANCZOS )
   end do
                                    ! Report 
   if ( iproc == 0 ) then 
      write(*,*) "BART: END  LANCZOS"  !debug
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Done Lanczos'
      close(FLOG) 
   end if

  ! Default value in the activation part is false.
  IN_MINIMUN = .False.
 
END SUBROUTINE check_min
