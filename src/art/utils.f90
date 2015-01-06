!> @file
!! This file contains a series of utilities that could be used by a
!! number of programs. They suppose very little.
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> This subroutine prints the initial details for an event (ART)
subroutine print_event( ievent_current, temperat )

  use defs
  implicit none

  !Arguments
  integer, intent(in) :: ievent_current
  real(kind=8), intent(in) :: temperat

  !Local variables
  integer :: ierror

  write(*,*) 'BART: Simulation : ', ievent_current
  write(*,*) 'BART: Attempt    : ', atp
  write(*,*) 'BART: Starting from minconf : ', mincounter
  write(*,*) 'BART: Reference Energy (eV) : ', ref_energy
  write(*,*) 'BART: Temperature : ', temperat

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) ' __________________________________________________'
  write(FLOG,'(1X,A34,I17)') ' - Simulation                   : ', ievent_current
  write(FLOG,'(1X,A34,I17)') ' - Attempt                      : ', atp
  write(FLOG,'(1X,A34,I17)') ' - Starting from minconf        : ', mincounter
  write(FLOG,'(1X,A34,(1p,e17.10,0p))') ' - Reference Energy (eV)        : ', ref_energy 
  write(FLOG,'(1X,A34,F17.6)') ' - Temperature                  : ', temperat
  close(FLOG)

END SUBROUTINE print_event


!> It (ART) computes the distance between two configurations and 
!! the number of particles having moved by more than a THRESHOLD
subroutine displacement( posa, posb, delr, npart )

  use defs
  implicit none

  !Arguments
  real(kind=8), dimension(vecsize), intent(in), target :: posa
  real(kind=8), dimension(vecsize), intent(in), target :: posb 
  real(kind=8), intent(out)                            :: delr
  integer, intent(out)                                 :: npart

  !Local variables
  integer :: i
  real(kind=8),dimension(3) :: invbox
  real(kind=8), parameter :: THRESHOLD = 0.1d0  ! In Angstroems
  real(kind=8), dimension(:), pointer :: xa, ya, za, xb, yb, zb
  real(kind=8) :: delx, dely, delz, dr, dr2, delr2

  ! We first set-up pointers for the x, y, z components for posa and posb
  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  xb => posb(1:NATOMS)
  yb => posb(NATOMS+1:2*NATOMS)
  zb => posb(2*NATOMS+1:3*NATOMS)

  delr2 = 0.0d0
  npart = 0

  invbox = 1.0d0/box
  do i = 1, NATOMS
     if ( boundary == 'P') then

        delx = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
        dely = ( ya(i) - yb(i) ) - box(2) * nint(( ya(i) - yb(i) )*invbox(2))
        delz = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))
         
     else if ( boundary == 'S') then

        !be carefull with boundaries if surface: In bigdft the surface is 
        !always perpendicular to "y"

        delx = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
        dely = ( ya(i) - yb(i) )
        delz = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))
      else if ( boundary == 'F' ) then

        delx = ( xa(i) - xb(i) )
        dely = ( ya(i) - yb(i) )
        delz = ( za(i) - zb(i) )
      end if

     dr2   = delx*delx + dely*dely + delz*delz
     delr2 = delr2 + dr2
     dr    = sqrt(dr2) 

     ! could comment this part if you are not interested in counting the moved atoms 
     if ( dr > THRESHOLD ) then 
        npart = npart + 1
     end if
  end do

  delr = sqrt(delr2)

END SUBROUTINE displacement


!> ART store
!!    This subroutine stores the configurations at minima and activated points
!!    By definition, it uses pos, box and scala
subroutine store( fname )

  use defs
  implicit none

  !Arguments
  character(len=7 ), intent(in) :: fname

  !Local variables
  integer :: i, ierror
  real(kind=8),dimension(3) :: boxl
  character(len=*), parameter :: extension = ".xyz"
  character(len=20) :: fnamexyz
  character(len=4), dimension(natoms) :: frzchain

  boxl = box * scala                  ! Update the box size

  write(*,*) ' Writing to file : ', fname
   
  open(unit=FCONF,file=fname,status='unknown',action='write',iostat=ierror)
  write(FCONF,*) 'run_id: ', mincounter
  write(FCONF,*) 'total_energy: ', total_energy
  write(FCONF,*) boxl
  do i=1, NATOMS
     write(FCONF,'(1x,i6,3(2x,f16.8))') typat(i), x(i), y(i), z(i)
  end do
  close(FCONF)

  ! Added by Fedwa El-Mellouhi July 2002, writes the configuration in .xyz format. 
  ! Modified by E. Machado-charry for v_sim and BigDFT.
!  if ( write_xyz ) then

     ! If there is a constraint over a given atom, is written in the geometry file.
     do i = 1, NATOMS
        if      ( constr(i)== 0) then
           frzchain(i)='    '
        else if (constr(i) == 1) then
           frzchain(i)='   f'
        else if (constr(i) == 2) then
           frzchain(i)='  fy'
        else if (constr(i) == 3) then
           frzchain(i)=' fxz'
        end if
     end do

     fnamexyz = trim(fname) // extension
     write(*,*) ' Writing to file : ', fnamexyz
     open(unit=XYZ,file=fnamexyz,status='unknown',action='write',iostat=ierror)
     write(XYZ,*) NATOMS , 'angstroemd0' 
     if (boundary == 'P') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))')'periodic', (boxl(i),i=1,3)
     else if (boundary == 'S') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))')'surface',(boxl(i),i=1,3)
     else
        write(XYZ,*)'free'
     end if
     do i=1, NATOMS
        write(XYZ,'(1x,A2,3(2x,f16.8),2x,a4)')   Atom(i), x(i), y(i), z(i), frzchain(i)
     end do
     close(XYZ)
!  end if

END SUBROUTINE store


!> ART save_intermediate
!!    It saves the configuration at every step in xyz format.
!!    The name of the file will look like this example:
!!    p_1001_05_030_K.xyz
!!    1001:: is the 'mincounter' of the event, 
!!    05  :: the attempt
!!    030 :: the step
!!    K is the argument 'stage' ( K= basin activation, L=lanczos, D= DIIS )
subroutine save_intermediate( stage )

  use defs
  implicit none

  !Arguments
  character(len=1), intent(in) :: stage

  !Local variables
  integer :: i, ierror
  real(kind=8), dimension(3) :: boxl
  character(len=40) :: fname
  character(len=4)  :: scounter, rcounter, pcounter
  character(len=4), dimension(natoms) :: frzchain

  boxl = box * scala                  ! Update the box size

                                      ! subroutines in utils.f90 
  if ( iproc == 0 ) then
                                      ! set up of xyz file.
     call convert_to_chain( mincounter, 4, scounter ) 
     call convert_to_chain( atp       , 2, rcounter )
     call convert_to_chain( pas       , 3, pcounter )
     fname = 'p_'//trim(scounter)//'_'//trim(rcounter)//'_'//trim(pcounter)//'_'//stage//".xyz"
     fname = trim(fname)

     ! If there is a constraint over a given atom, is written in the geometry file.
     do i = 1, NATOMS
        if      ( constr(i)== 0) then
           frzchain(i)='    '
        else if (constr(i) == 1) then
           frzchain(i)='   f'
        else if (constr(i) == 2) then
           frzchain(i)='  fy'
        else if (constr(i) == 3) then
           frzchain(i)=' fxz'
        end if
     end do

     open(unit=XYZ,file=fname,status='unknown',action='write',iostat=ierror)

     write(XYZ,*) NATOMS,  'angstroemd0' 

     if (boundary == 'P') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))') 'periodic', (boxl(i),i=1,3)
     else if (boundary == 'S') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p))') 'surface',  (boxl(i),i=1,3)
     else
        write(XYZ,*)'free'
     end if

     do i= 1, NATOMS
        write(XYZ,'(1x,A2,3(2x,f16.8),2x,a4)')   Atom(i), x(i), y(i), z(i), frzchain(i)
     end do

     write(XYZ,*) '# simulation ',scounter,", attempt ", atp,", step ", pas,", stage ", stage
     write(XYZ,'(a,(1p,e17.10,0p))') ' # total energy (eV) : ', total_energy

     close(XYZ)

  end if

END SUBROUTINE save_intermediate

!subroutine convert_to_chain(init_number,chain)
!  integer, intent(in) :: init_number
!  character(len=7), intent(out) :: chain
!
!  ! write to string
!  write(chain,'(I7)') init_number
!  ! flush left
!  chain = adjustl(chain)
!
!end subroutine convert_to_chain


!> ART convert_to_chain
!!    It takes an integer and transforms it into a chain of character. 
subroutine convert_to_chain( init_number, length, chain )

  implicit none

  !Arguments
  integer,          intent(in)  :: init_number
  integer,          intent(in)  :: length      !can be up to 4 
  character(len=4), intent(out) :: chain

  !Local variables
  character(len=10) :: digit = '0123456789'
  integer :: i, decades, divider, remainder, number, lm

  number = init_number
  if ( number == 0 ) then
     chain =''
     do i = 1, length
        chain = trim(chain)//"0"
     end do 
     return
  else
     decades = int(log10( 1.0d0 * real(number,kind=8))) + 1
  end if

  if ( decades > length ) then   ! WARNING: chain will be is zero. 
     chain =''
     do i = 1, length
        chain = trim(chain)//"0"
     end do 
     return 
  end if

  divider = 1
  do i = 2, decades
     divider =  divider * 10 
  end do

  lm = length - decades - 1 
  chain = ''
  do i = 1, decades
     remainder = number / divider  + 1
     chain = trim(chain)// digit(remainder:remainder)
     remainder = remainder -1
     number = number - remainder * divider
     divider = divider / 10
  end do
    
  do i= 1, length
     if ( len(trim(chain)) == length ) exit
     chain = "0"//trim(chain) 
  end do 

END SUBROUTINE convert_to_chain


!> ART write_refconfig
! THIS IS NOT NECESSARY IN BIGDFT
!!!    This subroutine writes the atomic positions and others to a "refconfig" file
!!!    which will be used a the reference point until a new events gets accepted.
!subroutine write_refconfig( )
!
!  use defs
!  implicit none
!
!  !Local variables
!  integer :: i, ierror
!  real(kind=8),dimension(3) :: boxl
!
!  boxl = box * scala                  ! Update the box size 
!
!                                      ! switch replace for unknown
!  open(unit=FREFCONFIG,file=REFCONFIG,status='unknown',action='write',iostat=ierror) 
!  write(FREFCONFIG,*) 'run_id: ', refcounter
!  write(FREFCONFIG,*) 'total_energy: ', total_energy
!  write(FREFCONFIG,*) boxl
!  do i = 1, NATOMS
!     write(FREFCONFIG,'(1x,i6, 3(2x,F16.8))') typat(i), x(i), y(i), z(i)
!  end do
!  close(FREFCONFIG)
!
!END SUBROUTINE write_refconfig


!> ART print_proj
subroutine print_proj( repetitions, stage, vector, eigenvalue, stepsize )

  use defs

  implicit none
  integer,          intent(in) :: repetitions
  character(len=1), intent(in) :: stage
  real(kind=8),     intent(in) :: eigenvalue
  real(kind=8),     intent(in), dimension(3*natoms) :: vector 
  real(kind=8),     intent(in) :: stepsize 

  !Local variables

  integer :: ierror
  integer :: i
  real(kind=8), dimension(3) :: boxl
  real(kind=8), allocatable  :: pc(:,:)
  character(len=40) :: fname
  character(len=4)  :: rcounter

  allocate(pc(3, NATOMS))
  do i = 1, NATOMS, 1
     pc(:, i) = (/ vector(i), vector(natoms + i), vector(2 * natoms + i) /) 
  end do

  if ( iproc == 0 ) then

     boxl = box * scala  ! Update the box size

     call convert_to_chain( repetitions, 2, rcounter )
     fname = 'proj_'//trim(rcounter)//'_'//stage//".xyz"
     fname = trim(fname)

     open(unit=XYZ,file=fname,status='unknown',action='write',iostat=ierror)

     write(XYZ,*) NATOMS,  'angstroem' 

     if (boundary == 'P') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p),2x,a,2x,F12.6)')'periodic', (boxl(i),i=1,3),'#', eigenvalue
     else if (boundary == 'S') then
        write(XYZ,'(a,3(1x,1p,e24.17,0p),2x,a,2x,F12.6)')'surface',  (boxl(i),i=1,3),'#', eigenvalue
     else
        write(XYZ,*)'free ',' # ', eigenvalue
     end if

     do i= 1, NATOMS
        write(XYZ,'(1x,A2,3(2x,f16.8),2x,A,2x,I3,3(2x,f12.8))') &
       &                               Atom(i), x(i) + stepsize*pc(1,i), &
       &                                        y(i) + stepsize*pc(2,i), &
       &                                        z(i) + stepsize*pc(3,i), &
       &                               '#',i, pc(1,i), pc(2,i),pc(3,i)
     end do

     close(XYZ)
  end if

  deallocate(pc)

END SUBROUTINE print_proj 


!> ART boundary_cond
subroutine boundary_cond ( posr, posa, posb )

  use defs, only : natoms, vecsize, box, boundary
  implicit none

  !Arguments
  real(kind=8), dimension(vecsize), intent(out), target :: posr
  real(kind=8), dimension(vecsize), intent(in), target  :: posa
  real(kind=8), dimension(vecsize), intent(in), target  :: posb 

  !Local variables
  integer :: i
  real(kind=8),dimension(3) :: invbox
  real(kind=8), dimension(:), pointer :: xr, yr, zr, xa, ya, za, xb, yb, zb

  ! We first set-up pointers for the x, y, z components for posr, posa and posb

  xr => posr(1:NATOMS)
  yr => posr(NATOMS+1:2*NATOMS)
  zr => posr(2*NATOMS+1:3*NATOMS)

  xa => posa(1:NATOMS)
  ya => posa(NATOMS+1:2*NATOMS)
  za => posa(2*NATOMS+1:3*NATOMS)

  xb => posb(1:NATOMS)
  yb => posb(NATOMS+1:2*NATOMS)
  zb => posb(2*NATOMS+1:3*NATOMS)

  invbox = 1.0d0/box

  do i = 1, natoms
     if ( boundary == 'P') then

        xr(i) = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
        yr(i) = ( ya(i) - yb(i) ) - box(2) * nint(( ya(i) - yb(i) )*invbox(2))
        zr(i) = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))
         
     else if ( boundary == 'S') then

        !be carefull with boundaries if surface: In bigdft the surface is 
        !always perpendicular to "y"

        xr(i) = ( xa(i) - xb(i) ) - box(1) * nint(( xa(i) - xb(i) )*invbox(1))
        yr(i) = ( ya(i) - yb(i) )
        zr(i) = ( za(i) - zb(i) ) - box(3) * nint(( za(i) - zb(i) )*invbox(3))
      else if ( boundary == 'F' ) then

        xr(i) = ( xa(i) - xb(i) )
        yr(i) = ( ya(i) - yb(i) )
        zr(i) = ( za(i) - zb(i) )
     end if
  enddo

END SUBROUTINE boundary_cond
