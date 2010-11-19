!!****f* BigDFT/utils
!! FUNCTION 
!!    This file contains a series of utilities that could be used by a
!!    number of program. They suppose very little.
!!
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!!****f* BigDFT/print_newevent
!! FUNCTION
!!    This subroutine prints the initial details for a new events
!! SOURCE
!!
subroutine print_newevent( ievent_current, temperat )

  use defs
  implicit none

  !Arguments
  integer, intent(in) :: ievent_current
  real(kind=8), intent(in) :: temperat

  !Local variables
  integer :: ierror

  write(*,*) 'BART: Simulation : ', ievent_current
  write(*,*) 'BART: Starting from minconf : ', mincounter
  write(*,*) 'BART: Reference Energy (eV) : ', ref_energy
  write(*,*) 'BART: Temperature : ', temperat

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) ' _______________________________________'
  write(FLOG,'(1X,A34,I17)') ' - Simulation                   : ', ievent_current
  write(FLOG,'(1X,A34,I17)') ' - Starting from minconf        : ', mincounter
  write(FLOG,'(1X,A34,(1p,e17.10,0p))') ' - Reference Energy (eV)        : ', ref_energy 
  write(FLOG,'(1X,A34,F17.6)') ' - Temperature                  : ', temperat
  close(FLOG)

END SUBROUTINE print_newevent
!!***


!!****f* BigDFT/convert_to_chain
!! FUNCTION
!!    This subroutine takes an integer and transforms it into a
!!    chain of character.
!! SOURCE
!!
subroutine convert_to_chain( init_number, chain )

  implicit none

  !Arguments
  integer,          intent(in)  :: init_number
  character(len=4), intent(out) :: chain

  !Local variables
  character(len=10), parameter :: digit = '0123456789'
  integer :: i, decades, divider, remainder, number

  number = init_number
  if ( number == 0 ) then
     chain = '0'
     return
  else
     decades = log10( 1.0d0 * number) + 1
  end if

  divider = 1
  do i = 2, decades
     divider =  divider * 10 
  end do
     
  chain = 'AB'
  do i = 1, decades
     remainder = number / divider  + 1
     chain = chain(1:i-1) // digit(remainder:remainder)
     remainder = remainder -1
     number = number - remainder * divider
     divider = divider / 10
  end do

END SUBROUTINE convert_to_chain 
!!***

!!****f* BigDFT/displacement
!! FUNCTION
!!    It computes the distance between two configurations and 
!!    the number of particles having moved by more than a THRESHOLD
!! SOURCE
!!
subroutine displacement( posa, posb, delr, npart )

  use defs
  implicit none

  !Arguments
  real(kind=8), dimension(vecsize), intent(in), target :: posa
  real(kind=8), dimension(vecsize), intent(in), target :: posb 
  real(kind=8), intent(out)                            :: delr
  integer, intent(out)                            :: npart

  !Local variables
  integer :: i, j
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

  do i = 1, NATOMS
     delx = ( xa(i) - xb(i) )
     dely = ( ya(i) - yb(i) )
     delz = ( za(i) - zb(i) )

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
!!***


!!****f* BigDFT/store
!! FUNCTION
!!    Subroutine store
!!    This subroutine stores the configurations at minima and activated points
!!    By definition, it uses pos, box and scala
!! SOURCE
!! 
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

  ! added by Fedwa El-Mellouhi July 2002, writes the configuration in jmol format 
  if ( write_jmol ) then
     fnamexyz = trim(fname) // extension
     write(*,*) ' Writing to file : ', fnamexyz
     open(unit=XYZ,file=fnamexyz,status='unknown',action='write',iostat=ierror)
     write(XYZ,*) NATOMS , 'angstroem' 
     if (boundary == 'P') then
        write(XYZ,'(a,3(1x,1pe24.17))')'periodic', (boxl(i),i=1,3)
     else if (boundary == 'S') then
        write(XYZ,'(a,3(1x,1pe24.17))')'surface',(boxl(i),i=1,3)
     else
        write(XYZ,*)'free'
     end if
     do i=1, NATOMS
        write(XYZ,'(1x,A2,3(2x,f16.8))')   Atom(i), x(i), y(i), z(i)
     end do
     close(XYZ)
  end if

END SUBROUTINE store
!!***


!!****f* BigDFT/write_refconfig
!! FUNCTION
!!    This subroutine writes the atomic positions and others to a "refconfig" file
!!    which will be used a the reference point until a new events gets accepted.
!! SOURCE
!!
subroutine write_refconfig( )

  use defs
  implicit none

  !Local variables
  integer :: i, ierror
  real(kind=8),dimension(3) :: boxl

  boxl = box * scala                  ! Update the box size 

                                      ! switch replace for unknown
  open(unit=FREFCONFIG,file=REFCONFIG,status='unknown',action='write',iostat=ierror) 
  write(FREFCONFIG,*) 'run_id: ', refcounter
  write(FREFCONFIG,*) 'total_energy: ', total_energy
  write(FREFCONFIG,*) boxl
  do i = 1, NATOMS
     write(FREFCONFIG,'(1x,i6, 3(2x,F16.8))') typat(i), x(i), y(i), z(i)
  end do
  close(FREFCONFIG)

END SUBROUTINE write_refconfig
!!***


!!****f* BigDFT/store_part
!! FUNCTION
!!    This subroutine stores partial configurations. 
!! SOURCE
!! 
subroutine store_part( fname, scounter, rcounter, stage )

  use defs
  implicit none

  !Arguments
  character(len=20), intent(in) :: fname
  character(len=4),  intent(in) :: scounter
  character(len=4),  intent(in) :: rcounter
  character(len=1),  intent(in) :: stage

  !Local variables
  integer :: ierror
  integer :: i
  real(kind=8), dimension(3) :: boxl
  logical :: exists_already
  character(len=24) :: fnamexyz
  character(len=*), parameter :: extension = ".xyz"
  character(len=9) :: digit = "123456789"
  character(len=150) :: commande
  character(len=100) :: temp

  boxl = box * scala  ! Update the box size

  fnamexyz = trim(fname)// extension

  temp=  fnamexyz
  do i = 1, 9 
     inquire(file=temp,exist=exists_already)
     if ( exists_already )  then
        temp = trim(fnamexyz) // "." // digit(i:i)
     else 
        if ( i > 1 ) then
           commande = "mv " // fnamexyz// "  " // temp
           call system(commande)
        end if
        exit
     end if 
  end do
   
  open(unit=XYZ,file=fnamexyz,status='unknown',action='write',iostat=ierror)

  write(XYZ,*) NATOMS,  'angstroem' 

  if (boundary == 'P') then
     write(XYZ,'(a,3(1x,1pe24.17))')'periodic', (boxl(i),i=1,3)
  else if (boundary == 'S') then
     write(XYZ,'(a,3(1x,1pe24.17))')'surface',(boxl(i),i=1,3)
  else
     write(XYZ,*)'free'
  end if

  do i= 1, NATOMS
     write(XYZ,'(1x,A2,3(2x,f16.8))')   Atom(i), x(i), y(i), z(i)
  end do

  write(XYZ,*) '# simulation ',scounter," stage ", stage, " iteration ",rcounter
  write(XYZ,'(a,(1p,e17.10,0p))') ' # total energy (eV) : ', total_energy

  close(XYZ)

END SUBROUTINE store_part
!!***


!!****f* BigDFT/convert_to_chain2
!! FUNCTION
!!    The subroutine convert_to_chain takes an integer and transforms it into a
!!    chain of character. It is the same convert_to_chain, but with a small
!!    modification.
!! SOURCE
!!
subroutine convert_to_chain_2( init_number, chain )

  implicit none

  !Arguments
  integer,          intent(in)  :: init_number
  character(len=4), intent(out) :: chain

  !Local variables
  character(len=10) :: digit = '0123456789'
  integer :: i, decades, divider, remainder, number, lm

  number = init_number
  if ( number == 0 ) then
     chain = '000'
     return
  else
     decades = log10( 1.0d0 * number) + 1
  end if

  divider = 1
  do i = 2, decades
     divider =  divider * 10 
  end do
     
  lm = 3 - decades - 1 
  chain = '0000'
  do i = 1, decades
     remainder = number / divider  + 1
     chain = chain(1:lm+i) // digit(remainder:remainder)
     remainder = remainder -1
     number = number - remainder * divider
     divider = divider / 10
  end do

END SUBROUTINE convert_to_chain_2
!!***


!!****f* BigDFT/save_intermediate
!! FUNCTION
!!    It saves the configuration at every step in xyz format.
!!    The name of the file will be:
!!     conf_1001_030_K.xyz
!!    where 1001 is the 'mincounter' is 1001,
!!    K is the argument 'stage' ( K= basin activation, L=lanczos, D= DIIS),
!!    and 030 is the step
!! SOURCE
!! 
subroutine save_intermediate( stage )

  use defs
  implicit none

  !Arguments
  character(len=1), intent(in) :: stage

  !Local variables
  character(len=20) :: fname
  character(len=4)  :: rcounter, scounter
                                      ! subroutines in utils.f90 
  if ( iproc == 0 ) then
     call convert_to_chain( mincounter, scounter )
     call convert_to_chain_2( pas, rcounter )
     fname = 'conf_'//trim(scounter)//'_'//trim(rcounter)//'_'//stage
     call store_part( fname, scounter, rcounter, stage )
  end if

END SUBROUTINE save_intermediate
!!***


!!****f* BigDFT/move_intermediate
!! FUNCTION
!!
!! SOURCE
!! 
subroutine move_intermediate( )

  use defs
  implicit none

  !Local variables
  integer :: i, j
  logical :: exists_already
  character(len=4)   :: scounter
  character(len=150) :: commande
  character(len=25)  :: fname
  character(len=9)   :: digit = "123456789"

  if ( iproc == 0 ) then
     call convert_to_chain( mincounter, scounter )

     fname = 'conf_'//trim(scounter)//"_000_K.xyz"

     do i = 1, 9 

        inquire( file = fname, exist = exists_already )

        if ( exists_already ) then
           fname = 'conf_'//trim(scounter)//"_000_K.xyz"//"."// digit(i:i)
        else 
           if ( i > 1 ) then

              j=i-1    
              commande = "ls -1 " //" conf_"//trim(scounter)//"*.xyz"//& 
            & " |sed -e 's|\(.*\)|mv \1 \1."//digit(j:j)//"|g'|sh"
              call system ( commande )

           end if
           exit
        end if 
      end do
  end if

END SUBROUTINE move_intermediate
!!***
