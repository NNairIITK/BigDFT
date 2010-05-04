!!****f* art/utils
!! DESCRIPTION 
!!   This file contains a series of utilities that could be used by a
!!   number of program. They suppose very little.
!!
!!****f* utils/convert_to_chain
!!   This subroutine takes an integer and transforms it into a
!!   chain of character.
!!
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
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

!!****f* utils/center
!! FUNCTION
!!   It places the center of mass of a 3D vector at (0,0,0) 
!!
!! SOURCE
!!
subroutine center( vector, vecsize )

  implicit none

  !Arguments
  integer, intent(in) :: vecsize
  real(kind=8), dimension(vecsize), intent(inout), target :: vector

  !Local variables
  integer :: i, natoms
  real(kind=8), dimension(:), pointer :: x, y, z     ! Pointers for coordinates
  real(kind=8) :: xtotal, ytotal, ztotal

  natoms = vecsize / 3

  ! We first set-up pointers for the x, y, z components 
  x => vector(1:natoms)
  y => vector(natoms+1:2*natoms)
  z => vector(2*natoms+1:3*natoms)

  xtotal = 0.0d0
  ytotal = 0.0d0
  ztotal = 0.0d0

  do i = 1, natoms
     xtotal = xtotal + x(i)
     ytotal = ytotal + y(i)
     ztotal = ztotal + z(i)
  enddo 

  xtotal = xtotal / natoms
  ytotal = ytotal / natoms
  ztotal = ztotal / natoms

  do i = 1, natoms
     x(i) = x(i) - xtotal
     y(i) = y(i) - ytotal
     z(i) = z(i) - ztotal
  end do

END SUBROUTINE center
!!***

!!****f* utils/displacement
!! FUNCTION
!!   It computes the distance between two configurations and 
!!   the number of particles having moved by more than a THRESHOLD
!!
!! SOURCE
!!
subroutine displacement( posa, posb, delr, npart )

  use defs
  implicit none

  !Arguments
  real(kind=8), dimension(vecsize), intent(in), target :: posa
  real(kind=8), dimension(vecsize), intent(in), target :: posb 
  real(kind=8), intent(out)                            :: delr
  integer, intent(out)                                 :: npart

  !Local variables
  real(kind=8), parameter :: THRESHOLD = 0.1  ! In Angstroems
  real(kind=8), dimension(:), pointer :: xa, ya, za, xb, yb, zb
  integer :: i, j
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


!!****f* utils/store
!! FUNCTION
!!   Subroutine store
!!   This subroutine stores the configurations at minima and activated points
!!   By definition, it uses pos, box and scala
!! SOURCE
!! 
subroutine store( fname )

  use defs
  implicit none

  !Arguments
  character(len=7 ), intent(in) :: fname

  !Local variables
  integer :: ierror
  integer :: i
  real(kind=8) :: boxl
  character(len=*), parameter :: extension = ".xyz"
  character(len=20) :: fnamexyz

  ! We first set-up pointers for the x, y, z components for posa and posb

  boxl = box * scala  ! Update the box size

! added by Fedwa El-Mellouhi July 2002, writes the configuration in jmol format 
  fnamexyz = trim(fname) // extension

  write(*,*) 'BART: Writing to file : ', fname
  write(*,*) 'BART: Writing to file : ', fnamexyz
   
  open(unit=FCONF,file=fname,status='unknown',action='write',iostat=ierror)
  open(unit=XYZ,file=fnamexyz,status='unknown',action='write',iostat=ierror)

  write(XYZ,*) NATOMS , 'angstroem' 
  write(XYZ,*) boxl

  write(FCONF,*) 'run_id: ', mincounter
  write(FCONF,*) 'total energy : ', total_energy
  write(FCONF,*) boxl

  do i=1, NATOMS
     write(XYZ,'(1x,A2,3(2x,f16.8))')   Atom(i), x(i), y(i), z(i)
     write(FCONF,'(1x,i6,3(2x,f16.8))') type(i), x(i), y(i), z(i)
  end do

  close(FCONF)
  close(XYZ)

END SUBROUTINE store
!!***


!!****f* utils/write_refconfig
!! FUNCTION
!!   This subroutine writes the atomic positions and others to a "refconfig" file
!!   which will be used a the reference point until a new events gets accepted.
!!
!! SOURCE
!!
subroutine write_refconfig( )

  use defs
  implicit none

  !Local variables
  integer :: i, ierror
  real(kind=8) :: boxl

  boxl = box * scala                  ! Update the box size 

                                      !switch replace for unknown
  open(unit=FREFCONFIG,file=REFCONFIG,status='unknown',action='write',iostat=ierror) 
  write(FREFCONFIG,*) 'run_id: ', mincounter
  write(FREFCONFIG,*) 'total energy : '
  write(FREFCONFIG,*) total_energy
  write(FREFCONFIG,*) boxl
  do i = 1, NATOMS
    write(FREFCONFIG,'(1x,i6, 3(2x,F16.8))') type(i), x(i), y(i), z(i)
  end do
  close(FREFCONFIG)

END SUBROUTINE write_refconfig
!!***


!!****f* utils/store_part
!! FUNCTION
!!   This subroutine stores partial configurations. 
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
  real(kind=8) :: boxl
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
  write(XYZ,*) boxl

  do i= 1, NATOMS
     write(XYZ,'(1x,A2,3(2x,f16.8))')   Atom(i), x(i), y(i), z(i)
  end do

  write(XYZ,*) '# simulation ',scounter," stage ", stage, " iteration ",rcounter
  write(XYZ,'(a,(1p,e17.10,0p))') ' # total energy (eV) : ', total_energy

  close(XYZ)

END SUBROUTINE store_part
!!***


!!****f* utils/convert_to_chain2
!!   The subroutine convert_to_chain takes an integer and transforms it into a
!!   chain of character. It is the same convert_to_chain, but with a small
!!   modification.
!!
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
     chain = '0'
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


!!****f* utils/save_intermediate
!! FUNCTION
!!   It saves the configuration at every step in xyz format.
!!   The name of the file will be:
!!    conf_1001_K_030.xyz
!!   if the 'mincounter' is 1001, 'stage' is 'K', and 'adv',
!!   i.e the step, is 30.
!!
!! SOURCE
!! 
subroutine save_intermediate( stage, adv )

  use defs
  implicit none

  !Arguments
  character(len=1), intent(in) :: stage
  integer,          intent(in) :: adv 

  !Local variables
  character(len=20) :: fname
  character(len=4)  :: rcounter, scounter

                                      ! subroutines in utils.f90 
  if ( iproc == 0 ) then 

     call convert_to_chain( mincounter, scounter )
     call convert_to_chain_2( adv, rcounter )
     fname = 'conf_'//trim(scounter)//'_'//stage//'_'//trim(rcounter)
     call store_part( fname, scounter, rcounter, stage )

  end if

END SUBROUTINE save_intermediate
!!***
