!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

module module_global_variables
    use module_base, only: gp !bigdft base module
    implicit none
    character(len = *), public, parameter :: mhgps_version   = '0.01'
    character(len = *), public, parameter :: inputdir   = 'input'
    character(len = *), public, parameter :: outputdir   = 'input'

!    !others
!    integer, save               :: nbond = 1
    integer, save, allocatable  :: iconnect(:,:) 
    real(gp), save, allocatable :: minmode(:,:)
!    character(len=60), save     :: saddle_filename='saddle.mon'
    logical, save               :: isForceField=.false.
    real(gp), save              :: ef_counter=0.d0
    character(len=8), save      :: currDir
!    character(len=3), parameter :: prefix='pos'
    character(len=5), save      :: isadc
    integer, save               :: isad
    integer, save               :: ntodo
    character(len=5), save      :: isadprobc
    integer , save              :: isadprob=0

    integer, save :: inputPsiId=0
    integer, save :: iproc=0,nproc=1,igroup=0,ngroups=1
    integer, save :: itermin=0
    real(gp), save :: frac_fluct=0.d0


end module
