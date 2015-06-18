!> @file
!!  Routines to read and print input variables
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 



!!!!> Create the directory output
!!!subroutine create_dir_output(iproc, in)
!!!  use yaml_output
!!!  use module_types
!!!  use module_base
!!!  implicit none
!!!  integer, intent(in) :: iproc
!!!  type(input_variables), intent(inout) :: in
!!!
!!!  integer(kind=4),parameter    :: dirlen=100
!!!  character(len=dirlen) :: dirname
!!!  integer :: ierror
!!!  integer(kind=4) :: i_stat, ierr
!!!  integer         :: ierrr
!!!
!!!  ! Create a directory to put the files in.
!!!  !dirname=repeat(' ',len(dirname))
!!!  call f_zero(dirname)
!!!  if (iproc == 0) then
!!!     call f_mkdir(in%dir_output,dirname)
!!!!!$     call getdir(in%dir_output, int(len_trim(in%dir_output),kind=4), dirname, dirlen, i_stat)
!!!!!$     if (i_stat /= 0) then
!!!!!$        call yaml_warning("Cannot create output directory '" // trim(in%dir_output) // "'.")
!!!!!$        call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierrr)
!!!!!$     end if
!!!  end if
!!!  call mpibcast(dirname,comm=bigdft_mpi%mpi_comm)
!!!  !call MPI_BCAST(dirname,len(dirname),MPI_CHARACTER,0,bigdft_mpi%mpi_comm,ierrr)
!!!  in%dir_output=dirname
!!!END SUBROUTINE create_dir_output
