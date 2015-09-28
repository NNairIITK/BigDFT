!> @file
!!    Modulefile for handling of the Simulation box of the Poisson Solver
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module PSbox
  use wrapper_MPI
  use PStypes, only: coulomb_operator, PSolver_energies
  use PSbase
  implicit none
  private

  interface PS_reduce
     module procedure reduce_scalar,reduce_array,reduce_energies
  end interface PS_reduce

  public :: PS_reduce,PS_gather

contains

  !>gather a distributed array to have a full array
  !!if only src is present this is assumed to be a full array
  !!otherwise it is assumed to be a distributed array
  subroutine PS_gather(src,kernel,dest)
    use dynamic_memory, only: f_memcpy
    implicit none
    !> input array. If dest is present, the values are assumed to be distributed
    !!otherwise the values are not modified and are gathered in the dest
    !!array
    real(dp), dimension(*), intent(inout) :: src
    type(coulomb_operator), intent(in) :: kernel
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(out), optional :: dest

    if (present(dest)) then
       if (kernel%mpi_env%nproc > 1) then
          call mpiallgather(src(1),recvbuf=dest(1,1,1),&
               recvcounts=kernel%counts,&
               displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
       else
          call f_memcpy(n=size(dest),src=src(1),dest=dest(1,1,1))
       end if
    else
       if (kernel%mpi_env%nproc > 1) then
          call mpiallgather(src(1),recvcounts=kernel%counts,&
               displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
       end if
    end if
  end subroutine PS_gather

  !> reduce all the given information 
  !! MPI_SUM is applied in the case of unspecified op
  subroutine reduce_scalar(val,kernel,op)
    implicit none
    real(dp), intent(inout) :: val
    type(coulomb_operator), intent(in) :: kernel
    integer, intent(in), optional :: op !< operation to be done
    !local variables
    integer :: mpi_op

    mpi_op=MPI_SUM
    if (present(op)) mpi_op=op

    if (kernel%mpi_env%nproc > 1) &
         call mpiallred(val,1,op=mpi_op,comm=kernel%mpi_env%mpi_comm)
    
  end subroutine reduce_scalar

  subroutine reduce_array(val,kernel,op)
    implicit none
    real(dp), dimension(:), intent(inout) :: val
    type(coulomb_operator), intent(in) :: kernel
    integer, intent(in), optional :: op !< operation to be done
    !local variables
    integer :: mpi_op

    mpi_op=MPI_SUM
    if (present(op)) mpi_op=op

    if (kernel%mpi_env%nproc > 1) &
         call mpiallred(val(1),size(val),op=mpi_op,comm=kernel%mpi_env%mpi_comm)

  end subroutine reduce_array

  !>this is of course to do the sum
  subroutine reduce_energies(e,kernel)
    type(PSolver_energies), intent(inout) :: e
    type(coulomb_operator), intent(in) :: kernel
    !local variables
    integer, parameter :: energsize=10
    real(gp), dimension(energsize) :: vals

    if (kernel%mpi_env%nproc > 1) then
       vals(1)   =e%hartree   
       vals(2)   =e%elec      
       vals(3)   =e%eVextra   
       vals(4)   =e%cavitation
       vals(5:10)=e%strten     
       call PS_reduce(vals,kernel)
       e%hartree   =vals(1)   
       e%elec      =vals(2)   
       e%eVextra   =vals(3)   
       e%cavitation=vals(4)   
       e%strten    =vals(5:10)
    end if

  end subroutine reduce_energies
end module PSbox
