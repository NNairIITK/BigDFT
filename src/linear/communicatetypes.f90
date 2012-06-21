!> @file 
!!   Routines to communicate types
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


module module_getbounds
   implicit none
   interface
      subroutine getbounds(iproc, root, array, is1, ie1, is2, ie2, is3, ie3)
         use module_base
         implicit none
         integer,intent(in):: iproc, root
         integer,dimension(:,:,:),pointer,intent(in):: array
         integer,intent(out):: is1, ie1, is2, ie2, is3, ie3
      END SUBROUTINE getbounds
   end interface
END MODULE module_getbounds


module module_communicatetypes

   implicit none

   interface

      subroutine communicate_locreg_descriptors(iproc, root, llr)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, root
         type(locreg_descriptors),intent(inout):: llr
      END SUBROUTINE communicate_locreg_descriptors

      subroutine communicate_grid_dimensions(iproc, root, d)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, root
         type(grid_dimensions),intent(inout):: d
      END SUBROUTINE communicate_grid_dimensions

      subroutine communicate_wavefunctions_descriptors(iproc, root, wfd)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, root
         type(wavefunctions_descriptors),intent(inout):: wfd
      END SUBROUTINE communicate_wavefunctions_descriptors

      subroutine communicate_convolutions_bounds(iproc, root, bounds)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, root
         type(convolutions_bounds),intent(inout):: bounds
      END SUBROUTINE communicate_convolutions_bounds

      subroutine communicate_kinetic_bounds(iproc, root, kb)
         use module_base
         use module_types
         use module_getbounds
         implicit none
         integer,intent(in):: iproc, root
         type(kinetic_bounds),intent(inout):: kb
      END SUBROUTINE communicate_kinetic_bounds

      subroutine communicate_shrink_bounds(iproc, root, sb)
         use module_base
         use module_types
         use module_getbounds
         implicit none
         integer,intent(in):: iproc, root
         type(shrink_bounds),intent(inout):: sb
      END SUBROUTINE communicate_shrink_bounds

      subroutine communicate_grow_bounds(iproc, root, gb)
         use module_base
         use module_types
         use module_getbounds
         implicit none
         integer,intent(in):: iproc, root
         type(grow_bounds),intent(inout):: gb
      END SUBROUTINE communicate_grow_bounds

   end interface
END MODULE module_communicatetypes




subroutine communicate_locreg_descriptors(iproc, root, llr)
   use module_base
   use module_types
   use module_communicatetypes, except_this_one => communicate_locreg_descriptors
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(locreg_descriptors),intent(inout):: llr

   ! Local variables
   integer:: ierr
   !integer:: commtype
   !integer,dimension(12):: types
   integer,dimension(10):: temparr_int
   real(8),dimension(4):: temparr_dbl

   

   !!types = (/ mpi_character, mpi_logical, mpi_integer, mpi_integer, mpi_integer, mpi_integer, &
   !!          mpi_integer, mpi_integer, mpi_integer, mpi_integer, mpi_double_precision, mpi_double_precision /)


   !!! First communicate all scalars and fixed-size arrays
   !!call mpi_bcast(llr%geocode, 1, mpi_character, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%hybrid_on, 1, mpi_logical, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%ns1, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%ns2, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%ns3, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%nsi1, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%nsi2, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%nsi3, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%localnorb, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%outofzone, 3, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%locregCenter, 3, mpi_double_precision, root, mpi_comm_world, ierr)
   !!call mpi_bcast(llr%locrad, 1, mpi_double_precision, root, mpi_comm_world, ierr)

   temparr_int(1)=llr%ns1
   temparr_int(2)=llr%ns2
   temparr_int(3)=llr%ns3
   temparr_int(4)=llr%nsi1
   temparr_int(5)=llr%nsi2
   temparr_int(6)=llr%nsi3
   temparr_int(7)=llr%localnorb
   temparr_int(8)=llr%outofzone(1)
   temparr_int(9)=llr%outofzone(2)
   temparr_int(10)=llr%outofzone(3)
   temparr_dbl(1)=llr%locregCenter(1)
   temparr_dbl(2)=llr%locregCenter(2)
   temparr_dbl(3)=llr%locregCenter(3)
   temparr_dbl(4)=llr%locrad
   call mpi_bcast(llr%geocode, 1, mpi_character, root, mpi_comm_world, ierr)
   call mpi_bcast(llr%hybrid_on, 1, mpi_logical, root, mpi_comm_world, ierr)
   call mpi_bcast(temparr_int(1), 10, mpi_integer, root, mpi_comm_world, ierr)
   call mpi_bcast(temparr_dbl(1), 4, mpi_double_precision, root, mpi_comm_world, ierr)
   llr%ns1=temparr_int(1)
   llr%ns2=temparr_int(2)
   llr%ns3=temparr_int(3)
   llr%nsi1=temparr_int(4)
   llr%nsi2=temparr_int(5)
   llr%nsi3=temparr_int(6)
   llr%localnorb=temparr_int(7)
   llr%outofzone(1)=temparr_int(8)
   llr%outofzone(2)=temparr_int(9)
   llr%outofzone(3)=temparr_int(10)
   llr%locregCenter(1)=temparr_dbl(1)
   llr%locregCenter(2)=temparr_dbl(2)
   llr%locregCenter(3)=temparr_dbl(3)
   llr%locrad=temparr_dbl(4)


   ! Now communicate the types
   call communicate_grid_dimensions(iproc, root, llr%d)
   call communicate_wavefunctions_descriptors(iproc, root, llr%wfd)
   if (llr%geocode == 'F') call communicate_convolutions_bounds(iproc, root, llr%bounds)

END SUBROUTINE communicate_locreg_descriptors



subroutine communicate_grid_dimensions(iproc, root, d)
   use module_base
   use module_types
   use module_communicatetypes, except_this_one => communicate_grid_dimensions
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(grid_dimensions),intent(inout):: d

   ! Local variables
   integer:: ierr
   integer,dimension(12):: temparr_int


   !!call mpi_bcast(d%n1, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%n2, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%n3, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%nfl1, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%nfu1, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%nfl2, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%nfu2, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%nfl3, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%nfu3, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%n1i, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%n2i, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(d%n3i, 1, mpi_integer, root, mpi_comm_world, ierr)

   temparr_int(1)=d%n1
   temparr_int(2)=d%n2
   temparr_int(3)=d%n3
   temparr_int(4)=d%nfl1
   temparr_int(5)=d%nfu1
   temparr_int(6)=d%nfl2
   temparr_int(7)=d%nfu2
   temparr_int(8)=d%nfl3
   temparr_int(9)=d%nfu3
   temparr_int(10)=d%n1i
   temparr_int(11)=d%n2i
   temparr_int(12)=d%n3i
   call mpi_bcast(temparr_int(1), 12, mpi_integer, root, mpi_comm_world, ierr)
   d%n1=temparr_int(1)
   d%n2=temparr_int(2)
   d%n3=temparr_int(3)
   d%nfl1=temparr_int(4)
   d%nfu1=temparr_int(5)
   d%nfl2=temparr_int(6)
   d%nfu2=temparr_int(7)
   d%nfl3=temparr_int(8)
   d%nfu3=temparr_int(9)
   d%n1i=temparr_int(10)
   d%n2i=temparr_int(11)
   d%n3i=temparr_int(12)

END SUBROUTINE communicate_grid_dimensions




subroutine communicate_wavefunctions_descriptors(iproc, root, wfd)
   use module_base
   use module_types
   use module_communicatetypes, except_this_one => communicate_wavefunctions_descriptors
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(wavefunctions_descriptors),intent(inout):: wfd

   ! Local variables
   integer:: ierr
   character(len=*),parameter:: subname='communicate_wavefunctions_descriptors'
   integer,dimension(4):: temparr_int

   ! First communicate all scalars
   !!call mpi_bcast(wfd%nvctr_c, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(wfd%nvctr_f, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(wfd%nseg_c, 1, mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(wfd%nseg_f, 1, mpi_integer, root, mpi_comm_world, ierr)

   temparr_int(1)=wfd%nvctr_c
   temparr_int(2)=wfd%nvctr_f
   temparr_int(3)=wfd%nseg_c
   temparr_int(4)=wfd%nseg_f
   call mpi_bcast(temparr_int(1), 4, mpi_integer, root, mpi_comm_world, ierr)
   wfd%nvctr_c=temparr_int(1)
   wfd%nvctr_f=temparr_int(2)
   wfd%nseg_c=temparr_int(3)
   wfd%nseg_f=temparr_int(4)

   ! Allocate the arrays
   if(iproc/=root) call allocate_wfd(wfd,subname)

   ! Communicate the arrays
   call mpi_bcast(wfd%keyglob, 2*(wfd%nseg_c+wfd%nseg_f), mpi_integer, root, mpi_comm_world, ierr)
   call mpi_bcast(wfd%keygloc, 2*(wfd%nseg_c+wfd%nseg_f), mpi_integer, root, mpi_comm_world, ierr)
   call mpi_bcast(wfd%keyvloc, wfd%nseg_c+wfd%nseg_f, mpi_integer, root, mpi_comm_world, ierr)
   call mpi_bcast(wfd%keyvglob, wfd%nseg_c+wfd%nseg_f, mpi_integer, root, mpi_comm_world, ierr)



END SUBROUTINE communicate_wavefunctions_descriptors




subroutine communicate_convolutions_bounds(iproc, root, bounds)
   use module_base
   use module_types
   use module_getbounds
   use module_communicatetypes, except_this_one => communicate_convolutions_bounds
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(convolutions_bounds),intent(inout):: bounds

   ! Local variables
   integer:: is1, ie1, is2, ie2, is3, ie3, istat, ierr
   character(len=*),parameter:: subname='communicate_convolutions_bounds'

   call communicate_kinetic_bounds(iproc, root, bounds%kb)
   call communicate_shrink_bounds(iproc, root, bounds%sb)
   call communicate_grow_bounds(iproc, root, bounds%gb)
   call getbounds(iproc, root, bounds%ibyyzz_r, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(bounds%ibyyzz_r(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, bounds%ibyyzz_r, 'bounds%ibyyzz_r,', subname)
   end if
   call mpi_bcast(bounds%ibyyzz_r, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)


END SUBROUTINE communicate_convolutions_bounds




subroutine communicate_kinetic_bounds(iproc, root, kb)
   use module_base
   use module_types
   use module_getbounds
   use module_communicatetypes, except_this_one => communicate_kinetic_bounds
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(kinetic_bounds),intent(inout):: kb

   ! Local variables
   integer:: ierr, istat, is1, ie1, is2, ie2, is3, ie3
   character(len=*),parameter:: subname='communicate_kinetic_bounds'

   call getbounds(iproc, root, kb%ibyz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibyz_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibyz_c, 'kb%ibyz_c', subname)
   end if
   call mpi_bcast(kb%ibyz_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, kb%ibxz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxz_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxz_c, 'kb%ibxz_c', subname)
   end if
   call mpi_bcast(kb%ibxz_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, kb%ibxy_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxy_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxy_c, 'kb%ibxy_c', subname)
   end if
   call mpi_bcast(kb%ibxy_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, kb%ibyz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibyz_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibyz_f, 'kb%ibyz_f', subname)
   end if
   call mpi_bcast(kb%ibyz_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, kb%ibxz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxz_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxz_f, 'kb%ibxz_f', subname)
   end if
   call mpi_bcast(kb%ibxz_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, kb%ibxy_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxy_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxy_f, 'kb%ibxy_f', subname)
   end if
   call mpi_bcast(kb%ibxy_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)



END SUBROUTINE communicate_kinetic_bounds


subroutine communicate_shrink_bounds(iproc, root, sb)
   use module_base
   use module_types
   use module_getbounds
   use module_communicatetypes, except_this_one => communicate_shrink_bounds
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(shrink_bounds),intent(inout):: sb

   ! Local variables
   integer:: ierr, istat, is1, ie1, is2, ie2, is3, ie3
   character(len=*),parameter:: subname='communicate_shrink_bounds'

   call getbounds(iproc, root, sb%ibzzx_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibzzx_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibzzx_c, 'sb%ibzzx_c', subname)
   end if

   call mpi_bcast(sb%ibzzx_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibyyzz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibyyzz_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibyyzz_c, 'sb%ibyyzz_c', subname)
   end if
   call mpi_bcast(sb%ibyyzz_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibxy_ff, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibxy_ff(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibxy_ff, 'sb%ibxy_ff', subname)
   end if
   call mpi_bcast(sb%ibxy_ff, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibzzx_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibzzx_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibzzx_f, 'sb%ibzzx_f', subname)
   end if
   call mpi_bcast(sb%ibzzx_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibyyzz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibyyzz_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibyyzz_f, 'sb%ibyyzz_f', subname)
   end if
   call mpi_bcast(sb%ibyyzz_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)


END SUBROUTINE communicate_shrink_bounds



subroutine communicate_grow_bounds(iproc, root, gb)
   use module_base
   use module_types
   use module_getbounds
   use module_communicatetypes, except_this_one => communicate_grow_bounds
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   type(grow_bounds),intent(inout):: gb

   ! Local variables
   integer:: ierr, istat, is1, ie1, is2, ie2, is3, ie3
   character(len=*),parameter:: subname='communicate_shrink_bounds'

   call getbounds(iproc, root, gb%ibzxx_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibzxx_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, gb%ibzxx_c, 'gb%ibzxx_c', subname)
   end if
   call mpi_bcast(gb%ibzxx_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, gb%ibxxyy_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibxxyy_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, gb%ibxxyy_c, 'gb%ibxxyy_c', subname)
   end if
   call mpi_bcast(gb%ibxxyy_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, gb%ibyz_ff, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibyz_ff(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, gb%ibyz_ff, 'gb%ibyz_ff', subname)
   end if
   call mpi_bcast(gb%ibyz_ff, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, gb%ibzxx_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibzxx_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, gb%ibzxx_f, 'gb%ibzxx_f', subname)
   end if
   call mpi_bcast(gb%ibzxx_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, gb%ibxxyy_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibxxyy_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, gb%ibxxyy_f, 'gb%ibxxyy_f', subname)
   end if
   call mpi_bcast(gb%ibxxyy_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

END SUBROUTINE communicate_grow_bounds



subroutine getbounds(iproc, root, array, is1, ie1, is2, ie2, is3, ie3)
   use module_base
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   integer,dimension(:,:,:),pointer,intent(in):: array
   integer,intent(out):: is1, ie1, is2, ie2, is3, ie3

   ! Local variables
   integer,dimension(6):: ind
   integer:: ierr

   if(iproc==root) then
      ind(1)=lbound(array,1)
      ind(2)=ubound(array,1)
      ind(3)=lbound(array,2)
      ind(4)=ubound(array,2)
      ind(5)=lbound(array,3)
      ind(6)=ubound(array,3)
   end if

   call mpi_bcast(ind, 6, mpi_integer, root, mpi_comm_world, ierr)

   is1=ind(1)
   ie1=ind(2)
   is2=ind(3)
   ie2=ind(4)
   is3=ind(5)
   ie3=ind(6)


END SUBROUTINE getbounds
