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
      subroutine getbounds_5(iproc, root, arr1, arr2, arr3, arr4, arr5, ise)
         use module_base
         implicit none
         integer,intent(in):: iproc, root
         integer,dimension(:,:,:),pointer,intent(in):: arr1, arr2, arr3, arr4, arr5
         integer,dimension(6,5),intent(out):: ise
      END SUBROUTINE getbounds_5
      subroutine getbounds_6(iproc, root, arr1, arr2, arr3, arr4, arr5, arr6, ise)
         use module_base
         implicit none
         integer,intent(in):: iproc, root
         integer,dimension(:,:,:),pointer,intent(in):: arr1, arr2, arr3, arr4, arr5, arr6
         integer,dimension(6,6),intent(out):: ise
      END SUBROUTINE getbounds_6
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

   !!! Local type
   !!type: local_scalars
   !!  character(len=1) :: geocode
   !!  logical :: hybrid_on
   !!  integer :: ns1, ns2, ns3, nsi1, nsi2, nsi3, Localnorb
   !!  integer,dimension(3) :: outofzone
   !!  real(8),dimension(3):: locregCenter
   !!  real(8):: locrad
   !!end type local_scalars

   ! Local variables
   integer:: ierr, ncount, mpi_tmptype
   integer,dimension(12):: blocklengths, types
   integer(kind=mpi_address_kind),dimension(12):: dspls
   integer(kind=mpi_address_kind):: addr_geocode, addr_hybrid_on, addr_ns1, addr_ns2, addr_ns3, addr_nsi1, addr_nsi2
   integer(kind=mpi_address_kind):: addr_nsi3, addr_localnorb, addr_outofzone, addr_locregCenter, addr_locrad, addr_llr

   
   !!! Copy data to datatype
   !!tmptype%geocode = llr%geocode
   !!tmptype%hybrid_on = llr%hybrid_on
   !!tmptype%ns1 = llr%ns1
   !!tmptype%ns2 = llr%ns2
   !!tmptype%ns3 = llr%ns3
   !!tmptype%nsi1 = llr%nsi1
   !!tmptype%nsi2 = llr%nsi2
   !!tmptype%nsi3 = llr%nsi3
   !!tmptype%localnorb = llr%localnorb
   !!tmptype%outofzone = llr%outofzone
   !!tmptype%locregCenter = llr%locregCenter
   !!tmptype%locrad = llr%locrad
   
   ! Build MPI datatype
   ncount=12
   blocklengths=(/1,1,1,1,1,1,1,1,1,3,3,1/)
   call mpi_get_address(llr, addr_llr, ierr)
   call mpi_get_address(llr%geocode, addr_geocode, ierr)
   call mpi_get_address(llr%hybrid_on, addr_hybrid_on, ierr)
   call mpi_get_address(llr%ns1,  addr_ns1, ierr)
   call mpi_get_address(llr%ns2,  addr_ns2, ierr)
   call mpi_get_address(llr%ns3,  addr_ns3, ierr)
   call mpi_get_address(llr%nsi1, addr_nsi1, ierr)
   call mpi_get_address(llr%nsi2, addr_nsi2, ierr)
   call mpi_get_address(llr%nsi3, addr_nsi3, ierr)
   call mpi_get_address(llr%localnorb, addr_localnorb, ierr)
   call mpi_get_address(llr%outofzone, addr_outofzone, ierr)
   call mpi_get_address(llr%locregCenter, addr_locregCenter, ierr)
   call mpi_get_address(llr%locrad, addr_locrad, ierr)
   
   dspls(1) = addr_geocode - addr_llr
   dspls(2) = addr_hybrid_on - addr_llr
   dspls(3) = addr_ns1 - addr_llr
   dspls(4) = addr_ns2 - addr_llr
   dspls(5) = addr_ns3 - addr_llr
   dspls(6) = addr_nsi1 - addr_llr
   dspls(7) = addr_nsi2 - addr_llr
   dspls(8) = addr_nsi3 - addr_llr
   dspls(9) = addr_localnorb - addr_llr
   dspls(10) = addr_outofzone - addr_llr
   dspls(11) = addr_locregCenter - addr_llr
   dspls(12) = addr_locrad - addr_llr

   types = (/mpi_character, mpi_logical, mpi_integer, mpi_integer, mpi_integer, mpi_integer, &
             mpi_integer, mpi_integer, mpi_integer, mpi_integer, mpi_double_precision, mpi_double_precision/)

   call mpi_type_create_struct(ncount, blocklengths, dspls, types, mpi_tmptype, ierr)
   !call mpi_type_create_struct(1, 1, 0, mpi_character, mpi_tmptype, ierr)
   call mpi_type_commit(mpi_tmptype, ierr)
   call mpi_bcast(llr, 1, mpi_tmptype, root, mpi_comm_world, ierr)
   call mpi_type_free(mpi_tmptype, ierr)


   ! Now communicate the types
   call communicate_grid_dimensions(iproc, root, llr%d)
   call communicate_wavefunctions_descriptors(iproc, root, llr%wfd)
   !!if (llr%geocode == 'F') call communicate_convolutions_bounds(iproc, root, llr%bounds)

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
   integer:: ierr, grid_dimensions_type

   call mpi_type_contiguous(12, mpi_integer, grid_dimensions_type, ierr)
   call mpi_type_commit(grid_dimensions_type, ierr)
   call mpi_bcast(d, 1, grid_dimensions_type, root, mpi_comm_world, ierr)
   call mpi_type_free(grid_dimensions_type, ierr)


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
   integer:: ierr, ncount, commtype
   character(len=*),parameter:: subname='communicate_wavefunctions_descriptors'
   integer,dimension(4):: blocklengths,types
   integer(kind=mpi_address_kind):: addr_wfd, addr_nvctr_c, addr_nvctr_f, addr_nseg_c, addr_nseg_f
   integer(kind=mpi_address_kind):: addr_keyglob, addr_keygloc, addr_keyvloc, addr_keyvglob
   integer(kind=mpi_address_kind),dimension(4):: dspls

   ncount=4 
   blocklengths=(/1,1,1,1/)
   call mpi_get_address(wfd, addr_wfd, ierr)
   call mpi_get_address(wfd%nvctr_c, addr_nvctr_c, ierr)
   call mpi_get_address(wfd%nvctr_f, addr_nvctr_f, ierr)
   call mpi_get_address(wfd%nseg_c, addr_nseg_c, ierr)
   call mpi_get_address(wfd%nseg_f, addr_nseg_f, ierr)

   dspls(1) = addr_nvctr_c - addr_wfd
   dspls(2) = addr_nvctr_f - addr_wfd
   dspls(3) = addr_nseg_c - addr_wfd
   dspls(4) = addr_nseg_f - addr_wfd

   types = (/mpi_integer, mpi_integer, mpi_integer, mpi_integer/)
   
   call mpi_type_create_struct(ncount, blocklengths, dspls, types, commtype, ierr)
   call mpi_type_commit(commtype, ierr)
   call mpi_bcast(wfd, 1, commtype, root, mpi_comm_world, ierr)
   call mpi_type_free(commtype, ierr)

   ! Allocate the arrays
   if(iproc/=root) call allocate_wfd(wfd,subname)

   ncount=4 
   blocklengths=(/2*(wfd%nseg_c+wfd%nseg_f), 2*(wfd%nseg_c+wfd%nseg_f), wfd%nseg_c+wfd%nseg_f, wfd%nseg_c+wfd%nseg_f/)
   call mpi_get_address(wfd, addr_wfd, ierr)
   call mpi_get_address(wfd%keyglob, addr_keyglob, ierr)
   call mpi_get_address(wfd%keygloc, addr_keygloc, ierr)
   call mpi_get_address(wfd%keyvloc, addr_keyvloc, ierr)
   call mpi_get_address(wfd%keyvglob, addr_keyvglob, ierr)

   dspls(1) = addr_keyglob - addr_wfd
   dspls(2) = addr_keygloc - addr_wfd
   dspls(3) = addr_keyvloc - addr_wfd
   dspls(4) = addr_keyvglob - addr_wfd

   types = (/mpi_integer, mpi_integer, mpi_integer, mpi_integer/)

   call mpi_type_create_struct(ncount, blocklengths, dspls, types, commtype, ierr)
   call mpi_type_commit(commtype, ierr)


   ! Communicate the arrays
   call mpi_bcast(wfd, 1, commtype, root, mpi_comm_world, ierr)
   call mpi_type_free(commtype, ierr)


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
   integer:: commtype
   integer,parameter:: ncount=6
   integer,dimension(ncount):: types, blocklengths
   integer,dimension(6,6):: ise
   integer(kind=mpi_address_kind):: addr_kb, addr_ibyz_c, addr_ibxz_c, addr_ibxy_c, addr_ibyz_f, addr_ibxz_f, addr_ibxy_f
   integer(kind=mpi_address_kind),dimension(ncount):: dspls

   call getbounds_6(iproc, root, kb%ibyz_c, kb%ibxz_c, kb%ibxy_c, kb%ibyz_f, kb%ibxz_f, kb%ibxy_f, ise)
   call mpi_barrier(mpi_comm_world, ierr)

   !call getbounds(iproc, root, kb%ibyz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibyz_c(ise(1,1):ise(2,1),ise(3,1):ise(4,1),ise(5,1):ise(6,1)), stat=istat)
      call memocc(istat, kb%ibyz_c, 'kb%ibyz_c', subname)
   end if
   blocklengths(1) = (ise(2,1)-ise(1,1)+1)*(ise(4,1)-ise(3,1)+1)*(ise(6,1)-ise(5,1)+1)
   types(1) = mpi_integer

   !call getbounds(iproc, root, kb%ibxz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxz_c(ise(1,2):ise(2,2),ise(3,2):ise(4,2),ise(5,2):ise(6,2)), stat=istat)
      call memocc(istat, kb%ibxz_c, 'kb%ibxz_c', subname)
   end if
   blocklengths(2) = (ise(2,2)-ise(1,2)+1)*(ise(4,2)-ise(3,2)+1)*(ise(6,2)-ise(5,2)+1)
   types(2) = mpi_integer


   !call getbounds(iproc, root, kb%ibxy_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxy_c(ise(1,3):ise(2,3),ise(3,3):ise(4,3),ise(5,3):ise(6,3)), stat=istat)
      call memocc(istat, kb%ibxy_c, 'kb%ibxy_c', subname)
   end if
   blocklengths(3) = (ise(2,3)-ise(1,3)+1)*(ise(4,3)-ise(3,3)+1)*(ise(6,3)-ise(5,3)+1)
   types(3) = mpi_integer


   !call getbounds(iproc, root, kb%ibyz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibyz_f(ise(1,4):ise(2,4),ise(3,4):ise(4,4),ise(5,4):ise(6,4)), stat=istat)
      call memocc(istat, kb%ibyz_f, 'kb%ibyz_f', subname)
   end if
   blocklengths(4) = (ise(2,4)-ise(1,4)+1)*(ise(4,4)-ise(3,4)+1)*(ise(6,4)-ise(5,4)+1)
   types(4) = mpi_integer

   !call getbounds(iproc, root, kb%ibxz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxz_f(ise(1,5):ise(2,5),ise(3,5):ise(4,5),ise(5,5):ise(6,5)), stat=istat)
      call memocc(istat, kb%ibxz_f, 'kb%ibxz_f', subname)
   end if
   blocklengths(5) = (ise(2,5)-ise(1,5)+1)*(ise(4,5)-ise(3,5)+1)*(ise(6,5)-ise(5,5)+1)
   types(5) = mpi_integer

   !call getbounds(iproc, root, kb%ibxy_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxy_f(ise(1,6):ise(2,6),ise(3,6):ise(4,6),ise(5,6):ise(6,6)), stat=istat)
      call memocc(istat, kb%ibxy_f, 'kb%ibxy_f', subname)
   end if
   blocklengths(6) = (ise(2,6)-ise(1,6)+1)*(ise(4,6)-ise(3,6)+1)*(ise(6,6)-ise(5,6)+1)
   types(6) = mpi_integer

   call mpi_get_address(kb, addr_kb, ierr)
   call mpi_get_address(kb%ibyz_c, addr_ibyz_c, ierr)
   dspls(1) = addr_ibyz_c - addr_kb
   call mpi_get_address(kb%ibxz_c, addr_ibxz_c, ierr)
   dspls(2) = addr_ibxz_c - addr_kb
   call mpi_get_address(kb%ibxy_c, addr_ibxy_c, ierr)
   dspls(3) = addr_ibxy_c - addr_kb
   call mpi_get_address(kb%ibyz_f, addr_ibyz_f, ierr)
   dspls(4) = addr_ibyz_f - addr_kb
   call mpi_get_address(kb%ibxz_f, addr_ibxz_f, ierr)
   dspls(5) = addr_ibxz_f - addr_kb
   call mpi_get_address(kb%ibxy_f, addr_ibxy_f, ierr)
   dspls(6) = addr_ibxy_f - addr_kb


   call mpi_type_create_struct(ncount, blocklengths, dspls, types, commtype, ierr)
   call mpi_type_commit(commtype, ierr)
   call mpi_bcast(kb, 1, commtype, root, mpi_comm_world, ierr)
   call mpi_type_free(commtype, ierr)

   !!call mpi_bcast(kb%ibyz_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(kb%ibxz_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(kb%ibxy_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(kb%ibyz_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(kb%ibxz_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)
   !!call mpi_bcast(kb%ibxy_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)



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
   integer:: ierr, istat, commtype
   character(len=*),parameter:: subname='communicate_shrink_bounds'
   integer,parameter:: ncount=5
   integer,dimension(ncount):: types, blocklengths
   integer,dimension(6,5):: ise
   integer(kind=mpi_address_kind):: addr_sb, addr_ibzzx_c, addr_ibyyzz_c, addr_ibxy_ff, addr_ibzzx_f, addr_ibyyzz_f
   integer(kind=mpi_address_kind),dimension(ncount):: dspls


   call getbounds_5(iproc, root, sb%ibzzx_c, sb%ibyyzz_c, sb%ibxy_ff, sb%ibzzx_f, sb%ibyyzz_f, ise)
   call mpi_barrier(mpi_comm_world, ierr)

   if(iproc/=root) then
      allocate(sb%ibzzx_c(ise(1,1):ise(2,1),ise(3,1):ise(4,1),ise(5,1):ise(6,1)), stat=istat)
      call memocc(istat, sb%ibzzx_c, 'sb%ibzzx_c', subname)
   end if
   blocklengths(1) = (ise(2,1)-ise(1,1)+1)*(ise(4,1)-ise(3,1)+1)*(ise(6,1)-ise(5,1)+1)
   types(1) = mpi_integer

   if(iproc/=root) then
      allocate(sb%ibyyzz_c(ise(1,2):ise(2,2),ise(3,2):ise(4,2),ise(5,2):ise(6,2)), stat=istat)
      call memocc(istat, sb%ibyyzz_c, 'sb%ibyyzz_c', subname)
   end if
   blocklengths(2) = (ise(2,2)-ise(1,2)+1)*(ise(4,2)-ise(3,2)+1)*(ise(6,2)-ise(5,2)+1)
   types(2) = mpi_integer

   if(iproc/=root) then
      allocate(sb%ibxy_ff(ise(1,3):ise(2,3),ise(3,3):ise(4,3),ise(5,3):ise(6,3)), stat=istat)
      call memocc(istat, sb%ibxy_ff, 'sb%ibxy_ff', subname)
   end if
   blocklengths(3) = (ise(2,3)-ise(1,3)+1)*(ise(4,3)-ise(3,3)+1)*(ise(6,3)-ise(5,3)+1)
   types(3) = mpi_integer

   if(iproc/=root) then
      allocate(sb%ibzzx_f(ise(1,4):ise(2,4),ise(3,4):ise(4,4),ise(5,4):ise(6,4)), stat=istat)
      call memocc(istat, sb%ibzzx_f, 'sb%ibzzx_f', subname)
   end if
   blocklengths(4) = (ise(2,4)-ise(1,4)+1)*(ise(4,4)-ise(3,4)+1)*(ise(6,4)-ise(5,4)+1)
   types(4) = mpi_integer

   if(iproc/=root) then
      allocate(sb%ibyyzz_f(ise(1,5):ise(2,5),ise(3,5):ise(4,5),ise(5,5):ise(6,5)), stat=istat)
      call memocc(istat, sb%ibyyzz_f, 'sb%ibyyzz_f', subname)
   end if
   blocklengths(5) = (ise(2,5)-ise(1,5)+1)*(ise(4,5)-ise(3,5)+1)*(ise(6,5)-ise(5,5)+1)
   types(5) = mpi_integer

   call mpi_get_address(sb, addr_sb, ierr)
   call mpi_get_address(sb%ibzzx_c, addr_ibzzx_c, ierr)
   dspls(1) = addr_ibzzx_c - addr_sb
   call mpi_get_address(sb%ibyyzz_c, addr_ibyyzz_c, ierr)
   dspls(2) = addr_ibyyzz_c - addr_sb
   call mpi_get_address(sb%ibxy_ff, addr_ibxy_ff, ierr)
   dspls(3) = addr_ibxy_ff - addr_sb
   call mpi_get_address(sb%ibzzx_f, addr_ibzzx_f, ierr)
   dspls(4) = addr_ibzzx_f - addr_sb
   call mpi_get_address(sb%ibyyzz_f, addr_ibyyzz_f, ierr)
   dspls(5) = addr_ibyyzz_f - addr_sb


   call mpi_type_create_struct(ncount, blocklengths, dspls, types, commtype, ierr)
   call mpi_type_commit(commtype, ierr)
   call mpi_bcast(sb, 1, commtype, root, mpi_comm_world, ierr)
   call mpi_type_free(commtype, ierr)


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
   integer:: ierr, istat
   character(len=*),parameter:: subname='communicate_shrink_bounds'
   integer:: commtype
   integer,parameter:: ncount=5
   integer,dimension(ncount):: types, blocklengths
   integer,dimension(6,5):: ise
   integer:: is1, ie1, is2, ie2, is3, ie3
   integer(kind=mpi_address_kind):: addr_gb, addr_ibzxx_c, addr_ibxxyy_c, addr_ibyz_ff, addr_ibzxx_f, addr_ibxxyy_f
   integer(kind=mpi_address_kind),dimension(ncount):: dspls

   call getbounds_5(iproc, root, gb%ibzxx_c, gb%ibxxyy_c, gb%ibyz_ff, gb%ibzxx_f, gb%ibxxyy_f, ise)
   call mpi_barrier(mpi_comm_world, ierr)

   !call getbounds(iproc, root, gb%ibzxx_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibzxx_c(ise(1,1):ise(2,1),ise(3,1):ise(4,1),ise(5,1):ise(6,1)), stat=istat)
      call memocc(istat, gb%ibzxx_c, 'gb%ibzxx_c', subname)
   end if
   blocklengths(1) = (ise(2,1)-ise(1,1)+1)*(ise(4,1)-ise(3,1)+1)*(ise(6,1)-ise(5,1)+1)
   types(1) = mpi_integer
   !call mpi_bcast(gb%ibzxx_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   !call getbounds(iproc, root, gb%ibxxyy_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibxxyy_c(ise(1,2):ise(2,2),ise(3,2):ise(4,2),ise(5,2):ise(6,2)), stat=istat)
      call memocc(istat, gb%ibxxyy_c, 'gb%ibxxyy_c', subname)
   end if
   blocklengths(2) = (ise(2,2)-ise(1,2)+1)*(ise(4,2)-ise(3,2)+1)*(ise(6,2)-ise(5,2)+1)
   types(2) = mpi_integer
   !call mpi_bcast(gb%ibxxyy_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   !call getbounds(iproc, root, gb%ibyz_ff, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibyz_ff(ise(1,3):ise(2,3),ise(3,3):ise(4,3),ise(5,3):ise(6,3)), stat=istat)
      call memocc(istat, gb%ibyz_ff, 'gb%ibyz_ff', subname)
   end if
   blocklengths(3) = (ise(2,3)-ise(1,3)+1)*(ise(4,3)-ise(3,3)+1)*(ise(6,3)-ise(5,3)+1)
   types(3) = mpi_integer
   !call mpi_bcast(gb%ibyz_ff, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   !call getbounds(iproc, root, gb%ibzxx_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibzxx_f(ise(1,4):ise(2,4),ise(3,4):ise(4,4),ise(5,4):ise(6,4)), stat=istat)
      call memocc(istat, gb%ibzxx_f, 'gb%ibzxx_f', subname)
   end if
   blocklengths(4) = (ise(2,4)-ise(1,4)+1)*(ise(4,4)-ise(3,4)+1)*(ise(6,4)-ise(5,4)+1)
   types(4) = mpi_integer
   !call mpi_bcast(gb%ibzxx_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   !call getbounds(iproc, root, gb%ibxxyy_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(gb%ibxxyy_f(ise(1,5):ise(2,5),ise(3,5):ise(4,5),ise(5,5):ise(6,5)), stat=istat)
      call memocc(istat, gb%ibxxyy_f, 'gb%ibxxyy_f', subname)
   end if
   blocklengths(5) = (ise(2,5)-ise(1,5)+1)*(ise(4,5)-ise(3,5)+1)*(ise(6,5)-ise(5,5)+1)
   types(5) = mpi_integer
   !call mpi_bcast(gb%ibxxyy_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)


   call mpi_get_address(gb, addr_gb, ierr)
   call mpi_get_address(gb%ibzxx_c, addr_ibzxx_c, ierr)
   dspls(1) = addr_ibzxx_c - addr_gb
   call mpi_get_address(gb%ibxxyy_c, addr_ibxxyy_c, ierr)
   dspls(2) = addr_ibxxyy_c - addr_gb
   call mpi_get_address(gb%ibyz_ff, addr_ibyz_ff, ierr)
   dspls(3) = addr_ibyz_ff - addr_gb
   call mpi_get_address(gb%ibzxx_f, addr_ibzxx_f, ierr)
   dspls(4) = addr_ibzxx_f - addr_gb
   call mpi_get_address(gb%ibxxyy_f, addr_ibxxyy_f, ierr)
   dspls(5) = addr_ibxxyy_f - addr_gb

   call mpi_type_create_struct(ncount, blocklengths, dspls, types, commtype, ierr)
   call mpi_type_commit(commtype, ierr)
   call mpi_bcast(gb, 1, commtype, root, mpi_comm_world, ierr)
   call mpi_type_free(commtype, ierr)

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

subroutine getbounds_5(iproc, root, arr1, arr2, arr3, arr4, arr5, ise)
   use module_base
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   integer,dimension(:,:,:),pointer,intent(in):: arr1, arr2, arr3, arr4, arr5
   integer,dimension(6,5),intent(out):: ise

   ! Local variables
   integer:: ierr

   if(iproc==root) then
      ise(1,1)=lbound(arr1,1)
      ise(2,1)=ubound(arr1,1)
      ise(3,1)=lbound(arr1,2)
      ise(4,1)=ubound(arr1,2)
      ise(5,1)=lbound(arr1,3)
      ise(6,1)=ubound(arr1,3)

      ise(1,2)=lbound(arr2,1)
      ise(2,2)=ubound(arr2,1)
      ise(3,2)=lbound(arr2,2)
      ise(4,2)=ubound(arr2,2)
      ise(5,2)=lbound(arr2,3)
      ise(6,2)=ubound(arr2,3)

      ise(1,3)=lbound(arr3,1)
      ise(2,3)=ubound(arr3,1)
      ise(3,3)=lbound(arr3,2)
      ise(4,3)=ubound(arr3,2)
      ise(5,3)=lbound(arr3,3)
      ise(6,3)=ubound(arr3,3)

      ise(1,4)=lbound(arr4,1)
      ise(2,4)=ubound(arr4,1)
      ise(3,4)=lbound(arr4,2)
      ise(4,4)=ubound(arr4,2)
      ise(5,4)=lbound(arr4,3)
      ise(6,4)=ubound(arr4,3)

      ise(1,5)=lbound(arr5,1)
      ise(2,5)=ubound(arr5,1)
      ise(3,5)=lbound(arr5,2)
      ise(4,5)=ubound(arr5,2)
      ise(5,5)=lbound(arr5,3)
      ise(6,5)=ubound(arr5,3)
   end if

   call mpi_bcast(ise(1,1), 30, mpi_integer, root, mpi_comm_world, ierr)


END SUBROUTINE getbounds_5


subroutine getbounds_6(iproc, root, arr1, arr2, arr3, arr4, arr5, arr6, ise)
   use module_base
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, root
   integer,dimension(:,:,:),pointer,intent(in):: arr1, arr2, arr3, arr4, arr5, arr6
   integer,dimension(6,6),intent(out):: ise

   ! Local variables
   integer:: ierr

   if(iproc==root) then
      ise(1,1)=lbound(arr1,1)
      ise(2,1)=ubound(arr1,1)
      ise(3,1)=lbound(arr1,2)
      ise(4,1)=ubound(arr1,2)
      ise(5,1)=lbound(arr1,3)
      ise(6,1)=ubound(arr1,3)

      ise(1,2)=lbound(arr2,1)
      ise(2,2)=ubound(arr2,1)
      ise(3,2)=lbound(arr2,2)
      ise(4,2)=ubound(arr2,2)
      ise(5,2)=lbound(arr2,3)
      ise(6,2)=ubound(arr2,3)

      ise(1,3)=lbound(arr3,1)
      ise(2,3)=ubound(arr3,1)
      ise(3,3)=lbound(arr3,2)
      ise(4,3)=ubound(arr3,2)
      ise(5,3)=lbound(arr3,3)
      ise(6,3)=ubound(arr3,3)

      ise(1,4)=lbound(arr4,1)
      ise(2,4)=ubound(arr4,1)
      ise(3,4)=lbound(arr4,2)
      ise(4,4)=ubound(arr4,2)
      ise(5,4)=lbound(arr4,3)
      ise(6,4)=ubound(arr4,3)

      ise(1,5)=lbound(arr5,1)
      ise(2,5)=ubound(arr5,1)
      ise(3,5)=lbound(arr5,2)
      ise(4,5)=ubound(arr5,2)
      ise(5,5)=lbound(arr5,3)
      ise(6,5)=ubound(arr5,3)

      ise(1,6)=lbound(arr6,1)
      ise(2,6)=ubound(arr6,1)
      ise(3,6)=lbound(arr6,2)
      ise(4,6)=ubound(arr6,2)
      ise(5,6)=lbound(arr6,3)
      ise(6,6)=ubound(arr6,3)
   end if

   call mpi_bcast(ise(1,1), 36, mpi_integer, root, mpi_comm_world, ierr)


END SUBROUTINE getbounds_6


subroutine get_convarrays_bounds(bounds, ab, mpi_ncounts)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(convolutions_bounds),intent(in):: bounds
  integer,dimension(6,6,4),intent(out):: ab
  integer(kind=mpi_address_kind),dimension(6,4),intent(out):: mpi_ncounts

  ab(1,1,1)=lbound(bounds%kb%ibyz_c,1)
  ab(2,1,1)=ubound(bounds%kb%ibyz_c,1)
  ab(3,1,1)=lbound(bounds%kb%ibyz_c,2)
  ab(4,1,1)=ubound(bounds%kb%ibyz_c,2)
  ab(5,1,1)=lbound(bounds%kb%ibyz_c,3)
  ab(6,1,1)=ubound(bounds%kb%ibyz_c,3)

  ab(1,2,1)=lbound(bounds%kb%ibxz_c,1)
  ab(2,2,1)=ubound(bounds%kb%ibxz_c,1)
  ab(3,2,1)=lbound(bounds%kb%ibxz_c,2)
  ab(4,2,1)=ubound(bounds%kb%ibxz_c,2)
  ab(5,2,1)=lbound(bounds%kb%ibxz_c,3)
  ab(6,2,1)=ubound(bounds%kb%ibxz_c,3)

  ab(1,3,1)=lbound(bounds%kb%ibxy_c,1)
  ab(2,3,1)=ubound(bounds%kb%ibxy_c,1)
  ab(3,3,1)=lbound(bounds%kb%ibxy_c,2)
  ab(4,3,1)=ubound(bounds%kb%ibxy_c,2)
  ab(5,3,1)=lbound(bounds%kb%ibxy_c,3)
  ab(6,3,1)=ubound(bounds%kb%ibxy_c,3)

  ab(1,4,1)=lbound(bounds%kb%ibyz_f,1)
  ab(2,4,1)=ubound(bounds%kb%ibyz_f,1)
  ab(3,4,1)=lbound(bounds%kb%ibyz_f,2)
  ab(4,4,1)=ubound(bounds%kb%ibyz_f,2)
  ab(5,4,1)=lbound(bounds%kb%ibyz_f,3)
  ab(6,4,1)=ubound(bounds%kb%ibyz_f,3)

  ab(1,5,1)=lbound(bounds%kb%ibxz_f,1)
  ab(2,5,1)=ubound(bounds%kb%ibxz_f,1)
  ab(3,5,1)=lbound(bounds%kb%ibxz_f,2)
  ab(4,5,1)=ubound(bounds%kb%ibxz_f,2)
  ab(5,5,1)=lbound(bounds%kb%ibxz_f,3)
  ab(6,5,1)=ubound(bounds%kb%ibxz_f,3)

  ab(1,6,1)=lbound(bounds%kb%ibxy_f,1)
  ab(2,6,1)=ubound(bounds%kb%ibxy_f,1)
  ab(3,6,1)=lbound(bounds%kb%ibxy_f,2)
  ab(4,6,1)=ubound(bounds%kb%ibxy_f,2)
  ab(5,6,1)=lbound(bounds%kb%ibxy_f,3)
  ab(6,6,1)=ubound(bounds%kb%ibxy_f,3)

  ab(1,1,2)=lbound(bounds%sb%ibzzx_c,1)
  ab(2,1,2)=ubound(bounds%sb%ibzzx_c,1)
  ab(3,1,2)=lbound(bounds%sb%ibzzx_c,2)
  ab(4,1,2)=ubound(bounds%sb%ibzzx_c,2)
  ab(5,1,2)=lbound(bounds%sb%ibzzx_c,3)
  ab(6,1,2)=ubound(bounds%sb%ibzzx_c,3)

  ab(1,2,2)=lbound(bounds%sb%ibyyzz_c,1)
  ab(2,2,2)=ubound(bounds%sb%ibyyzz_c,1)
  ab(3,2,2)=lbound(bounds%sb%ibyyzz_c,2)
  ab(4,2,2)=ubound(bounds%sb%ibyyzz_c,2)
  ab(5,2,2)=lbound(bounds%sb%ibyyzz_c,3)
  ab(6,2,2)=ubound(bounds%sb%ibyyzz_c,3)

  ab(1,3,2)=lbound(bounds%sb%ibxy_ff,1)
  ab(2,3,2)=ubound(bounds%sb%ibxy_ff,1)
  ab(3,3,2)=lbound(bounds%sb%ibxy_ff,2)
  ab(4,3,2)=ubound(bounds%sb%ibxy_ff,2)
  ab(5,3,2)=lbound(bounds%sb%ibxy_ff,3)
  ab(6,3,2)=ubound(bounds%sb%ibxy_ff,3)

  ab(1,4,2)=lbound(bounds%sb%ibzzx_f,1)
  ab(2,4,2)=ubound(bounds%sb%ibzzx_f,1)
  ab(3,4,2)=lbound(bounds%sb%ibzzx_f,2)
  ab(4,4,2)=ubound(bounds%sb%ibzzx_f,2)
  ab(5,4,2)=lbound(bounds%sb%ibzzx_f,3)
  ab(6,4,2)=ubound(bounds%sb%ibzzx_f,3)

  ab(1,5,2)=lbound(bounds%sb%ibyyzz_f,1)
  ab(2,5,2)=ubound(bounds%sb%ibyyzz_f,1)
  ab(3,5,2)=lbound(bounds%sb%ibyyzz_f,2)
  ab(4,5,2)=ubound(bounds%sb%ibyyzz_f,2)
  ab(5,5,2)=lbound(bounds%sb%ibyyzz_f,3)
  ab(6,5,2)=ubound(bounds%sb%ibyyzz_f,3)

  ab(1,1,3)=lbound(bounds%gb%ibzxx_c,1)
  ab(2,1,3)=ubound(bounds%gb%ibzxx_c,1)
  ab(3,1,3)=lbound(bounds%gb%ibzxx_c,2)
  ab(4,1,3)=ubound(bounds%gb%ibzxx_c,2)
  ab(5,1,3)=lbound(bounds%gb%ibzxx_c,3)
  ab(6,1,3)=ubound(bounds%gb%ibzxx_c,3)

  ab(1,2,3)=lbound(bounds%gb%ibxxyy_c,1)
  ab(2,2,3)=ubound(bounds%gb%ibxxyy_c,1)
  ab(3,2,3)=lbound(bounds%gb%ibxxyy_c,2)
  ab(4,2,3)=ubound(bounds%gb%ibxxyy_c,2)
  ab(5,2,3)=lbound(bounds%gb%ibxxyy_c,3)
  ab(6,2,3)=ubound(bounds%gb%ibxxyy_c,3)

  ab(1,3,3)=lbound(bounds%gb%ibyz_ff,1)
  ab(2,3,3)=ubound(bounds%gb%ibyz_ff,1)
  ab(3,3,3)=lbound(bounds%gb%ibyz_ff,2)
  ab(4,3,3)=ubound(bounds%gb%ibyz_ff,2)
  ab(5,3,3)=lbound(bounds%gb%ibyz_ff,3)
  ab(6,3,3)=ubound(bounds%gb%ibyz_ff,3)

  ab(1,4,3)=lbound(bounds%gb%ibzxx_f,1)
  ab(2,4,3)=ubound(bounds%gb%ibzxx_f,1)
  ab(3,4,3)=lbound(bounds%gb%ibzxx_f,2)
  ab(4,4,3)=ubound(bounds%gb%ibzxx_f,2)
  ab(5,4,3)=lbound(bounds%gb%ibzxx_f,3)
  ab(6,4,3)=ubound(bounds%gb%ibzxx_f,3)

  ab(1,5,3)=lbound(bounds%gb%ibxxyy_f,1)
  ab(2,5,3)=ubound(bounds%gb%ibxxyy_f,1)
  ab(3,5,3)=lbound(bounds%gb%ibxxyy_f,2)
  ab(4,5,3)=ubound(bounds%gb%ibxxyy_f,2)
  ab(5,5,3)=lbound(bounds%gb%ibxxyy_f,3)
  ab(6,5,3)=ubound(bounds%gb%ibxxyy_f,3)

  ab(1,1,4)=lbound(bounds%ibyyzz_r,1)
  ab(2,1,4)=ubound(bounds%ibyyzz_r,1)
  ab(3,1,4)=lbound(bounds%ibyyzz_r,2)
  ab(4,1,4)=ubound(bounds%ibyyzz_r,2)
  ab(5,1,4)=lbound(bounds%ibyyzz_r,3)
  ab(6,1,4)=ubound(bounds%ibyyzz_r,3)

  mpi_ncounts(1,1) = (ab(2,1,1)-ab(1,1,1)+1)*(ab(4,1,1)-ab(3,1,1)+1)*(ab(6,1,1)-ab(5,1,1)+1)
  mpi_ncounts(2,1) = (ab(2,2,1)-ab(1,2,1)+1)*(ab(4,2,1)-ab(3,2,1)+1)*(ab(6,2,1)-ab(5,2,1)+1)
  mpi_ncounts(3,1) = (ab(2,3,1)-ab(1,3,1)+1)*(ab(4,3,1)-ab(3,3,1)+1)*(ab(6,3,1)-ab(5,3,1)+1)
  mpi_ncounts(4,1) = (ab(2,4,1)-ab(1,4,1)+1)*(ab(4,4,1)-ab(3,4,1)+1)*(ab(6,4,1)-ab(5,4,1)+1)
  mpi_ncounts(5,1) = (ab(2,5,1)-ab(1,5,1)+1)*(ab(4,5,1)-ab(3,5,1)+1)*(ab(6,5,1)-ab(5,5,1)+1)
  mpi_ncounts(6,1) = (ab(2,6,1)-ab(1,6,1)+1)*(ab(4,6,1)-ab(3,6,1)+1)*(ab(6,6,1)-ab(5,6,1)+1)

  mpi_ncounts(1,2) = (ab(2,1,2)-ab(1,1,2)+1)*(ab(4,1,2)-ab(3,1,2)+1)*(ab(6,1,2)-ab(5,1,2)+1)
  mpi_ncounts(2,2) = (ab(2,2,2)-ab(1,2,2)+1)*(ab(4,2,2)-ab(3,2,2)+1)*(ab(6,2,2)-ab(5,2,2)+1)
  mpi_ncounts(3,2) = (ab(2,3,2)-ab(1,3,2)+1)*(ab(4,3,2)-ab(3,3,2)+1)*(ab(6,3,2)-ab(5,3,2)+1)
  mpi_ncounts(4,2) = (ab(2,4,2)-ab(1,4,2)+1)*(ab(4,4,2)-ab(3,4,2)+1)*(ab(6,4,2)-ab(5,4,2)+1)
  mpi_ncounts(5,2) = (ab(2,5,2)-ab(1,5,2)+1)*(ab(4,5,2)-ab(3,5,2)+1)*(ab(6,5,2)-ab(5,5,2)+1)

  mpi_ncounts(1,3) = (ab(2,1,3)-ab(1,1,3)+1)*(ab(4,1,3)-ab(3,1,3)+1)*(ab(6,1,3)-ab(5,1,3)+1)
  mpi_ncounts(2,3) = (ab(2,2,3)-ab(1,2,3)+1)*(ab(4,2,3)-ab(3,2,3)+1)*(ab(6,2,3)-ab(5,2,3)+1)
  mpi_ncounts(3,3) = (ab(2,3,3)-ab(1,3,3)+1)*(ab(4,3,3)-ab(3,3,3)+1)*(ab(6,3,3)-ab(5,3,3)+1)
  mpi_ncounts(4,3) = (ab(2,4,3)-ab(1,4,3)+1)*(ab(4,4,3)-ab(3,4,3)+1)*(ab(6,4,3)-ab(5,4,3)+1)
  mpi_ncounts(5,3) = (ab(2,5,3)-ab(1,5,3)+1)*(ab(4,5,3)-ab(3,5,3)+1)*(ab(6,5,3)-ab(5,5,3)+1)

  mpi_ncounts(1,4) = (ab(2,1,4)-ab(1,1,4)+1)*(ab(4,1,4)-ab(3,1,4)+1)*(ab(6,1,4)-ab(5,1,4)+1)

END SUBROUTINE get_convarrays_bounds




subroutine create_convolutions_windows_root(ab, mpi_ncounts, bounds, wins_bounds, wins_arrays)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(6,6,4),intent(in):: ab
  integer(kind=mpi_address_kind),dimension(6,4),intent(in):: mpi_ncounts
  type(convolutions_bounds),intent(in):: bounds
  integer,dimension(4),intent(out):: wins_bounds
  integer,dimension(6,4),intent(out):: wins_arrays

  ! Local variables
  integer(kind=mpi_address_kind):: size_of_integer, mpi_6, mpi_30, mpi_36
  integer:: ierr

  call mpi_type_size(mpi_integer, size_of_integer, ierr)
  mpi_6=6
  mpi_30=30
  mpi_36=36

  call mpi_win_create(ab(1,1,1), size_of_integer*mpi_36, size_of_integer, mpi_info_null, &
       mpi_comm_world, wins_bounds(1), ierr)
  call mpi_win_create(ab(1,1,2), size_of_integer*mpi_30, size_of_integer, mpi_info_null, &
       mpi_comm_world, wins_bounds(2), ierr)
  call mpi_win_create(ab(1,1,3), size_of_integer*mpi_30, size_of_integer, mpi_info_null, &
       mpi_comm_world, wins_bounds(3), ierr)
  call mpi_win_create(ab(1,1,4), size_of_integer*mpi_6, size_of_integer, mpi_info_null, &
       mpi_comm_world, wins_bounds(4), ierr)

  call mpi_win_create(bounds%kb%ibyz_c, size_of_integer*mpi_ncounts(1,1), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(1,1), ierr)
  call mpi_win_create(bounds%kb%ibxz_c, size_of_integer*mpi_ncounts(2,1), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(2,1), ierr)
  call mpi_win_create(bounds%kb%ibxy_c, size_of_integer*mpi_ncounts(3,1), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(3,1), ierr)
  call mpi_win_create(bounds%kb%ibyz_f, size_of_integer*mpi_ncounts(4,1), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(4,1), ierr)
  call mpi_win_create(bounds%kb%ibxz_f, size_of_integer*mpi_ncounts(5,1), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(5,1), ierr)
  call mpi_win_create(bounds%kb%ibxy_f, size_of_integer*mpi_ncounts(6,1), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(6,1), ierr)

  call mpi_win_create(bounds%sb%ibzzx_c, size_of_integer*mpi_ncounts(1,2), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(1,2), ierr)
  call mpi_win_create(bounds%sb%ibyyzz_c, size_of_integer*mpi_ncounts(2,2), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(2,2), ierr)
  call mpi_win_create(bounds%sb%ibxy_ff, size_of_integer*mpi_ncounts(3,2), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(3,2), ierr)
  call mpi_win_create(bounds%sb%ibzzx_f, size_of_integer*mpi_ncounts(4,2), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(4,2), ierr)
  call mpi_win_create(bounds%sb%ibyyzz_f, size_of_integer*mpi_ncounts(5,2), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(5,2), ierr)

  call mpi_win_create(bounds%gb%ibzxx_c, size_of_integer*mpi_ncounts(1,3), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(1,3), ierr)
  call mpi_win_create(bounds%gb%ibxxyy_c, size_of_integer*mpi_ncounts(2,3), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(2,3), ierr)
  call mpi_win_create(bounds%gb%ibyz_ff, size_of_integer*mpi_ncounts(3,3), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(3,3), ierr)
  call mpi_win_create(bounds%gb%ibzxx_f, size_of_integer*mpi_ncounts(4,3), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(4,3), ierr)
  call mpi_win_create(bounds%gb%ibxxyy_f, size_of_integer*mpi_ncounts(5,3), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(5,3), ierr)

  call mpi_win_create(bounds%ibyyzz_r, size_of_integer*mpi_ncounts(1,4), size_of_integer, &
       mpi_info_null, mpi_comm_world, wins_arrays(1,4), ierr)

END SUBROUTINE create_convolutions_windows_root



subroutine create_convolutions_windows_else(ab, ifake, wins_bounds, wins_arrays)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(6,6,4),intent(in):: ab
  integer,dimension(6,4),intent(in):: ifake
  integer,dimension(4),intent(out):: wins_bounds
  integer,dimension(6,4),intent(out):: wins_arrays

  ! Local variables
  integer(kind=mpi_address_kind):: mpi_0
  integer:: ierr

  mpi_0=0

  call mpi_win_create(ab(1,1,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_bounds(1), ierr)
  call mpi_win_create(ab(1,1,2), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_bounds(2), ierr)
  call mpi_win_create(ab(1,1,3), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_bounds(3), ierr)
  call mpi_win_create(ab(1,1,4), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_bounds(4), ierr)

  call mpi_win_create(ifake(1,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(1,1), ierr)
  call mpi_win_create(ifake(2,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(2,1), ierr)
  call mpi_win_create(ifake(3,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(3,1), ierr)
  call mpi_win_create(ifake(4,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(4,1), ierr)
  call mpi_win_create(ifake(5,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(5,1), ierr)
  call mpi_win_create(ifake(6,1), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(6,1), ierr)

  call mpi_win_create(ifake(1,2), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(1,2), ierr)
  call mpi_win_create(ifake(2,2), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(2,2), ierr)
  call mpi_win_create(ifake(3,2), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(3,2), ierr)
  call mpi_win_create(ifake(4,2), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(4,2), ierr)
  call mpi_win_create(ifake(5,2), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(5,2), ierr)

  call mpi_win_create(ifake(1,3), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(1,3), ierr)
  call mpi_win_create(ifake(2,3), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(2,3), ierr)
  call mpi_win_create(ifake(3,3), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(3,3), ierr)
  call mpi_win_create(ifake(4,3), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(4,3), ierr)
  call mpi_win_create(ifake(5,3), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(5,3), ierr)

  call mpi_win_create(ifake(1,4), mpi_0, 1, mpi_info_null, mpi_comm_world, wins_arrays(1,4), ierr)

END SUBROUTINE create_convolutions_windows_else


subroutine convolutions_bounds_fences(wins_bounds)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(4),intent(inout):: wins_bounds

  ! Local variables
  integer:: ierr

  call mpi_win_fence(0, wins_bounds(1), ierr)
  call mpi_win_fence(0, wins_bounds(2), ierr)
  call mpi_win_fence(0, wins_bounds(3), ierr)
  call mpi_win_fence(0, wins_bounds(4), ierr)

END SUBROUTINE convolutions_bounds_fences


subroutine convolutions_arrays_fences(wins_arrays)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(6,4),intent(inout):: wins_arrays

  ! Local variables
  integer:: ierr

  call mpi_win_fence(0, wins_arrays(1,1), ierr)
  call mpi_win_fence(0, wins_arrays(2,1), ierr)
  call mpi_win_fence(0, wins_arrays(3,1), ierr)
  call mpi_win_fence(0, wins_arrays(4,1), ierr)
  call mpi_win_fence(0, wins_arrays(5,1), ierr)
  call mpi_win_fence(0, wins_arrays(6,1), ierr)

  call mpi_win_fence(0, wins_arrays(1,2), ierr)
  call mpi_win_fence(0, wins_arrays(2,2), ierr)
  call mpi_win_fence(0, wins_arrays(3,2), ierr)
  call mpi_win_fence(0, wins_arrays(4,2), ierr)
  call mpi_win_fence(0, wins_arrays(5,2), ierr)

  call mpi_win_fence(0, wins_arrays(1,3), ierr)
  call mpi_win_fence(0, wins_arrays(2,3), ierr)
  call mpi_win_fence(0, wins_arrays(3,3), ierr)
  call mpi_win_fence(0, wins_arrays(4,3), ierr)
  call mpi_win_fence(0, wins_arrays(5,3), ierr)

  call mpi_win_fence(0, wins_arrays(1,4), ierr)

END SUBROUTINE convolutions_arrays_fences


subroutine allocate_convolutions_bounds(ab, subname, bounds)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(6,6,4),intent(in):: ab
  character(len=*),intent(in):: subname
  type(convolutions_bounds),intent(out):: bounds

  ! Local variables
  integer:: istat

  allocate(bounds%kb%ibyz_c(ab(1,1,1):ab(2,1,1),ab(3,1,1):ab(4,1,1),ab(5,1,1):ab(6,1,1)), stat=istat)
  call memocc(istat, bounds%kb%ibyz_c, 'bounds%kb%ibyz_c', subname)
  allocate(bounds%kb%ibxz_c(ab(1,2,1):ab(2,2,1),ab(3,2,1):ab(4,2,1),ab(5,2,1):ab(6,2,1)), stat=istat)
  call memocc(istat, bounds%kb%ibxz_c, 'bounds%kb%ibxz_c', subname)
  allocate(bounds%kb%ibxy_c(ab(1,3,1):ab(2,3,1),ab(3,3,1):ab(4,3,1),ab(5,3,1):ab(6,3,1)), stat=istat)
  call memocc(istat, bounds%kb%ibxy_c, 'bounds%kb%ibxy_c', subname)
  allocate(bounds%kb%ibyz_f(ab(1,4,1):ab(2,4,1),ab(3,4,1):ab(4,4,1),ab(5,4,1):ab(6,4,1)), stat=istat)
  call memocc(istat, bounds%kb%ibyz_f, 'bounds%kb%ibyz_f', subname)
  allocate(bounds%kb%ibxz_f(ab(1,5,1):ab(2,5,1),ab(3,5,1):ab(4,5,1),ab(5,5,1):ab(6,5,1)), stat=istat)
  call memocc(istat, bounds%kb%ibxz_f, 'bounds%kb%ibxz_f', subname)
  allocate(bounds%kb%ibxy_f(ab(1,6,1):ab(2,6,1),ab(3,6,1):ab(4,6,1),ab(5,6,1):ab(6,6,1)), stat=istat)
  call memocc(istat, bounds%kb%ibxy_f, 'bounds%kb%ibxy_f', subname)

  allocate(bounds%sb%ibzzx_c(ab(1,1,2):ab(2,1,2),ab(3,1,2):ab(4,1,2),ab(5,1,2):ab(6,1,2)), stat=istat)
  call memocc(istat, bounds%sb%ibzzx_c, 'bounds%sb%ibzzx_c', subname)
  allocate(bounds%sb%ibyyzz_c(ab(1,2,2):ab(2,2,2),ab(3,2,2):ab(4,2,2),ab(5,2,2):ab(6,2,2)), stat=istat)
  call memocc(istat, bounds%sb%ibyyzz_c, 'bounds%sb%ibyyzz_c', subname)
  allocate(bounds%sb%ibxy_ff(ab(1,3,2):ab(2,3,2),ab(3,3,2):ab(4,3,2),ab(5,3,2):ab(6,3,2)), stat=istat)
  call memocc(istat, bounds%sb%ibxy_ff, 'bounds%sb%ibxy_ff', subname)
  allocate(bounds%sb%ibzzx_f(ab(1,4,2):ab(2,4,2),ab(3,4,2):ab(4,4,2),ab(5,4,2):ab(6,4,2)), stat=istat)
  call memocc(istat, bounds%sb%ibzzx_f, 'bounds%sb%ibzzx_f', subname)
  allocate(bounds%sb%ibyyzz_f(ab(1,5,2):ab(2,5,2),ab(3,5,2):ab(4,5,2),ab(5,5,2):ab(6,5,2)), stat=istat)
  call memocc(istat, bounds%sb%ibyyzz_f, 'bounds%sb%ibyyzz_f', subname)

  allocate(bounds%gb%ibzxx_c(ab(1,1,3):ab(2,1,3),ab(3,1,3):ab(4,1,3),ab(5,1,3):ab(6,1,3)), stat=istat)
  call memocc(istat, bounds%gb%ibzxx_c, 'bounds%gb%ibzxx_c', subname)
  allocate(bounds%gb%ibxxyy_c(ab(1,2,3):ab(2,2,3),ab(3,2,3):ab(4,2,3),ab(5,2,3):ab(6,2,3)), stat=istat)
  call memocc(istat, bounds%gb%ibxxyy_c, 'bounds%gb%ibxxyy_c', subname)
  allocate(bounds%gb%ibyz_ff(ab(1,3,3):ab(2,3,3),ab(3,3,3):ab(4,3,3),ab(5,3,3):ab(6,3,3)), stat=istat)
  call memocc(istat, bounds%gb%ibyz_ff, 'bounds%gb%ibyz_ff', subname)
  allocate(bounds%gb%ibzxx_f(ab(1,4,3):ab(2,4,3),ab(3,4,3):ab(4,4,3),ab(5,4,3):ab(6,4,3)), stat=istat)
  call memocc(istat, bounds%gb%ibzxx_f, 'bounds%gb%ibzxx_f', subname)
  allocate(bounds%gb%ibxxyy_f(ab(1,5,3):ab(2,5,3),ab(3,5,3):ab(4,5,3),ab(5,5,3):ab(6,5,3)), stat=istat)
  call memocc(istat, bounds%gb%ibxxyy_f, 'bounds%gb%ibxxyy_f', subname)

  allocate(bounds%ibyyzz_r(ab(1,1,4):ab(2,1,4),ab(3,1,4):ab(4,1,4),ab(5,1,4):ab(6,1,4)), stat=istat)
  call memocc(istat, bounds%ibyyzz_r, 'bounds%ibyyzz_r', subname)

END SUBROUTINE allocate_convolutions_bounds


subroutine get_convolutions_bounds(root, ab, wins_bounds)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: root
  integer,dimension(6,6,4),intent(inout):: ab
  integer,dimension(4),intent(inout):: wins_bounds
  integer(kind=mpi_address_kind):: mpi_0

  ! Local variables
  integer:: ierr

  mpi_0=0

  call mpi_get(ab(1,1,1), 36, mpi_integer, root, mpi_0, 36, mpi_integer, wins_bounds(1), ierr)
  call mpi_get(ab(1,1,2), 30, mpi_integer, root, mpi_0, 30, mpi_integer, wins_bounds(2), ierr)
  call mpi_get(ab(1,1,3), 30, mpi_integer, root, mpi_0, 30, mpi_integer, wins_bounds(3), ierr)
  call mpi_get(ab(1,1,4), 6, mpi_integer, root, mpi_0, 6, mpi_integer, wins_bounds(4), ierr)

END SUBROUTINE get_convolutions_bounds



subroutine get_convolutions_arrays(root, ab, bounds, wins_arrays)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: root
  integer,dimension(6,6,4),intent(in):: ab
  type(convolutions_bounds),intent(inout):: bounds
  integer,dimension(6,4),intent(inout):: wins_arrays

  ! Local variables
  integer,dimension(6,4):: ncounts
  integer:: ierr
  integer(kind=mpi_address_kind):: mpi_0


  ncounts(1,1) = (ab(2,1,1)-ab(1,1,1)+1)*(ab(4,1,1)-ab(3,1,1)+1)*(ab(6,1,1)-ab(5,1,1)+1)
  ncounts(2,1) = (ab(2,2,1)-ab(1,2,1)+1)*(ab(4,2,1)-ab(3,2,1)+1)*(ab(6,2,1)-ab(5,2,1)+1)
  ncounts(3,1) = (ab(2,3,1)-ab(1,3,1)+1)*(ab(4,3,1)-ab(3,3,1)+1)*(ab(6,3,1)-ab(5,3,1)+1)
  ncounts(4,1) = (ab(2,4,1)-ab(1,4,1)+1)*(ab(4,4,1)-ab(3,4,1)+1)*(ab(6,4,1)-ab(5,4,1)+1)
  ncounts(5,1) = (ab(2,5,1)-ab(1,5,1)+1)*(ab(4,5,1)-ab(3,5,1)+1)*(ab(6,5,1)-ab(5,5,1)+1)
  ncounts(6,1) = (ab(2,6,1)-ab(1,6,1)+1)*(ab(4,6,1)-ab(3,6,1)+1)*(ab(6,6,1)-ab(5,6,1)+1)

  ncounts(1,2) = (ab(2,1,2)-ab(1,1,2)+1)*(ab(4,1,2)-ab(3,1,2)+1)*(ab(6,1,2)-ab(5,1,2)+1)
  ncounts(2,2) = (ab(2,2,2)-ab(1,2,2)+1)*(ab(4,2,2)-ab(3,2,2)+1)*(ab(6,2,2)-ab(5,2,2)+1)
  ncounts(3,2) = (ab(2,3,2)-ab(1,3,2)+1)*(ab(4,3,2)-ab(3,3,2)+1)*(ab(6,3,2)-ab(5,3,2)+1)
  ncounts(4,2) = (ab(2,4,2)-ab(1,4,2)+1)*(ab(4,4,2)-ab(3,4,2)+1)*(ab(6,4,2)-ab(5,4,2)+1)
  ncounts(5,2) = (ab(2,5,2)-ab(1,5,2)+1)*(ab(4,5,2)-ab(3,5,2)+1)*(ab(6,5,2)-ab(5,5,2)+1)

  ncounts(1,3) = (ab(2,1,3)-ab(1,1,3)+1)*(ab(4,1,3)-ab(3,1,3)+1)*(ab(6,1,3)-ab(5,1,3)+1)
  ncounts(2,3) = (ab(2,2,3)-ab(1,2,3)+1)*(ab(4,2,3)-ab(3,2,3)+1)*(ab(6,2,3)-ab(5,2,3)+1)
  ncounts(3,3) = (ab(2,3,3)-ab(1,3,3)+1)*(ab(4,3,3)-ab(3,3,3)+1)*(ab(6,3,3)-ab(5,3,3)+1)
  ncounts(4,3) = (ab(2,4,3)-ab(1,4,3)+1)*(ab(4,4,3)-ab(3,4,3)+1)*(ab(6,4,3)-ab(5,4,3)+1)
  ncounts(5,3) = (ab(2,5,3)-ab(1,5,3)+1)*(ab(4,5,3)-ab(3,5,3)+1)*(ab(6,5,3)-ab(5,5,3)+1)

  ncounts(1,4) = (ab(2,1,4)-ab(1,1,4)+1)*(ab(4,1,4)-ab(3,1,4)+1)*(ab(6,1,4)-ab(5,1,4)+1)

  mpi_0=0

  call mpi_get(bounds%kb%ibyz_c, ncounts(1,1), mpi_integer, root, mpi_0, ncounts(1,1), mpi_integer, wins_arrays(1,1), ierr)
  call mpi_get(bounds%kb%ibxz_c, ncounts(2,1), mpi_integer, root, mpi_0, ncounts(2,1), mpi_integer, wins_arrays(2,1), ierr)
  call mpi_get(bounds%kb%ibxy_c, ncounts(3,1), mpi_integer, root, mpi_0, ncounts(3,1), mpi_integer, wins_arrays(3,1), ierr)
  call mpi_get(bounds%kb%ibyz_f, ncounts(4,1), mpi_integer, root, mpi_0, ncounts(4,1), mpi_integer, wins_arrays(4,1), ierr)
  call mpi_get(bounds%kb%ibxz_f, ncounts(5,1), mpi_integer, root, mpi_0, ncounts(5,1), mpi_integer, wins_arrays(5,1), ierr)
  call mpi_get(bounds%kb%ibxy_f, ncounts(6,1), mpi_integer, root, mpi_0, ncounts(6,1), mpi_integer, wins_arrays(6,1), ierr)

  call mpi_get(bounds%sb%ibzzx_c, ncounts(1,2), mpi_integer, root, mpi_0, ncounts(1,2), mpi_integer, wins_arrays(1,2), ierr)
  call mpi_get(bounds%sb%ibyyzz_c, ncounts(2,2), mpi_integer, root, mpi_0, ncounts(2,2), mpi_integer, wins_arrays(2,2), ierr)
  call mpi_get(bounds%sb%ibxy_ff, ncounts(3,2), mpi_integer, root, mpi_0, ncounts(3,2), mpi_integer, wins_arrays(3,2), ierr)
  call mpi_get(bounds%sb%ibzzx_f, ncounts(4,2), mpi_integer, root, mpi_0, ncounts(4,2), mpi_integer, wins_arrays(4,2), ierr)
  call mpi_get(bounds%sb%ibyyzz_f, ncounts(5,2), mpi_integer, root, mpi_0, ncounts(5,2), mpi_integer, wins_arrays(5,2), ierr)

  call mpi_get(bounds%gb%ibzxx_c, ncounts(1,3), mpi_integer, root, mpi_0, ncounts(1,3), mpi_integer, wins_arrays(1,3), ierr)
  call mpi_get(bounds%gb%ibxxyy_c, ncounts(2,3), mpi_integer, root, mpi_0, ncounts(2,3), mpi_integer, wins_arrays(2,3), ierr)
  call mpi_get(bounds%gb%ibyz_ff, ncounts(3,3), mpi_integer, root, mpi_0, ncounts(3,3), mpi_integer, wins_arrays(3,3), ierr)
  call mpi_get(bounds%gb%ibzxx_f, ncounts(4,3), mpi_integer, root, mpi_0, ncounts(4,3), mpi_integer, wins_arrays(4,3), ierr)
  call mpi_get(bounds%gb%ibxxyy_f, ncounts(5,3), mpi_integer, root, mpi_0, ncounts(5,3), mpi_integer, wins_arrays(5,3), ierr)

  call mpi_get(bounds%ibyyzz_r, ncounts(1,4), mpi_integer, root, mpi_0, ncounts(1,4), mpi_integer, wins_arrays(1,4), ierr)


end subroutine get_convolutions_arrays





subroutine free_convolutions_bounds_windows(wins_bounds)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(4),intent(inout):: wins_bounds

  ! Local variables
  integer:: ierr

  call mpi_win_free(wins_bounds(1), ierr)
  call mpi_win_free(wins_bounds(2), ierr)
  call mpi_win_free(wins_bounds(3), ierr)
  call mpi_win_free(wins_bounds(4), ierr)

END SUBROUTINE 


subroutine free_convolutions_arrays_windows(wins_arrays)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,dimension(6,4),intent(inout):: wins_arrays

  ! Local variables
  integer:: ierr

  call mpi_win_free(wins_arrays(1,1), ierr)
  call mpi_win_free(wins_arrays(2,1), ierr)
  call mpi_win_free(wins_arrays(3,1), ierr)
  call mpi_win_free(wins_arrays(4,1), ierr)
  call mpi_win_free(wins_arrays(5,1), ierr)
  call mpi_win_free(wins_arrays(6,1), ierr)

  call mpi_win_free(wins_arrays(1,2), ierr)
  call mpi_win_free(wins_arrays(2,2), ierr)
  call mpi_win_free(wins_arrays(3,2), ierr)
  call mpi_win_free(wins_arrays(4,2), ierr)
  call mpi_win_free(wins_arrays(5,2), ierr)

  call mpi_win_free(wins_arrays(1,3), ierr)
  call mpi_win_free(wins_arrays(2,3), ierr)
  call mpi_win_free(wins_arrays(3,3), ierr)
  call mpi_win_free(wins_arrays(4,3), ierr)
  call mpi_win_free(wins_arrays(5,3), ierr)

  call mpi_win_free(wins_arrays(1,4), ierr)

END SUBROUTINE free_convolutions_arrays_windows
