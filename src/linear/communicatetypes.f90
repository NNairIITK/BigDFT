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
   integer:: addr_geocode, addr_hybrid_on, addr_ns1, addr_ns2, addr_ns3, addr_nsi1, addr_nsi2
   integer:: addr_nsi3, addr_localnorb, addr_outofzone, addr_locregCenter, addr_locrad, addr_llr
   integer,dimension(12):: blocklengths, dspls, types

   
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

   call mpi_type_struct(ncount, blocklengths, dspls, types, mpi_tmptype, ierr)
   !call mpi_type_struct(1, 1, 0, mpi_character, mpi_tmptype, ierr)
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
   integer:: ierr, ncount, commtype, addr_wfd, addr_nvctr_c, addr_nvctr_f, addr_nseg_c, addr_nseg_f
   integer:: addr_keyglob, addr_keygloc, addr_keyvloc, addr_keyvglob
   character(len=*),parameter:: subname='communicate_wavefunctions_descriptors'
   integer,dimension(4):: blocklengths, dspls, types

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
   
   call mpi_type_struct(ncount, blocklengths, dspls, types, commtype, ierr)
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

   call mpi_type_struct(ncount, blocklengths, dspls, types, commtype, ierr)
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
   integer,dimension(ncount):: types, blocklengths, dspls
   integer:: addr_kb, addr_ibyz_c, addr_ibxz_c, addr_ibxy_c, addr_ibyz_f, addr_ibxz_f, addr_ibxy_f

   call mpi_get_address(kb, addr_kb, ierr)

   call getbounds(iproc, root, kb%ibyz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibyz_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibyz_c, 'kb%ibyz_c', subname)
   end if
   blocklengths(1) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(kb%ibyz_c, addr_ibyz_c, ierr)
   dspls(1) = addr_ibyz_c - addr_kb
   types(1) = mpi_integer

   call getbounds(iproc, root, kb%ibxz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxz_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxz_c, 'kb%ibxz_c', subname)
   end if
   blocklengths(2) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(kb%ibxz_c, addr_ibxz_c, ierr)
   dspls(2) = addr_ibxz_c - addr_kb
   types(2) = mpi_integer


   call getbounds(iproc, root, kb%ibxy_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxy_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxy_c, 'kb%ibxy_c', subname)
   end if
   blocklengths(3) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(kb%ibxy_c, addr_ibxy_c, ierr)
   dspls(3) = addr_ibxy_c - addr_kb
   types(3) = mpi_integer


   call getbounds(iproc, root, kb%ibyz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibyz_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibyz_f, 'kb%ibyz_f', subname)
   end if
   blocklengths(4) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(kb%ibyz_f, addr_ibyz_f, ierr)
   dspls(4) = addr_ibyz_f - addr_kb
   types(4) = mpi_integer

   call getbounds(iproc, root, kb%ibxz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxz_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxz_f, 'kb%ibxz_f', subname)
   end if
   blocklengths(5) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(kb%ibxz_f, addr_ibxz_f, ierr)
   dspls(5) = addr_ibxz_f - addr_kb
   types(5) = mpi_integer

   call getbounds(iproc, root, kb%ibxy_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(kb%ibxy_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, kb%ibxy_f, 'kb%ibxy_f', subname)
   end if
   blocklengths(6) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(kb%ibxy_f, addr_ibxy_f, ierr)
   dspls(6) = addr_ibxy_f - addr_kb
   types(6) = mpi_integer

   call mpi_type_struct(ncount, blocklengths, dspls, types, commtype, ierr)
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
   integer:: ierr, istat, is1, ie1, is2, ie2, is3, ie3
   character(len=*),parameter:: subname='communicate_shrink_bounds'
   integer,parameter:: ncount=5
   integer,dimension(ncount):: types, dspls, blocklengths
   integer:: addr_sb, addr_ibzzx_c, addr_ibyyzz_c, addr_ibxy_ff, addr_ibzzx_f, addr_ibyyzz_f, commtype


   call mpi_get_address(sb, addr_sb, ierr)

   call getbounds(iproc, root, sb%ibzzx_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibzzx_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibzzx_c, 'sb%ibzzx_c', subname)
   end if
   blocklengths(1) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(sb%ibzzx_c, addr_ibzzx_c, ierr)
   dspls(1) = addr_ibzzx_c - addr_sb
   types(1) = mpi_integer
   !call mpi_bcast(sb%ibzzx_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibyyzz_c, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibyyzz_c(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibyyzz_c, 'sb%ibyyzz_c', subname)
   end if
   blocklengths(2) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(sb%ibyyzz_c, addr_ibyyzz_c)
   dspls(2) = addr_ibyyzz_c - addr_sb
   types(2) = mpi_integer
   !call mpi_bcast(sb%ibyyzz_c, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibxy_ff, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibxy_ff(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibxy_ff, 'sb%ibxy_ff', subname)
   end if
   blocklengths(3) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(sb%ibxy_ff, addr_ibxy_ff, ierr)
   dspls(3) = addr_ibxy_ff - addr_sb
   types(3) = mpi_integer
   !call mpi_bcast(sb%ibxy_ff, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibzzx_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibzzx_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibzzx_f, 'sb%ibzzx_f', subname)
   end if
   blocklengths(4) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(sb%ibzzx_f, addr_ibzzx_f, ierr)
   dspls(4) = addr_ibzzx_f - addr_sb
   types(4) = mpi_integer
   !call mpi_bcast(sb%ibzzx_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call getbounds(iproc, root, sb%ibyyzz_f, is1, ie1, is2, ie2, is3, ie3)
   if(iproc/=root) then
      allocate(sb%ibyyzz_f(is1:ie1,is2:ie2,is3:ie3), stat=istat)
      call memocc(istat, sb%ibyyzz_f, 'sb%ibyyzz_f', subname)
   end if
   blocklengths(5) = (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1)
   call mpi_get_address(sb%ibyyzz_f, addr_ibyyzz_f, ierr)
   dspls(5) = addr_ibyyzz_f - addr_sb
   types(5) = mpi_integer
   !call mpi_bcast(sb%ibyyzz_f, (ie1-is1+1)*(ie2-is2+1)*(ie3-is3+1), mpi_integer, root, mpi_comm_world, ierr)

   call mpi_type_struct(ncount, blocklengths, dspls, types, commtype, ierr)
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
