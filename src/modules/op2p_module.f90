!> @file
!!  File to define information used for the overlap point to point between wavefunctions
!! @author
!!    Copyright (C) 2011-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Modules which contains the defintions for overlap point to point
module overlap_point_to_point
   use module_base
   implicit none

   !By default variables are internal to the module
   private
   public :: initialize_OP2P_descriptors,OP2P_communication,OP2P_descriptors,free_OP2P_descriptors

   type OP2P_descriptors
      logical :: forsymop !< descriptor for symmetric operation
      integer :: ngroup,nsteps_max,ngroupp_max,ncomponents_max,ncomponents_results_max,iprocref
      integer, dimension(:), pointer :: ngroupp         !<= number of groups which belongs to each processor
      integer, dimension(:), pointer :: nprocgr         !<= number of processors which belongs to each group
      integer, dimension(:,:), pointer :: igrpr         !<= groups which are associated to each processor
      integer, dimension(:,:,:), pointer :: nvctr_par   !<= number of elements per group, for objects and results
      integer, dimension(:,:,:), pointer :: ioffp_group !<= global and local offsets of each group on the processor
      integer, dimension(:,:,:,:), pointer :: iprocpm1  !<= ascending and descending order for processors in the same group
      integer, dimension(:,:,:,:), pointer :: communication_schedule !<= processes to send and receive at each step
   end type OP2P_descriptors

   contains

   !<initialise the data needed to manage communication operations
   !objects_attributes(:,1) <= group to which the object belongs
   !objects_attributes(:,2) <= size in number of elements of the object
   !objects_attributes(:,3) <= size in number of elements of the results of the operations associated to the object
   !objects_repartition(:) <= distribution of the objects among processors
   subroutine initialize_OP2P_descriptors(symop,iproc,nproc,nobjects,objects_attributes,objects_repartition,OP2P)
      implicit none
      logical, intent(in) :: symop
      integer, intent(in) :: iproc,nproc,nobjects
      integer, dimension(0:nproc-1), intent(in) :: objects_repartition
      integer, dimension(nobjects,3), intent(in) :: objects_attributes
      type(OP2P_descriptors), intent(out) :: OP2P
      !local variables
      character(len=*), parameter :: subname='initialize_OP2P_communications'
      integer :: i_stat,ndebug=0  !< to be removed whenever with module_base
      integer :: igroup,iobject,jproc,ioffset_local,ioffset_global,iobj,nsize,igroup_previous
      integer :: nprocgr_max,iprocgrs,iprocgrr,icount,istep,kproc,nstepsm1
      integer :: iaddress_local,iaddress_results_local,nsize_results

      OP2P%forsymop=symop !< Assign the symmetry constraint to the structure

      !calculate total number of groups
      OP2P%ngroup=0
      do iobject=1,nobjects
         OP2P%ngroup=max(OP2P%ngroup,objects_attributes(iobject,1))
      end do
      !verify that all the groups exists
      parse_groups: do igroup=1,OP2P%ngroup
         do iobject=1,nobjects
            if (objects_attributes(iobject,1)==igroup) cycle parse_groups
         end do
         write(*,*)'ERROR: group id',igroup,' not found in the attributes'
         stop
      end do parse_groups

      !allocate arrays which depends of the number of groups
      allocate(OP2P%nvctr_par(0:nproc-1,OP2P%ngroup,2+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%nvctr_par,'OP2P%nvctr_par',subname)

      allocate(OP2P%ioffp_group(5,0:nproc-1,OP2P%ngroup+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%ioffp_group,'OP2P%ioffp_group',subname)

      !fill these arrays
      !initialise the arrays
      do igroup=1,OP2P%ngroup
         do jproc=0,nproc-1
            OP2P%nvctr_par(jproc,igroup,1)=0 !<number of elements for the objects
            OP2P%nvctr_par(jproc,igroup,2)=0 !<number of elements for the results
            OP2P%ioffp_group(1,jproc,igroup)=UNINITIALIZED(1) !< global object offset !fake value for cross-check initialisation
            OP2P%ioffp_group(2,jproc,igroup)=UNINITIALIZED(1) !< local object offset fake value for cross-check initialisation
            OP2P%ioffp_group(3,jproc,igroup)=1  !< local address position
            OP2P%ioffp_group(4,jproc,igroup)=0  !< local number of objects starting from local object
            OP2P%ioffp_group(5,jproc,igroup)=1  !< local address position of the result
         end do
      end do

      iobject=0
      ioffset_global=0 !< object to which a given couple group-processor is associated
      do jproc=0,nproc-1
         igroup_previous=UNINITIALIZED(1)
         ioffset_local=0
         iaddress_local=1
         iaddress_results_local=1
         do iobj=1,objects_repartition(jproc)
            iobject=iobject+1

            igroup=objects_attributes(iobject,1)
            nsize=objects_attributes(iobject,2)
            nsize_results=objects_attributes(iobject,3)
            !count the number of elements for each group and processor
            OP2P%nvctr_par(jproc,igroup,1)=OP2P%nvctr_par(jproc,igroup,1)+nsize
            OP2P%nvctr_par(jproc,igroup,2)=OP2P%nvctr_par(jproc,igroup,2)+nsize_results
            !any time that the group considered changes update the array of the offsets
            if (igroup_previous /= igroup) then
               !the array ioffp_group can be filled only once.
               !if it is not the case this means that the same group is 
               !not contiguous in a given processor. This is forbidden by the parallelisation scheme
               if (OP2P%ioffp_group(1,jproc,igroup) /= UNINITIALIZED(1)) then
                  write(*,*)'ERROR (global): the group',igroup,' is non contiguously distributed in processor ',jproc 
                  stop
               else
                  OP2P%ioffp_group(1,jproc,igroup)=ioffset_global 
               end if
               if (OP2P%ioffp_group(2,jproc,igroup) /= UNINITIALIZED(1)) then
                  write(*,*)'ERROR (local): the group',igroup,' is non contiguously distributed in processor ',jproc 
                  stop
               else
                  OP2P%ioffp_group(2,jproc,igroup)=ioffset_local
               end if
               OP2P%ioffp_group(3,jproc,igroup)=iaddress_local
               OP2P%ioffp_group(4,jproc,igroup)=1
               OP2P%ioffp_group(5,jproc,igroup)=iaddress_results_local
               igroup_previous=igroup
            else
               OP2P%ioffp_group(4,jproc,igroup)=OP2P%ioffp_group(4,jproc,igroup)+1
            end if
            !the offset is expressed in number of objects
            ioffset_global=ioffset_global+1
            ioffset_local=ioffset_local+1
            iaddress_local=iaddress_local+nsize
            iaddress_results_local=iaddress_results_local+nsize_results
         end do
      end do
      !check
      if (iobject /= nobjects) stop 'ERROR, objects_repartition'

      !clean the information for the offset, for the first part of the array
      do igroup=1,OP2P%ngroup
         do jproc=0,nproc-1
            if (OP2P%ioffp_group(1,jproc,igroup)==UNINITIALIZED(1)) OP2P%ioffp_group(1,jproc,igroup)=0
            if (OP2P%ioffp_group(2,jproc,igroup)==UNINITIALIZED(1)) OP2P%ioffp_group(2,jproc,igroup)=0
         end do
      end do

      allocate(OP2P%ngroupp(0:nproc-1+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%ngroupp,'OP2P%ngroupp',subname)

      !calculate the number of groups which belong to each processor
      !here we can allocate the working arrays giving the maximum
      !between the components for each group
      OP2P%ngroupp_max=0
      OP2P%ncomponents_max=0
      OP2P%ncomponents_results_max=0
      OP2P%iprocref=0 !processor which has the maximum number of groups
      do jproc=0,nproc-1
         OP2P%ngroupp(jproc)=0
         do igroup=1,OP2P%ngroup
            OP2P%ncomponents_max=max(OP2P%ncomponents_max,OP2P%nvctr_par(jproc,igroup,1))
            OP2P%ncomponents_results_max=max(OP2P%ncomponents_results_max,OP2P%nvctr_par(jproc,igroup,2))
            if (OP2P%nvctr_par(jproc,igroup,1) > 0) then !should be identical by considering objects or results
               OP2P%ngroupp(jproc)=OP2P%ngroupp(jproc)+1
            end if
         end do
         if (OP2P%ngroupp_max < OP2P%ngroupp(jproc)) OP2P%iprocref=jproc
         OP2P%ngroupp_max=max(OP2P%ngroupp_max,OP2P%ngroupp(jproc))
      end do

      !print *,'ngroupp_max,nproc,ngroupp',ngroupp_max,nproc,ngroupp

      allocate(OP2P%nprocgr(OP2P%ngroup+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%nprocgr,'OP2P%nprocgr',subname)

      !calculate the number of processors which belong to each group
      !here we can allocate the working arrays giving the maximum
      !between the components for each group, both for objects and results
      nprocgr_max=0
      do igroup=1,OP2P%ngroup
         OP2P%nprocgr(igroup)=0
         do jproc=0,nproc-1
            if (OP2P%nvctr_par(jproc,igroup,1) > 0) then !should be identical by considering objects or results
               OP2P%nprocgr(igroup)=OP2P%nprocgr(igroup)+1
            end if
         end do
         nprocgr_max=max(nprocgr_max,OP2P%nprocgr(igroup))
      end do
      if (OP2P%forsymop) then
         OP2P%nsteps_max=(nprocgr_max-1)/2+1 !here the number of steps should be changed for non-symmetric operation
      else
         OP2P%nsteps_max=(nprocgr_max-1)
      end if

      !determine the array of the groups which are of interest for this processor
      allocate(OP2P%igrpr(OP2P%ngroupp_max,0:nproc-1+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%igrpr,'OP2P%igrpr',subname)
      !processors lying above and below iproc in the list of communicating process
      allocate(OP2P%iprocpm1(2,0:nprocgr_max,OP2P%ngroupp_max,0:nproc-1+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%iprocpm1,'OP2P%iprocpm1',subname)

      !determine for each processor the groups which has to be used
      OP2P%iprocpm1(:,:,:,:)=UNINITIALIZED(1)
      do jproc=0,nproc-1
         !groups which are associated to a given processor
         icount=0
         do igroup=1,OP2P%ngroup
            if (OP2P%nvctr_par(jproc,igroup,1) > 0) then
               icount=icount+1
               OP2P%igrpr(icount,jproc)=igroup
            end if
         end do

         !calculate the processor which lies after and before the present in the list
         do igroup=1,OP2P%ngroupp(jproc)
            iprocgrs=-1 !<= distance from the processor for ascending order
            iprocgrr=-1 !<= distance from the processor for descending order
            !define the number of data to calculate in total
            do kproc=0,nproc-1
               if (OP2P%nvctr_par(modulo(jproc+kproc,nproc),OP2P%igrpr(igroup,jproc),1) > 0) then
                  iprocgrs=iprocgrs+1
                  OP2P%iprocpm1(1,iprocgrs,igroup,jproc)=modulo(jproc+kproc,nproc)
               end if
               if (OP2P%nvctr_par(modulo(jproc-kproc,nproc),OP2P%igrpr(igroup,jproc),1) > 0) then
                  iprocgrr=iprocgrr+1
                  OP2P%iprocpm1(2,iprocgrr,igroup,jproc)=modulo(jproc-kproc,nproc)
               end if
            end do
         end do
      end do

      !print *,'a',igrpr(:,:),iprocpm1,nvctr_par
      !stop

      !calculate the list of send-receive operations which have to be performed per group
      !allocate it at the maximum size needed (to be changed if symmetric
      allocate(OP2P%communication_schedule(4,0:OP2P%nsteps_max,OP2P%ngroupp_max,0:nproc-1+ndebug),stat=i_stat)
      call memocc(i_stat,OP2P%communication_schedule,'OP2P%communication_schedule',subname)

      !un-initalize array, no communication in any step
      OP2P%communication_schedule(:,:,:,:)=UNINITIALIZED(1)

      do jproc=0,nproc-1
         do igroup=1,OP2P%ngroupp(jproc)
            if (OP2P%forsymop) then
               nstepsm1=(OP2P%nprocgr(OP2P%igrpr(igroup,jproc))-1)/2-1 !here the number of steps should be changed for non-symmetric operation
            else
               nstepsm1=(OP2P%nprocgr(OP2P%igrpr(igroup,jproc))-1)
            end if
            !do not send anything if there is only one member in the group
            if (OP2P%nprocgr(OP2P%igrpr(igroup,jproc)) > 1) then
               do istep=0,nstepsm1
                  !define the arrays for send-receive of data
                  OP2P%communication_schedule(1,istep,igroup,jproc)=OP2P%iprocpm1(2,istep,igroup,jproc)
                  OP2P%communication_schedule(2,istep,igroup,jproc)=OP2P%iprocpm1(2,istep+1,igroup,jproc)
                  if (istep > 0) then
                     OP2P%communication_schedule(3,istep,igroup,jproc)=OP2P%iprocpm1(2,istep,igroup,jproc)
                     OP2P%communication_schedule(4,istep,igroup,jproc)=OP2P%iprocpm1(1,istep,igroup,jproc)
                  end if
               end do
               if (OP2P%forsymop) then
                  istep=nstepsm1+1
                  !the last step behaves differently if the group number is odd or even
                  if (modulo(OP2P%nprocgr(OP2P%igrpr(igroup,jproc)),2) == 0) then
                     OP2P%communication_schedule(1,istep,igroup,jproc)=OP2P%iprocpm1(2,istep,igroup,jproc)
                     OP2P%communication_schedule(2,istep,igroup,jproc)=OP2P%iprocpm1(2,istep+1,igroup,jproc)
                     if (istep > 0) then
                        OP2P%communication_schedule(3,istep,igroup,jproc)=OP2P%iprocpm1(2,istep,igroup,jproc)
                        OP2P%communication_schedule(4,istep,igroup,jproc)=OP2P%iprocpm1(1,istep,igroup,jproc)
                     end if
                  else
                     OP2P%communication_schedule(3,istep,igroup,jproc)=OP2P%iprocpm1(2,istep,igroup,jproc)
                     OP2P%communication_schedule(4,istep,igroup,jproc)=OP2P%iprocpm1(1,istep,igroup,jproc)
                  end if
               end if
            end if
         end do
         !if (iproc==0 .and. jproc==1) print *,'iproc',jproc,OP2P%communication_schedule(3:4,0:1,:,jproc)
      end do

      !after initialisation, the results of the groups can be printed
      if (iproc ==0) call print_group_schemes(6,nproc,OP2P)

   END SUBROUTINE initialize_OP2P_descriptors


   subroutine print_group_schemes(unit,nproc,OP2P)
      implicit none
      !Arguments
      integer, intent(in) :: nproc,unit
      type(OP2P_descriptors), intent(in) :: OP2P
      !local variables
      integer :: jproc,nvctrp,igroup,nsteps,i,istep

      write(unit,'(1x,a,a)')repeat('-',45),'Overlap point-to-point data repartition'
      write(unit,'(1x,8(a))')'| proc |',' No. Groups  |  Grp | #Obj- #Res ',&
         &   '|| N. Components | Steps|    Components   |'
      do jproc=0,nproc-1
         nvctrp=0
         do igroup=1,OP2P%ngroupp(jproc)
            nvctrp=nvctrp+&
               &   OP2P%nvctr_par(jproc,OP2P%igrpr(igroup,jproc),1)
         end do
         !print total number of orbitals and components
         write(unit,'(1x,a,i4,a,i8,a,i13,a)')'| ',jproc,' |',OP2P%ngroupp(jproc),&
            &   repeat(' ',5)//'|'//repeat('-',6)//'|'//repeat('-',12)//'||',&
            &   nvctrp,&
            &   repeat(' ',2)//'|'//repeat('-',6)//'|'//repeat('-',17)//'|'
         !change the values to zero if there is no orbital

         do igroup=1,OP2P%ngroupp(jproc)
            nsteps=0
            do istep=0,OP2P%nsteps_max
               do i=1,4
                  if (OP2P%communication_schedule(i,istep,igroup,jproc) /= UNINITIALIZED(1)) then
                     nsteps=nsteps+1
                     exit
                  end if
               end do
            end do
            write(unit,'(a,i4,a,i5,a,i5,a,i4,a,i8,a,i8,a)')&
               &   ' |'//repeat(' ',6)//'|'//repeat(' ',13)//'|',&
               &   OP2P%igrpr(igroup,jproc),'  |',OP2P%nvctr_par(jproc,OP2P%igrpr(igroup,jproc),1),&
               &   '-',OP2P%nvctr_par(jproc,OP2P%igrpr(igroup,jproc),2),&
               &   ' ||'//repeat(' ',15)//'|',&
               &   nsteps,'  |',1,'-',1,'|'
         end do
      end do
      write(unit,'(1x,a,a)')repeat('-',84)

   END SUBROUTINE print_group_schemes


   subroutine OP2P_communication(iproc,nproc,OP2P,objects_data,results_data,apply_symmetric_operator,&
         &   send_op,receive_op,wait_op) 
      implicit none
      interface
         subroutine send_op(istep,isendproc,irecvproc,ncount,itag,irequest,sendbuf)
            integer, intent(in) :: istep,isendproc,irecvproc,ncount,itag,irequest
            real(kind=8), intent(in) :: sendbuf
         END SUBROUTINE send_op

         subroutine receive_op(istep,isendproc,irecvproc,ncount,itag,irequest,recvbuf)
            integer, intent(in) :: istep,isendproc,irecvproc,ncount,itag,irequest
            real(kind=8), intent(in) :: recvbuf
         END SUBROUTINE receive_op

         subroutine wait_op(iproc,istep,nreq,requests)
            implicit none
            integer, intent(in) :: iproc,istep,nreq
            integer, dimension(nreq), intent(in) :: requests
         END SUBROUTINE wait_op

         !here a procedure can be added if the objects to be sent are not of double precision
         subroutine apply_symmetric_operator(istep,iproc,igroup,remote_result,idata_glob,jdata_glob,&
               &   idata_loc,jdata_loc,ndatai,ndataj,nvctri,nvctrj,nvctri_results,nvctrj_results,&
               &   objects_data_i,objects_data_j,&
               &   results_data_i,results_data_j)
            implicit none
            logical, intent(in) :: remote_result
            integer, intent(in) :: istep,iproc,igroup,idata_glob,jdata_glob,idata_loc,jdata_loc,ndatai,ndataj
            integer, intent(in) :: nvctri,nvctrj,nvctri_results,nvctrj_results
            real(kind=8), dimension(nvctri), intent(in) :: objects_data_i
            real(kind=8), dimension(nvctrj), intent(in) :: objects_data_j 
            real(kind=8), dimension(nvctri_results), intent(inout) :: results_data_i
            real(kind=8), dimension(nvctrj_results), intent(inout) :: results_data_j
         END SUBROUTINE apply_symmetric_operator

      end interface

      integer, intent(in) :: iproc,nproc
      type(OP2P_descriptors), intent(in) :: OP2P
      real(kind=8), dimension(*), intent(inout) :: objects_data !< the dimension of the initial objects is not specified
      real(kind=8), dimension(*), intent(out) :: results_data
      optional :: send_op,receive_op,wait_op

      !local variables
      character(len=*), parameter :: subname='op2p_communication'
      logical :: doit,remote_result
      integer :: isnow,ncommsstep,irnow,igroup,istep,isnow_results,irnow_results,ncommsstep_results,iaddress_local
      integer :: istart_results,i_all,i_stat,nvctri,nvctrj,nvctri_results,nvctrj_results
      integer :: isorb,jsorb,iorbs,jorbs,norbi,norbj,jproc,istart,jproc_to_send,jproc_to_recv,nelems_to_send,nelems_to_recv
      integer, dimension(:,:), allocatable :: mpireq
      real(kind=8), dimension(:,:,:), allocatable :: sendreceive_buffer,restemp_buffer

      !quick return if the number of groups is zero here
      if (OP2P%ngroupp(iproc) == 0) return

      !allocate the array of sendreceive buffers
      allocate(sendreceive_buffer(OP2P%ncomponents_max,2,OP2P%ngroupp(iproc)+ndebug),stat=i_stat)
      call memocc(i_stat,sendreceive_buffer,'sendreceive_buffer',subname)
      !allocate the results buffer only in the case of symmetric OP2P communication
      !if (OP2P%forsym) then
      allocate(restemp_buffer(OP2P%ncomponents_results_max,3,OP2P%ngroupp(iproc)+ndebug),stat=i_stat)
      call memocc(i_stat,restemp_buffer,'restemp_buffer',subname)
      !end if
      allocate(mpireq(2*OP2P%ngroupp_max,2+ndebug),stat=i_stat)
      call memocc(i_stat,mpireq,'mpireq',subname)
      !ncalls=0
      !real communication
      isnow=1
      isnow_results=1
      irnow_results=1
      ncommsstep_results=0
      do istep=0,OP2P%nsteps_max
         !print *,'istep,iproc',istep,iproc,nsteps_max
         irnow=3-isnow
         ncommsstep=0
         !sending receiving data
         do igroup=1,OP2P%ngroupp(iproc)
            if (OP2P%communication_schedule(1,istep,igroup,iproc) /= UNINITIALIZED(1)) then
               ncommsstep=ncommsstep+1 !just check the receives
               jproc_to_send=OP2P%iprocpm1(1,1,igroup,iproc)
               if (OP2P%ngroupp(jproc_to_send) == 0 ) then
                  print *,'ERROR,iproc',iproc
                  stop
               end if
               nelems_to_send=OP2P%nvctr_par(OP2P%communication_schedule(1,istep,igroup,iproc),OP2P%igrpr(igroup,iproc),1) !objects
               iaddress_local=OP2P%ioffp_group(3,iproc,OP2P%igrpr(igroup,iproc))
               if (istep == 0) then
                  if (present(send_op)) then
                     call send_op(istep,iproc,jproc_to_send,nelems_to_send,&
                        &   iproc+2*nproc*istep,mpireq(ncommsstep,1),objects_data(iaddress_local))
                  else
                     call send_mpi(istep,jproc_to_send,nelems_to_send,&
                        &   iproc+2*nproc*istep,mpireq(ncommsstep,1),objects_data(iaddress_local))
                  end if
               else
                  !multiply by two the sending objects for the accumulation array
                  if (present(send_op)) then
                     call send_op(istep,iproc,jproc_to_send,nelems_to_send,&
                        &   iproc+2*nproc*istep,mpireq(ncommsstep,1),sendreceive_buffer(1,isnow,igroup))
                  else
                     call send_mpi(istep,jproc_to_send,nelems_to_send,&
                        &   iproc+2*nproc*istep,mpireq(ncommsstep,1),sendreceive_buffer(1,isnow,igroup))
                  end if
               end if
            end if
            if (OP2P%communication_schedule(2,istep,igroup,iproc) /= UNINITIALIZED(1)) then
               ncommsstep=ncommsstep+1
               jproc_to_recv=OP2P%iprocpm1(2,1,igroup,iproc)
               if (OP2P%ngroupp(jproc_to_recv) == 0 ) then
                  print *,'ERROR(recv),iproc',iproc
                  stop
               end if
               nelems_to_recv=OP2P%nvctr_par(OP2P%communication_schedule(2,istep,igroup,iproc),OP2P%igrpr(igroup,iproc),1) !objects
               iaddress_local=OP2P%ioffp_group(3,iproc,OP2P%igrpr(igroup,iproc))
               if (present(receive_op)) then
                  call receive_op(istep,jproc_to_recv,iproc,nelems_to_recv,&
                     &   jproc_to_recv+2*nproc*istep,mpireq(ncommsstep,1),sendreceive_buffer(1,irnow,igroup))
               else
                  call receive_mpi(istep,jproc_to_recv,iproc,nelems_to_recv,&
                     &   jproc_to_recv+2*nproc*istep,mpireq(ncommsstep,1),sendreceive_buffer(1,irnow,igroup))
               end if
            end if
         end do
         !print *,'starting',iproc,istep,ncommsstep,nsteps_max
         !at this point the first pool of data should have been sent
         !the operation can be performed
         do igroup=1,OP2P%ngroupp(iproc)
            !decide if the operation has to be performed (if something has been received)
            doit = (istep == 0) 
            !to avoid boundary problems in the previous if statement
            if (.not. doit) doit = OP2P%communication_schedule(2,istep-1,igroup,iproc) /=UNINITIALIZED(1)
            if (doit) then
               !starting object (global)
               isorb=OP2P%ioffp_group(1,iproc,OP2P%igrpr(igroup,iproc))
               !starting object (local)
               iorbs=OP2P%ioffp_group(2,iproc,OP2P%igrpr(igroup,iproc))
               !number of objects to be treated locally for iproc
               norbi=OP2P%ioffp_group(4,iproc,OP2P%igrpr(igroup,iproc))
               !number of components associated to the object array
               nvctri=OP2P%nvctr_par(iproc,OP2P%igrpr(igroup,iproc),1)

               if (istep == 0) then
                  jproc=iproc
               else
                  jproc=OP2P%communication_schedule(2,istep-1,igroup,iproc)
               end if
               !starting address on the local array
               istart=OP2P%ioffp_group(3,iproc,OP2P%igrpr(igroup,iproc))

               !starting object (global)
               jsorb=OP2P%ioffp_group(1,jproc,OP2P%igrpr(igroup,iproc))
               !starting object (local)
               jorbs=OP2P%ioffp_group(2,jproc,OP2P%igrpr(igroup,iproc))
               !number of objects to be treated locally for iproc
               norbj=OP2P%ioffp_group(4,jproc,OP2P%igrpr(igroup,iproc))
               !number of components associated to the object array
               nvctrj=OP2P%nvctr_par(jproc,OP2P%igrpr(igroup,iproc),1)

               !starting address on the results_array
               istart_results=OP2P%ioffp_group(5,iproc,OP2P%igrpr(igroup,iproc))
               !number of components associated to the results array
               nvctri_results=OP2P%nvctr_par(iproc,OP2P%igrpr(igroup,iproc),2)
               !number of components associated to the results array
               nvctrj_results=OP2P%nvctr_par(jproc,OP2P%igrpr(igroup,iproc),2)


               !print *,'istep,iproc,istart_results',istep,iproc,istart_results
               !do remote operation for symmetric cases
               remote_result=OP2P%forsymop .and. &
                  &   (istep /=0 .and. OP2P%communication_schedule(3,istep,igroup,iproc) /= UNINITIALIZED(1))
               if (istep > 0) then
                  call apply_symmetric_operator(istep,iproc,OP2P%igrpr(igroup,iproc),remote_result,&
                     &   isorb,jsorb,iorbs,jorbs,norbi,norbj,nvctri,nvctrj,nvctri_results,nvctrj_results,&
                     &   objects_data(istart),sendreceive_buffer(1,isnow,igroup),&
                     &   results_data(istart_results),& !the result may have a different size than the wavefunction
                  restemp_buffer(1,3,igroup))
               else
                  call apply_symmetric_operator(istep,iproc,OP2P%igrpr(igroup,iproc),remote_result,&
                     &   isorb,jsorb,iorbs,jorbs,norbi,norbj,nvctri,nvctrj,nvctri_results,nvctrj_results,&
                     &   objects_data(istart),objects_data(istart),& 
                  results_data(istart_results),& !the result may have a different size than the wavefunction
                  restemp_buffer(1,3,igroup))
               end if
            end if
         end do
         !print *,'ending',iproc,istep

         !check the sending of the results array (this operation can also be performed at the end of the cycle
         if (ncommsstep_results > 0) then
            if (present(wait_op)) then
               call wait_op(iproc,istep,ncommsstep_results,mpireq(1,2))
            else
               call wait_mpi(iproc,istep,ncommsstep_results,mpireq(1,2))
            end if
            !copy the results which have been received (the messages sending are after)
            !to be modified via the creation of an array which follows psi
            do igroup=1,OP2P%ngroupp(iproc)
               !starting address on the local array (supposing that the array and the results have the same size)
               istart_results=OP2P%ioffp_group(5,iproc,OP2P%igrpr(igroup,iproc))
               if (OP2P%communication_schedule(4,istep-1,igroup,iproc) /= UNINITIALIZED(1)) then
                  !reduce temporary result
                  !this routine has to be defined outside from the module
                  call daxpy(OP2P%nvctr_par(iproc,OP2P%igrpr(igroup,iproc),2),1.0d0,restemp_buffer(1,irnow_results,igroup),1,&
                     &   results_data(istart_results),1)
                  !print *,'iproc,reduce',iproc,istep,results_data(istart_results),restemp_buffer(1,irnow_results,igroup)
               end if
            end do
         end if

         ncommsstep_results=0
         !meanwhile, we can receive the result from the processor which has the psi 
         irnow_results=3-isnow_results
         do igroup=1,OP2P%ngroupp(iproc)
            if (OP2P%communication_schedule(3,istep,igroup,iproc) /= UNINITIALIZED(1)) then
               ncommsstep_results=ncommsstep_results+1 !only receive is taken into account
               jproc_to_send=OP2P%communication_schedule(3,istep,igroup,iproc)
               nelems_to_send=OP2P%nvctr_par(OP2P%communication_schedule(3,istep,igroup,iproc),OP2P%igrpr(igroup,iproc),2) !results

               !copy data before sending (to be performed in case the previous send has not finished)
               call dcopy(nelems_to_send,restemp_buffer(1,3,igroup),1,&
                  &   restemp_buffer(1,isnow_results,igroup),1) !is this copy unavoidable? unfortunately yes...

               if (present(send_op)) then
                  call send_op(istep,iproc,jproc_to_send,nelems_to_send,&
                     &   iproc+nproc+2*nproc*istep,mpireq(ncommsstep_results,2),restemp_buffer(1,isnow_results,igroup))
               else
                  call send_mpi(istep,jproc_to_send,nelems_to_send,&
                     &   iproc+nproc+2*nproc*istep,mpireq(ncommsstep_results,2),restemp_buffer(1,isnow_results,igroup))
               end if
            end if
            if (OP2P%communication_schedule(4,istep,igroup,iproc) /= UNINITIALIZED(1)) then
               ncommsstep_results=ncommsstep_results+1
               jproc_to_recv=OP2P%communication_schedule(4,istep,igroup,iproc)
               nelems_to_recv=OP2P%nvctr_par(iproc,OP2P%igrpr(igroup,iproc),2) !results
               if (present(receive_op)) then
                  call receive_op(istep,jproc_to_recv,iproc,nelems_to_recv,&
                     &   jproc_to_recv+nproc+2*nproc*istep,mpireq(ncommsstep_results,2),&
                     &   restemp_buffer(1,irnow_results,igroup))
               else
                  call receive_mpi(istep,jproc_to_recv,iproc,nelems_to_recv,&
                     &   jproc_to_recv+nproc+2*nproc*istep,mpireq(ncommsstep_results,2),&
                     &   restemp_buffer(1,irnow_results,igroup))
               end if
            end if
         end do
         if (istep>1) isnow_results=3-isnow_results

         if (ncommsstep /=0) then
            if (present(wait_op)) then
               call wait_op(iproc,istep,ncommsstep,mpireq(1,1))
            else
               call wait_mpi(iproc,istep,ncommsstep,mpireq(1,1))
            end if
         end if
         isnow=3-isnow
         ncommsstep=0
         !print *,'iproc,istep',iproc,istep
      end do

      i_all=-product(shape(sendreceive_buffer))*kind(sendreceive_buffer)
      deallocate(sendreceive_buffer,stat=i_stat)
      call memocc(i_stat,i_all,'sendreceive_buffer',subname)
      !if (OP2P%forsymop) then
      i_all=-product(shape(restemp_buffer))*kind(restemp_buffer)
      deallocate(restemp_buffer,stat=i_stat)
      call memocc(i_stat,i_all,'restemp_buffer',subname)
      !end if
      i_all=-product(shape(mpireq))*kind(mpireq)
      deallocate(mpireq,stat=i_stat)
      call memocc(i_stat,i_all,'mpireq',subname)

   END SUBROUTINE OP2P_communication

   !> Deallocate everything for OP2P communications
   subroutine free_OP2P_descriptors(OP2P,subname)
      implicit none
      character(len=*), intent(in) :: subname
      type(OP2P_descriptors), intent(inout) :: OP2P
      !local variables
      integer :: i_all,i_stat
      !debugbprint *,'ecco5'

      i_all=-product(shape(OP2P%ngroupp))*kind(OP2P%ngroupp)
      deallocate(OP2P%ngroupp,stat=i_stat)
      call memocc(i_stat,i_all,'ngroupp',subname)
      i_all=-product(shape(OP2P%nprocgr))*kind(OP2P%nprocgr)
      deallocate(OP2P%nprocgr,stat=i_stat)
      call memocc(i_stat,i_all,'nprocgr',subname)
      i_all=-product(shape(OP2P%igrpr))*kind(OP2P%igrpr)
      deallocate(OP2P%igrpr,stat=i_stat)
      call memocc(i_stat,i_all,'igrpr',subname)
      i_all=-product(shape(OP2P%nvctr_par))*kind(OP2P%nvctr_par)
      deallocate(OP2P%nvctr_par,stat=i_stat)
      call memocc(i_stat,i_all,'nvctr_par',subname)
      i_all=-product(shape(OP2P%ioffp_group))*kind(OP2P%ioffp_group)
      deallocate(OP2P%ioffp_group,stat=i_stat)
      call memocc(i_stat,i_all,'ioffp_group',subname)
      i_all=-product(shape(OP2P%iprocpm1))*kind(OP2P%iprocpm1)
      deallocate(OP2P%iprocpm1,stat=i_stat)
      call memocc(i_stat,i_all,'iprocpm1',subname)
      i_all=-product(shape(OP2P%communication_schedule))*kind(OP2P%communication_schedule)
      deallocate(OP2P%communication_schedule,stat=i_stat)
      call memocc(i_stat,i_all,'communication_schedule',subname)

   END SUBROUTINE free_OP2P_descriptors


   subroutine wait_mpi(iproc,istep,nreq,requests)
      implicit none
      integer, intent(in) :: iproc,istep,nreq
      integer, dimension(nreq), intent(in) :: requests
      !local variables
      integer, dimension(MPI_STATUS_SIZE,4) :: mpistat
      !local variables
      integer :: ierr

      !verify that the messages have been passed
      call MPI_WAITALL(nreq,requests,mpistat,ierr)
      if (ierr /=0)  then
         write(*,*) 'ERROR WAITALL, iproc,step,ierr:',iproc,istep,ierr,mpistat !,MPI_STATUSES_IGNORE
      end if
   END SUBROUTINE wait_mpi


   !> fake receiving of the arrays
   subroutine receive_mpi(istep,isendproc,irecvproc,ncount,itag,irequest,recvbuf)
      implicit none
      integer, intent(in) :: istep,isendproc,irecvproc,ncount,itag,irequest
      real(kind=8), intent(in) :: recvbuf
      !local variables
      integer :: ierr

      !here we can add something to trap the IRECV call
      !print '(3(a,i4),i4)','NON_BLOCKING RECV, from',isendproc,' to',irecvproc,', step, elems:',istep,ncount

      call MPI_IRECV(recvbuf,ncount,MPI_DOUBLE_PRECISION,isendproc,&
         &   itag,MPI_COMM_WORLD,irequest,ierr)

      !output error signal
      if (ierr /=0) then
         write(*,*)'ERROR in IRECV, iproc, istep',irecvproc,istep
      end if

   END SUBROUTINE receive_mpi

   !> fake sending of the arrays
   subroutine send_mpi(istep,irecvproc,ncount,itag,irequest,sendbuf) 
      implicit none
      integer, intent(in) :: istep,irecvproc,ncount,itag,irequest
      real(kind=8), intent(in) :: sendbuf
      !local variables
      integer :: ierr

      !here we can add something to trap the ISEND call
      !print '(3(a,i4),i4)','NON_BLOCKING SEND, from',isendproc,' to',irecvproc,', step, elems:',istep,ncount

      call MPI_ISEND(sendbuf,ncount,MPI_DOUBLE_PRECISION,irecvproc,&
         &   itag,MPI_COMM_WORLD,irequest,ierr)

      !output error signal
      if (ierr /=0) then
         write(*,*)'ERROR in ISEND, iproc, istep',irecvproc,istep
      end if
   END SUBROUTINE send_mpi

END MODULE overlap_point_to_point
