!> @file
!!  File defining the routines to initialize the communications between processes
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining routines to initialize the communications
module communications_init
  use communications_base
  implicit none

  private

  public :: init_comms_linear
  public :: init_comms_linear_sumrho
  public :: initialize_communication_potential
  public :: orbitals_communicators


  contains

    subroutine init_comms_linear(iproc, nproc, imethod_overlap, npsidim_orbs, orbs, lzd, nspin, collcom)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, imethod_overlap, npsidim_orbs, nspin
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(comms_linear),intent(inout) :: collcom
      
      ! Local variables
      integer :: iorb, iiorb, ilr, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, i3, ii3, jj3min, jj3max, jj3
      integer :: ipt, nvalp_c, nvalp_f, i3s, n3p, ii, i, jjproc, jproc, ii3min, ii3max, np, n1p1, iseg, j0, j1, ierr
      integer :: j3start, j3end
      real(kind=8),dimension(:,:,:),allocatable :: weightppp_c, weightppp_f
      real(kind=8) :: weight_c_tot, weight_f_tot, weightp_c, weightp_f, tt
      integer,dimension(:,:),allocatable :: istartend_c, istartend_f
      integer,dimension(:,:,:),allocatable :: index_in_global_c, index_in_global_f
      real(kind=8),dimension(:,:,:),allocatable :: weightloc_c, weightloc_f
      integer :: window_c, window_f, i3start, i3end
      real(kind=8) :: weight_c_tot_check, weight_f_tot_check
      
      real(kind=4) :: tr0, tr1, trt0, trt1
      real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
      logical, parameter :: extra_timing=.false.

      call timing(iproc,'init_collcomm ','ON')
      if (extra_timing) call cpu_time(trt0)   
      call f_routine('init_comms_linear')

      ! method to calculate the overlap
      collcom%imethod_overlap = imethod_overlap

      ! Split up the z dimension in disjoint pieces.
      tt = real(lzd%glr%d%n3+1,kind=8)/real(nproc,kind=8)
      ii = floor(tt)
      n3p = ii
      jjproc = lzd%glr%d%n3+1 - nproc*ii
      if (iproc<=jjproc-1) n3p = n3p + 1
      i=1
      do jproc=0,nproc-1
          if (iproc==jproc) i3s = i
          i = i + ii
          if (jproc<=jjproc-1) i = i + 1
      end do


      i3start=1000000000
      i3end=-1000000000
      !j3start=1000000000
      !j3end=-1000000000
      do iorb=1,orbs%norbp
          iiorb = orbs%isorb+iorb
          !if (orbs%spinsgn(iiorb)<0.d0) cycle !consider only up orbitals
          ilr = orbs%inwhichlocreg(iiorb)
          i3start = min(i3start,lzd%llr(ilr)%ns3)
          i3end = max(i3end,lzd%llr(ilr)%ns3+lzd%llr(ilr)%d%n3)
          !j3start = min(j3start,modulo(lzd%llr(ilr)%ns3-lzd%llr(ilr)%ns3,lzd%llr(ilr)%d%n3)+1)
          !j3end = max(j3end,modulo(lzd%llr(ilr)%ns3+lzd%llr(ilr)%d%n3-lzd%llr(ilr)%ns3,lzd%llr(ilr)%d%n3)+1)
      end do
      if (orbs%norbp==0) then!.or.i3start==1000000000) then ! need to account for the case when norbp/=0 but all orbitals were down but should probably do in a better way
         i3end=1
         i3start=0
      end if

      !j3start = modulo(i3start-i3start,(lzd%glr%d%n3+1))+1 !should give 1
      !j3end = modulo(i3end-i3start,(lzd%glr%d%n3+1))+1
      j3start = i3start-i3start+1 !should give 1
      j3end = i3end-i3start+1
      ! Shrink if the extent is larger than the box
      j3end=min(j3end,lzd%glr%d%n3+1)

      !if (i3end-i3start/=j3end-j3start) stop 'i3end-i3start/=j3end-j3start'

      !write(*,'(a,9i8)') 'iproc, i3s, n3p, ii3min, ii3max, i3start, i3end, j3start, j3end', iproc, i3s, n3p, ii3min, ii3max, i3start, i3end, j3start, j3end
      !write(*,'(a,7i8)') 'iproc, i3s, n3p, i3start, i3end, j3start, j3end', iproc, i3s, n3p, i3start, i3end, j3start, j3end

      !!i3start=j3start
      !!i3end=j3end

      !weightloc_c = f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,1.to.(i3end-i3start+1)/),id='weightloc_c')
      !weightloc_f = f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,1.to.(i3end-i3start+1)/),id='weightloc_f')
      weightloc_c = f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,j3start.to.j3end/),id='weightloc_c')
      weightloc_f = f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,j3start.to.j3end/),id='weightloc_f')



      !!! Determine the maximal extent in the z direction that iproc has to handle
      !!ii3min = 1000000000
      !!ii3max = -1000000000
      !!jj3min = 1000000000
      !!jj3max = -1000000000
      !!do iorb=1,orbs%norbp
      !!    iiorb=orbs%isorb+iorb
      !!    ilr=orbs%inwhichlocreg(iiorb)
      !!    if (lzd%llr(ilr)%wfd%nseg_c>0) then
      !!        !n1p1=lzd%llr(ilr)%d%n1+1
      !!        n1p1=lzd%glr%d%n1+1
      !!        !np=n1p1*(lzd%llr(ilr)%d%n2+1)
      !!        np=n1p1*(lzd%glr%d%n2+1)
      !!        do iseg=1,lzd%llr(ilr)%wfd%nseg_c
      !!            !write(*,*) 'keygloc(1,iseg), keyglob(1,iseg)', lzd%llr(ilr)%wfd%keygloc(1,iseg), lzd%llr(ilr)%wfd%keyglob(1,iseg)
      !!            j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
      !!            j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
      !!            ii=j0-1
      !!            i3=ii/np
      !!            !ii3=i3+lzd%llr(ilr)%ns3
      !!            ii3=i3+lzd%glr%ns3
      !!            ii3min = min(ii3min,ii3)
      !!            ii3max = max(ii3max,ii3)
      !!            jj3 = modulo(ii3-lzd%llr(ilr)%ns3-1,lzd%glr%d%n3)+1
      !!            jj3min = min(jj3min,jj3)
      !!            jj3max = max(jj3max,jj3)
      !!        end do
      !!    end if
      !!end do
      
      !!write(*,*) 'ii3min, ii3max, jj3min, jj3max', ii3min, ii3max, jj3min, jj3max

      ii3min=i3start
      ii3max=i3end
    
      !!index_in_global_c=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,0.to.lzd%glr%d%n3/),id='index_in_global_c')
      !!index_in_global_f=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,0.to.lzd%glr%d%n3/),id='index_in_global_f')
      !index_in_global_c=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,ii3min.to.ii3max/),id='index_in_global_c')
      !index_in_global_f=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,ii3min.to.ii3max/),id='index_in_global_f')
      !index_in_global_c=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,i3start.to.i3end/),id='index_in_global_c')
      !index_in_global_f=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,i3start.to.i3end/),id='index_in_global_f')
      index_in_global_c=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,j3start.to.j3end/),id='index_in_global_c')
      index_in_global_f=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,j3start.to.j3end/),id='index_in_global_f')

      weightppp_c=f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,1.to.max(1,n3p)/),id='weightppp_c')
      weightppp_f=f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,1.to.max(1,n3p)/),id='weightppp_c')

      
    
      call get_weights(iproc, nproc, orbs, lzd, i3s, n3p, i3start, i3end, j3start, j3end, &
           weightloc_c, weightloc_f, window_c, window_f, weightppp_c, weightppp_f, &
           weight_c_tot_check, weight_f_tot_check)
    
      ! Assign the grid points to the processes such that the work is equally distributed
      istartend_c=f_malloc((/1.to.2,0.to.nproc-1/),id='istartend_c')
      istartend_f=f_malloc((/1.to.2,0.to.nproc-1/),id='istartend_f')

      if (extra_timing) call cpu_time(tr0)
      !call assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
      !     istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
      !     weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, nvalp_c, nvalp_f)
      call assign_weight_to_process(iproc, nproc, lzd, i3s, n3p, window_c, window_f, &
           weight_c_tot_check, weight_f_tot_check, &
           weightppp_c, weightppp_f, weight_c_tot, weight_f_tot, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, nvalp_c, nvalp_f)

      call f_free(weightloc_c)
      call f_free(weightloc_f)
     
      if (extra_timing) call cpu_time(tr1)
      if (extra_timing) time0=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0)
      !call assign_weight_to_process_new(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
      !     istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
      !     weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, nvalp_c, nvalp_f)
      if (extra_timing) call cpu_time(tr1)
      if (extra_timing) time1=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0)   
      ! Determine the index of a grid point i1,i2,i3 in the compressed array
      call get_index_in_global2(lzd%glr, j3start, j3end, i3start, index_in_global_c, index_in_global_f)
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time2=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0) 
      ! Determine values for mpi_alltoallv
      call allocate_MPI_communication_arrays(nproc, collcom)
      call determine_communication_arrays(iproc, nproc, npsidim_orbs, orbs, nspin, lzd, istartend_c, istartend_f, &
           j3start, j3end, i3start, index_in_global_c, index_in_global_f, nvalp_c, nvalp_f, &
           collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
           collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f)
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time3=real(tr1-tr0,kind=8)
    
      !Now set some integers in the collcomm structure
      collcom%ndimind_c = sum(collcom%nrecvcounts_c)
      collcom%ndimind_f = sum(collcom%nrecvcounts_f)
    
      ! Now rearrange the data on the process to communicate them
      collcom%ndimpsi_c=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          collcom%ndimpsi_c=collcom%ndimpsi_c+lzd%llr(ilr)%wfd%nvctr_c
      end do
      collcom%ndimpsi_f=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          collcom%ndimpsi_f=collcom%ndimpsi_f+lzd%llr(ilr)%wfd%nvctr_f
      end do
    
      call allocate_local_comms_cubic(collcom)
    
      if (extra_timing) call cpu_time(tr0) 
      call determine_num_orbs_per_gridpoint_new(iproc, nproc, lzd, i3s, n3p, weightppp_c, weightppp_f, &
           i3start, istartend_c, istartend_f, &
           istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, &
           collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time4=real(tr1-tr0,kind=8)
      call f_free(weightppp_c)
      call f_free(weightppp_f)
      if (extra_timing) call cpu_time(tr0)   
      call get_switch_indices(iproc, nproc, orbs, lzd, nspin, &
           collcom%nptsp_c, collcom%nptsp_f, collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f, &
           collcom%ndimpsi_c, collcom%ndimpsi_f, istartend_c, istartend_f, &
           istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%ndimind_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
           collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%ndimind_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f, &
           j3start, j3end, i3start, index_in_global_c, index_in_global_f, &
           weightp_c, weightp_f, collcom%isendbuf_c, collcom%irecvbuf_c, collcom%isendbuf_f, collcom%irecvbuf_f, &
           collcom%indexrecvorbital_c, collcom%iextract_c, collcom%iexpand_c, &
           collcom%indexrecvorbital_f, collcom%iextract_f, collcom%iexpand_f)
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time5=real(tr1-tr0,kind=8)
    
      ! These variables are used in various subroutines to speed up the code
      collcom%isptsp_c(1) = 0
      do ipt=2,collcom%nptsp_c
            collcom%isptsp_c(ipt) = collcom%isptsp_c(ipt-1) + collcom%norb_per_gridpoint_c(ipt-1)
      end do
      if (maxval(collcom%isptsp_c)>collcom%ndimind_c) stop 'maxval(collcom%isptsp_c)>collcom%ndimind_c'
    
      collcom%isptsp_f(1) = 0
      do ipt=2,collcom%nptsp_f
            collcom%isptsp_f(ipt) = collcom%isptsp_f(ipt-1) + collcom%norb_per_gridpoint_f(ipt-1)
      end do
      if (maxval(collcom%isptsp_f)>collcom%ndimind_f) stop 'maxval(collcom%isptsp_f)>collcom%ndimind_f'
    
      ! Not used any more, so deallocate...
      call f_free(istartend_c)
      call f_free(istartend_f)
    
      call f_free(index_in_global_c)
      call f_free(index_in_global_f)
        
      call f_release_routine()
      
      call timing(iproc,'init_collcomm ','OF')
      if (extra_timing) call cpu_time(trt1)   
      if (extra_timing) ttime=real(trt1-trt0,kind=8)

      if (extra_timing.and.iproc==0) print*,'time0,time1',time0,time1,time2,time3,time4,time5,&
           time0+time1+time2+time3+time4+time5,ttime
  
    end subroutine init_comms_linear


    subroutine get_weights(iproc, nproc, orbs, lzd, i3s, n3p, i3start, i3end, j3start, j3end, &
               weightloc_c, weightloc_f, window_c, window_f, weightppp_c, weightppp_f, &
               weight_c_tot_check, weight_f_tot_check)
      use module_base
      use module_types
      use locregs, only: get_extent_of_overlap
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, i3s, n3p, i3start, i3end, j3start, j3end
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      !real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,1:(i3end-i3start+1)),intent(inout) :: weightloc_c, weightloc_f
      !real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,i3start:i3end),intent(inout) :: weightloc_c, weightloc_f
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,j3start:j3end),intent(inout) :: weightloc_c, weightloc_f
      integer,intent(inout) :: window_c, window_f
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,1:max(1,n3p)),intent(out) :: weightppp_c, weightppp_f
      real(kind=8),intent(out) :: weight_c_tot_check, weight_f_tot_check
      
      ! Local variables
      integer :: iorb, iiorb, i0, i1, i2, i3, ii, iseg, ilr, istart, iend, i, j0, j1, ii1, ii2, ii3, n1p1, np
      integer :: i3e, ii3s, ii3e, is, ie, size_of_double, ierr, info, window, jproc, ncount, k, n, js, je, imin, imax
      integer :: request_c, request_f, jj3, isize
      real(kind=8),dimension(:),allocatable :: reducearr
      !real(kind=8),dimension(:,:,:),allocatable :: weightloc
      integer,dimension(:,:),allocatable :: i3startend
      real(kind=8) :: tt
      integer,dimension(2) :: ks, ke, nlen
      logical :: communicate
    
      call f_routine(id='get_weights')
    
      ii=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
      !weight_c_tot=0.d0
      !weight_f_tot=0.d0
    
   
!      orbs_it=>orbital_iterator(orbs)
!      do while(associated(orbs_it))
!        iorb=get_absolute_orbital(orbs_it)
!        ilr=get_orbital_locreg(orbs_it)
!        
!        [....]
!
!        orbs_it=>orbital_next(orbs_it)
!      end do

    
      i3startend = f_malloc0((/1.to.4,0.to.nproc-1/),id='i3startend')
      i3startend(1,iproc) = i3start+1
      i3startend(2,iproc) = i3end+1
      i3startend(3,iproc) = i3s
      i3startend(4,iproc) = i3s+n3p-1
      if (nproc>1) then
          call mpiallred(i3startend, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if


      !@NEW ##################################
      ! coarse part


      !!i3start=1000000000
      !!i3end=-1000000000
      !!do iorb=1,orbs%norbp
      !!    iiorb = orbs%isorb+iorb
      !!    if (orbs%spinsgn(iiorb)<0.d0) cycle !consider only up orbitals
      !!    ilr = orbs%inwhichlocreg(iiorb)
      !!    i3start = min(i3start,lzd%llr(ilr)%ns3)
      !!    i3end = max(i3end,lzd%llr(ilr)%ns3+lzd%llr(ilr)%d%n3)
      !!end do
      !!if (orbs%norbp==0) then
      !!   !want i3end-i3start+1=0
      !!   !i3start+1>lzd%glr%d%n3+1 or 1>i3end+1
      !!   i3end=0
      !!   i3start=1
      !!end if

      !weightloc = f_malloc0((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,1.to.(i3end-i3start+1)/),id='weightloc')
      ncount = (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      reducearr = f_malloc(ncount,id='reducearr')


      !call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p,weightppp_c(0,0,1))
      i3e=i3s+n3p-1
      isize = 0
      communicate = .false.
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          if (orbs%spinsgn(iiorb)<0.d0) cycle !consider only up orbitals
          ! If one is at least once beyind this cycle, the communication later on can be performed
          communicate = .true.
          ilr = orbs%inwhichlocreg(iiorb)
          isize = isize + lzd%llr(ilr)%wfd%nvctr_c
          ii3s = lzd%llr(ilr)%ns3
          ii3e = ii3s + lzd%llr(ilr)%d%n3
          !!write(*,'(a,6i8)') 'init: iproc, iorb, ii3s, ii3e, i3s, i3e', iproc, iorb, ii3s, ii3e, i3s, i3e
          !if (ii3s+1>i3e .or. ii3e+1<i3s) cycle !+1 since ns3 starts at 0, but is3 at 1

          !n1p1=lzd%llr(ilr)%d%n1+1
          !np=n1p1*(lzd%llr(ilr)%d%n2+1)
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)

          if (lzd%llr(ilr)%wfd%nseg_c>0) then
              !!$omp do
              imin=10000000
              imax=-10000000
              do iseg=1,lzd%llr(ilr)%wfd%nseg_c
                  j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/np
                  !ii3=i3+lzd%llr(ilr)%ns3
                  ii3=i3
                  if (ii3<imin) imin=ii3
                  if (ii3>imax) imax=ii3
                  !!if (ii3>i3end) stop 'strange 1'
                  !!if (ii3<i3start) stop 'strange 2'
                  !if (ii3+i3start>i3end) then
                  !    write(*,*) 'ii3, i3start, i3end', ii3, i3start, i3end
                  !    stop 'strange 1'
                  !end if
                  !if (ii3+i3start<i3start) then
                  !    write(*,*) 'ii3, i3start, i3end', ii3, i3start, i3end
                  !    stop 'strange 2'
                  !end if
                  jj3=modulo(ii3-i3start,(lzd%glr%d%n3+1))+1
                  if (jj3>j3end) then
                      write(*,'(a,5i8)') 'ii3, i3start, lzd%glr%d%n3, jj3, j3end', ii3, i3start, lzd%glr%d%n3, jj3, j3end
                      stop 'strange 2.1'
                  end if
                  if (jj3<j3start) then
                      write(*,'(a,5i8)') 'ii3, i3start, lzd%glr%d%n3, jj3, j3start', ii3, i3start, lzd%glr%d%n3, jj3, j3start
                      stop 'strange 2.2'
                  end if
                  !if (ii3+1<i3s) cycle
                  !if (ii3+1>i3e) exit
                  ii=ii-i3*np
                  i2=ii/n1p1
                  i0=ii-i2*n1p1
                  i1=i0+j1-j0
                  !!write(400,'(a,9i8)') 'j0, j1, ii, i0, i1, i2, i3, ii3, jj3',j0,j1,ii,i0,i1,i2,i3,ii3,jj3
                  !ii2=i2+lzd%llr(ilr)%ns2
                  ii2=i2
                  do i=i0,i1
                      !ii1=i+lzd%llr(ilr)%ns1
                      ii1=i
                      !weightppp_c(ii1,ii2,ii3+1-i3s+1)=weightppp_c(ii1,ii2,ii3+1-i3s+1)+1.d0
                      !weightloc_c(ii1,ii2,ii3-i3start+1)=weightloc_c(ii1,ii2,ii3-i3start+1)+1.d0
                      weightloc_c(ii1,ii2,jj3)=weightloc_c(ii1,ii2,jj3)+1.d0
                      !weight_c_tot=weight_c_tot+1.d0
                  end do
              end do
              !!write(*,'(a,4i8)') 'iproc, ilr, imin, imax', iproc, ilr, imin, imax
              !!$omp end do
          end if
      end do

      ! First a local check, then reduction for later
      weight_c_tot_check = sum(weightloc_c)
      if (nint(weight_c_tot_check)/=isize) then
          write(*,'(a,2i12)') 'weight_c_tot_check, isize', nint(weight_c_tot_check), isize
          stop 'weight_c_tot_check/=isize'
      end if

      if (nproc>1) then
          call mpiallred(weight_c_tot_check, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if



      !!tt = sum(weightloc_c)
      !!call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm)
      !!write(*,*) 'tt1',tt
      


      
      !do i3=1,lzd%glr%d%n3+1
      !    ! Check whether this slice has been (partially) calculated by iproc,
      !    ! otherwise fill with zero
      !    if (i3start+1<=i3 .and. i3<=i3end+1) then
      !        call vcopy(ncount, weightloc(0,0,i3-i3start), 1, reducearr(1), 1)
      !    else
      !        call f_zero(reducearr)
      !    end if

      !    ! Communicate the slice and the zeros (a bit wasteful...)
      !    if (nproc>1) then
      !        call mpiallred(reducearr(1), ncount, mpi_sum, bigdft_mpi%mpi_comm)
      !    end if

      !    ! Check whether iproc needs this slice
      !    if (i3s<=i3 .and. i3<=i3s+n3p-1) then
      !        call vcopy(ncount, reducearr(1), 1, weightppp_c(0,0,i3-i3s+1), 1)
      !    end if
      !end do

      !!do i3=1,max(1,n3p)
      !!    do i2=0,lzd%glr%d%n2
      !!        do i1=0,lzd%glr%d%n1
      !!            write(1000+iproc,*) i1, i2, i3, weightppp_c(i1,i2,i3)
      !!        end do
      !!    end do
      !!end do

      !@NEW #########################################
      !!weightppp_c = 0.d0
      !!call mpi_type_size(mpi_double_precision, size_of_double, ierr)
      !!call mpi_info_create(info, ierr)
      !!call mpi_info_set(info, "no_locks", "true", ierr)
      !!call mpi_win_create(weightppp_c(0,0,1), &
      !!     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p*size_of_double,kind=mpi_address_kind), size_of_double, &
      !!     info, bigdft_mpi%mpi_comm, window_c, ierr)
      !!call mpi_info_free(info, ierr)
      !!call mpi_win_fence(mpi_mode_noprecede, window_c, ierr)

      if (nproc>1) then
          window_c = mpiwindow((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p, weightppp_c(0,0,1), bigdft_mpi%mpi_comm)
      end if


      !write(*,*) 'sum(weightloc_c)', sum(weightloc_c)
      !do i3=j3start,j3end
      !  do i2=0,lzd%glr%d%n2
      !    do i1=0,lzd%glr%d%n1
      !      write(1200,*) 'i1,i2,i3,val',i1,i2,i3+i3start-1,weightloc_c(i1,i2,i3)
      !    end do
      !  end do
      !end do

      !if (nproc>1) then
      if (communicate) then
          do jproc=0,nproc-1
              !Check whether there is an overlap
              ! Start and end on task iproc, possibly out of box
              !!is = max(i3startend(1,iproc),i3startend(3,jproc))
              !!ie = min(i3startend(2,iproc),i3startend(4,jproc))
              !!write(*,'(a,3i8)') 'i3startend(1,iproc),i3startend(3,jproc), is', i3startend(1,iproc),i3startend(3,jproc), is
              !!write(*,'(a,3i8)') 'i3startend(2,iproc),i3startend(4,jproc), ie', i3startend(2,iproc),i3startend(4,jproc), ie
              ! The min is for cases where a task has more than the entire box
              is=modulo(i3startend(1,iproc)-1,lzd%glr%d%n3+1)+1
              if (i3startend(2,iproc)-i3startend(1,iproc)>lzd%glr%d%n3) then 
                  ie=modulo(min(is-1,i3startend(2,iproc))-1,lzd%glr%d%n3+1)+1
              else
                  ie=modulo(i3startend(2,iproc)-1,lzd%glr%d%n3+1)+1
              end if
              js=i3startend(3,jproc)
              je=i3startend(4,jproc)
              if (je>=js) then
                  call get_extent_of_overlap(is, ie, js, je, n, ks, ke, nlen)
                  !write(*,'(a,11i7)') 'is, ie, js, je, n, ks, ke, nlen', is, ie, js, je, n, ks, ke, nlen
                  do k=1,n
                      ! Undo the periodic wrap around if required
                      !if (ks(k)>i3end) then
                      !    ii=ks(k)-(lzd%glr%d%n3+1)
                      !else
                      !    ii=ks(k)
                      !end if
                      ii=modulo(ks(k)-i3start-1,(lzd%glr%d%n3+1))+1
                      !write(*,'(a,9i9)') 'k, ks(k), ke(k), nlen(k), i3start, ii, ks(k), i3startend(3,jproc), n3p', k, ks(k), ke(k), nlen(k), i3start, ii, ks(k), i3startend(3,jproc), n3p
                      if (nproc>1) then
                          call mpiaccumulate(weightloc_c(0,0,ii), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), &
                               jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ks(k)-i3startend(3,jproc)),kind=mpi_address_kind), &
                               (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), mpi_sum, window_c)
                      else
                          call axpy((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), 1.d0, weightloc_c(0,0,ii), 1, &
                               weightppp_c(0,0,1+(ks(k)-i3startend(3,jproc))), 1)
                      end if
                      !!call mpiaccumulate(weightloc_c(0,0,ks(k)-modulo(i3start-1,lzd%glr%d%n1+1)+1), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), &
                      !!     jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ks(k)-i3startend(3,jproc)),kind=mpi_address_kind), &
                      !!     (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), mpi_sum, window_c)
                  end do
              end if
              !if (ie-is>=0) then
              !    !!call mpi_accumulate(weightloc_c(0,0,is-i3start), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), &
              !    !!     mpi_double_precision, jproc, &
              !    !!     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(is-i3startend(3,jproc)),kind=mpi_address_kind), &
              !    !!     (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), mpi_double_precision, &
              !    !!     mpi_sum, window_c, ierr)
              !    call mpiaccumulate(weightloc_c(0,0,is-i3start), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), &
              !         jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(is-i3startend(3,jproc)),kind=mpi_address_kind), &
              !         (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), mpi_sum, window_c)
              !end if
          end do
      end if
      !else
      !    !!is = i3startend(1,iproc)
      !    !!ie = i3startend(2,iproc)
      !    !!call vcopy((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), weightloc_c(0,0,is-i3start), 1, &
      !    !!     weightppp_c(0,0,is-i3startend(3,iproc)+1), 1)
      !    call f_memcpy(n=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(j3end-j3start+1), src=weightloc_c(0,0,j3start), dest=weightppp_c(0,0,1))
      !end if


      !call f_free(i3startend)
      !call mpi_win_fence(0, window_c, ierr)
      !call mpi_win_free(window_c, ierr)
      !@END NEW #####################################
      !!do i3=1,max(1,n3p)
      !!    do i2=0,lzd%glr%d%n2
      !!        do i1=0,lzd%glr%d%n1
      !!            write(2000+iproc,*) i1, i2, i3, weightppp_c(i1,i2,i3)
      !!        end do
      !!    end do
      !!end do

      !!call mpi_type_size(mpi_double_precision, size_of_double, ierr)
      !!call mpi_info_create(info, ierr)
      !!call mpi_info_set(info, "no_locks", "true", ierr)
      !!call mpi_win_create(weightppp_c(0,0,1), &
      !!     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p*size_of_double,kind=mpi_address_kind), size_of_double, &
      !!     info, bigdft_mpi%mpi_comm, window, ierr)
      !!call mpi_info_free(info, ierr)
      !!call mpi_win_fence(mpi_mode_noprecede, window, ierr)
      !!dummybuf = f_malloc((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p,id='dummybuf')
      !!i3startend = f_malloc0((/1.to.2,0.to.nproc-1/),id='i3startend')
      !!call mpiallred(i3startend(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm)
      !!do jproc=0,nproc-1
      !!    !Check whether there is an overlap
      !!    is = max(i3startend(1,iproc),i3startend(1,jproc))
      !!    ie = min(i3startend(2,iproc),i3startend(2,jproc))
      !!    if (ie-is>=0) then
      !!        call mpi_fetch_and_op(weightloc_c(0,0,i2-i3startend(1,iproc)+1), dummybuf(1), &
      !!             mpi_double_precision, jproc, &
      !!             int(is-i3startend(1,jproc),kind=mpi_address_kind), mpi_sum, window, ierr)
      !!    end if
      !!end do
      !!call mpi_win_fence(0, window, ierr)
      !!call mpi_win_free(window, ierr)
      !!call f_free(dummybuf)

      !weight_c_tot = 0.d0
      !do i3=1,n3p
      !    do i2=0,lzd%glr%d%n2
      !        do i1=0,lzd%glr%d%n1
      !            weightppp_c(i1,i2,i3)=weightppp_c(i1,i2,i3)**2
      !            weight_c_tot = weight_c_tot + weightppp_c(i1,i2,i3)
      !        end do
      !    end do
      !end do
      !if (nproc>1) then
      !    call mpiallred(weight_c_tot, 1, mpi_sum, bigdft_mpi%mpi_comm)
      !end if



      ! fine part
      !if (i3end-i3start>=0) then
      !    call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(i3end-i3start+1), weightloc_f(0,0,1))
      !end if
      !call f_zero(weightloc)
      !call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p,weightppp_f(0,0,1))
      i3e=i3s+n3p-1
      isize = 0
      communicate = .false.
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          if (orbs%spinsgn(iiorb)<0.d0) cycle !consider only up orbitals
          ! If one is at least once beyind this cycle, the communication later on can be performed
          communicate = .true.
          ilr = orbs%inwhichlocreg(iiorb)
          isize = isize + lzd%llr(ilr)%wfd%nvctr_f
          ii3s = lzd%llr(ilr)%ns3
          ii3e = ii3s + lzd%llr(ilr)%d%n3
          !!write(*,'(a,6i8)') 'init: iproc, iorb, ii3s, ii3e, i3s, i3e', iproc, iorb, ii3s, ii3e, i3s, i3e
          !if (ii3s+1>i3e .or. ii3e+1<i3s) cycle !+1 since ns3 starts at 0, but is3 at 1

          !n1p1=lzd%llr(ilr)%d%n1+1
          !np=n1p1*(lzd%llr(ilr)%d%n2+1)
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)

          istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
          iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
          if (istart<=iend) then
              !!$omp do
              do iseg=istart,iend
                  j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/np
                  !ii3=i3+lzd%llr(ilr)%ns3
                  ii3=i3
                  !if (ii3>i3end) stop 'strange 1'
                  !if (ii3<i3start) stop 'strange 2'
                  !!if (ii3+i3start>i3end) stop 'strange 1'
                  !!if (ii3+i3start<i3start) stop 'strange 2'
                  !if (ii3+1<i3s) cycle
                  !if (ii3+1>i3e) exit
                  jj3=modulo(ii3-i3start,(lzd%glr%d%n3+1))+1
                  !!if (jj3>i3end) stop 'strange 1'
                  !!if (jj3<i3start) stop 'strange 2'
                  if (jj3>j3end) then
                      write(*,'(a,5i8)') 'ii3, i3start, lzd%glr%d%n3, jj3, j3end', ii3, i3start, lzd%glr%d%n3, jj3, j3end
                      stop 'strange 1.1'
                  end if
                  if (jj3<j3start) then
                      write(*,'(a,5i8)') 'ii3, i3start, lzd%glr%d%n3, jj3, j3start', ii3, i3start, lzd%glr%d%n3, jj3, j3start
                      stop 'strange 1.2'
                  end if
                  ii=ii-i3*np
                  i2=ii/n1p1
                  i0=ii-i2*n1p1
                  i1=i0+j1-j0
                  !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
                  !ii2=i2+lzd%llr(ilr)%ns2
                  ii2=i2
                  do i=i0,i1
                      !ii1=i+lzd%llr(ilr)%ns1
                      ii1=i
                      !weightppp_f(ii1,ii2,ii3+1-i3s+1)=weightppp_f(ii1,ii2,ii3+1-i3s+1)+1.d0
                      !weightloc_f(ii1,ii2,ii3-i3start+1)=weightloc_f(ii1,ii2,ii3-i3start+1)+1.d0
                      weightloc_f(ii1,ii2,jj3)=weightloc_f(ii1,ii2,jj3)+1.d0
                      !if (ii1==33 .and. ii2==33 .and. jj3==10) write(*,'(a,3i8,f11.1)') 'NONZERO, j0, j1, ii3, wl(ii1,ii2,jj3)', j0, j1, ii3, weightloc_f(ii1,ii2,jj3)
                      !if (ii1==33 .and. ii2==33 .and. ii3==125) write(*,'(a,3i8,f11.1)') 'NONZERO, j0, j1, jj3, wl(ii1,ii2,jj3)', j0, j1, jj3, weightloc_f(ii1,ii2,jj3)
                  end do
              end do
              !!$omp end do
          end if
      end do

      ! First a local check, then reduction for later
      weight_f_tot_check = sum(weightloc_f)
      if (nint(weight_f_tot_check)/=isize) then
          write(*,'(a,2i12)') 'weight_f_tot_check, isize', nint(weight_f_tot_check), isize
          stop 'weight_f_tot_check/=isize'
      end if

      if (nproc>1) then
          call mpiallred(weight_f_tot_check, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if


      !write(*,*) 'sum(weightloc_f)', sum(weightloc_f)
      !do i3=j3start,j3end
      !  jj3=modulo(i3-i3start,lzd%glr%d%n3)+1
      !  do i2=0,lzd%glr%d%n2
      !    do i1=0,lzd%glr%d%n1
      !      write(1201,'(a,4i8,f12.2)') 'i1,i2,i3,jj3,val',i1,i2,i3,jj3,weightloc_f(i1,i2,i3)
      !    end do
      !  end do
      !end do

      !do i3=1,lzd%glr%d%n3+1
      !    ! Check whether this slice has been (partially) calculated by iproc,
      !    ! otherwise fill with zero
      !    if (i3start+1<=i3 .and. i3<=i3end+1) then
      !        call vcopy(ncount, weightloc(0,0,i3-i3start), 1, reducearr(1), 1)
      !    else
      !        call to_zero(ncount, reducearr(1))
      !    end if

      !    ! Communicate the slice and the zeros (a bit wasteful...)
      !    if (nproc>1) then
      !        call mpiallred(reducearr(1), ncount, mpi_sum, bigdft_mpi%mpi_comm)
      !    end if

      !    ! Check whether iproc needs this slice
      !    if (i3s<=i3 .and. i3<=i3s+n3p-1) then
      !        call vcopy(ncount, reducearr(1), 1, weightppp_f(0,0,i3-i3s+1), 1)
      !    end if
      !end do

      !@NEW #########################################
      !call mpi_type_size(mpi_double_precision, size_of_double, ierr)
      !call mpi_info_create(info, ierr)
      !call mpi_info_set(info, "no_locks", "true", ierr)
      !call mpi_win_create(weightppp_f(0,0,1), &
      !     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p*size_of_double,kind=mpi_address_kind), size_of_double, &
      !     info, bigdft_mpi%mpi_comm, window_f, ierr)
      !call mpi_info_free(info, ierr)
      !call mpi_win_fence(mpi_mode_noprecede, window_f, ierr)

      if (nproc>1) then
          window_f = mpiwindow((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p, weightppp_f(0,0,1), bigdft_mpi%mpi_comm)
      end if

      !!i3startend = f_malloc0((/1.to.4,0.to.nproc-1/),id='i3startend')
      !!i3startend(1,iproc) = i3start+1
      !!i3startend(2,iproc) = i3end+1
      !!i3startend(3,iproc) = i3s
      !!i3startend(4,iproc) = i3s+n3p-1
      !!call mpiallred(i3startend(1,0), 4*nproc, mpi_sum, bigdft_mpi%mpi_comm)

      if (communicate) then
          do jproc=0,nproc-1
              !Check whether there is an overlap
              ! Start and end on task iproc, possibly out of box
              !!is = max(i3startend(1,iproc),i3startend(3,jproc))
              !!ie = min(i3startend(2,iproc),i3startend(4,jproc))
              !!write(*,'(a,3i8)') 'i3startend(1,iproc),i3startend(3,jproc), is', i3startend(1,iproc),i3startend(3,jproc), is
              !!write(*,'(a,3i8)') 'i3startend(2,iproc),i3startend(4,jproc), ie', i3startend(2,iproc),i3startend(4,jproc), ie
              ! The min is for cases where a task has more than the entire box
              is=modulo(i3startend(1,iproc)-1,lzd%glr%d%n3+1)+1
              if (i3startend(2,iproc)-i3startend(1,iproc)>lzd%glr%d%n3) then 
                  ie=modulo(min(is-1,i3startend(2,iproc))-1,lzd%glr%d%n3+1)+1
              else
                  ie=modulo(i3startend(2,iproc)-1,lzd%glr%d%n3+1)+1
              end if
              js=i3startend(3,jproc)
              je=i3startend(4,jproc)
              if (je>=js) then
                  call get_extent_of_overlap(is, ie, js, je, n, ks, ke, nlen)
                  do k=1,n
                      ! Undo the periodic wrap around if required
                      !if (ks(k)>i3end) then
                      !    ii=ks(k)-(lzd%glr%d%n3+1)
                      !else
                      !    ii=ks(k)
                      !end if
                      ii=modulo(ks(k)-i3start-1,(lzd%glr%d%n3+1))+1
                      !write(*,'(a,7i9)') 'k, ks(k), ke(k), nlen(k), i3start, ks(k)-i3startend(3,jproc), ii', k, ks(k), ke(k), nlen(k), i3start, ks(k)-i3startend(3,jproc), ii
                      if (nproc>1) then
                          call mpiaccumulate(weightloc_f(0,0,ii), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), &
                               jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ks(k)-i3startend(3,jproc)),kind=mpi_address_kind), &
                               (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), mpi_sum, window_f)
                      else
                          call axpy((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), 1.d0, weightloc_f(0,0,ii), 1, &
                               weightppp_f(0,0,1+(ks(k)-i3startend(3,jproc))), 1)
                      end if
                      !!call mpiaccumulate(weightloc_c(0,0,ks(k)-modulo(i3start-1,lzd%glr%d%n1+1)+1), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), &
                      !!     jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ks(k)-i3startend(3,jproc)),kind=mpi_address_kind), &
                      !!     (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*nlen(k), mpi_sum, window_c)
                  end do
              end if
              !if (ie-is>=0) then
              !    !!call mpi_accumulate(weightloc_c(0,0,is-i3start), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), &
              !    !!     mpi_double_precision, jproc, &
              !    !!     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(is-i3startend(3,jproc)),kind=mpi_address_kind), &
              !    !!     (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), mpi_double_precision, &
              !    !!     mpi_sum, window_c, ierr)
              !    call mpiaccumulate(weightloc_c(0,0,is-i3start), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), &
              !         jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(is-i3startend(3,jproc)),kind=mpi_address_kind), &
              !         (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), mpi_sum, window_c)
              !end if
              !!!Check whether there is an overlap
              !!is = max(i3startend(1,iproc),i3startend(3,jproc))
              !!ie = min(i3startend(2,iproc),i3startend(4,jproc))
              !!if (ie-is>=0) then
              !!    !call mpi_accumulate(weightloc_f(0,0,is-i3start), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), &
              !!    !     mpi_double_precision, jproc, &
              !!    !     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(is-i3startend(3,jproc)),kind=mpi_address_kind), &
              !!    !     (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), mpi_double_precision, &
              !!    !     mpi_sum, window_f, ierr)
              !!    call mpiaccumulate(weightloc_f(0,0,is-i3start), (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), &
              !!         jproc, int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(is-i3startend(3,jproc)),kind=mpi_address_kind), &
              !!         (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), mpi_sum, window_f)
              !!end if
          end do
      end if
      !!else
      !!    !!is = i3startend(1,iproc)
      !!    !!ie = i3startend(2,iproc)
      !!    !!call vcopy((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(ie-is+1), weightloc_f(0,0,is-i3start), 1, &
      !!    !!     weightppp_f(0,0,is-i3startend(3,iproc)+1), 1)
      !!    call f_memcpy(n=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(j3end-j3start+1), src=weightloc_f(0,0,j3start), dest=weightppp_f(0,0,1))
      !!end if
      call f_free(i3startend)
      !call mpi_win_fence(0, window, ierr)
      !call mpi_win_free(window, ierr)
      !@END NEW #####################################

      call f_free(reducearr)
      !call f_free(weightloc)


      !weight_f_tot = 0.d0
      !do i3=1,n3p
      !    do i2=0,lzd%glr%d%n2
      !        do i1=0,lzd%glr%d%n1
      !            weightppp_f(i1,i2,i3)=weightppp_f(i1,i2,i3)**2
      !            weight_f_tot = weight_f_tot + weightppp_f(i1,i2,i3)
      !        end do
      !    end do
      !end do
      !if (nproc>1) then
      !    call mpiallred(weight_f_tot, 1, mpi_sum, bigdft_mpi%mpi_comm)
      !end if
      !write(*,*) 'iproc, weight_f_tot', iproc, weight_f_tot
      !@ENDNEW ##################################

      !!! Wait for the local completion of the mpi_raccumulate calls
      !!call mpiwait(request_c)
      !!call mpiwait(request_f)


      call f_release_routine()
    
    end subroutine get_weights


    subroutine assign_weight_to_process(iproc, nproc, lzd, i3s, n3p, window_c, window_f, &
               weight_c_tot_check, weight_f_tot_check, &
               weightppp_c, weightppp_f, weight_tot_c, weight_tot_f, &
               istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
               weightp_c, weightp_f, nptsp_c, nptsp_f, nvalp_c, nvalp_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, i3s, n3p
      type(local_zone_descriptors),intent(in) :: lzd
      integer,intent(inout) :: window_c, window_f
      real(kind=8),intent(in) :: weight_c_tot_check, weight_f_tot_check
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,1:max(1,n3p)),intent(inout) :: weightppp_c, weightppp_f
      real(kind=8),intent(out) :: weight_tot_c, weight_tot_f
      integer,dimension(2,0:nproc-1),intent(out) :: istartend_c, istartend_f
      integer,intent(out) :: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
      real(kind=8),intent(out) :: weightp_c, weightp_f
      integer,intent(out) :: nptsp_c, nptsp_f
      integer,intent(out) :: nvalp_c, nvalp_f
      
      ! Local variables
      integer :: jproc, i1, i2, i3, ii, istart, iend, j0, j1, ii_c, ii_f, n1p1, np, jjproc, jjjproc
      !!$$integer :: ii2, iiseg, jprocdone
      integer :: i, iseg, i0, iitot, ii3, ierr
      real(kind=8) :: tt, tt2, weight_c_ideal, weight_f_ideal, ttt, weight_prev
      real(kind=8),dimension(:,:),allocatable :: weights_c_startend, weights_f_startend
      character(len=*),parameter :: subname='assign_weight_to_process'
      integer,dimension(:),allocatable :: points_per_process, nval_c, nval_f
      integer,dimension(:,:),allocatable :: istartendseg_c
      real(kind=8),dimension(:),allocatable :: weightpp_c, weight_per_process_c
      integer,dimension(:,:),allocatable :: istartendseg_f
      real(kind=8),dimension(:),allocatable :: weightpp_f, weight_per_process_f

      call f_routine(id='assign_weight_to_process')

      ! Wait for the completion of the mpi_accumulate call started in get_weights
      if (nproc>1) then
          !!call mpi_win_fence(0, window_c, ierr)
          !!call mpi_win_free(window_c, ierr)
          call mpi_fenceandfree(window_c)
      end if

      tt=sum(weightppp_c)
      if (nproc>1) then
          call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (tt/=weight_c_tot_check) then
          write(*,'(a,2es20.10)') 'tt, weight_c_tot_check', tt, weight_c_tot_check
          stop 'tt/=weight_c_tot_check'
      end if

      !write(*,*) 'sum(weightppp_c)', sum(weightppp_c)
      !do i3=1,n3p
      !  do i2=0,lzd%glr%d%n2
      !    do i1=0,lzd%glr%d%n1
      !      write(1210,*) 'i1,i2,i3,val',i1,i2,i3,weightppp_c(i1,i2,i3)
      !    end do
      !  end do
      !end do

      weight_tot_c = 0.d0
      do i3=1,n3p
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  weightppp_c(i1,i2,i3)=weightppp_c(i1,i2,i3)**2
                  weight_tot_c = weight_tot_c + weightppp_c(i1,i2,i3)
              end do
          end do
      end do
      if (nproc>1) then
          call mpiallred(weight_tot_c, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
    
      weights_c_startend = f_malloc((/ 1.to.2, 0.to.nproc-1 /),id='weights_c_startend')
      weights_f_startend = f_malloc((/ 1.to.2, 0.to.nproc-1 /),id='weights_f_startend')

      ! Ideal weight per process.
      weight_c_ideal=weight_tot_c/dble(nproc)
    
      tt=0.d0
      weights_c_startend(1,0)=0.d0
      do jproc=0,nproc-2
          tt=tt+weight_c_ideal
          weights_c_startend(2,jproc)=dble(floor(tt,kind=8))
          weights_c_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
      end do
      weights_c_startend(2,nproc-1)=weight_tot_c
    
      ! Iterate through all grid points and assign them to processes such that the
      ! load balancing is optimal.


          !@NEW #################################
          !call to_zero(2*nproc, istartend_c(1,0))
          call f_zero(istartend_c)
          weight_per_process_c = f_malloc0(0.to.nproc-1,id='weight_per_process_c')
          points_per_process = f_malloc0(0.to.nproc-1,id='points_per_process')
          istartendseg_c = f_malloc0((/1.to.2,0.to.nproc-1/),id='istartendseg_c')
          !LG: why these parts are not allocated until nproc-1 only?
          nval_c = f_malloc0(0.to.nproc,id='nval_c')
          weightpp_c = f_malloc0(0.to.nproc,id='weightpp_c')

          weight_per_process_c(iproc) = sum(weightppp_c)
          if (nproc>1) then
             !this array was allocated with npproc-1 before
              call mpiallred(weight_per_process_c, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          if (sum(weight_per_process_c)/=weight_tot_c) then
              write(*,'(a,2f16.2)') 'sum(weight_per_process_c), weight_tot_c', sum(weight_per_process_c), weight_tot_c
              stop 'sum(weight_per_process_c)/=weight_tot_c'
          end if
          if (iproc==0) then
              weight_prev = 0.d0
              jjproc = 0
          else
              weight_prev = sum(weight_per_process_c(0:iproc-1)) !total weight of process up to iproc-1
              jjproc = nproc-1
              do jproc=0,nproc-1
                  !write(*,'(a,2i5,3f10.1)') 'iproc, jproc, weight_prev, (weights_c_startend(:,jproc))', iproc, jproc, weight_prev, (weights_c_startend(:,jproc))
                  !if (weight_prev<weights_c_startend(1,jproc)) then
                  if (weights_c_startend(1,jproc)<=weight_prev .and. weight_prev<=weights_c_startend(2,jproc)) then
                      ! This process starts the assignment with process jjproc
                      jjproc = jproc
                      exit
                  end if
                  !if (weight_prev+1.d0<=weights_c_startend(1,jproc) .and. &
                  !     weight_prev+weight_per_process_c(iproc)>=weights_c_startend(1,jproc)) then
                  !    jjproc=max(jproc-1,0)
                  !    exit
                  !end if
              end do
          end if

          ! Determine the number of grid points handled by each process
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)
          do iseg=1,lzd%glr%wfd%nseg_c
              j0=lzd%glr%wfd%keyglob(1,iseg)
              j1=lzd%glr%wfd%keyglob(2,iseg)
              ii=j0-1
              i3=ii/np
              if (i3+1<i3s) cycle
              if (i3+1>i3s+n3p-1) exit
              ii3=i3-i3s+1
              ii=ii-i3*np
              i2=ii/n1p1
              i0=ii-i2*n1p1
              i1=i0+j1-j0
              do i=i0,i1
                  points_per_process(iproc) = points_per_process(iproc) + 1
              end do
          end do
          if (nproc>1) then
              call mpiallred(points_per_process, mpi_sum,comm=bigdft_mpi%mpi_comm)
          end if


          !write(*,*) 'n3p, sum(weightppp_c(:,:,1:n3p))', n3p, sum(weightppp_c(:,:,1:n3p))

          tt = weight_prev
          ! number of gris points handled by processes 0..iproc-1
          iitot = sum(points_per_process(0:iproc-1)) !total number of grid points up to iproc-1
          !!write(*,'(a,i7,f14.1,2i9)') 'start: iproc, tt, iitot, jjproc', iproc, tt, iitot, jjproc
          !!write(*,'(a,i5,100f12.1)') 'iproc, weights_c_startend', iproc, weights_c_startend
          ! only do this on task 0 due to the allreduce later
          if (iproc==0) then
              istartend_c(1,0) = 1
              istartendseg_c(1,0) = 1
          end if
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)
          jjjproc = 0
          do iseg=1,lzd%glr%wfd%nseg_c
              j0=lzd%glr%wfd%keyglob(1,iseg)
              j1=lzd%glr%wfd%keyglob(2,iseg)
              ii=j0-1
              i3=ii/np
              if (i3+1<i3s) cycle
              if (i3+1>i3s+n3p-1) exit
              ii3=i3+1-i3s+1
              ii=ii-i3*np
              i2=ii/n1p1
              i0=ii-i2*n1p1
              i1=i0+j1-j0
              do i=i0,i1
                  iitot = iitot + 1
                  tt = tt + weightppp_c(i,i2,ii3)
                  if (jjproc<nproc-1) then
                      if (tt>=weights_c_startend(1,jjproc+1)) then
                          !write(*,'(a,2i6,2f10.1)') 'iproc, jjproc, tt, weights_c_startend(1,jjproc+1)', iproc, jjproc, tt, weights_c_startend(1,jjproc+1)
                          jjproc = jjproc + 1
                          jjjproc = jjproc
                          istartend_c(1,jjproc) = iitot
                          istartendseg_c(1,jjproc) = iseg
                      end if
                  end if
                  if (weightppp_c(i,i2,ii3)>0.d0) then
                      nval_c(jjproc) = nval_c(jjproc) + nint(sqrt(weightppp_c(i,i2,ii3))) !total number of grid points to be handled by process jjproc
                      weightpp_c(jjproc) = weightpp_c(jjproc) + weightppp_c(i,i2,ii3) !total weight to be handled by process jjproc
                      !write(500,'(a,3i7,f12.2)') 'i, i2, ii3, weightpp_c(jjproc)',i, i2, ii3, weightpp_c(jjproc)
                      !!weightppp_c(i,i2,ii3)=0.d0
                  end if
              end do
          end do
          !write(*,*) 'AFTER SET TO ZERO: sum(weightppp_c)',sum(weightppp_c)
          !!do i3=1,n3p
          !!    write(*,*) 'i3, sum(weightppp_c(:,:,i3))', i3, sum(weightppp_c(:,:,i3))
          !!end do

          ! Communicate the data and assign the processor specific values
          if (nproc>1) then
              call mpiallred(istartend_c(1,0), 2*nproc, mpi_sum, comm=bigdft_mpi%mpi_comm) !a bit wasteful to communicate the zeros of the second entry...
              call mpiallred(istartendseg_c, mpi_sum, comm=bigdft_mpi%mpi_comm) !a bit wasteful to communicate the zeros of the second entry...
              call mpiallred(nval_c(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(weightpp_c(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(jjjproc, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
          end if
          ! jjjproc is the last task which has been assigned
          !do jproc=0,nproc-2
          do jproc=0,jjjproc-1
              istartend_c(2,jproc) = istartend_c(1,jproc+1)-1
              istartendseg_c(2,jproc) = istartendseg_c(1,jproc+1)
          end do
          ! Take the rest
          !!istartend_c(2,nproc-1) = lzd%glr%wfd%nvctr_c
          !!istartendseg_c(2,nproc-1) = lzd%glr%wfd%nseg_c
          istartend_c(2,jjjproc) = lzd%glr%wfd%nvctr_c
          istartendseg_c(2,jjjproc) = lzd%glr%wfd%nseg_c
          ! Fill with "empty" values
          do jproc=jjjproc+1,nproc-1
              istartend_c(1,jproc) = lzd%glr%wfd%nvctr_c + 1
              istartend_c(2,jproc) = lzd%glr%wfd%nvctr_c
              istartendseg_c(1,jproc) = lzd%glr%wfd%nseg_c + 1
              istartendseg_c(2,jproc) = lzd%glr%wfd%nseg_c
          end do
          istartp_seg_c = istartendseg_c(1,iproc)
          iendp_seg_c = istartendseg_c(2,iproc)
          nvalp_c = nval_c(iproc)
          weightp_c = weightpp_c(iproc)
          nptsp_c=istartend_c(2,iproc)-istartend_c(1,iproc)+1

          call f_free(weight_per_process_c)
          call f_free(points_per_process)
          call f_free(istartendseg_c)
          call f_free(nval_c)
          call f_free(weightpp_c)
          !write(*,'(a,i7,100i12)') 'new: iproc, istartend_c',iproc, istartend_c 
          !!write(*,'(a,i7,100i12)') 'new: iproc, istartp_seg_c', iproc, istartp_seg_c
          !!write(*,'(a,i7,100i12)') 'new: iproc, iendp_seg_c', iproc, iendp_seg_c
          !!write(*,'(a,i7,100i12)') 'new: iproc, nvalp_c', iproc, nvalp_c
          !!write(*,'(a,i7,100f12.1)') 'new: iproc, weightp_c', iproc, weightp_c
          !!write(*,'(a,i7,100i12)') 'new: iproc, nptsp_c', iproc, nptsp_c

          ! Some checks
          ii_c=istartend_c(2,iproc)-istartend_c(1,iproc)+1
          if (nproc > 1) then
            call mpiallred(ii_c, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          if(ii_c/=lzd%glr%wfd%nvctr_c) then
             write(*,*) 'ii_c/=lzd%glr%wfd%nvctr_c',ii_c,lzd%glr%wfd%nvctr_c
             stop
          end if
    
          if (nproc > 1) then
             call mpiallred(weightp_c,1,mpi_sum,comm=bigdft_mpi%mpi_comm,recvbuf=tt)
          else
              tt=weightp_c
          end if
          if(tt/=weight_tot_c) then
              write(*,*) 'tt, weight_tot_c', tt, weight_tot_c
              stop 'wrong partition of coarse weights'
          end if

          if (nproc > 1) then
             call mpiallred(nptsp_c, 1,mpi_sum, comm=bigdft_mpi%mpi_comm,recvbuf=ii)
          else
              ii=nptsp_c
          end if
          if(ii/=lzd%glr%wfd%nvctr_c) then
              write(*,*) 'ii, lzd%glr%wfd%nvctr_c', ii, lzd%glr%wfd%nvctr_c
              stop 'wrong partition of coarse grid points'
          end if
          !@END NEW #############################


      !!end if
    
      ! Same for fine region

      ! Wait for the completion of the mpi_accumulate call started in get_weights
      if (nproc>1) then
          !!call mpi_win_fence(0, window_f, ierr)
          !!call mpi_win_free(window_f, ierr)
          call mpi_fenceandfree(window_f)
      end if


      tt=sum(weightppp_f)
      if (nproc>1) then
          call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (tt/=weight_f_tot_check) then
          write(*,'(a,2es20.10)') 'tt, weight_f_tot_check', tt, weight_f_tot_check
          stop 'tt/=weight_f_tot_check'
      end if

      !do i3=1,n3p
      !  do i2=0,lzd%glr%d%n2
      !    do i1=0,lzd%glr%d%n1
      !      write(1211,*) 'i1,i2,i3,val',i1,i2,i3,weightppp_f(i1,i2,i3)
      !    end do
      !  end do
      !end do

      weight_tot_f = 0.d0
      do i3=1,n3p
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  weightppp_f(i1,i2,i3)=weightppp_f(i1,i2,i3)**2
                  weight_tot_f = weight_tot_f + weightppp_f(i1,i2,i3)
              end do
          end do
      end do
      if (nproc>1) then
          call mpiallred(weight_tot_f, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      ! Ideal weight per process.
      weight_f_ideal=weight_tot_f/dble(nproc)

      tt=0.d0
      weights_f_startend(1,0)=0.d0
      do jproc=0,nproc-2
          tt=tt+weight_f_ideal
          weights_f_startend(2,jproc)=dble(floor(tt,kind=8))
          weights_f_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
      end do
      weights_f_startend(2,nproc-1)=weight_tot_f
    



          !@NEW #################################
      !call to_zero(2*nproc, istartend_f(1,0))
          call f_zero(istartend_f)
          weight_per_process_f = f_malloc0(0.to.nproc,id='weight_per_process_f')
          points_per_process = f_malloc0(0.to.nproc-1,id='points_per_process')
          istartendseg_f = f_malloc0((/1.to.2,0.to.nproc-1/),id='istartendseg_f')
          nval_f = f_malloc0(0.to.nproc,id='nval_f')
          weightpp_f = f_malloc0(0.to.nproc,id='weightpp_f')

          weight_per_process_f(iproc) = sum(weightppp_f)
          if (nproc>1) then
              call mpiallred(weight_per_process_f(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          if (sum(weight_per_process_f)/=weight_tot_f) then
              write(*,'(a,2f16.2)') 'sum(weight_per_process_f), weight_tot_f', sum(weight_per_process_f), weight_tot_f
              stop 'sum(weight_per_process_f)/=weight_tot_f'
          end if
          if (iproc==0) then
              weight_prev = 0.d0
              jjproc = 0
          else
              weight_prev = sum(weight_per_process_f(0:iproc-1)) !total weight of process up to iproc-1
              jjproc = nproc-1
              do jproc=0,nproc-1
                  !write(*,'(a,2i5,3f10.1)') 'iproc, jproc, weight_prev, (weights_f_startend(:,jproc))', iproc, jproc, weight_prev, (weights_f_startend(:,jproc))
                  if (weights_f_startend(1,jproc)<=weight_prev .and. weight_prev<=weights_f_startend(2,jproc)) then
                      ! This process starts the assignment with process jjproc
                      jjproc = jproc
                      exit
                  end if
                  !if (weight_prev+1.d0<=weights_f_startend(1,jproc) .and. &
                  !    weight_prev+weight_per_process_f(iproc)>=weights_f_startend(1,jproc)) then
                  !    jjproc=max(jproc-1,0)
                  !    exit
                  !end if
              end do
          end if

          ! Determine the number of grid points handled by each process
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)
          istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
          iend=istart+lzd%glr%wfd%nseg_f-1
          if (istart<=iend) then
              do iseg=istart,iend
                  j0=lzd%glr%wfd%keyglob(1,iseg)
                  j1=lzd%glr%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/np
                  if (i3+1<i3s) cycle
                  if (i3+1>i3s+n3p-1) exit
                  ii3=i3-i3s+1
                  ii=ii-i3*np
                  i2=ii/n1p1
                  i0=ii-i2*n1p1
                  i1=i0+j1-j0
                  do i=i0,i1
                      points_per_process(iproc) = points_per_process(iproc) + 1
                  end do
              end do
          end if
          if (nproc>1) then
              call mpiallred(points_per_process(0), nproc, mpi_sum,comm=bigdft_mpi%mpi_comm)
          end if



          tt = weight_prev
          ! number of gris points handled by processes 0..iproc-1
          iitot = sum(points_per_process(0:iproc-1)) !total number of grid points up to iproc-1
          !!write(*,'(a,i7,f14.1,2i9)') 'start: iproc, tt, iitot, jjproc', iproc, tt, iitot, jjproc
          !!write(*,'(a,i5,100f12.1)') 'iproc, weights_f_startend', iproc, weights_f_startend
          ! only do this on task 0 due to the allreduce later
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)
          istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
          iend=istart+lzd%glr%wfd%nseg_f-1
          if (iproc==0) then
              istartend_f(1,0) = 1
              istartendseg_f(1,0) = istart
          end if
          jjjproc = 0
          if (istart<=iend) then
              do iseg=istart,iend
                  j0=lzd%glr%wfd%keyglob(1,iseg)
                  j1=lzd%glr%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/np
                  if (i3+1<i3s) cycle
                  if (i3+1>i3s+n3p-1) exit
                  ii3=i3+1-i3s+1
                  ii=ii-i3*np
                  i2=ii/n1p1
                  i0=ii-i2*n1p1
                  i1=i0+j1-j0
                  do i=i0,i1
                      iitot = iitot + 1
                      tt = tt + weightppp_f(i,i2,ii3)
                      if (jjproc<nproc-1) then
                          if (tt>=weights_f_startend(1,jjproc+1)) then
                              !write(*,'(a,2i6,2f10.1)') 'iproc, jjproc, tt, weights_f_startend(1,jjproc+1)', iproc, jjproc, tt, weights_f_startend(1,jjproc+1)
                              jjproc = jjproc + 1
                              jjjproc = jjproc
                              istartend_f(1,jjproc) = iitot
                              istartendseg_f(1,jjproc) = iseg
                          end if
                      end if
                      if (weightppp_f(i,i2,ii3)>0.d0) then
                          nval_f(jjproc) = nval_f(jjproc) + nint(sqrt(weightppp_f(i,i2,ii3))) !total number of grid points to be handled by process jjproc
                          weightpp_f(jjproc) = weightpp_f(jjproc) + weightppp_f(i,i2,ii3) !total weight to be handled by process jjproc
                      end if
                  end do
              end do
          !!!!end if

          ! Communicate the data and assign the processor specific values
          if (nproc>1) then
              call mpiallred(istartend_f(1,0), 2*nproc, mpi_sum,comm= bigdft_mpi%mpi_comm) !a bit wasteful to communicate the zeros of the second entry...
              call mpiallred(istartendseg_f(1,0), 2*nproc, mpi_sum, comm=bigdft_mpi%mpi_comm) !a bit wasteful to communicate the zeros of the second entry...
              call mpiallred(nval_f(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(weightpp_f(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(jjjproc, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
          end if
          ! jjjproc is the last task which has been assigned
          !!do jproc=0,nproc-2
          do jproc=0,jjjproc-1
              istartend_f(2,jproc) = istartend_f(1,jproc+1)-1
              istartendseg_f(2,jproc) = istartendseg_f(1,jproc+1)
          end do
          ! Take the rest
          !!istartend_f(2,nproc-1) = lzd%glr%wfd%nvctr_f
          !!istartendseg_f(2,nproc-1) = lzd%glr%wfd%nseg_c + lzd%glr%wfd%nseg_f
          istartend_f(2,jjjproc) = lzd%glr%wfd%nvctr_f
          istartendseg_f(2,jjjproc) = lzd%glr%wfd%nseg_c + lzd%glr%wfd%nseg_f
          ! Fill with "empty" values
          do jproc=jjjproc+1,nproc-1
              istartend_f(1,jproc) = lzd%glr%wfd%nvctr_f + 1
              istartend_f(2,jproc) = lzd%glr%wfd%nvctr_f
              istartendseg_f(1,jproc) = lzd%glr%wfd%nseg_c + lzd%glr%wfd%nseg_f + 1
              istartendseg_f(2,jproc) = lzd%glr%wfd%nseg_c + lzd%glr%wfd%nseg_f
          end do
          istartp_seg_f = istartendseg_f(1,iproc)
          iendp_seg_f = istartendseg_f(2,iproc)
          nvalp_f = nval_f(iproc)
          weightp_f = weightpp_f(iproc)
          nptsp_f=istartend_f(2,iproc)-istartend_f(1,iproc)+1

          call f_free(weight_per_process_f)
          call f_free(points_per_process)
          call f_free(istartendseg_f)
          call f_free(nval_f)
          call f_free(weightpp_f)
          !write(*,'(a,i7,100i12)') 'new: iproc, istartend_f',iproc, istartend_f 
          !!write(*,'(a,i7,100i12)') 'new: iproc, istartp_seg_f', iproc, istartp_seg_f
          !!write(*,'(a,i7,100i12)') 'new: iproc, iendp_seg_f', iproc, iendp_seg_f
          !write(*,'(a,i7,100i12)') 'new: iproc, nvalp_f', iproc, nvalp_f
          !!write(*,'(a,i7,100f12.1)') 'new: iproc, weightp_f', iproc, weightp_f
          !!write(*,'(a,i7,100i12)') 'new: iproc, nptsp_f', iproc, nptsp_f

          ! Some checks
          ii_f=istartend_f(2,iproc)-istartend_f(1,iproc)+1
          if (nproc > 1) then
            call mpiallred(ii_f, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          if(ii_f/=lzd%glr%wfd%nvctr_f) then
             write(*,*) 'ii_f/=lzd%glr%wfd%nvctr_f',ii_f,lzd%glr%wfd%nvctr_f
             stop
          end if
    
          if (nproc > 1) then
             call mpiallred(weightp_f,1,mpi_sum, comm=bigdft_mpi%mpi_comm,recvbuf=tt)
          else
              tt=weightp_f
          end if
          if(tt/=weight_tot_f) then
              write(*,*) 'tt, weight_tot_f', tt, weight_tot_f
              stop 'wrong partition of fine weights'
          end if

          if (nproc > 1) then
             call mpiallred(nptsp_f, 1,mpi_sum, comm=bigdft_mpi%mpi_comm,recvbuf=ii)
          else
              ii=nptsp_f
          end if
          if(ii/=lzd%glr%wfd%nvctr_f) then
              write(*,*) 'ii, lzd%glr%wfd%nvctr_f', ii, lzd%glr%wfd%nvctr_f
              stop 'wrong partition of coarse grid points'
          end if
          !@END NEW #############################



      end if
    
    
    
    
      call f_free(weights_c_startend)
      call f_free(weights_f_startend)
    

      call f_release_routine()
      
    end subroutine assign_weight_to_process

!!    subroutine assign_weight_to_process_new(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
!!               istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
!!               weightp_c, weightp_f, nptsp_c, nptsp_f, nvalp_c, nvalp_f)
!!      use module_base
!!      use module_types
!!      implicit none
!!      
!!      ! Calling arguments
!!      integer,intent(in) :: iproc, nproc
!!      type(local_zone_descriptors),intent(in) :: lzd
!!      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: weight_c, weight_f
!!      real(kind=8),intent(in) :: weight_tot_c, weight_tot_f
!!      integer,dimension(2,0:nproc-1),intent(out) :: istartend_c, istartend_f
!!      integer,intent(out) :: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
!!      real(kind=8),intent(out) :: weightp_c, weightp_f
!!      integer,intent(out) :: nptsp_c, nptsp_f
!!      integer,intent(out) :: nvalp_c, nvalp_f
!!      
!!      ! Local variables
!!      integer :: jproc, i1, i2, i3, ii, istart, iend, j0, j1, ii_c, ii_f, n1p1, np, it, npr
!!      integer :: i, iseg, i0, iitot, eproc, sproc, ith, nth, ierr, eseg, iitotseg, iitote
!!      real(kind=8) :: tt, ttseg, tt2, weight_c_ideal, weight_f_ideal, nproc_block, ttt
!!      real(kind=8),dimension(:,:),allocatable :: weights_c_startend, weights_f_startend
!!      real(kind=4) :: tr0, tr1
!!      real(kind=8) :: time1, time2
!!      !$ integer  :: omp_get_thread_num,omp_get_max_threads
!!    
!!      ! Ideal weight per process.
!!      weight_c_ideal=weight_tot_c/dble(nproc)
!!      weight_f_ideal=weight_tot_f/dble(nproc)
!!    
!!      weights_c_startend = f_malloc((/ 1.to.2, 0.to.nproc-1 /),id='weights_c_startend')
!!      weights_f_startend = f_malloc((/ 1.to.2, 0.to.nproc-1 /),id='weights_f_startend')
!!    
!!      tt=0.d0
!!      weights_c_startend(1,0)=0.d0
!!      do jproc=0,nproc-2
!!          tt=tt+weight_c_ideal
!!          weights_c_startend(2,jproc)=dble(floor(tt,kind=8))
!!          weights_c_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
!!      end do
!!      weights_c_startend(2,nproc-1)=weight_tot_c
!!
!!      !split into subroutine and divide into larger sections so each thread does one chunk of MPI procs
!!      istart=1
!!      iend=lzd%glr%wfd%nseg_c
!!
!!      if (nproc==1) then
!!         istartend_c(1,0)=1
!!         istartend_c(2,0)=lzd%glr%wfd%nvctr_c
!!         weightp_c = weight_tot_c 
!!         istartp_seg_c=istart
!!         iendp_seg_c=iend
!!         ttt=0.d0
!!         do i1=0,lzd%glr%d%n1
!!            do i2=0,lzd%glr%d%n2
!!               do i3=0,lzd%glr%d%n3
!!                  ttt = ttt+sqrt(weight_c(i1,i2,i3))
!!               end do
!!            end do
!!         end do
!!         nvalp_c=nint(ttt)
!!      else
!!         nth=1
!!         !$  nth = OMP_GET_max_threads()
!!         nproc_block=real(nproc,kind=8)/nth
!!
!!         !$omp parallel default(none) &
!!         !$omp private(sproc,eproc,iitot,iitote,iitotseg,ttseg,eseg,ith,time1,tr0,tr1,time2) &
!!         !$omp shared(nproc_block,istart,iend,lzd,weight_c,iproc,nproc,weights_c_startend,weight_tot_c) &
!!         !$omp shared(istartp_seg_c,iendp_seg_c,weightp_c,istartend_c,nth,nvalp_c)
!!         ith=0
!!         !$ ith = OMP_GET_THREAD_NUM()
!!         !check we have enough MPI tasks for each thread
!!         if (nproc_block>1) then
!!            sproc=nint(nproc_block*ith)
!!            if (ith==nth-1) then
!!               eproc=nproc-1
!!            else
!!               eproc=nint(nproc_block*(ith+1))-1
!!            end if
!!         else
!!            if (ith<nproc) then
!!               sproc=ith
!!               eproc=ith
!!            else
!!               sproc=-1
!!               eproc=-1
!!            end if
!!         end if
!!
!!         !call cpu_time(tr0)
!!         if (ith/=0.and.sproc/=-1) then
!!            call assign_weight_to_process_find_end_point(nproc, sproc, lzd, weight_c, &
!!                 weights_c_startend, istart, iend, ttseg, iitotseg, eseg)
!!         else
!!            ttseg=0.d0
!!            iitotseg=0
!!            eseg=istart
!!         end if
!!         !call cpu_time(tr1)
!!         !time2=real(tr1-tr0,kind=8)
!!
!!         if (sproc/=-1) call assign_weight_to_process_sub(iproc, eproc+1, lzd, weight_c, weight_tot_c, &
!!              istartend_c(1,sproc), istartp_seg_c, iendp_seg_c, weightp_c, weights_c_startend(1,sproc), &
!!              eseg, iend, ttseg, iitotseg, sproc, nvalp_c)
!!         !call cpu_time(tr0)
!!         !time1=real(tr0-tr1,kind=8)
!!         !if (iproc==0) print*,'thread times',iproc,ith,time2,time1,time1+time2
!!         !$omp end parallel
!!
!!         ! check
!!         !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!         !if (iproc==0) print*,''
!!         !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!         !do jproc=0,nproc-1!sproc,eproc
!!         !   if (iproc==jproc) then
!!         !      print*,istartend_c(1,jproc),istartend_c(2,jproc),nint(weightp_c),nint(weight_tot_c/dble(nproc)),nint(weightp_c-weight_tot_c/dble(nproc))
!!         !   end if
!!         !   call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!         !end do
!!      end if
!!
!!      ! Same for fine region
!!      tt=0.d0
!!      weights_f_startend(1,0)=0.d0
!!      do jproc=0,nproc-2
!!          tt=tt+weight_f_ideal
!!          weights_f_startend(2,jproc)=dble(floor(tt,kind=8))
!!          weights_f_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
!!      end do
!!      weights_f_startend(2,nproc-1)=weight_tot_f
!!    
!!      istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
!!      iend=istart+lzd%glr%wfd%nseg_f-1
!!
!!      if (nproc==1) then
!!         istartend_f(1,0)=1
!!         istartend_f(2,0)=lzd%glr%wfd%nvctr_f
!!         weightp_f = weight_tot_f
!!         istartp_seg_f=istart
!!         iendp_seg_f=iend
!!         ttt=0.d0
!!         do i1=0,lzd%glr%d%n1
!!            do i2=0,lzd%glr%d%n2
!!               do i3=0,lzd%glr%d%n3
!!                  ttt = ttt+sqrt(weight_f(i1,i2,i3))
!!               end do
!!            end do
!!         end do
!!         nvalp_f=nint(ttt)
!!      else
!!
!!         !$omp parallel default(none) &
!!         !$omp private(sproc,eproc,iitot,iitote,iitotseg,ttseg,eseg,ith) &
!!         !$omp shared(nproc_block,istart,iend,lzd,weight_f,iproc,nproc,weights_f_startend,weight_tot_f) &
!!         !$omp shared(istartp_seg_f,iendp_seg_f,weightp_f,istartend_f,nth,nvalp_f)
!!         ith=0
!!         !$ ith = OMP_GET_THREAD_NUM()
!!         !check we have enough MPI tasks for each thread
!!         if (nproc_block>1) then
!!            sproc=nint(nproc_block*ith)
!!            if (ith==nth-1) then
!!               eproc=nproc-1
!!            else
!!               eproc=nint(nproc_block*(ith+1))-1
!!            end if
!!         else
!!            if (ith<nproc) then
!!               sproc=ith
!!               eproc=ith
!!            else
!!               sproc=-1
!!               eproc=-1
!!            end if
!!         end if
!!
!!         if (ith/=0.and.sproc/=-1) then
!!            call assign_weight_to_process_find_end_point(nproc, sproc, lzd, weight_f, &
!!                 weights_f_startend, istart, iend, ttseg, iitotseg, eseg)
!!         else
!!            ttseg=0.d0
!!            iitotseg=0
!!            eseg=istart
!!         end if
!!
!!         if (sproc/=-1)call assign_weight_to_process_sub(iproc, eproc+1, lzd, weight_f, weight_tot_f, &
!!              istartend_f(1,sproc), istartp_seg_f, iendp_seg_f, weightp_f, weights_f_startend(1,sproc), &
!!              eseg, iend, ttseg, iitotseg, sproc, nvalp_f)
!!         !$omp end parallel
!!
!!         ! check
!!         !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!         !if (iproc==0) print*,''
!!         !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!         !do jproc=0,nproc-1!sproc,eproc
!!         !   if (iproc==jproc) then
!!         !      print*,istartend_f(1,jproc),istartend_f(2,jproc),nint(weightp_f),nint(weight_tot_f/dble(nproc)),nint(weightp_f-weight_tot_f/dble(nproc))
!!         !   end if
!!         !   call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!         !end do
!!      end if
!!
!!      call f_free(weights_c_startend)
!!      call f_free(weights_f_startend)
!!    
!!      nptsp_c=istartend_c(2,iproc)-istartend_c(1,iproc)+1
!!      nptsp_f=istartend_f(2,iproc)-istartend_f(1,iproc)+1
!!        
!!      ! some check
!!      ii_f=istartend_f(2,iproc)-istartend_f(1,iproc)+1
!!      if (nproc > 1) then
!!        call mpiallred(ii_f, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!      end if
!!      !if(ii_f/=lzd%glr%wfd%nvctr_f) stop 'assign_weight_to_process: ii_f/=lzd%glr%wfd%nvctr_f'
!!      if(ii_f/=lzd%glr%wfd%nvctr_f) then
!!         write(*,*) 'ii_f/=lzd%glr%wfd%nvctr_f',ii_f,lzd%glr%wfd%nvctr_f
!!         if (iproc==0) then
!!             do jproc=0,nproc-1
!!                 write(*,*) jproc, istartend_f(1,jproc), istartend_f(2,jproc)
!!             end do
!!         end if
!!         stop
!!      end if
!!     
!!      ii_c=istartend_c(2,iproc)-istartend_c(1,iproc)+1
!!      if (nproc > 1) then
!!        call mpiallred(ii_c, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!      end if
!!      if(ii_c/=lzd%glr%wfd%nvctr_c) then
!!         write(*,*) 'ii_c/=lzd%glr%wfd%nvctr_c',ii_c,lzd%glr%wfd%nvctr_c
!!         stop
!!      end if
!!    
!!      ! some checks
!!      if (nproc > 1) then
!!         call mpiallred(weightp_c,1,mpi_sum, comm=bigdft_mpi%mpi_comm,recvbuf=tt)
!!         !call mpi_allreduce(weightp_c, tt, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!      else
!!          tt=weightp_c
!!      end if
!!      if(tt/=weight_tot_c) then
!!         write(*,*) 'wrong partition of coarse weights',tt,weight_tot_c
!!         stop
!!      end if
!!      if (nproc > 1) then
!!         call mpiallred(weightp_f,1,mpi_sum,comm=bigdft_mpi%mpi_comm,recvbuf=tt)
!!         !call mpi_allreduce(weightp_f, tt, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!      else
!!          tt=weightp_f
!!      end if     
!!      if(tt/=weight_tot_f) then
!!         write(*,*) 'wrong partition of fine weights',tt,weight_tot_f
!!         stop
!!      end if
!!      if (nproc > 1) then
!!         call mpiallred(nptsp_c, 1,mpi_sum, comm=bigdft_mpi%mpi_comm,recvbuf=ii)
!!         !call mpi_allreduce(nptsp_c, ii, 1, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!      else
!!          ii=nptsp_c
!!      end if
!!      if(ii/=lzd%glr%wfd%nvctr_c) then
!!         write(*,*) 'wrong partition of coarse grid points',ii,lzd%glr%wfd%nvctr_c
!!         stop
!!      end if
!!      if (nproc > 1) then
!!         call mpiallred(nptsp_f, 1,mpi_sum,comm= bigdft_mpi%mpi_comm,recvbuf=ii)
!!         !call mpi_allreduce(nptsp_f, ii, 1, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!      else
!!          ii=nptsp_f
!!      end if
!!      if(ii/=lzd%glr%wfd%nvctr_f) then
!!         write(*,*) 'wrong partition of fine grid points',ii,lzd%glr%wfd%nvctr_f
!!         stop
!!      end if
!!      
!!    end subroutine assign_weight_to_process_new
!!
!!    !better name and change names to indicate coarse OR fine
!!    subroutine assign_weight_to_process_sub(iproc, nproc, lzd, weight_c, weight_tot_c, &
!!               istartend_c, istartp_seg_c, iendp_seg_c, weightp_c, weights_c_startend,&
!!               istart, iend, ttseg, iitotseg, jprocs, nvalp_c)
!!      use module_base
!!      use module_types
!!      implicit none
!!      
!!      ! Calling arguments
!!      integer,intent(in) :: iproc, nproc
!!      type(local_zone_descriptors),intent(in) :: lzd
!!      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: weight_c
!!      real(kind=8),intent(in) :: weight_tot_c
!!      integer,dimension(2,jprocs:nproc-1),intent(out) :: istartend_c
!!      integer,intent(inout) :: istartp_seg_c, iendp_seg_c
!!      real(kind=8),intent(inout) :: weightp_c
!!      real(kind=8),dimension(1:2,jprocs:nproc-1),intent(in) :: weights_c_startend
!!      integer, intent(in) :: istart, iend
!!      integer, intent(in) :: iitotseg
!!      real(kind=8), intent(in) :: ttseg
!!      integer, intent(in) :: jprocs
!!      integer, intent(out) :: nvalp_c
!!      
!!      ! Local variables
!!      integer ::  jproc, i1, i2, i3, ii, j0, j1, n1p1, np, i, iseg, i0, iitot
!!      real(kind=8) :: tt2, tt, ttt
!!    
!!      ! Iterate through all grid points and assign them to processes such that the
!!      ! load balancing is optimal.
!!      jproc=jprocs
!!      tt=ttseg
!!      tt2=0.d0
!!      ttt=0.d0
!!      iitot=iitotseg
!!      n1p1=lzd%glr%d%n1+1
!!      np=n1p1*(lzd%glr%d%n2+1)
!!      loop_nseg_c: do iseg=istart,iend
!!         j0=lzd%glr%wfd%keyglob(1,iseg)
!!         j1=lzd%glr%wfd%keyglob(2,iseg)
!!         ii=j0-1
!!         i3=ii/np
!!         ii=ii-i3*np
!!         i2=ii/n1p1
!!         i0=ii-i2*n1p1
!!         i1=i0+j1-j0
!!         do i=i0,i1
!!            tt=tt+weight_c(i,i2,i3)
!!            tt2=tt2+weight_c(i,i2,i3)
!!            ttt=ttt+sqrt(weight_c(i,i2,i3))
!!            iitot=iitot+1
!!            if (jproc==nproc) then
!!               if (tt>weights_c_startend(2,jproc-1)+1) exit loop_nseg_c
!!            else if (tt>weights_c_startend(1,jproc)) then
!!               if (jproc>jprocs) then
!!                  if (iproc==jproc) then
!!                     istartp_seg_c=iseg
!!                  else if (iproc==jproc-1) then
!!                     iendp_seg_c=iseg
!!                     weightp_c=tt2
!!                     nvalp_c=nint(ttt)
!!                  end if
!!                  tt2=0.d0
!!                  ttt=0.d0
!!                  istartend_c(1,jproc)=iitot+1
!!               else if (jproc==jprocs) then
!!                  if (jproc==0) then
!!                     istartend_c(1,jproc)=1
!!                     if (iproc==jproc) istartp_seg_c=istart
!!                  else
!!                     tt2=0.d0
!!                     ttt=0.d0
!!                     istartend_c(1,jproc)=iitot+1
!!                     if (iproc==jproc) istartp_seg_c=iseg
!!                  end if
!!               end if
!!               jproc=jproc+1
!!            end if
!!         end do
!!      end do loop_nseg_c
!!    
!!      do jproc=jprocs,nproc-2
!!         istartend_c(2,jproc)=istartend_c(1,jproc+1)-1
!!      end do
!!      istartend_c(2,nproc-1)=iitot    
!!
!!      if(iproc==nproc-1) then
!!         weightp_c=tt2
!!         iendp_seg_c=min(iseg,iend)
!!         nvalp_c=nint(ttt)
!!      end if
!!
!!    end subroutine assign_weight_to_process_sub
!!
!!    subroutine assign_weight_to_process_find_end_point(iproc, nproc, lzd, weight_c, &
!!               weights_c_startend, istart, iend, ttseg, iitot, iseg)
!!      use module_base
!!      use module_types
!!      implicit none
!!      
!!      ! Calling arguments
!!      integer,intent(in) :: iproc, nproc !technically nproc isn't full nproc...
!!      type(local_zone_descriptors),intent(in) :: lzd
!!      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: weight_c
!!      real(kind=8),dimension(1:2,0:nproc),intent(in) :: weights_c_startend
!!      integer, intent(in) :: istart, iend
!!      integer, intent(out) :: iitot, iseg
!!      real(kind=8), intent(out) :: ttseg
!!      
!!      ! Local variables
!!      integer ::  i1, i2, i3, ii, j0, j1, n1p1, np, i, i0
!!      real(kind=8) :: tt
!!
!!      tt=0.d0
!!      iitot=0
!!      n1p1=lzd%glr%d%n1+1
!!      np=n1p1*(lzd%glr%d%n2+1)
!!      loop_nseg_c: do iseg=istart,iend
!!         j0=lzd%glr%wfd%keyglob(1,iseg)
!!         j1=lzd%glr%wfd%keyglob(2,iseg)
!!         ii=j0-1
!!         i3=ii/np
!!         ii=ii-i3*np
!!         i2=ii/n1p1
!!         i0=ii-i2*n1p1
!!         i1=i0+j1-j0
!!         tt=tt+sum(weight_c(i0:i1,i2,i3))
!!         if (tt>weights_c_startend(1,nproc)) exit
!!         ttseg=tt
!!         iitot=iitot+i1-i0+1
!!      end do loop_nseg_c
!!
!!    end subroutine assign_weight_to_process_find_end_point

    subroutine get_index_in_global2(lr, ii3min, ii3max, jj3min, index_in_global_c, index_in_global_f)
    use module_base
    use module_types
    implicit none
    
    ! Calling arguments
    type(locreg_descriptors),intent(in) :: lr
    integer,intent(in) :: ii3min, ii3max, jj3min
    integer,dimension(0:lr%d%n1,0:lr%d%n2,ii3min:ii3max),intent(out) :: index_in_global_c, index_in_global_f
    
    ! Local variables
    integer :: iitot, iseg, j0, j1, ii, i1, i2, i3, i0, i, istart, iend, np, n1p1, jj3
    
    call f_routine(id='get_index_in_global2')

    ! Could optimize these loops by cycling and updating iitot as soon as
    ! (i3<ii3min .or. i3>ii3max)
    
    iitot=0
    n1p1=lr%d%n1+1
    np=n1p1*(lr%d%n2+1)
    do iseg=1,lr%wfd%nseg_c
       j0=lr%wfd%keyglob(1,iseg)
       j1=lr%wfd%keyglob(2,iseg)
       ii=j0-1
       i3=ii/np
       jj3=modulo(i3-jj3min,(lr%d%n3+1))+1
       ii=ii-i3*np
       i2=ii/n1p1
       i0=ii-i2*n1p1
       i1=i0+j1-j0
       do i=i0,i1
          iitot=iitot+1
          !!if (i3>=ii3min .and. i3<=ii3max) then
          !!    index_in_global_c(i,i2,i3)=iitot
          if (jj3>=ii3min .and. jj3<=ii3max) then
              index_in_global_c(i,i2,jj3)=iitot
          end if
       end do
    end do 
    
    
    iitot=0
    istart=lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)
    iend=istart+lr%wfd%nseg_f-1
    do iseg=istart,iend
       j0=lr%wfd%keyglob(1,iseg)
       j1=lr%wfd%keyglob(2,iseg)
       ii=j0-1
       i3=ii/np
       jj3=modulo(i3-jj3min,(lr%d%n3+1))+1
       ii=ii-i3*np
       i2=ii/n1p1
       i0=ii-i2*n1p1
       i1=i0+j1-j0
       do i=i0,i1
          iitot=iitot+1
          !!if (i3>=ii3min .and. i3<=ii3max) then
          !!    index_in_global_f(i,i2,i3)=iitot
          if (jj3>=ii3min .and. jj3<=ii3max) then
              index_in_global_f(i,i2,jj3)=iitot
          end if
       end do
    end do

    call f_release_routine()
    
    end subroutine get_index_in_global2


    subroutine determine_communication_arrays(iproc, nproc, npsidim_orbs, orbs, nspin, lzd, &
               istartend_c, istartend_f, ii3min, ii3max, jj3min, index_in_global_c, index_in_global_f, &
               nvalp_c, nvalp_f,  nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c, &
               nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs, nspin, ii3min, ii3max, jj3min
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(2,0:nproc-1),intent(in) :: istartend_c, istartend_f
      integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,ii3min:ii3max),intent(in) :: index_in_global_c, index_in_global_f
      integer,intent(in) :: nvalp_c, nvalp_f
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
      
      ! Local variables
      integer :: iorb, iiorb, i1, i2, i3, ii, jproc, jproctarget, ierr, ilr, j0, j1, i0, i, ind, n1p1, np
      integer :: ii1, ii2, ii3, iseg, istart, iend, jj3
      integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
      character(len=*),     parameter :: subname='determine_communication_arrays'

      call f_routine(id='determine_communication_arrays')
    
      ! Determine values for mpi_alltoallv
      ! first nsendcounts
      nsendcounts_c=0
      nsendcounts_f=0
    
      !$omp parallel default(private) shared(ilr,nproc,orbs,lzd,index_in_global_c,istartend_c,nsendcounts_c,nsendcounts_f) &
      !$omp shared(istartend_f,index_in_global_f,n1p1,np,jj3min)
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          if (lzd%llr(ilr)%wfd%nseg_c>0) then
              !!n1p1=lzd%llr(ilr)%d%n1+1
              !!np=n1p1*(lzd%llr(ilr)%d%n2+1)
              n1p1=lzd%glr%d%n1+1
              np=n1p1*(lzd%glr%d%n2+1)
              !$omp do firstprivate(ilr) reduction(+:nsendcounts_c)
              do iseg=1,lzd%llr(ilr)%wfd%nseg_c
                  j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/np
                  jj3=modulo(i3-jj3min,(lzd%glr%d%n3+1))+1
                  ii=ii-i3*np
                  i2=ii/n1p1
                  i0=ii-i2*n1p1
                  i1=i0+j1-j0
                  !ii2=i2+lzd%llr(ilr)%ns2
                  !ii3=i3+lzd%llr(ilr)%ns3
                  ii2=i2
                  ii3=i3
                  do i=i0,i1
                      !ii1=i+lzd%llr(ilr)%ns1
                      ii1=i
                      !ind=index_in_global_c(ii1,ii2,ii3)
                      ind=index_in_global_c(ii1,ii2,jj3)
                      jproctarget=-1
                      do jproc=0,nproc-1
                          if(ind>=istartend_c(1,jproc) .and. ind<=istartend_c(2,jproc)) then
                              jproctarget=jproc
                              exit
                          end if
                      end do
                      if (jproctarget /= -1) &
                           nsendcounts_c(jproctarget)=nsendcounts_c(jproctarget)+1
                  end do
              end do
              !$omp end do
          end if
      end do
    
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
          iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
          if (istart<iend) then
              !n1p1=lzd%llr(ilr)%d%n1+1
              !np=n1p1*(lzd%llr(ilr)%d%n2+1)
              n1p1=lzd%glr%d%n1+1
              np=n1p1*(lzd%glr%d%n2+1)
              !$omp do firstprivate(ilr) reduction(+:nsendcounts_f)
              do iseg=istart,iend
                  j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/np
                  jj3=modulo(i3-jj3min,(lzd%glr%d%n3+1))+1
                  ii=ii-i3*np
                  i2=ii/n1p1
                  i0=ii-i2*n1p1
                  i1=i0+j1-j0
                  !ii2=i2+lzd%llr(ilr)%ns2
                  !ii3=i3+lzd%llr(ilr)%ns3
                  ii2=i2
                  ii3=i3
                  do i=i0,i1
                      !ii1=i+lzd%llr(ilr)%ns1
                      ii1=i
                      !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'f', ind)
                      !ind=index_in_global_f(ii1,ii2,ii3)
                      ind=index_in_global_f(ii1,ii2,jj3)
                      jproctarget=-1
                      do jproc=0,nproc-1
                          if(ind>=istartend_f(1,jproc) .and. ind<=istartend_f(2,jproc)) then
                              jproctarget=jproc
                              exit
                          end if
                      end do
                      if (jproctarget /= -1) &
                           nsendcounts_f(jproctarget)=nsendcounts_f(jproctarget)+1
                  end do
              end do
              !$omp end do
          end if
       end do
       !$omp end parallel
    
    
      ! The first check is to make sure that there is no stop in case this process has no orbitals (in which case
      ! npsidim_orbs is 1 and not 0 as assumed by the check)
      if(npsidim_orbs>1 .and. sum(nsendcounts_c)+7*sum(nsendcounts_f)/=npsidim_orbs) then
          write(*,'(a,2i10)') 'sum(nsendcounts_c)+sum(nsendcounts_f)/=npsidim_orbs', &
                              sum(nsendcounts_c)+sum(nsendcounts_f), npsidim_orbs
          stop
      end if
      !write(*,'(a,2i10)') 'sum(nsendcounts_c), 7*sum(nsendcounts_f)', sum(nsendcounts_c), 7*sum(nsendcounts_f)
      
      ! now nsenddspls
      nsenddspls_c(0)=0
      do jproc=1,nproc-1
          nsenddspls_c(jproc)=nsenddspls_c(jproc-1)+nsendcounts_c(jproc-1)
      end do
      nsenddspls_f(0)=0
      do jproc=1,nproc-1
          nsenddspls_f(jproc)=nsenddspls_f(jproc-1)+nsendcounts_f(jproc-1)
      end do
   
      ! now nrecvcounts
      ! use an mpi_alltoallv to gather the data
      nsendcounts_tmp = f_malloc(0.to.nproc-1,id='nsendcounts_tmp')
      nsenddspls_tmp = f_malloc(0.to.nproc-1,id='nsenddspls_tmp')
      nrecvcounts_tmp = f_malloc(0.to.nproc-1,id='nrecvcounts_tmp')
      nrecvdspls_tmp = f_malloc(0.to.nproc-1,id='nrecvdspls_tmp')
      nsendcounts_tmp=1
      nrecvcounts_tmp=1
      do jproc=0,nproc-1
          nsenddspls_tmp(jproc)=jproc
          nrecvdspls_tmp(jproc)=jproc
      end do
      if(nproc>1) then
          call mpialltoallv(nsendcounts_c(0), nsendcounts_tmp, nsenddspls_tmp, &
               nrecvcounts_c(0), nrecvcounts_tmp, nrecvdspls_tmp, bigdft_mpi%mpi_comm)
          call mpialltoallv(nsendcounts_f(0), nsendcounts_tmp, nsenddspls_tmp, &
               nrecvcounts_f(0), nrecvcounts_tmp, nrecvdspls_tmp, bigdft_mpi%mpi_comm)
      else
          nrecvcounts_c=nsendcounts_c
          nrecvcounts_f=nsendcounts_f
      end if
      call f_free(nsendcounts_tmp)
      call f_free(nsenddspls_tmp)
      call f_free(nrecvcounts_tmp)
      call f_free(nrecvdspls_tmp)
    
      ! now recvdspls
      nrecvdspls_c(0)=0
      do jproc=1,nproc-1
          nrecvdspls_c(jproc)=nrecvdspls_c(jproc-1)+nrecvcounts_c(jproc-1)
      end do
      nrecvdspls_f(0)=0
      do jproc=1,nproc-1
          nrecvdspls_f(jproc)=nrecvdspls_f(jproc-1)+nrecvcounts_f(jproc-1)
      end do
    
      if(sum(nrecvcounts_c)/=nspin*nvalp_c) then
          write(*,*) 'sum(nrecvcounts_c), nspin*nvalp_c', sum(nrecvcounts_c), nspin*nvalp_c
          stop 'sum(nrecvcounts_c)/=nspin*nvalp_c'
      end if
      if(sum(nrecvcounts_f)/=nspin*nvalp_f) then
          write(*,*) 'sum(nrecvcounts_f), nspin*nvalp_f', sum(nrecvcounts_f), nspin*nvalp_f
          stop 'sum(nrecvcounts_f)/=nspin*nvalp_f'
      end if

      call f_release_routine()
    
    end subroutine determine_communication_arrays


    subroutine determine_num_orbs_per_gridpoint_new(iproc, nproc, lzd, i3s, n3p, weightppp_c, weightppp_f, &
               jj3min, istartend_c, istartend_f, &
               istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
               weightp_c, weightp_f, nptsp_c, nptsp_f, &
               norb_per_gridpoint_c, norb_per_gridpoint_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc, i3s, n3p, nptsp_c, nptsp_f, jj3min, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
      type(local_zone_descriptors),intent(in):: lzd
      integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
      real(kind=8),intent(in):: weightp_c, weightp_f
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,1:max(1,n3p)),intent(in),target :: weightppp_c, weightppp_f
      integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
      integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f
      
      ! Local variables
      integer :: ii, i1, i2, i3, iipt, iseg, jj, j0, j1, iitot, i, istart, iend, i0
      integer :: icheck_c,icheck_f,iiorb_c,iiorb_f, npgp_c,npgp_f,np,n1p1, jj3
      integer :: window, i3min_c, i3max_c, i3min_f, i3max_f, size_of_double, ierr, jproc, is, ie, info, ncount
      integer,dimension(:),allocatable :: i3s_par, n3_par
      real(kind=8),dimension(:,:,:),pointer :: workrecv_c, workrecv_f
      !!integer,dimension(:),allocatable:: iseg_start_c, iseg_start_f
    
      call f_routine(id='determine_num_orbs_per_gridpoint_new')
    
      icheck_c = 0
      icheck_f = 0
      iiorb_f=0
      iiorb_c=0
      iipt=0
    
      n1p1=lzd%glr%d%n1+1
      np=n1p1*(lzd%glr%d%n2+1)


      !@NEW ######################################
       

       ! Initialize the MPI window
       if (nproc>1) then
           ! These arrays start at 1 instead of 0
           i3min_c = (lzd%glr%wfd%keyglob(1,istartp_seg_c)-1)/np + 1
           i3max_c = (lzd%glr%wfd%keyglob(2,iendp_seg_c)-1)/np + 1

           workrecv_c = f_malloc_ptr((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,i3min_c.to.i3max_c/),id='workrecv')

           i3s_par = f_malloc0(0.to.nproc-1,id='i3s_par')
           i3s_par(iproc)=i3s
           call mpiallred(i3s_par,mpi_sum, comm=bigdft_mpi%mpi_comm)
           n3_par = f_malloc0(0.to.nproc-1,id='n3_par')
           n3_par(iproc)=n3p
           call mpiallred(n3_par, mpi_sum, comm=bigdft_mpi%mpi_comm)

           !!call mpi_type_size(mpi_double_precision, size_of_double, ierr)
           !!call mpi_info_create(info, ierr)
           !!call mpi_info_set(info, "no_locks", "true", ierr)
           !!call mpi_win_create(weightppp_c(0,0,1), &
           !!     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p*size_of_double,kind=mpi_address_kind), size_of_double, &
           !!     info, bigdft_mpi%mpi_comm, window, ierr)
           !!call mpi_info_free(info, ierr)
           !!call mpi_win_fence(mpi_mode_noprecede, window, ierr)
           call mpi_type_size(mpi_double_precision, size_of_double, ierr)
           window = mpiwindow((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p, weightppp_c(0,0,1), bigdft_mpi%mpi_comm)

           do jproc=0,nproc-1
               ! Check whether ther is an overlap
               is = max(i3min_c,i3s_par(jproc))
               ie = min(i3max_c,i3s_par(jproc)+n3_par(jproc)-1)
               if (ie-is>=0) then
                   ncount = (ie-is+1)*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
                   !write(*,'(9(a,i0),a)') 'process ',iproc,'(i3min=',i3min_c,',i3max=',i3max_c,') gets ',(ie-is+1), &
                   !                    ' lines at ',is,' from ',is-i3s_par(jproc)+1,' on process ', &
                   !                    jproc,'(i3s=',i3s_par(jproc),',n3p=',n3_par(jproc),')'
                   !if (iproc/=jproc) then
                       call mpi_get(workrecv_c(0,0,is), ncount, mpi_double_precision, jproc, &
                            int((is-i3s_par(jproc))*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1),kind=mpi_address_kind), &
                            ncount, mpi_double_precision, window, ierr)
                   !else
                   !    ncount = (lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
                   !    call vcopy(ncount, weightppp_c(0,0,is-i3s_par(jproc)+1), 1, workrecv_c(0,0,is), 1)
                   !end if
               end if
           end do
           ! Synchronize the communication
           !!call mpi_win_fence(0, window, ierr)
           !!call mpi_win_free(window, ierr)
           call mpi_fenceandfree(window)

       else
           workrecv_c => weightppp_c
       end if


       icheck_c = 0
       iiorb_c = 0
       do iseg=1,lzd%glr%wfd%nseg_c !istartp_seg_c,iendp_seg_c
           jj=lzd%glr%wfd%keyvglob(iseg)
           j0=lzd%glr%wfd%keyglob(1,iseg)
           j1=lzd%glr%wfd%keyglob(2,iseg)
           ii=j0-1
           i3=ii/np
           !jj3=modulo(i3-jj3min,(lzd%glr%d%n3+1))+1
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           do i=i0,i1
               iitot=jj+i-i0
               if(iitot>=istartend_c(1,iproc) .and. iitot<=istartend_c(2,iproc)) then
                   !write(1100+iproc,'(a,4i8,f10.1)') 'iproc, i, i2, i3, workrecv_c_c(i,i2,i3+1)', iproc, i, i2, i3, workrecv_c(i,i2,i3+1)
                   icheck_c = icheck_c + 1
                   iipt=jj-istartend_c(1,iproc)+i-i0+1
                   npgp_c = nint(sqrt(workrecv_c(i,i2,i3+1)))
                   iiorb_c=iiorb_c+nint(workrecv_c(i,i2,i3+1))
                   !!npgp_c = nint(sqrt(workrecv_c(i,i2,jj3)))
                   !!iiorb_c=iiorb_c+nint(workrecv_c(i,i2,jj3))
                   norb_per_gridpoint_c(iipt)=npgp_c
               end if
           end do
       end do

       if (nproc>1) then
           call f_free_ptr(workrecv_c)
           call f_free(i3s_par)
           call f_free(n3_par)
       end if


      !@ENDNEW ######################################
    
      if(icheck_c/=nptsp_c) stop 'icheck_c/=nptsp_c'
      if(iiorb_c/=nint(weightp_c)) then
          write(*,*) 'iiorb_c, nint(weightp_c)', iiorb_c, nint(weightp_c)
          stop 'iiorb_c/=weightp_c'
      end if
    
    

      !@NEW ######################################
       

       if (nproc>1) then
           ! These arrays start at one instead of 0
           i3min_f = (lzd%glr%wfd%keyglob(1,istartp_seg_f)-1)/np + 1
           i3max_f = (lzd%glr%wfd%keyglob(2,iendp_seg_f)-1)/np + 1

           workrecv_f = f_malloc_ptr((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,i3min_f.to.i3max_f/),id='workrecv')

           i3s_par = f_malloc0(0.to.nproc-1,id='i3s_par')
           i3s_par(iproc)=i3s
           call mpiallred(i3s_par, mpi_sum, comm=bigdft_mpi%mpi_comm)
           n3_par = f_malloc0(0.to.nproc-1,id='n3_par')
           n3_par(iproc)=n3p
           call mpiallred(n3_par, mpi_sum, comm=bigdft_mpi%mpi_comm)

           ! Initialize the MPI window
           !!call mpi_type_size(mpi_double_precision, size_of_double, ierr)
           !!call mpi_info_create(info, ierr)
           !!call mpi_info_set(info, "no_locks", "true", ierr)
           !!call mpi_win_create(weightppp_f(0,0,1), &
           !!     int((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p*size_of_double,kind=mpi_address_kind), size_of_double, &
           !!     info, bigdft_mpi%mpi_comm, window, ierr)
           !!call mpi_info_free(info, ierr)
           !!call mpi_win_fence(mpi_mode_noprecede, window, ierr)
           call mpi_type_size(mpi_double_precision, size_of_double, ierr)
           window = mpiwindow((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*n3p, weightppp_f(0,0,1), bigdft_mpi%mpi_comm)

           do jproc=0,nproc-1
               ! Check whether ther is an overlap
               is = max(i3min_f,i3s_par(jproc))
               ie = min(i3max_f,i3s_par(jproc)+n3_par(jproc)-1)
               if (ie-is>=0) then
                   ncount = (ie-is+1)*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
                   !!write(*,'(9(a,i0),a)') 'process ',iproc,'(i3min=',i3min_f,',i3max=',i3max_f,') gets ',(ie-is+1), &
                   !!                    ' lines at ',is,' from ',is-i3s_par(jproc)+1,' on process ', &
                   !!                    jproc,'(i3s=',i3s_par(jproc),',n3p=',n3_par(jproc),')'
                   call mpi_get(workrecv_f(0,0,is), ncount, mpi_double_precision, jproc, &
                        int((is-i3s_par(jproc))*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1),kind=mpi_address_kind), &
                        ncount, mpi_double_precision, window, ierr)
               end if
           end do
           ! Synchronize the communication
           !!call mpi_win_fence(0, window, ierr)
           !!call mpi_win_free(window, ierr)
           call mpi_fenceandfree(window)
       else
           workrecv_f => weightppp_f
       end if


       icheck_f = 0
       iiorb_f = 0
       !do iseg=lzd%glr%wfd%nseg_c+1,lzd%glr%wfd%nseg_c+lzd%glr%wfd%nseg_f !istartp_seg_f,iendp_seg_f
       do iseg=istartp_seg_f,iendp_seg_f
           jj=lzd%glr%wfd%keyvglob(iseg)
           j0=lzd%glr%wfd%keyglob(1,iseg)
           j1=lzd%glr%wfd%keyglob(2,iseg)
           ii=j0-1
           i3=ii/np
           !!jj3=modulo(i3-jj3min,(lzd%glr%d%n3+1))+1
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           do i=i0,i1
               iitot=jj+i-i0
               if(iitot>=istartend_f(1,iproc) .and. iitot<=istartend_f(2,iproc)) then
                   !!write(1100+iproc,'(a,4i8,f10.1)') 'iproc, i, i2, i3, workrecv_f(i,i2,i3+1)', iproc, i, i2, i3, workrecv_f(i,i2,i3+1)
                   icheck_f = icheck_f + 1
                   iipt=jj-istartend_f(1,iproc)+i-i0+1
                   npgp_f = nint(sqrt(workrecv_f(i,i2,i3+1)))
                   iiorb_f=iiorb_f+nint(workrecv_f(i,i2,i3+1))
                   !!npgp_f = nint(sqrt(workrecv_f(i,i2,jj3)))
                   !!iiorb_f=iiorb_f+nint(workrecv_f(i,i2,jj3))
                   norb_per_gridpoint_f(iipt)=npgp_f
                   !write(*,'(a,6i8)') 'iipt, i, i2, i3, jj3, npgp_f', iipt, i, i2, i3, jj3, npgp_f
                   !write(*,'(a,6i8)') 'iipt, i, i2, i3, npgp_f', iipt, i, i2, i3, npgp_f
               end if
           end do
       end do

       if (nproc>1) then
           call f_free_ptr(workrecv_f)
           call f_free(i3s_par)
           call f_free(n3_par)
       end if


      !@ENDNEW ######################################
    
      if(icheck_f/=nptsp_f) stop 'icheck_f/=nptsp_f'
      if(iiorb_f/=nint(weightp_f)) then
          write(*,*) 'iiorb_f, weightp_f', iiorb_f, weightp_f
          stop 'iiorb_f/=weightp_f'
      end if

      call f_release_routine()
    
    end subroutine determine_num_orbs_per_gridpoint_new



    subroutine get_switch_indices(iproc, nproc, orbs, lzd, nspin, &
               nptsp_c, nptsp_f, norb_per_gridpoint_c, norb_per_gridpoint_f, &
               ndimpsi_c, ndimpsi_f, istartend_c, istartend_f, &
               istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
               nsendcounts_c, nsenddspls_c, ndimind_c, nrecvcounts_c, nrecvdspls_c, &
               nsendcounts_f, nsenddspls_f, ndimind_f, nrecvcounts_f, nrecvdspls_f, &
               ii3min, ii3max, jj3min, index_in_global_c, index_in_global_f, &
               weightp_c, weightp_f,  isendbuf_c, irecvbuf_c, isendbuf_f, irecvbuf_f, &
               indexrecvorbital_c, iextract_c, iexpand_c, indexrecvorbital_f, iextract_f, iexpand_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin, nptsp_c, nptsp_f, ndimpsi_c, ndimpsi_f, ndimind_c,ndimind_f
      integer,intent(in) :: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, ii3min, ii3max, jj3min
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(nptsp_c),intent(in):: norb_per_gridpoint_c
      integer,dimension(nptsp_f),intent(in):: norb_per_gridpoint_f
      integer,dimension(2,0:nproc-1),intent(in) :: istartend_c, istartend_f
      integer,dimension(0:nproc-1),intent(in) :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
      integer,dimension(0:nproc-1),intent(in) :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
      integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,ii3min:ii3max),intent(in) :: index_in_global_c, index_in_global_f
      real(kind=8),intent(in) :: weightp_c, weightp_f
      integer,dimension(ndimpsi_c),intent(out) :: isendbuf_c, irecvbuf_c
      integer,dimension(ndimpsi_f),intent(out) :: isendbuf_f, irecvbuf_f
      integer,dimension(ndimind_c),intent(out) :: indexrecvorbital_c, iextract_c, iexpand_c
      integer,dimension(ndimind_f),intent(out) :: indexrecvorbital_f, iextract_f, iexpand_f
      
      ! Local variables
      integer :: i, iorb, iiorb, i1, i2, i3, ind, jproc, jproctarget, ii, ierr, iseg, iitot, ilr, n1p1, np
      integer :: i3min_c, i3max_c, i3min_f, i3max_f, jj3
      integer :: istart, iend, indglob, ii1, ii2, ii3, j1, i0, j0, ipt
      integer,dimension(:),allocatable :: nsend_c,nsend_f, indexsendorbital2, indexrecvorbital2
      integer,dimension(:),allocatable :: gridpoint_start_c, gridpoint_start_f, gridpoint_start_tmp_c, gridpoint_start_tmp_f
      real(kind=8),dimension(:,:,:),allocatable :: weight_c, weight_f
      integer,dimension(:),allocatable :: indexsendorbital_c, indexsendbuf_c, indexrecvbuf_c
      integer,dimension(:),allocatable :: indexsendorbital_f, indexsendbuf_f, indexrecvbuf_f
      character(len=*),parameter :: subname='get_switch_indices'
    
    
      call f_routine(id='get_switch_indices')
      
      indexsendorbital_c = f_malloc(max(ndimpsi_c,1),id='indexsendorbital_c')
      indexsendbuf_c = f_malloc(max(ndimpsi_c,1),id='indexsendbuf_c')
      indexrecvbuf_c = f_malloc(sum(nrecvcounts_c),id='indexrecvbuf_c')
      indexsendorbital_f = f_malloc(max(ndimpsi_f,1),id='indexsendorbital_f')
      indexsendbuf_f = f_malloc(max(ndimpsi_f,1),id='indexsendbuf_f')
      indexrecvbuf_f = f_malloc(sum(nrecvcounts_f),id='indexrecvbuf_f')
      gridpoint_start_c = f_malloc(istartend_c(1,iproc).to.istartend_c(2,iproc),id='gridpoint_start_c')
      gridpoint_start_f = f_malloc(istartend_f(1,iproc).to.istartend_f(2,iproc),id='gridpoint_start_f')
      if (nspin==2) then
          !!gridpoint_start_tmp_c = f_malloc((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1),id='gridpoint_start_tmp_c')
          !!gridpoint_start_tmp_f = f_malloc((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1),id='gridpoint_start_tmp_f')
          gridpoint_start_tmp_c = f_malloc(istartend_c(1,iproc).to.istartend_c(2,iproc),id='gridpoint_start_tmp_c')
          gridpoint_start_tmp_f = f_malloc(istartend_f(1,iproc).to.istartend_f(2,iproc),id='gridpoint_start_tmp_f')
      end if
      gridpoint_start_c=-1
      gridpoint_start_f=-1


      !!write(*,'(a,i7,4i9)') 'iproc, i3min_c, i3max_c, i3min_f, i3max_f', iproc, i3min_c, i3max_c, i3min_f, i3max_f
    
      !write(*,*) 'ndimpsi_f, sum(nrecvcounts_f)', ndimpsi_f, sum(nrecvcounts_f)
    
      nsend_c = f_malloc(0.to.nproc-1,id='nsend_c')
      nsend_f = f_malloc(0.to.nproc-1,id='nsend_f')
    
      nsend_c=0
      nsend_f=0
    
      !$omp parallel default(private) shared(orbs,lzd,index_in_global_c,index_in_global_f,istartend_c,istartend_f)&
      !$omp shared(nsend_c,nsend_f,nsenddspls_c,nsenddspls_f,ndimpsi_c,ndimpsi_f,nsendcounts_c,nsendcounts_f,nproc,jj3min) &
      !$omp shared(isendbuf_c,isendbuf_f,indexsendbuf_c,indexsendbuf_f,indexsendorbital_c,indexsendorbital_f)
    
      !$omp sections
      !$omp section
      iitot=0
     
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          !n1p1=lzd%llr(ilr)%d%n1+1
          !np=n1p1*(lzd%llr(ilr)%d%n2+1)
          n1p1=lzd%glr%d%n1+1
          np=n1p1*(lzd%glr%d%n2+1)
          do iseg=1,lzd%llr(ilr)%wfd%nseg_c
              !jj=lzd%llr(ilr)%wfd%keyvglob(iseg)
              j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
              j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
              ii=j0-1
              i3=ii/np
              jj3=modulo(i3-jj3min,(lzd%glr%d%n3+1))+1
              ii=ii-i3*np
              i2=ii/n1p1
              i0=ii-i2*n1p1
              i1=i0+j1-j0
              !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
              !ii2=i2+lzd%llr(ilr)%ns2
              !ii3=i3+lzd%llr(ilr)%ns3
              ii2=i2
              ii3=i3
              do i=i0,i1
                  !ii1=i+lzd%llr(ilr)%ns1
                  ii1=i
                  !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'c', indglob)
                  !indglob=index_in_global_c(ii1,ii2,ii3)
                  indglob=index_in_global_c(ii1,ii2,jj3)
                  iitot=iitot+1
               jproctarget=-1
               do jproc=0,nproc-1
                  if(indglob>=istartend_c(1,jproc) .and. indglob<=istartend_c(2,jproc)) then
                     jproctarget=jproc
                     exit
                  end if
               end do
               !write(600+iproc,'(a,2(i0,1x),i0,a,i0)') 'point ',ii1,ii2,ii3,' goes to process ',jproctarget
              
               if (jproctarget/=-1) then
                  nsend_c(jproctarget)=nsend_c(jproctarget)+1
                  ind=nsenddspls_c(jproctarget)+nsend_c(jproctarget)
                  isendbuf_c(iitot)=ind
                  indexsendbuf_c(ind)=indglob
                  indexsendorbital_c(iitot)=iiorb
               end if
               !indexsendorbital(ind)=iiorb
            end do
         end do
      end do
      !write(*,*) 'iitot,ndimpsi_c',iitot,ndimpsi_c
      if(iitot/=ndimpsi_c) stop 'iitot/=ndimpsi_c'

      ! Check only up to ndimpsi_c in case it was allocated with size 1 to avoid an allocation with size 0
      if (minval(indexsendorbital_c(1:ndimpsi_c))<1) then
          write(*,*) 'minval(indexsendorbital_c)',minval(indexsendorbital_c)
          stop 'minval(indexsendorbital_c)<1'
      end if
    
      !check
      do jproc=0,nproc-1
          if(nsend_c(jproc)/=nsendcounts_c(jproc)) stop 'nsend_c(jproc)/=nsendcounts_c(jproc)'
      end do
    
    
      !$omp section
      ! fine part
      iitot=0
      do iorb=1,orbs%norbp
         iiorb=orbs%isorb+iorb
         ilr=orbs%inwhichlocreg(iiorb)
         istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
         iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
         !n1p1=lzd%llr(ilr)%d%n1+1
         !np=n1p1*(lzd%llr(ilr)%d%n2+1)
         n1p1=lzd%glr%d%n1+1
         np=n1p1*(lzd%glr%d%n2+1)
         do iseg=istart,iend
            !jj=lzd%llr(ilr)%wfd%keyvglob(iseg)
            j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
            j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
            ii=j0-1
            i3=ii/np
            jj3=modulo(i3-jj3min,(lzd%glr%d%n3+1))+1
            ii=ii-i3*np
            i2=ii/n1p1
            i0=ii-i2*n1p1
            i1=i0+j1-j0
            !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
            !ii2=i2+lzd%llr(ilr)%ns2
            !ii3=i3+lzd%llr(ilr)%ns3
            ii2=i2
            ii3=i3
            do i=i0,i1
               !ii1=i+lzd%llr(ilr)%ns1
               ii1=i
               !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'f', indglob)
               !indglob=index_in_global_f(ii1,ii2,ii3)
               indglob=index_in_global_f(ii1,ii2,jj3)
               iitot=iitot+1
               jproctarget=-1
               do jproc=0,nproc-1
                  if(indglob>=istartend_f(1,jproc) .and. indglob<=istartend_f(2,jproc)) then
                     jproctarget=jproc
                     exit
                  end if
               end do
               if (jproctarget/=-1) then
                  nsend_f(jproctarget)=nsend_f(jproctarget)+1
                  ind=nsenddspls_f(jproctarget)+nsend_f(jproctarget)
                  isendbuf_f(iitot)=ind
                  indexsendbuf_f(ind)=indglob
                  !write(*,*) 'iiorb, ind, indglob', iiorb, ind, indglob
                  indexsendorbital_f(iitot)=iiorb
               end if
               !indexsendorbital(ind)=iiorb
            end do
         end do
     
      end do
      
      if(iitot/=ndimpsi_f) stop 'iitot/=ndimpsi_f'

      ! Check only up to ndimpsi_f in case it was allocated with size 1 to avoid an allocation with size 0
      if (minval(indexsendorbital_f(1:ndimpsi_f))<1) then
          write(*,*) 'minval(indexsendorbital_f)',minval(indexsendorbital_f)
          stop 'minval(indexsendorbital_f)<1'
      end if
    
      !$omp end sections
      !$omp end parallel
    
      !check
      do jproc=0,nproc-1
          !write(*,*) 'nsend(jproc), nsendcounts_f(jproc)', nsend(jproc), nsendcounts_f(jproc)
          if(nsend_f(jproc)/=nsendcounts_f(jproc)) stop 'nsend_f(jproc)/=nsendcounts_f(jproc)'
      end do
    
      indexsendorbital2 = f_malloc(max(1,ndimpsi_c),id='indexsendorbital2')
      call vcopy(ndimpsi_c, indexsendorbital_c(1), 1, indexsendorbital2(1), 1)
      do i=1,ndimpsi_c
          ind=isendbuf_c(i)
          indexsendorbital_c(ind)=indexsendorbital2(i)
      end do
    
      ! Inverse of isendbuf
      call get_reverse_indices(ndimpsi_c, isendbuf_c, irecvbuf_c)
    
      call f_free(indexsendorbital2)
      indexsendorbital2 = f_malloc(max(1,ndimpsi_f),id='indexsendorbital2')
      call vcopy(ndimpsi_f, indexsendorbital_f(1), 1, indexsendorbital2(1), 1)
      do i=1,ndimpsi_f
          ind=isendbuf_f(i)
          indexsendorbital_f(ind)=indexsendorbital2(i)
      end do
    
      ! Inverse of isendbuf
    
      call get_reverse_indices(ndimpsi_f, isendbuf_f, irecvbuf_f)
      call f_free(indexsendorbital2)
    
    
      if(nproc>1) then
          ! Communicate indexsendbuf
          call mpialltoallv(indexsendbuf_c(1), nsendcounts_c, nsenddspls_c, &
               indexrecvbuf_c(1), nrecvcounts_c, nrecvdspls_c, bigdft_mpi%mpi_comm)
          ! Communicate indexsendorbitals
          call mpialltoallv(indexsendorbital_c(1), nsendcounts_c, nsenddspls_c, &
               indexrecvorbital_c(1), nrecvcounts_c, nrecvdspls_c, bigdft_mpi%mpi_comm)
    
          ! Communicate indexsendbuf
          call mpialltoallv(indexsendbuf_f(1), nsendcounts_f, nsenddspls_f, &
               indexrecvbuf_f(1), nrecvcounts_f, nrecvdspls_f, bigdft_mpi%mpi_comm)
          ! Communicate indexsendorbitals
          call mpialltoallv(indexsendorbital_f(1), nsendcounts_f, nsenddspls_f, &
               indexrecvorbital_f(1), nrecvcounts_f, nrecvdspls_f, bigdft_mpi%mpi_comm)
       else
           indexrecvbuf_c=indexsendbuf_c
           indexrecvorbital_c=indexsendorbital_c
           indexrecvbuf_f=indexsendbuf_f
           indexrecvorbital_f=indexsendorbital_f
       end if
    
        

      ! gridpoint_start is the starting index of a given grid point in the overall array
      ii=1
      do ipt=1,nptsp_c
          i=ipt+istartend_c(1,iproc)-1
          if (norb_per_gridpoint_c(ipt)>0) then
              gridpoint_start_c(i)=ii
          else
              gridpoint_start_c(i)=0
          end if
          ii=ii+norb_per_gridpoint_c(ipt)
      end do

      ii=1
      do ipt=1,nptsp_f
          i=ipt+istartend_f(1,iproc)-1
          if (norb_per_gridpoint_f(ipt)>0) then
              gridpoint_start_f(i)=ii
          else
              !write(*,*) 'SET gridpoint_start_f TO ZERO: i',i
              gridpoint_start_f(i)=0
          end if
          ii=ii+norb_per_gridpoint_f(ipt)
      end do



      if (nspin==2) then
          gridpoint_start_tmp_c=gridpoint_start_c
          gridpoint_start_tmp_f=gridpoint_start_f
      end if
        
    
      if(maxval(gridpoint_start_c)>sum(nrecvcounts_c)) stop '1: maxval(gridpoint_start_c)>sum(nrecvcounts_c)'
      if(maxval(gridpoint_start_f)>sum(nrecvcounts_f)) stop '1: maxval(gridpoint_start_f)>sum(nrecvcounts_f)'

      ! Rearrange the communicated data
      if (nspin==1) then
          do i=1,sum(nrecvcounts_c)
              ii=indexrecvbuf_c(i)
              ind=gridpoint_start_c(ii)
              iextract_c(i)=ind
              gridpoint_start_c(ii)=gridpoint_start_c(ii)+1  
          end do
      else
          do i=1,sum(nrecvcounts_c)
              ii=indexrecvbuf_c(i)
              ind=gridpoint_start_c(ii)
              if (gridpoint_start_c(ii)-gridpoint_start_tmp_c(ii)+1>norb_per_gridpoint_c(ii-istartend_c(1,iproc)+1)) then
                  ! orbitals which fulfill this condition are down orbitals which should be put at the end
                  !ind = ind + ((ndimind_c+ndimind_f)/2-norb_per_gridpoint_c(ii-istartend_c(1,iproc)+1))
                  ind = ind + (ndimind_c/2-norb_per_gridpoint_c(ii-istartend_c(1,iproc)+1))
              end if
              iextract_c(i)=ind
              gridpoint_start_c(ii)=gridpoint_start_c(ii)+1  
          end do
      end if
      !write(*,'(a,2i12)') 'sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)', sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)
      !if(sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)) stop 'sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)'
      if(maxval(iextract_c)>sum(nrecvcounts_c)) then
          stop 'maxval(iextract_c)>sum(nrecvcounts_c)'
      end if
      if(minval(iextract_c)<1) stop 'minval(iextract_c)<1'
    
      ! Rearrange the communicated data
      if (nspin==1) then
          do i=1,sum(nrecvcounts_f)
              ii=indexrecvbuf_f(i)
              ind=gridpoint_start_f(ii)
              !if(ind<1) write(*,*) 'ERROR: ind<1, i, ii, ind',i, ii, ind
              iextract_f(i)=ind
              gridpoint_start_f(ii)=gridpoint_start_f(ii)+1  
          end do
      else
          do i=1,sum(nrecvcounts_f)
              ii=indexrecvbuf_f(i)
              ind=gridpoint_start_f(ii)
              if (gridpoint_start_f(ii)-gridpoint_start_tmp_f(ii)+1>norb_per_gridpoint_f(ii-istartend_f(1,iproc)+1)) then
                  ! orbitals which fulfill this condition are down orbitals which should be put at the end
                  ind = ind + (ndimind_f/2-norb_per_gridpoint_f(ii-istartend_f(1,iproc)+1))
              end if
              iextract_f(i)=ind
              gridpoint_start_f(ii)=gridpoint_start_f(ii)+1  
          end do
      end if
      if(maxval(iextract_f)>sum(nrecvcounts_f)) stop 'maxval(iextract_f)>sum(nrecvcounts_f)'
      if(minval(iextract_f)<1) stop 'minval(iextract_f)<1'
        
    
      ! Get the array to transfrom back the data
      call get_reverse_indices(sum(nrecvcounts_c), iextract_c, iexpand_c)
      call get_reverse_indices(sum(nrecvcounts_f), iextract_f, iexpand_f)
          
    
      indexrecvorbital2 = f_malloc(sum(nrecvcounts_c),id='indexrecvorbital2')
      indexrecvorbital2=indexrecvorbital_c
      do i=1,sum(nrecvcounts_c)
          ind=iextract_c(i)
          indexrecvorbital_c(ind)=indexrecvorbital2(i)
      end do
      call f_free(indexrecvorbital2)
      indexrecvorbital2 = f_malloc(sum(nrecvcounts_f),id='indexrecvorbital2')
      indexrecvorbital2=indexrecvorbital_f
      do i=1,sum(nrecvcounts_f)
          ind=iextract_f(i)
          indexrecvorbital_f(ind)=indexrecvorbital2(i)
      end do
      call f_free(indexrecvorbital2)
    
    
      if(minval(indexrecvorbital_c)<1) stop 'minval(indexrecvorbital_c)<1'
      if(maxval(indexrecvorbital_c)>orbs%norb) stop 'maxval(indexrecvorbital_c)>orbs%norb'
      if(minval(indexrecvorbital_f)<1) stop 'minval(indexrecvorbital_f)<1'
      if(maxval(indexrecvorbital_f)>orbs%norb) stop 'maxval(indexrecvorbital_f)>orbs%norb'
    

      call f_free(indexsendorbital_c)
      call f_free(indexsendbuf_c)
      call f_free(indexrecvbuf_c)
      call f_free(indexsendorbital_f)
      call f_free(indexsendbuf_f)
      call f_free(indexrecvbuf_f)
      call f_free(gridpoint_start_c)
      call f_free(gridpoint_start_f)
      if (nspin==2) then
          call f_free(gridpoint_start_tmp_c)
          call f_free(gridpoint_start_tmp_f)
      end if
      call f_free(nsend_c)
      call f_free(nsend_f)


      call f_release_routine()
    
    end subroutine get_switch_indices



    subroutine get_reverse_indices(n, indices, reverse_indices)
      use module_base
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: n
      integer,dimension(n),intent(in) :: indices
      integer,dimension(n),intent(out) :: reverse_indices
    
      ! Local variables
      integer :: i, j, m, j0, j1, j2, j3
    
      !$omp parallel default(private) &
      !$omp shared(n, m, indices, reverse_indices)
    
      m=mod(n,4)
      if (m/=0) then
          do i=1,m
              j=indices(i)
              reverse_indices(j)=i
          end do
      end if
    
      !$omp do
      do i=m+1,n,4
          j0=indices(i+0)
          reverse_indices(j0)=i+0
          j1=indices(i+1)
          reverse_indices(j1)=i+1
          j2=indices(i+2)
          reverse_indices(j2)=i+2
          j3=indices(i+3)
          reverse_indices(j3)=i+3
      end do
      !$omp end do
    
      !$omp end parallel
    
      !!do i=1,n
      !!    j=indices(i)
      !!    reverse_indices(j)=i
      !!end do
    
    end subroutine get_reverse_indices



    subroutine get_gridpoint_start(iproc, nproc, lzd, ndimind_c, nrecvcounts_c, ndimind_f, nrecvcounts_f, &
               indexrecvbuf_c, indexrecvbuf_f, &
               i3min_c, i3max_c, weight_c, i3min_f, i3max_f, weight_f, gridpoint_start_c, gridpoint_start_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc,ndimind_c,ndimind_f,i3min_c,i3max_c,i3min_f,i3max_f
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1),intent(in) :: nrecvcounts_c, nrecvcounts_f
      integer,dimension(ndimind_c),intent(in) :: indexrecvbuf_c
      integer,dimension(ndimind_f),intent(in) :: indexrecvbuf_f
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,i3min_c:i3max_c),intent(out) :: weight_c
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,i3min_f:i3max_f),intent(out) :: weight_f
      integer,dimension((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)),intent(out) :: gridpoint_start_c, gridpoint_start_f
      
      ! Local variables
      integer :: i, ii, jj, i1, i2, i3, n1p1, np
    
    
      !weight_c=0.d0
      !call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(i3max_c-i3min_c+1), weight_c(0,0,i3min_c))
      !call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(i3max_f-i3min_f+1), weight_f(0,0,i3min_f))
      call f_zero(weight_c)
      call f_zero(weight_f)

      n1p1=lzd%glr%d%n1+1
      np=n1p1*(lzd%glr%d%n2+1)

      !!$omp parallel default(private) shared(lzd,nrecvcounts_c,indexrecvbuf_c,weight_c,gridpoint_start_c) &
      !!$omp shared(nrecvcounts_f,indexrecvbuf_f,weight_f,gridpoint_start_f,np,n1p1)
    
      !!$omp sections
      !!$omp section
      do i=1,sum(nrecvcounts_c)
          ii=indexrecvbuf_c(i)
          !write(650+iproc,*) i, ii
          jj=ii-1
          i3=jj/np
          jj=jj-i3*np
          i2=jj/n1p1
          i1=jj-i2*n1p1
          weight_c(i1,i2,i3)=weight_c(i1,i2,i3)+1.d0
      end do
    
      !write(*,*) 'in get_gridpoint_start: maxval(weight_c)', maxval(weight_c)
    
      ii=1
      i=0
      !gridpoint_start_c=0
      do i3=0,lzd%glr%d%n3
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  i=i+1
                  if(weight_c(i1,i2,i3)>0.d0) then
                      gridpoint_start_c(i)=ii
                      ii=ii+nint(weight_c(i1,i2,i3))
                  else
                      gridpoint_start_c(i) = 0
                  end if
              end do
          end do
      end do
    
      !!$omp section
     
      do i=1,sum(nrecvcounts_f)
          ii=indexrecvbuf_f(i)
          write(*,*) 'i, ii', i, ii
          jj=ii-1
          i3=jj/np
          jj=jj-i3*np
          i2=jj/n1p1
          i1=jj-i2*n1p1
          weight_f(i1,i2,i3)=weight_f(i1,i2,i3)+1.d0
      end do
    
    
      ii=1
      i=0
      !gridpoint_start_f=0
      do i3=0,lzd%glr%d%n3
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  i=i+1
                  if(weight_f(i1,i2,i3)>0.d0) then
                      gridpoint_start_f(i)=ii
                      ii=ii+nint(weight_f(i1,i2,i3))
                  else
                      gridpoint_start_f(i)=0
                  end if
              end do
          end do
      end do
    
      !!$omp end sections
      !!$omp end parallel
    
    end subroutine get_gridpoint_start


    !> The sumrho routines
    subroutine init_comms_linear_sumrho(iproc, nproc, lzd, orbs, nspin, nscatterarr, collcom_sr)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      type(comms_linear),intent(inout) :: collcom_sr
    
      ! Local variables
      integer :: ipt, ii
      real(kind=8) :: weight_tot, weight_ideal
      integer(kind=8),dimension(:,:),allocatable :: istartend
      character(len=*),parameter :: subname='init_comms_linear_sumrho'
      real(kind=8),dimension(:),allocatable :: weights_per_slice, weights_per_zpoint
    
      ! Note: all weights are double precision to avoid integer overflow
      call timing(iproc,'init_collco_sr','ON')
    
      istartend = f_malloc((/ 1.to.2, 0.to.nproc-1 /),id='istartend')
      weights_per_slice = f_malloc(0.to.nproc-1,id='weights_per_slice')
      weights_per_zpoint = f_malloc(lzd%glr%d%n3i,id='weights_per_zpoint')
      call get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, weight_tot, weight_ideal, &
           weights_per_slice, weights_per_zpoint)

    
      call assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, weights_per_slice, &
           lzd, orbs, nscatterarr, istartend, collcom_sr%nptsp_c)

    
      call f_free(weights_per_slice)
    

      call allocate_MPI_communication_arrays(nproc, collcom_sr, only_coarse=.true.)

      call determine_communication_arrays_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, istartend, &
           collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, &
           collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, collcom_sr%ndimpsi_c)

      !Now set some integers in the collcomm structure
      collcom_sr%ndimind_c = sum(collcom_sr%nrecvcounts_c)

      call allocate_local_comms_cubic(collcom_sr, only_coarse=.true.)
    
      call determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, &
           istartend, weight_tot, weights_per_zpoint, collcom_sr%norb_per_gridpoint_c)
    
      ! Some check
      ii=sum(collcom_sr%norb_per_gridpoint_c)
      if (nspin*ii/=collcom_sr%ndimind_c) then
          write(*,*) 'nspin*ii/=collcom_sr%ndimind_c', ii, collcom_sr%ndimind_c
          stop 'nspin*ii/=collcom_sr%ndimind_c'
      end if
    
    
      collcom_sr%psit_c=f_malloc_ptr(collcom_sr%ndimind_c,id='collcom_sr%psit_c')
    
      call get_switch_indices_sumrho(iproc, nproc, collcom_sr%nptsp_c, collcom_sr%ndimpsi_c, collcom_sr%ndimind_c, lzd, &
           orbs, nspin, istartend, collcom_sr%norb_per_gridpoint_c, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, &
           collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, collcom_sr%isendbuf_c, collcom_sr%irecvbuf_c, &
           collcom_sr%iextract_c, collcom_sr%iexpand_c, collcom_sr%indexrecvorbital_c)
    
      ! These variables are used in various subroutines to speed up the code
      collcom_sr%isptsp_c(1) = 0
      do ipt=2,collcom_sr%nptsp_c
            collcom_sr%isptsp_c(ipt) = collcom_sr%isptsp_c(ipt-1) + collcom_sr%norb_per_gridpoint_c(ipt-1)
      end do
    
      !!call allocate_MPI_comms_cubic_repartition(nproc, collcom_sr)
    
      !!call communication_arrays_repartitionrho(iproc, nproc, lzd, nscatterarr, istartend, &
      !!     collcom_sr%nsendcounts_repartitionrho, collcom_sr%nsenddspls_repartitionrho, &
      !!     collcom_sr%nrecvcounts_repartitionrho, collcom_sr%nrecvdspls_repartitionrho)
    
      call communication_arrays_repartitionrho_general(iproc, nproc, lzd, nscatterarr, istartend, & 
           collcom_sr%ncomms_repartitionrho, collcom_sr%commarr_repartitionrho)
    
      call f_free(weights_per_zpoint)
      call f_free(istartend)
    
      call timing(iproc,'init_collco_sr','OF')
    
    end subroutine init_comms_linear_sumrho



    subroutine get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, &
               weight_tot, weight_ideal, weights_per_slice, weights_per_zpoint)
      use module_base
      use module_types
      use locregs, only: check_whether_bounds_overlap
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      real(kind=8),intent(out) :: weight_tot, weight_ideal
      real(kind=8),dimension(0:nproc-1),intent(out) :: weights_per_slice
      real(kind=8),dimension(lzd%glr%d%n3i),intent(out) :: weights_per_zpoint
    
      ! Local variables
      integer :: iorb, ilr, i3, i2, i1, is1, ie1, is2, ie2, is3, ie3, js3, je3, ii1, ii2
      real(kind=8) :: tt, zz
      real(kind=8),dimension(:,:),allocatable :: weight_xy
    
      call f_routine(id='get_weights_sumrho')
    
      !call to_zero(lzd%glr%d%n3i, weights_per_zpoint(1))
      call f_zero(weights_per_zpoint)
    
      weight_xy=f_malloc((/lzd%glr%d%n1i,lzd%glr%d%n2i/),id='weight_xy')
    
      !write(*,*) 'iproc, nscatterarr', iproc, nscatterarr(iproc,:)
    
      tt=0.d0
      weights_per_slice(:) = 0.0d0
      js3=nscatterarr(iproc,3)+1
      je3=nscatterarr(iproc,3)+nscatterarr(iproc,2)
      !do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,2)
      do i3=js3,je3
         !call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i, weight_xy(1,1))
          call f_zero(weight_xy)
          do iorb=1,orbs%norb
              if (orbs%spinsgn(iorb)<0.d0) cycle !consider only up orbitals
              ilr=orbs%inwhichlocreg(iorb)
              !is3=1+lzd%Llr(ilr)%nsi3
              !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
              ie3=modulo(lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
              !write(*,'(a,3i8,l3)') 'is3,ie3,i3,overlap',is3,ie3,i3,check_whether_bounds_overlap(is3,ie3,i3,i3)
              !if ((is3>i3 .or. i3>ie3).neqv..not.check_whether_bounds_overlap(is3,ie3,i3,i3)) then
              !    write(*,'(a,4i8,l3)') 'ERROR: is3,ie3,i3,i3,overlap',is3,ie3,i3,i3,check_whether_bounds_overlap(is3,ie3,i3,i3)
              !end if
              !if (is3>i3 .or. i3>ie3) cycle
              !write(*,*) 'is3, ie3, i3, res', is3, ie3, i3, check_whether_bounds_overlap(is3,ie3,i3,i3)
              if (.not.check_whether_bounds_overlap(is3,ie3,i3,i3)) cycle
              !is1=1+lzd%Llr(ilr)%nsi1
              !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              !is2=1+lzd%Llr(ilr)%nsi2
              !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
              !ie1=modulo(lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i-1,lzd%glr%d%n1i)+1
              ie1=is1+lzd%llr(ilr)%d%n1i-1
              is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
              !ie2=modulo(lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i-1,lzd%glr%d%n2i)+1
              ie2=is2+lzd%llr(ilr)%d%n2i-1
              !write(*,*) 'is2,ie2,n2i_l,n2i_g',is2,ie2,lzd%llr(ilr)%d%n2i,lzd%glr%d%n2i
              !!if (1+lzd%Llr(ilr)%nsi1/=is1) stop
              !!if (1+lzd%Llr(ilr)%nsi2/=is2) stop
              !!if (lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i/=ie1) stop
              !!if (lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i/=ie2) stop
              !!$omp parallel default(none) shared(is2, ie2, is1, ie1, weight_xy, lzd) private(i2, i1, ii2, ii1)
              !!$omp do
              do i2=is2,ie2
                  ii2=modulo(i2-1,lzd%glr%d%n2i)+1
                  !write(*,*) 'i3,i2,ii2',i3,i2,ii2
                  do i1=is1,ie1
                      ii1=modulo(i1-1,lzd%glr%d%n1i)+1
                      !weight_xy(i1,i2) = weight_xy(i1,i2)+1.d0
                      weight_xy(ii1,ii2) = weight_xy(ii1,ii2)+1.d0
                  end do
              end do
              !!$omp end do
              !!$omp end parallel
          end do
          zz=0.d0
          !$omp parallel default(none) shared(lzd, weight_xy, zz, tt) private(i2, i1)
          !$omp do reduction(+: tt, zz)
          do i2=1,lzd%glr%d%n2i
              do i1=1,lzd%glr%d%n1i
                 tt = tt + .5d0*(weight_xy(i1,i2)*(weight_xy(i1,i2)+1.d0))
                 zz = zz + .5d0*(weight_xy(i1,i2)*(weight_xy(i1,i2)))
              end do
          end do
          !$omp end do
          !$omp end parallel
          weights_per_zpoint(i3)=zz
      end do
      weights_per_slice(iproc)=tt
      if (nproc > 1) then
         call mpiallred(weights_per_slice(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
         call mpiallred(tt,1,mpi_sum, comm=bigdft_mpi%mpi_comm,recvbuf=weight_tot)
         !call mpi_allreduce(tt, weight_tot, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
         call mpiallred(weights_per_zpoint(1), lzd%glr%d%n3i, mpi_sum, comm=bigdft_mpi%mpi_comm)
      else
         weight_tot=tt
      end if

      call f_free(weight_xy)
    
      ! Ideal weight per process
      weight_ideal = weight_tot/dble(nproc)
    
      call f_release_routine()
    
    end subroutine get_weights_sumrho


    subroutine assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, weights_per_slice, &
               lzd, orbs, nscatterarr, istartend, nptsp)
      use module_base
      use module_types
      use locregs, only: check_whether_bounds_overlap
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      real(kind=8),intent(in) :: weight_tot, weight_ideal
      real(kind=8),dimension(0:nproc-1),intent(in) :: weights_per_slice
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer(kind=8),dimension(2,0:nproc-1),intent(out) :: istartend
      integer,intent(out) :: nptsp
    
      ! Local variables
      integer :: jproc, i1, i2, i3, iorb, ilr, is1, ie1, is2, ie2, is3, ie3, jproc_out, iii, ierr
      integer,dimension(:),allocatable :: recvcounts, displs
      real(kind=8),dimension(:,:),allocatable :: slicearr
      real(kind=8), dimension(:,:),allocatable :: weights_startend
      real(kind=8) :: tt
      integer(kind=8) :: ii, sendbuf
    
      call f_routine(id='assign_weight_to_process_sumrho')
    
      weights_startend=f_malloc((/1.to.2,0.to.nproc-1/),id='weights_startend')
    
      tt=0.d0
      weights_startend(1,0)=0.d0
      do jproc=0,nproc-2
          tt=tt+weight_ideal
          weights_startend(2,jproc)=dble(floor(tt,kind=8))
          weights_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
      end do
      weights_startend(2,nproc-1)=weight_tot
    
      ! Iterate through all grid points and assign them to processes such that the
      ! load balancing is optimal.
      if (nproc==1) then
          istartend(1,0)=int(1,kind=8)
          istartend(2,0)=int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n3i,kind=8)
      else
          slicearr=f_malloc((/lzd%glr%d%n1i,lzd%glr%d%n2i/),id='slicearr')
          istartend(1,:)=0
          istartend(2,:)=0
          tt=0.d0
          jproc=0
          ii=int(0,kind=8)
          outer_loop: do jproc_out=0,nproc-1
              if (tt+weights_per_slice(jproc_out)<weights_startend(1,iproc)) then
                  tt=tt+weights_per_slice(jproc_out)
                  ii=ii+int(nscatterarr(jproc_out,2),kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)
                  cycle outer_loop
              end if
              i3_loop: do i3=nscatterarr(jproc_out,3)+1,nscatterarr(jproc_out,3)+nscatterarr(jproc_out,2)
                 !call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i, slicearr(1,1))
                  call f_zero(slicearr)
                  do iorb=1,orbs%norb
                      if (orbs%spinsgn(iorb)<0.d0) cycle !consider only up orbitals
                      ilr=orbs%inwhichlocreg(iorb)
                      !is1=1+lzd%Llr(ilr)%nsi1
                      !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                      !is2=1+lzd%Llr(ilr)%nsi2
                      !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                      !is3=1+lzd%Llr(ilr)%nsi3
                      !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                      is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
                      ie1=modulo(lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i-1,lzd%glr%d%n1i)+1
                      is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
                      ie2=modulo(lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i-1,lzd%glr%d%n2i)+1
                      is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
                      ie3=modulo(lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
                      !if (is3>i3 .or. i3>ie3) cycle
                      if (.not.check_whether_bounds_overlap(is3,ie3,i3,i3)) cycle
                      !$omp parallel default(none) shared(lzd, slicearr, is1, ie1, is2, ie2) private(i1, i2)
                      !$omp do
                      do i2=1,lzd%glr%d%n2i
                          !if (is2>i2 .or. i2>ie2) cycle
                          if (.not.check_whether_bounds_overlap(is2,ie2,i2,i2)) cycle
                          do i1=1,lzd%glr%d%n1i
                              !if (is1<=i1 .and. i1<=ie1) then
                              if (check_whether_bounds_overlap(is1,ie1,i1,i1)) then
                                  slicearr(i1,i2)=slicearr(i1,i2)+1.d0
                              end if
                          end do
                      end do
                      !$omp end do
                      !$omp end parallel
                  end do
                  do i2=1,lzd%glr%d%n2i
                      do i1=1,lzd%glr%d%n1i
                          ii=ii+int(1,kind=8)
                          tt=tt+.5d0*slicearr(i1,i2)*(slicearr(i1,i2)+1.d0)
                          if (tt>=weights_startend(1,iproc)) then
                              istartend(1,iproc)=ii
                              exit outer_loop
                          end if
                      end do
                   end do
               end do i3_loop
            end do outer_loop
            call f_free(slicearr)
      end if
    
        
      if (nproc > 1) then
          ! call mpiallred(istartend(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm)
          recvcounts = f_malloc(0.to.nproc-1,id='recvcounts')
          displs = f_malloc(0.to.nproc-1,id='displs')
          do jproc=0,nproc-1
              recvcounts(jproc) = 1
              displs(jproc) = 2*jproc
          end do
          sendbuf = istartend(1,iproc)
          call mpi_allgatherv(sendbuf, 1, mpi_integer8, istartend(1,0), &
               recvcounts, displs, mpi_integer8, bigdft_mpi%mpi_comm, ierr)
          call f_free(recvcounts)
          call f_free(displs)
      end if
    
      do jproc=0,nproc-2
          istartend(2,jproc)=istartend(1,jproc+1)-1
      end do
      istartend(2,nproc-1)=int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n3i,kind=8)
    
      nptsp = int(istartend(2,iproc)-istartend(1,iproc),kind=4) + 1

      !!write(*,*) 'iproc, npts', iproc, lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
      !!write(*,*) 'iproc, istartend', iproc, istartend
      !!write(*,*) 'weight_tot', weight_tot
    
      call f_free(weights_startend)
    
      ! Some check
      tt=real(nptsp,kind=8)
      if (nproc > 1) then
        call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (tt/=real(lzd%glr%d%n1i,kind=8)*real(lzd%glr%d%n2i,kind=8)*real(lzd%glr%d%n3i,kind=8)) then
          write(*,'(a,2es24.14)') 'tt, real(lzd%glr%d%n1i,kind=8)*real(lzd%glr%d%n2i,kind=8)*real(lzd%glr%d%n3i,kind=8)', &
                      tt, real(lzd%glr%d%n1i,kind=8)*real(lzd%glr%d%n2i,kind=8)*real(lzd%glr%d%n3i,kind=8)
          stop 'tt/=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i'
      end if
    
      call f_release_routine()
    
    end subroutine assign_weight_to_process_sumrho



    subroutine determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, nptsp, lzd, orbs, &
               istartend, weight_tot, weights_per_zpoint, norb_per_gridpoint)
      use module_base
      use module_types
      use locregs, only: check_whether_bounds_overlap
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nptsp
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer(kind=8),dimension(2,0:nproc-1),intent(in) :: istartend
      real(kind=8),intent(in) :: weight_tot
      real(kind=8),dimension(lzd%glr%d%n3i),intent(in) :: weights_per_zpoint
      integer,dimension(nptsp),intent(out) :: norb_per_gridpoint
    
      ! Local variables
      integer :: i3, i2, i1, ipt, ilr, is1, ie1, is2, ie2, is3, ie3, iorb, i, jproc, j1, j2
      real(kind=8) :: tt, weight_check
      integer(kind=8) :: ii, ii2, ii3    
    
      if (nptsp>0) then
         !call to_zero(nptsp, norb_per_gridpoint(1))
          call f_zero(norb_per_gridpoint)
      end if
      do i3=1,lzd%glr%d%n3i
          !if (iproc==0) write(*,'(a,4i12)') 'v1, b1, v2, b2', int(i3,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8),istartend(1,iproc), &
          !        int(i3-1,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8),istartend(2,iproc)
          if (int(i3,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)<istartend(1,iproc) .or. &
              int(i3-1,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)+int(1,kind=8)>istartend(2,iproc)) then
              !if (iproc==0) write(*,'(a,i0)') 'cycle for i3=',i3
              cycle
          end if
          ii3=int(i3-1,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)
          if (weights_per_zpoint(i3)==0.d0) then
              cycle
          end if
          do iorb=1,orbs%norbu
              ilr=orbs%inwhichlocreg(iorb)
              !is3=1+lzd%Llr(ilr)%nsi3
              !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
              ie3=modulo(lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
              !if (is3>i3 .or. i3>ie3) cycle
              if (.not.check_whether_bounds_overlap(is3,ie3,i3,i3)) cycle
              !is2=1+lzd%Llr(ilr)%nsi2
              !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              !is1=1+lzd%Llr(ilr)%nsi1
              !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
              !ie2=modulo(lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i-1,lzd%glr%d%n2i)+1
              ie2=is2+lzd%llr(ilr)%d%n2i-1
              is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
              !ie1=modulo(lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i-1,lzd%glr%d%n1i)+1
              ie1=is1+lzd%llr(ilr)%d%n1i-1
              !$omp parallel default(none) &
              !$omp shared(i3, ii3, is2, ie2, is1, ie1, lzd, istartend, iproc, norb_per_gridpoint) &
              !$omp private(j2, j1, i2, i1, ii, ii2, ipt)
              !$omp do
              do i2=is2,ie2
                  j2=modulo(i2-1,lzd%glr%d%n2i)+1
                  !ii2=ii3+int(i2-1,kind=8)*int(lzd%glr%d%n1i,kind=8)
                  ii2=ii3+int(j2-1,kind=8)*int(lzd%glr%d%n1i,kind=8)
                  do i1=is1,ie1
                      j1=modulo(i1-1,lzd%glr%d%n1i)+1
                      ii=ii2+int(j1,kind=8)
                      if (ii>=istartend(1,iproc) .and. ii<=istartend(2,iproc)) then
                          ipt=int(ii-istartend(1,iproc),kind=4)+1
                          !write(1000+iproc,'(a,5i9)') 'i1, i2, i3, ipt, npg',i1, i2, i3, ipt, norb_per_gridpoint(ipt)
                          norb_per_gridpoint(ipt)=norb_per_gridpoint(ipt)+1
                      end if
                  end do
              end do
              !$omp end do
              !$omp end parallel
          end do
      end do
      !call mpi_finalize(i)
      !stop
    
      tt=0.d0
      !$omp parallel default(none) shared(tt, nptsp, norb_per_gridpoint) private(i)
      !$omp do reduction(+:tt)
      do i=1,nptsp
          tt=tt+.5d0*dble(norb_per_gridpoint(i)*(norb_per_gridpoint(i)+1))
      end do
      !$omp end do
      !$omp end parallel
      weight_check=tt
    
      ! Some check
      if (nproc > 1) then
        call mpiallred(weight_check, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (abs(weight_check-weight_tot) > 1.d-3) then
          write(*,*) 'ERROR: weight_check/=weight_tot', weight_check, weight_tot
          stop '2: weight_check/=weight_tot'
      else if (abs(weight_check-weight_tot) > 0.d0) then
         call yaml_warning('The total weight for density seems inconsistent! Ref:'//&
               trim(yaml_toa(weight_tot,fmt='(1pe25.17)'))//', Check:'//&
               trim(yaml_toa(weight_check,fmt='(1pe25.17)')))
      end if
    
    end subroutine determine_num_orbs_per_gridpoint_sumrho


    subroutine determine_communication_arrays_sumrho(iproc, nproc, nptsp, lzd, orbs, &
               istartend, nsendcounts, nsenddspls, nrecvcounts, &
               nrecvdspls, ndimpsi)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nptsp
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer(kind=8),dimension(2,0:nproc-1),intent(in) :: istartend
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      integer,intent(out) :: ndimpsi
    
      ! Local variables
      integer :: iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, jproc, i3, i2, i1, ii, ierr, ii0, j1, j2, j3
      integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
      character(len=*),parameter :: subname='determine_communication_arrays_sumrho'
      integer(kind=8) :: ind, ii2, ii3
    
    
      call f_zero(nsendcounts)
     
    
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          !is1=1+lzd%Llr(ilr)%nsi1
          !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
          !is2=1+lzd%Llr(ilr)%nsi2
          !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
          !is3=1+lzd%Llr(ilr)%nsi3
          !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
          is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
          !ie1=modulo(lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i-1,lzd%glr%d%n1i)+1
          ie1=is1+lzd%llr(ilr)%d%n1i-1
          is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
          !ie2=modulo(lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i-1,lzd%glr%d%n2i)+1
          ie2=is2+lzd%llr(ilr)%d%n2i-1
          is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
          !ie3=modulo(lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
          ie3=is3+lzd%llr(ilr)%d%n3i-1
          do jproc=0,nproc-1
              ii=0
              do i3=is3,ie3
                  j3=modulo(i3-1,lzd%glr%d%n3i)+1
                  if (int(j3,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)<istartend(1,jproc) .or. &
                      int(j3-1,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)+int(1,kind=8)>istartend(2,jproc)) then
                      cycle
                  end if
                  ii0=0
                  ii3=int(j3-1,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)
                  !$omp parallel default(none) &
                  !$omp shared(is2, ie2, is1, ie1, lzd, istartend, jproc, ii0, ii3) &
                  !$omp private(i2, i1, ind, ii2, j1, j2)
                  !$omp do reduction(+:ii0)
                  do i2=is2,ie2
                      j2=modulo(i2-1,lzd%glr%d%n2i)+1
                      ii2=ii3+int(j2-1,kind=8)*int(lzd%glr%d%n1i,kind=8)
                      do i1=is1,ie1
                        j1=modulo(i1-1,lzd%glr%d%n1i)+1
                        ind = ii2+int(j1,kind=8)
                        if (ind>=istartend(1,jproc) .and. ind<=istartend(2,jproc)) then
                            !nsendcounts(jproc)=nsendcounts(jproc)+1
                            ii0=ii0+1
                        end if
                      end do
                  end do
                  !$omp end do
                  !$omp end parallel
                  ii=ii+ii0
              end do
             nsendcounts(jproc)=nsendcounts(jproc)+ii
           end do
      end do
    
    
      ! Some check
      ii=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          ii = ii + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
      end do
      if (ii/=sum(nsendcounts)) then
          call f_err_throw('Error in determine_communication_arrays_sumrho: ii/=sum(nsendcounts); values are'//&
               trim(yaml_toa(ii))//' and '//&
               trim(yaml_toa(sum(nsendcounts))),&
               err_name='BIGDFT_RUNTIME_ERROR')
      end if
      ndimpsi=ii
    
    
      nsenddspls(0)=0
      do jproc=1,nproc-1
          nsenddspls(jproc)=nsenddspls(jproc-1)+nsendcounts(jproc-1)
      end do
    
      nsendcounts_tmp = f_malloc(0.to.nproc-1,id='nsendcounts_tmp')
      nsenddspls_tmp = f_malloc(0.to.nproc-1,id='nsenddspls_tmp')
      nrecvcounts_tmp = f_malloc(0.to.nproc-1,id='nrecvcounts_tmp')
      nrecvdspls_tmp = f_malloc(0.to.nproc-1,id='nrecvdspls_tmp')
      nsendcounts_tmp=1
      nrecvcounts_tmp=1
      do jproc=0,nproc-1
          nsenddspls_tmp(jproc)=jproc
          nrecvdspls_tmp(jproc)=jproc
      end do
      if(nproc>1) then
          call mpialltoallv(nsendcounts(0), nsendcounts_tmp, nsenddspls_tmp, &
               nrecvcounts(0), nrecvcounts_tmp, nrecvdspls_tmp, bigdft_mpi%mpi_comm)
      else
          nrecvcounts=nsendcounts
      end if
      call f_free(nsendcounts_tmp)
      call f_free(nsenddspls_tmp)
      call f_free(nrecvcounts_tmp)
      call f_free(nrecvdspls_tmp)
    
      !!ndimind = sum(nrecvcounts)
    
      !!! Some check
      !!ii=sum(norb_per_gridpoint)
      !!if (ii/=ndimind) stop 'ii/=sum(nrecvcounts)'
    
      nrecvdspls(0)=0
      do jproc=1,nproc-1
          nrecvdspls(jproc)=nrecvdspls(jproc-1)+nrecvcounts(jproc-1)
      end do
    
    end subroutine determine_communication_arrays_sumrho



    subroutine get_switch_indices_sumrho(iproc, nproc, nptsp, ndimpsi, ndimind, lzd, orbs, nspin, istartend, &
               norb_per_gridpoint, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, &
               isendbuf, irecvbuf, iextract, iexpand, indexrecvorbital)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nptsp, ndimpsi, ndimind, nspin
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer(kind=8),dimension(2,0:nproc-1),intent(in) :: istartend
      integer,dimension(nptsp),intent(in) :: norb_per_gridpoint
      integer,dimension(0:nproc-1),intent(in) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      integer,dimension(ndimpsi),intent(out) :: isendbuf, irecvbuf
      integer,dimension(ndimind),intent(out) :: iextract, iexpand, indexrecvorbital
    
      ! Local variables
      integer :: jproc, iitot, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, ind, ierr, ii, j1, j2, j3
      integer :: iorb, i, ipt, itotadd
      integer,dimension(:),allocatable :: nsend, indexsendorbital, indexsendorbital2, indexrecvorbital2
      integer,dimension(:),allocatable :: gridpoint_start, gridpoint_start_tmp
      character(len=*),parameter :: subname='get_switch_indices_sumrho'
      integer(kind=8) :: indglob3a, indglob3, indglob2, indglob
      integer(kind=8),dimension(:),allocatable :: indexsendbuf, indexrecvbuf
      integer(kind=8) :: iilong, ilong
    
    
      nsend = f_malloc(0.to.nproc-1,id='nsend')
      nsend=0
      indexsendbuf = f_malloc(max(1,ndimpsi),id='indexsendbuf')
      indexsendorbital = f_malloc(max(1,ndimpsi),id='indexsendorbital')
      !!allocate(isendbuf(ndimpsi), stat=istat)
      !!call memocc(istat, isendbuf, 'isendbuf', subname)
    
      iitot=0
      !!$omp parallel default(shared) &
      !!$omp private(iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, indglob, ind)
      !!$omp do lastprivate(iitot)
      do jproc=0,nproc-1
          iitot=0
          do iorb=1,orbs%norbp
              iiorb=orbs%isorb+iorb
              ilr=orbs%inwhichlocreg(iiorb)
              !is1=1+lzd%Llr(ilr)%nsi1
              !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              !is2=1+lzd%Llr(ilr)%nsi2
              !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              !is3=1+lzd%Llr(ilr)%nsi3
              !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
              !ie1=modulo(lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i-1,lzd%glr%d%n1i)+1
              ie1=is1+lzd%llr(ilr)%d%n1i-1
              is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
              !ie2=modulo(lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i-1,lzd%glr%d%n2i)+1
              ie2=is2+lzd%llr(ilr)%d%n2i-1
              is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
              !ie3=modulo(lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
              ie3=is3+lzd%llr(ilr)%d%n3i-1
              itotadd=(ie2-is2+1)*(ie1-is1+1)
              do i3=is3,ie3
                  j3=modulo(i3-1,lzd%glr%d%n3i)+1
                  indglob3a=int(j3,kind=8)*int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)
                  indglob3=indglob3a-int(lzd%glr%d%n1i,kind=8)*int(lzd%glr%d%n2i,kind=8)
                  if (indglob3a<istartend(1,jproc) .or. &
                      indglob3+int(1,kind=8)>istartend(2,jproc)) then
                      iitot=iitot+itotadd
                      cycle
                  end if
                  do i2=is2,ie2
                      j2=modulo(i2-1,lzd%glr%d%n2i)+1
                      indglob2=indglob3+int(j2-1,kind=8)*int(lzd%glr%d%n1i,kind=8)
                      do i1=is1,ie1
                          j1=modulo(i1-1,lzd%glr%d%n1i)+1
                          indglob = indglob2+int(j1,kind=8)
                          iitot=iitot+1
                          if (indglob>=istartend(1,jproc) .and. indglob<=istartend(2,jproc)) then
                              nsend(jproc)=nsend(jproc)+1
                              ind=nsenddspls(jproc)+nsend(jproc)
                              isendbuf(iitot)=ind
                              indexsendbuf(ind)=indglob
                              indexsendorbital(iitot)=iiorb
                              !exit
                          end if
                      end do
                  end do
              end do
          end do
      end do
      !!$omp end do
      !!$omp end parallel
    
    
      if(iitot/=ndimpsi) stop 'iitot/=ndimpsi'
    
      !check
      do jproc=0,nproc-1
          if(nsend(jproc)/=nsendcounts(jproc)) stop 'nsend(jproc)/=nsendcounts(jproc)'
      end do
    
    !!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    !!t2=mpi_wtime()
    !!tt=t2-t1
    !!if(iproc==0) write(*,*) 'time 5.1: iproc', iproc, tt
    
    
    
      !!allocate(irecvbuf(ndimpsi), stat=istat)
      !!call memocc(istat, irecvbuf, 'irecvbuf', subname)
    
      indexsendorbital2 = f_malloc(max(1,ndimpsi),id='indexsendorbital2')
      call vcopy(ndimpsi, indexsendorbital(1), 1, indexsendorbital2(1), 1)
      do i=1,ndimpsi
          ind=isendbuf(i)
          indexsendorbital(ind)=indexsendorbital2(i)
      end do
    
      ! Inverse of isendbuf
      call get_reverse_indices(ndimpsi, isendbuf, irecvbuf)
    
      call f_free(indexsendorbital2)
    
    
      indexrecvbuf = f_malloc(ndimind,id='indexrecvbuf')
      !!allocate(indexrecvorbital(ndimind), stat=istat)
      !!call memocc(istat, indexrecvorbital, 'indexrecvorbital', subname)
    
      if(nproc>1) then
          ! Communicate indexsendbuf
          call mpialltoallv(indexsendbuf(1), nsendcounts, nsenddspls, &
               indexrecvbuf(1), nrecvcounts, nrecvdspls, bigdft_mpi%mpi_comm)
          ! Communicate indexsendorbitals
          call mpialltoallv(indexsendorbital(1), nsendcounts, nsenddspls, &
               indexrecvorbital(1), nrecvcounts, nrecvdspls, bigdft_mpi%mpi_comm)
       else
           indexrecvbuf=indexsendbuf
           indexrecvorbital=indexsendorbital
       end if
    
      call f_free(indexsendbuf)
    
      call f_free(indexsendorbital)
    !!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    !!t2=mpi_wtime()
    !!tt=t2-t1
    !!if(iproc==0) write(*,*) 'time 5.2: iproc', iproc, tt
    
    
       gridpoint_start = f_malloc(0.to.int(istartend(2,iproc)-istartend(1,iproc),kind=4),id='gridpoint_start')
       gridpoint_start_tmp = f_malloc(0.to.int(istartend(2,iproc)-istartend(1,iproc),kind=4),id='gridpoint_start_tmp')
    
       ii=1
       do ipt=1,nptsp
           ilong=int(ipt,kind=8)+istartend(1,iproc)-int(1,kind=8)
           if (norb_per_gridpoint(ipt)>0) then
               gridpoint_start(ilong-istartend(1,iproc))=ii
           else
               gridpoint_start(ilong-istartend(1,iproc))=0
           end if
           ii=ii+norb_per_gridpoint(ipt)
       end do
    
       if (nspin*ii/=ndimind+nspin) then
           stop '(nspin*ii/=ndimind+nspin)'
       end if
       if(maxval(gridpoint_start)>ndimind) stop '1: maxval(gridpoint_start)>sum(nrecvcountc)'
       if(minval(indexrecvbuf)<istartend(1,iproc)) stop '1: minval(indexrecvbuf)<istartend(1,iproc)'
       if(maxval(indexrecvbuf)>istartend(2,iproc)) stop '1: maxval(indexrecvbuf)>istartend(2,iproc)'
    
       !!allocate(iextract(ndimind), stat=istat)
       !!call memocc(istat, iextract, 'iextract', subname)

       gridpoint_start_tmp = gridpoint_start
    
      ! Rearrange the communicated data
      do i=1,ndimind
          iilong=indexrecvbuf(i)
          ind=gridpoint_start(iilong-istartend(1,iproc))
          if (gridpoint_start(iilong-istartend(1,iproc))-gridpoint_start_tmp(iilong-istartend(1,iproc))+1 > &
                norb_per_gridpoint(int(iilong-istartend(1,iproc),kind=4)+1)) then
              ! orbitals which fulfill this condition are down orbitals which
              ! should be put at the end
              ind = ind + (ndimind/2-norb_per_gridpoint(int(iilong-istartend(1,iproc),kind=8)+1))
          end if
          iextract(i)=ind
          gridpoint_start(iilong-istartend(1,iproc))=gridpoint_start(iilong-istartend(1,iproc))+1
      end do

    
      if(maxval(iextract)>ndimind) then
          stop 'maxval(iextract)>ndimind'
      end if
      if(minval(iextract)<1) stop 'minval(iextract)<1'
    
      call f_free(indexrecvbuf)
    
    
      !! allocate(iexpand(ndimind), stat=istat)
      !! call memocc(istat, iexpand, 'iexpand', subname)
      ! Get the array to transfrom back the data
      call get_reverse_indices(ndimind, iextract, iexpand)
    
    !!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    !!t2=mpi_wtime()
    !!tt=t2-t1
    !!if(iproc==0) write(*,*) 'time 5.3: iproc', iproc, tt
    
      indexrecvorbital2 = f_malloc(ndimind,id='indexrecvorbital2')
    
      if (ndimind>0) then
          call vcopy(ndimind, indexrecvorbital(1), 1, indexrecvorbital2(1), 1)
      end if
    
      !$omp parallel default(none) &
      !$omp shared(ndimind, iextract, indexrecvorbital, indexrecvorbital2) private(i, ind)
      !$omp do
      do i=1,ndimind
          ind=iextract(i)
          indexrecvorbital(ind)=indexrecvorbital2(i)
      end do
      !$omp end do
      !$omp end parallel
    
      call f_free(indexrecvorbital2)
    
      if(minval(indexrecvorbital)<1) stop 'minval(indexrecvorbital)<1'
      if(maxval(indexrecvorbital)>orbs%norb) stop 'maxval(indexrecvorbital)>orbs%norb'
    
    
      call f_free(gridpoint_start)
      call f_free(gridpoint_start_tmp)
      call f_free(nsend)
    
    
    end subroutine get_switch_indices_sumrho



    subroutine communication_arrays_repartitionrho(iproc, nproc, lzd, nscatterarr, istartend, &
               nsendcounts_repartitionrho, nsenddspls_repartitionrho, &
               nrecvcounts_repartitionrho, nrecvdspls_repartitionrho)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer(kind=8),dimension(2,0:nproc-1),intent(in) :: istartend
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts_repartitionrho, nsenddspls_repartitionrho
      integer,dimension(0:nproc-1),intent(out) :: nrecvcounts_repartitionrho, nrecvdspls_repartitionrho
    
      ! Local variables
      integer :: jproc_send, jproc_recv, i3, i2, i1, jproc
      integer(kind=8) :: ii
    
      jproc_send=0
      jproc_recv=0
      ii=int(0,kind=8)
      nsendcounts_repartitionrho=0
      nrecvcounts_repartitionrho=0
      do i3=1,lzd%glr%d%n3i
          do i2=1,lzd%glr%d%n2i
              do i1=1,lzd%glr%d%n1i
                  ii=ii+int(1,kind=8)
                  if (ii>istartend(2,jproc_send)) then
                      jproc_send=jproc_send+1
                  end if
                  if (i3>nscatterarr(jproc_recv,3)+nscatterarr(jproc_recv,2)) then
                      jproc_recv=jproc_recv+1
                  end if
                  if (iproc==jproc_send) then
                      nsendcounts_repartitionrho(jproc_recv)=nsendcounts_repartitionrho(jproc_recv)+1
                  end if
                  if (iproc==jproc_recv) then
                      nrecvcounts_repartitionrho(jproc_send)=nrecvcounts_repartitionrho(jproc_send)+1
                  end if
              end do
          end do
      end do
    
      nsenddspls_repartitionrho(0)=0
      nrecvdspls_repartitionrho(0)=0
      do jproc=1,nproc-1
          nsenddspls_repartitionrho(jproc)=nsenddspls_repartitionrho(jproc-1)+&
                                                      nsendcounts_repartitionrho(jproc-1)
          nrecvdspls_repartitionrho(jproc)=nrecvdspls_repartitionrho(jproc-1)+&
                                                      nrecvcounts_repartitionrho(jproc-1)
      end do
    
    end subroutine communication_arrays_repartitionrho
    
    
    subroutine communication_arrays_repartitionrho_general(iproc, nproc, lzd, nscatterarr, istartend, &
               ncomms_repartitionrho, commarr_repartitionrho)
      use module_base
      use module_types
      use locregs, only: get_extent_of_overlap
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer(kind=8),dimension(2,0:nproc-1),intent(in) :: istartend
      integer,intent(out) :: ncomms_repartitionrho
      integer,dimension(:,:),pointer,intent(out) :: commarr_repartitionrho
      character(len=*),parameter :: subname='communication_arrays_repartitionrho_general'
    
      ! Local variables
      integer :: i1, i2, i3, jproc, jproc_send, iidest, nel, ioverlaps, iassign, is3, ie3, iis3, iie3, i
      logical :: started
      integer(kind=8) :: ii, iis, iie
      integer(kind=8),dimension(2) :: iiis, iiie, nlen
      integer(kind=8) :: is
      integer :: n, j
    
      call f_routine(id='communication_arrays_repartitionrho_general')

      !!write(*,'(a,4i8,3x,6i8)') 'n1, n2, n3, ntot, istartend',lzd%glr%d%n1i,lzd%glr%d%n2i,lzd%glr%d%n3i,lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i,istartend
      !!write(*,'(a,2i8,3x,6i8)') 'iproc, n3i, nscatterarr(iproc,:)', iproc, lzd%glr%d%n3i, nscatterarr(iproc,:)
    
      ! only do this if task iproc has to receive a part of the potential
      if (nscatterarr(iproc,1)>0) then
        
          !!!! First process from which iproc has to receive data
          !!!ncomms_repartitionrho=0
          !!!i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)
          !!!ii=int(i3,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)+int(1,kind=8)
          !!!do jproc=nproc-1,0,-1
          !!!    if (ii>=istartend(1,jproc)) then
          !!!        jproc_send=jproc
          !!!        write(*,'(a,5i8)') 'FIRST: iproc, i3, ii, jproc_send', iproc, i3, ii, jproc_send
          !!!        ncomms_repartitionrho=ncomms_repartitionrho+1
          !!!        exit
          !!!    end if
          !!!end do
        
          !!!! The remaining processes
          !!!iidest=0
          !!!nel=0
          !!!started=.false.
          !!!do i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+1,nscatterarr(iproc,3)-nscatterarr(iproc,4)+nscatterarr(iproc,1)
          !!!    write(*,'(a,3i8)') 'iproc, i3, lzd%glr%d%n3i', iproc, i3, lzd%glr%d%n3i
          !!!    ii=int(i3-1,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)
          !!!    do i2=1,lzd%glr%d%n2i
          !!!        do i1=1,lzd%glr%d%n1i
          !!!            ii=ii+int(1,kind=8)
          !!!            iidest=iidest+1
          !!!            if (ii>=istartend(1,jproc_send) .and. ii<=istartend(2,jproc_send)) then
          !!!                nel=nel+1
          !!!            else
          !!!                write(*,'(a,5i8)') 'iproc, i1, i2, i3, ii, jproc_send', iproc, i1, i2, i3, ii, jproc_send
          !!!                jproc_send=jproc_send+1
          !!!                ncomms_repartitionrho=ncomms_repartitionrho+1
          !!!            end if
          !!!        end do
          !!!    end do
          !!!end do

          !@NEW #########################
          ! Starting and ending point of the density required by task iproc (in z direction)
          is3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+1
          ie3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+nscatterarr(iproc,1)
          ! Due to perdiodic boundary conditions, iie3 might be smaller than iis3
          iis3=modulo(is3-1,lzd%glr%d%n3i)+1
          iie3=modulo(ie3-1,lzd%glr%d%n3i)+1
          
          ! Starting and ending point of the density required by task iproc (in global coordinates)
          iis=int(iis3-1,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)+int(1,kind=8)
          iie=int(iie3,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)

          ! Check whether there is an overlap between the density required on
          ! task iproc and the one calculated on task jproc
          !jproc_send=0
          ncomms_repartitionrho=0
          do jproc=0,nproc-1
              !if(check_whether_bounds_overlap(iis,iie,istartend(1,jproc),istartend(2,jproc))) then
              if (istartend(2,jproc)>=istartend(1,jproc)) then
                  call get_extent_of_overlap(iis,iie,istartend(1,jproc),istartend(2,jproc), n, iiis, iiie, nlen)
              end if
              !!write(*,'(a,11i11)') 'iproc, jproc, iis, iie, ise(1), ise(2), n, iiis, iiie', iproc, jproc, iis, iie, istartend(1,jproc), istartend(2,jproc), n, iiis, iiie
              !jproc_send=jproc_send+1
              ncomms_repartitionrho=ncomms_repartitionrho+n
              !end if
          end do
          !@END NEW #####################
        
        
          call allocate_MPI_comms_cubic_repartitionp2p(ncomms_repartitionrho, commarr_repartitionrho)
        
        
          !!! First process from which iproc has to receive data
          !!ioverlaps=0
          !!i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)
          !!ii=int(i3,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)+int(1,kind=8)
          !!do jproc=nproc-1,0,-1
          !!    if (ii>=istartend(1,jproc)) then
          !!        jproc_send=jproc
          !!        ioverlaps=ioverlaps+1
          !!        exit
          !!    end if
          !!end do
        
        
          !!! The remaining processes
          !!iassign=0
          !!iidest=0
          !!nel=0
          !!started=.false.
          !!do i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+1,nscatterarr(iproc,3)-nscatterarr(iproc,4)+nscatterarr(iproc,1)
          !!    ii=int(i3-1,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)
          !!    do i2=1,lzd%glr%d%n2i
          !!        do i1=1,lzd%glr%d%n1i
          !!            ii=ii+int(1,kind=8)
          !!            iidest=iidest+1
          !!            if (ii>=istartend(1,jproc_send) .and. ii<=istartend(2,jproc_send)) then
          !!                nel=nel+1
          !!            else
          !!                commarr_repartitionrho(4,ioverlaps)=nel
          !!                jproc_send=jproc_send+1
          !!                ioverlaps=ioverlaps+1
          !!                nel=1
          !!                started=.false.
          !!            end if
          !!            if (.not.started) then
          !!                if (jproc_send>=nproc) stop 'ERROR: jproc_send>=nproc'
          !!                commarr_repartitionrho(1,ioverlaps)=jproc_send
          !!                commarr_repartitionrho(2,ioverlaps)=int(ii-istartend(1,jproc_send),kind=8)+1
          !!                commarr_repartitionrho(3,ioverlaps)=iidest
          !!                started=.true.
          !!                iassign=iassign+1
          !!            end if
          !!        end do
          !!    end do
          !!end do
          !!commarr_repartitionrho(4,ioverlaps)=nel

          !!do i=1,ncomms_repartitionrho
          !!    write(*,'(a,i5,3x,4i8)') 'FIRST: iproc, cr(1:4,i)', iproc, commarr_repartitionrho(1:4,i)
          !!end do

          !@ NEW ###############################
          ! For each overlap, get the starting, ending point and extent
          !jproc_send=0
          ioverlaps=0
          do jproc=0,nproc-1
              !if(check_whether_bounds_overlap(iis,iie,istartend(1,jproc),istartend(2,jproc))) then
              if (istartend(2,jproc)>=istartend(1,jproc)) then
                  call get_extent_of_overlap(iis,iie,istartend(1,jproc),istartend(2,jproc), n, iiis, iiie, nlen)
                  !write(*,'(a,12i8)') 'iproc, iis,iie,istartend(:,jproc), n, iiis, iiie, nlen', iproc, iis,iie,istartend(:,jproc), n, iiis, iiie, nlen
                  ! Do nothing if n==0
                  do j=1,n
                      !jproc_send=jproc_send+1
                      ioverlaps=ioverlaps+1
                      !call get_extent_of_overlap(iis,iie,istartend(1,jproc),istartend(2,jproc), n, iiis, iiie, nlen)
                      !!if (n>1) then
                      !!    write(*,*) 'WARNING: THIS WRONG AND NEEDS A FIX'
                      !!    is=minval(iiis)
                      !!else
                      !!    is=iiis(1)
                      !!end if
                      commarr_repartitionrho(1,ioverlaps)=jproc
                      commarr_repartitionrho(2,ioverlaps)=int(iiis(j)-istartend(1,jproc),kind=4)+1
                      !i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+1
                      i3=modulo(nscatterarr(iproc,3)-nscatterarr(iproc,4)+1-1,lzd%glr%d%n3i)+1
                      ! iidest is the offest to the start of the density received by iproc
                      ii = int(i3-1,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8)
                      !iidest = int(iiis(j)-int(i3-1,kind=8)*int(lzd%glr%d%n2i,kind=8)*int(lzd%glr%d%n1i,kind=8),kind=4)
                      if (iiis(j)>ii) then
                          !standard case
                          iidest = int(iiis(j)-ii,kind=4)
                      else
                          !case with periodic wrap around
                          iidest = int(iiis(j)+lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i-ii,kind=4)
                      end if
                      !!write(*,'(a,12i11)') 'iproc, jproc, iis, iie, istartend(1,jproc), istartend(2,jproc), j, iiis(j), iiie(j), i3, iidest, iisrc', &
                      !!    iproc, jproc, iis, iie, istartend(1,jproc), istartend(2,jproc), j, iiis(j), iiie(j), i3, iidest, int(iiis(j)-istartend(1,jproc),kind=4)+1
                      commarr_repartitionrho(3,ioverlaps)=iidest
                      commarr_repartitionrho(4,ioverlaps)=int(nlen(j),kind=4)
                  end do
              end if
          end do
          !@ END NEW ###########################
          if (ioverlaps/=ncomms_repartitionrho) stop 'ERROR: ioverlaps/=ncomms_repartitionrho'
          !if (iassign/=ncomms_repartitionrho) stop 'ERROR: iassign/=ncomms_repartitionrho'

          !!do i=1,ncomms_repartitionrho
          !!    write(*,'(a,i5,3x,4i8)') 'SECOND: iproc, cr(1:4,i)', iproc, commarr_repartitionrho(1:4,i)
          !!end do
        
          ! some checks
          nel=0
          !nel_array=f_malloc0(0.to.nproc-1,id='nel_array')
          do ioverlaps=1,ncomms_repartitionrho
              nel=nel+commarr_repartitionrho(4,ioverlaps)
              !ii=commarr_repartitionrho(1,ioverlaps)
              !nel_array(ii)=nel_array(ii)+commarr_repartitionrho(4,ioverlaps)
          end do
          if (nel/=nscatterarr(iproc,1)*lzd%glr%d%n2i*lzd%glr%d%n1i) then
              write(*,'(a,3i12)') 'nel, nscatterarr(iproc,1), lzd%glr%d%n2i*lzd%glr%d%n1i', &
                   nel, nscatterarr(iproc,1), lzd%glr%d%n2i*lzd%glr%d%n1i
              stop 'nel/=nscatterarr(iproc,2)*lzd%glr%d%n2i*lzd%glr%d%n1i'
          end if
          !!call mpiallred(nel_array(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          !!if (nel_array(iproc)/=istartend(2,iproc)-istartend(1,iproc)+1) then
          !!    !stop 'nel_array(iproc)/=istartend(2,iproc)-istartend(1,iproc)+1'
          !!end if
          !!call f_free(nel_array)
    
      else
          ncomms_repartitionrho=0
          call allocate_MPI_comms_cubic_repartitionp2p(1, commarr_repartitionrho)
    
      end if
    
      call f_release_routine()
    
    end subroutine communication_arrays_repartitionrho_general




    !> Potential communication
    subroutine initialize_communication_potential(iproc, nproc, nscatterarr, orbs, lzd, nspin, comgp)
      use module_base
      use module_types
      use communications_base, only: p2pComms_null, bgq
      use locregs, only: get_extent_of_overlap, check_whether_bounds_overlap
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc
      integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      type(orbitals_data),intent(in):: orbs
      type(local_zone_descriptors),intent(in):: lzd
      integer,intent(in) :: nspin
      type(p2pComms),intent(out):: comgp
      
      ! Local variables
      integer:: is1, ie1, is2, ie2, is3, ie3, ilr, ii, iorb, iiorb, jproc, kproc, istsource, is, ie, iie3j
      integer:: ioverlap, is3j, ie3j, is3k, ie3k, mpidest, istdest, ioffsetx, ioffsety, ioffsetz, iel
      integer :: is3min, ie3max, tag, ncount, ierr, nmaxoverlap, nlen, iseg, j3, ileny, ioffset, isegx
      logical :: datatype_defined
      character(len=*),parameter:: subname='initialize_communication_potential'
      integer,dimension(6) :: ise
      integer,dimension(2) :: blocklengthsx, blocklengthsy, types, xyblock_type, nblocksy
      integer(kind=mpi_address_kind),dimension(2) :: displacementsx, displacementsy
      integer(kind=mpi_address_kind) :: lb, extent
      integer :: nsegx, nsegy, xline_type, size_of_double, size_datatype, n1, n2, n3, i
      integer,dimension(2) :: iis3, iie3, iis2, iie2, iis1, iie1, nlen1, nlen2, nlen3
      !integer,dimension(:),allocatable :: derived_types


      call timing(iproc,'init_commPot  ','ON')
      
      !call nullify_p2pComms(comgp)
      comgp = p2pComms_null()
    
      !allocate(comgp%ise(6,0:nproc-1), stat=istat)
      !call memocc(istat, comgp%ise, 'comgp%ise', subname)
      !!comgp%ise = f_malloc_ptr((/1.to.6,0.to.nproc-1/),id='comgp%ise')
      
      ! Determine the bounds of the potential that we need for
      ! the orbitals on this process.
      !iiorb=0
      is1=1000000000
      ie1=-1000000000
      is2=1000000000
      ie2=-1000000000
      is3=1000000000
      ie3=-1000000000
      !do iorb=1,orbs%norbu_par(jproc,0)
      do iorb=1,orbs%norbp
          
          iiorb=orbs%isorb+iorb 
          ilr=orbs%inwhichlocreg(iiorb)
      
          is=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
          if(is < is1) then
              is1=is
          end if
          !ii=lzd%Llr(ilr)%nsi1+lzd%Llr(ilr)%d%n1i
          ie=is+lzd%llr(ilr)%d%n1i-1
          if(ie > ie1) then
              ie1=ie
          end if
      
          is=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
          if(is < is2) then
              is2=is
          end if
          !ii=lzd%Llr(ilr)%nsi2+lzd%Llr(ilr)%d%n2i
          ie=is+lzd%llr(ilr)%d%n2i-1
          if(ie > ie2) then
              ie2=ie
          end if
      
          !ii=1+lzd%Llr(ilr)%nsi3
          is=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
          if(is < is3) then
              is3=is
          end if
          !ii=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i
          ie=is+lzd%llr(ilr)%d%n3i-1
          if(ie > ie3) then
              ie3=ie
          end if

          !!write(*,'(a,7i8)') 'ilr, lnsi1, lni1, gnsi1, gni1, is1, ie1', ilr, lzd%Llr(ilr)%nsi1, lzd%llr(ilr)%d%n1i, lzd%glr%nsi1, lzd%glr%d%n1i, is1, ie1
          !!write(*,'(a,7i8)') 'ilr, lnsi3, lni3, gnsi3, gni3, is, ie', ilr, lzd%Llr(ilr)%nsi3, lzd%llr(ilr)%d%n3i, lzd%glr%nsi3, lzd%glr%d%n3i, is, ie
      
      end do
      !!write(*,'(a,i4,3x,9i6)') 'iproc, is1, ie1, n1, is2, ie2, n2, is3, ie3, n3', iproc, is1, ie1, lzd%glr%d%n1i, is2, ie2, lzd%glr%d%n2i, is3, ie3, lzd%glr%d%n3i

      ! For non-free boundary conditions the values ie1, ie2, ie3 may lie outside of the box!
      ! Make sure that the wrapped around end is smaller than the beginning
      !if (ie1>lzd%glr%d%n1i) then
      !    ie1=min(modulo(ie1-1,lzd%glr%d%n1i)+1,is1-1)
      !end if
      !if (ie2>lzd%glr%d%n2i) then
      !    ie2=min(modulo(ie2-1,lzd%glr%d%n2i)+1,is2-1)
      !end if
      !if (ie3>lzd%glr%d%n3i) then
      !    ie3=min(modulo(ie3-1,lzd%glr%d%n3i)+1,is3-1)
      !end if
      if (ie1>lzd%glr%d%n1i) then
          ie1=min(modulo(ie1-1,lzd%glr%d%n1i)+1,is1-1)
          ie1=modulo(ie1-1,lzd%glr%d%n1i)+1
      end if
      if (ie2>lzd%glr%d%n2i) then
          ie2=min(modulo(ie2-1,lzd%glr%d%n2i)+1,is2-1)
          ie2=modulo(ie2-1,lzd%glr%d%n2i)+1
      end if
      if (ie3>lzd%glr%d%n3i) then
          ie3=min(modulo(ie3-1,lzd%glr%d%n3i)+1,is3-1)
          ie3=modulo(ie3-1,lzd%glr%d%n3i)+1
      end if
      !!write(*,'(a,i4,3x,9i6)') 'AFTER: iproc, is1, ie1, n1, is2, ie2, n2, is3, ie3, n3', iproc, is1, ie1, lzd%glr%d%n1i, is2, ie2, lzd%glr%d%n2i, is3, ie3, lzd%glr%d%n3i
      if (.not.bgq) then
          ! Communicate only the essential part, i.e. a subbox of the slices
          comgp%ise(1)=is1
          comgp%ise(2)=ie1
          comgp%ise(3)=is2
          comgp%ise(4)=ie2
      else
          ! Communicate always the entire slices
          comgp%ise(1)=1
          comgp%ise(2)=lzd%glr%d%n1i
          comgp%ise(3)=1
          comgp%ise(4)=lzd%glr%d%n2i
      end if
      comgp%ise(5)=is3
      comgp%ise(6)=ie3
    
      !!write(*,'(a,i5,6i6)') 'iproc, ise', iproc, comgp%ise
    
      
      ! Determine how many slices each process receives.
      !allocate(comgp%noverlaps(0:nproc-1), stat=istat)
      !call memocc(istat, comgp%noverlaps, 'comgp%noverlaps', subname)
      !comgp%noverlaps = f_malloc_ptr(0.to.nproc-1,id='comgp%noverlaps')
      !nmaxoverlap=0
      !do jproc=0,nproc-1
          is3j=comgp%ise(5)
          ie3j=comgp%ise(6)
          if (ie3j>lzd%glr%d%n3i) then
              ! Take modulo and make sure that it stays smaller than the beginning
              iie3j=modulo(ie3j-1,lzd%glr%d%n3i)+1
              iie3j=min(iie3j,is3j-1)
          else
              iie3j=ie3j
          end if
          mpidest=iproc
          ioverlap=0
          do kproc=0,nproc-1
              is3k=nscatterarr(kproc,3)+1
              ie3k=is3k+nscatterarr(kproc,2)-1
              if (ie3k<is3k) cycle !process kproc has no planes
              !if(is3j<=ie3k .and. ie3j>=is3k) then
              !!write(*,'(a,6i8,l6)') 'iproc, is3j, ie3j, iie3j, is3k, ie3k, overlap', iproc, is3j, ie3j, iie3j, is3k, ie3k, check_whether_bounds_overlap(is3j, iie3j, is3k, ie3k)
              !!if(check_whether_bounds_overlap(is3j, iie3j, is3k, ie3k)) then
              !!    ioverlap=ioverlap+1
              !!    !if(iproc==0) write(*,'(2(a,i0),a)') 'process ',jproc,' gets potential from process ',kproc,'.' 
              !!!TAKE INTO ACCOUNT THE PERIODICITY HERE
              !!!else if(ie3j > lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F') then
              !!!    stop 'periodicity here deprecated'
              !!!    ie3j = comgp%ise(6) - lzd%Glr%d%n3i
              !!!    if(ie3j>=is3k) then
              !!!       ioverlap=ioverlap+1
              !!!    end if
              !!!    if(is3j <= ie3k)then
              !!!       ioverlap=ioverlap+1
              !!!    end if
              !!end if
              call get_extent_of_overlap(is3j, iie3j, is3k, ie3k, n3, iis3, iie3, nlen3)
              ioverlap=ioverlap+n3
          end do
          !if (ioverlap>nmaxoverlap) nmaxoverlap=ioverlap
          comgp%noverlaps=ioverlap
          !!if(iproc==0) write(*,'(2(a,i0),a)') 'Process ',jproc,' gets ',ioverlap,' potential slices.'
      !end do
      
      ! Determine the parameters for the communications.
      !allocate(comgp%comarr(6,nmaxoverlap,0:nproc-1))
      !call memocc(istat, comgp%comarr, 'comgp%comarr', subname)
      !call to_zero(6*nmaxoverlap*nproc, comgp%comarr(1,1,0))
      comgp%comarr = f_malloc0_ptr((/1.to.6,1.to.comgp%noverlaps/),id='comgp%comarr')
      !allocate(comgp%mpi_datatypes(0:nmaxoverlap,0:nproc-1), stat=istat)
      !call memocc(istat, comgp%mpi_datatypes, 'comgp%mpi_datatypes', subname)
      !call to_zero((nmaxoverlap+1)*nproc, comgp%mpi_datatypes(0,0))
      comgp%mpi_datatypes = f_malloc0_ptr(0.to.comgp%noverlaps,id='comgp%mpi_datatypes')
      comgp%onedtypearr = f_malloc_ptr((/2,2*lzd%glr%d%n2i/),id='comgp%onedtypearr')
      comgp%nrecvBuf = 0
      is3min=0
      ie3max=0
    
      ! Only do this if we have more than one MPI task
      !nproc_if: if (nproc>1) then
          is3j=comgp%ise(5)
          ie3j=comgp%ise(6)
          if (ie3j>lzd%glr%d%n3i) then
              ! Take modulo and make sure that it stays smaller than the beginning
              iie3j=modulo(ie3j-1,lzd%glr%d%n3i)+1
              iie3j=min(iie3j,is3j-1)
          else
              iie3j=ie3j
          end if
          mpidest=iproc
          ioverlap=0
          istdest=1
          datatype_defined =.false.
          do kproc=0,nproc-1
              is3k=nscatterarr(kproc,3)+1
              ie3k=is3k+nscatterarr(kproc,2)-1
              if (ie3k<is3k) cycle !process kproc has no planes
              !SHOULD TAKE INTO ACCOUNT THE PERIODICITY HERE
              !Need to split the region
              !if(is3j<=ie3k .and. ie3j>=is3k) then
              if(check_whether_bounds_overlap(is3j, iie3j, is3k, ie3k)) then
                  !is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
                  !ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
                  call get_extent_of_overlap(is3j, iie3j, is3k, ie3k, n3, iis3, iie3, nlen3)
                  !!write(*,'(a,7i7)') 'iproc, kproc, is3j, iie3j, is3k, ie3k, n3', iproc, kproc, is3j, iie3j, is3k, ie3k, n3
                  do j3=1,n3
                      !!write(*,'(a,8i8)') 'iproc, kproc, is3j, iie3j, is3k, ie3k, iis3(j3), iie3(j3)', iproc, kproc, is3j, iie3j, is3k, ie3k, iis3(j3), iie3(j3)
                      ioffsetz=iis3(j3)-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
                      if (comgp%ise(4)>lzd%glr%d%n2i) then
                          ! Take modulo and make sure that it stays smaller than the beginning
                          ii=modulo(comgp%ise(4)-1,lzd%glr%d%n2i)+1
                          ii=min(ii,comgp%ise(3)-1)
                      else
                          ii=comgp%ise(4)
                      end if
                      call get_extent_of_overlap(comgp%ise(3), ii, 1, lzd%glr%d%n2i, n2, iis2, iie2, nlen2)
                      if (n2>1) then
                          is2=minval(iis2)
                          ie2=maxval(iie2)
                      end if
                      !if (n2/=1) stop 'n2/=1'
                      ioffsety = is2-1
                      if (comgp%ise(2)>lzd%glr%d%n1i) then
                          ! Take modulo and make sure that it stays smaller than the beginning
                          ii=modulo(comgp%ise(2)-1,lzd%glr%d%n1i)+1
                          ii=min(ii,comgp%ise(1)-1)
                      else
                          ii=comgp%ise(2)
                      end if
                      !write(*,*) 'iproc, comgp%ise(1), comgp%ise(2)', iproc, comgp%ise(1), comgp%ise(2)
                      call get_extent_of_overlap(comgp%ise(1), ii, 1, lzd%glr%d%n1i, n1, iis1, iie1, nlen1)
                      if (n1>1) then
                          is1=minval(iis1)
                          ie1=maxval(iie1)
                      end if
                      !if (n1/=1) stop 'n1/=1'
                      ioffsetx = iis1(1)
                      !if (comgp%ise(4)>lzd%glr%d%n2i) then
                      !    ! Take modulo and make sure that it stays smaller than the beginning
                      !    ioffsety=modulo(comgp%ise(4)-1,lzd%glr%d%n2i)+1
                      !    ioffsety=min(ioffsety,comgp%ise(3))-1
                      !else
                      !    ioffsety=comgp%ise(3)-1
                      !end if
                      !!ioffsetx=comgp%ise(1)
                      !if (comgp%ise(2)>lzd%glr%d%n1i) then
                      !    ! Take modulo and make sure that it stays smaller than the beginning
                      !    ioffsetx=modulo(comgp%ise(2)-1,lzd%glr%d%n1i)+1
                      !    ioffsetx=min(ioffsetx,comgp%ise(1))
                      !else
                      !    ioffsetx=comgp%ise(1)-1
                      !end if
                      !!write(*,'(a,6i8)') 'iproc, ioffsetx, ie1, ioffsety, ie2, ioffsetz', iproc, ioffsetx, ie1, ioffsety, ie2, ioffsetz
                      ioverlap=ioverlap+1

                      ! Check whether there are holes in the slices
                      if (comgp%ise(2)>lzd%glr%d%n1i) then
                          ! Take modulo and make sure that it stays smaller than the beginning
                          ii=modulo(comgp%ise(2)-1,lzd%glr%d%n1i)+1
                          ii=min(ii,comgp%ise(1)-1)
                      else
                          ii=comgp%ise(2)
                      end if
                      !if (comgp%ise(1)>is1 .and. ii<ie1) then
                      if (nproc>1) then
                          call mpi_type_size(mpi_double_precision, size_of_double, ierr)
                      else
                          size_of_double = 8
                      end if
                      !!write(*,'(a,5i8)') 'ii, is1, ie1, comgp%ise(1:2)', ii, is1, ie1, comgp%ise(1:2)
                      if (ii<comgp%ise(1) .and. ii>=is1 .and. comgp%ise(1)<ie1) then
                          !!write(*,'(a,5i8)') 'hole in x, iproc, is1, ie1, comgp%ise(1), ii', iproc, is1, ie1, comgp%ise(1), ii
                          nsegx=2
                          !!blocklengthsx(1)=comgp%ise(1)-is1+1
                          !!blocklengthsx(2)=ie1-ii+1
                          !!displacementsx(1)=int(0,kind=mpi_address_kind)
                          !!displacementsx(2)=int(ii-1,kind=mpi_address_kind)
                          blocklengthsx(1)=ii-is1+1
                          blocklengthsx(2)=ie1-comgp%ise(1)+1
                          displacementsx(1)=int(0*size_of_double,kind=mpi_address_kind)
                          displacementsx(2)=int((comgp%ise(1)-1)*size_of_double,kind=mpi_address_kind)
                      else
                          nsegx=1
                          blocklengthsx(1)=comgp%ise(2)-comgp%ise(1)+1
                          blocklengthsx(2)=0
                          displacementsx(1)=int(0*size_of_double,kind=mpi_address_kind)
                          displacementsx(2)=int(comgp%ise(2)*size_of_double,kind=mpi_address_kind)
                      end if
                      if (comgp%ise(4)>lzd%glr%d%n2i) then
                          ! Take modulo and make sure that it stays smaller than the beginning
                          ii=modulo(comgp%ise(4)-1,lzd%glr%d%n2i)+1
                          ii=min(ii,comgp%ise(3)-1)
                      else
                          ii=comgp%ise(4)
                      end if
                      !!write(*,'(a,6i8)') 'iproc, ii, is2, ie2, comgp%ise(3:4)', iproc, ii, is2, ie2, comgp%ise(3:4)
                      !if (comgp%ise(3)>is2 .and. ii<ie2) then
                      if (ii<comgp%ise(3) .and. ii>=is2 .and. comgp%ise(3)<ie2) then
                          !!write(*,*) 'iproc, hole in y', iproc
                          nsegy=2
                          !!blocklengthsy(1)=comgp%ise(3)-is2+1
                          !!blocklengthsy(2)=ie2-ii+1
                          !!displacementsy(1)=int(0,kind=mpi_address_kind)
                          !!displacementsy(2)=int(ii-1,kind=mpi_address_kind)
                          blocklengthsy(1)=ii-is2+1
                          blocklengthsy(2)=ie2-comgp%ise(3)+1
                          displacementsy(1)=int(0*lzd%glr%d%n1i*size_of_double,kind=mpi_address_kind)
                          displacementsy(2)=int((comgp%ise(3)-1)*lzd%glr%d%n1i*size_of_double,kind=mpi_address_kind)
                      else
                          nsegy=1
                          blocklengthsy(1)=comgp%ise(4)-comgp%ise(3)+1
                          blocklengthsy(2)=0
                          displacementsy(1)=int(0*lzd%glr%d%n1i*size_of_double,kind=mpi_address_kind)
                          displacementsy(2)=int(comgp%ise(4)*lzd%glr%d%n1i*size_of_double,kind=mpi_address_kind)
                      end if

                      !!!if(is3<is3min .or. ioverlap==1) then
                      !!!    is3min=is3
                      !!!end if
                      !!!if(ie3>ie3max .or. ioverlap==1) then
                      !!!    ie3max=ie3
                      !!!end if
                      !!if (.not.bgq) then
                          ! Communicate only the essential part, i.e. set the starting point to this part
                          istsource = ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i + ioffsety*lzd%glr%d%n1i + ioffsetx
                      !!else
                      !!    ! Communicate the entire slice, i.e. set the starting point to the beginning of a slice
                      !!    istsource = ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i
                      !!end if
                      ncount = 1
                      comgp%comarr(1,ioverlap)=kproc
                      comgp%comarr(2,ioverlap)=istsource
                      comgp%comarr(3,ioverlap)=iproc
                      comgp%comarr(4,ioverlap)=istdest
                      comgp%comarr(5,ioverlap)=iie3(j3)-iis3(j3)+1 !nlen3(j3) !iie3(j3)-iis3(j3)+1
                      comgp%comarr(6,ioverlap)=lzd%glr%d%n1i*lzd%glr%d%n2i
                      if (.not. datatype_defined) then
                          !!write(*,'(a,8i8)') 'iproc, nsegx, blocklengthsx, displacementsx, comgp%ise(1), comgp%ise(2)', &
                          !!                    iproc, nsegx, blocklengthsx, displacementsx, comgp%ise(1), comgp%ise(2)
                          types(:)=mpi_double_precision
                          if (nproc>1) then
                              call mpi_type_create_struct(nsegx, blocklengthsx, displacementsx, &
                                   types, xline_type, ierr)
                              call mpi_type_commit(xline_type, ierr)
                              call mpi_type_size(xline_type, ii, ierr)
                              call mpi_type_get_extent(xline_type, lb, extent, ierr)
                          end if
                          !!write(*,'(a,4i10)') 'iproc, size, lb, extent, of xline_type', iproc, ii, lb, extent
                          !write(*,*) 'iproc, size of xline_type', iproc, ii
                          !!call mpi_type_vector(comgp%ise(4)-comgp%ise(3)+1, comgp%ise(2)-comgp%ise(1)+1, &
                          !!     lzd%glr%d%n1i, mpi_double_precision, comgp%mpi_datatypes(0), ierr)
                          !derived_types = f_malloc(comgp%ise(4)-comgp%ise(3)+1,id='derived_types')
                          !derived_types(:)=xline_type
                          !!call mpi_type_create_hvector(comgp%ise(4)-comgp%ise(3)+1, 1, &
                          !!     int(size_of_double*lzd%glr%d%n1i,kind=mpi_address_kind), &
                          !!     xline_type, comgp%mpi_datatypes(0), ierr)
                          ! Now create a type describing one block
                          xyblock_type(:)=0 !just to initialize
                          iel = 0
                          do iseg=1,nsegy
                              !!write(*,*) 'iproc, iseg, blocklengthsy(iseg)', iproc, iseg, blocklengthsy(iseg)
                              if (nproc>1) then
                                  call mpi_type_create_hvector(blocklengthsy(iseg), 1, &
                                       int(size_of_double*lzd%glr%d%n1i,kind=mpi_address_kind), &
                                       xline_type, xyblock_type(iseg), ierr)
                                  call mpi_type_commit(xyblock_type(iseg), ierr)
                                  call mpi_type_size(xyblock_type(iseg), ii, ierr)
                                  call mpi_type_get_extent(xyblock_type(iseg), lb, extent, ierr)
                              end if
                              !!write(*,'(a,4i14)') 'iproc, size, lb, extent, of xyblock_type(iseg)', iproc, ii, lb, extent
                              types(iseg)=xyblock_type(iseg)
                              nblocksy(iseg)=1
                              do ileny=1,blocklengthsy(iseg)
                                  do isegx=1,nsegx
                                      iel = iel + 1
                                      ioffset = int(displacementsy(iseg)/size_of_double,kind=4) + &
                                                (ileny-1)*lzd%glr%d%n1i + &
                                                int(displacementsx(isegx)/size_of_double,kind=4)
                                      comgp%onedtypearr(1,iel) = ioffset
                                      comgp%onedtypearr(2,iel) = blocklengthsx(isegx)
                                      !write(*,*) 'isegx, blocklengthsx(isegx)', isegx, blocklengthsx(isegx)
                                      !write(*,*) 'iproc, ioverlap, comgp%noverlaps', iproc, ioverlap, comgp%noverlaps
                                      !write(*,'(a,3i5,2i12)') 'iproc, iel, ioverlap, 1darr', &
                                      !    iproc, iel, ioverlap, comgp%onedtypearr(:,iel,ioverlap)
                                  end do
                              end do
                          end do
                          comgp%onedtypeovrlp = iel
                          types(:)=xyblock_type
                          if (nproc>1) then
                              call mpi_type_create_struct(nsegy, nblocksy, displacementsy, &
                                   types, comgp%mpi_datatypes(0), ierr)
                              !call f_free(derived_types)
                              call mpi_type_commit(comgp%mpi_datatypes(0), ierr)
                              call mpi_type_size(comgp%mpi_datatypes(0), ii, ierr)
                              call mpi_type_get_extent(comgp%mpi_datatypes(0), lb, extent, ierr)
                              !!write(*,'(a,4i14)') 'iproc, size, lb, extent, of comgp%mpi_datatypes(0)', iproc, ii, lb, extent
                              do iseg=1,nsegy
                                  call mpi_type_free(xyblock_type(iseg), ierr)
                              end do
                              call mpi_type_free(xline_type, ierr)
                          end if
                          datatype_defined=.true.
                  end if
                  !!istdest = istdest + &
                  !!          (ie3-is3+1)*(comgp%ise(2)-comgp%ise(1)+1)*(comgp%ise(4)-comgp%ise(3)+1)
                  !!comgp%nrecvBuf = comgp%nrecvBuf + &
                  !!          (ie3-is3+1)*(comgp%ise(2)-comgp%ise(1)+1)*(comgp%ise(4)-comgp%ise(3)+1)
                  if (nproc>1) then
                      call mpi_type_size(comgp%mpi_datatypes(0), size_datatype, ierr)
                      size_datatype=size_datatype/size_of_double
                  else
                      size_datatype = 0
                      do i=1,comgp%onedtypeovrlp
                          size_datatype = size_datatype + comgp%onedtypearr(2,i)
                      end do
                  end if
                  istdest = istdest + nlen3(j3)*size_datatype
                  comgp%nrecvBuf = comgp%nrecvBuf + nlen3(j3)*size_datatype
                  !!write(*,'(a,4i9)') 'j3, nlen3(j3), size_datatype, comgp%nrecvBuf', j3, nlen3(j3), size_datatype, comgp%nrecvBuf
              !!else if(ie3j > lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F')then
              !!     stop 'WILL PROBABLY NOT WORK!'
              !!     ie3j = comgp%ise(6) - lzd%Glr%d%n3i
              !!     if(ie3j>=is3k) then
              !!         is3=max(0,is3k) ! starting index in z dimension for data to be sent
              !!         ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
              !!         ioffsetz=is3-0 ! starting index (in z direction) of data to be sent (actually it is the index -1)
              !!         ioverlap=ioverlap+1
              !!         !tag=tag+1
              !!         !!tag=p2p_tag(jproc)
              !!         if(is3<is3min .or. ioverlap==1) then
              !!             is3min=is3
              !!         end if
              !!         if(ie3>ie3max .or. ioverlap==1) then
              !!             ie3max=ie3
              !!         end if
              !!         !!call setCommunicationPotential(kproc, is3, ie3, ioffsetz, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
              !!         !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
              !!         istsource=ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i+1
              !!         !ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
              !!         ncount=lzd%glr%d%n1i*lzd%glr%d%n2i
              !!         call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap))
              !!         comgp%comarr(7,ioverlap)=(ie3-is3+1)
              !!         comgp%comarr(8,ioverlap)=lzd%glr%d%n1i*lzd%glr%d%n2i
              !!         istdest = istdest + (ie3-is3+1)*ncount
              !!         if(iproc==jproc) then
              !!             comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
              !!         end if
              !!     end if
              !!     if(is3j <= ie3k)then
              !!         is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
              !!         ie3=min(lzd%Glr%d%n3i,ie3k) ! ending index in z dimension for data to be sent
              !!         ioffsetz=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
              !!         ioverlap=ioverlap+1
              !!         !tag=tag+1
              !!         !!tag=p2p_tag(jproc)
              !!         if(is3<is3min .or. ioverlap==1) then
              !!             is3min=is3
              !!         end if
              !!         if(ie3>ie3max .or. ioverlap==1) then
              !!             ie3max=ie3
              !!         end if
              !!         !!call setCommunicationPotential(kproc, is3, ie3, ioffsetz, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
              !!         !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
              !!         istsource=ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i+1
              !!         !ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
              !!         ncount=lzd%glr%d%n1i*lzd%glr%d%n2i
              !!         call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap))
              !!         comgp%comarr(7,ioverlap)=ie3-is3+1
              !!         comgp%comarr(8,ioverlap)=lzd%glr%d%n1i*lzd%glr%d%n2i
              !!         istdest = istdest + (ie3-is3+1)*ncount
              !!         if(iproc==jproc) then
              !!             comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
              !!         end if
              !!     end if
              end do
              end if
          end do
          !!comgp%ise3(1,jproc)=is3min
          !!comgp%ise3(2,jproc)=ie3max
          !!if (iproc==0) write(*,*) 'is3min,comgp%ise(5,jproc)', is3min,comgp%ise(5,jproc)
          !!if (iproc==0) write(*,*) 'ie3max,comgp%ise(6,jproc)', ie3max,comgp%ise(6,jproc)
          !if (comgp%ise(5,jproc)/=is3min) stop 'ERROR 1'
          !if (comgp%ise(6,jproc)/=ie3max) stop 'ERROR 2'
          if(ioverlap/=comgp%noverlaps) stop 'ioverlap/=comgp%noverlaps'
    
      !else nproc_if ! monoproc

      !    !write(*,*) 'comgp%ise',comgp%ise
    
      !    !comgp%nrecvbuf = (comgp%ise(2)-comgp%ise(1)+1)*(comgp%ise(4)-comgp%ise(3)+1)*&
      !    !                 (comgp%ise(6)-comgp%ise(5)+1)
      !    ! Probably too much, but ok for the moment
      !    comgp%nrecvbuf = lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
      !
      !end if nproc_if
    
      ! This is the size of the communication buffer without spin
      comgp%nrecvbuf=max(comgp%nrecvbuf,1)

      ! Copy the spin value
      comgp%nspin=nspin
      
      ! To indicate that no communication is going on.
      comgp%communication_complete=.true.
      !!comgp%messages_posted=.false.
    
      call timing(iproc,'init_commPot  ','OF')
    
    end subroutine initialize_communication_potential



    !> Routines for the cubic version

    !> Partition the orbitals between processors to ensure load balancing
    !! the criterion will depend on GPU computation
    !! and/or on the sizes of the different localisation region.
    !!
    !! Calculate the number of elements to be sent to each process
    !! and the array of displacements.
    !! Cubic strategy: 
    !!    - the components are equally distributed among the wavefunctions
    !!    - each processor has all the orbitals in transposed form
    !!    - each wavefunction is equally distributed in its transposed form
    !!    - this holds for each k-point, which regroups different processors
    subroutine orbitals_communicators(iproc,nproc,lr,orbs,comms,basedist)
      use module_base
      use module_types
      use yaml_output, only: yaml_toa
      implicit none
      integer, intent(in) :: iproc,nproc
      type(locreg_descriptors), intent(in) :: lr
      type(orbitals_data), intent(inout) :: orbs
      type(comms_cubic), intent(out) :: comms
      integer, dimension(0:nproc-1,orbs%nkpts), intent(in), optional :: basedist
      !local variables
      character(len=*), parameter :: subname='orbitals_communicators'
      logical :: yesorb,yescomp
      integer :: jproc,nvctr_tot,ikpts,iorbp,jorb,norb_tot,ikpt
      integer :: nkptsp,ierr,kproc,jkpts,jkpte,jsorb,lubo,lubc,info,jkpt
      integer, dimension(:), allocatable :: mykpts
      logical, dimension(:), allocatable :: GPU_for_comp
      integer, dimension(:,:), allocatable :: nvctr_par,norb_par !<for all the components and orbitals (with k-pts)
      
      !check of allocation of important arrays
      if (.not. associated(orbs%norb_par)) then
         write(*,*)'ERROR: norb_par array not allocated'
         stop
      end if
   
      !Allocations of nvctr_par and norb_par
      nvctr_par = f_malloc((/ 0.to.nproc-1, 0.to.orbs%nkpts /),id='nvctr_par')
      norb_par = f_malloc((/ 0.to.nproc-1, 0.to.orbs%nkpts /),id='norb_par')
      mykpts = f_malloc(orbs%nkpts,id='mykpts')
    
      !initialise the arrays
      do ikpts=0,orbs%nkpts
         do jproc=0,nproc-1
            nvctr_par(jproc,ikpts)=0 
            norb_par(jproc,ikpts)=0 
         end do
      end do
    
      !calculate the same k-point distribution for the orbitals
      !assign the k-point to the given orbital, counting one orbital after each other
      jorb=1
      ikpts=1
      do jproc=0,nproc-1
         do iorbp=1,orbs%norb_par(jproc,0)
            norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
            if (mod(jorb,orbs%norb)==0) then
               ikpts=ikpts+1
            end if
            jorb=jorb+1
         end do
      end do
      !some checks
      if (orbs%norb /= 0) then
         !check the distribution
         do ikpts=1,orbs%nkpts
            !print *,'partition',ikpts,orbs%nkpts,'ikpts',norb_par(:,ikpts)
            norb_tot=0
            do jproc=0,nproc-1
               norb_tot=norb_tot+norb_par(jproc,ikpts)
            end do
            if(norb_tot /= orbs%norb) then
               call f_err_throw('Orbital partition is incorrect for k-point'//&
                    trim(yaml_toa(ikpts))//'; expected '//&
                    trim(yaml_toa(orbs%norb))//' orbitals, found'//&
                    trim(yaml_toa(norb_tot)),&
                    err_name='BIGDFT_RUNTIME_ERROR')
            end if
         end do
      end if
    
      !balance the components between processors
      !in the most symmetric way
      !here the components are taken into account for all the k-points
    
      !create an array which indicate which processor has a GPU associated 
      !from the viewpoint of the BLAS routines (deprecated, not used anymore)
      GPU_for_comp = f_malloc(0.to.nproc-1,id='GPU_for_comp')
    
      if (nproc > 1 .and. .not. GPUshare) then
         call MPI_ALLGATHER(GPUblas,1,MPI_LOGICAL,GPU_for_comp(0),1,MPI_LOGICAL,&
              bigdft_mpi%mpi_comm,ierr)
      else
         GPU_for_comp(0)=GPUblas
      end if
    
      call f_free(GPU_for_comp)
    
      !old k-point repartition
    !!$  !decide the repartition for the components in the same way as the orbitals
    !!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par)
    
    !!$  ikpts=1
    !!$  ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
    !!$  do jproc=0,nproc-1
    !!$     loop_comps: do
    !!$        if (nvctr_par(jproc,0) >= ncomp_res) then
    !!$           nvctr_par(jproc,ikpts)= ncomp_res
    !!$           ikpts=ikpts+1
    !!$           nvctr_par(jproc,0)=nvctr_par(jproc,0)-ncomp_res
    !!$           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
    !!$        else
    !!$           nvctr_par(jproc,ikpts)= nvctr_par(jproc,0)
    !!$           ncomp_res=ncomp_res-nvctr_par(jproc,0)
    !!$           nvctr_par(jproc,0)=0
    !!$           exit loop_comps
    !!$        end if
    !!$        if (nvctr_par(jproc,0) == 0 ) then
    !!$           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
    !!$           exit loop_comps
    !!$        end if
    !!$
    !!$     end do loop_comps
    !!$  end do
    
      !new k-point repartition
      if (present(basedist)) then
         do jkpt=1,orbs%nkpts
            do jproc=0,nproc-1
               nvctr_par(jproc,jkpt)=basedist(jproc,jkpt)
            end do
         end do
      else
         !first try the naive repartition
         call kpts_to_procs_via_obj(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par(0,1))
      end if
      !then silently check whether the distribution agree
      info=-1
      call check_kpt_distributions(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
           norb_par(0,1),nvctr_par(0,1),info,lubo,lubc)
      if (info/=0 .and. .not. present(basedist)) then !redo the distribution based on the orbitals scheme
         info=-1
         call components_kpt_distribution(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),norb_par(0,1),nvctr_par(0,1))
         call check_kpt_distributions(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
              norb_par(0,1),nvctr_par(0,1),info,lubo,lubc)
      end if
      if (info /=0) then
         !if (iproc==0) then
         !   write(*,*)'ERROR for nproc,nkpts,norb,nvctr',nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
         !   call print_distribution_schemes(nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
         !end if
         !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
         !stop
         if (iproc==0) call print_distribution_schemes(nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
         call f_err_throw('ERROR for nproc,nkpts,norb,nvctr' // &
              & trim(yaml_toa( (/ nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) /) )), &
              & err_id=BIGDFT_RUNTIME_ERROR)
      end if
    
    !write(*,'(a,i2,3x,8i7,i10)') 'iproc, nvctr_par(jproc), sum', iproc, (nvctr_par(jproc,1), jproc=0,nproc-1), sum(nvctr_par(:,1))
    !write(*,*) 'iproc, (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norbp', iproc, (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norbp
      !some checks
      !check the distribution
      do ikpts=1,orbs%nkpts
         !print *,'iproc,cpts:',lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nvctr_par(:,ikpts)
         nvctr_tot=0
         do jproc=0,nproc-1
            nvctr_tot=nvctr_tot+nvctr_par(jproc,ikpts)
         end do
         if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
            write(*,*)'ERROR: partition of components incorrect, kpoint:',ikpts
            stop
         end if
      end do
    
      !this function which associates a given k-point to a processor in the component distribution
      !the association is chosen such that each k-point is associated to only
      !one processor
      !if two processors treat the same k-point the processor which highest rank is chosen
      do ikpts=1,orbs%nkpts
         loop_jproc: do jproc=nproc-1,0,-1
            if (nvctr_par(jproc,ikpts) /= 0) then
               orbs%ikptproc(ikpts)=jproc
               exit loop_jproc
            end if
         end do loop_jproc
      end do
      
      !print*,'check',orbs%ikptproc(:)
    
    !write(*,*) 'orbs%norb_par',orbs%norb_par
    
      !calculate the number of k-points treated by each processor in both
      ! the component distribution and the orbital distribution.
      !to have a correct distribution, a k-point should be divided between the same processors
      nkptsp=0
      orbs%iskpts=-1
      do ikpts=1,orbs%nkpts
         if (nvctr_par(iproc,ikpts) /= 0 .or. norb_par(iproc,ikpts) /= 0) then
            if (orbs%iskpts == -1) orbs%iskpts=ikpts-1
            nkptsp=nkptsp+1
            mykpts(nkptsp) = ikpts
         end if
      end do
      orbs%nkptsp=nkptsp
    
    !!$  allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
    !!$  call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
    !!$  orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)
    
      !print the distribution scheme used for this set of orbital
      !in the case of multiple k-points
      if (iproc == 0 .and. verbose > 1 .and. orbs%nkpts > 1) then
         call print_distribution_schemes(nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
      end if
    
      !print *,iproc,orbs%nkptsp,orbs%norbp,orbs%norb,orbs%nkpts
      !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
      !call MPI_FINALIZE(ierr)
      !stop
      !check that for any processor the orbital k-point repartition is contained into the components
      if (orbs%norb /= 0) then
         do jproc=0,nproc-1
            jsorb=0
            do kproc=0,jproc-1
               jsorb=jsorb+orbs%norb_par(kproc,0)
            end do
            jkpts=min(jsorb/orbs%norb+1,orbs%nkpts)
            if (nvctr_par(jproc,jkpts) == 0 .and. orbs%norb_par(jproc,0) /=0 ) then
               if (iproc ==0) write(*,*)'ERROR, jproc: ',jproc,' the orbital k-points distribution starts before the components one'
               !print *,jsorb,jkpts,jproc,orbs%iskpts,nvctr_par(jproc,jkpts)
               stop
            end if
            jkpte=min((jsorb+orbs%norb_par(jproc,0)-1)/orbs%norb+1,orbs%nkpts)
            if (nvctr_par(jproc,jkpte) == 0 .and. orbs%norb_par(jproc,0) /=0) then
               if (iproc ==0) write(*,*)'ERROR, jproc: ',jproc,&
                    ' the orbital k-points distribution ends after the components one'
               print *,jsorb,jkpte,jproc,orbs%iskpts,orbs%nkptsp,nvctr_par(jproc,jkpte)
               stop
            end if
         end do
      end if
    
      !before printing the distribution schemes, check that the two distributions contain
      !the same k-points
      yesorb=.false.
      kpt_components: do ikpts=1,orbs%nkptsp
         ikpt=orbs%iskpts+ikpts
         do jorb=1,orbs%norbp
            if (orbs%iokpt(jorb) == ikpt) yesorb=.true.
         end do
         if (.not. yesorb .and. orbs%norbp /= 0) then
            write(*,*)' ERROR: processor ', iproc,' kpt ',ikpt,&
                 ' not found in the orbital distribution'
            call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
         end if
      end do kpt_components
    
      yescomp=.false.
      kpt_orbitals: do jorb=1,orbs%norbp
         ikpt=orbs%iokpt(jorb)   
         do ikpts=1,orbs%nkptsp
            if (orbs%iskpts+ikpts == ikpt) yescomp=.true.
         end do
         if (.not. yescomp) then
            write(*,*)' ERROR: processor ', iproc,' kpt,',ikpt,&
                 'not found in the component distribution'
            call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
         end if
      end do kpt_orbitals
    
      !print *,'AAAAiproc',iproc,orbs%iskpts,orbs%iskpts+orbs%nkptsp
    
      !allocate communication arrays
      comms%nvctr_par = f_malloc_ptr((/ 0.to.nproc-1, 0.to.orbs%nkpts /),id='comms%nvctr_par')
      comms%ncntd = f_malloc_ptr(0.to.nproc-1,id='comms%ncntd')
      comms%ncntt = f_malloc_ptr(0.to.nproc-1,id='comms%ncntt')
      comms%ndspld = f_malloc_ptr(0.to.nproc-1,id='comms%ndspld')
      comms%ndsplt = f_malloc_ptr(0.to.nproc-1,id='comms%ndsplt')
    
      !assign the partition of the k-points to the communication array
      !calculate the number of componenets associated to the k-point
      do jproc=0,nproc-1
         comms%nvctr_par(jproc,0)=0
         do ikpt=1,orbs%nkpts
            comms%nvctr_par(jproc,0)=comms%nvctr_par(jproc,0)+&
                 nvctr_par(jproc,ikpt) 
            comms%nvctr_par(jproc,ikpt)=nvctr_par(jproc,ikpt)
         end do
      end do
    !!$  do ikpts=1,orbs%nkptsp
    !!$     ikpt=orbs%iskpts+ikpts!orbs%ikptsp(ikpts)
    !!$     do jproc=0,nproc-1
    !!$        comms%nvctr_par(jproc,ikpts)=nvctr_par(jproc,ikpt) 
    !!$     end do
    !!$  end do
    
      !with this distribution the orbitals and the components are ordered following k-points
      !there must be no overlap for the components
      !here we will print out the k-points components distribution, in the transposed and in the direct way
    
      do jproc=0,nproc-1
         comms%ncntd(jproc)=0
         do ikpts=1,orbs%nkpts
            comms%ncntd(jproc)=comms%ncntd(jproc)+&
                 nvctr_par(jproc,ikpts)*norb_par(iproc,ikpts)*orbs%nspinor
         end do
      end do
      comms%ndspld(0)=0
      do jproc=1,nproc-1
         comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
      end do
      !receive buffer
      do jproc=0,nproc-1
         comms%ncntt(jproc)=0
         do ikpts=1,orbs%nkpts
            comms%ncntt(jproc)=comms%ncntt(jproc)+&
                 nvctr_par(iproc,ikpts)*norb_par(jproc,ikpts)*orbs%nspinor
         end do
      end do
      comms%ndsplt(0)=0
      do jproc=1,nproc-1
         comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)
      end do
    
      !print *,'iproc,comms',iproc,comms%ncntd,comms%ndspld,comms%ncntt,comms%ndsplt
    
      call f_free(nvctr_par)
      call f_free(norb_par)
      call f_free(mykpts)
    
      !calculate the dimension of the wavefunction
      !for the given processor (this is only the cubic strategy)
      orbs%npsidim_orbs=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor
      orbs%npsidim_comp=sum(comms%ncntt(0:nproc-1))
        
    !!$  orbs%npsidim=max((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor,&
    !!$       sum(comms%ncntt(0:nproc-1)))
    
    END SUBROUTINE orbitals_communicators



end module communications_init
