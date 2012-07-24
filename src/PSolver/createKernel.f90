!> @file
!!    Routines to create the kernel for Poisson solver
!! @author
!!    Copyright (C) 2002-2011 BigDFT group  (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

function pkernel_init(iproc,nproc,taskgroup_size,igpu,geocode,ndims,hgrids,itype_scf,&
     mu0_screening,angrad) result(kernel)
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: itype_scf,iproc,nproc,taskgroup_size,igpu
  character(len=1), intent(in) :: geocode
  integer, dimension(3), intent(in) :: ndims
  real(gp), dimension(3), intent(in) :: hgrids
  real(kind=8), intent(in), optional :: mu0_screening
  real(gp), dimension(3), intent(in), optional :: angrad
  type(coulomb_operator) :: kernel
  !local variables
  logical :: dump
  real(dp) :: alphat,betat,gammat,mu0t
  integer :: base_grp,group_id,temp_comm,grp,i,j,ierr,nthreads,group_size
  integer, dimension(nproc) :: group_list !using nproc instead of taskgroup_size
  !$ integer :: omp_get_max_threads

  group_size=taskgroup_size
  !nullification
  kernel=pkernel_null()

  if (present(angrad)) then
     kernel%angrad=angrad
  else
     alphat = 2.0_dp*datan(1.0_dp)
     betat = 2.0_dp*datan(1.0_dp)
     gammat = 2.0_dp*datan(1.0_dp)
     kernel%angrad=(/alphat,betat,gammat/)
  end if
  if (.not. present(mu0_screening)) then
     mu0t=0.0_gp
  else
     mu0t=mu0_screening
  end if
  kernel%mu=mu0t

  !geocode and ISF family
  kernel%geocode=geocode
  kernel%itype_scf=itype_scf

  !dimensions and grid spacings
  kernel%ndims=ndims
  kernel%hgrids=hgrids

  !part of decision of the communication procedure
  kernel%iproc_world=iproc
  kernel%iproc=iproc
  kernel%nproc=nproc
  kernel%mpi_comm=MPI_COMM_WORLD
  kernel%igpu=igpu  

  dump=iproc==0

  if (dump) then 
     if (mu0t==0.0_gp) then 
        call yaml_open_map('Poisson Kernel Initialization')
     else
        call yaml_open_map('Helmholtz Kernel Initialization')
         call yaml_map('Screening Length (AU)',1/mu0t,fmt='(g25.17)')
     end if
  end if

  if (nproc >1 .and. group_size > 0) then
     !create taskgroups if the number of processes is bigger than one and multiple of group_size
     !print *,'am i here',nproc >1 .and. group_size < nproc .and. mod(nproc,group_size)==0
!     print *,nproc,group_size,mod(nproc,group_size)
     
     if (nproc >1 .and. group_size < nproc .and. mod(nproc,group_size)==0) then
        group_id=iproc/group_size
        kernel%iproc=mod(iproc,group_size)
        kernel%nproc=group_size
        !take the base group
        call MPI_COMM_GROUP(MPI_COMM_WORLD,base_grp,ierr)
        if (ierr /=0) then
           call yaml_warning('Problem in group creation, ierr:'//yaml_toa(ierr))
           call MPI_ABORT(MPI_COMM_WORLD,ierr)
        end if
        do i=0,nproc/group_size-1
           !define the new groups and thread_id
           do j=0,group_size-1
              group_list(j+1)=i*group_size+j
           enddo
           call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
           if (ierr /=0) then
              call yaml_warning('Problem in group inclusion, ierr:'//yaml_toa(ierr))
              call MPI_ABORT(MPI_COMM_WORLD,ierr)
           end if
           call MPI_COMM_CREATE(MPI_COMM_WORLD,grp,temp_comm,ierr)
           if (ierr /=0) then
              call yaml_warning('Problem in communicator creator, ierr:'//yaml_toa(ierr))
              call MPI_ABORT(MPI_COMM_WORLD,ierr)
           end if
           !print *,'i,group_id,temp_comm',i,group_id,temp_comm
           if (i.eq.group_id) kernel%mpi_comm=temp_comm
        enddo
        if (dump) then
             call yaml_map('Total No. of Taskgroups created',nproc/kernel%nproc)
        end if
        
     end if
  end if

  !gpu can be used only for one nproc
  if (kernel%nproc > 1) kernel%igpu=0

  !-------------------
  nthreads=0
  if (kernel%iproc_world ==0) then
     !$ nthreads = omp_get_max_threads()
     call yaml_map('MPI tasks',kernel%nproc)
     if (nthreads /=0) call yaml_map('OpenMP threads per task',nthreads)
     if (kernel%igpu==1) call yaml_map('Kernel copied on GPU',.true.)
     call yaml_close_map() !kernel
  end if

end function pkernel_init

subroutine pkernel_free(kernel,subname)
  use module_base
  implicit none
  character(len=*), intent(in) :: subname
  type(coulomb_operator), intent(inout) :: kernel
  !local variables
  integer :: i_all,i_stat

  if (associated(kernel%kernel)) then
     i_all=-product(shape(kernel%kernel))*kind(kernel%kernel)
     deallocate(kernel%kernel,stat=i_stat)
     call memocc(i_stat,i_all,'kernel',subname)
  end if

  !free GPU data
  if (kernel%igpu == 1) then
     call cudafree(kernel%work1_GPU)
     call cudafree(kernel%work2_GPU)
     call cudafree(kernel%k_GPU)
  end if
  

  !cannot yet free the communicators of the poisson kernel

end subroutine pkernel_free


!> Allocate a pointer which corresponds to the zero-padded FFT slice needed for
!! calculating the convolution with the kernel expressed in the interpolating scaling
!! function basis. The kernel pointer is unallocated on input, allocated on output.
!! SYNOPSIS
!!    @param geocode  Indicates the boundary conditions (BC) of the problem:
!!              - 'F' free BC, isolated systems.
!!                    The program calculates the solution as if the given density is
!!                    "alone" in R^3 space.
!!              - 'S' surface BC, isolated in y direction, periodic in xz plane                
!!                    The given density is supposed to be periodic in the xz plane,
!!                    so the dimensions in these direction mus be compatible with the FFT
!!                    Beware of the fact that the isolated direction is y!
!!              - 'P' periodic BC.
!!                    The density is supposed to be periodic in all the three directions,
!!                    then all the dimensions must be compatible with the FFT.
!!                    No need for setting up the kernel.
!!              - 'W' Wires BC.
!!                    The density is supposed to be periodic in z direction, 
!!                    which has to be compatible with the FFT.
!!    @param iproc,nproc number of process, number of processes
!!    @param n01,n02,n03 dimensions of the real space grid to be hit with the Poisson Solver
!!    @param itype_scf   order of the interpolating scaling functions used in the decomposition
!!    @param hx,hy,hz grid spacings. For the isolated BC case for the moment they are supposed to 
!!                    be equal in the three directions
!!    @param kernel   pointer for the kernel FFT. Unallocated on input, allocated on output.
!!                    Its dimensions are equivalent to the region of the FFT space for which the
!!                    kernel is injective. This will divide by two each direction, 
!!                    since the kernel for the zero-padded convolution is real and symmetric.
!!    @param gpu      tag for CUDA gpu   0: CUDA GPU is disabled
!!
!! @warning
!!    Due to the fact that the kernel dimensions are unknown before the calling, the kernel
!!    must be declared as pointer in input of this routine.
!!    To avoid that, one can properly define the kernel dimensions by adding 
!!    the nd1,nd2,nd3 arguments to the PS_dim4allocation routine, then eliminating the pointer
!!    declaration.
subroutine pkernel_set(kernel,wrtmsg) !optional arguments
  use module_base, only: ndebug
  use yaml_output
  implicit none
  logical, intent(in) :: wrtmsg
  type(coulomb_operator), intent(inout) :: kernel
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !local variables
  logical :: dump
  character(len=*), parameter :: subname='createKernel'
  integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,i_stat
  integer :: jproc,nlimd,nlimk,jfd,jhd,jzd,jfk,jhk,jzk,npd,npk
  real(kind=8) :: alphat,betat,gammat,mu0t
  real(kind=8), dimension(:), pointer :: pkernel2
  integer :: i1,i2,i3,j1,j2,j3,ind,indt,switch_alg,size2,sizek,i_all
  integer,dimension(3) :: n
  !$ integer :: omp_get_max_threads

  call timing(kernel%iproc_world,'PSolvKernel   ','ON')

  dump=wrtmsg .and. kernel%iproc_world==0

  mu0t=kernel%mu
  alphat=kernel%angrad(1)
  betat=kernel%angrad(2)
  gammat=kernel%angrad(3)

  if (dump) then 
     if (mu0t==0.0_gp) then 
        call yaml_open_map('Poisson Kernel Creation')
     else
        call yaml_open_map('Helmholtz Kernel Creation')
        call yaml_map('Screening Length (AU)',1/mu0t,fmt='(g25.17)')
     end if
  end if



  if (kernel%geocode == 'P') then
     
     if (dump) then
        call yaml_map('Boundary Conditions','Periodic')
     end if
     call P_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,kernel%nproc)

     allocate(kernel%kernel(nd1*nd2*nd3/kernel%nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel%kernel,'kernel',subname)

     call Periodic_Kernel(n1,n2,n3,nd1,nd2,nd3,&
          kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),&
          kernel%itype_scf,kernel%kernel,kernel%iproc,kernel%nproc,mu0t,alphat,betat,gammat)

     nlimd=n2
     nlimk=n3/2+1

  else if (kernel%geocode == 'S') then
     
     if (dump) then
        call yaml_map('Boundary Conditions','Surface')
     end if
     !Build the Kernel
     call S_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,kernel%nproc,kernel%igpu)

     allocate(kernel%kernel(nd1*nd2*nd3/kernel%nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel%kernel,'kernel',subname)

     !the kernel must be built and scattered to all the processes
     call Surfaces_Kernel(kernel%iproc,kernel%nproc,kernel%mpi_comm,n1,n2,n3,m3,nd1,nd2,nd3,&
          kernel%hgrids(1),kernel%hgrids(3),kernel%hgrids(2),&
          kernel%itype_scf,kernel%kernel,mu0t,alphat,betat,gammat)

     !last plane calculated for the density and the kernel
     nlimd=n2
     nlimk=n3/2+1

  else if (kernel%geocode == 'F') then

     if (dump) then
        call yaml_map('Boundary Conditions','Free')
     end if
!     print *,'debug',kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3)
     !Build the Kernel
     call F_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
          md1,md2,md3,nd1,nd2,nd3,kernel%nproc,kernel%igpu)

     allocate(kernel%kernel(nd1*nd2*nd3/kernel%nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel%kernel,'kernel',subname)

     !the kernel must be built and scattered to all the processes
     call Free_Kernel(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          n1,n2,n3,nd1,nd2,nd3,kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),&
          kernel%itype_scf,kernel%iproc,kernel%nproc,kernel%kernel,mu0t)

     !last plane calculated for the density and the kernel
     nlimd=n2/2
     nlimk=n3/2+1
     
  else if (kernel%geocode == 'W') then

     if (dump) then
        call yaml_map('Boundary Conditions','Wire')
     end if
     call W_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,kernel%nproc,kernel%igpu)

     allocate(kernel%kernel(nd1*nd2*nd3/kernel%nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel%kernel,'kernel',subname)

     call Wires_Kernel(kernel%iproc,kernel%nproc,&
          kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),&
          n1,n2,n3,nd1,nd2,nd3,kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),&
          kernel%itype_scf,kernel%kernel,mu0t)

     nlimd=n2
     nlimk=n3/2+1
              
  else
     
     !if (iproc==0) 
     write(*,'(1x,a,3a)')'createKernel, geocode not admitted',kernel%geocode

     stop
  end if
!print *,'thereAAA',iproc,nproc,kernel%iproc,kernel%nproc,kernel%mpi_comm
!call MPI_BARRIER(kernel%mpi_comm,ierr)

  if (dump) then
     call yaml_open_map('Memory Requirements per MPI task')
       call yaml_map('Density (MB)',8.0_gp*real(md1*md3,gp)*real(md2/kernel%nproc,gp)/(1024.0_gp**2),fmt='(f8.2)')
       call yaml_map('Kernel (MB)',8.0_gp*real(nd1*nd3,gp)*real(nd2/kernel%nproc,gp)/(1024.0_gp**2),fmt='(f8.2)')
       call yaml_map('Full Grid Arrays (MB)',&
            8.0_gp*real(kernel%ndims(1)*kernel%ndims(2),gp)*real(kernel%ndims(3),gp)/(1024.0_gp**2),fmt='(f8.2)')
       !print the load balancing of the different dimensions on screen
     if (kernel%nproc > 1) then
        call yaml_open_map('Load Balancing of calculations')
        !write(*,'(1x,a)')&
        !     'Load Balancing for Poisson Solver related operations:'
        jhd=10000
        jzd=10000
        npd=0
        load_balancing: do jproc=0,kernel%nproc-1
           !print *,'jproc,jfull=',jproc,jproc*md2/nproc,(jproc+1)*md2/nproc
           if ((jproc+1)*md2/kernel%nproc <= nlimd) then
              jfd=jproc
           else if (jproc*md2/kernel%nproc <= nlimd) then
              jhd=jproc
              npd=nint(real(nlimd-(jproc)*md2/kernel%nproc,kind=8)/real(md2/kernel%nproc,kind=8)*100.d0)
           else
              jzd=jproc
              exit load_balancing
           end if
        end do load_balancing
        call yaml_open_map('Density')
         call yaml_map('MPI tasks 0-'//trim(yaml_toa(jfd,fmt='(i5)')),'100%')
         if (jfd < kernel%nproc-1) &
              call yaml_map('MPI task'//trim(yaml_toa(jhd,fmt='(i5)')),trim(yaml_toa(npd,fmt='(i5)'))//'%')
         if (jhd < kernel%nproc-1) &
              call yaml_map('MPI tasks'//trim(yaml_toa(jhd,fmt='(i5)'))//'-'//&
              yaml_toa(kernel%nproc-1,fmt='(i3)'),'0%')
        call yaml_close_map()
        jhk=10000
        jzk=10000
        npk=0
       ! if (geocode /= 'P') then
           load_balancingk: do jproc=0,kernel%nproc-1
              !print *,'jproc,jfull=',jproc,jproc*nd3/kernel%nproc,(jproc+1)*nd3/kernel%nproc
              if ((jproc+1)*nd3/kernel%nproc <= nlimk) then
                 jfk=jproc
              else if (jproc*nd3/kernel%nproc <= nlimk) then
                 jhk=jproc
                 npk=nint(real(nlimk-(jproc)*nd3/kernel%nproc,kind=8)/real(nd3/kernel%nproc,kind=8)*100.d0)
              else
                 jzk=jproc
                 exit load_balancingk
              end if
           end do load_balancingk
           call yaml_open_map('Kernel')
           call yaml_map('MPI tasks 0-'//trim(yaml_toa(jfk,fmt='(i5)')),'100%')
!           print *,'here,npk',npk
           if (jfk < kernel%nproc-1) &
                call yaml_map('MPI task'//trim(yaml_toa(jhk,fmt='(i5)')),trim(yaml_toa(npk,fmt='(i5)'))//'%')
           if (jhk < kernel%nproc-1) &
                call yaml_map('MPI tasks'//trim(yaml_toa(jhk,fmt='(i5)'))//'-'//&
                yaml_toa(kernel%nproc-1,fmt='(i3)'),'0%')
           call yaml_close_map()
        call yaml_map('Complete LB per task','1/3 LB_density + 2/3 LB_kernel')
        call yaml_close_map()
     end if
     call yaml_close_map() !memory

  end if

  if (kernel%igpu == 1) then

    size2=2*n1*n2*n3
    sizek=(n1/2+1)*n2*n3

    call cudamalloc(size2,kernel%work1_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size2,kernel%work2_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizek,kernel%k_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat

    allocate(pkernel2((n1/2+1)*n2*n3+ndebug),stat=i_stat)
    call memocc(i_stat,pkernel2,'pkernel2',subname)

    ! transpose kernel for GPU
    do i3=1,n3
       j3=i3+(i3/(n3/2+2))*(n3+2-2*i3)!injective dimension
       do i2=1,n2
          j2=i2+(i2/(n2/2+2))*(n2+2-2*i2)!injective dimension
          do i1=1,n1
             j1=i1+(i1/(n1/2+2))*(n1+2-2*i1)!injective dimension
             !injective index
             ind=j1+(j2-1)*nd1+(j3-1)*nd1*nd2
             !unfolded index
             indt=i2+(j1-1)*n2+(i3-1)*nd1*n2
             pkernel2(indt)=kernel%kernel(ind)
          end do
       end do
    end do
    !offset to zero
    if (kernel%geocode == 'P') pkernel2(1)=0.0_dp

    if(kernel%geocode == 'P') then
     kernel%geo(1)=1
     kernel%geo(2)=1
     kernel%geo(3)=1
    else if (kernel%geocode == 'S') then
     kernel%geo(1)=1
     kernel%geo(2)=0
     kernel%geo(3)=1
    else if (kernel%geocode == 'F') then
     kernel%geo(1)=0
     kernel%geo(2)=0
     kernel%geo(3)=0
    else if (kernel%geocode == 'W') then
     kernel%geo(1)=0
     kernel%geo(2)=0
     kernel%geo(3)=1
    end if

    call reset_gpu_data((n1/2+1)*n2*n3,pkernel2,kernel%k_GPU)

    n(1)=n1!kernel%ndims(1)*(2-kernel%geo(1))
    n(2)=n3!kernel%ndims(2)*(2-kernel%geo(2))
    n(3)=n2!kernel%ndims(3)*(2-kernel%geo(3))

    call cuda_3d_psolver_general_plan(n,kernel%plan,switch_alg,kernel%geo)

    i_all=-product(shape(pkernel2))*kind(pkernel2)
    deallocate(pkernel2,stat=i_stat)
    call memocc(i_stat,i_all,'pkernel2',subname)
    
    if (dump) call yaml_map('Kernel Copied on GPU',.true.)

 endif
  
!print *,'there',iproc,nproc,kernel%iproc,kernel%nproc,kernel%mpi_comm
!call MPI_BARRIER(kernel%mpi_comm,ierr)
!print *,'okcomm',kernel%mpi_comm,kernel%iproc
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (dump) call yaml_close_map() !kernel

  call timing(kernel%iproc_world,'PSolvKernel   ','OF')

END SUBROUTINE pkernel_set
