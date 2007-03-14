subroutine PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
     rhopot,karray,pot_ion,eh,exc,vxc,offset)
  implicit none
  include 'mpif.h'
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc
  real(kind=8), intent(in) :: hx,hy,hz,offset
  real(kind=8), dimension(*), intent(in) :: karray,pot_ion
  real(kind=8), intent(out) :: eh,exc,vxc
  real(kind=8), dimension(n01,n02,*), intent(inout) :: rhopot
  !local variables
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: i_allocated,i_stat,ierr
  integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh
  integer :: nxc,nwbl,nwbr,nxt,nwb,nxcl,nxcr,nlim
  real(kind=8) :: ehartreeLOC,eexcuLOC,vexcuLOC
  real(kind=8) :: hgrid,scal,newoffset,correction,pot,factor
  real(kind=8), dimension(:,:,:), allocatable :: zf,zfionxc
  integer, dimension(:,:), allocatable :: gather_arr
  real(kind=8), dimension(:), allocatable :: energies_mpi


  !calculate the dimensions wrt the geocode
  if (geocode == 'P') then
     if (iproc==0) &
          write(*,'(1x,a,3(i5),a,i3,a,i3,a)',advance='no')&
          'PSolver, periodic BC, dimensions: ',n01,n02,n03,'  nproc',nproc,'   ixc:',ixc,' ...'
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'S') then
     if (iproc==0) &
          write(*,'(1x,a,3(i5),a,i3,a,i3,a)',advance='no')&
          'PSolver, surfaces BC, dimensions: ',n01,n02,n03,'  nproc',nproc,'   ixc:',ixc,' ...'
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'F') then
     if (iproc==0) &
          write(*,'(1x,a,3(i5),a,i3,a,i3,a)',advance='no')&
          'PSolver, free  BC, dimensions: ',n01,n02,n03,'  nproc',nproc,'   ixc:',ixc,' ...'
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else
     stop 'PSolver: geometry code not admitted'
  end if

  !array allocations
  i_allocated=0
  allocate(zf(md1,md3,md2/nproc),stat=i_stat)
  i_allocated=i_allocated+i_stat
  allocate(zfionxc(md1,md3,md2/nproc),stat=i_stat)
  i_allocated=i_allocated+i_stat
  if (nproc > 1) then
     if (datacode == 'G') then
        allocate(gather_arr(0:nproc-1,2),stat=i_stat)
        i_allocated=i_allocated+i_stat
     end if
     allocate(energies_mpi(6),stat=i_stat)
     i_allocated=i_allocated+i_stat
  end if
  if (i_allocated /= 0) then
     print *,"PSolver: Problem of memory allocation"
     stop
  end if


  call timing(iproc,'Exchangecorr  ','ON')
  !dimension for exchange-correlation (different in the global or distributed case)
  !let us calculate the dimension of the portion of the rhopot array to be passed 
  !to the xc routine
  !this portion will depend on the need of calculating the gradient or not, 
  !and whether the White-Bird correction must be inserted or not 
  !(absent only in the LB ixc=13 case)
  
  !nxc is the effective part of the third dimension that is being processed
  !nxt is the dimension of the part of rhopot that must be passed to the gradient routine
  !nwb is the dimension of the part of rhopot in the wb-postprocessing routine
  !note: nxc <= nwb <= nxt
  !the dimension are related by the values of nwbl and nwbr
  !      nxc+nxcl+nxcr-2 = nwb
  !      nwb+nwbl+nwbr = nxt
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)
  nxc=iend-istart
  if (ixc >= 11 .and. ixc <= 16 .and. geocode == 'F') then
     if (ixc==13) then
        nwbl=min(istart,nordgr)
        nwbr=min(m2-iend,nordgr)
        nxcl=1
        nxcr=1
     else
        if(istart<=nordgr) then
           nxcl=istart+1
           nwbl=0
        else
           nxcl=nordgr+1
           nwbl=min(nordgr,istart-nordgr)
        end if
        if(iend>=m2-nordgr+1) then
           nxcr=m2-iend+1
           nwbr=0
        else
           nxcr=nordgr+1
           nwbr=min(nordgr,m2-nordgr-iend)
        end if
     end if
  else !(for the moment GGA is not implemented in the non free BC)
     nwbl=0
     nwbr=0
     nxcl=1
     nxcr=1
  end if
  nwb=nxcl+nxc+nxcr-2
  nxt=nwbr+nwb+nwbl

  if (datacode=='G') then
     !starting address of rhopot in the case of global i/o
     i3start=istart+2-nxcl-nwbl
  else if (datacode == 'D') then
     !distributed i/o
     i3start=1
  else
     stop 'PSolver: datacode not admitted'
  end if

  !calculate the actual limit of the array for the zero padded FFT
  if (geocode == 'P') then
     nlim=n2
  else if (geocode == 'S') then
     nlim=n2
  else if (geocode == 'F') then
     nlim=n2/2
  end if

!!$  print *,'density must go from',min(istart+1,m2),'to',iend,'with n2/2=',n2/2
!!$  print *,'        it goes from',i3start+nwbl+nxcl-1,'to',i3start+nxc-1

  if (istart+1 <= m2) then 
     call xc_energy(geocode,m1,m2,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,nxcl,nxcr,&
          ixc,hx,hy,hz,rhopot(1,1,i3start),pot_ion,zf,zfionxc,eexcuLOC,vexcuLOC,iproc,nproc)
  else if (istart+1 <= nlim) then !this condition assures that we have perform good zero padding
     do i2=istart+1,min(nlim,istart+md2/nproc)
        j2=i2-istart
        do i3=1,md3
           do i1=1,md1
              zf(i1,i3,j2)=0.d0
              zfionxc(i1,i3,j2)=0.d0 !this is not needed, only if pot is updated in Solver
           end do
        end do
     end do
     eexcuLOC=0.d0
     vexcuLOC=0.d0
  else
     eexcuLOC=0.d0
     vexcuLOC=0.d0
  end if
  call timing(iproc,'Exchangecorr  ','OF')

  !this routine builds the values for each process of the potential (zf), multiplying by scal 
  if(geocode == 'P') then
     !no powers of hgrid because they are incorporated in the plane wave treatment
     scal=1.d0/real(n1*n2*n3,kind=8)
     call P_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,zf,&
          scal,hx,hy,hz)

     !offset correction for the periodic treatment
     if (iproc == 0) newoffset=zf(1,1,1)
     !send the value of the offset to the other processes
     call timing(iproc,'PSolv_commun  ','ON')
     call MPI_BCAST(newoffset,1,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
     call timing(iproc,'PSolv_commun  ','OF')
     correction=offset-newoffset
     factor=0.5d0*hx*hy*hz

  else if (geocode == 'S') then
     !only one power of hgrid 
     scal=hy/real(n1*n2*n3,kind=8)
     call S_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf,&
          scal,hx,hy,hz)!,ehartreeLOC)
     correction=0.d0
     factor=0.5d0*hx*hy*hz
  else if (geocode == 'F') then
     hgrid=max(hx,hy,hz)
     scal=hgrid**3/real(n1*n2*n3,kind=8)
     call F_PoissonSolver(n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nproc,iproc,karray,zf,&
          scal,hgrid)!,ehartreeLOC)
     correction=0.d0
     factor=0.5d0*hgrid**3

  end if
  
  call timing(iproc,'PSolv_comput  ','ON')

  !the value of the shift depends on the distributed i/o or not
  if (datacode=='G') then
     i3xcsh=istart
  else if (datacode=='D') then
     i3xcsh=nxcl+nwbl-1
  end if

  ehartreeLOC=0.d0
  do j2=1,nxc!i2=istart+1,iend !the index i2 must be changed in the distributed i/o case
     i2=j2+i3xcsh!j2=i2-istart
     do i3=1,m3
        do i1=1,m1
           pot=zf(i1,i3,j2)+correction
           ehartreeLOC=ehartreeLOC+rhopot(i1,i3,i2)*pot
           rhopot(i1,i3,i2)=pot+zfionxc(i1,i3,j2)
        end do
     end do
  end do

  if (ixc==0) then
     ehartreeLOC=ehartreeLOC*factor
  else
     ehartreeLOC=ehartreeLOC*2.d0*factor
  end if

  deallocate(zf,zfionxc)

  call timing(iproc,'PSolv_comput  ','OF')

  !gathering the data to obtain the distribution array
  !evaluating the total ehartree,eexcu,vexcu
  if (nproc.gt.1) then
     call timing(iproc,'PSolv_commun  ','ON')
     energies_mpi(1)=ehartreeLOC
     energies_mpi(2)=eexcuLOC
     energies_mpi(3)=vexcuLOC
     call MPI_ALLREDUCE(energies_mpi(1),energies_mpi(4),3,MPI_double_precision,  &
          MPI_SUM,MPI_COMM_WORLD,ierr)
     eh=energies_mpi(4)
     exc=energies_mpi(5)
     vxc=energies_mpi(6)
     call timing(iproc,'PSolv_commun  ','OF')

     if (datacode == 'G') then
        !building the array of the data to be sent from each process
        !and the array of the displacement
        call timing(iproc,'PSolv_comput  ','ON')
        do jproc=0,nproc-1
           istart=min(jproc*(md2/nproc),m2-1)
           jend=max(min(md2/nproc,m2-md2/nproc*jproc),0)
           gather_arr(jproc,1)=m1*m3*jend
           gather_arr(jproc,2)=m1*m3*istart
        end do
        call timing(iproc,'PSolv_comput  ','OF')

        !gather all the results in the same rhopot array
        istart=min(iproc*(md2/nproc),m2-1)
        call timing(iproc,'PSolv_commun  ','ON')
        call MPI_ALLGATHERV(rhopot(1,1,istart+1),gather_arr(iproc,1),MPI_double_precision,&
             rhopot,gather_arr(:,1),gather_arr(:,2),MPI_double_precision,MPI_COMM_WORLD,ierr)
        call timing(iproc,'PSolv_commun  ','OF')
     end if

  else
     eh=ehartreeLOC
     exc=eexcuLOC
     vxc=vexcuLOC
  end if

  if (allocated(gather_arr)) deallocate(gather_arr)
  if (allocated(energies_mpi)) deallocate(energies_mpi)

  if (iproc==0) write(*,*)'done.'

end subroutine PSolver



subroutine PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s)
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc
  integer, intent(out) :: n3d,n3p,n3pi,i3xcsh,i3s
  !local variables
  integer, parameter :: nordgr=4
  integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
  integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr


  !calculate the dimensions wrt the geocode
  if (geocode == 'P') then
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'S') then
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else if (geocode == 'F') then
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
  else
     stop 'PS_dim4allocation: geometry code not admitted'
  end if

  !formal start and end of the slice
  istart=iproc*(md2/nproc)
  iend=min((iproc+1)*md2/nproc,m2)

  if (datacode == 'D') then
     if (istart <= m2-1) then
        nxc=iend-istart
        if (ixc >= 11 .and. ixc <= 16 .and. geocode == 'F') then
           if (ixc==13) then
              !now the dimension of the part required for the gradient
              nwbl=min(istart,nordgr)
              nwbr=min(m2-iend,nordgr)
              nxcl=1
              nxcr=1
           else
              !now the dimension of the part required for the gradient
              if(istart<=nordgr) then
                 nxcl=istart+1
                 nwbl=0
              else
                 nxcl=nordgr+1
                 nwbl=min(nordgr,istart-nordgr)
              end if
              if(iend>=m2-nordgr+1) then
                 nxcr=m2-iend+1
                 nwbr=0
              else
                 nxcr=nordgr+1
                 nwbr=min(nordgr,m2-nordgr-iend)
              end if
           end if
        else !(for the moment GGA is not implemented in the non free BC)
           nwbl=0
           nwbr=0
           nxcl=1
           nxcr=1
        end if
        nwb=nxcl+nxc+nxcr-2
        nxt=nwbr+nwb+nwbl

        i3xcsh=nxcl+nwbl-1
        i3s=istart+1-i3xcsh
     else
        nxc=0
        nxt=0
        i3xcsh=0
        i3s=m2
     end if
     n3p=nxc
     n3d=nxt
     n3pi=n3p
  else if (datacode == 'G') then
     n3d=n03
     n3p=n03
     i3xcsh=0
     i3s=min(istart,m2-1)+1
     n3pi=max(iend-istart,0)
  else
     print *,datacode
     stop 'PS_dim4allocation: data code not admitted'
  end if

!!$  print *,'P4,iproc',iproc,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!$       'ixc,n3d,n3p,i3xcsh,i3s',ixc,n3d,n3p,i3xcsh,i3s

end subroutine PS_dim4allocation


!!***
!!****h* BigDFT/P_FFT_dimensions
!! NAME
!!   P_FFT_dimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    convolution for the periodic system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with m2 and m3 exchanged
!!
!!    n1,n2,n3 the first FFT dimensions, for the moment supposed to be even
!!
!!    md1,md2,md3 the n1,n2,n3 dimensions. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the other geometries of the Poisson Solver.
!!    The dimensions 2 and 3 are exchanged.
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3,l1A,l3A

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 !beware of the half dimension
    call fourier_dim(l1,n1)
    if (n1 == m1) then
    else
       print *,'the FFT in the x direction is not allowed'
       print *,'n01 dimension',n01
       stop
    end if
    l1=l1+1
    call fourier_dim(l2,n2)
    if (n2 == m2) then
    else
       print *,'the FFT in the z direction is not allowed'
       print *,'n03 dimension',n03
       stop
    end if
    call fourier_dim(l3,n3)
    if (n3 == m3) then
    else
       print *,'the FFT in the y direction is not allowed'
       print *,'n02 dimension',n02
       stop
    end if

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3
151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1 
 nd2=n2/2+1
 nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine P_FFT_dimensions


!!****h* BigDFT/S_FFT_dimensions
!! NAME
!!   S_FFT_dimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    convolution for the surface system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    October 2006
!!
!! SOURCE
!!
subroutine S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3,l1A,l3A

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 !beware of the half dimension
    call fourier_dim(l1,n1)
    if (n1 == m1) then
    else
       print *,'the FFT in the x direction is not allowed'
       print *,'n01 dimension',n01
       stop
    end if
    l1=l1+1
    call fourier_dim(l2,n2)
    if (n2 == m2) then
    else
       print *,'the FFT in the z direction is not allowed'
       print *,'n03 dimension',n03
       stop
    end if
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3/2
151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    goto 151
 endif


 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1

 nd3=n3/2+1
250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine S_FFT_dimensions
!!***

!!****h* BigDFT/F_FFT_dimensions
!! NAME
!!   F_FFT_pardimensions
!!
!! FUNCTION
!!    Calculate four sets of dimension needed for the calculation of the
!!    zero-padded convolution
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
!!
!!    n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!!                properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! WARNING
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!!
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2006
!!
!! SOURCE
!!
subroutine F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)
 implicit none
 integer, intent(in) :: n01,n02,n03,nproc
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3,l1A,l3A

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 l3=m3 !beware of the half dimension
 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do
 do
    call fourier_dim(l2,n2)
    if (modulo(n2,2) == 0) then
       exit
    end if
    l2=l2+1
 end do
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
151 if (nproc*(md2/nproc).lt.n2/2) then
    md2=md2+1
    goto 151
 endif

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1

250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    goto 250
 endif

end subroutine F_FFT_dimensions
!!***
