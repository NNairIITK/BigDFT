!> Calculates the FFT of the distributed kernel
!! (Based on suitable modifications of S.Goedecker routines)
!!
!! SYNOPSIS
!!   @param  zf:          Real kernel (input)
!!                        zf(i1,i2,i3)
!!   @param  zr:          Distributed Kernel FFT 
!!                        zr(2,i1,i2,i3)
!!   @param  nproc:       number of processors used as returned by MPI_COMM_SIZE
!!   @param  iproc:       [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!!   @param  n1,n2,n3:    logical dimension of the transform. As transform lengths 
!!                        most products of the prime factors 2,3,5 are allowed.
!!                        The detailed table with allowed transform lengths can 
!!                        be found in subroutine ctrig_sg
!!   @param  nd1,nd2,nd3: Dimensions of work arrays
subroutine kernelfft(n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,mpi_comm,nproc,iproc,zf,zr)
  use Poisson_Solver, only: dp
  use wrapper_mpi
  implicit none
  include 'perfdata.inc'
  !Arguments
  integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,nk1,nk2,nk3,nproc,iproc,mpi_comm
  real(dp), dimension(n1/2+1,n3/2+1,nd2/nproc), intent(in) :: zf
  real(dp), dimension(nk1,nk2,nk3/nproc), intent(inout) :: zr
  !Local variables
  character(len=*), parameter :: subname='kernelfft'
  !Maximum number of points for FFT (should be same number in fft3d routine)
  integer :: ncache,lzt,lot,ma,mb,nfft,ic1,ic2,ic3,Jp2st,J2st
  integer :: j2,j3,i1,i3,i,j,inzee,ierr,i_all,i_stat,ntrig
  real(dp) :: twopion
  !work arrays for transpositions
  real(dp), dimension(:,:,:), allocatable :: zt
  !work arrays for MPI
  real(dp), dimension(:,:,:,:,:), allocatable :: zmpi1
  real(dp), dimension(:,:,:,:), allocatable :: zmpi2
  !cache work array
  real(dp), dimension(:,:,:), allocatable :: zw
  !FFT work arrays
  real(dp), dimension(:,:), allocatable :: trig1,trig2,trig3,cosinarr
  integer, dimension(:), allocatable :: after1,now1,before1, & 
       after2,now2,before2,after3,now3,before3

  !Body

  !check input
  if (nd1.lt.n1) stop 'ERROR:nd1'
  if (nd2.lt.n2) stop 'ERROR:nd2'
  if (nd3.lt.n3/2+1) stop 'ERROR:nd3'
  if (mod(nd3,nproc).ne.0) stop 'ERROR:nd3'
  if (mod(nd2,nproc).ne.0) stop 'ERROR:nd2'
  
  !defining work arrays dimensions
  ncache=ncache_optimal
  if (ncache <= max(n1,n2,n3/2)*4) ncache=max(n1,n2,n3/2)*4
  lzt=n2
  if (mod(n2,2).eq.0) lzt=lzt+1
  if (mod(n2,4).eq.0) lzt=lzt+1
  
  ntrig=max(n1,n2,n3/2)

  !Allocations
  allocate(trig1(2,ntrig+ndebug),stat=i_stat)
  call memocc(i_stat,trig1,'trig1',subname)
  allocate(after1(7+ndebug),stat=i_stat)
  call memocc(i_stat,after1,'after1',subname)
  allocate(now1(7+ndebug),stat=i_stat)
  call memocc(i_stat,now1,'now1',subname)
  allocate(before1(7+ndebug),stat=i_stat)
  call memocc(i_stat,before1,'before1',subname)
  allocate(trig2(2,ntrig+ndebug),stat=i_stat)
  call memocc(i_stat,trig2,'trig2',subname)
  allocate(after2(7+ndebug),stat=i_stat)
  call memocc(i_stat,after2,'after2',subname)
  allocate(now2(7+ndebug),stat=i_stat)
  call memocc(i_stat,now2,'now2',subname)
  allocate(before2(7+ndebug),stat=i_stat)
  call memocc(i_stat,before2,'before2',subname)
  allocate(trig3(2,ntrig+ndebug),stat=i_stat)
  call memocc(i_stat,trig3,'trig3',subname)
  allocate(after3(7+ndebug),stat=i_stat)
  call memocc(i_stat,after3,'after3',subname)
  allocate(now3(7+ndebug),stat=i_stat)
  call memocc(i_stat,now3,'now3',subname)
  allocate(before3(7+ndebug),stat=i_stat)
  call memocc(i_stat,before3,'before3',subname)
  allocate(zw(2,ncache/4,2+ndebug),stat=i_stat)
  call memocc(i_stat,zw,'zw',subname)
  allocate(zt(2,lzt,n1+ndebug),stat=i_stat)
  call memocc(i_stat,zt,'zt',subname)
  allocate(zmpi2(2,n1,nd2/nproc,nd3+ndebug),stat=i_stat)
  call memocc(i_stat,zmpi2,'zmpi2',subname)
  allocate(cosinarr(2,n3/2+ndebug),stat=i_stat)
  call memocc(i_stat,cosinarr,'cosinarr',subname)
  if (nproc.gt.1) then
     allocate(zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc+ndebug),stat=i_stat)
     call memocc(i_stat,zmpi1,'zmpi1',subname)
  end if

  
  !calculating the FFT work arrays (beware on the HalFFT in n3 dimension)
  call ctrig_sg(n3/2,ntrig,trig3,after3,before3,now3,1,ic3)
  call ctrig_sg(n1,ntrig,trig1,after1,before1,now1,1,ic1)
  call ctrig_sg(n2,ntrig,trig2,after2,before2,now2,1,ic2)
  
  !Calculating array of phases for HalFFT decoding
  twopion=8._dp*datan(1._dp)/real(n3,dp)
  do i3=1,n3/2
     cosinarr(1,i3)= dcos(twopion*real(i3-1,dp))
     cosinarr(2,i3)=-dsin(twopion*real(i3-1,dp))
  end do
  
  !transform along z axis

  lot=ncache/(2*n3)
  if (lot.lt.1) stop 'kernelfft:enlarge ncache for z'
  
  do j2=1,nd2/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd2/nproc)+j2.le.n2) then
        do i1=1,n1,lot
           ma=i1
           mb=min(i1+(lot-1),n1)
           nfft=mb-ma+1

           !inserting real data into complex array of half lenght
           !input: I1,I3,J2,(Jp2)

           call inserthalf(n1,n3,lot,nfft,i1,zf(1,1,j2),zw(1,1,1))

           !performing FFT
           inzee=1
           do i=1,ic3
              call fftstp_sg(lot,nfft,n3/2,lot,n3/2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ntrig,trig3,after3(i),now3(i),before3(i),1)
              inzee=3-inzee
           enddo
           !output: I1,i3,J2,(Jp2)

           !unpacking FFT in order to restore correct result, 
           !while exchanging components
           !input: I1,i3,J2,(Jp2)
           call scramble_unpack(i1,j2,lot,nfft,n1,n3,nd2,nproc,nd3,zw(1,1,inzee),zmpi2,cosinarr)
           !output: I1,J2,i3,(Jp2)
        end do
     endif
  end do

  !Interprocessor data transposition
  !input: I1,J2,j3,jp3,(Jp2)
  if (nproc.gt.1) then
     !communication scheduling
     call MPI_ALLTOALL(zmpi2,2*n1*(nd2/nproc)*(nd3/nproc), &
          MPI_double_precision, &
          zmpi1,2*n1*(nd2/nproc)*(nd3/nproc), &
          MPI_double_precision,mpi_comm,ierr)
     ! output: I1,J2,j3,Jp2,(jp3)
  endif


  do j3=1,nd3/nproc
     !this condition ensures that we manage only the interesting part for the FFT
     if (iproc*(nd3/nproc)+j3.le.n3/2+1) then
        Jp2st=1
        J2st=1
        
        !transform along x axis
        lot=ncache/(4*n1)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for x'
        
        do j=1,n2,lot
           ma=j
           mb=min(j+(lot-1),n2)
           nfft=mb-ma+1

           !reverse ordering
           !input: I1,J2,j3,Jp2,(jp3)
           if (nproc.eq.1) then
              call mpiswitchPS(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi2,zw(1,1,1))
           else
              call mpiswitchPS(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw(1,1,1))
           endif
           !output: J2,Jp2,I1,j3,(jp3)

           !performing FFT
           !input: I2,I1,j3,(jp3)          
           inzee=1
           do i=1,ic1-1
              call fftstp_sg(lot,nfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ntrig,trig1,after1(i),now1(i),before1(i),1)
              inzee=3-inzee
           enddo
           !storing the last step into zt
           i=ic1
           call fftstp_sg(lot,nfft,n1,lzt,n1,zw(1,1,inzee),zt(1,j,1), & 
                ntrig,trig1,after1(i),now1(i),before1(i),1)
           !output: I2,i1,j3,(jp3)
        end do

        !transform along y axis, and taking only the first half
        lot=ncache/(4*n2)
        if (lot.lt.1) stop 'kernelfft:enlarge ncache for y'

        do j=1,nk1,lot
           ma=j
           mb=min(j+(lot-1),nk1)
           nfft=mb-ma+1

           !reverse ordering
           !input: I2,i1,j3,(jp3)
           call switchPS(nfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
           !output: i1,I2,j3,(jp3)

           !performing FFT
           !input: i1,I2,j3,(jp3)
           inzee=1
           do i=1,ic2
              call fftstp_sg(lot,nfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                   ntrig,trig2,after2(i),now2(i),before2(i),1)
              inzee=3-inzee
           enddo

           call realcopy(lot,nfft,n2,nk1,nk2,zw(1,1,inzee),zr(j,1,j3))
          
        end do
        !output: i1,i2,j3,(jp3)
     endif
  end do

  !De-allocations
  i_all=-product(shape(trig1))*kind(trig1)
  deallocate(trig1,stat=i_stat)
  call memocc(i_stat,i_all,'trig1',subname)
  i_all=-product(shape(after1))*kind(after1)
  deallocate(after1,stat=i_stat)
  call memocc(i_stat,i_all,'after1',subname)
  i_all=-product(shape(now1))*kind(now1)
  deallocate(now1,stat=i_stat)
  call memocc(i_stat,i_all,'now1',subname)
  i_all=-product(shape(before1))*kind(before1)
  deallocate(before1,stat=i_stat)
  call memocc(i_stat,i_all,'before1',subname)
  i_all=-product(shape(trig2))*kind(trig2)
  deallocate(trig2,stat=i_stat)
  call memocc(i_stat,i_all,'trig2',subname)
  i_all=-product(shape(after2))*kind(after2)
  deallocate(after2,stat=i_stat)
  call memocc(i_stat,i_all,'after2',subname)
  i_all=-product(shape(now2))*kind(now2)
  deallocate(now2,stat=i_stat)
  call memocc(i_stat,i_all,'now2',subname)
  i_all=-product(shape(before2))*kind(before2)
  deallocate(before2,stat=i_stat)
  call memocc(i_stat,i_all,'before2',subname)
  i_all=-product(shape(trig3))*kind(trig3)
  deallocate(trig3,stat=i_stat)
  call memocc(i_stat,i_all,'trig3',subname)
  i_all=-product(shape(after3))*kind(after3)
  deallocate(after3,stat=i_stat)
  call memocc(i_stat,i_all,'after3',subname)
  i_all=-product(shape(now3))*kind(now3)
  deallocate(now3,stat=i_stat)
  call memocc(i_stat,i_all,'now3',subname)
  i_all=-product(shape(before3))*kind(before3)
  deallocate(before3,stat=i_stat)
  call memocc(i_stat,i_all,'before3',subname)
  i_all=-product(shape(zmpi2))*kind(zmpi2)
  deallocate(zmpi2,stat=i_stat)
  call memocc(i_stat,i_all,'zmpi2',subname)
  i_all=-product(shape(zw))*kind(zw)
  deallocate(zw,stat=i_stat)
  call memocc(i_stat,i_all,'zw',subname)
  i_all=-product(shape(zt))*kind(zt)
  deallocate(zt,stat=i_stat)
  call memocc(i_stat,i_all,'zt',subname)
  i_all=-product(shape(cosinarr))*kind(cosinarr)
  deallocate(cosinarr,stat=i_stat)
  call memocc(i_stat,i_all,'cosinarr',subname)
  if (nproc.gt.1) then
     i_all=-product(shape(zmpi1))*kind(zmpi1)
     deallocate(zmpi1,stat=i_stat)
     call memocc(i_stat,i_all,'zmpi1',subname)
  end if

END SUBROUTINE kernelfft
subroutine realcopy(lot,nfft,n2,nk1,nk2,zin,zout)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: nfft,lot,n2,nk1,nk2
  real(dp), dimension(2,lot,n2), intent(in) :: zin
  real(dp), dimension(nk1,nk2), intent(out) :: zout
  !local variables
  integer :: i,j
 
  do i=1,nk2
     do j=1,nfft
        zout(j,i)=zin(1,j,i)
     end do
  end do

END SUBROUTINE realcopy

subroutine copyreal(n1,nk1,nfft,halfft,kernelfour)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: n1,nk1,nfft
  real(dp), dimension(2,n1,nfft), intent(in) :: halfft
  real(dp), dimension(nk1,nfft), intent(out) :: kernelfour
  !local varaibles
  integer :: ifft
  
  do ifft=1,nfft
     !if (ifft==39) write(130,'(1pe24.17)')halfft(1,:,ifft)
    call vcopy(nk1,halfft(1,1,ifft),2,kernelfour(1,ifft),1)  
  enddo
END SUBROUTINE copyreal

subroutine switchPS(nfft,n2,lot,n1,lzt,zt,zw)
   use Poisson_Solver, only: dp
   implicit none
   !Arguments
   integer, intent(in) :: nfft,n2,lot,n1,lzt
   real(dp) :: zw(2,lot,n2),zt(2,lzt,n1)
   !Local variables
   integer :: i,j
   
   do j=1,nfft
      do i=1,n2
         zw(1,j,i)=zt(1,i,j)
         zw(2,j,i)=zt(2,i,j)
      end do
   end do

END SUBROUTINE switchPS


subroutine mpiswitchPS(j3,nfft,Jp2st,J2st,lot,n1,nd2,nd3,nproc,zmpi1,zw)
   use Poisson_Solver, only: dp
   implicit none
   !Arguments
   integer, intent(in) :: j3,nfft,lot,n1,nd2,nd3,nproc
   integer, intent(inout) :: Jp2st,J2st
   real(dp) :: zmpi1(2,n1,nd2/nproc,nd3/nproc,nproc),zw(2,lot,n1)
   !Local variables
   integer :: mfft,Jp2,J2,I1
   mfft=0
   do Jp2=Jp2st,nproc
      do J2=J2st,nd2/nproc
         mfft=mfft+1
         if (mfft.gt.nfft) then
            Jp2st=Jp2
            J2st=J2
            return
         endif
         do I1=1,n1
            zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
            zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
         end do
      end do
      J2st=1
   end do

END SUBROUTINE mpiswitchPS


subroutine fill_halfft(nreal,n1,n_range,nfft,kernelreal,halfft)
  use Poisson_Solver, only: dp
  implicit none
  integer, intent(in) :: n1,n_range,nfft,nreal
  real(dp), dimension(nreal,nfft), intent(in) :: kernelreal
  real(dp), dimension(2,n1,nfft), intent(out) :: halfft
  !local variables
  integer :: i,ifft,ntot

  ntot=min(n1/2,n_range)
  halfft(:,:,:)=0.0_dp
  !if (n1 > 2*n_range+1) stop 'n1 > 2*n_range+1'
  do ifft=1,nfft
     do i=0,ntot
        !print *,n1/2+i+1,n1,n_range
        halfft(1,n1/2+i+1,ifft)=kernelreal(i+1,ifft)/real(n1,dp)
     end do
     do i=1,ntot
        !print *,n1/2-i+1,n1,n_range
        halfft(1,n1/2-i+1,ifft)=kernelreal(i+1,ifft)/real(n1,dp)
     end do
     if (ifft==39) halfft(1,:,ifft)=11.0_dp
  end do

END SUBROUTINE fill_halfft
