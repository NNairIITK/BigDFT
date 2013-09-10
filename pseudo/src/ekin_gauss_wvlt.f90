!> @file
!! Part of the pseudo program (pseudopotential generation)
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by Matthias Krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> This procedure is supposed to call a minimum of
!! routines taken from an older version of BigDFT
!! to calculate the kinetic energy in a wavelet
!! basis. Input shall be a given gaussian
!! representation of the wavefunction.
subroutine  ekin_wvlt(verbose,iproc,nproc,ng,ngmx,&
   noccmx,lmx,nspin,nsmx,&
   nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,&
   xp_in,psi_in,occup_in, ekin_pen,time) 

   implicit none

   !Arguments
   integer, intent(in) :: iproc,nproc,ng,ngmx,noccmx,lmx,nspin,nsmx,nhpow
   integer, intent(out) :: nhgrid
   !> The shape of the psi and xp array as defined in the main program pseudo.f90
   real(kind=8), dimension(0:ngmx,noccmx,lmx,nsmx), intent(in) :: psi_in
   real(kind=8), dimension(0:ngmx), intent(in) :: xp_in
   real(kind=8), dimension(noccmx,lmx,nsmx), intent(in) :: occup_in
   real(kind=8), dimension(3), intent(inout) :: time
   real(kind=8), intent(in) :: hgridmin, hgridmax,crmult,frmult,ampl
   real(kind=8), intent(out) :: ekin_pen
   logical, intent(in) :: verbose
   !Local variables
   real(kind=8), parameter :: eps_mach=1.d-12,onem=1.d0-eps_mach
   integer, dimension(lmx,nsmx) :: nl
   real(kind=8), dimension(0:ngmx) :: xp
   real(kind=8), dimension(:,:,:), allocatable :: ekin_plot
   character(len=20) :: frmt

   !Arguments of createWavefunctionsDescriptors
   integer :: n1,n2,n3,norb
   logical :: parallel 
   real(kind=8) :: hgrid
   integer, dimension(:,:,:), allocatable :: ibyz_c, ibxz_c, ibxy_c
   integer, dimension(:,:,:), allocatable :: ibyz_f, ibxz_f, ibxy_f
   ! wavefunction 
   ! real(kind=8), allocatable :: psi(:,:)
   ! wavefunction gradients
   ! arrays for DIIS convergence accelerator
   real(kind=8) :: t
   logical, dimension(:,:,:), allocatable :: logrid_c, logrid_f

   !Arguments of createAtomicOrbitals
   integer :: norbe, norbep, ngx
   integer :: i,ihgrid,iocc,iorb,ispin,iprocs,iserial
   integer :: ierr,l,nprocs,nfull,ncoeff
   integer :: nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
   
   real(kind=8), dimension(3) :: rxyz
   real(kind=8) :: alat1, alat2, alat3,cxmin,cxmax,cymin,cymax,czmin,czmax
   real(kind=8) :: drxyz, ekin_mypen,rad,wghth,wghthsum
   
   include 'mpif.h'
   parallel=(nproc>1)
   ngx=ngmx
   !write(6,*)'parallel=',parallel
   
    
   ! Each process needs the wfn of configuration 1 from process zero 
   if (parallel)then
      ncoeff=product(shape(psi_in))
      call cpu_time(t)
      time(2)=time(2)-t
      call MPI_BCAST(psi_in,ncoeff,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call cpu_time(t)
      time(2)=time(2)+t
      if (ierr.ne.0) write(*,*)'MPI_BCAST of PSI: ierr',ierr
   end if
   
   ! Don't read atomic positions from posinp.xyz.
   ! Nonzero coordinates are available with random shifts.
   rxyz=0d0


   !get number of (semicore and) valence orbitals from occup that,
   !as we only want to enforce softness for the occupied states.
   !Also count nl, the occupied number of orbitals per l&s channel

   !To do so, take occup from process zero = configuration one := ground state
   !and broadcast the result.
   nl=0
   norb=0
   if (iproc==0)then
      do l=1,lmx
         do ispin=1,nspin
            do iocc=1,noccmx
               !convention from atom.f : occupation numbers are not written out
               !as zeroes, but such that multiplying by 1d20 occupies the orbital 
               if(occup_in(iocc,l,ispin)>1d-10)then
                  norb=norb+1
                  nl(l,ispin)=nl(l,ispin)+1
               end if
            end do
            !write(6,*)'debug: l,s,nl(l,s)',l,ispin,nl(l,ispin)
        end do
     end do
  end if

if(parallel)then 
    ncoeff=product(shape(nl))      
         call cpu_time(t)
         time(2)=time(2)-t
    call MPI_BCAST(nl,ncoeff,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call cpu_time(t)
         time(2)=time(2)+t
end if
norb=sum(nl)
norbe=norb! we do not treat unoccupied orbitals, thus norbe=norb
!write(6,*)'debug: occupied orbitals and norbe:',norb,norbe


if(verbose)then
   write(6,'(1x,a)') 'Penalty contributions from wavelet transformtaions' 
   write(6,'(1x,a)') '__________________________________________________'
   !output some additional information. This is not really needed,
   !only done once after the amoeba has converged/exited.
   write(6,'(/,1x,a)') 'Occupation numbers and orbital indices for wavelet transformation'
   write(6,'(1x,a)')   '          l          inl          s   -> orbital   occupation number' 
   iorb=0
   do l=1,lmx
      do ispin=1,nspin
          do iocc=1,noccmx
              if(occup_in(iocc,l,ispin)>1d-10)then
              iorb=iorb+1
              write(6,*)l, iocc, ispin, iorb, occup_in(iocc,l,ispin)
              end if
          end do
      end do
   end do
   write(6,*)
end if

! For time profiling, we consider the following section as calls to wavelet libraries

call cpu_time(t)
time(3)=time(3)-t


   ! Treat nhgrid different values for the grid spacing hgrid in the
   ! range from hgridmin to hgridmax. Do this in parallel and write 
   ! a warning if the number of processes does not make sense.

   ! Here we have different conventions how to distribute the work.
   ! If there is exactly one process per transformed orbital,
   ! the hgrid samples are computed in a loop, just as in the serial
   ! case. Otherwise, all hgrid samples are treated in parallel,
   ! assigning a part of the processes to each sample. 

if(parallel.and.norb.ne.nproc)then
   
   if(nproc<nhgrid)then
       write(6,*)'NOTICE: Not enough MPI processes to treat',&
                 nhgrid,'grid samples. Setting nhgrid=nproc=',nproc
       nhgrid=nproc
   end if
   
   if(nproc>nhgrid*norb)then
       write(6,*)'NOTICE: Found more MPI processes than can be used for',&
                 nhgrid,'grid samples and',norb,'orbitals.'
       nhgrid=nproc/norb
       write(6,*)'Raising the number of samples to',nhgrid
   end if
   
   ! Assign at least one process per hgrid sample.
   ! Distrubute the remaining processes such that 
   ! the finest grid samples, which are computationally
   ! costly, have norb processes each. 
   ! 
   
   if(iproc<nhgrid)then
       !first at least one process goes to each  sample
       ihgrid = iproc 
   else
       !then assign the remaining processeses to increasing values of hgrid  
       ihgrid = (iproc-nhgrid)/(norb-1) ! 
       !note: norb is never one when nproc > nhgrid, because of auto adjust
   end if
   
   
   if(ihgrid<0.or.ihgrid>nhgrid-1)&
      write(6,*)'WARNING, BUG detected: Invalid value for ihgrid:',ihgrid
   
   ! Linear interpolation of the sample value
   hgrid = hgridmin + ihgrid/(nhgrid-1d0)*(hgridmax-hgridmin)   
   if(nhgrid==1)hgrid=hgridmin
   
   
   ! We have a 2nd level of parallelization, where several processes but not all
   ! perform the transformation of the gaussians to the wavelet basis for a given
   ! hgrid sample. Thus, we need new variables equivalent to nproc and iproc:
   ! Assign nprocs and iprocs among those processes that treat this sample.
   
   ! (modulo function is slow, but this seems easier to read)
   
   if(iproc<nhgrid)then
      iprocs = 0
   else
      iprocs = 1+mod((iproc-nhgrid),norb-1)
   end if
   
   if(norb==1)then
      nprocs=1
   else
      !get the number of samples that are assigned norb processes each
      nfull = (nproc-nhgrid)/(norb-1)
      !thus ihgrid==nfull is the only sample with nprocs other than 1 or norb
   
      !if(ihgrid<nfull)&
      nprocs = norb
      if(ihgrid==nfull) nprocs = 1+ mod((nproc-nhgrid),(norb-1))
      if(ihgrid > nfull) nprocs = 1
   end if

else
   !Either the serial case or norb is equal to nproc
   !In both cases, the hgrid are sampled in a loop
   iprocs=iproc
   nprocs=nproc
   ihgrid=0
   hgrid=hgridmin
end if

allocate(ekin_plot(norbe,nhgrid,2))
ekin_plot=0d0
ekin_pen=0d0

!write(6,*)'DEBUG: ihgrid,hgrid,iprocs,nprocs', ihgrid,hgrid,iprocs,nprocs
 
do iserial=1,nhgrid
 !this loop is for the treatment of several hgrid samples in the serial case 
 !or when nproc = norb. Otherwise, exit this loop after the first iteration.
 


! determine size alat of overall simulation cell
  call system_size(rxyz,2.0*crmult, &
       cxmin,cxmax,cymin,cymax,czmin,czmax)
  alat1=(cxmax-cxmin+ampl)
  alat2=(cymax-cymin+ampl)
  alat3=(czmax-czmin+ampl)

!!write(16,*)'DEBUG: box and radii',alat1,alat2,alat3,crmult,frmult

! shift atomic positions such that molecule is inside cell
     rxyz(1)=rxyz(1)-cxmin
     rxyz(2)=rxyz(2)-cymin
     rxyz(3)=rxyz(3)-czmin

        if(ampl>0)then
           call random_number(drxyz)
           rxyz= rxyz+ampl*drxyz
        end if
     

! set up the grid for the wavelet basis.


! grid sizes n1,n2,n3
  n1=int(alat1/hgrid)
  if (mod(n1+1,4).eq.0) n1=n1+1
  n2=int(alat2/hgrid)
  if (mod(n2+1,8).eq.0) n2=n2+1
  n3=int(alat3/hgrid)
  alat1=n1*hgrid ; alat2=n2*hgrid ; alat3=n3*hgrid
!!  if (iproc.eq.0) then
!!     write(*,'(1x,a,3(1x,i0))') 'n1,n2,n3',n1,n2,n3
!!     write(*,'(1x,a,3(1x,i0))') 'total number of grid points',(n1+1)*(n2+1)*(n3+1)
!!     write(*,'(1x,a,3(1x,1pe12.5))') 'simulation cell',alat1,alat2,alat3
!!  endif

! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 ; nfl2=n2 ; nfl3=n3
  nfu1=0 ; nfu2=0 ; nfu3=0
!  do iat=1,nat
     rad=frmult
     nfl1=min(nfl1,int(onem+(rxyz(1)-rad)/hgrid)) ; nfu1=max(nfu1,int((rxyz(1)+rad)/hgrid))
     nfl2=min(nfl2,int(onem+(rxyz(2)-rad)/hgrid)) ; nfu2=max(nfu2,int((rxyz(2)+rad)/hgrid))
     nfl3=min(nfl3,int(onem+(rxyz(3)-rad)/hgrid)) ; nfu3=max(nfu3,int((rxyz(3)+rad)/hgrid))
!  enddo
!!!  if (iproc.eq.0) then
!!!     write(*,'(1x,a,2(1x,i0))') 'nfl1,nfu1 ',nfl1,nfu1
!!!     write(*,'(1x,a,2(1x,i0))') 'nfl2,nfu2 ',nfl2,nfu2
!!!     write(*,'(1x,a,2(1x,i0))') 'nfl3,nfu3 ',nfl3,nfu3
!!!  endif

! Create wavefunctions descriptors and allocate them
! To be  more precise: We only need the bounds for the convolution,
! as the compressed wavelet form of the wfns is never used.
  allocate(ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2))
  allocate(ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2))
!
!
  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3))
  allocate(logrid_f(0:n1,0:n2,0:n3))

  ! coarse grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,rxyz, &
       2*crmult,hgrid,logrid_c)
  !     2*rprb,crmult,hgrid,logrid_c)
  call bounds(n1,n2,n3,logrid_c,ibyz_c,ibxz_c,ibxy_c)

  ! fine grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,rxyz, &
       frmult,hgrid,logrid_f)
!      rprb*frmult,hgrid,logrid_f)
!   if (iproc.eq.0) write(*,'(1x,a,2(1x,i10))') 'orbitals have fine   segment, elements',nseg_f,7*nvctr_f
  call bounds(n1,n2,n3,logrid_f,ibyz_f,ibxz_f,ibxy_f)


  deallocate(logrid_c,logrid_f)


!  if(norbe/=norb)write(6,*)'(norbe/=norb)'



! atomkin uses a different form of the gaussians than gatom does
  do i=0,ng
     xp(i)=sqrt(0.5d0/xp_in(i))
  end do


! Call the routine that performs the wavelet transformation and convolution for the kinetic energy
 call createAtomicOrbitals(iprocs, nprocs, nspin,  &
         rxyz, norbe, norbep,  ngx, xp, psi_in, ng+1, nl, &
         n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, &
         ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
         ekin_mypen,ekin_plot(1:norbe,ihgrid+1,1))


  !multiply the contribution from this hgrid value according to the power law from input.dat
  !sum up the hgrid**(-nhpow) to normalize these weights.
  wghthsum=0d0
  if(nhgrid>1)then
     do i=0,nhgrid-1
        wghthsum=wghthsum+1d0/(hgridmin + i/(nhgrid-1d0)*(hgridmax-hgridmin))**nhpow
     end do
     wghth=1d0/(hgrid**nhpow)/wghthsum
  else
     wghth=1d0
  end if

  ekin_pen=ekin_pen+ekin_mypen*wghth
!!  write(16,*)'DEBUG ekin_mypen,wght,wghtsum,ekin_pen',ekin_mypen,wghth,wghtsum,ekin_pen
  !Note: So far, this is the contribution from the hgrid and orbitals assigned to this
  !process. The sum over all hgrid and orbitals is done in the calling penalty subroutine.
   
  deallocate(ibxz_c,ibyz_c,ibxy_c,ibxz_f,ibyz_f,ibxy_f)

  if(parallel.and.norb.ne.nproc)exit
  !assign new hgrid value in the serial case
  ihgrid= iserial
  hgrid = hgridmin + ihgrid/(nhgrid-1d0)*(hgridmax-hgridmin)   
  if(verbose)write(6,'(5x,2(a,i3),a)')'(hgrid ',iserial,' of ',nhgrid,' done.)'
  
end do ! serial case: loop over hgrid samples


! for time profiling, we consider the calls to wavelet libraries done here
call cpu_time(t)
time(3)=time(3)+t

if(verbose)then
   !write out some info and do the plotting. 
   !MPI procedure to gather all information from other processes: ekin_plot(orbitals, hgrid, buff)
   if(parallel)then
      call cpu_time(t)
      time(2)=time(2)-t
      call MPI_ALLREDUCE(ekin_plot(1,1,1),ekin_plot(1,1,2),norb*nhgrid,&
                         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call cpu_time(t)
      time(2)=time(2)+t
   else
      ekin_plot(:,:,2)=ekin_plot(:,:,1)
   end if
 
   write(6,*)
   if (iproc==0)then
       open(101,file='dEkin.orb.dat')
       write(101,'(a)',advance='no')'# hgrid            dE_kin(gauss-wvlt) weighted by hgrid '
   end if
   write(6,      '(a)',advance='no')'  hgrid            dE_kin(gauss-wvlt) weighted by hgrid '
   do iorb=1,norb
      if(iproc==0)&
      write(101,'(a,i3,8x)',advance='no') 'orbital',iorb
      write(  6,'(a,i3,8x)',advance='no') 'orbital',iorb
   end do
   if(iproc==0)&
   write(101,*)
   write(  6,*)

   do i=0,nhgrid-1
      if(nhgrid==1)then
         hgrid=hgridmin
         wghth=1d0
      else
         hgrid=hgridmin + i/(nhgrid-1d0)*(hgridmax-hgridmin)
         wghth=1d0/(hgrid**nhpow)/wghthsum
      end if
      ekin_mypen = sqrt(sum(ekin_plot(:,i+1,2)))
      write(frmt,'(a,i2,a)')'(',norb+3,'e18.10)'
      write(6,  frmt)hgrid, ekin_mypen,ekin_mypen*wghth, sqrt(ekin_plot(:,i+1,2))
      if(iproc==0)&
      write(101,frmt)hgrid, ekin_mypen,ekin_mypen*wghth, sqrt(ekin_plot(:,i+1,2))
   end do
   if(iproc==0) close(101)
end if
      
deallocate(ekin_plot)

end subroutine ekin_wvlt


subroutine createAtomicOrbitals(iproc, nproc, nspin,&
       rxyz, norbe, norbep, ngx, xp,psi_in, ng, nl, &
       n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,&
       ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
       ekin_mypen,ekin_plot)

  implicit none
  !Arguments
  integer, intent(in) :: norbe, norbep, ngx, iproc, nproc
  integer, intent(in) :: n1, n2, n3 
  integer, intent(in) :: nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
  integer, intent(in) :: ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  integer, intent(in) :: ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  real(kind=8), intent(in) :: hgrid 
  !character(len = 20) :: pspatomnames(npsp)
  !character(len = 20) :: atomnames
  !integer :: ng, nl(4)
  integer, intent(in) :: ng, nl(5,2)
  real(kind=8), dimension(3), intent(in) :: rxyz
  real(kind=8), dimension(ngx), intent(in) :: xp
  !Local variables
  integer, parameter :: nterm_max=3
  integer, dimension(nterm_max) :: lx,ly,lz
  real(kind=8), dimension(nterm_max) :: fac_arr
  integer :: ispin,nspin, iorb, jorb,  i,  inl, l, m,  nterm
  real(kind=8) :: rx, ry, rz,  scpr, ekin, ekgauss
  real(kind=8), dimension(0:32,7,5,2) :: psi_in
  logical :: myorbital
  real(kind=8), dimension(:), allocatable :: psig, psigp, psiatn, psiat
  real(kind=8), dimension(norbe) :: ekin_plot
  real(kind=8) :: ekin_mypen
  real(kind=8), external :: dnrm2

  allocate (psiatn(ngx), psiat(ngx))
  allocate(psig (8*(1+n1)*(1+n2)*(1+n3)),psigp(8*(1+n1)*(1+n2)*(1+n3) )) 
  ekin_mypen=0.d0
  ekin_plot=0d0

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)

     iorb=0! jorb starts from iproc*norbp+1 

     rx=rxyz(1)
     ry=rxyz(2)
     rz=rxyz(3)

     do l=1,4
      do ispin=1,nspin
!       write(6,*)'nl(l,ispin)',l,ispin,nl(l,ispin)
        do inl=1,nl(l,ispin)!  zero for unoccupied l/s - channels
            !pick coefficients of psi from gatom to form psiat
            ! psi(gridpoint, iocc=inl, l, s) --> psiat(gridpoint)
          
            !psiat=0d0
            do i=1,ng+1
               psiat(i)=psi_in(i-1,inl,l,ispin)
            end do


            call atomkin(l-1,ng,xp(1),psiat,psiatn,ekgauss)
            do m=1,1! m=1, 2*l-1  ! for now, do not treat several m
              iorb=iorb+1
!              write(6,*)'DEBUG:  atomkin, l, s, iocc',eki,l,ispin,inl
              jorb=iorb-iproc*norbep ! index among orbitals for this process
!             only perform the costly transform to 3d space for orbitals of this process
              if (myorbital(iorb,norbe,iproc,nproc)) then
                 !this will calculate the proper spherical harmonics
                 call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
                 !fac_arr=1.d0
!                !perform the 3d transform to daubechie wavelets
                 call crtonewave(n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,xp(1),psiatn,&
                      rx,ry,rz,hgrid, &
                      0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,psig)
!                Normalize the resulting wavelets
!                We do not use the compressed form, thus we can call BLAS routines
!                instead of their  counterparts for the compressed wavelet representation
                 scpr=dnrm2(8*(n1+1)*(n2+1)*(n3+1),psig,1)
!                 call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr)
!                write(*,'(1x,a,2(a3,i1),a16,i4,i4,1x,e14.7)')&
!                     'my ATOMIC INPUT ORBITAL ',&
!                     'l=',l,'spin=',ispin,'iorb,jorb,1-norm',iorb,jorb,1d0-scpr
                 !scpr=1.d0/sqrt(scpr)
                 scpr=1.d0/scpr
!                call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
                 call dscal(8*(n1+1)*(n2+1)*(n3+1),scpr,psig,1)
!                call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr)
                 scpr=dnrm2(8*(n1+1)*(n2+1)*(n3+1),psig,1)
                 if(abs(scpr-1d0)>1d-8)write(6,*) 'Bad normalization for orbital ',iorb,scpr
!                calculate the kinetic energy in the wavelet representation here
                 call ConvolkineticP(n1,n2,n3,  &
                      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,ibyz_c,ibxz_c,&
                      ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
                      psig,ekin)
                 ekin_plot(iorb)=(ekin-ekgauss)**2
                 ekin_mypen=ekin_mypen+ekin_plot(iorb)
!                DEBUG Line for high precision output
!                write(6,'(a,i3,3e24.16)')' iorb, ekin_gauss, ekin_wvlt,dE',&
!                                           iorb,ekgauss,ekin,ekgauss-ekin
              endif
            end do !m, currently one set sonstant at 1
        end do ! nl(l) 
      end do ! spin
     end do ! l-chanel

  deallocate(psig ,psigp,psiat,psiatn)

!  end do

END SUBROUTINE createAtomicOrbitals


!> Returns an input guess orbital that is a Gaussian centered at a Wannier center
!! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
!! in the arrays psi_c, psi_f
!! @todo Use the routine from bigdft program.
subroutine crtonewave(n1,n2,n3,nterm,ntp,lx,ly,lz,fac_arr,xp,psiat,rx,ry,rz,hgrid, &
                   nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,psig)
        implicit real(kind=8) (a-h,o-z)
        parameter(nw=16000)
        dimension xp(nterm),psiat(nterm),fac_arr(ntp)
        dimension lx(ntp),ly(ntp),lz(ntp)
        real(kind=8), allocatable, dimension(:,:) :: wprojx, wprojy, wprojz
        real(kind=8), allocatable, dimension(:,:) :: work
        !new variable
        real(kind=8)  ::psig(0:n1,2,0:n2,2,0:n3,2)
        psig=0.d0

        allocate(wprojx(0:n1,2),wprojy(0:n2,2),wprojz(0:n3,2),work(0:nw,2))

!       write(6,*)'DEBUG:n1,n2,n3'
!       write(6,*)       n1,n2,n3 
!       write(6,*)'DEBUG:nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f'
!       write(6,*)       nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f 

      iterm=1
      itp=1
        gau_a=xp(iterm)
        n_gau=lx(itp)
        CALL GAUSS_TO_DAUB(hgrid,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        n_gau=ly(itp)
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        n_gau=lz(itp)
        CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)


!      IMPORTANT modification: Here was an implicit compress step of the wfn which has been eliminated
!                              by storing the representation which would result from the following 
!                              uncompres step in applylocpotkinone

! First term: coarse projector components
        do i3=nl3_c,nu3_c
           do i2=nl2_c,nu2_c
              do i1=nl1_c,nu1_c
                 psig(i1,1,i2,1,i3,1)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
              enddo
           enddo
        enddo

! First term: fine projector components
        do i3=nl3_f,nu3_f
           do i2=nl2_f,nu2_f
              do i1=nl1_f,nu1_f
                 psig(i1,2,i2,1,i3,1)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
                 psig(i1,1,i2,2,i3,1)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
                 psig(i1,2,i2,2,i3,1)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
                 psig(i1,1,i2,1,i3,2)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
                 psig(i1,2,i2,1,i3,2)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
                 psig(i1,1,i2,2,i3,2)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
                 psig(i1,2,i2,2,i3,2)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
              enddo
           enddo
        enddo

        do iterm=2,nterm
           gau_a=xp(iterm)
           n_gau=lx(itp)
           CALL GAUSS_TO_DAUB(hgrid,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
           n_gau=ly(itp)
           CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
           n_gau=lz(itp)
           CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

           ! First term: coarse projector components
           do i3=nl3_c,nu3_c
              do i2=nl2_c,nu2_c
                 do i1=nl1_c,nu1_c
                    psig(i1,1,i2,1,i3,1)=psig(i1,1,i2,1,i3,1)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
                 enddo
              enddo
           enddo

! First term: fine projector components
           do i3=nl3_f,nu3_f
              do i2=nl2_f,nu2_f
                 do i1=nl1_f,nu1_f
                    psig(i1,2,i2,1,i3,1)=psig(i1,2,i2,1,i3,1)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
                    psig(i1,1,i2,2,i3,1)=psig(i1,1,i2,2,i3,1)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
                    psig(i1,2,i2,2,i3,1)=psig(i1,2,i2,2,i3,1)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
                    psig(i1,1,i2,1,i3,2)=psig(i1,1,i2,1,i3,2)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
                    psig(i1,2,i2,1,i3,2)=psig(i1,2,i2,1,i3,2)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
                    psig(i1,1,i2,2,i3,2)=psig(i1,1,i2,2,i3,2)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
                    psig(i1,2,i2,2,i3,2)=psig(i1,2,i2,2,i3,2)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
                 enddo
              enddo
           enddo

        end do

        do itp=2,ntp

        do iterm=1,nterm
           gau_a=xp(iterm)
           n_gau=lx(itp)
           CALL GAUSS_TO_DAUB(hgrid,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
           n_gau=ly(itp)
           CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
           n_gau=lz(itp)
           CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

           ! First term: coarse projector components
           do i3=nl3_c,nu3_c
              do i2=nl2_c,nu2_c
                 do i1=nl1_c,nu1_c
                    psig(i1,1,i2,1,i3,1)=psig(i1,1,i2,1,i3,1)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
                 enddo
              enddo
           enddo

! First term: fine projector components
           do i3=nl3_f,nu3_f
              do i2=nl2_f,nu2_f
                 do i1=nl1_f,nu1_f
                    psig(i1,2,i2,1,i3,1)=psig(i1,2,i2,1,i3,1)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
                    psig(i1,1,i2,2,i3,1)=psig(i1,1,i2,2,i3,1)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
                    psig(i1,2,i2,2,i3,1)=psig(i1,2,i2,2,i3,1)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
                    psig(i1,1,i2,1,i3,2)=psig(i1,1,i2,1,i3,2)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
                    psig(i1,2,i2,1,i3,2)=psig(i1,2,i2,1,i3,2)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
                    psig(i1,1,i2,2,i3,2)=psig(i1,1,i2,2,i3,2)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
                    psig(i1,2,i2,2,i3,2)=psig(i1,2,i2,2,i3,2)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
                 enddo
              enddo
           enddo

        end do


        end do


!This part is eliminated
!!wavefunction compression
!! coarse part
!    do iseg=1,nseg_c
!          jj=keyv_c(iseg)
!          j0=keyg_c(1,iseg)
!          j1=keyg_c(2,iseg)
!             ii=j0-1
!             i3=ii/((n1+1)*(n2+1))
!             ii=ii-i3*(n1+1)*(n2+1)
!             i2=ii/(n1+1)
!             i0=ii-i2*(n1+1)
!             i1=i0+j1-j0
!      do i=i0,i1
!            psi_c(i-i0+jj)=psig_c(i,i2,i3)
!          enddo
!        enddo
!
!! fine part
!    do iseg=1,nseg_f
!          jj=keyv_f(iseg)
!          j0=keyg_f(1,iseg)
!          j1=keyg_f(2,iseg)
!             ii=j0-1
!             i3=ii/((n1+1)*(n2+1))
!             ii=ii-i3*(n1+1)*(n2+1)
!             i2=ii/(n1+1)
!             i0=ii-i2*(n1+1)
!             i1=i0+j1-j0
!      do i=i0,i1
!            psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
!            psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
!            psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
!            psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
!            psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
!            psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
!            psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
!          enddo
!        enddo

          deallocate(wprojx,wprojy,wprojz,work)

END SUBROUTINE



!> y = y + (kinetic energy operator)x 
!! This version does not compute the y array, as only the scalar ekin
!! is needed and a bit of time can be saved this way. Remove all !y
!! comments to undo this modification
subroutine ConvolkineticP(n1,n2,n3, &
               nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
               hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x,ekin)
    implicit real(kind=8) (a-h,o-z)
    logical :: firstcall=.true.
    integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
    dimension x(0:n1,2,0:n2,2,0:n3,2)
    dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
    dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)


    parameter(lowfil=-14,lupfil=14)
    dimension a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
    scale=-.5d0/hgrid**2

!---------------------------------------------------------------------------
! second derivative filters for Daubechies 16
!  <phi|D^2|phi_i>
    a(0)=   -3.5536922899131901941296809374d0*scale
    a(1)=    2.2191465938911163898794546405d0*scale
    a(2)=   -0.6156141465570069496314853949d0*scale
    a(3)=    0.2371780582153805636239247476d0*scale
    a(4)=   -0.0822663999742123340987663521d0*scale
    a(5)=    0.02207029188482255523789911295638968409d0*scale
    a(6)=   -0.409765689342633823899327051188315485d-2*scale
    a(7)=    0.45167920287502235349480037639758496d-3*scale
    a(8)=   -0.2398228524507599670405555359023135d-4*scale
    a(9)=    2.0904234952920365957922889447361d-6*scale
    a(10)=  -3.7230763047369275848791496973044d-7*scale
    a(11)=  -1.05857055496741470373494132287d-8*scale
    a(12)=  -5.813879830282540547959250667d-11*scale
    a(13)=   2.70800493626319438269856689037647576d-13*scale
    a(14)=  -6.924474940639200152025730585882d-18*scale
    do i=1,14
        a(-i)=a(i)
    enddo
!  <phi|D^2|psi_i>
    c(-14)=     -3.869102413147656535541850057188d-18*scale
    c(-13)=      1.5130616560866154733900029272077362d-13*scale
    c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
    c(-11)=     -5.96264938781402337319841002642d-9*scale
    c(-10)=     -2.1656830629214041470164889350342d-7*scale
    c(-9 )=      8.7969704055286288323596890609625d-7*scale
    c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
    c(-7 )=      0.00021710795484646138591610188464622454d0*scale
    c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
    c(-5 )=      0.00713761218453631422925717625758502986d0*scale
    c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
    c(-3 )=      0.14327329352510759457155821037742893841d0*scale
    c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
    c(-1 )=      0.65703074007121357894896358254040272157d0*scale
    c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
    c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
    c( 2 )=      0.63457035267892488185929915286969303251d0*scale
    c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
    c( 4 )=      0.23370490631751294307619384973520033236d0*scale
    c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
    c( 6 )=      0.0080872029411844780634067667008050127d0*scale
    c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
    c( 8 )=      0.00005075807947289728306309081261461095d0*scale
    c( 9 )=     -4.62561497463184262755416490048242d-6*scale
    c( 10)=      6.3919128513793415587294752371778d-7*scale
    c( 11)=      1.87909235155149902916133888931d-8*scale
    c( 12)=      1.04757345962781829480207861447155543883d-10*scale
    c( 13)=     -4.84665690596158959648731537084025836d-13*scale
    c( 14)=      1.2392629629188986192855777620877d-17*scale
!  <psi|D^2|phi_i>
    do i=-14,14
        b(i)=c(-i)
    enddo
    !<psi|D^2|psi_i>
    e(0)=   -24.875846029392331358907766562d0*scale
    e(1)=   -7.1440597663471719869313377994d0*scale
    e(2)=   -0.04251705323669172315864542163525830944d0*scale
    e(3)=   -0.26995931336279126953587091167128839196d0*scale
    e(4)=    0.08207454169225172612513390763444496516d0*scale
    e(5)=   -0.02207327034586634477996701627614752761d0*scale
    e(6)=    0.00409765642831595181639002667514310145d0*scale
    e(7)=   -0.00045167920287507774929432548999880117d0*scale
    e(8)=    0.00002398228524507599670405555359023135d0*scale
    e(9)=   -2.0904234952920365957922889447361d-6*scale
    e(10)=   3.7230763047369275848791496973044d-7*scale
    e(11)=   1.05857055496741470373494132287d-8*scale
    e(12)=   5.8138798302825405479592506674648873655d-11*scale
    e(13)=  -2.70800493626319438269856689037647576d-13*scale
    e(14)=   6.924474940639200152025730585882d-18*scale
    do i=1,14
        e(-i)=e(i)
    enddo


  if (firstcall) then

! (1/2) d^2/dx^2
    mflop1=0
    do i3=0,n3
    do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                do l=max(-i1,lowfil),min(lupfil,n1-i1)
                    mflop1=mflop1+4
                enddo
                    mflop1=mflop1+4
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                do l=max(-i1,lowfil),min(lupfil,n1-i1)
                    mflop1=mflop1+3
                enddo
        enddo
    enddo
    enddo

! + (1/2) d^2/dy^2
    mflop2=0
    do i3=0,n3
    do i1=0,n1
        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
                do l=max(-i2,lowfil),min(lupfil,n2-i2)
                    mflop2=mflop2+4
                enddo
                    mflop2=mflop2+4
        enddo

        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
                do l=max(-i2,lowfil),min(lupfil,n2-i2)
                    mflop2=mflop2+3
                enddo
        enddo
    enddo
    enddo


! + (1/2) d^2/dz^2
    mflop3=0
    do i2=0,n2
    do i1=0,n1
        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                do l=max(-i3,lowfil),min(lupfil,n3-i3)
                    mflop3=mflop3+4
                enddo
                    mflop3=mflop3+4
        enddo

        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
                do l=max(-i3,lowfil),min(lupfil,n3-i3)
                    mflop3=mflop3+3
                enddo
        enddo
    enddo
    enddo

! wavelet part
 ! (1/2) d^2/dx^2
    nflop1=0
    do i3=nfl3,nfu3
        do i2=nfl2,nfu2
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                    nflop1=nflop1+26
                enddo
                    nflop1=nflop1+21
            enddo
        enddo
    enddo

 ! + (1/2) d^2/dy^2
    nflop2=0
    do i3=nfl3,nfu3
    do i1=nfl1,nfu1
       do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
                do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                    nflop2=nflop2+26
                enddo
                    nflop2=nflop2+21
       enddo
    enddo
    enddo

 ! + (1/2) d^2/dz^2
    nflop3=0
    do i2=nfl2,nfu2
    do i1=nfl1,nfu1
       do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
                do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                    nflop3=nflop3+26
                enddo
                    nflop3=nflop3+21
       enddo
    enddo
    enddo

    firstcall=.false.
    endif


!---------------------------------------------------------------------------

     ekin=0.d0

! Scaling function part

       call system_clock(ncount0,ncount_rate,ncount_max)


! (1/2) d^2/dx^2
    do i3=0,n3
    do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                t111=0.d0 ; s111=0.d0
                do l=max(-i1,lowfil),min(lupfil,n1-i1)
                    t111=t111 + x(i1+l,1,i2,1,i3,1)*a(l)
                    s111=s111 + x(i1+l,2,i2,1,i3,1)*b(l)
                enddo
!y              y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
                ekin=ekin+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                t211=0.d0
                do l=max(-i1,lowfil),min(lupfil,n1-i1)
                    t211=t211 + x(i1+l,1,i2,1,i3,1)*c(l)
                enddo
!y              y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
                ekin=ekin+t211*x(i1,2,i2,1,i3,1)
        enddo
    enddo
    enddo

       call system_clock(ncount1,ncount_rate,ncount_max)
       tel=dble(ncount1-ncount0)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:x',tel,1.d-6*mflop1/tel

! + (1/2) d^2/dy^2
    do i3=0,n3
    do i1=0,n1
        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
                t111=0.d0 ; s111=0.d0
                do l=max(-i2,lowfil),min(lupfil,n2-i2)
                    t111=t111 + x(i1,1,i2+l,1,i3,1)*a(l)
                    s111=s111 + x(i1,1,i2+l,2,i3,1)*b(l)
                enddo
!y              y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
                ekin=ekin+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
                t121=0.d0
                do l=max(-i2,lowfil),min(lupfil,n2-i2)
                    t121=t121 + x(i1,1,i2+l,1,i3,1)*c(l)
                enddo
!y              y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
                ekin=ekin+t121*x(i1,1,i2,2,i3,1)
        enddo
    enddo
    enddo


       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       !write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:y',tel,1.d-6*mflop2/tel

! + (1/2) d^2/dz^2
    do i2=0,n2
    do i1=0,n1
        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
                t111=0.d0 ; s111=0.d0
                do l=max(-i3,lowfil),min(lupfil,n3-i3)
                    t111=t111 + x(i1,1,i2,1,i3+l,1)*a(l)
                    s111=s111 + x(i1,1,i2,1,i3+l,2)*b(l)
                enddo
!y              y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
                ekin=ekin+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
                t112=0.d0
                do l=max(-i3,lowfil),min(lupfil,n3-i3)
                    t112=t112 + x(i1,1,i2,1,i3+l,1)*c(l)
                enddo
!y              y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
                ekin=ekin+t112*x(i1,1,i2,1,i3,2)
        enddo
    enddo
    enddo

       call system_clock(ncount3,ncount_rate,ncount_max)
       tel=dble(ncount3-ncount2)/dble(ncount_rate)
       !write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:z',tel,1.d-6*mflop3/tel

! wavelet part
 ! (1/2) d^2/dx^2
    do i3=nfl3,nfu3
        do i2=nfl2,nfu2
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
                do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                    t112=t112 + x(i1+l,1,i2,1,i3,2)*a(l) + x(i1+l,2,i2,1,i3,2)*b(l)
                    t121=t121 + x(i1+l,1,i2,2,i3,1)*a(l) + x(i1+l,2,i2,2,i3,1)*b(l)
                    t122=t122 + x(i1+l,1,i2,2,i3,2)*a(l) + x(i1+l,2,i2,2,i3,2)*b(l)
                    t212=t212 + x(i1+l,1,i2,1,i3,2)*c(l) + x(i1+l,2,i2,1,i3,2)*e(l)
                    t221=t221 + x(i1+l,1,i2,2,i3,1)*c(l) + x(i1+l,2,i2,2,i3,1)*e(l)
                    t222=t222 + x(i1+l,1,i2,2,i3,2)*c(l) + x(i1+l,2,i2,2,i3,2)*e(l)
                    t211=t211 + x(i1+l,2,i2,1,i3,1)*e(l)
                enddo

!y              y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
!y              y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
!y              y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
!y              y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
!y              y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
!y              y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
!y              y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
                ekin=ekin+x(i1,1,i2,1,i3,2)*t112
                ekin=ekin+x(i1,1,i2,2,i3,1)*t121
                ekin=ekin+x(i1,2,i2,1,i3,1)*t211
                ekin=ekin+x(i1,1,i2,2,i3,2)*t122
                ekin=ekin+x(i1,2,i2,1,i3,2)*t212
                ekin=ekin+x(i1,2,i2,2,i3,1)*t221
                ekin=ekin+x(i1,2,i2,2,i3,2)*t222
            enddo
        enddo
    enddo

       call system_clock(ncount4,ncount_rate,ncount_max)
       tel=dble(ncount4-ncount3)/dble(ncount_rate)
       !write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:x',tel,1.d-6*nflop1/tel


 ! + (1/2) d^2/dy^2
    do i3=nfl3,nfu3
    do i1=nfl1,nfu1
       do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
                do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                    t112=t112 + x(i1,1,i2+l,1,i3,2)*a(l) + x(i1,1,i2+l,2,i3,2)*b(l)
                    t211=t211 + x(i1,2,i2+l,1,i3,1)*a(l) + x(i1,2,i2+l,2,i3,1)*b(l)
                    t122=t122 + x(i1,1,i2+l,1,i3,2)*c(l) + x(i1,1,i2+l,2,i3,2)*e(l)
                    t212=t212 + x(i1,2,i2+l,1,i3,2)*a(l) + x(i1,2,i2+l,2,i3,2)*b(l)
                    t221=t221 + x(i1,2,i2+l,1,i3,1)*c(l) + x(i1,2,i2+l,2,i3,1)*e(l)
                    t222=t222 + x(i1,2,i2+l,1,i3,2)*c(l) + x(i1,2,i2+l,2,i3,2)*e(l)
                    t121=t121 + x(i1,1,i2+l,2,i3,1)*e(l)
                enddo

!y              y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
!y              y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
!y              y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
!y              y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
!y              y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
!y              y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
!y              y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
                ekin=ekin+x(i1,1,i2,1,i3,2)*t112
                ekin=ekin+x(i1,1,i2,2,i3,1)*t121
                ekin=ekin+x(i1,2,i2,1,i3,1)*t211
                ekin=ekin+x(i1,1,i2,2,i3,2)*t122
                ekin=ekin+x(i1,2,i2,1,i3,2)*t212
                ekin=ekin+x(i1,2,i2,2,i3,1)*t221
                ekin=ekin+x(i1,2,i2,2,i3,2)*t222
       enddo
    enddo
    enddo

       call system_clock(ncount5,ncount_rate,ncount_max)
       tel=dble(ncount5-ncount4)/dble(ncount_rate)
       !write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:y',tel,1.d-6*nflop2/tel

 ! + (1/2) d^2/dz^2
    do i2=nfl2,nfu2
    do i1=nfl1,nfu1
       do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
                do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                    t121=t121 + x(i1,1,i2,2,i3+l,1)*a(l) + x(i1,1,i2,2,i3+l,2)*b(l)
                    t211=t211 + x(i1,2,i2,1,i3+l,1)*a(l) + x(i1,2,i2,1,i3+l,2)*b(l)
                    t122=t122 + x(i1,1,i2,2,i3+l,1)*c(l) + x(i1,1,i2,2,i3+l,2)*e(l)
                    t212=t212 + x(i1,2,i2,1,i3+l,1)*c(l) + x(i1,2,i2,1,i3+l,2)*e(l)
                    t221=t221 + x(i1,2,i2,2,i3+l,1)*a(l) + x(i1,2,i2,2,i3+l,2)*b(l)
                    t222=t222 + x(i1,2,i2,2,i3+l,1)*c(l) + x(i1,2,i2,2,i3+l,2)*e(l)
                    t112=t112 + x(i1,1,i2,1,i3+l,2)*e(l)
                enddo

!y              y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
!y              y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
!y              y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
!y              y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
!y              y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
!y              y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
!y              y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
                ekin=ekin+x(i1,1,i2,1,i3,2)*t112
                ekin=ekin+x(i1,1,i2,2,i3,1)*t121
                ekin=ekin+x(i1,2,i2,1,i3,1)*t211
                ekin=ekin+x(i1,1,i2,2,i3,2)*t122
                ekin=ekin+x(i1,2,i2,1,i3,2)*t212
                ekin=ekin+x(i1,2,i2,2,i3,1)*t221
                ekin=ekin+x(i1,2,i2,2,i3,2)*t222

       enddo
    enddo
    enddo

       call system_clock(ncount6,ncount_rate,ncount_max)
       tel=dble(ncount6-ncount5)/dble(ncount_rate)
      ! write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:z',tel,1.d-6*nflop3/tel

       tel=dble(ncount6-ncount0)/dble(ncount_rate)
      ! write(99,'(a40,2(1x,e10.3))') 'P:ALL   PART',  &
      !      tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

end subroutine ConvolkineticP


!> Calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
subroutine system_size(rxyz,rad, &
                   cxmin,cxmax,cymin,cymax,czmin,czmax)


   implicit real(kind=8) (a-h,o-z)
   parameter(eps_mach=1.d-12)
   dimension rxyz(3)



   cxmax=-1.d100 ; cxmin=1.d100
   cymax=-1.d100 ; cymin=1.d100
   czmax=-1.d100 ; czmin=1.d100

   cxmax=max(cxmax,rxyz(1)+rad) ; cxmin=min(cxmin,rxyz(1)-rad)
   cymax=max(cymax,rxyz(2)+rad) ; cymin=min(cymin,rxyz(2)-rad)
   czmax=max(czmax,rxyz(3)+rad) ; czmin=min(czmin,rxyz(3)-rad)

   cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
   cymax=cymax-eps_mach ; cymin=cymin+eps_mach
   czmax=czmax-eps_mach ; czmin=czmin+eps_mach

END SUBROUTINE system_size




logical function myorbital(iorb,norbe,iproc,nproc)
        implicit real(kind=8) (a-h,o-z)
        parameter(eps_mach=1.d-12)

        tt=dble(norbe)/dble(nproc)
        norbep=int((1.d0-eps_mach*tt) + tt)
        if (iorb .ge. iproc*norbep+1 .and. iorb .le. min((iproc+1)*norbep,norbe)) then
        myorbital=.true.
        else
        myorbital=.false.
        endif

        return
end function myorbital



!> Calculates the kinetic energy of an atomic wavefunction expressed in Gaussians
!! the output psiatn is a normalized version of psiat
subroutine atomkin(l,ng,xp,psiat,psiatn,ek)
        implicit real(kind=8) (a-h,o-z)
        dimension xp(ng),psiat(ng),psiatn(ng)
!        dimension xp(31),psiat(31),psiatn(31)
!!!write(*,*)'entered ATOMKIN, debug lines:'
!!!write(*,*)'l,ng,xp,psiat,psiatn,ek'
!!!write(*,*)l,ng,xp,psiat,psiatn,ek
!!!write(*,*)'-----------------------'
!        gml=.5d0*gamma(.5d0+l)
        gml = 0.d0
        if (l.eq.0) then 
            gml=0.88622692545275801365d0
        else if (l.eq.1) then 
            gml=0.44311346272637900682d0
        else if (l.eq.2) then 
            gml=0.66467019408956851024d0
        else if (l.eq.3) then 
            gml=1.6616754852239212756d0
        else
          stop 'atomkin'
        endif



        ek=0.d0
        tt=0.d0
        do i=1,ng
        xpi=.5d0/xp(i)**2
        do j=1,ng
!        do j=i,ng
        xpj=.5d0/xp(j)**2
        d=xpi+xpj
        sxp=1.d0/d
        const=gml*sqrt(sxp)**(2*l+1)
! kinetic energy  matrix element hij
        hij=.5d0*const*sxp**2* ( 3.d0*xpi*xpj +                  &
                     l*(6.d0*xpi*xpj-xpi**2-xpj**2) -        &
                     l**2*(xpi-xpj)**2  ) + .5d0*l*(l+1.d0)*const
        sij=const*sxp*(l+.5d0)
        ek=ek+hij*psiat(i)*psiat(j)
        tt=tt+sij*psiat(i)*psiat(j)
        enddo
        enddo

        if (abs(tt-1.d0).gt.1.d-2) write(*,*) 'presumably wrong inguess data',l,tt
! energy expectation value
        ek=ek/tt
!        write(*,*) 'ek=',ek,tt,l,ng
! scale atomic wavefunction
        tt=sqrt(1.d0/tt)
!!$        if (l.eq.0) then  ! multiply with 1/sqrt(4*pi)
!!$        tt=tt*0.28209479177387814347d0
!!$        else if (l.eq.1) then  ! multiply with sqrt(3/(4*pi))
!!$        tt=tt*0.48860251190291992159d0
!!$        !decide the value of the normalization to be used
!!$        endif
        do i=1,ng
        psiatn(i)=psiat(i)*tt
        enddo

        return
        END SUBROUTINE



subroutine calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
  
  implicit none
  integer :: l,m,nterm_max
  integer :: nterm
  integer, dimension(nterm_max) :: lx,ly,lz
  real(kind=8), dimension(nterm_max) :: fac_arr

  if (l.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.28209479177387814347d0

  else if (l.eq.2  .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.48860251190291992159d0
  else if (l.eq.2  .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=0.48860251190291992159d0
  else if (l.eq.2  .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=0.48860251190291992159d0

  else if (l.eq.3  .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.092548430592079d0
  else if (l.eq.3  .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.092548430592079d0
  else if (l.eq.3  .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.092548430592079d0
  else if (l.eq.3  .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.5462742152960396d0
     fac_arr(2)=-0.5462742152960396d0
  else if (l.eq.3  .and. m.eq.5) then !to be controlled, non normalized
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.3153915652525201d0
     fac_arr(2)=-0.3153915652525201d0
     fac_arr(3)=2.d0*0.3153915652525201d0 !!! this was wrong as 3d0*
!!$     nterm=2
!!$     lx(1)=0 ; ly(1)=0 ; lz(1)=0
!!$     lx(2)=0 ; ly(2)=0 ; lz(2)=2
!!$     fac_arr(1)=-0.3153915652525201d0
!!$     fac_arr(2)=3.d0*0.3153915652525201d0

  else if (l.eq.4  .and. m.eq.1) then
     nterm=2
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=0 ; lz(2)=2
     fac_arr(1)=-0.4570457994644658d0
     fac_arr(2)=5.d0*0.4570457994644658d0
  else if (l.eq.4  .and. m.eq.2) then
     nterm=2
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=1 ; lz(2)=2
     fac_arr(1)=-0.4570457994644658d0
     fac_arr(2)=5.d0*0.4570457994644658d0
  else if (l.eq.4  .and. m.eq.3) then
     nterm=2
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=0 ; lz(2)=3
     fac_arr(1)=-3.d0*0.3731763325901154d0
     fac_arr(2)=5.d0*0.3731763325901154d0
  else if (l.eq.4  .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.5900435899266436d0
     fac_arr(2)=-3.d0*0.5900435899266436d0
  else if (l.eq.4  .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-3.d0*0.5900435899266436d0
     fac_arr(2)=0.5900435899266436d0
  else if (l.eq.4  .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.445305721320277d0
     fac_arr(2)=-1.445305721320277d0
  else if (l.eq.4  .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.890611442640554d0
  else
     stop 'input guess format error'
  endif
  
END SUBROUTINE calc_coeff_inguess


        subroutine bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
        implicit real(kind=8) (a-h,o-z)
        logical logrid
        dimension logrid(0:n1,0:n2,0:n3)
        dimension ibyz(2,0:n2,0:n3),ibxz(2,0:n1,0:n3),ibxy(2,0:n1,0:n2)


        do 100, i3=0,n3
        do 100, i2=0,n2
        ibyz(1,i2,i3)= 1000
        ibyz(2,i2,i3)=-1000

        do i1=0,n1
         if (logrid(i1,i2,i3)) then
            ibyz(1,i2,i3)=i1
            goto 10
         endif
        enddo
10      continue
        do i1=n1,0,-1
         if (logrid(i1,i2,i3)) then
            ibyz(2,i2,i3)=i1
            goto 11
         endif
        enddo
11      continue

100     continue

        do 200,i3=0,n3
        do 200,i1=0,n1
        ibxz(1,i1,i3)= 1000
        ibxz(2,i1,i3)=-1000

        do i2=0,n2
         if (logrid(i1,i2,i3)) then
             ibxz(1,i1,i3)=i2
             goto 20
         endif
        enddo
20      continue
        do i2=n2,0,-1
         if (logrid(i1,i2,i3)) then
             ibxz(2,i1,i3)=i2
             goto 21
         endif
        enddo
21      continue

200     continue


        do 300, i2=0,n2
        do 300, i1=0,n1
        ibxy(1,i1,i2)= 1000
        ibxy(2,i1,i2)=-1000

        do i3=0,n3
         if (logrid(i1,i2,i3)) then
             ibxy(1,i1,i2)=i3
             goto 30
         endif
        enddo
30      continue
        do i3=n3,0,-1
         if (logrid(i1,i2,i3)) then
             ibxy(2,i1,i2)=i3
             goto 31
         endif
        enddo
31      continue

300     continue

        return
        END SUBROUTINE

       subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,  &
                               rxyz,rad,hgrid,logrid)
! set up an array logrid(i1,i2,i3) that specifies whether the grid point
! i1,i2,i3 is the center of a scaling function/wavelet
        implicit real(kind=8) (a-h,o-z)
        logical logrid
        parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
        dimension rxyz(3)!
        dimension logrid(0:n1,0:n2,0:n3)


        do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
         logrid(i1,i2,i3)=.false.
        enddo ; enddo ; enddo

!        rad=radii*rmult+nbuf*hgrid ! add/subtract nbuf below instead
        ml1=int(onem+(rxyz(1)-rad)/hgrid-nbuf)  ; mu1=int((rxyz(1)+rad)/hgrid+nbuf)
        ml2=int(onem+(rxyz(2)-rad)/hgrid-nbuf)  ; mu2=int((rxyz(2)+rad)/hgrid+nbuf)
        ml3=int(onem+(rxyz(3)-rad)/hgrid-nbuf)  ; mu3=int((rxyz(3)+rad)/hgrid+nbuf)
!NOTE buffer nbuf and lower bounds nl are zero, upper bounds nu are n1 n2 n3
        if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
        if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
        if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
        do i3=ml3,mu3
        dz2=(i3*hgrid-rxyz(3))**2
        do i2=ml2,mu2
        dy2=(i2*hgrid-rxyz(2))**2
        do i1=ml1,mu1
        dx=i1*hgrid-rxyz(1)
        if (dx**2+(dy2+dz2).lt.rad**2) then
              logrid(i1,i2,i3)=.true.
        endif
        enddo ; enddo ; enddo
!      enddo

        return
        END SUBROUTINE



!=================== gauss to daub ==================================
!!***
 !       PROGRAM MAIN
 !       implicit real(kind=8) (a-h,o-z)          
 !       PARAMETER(NMAX=1000,NWORK=10000)
 !       DIMENSION C(0:NMAX,2)
!!
!!       NOW THE DIMENSION OF C IS (0:NMAX,2) INSTEAD OF (-N_INTVX:N_INTVX)
!!
 !       DIMENSION WW(0:NWORK,2)  

 !           GAU_A=1.D0
 !           GAU_CEN=20.D0
 !           N_GAU=0
 !          FACTOR=1.D0
 !
 !        HGRID=1.D0

 !       CALL GAUSS_TO_DAUB(HGRID,FACTOR,GAU_CEN,GAU_A,N_GAU,&!NO ERR, ERRSUC
 !            NMAX,N_LEFT,N_RIGHT,C,ERR_NORM,&               !NO ERR_WAV. NMAX INSTEAD OF N_INTVX
 !            WW,NWORK)                                            !ADDED WORK ARRAY WW(:,:)  


 !           WRITE(*,*)'ERROR=',ERR_NORM

 !       END


         SUBROUTINE GAUSS_TO_DAUB(HGRID,FACTOR,GAU_CEN,GAU_A,N_GAU,&!NO ERR, ERRSUC
              NMAX,N_LEFT,N_RIGHT,C,ERR_NORM,&              !NO ERR_WAV. NMAX INSTEAD OF N_INTVX
              WW,NWORK)                             !ADDED WORK ARRAYS WW WITH DIMENSION  NWORK
! Gives the expansion coefficients of exp(-(1/2)*(x/gau_a)**2)
!!  INPUT: hgrid
!          FACTOR
!!         gau_cen
!!          gau_a
!!          n_gau
!           NMAX
!! OUTPUT: N_LEFT,N_RIGHT: Intervall where the GAUSSIAN IS LARGER THAN
!than thE MACHINE PRECISION
!!         C(:,1) array of scaling function coefficients:
!!         C(:,2) array of wavelet coefficients:
!!         WW(:,1),WW(:,2): work arrays that have to be 16 times larger than C
            implicit real(kind=8) (a-h,o-z)
            INTEGER LEFTS(0:4),RIGHTS(0:4),RIGHTX,LEFTX,RIGHT_T
            DIMENSION C(0:NMAX,2)
            DIMENSION WW(0:NWORK,2)
         INCLUDE 'recs16.inc'! MAGIC FILTER  
         INCLUDE 'intots.inc'! HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
         INCLUDE 'sym_16.inc'! WAVELET FILTERS

!
!            RESCALE THE PARAMETERS SO THAT HGRID GOES TO 1.D0  
!                   
             A=GAU_A/HGRID
             I0=NINT(GAU_CEN/HGRID) ! THE ARRAY IS CENTERED AT I0
             Z0=GAU_CEN/HGRID-I0
         
             H=.125D0*.5d0
!
!            CALCULATE THE ARRAY SIZES;
!            AT LEVEL 0, POSITIONS SHIFTED BY I0 
!
             RIGHT_T= CEILING(15.D0*A)

!             WRITE(*,*)'RIGHT_T=',RIGHT_T,'A=',A,'HGRID=',HGRID,'NMAX=',NMAX

!             IF (2*RIGHT_T.GT.NMAX) &
             IF (2*RIGHT_T.GT.NWORK) &
!             STOP 'A GAUSSIAN IS GREATER THAN THE CELL'
             STOP 'INCREASE THE NWORK IN SUBROUTINE gauss_to_daub.f90'

             LEFTS( 0)=MAX(I0-RIGHT_T,   0)
             RIGHTS(0)=MIN(I0+RIGHT_T,NMAX)

             N_LEFT=LEFTS(0)
             N_RIGHT=RIGHTS(0)

             DO K=1,4
               RIGHTS(K)=2*RIGHTS(K-1)+M
               LEFTS( K)=2*LEFTS( K-1)-M
             ENDDO 

             LEFTX = LEFTS(4)-N
             RIGHTX=RIGHTS(4)+N  
!
!            EIGENTLICH, CALCULATE THE EXPANSION COEFFICIENTS
!            AT LEVEL 4, POSITIONS SHIFTED BY 16*I0 
!         
             DO I=LEFTX,RIGHTX
               WW(I-LEFTX,1)=PSI_G((I-I0*16)*H,A,Z0,N_GAU)
             ENDDO 

             CALL APPLY_W(WW(:,1),WW(:,2),&
                                LEFTX   ,RIGHTX   ,LEFTS(4),RIGHTS(4),H)

             CALL FORWARD_C(WW(:,2),WW(:,1),&
                                LEFTS(4),RIGHTS(4),LEFTS(3),RIGHTS(3)) 
             CALL FORWARD_C(WW(:,1),WW(:,2),&
                                LEFTS(3),RIGHTS(3),LEFTS(2),RIGHTS(2)) 
             CALL FORWARD_C(WW(:,2),WW(:,1),&
                                LEFTS(2),RIGHTS(2),LEFTS(1),RIGHTS(1)) 

             CALL FORWARD(  WW(:,1),WW(:,2),&
                                LEFTS(1),RIGHTS(1),LEFTS(0),RIGHTS(0)) 

             C=0.D0
             LENGTH=N_RIGHT-N_LEFT+1
             DO I=0,LENGTH-1
               C(I+N_LEFT,1)=WW(I       ,2) !N_LEFT..N_RIGHT <->    0  ..  LENGTH-1
               C(I+N_LEFT,2)=WW(I+LENGTH,2) !N_LEFT..N_RIGHT <-> LENGTH..2*LENGTH-1
             ENDDO 
!
!            CALCULATE THE (RELATIVE) ERROR
!
             CN2=0.D0
             DO I=0,LENGTH*2-1
               CN2=CN2+WW(I,2)**2
             ENDDO 
             
             THEOR_NORM2=VALINTS(N_GAU)*A**(2*N_GAU+1)

             ERROR=SQRT(ABS(1-CN2/THEOR_NORM2))
!
!            RESCALE BACK THE COEFFICIENTS AND THE ERROR
!
             FAC= HGRID**N_GAU*SQRT(HGRID)*FACTOR
             C=C*FAC
             ERR_NORM=ERROR*FAC
!
!            CALCULATE THE OUTPUT ARRAY DIMENSIONS
!

!             WRITE(*,*)'N_LEFT=',N_LEFT,'        N_RIGHT=',N_RIGHT

         END SUBROUTINE GAUSS_TO_DAUB
         
        SUBROUTINE APPLY_W(CX,C,LEFTX,RIGHTX,LEFT,RIGHT,H)
!
!       APPLYING THE MAGIC FILTER ("SHRINK") 
!
        implicit real(kind=8) (a-h,o-z)
        INTEGER RIGHT,RIGHTX
        DIMENSION CX(LEFTX:RIGHTX),C(LEFT:RIGHT)
        INCLUDE 'recs16.inc'

        SQH=SQRT(H)
        
        DO I=LEFT,RIGHT
          CI=0.D0
          DO J=-N,N
            CI=CI+CX(I+J)*W(J)         
          ENDDO
          C(I)=CI*SQH
        ENDDO
        
        END SUBROUTINE APPLY_W


      SUBROUTINE FORWARD_C(C,C_1,LEFT,RIGHT,LEFT_1,RIGHT_1)
!
!      FORWARD WAVELET TRANSFORM WITHOUT WAVELETS ("SHRINK")
!
       implicit real(kind=8) (a-h,o-z)
       INTEGER RIGHT,RIGHT_1
       DIMENSION C(LEFT:RIGHT)
       DIMENSION C_1(LEFT_1:RIGHT_1)

       INCLUDE 'sym_16.inc'

!
!      GET THE COARSE SCFUNCTIONS AND WAVELETS
!
       DO I=LEFT_1,RIGHT_1
         I2=2*I
         CI=0.D0
         DO J=-M,M
           CI=CI+CHT(J)*C(J+I2)
         ENDDO
         C_1(I)=CI
       ENDDO

       END SUBROUTINE FORWARD_C

      SUBROUTINE FORWARD(C,CD_1,LEFT,RIGHT,LEFT_1,RIGHT_1)
!
!      CONVENTIONAL FORWARD WAVELET TRANSFORM ("SHRINK")
!
       implicit real(kind=8) (a-h,o-z)
       INTEGER RIGHT,RIGHT_1
       DIMENSION C(LEFT:RIGHT)
       DIMENSION CD_1(LEFT_1:RIGHT_1,2)

       INCLUDE 'sym_16.inc'

!
!      GET THE COARSE SCFUNCTIONS AND WAVELETS
!
       DO I=LEFT_1,RIGHT_1
         I2=2*I
         CI=0.D0
         DI=0.D0
         DO J=-M,M
           CI=CI+CHT(J)*C(J+I2)
           DI=DI+CGT(J)*C(J+I2)
         ENDDO
         CD_1(I,1)=CI
         CD_1(I,2)=DI
       ENDDO
 
       END SUBROUTINE FORWARD

       function psi_g(x,GAU_A,GAU_CEN,N_GAU)
       implicit real(kind=8) (a-h,o-z)
         psi_g=(X-GAU_CEN)**N_GAU*exp(-0.5d0*((X-GAU_CEN)/GAU_A)**2)
       end function psi_g

