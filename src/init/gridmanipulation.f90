!!****f* BigDFT/system_size
!! FUNCTION
!!   Calculates the overall size of the simulation cell 
!!   and shifts the atoms such that their position is the most symmetric possible.
!!   Assign these values to the global localisation region descriptor.
!!
!! COPYRIGHT
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: iproc
  real(gp), intent(in) :: crmult,frmult
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(gp), intent(inout) :: hx,hy,hz
  type(locreg_descriptors), intent(out) :: Glr
  real(gp), dimension(3), intent(out) :: shift
  !local variables
  integer, parameter :: lupfil=14
  real(gp), parameter ::eps_mach=1.e-12_gp
  integer :: iat,j,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  real(gp) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3

  !check the geometry code with the grid spacings
  if (atoms%geocode == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
     write(*,'(1x,a)')'ERROR: The values of the grid spacings must be equal in the Free BC case'
     stop
  end if

  !calculate the extremes of the boxes taking into account the spheres around the atoms
  cxmax=-1.e10_gp 
  cxmin=1.e10_gp

  cymax=-1.e10_gp 
  cymin=1.e10_gp

  czmax=-1.e10_gp 
  czmin=1.e10_gp

  do iat=1,atoms%nat

     rad=radii_cf(atoms%iatype(iat),1)*crmult

     cxmax=max(cxmax,rxyz(1,iat)+rad) 
     cxmin=min(cxmin,rxyz(1,iat)-rad)

     cymax=max(cymax,rxyz(2,iat)+rad) 
     cymin=min(cymin,rxyz(2,iat)-rad)
     
     czmax=max(czmax,rxyz(3,iat)+rad) 
     czmin=min(czmin,rxyz(3,iat)-rad)
  enddo

!eliminate epsilon form the grid size calculation
!!  cxmax=cxmax+eps_mach 
!!  cymax=cymax+eps_mach  
!!  czmax=czmax+eps_mach  
!!
!!  cxmin=cxmin-eps_mach
!!  cymin=cymin-eps_mach
!!  czmin=czmin-eps_mach


  !define the box sizes for free BC, and calculate dimensions for the fine grid with ISF
  if (atoms%geocode == 'F') then
     atoms%alat1=(cxmax-cxmin)
     atoms%alat2=(cymax-cymin)
     atoms%alat3=(czmax-czmin)

     ! grid sizes n1,n2,n3
     n1=int(atoms%alat1/hx)
!if (mod(n1,2)==1) n1=n1+1
     n2=int(atoms%alat2/hy)
!if (mod(n2,2)==1) n2=n2+1
     n3=int(atoms%alat3/hz)
!if (mod(n3,2)==1) n3=n3+1
     alatrue1=real(n1,gp)*hx
     alatrue2=real(n2,gp)*hy
     alatrue3=real(n3,gp)*hz

     n1i=2*n1+31
     n2i=2*n2+31
     n3i=2*n3+31

  else if (atoms%geocode == 'P') then 
     !define the grid spacings, controlling the FFT compatibility
     call correct_grid(atoms%alat1,hx,n1)
     call correct_grid(atoms%alat2,hy,n2)
     call correct_grid(atoms%alat3,hz,n3)
     alatrue1=(cxmax-cxmin)
     alatrue2=(cymax-cymin)
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+2
     n3i=2*n3+2

  else if (atoms%geocode == 'S') then
     call correct_grid(atoms%alat1,hx,n1)
     atoms%alat2=(cymax-cymin)
     call correct_grid(atoms%alat3,hz,n3)

     alatrue1=(cxmax-cxmin)
     n2=int(atoms%alat2/hy)
     alatrue2=real(n2,gp)*hy
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+31
     n3i=2*n3+2

  end if

  !balanced shift taking into account the missing space
  cxmin=cxmin+0.5_gp*(atoms%alat1-alatrue1)
  cymin=cymin+0.5_gp*(atoms%alat2-alatrue2)
  czmin=czmin+0.5_gp*(atoms%alat3-alatrue3)

  !correct the box sizes for the isolated case
  if (atoms%geocode == 'F') then
     atoms%alat1=alatrue1
     atoms%alat2=alatrue2
     atoms%alat3=alatrue3
  else if (atoms%geocode == 'S') then
     cxmin=0.0_gp
     atoms%alat2=alatrue2
     czmin=0.0_gp
  else if (atoms%geocode == 'P') then
     !for the moment we do not put the shift, at the end it will be tested
     !here we should put the center of mass
     cxmin=0.0_gp
     cymin=0.0_gp
     czmin=0.0_gp
  end if

  !assign the shift to the atomic positions
  shift(1)=cxmin
  shift(2)=cymin
  shift(3)=czmin

  !here we can put a modulo operation for periodic directions
  do iat=1,atoms%nat
     rxyz(1,iat)=rxyz(1,iat)-shift(1)
     rxyz(2,iat)=rxyz(2,iat)-shift(2)
     rxyz(3,iat)=rxyz(3,iat)-shift(3)
  enddo

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 
  nfl2=n2 
  nfl3=n3

  nfu1=0 
  nfu2=0 
  nfu3=0

  do iat=1,atoms%nat
     rad=radii_cf(atoms%iatype(iat),2)*frmult
     if (rad > 0.0_gp) then
        nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hx - eps_mach))
        nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hx + eps_mach))
        
        nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hy - eps_mach))
        nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hy + eps_mach))
        
        nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hz - eps_mach)) 
        nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hz + eps_mach))
     end if
  enddo

  !correct the values of the delimiter if they go outside the box
  if (nfl1 < 0 .or. nfu1 > n1) then
     nfl1=0
     nfu1=n1
  end if
  if (nfl2 < 0 .or. nfu2 > n2) then
     nfl2=0
     nfu2=n2
  end if
  if (nfl3 < 0 .or. nfu3 > n3) then
     nfl3=0
     nfu3=n3
  end if

  !correct the values of the delimiter if there are no wavelets
  if (nfl1 == n1 .and. nfu1 == 0) then
     nfl1=n1/2
     nfu1=n1/2
  end if
  if (nfl2 == n2 .and. nfu2 == 0) then
     nfl2=n2/2
     nfu2=n2/2
  end if
  if (nfl3 == n3 .and. nfu3 == 0) then
     nfl3=n3/2
     nfu3=n3/2
  end if


  if (iproc == 0) then
     write(*,'(1x,a,19x,a)') 'Shifted atomic positions, Atomic Units:','grid spacing units:'
     do iat=1,atoms%nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5),3x,3(1x,0pf9.3))') &
             iat,trim(atoms%atomnames(atoms%iatype(iat))),&
             (rxyz(j,iat),j=1,3),rxyz(1,iat)/hx,rxyz(2,iat)/hy,rxyz(3,iat)/hz
     enddo
     write(*,'(1x,a,3(1x,1pe12.5),a,3(1x,0pf7.4))') &
          '   Shift of=',-cxmin,-cymin,-czmin,' H grids=',hx,hy,hz
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '  Box Sizes=',atoms%alat1,atoms%alat2,atoms%alat3,n1,n2,n3
     write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
          '      Extremes for the high resolution grid points:',&
          nfl1,'<',nfu1,nfl2,'<',nfu2,nfl3,'<',nfu3
  endif

  !assign the values
  Glr%d%n1  =n1  
  Glr%d%n2  =n2  
  Glr%d%n3  =n3  
  Glr%d%n1i =n1i 
  Glr%d%n2i =n2i 
  Glr%d%n3i =n3i 
  Glr%d%nfl1=nfl1
  Glr%d%nfl2=nfl2
  Glr%d%nfl3=nfl3
  Glr%d%nfu1=nfu1
  Glr%d%nfu2=nfu2
  Glr%d%nfu3=nfu3

  Glr%ns1=0
  Glr%ns2=0
  Glr%ns3=0

  !while using k-points this condition should be disabled
  !evaluate if the conditiond for the hybrid evaluation if periodic BC hold
  Glr%hybrid_on=                   (nfu1-nfl1+lupfil < n1+1)
  Glr%hybrid_on=(Glr%hybrid_on.and.(nfu2-nfl2+lupfil < n2+1))
  Glr%hybrid_on=(Glr%hybrid_on.and.(nfu3-nfl3+lupfil < n3+1))

  if (Glr%hybrid_on) then
     if (iproc == 0) write(*,*)'wavelet localization is ON'
  else
     if (iproc == 0) write(*,*)'wavelet localization is OFF'
  endif

END SUBROUTINE system_size
!!***


!!****f* BigDFT/correct_grid
!! FUNCTION
!!   Here the dimensions should be corrected in order to 
!!   allow the fft for the preconditioner and for Poisson Solver
!! SOURCE
!!
subroutine correct_grid(a,h,n)
  use module_base
  use Poisson_Solver
  implicit none
  real(gp), intent(in) :: a
  integer, intent(inout) :: n
  real(gp), intent(inout) :: h
  !local variables
  integer :: m,m2,nt

  n=ceiling(a/h)-1
  nt=n+1
  do
     !correct the direct dimension
     call fourier_dim(nt,m)

     !control if the double of this dimension is compatible with the FFT
     call fourier_dim(2*m,m2)
     !if this check is passed both the preconditioner and the PSolver works
     if (m2==2*m .and. mod(m,2) ==0) exit !only even dimensions are considered so far

     nt=m+1
  end do
  n=m-1

!!!  !here the dimensions should be corrected in order to 
!!!  !allow the fft for the preconditioner
!!!  m=2*n+2
!!!  do 
!!!     call fourier_dim(m,m)
!!!     if ((m/2)*2==m) then
!!!        n=(m-2)/2
!!!        exit
!!!     else
!!!        m=m+1
!!!     end if
!!!  end do

  h=a/real(n+1,gp)
  
END SUBROUTINE correct_grid
!!***


!!****f* BigDFT/num_segkeys
!! FUNCTION
!!   Calculates the length of the keys describing a wavefunction data structure
!! SOURCE
!!
subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid 
  integer, intent(out) :: mseg,mvctr
  !local variables
  logical :: plogrid
  integer :: i1,i2,i3,nsrt,nend,nsrti,nendi,mvctri
  mvctr=0
  nsrt=0
  nend=0
!$omp parallel default(private) shared(nl3,nu3,nl2,nu2,nl1,nu1,logrid,mvctr,nsrt,nend)
  mvctri=0
  nsrti=0
  nendi=0
!$omp do  
  do i3=nl3,nu3 
     do i2=nl2,nu2
        plogrid=.false.
        do i1=nl1,nu1
           if (logrid(i1,i2,i3)) then
              mvctri=mvctri+1
              if (.not. plogrid) then
                 nsrti=nsrti+1
              endif
           else
              if (plogrid) then
                 nendi=nendi+1
              endif
           endif
           plogrid=logrid(i1,i2,i3)
        enddo
        if (plogrid) then
           nendi=nendi+1
        endif
     enddo
  enddo
!$omp enddo
!$omp critical
mvctr=mvctr+mvctri
nsrt=nsrt+nsrti
nend=nend+nendi
!$omp end critical
!$omp end parallel
  if (nend /= nsrt) then 
     write(*,*)' ERROR: nend <> nsrt',nend,nsrt
     stop 
  endif
  mseg=nend
  
END SUBROUTINE num_segkeys
!!***


!!****f* BigDFT/segkeys
!! FUNCTION
!!   Calculates the keys describing a wavefunction data structure
!! SOURCE
!!
subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
  !implicit real(kind=8) (a-h,o-z)
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,mseg
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid  
  integer, dimension(mseg), intent(out) :: keyv
  integer, dimension(2,mseg), intent(out) :: keyg
  !local variables
  logical :: plogrid
  integer :: mvctr,nsrt,nend,i1,i2,i3,ngridp

  mvctr=0
  nsrt=0
  nend=0
  do i3=nl3,nu3 
     do i2=nl2,nu2
        plogrid=.false.
        do i1=nl1,nu1
           ngridp=i3*((n1+1)*(n2+1)) + i2*(n1+1) + i1+1
           if (logrid(i1,i2,i3)) then
              mvctr=mvctr+1
              if (.not. plogrid) then
                 nsrt=nsrt+1
                 keyg(1,nsrt)=ngridp
                 keyv(nsrt)=mvctr
              endif
           else
              if (plogrid) then
                 nend=nend+1
                 keyg(2,nend)=ngridp-1
              endif
           endif
           plogrid=logrid(i1,i2,i3)
        enddo
        if (plogrid) then
           nend=nend+1
           keyg(2,nend)=ngridp
        endif
     enddo
  enddo
  if (nend /= nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif
  !mseg=nend
END SUBROUTINE segkeys
!!***


!!****f* BigDFT/fill_logrid
!! FUNCTION
!!   set up an array logrid(i1,i2,i3) that specifies whether the grid point
!!   i1,i2,i3 is the center of a scaling function/wavelet
!!
!! SOURCE
!!
subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
  real(gp), intent(in) :: rmult,hx,hy,hz
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(ntypes), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
  real(gp) :: dx,dy2,dz2,rad

  !some checks
  if (geocode /='F') then
     !the nbuf value makes sense only in the case of free BC
     if (nbuf /=0) then
        write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
        stop
     end if
     !the grid spacings must be the same
     if (hx/= hy .or. hy /=hz .or. hx/=hz) then
!        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
     end if
  end if

  if (geocode == 'F') then
     do i3=nl3,nu3 
        do i2=nl2,nu2 
           do i1=nl1,nu1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  else !
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  end if

  do iat=1,nat
     rad=radii(iatype(iat))*rmult+real(nbuf,gp)*hx
     if (rad /= 0.0_gp) then
        ml1=ceiling((rxyz(1,iat)-rad)/hx - eps_mach)  
        ml2=ceiling((rxyz(2,iat)-rad)/hy - eps_mach)   
        ml3=ceiling((rxyz(3,iat)-rad)/hz - eps_mach)   
        mu1=floor((rxyz(1,iat)+rad)/hx + eps_mach)
        mu2=floor((rxyz(2,iat)+rad)/hy + eps_mach)
        mu3=floor((rxyz(3,iat)+rad)/hz + eps_mach)
        !for Free BC, there must be no incoherences with the previously calculated delimiters
        if (geocode == 'F') then
           if (ml1 < nl1) stop 'ml1 < nl1'
           if (ml2 < nl2) stop 'ml2 < nl2'
           if (ml3 < nl3) stop 'ml3 < nl3'

           if (mu1 > nu1) stop 'mu1 > nu1'
           if (mu2 > nu2) stop 'mu2 > nu2'
           if (mu3 > nu3) stop 'mu3 > nu3'
        end if
        !what follows works always provided the check before
!$omp parallel default(shared) private(i3,dz2,j3,i2,dy2,j2,i1,j1,dx)
!$omp do
        do i3=max(ml3,-n3/2-1),min(mu3,n3+n3/2+1)
           dz2=(real(i3,gp)*hz-rxyz(3,iat))**2
           j3=modulo(i3,n3+1)
           do i2=max(ml2,-n2/2-1),min(mu2,n2+n2/2+1)
              dy2=(real(i2,gp)*hy-rxyz(2,iat))**2
              j2=modulo(i2,n2+1)
              do i1=max(ml1,-n1/2-1),min(mu1,n1+n1/2+1)
                 j1=modulo(i1,n1+1)
                 dx=real(i1,gp)*hx-rxyz(1,iat)
                 if (dx**2+(dy2+dz2) <= rad**2) then 
                    logrid(j1,j2,j3)=.true.
                 endif
              enddo
           enddo
        enddo
!$omp enddo
!$omp end parallel
     end if
  enddo

END SUBROUTINE fill_logrid
!!***


!!****f* BigDFT/make_bounds
!! FUNCTION
!!
!! SOURCE
!!
subroutine make_bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
  implicit none
  integer, intent(in) :: n1,n2,n3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid
  integer, dimension(2,0:n2,0:n3), intent(out) :: ibyz
  integer, dimension(2,0:n1,0:n3), intent(out) :: ibxz
  integer, dimension(2,0:n1,0:n2), intent(out) :: ibxy
  !local variables
  integer :: i1,i2,i3

  do i3=0,n3 
     do i2=0,n2 
        ibyz(1,i2,i3)= 1000
        ibyz(2,i2,i3)=-1000

        loop_i1s: do i1=0,n1
           if (logrid(i1,i2,i3)) then 
              ibyz(1,i2,i3)=i1
              exit loop_i1s
           endif
        enddo loop_i1s

        loop_i1e: do i1=n1,0,-1
           if (logrid(i1,i2,i3)) then 
              ibyz(2,i2,i3)=i1
              exit loop_i1e
           endif
        enddo loop_i1e
     end do
  end do


  do i3=0,n3 
     do i1=0,n1
        ibxz(1,i1,i3)= 1000
        ibxz(2,i1,i3)=-1000

        loop_i2s: do i2=0,n2 
           if (logrid(i1,i2,i3)) then 
              ibxz(1,i1,i3)=i2
              exit loop_i2s
           endif
        enddo loop_i2s

        loop_i2e: do i2=n2,0,-1
           if (logrid(i1,i2,i3)) then 
              ibxz(2,i1,i3)=i2
              exit loop_i2e
           endif
        enddo loop_i2e

     end do
  end do

  do i2=0,n2 
     do i1=0,n1 
        ibxy(1,i1,i2)= 1000
        ibxy(2,i1,i2)=-1000

        loop_i3s: do i3=0,n3
           if (logrid(i1,i2,i3)) then 
              ibxy(1,i1,i2)=i3
              exit loop_i3s
           endif
        enddo loop_i3s

        loop_i3e: do i3=n3,0,-1
           if (logrid(i1,i2,i3)) then 
              ibxy(2,i1,i2)=i3
              exit loop_i3e
           endif
        enddo loop_i3e
     end do
  end do

END SUBROUTINE make_bounds
!!***






subroutine initializeLocRegLIN(iproc, nproc, lr, orbsLIN, at, input, rxyz, radii_cf)
!
! Purpose:
! ========
!   Determines the localization regions for each orbital.
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces
implicit none 

! Calling arguments
integer:: iproc, nproc
type(locreg_descriptors):: lr
type(orbitals_data),intent(in):: orbsLIN
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat):: rxyz
real(8),dimension(at%ntypes,3):: radii_cf
type(communications_arrays):: commsLIN

! Local variables
integer:: iorb, iiAt, iitype, istat, ierr, ii, ngroups, norbPerGroup, jprocStart, jprocEnd, norbtot1
integer:: norbtot2, igroup, jproc, norbpMax, lproc, uproc, wholeGroup
real(8):: radius, radiusCut, idealSplit
logical,dimension(:,:,:,:),allocatable:: logridCut_c
logical,dimension(:,:,:,:),allocatable:: logridCut_f
character(len=*),parameter:: subname='initializeLocRegLIN'
integer,dimension(:),allocatable:: norbPerGroupArr, newID, newGroup, newComm
integer,dimension(:,:),allocatable:: procsInGroup
integer,dimension(:,:,:),allocatable:: tempArr
real(8),dimension(:),pointer:: psi, psiWork


! First check wheter we have free boundary conditions.
if(lr%geocode/='F' .and. lr%geocode/='f') then
    if(iproc==0) write(*,'(a)') 'initializeLocRegLIN only implemented for free boundary conditions!'
    call mpi_barrier(mpi_comm_world, ierr)
    stop
end if

! logridCut is a logical array that is trie for a given grid point if this point is within the
! localization radius and false otherwise. It, of course, different for each orbital.
! Maybe later this array can be changed such that it does not cover the whole simulation box,
! but only a part of it.
allocate(logridCut_c(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3,orbsLIN%norbp), stat=istat)
allocate(logridCut_f(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3,orbsLIN%norbp), stat=istat)

norbpMax=maxval(orbsLIN%norb_par)
allocate(lr%wfdLIN(norbpMax,0:nproc-1), stat=istat)
do jproc=0,nproc-1
    do iorb=1,norbpMax
        lr%wfdLIN(iorb,jproc)%nseg_c=0
        lr%wfdLIN(iorb,jproc)%nseg_f=0
        lr%wfdLIN(iorb,jproc)%nvctr_c=0
        lr%wfdLIN(iorb,jproc)%nvctr_f=0
    end do
end do
write(*,*) 'iproc, size(lr%wfdLIN)', iproc, size(lr%wfdLIN)

! Fill logridCut for each orbital
radiusCut=10.d0
do iorb=1,orbsLIN%norbp

    iiAt=orbsLIN%onWhichAtom(iorb)
    iitype=at%iatype(iiAt)
    radius=radii_cf(1,iitype)
    call fill_logridCUT(lr%geocode, lr%d%n1, lr%d%n2, lr%d%n3, 0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, 0, 1, &
         1, 1, rxyz(1,iiAt), radius, radiusCut, input%hx, input%hy, input%hz, logridCut_c(0,0,0,iorb))
    !!subroutine fill_logridCut(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
    !!     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)

    ! Calculate the number of segments and the number of grid points for each orbital.
    call num_segkeys(lr%d%n1, lr%d%n2, lr%d%n3, 0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, logridCut_c(0,0,0,iorb), &
        lr%wfdLIN(iorb,iproc)%nseg_c, lr%wfdLIN(iorb,iproc)%nvctr_c)
    !!subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
    !write(*,'(a,2i4,2x,2i4,2i8)') 'iproc, iorb, iiAt, iitype, lr%wfdLIN(iorb)%nseg_c, lr%wfdLIN(iorb)%nvctr_c', iproc, iorb, iiAt, iitype, lr%wfdLIN(iorb)%nseg_c, lr%wfdLIN(iorb)%nvctr_c

    ! Now the same procedure for the fine radius.
    radius=radii_cf(2,iitype)
    call fill_logridCUT(lr%geocode, lr%d%n1, lr%d%n2, lr%d%n3, 0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, 0, 1, &
         1, 1, rxyz(1,iiAt), radius, radiusCut, input%hx, input%hy, input%hz, logridCut_f(0,0,0,iorb))
    !!subroutine fill_logridCut(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
    !!     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)

    ! Calculate the number of segments and the number of grid points for each orbital.
    call num_segkeys(lr%d%n1, lr%d%n2, lr%d%n3, 0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, logridCut_f(0,0,0,iorb), &
        lr%wfdLIN(iorb,iproc)%nseg_f, lr%wfdLIN(iorb,iproc)%nvctr_f)
    !!subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
    !write(*,'(a,2i4,2x,2i4,2i8)') 'iproc, iorb, iiAt, iitype, lr%wfdLIN(iorb)%nseg_f, lr%wfdLIN(iorb)%nvctr_f', iproc, iorb, iiAt, iitype, lr%wfdLIN(iorb)%nseg_f, lr%wfdLIN(iorb)%nvctr_f



    ! Now fill the descriptors.
    call allocate_wfd(lr%wfdLIN(iorb,iproc),subname)
    ! First the coarse part
    call segkeys(lr%d%n1, lr%d%n2, lr%d%n3, 0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, logridCut_c(0,0,0,iorb), &
        lr%wfdLIN(iorb,iproc)%nseg_c, lr%wfdLIN(iorb,iproc)%keyg(1,1), lr%wfdLIN(iorb,iproc)%keyv(1))
    ! And then the fine part
    ii=lr%wfdLIN(iorb,iproc)%nseg_c+1
    call segkeys(lr%d%n1, lr%d%n2, lr%d%n3, 0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, logridCut_f(0,0,0,iorb), &
        lr%wfdLIN(iorb,iproc)%nseg_f, lr%wfdLIN(iorb,iproc)%keyg(1,ii), lr%wfdLIN(iorb,iproc)%keyv(ii))
    !!subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
end do


allocate(tempArr(1:norbpMax,0:nproc-1,2))
tempArr=0
do jproc=0,nproc-1
    do iorb=1,orbsLIN%norb_par(jproc)
        tempArr(iorb,jproc,2)=lr%wfdLIN(iorb,jproc)%nvctr_c
    end do
end do
call mpi_allreduce(tempArr(1,0,2), tempArr(1,0,1), norbpMax*nproc, mpi_integer, mpi_sum, mpi_comm_world, ierr)
do jproc=0,nproc-1
    do iorb=1,orbsLIN%norb_par(jproc)
        lr%wfdLIN(iorb,jproc)%nvctr_c=tempArr(iorb,jproc,1)
        tempArr(iorb,jproc,1)=0
    end do
end do

do jproc=0,nproc-1
    do iorb=1,orbsLIN%norb_par(jproc)
        tempArr(iorb,jproc,2)=lr%wfdLIN(iorb,jproc)%nvctr_f
    end do
end do
call mpi_allreduce(tempArr(1,0,2), tempArr(1,0,1), norbpMax*nproc, mpi_integer, mpi_sum, mpi_comm_world, ierr)
do jproc=0,nproc-1
    do iorb=1,orbsLIN%norb_par(jproc)
        lr%wfdLIN(iorb,jproc)%nvctr_f=tempArr(iorb,jproc,1)
    end do
end do




! Now divide the system in parts for the orthonormalization.
! At the moment just make ngroups groups.
norbPerGroup=30
ngroups=nint(dble(orbsLIN%norb)/dble(norbPerGroup))
if(ngroups>nproc) then
    if(iproc==0) write(*,'(a,i0,a,i0,a)') 'WARNING: change ngroups from ', ngroups, ' to ', nproc,'!'
    ngroups=nproc
end if
ngroups=min(ngroups,nproc)
idealSplit=dble(orbsLIN%norb)/dble(ngroups)
allocate(norbPerGroupArr(ngroups), stat=istat)
allocate(procsInGroup(2,ngroups))
! Now distribute the orbitals to the groups. Do not split MPI processes, i.e.
! all orbitals for one process will remain with this proccess.
norbtot1=0
norbtot2=0
igroup=1
jprocStart=0
jprocEnd=0
do jproc=0,nproc-1
    if(igroup==ngroups) then
        ! This last group has to take all the rest
        do ii=jproc,nproc-1
            norbtot1=norbtot1+orbsLIN%norb_par(jproc)
            jprocEnd=jprocEnd+1
        end do
    else
        norbtot1=norbtot1+orbsLIN%norb_par(jproc)
        if(jproc<nproc-1) norbtot2=norbtot1+orbsLIN%norb_par(jproc+1)
        jprocEnd=jprocEnd+1
    end if
    if(abs(dble(norbtot1)-idealSplit)<abs(dble(norbtot2-idealSplit)) .or. igroup==ngroups) then
        ! Here is the split between two groups
        norbPerGroupArr(igroup)=norbtot1
        procsInGroup(1,igroup)=jprocStart
        procsInGroup(2,igroup)=jprocEnd-1
        norbtot1=0
        norbtot2=0
        jprocStart=jproc+1
        if(igroup==ngroups) exit
        igroup=igroup+1
    end if
end do

do igroup=1,ngroups
    if(iproc==0) write(*,'(a,3i12)') 'norbPerGroupArr(igroup), procsInGroup(1,igroup), procsInGroup(2,igroup)',&
        norbPerGroupArr(igroup), procsInGroup(1,igroup), procsInGroup(2,igroup)
end do



! Create the new communicators.
allocate(newID(0:nproc), stat=istat)
allocate(newGroup(1:ngroups))
allocate(newComm(1:ngroups))
do jproc=0,nproc-1
    newID(jproc)=jproc
end do
call mpi_comm_group(mpi_comm_world, wholeGroup, ierr)
do igroup=1,ngroups
    call mpi_group_incl(wholeGroup, newID(procsInGroup(2,igroup))-newID(procsInGroup(1,igroup))+1,&
        newID(procsInGroup(1,igroup)), newGroup(igroup), ierr)
    call mpi_comm_create(mpi_comm_world, newGroup(igroup), newComm(igroup), ierr)
end do




! Now create the parameters for the transposition.
do igroup=1,ngroups
    lproc=procsInGroup(1,igroup)
    uproc=procsInGroup(2,igroup)
    if(iproc>=lproc .and. iproc<=uproc) then
        call orbitals_communicatorsLIN_group(iproc, lproc, uproc, lr, orbsLIN, commsLIN, newComm(igroup), norbPerGroupArr(igroup))
    end if
end do


do igroup=1,ngroups
    lproc=procsInGroup(1,igroup)
    uproc=procsInGroup(2,igroup)
    do jproc=lproc,uproc
        if(iproc>=lproc .and. iproc<=uproc) write(*,'(a,3i5,4i12)') 'iproc, igroup, jproc, commsLIN%ncntdLIN(jproc), &
            & commsLIN%ndspldLIN(jproc), commsLIN%ncnttLIN(jproc), commsLIN%ndspltLIN(jproc)', &
            iproc, igroup, jproc, commsLIN%ncntdLIN(jproc), commsLIN%ndspldLIN(jproc), &
            commsLIN%ncnttLIN(jproc), commsLIN%ndspltLIN(jproc)
    end do
end do


!!! Test the transposition.
!!allocate(psi(orbsLIN%npsidim), stat=istat)
!!allocate(psiWork(orbsLIN%npsidim), stat=istat)
!!call random_number(psi)
!!write(100+iproc,*) psi
!!do igroup=1,ngroups
!!    lproc=procsInGroup(1,igroup)
!!    uproc=procsInGroup(2,igroup)
!!    if(iproc>=lproc .and. iproc<=uproc) then
!!        call transpose_vLIN(iproc, lproc, uproc, norbPerGroupArr(igroup), orbsLIN, commsLIN, psi, lr, newComm(igroup), work=psiWork)
!!        write(200+iproc,*) psi
!!        call orthogonalizeLIN(iproc, lproc, uproc, norbPerGroupArr(igroup), orbsLIN, commsLIN, psi, input, newComm(igroup))
!!        call untranspose_vLIN(iproc, lproc, uproc, norbPerGroupArr(igroup), orbsLIN, commsLIN, psi, lr, newComm(igroup), work=psiWork)
!!    end if
!!end do
!!write(300+iproc,*) psi
!!nullify(psi)
!!nullify(psiWork)




end subroutine initializeLocRegLIN




!!****f* BigDFT/fill_logrid
!! FUNCTION
!!   set up an array logrid(i1,i2,i3) that specifies whether the grid point
!!   i1,i2,i3 is the center of a scaling function/wavelet
!!
!! SOURCE
!!
subroutine fill_logridCut(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
  real(gp), intent(in) :: rmult,hx,hy,hz
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(ntypes), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
  real(gp) :: dx,dy2,dz2,rad

!write(*,'(a,3i15,4x,3i15,3es12.3)') 'at start fill_logridCut: nl1, nu1, nl2, nu2, nl3, nu3, hx, hy, hz', nl1, nu1, nl2, nu2, nl3, nu3, hx, hy, hz
!write(*,'(a,3es15.5)') 'rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)', rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)

  !some checks
  if (geocode /='F') then
     !the nbuf value makes sense only in the case of free BC
     if (nbuf /=0) then
        write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
        stop
     end if
     !the grid spacings must be the same
     if (hx/= hy .or. hy /=hz .or. hx/=hz) then
!        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
     end if
  end if

  if (geocode == 'F') then
     do i3=nl3,nu3
        do i2=nl2,nu2
           do i1=nl1,nu1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  else !
     do i3=0,n3
        do i2=0,n2
           do i1=0,n1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  end if

  do iat=1,nat
     rad=radii(iatype(iat))*rmult+real(nbuf,gp)*hx
     if (rad /= 0.0_gp) then
        ml1=max(ceiling((rxyz(1,iat)-rad)/hx - eps_mach), nl1)
        ml2=max(ceiling((rxyz(2,iat)-rad)/hy - eps_mach), nl2)
        ml3=max(ceiling((rxyz(3,iat)-rad)/hz - eps_mach), nl3)
        mu1=min(floor((rxyz(1,iat)+rad)/hx + eps_mach), nu1)
        mu2=min(floor((rxyz(2,iat)+rad)/hy + eps_mach), nu2)
        mu3=min(floor((rxyz(3,iat)+rad)/hz + eps_mach), nu3)
        !write(*,'(a,3es15.5)') 'rxyz(1,iat)+rad, rxyz(2,iat)+rad, rxyz(3,iat)+rad', rxyz(1,iat)+rad, rxyz(2,iat)+rad, rxyz(3,iat)+rad
        !write(*,'(a,f9.2,3es12.4)') 'rxyz(1,iat), rxyz(2,iat), rxyz(3,iat), rad', rxyz(1,iat), rxyz(2,iat), rxyz(3,iat), rad
        !for Free BC, there must be no incoherences with the previously calculated delimiters
        if (geocode == 'F') then
           if (ml1 < nl1) stop 'ml1 < nl1'
           if (ml2 < nl2) stop 'ml2 < nl2'
           if (ml3 < nl3) stop 'ml3 < nl3'

           if (mu1 > nu1) stop 'mu1 > nu1'
           if (mu2 > nu2) stop 'mu2 > nu2'
           if (mu3 > nu3) stop 'mu3 > nu3'
        end if
        !write(*,'(a,f6.1,3i15,4x,3i15)') 'rad, ml1, ml2, ml3, mu1, mu2, mu3', rad, ml1, ml2, ml3, mu1, mu2, mu3
        !what follows works always provided the check before
!$omp parallel default(shared) private(i3,dz2,j3,i2,dy2,j2,i1,j1,dx)
!$omp do
        do i3=max(ml3,-n3/2-1),min(mu3,n3+n3/2+1)
           dz2=(real(i3,gp)*hz-rxyz(3,iat))**2
           j3=modulo(i3,n3+1)
           do i2=max(ml2,-n2/2-1),min(mu2,n2+n2/2+1)
              dy2=(real(i2,gp)*hy-rxyz(2,iat))**2
              j2=modulo(i2,n2+1)
              do i1=max(ml1,-n1/2-1),min(mu1,n1+n1/2+1)
                 j1=modulo(i1,n1+1)
                 dx=real(i1,gp)*hx-rxyz(1,iat)
                 if (dx**2+(dy2+dz2) <= rad**2) then
                    !if(j1<0 .or. j1>n1 .or. j2<0 .or. j2>n2 .or. j3<0 .or. j3>n3) write(*,'(a,f6.1,3i6,2x,3i6)') 'dx**2+(dy2+dz2), j1, j2, j3, n1, n2, n3', dx**2+(dy2+dz2), j1, j2, j3, n1, n2, n3
                    logrid(j1,j2,j3)=.true.
                 endif
              enddo
           enddo
        enddo
!$omp enddo
!$omp end parallel
  end if
  enddo

END SUBROUTINE fill_logridCut
!!***









!!****f* BigDFT/orbitals_communicatorsLIN_group
!! FUNCTION
!!   Partition the orbitals between processors to ensure load balancing
!!   the criterion will depend on GPU computation
!!   and/or on the sizes of the different localisation region
!! DESCRIPTION
!!   Calculate the number of elements to be sent to each process
!!   and the array of displacements
!!   Cubic strategy: 
!!      - the components are equally distributed among the wavefunctions
!!      - each processor has all the orbitals in transposed form
!!      - each wavefunction is equally distributed in its transposed form
!!      - this holds for each k-point, which regroups different processors
!!
!! SOURCE
!!
subroutine orbitals_communicatorsLIN_group(iproc, lproc, uproc, lr, orbs, comms, newComm, norbPerGroup)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, lproc, uproc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(out) :: comms
integer:: newComm
integer:: norbPerGroup
  !local variables
  character(len=*), parameter :: subname='orbitals_communicators'
  logical :: yesorb,yescomp
  integer :: jproc,nvctr_tot,ikpts,iorbp,jorb,norb_tot,ikpt,i_stat,i_all
  integer :: ncomp_res,nkptsp,ierr
  integer, dimension(:), allocatable :: mykpts
  logical, dimension(:), allocatable :: GPU_for_comp
  integer, dimension(:,:), allocatable :: nvctr_par,norb_par !for all the components and orbitals (with k-pts)
  integer, dimension(:,:,:,:), allocatable :: nvctr_parLIN, nvctr_parLIN2 !for all the components and orbitals (with k-pts)
integer:: iorb, ii, outproc, nproc
integer,dimension(:,:),allocatable:: ncntdLIN, ndspldLIN, ncnttLIN, ndspltLIN
  
    nproc=uproc-lproc+1


  !calculate the number of elements to be sent to each process
  !and the array of displacements
  !cubic strategy: -the components are equally distributed among the wavefunctions
  !                -each processor has all the orbitals in transposed form
  !                -each wavefunction is equally distributed in its transposed form
  !                -this holds for each k-point, which regroups different processors

  !check of allocation of important arrays
  if (.not. associated(orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if
  

  !allocate(nvctr_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  allocate(nvctr_par(lproc:uproc,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,nvctr_par,'nvctr_par',subname)
  !allocate(norb_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  allocate(norb_par(lproc:uproc,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)
  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,mykpts,'mykpts',subname)

  !initialise the arrays
  do ikpts=0,orbs%nkpts
     !do jproc=0,nproc-1
     do jproc=lproc,uproc
        nvctr_par(jproc,ikpts)=0 
        norb_par(jproc,ikpts)=0 
     end do
  end do


 !! THIS PART IS COPIED FROM BELOW !!
  jorb=1
  ikpts=1
  !print *,'here',orbs%norb_par(:)
  !do jproc=0,nproc-1
  do jproc=lproc,uproc
     do iorbp=1,orbs%norb_par(jproc)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (mod(jorb,orbs%norb)==0) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do

  !allocate(nvctr_parLIN2(maxval(norb_par(:,1)),0:nproc-1,0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  allocate(nvctr_parLIN2(maxval(norb_par(:,1)),lproc:uproc,lproc:uproc,0:orbs%nkpts+ndebug),stat=i_stat)
  !allocate(nvctr_parLIN(maxval(norb_par(:,1)),0:nproc-1,0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  allocate(nvctr_parLIN(maxval(norb_par(:,1)),lproc:uproc,lproc:uproc,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,nvctr_par,'nvctr_parLIN',subname)
  nvctr_parLIN2=0
  nvctr_parLIN=0

  !balance the components between processors
  !in the most symmetric way
  !here the components are taken into account for all the k-points

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines
  allocate(GPU_for_comp(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,GPU_for_comp,'GPU_for_comp',subname)

  if (nproc > 1 .and. .not. GPUshare) then
     call MPI_ALLGATHER(GPUblas,1,MPI_LOGICAL,GPU_for_comp(0),1,MPI_LOGICAL,&
          MPI_COMM_WORLD,ierr)
  else
     GPU_for_comp(0)=GPUblas
  end if

  i_all=-product(shape(GPU_for_comp))*kind(GPU_for_comp)
  deallocate(GPU_for_comp,stat=i_stat)
  call memocc(i_stat,i_all,'GPU_for_comp',subname)


  !decide the repartition for the components in the same way as the orbitals
  !call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par)
  !do jproc=0,nproc-1
     do iorb=1,orbs%norbp
         !!! lr%wfdCut(iorb)%nvctr_c and lr%wfdCut(iorb)%nvctr_f NOT YET ALLOCATED AT CALLING TIME !!
         !call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par)
         write(*,*) 'calling parallel_repartition_with_kpoints, iproc, iorb',iproc, iorb
         !call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfdCut(iorb,iproc)%nvctr_c+7*lr%wfdCut(iorb,iproc)%nvctr_f),nvctr_par)
         !call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfdLIN(iorb,iproc)%nvctr_c+7*lr%wfdLIN(iorb,iproc)%nvctr_f),nvctr_par)
         call parallel_repartition_with_kpoints(nproc,orbs%nkpts,&
             (lr%wfdLIN(iorb,iproc)%nvctr_c+7*lr%wfdLIN(iorb,iproc)%nvctr_f),nvctr_par)
         !do jproc=0,nproc-1
         do jproc=lproc,uproc
             !nvctr_parLIN(iorb,0:nproc-1,0)=nvctr_par(jproc,0)
             nvctr_parLIN2(iorb,iproc,jproc,0)=nvctr_par(jproc,0)
         end do
     end do
  !end do
!call mpi_allreduce(nvctr_parLIN2(1,0,0,0), nvctr_parLIN(1,0,0,0), maxval(norb_par(:,1))*nproc*nproc*(orbs%nkpts+1), mpi_integer, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(nvctr_parLIN2(1,lproc,lproc,0), nvctr_parLIN(1,lproc,lproc,0), &
     maxval(norb_par(:,1))*nproc*nproc*(orbs%nkpts+1), mpi_integer, mpi_sum, newComm, ierr)
call mpi_barrier(newComm, ierr)
write(*,*) 'after loop calling parallel_repartition_with_kpoints, iproc', iproc

!!!!do outproc=0,nproc-1
!!!do outproc=lproc,uproc
!!!  do iorb=1,orbs%norbp
!!!      write(*,'(a,i3,3x,8i8)') 'LIN: iproc, nvctr_parLIN(iorb,outproc,0:nproc-1,0)',iproc, nvctr_parLIN(iorb,outproc,0:nproc-1,0)
!!!  end do
!!!end do

!!$  i=1
!!$  j=1
!!$  loop_components: do 
!!$     jproc=mod(i-1,nproc)
!!$     if (.true.) then !here there is the criterion for filling a processor
!!$        nvctr_par(jproc,0)=nvctr_par(jproc,0)+1
!!$        j=j+1
!!$     end if
!!$     if (j > (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nkpts) exit loop_components
!!$     i=i+1
!!$  end do loop_components


  !write(*,'(a,2i7)') 'size(nvctr_parLIN, 2), nproc', size(nvctr_parLIN, 2), nproc
  !do outproc=0,nproc-1
  do outproc=lproc,uproc
    !do iorb=1,orbs%norbp
    do iorb=1,norb_par(outproc,1) ! index 1 for k-points
        ikpts=1
        write(*,'(a,2i6)') 'increasing iorb: iproc, iorb', iproc, iorb
        !ncomp_res=(lr%wfdCut(iorb,outproc)%nvctr_c+7*lr%wfdCut(iorb,outproc)%nvctr_f)
        ncomp_res=(lr%wfdLIN(iorb,outproc)%nvctr_c+7*lr%wfdLIN(iorb,outproc)%nvctr_f)
        !do jproc=0,nproc-1
        do jproc=lproc,uproc
           loop_comps: do
              !write(*,'(a,4i12)') 'iproc, nvctr_parLIN(iorb,jproc,0), ncomp_res, ikpts, ', iproc, nvctr_parLIN(iorb,jproc,0), ncomp_res, ikpts
              if(ikpts>1) write(*,*) 'ATTENTION: iproc',iproc
              if (nvctr_parLIN(iorb,outproc,jproc,0) >= ncomp_res) then
                 nvctr_parLIN(iorb,outproc,jproc,ikpts)= ncomp_res
                 ikpts=ikpts+1
                 nvctr_parLIN(iorb,outproc,jproc,0)=nvctr_parLIN(iorb,outproc,jproc,0)-ncomp_res
                 ncomp_res=(lr%wfdLIN(iorb,outproc)%nvctr_c+7*lr%wfdLIN(iorb,outproc)%nvctr_f)
              else
                 nvctr_parLIN(iorb,outproc,jproc,ikpts)= nvctr_parLIN(iorb,outproc,jproc,0)
             if(nvctr_parLIN(iorb,outproc,jproc,ikpts)==0) write(*,'(a,3i6)') 'ATTENTION: iorb, outproc, jproc', iorb, &
                     outproc, jproc
                 ncomp_res=ncomp_res-nvctr_parLIN(iorb,outproc,jproc,0)
                 nvctr_parLIN(iorb,outproc,jproc,0)=0
                 exit loop_comps
              end if
              if (nvctr_parLIN(iorb,outproc,jproc,0) == 0 ) then
                 ncomp_res=(lr%wfdLIN(iorb,outproc)%nvctr_c+7*lr%wfdLIN(iorb,outproc)%nvctr_f)
                 exit loop_comps
              end if
    
           end do loop_comps
      end do
   end do
 end do
call mpi_barrier(newComm, ierr)
!do outproc=0,nproc-1
do outproc=lproc,uproc
  do iorb=1,orbs%norbp
      !write(*,'(a,i3,3x,8i8)') 'LIN: iproc, nvctr_par(iorb,outproc,lproc:uproc,1)',iproc, nvctr_parLIN(iorb,outproc,lproc:uproc,1)
  end do
end do
write(*,*) 'PERFORMING CHECKS..., iproc', iproc
 
   !some checks
   !check the distribution
   do ikpts=1,orbs%nkpts
      !print *,'iproc,cpts:',lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nvctr_par(:,ikpts)
     !do outproc=0,nproc-1
     do outproc=lproc,uproc
       do iorb=1,orbs%norbp
          nvctr_tot=0
          !do jproc=0,nproc-1
          do jproc=lproc,uproc
              nvctr_tot=nvctr_tot+nvctr_parLIN(iorb,outproc,jproc,ikpts)
          end do
          if(nvctr_tot /= lr%wfdLIN(iorb,outproc)%nvctr_c+7*lr%wfdLIN(iorb,outproc)%nvctr_f) then
             write(*,*)'ERROR: partition of components incorrect, iorb, kpoint:',iorb, ikpts
             stop
          end if
       end do
     end do
   end do
 
 
   !!  DO THIS LATER ... OR ALREADY DONE BY NORMAL SUBROUTINE? !!!
   !this function which associates a given k-point to a processor in the component distribution
   !the association is chosen such that each k-point is associated to only
   !one processor
   !if two processors treat the same k-point the processor which highest rank is chosen
   do ikpts=1,orbs%nkpts
      !loop_jproc: do jproc=nproc-1,0,-1
      loop_jproc: do jproc=uproc,lproc,-1
         if (nvctr_par(jproc,ikpts) /= 0) then
            orbs%ikptproc(ikpts)=jproc
            exit loop_jproc
         end if
      end do loop_jproc
   end do
   
   !print *,'check',orbs%ikptproc(:)
 
 
   !!  DO THIS LATER ... OR ALREADY DONE BY NORMAL SUBROUTINE? !!!
   !calculate the same k-point distribution for the orbitals
   !assign the k-point to the given orbital, counting one orbital after each other

   !! ALREADY DONE ABOVE !!
   !!jorb=1
   !!ikpts=1
   !!!print *,'here',orbs%norb_par(:)
   !!do jproc=0,nproc-1
   !!   do iorbp=1,orbs%norb_par(jproc)
   !!      norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
   !!      if (mod(jorb,orbs%norb)==0) then
   !!         ikpts=ikpts+1
   !!      end if
   !!      jorb=jorb+1
   !!   end do
   !!end do
   !some checks
   if (orbs%norb /= 0) then
      !check the distribution
      do ikpts=1,orbs%nkpts
         !print *,'partition',ikpts,orbs%nkpts,'ikpts',norb_par(:,ikpts)
         norb_tot=0
         !do jproc=0,nproc-1
         do jproc=lproc,uproc
            norb_tot=norb_tot+norb_par(jproc,ikpts)
         end do
         !if(norb_tot /= orbs%norb) then
         if(norb_tot /= norbPerGroup) then
            write(*,*)'ERROR: partition of orbitals incorrect, kpoint:',ikpts
            write(*,*) 'norb_tot, orbs%norb', norb_tot, orbs%norb
            stop
         end if
      end do
   end if
 
 
   !!  DO THIS LATER ... OR ALREADY DONE BY NORMAL SUBROUTINE? !!!
   !calculate the number of k-points treated by each processor in both
   ! the component distribution and the orbital distribution.
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
 
   !!  DO THIS LATER ... OR ALREADY DONE BY NORMAL SUBROUTINE? !!!
   !print the distribution scheme ussed for this set of orbital
   !in the case of multiple k-points
   if (iproc == 0 .and. verbose > 1 .and. orbs%nkpts > 1) then
      call print_distribution_schemes(nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
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
         stop
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
         stop
      end if
   end do kpt_orbitals


 
   !allocate communication arrays
   !allocate(comms%nvctr_par(0:nproc-1,orbs%nkptsp+ndebug),stat=i_stat)
   !call memocc(i_stat,comms%nvctr_par,'nvctr_par',subname)
 
   !allocate(comms%nvctr_parLIN(0:nproc-1,1:orbs%norb,orbs%nkptsp+ndebug),stat=i_stat)
   !call memocc(i_stat,comms%nvctr_parLIN,'nvctr_parLIN',subname)
   !allocate(comms%nvctr_parLIN(1:orbs%norb,0:nproc-1,0:nproc-1,orbs%nkptsp+ndebug),stat=i_stat)
   allocate(comms%nvctr_parLIN(1:orbs%norb,lproc:uproc,lproc:uproc,orbs%nkptsp+ndebug),stat=i_stat)
   call memocc(i_stat,comms%nvctr_parLIN,'nvctr_parLIN',subname)
 
   !!allocate(comms%ncntd(0:nproc-1+ndebug),stat=i_stat)
   !!call memocc(i_stat,comms%ncntd,'ncntd',subname)
 
   !!allocate(comms%ncntt(0:nproc-1+ndebug),stat=i_stat)
   !!call memocc(i_stat,comms%ncntt,'ncntt',subname)
   !!allocate(comms%ndspld(0:nproc-1+ndebug),stat=i_stat)
   !!call memocc(i_stat,comms%ndspld,'ndspld',subname)
   !!allocate(comms%ndsplt(0:nproc-1+ndebug),stat=i_stat)
   !!call memocc(i_stat,comms%ndsplt,'ndsplt',subname)
 
   !assign the partition of the k-points to the communication array
   do ikpts=1,orbs%nkptsp
      ikpt=orbs%iskpts+ikpts!orbs%ikptsp(ikpts)
      !do jproc=0,nproc-1
      do jproc=lproc,uproc
         !do outproc=0,nproc-1
         do outproc=lproc,uproc
            do iorb=1,norb_par(outproc,ikpt)
               comms%nvctr_parLIN(iorb,outproc,jproc,ikpt)=nvctr_parLIN(iorb,outproc,jproc,ikpt) 
            end do
         end do
      end do
   end do
 

   !with this distribution the orbitals and the components are ordered following k-points
   !there must be no overlap for the components
   !here we will print out the k-points components distribution, in the transposed and in the direct way
 
 
   ! send counts for mpi_alltoallv
   !allocate(comms%ncntdLIN(0:nproc-1,0:nproc-1), stat=i_stat)
   !allocate(comms%ncntdLIN(0:nproc-1), stat=i_stat)
   allocate(comms%ncntdLIN(lproc:uproc), stat=i_stat)
   comms%ncntdLIN=0
   !allocate(ncntdLIN(0:nproc-1,0:nproc-1), stat=i_stat)
   allocate(ncntdLIN(lproc:uproc,lproc:uproc), stat=i_stat)
   ncntdLIN=0
   call memocc(i_stat,comms%ncntdLIN,'comms%ncntdLIN',subname)
   !do jproc=0,nproc-1
   do jproc=lproc,uproc
      !do iorb=1,orbs%norb
         do ikpts=1,orbs%nkpts
            ii=0
            do jorb=1,norb_par(iproc,ikpts)
                ii=ii+nvctr_parLIN(jorb,iproc,jproc,ikpts)*orbs%nspinor
            end do
            !ncntdLIN(jproc,iproc)=ncntdLIN(jproc,iproc)+ii
            comms%ncntdLIN(jproc)=comms%ncntdLIN(jproc)+ii
            !comms%ncntdLIN(iorb,jproc)=comms%ncntdLIN(iorb,jproc)+&
            !     nvctr_parLIN(iorb,jproc,ikpts)*norb_par(iproc,ikpts)*orbs%nspinor
         end do
      !end do
   end do
 
   ! send displacements for mpi_alltoallv
   !allocate(comms%ndspldLIN(0:nproc-1,0:nproc-1))
   !allocate(comms%ndspldLIN(0:nproc-1))
   allocate(comms%ndspldLIN(lproc:uproc))
   comms%ndspldLIN=0
   !allocate(ndspldLIN(0:nproc-1,0:nproc-1))
   !allocate(ndspldLIN(lproc:uproc))
   !ndspldLIN=0
   call memocc(i_stat,comms%ndspldLIN,'comms%ndspldLIN',subname)
   !do iorb=1,orbs%norb
      !do jproc=1,nproc-1
      do jproc=lproc+1,uproc
         !ndspldLIN(jproc,iproc)=ndspldLIN(jproc-1,iproc)+ncntdLIN(jproc-1,iproc)
         comms%ndspldLIN(jproc)=comms%ndspldLIN(jproc-1)+comms%ncntdLIN(jproc-1)
      end do
   !end do
 
 
 !write(*,'(a,i4,4x,8i5)') 'iproc, norb_par(0:nproc-1,1)',iproc, norb_par(0:nproc-1,1)
   !receive buffer
   ! receive counts for mpi_alltoallv
   !allocate(comms%ncnttLIN(0:nproc-1,0:nproc-1), stat=i_stat)
   !allocate(comms%ncnttLIN(0:nproc-1), stat=i_stat)
   allocate(comms%ncnttLIN(lproc:uproc), stat=i_stat)
   comms%ncnttLIN=0
   !allocate(ncnttLIN(0:nproc-1,0:nproc-1), stat=i_stat)
   allocate(ncnttLIN(lproc:uproc,lproc:uproc), stat=i_stat)
   ncnttLIN=0
   call memocc(i_stat,comms%ncnttLIN,'comms%ncnttLIN',subname)
   !do jproc=0,nproc-1
   do jproc=lproc,uproc
      !do iorb=1,orbs%norb
         do ikpts=1,orbs%nkpts
            ii=0
            do jorb=1,norb_par(jproc,ikpts)
                ii=ii+nvctr_parLIN(jorb,jproc,iproc,ikpts)*orbs%nspinor
            end do
            !ncnttLIN(jproc,iproc)=ncnttLIN(jproc,iproc)+ii
            comms%ncnttLIN(jproc)=comms%ncnttLIN(jproc)+ii
            !comms%ncnttLIN(iorb,jproc)=comms%ncnttLIN(iorb,jproc)+&
            !     nvctr_parLIN(iorb,iproc,ikpts)*norb_par(jproc,ikpts)*orbs%nspinor
         end do
      !end do
   end do
   
   ! receive displacements for mpi_alltoallv
   !allocate(comms%ndspltLIN(0:nproc-1,0:nproc-1), stat=i_stat)
   !allocate(comms%ndspltLIN(0:nproc-1), stat=i_stat)
   allocate(comms%ndspltLIN(lproc:uproc), stat=i_stat)
   comms%ndspltLIN=0
   !allocate(ndspltLIN(0:nproc-1,0:nproc-1), stat=i_stat)
   allocate(ndspltLIN(lproc,uproc), stat=i_stat)
   ndspltLIN=0
   call memocc(i_stat,comms%ndspltLIN,'comms%ndspltLIN',subname)
   !do iorb=1,orbs%norb
      !do jproc=1,nproc-1
      do jproc=lproc+1,uproc
         !ndspltLIN(jproc,iproc)=ndspltLIN(jproc-1,iproc)+ncnttLIN(jproc-1,iproc)
         comms%ndspltLIN(jproc)=comms%ndspltLIN(jproc-1)+comms%ncnttLIN(jproc-1)
      end do
   !end do
 
  !!! collect parameters for mpi_alltoallv from all processes
  !!call mpi_allreduce(ncntdLIN(0,0), comms%ncntdLIN(0,0), nproc**2, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  !!call mpi_allreduce(ndspldLIN(0,0), comms%ndspldLIN(0,0), nproc**2, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  !!call mpi_allreduce(ncnttLIN(0,0), comms%ncnttLIN(0,0), nproc**2, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  !!call mpi_allreduce(ndspltLIN(0,0), comms%ndspltLIN(0,0), nproc**2, mpi_integer, mpi_sum, mpi_comm_world, ierr)

  !comms%ncntdLIN(0:nproc-1)=ncntdLIN(0:nproc-1,iproc)
  !comms%ndspldLIN(0:nproc-1)=ndspldLIN(0:nproc-1,iproc)
  !comms%ncnttLIN(0:nproc-1)=ncnttLIN(0:nproc-1,iproc)
  !comms%ndspltLIN(0:nproc-1)=ndspltLIN(0:nproc-1,iproc)
 
 
   !print *,'iproc,comms',iproc,comms%ncntd,comms%ndspld,comms%ncntt,comms%ndsplt
 
   i_all=-product(shape(nvctr_par))*kind(nvctr_par)
   deallocate(nvctr_par,stat=i_stat)
   call memocc(i_stat,i_all,'nvctr_par',subname)
   i_all=-product(shape(norb_par))*kind(norb_par)
   deallocate(norb_par,stat=i_stat)
   call memocc(i_stat,i_all,'norb_par',subname)
   i_all=-product(shape(mykpts))*kind(mykpts)
   deallocate(mykpts,stat=i_stat)
   call memocc(i_stat,i_all,'mykpts',subname)
 
   !calculate the dimension of the wavefunction
   !for the given processor
   !take into account max one k-point per processor
   !!if(.not. associated(orbs%npsidimLIN)) then
   !!    allocate(orbs%npsidimLIN(1:orbs%norb))
   !!end if
 
   !!! DO THIS LATER !!!
   !!do iorb=1,orbs%norbp
   !!    orbs%npsidimLIN(iorb)=max((lr%wfdCut(iorb)%nvctr_c+7*lr%wfdCut(iorb)%nvctr_f)*orbs%norb_par(iproc)*orbs%nspinor,&
   !!         sum(comms%ncnttLIN(iorb,0:nproc-1)))
     !orbs%npsidimLIN=0
     write(*,*) 'ATTENTION: changing orbs%npsidim!!'
     orbs%npsidim=0
     do iorb=1,orbs%norbp
         !orbs%npsidimLIN=orbs%npsidimLIN+(lr%wfdCut(iorb,iproc)%nvctr_c+7*lr%wfdCut(iorb,iproc)%nvctr_f)*orbs%nspinor
         !orbs%npsidimLIN=orbs%npsidimLIN+(lr%wfdLIN(iorb,iproc)%nvctr_c+7*lr%wfdLIN(iorb,iproc)%nvctr_f)*orbs%nspinor
         orbs%npsidim=orbs%npsidim+(lr%wfdLIN(iorb,iproc)%nvctr_c+7*lr%wfdLIN(iorb,iproc)%nvctr_f)*orbs%nspinor
     end do
     !orbs%npsidimLIN=max(orbs%npsidimLIN,sum(comms%ncnttLIN(0:nproc-1,iproc)))
     !orbs%npsidimLIN=max(orbs%npsidimLIN,sum(comms%ncnttLIN(0:nproc-1)))
     !orbs%npsidimLIN=max(orbs%npsidimLIN,sum(comms%ncnttLIN(lproc:uproc)))
     orbs%npsidim=max(orbs%npsidim,sum(comms%ncnttLIN(lproc:uproc)))
   !!end do
 
   if (iproc == 0) write(*,'(1x,a,i0)') &
        'Wavefunctions memory occupation for root processor (Bytes): ',&
        orbs%npsidim*8
   write(*,'(1x,a,i0,4x,i5)') &
        'LIN: Wavefunctions memory occupation for root processor (Bytes), iproc ',&
        !orbs%npsidimLIN*8, iproc
        orbs%npsidim*8, iproc
 
 
 call mpi_barrier(newComm, ierr)
 write(*,*) 'at end orbitals_communicatorsLIN, iproc', iproc
 call mpi_barrier(newComm, ierr)

END SUBROUTINE orbitals_communicatorsLIN_group
!!***
