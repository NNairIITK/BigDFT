module libBigDFT

  private

!- High level methods.
  !- Do a minimisation loop to find forces and energy.
  public :: cluster
  !- Do a DIIS step on wavefunctions.
  public :: diisstp

  !- Compute the density from the wavefunctions.
  public :: sumrho
  !- Apply the kinetic operator and the local part of potential
  public :: applylocpotkinone
  public :: applylocpotkinall
  !- Apply the non local part of potential
  public :: applyprojectorsone
  public :: applyprojectorsall
  !- Orthogonalize wavefunctions (serial and parallel MPI routine)
  public :: orthoconstraint, orthoconstraint_p
  public :: orthon, orthon_p
  !- Precondition on wavefunctions
!  public :: preconditionall
  !- Make some Loewedin orthogonalisations (serial and parallel MPI).
  public :: loewe, loewe_p
  !- Control orthogonality (serial and parallel MPI).
  public :: checkortho, checkortho_p
  !- Solve a KS system (compute an hamiltonian, diagonalise it and return eigen vectors)
  public :: KStrans, KStrans_p, solveKS
  
  !- Initialisation methods.
  !- Create and allocate access arrays for wavefunctions.
  public :: createWavefunctionsDescriptors
  !- Create and allocate projectors (and their access arrays).
  public :: createProjectorsArrays
  public :: crtproj
  !- Create a ionic density (without gaussian parts)
  public :: input_rho_ion
  !- Add gaussian part to local part of ionic potential
  public :: addlocgauspsp
  !- Create the kernel for Poisson solver.
  public :: createKernel
  !- Compute input guess wavefunctions from aatomic orbitals.
  public :: input_wf_diag
  
  !- Wavefunctions related methods.
  !- Create wavefunctions from atomic orbitals
  public :: createAtomicOrbitals
  !- Read a file with description of atomic orbitals for all elements
  public :: readAtomicOrbitals
  !- Read wavefunctions from disk
  public :: readmywaves
  !- Atomic read for a given wavefunction
  public :: readonewave
  !- Write wavefunctions to disk using the BigDFT format
  public :: writemywaves
  !- Atomic write for a given wavefunction
  public :: writeonewave
  !- Transform wavefunctions from old grid to new grid.
  public :: reformatmywaves
  !- MPI communications for wavefunctions
  public :: transallwaves
  public :: untransallwaves

  !- Geometry method
  !- Get a box that contains all atoms and their active grid points.
  public :: system_size
  
  !- Linear algebra method
  !- Compute a scalar product
  public :: wnrm
  public :: wpdot

contains

subroutine cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames, rxyz, energy, fxyz,  &
     & psi, keyg, keyv, nvctr_c, nvctr_f, nseg_c, nseg_f, norbp, norb, eval, inputPsiId, &
     & output_grid, output_wf, n1, n2, n3, hgrid, rxyz_old)
  ! inputPsiId = 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
  ! inputPsiId = 1 : read waves from argument psi, using n1, n2, n3, hgrid and rxyz_old
  !                  as definition of the previous system.
  ! inputPsiId = 2 : read waves from disk
  ! does an electronic structure calculation. Output is the total energy and the forces 
  ! psi, keyg, keyv and eval should be freed after use outside of the routine.
  implicit real*8 (a-h,o-z)
  character*30 label ; character*27 filename ; character*20 atomnames
  logical logrid_c,logrid_f,parallel,calc_tail,output_wf,output_grid
  parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
  ! work array for ALLREDUCE
  dimension wrkallred(5,2) 
  ! atomic coordinates
  dimension rxyz(3,nat),rxyz_old(3,nat),fxyz(3,nat),iatype(nat),atomnames(100)
  allocatable :: gxyz(:,:)
  ! active grid points, segments of real space grid
  allocatable :: logrid_c(:,:,:) ,  logrid_f(:,:,:)
  allocatable :: ibyz_c(:,:,:),ibxz_c(:,:,:),ibxy_c(:,:,:),  & 
       ibyz_f(:,:,:),ibxz_f(:,:,:),ibxy_f(:,:,:)
  ! occupation numbers, eigenvalues
  allocatable :: occup(:)
  real*8, pointer :: eval(:),eval_old(:)

  ! wavefunction segments
  integer, pointer :: keyv(:)
  ! wavefunction segments on real space grid
  integer, pointer :: keyg(:,:)
  ! wavefunction 
  real*8, pointer :: psi(:,:)
  real*8, pointer :: psit(:,:)
  ! wavefunction gradients
  real*8, pointer :: hpsi(:,:),hpsit(:,:)

  ! Pointers and variables to store the last psi
  ! before reformating if useFormattedInput is .true.
  integer :: nseg_c_old, nseg_f_old, nvctr_c_old, nvctr_f_old
  integer, pointer :: keyg_old(:,:), keyv_old(:)
  real*8, pointer :: psi_old(:,:)

  ! Charge density/potential,ionic potential, pkernel
  allocatable :: rhopot(:,:,:),pot_ion(:)
  real*8, pointer     :: pkernel(:)

  ! projector segments on real space grid
  pointer :: keyg_p(:,:), keyv_p(:)
  allocatable :: nvctr_p(:), nseg_p(:)
  ! projectors 
  real*8, pointer :: proj(:)
  ! Parameters for the boxes containing the projectors
  allocatable :: nboxp_c(:,:,:),nboxp_f(:,:,:)

  ! pseudopotential parameters
  allocatable :: psppar(:,:,:),nelpsp(:),radii_cf(:,:),npspcode(:)

  ! arrays for DIIS convergence accelerator
  real*8, pointer :: ads(:,:,:),psidst(:,:,:),hpsidst(:,:,:)

  ! arrays for correction to kinetic energy
  allocatable :: ibbyz_c(:,:,:),ibbxz_c(:,:,:),ibbxy_c(:,:,:),  & 
       ibbyz_f(:,:,:),ibbxz_f(:,:,:),ibbxy_f(:,:,:)
  allocatable :: keybv(:),keybg(:,:),txyz(:,:)
  allocatable :: psib(:),hpsib(:)
  allocatable :: psifscf(:),psir(:),psig(:),psigp(:)
  allocatable :: rho(:),hartpot(:),rhopotb(:,:,:)

  integer :: ierror

  include 'mpif.h'

  if (iproc.eq.0) write(*,*) 'CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER',inputPsiId
  if (parallel) then
     call timing(iproc,'parallel     ','IN')
  else
     call timing(iproc,'             ','IN')
  end if
  if (iproc.eq.0) open(unit=79,file='malloc')
  call cpu_time(tcpu0)
  call system_clock(ncount0,ncount_rate,ncount_max)

  ! We save the variables that defined the previous psi if
  ! restartOnPsi is .true.
  if (inputPsiId == 1) then
     hgrid_old   = hgrid
     n1_old      = n1
     n2_old      = n2
     n3_old      = n3
     nvctr_c_old = nvctr_c
     nvctr_f_old = nvctr_f
     nseg_c_old  = nseg_c
     nseg_f_old  = nseg_f
     allocate(keyg_old(2,nseg_c_old+nseg_f_old))
     allocate(keyv_old(nseg_c_old+nseg_f_old))
     allocate(psi_old(nvctr_c_old+7*nvctr_f_old,norbp))
     allocate(eval_old(norb))
     do iseg=1,nseg_c_old+nseg_f_old
        keyg_old(1,iseg)    = keyg(1,iseg)
        keyg_old(2,iseg)    = keyg(2,iseg)
        keyv_old(iseg)      = keyv(iseg)
     enddo
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        tt=0.d0
        do j=1,nvctr_c_old+7*nvctr_f_old
           psi_old(j,iorb-iproc*norbp)     = psi(j,iorb-iproc*norbp)
           tt=tt+psi(j,iorb-iproc*norbp)**2
        enddo
        tt=sqrt(tt)
        if (abs(tt-1.d0).gt.1.d-8) stop 'wrong psi_old'
        eval_old(iorb) = eval(iorb)
     enddo
     deallocate(keyg, keyv)
     deallocate(psi, eval)
  end if

  ! Read the input variables.
  open(unit=1,file='input.dat',status='old')
  !First line for the main routine (the program)
  read(1,*)
  !Parameters 
  read(1,*) hgrid
  read(1,*) crmult
  read(1,*) frmult
  read(1,*) cpmult
  read(1,*) fpmult
  if (fpmult.gt.frmult) write(*,*) 'NONSENSE: fpmult > frmult'
  read(1,*) gnrm_cv
  read(1,*) itermax
  read(1,*) ncong
  read(1,*) idsx
  read(1,*) calc_tail
  read(1,*) rbuf
  read(1,*) ncongt
  close(1)

  if (iproc.eq.0) then 
     write(*,*) 'hgrid=',hgrid
     write(*,*) 'crmult=',crmult
     write(*,*) 'frmult=',frmult
     write(*,*) 'cpmult=',cpmult
     write(*,*) 'fpmult=',fpmult
     write(*,*) 'gnrm_cv=',gnrm_cv
     write(*,*) 'itermax=',itermax
     write(*,*) 'ncong=',ncong
     write(*,*) 'idsx=',idsx
     write(*,*) 'calc_tail',calc_tail
     write(*,*) 'rbuf=',rbuf
     write(*,*) 'ncongt=',ncongt
  endif

! grid spacing (same in x,y and z direction)
  hgridh=.5d0*hgrid

! store PSP parameters
! modified to accept both GTH and HGH pseudopotential types
  allocate(psppar(0:4,0:4,ntypes),nelpsp(ntypes),radii_cf(ntypes,2),npspcode(ntypes))
  do ityp=1,ntypes
     filename = 'psppar.'//atomnames(ityp)
     ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
     open(unit=11,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*) 'iproc=',iproc,': Failed to open the file (it must be in ABINIT format!) "',&
             trim(filename),'"'
        stop
     end if
     read(11,*)
     read(11,*) n_abinitzatom,nelpsp(ityp)
     read(11,*) npspcode(ityp)
     if (iproc.eq.0) write(*,'(a,1x,a,a,i3,a,i3)') 'atom type ',atomnames(ityp), & 
          ' is described by a ',nelpsp(ityp),' electron, with pspcode=',npspcode(ityp)
     psppar(:,:,ityp)=0.d0
     read(11,*) (psppar(0,j,ityp),j=0,4)
     if (npspcode(ityp) == 2) then !GTH case
        do i=1,2
           read(11,*) (psppar(i,j,ityp),j=0,3-i)
        enddo
        read(11,*) radii_cf(ityp,1),radii_cf(ityp,2)
     else if (npspcode(ityp) == 3) then !HGH case
        read(11,*) (psppar(1,j,ityp),j=0,3)
        do i=2,4
           read(11,*) (psppar(i,j,ityp),j=0,3)
           read(11,*) !k coefficients, not used (no spin-orbit coupling)
        enddo
        read(11,*) radii_cf(ityp,1),radii_cf(ityp,2)
     else
        stop 'unrecognized pspcode (accepts only GTH and HGH pseudopotentials in ABINIT format)'
     end if
     close(11)
  enddo

! Number of orbitals and their occupation number
  norb_vir=0

  nelec=0
  do iat=1,nat
     ityp=iatype(iat)
     nelec=nelec+nelpsp(ityp)
  enddo
  if (iproc.eq.0) then
     write(*,*) 'number of electrons',nelec
     if (mod(nelec,2).ne.0) write(*,*) 'WARNING: odd number of electrons, no closed shell system'
  end if
  norb=(nelec+1)/2+norb_vir

  allocate(occup(norb),eval(norb))

  nt=0
  do iorb=1,norb
     it=min(2,nelec-nt)
     occup(iorb)=it
     nt=nt+it
  enddo

  if (iproc.eq.0) then 
     write(*,*) 'number of orbitals',norb
     do iorb=1,norb
        write(*,*) 'occup(',iorb,')=',occup(iorb)
     enddo
  endif

! determine size alat of overall simulation cell
  call system_size(nat,rxyz,radii_cf(1,1),crmult,iatype,ntypes, &
       cxmin,cxmax,cymin,cymax,czmin,czmax)
  alat1=(cxmax-cxmin)
  alat2=(cymax-cymin)
  alat3=(czmax-czmin)

! shift atomic positions such that molecule is inside cell
  if (iproc.eq.0) write(*,'(a,3(1x,e14.7))') 'atomic positions shifted by',-cxmin,-cymin,-czmin
  do iat=1,nat
     rxyz(1,iat)=rxyz(1,iat)-cxmin
     rxyz(2,iat)=rxyz(2,iat)-cymin
     rxyz(3,iat)=rxyz(3,iat)-czmin
  enddo

  if (iproc.eq.0) then
     write(*,*) 'Shifted atomic positions'
     do iat=1,nat
        write(*,'(i3,3(1x,e14.7))') iat,(rxyz(j,iat),j=1,3)
     enddo
  endif

! grid sizes n1,n2,n3
  n1=int(alat1/hgrid)
  n2=int(alat2/hgrid)
  n3=int(alat3/hgrid)
  alat1=n1*hgrid ; alat2=n2*hgrid ; alat3=n3*hgrid
  if (iproc.eq.0) then 
     write(*,*) 'n1,n2,n3',n1,n2,n3
     write(*,*) 'total number of grid points',(n1+1)*(n2+1)*(n3+1)
     write(*,'(a,3(1x,e12.5))') 'simulation cell',alat1,alat2,alat3
     !           write(*,'(a,3(1x,e24.17))') 'simulation cell',alat1,alat2,alat3
  endif

! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 ; nfl2=n2 ; nfl3=n3
  nfu1=0 ; nfu2=0 ; nfu3=0
  do iat=1,nat
     rad=radii_cf(iatype(iat),2)*frmult
     nfl1=min(nfl1,int(onem+(rxyz(1,iat)-rad)/hgrid)) ; nfu1=max(nfu1,int((rxyz(1,iat)+rad)/hgrid))
     nfl2=min(nfl2,int(onem+(rxyz(2,iat)-rad)/hgrid)) ; nfu2=max(nfu2,int((rxyz(2,iat)+rad)/hgrid))
     nfl3=min(nfl3,int(onem+(rxyz(3,iat)-rad)/hgrid)) ; nfu3=max(nfu3,int((rxyz(3,iat)+rad)/hgrid))
  enddo
  if (iproc.eq.0) then
     write(*,*) 'nfl1,nfu1 ',nfl1,nfu1
     write(*,*) 'nfl2,nfu2 ',nfl2,nfu2
     write(*,*) 'nfl3,nfu3 ',nfl3,nfu3
  endif

! Create wavefunctions descriptors and allocate them
  allocate(ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2))
  allocate(ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2))

! Create the file grid.ascii to visualize the grid of functions
  if (iproc.eq.0 .and. output_grid) then
     open(unit=22,file='grid.ascii',status='unknown')
     write(22,*) nat
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     do iat=1,nat
        write(22,'(3(1x,e12.5),3x,a20)') rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
     enddo
  endif

  call createWavefunctionsDescriptors(parallel, iproc, nproc, idsx, n1, n2, n3, output_grid, hgrid, &
       & nat, ntypes, iatype, atomnames, rxyz, radii_cf, crmult, frmult, ibyz_c,ibxz_c,ibxy_c, &
       & ibyz_f, ibxz_f, ibxy_f, nseg_c, nseg_f, nvctr_c, nvctr_f, nvctrp, keyg,  keyv, &
       & norb,norbp,psi,hpsi,psit,psidst,hpsidst,ads)

! Calculate all projectors
  allocate(nseg_p(0:2*nat))
  allocate(nvctr_p(0:2*nat))
  allocate(nboxp_c(2,3,nat),nboxp_f(2,3,nat))
  call createProjectorsArrays(iproc, n1, n2, n3, rxyz, nat, ntypes, iatype, atomnames, &
       & psppar, npspcode, radii_cf, cpmult, fpmult, hgrid, nvctr_p, nseg_p, &
       & keyg_p, keyv_p, nproj, nprojel, istart, nboxp_c, nboxp_f, proj)


  if (iproc.eq.0) write(*,*) 'Size of real space grids',(2*n1+31),(2*n2+31),(2*n3+31)
! Charge density, Potential in real space
  if (iproc.eq.0) write(79,'(a40,i10)') 'words for rhopot and pot_ion ',2*(2*n1+31)*(2*n2+31)*(2*n3+31)
  allocate(rhopot((2*n1+31),(2*n2+31),(2*n3+31)),pot_ion((2*n1+31)*(2*n2+31)*(2*n3+31)))
  call razero((2*n1+31)*(2*n2+31)*(2*n3+31),pot_ion)
  if (iproc.eq.0) write(79,*) 'allocation done'
! Allocate and calculate the 1/|r-r'| kernel for the solution of Poisson's equation and test it
  ndegree_ip=14

  call createKernel(parallel, nfft1, nfft2, nfft3, n1, n2, n3, hgridh, &
       & ndegree_ip, iproc, nproc, pkernel)
  !call PARtest_kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,pot_ion,rhopot,iproc,nproc)

! Precalculate ionic potential from PSP charge densities and local Gaussian terms
  call input_rho_ion(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,n1,n2,n3,hgrid,pot_ion,eion)
  if (iproc.eq.0) write(*,*) 'ion-ion interaction energy',eion
  if (parallel) then
     call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,0,  &
          rhopot,pot_ion,ehart,eexcu,vexcu,iproc,nproc)
  else
     call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,0,  &
          rhopot,pot_ion,ehart,eexcu,vexcu)
  end if

  call addlocgauspsp(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,n1,n2,n3,hgrid,pot_ion)

     ! INPUT WAVEFUNCTIONS
  if (inputPsiId == 0) then
     call input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          nat,norb,norbp,n1,n2,n3,nfft1,nfft2,nfft3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
          rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
          atomnames,ntypes,iatype,pkernel,psppar,npspcode,psi,eval,accurex)
     if (iproc.eq.0) then
        write(*,*) 'expected accuracy in total energy due to grid size',accurex
        write(*,*) 'suggested value for gnrm_cv ',accurex
     endif
     if (iproc.eq.0) write(*,*) 'input wavefunction has been calculated'

  else if (inputPsiId == 1 ) then
     if (iproc.eq.0) write(*,*) 'START reformatting psi from old psi'
     call reformatmywaves(iproc, norb, norbp, nat, &
          & hgrid_old, nvctr_c_old, nvctr_f_old, n1_old, n2_old, n3_old, rxyz_old, &
          & nseg_c_old, nseg_f_old, keyg_old, keyv_old, psi_old, &
          & hgrid, nvctr_c, nvctr_f, n1, n2, n3, rxyz, &
          & nseg_c, nseg_f, keyg, keyv, psi)
            eval=eval_old
     deallocate(keyg_old, keyv_old)
     deallocate(psi_old,eval_old)

  else if (inputPsiId == 2) then
     call readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
  end if

  if (parallel) then
     call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
     call orthon_p(iproc,nproc,norb,norbp,nvctrp,psit)
     !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
     call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
  else
     call orthon(norb,norbp,nvctrp,psi)
     call checkortho(norb,norbp,nvctrp,psi)
  endif

  alpha=1.d0
  energy=1.d100
  gnrm=1.d100
  ekin_sum=0.d0 ; epot_sum=0.d0 ; eproj_sum=0.d0
! loop for wavefunction minimization
  do 1000, iter=1,itermax
     if (idsx.gt.0) mids=mod(iter-1,idsx)+1
     if (iproc.eq.0) then 
        write(*,*) '-------------------------------------- iter= ',iter
        if (gnrm.le.gnrm_cv) then
           write(*,'(a,i3,3(1x,1pe18.11))') 'iproc,ehart,eexcu,vexcu',iproc,ehart,eexcu,vexcu
           write(*,'(a,3(1x,1pe18.11))') 'final ekin_sum,epot_sum,eproj_sum',ekin_sum,epot_sum,eproj_sum
           write(*,'(a,3(1x,1pe18.11))') 'final ehart,eexcu,vexcu',ehart,eexcu,vexcu
           write(*,'(a,i6,2x,1pe19.12,1x,1pe9.2)') 'FINAL iter,total energy,gnrm',iter,energy,gnrm
        endif
     endif

! Potential from electronic charge density
     call sumrho(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)

!     ixc=12   ! PBE functional
     ixc=1   ! LDA functional
     if (parallel) then
        call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,ixc, & 
             pot_ion,rhopot,ehart,eexcu,vexcu,iproc,nproc)
     else
        call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,ixc, & 
             pot_ion,rhopot,ehart,eexcu,vexcu)
     end if

! local potential and kinetic energy for all orbitals belonging to iproc
     call applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          hgrid,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
          ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          psi,rhopot,hpsi,epot_sum,ekin_sum)

! apply all PSP projectors for all orbitals belonging to iproc
     call applyprojectorsall(iproc,ntypes,nat,iatype,psppar,npspcode,occup, &
          nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
          norb,norbp,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)

     if (parallel) then
        wrkallred(1,2)=ekin_sum ; wrkallred(2,2)=epot_sum ; wrkallred(3,2)=eproj_sum
        call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        ekin_sum=wrkallred(1,1) ; epot_sum=wrkallred(2,1) ; eproj_sum=wrkallred(3,1) 
     endif
     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-ehart+eexcu-vexcu+eion

     !check for convergence or whether max. numb. of iterations exceeded
     if (gnrm.le.gnrm_cv .or. iter.eq.itermax) then 
        if (iproc.eq.0) then 
           write(*,*) iter,' minimization iterations required'
        end if
        goto 1010
     endif


! Apply  orthogonality constraints to all orbitals belonging to iproc
     if (parallel) then
        allocate(hpsit(nvctrp,norbp*nproc))
        call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
        call  orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
        call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsit,hpsi)
        deallocate(hpsit)
     else
        call  orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
     endif
     
! norm of gradient
     gnrm=0.d0
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        scpr=dnrm2(nvctr_c+7*nvctr_f,hpsi(1,iorb-iproc*norbp),1) 
        gnrm=gnrm+scpr**2
     enddo
     if (parallel) then
        tt=gnrm
        call MPI_ALLREDUCE(tt,gnrm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     endif
     gnrm=sqrt(gnrm/norb)

! Preconditions all orbitals belonging to iproc
     call preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid, &
          ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)

     !       call plot_wf(10,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f, &
     !                       rxyz(1,1),rxyz(2,1),rxyz(3,1),psi)
     !       call plot_wf(20,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f, &
     !                       rxyz(1,1),rxyz(2,1),rxyz(3,1),hpsi)


     if (parallel) then
        if (idsx.gt.0) then
           call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsidst(:,:,mids))
           do iorb=1,norb
              do k=1,nvctrp
                 psidst(k,iorb,mids)= psit(k,iorb) 
              enddo
           enddo
        else
           allocate(hpsit(nvctrp,norbp*nproc))
           call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
        endif
     else
        if (idsx.gt.0) then
           do iorb=1,norb
              do k=1,nvctrp
                 psidst(k,iorb,mids)= psi(k,iorb)
                 hpsidst(k,iorb,mids)=hpsi(k,iorb)
              enddo
           enddo
        endif
     endif


     if (idsx.gt.0) then 

        if (parallel) then
           call diisstp(parallel,norb,norbp,nproc,iproc,  &
                ads,iter,mids,idsx,nvctrp,psit,psidst,hpsidst)
        else
           call diisstp(parallel,norb,norbp,nproc,iproc,  &
                ads,iter,mids,idsx,nvctrp,psi,psidst,hpsidst)
        endif

     else

! update all wavefunctions with the preconditioned gradient
        if (energy.gt.energy_old) then
           alpha=max(.125d0,.5d0*alpha)
           if (alpha.eq..125d0) write(*,*) 'Convergence problem or limit'
        else
           alpha=min(1.05d0*alpha,1.d0)
        endif
        if (iproc.eq.0) write(*,*) 'alpha=',alpha

        if (parallel) then
           do iorb=1,norb
              call DAXPY(nvctrp,-alpha,hpsit(1,iorb),1,psit(1,iorb),1)
           enddo
           deallocate(hpsit)
        else
           do iorb=1,norb
              call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psi(1,iorb),1)
           enddo
        endif

     endif
     
     if (parallel) then
        call orthon_p(iproc,nproc,norb,norbp,nvctrp,psit)
        !       call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
     else
        call orthon(norb,norbp,nvctrp,psi)
        !       call checkortho(norb,norbp,nvctrp,psi)
     endif


     tt=energybs-scprsum
     if (abs(tt).gt.1.d-8) then 
        write(*,*) 'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
     endif
     if (iproc.eq.0) then
        write(*,'(a,3(1x,e18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             ekin_sum,epot_sum,eproj_sum
        write(*,'(a,3(1x,e18.11))') 'ehart,eexcu,vexcu',ehart,eexcu,vexcu
        write(*,'(a,i6,2x,e19.12,1x,e9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
     endif

     if (parallel) then
        call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
     endif

1000 continue
  write(*,*) 'No convergence within the allowed number of minimization steps'
1010 continue
  if (idsx.gt.0) then
     deallocate(psidst,hpsidst,ads)
  end if

!------------------------------------------------------------------------
! transform to KS orbitals
  if (parallel) then
     allocate(hpsit(nvctrp,norbp*nproc))
     call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
     call KStrans_p(iproc,nproc,norb,norbp,nvctrp,occup,hpsit,psit,evsum,eval)
     deallocate(hpsit)
     call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
  else
     call KStrans(norb,norbp,nvctrp,occup,hpsi,psi,evsum,eval)
  endif
  deallocate(hpsi)
  if (parallel) deallocate(psit)
  if (abs(evsum-energybs).gt.1.d-8) write(*,*) 'Difference:evsum,energybs',evsum,energybs

!  write all the wavefunctions into files
  if (output_wf) then
     call  writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
              nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
     write(*,*) iproc,' finished writing waves'
  end if


!------------------------------------------------------------------------
! here we start the calculation of the forces
  if (iproc.eq.0) write(*,*)'calculation of forces'

! Selfconsistent potential is saved in rhopot, new arrays rho, hartpot for calculation of forces
! ground state electronic density
  deallocate(pot_ion)
  allocate(rho((2*n1+31)*(2*n2+31)*(2*n3+31)))
  call sumrho(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  &
              nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho)


  allocate(hartpot((2*n1+31)*(2*n2+31)*(2*n3+31)))
  call DCOPY((2*n1+31)*(2*n2+31)*(2*n3+31),rho,1,hartpot,1) 
  if (parallel) then
     call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,0, &
          hartpot,hartpot,ehart_fake,eexcu_fake,vexcu_fake,iproc,nproc)
  else
     call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,0, &
          hartpot,hartpot,ehart_fake,eexcu_fake,vexcu_fake)
  end if
  deallocate(pkernel)

  if (iproc.eq.0) write(*,*)'electronic potential calculated'
  allocate(gxyz(3,nat))

        call timing(iproc,'Forces        ','ON')
! calculate local part of the forces gxyz
   call local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
                     n1,n2,n3,rho,hartpot,gxyz)
   deallocate(rho,hartpot)

! Add the nonlocal part of the forces to gxyz
! calculating derivatives of the projectors (for the moment recalculate projectors)
  call nonlocal_forces(iproc,nproc,n1,n2,n3,nboxp_c,nboxp_f, &
     ntypes,nat,norb,norbp,istart,nprojel,nproj,&
     iatype,psppar,npspcode,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,  &
     keyg,keyv,keyg_p,keyv_p,psi,rxyz,radii_cf,cpmult,fpmult,hgrid,gxyz)

! Add up all the force contributions
  if (parallel) then
     call MPI_ALLREDUCE(gxyz,fxyz,3*nat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     do iat=1,nat
        fxyz(1,iat)=gxyz(1,iat)
        fxyz(2,iat)=gxyz(2,iat)
        fxyz(3,iat)=gxyz(3,iat)
     enddo
  end if

  deallocate(gxyz)

  call timing(iproc,'Forces        ','OF')

  !------------------------------------------------------------------------
  if (calc_tail) then
     call timing(iproc,'Tail          ','ON')
!    Calculate kinetic energy correction due to boundary conditions
     nbuf=nint(rbuf/hgrid)
     if (iproc.eq.0) write(*,*) 'tail requires ',nbuf,' additional grid points around cell'
!    --- new grid sizes n1,n2,n3
     nb1=n1+2*nbuf
     nb2=n2+2*nbuf
     nb3=n3+2*nbuf
     alatb1=nb1*hgrid ; alatb2=nb2*hgrid ; alatb3=nb3*hgrid
     if (iproc.eq.0) then 
        write(*,*) 'BIG: n1,n2,n3',nb1,nb2,nb3
        write(*,*) 'BIG: total number of grid points',(nb1+1)*(nb2+1)*(nb3+1)
        write(*,'(a,3(1x,e12.5))') 'BIG: simulation cell',alatb1,alatb2,alatb3
     endif
!    ---reformat potential
     allocate(rhopotb((2*nb1+31),(2*nb2+31),(2*nb3+31)))
     call razero((2*nb1+31)*(2*nb2+31)*(2*nb3+31),rhopotb)
     do i3=1+2*nbuf,2*n3+31+2*nbuf
     do i2=1+2*nbuf,2*n2+31+2*nbuf
     do i1=1+2*nbuf,2*n1+31+2*nbuf
        rhopotb(i1,i2,i3)=rhopot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf)
     enddo ; enddo ; enddo
     deallocate(rhopot)

!    ---reformat keyg_p
     do iseg=1,nseg_p(2*nat)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/((n1+1)*(n2+1))
        ii=ii-i3*(n1+1)*(n2+1)
        i2=ii/(n1+1)
        i0=ii-i2*(n1+1)
        i1=i0+j1-j0
        i3=i3+nbuf
        i2=i2+nbuf
        i1=i1+nbuf
        i0=i0+nbuf
        j0=i3*((nb1+1)*(nb2+1)) + i2*(nb1+1) + i0+1
        j1=i3*((nb1+1)*(nb2+1)) + i2*(nb1+1) + i1+1
        keyg_p(1,iseg)=j0
        keyg_p(2,iseg)=j1
     enddo

!    ---reformat wavefunctions

!    fine grid size (needed for creation of input wavefunction, preconditioning)
     nbfl1=nfl1+nbuf ; nbfl2=nfl2+nbuf ; nbfl3=nfl3+nbuf
     nbfu1=nfu1+nbuf ; nbfu2=nfu2+nbuf ; nbfu3=nfu3+nbuf
     if (iproc.eq.0) then
        write(*,*) 'BIG:nfl1,nfu1 ',nbfl1,nbfu1
        write(*,*) 'BIG:nfl2,nfu2 ',nbfl2,nbfu2
        write(*,*) 'BIG:nfl3,nfu3 ',nbfl3,nbfu3
     endif

     allocate(txyz(3,nat))
     do iat=1,nat
        txyz(1,iat)=rxyz(1,iat)+nbuf*hgrid
        txyz(2,iat)=rxyz(2,iat)+nbuf*hgrid
        txyz(3,iat)=rxyz(3,iat)+nbuf*hgrid
     enddo

!    determine localization region for all orbitals, but do not yet fill the descriptor arrays
     allocate(logrid_c(0:nb1,0:nb2,0:nb3),logrid_f(0:nb1,0:nb2,0:nb3))
     allocate(ibbyz_c(2,0:nb2,0:nb3),ibbxz_c(2,0:nb1,0:nb3),ibbxy_c(2,0:nb1,0:nb2))
     allocate(ibbyz_f(2,0:nb2,0:nb3),ibbxz_f(2,0:nb1,0:nb3),ibbxy_f(2,0:nb1,0:nb2))

!    coarse grid quantities
     call fill_logrid(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,nbuf,nat,ntypes,iatype,txyz, & 
                     radii_cf(1,1),crmult,hgrid,logrid_c)
     if (iproc.eq.0 .and. output_grid) then
        open(unit=22,file='grid.ascii',status='unknown')
        write(22,*) nat
        write(22,*) alat1,' 0. ',alat2
        write(22,*) ' 0. ',' 0. ',alat3
        do iat=1,nat
           write(22,'(3(1x,e12.5),3x,a20)') txyz(1,iat),txyz(2,iat),txyz(3,iat),atomnames(iatype(iat))
        enddo
        do i3=0,nb3 ; do i2=0,nb2 ; do i1=0,nb1
           if (logrid_c(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
        enddo ; enddo ; enddo 
     endif
     call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,nvctrb_c)
     if (iproc.eq.0) write(*,*) 'BIG:orbitals have coarse segment, elements',nsegb_c,nvctrb_c
     call bounds(nb1,nb2,nb3,logrid_c,ibbyz_c,ibbxz_c,ibbxy_c)

!    fine grid quantities
     call fill_logrid(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,0,nat,ntypes,iatype,txyz, & 
                     radii_cf(1,2),frmult,hgrid,logrid_f)
     if (iproc.eq.0 .and. output_grid) then
        do i3=0,nb3 ; do i2=0,nb2 ; do i1=0,nb1
           if (logrid_f(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
          enddo ; enddo ; enddo 
         endif
	 call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,nvctrb_f)
        if (iproc.eq.0) write(*,*) 'BIG:orbitals have fine   segment, elements',nsegb_f,7*nvctrb_f
        call bounds(nb1,nb2,nb3,logrid_f,ibbyz_f,ibbxz_f,ibbxy_f)

        if (iproc.eq.0) close(22)

! now fill the wavefunction descriptor arrays
        allocate(keybg(2,nsegb_c+nsegb_f),keybv(nsegb_c+nsegb_f))
! coarse grid quantities
        call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,keybg(1,1),keybv(1))

! fine grid quantities
        call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,keybg(1,nsegb_c+1),keybv(nsegb_c+1))

        deallocate(logrid_c,logrid_f)
! allocations for arrays holding the wavefunction
        if (iproc.eq.0) write(79,'(a40,i10)') 'words for psib and hpsib ',2*(nvctrb_c+7*nvctrb_f)
        allocate(psib(nvctrb_c+7*nvctrb_f),hpsib(nvctrb_c+7*nvctrb_f))
        if (iproc.eq.0) write(79,*) 'allocation done'

! work arrays applylocpotkin
        allocate(psig(8*(nb1+1)*(nb2+1)*(nb3+1)) )
        allocate(psigp(8*(nb1+1)*(nb2+1)*(nb3+1)) )
        allocate(psifscf((2*nb1+31)*(2*nb2+31)*(2*nb3+16)) )
        allocate(psir((2*nb1+31)*(2*nb2+31)*(2*nb3+31)))

      ekin_sum=0.d0
      epot_sum=0.d0
      eproj_sum=0.d0
    do 2500 iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call uncompress(n1,n2,n3,-nbuf,n1+nbuf,-nbuf,n2+nbuf,-nbuf,n3+nbuf, & 
                    nseg_c,nvctr_c,keyg,              keyv,   &
                    nseg_f,nvctr_f,keyg(:,nseg_c+1),keyv(nseg_c+1:),   &
                    psi(:,iorb-iproc*norbp),psi(nvctr_c+1:,iorb-iproc*norbp),psig)
        call compress(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,  &
                    nsegb_c,nvctrb_c,keybg(1,1),       keybv(1),   &
                    nsegb_f,nvctrb_f,keybg(1,nsegb_c+1),keybv(nsegb_c+1),   &
                    psig,psib(1),psib(nvctrb_c+1))

        npt=2
        do 1500 ipt=1,npt

!          calculate gradient
           call applylocpotkinone(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3, & 
                   hgrid,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f,keybg,keybv,  & 
                   ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f, & 
                   psig,psigp,psifscf,psir,  &
                   psib,rhopotb,hpsib,epot,ekin)
           call applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    nsegb_c,nsegb_f,keybg,keybv,nvctrb_c,nvctrb_f,  & 
                    psib,hpsib,eproj)
           tt=0.d0
           do i=1,nvctrb_c+7*nvctrb_f
              hpsib(i)=hpsib(i)-eval(iorb)*psib(i)
              tt=tt+hpsib(i)**2
           enddo
           tt=sqrt(tt)
           write(*,'(a,i3,3(1x,e21.14),1x,e10.3)') 'BIG:iorb,ekin,epot,eproj,gnrm',iorb,ekin,epot,eproj,tt
           if (ipt.eq.npt) goto 1600
!          calculate tail
           cprecr=-eval(iorb)
           call timing(iproc,'Tail          ','OF')
           call precong(iorb,nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3, &
                     nsegb_c,nvctrb_c,nsegb_f,nvctrb_f,keybg,keybv, &
                     ncongt,cprecr,hgrid,ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,hpsib)
           call timing(iproc,'Tail          ','ON')
!          call plot_wf(10,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keybg,keybv,nsegb_f,nvctrb_f,  & 
!                    txyz(1,1),txyz(2,1),txyz(3,1),psib)
!          add tail
           sum=0.d0
           do i=1,nvctrb_c+7*nvctrb_f
              psib(i)=psib(i)-hpsib(i)
              sum=sum+psib(i)**2
           enddo
           sum=sqrt(sum)
           write(*,*) 'norm orbital + tail',iorb,sum
!          call plot_wf(20,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keybg,keybv,nsegb_f,nvctrb_f,  & 
!                    txyz(1,1),txyz(2,1),txyz(3,1),psib)

           sum=1.d0/sum
           do i=1,nvctrb_c+7*nvctrb_f
              psib(i)=psib(i)*sum
           enddo

1500    continue
1600    continue

        ekin_sum=ekin_sum+ekin*occup(iorb)
        epot_sum=epot_sum+epot*occup(iorb)
        eproj_sum=eproj_sum+eproj*occup(iorb)


2500 continue
     if (parallel) then
        wrkallred(1,2)=ekin_sum ; wrkallred(2,2)=epot_sum ; wrkallred(3,2)=eproj_sum 
        call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        ekin_sum=wrkallred(1,1) ; epot_sum=wrkallred(2,1) ; eproj_sum=wrkallred(3,1)  
     endif

     energybs=ekin_sum+epot_sum+eproj_sum
     energy=energybs-ehart+eexcu-vexcu+eion

     if (iproc.eq.0) then
        write(*,*) 'total energy with tail correction',energy
        write(*,*) 'ekin,epot,eproj with tail correction',ekin_sum,epot_sum,eproj_sum
     endif

     deallocate(txyz)
     deallocate(rhopotb,psig,psigp,psifscf,psir)

     deallocate(psib,hpsib)
     deallocate(keybg,keybv)
     deallocate(ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f)
     call timing(iproc,'Tail          ','OF')
  else
!    No tail calculation
     deallocate(rhopot)
  endif
  ! --- End if of tail calculation

  deallocate(ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f)

  deallocate(keyg_p,keyv_p,proj)
  deallocate(occup)
  deallocate(nvctr_p,nseg_p)
  deallocate(psppar,nelpsp,radii_cf,npspcode)

  call timing(iproc,'              ','RE')
  call cpu_time(tcpu1)
  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  write(*,'(a,1x,i4,2(1x,f12.2))') 'iproc, elapsed, CPU time ', iproc,tel,tcpu1-tcpu0

        if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE cluster


subroutine transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
   implicit real*8 (a-h,o-z)
   logical, parameter :: parallel=.true.
   integer recvcount,sendcount
   dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
   real*8, allocatable :: psiw(:,:,:)
   include 'mpif.h'
 
   call timing(iproc,'Un-Transall   ','ON')
   allocate(psiw(nvctrp,norbp,nproc))
 
   sendcount=nvctrp*norbp
   recvcount=nvctrp*norbp

!  reformatting: psiw(i,iorb,j,jorb) <- psi(ij,iorb,jorb)
   do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
      ij=1
      do j=1,nproc
         do i=1,nvctrp
            if (ij .le. nvctr_c+7*nvctr_f) then
               psiw(i,iorb-iproc*norbp,j)=psi(ij,iorb-iproc*norbp)
            else
               psiw(i,iorb-iproc*norbp,j)=0.d0
            endif
            ij=ij+1
         enddo
      enddo
   enddo

!  transposition: psit(i,iorb,jorb,j) <- psiw(i,iorb,j,jorb) 
   call MPI_ALLTOALL(psiw,sendcount,MPI_DOUBLE_PRECISION,  &
                     psit,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
 
   call timing(iproc,'Un-Transall   ','OF')

   deallocate(psiw)

END SUBROUTINE transallwaves


subroutine untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
   implicit real*8 (a-h,o-z)
   logical, parameter :: parallel=.true.
   integer recvcount,sendcount
   dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
   real*8, allocatable :: psiw(:,:,:)
   include 'mpif.h'
 
   call timing(iproc,'Un-Transall   ','ON')

   allocate(psiw(nvctrp,norbp,nproc))

   sendcount=nvctrp*norbp
   recvcount=nvctrp*norbp
 
!  transposition: psiw(i,iorb,j,jorb) <- psit(i,iorb,jorb,j) 
   call MPI_ALLTOALL(psit,sendcount,MPI_DOUBLE_PRECISION,  &
                     psiw,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
 
!  reformatting: psi(ij,iorb,jorb) <- psiw(i,iorb,j,jorb)
   loop_iorb: do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
      ij=1
      do j=1,nproc
         do i=1,nvctrp
            psi(ij,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
            ij=ij+1
            if (ij.gt. nvctr_c+7*nvctr_f) cycle loop_iorb
         end do
      end do
   end do loop_iorb
 
   deallocate(psiw)
   call timing(iproc,'Un-Transall   ','OF')
 
END SUBROUTINE untransallwaves

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
        
subroutine input_rho_ion(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar, &
     & nelpsp,n1,n2,n3,hgrid,rho,eion)
  ! Creates charge density arising from the ionoc PSP cores
  implicit real*8 (a-h,o-z)
  character*20 :: atomnames(100)
  dimension psppar(0:4,0:4,ntypes),rxyz(3,nat),iatype(nat),nelpsp(ntypes)
  dimension rho(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

  call timing(iproc,'CrtLocPot     ','ON')

  hgridh=hgrid*.5d0 
  pi=4.d0*atan(1.d0)
  call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

  ! Ionic charge 
  rholeaked=0.d0
  eion=0.d0
  do iat=1,nat
     ityp=iatype(iat)
     rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
     ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
     !    ion-ion interaction
     do jat=1,iat-1
        dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
        jtyp=iatype(jat)
        eion=eion+nelpsp(jtyp)*nelpsp(ityp)/dist
     enddo

     rloc=psppar(0,0,ityp)
!     if (iproc.eq.0) write(*,'(a,i3,a,a,a,e10.3)') 'atom ',iat,' of type ',atomnames(ityp),' has an ionic charge whose rloc is ',rloc
     charge=nelpsp(ityp)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

     do i3=iz-ii,iz+ii
        do i2=iy-ii,iy+ii
           do i1=ix-ii,ix+ii
              x=i1*hgridh-rx
              y=i2*hgridh-ry
              z=i3*hgridh-rz
              r2=x**2+y**2+z**2
              arg=r2/rloc**2
              xp=exp(-.5d0*arg)
              if (i3.ge.-14 .and. i3.le.2*n3+16  .and.  & 
                   i2.ge.-14 .and. i2.le.2*n2+16  .and.  & 
                   i1.ge.-14 .and. i1.le.2*n1+16 ) then
                 rho(i1,i2,i3)=rho(i1,i2,i3)-xp*charge
              else
                 rholeaked=rholeaked+xp*charge
              endif
           enddo
        enddo
     enddo

  enddo
  
  ! Check
  tt=0.d0
  do i3= -14,2*n3+16
     do i2= -14,2*n2+16
        do i1= -14,2*n1+16
           tt=tt+rho(i1,i2,i3)
        enddo
     enddo
  enddo
  tt=tt*hgridh**3
  rholeaked=rholeaked*hgridh**3
  if (iproc.eq.0) write(*,'(a,e21.14,1x,e10.3)') 'total ionic charge,leaked charge: ',tt,rholeaked

  call timing(iproc,'CrtLocPot     ','OF')

  return
END SUBROUTINE 


        subroutine addlocgauspsp(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,n1,n2,n3,hgrid,pot)
! Add local Gaussian terms of the PSP to pot 
        implicit real*8 (a-h,o-z)
        character*20 :: atomnames(100)
        dimension psppar(0:4,0:4,ntypes),rxyz(3,nat),iatype(nat)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 

   do iat=1,nat
   ityp=iatype(iat)
   rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
   ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
! determine number of local terms
     nloc=0
     do iloc=1,4
     if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo

     rloc=psppar(0,0,ityp)
!     if (iproc.eq.0) write(*,'(a,i4,a,a,a,i3,a,1x,e9.2)')  & 
!     'atom ',iat,' is of type ',atomnames(ityp),' and has ',nloc,' local PSP terms with rloc',rloc
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

      do i3=max(-14,iz-ii),min(2*n3+16,iz+ii)
      do i2=max(-14,iy-ii),min(2*n2+16,iy+ii)
      do i1=max(-14,ix-ii),min(2*n1+16,ix+ii)
         x=i1*hgridh-rx
         y=i2*hgridh-ry
         z=i3*hgridh-rz
         r2=x**2+y**2+z**2
         arg=r2/rloc**2
        xp=exp(-.5d0*arg)
        tt=psppar(0,nloc,ityp)
        do iloc=nloc-1,1,-1
        tt=arg*tt+psppar(0,iloc,ityp)
        enddo
        pot(i1,i2,i3)=pot(i1,i2,i3)+xp*tt
      enddo
      enddo
      enddo

!! For testing only: Add erf part (in the final version that should be part of Hartree pot)
!      do i3=-14,2*n3+16
!      do i2=-14,2*n2+16
!      do i1=-14,2*n1+16
!         x=(i1-ix)*hgridh
!         y=(i2-iy)*hgridh
!         z=(i3-iz)*hgridh
!         r2=x**2+y**2+z**2
!         r=sqrt(r2)
!         arg=r*(sqrt(.5d0)/rloc)
!         if (arg.lt.1.d-7) then 
!! Taylor expansion
!         x=arg**2
!         tt=   -0.37612638903183752463d0*x + 1.1283791670955125739d0
!         tt=tt*(sqrt(.5d0)/rloc)
!         else
!          tt=derf(arg)/r
!         endif
!        pot(i1,i2,i3)=pot(i1,i2,i3)+nelpsp(ityp)*tt
!      enddo
!      enddo
!      enddo


    enddo

        return
	END SUBROUTINE



        subroutine applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                   hgrid,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
                   psi,pot,hpsi,epot_sum,ekin_sum)
!  Applies the local potential and kinetic energy operator to all wavefunctions belonging to processor
! Input: pot,psi
! Output: hpsi,epot,ekin
        implicit real*8 (a-h,o-z)
        dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
        dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
        dimension occup(norb),pot((2*n1+31)*(2*n2+31)*(2*n3+31))
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension  psi(nvctr_c+7*nvctr_f,norbp)
        dimension hpsi(nvctr_c+7*nvctr_f,norbp)
        real*8, allocatable, dimension(:) :: psifscf,psir,psig,psigp

      call timing(iproc,'ApplyLocPotKin','ON')

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate( psig(8*(n1+1)*(n2+1)*(n3+1)) )
        allocate(psigp(8*(n1+1)*(n2+1)*(n3+1)) )
        allocate(psifscf(max( (2*n1+31)*(2*n2+31)*(2*n3+16),&
                        &     (2*n1+16)*(2*n2+31)*(2*n3+31))) )
! Wavefunction in real space
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))

        ekin_sum=0.d0
        epot_sum=0.d0
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                   hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
                   ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
                   psig,psigp,psifscf,psir,  &
                   psi(1,iorb-iproc*norbp),pot,hpsi(1,iorb-iproc*norbp),epot,ekin)
        ekin_sum=ekin_sum+occup(iorb)*ekin
        epot_sum=epot_sum+occup(iorb)*epot
!        write(*,'(a,i5,2(1x,e17.10))') 'iorb,ekin,epot',iorb,ekin,epot

     enddo


        deallocate(psig,psigp,psifscf,psir)

      call timing(iproc,'ApplyLocPotKin','OF')

        END SUBROUTINE


        subroutine applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                   hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
                   ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
                   psig,psigp,psifscf,psir,  &
                   psi,pot,hpsi,epot,ekin)
!  Applies the local potential and kinetic energy operator to one wavefunction 
! Input: pot,psi
! Output: hpsi,epot,ekin
        implicit real*8 (a-h,o-z)
        dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
        dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
        dimension pot((2*n1+31)*(2*n2+31)*(2*n3+31))
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension  psi(nvctr_c+7*nvctr_f)
        dimension hpsi(nvctr_c+7*nvctr_f)
        dimension   psig(8*(n1+1)*(n2+1)*(n3+1))
        dimension psigp(8*(n1+1)*(n2+1)*(n3+1)) 
        dimension psifscf(max( (2*n1+31)*(2*n2+31)*(2*n3+16), &
                        &      (2*n1+16)*(2*n2+31)*(2*n3+31)) ) 
        dimension psir((2*n1+31)*(2*n2+31)*(2*n3+31))


        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1),psi(nvctr_c+1),psig)

        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n1+2))

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=psifscf and ww(((2*n1+31)*(2*n2+31)*(2*n1+16))
        epot=0.d0
        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
          tt=pot(i)*psir(i)
          epot=epot+tt*psir(i)
          psir(i)=tt
        enddo

        call convolut_magic_t(2*n1+15,2*n2+15,2*n3+15,psir,psifscf) !psifscf=ww((2*n1+16)*(2*n2+31)*(2*n1+31))

        call analyse_shrink(n1,n2,n3,psir,psifscf,psigp)  !psir=ww(((2*n1+31)*(2*n2+31)*(2*n1+16))

         call ConvolkineticP(n1,n2,n3,  &
              nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
              psig,psigp,ekin)

        call compress(n1,n2,n3,0,n1,0,n2,0,n3,  &
                    nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psigp,hpsi(1),hpsi(nvctr_c+1))

        return
    END SUBROUTINE

    
        subroutine uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                              mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                              mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                              psi_c,psi_f,psig)
! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)

        call razero(8*(nu1-nl1+1)*(nu2-nl2+1)*(nu3-nl3+1),psig)

! coarse part
        do iseg=1,mseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
          ii=ii+1
            psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
          enddo
         enddo

! fine part
        do iseg=1,mseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
            psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
            psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
            psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
            psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
            psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
            psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
            psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
          enddo
         enddo

        END SUBROUTINE

    
        subroutine compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                            mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                            mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                            psig,psi_c,psi_f)
! Compresses a psig wavefunction into psi_c,psi_f form
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)
        
! coarse part
        do iseg=1,mseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
            psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
          enddo
        enddo

! fine part
        do iseg=1,mseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
            psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
            psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
            psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
            psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
            psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
            psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
            psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
          enddo
        enddo

        END SUBROUTINE



        subroutine convolut_magic_n(n1,n2,n3,x,y)
! Applies the magic filter matrix ( no transposition) ; data set grows
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
        dimension x(0:n1,0:n2,0:n3)
        dimension y(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil)
 
!  (i1,i2*i3) -> (i2*i3,I1)
        ndat=(n2+1)*(n3+1)
        call convrot_grow(n1,ndat,x,y)
!  (i2,i3*I1) -> (i3*i1,I2)
        ndat=(n3+1)*(n1+1+lupfil-lowfil)
        call convrot_grow(n2,ndat,y,x)
!  (i3,I1*I2) -> (iI*I2,I3)
        ndat=(n1+1+lupfil-lowfil)*(n2+1+lupfil-lowfil)
        call convrot_grow(n3,ndat,x,y)


        END SUBROUTINE convolut_magic_n


        subroutine convolut_magic_t(n1,n2,n3,x,y)
! Applies the magic filter matrix transposed ; data set shrinks
! The input array x is overwritten
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
        dimension x(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil),y(0:n1,0:n2,0:n3)

!  (I1,I2*I3) -> (I2*I3,i1)
        ndat=(n2+1+lupfil-lowfil)*(n3+1+lupfil-lowfil)
        call convrot_shrink(n1,ndat,x,y)
!  (I2,I3*i1) -> (I3*i1,i2)
        ndat=(n3+1+lupfil-lowfil)*(n1+1)
        call convrot_shrink(n2,ndat,y,x)
!  (I3,i1*i2) -> (i1*i2,i3)
        ndat=(n1+1)*(n2+1)
        call convrot_shrink(n3,ndat,x,y)

        END SUBROUTINE convolut_magic_t


        subroutine synthese_grow(n1,n2,n3,ww,x,y)
! A synthesis wavelet transformation where the size of the data is allowed to grow
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension x(0:n1,2,0:n2,2,0:n3,2)
        dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
        dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

! i1,i2,i3 -> i2,i3,I1
        nt=(2*n2+2)*(2*n3+2)
        call  syn_rot_grow(n1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(2*n3+2)*(2*n1+16)
        call  syn_rot_grow(n2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(2*n1+16)*(2*n2+16)
        call  syn_rot_grow(n3,nt,ww,y)

        END SUBROUTINE



        subroutine analyse_shrink(n1,n2,n3,ww,y,x)
! A analysis wavelet transformation where the size of the data is forced to shrink
! The input array y is overwritten
        implicit real*8 (a-h,o-z)
        dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
        dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
        dimension x(0:n1,2,0:n2,2,0:n3,2)

! I1,I2,I3 -> I2,I3,i1
        nt=(2*n2+16)*(2*n3+16)
        call  ana_rot_shrink(n1,nt,y,ww)
! I2,I3,i1 -> I3,i1,i2
        nt=(2*n3+16)*(2*n1+2)
        call  ana_rot_shrink(n2,nt,ww,y)
! I3,i1,i2 -> i1,i2,i3
        nt=(2*n1+2)*(2*n2+2)
        call  ana_rot_shrink(n3,nt,y,x)

        return
        END SUBROUTINE




!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine sumrho(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
                              nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho)
! Calculates the charge density by summing the square of all orbitals
! Input: psi
! Output: rho
        implicit real*8 (a-h,o-z)
        logical parallel
        dimension rho((2*n1+31)*(2*n2+31)*(2*n3+31)),occup(norb)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp)
        real*8, allocatable :: psig(:,:,:,:,:,:),psifscf(:),psir(:),rho_p(:)
        include 'mpif.h'


        hgridh=hgrid*.5d0 

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psig(0:n1,2,0:n2,2,0:n3,2))
        allocate( psifscf((2*n1+31)*(2*n2+31)*(2*n3+16)) )
! Wavefunction in real space
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))

 if (parallel) then
      call timing(iproc,'Rho_comput    ','ON')
        allocate(rho_p((2*n1+31)*(2*n2+31)*(2*n3+31)))

   !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenmminustwenty((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)
        !call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psig)
        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n1+2))

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=ww(((2*n1+31)*(2*n2+31)*(2*n1+16))

       do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         rho_p(i)=rho_p(i)+(occup(iorb)/hgridh**3)*psir(i)**2
        enddo

      enddo

      call timing(iproc,'Rho_comput    ','OF')
      call timing(iproc,'Rho_commun    ','ON')
        call MPI_ALLREDUCE(rho_p,rho,(2*n1+31)*(2*n2+31)*(2*n3+31),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call timing(iproc,'Rho_commun    ','OF')

        deallocate(rho_p)
 else

      call timing(iproc,'Rho_comput    ','ON')
    !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenmminustwenty((2*n1+31)*(2*n2+31)*(2*n3+31),rho)
        !call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

     do iorb=1,norb

        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psig)
        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n1+2))`

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=ww(((2*n1+31)*(2*n2+31)*(2*n1+16))

       do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         rho(i)=rho(i)+(occup(iorb)/hgridh**3)*psir(i)**2
        enddo

     enddo
      call timing(iproc,'Rho_comput    ','OF')
 endif

! Check
        tt=0.d0
        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         tt=tt+rho(i)
        enddo
        !factor of two to restore the total charge
        tt=tt*hgridh**3
        if (iproc.eq.0) write(*,*) 'Total charge from routine chargedens',tt,iproc


        deallocate(psig,psifscf,psir)


        END SUBROUTINE

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
        
!!****f* BigDFT/tenmminustwenty
!! NAME
!!   tenmminustwenty
!!
!! FUNCTION
!!   Set to 10^-20 an array x(n)
!!
!! SOURCE
!!
subroutine tenmminustwenty(n,x)
  implicit none
! Arguments
  integer :: n
  real*8 :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-20
  end do
END SUBROUTINE 
!!***



        subroutine numb_proj(ityp,ntypes,psppar,npspcode,mproj)
! Determines the number of projectors (valid for GTH and HGH pseudopotentials)
        implicit real*8 (a-h,o-z)
        dimension psppar(0:4,0:4,ntypes),npspcode(ntypes)

        mproj=0
        if (npspcode(ityp) == 2) then !GTH
           do l=1,2 
              do i=1,2 
                 if (psppar(l,i,ityp).ne.0.d0) mproj=mproj+2*l-1
              enddo
           enddo
        else if (npspcode(ityp) == 3) then !HGH
              do l=1,4 
                 do i=1,3 
                 if (psppar(l,i,ityp).ne.0.d0) mproj=mproj+2*l-1
              enddo
           enddo
        end if
        return
        END SUBROUTINE


        subroutine applyprojectorsall(iproc,ntypes,nat,iatype,psppar,npspcode,occup, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norb,norbp,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)
! Applies all the projectors onto a wavefunction
! Input: psi_c,psi_f
! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
        implicit real*8 (a-h,o-z)
        dimension psppar(0:4,0:4,ntypes),iatype(nat),npspcode(ntypes)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp)
        dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
        dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
        dimension proj(nprojel),occup(norb)

      call timing(iproc,'ApplyProj     ','ON')

  eproj_sum=0.d0
! loop over all my orbitals
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
      call applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,  & 
                    psi(1,iorb-iproc*norbp),hpsi(1,iorb-iproc*norbp),eproj)
     eproj_sum=eproj_sum+occup(iorb)*eproj
!     write(*,*) 'iorb,eproj',iorb,eproj
  enddo

      call timing(iproc,'ApplyProj     ','OF')

         END SUBROUTINE

        subroutine applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj)
! Applies all the projectors onto a single wavefunction
! Input: psi_c,psi_f
! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
        implicit real*8 (a-h,o-z)
        dimension psppar(0:4,0:4,ntypes),iatype(nat),npspcode(ntypes)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f),hpsi(nvctr_c+7*nvctr_f)
        dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
        dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
        dimension proj(nprojel)

! loop over all projectors
    iproj=0
    eproj=0.d0
    istart_c=1
    do iat=1,nat
        mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
        mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
        jseg_c=nseg_p(2*iat-2)+1
        jseg_f=nseg_p(2*iat-1)+1
        mbvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
        mbvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
     ityp=iatype(iat)
     !GTH and HGH pseudopotentials
     do l=1,4
        do i=1,3
           if (psppar(l,i,ityp).ne.0.d0) then
              do m=1,2*l-1
                 iproj=iproj+1
                 istart_f=istart_c+mbvctr_c
                 call wpdot(  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),scpr)

                 ! test (will sometimes give wrong result)
                 call wpdot(  &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),tcpr)
                 if (scpr.ne.tcpr) then
                    print *,'projectors: scpr.ne.tcpr'
                    print *,'l,i,m,h_i^l=',l,i,m,psppar(l,i,ityp)
                    stop 
                 end if
                 ! testend

                 scprp=scpr*psppar(l,i,ityp)
                 eproj=eproj+scprp*scpr

                 call waxpy(&
                      scprp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                 istart_c=istart_f+7*mbvctr_f
              enddo
              if (npspcode(ityp) == 3 .and. l/=4 .and. i/=3) then !HGH case, offdiagonal terms
                 loop_j: do j=i+1,3
                    if (psppar(l,j,ityp) .eq. 0.d0) exit loop_j
                    !calculate the coefficients for the off-diagonal terms
                    if (l==1) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(3.d0/5.d0)
                          if (j==3) offdiagcoeff=0.5d0*sqrt(5.d0/21.d0)
                       else
                          offdiagcoeff=-0.5d0*sqrt(100.d0/63.d0)
                       end if
                    else if (l==2) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(5.d0/7.d0)
                          if (j==3) offdiagcoeff=1.d0/6.d0*sqrt(35.d0/11.d0)
                       else
                          offdiagcoeff=-7.d0/3.d0*sqrt(1.d0/11.d0)
                       end if
                    else if (l==3) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(7.d0/9.d0)
                          if (j==3) offdiagcoeff=0.5d0*sqrt(63.d0/143.d0)
                       else
                          offdiagcoeff=-9.d0*sqrt(1.d0/143.d0)
                       end if
                    end if
                    istart_c_i=istart_c-(2*l-1)*(mbvctr_c+7*mbvctr_f)
                    istart_c_j=istart_c_i+(j-i)*(2*l-1)*(mbvctr_c+7*mbvctr_f)
                    do m=1,2*l-1
                       !starting addresses of the projectors
                       istart_f_j=istart_c_j+mbvctr_c
                       istart_f_i=istart_c_i+mbvctr_c
                       call wpdot(&
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_j),proj(istart_f_j),scpr_j)

                       call wpdot(&
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),psi(1),psi(nvctr_c+1),  &
                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_i),proj(istart_f_i),scpr_i)


                       scprp_j=scpr_j*offdiagcoeff*psppar(l,j,ityp)
                       scprp_i=scpr_i*offdiagcoeff*psppar(l,j,ityp)
                       !scpr_i*h_ij*scpr_j+scpr_j*h_ij*scpr_i
                       eproj=eproj+2.d0*scprp_j*scpr_i

                       !|hpsi>=|hpsi>+h_ij (<p_i|psi>|p_j>+<p_j|psi>|p_i>)
                       call waxpy(&
                            scprp_j,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                            keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_i),proj(istart_f_i),  &
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                       call waxpy(&
                            scprp_i,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                            keyv_p(jseg_c),keyv_p(jseg_f),  &
                            keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                            proj(istart_c_j),proj(istart_f_j),  &
                            nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                            keyg(1,1),keyg(1,nseg_c+1),hpsi(1),hpsi(nvctr_c+1))

                       istart_c_j=istart_f_j+7*mbvctr_f
                       istart_c_i=istart_f_i+7*mbvctr_f
                    enddo
                 end do loop_j
              end if
           end if
        enddo
     enddo
  enddo
     if (iproj.ne.nproj) stop '1:applyprojectorsone'
     if (istart_c-1.ne.nprojel) stop '2:applyprojectorsone'
         return
         END SUBROUTINE


        subroutine crtproj(iproc,nterm,n1,n2,n3, & 
                     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                     radius_f,cpmult,fpmult,hgrid,gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
                     mvctr_c,mvctr_f,proj_c,proj_f)
! returns the compressed form of a Gaussian projector 
! x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays proj_c, proj_f
        implicit real*8 (a-h,o-z)
        parameter(ntermx=3,nw=16000)
        dimension lx(nterm),ly(nterm),lz(nterm)
        dimension fac_arr(nterm)
        dimension proj_c(mvctr_c),proj_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work

        allocate(wprojx(0:n1,2,nterm),wprojy(0:n2,2,nterm),wprojz(0:n3,2,nterm),work(0:nw,2))


        rad_c=radius_f*cpmult
        rad_f=radius_f*fpmult

! make sure that the coefficients returned by CALL GAUSS_TO_DAUB are zero outside [ml:mr] 
        err_norm=0.d0 
      do 100,iterm=1,nterm
        factor=fac_arr(iterm)
        n_gau=lx(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,factor,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
        n_gau=ly(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
        n_gau=lz(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
       if (iproc.eq.0)  then
!       write(*,'(a,6(1x,i4))') 'Proj. box X: nl1_c,nl1_f,ml1,mu1,nu1_f,nu1_c',nl1_c,nl1_f,ml1,mu1,nu1_f,nu1_c
!       write(*,'(a,6(1x,i4))') 'Proj. box Y: nl2_c,nl2_f,ml2,mu2,nu2_f,nu2_c',nl2_c,nl2_f,ml2,mu2,nu2_f,nu2_c
!       write(*,'(a,6(1x,i4))') 'Proj. box Z: nl3_c,nl3_f,ml3,mu3,nu3_f,nu3_c',nl3_c,nl3_f,ml3,mu3,nu3_f,nu3_c
       if (ml1.gt.min(nl1_c,nl1_f)) write(*,*) 'Projector box larger than needed: ml1'
       if (ml2.gt.min(nl2_c,nl2_f)) write(*,*) 'Projector box larger than needed: ml2'
       if (ml3.gt.min(nl3_c,nl3_f)) write(*,*) 'Projector box larger than needed: ml3'
       if (mu1.lt.max(nu1_c,nu1_f)) write(*,*) 'Projector box larger than needed: mu1'
       if (mu2.lt.max(nu2_c,nu2_f)) write(*,*) 'Projector box larger than needed: mu2'
       if (mu3.lt.max(nu3_c,nu3_f)) write(*,*) 'Projector box larger than needed: mu3'
       endif
100   continue
!        if (iproc.eq.0) write(*,*) 'max err_norm ',err_norm

! First term: coarse projector components
!          write(*,*) 'rad_c=',rad_c
          mvctr=0
          do i3=nl3_c,nu3_c
          dz2=(i3*hgrid-rz)**2
          do i2=nl2_c,nu2_c
          dy2=(i2*hgrid-ry)**2
          do i1=nl1_c,nu1_c
          dx=i1*hgrid-rx
          if (dx**2+(dy2+dz2).lt.rad_c**2) then
            mvctr=mvctr+1
            proj_c(mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
          endif
          enddo ; enddo ; enddo
          if (mvctr.ne.mvctr_c) stop 'mvctr >< mvctr_c'

! First term: fine projector components
          mvctr=0
          do i3=nl3_f,nu3_f
          dz2=(i3*hgrid-rz)**2
          do i2=nl2_f,nu2_f
          dy2=(i2*hgrid-ry)**2
          do i1=nl1_f,nu1_f
          dx=i1*hgrid-rx
          if (dx**2+(dy2+dz2).lt.rad_f**2) then
            mvctr=mvctr+1
            proj_f(1,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
            proj_f(2,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
            proj_f(3,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
            proj_f(4,mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
            proj_f(5,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
            proj_f(6,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
            proj_f(7,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
          endif
          enddo ; enddo ; enddo
          if (mvctr.ne.mvctr_f) stop 'mvctr >< mvctr_f'
  

         do iterm=2,nterm

! Other terms: coarse projector components
         mvctr=0
         do i3=nl3_c,nu3_c
         dz2=(i3*hgrid-rz)**2
         do i2=nl2_c,nu2_c
         dy2=(i2*hgrid-ry)**2
         do i1=nl1_c,nu1_c
         dx=i1*hgrid-rx
         if (dx**2+(dy2+dz2).lt.rad_c**2) then
           mvctr=mvctr+1
           proj_c(mvctr)=proj_c(mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
         endif
         enddo ; enddo ; enddo

! Other terms: fine projector components
         mvctr=0
         do i3=nl3_f,nu3_f
         dz2=(i3*hgrid-rz)**2
         do i2=nl2_f,nu2_f
         dy2=(i2*hgrid-ry)**2
         do i1=nl1_f,nu1_f
         dx=i1*hgrid-rx
         if (dx**2+(dy2+dz2).lt.rad_f**2) then
           mvctr=mvctr+1
           proj_f(1,mvctr)=proj_f(1,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
           proj_f(2,mvctr)=proj_f(2,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
           proj_f(3,mvctr)=proj_f(3,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
           proj_f(4,mvctr)=proj_f(4,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
           proj_f(5,mvctr)=proj_f(5,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
           proj_f(6,mvctr)=proj_f(6,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
           proj_f(7,mvctr)=proj_f(7,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
         endif
         enddo ; enddo ; enddo
          

          enddo
  
          deallocate(wprojx,wprojy,wprojz,work)

    return
    END SUBROUTINE
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    subroutine wpdot(  &
         mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
! Warning: the subroutine assumes that bpsi has only one segment along each line,
! whereas apsi can have several segments. This assumption is true if bpsi is a projector 
        implicit real*8 (a-h,o-z)
        dimension keyav_c(maseg_c),keyag_c(2,maseg_c),keyav_f(maseg_f),keyag_f(2,maseg_f)
        dimension keybv_c(mbseg_c),keybg_c(2,mbseg_c),keybv_f(mbseg_f),keybg_f(2,mbseg_f)
        dimension apsi_c(mavctr_c),apsi_f(7,mavctr_f),bpsi_c(mbvctr_c),bpsi_f(7,mbvctr_f)

!        llc=0
        scpr=0.d0
! coarse part
        ibseg=1
        do iaseg=1,maseg_c
          jaj=keyav_c(iaseg)
          ja0=keyag_c(1,iaseg)
          ja1=keyag_c(2,iaseg)

100       jb1=keybg_c(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_c) goto 111
             goto 100
          endif
          jb0=keybg_c(1,ibseg)
          jbj=keybv_c(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
!           write(*,*) 'ja0,ja1,jb0,jb1',ja0,ja1,jb0,jb1,length
!          write(*,'(5(a,i5))') 'C:from ',jaj+iaoff,' to ',jaj+iaoff+length,' and from ',jbj+iboff,' to ',jbj+iboff+length
          do i=0,length
!          llc=llc+1
          scpr=scpr+apsi_c(jaj+iaoff+i)*bpsi_c(jbj+iboff+i) 
          enddo
        enddo
111     continue


!        llf=0
        scpr1=0.d0
        scpr2=0.d0
        scpr3=0.d0
        scpr4=0.d0
        scpr5=0.d0
        scpr6=0.d0
        scpr7=0.d0
! fine part
        ibseg=1
        do iaseg=1,maseg_f
          jaj=keyav_f(iaseg)
          ja0=keyag_f(1,iaseg)
          ja1=keyag_f(2,iaseg)

200       jb1=keybg_f(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_f) goto 222
             goto 200
          endif
          jb0=keybg_f(1,ibseg)
          jbj=keybv_f(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
          do i=0,length
!          llf=llf+1
          scpr1=scpr1+apsi_f(1,jaj+iaoff+i)*bpsi_f(1,jbj+iboff+i) 
          scpr2=scpr2+apsi_f(2,jaj+iaoff+i)*bpsi_f(2,jbj+iboff+i) 
          scpr3=scpr3+apsi_f(3,jaj+iaoff+i)*bpsi_f(3,jbj+iboff+i) 
          scpr4=scpr4+apsi_f(4,jaj+iaoff+i)*bpsi_f(4,jbj+iboff+i) 
          scpr5=scpr5+apsi_f(5,jaj+iaoff+i)*bpsi_f(5,jbj+iboff+i) 
          scpr6=scpr6+apsi_f(6,jaj+iaoff+i)*bpsi_f(6,jbj+iboff+i) 
          scpr7=scpr7+apsi_f(7,jaj+iaoff+i)*bpsi_f(7,jbj+iboff+i) 
          enddo
        enddo
222     continue

        scpr=scpr+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
!        write(*,*) 'llc,llf',llc,llf

    return
    END SUBROUTINE



    subroutine waxpy(  & 
        scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f, & 
        mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f)
! rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
! The update is only done in the localization region of apsi
        implicit real*8 (a-h,o-z)
        dimension keyav_c(maseg_c),keyag_c(2,maseg_c),keyav_f(maseg_f),keyag_f(2,maseg_f)
        dimension keybv_c(mbseg_c),keybg_c(2,mbseg_c),keybv_f(mbseg_f),keybg_f(2,mbseg_f)
        dimension apsi_c(mavctr_c),apsi_f(7,mavctr_f),bpsi_c(mbvctr_c),bpsi_f(7,mbvctr_f)

!        llc=0
! coarse part
        ibseg=1
        do iaseg=1,maseg_c
          jaj=keyav_c(iaseg)
          ja0=keyag_c(1,iaseg)
          ja1=keyag_c(2,iaseg)

100       jb1=keybg_c(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_c) goto 111
             goto 100
          endif
          jb0=keybg_c(1,ibseg)
          jbj=keybv_c(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
          do i=0,length
!          llc=llc+1
          apsi_c(jaj+iaoff+i)=apsi_c(jaj+iaoff+i)+scpr*bpsi_c(jbj+iboff+i) 
          enddo
        enddo
111     continue

!        llf=0
! fine part
        ibseg=1
        do iaseg=1,maseg_f
          jaj=keyav_f(iaseg)
          ja0=keyag_f(1,iaseg)
          ja1=keyag_f(2,iaseg)

200       jb1=keybg_f(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_f) goto 222
             goto 200
          endif
          jb0=keybg_f(1,ibseg)
          jbj=keybv_f(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
          do i=0,length
!          llf=llf+1
          apsi_f(1,jaj+iaoff+i)=apsi_f(1,jaj+iaoff+i)+scpr*bpsi_f(1,jbj+iboff+i) 
          apsi_f(2,jaj+iaoff+i)=apsi_f(2,jaj+iaoff+i)+scpr*bpsi_f(2,jbj+iboff+i) 
          apsi_f(3,jaj+iaoff+i)=apsi_f(3,jaj+iaoff+i)+scpr*bpsi_f(3,jbj+iboff+i) 
          apsi_f(4,jaj+iaoff+i)=apsi_f(4,jaj+iaoff+i)+scpr*bpsi_f(4,jbj+iboff+i) 
          apsi_f(5,jaj+iaoff+i)=apsi_f(5,jaj+iaoff+i)+scpr*bpsi_f(5,jbj+iboff+i) 
          apsi_f(6,jaj+iaoff+i)=apsi_f(6,jaj+iaoff+i)+scpr*bpsi_f(6,jbj+iboff+i) 
          apsi_f(7,jaj+iaoff+i)=apsi_f(7,jaj+iaoff+i)+scpr*bpsi_f(7,jbj+iboff+i) 
          enddo
        enddo
222     continue
!        write(*,*) 'waxpy,llc,llf',llc,llf

    return
    END SUBROUTINE



    subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
! calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
        implicit real*8 (a-h,o-z)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

        scpr=0.d0
    do i=1,mvctr_c
           scpr=scpr+psi_c(i)**2
        enddo
        scpr1=0.d0
        scpr2=0.d0
        scpr3=0.d0
        scpr4=0.d0
        scpr5=0.d0
        scpr6=0.d0
        scpr7=0.d0
    do i=1,mvctr_f
           scpr1=scpr1+psi_f(1,i)**2
           scpr2=scpr2+psi_f(2,i)**2
           scpr3=scpr3+psi_f(3,i)**2
           scpr4=scpr4+psi_f(4,i)**2
           scpr5=scpr5+psi_f(5,i)**2
           scpr6=scpr6+psi_f(6,i)**2
           scpr7=scpr7+psi_f(7,i)**2
        enddo
        scpr=scpr+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

    return
    END SUBROUTINE



    subroutine wscal(mvctr_c,mvctr_f,scal,psi_c,psi_f)
! multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
        implicit real*8 (a-h,o-z)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

    do i=1,mvctr_c
           psi_c(i)=psi_c(i)*scal
        enddo
    do i=1,mvctr_f
           psi_f(1,i)=psi_f(1,i)*scal
           psi_f(2,i)=psi_f(2,i)*scal
           psi_f(3,i)=psi_f(3,i)*scal
           psi_f(4,i)=psi_f(4,i)*scal
           psi_f(5,i)=psi_f(5,i)*scal
           psi_f(6,i)=psi_f(6,i)*scal
           psi_f(7,i)=psi_f(7,i)*scal
        enddo

    return
    END SUBROUTINE


    subroutine wzero(mvctr_c,mvctr_f,psi_c,psi_f)
! initializes a wavefunction to zero
        implicit real*8 (a-h,o-z)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

    do i=1,mvctr_c
           psi_c(i)=0.d0
        enddo
    do i=1,mvctr_f
           psi_f(1,i)=0.d0
           psi_f(2,i)=0.d0
           psi_f(3,i)=0.d0
           psi_f(4,i)=0.d0
           psi_f(5,i)=0.d0
           psi_f(6,i)=0.d0
           psi_f(7,i)=0.d0
        enddo

    return
    END SUBROUTINE


        subroutine orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
!Effect of orthogonality constraints on gradient 
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.true.
        dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc),occup(norb)
        allocatable :: alag(:,:,:)
        include 'mpif.h'

      call timing(iproc,'LagrM_comput  ','ON')

      allocate(alag(norb,norb,2))
!     alag(jorb,iorb,2)=+psit(k,jorb)*hpsit(k,iorb)
      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,alag(1,1,2),norb)

      call timing(iproc,'LagrM_comput  ','OF')
      call timing(iproc,'LagrM_commun  ','ON')
     call MPI_ALLREDUCE(alag(1,1,2),alag(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call timing(iproc,'LagrM_commun  ','OF')
      call timing(iproc,'LagrM_comput  ','ON')
!        if (iproc.eq.0) then
!        write(*,*) 'ALAG',iproc
!        do iorb=1,norb
!        write(*,'(10(1x,e10.3))') (alag(iorb,jorb,1),jorb=1,norb)
!        enddo
!        endif

     scprsum=0.d0
     do iorb=1,norb
         scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
     enddo

! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psit,nvctrp,alag,norb,1.d0,hpsit,nvctrp)
     deallocate(alag)
 
      call timing(iproc,'LagrM_comput  ','OF')

     END SUBROUTINE



        subroutine orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
!Effect of orthogonality constraints on gradient 
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.false.
        dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp),occup(norb)
        allocatable :: alag(:,:,:)

      call timing(iproc,'LagrM_comput  ','ON')

     allocate(alag(norb,norb,2))
 
!     alag(jorb,iorb,2)=+psi(k,jorb)*hpsi(k,iorb)
      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,alag(1,1,1),norb)

     scprsum=0.d0
     do iorb=1,norb
         scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
     enddo

! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psi,nvctrp,alag,norb,1.d0,hpsi,nvctrp)

     deallocate(alag)

      call timing(iproc,'LagrM_comput  ','OF')

     END SUBROUTINE



!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine system_size(nat,rxyz,radii,rmult,iatype,ntypes, &
                   cxmin,cxmax,cymin,cymax,czmin,czmax)
! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),radii(ntypes),iatype(nat)

        cxmax=-1.d100 ; cxmin=1.d100
        cymax=-1.d100 ; cymin=1.d100
        czmax=-1.d100 ; czmin=1.d100
        do iat=1,nat
            rad=radii(iatype(iat))*rmult
            cxmax=max(cxmax,rxyz(1,iat)+rad) ; cxmin=min(cxmin,rxyz(1,iat)-rad)
            cymax=max(cymax,rxyz(2,iat)+rad) ; cymin=min(cymin,rxyz(2,iat)-rad)
            czmax=max(czmax,rxyz(3,iat)+rad) ; czmin=min(czmin,rxyz(3,iat)-rad)
        enddo
  
      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
      czmax=czmax-eps_mach ; czmin=czmin+eps_mach

        return
        END SUBROUTINE


        subroutine pregion_size(rxyz,radii,rmult,iatype,ntypes, &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
! finds the size of the smallest subbox that contains a localization region made 
! out of atom centered spheres
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension rxyz(3),radii(ntypes)

            rad=radii(iatype)*rmult
            cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
            cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
            czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
!        write(*,*) radii(iatype(iat)),rmult
!        write(*,*) rxyz(1),rxyz(2),rxyz(3)
  
      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
      czmax=czmax-eps_mach ; czmin=czmin+eps_mach
      onem=1.d0-eps_mach
      nl1=int(onem+cxmin/hgrid)   
      nl2=int(onem+cymin/hgrid)   
      nl3=int(onem+czmin/hgrid)   
      nu1=int(cxmax/hgrid)  
      nu2=int(cymax/hgrid)  
      nu3=int(czmax/hgrid)  
!        write(*,'(a,6(i4))') 'projector region ',nl1,nu1,nl2,nu2,nl3,nu3
!        write(*,*) ' projector region size ',cxmin,cxmax
!        write(*,*) '             ',cymin,cymax
!        write(*,*) '             ',czmin,czmax
      if (nl1.lt.0)   stop 'nl1: projector region outside cell'
      if (nl2.lt.0)   stop 'nl2: projector region outside cell'
      if (nl3.lt.0)   stop 'nl3: projector region outside cell'
      if (nu1.gt.n1)   stop 'nu1: projector region outside cell'
      if (nu2.gt.n2)   stop 'nu2: projector region outside cell'
      if (nu3.gt.n3)   stop 'nu3: projector region outside cell'

        return
        END SUBROUTINE


    subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
! Calculates the length of the keys describing a wavefunction data structure
        implicit real*8 (a-h,o-z)
        logical logrid,plogrid
        dimension logrid(0:n1,0:n2,0:n3)

        mvctr=0
        nsrt=0
        nend=0
        do i3=nl3,nu3 ; do i2=nl2,nu2

        plogrid=.false.
        do i1=nl1,nu1
         if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
             nsrt=nsrt+1
           endif
         else
           if (plogrid .eqv. .true.) then
             nend=nend+1
           endif
         endif
         plogrid=logrid(i1,i2,i3)
        enddo 
           if (plogrid .eqv. .true.) then
             nend=nend+1
           endif
        enddo ; enddo
        if (nend.ne.nsrt) then 
           write(*,*) 'nend , nsrt',nend,nsrt
           stop 'nend <> nsrt'
        endif
        mseg=nend

    return
        END SUBROUTINE


    subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
! Calculates the keys describing a wavefunction data structure
        implicit real*8 (a-h,o-z)
        logical logrid,plogrid
        dimension logrid(0:n1,0:n2,0:n3),keyg(2,mseg),keyv(mseg)

        mvctr=0
        nsrt=0
        nend=0
        do i3=nl3,nu3 ; do i2=nl2,nu2

        plogrid=.false.
        do i1=nl1,nu1
         ngridp=i3*((n1+1)*(n2+1)) + i2*(n1+1) + i1+1
         if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
             nsrt=nsrt+1
             keyg(1,nsrt)=ngridp
             keyv(nsrt)=mvctr
           endif
         else
           if (plogrid .eqv. .true.) then
             nend=nend+1
             keyg(2,nend)=ngridp-1
           endif
         endif
         plogrid=logrid(i1,i2,i3)
        enddo 
           if (plogrid .eqv. .true.) then
             nend=nend+1
             keyg(2,nend)=ngridp
           endif
        enddo ; enddo
        if (nend.ne.nsrt) then 
           write(*,*) 'nend , nsrt',nend,nsrt
           stop 'nend <> nsrt'
        endif
        mseg=nend

    return
        END SUBROUTINE




       subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
                               ntypes,iatype,rxyz,radii,rmult,hgrid,logrid)
! set up an array logrid(i1,i2,i3) that specifies whether the grid point
! i1,i2,i3 is the center of a scaling function/wavelet
        implicit real*8 (a-h,o-z)
        logical logrid
        parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
        dimension rxyz(3,nat),iatype(nat),radii(ntypes)
        dimension logrid(0:n1,0:n2,0:n3)

        do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
         logrid(i1,i2,i3)=.false.
        enddo ; enddo ; enddo

      do iat=1,nat
        rad=radii(iatype(iat))*rmult+nbuf*hgrid
!        write(*,*) 'iat,nat,rad',iat,nat,rad
        ml1=int(onem+(rxyz(1,iat)-rad)/hgrid)  ; mu1=int((rxyz(1,iat)+rad)/hgrid)
        ml2=int(onem+(rxyz(2,iat)-rad)/hgrid)  ; mu2=int((rxyz(2,iat)+rad)/hgrid)
        ml3=int(onem+(rxyz(3,iat)-rad)/hgrid)  ; mu3=int((rxyz(3,iat)+rad)/hgrid)
        if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
        if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
        if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
        do i3=ml3,mu3
        dz2=(i3*hgrid-rxyz(3,iat))**2
        do i2=ml2,mu2
        dy2=(i2*hgrid-rxyz(2,iat))**2
        do i1=ml1,mu1
        dx=i1*hgrid-rxyz(1,iat)
        if (dx**2+(dy2+dz2).lt.rad**2) then 
              logrid(i1,i2,i3)=.true.
        endif
        enddo ; enddo ; enddo
      enddo

        return
        END SUBROUTINE


        subroutine bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
        implicit real*8 (a-h,o-z)
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine orthon_p(iproc,nproc,norb,norbp,nvctrp,psit)
! Gram-Schmidt orthogonalisation
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.true.
        dimension psit(nvctrp,norbp*nproc)
        real*8, allocatable :: ovrlp(:,:,:)
        include 'mpif.h'

      call timing(iproc,'GramS_comput  ','ON')

 if (norb.eq.1) stop 'more than one orbital needed for a parallel run'

        allocate(ovrlp(norb,norb,2))

! Upper triangle of overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
        call DSYRK('L','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)
        
       call timing(iproc,'GramS_comput  ','OF')
       call timing(iproc,'GramS_commun  ','ON')
        call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       call timing(iproc,'GramS_commun  ','OF')
       call timing(iproc,'GramS_comput  ','ON')

!  write(*,*) 'parallel ovrlp'
!  do i=1,norb
!  write(*,'(10(1x,e10.3))') (ovrlp(i,j,1),j=1,norb)
!  enddo

       
       ! Cholesky factorization
        call dpotrf( 'L', norb, ovrlp, norb, info )
        if (info.ne.0) write(6,*) 'info Cholesky factorization', info

! calculate L^{-1}
        call DTRTRI( 'L', 'N', norb, ovrlp, norb, info )
        if (info.ne.0) write(6,*) 'info L^-1', info

! new vectors   
      call DTRMM ('R', 'L', 'T', 'N', nvctrp, norb, 1.d0, ovrlp, norb, psit, nvctrp)

        deallocate(ovrlp)

       call timing(iproc,'GramS_comput  ','OF')

    END SUBROUTINE orthon_p


        
       subroutine orthon(norb,norbp,nvctrp,psi)
! Gram-Schmidt orthogonalisation
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.false.
        dimension psi(nvctrp,norbp)
        real*8, allocatable :: ovrlp(:,:)

       call timing(iproc,'GramS_comput  ','ON')

 if (norb.eq.1) then
        tt=dnrm2(nvctrp,psi,1)
        tt=1.d0/tt
        call dscal(nvctrp,tt,psi,1)

 else

        allocate(ovrlp(norb,norb))

! Overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
        call DSYRK('L','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp,norb)

!  write(*,*) 'ovrlp'
!  do i=1,norb
!  write(*,'(10(1x,e10.3))') (ovrlp(i,j),j=1,norb)
!  enddo

! Cholesky factorization
        call dpotrf( 'L', norb, ovrlp, norb, info )
        if (info.ne.0) write(6,*) 'info Cholesky factorization', info

! calculate L^{-1}
        call DTRTRI( 'L', 'N', norb, ovrlp, norb, info )
        if (info.ne.0) write(6,*) 'info L^-1', info

! new vectors   
      call DTRMM ('R', 'L', 'T', 'N', nvctrp, norb, 1.d0, ovrlp, norb, psi, nvctrp)

      deallocate(ovrlp)

 endif

       call timing(iproc,'GramS_comput  ','OF')

end  subroutine orthon

subroutine createWavefunctionsDescriptors(parallel, iproc, nproc, idsx, n1, n2, n3, output_grid, &
     & hgrid, nat, ntypes, iatype, atomnames, rxyz, radii_cf, crmult, frmult, &
     & ibyz_c,ibxz_c,ibxy_c, ibyz_f, ibxz_f, ibxy_f, nseg_c, nseg_f, nvctr_c, nvctr_f, nvctrp, &
     & keyg, keyv,norb,norbp,psi,hpsi,psit,psidst,hpsidst,ads)
!calculates the descriptor arrays keyg and keyv as well as nseg_c, nseg_f, nvctr_c, nvctr_f, nvctrp
!calculates also the arrays ibyz_c,ibxz_c,ibxy_c, ibyz_f, ibxz_f, ibxy_f needed for convolut_standard
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,idsx,n1,n2,n3,nat,ntypes,norb
  integer, intent(in) :: nseg_c,nseg_f,nvctr_c,nvctr_f
  integer, intent(out) :: norbp,nvctrp
  logical, intent(in) :: parallel, output_grid
  integer, intent(in) :: iatype(nat)
  real*8, intent(in) :: hgrid,crmult,frmult
  integer, intent(in) :: ibyz_c(2,0:n2,0:n3), ibxz_c(2,0:n1,0:n3), ibxy_c(2,0:n1,0:n2)
  integer, intent(in) :: ibyz_f(2,0:n2,0:n3), ibxz_f(2,0:n1,0:n3), ibxy_f(2,0:n1,0:n2)
  real*8 :: rxyz(3, nat), radii_cf(ntypes, 2)
  character(len=20), intent(in) :: atomnames(100)
  integer, pointer :: keyg(:,:), keyv(:)
  ! wavefunction 
  real*8, pointer :: psi(:,:)
  real*8, pointer :: psit(:,:)
  ! wavefunction gradients
  real*8, pointer :: hpsi(:,:),hpsit(:,:)
  ! arrays for DIIS convergence accelerator
  real*8, pointer :: ads(:,:,:),psidst(:,:,:),hpsidst(:,:,:)
  !Local variables
  real*8, parameter :: eps_mach=1.d-12,onem=1.d0-eps_mach
  integer :: iat,i1,i2,i3,norbme
  real*8 :: tt
  logical, allocatable :: logrid_c(:,:,:), logrid_f(:,:,:)

  call timing(iproc,'CrtDescriptors','ON')

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3))
  allocate(logrid_f(0:n1,0:n2,0:n3))

  ! coarse grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,1),crmult,hgrid,logrid_c)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3 ; do i2=0,n2 ; do i1=0,n1
        if (logrid_c(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
     enddo; enddo ; enddo 
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,nvctr_c)
  if (iproc.eq.0) write(*,*) 'orbitals have coarse segment, elements',nseg_c,nvctr_c
  call bounds(n1,n2,n3,logrid_c,ibyz_c,ibxz_c,ibxy_c)

  ! fine grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,2),frmult,hgrid,logrid_f)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3 ; do i2=0,n2 ; do i1=0,n1
        if (logrid_f(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
     enddo; enddo ; enddo 
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,nvctr_f)
  if (iproc.eq.0) write(*,*) 'orbitals have fine   segment, elements',nseg_f,7*nvctr_f
  call bounds(n1,n2,n3,logrid_f,ibyz_f,ibxz_f,ibxy_f)

  if (iproc.eq.0 .and. output_grid) close(22)

  ! allocations for arrays holding the wavefunctions and their data descriptors
  allocate(keyg(2,nseg_c+nseg_f))
  allocate(keyv(nseg_c+nseg_f))

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,keyg(:,1:nseg_c),keyv(1:nseg_c))

  ! fine grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,keyg(:,nseg_c+1:nseg_c+nseg_f), &
    & keyv(nseg_c+1:nseg_c+nseg_f))

  deallocate(logrid_c,logrid_f)

! allocate wavefunction arrays
  tt=dble(norb)/dble(nproc)
  norbp=int((1.d0-eps_mach*tt) + tt)
  write(*,*) 'norbp=',norbp
  allocate(psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp))
  norbme=max(min((iproc+1)*norbp,norb)-iproc*norbp,0)
  write(*,*) 'iproc ',iproc,' treats ',norbme,' orbitals '

  tt=dble(nvctr_c+7*nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)
  if (parallel) then
     if (iproc.eq.0) write(79,'(a40,i10)') 'words for psit',nvctrp*norbp*nproc
     allocate(psit(nvctrp,norbp*nproc))
     if (iproc.eq.0) write(79,*) 'allocation done'
  endif

! allocate arrays necessary for DIIS convergence acceleration
  if (idsx.gt.0) then
     if (iproc.eq.0) write(79,'(a40,i10)') 'words for psidst and hpsidst',2*nvctrp*norbp*nproc*idsx
     allocate( psidst(nvctrp,norbp*nproc,idsx))
     allocate(hpsidst(nvctrp,norbp*nproc,idsx))
     if (iproc.eq.0) write(79,*) 'allocation done'
     allocate(ads(idsx+1,idsx+1,3))
     call razero(3*(idsx+1)**2,ads)
  endif

  call timing(iproc,'CrtDescriptors','OF')

END SUBROUTINE createWavefunctionsDescriptors

subroutine createKernel(parallel, nfft1, nfft2, nfft3, n1, n2, n3, hgridh, &
     & ndegree_ip, iproc, nproc, pkernel)
  implicit real*8 (a-h,o-z)
  real*8, pointer :: pkernel(:)
  logical :: parallel
  include "mpif.h"

  call timing(iproc,'PSolvKernel   ','ON')

  if (parallel) then
     call calculate_pardimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,nf1,nf2,nf3,md1,md2,md3,nfft1,nfft2,nfft3,nproc)
     !call Dimensions_FFT(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3)
     if (iproc.eq.0) then
        write(*,*) 'dimension of FFT grid',nf1,nf2,nf3
        write(*,*) 'dimension of kernel',nfft1,nfft2,nfft3/nproc
        if (iproc.eq.0) write(79,'(a40,i10)') 'words for kernel ',nfft1*nfft2*nfft3/nproc
     endif
     allocate(pkernel(nfft1*nfft2*nfft3/nproc))
     if (iproc.eq.0) write(79,*) 'allocation done'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call ParBuild_Kernel(2*n1+31,2*n2+31,2*n3+31,nf1,nf2,nf3,nfft1,nfft2,nfft3, &
          hgridh,ndegree_ip,iproc,nproc,pkernel)
     if (iproc.eq.0) print *,"Poisson Solver Kernel built!"
     !          call PARtest_kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,pot_ion,rhopot,iproc,nproc) 

  else
     call Dimensions_FFT(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3)
     write(*,*) 'dimension of FFT grid',nfft1,nfft2,nfft3
     write(*,*) 'dimension of kernel',nfft1/2+1,nfft2/2+1,nfft3/2+1
     if (iproc.eq.0) write(79,'(a40,i10)') 'words for kernel ',(nfft1/2+1)*(nfft2/2+1)*(nfft3/2+1)
     allocate(pkernel((nfft1/2+1)*(nfft2/2+1)*(nfft3/2+1)))
     if (iproc.eq.0) write(79,*) 'allocation done'
     call Build_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3, &
          hgridh,ndegree_ip,pkernel)

     !          call test_kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,pot_ion,rhopot)
  end if

  call timing(iproc,'PSolvKernel   ','OF')

END SUBROUTINE createKernel

  subroutine createProjectorsArrays(iproc, n1, n2, n3, rxyz, nat, ntypes, iatype, atomnames, &
       & psppar, npspcode, radii_cf, cpmult, fpmult, hgrid, nvctr_p, nseg_p, &
       & keyg_p, keyv_p, nproj, nprojel, istart, nboxp_c, nboxp_f, proj)
    implicit real*8 (a-h,o-z)
    character*20 :: atomnames(100)
    dimension rxyz(3,nat),iatype(nat),radii_cf(ntypes,2),psppar(0:4,0:4,ntypes),npspcode(ntypes)
    integer :: nvctr_p(0:2*nat), nseg_p(0:2*nat)
    integer :: nboxp_c(2,3,nat), nboxp_f(2,3,nat)
    real*8, pointer :: proj(:)
    integer, pointer :: keyg_p(:,:), keyv_p(:)
    real(kind=8), dimension(:), allocatable :: fac_arr
    integer, dimension(:), allocatable :: lx,ly,lz

    logical, allocatable :: logrid(:,:,:)

  call timing(iproc,'CrtProjectors ','ON')

    ! determine localization region for all projectors, but do not yet fill the descriptor arrays
    allocate(logrid(0:n1,0:n2,0:n3))

    nseg_p(0)=0 
    nvctr_p(0)=0 

    istart=1
    nproj=0
    do iat=1,nat

       if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++ iat=',iat

       call numb_proj(iatype(iat),ntypes,psppar,npspcode,mproj)
       if (mproj.ne.0) then 

          if (iproc.eq.0) write(*,*) 'projector descriptors for atom with mproj ',iat,mproj
          nproj=nproj+mproj

          ! coarse grid quantities
          call  pregion_size(rxyz(1,iat),radii_cf(1,2),cpmult,iatype(iat),ntypes, &
               hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
!          if (iproc.eq.0) write(*,'(a,6(i4))') 'coarse grid',nl1,nu1,nl2,nu2,nl3,nu3
          nboxp_c(1,1,iat)=nl1 ; nboxp_c(2,1,iat)=nu1
          nboxp_c(1,2,iat)=nl2 ; nboxp_c(2,2,iat)=nu2
          nboxp_c(1,3,iat)=nl3 ; nboxp_c(2,3,iat)=nu3
          call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)
          call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
          if (iproc.eq.0) write(*,*) 'mseg,mvctr,coarse projectors ',mseg,mvctr
          nseg_p(2*iat-1)=nseg_p(2*iat-2) + mseg
          nvctr_p(2*iat-1)=nvctr_p(2*iat-2) + mvctr
          istart=istart+mvctr*mproj

          ! fine grid quantities
          call  pregion_size(rxyz(1,iat),radii_cf(1,2),fpmult,iatype(iat),ntypes, &
               hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
!          if (iproc.eq.0) write(*,'(a,6(i4))') 'fine   grid',nl1,nu1,nl2,nu2,nl3,nu3
          nboxp_f(1,1,iat)=nl1 ; nboxp_f(2,1,iat)=nu1
          nboxp_f(1,2,iat)=nl2 ; nboxp_f(2,2,iat)=nu2
          nboxp_f(1,3,iat)=nl3 ; nboxp_f(2,3,iat)=nu3
          call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
          call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
          if (iproc.eq.0) write(*,*) 'mseg,mvctr, fine  projectors ',mseg,mvctr
          nseg_p(2*iat)=nseg_p(2*iat-1) + mseg
          nvctr_p(2*iat)=nvctr_p(2*iat-1) + mvctr
          istart=istart+7*mvctr*mproj

       else  !(atom has no nonlocal PSP, e.g. H)
          nseg_p(2*iat-1)=nseg_p(2*iat-2) 
          nvctr_p(2*iat-1)=nvctr_p(2*iat-2) 
          nseg_p(2*iat)=nseg_p(2*iat-1) 
          nvctr_p(2*iat)=nvctr_p(2*iat-1) 
       endif
    enddo

    if (iproc.eq.0) write(*,*) 'total number of projectors',nproj
    ! allocations for arrays holding the projectors and their data descriptors
    allocate(keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat)))
    nprojel=istart-1
    if (iproc.eq.0) write(79,'(a40,i10)') 'words for proj ',nprojel
    allocate(proj(nprojel))
    if (iproc.eq.0) write(79,*) 'allocation done'


    if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    ! After having determined the size of the projector descriptor arrays fill them
    istart_c=1
    do iat=1,nat
       call numb_proj(iatype(iat),ntypes,psppar,npspcode,mproj)
       if (mproj.ne.0) then 

          ! coarse grid quantities
          nl1=nboxp_c(1,1,iat) ; nu1=nboxp_c(2,1,iat)
          nl2=nboxp_c(1,2,iat) ; nu2=nboxp_c(2,2,iat)
          nl3=nboxp_c(1,3,iat) ; nu3=nboxp_c(2,3,iat)
          call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)

          iseg=nseg_p(2*iat-2)+1
          mseg=nseg_p(2*iat-1)-nseg_p(2*iat-2)
          call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
               logrid,mseg,keyg_p(:,iseg:iseg+mseg-1),keyv_p(iseg:iseg+mseg-1))

          ! fine grid quantities
          nl1=nboxp_f(1,1,iat) ; nu1=nboxp_f(2,1,iat)
          nl2=nboxp_f(1,2,iat) ; nu2=nboxp_f(2,2,iat)
          nl3=nboxp_f(1,3,iat) ; nu3=nboxp_f(2,3,iat)
          call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
          iseg=nseg_p(2*iat-1)+1
          mseg=nseg_p(2*iat)-nseg_p(2*iat-1)
          call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
               logrid,mseg,keyg_p(:,iseg:iseg+mseg-1),keyv_p(iseg:iseg+mseg-1))

       endif
    enddo

    if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    !allocate these vectors up to the maximum size we can get
    nterm_max=10 !if GTH nterm_max=3
    allocate(fac_arr(nterm_max))
    allocate(lx(nterm_max),ly(nterm_max),lz(nterm_max))

    iproj=0
    fpi=(4.d0*atan(1.d0))**(-.75d0)
    do iat=1,nat
       rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
       ityp=iatype(iat)

       !decide the loop bounds
       do l=1,4 !generic case, also for HGH (for GTH it will stop at l=2)
          do i=1,3 !generic case, also for HGH (for GTH it will stop at i=2)
             if (psppar(l,i,ityp).ne.0.d0) then
                gau_a=psppar(l,0,ityp)
                factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
                do m=1,2*l-1
                   mvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
                   mvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
                   istart_f=istart_c+mvctr_c
                   nl1_c=nboxp_c(1,1,iat) ; nu1_c=nboxp_c(2,1,iat)
                   nl2_c=nboxp_c(1,2,iat) ; nu2_c=nboxp_c(2,2,iat)
                   nl3_c=nboxp_c(1,3,iat) ; nu3_c=nboxp_c(2,3,iat)
                   nl1_f=nboxp_f(1,1,iat) ; nu1_f=nboxp_f(2,1,iat)
                   nl2_f=nboxp_f(1,2,iat) ; nu2_f=nboxp_f(2,2,iat)
                   nl3_f=nboxp_f(1,3,iat) ; nu3_f=nboxp_f(2,3,iat)

                   call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)

                   fac_arr(1:nterm)=factor*fac_arr(1:nterm)

                   call crtproj(iproc,nterm,n1,n2,n3,nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c, &
                     & nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,radii_cf(iatype(iat),2), & 
                     & cpmult,fpmult,hgrid,gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
                     & mvctr_c,mvctr_f,proj(istart_c:istart_c+mvctr_c-1), &
                     & proj(istart_f:istart_f+7*mvctr_f-1))

                   iproj=iproj+1
                   ! testing
                   call wnrm(mvctr_c,mvctr_f,proj(istart_c:istart_c+mvctr_c-1), &
                     & proj(istart_f:istart_f + 7 * mvctr_f - 1),scpr)
                   if (abs(1.d0-scpr).gt.1.d-1) then
                      !print *,'norm projector','iproc=',iproc,i,j,scpr,nadd_arr,nterm
                      stop 'norm projector'
                   end if
                   ! testing end
                   istart_c=istart_f+7*mvctr_f
                   if (istart_c.gt.istart) stop 'istart_c > istart'

                   do iterm=1,nterm
                      if (iproc.eq.0) write(*,'(a,i3,1x,a,e10.3,3(i2))') 'projector: iat,atomname,gau_a,lx,ly,lz,err', & 
                           iat,atomnames(iatype(iat)),gau_a,lx(iterm),ly(iterm),lz(iterm)
                   enddo


                enddo
             endif
          enddo
       enddo
    enddo
    if (iproj.ne.nproj) stop 'incorrect number of projectors created'
    ! projector part finished
    deallocate(logrid)
    deallocate(fac_arr,lx,ly,lz)
  call timing(iproc,'CrtProjectors ','OF')

END SUBROUTINE 


subroutine calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)
  
  implicit none
  integer, intent(in) :: l,i,m,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(kind=8), dimension(nterm_max), intent(out) :: fac_arr

  if (l.eq.1 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.7071067811865475244008444d0
  else if (l.eq.1 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3651483716701107423046465d0
     fac_arr(2)=0.3651483716701107423046465d0
     fac_arr(3)=0.3651483716701107423046465d0
  else if (l.eq.1 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.09200874124564722903948358d0
     fac_arr(2)=0.1840174824912944580789672d0
     fac_arr(3)=0.09200874124564722903948358d0
     fac_arr(4)=0.1840174824912944580789672d0
     fac_arr(5)=0.1840174824912944580789672d0
     fac_arr(6)=0.09200874124564722903948358d0
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000d0
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000d0
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.000000000000000000000000d0
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473d0
     fac_arr(2)=0.3380617018914066310038473d0
     fac_arr(3)=0.3380617018914066310038473d0
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473d0
     fac_arr(2)=0.3380617018914066310038473d0
     fac_arr(3)=0.3380617018914066310038473d0
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3380617018914066310038473d0
     fac_arr(2)=0.3380617018914066310038473d0
     fac_arr(3)=0.3380617018914066310038473d0
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187d0
     fac_arr(2)=0.1359059177167001452365437d0
     fac_arr(3)=0.06795295885835007261827187d0
     fac_arr(4)=0.1359059177167001452365437d0
     fac_arr(5)=0.1359059177167001452365437d0
     fac_arr(6)=0.06795295885835007261827187d0
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187d0
     fac_arr(2)=0.1359059177167001452365437d0
     fac_arr(3)=0.06795295885835007261827187d0
     fac_arr(4)=0.1359059177167001452365437d0
     fac_arr(5)=0.1359059177167001452365437d0
     fac_arr(6)=0.06795295885835007261827187d0
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.06795295885835007261827187d0
     fac_arr(2)=0.1359059177167001452365437d0
     fac_arr(3)=0.06795295885835007261827187d0
     fac_arr(4)=0.1359059177167001452365437d0
     fac_arr(5)=0.1359059177167001452365437d0
     fac_arr(6)=0.06795295885835007261827187d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.414213562373095048801689d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.7071067811865475244008444d0
     fac_arr(2)=-0.7071067811865475244008444d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.5) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.4082482904638630163662140d0
     fac_arr(2)=-0.4082482904638630163662140d0
     fac_arr(3)=0.8164965809277260327324280d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=1
     lx(2)=0 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046d0
     fac_arr(2)=0.3563483225498991795794046d0
     fac_arr(3)=0.3563483225498991795794046d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=1
     lx(2)=1 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046d0
     fac_arr(2)=0.3563483225498991795794046d0
     fac_arr(3)=0.3563483225498991795794046d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=0
     lx(2)=1 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3563483225498991795794046d0
     fac_arr(2)=0.3563483225498991795794046d0
     fac_arr(3)=0.3563483225498991795794046d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.4) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=4 ; lz(2)=0
     lx(3)=2 ; ly(3)=0 ; lz(3)=2
     lx(4)=0 ; ly(4)=2 ; lz(4)=2
     fac_arr(1)=0.1781741612749495897897023d0
     fac_arr(2)=-0.1781741612749495897897023d0
     fac_arr(3)=0.1781741612749495897897023d0
     fac_arr(4)=-0.1781741612749495897897023d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.5) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=-0.1028688999747279401740630d0
     fac_arr(2)=-0.2057377999494558803481260d0
     fac_arr(3)=-0.1028688999747279401740630d0
     fac_arr(4)=0.1028688999747279401740630d0
     fac_arr(5)=0.1028688999747279401740630d0
     fac_arr(6)=0.2057377999494558803481260d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=1
     lx(2)=2 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=5 ; lz(3)=1
     lx(4)=2 ; ly(4)=1 ; lz(4)=3
     lx(5)=0 ; ly(5)=3 ; lz(5)=3
     lx(6)=0 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993d0
     fac_arr(2)=0.1191973750047131197905399d0
     fac_arr(3)=0.05959868750235655989526993d0
     fac_arr(4)=0.1191973750047131197905399d0
     fac_arr(5)=0.1191973750047131197905399d0
     fac_arr(6)=0.05959868750235655989526993d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=1
     lx(2)=3 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=4 ; lz(3)=1
     lx(4)=3 ; ly(4)=0 ; lz(4)=3
     lx(5)=1 ; ly(5)=2 ; lz(5)=3
     lx(6)=1 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993d0
     fac_arr(2)=0.1191973750047131197905399d0
     fac_arr(3)=0.05959868750235655989526993d0
     fac_arr(4)=0.1191973750047131197905399d0
     fac_arr(5)=0.1191973750047131197905399d0
     fac_arr(6)=0.05959868750235655989526993d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=0
     lx(2)=3 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=5 ; lz(3)=0
     lx(4)=3 ; ly(4)=1 ; lz(4)=2
     lx(5)=1 ; ly(5)=3 ; lz(5)=2
     lx(6)=1 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.05959868750235655989526993d0
     fac_arr(2)=0.1191973750047131197905399d0
     fac_arr(3)=0.05959868750235655989526993d0
     fac_arr(4)=0.1191973750047131197905399d0
     fac_arr(5)=0.1191973750047131197905399d0
     fac_arr(6)=0.05959868750235655989526993d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.4) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=4 ; ly(5)=0 ; lz(5)=2
     lx(6)=0 ; ly(6)=4 ; lz(6)=2
     lx(7)=2 ; ly(7)=0 ; lz(7)=4
     lx(8)=0 ; ly(8)=2 ; lz(8)=4
     fac_arr(1)=0.02979934375117827994763496d0
     fac_arr(2)=0.02979934375117827994763496d0
     fac_arr(3)=-0.02979934375117827994763496d0
     fac_arr(4)=-0.02979934375117827994763496d0
     fac_arr(5)=0.05959868750235655989526993d0
     fac_arr(6)=-0.05959868750235655989526993d0
     fac_arr(7)=0.02979934375117827994763496d0
     fac_arr(8)=-0.02979934375117827994763496d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.5) then
     nterm=7
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=2 ; ly(5)=0 ; lz(5)=4
     lx(6)=0 ; ly(6)=2 ; lz(6)=4
     lx(7)=0 ; ly(7)=0 ; lz(7)=6
     fac_arr(1)=-0.01720465913641697233541246d0
     fac_arr(2)=-0.05161397740925091700623738d0
     fac_arr(3)=-0.05161397740925091700623738d0
     fac_arr(4)=-0.01720465913641697233541246d0
     fac_arr(5)=0.05161397740925091700623738d0
     fac_arr(6)=0.05161397740925091700623738d0
     fac_arr(7)=0.03440931827283394467082492d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894d0
     fac_arr(2)=0.3162277660168379331998894d0
     fac_arr(3)=-1.264911064067351732799557d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894d0
     fac_arr(2)=0.3162277660168379331998894d0
     fac_arr(3)=-1.264911064067351732799557d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.7745966692414833770358531d0
     fac_arr(2)=0.7745966692414833770358531d0
     fac_arr(3)=-0.5163977794943222513572354d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.4082482904638630163662140d0
     fac_arr(2)=-1.224744871391589049098642d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-1.224744871391589049098642d0
     fac_arr(2)=0.4082482904638630163662140d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.000000000000000000000000d0
     fac_arr(2)=-1.000000000000000000000000d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.000000000000000000000000d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506d0
     fac_arr(2)=0.1271283452327456420595701d0
     fac_arr(3)=0.06356417261637282102978506d0
     fac_arr(4)=-0.1906925178491184630893552d0
     fac_arr(5)=-0.1906925178491184630893552d0
     fac_arr(6)=-0.2542566904654912841191402d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506d0
     fac_arr(2)=0.1271283452327456420595701d0
     fac_arr(3)=0.06356417261637282102978506d0
     fac_arr(4)=-0.1906925178491184630893552d0
     fac_arr(5)=-0.1906925178491184630893552d0
     fac_arr(6)=-0.2542566904654912841191402d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.1556997888323045941832351d0
     fac_arr(2)=0.3113995776646091883664703d0
     fac_arr(3)=0.1556997888323045941832351d0
     fac_arr(4)=0.05189992961076819806107838d0
     fac_arr(5)=0.05189992961076819806107838d0
     fac_arr(6)=-0.1037998592215363961221568d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.4) then
     nterm=5
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     fac_arr(1)=0.08206099398622182182282711d0
     fac_arr(2)=-0.1641219879724436436456542d0
     fac_arr(3)=-0.2461829819586654654684813d0
     fac_arr(4)=0.08206099398622182182282711d0
     fac_arr(5)=-0.2461829819586654654684813d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.5) then
     nterm=5
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     fac_arr(1)=-0.2461829819586654654684813d0
     fac_arr(2)=-0.1641219879724436436456542d0
     fac_arr(3)=0.08206099398622182182282711d0
     fac_arr(4)=-0.2461829819586654654684813d0
     fac_arr(5)=0.08206099398622182182282711d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.6) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=4 ; lz(2)=1
     lx(3)=2 ; ly(3)=0 ; lz(3)=3
     lx(4)=0 ; ly(4)=2 ; lz(4)=3
     fac_arr(1)=0.2010075630518424150978747d0
     fac_arr(2)=-0.2010075630518424150978747d0
     fac_arr(3)=0.2010075630518424150978747d0
     fac_arr(4)=-0.2010075630518424150978747d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.7) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=1
     lx(2)=1 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.4020151261036848301957494d0
     fac_arr(2)=0.4020151261036848301957494d0
     fac_arr(3)=0.4020151261036848301957494d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.1) then
     nterm=10
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     lx(10)=1 ; ly(10)=0 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687d0
     fac_arr(2)=0.02731154967995675489524106d0
     fac_arr(3)=0.02731154967995675489524106d0
     fac_arr(4)=0.009103849893318918298413687d0
     fac_arr(5)=-0.01820769978663783659682737d0
     fac_arr(6)=-0.03641539957327567319365475d0
     fac_arr(7)=-0.01820769978663783659682737d0
     fac_arr(8)=-0.06372694925323242808889581d0
     fac_arr(9)=-0.06372694925323242808889581d0
     fac_arr(10)=-0.03641539957327567319365475d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.2) then
     nterm=10
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     lx(10)=0 ; ly(10)=1 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687d0
     fac_arr(2)=0.02731154967995675489524106d0
     fac_arr(3)=0.02731154967995675489524106d0
     fac_arr(4)=0.009103849893318918298413687d0
     fac_arr(5)=-0.01820769978663783659682737d0
     fac_arr(6)=-0.03641539957327567319365475d0
     fac_arr(7)=-0.01820769978663783659682737d0
     fac_arr(8)=-0.06372694925323242808889581d0
     fac_arr(9)=-0.06372694925323242808889581d0
     fac_arr(10)=-0.03641539957327567319365475d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.3) then
     nterm=10
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=2 ; ly(6)=2 ; lz(6)=3
     lx(7)=0 ; ly(7)=4 ; lz(7)=3
     lx(8)=2 ; ly(8)=0 ; lz(8)=5
     lx(9)=0 ; ly(9)=2 ; lz(9)=5
     lx(10)=0 ; ly(10)=0 ; lz(10)=7
     fac_arr(1)=0.02229978693352242055222348d0
     fac_arr(2)=0.06689936080056726165667044d0
     fac_arr(3)=0.06689936080056726165667044d0
     fac_arr(4)=0.02229978693352242055222348d0
     fac_arr(5)=0.02973304924469656073629797d0
     fac_arr(6)=0.05946609848939312147259594d0
     fac_arr(7)=0.02973304924469656073629797d0
     fac_arr(8)=-0.007433262311174140184074493d0
     fac_arr(9)=-0.007433262311174140184074493d0
     fac_arr(10)=-0.01486652462234828036814899d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.4) then
     nterm=9
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     fac_arr(1)=0.01175301967439877980816756d0
     fac_arr(2)=-0.01175301967439877980816756d0
     fac_arr(3)=-0.05876509837199389904083778d0
     fac_arr(4)=-0.03525905902319633942450267d0
     fac_arr(5)=0.02350603934879755961633511d0
     fac_arr(6)=-0.04701207869759511923267022d0
     fac_arr(7)=-0.07051811804639267884900533d0
     fac_arr(8)=0.01175301967439877980816756d0
     fac_arr(9)=-0.03525905902319633942450267d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.5) then
     nterm=9
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     fac_arr(1)=-0.03525905902319633942450267d0
     fac_arr(2)=-0.05876509837199389904083778d0
     fac_arr(3)=-0.01175301967439877980816756d0
     fac_arr(4)=0.01175301967439877980816756d0
     fac_arr(5)=-0.07051811804639267884900533d0
     fac_arr(6)=-0.04701207869759511923267022d0
     fac_arr(7)=0.02350603934879755961633511d0
     fac_arr(8)=-0.03525905902319633942450267d0
     fac_arr(9)=0.01175301967439877980816756d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.6) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=0 ; ly(6)=4 ; lz(6)=3
     lx(7)=2 ; ly(7)=0 ; lz(7)=5
     lx(8)=0 ; ly(8)=2 ; lz(8)=5
     fac_arr(1)=0.02878890113916869875409405d0
     fac_arr(2)=0.02878890113916869875409405d0
     fac_arr(3)=-0.02878890113916869875409405d0
     fac_arr(4)=-0.02878890113916869875409405d0
     fac_arr(5)=0.05757780227833739750818811d0
     fac_arr(6)=-0.05757780227833739750818811d0
     fac_arr(7)=0.02878890113916869875409405d0
     fac_arr(8)=-0.02878890113916869875409405d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.7) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=1
     lx(2)=3 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=5 ; lz(3)=1
     lx(4)=3 ; ly(4)=1 ; lz(4)=3
     lx(5)=1 ; ly(5)=3 ; lz(5)=3
     lx(6)=1 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05757780227833739750818811d0
     fac_arr(2)=0.1151556045566747950163762d0
     fac_arr(3)=0.05757780227833739750818811d0
     fac_arr(4)=0.1151556045566747950163762d0
     fac_arr(5)=0.1151556045566747950163762d0
     fac_arr(6)=0.05757780227833739750818811d0

  else
     stop 'PSP format error'
  endif
  
END SUBROUTINE calc_coeff_proj


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine loewe_p(iproc,nproc,norb,norbp,nvctrp,psit)
! loewdin orthogonalisation
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.true.
        dimension psit(nvctrp,norbp*nproc)
        real*8, allocatable :: ovrlp(:,:,:),evall(:),psitt(:,:)
        include 'mpif.h'

 if (norb.eq.1) stop 'more than one orbital needed for a parallel run'

        allocate(ovrlp(norb,norb,3),evall(norb))

! Upper triangle of overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
        call DSYRK('U','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

! Full overlap matrix using  BLAS
!     ovrlap(jorb,iorb,2)=+psit(k,jorb)*psit(k,iorb)
!      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,nvctrp,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

        call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!       write(*,*) 'OVERLAP',iproc
!       do i=1,norb
!       write(*,'(10(x,e17.10))') (ovrlp(i,j,1),j=1,norb)
!       enddo

! LAPACK
        call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
        if (info.ne.0) write(6,*) 'info loewe', info
!        if (iproc.eq.0) then 
!          write(6,*) 'overlap eigenvalues'
!77        format(8(1x,e10.3))
!          if (norb.le.16) then
!          write(6,77) evall
!          else
!          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
!          endif
!        endif

! calculate S^{-1/2} ovrlp(*,*,3)
        do 2935,lorb=1,norb
        do 2935,jorb=1,norb
2935    ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
!        do 3985,j=1,norb
!        do 3985,i=1,norb
!        ovrlp(i,j,3)=0.d0
!        do 3985,l=1,norb
!3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
! BLAS:
        call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

        allocate(psitt(nvctrp,norbp*nproc))
! new eigenvectors
!   psitt(i,iorb)=psit(i,jorb)*ovrlp(jorb,iorb,3)
      call DGEMM('N','N',nvctrp,norb,norb,1.d0,psit,nvctrp,ovrlp(1,1,3),norb,0.d0,psitt,nvctrp)
      call DCOPY(nvctrp*norbp*nproc,psitt,1,psit,1)
      deallocate(psitt)

        deallocate(ovrlp,evall)

        END SUBROUTINE



        subroutine loewe(norb,norbp,nvctrp,psi)
! loewdin orthogonalisation
        implicit real*8 (a-h,o-z)
        dimension psi(nvctrp,norbp)
        real*8, allocatable :: ovrlp(:,:,:),evall(:),tpsi(:,:)

 if (norb.eq.1) then
        tt=0.d0
        do i=1,nvctrp
        tt=tt+psi(i,1)**2
        enddo
        tt=1.d0/sqrt(tt)
        do i=1,nvctrp
        psi(i,1)=psi(i,1)*tt
        enddo

 else

        allocate(ovrlp(norb,norb,3),evall(norb))

! Overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
        call DSYRK('U','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp(1,1,1),norb)

!       write(*,*) 'OVERLAP'
!       do i=1,norb
!       write(*,'(10(x,e17.10))') (ovrlp(i,j,1),j=1,norb)
!       enddo


! LAPACK
        call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
        if (info.ne.0) write(6,*) 'info loewe', info
!          write(6,*) 'overlap eigenvalues'
!77        format(8(1x,e10.3))
!          if (norb.le.16) then
!          write(6,77) evall
!          else
!          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
!          endif

! calculate S^{-1/2} ovrlp(*,*,3)
        do 3935,lorb=1,norb
        do 3935,jorb=1,norb
3935    ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
!        do 3985,j=1,norb
!        do 3985,i=1,norb
!        ovrlp(i,j,3)=0.d0
!        do 3985,l=1,norb
!3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
! BLAS:
        call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

! new eigenvectors
      allocate(tpsi(nvctrp,norb))
!   tpsi(i,iorb)=psi(i,jorb)*ovrlp(jorb,iorb,3)
      call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi,nvctrp,ovrlp(1,1,3),norb,0.d0,tpsi,nvctrp)
      call DCOPY(nvctrp*norb,tpsi,1,psi,1)
      deallocate(tpsi)

        deallocate(ovrlp,evall)

 endif

        END SUBROUTINE


        subroutine checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
        implicit real*8 (a-h,o-z)
        dimension psit(nvctrp,norbp*nproc)
        real*8, allocatable :: ovrlp(:,:,:)
        include 'mpif.h'

        allocate(ovrlp(norb,norb,2))

     do 100,iorb=1,norb
        do 100,jorb=1,norb
        ovrlp(iorb,jorb,2)=ddot(nvctrp,psit(1,iorb),1,psit(1,jorb),1)
100    continue

     call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        toler=1.d-10
        dev=0.d0
     do 110,iorb=1,norb
        do 110,jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
        dev=dev+(scpr-1.d0)**2
        else
        dev=dev+scpr**2
        endif
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
110     continue

        if (dev.gt.1.d-10) write(*,*) 'Deviation from orthogonality ',iproc,dev

        deallocate(ovrlp)

        return
        END SUBROUTINE


        subroutine checkortho(norb,norbp,nvctrp,psi)
        implicit real*8 (a-h,o-z)
        dimension psi(nvctrp,norbp)
        real*8, allocatable :: ovrlp(:,:,:)

        allocate(ovrlp(norb,norb,1))

     do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,1)=ddot(nvctrp,psi(1,iorb),1,psi(1,jorb),1)
     enddo ; enddo

        toler=1.d-10
        dev=0.d0
     do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
        dev=dev+(scpr-1.d0)**2
        else
        dev=dev+scpr**2
        endif
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
     enddo  ; enddo

        if (dev.gt.1.d-10) write(*,*) 'Deviation from orthogonality ',0,dev

        deallocate(ovrlp)


        return
        END SUBROUTINE

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine readAtomicOrbitals(iproc, filename, ngx, npsp, xp, psiat, occupat, ng, &
     & ns, np, pspatomnames, norbe, atomnames, iatype, nat)
  character(len = *), intent(in) :: filename
  integer, intent(in) :: ngx, npsp, iproc
  real*8, intent(out) :: xp(ngx, npsp), psiat(ngx, 3, npsp), occupat(3, npsp)
  integer, intent(out) :: ng(npsp), ns(npsp), np(npsp)
  character(len = 20), intent(out) :: pspatomnames(npsp)
  character(len = 20), intent(in) :: atomnames(100)
  integer, intent(out) :: norbe
  integer, intent(in) :: nat
  integer, intent(in) :: iatype(nat)

  integer :: ity, iat, i, j, ipsp

  ! Read the data file.
  ns(:) = 0
  np(:) = 0
  ng(:) = 0
  xp(:, :) = 0.d0
  psiat(:, :, :) = 0.d0
  open(unit=24,file=trim(filename),form='formatted',status='unknown')
  do ity=1,npsp
33   format(30(e12.5))
19   format(a)
     read(24,19) pspatomnames(ity)
     read(24,*) ns(ity),(occupat(i,ity),i=1,ns(ity)),  &
          np(ity),(occupat(i,ity),i=ns(ity)+1,ns(ity)+np(ity))
     if (ns(ity)+np(ity).gt.3) stop 'error ns+np'
     read(24,*) ng(ity)
     if (ng(ity).gt.ngx) stop 'enlarge ngx'
     read(24,33) (xp(i,ity)  ,i=1,ng(ity))
     do i=1,ng(ity) 
        read(24,*) (psiat(i,j,ity),j=1,ns(ity)+np(ity))
     enddo
  enddo
  close(unit=24)

  ! number of orbitals 
  norbe=0
  ipsp = 1
  do iat=1,nat
     ity=iatype(iat)
     do i=1,npsp
        if (pspatomnames(i).eq.atomnames(ity)) then
           ipsp=i
           goto 333
        endif
     enddo
     if (iproc.eq.0) write(*,*) 'no PSP for ',atomnames(ity)
     stop 
333  continue
     if (iproc.eq.0) write(*,*) 'found input wavefunction data for ',iat,atomnames(ity)
     norbe=norbe+ns(ipsp)+3*np(ipsp)
     !        if (iproc.eq.0) write(*,*) 'ns(ipsp),np(ipsp)',ns(ipsp),np(ipsp)
  enddo
  if (iproc.eq.0) write(*,*) 'number of orbitals used in the construction of input guess ',norbe

END SUBROUTINE 

subroutine createAtomicOrbitals(iproc, nproc, npsp, pspatomnames, atomnames, nat, rxyz, &
     & norbe, norbep, occupe, occupat, ngx, xp, psiat, psiatn, ng, ns, np, nvctr_c, nvctr_f, &
     & n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, nseg_c, nseg_f, &
     & keyg, keyv, iatype, psi, eks)

  implicit none

  integer, intent(in) :: npsp, nat, norbe, norbep, ngx, iproc, nproc
  integer, intent(in) :: nvctr_c, nvctr_f, n1, n2, n3, nseg_c, nseg_f
  integer, intent(in) :: nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
  integer, intent(in) :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
  integer, intent(in) :: iatype(nat)
  real*8, intent(in) :: hgrid
  real*8, intent(out) :: eks
  character(len = 20), intent(in) :: pspatomnames(npsp)
  character(len = 20), intent(in) :: atomnames(100)
  integer, intent(in) :: ng(npsp), ns(npsp), np(npsp)
  real*8, intent(in) :: rxyz(3, nat)
  real*8, intent(in) :: xp(ngx, npsp), psiat(ngx, 3, npsp)
  real*8, intent(in) :: occupat(3, npsp), psiatn(ngx)
  real*8, intent(out) :: psi(nvctr_c + 7 * nvctr_f, norbep), occupe(norbe)

  integer :: iorb, jorb, iat, ity, ipsp, i
  real*8 :: rx, ry, rz, ek, scpr

  
  eks=0.d0
  iorb=0
  ipsp = 1
  do iat=1,nat

     rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)

     ity=iatype(iat)
     do i=1,npsp
        if (pspatomnames(i).eq.atomnames(ity)) then
           ipsp=i
           goto 444
        endif
     enddo
     if (iproc.eq.0) write(*,*) 'no PSP for ',atomnames(ity)
     stop 
444  continue

     ! H:
     if (ns(ipsp).eq.1 .and. np(ipsp).eq.0) then
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        if (iorb.gt.norbe) stop 'transgpw occupe'
        occupe(iorb)=occupat(1,ipsp)
        call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'H, s, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif

        ! Li,Be 
     else if (ns(ipsp).eq.2 .and. np(ipsp).eq.0) then
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        if (iorb+1.gt.norbe) stop 'transgpw occupe'
 !      Semicore s
        occupe(iorb)=occupat(1,ipsp)
        call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'Li, Be, semicore sATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !    valence s
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(2,ipsp)
        call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,2,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) ; scpr=1.d0/sqrt(scpr)
           write(*,'(a50,i4,1x,e14.7)') 'Li, Be, s , ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr 
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif

        ! Na,Mg
     else if (ns(ipsp).eq.2 .and. np(ipsp).eq.1) then
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        if (iorb+1.gt.norbe) stop 'transgpw occupe'
!       Semicore s
        occupe(iorb)=occupat(1,ipsp)
        call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'Na,Mg, semicore s, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !    valence s
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(2,ipsp)
        call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,2,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) ; scpr=1.d0/sqrt(scpr)
           write(*,'(a50,i4,1x,e14.7)') 'Na,Mg, s, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr 
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !    semicore px
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(3,ipsp)*(1.d0/3.d0)
        call atomkin(1,ng(ipsp),xp(1,ipsp),psiat(1,3,ipsp),psiatn,ek) ; eks=eks+ek*3.d0*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),1,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'Na,Mg, semicore px, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !    semicore py
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(3,ipsp)*(1.d0/3.d0)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,1,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'Na,Mg, semicore py, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !    semicore pz
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(3,ipsp)*(1.d0/3.d0)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,1,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'Na,Mg, semicore pz, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif

        ! B,C,N,O,F , Al,Si,P,S,Cl
     else if (ns(ipsp).eq.1 .and. np(ipsp).eq.1) then
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        if (iorb+3.gt.norbe) stop 'transgpw occupe'
        !    Valence s
        occupe(iorb)=occupat(1,ipsp)
        call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'valence s, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !   valence px
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(2,ipsp)*(1.d0/3.d0)
        call atomkin(1,ng(ipsp),xp(1,ipsp),psiat(1,2,ipsp),psiatn,ek) ; eks=eks+3.d0*ek*occupe(iorb)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),1,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'valence px, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !   valence py
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(2,ipsp)*(1.d0/3.d0)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,1,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'valence py, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif
        !   valence pz
        iorb=iorb+1
        jorb=iorb-iproc*norbep
        occupe(iorb)=occupat(2,ipsp)*(1.d0/3.d0)
        if (myorbital(iorb,norbe,iproc,nproc)) then
           call crtonewave(n1,n2,n3,ng(ipsp),0,0,1,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
                0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                psi(1,jorb),psi(nvctr_c+1,jorb))
           call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
           write(*,'(a50,i4,1x,e14.7)') 'valence pz, ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
           call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
        endif

        ! ERROR
     else 
        stop 'error in inguess.dat'
     endif
  end do

END SUBROUTINE 


subroutine input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nat,norb,norbp,n1,n2,n3,nfft1,nfft2,nfft3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
     rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
     atomnames,ntypes,iatype,pkernel,psppar,npspcode,ppsi,eval,accurex)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors writes its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  implicit real*8 (a-h,o-z)
  logical parallel
  character*20 atomnames(100)
  character*20 pspatomnames(15)
  parameter(eps_mach=1.d-12)
  parameter (ngx=31)
  parameter (npsp=15)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension  xp(ngx,15),psiat(ngx,3,15),occupat(3,15),ng(15),ns(15),np(15),psiatn(ngx)
  dimension rxyz(3,nat),iatype(nat),eval(norb)
  dimension rhopot((2*n1+31)*(2*n2+31)*(2*n3+31)),pot_ion((2*n1+31)*(2*n2+31)*(2*n3+31))
  dimension pkernel(*)
  dimension psppar(0:4,0:4,ntypes),npspcode(ntypes)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
  dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
  dimension proj(nprojel)
  dimension ppsi(nvctr_c+7*nvctr_f,norbp)

  allocatable :: psi(:,:),hpsi(:,:),psit(:,:),hpsit(:,:),ppsit(:,:),occupe(:)
  allocatable :: hamovr(:,:,:),evale(:),work_lp(:)
  include 'mpif.h'

  ! Read the inguess.dat file.
  call readAtomicOrbitals(iproc, 'inguess.dat', ngx, npsp, xp, psiat, occupat, ng, &
       & ns, np, pspatomnames, norbe, atomnames, iatype, nat)

  !  allocate wavefunctions and their occupation numbers
  allocate(occupe(norbe))
  tt=dble(norbe)/dble(nproc)
  norbep=int((1.d0-eps_mach*tt) + tt)
  if (iproc.eq.0) write(79,'(a40,i10)') 'words for (h)psi inguess',2*(nvctr_c+7*nvctr_f)*norbep
  allocate(psi(nvctr_c+7*nvctr_f,norbep))
  allocate(hpsi(nvctr_c+7*nvctr_f,norbep))
  if (iproc.eq.0) write(79,*) 'allocation done'
  norbeme=max(min((iproc+1)*norbep,norbe)-iproc*norbep,0)
  write(*,*) 'iproc ',iproc,' treats ',norbeme,' inguess orbitals '

  hgridh=.5d0*hgrid

  ! Create input guess orbitals
  call createAtomicOrbitals(iproc, nproc, npsp, pspatomnames, atomnames, nat, rxyz, &
       & norbe, norbep, occupe, occupat, ngx, xp, psiat, psiatn, ng, ns, np, nvctr_c, nvctr_f, &
       & n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, nseg_c, nseg_f, &
       & keyg, keyv, iatype, psi, eks)

! resulting charge density and potential
       call sumrho(parallel,iproc,norbe,norbep,n1,n2,n3,hgrid,occupe,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)
      ixc=1   ! LDA functional
      if (parallel) then
          call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,ixc, &
               pot_ion,rhopot,ehart,eexcu,vexcu,iproc,nproc)
       else
          call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,ixc, &
               pot_ion,rhopot,ehart,eexcu,vexcu)
       end if

! set up subspace Hamiltonian 
        allocate(hamovr(norbe,norbe,4))

        call applylocpotkinall(iproc,norbe,norbep,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
                 hgrid,occupe,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
                 psi,rhopot,hpsi,epot_sum,ekin_sum)

       if (parallel) then
       tt=ekin_sum
       call MPI_ALLREDUCE(tt,ekin_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       endif

       accurex=abs(eks-ekin_sum)
       write(*,*) 'ekin_sum,eks',ekin_sum,eks

        call applyprojectorsall(iproc,ntypes,nat,iatype,psppar,npspcode,occupe, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norbe,norbep,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)

 if (parallel) then
        if (iproc.eq.0) write(79,'(a40,i10)') 'words for psit inguess',nvctrp*norbep*nproc
        allocate(psit(nvctrp,norbep*nproc))
        if (iproc.eq.0) write(79,*) 'allocation done'

        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,psi,psit)

        deallocate(psi)

        if (iproc.eq.0) write(79,'(a40,i10)') 'words for hpsit inguess',2*nvctrp*norbep*nproc
        allocate(hpsit(nvctrp,norbep*nproc))
        if (iproc.eq.0) write(79,*) 'allocation done'

        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)

        deallocate(hpsi)

!       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
!       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,hamovr(1,1,3),norbe)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp, psit,nvctrp,0.d0,hamovr(1,1,4),norbe)
        deallocate(hpsit)

        call MPI_ALLREDUCE (hamovr(1,1,3),hamovr(1,1,1),2*norbe**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

! calculate  KS orbitals
!      if (iproc.eq.0) then
!        write(*,*) 'KS Hamiltonian',iproc
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
!        enddo
!        write(*,*) 'Overlap',iproc
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
!        enddo
!     endif

        n_lp=5000
        allocate(work_lp(n_lp),evale(norbe))
        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
        if (iproc.eq.0) then
        do iorb=1,norbe
        write(*,*) 'evale(',iorb,')=',evale(iorb)
        enddo
        endif
        do iorb=1,norb
        eval(iorb)=evale(iorb)
        enddo
        deallocate(work_lp,evale)

        if (iproc.eq.0) write(79,'(a40,i10)') 'words for ppsit ',nvctrp*norbp*nproc
        allocate(ppsit(nvctrp,norbp*nproc))
        if (iproc.eq.0) write(79,*) 'allocation done'

! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psit,nvctrp,hamovr,norbe,0.d0,ppsit,nvctrp)
       call  untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,ppsit,ppsi)
        deallocate(psit,ppsit)

        if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  else !serial case
!       hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,hamovr(1,1,1),norbe)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp, psi,nvctrp,0.d0,hamovr(1,1,2),norbe)
        deallocate(hpsi)

! calculate  KS orbitals
!        write(*,*) 'KS Hamiltonian'
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
!        enddo
!        write(*,*) 'Overlap'
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
!        enddo

        n_lp=5000
        allocate(work_lp(n_lp),evale(norbe))
        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
        if (iproc.eq.0) then
        do iorb=1,norbe
        write(*,*) 'evale(',iorb,')=',evale(iorb)
        enddo
        endif
        do iorb=1,norb
        eval(iorb)=evale(iorb)
        enddo
        deallocate(work_lp,evale)

! ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
        call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psi,nvctrp,hamovr,norbe,0.d0,ppsi,nvctrp)
        deallocate(psi)

  endif

        deallocate(hamovr,occupe)

    return
    END SUBROUTINE

subroutine solveKS(parallel, iproc, nproc, norb, norbp, norbe, norbep, nvctr_c, &
  & nvctr_f, nvctrp, psi, hpsi, ppsi, eval)

implicit real*8 (a-h,o-z)

logical, intent(in) :: parallel
integer, intent(in) :: nvctrp, norbe, norbep, nproc, iproc
integer, intent(in) :: nvctr_c, nvctr_f, norb, norbp
dimension psi(nvctr_c+7*nvctr_f,norbep), hpsi(nvctr_c+7*nvctr_f,norbep)
real*8, intent(out) :: ppsi(nvctr_c+7*nvctr_f,norbp), eval(norb)
        include 'mpif.h'

allocatable :: ppsit(:,:), psit(:,:), hpsit(:,:), hamovr(:,:,:),work_lp(:),evale(:)

        allocate(hamovr(norbe,norbe,4))
 if (parallel) then
        write(79,'(a40,i10)') 'words for psit inguess',nvctrp*norbep*nproc
        allocate(psit(nvctrp,norbep*nproc))
        write(79,*) 'allocation done'

        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,psi,psit)

        write(79,'(a40,i10)') 'words for hpsit inguess',2*nvctrp*norbep*nproc
        allocate(hpsit(nvctrp,norbep*nproc))
        write(79,*) 'allocation done'

        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)

!       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
!       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,hamovr(1,1,3),norbe)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp, psit,nvctrp,0.d0,hamovr(1,1,4),norbe)
        deallocate(hpsit)

        call MPI_ALLREDUCE (hamovr(1,1,3),hamovr(1,1,1),2*norbe**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

! calculate  KS orbitals
!      if (iproc.eq.0) then
!        write(*,*) 'KS Hamiltonian',iproc
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
!        enddo
!        write(*,*) 'Overlap',iproc
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
!        enddo
!     endif

        n_lp=5000
        allocate(work_lp(n_lp),evale(norbe))
        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
        if (iproc.eq.0) then
        do iorb=1,norbe
        write(*,*) 'evale(',iorb,')=',evale(iorb)
        enddo
        endif
        eval(1:norb) = evale(1:norb)
        deallocate(work_lp,evale)

        write(79,'(a40,i10)') 'words for ppsit ',nvctrp*norbep*nproc
        allocate(ppsit(nvctrp,norbp*nproc))
        write(79,*) 'allocation done'

! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psit,nvctrp,hamovr,norbe,0.d0,ppsit,nvctrp)

       call  untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,ppsit,ppsi)

        deallocate(psit,ppsit)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  else !serial case
!       hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,hamovr(1,1,1),norbe)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp, psi,nvctrp,0.d0,hamovr(1,1,2),norbe)

! calculate  KS orbitals
!        write(*,*) 'KS Hamiltonian'
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
!        enddo
!        write(*,*) 'Overlap'
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
!        enddo

        n_lp=5000
        allocate(work_lp(n_lp),evale(norbe))
        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
        if (iproc.eq.0) then
        do iorb=1,norbe
        write(*,*) 'evale(',iorb,')=',evale(iorb)
        enddo
        endif
        eval(1:norb) = evale(1:norb)
        deallocate(work_lp,evale)

! ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
        call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psi,nvctrp,hamovr,norbe,0.d0,ppsi,nvctrp)

  endif
  deallocate(hamovr)
END SUBROUTINE


        logical function myorbital(iorb,norbe,iproc,nproc)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)

        tt=dble(norbe)/dble(nproc)
        norbep=int((1.d0-eps_mach*tt) + tt)
        if (iorb .ge. iproc*norbep+1 .and. iorb .le. min((iproc+1)*norbep,norbe)) then
        myorbital=.true.
        else
        myorbital=.false.
        endif

        return
        end function


    subroutine KStrans_p(iproc,nproc,norb,norbp,nvctrp,occup,  & 
                           hpsit,psit,evsum,eval)
! at the start each processor has all the Psi's but only its part of the HPsi's
! at the end each processor has only its part of the Psi's
        implicit real*8 (a-h,o-z)
        dimension occup(norb),eval(norb)
        dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc)
! arrays for KS orbitals
        allocatable :: hamks(:,:,:),work_lp(:),psitt(:,:)
        include 'mpif.h'

! set up Hamiltonian matrix
        allocate(hamks(norb,norb,2))
        do jorb=1,norb
        do iorb=1,norb
        hamks(iorb,jorb,2)=0.d0
    enddo
    enddo
        do iorb=1,norb
        do jorb=1,norb
        scpr=ddot(nvctrp,psit(1,jorb),1,hpsit(1,iorb),1)
        hamks(iorb,jorb,2)=scpr
        enddo
        enddo


        call MPI_ALLREDUCE (hamks(1,1,2),hamks(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        
!        write(*,*) 'KS Hamiltonian',iproc
!        do iorb=1,norb
!        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
!        enddo

        n_lp=1000
        allocate(work_lp(n_lp))
        call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
        evsum=0.d0
        do iorb=1,norb
        evsum=evsum+eval(iorb)*occup(iorb)
        if (iproc.eq.0) write(*,*) 'eval(',iorb,')=',eval(iorb)
        enddo
        deallocate(work_lp)
        if (info.ne.0) write(*,*) 'DSYEV ERROR',info

        allocate(psitt(nvctrp,norbp*nproc))
! Transform to KS orbitals
      do iorb=1,norb
        call razero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        enddo
      enddo
        deallocate(hamks)

        call DCOPY(nvctrp*norbp*nproc,psitt,1,psit,1)
        deallocate(psitt)

    return
        END SUBROUTINE



    subroutine KStrans(norb,norbp,nvctrp,occup,hpsi,psi,evsum,eval)
! at the start each processor has all the Psi's but only its part of the HPsi's
! at the end each processor has only its part of the Psi's
        implicit real*8 (a-h,o-z)
        dimension occup(norb),eval(norb)
        dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp)
! arrays for KS orbitals
        allocatable :: hamks(:,:,:),work_lp(:),psitt(:,:)

! set up Hamiltonian matrix
        allocate(hamks(norb,norb,2))
        do jorb=1,norb
        do iorb=1,norb
        hamks(iorb,jorb,2)=0.d0
    enddo
    enddo
        do iorb=1,norb
        do jorb=1,norb
        scpr=ddot(nvctrp,psi(1,jorb),1,hpsi(1,iorb),1)
        hamks(iorb,jorb,1)=scpr
        enddo
        enddo


!        write(*,*) 'KS Hamiltonian',0
!        do iorb=1,norb
!        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
!        enddo

        n_lp=1000
        allocate(work_lp(n_lp))
        call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
        evsum=0.d0
        do iorb=1,norb
        evsum=evsum+eval(iorb)*occup(iorb)
        write(*,*) 'eval(',iorb,')=',eval(iorb)
        enddo
        deallocate(work_lp)
        if (info.ne.0) write(*,*) 'DSYEV ERROR',info

        allocate(psitt(nvctrp,norbp))
! Transform to KS orbitals
      do iorb=1,norb
        call razero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psi(1,jorb),1,psitt(1,iorb),1)
        enddo
      enddo
        deallocate(hamks)

        call DCOPY(nvctrp*norbp,psitt,1,psi,1)
        deallocate(psitt)

    return
        END SUBROUTINE



        subroutine crtonewave(n1,n2,n3,nterm,lx,ly,lz,xp,psiat,rx,ry,rz,hgrid, & 
                   nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                   nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f)
! returns an input guess orbital that is a Gaussian centered at a Wannier center
! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays psi_c, psi_f
        implicit real*8 (a-h,o-z)
        parameter(nw=16000)
        dimension xp(nterm),psiat(nterm)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work
        real*8, allocatable :: psig_c(:,:,:), psig_f(:,:,:,:)

        allocate(wprojx(0:n1,2),wprojy(0:n2,2),wprojz(0:n3,2),work(0:nw,2))
        allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c))
        allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f))

      iterm=1
        gau_a=xp(iterm)
        n_gau=lx
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        n_gau=ly
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        n_gau=lz
        CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

! First term: coarse projector components
          do i3=nl3_c,nu3_c
          do i2=nl2_c,nu2_c
          do i1=nl1_c,nu1_c
            psig_c(i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
          enddo ; enddo ; enddo

! First term: fine projector components
          do i3=nl3_f,nu3_f
          do i2=nl2_f,nu2_f
          do i1=nl1_f,nu1_f
            psig_f(1,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
            psig_f(2,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(3,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(4,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(5,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(6,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
            psig_f(7,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
          enddo ; enddo ; enddo

      do 100,iterm=2,nterm
        gau_a=xp(iterm)
        n_gau=lx
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        n_gau=ly
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        n_gau=lz
        CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

! First term: coarse projector components
          do i3=nl3_c,nu3_c
          do i2=nl2_c,nu2_c
          do i1=nl1_c,nu1_c
            psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
          enddo ; enddo ; enddo

! First term: fine projector components
          do i3=nl3_f,nu3_f
          do i2=nl2_f,nu2_f
          do i1=nl1_f,nu1_f
            psig_f(1,i1,i2,i3)=psig_f(1,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
            psig_f(2,i1,i2,i3)=psig_f(2,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(3,i1,i2,i3)=psig_f(3,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(4,i1,i2,i3)=psig_f(4,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(5,i1,i2,i3)=psig_f(5,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(6,i1,i2,i3)=psig_f(6,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
            psig_f(7,i1,i2,i3)=psig_f(7,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
          enddo ; enddo ; enddo

100    continue

! coarse part
    do iseg=1,nseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
      do i=i0,i1
            psi_c(i-i0+jj)=psig_c(i,i2,i3)
          enddo
        enddo

! fine part
    do iseg=1,nseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
      do i=i0,i1
            psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
            psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
            psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
            psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
            psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
            psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
            psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
          enddo
        enddo
  
          deallocate(wprojx,wprojy,wprojz,work,psig_c,psig_f)

    return
    END SUBROUTINE
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 subroutine reformatmywaves(iproc, norb, norbp, nat, &
      & hgrid_old, nvctr_c_old, nvctr_f_old, n1_old, n2_old, n3_old, rxyz_old, &
      & nseg_c_old, nseg_f_old, keyg_old, keyv_old, psi_old, &
      & hgrid, nvctr_c, nvctr_f, n1, n2, n3, rxyz, &
      & nseg_c, nseg_f, keyg, keyv, psi)
   implicit real*8 (a-h,o-z)
   dimension :: rxyz(3,nat), rxyz_old(3,nat), center(3), center_old(3)
   dimension :: keyg_old(2, nseg_c_old + nseg_f_old), keyv_old(nseg_c_old + nseg_f_old)
   dimension :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
   dimension :: psi_old(nvctr_c_old + 7 * nvctr_f_old, norbp), psi(nvctr_c + 7 * nvctr_f, norbp)

   allocatable :: psifscf(:,:,:), psigold(:,:,:,:,:,:)

   allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8))

! calculate center of molecule
        c1=0.d0 ; c2=0.d0 ; c3=0.d0
        do iat=1,nat
        c1=c1+rxyz(1,iat) ; c2=c2+rxyz(2,iat) ; c3=c3+rxyz(3,iat)
        enddo
        center(1)=c1/nat ; center(2)=c2/nat ; center(3)=c3/nat
        c1=0.d0 ; c2=0.d0 ; c3=0.d0
        do iat=1,nat
        c1=c1+rxyz_old(1,iat) ; c2=c2+rxyz_old(2,iat) ; c3=c3+rxyz_old(3,iat)
        enddo
        center_old(1)=c1/nat ; center_old(2)=c2/nat ; center_old(3)=c3/nat


   do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

      if (hgrid_old.eq. hgrid .and. nvctr_c_old.eq.nvctr_c .and. nvctr_f_old.eq.nvctr_f  & 
           .and. n1_old.eq.n1  .and. n2_old.eq.n2 .and. n3_old.eq.n3  .and.  &
             abs(center(1)-center_old(1)).lt.1.d-3 .and. &
             abs(center(2)-center_old(2)).lt.1.d-3 .and. &
             abs(center(3)-center_old(3)).lt.1.d-3  ) then


         write(*,*) 'wavefunction ',iorb,' needs NO reformatting on processor',iproc
         do j=1,nvctr_c_old
            psi(j,iorb-iproc*norbp)=psi_old(j, iorb - iproc * norbp)
         enddo
         do j=1,7*nvctr_f_old-6,7
            psi(nvctr_c+j+0,iorb-iproc*norbp)=psi_old(nvctr_c+j+0,iorb-iproc*norbp)
            psi(nvctr_c+j+1,iorb-iproc*norbp)=psi_old(nvctr_c+j+1,iorb-iproc*norbp)
            psi(nvctr_c+j+2,iorb-iproc*norbp)=psi_old(nvctr_c+j+2,iorb-iproc*norbp)
            psi(nvctr_c+j+3,iorb-iproc*norbp)=psi_old(nvctr_c+j+3,iorb-iproc*norbp)
            psi(nvctr_c+j+4,iorb-iproc*norbp)=psi_old(nvctr_c+j+4,iorb-iproc*norbp)
            psi(nvctr_c+j+5,iorb-iproc*norbp)=psi_old(nvctr_c+j+5,iorb-iproc*norbp)
            psi(nvctr_c+j+6,iorb-iproc*norbp)=psi_old(nvctr_c+j+6,iorb-iproc*norbp)
         enddo

      else
         write(*,*) 'wavefunction ',iorb,' needs reformatting on processor',iproc
         if (hgrid_old.ne.hgrid) then 
           write(*,*) 'because hgrid_old >< hgrid',hgrid_old, hgrid
         else if (nvctr_c_old.ne.nvctr_c) then
           write(*,*) 'because nvctr_c_old >< nvctr_c',nvctr_c_old,nvctr_c
         else if (nvctr_f_old.ne.nvctr_f)  then
           write(*,*) 'because nvctr_f_old >< nvctr_f',nvctr_f_old,nvctr_f
         else if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 )  then  
           write(*,*) 'because cell size has changed',n1_old,n1  , n2_old,n2 , n3_old,n3
         else
           write(*,*) 'molecule was shifted' 
         endif

         allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2))

         call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)


         ! coarse part
         do iseg=1,nseg_c_old
            jj=keyv_old(iseg)
            j0=keyg_old(1,iseg)
            j1=keyg_old(2,iseg)
            ii=j0-1
            i3=ii/((n1_old+1)*(n2_old+1))
            ii=ii-i3*(n1_old+1)*(n2_old+1)
            i2=ii/(n1_old+1)
            i0=ii-i2*(n1_old+1)
            i1=i0+j1-j0
            do i=i0,i1
               psigold(i, 1, i2, 1, i3, 1) = psi_old(i - i0 + jj, iorb - iproc * norbp)
            enddo
         enddo

         ! fine part
         do iseg=1,nseg_f_old
            jj=keyv_old(nseg_c_old + iseg)
            j0=keyg_old(1,nseg_c_old + iseg)
            j1=keyg_old(2,nseg_c_old + iseg)
            ii=j0-1
            i3=ii/((n1_old+1)*(n2_old+1))
            ii=ii-i3*(n1_old+1)*(n2_old+1)
            i2=ii/(n1_old+1)
            i0=ii-i2*(n1_old+1)
            i1=i0+j1-j0
            do i=i0,i1
               psigold(i,2,i2,1,i3,1)=psi_old(nvctr_c_old + 1 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
               psigold(i,1,i2,2,i3,1)=psi_old(nvctr_c_old + 2 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
               psigold(i,2,i2,2,i3,1)=psi_old(nvctr_c_old + 3 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
               psigold(i,1,i2,1,i3,2)=psi_old(nvctr_c_old + 4 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
               psigold(i,2,i2,1,i3,2)=psi_old(nvctr_c_old + 5 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
               psigold(i,1,i2,2,i3,2)=psi_old(nvctr_c_old + 6 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
               psigold(i,2,i2,2,i3,2)=psi_old(nvctr_c_old + 7 + 7 * (i-i0+jj - 1), iorb - iproc * norbp)
            enddo
         enddo

         call reformatonewave(iproc, hgrid_old, &
              & n1_old, n2_old, n3_old, center_old, psigold, hgrid, &
              & nvctr_c, nvctr_f, n1, n2, n3, center, nseg_c, nseg_f, keyg, keyv, psifscf, & 
              & psi(1,iorb - iproc * norbp))

         deallocate(psigold)
      end if
   end do
   
   deallocate(psifscf)  

 END SUBROUTINE 

 subroutine reformatonewave(iproc, hgrid_old, n1_old, n2_old, n3_old, &
      & center_old, psigold, hgrid, nvctr_c, nvctr_f, n1, n2, n3, center, nseg_c, nseg_f, &
      & keyg, keyv, psifscf, psi)
   implicit real*8 (a-h,o-z)
   logical cif1,cif2,cif3
   dimension xya(-1:1,-1:1),xa(-1:1)
   dimension :: center(3), center_old(3)
   dimension :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
   dimension :: psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2), psi(nvctr_c + 7 * nvctr_f)
   dimension :: psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

   allocatable :: psifscfold(:,:,:),psifscfoex(:,:,:),psig(:,:,:,:,:,:),ww(:)

   allocate(psifscfold(-7:2*n1_old+8,-7:2*n2_old+8,-7:2*n3_old+8), &
        psifscfoex(-8:2*n1_old+9,-8:2*n2_old+9,-8:2*n3_old+9))

   ! calculate fine scaling functions, psifscfoex=wwold((2*n1_old+16)*(2*n2_old+16)*(2*n3_old+16))
   call synthese_grow(n1_old,n2_old,n3_old,psifscfoex,psigold,psifscfold) 

   do i3=-7,2*n3_old+8
      do i2=-7,2*n2_old+8
         i1=-8
         psifscfoex(i1,i2,i3)=0.d0
         do i1=-7,2*n1_old+8
            psifscfoex(i1,i2,i3)=psifscfold(i1,i2,i3)
         enddo
         i1=2*n1_old+9
         psifscfoex(i1,i2,i3)=0.d0
      enddo
   enddo

   i3=-8
   do i2=-8,2*n2_old+9
      do i1=-8,2*n1_old+9
         psifscfoex(i1,i2,i3)=0.d0
      enddo
   enddo
   i3=2*n3_old+9
   do i2=-8,2*n2_old+9
      do i1=-8,2*n1_old+9
         psifscfoex(i1,i2,i3)=0.d0
      enddo
   enddo

   i2=-8
   do i3=-8,2*n3_old+9
      do i1=-8,2*n1_old+9
         psifscfoex(i1,i2,i3)=0.d0
      enddo
   enddo
   i2=2*n2_old+9
   do i3=-8,2*n3_old+9
      do i1=-8,2*n1_old+9
         psifscfoex(i1,i2,i3)=0.d0
      enddo
   enddo

   ! transform to new structure    
   dx=center(1)-center_old(1)
   dy=center(2)-center_old(2)
   dz=center(3)-center_old(3)
!   write(*,*) 'dxyz',dx,dy,dz
   hgridh=.5d0*hgrid
   hgridh_old=.5d0*hgrid_old
   call razero((2*n1+16)*(2*n2+16)*(2*n3+16),psifscf)
   do i3=-7,2*n3+8
      z=i3*hgridh
      j3=nint((z-dz)/hgridh_old)
      cif3=(j3.ge.-7 .and. j3.le.2*n3_old+8)
      do i2=-7,2*n2+8
         y=i2*hgridh
         j2=nint((y-dy)/hgridh_old)
         cif2=(j2.ge.-7 .and. j2.le.2*n2_old+8)
         do i1=-7,2*n1+8
            x=i1*hgridh
            j1=nint((x-dx)/hgridh_old)
            cif1=(j1.ge.-7 .and. j1.le.2*n1_old+8)

            !        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfold(j1,j2,j3)
            !        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfoex(j1,j2,j3)

            if (cif1 .and. cif2 .and. cif3) then 
               zr = ((z-dz)-j3*hgridh_old)/hgridh_old
               do l2=-1,1
                  do l1=-1,1
                     ym1=psifscfoex(j1+l1,j2+l2,j3-1)
                     y00=psifscfoex(j1+l1,j2+l2,j3  )
                     yp1=psifscfoex(j1+l1,j2+l2,j3+1)
                     xya(l1,l2)=ym1 + (1.d0 + zr)*(y00 - ym1 + zr*(.5d0*ym1 - y00  + .5d0*yp1))
                  enddo
               enddo

               yr = ((y-dy)-j2*hgridh_old)/hgridh_old
               do l1=-1,1
                  ym1=xya(l1,-1)
                  y00=xya(l1,0)
                  yp1=xya(l1,1)
                  xa(l1)=ym1 + (1.d0 + yr)*(y00 - ym1 + yr*(.5d0*ym1 - y00  + .5d0*yp1))
               enddo

               xr = ((x-dx)-j1*hgridh_old)/hgridh_old
               ym1=xa(-1)
               y00=xa(0)
               yp1=xa(1)
               psifscf(i1,i2,i3)=ym1 + (1.d0 + xr)*(y00 - ym1 + xr*(.5d0*ym1 - y00  + .5d0*yp1))

            endif

         enddo
      enddo
   enddo

   deallocate(psifscfold,psifscfoex)
   allocate(psig(0:n1,2,0:n2,2,0:n3,2),ww((2*n1+16)*(2*n2+16)*(2*n3+16)))

        call analyse_shrink(n1,n2,n3,ww,psifscf,psig)
        call compress(n1,n2,n3,0,n1,0,n2,0,n3,  &
                    nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psig,psi(1),psi(nvctr_c+1))

   deallocate(psig,ww)
 END SUBROUTINE 




 subroutine readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,  & 
      nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
   ! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
   implicit real*8 (a-h,o-z)
   character*50 filename
   character*4 f4
   dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
   dimension psi(nvctr_c+7*nvctr_f,norbp)
   dimension rxyz(3,nat),eval(norb),center(3)
   allocatable :: psifscf(:,:,:)

   call cpu_time(tr0)
   call system_clock(ncount1,ncount_rate,ncount_max)

   allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8))

! calculate center of molecule
        c1=0.d0 ; c2=0.d0 ; c3=0.d0
        do iat=1,nat
        c1=c1+rxyz(1,iat) ; c2=c2+rxyz(2,iat) ; c3=c3+rxyz(3,iat)
        enddo
        center(1)=c1/nat ; center(2)=c2/nat ; center(3)=c3/nat

   do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

      write(f4,'(i4.4)') iorb
      filename = 'wavefunction.'//f4
      open(unit=99,file=filename,status='unknown')

      call readonewave(99, .true., iorb,iproc,n1,n2,n3, &
           & hgrid,center,nseg_c,nseg_f, nvctr_c,nvctr_f,keyg,keyv,psi(1,iorb-iproc*norbp),eval(iorb),psifscf)

      close(99)
   end do

   deallocate(psifscf)

   call cpu_time(tr1)
   call system_clock(ncount2,ncount_rate,ncount_max)
   tel=dble(ncount2-ncount1)/dble(ncount_rate)
   write(*,'(a40,i4,2(1x,e10.3))') 'READING WAVES TIME',iproc,tr1-tr0,tel
 END SUBROUTINE 


 subroutine readonewave(unitwf, useFormattedInput, iorb,iproc,n1,n2,n3, &
      & hgrid,center,nseg_c,nseg_f, nvctr_c,nvctr_f,keyg,keyv,psi,eval,psifscf)
   implicit real*8 (a-h,o-z)
   dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
   dimension psi(nvctr_c+7*nvctr_f)
   dimension center(3),center_old(3)
   dimension :: psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
   allocatable :: psigold(:,:,:,:,:,:)
   integer :: unitwf
   logical :: useFormattedInput

   if (useFormattedInput) then
      read(unitwf,*) iorb_old,eval
   else
      read(unitwf) iorb_old,eval
   end if
   if (iorb_old.ne.iorb) stop 'readonewave'
   if (useFormattedInput) then
      read(unitwf,*) hgrid_old
      read(unitwf,*) n1_old,n2_old,n3_old
      read(unitwf,*) (center_old(j),j=1,3)
   else
      read(unitwf) hgrid_old
      read(unitwf) n1_old,n2_old,n3_old
      read(unitwf) (center_old(j),j=1,3)
   end if
   write(*,'(i2,6(1x,e14.7))') iproc,(center(j),j=1,3),(center_old(j),j=1,3)
   if (useFormattedInput) then
      read(unitwf,*) nvctr_c_old, nvctr_f_old
   else
      read(unitwf) nvctr_c_old, nvctr_f_old
   end if

!           write(*,*) iorb,' hgrid_old,hgrid ',hgrid_old,hgrid
!           write(*,*) iorb,' nvctr_c_old,nvctr_c ',nvctr_c_old,nvctr_c
!           write(*,*) iorb,' nvctr_f_old,nvctr_f ',nvctr_f_old,nvctr_f
!           write(*,*) iorb,' n1_old,n1 ',n1_old,n1
!           write(*,*) iorb,' n2_old,n2 ',n2_old,n2
!           write(*,*) iorb,' n3_old,n3 ',n3_old,n3

   if (hgrid_old.eq. hgrid .and. nvctr_c_old.eq.nvctr_c .and. nvctr_f_old.eq.nvctr_f  & 
        .and. n1_old.eq.n1  .and. n2_old.eq.n2 .and. n3_old.eq.n3 ) then

      write(*,*) 'wavefunction ',iorb,' needs NO reformatting on processor',iproc
      do j=1,nvctr_c_old
         if (useFormattedInput) then
            read(unitwf,*) i1,i2,i3,tt
         else
            read(unitwf) i1,i2,i3,tt
         end if
         psi(j)=tt
      enddo
      do j=1,7*nvctr_f_old-6,7
         if (useFormattedInput) then
            read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         else
            read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         end if
         psi(nvctr_c+j+0)=t1
         psi(nvctr_c+j+1)=t2
         psi(nvctr_c+j+2)=t3
         psi(nvctr_c+j+3)=t4
         psi(nvctr_c+j+4)=t5
         psi(nvctr_c+j+5)=t6
         psi(nvctr_c+j+6)=t7
      enddo

   else
      write(*,*) 'wavefunction ',iorb,' needs reformatting on processor',iproc
      if (hgrid_old.ne.hgrid) write(*,*) 'because hgrid_old >< hgrid',hgrid_old,hgrid
      if (nvctr_c_old.ne.nvctr_c) write(*,*) 'because nvctr_c_old >< nvctr_c',nvctr_c_old,nvctr_c
      if (nvctr_f_old.ne.nvctr_f) write(*,*) 'because nvctr_f_old >< nvctr_f',nvctr_f_old,nvctr_f
      if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 ) &
           write(*,*) 'because cell size has changed',n1_old,n1  , n2_old,n2 , n3_old,n3
      
      allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2))

      call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)
      do iel=1,nvctr_c_old
         if (useFormattedInput) then
            read(unitwf,*) i1,i2,i3,tt
         else
            read(unitwf) i1,i2,i3,tt
         end if
         psigold(i1,1,i2,1,i3,1)=tt
      enddo
      do iel=1,nvctr_f_old
         if (useFormattedInput) then
            read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         else
            read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         end if
         psigold(i1,2,i2,1,i3,1)=t1
         psigold(i1,1,i2,2,i3,1)=t2
         psigold(i1,2,i2,2,i3,1)=t3
         psigold(i1,1,i2,1,i3,2)=t4
         psigold(i1,2,i2,1,i3,2)=t5
         psigold(i1,1,i2,2,i3,2)=t6
         psigold(i1,2,i2,2,i3,2)=t7
      enddo

      ! I put nat = 1 here, since only one position is saved in wavefunction files.
      call reformatonewave(iproc, hgrid_old, n1_old, n2_old, n3_old, &
           & center_old, psigold, hgrid, nvctr_c, nvctr_f, n1, n2, n3, center, nseg_c, nseg_f, &
           & keyg, keyv, psifscf, psi)

      deallocate(psigold)

   endif
 END SUBROUTINE 



        subroutine writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
                   nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
! write all my wavefunctions in files by calling writeonewave
        implicit real*8 (a-h,o-z)
        character*4 f4
        character*50 filename
        dimension rxyz(3,nat),eval(norb),center(3)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp)

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)

! calculate center of molecule
        c1=0.d0 ; c2=0.d0 ; c3=0.d0
        do iat=1,nat
        c1=c1+rxyz(1,iat) ; c2=c2+rxyz(2,iat) ; c3=c3+rxyz(3,iat)
        enddo
        center(1)=c1/nat ; center(2)=c2/nat ; center(3)=c3/nat

       do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        write(f4,'(i4.4)') iorb
        filename = 'wavefunction.'//f4
        write(*,*) 'opening ',filename
        open(unit=99,file=filename,status='unknown')

       call writeonewave(99, .true., iorb,n1,n2,n3,hgrid,center,  & 
                         nseg_c,nvctr_c,keyg(1,1),keyv(1)  & 
                        ,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                        psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),norb,eval)
          close(99)

       enddo

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(*,'(a40,i4,2(1x,e10.3))') 'WRITE WAVES TIME',iproc,tr1-tr0,tel


       return
       END SUBROUTINE



        subroutine writeonewave(unitwf, useFormattedOutput, iorb,n1,n2,n3,hgrid,center,  & 
                           nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                           nseg_f,nvctr_f,keyg_f,keyv_f, & 
                              psi_c,psi_f,norb,eval)
        implicit real*8 (a-h,o-z)
        logical :: useFormattedOutput
        integer :: unitwf
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f),center(3),eval(norb)


        if (useFormattedOutput) then
         write(unitwf,*) iorb,eval(iorb)
         write(unitwf,*) hgrid
         write(unitwf,*) n1,n2,n3
         write(unitwf,'(3(1x,e24.17))') (center(j),j=1,3)
         write(unitwf,*) nvctr_c, nvctr_f
        else
         write(unitwf) iorb,eval(iorb)
         write(unitwf) hgrid
         write(unitwf) n1,n2,n3
         write(unitwf) (center(j),j=1,3)
         write(unitwf) nvctr_c, nvctr_f
        end if

! coarse part
        do iseg=1,nseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
            tt=psi_c(i-i0+jj) 
            if (useFormattedOutput) then
              write(unitwf,'(3(i4),1x,e19.12)') i,i2,i3,tt
            else
              write(unitwf) i,i2,i3,tt
            end if
          enddo
         enddo
                                                                                                                             
! fine part
        do iseg=1,nseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
            t1=psi_f(1,i-i0+jj)
            t2=psi_f(2,i-i0+jj)
            t3=psi_f(3,i-i0+jj)
            t4=psi_f(4,i-i0+jj)
            t5=psi_f(5,i-i0+jj)
            t6=psi_f(6,i-i0+jj)
            t7=psi_f(7,i-i0+jj)
            if (useFormattedOutput) then
              write(unitwf,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
              write(unitwf) i,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
          enddo
         enddo

    write(*,*) iorb,'th wavefunction written'


    END SUBROUTINE
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,rx,ry,rz,psi)
        implicit real*8 (a-h,o-z)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f)
        real*8, allocatable :: psifscf(:),psir(:),psig(:,:,:,:,:,:)

        allocate(psig(0:n1,2,0:n2,2,0:n3,2))
        allocate( psifscf((2*n1+31)*(2*n2+31)*(2*n3+16)) )
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))

        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1),psi(nvctr_c+1),psig)
        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n1+2))
        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

        call plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,psir)

        deallocate(psig,psifscf,psir)
        return
    END SUBROUTINE




        subroutine plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,pot)
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=.5d0*hgrid
        open(iounit) 
        open(iounit+1) 
        open(iounit+2) 

        i3=nint(rz/hgridh)
        i2=nint(ry/hgridh)
        write(*,*) 'plot_p, i2,i3,n2,n3 ',i2,i3,n2,n3
        do i1=-14,2*n1+16
        write(iounit,*) i1*hgridh,pot(i1,i2,i3)
        enddo

        i1=nint(rx/hgridh)
        i2=nint(ry/hgridh)
        write(*,*) 'plot_p, i1,i2 ',i1,i2
        do i3=-14,2*n3+16
        write(iounit+1,*) i3*hgridh,pot(i1,i2,i3)
        enddo

        i1=nint(rx/hgridh)
        i3=nint(rz/hgridh)
        write(*,*) 'plot_p, i1,i3 ',i1,i3
        do i2=-14,2*n2+16
        write(iounit+2,*) i2*hgridh,pot(i1,i2,i3)
        enddo

        close(iounit) 
        close(iounit+1) 
        close(iounit+2) 

        return
        END SUBROUTINE



        subroutine plot_psifscf(iunit,hgrid,n1,n2,n3,psifscf)
        implicit real*8 (a-h,o-z)
        dimension psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

    hgridh=.5d0*hgrid

! along x-axis
    i3=n3
    i2=n2
    do i1=-7,2*n1+8
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
    enddo 

! 111 diagonal
    do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=i ; i2=i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
    enddo 

! 1-1-1 diagonal
    do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=i ; i2=-i ; i3=-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
    enddo 

! -11-1 diagonal
    do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=-i ; i2=i ; i3=-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
    enddo 

! -1-11 diagonal
    do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=-i ; i2=-i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
    enddo 

        return
        END SUBROUTINE

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      subroutine atomkin(l,ng,xp,psiat,psiatn,ek)
! calculates the kinetic energy of an atomic wavefunction expressed in Gaussians
! the output psiatn is a normalized version of psiat
        implicit real*8 (a-h,o-z)
        dimension xp(ng),psiat(ng),psiatn(ng)

!        gml=.5d0*gamma(.5d0+l)
        gml = 0.d0
        if (l.eq.0) then 
            gml=0.88622692545275801365d0
        else if (l.eq.1) then 
            gml=0.44311346272637900682d0
!        else if (l.eq.2) then 
!            gml=0.66467019408956851024d0
!        else if (l.eq.3) then 
!            gml=1.6616754852239212756d0
        else
          stop 'atomkin'
        endif

        ek=0.d0
        tt=0.d0
        do i=1,ng
        xpi=.5d0/xp(i)**2
        do j=1,ng
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
        if (l.eq.0) then  ! multiply with 1/sqrt(4*pi)
        tt=tt*0.28209479177387814347d0
        else if (l.eq.1) then  ! multiply with sqrt(3/(4*pi))
        tt=tt*0.48860251190291992159d0
        endif
        do i=1,ng
        psiatn(i)=psiat(i)*tt
        enddo

        return
        END SUBROUTINE



        subroutine diisstp(parallel,norb,norbp,nproc,iproc,  & 
                   ads,ids,mids,idsx,nvctrp,psit,psidst,hpsidst)
! diis subroutine:
! calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
! using  the previous iteration points phidst and the associated error 
! vectors (preconditione gradients) hpsidst
        implicit real*8 (a-h,o-z)
        include 'mpif.h'
        logical parallel
        dimension psit(nvctrp,norbp*nproc),ads(idsx+1,idsx+1,3), &
        psidst(nvctrp,norbp*nproc,idsx),hpsidst(nvctrp,norbp*nproc,idsx)
        allocatable :: ipiv(:),rds(:)

        call timing(iproc,'Diis          ','ON')

        allocate(ipiv(idsx+1),rds(idsx+1))

! set up DIIS matrix (upper triangle)
        if (ids.gt.idsx) then
! shift left up matrix
        do 3079,i=1,idsx-1
        do 3079,j=1,i
3079    ads(j,i,1)=ads(j+1,i+1,1)
        endif

! calculate new line, use rds as work array for summation
        call razero(idsx,rds)
        ist=max(1,ids-idsx+1)
        do i=ist,ids
           mi=mod(i-1,idsx)+1
           do iorb=1,norb
              tt=DDOT(nvctrp,hpsidst(1,iorb,mids),1,hpsidst(1,iorb,mi),1)
              rds(i-ist+1)=rds(i-ist+1)+tt
           end do
        end do

        if (parallel) then
           call MPI_ALLREDUCE(rds,ads(1,min(idsx,ids),1),min(ids,idsx),  & 
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        else
           do i=1,min(ids,idsx)
              ads(i,min(idsx,ids),1)=rds(i)
           end do
        endif


! copy to work array, right hand side, boundary elements
        do 3983,j=1,min(idsx,ids)
        ads(j,min(idsx,ids)+1,2)=1.d0
        rds(j)=0.d0
        do 3983,i=j,min(idsx,ids)
        ads(j,i,2)=ads(j,i,1)
3983    continue
        ads(min(idsx,ids)+1,min(idsx,ids)+1,2)=0.d0
        rds(min(idsx,ids)+1)=1.d0

!        write(6,*) 'DIIS matrix'
!        do i=1,min(idsx,ids)+1
!        write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(idsx,ids)+1),rds(i)
!        enddo
        if (ids.gt.1) then
! solve linear system:(LAPACK)
        call DSYSV('U',min(idsx,ids)+1,1,ads(1,1,2),idsx+1,  & 
                   ipiv,rds,idsx+1,ads(1,1,3),(idsx+1)**2,info)
        if (info.ne.0) print*, 'DGESV',info
        if (info.ne.0) stop 'DGESV'
        else
        rds(1)=1.d0
        endif
        if (iproc.eq.0) then 
           write(6,*) 'DIIS weights'
           write(6,'(12(1x,e9.2))') (rds(j),j=1,min(idsx,ids)+1)
        endif

! new guess
        do 6633,iorb=1,norb
        call razero(nvctrp,psit(1,iorb))

        jst=max(1,ids-idsx+1)
        jj=0
        do 6612,j=jst,ids
        jj=jj+1
        mj=mod(j-1,idsx)+1
        do 6612,k=1,nvctrp
        psit(k,iorb)=psit(k,iorb)+rds(jj)*(psidst(k,iorb,mj)-hpsidst(k,iorb,mj))
6612    continue
6633    continue

        deallocate(ipiv,rds)
        call timing(iproc,'Diis          ','OF')

        return
        END SUBROUTINE

end module
