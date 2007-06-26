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
  !- Create the ionic potential from the distribution of charges.
  public :: createIonicPotential
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
  public :: reformatonewave
  !- MPI communications for wavefunctions
  public :: switch_waves
  public :: unswitch_waves

  !- SCF handling
  public :: HamiltonianApplication
  public :: hpsitopsi

  !- Geometry method
  !- Get a box that contains all atoms and their active grid points.
  public :: system_size
  
  !- Linear algebra method
  !- Compute a scalar product
!  public :: wnrm
!  public :: wpdot

  !estimation of the memory
  public :: MemoryEstimator

contains

subroutine cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames, rxyz, energy, fxyz,  &
     & psi, keyg, keyv, nvctr_c, nvctr_f, nseg_c, nseg_f, norbp, norb, eval, inputPsiId, &
     & output_grid, output_wf, n1, n2, n3, hgrid, rxyz_old, infocode)
  ! inputPsiId = 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
  ! inputPsiId = 1 : read waves from argument psi, using n1, n2, n3, hgrid and rxyz_old
  !                  as definition of the previous system.
  ! inputPsiId = 2 : read waves from disk
  ! does an electronic structure calculation. Output is the total energy and the forces 
  ! psi, keyg, keyv and eval should be freed after use outside of the routine.
  ! infocode -> encloses some information about the status of the run
  !          =0 run succesfully succeded
  !          =1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
  !             forces may be meaningless   
  !          =2 (present only for inputPsiId=1) gnrm of the first iteration > 1 AND growing in
  !             the second iteration OR grnm 1st >2.
  !             Input wavefunctions need to be recalculated. Routine exits.
  !          =3 (present only for inputPsiId=0) gnrm > 4. SCF error. Routine exits.

  use Poisson_Solver

  implicit real*8 (a-h,o-z)
  character*30 label
  character*27 filename
  character*20 atomnames
  character*2 symbol
  character*1 datacode
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
  real*8, pointer :: hpsi(:,:)

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
  allocatable :: psppar(:,:,:),nelpsp(:),radii_cf(:,:),npspcode(:),nzatom(:),iasctype(:)
  allocatable :: derproj(:)
  ! arrays for DIIS convergence accelerator
  real*8, pointer :: ads(:,:,:),psidst(:,:,:),hpsidst(:,:,:)

  ! arrays for calculation of forces and tail correction to kinetic energy
  allocatable :: rho(:),pot(:,:,:)
  allocatable :: neleconf(:,:),nscatterarr(:,:),ngatherarr(:,:)

  integer :: ierror

  include 'mpif.h'

  if (iproc.eq.0) write(*,'(1x,a,1x,i0)') 'CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER',inputPsiId
  if (parallel) then
     call timing(iproc,'parallel     ','IN')
  else
     call timing(iproc,'             ','IN')
  end if
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
     !add the number of distributed point for the compressed wavefunction
     tt=dble(nvctr_c_old+7*nvctr_f_old)/dble(nproc)
     nvctrp_old=int((1.d0-eps_mach*tt) + tt)

     !allocations
     allocate(keyg_old(2,nseg_c_old+nseg_f_old),stat=i_stat)
     call memocc(i_stat,product(shape(keyg_old))*kind(keyg_old),'keyg_old','cluster')
     allocate(keyv_old(nseg_c_old+nseg_f_old),stat=i_stat)
     call memocc(i_stat,product(shape(keyv_old))*kind(keyv_old),'keyv_old','cluster')
     allocate(psi_old(nvctr_c_old+7*nvctr_f_old,norbp),stat=i_stat)
     call memocc(i_stat,product(shape(psi_old))*kind(psi_old),'psi_old','cluster')
     allocate(eval_old(norb),stat=i_stat)
     call memocc(i_stat,product(shape(eval_old))*kind(eval_old),'eval_old','cluster')

     do iseg=1,nseg_c_old+nseg_f_old
        keyg_old(1,iseg)    = keyg(1,iseg)
        keyg_old(2,iseg)    = keyg(2,iseg)
        keyv_old(iseg)      = keyv(iseg)
     enddo
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        tt=0.d0
        !the control between the different psi
        !must be such that the allocation dimensions are respected
        !(remember that the allocations are slightly enlarged to avoid
        !allocation of work arrays for transposition in the parallel case)
        do j=1,nvctr_c_old+7*nvctr_f_old
           ind=j+(nvctr_c_old+7*nvctr_f_old)*(iorb-iproc*norbp-1)
           i1=mod(ind,nvctrp_old)
           i2=(ind-i1)/nvctrp_old+1
           psi_old(j,iorb-iproc*norbp)     = psi(i1,i2)
           tt=tt+psi(i1,i2)**2
        enddo
        tt=sqrt(tt)
        if (abs(tt-1.d0).gt.1.d-8) then
           write(*,*)'wrong psi_old',iorb,tt
           stop 
        end if
        eval_old(iorb) = eval(iorb)
     enddo
     !deallocation
     i_all=-product(shape(keyg))*kind(keyg)
     deallocate(keyg,stat=i_stat)
     call memocc(i_stat,i_all,'keyg','cluster')
     i_all=-product(shape(keyv))*kind(keyv)
     deallocate(keyv,stat=i_stat)
     call memocc(i_stat,i_all,'keyv','cluster')
     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi','cluster')
     i_all=-product(shape(eval))*kind(eval)
     deallocate(eval,stat=i_stat)
     call memocc(i_stat,i_all,'eval','cluster')

  end if

  !temporary flag, added for debugging purposes
  !in case of doubts put these flags to the "robust" position 'G'
  datacode='D'
  if (.not. parallel) datacode='G'

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
  read(1,*) ixc
  read(1,*) ncharge,elecfield
  read(1,*) gnrm_cv
  read(1,*) itermax
  read(1,*) ncong
  read(1,*) idsx
  read(1,*) calc_tail
  read(1,*) rbuf
  read(1,*) ncongt
  close(1)
 
  if (iproc.eq.0) then 
     write(*,'(1x,a)')&
          '------------------------------------------------------------------- Input Parameters'
     write(*,'(1x,a)')&
          '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
     write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
          'Grid spacing=',hgrid,    '|  Coarse Wfs.=',crmult,'| Wavefns Conv.=',gnrm_cv,&
          '| Calculate=',calc_tail
     write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,f4.1)')&
          '       XC id=',ixc,      '|    Fine Wfs.=',frmult,'| Max. N. Iter.=',itermax,&
          '| Extension=',rbuf
     write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,i4)')&
          'total charge=',ncharge,  '| Coarse Proj.=',cpmult,'| CG Prec.Steps=',ncong,&
          '|  CG Steps=',ncongt
     write(*,'(1x,a,1pe7.1,1x,a,0pf5.2,1x,a,i8)')&
          ' elec. field=',elecfield,'|   Fine Proj.=',fpmult,'| DIIS Hist. N.=',idsx

!!$     write(*,'(1x,a,f6.3)')    'hgrid=',hgrid
!!$     write(*,'(1x,a,f6.3)')    'crmult=',crmult
!!$     write(*,'(1x,a,f6.3)')    'frmult=',frmult
!!$     write(*,'(1x,a,f6.3)')    'cpmult=',cpmult
!!$     write(*,'(1x,a,f6.3)')    'fpmult=',fpmult
!!$     write(*,'(1x,a,i0)')      'ixc= ',ixc
!!$     write(*,'(1x,a,i0)')      'ncharge= ',ncharge
!!$     write(*,'(1x,a,1pe9.2)')  'electric_field=',elecfield
!!$     write(*,'(1x,a,1pe9.2)')  'gnrm_cv=',gnrm_cv
!!$     write(*,'(1x,a,i0)')      'itermax= ',itermax
!!$     write(*,'(1x,a,i0)')      'ncong= ',ncong
!!$     write(*,'(1x,a,i0)')      'idsx= ',idsx
!!$     write(*,*)                'calc_tail',calc_tail
!!$     write(*,'(1x,a,f6.3)')    'rbuf=',rbuf
!!$     write(*,'(1x,a,i0)')      'ncongt= ',ncongt
  endif


  ! grid spacing (same in x,y and z direction)
  hgridh=.5d0*hgrid

  ! store PSP parameters
  ! modified to accept both GTH and HGH pseudopotential types
  !allocation
  allocate(psppar(0:4,0:4,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(psppar))*kind(psppar),'psppar','cluster')
  allocate(nelpsp(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nelpsp))*kind(nelpsp),'nelpsp','cluster')
  allocate(radii_cf(ntypes,2),stat=i_stat)
  call memocc(i_stat,product(shape(radii_cf))*kind(radii_cf),'radii_cf','cluster')
  allocate(npspcode(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(npspcode))*kind(npspcode),'npspcode','cluster')
  allocate(nzatom(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nzatom))*kind(nzatom),'nzatom','cluster')
  allocate(iasctype(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(iasctype))*kind(iasctype),'iasctype','cluster')
  allocate(neleconf(6,0:3),stat=i_stat)
  call memocc(i_stat,product(shape(neleconf))*kind(neleconf),'neleconf','cluster')

  if (iproc==0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------------ System Properties'
     write(*,'(1x,a)')&
          'Atom Name   Ext.Electrons  PSP Code  Radii: Coarse     Fine   Calculated   From File'

  end if

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
     read(11,*) nzatom(ityp),nelpsp(ityp)
     read(11,*) npspcode(ityp)
     psppar(:,:,ityp)=0.d0
     read(11,*) (psppar(0,j,ityp),j=0,4)
     if (npspcode(ityp) == 2) then !GTH case
        do i=1,2
           read(11,*) (psppar(i,j,ityp),j=0,3-i)
        enddo
     else if (npspcode(ityp) == 3) then !HGH case
        read(11,*) (psppar(1,j,ityp),j=0,3)
        do i=2,4
           read(11,*) (psppar(i,j,ityp),j=0,3)
           read(11,*) !k coefficients, not used (no spin-orbit coupling)
        enddo
     else
        if (iproc == 0) then
           write(*,'(1x,a,a)')trim(atomnames(ityp)),&
                'unrecognized pspcode (accepts only GTH & HGH pseudopotentials in ABINIT format)'
        end if
        stop
     end if
     !see whether the atom is semicore or not
     call eleconf(nzatom(ityp),nelpsp(ityp),symbol,rcov,rprb,ehomo,neleconf,iasctype(ityp))
     !if you want no semicore electrons, uncomment the following line
     !iasctype(ityp)=0

     !old way of calculating the radii, requires modification of the PSP files
     read(11,*,iostat=ierror) radii_cf(ityp,1),radii_cf(ityp,2)
     if (ierror.eq.0) then
        if (iproc==0) write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
             trim(atomnames(ityp)),nelpsp(ityp),npspcode(ityp),&
             radii_cf(ityp,1),radii_cf(ityp,2),&
             '                   X    '
        !if (iproc.eq.0) write(*,'(1x,a,a)') &
        !   'radii is given in the psp file for atom type ',trim(atomnames(ityp))
     else
        !new method for assigning the radii
        radii_cf(ityp,1)=1.d0/sqrt(abs(2.d0*ehomo))
        radfine=100.d0
        do i=0,4
           if (psppar(i,0,ityp)/=0.d0) then
              radfine=min(radfine,psppar(i,0,ityp))
           end if
        end do
        radii_cf(ityp,2)=radfine
        if (iproc==0) write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
             trim(atomnames(ityp)),nelpsp(ityp),npspcode(ityp),&
             radii_cf(ityp,1),radii_cf(ityp,2),&
             '       X                '
     end if
     close(11)
     !if (iproc.eq.0) write(*,'(1x,a,a,a,i0,a,i0,a,2(f8.5))') 'atom type ',trim(atomnames(ityp)), & 
     !     ' is described by ',nelpsp(ityp),' electrons, with pspcode= ',npspcode(ityp),&
     !     ' and radii=',radii_cf(ityp,1),radii_cf(ityp,2)
  enddo

  !deallocation
  i_all=-product(shape(neleconf))*kind(neleconf)
  deallocate(neleconf,stat=i_stat)
  call memocc(i_stat,i_all,'neleconf','cluster')


! Number of orbitals and their occupation number
  norb_vir=0!2 !modified for testing purposes

! Number of electrons and number of semicore atoms
  nelec=0
  natsc=0
  do iat=1,nat
     ityp=iatype(iat)
     nelec=nelec+nelpsp(ityp)
     if (iasctype(ityp) /= 0) natsc=natsc+1
  enddo
  nelec=nelec-ncharge
  if (iproc.eq.0) then
     write(*,'(1x,a,i8)') &
          'Total Number of Electrons ',nelec
     if (mod(nelec,2).ne.0) write(*,*) &
          'WARNING: odd number of electrons, no closed shell system'
  end if
  norb=(nelec+1)/2+norb_vir

  allocate(occup(norb),stat=i_stat)
  call memocc(i_stat,product(shape(occup))*kind(occup),'occup','cluster')
  allocate(eval(norb),stat=i_stat)
  call memocc(i_stat,product(shape(eval))*kind(eval),'eval','cluster')

!!$  occup(1)=2.d0 !added for testing purposes
!!$  do iorb=2,norb
!!$     occup(iorb)=2.d0/3.d0
!!$  enddo
!!$  nt=4

  nt=0
  do iorb=1,norb
     it=min(2,nelec-nt)
     occup(iorb)=it
     nt=nt+it
  enddo

  if (iproc.eq.0) then 
     if (norb_vir /=0) write(*,'(1x,a,i8)') &
          '         Virtual Orbitals ',norb_vir
     write(*,'(1x,a,i8)') &
          'Total Number of  Orbitals ',norb
     !write(*,'(1x,a,i0)') 'number of orbitals ',norb
     iorb1=1
     rocc=occup(1)
     do iorb=1,norb
        if (occup(iorb) /= rocc) then
           if (iorb1 == iorb-1) then
              write(*,'(1x,a,i0,a,f3.1)') 'occup(',iorb1,')= ',rocc
           else
              write(*,'(1x,a,i0,a,i0,a,f3.1)') 'occup(',iorb1,':',iorb-1,')= ',rocc
           end if
           rocc=occup(iorb)
           iorb1=iorb
        end if
     enddo
     if (iorb1 == norb) then
        write(*,'(1x,a,i0,a,f3.1)') 'occup(',norb,')= ',occup(norb)
     else
        write(*,'(1x,a,i0,a,i0,a,f3.1)') 'occup(',iorb1,':',norb,')= ',occup(norb)
     end if
  endif

! determine size alat of overall simulation cell
  call system_size(nat,rxyz,radii_cf(1,1),crmult,iatype,ntypes, &
       cxmin,cxmax,cymin,cymax,czmin,czmax)
  alat1=(cxmax-cxmin)
  alat2=(cymax-cymin)
  alat3=(czmax-czmin)

!!$! grid sizes n1,n2,n3 !added for testing purposes
!!$  n1=int(alat1/hgrid)
!!$  if (mod(n1,2).eq.0) n1=n1+1
!!$  n2=int(alat2/hgrid)
!!$  if (mod(n2,2).eq.0) n2=n2+1
!!$  n3=int(alat3/hgrid)
!!$  if (mod(n3,2).eq.0) n3=n3+1
!!$  alat1=n1*hgrid 
!!$  alat2=n2*hgrid 
!!$  alat3=n3*hgrid
!!$  do iat=1,nat
!!$     rxyz(1,iat)=(real(n1/2,kind=8)+0.5)*hgrid 
!!$     rxyz(2,iat)=(real(n1/2,kind=8)+0.5)*hgrid 
!!$     rxyz(3,iat)=(real(n1/2,kind=8)+0.5)*hgrid 
!!$  enddo


! shift atomic positions such that molecule is inside cell
  !if (iproc.eq.0) write(*,'(1x,a,3(1x,1pe14.7))') 'Atomic positions shifted by',-cxmin,-cymin,-czmin

  do iat=1,nat
     rxyz(1,iat)=rxyz(1,iat)-cxmin
     rxyz(2,iat)=rxyz(2,iat)-cymin
     rxyz(3,iat)=rxyz(3,iat)-czmin
  enddo

  if (iproc.eq.0) then
     write(*,'(1x,a,19x,a)') 'Shifted atomic positions, Atomic Units:','grid spacing units:'
     do iat=1,nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5),3x,3(1x,0pf9.3))') &
             iat,trim(atomnames(iatype(iat))),&
             (rxyz(j,iat),j=1,3),rxyz(1,iat)/hgrid,rxyz(2,iat)/hgrid,rxyz(3,iat)/hgrid
     enddo
  endif

! grid sizes n1,n2,n3
  n1=int(alat1/hgrid)
  !if (mod(n1+1,4).eq.0) n1=n1+1
  n2=int(alat2/hgrid)
  !if (mod(n2+1,8).eq.0) n2=n2+1
  n3=int(alat3/hgrid)
  alat1=n1*hgrid ; alat2=n2*hgrid ; alat3=n3*hgrid
  if (iproc.eq.0) then 
     write(*,'(1x,a,3(1x,1pe12.5))') &
          '   Shift of=',-cxmin,-cymin,-czmin
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '  Box Sizes=',alat1,alat2,alat3,n1,n2,n3
     !write(*,'(1x,a,3(1x,i0))') 'n1,n2,n3',n1,n2,n3
     !write(*,'(1x,a,3(1x,i0))') 'total number of grid points',(n1+1)*(n2+1)*(n3+1)
     !write(*,'(1x,a,3(1x,1pe12.5))') 'simulation cell',alat1,alat2,alat3
  endif

! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 ; nfl2=n2 ; nfl3=n3
  nfu1=0 ; nfu2=0 ; nfu3=0
  do iat=1,nat
     rad=radii_cf(iatype(iat),2)*frmult
     nfl1=min(nfl1,int(onem+(rxyz(1,iat)-rad)/hgrid))
     nfu1=max(nfu1,int((rxyz(1,iat)+rad)/hgrid))

     nfl2=min(nfl2,int(onem+(rxyz(2,iat)-rad)/hgrid))
     nfu2=max(nfu2,int((rxyz(2,iat)+rad)/hgrid))

     nfl3=min(nfl3,int(onem+(rxyz(3,iat)-rad)/hgrid)) 
     nfu3=max(nfu3,int((rxyz(3,iat)+rad)/hgrid))
  enddo
  if (iproc.eq.0) then
     write(*,'(1x,a,3x,3(2x,i4,a1,i0))')&
          '      Extremes for the high resolution grid points:',&
          nfl1,'<',nfu1,nfl2,'<',nfu2,nfl3,'<',nfu3
     !write(*,'(1x,a,2(1x,i0))') 'nfl1,nfu1 ',nfl1,nfu1
     !write(*,'(1x,a,2(1x,i0))') 'nfl2,nfu2 ',nfl2,nfu2
     !write(*,'(1x,a,2(1x,i0))') 'nfl3,nfu3 ',nfl3,nfu3
  endif


  !memory estimation
  if (iproc==0) then
     call MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
          rxyz,radii_cf,crmult,frmult,norb,atomnames,.false.)
  end if

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=14
  call createKernel('F',2*n1+31,2*n2+31,2*n3+31,hgridh,hgridh,hgridh,ndegree_ip,&
       iproc,nproc,pkernel)

! Create wavefunctions descriptors and allocate them
  allocate(ibyz_c(2,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(ibyz_c))*kind(ibyz_c),'ibyz_c','cluster')
  allocate(ibxz_c(2,0:n1,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(ibxz_c))*kind(ibxz_c),'ibxz_c','cluster')
  allocate(ibxy_c(2,0:n1,0:n2),stat=i_stat)
  call memocc(i_stat,product(shape(ibxy_c))*kind(ibxy_c),'ibxy_c','cluster')
  allocate(ibyz_f(2,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(ibyz_f))*kind(ibyz_f),'ibyz_f','cluster')
  allocate(ibxz_f(2,0:n1,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(ibxz_f))*kind(ibxz_f),'ibxz_f','cluster')
  allocate(ibxy_f(2,0:n1,0:n2),stat=i_stat)
  call memocc(i_stat,product(shape(ibxy_f))*kind(ibxy_f),'ibxy_f','cluster')
  
  call createWavefunctionsDescriptors(iproc,nproc,idsx,n1,n2,n3,output_grid,hgrid,&
       nat,ntypes,iatype,atomnames,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
       ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp,keyg,keyv,&
       norb,norbp)

! Calculate all projectors
  allocate(nseg_p(0:2*nat),stat=i_stat)
  call memocc(i_stat,product(shape(nseg_p))*kind(nseg_p),'nseg_p','cluster')
  allocate(nvctr_p(0:2*nat),stat=i_stat)
  call memocc(i_stat,product(shape(nvctr_p))*kind(nvctr_p),'nvctr_p','cluster')
  allocate(nboxp_c(2,3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(nboxp_c))*kind(nboxp_c),'nboxp_c','cluster')
  allocate(nboxp_f(2,3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(nboxp_f))*kind(nboxp_f),'nboxp_f','cluster')

  call createProjectorsArrays(iproc, n1, n2, n3, rxyz, nat, ntypes, iatype, atomnames, &
       & psppar, npspcode, radii_cf, cpmult, fpmult, hgrid, nvctr_p, nseg_p, &
       & keyg_p, keyv_p, nproj, nprojel, istart, nboxp_c, nboxp_f, proj)


  !if (iproc.eq.0) write(*,'(1x,a,3(1x,i0))') &
  !     'Size of real space grids',(2*n1+31),(2*n2+31),(2*n3+31)
    
  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(nscatterarr(0:nproc-1,4),stat=i_stat)
  call memocc(i_stat,product(shape(nscatterarr))*kind(nscatterarr),'nscatterarr','cluster')
  if (datacode == 'D') then
     do jproc=0,iproc-1
        call PS_dim4allocation('F',datacode,jproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d            !number of planes for the density
        nscatterarr(jproc,2)=n3p            !number of planes for the potential
        nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
        nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
     end do
     do jproc=iproc+1,nproc-1
        call PS_dim4allocation('F',datacode,jproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d
        nscatterarr(jproc,2)=n3p
        nscatterarr(jproc,3)=i3s+i3xcsh-1
        nscatterarr(jproc,4)=i3xcsh
     end do
  end if

  call PS_dim4allocation('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s)
  nscatterarr(iproc,1)=n3d
  nscatterarr(iproc,2)=n3p
  nscatterarr(iproc,3)=i3s+i3xcsh-1
  nscatterarr(iproc,4)=i3xcsh

  !allocate array for the communications of the potential
  allocate(ngatherarr(0:nproc-1,2),stat=i_stat)
  call memocc(i_stat,product(shape(ngatherarr))*kind(ngatherarr),'ngatherarr','cluster')
  ngatherarr(:,1)=(2*n1+31)*(2*n2+31)*nscatterarr(:,2)
  ngatherarr(:,2)=(2*n1+31)*(2*n2+31)*nscatterarr(:,3)


  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion((2*n1+31)*(2*n2+31)*n3pi),stat=i_stat)
     call memocc(i_stat,product(shape(pot_ion))*kind(pot_ion),'pot_ion','cluster')
  else
     allocate(pot_ion(1),stat=i_stat)
     call memocc(i_stat,product(shape(pot_ion))*kind(pot_ion),'pot_ion','cluster')
  end if

  call createIonicPotential(iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,hgrid,&
     elecfield,n1,n2,n3,n3pi,i3s+i3xcsh,pkernel,pot_ion,eion)

  !Allocate Charge density, Potential in real space
  if (n3d >0) then
     allocate(rhopot((2*n1+31),(2*n2+31),n3d),stat=i_stat)
     call memocc(i_stat,product(shape(rhopot))*kind(rhopot),'rhopot','cluster')
  else
     allocate(rhopot(1,1,1),stat=i_stat)
     call memocc(i_stat,product(shape(rhopot))*kind(rhopot),'rhopot','cluster')
  end if

     ! INPUT WAVEFUNCTIONS
  if (inputPsiId == -1) then
     call input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          nat,natsc,norb,norbp,n1,n2,n3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
          rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
          atomnames,ntypes,iatype,iasctype,pkernel,nzatom,nelpsp,psppar,npspcode,&
          ixc,psi,psit,eval,accurex,datacode,nscatterarr,ngatherarr)
     if (iproc.eq.0) then
        write(*,'(1x,a,1pe9.2)') 'expected accuracy in kinetic energy due to grid size',accurex
        write(*,'(1x,a,1pe9.2)') 'suggested value for gnrm_cv ',accurex/norb
     endif

     if (parallel) then
        !allocate hpsi array (used also as transposed)
        !allocated in the transposed way such as 
        !it can also be used as the transposed hpsi
        allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','cluster')
     else
        !allocate hpsi array
        allocate(hpsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','cluster')
     endif

  else 

     !allocate principal wavefunction
     if (parallel) then
        !allocated in the transposed way such as 
        !it can also be used as a work array for transposition
        allocate(psi(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     else
        allocate(psi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     end if

     if (inputPsiId == 0 ) then !WARNING TO BE CHANGED
        !apply the results of the gaussian wavefunctions
        call gautowav(iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
             nvctr_c,nvctr_f,nseg_c,nseg_f,keyg,keyv,iatype,rxyz,hgrid,psi)

     else if (inputPsiId == 1 ) then

        if (iproc.eq.0) write(*,*) 'START reformatting psi from old psi'
        call reformatmywaves(iproc, norb, norbp, nat, &
             & hgrid_old, nvctr_c_old, nvctr_f_old, n1_old, n2_old, n3_old, rxyz_old, &
             & nseg_c_old, nseg_f_old, keyg_old, keyv_old, psi_old, &
             & hgrid, nvctr_c, nvctr_f, n1, n2, n3, rxyz, &
             & nseg_c, nseg_f, keyg, keyv, psi)
        eval=eval_old
        i_all=-product(shape(keyg_old))*kind(keyg_old)
        deallocate(keyg_old,stat=i_stat)
        call memocc(i_stat,i_all,'keyg_old','cluster')
        i_all=-product(shape(keyv_old))*kind(keyv_old)
        deallocate(keyv_old,stat=i_stat)
        call memocc(i_stat,i_all,'keyv_old','cluster')
        i_all=-product(shape(psi_old))*kind(psi_old)
        deallocate(psi_old,stat=i_stat)
        call memocc(i_stat,i_all,'psi_old','cluster')
        i_all=-product(shape(eval_old))*kind(eval_old)
        deallocate(eval_old,stat=i_stat)
        call memocc(i_stat,i_all,'eval_old','cluster')

        !initialise control value for gnrm in the case of a restart
        gnrm_check=0.d0
       
     else if (inputPsiId == 2) then
        
        call readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,&
             nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
        
     end if

     if (parallel) then
        !allocate hpsi array (used also as transposed)
        !allocated in the transposed way such as 
        !it can also be used as the transposed hpsi
        allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','cluster')

        !transpose the psi wavefunction
        call timing(iproc,'Un-Transall   ','ON')
        !here hpsi is used as a work array
        call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,hpsi)
        !allocate transposed principal wavefunction
        allocate(psit(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psit))*kind(psit),'psit','cluster')
        call MPI_ALLTOALL(hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-Transall   ','OF')
        !end of transposition

        call orthon_p(iproc,nproc,norb,norbp,nvctrp,psit)
        !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

        !retranspose the psit wavefunction into psi
        call timing(iproc,'Un-Transall   ','ON')
        !here hpsi is used as a work array
        call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
        call timing(iproc,'Un-Transall   ','OF')
        !end of retransposition
     else
        call orthon(norb,norbp,nvctrp,psi)
        call checkortho(norb,norbp,nvctrp,psi)
        !allocate hpsi array
        allocate(hpsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','cluster')
     endif

  end if

  !no need of using nzatom array
  i_all=-product(shape(nzatom))*kind(nzatom)
  deallocate(nzatom,stat=i_stat)
  call memocc(i_stat,i_all,'nzatom','cluster')

!!$  !plot the initial wavefunctions in the different orbitals
!!$  do i=2*iproc+1,2*iproc+2
!!$     iounit=27+3*(i-1)
!!$     print *,'iounit',iounit,'-',iounit+2
!!$     call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$          rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc))
!!$  end do

! allocate arrays necessary for DIIS convergence acceleration
  if (idsx.gt.0) then
     allocate(psidst(nvctrp,norbp*nproc,idsx),stat=i_stat)
     call memocc(i_stat,product(shape(psidst))*kind(psidst),'psidst','cluster')
     allocate(hpsidst(nvctrp,norbp*nproc,idsx),stat=i_stat)
     call memocc(i_stat,product(shape(hpsidst))*kind(hpsidst),'hpsidst','cluster')
     allocate(ads(idsx+1,idsx+1,3),stat=i_stat)
     call memocc(i_stat,product(shape(ads))*kind(ads),'ads','cluster')
     call razero(3*(idsx+1)**2,ads)
  endif

  alpha=1.d0
  energy=1.d100
  gnrm=1.d100
  ekin_sum=0.d0 ; epot_sum=0.d0 ; eproj_sum=0.d0
! loop for wavefunction minimization
  do 1000, iter=1,itermax
     if (idsx.gt.0) mids=mod(iter-1,idsx)+1
     if (iproc.eq.0) then 
        write(*,'(1x,a,i0)')&
         '---------------------------------------------------------------------------- iter= ',&
         iter
        if (gnrm.le.gnrm_cv) then
!!$           write(*,'(1x,a,i3,3(1x,1pe18.11))') 'iproc,ehart,eexcu,vexcu',iproc,ehart,eexcu,vexcu
!!$           write(*,'(1x,a,3(1x,1pe18.11))') 'final ekin_sum,epot_sum,eproj_sum',ekin_sum,epot_sum,eproj_sum
!!$           write(*,'(1x,a,3(1x,1pe18.11))') 'final ehart,eexcu,vexcu',ehart,eexcu,vexcu
!!$           write(*,'(1x,a,i6,2x,1pe19.12,1x,1pe9.2)') 'FINAL iter,total energy,gnrm',iter,energy,gnrm
           !write(61,*)hgrid,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu

        endif
     endif

     ! Potential from electronic charge density
     if (datacode=='G') then
        call sumrho_old(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
             nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)
     else
        call sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
             nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot, &
             (2*n1+31)*(2*n2+31)*n3d,nscatterarr)
     end if

!     ixc=12  ! PBE functional
!     ixc=1   ! LDA functional
     call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.)

     call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
     psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,ngatherarr,n3p,&
     rhopot(1,1,1+i3xcsh),psi,hpsi,ekin_sum,epot_sum,eproj_sum)

     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-ehart+eexcu-vexcu+eion

     !check for convergence or whether max. numb. of iterations exceeded
     if (gnrm.le.gnrm_cv .or. iter.eq.itermax) then 
        if (iproc.eq.0) then 
           write(*,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write(*,'(1x,a,i3,3(1x,1pe18.11))') &
                'iproc,ehart,eexcu,vexcu',iproc,ehart,eexcu,vexcu
           write(*,'(1x,a,3(1x,1pe18.11))') &
                'final ekin_sum,epot_sum,eproj_sum',ekin_sum,epot_sum,eproj_sum
           write(*,'(1x,a,3(1x,1pe18.11))') &
                'final ehart,eexcu,vexcu',ehart,eexcu,vexcu
           write(*,'(1x,a,i6,2x,1pe19.12,1x,1pe9.2)') &
                'FINAL iter,total energy,gnrm',iter,energy,gnrm

        end if
        infocode=0
        goto 1010
     endif

     call hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp,nseg_c,nseg_f,&
     keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
     psi,psit,hpsi,psidst,hpsidst)

     if (inputPsiId == 0) then
        if (gnrm > 4.d0) then
           if (iproc == 0) write(*,'(1x,a)')&
                'Error: the norm of the residue is too large also with input wavefunctions.'
           infocode=3
           return
        end if
     else if (inputPsiId == 1) then
        if (gnrm > 2.d0) then
           if (iproc == 0) write(*,'(1x,a)')&
                'The norm of the residue is too large, need to recalculate input wavefunctions'
           infocode=2
           return
        else if (iter == 1 .and. gnrm > 1.d0) then
        !check the value of the first gnrm to see whether it is the case
        !to recalculate input guess
           gnrm_check=gnrm
        else if (iter == 2 .and. gnrm_check > 1.d0) then
           !control whether it is the case to exit the program
            if (gnrm >= gnrm_check) then
               if (iproc == 0) write(*,'(1x,a)')&
                   'The norm of the residue is growing, need to recalculate input wavefunctions'
              infocode=2
              return
           end if
        end if
     end if


     tt=energybs-scprsum
     if (abs(tt).gt.1.d-8) then 
        write(*,'(1x,a,3(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
     endif
     if (iproc.eq.0) then
        write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             ekin_sum,epot_sum,eproj_sum
        write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
        write(*,'(1x,a,i6,2x,1pe19.12,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
     endif

1000 continue
     write(*,'(1x,a)')'No convergence within the allowed number of minimization steps'
     infocode=1
1010 continue
  if (idsx.gt.0) then
       i_all=-product(shape(psidst))*kind(psidst)
       deallocate(psidst,stat=i_stat)
       call memocc(i_stat,i_all,'psidst','cluster')
       i_all=-product(shape(hpsidst))*kind(hpsidst)
       deallocate(hpsidst,stat=i_stat)
       call memocc(i_stat,i_all,'hpsidst','cluster')
       i_all=-product(shape(ads))*kind(ads)
       deallocate(ads,stat=i_stat)
       call memocc(i_stat,i_all,'ads','cluster')
  end if

!------------------------------------------------------------------------
! transform to KS orbitals
  if (parallel) then

     !transpose the hpsi wavefunction
     call timing(iproc,'Un-Transall   ','ON')
     !here psi is used as a work array
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     !here hpsi is the transposed array
     call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-Transall   ','OF')
     !end of transposition

     call KStrans_p(iproc,nproc,norb,norbp,nvctrp,occup,hpsi,psit,evsum,eval)

     !retranspose the psit wavefunction into psi
     call timing(iproc,'Un-Transall   ','ON')
     !here hpsi is used as a work array
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-Transall   ','OF')
     !end of retransposition

     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit','cluster')

  else
     call KStrans(norb,norbp,nvctrp,occup,hpsi,psi,evsum,eval)
  endif
  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi','cluster')
  if (abs(evsum-energybs).gt.1.d-8 .and. iproc==0) write(*,'(1x,a,2(1x,1pe20.13))')&
       'Difference:evsum,energybs',evsum,energybs

!!$  !plot the converged wavefunctions in the different orbitals
!!$  do i=2*iproc+1,2*iproc+2
!!$     iounit=39+3*(i-1)
!!$     call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$          rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc))
!!$  end do


!  write all the wavefunctions into files
  if (output_wf) then
     call  writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
              nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
     write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
  end if


!------------------------------------------------------------------------
! here we start the calculation of the forces
  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------------- Forces Calculation'
  end if

! Selfconsistent potential is saved in rhopot, 
! new arrays rho,pot for calculation of forces ground state electronic density

  ! Potential from electronic charge density
  if (datacode=='G') then
     allocate(rho((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
     call memocc(i_stat,product(shape(rho))*kind(rho),'rho','cluster')
     call sumrho_old(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho)
  else

     !manipulate scatter array for avoiding the GGA shift
     do jproc=0,nproc-1
        !n3d=n3p
        nscatterarr(jproc,1)=nscatterarr(jproc,2)
        !i3xcsh=0
        nscatterarr(jproc,4)=0
     end do

     if (n3p>0) then
        allocate(rho((2*n1+31)*(2*n2+31)*n3p),stat=i_stat)
        call memocc(i_stat,product(shape(rho))*kind(rho),'rho','cluster')
     else
        allocate(rho(1),stat=i_stat)
        call memocc(i_stat,product(shape(rho))*kind(rho),'rho','cluster')
     end if

     !use pot_ion array for building total rho
     call sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho,&
          (2*n1+31)*(2*n2+31)*n3p,nscatterarr)
     !end of switch

  end if


  if (iproc.eq.0 .and. output_grid) then
     open(unit=22,file='density.pot',status='unknown')
     write(22,*)'density'
     write(22,*) 2*n1,2*n2,2*n3
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     write(22,*)'xyz'
     do i3=1,2*n3
        do i2=1,2*n2
           do i1=1,2*n1
              ind=i1+14+(i2+13)*(2*n1+31)+(i3+13)*(2*n1+31)*(2*n2+31)
              write(22,*)rho(ind)
           end do
        end do
     end do
     close(22)
  endif


  if (datacode == 'D') then

     !switch between the old and the new forces calculation
     i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion','cluster')

     if (n3p>0) then
        allocate(pot((2*n1+31),(2*n2+31),n3p),stat=i_stat)
        call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')
     else
        allocate(pot(1,1,1),stat=i_stat)
        call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')
     end if
     call DCOPY((2*n1+31)*(2*n2+31)*n3p,rho,1,pot,1) 
     call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,0,hgridh,hgridh,hgridh,&
          pot,pkernel,pot,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.)

  else
     allocate(pot((2*n1+31),(2*n2+31),(2*n3+31)),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')
     call DCOPY((2*n1+31)*(2*n2+31)*(2*n3+31),rho,1,pot,1) 
     
     call PSolver('F','G',iproc,nproc,2*n1+31,2*n2+31,2*n3+31,0,hgridh,hgridh,hgridh,&
          pot,pkernel,pot_ion,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.)
     i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion','cluster')
  end if


  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'pkernel','cluster')

  allocate(gxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(gxyz))*kind(gxyz),'gxyz','cluster')

  call timing(iproc,'Forces        ','ON')
  ! calculate local part of the forces gxyz
  call local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
                     n1,n2,n3,n3p,i3s+i3xcsh,rho,pot,gxyz)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho','cluster')
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot','cluster')

  allocate(derproj(3*nprojel),stat=i_stat)
  call memocc(i_stat,product(shape(derproj))*kind(derproj),'derproj','cluster')

  if (iproc == 0) write(*,'(1x,a)',advance='no')'Calculate projectors derivatives...'

  !the calculation of the derivatives of the projectors has been decoupled
  !from the one of nonlocal forces, in this way forces can be calculated
  !diring the minimization if needed
  call projectors_derivatives(iproc,n1,n2,n3,nboxp_c,nboxp_f, & 
     ntypes,nat,norb,nprojel,nproj,&
     iatype,psppar,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,  &
     keyg,keyv,keyg_p,keyv_p,rxyz,radii_cf,cpmult,fpmult,hgrid,derproj)

  if (iproc == 0) write(*,'(1x,a)',advance='no')'done, calculate nonlocal forces...'

  call nonlocal_forces(iproc,ntypes,nat,norb,norbp,nprojel,nproj,&
       iatype,psppar,npspcode,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,derproj,  &
       keyg,keyv,keyg_p,keyv_p,psi,gxyz)

  if (iproc == 0) write(*,'(1x,a)')'done.'
  
  i_all=-product(shape(derproj))*kind(derproj)
  deallocate(derproj,stat=i_stat)
  call memocc(i_stat,i_all,'derproj','cluster')

  i_all=-product(shape(nboxp_c))*kind(nboxp_c)
  deallocate(nboxp_c,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_c','cluster')
  i_all=-product(shape(nboxp_f))*kind(nboxp_f)
  deallocate(nboxp_f,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_f','cluster')

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

  i_all=-product(shape(gxyz))*kind(gxyz)
  deallocate(gxyz,stat=i_stat)
  call memocc(i_stat,i_all,'gxyz','cluster')

  call timing(iproc,'Forces        ','OF')

  !------------------------------------------------------------------------
  if (calc_tail) then
     call timing(iproc,'Tail          ','ON')
!    Calculate energy correction due to finite size effects
     !    ---reformat potential
     allocate(pot((2*n1+31),(2*n2+31),(2*n3+31)),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')
 
     if (datacode=='D') then
        call MPI_ALLGATHERV(rhopot(1,1,1+i3xcsh),(2*n1+31)*(2*n2+31)*n3p,MPI_DOUBLE_PRECISION, &
             pot,ngatherarr(0,1),ngatherarr(0,2), & 
             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        i_all=-product(shape(nscatterarr))*kind(nscatterarr)
        deallocate(nscatterarr,stat=i_stat)
        call memocc(i_stat,i_all,'nscatterarr','cluster')
        i_all=-product(shape(ngatherarr))*kind(ngatherarr)
        deallocate(ngatherarr,stat=i_stat)
        call memocc(i_stat,i_all,'ngatherarr','cluster')
     else
        do i3=1,2*n3+31
           do i2=1,2*n2+31
              do i1=1,2*n1+31
                 pot(i1,i2,i3)=rhopot(i1,i2,i3)
              enddo
           enddo
        enddo
     end if
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot','cluster')
     call timing(iproc,'Tail          ','OF')

     call CalculateTailCorrection(iproc,nproc,n1,n2,n3,rbuf,norb,norbp,nat,ntypes,&
     nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nproj,nprojel,ncongt,&
     keyv,keyg,nseg_p,keyv_p,keyg_p,nvctr_p,psppar,npspcode,eval,&
     pot,hgrid,rxyz,radii_cf,crmult,frmult,iatype,atomnames,&
     proj,psi,occup,output_grid,parallel,ekin_sum,epot_sum,eproj_sum)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot','cluster')

     if (iproc==0) then
        open(61)
        write(61,'(4(f9.3),1x,7(1pe19.11))',advance='no')&
             hgrid,alat1,alat2,alat3,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
     end if
    
     energybs=ekin_sum+epot_sum+eproj_sum
     energy=energybs-ehart+eexcu-vexcu+eion

     if (iproc==0) then
        write(61,'(1pe19.11)')energy
        close(61)
     end if

     if (iproc.eq.0) then
        write(*,'(1x,a,3(1x,1pe18.11))')&
             '  Corrected ekin,epot,eproj',ekin_sum,epot_sum,eproj_sum
        write(*,'(1x,a,1x,1pe19.12)')&
             'Total energy with tail correction',energy
     endif


  else
!    No tail calculation
     if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot','cluster')
     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
     deallocate(nscatterarr,stat=i_stat)
     call memocc(i_stat,i_all,'nscatterarr','cluster')
     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr','cluster')
  endif
! --- End if of tail calculation

  i_all=-product(shape(ibyz_c))*kind(ibyz_c)
  deallocate(ibyz_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibyz_c','cluster')
  i_all=-product(shape(ibxz_c))*kind(ibxz_c)
  deallocate(ibxz_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibxz_c','cluster')
  i_all=-product(shape(ibxy_c))*kind(ibxy_c)
  deallocate(ibxy_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibxy_c','cluster')
  i_all=-product(shape(ibyz_f))*kind(ibyz_f)
  deallocate(ibyz_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibyz_f','cluster')
  i_all=-product(shape(ibxz_f))*kind(ibxz_f)
  deallocate(ibxz_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibxz_f','cluster')
  i_all=-product(shape(ibxy_f))*kind(ibxy_f)
  deallocate(ibxy_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibxy_f','cluster')

  i_all=-product(shape(keyg_p))*kind(keyg_p)
  deallocate(keyg_p,stat=i_stat)
  call memocc(i_stat,i_all,'keyg_p','cluster')
  i_all=-product(shape(keyv_p))*kind(keyv_p)
  deallocate(keyv_p,stat=i_stat)
  call memocc(i_stat,i_all,'keyv_p','cluster')
  i_all=-product(shape(proj))*kind(proj)
  deallocate(proj,stat=i_stat)
  call memocc(i_stat,i_all,'proj','cluster')
  i_all=-product(shape(occup))*kind(occup)
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup','cluster')
  i_all=-product(shape(nvctr_p))*kind(nvctr_p)
  deallocate(nvctr_p,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_p','cluster')
  i_all=-product(shape(nseg_p))*kind(nseg_p)
  deallocate(nseg_p,stat=i_stat)
  call memocc(i_stat,i_all,'nseg_p','cluster')
  i_all=-product(shape(psppar))*kind(psppar)
  deallocate(psppar,stat=i_stat)
  call memocc(i_stat,i_all,'psppar','cluster')
  i_all=-product(shape(nelpsp))*kind(nelpsp)
  deallocate(nelpsp,stat=i_stat)
  call memocc(i_stat,i_all,'nelpsp','cluster')
  i_all=-product(shape(iasctype))*kind(iasctype)
  deallocate(iasctype,stat=i_stat)
  call memocc(i_stat,i_all,'iasctype','cluster')
  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf','cluster')
  i_all=-product(shape(npspcode))*kind(npspcode)
  deallocate(npspcode,stat=i_stat)
  call memocc(i_stat,i_all,'npspcode','cluster')

  call timing(iproc,'              ','RE')
  call cpu_time(tcpu1)
  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  write(*,'(a,1x,i4,2(1x,f12.2))') '- iproc, elapsed, CPU time ', iproc,tel,tcpu1-tcpu0

END SUBROUTINE cluster


subroutine HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
     psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,ngatherarr,n3p,&
     potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum)
  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,nproj,nprojel,n3p
  integer, intent(in) :: nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f
  real(kind=8), intent(in) :: hgrid
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
  integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
  integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
  integer, dimension(ntypes), intent(in) :: npspcode
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(*), intent(in) :: potential
  real(kind=8), dimension(nprojel), intent(in) :: proj
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: hpsi
  real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
  !local variables
  integer :: i_all,i_stat,ierr
  real(kind=8), dimension(:), allocatable :: pot
  real(kind=8), dimension(:,:), allocatable :: wrkallred

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  if (datacode=='D') then
     !allocate full potential
     allocate(pot((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','hamiltonianapplication')

     call timing(iproc,'ApplyLocPotKin','ON')

     call MPI_ALLGATHERV(potential,(2*n1+31)*(2*n2+31)*n3p,&
          MPI_DOUBLE_PRECISION,pot,ngatherarr(0,1),&
          ngatherarr(0,2),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     call timing(iproc,'ApplyLocPotKin','OF')

     call applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,0, &
          hgrid,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,&
          ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          psi,pot,hpsi,epot_sum,ekin_sum)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot','hamiltonianapplication')

  else

     call applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,0, &
          hgrid,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,&
          ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          psi,potential,hpsi,epot_sum,ekin_sum)

  end if

  ! apply all PSP projectors for all orbitals belonging to iproc
  call applyprojectorsall(iproc,ntypes,nat,iatype,psppar,npspcode,occup, &
       nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
       norb,norbp,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)

  if (parallel) then
     allocate(wrkallred(3,2),stat=i_stat)
     call memocc(i_stat,product(shape(wrkallred))*kind(wrkallred),'wrkallred','hamiltonianapplication')

     wrkallred(1,2)=ekin_sum 
     wrkallred(2,2)=epot_sum 
     wrkallred(3,2)=eproj_sum
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     ekin_sum=wrkallred(1,1)
     epot_sum=wrkallred(2,1)
     eproj_sum=wrkallred(3,1) 

     i_all=-product(shape(wrkallred))*kind(wrkallred)
     deallocate(wrkallred,stat=i_stat)
     call memocc(i_stat,i_all,'wrkallred','hamiltonianapplication')
  endif

end subroutine HamiltonianApplication


subroutine hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp,nseg_c,nseg_f,&
     keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
     psi,psit,hpsi,psidst,hpsidst)
  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
  integer, intent(in) :: nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp
  real(kind=8), intent(in) :: hgrid,energy,energy_old
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  real(kind=8), dimension(norb), intent(in) :: occup,eval
  real(kind=8), intent(inout) :: alpha
  real(kind=8), intent(inout) :: gnrm,scprsum
  real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
  real(kind=8), dimension(:,:,:), pointer :: psidst,hpsidst,ads
  !local variables
  integer :: ierr,ind,i1,i2,iorb,k
  real(kind=8) :: tt,scpr,dnrm2

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, orthoconstraint...'
  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  if (parallel) then
     !transpose the hpsi wavefunction
     call timing(iproc,'Un-Transall   ','ON')
     !here psi is used as a work array
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     !here hpsi is the transposed array
     call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-Transall   ','OF')
     !end of transposition

     call  orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsi,scprsum)

     !retranspose the hpsi wavefunction
     call timing(iproc,'Un-Transall   ','ON')
     !here psi is used as a work array
     call MPI_ALLTOALL(hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     !here hpsi is the direct array
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,hpsi)
     call timing(iproc,'Un-Transall   ','OF')
     !end of retransposition
  else
     call orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
  endif

  ! norm of gradient
  gnrm=0.d0
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     !calculate the address to start from for calculating the 
     !norm of the residue if hpsi is allocated in the transposed way
     !it can be eliminated when including all this procedure in a subroutine
     if (parallel) then
        ind=1+(nvctr_c+7*nvctr_f)*(iorb-iproc*norbp-1)
        i1=mod(ind,nvctrp)
        i2=(ind-i1)/nvctrp+1
     else
        i1=1
        i2=iorb-iproc*norbp
     end if
     scpr=dnrm2(nvctr_c+7*nvctr_f,hpsi(i1,i2),1)
     !scpr=dnrm2(nvctr_c+7*nvctr_f,hpsi(1,iorb-iproc*norbp),1)
     !lines for writing the residue following the orbitals
     !if (iorb <=5) write(83,*)iter,iorb,scpr
     gnrm=gnrm+scpr**2
  enddo
  if (parallel) then
     tt=gnrm
     call MPI_ALLREDUCE(tt,gnrm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  gnrm=sqrt(gnrm/norb)

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, preconditioning...'
  end if

  ! Preconditions all orbitals belonging to iproc
  call preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid, &
       ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,&
       ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)

  !       call plot_wf(10,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f, &
  !                       rxyz(1,1),rxyz(2,1),rxyz(3,1),psi)
  !       call plot_wf(20,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f, &
  !                       rxyz(1,1),rxyz(2,1),rxyz(3,1),hpsi)


  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if

  !apply the minimization method (DIIS or steepest descent)
  if (idsx.gt.0) then
     if (parallel) then
        !transpose the hpsi wavefunction into the diis array
        call timing(iproc,'Un-Transall   ','ON')
        !here psi is used as a work array
        call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
        call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsidst(:,:,mids),nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-Transall   ','OF')
        !end of transposition

        do iorb=1,norb
           do k=1,nvctrp
              psidst(k,iorb,mids)= psit(k,iorb) 
           enddo
        enddo

        call diisstp(parallel,norb,norbp,nproc,iproc,  &
             ads,iter,mids,idsx,nvctrp,psit,psidst,hpsidst)
     else
        do iorb=1,norb
           do k=1,nvctrp
              psidst(k,iorb,mids)= psi(k,iorb)
              hpsidst(k,iorb,mids)=hpsi(k,iorb)
           enddo
        enddo

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
     if (iproc.eq.0) write(*,'(1x,a,1pe11.3)') 'alpha=',alpha

     if (parallel) then
        !transpose the hpsi wavefunction
        call timing(iproc,'Un-Transall   ','ON')
        !here psi is used as a work array
        call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
        !here hpsi is the transposed array
        call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-Transall   ','OF')
        !end of transposition

        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psit(1,iorb),1)
        enddo
     else
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psi(1,iorb),1)
        enddo
     endif


  endif

 if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if


  if (parallel) then
     call orthon_p(iproc,nproc,norb,norbp,nvctrp,psit)
     !       call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

     !retranspose the psit wavefunction into psi
     call timing(iproc,'Un-Transall   ','ON')
     !here hpsi is used as a work array
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-Transall   ','OF')
     !end of retransposition
  else
     call orthon(norb,norbp,nvctrp,psi)
     !       call checkortho(norb,norbp,nvctrp,psi)
  endif

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if

end subroutine hpsitopsi

subroutine CalculateTailCorrection(iproc,nproc,n1,n2,n3,rbuf,norb,norbp,nat,ntypes,&
     nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nproj,nprojel,ncongt,&
     keyv,keyg,nseg_p,keyv_p,keyg_p,nvctr_p,psppar,npspcode,eval,&
     pot,hgrid,rxyz,radii_cf,crmult,frmult,iatype,atomnames,&
     proj,psi,occup,output_grid,parallel,ekin_sum,epot_sum,eproj_sum)
implicit none
include 'mpif.h'
logical, intent(in) :: output_grid,parallel
character(len=20), dimension(100), intent(in) :: atomnames
integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,ncongt
integer, intent(in) :: nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nproj,nprojel
real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
integer, dimension(2,nseg_p(2*nat)), intent(inout) :: keyg_p
integer, dimension(ntypes), intent(in) :: npspcode
integer, dimension(nat), intent(in) :: iatype
real(kind=8), dimension(norb), intent(in) :: occup,eval
real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
real(kind=8), dimension(3,nat), intent(in) :: rxyz
real(kind=8), dimension(2*n1+31,2*n2+31,2*n3+31), intent(in) :: pot
real(kind=8), dimension(nprojel), intent(in) :: proj
real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
!local variables
logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f
integer :: iseg,i0,j0,i1,j1,i2,i3,ii,iat,iorb,npt,ipt,i,ierr,i_all,i_stat
integer :: nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,nbuf,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f
integer, dimension(:,:,:), allocatable :: ibbyz_c,ibbyz_f,ibbxz_c,ibbxz_f,ibbxy_c,ibbxy_f
integer, dimension(:), allocatable :: keybv
integer, dimension(:,:), allocatable :: keybg
real(kind=8) :: alatb1,alatb2,alatb3,ekin,epot,eproj,tt,cprecr,sum_tail,ekin1,epot1,eproj1
real(kind=8), dimension(:,:), allocatable :: txyz,wrkallred
real(kind=8), dimension(:), allocatable :: psib,hpsib,psig,psigp,psifscf,psir

call timing(iproc,'Tail          ','ON')

nbuf=nint(rbuf/hgrid)
!    --- new grid sizes n1,n2,n3
nb1=n1+2*nbuf
nb2=n2+2*nbuf
nb3=n3+2*nbuf
alatb1=nb1*hgrid 
alatb2=nb2*hgrid 
alatb3=nb3*hgrid

if (iproc.eq.0) then
   write(*,'(1x,a)')&
        '---------------------------------------------- Estimation of Finite-Size Corrections'
   write(*,'(1x,a,f6.3,a)') &
        'F-S Correction for an effective space of ',rbuf,' AU more around each external atom'
   write(*,'(1x,a,i6,a)') &
        '                  requires the adding of ',nbuf,' additional grid points around cell'
   write(*,'(1x,a,19x,a)') &
        '   Effective box size,   Atomic Units:','grid spacing units:'
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '            ',alatb1,alatb2,alatb3,nb1,nb2,nb3
   !write(*,'(1x,a,i0,a)') &
   !     'BIG: tail requires ',nbuf,' additional grid points around cell'
   !write(*,'(1x,a,3(1x,i0))')     &
   !     'BIG: n1,n2,n3',nb1,nb2,nb3
   !write(*,'(1x,a,1x,i0)')        &
   !     'BIG: total number of grid points',(nb1+1)*(nb2+1)*(nb3+1)
   !write(*,'(1x,a,3(1x,1pe12.5))')&
   !     'BIG: simulation cell',alatb1,alatb2,alatb3
end if


!---reformat keyg_p
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
end do

!---reformat wavefunctions

! fine grid size (needed for creation of input wavefunction, preconditioning)
nbfl1=nfl1+nbuf ; nbfl2=nfl2+nbuf ; nbfl3=nfl3+nbuf
nbfu1=nfu1+nbuf ; nbfu2=nfu2+nbuf ; nbfu3=nfu3+nbuf
if (iproc.eq.0) then
     write(*,'(1x,a,3x,3(2x,i4,a1,i0))')&
          '  Extremes for the new high resolution grid points:',&
          nbfl1,'<',nbfu1,nbfl2,'<',nbfu2,nbfl3,'<',nbfu3
   !write(*,'(1x,a,2(1x,i0))') 'BIG: nfl1,nfu1',nbfl1,nbfu1
   !write(*,'(1x,a,2(1x,i0))') 'BIG: nfl2,nfu2',nbfl2,nbfu2
   !write(*,'(1x,a,2(1x,i0))') 'BIG: nfl3,nfu3',nbfl3,nbfu3
endif

! change atom coordinates according to the enlarged box
allocate(txyz(3,nat),stat=i_stat)
call memocc(i_stat,product(shape(txyz))*kind(txyz),'txyz','calculatetailcorrection')
do iat=1,nat
   txyz(1,iat)=rxyz(1,iat)+nbuf*hgrid
   txyz(2,iat)=rxyz(2,iat)+nbuf*hgrid
   txyz(3,iat)=rxyz(3,iat)+nbuf*hgrid
enddo

! determine localization region for all orbitals, but do not yet fill the descriptor arrays
allocate(logrid_c(0:nb1,0:nb2,0:nb3),stat=i_stat)
call memocc(i_stat,product(shape(logrid_c))*kind(logrid_c),'logrid_c','calculatetailcorrection')
allocate(logrid_f(0:nb1,0:nb2,0:nb3),stat=i_stat)
call memocc(i_stat,product(shape(logrid_f))*kind(logrid_f),'logrid_f','calculatetailcorrection')
allocate(ibbyz_c(2,0:nb2,0:nb3),stat=i_stat)
call memocc(i_stat,product(shape(ibbyz_c))*kind(ibbyz_c),'ibbyz_c','calculatetailcorrection')
allocate(ibbxz_c(2,0:nb1,0:nb3),stat=i_stat)
call memocc(i_stat,product(shape(ibbxz_c))*kind(ibbxz_c),'ibbxz_c','calculatetailcorrection')
allocate(ibbxy_c(2,0:nb1,0:nb2),stat=i_stat)
call memocc(i_stat,product(shape(ibbxy_c))*kind(ibbxy_c),'ibbxy_c','calculatetailcorrection')
allocate(ibbyz_f(2,0:nb2,0:nb3),stat=i_stat)
call memocc(i_stat,product(shape(ibbyz_f))*kind(ibbyz_f),'ibbyz_f','calculatetailcorrection')
allocate(ibbxz_f(2,0:nb1,0:nb3),stat=i_stat)
call memocc(i_stat,product(shape(ibbxz_f))*kind(ibbxz_f),'ibbxz_f','calculatetailcorrection')
allocate(ibbxy_f(2,0:nb1,0:nb2),stat=i_stat)
call memocc(i_stat,product(shape(ibbxy_f))*kind(ibbxy_f),'ibbxy_f','calculatetailcorrection')

! coarse grid quantities
call fill_logrid(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,nbuf,nat,ntypes,iatype,txyz, & 
     radii_cf(1,1),crmult,hgrid,logrid_c)
if (iproc.eq.0 .and. output_grid) then
   write(*,'(1x,a)')&
        'Writing the file describing the new atomic positions of the effective system'
   open(unit=22,file='grid_tail.ascii',status='unknown')
   write(22,*) nat
   write(22,*) alatb1,' 0. ',alatb2
   write(22,*) ' 0. ',' 0. ',alatb3
   do iat=1,nat
      write(22,'(3(1x,e12.5),3x,a20)') txyz(1,iat),txyz(2,iat),txyz(3,iat),atomnames(iatype(iat))
   end do
   do i3=0,nb3
      do i2=0,nb2
         do i1=0,nb1
            if (logrid_c(i1,i2,i3)) then
               write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
            end if
         enddo
      enddo
   end do
endif
call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,nvctrb_c)

if (iproc.eq.0) then
   write(*,'(2(1x,a,i10))') &
        'Coarse resolution grid: Number of segments= ',nsegb_c,'points=',nvctrb_c
   !write(*,'(1x,a,2(1x,i10))') 'BIG: orbitals have coarse segment, elements',nsegb_c,nvctrb_c
end if
call bounds(nb1,nb2,nb3,logrid_c,ibbyz_c,ibbxz_c,ibbxy_c)

! fine grid quantities
call fill_logrid(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,0,nat,ntypes,iatype,txyz, & 
     radii_cf(1,2),frmult,hgrid,logrid_f)
if (iproc.eq.0 .and. output_grid) then
   do i3=0,nb3 
      do i2=0,nb2 
         do i1=0,nb1
            if (logrid_f(i1,i2,i3)) then
               write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
            end if
         enddo
      enddo
   enddo
   close(22)
endif
call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,nvctrb_f)
if (iproc.eq.0) then
   write(*,'(2(1x,a,i10))') &
        '  Fine resolution grid: Number of segments= ',nsegb_f,'points=',nvctrb_f
   !write(*,'(1x,a,2(1x,i10))') 'BIG: orbitals have fine   segment, elements',nsegb_f,7*nvctrb_f
end if
call bounds(nb1,nb2,nb3,logrid_f,ibbyz_f,ibbxz_f,ibbxy_f)

! now fill the wavefunction descriptor arrays
allocate(keybg(2,nsegb_c+nsegb_f),stat=i_stat)
call memocc(i_stat,product(shape(keybg))*kind(keybg),'keybg','calculatetailcorrection')
allocate(keybv(nsegb_c+nsegb_f),stat=i_stat)
call memocc(i_stat,product(shape(keybv))*kind(keybv),'keybv','calculatetailcorrection')
! coarse grid quantities
call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,keybg(1,1),keybv(1))

! fine grid quantities
call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,keybg(1,nsegb_c+1),keybv(nsegb_c+1))

i_all=-product(shape(logrid_c))*kind(logrid_c)
deallocate(logrid_c,stat=i_stat)
call memocc(i_stat,i_all,'logrid_c','calculatetailcorrection')
i_all=-product(shape(logrid_f))*kind(logrid_f)
deallocate(logrid_f,stat=i_stat)
call memocc(i_stat,i_all,'logrid_f','calculatetailcorrection')

! allocations for arrays holding the wavefunction
!if (iproc.eq.0) &
!     write(*,'(1x,a,i0)') 'Allocate words for psib and hpsib ',2*(nvctrb_c+7*nvctrb_f)
allocate(psib(nvctrb_c+7*nvctrb_f),stat=i_stat)
call memocc(i_stat,product(shape(psib))*kind(psib),'psib','calculatetailcorrection')
allocate(hpsib(nvctrb_c+7*nvctrb_f),stat=i_stat)
call memocc(i_stat,product(shape(hpsib))*kind(hpsib),'hpsib','calculatetailcorrection')
!if (iproc.eq.0) write(*,*) 'Allocation done'

! work arrays applylocpotkin
allocate(psig(8*(nb1+1)*(nb2+1)*(nb3+1)),stat=i_stat)
call memocc(i_stat,product(shape(psig))*kind(psig),'psig','calculatetailcorrection')
allocate(psigp(8*(nb1+1)*(nb2+1)*(nb3+1)),stat=i_stat)
call memocc(i_stat,product(shape(psigp))*kind(psigp),'psigp','calculatetailcorrection')
allocate(psifscf(max((2*nb1+31)*(2*nb2+31)*(2*nb3+16),(2*nb1+16)*(2*nb2+31)*(2*nb3+31))),stat=i_stat)
call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','calculatetailcorrection')
allocate(psir((2*nb1+31)*(2*nb2+31)*(2*nb3+31)),stat=i_stat)
call memocc(i_stat,product(shape(psir))*kind(psir),'psir','calculatetailcorrection')

if (iproc.eq.0) then
   write(*,'(1x,a,i0)') &
        'Wavefunction memory occupation in the extended grid (Bytes): ',&
        (nvctrb_c+7*nvctrb_f)*8
   write(*,'(1x,a,i0,a)') &
        'Calculating tail corrections on ',norbp,' orbitals per processor.'
   write(*,'(1x,a)',advance='no') &
        '     orbitals are processed separately'
end if

ekin_sum=0.d0
epot_sum=0.d0
eproj_sum=0.d0
do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

   !build the compressed wavefunction in the enlarged box
   call uncompress(n1,n2,n3,-nbuf,n1+nbuf,-nbuf,n2+nbuf,-nbuf,n3+nbuf, & 
        nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
        nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
        psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psig)
   call compress(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,  &
        nsegb_c,nvctrb_c,keybg(1,1),keybv(1),   &
        nsegb_f,nvctrb_f,keybg(1,nsegb_c+1),keybv(nsegb_c+1),   &
        psig,psib(1),psib(nvctrb_c+1))

   npt=2
   tail_adding: do ipt=1,npt
      
      !calculate gradient
      call applylocpotkinone(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,nbuf, & 
           hgrid,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f,keybg,keybv,  & 
           ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f, & 
           psig,psigp,psifscf,psir,  &
           psib,pot,hpsib,epot,ekin)
      call applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
           nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
           nsegb_c,nsegb_f,keybg,keybv,nvctrb_c,nvctrb_f,  & 
           psib,hpsib,eproj)
      !calculate residue for the single orbital
      tt=0.d0
      do i=1,nvctrb_c+7*nvctrb_f
         hpsib(i)=hpsib(i)-eval(iorb)*psib(i)
         tt=tt+hpsib(i)**2
      enddo
      tt=sqrt(tt)

      if (ipt.eq.npt) exit tail_adding

      !calculate tail using the preconditioner as solver for the green function application
      cprecr=-eval(iorb)
      call timing(iproc,'Tail          ','OF')
      call precong(iorb,nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3, &
           nsegb_c,nvctrb_c,nsegb_f,nvctrb_f,keybg,keybv, &
           ncongt,cprecr,hgrid,ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,hpsib)
      call timing(iproc,'Tail          ','ON')
      !call plot_wf(10,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keybg,keybv,nsegb_f,nvctrb_f,  & 
      !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

      ! add tail to the bulk wavefunction
      sum_tail=0.d0
      do i=1,nvctrb_c+7*nvctrb_f
         psib(i)=psib(i)-hpsib(i)
         sum_tail=sum_tail+psib(i)**2
      enddo
      sum_tail=sqrt(sum_tail)
      !write(*,'(1x,a,i3,3(1x,1pe13.6),1x,1pe9.2)') &
      !     'BIG: iorb,ekin,epot,eproj,gnrm',iorb,ekin,epot,eproj,tt
      !values of the energies before tail application
      ekin1=ekin
      epot1=epot
      eproj1=eproj
      !write(*,'(1x,a,1x,i0,f18.14)') 'norm orbital + tail',iorb,sum_tail
      !call plot_wf(20,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keybg,keybv,nsegb_f,nvctrb_f,  & 
      !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

      sum_tail=1.d0/sum_tail
      do i=1,nvctrb_c+7*nvctrb_f
         psib(i)=psib(i)*sum_tail
      enddo

   end do tail_adding

   !write(*,'(1x,a,i3,3(1x,1pe13.6),2(1x,1pe9.2))') &
   !     'BIG: iorb,denergies,gnrm,dnorm',&
   !     iorb,ekin-ekin1,epot-epot1,eproj-eproj1,tt,sum_tail-1.d0

   if (iproc == 0) then
         write(*,'(a)',advance='no') &
           repeat('.',(iorb*40)/norbp-((iorb-1)*40)/norbp)
   end if
   ekin_sum=ekin_sum+ekin*occup(iorb)
   epot_sum=epot_sum+epot*occup(iorb)
   eproj_sum=eproj_sum+eproj*occup(iorb)
end do

if (iproc == 0) then
   write(*,'(1x,a)')'done.'
end if


i_all=-product(shape(txyz))*kind(txyz)
deallocate(txyz,stat=i_stat)
call memocc(i_stat,i_all,'txyz','calculatetailcorrection')
i_all=-product(shape(psig))*kind(psig)
deallocate(psig,stat=i_stat)
call memocc(i_stat,i_all,'psig','calculatetailcorrection')
i_all=-product(shape(psigp))*kind(psigp)
deallocate(psigp,stat=i_stat)
call memocc(i_stat,i_all,'psigp','calculatetailcorrection')
i_all=-product(shape(psifscf))*kind(psifscf)
deallocate(psifscf,stat=i_stat)
call memocc(i_stat,i_all,'psifscf','calculatetailcorrection')
i_all=-product(shape(psir))*kind(psir)
deallocate(psir,stat=i_stat)
call memocc(i_stat,i_all,'psir','calculatetailcorrection')
i_all=-product(shape(psib))*kind(psib)
deallocate(psib,stat=i_stat)
call memocc(i_stat,i_all,'psib','calculatetailcorrection')
i_all=-product(shape(hpsib))*kind(hpsib)
deallocate(hpsib,stat=i_stat)
call memocc(i_stat,i_all,'hpsib','calculatetailcorrection')
i_all=-product(shape(keybg))*kind(keybg)
deallocate(keybg,stat=i_stat)
call memocc(i_stat,i_all,'keybg','calculatetailcorrection')
i_all=-product(shape(keybv))*kind(keybv)
deallocate(keybv,stat=i_stat)
call memocc(i_stat,i_all,'keybv','calculatetailcorrection')
i_all=-product(shape(ibbyz_c))*kind(ibbyz_c)
deallocate(ibbyz_c,stat=i_stat)
call memocc(i_stat,i_all,'ibbyz_c','calculatetailcorrection')
i_all=-product(shape(ibbxz_c))*kind(ibbxz_c)
deallocate(ibbxz_c,stat=i_stat)
call memocc(i_stat,i_all,'ibbxz_c','calculatetailcorrection')
i_all=-product(shape(ibbxy_c))*kind(ibbxy_c)
deallocate(ibbxy_c,stat=i_stat)
call memocc(i_stat,i_all,'ibbxy_c','calculatetailcorrection')
i_all=-product(shape(ibbyz_f))*kind(ibbyz_f)
deallocate(ibbyz_f,stat=i_stat)
call memocc(i_stat,i_all,'ibbyz_f','calculatetailcorrection')
i_all=-product(shape(ibbxz_f))*kind(ibbxz_f)
deallocate(ibbxz_f,stat=i_stat)
call memocc(i_stat,i_all,'ibbxz_f','calculatetailcorrection')
i_all=-product(shape(ibbxy_f))*kind(ibbxy_f)
deallocate(ibbxy_f,stat=i_stat)
call memocc(i_stat,i_all,'ibbxy_f','calculatetailcorrection')


if (parallel) then
   !if (iproc.eq.0) then
   !   write(*,'(1x,a,f27.14)')'Tail calculation ended'
   !endif
   allocate(wrkallred(3,2),stat=i_stat)
   call memocc(i_stat,product(shape(wrkallred))*kind(wrkallred),'wrkallred','calculatetailcorrection')
   wrkallred(1,2)=ekin_sum
   wrkallred(2,2)=epot_sum 
   wrkallred(3,2)=eproj_sum 
   call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   ekin_sum=wrkallred(1,1) 
   epot_sum=wrkallred(2,1) 
   eproj_sum=wrkallred(3,1)
   i_all=-product(shape(wrkallred))*kind(wrkallred)
   deallocate(wrkallred,stat=i_stat)
   call memocc(i_stat,i_all,'wrkallred','calculatetailcorrection')
endif

call timing(iproc,'Tail          ','OF')

end subroutine CalculateTailCorrection


subroutine transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
   implicit real*8 (a-h,o-z)
   logical, parameter :: parallel=.true.
   integer recvcount,sendcount
   dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
   real*8, allocatable :: psiw(:,:,:)
   include 'mpif.h'
 
   call timing(iproc,'Un-Transall   ','ON')
   allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
   call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','transallwaves')
 
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

   i_all=-product(shape(psiw))*kind(psiw)
   deallocate(psiw,stat=i_stat)
   call memocc(i_stat,i_all,'psiw','transallwaves')

END SUBROUTINE transallwaves

subroutine switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psiw)
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
  real(kind=8), dimension(nvctrp,norbp,nproc), intent(out) :: psiw
  !local variables
  integer :: iorb,i,j,ij

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

end subroutine switch_waves

subroutine unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,psi)
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp
  real(kind=8), dimension(nvctrp,norbp,nproc), intent(in) :: psiw
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
  !local variables
  integer :: iorb,i,j,ij

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     ij=1
     loop: do j=1,nproc
        do i=1,nvctrp
           psi(ij,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
           ij=ij+1
           if (ij.gt. nvctr_c+7*nvctr_f) exit loop
        enddo
     enddo loop
  enddo

end subroutine unswitch_waves

subroutine untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
  implicit real*8 (a-h,o-z)
  logical, parameter :: parallel=.true.
  integer recvcount,sendcount
  dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
  real*8, allocatable :: psiw(:,:,:)
  include 'mpif.h'

  call timing(iproc,'Un-Transall   ','ON')

  allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
  call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','untransallwaves')

  sendcount=nvctrp*norbp
  recvcount=nvctrp*norbp

  ! transposition: psiw(i,iorb,j,jorb) <- psit(i,iorb,jorb,j) 
  call MPI_ALLTOALL(psit,sendcount,MPI_DOUBLE_PRECISION,  &
       psiw,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  ! reformatting: psi(ij,iorb,jorb) <- psiw(i,iorb,j,jorb)
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     ij=1
     do j=1,nproc
        do i=1,nvctrp
           psi(ij,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
           ij=ij+1
           if (ij.gt. nvctr_c+7*nvctr_f) goto 333
        enddo
     enddo
333  continue
  enddo

  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw','untransallwaves')
  call timing(iproc,'Un-Transall   ','OF')


END SUBROUTINE untransallwaves

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine createIonicPotential(iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,hgrid,&
     elecfield,n1,n2,n3,n3pi,i3s,pkernel,pot_ion,eion)

  use Poisson_Solver
  
  implicit none
  integer, intent(in) :: iproc,nproc,nat,ntypes,n1,n2,n3,n3pi,i3s
  real(kind=8), intent(in) :: hgrid,elecfield
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), intent(out) :: eion
  real(kind=8), dimension(*), intent(out) :: pot_ion
  !local variables
  real(kind=8) :: hgridh,ehart,eexcu,vexcu

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------- Ionic Potential Creation'
  end if


  hgridh=0.5d0*hgrid

  ! Precalculate ionic potential from PSP charge densities and local Gaussian terms
  call input_rho_ion(iproc,nproc,ntypes,nat,iatype,rxyz,psppar, &
       & nelpsp,n1,n2,n3,n3pi,i3s,hgrid,pot_ion,eion)
  if (iproc.eq.0) write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion

  !here the value of the datacode must be kept fixed
  call PSolver('F','D',iproc,nproc,2*n1+31,2*n2+31,2*n3+31,0,hgridh,hgridh,hgridh,&
       pot_ion,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.false.)

  !print *,'ehartree',ehart
  if (n3pi > 0) then
     call addlocgauspsp(iproc,ntypes,nat,iatype,rxyz,psppar,&
          n1,n2,n3,n3pi,i3s,hgrid,pot_ion)
  end if

  !use rhopot to calculate the potential from a constant electric field along x direction
  if (elecfield /= 0.d0) then
     if (iproc.eq.0) write(*,'(1x,a,1pe10.2)') &
          'Adding constant electric field of intensity',elecfield,&
          'Ha*Bohr'

     if (n3pi > 0) call pot_constantfield(iproc,n1,n2,n3,n3pi,pot_ion,hgrid,elecfield)

  end if

end subroutine createIonicPotential


        
subroutine input_rho_ion(iproc,nproc,ntypes,nat,iatype,rxyz,psppar, &
     & nelpsp,n1,n2,n3,n3pi,i3s,hgrid,rho,eion)
  !Creates charge density arising from the ionic PSP cores
  implicit none
  include 'mpif.h'
  integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3,n3pi,i3s
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(out) :: eion
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(inout) :: rho
  !local variables
  integer :: iat,jat,i1,i2,i3,j3,ii,ix,iy,iz,i3start,i3end,ierr,ityp,jtyp,ind,i_all,i_stat
  real(kind=8) :: hgridh,pi,rholeaked,dist,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot
  real(kind=8), dimension(:), allocatable :: charges_mpi
  
  call timing(iproc,'CrtLocPot     ','ON')

  hgridh=hgrid*.5d0 
  pi=4.d0*atan(1.d0)
  ! Ionic charge 
  rholeaked=0.d0
  eion=0.d0

  if (n3pi >0 ) then
     call razero((2*n1+31)*(2*n2+31)*n3pi,rho)

     do iat=1,nat
        ityp=iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        ix=nint(rx/hgridh) 
        iy=nint(ry/hgridh) 
        iz=nint(rz/hgridh)
        !    ion-ion interaction
        do jat=1,iat-1
           dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
           jtyp=iatype(jat)
           eion=eion+nelpsp(jtyp)*nelpsp(ityp)/dist
        enddo

        rloc=psppar(0,0,ityp)
        charge=nelpsp(ityp)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc
        ii=nint(cutoff/hgridh)

        !calculate start and end of the distributed pot
        i3start=max(max(-14,iz-ii),i3s-15)
        i3end=min(min(2*n3+16,iz+ii),i3s+n3pi-16)

        do i3=iz-ii,iz+ii
           j3=i3+15-i3s+1
           do i2=iy-ii,iy+ii
              do i1=ix-ii,ix+ii
                 x=i1*hgridh-rx
                 y=i2*hgridh-ry
                 z=i3*hgridh-rz
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (i3.ge.i3start .and. i3.le.i3end  .and.  & 
                      i2.ge.-14 .and. i2.le.2*n2+16  .and.  & 
                      i1.ge.-14 .and. i1.le.2*n1+16 ) then
                    ind=i1+15+(i2+14)*(2*n1+31)+(j3-1)*(2*n1+31)*(2*n2+31)
                    rho(ind)=rho(ind)-xp*charge
                 else if (i3.lt.-14 .or. i3.gt.2*n3+16 ) then
                    rholeaked=rholeaked+xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if
  ! Check
  tt=0.d0
  do j3= 1,n3pi!i3start,i3end
     !j3=i3+15-i3s+1
     do i2= -14,2*n2+16
        do i1= -14,2*n1+16
           ind=i1+15+(i2+14)*(2*n1+31)+(j3-1)*(2*n1+31)*(2*n2+31)
           tt=tt+rho(ind)
        enddo
     enddo
  enddo

  tt=tt*hgridh**3
  rholeaked=rholeaked*hgridh**3

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt

  if (nproc > 1) then
     allocate(charges_mpi(4),stat=i_stat)
     call memocc(i_stat,product(shape(charges_mpi))*kind(charges_mpi),'charges_mpi','input_rho_ion')
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked
     call MPI_ALLREDUCE(charges_mpi(1),charges_mpi(3),2,MPI_double_precision,  &
          MPI_SUM,MPI_COMM_WORLD,ierr)
     tt_tot=charges_mpi(3)
     rholeaked_tot=charges_mpi(4)
     i_all=-product(shape(charges_mpi))*kind(charges_mpi)
     deallocate(charges_mpi,stat=i_stat)
     call memocc(i_stat,i_all,'charges_mpi','input_rho_ion')
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc.eq.0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot


  call timing(iproc,'CrtLocPot     ','OF')

end subroutine input_rho_ion

subroutine pot_constantfield(iproc,n1,n2,n3,n3pi,pot,hgrid,elecfield)
  !Creates charge density arising from the ionic PSP cores
  implicit none
  include 'mpif.h'
  integer, intent(in) :: iproc,n1,n2,n3,n3pi
  real(kind=8), intent(in) :: hgrid,elecfield
  real(kind=8), dimension(*), intent(inout) :: pot
  !local variables
  integer :: i1,i2,i3,ind
  
  call timing(iproc,'CrtLocPot     ','ON')

  do i3=1,n3pi
     do i2= -14,2*n2+16
        do i1= -14,2*n1+16
           ind=i1+15+(i2+14)*(2*n1+31)+(i3-1)*(2*n1+31)*(2*n2+31)
           pot(ind)=pot(ind)+0.25d0*elecfield*hgrid*real(i1-n1,kind=8)
        enddo
     enddo
  enddo

  call timing(iproc,'CrtLocPot     ','OF')

end subroutine pot_constantfield



subroutine addlocgauspsp(iproc,ntypes,nat,iatype,rxyz,psppar,&
     n1,n2,n3,n3pi,i3s,hgrid,pot)
  ! Add local Gaussian terms of the PSP to pot, where pot is distributed 
  implicit none
  integer, intent(in) :: ntypes,nat,n1,n2,n3,n3pi,iproc,i3s
  real(kind=8), intent(in) :: hgrid
  integer, dimension(nat), intent(in) :: iatype
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(-14:2*n1+16,-14:2*n2+16,n3pi), intent(inout) :: pot
  !local variables
  integer :: iat,i1,i2,i3,ii,ix,iy,iz,ityp,iloc,nloc,i3start,i3end,j3
  real(kind=8) :: hgridh,rloc,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  
  hgridh=hgrid*.5d0 

  do iat=1,nat
     ityp=iatype(iat)

     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     ix=nint(rx/hgridh)
     iy=nint(ry/hgridh)
     iz=nint(rz/hgridh)

     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo
     rloc=psppar(0,0,ityp)
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

     !calculate start and end of the distributed pot
     i3start=max(max(-14,iz-ii),i3s-15)
     i3end=min(min(2*n3+16,iz+ii),i3s+n3pi-16)

     do i3=i3start,i3end
        j3=i3+15-i3s+1
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
              pot(i1,i2,j3)=pot(i1,i2,j3)+xp*tt
           enddo
        enddo
     enddo

   enddo

   return
 end subroutine addlocgauspsp


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
        

subroutine applylocpotkinall(iproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,&
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     psi,pot,hpsi,epot_sum,ekin_sum)
  !  Applies the local potential and kinetic energy operator to all wavefunctions belonging to processor
  ! Input: pot,psi
  ! Output: hpsi,epot,ekin
  implicit real*8 (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension occup(norb),pot(*)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension  psi(nvctr_c+7*nvctr_f,norbp)
  dimension hpsi(nvctr_c+7*nvctr_f,norbp)
  real*8, allocatable, dimension(:) :: psifscf,psir,psig,psigp

  call timing(iproc,'ApplyLocPotKin','ON')

  ! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
  allocate(psig(8*(n1+1)*(n2+1)*(n3+1)),stat=i_stat)
  call memocc(i_stat,product(shape(psig))*kind(psig),'psig','applylocpotkinall')
  allocate(psigp(8*(n1+1)*(n2+1)*(n3+1)),stat=i_stat)
  call memocc(i_stat,product(shape(psigp))*kind(psigp),'psigp','applylocpotkinall')
  allocate(psifscf(max( (2*n1+31)*(2*n2+31)*(2*n3+16),            (2*n1+16)*(2*n2+31)*(2*n3+31))),stat=i_stat)
  call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','applylocpotkinall')
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','applylocpotkinall')

  ekin_sum=0.d0
  epot_sum=0.d0
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     
     call applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
          hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
          ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
          psig,psigp,psifscf,psir,  &
          psi(1,iorb-iproc*norbp),pot,hpsi(1,iorb-iproc*norbp),epot,ekin)

     ekin_sum=ekin_sum+occup(iorb)*ekin
     epot_sum=epot_sum+occup(iorb)*epot
     !        write(*,'(a,i0,2(1x,1pe17.10))') 'iorb,ekin,epot',iorb,ekin,epot
        
  enddo


  i_all=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=i_stat)
  call memocc(i_stat,i_all,'psig','applylocpotkinall')
  i_all=-product(shape(psigp))*kind(psigp)
  deallocate(psigp,stat=i_stat)
  call memocc(i_stat,i_all,'psigp','applylocpotkinall')
  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf','applylocpotkinall')
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','applylocpotkinall')

  call timing(iproc,'ApplyLocPotKin','OF')

END SUBROUTINE applylocpotkinall
 

        subroutine applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
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
        dimension pot((2*n1-4*nbuf+31)*(2*n2-4*nbuf+31)*(2*n3-4*nbuf+31))
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f)
        dimension hpsi(nvctr_c+7*nvctr_f)
        dimension psig(8*(n1+1)*(n2+1)*(n3+1))
        dimension psigp(8*(n1+1)*(n2+1)*(n3+1)) 
        dimension psifscf(max( (2*n1+31)*(2*n2+31)*(2*n3+16), &
                        &      (2*n1+16)*(2*n2+31)*(2*n3+31)) ) 
        dimension psir((2*n1+31)*(2*n2+31)*(2*n3+31))


        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1),psi(nvctr_c+1),psig)

        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n3+2))

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=psifscf and ww(((2*n1+31)*(2*n2+31)*(2*n3+16))
     if (nbuf.eq.0) then
        epot=0.d0
        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
          tt=pot(i)*psir(i)
          epot=epot+tt*psir(i)
          psir(i)=tt
        enddo
     else
        call psipotbuf(n1,n2,n3,nbuf,pot,psir,epot)
     endif

        call convolut_magic_t(2*n1+15,2*n2+15,2*n3+15,psir,psifscf) !psifscf=ww((2*n1+16)*(2*n2+31)*(2*n3+31))

        call analyse_shrink(n1,n2,n3,psir,psifscf,psigp)  !psir=ww(((2*n1+2)*(2*n2+16)*(2*n3+16))

         call ConvolkineticP(n1,n2,n3,  &
              nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
              psig,psigp,ekin)

        call compress(n1,n2,n3,0,n1,0,n2,0,n3,  &
                    nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psigp,hpsi(1),hpsi(nvctr_c+1))

        return
    END SUBROUTINE

     subroutine psipotbuf(n1,n2,n3,nbuf,pot,psir,epot)
     implicit real*8 (a-h,o-z)
     dimension pot(2*n1-4*nbuf+31,2*n2-4*nbuf+31,2*n3-4*nbuf+31),psir(2*n1+31,2*n2+31,2*n3+31)

! if nbuf.ne.0 n1=nb1, np1=n1
    epot=0.d0
     np1=n1-2*nbuf
     np2=n2-2*nbuf
     np3=n3-2*nbuf
!     write(*,*) 'size pot',2*n1-4*nbuf+31,2*n2-4*nbuf+31,2*n3-4*nbuf+31
!     write(*,*) 'size psir',2*n1+31,2*n2+31,2*n3+31
!     write(*,*) 'nbuf',nbuf
!     write(*,*) 'np1,np2,np3',np1,np2,np3
!     write(*,*) 'nb1,nb2,nb3',n1,n2,n3
!     write(*,*) 'limits psir',1+2*nbuf,2*np3+31+2*nbuf
     do i3=1,2*np3+31
     do i2=1,2*np2+31
     do i1=1,2*np1+31
          tt=pot(i1,i2,i3)*psir(i1+2*nbuf,i2+2*nbuf,i3+2*nbuf)
          epot=epot+tt*psir(i1+2*nbuf,i2+2*nbuf,i3+2*nbuf)
          psir(i1+2*nbuf,i2+2*nbuf,i3+2*nbuf)=tt
     enddo
     enddo
     enddo
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
        dimension x(*)
        dimension y(*)
 
!  (i1,i2*i3) -> (i2*i3,I1)
        ndat=(n2+1)*(n3+1)
        call convrot_grow(n1,ndat,x,y)
!  (i2,i3*I1) -> (i3*I1,I2)
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
        dimension x(*),y(*)

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

subroutine sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
     nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho,nrho,nscatterarr)
! Calculates the charge density by summing the square of all orbitals
! Input: psi
! Output: rho
  implicit real*8 (a-h,o-z)
  logical parallel,withmpi2
  dimension rho(nrho),occup(norb)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f,norbp)
  dimension nscatterarr(0:nproc-1,4)!n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real*8, allocatable :: psig(:,:,:,:,:,:),psifscf(:),psir(:),rho_p(:)
  include 'mpif.h'
  
  hgridh=hgrid*.5d0 
  
  ! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
  allocate(psig(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
  call memocc(i_stat,product(shape(psig))*kind(psig),'psig','sumrho')
  allocate(psifscf((2*n1+31)*(2*n2+31)*(2*n3+16)),stat=i_stat)
  call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','sumrho')
  ! Wavefunction in real space
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','sumrho')
  
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if
 
 if (parallel) then
    call timing(iproc,'Rho_comput    ','ON')
    !calculate dimensions of the complete array to be allocated before the reduction procedure
    nrhotot=0
    do jproc=0,nproc-1
       nrhotot=nrhotot+nscatterarr(jproc,1)
    end do
    allocate(rho_p((2*n1+31)*(2*n2+31)*nrhotot),stat=i_stat)
    call memocc(i_stat,product(shape(rho_p))*kind(rho_p),'rho_p','sumrho')

    !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
    call tenmminustwenty((2*n1+31)*(2*n2+31)*nrhotot,rho_p)
    !call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

    do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

       call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
            nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
            nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
            psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psig)
       call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n3+2))

       call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=ww(((2*n1+31)*(2*n2+31)*(2*n3+16))

       !sum different slices by taking into account the overlap
       i3s=0
       loop_xc_overlap: do jproc=0,nproc-1
          i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
          n3d=nscatterarr(jproc,1)
          if (n3d==0) exit loop_xc_overlap
          do i3=i3off+1,i3off+n3d
             i3s=i3s+1
             ind3=(i3-1)*(2*n1+31)*(2*n2+31)
             ind3s=(i3s-1)*(2*n1+31)*(2*n2+31)
             do i2=1,2*n2+31
                ind2=(i2-1)*(2*n1+31)+ind3
                ind2s=(i2-1)*(2*n1+31)+ind3s
                do i1=1,2*n1+31
                   ind1=i1+ind2
                   ind1s=i1+ind2s
                   !do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
                   rho_p(ind1s)=rho_p(ind1s)+(occup(iorb)/hgridh**3)*psir(ind1)**2
                end do
             end do
          end do
       end do loop_xc_overlap

       if (i3s /= nrhotot) then
          print *,'problem with rhopot array in sumrho,i3s,nrhotot,',i3s,nrhotot
          stop
       end if

    enddo

    call timing(iproc,'Rho_comput    ','OF')
    call timing(iproc,'Rho_commun    ','ON')
    call MPI_REDUCE_SCATTER(rho_p,rho,(2*n1+31)*(2*n2+31)*nscatterarr(:,1),&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    call timing(iproc,'Rho_commun    ','OF')

    ! Check
    call timing(iproc,'Rho_comput    ','ON')
    tt=0.d0
    i3off=(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)
    do i=1,(2*n1+31)*(2*n2+31)*nscatterarr(iproc,2)
       tt=tt+rho(i+i3off)
    enddo
    call timing(iproc,'Rho_comput    ','OF')
    call timing(iproc,'Rho_commun    ','ON')
    call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
         'done. Total electronic charge=',charge*hgridh**3
    call timing(iproc,'Rho_commun    ','OF')
    i_all=-product(shape(rho_p))*kind(rho_p)
    deallocate(rho_p,stat=i_stat)
    call memocc(i_stat,i_all,'rho_p','sumrho')

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
       call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n3+2))`

       call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=ww(((2*n1+31)*(2*n2+31)*(2*n3+16))

       do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
          rho(i)=rho(i)+(occup(iorb)/hgridh**3)*psir(i)**2
       enddo

    enddo
    ! Check
    tt=0.d0
    do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
       tt=tt+rho(i)
    enddo
    tt=tt*hgridh**3
    if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
         'done. Total electronic charge=',tt

    call timing(iproc,'Rho_comput    ','OF')
 endif

 i_all=-product(shape(psig))*kind(psig)
 deallocate(psig,stat=i_stat)
 call memocc(i_stat,i_all,'psig','sumrho')
 i_all=-product(shape(psifscf))*kind(psifscf)
 deallocate(psifscf,stat=i_stat)
 call memocc(i_stat,i_all,'psifscf','sumrho')
 i_all=-product(shape(psir))*kind(psir)
 deallocate(psir,stat=i_stat)
 call memocc(i_stat,i_all,'psir','sumrho')


END SUBROUTINE



        subroutine sumrho_old(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
             nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho)
! Calculates the charge density by summing the square of all orbitals
! Input: psi
! Output: rho
        implicit real*8 (a-h,o-z)
        logical parallel,withmpi2
        dimension rho((2*n1+31)*(2*n2+31)*(2*n3+31)),occup(norb)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp)
        real*8, allocatable :: psig(:,:,:,:,:,:),psifscf(:),psir(:),rho_p(:)
        include 'mpif.h'
        !flag indicating the MPI libraries used
        withmpi2=.true.

        hgridh=hgrid*.5d0 

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psig(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
        call memocc(i_stat,product(shape(psig))*kind(psig),'psig','sumrho_old')
        allocate(psifscf((2*n1+31)*(2*n2+31)*(2*n3+16)),stat=i_stat)
        call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','sumrho_old')
! Wavefunction in real space
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
        call memocc(i_stat,product(shape(psir))*kind(psir),'psir','sumrho_old')

        if (iproc==0) then
           write(*,'(1x,a)',advance='no')&
                'Calculation of charge density...'
        end if


 if (parallel) then
    if (withmpi2) then
      call timing(iproc,'Rho_comput    ','ON')
      !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
      call tenmminustwenty((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psig)
        call synthese_grow(n1,n2,n3,psir,psig,psifscf) 

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         rho(i)=rho(i)+(occup(iorb)/hgridh**3)*psir(i)**2
        enddo

      enddo

      call timing(iproc,'Rho_comput    ','OF')
      call timing(iproc,'Rho_commun    ','ON')
      call MPI_ALLREDUCE(MPI_IN_PLACE,rho,(2*n1+31)*(2*n2+31)*(2*n3+31),&
           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call timing(iproc,'Rho_commun    ','OF')

    else
      call timing(iproc,'Rho_comput    ','ON')
        allocate(rho_p((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
        call memocc(i_stat,product(shape(rho_p))*kind(rho_p),'rho_p','sumrho_old')

   !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        call tenmminustwenty((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)
        !call razero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psig)
        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n3+2))

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=ww(((2*n1+31)*(2*n2+31)*(2*n3+16))

        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         rho_p(i)=rho_p(i)+(occup(iorb)/hgridh**3)*psir(i)**2
        enddo

      enddo

      call timing(iproc,'Rho_comput    ','OF')
      call timing(iproc,'Rho_commun    ','ON')
        call MPI_ALLREDUCE(rho_p,rho,(2*n1+31)*(2*n2+31)*(2*n3+31),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call timing(iproc,'Rho_commun    ','OF')

        i_all=-product(shape(rho_p))*kind(rho_p)
        deallocate(rho_p,stat=i_stat)
        call memocc(i_stat,i_all,'rho_p','sumrho_old')
     end if
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
        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n3+2))`

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) !psifscf=ww(((2*n1+31)*(2*n2+31)*(2*n3+16))

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
        if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
             'done. Total electronic charge=',tt


        i_all=-product(shape(psig))*kind(psig)
        deallocate(psig,stat=i_stat)
        call memocc(i_stat,i_all,'psig','sumrho_old')
        i_all=-product(shape(psifscf))*kind(psifscf)
        deallocate(psifscf,stat=i_stat)
        call memocc(i_stat,i_all,'psifscf','sumrho_old')
        i_all=-product(shape(psir))*kind(psir)
        deallocate(psir,stat=i_stat)
        call memocc(i_stat,i_all,'psir','sumrho_old')


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

        allocate(wprojx(0:n1,2,nterm),stat=i_stat)
        call memocc(i_stat,product(shape(wprojx))*kind(wprojx),'wprojx','crtproj')
        allocate(wprojy(0:n2,2,nterm),stat=i_stat)
        call memocc(i_stat,product(shape(wprojy))*kind(wprojy),'wprojy','crtproj')
        allocate(wprojz(0:n3,2,nterm),stat=i_stat)
        call memocc(i_stat,product(shape(wprojz))*kind(wprojz),'wprojz','crtproj')
        allocate(work(0:nw,2),stat=i_stat)
        call memocc(i_stat,product(shape(work))*kind(work),'work','crtproj')


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
  
          i_all=-product(shape(wprojx))*kind(wprojx)
          deallocate(wprojx,stat=i_stat)
          call memocc(i_stat,i_all,'wprojx','crtproj')
          i_all=-product(shape(wprojy))*kind(wprojy)
          deallocate(wprojy,stat=i_stat)
          call memocc(i_stat,i_all,'wprojy','crtproj')
          i_all=-product(shape(wprojz))*kind(wprojz)
          deallocate(wprojz,stat=i_stat)
          call memocc(i_stat,i_all,'wprojz','crtproj')
          i_all=-product(shape(work))*kind(work)
          deallocate(work,stat=i_stat)
          call memocc(i_stat,i_all,'work','crtproj')

    return
    END SUBROUTINE
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
!Effect of orthogonality constraints on gradient 
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.true.
        dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc),occup(norb)
        allocatable :: alag(:,:,:)
        include 'mpif.h'

      call timing(iproc,'LagrM_comput  ','ON')

      allocate(alag(norb,norb,2),stat=i_stat)
      call memocc(i_stat,product(shape(alag))*kind(alag),'alag','orthoconstraint_p')
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
!        write(*,'(10(1x,1pe10.3))') (alag(iorb,jorb,1),jorb=1,norb)
!        enddo
!        endif

     scprsum=0.d0
     do iorb=1,norb
         scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
     enddo

! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psit,nvctrp,alag,norb,1.d0,hpsit,nvctrp)
     i_all=-product(shape(alag))*kind(alag)
     deallocate(alag,stat=i_stat)
     call memocc(i_stat,i_all,'alag','orthoconstraint_p')
 
      call timing(iproc,'LagrM_comput  ','OF')

     END SUBROUTINE



        subroutine orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
!Effect of orthogonality constraints on gradient 
        implicit real*8 (a-h,o-z)
        logical, parameter :: parallel=.false.
        dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp),occup(norb)
        allocatable :: alag(:,:,:)

      call timing(iproc,'LagrM_comput  ','ON')

     allocate(alag(norb,norb,2),stat=i_stat)
     call memocc(i_stat,product(shape(alag))*kind(alag),'alag','orthoconstraint')
 
!     alag(jorb,iorb,2)=+psi(k,jorb)*hpsi(k,iorb)
      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,alag(1,1,1),norb)

     scprsum=0.d0
     do iorb=1,norb
         scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
     enddo

! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psi,nvctrp,alag,norb,1.d0,hpsi,nvctrp)

     i_all=-product(shape(alag))*kind(alag)
     deallocate(alag,stat=i_stat)
     call memocc(i_stat,i_all,'alag','orthoconstraint')

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
  
!!$      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
!!$      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
!!$      czmax=czmax-eps_mach ; czmin=czmin+eps_mach
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
  END SUBROUTINE num_segkeys


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
    enddo; enddo
    if (nend.ne.nsrt) then 
       write(*,*) 'nend , nsrt',nend,nsrt
       stop 'nend <> nsrt'
    endif
    mseg=nend

    return
  END SUBROUTINE segkeys




  subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
       ntypes,iatype,rxyz,radii,rmult,hgrid,logrid)
    ! set up an array logrid(i1,i2,i3) that specifies whether the grid point
    ! i1,i2,i3 is the center of a scaling function/wavelet
    implicit real*8 (a-h,o-z)
    logical logrid
    parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
    dimension rxyz(3,nat),iatype(nat),radii(ntypes)
    dimension logrid(0:n1,0:n2,0:n3)

    do i3=nl3,nu3 
       do i2=nl2,nu2 
          do i1=nl1,nu1
             logrid(i1,i2,i3)=.false.
          enddo
       enddo
    enddo

    do iat=1,nat
       rad=radii(iatype(iat))*rmult+nbuf*hgrid
       !        write(*,*) 'iat,nat,rad',iat,nat,rad
       ml1=int(onem+(rxyz(1,iat)-rad)/hgrid)  ; mu1=int((rxyz(1,iat)+rad)/hgrid)
       ml2=int(onem+(rxyz(2,iat)-rad)/hgrid)  ; mu2=int((rxyz(2,iat)+rad)/hgrid)
       ml3=int(onem+(rxyz(3,iat)-rad)/hgrid)  ; mu3=int((rxyz(3,iat)+rad)/hgrid)

!!$       print *,'values of the mls',ml1,nl1,ml2,nl2,ml3,nl3,nbuf,rmult,rad
!!$       print *,'values of the mus',mu1,nu1,mu2,nu2,mu3,nu3,nbuf,hgrid,rxyz(1,iat)

       if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
       if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
       if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
       do i3=ml3,mu3
          dz2=(i3*hgrid-rxyz(3,iat))**2
          do i2=ml2,mu2
             dy2=(i2*hgrid-rxyz(2,iat))**2
             do i1=ml1,mu1
                dx=i1*hgrid-rxyz(1,iat)
                if (dx**2+(dy2+dz2).le.rad**2) then 
                   logrid(i1,i2,i3)=.true.
                endif
             enddo
          enddo
       enddo
    enddo

    return
  END SUBROUTINE fill_logrid


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

      if (norb.eq.1) then

!!$         if (iproc .eq. 0) then
!!$            tt=dnrm2(nvctrp,psi,1)
!!$            tt=1.d0/tt
!!$            call dscal(nvctrp,tt,psi,1)
!!$         end if

         stop 'more than one orbital needed for a parallel run'

      else

        allocate(ovrlp(norb,norb,2),stat=i_stat)
        call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','orthon_p')

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

        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp,stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp','orthon_p')

     end if

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

        allocate(ovrlp(norb,norb),stat=i_stat)
        call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','orthon')

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

      i_all=-product(shape(ovrlp))*kind(ovrlp)
      deallocate(ovrlp,stat=i_stat)
      call memocc(i_stat,i_all,'ovrlp','orthon')

 endif

       call timing(iproc,'Grams_comput  ','OF')

end subroutine orthon

subroutine MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
     rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid)
  use Poisson_Solver
  implicit none
  !Arguments
  logical, intent(in) :: output_grid
  integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb
  integer, dimension(nat), intent(in) :: iatype
  character(len=20), dimension(100), intent(in) :: atomnames
  real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  character(len=1) :: geocode
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f,norbp,nvctrp,i_all,i_stat
  integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,iat,i1,i2,i3
  real(kind=8) :: omemwf,omemker,omemden,omempot
  real(kind=8) :: tt,tmemker,tmemden,tmemps,tmemha,tminamount
  logical, dimension(:,:,:), allocatable :: logrid

! Create the file grid.ascii to visualize the grid of functions
  if (output_grid) then
     open(unit=22,file='grid.ascii',status='unknown')
     write(22,*) nat
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     do iat=1,nat
        write(22,'(3(1x,e12.5),3x,a20)') &
             rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
     enddo
  endif

  ! determine localization region for all orbitals
  allocate(logrid(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid))*kind(logrid),'logrid','memoryestimator')

  ! coarse grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,1),crmult,hgrid,logrid)
  if (output_grid) then
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
           enddo
        enddo
     end do
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid,nseg_c,nvctr_c)

  ! fine grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,2),frmult,hgrid,logrid)
  if (output_grid) then
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
           enddo
        enddo
     enddo
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid,nseg_f,nvctr_f)

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid','memoryestimator')

  tt=dble(norb)/dble(nproc)
  norbp=int((1.d0-eps_mach*tt) + tt)
  tt=dble(nvctr_c+7*nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)

  !wavefunction memory per orbitals
  omemwf=real(nvctrp*nproc*8,kind=8)
  
  geocode='F'

  if (geocode == 'P') then
     call F_FFT_dimensions(n1,n2,n3,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc)
     n01=n1
     n02=n2
     n03=n3
     tt=1.d0
  else if (geocode == 'S') then
     call S_FFT_dimensions(n1,2*n2+31,n3,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc)
     n01=n1
     n02=2*n2+31
     n03=n3
     tt=real(2*n2,kind=8)/real(n02,kind=8)
  else if (geocode == 'F') then
     call F_FFT_dimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc)
     n01=2*n1+31
     n02=2*n2+31
     n03=2*n3+31
     tt=real(8*n1*n2*n3,kind=8)/real(n01*n02*n03,kind=8)
  end if

  !density memory
  omemden=real(md1*md3*md2/nproc*8,kind=8)
  !kernel memory
  omemker=real(nd1*nd2*nd3/nproc*8,kind=8)
  !memory of full grid arrays
  omempot=real(n01*n02*n03*8,kind=8)

  write(*,'(1x,a)')&
       '------------------------------------------------------------------ Memory Estimation'
  write(*,'(1x,a,i5,a,i6,a,3(i5))')&
       'Number of atoms=',nat,' Number of orbitals=',norb,' Sim. Box Dimensions= ', n1,n2,n3
  write(*,'(1x,a,i0,a)')&
       'Estimation performed for ',nproc,' processors.'
  write(*,'(1x,a)')&
       'Memory occupation for principal arrays:'
  write(*,'(1x,a,2(i6,a3))')&
       '              Poisson Solver Kernel (K):',int(omemker/1048576.d0),'MB',&
       ceiling((omemker-aint(omemker/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  write(*,'(1x,a,2(i6,a3))')&
       '             Poisson Solver Density (D):',int(omemden/1048576.d0),'MB',&
       ceiling((omemden-aint(omemden/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  write(*,'(1x,a,2(i6,a3))')&
       '    Single Wavefunction for one orbital:',int(omemwf/1048576.d0),'MB',&
       ceiling((omemwf-aint(omemwf/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  if(nproc > 1 ) omemwf=24.d0*real(norbp*nvctrp*nproc,kind=8)
  if(nproc == 1 ) omemwf=16.d0*real(norbp*nvctrp*nproc,kind=8)
  write(*,'(1x,a,2(i6,a3))')&
       '   All Wavefunctions for each processor:',int(omemwf/1048576.d0),'MB',&
       ceiling((omemwf-aint(omemwf/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  if(nproc > 1 ) omemwf=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
  if(nproc == 1 ) omemwf=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
  write(*,'(1x,a,2(i6,a3))')&
       '      Wavefunctions + DIIS per proc (W):',int(omemwf/1048576.d0),'MB',&
       ceiling((omemwf-aint(omemwf/1048576.d0)*1048576.d0)/1024.d0),'KB'  
  write(*,'(1x,a,2(i6,a3))')&
       '   Arrays of full uncompressed grid (U):',int(omempot/1048576.d0),'MB',&
       ceiling((omempot-aint(omempot/1048576.d0)*1048576.d0)/1024.d0),'KB'  

  write(*,'(1x,a)')&
       'Estimation of Memory requirements for principal code sections:'
  write(*,'(1x,a)')&
       ' Kernel calculation | Density Construction | Poisson Solver | Hamiltonian application'
  if (nproc > 1) then 
     write(*,'(1x,a)')&
       '      ~19*K         |      W+(~4)*U+D+K    |    ~12*D+K+W   |     ~*W+(~5)*U+D+K '
     tmemker=19.d0*omemker
     tmemden=omemwf+(3.d0+tt)*omempot+omemker
     tmemps=12.d0*omemden+omemwf+omemker
     tmemha=(3.d0+2*tt)*omempot+omemwf+omemker
  else
     write(*,'(1x,a)')&
       '      ~11*K         |       ~W+(~3)*U      |     ~8*D+*W     |       ~W+(~3)*U '
     tmemker=11.d0*omemker
     tmemden=omemwf+(2.d0+tt)*omempot+omemker
     tmemps=8.d0*omemden+omemwf+omemker
     tmemha=(2.d0+2*tt)*omempot+omemwf+omemker
  end if
  write(*,'(1x,4(1x,i8,a))')&
       int(tmemker/1048576.d0),'MB         | ',int(tmemden/1048576.d0),'MB          |',&
       int(tmemps/1048576.d0),'MB     |     ',int(tmemha/1048576.d0),'MB'
  write(*,'(1x,a,i0,a)')&
       'The overall memory requirement needed for this calculation is thus: ',&
       int(max(tmemker,tmemden,tmemps,tmemha)/1048576.d0),' MB'
  tminamount=real(3*(nvctr_c+7*nvctr_f)*8,kind=8)+3.d0*n01*n02+&
       (3.d0+2*tt)*omempot
  write(*,'(1x,a)')&
       'By reducing the DIIS history and/or decreasing the number of processors the amount of'
  write(*,'(1x,a,i0,a)')&
       ' memory can be reduced but for this system it will never be less than ',&
       int(tminamount/1048576.d0),' MB'

end subroutine MemoryEstimator


subroutine createWavefunctionsDescriptors(iproc,nproc,idsx,n1,n2,n3,output_grid,&
     & hgrid,nat,ntypes,iatype,atomnames,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
     & ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp,&
     & keyg,keyv,norb,norbp)
!calculates the descriptor arrays keyg and keyv as well as nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp
!calculates also the arrays ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f needed for convolut_standard
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,idsx,n1,n2,n3,nat,ntypes,norb
  integer, intent(out) :: nseg_c,nseg_f,nvctr_c,nvctr_f
  integer, intent(out) :: norbp,nvctrp
  logical, intent(in) :: output_grid
  integer, intent(in) :: iatype(nat)
  real*8, intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
  integer, intent(out) :: ibyz_c(2,0:n2,0:n3), ibxz_c(2,0:n1,0:n3), ibxy_c(2,0:n1,0:n2)
  integer, intent(out) :: ibyz_f(2,0:n2,0:n3), ibxz_f(2,0:n1,0:n3), ibxy_f(2,0:n1,0:n2)
  real*8 :: rxyz(3, nat), radii_cf(ntypes, 2)
  character(len=20), intent(in) :: atomnames(100)
  integer, pointer :: keyg(:,:), keyv(:)

  !Local variables
  real*8, parameter :: eps_mach=1.d-12,onem=1.d0-eps_mach
  integer :: iat,i1,i2,i3,norbme,norbyou,jpst,jproc,i_all,i_stat
  real*8 :: tt
  logical, allocatable :: logrid_c(:,:,:), logrid_f(:,:,:)

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------- Wavefunctions Descriptors Creation'
  end if


  call timing(iproc,'CrtDescriptors','ON')

! Create the file grid.ascii to visualize the grid of functions
  if (iproc.eq.0 .and. output_grid) then
     open(unit=22,file='grid.ascii',status='unknown')
     write(22,*) nat
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     do iat=1,nat
        write(22,'(3(1x,e12.5),3x,a20)') rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
        !write(*,'(3(1x,e12.5),3x,a20)') rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
     enddo
  endif

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_c))*kind(logrid_c),'logrid_c','createwavefunctionsdescriptors')
  allocate(logrid_f(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_f))*kind(logrid_f),'logrid_f','createwavefunctionsdescriptors')

  ! coarse grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,1),crmult,hgrid,logrid_c)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
           enddo
        enddo
     end do
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,nvctr_c)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       'Coarse resolution grid: Number of segments= ',nseg_c,'points=',nvctr_c
  !if (iproc.eq.0) write(*,'(1x,a,2(1x,i10))') &
  !     'orbitals have coarse segment, elements',nseg_c,nvctr_c
  call bounds(n1,n2,n3,logrid_c,ibyz_c,ibxz_c,ibxy_c)

  ! fine grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,2),frmult,hgrid,logrid_f)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
           enddo
        enddo
     enddo
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,nvctr_f)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       '  Fine resolution grid: Number of segments= ',nseg_f,'points=',nvctr_f
  !if (iproc.eq.0) write(*,'(1x,a,2(1x,i10))') &
  !     'orbitals have fine   segment, elements',nseg_f,7*nvctr_f
  call bounds(n1,n2,n3,logrid_f,ibyz_f,ibxz_f,ibxy_f)

  if (iproc.eq.0 .and. output_grid) close(22)

  ! allocations for arrays holding the wavefunctions and their data descriptors
  allocate(keyg(2,nseg_c+nseg_f),stat=i_stat)
  call memocc(i_stat,product(shape(keyg))*kind(keyg),'keyg','createwavefunctionsdescriptors')
  allocate(keyv(nseg_c+nseg_f),stat=i_stat)
  call memocc(i_stat,product(shape(keyv))*kind(keyv),'keyv','createwavefunctionsdescriptors')

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,keyg(:,1:nseg_c),keyv(1:nseg_c))

  ! fine grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,keyg(:,nseg_c+1:nseg_c+nseg_f), &
    & keyv(nseg_c+1:nseg_c+nseg_f))

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c','createwavefunctionsdescriptors')
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f','createwavefunctionsdescriptors')

! allocate wavefunction arrays
  tt=dble(norb)/dble(nproc)
  norbp=int((1.d0-eps_mach*tt) + tt)
  !if (iproc.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp
  norbme=max(min((iproc+1)*norbp,norb)-iproc*norbp,0)
  !write(*,'(a,i0,a,i0,a)') '- iproc ',iproc,' treats ',norbme,' orbitals '
  if (iproc == 0 .and. nproc>1) then
     jpst=0
     do jproc=0,nproc-2
        norbme=max(min((jproc+1)*norbp,norb)-jproc*norbp,0)
        norbyou=max(min((jproc+2)*norbp,norb)-(jproc+1)*norbp,0)
        if (norbme /= norbyou) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' orbitals '
           jpst=jproc+1
        end if
     end do
     write(*,'(3(a,i0),a)')&
          ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' orbitals '
  end if

  tt=dble(nvctr_c+7*nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)

  if (iproc.eq.0) write(*,'(1x,a,i0)') &
       'Wavefunction memory occupation per orbital (Bytes): ',&
       nvctrp*nproc*8

  call timing(iproc,'CrtDescriptors','OF')

END SUBROUTINE createWavefunctionsDescriptors

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

    if (iproc.eq.0) then
       write(*,'(1x,a)')&
         '------------------------------------------------------------ PSP Projectors Creation'
       write(*,'(1x,a4,4x,a4,1x,a)')&
            'Atom','Name','Number of projectors'
    end if

    call timing(iproc,'CrtProjectors ','ON')


    ! determine localization region for all projectors, but do not yet fill the descriptor arrays
    allocate(logrid(0:n1,0:n2,0:n3),stat=i_stat)
    call memocc(i_stat,product(shape(logrid))*kind(logrid),'logrid','createprojectorsarrays')

    nseg_p(0)=0 
    nvctr_p(0)=0 

    istart=1
    nproj=0
    do iat=1,nat

       call numb_proj(iatype(iat),ntypes,psppar,npspcode,mproj)
       if (mproj.ne.0) then 

          if (iproc.eq.0) write(*,'(1x,i4,2x,a6,1x,i20)')&
               iat,trim(atomnames(iatype(iat))),mproj


          !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))')&
          !     'projector descriptors for atom with mproj ',iat,mproj
          nproj=nproj+mproj

          ! coarse grid quantities
          call  pregion_size(rxyz(1,iat),radii_cf(1,2),cpmult,iatype(iat),ntypes, &
               hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
          !if (iproc.eq.0) write(*,'(a,6(i4))') 'coarse grid',nl1,nu1,nl2,nu2,nl3,nu3
          nboxp_c(1,1,iat)=nl1 ; nboxp_c(2,1,iat)=nu1
          nboxp_c(1,2,iat)=nl2 ; nboxp_c(2,2,iat)=nu2
          nboxp_c(1,3,iat)=nl3 ; nboxp_c(2,3,iat)=nu3
          call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)
          call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
          !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr,coarse projectors ',mseg,mvctr
          nseg_p(2*iat-1)=nseg_p(2*iat-2) + mseg
          nvctr_p(2*iat-1)=nvctr_p(2*iat-2) + mvctr
          istart=istart+mvctr*mproj

          ! fine grid quantities
          call  pregion_size(rxyz(1,iat),radii_cf(1,2),fpmult,iatype(iat),ntypes, &
               hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
          !if (iproc.eq.0) write(*,'(a,6(i4))') 'fine   grid',nl1,nu1,nl2,nu2,nl3,nu3
          nboxp_f(1,1,iat)=nl1 ; nboxp_f(2,1,iat)=nu1
          nboxp_f(1,2,iat)=nl2 ; nboxp_f(2,2,iat)=nu2
          nboxp_f(1,3,iat)=nl3 ; nboxp_f(2,3,iat)=nu3
          call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
          call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
          !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr
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

    if (iproc.eq.0) then
       write(*,'(28x,a)') '------'
       write(*,'(1x,a,i5)') 'Total number of projectors =',nproj
    end if

    ! allocations for arrays holding the projectors and their data descriptors
    allocate(keyg_p(2,nseg_p(2*nat)),stat=i_stat)
    call memocc(i_stat,product(shape(keyg_p))*kind(keyg_p),'keyg_p','createprojectorsarrays')
    allocate(keyv_p(nseg_p(2*nat)),stat=i_stat)
    call memocc(i_stat,product(shape(keyv_p))*kind(keyv_p),'keyv_p','createprojectorsarrays')
    nprojel=istart-1
    allocate(proj(nprojel),stat=i_stat)
    call memocc(i_stat,product(shape(proj))*kind(proj),'proj','createprojectorsarrays')


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

    !if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


    if (iproc.eq.0) write(*,'(1x,a)',advance='no') &
         'Calculating wavelets expansion of projectors...'
    !allocate these vectors up to the maximum size we can get
    nterm_max=10 !if GTH nterm_max=3
    allocate(fac_arr(nterm_max),stat=i_stat)
    call memocc(i_stat,product(shape(fac_arr))*kind(fac_arr),'fac_arr','createprojectorsarrays')
    allocate(lx(nterm_max),stat=i_stat)
    call memocc(i_stat,product(shape(lx))*kind(lx),'lx','createprojectorsarrays')
    allocate(ly(nterm_max),stat=i_stat)
    call memocc(i_stat,product(shape(ly))*kind(ly),'ly','createprojectorsarrays')
    allocate(lz(nterm_max),stat=i_stat)
    call memocc(i_stat,product(shape(lz))*kind(lz),'lz','createprojectorsarrays')

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
                      print *,'norm projector for atom ',trim(atomnames(iatype(iat))),&
                           'iproc,l,i,rl,scpr=',iproc,l,i,gau_a,scpr
                      stop 'norm projector'
                   end if

!!$                   !plot the p projector
!!$                   if (l==2 .and. m==1) then
!!$                      mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
!!$                      mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
!!$                      jseg_c=nseg_p(2*iat-2)+1
!!$                      call plot_wf(51,n1,n2,n3,hgrid,mbseg_c,mvctr_c,&
!!$                           keyg_p(:,jseg_c:jseg_c+mbseg_c+mbseg_f-1),&
!!$                           keyv_p(jseg_c:jseg_c+mbseg_c+mbseg_f-1),mbseg_f,nvctr_f, &
!!$                           rx,ry,rz,proj(istart_c:istart_f+7*mvctr_f-1))
!!$                   end if

                   ! testing end
                   istart_c=istart_f+7*mvctr_f
                   if (istart_c.gt.istart) stop 'istart_c > istart'

                   !do iterm=1,nterm
                   !   if (iproc.eq.0) write(*,'(1x,a,i0,1x,a,1pe10.3,3(1x,i0))') &
                   !        'projector: iat,atomname,gau_a,lx,ly,lz ', & 
                   !        iat,trim(atomnames(iatype(iat))),gau_a,lx(iterm),ly(iterm),lz(iterm)
                   !enddo


                enddo
             endif
          enddo
       enddo
    enddo
    if (iproj.ne.nproj) stop 'incorrect number of projectors created'
    ! projector part finished
    if (iproc ==0) write(*,'(1x,a)')'done.'

    i_all=-product(shape(logrid))*kind(logrid)
    deallocate(logrid,stat=i_stat)
    call memocc(i_stat,i_all,'logrid','createprojectorsarrays')
    i_all=-product(shape(fac_arr))*kind(fac_arr)
    deallocate(fac_arr,stat=i_stat)
    call memocc(i_stat,i_all,'fac_arr','createprojectorsarrays')
    i_all=-product(shape(lx))*kind(lx)
    deallocate(lx,stat=i_stat)
    call memocc(i_stat,i_all,'lx','createprojectorsarrays')
    i_all=-product(shape(ly))*kind(ly)
    deallocate(ly,stat=i_stat)
    call memocc(i_stat,i_all,'ly','createprojectorsarrays')
    i_all=-product(shape(lz))*kind(lz)
    deallocate(lz,stat=i_stat)
    call memocc(i_stat,i_all,'lz','createprojectorsarrays')
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

        allocate(ovrlp(norb,norb,3),stat=i_stat)
        call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','loewe_p')
        allocate(evall(norb),stat=i_stat)
        call memocc(i_stat,product(shape(evall))*kind(evall),'evall','loewe_p')

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

        allocate(psitt(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psitt))*kind(psitt),'psitt','loewe_p')
! new eigenvectors
!   psitt(i,iorb)=psit(i,jorb)*ovrlp(jorb,iorb,3)
      call DGEMM('N','N',nvctrp,norb,norb,1.d0,psit,nvctrp,ovrlp(1,1,3),norb,0.d0,psitt,nvctrp)
      call DCOPY(nvctrp*norbp*nproc,psitt,1,psit,1)
      i_all=-product(shape(psitt))*kind(psitt)
      deallocate(psitt,stat=i_stat)
      call memocc(i_stat,i_all,'psitt','loewe_p')

        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp,stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp','loewe_p')
        i_all=-product(shape(evall))*kind(evall)
        deallocate(evall,stat=i_stat)
        call memocc(i_stat,i_all,'evall','loewe_p')

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

        allocate(ovrlp(norb,norb,3),stat=i_stat)
        call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','loewe')
        allocate(evall(norb),stat=i_stat)
        call memocc(i_stat,product(shape(evall))*kind(evall),'evall','loewe')

! Overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
        call DSYRK('U','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp(1,1,1),norb)

!       write(*,*) 'OVERLAP'
!       do i=1,norb
!       write(*,'(10(1x,1pe17.10))') (ovrlp(i,j,1),j=1,norb)
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
      allocate(tpsi(nvctrp,norb),stat=i_stat)
      call memocc(i_stat,product(shape(tpsi))*kind(tpsi),'tpsi','loewe')
!   tpsi(i,iorb)=psi(i,jorb)*ovrlp(jorb,iorb,3)
      call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi,nvctrp,ovrlp(1,1,3),norb,0.d0,tpsi,nvctrp)
      call DCOPY(nvctrp*norb,tpsi,1,psi,1)
      i_all=-product(shape(tpsi))*kind(tpsi)
      deallocate(tpsi,stat=i_stat)
      call memocc(i_stat,i_all,'tpsi','loewe')

        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp,stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp','loewe')
        i_all=-product(shape(evall))*kind(evall)
        deallocate(evall,stat=i_stat)
        call memocc(i_stat,i_all,'evall','loewe')

 endif

        END SUBROUTINE


        subroutine checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
        implicit real*8 (a-h,o-z)
        dimension psit(nvctrp,norbp*nproc)
        real*8, allocatable :: ovrlp(:,:,:)
        include 'mpif.h'

        allocate(ovrlp(norb,norb,2),stat=i_stat)
        call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','checkortho_p')

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
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,'(1x,a,3(1x,i0))') 'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)       write(*,'(1x,a,3(1x,i0))') 'ERROR ORTHO',iorb,jorb,scpr
110     continue

        if (dev.gt.1.d-10) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',iproc,dev

        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp,stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp','checkortho_p')

        return
        END SUBROUTINE


        subroutine checkortho(norb,norbp,nvctrp,psi)
        implicit real*8 (a-h,o-z)
        dimension psi(nvctrp,norbp)
        real*8, allocatable :: ovrlp(:,:,:)

        allocate(ovrlp(norb,norb,1),stat=i_stat)
        call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','checkortho')

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
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,'(1x,a,3(1x,i0))') 'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)       write(*,'(1x,a,3(1x,i0))') 'ERROR ORTHO',iorb,jorb,scpr
     enddo  ; enddo

        if (dev.gt.1.d-10) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',0,dev

        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp,stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp','checkortho')


        return
        END SUBROUTINE

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,nzatom,nelpsp,&
     & psppar,npspcode,norbe,norbsc,atomnames,ntypes,iatype,iasctype,nat,natsc,&
     & scorb,norbsc_arr)
  implicit none
! character(len = *), intent(in) :: filename
  integer, intent(in) :: ngx, iproc, ntypes
  integer, intent(in) :: nzatom(ntypes), nelpsp(ntypes)
  real*8, intent(in) :: psppar(0:4,0:4,ntypes)
  integer, intent(in) :: npspcode(ntypes),iasctype(ntypes)
  real*8, intent(out) :: xp(ngx, ntypes), psiat(ngx, 5, ntypes), occupat(5, ntypes)
  integer, intent(out) :: ng(ntypes), nl(4,ntypes)
  character(len = 20), intent(in) :: atomnames(100)
  integer, intent(out) :: norbe,norbsc
  integer, intent(in) :: nat,natsc
  integer, intent(in) :: iatype(nat)
  logical, dimension(4,natsc), intent(out) :: scorb
  integer, dimension(natsc+1), intent(out) :: norbsc_arr

  character(len = 20) :: pspatomname
  logical :: exists,found
  integer :: ity,i,j,l,ipsp,ifile,ng_fake,ierror,iatsc,iat,ipow,lsc,inorbsc,iorbsc_count
  real(kind=8) :: sccode

! Read the data file.
  nl(1:4,1:ntypes) = 0
  ng(1:ntypes) = 0
  xp(1:ngx,1:ntypes) = 0.d0
  psiat(1:ngx,1:5,1:ntypes) = 0.d0
  occupat(1:5,1:ntypes)= 0.d0

! Test if the file 'inguess.dat exists
  inquire(file='inguess.dat',exist=exists)
  if (exists) then
     open(unit=24,file='inguess.dat',form='formatted',action='read',status='old')
  end if

  loop_assign: do ity=1,ntypes

     if (exists) then
        rewind(24)
     end if
     found = .false.

     loop_find: do
        if (.not.exists) then
!          The file 'inguess.dat' does not exist: automatic generation
           exit loop_find
        end if
        read(24,'(a)',iostat=ierror) pspatomname
        if (ierror /= 0) then
!          Read error or end of file
           exit loop_find
        end if

        if (pspatomname .eq. atomnames(ity)) then
           if (iproc.eq.0) then
              write(*,'(1x,a,a,a)') 'input wavefunction data for atom ',trim(atomnames(ity)),&
                   ' found'
           end if
           found = .true.
           read(24,*) nl(1,ity),(occupat(i,ity),i=1,nl(1,ity)),  &
                nl(2,ity),(occupat(i,ity),i=1+nl(1,ity),nl(2,ity)+nl(1,ity)) ,&
                nl(3,ity),(occupat(i,ity),i=1+nl(2,ity)+nl(1,ity),nl(3,ity)+nl(2,ity)+nl(1,ity)) ,&
                nl(4,ity),(occupat(i,ity),&
                i=1+nl(3,ity)+nl(2,ity)+nl(1,ity),nl(4,ity)+nl(3,ity)+nl(2,ity)+nl(1,ity))
           !print *,nl(:,ity),occupat(:,ity)
           if (nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity).gt.5) then
              print *,'error: number of valence orbitals too big'
              print *,nl(:,ity),occupat(:,ity)
              stop
           end if
           read(24,*) ng(ity)
           !print *, pspatomnames(ity),(nl(l,ity),l=1,4),ng(ity),ngx,npsp
           if (ng(ity).gt.ngx) stop 'enlarge ngx'
           !read(24,'(30(e12.5))') (xp(i,ity)  ,i=1,ng(ity))
           read(24,*) (xp(i,ity)  ,i=1,ng(ity))
           do i=1,ng(ity) 
              read(24,*) (psiat(i,j,ity),j=1,nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity))
           enddo

           exit loop_find

        else

           read(24,*)
           read(24,*)ng_fake
           read(24,*) 
           do i=1,ng_fake
              read(24,*)
           enddo

        end if

     enddo loop_find

     if (.not.found) then

        if (iproc.eq.0) then
           write(*,'(1x,a,a6,a)',advance='no')&
                'Input wavefunction data for atom ',trim(atomnames(ity)),&
                ' NOT found, automatic generation...'
        end if

        !the default value for the gaussians is chosen to be 21
        ng(ity)=21
        call iguess_generator(iproc, nzatom(ity), nelpsp(ity),psppar(0,0,ity),npspcode(ity),&
             ng(ity)-1,nl(1,ity),5,occupat(1:5,ity),xp(1:ng(ity),ity),psiat(1:ng(ity),1:5,ity))

!values obtained from the input guess generator in iguess.dat format
        !write these values on a file
        if (iproc .eq. 0) then
           open(unit=12,file='inguess.new',status='unknown')

           !write(*,*)' --------COPY THESE VALUES INSIDE inguess.dat--------'
           write(12,*)trim(atomnames(ity))//' (remove _lda)'
           write(12,*)nl(1,ity),(occupat(i,ity),i=1,nl(1,ity)),&
                nl(2,ity),(occupat(i+nl(1,ity),ity),i=1,nl(2,ity)),&
                nl(3,ity),(occupat(i+nl(1,ity)+nl(2,ity),ity),i=1,nl(3,ity)),&
                nl(4,ity),(occupat(i+nl(1,ity)+nl(2,ity)+nl(3,ity),ity),i=1,nl(4,ity))
           write(12,*)ng(ity)
           write(12,'(30(e12.5))')xp(1:ng(ity),ity)
           do j=1,ng(ity)
              write(12,*)(psiat(j,i,ity),i=1,nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity))
           end do
           !print *,' --------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--------'
           write(*,'(1x,a)')'done.'
        end if

     end if


  end do loop_assign

  close(unit=24)

  ! number of orbitals, total and semicore
  norbe=0
  norbsc=0
  iatsc=0
  do iat=1,nat
     ity=iatype(iat)
     norbe=norbe+nl(1,ity)+3*nl(2,ity)+5*nl(3,ity)+7*nl(4,ity)
      if (iasctype(ity)/=0) then !the atom has some semicore orbitals
        iatsc=iatsc+1
        sccode=real(iasctype(ity),kind=8)
        inorbsc=ceiling(dlog(sccode)/dlog(10.d0))
        if (sccode==1.d0) inorbsc=1
        ipow=inorbsc-1
        scorb(:,iatsc)=.false.
        !count the semicore orbitals for this atom
        iorbsc_count=0
        do i=1,inorbsc
           lsc=floor(sccode/10.d0**ipow)
           iorbsc_count=iorbsc_count+(2*lsc-1)
           scorb(lsc,iatsc)=.true.
           sccode=sccode-real(lsc,kind=8)*10.d0**ipow
           ipow=ipow-1
           !print *,iasctype(ity),inorbsc,lsc
        end do
        norbsc_arr(iatsc)=iorbsc_count
        norbsc=norbsc+iorbsc_count
     end if
  end do
  !orbitals which are non semicore
  norbsc_arr(natsc+1)=norbe-norbsc

  if (iproc ==0) then
     write(*,'(1x,a,i0,a)')'Generating ',norbe,' Atomic Input Orbitals'
     if (norbsc /=0)   write(*,'(1x,a,i0,a)')'  of which ',norbsc,' are semicore orbitals'
  end if


END SUBROUTINE 

subroutine input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nat,natsc,norb,norbp,n1,n2,n3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
     rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
     atomnames,ntypes,iatype,iasctype,pkernel,nzatom,nelpsp,psppar,npspcode,ixc,&
     ppsi,ppsit,eval,accurex,datacode,nscatterarr,ngatherarr)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors writes its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave

  use Poisson_Solver

  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  character(len=20), dimension(100), intent(in) :: atomnames
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,nat,natsc,ntypes,norb,norbp,n1,n2,n3,nprojel,nproj,ixc
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp,nseg_c,nseg_f
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(out) :: accurex
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: iasctype,npspcode,nzatom,nelpsp
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
  integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
  integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(nprojel), intent(in) :: proj
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:), pointer :: ppsi,ppsit
  !real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: ppsi
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer, parameter :: ngx=31
  integer :: i,iorb,iorbsc,imatrsc,iorbst,imatrst,i_stat,i_all,ierr,info,jproc,jpst,norbeyou
  integer :: norbe,norbep,norbi,norbj,norbeme,ndim_hamovr,n_lp,norbi_max,norbsc
  real(kind=8) :: hgridh,tt,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum
  logical, dimension(:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: norbsc_arr,ng
  integer, dimension(:,:), allocatable :: nl
  real(kind=8), dimension(:), allocatable :: work_lp,pot,evale,occupe
  real(kind=8), dimension(:,:), allocatable :: xp,occupat,hamovr,psi,hpsi
  real(kind=8), dimension(:,:,:), allocatable :: psiw,psiat

  allocate(xp(ngx,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(xp))*kind(xp),'xp','input_wf_diag')
  allocate(psiat(ngx,5,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(psiat))*kind(psiat),'psiat','input_wf_diag')
  allocate(occupat(5,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(occupat))*kind(occupat),'occupat','input_wf_diag')
  allocate(ng(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(ng))*kind(ng),'ng','input_wf_diag')
  allocate(nl(4,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nl))*kind(nl),'nl','input_wf_diag')
  allocate(scorb(4,natsc),stat=i_stat)
  call memocc(i_stat,product(shape(scorb))*kind(scorb),'scorb','input_wf_diag')

  allocate(norbsc_arr(natsc+1),stat=i_stat)
  call memocc(i_stat,product(shape(norbsc_arr))*kind(norbsc_arr),'norbsc_arr','input_wf_diag')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  ! Read the inguess.dat file or generate the input guess via the inguess_generator
  call readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,nzatom,nelpsp,psppar,&
       & npspcode,norbe,norbsc,atomnames,ntypes,iatype,iasctype,nat,natsc,scorb,&
       & norbsc_arr)

  !  allocate wavefunctions and their occupation numbers
  allocate(occupe(norbe),stat=i_stat)
  call memocc(i_stat,product(shape(occupe))*kind(occupe),'occupe','input_wf_diag')
  tt=dble(norbe)/dble(nproc)
  norbep=int((1.d0-eps_mach*tt) + tt)

  if (iproc == 0 .and. nproc>1) then
     jpst=0
     do jproc=0,nproc-2
        norbeme=max(min((jproc+1)*norbep,norbe)-jproc*norbep,0)
        norbeyou=max(min((jproc+2)*norbep,norbe)-(jproc+1)*norbep,0)
        if (norbeme /= norbeyou) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbeme,' inguess orbitals '
           jpst=jproc+1
        end if
     end do
     write(*,'(3(a,i0),a)')&
          ' Processes from ',jpst,' to ',nproc-1,' treat ',norbeyou,' inguess orbitals '
  end if

  hgridh=.5d0*hgrid

  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(psi(nvctrp,norbep*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','input_wf_diag')
  else
     allocate(psi(nvctr_c+7*nvctr_f,norbep),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','input_wf_diag')
  end if

  ! Create input guess orbitals
  call createAtomicOrbitals(iproc,nproc,atomnames,&
       & nat,rxyz,norbe,norbep,norbsc,occupe,occupat,ngx,xp,psiat,ng,nl,&
       & nvctr_c,nvctr_f,n1,n2,n3,hgrid,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       & nseg_c,nseg_f,keyg,keyv,iatype,ntypes,iasctype,natsc,psi,eks,scorb)

!!$  !plot the initial LCAO wavefunctions
!!$  do i=2*iproc+1,2*iproc+2
!!$     iounit=15+3*(i-1)
!!$     print *,'iounit',iounit,'-',iounit+2
!!$     call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$          rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc))
!!$  end do


  i_all=-product(shape(scorb))*kind(scorb)
  deallocate(scorb,stat=i_stat)
  call memocc(i_stat,i_all,'scorb','input_wf_diag')
  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp','input_wf_diag')
  i_all=-product(shape(psiat))*kind(psiat)
  deallocate(psiat,stat=i_stat)
  call memocc(i_stat,i_all,'psiat','input_wf_diag')
  i_all=-product(shape(occupat))*kind(occupat)
  deallocate(occupat,stat=i_stat)
  call memocc(i_stat,i_all,'occupat','input_wf_diag')
  i_all=-product(shape(ng))*kind(ng)
  deallocate(ng,stat=i_stat)
  call memocc(i_stat,i_all,'ng','input_wf_diag')
  i_all=-product(shape(nl))*kind(nl)
  deallocate(nl,stat=i_stat)
  call memocc(i_stat,i_all,'nl','input_wf_diag')


  ! resulting charge density and potential
  if (datacode=='G') then
     call sumrho_old(parallel,iproc,norbe,norbep,n1,n2,n3,hgrid,occupe,  & 
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)
  else
     call sumrho(parallel,iproc,nproc,norbe,norbep,n1,n2,n3,hgrid,occupe,  & 
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot,&
          (2*n1+31)*(2*n2+31)*nscatterarr(iproc,1),nscatterarr)
  end if
  !      ixc=1   ! LDA functional
  call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.)


  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(hpsi(nvctrp,norbep*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','input_wf_diag')
  else
     allocate(hpsi(nvctr_c+7*nvctr_f,norbep),stat=i_stat)
     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','input_wf_diag')
  end if


  call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
       psppar,npspcode,norbe,norbep,occupe,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
       nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,ngatherarr,nscatterarr(iproc,2),&
       rhopot(1+(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum)

  i_all=-product(shape(occupe))*kind(occupe)
  deallocate(occupe,stat=i_stat)
  call memocc(i_stat,i_all,'occupe','input_wf_diag')


  accurex=abs(eks-ekin_sum)
  if (iproc.eq.0) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks


  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted

  if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
       'Input Wavefunctions Orthogonalization:'

  norbi_max=maxval(norbsc_arr)
  
  !calculate the dimension of the overlap matrix
  ndim_hamovr=0
  do i=1,natsc+1
     ndim_hamovr=ndim_hamovr+norbsc_arr(i)**2
  end do

  if (parallel) then

     !transpose all the wavefunctions for having a piece of all the orbitals 
     !for each processor
     call timing(iproc,'Un-Transall   ','ON')
     allocate(psiw(nvctrp,norbep,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','input_wf_diag')

     call switch_waves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,psi,psiw)
     call MPI_ALLTOALL(psiw,nvctrp*norbep,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbep,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     call switch_waves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,hpsi,psiw)
     call MPI_ALLTOALL(psiw,nvctrp*norbep,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbep,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw','input_wf_diag')
     call timing(iproc,'Un-Transall   ','OF')
     !end of transposition

     allocate(hamovr(ndim_hamovr,4),stat=i_stat)
     call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','input_wf_diag')

     !calculate the overlap matrix for each group of the semicore atoms
!       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
!       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
     iorbst=1
     imatrst=1
     !print *,'norbi',norbi,natsc,norbsc_arr(natsc+1)

     if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
          'Overlap Matrix...'

     do i=1,natsc
        norbi=norbsc_arr(i)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,3),norbi)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,4),norbi)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
     norbi=norbsc_arr(natsc+1)
     call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
          0.d0,hamovr(imatrst,3),norbi)

     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi','input_wf_diag')

     call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
          0.d0,hamovr(imatrst,4),norbi)
     
     !reduce the overlap matrix between all the processors
     call MPI_ALLREDUCE(hamovr(1,3),hamovr(1,1),2*ndim_hamovr,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     !found the eigenfunctions for each group
     n_lp=max(10,4*norbi_max)
     allocate(work_lp(n_lp),stat=i_stat)
     call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','input_wf_diag')
     allocate(evale(norbi_max),stat=i_stat)
     call memocc(i_stat,product(shape(evale))*kind(evale),'evale','input_wf_diag')
     
     if (iproc.eq.0) write(*,'(1x,a)')'Linear Algebra...'

     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DSYGV(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
             norbi,evale,work_lp,n_lp,info)

        if (info.ne.0) write(*,*) 'DSYGV ERROR',info,i,natsc+1
!!$        !write the matrices on a file
!!$        open(33+2*(i-1))
!!$        do jjorb=1,norbi
!!$           write(33+2*(i-1),'(2000(1pe10.2))')&
!!$                (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        end do
!!$        close(33+2*(i-1))
!!$        open(34+2*(i-1))
!!$        do jjorb=1,norbi
!!$           write(34+2*(i-1),'(2000(1pe10.2))')&
!!$                (hamovr(imatrst-1+jjorb+(jiorb-1)*norbi,1),jiorb=1,norbi)
!!$        end do
!!$        close(34+2*(i-1))

        if (iproc.eq.0) then
        do iorb=1,norbi
        write(*,'(1x,a,i0,a,1x,1pe21.14)') 'evale(',iorb+iorbst-1,')=',evale(iorb)
        enddo
        endif
        do iorb=iorbst,min(norbi+iorbst-1,norb)
           eval(iorb)=evale(iorb-iorbst+1)
        enddo
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do

     i_all=-product(shape(work_lp))*kind(work_lp)
     deallocate(work_lp,stat=i_stat)
     call memocc(i_stat,i_all,'work_lp','input_wf_diag')
     i_all=-product(shape(evale))*kind(evale)
     deallocate(evale,stat=i_stat)
     call memocc(i_stat,i_all,'evale','input_wf_diag')

     !allocate the transposed wavefunction
     allocate(ppsit(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(ppsit))*kind(ppsit),'ppsit','input_wf_diag')

     if (iproc.eq.0) write(*,'(1x,a)',advance='no')'Building orthogonal Input Wavefunctions...'

     !perform the vector-matrix multiplication for building the input wavefunctions
     ! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
     iorbst=1
     imatrst=1
     do i=1,natsc
        norbi=norbsc_arr(i)
        call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
             hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst),nvctrp)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
     norbi=norbsc_arr(natsc+1)
     norbj=norb-norbsc
     call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
          hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst),nvctrp)

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi','input_wf_diag')
     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr','input_wf_diag')
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr','input_wf_diag')
     
     !retranspose the wavefunctions
     call timing(iproc,'Un-Transall   ','ON')
     allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','input_wf_diag')

     call MPI_ALLTOALL(ppsit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psiw,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!!$     i_all=-product(shape(ppsit))*kind(ppsit)
!!$     deallocate(ppsit,stat=i_stat)
!!$     call memocc(i_stat,i_all,'ppsit','input_wf_diag')

     !allocate the direct wavefunction
     allocate(ppsi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','input_wf_diag')

     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,ppsi)

     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw','input_wf_diag')

     call timing(iproc,'Un-Transall   ','OF')

  else !serial case

     write(*,'(1x,a)',advance='no')'Overlap Matrix...'

     allocate(hamovr(ndim_hamovr,4),stat=i_stat)
     call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','input_wf_diag')
     !hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,1),norbi)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,2),norbi)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do

     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi','input_wf_diag')


     n_lp=max(10,4*norbi_max)
     allocate(work_lp(n_lp),stat=i_stat)
     call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','input_wf_diag')
     allocate(evale(norbe),stat=i_stat)
     call memocc(i_stat,product(shape(evale))*kind(evale),'evale','input_wf_diag')

     write(*,'(1x,a)')'Linear Algebra...'
     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DSYGV(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
             norbi,evale(iorbst),work_lp,n_lp,info)

        if (info.ne.0) write(*,*) 'DSYGV ERROR',info,i,natsc+1
        if (iproc.eq.0) then
        do iorb=1,norbi
        write(*,'(1x,a,i0,a,1x,1pe21.14)') 'evale(',iorb+iorbst-1,')=',evale(iorb)
        enddo
        endif
        do iorb=iorbst,min(norbi+iorbst-1,norb)
           eval(iorb)=evale(iorb-iorbst+1)
        enddo
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
     i_all=-product(shape(work_lp))*kind(work_lp)
     deallocate(work_lp,stat=i_stat)
     call memocc(i_stat,i_all,'work_lp','input_wf_diag')
     i_all=-product(shape(evale))*kind(evale)
     deallocate(evale,stat=i_stat)
     call memocc(i_stat,i_all,'evale','input_wf_diag')

     write(*,'(1x,a)',advance='no')'Building orthogonal Input Wavefunctions...'

     !allocate the wavefunction
     allocate(ppsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
     call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','input_wf_diag')

     !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
     iorbst=1
     imatrst=1
     norbsc=0
     do i=1,natsc
        norbi=norbsc_arr(i)
        norbsc=norbsc+norbi
        call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
             hamovr(imatrst,1),norbi,0.d0,ppsi(1,iorbst),nvctrp)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
     norbi=norbsc_arr(natsc+1)
     norbj=norb-norbsc
     call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
          hamovr(imatrst,1),norbi,0.d0,ppsi(1,iorbst),nvctrp)

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi','input_wf_diag')
     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr','input_wf_diag')
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr','input_wf_diag')

  endif

  if (iproc.eq.0) write(*,'(1x,a)')'done.'

  return
END SUBROUTINE

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
        allocate(hamks(norb,norb,2),stat=i_stat)
        call memocc(i_stat,product(shape(hamks))*kind(hamks),'hamks','kstrans_p')
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

        n_lp=max(4*norb,1000)
        allocate(work_lp(n_lp),stat=i_stat)
        call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','kstrans_p')
        call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
        evsum=0.d0
        do iorb=1,norb
        evsum=evsum+eval(iorb)*occup(iorb)
        if (iproc.eq.0) write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        enddo
        i_all=-product(shape(work_lp))*kind(work_lp)
        deallocate(work_lp,stat=i_stat)
        call memocc(i_stat,i_all,'work_lp','kstrans_p')
        if (info.ne.0) write(*,*) 'DSYEV ERROR',info

        allocate(psitt(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psitt))*kind(psitt),'psitt','kstrans_p')
! Transform to KS orbitals
      do iorb=1,norb
        call razero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        enddo
      enddo
        i_all=-product(shape(hamks))*kind(hamks)
        deallocate(hamks,stat=i_stat)
        call memocc(i_stat,i_all,'hamks','kstrans_p')

        call DCOPY(nvctrp*norbp*nproc,psitt,1,psit,1)
        i_all=-product(shape(psitt))*kind(psitt)
        deallocate(psitt,stat=i_stat)
        call memocc(i_stat,i_all,'psitt','kstrans_p')

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
        allocate(hamks(norb,norb,2),stat=i_stat)
        call memocc(i_stat,product(shape(hamks))*kind(hamks),'hamks','kstrans')
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

        n_lp=max(4*norb,1000)
        allocate(work_lp(n_lp),stat=i_stat)
        call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','kstrans')
        call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
        evsum=0.d0
        do iorb=1,norb
        evsum=evsum+eval(iorb)*occup(iorb)
        write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        enddo
        i_all=-product(shape(work_lp))*kind(work_lp)
        deallocate(work_lp,stat=i_stat)
        call memocc(i_stat,i_all,'work_lp','kstrans')
        if (info.ne.0) write(*,*) 'DSYEV ERROR',info

        allocate(psitt(nvctrp,norbp),stat=i_stat)
        call memocc(i_stat,product(shape(psitt))*kind(psitt),'psitt','kstrans')
! Transform to KS orbitals
      do iorb=1,norb
        call razero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psi(1,jorb),1,psitt(1,iorb),1)
        enddo
      enddo
        i_all=-product(shape(hamks))*kind(hamks)
        deallocate(hamks,stat=i_stat)
        call memocc(i_stat,i_all,'hamks','kstrans')

        call DCOPY(nvctrp*norbp,psitt,1,psi,1)
        i_all=-product(shape(psitt))*kind(psitt)
        deallocate(psitt,stat=i_stat)
        call memocc(i_stat,i_all,'psitt','kstrans')

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

   allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8),stat=i_stat)
   call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','reformatmywaves')

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

         allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2),stat=i_stat)
         call memocc(i_stat,product(shape(psigold))*kind(psigold),'psigold','reformatmywaves')

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

         i_all=-product(shape(psigold))*kind(psigold)
         deallocate(psigold,stat=i_stat)
         call memocc(i_stat,i_all,'psigold','reformatmywaves')
      end if
   end do
   
   i_all=-product(shape(psifscf))*kind(psifscf)
   deallocate(psifscf,stat=i_stat)
   call memocc(i_stat,i_all,'psifscf','reformatmywaves')

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

   allocate(psifscfold(-7:2*n1_old+8,-7:2*n2_old+8,-7:2*n3_old+8),stat=i_stat)
   call memocc(i_stat,product(shape(psifscfold))*kind(psifscfold),'psifscfold','reformatonewave')
   allocate(psifscfoex(-8:2*n1_old+9,-8:2*n2_old+9,-8:2*n3_old+9),stat=i_stat)
   call memocc(i_stat,product(shape(psifscfoex))*kind(psifscfoex),'psifscfoex','reformatonewave')

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

   i_all=-product(shape(psifscfold))*kind(psifscfold)
   deallocate(psifscfold,stat=i_stat)
   call memocc(i_stat,i_all,'psifscfold','reformatonewave')
   i_all=-product(shape(psifscfoex))*kind(psifscfoex)
   deallocate(psifscfoex,stat=i_stat)
   call memocc(i_stat,i_all,'psifscfoex','reformatonewave')
   allocate(psig(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
   call memocc(i_stat,product(shape(psig))*kind(psig),'psig','reformatonewave')
   allocate(ww((2*n1+16)*(2*n2+16)*(2*n3+16)),stat=i_stat)
   call memocc(i_stat,product(shape(ww))*kind(ww),'ww','reformatonewave')

        call analyse_shrink(n1,n2,n3,ww,psifscf,psig)
        call compress(n1,n2,n3,0,n1,0,n2,0,n3,  &
                    nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psig,psi(1),psi(nvctr_c+1))

   i_all=-product(shape(psig))*kind(psig)
   deallocate(psig,stat=i_stat)
   call memocc(i_stat,i_all,'psig','reformatonewave')
   i_all=-product(shape(ww))*kind(ww)
   deallocate(ww,stat=i_stat)
   call memocc(i_stat,i_all,'ww','reformatonewave')
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

   allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8),stat=i_stat)
   call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','readmywaves')

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

   i_all=-product(shape(psifscf))*kind(psifscf)
   deallocate(psifscf,stat=i_stat)
   call memocc(i_stat,i_all,'psifscf','readmywaves')

   call cpu_time(tr1)
   call system_clock(ncount2,ncount_rate,ncount_max)
   tel=dble(ncount2-ncount1)/dble(ncount_rate)
   write(*,'(a,i4,2(1x,e10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
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
   write(*,'(1x,i2,6(1x,e14.7))') iproc,(center(j),j=1,3),(center_old(j),j=1,3)
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
      
      allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2),stat=i_stat)
      call memocc(i_stat,product(shape(psigold))*kind(psigold),'psigold','readonewave')

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

      i_all=-product(shape(psigold))*kind(psigold)
      deallocate(psigold,stat=i_stat)
      call memocc(i_stat,i_all,'psigold','readonewave')

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
       write(*,'(a,i4,2(1x,e10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel


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

    write(*,'(1x,i0,a)') iorb,'th wavefunction written'


    END SUBROUTINE
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,rx,ry,rz,psi)
        implicit real*8 (a-h,o-z)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f)
        real*8, allocatable :: psifscf(:),psir(:),psig(:,:,:,:,:,:)

        allocate(psig(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
        call memocc(i_stat,product(shape(psig))*kind(psig),'psig','plot_wf')
        allocate(psifscf((2*n1+31)*(2*n2+31)*(2*n3+16)),stat=i_stat)
        call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','plot_wf')
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
        call memocc(i_stat,product(shape(psir))*kind(psir),'psir','plot_wf')

        call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
                    nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1),psi(nvctr_c+1),psig)
        call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n1+2))
        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

        call plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,psir)

        i_all=-product(shape(psig))*kind(psig)
        deallocate(psig,stat=i_stat)
        call memocc(i_stat,i_all,'psig','plot_wf')
        i_all=-product(shape(psifscf))*kind(psifscf)
        deallocate(psifscf,stat=i_stat)
        call memocc(i_stat,i_all,'psifscf','plot_wf')
        i_all=-product(shape(psir))*kind(psir)
        deallocate(psir,stat=i_stat)
        call memocc(i_stat,i_all,'psir','plot_wf')
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

        allocate(ipiv(idsx+1),stat=i_stat)
        call memocc(i_stat,product(shape(ipiv))*kind(ipiv),'ipiv','diisstp')
        allocate(rds(idsx+1),stat=i_stat)
        call memocc(i_stat,product(shape(rds))*kind(rds),'rds','diisstp')

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
           !write(*,*) 'DIIS weights'
           write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j),j=1,min(idsx,ids)+1)
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

        i_all=-product(shape(ipiv))*kind(ipiv)
        deallocate(ipiv,stat=i_stat)
        call memocc(i_stat,i_all,'ipiv','diisstp')
        i_all=-product(shape(rds))*kind(rds)
        deallocate(rds,stat=i_stat)
        call memocc(i_stat,i_all,'rds','diisstp')
        call timing(iproc,'Diis          ','OF')

        return
        END SUBROUTINE

subroutine iguess_generator(iproc,izatom,ielpsp,psppar,npspcode,ng,nl,nmax_occ,occupat,expo,psiat)
  implicit none
  integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,nmax_occ
  real(kind=8), dimension(0:4,0:4), intent(in) :: psppar
  integer, dimension(4), intent(out) :: nl
  real(kind=8), dimension(ng+1), intent(out) :: expo
  real(kind=8), dimension(nmax_occ), intent(out) :: occupat
  real(kind=8), dimension(ng+1,nmax_occ), intent(out) :: psiat
  
  !local variables
  character(len=27) :: string 
  character(len=2) :: symbol
  integer, parameter :: lmax=3,nint=100,noccmax=2
  real(kind=8), parameter :: fact=4.d0
  real(kind=8), dimension(:), allocatable :: xp,gpot,alps,ott
  real(kind=8), dimension(:,:), allocatable :: aeval,chrg,res,vh,hsep,occup,ofdcoef
  real(kind=8), dimension(:,:,:), allocatable :: psi
  real(kind=8), dimension(:,:,:,:), allocatable :: rmt
  integer, dimension(:,:), allocatable :: neleconf
  logical :: exists
  integer :: lpx,ncount,nsccode
  integer :: l,i,j,iocc,il,lwrite,i_all,i_stat
  real(kind=8) :: alpz,alpl,rcov,rprb,zion,rij,a,a0,a0in,tt,ehomo

  !filename = 'psppar.'//trim(atomname)

  lpx=0
  lpx_determination: do i=1,4
     if (psppar(i,0) == 0.d0) then
     exit lpx_determination
     else
        lpx=i-1
     end if
  end do lpx_determination

  allocate(gpot(3),stat=i_stat)
  call memocc(i_stat,product(shape(gpot))*kind(gpot),'gpot','iguess_generator')
  allocate(alps(lpx+1),stat=i_stat)
  call memocc(i_stat,product(shape(alps))*kind(alps),'alps','iguess_generator')
  allocate(hsep(6,lpx+1),stat=i_stat)
  call memocc(i_stat,product(shape(hsep))*kind(hsep),'hsep','iguess_generator')
  allocate(ott(6),stat=i_stat)
  call memocc(i_stat,product(shape(ott))*kind(ott),'ott','iguess_generator')
  allocate(occup(noccmax,lmax+1),stat=i_stat)
  call memocc(i_stat,product(shape(occup))*kind(occup),'occup','iguess_generator')
  allocate(ofdcoef(3,4),stat=i_stat)
  call memocc(i_stat,product(shape(ofdcoef))*kind(ofdcoef),'ofdcoef','iguess_generator')
  allocate(neleconf(6,4),stat=i_stat)
  call memocc(i_stat,product(shape(neleconf))*kind(neleconf),'neleconf','iguess_generator')

  !assignation of radii and coefficients of the local part
  alpz=psppar(0,0)
  alpl=psppar(0,0)
  alps(1:lpx+1)=psppar(1:lpx+1,0)
  gpot(1:3)=psppar(0,1:3)

  !assignation of the coefficents for the nondiagonal terms
  if (npspcode == 2) then !GTH case
     ofdcoef(:,:)=0.d0
  else if (npspcode == 3) then !HGH case
     ofdcoef(1,1)=-0.5d0*sqrt(3.d0/5.d0) !h2
     ofdcoef(2,1)=0.5d0*sqrt(5.d0/21.d0) !h4
     ofdcoef(3,1)=-0.5d0*sqrt(100.d0/63.d0) !h5

     ofdcoef(1,2)=-0.5d0*sqrt(5.d0/7.d0) !h2
     ofdcoef(2,2)=1.d0/6.d0*sqrt(35.d0/11.d0) !h4
     ofdcoef(3,2)=-7.d0/3.d0*sqrt(1.d0/11.d0) !h5

     ofdcoef(1,3)=-0.5d0*sqrt(7.d0/9.d0) !h2
     ofdcoef(2,3)=0.5d0*sqrt(63.d0/143.d0) !h4
     ofdcoef(3,3)=-9.d0*sqrt(1.d0/143.d0) !h5

     ofdcoef(1,4)=0.d0 !h2
     ofdcoef(2,4)=0.d0 !h4
     ofdcoef(3,4)=0.d0 !h5
  end if
  !define the values of hsep starting from the pseudopotential file
  do l=1,lpx+1
     hsep(1,l)=psppar(l,1)
     hsep(2,l)=psppar(l,2)*ofdcoef(1,l)
     hsep(3,l)=psppar(l,2)
     hsep(4,l)=psppar(l,3)*ofdcoef(2,l)
     hsep(5,l)=psppar(l,3)*ofdcoef(3,l)
     hsep(6,l)=psppar(l,3)
  end do

  !Now the treatment of the occupation number
  call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode)

  occup(:,:)=0.d0
   do l=0,lmax
     iocc=0
     do i=1,6
        ott(i)=real(neleconf(i,l+1),kind=8)
        if (ott(i).gt.0.d0) then
           iocc=iocc+1
            if (iocc.gt.noccmax) stop 'iguess_generator: noccmax too small'
           occup(iocc,l+1)=ott(i)
        endif
     end do
     nl(l+1)=iocc
  end do

  !allocate arrays for the gatom routine
  allocate(aeval(noccmax,lmax+1),stat=i_stat)
  call memocc(i_stat,product(shape(aeval))*kind(aeval),'aeval','iguess_generator')
  allocate(chrg(noccmax,lmax+1),stat=i_stat)
  call memocc(i_stat,product(shape(chrg))*kind(chrg),'chrg','iguess_generator')
  allocate(res(noccmax,lmax+1),stat=i_stat)
  call memocc(i_stat,product(shape(res))*kind(res),'res','iguess_generator')
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2),stat=i_stat)
  call memocc(i_stat,product(shape(vh))*kind(vh),'vh','iguess_generator')
  allocate(psi(0:ng,noccmax,lmax+1),stat=i_stat)
  call memocc(i_stat,product(shape(psi))*kind(psi),'psi','iguess_generator')
  allocate(xp(0:ng),stat=i_stat)
  call memocc(i_stat,product(shape(xp))*kind(xp),'xp','iguess_generator')
  allocate(rmt(nint,0:ng,0:ng,lmax+1),stat=i_stat)
  call memocc(i_stat,product(shape(rmt))*kind(rmt),'rmt','iguess_generator')

  zion=real(ielpsp,kind=8)

  !can be switched on for debugging
  !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
  !     'Input Guess Generation for atom',trim(atomname),&
  !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb

  rij=3.d0
  ! exponents of gaussians
  a0in=alpz
  a0=a0in/rij
  !       tt=sqrt(sqrt(2.d0))
  tt=2.d0**.3d0
  do i=0,ng
     a=a0*tt**i
     xp(i)=.5d0/a**2
  end do

  ! initial guess
  do l=0,lmax
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.d0
        end do
     end do
  end do

  call crtvh(ng,lmax,xp,vh,rprb,fact,nint,rmt)

  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nint,&
       aeval,ng,psi,res,chrg)

  !postreatment of the inguess data
  do i=1,ng+1
     expo(i)=sqrt(0.5/xp(i-1))
  end do

  i=0
  do l=1,4
     do iocc=1,nl(l)
        i=i+1
        occupat(i)=occup(iocc,l)
        do j=1,ng+1
           psiat(j,i)=psi(j-1,iocc,l)
        end do
     end do
  end do

  i_all=-product(shape(aeval))*kind(aeval)
  deallocate(aeval,stat=i_stat)
  call memocc(i_stat,i_all,'aeval','iguess_generator')
  i_all=-product(shape(chrg))*kind(chrg)
  deallocate(chrg,stat=i_stat)
  call memocc(i_stat,i_all,'chrg','iguess_generator')
  i_all=-product(shape(res))*kind(res)
  deallocate(res,stat=i_stat)
  call memocc(i_stat,i_all,'res','iguess_generator')
  i_all=-product(shape(vh))*kind(vh)
  deallocate(vh,stat=i_stat)
  call memocc(i_stat,i_all,'vh','iguess_generator')
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi','iguess_generator')
  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp','iguess_generator')
  i_all=-product(shape(rmt))*kind(rmt)
  deallocate(rmt,stat=i_stat)
  call memocc(i_stat,i_all,'rmt','iguess_generator')
  i_all=-product(shape(gpot))*kind(gpot)
  deallocate(gpot,stat=i_stat)
  call memocc(i_stat,i_all,'gpot','iguess_generator')
  i_all=-product(shape(hsep))*kind(hsep)
  deallocate(hsep,stat=i_stat)
  call memocc(i_stat,i_all,'hsep','iguess_generator')
  i_all=-product(shape(alps))*kind(alps)
  deallocate(alps,stat=i_stat)
  call memocc(i_stat,i_all,'alps','iguess_generator')
  i_all=-product(shape(ott))*kind(ott)
  deallocate(ott,stat=i_stat)
  call memocc(i_stat,i_all,'ott','iguess_generator')
  i_all=-product(shape(occup))*kind(occup)
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup','iguess_generator')
  i_all=-product(shape(ofdcoef))*kind(ofdcoef)
  deallocate(ofdcoef,stat=i_stat)
  call memocc(i_stat,i_all,'ofdcoef','iguess_generator')
  i_all=-product(shape(neleconf))*kind(neleconf)
  deallocate(neleconf,stat=i_stat)
  call memocc(i_stat,i_all,'neleconf','iguess_generator')

end subroutine iguess_generator

subroutine gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
     zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
     aeval,ng,psi,res,chrg)
     implicit real*8 (a-h,o-z)
     logical noproj
     parameter(nint=100)
     dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
       aux(2*ng+2),&
       gpot(3),hsep(6,lpx+1),rmt(nint,0:ng,0:ng,lmax+1),&
       pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
       potgrd(nint),&
       rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(nint),&
       occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
       vh(0:ng,0:ng,4,0:ng,0:ng,4),&
       res(noccmax,lmax+1),xp(0:ng)
     if (nintp.ne.nint) stop 'nint><nintp'

     do l=0,lmax
        if (occup(1,l+1).gt.0.d0) lcx=l
     end do
     !write(6,*) 'lcx',lcx
 
     noproj=.true.
     do l=1,lpx+1
        noproj = noproj .and. (alps(l) .eq. 0.d0)
     end do


!    projectors, just in case
     if (.not. noproj) then
        do l=0,lpx
           gml1=sqrt( gamma(l+1.5d0) / (2.d0*alps(l+1)**(2*l+3)) )
           gml2=sqrt( gamma(l+3.5d0) / (2.d0*alps(l+1)**(2*l+7)) )&
               /(l+2.5d0)
           gml3=sqrt( gamma(l+5.5d0) / (2.d0*alps(l+1)**(2*l+11)) )&
               /((l+3.5d0)*(l+4.5d0))
           tt=1.d0/(2.d0*alps(l+1)**2)
           do i=0,ng
              ttt=1.d0/(xp(i)+tt)
              pp1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
              pp2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
              pp3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
           end do
        end do
     else
        pp1(:,:)=0.d0
        pp2(:,:)=0.d0
        pp3(:,:)=0.d0
     end if

     do l=0,lmax
        do j=0,ng
           do i=0,ng
              rho(i,j,l+1)=0.d0
           end do
        end do
     end do

     evsum=1.d30
     do 2000,it=1,50
        evsumold=evsum
        evsum=0.d0
        
        ! coefficients of charge density
        do l=0,lmax
           do j=0,ng
              do i=0,ng
                 rhoold(i,j,l+1)=rho(i,j,l+1)
                 rho(i,j,l+1)=0.d0        
              end do
           end do
        end do

        do l=0,lmax
           do iocc=1,noccmax
              if (occup(iocc,l+1).gt.0.d0) then
                 do j=0,ng
                    do i=0,ng
                       rho(i,j,l+1)=rho(i,j,l+1) + &
                            psi(i,iocc,l+1)*psi(j,iocc,l+1)*occup(iocc,l+1)
                    end do
                 end do
              end if
           end do
        end do


        rmix=.5d0
        if (it.eq.1) rmix=1.d0
        do 2834,l=0,lmax
        do 2834,j=0,ng
        do 2834,i=0,ng
        tt=rmix*rho(i,j,l+1) + (1.d0-rmix)*rhoold(i,j,l+1)
        rho(i,j,l+1)=tt
2834        continue

! XC potential on grid
!        do 3728,k=1,nint
!3728        xcgrd(k)=0.d0        
!        do 3328,l=0,lmax
!        do 3328,j=0,ng
!        do 3328,i=0,ng
!        do 3328,k=1,nint
!3328        xcgrd(k)=xcgrd(k)+rmt(k,i,j,l+1)*rho(i,j,l+1)
        call DGEMV('N',nint,(lcx+1)*(ng+1)**2,1.d0,&
                   rmt,nint,rho,1,0.d0,xcgrd,1)

        dr=fact*rprb/nint
        do 3167,k=1,nint
        r=(k-.5d0)*dr
! divide by 4 pi
        tt=xcgrd(k)*0.07957747154594768d0
! multiply with r^2 to speed up calculation of matrix elements
        xcgrd(k)=emuxc(tt)*r**2
3167        continue

        do 1000,l=0,lmax
        gml=.5d0*gamma(.5d0+l)

!  lower triangles only
        do 100,i=0,ng
        do 100,j=0,i
        d=xp(i)+xp(j)
        sxp=1.d0/d
        const=gml*sqrt(sxp)**(2*l+1)
! overlap        
        ss(i,j)=const*sxp*(l+.5d0)
! kinetic energy
        hh(i,j)=.5d0*const*sxp**2* ( 3.d0*xp(i)*xp(j) +&
             l*(6.d0*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -&
             l**2*(xp(i)-xp(j))**2  ) + .5d0*l*(l+1.d0)*const
! potential energy from parabolic potential
        hh(i,j)=hh(i,j) +&
             .5d0*const*sxp**2*(l+.5d0)*(l+1.5d0)/rprb**4 
! hartree potential from ionic core charge
        tt=sqrt(1.d0+2.d0*alpz**2*d)
        if (l.eq.0) then
        hh(i,j)=hh(i,j) -zion/(2.d0*d*tt)
        else if (l.eq.1) then
        hh(i,j)=hh(i,j) -zion* &
             (1.d0 + 3.d0*alpz**2*d)/(2.d0*d**2*tt**3)
        else if (l.eq.2) then
        hh(i,j)=hh(i,j) -zion* &
             (2.d0 + 10.d0*alpz**2*d + 15.d0*alpz**4*d**2)/(2.d0*d**3*tt**5)
        else if (l.eq.3) then
        hh(i,j)=hh(i,j) -zion*3.d0* &
             (2.d0 +14.d0*alpz**2*d +35.d0*alpz**4*d**2 +35.d0*alpz**6*d**3)/&
             (2.d0*d**4*tt**7)
        else 
        stop 'l too big'
        end if
! potential from repulsive gauss potential
        tt=alpl**2/(.5d0+d*alpl**2)
        hh(i,j)=hh(i,j)  + gpot(1)*.5d0*gamma(1.5d0+l)*tt**(1.5d0+l)&
             + (gpot(2)/alpl**2)*.5d0*gamma(2.5d0+l)*tt**(2.5d0+l)&
             + (gpot(3)/alpl**4)*.5d0*gamma(3.5d0+l)*tt**(3.5d0+l)
        ! separable terms
        if (l.le.lpx) then
           hh(i,j)=hh(i,j) + pp1(i,l+1)*hsep(1,l+1)*pp1(j,l+1)&
                + pp1(i,l+1)*hsep(2,l+1)*pp2(j,l+1)&
                + pp2(i,l+1)*hsep(2,l+1)*pp1(j,l+1)&
                + pp2(i,l+1)*hsep(3,l+1)*pp2(j,l+1)&
                + pp1(i,l+1)*hsep(4,l+1)*pp3(j,l+1)&
                + pp3(i,l+1)*hsep(4,l+1)*pp1(j,l+1)&
                + pp2(i,l+1)*hsep(5,l+1)*pp3(j,l+1)&
                + pp3(i,l+1)*hsep(5,l+1)*pp2(j,l+1)&
                + pp3(i,l+1)*hsep(6,l+1)*pp3(j,l+1)
        end if
! hartree potential from valence charge distribution
!        tt=0.d0
!        do 4982,lp=0,lcx
!        do 4982,jp=0,ng
!        do 4982,ip=0,ng
!        tt=tt + vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)
!4982        continue
        tt=DDOT((lcx+1)*(ng+1)**2,vh(0,0,1,i,j,l+1),1,rho(0,0,1),1)
        hh(i,j)=hh(i,j) + tt
! potential from XC potential
        dr=fact*rprb/nint
!        tt=0.d0
!        do 8049,k=1,nint
!8049        tt=tt+xcgrd(k)*rmt(k,i,j,l+1)
        tt=DDOT(nint,rmt(1,i,j,l+1),1,xcgrd(1),1)
        hh(i,j)=hh(i,j)+tt*dr
100        continue
 
! ESSL
!        call DSYGV(1,hh,ng+1,ss,ng+1,eval,evec,ng+1,ng+1,aux,2*ng+2)
! LAPACK
        call DSYGV(1,'V','L',&
             ng+1,hh,ng+1,ss,ng+1,eval,evec,(ng+1)**2,info)
        if (info.ne.0) write(6,*) 'LAPACK',info
        do 334,iocc=0,noccmax-1
        do 334,i=0,ng
334        evec(i,iocc)=hh(i,iocc)
! end LAPACK
        do 9134,iocc=1,noccmax
        evsum=evsum+eval(iocc-1)
        aeval(iocc,l+1)=eval(iocc-1)
        do 9134,i=0,ng
9134        psi(i,iocc,l+1)=evec(i,iocc-1)

!        write(6,*) 'eval',l
!55        format(5(e14.7))
!        write(6,55) eval 
!        write(6,*) 'diff eval'
!        write(6,55) (eval(i)-eval(i-1),i=1,ng)
        
!        write(6,*) 'evec',l
!        do i=0,ng
!33        format(10(e9.2))
!        write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
!        end do

1000        continue

        tt=abs(evsum-evsumold)
!        write(6,*) 'evdiff',it,tt
        if (tt.lt.1.d-12) goto 3000
2000        continue
3000        continue
        call resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,ng,res,&
             zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,nint,&
             potgrd,xcgrd)

! charge up to radius rcov
        if (lmax.gt.3) stop 'cannot calculate chrg'
        do 3754,l=0,lmax
        do 3754,iocc=1,noccmax
3754        chrg(iocc,l+1)=0.d0

        do 3761,iocc=1,noccmax
        do 3761,j=0,ng
        do 3761,i=0,ng
        d=xp(i)+xp(j)
        sd=sqrt(d)
        terf=erf(sd*rcov) 
        texp=exp(-d*rcov**2)

        tt=0.4431134627263791d0*terf/sd**3 - 0.5d0*rcov*texp/d
        chrg(iocc,1)=chrg(iocc,1) + psi(i,iocc,1)*psi(j,iocc,1)*tt
        if (lmax.eq.0) goto 3761
        tt=0.6646701940895686d0*terf/sd**5 + &
             (-0.75d0*rcov*texp - 0.5d0*d*rcov**3*texp)/d**2
        chrg(iocc,2)=chrg(iocc,2) + psi(i,iocc,2)*psi(j,iocc,2)*tt
        if (lmax.eq.1) goto 3761
        tt=1.661675485223921d0*terf/sd**7 + &
             (-1.875d0*rcov*texp-1.25d0*d*rcov**3*texp-.5d0*d**2*rcov**5*texp) &
             /d**3
        chrg(iocc,3)=chrg(iocc,3) + psi(i,iocc,3)*psi(j,iocc,3)*tt
        if (lmax.eq.2) goto 3761
        tt=5.815864198283725d0*terf/sd**9 + &
             (-6.5625d0*rcov*texp - 4.375d0*d*rcov**3*texp - &
             1.75d0*d**2*rcov**5*texp - .5d0*d**3*rcov**7*texp)/d**4
        chrg(iocc,4)=chrg(iocc,4) + psi(i,iocc,4)*psi(j,iocc,4)*tt
3761        continue


        !writing lines suppressed
!!$        write(66,*)  lmax+1
!!$        write(66,*) ' #LINETYPE{1324}' 
!!$        write(66,*) ' $' 
        do l=0,lmax
!!$           write(66,*) ' 161'
           r=0.d0
           do
              tt= wave(ng,l,xp,psi(0,1,l+1),r)
!!$              write(66,*) r,tt
              r=r+.025d0
              if(r > 4.00001d0) exit
           end do
        end do
        !writing lines suppressed
!!$        write(67,*) min(lmax+1,3)
!!$        write(67,*) ' #LINETYPE{132}'
!!$        write(67,*) ' #TITLE{FOURIER}' 
!!$        write(67,*) ' $'
        dr=6.28d0/rprb/200.d0
!!$        write(67,*) ' 200'
        rk=0.d0
        loop_rk1: do 
           tt=0.d0
           do i=0,ng
              texp=exp(-.25d0*rk**2/xp(i))
!              texp=exp(-.5d0*energy/xp(i))
              sd=sqrt(xp(i))
              tt=tt+psi(i,1,1)*0.4431134627263791d0*texp/sd**3
           end do
!!$           write(67,*) rk,tt
           rk=rk+dr
           if(rk > 6.28d0/rprb-.5d0*dr) exit loop_rk1
        end do loop_rk1
        if (lmax.ge.1) then
!!$           write(67,*) ' 200'
           rk=0.d0
           loop_rk2: do 
              tt=0.d0
              do i=0,ng
                 texp=exp(-.25d0*rk**2/xp(i))
                 sd=sqrt(xp(i))
                 tt=tt+psi(i,1,2)*0.2215567313631895d0*rk*texp/sd**5
              end do
!!$              write(67,*) rk,tt
              rk=rk+dr
              if (rk > 6.28d0/rprb-.5d0*dr) exit loop_rk2
           end do loop_rk2
        end if
        if (lmax.ge.2) then
!!$           write(67,*) ' 200'
           rk=0.d0
           do 
              tt=0.d0
              do i=0,ng
                 texp=exp(-.25d0*rk**2/xp(i))
                 sd=sqrt(xp(i))
              tt=tt+psi(i,1,3)*0.1107783656815948d0*rk**2*texp/sd**7
              end do
!!$              write(67,*) rk,tt
              rk=rk+dr
              if (rk > 6.28d0/rprb-.5d0*dr) exit
           end do
        end if


      end subroutine gatom



      subroutine resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,&
           ng,res,zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,nint,&
           potgrd,xcgrd)
        implicit real*8 (a-h,o-z)
        dimension psi(0:ng,noccmax,lmax+1),rho(0:ng,0:ng,lmax+1),&
             gpot(3),pp1(0:ng,lmax+1),pp2(0:ng,lmax+1),pp3(0:ng,lmax+1),&
             alps(lmax+1),hsep(6,lmax+1),res(noccmax,lmax+1),xp(0:ng),&
             xcgrd(nint),aeval(noccmax,lmax+1),potgrd(nint)
        
!   potential on grid 
        dr=fact*rprb/nint
        do 9873,k=1,nint
        r=(k-.5d0)*dr
        potgrd(k)= .5d0*(r/rprb**2)**2 - &
             zion*erf(r/(sqrt(2.d0)*alpz))/r &
             + exp(-.5d0*(r/alpl)**2)*&
             ( gpot(1) + gpot(2)*(r/alpl)**2 + gpot(3)*(r/alpl)**4 )&
             + xcgrd(k)/r**2
        do 2487,j=0,ng
        do 2487,i=0,ng
        spi=1.772453850905516d0
        d=xp(i)+xp(j)
        sd=sqrt(d)
        tx=exp(-d*r**2)
        tt=spi*erf(sd*r)
           ud0=tt/(4.d0*sd**3*r)
        potgrd(k)=potgrd(k)+ud0*rho(i,j,1)
           ud1=-tx/(4.d0*d**2) + 3.d0*tt/(8.d0*sd**5*r)
        if (lmax.ge.1) potgrd(k)=potgrd(k)+ud1*rho(i,j,2)
        ud2=-tx*(7.d0 + 2.d0*d*r**2)/(8.d0*d**3) +&
             15.d0*tt/(16.d0*sd**7*r)
        if (lmax.ge.2) potgrd(k)=potgrd(k)+ud2*rho(i,j,3)
        ud3=-tx*(57.d0+22.d0*d*r**2+4.d0*d**2*r**4)/(16.d0*d**4) + &
             105.d0*tt/(32.d0*sd**9*r)
        if (lmax.ge.3) potgrd(k)=potgrd(k)+ud3*rho(i,j,4)
2487    continue
9873        continue

        do 1500,ll=0,lmax
        if (ll.le.lpx) then
        rnrm1=1.d0/sqrt(.5d0*gamma(ll+1.5d0)*alps(ll+1)**(2*ll+3))
        rnrm2=1.d0/sqrt(.5d0*gamma(ll+3.5d0)*alps(ll+1)**(2*ll+7))
        rnrm3=1.d0/sqrt(.5d0*gamma(ll+5.5d0)*alps(ll+1)**(2*ll+11))
        end if
        do 1500,iocc=1,noccmax
! separabel part
        if (ll.le.lpx) then
        scpr1=DDOT(ng+1,psi(0,iocc,ll+1),1,pp1(0,ll+1),1)
        scpr2=DDOT(ng+1,psi(0,iocc,ll+1),1,pp2(0,ll+1),1)
        scpr3=DDOT(ng+1,psi(0,iocc,ll+1),1,pp3(0,ll+1),1)
        end if
        res(iocc,ll+1)=0.d0
        do 1500,j=1,nint
!  wavefunction on grid
        r=(j-.5d0)*dr
        psigrd = wave(ng,ll,xp,psi(0,iocc,ll+1),r)
!   kinetic energy        
        rkin=0.d0
        do 5733,i=0,ng
        rkin=rkin + psi(i,iocc,ll+1) *  (&
             xp(i)*(3.d0+2.d0*ll-2.d0*xp(i)*r**2)*exp(-xp(i)*r**2) )
5733        continue
        rkin=rkin*r**ll
!   separabel part
        if (ll.le.lpx) then
           sep =& 
                (scpr1*hsep(1,ll+1) + scpr2*hsep(2,ll+1) + scpr3*hsep(4,ll+1))&
                *rnrm1*r**ll*exp(-.5d0*(r/alps(ll+1))**2)   +&
                (scpr1*hsep(2,ll+1) + scpr2*hsep(3,ll+1) + scpr3*hsep(5,ll+1))&
                *rnrm2*r**(ll+2)*exp(-.5d0*(r/alps(ll+1))**2)   +&
                (scpr1*hsep(4,ll+1) + scpr2*hsep(5,ll+1) + scpr3*hsep(6,ll+1))&
                *rnrm3*r**(ll+4)*exp(-.5d0*(r/alps(ll+1))**2)
        else
        sep=0.d0
        end if
! resdidue
        tt=rkin+sep+(potgrd(j)-aeval(iocc,ll+1))*psigrd
!384        format(6(e12.5))
!12        format(i2,i2,e9.2,3(e12.5),e10.3)
        res(iocc,ll+1)=res(iocc,ll+1) + tt**2*dr
1500        continue
!        do 867,l=0,lmax
!        do 867,iocc=1,noccmax
!867        write(6,*) 'res',l,iocc,res(iocc,l+1)
        return
        end subroutine resid



        subroutine crtvh(ng,lmax,xp,vh,rprb,fact,nint,rmt)
        implicit real*8 (a-h,o-z)
        dimension vh(0:ng,0:ng,0:3,0:ng,0:ng,0:3),xp(0:ng),&
             rmt(nint,0:ng,0:ng,lmax+1)
        if (lmax.gt.3) stop 'crtvh'

        dr=fact*rprb/nint
        do 8049,l=0,lmax
        do 8049,k=1,nint
        r=(k-.5d0)*dr
        do 8049,j=0,ng
        do 8049,i=0,ng
8049        rmt(k,i,j,l+1)=(r**2)**l*exp(-(xp(i)+xp(j))*r**2)

        do 100,j=0,ng
        do 100,i=0,ng
        c=xp(i)+xp(j)
        do 100,jp=0,ng
        do 100,ip=0,ng
        d=xp(ip)+xp(jp)
        scpd=sqrt(c+d)
        vh(ip,jp,0,i,j,0)=0.2215567313631895d0/(c*d*scpd)
        vh(ip,jp,1,i,j,0)=&
             .1107783656815948d0*(2.d0*c+3.d0*d)/(c*d**2*scpd**3)
        vh(ip,jp,2,i,j,0)=.05538918284079739d0*&
             (8.d0*c**2+20.d0*c*d+15.d0*d**2)/(c*d**3*scpd**5)
        vh(ip,jp,3,i,j,0)=.0830837742611961d0*&
        (16.d0*c**3+56.d0*c**2*d+70.d0*c*d**2+35.d0*d**3)/&
             (c*d**4*scpd**7)


        vh(ip,jp,0,i,j,1)=&
             .1107783656815948d0*(3.d0*c+2.d0*d)/(c**2*d*scpd**3)
        vh(ip,jp,1,i,j,1)=&
             .05538918284079739d0*(6.d0*c**2+15.d0*c*d+6.d0*d**2)/&
             (c**2*d**2*scpd**5)
        vh(ip,jp,2,i,j,1)=.02769459142039869d0*&
             (24.d0*c**3+84.d0*c**2*d+105.d0*c*d**2+30.d0*d**3)/&
             (c**2*d**3*scpd**7)
        vh(ip,jp,3,i,j,1)=0.04154188713059803d0*&
             (48.d0*c**4+216.d0*c**3*d+378.d0*c**2*d**2+&
             315.d0*c*d**3+70.d0*d**4)/(c**2*d**4*scpd**9)

        vh(ip,jp,0,i,j,2)=&
             .05538918284079739d0*(15.d0*c**2+20.d0*c*d+8.d0*d**2)/&
             (c**3*d*scpd**5)
        vh(ip,jp,1,i,j,2)=.02769459142039869d0*&
             (30.d0*c**3+105.d0*c**2*d+84.d0*c*d**2+24.d0*d**3)/&
             (c**3*d**2*scpd**7)
        vh(ip,jp,2,i,j,2)=&
             .2077094356529901d0*(8.d0*c**4+36.d0*c**3*d+63.d0*c**2*d**2+&
             36.d0*c*d**3+8.d0*d**4)/(c**3*d**3*scpd**9)
        vh(ip,jp,3,i,j,2)=&
             .1038547178264951d0*(48.d0*c**5+264.d0*c**4*d+594.d0*c**3*d**2+&
             693.d0*c**2*d**3+308.d0*c*d**4+56.d0*d**5)/&
             (c**3*d**4*scpd**11)

        vh(ip,jp,0,i,j,3)=.0830837742611961d0*&
             (35.d0*c**3+70.d0*c**2*d+56.d0*c*d**2+16.d0*d**3)/&
             (c**4*d*scpd**7)
        vh(ip,jp,1,i,j,3)=&
             .04154188713059803d0*(70.d0*c**4+315.d0*c**3*d+378.d0*c**2*d**2+&
             216.d0*c*d**3+48.d0*d**4)/(c**4*d**2*scpd**9)
        vh(ip,jp,2,i,j,3)=&
             .1038547178264951d0*(56.d0*c**5+308.d0*c**4*d+693.d0*c**3*d**2+&
             594.d0*c**2*d**3+264.d0*c*d**4+48.d0*d**5)/&
             (c**4*d**3*scpd**11)
        vh(ip,jp,3,i,j,3)=&
             1.090474537178198d0*(16.d0*c**6+104.d0*c**5*d+286.d0*c**4*d**2+&
             429.d0*c**3*d**3+286.d0*c**2*d**4+104.d0*c*d**5+16.d0*d**6)/&
             (c**4*d**4*scpd**13)
100        continue
        return
        end subroutine crtvh

        real*8 function wave(ng,ll,xp,psi,r)
        implicit real*8 (a-h,o-z)
        dimension psi(0:ng),xp(0:ng)

        wave=0.d0
        do 9373,i=0,ng
9373    wave=wave + psi(i)*exp(-xp(i)*r**2)
        wave=wave*r**ll
        return
        end function wave


      real*8 function emuxc(rho)
      implicit real*8 (a-h,o-z)
      parameter (a0p=.4581652932831429d0,&
           a1p=2.217058676663745d0,&
           a2p=0.7405551735357053d0,&
           a3p=0.01968227878617998d0)
      parameter (b1p=1.0d0,&
           b2p=4.504130959426697d0,&
           b3p=1.110667363742916d0,&
           b4p=0.02359291751427506d0)
      parameter (rsfac=.6203504908994000d0,ot=1.d0/3.d0)
      parameter (c1=4.d0*a0p*b1p/3.0d0,  c2=5.0d0*a0p*b2p/3.0d0+a1p*b1p,&
           c3=2.0d0*a0p*b3p+4.0d0*a1p*b2p/3.0d0+2.0d0*a2p*b1p/3.0d0,&
           c4=7.0d0*a0p*b4p/3.0d0+5.0d0*a1p*b3p/3.0d0+a2p*b2p+a3p*b1p/3.0d0,&
           c5=2.0d0*a1p*b4p+4.0d0*a2p*b3p/3.0d0+2.0d0*a3p*b2p/3.0d0,&
           c6=5.0d0*a2p*b4p/3.0d0+a3p*b3p,c7=4.0d0*a3p*b4p/3.0d0)
      if(rho.lt.1.d-24) then
        emuxc=0.d0
      else
        if(rho.lt.0.d0) write(6,*) ' rho less than zero',rho
        rs=rsfac*rho**(-ot)
        top=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
        bot=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
        emuxc=top/(bot*bot)
      end if
      end function emuxc

      real*8 function gamma(x)
!     restricted version of the Gamma function
      implicit real*8 (a-h,o-z)

      if (x.le.0.d0) stop 'wrong argument for gamma'
      if (mod(x,1.d0).eq.0.d0) then
         ii=x
         do i=2,ii
            gamma=gamma*(i-1)
         end do
      else if (mod(x,.5d0).eq.0.d0) then
         ii=x-.5d0
!         gamma=sqrt(3.14159265358979d0)
         gamma=1.772453850905516027d0
         do i=1,ii
            gamma=gamma*(i-.5d0)
         end do
      else
         stop 'wrong argument for gamma'
      end if
      end function gamma


end module
