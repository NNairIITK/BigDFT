module libBigDFT
  
  private

  !- Do a minimisation loop to find forces and energy.
  public :: cluster

  !input variables
  public :: read_input_variables
  
  !- Initialisation methods.
  !- Create and allocate access arrays for wavefunctions.
  public :: createWavefunctionsDescriptors
  !- Create and allocate data descriptors for nonlocal PSP projectors
  public :: createProjectorsArrays
  !- Create the descriptors for the density and the potential
  public :: createDensPotDescriptors
  !- Import CP2K files
  public :: import_gaussians
  !- Compute input guess wavefunctions from atomic orbitals.
  public :: input_wf_diag
  
  !- Orthogonalise after restart (reformat or read)
  public :: first_orthon
  !- KS transformation after minimisation
  public :: last_orthon
  !- SCF handling
  public :: hpsitopsi

  !-Input variable structure
  !structure of the variables read by input.dat file
  type, public :: input_variables
     integer :: ncount_cluster_x
     real(kind=8) :: frac_fluct,randdis,betax
     integer :: ixc,ncharge,itermax,ncong,idsx,ncongt,inputPsiId,nspin,mpol
     real(kind=8) :: hgrid,crmult,frmult,cpmult,fpmult,elecfield,gnrm_cv,rbuf
     logical :: output_grid,output_wf,calc_tail
  end type input_variables

  !- Interfaces for all outside public routines.
  include "input/interface.f90"
  include "profiling/interface.f90"
  include "interface.f90"

contains
 
  subroutine cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz,&
       psi,keyg,keyv,nvctr_c,nvctr_f,nseg_c,nseg_f,norbp,norb,eval,&
       n1,n2,n3,rxyz_old,in,infocode)
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

  implicit real(kind=8) (a-h,o-z)
  character(len=30) :: label
  character(len=20), dimension(100) :: atomnames
  character(len=1) :: datacode
  logical logrid_c,logrid_f,parallel,calc_tail,output_wf,output_grid
  real(kind=8), parameter :: eps_mach=1.d-12
  ! work array for ALLREDUCE
  dimension wrkallred(5,2) 
  ! atomic coordinates
  dimension rxyz(3,nat),rxyz_old(3,nat),fxyz(3,nat),iatype(nat)
  allocatable :: gxyz(:,:)
  ! active grid points, segments of real space grid
  allocatable :: logrid_c(:,:,:) ,  logrid_f(:,:,:)
  allocatable :: ibyz_c(:,:,:),ibxz_c(:,:,:),ibxy_c(:,:,:),  & 
       ibyz_f(:,:,:),ibxz_f(:,:,:),ibxy_f(:,:,:)
  ! occupation numbers, eigenvalues
  real(kind=8), allocatable :: occup(:),spinar(:),spinar_foo(:)
  real(kind=8), pointer :: eval(:),eval_old(:)

  ! wavefunction segments
  integer, pointer :: keyv(:)
  ! wavefunction segments on real space grid
  integer, pointer :: keyg(:,:)
  ! wavefunction 
  real(kind=8), pointer :: psi(:,:)
  real(kind=8), pointer :: psit(:,:)
  ! wavefunction gradients
  real(kind=8), pointer :: hpsi(:,:)

  ! Pointers and variables to store the last psi
  ! before reformating if useFormattedInput is .true.
  integer :: nseg_c_old, nseg_f_old, nvctr_c_old, nvctr_f_old
  integer, pointer :: keyg_old(:,:), keyv_old(:)
  real(kind=8), pointer :: psi_old(:,:)

  ! Charge density/potential,ionic potential, pkernel
  real(kind=8), allocatable :: rhopot(:,:,:,:),pot_ion(:)
  real(kind=8), pointer     :: pkernel(:)

  ! projector segments on real space grid
  pointer :: keyg_p(:,:), keyv_p(:)
  allocatable :: nvctr_p(:), nseg_p(:)
  ! projectors 
  real(kind=8), pointer :: proj(:)
  ! Parameters for the boxes containing the projectors
  allocatable :: nboxp_c(:,:,:),nboxp_f(:,:,:)

  ! pseudopotential parameters
  allocatable :: psppar(:,:,:),nelpsp(:),radii_cf(:,:),npspcode(:),nzatom(:),iasctype(:)
  allocatable :: derproj(:)
  ! arrays for DIIS convergence accelerator
  real(kind=8), pointer :: ads(:,:,:),psidst(:,:,:),hpsidst(:,:,:)

  ! arrays for calculation of forces and tail correction to kinetic energy
  allocatable :: rho(:),pot(:,:,:,:)
  allocatable :: neleconf(:,:),nscatterarr(:,:),ngatherarr(:,:)

!*****************************Alexey************************************************************
!for shrink:
    integer,allocatable,dimension(:,:,:)::ibzzx_c,ibyyzz_c
    integer,allocatable,dimension(:,:,:)::ibxy_ff,ibzzx_f,ibyyzz_f

!for grow:
    integer,allocatable,dimension(:,:,:)::ibzxx_c,ibxxyy_c
    integer,allocatable,dimension(:,:,:)::ibyz_ff,ibzxx_f,ibxxyy_f
!***********************************************************************************************
    integer,allocatable,dimension(:,:,:)::ibyyzz_r ! real space border
!*************Alexey***************************************************************************
!    real(kind=8),allocatable,dimension(:,:,:)::xc!input 
!    real(kind=8),allocatable::xf(:,:,:,:)! input
!    real(kind=8),allocatable,dimension(:):: w1,w2
!***********************************************************************************************
  integer :: ierror

  !input variables
  type(input_variables), intent(in) :: in

  include 'mpif.h'

  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be 
  !inserted inside the structure

  hgrid=in%hgrid
  crmult=in%crmult
  frmult=in%frmult
  cpmult=in%cpmult
  fpmult=in%fpmult
  ixc=in%ixc
  ncharge=in%ncharge
  elecfield=in%elecfield
  gnrm_cv=in%gnrm_cv
  itermax=in%itermax
  ncong=in%ncong
  idsx=in%idsx
  calc_tail=in%calc_tail
  rbuf=in%rbuf
  ncongt=in%ncongt
  nspin=in%nspin
  mpol=in%mpol

  inputPsiId=in%inputPsiId
  output_grid=in%output_grid
  output_wf=in%output_wf


  if (iproc.eq.0) write(*,'(1x,a,1x,i0)') &
       '                               BigDFT Wavefunction Optimization',inputPsiId
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
           i1=mod(ind-1,nvctrp_old)+1
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

  if(nspin<1.or.nspin>2) nspin=1
  if(nspin==1) mpol=0

  ! grid spacing (same in x,y and z direction)
  hgridh=.5d0*hgrid

  if (iproc==0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------------ System Properties'
  end if

  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(psppar(0:4,0:6,ntypes),stat=i_stat)
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

  call read_system_variables(iproc,nproc,nat,ntypes,nspin,ncharge,mpol,atomnames,iatype,&
       psppar,radii_cf,npspcode,iasctype,nelpsp,nzatom,nelec,natsc,norb,norbu,norbd,norbp,iunit)

  allocate(occup(norb),stat=i_stat)
  call memocc(i_stat,product(shape(occup))*kind(occup),'occup','cluster')
  allocate(spinar(norb),stat=i_stat)
  call memocc(i_stat,product(shape(spinar))*kind(spinar),'occup','cluster')
  
! Occupation numbers
  call input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,occup,spinar)

! Determine size alat of overall simulation cell and shift atom positions
! then calculate the size in units of the grid space
  call system_size(iproc,nat,ntypes,rxyz,radii_cf,crmult,frmult,hgrid,iatype,atomnames, &
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)

  !save the new atomic positions in the rxyz_old array
   do iat=1,nat
     rxyz_old(1,iat)=rxyz(1,iat)
     rxyz_old(2,iat)=rxyz(2,iat)
     rxyz_old(3,iat)=rxyz(3,iat)
   enddo

  !memory estimation
  if (iproc==0) then
     call MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
          rxyz,radii_cf,crmult,frmult,norb,atomnames,.false.,nspin) 
  end if

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=14 !default value
  call createKernel('F',2*n1+31,2*n2+31,2*n3+31,hgridh,hgridh,hgridh,ndegree_ip,&
       iproc,nproc,pkernel)

! Create wavefunctions descriptors and allocate them
  call timing(iproc,'CrtDescriptors','ON')

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
  
!*********************************Alexey*********************************************************
  !   allocate for grow
  allocate(ibzxx_c(2,0:n3,-14:2*n1+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibzxx_c))*kind(ibzxx_c),'ibzxx_c','cluster')
  allocate(ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibxxyy_c))*kind(ibxxyy_c),'ibxxyy_c','cluster')
  allocate(ibyz_ff(2,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
  call memocc(i_stat,product(shape(ibyz_ff))*kind(ibyz_ff),'ibyz_ff','cluster')
  allocate(ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibzxx_f))*kind(ibzxx_f),'ibzxx_f','cluster')
  allocate(ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibxxyy_f))*kind(ibxxyy_f),'ibxxyy_f','cluster')

  !allocate for shrink
  allocate(ibzzx_c(2,-14:2*n3+16,0:n1),stat=i_stat)
  call memocc(i_stat,product(shape(ibzzx_c))*kind(ibzzx_c),'ibzzx_c','cluster')
  allocate(ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibyyzz_c))*kind(ibyyzz_c),'ibyyzz_c','cluster')
  allocate(ibxy_ff(2,nfl1:nfu1,nfl2:nfu2),stat=i_stat)
  call memocc(i_stat,product(shape(ibxy_ff))*kind(ibxy_ff),'ibxy_ff','cluster')
  allocate(ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1),stat=i_stat)
  call memocc(i_stat,product(shape(ibzzx_f))*kind(ibzzx_f),'ibzzx_f','cluster')
  allocate(ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibyyzz_f))*kind(ibyyzz_f),'ibyyzz_f','cluster')

  !allocate for real space
  allocate(ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16),stat=i_stat)
  call memocc(i_stat,product(shape(ibyyzz_r))*kind(ibyyzz_r),'ibyyzz_r','cluster')
!***********************************************************************************************

  call createWavefunctionsDescriptors(iproc,nproc,idsx,n1,n2,n3,output_grid,hgrid,&
       & nat,ntypes,iatype,atomnames,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
       ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp,&
       keyg,keyv,norb,norbp,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
       ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
  call timing(iproc,'CrtDescriptors','OF')

! Calculate all projectors
  call timing(iproc,'CrtProjectors ','ON')

  allocate(nseg_p(0:2*nat),stat=i_stat)
  call memocc(i_stat,product(shape(nseg_p))*kind(nseg_p),'nseg_p','cluster')
  allocate(nvctr_p(0:2*nat),stat=i_stat)
  call memocc(i_stat,product(shape(nvctr_p))*kind(nvctr_p),'nvctr_p','cluster')
  allocate(nboxp_c(2,3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(nboxp_c))*kind(nboxp_c),'nboxp_c','cluster')
  allocate(nboxp_f(2,3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(nboxp_f))*kind(nboxp_f),'nboxp_f','cluster')

  call createProjectorsArrays(iproc,n1,n2,n3,rxyz,nat,ntypes,iatype,atomnames,&
       & psppar,npspcode,radii_cf,cpmult,fpmult,hgrid,nvctr_p,nseg_p,&
       & keyg_p,keyv_p,nproj,nprojel,istart,nboxp_c,nboxp_f,proj)
  call timing(iproc,'CrtProjectors ','OF')
    
  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(nscatterarr(0:nproc-1,4),stat=i_stat)
  call memocc(i_stat,product(shape(nscatterarr))*kind(nscatterarr),'nscatterarr','cluster')
  !allocate array for the communications of the potential
  allocate(ngatherarr(0:nproc-1,2),stat=i_stat)
  call memocc(i_stat,product(shape(ngatherarr))*kind(ngatherarr),'ngatherarr','cluster')

  !create the descriptors for the density and the potential
  call createDensPotDescriptors(iproc,nproc,'F',datacode,n1,n2,n3,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

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
     allocate(rhopot((2*n1+31),(2*n2+31),n3d,nspin),stat=i_stat)
     call memocc(i_stat,product(shape(rhopot))*kind(rhopot),'rhopot','cluster')
  else
     allocate(rhopot(1,1,1,nspin),stat=i_stat)
     call memocc(i_stat,product(shape(rhopot))*kind(rhopot),'rhopot','cluster')
  end if

  allocate(eval(norb),stat=i_stat)
  call memocc(i_stat,product(shape(eval))*kind(eval),'eval','cluster')

! INPUT WAVEFUNCTIONS
  if (inputPsiId == -1) then

     !import gaussians form CP2K (data in files def_gaubasis.dat and gaucoeff.dat)
     !and calculate eigenvalues
     call import_gaussians(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          nat,norb,norbp,occup,n1,n2,n3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
          rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
          atomnames,ntypes,iatype,pkernel,psppar,npspcode,ixc,&
          psi,psit,hpsi,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar,&
          ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
          ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  else if (inputPsiId == 0) then

     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          nat,natsc,norb,norbp,n1,n2,n3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
          rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
          atomnames,ntypes,iatype,iasctype,pkernel,nzatom,nelpsp,psppar,npspcode,&
          ixc,psi,psit,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar,&
          ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
          ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

     if (iproc.eq.0) then
        write(*,'(1x,a,1pe9.2)') 'expected accuracy in kinetic energy due to grid size',accurex
        write(*,'(1x,a,1pe9.2)') 'suggested value for gnrm_cv ',accurex/real(norb,kind=8)
     endif

     if (parallel) then
        !allocate hpsi array (used also as transposed)
        !allocated in the transposed way such as 
        !it can also be used as the transposed hpsi
        allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','cluster')
     else
        !allocate hpsi array
        allocate(hpsi(nvctr_c+7*nvctr_f,norb),stat=i_stat)
        call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','cluster')
     endif
     
  else if (inputPsiId == 1 ) then 
     !restart from previously calculated wavefunctions, in memory

     !allocate principal wavefunction
     if (parallel) then
        !allocated in the transposed way such as 
        !it can also be used as a work array for transposition
        allocate(psi(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     else
        allocate(psi(nvctr_c+7*nvctr_f,norb),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     end if
     
     if (iproc.eq.0) then
        write(*,'(1x,a)')&
          '-------------------------------------------------------------- Wavefunctions Restart'
     end if
     call reformatmywaves(iproc,norb,norbp,nat, &
          & hgrid_old,nvctr_c_old,nvctr_f_old,n1_old,n2_old,n3_old,rxyz_old, &
          & nseg_c_old,nseg_f_old,keyg_old,keyv_old,psi_old, &
          & hgrid,nvctr_c,nvctr_f,n1,n2,n3,rxyz, &
          & nseg_c,nseg_f,keyg,keyv,psi)
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

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
          nspin,psi,hpsi,psit)
     
  else if (inputPsiId == 2 ) then 
     !restart from previously calculated wavefunctions, on disk

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

     if (iproc.eq.0) then
        write(*,'(1x,a)')&
          '---------------------------------------------------- Reading Wavefunctions from disk'
     end if


     call readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,&
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
  
     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0
 
     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
          nspin,psi,hpsi,psit)

  else
     
     if (iproc == 0) then
        write(*,'(1x,a)')'The supported values of inputPsiId are integers from -1 to 2'
        write(*,'(1x,a,i0)')'                        while we found',inputPsiId
     end if
     stop
    
  end if

  !no need of using nzatom array, semicores useful only for the input guess
  i_all=-product(shape(nzatom))*kind(nzatom)
  deallocate(nzatom,stat=i_stat)
  call memocc(i_stat,i_all,'nzatom','cluster')
  i_all=-product(shape(iasctype))*kind(iasctype)
  deallocate(iasctype,stat=i_stat)
  call memocc(i_stat,i_all,'iasctype','cluster')


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
  energy=1.d10
  gnrm=1.d10
  ekin_sum=0.d0 
  epot_sum=0.d0 
  eproj_sum=0.d0
  !set the infocode to the value it would have in the case of no convergence
  infocode=1
! loop for wavefunction minimization
  wfn_loop: do iter=1,itermax
     if (idsx.gt.0) mids=mod(iter-1,idsx)+1
     if (iproc.eq.0) then 
        write(*,'(1x,a,i0)')&
         '---------------------------------------------------------------------------- iter= ',&
         iter
     endif

     ! Potential from electronic charge density
     call sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
             nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot, &
             (2*n1+31)*(2*n2+31)*n3d,nscatterarr,nspin,spinar,&
             nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
             ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

!     ixc=11  ! PBE functional
!     ixc=1   ! LDA functional

     call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin) !add NSPIN

     call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
             psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
             nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
             nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,ngatherarr,n3p,&
             rhopot(1,1,1+i3xcsh,1),psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,spinar,&
             ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
             ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

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
           !write(61,*)hgrid,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
        end if
        infocode=0
        exit wfn_loop 
     endif

     call hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp,nseg_c,nseg_f,&
     keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
     psi,psit,hpsi,psidst,hpsidst,nspin,spinar)! add NSPIN

     tt=energybs-scprsum
     if (abs(tt).gt.1.d-8 .and. iproc==0) then 
        write(*,'(1x,a,3(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
     endif
     if (iproc.eq.0) then
        write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             ekin_sum,epot_sum,eproj_sum
        write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
        write(*,'(1x,a,i6,2x,1pe19.12,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
     endif

     if (inputPsiId == 0) then
        if (gnrm > 4.d0) then
           if (iproc == 0) then
              write(*,'(1x,a)')&
                'Error: the norm of the residue is too large also with input wavefunctions.'
           end if
           infocode=3
           call deallocate_before_exiting
           return
        end if
     else if (inputPsiId == 1) then
        if (gnrm > 1.d0) then
           if (iproc == 0) then
              write(*,'(1x,a)')&
                'The norm of the residue is too large, need to recalculate input wavefunctions'
           end if
           infocode=2
           if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           call deallocate_before_exiting
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
              if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              call deallocate_before_exiting
              return
           end if
        end if
     end if

  end do wfn_loop
     if (iter == itermax+1 .and. iproc == 0 ) &
          write(*,'(1x,a)')'No convergence within the allowed number of minimization steps'

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

! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
  call last_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
       nspin,psi,hpsi,psit,occup,evsum,eval)

  if (abs(evsum-energybs).gt.1.d-8 .and. iproc==0) write(*,'(1x,a,2(1x,1pe20.13))')&
       'Difference:evsum,energybs',evsum,energybs

!!$  !plot the converged wavefunctions in the different orbitals
!!$  do i=2*iproc+1,2*iproc+2
!!$     iounit=39+3*(i-1)
!!$     call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$          rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc),&
!!$          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
!!$          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
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

  allocate(spinar_foo(norb),stat=i_stat)
  call memocc(i_stat,product(shape(spinar_foo))*kind(spinar_foo),'spinar_foo','cluster')
  spinar_foo(:)=1.0d0
  ! Potential from electronic charge density
  
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
       (2*n1+31)*(2*n2+31)*n3p,nscatterarr,1,spinar_foo,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  i_all=-product(shape(spinar_foo))*kind(spinar_foo)
  deallocate(spinar_foo,stat=i_stat)
  call memocc(i_stat,i_all,'spinar_foo','cluster')
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
  
  !switch between the old and the new forces calculation
  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion','cluster')

  if (n3p>0) then
     allocate(pot((2*n1+31),(2*n2+31),n3p,1),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')
  else
     allocate(pot(1,1,1,1),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')
  end if
  call DCOPY((2*n1+31)*(2*n2+31)*n3p,rho,1,pot,1) 
  call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,0,hgridh,hgridh,hgridh,&
          pot,pkernel,pot,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.,1)
     !here nspin=1 since ixc=0


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
  !diring the wavefunction minimization if needed
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
     allocate(pot((2*n1+31),(2*n2+31),(2*n3+31),nspin),stat=i_stat)
     call memocc(i_stat,product(shape(pot))*kind(pot),'pot','cluster')

     if (datacode=='D') then
        call MPI_ALLGATHERV(rhopot(1,1,1+i3xcsh,1),(2*n1+31)*(2*n2+31)*n3p,&
             MPI_DOUBLE_PRECISION,pot(1,1,1,1),ngatherarr(0,1),ngatherarr(0,2), & 
             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        !print '(a,2f12.6)','RHOup',sum(abs(rhopot(:,:,:,1))),sum(abs(pot(:,:,:,1)))
        if(nspin==2) then
           !print '(a,2f12.6)','RHOdw',sum(abs(rhopot(:,:,:,2))),sum(abs(pot(:,:,:,2)))
           if (n3d /= n3p) then
              i03=1+i3xcsh+n3p
              i04=1
           else
              i03=1
              i04=2
           end if
           call MPI_ALLGATHERV(rhopot(1,1,i03,i04),(2*n1+31)*(2*n2+31)*n3p,&
                MPI_DOUBLE_PRECISION,pot(1,1,1,2),ngatherarr(0,1),ngatherarr(0,2), & 
                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        end if
     else
        do ispin=1,nspin
        !here one could have not allocated pot and: call move_alloc(rhopot,pot) 
        !(but it is a Fortran 95/2003 spec)
           do i3=1,2*n3+31
              do i2=1,2*n2+31
                 do i1=1,2*n1+31
                    pot(i1,i2,i3,ispin)=rhopot(i1,i2,i3,ispin)
                 enddo
              enddo
           enddo
        end do
     end if
     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
     deallocate(nscatterarr,stat=i_stat)
     call memocc(i_stat,i_all,'nscatterarr','cluster')
     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr','cluster')
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot','cluster')

     call CalculateTailCorrection(iproc,nproc,n1,n2,n3,rbuf,norb,norbp,nat,ntypes,&
     nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nproj,nprojel,ncongt,&
     keyv,keyg,nseg_p,keyv_p,keyg_p,nvctr_p,psppar,npspcode,eval,&
     pot,hgrid,rxyz,radii_cf,crmult,frmult,iatype,atomnames,nspin,spinar,&
     proj,psi,occup,output_grid,parallel,ekin_sum,epot_sum,eproj_sum)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot','cluster')

     !if (iproc==0) then
     !   open(61)
     !   write(61,'(4(f9.3),1x,7(1pe19.11))',advance='no')&
     !        hgrid,alat1,alat2,alat3,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
     !end if
    
     energybs=ekin_sum+epot_sum+eproj_sum
     energy=energybs-ehart+eexcu-vexcu+eion

     !if (iproc==0) then
     !   write(61,'(1pe19.11)')energy
     !   close(61)
     !end if

     if (iproc.eq.0) then
        write(*,'(1x,a,3(1x,1pe18.11))')&
             '  Corrected ekin,epot,eproj',ekin_sum,epot_sum,eproj_sum
        write(*,'(1x,a,1x,1pe19.12)')&
             'Total energy with tail correction',energy
     endif

     call timing(iproc,'Tail          ','OF')
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

  call deallocate_before_exiting

contains

  !routine which deallocate the pointers and the arrays before exiting 
  subroutine deallocate_before_exiting
    implicit real(kind=8) (a-h,o-z)

    !when this condition is verified we are in the middle of the SCF cycle
    if (infocode /=0 .and. infocode /=1) then

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

       if (parallel) then
          i_all=-product(shape(psit))*kind(psit)
          deallocate(psit,stat=i_stat)
          call memocc(i_stat,i_all,'psit','cluster')
       end if

       i_all=-product(shape(hpsi))*kind(hpsi)
       deallocate(hpsi,stat=i_stat)
       call memocc(i_stat,i_all,'hpsi','cluster')

       i_all=-product(shape(pot_ion))*kind(pot_ion)
       deallocate(pot_ion,stat=i_stat)
       call memocc(i_stat,i_all,'pot_ion','cluster')

       i_all=-product(shape(pkernel))*kind(pkernel)
       deallocate(pkernel,stat=i_stat)
       call memocc(i_stat,i_all,'pkernel','cluster')

       i_all=-product(shape(nboxp_c))*kind(nboxp_c)
       deallocate(nboxp_c,stat=i_stat)
       call memocc(i_stat,i_all,'nboxp_c','cluster')
       i_all=-product(shape(nboxp_f))*kind(nboxp_f)
       deallocate(nboxp_f,stat=i_stat)
       call memocc(i_stat,i_all,'nboxp_f','cluster')

       ! calc_tail false
       i_all=-product(shape(rhopot))*kind(rhopot)
       deallocate(rhopot,stat=i_stat)
       call memocc(i_stat,i_all,'rhopot','cluster')
       i_all=-product(shape(nscatterarr))*kind(nscatterarr)
       deallocate(nscatterarr,stat=i_stat)
       call memocc(i_stat,i_all,'nscatterarr','cluster')
       i_all=-product(shape(ngatherarr))*kind(ngatherarr)
       deallocate(ngatherarr,stat=i_stat)
       call memocc(i_stat,i_all,'ngatherarr','cluster')

    end if

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

    !*****************************Alexey*************************
    i_all=-product(shape(ibzzx_c))*kind(ibzzx_c)
    deallocate(ibzzx_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibzzx_c','cluster')
    i_all=-product(shape(ibyyzz_c))*kind(ibyyzz_c)
    deallocate(ibyyzz_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_c','cluster')
    i_all=-product(shape(ibxy_ff))*kind(ibxy_ff)
    deallocate(ibxy_ff,stat=i_stat)
    call memocc(i_stat,i_all,'ibxy_ff','cluster')
    i_all=-product(shape(ibzzx_f))*kind(ibzzx_f)
    deallocate(ibzzx_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibzzx_f','cluster')
    i_all=-product(shape(ibyyzz_f))*kind(ibyyzz_f)
    deallocate(ibyyzz_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_f','cluster')
    i_all=-product(shape(ibzxx_c))*kind(ibzxx_c)
    deallocate(ibzxx_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibzxx_c','cluster')
    i_all=-product(shape(ibxxyy_c))*kind(ibxxyy_c)
    deallocate(ibxxyy_c,stat=i_stat)
    call memocc(i_stat,i_all,'ibxxyy_c','cluster')
    i_all=-product(shape(ibyz_ff))*kind(ibyz_ff)
    deallocate(ibyz_ff,stat=i_stat)
    call memocc(i_stat,i_all,'ibyz_ff','cluster')
    i_all=-product(shape(ibzxx_f))*kind(ibzxx_f)
    deallocate(ibzxx_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibzxx_f','cluster')
    i_all=-product(shape(ibxxyy_f))*kind(ibxxyy_f)
    deallocate(ibxxyy_f,stat=i_stat)
    call memocc(i_stat,i_all,'ibxxyy_f','cluster')
    i_all=-product(shape(ibyyzz_r))*kind(ibyyzz_r)
    deallocate(ibyyzz_r,stat=i_stat)
    call memocc(i_stat,i_all,'ibyyzz_r','cluster')
    !************************************************************

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
    i_all=-product(shape(spinar))*kind(spinar)
    deallocate(spinar,stat=i_stat)
    call memocc(i_stat,i_all,'spinar','cluster')
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
    if (iproc == 0) &
         write(*,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time for root process ', iproc,tel,tcpu1-tcpu0

  end subroutine deallocate_before_exiting

END SUBROUTINE cluster

subroutine hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp,nseg_c,nseg_f,&
     keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
     psi,psit,hpsi,psidst,hpsidst,nspin,spinar)
  implicit none
  include 'mpif.h'
  real(kind=8), parameter :: eps_mach=1.d-12
  logical, intent(in) :: parallel
  integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
  integer, intent(in) :: nseg_c,nseg_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp
  real(kind=8), intent(in) :: hgrid,energy,energy_old
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  real(kind=8), dimension(norb), intent(in) :: occup,eval,spinar
  real(kind=8), intent(inout) :: alpha
  real(kind=8), intent(inout) :: gnrm,scprsum
  real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
! real(kind=8), dimension(:,:), pointer :: psit
! real(kind=8), dimension(nvctr_c+7*nvctr_f,norb) :: psi,hpsi
  real(kind=8), dimension(:,:,:), pointer :: psidst,hpsidst,ads
  !local variables
  integer :: ierr,ind,i1,i2,iorb,k,norbu,norbd,nspin
  real(kind=8) :: tt,scpr,dnrm2,scprpart

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, orthoconstraint...'
  end if
  
  !Calculate no. up and dw orbitals for spin-polarized starting guess
  norbu=0
  norbd=0
  do iorb=1,norb
     if(spinar(iorb)>0.0d0) norbu=norbu+1
     if(spinar(iorb)<0.0d0) norbd=norbd+1
  end do
  !write(*,'(1x,a,3i4,30f6.2)')'Spins: ',norb,norbu,norbd,(spinar(iorb),iorb=1,norb)

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  if (parallel) then
     !transpose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !here hpsi is the transposed array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition


     if(nspin==1) then
        call  orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsi,scprsum)
     else
        call  orthoconstraint_p(iproc,nproc,norbu,occup,nvctrp,psit,hpsi,scprsum)
        scprpart=0.0d0
        if(norbd>0) then
           scprpart=scprsum 
           call  orthoconstraint_p(iproc,nproc,norbd,occup,nvctrp,psit(1,norbu+1),hpsi(1,norbu+1),scprsum)
        end if
        scprsum=scprsum+scprpart
     end if

     !retranspose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !here hpsi is the direct array
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,hpsi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
        if(nspin==1) then
           call orthoconstraint(norb,occup,nvctrp,psi,hpsi,scprsum)
        else
           call orthoconstraint(norbu,occup,nvctrp,psi,hpsi,scprsum)
           scprpart=0.0d0
           if(norbd>0) then
              scprpart=scprsum 
              call orthoconstraint(norbd,occup,nvctrp,psi(1,norbu+1),hpsi(1,norbu+1),scprsum)
           end if
           scprsum=scprsum+scprpart
        end if
  endif


  ! norm of gradient
  gnrm=0.d0
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     !calculate the address to start from for calculating the 
     !norm of the residue if hpsi is allocated in the transposed way
     !it can be eliminated when including all this procedure in a subroutine
     if (parallel) then
        ind=1+(nvctr_c+7*nvctr_f)*(iorb-iproc*norbp-1)
        i1=mod(ind-1,nvctrp)+1
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
  gnrm=sqrt(gnrm/real(norb,kind=8))

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'done, preconditioning...'
  end if

  !   write(*,'(10f10.6)') (eval(iorb),iorb=1,norb)
  ! Preconditions all orbitals belonging to iproc
  call preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid, &
       ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,&
       ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)



  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')


!!$  if (iproc==0) then
!!$       do iorb=1,norb
!!$        do i1=1,nvctr_c+7*nvctr_f
!!$!          write(200,'(f12.8)') psi(i1,iorb)
!!$!          write(220,'(f12.8)') hpsi(i1,iorb)
!!$        end do
!!$     end do
!!$  end if

  !apply the minimization method (DIIS or steepest descent)
  if (idsx.gt.0) then
     if (parallel) then
        !transpose the hpsi wavefunction into the diis array

        !here psi is used as a work array
        call timing(iproc,'Un-TransSwitch','ON')
        call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
        call timing(iproc,'Un-TransSwitch','OF')
        call timing(iproc,'Un-TransComm  ','ON')
        call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsidst(1,1,mids),nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-TransComm  ','OF')
        !end of transposition

        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           do k=1,nvctrp
              psidst(k,iorb,mids)= psit(k,iorb) 
           enddo
        enddo

        call diisstp(parallel,norb,norbp,nproc,iproc,  &
             ads,iter,mids,idsx,nvctrp,psit,psidst,hpsidst)
     else
        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           do k=1,nvctrp
              psidst(k,iorb,mids)= psi(k,iorb)
              hpsidst(k,iorb,mids)=hpsi(k,iorb)
           enddo
        enddo

        call diisstp(parallel,norb,norbp,nproc,iproc,  &
             ads,iter,mids,idsx,nvctrp,psi,psidst,hpsidst)

     endif
!    write(*,*) 'psi update done',iproc

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

        !here psi is used as a work array
        call timing(iproc,'Un-TransSwitch','ON')
        call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
        call timing(iproc,'Un-TransSwitch','OF')
        !here hpsi is the transposed array
        call timing(iproc,'Un-TransComm  ','ON')
        call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
             hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call timing(iproc,'Un-TransComm  ','OF')
        !end of transposition

        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psit(1,iorb),1)
        enddo
     else
        call timing(iproc,'Diis          ','ON')
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psi(1,iorb),1)
        enddo
     endif

  endif

  call timing(iproc,'Diis          ','OF')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Orthogonalization...'
  end if
       do iorb=1,norb
        do i1=1,nvctr_c+7*nvctr_f
!          write(201,'(f12.8)') psi(i1,iorb)
!          write(210,'(f12.8)') hpsi(i1,iorb)
        end do
     end do

  if (parallel) then
     if(nspin==1) then
        call orthon_p(iproc,nproc,norb,nvctrp,nvctr_c+7*nvctr_f,psit)
     else
        call orthon_p(iproc,nproc,norbu,nvctrp,nvctr_c+7*nvctr_f,psit)
        if(norbd>0) then
           call orthon_p(iproc,nproc,norbd,nvctrp,nvctr_c+7*nvctr_f,psit(1,norbu+1))
        end if
     end if
     !       call checkortho_p(iproc,nproc,norb,nvctrp,psit)

     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     if(nspin==1) then
        call orthon(norb,nvctrp,psi)
        !          call checkortho(norb,nvctrp,psi)
     else
        call orthon(norbu,nvctrp,psi)
        !          call checkortho(norbu,nvctrp,psi)
        if(norbd>0) then
           call orthon(norbd,nvctrp,psi(1,norbu+1))
        end if
        !          call checkortho(norbd,nvctrp,psi(1,norbu+1))
     end if
  endif
  
       do iorb=1,norb
        do i1=1,nvctr_c+7*nvctr_f
!           write(200,'(f12.8)') psi(i1,iorb)
!          write(210,'(f12.8)') hpsi(i1,iorb)
        end do
     end do

  if (iproc==0) then
     write(*,'(1x,a)')&
          'done.'
  end if

end subroutine hpsitopsi

subroutine first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
     nspin,psi,hpsi,psit)
  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
  real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
  !local variables
  integer :: i_all,i_stat,ierr

  if (parallel) then
     !allocate hpsi array (used also as transposed)
     !allocated in the transposed way such as 
     !it can also be used as the transposed hpsi
     allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','first_orthon')

     !transpose the psi wavefunction
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,hpsi)
     call timing(iproc,'Un-TransSwitch','OF')
     !allocate transposed principal wavefunction
     allocate(psit(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psit))*kind(psit),'psit','first_orthon')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition

     if(nspin==1) then
        call orthon_p(iproc,nproc,norb,nvctrp,nvctr_c+7*nvctr_f,psit) 
     else
        call orthon_p(iproc,nproc,norbu,nvctrp,nvctr_c+7*nvctr_f,psit) 
        if(norbd>0) then
           call orthon_p(iproc,nproc,norbd,nvctrp,nvctr_c+7*nvctr_f,psit(1,norbu+1)) 
        end if
     end if
     !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

     !retranspose the psit wavefunction into psi
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition
  else
     if(nspin==1) then
        call orthon(norb,nvctrp,psi)
        !call checkortho(norb,norbp,nvctrp,psi)
     else
        call orthon(norbu,nvctrp,psi)
        !call checkortho(norbu,nvctrp,psi)
        if(norbd>0) then
           call orthon(norbd,nvctrp,psi(1,norbu+1))
           !call checkortho(norbd,nvctrp,psi(1,norbu+1))
        end if
     end if
     !allocate hpsi array
     allocate(hpsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'hpsi','first_orthon')
  endif

end subroutine first_orthon

! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
     nspin,psi,hpsi,psit,occup,evsum,eval)
  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), intent(out) :: evsum
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
  !local variables
  integer :: i_all,i_stat,ierr,iorb,jorb
  real(kind=8) :: evpart

  if (parallel) then

     !transpose the hpsi wavefunction
     !here psi is used as a work array
     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !here hpsi is the transposed array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     !end of transposition

     if(nspin==1) then
        call KStrans_p(iproc,nproc,norb,norbp*nproc,nvctrp,occup,hpsi,psit,evsum,eval)
     else
        call KStrans_p(iproc,nproc,norbu,norbu,nvctrp,occup,hpsi,psit,evsum,eval)
        evpart=evsum
        if(norbd>0) then
           call KStrans_p(iproc,nproc,norbd,norbd,nvctrp,occup,&
                hpsi(1,norbu+1),psit(1,norbu+1),evsum,eval(norbu+1))
           evsum=evsum+evpart
        end if
     end if

     !retranspose the psit wavefunction into psi
     !here hpsi is used as a work array
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')
     !end of retransposition

     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit','last_orthon')

  else
     if(nspin==1) then
        call KStrans(norb,nvctrp,occup,hpsi,psi,evsum,eval)
     else
        call KStrans(norbu,nvctrp,occup,hpsi,psi,evsum,eval)
        evpart=evsum
        if(norbd>0) then
           call KStrans(norbd,nvctrp,occup,hpsi(1,norbu+1),psi(1,norbu+1),&
                evsum,eval(norbu+1))
           evsum=evsum+evpart
        end if
     end if
  endif

  !print the found eigenvalues
  if (iproc == 0) then
     write(*,'(1x,a)')&
          '-------------------------------------------------------------- Kohn-Sham Eigenvalues'
     if (nspin==1) then
        do iorb=1,norb!/2
!!$           jorb=norb/2+iorb
!!$           write(*,'(1x,a,i4,a,1x,1pe21.14,17x,a,i4,a,1x,1pe21.14)') &
!!$                'eval(',iorb,')=',eval(iorb),'eval(',jorb,')=',eval(jorb)
!!$        end do
!!$        if (2*norb/2 /= norb) then
!!$           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',norb/2+1,')=',eval(norb/2+1)
!!$        end if
           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        end do
     else
        do iorb=1,min(norbu,norbd)
           jorb=norbu+iorb
           write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
                'eval(',iorb,',u)=',eval(iorb),'eval(',iorb,',d)=',eval(jorb)
        end do
        if (norbu > norbd) then
           do iorb=norbd+1,norbu
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'eval(',iorb,',u)=',eval(iorb)
           end do
        else if (norbd > norbu) then
           do iorb=2*norbu+1,norbu+norbd
              write(*,'(48x,a,i0,a,1x,1pe21.14)') 'eval(',iorb-norbu,',d)=',eval(iorb)
           end do
        end if
     end if
  end if

  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi','last_orthon')

end subroutine last_orthon


!create the descriptors for the density and the potential
subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1,n2,n3,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  use Poisson_Solver

  implicit none
  character(len=1), intent(in) :: geocode,datacode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,ixc
  integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
  integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
  !local variables
  integer :: jproc

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

  ngatherarr(:,1)=(2*n1+31)*(2*n2+31)*nscatterarr(:,2)
  ngatherarr(:,2)=(2*n1+31)*(2*n2+31)*nscatterarr(:,3)

end subroutine createDensPotDescriptors


subroutine createWavefunctionsDescriptors(iproc,nproc,idsx,n1,n2,n3,output_grid,&
     hgrid,nat,ntypes,iatype,atomnames,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp,&
     keyg, keyv,norb,norbp,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
!calculates the descriptor arrays keyg and keyv as well as nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp
!calculates also the arrays ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f needed for convolut_standard
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,idsx,n1,n2,n3,nat,ntypes,norb,norbp
  integer, intent(out) :: nseg_c,nseg_f,nvctr_c,nvctr_f
  integer, intent(out) :: nvctrp
  logical, intent(in) :: output_grid
  integer, intent(in) :: iatype(nat)
  real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
  integer, intent(out) :: ibyz_c(2,0:n2,0:n3), ibxz_c(2,0:n1,0:n3), ibxy_c(2,0:n1,0:n2)
  integer, intent(out) :: ibyz_f(2,0:n2,0:n3), ibxz_f(2,0:n1,0:n3), ibxy_f(2,0:n1,0:n2)
  real(kind=8) :: rxyz(3, nat), radii_cf(ntypes, 2)
  character(len=20), intent(in) :: atomnames(100)
  integer, pointer :: keyg(:,:), keyv(:)
  !**********************Alexey*************************************************************
  integer,intent(in):: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  !**********************Alexey*************************************************************
  !for shrink:
  integer ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  !for grow:
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  !for real space:
  integer,intent(out):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
  !*****************************************************************************************


  !Local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: iat,i1,i2,i3,norbme,norbyou,jpst,jproc,i_all,i_stat
  real(kind=8) :: tt
  logical, allocatable :: logrid_c(:,:,:), logrid_f(:,:,:)

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------- Wavefunctions Descriptors Creation'
  end if

! Create the file grid.ascii to visualize the grid of functions
  if (iproc.eq.0 .and. output_grid) then
     open(unit=22,file='grid.ascii',status='unknown')
     write(22,*) nat
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     do iat=1,nat
        write(22,'(3(1x,e12.5),3x,a20)') &
             rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
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
                   write(22,'(3(1x,e10.3),1x,a4)') &
                        real(i1,kind=8)*hgrid,real(2,kind=8)*hgrid,real(i3,kind=8)*hgrid,'  g '
           enddo
        enddo
     end do
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,nvctr_c)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       'Coarse resolution grid: Number of segments= ',nseg_c,'points=',nvctr_c
  call bounds(n1,n2,n3,logrid_c,ibyz_c,ibxz_c,ibxy_c)

  ! fine grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,2),frmult,hgrid,logrid_f)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') &
                        real(i1,kind=8)*hgrid,real(i2,kind=8)*hgrid,real(i3,kind=8)*hgrid,'  G '
           enddo
        enddo
     enddo
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,nvctr_f)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       '  Fine resolution grid: Number of segments= ',nseg_f,'points=',nvctr_f
  call bounds(n1,n2,n3,logrid_f,ibyz_f,ibxz_f,ibxy_f)

  if (iproc.eq.0 .and. output_grid) close(22)

  ! allocations for arrays holding the wavefunctions and their data descriptors
  allocate(keyg(2,nseg_c+nseg_f),stat=i_stat)
  call memocc(i_stat,product(shape(keyg))*kind(keyg),'keyg','createwavefunctionsdescriptors')
  allocate(keyv(nseg_c+nseg_f),stat=i_stat)
  call memocc(i_stat,product(shape(keyv))*kind(keyv),'keyv','createwavefunctionsdescriptors')

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,nseg_c,keyg(1,1),keyv(1))

  ! fine grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,nseg_f,keyg(1,nseg_c+1), &
    & keyv(nseg_c+1))

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c','createwavefunctionsdescriptors')
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f','createwavefunctionsdescriptors')

! allocate wavefunction arrays
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

!*********Alexey******************************************************************************
 
  call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
       ibyz_c,ibzxx_c,ibxxyy_c,ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

!***********************************************************************************************
END SUBROUTINE createWavefunctionsDescriptors

subroutine createProjectorsArrays(iproc, n1, n2, n3, rxyz, nat, ntypes, iatype, atomnames, &
     & psppar, npspcode, radii_cf, cpmult, fpmult, hgrid, nvctr_p, nseg_p, &
     & keyg_p, keyv_p, nproj, nprojel, istart, nboxp_c, nboxp_f, proj)
  implicit real(kind=8) (a-h,o-z)
  character(len=20) :: atomnames(100)
  dimension rxyz(3,nat),iatype(nat),radii_cf(ntypes,2),psppar(0:4,0:6,ntypes),npspcode(ntypes)
  integer :: nvctr_p(0:2*nat), nseg_p(0:2*nat)
  integer :: nboxp_c(2,3,nat), nboxp_f(2,3,nat)
  real(kind=8), pointer :: proj(:)
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
     do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
        do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
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
                 if (abs(1.d0-scpr).gt.1.d-2) then
                    print *,'norm projector for atom ',trim(atomnames(iatype(iat))),&
                         'iproc,l,i,rl,scpr=',iproc,l,i,gau_a,scpr
                    stop 'norm projector'
                 end if

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

END SUBROUTINE createProjectorsArrays

subroutine import_gaussians(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nat,norb,norbp,occup,n1,n2,n3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
     rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
     atomnames,ntypes,iatype,pkernel,psppar,npspcode,ixc,&
     psi,psit,hpsi,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar, &
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  use Poisson_Solver

  implicit none
  include 'mpif.h'
  logical, intent(in) :: parallel
  character(len=20), dimension(100), intent(in) :: atomnames
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,nprojel,nproj,ixc
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctr_c,nvctr_f,nvctrp,nseg_c,nseg_f
  integer, intent(in) :: nspin
  real(kind=8), dimension(norb), intent(in) :: spinar
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(out) :: accurex
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: npspcode
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
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(nprojel), intent(in) :: proj
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
  !real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: ppsi
  !********************Alexey***************************************************************
  !for shrink:
  integer ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  !for grow:
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  !for real space:
  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
  !*****************************************************************************************

  !local variables
  integer :: i,iorb,i_stat,i_all,ierr,info,jproc,n_lp,jorb
  real(kind=8) :: hgridh,tt,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum
  real(kind=8), dimension(:), allocatable :: work_lp,pot,ones
  real(kind=8), dimension(:,:), allocatable :: hamovr

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------------- Import Gaussians from CP2K'
  end if

  hgridh=.5d0*hgrid

  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(psi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','import_gaussians')
  else
     allocate(psi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','import_gaussians')
  end if

  !read the values for the gaussian code and insert them on psi 
 call gautowav(iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       nvctr_c,nvctr_f,nseg_c,nseg_f,keyg,keyv,iatype,occup,rxyz,hgrid,psi,eks)

!!$  !!plot the initial gaussian wavefunctions
!!$  !do i=2*iproc+1,2*iproc+2
!!$  !   iounit=15+3*(i-1)
!!$  !   print *,'iounit',iounit,'-',iounit+2
!!$  !   call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$  !        rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc), & 
!!$          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
!!$          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!!$  !end do
 

 ! resulting charge density and potential
 allocate(ones(norb),stat=i_stat)
 call memocc(i_stat,product(shape(ones))*kind(ones),'ones','import_gaussians')
 ones(:)=1.0d0

 call sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
      nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot,&
      (2*n1+31)*(2*n2+31)*nscatterarr(iproc,1),nscatterarr,1,ones,&
      nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

 call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
      rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,1)


 if (parallel) then
    !allocate the wavefunction in the transposed way to avoid allocations/deallocations
    allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
    call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','import_gaussians')
 else
    allocate(hpsi(nvctr_c+7*nvctr_f,norbp),stat=i_stat)
    call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','import_gaussians')
 end if


 call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
      psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
      nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
      nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,ngatherarr,nscatterarr(iproc,2),&
      rhopot(1+(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)),&
      psi,hpsi,ekin_sum,epot_sum,eproj_sum,1,ones,&
      ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
      ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  i_all=-product(shape(ones))*kind(ones)
  deallocate(ones,stat=i_stat)
  call memocc(i_stat,i_all,'ones','import_gaussians')

 
 accurex=abs(eks-ekin_sum)
 if (iproc.eq.0) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks


 !after having applied the hamiltonian to all the atomic orbitals
 !we split the semicore orbitals from the valence ones
 !this is possible since the semicore orbitals are the first in the 
 !order, so the linear algebra on the transposed wavefunctions 
 !may be splitted

 if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
      'Imported Wavefunctions Orthogonalization:'

 if (parallel) then

    !transpose all the wavefunctions for having a piece of all the orbitals 
    !for each processor
    !allocate the wavefunction in the transposed way to avoid allocations/deallocations
    allocate(psit(nvctrp,norbp*nproc),stat=i_stat)
    call memocc(i_stat,product(shape(psit))*kind(psit),'psit','import_gaussians')

    call timing(iproc,'Un-TransSwitch','ON')
    call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
    call timing(iproc,'Un-TransSwitch','OF')
    call timing(iproc,'Un-TransComm  ','ON')
    call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
         psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call timing(iproc,'Un-TransComm  ','OF')

    call timing(iproc,'Un-TransSwitch','ON')
    call switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psit)
    call timing(iproc,'Un-TransSwitch','OF')
    call timing(iproc,'Un-TransComm  ','ON')
    call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
         hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call timing(iproc,'Un-TransComm  ','OF')

    !end of transposition

    allocate(hamovr(norb**2,4),stat=i_stat)
    call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','import_gaussians')

    !calculate the overlap matrix for each group of the semicore atoms
    !       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
    !       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)

    if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
         'Overlap Matrix...'

    call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,&
         0.d0,hamovr(1,3),norb)

    call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,psi,nvctrp,&
         0.d0,hamovr(1,4),norb)

    !reduce the overlap matrix between all the processors
    call MPI_ALLREDUCE(hamovr(1,3),hamovr(1,1),2*norb**2,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    !print the overlap matrix in the wavelet case
    print *,norb
    open(33)
    do iorb=1,norb
       write(33,'(2000(1pe10.2))')&
            (hamovr(jorb+(iorb-1)*norb,2),jorb=1,norb)
    end do
    close(33)

    !stop

    !found the eigenfunctions for each group
    n_lp=max(10,4*norb)
    allocate(work_lp(n_lp),stat=i_stat)
    call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','import_gaussians')

    if (iproc.eq.0) write(*,'(1x,a)')'Linear Algebra...'

    call DSYGV(1,'V','U',norb,hamovr(1,1),norb,hamovr(1,2),norb,eval,work_lp,n_lp,info)

    if (info.ne.0) write(*,*) 'DSYGV ERROR',info
!!$        !!write the matrices on a file
!!$        !open(33+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(33+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(33+2*(i-1))
!!$        !open(34+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(34+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jjorb+(jiorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(34+2*(i-1))

    if (iproc.eq.0) then
       do iorb=1,norb
          write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
       enddo
    endif

    i_all=-product(shape(work_lp))*kind(work_lp)
    deallocate(work_lp,stat=i_stat)
    call memocc(i_stat,i_all,'work_lp','import_gaussians')

    if (iproc.eq.0) write(*,'(1x,a)',advance='no')'Building orthogonal Imported Wavefunctions...'

    !perform the vector-matrix multiplication for building the input wavefunctions
    ! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)

    call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi,nvctrp,&
         hamovr(1,1),norb,0.d0,psit,nvctrp)

    i_all=-product(shape(hamovr))*kind(hamovr)
    deallocate(hamovr,stat=i_stat)
    call memocc(i_stat,i_all,'hamovr','import_gaussians')

    !retranspose the wavefunctions
    call timing(iproc,'Un-TransComm  ','ON')
    call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
         hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call timing(iproc,'Un-TransComm  ','OF')
    call timing(iproc,'Un-TransSwitch','ON')
    call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,psi)
    call timing(iproc,'Un-TransSwitch','OF')

 else !serial case

    write(*,'(1x,a)',advance='no')'Overlap Matrix...'

    allocate(hamovr(norb**2,4),stat=i_stat)
    call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','import_gaussians')
    !hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
    call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,&
         0.d0,hamovr(1,1),norb)
    call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,psi,nvctrp,&
         0.d0,hamovr(1,2),norb)

    n_lp=max(10,4*norb)
    allocate(work_lp(n_lp),stat=i_stat)
    call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','import_gaussians')

    write(*,'(1x,a)')'Linear Algebra...'
    call DSYGV(1,'V','U',norb,hamovr(1,1),norb,hamovr(1,2),norb,eval,work_lp,n_lp,info)

    if (info.ne.0) write(*,*) 'DSYGV ERROR',info
    if (iproc.eq.0) then
       do iorb=1,norb
          write(*,'(1x,a,i0,a,1x,1pe21.14)') 'evale(',iorb,')=',eval(iorb)
       enddo
    endif

    i_all=-product(shape(work_lp))*kind(work_lp)
    deallocate(work_lp,stat=i_stat)
    call memocc(i_stat,i_all,'work_lp','import_gaussians')

    write(*,'(1x,a)',advance='no')'Building orthogonal Imported Wavefunctions...'

    !copy the values into hpsi
    do iorb=1,norb
       do i=1,nvctr_c+7*nvctr_f
          hpsi(i,iorb)=psi(i,iorb)
       end do
    end do
    !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
    call DGEMM('N','N',nvctrp,norb,norb,1.d0,hpsi,nvctrp,hamovr(1,1),norb,0.d0,psi,nvctrp)

    i_all=-product(shape(hamovr))*kind(hamovr)
    deallocate(hamovr,stat=i_stat)
    call memocc(i_stat,i_all,'hamovr','import_gaussians')

 endif

 if (iproc.eq.0) write(*,'(1x,a)')'done.'

END SUBROUTINE import_gaussians

subroutine input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nat,natsc,norb,norbp,n1,n2,n3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
     rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
     nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
     atomnames,ntypes,iatype,iasctype,pkernel,nzatom,nelpsp,psppar,npspcode,ixc,&
     ppsi,ppsit,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar, &
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
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
  integer, intent(in) :: nspin
  real(kind=8), dimension(norb), intent(in) :: spinar
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
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(nprojel), intent(in) :: proj
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:), pointer :: ppsi,ppsit
  !real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: ppsi
  !********************Alexey***************************************************************
  !for shrink:
  integer ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  !for grow:
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  !for real space:
  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
!*****************************************************************************************

  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer, parameter :: ngx=31
  integer :: i,iorb,iorbsc,imatrsc,iorbst,imatrst,i_stat,i_all,ierr,info,jproc,jpst,norbeyou
  integer :: norbe,norbep,norbi,norbj,norbeme,ndim_hamovr,n_lp,norbi_max,norbsc
  integer :: ispin,norbu,norbd,iorbst2
  real(kind=8) :: hgridh,tt,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum
  logical, dimension(:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: norbsc_arr,ng
  integer, dimension(:,:), allocatable :: nl
  real(kind=8), dimension(:), allocatable :: work_lp,pot,evale,occupe,ones
  real(kind=8), dimension(:,:), allocatable :: xp,occupat,hamovr,psi,hpsi
  real(kind=8), dimension(:,:,:), allocatable :: psiw,psiat

  !Calculate no. up and dw orbitals for spin-polarized starting guess
  norbu=0
  norbd=0
  do iorb=1,norb
     if(spinar(iorb)>0.0d0) norbu=norbu+1
     if(spinar(iorb)<0.0d0) norbd=norbd+1
  end do

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
       nat,rxyz,norbe,norbep,norbsc,occupe,occupat,ngx,xp,psiat,ng,nl,&
       nvctr_c,nvctr_f,n1,n2,n3,hgrid,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       nseg_c,nseg_f,keyg,keyv,iatype,ntypes,iasctype,natsc,psi,eks,scorb)

!!$  !!plot the initial LCAO wavefunctions
!!$  !do i=2*iproc+1,2*iproc+2
!!$  !   iounit=15+3*(i-1)
!!$  !   print *,'iounit',iounit,'-',iounit+2
!!$  !   call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$  !        rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc), &
!!$          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
!!$          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!!$  !end do


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
  allocate(ones(norbe),stat=i_stat)
  call memocc(i_stat,product(shape(ones))*kind(ones),'ones','input_wf_diag')
  ones(:)=1.0d0

  call sumrho(parallel,iproc,nproc,norbe,norbep,n1,n2,n3,hgrid,occupe,  & 
          nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot,&
          (2*n1+31)*(2*n2+31)*nscatterarr(iproc,1),nscatterarr,1,ones, &
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,1)

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
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,1,ones,&
       ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
       ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  i_all=-product(shape(ones))*kind(ones)
  deallocate(ones,stat=i_stat)
  call memocc(i_stat,i_all,'ones','input_wf_diag')

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
     !here the timing is related to the input guess part
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
!!$        !open(33+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(33+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(33+2*(i-1))
!!$        !open(34+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(34+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jjorb+(jiorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(34+2*(i-1))

        if (iproc.eq.0) then
        do iorb=1,norbi
           if (nspin==1) then
              if (iorb+iorbst-1 == norb) then
                 write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                      'evale(',iorb+iorbst-1,')=',evale(iorb),&
                      ' <- Last eigenvalue for input wavefncts'
              else
                 write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                      'evale(',iorb+iorbst-1,')=',evale(iorb)
              end if
           else
              if (norbu==norbd) then
                 if (iorb+iorbst-1 == norbu) then
                    write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb),&
                         ' <- Last eigenvalue for input wavefncts of both spins'
                 else
                    write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb)
                 end if
              else
                 if (iorb+iorbst-1 == norbu) then
                    write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb),&
                         ' <- Last eigenvalue for input wavefncts of spin up'
                 else if (iorb+iorbst-1 == norbd) then
                    write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb),&
                         ' <- Last eigenvalue for input wavefncts of spin down'
                 else
                    write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb)
                 end if
              end if
           end if
        enddo
        endif
        do iorb=iorbst,min(norbi+iorbst-1,norb)
           eval(iorb)=evale(iorb-iorbst+1)
        enddo
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do

!     write(*,*) "NORBO",norb,norbe,norbu,norbd,norbu+norbd,norbp
     ! Copy eigenvalues from NM to spin-polarized channels
     if(nspin>1) then
!        do iorb=1,norbu
!!           write(*,*) 'jorb:',iorb,norbi_max,norbe,norb,product(shape(eval)),product(shape(evale))
!           evale(iorb)=eval(iorb)
!        end do
        do iorb=1,norbd
!           write(*,*) 'korb:',iorb,iorb+norbu,norbi_max,norbe,norb,product(shape(evale))
           eval(iorb+norbu)=eval(iorb)
        end do
     end if

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
!!     iorbst=1
!!     imatrst=1
!!     do i=1,natsc
!!        norbi=norbsc_arr(i)
!!        call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
!!             hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst),nvctrp)
!!        iorbst=iorbst+norbi
!!        imatrst=imatrst+norbi**2
!!     end do
!!     norbi=norbsc_arr(natsc+1)
!!     norbj=norb-norbsc
!!     call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
!!          hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst),nvctrp)

     !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
     do ispin=1,nspin
        iorbst=1
        iorbst2=0
        if(ispin==2) iorbst2=norbu
        imatrst=1
        norbsc=0
        do i=1,natsc
           norbi=norbsc_arr(i)
           norbsc=norbsc+norbi
           call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst2+iorbst),nvctrp)
           iorbst=iorbst+norbi
           imatrst=imatrst+norbi**2
        end do
        norbi=norbsc_arr(natsc+1)
        if(ispin==1) norbj=norbu-norbsc
        if(ispin==2) norbj=norbd-norbsc
!        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
        !        norbj=norb-norbsc
        if(norbj>0) then
           call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst2+iorbst),nvctrp)
        end if
     end do


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
     allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','input_wf_diag')

     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(ppsit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psiw,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')

     !allocate the direct wavefunction
     allocate(ppsi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','input_wf_diag')

     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,ppsi)
     call timing(iproc,'Un-TransSwitch','OF')

     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw','input_wf_diag')

  else !serial case

     write(*,'(1x,a)',advance='no')'Overlap Matrix...'

     allocate(hamovr(ndim_hamovr,2),stat=i_stat)
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
     allocate(evale(max(norbi_max,norb)),stat=i_stat)
     call memocc(i_stat,product(shape(evale))*kind(evale),'evale','input_wf_diag')

     write(*,'(1x,a)')'Linear Algebra...'
     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DSYGV(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
             norbi,evale,work_lp,n_lp,info)

        if (info.ne.0) write(*,*) 'DSYGV ERROR',info,i,natsc+1
        if (iproc.eq.0) then
           do iorb=1,norbi
              if (nspin==1) then
                 if (iorb+iorbst-1 == norb) then
                    write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb),&
                         ' <- Last eigenvalue for input wavefncts'
                 else
                    write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb)
                 end if
              else
                 if (norbu==norbd) then
                    if (iorb+iorbst-1 == norbu) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of both spins'
                    else
                       write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb)
                    end if
                 else
                    if (iorb+iorbst-1 == norbu) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of spin up'
                    else if (iorb+iorbst-1 == norbd) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of spin down'
                    else
                       write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb)
                    end if
                 end if
              end if
           enddo
        endif
        do iorb=iorbst,min(norbi+iorbst-1,norb)
!           write(*,*) 'iorb:',iorb,iorb-iorbst+1,norbi_max,norbe,norb,product(shape(evale))
           eval(iorb)=evale(iorb-iorbst+1)
        enddo
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
!     write(*,*) "NORBO",norb,norbe,norbu,norbd,norbu+norbd,norbp
     ! Copy eigenvalues from NM to spin-polarized channels
     if(nspin>1) then
!        do iorb=1,norbu
!!           write(*,*) 'jorb:',iorb,norbi_max,norbe,norb,product(shape(eval)),product(shape(evale))
!           evale(iorb)=eval(iorb)
!        end do
        do iorb=1,norbd
!           write(*,*) 'korb:',iorb,iorb+norbu,norbi_max,norbe,norb,product(shape(evale))
           eval(iorb+norbu)=eval(iorb)
        end do
     end if
     !
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
     do ispin=1,nspin
        iorbst=1
        iorbst2=0
        if(ispin==2) iorbst2=norbu
        imatrst=1
        norbsc=0
        do i=1,natsc
           norbi=norbsc_arr(i)
           norbsc=norbsc+norbi
           call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsi(1,iorbst2+iorbst),nvctrp)
           iorbst=iorbst+norbi
           imatrst=imatrst+norbi**2
        end do
        norbi=norbsc_arr(natsc+1)
        if(ispin==1) norbj=norbu-norbsc
        if(ispin==2) norbj=norbd-norbsc
!        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
        !        norbj=norb-norbsc
        if(norbj>0) then
           call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsi(1,iorbst2+iorbst),nvctrp)
        end if
     end do

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

END SUBROUTINE input_wf_diag

subroutine diisstp(parallel,norb,norbp,nproc,iproc,  & 
                   ads,ids,mids,idsx,nvctrp,psit,psidst,hpsidst)
! diis subroutine:
! calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
! using  the previous iteration points phidst and the associated error 
! vectors (preconditione gradients) hpsidst
  implicit real(kind=8) (a-h,o-z)
  include 'mpif.h'
  logical parallel
  dimension psit(nvctrp,norbp*nproc),ads(idsx+1,idsx+1,3), &
  psidst(nvctrp,norbp*nproc,idsx),hpsidst(nvctrp,norbp*nproc,idsx)
  allocatable :: ipiv(:),rds(:)

  allocate(ipiv(idsx+1),stat=i_stat)
  call memocc(i_stat,product(shape(ipiv))*kind(ipiv),'ipiv','diisstp')
  allocate(rds(idsx+1),stat=i_stat)
  call memocc(i_stat,product(shape(rds))*kind(rds),'rds','diisstp')

! set up DIIS matrix (upper triangle)
  if (ids.gt.idsx) then
! shift left up matrix
     do i=1,idsx-1
        do j=1,i
           ads(j,i,1)=ads(j+1,i+1,1)
        end do
     end do
  end if

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
  do j=1,min(idsx,ids)
     ads(j,min(idsx,ids)+1,2)=1.d0
     rds(j)=0.d0
     do i=j,min(idsx,ids)
        ads(j,i,2)=ads(j,i,1)
     end do
  end do
  ads(min(idsx,ids)+1,min(idsx,ids)+1,2)=0.d0
  rds(min(idsx,ids)+1)=1.d0

!  write(6,*) 'DIIS matrix'
!  do i=1,min(idsx,ids)+1
!  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(idsx,ids)+1),rds(i)
!  enddo
  if (ids.gt.1) then
           ! solve linear system:(LAPACK)
           call DSYSV('U',min(idsx,ids)+1,1,ads(1,1,2),idsx+1,  & 
                ipiv,rds,idsx+1,ads(1,1,3),(idsx+1)**2,info)
           if (info.ne.0) then
              print*, 'DGESV',info
              stop
           end if
  else
     rds(1)=1.d0
  endif
  if (iproc.eq.0) then 
     !write(*,*) 'DIIS weights'
     write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j),j=1,min(idsx,ids)+1)
  endif

! new guess
  do iorb=1,norb
     call razero(nvctrp,psit(1,iorb))

     jst=max(1,ids-idsx+1)
     jj=0
     do j=jst,ids
        jj=jj+1
        mj=mod(j-1,idsx)+1
        do k=1,nvctrp
           psit(k,iorb)=psit(k,iorb)+rds(jj)*(psidst(k,iorb,mj)-hpsidst(k,iorb,mj))
        end do
     end do
  end do

  i_all=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv,stat=i_stat)
  call memocc(i_stat,i_all,'ipiv','diisstp')
  i_all=-product(shape(rds))*kind(rds)
  deallocate(rds,stat=i_stat)
  call memocc(i_stat,i_all,'rds','diisstp')

END SUBROUTINE diisstp

subroutine read_input_variables(iproc,in)

  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  integer :: ierror

  ! Read the input variables.
  open(unit=1,file='input.dat',status='old')
  read(1,*,iostat=ierror) in%ncount_cluster_x
  read(1,*,iostat=ierror) in%frac_fluct
  read(1,*,iostat=ierror) in%randdis
  read(1,*,iostat=ierror) in%betax
  read(1,*,iostat=ierror) in%hgrid
  read(1,*,iostat=ierror) in%crmult
  read(1,*,iostat=ierror) in%frmult
  read(1,*,iostat=ierror) in%cpmult
  read(1,*,iostat=ierror) in%fpmult
  if (in%fpmult.gt.in%frmult) write(*,*) ' NONSENSE: fpmult > frmult'
  read(1,*,iostat=ierror) in%ixc
  read(1,*,iostat=ierror) in%ncharge,in%elecfield
  read(1,*,iostat=ierror) in%gnrm_cv
  read(1,*,iostat=ierror) in%itermax
  read(1,*,iostat=ierror) in%ncong
  read(1,*,iostat=ierror) in%idsx
  read(1,*,iostat=ierror) in%calc_tail
  read(1,*,iostat=ierror) in%rbuf
  read(1,*,iostat=ierror) in%ncongt
  read(1,*,iostat=ierror) in%nspin,in%mpol

  if (ierror/=0) then
     write(*,'(1x,a)') 'Error while reading the file "input.dat"'
     stop
  end if

  close(1,iostat=ierror)

  !these values are hard-coded for the moment but they can be entered in the input file
  !in case of need also other variables can be entered without any changements
  in%output_grid=.false. 
  in%inputPsiId=0
  in%output_wf=.false. 

  if (iproc == 0) then
     write(*,'(1x,a,i0)') 'Max. number of wavefnctn optim ',in%ncount_cluster_x
     write(*,'(1x,a,1pe10.2)') 'Convergence criterion for forces: fraction of noise ',&
          in%frac_fluct
     write(*,'(1x,a,1pe10.2)') 'Random displacement amplitude ',in%randdis
     write(*,'(1x,a,1pe10.2)') 'Steepest descent step ',in%betax

     write(*,'(1x,a)')&
          '------------------------------------------------------------------- Input Parameters'
     write(*,'(1x,a)')&
          '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
     write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
          'Grid spacing=',in%hgrid,   '|  Coarse Wfs.=',in%crmult,'| Wavefns Conv.=',in%gnrm_cv,&
          '| Calculate=',in%calc_tail
     write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,f4.1)')&
          '       XC id=',in%ixc,     '|    Fine Wfs.=',in%frmult,'| Max. N. Iter.=',in%itermax,&
          '| Extension=',in%rbuf
     write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,i4)')&
          'total charge=',in%ncharge, '| Coarse Proj.=',in%cpmult,'| CG Prec.Steps=',in%ncong,&
          '|  CG Steps=',in%ncongt
     write(*,'(1x,a,1pe7.1,1x,a,0pf5.2,1x,a,i8)')&
          ' elec. field=',in%elecfield,'|   Fine Proj.=',in%fpmult,'| DIIS Hist. N.=',in%idsx
  end if

end subroutine read_input_variables


end module libBigDFT
