 subroutine call_cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
     psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(input_variables),intent(inout) :: in
  type(wavefunctions_descriptors), intent(inout) :: wfd
  type(atoms_data), intent(inout) :: atoms
  integer, intent(inout) :: infocode,n1,n2,n3,norbp,norb
  real(kind=8), intent(out) :: energy
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz,rxyz_old
  real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:), pointer :: psi
  !local variables
  character(len=*), parameter :: subname='call_cluster'
  integer :: i_stat,i_all,ierr,inputPsiId_orig
  !temporary interface
  interface
     subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
       use module_types
       implicit none
       integer, intent(in) :: nproc,iproc
       integer, intent(inout) :: n1,n2,n3,norbp,norb
       integer, intent(out) :: infocode
       type(input_variables), intent(in) :: in
       type(wavefunctions_descriptors), intent(inout) :: wfd
       type(atoms_data), intent(inout) :: atoms
       real(kind=8), intent(out) :: energy
       real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz,rxyz_old
       real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz
       real(kind=8), dimension(:), pointer :: eval
       real(kind=8), dimension(:), pointer :: psi
     end subroutine cluster
  end interface

  inputPsiId_orig=in%inputPsiId

  loop_cluster: do

     if (in%inputPsiId == 0 .and. associated(psi)) then
        i_all=-product(shape(psi))*kind(psi)
        deallocate(psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(eval))*kind(eval)
        deallocate(eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(wfd,'call_cluster')
     end if

     call cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)

     if (in%inputPsiId==1 .and. infocode==2) then
        in%inputPsiId=0
     else if (in%inputPsiId==1 .and. infocode==1) then
        in%inputPsiId=0
     else if (in%inputPsiId == 0 .and. infocode==3) then
        if (iproc.eq.0) then
           write(*,'(1x,a)')'Convergence error, cannot proceed.'
           write(*,'(1x,a)')' writing positions in file posout_999.xyz then exiting'

           call wtposout(999,energy,rxyz,atoms)

        end if

        i_all=-product(shape(psi))*kind(psi)
        deallocate(psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(eval))*kind(eval)
        deallocate(eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(wfd,'call_cluster')
        !finalize memory counting (there are still the positions and the forces allocated)
        call memocc(0,0,'count','stop')

        if (nproc > 1) call MPI_FINALIZE(ierr)

        stop
     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  in%inputPsiId=inputPsiId_orig

end subroutine call_cluster

subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
     psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
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
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  implicit none
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: n1,n2,n3,norbp,norb
  integer, intent(out) :: infocode
  type(input_variables), intent(in) :: in
  type(wavefunctions_descriptors), intent(inout) :: wfd
  type(atoms_data), intent(inout) :: atoms
  real(kind=8), intent(out) :: energy
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz,rxyz_old
  real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:), pointer :: psi
  !local variables
  character(len=*), parameter :: subname='cluster'
  character(len=1) :: geocode
  character(len=10) :: orbname
  logical :: calc_tail,switchSD
  integer :: ixc,ncharge,ncong,idsx,ncongt,nspin,mpol,itermax,idsx_actual,nvirte,nvirtep,nvirt
  integer :: nelec,norbu,norbd,ndegree_ip,nvctrp,mids,iorb,iounit,ids,idiistol,j
  integer :: n1_old,n2_old,n3_old,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n3d,n3p,n3pi,i3xcsh,i3s
  integer :: ncount0,ncount1,ncount_rate,ncount_max,iunit,n1i,n2i,n3i,nl1,nl2,nl3
  integer :: i1,i2,i3,ind,iat,ierror,i_all,i_stat,iter,ierr,i03,i04,jproc,ispin,nspinor,nplot
  real :: tcpu0,tcpu1
  real(kind=8) :: hgrid,crmult,frmult,cpmult,fpmult,elecfield,gnrm_cv,rbuf,hx,hy,hz,hxh,hyh,hzh
  real(kind=8) :: peakmem,alat1,alat2,alat3,gnrm_check,hgrid_old,energy_old,sumz
  real(kind=8) :: eion,epot_sum,ekin_sum,eproj_sum,ehart,eexcu,vexcu,alpha,gnrm,evsum,sumx,sumy
  real(kind=8) :: scprsum,energybs,tt,tel,eexcu_fake,vexcu_fake,ehart_fake,energy_min,psoffset
  real(kind=8) :: factor,rhon,rhos,ttsum,hx_old,hy_old,hz_old
  type(wavefunctions_descriptors) :: wfd_old
  type(convolutions_bounds) :: bounds
  type(nonlocal_psp_descriptors) :: nlpspd
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:), allocatable :: occup,spinsgn,spinsgn_foo,rho
  real(kind=8), dimension(:,:), allocatable :: radii_cf,gxyz,fion!,derproj
  ! Charge density/potential,ionic potential, pkernel
  real(kind=8), dimension(:), allocatable :: pot_ion
  real(kind=8), dimension(:,:,:,:), allocatable :: rhopot,pot,rho_diag
  real(kind=8), dimension(:,:,:), allocatable :: m_norm
  real(kind=8), dimension(:), pointer :: pkernel
  real(kind=8), dimension(:), pointer :: eval_old !should be removed from copy_old_wavefunctions
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:), pointer :: hpsi,psit,psi_old,psivirt,psidst,hpsidst
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj
  ! arrays for DIIS convergence accelerator
  real(kind=8), dimension(:,:,:), pointer :: ads
  ! tmp debug array
  real(kind=8), dimension(:,:), allocatable :: tmred
  

  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be 
  !reassigned inside the structure

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
  if(nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  mpol=in%mpol

  nvirt=in%nvirt
  nplot=in%nplot

  hx=in%hgrid
  hy=in%hgrid
  hz=in%hgrid

  geocode=in%geocode

  alat1=in%alat1
  alat2=in%alat2
  alat3=in%alat3

  if (iproc.eq.0) then
     write(*,'(1x,a,1x,i0)') &
       '===================== BigDFT Wavefunction Optimization =============== inputPsiId=',&
       in%inputPsiId
     call print_input_parameters(in)
  end if
  if (nproc > 1) then
     call timing(iproc,'parallel     ','IN')
  else
     call timing(iproc,'             ','IN')
  end if
  call cpu_time(tcpu0)
  call system_clock(ncount0,ncount_rate,ncount_max)



  ! We save the variables that defined the previous psi if
  ! restartOnPsi is .true.
  if (in%inputPsiId == 1) then
     !regeneerate grid spacings
     if (geocode == 'P') then
        call correct_grid(in%alat1,hx,n1)
        call correct_grid(in%alat2,hy,n2)
        call correct_grid(in%alat3,hz,n3)
     end if
     call copy_old_wavefunctions(iproc,nproc,norb,norbp,nspinor,hx,hy,hz,n1,n2,n3,wfd,psi,&
          hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,wfd_old,psi_old)
  end if

!!$  !datacodes for the poisson solver, depending on the implementation
!!$  datacode='D'
!!$  !leaving always datacode to D
!!$  !if () datacode='G'

  if(nspin/=1.and.nspin/=2.and.nspin/=4) nspin=1
  if(nspin==1) mpol=0

  ! grid spacing (same in x,y and z direction)

  if (iproc==0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------------ System Properties'
  end if

  allocate(radii_cf(atoms%ntypes,2+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call read_system_variables(iproc,nproc,in,atoms,radii_cf,nelec,norb,norbu,norbd,norbp,iunit)

  allocate(occup(norb+ndebug),stat=i_stat)
  call memocc(i_stat,occup,'occup',subname)
  allocate(spinsgn(norb+ndebug),stat=i_stat)
  call memocc(i_stat,spinsgn,'spinsgn',subname)

  ! Occupation numbers
  call input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinsgn)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,geocode,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,&
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)

  hxh=0.5d0*hx
  hyh=0.5d0*hy
  hzh=0.5d0*hz

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=16 !default value to be put to 16 and update references for test
  call createKernel(geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,iproc,nproc,pkernel)

  ! Create wavefunctions descriptors and allocate them
  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,nproc,geocode,n1,n2,n3,in%output_grid,hx,hy,hz,&
       atoms,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,wfd,&
       nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,nspinor)
  call timing(iproc,'CrtDescriptors','OF')

  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call timing(iproc,'CrtProjectors ','ON')
  call createProjectorsArrays(geocode,iproc,n1,n2,n3,rxyz,atoms,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  call timing(iproc,'CrtProjectors ','OF')

  !memory estimation
  if (iproc==0) then
     call MemoryEstimator(geocode,nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hx,hy,hz,&
          atoms%nat,atoms%ntypes,atoms%iatype,rxyz,radii_cf,crmult,frmult,norb,nlpspd%nprojel,&
          atoms%atomnames,.false.,nspin,peakmem) 
  end if

  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  !allocate array for the communications of the potential
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)

  !create the descriptors for the density and the potential
  call createDensPotDescriptors(iproc,nproc,geocode,'D',n1i,n2i,n3i,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion(n1i*n2i*n3pi+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  else
     allocate(pot_ion(1+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  end if

  !here calculate the ionic energy and forces accordingly
  allocate(fion(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fion,'fion',subname)

  call IonicEnergyandForces(geocode,iproc,nproc,atoms,hxh,hyh,hzh,alat1,alat2,alat3,rxyz,eion,fion,&
       psoffset,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)

  !can pass the atoms data structure as argument
  call createIonicPotential(geocode,iproc,nproc,atoms%nat,atoms%ntypes,atoms%iatype,&
       atoms%psppar,atoms%nelpsp,rxyz,hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s+i3xcsh,&
       n1i,n2i,n3i,pkernel,pot_ion,eion,psoffset)

  !Allocate Charge density, Potential in real space
  if (n3d >0) then
     allocate(rhopot(n1i,n2i,n3d,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  else
     allocate(rhopot(1,1,1,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopot,'rhopot',subname)
  end if

  if (in%inputPsiId /= 1) then
     allocate(eval(norb+ndebug),stat=i_stat)
     call memocc(i_stat,eval,'eval',subname)
  end if

  ! INPUT WAVEFUNCTIONS, added also random input guess
  if (in%inputPsiId == -2) then

     if (iproc.eq.0) then
        write(*,'(1x,a)')&
             '------------------------------------------------ Random wavefunctions initialization'
     end if

     !random initialisation of the wavefunctions
     allocate(psi(nvctrp*nspinor*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     psi=0.0d0
     ttsum=0.0d0
     do iorb=1,norbp*nproc*max(1,nspinor)
        if(mod(iorb-1,nspinor)==0) then
           do i1=1,nvctrp
              do j=0,iproc-1
                 call random_number(tt)
              end do
              call random_number(tt)
              psi(i1+nvctrp*(iorb-1))=real(tt,kind=8)*0.01d0
              ttsum=ttsum+psi(i1+nvctrp*(iorb-1))
              do j=iproc+1,nproc
                 call random_number(tt)
              end do
           end do
        end if
     end do
     !write(*,'(a,30f10.4)') 'Rand Check',ttsum,(sum(psi(:,iorb)),iorb=1,norbp*nproc*nspinor)
 
     eval(:)=-0.5d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else if (in%inputPsiId == -1) then

     !import gaussians form CP2K (data in files gaubasis.dat and gaucoeff.dat)
     !and calculate eigenvalues
     call import_gaussians(geocode,iproc,nproc,cpmult,fpmult,radii_cf,atoms,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          norb,norbp,occup,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj, &
          pkernel,ixc,psi,psit,hpsi,eval,nscatterarr,ngatherarr,nspin,spinsgn)

  else if (in%inputPsiId == 0) then

     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(geocode,iproc,nproc,cpmult,fpmult,radii_cf,atoms,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          norb,norbp,nvirte,nvirtep,nvirt,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,&
          wfd,bounds,nlpspd,proj,pkernel,ixc,psi,hpsi,psit,psivirt,eval,&
          nscatterarr,ngatherarr,nspin,spinsgn)
  
  else if (in%inputPsiId == 1 ) then 
     !these parts should be reworked for the non-collinear spin case

     !restart from previously calculated wavefunctions, in memory

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition

     allocate(psi(nvctrp*norbp*nproc*nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     if (iproc.eq.0) then
        write(*,'(1x,a)')&
             '-------------------------------------------------------------- Wavefunctions Restart'
     end if

     call reformatmywaves(iproc,norb*nspinor,norbp*nspinor,atoms%nat,hx_old,hy_old,hz_old,&
          n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
!!$     eval=eval_old

     call deallocate_wfd(wfd_old,'cluster')

     i_all=-product(shape(psi_old))*kind(psi_old)
     deallocate(psi_old,stat=i_stat)
     call memocc(i_stat,i_all,'psi_old',subname)

!!$     i_all=-product(shape(eval_old))*kind(eval_old)
!!$     deallocate(eval_old,stat=i_stat)
!!$     call memocc(i_stat,i_all,'eval_old',subname)

     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else if (in%inputPsiId == 2 ) then 
     !restart from previously calculated wavefunctions, on disk

     !allocate principal wavefunction
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(psi(nvctrp*norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)

     if (iproc.eq.0) then
        write(*,'(1x,a)')&
             '---------------------------------------------------- Reading Wavefunctions from disk'
     end if

     call readmywaves(iproc,norb,norbp,n1,n2,n3,hx,hy,hz,atoms%nat,rxyz,wfd,psi,eval)

     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit)

  else

     if (iproc == 0) then
        write(*,'(1x,a)')'The supported values of inputPsiId are integers from -2 to 2'
        write(*,'(1x,a,i0)')'                        while we found',in%inputPsiId
     end if
     stop

  end if

  !this rearrange the value of norbu, to be changed
  if(nspinor==4) then
     norbu=norb
     norbd=0
  end if

  !save the new atomic positions in the rxyz_old array
  do iat=1,atoms%nat
     rxyz_old(1,iat)=rxyz(1,iat)
     rxyz_old(2,iat)=rxyz(2,iat)
     rxyz_old(3,iat)=rxyz(3,iat)
  enddo

  !no need of using nzatom array, semicores useful only for the input guess
  i_all=-product(shape(atoms%nzatom))*kind(atoms%nzatom)
  deallocate(atoms%nzatom,stat=i_stat)
  call memocc(i_stat,i_all,'nzatom',subname)
  i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
  deallocate(atoms%iasctype,stat=i_stat)
  call memocc(i_stat,i_all,'iasctype',subname)

  ! allocate arrays necessary for DIIS convergence acceleration
  if (idsx.gt.0) then
     allocate(psidst(nvctrp*nspinor*norbp*nproc*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,psidst,'psidst',subname)
     allocate(hpsidst(nvctrp*nspinor*norbp*nproc*idsx+ndebug),stat=i_stat)
     call memocc(i_stat,hpsidst,'hpsidst',subname)
     allocate(ads(idsx+1,idsx+1,3+ndebug),stat=i_stat)
     call memocc(i_stat,ads,'ads',subname)
     call razero(3*(idsx+1)**2,ads)
  endif

  alpha=2.d0
  energy=1.d10
  gnrm=1.d10
  ekin_sum=0.d0 
  epot_sum=0.d0 
  eproj_sum=0.d0
  !minimum value of the energy during the minimisation procedure
  energy_min=1.d10
  !set the infocode to the value it would have in the case of no convergence
  infocode=1
  !logical control variable for switch DIIS-SD
  switchSD=.false.
  !local variable for the diis history
  idsx_actual=idsx

  ids=0
  idiistol=0

  !end of the initialization part
  call timing(iproc,'INIT','PR')

  ! loop for wavefunction minimization
  wfn_loop: do iter=1,itermax
     if (idsx.gt.0) then
        mids=mod(ids,idsx)+1
        ids=ids+1
     end if
     if (iproc.eq.0) then 
        write(*,'(1x,a,i0)')&
             '---------------------------------------------------------------------------- iter= ',&
             iter
     endif

     !control whether the minimisation iterations ended
     if (gnrm <= gnrm_cv .or. iter == itermax) call timing(iproc,'WFN_OPT','PR')

     ! Potential from electronic charge density
     call sumrho(geocode,iproc,nproc,norb,norbp,ixc,n1,n2,n3,hxh,hyh,hzh,occup,  & 
          wfd,psi,rhopot,n1i*n2i*n3d,nscatterarr,nspin,nspinor,spinsgn,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)

     !this section should be inserted into sumrho
     if(nspinor==4) then
        !        call calc_moments(iproc,nproc,norb,norbp,nvctrp*nproc,nspinor,psi)
        allocate(tmred(nspin+1,2+ndebug),stat=i_stat)
        call memocc(i_stat,tmred,'tmred',subname)
        tmred=0.0d0
        do ispin=1,nspin
           tmred(ispin,1)=sum(rhopot(:,:,:,ispin))
        end do
        tmred(nspin+1,1)=sum(rhopot(:,:,:,:))
        if (nproc>1) then
           call MPI_ALLREDUCE(tmred(:,1),tmred(:,2),nspin+1,MPI_DOUBLE_PRECISION,&
                MPI_SUM,MPI_COMM_WORLD,ierr)
           tt=sqrt(tmred(2,2)**2+tmred(3,2)**2+tmred(4,2)**2)
           if(iproc==0.and.tt>0.0d0) write(*,'(a,5f10.4)') '  Magnetic density orientation:', &
                (tmred(ispin,2)/tmred(1,2),ispin=2,nspin)
        else
           tt=sqrt(tmred(2,1)**2+tmred(3,1)**2+tmred(4,1)**2)
           if(iproc==0.and.tt>0.0d0) write(*,'(a,5f10.4)') '  Magnetic density orientation:',&
                (tmred(ispin,1)/tmred(1,1),ispin=2,nspin)
        end if
        i_all=-product(shape(tmred))*kind(tmred)
        deallocate(tmred,stat=i_stat)
        call memocc(i_stat,i_all,'tmred',subname)
     end if
     
     if(nspinor==4) then
        !this wrapper can be inserted inside the poisson solver or in sumrho
         call PSolverNC(geocode,'D',iproc,nproc,n1i,n2i,n3i,n3d,ixc,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin)
     else
              
        call PSolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,ixc,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin)
        
     end if

     call HamiltonianApplication(geocode,iproc,nproc,atoms,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
          norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          wfd,bounds,nlpspd,proj,ngatherarr,n1i*n2i*n3p,&
          rhopot(1,1,1+i3xcsh,1),psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,nspinor,spinsgn)

     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-ehart+eexcu-vexcu+eion
     energy_min=min(energy_min,energy)

     !check for convergence or whether max. numb. of iterations exceeded
     if (gnrm <= gnrm_cv .or. iter == itermax) then 
        if (iproc.eq.0) then 
           write(*,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write(*,'(1x,a)') &
                '--------------------------------------------------- End of Wavefunction Optimisation'
           write(*,'(1x,a,3(1x,1pe18.11))') &
                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
           write(*,'(1x,a,3(1x,1pe18.11))') &
                'final ehart, eexcu,  vexcu ',ehart,eexcu,vexcu
           write(*,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') &
                'FINAL iter,total energy,gnrm',iter,energy,gnrm
           !write(61,*)hgrid,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
           if (energy > energy_min) write(*,'(1x,a,1pe9.2)')&
                'WARNING: Found an energy value lower than the FINAL energy, delta:',energy-energy_min
        end if
        if (gnrm <= gnrm_cv) infocode=0
        exit wfn_loop 
     endif

     call hpsitopsi(geocode,ids,iproc,nproc,norb,norbp,occup,hx,hy,hz,n1,n2,n3,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,bounds%kb,&
          eval,ncong,mids,idsx_actual,ads,energy,energy_old,alpha,gnrm,scprsum,&
          psi,psit,hpsi,psidst,hpsidst,nspin,nspinor,spinsgn)

     tt=(energybs-scprsum)/scprsum
     if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
          (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
        write(*,'(1x,a,1pe9.2,2(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
     endif
     if (iproc.eq.0) then
        write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             ekin_sum,epot_sum,eproj_sum
        write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
        write(*,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
     endif

     if (in%inputPsiId == 0) then
        if ((gnrm > 4.d0 .and. norbu/=norbd) .or. (norbu==norbd .and. gnrm > 10.d0)) then
           if (iproc == 0) then
              write(*,'(1x,a)')&
                   'Error: the norm of the residue is too large also with input wavefunctions.'
           end if
           infocode=3
           call deallocate_before_exiting
           return
        end if
     else if (in%inputPsiId == 1) then
        if (gnrm > 1.d0) then
           if (iproc == 0) then
              write(*,'(1x,a)')&
                   'The norm of the residue is too large, need to recalculate input wavefunctions'
           end if
           infocode=2
           if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
              if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              call deallocate_before_exiting
              return
           end if
        end if
     end if
 
     !add switch between DIIS and SD
     !this section should be inserted into hpsitopsi
     if (energy == energy_min .and. .not. switchSD) idiistol=0
     if (energy > energy_min .and. idsx >0 .and. .not. switchSD) then
        idiistol=idiistol+1
     end if
     if (idiistol > idsx .and. .not. switchSD) then
        !the energy has not decreasing for too much steps, switching to SD for next steps
        if (iproc ==0) write(*,'(1x,a,1pe9.2,a)')&
                'WARNING: The energy value is growing (delta=',energy-energy_min,') switch to SD'
        switchSD=.true.
        i_all=-product(shape(psidst))*kind(psidst)
        deallocate(psidst,stat=i_stat)
        call memocc(i_stat,i_all,'psidst',subname)
        i_all=-product(shape(hpsidst))*kind(hpsidst)
        deallocate(hpsidst,stat=i_stat)
        call memocc(i_stat,i_all,'hpsidst',subname)
        i_all=-product(shape(ads))*kind(ads)
        deallocate(ads,stat=i_stat)
        call memocc(i_stat,i_all,'ads',subname)
        idsx_actual=0
        idiistol=0
     end if

     if ((energy == energy_min) .and. switchSD) then
        idiistol=idiistol+1
     end if
     if (idiistol > idsx .and. switchSD) then
        !restore the original DIIS
        if (iproc ==0) write(*,'(1x,a,1pe9.2)')&
                'WARNING: The energy value is now decreasing again, coming back to DIIS'
        switchSD=.false.
        idsx_actual=idsx
        ids=0
        idiistol=0

        allocate(psidst(nvctrp*nspinor*norbp*nproc*idsx+ndebug),stat=i_stat)
        call memocc(i_stat,psidst,'psidst',subname)
        allocate(hpsidst(nvctrp*nspinor*norbp*nproc*idsx+ndebug),stat=i_stat)
        call memocc(i_stat,hpsidst,'hpsidst',subname)
        allocate(ads(idsx+1,idsx+1,3+ndebug),stat=i_stat)
        call memocc(i_stat,ads,'ads',subname)
        call razero(3*(idsx+1)**2,ads)
     end if


  end do wfn_loop
  if (iter == itermax .and. iproc == 0 ) &
       write(*,'(1x,a)')'No convergence within the allowed number of minimization steps'

  if (idsx_actual > 0) then
     i_all=-product(shape(psidst))*kind(psidst)
     deallocate(psidst,stat=i_stat)
     call memocc(i_stat,i_all,'psidst',subname)
     i_all=-product(shape(hpsidst))*kind(hpsidst)
     deallocate(hpsidst,stat=i_stat)
     call memocc(i_stat,i_all,'hpsidst',subname)
     i_all=-product(shape(ads))*kind(ads)
     deallocate(ads,stat=i_stat)
     call memocc(i_stat,i_all,'ads',subname)
  end if

  ! transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
  call last_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,nspin,psi,hpsi,psit,&
       occup,evsum,eval)

  if (abs(evsum-energybs).gt.1.d-8 .and. iproc==0) write(*,'(1x,a,2(1x,1pe20.13))')&
       'Difference:evsum,energybs',evsum,energybs
 
  if (nvirt > 0 .and. in%inputPsiId == 0) then
     call davidson(geocode,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i,atoms,&
          cpmult,fpmult,radii_cf,&
          norb,norbu,norbp,nvirte,nvirtep,nvirt,gnrm_cv,nplot,n1,n2,n3,nvctrp,&
          hx,hy,hz,rxyz,rhopot,occup,i3xcsh,n3p,itermax,wfd,bounds,nlpspd,proj,  &
          pkernel,ixc,psi,psivirt,eval,ncong,nscatterarr,ngatherarr)
  end if

  !  write all the wavefunctions into files
  if (in%output_wf) then
     !add flag for writing waves in the gaussian basis form
     if (in%inputPsiId >= 10) then
        !extract the gaussian basis from the pseudopotential

        !extract the gaussian basis from the wavefunctions
        
     else
        call  writemywaves(iproc,norb,norbp,n1,n2,n3,hx,hy,hz,atoms%nat,rxyz,wfd,psi,eval)
        write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
     end if
  end if


  !------------------------------------------------------------------------
  ! here we start the calculation of the forces
  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------------- Forces Calculation'
  end if

  ! Selfconsistent potential is saved in rhopot, 
  ! new arrays rho,pot for calculation of forces ground state electronic density

  allocate(spinsgn_foo(norb+ndebug),stat=i_stat)
  call memocc(i_stat,spinsgn_foo,'spinsgn_foo',subname)
  spinsgn_foo(:)=1.0d0
  ! Potential from electronic charge density

  !manipulate scatter array for avoiding the GGA shift
  do jproc=0,nproc-1
     !n3d=n3p
     nscatterarr(jproc,1)=nscatterarr(jproc,2)
     !i3xcsh=0
     nscatterarr(jproc,4)=0
  end do

  !here there are the spinor which must be taken into account
  if(nproc>1) then
     if (n3p>0) then
        allocate(rho(n1i*n2i*n3p+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     else
        allocate(rho(1+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     end if
  else
     if (n3p>0) then
        allocate(rho(n1i*n2i*n3p*nspinor+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     else
        allocate(rho(1+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     end if
  end if

  call sumrho(geocode,iproc,nproc,norb,norbp,0,n1,n2,n3,hxh,hyh,hzh,occup,  & 
       wfd,psi,rho,n1i*n2i*n3p,nscatterarr,1,nspinor,spinsgn_foo,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)

  i_all=-product(shape(spinsgn_foo))*kind(spinsgn_foo)
  deallocate(spinsgn_foo,stat=i_stat)
  call memocc(i_stat,i_all,'spinsgn_foo',subname)

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)

  if (in%output_grid) then
     if (iproc == 0) then
        open(unit=22,file='density.pot',status='unknown')
        write(22,*)'normalised density'
        write(22,*) 2*n1+2,2*n2+2,2*n3+2
        write(22,*) alat1,' 0. ',alat2
        write(22,*) ' 0. ',' 0. ',alat3
        write(22,*)'xyz   periodic'
     end if

     !conditions for periodicity in the three directions
     !value of the buffer in the x and z direction
     if (geocode /= 'F') then
        nl1=1
        nl3=1
     else
        nl1=14
        nl3=14
     end if
     !value of the buffer in the y direction
     if (geocode == 'P') then
        nl2=1
     else
        nl2=14
     end if


     if (nproc > 1) then
        !allocate full density in pot_ion directory
        allocate(pot_ion(n1i*n2i*n3i+ndebug),stat=i_stat)
        call memocc(i_stat,pot_ion,'pot_ion',subname)

        call MPI_ALLGATHERV(rho,n1i*n2i*n3p,&
             MPI_DOUBLE_PRECISION,pot_ion,ngatherarr(0,1),&
             ngatherarr(0,2),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        if (iproc == 0) then
           do i3=0,2*n3+1
              do i2=0,2*n2+1
                 do i1=0,2*n1+1
                    ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                    write(22,*)pot_ion(ind)/real(nelec,dp)
                 end do
              end do
           end do
        end if

        i_all=-product(shape(pot_ion))*kind(pot_ion)
        deallocate(pot_ion,stat=i_stat)
        call memocc(i_stat,i_all,'pot_ion',subname)
     else
        do i3=0,2*n3+1
           do i2=0,2*n2+1
              do i1=0,2*n1+1
                 ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                 write(22,*)rho(ind)/real(nelec,dp)
              end do
           end do
        end do
     end if
     if (iproc == 0 )close(22)
  endif

  if (n3p>0) then
     allocate(pot(n1i,n2i,n3p,1+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
  else
     allocate(pot(1,1,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
  end if

  call DCOPY(n1i*n2i*n3p,rho,1,pot,1) 
  call PSolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
       pot,pkernel,pot,ehart_fake,eexcu_fake,vexcu_fake,0.d0,.false.,1)
  !here nspin=1 since ixc=0

  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'pkernel',subname)

  allocate(gxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gxyz,'gxyz',subname)

  call timing(iproc,'Forces        ','ON')
  ! calculate local part of the forces gxyz
  call local_forces(geocode,iproc,nproc,atoms,rxyz,hxh,hyh,hzh,&
       n1,n2,n3,n3p,i3s+i3xcsh,n1i,n2i,n3i,rho,pot,gxyz)

!!$  sumx=0.d0
!!$  sumy=0.d0
!!$  sumz=0.d0
!!$  write(*,'(1x,a,19x,a)') 'Final values of the Local Forces for each atom'
!!$  do iat=1,atoms%nat
!!$     write(*,'(1x,a,i5,i5,1x,a6,3(1x,1pe12.5))') &
!!$          'L',iproc,iat,trim(atoms%atomnames(atoms%iatype(iat))),(gxyz(j,iat),j=1,3)
!!$     sumx=sumx+gxyz(1,iat)
!!$     sumy=sumy+gxyz(2,iat)
!!$     sumz=sumz+gxyz(3,iat)
!!$  enddo
!!$  write(*,'(1x,a)')'the sum of the forces is'
!!$  write(*,'(1x,a,3x,i5,1pe16.8)')'x direction(L)',iproc,sumx
!!$  write(*,'(1x,a,3x,i5,1pe16.8)')'y direction(L)',iproc,sumy
!!$  write(*,'(1x,a,3x,i5,1pe16.8)')'z direction(L)',iproc,sumz

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)

!!$  allocate(derproj(nlpspd%nprojel,3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,derproj,'derproj',subname)

  if (iproc == 0) write(*,'(1x,a)',advance='no')'Calculate projectors derivatives...'

  !the calculation of the derivatives of the projectors has been decoupled
  !from the one of nonlocal forces, in this way forces can be calculated
  !during the wavefunction minimization if needed
!!$  call projectors_derivatives(geocode,iproc,atoms,n1,n2,n3,norb,nlpspd,proj,&
!!$       rxyz,radii_cf,cpmult,fpmult,hx,hy,hz,derproj)

  if (iproc == 0) write(*,'(1x,a)',advance='no')'done, calculate nonlocal forces...'

  call nonlocal_forces(geocode,iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,atoms,rxyz,radii_cf,&
     norb,norbp,nspinor,occup,nlpspd,proj,wfd,psi,gxyz,calc_tail) !refill projectors for tails
!!$
!!$  do i1=1,3
!!$     call fill_projectors(geocode,iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,atoms,rxyz,radii_cf,&
!!$          nlpspd,derproj(1,i1),i1)
!!$  end do
!!$
!!$  call nonlocal_forcesold(iproc,atoms,norb,norbp,occup,nlpspd,proj,derproj,wfd,psi,gxyz,nspinor)

  if (iproc == 0) write(*,'(1x,a)')'done.'

!!$  sumx=0.d0
!!$  sumy=0.d0
!!$  sumz=0.d0
!!$  write(*,'(1x,a,19x,a)') 'Final values of the NonLocal+LOCAL Forces for each atom'
!!$  do iat=1,atoms%nat
!!$     write(*,'(1x,a,i5,i5,1x,a6,3(1x,1pe12.5))') &
!!$          'NL',iproc,iat,trim(atoms%atomnames(atoms%iatype(iat))),(gxyz(j,iat),j=1,3)
!!$     sumx=sumx+gxyz(1,iat)
!!$     sumy=sumy+gxyz(2,iat)
!!$     sumz=sumz+gxyz(3,iat)
!!$  enddo
!!$  write(*,'(1x,a)')'the sum of the forces is'
!!$  write(*,'(1x,a,3x,i5,1pe16.8)')'x direction(NL)',iproc,sumx
!!$  write(*,'(1x,a,3x,i5,1pe16.8)')'y direction(NL)',iproc,sumy
!!$  write(*,'(1x,a,3x,i5,1pe16.8)')'z direction(NL)',iproc,sumz


!!$  i_all=-product(shape(derproj))*kind(derproj)
!!$  deallocate(derproj,stat=i_stat)
!!$  call memocc(i_stat,i_all,'derproj',subname)

  ! Add up all the force contributions
  if (nproc > 1) then
     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     do iat=1,atoms%nat
        fxyz(1,iat)=gxyz(1,iat)
        fxyz(2,iat)=gxyz(2,iat)
        fxyz(3,iat)=gxyz(3,iat)
     enddo
  end if

  !add to the forces the ionic contribution 
  fxyz(:,:)=fxyz(:,:)+fion(:,:)

  i_all=-product(shape(fion))*kind(fion)
  deallocate(fion,stat=i_stat)
  call memocc(i_stat,i_all,'fion',subname)
  i_all=-product(shape(gxyz))*kind(gxyz)
  deallocate(gxyz,stat=i_stat)
  call memocc(i_stat,i_all,'gxyz',subname)

  call timing(iproc,'Forces        ','OF')

  !------------------------------------------------------------------------
  if (calc_tail .and. geocode == 'F') then
     call timing(iproc,'Tail          ','ON')
     !    Calculate energy correction due to finite size effects
     !    ---reformat potential
     allocate(pot(n1i,n2i,n3i,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)

     if (nproc > 1) then
        call MPI_ALLGATHERV(rhopot(1,1,1+i3xcsh,1),n1i*n2i*n3p,&
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
           call MPI_ALLGATHERV(rhopot(1,1,i03,i04),n1i*n2i*n3p,&
                MPI_DOUBLE_PRECISION,pot(1,1,1,2),ngatherarr(0,1),ngatherarr(0,2), & 
                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        end if
     else
        do ispin=1,nspin
           !here one could have not allocated pot and: call move_alloc(rhopot,pot) 
           !(but it is a Fortran 95/2003 spec)
           do i3=1,n3i
              do i2=1,n2i
                 do i1=1,n1i
                    pot(i1,i2,i3,ispin)=rhopot(i1,i2,i3,ispin)
                 enddo
              enddo
           enddo
        end do
     end if
     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
     deallocate(nscatterarr,stat=i_stat)
     call memocc(i_stat,i_all,'nscatterarr',subname)
     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr',subname)
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot',subname)

     !pass hx instead of hgrid since we are only in free BC
     call CalculateTailCorrection(iproc,nproc,atoms,n1,n2,n3,rbuf,norb,norbp,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,eval,&
          pot,hx,rxyz,radii_cf,crmult,frmult,cpmult,fpmult,nspin,spinsgn,&
          proj,psi,occup,in%output_grid,ekin_sum,epot_sum,eproj_sum)

     i_all=-product(shape(pot))*kind(pot)
     deallocate(pot,stat=i_stat)
     call memocc(i_stat,i_all,'pot',subname)

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
        write(*,'(1x,a,1x,1pe24.17)')&
             'Total energy with tail correction',energy
     endif

     call timing(iproc,'Tail          ','OF')
  else
     !    No tail calculation
     if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     i_all=-product(shape(rhopot))*kind(rhopot)
     deallocate(rhopot,stat=i_stat)
     call memocc(i_stat,i_all,'rhopot',subname)
     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
     deallocate(nscatterarr,stat=i_stat)
     call memocc(i_stat,i_all,'nscatterarr',subname)
     i_all=-product(shape(ngatherarr))*kind(ngatherarr)
     deallocate(ngatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'ngatherarr',subname)
  endif
  ! --- End if of tail calculation

  call deallocate_before_exiting

contains

  !routine which deallocate the pointers and the arrays before exiting 
  subroutine deallocate_before_exiting

    !when this condition is verified we are in the middle of the SCF cycle
    if (infocode /=0 .and. infocode /=1) then

       if (idsx_actual > 0) then
          i_all=-product(shape(psidst))*kind(psidst)
          deallocate(psidst,stat=i_stat)
          call memocc(i_stat,i_all,'psidst',subname)
          i_all=-product(shape(hpsidst))*kind(hpsidst)
          deallocate(hpsidst,stat=i_stat)
          call memocc(i_stat,i_all,'hpsidst',subname)
          i_all=-product(shape(ads))*kind(ads)
          deallocate(ads,stat=i_stat)
          call memocc(i_stat,i_all,'ads',subname)
       end if

       if (nproc > 1) then
          i_all=-product(shape(psit))*kind(psit)
          deallocate(psit,stat=i_stat)
          call memocc(i_stat,i_all,'psit',subname)
       end if

       i_all=-product(shape(hpsi))*kind(hpsi)
       deallocate(hpsi,stat=i_stat)
       call memocc(i_stat,i_all,'hpsi',subname)

       i_all=-product(shape(pot_ion))*kind(pot_ion)
       deallocate(pot_ion,stat=i_stat)
       call memocc(i_stat,i_all,'pot_ion',subname)

       i_all=-product(shape(pkernel))*kind(pkernel)
       deallocate(pkernel,stat=i_stat)
       call memocc(i_stat,i_all,'pkernel',subname)

       ! calc_tail false
       i_all=-product(shape(rhopot))*kind(rhopot)
       deallocate(rhopot,stat=i_stat)
       call memocc(i_stat,i_all,'rhopot',subname)
       i_all=-product(shape(nscatterarr))*kind(nscatterarr)
       deallocate(nscatterarr,stat=i_stat)
       call memocc(i_stat,i_all,'nscatterarr',subname)
       i_all=-product(shape(ngatherarr))*kind(ngatherarr)
       deallocate(ngatherarr,stat=i_stat)
       call memocc(i_stat,i_all,'ngatherarr',subname)

       i_all=-product(shape(fion))*kind(fion)
       deallocate(fion,stat=i_stat)
       call memocc(i_stat,i_all,'fion',subname)

    end if
    !deallocate wavefunction for virtual orbitals
    if (in%nvirt > 0) then
       i_all=-product(shape(psivirt))*kind(psivirt)
       deallocate(psivirt)
       call memocc(i_stat,i_all,'psivirt','cluster')
    end if

    if (geocode == 'F') then
       call deallocate_bounds(bounds,'cluster')
    end if

    i_all=-product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c)
    deallocate(nlpspd%nboxp_c,stat=i_stat)
    call memocc(i_stat,i_all,'nboxp_c',subname)
    i_all=-product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f)
    deallocate(nlpspd%nboxp_f,stat=i_stat)
    call memocc(i_stat,i_all,'nboxp_f',subname)
    i_all=-product(shape(nlpspd%keyg_p))*kind(nlpspd%keyg_p)
    deallocate(nlpspd%keyg_p,stat=i_stat)
    call memocc(i_stat,i_all,'keyg_p',subname)
    i_all=-product(shape(nlpspd%keyv_p))*kind(nlpspd%keyv_p)
    deallocate(nlpspd%keyv_p,stat=i_stat)
    call memocc(i_stat,i_all,'keyv_p',subname)
    i_all=-product(shape(nlpspd%nvctr_p))*kind(nlpspd%nvctr_p)
    deallocate(nlpspd%nvctr_p,stat=i_stat)
    call memocc(i_stat,i_all,'nvctr_p',subname)
    i_all=-product(shape(nlpspd%nseg_p))*kind(nlpspd%nseg_p)
    deallocate(nlpspd%nseg_p,stat=i_stat)
    call memocc(i_stat,i_all,'nseg_p',subname)

    i_all=-product(shape(proj))*kind(proj)
    deallocate(proj,stat=i_stat)
    call memocc(i_stat,i_all,'proj',subname)

    i_all=-product(shape(occup))*kind(occup)
    deallocate(occup,stat=i_stat)
    call memocc(i_stat,i_all,'occup',subname)
    i_all=-product(shape(spinsgn))*kind(spinsgn)
    deallocate(spinsgn,stat=i_stat)
    call memocc(i_stat,i_all,'spinsgn',subname)
    i_all=-product(shape(atoms%psppar))*kind(atoms%psppar)
    deallocate(atoms%psppar,stat=i_stat)
    call memocc(i_stat,i_all,'psppar',subname)
    i_all=-product(shape(atoms%nelpsp))*kind(atoms%nelpsp)
    deallocate(atoms%nelpsp,stat=i_stat)
    call memocc(i_stat,i_all,'nelpsp',subname)
    i_all=-product(shape(radii_cf))*kind(radii_cf)
    deallocate(radii_cf,stat=i_stat)
    call memocc(i_stat,i_all,'radii_cf',subname)
    i_all=-product(shape(atoms%npspcode))*kind(atoms%npspcode)
    deallocate(atoms%npspcode,stat=i_stat)
    call memocc(i_stat,i_all,'npspcode',subname)

    !end of wavefunction minimisation
    call timing(iproc,'LAST','PR')
    call timing(iproc,'              ','RE')
    call cpu_time(tcpu1)
    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)
    if (iproc == 0) &
         write(*,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time for root process ', iproc,tel,tcpu1-tcpu0

  end subroutine deallocate_before_exiting

END SUBROUTINE cluster


