
 subroutine call_cluster(parallel,nproc,iproc,atoms,rxyz,energy,fxyz,&
     psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
  use module_types
  implicit none
  logical, intent(in) :: parallel
  integer, intent(in) :: iproc,nproc
  type(input_variables),intent(inout) :: in
  type(wavefunctions_descriptors), intent(inout) :: wfd
  type(atoms_data), intent(inout) :: atoms
  integer, intent(inout) :: infocode,n1,n2,n3,norbp,norb
  real(kind=8), intent(out) :: energy
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz,rxyz_old
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:,:), pointer :: psi
  !local variables
  integer :: i_stat,i_all,ierr,inputPsiId_orig
  !temporary interface
  interface
     subroutine cluster(parallel,nproc,iproc,atoms,rxyz,energy,fxyz,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
       use module_types
       use module_interfaces
       use Poisson_Solver
       implicit none
       logical, intent(in) :: parallel
       integer, intent(in) :: nproc,iproc
       integer, intent(inout) :: n1,n2,n3,norbp,norb
       integer, intent(out) :: infocode
       type(input_variables), intent(in) :: in
       type(wavefunctions_descriptors), intent(inout) :: wfd
       type(atoms_data), intent(inout) :: atoms
       real(kind=8), intent(out) :: energy
       real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
       real(kind=8), dimension(3,atoms%nat), intent(out) :: rxyz_old,fxyz
       real(kind=8), dimension(:), pointer :: eval
       real(kind=8), dimension(:,:), pointer :: psi
     end subroutine cluster
  end interface

  inputPsiId_orig=in%inputPsiId

  loop_cluster: do

     if (in%inputPsiId == 0 .and. associated(psi)) then
        i_all=-product(shape(psi))*kind(psi)
        deallocate(psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi','call_cluster')
        i_all=-product(shape(eval))*kind(eval)
        deallocate(eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval','call_cluster')

        call deallocate_wfd(wfd,'call_cluster')
     end if

     call cluster(parallel,nproc,iproc,atoms,rxyz,energy,fxyz,&
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
        call memocc(i_stat,i_all,'psi','call_cluster')
        i_all=-product(shape(eval))*kind(eval)
        deallocate(eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval','call_cluster')

        call deallocate_wfd(wfd,'call_cluster')
        !finalize memory counting (there are still the positions and the forces allocated)
        call memocc(0,0,'count','stop')

        if (parallel) call MPI_FINALIZE(ierr)

        stop
     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  inputPsiId_orig=in%inputPsiId

end subroutine call_cluster

subroutine cluster(parallel,nproc,iproc,atoms,rxyz,energy,fxyz,&
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

  use module_types
  use module_interfaces
  use Poisson_Solver
  implicit none
  include 'mpif.h'
  !implicit real(kind=8) (a-h,o-z)
  logical, intent(in) :: parallel
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: n1,n2,n3,norbp,norb
  integer, intent(out) :: infocode
  type(input_variables), intent(in) :: in
  type(wavefunctions_descriptors), intent(inout) :: wfd
  type(atoms_data), intent(inout) :: atoms
  real(kind=8), intent(out) :: energy
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(kind=8), dimension(3,atoms%nat), intent(out) :: rxyz_old,fxyz
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:,:), pointer :: psi
  !local variables
  character(len=1) :: datacode,geocode
  logical :: calc_tail,output_wf,output_grid
  integer :: ixc,ncharge,ncong,idsx,ncongt,nspin,mpol,inputPsiId,itermax
  integer :: nelec,norbu,norbd,ndegree_ip,nvctrp,mids,iorb,iounit
  integer :: n1_old,n2_old,n3_old,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n3d,n3p,n3pi,i3xcsh,i3s
  integer :: ncount0,ncount1,ncount_rate,ncount_max,iunit,n1i,n2i,n3i
  integer :: i1,i2,i3,ind,iat,ierror,i_all,i_stat,iter,ierr,i03,i04,jproc,ispin
  real :: tcpu0,tcpu1
  real(kind=8) :: hgrid,crmult,frmult,cpmult,fpmult,elecfield,gnrm_cv,rbuf,hx,hy,hz,hxh,hyh,hzh
  real(kind=8) :: hgridh,peakmem,alat1,alat2,alat3,accurex,gnrm_check,hgrid_old,energy_old
  real(kind=8) :: eion,epot_sum,ekin_sum,eproj_sum,ehart,eexcu,vexcu,alpha,gnrm,evsum
  real(kind=8) :: scprsum,energybs,tt,tel,eexcu_fake,vexcu_fake,ehart_fake,energy_min
  type(wavefunctions_descriptors) :: wfd_old
  type(convolutions_bounds) :: bounds
  type(nonlocal_psp_descriptors) :: nlpspd
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(kind=8), dimension(:), allocatable :: occup,spinar,spinar_foo,derproj,rho
  real(kind=8), dimension(:,:), allocatable :: radii_cf,gxyz
  ! Charge density/potential,ionic potential, pkernel
  real(kind=8), dimension(:), allocatable :: pot_ion
  real(kind=8), dimension(:,:,:,:), allocatable :: rhopot,pot
  real(kind=8), dimension(:), pointer :: pkernel
  real(kind=8), dimension(:), pointer :: eval_old !should be removed from copy_old_wavefunctions
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(kind=8), dimension(:,:), pointer :: hpsi,psit,psi_old
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: proj
  ! arrays for DIIS convergence accelerator
  real(kind=8), dimension(:,:,:), pointer :: ads,psidst,hpsidst

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
  mpol=in%mpol

  inputPsiId=in%inputPsiId
  output_grid=in%output_grid
  output_wf=in%output_wf

  hx=in%hgrid
  hy=in%hgrid
  hz=in%hgrid


  !define the geometry code: hard coded to free BC for the moment
  geocode='F'

  if (iproc.eq.0) then
     write(*,'(1x,a,1x,i0)') &
       '===================== BigDFT Wavefunction Optimization =============== inputPsiId=',inputPsiId
     call print_input_parameters(in)
  end if
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
     call copy_old_wavefunctions(iproc,nproc,norb,norbp,hgrid,n1,n2,n3,eval,wfd,psi,&
          hgrid_old,n1_old,n2_old,n3_old,eval_old,wfd_old,psi_old)
  end if

  !datacodes for the poisson solver, depending on the implementation
  datacode='D'
  if (.not. parallel) datacode='G'

  if(nspin<1.or.nspin>2) nspin=1
  if(nspin==1) mpol=0

  ! grid spacing (same in x,y and z direction)
  !modify the grid spacings such that it allows for inhomogeneous values of hgrid
  hgridh=.5d0*hgrid

  if (iproc==0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------------ System Properties'
  end if

  allocate(radii_cf(atoms%ntypes,2),stat=i_stat)
  call memocc(i_stat,product(shape(radii_cf))*kind(radii_cf),'radii_cf','cluster')


  call read_system_variables(iproc,nproc,in,atoms,radii_cf,nelec,norb,norbu,norbd,norbp,iunit)

  allocate(occup(norb),stat=i_stat)
  call memocc(i_stat,product(shape(occup))*kind(occup),'occup','cluster')
  allocate(spinar(norb),stat=i_stat)
  call memocc(i_stat,product(shape(spinar))*kind(spinar),'spinar','cluster')

  ! Occupation numbers
  call input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinar)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,geocode,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,&
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)

  hxh=0.5d0*hx
  hyh=0.5d0*hy
  hzh=0.5d0*hz

  !memory estimation
  if (iproc==0) then
     !for the moment leave the atoms data as they are due to the module dependence
     call MemoryEstimator(geocode,nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hx,hy,hz,&
          atoms%nat,atoms%ntypes,atoms%iatype,rxyz,radii_cf,crmult,frmult,norb,&
          atoms%atomnames,.false.,nspin,peakmem) 
  end if

  !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
  ndegree_ip=16 !default value to be put to 16 and update references for test
  call createKernel(geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,iproc,nproc,pkernel)

  ! Create wavefunctions descriptors and allocate them
  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,nproc,geocode,n1,n2,n3,output_grid,hx,hy,hz,&
       atoms,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,wfd,&
       nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
  call timing(iproc,'CrtDescriptors','OF')

  ! Calculate all projectors
  call timing(iproc,'CrtProjectors ','ON')
  call createProjectorsArrays(geocode,iproc,n1,n2,n3,rxyz,atoms,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  call timing(iproc,'CrtProjectors ','OF')

  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(nscatterarr(0:nproc-1,4),stat=i_stat)
  call memocc(i_stat,product(shape(nscatterarr))*kind(nscatterarr),'nscatterarr','cluster')
  !allocate array for the communications of the potential
  allocate(ngatherarr(0:nproc-1,2),stat=i_stat)
  call memocc(i_stat,product(shape(ngatherarr))*kind(ngatherarr),'ngatherarr','cluster')

  !create the descriptors for the density and the potential
  call createDensPotDescriptors(iproc,nproc,geocode,datacode,n1i,n2i,n3i,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion(n1i*n2i*n3pi),stat=i_stat)
     call memocc(i_stat,product(shape(pot_ion))*kind(pot_ion),'pot_ion','cluster')
  else
     allocate(pot_ion(1),stat=i_stat)
     call memocc(i_stat,product(shape(pot_ion))*kind(pot_ion),'pot_ion','cluster')
  end if

  !do not pass the atoms data structure as arguments for dependence problem in compiling
  call createIonicPotential(geocode,iproc,nproc,atoms%nat,atoms%ntypes,atoms%iatype,&
       atoms%psppar,atoms%nelpsp,rxyz,hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s+i3xcsh,&
       n1i,n2i,n3i,pkernel,pot_ion,eion)

  !Allocate Charge density, Potential in real space
  if (n3d >0) then
     allocate(rhopot(n1i,n2i,n3d,nspin),stat=i_stat)
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
     call import_gaussians(parallel,iproc,nproc,atoms,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          norb,norbp,occup,n1,n2,n3,nvctrp,hgrid,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj,  &
          pkernel,ixc,psi,psit,hpsi,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar)

  else if (inputPsiId == 0) then

     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(geocode,iproc,nproc,atoms,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          norb,norbp,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj, &
          pkernel,ixc,psi,psit,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar)

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
        allocate(hpsi(nvctrp,norb),stat=i_stat)
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
        allocate(psi(nvctrp,norb),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     end if

     if (iproc.eq.0) then
        write(*,'(1x,a)')&
             '-------------------------------------------------------------- Wavefunctions Restart'
     end if
     call reformatmywaves(iproc,norb,norbp,atoms%nat,hgrid_old,n1_old,n2_old,n3_old,&
          rxyz_old,wfd_old,psi_old,hgrid,n1,n2,n3,rxyz,wfd,psi)
     eval=eval_old

     call deallocate_wfd(wfd_old,'cluster')

     i_all=-product(shape(psi_old))*kind(psi_old)
     deallocate(psi_old,stat=i_stat)
     call memocc(i_stat,i_all,'psi_old','cluster')
     i_all=-product(shape(eval_old))*kind(eval_old)
     deallocate(eval_old,stat=i_stat)
     call memocc(i_stat,i_all,'eval_old','cluster')

     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,&
          wfd%nvctr_c,wfd%nvctr_f,nvctrp,nspin,psi,hpsi,psit)

  else if (inputPsiId == 2 ) then 
     !restart from previously calculated wavefunctions, on disk

     !allocate principal wavefunction
     if (parallel) then
        !allocated in the transposed way such as 
        !it can also be used as a work array for transposition
        allocate(psi(nvctrp,norbp*nproc),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     else
        allocate(psi(nvctrp,norb),stat=i_stat)
        call memocc(i_stat,product(shape(psi))*kind(psi),'psi','cluster')
     end if

     if (iproc.eq.0) then
        write(*,'(1x,a)')&
             '---------------------------------------------------- Reading Wavefunctions from disk'
     end if


     call readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,atoms%nat,rxyz,wfd,psi,eval)

     !initialise control value for gnrm in the case of a restart
     gnrm_check=0.d0

     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,&
          wfd%nvctr_c,wfd%nvctr_f,nvctrp,nspin,psi,hpsi,psit)

  else

     if (iproc == 0) then
        write(*,'(1x,a)')'The supported values of inputPsiId are integers from -1 to 2'
        write(*,'(1x,a,i0)')'                        while we found',inputPsiId
     end if
     stop

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
  call memocc(i_stat,i_all,'nzatom','cluster')
  i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
  deallocate(atoms%iasctype,stat=i_stat)
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
  !minimum value of the energy during the minimisation procedure
  energy_min=1.d10
  !set the infocode to the value it would have in the case of no convergence
  infocode=1

  !end of the initialization part
  call timing(iproc,'INIT','PR')

  ! loop for wavefunction minimization
  wfn_loop: do iter=1,itermax
     if (idsx.gt.0) mids=mod(iter-1,idsx)+1
     if (iproc.eq.0) then 
        write(*,'(1x,a,i0)')&
             '---------------------------------------------------------------------------- iter= ',&
             iter
     endif

     !control whether the minimisation iterations ended
     if (gnrm.le.gnrm_cv .or. iter.eq.itermax) call timing(iproc,'WFN_OPT','PR')

     ! Potential from electronic charge density
     call sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
          wfd,psi,rhopot,(2*n1+31)*(2*n2+31)*n3d,nscatterarr,nspin,spinar,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)

     !     ixc=11  ! PBE functional
     !     ixc=1   ! LDA functional

     call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin) !add NSPIN

     call HamiltonianApplication(parallel,datacode,iproc,nproc,atoms,hgrid,&
          norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          wfd,bounds,nlpspd,proj,ngatherarr,n3p,&
          rhopot(1,1,1+i3xcsh,1),psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,spinar)
     energybs=ekin_sum+epot_sum+eproj_sum
     energy_old=energy
     energy=energybs-ehart+eexcu-vexcu+eion
     energy_min=min(energy_min,energy)

     !check for convergence or whether max. numb. of iterations exceeded
     if (gnrm.le.gnrm_cv .or. iter.eq.itermax) then 
        if (iproc.eq.0) then 
           write(*,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write(*,'(1x,a)') &
                '--------------------------------------------------- End of Wavefunction Optimisation'
           write(*,'(1x,a,3(1x,1pe18.11))') &
                'final  ekin,  epot,  eproj ',ekin_sum,epot_sum,eproj_sum
           write(*,'(1x,a,3(1x,1pe18.11))') &
                'final ehart, eexcu,  vexcu ',ehart,eexcu,vexcu
           write(*,'(1x,a,i6,2x,1pe19.12,1x,1pe9.2)') &
                'FINAL iter,total energy,gnrm',iter,energy,gnrm
           !write(61,*)hgrid,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
           if (energy > energy_min) write(*,'(1x,a,1pe9.2)')&
                'WARNING: Found an energy value lower than the FINAL energy, delta:',energy-energy_min
        end if
        infocode=0
        exit wfn_loop 
     endif

     call hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,bounds%kb,&
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
        if ((gnrm > 4.d0 .and. norbu/=norbd) .or. (norbu==norbd .and. gnrm > 10.d0)) then
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
  call last_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,&
       nspin,psi,hpsi,psit,occup,evsum,eval)

  if (abs(evsum-energybs).gt.1.d-8 .and. iproc==0) write(*,'(1x,a,2(1x,1pe20.13))')&
       'Difference:evsum,energybs',evsum,energybs

  !plot the converged wavefunctions in the different orbitals
!!$  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
!!$     iounit=iorb
!!$     call plot_wf(iounit,n1,n2,n3,hgrid,wfd%nseg_c,wfd%nvctr_c,wfd%keyg,wfd%keyv,wfd%nseg_f,wfd%nvctr_f,  & 
!!$          rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(1,iorb-iproc*norbp),&
!!$          bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,&
!!$          bounds%gb%ibxxyy_f,bounds%ibyyzz_r,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!!$  end do

  !  write all the wavefunctions into files
  if (output_wf) then
     call  writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,atoms%nat,rxyz,wfd,psi,eval)
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
       wfd,psi,rho,(2*n1+31)*(2*n2+31)*n3p,nscatterarr,1,spinar_foo,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)

  i_all=-product(shape(spinar_foo))*kind(spinar_foo)
  deallocate(spinar_foo,stat=i_stat)
  call memocc(i_stat,i_all,'spinar_foo','cluster')

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion','cluster')

  if (iproc.eq.0 .and. output_grid) then

     open(unit=22,file='density.pot',status='unknown')
     write(22,*)'density'
     write(22,*) 2*n1+1,2*n2+1,2*n3+1
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     write(22,*)'xyz'

     if (nproc > 1) then
        !allocate full density in pot_ion directory
        allocate(pot_ion((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
        call memocc(i_stat,product(shape(pot))*kind(pot),'pot','hamiltonianapplication')

        call MPI_ALLGATHERV(pot_ion,(2*n1+31)*(2*n2+31)*n3p,&
             MPI_DOUBLE_PRECISION,rho,ngatherarr(0,1),&
             ngatherarr(0,2),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        do i3=1,2*n3+1
           do i2=1,2*n2+1
              do i1=1,2*n1+1
                 ind=i1+14+(i2+13)*(2*n1+31)+(i3+13)*(2*n1+31)*(2*n2+31)
                 write(22,*)rho(ind)
              end do
           end do
        end do
        i_all=-product(shape(pot_ion))*kind(pot_ion)
        deallocate(pot_ion,stat=i_stat)
        call memocc(i_stat,i_all,'pot_ion','cluster')
     else
        do i3=1,2*n3+1
           do i2=1,2*n2+1
              do i1=1,2*n1+1
                 ind=i1+14+(i2+13)*(2*n1+31)+(i3+13)*(2*n1+31)*(2*n2+31)
                 write(22,*)rho(ind)
              end do
           end do
        end do
     end if
     close(22)
  endif

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

  allocate(gxyz(3,atoms%nat),stat=i_stat)
  call memocc(i_stat,product(shape(gxyz))*kind(gxyz),'gxyz','cluster')

  call timing(iproc,'Forces        ','ON')
  ! calculate local part of the forces gxyz
  call local_forces(iproc,nproc,atoms,rxyz,hgrid,n1,n2,n3,n3p,i3s+i3xcsh,rho,pot,gxyz)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho','cluster')
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot','cluster')

  allocate(derproj(3*nlpspd%nprojel),stat=i_stat)
  call memocc(i_stat,product(shape(derproj))*kind(derproj),'derproj','cluster')

  if (iproc == 0) write(*,'(1x,a)',advance='no')'Calculate projectors derivatives...'

  !the calculation of the derivatives of the projectors has been decoupled
  !from the one of nonlocal forces, in this way forces can be calculated
  !diring the wavefunction minimization if needed
  call projectors_derivatives(geocode,iproc,atoms,n1,n2,n3,norb,nlpspd,proj,&
       rxyz,radii_cf,cpmult,fpmult,hx,hy,hz,derproj)

  if (iproc == 0) write(*,'(1x,a)',advance='no')'done, calculate nonlocal forces...'

  call nonlocal_forces(iproc,atoms,norb,norbp,occup,nlpspd,proj,derproj,wfd,psi,gxyz)

  if (iproc == 0) write(*,'(1x,a)')'done.'

  i_all=-product(shape(derproj))*kind(derproj)
  deallocate(derproj,stat=i_stat)
  call memocc(i_stat,i_all,'derproj','cluster')

  i_all=-product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c)
  deallocate(nlpspd%nboxp_c,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_c','cluster')
  i_all=-product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f)
  deallocate(nlpspd%nboxp_f,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_f','cluster')

  ! Add up all the force contributions
  if (parallel) then
     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     do iat=1,atoms%nat
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
  if (calc_tail .and. geocode == 'F') then
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

     !pass hx instead of hgrid since we are only in free BC
     call CalculateTailCorrection(iproc,nproc,atoms,n1,n2,n3,rbuf,norb,norbp,&
          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,eval,&
          pot,hx,rxyz,radii_cf,crmult,frmult,nspin,spinar,&
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

       i_all=-product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c)
       deallocate(nlpspd%nboxp_c,stat=i_stat)
       call memocc(i_stat,i_all,'nboxp_c','cluster')
       i_all=-product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f)
       deallocate(nlpspd%nboxp_f,stat=i_stat)
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

    call deallocate_bounds(bounds,'cluster')

    i_all=-product(shape(nlpspd%keyg_p))*kind(nlpspd%keyg_p)
    deallocate(nlpspd%keyg_p,stat=i_stat)
    call memocc(i_stat,i_all,'keyg_p','cluster')
    i_all=-product(shape(nlpspd%keyv_p))*kind(nlpspd%keyv_p)
    deallocate(nlpspd%keyv_p,stat=i_stat)
    call memocc(i_stat,i_all,'keyv_p','cluster')
    i_all=-product(shape(nlpspd%nvctr_p))*kind(nlpspd%nvctr_p)
    deallocate(nlpspd%nvctr_p,stat=i_stat)
    call memocc(i_stat,i_all,'nvctr_p','cluster')
    i_all=-product(shape(nlpspd%nseg_p))*kind(nlpspd%nseg_p)
    deallocate(nlpspd%nseg_p,stat=i_stat)
    call memocc(i_stat,i_all,'nseg_p','cluster')
    i_all=-product(shape(proj))*kind(proj)
    deallocate(proj,stat=i_stat)
    call memocc(i_stat,i_all,'proj','cluster')
    i_all=-product(shape(occup))*kind(occup)
    deallocate(occup,stat=i_stat)
    call memocc(i_stat,i_all,'occup','cluster')
    i_all=-product(shape(spinar))*kind(spinar)
    deallocate(spinar,stat=i_stat)
    call memocc(i_stat,i_all,'spinar','cluster')
    i_all=-product(shape(atoms%psppar))*kind(atoms%psppar)
    deallocate(atoms%psppar,stat=i_stat)
    call memocc(i_stat,i_all,'psppar','cluster')
    i_all=-product(shape(atoms%nelpsp))*kind(atoms%nelpsp)
    deallocate(atoms%nelpsp,stat=i_stat)
    call memocc(i_stat,i_all,'nelpsp','cluster')
    i_all=-product(shape(radii_cf))*kind(radii_cf)
    deallocate(radii_cf,stat=i_stat)
    call memocc(i_stat,i_all,'radii_cf','cluster')
    i_all=-product(shape(atoms%npspcode))*kind(atoms%npspcode)
    deallocate(atoms%npspcode,stat=i_stat)
    call memocc(i_stat,i_all,'npspcode','cluster')

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


