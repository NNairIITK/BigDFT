!> @file
!! Main file for XANES calculation
!! @author Copyright (C) 2009-2011 BigDFT group (AM, ESRF)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>  Main program for XANES calculation (absorption calculation)
program abscalc_main

   use module_base
   use module_types
   use module_interfaces
   use m_ab6_symmetry
   !  use minimization, only: parameterminimization 

   implicit none
   character(len=*), parameter :: subname='abscalc_main'
   integer :: iproc,nproc,iat,j,i_stat,i_all,ierr,infocode
   real(gp) :: etot,sumx,sumy,sumz
   logical :: exist_list
   !input variables
   type(atoms_data) :: atoms
   type(input_variables) :: inputs
   type(restart_objects) :: rst
   character(len=50), dimension(:), allocatable :: arr_posinp
   character(len=60) :: radical
   ! atomic coordinates, forces
   real(gp), dimension(:,:), allocatable :: fxyz
   real(gp), dimension(:,:), pointer :: rxyz
   integer :: iconfig,nconfig,istat
   logical :: exists

   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

   ! Read a possible radical format argument.
   call get_command_argument(1, value = radical, status = istat)
   if (istat > 0) then
      write(radical, "(A)") "input"
   end if

   ! find out which input files will be used
   inquire(file="list_posinp",exist=exist_list)
   if (exist_list) then
      open(54,file="list_posinp")
      read(54,*) nconfig
      if (nconfig > 0) then 
         !allocation not referenced since memocc count not initialised
         allocate(arr_posinp(1:nconfig))

         do iconfig=1,nconfig
            read(54,*) arr_posinp(iconfig)
         enddo
      else
         !normal case
         nconfig=1
         allocate(arr_posinp(1:1))
         if (istat > 0) then
            arr_posinp(1)='posinp'
         else
            arr_posinp(1)=trim(radical)
         end if
      endif
      close(54)
   else
      nconfig=1
      allocate(arr_posinp(1:1))
      if (istat > 0) then
         arr_posinp(1)='posinp'
      else
         arr_posinp(1)=trim(radical)
      end if
   end if

   do iconfig=1,nconfig

      !Welcome screen
      if (iproc==0) call print_logo()

      ! Read all input files.
      !standard names
      call standard_inputfile_names(inputs,radical,nproc)
      call read_input_variables(iproc,trim(arr_posinp(iconfig)),inputs, atoms, rxyz)

      !Initialize memory counting
      !call memocc(0,iproc,'count','start')

      !Read absorption-calculation input variables
      !inquire for the needed file 
      !if not present, set default (no absorption calculation)

      inquire(file=trim(radical)//".abscalc",exist=exists)
      if (.not. exists) then
         if (iproc == 0) write(*,*) 'ERROR: need file input.abscalc for x-ray absorber treatment.'
         if(nproc/=0)   call MPI_FINALIZE(ierr)
         stop
      end if
      call abscalc_input_variables(iproc,trim(radical)//".abscalc",inputs)
      if( inputs%iat_absorber <1 .or. inputs%iat_absorber > atoms%nat) then
         if (iproc == 0) write(*,*)'ERROR: inputs%iat_absorber  must .ge. 1 and .le. number_of_atoms '
         if(nproc/=0)   call MPI_FINALIZE(ierr)
         stop
      endif


      !Allocations
      allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
      call memocc(i_stat,fxyz,'fxyz',subname)

      call init_restart_objects(iproc,inputs%iacceleration,atoms,rst,subname)

      call call_abscalc(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,rst,infocode)

      if (iproc == 0) then
         sumx=0.d0
         sumy=0.d0
         sumz=0.d0
         write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom'
         do iat=1,atoms%nat
            write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
               &   iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
            sumx=sumx+fxyz(1,iat)
            sumy=sumy+fxyz(2,iat)
            sumz=sumz+fxyz(3,iat)
         enddo
         if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
            write(*,'(1x,a)')'the sum of the forces is'
            write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
            write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
            write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
         end if
      endif

      !De-allocations
      call deallocate_abscalc_input(inputs, subname)
      call deallocate_atoms(atoms,subname) 
!      call deallocate_local_zone_descriptors(rst%Lzd, subname)

      call free_restart_objects(rst,subname)

      i_all=-product(shape(rxyz))*kind(rxyz)
      deallocate(rxyz,stat=i_stat)
      call memocc(i_stat,i_all,'rxyz',subname)
      i_all=-product(shape(fxyz))*kind(fxyz)
      deallocate(fxyz,stat=i_stat)
      call memocc(i_stat,i_all,'fxyz',subname)

      call free_input_variables(inputs)

      !finalize memory counting
      call memocc(0,0,'count','stop')

      !     call sg_end()

   enddo !loop over iconfig


   !No referenced by memocc!
   deallocate(arr_posinp)

   call MPI_FINALIZE(ierr)

END PROGRAM abscalc_main


!> Routines to use abscalc as a blackbox
subroutine call_abscalc(nproc,iproc,atoms,rxyz,in,energy,fxyz,rst,infocode)
   use module_base
   use module_types
   use module_interfaces
   implicit none
   !Arguments
   integer, intent(in) :: iproc,nproc
   type(input_variables),intent(inout) :: in
   type(atoms_data), intent(inout) :: atoms
   type(restart_objects), intent(inout) :: rst
   integer, intent(inout) :: infocode
   real(gp), intent(out) :: energy !< only iproc has the right value
   !Local variables
   real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
   real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
   !local variables
   character(len=*), parameter :: subname='call_bigdft'
   character(len=40) :: comment
   integer :: i_stat,i_all,ierr,inputPsiId_orig,icycle

   !temporary interface
   interface
      subroutine abscalc(nproc,iproc,atoms,rxyz,&
          psi,Lzd,orbs,hx_old,hy_old,hz_old,in,GPU,infocode)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: nproc,iproc
         integer, intent(out) :: infocode
         real(gp), intent(inout) :: hx_old,hy_old,hz_old
         type(input_variables), intent(in) :: in
       type(local_zone_descriptors), intent(inout) :: Lzd
         type(atoms_data), intent(inout) :: atoms
         type(orbitals_data), intent(inout) :: orbs
         type(GPU_pointers), intent(inout) :: GPU
         real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
         real(wp), dimension(:), pointer :: psi
      END SUBROUTINE abscalc 
   end interface

   !put a barrier for all the processes
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   !assign the verbosity of the output
   !the verbose variables is defined in module_base
   verbose=in%verbosity

   !Assign a value for energy to avoid compiler warning and to check the calculation
   energy = huge(1.d0)

   inputPsiId_orig=in%inputPsiId

   loop_cluster: do icycle=1,in%nrepmax

      if (in%inputPsiId == 0 .and. associated(rst%KSwfn%psi)) then
         i_all=-product(shape(rst%KSwfn%psi))*kind(rst%KSwfn%psi)
         deallocate(rst%KSwfn%psi,stat=i_stat)
         call memocc(i_stat,i_all,'psi',subname)
         i_all=-product(shape(rst%KSwfn%orbs%eval))*kind(rst%KSwfn%orbs%eval)
         deallocate(rst%KSwfn%orbs%eval,stat=i_stat)
         call memocc(i_stat,i_all,'eval',subname)
         nullify(rst%KSwfn%orbs%eval)

        call deallocate_wfd(rst%KSwfn%Lzd%Glr%wfd,subname)
      end if

      if(.not. in%c_absorbtion) then 

         stop 'ERROR'
      else

         call abscalc(nproc,iproc,atoms,rxyz,&
             rst%KSwfn%psi,rst%KSwfn%Lzd,rst%KSwfn%orbs,&
             rst%hx_old,rst%hy_old,rst%hz_old,in,rst%GPU,infocode)
         fxyz(:,:) = 0.d0
      endif

      if (in%inputPsiId==1 .and. infocode==2) then
         if (in%gaussian_help) then
            in%inputPsiId=11
         else
            in%inputPsiId=0
         end if
      else if ((in%inputPsiId==1 .or. in%inputPsiId==0) .and. infocode==1) then
         !in%inputPsiId=0 !better to diagonalise that to restart an input guess
         in%inputPsiId=1
         if(iproc==0) then
            write(*,*)&
               &   ' WARNING: Wavefunctions not converged after cycle',icycle
            write(*,*)' restart after diagonalisation'
         end if

      else if (in%inputPsiId == 0 .and. infocode==3) then
         if (iproc == 0) then
            write( *,'(1x,a)')'Convergence error, cannot proceed.'
            write( *,'(1x,a)')' writing positions in file posfail.xyz then exiting'
            write(comment,'(a)')'UNCONVERGED WF '
            !call wtxyz('posfail',energy,rxyz,atoms,trim(comment))

            call write_atomic_file("posfail",energy,rxyz,atoms,trim(comment))

         end if 

         i_all=-product(shape(rst%KSwfn%psi))*kind(rst%KSwfn%psi)
         deallocate(rst%KSwfn%psi,stat=i_stat)
         call memocc(i_stat,i_all,'psi',subname)
         i_all=-product(shape(rst%KSwfn%orbs%eval))*kind(rst%KSwfn%orbs%eval)
         deallocate(rst%KSwfn%orbs%eval,stat=i_stat)
         call memocc(i_stat,i_all,'eval',subname)
         nullify(rst%KSwfn%orbs%eval)

        call deallocate_wfd(rst%KSwfn%Lzd%Glr%wfd,subname)
         !finalize memory counting (there are still the positions and the forces allocated)
         call memocc(0,0,'count','stop')

         if (nproc > 1) call MPI_FINALIZE(ierr)

         stop 'normal end'
      else
         exit loop_cluster
      end if

   end do loop_cluster

   !preserve the previous value
   in%inputPsiId=inputPsiId_orig

   !put a barrier for all the processes
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE call_abscalc


!>   Absorption (XANES) calculation
!!   @param psi should be freed after use outside of the routine.
!!   @param infocode -> encloses some information about the status of the run
!!          - 0 run succesfully succeded
!!          - 1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
!!               forces may be meaningless   
!!          - 2 (present only for inputPsiId=1) gnrm of the first iteration > 1 AND growing in
!!               the second iteration OR grnm 1st >2.
!!               Input wavefunctions need to be recalculated. Routine exits.
!!          - 3 (present only for inputPsiId=0) gnrm > 4. SCF error. Routine exits.
!!
subroutine abscalc(nproc,iproc,atoms,rxyz,&
     psi,Lzd,orbsAO,hx_old,hy_old,hz_old,in,GPU,infocode)
   use module_base
   use module_types
   use module_interfaces
   use Poisson_Solver
   use module_xc
   use vdwcorrection
   use esatto
   use m_ab6_symmetry
   use m_ab6_mixing
   use m_ab6_kpoints
   implicit none
   integer, intent(in) :: nproc,iproc
   real(gp), intent(inout) :: hx_old,hy_old,hz_old
   type(input_variables), intent(in) :: in
  type(local_zone_descriptors), intent(inout) :: Lzd
   type(atoms_data), intent(inout) :: atoms
   type(orbitals_data), intent(inout) :: orbsAO
   type(GPU_pointers), intent(inout) :: GPU
   real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
   integer, intent(out) :: infocode
   real(wp), dimension(:), pointer :: psi
   !local variables
   type(orbitals_data) :: orbs
   character(len=*), parameter :: subname='abscalc'
   character(len=3) :: PSquiet
   integer :: ixc,ncong,idsx,ncongt,nspin,itermax
   integer :: nvirt,nsym
   integer :: nelec,ndegree_ip,j,n1,n2,n3
!   integer :: n3d,n3p,n3pi,i3xcsh,i3s
   integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i
   integer :: iat,i_all,i_stat,ierr,inputpsi
   real :: tcpu0,tcpu1
   real(gp), dimension(3) :: shift
   real(kind=8) :: crmult,frmult,cpmult,fpmult,gnrm_cv,rbuf,hxh,hyh,hzh,hx,hy,hz
   real(kind=8) :: peakmem
!   real(kind=8) :: eion,epot_sum,ekin_sum,eproj_sum
   real(kind=8) :: tel,psoffset
   !real(gp) :: edisp ! Dispersion energy
   type(nonlocal_psp_descriptors) :: nlpspd
   type(communications_arrays) :: comms
   type(gaussian_basis) :: Gvirt
   type(rho_descriptors)  :: rhodsc
   type(energy_terms) :: energs

   !integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
   real(kind=8), dimension(:,:), allocatable :: radii_cf
   !real(kind=8), dimension(:,:), allocatable :: gxyz
   real(gp), dimension(:,:),pointer :: fdisp,fion
   ! Charge density/potential,ionic potential, pkernel
   real(kind=8), dimension(:), allocatable :: pot_ion

   real(kind=8), dimension(:,:,:,:), allocatable, target :: rhopot, rhopotTOTO, rhoXanes
   real(kind=8), dimension(:,:,:,:), pointer ::  rhopottmp, rhopotExtra, rhotarget
   integer :: b2Bcounter, b2BN
   character(len=100) :: filename
   real(kind=8), dimension(:), pointer :: pkernel

   !wavefunction gradients, hamiltonian on vavefunction
   !transposed  wavefunction
   ! Pointers and variables to store the last psi
   ! before reformatting if useFormattedInput is .true.
   real(kind=8), dimension(:), pointer :: hpsi,psit,psivirt
   real(dp), dimension(:,:,:,:), pointer :: rhocore
   !real(kind=8), dimension(:), pointer :: psidst,hpsidst
   ! PSP projectors 
   real(kind=8), dimension(:), pointer :: proj
   ! arrays for DIIS convergence accelerator
   !real(kind=8), dimension(:,:,:), pointer :: ads
   ! Arrays for the symmetrisation, not used here...
   type(symmetry_data) :: symObj
   type(denspot_distribution) :: dpcom
   character(len=5) :: gridformat

   !for xabsorber
   integer :: ix, iy, iz , ixnl, iynl, iznl !n(c) lpot_a
   real(gp) :: rpot_a,spot_a,hpot_a,espo,harmo,r,rx,ry,rz,minr  
   real(gp), pointer :: radpot(:,:)
   integer :: radpotcount, igrid
   real(gp), dimension(6) :: ewaldstr
   real(gp), dimension(:,:), pointer :: rxyz_b2B
   integer, dimension(:), pointer :: iatype_b2B, znucl_b2B
   real(gp) :: shift_b2B(3)
   integer :: itype, nd
   integer :: n1i_bB,n2i_bB,n3i_bB
   real(gp), dimension(:,:,:,:), pointer :: pot_bB
   real(gp), dimension(:), pointer ::  intfunc_x, intfunc_y
   real(gp) :: factx, facty, factz
   integer :: idelta
   integer :: ix_bB, iy_bB, iz_bB
   integer :: maxX_B, maxY_B, maxZ_B
   integer :: minX_B, minY_B, minZ_B
   real(gp) :: rx_bB, ry_bB, rz_bB
   integer :: nrange
   real(gp), pointer :: auxint(:)

   logical  :: exists 
   integer  :: nat_b2B
   integer  :: Nreplicas, ireplica, replicaoffset
   real(gp) :: dumvect3D(3)
   real(gp) :: shiftdiff
   real(gp) :: potmodified_maxr, potmodified_shift

   type(atoms_data) :: atoms_clone
   integer :: nsp, nspinor !n(c) noncoll
   integer, parameter :: nelecmax=32,lmax=4 !n(c) nmax=6
   integer, parameter :: noccmax=2

   !! to apply pc_projector
   type(pcproj_data_type) ::PPD
   !! to apply paw projectors
   type(PAWproj_data_type) ::PAWD


   if (in%potshortcut==0) then
      if(nproc>1) call MPI_Finalize(ierr)
      stop '   in%potshortcut==0 calculating spectra. Use rather box2Box option      '
   endif

   crmult=in%crmult
   frmult=in%frmult
   cpmult=in%frmult
   fpmult=in%frmult
   ixc=in%ixc
   gnrm_cv=in%gnrm_cv
   itermax=in%itermax
   ncong=in%ncong
   idsx=in%idsx
   rbuf=in%rbuf
   ncongt=in%ncongt
   nspin=in%nspin

   nvirt=in%nvirt

   hx=in%hx
   hy=in%hy
   hz=in%hz

   write(gridformat, "(A)") ""
   select case (in%output_denspot_format)
   case (output_denspot_FORMAT_ETSF)
      write(gridformat, "(A)") ".etsf"
   case (output_denspot_FORMAT_CUBE)
      write(gridformat, "(A)") ".bin"
   end select

   if (iproc == 0) then
      write( *,'(1x,a,1x,i0)') &
         &   '===================== BigDFT XANE calculation =============== inputPsiId=',&
         &   in%inputPsiId
      call print_dft_parameters(in,atoms)
   end if
   !time initialization
   call timing(nproc,trim(in%dir_output)//'time.prc','IN')
   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max)
 
   if(nspin/=1 .and. nspin/=2 .and. nspin/=4) nspin=1

   ! grid spacing (same in x,y and z direction)

   !these routines can be regrouped in one

   allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)


   if (iproc==0) then
      write( *,'(1x,a)')&
           &   '------------------------------------------------------------------ System Properties'
   end if

   if (ixc < 0) then
      call xc_init(ixc, XC_MIXED, nspin)
   else
      call xc_init(ixc, XC_ABINIT, nspin)
   end if

   !character string for quieting the Poisson solver
   if (verbose >1) then
      PSquiet='NO'
   else
      PSquiet='YES'
   end if

   call system_properties(iproc,nproc,in,atoms,orbsAO,radii_cf,nelec)

   call nullify_locreg_descriptors(Lzd%Glr)

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space

   call system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Lzd%Glr,shift)

   if ( orbsAO%nspinor.gt.1) then
      !!  hybrid_on is not compatible with kpoints
     Lzd%Glr%hybrid_on=.false.
   endif

   ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
   call createWavefunctionsDescriptors(iproc,hx,hy,hz,&
       atoms,rxyz,radii_cf,crmult,frmult,Lzd%Glr)

   Lzd%hgrids(1)=hx
   Lzd%hgrids(2)=hy
   Lzd%hgrids(3)=hz

   !variables substitution for the PSolver part
   hxh=0.5d0*hx
   hyh=0.5d0*hy
   hzh=0.5d0*hz
   n1i=Lzd%Glr%d%n1i
   n2i=Lzd%Glr%d%n2i
   n3i=Lzd%Glr%d%n3i

   n1=Lzd%Glr%d%n1
   n2=Lzd%Glr%d%n2
   n3=Lzd%Glr%d%n3

   ! Calculate all projectors, or allocate array for on-the-fly calculation

   !de-allocate orbs and recreate it with one orbital only

   !allocate communications arrays (allocate it before Projectors because of the definition
   !of iskpts and nkptsp)

   call orbitals_descriptors(iproc,nproc,1,1,0,in%nspin,1,in%nkpt,in%kpt,in%wkpt,orbs,.false.)
   call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)  

   call createProjectorsArrays(iproc,Lzd%Glr,rxyz,atoms,orbs,&
        radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)

   call check_linear_and_create_Lzd(iproc,nproc,in%linear,Lzd,atoms,orbs,in%nspin,rxyz)

   !calculate the partitioning of the orbitals between the different processors
   !memory estimation
   if (iproc==0 .and. verbose > 0) then
     call MemoryEstimator(nproc,idsx,Lzd%Glr,&
         &   atoms%nat,orbs%norb,orbs%nspinor,orbs%nkpts,nlpspd%nprojel,&
         &   in%nspin,in%itrpmax,in%iscf,peakmem)
   end if

  !calculate the descriptors for rho and the potentials.
   call denspot_communications(iproc,nproc,Lzd%Glr%d,hxh,hyh,hzh,in,&
        atoms,rxyz,radii_cf,dpcom,rhodsc)

!!$   !these arrays should be included in the comms descriptor
!!$   !allocate values of the array for the data scattering in sumrho
!!$   !its values are ignored in the datacode='G' case
!!$   allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
!!$   call memocc(i_stat,nscatterarr,'nscatterarr',subname)
!!$   !allocate array for the communications of the potential
!!$   allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
!!$   call memocc(i_stat,ngatherarr,'ngatherarr',subname)
!!$   !create the descriptors for the density and the potential
!!$   !these descriptors should take into account the localisation regions
!!$  call createDensPotDescriptors(iproc,nproc,atoms,Lzd%Glr%d,hxh,hyh,hzh,&
!!$       rxyz,in%crmult,in%frmult,radii_cf,in%nspin,'D',ixc,in%rho_commun,&
!!$       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)

   !allocate ionic potential
   if (iproc == 0) write(*,*) " allocate ionic potential " 
   if (dpcom%n3pi > 0) then
      allocate(pot_ion(n1i*n2i*dpcom%n3pi+ndebug),stat=i_stat)
      call memocc(i_stat,pot_ion,'pot_ion',subname)
   else
      allocate(pot_ion(1+ndebug),stat=i_stat)
      call memocc(i_stat,pot_ion,'pot_ion',subname)
   end if

   ! A message about dispersion forces.
   call vdwcorrection_initializeparams(in%ixc, in%dispersion)
   if (iproc == 0) call vdwcorrection_warnings(atoms, in)

   !calculation of the Poisson kernel anticipated to reduce memory peak for small systems
   ndegree_ip=16 !default value 
  call createKernel(iproc,nproc,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,pkernel,&
      &   (verbose > 1))

   !calculate the irreductible zone for this region, if necessary.
   call symmetry_set_irreductible_zone(atoms%sym,atoms%geocode,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i, in%nspin)

!!$   !calculate the irreductible zone for this region, if necessary.
!!$   if (atoms%sym%symObj >= 0) then
!!$      call symmetry_get_n_sym(atoms%sym%symObj, nsym, i_stat)
!!$      if (nsym > 1) then
!!$         ! Current third dimension is set to 1 always
!!$         ! since nspin == nsppol always in BigDFT
!!$         allocate(atoms%sym%irrzon(n1i*n2i*n3i,2,1+ndebug),stat=i_stat)
!!$         call memocc(i_stat,atoms%sym%irrzon,'irrzon',subname)
!!$         allocate(atoms%sym%phnons(2,n1i*n2i*n3i,1+ndebug),stat=i_stat)
!!$         call memocc(i_stat,atoms%sym%phnons,'phnons',subname)
!!$         call kpoints_get_irreductible_zone(atoms%sym%irrzon, atoms%sym%phnons, &
!!$              &   n1i, n2i, n3i, in%nspin, in%nspin, atoms%sym%symObj, i_stat)
!!$      end if
!!$   end if
!!$   if (.not. associated(atoms%sym%irrzon)) then
!!$      ! Allocate anyway to small size otherwise the bounds check does not pass.
!!$      allocate(atoms%sym%irrzon(1,2,1+ndebug),stat=i_stat)
!!$      call memocc(i_stat,atoms%sym%irrzon,'irrzon',subname)
!!$      allocate(atoms%sym%phnons(2,1,1+ndebug),stat=i_stat)
!!$      call memocc(i_stat,atoms%sym%phnons,'phnons',subname)
!!$   end if


   if(sum(atoms%paw_NofL).gt.0) then
      ! Calculate all paw_projectors, or allocate array for on-the-fly calculation
      call timing(iproc,'CrtPawProjects ','ON')
      PAWD%DistProjApply =  .false. !! .true.
      ! the following routine calls a specialized version of localize_projectors
      ! which does not interfere with the global DistProjApply
      call createPawProjectorsArrays(iproc,n1,n2,n3,rxyz,atoms,orbs,&
           radii_cf,cpmult,fpmult,hx,hy,hz, &
           PAWD, Lzd%Glr )
      call timing(iproc,'CrtPawProjects ','OF')
   endif

   if (in%iabscalc_type==3) then
      ! Calculate all pc_projectors, or allocate array for on-the-fly calculation
      call timing(iproc,'CrtPcProjects ','ON')
      PPD%DistProjApply  =  DistProjApply
      ! the following routine calls  localize_projectors again
      ! but this should be in coherence with the previous call for psp projectos 
      call createPcProjectorsArrays(iproc,n1,n2,n3,rxyz,atoms,orbs,&
           radii_cf,cpmult,fpmult,hx,hy,hz,-0.1_gp, &
           PPD, Lzd%Glr  )
      call timing(iproc,'CrtPcProjects ','OF')
   endif

   call IonicEnergyandForces(iproc,nproc,atoms,hxh,hyh,hzh,in%elecfield,rxyz,&
        energs%eion,fion,in%dispersion,energs%edisp,fdisp,ewaldstr,&
        psoffset,n1,n2,n3,n1i,n2i,n3i,dpcom%i3s+dpcom%i3xcsh,dpcom%n3pi,pot_ion,pkernel)

   call createIonicPotential(atoms%geocode,iproc,nproc, (iproc == 0), atoms,rxyz,hxh,hyh,hzh,&
        in%elecfield,n1,n2,n3,dpcom%n3pi,dpcom%i3s+dpcom%i3xcsh,n1i,n2i,n3i,pkernel,pot_ion,psoffset)



   !Allocate Charge density, Potential in real space
   if (dpcom%n3d >0) then
      allocate(rhopot(n1i,n2i,dpcom%n3d,in%nspin+ndebug),stat=i_stat)
      call memocc(i_stat,rhopot,'rhopot',subname)
   else
      allocate(rhopot(1,1,1,in%nspin+ndebug),stat=i_stat)
      call memocc(i_stat,rhopot,'rhopot',subname)
   end if

   if( iand( in%potshortcut,16)>0) then
      allocate(rhoXanes(n1i,n2i,n3i,in%nspin+ndebug),stat=i_stat)
      call memocc(i_stat,rhoXanes,'rhoXanes',subname)
   else
      allocate(rhoXanes(1,1,1,1+ndebug),stat=i_stat)
      call memocc(i_stat,rhoXanes,'rhoXanes',subname)     
   endif


   if(  iand( in%potshortcut, 2)  > 0 ) then
      allocate(rhopotTOTO(n1i,n2i,n3i,  in%nspin+ndebug),  stat=i_stat  )
      call memocc(i_stat,rhopotTOTO,'rhopotTOTO',subname)
   endif


   nullify(rhocore)
   !check the communication distribution
  !call check_communications(iproc,nproc,orbs,Lzd%Glr,comms)


   if( iand( in%potshortcut,4)  .gt. 0 ) then
      if (dpcom%n3d >0) then
         allocate(rhopotExtra(n1i,n2i,dpcom%n3d,in%nspin+ndebug),stat=i_stat)
         call memocc(i_stat,rhopotExtra,'rhopotExtra',subname)
      else
         allocate(rhopotExtra(1,1,1,in%nspin+ndebug),stat=i_stat)
         call memocc(i_stat,rhopotExtra,'rhopotExtra',subname)
      end if

      atoms_clone = atoms
      nullify(atoms_clone%aocc)
      nullify(atoms_clone%iasctype)


      allocate(atoms_clone%aocc(lbound(atoms%aocc,1 ):ubound(atoms%aocc,1),&
         &   lbound(atoms%aocc,2):ubound(atoms%aocc,2)),stat=i_stat)
      call memocc(i_stat,atoms%aocc,'atoms_clone%aocc',subname)

      allocate(atoms_clone%iasctype(lbound(atoms%iasctype,1 ):ubound(atoms%iasctype,1)),stat=i_stat)
      call memocc(i_stat,atoms%iasctype,'atoms_clone%iasctype',subname)


      atoms_clone%aocc=0.0_gp
      atoms_clone%iasctype=0


      read(in%extraOrbital,*,iostat=ierr)iat
      !control the spin
      select case(in%nspin)
      case(1)
         nsp=1
         nspinor=1
         !n(c) noncoll=1
      case(2)
         nsp=2
         nspinor=1
         !n(c) noncoll=1
      case(4)
         nsp=1
         nspinor=4
         !n(c) noncoll=2
      case default
         if (iproc == 0) write(*,*)' ERROR: nspin not valid:',nspin
         stop
      end select

      print *, " Going to create extra potential for orbital "
      print *, in%extraOrbital
      print *, "using hard-coded parameters "
      print *, "noccmax, nelecmax,lmax ", noccmax, nelecmax,lmax

      call read_eleconf(in%extraOrbital ,nsp,nspinor,noccmax, nelecmax,lmax, &
         &   atoms_clone%aocc(1,iat), atoms_clone%iasctype(iat))

      nspin=in%nspin
      symObj%symObj = -1

      !-- calculate input guess from non-diagonalised of LCAO basis (written in wavelets)
      !-- if spectra calculation is energy dependent  input_wf_diag will write
      !    the density to the file electronic_density.cube
      !   To tell  input_wf_diag to do this ne must set the 5th bit on in in%potshortcut
      if( iand( in%potshortcut,16)>0) then
         if( in%iabscalc_type<3) then
            if(iproc==0) write(*,*)  ' Energy dependent potential has been asked in abscalc but  iabscalc_type<3 '
            if(nproc>1) call MPI_Finalize(ierr)
            stop '      Energy dependent potential has been asked in abscalc but  iabscalc_type<3    '
         endif
      endif

      call extract_potential_for_spectra(iproc,nproc,atoms_clone,rhodsc,dpcom,&
          orbsAO,nvirt,comms,Lzd,hx,hy,hz,rxyz,rhopotExtra,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernel,ixc,psi,hpsi,psit,Gvirt,&
          nspin, in%potshortcut, symObj, GPU,in)
      
      if( iand( in%potshortcut,32)  .gt. 0 .and. in%iabscalc_type==3 ) then
         print *, " ============== TESTING PC_PROJECTORS =========== "
         allocate(hpsi(max(orbsAO%npsidim_orbs,orbsAO%npsidim_comp)+ndebug),stat=i_stat)
         hpsi=0.0_wp
         PPD%iproj_to_factor(1:PPD%mprojtot) = 2.0_gp
        call applyPCprojectors(orbsAO,atoms,hx,hy,hz,Lzd%Glr,PPD,psi,hpsi, .true.)
         deallocate(hpsi)
      end if

      if( iand( in%potshortcut,16)>0) then

         if(iproc==0) write(*,*) "re-reading electronic_density for Xanes energy dependent potential "
         STOP " this part has to be rearranged to keep into account distributed potentials "
         call read_density_cube_old("electronic_density",&
              n1i,n2i,n3i,1, hx ,hy ,hz, atoms%nat, rxyz_b2B, pot_bB)
         rhoXanes=0.0_gp
         do iz = 1,n3i
            do iy=1,n2i
               do ix=1,n1i
                  rhopottmp(ix,iy,iz +dpcom%i3xcsh,1) =  pot_bB(ix,iy,iz,1)  !pot_bB(ix+ (iy-1)*n1i  + (iz-1)*n1i*n2i,1)  
               enddo
            enddo
         enddo

         i_all=-product(shape(pot_bB))*kind(pot_bB)
         deallocate(pot_bB,stat=i_stat)
         call memocc(i_stat,i_all,'pot_bB',subname)

         i_all=-product(shape(rxyz_b2B))*kind(rxyz_b2B)
         deallocate(rxyz_b2B,stat=i_stat)
         call memocc(i_stat,i_all,'rxyz',subname)


      endif



      i_all=-product(shape(psi))*kind(psi)
      deallocate(psi,stat=i_stat)
      call memocc(i_stat,i_all,'psi',subname)

      i_all=-product(shape(atoms_clone%aocc))*kind(atoms_clone%aocc)
      deallocate(atoms_clone%aocc,stat=i_stat)
      call memocc(i_stat,i_all,'atoms_clone%aocc',subname)
      nullify(atoms_clone%aocc)
      i_all=-product(shape(atoms_clone%iasctype))*kind(atoms_clone%iasctype)
      deallocate(atoms_clone%iasctype,stat=i_stat)
      call memocc(i_stat,i_all,'atoms_clone%iasctype',subname)
      nullify(atoms_clone%iasctype)

   endif


   if( iand( in%potshortcut,1)  .gt. 0 ) then

      inputpsi=in%inputPsiId

      nspin=in%nspin

      symObj%symObj = -1


      !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
      call extract_potential_for_spectra(iproc,nproc,atoms,rhodsc,dpcom,&
          orbsAO,nvirt,comms,Lzd,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernel,ixc,psi,hpsi,psit,Gvirt,&
          nspin, in%potshortcut, symObj, GPU, in)

      i_all=-product(shape(psi))*kind(psi)
      deallocate(psi,stat=i_stat)
      call memocc(i_stat,i_all,'psi',subname)

   end if

   nullify(psit)

   i_all=-product(shape(pot_ion))*kind(pot_ion)
   deallocate(pot_ion,stat=i_stat)
   call memocc(i_stat,i_all,'pot_ion',subname)

   i_all=-product(shape(pkernel))*kind(pkernel)
   deallocate(pkernel,stat=i_stat)
   call memocc(i_stat,i_all,'kernel',subname)

   ! needs something to let to  bigdft to deallocate
   allocate(psi(2+ndebug),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)

   allocate(orbsAO%eval(2+ndebug),stat=i_stat)
   call memocc(i_stat, orbsAO%eval,'eval',subname)


   if ( in%c_absorbtion ) then


      !!$
      !!$     rhopot(10,9,8+i3xcsh,1)=100.0

      if (in%output_denspot == output_denspot_DENSPOT) then
         if (in%output_denspot_format == output_denspot_FORMAT_TEXT) then
            if (iproc == 0) write(*,*) 'writing local_potential'
            call plot_density(iproc,nproc,'local_potentialb2B' // gridformat,&
                 atoms,rxyz,dpcom,in%nspin,rhopot(1,1,1,1))
         else
            call plot_density_cube_old('local_potentialb2B',iproc,nproc,&
               &   n1,n2,n3,n1i,n2i,n3i,dpcom%n3p,&
               &   in%nspin,hxh,hyh,hzh,atoms,rxyz,dpcom%ngatherarr,rhopot(1,1,1,1))
         endif
      end if

      !!$     call  read_potfile4b2B("local_potential.pot",n1i_bB,n2i_bB,n3i_bB, pot_bB, alat1_bB, alat2_bB, alat3_bB)
      !!$     print *, pot_bB(10  + (9-1)*n1i_bB  + (8-1)*n1i_bB*n2i_bB)
      !!$     stop

      if(iproc==0) print *, " going to calculate spectra "

      if(  iand( in%potshortcut, 2)  > 0 ) then

         if(  iand( in%potshortcut, 16)  > 0 ) then
            b2BN=2
         else
            b2BN=1
         endif

         !Big loop (b2Bcounter)
         do b2Bcounter=1,b2BN

            if(b2Bcounter==1) then
               write(filename,'(A)' ) 'b2B_xanes'
               rhotarget=>rhopotTOTO
            else
               STOP " reimplement allocation of rhoXanesTOTO" 
               write(filename,'(A)') 'b2B_rho'
               !              rhotarget=>rhoXanesTOTO
            endif


            inquire(file=trim(trim(filename)//'.cube'),exist=exists)
            print *, "check ",  trim(filename)//'.cube', exists

            if(exists) then

               nullify(pot_bB);
               call read_cube(trim(filename),atoms%geocode,n1i_bB,n2i_bB,n3i_bB, &
                  &   nspin , hx_old ,hy_old ,hz_old ,pot_bB, nat_b2B, rxyz_b2B, iatype_b2B, znucl_b2B)
               !call read_density_cube_old(trim(filename), n1i_bB,n2i_bB,n3i_bB, 1 , hx_old ,hy_old ,hz_old , nat_b2B, rxyz_b2B, pot_bB )
               hx_old=hx_old*2
               hy_old=hy_old*2
               hz_old=hz_old*2


               if( (atoms%nat/nat_b2B)*nat_b2B /=  atoms%nat ) then
                  if(iproc==0) write(*,*)  "   b2B_xanes cube  is not compatible with actual positions" 
                  if(nproc>1) call MPI_Finalize(ierr)
                  stop '      b2B_xanes cube  is not compatible with actual positions          '
               end if

            else
               print  *, " reading potential from file ","b2B_xanes.pot"
               write(6,*) " reading pot for b2B from potfile has been disactivated. Use cube format instead  " 
               stop  " reading pot for b2B from potfile has been disactivated. Use cube format instead  " 

            endif


            allocate(rhopottmp( max(n1i_bB,n1i),max(n2i_bB,n2i),max(n3i_bB,n3i),in%nspin+ndebug),stat=i_stat)
            !! allocate(rhopottmp( max(n1i_bB,n1i),max(n2i_bB,n2i),max(n3i_bB,n3d),in%nspin+ndebug),stat=i_stat)
            call memocc(i_stat,rhopottmp,'rhopottmp',subname)


            rhotarget=0.0_gp


            itype=16
            nd=2**20

            allocate(  intfunc_x(0:nd+ndebug),stat=i_stat )
            call memocc(i_stat,intfunc_x,'intfunc_x',subname)
            allocate( intfunc_y(0:nd+ndebug) ,stat=i_stat )
            call memocc(i_stat,intfunc_y,'intfunc_y',subname)

            print *, " scaling function for interpolation "

            call scaling_function4b2B(itype,nd,nrange,intfunc_x,intfunc_y)  ! intervallo di 32 con 2**20 punti
            if( abs(intfunc_y(nd/2)-1)>1.0e-10 ) then
               stop " wrong scaling function 4b2B: not a centered one "
            endif

            i_all=-product(shape(intfunc_x))*kind(intfunc_x)
            deallocate(intfunc_x,stat=i_stat)
            call memocc(i_stat,i_all,'intfunc_x',subname)

            allocate(auxint(n1i+n2i+n3i+ndebug),stat=i_stat)
            call memocc(i_stat,auxint,'auxint',subname)

            Nreplicas = atoms%nat / nat_b2B
            dumvect3d(1)=atoms%alat1
            dumvect3d(2)=atoms%alat2
            dumvect3d(3)=atoms%alat3


            !Loop over ireplica
            do ireplica=0, Nreplicas-1

               replicaoffset = (ireplica)*(   nat_b2B      )

               do j=1,3
                  shift_b2B(j) = rxyz(j,1+replicaoffset) - rxyz_b2B(j,1) 

                  do iat=2+ireplica*nat_b2B, (ireplica+1)*nat_b2B

                     shiftdiff = shift_b2B(j) - (rxyz(j,iat) -rxyz_b2B(j,iat-replicaoffset)   )

                     if( abs( shiftdiff )>1.0e-4 .and.  abs(abs(shiftdiff)-dumvect3d(j))  >1.0e-4) then
                        if(iproc==0) write(*,*)  "   b2B_xanes  positions are not compatible with actual positions" 
                        if(nproc>1) call MPI_Finalize(ierr)
                        stop '      b2B_xanes positions are not compatible with actual positions          '
                     end if
                  enddo
               enddo

               if (iproc == 0) write(*,'(a,1x,i2,1x,a,1x,3ES13.6)')  "for replica ", ireplica,  "SHIFT " , shift_b2B

               rhopottmp=0.0_gp
               do iz_bB = 1,n3i_bB
                  do iy_bB=1,n2i_bB
                     do ix_bB=1,n1i_bB
                        rhopottmp(ix_bB,iy_bB,iz_bB ,1) =  pot_bB(ix_bB,iy_bB,iz_bB,1)
                                                          !pot_bB(ix_bB  + (iy_bB-1)*n1i_bB  + (iz_bB-1)*n1i_bB*n2i_bB,1)
                        !! rhopottmp(ix_bB,iy_bB,iz_bB ,1) =  pot_bB(ix_bB  + (iy_bB-1)*n1i_bB  + (iz_bB-1)*n1i_bB*n2i_bB,1)
                     enddo
                  enddo
               enddo

               if (ireplica==Nreplicas-1) then
                  i_all=-product(shape(pot_bB))*kind(pot_bB)
                  deallocate(pot_bB,stat=i_stat)
                  call memocc(i_stat,i_all,'pot_bB',subname)
               endif


               do iz_bB = 1,n3i_bB
                  do iy_bB=1,n2i_bB
                     auxint = 0.0_gp
                     do ix_bB=1,n1i_bB
                        rx_bB = hx_old*(ix_bB-1)           /2.0   +  shift_b2B(1)
                        minX_B  =  max(1,NINT((rx_bB -8*hx_old/2)/(hx/2.0)))
                        maxX_B  =  min(n1i,NINT((rx_bB +8*hx_old/2)/(hx/2.0)))

                        minX_B  =  NINT((rx_bB -8*hx_old/2)/(hx/2.0))
                        maxX_B  =  NINT((rx_bB +8*hx_old/2)/(hx/2.0))

                        do ixnl= minX_B , maxX_B 
                           ix = mod(ixnl-1 + n1i , n1i) +1

                           rx = hx*(ix-1  )/2.0  

                           shiftdiff = (rx-rx_bB)
                           if ( abs(shiftdiff -atoms%alat1) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat1
                           if ( abs(shiftdiff +atoms%alat1) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat1

                           idelta = NINT( shiftdiff *2**15/(hx_old/2))  
                           factx = intfunc_y(nd/2+idelta)
                           auxint(ix) = auxint(ix) + &
                              &   factx * rhopottmp(ix_bB,iy_bB,iz_bB ,1)
                        enddo
                     enddo
                     rhopottmp(1:n1i,iy_bB,iz_bB,1)=auxint(1:n1i)
                  enddo
               enddo

               do iz_bB = 1,n3i_bB
                  do ix_bB=1,n1i
                     auxint = 0.0_gp
                     do iy_bB=1,n2i_bB
                        ry_bB = hy_old*(iy_bB-1)/2.0   +  shift_b2B(2)
                        minY_B  =  max(1  ,NINT((ry_bB -8*hy_old/2)/(hy/2.0)))
                        maxY_B  =  min(n2i,NINT((ry_bB +8*hy_old/2)/(hy/2.0)))

                        minY_B  =  NINT((ry_bB -8*hy_old/2)/(hy/2.0))
                        maxY_B  =  NINT((ry_bB +8*hy_old/2)/(hy/2.0))

                        do iynl= minY_B , maxY_B 
                           iy = mod(iynl-1 + n2i , n2i) +1

                           ry = hy*(iy-1  )/2.0  

                           shiftdiff = (ry-ry_bB)
                           if ( abs(shiftdiff -atoms%alat2) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat2
                           if ( abs(shiftdiff +atoms%alat2) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat2


                           idelta = NINT(shiftdiff *2**15/(hy_old/2))
                           facty = intfunc_y(nd/2+idelta)
                           auxint(iy) = auxint(iy) + &
                              &   facty * rhopottmp(ix_bB,iy_bB,iz_bB,1)
                        enddo
                     enddo
                     rhopottmp(ix_bB ,1:n2i,iz_bB,1)=auxint(1:n2i)
                  enddo
               enddo

               do ix_bB=1,n1i
                  do iy_bB=1,n2i
                     auxint = 0.0_gp
                     do iz_bB = 1,n3i_bB
                        rz_bB = hz_old*(iz_bB-1)           /2.0   +  shift_b2B(3)


                        minZ_B  =    NINT((rz_bB -8*hz_old/2)/(hz/2.0))
                        maxZ_B  =    NINT((rz_bB +8*hz_old/2)/(hz/2.0))

                        do iznl= minZ_B , maxZ_B 

                           iz = mod(iznl-1 + n3i , n3i) +1

                           rz = hz*(iz-1  )/2.0  

                           shiftdiff = (rz-rz_bB)
                           if ( abs(shiftdiff -atoms%alat3) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat3
                           if ( abs(shiftdiff +atoms%alat3) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat3

                           idelta = NINT( shiftdiff *2**15/(hz_old/2.0))     
                           factz = intfunc_y(nd/2+idelta)
                           auxint(iz) = auxint(iz) + &
                              &   factz * rhopottmp(ix_bB,iy_bB,iz_bB,1)
                        enddo
                     enddo
                     rhotarget(ix_bB ,iy_bB, 1:n3i ,1)= rhotarget(ix_bB ,iy_bB, : ,1)+auxint(1:n3i)
                  enddo
               enddo
            enddo !End of loop over ireplica

            !De-allocations
            i_all=-product(shape(rxyz_b2B))*kind(rxyz_b2B)
            deallocate(rxyz_b2B,stat=i_stat)
            call memocc(i_stat,i_all,'rxyz',subname)
            i_all=-product(shape(iatype_b2B))*kind(iatype_b2B)
            deallocate(iatype_b2B,stat=i_stat)
            call memocc(i_stat,i_all,'iatype',subname)
            i_all=-product(shape(znucl_b2B))*kind(znucl_b2B)
            deallocate(znucl_b2B,stat=i_stat)
            call memocc(i_stat,i_all,'znucl',subname)
            i_all=-product(shape(rhopottmp))*kind(rhopottmp)
            deallocate(rhopottmp,stat=i_stat)
            call memocc(i_stat,i_all,'rhopottmp',subname)
            i_all=-product(shape(auxint))*kind(auxint)
            deallocate(auxint,stat=i_stat)
            call memocc(i_stat,i_all,'auxint',subname)
            i_all=-product(shape(intfunc_y))*kind(intfunc_y)
            deallocate(intfunc_y,stat=i_stat)
            call memocc(i_stat,i_all,'intfunc_y',subname)
         enddo  !End of loop of B2counter




         ! if (iproc == 0) write(*,*) 'writing NEW local_potential.pot'
         ! call plot_density('local_potentialb2BNEW',iproc,nproc,&
               !      n1,n2,n3,n1i,n2i,n3i,n3p,&
               !      in%nspin,hxh,hyh,hzh,&
               !      atoms,rxyz,ngatherarr,rhopot(1,1,1,1))



         do iz = 1,dpcom%n3p
            do iy = 1,n2i
               do ix = 1,n1i
                  rhopot(ix ,iy , iz  ,1)= rhopotTOTO(ix ,iy,iz+ dpcom%i3s+dpcom%i3xcsh-1  ,1)
               end do
            end do
         end do




      endif !End of if ( iand( in%potshortcut, 2)  > 0)


      if( iand( in%potshortcut,4)  .gt. 0 ) then
         do ix=1,n1i
            do iy=1,n2i
               do iz = 1,n3i
                  rhopot(ix ,iy, iz ,1)= rhopot(ix ,iy, iz ,1)+rhopotExtra(ix ,iy, iz ,1)
               enddo
            enddo
         enddo
         i_all=-product(shape(rhopotExtra))*kind(rhopotExtra)
         deallocate(rhopotExtra,stat=i_stat)
         call memocc(i_stat,i_all,'rhopotExtra',subname)
      endif

      alteration: if(  in%abscalc_alterpot) then
         ! Attention :  modification of the  potential for the  
         ! exactly resolvable case 

         !n(c) lpot_a=1
         rpot_a = 7.5d0
         spot_a = 0.8d0
         hpot_a = 3.0d0

         allocate(radpot(60000 ,2+ndebug ))
         radpotcount=60000

         open(unit=22,file='pot.dat', status='old')
         do igrid=1, radpotcount
            read(22,*)  radpot(igrid ,1 ),  radpot(igrid , 2 )
         enddo
         close(unit=22)

         minr=1000.0
         potmodified_maxr=0
         potmodified_shift=0



         do ix=1,n1i
            do iy=1,n2i
               do iz = 1,dpcom%n3p
                  rx = hx*(ix-1)           /2.0  -  rxyz(1,in%iat_absorber )
                  ry = hy*(iy-1)           /2.0  -  rxyz(2,in%iat_absorber )
                  rz = hz*(iz-1 +dpcom%i3xcsh + dpcom%i3s -1 )/2.0  -  rxyz(3,in%iat_absorber )

                  r  = sqrt( rx*rx+ry*ry+rz*rz)

                  if(r>3.5) then

                     if( r>29) then
                        rhopot(ix,iy,iz,1)=0.0
                     else
                        igrid = binary_search( r, radpot, radpotcount )
                        rhopot(ix,iy,iz,1) = &
                           &   ( radpot(igrid,2)*(radpot(igrid+1,1)-R) + radpot(igrid+1,2)*(R-radpot(igrid,1)) )/&
                           &   ( radpot(igrid+1,1) -radpot(igrid,1) )
                     endif
                  else
                     if(potmodified_maxr<r) then 
                        potmodified_maxr=r
                        igrid = binary_search( r, radpot, radpotcount )
                        potmodified_shift =&
                           &   ( radpot(igrid,2)*(radpot(igrid+1,1)-R) + radpot(igrid+1,2)*(R-radpot(igrid,1)) )/&
                           &   ( radpot(igrid+1,1) -radpot(igrid,1) ) &
                           &   -rhopot(ix,iy,iz,1) 
                     endif
                  endif

                  if(r<minr) minr=r

                  if( r.ge.3.5) then
                     ! harmo = (rx+2*ry+3*rz)/sqrt(14.0)/r *sqrt( 3.0/4.0/3.1415926535)
                     !! harmo = sqrt( 1.0/4.0/3.1415926535)
                     ! harmo = (rz)/sqrt(1.0)/r *sqrt( 3.0/4.0/3.1415926535)
                     harmo = (rx)/sqrt(1.0)/r *sqrt( 3.0/4.0/3.1415926535)
                  else
                     harmo=0.0_gp
                  endif

                  espo  = ((r-rpot_a)**2)/spot_a/spot_a/2.0
                  if(espo<100) then
                     rhopot(ix,iy,iz,1) = rhopot(ix,iy,iz,1) +  hpot_a * exp(-espo) *harmo
                  endif
               enddo
            enddo
         enddo
         do ix=1,n1i
            do iy=1,n2i
               do iz = 1,dpcom%n3p
                  rx = hx*(ix-1)           /2.0  -  rxyz(1,in%iat_absorber )
                  ry = hy*(iy-1)           /2.0  -  rxyz(2,in%iat_absorber )
                  rz = hz*(iz-1 +dpcom%i3xcsh + dpcom%i3s -1 )/2.0  -  rxyz(3,in%iat_absorber )

                  r  = sqrt( rx*rx+ry*ry+rz*rz)

                  if(r<=3.5) then
                     rhopot(ix,iy,iz,1)=rhopot(ix,iy,iz,1)+potmodified_shift*0
                  endif
               enddo
            enddo
         enddo
         print *, "  potmodified_shift =", potmodified_shift
      end if alteration

      infocode=0

      if (in%iabscalc_type==2) then
         call xabs_lanczos(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Lzd,dpcom,&
             rhopot(1,1,1,1),energs,in%nspin,GPU,&
             in%iat_absorber,in,PAWD,orbs)

      else if (in%iabscalc_type==1) then
         call xabs_chebychev(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Lzd,dpcom,&
            &   rhopot(1,1,1,1) ,energs,in%nspin,GPU &
            &   , in%iat_absorber, in, PAWD, orbs)
      else if (in%iabscalc_type==3) then
         call xabs_cg(iproc,nproc,atoms,hx,hy,hz,rxyz,&
             radii_cf,nlpspd,proj,Lzd,dpcom,&
            &   rhopot(1,1,1,1) ,energs,in%nspin,GPU &
            &   , in%iat_absorber, in, rhoXanes(1,1,1,1), PAWD, PPD, orbs)
      else
         if (iproc == 0) write(*,*)' iabscalc_type not known, does not perform calculation'
      endif

   end if

   !    No tail calculation
   if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   call deallocate_before_exiting
!   call deallocate_local_zone_descriptors(lzd, subname)

   contains


   !routine which deallocate the pointers and the arrays before exiting 
   subroutine deallocate_before_exiting

      !when this condition is verified we are in the middle of the SCF cycle

      !! if (infocode /=0 .and. infocode /=1) then
      if (.true.) then

         !!$       if (idsx_actual > 0) then
         !!$          i_all=-product(shape(psidst))*kind(psidst)
         !!$          deallocate(psidst,stat=i_stat)
         !!$          call memocc(i_stat,i_all,'psidst',subname)
         !!$          i_all=-product(shape(hpsidst))*kind(hpsidst)
         !!$          deallocate(hpsidst,stat=i_stat)
         !!$          call memocc(i_stat,i_all,'hpsidst',subname)
         !!$          i_all=-product(shape(ads))*kind(ads)
         !!$          deallocate(ads,stat=i_stat)
         !!$          call memocc(i_stat,i_all,'ads',subname)
         !!$       end if

         !!$       if (nproc > 1) then
         !!$          i_all=-product(shape(psit))*kind(psit)
         !!$          deallocate(psit,stat=i_stat)
         !!$          call memocc(i_stat,i_all,'psit',subname)
         !!$       end if
         !!$       
         !!$       i_all=-product(shape(hpsi))*kind(hpsi)
         !!$       deallocate(hpsi,stat=i_stat)
         !!$       call memocc(i_stat,i_all,'hpsi',subname)




         !!$       i_all=-product(shape(pot_ion))*kind(pot_ion)
         !!$       deallocate(pot_ion,stat=i_stat)
         !!$       call memocc(i_stat,i_all,'pot_ion',subname)
         !!$       
         !!$       i_all=-product(shape(pkernel))*kind(pkernel)
         !!$       deallocate(pkernel,stat=i_stat)
         !!$       call memocc(i_stat,i_all,'pkernel',subname)
         !!$

         !!$       if (in%read_ref_den) then
         !!$          i_all=-product(shape(pkernel_ref))*kind(pkernel_ref)
         !!$          deallocate(pkernel_ref,stat=i_stat)
         !!$          call memocc(i_stat,i_all,'pkernel_ref',subname)
         !!$       end if

         ! calc_tail false
         i_all=-product(shape(rhopot))*kind(rhopot)
         deallocate(rhopot,stat=i_stat)
         call memocc(i_stat,i_all,'rhopot',subname)

         if(  iand( in%potshortcut, 2)  > 0 ) then
            i_all=-product(shape(rhopotTOTO))*kind(rhopotTOTO)
            deallocate(rhopotTOTO,stat=i_stat)
            call memocc(i_stat,i_all,'rhopotTOTO',subname)
         endif

         i_all=-product(shape(rhoXanes))*kind(rhoXanes)
         deallocate(rhoXanes,stat=i_stat)
         call memocc(i_stat,i_all,'rhoXanes',subname)

         !!$       if (in%read_ref_den) then
         !!$          i_all=-product(shape(rhoref))*kind(rhoref)
         !!$          deallocate(rhoref,stat=i_stat)
         !!$          call memocc(i_stat,i_all,'rhoref',subname)
         !!$       end if

         i_all=-product(shape(dpcom%nscatterarr))*kind(dpcom%nscatterarr)
         deallocate(dpcom%nscatterarr,stat=i_stat)
         call memocc(i_stat,i_all,'dpcom%nscatterarr',subname)

         i_all=-product(shape(dpcom%ngatherarr))*kind(dpcom%ngatherarr)
         deallocate(dpcom%ngatherarr,stat=i_stat)
         call memocc(i_stat,i_all,'dpcom%ngatherarr',subname)

         i_all=-product(shape(fion))*kind(fion)
         deallocate(fion,stat=i_stat)
         call memocc(i_stat,i_all,'fion',subname)

         i_all=-product(shape(fdisp))*kind(fdisp)
         deallocate(fdisp,stat=i_stat)
         call memocc(i_stat,i_all,'fdisp',subname)


      end if

      !deallocate wavefunction for virtual orbitals
      !if it is the case
      if (in%nvirt > 0) then
         !call deallocate_gwf(Gvirt,subname)
         i_all=-product(shape(psivirt))*kind(psivirt)
         deallocate(psivirt,stat=i_stat)
         call memocc(i_stat,i_all,'psivirt',subname)
      end if

      !De-allocations
      call deallocate_bounds(atoms%geocode,Lzd%Glr%hybrid_on,&
           Lzd%Glr%bounds,subname)
      call deallocate_Lzd_except_Glr(Lzd, subname)
      i_all=-product(shape(Lzd%Glr%projflg))*kind(Lzd%Glr%projflg)
      deallocate(Lzd%Glr%projflg,stat=i_stat)
      call memocc(i_stat,i_all,'Lzd%Glr%projflg',subname)  

      call deallocate_comms(comms,subname)

      call deallocate_orbs(orbs,subname)
      call deallocate_orbs(orbsAO,subname)

      call deallocate_proj_descr(nlpspd,subname)

      i_all=-product(shape(proj))*kind(proj)
      deallocate(proj,stat=i_stat)
      call memocc(i_stat,i_all,'proj',subname)

      i_all=-product(shape(radii_cf))*kind(radii_cf)
      deallocate(radii_cf,stat=i_stat)
      call memocc(i_stat,i_all,'radii_cf',subname)

      call deallocate_rho_descriptors(rhodsc,subname)

      if( in%iabscalc_type==3) then
         call deallocate_pcproj_data(PPD,subname)
      endif
      if(sum(atoms%paw_NofL).gt.0) then
         call deallocate_pawproj_data(PAWD,subname)       
      endif
      !! this is included in deallocate_atomdatapaw
      !! call deallocate_atomdatapaw(atoms,subname)
     
      ! Free the libXC stuff if necessary.
      call xc_end()

      !end of wavefunction minimisation
      call timing(iproc,'LAST','PR')
      call timing(iproc,'              ','RE')
      call cpu_time(tcpu1)
      call system_clock(ncount1,ncount_rate,ncount_max)
      tel=dble(ncount1-ncount0)/dble(ncount_rate)
      if (iproc == 0) &
         &   write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0

   END SUBROUTINE deallocate_before_exiting

END SUBROUTINE abscalc


subroutine zero4b2B(n,x)
   implicit none
   !Arguments
   integer, intent(in) :: n
   real(kind=8), intent(out) :: x(n)
   !Local variables
   integer :: i
   do i=1,n
      x(i)=0.d0
   end do
END SUBROUTINE zero4b2B


!> Backward wavelet transform
!! @param nd length of data set
!! @param nt length of data in data set to be transformed
!! @param m filter length (m has to be even!)
!! @param x input data, y output data
subroutine back_trans_14_4b2B(nd,nt,x,y)
   implicit none
   !Arguments
   integer, intent(in) :: nd,nt
   real(kind=8), intent(in) :: x(0:nd-1)
   real(kind=8), intent(out) :: y(0:nd-1)
   !Local variables
   integer :: i,j,ind

   include 'lazy_16.inc'

   do i=0,nt/2-1
      y(2*i+0)=0.d0
      y(2*i+1)=0.d0

      do j=-m/2,m/2-1

         ! periodically wrap index if necessary
         ind=i-j
         loop99: do
            if (ind.lt.0) then 
               ind=ind+nt/2
               cycle loop99
            end if
            if (ind.ge.nt/2) then 
               ind=ind-nt/2
               cycle loop99
            end if
            exit loop99
         end do loop99

         y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
         y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
      end do
   end do

END SUBROUTINE back_trans_14_4b2B


subroutine scaling_function4b2B(itype,nd,nrange,a,x)
   use module_base
   implicit none
   !Arguments
   !Type of interpolating functions
   integer, intent(in) :: itype
   !Number of points: must be 2**nex
   integer, intent(in) :: nd
   integer, intent(out) :: nrange
   real(kind=8), dimension(0:nd), intent(out) :: a,x
   !Local variables
   character(len=*), parameter :: subname='scaling_function4b2B'
   real(kind=8), dimension(:), allocatable :: y
   integer :: i,nt,ni,i_all,i_stat  

   !Only itype=8,14,16,20,24,30,40,50,60,100
   select case(itype)
   case(8,14,16,20,24,30,40,50,60,100)
      !O.K.
   case default
      print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
      stop
   end select
   !!$  write(unit=*,fmt="(1x,a,i0,a)") &
   !!$       "Use interpolating scaling functions of ",itype," order"

   !Give the range of the scaling function
   !from -itype to itype
   ni=2*itype
   nrange = ni
   allocate(y(0:nd+ndebug),stat=i_stat)
   call memocc(i_stat,y,'y',subname)

   ! plot scaling function
   call zero4b2B(nd+1,x)
   call zero4b2B(nd+1,y)
   nt=ni
   x(nt/2)=1.d0
   loop1: do
      nt=2*nt
      ! write(6,*) 'nd,nt',nd,nt
      select case(itype)
      case(8)
         stop
      case(14)
         stop
      case(16)
         call back_trans_14_4b2B(nd,nt,x,y)
      case(20)
         stop
      case(24)
         stop
      case(30)
         stop
      case(40)
         stop
      case(50)
         stop
      case(60)
         stop
      case(100)
         stop
      end select

      do i=0,nt-1
         x(i)=y(i)
      end do
      if (nt.eq.nd) then
         exit loop1
      end if
   end do loop1

   !open (unit=1,file='scfunction',status='unknown')
   do i=0,nd
      a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
      !write(1,*) a(i),x(i)
   end do
   !close(1)

   i_all=-product(shape(y))*kind(y)
   deallocate(y,stat=i_stat)
   call memocc(i_stat,i_all,'y',subname)
END SUBROUTINE scaling_function4b2B


subroutine read_potfile4b2B(filename,n1i,n2i,n3i, rho, alat1, alat2, alat3)
   use module_base
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(out) :: n1i,n2i,n3i
   real(gp) alat1, alat2, alat3, dum, dum1
   ! real(dp), dimension(n1i*n2i*n3d), intent(out) :: rho
   real(gp), pointer :: rho(:)
   !local variables
   integer :: nl1,nl2,nl3,i_stat,i1,i2,i3,ind
   real(gp) :: value
   character(len=*), parameter :: subname='read_potfile4b2B'

   open(unit=22,file=filename,status='unknown')
   read(22,*)!'normalised density'
   read(22,*) n1i,n2i,n3i
   read(22,*) alat1,dum ,alat2
   read(22,*)  dum, dum1, alat3
   read(22,*)!xyz   periodic' !not true in general but needed in the case

   !conditions for periodicity in the three directions
   !value of the buffer in the x and z direction
   nl1=1
   nl3=1
   nl2=1

   print *, " allocation for rho for  n1i,n2i,n3i ",  n1i,n2i,n3i

   allocate( rho( n1i*n2i*n3i+ndebug) , stat=i_stat )
   call memocc(i_stat,rho,'rho',subname)

   print *, " going to read all pot points " 
   do i3=0,n3i-1
      do i2=0,n2i-1
         do i1=0,n1i-1
            ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
            read(22,*)value
            rho(ind)=value
         end do
      end do
   end do
   print *, " closing file  " 
   close(22)

END SUBROUTINE read_potfile4b2B


subroutine extract_potential_for_spectra(iproc,nproc,at,rhod,dpcom,&
     orbs,nvirt,comms,Lzd,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
     nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
     nspin,potshortcut,symObj,GPU,input)
   use module_base
   use module_interfaces, except_this_one => extract_potential_for_spectra
   use module_types
   use Poisson_Solver
   use libxc_functionals
   implicit none
   !Arguments
   integer, intent(in) :: iproc,nproc,ixc
   integer, intent(inout) :: nspin,nvirt
   real(gp), intent(in) :: hx,hy,hz
   type(atoms_data), intent(inout) :: at
   type(rho_descriptors),intent(in) :: rhod
   type(denspot_distribution), intent(in) :: dpcom
   type(orbitals_data), intent(inout) :: orbs
   type(nonlocal_psp_descriptors), intent(in) :: nlpspd
   type(local_zone_descriptors), intent(inout) :: Lzd
   type(communications_arrays), intent(in) :: comms
   type(GPU_pointers), intent(inout) :: GPU
   type(input_variables):: input
   type(symmetry_data), intent(in) :: symObj
   !integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
   !integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
   real(gp), dimension(3,at%nat), intent(in) :: rxyz
   real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
   real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
   type(gaussian_basis), intent(out) :: G !basis for davidson IG
   real(wp), dimension(:), pointer :: psi,hpsi,psit
   real(wp), dimension(:,:,:,:), pointer :: rhocore
   real(dp), dimension(:), pointer :: pkernel,pkernelseq
   integer, intent(in) ::potshortcut

  !local variables
  character(len=*), parameter :: subname='extract_potential_for_spectra'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0,ncplx,nrhodim,irhotot_add,irho_add,ispin
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,ehart,eexcu,vexcu,etol,accurex
  type(orbitals_data) :: orbse
  type(communications_arrays) :: commse
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc,passmat
  !real(wp), dimension(:,:,:), allocatable :: mom_vec
  real(gp), dimension(:), allocatable :: locrad
  !   real(wp), dimension(:), pointer :: pot,pot1
  real(dp), dimension(:,:), pointer :: rho_p
  real(wp), dimension(:,:,:), pointer :: psigau
  type(confpot_data), dimension(:), allocatable :: confdatarr
  ! #### Linear Scaling Variables
  integer :: ilr,ityp
  logical :: calc                           
  real(dp),dimension(:),pointer:: Lpsi,Lhpsi
  logical :: linear, withConfinement
  integer,dimension(:),allocatable :: norbsc
  logical,dimension(:),allocatable:: calculateBounds
  integer :: i!,iorb,jorb,icplx
  real(gp), dimension(:), allocatable :: ovrlp
  real(gp), dimension(:,:), allocatable :: smat,tmp
  real(gp), dimension(6) :: xcstr
  type(local_zone_descriptors) :: Lzde

  allocate(norbsc_arr(at%natsc+1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          &   '------------------------------------------------------- Input Wavefunctions Creation'
     !yaml_output
     !      write(70,'(a)')repeat(' ',yaml_indent)//'- Input Hamiltonian: { '
!     yaml_indent=yaml_indent+2 !list element
  end if
  !spin for inputguess orbitals
  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin_ig,&
       &   orbs,orbse,norbsc_arr,locrad,G,psigau,eks)

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbse,commse,subname)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbse,commse,basedist=comms%nvctr_par(0:,1:))  

  !use the eval array of orbse structure to save the original values
  allocate(orbse%eval(orbse%norb*orbse%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbse%eval,'orbse%eval',subname)

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  if(potshortcut<=0) then
     call nullify_local_zone_descriptors(Lzde)
     call create_LzdLIG(iproc,nproc,orbs%nspin,input%linear,hx,hy,hz,Lzd%Glr,at,orbse,rxyz,Lzde)
  else
     call nullify_local_zone_descriptors(Lzde)
     Lzde = Lzd
  end if

  ! determine the wavefunction dimension
  call wavefunction_dimension(Lzde,orbse)

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(psi(max(orbse%npsidim_orbs,orbse%npsidim_comp)+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  !allocate arrays for the GPU if a card is present
  switchGPUconv=.false.
  switchOCLconv=.false.
  if (GPUconv .and. potshortcut ==0 ) then
     call prepare_gpu_for_locham(Lzde%Glr%d%n1,Lzde%Glr%d%n2,Lzde%Glr%d%n3,nspin_ig,&
          hx,hy,hz,Lzd%Glr%wfd,orbse,GPU)
  else if (OCLconv .and. potshortcut ==0) then
     call allocate_data_OCL(Lzde%Glr%d%n1,Lzde%Glr%d%n2,Lzde%Glr%d%n3,at%geocode,&
          nspin_ig,Lzde%Glr%wfd,orbse,GPU)
     if (iproc == 0) write(*,*)&
          'GPU data allocated'
  else if (GPUconv .and. potshortcut >0 ) then
     switchGPUconv=.true.
     GPUconv=.false.
  else if (OCLconv .and. potshortcut >0 ) then
     switchOCLconv=.true.
     OCLconv=.false.
  end if

  call timing(iproc,'wavefunction  ','ON')   
  !use only the part of the arrays for building the hamiltonian matrix
  call gaussians_to_wavelets_new(iproc,nproc,Lzde,orbse,G,&
       psigau(1,1,min(orbse%isorb+1,orbse%norb)),psi)
  call timing(iproc,'wavefunction  ','OF')
  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  !spin adaptation for the IG in the spinorial case
  orbse%nspin=nspin
  call sumrho(iproc,nproc,orbse,Lzde,hxh,hyh,hzh,dpcom%nscatterarr,&
       GPU,symObj,rhod,psi,rho_p)
  call communicate_density(iproc,nproc,orbse%nspin,hxh,hyh,hzh,Lzde,&
       rhod,dpcom%nscatterarr,rho_p,rhopot,.false.)
  orbse%nspin=nspin_ig

  !-- if spectra calculation uses a energy dependent potential
  !    input_wf_diag will write (to be used in abscalc)
  !    the density to the file electronic_density.cube
  !  The writing is activated if  5th bit of  in%potshortcut is on.
  if( iand( potshortcut,16)==0 .and. potshortcut /= 0) then
     call plot_density_cube_old('electronic_density',&
          iproc,nproc,Lzde%Glr%d%n1,Lzde%Glr%d%n2,Lzde%Glr%d%n3,&
          Lzde%Glr%d%n1i,Lzde%Glr%d%n2i,Lzde%Glr%d%n3i,dpcom%nscatterarr(iproc,2),  & 
          nspin,hxh,hyh,hzh,at,rxyz,dpcom%ngatherarr,&
          rhopot(1+dpcom%nscatterarr(iproc,4)*Lzde%Glr%d%n1i*Lzd%Glr%d%n2i))
  endif
  !---

  if(orbs%nspinor==4) then
     !this wrapper can be inserted inside the poisson solver 
     call PSolverNC(at%geocode,'D',iproc,nproc,Lzde%Glr%d%n1i,Lzde%Glr%d%n2i,Lzde%Glr%d%n3i,&
          dpcom%nscatterarr(iproc,1),& !this is n3d
          ixc,hxh,hyh,hzh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
  else
     !Allocate XC potential
     if (dpcom%nscatterarr(iproc,2) >0) then
        allocate(potxc(Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*dpcom%nscatterarr(iproc,2)*nspin+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     else
        allocate(potxc(1+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     end if

     call XC_potential(at%geocode,'D',iproc,nproc,&
          Lzde%Glr%d%n1i,Lzde%Glr%d%n2i,Lzde%Glr%d%n3i,ixc,hxh,hyh,hzh,&
          rhopot,eexcu,vexcu,nspin,rhocore,potxc,xcstr)
     if( iand(potshortcut,4)==0) then
        call H_potential(at%geocode,'D',iproc,nproc,&
             Lzde%Glr%d%n1i,Lzde%Glr%d%n2i,Lzde%Glr%d%n3i,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
     endif

     !sum the two potentials in rhopot array
     !fill the other part, for spin, polarised
     if (nspin == 2) then
        call dcopy(Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*dpcom%nscatterarr(iproc,2),rhopot(1),1,&
             rhopot(Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*dpcom%nscatterarr(iproc,2)+1),1)
     end if
     !spin up and down together with the XC part
     call axpy(Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*dpcom%nscatterarr(iproc,2)*nspin,1.0_dp,potxc(1),1,&
          rhopot(1),1)


     i_all=-product(shape(potxc))*kind(potxc)
     deallocate(potxc,stat=i_stat)
     call memocc(i_stat,i_all,'potxc',subname)

  end if

  if (switchGPUconv) then
     GPUconv=.true.
  end if
  if (switchOCLconv) then
     OCLconv=.true.
  end if

  call deallocate_orbs(orbse,subname)
  i_all=-product(shape(orbse%eval))*kind(orbse%eval)
  deallocate(orbse%eval,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%eval',subname)


  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)
  if(potshortcut<=0) call deallocate_local_zone_descriptors(Lzde, subname)  

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)
  call deallocate_comms(commse,subname)
  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)


end subroutine extract_potential_for_spectra
