!> @file
!!  Routines to do diagonalisation with Davidson algorithm
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Naive subroutine which performs a direct minimization of the energy 
!!  for a given hamiltonian
subroutine direct_minimization(iproc,nproc,in,at,nvirt,rxyz,rhopot,nlpsp, &
     pkernel,dpcom,GPU,KSwfn,VTwfn)
   use module_base
   use module_types
   use module_interfaces, except_this_one => direct_minimization
   use module_xc
   use yaml_output
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc,nvirt
   type(input_variables), intent(in) :: in
   type(atoms_data), intent(in) :: at
   type(DFT_PSP_projectors), intent(inout) :: nlpsp
   type(denspot_distribution), intent(in) :: dpcom
   type(DFT_wavefunction), intent(inout) :: KSwfn,VTwfn
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   type(coulomb_operator), intent(in) :: pkernel
   real(dp), dimension(*), intent(in), target :: rhopot
   type(GPU_pointers), intent(inout) :: GPU
   !local variables
   character(len=*), parameter :: subname='direct_minimization'
   logical :: msg,exctX,occorbs,endloop !extended output
   integer :: nrhodim,i3rho_add !n(c) occnorb, occnorbu, occnorbd
   integer :: i_stat,i_all,iter,ikpt,idsx_actual_before,ndiis_sd_sw
   real(gp) :: gnrm,gnrm_zero
   type(energy_terms) :: energs
   real(wp), dimension(:), pointer :: psiw,psirocc,pot

   !supplementary messages
   msg=.false.

   !logical flag which control to othogonalise wrt the occupied orbitals or not
   if (KSwfn%orbs%nkpts /= VTwfn%orbs%nkpts) then
      occorbs=.false.
   else
      occorbs=.true.
      do ikpt = 1, KSwfn%orbs%nkpts
         if (abs(maxval(KSwfn%orbs%kpts(:,ikpt) - VTwfn%orbs%kpts(:,ikpt))) > 1.d-6) then
            occorbs=.false.
            exit
         end if
      end do
   end if
   !n(c) if (occorbs) then
   !n(c)   occnorb = 0
   !n(c)   occnorbu = 0
   !n(c)   occnorbd = 0
   !n(c) else
   !n(c)   occnorb = orbs%norb
   !n(c)   occnorbu = orbs%norbu
   !n(c)   occnorbd = orbs%norbd
   !n(c) end if

   !in the GPU case, the wavefunction should be copied to the card 
   !at each HamiltonianApplication
   !rebind the GPU pointers to the orbsv structure
   if (GPUconv) then
      call free_gpu(GPU,KSwfn%orbs%norbp)
      call prepare_gpu_for_locham(VTwfn%Lzd%Glr%d%n1,VTwfn%Lzd%Glr%d%n2,VTwfn%Lzd%Glr%d%n3,in%nspin,&
           VTwfn%Lzd%hgrids(1), VTwfn%Lzd%hgrids(2), VTwfn%Lzd%hgrids(3),VTwfn%Lzd%Glr%wfd,VTwfn%orbs,GPU)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,KSwfn%orbs,in%nspin)    
      call allocate_data_OCL(VTwfn%Lzd%Glr%d%n1,VTwfn%Lzd%Glr%d%n2,VTwfn%Lzd%Glr%d%n3,at%astruct%geocode,&
         &   in%nspin,VTwfn%Lzd%Glr%wfd,VTwfn%orbs,GPU)
      if (iproc == 0) call yaml_map('GPU data allocated',.true.)
      !if (iproc == 0) write(*,*) 'GPU data allocated'
   end if

   GPU%full_locham=.true.
   !verify whether the calculation of the exact exchange term
   !should be performed
   energs%eexctX=0.0_gp
   exctX = xc_exctXfac() /= 0.0_gp
   if (in%exctxpar == 'OP2P') energs%eexctX = UNINITIALIZED(1.0_gp)
   !check the size of the rhopot array related to NK SIC
   nrhodim=in%nspin
   i3rho_add=0
   if (in%SIC%approach=='NK') then
      nrhodim=2*nrhodim
      i3rho_add=VTwfn%Lzd%Glr%d%n1i*VTwfn%Lzd%Glr%d%n2i*dpcom%nscatterarr(iproc,4)+1
   end if

   if(iproc==0) then
      call yaml_comment('Iterative subspace diagonalization of virtual orbitals (Direct Minimization)',hfill='-')
      !write(*,'(1x,a)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      !write(*,'(1x,a)') "Iterative subspace diagonalization of virtual orbitals (Direct Minimization)."
   end if


   !before transposition, create the array of the occupied
   !wavefunctions in real space, for exact exchange calculations
   !still the exact exchange with occorbs=.false. has to be verified
   if (exctX) then
      allocate(psirocc(max(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*KSwfn%orbs%norbp,&
         &   dpcom%ngatherarr(0,1)*KSwfn%orbs%norb),1)+ndebug),stat=i_stat)
      call memocc(i_stat,psirocc,'psirocc',subname)

      call prepare_psirocc(iproc,nproc,KSwfn%Lzd%Glr,KSwfn%orbs,dpcom%nscatterarr(iproc,2),dpcom%ngatherarr(0,1),KSwfn%psi,psirocc)
   else if (in%SIC%approach=='NK') then
      allocate(psirocc(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*2*KSwfn%orbs%nspin+ndebug),stat=i_stat)
      call memocc(i_stat,psirocc,'psirocc',subname)
   else
      nullify(psirocc)
   end if

   !n2virt=2*KSwfn%orbsv%norb! the dimension of the subspace

   if (occorbs) then
      !disassociate work array for transposition in serial
      if (nproc > 1) then
         allocate(psiw(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp)+ndebug),stat=i_stat)
         call memocc(i_stat,psiw,'psiw',subname)
      else
         psiw => null()
      endif

      !transpose the wavefunction psi 
      call transpose_v(iproc,nproc,KSwfn%orbs,KSwfn%lzd%glr%wfd,KSwfn%comms,KSwfn%psi(1),psiw(1))

      if (nproc > 1) then
         i_all=-product(shape(psiw))*kind(psiw)
         deallocate(psiw,stat=i_stat)
         call memocc(i_stat,i_all,'psiw',subname)
      end if
   end if

   allocate(VTwfn%orbs%eval(VTwfn%orbs%norb*VTwfn%orbs%nkpts+ndebug),stat=i_stat)
   call memocc(i_stat,VTwfn%orbs%eval,'eval',subname)

   VTwfn%orbs%eval(1:VTwfn%orbs%norb*VTwfn%orbs%nkpts)=-0.5d0

   !prepare the v array starting from a set of gaussians
   call psivirt_from_gaussians(iproc,nproc,at,VTwfn%orbs,VTwfn%Lzd,VTwfn%comms,rxyz,&
        VTwfn%Lzd%hgrids(1),VTwfn%Lzd%hgrids(2),VTwfn%Lzd%hgrids(3),in%nspin,&
        VTwfn%psi)

   !if(iproc==0) call yaml_map('Orthogonality to occupied psi',.true.)
   !if(iproc==0) write(*,'(1x,a)',advance="no") "Orthogonality to occupied psi..."

   !project v such that they are orthogonal to all occupied psi
   !Orthogonalize before and afterwards.

   !here nvirte=VTwfn%orbs%norb
   !     nvirtep=VTwfn%orbs%norbp

   !this is the same also in serial
   call orthogonalize(iproc,nproc,VTwfn%orbs,VTwfn%comms,VTwfn%psi,in%orthpar)

   if (occorbs) then
      call orthon_virt_occup(iproc,nproc,KSwfn%orbs,VTwfn%orbs,KSwfn%comms,VTwfn%comms,KSwfn%psi,VTwfn%psi,msg)
      !and orthonormalize them using "gram schmidt"  (conserve orthogonality to psi)
      call orthogonalize(iproc,nproc,VTwfn%orbs,VTwfn%comms,VTwfn%psi,in%orthpar)
   end if

   !retranspose v
   if(nproc > 1)then
      !reallocate the work array with the good size
      allocate(psiw(max(VTwfn%orbs%npsidim_orbs,VTwfn%orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,psiw,'psiw',subname)
   end if

   call untranspose_v(iproc,nproc,VTwfn%orbs,VTwfn%Lzd%Glr%wfd,VTwfn%comms,VTwfn%psi(1),psiw(1))

   ! 1st Hamilton application on psivirt
   !if(iproc==0)write(*,'(1x,a)')"done."

   allocate(VTwfn%hpsi(max(VTwfn%orbs%npsidim_orbs,VTwfn%orbs%npsidim_comp)+ndebug),stat=i_stat)
   call memocc(i_stat,VTwfn%hpsi,'VTwfn%hpsi',subname)
   if (nproc > 1) then
      allocate(VTwfn%psit(max(VTwfn%orbs%npsidim_orbs,VTwfn%orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,VTwfn%psit,'VTwfn%psit',subname)
      !transpose the psivirt 
      call transpose_v(iproc,nproc,VTwfn%orbs,VTwfn%lzd%glr%wfd,VTwfn%comms,VTwfn%psi(1),psiw(1),out_add=VTwfn%psit(1))
   else
      nullify(VTwfn%psit)
   end if

   !allocate the potential in the full box
   call full_local_potential(iproc,nproc,VTwfn%orbs,VTwfn%Lzd,0,dpcom,rhopot,pot)
   !iproc,nproc,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,2),&
   !     Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
   !     in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,1)*nrhodim,i3rho_add,&
   !     VTwfn%orbs,Lzd,0,ngatherarr,rhopot,pot)

   call local_potential_dimensions(iproc,VTwfn%Lzd,VTwfn%orbs,dpcom%ngatherarr(0,1))

   !in the case of NK SIC, put the total density in the psirocc pointer, so that it could be reused for building the 
   !Hamiltonian Application
   if (in%SIC%approach=='NK') then
      !put the wxd term in the psirocc array
      !!$     call NK_SIC_potential(Lzd%Glr,orbs,in%ixc,0.5_gp,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,&
      !!$          psi,pot,eSIC_DC,wxdsave=psirocc)
      !put the density in the *second* part of psirocc (off diangonal term presence should be verified still)
      call vcopy(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*KSwfn%orbs%nspin,&
           pot(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*KSwfn%orbs%nspin+1),1,&
         &   psirocc(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*KSwfn%orbs%nspin+1),1)
      call to_zero(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*KSwfn%orbs%nspin,psirocc(1))
   end if


   !-----------starting point of the routine of direct minimisation

   ! allocate arrays necessary for DIIS convergence acceleration
   !the allocation with npsidim is not necessary here since DIIS arrays
   !are always calculated in the transpsed form
   call allocate_diis_objects(in%idsx,in%alphadiis,sum(VTwfn%comms%ncntt(0:nproc-1)),&
      &   VTwfn%orbs%nkptsp,VTwfn%orbs%nspinor,VTwfn%diis,subname)  
   !print *,'check',in%idsx,sum(VTwfn%comms%ncntt(0:nproc-1)),VTwfn%orbs%nkptsp

   energs%eKS=1.d10
   gnrm=1.d10
   gnrm_zero=0.0_gp
!!$   ekin_sum=0.d0 
!!$   epot_sum=0.d0 
!!$   eproj_sum=0.d0
!!$   eSIC_DC=0.0_gp

   !number of switching betweed DIIS and SD during self-consistent loop
   ndiis_sd_sw=0
   !previous value of idsx_actual to control if switching has appeared
   idsx_actual_before=VTwfn%diis%idsx

   if (iproc == 0) call yaml_open_sequence('Optimization of virtual orbitals')

   wfn_loop: do iter=1,in%itermax+100

      if (iproc == 0 .and. verbose > 0) then 
         call yaml_comment('iter=' // trim(yaml_toa(iter)),hfill='-')
         call yaml_sequence(advance='no')
         call yaml_open_map(flow=.true.)
         !write( *,'(1x,a,i0)') repeat('~',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
      endif
      !control whether the minimisation iterations ended
      endloop= gnrm <= in%gnrm_cv .or. iter == in%itermax+100

      !control how many times the DIIS has switched into SD
      if (VTwfn%diis%idsx /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

      !leave SD if the DIIS did not work the second time
      if (ndiis_sd_sw > 1) then
         VTwfn%diis%switchSD=.false.
      end if

      !terminate SCF loop if forced to switch more than once from DIIS to SD
      !endloop=endloop .or. ndiis_sd_sw > 2

      call FullHamiltonianApplication(iproc,nproc,at,VTwfn%orbs,rxyz,&
           VTwfn%Lzd,nlpsp,VTwfn%confdatarr,dpcom%ngatherarr,pot,VTwfn%psi,VTwfn%hpsi,&
           energs,in%SIC,GPU,&
           pkernel,KSwfn%orbs,psirocc)

      energs%ebs=energs%ekin+energs%epot+energs%eproj
      !n(c) energy_old=energy
      energs%eKS=energs%ebs-energs%eexctX

      !check for convergence or whether max. numb. of iterations exceeded
      if (endloop) then 
         if (iproc == 0) then 
            !if (verbose > 1) call yaml_map('Minimization iterations required',iter)
            call write_energies(iter,0,energs,gnrm,0.d0,' ')
            call yaml_close_map()
            call yaml_comment('End of Virtual Wavefunction Optimisation',hfill='-')
            if (VTwfn%diis%energy > VTwfn%diis%energy_min) then
               call yaml_warning('Found an energy value lower than the FINAL energy, delta' // &
                    & trim(yaml_toa(VTwfn%diis%energy-VTwfn%diis%energy_min,fmt='(1pe9.2)')))
            end if
         end if
         exit wfn_loop 
      endif

      !evaluate the functional of the wavefucntions and put it into the diis structure
      !the energy values should be printed out here
      call total_energies(energs, iter, iproc)
      call calculate_energy_and_gradient(iter,iproc,nproc,GPU,in%ncong,in%iscf,energs,&
           VTwfn,gnrm,gnrm_zero)

      !control the previous value of idsx_actual
      idsx_actual_before=VTwfn%diis%idsx

      call hpsitopsi(iproc,nproc,iter,in%idsx,VTwfn,at,nlpsp)

      if (occorbs) then
         !if this is true the transposition for psivirt which is done in hpsitopsi
         !is useless, but we leave it for simplicity

         if (nproc == 1) VTwfn%psit => VTwfn%psi
         !project psivirt such that they are orthogonal to all occupied psi
         call orthon_virt_occup(iproc,nproc,KSwfn%orbs,VTwfn%orbs,KSwfn%comms,VTwfn%comms,KSwfn%psi,VTwfn%psit,msg)
         call orthogonalize(iproc,nproc,VTwfn%orbs,VTwfn%comms,VTwfn%psit,in%orthpar)
         !retranspose the psivirt
         call untranspose_v(iproc,nproc,VTwfn%orbs,VTwfn%Lzd%Glr%wfd,VTwfn%comms,VTwfn%psit(1),&
            &   psiw(1),out_add=VTwfn%psi(1))
      end if

      if (iproc == 0) call yaml_close_map()

   end do wfn_loop

   if (iproc == 0) then
      call yaml_close_sequence() !wfn iterations
      if (iter == in%itermax+100) then
         call yaml_warning('No convergence within the allowed number of minimization steps')
      else if (verbose > 1) then
         call yaml_map('Minimization iterations required',iter)
      end if
   end if

   !deallocate real array of wavefunctions
   if(exctX .or. in%SIC%approach=='NK')then
      i_all=-product(shape(psirocc))*kind(psirocc)
      deallocate(psirocc,stat=i_stat)
      call memocc(i_stat,i_all,'psirocc',subname)
   end if

   call deallocate_diis_objects(VTwfn%diis,subname)

   !this deallocates also hpsivirt and psitvirt
   call last_orthon(iproc,nproc,iter,VTwfn,energs%evsum)

   !resize work array before final transposition
   if(nproc > 1)then
      i_all=-product(shape(psiw))*kind(psiw)
      deallocate(psiw,stat=i_stat)
      call memocc(i_stat,i_all,'psiw',subname)

      allocate(psiw(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,psiw,'psiw',subname)
   end if

   call untranspose_v(iproc,nproc,KSwfn%orbs,KSwfn%Lzd%Glr%wfd,KSwfn%comms,KSwfn%psi(1),psiw(1))

   if(nproc > 1) then
      i_all=-product(shape(psiw))*kind(psiw)
      deallocate(psiw,stat=i_stat)
      call memocc(i_stat,i_all,'psiw',subname)
   end if
   !!!!! end point of the direct minimisation procedure

   !deallocate potential
   call free_full_potential(dpcom%mpi_env%nproc,0,pot,subname)

   if (GPUconv) then
      call free_gpu(GPU,VTwfn%orbs%norbp)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,VTwfn%orbs,in%nspin)
   end if

   call calculate_HOMO_LUMO_gap(iproc,KSwfn%orbs,VTwfn%orbs)

   !the plotting should be added here (perhaps build a common routine?)
   call write_eigen_objects(iproc,occorbs,in%nspin,nvirt,in%nplot,VTwfn%Lzd%hgrids(1),VTwfn%Lzd%hgrids(2),VTwfn%Lzd%hgrids(3),&
        at,rxyz,KSwfn%Lzd%Glr,KSwfn%orbs,VTwfn%orbs,KSwfn%psi,VTwfn%psi,in%output_wf_format)

 END SUBROUTINE direct_minimization


!> Davidsons method for iterative diagonalization of virtual Kohn Sham orbitals
!!   under orthogonality constraints to occupied orbitals psi. The nvirt input
!!   variable gives the number of unoccupied orbitals for which the exit criterion
!!   for the gradients norm holds. nvirte = norbe - norb >= nvirt is the number of
!!   virtual orbitals processed by the method. The dimension of the subspace for
!!   diagonalization is 2*nvirte = n2virt
!!                                                                   Alex Willand
!!   Algorithm
!!   _________
!!   (parallel)
!!
!!   (transpose psi, v is already transposed)\n
!!   orthogonality of v to psi\n
!!   orthogonalize v\n
!!   (retranspose v)\n
!!   Hamilton(v) --> hv\n
!!   transpose v and hv\n
!!   Rayleigh quotients  e\n
!!   do\n
!!      gradients g= e*v -hv\n
!!      exit condition gnrm\n
!!      orthogonality of g to psi\n
!!      (retranspose g)\n
!!      preconditioning of g\n
!!      (transpose g again)\n
!!      orthogonality of g to psi\n
!!      (retranspose g)\n
!!      Hamilton(g) --> hg\n
!!      (transpose g and hg)\n
!!      subspace matrices H and S\n
!!      DSYGV(H,e,S)  --> H\n
!!      update v with eigenvectors H\n
!!      orthogonality of v to psi\n
!!      orthogonalize v\n
!!      (retranspose v)\n
!!      Hamilton(v) --> hv\n
!!      (transpose v and hv)\n
!!   end do\n
!!   (retranspose v and psi)\n
subroutine davidson(iproc,nproc,in,at,&
     & orbs,orbsv,nvirt,Lzd,comms,commsv,&
     & rxyz,rhopot,nlpsp,pkernel,psi,v,dpcom,GPU)
   use module_base
   use module_types
   use module_interfaces, except_this_one => davidson
   use module_xc
   use yaml_output
   use communications_base, only: comms_cubic
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc
   integer, intent(in) :: nvirt
   type(input_variables), intent(in) :: in
   type(atoms_data), intent(in) :: at
   type(DFT_PSP_projectors), intent(inout) :: nlpsp
   type(local_zone_descriptors), intent(inout) :: Lzd
   type(orbitals_data), intent(inout) :: orbs !<could be modify in calculate_HOMO_LUMO_gap
   type(comms_cubic), intent(in) :: comms, commsv
   type(denspot_distribution), intent(in) :: dpcom
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   type(coulomb_operator), intent(in) :: pkernel
   real(dp), dimension(*), intent(in) :: rhopot
   type(orbitals_data), intent(inout) :: orbsv
   type(GPU_pointers), intent(inout) :: GPU
   real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
   !v, that is psivirt, is transposed on input and direct on output
   !local variables
   character(len=*), parameter :: subname='davidson',print_precise='1pe20.12',print_rough='1pe12.4 '
   character(len=8) :: prteigu,prteigd !format for eigenvalues printing
   logical :: msg,exctX,occorbs !extended output
   integer :: nrhodim,i3rho_add !n(c) occnorb, occnorbu, occnorbd
   integer :: ierr,i_stat,i_all,iorb,jorb,iter,nwork,norb,nspinor
   integer :: ise,j,ispsi,ikpt,ikptp,nvctrp,ncplx,ncomp,norbs,ispin,ish1,ish2,nspin
   real(gp) :: tt,gnrm,gnrm_fake
   integer, dimension(:,:), allocatable :: ndimovrlp
   real(wp), dimension(:), allocatable :: work,work_rp,hamovr
   real(wp), dimension(:), allocatable :: hv,g,hg,ew
   real(wp), dimension(:,:,:), allocatable :: e
   real(wp), dimension(:), pointer :: psiw,psirocc,pot
   type(confpot_data), dimension(:), allocatable :: confdatarr
   type(energy_terms) :: energs

   !logical flag which control to othogonalise wrt the occupied orbitals or not
   if (orbs%nkpts /= orbsv%nkpts) then
      occorbs=.false.
   else
      occorbs=.true.
      do ikpt = 1, orbs%nkpts
         if (abs(maxval(orbs%kpts(:,ikpt) - orbsv%kpts(:,ikpt))) > 1.d-6) then
            occorbs=.false.
            exit
         end if
      end do
   end if
   !n(c) if (occorbs) then
   !n(c)   occnorb = 0
   !n(c)   occnorbu = 0
   !n(c)   occnorbd = 0
   !n(c) else
   !n(c)   occnorb = orbs%norb
   !n(c)   occnorbu = orbs%norbu
   !n(c)   occnorbd = orbs%norbd
   !n(c) end if

   !in the GPU case, the wavefunction should be copied to the card 
   !at each HamiltonianApplication
   !rebind the GPU pointers to the orbsv structure
   if (GPUconv) then
      call free_gpu(GPU,orbs%norbp)
      call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,in%nspin,&
           Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),Lzd%Glr%wfd,orbsv,GPU)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,orbs,in%nspin)    
      call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
           in%nspin,Lzd%Glr%wfd,orbsv,GPU)
   end if

   GPU%full_locham=.true.
   !verify whether the calculation of the exact exchange term
   !should be performed
   energs%eexctX=0.0_gp
   exctX = xc_exctXfac() /= 0.0_gp
   if (in%exctxpar == 'OP2P') energs%eexctX = UNINITIALIZED(1.0_gp)

   !check the size of the rhopot array related to NK SIC
   nrhodim=in%nspin
   i3rho_add=0
   if (trim(in%SIC%approach)=='NK') then
      nrhodim=2*nrhodim
      i3rho_add=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*dpcom%nscatterarr(iproc,4)+1
   end if

   !last index of e and hamovr are for mpi_alLzd%Glreduce. 
   !e (eigenvalues) is also used as 2 work arrays

   msg=verbose > 2 .and. iproc ==0! no extended output
   !msg =(iproc==0)!extended output

   if (iproc==0) then
      call yaml_comment('Iterative subspace diagonalization of virtual orbitals',hfill='-')
      !write(*,'(1x,a)')"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      !write(*,'(1x,a)')"Iterative subspace diagonalization of virtual orbitals."
   end if

   !if(msg)write(*,*)'shape(v)',shape(v),'size(v)',size(v)


   !before transposition, create the array of the occupied
   !wavefunctions in real space, for exact exchange calculations
   if (exctX) then
      allocate(psirocc(max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%norbp,&
         &   dpcom%ngatherarr(0,1)*orbs%norb),1)+ndebug),stat=i_stat)
      call memocc(i_stat,psirocc,'psirocc',subname)

      call prepare_psirocc(iproc,nproc,Lzd%Glr,orbs,dpcom%nscatterarr(iproc,2),dpcom%ngatherarr(0,1),psi,psirocc)
   else if (in%SIC%approach=='NK') then
      allocate(psirocc(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*2*orbs%nspin+ndebug),stat=i_stat)
      call memocc(i_stat,psirocc,'psirocc',subname)
   else
      nullify(psirocc)
   end if

   !n2virt=2*orbsv%norb! the dimension of the subspace

   if (occorbs) then
      !disassociate work array for transposition in serial
      if (nproc > 1) then
         allocate(psiw(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
         call memocc(i_stat,psiw,'psiw',subname)
      else
         psiw => null()
      endif

      !transpose the wavefunction psi 
      call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),psiw(1))

      if (nproc > 1) then
         i_all=-product(shape(psiw))*kind(psiw)
         deallocate(psiw,stat=i_stat)
         call memocc(i_stat,i_all,'psiw',subname)
      end if
   end if

   allocate(orbsv%eval(orbsv%norb*orbsv%nkpts+ndebug),stat=i_stat)
   call memocc(i_stat,orbsv%eval,'eval',subname)

   orbsv%eval(1:orbsv%norb*orbsv%nkpts)=-0.5d0

   !prepare the v array starting from a set of gaussians
   call psivirt_from_gaussians(iproc,nproc,at,orbsv,Lzd,commsv,rxyz,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),in%nspin,v)

   !if(iproc==0) call yaml_open_map('Orthogonality to occupied psi',flow=.true.)
   !if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."

   !project v such that they are orthogonal to all occupied psi
   !Orthogonalize before and afterwards.

   !here nvirte=orbsv%norb
   !     nvirtep=orbsv%norbp

   !this is the same also in serial
   call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)

   if (occorbs) then
      call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
      !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
      call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)
   end if

   !retranspose v
   if(nproc > 1)then
      !reallocate the work array with the good size
      allocate(psiw(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,psiw,'psiw',subname)
   end if

   call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,v(1),psiw(1))

   ! 1st Hamilton application on psivirt
   !if(iproc==0) then
   !   call yaml_map('first','done')
   !   call yaml_close_map()
   !end if
   !if(iproc==0)vwrite(*,'(1x,a)',advance="no")"done. first "

   allocate(hv(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
   call memocc(i_stat,hv,'hv',subname)

   !allocate the potential in the full box
   call full_local_potential(iproc,nproc,orbsv,Lzd,0,dpcom,rhopot,pot)
   !(iproc,nproc,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,2),&
   !     Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
   !     in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,1)*nrhodim,i3rho_add,&
   !     orbsv,Lzd,0,ngatherarr,rhopot,pot)

   call local_potential_dimensions(iproc,Lzd,orbsv,dpcom%ngatherarr(0,1))
   allocate(confdatarr(orbsv%norbp))
   call default_confinement_data(confdatarr,orbsv%norbp)


   !in the case of NK SIC, put the total density in the psirocc pointer, so that it could be reused for building the 
   !Hamiltonian Application
   if (in%SIC%approach=='NK') then
      !put the density in the *second* part of psirocc
      call vcopy(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin,&
           pot(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+1),1,&
           psirocc(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+1),1)
      call to_zero(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin,psirocc(1))
      !!$     call NK_SIC_potential(Lzd%Glr,orbs,in%SIC%ixc,in%SIC%fref,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,&
      !!$          psi,pot(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+1:),eSIC_DC,wxdsave=psirocc)
   end if
   !experimental: add confining potential to the hamiltonian
      !should already be guaranteed by the crmult terms
      !call add_confining_potential(Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,orbs%nspin,1.e-10_gp,1.e-14_gp,-0.5_gp,&
      !     pot(1),pot(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+1))
   

   !experimental: add parabolic potential to the hamiltonian
   !call add_parabolic_potential(at%astruct%geocode,at%astruct%nat,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,12.0_gp,rxyz,pot)

   call FullHamiltonianApplication(iproc,nproc,at,orbsv,rxyz,&
        Lzd,nlpsp,confdatarr,dpcom%ngatherarr,pot,v,hv,&
        energs,in%SIC,GPU,&
        pkernel,orbs,psirocc)

   !if(iproc==0)write(*,'(1x,a)',advance="no")"done. Rayleigh quotients..."

   allocate(e(orbsv%norb,orbsv%nkpts,2+ndebug),stat=i_stat)
   call memocc(i_stat,e,'e',subname)

   !transpose  v and hv
   call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,v(1),psiw(1))
   call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,hv(1),psiw(1))

   call timing(iproc,'Davidson      ','ON')
   !Timing excludes transposition, hamilton application and preconditioning
   call to_zero(orbsv%norb*2*orbsv%nkpts,e)
   ! Rayleigh quotients.

   !probably this loop can be rewritten using GEMMs
   ispsi=1
   do ikptp=1,orbsv%nkptsp
      ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
      nvctrp=commsv%nvctr_par(iproc,ikpt)
      if (nvctrp == 0) cycle

      nspinor=orbsv%nspinor
      do iorb=1,orbsv%norb ! temporary variables 
         !for complex wavefunctions the diagonal is always real
         e(iorb,ikpt,1)= dot(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1,&
            &   hv(ispsi+nvctrp*nspinor*(iorb-1)),1)          != <psi|H|psi> 
         e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1)**2   != <psi|psi> 
      end do
      ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
   end do

   if(nproc > 1)then
      !sum up the contributions of nproc sets with 
      !commsv%nvctr_par(iproc,1) wavelet coefficients each
      call mpiallred(e(1,1,1),2*orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

   end if

   !if(iproc==0)write(*,'(1x,a)')"done."
   !if(iproc==0)write(*,'(1x,a)')"      1-sqnorm   Rayleigh quotient"
   if(iproc==0) call yaml_open_sequence('L2 Norm - 1 and Rayleigh quotient (Davidson)')

   do ikpt=1,orbsv%nkpts
      !if (orbsv%nkpts > 1 .and. iproc == 0) then
      if (iproc == 0) then
         call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)')))&
              // ' BZ coord. = ' // &
              trim(yaml_toa(orbs%kpts(:, ikpt),fmt='(f12.6)')))
         !call yaml_sequence(advance='no')
         !write(*,"(1x,A,I3.3,A,3F12.6)") " Kpt #", ikpt, " BZ coord. = ", orbsv%kpts(:, ikpt)
      end if
      do iorb=1,orbsv%norb
         !e(:,1,1) = <psi|H|psi> / <psi|psi>
         e(iorb,ikpt,1)=e(iorb,ikpt,1)/e(iorb,ikpt,2)
         if (iproc == 0) then
            call yaml_sequence(trim(yaml_toa((/ 1.0_gp-e(iorb,ikpt,2),e(iorb,ikpt,1) /),fmt='(1pe13.6)')),&
                               advance='no')
            !call yaml_map('Orbitals',&
            !     ,advance='no')
            call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
         end if
         !if(iproc==0) write(*,'(1x,i3,1x,1pe13.6,1x,1pe12.5)') iorb,1.0_gp-e(iorb,ikpt,2),e(iorb,ikpt,1)
      end do
   end do
   if(iproc==0) call yaml_close_sequence()

   !if(msg)then
   !write(*,*)"******** transposed v,hv 1st elements"
   !do iorb=1,10
   !  write(*,*)v(iorb,1),hv(iorb,1)
   !end do
   !write(*,*)"**********"
   !end if

   !calculate the dimension of the overlap matrix for each k-point
   if (orbsv%norbd > 0) then
      nspin=2
   else
      nspin=1
   end if

   !number of components for the overlap matrix in wp-kind real numbers

   allocate(ndimovrlp(nspin,0:orbsv%nkpts+ndebug),stat=i_stat)
   call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

   call dimension_ovrlp(nspin,orbsv,ndimovrlp)

   !the dimension should be chosen with the max between k-points
   !allocate(hamovr(n2virt,n2virt,2,orbsv%nkpts+ndebug),stat=i_stat)
   allocate(hamovr(8*ndimovrlp(nspin,orbsv%nkpts)+ndebug),stat=i_stat)
   call memocc(i_stat,hamovr,'hamovr',subname)

   !put to zero all the k-points which are not needed
   call to_zero(8*ndimovrlp(nspin,orbsv%nkpts),hamovr)

   if (orbsv%nspinor > 1) then
      ncplx=2
   else
      ncplx=1
   end if
   nwork=max(10,16*orbsv%norb)
   allocate(work(ncplx*nwork+ndebug),stat=i_stat)
   call memocc(i_stat,work,'work',subname)

   allocate(ew(2*orbsv%norb+ndebug),stat=i_stat)
   call memocc(i_stat,ew,'ew',subname)


   !itermax=... use the input variable instead
   iter=1
   davidson_loop: do 

      if(iproc==0) call yaml_comment(' iter= ' // trim(yaml_toa(iter)),hfill='-')
      if(msg) call yaml_open_sequence('squared norm of the (nvirt) gradients')
      !if(iproc==0) write( *,'(1x,a,i0)') repeat('~',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
      !if(msg) write(*,'(1x,a)')"squared norm of the (nvirt) gradients"

      allocate(g(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,g,'g',subname)

      call vcopy(max(orbsv%npsidim_orbs,orbsv%npsidim_comp),hv(1),1,g(1),1)! don't overwrite hv

      call to_zero(orbsv%norb*orbsv%nkpts,e(1,1,2))
      !also these operations are presumably GEMMs
      !here we should add the ncomp term for non-collinear case
      ispsi=1
      do ikptp=1,orbsv%nkptsp
         ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
         nvctrp=commsv%nvctr_par(iproc,ikpt)
         if (nvctrp == 0) cycle

         nspinor=orbsv%nspinor
         do iorb=1,orbsv%norb
            !gradient = hv-e*v
            call axpy(nvctrp*nspinor,-e(iorb,ikpt,1),v(ispsi+nvctrp*nspinor*(iorb-1)),1,&
               &   g(ispsi+nvctrp*nspinor*(iorb-1)),1)

            !local contribution to the square norm
            e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)**2
         end do
         ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
      end do

      if(nproc > 1)then
         !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
         call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
      end if

      gnrm=0._dp
      do ikpt=1,orbsv%nkpts
         if (msg) then
            call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))) // ' BZ coord. = ' // &
            & trim(yaml_toa(orbs%kpts(:, ikpt),fmt='(f12.6)')))
            call yaml_sequence(advance='no')
            call yaml_open_map('tt',flow=.true.)
         end if
         do iorb=1,orbsv%norb
            tt=real(e(iorb,ikpt,2)*orbsv%kwgts(ikpt),dp)
            if (msg) call yaml_map('Orbital No.'//trim(yaml_toa(iorb)), tt, fmt='(1pe21.14)')
            !if(msg) write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
            if (iorb <= nvirt) gnrm=gnrm+tt
         end do
         if (nspin == 2) then
            do iorb=1,orbsv%norbu
               tt=real(e(iorb+orbsv%norbu,ikpt,2)*orbsv%kwgts(ikpt),dp)
               if (msg) call yaml_map('Orbital No.'//trim(yaml_toa(iorb+orbsv%norbu)), tt, fmt='(1pe21.14)')
               !if(msg) write(*,'(1x,i3,1x,1pe21.14)')iorb+orbsv%norbu,tt
               if (iorb <= nvirt) gnrm=gnrm+tt
            end do
         end if
         if (msg) call yaml_close_map()
      end do
      if (msg) call yaml_close_sequence()

      !the gnrm defined should be the average of the active gnrms
      gnrm=dsqrt(gnrm/real(nvirt,dp))

      if (iproc == 0) then
         call yaml_open_map('Gradient Norm',flow=.true.)
         call yaml_map('Value',gnrm,fmt='(1pe12.5)')
         call yaml_map('Exit criterion',in%gnrm_cv,fmt='(1pe9.2)')
         call yaml_close_map()
         !write(*,'(1x,a,2(1x,1pe12.5))') "|gradient|=gnrm and exit criterion ",gnrm,in%gnrm_cv
      end if

      if(gnrm < in%gnrm_cv) then
         i_all=-product(shape(g))*kind(g)
         deallocate(g,stat=i_stat)
         call memocc(i_stat,i_all,'g',subname)
         exit davidson_loop! iteration loop
      end if
      call timing(iproc,'Davidson      ','OF')

      !if(iproc==0) call yaml_map('Orthogonality of gradients to occupied psi...',.true.)
      !if(iproc==0) write(*,'(1x,a)',advance="no") "Orthogonality of gradients to occupied psi..."

      !project g such that they are orthogonal to all occupied psi. 
      !Gradients do not need orthogonality.
      if (occorbs) then
         call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
      end if

      call timing(iproc,'Davidson      ','ON')
      !if(iproc==0)write(*,'(1x,a)',advance="no")"done."

      if(msg) then
         call to_zero(orbsv%norb*orbsv%nkpts,e(1,1,2))
         call yaml_open_sequence('squared norm of all gradients after projection')
         !write(*,'(1x,a)')"squared norm of all gradients after projection"
         ispsi=1
         do ikptp=1,orbsv%nkptsp
            ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
            nvctrp=commsv%nvctr_par(iproc,ikpt)
            if (nvctrp == 0) cycle

            nspinor=orbsv%nspinor
            do iorb=1,orbsv%norb
               e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)**2
            end do

            ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
         end do

         if(nproc > 1)then
            !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
            call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
         end if

         gnrm=0._dp
         do ikpt=1,orbsv%nkpts
            call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))) // ' BZ coord. = ' // &
            & trim(yaml_toa(orbs%kpts(:, ikpt),fmt='(f12.6)')))
            call yaml_sequence(advance='no')
            call yaml_open_map('tt',flow=.true.)
            do iorb=1,orbsv%norb
               tt=real(e(iorb,ikpt,2)*orbsv%kwgts(ikpt),dp)
               call yaml_map('Orbital No.'//trim(yaml_toa(iorb)), tt, fmt='(1pe21.14)')
               !write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
               gnrm=gnrm+tt
            end do
            if (nspin == 2) then
               do iorb=1,orbsv%norb
                  tt=real(e(iorb+orbsv%norbu,ikpt,2)*orbsv%kwgts(ikpt),dp)
                  call yaml_map('Orbital No.'//trim(yaml_toa(iorb+orbsv%norbu)), tt, fmt='(1pe21.14)')
                  !write(*,'(1x,i3,1x,1pe21.14)') iorb+orbsv%norbu,tt
                  gnrm=gnrm+tt
               end do
            end if
            call yaml_close_map()
         end do
         gnrm=sqrt(gnrm/real(orbsv%norb,dp))

         call yaml_map('gnrm of all',gnrm,fmt='(1pe21.14)')
         !write(*,'(1x,a,2(1x,1pe21.14))')"gnrm of all ",gnrm
         call yaml_close_sequence()
      end if

      !if (iproc==0) write(*,'(1x,a)',advance='no')'Preconditioning...'

      call timing(iproc,'Davidson      ','OF')

      !retranspose the gradient g 
      call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,g(1),psiw(1))

      ! Here the gradients norm could be calculated in the direct form instead,
      ! as it is done in hpsiortho before preconditioning. 
      ! However, this should not make a difference and is not really simpler 

      call timing(iproc,'Precondition  ','ON')

      !we fill the values of the eval for the orbitals used in the preconditioner
      do ikpt=1,orbsv%nkpts
         do iorb=1,orbsv%norb
            !write(*,*) 'iorb,e(iorb,ikpt,1)',iorb,e(iorb,ikpt,1)
            !orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=min(e(iorb,ikpt,1),-.5d0)
            orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=e(iorb,ikpt,1)
         end do
      end do
      !we use for preconditioning the eval from the lowest value of the KS wavefunctions
      !WARNING this pointer association is nor coherent: 
      !        it will give errors if orbs%norb < orbsv%norb
      !     orbsv%eval=>orbs%eval
      !     if (orbs%norb < orbsv%norb) then
      !        write(*,*)'ERROR: too many virtual orbitals demanded, cannot proceed.'
      !        write(*,*)'       change the eigenvalues array to fix this problem'
      !        stop
      !     end if

      call preconditionall(orbsv,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),in%ncong,g,gnrm_fake,gnrm_fake)

      call timing(iproc,'Precondition  ','OF')
      !if (iproc==0)write(*,'(1x,a)')'done.'

      !if(iproc==0) call yaml_map('Orthogonality of preconditioned gradients to occupied psi',.true.)
      !if(iproc==0)write(*,'(1x,a)',advance="no") "Orthogonality of preconditioned gradients to occupied psi..."

      if (occorbs) then
         !transpose  g 
         call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,g(1),psiw(1))
         !project g such that they are orthogonal to all occupied psi
         call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
         !retranspose the gradient g
         call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,g(1),psiw(1))
      end if

      !if(iproc==0)write(*,'(1x,a)')"done."

      allocate(hg(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,hg,'hg',subname)

      call FullHamiltonianApplication(iproc,nproc,at,orbsv,rxyz,&
           Lzd,nlpsp,confdatarr,dpcom%ngatherarr,pot,g,hg,&
           energs,in%SIC,GPU,&
           pkernel,orbs,psirocc)

      !transpose  g and hg
      call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,g(1),psiw(1))
      call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,hg(1),psiw(1))

      call timing(iproc,'Davidson      ','ON')
      !if(iproc==0)write(*,'(1x,a)',advance="no")"done."

      if(msg) then
         call to_zero(orbsv%norb*orbsv%nkpts,e(1,1,2))
         write(*,'(1x,a)')"Norm of all preconditioned gradients"
         ispsi=1
         do ikptp=1,orbsv%nkptsp
            ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
            nvctrp=commsv%nvctr_par(iproc,ikpt)
            if (nvctrp == 0) cycle

            nspinor=orbsv%nspinor
            do iorb=1,orbsv%norb
               e(iorb,ikpt,2)=nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)**2
            end do
            ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
         end do

         if(nproc > 1)then
            !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
            call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
         end if

         gnrm=0.0_dp
         do ikpt=1,orbsv%nkpts
            do iorb=1,orbsv%norb
               tt=real(e(iorb,ikpt,2)*orbsv%kwgts(ikpt),dp)
               if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
               gnrm=gnrm+tt
            end do
         end do
         gnrm=sqrt(gnrm/real(orbsv%norb,dp))

         write(*,'(1x,a,2(1x,1pe21.14))')"gnrm of all",gnrm
      end if

      !if(iproc==0)write(*,'(1x,a)',advance="no")"Expanding subspace matrices..."

      !                 <vi | hvj>      <vi | hgj-n>                   <vi | vj>      <vi | gj-n>
      ! hamovr(i,j,1)=                               ;  hamovr(i,j,2)=  
      !                 <gi-n | hvj>  <gi-n | hgj-n>                   <gi-n | vj>  <gi-n | gj-n>
      !put to zero all the k-points which are not needed
      call to_zero(8*ndimovrlp(nspin,orbsv%nkpts),hamovr)


      ! store upper triangular part of these matrices only
      ! therefore, element (iorb+nvirte,jorb) is transposed to (j,nvirt+iorb)
      ispsi=1
      do ikptp=1,orbsv%nkptsp
         ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)

         do ispin=1,nspin

            call orbitals_and_components(iproc,ikpt,ispin,orbsv,commsv,&
               &   nvctrp,norb,norbs,ncomp,nspinor)
            if (nvctrp == 0) cycle
            !print *,iproc,ikpt,ispin,norb,nspinor,ncplx,nvctrp,8*ndimovrlp(ispin,ikpt-1)+1,8*ndimovrlp(nspin,orbsv%nkpts)
            call Davidson_subspace_hamovr(norb,nspinor,ncplx,nvctrp,&
               &   hamovr(8*ndimovrlp(ispin,ikpt-1)+1),&
               &   v(ispsi:),g(ispsi),hv(ispsi),hg(ispsi))

            ispsi=ispsi+nvctrp*norb*nspinor
         end do
      end do

      i_all=-product(shape(hg))*kind(hg)
      deallocate(hg,stat=i_stat)
      call memocc(i_stat,i_all,'hg',subname)

      if(nproc > 1)then
         !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
         call mpiallred(hamovr(1),8*ndimovrlp(nspin,orbsv%nkpts),&
            &   MPI_SUM,bigdft_mpi%mpi_comm,ierr)
      end if

      !if(iproc==0)write(*,'(1x,a)')"done."

      ispsi=1
      do ikptp=1,orbsv%nkptsp
         ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)

         !if(msg .or. (iproc==0 .and. ikpt == 1)) write(*,'(1x,a)',advance='no')"Diagonalization..."

         ise=0
         do ispin=1,nspin
            call orbitals_and_components(iproc,ikpt,ispin,orbsv,commsv,&
               &   nvctrp,norb,norbs,ncomp,nspinor)
            if (nvctrp == 0) cycle

            if (nspinor /= 1) then
               allocate(work_rp(6*norb+1+ndebug),stat=i_stat)
               call memocc(i_stat,work_rp,'work_rp',subname)
            end if

            ish1=8*ndimovrlp(ispin,ikpt-1)+1
            ish2=8*ndimovrlp(ispin,ikpt-1)+4*norbs*norb+1

            if(msg)then
               write(*,*)"subspace matrices, upper triangular (diagonal elements first)"
               write(*,'(1x)')
               write(*,*)"subspace H "
               do iorb=1,2*norbs
                  write(*,'(100(1x,1pe12.5))')&
                     &   (hamovr(ish1-1+iorb+(jorb-1)*2*norbs),jorb=1,2*norb)
                  write(*,*)
               end do
               write(*,*)"subspace S"
               write(*,*)
               do iorb=1,2*norbs
                  write(*,'(100(1x,1pe12.5))')&
                     &   (hamovr(ish2-1+iorb+(jorb-1)*2*norbs),jorb=1,2*norb)
                  write(*,*)
               end do
            end if

            if (nspinor == 1) then
               call sygv(1,'V','U',2*norb,hamovr(ish1),2*norb,hamovr(ish2),2*norb,&
                  &   ew(1),work(1),nwork,i_stat)! Lapack GEVP
               if (i_stat /= 0) write(*,*) &
                  &   'Error in SYGV on process ',iproc,', infocode ', i_stat

            else
               call hegv(1,'V','U',2*norb,hamovr(ish1),2*norb,hamovr(ish2),2*norb,&
                  &   ew(1),work(1),nwork,work_rp(1),i_stat)! Lapack GEVP

               if (i_stat /= 0) write(*,*) &
                  &   'Error in HEGV on process ',iproc,', infocode ', i_stat

            end if

            do iorb=1,norb
               e(iorb+ise,ikpt,1)=ew(iorb)
               if (msg) e(iorb+ise,ikpt,2)=ew(iorb+norb)
            end do
            ise=norb

            if (nspinor /= 1) then
               i_all=-product(shape(work_rp))*kind(work_rp)
               deallocate(work_rp,stat=i_stat)
               call memocc(i_stat,i_all,'work_rp',subname)
            end if

            if(msg)then
               write(*,'(1x,a)')'    e(update)           e(not used)'
               do iorb=1,orbsv%norb
                  write(*,'(1x,i3,2(1pe21.14))')iorb, (e(iorb,ikpt,j),j=1,2)
               end do
               write(*,*)
               write(*,*)"and the eigenvectors are"
               write(*,*)
               do iorb=1,2*norbs
                  write(*,'(100(1x,1pe12.5))')&
                     &   (hamovr(ish1-1+iorb+(jorb-1)*2*norbs),jorb=1,2*norb)
                  write(*,*)
               end do
            end if


            !if(msg .or. (iproc==0 .and. ikpt == 1) .and. ispin==1) &
            !  & write(*,'(1x,a)',advance="no")"done. Update v with eigenvectors..."

            !!$     !Update v, that is the wavefunction, using the eigenvectors stored in hamovr(:,:,1)
            !!$     !Lets say we have 4 quarters top/bottom left/right, then
            !!$     !v = matmul(v, hamovr(topleft)  ) + matmul(g, hamovr(bottomleft)  )     needed    
            !!$     !g=  matmul(v, hamovr(topright) ) + matmul(g, hamovr(bottomright) ) not needed
            !!$     !use hv as work arrray
            !!$
            !!$     ispsi=1
            !!$     do ikptp=1,orbsv%nkptsp
            !!$        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
            !!$        nvctrp=commsv%nvctr_par(iproc,ikptp)
            !!$
            !!$        do jorb=1,orbsv%norb! v to update
            !!$           call to_zero(nvctrp,hv(ispsi+nvctrp*(jorb-1)))
            !!$           do iorb=1,orbsv%norb ! sum over v and g
            !!$              tt=hamovr(iorb,jorb,ikpt,1)
            !!$              call axpy(nvctrp,tt,v(ispsi+nvctrp*(iorb-1)),1,&
            !!$                   hv(ispsi+nvctrp*(jorb-1)),1)
            !!$              tt=hamovr(iorb+orbsv%norb,jorb,ikpt,1)
            !!$              call axpy(nvctrp,tt,g(ispsi+nvctrp*(iorb-1)),1,&
            !!$                   hv(ispsi+nvctrp*(jorb-1)),1)
            !!$           enddo
            !!$        enddo
            !!$
            !!$        call vcopy(nvctrp*orbsv%norb,hv(ispsi),1,v(ispsi),1)
            !!$
            !!$        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
            !!$     end do

            call update_psivirt(norb,nspinor,ncplx,nvctrp,&
               &   hamovr(ish1),v(ispsi:),g(ispsi),hv(ispsi))

            ispsi=ispsi+nvctrp*norb*nspinor
         end do

         if(msg .or. (iproc==0 .and. ikpt == 1)) then
            call yaml_open_sequence('Eigenvalues and eigenstate residue')
            !write(*,'(1x,a)')'done. Eigenvalues, gnrm'
            if (nspin ==1) then
               do iorb=1,orbsv%norb
                  !show the eigenvalue in full form only if it has reached convergence
                  if (sqrt(e(iorb,ikpt,2)) <= in%gnrm_cv) then
                     prteigu=print_precise
                  else
                     prteigu=print_rough
                  end if
                     call yaml_sequence(trim(yaml_toa((/ e(iorb,ikpt,1),sqrt(e(iorb,ikpt,2)) /),fmt='('//prteigu//')')),&
                          advance='no')
                     call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
                  !write(*,'(1x,i5,'//prteigu//',1pe9.2)')iorb,e(iorb,ikpt,1),sqrt(e(iorb,ikpt,2))
               end do
            else if (nspin == 2) then
               do iorb=1,min(orbsv%norbu,orbsv%norbd) !they should be equal
                  if (sqrt(e(iorb,ikpt,2)) <= in%gnrm_cv) then
                     prteigu=print_precise
                  else
                     prteigu=print_rough
                  end if
                  if (sqrt(e(iorb+orbsv%norbu,ikpt,2)) <= in%gnrm_cv) then
                     prteigd=print_precise
                  else
                     prteigd=print_rough
                  end if
                  call yaml_sequence(trim(yaml_toa((/ &
                       & e(iorb,ikpt,1),sqrt(e(iorb,ikpt,2)), &
                          & e(iorb+orbsv%norbu,ikpt,1),sqrt(e(iorb+orbsv%norbu,ikpt,2)) /),fmt='('//prteigu//')')),&
                          advance='no')
                     call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
                  !write(*,'(1x,i5,'//prteigu//',1pe9.2,t50,'//prteigd//',1pe9.2)')&
                  !   &   iorb,e(iorb,ikpt,1),sqrt(e(iorb,ikpt,2)),e(iorb+orbsv%norbu,ikpt,1),sqrt(e(iorb+orbsv%norbu,ikpt,2))
               end do
            end if
            call yaml_close_sequence()
         end if

      end do

      i_all=-product(shape(g))*kind(g)
      deallocate(g,stat=i_stat)
      call memocc(i_stat,i_all,'g',subname)

      !if(iproc==0)write(*,'(1x,a)')"done."
      !if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."
      !project v such that they are orthogonal to all occupied psi
      !Orthogonalize before and afterwards.

      call timing(iproc,'Davidson      ','OF')

      !these routines should work both in parallel or in serial
      call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)

      if (occorbs) then
         call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
         !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
         call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)
      end if

      !retranspose v
      call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,v(1),psiw(1))

      ! Hamilton application on v
      !if(iproc==0)write(*,'(1x,a)',advance="no")"done."

      call FullHamiltonianApplication(iproc,nproc,at,orbsv,rxyz,&
           Lzd,nlpsp,confdatarr,dpcom%ngatherarr,pot,v,hv,&
           energs,in%SIC,GPU,&
           pkernel,orbs,psirocc)

      !transpose  v and hv
      call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,v(1),psiw(1))
      call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,hv(1),psiw(1))

      !if(iproc==0 .and. verbose > 1) write(*,'(1x,a)')"done. "
      call timing(iproc,'Davidson      ','ON')
      iter=iter+1
      if(iter>in%itermax+100)then !an input variable should be put
         if(iproc==0)write(*,'(1x,a)')&
            &   'No convergence within the allowed number of minimization steps (itermax + 100)'
         exit davidson_loop
      end if

   end do davidson_loop

   !deallocate potential
   call free_full_potential(dpcom%mpi_env%nproc,0,pot,subname)

   i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
   deallocate(ndimovrlp,stat=i_stat)
   call memocc(i_stat,i_all,'ndimovrlp',subname)

   i_all=-product(shape(hamovr))*kind(hamovr)
   deallocate(hamovr,stat=i_stat)
   call memocc(i_stat,i_all,'hamovr',subname)

   i_all=-product(shape(work))*kind(work)
   deallocate(work,stat=i_stat)
   call memocc(i_stat,i_all,'work',subname)

   i_all=-product(shape(ew))*kind(ew)
   deallocate(ew,stat=i_stat)
   call memocc(i_stat,i_all,'ew',subname)


   !deallocate real array of wavefunctions
   if(exctX .or. in%SIC%approach=='NK')then
      i_all=-product(shape(psirocc))*kind(psirocc)
      deallocate(psirocc,stat=i_stat)
      call memocc(i_stat,i_all,'psirocc',subname)
   end if


   if(iter <=in%itermax+100) then
      if(iproc==0) call yaml_map('Iteration for Davidson convergence',iter-1)
      !if(iproc==0) write(*,'(1x,a,i3,a)') "Davidson's method: Convergence after ",iter-1,' iterations.'
   end if
   !finalize: Retranspose, deallocate

   ! Send all eigenvalues to all procs.
   call broadcast_kpt_objects(nproc,orbsv%nkpts,orbsv%norb,e(1,1,1),orbsv%ikptproc)

   call timing(iproc,'Davidson      ','OF')

   !retranspose v and psi
   call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,v(1),psiw(1))

   !resize work array before final transposition
   if(nproc > 1)then
      i_all=-product(shape(psiw))*kind(psiw)
      deallocate(psiw,stat=i_stat)
      call memocc(i_stat,i_all,'psiw',subname)

      allocate(psiw(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,psiw,'psiw',subname)
   end if

   call untranspose_v(iproc,nproc,orbs,Lzd%Glr%wfd,comms,psi(1),psiw(1))

   if(nproc > 1) then
      i_all=-product(shape(psiw))*kind(psiw)
      deallocate(psiw,stat=i_stat)
      call memocc(i_stat,i_all,'psiw',subname)
   end if

   i_all=-product(shape(hv))*kind(hv)
   deallocate(hv,stat=i_stat)
   call memocc(i_stat,i_all,'hv',subname)

   !copy the values in the eval array of the davidson procedure
   do ikpt=1,orbsv%nkpts
      do iorb=1,orbsv%norb
         !write(*,*) 'iorb,e(iorb,ikpt,1)',iorb,e(iorb,ikpt,1)
         !orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=min(e(iorb,ikpt,1),-.5d0)
         orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=e(iorb,ikpt,1)
      end do
   end do

   i_all=-product(shape(e))*kind(e)
   deallocate(e,stat=i_stat)
   call memocc(i_stat,i_all,'e',subname)

   !calculate gap
   call calculate_HOMO_LUMO_gap(iproc,orbs,orbsv)

   !write the results on the screen
   call write_eigen_objects(iproc,occorbs,nspin,nvirt,in%nplot,&
        Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),at,rxyz,Lzd%Glr,&
        orbs,orbsv,psi,v,in%output_wf_format)

   deallocate(confdatarr)

   if (GPUconv) then
      call free_gpu(GPU,orbsv%norbp)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,orbsv,in%nspin)
   end if

END SUBROUTINE davidson

!>   Generate upper triangular matrix in the subspace of Davidson algorithm
subroutine Davidson_subspace_hamovr(norb,nspinor,ncplx,nvctrp,hamovr,v,g,hv,hg)
   use module_base
   implicit none
   integer, intent(in) :: norb,nvctrp,nspinor,ncplx
   real(wp), dimension(nspinor*nvctrp*norb), intent(in) :: v,g,hv,hg
   real(wp), dimension(ncplx,2*norb,2*norb,2), intent(out) :: hamovr
   !local variables
   !n(c) character(len=*), parameter :: subname='Davidson_subspace_hamovr'
   integer :: iorb,jorb,icplx,ncomp

   if (nspinor == 4) then
      ncomp=2
   else
      ncomp=1
   end if

   !                 <vi | hvj>      <vi | hgj-n>                   <vi | vj>      <vi | gj-n>
   ! hamovr(i,j,1)=                               ;  hamovr(i,j,2)=  
   !                 <gi-n | hvj>  <gi-n | hgj-n>                   <gi-n | vj>  <gi-n | gj-n>

   !  do iorb=1,norb
   !     do jorb=iorb,norb
   !        hamovr(1,iorb,jorb,1)=               dot(nvctrp,v(1,iorb),1,hv(1,jorb),1)
   !        hamovr(1,jorb,iorb+norb,1)=        dot(nvctrp,g(1,iorb),1,hv(1,jorb),1)
   !        hamovr(1,iorb,jorb+norb,1)=        dot(nvctrp,v(1,iorb),1,hg(1,jorb),1)
   !        hamovr(1,iorb+norb,jorb+norb,1)= dot(nvctrp,g(1,iorb),1,hg(1,jorb),1)
   !               
   !        hamovr(1,iorb,jorb,2)=               dot(nvctrp,v(1,iorb),1, v(1,jorb),1)
   !        hamovr(1,jorb,iorb+norb,2)=       dot(nvctrp,g(1,iorb),1, v(1,jorb),1)
   !        hamovr(1,iorb,jorb+norb,2)=        dot(nvctrp,v(1,iorb),1, g(1,jorb),1)
   !        hamovr(1,iorb+norb,jorb+norb,2)= dot(nvctrp,g(1,iorb),1, g(1,jorb),1)
   !     enddo
   !  enddo
   !
   !use lapack operations to generalise the calculation to different nspinor

   !4 gemm + 2 dsyrk operations

   !<vi | hvj> 
   if(nspinor==1) then
      call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,v(1),&
         &   max(1,nvctrp),hv(1),max(1,nvctrp),0.0_wp,&
         &   hamovr(1,1,1,1),2*norb)
   else
      call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),v(1),&
         &   max(1,ncomp*nvctrp), &
         &   hv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
         &   hamovr(1,1,1,1),2*norb)
   end if

   !<gi | hvj> 
   if(nspinor==1) then
      call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
         &   max(1,nvctrp),hv(1),max(1,nvctrp),0.0_wp,&
         &   hamovr(1,norb+1,1,1),2*norb)
   else
      call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
         &   max(1,ncomp*nvctrp), &
         &   hv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
         &   hamovr(1,norb+1,1,1),2*norb)
   end if

   !<gi | hgj>
   if(nspinor==1) then
      call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
         &   max(1,nvctrp),hg(1),max(1,nvctrp),0.0_wp,&
         &   hamovr(1,norb+1,norb+1,1),2*norb)
   else
      call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
         &   max(1,ncomp*nvctrp), &
         &   hg(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
         &   hamovr(1,norb+1,norb+1,1),2*norb)
   end if

   !<vi | vj> 
   if(nspinor==1) then
      call syrk('U','T',norb,nvctrp,1.0_wp,v(1),max(1,nvctrp),&
         &   0.0_wp,hamovr(1,1,1,2),2*norb)
   else
      call herk('U','C',norb,ncomp*nvctrp,1.0_wp,v(1),max(1,ncomp*nvctrp),&
         &   0.0_wp,hamovr(1,1,1,2),2*norb)
   end if

   !<gi | vj> => hsub(:,:,:,5)
   if(nspinor==1) then
      call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
         &   max(1,nvctrp),v(1),max(1,nvctrp),0.0_wp,&
         &   hamovr(1,norb+1,1,2),2*norb)
   else
      call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
         &   max(1,ncomp*nvctrp), &
         &   v(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
         &   hamovr(1,norb+1,1,2),2*norb)
   end if

   !<gi | gj>
   if(nspinor==1) then
      call syrk('U','T',norb,nvctrp,1.0_wp,g(1),max(1,nvctrp),&
         &   0.0_wp,hamovr(1,norb+1,norb+1,2),2*norb)
   else
      call herk('U','C',norb,ncomp*nvctrp,1.0_wp,g(1),max(1,ncomp*nvctrp),&
         &   0.0_wp,hamovr(1,norb+1,norb+1,2),2*norb)
   end if


   !fill hamovr (Upper triangular)
   do jorb=1,norb
      do iorb=1,norb
         do icplx=1,ncplx
            hamovr(icplx,iorb,jorb+norb,1) = (-1)**(icplx-1)*hamovr(icplx,jorb+norb,iorb,1)  
            hamovr(icplx,iorb,jorb+norb,2) = (-1)**(icplx-1)*hamovr(icplx,jorb+norb,iorb,2)  
         end do
      enddo
   enddo

END SUBROUTINE Davidson_subspace_hamovr


subroutine update_psivirt(norb,nspinor,ncplx,nvctrp,hamovr,v,g,work)
   use module_base
   implicit none
   integer, intent(in) :: norb,nvctrp,nspinor,ncplx
   real(wp), dimension(nspinor*nvctrp*norb), intent(in) :: g
   real(wp), dimension(ncplx,2*norb,2*norb), intent(in) :: hamovr
   real(wp), dimension(nspinor*nvctrp*norb), intent(inout) :: v
   real(wp), dimension(nspinor*nvctrp*norb), intent(inout) :: work
   !local variables
   !n(c) character(len=*), parameter :: subname='update_psivirt'
   integer :: ncomp

   if (nspinor == 4) then
      ncomp=2
   else
      ncomp=1
   end if

   !Update v, that is the wavefunction, using eigenvectors stored in hamovr(:,:,1)
   !Lets say we have 4 quarters top/bottom left/right, then
   !v = matmul(v, hamovr(topleft)  ) + matmul(g, hamovr(bottomleft)  )  needed    
   !g=  matmul(v, hamovr(topright) ) + matmul(g, hamovr(bottomright) ) not needed
   !use hv as work arrray

   !Note: The previous data layout allowed level 3 BLAS
   !call DGEMM('N','N',nvctrp,nvirte,n2virt,1.d0,v(1,1),nvctrp,hamovr(1,1,1),n2virt,0.d0,hv(1,1),nvctrp)
   !    dimensions    =m      =n   =k          m,k        k,n                   m,n             
   !call vcopy(nvctrp*nvirte,hv(1,1),1,v(1,1),1)

   if(nspinor==1) then
      call gemm('N','N',nvctrp,norb,norb,1.0_wp,v(1),&
         &   max(1,nvctrp),hamovr(1,1,1),max(1,2*norb),0.0_wp,&
         &   work(1),nvctrp)
      call gemm('N','N',nvctrp,norb,norb,1.0_wp,g(1),&
         &   max(1,nvctrp),hamovr(1,norb+1,1),max(1,2*norb),1.0_wp,&
         &   work(1),nvctrp)

   else
      call c_gemm('N','N',ncomp*nvctrp,norb,norb,(1.0_wp,0.0_wp),v(1),&
         &   max(1,ncomp*nvctrp),hamovr(1,1,1),max(1,2*norb),(0.0_wp,0.0_wp),&
         &   work(1),ncomp*nvctrp)
      call c_gemm('N','N',ncomp*nvctrp,norb,norb,(1.0_wp,0.0_wp),g(1),&
         &   max(1,ncomp*nvctrp),hamovr(1,norb+1,1),max(1,2*norb),(1.0_wp,0.0_wp),&
         &   work(1),ncomp*nvctrp)
   end if

   call vcopy(nspinor*nvctrp*norb,work(1),1,v(1),1)

END SUBROUTINE update_psivirt


subroutine psivirt_from_gaussians(iproc,nproc,at,orbs,Lzd,comms,rxyz,hx,hy,hz,nspin,psivirt)
   use module_base
   use module_types
   use module_interfaces
   use gaussians
   use communications_base, only: comms_cubic
   use communications, only: transpose_v
   implicit none
   integer, intent(in) :: iproc,nproc,nspin
   real(gp), intent(in) :: hx,hy,hz
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd
   type(comms_cubic), intent(in) :: comms
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), dimension(orbs%npsidim_orbs), intent(out) :: psivirt
   !local variables
   character(len=*), parameter :: subname='psivirt_from_gaussians'
   logical ::  randinp
   integer :: iorb,icoeff,i_all,i_stat,nwork,info,jorb,ikpt,korb
   integer :: iseg,i0,i1,i2,i3,jj,ispinor,i,ind_c,ind_f,jcoeff
   real(wp) :: rfreq,gnrm_fake
   real(wp), dimension(:,:,:), allocatable :: gaucoeffs
   real(gp), dimension(:), allocatable :: work,ev
   real(gp), dimension(:,:), allocatable :: ovrlp
   type(gaussian_basis) :: G
   real(wp), dimension(:), pointer :: gbd_occ,psiw


   !initialise some coefficients in the gaussian basis
   !nullify the G%rxyz pointer
   nullify(G%rxyz)
   !extract the gaussian basis from the pseudowavefunctions
   !use a better basis than the input guess
   call gaussian_pswf_basis(11,.false.,iproc,nspin,at,rxyz,G,gbd_occ)

   allocate(gaucoeffs(G%ncoeff,orbs%nspinor,orbs%norbp+ndebug),stat=i_stat)
   call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

   !the kinetic overlap is correctly calculated only with Free BC
   randinp =.true.!.false.!lr%geocode /= 'F'

   if (randinp) then
      call to_zero(G%ncoeff*orbs%norbp*orbs%nspinor,gaucoeffs)
      if (G%ncoeff >= orbs%norb) then
         do icoeff=1,G%ncoeff
            !choose the orbital which correspond to this coefficient
            jorb=modulo(icoeff-1,orbs%norb)+1
            !fo any of the k-points associated to the processor fill
            do iorb=1,orbs%norbp
               !orbital at the net of k-point
               ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
               korb=orbs%isorb+iorb-(ikpt-1)*orbs%norb
               if (korb==jorb) then
                  gaucoeffs(icoeff,1,iorb)=cos(real(jorb+icoeff,wp))
                  if (orbs%nspinor == 4) then
                     gaucoeffs(icoeff,3,iorb)=sin(real(jorb+icoeff,wp))
                  end if
               end if
            end do
         end do
      else
         do iorb=1,orbs%norbp
            !orbital at the net of k-point
            ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
            korb=orbs%isorb+iorb-(ikpt-1)*orbs%norb
            !choose the coefficients which are associated to this orbital
            jcoeff=modulo(korb-1,G%ncoeff)+1
            do icoeff=1,G%ncoeff
               if (icoeff==jcoeff) then
                  gaucoeffs(icoeff,1,iorb)=1.0_gp!cos(real(korb+icoeff,wp))
                  if (orbs%nspinor == 4) then
                     gaucoeffs(icoeff,3,iorb)=sin(real(korb+icoeff,wp))
                  end if
               end if
            end do
         end do
         !write(*,*)'ERROR, not enough gaussian coefficients',G%ncoeff,orbs%norb
         !stop
      end if
      !!$     !fill randomly the gaussian coefficients for the orbitals considered
      !!$     do icoeff=1,G%ncoeff !reversed loop
      !!$        !be sure to call always a different random number, per orbital
      !!$        do jorb=1,orbs%isorb*orbs%nspinor
      !!$           tt=builtin_rand(idum) !call random_number(tt)
      !!$        end do
      !!$        do iorb=1,orbs%norbp*orbs%nspinor
      !!$           !do jproc=0,iproc-1
      !!$           !   tt=builtin_rand(idum) !call random_number(tt)
      !!$           !end do
      !!$           tt=builtin_rand(idum) !call random_number(tt)
      !!$           gaucoeffs(icoeff,iorb)=real(tt,wp)
      !!$           !do jproc=iproc+1,nproc-1
      !!$           !   tt=builtin_rand(idum) !call random_number(tt)
      !!$           !end do
      !!$        end do
      !!$        do iorb=(orbs%isorb+orbs%norbp)*orbs%nspinor+1,orbs%norb*orbs%nspinor
      !!$           tt=builtin_rand(idum) !call random_number(tt)
      !!$        end do
      !!$     end do

      !othogonalise the gaussian basis (wrong with k-points)
      !call gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)

   else
      !as an alternative strategy we may take the eigenvectors of the kinetic+k hamiltonian

      !in view of complete gaussian calculation
      allocate(ovrlp(G%ncoeff,G%ncoeff),stat=i_stat)
      call memocc(i_stat,ovrlp,'ovrlp',subname)

      !overlap calculation of the kinetic operator, upper triangular part
      !call kinetic_overlap(G,G,ovrlp)
      call gaussian_overlap(G,G,ovrlp)
      nwork=3*G%ncoeff+1
      allocate(work(nwork+ndebug),stat=i_stat)
      call memocc(i_stat,work,'work',subname)
      allocate(ev(G%ncoeff+ndebug),stat=i_stat)
      call memocc(i_stat,ev,'ev',subname)

      !!$  if (iproc == 0) then
      !!$     do iat=1,G%ncoeff
      !!$        write(*,'(a,i0,10(1pe15.7))')'T',iat,ovrlp(1:iat,iat)
      !!$     end do
      !!$  end if

      !print *'nwork',nwork,3*nbasis-1
      call dsyev('V','U',G%ncoeff,ovrlp(1,1),G%ncoeff,ev(1),work(1),nwork,info)
      if (info /= 0) then
         if (iproc == 0) then
            write(*,*)'DSyev Error',info
         end if
         stop
      end if

      !!$  if (iproc == 0) then
      !!$     do iat=1,G%ncoeff
      !!$        write(*,'(a,i0,10(1pe15.7))')'Ev',iat,ovrlp(:,iat)
      !!$     end do
      !!$     do iat=1,G%ncoeff
      !!$        write(*,'(a,i0,10(1pe15.7))')'K',iat,ev(iat)
      !!$     end do
      !!$  end if

      !copy the eigenvectors to the matrix
      call to_zero(G%ncoeff*orbs%norbp*orbs%nspinor,gaucoeffs)
      if (orbs%norb > G%ncoeff) stop 'wrong gaussian basis'
      jorb=mod(orbs%isorb,orbs%norb)
      do iorb=1,orbs%norbp
         jorb=jorb+1
         if (jorb == orbs%norb+1) jorb=1 !for k-points calculation
         call vcopy(G%ncoeff,ovrlp(1,jorb),1,gaucoeffs(1,1,iorb),orbs%nspinor)
      end do


      i_all=-product(shape(ovrlp))*kind(ovrlp)
      deallocate(ovrlp,stat=i_stat)
      call memocc(i_stat,i_all,'ovrlp',subname)
      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work',subname)
      i_all=-product(shape(ev))*kind(ev)
      deallocate(ev,stat=i_stat)
      call memocc(i_stat,i_all,'ev',subname)

      !call MPI_BARRIER(bigdft_mpi%mpi_comm,info)
      !stop

   end if

   call gaussians_to_wavelets_new(iproc,nproc,Lzd,orbs,G,gaucoeffs,psivirt)

   !deallocate the gaussian basis descriptors
   call deallocate_gwf(G,subname)

   !deallocate gaussian array
   i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
   deallocate(gaucoeffs,stat=i_stat)
   call memocc(i_stat,i_all,'gaucoeffs',subname)
   i_all=-product(shape(gbd_occ))*kind(gbd_occ)
   deallocate(gbd_occ,stat=i_stat)
   call memocc(i_stat,i_all,'gbd_occ',subname)

   !add random background to the wavefunctions
   if (randinp .and. G%ncoeff >= orbs%norb) then
      !call to_zero(orbs%npsidim,psivirt)
      do iorb=1,orbs%norbp
         jorb=iorb+orbs%isorb
         do ispinor=1,orbs%nspinor
            !pseudo-random frequency (from 0 to 10*2pi)
            rfreq=real(jorb,wp)/real(orbs%norb*orbs%nkpts,wp)*62.8318530717958648_wp
            do iseg=1,Lzd%Glr%wfd%nseg_c
               call segments_to_grid(Lzd%Glr%wfd%keyvloc(iseg),Lzd%Glr%wfd%keygloc(1,iseg),Lzd%Glr%d,i0,i1,i2,i3,jj)
               do i=i0,i1
                  ind_c=i-i0+jj+((iorb-1)*orbs%nspinor+ispinor-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)
                  psivirt(ind_c)=psivirt(ind_c)+0.5_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
               end do
            end do
            do iseg=Lzd%Glr%wfd%nseg_c+1,Lzd%Glr%wfd%nseg_c+Lzd%Glr%wfd%nseg_f
               call segments_to_grid(Lzd%Glr%wfd%keyvloc(iseg),Lzd%Glr%wfd%keygloc(1,iseg),Lzd%Glr%d,i0,i1,i2,i3,jj)
               do i=i0,i1
                  ind_f=Lzd%Glr%wfd%nvctr_c+7*(i-i0+jj-1)+&
                     &   ((iorb-1)*orbs%nspinor+ispinor-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)
                  psivirt(ind_f+1)=psivirt(ind_f+1)+0.4_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
                  psivirt(ind_f+2)=psivirt(ind_f+2)+0.35_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
                  psivirt(ind_f+3)=psivirt(ind_f+3)+0.3_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
                  psivirt(ind_f+4)=psivirt(ind_f+4)+0.25_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
                  psivirt(ind_f+5)=psivirt(ind_f+5)+0.2_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
                  psivirt(ind_f+6)=psivirt(ind_f+6)+0.15_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
                  psivirt(ind_f+7)=psivirt(ind_f+7)+0.1_wp*&
                     &   sin(rfreq*(i1+real(jorb,wp)))*sin(rfreq*(i2+real(jorb,wp)))*sin(rfreq*(i3+real(jorb,wp)))
               end do
            end do
         end do
      end do
      !after having added random background, precondition the wavefunctions with an ncong of 10
      call preconditionall(orbs,Lzd%Glr,hx,hy,hz,10,psivirt,gnrm_fake,gnrm_fake)
   end if

   !transpose v
   if(nproc > 1)then
      !reallocate the work array with the good size
      allocate(psiw(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,psiw,'psiw',subname)
   end if

   !transpose the wavefunction in wavelet basis
   call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psivirt(1),psiw(1))

   !here one has to decide whether leave things like that or
   !multiply the transposed wavefunctions by the matrix of the coefficients

   if(nproc > 1)then
      i_all=-product(shape(psiw))*kind(psiw)
      deallocate(psiw,stat=i_stat)
      call memocc(i_stat,i_all,'psiw',subname)
   end if


END SUBROUTINE psivirt_from_gaussians


!> Write eigenvalues and related quantities
!! @todo: must add the writing directory to the files
subroutine write_eigen_objects(iproc,occorbs,nspin,nvirt,nplot,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt,output_wf_format)
   use module_base
   use module_types
   use yaml_output
   implicit none
   logical, intent(in) :: occorbs
   integer, intent(in) :: iproc,nspin,nvirt,nplot,output_wf_format
   real(gp), intent(in) :: hx,hy,hz
   type(atoms_data), intent(in) :: at
   type(locreg_descriptors), intent(in) :: lr
   type(orbitals_data), intent(in) :: orbs,orbsv
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), dimension(:), pointer :: psi,psivirt
   !local variables
   character(len=10) :: comment
   character(len=11) :: orbname,denname
   integer :: iorb,ikpt,jorb,ind,occnorb,occnorbu,occnorbd
   real(gp) :: valu,vald,val

   if (occorbs) then
      occnorb = 0
      occnorbu = 0
      occnorbd = 0
   else
      occnorb = orbs%norb
      occnorbu = orbs%norbu
      occnorbd = orbs%norbd
   end if

   if(iproc==0)then
      call yaml_open_sequence('Complete list of energy eigenvalues')
      if (nspin==1) then
         !write(*,'(1x,a)')'Complete list of energy eigenvalues'
         do ikpt=1,orbsv%nkpts
            call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))) // ' BZ coord. = ' // &
            & trim(yaml_toa(orbs%kpts(:, ikpt),fmt='(f12.6)')))
            !if (orbsv%nkpts > 1) write(*,"(1x,A,I3.3,A,3F12.6)") "Kpt #", ikpt, " BZ coord. = ", orbsv%kpts(:, ikpt)
            do iorb=1,orbs%norb
               if (occorbs) then
                  val = orbs%eval(iorb+(ikpt-1)*orbs%norb)
               else
                  val = orbsv%eval(iorb+(ikpt-1)*orbsv%norb)!(iorb, ikpt, 1)
               end if
               call yaml_sequence(advance='no')
               call yaml_map('e_occupied',val, fmt='(1pe21.14)',advance='no')
               call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)'))) 
               !write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_occupied(',iorb,')=',val
            end do
            call yaml_sequence(advance='no')
            call yaml_map('HOMO LUMO gap (AU, eV)', &
               &   (/ orbsv%eval(1+occnorb+(ikpt-1)*orbsv%norb)-val,&
               &   Ha_eV*(orbsv%eval(1+occnorb+(ikpt-1)*orbsv%norb)-val) /), fmt='(1pe21.14)')
            !write(*,'(1x,a,1pe21.14,a,0pf8.4,a)')&
            !   &   'HOMO LUMO gap   =',orbsv%eval(1+occnorb+(ikpt-1)*orbsv%norb)-val,&
            !   &   ' (',Ha_eV*(orbsv%eval(1+occnorb+(ikpt-1)*orbsv%norb)-val),&
            !   &   ' eV)'
            do iorb=1,orbsv%norb - occnorb
               call yaml_sequence(advance='no')
               call yaml_map('e_virtual',&
                    & orbsv%eval(iorb+occnorb+(ikpt-1)*orbsv%norb), fmt='(1pe21.14)',advance='no')
               call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
               !write(*,'(1x,a,i4,a,1x,1pe21.14)') &
               !   &   'e_virtual(',iorb,')=',orbsv%eval(iorb+occnorb+(ikpt-1)*orbsv%norb)!e(iorb+occnorb,ikpt,1)
            end do
         end do
         call yaml_close_sequence()
      else
         do ikpt=1,orbsv%nkpts
            call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))) // ' BZ coord. = ' // &
            & trim(yaml_toa(orbs%kpts(:, ikpt),fmt='(f12.6)')))
            !write(*,'(1x,a)')'Complete list of energy eigenvalues'
            do iorb=1,min(orbs%norbu,orbs%norbd)
               if (occorbs) then
                  valu = orbs%eval(iorb+(ikpt-1)*orbs%norb)
                  vald = orbs%eval(iorb+orbs%norbu+(ikpt-1)*orbs%norb)
               else
                  valu = orbsv%eval(iorb+(ikpt-1)*orbsv%norb)!e(iorb, ikpt, 1)
                  vald = orbsv%eval(iorb+orbsv%norbu+(ikpt-1)*orbsv%norb)!e(iorb+orbsv%norbu, ikpt, 1)
               end if
               call yaml_sequence(advance='no')
               call yaml_map('e_occ',(/ valu,vald /), fmt='(1pe21.14)',advance='no')
               call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
               !write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
               !   &   'e_occ(',iorb,',u)=',valu,'e_occ(',iorb,',d)=',vald
            end do
            if (orbs%norbu > orbs%norbd) then
               do iorb=orbs%norbd+1,orbs%norbu
                  if (occorbs) then
                     valu = orbs%eval(iorb+(ikpt-1)*orbs%norb)
                  else
                     valu = orbsv%eval(iorb+(ikpt-1)*orbsv%norb)!e(iorb, ikpt, 1)
                  end if
                  call yaml_sequence(advance='no')
                  call yaml_map('e_occ',valu, fmt='(1pe21.14)',advance='no')
                  call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
                  !write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_occ(',iorb,',u)=',valu
               end do
            else if (orbs%norbd > orbs%norbu) then
               do iorb=orbs%norbu+1,orbs%norbd
                  if (occorbs) then
                     vald = orbs%eval(iorb+orbs%norbu+(ikpt-1)*orbs%norb)
                  else
                     vald = orbsv%eval(iorb+orbsv%norbu+(ikpt-1)*orbsv%norb)!e(iorb+orbsv%norbu, ikpt, 1)
                  end if
                  call yaml_sequence(advance='no')
                  call yaml_map('e_occ',vald, fmt='(1pe21.14)',advance='no')
                  call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
                  !write(*,'(50x,a,i4,a,1x,1pe21.14)') 'e_occ(',iorb,',d)=',vald
               end do
            end if
            call yaml_sequence(advance='no')
            call yaml_map('HOMO LUMO gap (AU, eV)', &
               &   (/ orbsv%eval(1+occnorbu+(ikpt-1)*orbsv%norb)-valu, &
               &      orbsv%eval(orbsv%norbu+1+occnorbd+(ikpt-1)*orbsv%norb)-vald, &
               &      Ha_eV*(orbsv%eval(1+occnorbu+(ikpt-1)*orbsv%norb)-valu), &
               &      Ha_eV*(orbsv%eval(orbsv%norbu+1+occnorbd+(ikpt-1)*orbsv%norb)-vald) /), fmt='(1pe21.14)')
            !write(*,'(1x,a,1x,1pe21.14,a,0pf8.4,a,a,1x,1pe21.14,a,0pf8.4,a)') &
            !   &   'HOMO LUMO gap, u =', orbsv%eval(1+occnorbu+(ikpt-1)*orbsv%norb)-valu,&
            !   &   ' (',Ha_eV*(orbsv%eval(1+occnorbu+(ikpt-1)*orbsv%norb)-valu),' eV)',&
            !   &   ',d =',orbsv%eval(orbsv%norbu+1+occnorbd+(ikpt-1)*orbsv%norb)-vald,&
            !   &   ' (',Ha_eV*(orbsv%eval(orbsv%norbu+1+occnorbd+(ikpt-1)*orbsv%norb)-vald),' eV)'
            do iorb=1,min(orbsv%norbu-occnorbu,orbsv%norbd-occnorbd)
               jorb=orbsv%norbu+iorb+occnorbd
               call yaml_sequence(advance='no')
               call yaml_map('e_virt', (/ &
               &  orbsv%eval(iorb+(ikpt-1)*orbsv%norb), &
               &  orbsv%eval(jorb+(ikpt-1)*orbsv%norb) /), fmt='(1pe21.14)', advance='no')
               call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
               !write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
               !   &   'e_vrt(',iorb,',u)=',orbsv%eval(iorb+(ikpt-1)*orbsv%norb),&!e(iorb,ikpt,1),&
               !   &   'e_vrt(',iorb,',d)=',orbsv%eval(jorb+(ikpt-1)*orbsv%norb)!e(jorb,ikpt,1)
            end do
            call yaml_close_map()
            if (orbsv%norbu-occnorbu > orbsv%norbd-occnorbd) then
               do iorb=orbsv%norbd+1-occnorbu,orbsv%norbu-occnorbd
                  call yaml_sequence(advance='no')
                  call yaml_map('e_vrt u', &
                  &  orbsv%eval(iorb+(ikpt-1)*orbsv%norb), fmt='(1pe21.14)',advance='no')
                  call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
                  !write(*,'(1x,a,i4,a,1x,1pe21.14)') &
                  !   &   'e_vrt(',iorb,',u)=',orbsv%eval(iorb+(ikpt-1)*orbsv%norb)!e(iorb,ikpt,1)
               end do
            else if (orbsv%norbd-occnorbd > orbsv%norbu-occnorbu) then
               do iorb=2*orbsv%norbu+1-occnorbu,orbsv%norbu-occnorbu+orbsv%norbd-occnorbd
                  call yaml_sequence(advance='no')
                  call yaml_map('e_vrt d', &
                  &  orbsv%eval(iorb+(ikpt-1)*orbsv%norb), fmt='(1pe21.14)',advance='no')
                  call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)')))
                  !write(*,'(50x,a,i4,a,1x,1pe21.14)') &
                  !   &   'e_vrt(',iorb-orbsv%norbu-occnorbu,',d)=',orbsv%eval(iorb+(ikpt-1)*orbsv%norb)!e(iorb,ikpt,1)
               end do
            end if
         end do
         call yaml_close_sequence()
      end if
   end if

END SUBROUTINE write_eigen_objects


! PLOTTING

!>plot the converged wavefunctions in the different orbitals.
!nplot is the requested total of orbitals to plot, where
!states near the HOMO/LUMO gap are given higher priority.
!Occupied orbitals are only plotted when nplot>nvirt,
!otherwise a comment is given in the out file.
subroutine dump_eigenfunctions(dir_output,nplot,at,hgrids,lr,orbs,orbsv,rxyz,psi,psivirt)
  use module_base, only: gp,wp
  use locregs, only: locreg_descriptors
  use module_types, only: atoms_data,orbitals_data
  implicit none
  integer, intent(in) :: nplot !<number of eigenfuncitions to be plotted close to the fermi level
  type(atoms_data), intent(in) :: at !<descriptor of atomic properties
  type(orbitals_data), intent(in) :: orbs,orbsv !<orbitals, occupied and virtual respectively
  type(locreg_descriptors), intent(in) :: lr !<localization regions of the wavefunctions
  character(len=*), intent(in) :: dir_output !<directory where the data have to be put in
  real(gp), dimension(3), intent(in) :: hgrids !<grid spacings of the simulation domain
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz !<atomic positions
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi,psivirt !<occupied and virtual eigenfunctions
  !local variables
  integer :: ind,iorb
  real(gp) :: hx,hy,hz
  character(len=300) :: orbname,denname
  
  hx=hgrids(1)
  hy=hgrids(2)
  hz=hgrids(3)

  !add a modulo operator to get rid of the particular k-point
  do iorb=1,orbsv%norbp!requested: nvirt of nvirte orbitals

     if(modulo(iorb+orbsv%isorb-1,orbsv%norb)+1 > abs(nplot)) then
        exit 
        !if(iproc == 0 .and. abs(nplot) > 0) write(*,'(A)')'No plots of occupied orbitals requested.'
     end if

     ind=1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*(iorb-1)
     !plot the orbital and the density
     write(orbname,'(A,i4.4)')trim(dir_output)//'virtual',iorb+orbsv%isorb
     write(denname,'(A,i4.4)')trim(dir_output)//'denvirt',iorb+orbsv%isorb
     !write(comment,'(1pe10.3)')orbsv%eval(iorb+orbsv%isorb)!e(modulo(iorb+orbsv%isorb-1,orbsv%norb)+1,orbsv%iokpt(iorb),1)

     call plot_wf(trim(orbname),1,at,1.0_wp,lr,hx,hy,hz,rxyz,psivirt(ind:))
     call plot_wf(trim(denname),2,at,1.0_wp,lr,hx,hy,hz,rxyz,psivirt(ind:))

  end do

  do iorb=orbs%norbp,1,-1 ! sweep over highest occupied orbitals
     if(modulo(orbs%norb-iorb-orbs%isorb-0,orbs%norb)+1 <=  abs(nplot)) then  ! SG 
        !address
        ind=1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*(iorb-1)
        write(orbname,'(A,i4.4)')trim(dir_output)//'orbital',iorb+orbs%isorb
        write(denname,'(A,i4.4)')trim(dir_output)//'densocc',iorb+orbs%isorb
        !write(comment,'(1pe10.3)')orbs%eval(iorb+orbs%isorb)

        call plot_wf(trim(orbname),1,at,1.0_wp,lr,hx,hy,hz,rxyz,psi(ind:))
        call plot_wf(trim(denname),2,at,1.0_wp,lr,hx,hy,hz,rxyz,psi(ind:))

     endif
  end do
  ! END OF PLOTTING
end subroutine dump_eigenfunctions


!> Calculate the gap and fill the value in the orbs structure
subroutine calculate_HOMO_LUMO_gap(iproc,orbs,orbsv)
   use module_base
   use module_types
   use yaml_output
   implicit none
   integer, intent(in) :: iproc
   type(orbitals_data), intent(in) :: orbsv
   type(orbitals_data), intent(inout) :: orbs
   !local variables
   integer :: ikpt

   if (orbs%nkpts /= orbsv%nkpts) then
      return
      stop 'HL gap with Band structure not implemented yet'
   end if
   !depending on nspin
   orbs%HLgap=UNINITIALIZED(orbs%HLgap)
   if (orbs%nspin==1) then
      ikpt=1
       orbs%HLgap=orbsv%eval(1+(ikpt-1)*orbsv%norb)-orbs%eval(orbs%norb+(ikpt-1)*orbs%norb)
      !the minimum wrt all the k-points
      do ikpt=2,orbs%nkpts
         orbs%HLgap=min(orbs%HLgap,orbsv%eval(1+(ikpt-1)*orbsv%norb)&
              -orbs%eval(orbs%norb+(ikpt-1)*orbs%norb))
      end do
   else if (orbs%nspin==2) then
      ikpt=1
      orbs%HLgap=min(orbsv%eval(1+(ikpt-1)*orbsv%norb)-orbs%eval(orbs%norbu+(ikpt-1)*orbs%norb),&
           orbsv%eval(orbsv%norbu+1+(ikpt-1)*orbsv%norb)-orbs%eval(orbs%norbd+orbs%norbu+(ikpt-1)*orbs%norb))
      do ikpt=2,orbs%nkpts
         orbs%HLgap=min(orbs%HLgap,orbsv%eval(1+(ikpt-1)*orbsv%norb)-orbs%eval(orbs%norbu+(ikpt-1)*orbs%norb),&
              orbsv%eval(orbsv%norbu+1+(ikpt-1)*orbsv%norb)-orbs%eval(orbs%norbd+orbs%norbu+(ikpt-1)*orbs%norb))
      end do
   end if

   !warning if gap is negative
   if (orbs%HLgap < 0.0_gp .and. orbs%HLgap/=uninitialized(orbs%HLgap)) then
      if (iproc==0) call yaml_warning('HLgap is negative, convergence problem?')
      !if (iproc==0) write(*,*)'WARNING!! HLgap is negative, convergence problem?' 
   end if

END SUBROUTINE calculate_HOMO_LUMO_gap

!> Add a potential to the local potential which has the function of confining the 
!! Solutions to a given value
subroutine add_confining_potential(n1i,n2i,n3i,nspin,eps,dencutoff,rpow,pot,rho)
   use module_base
   implicit none
   integer, intent(in) :: n1i,n2i,n3i,nspin
   real(gp) , intent(in) :: rpow,eps,dencutoff
   real(dp), dimension(n1i,n2i,n3i,nspin), intent(in) :: rho
   real(wp), dimension(n1i,n2i,n3i,nspin), intent(inout) :: pot
   !local variables
   integer :: i1,i2,i3,ispin
   real(dp) :: density

   do ispin=1,nspin
      do i3=1,n3i
         do i2=1,n2i
            do i1=1,n1i
               !charge density value (not optimized)
               if (nspin==2) then
                  density=rho(i1,i2,i3,1)+rho(i1,i2,i3,2)
               else
                  density=rho(i1,i2,i3,1)
               end if
               pot(i1,i2,i3,ispin)=pot(i1,i2,i3,ispin)+eps*((density+dencutoff)**rpow)
            end do
         end do
      end do
   end do
END SUBROUTINE add_confining_potential

subroutine add_parabolic_potential(geocode,nat,n1i,n2i,n3i,hxh,hyh,hzh,rlimit,rxyz,pot)
   use module_base
   implicit none
   character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
   integer, intent(in) :: n1i,n2i,n3i,nat
   real(gp), intent(in) :: rlimit,hxh,hyh,hzh
   real(gp), dimension(3,nat), intent(in) :: rxyz
   real(wp), dimension(n1i*n2i*n3i), intent(inout) :: pot
   !local variables
   logical :: perx,pery,perz
   integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,i,i1,i2,i3,isx,iex,isy,iey,isz,iez
   integer :: iat,ind
   real(gp) :: x,y,z,r2,rx,ry,rz
   real(gp), dimension(3) :: cxyz

   !conditions for periodicity in the three directions
   perx=(geocode /= 'F')
   pery=(geocode == 'P')
   perz=(geocode /= 'F')

   call ext_buffers(perx,nbl1,nbr1)
   call ext_buffers(pery,nbl2,nbr2)
   call ext_buffers(perz,nbl3,nbr3)

   !calculate the center of the molecule
   call to_zero(3,cxyz)
   do iat=1,nat
      do i=1,3
         cxyz(i)=cxyz(i)+rxyz(i,iat)
      end do
   end do
   do i=1,3
      cxyz(i)=cxyz(i)/real(nat,gp)
   end do

   rx=cxyz(1) 
   ry=cxyz(2)
   rz=cxyz(3)

   isx=-nbl1
   isy=-nbl2
   isz=-nbl3

   iex=n1i-nbl1-1
   iey=n2i-nbl2-1
   iez=n3i-nbl3-1

   do i3=isz,iez
      z=real(i3,gp)*hzh-rz
      do i2=isy,iey
         y=real(i2,gp)*hyh-ry
         do i1=isx,iex
            x=real(i1,gp)*hxh-rx
            r2=x**2+y**2+z**2
            !add the parabolic correction to the potential
            if (r2 > rlimit**2) then
               ind=i1+1+nbl1+(i2+nbl2)*n1i+(i3+nbl3)*n1i*n2i
               pot(ind)=pot(ind)+0.1_gp*(sqrt(r2)-rlimit)**2
            endif
         enddo
      enddo
   enddo

END SUBROUTINE add_parabolic_potential
