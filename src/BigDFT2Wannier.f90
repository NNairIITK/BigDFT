program BigDFT2Wannier

   use BigDFT_API
   use Poisson_Solver
   implicit none
   character :: filetype*4
   !etsf
   type(locreg_descriptors) :: Glr
   type(orbitals_data) :: orbs,orbs2,orbsv
   type(atoms_data) :: atoms
   type(input_variables) :: input 
   type(workarr_sumrho) :: w
   type(communications_arrays), target :: comms, comms2,commsv
   integer :: iproc, nproc, i_stat, nelec, ind, ierr, nvp, nbp, npsidim, npsidim2, npsidimv
   integer :: jproc, norb_tot,isvirt,isband,nvctrp,npp,nvirtu,nvirtd,nvirt
   integer, dimension(:), allocatable :: n_bands_par,n_virt_tot_par,isvirt_par,isband_par,isproj_par,n_proj_par,n_vctr_par
   real(gp), dimension(:,:), pointer :: rxyz, rxyz_old
   real(gp), dimension(:,:), allocatable :: radii_cf
   real(gp), dimension(3) :: shift
   real(wp), allocatable :: psi_etsf(:,:),psi_etsfv(:,:),sph_har_etsf(:),sph_daub(:),psir(:),psir_re(:),psir_im(:)
   real(wp), allocatable :: psi_daub_im(:),psi_daub_re(:),psi_etsf2(:),psi_etsf3(:)
   real(wp), allocatable :: mmnk_v_re(:), mmnk_v_im(:)
   real(wp), pointer :: pwork(:)
   character :: filename*20
   !cube
   integer :: nx, ny, nz, nb, nb1, nb2, nk, inn
   integer, allocatable, dimension(:) :: Z
   real(kind=8) :: bx(3), by(3), bz(3), b1, b2, b3, r0x, r0y, r0z 
   real(kind=8) :: xx, yy, zz
   real(kind=8), allocatable :: ylm(:,:,:), func_r(:,:,:), sph_har(:,:,:,:), psi(:,:,:,:)
   real(kind=8), allocatable :: virt (:,:,:), orb(:,:,:)
   real(kind=8), allocatable :: at_pos(:,:)
   real(kind=8), allocatable :: amnk(:,:), amnk_tot(:), amnk_guess(:), amnk_guess_sorted(:)
   real(kind=8), allocatable :: mmnk_re(:,:,:), mmnk_im(:,:,:), mmnk_tot(:,:)
   integer :: i, j, k, np
   character :: subname*10, seedname*16
   logical :: calc_only_A
   real, dimension(3,3) :: real_latt, recip_latt
   integer :: n_kpts, n_proj, n_nnkpts, n_excb, n_at, n_bands, s
   integer :: n_occ, n_virt, n_virt_tot
   logical :: w_sph, w_ang, w_rad, pre_check
   real, allocatable, dimension (:,:) :: kpts
   real, allocatable, dimension (:,:) :: ctr_proj, x_proj, z_proj
   integer, allocatable, dimension (:) :: l, mr, rvalue
   real, allocatable, dimension (:) :: zona
   integer, allocatable, dimension (:,:) :: k_plus_b
   integer, allocatable, dimension (:,:) :: G_vec
   integer, allocatable, dimension (:) :: excb
   integer, allocatable, dimension (:) :: virt_list, amnk_bands_sorted
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   ! possible junk ?
   integer :: ncplx, j_end,pshft
   integer, allocatable, dimension (:) :: orb_list


   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

   if (nproc > 1) then
      call timing(iproc,'parallel     ','IN')
   else
      call timing(iproc,'             ','IN')
   end if

! Read input.inter file
! It defines the name of the system studied and some important integers :
! the number of occupied and unoccupied bands to compute mmn and mmn matrix, and the total number of unoccupied states to use in the precheck mode.
! pre_check defines if the pre-check calculations have to be done.
! Other logicals define if spherical harmonics and/or their angular and radial parts have to be written.
   call timing(iproc,'Precondition  ','ON')
   call read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_sph, w_ang, w_rad)


if (filetype == 'etsf' .or. filetype == 'ETSF') then

   ! Initalise the variables for the calculation
   call standard_inputfile_names(input)
   call read_input_variables(iproc,'posinp',input, atoms, rxyz)

   if (iproc == 0) call print_general_parameters(nproc,input,atoms)

   allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)
   call system_properties(iproc,nproc,input,atoms,orbs,radii_cf,nelec)
   call orbitals_descriptors(iproc,nproc,orbs%norb,orbs%norbu,orbs%norbd,orbs%nspin,orbs%nspinor,&
        orbs%nkpts,orbs%kpts,orbs%kwgts,orbs2)

   ! Set-up number of states and shifting values.
   nvirtu = abs(input%norbv)
   nvirtd = 0
   if (input%nspin==2) nvirtd=nvirtu
   call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
       & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv)


   if(orbs%nkpts > 1) stop 'BigDFT2Wannier does not work for nkpts > 1'
 

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space
   call system_size(iproc,atoms,rxyz,radii_cf,input%crmult,input%frmult,input%hx,input%hy,input%hz,&
        Glr,shift)


   ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
   call createWavefunctionsDescriptors(iproc,input%hx,input%hy,input%hz,&
        atoms,rxyz,radii_cf,input%crmult,input%frmult,Glr)


   ! Allocate communications arrays (allocate it before Projectors because of the definition of iskpts and nkptsp)
   call orbitals_communicators(iproc,nproc,Glr,orbs,comms)
   if(orbs%nspinor > 1) STOP 'BigDFT2Wannier does not work for nspinor > 1'
   call orbitals_communicators(iproc,nproc,Glr,orbsv,commsv)  
   if(orbsv%nspinor > 1) STOP 'BigDFT2Wannier does not work for nspinor > 1'


   ! Read wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
   allocate(psi_etsf(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp*orbs%nspinor),stat=i_stat)
   call memocc(i_stat,psi_etsf,'psi_etsf',subname)
   allocate(psi_etsfv(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbsv%norbp*orbsv%nspinor),stat=i_stat)
   call memocc(i_stat,psi_etsfv,'psi_etsfv',subname)
   allocate(rxyz_old(3,atoms%nat),stat=i_stat)
   call memocc(i_stat,rxyz_old,'rxyz_old',subname)

   if(associated(orbs%eval)) nullify(orbs%eval)
   allocate(orbs%eval(orbs%norb*orbs%nkpts), stat=i_stat)
   filename='wavefunction.etsf'
   call read_waves_etsf(iproc,filename,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
      Glr%wfd,psi_etsf)

   if(associated(orbsv%eval)) nullify(orbsv%eval)
   allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
   filename='virtuals.etsf'
   call read_waves_etsf(iproc,filename,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
      Glr%wfd,psi_etsfv)

   
   ! Read integers in order to allocate tables used to store wavefunctions, l, mr, rvalue, zona, ...
   call read_nnkp_int_alloc(iproc,seedname, n_kpts, n_proj, n_nnkpts, n_excb)
   allocate(kpts(n_kpts,3))
   allocate(ctr_proj(n_proj,3))
   allocate(x_proj(n_proj,3))
   allocate(z_proj(n_proj,3))
   allocate(l(n_proj))
   allocate(mr(n_proj))
   allocate(rvalue(n_proj))
   allocate(zona(n_proj))
   allocate(k_plus_b(n_kpts*n_nnkpts,2))
   allocate(G_vec(n_kpts*n_nnkpts,3))
   allocate(excb(n_excb))


   ! Read Wannier90 .nnkp file.
   ! The most important informations to be read are : 
   !  - ctr_proj : the coordinates of the center of the projections
   !  - l, mr, rvalue and zona = Z/a : the parameters used to build spherical harmonics
   !  - k_plus_b and G_vec : the parameters used to build the nearest neighbours k-points
   call read_nnkp(iproc,seedname, calc_only_A, real_latt, recip_latt, n_kpts, n_proj, n_nnkpts, &
                n_excb, kpts, ctr_proj, x_proj, z_proj, l, mr, rvalue, &
                zona, k_plus_b, G_vec, excb)

   ! Initialise the arrays n_proj_par and isproj_par
   allocate(n_proj_par(0:nproc-1))
   allocate(isproj_par(0:nproc-1))
   call parallel_repartition_with_kpoints(nproc,1,n_proj,n_proj_par)
   !check the distribution
   norb_tot=0
   do jproc=0,nproc-1
      isproj_par(jproc)=norb_tot
      norb_tot=norb_tot+n_proj_par(jproc)
   end do

   ! Initialise the array n_vctr_par
   allocate(n_vctr_par(0:nproc-1))
   call parallel_repartition_with_kpoints(nproc,n_kpts,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,n_vctr_par)

   ! Simplification of the notations
   nx=Glr%d%n1i
   ny=Glr%d%n2i
   nz=Glr%d%n3i
   n_at=atoms%nat
   call initialize_work_arrays_sumrho(Glr,w)

   ! Allocations for Amnk calculation
   ! Define orbs2 and comms2, used for the spherical harmonics transposition.
   orbs2%isorb=isproj_par(iproc)
   orbs2%nkpts=1
   orbs2%nspinor=1
   orbs2%iskpts=0
   orbs2%norb=n_proj
   orbs2%norb_par(:)=n_proj_par(:)
   deallocate(orbs2%iokpt)
   allocate(orbs2%iokpt(n_proj_par(iproc)))
   orbs2%iokpt=1
   orbs2%norbp=n_proj_par(iproc)
   call orbitals_communicators(iproc,nproc,Glr,orbs2,comms2)
   npsidim2=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*n_proj_par(iproc),sum(comms2%ncntt(0:nproc-1)))
   allocate(ylm(nx,ny,nz))
   allocate(func_r(nx,ny,nz))
   allocate(sph_har_etsf(nx*ny*nz))
   allocate(amnk_bands_sorted(n_virt))
   call timing(iproc,'Precondition  ','OF')
   
   if (pre_check .eqv. .true.) then

      call timing(iproc,'CrtProjectors ','ON')
      ! Initialise the arrays n_virt_par, isvirt_par
      allocate(n_virt_tot_par(0:nproc-1))
      allocate(isvirt_par(0:nproc-1))
      call parallel_repartition_with_kpoints(nproc,n_kpts,n_virt_tot,n_virt_tot_par)
      !check the distribution
      norb_tot=0
      do jproc=0,nproc-1
         isvirt_par(jproc)=norb_tot
         norb_tot=norb_tot+n_virt_tot_par(jproc)
      end do
      !check the distribution
      norb_tot=0
      do jproc=0,iproc-1
         norb_tot=norb_tot+n_virt_tot_par(jproc)
      end do
      !reference orbital for process
      isvirt=norb_tot
      nvp=n_virt_tot_par(iproc)


      ! Initialisation :
      ! - define a new psi_etsf that can be transposed in order to keep psi_etsf in memory (psi_etsf will be used during the Mmnk matrix calculation,)
      npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsv%norbp*orbsv%nspinor,sum(commsv%ncntt(0:nproc-1)))
      allocate(psi_etsf2(npsidim),stat=i_stat)
      call memocc(i_stat,psi_etsf2,'psi_etsf2',subname)
      k=0
      do i=1,orbsv%norbp
         do j=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
            k=k+1
            psi_etsf2(k)=psi_etsfv(j,i)
         end do
      end do

      ! - b1, b2 and b3 are the norm of the lattice parameters.
      b1=atoms%alat1
      b2=atoms%alat2
      b3=atoms%alat3
      ! - Allocations
      allocate(amnk(n_virt_tot,n_proj))
      amnk=0.0
      allocate(amnk_guess(n_virt_tot))
      allocate(sph_daub(npsidim2), stat=i_stat)
      call razero(npsidim2, sph_daub)
      ! Begining of the algorithm to compute the scalar product in order to find the best unoccupied orbitals to use to compute the actual Amnk matrix :
      if (iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '! Calculating amnk=<virt|sph_har>  !'
         write(*,*) '!       in pre-check mode :        !'
         write(*,*) '!==================================!'
         write(*,'(A12,4x,A15)') 'Virtual band', 'amnk_guess(nb)='
      end if

      ! Calculation of the spherical harmonics in parallel.
      ! It is done in the real space and then converted in the Daubechies representation.
      pshft = 0
      do npp=1, orbs2%norbp
         np = npp + orbs2%isorb
         ! Convolution buffer : n1i=2*n1+31 -> explains the '15*input%hx*0.5' term
         r0x=ctr_proj(np,1)*b1+15*input%hx*0.5
         r0y=ctr_proj(np,2)*b2+15*input%hy*0.5
         r0z=ctr_proj(np,3)*b3+15*input%hz*0.5
         do k=1,nz
            zz=(k-1)*input%hz*0.5-r0z
            do j=1,ny
               yy=(j-1)*input%hy*0.5-r0y
               do i=1,nx
                  ind=(k-1)*ny*nx+(j-1)*nx+i
                  xx=(i-1)*input%hx*0.5-r0x
                  call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                        xx, yy, zz, n_proj, ylm)
                  call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                        xx, yy, zz, n_proj, func_r)
                  ! The 'sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)' term is here to normalize spherical harmonics
                  sph_har_etsf(ind)=func_r(i,j,k)*ylm(i,j,k)*sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)
               end do
            end do
         end do
         call isf_to_daub(Glr,w,sph_har_etsf(1),sph_daub(1+pshft))
         pshft=pshft + max(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,comms2%ncntt(iproc)/n_proj_par(iproc))
      end do

      ! Tranposition of the distribution of the BigDFT wavefunctions : orbitals -> components.
      call timing(iproc,'CrtProjectors ','OF')
      allocate(pwork(npsidim))
      call transpose_v(iproc,nproc,orbsv,Glr%wfd,commsv,psi_etsf2,work=pwork)
      deallocate(pwork)

      ! Tranposition of the distribution of the spherical harmonics: orbitals -> components.
      allocate(pwork(npsidim2))
      call transpose_v(iproc,nproc,orbs2,Glr%wfd,comms2,sph_daub,work=pwork)
      deallocate(pwork)
      call timing(iproc,'ApplyProj     ','ON')

      ! Scalar product of amnk=<sph_daub|psi> in parallel.
      call razero(n_proj*n_virt_tot,amnk)
      nvctrp=n_vctr_par(iproc)
      call gemm('T','N',n_virt_tot,n_proj,nvctrp,1.0_wp,psi_etsf2(1),max(1,nvctrp),&
           sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),n_virt_tot)

      ! Construction of the whole Amnk_guess matrix.
      call mpiallred(amnk(1,1),n_proj*n_virt_tot,MPI_SUM,MPI_COMM_WORLD,ierr)

      ! For each unoccupied orbitals, check how they project on spherical harmonics.
      ! The greater amnk_guess(nb) is, the more they project on spherical harmonics.
      do nb=1,n_virt_tot
         amnk_guess(nb)=0.0
         do np=1,n_proj
            amnk_guess(nb)=amnk_guess(nb)+(amnk(nb,np))**2
         end do
         if (iproc==0) write(*,'(I4,11x,F12.6)') nb, amnk_guess(nb)
      end do

      ! Choice of the unoccupied orbitals to calculate the Amnk matrix
      if (iproc==0) then
         write(*,*) 
         write(*,'(1a)') 'These are the virtual bands to use to construct the actual Amn and Mmn matrices :'
         write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_guess)='
      end if
      allocate(amnk_guess_sorted(n_virt))
      do nb=1,n_virt
         amnk_guess_sorted(nb)=maxval(amnk_guess,n_virt_tot)
         amnk_bands_sorted(nb)=maxloc(amnk_guess,n_virt_tot)
         amnk_guess(amnk_bands_sorted(nb))=0.d0
      if (iproc==0) write(*,'(I4,3x,F12.6)') amnk_bands_sorted(nb), amnk_guess_sorted(nb)
      end do

      ! End of the pre-check mode
      allocate(pwork(npsidim2))
      call timing(iproc,'ApplyProj     ','OF')
      call untranspose_v(iproc,nproc,orbs2,Glr%wfd,comms2,sph_daub,work=pwork)
      call timing(iproc,'ApplyProj     ','ON')
      deallocate(pwork)
      deallocate(sph_daub)
      deallocate(psi_etsf2)
      deallocate(amnk)
      deallocate(amnk_guess_sorted)
      deallocate(amnk_guess)
      if (iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '! Calculating amnk=<virt|sph_har>  !'
         write(*,*) '!     in pre-check mode done       !'
         write(*,*) '!==================================!'
         write(*,*)
         write(*,*)
      end if

!     ! Rewrite the input.inter file in order to leave the pre-check mode the next time the program is run
!     if (iproc==0) call write_inter(n_virt, amnk_bands_sorted)
      call timing(iproc,'ApplyProj     ','OF')
   end if

   call timing(iproc,'CrtDescriptors','ON')
   ! Define which unoccupied orbitals have to be read.
   ! Here, virt_list nominates the virtual orbitals chosen in pre-check mode.
   if (n_virt .ne. 0) then
      allocate(virt_list(n_virt))
      if (pre_check .eqv. .true.) then 
         do i=1,n_virt
            virt_list(i)=amnk_bands_sorted(i)!+n_occ
         end do
      else if (pre_check .eqv. .false.) then
         call read_inter_list(iproc, n_virt, virt_list)
      end if
   end if

   ! Initialise the arrays n_bands_par, isband_par
   n_bands=n_occ+n_virt
   allocate(n_bands_par(0:nproc-1))
   allocate(isband_par(0:nproc-1))
   call parallel_repartition_with_kpoints(nproc,n_kpts,n_bands,n_bands_par)
   !check the distribution
   norb_tot=0
   do jproc=0,nproc-1
      isband_par(jproc)=norb_tot
      norb_tot=norb_tot+n_bands_par(jproc)
   end do
   !check the distribution
   norb_tot=0
   do jproc=0,iproc-1
      norb_tot=norb_tot+n_bands_par(jproc)
   end do
   !reference orbital for process
   isband=norb_tot
   nbp=n_bands_par(iproc)



   ! Algorithm to compute the scalar product of the input guess:
   ! The term 'sqrt(bx(1)*by(2)*bz(3))' is there to normalize spherical harmonics.
   ! Wavefunctions calculated by BigDFT already are normalized.
   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!  Calculating amnk=<psi|sph_har>  !'
      write(*,*) '!==================================!'
   end if

   call timing(iproc,'CrtDescriptors','OF')
   call timing(iproc,'CrtProjectors ','ON')
   ! Initialisation :
   ! - define a new psi_etsf that can be transposed in order to keep psi_etsf in memory (psi_etsf will be used during the Mmnk matrix calculation,)
   npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*nbp,sum(comms%ncntt(0:nproc-1)))
   allocate(psi_etsf2(npsidim),stat=i_stat)
   k=0
   do i=1,nbp
      do j=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
         k=k+1
         if (i + isband_par(iproc) > n_occ ) then
            psi_etsf2(k)=psi_etsfv(j,virt_list(i+isband_par(iproc)-n_occ))
         else
            psi_etsf2(k)=psi_etsf(j,i)
         end if
      end do
   end do
  
   ! - b1, b2 and b3 are the norm of the lattice parameters.
   b1=atoms%alat1
   b2=atoms%alat2
   b3=atoms%alat3
   ! - Allocations
   allocate(amnk(n_bands,n_proj))
   allocate(amnk_tot(n_bands))
   allocate(sph_daub(npsidim2), stat=i_stat)

   ! Calculation of the spherical harmonics in parallel.
   ! It is done in the real space and then converted in the Daubechies representation.
   do npp=1, n_proj_par(iproc)
      np=npp+isproj_par(iproc)
      ! Convolution buffer : n1i=2*n1+31 -> explains the '15*input%hx*0.5' term
      r0x=ctr_proj(np,1)*b1+15*input%hx*0.5
      r0y=ctr_proj(np,2)*b2+15*input%hy*0.5
      r0z=ctr_proj(np,3)*b3+15*input%hz*0.5
      do k=1,nz
         zz=(k-1)*input%hz*0.5-r0z
         do j=1,ny
            yy=(j-1)*input%hy*0.5-r0y
            do i=1,nx
               ind=(k-1)*ny*nx+(j-1)*nx+i
               xx=(i-1)*input%hx*0.5-r0x
               call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                     xx, yy, zz, n_proj, ylm)
               call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                     xx, yy, zz, n_proj, func_r)
               ! The 'sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)' term is here to normalize spherical harmonics
               sph_har_etsf(ind)=func_r(i,j,k)*ylm(i,j,k)*sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)
            end do
         end do
      end do
      call isf_to_daub(Glr,w,sph_har_etsf(1),sph_daub(1+(npp-1)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)))
   end do

   ! Tranposition of distribution : orbitals -> components.
   call timing(iproc,'CrtProjectors ','OF')
   allocate(pwork(npsidim))
   call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,psi_etsf2,work=pwork)
   deallocate(pwork)

   ! Tranposition of distribution : orbitals -> components.
   allocate(pwork(npsidim2))
   call transpose_v(iproc,nproc,orbs2,Glr%wfd,comms2,sph_daub,work=pwork)
   deallocate(pwork)
   call timing(iproc,'ApplyProj     ','ON')

   ! Scalar product of amnk=<sph_daub|psi> in parallel.
   nvctrp=n_vctr_par(iproc)
   call gemm('T','N',n_bands,n_proj,nvctrp,1.0_wp,psi_etsf2(1),max(1,nvctrp),&
        sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),n_bands)

   ! Construction of the whole Amnk matrix.
   call mpiallred(amnk(1,1),n_proj*n_bands,MPI_SUM,MPI_COMM_WORLD,ierr)

   if (iproc==0) then
      ! Check normalisation (the amnk_tot value must tend to 1).
      write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_tot)='
      do nb=1,n_bands
         amnk_tot(nb)=sum(amnk(nb,:)**2)
         write(*,'(I4,3x,F12.6)') nb, amnk_tot(nb)
      end do
      write(*,*) '!==================================!'
      write(*,*) '!  Calculating amnk=<psi|sph_har>  !'
      write(*,*) '!               done               !'
      write(*,*) '!==================================!'
      write(*,*)
      write(*,*)
   end if

   allocate(pwork(npsidim2))
   call timing(iproc,'CrtDescriptors','OF')
   call untranspose_v(iproc,nproc,orbs2,Glr%wfd,comms2,sph_daub,work=pwork)
   call timing(iproc,'CrtDescriptors','ON')
   deallocate(pwork)
   allocate(pwork(npsidim))
   call timing(iproc,'CrtDescriptors','OF')
   call untranspose_v(iproc,nproc,orbs,Glr%wfd,comms,psi_etsf2,work=pwork)
   call timing(iproc,'CrtDescriptors','ON')
   deallocate(pwork)
   deallocate(func_r)
   deallocate(ylm)
   deallocate(sph_har_etsf)
   deallocate(sph_daub)
   deallocate(amnk_tot)
   deallocate(amnk_bands_sorted)
   deallocate(psi_etsf2)


   ! Write the .amn file
   if (iproc == 0) call write_amn(seedname, n_bands, n_kpts, n_proj, amnk)
   deallocate(amnk)


   call timing(iproc,'ApplyProj     ','OF')
   call timing(iproc,'Input_comput  ','ON')


   ! Calculate the overlap matrix (scalar product of wavefunctions psi calculated by BigDFT and themselves)
   ! This values are written in a .mmn file that will be then used by Wannier90.
   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Calculating mmnk=<psi|psi> :   !'
      write(*,*) '!==================================!'
      write(*,*) 'The values of sqrt(mmnk_tot) check the normalization in each band.'
      write(*,*) 'They must tend to their lower limit value, which is equal to 1 :'
   end if
   
   ! Allocations
   allocate(mmnk_re(n_bands,n_bands,n_kpts*n_nnkpts))
   allocate(mmnk_im(n_bands,n_bands,n_kpts*n_nnkpts))
   allocate(mmnk_tot(n_bands,n_kpts*n_nnkpts))
   allocate(psir(nx*ny*nz))
   allocate(psir_re(nx*ny*nz))
   allocate(psir_im(nx*ny*nz))
   allocate(psi_daub_re(npsidim))
   allocate(psi_daub_im(npsidim))

   ! Copy psi_etsf into psi_etsf2
   allocate(psi_etsf2(npsidim),stat=i_stat)
   call memocc(i_stat,psi_etsf2,'psi_etsf2',subname)
   k=0
   do i=1,nbp
      do j=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
         k=k+1
         if (i > n_occ ) then
            psi_etsf2(k)=psi_etsfv(j,virt_list(i-n_occ))
         else
            psi_etsf2(k)=psi_etsf(j,i)
         end if
      end do
   end do

   ! Tranposition of distribution : orbitals -> components
   call timing(iproc,'Input_comput  ','OF')
   allocate(pwork(npsidim))
   call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,psi_etsf2,work=pwork)
   deallocate(pwork)
   call timing(iproc,'Input_comput  ','ON')

   ! Algorithm to compute the scalar product :
   do inn=1,n_kpts*n_nnkpts
      if (iproc==0) then
         write(*,*)
         write(*,'(A21,3(I4,1x))') 'k-point coordinates :', (G_vec(inn,np), np=1,3)
         write(*,'(A4,4x,A15)') 'Band', 'sqrt(mmnk_tot)='
      end if
   
      ! The scalar product to calculate is <psi|psi>, and it gives a complex result, 
      ! so it is required to calculate both real and imaginary parts. It is done by :
      ! 1- converting the Daubechies wavefunctions into real space, 
      ! 2- multiply psi by the cos(.) and sin(.) factor at each point of the real space to get real and imaginary parts,
      ! 3- convert back to the Daubechies representation for real and imaginary parts.
      do nb1=1,nbp
         if (nb1<=n_occ .and. nproc==1) call daub_to_isf(Glr,w,psi_etsf(1,nb1),psir)
         if (nb1>n_occ .and. nproc==1) call daub_to_isf(Glr,w,psi_etsfv(1,virt_list(nb1-n_occ)),psir)
         do k=1,nz
            zz=(k-1)*input%hz*0.5
            do j=1,ny
               yy=(j-1)*input%hy*0.5
               do i=1,nx
                  xx=(i-1)*input%hx*0.5
                  ind=(k-1)*ny*nx+(j-1)*nx+i
                  psir_re(ind)= psir(ind) * cos( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                  psir_im(ind)=-psir(ind) * sin( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
               end do
            end do
         end do
         call isf_to_daub(Glr,w,psir_re(1),psi_daub_re(1+(nb1-1)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)))
         call isf_to_daub(Glr,w,psir_im(1),psi_daub_im(1+(nb1-1)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)))
      end do

      ! Tranposition of distribution : orbitals -> components
      call timing(iproc,'Input_comput  ','OF')
      allocate(pwork(npsidim))
      call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,psi_daub_re,work=pwork)
      call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,psi_daub_im,work=pwork)
      deallocate(pwork)
      call timing(iproc,'Input_comput  ','ON')

      ! Scalar product to compute the overlap matrix
      nvctrp=comms%nvctr_par(iproc,1)
      allocate(mmnk_v_re(n_bands*n_bands))
      allocate(mmnk_v_im(n_bands*n_bands))
      call gemm('T','N',n_bands,n_bands,nvctrp,1.0_wp,psi_daub_re(1),max(1,nvctrp),&
           psi_etsf2(1),max(1,nvctrp),0.0_wp,mmnk_v_re(1),n_bands)
      call gemm('T','N',n_bands,n_bands,nvctrp,1.0_wp,psi_daub_im(1),max(1,nvctrp),&
           psi_etsf2(1),max(1,nvctrp),0.0_wp,mmnk_v_im(1),n_bands)

      ! Reduce the overlap matrix between all the processors
      if (nproc > 1) then
         call mpiallred(mmnk_v_re(1),n_bands*n_bands,MPI_SUM,MPI_COMM_WORLD,ierr)
         call mpiallred(mmnk_v_im(1),n_bands*n_bands,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if

      ! Reshape the overlap matrix elements into a more manageable disposition
      mmnk_re(:,:,inn)=reshape(mmnk_v_re,(/n_bands,n_bands/))
      mmnk_im(:,:,inn)=reshape(mmnk_v_im,(/n_bands,n_bands/))
      deallocate(mmnk_v_re)
      deallocate(mmnk_v_im)

      ! Check the normalisation
      do nb1=1, n_bands
         mmnk_tot(nb1,inn)=sum(mmnk_re(nb1,:,inn)**2+mmnk_im(nb1,:,inn)**2)
         if (iproc==0) write(*,'(I4,3x,F12.6)') nb1, sqrt(mmnk_tot(nb1,inn))
      end do
   end do

   deallocate(mmnk_tot)
   deallocate(psi_etsf2)
   deallocate(psir)
   deallocate(psir_re)
   deallocate(psir_im)
   deallocate(psi_daub_re)
   deallocate(psi_daub_im)

   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '! Calculating mmnk=<psi|psi> done  !'
      write(*,*) '!==================================!'
      write(*,*)
      write(*,*)
   end if

   ! Write the .mmn file
   if (iproc==0) call write_mmn(seedname, n_bands, n_kpts, n_nnkpts, k_plus_b, G_vec, mmnk_re, mmnk_im)
   deallocate(mmnk_re)
   deallocate(mmnk_im)

   call timing(iproc,'Input_comput  ','OF')

   if(nproc > 1) then
      call timing(iproc,'parallel      ','RE')
   else
      call timing(iproc,'             ','RE')
   end if

!######################################################################################################################################################################
! CUBE TREATMENT BEGINS HERE
! TO DO: Redefine the timing
!        Do the Parallelization
!######################################################################################################################################################################

else if ( (filetype == 'cube' .or. filetype == 'CUBE') .and. nproc==1 ) then


! Read integers in order to allocate tables used to store
! wavefunctions, l, mr, rvalue, zona, ...

   call read_cube_header_1(nx, ny, nz, n_at, bx, by, bz)
   allocate(Z(n_at)) 
   allocate(at_pos(n_at,3))
   call read_cube_header_2(n_at, Z, at_pos)
   call read_nnkp_int_alloc(iproc,seedname, n_kpts, n_proj, n_nnkpts, n_excb)


!! Find if the values read do not reach limitations of the program
!! Stops the program and send an error message in case there are problems
!
!   call limitations(seedname, n_proj, n_occ, n_virt_tot, n_bands, n_kpts)


! Allocations

   allocate(kpts(n_kpts,3))
   allocate(ctr_proj(n_proj,3))
   allocate(x_proj(n_proj,3))
   allocate(z_proj(n_proj,3))
   allocate(l(n_proj))
   allocate(mr(n_proj))
   allocate(rvalue(n_proj))
   allocate(zona(n_proj))
   allocate(k_plus_b(n_kpts*n_nnkpts,2))
   allocate(G_vec(n_kpts*n_nnkpts,3))
   allocate(excb(n_excb))


! Read Wannier90 .nnkp file.
! The most important informations to be read are : 
!  - ctr_proj : the coordinates of the center of the projections
!  - l, mr, rvalue and zona = Z/a : the parameters used to build spherical harmonics
!  - k_plus_b and G_vec : the parameters used to build the nearest neighbours k-points

   call read_nnkp(iproc,seedname, calc_only_A, real_latt, recip_latt, n_kpts, n_proj, n_nnkpts, &
                n_excb, kpts, ctr_proj, x_proj, z_proj, l, mr, rvalue, &
                zona, k_plus_b, G_vec, excb)

   if (pre_check .eqv. .true.) then 
! If pre-check mode is chosen, calculation of the the scalar product of all unoccupied wavefunctions calculated by BigDFT and spherical harmonics.

      ! Define which unoccupied orbitals have to be read in pre-check mode.
      ! virt_list nominates the number of all the virtual orbitals up to n_virt_tot (integer defined in input.inter)
      n_bands=n_occ+n_virt_tot
      if (n_virt_tot .ne. 0) then
         allocate(virt_list(n_virt_tot))
         do i=1,n_virt_tot
            virt_list(i)=i
         end do
      else 
         STOP
      end if

      ! Algorithm to compute the scalar product :
      ! The term 'sqrt(bx(1)*by(2)*bz(3))' is there to normalize spherical harmonics.
      ! Wavefunctions calculated by BigDFT already are normalized.
      write(*,*) '!==================================!'
      write(*,*) '! Calculating amnk=<virt|sph_har>  !'
      write(*,*) '!       in pre-check mode :        !'
      write(*,*) '!==================================!'
      write(*,'(A12,4x,A15)') 'Virtual band', 'amnk_guess(nb)='
      ! b1, b2 and b3 are the norm of the lattice parameters
      ! bx(:), by(:) and bz(:) define the basis used in BigDFT calculation (in Cartesian)
      ! nx, ny, and nz define the number of points where wavefunctions are calculated, following each cartesian axis
      allocate(amnk(n_virt_tot,n_proj))
      allocate(amnk_guess(n_virt_tot))
      allocate(amnk_guess_sorted(n_virt))
      allocate(amnk_bands_sorted(n_virt))
      b1=sqrt(nx*bx(1)*nx*bx(1)+ny*bx(2)*ny*bx(2)+nz*bx(3)*nz*bx(3))
      b2=sqrt(nx*by(1)*nx*by(1)+ny*by(2)*ny*by(2)+nz*by(3)*nz*by(3))
      b3=sqrt(nx*bz(1)*nx*bz(1)+ny*bz(2)*ny*bz(2)+nz*bz(3)*nz*bz(3))
      do nb=1,n_virt_tot
         amnk_guess(nb)=0.d0
         allocate(virt(nx,ny,nz))
         call read_virtcube_1(nx, ny, nz, n_at, virt, nb+n_occ, virt_list, n_virt_tot, n_occ)
         do np=1, n_proj
            amnk(nb,np)=0.d0
            r0x=ctr_proj(np,1)*b1
            r0y=ctr_proj(np,2)*b2
            r0z=ctr_proj(np,3)*b3
            do k=1,nz
               if ( nb == 1 ) then
                  zz=k*bz(3)-r0z
               end if
               do j=1,ny
                  if ( nb == 1 ) then
                     yy=j*by(2)-r0y
                  end if
                  do i=1,nx
                     if ( nb == 1 ) then
                        if ( (np==1) .and. (i==1) .and. (j==1) .and. (k==1) ) then
                           allocate(ylm(nx,ny,nz))
                           allocate(func_r(nx,ny,nz))
                           allocate(sph_har(nx,ny,nz,n_proj))
                        end if
                        xx=i*bx(1)-r0x
                        call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                           xx, yy, zz, n_proj, ylm)
                        call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                           xx, yy, zz, n_proj, func_r)
                        sph_har(i,j,k,np)=func_r(i,j,k)*ylm(i,j,k)
                        if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then
                           call write_cube(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, n_proj, &
                                    nx, ny, nz, n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
                           if ( np==n_proj) then
                              deallocate(func_r)
                              deallocate(ylm)
                           end if
                        end if
                     end if
                     amnk(nb,np)=amnk(nb,np)+virt(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
                  end do
               end do
            end do
            ! sqrt(amnk_tot) verifies the normalization in each band (it must tend to 1)
            amnk_guess(nb)=amnk_guess(nb)+(amnk(nb,np))**2
         end do
         deallocate(virt)
         write(*,'(I4,11x,F12.6)') nb, amnk_guess(nb)
      end do
      deallocate(virt_list)
      write(*,*)
      write(*,'(1a)') 'These are the virtual bands to use to construct the actual Amn and Mmn matrices :'
      write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_guess)='
      do nb1=1,n_virt
         amnk_guess_sorted(nb1)=maxval(amnk_guess,n_virt_tot)
         amnk_bands_sorted(nb1)=maxloc(amnk_guess,n_virt_tot)
         amnk_guess(amnk_bands_sorted(nb1))=0.d0
      write(*,'(I4,3x,F12.6)') amnk_bands_sorted(nb1), amnk_guess_sorted(nb1)
      end do
      deallocate(amnk)
      deallocate(amnk_guess)
      deallocate(amnk_guess_sorted)
      write(*,*) '!==================================!'
      write(*,*) '! Calculating amnk=<virt|sph_har>  !'
      write(*,*) '!     in pre-check mode done       !'
      write(*,*) '!==================================!'
      write(*,*)
      write(*,*)

   ! Rewrite the input.inter file in order to leave the pre-check mode the next time the program is run
      call write_inter(n_virt, amnk_bands_sorted)

   end if

! In this part of the program, the pre-check has already been done.
! Calculate the scalar product of chosen wavefunctions calculated by BigDFT and spherical harmonics computed above.
! This values are written in a .amn file that will be then used by Wannier90.

   ! Define which unoccupied orbitals have to be read.
   ! Here, virt_list nominates the virtual orbitals chosen in pre-check mode.
   n_bands=n_occ+n_virt
   if (n_virt .ne. 0) then
      allocate(virt_list(n_virt))
      if (pre_check .eqv. .true.) virt_list(:)=amnk_bands_sorted(:)
      if (pre_check .eqv. .false.) call read_inter_list(iproc, n_virt, virt_list)
   else 
      STOP
   end if

   call timing(iproc,'CrtProjectors ','ON')

   ! Quick algorithm to compute the amnk matrix :
   ! The term 'sqrt(bx(1)*by(2)*bz(3))' is there to normalize spherical harmonics.
   ! Wavefunctions already are normalized.
   write(*,*) '!==================================!'
   write(*,*) '! Calculating amnk=<psi|sph_har> : !'
   write(*,*) '!==================================!'
   write(*,*) 'The values of sqrt(amnk_tot) check the normalization in each band.'
   write(*,*) 'They must tend to their lower limit value, which is equal to 1 :'
   write(*,'(A4,4x,A15)') 'Band', 'sqrt(amnk_tot)='
   ! b1, b2 and b3 are the norm of the lattice parameters
   ! bx(:), by(:) and bz(:) define the basis used in BigDFT calculation (in Cartesian)
   ! nx, ny, and nz define the number of points where wavefunctions are calculated, following each cartesian axis
   allocate(amnk(n_bands,n_proj))
   if (pre_check .eqv. .false.) allocate(amnk_bands_sorted(n_virt))
   allocate(amnk_tot(n_bands))
   b1=sqrt(nx*bx(1)*nx*bx(1)+ny*bx(2)*ny*bx(2)+nz*bx(3)*nz*bx(3))
   b2=sqrt(nx*by(1)*nx*by(1)+ny*by(2)*ny*by(2)+nz*by(3)*nz*by(3))
   b3=sqrt(nx*bz(1)*nx*bz(1)+ny*bz(2)*ny*bz(2)+nz*bz(3)*nz*bz(3))
   if (n_virt == n_occ) then
      do nb=1,n_occ
         amnk_tot(nb)=0.d0
         amnk_tot(nb+n_occ)=0.d0
         allocate(orb(nx,ny,nz))
         call read_orbcube_1(nx, ny, nz, n_at, orb, nb) 
         allocate(virt(nx,ny,nz))
         call read_virtcube_1(nx, ny, nz, n_at, virt, nb+n_occ, virt_list, n_virt, n_occ)
         do np=1,n_proj
            amnk(nb,np)=0.d0
            amnk(nb+n_occ,np)=0.d0
            r0x=ctr_proj(np,1)*b1
            r0y=ctr_proj(np,2)*b2
            r0z=ctr_proj(np,3)*b3
            do k=1,nz
               if ( nb == 1 ) zz=k*bz(3)-r0z
               do j=1,ny
                  if ( nb == 1 ) yy=j*by(2)-r0y
                  do i=1,nx
                     if ( nb == 1 .and. (pre_check .eqv. .false.) ) then
                        if ( (np==1) .and. (i==1) .and. (j==1) .and. (k==1) ) then
                           allocate(ylm(nx,ny,nz))
                           allocate(func_r(nx,ny,nz))
                           allocate(sph_har(nx,ny,nz,n_proj))
                        end if
                        xx=i*bx(1)-r0x
                        call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                           xx, yy, zz, n_proj, ylm)
                        call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                           xx, yy, zz, n_proj, func_r)
                        sph_har(i,j,k,np)=func_r(i,j,k)*ylm(i,j,k)
                        if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then
                           call write_cube(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, n_proj, &
                                    nx, ny, nz, n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
                           if ( np==n_proj ) then
                              deallocate(func_r)
                              deallocate(ylm)
                           end if
                        end if
                     end if
                     amnk(nb,np)=amnk(nb,np)+orb(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
                     amnk(nb+n_occ,np)=amnk(nb+n_occ,np)+virt(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
                  end do
               end do
            end do
            ! sqrt(amnk_tot) verifies the normalization in each band (it must tend to 1)
            amnk_tot(nb)=amnk_tot(nb)+(amnk(nb,np))**2
            amnk_tot(nb+n_occ)=amnk_tot(nb+n_occ)+(amnk(nb+n_occ,np))**2
         end do
      if (allocated(orb)) deallocate(orb)
      if (allocated(virt)) deallocate(virt)
      end do
      do nb=1,n_bands
         write(*,'(I4,3x,F12.6)') nb, sqrt(amnk_tot(nb))
      end do
      deallocate(amnk_tot)
      deallocate(sph_har)
   else
      do nb=1,n_bands
         amnk_tot(nb)=0.d0
         if ( (nb .gt. 0) .and. (nb .le. n_occ) ) then
            allocate(orb(nx,ny,nz))
            call read_orbcube_1(nx, ny, nz, n_at, orb, nb)        
         else
            allocate(virt(nx,ny,nz))
            call read_virtcube_1(nx, ny, nz, n_at, virt, nb, virt_list, n_virt, n_occ)
         end if
         do np=1, n_proj
            amnk(nb,np)=0.d0
            r0x=ctr_proj(np,1)*b1
            r0y=ctr_proj(np,2)*b2
            r0z=ctr_proj(np,3)*b3
            do k=1,nz
               if ( nb == 1 ) then
                  zz=k*bz(3)-r0z
               end if
               do j=1,ny
                  if ( nb == 1 ) then
                     yy=j*by(2)-r0y
                  end if
                  do i=1,nx
                     if ( nb == 1 .and. (pre_check .eqv. .false.) ) then
                        if ( (np==1) .and. (i==1) .and. (j==1) .and. (k==1) ) then
                           allocate(ylm(nx,ny,nz))
                           allocate(func_r(nx,ny,nz))
                           allocate(sph_har(nx,ny,nz,n_proj))
                        end if
                        xx=i*bx(1)-r0x
                        call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                           xx, yy, zz, n_proj, ylm)
                        call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                           xx, yy, zz, n_proj, func_r)
                        sph_har(i,j,k,np)=func_r(i,j,k)*ylm(i,j,k)
                        if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then
                           call write_cube(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, n_proj, &
                                    nx, ny, nz, n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
                           if ( np==n_proj) then
                              deallocate(func_r)
                              deallocate(ylm)
                           end if
                        end if
                     end if
                     if ( (nb .gt. 0) .and. (nb .le. n_occ) ) then
                        amnk(nb,np)=amnk(nb,np)+orb(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
                     else
                        amnk(nb,np)=amnk(nb,np)+virt(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
                     end if
                  end do
               end do
            end do
            ! sqrt(amnk_tot) verifies the normalization in each band (it must tend to 1)
            amnk_tot(nb)=amnk_tot(nb)+(amnk(nb,np))**2
         end do
      if (allocated(orb)) deallocate(orb)
      if (allocated(virt)) deallocate(virt)
      write(*,'(I4,3x,F12.6)') nb, sqrt(amnk_tot(nb))
      end do
      deallocate(amnk_tot)
      deallocate(sph_har)
   end if
   write(*,*) '!==================================!'
   write(*,*) '!  Calculating amnk=<psi|sph_har>  !'
   write(*,*) '!               done               !'
   write(*,*) '!==================================!'
   write(*,*)
   write(*,*)


   ! Write the .amn file
   call write_amn(seedname, n_bands, n_kpts, n_proj, amnk)
   deallocate(amnk)

   call timing(iproc,'CrtProjectors ','OF')
   call timing(iproc,'CrtDescriptors','ON')

   ! Calculate the scalar product of wavefunctions psi calculated by BigDFT and themselves.
   ! This values are written in a .mmn file that will be then used by Wannier90.
   ! Quick algorithm to calculate mmnk matrix
   write(*,*) '!==================================!'
   write(*,*) '!   Calculating mmnk=<psi|psi> :   !'
   write(*,*) '!==================================!'
   write(*,*) 'The values of sqrt(mmnk_tot) check the normalization in each band.'
   write(*,*) 'They must tend to their lower limit value, which is equal to 1 :'
   allocate(mmnk_re(n_bands,n_bands,n_kpts*n_nnkpts))
   allocate(mmnk_im(n_bands,n_bands,n_kpts*n_nnkpts))
   allocate(mmnk_tot(n_bands,n_kpts*n_nnkpts))
   ! Algorithm to compute the scalar product :
   do inn=1,n_kpts*n_nnkpts
      write(*,*)
      write(*,'(A21,3(I4,1x))') 'k-point coordinates :', (G_vec(inn,np), np=1,3)!G_vec(inn,1), G_vec(inn,2), G_vec(inn,3)
      write(*,'(A4,4x,A15)') 'Band', 'sqrt(mmnk_tot)='
      if (n_occ == n_virt) then 
         do nb1=1,n_occ
            mmnk_tot(nb1,inn)=0.d0
            mmnk_tot(nb1+n_occ,inn)=0.d0
            if ( (nb1 == 1) .and. (inn == 1) ) allocate(psi(nx,ny,nz,n_bands))
            do nb2=1,n_occ
               mmnk_re(nb1,nb2,inn)=0.d0
               mmnk_im(nb1,nb2,inn)=0.d0
               mmnk_re(nb1+n_occ,nb2,inn)=0.d0
               mmnk_im(nb1+n_occ,nb2,inn)=0.d0
               mmnk_re(nb1,nb2+n_occ,inn)=0.d0
               mmnk_im(nb1,nb2+n_occ,inn)=0.d0
               mmnk_re(nb1+n_occ,nb2+n_occ,inn)=0.d0
               mmnk_im(nb1+n_occ,nb2+n_occ,inn)=0.d0
               do k=1,nz
                  zz=k*bz(3)
                  do j=1,ny
                     yy=j*by(2)
                     do i=1,nx
                        xx=i*bx(1)
                        if ( (inn==1) .and. (nb1 == 1) &
                               .and. (i==1) .and. (j==1) .and. (k==1)) then
                           call read_orbcube(nx, ny, nz, n_at, psi, n_bands, nb2)
                           call read_virtcube(nx, ny, nz, n_at, psi, nb2+n_occ, virt_list, n_virt, n_bands, n_occ)
                        end if
                        mmnk_re(nb1,nb2,inn)=mmnk_re(nb1,nb2,inn)+psi(i,j,k,nb1)*psi(i,j,k,nb2)*cos( 2*pi* &
                           (xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                        mmnk_im(nb1,nb2,inn)=mmnk_im(nb1,nb2,inn)-psi(i,j,k,nb1)*psi(i,j,k,nb2)*sin( 2*pi* &
                           (xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                        mmnk_re(nb1+n_occ,nb2,inn)=mmnk_re(nb1+n_occ,nb2,inn)+psi(i,j,k,nb1+n_occ)* &
                           psi(i,j,k,nb2)*cos( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+ &
                           zz*G_vec(inn,3)/b3) )
                        mmnk_im(nb1+n_occ,nb2,inn)=mmnk_im(nb1+n_occ,nb2,inn)-psi(i,j,k,nb1+n_occ)* &
                           psi(i,j,k,nb2)*sin( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+ &
                           zz*G_vec(inn,3)/b3) )
                        mmnk_re(nb1,nb2+n_occ,inn)=mmnk_re(nb1,nb2+n_occ,inn)+psi(i,j,k,nb1)* &
                           psi(i,j,k,nb2+n_occ)*cos( 2*pi*(xx*G_vec(inn,1)/b1+ &
                           yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                        mmnk_im(nb1,nb2+n_occ,inn)=mmnk_im(nb1,nb2+n_occ,inn)-psi(i,j,k,nb1)* &
                           psi(i,j,k,nb2+n_occ)*sin( 2*pi*(xx*G_vec(inn,1)/b1+ &
                           yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                        mmnk_re(nb1+n_occ,nb2+n_occ,inn)=mmnk_re(nb1+n_occ,nb2+n_occ,inn)+psi(i,j,k,nb1+n_occ)* &
                           psi(i,j,k,nb2+n_occ)*cos( 2*pi*(xx*G_vec(inn,1)/b1+ &
                           yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                        mmnk_im(nb1+n_occ,nb2+n_occ,inn)=mmnk_im(nb1+n_occ,nb2+n_occ,inn)-psi(i,j,k,nb1+n_occ)* &
                           psi(i,j,k,nb2+n_occ)* sin( 2*pi*(xx*G_vec(inn,1)/b1+ &
                           yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                     end do
                  end do
               end do
               ! sqrt(mmnk_tot) verifies the normalization in each band (it must tend to 1)
               mmnk_tot(nb1,inn)=mmnk_tot(nb1,inn)+mmnk_re(nb1,nb2,inn)**2+mmnk_im(nb1,nb2,inn)**2
               mmnk_tot(nb1+n_occ,inn)=mmnk_tot(nb1+n_occ,inn)+mmnk_re(nb1+n_occ,nb2,inn)**2+ &
                  mmnk_im(nb1+n_occ,nb2,inn)**2
               mmnk_tot(nb1,inn)=mmnk_tot(nb1,inn)+mmnk_re(nb1,nb2+n_occ,inn)**2+ &
                  mmnk_im(nb1,nb2+n_occ,inn)**2
               mmnk_tot(nb1+n_occ,inn)=mmnk_tot(nb1+n_occ,inn)+mmnk_re(nb1+n_occ,nb2+n_occ,inn)**2 &
                  +mmnk_im(nb1+n_occ,nb2+n_occ,inn)**2
            end do
         end do
         do nb1=1,n_bands
            write(*,'(I4,3x,F12.6)') nb1, sqrt(mmnk_tot(nb1,inn))       
         end do
      else
         do nb1=1,n_bands
            mmnk_tot(nb1,inn)=0.d0
            if ( (nb1 == 1) .and. (inn == 1) ) then
               allocate(psi(nx,ny,nz,n_bands))
               call read_orbcube(nx, ny, nz, n_at, psi, n_bands, nb1)
            end if
            do nb2=1,n_bands
               mmnk_re(nb1,nb2,inn)=0.d0
               mmnk_im(nb1,nb2,inn)=0.d0
               do k=1,nz
                  zz=k*bz(3)
                  do j=1,ny
                     yy=j*by(2)
                     do i=1,nx
                        xx=i*bx(1)
                        if ( (inn==1) .and. (nb1 == 1) .and. (nb2 .ne. 1) .and. (nb2 .le. n_occ) &
                                                      .and. (i==1) .and. (j==1) .and. (k==1)) then
                           call read_orbcube(nx, ny, nz, n_at, psi, n_bands, nb2)
                        end if
                        if ( (inn==1) .and. (nb1 == 1) .and. (nb2 .gt. n_occ) &
                                   .and. (i==1) .and. (j==1) .and. (k==1)) then
                           call read_virtcube(nx, ny, nz, n_at, psi, nb2, virt_list, n_virt, n_bands, n_occ)
                        end if
                        mmnk_re(nb1,nb2,inn)=mmnk_re(nb1,nb2,inn)+psi(i,j,k,nb1)*psi(i,j,k,nb2)* &
                           cos( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                        mmnk_im(nb1,nb2,inn)=mmnk_im(nb1,nb2,inn)-psi(i,j,k,nb1)*psi(i,j,k,nb2)* &
                           sin( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                     end do
                  end do
               end do
            ! sqrt(mmnk_tot) verifies the normalization in each band (it must tend to 1)
            mmnk_tot(nb1,inn)=mmnk_tot(nb1,inn)+mmnk_re(nb1,nb2,inn)**2+mmnk_im(nb1,nb2,inn)**2
            end do
         write(*,'(I4,3x,F12.6)') nb1, sqrt(mmnk_tot(nb1,inn))
         end do
      end if
   end do
   deallocate(mmnk_tot)
   write(*,*) '!==================================!'
   write(*,*) '! Calculating mmnk=<psi|psi> done  !'
   write(*,*) '!==================================!'
   write(*,*)
   write(*,*)


   ! Write the .mmn file
   call write_mmn(seedname, n_bands, n_kpts, n_nnkpts, k_plus_b, G_vec, mmnk_re, mmnk_im)
   deallocate(mmnk_re)
   deallocate(mmnk_im)

   call timing(iproc,'CrtDescriptors','OF')
   call timing(iproc,'            ','RE')


   ! Write UNKnk.s files containing all the Bloch states calculated by BigDFT. 
   ! These files are used for plotting.
   ! Beware : it is not implemented in the case where there are more than 999 k-points
   do nk=1, n_kpts
      ! s is the spin, set by default to 1
      s=1
      call write_unk(n_bands, nx, ny, nz, nk, s, psi)
   end do


   ! Deallocations
   if(allocated(psi))                   deallocate(psi)
   if(allocated(Z))                     deallocate(Z)
   if(allocated(at_pos))                deallocate(at_pos)
   if(allocated(kpts))                  deallocate(kpts)
   if(allocated(ctr_proj))              deallocate(ctr_proj)
   if(allocated(x_proj))                deallocate(x_proj)
   if(allocated(z_proj))                deallocate(z_proj)
   if(allocated(l))                     deallocate(l)
   if(allocated(mr))                    deallocate(mr)
   if(allocated(rvalue))                deallocate(rvalue)
   if(allocated(zona))                  deallocate(zona)
   if(allocated(k_plus_b))              deallocate(k_plus_b)
   if(allocated(G_vec))                 deallocate(G_vec)
   if(allocated(excb))                  deallocate(excb)
   if(allocated(amnk_bands_sorted))     deallocate(amnk_bands_sorted)

else
   if (iproc==0) write(*,*) 'Cubic code not parallelized'
end if
   
  ! Barrier suggested by support for titane.ccc.cea.fr, before finalise.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call MPI_FINALIZE(ierr)
stop

end program BigDFT2Wannier


!>
subroutine read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_sph, w_ang, w_rad)
   
   ! This routine reads the first lines of a .inter file

   implicit none

   ! I/O variables
   integer, intent(in) :: iproc
   character, intent(out) :: seedname*16, filetype*4
   integer, intent(out) :: n_occ, n_virt, n_virt_tot
   logical, intent(out) :: w_sph, w_ang, w_rad, pre_check

   ! Local variables
   character :: char1*1, char2*1, char3*1


   OPEN(11, FILE='input.inter', STATUS='OLD')

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Reading input.inter header :   !'
      write(*,*) '!==================================!'
   end if
   w_sph=.false.
   w_ang=.false.
   w_rad=.false.
   read(11,*) seedname
   if(iproc==0)write(*,*) 'System studied : ', trim(seedname)
   read(11,*) filetype
   if(iproc==0)write(*,*) 'file type : ', filetype
   read(11,*) n_occ
   if(iproc==0)write(*,'(A30,I4)') 'Number of occupied orbitals :', n_occ
   read(11,*) char1, n_virt_tot, n_virt
   if (char1=='T') then
      pre_check=.true.
      if(iproc==0)write(*,*) 'Pre-check before calculating Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A38,I4)') 'Total number of unnocupied orbitals :', n_virt_tot
   else
      pre_check=.false.
      if(iproc==0)write(*,*) 'Calculation of Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A39,I4)') 'Number of chosen unnocupied orbitals :', n_virt
   end if
   read(11,*) char1, char2, char3
   if (char1=='T') then
      w_sph=.true.
      if(iproc==0) write(*,*) 'You want to write .cube files for spherical harmonics'
   end if
   if (char2=='T') then
      w_ang=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for angular parts of the spherical harmonics'
   end if
   if (char3=='T') then
      w_rad=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for radial parts of the spherical harmonics'
   end if
   CLOSE(11)

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '! Reading input.inter header done  !'
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_inter_header


!>
subroutine limitations(seedname, n_proj, n_occ, n_virt, n_bands, n_kpts)

   ! This routine sends a message in case there are problems with some values

   implicit none

   ! I/O variables
   character, intent(in) :: seedname*16
   integer, intent(in) :: n_proj, n_occ, n_virt, n_bands, n_kpts


   ! Beware : it is not implemented in the case where there are more than 9999 occupied bands 
   if (n_occ>9999) then
      write(*,*) 'There are too many occupied bands'
      STOP
   else
      if (n_occ<=0) then
         write(*,*) 'Wrong number of occupied bands (lower than 0)'
         STOP
      end if
   end if

   ! Beware : it is not implemented in the case where there are more than 9999 unoccupied bands 
   if (n_virt>9999) then
      write(*,*) 'There are too many virtual bands'
      STOP
   else
      if (n_virt<0) then
         write(*,*) 'Wrong number of virtual bands (lower than 0)'
         STOP
      end if
   end if

   ! Beware : it is not implemented in the case where there are more than 9999 projections 
   if (n_proj>9999) then
      write(*,*) 'There are too many projections'
      STOP
   else
      if (n_proj<0) then
         write(*,*) 'Wrong number of projections (lower than 0)'
         STOP
      end if
   end if

   ! Beware : it is not implemented in the case where there are more than 9999 k-points 
   if (n_kpts>9999) then
      write(*,*) 'There are too many k-points '
      STOP
   else
      if (n_kpts<0) then
         write(*,*) 'Wrong number of k-points (lower than 0)'
         STOP
      end if
   end if

   ! Beware : there cannot be more projections than bands
   if (n_proj > n_bands) then
      write(*,*) 'There cannot be more projections than bands'
      write(*,*) 'Possible reasons for problems :'
      write(*,*) '- Wrong specifications for the number of occupied and unoccupied orbitals in input.inter'
      write(*,*) '- Wrong specifications for the projections in ',  trim(seedname)//'.win'
      write(*,*) 'In this second case, do not forget to restart Wannier90 in the Post-processing mode before restarting&
       & this interface program'
      STOP
   end if

end subroutine limitations


!>
subroutine read_inter_list(iproc,n_virt, virt_list)
   
   ! This routine reads the list of virtual orbitals needed

   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt,iproc
   integer, dimension(n_virt), intent(out) :: virt_list

   ! Local variables
   integer :: i,j


   OPEN(11, FILE='input.inter', STATUS='OLD')

!   write(*,*) '!==================================!'
!   write(*,*) '!  Reading virtual orbitals list : !'
!   write(*,*) '!==================================!'

   do i=1,6
      read(11,*) ! Skip first lines
   end do
   read(11,*) (virt_list(j), j=1,n_virt)
   CLOSE(11)

   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!Reading virtual orbitals list done!'
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_inter_list


!>
subroutine read_cube_header_1(nx, ny, nz, n_at, bx, by, bz)
   
   ! This routine reads the first lines of a .cube file

   implicit none

   ! I/O variables
   integer, intent(out) :: nx, ny, nz, n_at
   real(kind=8), intent(out) :: bx(3), by(3), bz(3)

   ! Local variables
   integer :: i


   OPEN(11, FILE='orbital0001.cube', STATUS='OLD')

   write(*,*) '!==================================!'
   write(*,*) '!       Reading .cube header :     !'
   write(*,*) '!==================================!'

   read(11,*) ! skip first line
   read(11,*) ! skip second line
   read(11,*) n_at
   read(11,*) nx, (bx(i), i=1,3)
   read(11,*) ny, (by(i), i=1,3)
   read(11,*) nz, (bz(i), i=1,3)
   print *, 'Values read :'
   write(*,'(A5,4(I4,A5))') 'n_at=',  n_at, '; nx=', nx, '; ny=', ny, '; nz=', nz
   CLOSE(11)

end subroutine read_cube_header_1


!>
subroutine read_cube_header_2(n_at, Z, at_pos)
   
   ! This routine reads the first lines of a .cube file

   implicit none

   ! I/O variables
   integer, intent(in) :: n_at
   integer, dimension(n_at), intent(out) :: Z
   real(kind=8), dimension(n_at,3), intent(out) :: at_pos

   ! Local variables
   integer :: i,j
   real :: zero


   OPEN(11, FILE='orbital0001.cube', STATUS='OLD')

   do i=1,6
       read(11,*) ! skip first lines
   end do
   read(11,*)   (Z(i), zero, (at_pos(i,j), j=1,3), i=1,n_at)
   CLOSE(11)

   write(*,*) '!==================================!'
   write(*,*) '!    Reading .cube header done     !'
   write(*,*) '!==================================!'
   print *
   print *

end subroutine read_cube_header_2


!>
subroutine read_orbcube_1(nx, ny, nz, n_at, psi, io)
   
   ! This routine reads a .cube file

   implicit none

   ! I/O variables
   integer, intent(in) :: nx, ny, nz, n_at, io
   real(kind=8), dimension(nx, ny, nz), intent(out) :: psi

   ! Local variables
   integer :: i, j, k
   character :: subname*11, io_c*3


   ! Finds the name of the file to be read
   if (io>0 .and. io<10) then 
      write(io_c, '(i1)') io
      subname='orbital000'//io_c
   else  
      if (io>9 .and. io<100) then
         write(io_c, '(i2)') io
         subname='orbital00'//io_c
      else
         if (io>99 .and. io<1000) then
            write(io_c, '(i3)') io
            subname='orbital0'//io_c
         else
            if (io>999 .and. io<10000) then
               write(io_c, '(i4)') io
               subname='orbital'//io_c
            end if
         end if
      end if
   end if

   ! Reading of the file
   OPEN(11, FILE=subname//'.cube', STATUS='OLD')
   do i=1,6+n_at
      read(11,*)   ! skip header
   end do
   read(11,*)   (((psi(i,j,k), k=1,nz), j=1,ny), i=1,nx)
   CLOSE(11)

end subroutine read_orbcube_1


!>
subroutine read_orbcube(nx, ny, nz, n_at, psi, n_bands, io)
   
   ! This routine reads a .cube file

   implicit none

   ! I/O variables
   integer, intent(in) :: nx, ny, nz, n_at, io, n_bands
   real(kind=8), dimension(nx, ny, nz, n_bands), intent(out) :: psi

   ! Local variables
   integer :: i, j, k
   character :: subname*11, io_c*3


   ! Finds the name of the file to be read
   if (io>0 .and. io<10) then 
      write(io_c, '(i1)') io
      subname='orbital000'//io_c
   else  
      if (io>9 .and. io<100) then
         write(io_c, '(i2)') io
         subname='orbital00'//io_c
      else
         if (io>99 .and. io<1000) then
            write(io_c, '(i3)') io
            subname='orbital0'//io_c
         else
            if (io>999 .and. io<10000) then
               write(io_c, '(i4)') io
               subname='orbital'//io_c
            end if
         end if
      end if
   end if

   ! Reading of the file
   OPEN(11, FILE=subname//'.cube', STATUS='OLD')

   do i=1,6+n_at
      read(11,*)   ! skip header
   end do
   read(11,*)   (((psi(i,j,k,io), k=1,nz), j=1,ny), i=1,nx)

   CLOSE(11)

end subroutine read_orbcube


!>
subroutine read_virtcube_1(nx, ny, nz, n_at, psi, nb, virt_list, n_virt_tot, n_occ)
   
   ! This routine reads a .cube file

   implicit none

   ! I/O variables
   integer, intent(in) :: nx, ny, nz, n_at, nb, n_virt_tot, n_occ
   integer, dimension(n_virt_tot), intent (in) :: virt_list
   real(kind=8), dimension(nx, ny, nz), intent(out) :: psi

   ! Local variables
   integer :: i, j, k, iv
   character :: subname*11, vl_c*3

   iv=nb-n_occ


   ! Finds the name of the file to be read
   if (virt_list(iv)>0 .and. virt_list(iv)<10) then 
      write(vl_c, '(i1)') virt_list(iv)
      subname='virtual000'//vl_c
   else  
      if (virt_list(iv)>9 .and. virt_list(iv)<100) then
         write(vl_c, '(i2)') virt_list(iv)
         subname='virtual00'//vl_c
      else
         if (virt_list(iv)>99 .and. virt_list(iv)<1000) then
            write(vl_c, '(i3)') virt_list(iv)
            subname='virtual0'//vl_c
         else
            if (virt_list(iv)>999 .and. virt_list(iv)<10000) then
               write(vl_c, '(i4)') virt_list(iv)
               subname='virtual'//vl_c
            end if
         end if
      end if
   end if


   ! Reading of the file
   OPEN(11, FILE=subname//'.cube', STATUS='OLD')

   do i=1,6+n_at
      read(11,*)   ! skip header
   end do
   read(11,*)   (((psi(i,j,k), k=1,nz), j=1,ny), i=1,nx)

   CLOSE(11)

end subroutine read_virtcube_1


!>
subroutine read_virtcube(nx, ny, nz, n_at, psi, nb, virt_list, n_virt, n_bands, n_occ)
   
   ! This routine reads a .cube file

   implicit none

   ! I/O variables
   integer, intent(in) :: nx, ny, nz, n_at, nb, n_virt, n_bands, n_occ
   integer, dimension(n_virt), intent (in) :: virt_list
   real(kind=8), dimension(nx, ny, nz, n_bands), intent(out) :: psi

   ! Local variables
   integer :: i, j, k, iv
   character :: subname*11, vl_c*3

   iv=nb-n_occ

   ! Finds the name of the file to be read
   if (virt_list(iv)>0 .and. virt_list(iv)<10) then 
      write(vl_c, '(i1)') virt_list(iv)
      subname='virtual000'//vl_c
   else  
      if (virt_list(iv)>9 .and. virt_list(iv)<100) then
         write(vl_c, '(i2)') virt_list(iv)
         subname='virtual00'//vl_c
      else
         if (virt_list(iv)>99 .and. virt_list(iv)<1000) then
            write(vl_c, '(i3)') virt_list(iv)
            subname='virtual0'//vl_c
         else
            if (virt_list(iv)>999 .and. virt_list(iv)<10000) then
               write(vl_c, '(i4)') virt_list(iv)
               subname='virtual'//vl_c
            end if
         end if
      end if
   end if


   ! Reading of the file
   OPEN(11, FILE=subname//'.cube', STATUS='OLD')

   do i=1,6+n_at
      read(11,*)   ! skip header
   end do
   read(11,*)   (((psi(i,j,k,nb), k=1,nz), j=1,ny), i=1,nx)

   CLOSE(11)

end subroutine read_virtcube


!>
subroutine read_nnkp_int_alloc(iproc, seedname, n_kpts, n_proj, n_nnkpts, n_excb)
   
   ! This routine reads integers used to allocate and also verifies the nnkp file is well written

   implicit none

   ! I/O variables
   integer, intent(in) :: iproc
   integer, intent(out) :: n_kpts, n_proj, n_nnkpts, n_excb
   character *16 :: seedname

   ! Local variables
   integer :: i
   character *16 :: char1, char2, char3, char4, char5


   OPEN(11, FILE=trim(seedname)//'.nnkp', STATUS='OLD')

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .nnkp integers :     !' 
      write(*,*) '!==================================!'
   end if

   read(11,*) ! skip first line

!=====calc_only_A=====!
   read(11,*) char1, char2, char3
   if ( (char1 .ne. 'calc_only_A') .or. ( (char3 .ne. 'T') .and. (char3 .ne. 'F') ) ) STOP

!=====real_lattice=====!
   read(11,*) char1, char2
   if ( (char1 .eq. 'begin') .and. (char2 .eq. 'real_lattice') ) then
      do i=1,3
         read(11,*)   ! skip real lattice coordinates
      end do
      read(11,*) char3, char4
      if ( (char3 .ne. 'end') .or. (char4 .ne. 'real_lattice') ) STOP
   else
      STOP
   end if

!=====recip_lattice=====!
   read(11,*) char1, char2
   if ( (char1 .eq. 'begin') .and. (char2.eq.'recip_lattice') ) then
      do i=1,3
         read(11,*)   ! skip reciprocal lattice coordinates
      end do
      read(11,*) char3, char4
      if ( (char3 .ne. 'end') .or. (char4 .ne. 'recip_lattice') ) STOP
   else
      STOP
   end if

!=====kpoints=====!
   read(11,*) char1, char2
   if ( (char1 .eq. 'begin') .and. (char2.eq.'kpoints') ) then
      read(11,*) char3
      if (char3 .ne. 'end') then ! verify that there are kpoints
         BACKSPACE 11
         read(11,*) n_kpts
         do i=1, n_kpts
            read(11,*)  ! skip kpoints coordinates
         end do
         read(11,*) char4, char5
         if ( (char4 .ne. 'end') .or. (char5 .ne. 'kpoints') ) STOP
      else 
         n_kpts=0
      end if
   else
      STOP
   end if

!=====projections=====!
   read(11,*) char1, char2
   if ( (char1 .eq. 'begin') .and. (char2 .eq. 'projections') ) then
      read(11,*) char3
      if (char3 .ne. 'end') then ! verify that there are projections
         BACKSPACE 11
         read(11,*) n_proj
         do i=1,2*n_proj
            read(11,*)   ! skip projection arguments
         end do
         read(11,*) char4, char5
         if ( (char4 .ne. 'end') .or. (char5 .ne. 'projections') ) STOP
      else 
         n_proj=0
      end if
   else
      STOP
   end if

!=====nnkpts=====!
   read(11,*) char1, char2
   if ( (char1 .eq. 'begin') .and. (char2 .eq. 'nnkpts') ) then
      read(11,*) char3
      if (char3 .ne. 'end') then ! verify that there are nnkpts
         BACKSPACE 11
         read(11,*) n_nnkpts
         do i=1,n_kpts*n_nnkpts
            read(11,*)   ! skip nearest neighours arguments
         end do
         read(11,*) char4, char5
         if ( (char4 .ne. 'end') .or. (char5 .ne. 'nnkpts') ) STOP
      else
         n_nnkpts=0
      end if
   else 
      STOP
   end if

!=====exclude_bands=====!
   read(11,*) char1, char2
   if ( (char1 .eq. 'begin') .and. (char2.eq.'exclude_bands') ) then
      read(11,*) n_excb
      if (n_excb .ne. 0) then ! verify that there are exclude_bands
         do i=1,n_excb
            read(11,*)   ! skip exclude bands
         end do
         read(11,*) char3, char4
         if ( (char3 .ne. 'end') .or. (char4 .ne. 'exclude_bands') ) STOP
      else 
         n_excb=0
      end if
   else
      STOP
   end if

   close(11)
   if(iproc==0) then
      print *, 'Values read :'
      write(*,'(4(A12,I4))') '     n_kpts=',  n_kpts, ';    n_proj=', n_proj, ';  n_nnkpts=', n_nnkpts, ';    n_excb=', n_excb
      write(*,*) '!==================================!'
      write(*,*) '!   Reading .nnkp integers  done   !' 
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_nnkp_int_alloc


!>
subroutine read_nnkp(iproc,seedname, calc_only_A, real_latt, recip_latt, n_kpts, n_proj, &
              n_nnkpts, n_excb, kpts, ctr_proj, x_proj, z_proj, &
              l, mr, rvalue, zona, k_plus_b, G_vec, excb)

   ! This routine reads an .nnkp file

   implicit none

   ! I/O variables
   character *16 :: seedname
   logical, intent(out) :: calc_only_A
   real, intent(out) :: real_latt(3,3), recip_latt(3,3)
   integer, intent(in) :: n_kpts, n_proj, n_nnkpts, n_excb,iproc
   real, dimension(n_kpts,3), intent(out) :: kpts
   real, dimension(n_proj,3), intent(out) :: ctr_proj, x_proj, z_proj
   real, dimension(n_proj), intent(out) :: zona
   integer, dimension(n_proj), intent(out) :: l, mr, rvalue
   integer, dimension(n_nnkpts,2), intent(out) :: k_plus_b
   integer, dimension(n_nnkpts,3), intent(out) :: G_vec
   integer, dimension(n_excb), intent(out) :: excb

   ! Local variables
   integer :: i, j
   character *16 :: char1, char2, char3, char4


   OPEN(11, FILE=trim(seedname)//'.nnkp', STATUS='OLD')

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!       Reading .nnkp file :       !'
      write(*,*) '!==================================!'
   end if

   READ(11,*) ! skip first line


!=====calc_only_A=====!
   READ(11,*) char1, char2, char3
   if (char3 .eq. 'T') then 
      calc_only_A=.TRUE.
      if(iproc==0) write(*,*) 'There needs to only calculate A'
   else 
      if (char3 .eq. 'F') then 
         calc_only_A=.FALSE.
         if(iproc==0) write(*,*) 'Both A and M matrices must be calculated'
      end if
   end if

!=====real_lattice=====!
   READ(11,*) char1, char2 ! skip "begin real_lattice"
   READ(11,*) ((real_latt(i,j), j=1,3), i=1,3)
   READ(11,*) char3, char4 ! skip "end real_lattice"
   if(iproc==0) write(*,*) '3 lines read, corresponding to real lattice coordinates'

!=====recip_lattice=====!
   READ(11,*) char1, char2 ! skip "begin recip_lattice"
   READ(11,*) ((recip_latt(i,j), j=1,3), i=1,3)
   READ(11,*) char3, char4 ! skip "end recip_lattice"
   if(iproc==0) write(*,*) '3 lines read, corresponding to reciprocal lattice coordinates'

!=====kpoints=====!
   READ(11,*) char1, char2 ! skip "begin kpoints"
   READ(11,*) char3
   if (char3 .ne. 'end') then ! verify that there are kpoints
      BACKSPACE 11
      READ(11,*) ! skip n_kpts
      READ(11,*) ((kpts(i,j), j=1,3), i=1,n_kpts)
      READ(11,*) char3, char4 ! skip "end kpoints"
      if(iproc==0) write(*,'(I4,A50)') n_kpts, 'lines read, corresponding to kpoints coordinates'
   end if

!=====projections=====!
   READ(11,*) char1, char2 ! skip "begin projections"
   READ(11,*) char3
   if (char3 .ne. 'end') then ! verify that there are projections
      BACKSPACE 11
      READ(11,*) ! skip n_proj
      READ(11,*) ((ctr_proj(i,j), j=1,3), l(i), mr(i), rvalue(i), (z_proj(i,j), j=1,3), (x_proj(i,j), j=1,3), zona(i), i=1,n_proj)
      READ(11,*) char3, char4 ! skip "end projections"
      if(iproc==0) write(*,'(I4,A52)') 2*n_proj, 'lines read, corresponding to projections arguments'
   end if

!=====nnkpts=====!
   READ(11,*) char1, char2 ! skip "begin nnkpts"
   READ(11,*) char3
   if (char3 .ne. 'end') then ! verify that there are nnkpts
      BACKSPACE 11
      READ(11,*) ! skip n_nnkpts
      READ(11,*) ((k_plus_b(i,j), j=1,2), (G_vec(i,j), j=1,3), i=1,n_kpts*n_nnkpts)
      READ(11,*) char3, char4 ! skip "end nnkpts"
      if(iproc==0) write(*,'(I4,A59)') n_nnkpts, 'lines read, corresponding to nearest neighbours arguments'
   end if

!=====exclude_bands=====!
   READ(11,*) char1, char2 ! skip "begin exclude_bands"
   READ(11,*) ! skip n_excb
   if (n_excb .ne. 0) then ! verify that there are exclude_bands
      READ(11,*) (excb(i), i=1,n_excb)
      READ(11,*) char3, char4 ! skip "end exclude_bands"
      if(iproc==0) write(*,'(I4,A47)') n_excb, 'lines read, corresponding to the exclude band'
   end if

   close(11)
   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .nnkp file done      !'
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_nnkp


!>
subroutine angularpart(l, mr, np, nx, ny, nz, ix, iy, iz, &
                    xx, yy, zz, n_proj, ylm)

   ! This routine returns the angular part of the spherical harmonic identified by indices (l,mr)
   ! Calculations are made in spherical coordinates

   implicit none

   ! I/O variables
   integer, intent(in) :: l(n_proj), mr(n_proj)
   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
   real(kind=8), intent(in) :: xx, yy, zz
   real(kind=8), dimension(nx,ny,nz), intent(out) :: ylm

   ! local variables
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8), parameter :: eps8  = 1.0e-8
   real(kind=8), external :: s, pz, px, py, dz2, dxz, dyz, dx2my2, dxy, &
                          fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(kind=8) :: rr, cost, phi
   real(kind=8) :: bs2, bs3, bs6, bs12

   bs2 = 1.d0/sqrt(2.d0)
   bs3 = 1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)

   if (l(np) > 3 .OR. l(np) < -5 ) then 
      write(*,*) 'error, l out of range '
   else
      if (l(np)>=0) then
         if (mr(np) < 1 .OR. mr(np) > 2*l(np)+1) then
            write(*,*) 'error, mr out of range'
         end if
      else
         if (mr(np) < 1 .OR. mr(np) > abs(l(np))+1 ) then 
            write(*,*) 'error, mr out of range'
         end if
      end if
   end if

   rr=sqrt(xx*xx+yy*yy+zz*zz)

   cost = zz / rr

   if (xx > eps8) then
      phi = atan( yy/xx )
   else if (xx < -eps8 ) then
      phi = atan( yy/xx ) + pi
   else
      phi = sign( pi/2.d0,yy )
   end if

   if (l(np)==0) then   ! s orbital
      ylm(ix,iy,iz) = s(cost,phi)  
   end if
    
   if (l(np)==1) then   ! p orbitals
      if (mr(np)==1) ylm(ix,iy,iz) = pz(cost,phi) 
      if (mr(np)==2) ylm(ix,iy,iz) = px(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = py(cost,phi)
   end if

   if (l(np)==2) then   ! d orbitals
      if (mr(np)==1) ylm(ix,iy,iz) = dz2(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = dxz(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = dyz(cost,phi)
      if (mr(np)==4) ylm(ix,iy,iz) = dx2my2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) = dxy(cost,phi)
   endif

   if (l(np)==3) then   ! f orbitals
      if (mr(np)==1) ylm(ix,iy,iz) = fz3(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = fxz2(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = fyz2(cost,phi)
      if (mr(np)==4) ylm(ix,iy,iz) = fzx2my2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) = fxyz(cost,phi)
      if (mr(np)==6) ylm(ix,iy,iz) = fxx2m3y2(cost,phi)
      if (mr(np)==7) ylm(ix,iy,iz) = fy3x2my2(cost,phi)
   endif

   if (l(np)==-1) then  !  sp hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = bs2 * ( s(cost,phi) + px(cost,phi) ) 
      if (mr(np)==2) ylm(ix,iy,iz) = bs2 * ( s(cost,phi) - px(cost,phi) ) 
   end if

   if (l(np)==-2) then  !  sp2 hybrids 
      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi) 
   end if

   if (l(np)==-3) then  !  sp3 hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+pz(cost,phi))
      if (mr(np)==2) ylm(ix,iy,iz) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-pz(cost,phi))
      if (mr(np)==3) ylm(ix,iy,iz) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-pz(cost,phi))
      if (mr(np)==4) ylm(ix,iy,iz) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+pz(cost,phi))
   end if

   if (l(np)==-4) then  !  sp3d hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi) 
      if (mr(np)==4) ylm(ix,iy,iz) = bs2*pz(cost,phi)+bs2*dz2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) =-bs2*pz(cost,phi)+bs2*dz2(cost,phi)
   end if

   if (l(np)==-5) then  ! sp3d2 hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
      if (mr(np)==4) ylm(ix,iy,iz) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) = bs6*s(cost,phi)-bs2*pz(cost,phi)+bs3*dz2(cost,phi)
      if (mr(np)==6) ylm(ix,iy,iz) = bs6*s(cost,phi)+bs2*pz(cost,phi)+bs3*dz2(cost,phi)
   end if

   return

end subroutine angularpart


!======== l = 0 =====================================================================
function s(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: s, cost, phi
   s = 1.d0/ sqrt(4*pi)
   return
end function s


!======== l = 1 =====================================================================
function pz(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: pz, cost, phi
   pz =  sqrt(3.d0/(4*pi)) * cost
   return
end function pz

function px(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/(4*pi)) * sint * cos(phi)
   return
end function px

function py(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/(4*pi)) * sint * sin(phi)
   return
end function py


!======== l = 2 =====================================================================
function dz2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: dz2, cost, phi
   dz2 =  sqrt(1.25d0/(4*pi)) * (3.d0* cost*cost-1.d0)
   return
end function dz2

function dxz(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) ::  dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/(4*pi)) * sint*cost * cos(phi)
   return
end function dxz

function dyz(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/(4*pi)) * sint*cost * sin(phi)
   return
end function dyz

function dx2my2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/(4*pi)) * sint*sint * cos(2.d0*phi)
   return
end function dx2my2

function dxy(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/(4*pi)) * sint*sint * sin(2.d0*phi)
   return
end function dxy


!======== l = 3 =====================================================================
function fz3(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fz3, cost, phi
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   return
end function fz3

function fxz2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   return
end function fxz2

function fyz2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   return
end function fyz2

function fzx2my2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   return
end function fzx2my2

function fxyz(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   return
end function fxyz

function fxx2m3y2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   return
end function fxx2m3y2

function fy3x2my2(cost,phi)
   implicit none
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8) :: fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   return
end function fy3x2my2


!!>
!subroutine angularpart(l, mr, np, nx, ny, nz, ix, iy, iz, &
!                    xx, yy, zz, n_proj, ylm)
!
!   ! This routine returns the angular part of the spherical harmonic identified by indices (l,mr)
!   ! Calcutations are made in Cartesian coordinates
!
!   implicit none
!
!   ! I/O variables
!   integer, intent(in) :: l(n_proj), mr(n_proj)
!   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
!   real(kind=8), intent(in) :: xx, yy, zz
!   real(kind=8), dimension(nx,ny,nz), intent(out) :: ylm
!
!   ! local variables
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8), parameter :: eps8  = 1.0e-8
!   real(kind=8), external :: s, pz, px, py, dz2, dxz, dyz, dx2my2, dxy, &
!                          fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
!   real(kind=8) :: rr
!   real(kind=8) :: bs2, bs3, bs6, bs12
!
!   bs2 = 1.d0/sqrt(2.d0)
!   bs3 = 1.d0/sqrt(3.d0)
!   bs6 = 1.d0/sqrt(6.d0)
!   bs12 = 1.d0/sqrt(12.d0)
!
!
!   if (l(np) > 3 .OR. l(np) < -5 ) then 
!      write(*,*) 'error, l out of range '
!   else
!      if (l(np)>=0) then
!         if (mr(np) < 1 .OR. mr(np) > 2*l(np)+1) then
!	    write(*,*) 'error, mr out of range'
!	 end if
!      else
!         if (mr(np) < 1 .OR. mr(np) > abs(l(np))+1 ) then 
!	    write(*,*) 'error, mr out of range'
!         end if
!      end if
!   end if
!
!   rr = sqrt( xx*xx + yy*yy + zz*zz )
!    
!   if (l(np)==0) then   ! s orbital
!      ylm(ix,iy,iz) = s(xx,yy,zz,rr)  
!   end if
!
!   if (l(np)==1) then   ! p orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = pz(xx,yy,zz,rr) 
!      if (mr(np)==2) ylm(ix,iy,iz) = px(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = py(xx,yy,zz,rr)
!   end if
!
!   if (l(np)==2) then   ! d orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = dz2(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = dxz(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = dyz(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = dx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = dxy(xx,yy,zz,rr)
!   endif
!
!   if (l(np)==3) then   ! f orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = fz3(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = fxz2(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = fyz2(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = fzx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = fxyz(xx,yy,zz,rr)
!      if (mr(np)==6) ylm(ix,iy,iz) = fxx2m3y2(xx,yy,zz,rr)
!      if (mr(np)==7) ylm(ix,iy,iz) = fy3x2my2(xx,yy,zz,rr)
!   endif
!
!   if (l(np)==-1) then  !  sp hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs2 * ( s(xx,yy,zz,rr) + px(xx,yy,zz,rr) ) 
!      if (mr(np)==2) ylm(ix,iy,iz) = bs2 * ( s(xx,yy,zz,rr) - px(xx,yy,zz,rr) ) 
!   end if
!
!   if (l(np)==-2) then  !  sp2 hybrids 
!      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr) +2.d0*bs6*px(xx,yy,zz,rr) 
!   end if
!
!   if (l(np)==-3) then  !  sp3 hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)+px(xx,yy,zz,rr)+py(xx,yy,zz,rr)+pz(xx,yy,zz,rr))
!      if (mr(np)==2) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)+px(xx,yy,zz,rr)-py(xx,yy,zz,rr)-pz(xx,yy,zz,rr))
!      if (mr(np)==3) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)-px(xx,yy,zz,rr)+py(xx,yy,zz,rr)-pz(xx,yy,zz,rr))
!      if (mr(np)==4) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)-px(xx,yy,zz,rr)-py(xx,yy,zz,rr)+pz(xx,yy,zz,rr))
!   end if
!
!   if (l(np)==-4) then  !  sp3d hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr) +2.d0*bs6*px(xx,yy,zz,rr) 
!      if (mr(np)==4) ylm(ix,iy,iz) = bs2*pz(xx,yy,zz,rr)+bs2*dz2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) =-bs2*pz(xx,yy,zz,rr)+bs2*dz2(xx,yy,zz,rr)
!   end if
!
!   if (l(np)==-5) then  ! sp3d2 hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*px(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)+.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*px(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)+.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)-.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)-.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*pz(xx,yy,zz,rr)+bs3*dz2(xx,yy,zz,rr)
!      if (mr(np)==6) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*pz(xx,yy,zz,rr)+bs3*dz2(xx,yy,zz,rr)
!   end if
!
!   return
!
!end subroutine angularpart
!
!
!! The following functions are used to calculate angular parts of the spherical harmonics
!!======== l = 0 =====================================================================
!function s(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) :: s, xx, yy, zz, rr
!   s = 1.d0/ sqrt(4*pi)
!   return
!end function s
!
!
!!======== l = 1 =====================================================================
!function pz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::pz, xx, yy, zz, rr
!   pz =  sqrt(3.d0/(4*pi)) * (zz/rr)
!   return
!end function pz
!
!function px(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::px, xx, yy, zz, rr
!   px =  sqrt(3.d0/(4*pi)) * (xx/rr)
!   return
!end function px
!
!function py(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::py, xx, yy, zz, rr
!   py =  sqrt(3.d0/(4*pi)) * (yy/rr)
!   return
!end function py
!
!
!!======== l = 2 =====================================================================
!function dz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dz2, xx, yy, zz, rr
!   dz2 =  sqrt(1.25d0/(4*pi)) * (-xx*xx-yy*yy+2.d0*zz*zz)/(rr*rr)
!   return
!end function dz2
!
!function dxz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dxz, xx, yy, zz, rr
!   dxz =  sqrt(15.d0/(4*pi)) * (xx*zz)/(rr*rr)
!   return
!end function dxz
!
!function dyz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dyz, xx, yy, zz, rr
!   dyz =  sqrt(15.d0/(4*pi)) * (yy*zz)/(rr*rr)
!   return
!end function dyz
!
!function dx2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dx2my2, xx, yy, zz, rr
!   dx2my2 =  sqrt(3.75d0/(4*pi)) * (xx*xx-yy*yy)/(rr*rr)
!   return
!end function dx2my2
!
!function dxy(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dxy, xx, yy, zz, rr
!   dxy =  sqrt(3.75d0/(4*pi)) * (xx*yy)/(rr*rr)
!   return
!end function dxy
!
!
!!======== l = 3 =====================================================================
!function fz3(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fz3, xx, yy, zz, rr
!   fz3 =  0.25d0*sqrt(7.d0/pi) * (zz*(2.d0*zz*zz-3.d0*xx*xx-3.d0*yy*yy)) / (rr*rr*rr)
!   return
!end function fz3
!
!function fxz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxz2, xx, yy, zz, rr
!   fxz2 =  0.25d0*sqrt(10.5d0/pi) * (xx*(4.d0*zz*zz-xx*xx-yy*yy)) / (rr*rr*rr)
!   return
!end function fxz2
!
!function fyz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fyz2, xx, yy, zz, rr
!   fyz2 =  0.25d0*sqrt(10.5d0/pi) * (yy*(4.d0*zz*zz-xx*xx-yy*yy)) / (rr*rr*rr)
!   return
!end function fyz2
!
!function fzx2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fzx2my2, xx, yy, zz, rr
!   fzx2my2 =  0.25d0*sqrt(105d0/pi) * (zz*(xx*xx-yy*yy)) / (rr*rr*rr)
!   return
!end function fzx2my2
!
!function fxyz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxyz, xx, yy, zz, rr
!   fxyz =  0.25d0*sqrt(105d0/pi) * (xx*yy*zz) / (rr*rr*rr)
!   return
!end function fxyz
!
!function fxx2m3y2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxx2m3y2, xx, yy, zz, rr
!   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * (xx*(xx*xx-3.d0*yy*yy)) / (rr*rr*rr)
!   return
!end function fxx2m3y2
!
!function fy3x2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fy3x2my2, xx, yy, zz, rr
!   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * (yy*(3.d0*xx*xx-yy*yy)) / (rr*rr*rr)
!   return
!end function fy3x2my2


!>
subroutine radialpart(rvalue, zona, np, nx, ny, nz, ix, iy, iz, &
                    xx, yy, zz, n_proj, func_r)

   ! This routine returns the radial part of the spherical harmonic identified by the indice rvalue.

   implicit none

   ! I/O variables
   integer, intent(in) :: rvalue(n_proj)
   real, intent(in) :: zona(n_proj)
   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
   real(kind=8), intent(in) :: xx, yy, zz
   real(kind=8), dimension(nx,ny,nz), intent(out) :: func_r

   ! local variables 
   real(kind=8), parameter :: eps8  = 1.0e-8
   real(kind=8) :: rr

   rr = sqrt( xx*xx + yy*yy + zz*zz )

   if (rr < eps8) then
      write(*,*) 'rr too small '
   end if
   if (rvalue(np)==1) func_r(ix,iy,iz) = 2.d0 * (zona(np)/0.529177208)**(3.d0/2.d0) * exp(-(zona(np)/0.529177208)*rr)
   if (rvalue(np)==2) func_r(ix,iy,iz) = 1.d0/sqrt(8.d0) * (zona(np)/0.529177208)**(3.d0/2.d0) * & 
      (2.0d0 - (zona(np)/0.529177208)*rr) * exp(-(zona(np)/0.529177208)*rr*0.5d0)
   if (rvalue(np)==3) func_r(ix,iy,iz) = sqrt(4.d0/27.d0) * (zona(np)/0.529177208)**(3.0d0/2.0d0) * &
      (1.d0 - (2.0d0/3.0d0)*(zona(np)/0.529177208)*rr + 2.d0*((zona(np)/0.529177208)*rr)**2/27.d0) * &
      exp(-(zona(np)/0.529177208)*rr/3.0d0)

   return

end subroutine radialpart


!>
subroutine write_cube(w_sph, w_ang, w_rad, fn1, fn2, fn3, np, n_proj, nx, ny, nz, &
           n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)

   ! This routine writes .cube files for radial part, angular part and
   ! the spherical harmonic given in argument

   implicit none

   ! I/O variables
   logical, intent(in) :: w_sph, w_ang, w_rad
   character, intent(in) :: fn1*7, fn2*6, fn3*3
   integer, intent(in) :: np, nx, ny, nz, n_at, n_proj
   real(kind=8), dimension(3), intent(in) :: bx, by, bz
   integer, dimension(n_at), intent(in) :: Z
   real(kind=8), dimension(n_at,3), intent(in) :: at_pos
   real(kind=8), dimension(nx,ny,nz,n_proj), intent(in) :: sph_har
   real(kind=8), dimension(nx,ny,nz), intent(in) :: func_r, ylm

   ! Local variables
   character :: subname1*11, subname2*10, subname3*7, np_c*3
   integer :: i, j, rem, ix, iy, iz


   ! Concatenations to write the names of .cube files
   if (np>0 .and. np<10) then 
      write(np_c, '(i1)') np
      subname1=fn1//'000'//np_c
      subname2=fn2//'000'//np_c
      subname3=fn3//'000'//np_c
   else  
      if (np>9 .and. np<100) then
         write(np_c, '(i2)') np
         subname1=fn1//'00'//np_c
         subname2=fn2//'00'//np_c
         subname3=fn3//'00'//np_c
      else
         if (np>99 .and. np<1000) then
            write(np_c, '(i3)') np
            subname1=fn1//'0'//np_c
            subname2=fn2//'0'//np_c
            subname3=fn3//'0'//np_c
         else
            if (np>999 .and. np<10000) then
               write(np_c, '(i4)') np
               subname1=fn1//np_c
               subname2=fn2//np_c
               subname3=fn3//np_c
            end if
         end if
      end if
   end if

   rem=nz-floor(nz/6.d0)*6

   ! Write the sph_harxxx.cube files
   if (w_sph .eqv. .true.) then
      OPEN(12, FILE=subname1//'.cube', STATUS='unknown')
      write(12,*) ' CUBE file for ISF field'
      write(12,*) ' Case for'
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') n_at, real(0.d0), real(0.d0), real(0.d0)
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nx, bx(1), bx(2), bx(3)
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') ny, by(1), by(2), by(3)
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nz, bz(1), bz(2), bz(3)
      do i=1, n_at
         write(12,'(I4,1X,F12.6,3(1X,F12.6))') Z(i), real(0.d0), (real(at_pos(i,j)), j=1,3)
      end do
      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
      do ix=nx,1,-1
         do iy=ny,1,-1
            do iz=nz,1,-1
               write(12,'(E14.6)',advance='no') real(sph_har(ix,iy,iz,np))
               if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. nz) ) .or. (iz .eq. 1) ) then
                  write(12,'(a)') ''
               end if
            end do
         end do
      end do
      CLOSE(12)
   end if

   ! Write the func_rxxx.cube file
   if (w_rad .eqv. .true.) then
      OPEN(13, FILE=subname2//'.cube', STATUS='unknown')
      write(13,*) ' CUBE file for ISF field'
      write(13,*) ' Case for'
      write(13,'(I4,1X,F12.6,2(1X,F12.6))') n_at, real(0.d0), real(0.d0), real(0.d0)
      write(13,'(I4,1X,F12.6,2(1X,F12.6))') nx, bx(1), bx(2), bx(3)
      write(13,'(I4,1X,F12.6,2(1X,F12.6))') ny, by(1), by(2), by(3)
      write(13,'(I4,1X,F12.6,2(1X,F12.6))') nz, bz(1), bz(2), bz(3)
      do i=1, n_at
         write(13,'(I4,1X,F12.6,3(1X,F12.6))') Z(i), real(0.d0), (real(at_pos(i,j)), j=1,3)
      end do
      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
      do ix=nx,1,-1
         do iy=ny,1,-1
            do iz=nz,1,-1
               write(13,'(E14.6)',advance='no') real(func_r(ix,iy,iz))
               if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. nz) ) .or. (iz .eq. 1) ) then
                  write(13,'(a)') ''
               end if
            end do
         end do
      end do
      CLOSE(13)
   end if

   ! Write the ylmxxx.cube file
   if (w_ang .eqv. .true.) then
      OPEN(14, FILE=subname3//'.cube', STATUS='unknown')
      write(14,*) ' CUBE file for ISF field'
      write(14,*) ' Case for'
      write(14,'(I4,1X,F12.6,2(1X,F12.6))') n_at, real(0.d0), real(0.d0), real(0.d0)
      write(14,'(I4,1X,F12.6,2(1X,F12.6))') nx, bx(1), bx(2), bx(3)
      write(14,'(I4,1X,F12.6,2(1X,F12.6))') ny, by(1), by(2), by(3)
      write(14,'(I4,1X,F12.6,2(1X,F12.6))') nz, bz(1), bz(2), bz(3)
      do i=1, n_at
         write(14,'(I4,1X,F12.6,3(1X,F12.6))') Z(i), real(0.d0), (real(at_pos(i,j)), j=1,3)
      end do
      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
      do ix=nx,1,-1
         do iy=ny,1,-1
            do iz=nz,1,-1
               write(14,'(E14.6)',advance='no') real(ylm(ix,iy,iz))
               if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. nz) ) .or. (iz .eq. 1) ) then
                  write(14,'(a)') ''
               end if
            end do
         end do
      end do
      CLOSE(14)
  end if

end subroutine write_cube


!>
subroutine write_inter(n_virt, amnk_bands_sorted)

   ! This routine write a input.inter file
   
   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt
   integer, dimension(n_virt), intent(inout) :: amnk_bands_sorted
  
   ! Local variables
   integer :: i
   character :: seedname*20, pre_check_mode_old*1, filetype*4
   integer :: n_occ, n_virt_tot, ng(3)

   ! Read data to keep
   OPEN(11, FILE='input.inter', STATUS='OLD')
!   write(*,*) '!==================================!'
!   write(*,*) '!   Writing an input.inter file :  !'
!   write(*,*) '!==================================!'
   read(11,*) seedname
   read(11,*) filetype
   read(11,*) n_occ
   read(11,*) pre_check_mode_old, n_virt_tot
   read(11,*) ! skip line
   read(11,*) ng(1), ng(2), ng(3)
   if (pre_check_mode_old == 'F') then
      read(11,*) (amnk_bands_sorted(i), i=1,n_virt)
   end if
   CLOSE(11)

   ! Write a new input.inter file (change the value of the logical pre_check)
   OPEN(11, FILE='input.inter', STATUS='OLD')
   write(11,'(1a,1x,1a)') seedname, '# name of the .win file'
   write(11,'(1a,1x,1a)') filetype, '# cube or etsf'
   write(11,'(I4,13x,1a)') n_occ, '# no. of occupied orbitals'
   if (pre_check_mode_old == 'T') then
      write(11,'(a1,2(1x,I4),6x,1a)') 'F', n_virt_tot, n_virt, '# Pre-check, no. of virtual orbitals&
      & for pre-check, no. of virtual orbitals used to write Amn and Mmn matrices'
   else
      write(11,'(a1,2(1x,I4),6x,1a)') 'T', n_virt_tot, n_virt, '# Pre-check, no. of virtual orbitals&
      & for pre-check, no. of virtual orbitals used to write Amn and Mmn matrices'
   end if
   write(11,'(1a)') 'F    F    F      # allow the writing of .cube files for spherical& 
   & harmonics, their angular parts and their radial parts'
   write(11,'(3(I4,x),a)') (ng(i), i=1,3), '# number of points for each axis in the cubic representation&
   & used by BigDFT (needed to use Wannier90)'
   do i=1,n_virt
       write(11,'(I4, 1x)',advance='no') amnk_bands_sorted(i)
!      if (filetype=='etsf' .or. filetype=='ETSF') write(11,'(I4, 1x)',advance='no') amnk_bands_sorted(i)!+n_occ
!      if (filetype=='cube' .or. filetype=='CUBE') write(11,'(I4, 1x)',advance='no') amnk_bands_sorted(i)
      if ( (mod(i,15) .eq. 0) .and. (i .ne. n_virt) ) then
         write(11,'(a)') ''
      end if
   end do
   CLOSE(11)
   write(*,*) '!==================================!'
   write(*,*) '! Writing an input.inter file done !'
   write(*,*) '!==================================!'
   write(*,*)
   write(*,*)

end subroutine write_inter


!>
subroutine write_amn(seedname, n_bands, n_kpts, n_proj, amnk)

   ! This routine writes a .amn file that can be read by Wannier90

   implicit none

   ! I/O variables
   character, intent(in) :: seedname*16
   integer, intent(in) :: n_bands, n_kpts, n_proj
   real(kind=8), dimension(n_bands,n_proj), intent(in) :: amnk

   ! Local variables
   integer :: nb, nk, np


   ! Writing the .amn file
   OPEN(12, FILE=trim(seedname)//'.amn', STATUS='unknown')
!   write(*,*) '!==================================!'
!   write(*,*) '!       Writing a .amn file :      !'
!   write(*,*) '!==================================!'
   write(12,*) 'File Created on'
   write(12,'(I4,2(1X,I4))') n_bands, n_kpts, n_proj
   do nk=1, n_kpts
      do np=1, n_proj
         do nb=1, n_bands
            write(12,'(3(I4,1X),E13.6,1X,E13.6)') nb, np, nk, amnk(nb,np), 0.d0
         end do
      end do
   end do
   CLOSE(12)

   write(*,*) '!==================================!'
   write(*,*) '!     Writing a .amn file done     !'
   write(*,*) '!==================================!'
   write(*,*)
   write(*,*)

end subroutine write_amn


!>
subroutine write_mmn(seedname, n_bands, n_kpts, n_nnkpts, k_plus_b, G_vec, mmnk_re, mmnk_im)

   ! This routine writes a .mmn file that can be read by Wannier90

   implicit none

   ! I/O variables
   character, intent(in) :: seedname*16
   integer, intent(in) :: n_bands, n_kpts, n_nnkpts
   integer, dimension(n_kpts*n_nnkpts,2), intent(in) :: k_plus_b
   integer, dimension(n_kpts*n_nnkpts,3), intent(in) :: G_vec
   real(kind=8), dimension(n_bands,n_bands,n_kpts*n_nnkpts), intent(in) :: mmnk_re
   real(kind=8), dimension(n_bands,n_bands,n_kpts*n_nnkpts), intent(in) :: mmnk_im

   ! Local variables
   integer :: nb1, nb2, n, i, j


   ! Writing the .mmn file
   OPEN(12, FILE=trim(seedname)//'.mmn', STATUS='unknown')
!   write(*,*) '!==================================!'
!   write(*,*) '!       Writing a .mmn file :      !'
!   write(*,*) '!==================================!'
   write(12,*) 'File Created on'
   write(12,'(I4,2(1X,I4))') n_bands, n_kpts, n_nnkpts
   do n=1, n_kpts*n_nnkpts
      write(12,'(I4,4(1X,I4))') (k_plus_b(n,i), i=1,2), (G_vec(n,j), j=1,3)
      do nb2=1, n_bands
         do nb1=1, n_bands
            write(12,'(E13.6,1X,E13.6)') mmnk_re(nb1,nb2,n), mmnk_im(nb1,nb2,n)
         end do
      end do
   end do
   CLOSE(12)
   write(*,*) '!==================================!'
   write(*,*) '!     Writing a .mmn file done     !'
   write(*,*) '!==================================!'
   write(*,*)
   write(*,*)

end subroutine write_mmn


!>
subroutine write_unk(n_bands, nx, ny, nz, nk, s, psi)

   ! This routine writes a UNKnk.s used for plotting

   implicit none

   ! I/O variables
   integer, intent(in) :: n_bands, nx, ny, nz, nk, s
   real(kind=8), dimension(nx,ny,nz,n_bands), intent(in) :: psi

   ! Local variables
   integer :: nb, i, j, k
   character :: s_c*1, nk_c*3, seedname*10


   ! Concatenations to write the names of UNKnk.s files.
   if (nk>0 .and. nk<10) then 
      write(nk_c, '(i1)') nk
      write(s_c, '(i1)') s
      seedname=(('UNK0000'//trim(nk_c))//'.')//s_c
   else  
      if (nk>9 .and. nk<100) then
         write(nk_c, '(i2)') nk
         write(s_c, '(i1)') s
         seedname=(('UNK000'//trim(nk_c))//'.')//s_c
      else
         if (nk>99 .and. nk<1000) then
            write(nk_c, '(i3)') nk
            write(s_c, '(i1)') s
            seedname=(('UNK00'//trim(nk_c))//'.')//s_c
         else
            if (nk>999 .and. nk<10000) then
               write(nk_c, '(i4)') nk
               write(s_c, '(i1)') s
               seedname=(('UNK0'//trim(nk_c))//'.')//s_c
            end if
         end if
      end if
   end if

   ! Writing the UNKnk.s file
   OPEN(12, FILE=seedname, STATUS='unknown')
!   write(*,*) '!==================================!'
!   write(*,*) '!     Writing a UNKnk.s file :     !'
!   write(*,*) '!==================================!'
   write(12,'(I4,4(1X,I4))') nx, ny, nz, nk, n_bands
   do nb=1, n_bands
      do k=1, nz
         do j=1, ny
            do i=1, nx
               write(12,'(E13.6, 1X, E13.6)') psi(i,j,k,nb), 0.d0
            end do
         end do
      end do
   end do
   CLOSE(12)
   write(*,*) '!==================================!'
   write(*,*) '!    Writing a UNKnk.s file done   !'
   write(*,*) '!==================================!'
   write(*,*)
   write(*,*)

end subroutine write_unk
