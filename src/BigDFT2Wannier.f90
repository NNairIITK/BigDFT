!> @file
!!   Interface between BigDFT and Wannier90
!! @author
!!   Copyright (C) 2010-2013 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Program to calculate all quantities needed by Wannier90
program BigDFT2Wannier

   use BigDFT_API
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   use module_interfaces
   use yaml_output
   use module_input_dicts
   implicit none
   character :: filetype*4
   !etsf
   character(len=*), parameter :: subname='BigDFT2Wannier'
   type(locreg_descriptors) :: Glr
   type(orbitals_data) :: orbs  !< describes the occupied orbitals
   type(orbitals_data) :: orbsp !< describes the projectors
   type(orbitals_data) :: orbsv !< describes only the non-occupied orbitals
   type(orbitals_data) :: orbsb !< describes the distribution of all orbitals
   type(atoms_data) :: atoms
   type(input_variables) :: input 
   type(workarr_sumrho) :: w
   type(communications_arrays), target :: comms, commsp,commsv,commsb
   integer, parameter :: WF_FORMAT_CUBE = 4
   integer :: iproc, nproc, nproctiming, i_stat, ind, ierr, npsidim, npsidim2
   integer :: n_proj,nvctrp,npp,nvirtu,nvirtd,pshft,nbl1,nbl2,nbl3,iformat,info
   integer :: ncount0,ncount1,ncount_rate,ncount_max,nbr1,nbr2,nbr3,shft,wshft,lwork
   real :: tcpu0,tcpu1
   real(kind=8) ::tel
   real(kind=8) :: znorm,xnorm,ortho,ddot
   real(kind=8),parameter :: eps6=1.0d-6!, eps8=1.0d-8
   real(gp), dimension(:,:), pointer :: rxyz_old
   real(gp), dimension(:,:), allocatable :: radii_cf
   real(gp), dimension(3) :: shift
   real(wp), allocatable :: psi_etsf(:,:),psi_etsfv(:),sph_har_etsf(:),psir(:),psir_re(:),psir_im(:),sph_daub(:)
   real(wp), allocatable :: psi_daub_im(:),psi_daub_re(:),psi_etsf2(:) !!,pvirt(:)
   real(wp), allocatable :: mmnk_v_re(:), mmnk_v_im(:)
   real(wp), pointer :: pwork(:)!,sph_daub(:)
   character(len=60) :: filename, run_id
   logical :: perx, pery,perz, residentity,write_resid
   integer :: nx, ny, nz, nb, nb1, nk, inn
   real(kind=8) :: b1, b2, b3, r0x, r0y, r0z
   real(kind=8) :: xx, yy, zz
   real(kind=8), allocatable :: ylm(:,:,:), func_r(:,:,:)
   real(kind=8), allocatable :: amnk(:,:), amnk_tot(:), amnk_guess(:), amnk_guess_sorted(:),overlap_proj(:,:)
   real(kind=8), allocatable :: mmnk_re(:,:,:), mmnk_im(:,:,:), mmnk_tot(:,:)
   integer :: i, j, k, np,i_all
   character :: seedname*16
   logical :: calc_only_A 
   real, dimension(3,3) :: real_latt, recip_latt
   integer :: n_kpts, n_nnkpts, n_excb, n_at, s
   integer :: n_occ, n_virt, n_virt_tot,nconfig
   logical :: w_unk, w_sph, w_ang, w_rad, pre_check
   real, allocatable, dimension (:,:) :: kpts
   real(kind=8), allocatable, dimension (:,:) :: ctr_proj, x_proj, y_proj, z_proj
   integer, allocatable, dimension (:) :: l, mr, rvalue
   real, allocatable, dimension (:) :: zona
   integer, allocatable, dimension (:,:) :: k_plus_b
   integer, allocatable, dimension (:,:) :: G_vec
   integer, allocatable, dimension (:) :: excb,ipiv
   integer, allocatable, dimension (:) :: virt_list, amnk_bands_sorted
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   integer, dimension(4) :: mpi_info
   type(dictionary), pointer :: user_inputs

   call f_lib_initialize()
   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(mpi_info,nconfig,run_id,ierr)

   !just for backward compatibility
   iproc=mpi_info(1)
   nproc=mpi_info(2)

!!$   ! Start MPI in parallel version
!!$   !in the case of MPIfake libraries the number of processors is automatically adjusted
!!$   call MPI_INIT(ierr)
!!$   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
!!$   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!!$
!!$   call mpi_environment_set(bigdft_mpi,iproc,nproc,MPI_COMM_WORLD,0)
!!$
!!$   call memocc_set_memory_limit(memorylimit)

   if (nconfig < 0) stop 'runs-file not supported for BigDFT2Wannier executable'

   call user_dict_from_files(user_inputs, trim(run_id)//trim(bigdft_run_id_toa()), &
        & 'posinp'//trim(bigdft_run_id_toa()), bigdft_mpi)
   call inputs_from_dict(input, atoms, user_inputs, .true.)
   if (iproc == 0) then
      call print_general_parameters(input,atoms)
   end if
   call dict_free(user_inputs)

   if (input%verbosity > 2) then
      nproctiming=-nproc !timing in debug mode                                                                                                                                                                  
   else
      nproctiming=nproc
   end if

   call timing(nproctiming,'b2w_time.prc','IN')

   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max) 

   ! Read input.inter file
   call timing(iproc,'Precondition  ','ON')
   call read_inter_header(iproc,seedname, filetype, residentity, write_resid, n_occ, pre_check, n_virt_tot,&
         n_virt, w_unk, w_sph, w_ang, w_rad)

   if(n_virt_tot < n_virt) then
      if (iproc == 0) then
         call yaml_warning('The total number of virtual states,' // trim(yaml_toa(n_virt_tot)))
         call yaml_comment('is smaller than number of desired states' // trim(yaml_toa(n_virt)))
         call yaml_comment('CORRECTION: Increase total number of virtual states')
         call yaml_comment('or decrease the number of desired states')
         !write(*,'(A,1x,I4)') 'ERROR: total number of virtual states :',n_virt_tot
         !write(*,'(A,1x,I4)') 'smaller than number of desired states:',n_virt
         !write(*,'(A)') 'CORRECTION: Increase total number of virtual states'
         !write(*,'(A)') 'or decrease the number of desired states'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   ! assign the input_wf_format
   iformat = WF_FORMAT_NONE
   select case (filetype)
   case ("ETSF","etsf")
      iformat = WF_FORMAT_ETSF
   case ("BIN","bin")
      iformat = WF_FORMAT_BINARY
   case ("FORM",'form')
      iformat = WF_FORMAT_PLAIN
   case ("cube","CUBE")
      if(iproc==0) call yaml_warning('Cube files are no longer implemented. TO DO!')
      !if(iproc==0) write(*,*)'Cube files are no longer implemented. TO DO!'
      iformat = WF_FORMAT_CUBE
      call mpi_finalize(ierr)
      stop
   case default
      if (iproc == 0) call yaml_warning('Missing specification of wavefunction files')
      !if (iproc == 0) write(*,*)' WARNING: Missing specification of wavefunction files'
      stop
   end select

!!$   posinp='posinp'
!!$
!!$   ! Initalise the variables for the calculation
!!$   call standard_inputfile_names(input,radical,nproc)
!!$   call read_input_variables(iproc,nproc,posinp,input, atoms, rxyz,1,radical,1)
!!$
!!$   if (iproc == 0) call print_general_parameters(input,atoms)

   allocate(radii_cf(atoms%astruct%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)

   call system_properties(iproc,nproc,input,atoms,orbs,radii_cf)

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space
   call system_size(atoms,atoms%astruct%rxyz,radii_cf,input%crmult,input%frmult,input%hx,input%hy,input%hz,&
      &   Glr,shift)
   if (iproc == 0) &
        & call print_atoms_and_grid(Glr, atoms, atoms%astruct%rxyz, shift, input%hx,input%hy,input%hz)

   ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
   call createWavefunctionsDescriptors(iproc,input%hx,input%hy,input%hz,&
      & atoms,atoms%astruct%rxyz,radii_cf,input%crmult,input%frmult,Glr)
   if (iproc == 0) call print_wfd(Glr%wfd)

   ! don't need radii_cf anymore
   i_all = -product(shape(radii_cf))*kind(radii_cf)
   deallocate(radii_cf,stat=i_stat)
   call memocc(i_stat,i_all,'radii_cf',subname)

   ! Allocate communications arrays (allocate it before Projectors because of the definition of iskpts and nkptsp)
   call orbitals_communicators(iproc,nproc,Glr,orbs,comms)

   if(orbs%nspinor > 1) STOP 'BigDFT2Wannier does not work for nspinor > 1'
   if(orbs%nkpts > 1) stop 'BigDFT2Wannier does not work for nkpts > 1'

   ! Read integers in order to allocate tables related to projectors 
   call read_nnkp_int_alloc(iproc,seedname, n_kpts, n_proj, n_nnkpts, n_excb)
   call allocate_initial()
   call orbitals_descriptors(iproc,nproc,n_proj,n_proj,0,1,1,&
        1,(/ 0.0_dp,0.0_dp,0.0_dp /),(/0.0_dp/),orbsp,.false.) 
   if(residentity) n_virt = n_proj

   ! Set-up number of virtual states
   if(pre_check .and. n_virt_tot > 0) then
     nvirtu = abs(input%norbv)
     nvirtd = 0
   else if(residentity) then
     nvirtu = n_proj
     nvirtd = 0
   else
     nvirtu = n_virt
     nvirtd = 0
   end if
   if (input%nspin==2) nvirtd=nvirtu
   call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
      &   orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv,.false.)

   ! Read Wannier90 .nnkp file.
   ! The most important informations to be read are : 
   !  - ctr_proj : the coordinates of the center of the projections
   !  - l, mr, rvalue and zona = Z/a : the parameters used to build spherical harmonics
   !  - k_plus_b and G_vec : the parameters used to build the nearest neighbours k-points
   call read_nnkp(iproc,seedname, calc_only_A, real_latt, recip_latt, n_kpts, n_proj, n_nnkpts, &
      &   n_excb, kpts, ctr_proj, x_proj, z_proj, l, mr, rvalue, &
      &   zona, k_plus_b, G_vec, excb)

   ! Check that z_proj and x_proj are orthonormal and build y_proj
   do np = 1, n_proj  
      znorm = sqrt(z_proj(np,1)**2+z_proj(np,2)**2+z_proj(np,3)**2)
      xnorm = sqrt(x_proj(np,1)**2+x_proj(np,2)**2+x_proj(np,3)**2)
      ortho = z_proj(np,1)*x_proj(np,1) + z_proj(np,2)*x_proj(np,2) + z_proj(np,3)*x_proj(np,3)
      if(abs(znorm - 1.d0) > eps6 .or. abs(znorm - 1.d0) > eps6 .or. abs(ortho) > eps6) then
         if(iproc == 0) then
            call yaml_warning('Check orthonormality of z_proj and x_proj:')
            call yaml_comment('z norm: ',trim(yaml_toa(znorm,fmt='(e9.7)')))
            call yaml_comment('x norm: ',trim(yaml_toa(xnorm,fmt='(e9.7)')))
            call yaml_comment('x dot z: ',trim(yaml_toa(ortho,fmt='(e9.7)')))
            !write(*,'(a)') 'check orthonormality of z_proj and x_proj:'
            !write(*,'(a,e9.7)') 'z norm: ',znorm
            !write(*,'(a,e9.7)') 'x norm: ',xnorm
            !write(*,'(a,e9.7)') 'x dot z: ',ortho
            stop
         end if
      end if
      !Now we need to calculate the y direction
      y_proj(np,1) = x_proj(np,3)*z_proj(np,2) - x_proj(np,2)*z_proj(np,3)
      y_proj(np,2) = x_proj(np,1)*z_proj(np,3) - x_proj(np,3)*z_proj(np,1)
      y_proj(np,3) = x_proj(np,2)*z_proj(np,1) - x_proj(np,1)*z_proj(np,2)
      !print *,'np,norms,ortho,y',np,znorm,xnorm,ortho,y_proj(np,:) 
   end do

   !distribute the projectors on the processes (contained in orbsp: norb,norbp,isorb,...)
   call split_vectors_for_parallel(iproc,nproc,n_proj,orbsp)
   call orbitals_communicators(iproc,nproc,Glr,orbsp,commsp)

   ! Simplification of the notations
   nx=Glr%d%n1i
   ny=Glr%d%n2i
   nz=Glr%d%n3i
   n_at=atoms%astruct%nat
   call initialize_work_arrays_sumrho(Glr,w)

   ! Allocations for Amnk calculation
   npsidim2=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsp%norbp,sum(commsp%ncntt(0:nproc-1)))
   allocate(ylm(nx,ny,nz),stat=i_stat)
   call memocc(i_stat,ylm,'ylm',subname)
   allocate(func_r(nx,ny,nz),stat=i_stat)
   call memocc(i_stat,func_r,'func_r',subname)
   allocate(sph_har_etsf(nx*ny*nz),stat=i_stat)
   call memocc(i_stat,sph_har_etsf,'sph_har_etsf',subname)
   allocate(amnk_bands_sorted(n_virt),stat=i_stat)
   call memocc(i_stat,amnk_bands_sorted,'amnk_bands_sorted',subname)


   call timing(iproc,'Precondition  ','OF')

      if (pre_check .and. n_virt_tot > 0 .and. .not. residentity) then

         call split_vectors_for_parallel(iproc,nproc,n_virt_tot,orbsv)
         call orbitals_communicators(iproc,nproc,Glr,orbsv,commsv) 

         call timing(iproc,'CrtProjectors ','ON')

         ! Read wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
         npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsv%norbp*orbsv%nspinor,sum(commsv%ncntt(0:nproc-1)))
         allocate(psi_etsfv(npsidim),stat=i_stat)
         call memocc(i_stat,psi_etsfv,'psi_etsfv',subname)
         if(associated(orbsv%eval)) nullify(orbsv%eval)
         allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
         call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)

         filename= trim(input%dir_output) // 'virtuals'
         call readmywaves(iproc,filename,iformat,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3, &
              & input%hx,input%hy,input%hz,atoms,rxyz_old,atoms%astruct%rxyz,  & 
              Glr%wfd,psi_etsfv)
         i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
         deallocate(orbsv%eval,stat=i_stat)
         nullify(orbsv%eval)
         call memocc(i_stat,i_all,'orbsv%eval',subname)

         if(nproc > 1) then
            allocate(pwork(npsidim),stat=i_stat)
            call memocc(i_stat,pwork,'pwork',subname)
            call transpose_v(iproc,nproc,orbsv,Glr%wfd,commsv,psi_etsfv(1),work=pwork)
            i_all = -product(shape(pwork))*kind(pwork)
            deallocate(pwork,stat=i_stat)
            call memocc(i_stat,i_all,'pwork',subname)
         end if

         ! - b1, b2 and b3 are the norm of the lattice parameters.
         b1=atoms%astruct%cell_dim(1)
         b2=atoms%astruct%cell_dim(2)
         b3=atoms%astruct%cell_dim(3)
         ! - Allocations
         allocate(amnk(orbsv%norb,orbsp%norb),stat=i_stat)
         call memocc(i_stat,amnk,'amnk',subname)
         allocate(amnk_guess(orbsv%norb),stat=i_stat)
         call memocc(i_stat,amnk_guess,'amnk_guess',subname)
         allocate(sph_daub(npsidim2), stat=i_stat)
         call memocc(i_stat,sph_daub,'sph_daub',subname)
         call to_zero(npsidim2,sph_daub(1))

         ! Begining of the algorithm to compute the scalar product in order to find the best unoccupied orbitals
         ! to use to compute the actual Amnk matrix :
         if (iproc==0) then
            call yaml_comment('',hfill='=')
            call yaml_comment('Calculating amnk=<virt|sph_har> in pre-check mode')
            call yaml_comment('',hfill='=')
            !write(*,*) '!==================================!'
            !write(*,*) '! Calculating amnk=<virt|sph_har>  !'
            !write(*,*) '!       in pre-check mode :        !'
            !write(*,*) '!==================================!'
            !write(*,'(A12,4x,A15)') 'Virtual band', 'amnk_guess(nb)='
         end if

         !calculate buffer shifts
         perx=(Glr%geocode /= 'F')
         pery=(Glr%geocode == 'P')
         perz=(Glr%geocode /= 'F')
         call ext_buffers(perx,nbl1,nbr1)
         call ext_buffers(pery,nbl2,nbr2)
         call ext_buffers(perz,nbl3,nbr3)

         ! Calculation of the spherical harmonics in parallel.
         ! It is done in the real space and then converted in the Daubechies representation.
         pshft = 0
         do npp=1, orbsp%norbp
            np = npp + orbsp%isorb
            ! Convolution buffer : n1i=2*n1+31 -> explains the '13*input%hx*0.5' term
            r0x=ctr_proj(np,1)*b1+nbl1*input%hx*0.5
            r0y=ctr_proj(np,2)*b2+nbl2*input%hy*0.5
            r0z=ctr_proj(np,3)*b3+nbl3*input%hz*0.5
            do k=1+nbl3,nz-nbr3
               zz=(k-1)*input%hz*0.5-r0z
               do j=1+nbl2,ny-nbr2
                  yy=(j-1)*input%hy*0.5-r0y
                  do i=1+nbl1,nx-nbl1
                     ind=(k-1)*ny*nx+(j-1)*nx+i
                     xx=(i-1)*input%hx*0.5-r0x
                     call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                        &   xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
                     call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                        &   xx, yy, zz, n_proj, func_r)
                     ! The 'sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)' term is here to normalize spherical harmonics
                     sph_har_etsf(ind)=func_r(i,j,k)*ylm(i,j,k)*sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)
                  end do
               end do
            end do
            xnorm = ddot(nz*ny*nx,sph_har_etsf(1),1,sph_har_etsf(1),1)
            call dscal(nz*ny*nx,1.0_dp/sqrt(xnorm),sph_har_etsf(1),1)
            if(w_sph .or. w_ang .or. w_rad) then
               call write_functions(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, Glr, &
               &    0.5_dp*input%hx, 0.5_dp*input%hy, 0.5_dp*input%hz, atoms, atoms%astruct%rxyz, sph_har_etsf, func_r, ylm)
            end if
            call isf_to_daub(Glr,w,sph_har_etsf(1),sph_daub(1+pshft))
            pshft=pshft + Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f 
         end do

         call deallocate_projectors()

         call timing(iproc,'CrtProjectors ','OF')

         ! Tranposition of the distribution of the spherical harmonics: orbitals -> components.
         allocate(pwork(npsidim2),stat=i_stat)
         call memocc(i_stat,pwork,'pwork',subname)
         call transpose_v(iproc,nproc,orbsp,Glr%wfd,commsp,sph_daub(1),work=pwork)
         i_all = -product(shape(pwork))*kind(pwork)
         deallocate(pwork,stat=i_stat)
         call memocc(i_stat,i_all,'pwork',subname)
         call timing(iproc,'ApplyProj     ','ON')

         ! calculate the inverse overlap of the projector to build proper amnk_tot
         allocate(overlap_proj(orbsp%norb,orbsp%norb))
         nvctrp=commsp%nvctr_par(iproc,1)
         call gemm('T','N',orbsp%norb,orbsp%norb,nvctrp,1.0_wp,sph_daub(1),max(1,nvctrp),&
            &   sph_daub(1),max(1,nvctrp),0.0_wp,overlap_proj(1,1),orbsp%norb)
         if(nproc > 1) then
            call mpiallred(overlap_proj(1,1),orbsp%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
         end if
         !print *,'overlap_proj',overlap_proj
         allocate(ipiv(orbsp%norb),stat=i_stat)
         call memocc(i_stat,ipiv,'ipiv',subname)
         call dgetrf( orbsp%norb, orbsp%norb, overlap_proj, orbsp%norb, ipiv, info )
         allocate(pwork(1),stat=i_stat)
         call memocc(i_stat,pwork,'pwork',subname)
         call dgetri( orbsp%norb, overlap_proj, orbsp%norb, ipiv, pwork, -1, info )
         lwork = int(pwork(1))
         i_all = -product(shape(pwork))*kind(pwork)
         deallocate(pwork,stat=i_stat)
         call memocc(i_stat,i_all,'pwork',subname)
         allocate(pwork(lwork),stat=i_stat)
         call memocc(i_stat,pwork,'pwork',subname)
         call dgetri( orbsp%norb, overlap_proj, orbsp%norb, ipiv, pwork, lwork, info )
         i_all = -product(shape(pwork))*kind(pwork)
         deallocate(pwork,stat=i_stat)
         call memocc(i_stat,i_all,'pwork',subname)
         i_all = -product(shape(ipiv))*kind(ipiv)
         deallocate(ipiv,stat=i_stat)
         call memocc(i_stat,i_all,'ipiv',subname)
         !print *,'inverse overlap_proj',overlap_proj

         ! Scalar product of amnk=<sph_daub|psi> in parallel.
         call razero(orbsp%norb*orbsv%norb,amnk)
         nvctrp=commsv%nvctr_par(iproc,1)
         call gemm('T','N',orbsv%norb,orbsp%norb,nvctrp,1.0_wp,psi_etsfv(1),max(1,nvctrp),&
            &   sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),orbsv%norb)

         ! Construction of the whole Amnk_guess matrix.
         if(nproc > 1) then
            call mpiallred(amnk(1,1),orbsv%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
         end if

         ! For each unoccupied orbitals, check how they project on spherical harmonics.
         ! The greater amnk_guess(nb) is, the more they project on spherical harmonics.
         do nb=1,orbsv%norb
            amnk_guess(nb)=0.0
            do np=1,orbsp%norb
               do j=1,orbsp%norb
                  amnk_guess(nb)= amnk_guess(nb) + amnk(nb,np)*amnk(nb,j)*overlap_proj(np,j)
               end do
            end do
            if (iproc==0) then
               call yaml_map('Virtual band',nb)
               call yaml_map('amnk_guess(nb)',sqrt(amnk_guess(nb)))
            end if
            !if (iproc==0) write(*,'(I4,11x,F12.6)') nb, sqrt(amnk_guess(nb))
         end do

         ! Choice of the unoccupied orbitals to calculate the Amnk matrix
         if (iproc==0) then
            write(*,*) 
            write(*,'(1a)') 'These are the virtual bands to use to construct the actual Amn and Mmn matrices :'
            write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_guess)='
         end if
         allocate(amnk_guess_sorted(n_virt),stat=i_stat)
         call memocc(i_stat,amnk_guess_sorted,'amnk_guess_sorted',subname)
         do nb=1,n_virt
            amnk_guess_sorted(nb)=maxval(amnk_guess,1)
            amnk_bands_sorted(nb)=maxloc(amnk_guess,1)
            amnk_guess(amnk_bands_sorted(nb))=0.d0
            if (iproc==0) write(*,'(I4,3x,F12.6)') amnk_bands_sorted(nb), sqrt(amnk_guess_sorted(nb))
         end do

         ! End of the pre-check mode
         i_all = -product(shape(psi_etsfv))*kind(psi_etsfv)
         deallocate(psi_etsfv,stat=i_stat)
         call memocc(i_stat,i_all,'psi_etsfv',subname)
         i_all = -product(shape(amnk))*kind(amnk)
         deallocate(amnk,stat=i_stat)
         call memocc(i_stat,i_all,'amnk',subname)
         i_all = -product(shape(amnk_guess_sorted))*kind(amnk_guess_sorted)
         deallocate(amnk_guess_sorted,stat=i_stat)
         call memocc(i_stat,i_all,'amnk_guess_sorted',subname)
         i_all = -product(shape(amnk_guess))*kind(amnk_guess)
         deallocate(amnk_guess,stat=i_stat)
         call memocc(i_stat,i_all,'amnk_guess',subname)
         call deallocate_comms(commsv,subname)

         if (iproc==0) then
            write(*,*) '!==================================!'
            write(*,*) '! Calculating amnk=<virt|sph_har>  !'
            write(*,*) '!     in pre-check mode done       !'
            write(*,*) '!==================================!'
            write(*,*)
            write(*,*)
         end if

         ! Rewrite the input.inter file to add the chosen unoccupied states.
         if (iproc==0) call write_inter(n_virt, nx, ny, nz, amnk_bands_sorted)

         call timing(iproc,'ApplyProj     ','OF')

      end if

      call timing(iproc,'CrtDescriptors','ON')

      ! Define which unoccupied orbitals have to be read.
      ! Here, virt_list nominates the virtual orbitals chosen in pre-check mode.
      allocate(virt_list(n_virt),stat=i_stat)
      call memocc(i_stat,virt_list,'virt_list',subname)

      if (n_virt .ne. 0) then
         if (pre_check) then 
            do i=1,n_virt
               virt_list(i)=amnk_bands_sorted(i)!+n_occ
            end do
         else if(residentity) then
            virt_list = 0
         else
            call read_inter_list(iproc, n_virt, virt_list)
         end if
      end if

      !Setup the description of the new subspace (they are similar to orbitals)
      call orbitals_descriptors(iproc,nproc,orbs%norb,orbs%norbu,orbs%norbd,orbs%nspin,orbs%nspinor,&
           orbs%nkpts,orbs%kpts,orbs%kwgts,orbsb,.false.)

      ! Initialise the arrays n_bands_par, isband_par
      call split_vectors_for_parallel(iproc,nproc,n_virt,orbsv)
      call split_vectors_for_parallel(iproc,nproc,n_occ+n_virt,orbsb)
      call orbitals_communicators(iproc,nproc,Glr,orbsb,commsb)


      call timing(iproc,'CrtDescriptors','OF')
      call timing(iproc,'CrtProjectors ','ON')

      npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsb%norbp*orbsb%nspinor,sum(commsb%ncntt(0:nproc-1)))
      allocate(psi_etsf(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,max(orbsb%norbp*orbsb%nspinor,1)),stat=i_stat)
      call memocc(i_stat,psi_etsf,'psi_etsf',subname)

      ! For the occupied orbitals, need to modifify norbp,isorb to match the total distributed scheme
      orbs%norbp = n_occ - orbsb%isorb
      if (orbsb%isorb + orbsb%norbp < n_occ ) orbs%norbp = orbsb%norbp
      if(orbsb%isorb > n_occ) orbs%norbp = 0
      orbs%isorb = orbsb%isorb
      if(associated(orbs%iokpt)) then
         i_all = -product(shape(orbs%iokpt))*kind(orbs%iokpt)
         deallocate(orbs%iokpt,stat=i_stat)
         call memocc(i_stat,i_all,'orbs%iokpt',subname)
      end if
      allocate(orbs%iokpt(orbs%norbp),stat=i_stat)
      call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)
      orbs%iokpt=1

      if(associated(orbs%eval)) nullify(orbs%eval)
      allocate(orbs%eval(orbs%norb*orbs%nkpts), stat=i_stat)
      call memocc(i_stat,orbs%eval,'orbs%eval',subname)
      call razero(orbs%norb*orbs%nkpts,orbs%eval)
      if(orbs%norbp > 0) then
            filename=trim(input%dir_output) // 'wavefunction'
            call readmywaves(iproc,filename,iformat,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
            input%hx,input%hy,input%hz,atoms,rxyz_old,atoms%astruct%rxyz,Glr%wfd,psi_etsf(1,1))
      end if
      ! For bin files, the eigenvalues are distributed, so reduce them
      if((filetype == 'bin' .or. filetype == 'BIN') .and. nproc > 0) then
        call mpiallred(orbs%eval(1),orbs%norb*orbs%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
      ! Write the eigenvalues into a file to output the hamiltonian matrix elements in Wannier functions
      if(iproc==0) then     
        open(15, file=trim(seedname)//'.eig', status='unknown')
        do nb = 1, orbs%norb            !TO DO: ADD KPTS
           write(15,'(I4,2x,I4,2x,E17.9)') nb, 1, orbs%eval(nb)
        end do
      end if
      i_all = -product(shape(orbs%eval))*kind(orbs%eval)
      deallocate(orbs%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbs%eval',subname)

      ! For the non-occupied orbitals, need to change norbp,isorb
      orbsv%norbp = orbsb%isorb + orbsb%norbp - n_occ
      if (orbsb%isorb + orbsb%norbp < n_occ ) orbsv%norbp = 0
      if (orbsb%isorb > n_occ) orbsv%norbp = orbsb%norbp
      orbsv%isorb = 0
      if(orbsb%isorb >= n_occ) orbsv%isorb = orbsb%isorb - n_occ
      if(associated(orbsv%iokpt)) then
         i_all = -product(shape(orbsv%iokpt))*kind(orbsv%iokpt)
         deallocate(orbsv%iokpt,stat=i_stat)
         call memocc(i_stat,i_all,'orbsv%iokpt',subname)
      end if
      allocate(orbsv%iokpt(orbsv%norbp),stat=i_stat)
      call memocc(i_stat,orbsv%iokpt,'orbsv%iokpt',subname)
      orbsv%iokpt=1
      !orbsv%spinsgn= 1.0

         ! read unoccupied wavefunctions
      if(associated(orbsv%eval)) nullify(orbsv%eval)
      allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
      call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)
      call razero(orbsv%norb*orbsv%nkpts,orbsv%eval)
      if(orbsv%norbp > 0 .and. .not. residentity) then
         filename=trim(input%dir_output) // 'virtuals'
            call readmywaves(iproc,filename,iformat,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
            input%hx,input%hy,input%hz,atoms,rxyz_old,atoms%astruct%rxyz, Glr%wfd,psi_etsf(1,1+orbs%norbp),virt_list)
      end if
      ! For bin files, the eigenvalues are distributed, so reduce them
      if(residentity)then
         orbsv%eval = 99.0_dp  !What to put for the energy?
      else if((filetype == 'bin' .or. filetype == 'BIN') .and.  nproc > 0 .and. orbsv%norb>0) then
         call mpiallred(orbsv%eval(1),orbsv%norb*orbsv%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
      ! Write the eigenvalues into a file to output the hamiltonian matrix elements in Wannier functions
      if(iproc==0) then
        do nb = 1, orbsv%norb            !TO DO: ADD KPTS
           write(15,'(I4,2x,I4,2x,E17.9)') nb+n_occ, 1, orbsv%eval(nb)
        end do
        close(15)
      end if
      i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
      deallocate(orbsv%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbsv%eval',subname)
      i_all = -product(shape(rxyz_old))*kind(rxyz_old)
      deallocate(rxyz_old,stat=i_stat)
      call memocc(i_stat,i_all,'rxyz_old',subname)

      ! Algorithm to compute the scalar product of the input guess:
      ! The term 'sqrt(bx(1)*by(2)*bz(3))' is there to normalize spherical harmonics.
      ! Wavefunctions calculated by BigDFT already are normalized.
      if (iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '!  Calculating amnk=<psi|sph_har>  !'
         write(*,*) '!==================================!'
      end if

      ! - b1, b2 and b3 are the norm of the lattice parameters.
      b1=atoms%astruct%cell_dim(1)
      b2=atoms%astruct%cell_dim(2)
      b3=atoms%astruct%cell_dim(3)
      ! - Allocations
      allocate(amnk(orbsb%norb,orbsp%norb),stat=i_stat)
      call memocc(i_stat,amnk,'amnk',subname)
      call to_zero(orbsb%norb*orbsp%norb,amnk(1,1))
      allocate(amnk_tot(orbsb%norb),stat=i_stat)
      call memocc(i_stat,amnk_tot,'amnk_tot',subname)

      call timing(iproc,'CrtProjectors ','OF')

      ! Calculation of the spherical harmonics in parallel (only if not done in precheck).
      ! It is done in the real space and then converted in the Daubechies representation.
      if(.not. pre_check .or. n_virt_tot == 0 .or. residentity) then
         call timing(iproc,'CrtProjectors ','ON')
         allocate(sph_daub(npsidim2), stat=i_stat)
         call memocc(i_stat,sph_daub,'sph_daub',subname)
         if(npsidim2 > 0) call to_zero(npsidim2,sph_daub(1))
         !calculate buffer shifts
         perx=(Glr%geocode /= 'F')
         pery=(Glr%geocode == 'P')
         perz=(Glr%geocode /= 'F')
         call ext_buffers(perx,nbl1,nbr1)
         call ext_buffers(pery,nbl2,nbr2)
         call ext_buffers(perz,nbl3,nbr3)

         pshft = 0
         do npp=1, orbsp%norbp
            np=npp+orbsp%isorb
            ! Convolution buffer : n1i=2*n1+31 -> explains the '13*input%hx*0.5' term
            r0x=ctr_proj(np,1)*b1+nbl1*input%hx*0.5
            r0y=ctr_proj(np,2)*b2+nbl2*input%hy*0.5
            r0z=ctr_proj(np,3)*b3+nbl3*input%hz*0.5
            call to_zero(nz*ny*nx,sph_har_etsf(1))
            do k=1+nbl3,nz-nbr3
               zz=(k-1)*input%hz*0.5-r0z
               do j=1+nbl2,ny-nbr2
                  yy=(j-1)*input%hy*0.5-r0y
                  do i=1+nbl1,nx-nbr1
                     ind=(k-1)*ny*nx+(j-1)*nx+i
                     xx=(i-1)*input%hx*0.5-r0x
                     call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
                        &   xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
                     call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
                        &   xx, yy, zz, n_proj, func_r)
                     ! The 'sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)' term is here to normalize spherical harmonics
                     sph_har_etsf(ind)=func_r(i,j,k)*ylm(i,j,k)*sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)
                  end do
               end do
            end do
            xnorm = ddot(nz*ny*nx,sph_har_etsf(1),1,sph_har_etsf(1),1)
            call dscal(nz*ny*nx,1.0_dp/sqrt(xnorm),sph_har_etsf(1),1)
            if(w_sph .or. w_ang .or. w_rad) then
               call write_functions(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, Glr, &
               &    0.5_dp*input%hx, 0.5_dp*input%hy, 0.5_dp*input%hz, atoms, atoms%astruct%rxyz, sph_har_etsf, func_r, ylm)
            end if
            !print *,'before isf_to_daub',sqrt(ddot(nz*ny*nx,sph_har_etsf(1),1,sph_har_etsf(1),1))
            call isf_to_daub(Glr,w,sph_har_etsf(1),sph_daub(1+pshft))
            !print *,'after isf_to_daub',sqrt(ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,sph_daub(1+pshft),1,sph_daub(1+pshft),1))
            pshft=pshft + Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f 
         end do

         call timing(iproc,'CrtProjectors ','OF')
         ! Tranposition of distribution : orbitals -> components.
         if(nproc>1) then
            allocate(pwork(npsidim2),stat=i_stat)
            call memocc(i_stat,pwork,'pwork',subname)
            call transpose_v(iproc,nproc,orbsp,Glr%wfd,commsp,sph_daub,work=pwork)
            i_all = -product(shape(pwork))*kind(pwork)
            deallocate(pwork,stat=i_stat)
            call memocc(i_stat,i_all,'pwork',subname)
         end if
         call deallocate_projectors()

         ! calculate the inverse overlap of the projector to build proper amnk_tot
         allocate(overlap_proj(orbsp%norb,orbsp%norb))
         nvctrp=commsp%nvctr_par(iproc,1)
         call gemm('T','N',orbsp%norb,orbsp%norb,nvctrp,1.0_wp,sph_daub(1),max(1,nvctrp),&
            &   sph_daub(1),max(1,nvctrp),0.0_wp,overlap_proj(1,1),orbsp%norb)
         if(nproc > 1) then
            call mpiallred(overlap_proj(1,1),orbsp%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
         end if
         !print *,'overlap_proj',overlap_proj
         allocate(ipiv(orbsp%norb),stat=i_stat)
         call memocc(i_stat,ipiv,'ipiv',subname)
         call dgetrf( orbsp%norb, orbsp%norb, overlap_proj, orbsp%norb, ipiv, info )
         allocate(pwork(1),stat=i_stat)
         call memocc(i_stat,pwork,'pwork',subname)
         call dgetri( orbsp%norb, overlap_proj, orbsp%norb, ipiv, pwork, -1, info )
         lwork = int(pwork(1))
         i_all = -product(shape(pwork))*kind(pwork)
         deallocate(pwork,stat=i_stat)
         call memocc(i_stat,i_all,'pwork',subname)
         allocate(pwork(lwork),stat=i_stat)
         call memocc(i_stat,pwork,'pwork',subname)
         call dgetri( orbsp%norb, overlap_proj, orbsp%norb, ipiv, pwork, lwork, info )
         i_all = -product(shape(pwork))*kind(pwork)
         deallocate(pwork,stat=i_stat)
         call memocc(i_stat,i_all,'pwork',subname)
         i_all = -product(shape(ipiv))*kind(ipiv)
         deallocate(ipiv,stat=i_stat)
         call memocc(i_stat,i_all,'ipiv',subname)

         !print *,'inverse overlap_proj',overlap_proj

      end if

      if(.not. residentity)then
         ! Tranposition of the distribution of the BigDFT wavefunctions : orbitals -> components.
         npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsb%norbp*orbsb%nspinor,sum(commsb%ncntt(0:nproc-1)))
         allocate(psi_etsf2(npsidim),stat=i_stat) !!doing this because psi_etsfv does not incorporate enough space for transpose
         call memocc(i_stat,psi_etsf2,'psi_etsf2',subname)
         call razero(npsidim,psi_etsf2)
         if(nproc > 1) then
            allocate(pwork(npsidim),stat=i_stat)
            call memocc(i_stat,pwork,'pwork',subname)
            call transpose_v(iproc,nproc,orbsb,Glr%wfd,commsb,psi_etsf(1,1),work=pwork,outadd=psi_etsf2(1))
            i_all = -product(shape(pwork))*kind(pwork)
            deallocate(pwork,stat=i_stat)
            call memocc(i_stat,i_all,'pwork',subname)
         else
            call vcopy(orbsb%norb*orbsb%nspinor*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),psi_etsf(1,1),1,psi_etsf2(1),1)
         end if
      else
         ! Tranposition of the distribution of the BigDFT occupied wavefunctions : orbitals -> components.
         npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsb%norbp*orbsb%nspinor,sum(commsb%ncntt(0:nproc-1)))
         allocate(psi_etsf2(npsidim),stat=i_stat) 
         call memocc(i_stat,psi_etsf2,'psi_etsf2',subname)
         call razero(npsidim,psi_etsf2)
         if(nproc > 1) then
            allocate(pwork(npsidim),stat=i_stat)
            call memocc(i_stat,pwork,'pwork',subname)
            call transpose_v(iproc,nproc,orbsb,Glr%wfd,commsb,psi_etsf(1,1),work=pwork,outadd=psi_etsf2(1))
            i_all = -product(shape(pwork))*kind(pwork)
            deallocate(pwork,stat=i_stat)
            call memocc(i_stat,i_all,'pwork',subname)
         else
            call vcopy(orbsb%norb*orbs%nspinor*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),psi_etsf(1,1),1,psi_etsf2(1),1)
         end if

         ! Scalar product of amnk=<sph_daub|occ> in parallel.
         nvctrp=commsb%nvctr_par(iproc,1)
         call gemm('T','N',orbsb%norb,orbsp%norb,nvctrp,1.0_wp,psi_etsf2(1),max(1,nvctrp),&
            &   sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),orbsb%norb)

         ! Construction of the occupied Amnk submatrix.
         if(nproc > 1) then
            call mpiallred(amnk(1,1),orbsb%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
         end if

         ! Now build the new states corresponding to: sph_daub - sum{amnk occ} and place them at the virtual states
         ! We are working in transpose space
         pshft = 0
         shft  = nvctrp*orbs%norb
         do np=1, orbsp%norb
            !np=npp+orbsp%isorb 
            call vcopy(nvctrp,sph_daub(1+pshft),1,psi_etsf2(1+shft),1)
            wshft = 0
            do nb = 1,orbs%norb*orbs%nspinor    !Should be on all occupied bands
               do j=1,nvctrp
                  psi_etsf2(j+shft) = psi_etsf2(j+shft) - amnk(nb,np)*psi_etsf2(j+wshft) 
               end do
               wshft = wshft + nvctrp
            end do
            pshft = pshft +  nvctrp
            shft = shft + nvctrp
         end do

         ! Now untranspose to make the psi
         if(nproc > 1) then
            allocate(pwork(npsidim),stat=i_stat)
            call memocc(i_stat,pwork,'pwork',subname)
            call untranspose_v(iproc,nproc,orbsb,Glr%wfd,commsb,psi_etsf2,pwork,outadd=psi_etsf(1,1))
            i_all = -product(shape(pwork))*kind(pwork)
            deallocate(pwork,stat=i_stat)
            call memocc(i_stat,i_all,'pwork',subname)
         else
            call vcopy(orbsb%norb*orbsb%nspinor*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),psi_etsf2(1),1,psi_etsf(1,1),1)
         end if
         ! Should write the symmetrized projectors to file
         if(write_resid .and. orbsv%norbp > 0)then
            allocate(orbsv%eval(orbsv%norb))
            orbsv%eval = 99.0_dp
            call writemywaves(iproc,trim(input%dir_output) // "virtuals",iformat,orbsv,Glr%d%n1,Glr%d%n2,&
              Glr%d%n3,input%hx,input%hy,input%hz,atoms,atoms%astruct%rxyz,Glr%wfd,psi_etsf(1,1+orbs%norbp))
            deallocate(orbsv%eval)
         end if
         
      end if

      ! Scalar product of amnk=<sph_daub|psi> in parallel.
      nvctrp=commsb%nvctr_par(iproc,1)
      call gemm('T','N',orbsb%norb,orbsp%norb,nvctrp,1.0_wp,psi_etsf2(1),max(1,nvctrp),&
         &   sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),orbsb%norb)

      ! Construction of the whole Amnk matrix.
      if(nproc > 0) then
         call mpiallred(amnk(1,1),orbsb%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if

      call to_zero(orbsb%norb,amnk_tot(1))
      if (iproc==0) then
         ! Check normalisation (the amnk_tot value must tend to 1).
         write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_tot)='
         do nb=1,orbsb%norb
            do np=1,orbsp%norb
               do j=1,orbsp%norb
                  amnk_tot(nb)= amnk_tot(nb) + amnk(nb,np)*amnk(nb,j)*overlap_proj(np,j)
               end do
            end do
            write(*,'(I4,3x,F12.6)') nb, sqrt(amnk_tot(nb))
         end do
         write(*,*) '!==================================!'
         write(*,*) '!  Calculating amnk=<psi|sph_har>  !'
         write(*,*) '!               done               !'
         write(*,*) '!==================================!'
         write(*,*)
         write(*,*)
      end if
 
     ! Write the .amn file
      if (iproc == 0) call write_amn(seedname, orbsb%norb, n_kpts, orbsp%norb, amnk)

      !TEST NEW SCHEME
      !!allocate(pwork(npsidim2),stat=i_stat)
      !!call memocc(i_stat,pwork,'pwork',subname)
      !!call untranspose_v(iproc,nproc,orbsp,Glr%wfd,commsp,sph_daub,work=pwork)
      !!allocate(pvirt(npsidim2))
      !!pshft = 0
      !!do npp=1, orbsp%norbp
      !!   np=npp+orbsp%isorb 
      !!   call vcopy(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,sph_daub(1+pshft),1,pvirt(1+pshft),1)
      !!   do i = 1,orbsb%norbp*orbsb%nspinor    !Should be on all bands
      !!      nb = i + orbsb%isorb
      !!      !print *,'amnk(nb,np)',nb,np,amnk(nb,np),ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,psi_etsf(1,i),1,psi_etsf(1,i),1)
      !!      do j=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      !!         pvirt(j+pshft) = pvirt(j+pshft) - amnk(nb,np)*psi_etsf(j,i) 
      !!      end do
      !!   end do
      !!   print *,'Norm of proj',np,ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,sph_daub(1+pshft),1,sph_daub(1+pshft),1)
      !!   print *,'Norm of Symproj',np,ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,pvirt(1+pshft),1,pvirt(1+pshft),1)
      !!   pshft = pshft +  Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      !!end do
      !!!call transpose_v(iproc,nproc,orbsp,Glr%wfd,commsp,pvirt,work=pwork)
      !!i_all = -product(shape(pwork))*kind(pwork)
      !!deallocate(pwork,stat=i_stat)
      !!call memocc(i_stat,i_all,'pwork',subname)
      !!allocate(orbsp%eval(orbsp%norb))
      !!orbsp%eval = 0.5_dp
      !!call writemywaves(iproc,trim(input%dir_output) // "virtuals",2,orbsp,Glr%d%n1,Glr%d%n2,&
      !!     Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz,Glr%wfd,pvirt)
      !!print *,'Generated the symproj'
      !!call mpi_finalize(ierr)
      !!stop
      !END TEST

      call deallocate_amnk_calculation()

      call timing(iproc,'ApplyProj     ','OF')
      call timing(iproc,'Input_comput  ','ON')

      if (iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '!   Calculating mmnk=<psi|psi> :   !'
         write(*,*) '!==================================!'
         write(*,*) 'The values of sqrt(mmnk_tot) check the normalization in each band.'
         write(*,*) 'They must tend to their lower limit value, which is equal to 1 :'
      end if

      call mmnk_calculation_allocation()

      !calculate buffer shifts
      perx=(Glr%geocode /= 'F')
      pery=(Glr%geocode == 'P')
      perz=(Glr%geocode /= 'F')
      call ext_buffers(perx,nbl1,nbr1)
      call ext_buffers(pery,nbl2,nbr2)
      call ext_buffers(perz,nbl3,nbr3)

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
         pshft = 0
         call to_zero(nx*ny*nz,psir_re(1))
         call to_zero(nx*ny*nz,psir_im(1))
         call to_zero(npsidim,psi_daub_re(1))
         call to_zero(npsidim,psi_daub_im(1))
         do nb1=1,orbsb%norbp
            call daub_to_isf(Glr,w,psi_etsf(1,nb1),psir)
            do k=1,nz
               zz=(k-nbl3)*input%hz*0.5
               do j=1,ny
                  yy=(j-nbl2)*input%hy*0.5
                  do i=1,nx
                     xx=(i-nbl1)*input%hx*0.5
                     ind=(k-1)*ny*nx+(j-1)*nx+i
                     psir_re(ind)= psir(ind) * cos( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                     psir_im(ind)=-psir(ind) * sin( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
                  end do
               end do
            end do
            call isf_to_daub(Glr,w,psir_re(1),psi_daub_re(1+pshft))
            call isf_to_daub(Glr,w,psir_im(1),psi_daub_im(1+pshft))
            !pshft = pshft + max(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,commsb%ncntt(iproc)/orbsb%norbp)
            pshft = pshft + Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
         end do

         ! Tranposition of distribution : orbitals -> components
         if(nproc>0 .or. orbsb%nspinor /=1) then
            call timing(iproc,'Input_comput  ','OF')
            allocate(pwork(npsidim),stat=i_stat)
            call memocc(i_stat,pwork,'pwork',subname)
            call transpose_v(iproc,nproc,orbsb,Glr%wfd,commsb,psi_daub_re,work=pwork)
            call transpose_v(iproc,nproc,orbsb,Glr%wfd,commsb,psi_daub_im,work=pwork)
            i_all = -product(shape(pwork))*kind(pwork)
            deallocate(pwork,stat=i_stat)
            call memocc(i_stat,i_all,'pwork',subname)
            call timing(iproc,'Input_comput  ','ON')
         end if

         ! Scalar product to compute the overlap matrix
         nvctrp=commsb%nvctr_par(iproc,1)
         allocate(mmnk_v_re(orbsb%norb*orbsb%norb),stat=i_stat)
         call memocc(i_stat,mmnk_v_re,'mmnk_v_re',subname)
         allocate(mmnk_v_im(orbsb%norb*orbsb%norb),stat=i_stat)
         call memocc(i_stat,mmnk_v_im,'mmnk_v_im',subname)

         call gemm('T','N',orbsb%norb,orbsb%norb,nvctrp,1.0_wp,psi_daub_re(1),max(1,nvctrp),&
            &   psi_etsf2(1),max(1,nvctrp),0.0_wp,mmnk_v_re(1),orbsb%norb)
         call gemm('T','N',orbsb%norb,orbsb%norb,nvctrp,1.0_wp,psi_daub_im(1),max(1,nvctrp),&
            &   psi_etsf2(1),max(1,nvctrp),0.0_wp,mmnk_v_im(1),orbsb%norb)

         ! Reduce the overlap matrix between all the processors
         if (nproc > 1) then
            call mpiallred(mmnk_v_re(1),orbsb%norb*orbsb%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
            call mpiallred(mmnk_v_im(1),orbsb%norb*orbsb%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
         end if

         ! Reshape the overlap matrix elements into a more manageable disposition
         mmnk_re(:,:,inn)=reshape(mmnk_v_re,(/orbsb%norb,orbsb%norb/))
         mmnk_im(:,:,inn)=reshape(mmnk_v_im,(/orbsb%norb,orbsb%norb/))

         i_all = -product(shape(mmnk_v_re))*kind(mmnk_v_re)
         deallocate(mmnk_v_re,stat=i_stat)
         call memocc(i_stat,i_all,'mmnk_v_re',subname)
         i_all = -product(shape(mmnk_v_im))*kind(mmnk_v_im)
         deallocate(mmnk_v_im,stat=i_stat)
         call memocc(i_stat,i_all,'mmnk_v_im',subname)

         ! Check the normalisation
         do nb1=1, orbsb%norb
            mmnk_tot(nb1,inn)=sum(mmnk_re(nb1,:,inn)**2+mmnk_im(nb1,:,inn)**2)
            if (iproc==0) write(*,'(I4,3x,F12.6)') nb1, sqrt(mmnk_tot(nb1,inn))
         end do
      end do

      if (iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '! Calculating mmnk=<psi|psi> done  !'
         write(*,*) '!==================================!'
         write(*,*)
         write(*,*)
      end if

      ! Write the .mmn file
      if (iproc==0) call write_mmn(seedname, orbsb%norb, n_kpts, n_nnkpts, k_plus_b, G_vec, mmnk_re, mmnk_im)

      ! Write UNK file
      do nk=1, n_kpts
         s=1 ! s is the spin, set by default to 1
         if (w_unk) then 
            if (iproc==0) call write_unk_bin(Glr,orbs,orbsv,orbsb,input,atoms, &
               &   atoms%astruct%rxyz,n_occ,n_virt,virt_list,nx,ny,nz,nk,s,iformat)
         end if
      end do
   
      call final_deallocations()
      call timing(iproc,'Input_comput  ','OF')


call timing(bigdft_mpi%mpi_comm,'             ','RE')

call cpu_time(tcpu1)
call system_clock(ncount1,ncount_rate,ncount_max)
tel=dble(ncount1-ncount0)/dble(ncount_rate)
if (iproc == 0) &
   &   write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0 


!!$!finalize memory counting
!!$call memocc(0,0,'count','stop')


call bigdft_finalize(ierr)

!!$! Barrier suggested by support for titane.ccc.cea.fr, before finalise.
!!$call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!$
!!$call MPI_FINALIZE(ierr)
call f_lib_finalize()
contains

subroutine allocate_initial()

  allocate(kpts(n_kpts,3),stat=i_stat)
  call memocc(i_stat,kpts,'kpts',subname)
  allocate(ctr_proj(n_proj,3),stat=i_stat)
  call memocc(i_stat,ctr_proj,'ctr_proj',subname)
  allocate(x_proj(n_proj,3),stat=i_stat)
  call memocc(i_stat,x_proj,'x_proj',subname)
  allocate(y_proj(n_proj,3),stat=i_stat)
  call memocc(i_stat,y_proj,'y_proj',subname)
  allocate(z_proj(n_proj,3),stat=i_stat)
  call memocc(i_stat,z_proj,'z_proj',subname)
  allocate(l(n_proj),stat=i_stat)
  call memocc(i_stat,l,'l',subname)
  allocate(mr(n_proj),stat=i_stat)
  call memocc(i_stat,mr,'mr',subname)
  allocate(rvalue(n_proj),stat=i_stat)
  call memocc(i_stat,rvalue,'rvalue',subname)
  allocate(zona(n_proj),stat=i_stat)
  call memocc(i_stat,zona,'zona',subname)
  allocate(k_plus_b(n_kpts*n_nnkpts,2),stat=i_stat)
  call memocc(i_stat,k_plus_b,'k_plus_b',subname)
  allocate(G_vec(n_kpts*n_nnkpts,3),stat=i_stat)
  call memocc(i_stat,G_vec,'G_vec',subname)
  allocate(excb(n_excb),stat=i_stat)
  call memocc(i_stat,excb,'excb',subname)
  allocate(rxyz_old(3,atoms%astruct%nat),stat=i_stat)
  call memocc(i_stat,rxyz_old,'rxyz_old',subname)

END SUBROUTINE allocate_initial

subroutine mmnk_calculation_allocation()

      allocate(mmnk_re(orbsb%norb,orbsb%norb,n_kpts*n_nnkpts),stat=i_stat)
      call memocc(i_stat,mmnk_re,'mmnk_re',subname)
      allocate(mmnk_im(orbsb%norb,orbsb%norb,n_kpts*n_nnkpts),stat=i_stat)
      call memocc(i_stat,mmnk_im,'mmnk_im',subname)
      allocate(mmnk_tot(orbsb%norb,n_kpts*n_nnkpts),stat=i_stat)
      call memocc(i_stat,mmnk_tot,'mmnk_tot',subname)
      allocate(psir(nx*ny*nz),stat=i_stat)
      call memocc(i_stat,psir,'psir',subname)
      allocate(psir_re(nx*ny*nz),stat=i_stat)
      call memocc(i_stat,psir_re,'psir_re',subname)
      allocate(psir_im(nx*ny*nz),stat=i_stat)
      call memocc(i_stat,psir_im,'psir_im',subname)
      allocate(psi_daub_re(npsidim),stat=i_stat)
      call memocc(i_stat,psi_daub_re,'psi_daub_re',subname)
      allocate(psi_daub_im(npsidim),stat=i_stat)
      call memocc(i_stat,psi_daub_im,'psi_daub_im',subname)

END SUBROUTINE mmnk_calculation_allocation

subroutine deallocate_projectors()

   i_all = -product(shape(zona))*kind(zona)
   deallocate(zona,stat=i_stat)
   call memocc(i_stat,i_all,'zona',subname)
   i_all = -product(shape(rvalue))*kind(rvalue)
   deallocate(rvalue,stat=i_stat)
   call memocc(i_stat,i_all,'rvalue',subname)
   i_all = -product(shape(l))*kind(l)
   deallocate(l,stat=i_stat)
   call memocc(i_stat,i_all,'l',subname)
   i_all = -product(shape(x_proj))*kind(x_proj)
   deallocate(x_proj,stat=i_stat) 
   call memocc(i_stat,i_all,'x_proj',subname)
   i_all = -product(shape(y_proj))*kind(y_proj)
   deallocate(y_proj,stat=i_stat) 
   call memocc(i_stat,i_all,'y_proj',subname)
   i_all = -product(shape(z_proj))*kind(z_proj)
   deallocate(z_proj,stat=i_stat) 
   call memocc(i_stat,i_all,'z_proj',subname)
   i_all = -product(shape(mr))*kind(mr)
   deallocate(mr,stat=i_stat)
   call memocc(i_stat,i_all,'mr',subname)
   i_all = -product(shape(ctr_proj))*kind(ctr_proj)
   deallocate(ctr_proj,stat=i_stat) 
   call memocc(i_stat,i_all,'ctr_proj',subname) 

end subroutine deallocate_projectors

subroutine deallocate_amnk_calculation()

  i_all = -product(shape(sph_daub))*kind(sph_Daub)
  deallocate(sph_daub,stat=i_stat)
  call memocc(i_stat,i_all,'sph_daub',subname)
  i_all = -product(shape(func_r))*kind(func_r)
  deallocate(func_r,stat=i_stat)
  call memocc(i_stat,i_all,'func_r',subname)
  i_all = -product(shape(ylm))*kind(ylm)
  deallocate(ylm,stat=i_stat)
  call memocc(i_stat,i_all,'ylm',subname)
  i_all = -product(shape(sph_har_etsf))*kind(sph_har_etsf)
  deallocate(sph_har_etsf,stat=i_stat)
  call memocc(i_stat,i_all,'sph_har_etsf',subname)
  i_all = -product(shape(amnk_tot))*kind(amnk_tot)
  deallocate(amnk_tot,stat=i_stat)
  call memocc(i_stat,i_all,'amnk_tot',subname)
  i_all = -product(shape(amnk_bands_sorted))*kind(amnk_bands_sorted)
  deallocate(amnk_bands_sorted,stat=i_stat)
  call memocc(i_stat,i_all,'amnk_bands_sorted',subname)
  i_all = -product(shape(amnk))*kind(amnk)
  deallocate(amnk,stat=i_stat)
  call memocc(i_stat,i_all,'amnk',subname)

END SUBROUTINE deallocate_amnk_calculation

subroutine final_deallocations()

  call deallocate_work_arrays_sumrho(w)
  i_all = -product(shape(psi_etsf))*kind(psi_etsf)
  deallocate(psi_etsf,stat=i_stat)
  call memocc(i_stat,i_all,'psi_etsf',subname)
  i_all = -product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  i_all = -product(shape(psir_re))*kind(psir_re)
  deallocate(psir_re,stat=i_stat)
  call memocc(i_stat,i_all,'psir_re',subname)
  i_all = -product(shape(psir_im))*kind(psir_im)
  deallocate(psir_im,stat=i_stat)
  call memocc(i_stat,i_all,'psir_im',subname)
  i_all = -product(shape(mmnk_tot))*kind(mmnk_tot)
  deallocate(mmnk_tot,stat=i_stat)
  call memocc(i_stat,i_all,'mmnk_tot',subname)
  i_all = -product(shape(psi_etsf2))*kind(psi_etsf2)
  deallocate(psi_etsf2,stat=i_stat)
  call memocc(i_stat,i_all,'psi_etsf2',subname)
  i_all = -product(shape(psi_daub_re))*kind(psi_daub_re)
  deallocate(psi_daub_re,stat=i_stat)
  call memocc(i_stat,i_all,'psi_daub_re',subname)
  i_all = -product(shape(psi_daub_im))*kind(psi_daub_im)
  deallocate(psi_daub_im,stat=i_stat)
  call memocc(i_stat,i_all,'psi_daub_im',subname)
  i_all = -product(shape(virt_list))*kind(virt_list)
  deallocate(virt_list,stat=i_stat)
  call memocc(i_stat,i_all,'virt_list',subname)
  i_all = -product(shape(mmnk_re))*kind(mmnk_re)
  deallocate(mmnk_re,stat=i_stat)
  call memocc(i_stat,i_all,'mmnk_re',subname)
  i_all = -product(shape(mmnk_im))*kind(mmnk_im)
  deallocate(mmnk_im,stat=i_stat)
  call memocc(i_stat,i_all,'mmnk_im',subname)
  i_all = -product(shape(G_vec))*kind(G_vec)
  deallocate(G_vec,stat=i_stat)
  call memocc(i_stat,i_all,'G_vec',subname)
  i_all = -product(shape(k_plus_b))*kind(k_plus_b)
  deallocate(k_plus_b,stat=i_stat)
  call memocc(i_stat,i_all,'k_plus_b',subname) 
  i_all = -product(shape(kpts))*kind(kpts)
  deallocate(kpts,stat=i_stat)
  call memocc(i_stat,i_all,'kpts',subname)
  i_all = -product(shape(excb))*kind(excb)
  deallocate(excb,stat=i_stat)
  call memocc(i_stat,i_all,'excb',subname)

  call deallocate_lr(Glr,subname)
  call deallocate_orbs(orbs,subname)
  call deallocate_comms(comms,subname)
  call deallocate_orbs(orbsv,subname)
  call deallocate_orbs(orbsp,subname)
  call deallocate_comms(commsp,subname) 
  call deallocate_orbs(orbsb,subname)
  call deallocate_comms(commsb,subname) 
  !call deallocate_atoms_scf(atoms,subname)
  call deallocate_atoms(atoms,subname)
!  call free_input_variables(input)
  call bigdft_free_input(input)

END SUBROUTINE final_deallocations
END PROGRAM BigDFT2Wannier

subroutine read_inter_header(iproc,seedname, filetype, residentity, write_resid, n_occ, pre_check,&
           n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)

   ! This routine reads the first lines of a .inter file

   implicit none

   ! I/O variables
   integer, intent(in) :: iproc
   character, intent(out) :: seedname*16, filetype*4
   integer, intent(out) :: n_occ, n_virt, n_virt_tot
   logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check,residentity,write_resid

   ! Local variables
   character :: char1*1, char2*1, char3*1, char4*1
   logical :: file_exist
   integer :: ierr

   ! Should check if it exists, if not, make a nice output message
   inquire(file="input.inter",exist=file_exist)
   if (.not. file_exist) then
      if(iproc == 0) then
         write(*,'(A)') 'ERROR : Input file, input.inter, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   OPEN(11, FILE='input.inter', STATUS='OLD')

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Reading input.inter header :   !'
      write(*,*) '!==================================!'
   end if

   w_unk=.false.
   w_sph=.false.
   w_ang=.false.
   w_rad=.false.

   ! First line
   read(11,*) seedname
   if(iproc==0)write(*,*) 'System studied : ', trim(seedname)

   ! Second line
   read(11,*) filetype
   if(iproc==0)write(*,*) 'file type : ', filetype

   ! Third line
   read(11,*) char1, char2, n_occ
   if(iproc==0)write(*,'(A30,I4)') 'Number of occupied orbitals :', n_occ
   if(char1=='T') then
     residentity = .true.
     if(iproc==0) write(*,*) 'Will use resolution of the identity to construct virtual states'
   else
     residentity = .false.
   end if
   if(residentity .and. char2=='T')then
     write_resid = .true.
     if(iproc==0) write(*,*) 'The constructed virtual states will be written to file.'
   else
     write_resid = .false.
   end if

   ! Fourth line
   read(11,*) char1, n_virt_tot, n_virt
   if (char1=='T' .and. .not. residentity) then
      pre_check=.true.
      if(iproc==0)write(*,*) 'Pre-check before calculating Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A38,I4)') 'Total number of unnocupied orbitals :', n_virt_tot
   else if(.not. residentity)then
      pre_check=.false.
      if(iproc==0)write(*,*) 'Calculation of Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A39,I4)') 'Number of chosen unnocupied orbitals :', n_virt
   else
      pre_check = .false.
      if(iproc==0)write(*,*) 'Calculation of Amnk and Mmnk matrices'
   end if

   ! Fifth line
   read(11,*) char1, char2, char3, char4
   if (char1=='T') then
      w_unk=.true.
      if(iproc==0) write(*,*) 'You want to write a UNKp.s file'
   else if (char1 /= 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_unk'
      STOP
   end if
   if (char2=='T') then
      w_sph=.true.
      if(iproc==0) write(*,*) 'You want to write .cube files for spherical harmonics'
   else if (char2 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_sph'
      STOP
   end if
   if (char3=='T') then
      w_ang=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for angular parts of the spherical harmonics'
   else if (char3 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_ang'
      STOP
   end if
   if (char4=='T') then
      w_rad=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for radial parts of the spherical harmonics'
   else if (char4 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_rad'
      STOP
   end if

   CLOSE(11)

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '! Reading input.inter header done  !'
      write(*,*) '!==================================!'
      print *
      print *
   end if

END SUBROUTINE read_inter_header

!>
subroutine read_inter_list(iproc,n_virt, virt_list)

   ! This routine reads the list of virtual orbitals needed

   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt,iproc
   integer, dimension(n_virt), intent(out) :: virt_list

   ! Local variables
   integer :: i,j,ierr

   open(11, file='input.inter', status='old')

   !   write(*,*) '!==================================!'
   !   write(*,*) '!  Reading virtual orbitals list : !'
   !   write(*,*) '!==================================!'

   do i=1,6
      read(11,*,iostat=ierr) ! Skip first lines                                                                                                                                                               
   end do
   read(11,*,iostat=ierr) (virt_list(j), j=1,n_virt)

   if(ierr < 0) then  !reached the end of file and no virt_list, so generate the trivial one
      do j= 1, n_virt
         virt_list(j) = j
      end do
   end if
   close(11)

   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!Reading virtual orbitals list done!'
      write(*,*) '!==================================!'
      print *
      print *
   end if

END SUBROUTINE read_inter_list

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
   logical :: file_exist

   ! Should check if it exists, if not, make a nice output message
   inquire(file=trim(seedname)//'.nnkp',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'.nnkp, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

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

END SUBROUTINE read_nnkp_int_alloc


!>
subroutine read_nnkp(iproc,seedname, calc_only_A, real_latt, recip_latt, n_kpts, n_proj, &
      &   n_nnkpts, n_excb, kpts, ctr_proj, x_proj, z_proj, &
      &   l, mr, rvalue, zona, k_plus_b, G_vec, excb)

   ! This routine reads an .nnkp file

   implicit none

   ! I/O variables
   character *16 :: seedname
   logical, intent(out) :: calc_only_A
   real, intent(out) :: real_latt(3,3), recip_latt(3,3)
   integer, intent(in) :: n_kpts, n_proj, n_nnkpts, n_excb,iproc
   real, dimension(n_kpts,3), intent(out) :: kpts
   real(kind=8), dimension(n_proj,3), intent(out) :: ctr_proj, x_proj, z_proj
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

END SUBROUTINE read_nnkp


!> This routine returns the angular part of the spherical harmonic identified by indices (l,mr)
!! Calculations are made in spherical coordinates
subroutine angularpart(l, mr, np, nx, ny, nz, ix, iy, iz, &
      &   xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)

   implicit none

   ! I/O variables
   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
   integer, intent(in) :: l(n_proj), mr(n_proj)
   real(kind=8), intent(in) :: xx, yy, zz
   real(kind=8), dimension(nx,ny,nz), intent(out) :: ylm
   real(kind=8), dimension (n_proj,3), intent(in) :: x_proj, y_proj, z_proj

   ! local variables
   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
   real(kind=8), parameter :: eps8  = 1.0e-8
   real(kind=8) :: rr, cost, phi, xdir, ydir, zdir
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

   if(rr < eps8)then
      ylm(ix,iy,iz) = 1.d0/ sqrt(4*pi)
      return
   end if      

   !Here, must define theta in a general fashion, taking into account
   !the z axis: z_proj. To do this, we must project the r vector on z_proj
   ! NOTE: xx,yy,zz are already expressed wrt the projection center
   zdir = z_proj(np,1)*xx + z_proj(np,2)*yy + z_proj(np,3)*zz
   cost = zdir / rr

   !Finally calculate the new x and y
   xdir = x_proj(np,1)*xx + x_proj(np,2)*yy + x_proj(np,3)*zz
   ydir = y_proj(np,1)*xx + y_proj(np,2)*yy + y_proj(np,3)*zz

   if (xdir > eps8) then
      phi = atan( ydir/xdir )
   else if (xdir < -eps8 ) then
      phi = atan( ydir/xdir ) + pi
   else
      phi = sign( pi/2.d0,ydir )
   end if

   if (l(np)==0) then   ! s orbital
      ylm(ix,iy,iz) = s()  
   end if

   if (l(np)==1) then   ! p orbitals
      if (mr(np)==1) ylm(ix,iy,iz) = pz(cost    ) 
      if (mr(np)==2) ylm(ix,iy,iz) = px(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = py(cost,phi)
   end if

   if (l(np)==2) then   ! d orbitals
      if (mr(np)==1) ylm(ix,iy,iz) = dz2(cost    )
      if (mr(np)==2) ylm(ix,iy,iz) = dxz(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = dyz(cost,phi)
      if (mr(np)==4) ylm(ix,iy,iz) = dx2my2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) = dxy(cost,phi)
   endif

   if (l(np)==3) then   ! f orbitals
      if (mr(np)==1) ylm(ix,iy,iz) = fz3(cost)
      if (mr(np)==2) ylm(ix,iy,iz) = fxz2(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = fyz2(cost,phi)
      if (mr(np)==4) ylm(ix,iy,iz) = fzx2my2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) = fxyz(cost,phi)
      if (mr(np)==6) ylm(ix,iy,iz) = fxx2m3y2(cost,phi)
      if (mr(np)==7) ylm(ix,iy,iz) = fy3x2my2(cost,phi)
   endif

   if (l(np)==-1) then  !  sp hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = bs2 * ( s() + px(cost,phi) ) 
      if (mr(np)==2) ylm(ix,iy,iz) = bs2 * ( s() - px(cost,phi) ) 
   end if

   if (l(np)==-2) then  !  sp2 hybrids 
      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s() +2.d0*bs6*px(cost,phi) 
   end if

   if (l(np)==-3) then  !  sp3 hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = 0.5d0*(s()+px(cost,phi)+py(cost,phi)+pz(cost))
      if (mr(np)==2) ylm(ix,iy,iz) = 0.5d0*(s()+px(cost,phi)-py(cost,phi)-pz(cost))
      if (mr(np)==3) ylm(ix,iy,iz) = 0.5d0*(s()-px(cost,phi)+py(cost,phi)-pz(cost))
      if (mr(np)==4) ylm(ix,iy,iz) = 0.5d0*(s()-px(cost,phi)-py(cost,phi)+pz(cost))
   end if

   if (l(np)==-4) then  !  sp3d hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s() +2.d0*bs6*px(cost,phi) 
      if (mr(np)==4) ylm(ix,iy,iz) = bs2*pz(cost)+bs2*dz2(cost)
      if (mr(np)==5) ylm(ix,iy,iz) =-bs2*pz(cost)+bs2*dz2(cost)
   end if

   if (l(np)==-5) then  ! sp3d2 hybrids
      if (mr(np)==1) ylm(ix,iy,iz) = bs6*s()-bs2*px(cost,phi)-bs12*dz2(cost)+.5d0*dx2my2(cost,phi)
      if (mr(np)==2) ylm(ix,iy,iz) = bs6*s()+bs2*px(cost,phi)-bs12*dz2(cost)+.5d0*dx2my2(cost,phi)
      if (mr(np)==3) ylm(ix,iy,iz) = bs6*s()-bs2*py(cost,phi)-bs12*dz2(cost)-.5d0*dx2my2(cost,phi)
      if (mr(np)==4) ylm(ix,iy,iz) = bs6*s()+bs2*py(cost,phi)-bs12*dz2(cost)-.5d0*dx2my2(cost,phi)
      if (mr(np)==5) ylm(ix,iy,iz) = bs6*s()-bs2*pz(cost    )+bs3*dz2(cost)
      if (mr(np)==6) ylm(ix,iy,iz) = bs6*s()+bs2*pz(cost    )+bs3*dz2(cost)
   end if

   contains

   !======== l = 0 =====================================================================
   function s()
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: s
      s = 1.d0/ sqrt(4*pi)
   END FUNCTION s


   !======== l = 1 =====================================================================
   function pz(cost)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: pz, cost
      pz =  sqrt(3.d0/(4*pi)) * cost
   END FUNCTION pz

   function px(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: px, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      px =  sqrt(3.d0/(4*pi)) * sint * cos(phi)
   END FUNCTION px

   function py(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: py, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      py =  sqrt(3.d0/(4*pi)) * sint * sin(phi)
   END FUNCTION py


   !======== l = 2 =====================================================================
   function dz2(cost)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: dz2, cost
      dz2 =  sqrt(1.25d0/(4*pi)) * (3.d0* cost*cost-1.d0)
   END FUNCTION dz2

   function dxz(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) ::  dxz, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      dxz =  sqrt(15.d0/(4*pi)) * sint*cost * cos(phi)
   END FUNCTION dxz

   function dyz(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: dyz, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      dyz =  sqrt(15.d0/(4*pi)) * sint*cost * sin(phi)
   END FUNCTION dyz

   function dx2my2(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: dx2my2, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      dx2my2 =  sqrt(3.75d0/(4*pi)) * sint*sint * cos(2.d0*phi)
   END FUNCTION dx2my2

   function dxy(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: dxy, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      dxy =  sqrt(3.75d0/(4*pi)) * sint*sint * sin(2.d0*phi)
   END FUNCTION dxy


   !======== l = 3 =====================================================================
   function fz3(cost)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fz3, cost
      fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   END FUNCTION fz3

   function fxz2(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fxz2, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   END FUNCTION fxz2

   function fyz2(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fyz2, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   END FUNCTION fyz2

   function fzx2my2(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fzx2my2, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   END FUNCTION fzx2my2

   function fxyz(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fxyz, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   END FUNCTION fxyz

   function fxx2m3y2(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fxx2m3y2, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   END FUNCTION fxx2m3y2

   function fy3x2my2(cost,phi)
      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
      real(kind=8) :: fy3x2my2, cost, phi, sint
      sint = sqrt(abs(1.d0 - cost*cost))
      fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   END FUNCTION fy3x2my2


END SUBROUTINE angularpart

!>
subroutine radialpart(rvalue, zona, np, nx, ny, nz, ix, iy, iz, &
      &   xx, yy, zz, n_proj, func_r)

   ! This routine returns the radial part of the spherical harmonic identified by the indice rvalue.

   implicit none

   ! I/O variables
   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
   integer, intent(in) :: rvalue(n_proj)
   real, intent(in) :: zona(n_proj)
   real(kind=8), intent(in) :: xx, yy, zz
   real(kind=8), dimension(nx,ny,nz), intent(out) :: func_r

   ! local variables 
   !real(kind=8), parameter :: eps8  = 1.0e-8
   real(kind=8) :: rr

   rr = sqrt( xx*xx + yy*yy + zz*zz )

   !   if (rr < eps8) then
   !      write(*,*) 'rr too small '
   !   end if
   if (rvalue(np)==1) func_r(ix,iy,iz) = 2.d0 * (zona(np)/0.529177208)**(3.d0/2.d0) * exp(-(zona(np)/0.529177208)*rr)
   if (rvalue(np)==2) func_r(ix,iy,iz) = 1.d0/sqrt(8.d0) * (zona(np)/0.529177208)**(3.d0/2.d0) * & 
   (2.0d0 - (zona(np)/0.529177208)*rr) * exp(-(zona(np)/0.529177208)*rr*0.5d0)
   if (rvalue(np)==3) func_r(ix,iy,iz) = sqrt(4.d0/27.d0) * (zona(np)/0.529177208)**(3.0d0/2.0d0) * &
      &   (1.d0 - (2.0d0/3.0d0)*(zona(np)/0.529177208)*rr + 2.d0*((zona(np)/0.529177208)*rr)**2/27.d0) * &
      &   exp(-(zona(np)/0.529177208)*rr/3.0d0)

END SUBROUTINE radialpart

! This routine writes .cube files for radial part, angular part and
! the spherical harmonic given in argument
subroutine write_functions(w_sph, w_ang, w_rad, fn1, fn2, fn3, np, Glr, &
      &   hxh, hyh, hzh, atoms, rxyz, sph_har, func_r, ylm)
   use module_types
   implicit none

   ! I/O variables
   logical, intent(in) :: w_sph, w_ang, w_rad                     !logicals controlling the plottings
   character(len=7), intent(in) :: fn1                            !basic name for the spherical harmonics
   character(len=6), intent(in) :: fn2                            !basic name for the radial functions
   character(len=3), intent(in) :: fn3                            !basic name for the ylm
   integer, intent(in) :: np                                      !number of the projector
   real(kind=8), intent(in) :: hxh, hyh, hzh                      !grid spacing
   type(locreg_descriptors), intent(in) :: Glr
   type(atoms_data), intent(in) :: atoms
   real(kind=8), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
   real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i), intent(in) :: sph_har, func_r, ylm
   ! Local variables
   character(len=13) :: subname1
   character(len=12) :: subname2
   character(len=9) :: subname3
   character(len=3) :: np_c
   integer :: i, j, rem, ix, iy, iz
   integer :: nl1, nl2, nl3, nbx, nby, nbz, nc1, nc2, nc3

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

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (atoms%astruct%geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
     nc1=Glr%d%n1i
     nc3=Glr%d%n3i
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
     nc1=Glr%d%n1i-31
     nc3=Glr%d%n3i-31
  end if
  !value of the buffer in the y direction
  if (atoms%astruct%geocode == 'P') then
     nl2=1
     nby = 1
     nc2=Glr%d%n2i
  else
     nl2=15
     nby = 0
     nc2=Glr%d%n2i-31
  end if


   rem=nc3-floor(nc3/6.d0)*6

   ! Write the sph_harxxx.cube files
   if (w_sph) then
      OPEN(12, FILE=trim(subname1)//'.cube', STATUS='unknown')
      write(12,*) ' CUBE file for ISF field'
      write(12,*) ' Case for'
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') atoms%astruct%nat, real(0.d0), real(0.d0), real(0.d0)
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc1, hxh, 0.0_gp, 0.0_gp
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc2, 0.0_gp, hyh, 0.0_gp
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc3, 0.0_gp, 0.0_gp, hzh
      do i=1, atoms%astruct%nat
         write(12,'(I4,1X,F12.6,3(1X,F12.6))') atoms%nzatom(atoms%astruct%iatype(i)), real(0.d0), (real(rxyz(j,i)), j=1,3)
      end do
      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
      do ix=0,nc1-1
         do iy=0,nc2-1 
            do iz=0,nc3-1 
               write(12,'(E14.6)',advance='no') real(sph_har(ix+nl1,iy+nl2,iz+nl3))
               if ( ( (mod(iz+6-rem,6) .eq. 0) .and. (iz+1 .ne. nc3) ) .or. (iz+1 .eq. 1) ) then
                  write(12,'(a)') ''
               end if
            end do
         end do
      end do
      CLOSE(12)
   end if

   ! Write the func_rxxx.cube file
   if (w_rad) then
      OPEN(13, FILE=trim(subname2)//'.cube', STATUS='unknown')
      write(12,*) ' CUBE file for ISF field'
      write(12,*) ' Case for'
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') atoms%astruct%nat, real(0.d0), real(0.d0), real(0.d0)
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc1, hxh, 0.0_gp, 0.0_gp
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc2, 0.0_gp, hyh, 0.0_gp
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc3, 0.0_gp, 0.0_gp, hzh
      do i=1, atoms%astruct%nat
         write(12,'(I4,1X,F12.6,3(1X,F12.6))') atoms%nzatom(atoms%astruct%iatype(i)), real(0.d0), (real(rxyz(j,i)), j=1,3)
      end do
      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
      do ix=0, nc1-1
         do iy=0, nc2-1
            do iz=0, nc3-1
               write(13,'(E14.6)',advance='no') real(func_r(ix+nl1,iy+nl2,iz+nl3))
               if ( ( (mod(iz+6-rem,6) .eq. 0) .and. (iz+1 .ne. nc3) ) .or. (iz+1 .eq. 1) ) then
                  write(13,'(a)') ''
               end if
            end do
         end do
      end do
      CLOSE(13)
   end if

   ! Write the ylmxxx.cube file
   if (w_ang) then
      OPEN(14, FILE=trim(subname3)//'.cube', STATUS='unknown')
      write(12,*) ' CUBE file for ISF field'
      write(12,*) ' Case for'
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') atoms%astruct%nat, real(0.d0), real(0.d0), real(0.d0)
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc1, hxh, 0.0_gp, 0.0_gp
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc2, 0.0_gp, hyh, 0.0_gp
      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nc3, 0.0_gp, 0.0_gp, hzh
      do i=1, atoms%astruct%nat
         write(12,'(I4,1X,F12.6,3(1X,F12.6))') atoms%nzatom(atoms%astruct%iatype(i)), real(0.d0), (real(rxyz(j,i)), j=1,3)
      end do
      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
      do ix=0, nc1-1
         do iy=0, nc2-1
            do iz=0, nc3-1
               write(14,'(E14.6)',advance='no') real(ylm(ix+nl1,iy+nl2,iz+nl3))
               if ( ( (mod(iz+6-rem,6) .eq. 0) .and. (iz+1 .ne. nc3) ) .or. (iz+1 .eq. 1) ) then
                  write(14,'(a)') ''
               end if
            end do
         end do
      end do
      CLOSE(14)
   end if

END SUBROUTINE write_functions




!>
subroutine write_inter(n_virt, nx, ny, nz, amnk_bands_sorted)

   ! This routine writes a new input.inter file from an older one.
   ! It reads the previous input.inter file, and it adds the virt_list at the end of the file.
   ! It automatically switches w_sph, w_ang and w_rad to false

   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt, nx, ny, nz
   integer, dimension(n_virt), intent(inout) :: amnk_bands_sorted

   ! Local variables
   integer :: i
   character :: seedname*20, pre_check_mode*1, filetype*4
   integer :: n_occ, n_virt_tot, ng(3)
   character :: char1*1,char2*1

   ! Read data to keep
   OPEN(11, FILE='input.inter', STATUS='OLD')
   !   write(*,*) '!==================================!'
   !   write(*,*) '!   Writing an input.inter file :  !'
   !   write(*,*) '!==================================!'
   read(11,*) seedname
   read(11,*) filetype
   read(11,*) char1, char2, n_occ
   read(11,*) pre_check_mode, n_virt_tot
   read(11,*) char1
   read(11,*) ng(1), ng(2), ng(3)
   if (pre_check_mode == 'F') then
      read(11,*) (amnk_bands_sorted(i), i=1,n_virt)
   end if
   CLOSE(11)

   ! Write a new input.inter file (change the value of the logical pre_check)
   OPEN(11, FILE='input.inter', STATUS='OLD')
   write(11,'(1a,1x,1a)') seedname, '# Name of the .win file'
   write(11,'(1a,17x,1a)') filetype, '# Format : cube or etsf'
   write(11,'(1a,3x,1a,3x,I4,17x,1a)') 'F','F',n_occ, '# No. of occupied orbitals'
   write(11,'(a1,2(1x,I4),10x,1a)') pre_check_mode, n_virt_tot, n_virt, '# Pre-check, n_virt_tot ' // &
      &   ', n_virt'
   !   if (pre_check_mode == 'T') then
   !      write(11,'(a1,2(1x,I4),10x,1a)') 'F', n_virt_tot, n_virt, '# Pre-check, no. of virtual orbitals&
   !      & for pre-check, no. of virtual orbitals used to write A and M matrices'
   !   else
   !      write(11,'(a1,2(1x,I4),9x,1a)') 'T', n_virt_tot, n_virt, '# Pre-check, no. of virtual orbitals&
   !      & for pre-check, no. of virtual orbitals used to write A and M matrices'
   !   end if
   write(11,'(1a,4x,1a)') char1, 'F    F    F     # Write_UNKp.s, write_spherical_harmonics, ' // &
      &   'write_angular_parts, write_radial_parts'
   write(11,'(3(I4,1x),6x,a)') nx, ny, nz, '# Number of points for each axis in the cubic ' // &
      &   'BigDFT representation (information needed by Wannier90)'
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

END SUBROUTINE write_inter


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

END SUBROUTINE write_amn


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

END SUBROUTINE write_mmn

!>
subroutine write_unk_bin(Glr,orbs,orbsv,orbsb,input,atoms,rxyz,n_occ,n_virt,virt_list,nx,ny,nz,nk,s,iformat)

   use BigDFT_API
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   implicit none
   ! I/O variables
   type(locreg_descriptors), intent(in) :: Glr
   type(orbitals_data), intent(inout) :: orbs,orbsv,orbsb 
   type(atoms_data), intent(in) :: atoms
   type(input_variables), intent(in) :: input
   integer, intent(in) :: n_occ,n_virt,nx,ny,nz,nk,s,iformat
   real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
   integer, dimension (n_virt), intent(in) :: virt_list
   ! Local variables
   logical :: perx,pery,perz
   integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3
   integer :: nb, i, j, k, n_bands, i_stat, ind, i_all
   character :: s_c*1, nk_c*3, seedname*10,filename*60
   character(len=*), parameter :: subname='write_unk_bin'
   real(wp), dimension(nx*ny*nz) :: psir
   real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,n_occ+n_virt) :: psi_etsf
   real(gp), dimension(3,atoms%astruct%nat) :: rxyz_old
   type(workarr_sumrho) :: w

   n_bands=n_occ+n_virt

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

   call split_vectors_for_parallel(0,1,n_occ,orbs)
   call split_vectors_for_parallel(0,1,n_virt+n_occ,orbsb)
   call split_vectors_for_parallel(0,1,n_virt,orbsv)

   call initialize_work_arrays_sumrho(Glr,w)

   ! Read occupied orbitals
   if(n_occ > 0) then
      if(associated(orbs%eval)) nullify(orbs%eval)
      allocate(orbs%eval(n_occ*orbs%nkpts), stat=i_stat)
      call memocc(i_stat,orbs%eval,'orbs%eval',subname)
      filename=trim(input%dir_output) // 'wavefunction'
      call readmywaves(0,filename,iformat,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
      Glr%wfd,psi_etsf(1,1))
      i_all = -product(shape(orbs%eval))*kind(orbs%eval)
      deallocate(orbs%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbs%eval',subname) 
   end if

   ! Read virtual orbitals chosen in pre-check mode 
   if(n_virt > 0) then
      filename=trim(input%dir_output) // 'virtuals'
      if(associated(orbsv%eval)) then
         i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
         deallocate(orbsv%eval,stat=i_stat)
         call memocc(i_stat,i_all,'orbsv%eval',subname)
      end if
      allocate(orbsv%eval(n_virt*orbsv%nkpts), stat=i_stat)
      call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)
      call readmywaves(0,filename,iformat,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
      Glr%wfd,psi_etsf(1,1+n_occ),virt_list)
      i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
      deallocate(orbsv%eval,stat=i_stat)
      call memocc(i_stat,i_all,'orbsv%eval',subname)
   end if

   !calculate buffer shifts
   perx=(Glr%geocode /= 'F')
   pery=(Glr%geocode == 'P')
   perz=(Glr%geocode /= 'F')
   call ext_buffers(perx,nbl1,nbr1)
   call ext_buffers(pery,nbl2,nbr2)
   call ext_buffers(perz,nbl3,nbr3)
   if(nbr1 > 0) nbr1 = nbr1 + 2
   if(nbr2 > 0) nbr2 = nbr2 + 2
   if(nbr3 > 0) nbr3 = nbr3 + 2

   ! Writing the UNKnk.s file
   OPEN(12, FILE=seedname, STATUS='unknown')
   write(*,*) '!==================================!'
   write(*,*) '!      Writing a UNKnk.s file      !'
   write(*,*) '!==================================!'
   write(12,'(I4,4(1X,I4))') nx-(nbl1+nbr1), ny-(nbl2+nbr2), nz-(nbl3+nbr3), nk, n_bands
   do nb=1, n_bands
      ! Convert from Daubechies to cube
      call daub_to_isf(Glr,w,psi_etsf(1,nb),psir)
      do k=nbl3+1, nz-nbr3
         do j=nbl2+1, ny-nbr2
            do i=nbl1+1, nx-nbr1
               ind=(k-1)*ny*nx+(j-1)*nx+i
               write(12,'(E13.6, 1X, E13.6)') psir(ind), 0.d0
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

   call deallocate_work_arrays_sumrho(w)

END SUBROUTINE write_unk_bin

subroutine split_vectors_for_parallel(iproc,nproc,nvctr,orbs)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iproc,nproc
   integer, intent(in) :: nvctr
   type(orbitals_data), intent(inout) :: orbs
   !local variables
   integer :: ntot,jproc,i_stat,i_all
   character(len=*), parameter :: subname='split_vectors_for_parallel'
   integer, dimension(:), allocatable :: nvctr_par,isvctr_par

   ! Initialise the arrays n_proj_par and isproj_par
   allocate(nvctr_par(0:nproc-1),stat=i_stat)
   call memocc(i_stat,nvctr_par,'nvctr_par',subname)
   allocate(isvctr_par(0:nproc-1),stat=i_stat)
   call memocc(i_stat,isvctr_par,'isvctr_par',subname)

   call parallel_repartition_with_kpoints(nproc,1,nvctr,nvctr_par)
   !  call kpts_to_procs_via_obj(nproc,nkpts,nvctr,nvctr_par) 

   ! check the distribution
   ntot=0
   do jproc=0,nproc-1
      isvctr_par(jproc)=ntot
      ntot=ntot+nvctr_par(jproc)
   end do

   orbs%norb = nvctr
   orbs%isorb = isvctr_par(iproc) 
   orbs%norbp = nvctr_par(iproc)

   do jproc=0,nproc-1
      orbs%norb_par(jproc,0) = nvctr_par(jproc)
      !     orbs%isorb_par(jproc) = isvctr_par(jproc) 
   end do

   ! For now, don't worry about kpoints: MUST CHANGE THIS?
   orbs%nkpts=1
   orbs%nspinor=1
   orbs%iskpts=0
   if(associated(orbs%iokpt)) then
      i_all = -product(shape(orbs%iokpt))*kind(orbs%iokpt)
      deallocate(orbs%iokpt,stat=i_stat)
      call memocc(i_stat,i_all,'orbs%iokpt',subname)
   end if
   allocate(orbs%iokpt(orbs%norbp),stat=i_stat)
   call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)
   orbs%iokpt=1

   ! For now, also don't consider spin
   orbs%norbu = nvctr
   orbs%norbd = 0

   i_all = -product(shape(nvctr_par))*kind(nvctr_par)
   deallocate(nvctr_par,stat=i_stat)
   call memocc(i_stat,i_all,'nvctr_par',subname)
   i_all = -product(shape(isvctr_par))*kind(isvctr_par) 
   deallocate(isvctr_par,stat=i_stat)
   call memocc(i_stat,i_all,'isvctr_par',subname)

END SUBROUTINE split_vectors_for_parallel


