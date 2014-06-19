!> @file 
!!   Routines to do diagonalisation (Davidson)
!! @author
!!   Copyright (C) 2010-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Davidsons method for iterative diagonalization of virtual Kohn Sham orbitals
!!   under orthogonality constraints to occupied orbitals psi. The nvirt input
!!   variable gives the number of unoccupied orbitals for which the exit criterion
!!   for the gradients norm holds. nvirte = norbe - norb >= nvirt is the number of
!!   virtual orbitals processed by the method. The dimension of the subspace for
!!   diagonalization is 2*nvirt
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
subroutine constrained_davidson(iproc,nproc,in,at,& 
     orbs,orbsv,nvirt,Lzd,comms,commsv,&
     hx,hy,hz,rxyz,rhopot,psi,v,dpcom,xc,GPU)
  use module_base
  use module_types
  use module_interfaces, except_this_one => constrained_davidson
  use module_xc
  use yaml_output
  use communications, only: transpose_v, untranspose_v
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: nvirt
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: at
  type(local_zone_descriptors), intent(in) :: Lzd
  type(orbitals_data), intent(in) :: orbs
  type(comms_cubic), intent(in) :: comms, commsv
  type(denspot_distribution), intent(in) :: dpcom
  type(xc_info), intent(in) :: xc
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: rhopot
  type(orbitals_data), intent(inout) :: orbsv
  type(GPU_pointers), intent(inout) :: GPU
  real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
                        !v, that is psivirt, is transposed on input and direct on output
  !local variables
  character(len=*), parameter :: subname='davidson'
  logical :: msg,exctX,occorbs !extended output
  integer :: nrhodim,i3rho_add !n(c) occnorb, occnorbu, occnorbd
  integer :: ierr,i_stat,i_all,iorb,jorb,iter,nwork,norb,nspinor,imin
  integer :: ise,ispsi,ikpt,ikptp,nvctrp,ncplx,ncomp,norbs,ispin,ish1,ish2,nspin
  real(gp) :: tt,gnrm,gnrm_fake,emin,diff_max,this_e
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: work,work_rp,hamovr
  real(wp), dimension(:), allocatable :: hv,g,hg,ew  !,Pv,Pg
  real(wp), dimension(:,:,:), allocatable :: e,eg,e_tmp,eg_tmp
  real(wp), dimension(:), pointer :: psiw,psirocc,pot
  real(wp), dimension(:), allocatable :: ALPHAR,ALPHAI,BETA,VR,VL
  type(paw_objects) ::paw !dummy herem, only used for PAW
  
  paw%usepaw=0 !Not using PAW
  call nullify_paw_objects(paw)
  
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
     !n(c) occnorb = 0
     !n(c) occnorbu = 0
     !n(c) occnorbd = 0
  !n(c) else
     !n(c) occnorb = orbs%norb
     !n(c) occnorbu = orbs%norbu
     !n(c) occnorbd = orbs%norbd
  !n(c) end if

  !in the GPU case, the wavefunction should be copied to the card 
  !at each HamiltonianApplication
  !rebind the GPU pointers to the orbsv structure
  if (GPUconv) then
     call free_gpu(GPU,orbs%norbp)
     call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,in%nspin,&
          hx,hy,hz,Lzd%Glr%wfd,orbsv,GPU)
  else if (GPU%OCLconv) then
     call free_gpu_OCL(GPU,orbs,in%nspin)   
     call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
          in%nspin,Lzd%Glr%wfd,orbsv,GPU) 
  end if
 
  GPU%full_locham=.true.
  !verify whether the calculation of the exact exchange term
  !should be performed
  exctX = (xc_exctXfac(xc) /= 0.0_gp)

  !check the size of the rhopot array related to NK SIC
  nrhodim=in%nspin
  i3rho_add=0
  if (in%SIC%approach=='NK') then
     nrhodim=2*nrhodim
     i3rho_add=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*dpcom%nscatterarr(iproc,4)+1
  end if

  !last index of e and hamovr are for mpi_alLzd%Glreduce. 
  !e (eigenvalues) is also used as 2 work arrays
  
  msg=verbose > 2 .and. iproc ==0! no extended output
  !msg =(iproc==0)!extended output

  if(iproc==0)write(*,'(1x,a)')"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  if(iproc==0)write(*,'(1x,a)')"Iterative subspace diagonalization of constrained virtual orbitals."


  !before transposition, create the array of the occupied
  !wavefunctions in real space, for exact exchange calculations
  if (exctX) then
     allocate(psirocc(max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%norbp,&
          dpcom%ngatherarr(0,1)*orbs%norb),1)+ndebug),stat=i_stat)
     call memocc(i_stat,psirocc,'psirocc',subname)

     call prepare_psirocc(iproc,nproc,Lzd%Glr,orbs,dpcom%nscatterarr(iproc,2),dpcom%ngatherarr(0,1),psi,psirocc)
  end if

  
  ! **********************************************
  ! some memory work memory 
  !
  !allocate the work array for transpositions
  if(nproc > 1)then
     allocate(psiw(max(orbs%npsidim_comp,orbs%npsidim_orbs)+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  else
     psiw => null()
  end if
  
  ! allocate and init eval array
  orbsv%eval = f_malloc_ptr(orbsv%norb*orbsv%nkpts,id='orbsv%eval')
  orbsv%eval(1:orbsv%norb*orbsv%nkpts)=-0.5d0
  ! end memory work memory 
  ! **********************************************

  
  !transpose the wavefunction psi if any 
  if (occorbs) then
     call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),psiw(1))
  end if


  ! prepare the v array starting from a set of gaussians
  call psivirt_from_gaussians(iproc,nproc,at,orbsv,Lzd,commsv,rxyz,hx,hy,hz,in%nspin,&
       & v, max(orbsv%npsidim_orbs, orbsv%npsidim_comp))


  ! allocate the potential in the full box
   call full_local_potential(iproc,nproc,orbsv,Lzd,0,dpcom,xc,rhopot,pot)
!!$   call full_local_potential(iproc,nproc,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*dpcom%nscatterarr(iproc,2),&
!!$        Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*dpcom%nscatterarr(iproc,1)*nrhodim,i3rho_add,&
!!$        orbsv,Lzd,0,dpcom%ngatherarr,rhopot,pot)
   
  
  ! **********************************************
  ! some memory work memory 
  allocate(hv(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
  call memocc(i_stat,hv,'hv',subname)

  allocate(e(orbsv%norb,orbsv%nkpts,2+ndebug),stat=i_stat)
  call memocc(i_stat,e,'e',subname)
  
  allocate(eg(orbsv%norb,orbsv%nkpts,2+ndebug),stat=i_stat)
  call memocc(i_stat,eg,'eg',subname)
  
  allocate(e_tmp(orbsv%norb,orbsv%nkpts,2+ndebug),stat=i_stat)
  call memocc(i_stat,eg,'e_tmp',subname)

  allocate(eg_tmp(orbsv%norb,orbsv%nkpts,2+ndebug),stat=i_stat)
  call memocc(i_stat,eg,'eg_tmp',subname)

  if (orbsv%nspinor > 1) then
     ncplx=2
  else
     ncplx=1
  end if
  
  nwork=max(10,16*orbsv%norb)
  allocate(work(ncplx*nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)
  
  allocate(ALPHAR(2*orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,ALPHAR,'ALPHAR',subname)

  allocate(ALPHAI(2*orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,ALPHAI,'ALPHAI',subname)
  
  allocate(BETA(2*orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,BETA,'BETA',subname)
 
  allocate(VR(4*orbsv%norb*orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,VR,'VR',subname)

  allocate(ew(2*orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,ew,'ew',subname)

  allocate(g(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
  call memocc(i_stat,g,'g',subname)

  allocate(hg(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
  call memocc(i_stat,hg,'hg',subname)
  ! end memory work memory 
  ! **********************************************


  
  ! ****************************************
  ! Orthogonality cycle: 
  !
  !   let v orthogonal to occupied state and <vj|P|vi>=delta_ij 
  !
  ! inform
  !
  if (iproc==0) write (*,'(1x,a)',advance="no") "Orthogonality..."
  
  call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)
  !
  ! set v orthogonal to all occupied psi
  !
  if (occorbs) then
     call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
  end if
  !
  ! orthonormalize v through the projections
  !
  call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)
  !
  ! untranspose v 
  !
  call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,v(1),psiw(1))
  !
  ! inform
  !
  !if(iproc==0)write(*,'(1x,a)',advance="no")"done."
  !
  ! End orthogonality cycle: 
  ! ****************************************

 
  ! ****************************************
  ! Hamiltonian application:
  !
  !   compute H|v> => hv, <v|H|v> => e(:,1) and <v|P|v> => e(:,2)
  !
  !stop 'Update HamiltonianApplication call'
!!$  call FullHamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
!!$       proj,Lzd,nlpspd,confdatarr,dpcom%ngatherarr,pot,v,hv,&
!!$       ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,in%SIC,GPU,&
!!$       pkernel,orbs,psirocc)

!!$  call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
!!$       nlpspd,proj,Lzd%Glr,ngatherarr,pot,v,hv,ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,in%SIC,GPU,&
!!$       pkernel,orbs,psirocc) ! optional arguments
  ! 
  !transpose  v and hv
  !
  call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,v(1),psiw(1))
  call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,hv(1),psiw(1))
  !
  ! reset e 
  !
  call to_zero(orbsv%norb*2*orbsv%nkpts,e)
  !
  ! compute rayleigh quotients.
  !
  ! starting index
  ispsi=1
  ! number of components
  nspinor=orbsv%nspinor
  ! loop over kpoint
  do ikptp=1,orbsv%nkptsp
     ! find kpoint starting index
     ikpt=orbsv%iskpts+ikptp
     ! number of coeff for this k point
     nvctrp=commsv%nvctr_par(iproc,ikpt)
     ! cycle if nothing to be done
     if (nvctrp == 0) cycle
     ! otherwise loop on orbitals 
     do iorb=1,orbsv%norb ! temporary variables 
        ! <v|H|v>
        e(iorb,ikpt,1)= dot(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1,hv(ispsi+nvctrp*nspinor*(iorb-1)),1)
        ! <v|P|v>   
        e(iorb,ikpt,2)= dot(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1,v(ispsi+nvctrp*nspinor*(iorb-1)),1)
     end do
     ! increment wf starting index
     ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
  end do
  !
  ! reduce e if necessary
  !
  if(nproc > 1)then
     !sum up the contributions of nproc sets with 
     !commsv%nvctr_par(iproc,1) wavelet coefficients each
     call mpiallred(e(1,1,1),2*orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm)
  end if
  !
  ! inform
  ! 
  !if(iproc==0)write(*,'(1x,a)')"done."
  !
  ! Hamiltonian application:
  ! ****************************************
  
  
  ! ****************************************
  ! Inform
  !
  if(iproc==0)write(*,'(1x,a)')"      1-sqnorm   Rayleigh quotient"
 
  do ikpt=1,orbsv%nkpts
     if (orbsv%nkpts > 1 .and.iproc == 0) write(*,"(1x,A,I3.3,A,3F12.6)") &
          & "Kpt #", ikpt, " BZ coord. = ", orbsv%kpts(:, ikpt)
     do iorb=1,orbsv%norb
        !e(:,1,1) = <psi|H|psi> / <psi|psi>
        e(iorb,ikpt,1)=e(iorb,ikpt,1)/e(iorb,ikpt,2)
        if(iproc==0) write(*,'(1x,i3,1x,1pe13.6,1x,1pe12.5)')&
             iorb,1.0_gp-e(iorb,ikpt,2),e(iorb,ikpt,1)
     end do
  end do
  !
  ! End inform
  ! ****************************************
  


  ! **********************************************
  ! Interaction/Overlap hamovr matrix allocation
  !
  ! calculate the dimension of the overlap matrix for each k-point
  !
  if (orbsv%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if
  !
  ! number of components for the overlap matrix in wp-kind real numbers
  !
  allocate(ndimovrlp(nspin,0:orbsv%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)
  !
  ! fill ndimovrlp
  !
  call dimension_ovrlp(nspin,orbsv,ndimovrlp)
  !
  ! the dimension should be chosen with the max between k-points
  !
  allocate(hamovr(8*ndimovrlp(nspin,orbsv%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,hamovr,'hamovr',subname)
  !
  ! put to zero all the k-points which are not needed
  !
  call to_zero(8*ndimovrlp(nspin,orbsv%nkpts),hamovr)
  !
  ! End interaction/overlap hamovr matrix allocation
  ! **********************************************
  
  

  ! **********************************************
  ! Davidson optimization cycle
  !
  iter=1
  davidson_loop: do 
     ! inform where we are
     if(iproc==0) write( *,'(1x,a,i0)') repeat('~',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
     
     
          
     ! **********************************************
     ! Gradients:
     !
     !   compute the gradients as: g=H|v>*<v|P|v>-P|v>*<v|H|v>
     !   precondition then orthogonalize gradients with occupied states
     !   compute the squared norm of the projected gradients P|g> => e(:,2)
     !
     ! copy hv => g, as we keep hv for later
     !
     call vcopy(max(orbsv%npsidim_orbs,orbsv%npsidim_comp),hv(1),1,g(1),1)  
     !
     ! substract g = g - v*<v|H|v>/<v|P|v>
     !
     ! starting index
     ispsi=1
     ! loop over kpoints
     do ikptp=1,orbsv%nkptsp
        ! this kpoint starting index
        ikpt=orbsv%iskpts+ikptp
        ! number of coeff for this kpoint
        nvctrp=commsv%nvctr_par(iproc,ikpt)
        ! cycle if nothing to do
        if (nvctrp == 0) cycle
        ! otherwise, loop on orbitals
        do iorb=1,orbsv%norb
           !gradient = hv-e*v
           call axpy(nvctrp*nspinor,-e(iorb,ikpt,1)/e(iorb,ikpt,2),v(ispsi+nvctrp*nspinor*(iorb-1)),1,&
                g(ispsi+nvctrp*nspinor*(iorb-1)),1)
        end do
        ! increment starting index
        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do
     !
     ! orthogonalize with respect to occupied states before computing gradients norm
     !
     if (occorbs) then
        ! set g such that they are orthogonal to all occupied psi
        call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
     end if
     !
     ! compute squared norms of the gradients
     !
     ! starting index
     ispsi=1
     ! loop over kpoints
     do ikptp=1,orbsv%nkptsp
        ! this kpoint starting index
        ikpt=orbsv%iskpts+ikptp
        ! number of coeff for this kpoint
        nvctrp=commsv%nvctr_par(iproc,ikpt)
        ! cycle if nothing to do
        if (nvctrp == 0) cycle
        ! otherwise, loop on orbitals
        do iorb=1,orbsv%norb
           !local contribution to the square norm
           e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)!**2
        end do
        ! increment starting index
        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do
     ! reduce if necessary
     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm)
     end if
     !
     ! untranspose gradients for preconditionning
     !
     call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,g(1),psiw(1))
     !
     ! use the values of the eval for the orbitals used in the preconditioner
     !
     do ikpt=1,orbsv%nkpts
        do iorb=1,orbsv%norb
           orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=e(iorb,ikpt,1)
        end do
     end do
     !
     ! apply preconditionner
     !
     call preconditionall(orbsv,Lzd%Glr,hx,hy,hz,in%ncong,g,gnrm_fake,gnrm_fake)
     !
     ! transpose gradients for orthogonalization and norm computation
     !
     call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,g(1),psiw(1))
     !
     ! orthogonalize with respect to occupied states again
     !
     if (occorbs) then
        ! set g such that they are orthogonal to all occupied psi
        call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
     end if
     !
     ! orthogonalize gradient directions (a bit more stable)
     !
     call orthogonalize(iproc,nproc,orbsv,commsv,g,in%orthpar)
     !
     ! untranspose gradients
     !
     call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,g(1),psiw(1))
     !
     ! End gradients
     ! **********************************************


     
     ! **********************************************
     ! Exit criterion on gradient 
     !
     ! compute maximum gradient
     !
     gnrm=0._dp
     if(msg) call yaml_open_sequence('Squared norm of the (nvirt) gradients')
     !if(msg) write(*,'(1x,a)')"squared norm of the (nvirt) gradients"
     ! loop on kpoints
     do ikpt=1,orbsv%nkpts
        ! single spin case
        !if(msg)write(*,'(1x,a)')'done. Eigenvalues, gnrm'
        if (iproc==0) call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))) // ' BZ coord. = ' // &
             & trim(yaml_toa(orbs%kpts(:, ikpt),fmt='(f12.6)')))
        if (nspin == 1) then
           do iorb=1,orbsv%norb!nvirt
              !to understand whether the sqrt should be placed or not
              !if(iproc==0)write(*,'(1x,i5,1pe22.14,1pe9.2)')iorb,e(iorb,ikpt,1),(e(iorb,ikpt,2))
              if (iproc ==0) then
                 call yaml_sequence(advance='no')
                 call yaml_map('Orbitals',(/ e(iorb,ikpt,1),e(iorb,ikpt,2) /), fmt='(1pe22.14)',advance='no')
                 call yaml_comment(trim(yaml_toa(iorb,fmt='(i4.4)'))) 
              end if
              !if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
              tt=real(e(iorb,ikpt,2),dp)!*orbsv%kwgts(ikpt)
              if (iorb <= nvirt) gnrm=max(gnrm,tt)
           end do
        else
           do iorb=1,min(orbsv%norbu,orbsv%norbd) !nvirt
              if(iproc==0)write(*,'(1x,i5,2(1pe22.14,1pe9.2,2x))')&
                   iorb,e(iorb,ikpt,1),(e(iorb,ikpt,2)),e(iorb+orbsv%norbu,ikpt,1),(e(iorb+orbsv%norbu,ikpt,2))
              tt=(real(e(iorb,ikpt,2),dp)+real(e(iorb+orbsv%norbu,ikpt,2),dp))!*orbsv%kwgts(ikpt)
              if (iorb <= nvirt) gnrm=max(gnrm,tt)
              !if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb+orbsv%norbu,tt
           end do
        end if
     end do
     !
     ! inform
     !    
     if (iproc == 0) write (*,'(1x,a,2(1x,1pe12.5))') "maximum gradient and exit criterion ",gnrm,in%gnrm_cv
     !
     ! exit if convergence has been reached
     ! 
     if(gnrm < in%gnrm_cv) then
        exit davidson_loop! iteration loop
     end if
     !
     ! End exit criterion on gradient 
     ! **********************************************
     
     
     
     ! **********************************************
     ! Subspace expansion
     !
     ! apply hamiltonian on gradients
     !
     !stop 'Luigi should work here.'
!!$     call FullHamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
!!$          proj,Lzd,nlpspd,confdatarr,dpcom%ngatherarr,pot,g,hg,&
!!$          ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,in%SIC,GPU,&
!!$          pkernel,orbs,psirocc)

     !call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
     !     nlpspd,proj,Lzd%Glr,ngatherarr,pot,g,hg,ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,in%SIC,GPU,&
     !     pkernel,orbs,psirocc) 
     !
     ! transpose  g and hg and Pg (v, hv and Pv are aLzd%Glready transposed)
     !
     call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,g(1),psiw(1))
     call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,hg(1),psiw(1))
     !
     ! reset expanded hamiltonian/overlap matrices
     !
     call to_zero(8*ndimovrlp(nspin,orbsv%nkpts),hamovr)
     !
     ! compute expanded hamiltonian/overlap matrices   
     ! store upper triangular part of these matrices only
     ! therefore, element (iorb+nvirte,jorb) is transposed to (j,nvirt+iorb)
     !
     ! staring index
     ispsi=1
     ! loop over k points
     do ikptp=1,orbsv%nkptsp
        ! this k point starting index
        ikpt=orbsv%iskpts+ikptp
        ! loop on spin
        do ispin=1,nspin
           ! get dimensions 
           call orbitals_and_components(iproc,ikpt,ispin,orbsv,commsv,&
                nvctrp,norb,norbs,ncomp,nspinor)
           ! cycle if nothing to be done
           if (nvctrp == 0) cycle
           ! otherwise, fill matrices
           call Davidson_constrained_subspace_hamovr(norb,nspinor,ncplx,nvctrp,&
                hamovr(8*ndimovrlp(ispin,ikpt-1)+1),&
                v(ispsi:),g(ispsi),hv(ispsi),hg(ispsi),v(ispsi:),g(ispsi))
           
!           call Davidson_subspace_hamovr(norb,nspinor,ncplx,nvctrp,&
!                hamovr(8*ndimovrlp(ispin,ikpt-1)+1),&
!                v(ispsi),g(ispsi),hv(ispsi),hg(ispsi))
           
           ! increment starting index
           ispsi=ispsi+nvctrp*norb*nspinor
        end do
     end do
     !
     ! reduce result if necessary 
     !
     if(nproc > 1)then
        call mpiallred(hamovr(1),8*ndimovrlp(nspin,orbsv%nkpts),MPI_SUM,bigdft_mpi%mpi_comm)
     end if
     !
     ! check asymmetry
     !
     diff_max=0.0_gp
     do ikptp=1,orbsv%nkptsp
        ! this k point starting index
        ikpt=orbsv%iskpts+ikptp
        ! loop on spin
        do ispin=1,nspin
           do iorb=1,2*norb
              do jorb=1,iorb
                 ! check for max assymetry
                 diff_max=max(diff_max,abs(hamovr(8*ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*2*norb)&
                      -hamovr(8*ndimovrlp(ispin,ikpt-1)+jorb+(iorb-1)*2*norb)))
                 diff_max=max(diff_max,abs(hamovr(8*ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*2*norb+4*norb*norb)&
                      -hamovr(8*ndimovrlp(ispin,ikpt-1)+jorb+(iorb-1)*2*norb+4*norb*norb)))
                 ! enforce symmetry
                 hamovr(8*ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*2*norb)=&
                      hamovr(8*ndimovrlp(ispin,ikpt-1)+jorb+(iorb-1)*2*norb)
                 hamovr(8*ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*2*norb+4*norb*norb)=&
                      hamovr(8*ndimovrlp(ispin,ikpt-1)+jorb+(iorb-1)*2*norb+4*norb*norb)                         
              end do                 
           end do
        end do
     end do
     !
     ! inform
     !
     !if(iproc==0)write(*,'(1x,a)')"done."
     !
     if ( iproc==0 .and. diff_max>1.0E-12_gp ) then
        print *,'WARNING: Important asymmetry of subspace expansion found:',diff_max
     end if
     !
     ! End subspace expansion
     ! **********************************************



     ! **********************************************
     ! Subspace diagonalization
     !
     ! loop on kpoints
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ! this kpoint starting index
        ikpt=orbsv%iskpts+ikptp
        ! offset in e array 
        ise=0
        ! loop on spin
        do ispin=1,nspin
           ! get dimensions
           call orbitals_and_components(iproc,ikpt,ispin,orbsv,commsv,&
                nvctrp,norb,norbs,ncomp,nspinor)
           ! loop if noting to be done
           if (nvctrp == 0) cycle
           !
           ! starting indices in hamovr arrays
           !
           ish1=8*ndimovrlp(ispin,ikpt-1)+1
           ish2=8*ndimovrlp(ispin,ikpt-1)+4*norbs*norb+1
           !
           ! diagonalization:
           ! 
           ! real case
           if (nspinor == 1) then
              ! diagonalize subspace
              call DGGEV( 'N', 'V', 2*norb, hamovr(ish1), 2*norb, hamovr(ish2), 2*norb, ALPHAR, ALPHAI,&
                         BETA, VL, 2*norb, VR, 2*norb, work(1), nwork, i_stat )
              ! check error status     
              if (i_stat /= 0) write(*,*) &
                   'Error in DGGEV on process ',iproc,', infocode ', i_stat
              !
              ! sort eigen energies and copy eigen states              
              !
              do iorb=1,norb
                !
                ! search for min energy
                !
                imin=1
                emin=1.0E32_wp
                do jorb=1,2*norb
                   ! if real finite eigen energy
                   if ( ALPHAI(jorb)==0.0_gp .and. BETA(jorb)>0.0_gp ) then
                      ! compute this energy
                      this_e=ALPHAR(jorb)/BETA(jorb)
                      ! compare with min
                      if (emin>this_e) then
                         emin=this_e
                         imin=jorb
                      end if
                   end if
                end do
                !
                ! copy best case in first part of hamovr
                !
                call vcopy(2*norb,VR(1+(imin-1)*2*norb),1,hamovr(ish1+(iorb-1)*2*norb),1)
                !
                ! store energy in e and scratch ew(imin)
                ! 
                e(iorb,ikpt,1)=emin
                BETA(imin)=0.0_gp
              end do
           ! complex case
           else
              ! allocate work array for complex arithmetic
              allocate(work_rp(6*norb+1+ndebug),stat=i_stat)
              call memocc(i_stat,work_rp,'work_rp',subname)
              ! call lapack GEVP
              call hegv(1,'V','U',2*norb,hamovr(ish1),2*norb,hamovr(ish2),2*norb,&
                   ew(1),work(1),nwork,work_rp(1),i_stat)! Lapack GEVP
              ! check error status
              if (i_stat /= 0) write(*,*) &
                   'Error in HEGV on process ',iproc,', infocode ', i_stat
              ! free memory
              i_all=-product(shape(work_rp))*kind(work_rp)
              deallocate(work_rp,stat=i_stat)
              call memocc(i_stat,i_all,'work_rp',subname)
           end if
           !
           ! keep eigen values
           !
           do iorb=1,norb
              e(iorb+ise,ikpt,1)=ew(iorb)
              if (msg) e(iorb+ise,ikpt,2)=ew(iorb+norb)
           end do
           !
           ! increment e array offset
           !
           ise=norb
        end do
     end do
     !
     ! End subspace diagonalization
     ! **********************************************
     

      
     ! **********************************************
     ! Form new eigen states and gradients
     !
     ! starting index in psi array
     ispsi=1
     ! loop over this proc kpoints
     do ikptp=1,orbsv%nkptsp
        ! index of this kpoint 
        ikpt=orbsv%iskpts+ikptp
        ! number of coeff
        nvctrp=commsv%nvctr_par(iproc,ikpt)
        ! loop on spin
        do ispin=1,nspin
           !
           ! update_psivirt works for kpoint too
           !          
           call update_psivirt(norb,nspinor,ncplx,nvctrp,&
                hamovr(ish1),v(ispsi),g(ispsi),hv(ispsi))
           ! increment starting index
           ispsi=ispsi+nvctrp*norb*nspinor       
        end do
     end do
     !
     ! End form new eigen states
     ! **********************************************
           

     
     ! ****************************************
     ! Orthogonality cycle: 
     !
     ! inform
     !
     if (iproc==0) write (*,'(1x,a)',advance="no") "Orthogonality..."
     !
     ! set v orthogonal to all occupied psi
     !
     if (occorbs) then
        call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
     end if
     !
     ! orthonormalize v set
     !
     call orthogonalize(iproc,nproc,orbsv,commsv,v,in%orthpar)
     !
     ! untranspose v 
     !
     call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,v(1),psiw(1))
     !
     ! inform
     !
     !if(iproc==0)write(*,'(1x,a)')"done."
     !
     ! End orthogonality cycle: 
     ! ****************************************



     ! ****************************************
     ! Hamiltonian application:
     !
     !   compute H|v> => hv 
     !
     !stop 'Here again'
!!$     call FullHamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
!!$          proj,Lzd,nlpspd,confdatarr,dpcom%ngatherarr,pot,v,hv,&
!!$          ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,in%SIC,GPU,&
!!$          pkernel,orbs,psirocc)

     !call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
     !     nlpspd,proj,Lzd%Glr,ngatherarr,pot,v,hv,ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,in%SIC,GPU,&
     !     pkernel,orbs,psirocc)
     !if(iproc==0)write(*,'(1x,a)')"done."
     ! 
     !transpose  v and hv
     !
     call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,v(1),psiw(1))
     call transpose_v(iproc,nproc,orbsv,lzd%glr%wfd,commsv,hv(1),psiw(1))
     !
     ! compute rayleigh quotients. 
     !
     call to_zero(orbsv%norb*2*orbsv%nkpts,e)
     ! starting index
     ispsi=1
     ! number of components
     nspinor=orbsv%nspinor
     ! loop over kpoint
     do ikptp=1,orbsv%nkptsp
        ! find kpoint starting index
        ikpt=orbsv%iskpts+ikptp
        ! number of coeff for this k point
        nvctrp=commsv%nvctr_par(iproc,ikpt)
        ! cycle if nothing to be done
        if (nvctrp == 0) cycle
        ! otherwise loop on orbitals 
        do iorb=1,orbsv%norb ! temporary variables    
           ! <v|H|v>
           e(iorb,ikpt,1) = dot(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1,hv(ispsi+nvctrp*nspinor*(iorb-1)),1)
           ! <v|P|v>   
           e(iorb,ikpt,2) = dot(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1,v(ispsi+nvctrp*nspinor*(iorb-1)),1)
        end do
        ! increment wf starting index
        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do
     !
     ! reduce e if necessary
     !
     if(nproc > 1)then
        !sum up the contributions of nproc sets with 
        !commsv%nvctr_par(iproc,1) wavelet coefficients each
        call mpiallred( e(1,1,1),2*orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm)
        call mpiallred(eg(1,1,1),2*orbsv%norb*orbsv%nkpts,MPI_SUM,bigdft_mpi%mpi_comm)
     end if
     !
     ! End Hamiltonian application:
     ! ****************************************
          
               
          
     ! ****************************************
     ! Inform 
     !
     if(iproc==0)write(*,'(1x,a)')"      1-sqnorm   Rayleigh quotient"

     do ikpt=1,orbsv%nkpts
        if (orbsv%nkpts > 1 .and.iproc == 0) write(*,"(1x,A,I3.3,A,3F12.6)") &
             & "Kpt #", ikpt, " BZ coord. = ", orbsv%kpts(:, ikpt)
        do iorb=1,orbsv%norb
           !e(:,1,1) = <psi|H|psi> / <psi|psi>
           e(iorb,ikpt,1)=e(iorb,ikpt,1)/e(iorb,ikpt,2)
           if(iproc==0) write(*,'(1x,i3,1x,1pe13.6,1x,1pe12.5)')&
                iorb,1.0_gp-e(iorb,ikpt,2),e(iorb,ikpt,1)
        end do
     end do
     !
     ! End hamiltonian application
     ! ****************************************


     
     ! **********************************************
     ! Exit criterion on iteration number 
     !
     iter=iter+1
     if (iter>in%itermax) then 
        if(iproc==0) call yaml_warning( &
          &   'No convergence within the allowed number of minimization steps (itermax)')
        !if(iproc==0) write(*,'(1x,a)')&
        !     'No convergence within the allowed number of minimization steps (itermax)'
        exit davidson_loop
     end if
     !
     ! Exit criterion on iteration number 
     ! **********************************************

  end do davidson_loop


  !retranspose v and psi
  call untranspose_v(iproc,nproc,orbsv,Lzd%Glr%wfd,commsv,v(1),psiw(1))
  call untranspose_v(iproc,nproc,orbs,Lzd%Glr%wfd,comms,psi(1),psiw(1))


  ! inform
  if (iter <=in%itermax) then
     if(iproc==0) call yaml_map('Iteration for Davidson convergence',iter-1)
     !if(iproc==0) write(*,'(1x,a,i3,a)') "Davidson's method: Convergence after ",iter-1,' iterations.'
  end if


  ! Send all eigenvalues to all procs.
  call broadcast_kpt_objects(nproc, orbsv%nkpts, orbsv%norb, e(1,1,1), orbsv%ikptproc)
  

  !copy the values in the eval array of the davidson procedure
  do ikpt=1,orbsv%nkpts
    do iorb=1,orbsv%norb
      orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=e(iorb,ikpt,1)
    end do
  end do


  !write the results on the screen
  call write_eigen_objects(iproc,occorbs,nspin,nvirt,in%nplot,hx,hy,hz,at,rxyz,Lzd%Glr,orbs,orbsv,psi,v,in%output_wf_format)


  ! ******************************************************
  ! memory work  
  ! deallocate potential
  i_all=-product(shape(hg))*kind(hg)
  deallocate(hg,stat=i_stat)
  call memocc(i_stat,i_all,'hg',subname)

  i_all=-product(shape(g))*kind(g)
  deallocate(g,stat=i_stat)
  call memocc(i_stat,i_all,'g',subname)

  call free_full_potential(dpcom%mpi_env%nproc,0,xc,pot,subname)

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
  if(exctX)then
     i_all=-product(shape(psirocc))*kind(psirocc)
     deallocate(psirocc,stat=i_stat)
     call memocc(i_stat,i_all,'psirocc',subname)
  end if

  if(nproc > 1) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  end if

  i_all=-product(shape(hv))*kind(hv)
  deallocate(hv,stat=i_stat)
  call memocc(i_stat,i_all,'hv',subname)

  i_all=-product(shape(e))*kind(e)
  deallocate(e,stat=i_stat)
  call memocc(i_stat,i_all,'e',subname)

  i_all=-product(shape(eg))*kind(eg)
  deallocate(eg,stat=i_stat)
  call memocc(i_stat,i_all,'eg',subname)

  i_all=-product(shape(e_tmp))*kind(e_tmp)
  deallocate(e_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'e_tmp',subname)

  i_all=-product(shape(eg_tmp))*kind(eg_tmp)
  deallocate(eg_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'eg_tmp',subname)
  
  
  i_all=-product(shape(ALPHAR))*kind(ALPHAR)
  deallocate(ALPHAR,stat=i_stat)
  call memocc(i_stat,i_all,'ALPHAR',subname)
  
  i_all=-product(shape(ALPHAI))*kind(ALPHAI)
  deallocate(ALPHAI,stat=i_stat)
  call memocc(i_stat,i_all,'ALPHAI',subname)
  
  i_all=-product(shape(BETA))*kind(BETA)
  deallocate(BETA,stat=i_stat)
  call memocc(i_stat,i_all,'BETA',subname)
  
  i_all=-product(shape(VR))*kind(VR)
  deallocate(VR,stat=i_stat)
  call memocc(i_stat,i_all,'VR',subname)

  ! end memory work  
  ! ******************************************************


  if (GPUconv) then
     call free_gpu(GPU,orbsv%norbp)
  else if (GPU%OCLconv) then
     call free_gpu_OCL(GPU,orbsv,in%nspin)
  end if

END SUBROUTINE constrained_davidson


!>   Generate upper triangular matrix in the subspace of Davidson algorithm
subroutine Davidson_constrained_subspace_hamovr(norb,nspinor,ncplx,nvctrp,hamovr,v,g,hv,hg,Pv,Pg)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvctrp,nspinor,ncplx
  real(wp), dimension(nspinor*nvctrp*norb), intent(in) :: v,g,hv,hg,Pv,Pg
  real(wp), dimension(ncplx,2*norb,2*norb,2), intent(out) :: hamovr

  !local variables
  !n(c) character(len=*), parameter :: subname='Davidson_subspace_hamovr'
  integer :: iorb,jorb,icplx,ncomp
  
  if (nspinor == 4) then
     ncomp=2
  else
     ncomp=1
  end if

  !                 <vi | hvj>      <vi | hgj-n>                   <vi |P| vj>      <vi |P| gj-n>
  ! hamovr(i,j,1)=                               ;  hamovr(i,j,2)=  
  !                 <gi-n | hvj>  <gi-n | hgj-n>                   <gi-n |P| vj>  <gi-n |P| gj-n>

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
          max(1,nvctrp),hv(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,1,1,1),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),v(1),&
          max(1,ncomp*nvctrp), &
          hv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,1,1,1),2*norb)
  end if

  !<gi | hvj> 
  if(nspinor==1) then
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),hv(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,1,1),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          hv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,1,1),2*norb)
  end if

  !<gi | hgj>
  if(nspinor==1) then
     call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),hg(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,norb+1,1),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          hg(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,norb+1,1),2*norb)
  end if
    
  !<vi | Pvj> 
  if(nspinor==1) then
     call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,v(1),&
          max(1,nvctrp),Pv(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,1,1,2),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),v(1),&
          max(1,ncomp*nvctrp), &
          Pv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,1,1,2),2*norb)
  end if
     
  !<gi | Pvj> => hsub(:,:,:,5)
  if(nspinor==1) then
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),Pv(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,1,2),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          Pv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,1,2),2*norb)
  end if

  !<gi | Pgj>
  if(nspinor==1) then
     call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),Pg(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,norb+1,2),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          Pg(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,norb+1,2),2*norb)
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

END SUBROUTINE Davidson_constrained_subspace_hamovr
