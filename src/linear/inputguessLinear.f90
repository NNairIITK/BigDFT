!> @file 
!!   Input guess for the linear version
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Input guess wavefunction diagonalization
!! Input wavefunctions are found by a diagonalization in a minimal basis set
!! Each processors write its initial wavefunctions into the wavefunction file
!! The files are then read by readwave
subroutine inputguessConfinement(iproc, nproc, at, input, hx, hy, hz, &
     rxyz, nlpspd, proj, GPU, orbs, tmb, tmblarge, denspot, rhopotold, energs)
  use module_base
  use module_interfaces, exceptThisOne => inputguessConfinement
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx, hy, hz
  type(atoms_data), intent(inout) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(GPU_pointers), intent(inout) :: GPU
  type(input_variables),intent(in) :: input
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  type(orbitals_data),intent(inout) :: orbs
  type(DFT_wavefunction),intent(inout) :: tmb
  type(DFT_wavefunction),intent(inout) :: tmblarge
  type(DFT_local_fields), intent(inout) :: denspot
  real(dp), dimension(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout) :: rhopotold
  type(energy_terms),intent(inout) :: energs

  ! Local variables
  type(gaussian_basis) :: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  integer :: istat,iall,iat,nspin_ig,iorb,nvirt,norbat,iseg,jjorb
  real(gp) :: hxh,hyh,hzh,eks,fnrm,V3prb,x0
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: locrad, density_kernel_compr
  real(wp), dimension(:,:,:), pointer :: psigau
  integer, dimension(:),allocatable :: norbsPerAt, mapping, inversemapping
  logical,dimension(:),allocatable :: covered
  integer, parameter :: nmax=6,lmax=3
  integer :: ist,jorb,iadd,ii,jj,ityp
  integer :: jlr,iiorb
  integer :: infoCoeff, it_coeff_opt
  type(orbitals_data) :: orbs_gauss
  type(GPU_pointers) :: GPUe
  character(len=2) :: symbol
  real(kind=8) :: rcov,rprb,ehomo,amu,pnrm
  real(kind=8) :: neleconf(nmax,0:lmax)                                        
  integer :: nsccode,mxpl,mxchg
  real(8),dimension(:),allocatable :: ham_compr, ovrlp_compr
  type(mixrhopotDIISParameters) :: mixdiis

  call nullify_orbitals_data(orbs_gauss)

  ! Allocate some arrays we need for the input guess.
  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=istat)
  call memocc(istat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=istat)
  call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)
  allocate(mapping(tmb%orbs%norb), stat=istat)
  call memocc(istat, mapping, 'mapping', subname)
  allocate(covered(tmb%orbs%norb), stat=istat)
  call memocc(istat, covered, 'covered', subname)
  allocate(inversemapping(tmb%orbs%norb), stat=istat)
  call memocc(istat, inversemapping, 'inversemapping', subname)


  GPUe = GPU

  ! Spin for inputguess orbitals
  if (input%nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=input%nspin
  end if

  ! Determine how many atomic orbitals we have. Maybe we have to increase this number to more than
  ! its 'natural' value.
  norbat=0
  ist=0
  do iat=1,at%nat
      ii=input%lin%norbsPerType(at%iatype(iat))
      iadd=0
      do 
          ! Count the number of atomic orbitals and increase the number if necessary until we have more
          ! (or equal) atomic orbitals than basis functions per atom.
          jj=1*nint(at%aocc(1,iat))+3*nint(at%aocc(3,iat))+&
               5*nint(at%aocc(7,iat))+7*nint(at%aocc(13,iat))
          if(jj>=ii) then
              ! we have enough atomic orbitals
              exit
          else
              ! add additional orbitals
              iadd=iadd+1
              select case(iadd)
                  case(1) 
                      at%aocc(1,iat)=1.d0
                  case(2) 
                      at%aocc(3,iat)=1.d0
                  case(3) 
                      at%aocc(7,iat)=1.d0
                  case(4) 
                      at%aocc(13,iat)=1.d0
                  case default 
                      write(*,'(1x,a)') 'ERROR: more than 16 basis functions per atom are not possible!'
                      stop
              end select
          end if
      end do
      norbsPerAt(iat)=jj
      norbat=norbat+norbsPerAt(iat)
  end do

  ! This array gives a mapping from the 'natural' orbital distribution (i.e. simply counting up the atoms) to
  ! our optimized orbital distribution (determined by in orbs%inwhichlocreg).
  iiorb=0
  covered=.false.
  do iat=1,at%nat
      do iorb=1,norbsPerAt(iat)
          iiorb=iiorb+1
          ! Search the corresponding entry in inwhichlocreg
          do jorb=1,tmb%orbs%norb
              if(covered(jorb)) cycle
              jlr=tmb%orbs%inwhichlocreg(jorb)
              if( tmb%lzd%llr(jlr)%locregCenter(1)==rxyz(1,iat) .and. &
                  tmb%lzd%llr(jlr)%locregCenter(2)==rxyz(2,iat) .and. &
                  tmb%lzd%llr(jlr)%locregCenter(3)==rxyz(3,iat) ) then
                  covered(jorb)=.true.
                  mapping(iiorb)=jorb
                  exit
              end if
          end do
      end do
  end do

  ! Inverse mapping
  do iorb=1,tmb%orbs%norb
      do jorb=1,tmb%orbs%norb
          if(mapping(jorb)==iorb) then
              inversemapping(iorb)=jorb
              exit
          end if
      end do
  end do

  nvirt=0

  do ityp=1,at%ntypes
     call eleconf(at%nzatom(ityp),at%nelpsp(ityp),symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
     if(4.d0*rprb>input%lin%locrad_type(ityp)) then
         if(iproc==0) write(*,'(3a,es10.2)') 'WARNING: locrad for atom type ',trim(symbol), &
                      ' is too small; minimal value is ',4.d0*rprb
     end if
     if(input%lin%potentialPrefac_lowaccuracy(ityp)>0.d0) then
         x0=(70.d0/input%lin%potentialPrefac_lowaccuracy(ityp))**.25d0
         if(iproc==0) write(*,'(a,a,2es11.2,es12.3)') 'type, 4.d0*rprb, x0, input%lin%locrad_type(ityp)', &
                      trim(symbol),4.d0*rprb, x0, input%lin%locrad_type(ityp)
         V3prb=input%lin%potentialPrefac_lowaccuracy(ityp)*(4.d0*rprb)**4
         if(iproc==0) write(*,'(a,es14.4)') 'V3prb',V3prb
     end if
  end do

  call inputguess_gaussian_orbitals_forLinear(iproc,nproc,tmb%orbs%norb,at,rxyz,nvirt,nspin_ig,&
       at%nat, norbsPerAt, mapping, &
       tmb%orbs,orbs_gauss,norbsc_arr,locrad,G,psigau,eks,input%lin%potentialPrefac_lowaccuracy)

  ! Take inwhichlocreg from tmb (otherwise there might be problems after the restart...
  !do iorb=1,tmb%orbs%norb
  !    orbs_gauss%inwhichlocreg(iorb)=tmb%orbs%onwhichatom(iorb)
  !end do


  ! Grid spacing on fine grid.
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  ! Transform the atomic orbitals to the wavelet basis.
  orbs_gauss%inwhichlocreg=tmb%orbs%inwhichlocreg
  call wavefunction_dimension(tmb%lzd,orbs_gauss)
  call to_zero(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp), tmb%psi(1))
  call gaussians_to_wavelets_new(iproc,nproc,tmb%lzd,orbs_gauss,G,&
       psigau(1,1,min(tmb%orbs%isorb+1,tmb%orbs%norb)),tmb%psi)

  iall=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=istat)
  call memocc(istat,iall,'psigau',subname)

  call deallocate_gwf(G,subname)

  ! Deallocate locrad, which is not used any longer.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)

  ! Create the potential. First calculate the charge density.
  do iorb=1,tmb%orbs%norb
      tmb%orbs%occup(iorb)=orbs_gauss%occup(iorb)
  end do

  !!call sumrho(denspot%dpbox,tmb%orbs,tmb%lzd,GPUe,at%sym,denspot%rhod,&
  !!     tmb%psi,denspot%rho_psi,inversemapping)
  !!call communicate_density(denspot%dpbox,input%nspin,&!hxh,hyh,hzh,tmbgauss%lzd,&
  !!     denspot%rhod,denspot%rho_psi,denspot%rhov,.false.)

  !Put the Density kernel to identity for now
  !!call to_zero(tmb%orbs%norb**2, tmb%wfnmd%density_kernel(1,1))
  !!do ii = 1, tmb%orbs%norb
  !!   tmb%wfnmd%density_kernel(ii,ii) = 1.d0*tmb%orbs%occup(inversemapping(ii))
  !!end do 

  allocate(density_kernel_compr(tmblarge%mad%nvctr), stat=istat)
  call memocc(istat, density_kernel_compr, 'density_kernel_compr', subname)

  ii=0
  do iseg=1,tmblarge%mad%nseg
      do jorb=tmblarge%mad%keyg(1,iseg),tmblarge%mad%keyg(2,iseg)
          ii=ii+1
          iiorb = (jorb-1)/tmblarge%orbs%norb + 1
          jjorb = jorb - (iiorb-1)*tmblarge%orbs%norb
          if(iiorb==jjorb) then
              density_kernel_compr(ii)=1.d0*tmb%orbs%occup(inversemapping(iiorb))
          else
              density_kernel_compr(ii)=0.0d0
          end if
      end do
  end do


  !Calculate the density in the new scheme
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, tmb%Lzd%hgrids(1), tmb%Lzd%hgrids(2), tmb%Lzd%hgrids(3), &
       tmb%orbs, tmblarge%mad, tmb%collcom_sr, density_kernel_compr, &
       tmb%Lzd%Glr%d%n1i*tmb%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

  !!do istat=1,size(denspot%rhov)
  !!    write(300+iproc,*) istat, denspot%rhov(istat)
  !!end do 
  !!call mpi_finalize(istat)
  !!stop

  iall=-product(shape(density_kernel_compr))*kind(density_kernel_compr)
  deallocate(density_kernel_compr, stat=istat)
  call memocc(istat, iall, 'density_kernel_compr', subname)


  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE &
           .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
      end if
  end if

  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
          call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
      end if
  end if

  if (input%exctxpar == 'OP2P') energs%eexctX = uninitialized(energs%eexctX)


!!  ! Give tmblarge%mad since this is the correct matrix description
!!  call orthonormalizeLocalized(iproc, nproc, -1, tmb%orbs, tmb%lzd, tmblarge%mad, tmb%collcom, &
!!       tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)

  allocate(ham_compr(tmblarge%mad%nvctr), stat=istat)
  call memocc(istat, ham_compr, 'ham_compr', subname)
  allocate(ovrlp_compr(tmblarge%mad%nvctr), stat=istat)
  call memocc(istat, ovrlp_compr, 'ovrlp_compr', subname)

  if (input%lin%scf_mode==LINEAR_FOE) then
      call get_coeff(iproc,nproc,LINEAR_FOE,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
           input%SIC,tmb,fnrm,.true.,.false.,tmblarge, ham_compr, ovrlp_compr, .true., it_coeff_opt)
  else
      call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
           input%SIC,tmb,fnrm,.true.,.false.,tmblarge, ham_compr, ovrlp_compr, .true., it_coeff_opt)
  end if

  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, tmb%Lzd%hgrids(1), tmb%Lzd%hgrids(2), tmb%Lzd%hgrids(3), &
       tmb%orbs, tmblarge%mad, tmb%collcom_sr, tmb%wfnmd%density_kernel_compr, &
       tmb%Lzd%Glr%d%n1i*tmb%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

  call plot_density(iproc,nproc,'initial',at,rxyz,denspot%dpbox,input%nspin,denspot%rhov)


  ! Mix the density.
  if (input%lin%mixing_after_inputguess) then
      if (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE &
           .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
         call mix_main(iproc, nproc, 0, input, tmb%Lzd%Glr, input%lin%alpha_mix_lowaccuracy, &
              denspot, mixdiis, rhopotold, pnrm)
      end if
  end if

  if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE &
       .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
      call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if

  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

  ! Mix the potential.
  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
         call mix_main(iproc, nproc, 0, input, tmb%Lzd%Glr, input%lin%alpha_mix_lowaccuracy, &
              denspot, mixdiis, rhopotold, pnrm)
      end if
  end if

  if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if


  iall=-product(shape(ham_compr))*kind(ham_compr)
  deallocate(ham_compr, stat=istat)
  call memocc(istat, iall, 'ham_compr', subname)
  iall=-product(shape(ovrlp_compr))*kind(ovrlp_compr)
  deallocate(ovrlp_compr, stat=istat)
  call memocc(istat, iall, 'ovrlp_compr', subname)


  ! Important: Don't use for the rest of the code
  tmblarge%can_use_transposed = .false.

  if(associated(tmblarge%psit_c)) then
      iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
      deallocate(tmblarge%psit_c, stat=istat)
      call memocc(istat, iall, 'tmblarge%psit_c', subname)
  end if
  if(associated(tmblarge%psit_f)) then
      iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
      deallocate(tmblarge%psit_f, stat=istat)
      call memocc(istat, iall, 'tmblarge%psit_f', subname)
  end if
  
  if(iproc==0) write(*,'(1x,a)') '------------------------------------------------------------- Input guess generated.'
  
  ! Deallocate all local arrays.

  ! Deallocate all types that are not needed any longer.
  call deallocate_orbitals_data(orbs_gauss, subname)

  ! Deallocate all remaining local arrays.
  iall=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=istat)
  call memocc(istat,iall,'norbsc_arr',subname)

  iall=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt, stat=istat)
  call memocc(istat, iall, 'norbsPerAt',subname)

  iall=-product(shape(mapping))*kind(mapping)
  deallocate(mapping, stat=istat)
  call memocc(istat, iall, 'mapping',subname)

  iall=-product(shape(covered))*kind(covered)
  deallocate(covered, stat=istat)
  call memocc(istat, iall, 'covered',subname)

  iall=-product(shape(inversemapping))*kind(inversemapping)
  deallocate(inversemapping, stat=istat)
  call memocc(istat, iall, 'inversemapping',subname)

END SUBROUTINE inputguessConfinement
