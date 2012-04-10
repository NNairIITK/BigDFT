!!!subroutine calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, atoms, in, hx, hy, hz,&
!!!    comms, lin, nlpspd, proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, phi,&
!!!    coeff, rhopot, fxyz, fnoise, radii_cf)
!!!! Purpose:
!!!! ========
!!!!   Calculates the forces we get with psi. It is copied from cluster, with an additional
!!!!   write statement to print the forces.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   -----------------
!!!!     iproc       process ID
!!!!     nproc       total number of processes
!!!!     n3p         ??
!!!!     i3s         ??
!!!!     i3xcsh      ??
!!!!     Glr         type describing the localization region
!!!!     orbs        type describing the physical orbitals psi
!!!!     atoms       type containing the parameters for the atoms
!!!!     in          type  containing some very general parameters
!!!!     lin         type containing parameters for the linear version
!!!!     nlpspd      nonlocal peudopotential descriptors
!!!!     ngatherarr  ??
!!!!     GPU         parameters for GPUs?
!!!!     irrzon      ??
!!!!     phnons      ??
!!!!     pkernel     ??
!!!!     rxyz        atomic positions
!!!!     fion        ionic forces
!!!!     fdisp       dispersion forces
!!!!   Input / Output arguments
!!!!   ------------------------
!!!!     proj        ??
!!!!     nscatterarr ??
!!!!     GPU         parameters for GPUs?
!!!!   Output arguments:
!!!!   -----------------
!!!!     fxyz        the forces
!!!!     fnoise      noise of the forces
!!!!
!!!use module_base
!!!use module_types
!!!use Poisson_Solver
!!!use module_interfaces, exceptThisOne => calculateForcesSub
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh
!!!real(gp),intent(in):: hx, hy, hz
!!!type(locreg_descriptors),intent(in):: Glr
!!!type(orbitals_data),intent(in):: orbs
!!!type(atoms_data),intent(in):: atoms
!!!type(input_variables),intent(in):: in
!!!type(communications_arrays),intent(in):: comms
!!!type(linearParameters),intent(inout):: lin
!!!type(nonlocal_psp_descriptors),intent(in) :: nlpspd
!!!real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
!!!integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr   !!! NOT NEEDED
!!!integer,dimension(0:nproc-1,4),intent(inout) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!type(GPU_pointers),intent(inout):: GPU
!!!integer,dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
!!!real(dp),dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
!!!real(dp),dimension(lin%as%size_pkernel),intent(in):: pkernel
!!!real(8),dimension(3,atoms%nat),intent(in):: rxyz, fion, fdisp
!!!real(8),dimension(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1)),intent(in):: rhopot
!!!real(8),dimension(3,atoms%nat),intent(out):: fxyz
!!!real(8),intent(out):: fnoise
!!!real(8),dimension(max(lin%gorbs%npsidim_orbs,lin%gorbs%npsidim_comp)),intent(inout):: phi
!!!!real(8),dimension(lin%gorbs%npsidim),intent(inout):: phi
!!!real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
!!!real(gp), dimension(atoms%ntypes,3+ndebug), intent(in) :: radii_cf
!!!! Local variables
!!!integer:: jproc, i_stat, i_all, iat, ierr, j
!!!real(8):: hxh, hyh, hzh, ehart_fake
!!!real(kind=8), dimension(:), allocatable :: rho
!!!real(gp), dimension(:,:), allocatable :: gxyz, fxyzConf
!!!real(kind=8), dimension(:,:,:,:), allocatable :: pot
!!!character(len=*),parameter:: subname='calculateForcesSub'
!!!logical:: refill_proj
!!!integer :: iels, ilr, ii, iorb, jorb
!!!real(wp) :: sum_psi
!!!
!!!  hxh=0.5d0*hx
!!!  hyh=0.5d0*hy
!!!  hzh=0.5d0*hz
!!!
!!!  if (iproc==0) then
!!!     write( *,'(1x,a)')&
!!!          '----------------------------------------------------------------- Forces Calculation'
!!!  end if
!!!
!!!  ! Selfconsistent potential is saved in rhopot, 
!!!  ! new arrays rho,pot for calculation of forces ground state electronic density
!!!
!!!  ! Potential from electronic charge density
!!!
!!!  !manipulate scatter array for avoiding the GGA shift
!!!  do jproc=0,nproc-1
!!!     !n3d=n3p
!!!     nscatterarr(jproc,1)=nscatterarr(jproc,2)
!!!     !i3xcsh=0
!!!     nscatterarr(jproc,4)=0
!!!  end do
!!!
!!!
!!!  ! The charge density has already been calculated and is in rhopot.
!!!  !!!if (n3p>0) then
!!!  !!   allocate(rho(Glr%d%n1i*Glr%d%n2i*n3p*in%nspin+ndebug),stat=i_stat)
!!!  !!   call memocc(i_stat,rho,'rho',subname)
!!!  !!else
!!!  !!   allocate(rho(1+ndebug),stat=i_stat)
!!!  !!   call memocc(i_stat,rho,'rho',subname)
!!!  !!end if
!!!  !!call sumrho(iproc,nproc,orbs,Glr,0,hxh,hyh,hzh,psi,rho,Glr%d%n1i*Glr%d%n2i*n3p,&
!!!  !!        nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!!
!!!  !calculate the total density in the case of nspin==2
!!!  if (in%nspin==2) then
!!!     call axpy(Glr%d%n1i*Glr%d%n2i*n3p,1.0_dp,rho(1+Glr%d%n1i*Glr%d%n2i*n3p),1,rho(1),1)
!!!  end if
!!!  if (n3p>0) then
!!!     allocate(pot(Glr%d%n1i,Glr%d%n2i,n3p,1+ndebug),stat=i_stat)
!!!     call memocc(i_stat,pot,'pot',subname)
!!!  else
!!!     allocate(pot(1,1,1,1+ndebug),stat=i_stat)
!!!     call memocc(i_stat,pot,'pot',subname)
!!!  end if
!!!
!!!  !calculate electrostatic potential
!!!  !call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rho,1,pot,1)
!!!  call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rhopot,1,pot,1)
!!!  call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!       Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,pot,pkernel,pot,ehart_fake,0.0_dp,.false.)
!!!
!!!
!!!  allocate(gxyz(3,atoms%nat+ndebug),stat=i_stat)
!!!  call memocc(i_stat,gxyz,'gxyz',subname)
!!!
!!!  call timing(iproc,'Forces        ','ON')
!!!  ! calculate local part of the forces gxyz
!!!  call local_forces(iproc,atoms,rxyz,hxh,hyh,hzh,&
!!!       Glr%d%n1,Glr%d%n2,Glr%d%n3,n3p,i3s+i3xcsh,Glr%d%n1i,Glr%d%n2i,rhopot,pot,gxyz)
!!!  !call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!
!!!  !!i_all=-product(shape(rho))*kind(rho)
!!!  !!deallocate(rho,stat=i_stat)
!!!  !!call memocc(i_stat,i_all,'rho',subname)
!!!  i_all=-product(shape(pot))*kind(pot)
!!!  deallocate(pot,stat=i_stat)
!!!  call memocc(i_stat,i_all,'pot',subname)
!!!
!!!  if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)',advance='no')'Calculate nonlocal forces...'
!!!
!!!  !refill projectors for tails, davidson
!!!  !refill_proj=(in%calc_tail .or. DoDavidson) .and. DoLastRunThings
!!!  refill_proj=.false.  !! IS THIS CORRECT??
!!!  !gxyz = 0.0_wp
!!!  !fxyz = 0.0_wp
!!!  !call nonlocal_forces(iproc,Glr,hx,hy,hz,atoms,rxyz,&
!!!  !     orbs,nlpspd,proj,Glr%wfd,psi,gxyz,refill_proj)
!!!
!!!  ! ATTENTION: passing phi (after proj, before gxyz) is just to pass something of the right size
!!!  call Linearnonlocal_forces(iproc, nproc, lin%lzd, nlpspd, hx, hy, hz, atoms, rxyz, orbs, &
!!!       proj, phi, gxyz, .false., lin%orbs, coeff, phi)
!!!
!!!!#####################################################################
!!!!DEBUG
!!!!!     
!!!!!sum_psi = sum(psi)
!!!!!call mpiallred(sum_psi,1,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!!if(iproc==0) print *,'sum(psi)',sum_psi
!!!!!
!!!!!     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!!    call mpiallred(gxyz(1,1),3*atoms%nat,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!
!!!!!    if(iproc==0) then
!!!!!    do i_all=1,atoms%nat
!!!!!    open(44,file='Force_ref.dat',status='unknown')
!!!!!    write(*,'(a,i0,a,i0,2x,3(3x,es15.6))') '(C)Forces on atom ',i_all,' :',iproc,fxyz(:,i_all)
!!!!!    write(44,*)'Forces on atom',i_all,' :',iproc,fxyz(:,i_all)
!!!!!    end do
!!!!!    end if
!!!!!
!!!!!    allocate(lin%Lzd%Lnlpspd(lin%Lzd%nlr),stat=i_stat)
!!!!!    do i_all=1,lin%Lzd%nlr
!!!!!       ! allocate projflg
!!!!!       allocate(lin%Lzd%Llr(i_all)%projflg(atoms%nat),stat=i_stat)
!!!!!       call memocc(i_stat,lin%Lzd%Llr(i_all)%projflg,'Lzd%Llr(ilr)%projflg',subname)
!!!!!       call nlpspd_to_locreg(in,iproc,lin%Lzd%Glr,lin%Lzd%Llr(i_all),rxyz,atoms,orbs,&
!!!!!        &      radii_cf,in%frmult,in%frmult,hx,hy,hz,lin%Lzd%Gnlpspd,lin%Lzd%Lnlpspd(i_all),lin%Lzd%Llr(i_all)%projflg)
!!!!!    end do
!!!!!
!!!!!!    sum_psi = 0.0 
!!!!!!    do iorb = 1,orbs%norb
!!!!!!       iels=1
!!!!!!       do jorb=1,lin%orbs%norbp
!!!!!!          ilr=lin%orbs%inwhichlocreg(jorb+lin%orbs%isorb)
!!!!!!          do ii=1,lin%Lzd%LLr(ilr)%wfd%nvctr_c+7*lin%Lzd%LLr(ilr)%wfd%nvctr_f
!!!!!!             sum_psi = sum_psi + coeff(jorb+lin%orbs%isorb,iorb)*phi(iels)
!!!!!!             iels=iels+1
!!!!!!          end do
!!!!!!       end do
!!!!!       !print *,'iproc,iorb,iels,size(phi)',iproc,iorb,iels-1,size(phi)
!!!!!!    end do
!!!!!    !print *,'orbs%npsidim,lin%orbs%npsidim',orbs%npsidim,lin%orbs%npsidim,lin%Lzd%Glr%wfd%nvctr_c+7*lin%Lzd%Glr%wfd%nvctr_f
!!!!!!    call mpiallred(sum_psi,1,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!!    if(iproc==0) print *,'linTMO:sum(psi)',sum_psi
!!!!!
!!!!!    proj = 0.0_wp
!!!!!    gxyz = 0.0_wp
!!!!!    fxyz = 0.0_wp
!!!!!    call Linearnonlocal_forces(iproc,nproc,lin%lzd,hx,hy,hz,atoms,rxyz,orbs,proj,psi,gxyz,.false.,lin%orbs,coeff,phi)
!!!!!!    call mpiallred(gxyz(1,1),3*atoms%nat,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!    call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!
!!!!!    if(iproc==0) then
!!!!!    open(44,file='Force.dat',status='unknown')
!!!!!    do i_all=1,atoms%nat
!!!!!    write(*,'(a,i0,a,i0,2x,3(3x,es15.6))') '(L)Forces on atom ',i_all,' :',iproc,fxyz(:,i_all)
!!!!!    write(44,*)'Forces on atom',i_all,' :',iproc,fxyz(:,i_all)
!!!!!    end do
!!!!!    end if
!!!!!    call mpi_finalize(ierr)
!!!!!    stop
!!!!!
!!!!END DEBUG
!!!!#############################################################################
!!!
!!!  if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)')'done.'
!!!
!!!  ! Add up all the force contributions
!!!  if (nproc > 1) then
!!!     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!     !do iat=1,atoms%nat
!!!     !   fxyz(1,iat)=gxyz(1,iat)
!!!     !   fxyz(2,iat)=gxyz(2,iat)
!!!     !   fxyz(3,iat)=gxyz(3,iat)
!!!     !enddo
!!!  else
!!!     do iat=1,atoms%nat
!!!        fxyz(1,iat)=gxyz(1,iat)
!!!        fxyz(2,iat)=gxyz(2,iat)
!!!        fxyz(3,iat)=gxyz(3,iat)
!!!     enddo
!!!  end if
!!!
!!!  !!$  if (iproc == 0) then
!!!  !!$     sumx=0.d0 ; sumy=0.d0 ; sumz=0.d0
!!!  !!$     fumx=0.d0 ; fumy=0.d0 ; fumz=0.d0
!!!  !!$     do iat=1,atoms%nat
!!!  !!$        sumx=sumx+fxyz(1,iat) ; sumy=sumy+fxyz(2,iat) ; sumz=sumz+fxyz(3,iat)
!!!  !!$        fumx=fumx+fion(1,iat) ; fumy=fumy+fion(2,iat) ; fumz=fumz+fion(3,iat)
!!!  !!$     enddo
!!!  !!$     write(77,'(a30,3(1x,e10.3))') 'translat. force total pot ',sumx,sumy,sumz
!!!  !!$     write(77,'(a30,3(1x,e10.3))') 'translat. force ionic pot ',fumx,fumy,fumz
!!!  !!$  endif
!!!
!!!    !add to the forces the ionic and dispersion contribution 
!!!    do iat=1,atoms%nat
!!!       fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)+fdisp(1,iat)
!!!       fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)+fdisp(2,iat)
!!!       fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)+fdisp(3,iat)
!!!    enddo
!!!
!!!!!!!!!!! TEST !!!
!!!!!!!!  call timing(iproc,'Forces        ','OF')
!!!!!!!!  !!call confinementCorrection()
!!!!!!!!  !!fxyz=fxyz+fxyzConf
!!!!!!!!  !call pulayCorrection()
!!!!!!!!  call timing(iproc,'Forces        ','ON')
!!!!!!!!!!!!!!!!!!!!
!!!
!!!    !i_all=-product(shape(fion))*kind(fion)
!!!    !deallocate(fion,stat=i_stat)
!!!    !call memocc(i_stat,i_all,'fion',subname)
!!!    !i_all=-product(shape(fdisp))*kind(fdisp)
!!!    !deallocate(fdisp,stat=i_stat)
!!!    !call memocc(i_stat,i_all,'fdisp',subname)
!!!    i_all=-product(shape(gxyz))*kind(gxyz)
!!!    deallocate(gxyz,stat=i_stat)
!!!    call memocc(i_stat,i_all,'gxyz',subname)
!!!
!!!
!!!  !!do iat=1,atoms%nat
!!!  !!   if(iproc==0) write(*,'(a,i0,3es14.5)') 'forces for atom ',iat, fxyz(1,iat), fxyz(2,iat), fxyz(3,iat)
!!!  !!end do
!!!
!!!  !subtraction of zero of the forces, disabled for the moment
!!!  !the zero of the forces depends on the atomic positions
!!!  !if (in%gaussian_help .and. .false.) then
!!!  call clean_forces(iproc,atoms,rxyz,fxyz,fnoise)
!!!  !end if
!!!
!!!  if(iproc==0) then
!!!      write(*,'(1x,a)') 'Force values for all atoms in x, y, z direction.'
!!!      do iat=1,atoms%nat
!!!         write(*,'(3x,i0,1x,a6,1x,3(1x,es12.5))') &
!!!              iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
!!!      end do
!!!  end if
!!!
!!!  call timing(iproc,'Forces        ','OF')
!!!
!!!
!!!
!!!
!!!!!contains
!!!!!
!!!!!
!!!!!  subroutine confinementCorrection
!!!!!  implicit none
!!!!!  integer:: ix, iy, iz, ist, istart, jstart, istat, iall, iorb, iiAt, ii, jj, jx, jy, jz, jx0, jy0, jz0
!!!!!  integer:: jorb, nvctrp, korb, lorb, lproc
!!!!!  real(8):: dx, dy, dz, dr2, prefac, tt, phirsq, ttmax, ttmax2, jj2, dnrm2, kx, ky, kz, ddot
!!!!!  real(8),dimension(:),allocatable:: phir, hphirx, hphiry, hphirz, hphix, hphiy, hphiz
!!!!!  real(8),dimension(:,:),allocatable:: gxyz
!!!!!  real(8),dimension(:,:,:),allocatable:: matx, maty, matz
!!!!!  type(workarr_sumrho):: w
!!!!!  type(workarr_locham):: w_lh
!!!!!  real(8),dimension(:),pointer:: phiWork
!!!!!
!!!!!
!!!!!    allocate(phiWork(max(size(phi),size(psi))), stat=istat)
!!!!!    call memocc(istat, phiWork, 'phiWork', subname)
!!!!!
!!!!!    !write(*,*) 'Glr%d%n1,Glr%d%n2,Glr%d%n3', Glr%d%n1,Glr%d%n2,Glr%d%n3
!!!!!    !write(*,*) 'Glr%d%n1i,Glr%d%n2i,Glr%d%n3i', Glr%d%n1i,Glr%d%n2i,Glr%d%n3i
!!!!!    allocate(fxyzConf(3,atoms%nat), stat=istat)
!!!!!    allocate(gxyz(3,atoms%nat), stat=istat)
!!!!!    fxyzConf=0.d0
!!!!!    gxyz=0.d0
!!!!!
!!!!!    allocate(matx(lin%orbs%norb,lin%orbs%norb,2), stat=istat)
!!!!!    allocate(maty(lin%orbs%norb,lin%orbs%norb,2), stat=istat)
!!!!!    allocate(matz(lin%orbs%norb,lin%orbs%norb,2), stat=istat)
!!!!!    matx=0.d0
!!!!!    maty=0.d0
!!!!!    matz=0.d0
!!!!!
!!!!!    allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
!!!!!    allocate(hphirx(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
!!!!!    allocate(hphiry(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
!!!!!    allocate(hphirz(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
!!!!!    allocate(hphix(lin%orbs%npsidim), stat=istat)
!!!!!    allocate(hphiy(lin%orbs%npsidim), stat=istat)
!!!!!    allocate(hphiz(lin%orbs%npsidim), stat=istat)
!!!!!    hphix=0.d0
!!!!!    hphiy=0.d0
!!!!!    hphiz=0.d0
!!!!!    call initialize_work_arrays_sumrho(Glr,w)
!!!!!    istart=1
!!!!!    orbLoop: do iorb=1,lin%orbs%norbp
!!!!!      call deallocate_work_arrays_sumrho(w)
!!!!!      call initialize_work_arrays_sumrho(Glr,w)
!!!!!      phir=0.d0
!!!!!      call daub_to_isf(Glr,w,phi(istart),phir(1))
!!!!!      iiAt=lin%onWhichAtom(iorb)
!!!!!      prefac=lin%potentialPrefac(atoms%iatype(iiAt))
!!!!!
!!!!!      dr2=0.d0
!!!!!      ist=0
!!!!!      !do iz=1,Glr%d%n3i
!!!!!      do iz=1+15,Glr%d%n3i-15
!!!!!        !do iy=1,Glr%d%n2i
!!!!!        do iy=1+15,Glr%d%n2i-15
!!!!!          !do ix=1,Glr%d%n1i
!!!!!          do ix=1+15,Glr%d%n1i-15
!!!!!            ist=ist+1
!!!!!            dx=hxh*ix-rxyz(1,iiAt)
!!!!!            dy=hyh*iy-rxyz(2,iiAt)
!!!!!            dz=hzh*iz-rxyz(3,iiAt)
!!!!!            !dr2=dr2+dx**2+dy**2+dz**2
!!!!!            dr2=dx**2+dy**2+dz**2
!!!!!            tt=4.d0*prefac*dr2
!!!!!            !phirsq=phir(ist)*phir(ist)
!!!!!            !gxyz(1,iiAt)=gxyz(1,iiAt)+tt*phirsq*dx
!!!!!            !gxyz(2,iiAt)=gxyz(2,iiAt)+tt*phirsq*dy
!!!!!            !gxyz(3,iiAt)=gxyz(3,iiAt)+tt*phirsq*dz
!!!!!            hphirx(ist)=tt*phir(ist)*dx
!!!!!            hphiry(ist)=tt*phir(ist)*dy
!!!!!            hphirz(ist)=tt*phir(ist)*dz
!!!!!          end do
!!!!!        end do
!!!!!      end do
!!!!!      kx=lin%orbs%kpts(1,lin%orbs%iokpt(iorb))
!!!!!      ky=lin%orbs%kpts(2,lin%orbs%iokpt(iorb))
!!!!!      kz=lin%orbs%kpts(3,lin%orbs%iokpt(iorb))
!!!!!      call initialize_work_arrays_locham(Glr,orbs%nspinor,w_lh)
!!!!!      call isf_to_daub(hx, hy, hz, kx, ky, kz, orbs%nspinor, Glr, w_lh, hphirx(1), hphix(istart), tt)
!!!!!      call isf_to_daub(hx, hy, hz, kx, ky, kz, orbs%nspinor, Glr, w_lh, hphiry(1), hphiy(istart), tt)
!!!!!      call isf_to_daub(hx, hy, hz, kx, ky, kz, orbs%nspinor, Glr, w_lh, hphirz(1), hphiz(istart), tt)
!!!!!      call deallocate_work_arrays_locham(Glr,w_lh)
!!!!!
!!!!!      istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor
!!!!!    end do orbLoop
!!!!!
!!!!!    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
!!!!!    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphix, work=phiWork)
!!!!!    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphiy, work=phiWork)
!!!!!    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphiz, work=phiWork)
!!!!!
!!!!!    nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
!!!!!    jstart=1
!!!!!    do jorb=1,lin%orbs%norb
!!!!!        istart=1
!!!!!        do iorb=1,lin%orbs%norb
!!!!!            matx(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphix(jstart), 1)
!!!!!            maty(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphiy(jstart), 1)
!!!!!            matz(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphiz(jstart), 1)
!!!!!            istart=istart+nvctrp
!!!!!        end do
!!!!!        jstart=jstart+nvctrp
!!!!!    end do
!!!!!    call mpi_allreduce(matx(1,1,2), matx(1,1,1), lin%orbs%norb**2, mpi_double_precision, &
!!!!!        mpi_sum, mpi_comm_world, ierr)
!!!!!    call mpi_allreduce(maty(1,1,2), maty(1,1,1), lin%orbs%norb**2, mpi_double_precision, &
!!!!!        mpi_sum, mpi_comm_world, ierr)
!!!!!    call mpi_allreduce(matz(1,1,2), matz(1,1,1), lin%orbs%norb**2, mpi_double_precision, &
!!!!!        mpi_sum, mpi_comm_world, ierr)
!!!!!
!!!!!    call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
!!!!!
!!!!!
!!!!!    gxyz=0.d0
!!!!!    do iat=1,atoms%nat
!!!!!        do iorb=1,orbs%norb
!!!!!            do korb=1,lin%orbs%norb
!!!!!                do lproc=0,nproc-1
!!!!!                    do lorb=1,lin%orbs%norb_par(lproc)
!!!!!                        if(iproc==lproc) then
!!!!!                            if(lin%onWhichAtom(lorb)==iat) then
!!!!!                                gxyz(1,iat)=gxyz(1,iat)+coeff(korb,iorb)*coeff(lorb,iorb)*matx(korb,lorb,1)
!!!!!                                gxyz(2,iat)=gxyz(2,iat)+coeff(korb,iorb)*coeff(lorb,iorb)*maty(korb,lorb,1)
!!!!!                                gxyz(3,iat)=gxyz(3,iat)+coeff(korb,iorb)*coeff(lorb,iorb)*matz(korb,lorb,1)
!!!!!                            end if
!!!!!                        end if
!!!!!                    end do
!!!!!                end do
!!!!!            end do
!!!!!        end do
!!!!!    end do
!!!!!
!!!!!    fxyzConf=0.d0
!!!!!    call mpi_allreduce(gxyz(1,1), fxyzConf(1,1), 3*atoms%nat, mpi_double_precision, &
!!!!!        mpi_sum, mpi_comm_world, ierr)
!!!!!
!!!!!
!!!!!    if(iproc==0) write(*,*) 'OLD FORCES'
!!!!!    do iat=1,atoms%nat
!!!!!      if(iproc==0) write(*,'(3es15.6)') fxyz(1,iat), fxyz(2,iat), fxyz(3,iat)
!!!!!    end do
!!!!!
!!!!!    if(iproc==0) write(*,*) 'NEW FORCES'
!!!!!    do iat=1,atoms%nat
!!!!!      if(iproc==0) write(*,'(3es15.6)') fxyzConf(1,iat), fxyzConf(2,iat), fxyzConf(3,iat)
!!!!!    end do
!!!!!
!!!!!    call deallocate_work_arrays_sumrho(w)
!!!!!
!!!!!    iall=-product(shape(phiWork))*kind(phiWork)
!!!!!    deallocate(phiWork, stat=istat)
!!!!!    call memocc(istat, iall, 'phiWork', subname)
!!!!!
!!!!!  end subroutine confinementCorrection
!!!
!!!
!!!
!!!
!!!
!!!!!*************************************************************************************************************88
!!!!!*************************************************************************************************************88
!!!!!*************************************************************************************************************88
!!!
!!!
!!!!!subroutine pulayCorrection()
!!!!!implicit none
!!!!!
!!!!!real(8),dimension(:),allocatable:: dpsix, dpsiy, dpsiz, drhox, drhoy, drhoz, dpot_ionx, dpot_iony, dpot_ionz
!!!!!real(8),dimension(:),allocatable:: dpotx, dpoty, dpotz, rho, pot
!!!!!real(8),dimension(:,:,:,:),allocatable:: dpotxcx, dpotxcy, dpotxcz, potxc
!!!!!real(8),dimension(:,:,:),allocatable:: ovrlpx, ovrlpy, ovrlpz
!!!!!integer,dimension(:,:),allocatable:: ndimovrlp
!!!!!real(8),dimension(:),pointer:: psiWork
!!!!!integer:: iorb, jorb, istat, nvctrp, istart, jstart, ierr, i
!!!!!real(8):: tt, tt2, ddot, deexcux, dvexcux, deexcuy, dvexcuy, deexcuz, dvexcuz, dehartx, deharty, dehartz, psoffset
!!!!!real(8):: eexcu, vexcu, ehart
!!!!!real(8),dimension(:),pointer:: rhocore
!!!!!character(len=3),parameter:: PSquiet='yes'
!!!!!
!!!!!
!!!!!allocate(dpsix(orbs%npsidim), stat=istat)
!!!!!call memocc(istat, dpsix, 'dpsix', subname)
!!!!!allocate(dpsiy(orbs%npsidim), stat=istat)
!!!!!call memocc(istat, dpsiy, 'dpsiy', subname)
!!!!!allocate(dpsiz(orbs%npsidim), stat=istat)
!!!!!call memocc(istat, dpsiz, 'dpsiz', subname)
!!!!!allocate(psiWork(orbs%npsidim), stat=istat)
!!!!!call memocc(istat, psiWork, 'psiWork', subname)
!!!!!
!!!!!allocate(ovrlpx(orbs%norb,orbs%norb,2), stat=istat)
!!!!!call memocc(istat, ovrlpx, 'ovrlpx', subname)
!!!!!allocate(ovrlpy(orbs%norb,orbs%norb,2), stat=istat)
!!!!!call memocc(istat, ovrlpy, 'ovrlpy', subname)
!!!!!allocate(ovrlpz(orbs%norb,orbs%norb,2), stat=istat)
!!!!!call memocc(istat, ovrlpz, 'ovrlpz', subname)
!!!!!
!!!!!allocate(drhox(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, drhox, 'drhox', subname)
!!!!!allocate(drhoy(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, drhoy, 'drhoy', subname)
!!!!!allocate(drhoz(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, drhoz, 'drhoz', subname)
!!!!!allocate(rho(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, rho, 'rho', subname)
!!!!!
!!!!!allocate(dpotx(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, dpotx, 'dpotx', subname)
!!!!!allocate(dpoty(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, dpoty, 'dpoty', subname)
!!!!!allocate(dpotz(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, dpotz, 'dpotz', subname)
!!!!!allocate(pot(lin%as%size_rhopot), stat=istat)
!!!!!call memocc(istat, pot, 'pot', subname)
!!!!!
!!!!!allocate(dpotxcx(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)), stat=istat)
!!!!!call memocc(istat, dpotxcx, 'dpotxcx', subname)
!!!!!allocate(dpotxcy(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)), stat=istat)
!!!!!call memocc(istat, dpotxcy, 'dpotxcy', subname)
!!!!!allocate(dpotxcz(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)), stat=istat)
!!!!!call memocc(istat, dpotxcz, 'dpotxcz', subname)
!!!!!allocate(potxc(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)), stat=istat)
!!!!!call memocc(istat, potxc, 'potxc', subname)
!!!!!
!!!!!allocate(dpot_ionx(lin%as%size_pot_ion),stat=i_stat)
!!!!!call memocc(i_stat,dpot_ionx,'dpot_ionx',subname)
!!!!!allocate(dpot_iony(lin%as%size_pot_ion),stat=i_stat)
!!!!!call memocc(i_stat,dpot_iony,'dpot_iony',subname)
!!!!!allocate(dpot_ionz(lin%as%size_pot_ion),stat=i_stat)
!!!!!call memocc(i_stat,dpot_ionz,'dpot_ionz',subname)
!!!!!
!!!!!
!!!!!! Create an input guess for dpsi. At the moment only random.
!!!!!call random_number(dpsix)
!!!!!call random_number(dpsiy)
!!!!!call random_number(dpsiz)
!!!!!dpsix=psi
!!!!!dpsiy=psi
!!!!!dpsiz=psi
!!!!!do i=1,orbs%npsidim
!!!!!    call random_number(tt)
!!!!!    tt=tt-.5d0
!!!!!    tt=tt*.1d0
!!!!!    dpsix(i)=dpsix(i)*tt
!!!!!    call random_number(tt)
!!!!!    tt=tt-.5d0
!!!!!    tt=tt*.1d0
!!!!!    dpsiy(i)=dpsiy(i)*tt
!!!!!    call random_number(tt)
!!!!!    tt=tt-.5d0
!!!!!    tt=tt*.1d0
!!!!!    dpsiz(i)=dpsiz(i)*tt
!!!!!end do
!!!!!
!!!!!
!!!!!!!allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=istat)
!!!!!!!call memocc(istat, ndimovrlp, 'ndimovrlp', subname)
!!!!!!!! Allocate the overlap matrix
!!!!!!!allocate(ovrlp(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=istat)
!!!!!!!call memocc(i_stat, ovrlp, 'ovrlp', subname)
!!!!!
!!!!!! Othogonalize dpsi to psi
!!!!!call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=psiWork)
!!!!!nvctrp=comms%nvctr_par(iproc,1) ! 1 for k-point
!!!!!istart=1
!!!!!do iorb=1,orbs%norb
!!!!!    jstart=1
!!!!!    do jorb=1,orbs%norb
!!!!!        ovrlpx(jorb,iorb,2)=ddot(nvctrp, dpsix(istart), 1, psi(jstart), 1)
!!!!!        !call daxpy(nvctrp, -tt, psi(jstart), 1, dpsix(istart), 1)
!!!!!        ovrlpy(jorb,iorb,2)=ddot(nvctrp, dpsiy(istart), 1, psi(jstart), 1)
!!!!!        !call daxpy(nvctrp, -tt, psi(jstart), 1, dpsix(istart), 1)
!!!!!        ovrlpz(jorb,iorb,2)=ddot(nvctrp, dpsiz(istart), 1, psi(jstart), 1)
!!!!!        !call daxpy(nvctrp, -tt, psi(jstart), 1, dpsix(istart), 1)
!!!!!        jstart=jstart+nvctrp
!!!!!    end do
!!!!!    istart=istart+nvctrp
!!!!!end do
!!!!!
!!!!!call mpi_allreduce(ovrlpx(1,1,2), ovrlpx(1,1,1), orbs%norb**2, mpi_double_precision, &
!!!!!    mpi_sum, mpi_comm_world, ierr)
!!!!!call mpi_allreduce(ovrlpy(1,1,2), ovrlpy(1,1,1), orbs%norb**2, mpi_double_precision, &
!!!!!    mpi_sum, mpi_comm_world, ierr)
!!!!!call mpi_allreduce(ovrlpz(1,1,2), ovrlpz(1,1,1), orbs%norb**2, mpi_double_precision, &
!!!!!    mpi_sum, mpi_comm_world, ierr)
!!!!!
!!!!!istart=1
!!!!!do iorb=1,orbs%norb
!!!!!    jstart=1
!!!!!    do jorb=1,orbs%norb
!!!!!        call daxpy(nvctrp, -ovrlpx(jorb,iorb,1), psi(jstart), 1, dpsix(istart), 1)
!!!!!        call daxpy(nvctrp, -ovrlpy(jorb,iorb,1), psi(jstart), 1, dpsiy(istart), 1)
!!!!!        call daxpy(nvctrp, -ovrlpz(jorb,iorb,1), psi(jstart), 1, dpsiz(istart), 1)
!!!!!        jstart=jstart+nvctrp
!!!!!    end do
!!!!!    istart=istart+nvctrp
!!!!!end do
!!!!!
!!!!!
!!!!!! Orthonormalize dpsi
!!!!!call orthogonalize(iproc, nproc, orbs, comms, Glr%wfd, dpsix, in)
!!!!!call orthogonalize(iproc, nproc, orbs, comms, Glr%wfd, dpsiy, in)
!!!!!call orthogonalize(iproc, nproc, orbs, comms, Glr%wfd, dpsiz, in)
!!!!!call dscal(orbs%norb*nvctrp, .1d0, dpsix, 1)
!!!!!call dscal(orbs%norb*nvctrp, .1d0, dpsiy, 1)
!!!!!call dscal(orbs%norb*nvctrp, .1d0, dpsiz, 1)
!!!!!
!!!!!!!! CHECK
!!!!!istart=1
!!!!!do iorb=1,orbs%norb
!!!!!    jstart=1
!!!!!    do jorb=1,orbs%norb
!!!!!        ovrlpx(jorb,iorb,2)=ddot(nvctrp, dpsix(istart), 1, dpsix(jstart), 1)
!!!!!        !call daxpy(nvctrp, -tt, psi(jstart), 1, dpsix(istart), 1)
!!!!!        ovrlpy(jorb,iorb,2)=ddot(nvctrp, dpsiy(istart), 1, dpsiy(jstart), 1)
!!!!!        !call daxpy(nvctrp, -tt, psi(jstart), 1, dpsix(istart), 1)
!!!!!        ovrlpz(jorb,iorb,2)=ddot(nvctrp, dpsiz(istart), 1, dpsiz(jstart), 1)
!!!!!        !call daxpy(nvctrp, -tt, psi(jstart), 1, dpsix(istart), 1)
!!!!!        jstart=jstart+nvctrp
!!!!!    end do
!!!!!    istart=istart+nvctrp
!!!!!end do
!!!!!call mpi_allreduce(ovrlpx(1,1,2), ovrlpx(1,1,1), orbs%norb**2, mpi_double_precision, &
!!!!!    mpi_sum, mpi_comm_world, ierr)
!!!!!call mpi_allreduce(ovrlpy(1,1,2), ovrlpy(1,1,1), orbs%norb**2, mpi_double_precision, &
!!!!!    mpi_sum, mpi_comm_world, ierr)
!!!!!call mpi_allreduce(ovrlpz(1,1,2), ovrlpz(1,1,1), orbs%norb**2, mpi_double_precision, &
!!!!!    mpi_sum, mpi_comm_world, ierr)
!!!!!do iorb=1,orbs%norb
!!!!!    do jorb=1,orbs%norb
!!!!!        if(iproc==0) write(*,*) 'iorb, jorb, ovrlp_x', iorb, jorb, ovrlpx(jorb,iorb,1)
!!!!!        if(iproc==0) write(*,*) 'iorb, jorb, ovrlp_x', iorb, jorb, ovrlpy(jorb,iorb,1)
!!!!!        if(iproc==0) write(*,*) 'iorb, jorb, ovrlp_x', iorb, jorb, ovrlpz(jorb,iorb,1)
!!!!!    end do
!!!!!end do
!!!!!
!!!!!
!!!!!!!!!! ATTENTION !!!
!!!!!!!! DEBUG
!!!!!!!call orthogonalize(iproc, nproc, orbs, comms, Glr%wfd, psi, in)
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=psiWork)
!!!!!call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, dpsix, work=psiWork)
!!!!!call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, dpsiy, work=psiWork)
!!!!!call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, dpsiz, work=psiWork)
!!!!!
!!!!!
!!!!!
!!!!!! Calculate the right hand side of the linear system. First apply the modified Hamiltonian to psi.
!!!!!
!!!!!!!!! ATTENTION
!!!!!!! DEBUG
!!!!!!dpsix=psi
!!!!!!!!!!!!!!!!!!!
!!!!!
!!!!!
!!!!!call sumrho(iproc,nproc,orbs,Glr,in%ixc,hxh,hyh,hzh,psi,rho,&
!!!!!     Glr%d%n1i*Glr%d%n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!!!!if(orbs%nspinor==4) then
!!!!!   !this wrapper can be inserted inside the poisson solver 
!!!!!   stop 'not yet implemented!!!'
!!!!!   !!call PSolverNC(atoms%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3d,&
!!!!!   !!     in%ixc,hxh,hyh,hzh,&
!!!!!   !!     rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
!!!!!else
!!!!!    ! rhocore is not associated... not used?
!!!!!    call XC_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!         Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,hxh,hyh,hzh,&
!!!!!         rho,eexcu,vexcu,in%nspin,rhocore,potxc)
!!!!!
!!!!!     call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!          Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!          rho,pkernel,pot,ehart,0.0_dp,.false.,&
!!!!!          quiet=PSquiet) !optional argument
!!!!!end if
!!!!!
!!!!!do iat=1,atoms%nat
!!!!!
!!!!!    ! First calculate the charge density.
!!!!!    call sumrho(iproc,nproc,orbs,Glr,in%ixc,hxh,hyh,hzh,dpsix,drhox,&
!!!!!         Glr%d%n1i*Glr%d%n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!!!!    call sumrho(iproc,nproc,orbs,Glr,in%ixc,hxh,hyh,hzh,dpsiy,drhoy,&
!!!!!         Glr%d%n1i*Glr%d%n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!!!!    call sumrho(iproc,nproc,orbs,Glr,in%ixc,hxh,hyh,hzh,dpsiz,drhoz,&
!!!!!         Glr%d%n1i*Glr%d%n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!!!!
!!!!!    dpotx=drhox
!!!!!    dpoty=drhoy
!!!!!    dpotz=drhoz
!!!!!    
!!!!!    if(orbs%nspinor==4) then
!!!!!       !this wrapper can be inserted inside the poisson solver 
!!!!!       stop 'not yet implemented!!!'
!!!!!       !!call PSolverNC(atoms%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3d,&
!!!!!       !!     in%ixc,hxh,hyh,hzh,&
!!!!!       !!     rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
!!!!!    else
!!!!!        ! rhocore is not associated... not used?
!!!!!        call XC_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,hxh,hyh,hzh,&
!!!!!             drhox,deexcux,dvexcux,in%nspin,rhocore,dpotxcx)
!!!!!        call XC_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,hxh,hyh,hzh,&
!!!!!             drhoy,deexcuy,dvexcuy,in%nspin,rhocore,dpotxcy)
!!!!!        call XC_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,hxh,hyh,hzh,&
!!!!!             drhoz,deexcux,dvexcux,in%nspin,rhocore,dpotxcz)
!!!!!
!!!!!         call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!              Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!              drhox,pkernel,dpotx,dehartx,0.0_dp,.false.,&
!!!!!              quiet=PSquiet) !optional argument
!!!!!         call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!              Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!              drhoy,pkernel,dpoty,deharty,0.0_dp,.false.,&
!!!!!              quiet=PSquiet) !optional argument
!!!!!         call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!              Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!              drhoz,pkernel,dpotz,dehartz,0.0_dp,.false.,&
!!!!!              quiet=PSquiet) !optional argument
!!!!!    
!!!!!        !! POT_ION IS NOT CORRECT
!!!!!        ! We need the derivative of pot_ion with respect to the atomic
!!!!!        ! coordinates.
!!!!!
!!!!!         ! IS THIS CORRECT?
!!!!!         psoffset=0.d0
!!!!!
!!!!!         !call createIonicPotentialModified(atoms%geocode,iproc,nproc,atoms,lin, iat, rxyz,&
!!!!!         !hxh,hyh,hzh,in%elecfield,Glr%d%n1,Glr%d%n2,Glr%d%n3,n3pi,i3s,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,pkernel, &
!!!!!         !drhox, drhoy, drhoz, dpotx, dpoty, dpotz, dpot_ionx, dpot_iony, dpot_ionz, &
!!!!!         !psoffset,in%nvacancy,&
!!!!!         !in%correct_offset)
!!!!!         call createIonicPotentialModified(atoms%geocode,iproc,nproc,atoms,lin, iat, rxyz,&
!!!!!         hxh,hyh,hzh,in%elecfield,Glr%d%n1,Glr%d%n2,Glr%d%n3,n3pi,i3s,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,pkernel, &
!!!!!         rho, rho, rho, pot, pot, pot, dpot_ionx, dpot_iony, dpot_ionz, &
!!!!!         psoffset,in%nvacancy,&
!!!!!         in%correct_offset)
!!!!!
!!!!!    
!!!!!    
!!!!!         call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!              Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!              drhox,pkernel,dpot_ionx,dehartx,0.0_dp,.true.,&
!!!!!              quiet=PSquiet) !optional argument
!!!!!         call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!              Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!              drhoy,pkernel,dpot_iony,deharty,0.0_dp,.true.,&
!!!!!              quiet=PSquiet) !optional argument
!!!!!         call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!              Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!              drhoz,pkernel,dpot_ionz,dehartz,0.0_dp,.true.,&
!!!!!              quiet=PSquiet) !optional argument
!!!!!    end if
!!!!!
!!!!!end do
!!!!!
!!!!!
!!!!!     !!if(orbs%nspinor==4) then
!!!!!     !!   !this wrapper can be inserted inside the poisson solver 
!!!!!     !!   call PSolverNC(atoms%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3d,&
!!!!!     !!        in%ixc,hxh,hyh,hzh,&
!!!!!     !!        rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
!!!!!     !!else
!!!!!     !!   call XC_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!     !!        Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,hxh,hyh,hzh,&
!!!!!     !!        rhopot,eexcu,vexcu,in%nspin,rhocore,potxc)
!!!!!
!!!!!     !!   call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!!!     !!        Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!     !!        rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.,&
!!!!!     !!        quiet=PSquiet) !optional argument
!!!!!
!!!!!     !!   !sum the two potentials in rhopot array
!!!!!     !!   !fill the other part, for spin, polarised
!!!!!     !!   if (in%nspin == 2) then
!!!!!     !!      call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rhopot(1),1,&
!!!!!     !!           rhopot(1+Glr%d%n1i*Glr%d%n2i*n3p),1)
!!!!!     !!   end if
!!!!!     !!   !spin up and down together with the XC part
!!!!!     !!   call axpy(Glr%d%n1i*Glr%d%n2i*n3p*in%nspin,1.0_dp,potxc(1,1,1,1),1,&
!!!!!     !!        rhopot(1),1)
!!!!!
!!!!!     !!end if
!!!!!end subroutine pulayCorrection
!!!
!!!
!!!
!!!
!!!
!!!
!!!end subroutine calculateForcesSub
!!!


!!!subroutine createIonicPotentialModified(geocode,iproc,nproc,at,lin, iiAt, rxyz,&
!!!     hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel, &
!!!     rhox, rhoy, rhoz, potx, poty, potz, pot_ionx, pot_iony, pot_ionz, &
!!!     psoffset,nvacancy,&
!!!     correct_offset)
!!!  use module_base
!!!  use module_types
!!!!  use module_interfaces, except_this_one => createIonicPotential
!!!  use Poisson_Solver
!!!  implicit none
!!!  character(len=1), intent(in) :: geocode
!!!  logical,intent(in) :: correct_offset
!!!  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,nvacancy, iiAt
!!!  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
!!!  type(atoms_data), intent(in) :: at
!!!  type(linearParameters),intent(in):: lin
!!!  real(gp), intent(in) :: elecfield
!!!  real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!!  real(dp), dimension(*), intent(in) :: pkernel
!!!  real(wp), dimension(lin%as%size_pot_ion), intent(inout) :: pot_ionx, pot_iony, pot_ionz
!!!  real(8),dimension(lin%as%size_rhopot):: rhox, rhoy, rhoz, potx, poty, potz
!!!  !local variables
!!!  character(len=*), parameter :: subname='createIonicPotentialModified'
!!!  logical :: perx,pery,perz,gox,goy,goz,htoobig=.false.,efwrite,check_potion=.false.
!!!  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp,nspin
!!!  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc
!!!  real(kind=8) :: pi,rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
!!!  real(kind=8) :: tt_tot,rholeaked_tot,potxyz,offset, prefactor
!!!  real(8):: tt1, tt2, tt3, potleakedx, potleakedy, potleakedz
!!!  real(8):: rhoelx, rhoely, rhoelz, Velx, Vely, Velz
!!!  real(8):: tt1_tot, tt2_tot, tt3_tot, potleakedx_tot, potleakedy_tot, potleakedz_tot
!!!  real(wp) :: maxdiff
!!!  real(gp) :: ehart
!!!  real(dp), dimension(2) :: charges_mpi
!!!  integer, dimension(:,:), allocatable :: ngatherarr
!!!  real(dp), dimension(:), allocatable :: potion_corr
!!!  real(dp), dimension(:), pointer :: pkernel_ref
!!!  real(kind=8), dimension(4) :: cprime
!!!
!!!  call timing(iproc,'CrtLocPot     ','ON')
!!!
!!!  if (iproc.eq.0) then
!!!     write(*,'(1x,a)')&
!!!          '----------------------------------------------------------- Ionic Potential Creation'
!!!  end if
!!!
!!!  pi=4.d0*atan(1.d0)
!!!  ! Ionic charge (must be calculated for the PS active processes)
!!!  rholeaked=0.d0
!!!  ! Ionic energy (can be calculated for all the processors)
!!!
!!!  !Creates charge density arising from the ionic PSP cores
!!!  call razero(n1i*n2i*n3pi,pot_ionx)
!!!  call razero(n1i*n2i*n3pi,pot_iony)
!!!  call razero(n1i*n2i*n3pi,pot_ionz)
!!!
!!!  !conditions for periodicity in the three directions
!!!  perx=(geocode /= 'F')
!!!  pery=(geocode == 'P')
!!!  perz=(geocode /= 'F')
!!!
!!!  call ext_buffers(perx,nbl1,nbr1)
!!!  call ext_buffers(pery,nbl2,nbr2)
!!!  call ext_buffers(perz,nbl3,nbr3)
!!!
!!!  if (n3pi >0 .and. .not. htoobig) then
!!!
!!!     !do iat=1,at%nat
!!!     do iat=iiAt,iiAt
!!!        ityp=at%iatype(iat)
!!!        rx=rxyz(1,iat) 
!!!        ry=rxyz(2,iat)
!!!        rz=rxyz(3,iat)
!!!
!!!        !building array of coefficients of the derivative of the gaussian part
!!!        cprime(1)=2.d0*at%psppar(0,2,ityp)-at%psppar(0,1,ityp)
!!!        cprime(2)=4.d0*at%psppar(0,3,ityp)-at%psppar(0,2,ityp)
!!!        cprime(3)=6.d0*at%psppar(0,4,ityp)-at%psppar(0,3,ityp)
!!!        cprime(4)=-at%psppar(0,4,ityp)
!!!
!!!        ! determine number of local terms
!!!        nloc=0
!!!        do iloc=1,4
!!!           if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
!!!        enddo
!!!
!!!
!!!
!!!        rloc=at%psppar(0,0,ityp)
!!!        prefactor=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**5)
!!!        charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
!!!        cutoff=10.d0*rloc
!!!
!!!        isx=floor((rx-cutoff)/hxh)
!!!        isy=floor((ry-cutoff)/hyh)
!!!        isz=floor((rz-cutoff)/hzh)
!!!
!!!        iex=ceiling((rx+cutoff)/hxh)
!!!        iey=ceiling((ry+cutoff)/hyh)
!!!        iez=ceiling((rz+cutoff)/hzh)
!!!
!!!        do i3=isz,iez
!!!           z=real(i3,kind=8)*hzh-rz
!!!           call ind_positions(perz,i3,n3,j3,goz) 
!!!           j3=j3+nbl3+1
!!!           do i2=isy,iey
!!!              y=real(i2,kind=8)*hyh-ry
!!!              call ind_positions(pery,i2,n2,j2,goy)
!!!              do i1=isx,iex
!!!                 x=real(i1,kind=8)*hxh-rx
!!!                 call ind_positions(perx,i1,n1,j1,gox)
!!!                 r2=x**2+y**2+z**2
!!!                 arg=r2/rloc**2
!!!                 xp=exp(-.5d0*arg)
!!!                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!!                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
!!!                    !gaussian part
!!!                    tt=0.d0
!!!                    if (nloc /= 0) then
!!!                       !derivative of the polynomial
!!!                       tt=cprime(nloc)
!!!                       do iloc=nloc-1,1,-1
!!!                          tt=arg*tt+cprime(iloc)
!!!                       enddo
!!!                       rhoelx=rhox(ind)
!!!                       rhoely=rhoy(ind)
!!!                       rhoelz=rhoz(ind)
!!!                       tt1=xp*tt*rhoelx
!!!                       tt2=xp*tt*rhoely
!!!                       tt3=xp*tt*rhoelz
!!!                       pot_ionx(ind)=pot_ionx(ind)+tt1*x
!!!                       pot_iony(ind)=pot_iony(ind)+tt2*y
!!!                       pot_ionz(ind)=pot_ionz(ind)+tt3*z
!!!                    end if
!!!                    !error function part
!!!                    Velx=potx(ind)
!!!                    Vely=poty(ind)
!!!                    Velz=potz(ind)
!!!                    pot_ionx(ind)=pot_ionx(ind)+prefactor*xp*Velx*x
!!!                    pot_iony(ind)=pot_iony(ind)+prefactor*xp*Vely*y
!!!                    pot_ionz(ind)=pot_ionz(ind)+prefactor*xp*Velz*z
!!!
!!!                 else if (.not. goz ) then
!!!                    potleakedx=potleakedx+xp*tt*rhox(1)
!!!                    potleakedy=potleakedy+xp*tt*rhoy(1)
!!!                    potleakedz=potleakedz+xp*tt*rhoz(1)
!!!                 endif
!!!              enddo
!!!           enddo
!!!        enddo
!!!
!!!     enddo
!!!
!!!  end if
!!!
!!!  ! IS THIS CORRECT?
!!!  pot_ionx=pot_ionx*hxh*hyh*hzh
!!!  pot_iony=pot_iony*hxh*hyh*hzh
!!!  pot_ionz=pot_ionz*hxh*hyh*hzh
!!!  
!!!
!!!  ! Check
!!!  tt1=0.d0
!!!  tt2=0.d0
!!!  tt3=0.d0
!!!  do j3=1,n3pi
!!!     do i2= -nbl2,2*n2+1+nbr2
!!!        do i1= -nbl1,2*n1+1+nbr1
!!!           ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
!!!           tt1=tt1+pot_ionx(ind)
!!!           tt2=tt2+pot_iony(ind)
!!!           tt3=tt3+pot_ionz(ind)
!!!        enddo
!!!     enddo
!!!  enddo
!!!
!!!  tt1=tt1*hxh*hyh*hzh
!!!  tt2=tt2*hxh*hyh*hzh
!!!  tt3=tt3*hxh*hyh*hzh
!!!  potleakedx=potleakedx*hxh*hyh*hzh
!!!  potleakedy=potleakedy*hxh*hyh*hzh
!!!  potleakedz=potleakedz*hxh*hyh*hzh
!!!
!!!  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt
!!!
!!!  if (nproc > 1) then
!!!     charges_mpi(1)=tt1
!!!     charges_mpi(2)=potleakedx
!!!     call mpiallred(charges_mpi(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!     tt1_tot=charges_mpi(1)
!!!     potleakedx_tot=charges_mpi(2)
!!!
!!!     charges_mpi(1)=tt2
!!!     charges_mpi(2)=potleakedy
!!!     call mpiallred(charges_mpi(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!     tt2_tot=charges_mpi(1)
!!!     potleakedy_tot=charges_mpi(2)
!!!
!!!     charges_mpi(1)=tt3
!!!     charges_mpi(2)=potleakedz
!!!     call mpiallred(charges_mpi(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!     tt3_tot=charges_mpi(1)
!!!     potleakedz_tot=charges_mpi(2)
!!!  else
!!!     tt1_tot=tt
!!!     potleakedx_tot=potleakedx
!!!
!!!     tt2_tot=tt
!!!     potleakedy_tot=potleakedy
!!!
!!!     tt3_tot=tt
!!!     potleakedz_tot=potleakedy
!!!  end if
!!!
!!!  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
!!!       'x: total ionic charge, leaked charge ',tt1_tot,potleakedx_tot
!!!  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
!!!       'y: total ionic charge, leaked charge ',tt2_tot,potleakedy_tot
!!!  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
!!!       'z: total ionic charge, leaked charge ',tt3_tot,potleakedz_tot
!!!
!!!  if (.not. htoobig) then
!!!     call timing(iproc,'CrtLocPot     ','OF')
!!!     !here the value of the datacode must be kept fixed
!!!     nspin=1
!!!
!!!     call H_potential(geocode,'D',iproc,nproc,&
!!!          n1i,n2i,n3i,hxh,hyh,hzh,&
!!!          pot_ionx,pkernel,pot_ionx,ehart,-psoffset,.false.)
!!!     call H_potential(geocode,'D',iproc,nproc,&
!!!          n1i,n2i,n3i,hxh,hyh,hzh,&
!!!          pot_iony,pkernel,pot_iony,ehart,-psoffset,.false.)
!!!     call H_potential(geocode,'D',iproc,nproc,&
!!!          n1i,n2i,n3i,hxh,hyh,hzh,&
!!!          pot_ionz,pkernel,pot_ionz,ehart,-psoffset,.false.)
!!!
!!!     call timing(iproc,'CrtLocPot     ','ON')
!!!     
!!!  ! Commented since check_potion is set to false...
!!!     !!!if (check_potion) then
!!!     !!!   if (iproc == 0) write(*,'(1x,a)',advance='no') &
!!!     !!!        'Check the ionic potential...'
!!!     !!!     
!!!     !!!   allocate(potion_corr(n1i*n2i*n3pi+ndebug),stat=i_stat)
!!!     !!!   call memocc(i_stat,potion_corr,'potion_corr',subname)
!!!
!!!     !!!   call razero(n1i*n2i*n3pi,potion_corr)
!!!
!!!     !!!   !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
!!!     !!!   !for the moment works only in the isolated BC case
!!!     !!!   do i3=1,n3pi
!!!     !!!      z=real(i3+i3s-1-nbl3-1,gp)*hzh
!!!     !!!      do i2=1,n2i
!!!     !!!         y=real(i2-nbl2-1,gp)*hyh
!!!     !!!         do i1=1,n1i
!!!     !!!            x=real(i1-nbl1-1,gp)*hxh
!!!     !!!            ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!!     !!!            !if (i1==49 .and. i2==46 .and. i3==44) then
!!!     !!!               call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
!!!     !!!            !   stop
!!!     !!!            !end if
!!!     !!!            potion_corr(ind)=potion_corr(ind)+potxyz
!!!     !!!            !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
!!!     !!!         end do
!!!     !!!      end do
!!!     !!!   end do
!!!
!!!     !!!   !then calculate the maximum difference in the sup norm
!!!     !!!   maxdiff=0.0_wp
!!!     !!!   do i3=1,n3pi
!!!     !!!      do i2=1,n2i
!!!     !!!         do i1=1,n1i
!!!     !!!            ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!!     !!!            maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
!!!     !!!            !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
!!!     !!!         end do
!!!     !!!      end do
!!!     !!!   end do
!!!
!!!     !!!   call mpiallred(maxdiff,1,MPI_MAX,MPI_COMM_WORLD,ierr)
!!!
!!!     !!!   if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff
!!!
!!!     !!!   stop
!!!
!!!     !!!   i_all=-product(shape(potion_corr))*kind(potion_corr)
!!!     !!!   deallocate(potion_corr,stat=i_stat)
!!!     !!!   call memocc(i_stat,i_all,'potion_corr',subname)
!!!
!!!     !!!end if
!!!
!!!  end if
!!!
!!!
!!!!!!  !calculate the value of the offset to be put
!!!!!!  tt_tot=0.d0
!!!!!!  do ind=1,n1i*n2i*n3i
!!!!!!     tt_tot=tt_tot+pot_ion(ind)
!!!!!!  end do
!!!!!!  print *,'previous offset',tt_tot*hxh*hyh*hzh
!!!
!!!
!!!
!!! ! Already done above.. ?
!!!  !!if (n3pi > 0) then
!!!  !!   do iat=1,at%nat
!!!  !!      ityp=at%iatype(iat)
!!!
!!!  !!      rx=rxyz(1,iat)
!!!  !!      ry=rxyz(2,iat)
!!!  !!      rz=rxyz(3,iat)
!!!
!!!  !!      ! determine number of local terms
!!!  !!      nloc=0
!!!  !!      do iloc=1,4
!!!  !!         if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
!!!  !!      enddo
!!!  !!      rloc=at%psppar(0,0,ityp)
!!!  !!      cutoff=10.d0*rloc
!!!
!!!  !!      isx=floor((rx-cutoff)/hxh)
!!!  !!      isy=floor((ry-cutoff)/hyh)
!!!  !!      isz=floor((rz-cutoff)/hzh)
!!!
!!!  !!      iex=ceiling((rx+cutoff)/hxh)
!!!  !!      iey=ceiling((ry+cutoff)/hyh)
!!!  !!      iez=ceiling((rz+cutoff)/hzh)
!!!  !!      
!!!  !!      !do not add the local part for the vacancy
!!!  !!      if (nloc /= 0) then
!!!
!!!  !!         do i3=isz,iez
!!!  !!            z=real(i3,kind=8)*hzh-rz
!!!  !!            call ind_positions(perz,i3,n3,j3,goz) 
!!!  !!            j3=j3+nbl3+1
!!!  !!            if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
!!!  !!               do i2=isy,iey
!!!  !!                  y=real(i2,kind=8)*hyh-ry
!!!  !!                  call ind_positions(pery,i2,n2,j2,goy)
!!!  !!                  if (goy) then
!!!  !!                     do i1=isx,iex
!!!  !!                        x=real(i1,kind=8)*hxh-rx
!!!  !!                        call ind_positions(perx,i1,n1,j1,gox)
!!!  !!                        if (gox) then
!!!  !!                           r2=x**2+y**2+z**2
!!!  !!                           arg=r2/rloc**2
!!!  !!                           xp=exp(-.5d0*arg)
!!!  !!                           tt=at%psppar(0,nloc,ityp)
!!!  !!                           do iloc=nloc-1,1,-1
!!!  !!                              tt=arg*tt+at%psppar(0,iloc,ityp)
!!!  !!                           enddo
!!!  !!                           ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
!!!  !!                           pot_ion(ind)=pot_ion(ind)+xp*tt
!!!  !!                        end if
!!!  !!                     enddo
!!!  !!                  end if
!!!  !!               enddo
!!!  !!            end if
!!!  !!         end do
!!!
!!!  !!      end if
!!!
!!!  !!   enddo
!!!
!!!  !!   if (htoobig) then
!!!  !!      !add to pot_ion an explicit error function to correct in the case of big grid spacing
!!!  !!      !for the moment works only in the isolated BC case
!!!  !!      do i3=1,n3pi
!!!  !!         z=real(i3+i3s-1-nbl3-1,gp)*hzh
!!!  !!         do i2=1,n2i
!!!  !!            y=real(i2-nbl2-1,gp)*hyh
!!!  !!            do i1=1,n1i
!!!  !!               x=real(i1-nbl1-1,gp)*hxh
!!!  !!               ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!!  !!               call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
!!!  !!               pot_ion(ind)=pot_ion(ind)+potxyz
!!!  !!            end do
!!!  !!         end do
!!!  !!      end do
!!!  !!   end if
!!!  !!   
!!!  !!end if
!!!  
!!!
!!!! nvacancy==0 is always true
!!! !! if (nvacancy /= 0) then
!!! !!    !for a vacancy reference calculation, save the ionic potential to be used
!!! !!    !in the following run
!!!
!!! !!    !first calculate the kernel in isolated BC
!!! !!    call timing(iproc,'CrtLocPot     ','OF')
!!! !!    call createKernel(iproc,nproc,'F',n1i,n2i,n3i,hxh,hyh,hzh,16,pkernel_ref)
!!! !!    call timing(iproc,'CrtLocPot     ','ON')
!!!
!!!
!!! !!    !calculate the ionic potential correction in the global data distribution
!!! !!    allocate(potion_corr(n1i*n2i*n3i+ndebug),stat=i_stat)
!!! !!    call memocc(i_stat,potion_corr,'potion_corr',subname)
!!!
!!! !!    call razero(n1i*n2i*n3i,potion_corr)
!!!
!!! !!    iat=nvacancy
!!! !!    ityp=at%iatype(iat)
!!! !!    rx=rxyz(1,iat) 
!!! !!    ry=rxyz(2,iat)
!!! !!    rz=rxyz(3,iat)
!!!
!!! !!    rloc=at%psppar(0,0,ityp)
!!! !!    charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
!!! !!    cutoff=10.d0*rloc
!!!
!!! !!    isx=floor((rx-cutoff)/hxh)
!!! !!    isy=floor((ry-cutoff)/hyh)
!!! !!    isz=floor((rz-cutoff)/hzh)
!!!
!!! !!    iex=ceiling((rx+cutoff)/hxh)
!!! !!    iey=ceiling((ry+cutoff)/hyh)
!!! !!    iez=ceiling((rz+cutoff)/hzh)
!!!
!!! !!    do i3=isz,iez
!!! !!       z=real(i3,kind=8)*hzh-rz
!!! !!       call ind_positions(perz,i3,n3,j3,goz) 
!!! !!       j3=j3+nbl3+1
!!! !!       do i2=isy,iey
!!! !!          y=real(i2,kind=8)*hyh-ry
!!! !!          call ind_positions(pery,i2,n2,j2,goy)
!!! !!          do i1=isx,iex
!!! !!             x=real(i1,kind=8)*hxh-rx
!!! !!             call ind_positions(perx,i1,n1,j1,gox)
!!! !!             r2=x**2+y**2+z**2
!!! !!             arg=r2/rloc**2
!!! !!             xp=exp(-.5d0*arg)
!!! !!             if (goz  .and. goy  .and. gox ) then
!!! !!                ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-1)*n1i*n2i
!!! !!                potion_corr(ind)=xp*charge !the sign is inverted here
!!! !!             endif
!!! !!          enddo
!!! !!       enddo
!!! !!    enddo
!!!
!!!!!!!$     !plot the ionic potential in a .pot file
!!!!!!!$     !allocate the arrays for plotting
!!!!!!!$     !its values are ignored in the datacode='G' case
!!!!!!!$     allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
!!!!!!!$     call memocc(i_stat,nscatterarr,'nscatterarr',subname)
!!!!!!!$     allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
!!!!!!!$     call memocc(i_stat,ngatherarr,'ngatherarr',subname)
!!!!!!!$     !create the descriptors for the density and the potential
!!!!!!!$     !these descriptors should take into account the localisation regions
!!!!!!!$     call createDensPotDescriptors(iproc,nproc,at%geocode,'D',n1i,n2i,n3i,0,&
!!!!!!!$          n3d_fake,n3p_fake,n3pi_fake,i3xcsh_fake,i3s_fake,nscatterarr,ngatherarr)
!!!!!!!$
!!!!!!!$     i_all=-product(shape(nscatterarr))*kind(nscatterarr)
!!!!!!!$     deallocate(nscatterarr,stat=i_stat)
!!!!!!!$     call memocc(i_stat,i_all,'nscatterarr',subname)
!!!!!!!$
!!!!!!!$
!!!!!!!$     call plot_density(at%geocode,'gaupotion.pot',iproc,1,n1,n2,n3,n1i,n2i,n3i,n3i,&
!!!!!!!$          at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)
!!!
!!!
!!!
!!! !!    call timing(iproc,'CrtLocPot     ','OF')
!!! !!    !here the value of the datacode must be kept fixed
!!! !!    call H_potential('F','G',iproc,nproc,&
!!! !!         n1i,n2i,n3i,hxh,hyh,hzh,&
!!! !!         potion_corr,pkernel_ref,potion_corr,ehart,0.0_gp,.false.)
!!!
!!!!!!!$     call PSolver('F','G',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
!!!!!!!$          potion_corr,pkernel_ref,potion_corr,ehart,eexcu,vexcu,0.0_gp,.false.,1)
!!! !!    call timing(iproc,'CrtLocPot     ','ON')
!!!
!!!
!!! !!    i_all=-product(shape(pkernel_ref))*kind(pkernel_ref)
!!! !!    deallocate(pkernel_ref,stat=i_stat)
!!! !!    call memocc(i_stat,i_all,'pkernel_ref',subname)
!!!
!!!
!!!!!!!!     call plot_density(at%geocode,'deltapotion.pot',iproc,1,n1,n2,n3,n1i,n2i,n3i,n3i,&
!!!!!!!!          at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)
!!!
!!!
!!! !!    iat=nvacancy
!!! !!    ityp=at%iatype(iat)
!!!
!!! !!    rx=rxyz(1,iat)
!!! !!    ry=rxyz(2,iat)
!!! !!    rz=rxyz(3,iat)
!!!
!!! !!    ! determine number of local terms
!!! !!    nloc=0
!!! !!    do iloc=1,4
!!! !!       if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
!!! !!    enddo
!!! !!    rloc=at%psppar(0,0,ityp)
!!! !!    cutoff=10.d0*rloc
!!!
!!! !!    isx=floor((rx-cutoff)/hxh)
!!! !!    isy=floor((ry-cutoff)/hyh)
!!! !!    isz=floor((rz-cutoff)/hzh)
!!!
!!! !!    iex=ceiling((rx+cutoff)/hxh)
!!! !!    iey=ceiling((ry+cutoff)/hyh)
!!! !!    iez=ceiling((rz+cutoff)/hzh)
!!!
!!! !!    !do not add the local part for the vacancy
!!! !!    if (nloc /= 0) then
!!!
!!! !!       do i3=isz,iez
!!! !!          z=real(i3,kind=8)*hzh-rz
!!! !!          call ind_positions(perz,i3,n3,j3,goz) 
!!! !!          j3=j3+nbl3+1
!!! !!          if (goz) then
!!! !!             do i2=isy,iey
!!! !!                y=real(i2,kind=8)*hyh-ry
!!! !!                call ind_positions(pery,i2,n2,j2,goy)
!!! !!                if (goy) then
!!! !!                   do i1=isx,iex
!!! !!                      x=real(i1,kind=8)*hxh-rx
!!! !!                      call ind_positions(perx,i1,n1,j1,gox)
!!! !!                      if (gox) then
!!! !!                         r2=x**2+y**2+z**2
!!! !!                         arg=r2/rloc**2
!!! !!                         xp=exp(-.5d0*arg)
!!! !!                         tt=at%psppar(0,nloc,ityp)
!!! !!                         do iloc=nloc-1,1,-1
!!! !!                            tt=arg*tt+at%psppar(0,iloc,ityp)
!!! !!                         enddo
!!! !!                         ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-1)*n1i*n2i
!!! !!                         potion_corr(ind)=&
!!! !!                              potion_corr(ind)-xp*tt ! the sign has changed here
!!! !!                      end if
!!! !!                   enddo
!!! !!                end if
!!! !!             enddo
!!! !!          end if
!!! !!       end do
!!!
!!! !!    end if
!!!
!!!
!!! !!    !call plot_density(at%geocode,'deltapotion_final.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3pi,&
!!! !!    !     at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)
!!!
!!! !!    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!! !!    !add the periodic pot_ion
!!! !!    ind=1+(i3s-1)*n1i*n2i
!!! !!    call axpy(n1i*n2i*n3pi,1.0_dp,pot_ion(1),1,potion_corr(ind),1)
!!!
!!! !!    if (correct_offset) then
!!! !!       !calculate the offset
!!! !!       tt=0.d0
!!! !!       do i1=1,n3pi*n2i*n1i
!!! !!          tt=tt+potion_corr(ind-1+i1)
!!! !!       enddo
!!! !!       !tt=tt*hxh*hyh*hzh
!!! !!       if (nproc > 1) then
!!! !!          call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
!!! !!               MPI_SUM,MPI_COMM_WORLD,ierr)
!!! !!       else
!!! !!          offset=tt
!!! !!       end if
!!!
!!! !!       if (iproc==0) print *,'offset for potion',offset
!!!
!!! !!       !now potion_corr has zero integral
!!! !!       potion_corr=potion_corr-offset/real(n1i*n2i*n3i,dp)
!!!
!!! !!       !calculate the offset
!!! !!       tt=0.d0
!!! !!       do i1=1,n3pi*n2i*n1i
!!! !!          tt=tt+potion_corr(ind-1+i1)
!!! !!       enddo
!!! !!       !tt=tt*hxh*hyh*hzh
!!! !!       if (nproc > 1) then
!!! !!          call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
!!! !!               MPI_SUM,MPI_COMM_WORLD,ierr)
!!! !!       else
!!! !!          offset=tt
!!! !!       end if
!!!
!!! !!       if (iproc==0) print *,'offset recheck',offset
!!! !!    end if
!!!
!!! !!    !here put nproc=1 
!!! !!    call plot_density('potion_corr.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3pi,&
!!! !!         at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr(ind))
!!!
!!!
!!!
!!!!!!!!     !reread file from disk
!!!!!!!!     !overwrite pot_ion with the potential previously created
!!!!!!!!     call read_potfile(at%geocode,'potion_corr.pot',n1,n2,n3,n1i,n2i,n3i,n3pi,i3s,1,potion_corr)
!!!!!!!!
!!!!!!!!     !calculate the offset
!!!!!!!!     tt=0.d0
!!!!!!!!     do ind=1,n3pi*n2i*n1i
!!!!!!!!        tt=tt+potion_corr(ind)
!!!!!!!!     enddo
!!!!!!!!     !tt=tt*hxh*hyh*hzh
!!!!!!!!     if (nproc > 1) then
!!!!!!!!        call MPI_ALLREDUCE(tt,offset,1,mpidtypd, &
!!!!!!!!             MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!!!!     else
!!!!!!!!        offset=tt
!!!!!!!!     end if
!!!!!!!!
!!!!!!!!     if (iproc==0) print *,'offset reread',offset
!!!!!!!!
!!!!!!!!     call plot_density(at%geocode,'potion_corr_2.pot',iproc,nproc,n1,n2,n3,n1i,n2i,n3i,&
!!!!!!!!          n3pi,at%alat1,at%alat2,at%alat3,ngatherarr,potion_corr)
!!! !!         
!!!
!!! !!    i_all=-product(shape(ngatherarr))*kind(ngatherarr)
!!! !!    deallocate(ngatherarr,stat=i_stat)
!!! !!    call memocc(i_stat,i_all,'ngatherarr',subname)
!!! !!    i_all=-product(shape(potion_corr))*kind(potion_corr)
!!! !!    deallocate(potion_corr,stat=i_stat)
!!! !!    call memocc(i_stat,i_all,'potion_corr',subname)
!!!
!!!
!!! !! end if
!!!
!!!!!!  !calculate the value of the offset to be put
!!!!!!  tt_tot=0.d0
!!!!!!  do ind=1,n1i*n2i*n3i
!!!!!!     tt_tot=tt_tot+pot_ion(ind)
!!!!!!  end do
!!!!!!  print *,'actual offset',tt_tot*hxh*hyh*hzh
!!!
!!!  !use rhopot to calculate the potential from a constant electric field along y direction
!!!  if (elecfield /= 0.0_gp) then
!!!      stop 'not yet implemented'
!!!
!!!!!!!!%%     !constant electric field allowed only for free BC
!!!!!!!!%%     if (geocode == 'P') then
!!!!!!!!%%     !if (iproc == 0) 
!!!!!!!!%%           write(*,'(1x,a)') &
!!!!!!!!%%          'The constant electric field is allowed only for Free and Surfaces BC'
!!!!!!!!%%     stop
!!!!!!!!%%     end if
!!!!!!!!%%     if (iproc == 0) write(*,'(1x,3(a,1pe10.2))') &
!!!!!!!!%%          'Constant electric field of',elecfield,' Ha/Bohr for:'
!!!!!!!!%%!or         'Parabolic confining potential: rprb=',elecfield,&
!!!!!!!!%%!           ';  v_conf(r)= 1/(2*rprb**4) * r**2'
!!!!!!!!%%
!!!!!!!!%%     !write or not electric field in a separate file
!!!!!!!!%%     efwrite=.true.
!!!!!!!!%%
!!!!!!!!%%     if (n3pi > 0) then
!!!!!!!!%%        do i3=1,n3pi
!!!!!!!!%%           !z=real(i3+i3s-1-nbl3-1,gp)*hzh
!!!!!!!!%%           do i2=1,n2i
!!!!!!!!%%              y=real(i2-nbl2-1,gp)*hyh
!!!!!!!!%%                 do i1=1,n1i
!!!!!!!!%%                    !x=real(i1-nbl1-1,gp)*hxh
!!!!!!!!%%                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!!!!!!!%%                    pot_ion(ind)=pot_ion(ind)+elecfield*y
!!!!!!!!%%!                    parabola: these two lines replace the above line 
!!!!!!!!%%!                              comment out the if case and calculate x, z
!!!!!!!!%%!                    r2=(x-rx)**2+(y-ry)**2+(z-rz)**2
!!!!!!!!%%!                    pot_ion(ind)=pot_ion(ind)+0.5_gp/(elecfield**4)*r2
!!!!!!!!%%                 end do
!!!!!!!!%%           end do
!!!!!!!!%%        end do
!!!!!!!!%%
!!!!!!!!%%        if (efwrite .and. iproc == 0) then
!!!!!!!!%%           open(unit=17,file='elecpotential_y',status='unknown')
!!!!!!!!%%           do i2=nbl2+1,n2i-nbr2-1
!!!!!!!!%%              y=real(i2-nbl2-1,gp)*hyh
!!!!!!!!%%                 write(17,*)i2,y,elecfield*y
!!!!!!!!%%           end do
!!!!!!!!%%           close(17)
!!!!!!!!%%        end if
!!!!!!!!%%
!!!!!!!!%%     end if
!!!  end if
!!!
!!!  call timing(iproc,'CrtLocPot     ','OF')
!!!
!!!END SUBROUTINE createIonicPotentialModified

!!!subroutine calculateForcesLinear(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, atoms, in, hx, hy, hz, &
!!!    comms, lin, nlpspd, proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, rho, psi, fxyz, fnoise)
!!!use module_base
!!!use module_types
!!!use Poisson_Solver
!!!use module_interfaces, exceptThisOne => calculateForcesLinear
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh
!!!real(gp), intent(in):: hx, hy, hz
!!!type(locreg_descriptors),intent(in):: Glr
!!!type(orbitals_data),intent(in):: orbs
!!!type(atoms_data),intent(in):: atoms
!!!type(input_variables),intent(in):: in
!!!type(communications_arrays),intent(in):: comms
!!!type(linearParameters),intent(inout):: lin
!!!type(nonlocal_psp_descriptors),intent(in) :: nlpspd
!!!real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
!!!integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr   !!! NOT NEEDED
!!!integer,dimension(0:nproc-1,4),intent(inout) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!type(GPU_pointers),intent(inout):: GPU
!!!integer,dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
!!!real(dp),dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
!!!real(dp),dimension(lin%as%size_pkernel),intent(in):: pkernel
!!!real(8),dimension(3,atoms%nat),intent(in):: rxyz, fion, fdisp
!!!real(8),dimension(3,atoms%nat),intent(out):: fxyz
!!!real(8),intent(out):: fnoise
!!!real(8),dimension(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1)),intent(in):: rho
!!!real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi
!!!
!!!! Local variables
!!!integer:: jproc, i_stat, i_all, iat, ierr, j
!!!real(8):: hxh, hyh, hzh, ehart_fake
!!!real(gp), dimension(:,:), allocatable :: gxyz, fxyzConf
!!!real(kind=8), dimension(:,:,:,:), allocatable :: pot
!!!character(len=*),parameter:: subname='calculateForcesSub'
!!!logical:: refill_proj
!!!
!!!  hxh=0.5d0*hx
!!!  hyh=0.5d0*hy
!!!  hzh=0.5d0*hz
!!!
!!!  if (iproc==0) then
!!!     write( *,'(1x,a)')&
!!!          '----------------------------------------------------------------- Forces Calculation'
!!!  end if
!!!
!!!  ! Selfconsistent potential is saved in rhopot, 
!!!  ! new arrays rho,pot for calculation of forces ground state electronic density
!!!
!!!  ! Potential from electronic charge density
!!!
!!!  !manipulate scatter array for avoiding the GGA shift
!!!  do jproc=0,nproc-1
!!!     !n3d=n3p
!!!     nscatterarr(jproc,1)=nscatterarr(jproc,2)
!!!     !i3xcsh=0
!!!     nscatterarr(jproc,4)=0
!!!  end do
!!!
!!!  ! The charge density has already been calculated and is in rhopot.
!!!  !!if (n3p>0) then
!!!  !!   allocate(rho(Glr%d%n1i*Glr%d%n2i*n3p*in%nspin+ndebug),stat=i_stat)
!!!  !!   call memocc(i_stat,rho,'rho',subname)
!!!  !!else
!!!  !!   allocate(rho(1+ndebug),stat=i_stat)
!!!  !!   call memocc(i_stat,rho,'rho',subname)
!!!  !!end if
!!!  !!!!call sumrho(iproc,nproc,orbs,Glr,0,hxh,hyh,hzh,psi,rho,Glr%d%n1i*Glr%d%n2i*n3p,&
!!!  !!!!        nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)
!!!  !!call sumrhoForLocalizedBasis2(iproc, nproc, lin%orbs, Glr, in, lin, coeff, phi, Glr%d%n1i*Glr%d%n2i*n3d, &
!!!  !!     rho, atoms, nscatterarr)
!!!
!!!  !calculate the total density in the case of nspin==2
!!!  !!if (in%nspin==2) then
!!!  !!   call axpy(Glr%d%n1i*Glr%d%n2i*n3p,1.0_dp,rho(1+Glr%d%n1i*Glr%d%n2i*n3p),1,rho(1),1)
!!!  !!end if
!!!  if (n3p>0) then
!!!     allocate(pot(Glr%d%n1i,Glr%d%n2i,n3p,1+ndebug),stat=i_stat)
!!!     call memocc(i_stat,pot,'pot',subname)
!!!  else
!!!     allocate(pot(1,1,1,1+ndebug),stat=i_stat)
!!!     call memocc(i_stat,pot,'pot',subname)
!!!  end if
!!!
!!!  !calculate electrostatic potential
!!!  call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rho,1,pot,1)
!!!  call H_potential(atoms%geocode,'D',iproc,nproc,&
!!!       Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,pot,pkernel,pot,ehart_fake,0.0_dp,.false.)
!!!
!!!
!!!  allocate(gxyz(3,atoms%nat+ndebug),stat=i_stat)
!!!  call memocc(i_stat,gxyz,'gxyz',subname)
!!!
!!!  call timing(iproc,'Forces        ','ON')
!!!  ! calculate local part of the forces gxyz
!!!  call local_forces(iproc,atoms,rxyz,hxh,hyh,hzh,&
!!!       Glr%d%n1,Glr%d%n2,Glr%d%n3,n3p,i3s+i3xcsh,Glr%d%n1i,Glr%d%n2i,rho,pot,gxyz)
!!!
!!!  !i_all=-product(shape(rho))*kind(rho)
!!!  !deallocate(rho,stat=i_stat)
!!!  !call memocc(i_stat,i_all,'rho',subname)
!!!  i_all=-product(shape(pot))*kind(pot)
!!!  deallocate(pot,stat=i_stat)
!!!  call memocc(i_stat,i_all,'pot',subname)
!!!
!!!  if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)',advance='no')'Calculate nonlocal forces...'
!!!
!!!  !refill projectors for tails, davidson
!!!  !refill_proj=(in%calc_tail .or. DoDavidson) .and. DoLastRunThings
!!!  refill_proj=.false.  !! IS THIS CORRECT??
!!!
!!!  call nonlocal_forces(iproc,Glr,hx,hy,hz,atoms,rxyz,&
!!!       orbs,nlpspd,proj,Glr%wfd,psi,gxyz,refill_proj)
!!!
!!!  if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)')'done.'
!!!
!!!  ! Add up all the force contributions
!!!  if (nproc > 1) then
!!!     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!  else
!!!     do iat=1,atoms%nat
!!!        fxyz(1,iat)=gxyz(1,iat)
!!!        fxyz(2,iat)=gxyz(2,iat)
!!!        fxyz(3,iat)=gxyz(3,iat)
!!!     enddo
!!!  end if
!!!
!!!  !!$  if (iproc == 0) then
!!!  !!$     sumx=0.d0 ; sumy=0.d0 ; sumz=0.d0
!!!  !!$     fumx=0.d0 ; fumy=0.d0 ; fumz=0.d0
!!!  !!$     do iat=1,atoms%nat
!!!  !!$        sumx=sumx+fxyz(1,iat) ; sumy=sumy+fxyz(2,iat) ; sumz=sumz+fxyz(3,iat)
!!!  !!$        fumx=fumx+fion(1,iat) ; fumy=fumy+fion(2,iat) ; fumz=fumz+fion(3,iat)
!!!  !!$     enddo
!!!  !!$     write(77,'(a30,3(1x,e10.3))') 'translat. force total pot ',sumx,sumy,sumz
!!!  !!$     write(77,'(a30,3(1x,e10.3))') 'translat. force ionic pot ',fumx,fumy,fumz
!!!  !!$  endif
!!!
!!!    !add to the forces the ionic and dispersion contribution 
!!!    do iat=1,atoms%nat
!!!       fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)+fdisp(1,iat)
!!!       fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)+fdisp(2,iat)
!!!       fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)+fdisp(3,iat)
!!!    enddo
!!!
!!!!!! TEST !!!
!!!  call timing(iproc,'Forces        ','OF')
!!!  !!call confinementCorrection()
!!!  !!fxyz=fxyz+fxyzConf
!!!  !call pulayCorrection()
!!!  call timing(iproc,'Forces        ','ON')
!!!!!!!!!!!!!!!
!!!
!!!    !i_all=-product(shape(fion))*kind(fion)
!!!    !deallocate(fion,stat=i_stat)
!!!    !call memocc(i_stat,i_all,'fion',subname)
!!!    !i_all=-product(shape(fdisp))*kind(fdisp)
!!!    !deallocate(fdisp,stat=i_stat)
!!!    !call memocc(i_stat,i_all,'fdisp',subname)
!!!    i_all=-product(shape(gxyz))*kind(gxyz)
!!!    deallocate(gxyz,stat=i_stat)
!!!    call memocc(i_stat,i_all,'gxyz',subname)
!!!
!!!
!!!  !!do iat=1,atoms%nat
!!!  !!   if(iproc==0) write(*,'(a,i0,3es14.5)') 'forces for atom ',iat, fxyz(1,iat), fxyz(2,iat), fxyz(3,iat)
!!!  !!end do
!!!
!!!  !subtraction of zero of the forces, disabled for the moment
!!!  !the zero of the forces depends on the atomic positions
!!!  !if (in%gaussian_help .and. .false.) then
!!!  call clean_forces(iproc,atoms,rxyz,fxyz,fnoise)
!!!  !end if
!!!
!!!  if(iproc==0) then
!!!      write(*,'(1x,a)') 'Force values for all atoms in x, y, z direction.'
!!!      do iat=1,atoms%nat
!!!         write(*,'(3x,i0,1x,a6,1x,3(1x,es17.10))') &
!!!              iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
!!!      end do
!!!  end if
!!!
!!!  call timing(iproc,'Forces        ','OF')
!!!
!!!
!!!end subroutine calculateForcesLinear
