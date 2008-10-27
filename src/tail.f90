!!****f* BigDFT/CalculateTailCorrection
!! NAME
!!   CalculateTailCorrection
!!
!! FUNCTION
!!  Calculate the finite size corrections over wavefunctions
!!  Conceived only for isolated Boundary Conditions
!!  
!! COPYRIGHT
!!    Copyright (C) 2007 UNIBAS,CEA
!!
!! SOURCE
!!
subroutine CalculateTailCorrection(iproc,nproc,at,n1,n2,n3,rbuf,norb,norbp,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,eval,&
     pot,hgrid,rxyz,radii_cf,crmult,frmult,cpmult,fpmult,nspin,spinsgn,&
     proj,psi,occup,output_grid,ekin_sum,epot_sum,eproj_sum)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
  logical, intent(in) :: output_grid
  integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,ncongt,nspin
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf,cpmult,fpmult
  real(kind=8), dimension(norb), intent(in) :: occup,eval,spinsgn
  real(kind=8), dimension(at%ntypes,3), intent(in) :: radii_cf
  real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
  real(kind=8), dimension(2*n1+31,2*n2+31,2*n3+31,nspin), intent(in) :: pot
  real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
  real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
  !local variables
  character(len=*), parameter :: subname='CalculateTailCorrection'
  integer :: iseg,i0,j0,i1,j1,i2,i3,ii,iat,iorb,npt,ipt,i,ierr,i_all,i_stat,nbuf,ispin
  integer :: nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f
  real(kind=8) :: alatb1,alatb2,alatb3,ekin,epot,eproj,tt,cprecr,sum_tail,ekin1,epot1,eproj1
  type(wavefunctions_descriptors) :: wfdb
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f
  integer, dimension(:,:,:), allocatable :: ibbyz_c,ibbyz_f,ibbxz_c,ibbxz_f,ibbxy_c,ibbxy_f
  real(kind=8), dimension(:,:), allocatable :: txyz,wrkallred
  real(kind=8), dimension(:), allocatable :: psib,hpsib,psir
  integer, dimension(:), pointer :: keybv
  integer, dimension(:,:), pointer :: keybg

  !for shrink:
  integer, allocatable, dimension(:,:,:) :: ibbzzx_c,ibbyyzz_c
  integer, allocatable, dimension(:,:,:) :: ibbxy_ff,ibbzzx_f,ibbyyzz_f

  !for grow:
  integer, allocatable, dimension(:,:,:) :: ibbzxx_c,ibbxxyy_c
  integer, allocatable, dimension(:,:,:) :: ibbyz_ff,ibbzxx_f,ibbxxyy_f

  !real space border:
  integer, allocatable, dimension(:,:,:) :: ibbyyzz_r 

  integer nw1,nw2

  real(kind=8), dimension(:,:,:), allocatable::x_c!input 
  real(kind=8), dimension(:,:,:,:), allocatable :: x_f ! input
  real(kind=8), dimension(:,:,:), allocatable :: x_f1,x_f2,x_f3 ! internal
  real(kind=8), dimension(:), allocatable :: w1,w2
  real(kind=8), dimension(:,:,:), allocatable::y_c!output 
  real(kind=8), dimension(:,:,:,:), allocatable :: y_f! output

  nbuf=nint(rbuf/hgrid)
  !    --- new grid sizes n1,n2,n3
  nb1=n1+2*nbuf
  nb2=n2+2*nbuf
  nb3=n3+2*nbuf
  alatb1=real(nb1,kind=8)*hgrid 
  alatb2=real(nb2,kind=8)*hgrid 
  alatb3=real(nb3,kind=8)*hgrid

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '---------------------------------------------- Estimation of Finite-Size Corrections'
     write(*,'(1x,a,f6.3,a)') &
          'F-S Correction for an effective space of ',rbuf,' AU more around each external atom'
     write(*,'(1x,a,i6,a)') &
          '                  requires the adding of ',nbuf,' additional grid points around cell'
     write(*,'(1x,a,19x,a)') &
          '   Effective box size,   Atomic Units:','grid spacing units:'
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '            ',alatb1,alatb2,alatb3,nb1,nb2,nb3
  end if


  !---reformat keyg_p
  do iseg=1,nlpspd%nseg_p(2*at%nat)
     j0=nlpspd%keyg_p(1,iseg)
     j1=nlpspd%keyg_p(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     i3=i3+nbuf
     i2=i2+nbuf
     i1=i1+nbuf
     i0=i0+nbuf
     j0=i3*((nb1+1)*(nb2+1)) + i2*(nb1+1) + i0+1
     j1=i3*((nb1+1)*(nb2+1)) + i2*(nb1+1) + i1+1
     nlpspd%keyg_p(1,iseg)=j0
     nlpspd%keyg_p(2,iseg)=j1
  end do

  !---reformat wavefunctions

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nbfl1=nfl1+nbuf ; nbfl2=nfl2+nbuf ; nbfl3=nfl3+nbuf
  nbfu1=nfu1+nbuf ; nbfu2=nfu2+nbuf ; nbfu3=nfu3+nbuf
  if (iproc.eq.0) then
     write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
          '  Extremes for the new high resolution grid points:',&
          nbfl1,'<',nbfu1,nbfl2,'<',nbfu2,nbfl3,'<',nbfu3
  endif

  ! change atom coordinates according to the enlarged box
  allocate(txyz(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,txyz,'txyz',subname)
  do iat=1,at%nat
     txyz(1,iat)=rxyz(1,iat)+real(nbuf,kind=8)*hgrid
     txyz(2,iat)=rxyz(2,iat)+real(nbuf,kind=8)*hgrid
     txyz(3,iat)=rxyz(3,iat)+real(nbuf,kind=8)*hgrid
  enddo

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:nb1,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_c,'logrid_c',subname)
  allocate(logrid_f(0:nb1,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_f,'logrid_f',subname)
  allocate(ibbyz_c(2,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,ibbyz_c,'ibbyz_c',subname)
  allocate(ibbxz_c(2,0:nb1,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxz_c,'ibbxz_c',subname)
  allocate(ibbxy_c(2,0:nb1,0:nb2+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxy_c,'ibbxy_c',subname)
  allocate(ibbyz_f(2,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,ibbyz_f,'ibbyz_f',subname)
  allocate(ibbxz_f(2,0:nb1,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxz_f,'ibbxz_f',subname)
  allocate(ibbxy_f(2,0:nb1,0:nb2+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxy_f,'ibbxy_f',subname)

  !!*********************************Alexey*********************************************************
  !   allocate for grow
  allocate(ibbzxx_c(2,0:nb3,-14:2*nb1+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbzxx_c,'ibbzxx_c',subname)
  allocate(ibbxxyy_c(2,-14:2*nb1+16,-14:2*nb2+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxxyy_c,'ibbxxyy_c',subname)
  allocate(ibbyz_ff(2,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,ibbyz_ff,'ibbyz_ff',subname)
  allocate(ibbzxx_f(2,0:nb3,-14:2*nb1+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbzxx_f,'ibbzxx_f',subname)
  allocate(ibbxxyy_f(2,-14:2*nb1+16,-14:2*nb2+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxxyy_f,'ibbxxyy_f',subname)

  !allocate for shrink
  allocate(ibbzzx_c(2,-14:2*nb3+16,0:nb1+ndebug),stat=i_stat)
  call memocc(i_stat,ibbzzx_c,'ibbzzx_c',subname)
  allocate(ibbyyzz_c(2,-14:2*nb2+16,-14:2*nb3+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbyyzz_c,'ibbyyzz_c',subname)
  allocate(ibbxy_ff(2,0:nb1,0:nb2+ndebug),stat=i_stat)
  call memocc(i_stat,ibbxy_ff,'ibbxy_ff',subname)
  allocate(ibbzzx_f(2,-14:2*nb3+16,0:nb1+ndebug),stat=i_stat)
  call memocc(i_stat,ibbzzx_f,'ibbzzx_f',subname)
  allocate(ibbyyzz_f(2,-14:2*nb2+16,-14:2*nb3+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbyyzz_f,'ibbyyzz_f',subname)

  !allocate for real space
  allocate(ibbyyzz_r(2,-14:2*nb2+16,-14:2*nb3+16+ndebug),stat=i_stat)
  call memocc(i_stat,ibbyyzz_r,'ibbyyzz_r',subname)
  !***********************************************************************************************
  ! coarse grid quantities
  call fill_logrid('F',nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,nbuf,at%nat,at%ntypes,at%iatype,txyz, & 
       radii_cf(1,1),crmult,hgrid,hgrid,hgrid,logrid_c)
  call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,nvctrb_c)

  if (iproc.eq.0) then
     write(*,'(2(1x,a,i10))') &
          'Coarse resolution grid: Number of segments= ',nsegb_c,'points=',nvctrb_c
     !write(*,'(1x,a,2(1x,i10))') 'BIG: orbitals have coarse segment, elements',nsegb_c,nvctrb_c
  end if
  call make_bounds(nb1,nb2,nb3,logrid_c,ibbyz_c,ibbxz_c,ibbxy_c)

  ! fine grid quantities
  call fill_logrid('F',nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,0,at%nat,at%ntypes,at%iatype,txyz, & 
       radii_cf(1,2),frmult,hgrid,hgrid,hgrid,logrid_f)
  call num_segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,nvctrb_f)
  if (iproc.eq.0) then
     write(*,'(2(1x,a,i10))') &
          '  Fine resolution grid: Number of segments= ',nsegb_f,'points=',nvctrb_f
     !write(*,'(1x,a,2(1x,i10))') 'BIG: orbitals have fine   segment, elements',nsegb_f,7*nvctrb_f
  end if
  call make_bounds(nb1,nb2,nb3,logrid_f,ibbyz_f,ibbxz_f,ibbxy_f)

! Create the file grid.xyz to visualize the grid of functions
  if (iproc ==0 .and. output_grid) then
     write(*,'(1x,a)')&
          'Writing the file describing the new atomic positions of the effective system'
     open(unit=22,file='grid_tail.xyz',status='unknown')
     write(22,*) nvctrb_c+nvctrb_f,' atomic' 
     write(22,*)'complete simulation grid for the tail correction'
     do iat=1,at%nat
        write(22,'(a6,2x,3(1x,e12.5),3x)') &
             trim(at%atomnames(at%iatype(iat))),txyz(1,iat),txyz(2,iat),txyz(3,iat)
     enddo
     do i3=0,nb3  
        do i2=0,nb2  
           do i1=0,nb1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  g ',real(i1,kind=8)*hgrid,real(i2,kind=8)*hgrid,real(i3,kind=8)*hgrid
           enddo
        enddo
     end do
     do i3=0,nb3 
        do i2=0,nb2 
           do i1=0,nb1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  G ',real(i1,kind=8)*hgrid,real(i2,kind=8)*hgrid,real(i3,kind=8)*hgrid
           enddo
        enddo
     enddo
     close(22)
  endif

  call make_all_ib(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,&
       ibbxy_c,ibbzzx_c,ibbyyzz_c,ibbxy_f,ibbxy_ff,ibbzzx_f,ibbyyzz_f,&
       ibbyz_c,ibbzxx_c,ibbxxyy_c,ibbyz_f,ibbyz_ff,ibbzxx_f,ibbxxyy_f,ibbyyzz_r)

  ! now fill the wavefunction descriptor arrays
  allocate(keybg(2,nsegb_c+nsegb_f+ndebug),stat=i_stat)
  call memocc(i_stat,keybg,'keybg',subname)
  allocate(keybv(nsegb_c+nsegb_f+ndebug),stat=i_stat)
  call memocc(i_stat,keybv,'keybv',subname)
  ! coarse grid quantities
  call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,keybg(1,1),keybv(1))

  ! fine grid quantities
  call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,keybg(1,nsegb_c+1),keybv(nsegb_c+1))

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c',subname)
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f',subname)


  !assign the values of the big wavefunction descriptors (used for on-the-fly projectors calc)
  if (DistProjApply) then
     wfdb%nvctr_c=nvctrb_c
     wfdb%nvctr_f=nvctrb_f
     wfdb%nseg_c=nsegb_c
     wfdb%nseg_f=nsegb_f
     wfdb%keyv => keybv
     wfdb%keyg => keybg
  end if


  ! allocations for arrays holding the wavefunction
  !if (iproc.eq.0) &
  !     write(*,'(1x,a,i0)') 'Allocate words for psib and hpsib ',2*(nvctrb_c+7*nvctrb_f)
  allocate(psib(nvctrb_c+7*nvctrb_f+ndebug),stat=i_stat)
  call memocc(i_stat,psib,'psib',subname)
  allocate(hpsib(nvctrb_c+7*nvctrb_f+ndebug),stat=i_stat)
  call memocc(i_stat,hpsib,'hpsib',subname)
  !if (iproc.eq.0) write(*,*) 'Allocation done'

  ! work arrays applylocpotkin
  allocate(psir((2*nb1+31)*(2*nb2+31)*(2*nb3+31)+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  if (iproc.eq.0) then
     write(*,'(1x,a,i0)') &
          'Wavefunction memory occupation in the extended grid (Bytes): ',&
          (nvctrb_c+7*nvctrb_f)*8
     write(*,'(1x,a,i0,a)') &
          'Calculating tail corrections on ',norbp,' orbitals per processor.'
     write(*,'(1x,a)',advance='no') &
          '     orbitals are processed separately'
  end if

  nw1=max(2*(nb3+1)*(2*nb1+31)*(2*nb2+31),&   ! shrink convention: nw1>nw2
       2*(nb1+1)*(2*nb2+31)*(2*nb3+31))
  nw2=max(4*(nb2+1)*(nb3+1)*(2*nb1+31),&
       4*(nb1+1)*(nb2+1)*(2*nb3+31))

  allocate(x_c(0:nb1,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,x_c,'x_c',subname)
  allocate(y_c(0:nb1,0:nb2,0:nb3+ndebug),stat=i_stat)
  call memocc(i_stat,y_c,'y_c',subname)
  allocate(x_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f,'x_f',subname)
  allocate(y_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3+ndebug),stat=i_stat)
  call memocc(i_stat,y_f,'y_f',subname)
  allocate(w1(nw1+ndebug),stat=i_stat)
  call memocc(i_stat,w1,'w1',subname)
  allocate(w2(nw2+ndebug),stat=i_stat)
  call memocc(i_stat,w2,'w2',subname)

  allocate(x_f1(nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f1,'x_f1',subname)
  allocate(x_f2(nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f2,'x_f2',subname)
  allocate(x_f3(nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f3,'x_f3',subname)
  !put to zero the arrays for the hamiltonian procedure
  call razero((nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f1)
  call razero((nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f2)
  call razero((nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f3)
  call razero((nb1+1)*(nb2+1)*(nb3+1),x_c)
  call razero(7*(nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),x_f)
  call razero((nb1+1)*(nb2+1)*(nb3+1),y_c)
  call razero(7*(nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),y_f)
  call razero((2*nb1+31)*(2*nb2+31)*(2*nb3+31),psir)
  ekin_sum=0.d0
  epot_sum=0.d0
  eproj_sum=0.d0

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

     !build the compressed wavefunction in the enlarged box
     call transform_fortail(n1,n2,n3,nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,&
          wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),&
          wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),  &
          nsegb_c,nvctrb_c,keybg(1,1),keybv(1),nsegb_f,nvctrb_f,&
          keybg(1,nsegb_c+1),keybv(nsegb_c+1),&
          nbuf,psi(1,iorb-iproc*norbp),psi(wfd%nvctr_c+1,iorb-iproc*norbp),  & 
          x_c,x_f,psib(1),psib(nvctrb_c+1))

     !write(*,*) 'transform_fortail finished',iproc,iorb


     npt=2
     tail_adding: do ipt=1,npt

        if(spinsgn(iorb)>0.0d0) then
           ispin=1
        else
           ispin=2
        end if

        !calculate gradient
        call applylocpotkinone(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,nbuf, &
             hgrid,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f,keybg,keybv,  &
             ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,y_c,y_f,psir,  &
             psib,pot(1,1,1,ispin),hpsib,epot,ekin, &
             x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
             ibbzzx_c,ibbyyzz_c,ibbxy_ff,ibbzzx_f,ibbyyzz_f,&
             ibbzxx_c,ibbxxyy_c,ibbyz_ff,ibbzxx_f,ibbxxyy_f,nw1,nw2,ibbyyzz_r,1)
        !write(*,'(a,3i3,2f12.8)') 'applylocpotkinone finished',iproc,iorb,ipt,epot,ekin

        if (DistProjApply) then
           !here the box do not change, only the projectors descriptors
           call applyprojectorsonthefly(0,1,1,1,(/1.0_gp/),at,n1,n2,n3,&
                rxyz,hgrid,hgrid,hgrid,cpmult,fpmult,radii_cf,wfdb,nlpspd,proj,psib,hpsib,eproj)
           !only the wavefunction descriptors must change
        else
           call applyprojectorsone(at%ntypes,at%nat,at%iatype,at%psppar,at%npspcode, &
                nlpspd%nprojel,nlpspd%nproj,nlpspd%nseg_p,&
                nlpspd%keyg_p,nlpspd%keyv_p,nlpspd%nvctr_p,proj,  &
                nsegb_c,nsegb_f,keybg,keybv,nvctrb_c,nvctrb_f,  & 
                psib,hpsib,eproj)
           !write(*,'(a,2i3,2f12.8)') 'applyprojectorsone finished',iproc,iorb,eproj,sum_tail
        end if


        !calculate residue for the single orbital
        tt=0.d0
        do i=1,nvctrb_c+7*nvctrb_f
           hpsib(i)=hpsib(i)-eval(iorb)*psib(i)
           tt=tt+hpsib(i)**2
        enddo
        tt=sqrt(tt)

        if (ipt.eq.npt) exit tail_adding

        !calculate tail using the preconditioner as solver for the green function application
        cprecr=-eval(iorb)
        call precong(iorb,nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3, &
             nsegb_c,nvctrb_c,nsegb_f,nvctrb_f,keybg,keybv, &
             ncongt,cprecr,hgrid,ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,hpsib)
        !call plot_wf(10,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keybg,keybv,nsegb_f,nvctrb_f,  & 
        !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

        ! add tail to the bulk wavefunction
        sum_tail=0.d0
        do i=1,nvctrb_c+7*nvctrb_f
           psib(i)=psib(i)-hpsib(i)
           sum_tail=sum_tail+psib(i)**2
        enddo
        sum_tail=sqrt(sum_tail)
        !write(*,'(1x,a,i3,3(1x,1pe13.6),1x,1pe9.2)') &
        !     'BIG: iorb,ekin,epot,eproj,gnrm',iorb,ekin,epot,eproj,tt

        !values of the energies before tail application
        ekin1=ekin
        epot1=epot
        eproj1=eproj
        !write(*,'(1x,a,1x,i0,f18.14)') 'norm orbital + tail',iorb,sum_tail
        !call plot_wf(20,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keybg,keybv,nsegb_f,nvctrb_f,  & 
        !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

        sum_tail=1.d0/sum_tail
        do i=1,nvctrb_c+7*nvctrb_f
           psib(i)=psib(i)*sum_tail
        enddo

     end do tail_adding

!!$     write(*,'(1x,a,i3,3(1x,1pe13.6),2(1x,1pe9.2))') &
!!$          'BIG: iorb,denergies,gnrm,dnorm',&
!!$          iorb,ekin-ekin1,epot-epot1,eproj-eproj1,tt,sum_tail-1.d0

     if (iproc == 0) then
        write(*,'(a)',advance='no') &
             repeat('.',(iorb*40)/norbp-((iorb-1)*40)/norbp)
     end if
     ekin_sum=ekin_sum+ekin*occup(iorb)
     epot_sum=epot_sum+epot*occup(iorb)
     eproj_sum=eproj_sum+eproj*occup(iorb)

  end do

  if (iproc == 0) then
     write(*,'(1x,a)')'done.'
  end if

  i_all=-product(shape(txyz))*kind(txyz)
  deallocate(txyz,stat=i_stat)
  call memocc(i_stat,i_all,'txyz',subname)
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  i_all=-product(shape(psib))*kind(psib)
  deallocate(psib,stat=i_stat)
  call memocc(i_stat,i_all,'psib',subname)
  i_all=-product(shape(hpsib))*kind(hpsib)
  deallocate(hpsib,stat=i_stat)
  call memocc(i_stat,i_all,'hpsib',subname)

  if (DistProjApply) then
     call deallocate_wfd(wfdb,subname)
  else
     i_all=-product(shape(keybg))*kind(keybg)
     deallocate(keybg,stat=i_stat)
     call memocc(i_stat,i_all,'keybg',subname)
     i_all=-product(shape(keybv))*kind(keybv)
     deallocate(keybv,stat=i_stat)
     call memocc(i_stat,i_all,'keybv',subname)
  end if

  i_all=-product(shape(ibbyz_c))*kind(ibbyz_c)
  deallocate(ibbyz_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbyz_c',subname)
  i_all=-product(shape(ibbxz_c))*kind(ibbxz_c)
  deallocate(ibbxz_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxz_c',subname)
  i_all=-product(shape(ibbxy_c))*kind(ibbxy_c)
  deallocate(ibbxy_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxy_c',subname)
  i_all=-product(shape(ibbyz_f))*kind(ibbyz_f)
  deallocate(ibbyz_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbyz_f',subname)
  i_all=-product(shape(ibbxz_f))*kind(ibbxz_f)
  deallocate(ibbxz_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxz_f',subname)
  i_all=-product(shape(ibbxy_f))*kind(ibbxy_f)
  deallocate(ibbxy_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxy_f',subname)

  !!*****************************Alexey************************************************************
  i_all=-product(shape(x_c))*kind(x_c)
  deallocate(x_c,stat=i_stat)
  call memocc(i_stat,i_all,'x_c',subname)

  i_all=-product(shape(y_c))*kind(y_c)
  deallocate(y_c,stat=i_stat)
  call memocc(i_stat,i_all,'y_c',subname)

  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1',subname)

  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2',subname)
  
  i_all=-product(shape(x_f1))*kind(x_f1)
  deallocate(x_f1,stat=i_stat)
  call memocc(i_stat,i_all,'x_f1',subname)
  i_all=-product(shape(x_f2))*kind(x_f2)
  deallocate(x_f2,stat=i_stat)
  call memocc(i_stat,i_all,'x_f2',subname)
  i_all=-product(shape(x_f3))*kind(x_f3)
  deallocate(x_f3,stat=i_stat)
  call memocc(i_stat,i_all,'x_f3',subname)
  
  i_all=-product(shape(x_f))*kind(x_f)
  deallocate(x_f,stat=i_stat)
  call memocc(i_stat,i_all,'x_f',subname)

  i_all=-product(shape(y_f))*kind(y_f)
  deallocate(y_f,stat=i_stat)
  call memocc(i_stat,i_all,'y_f',subname)

  i_all=-product(shape(ibbzzx_c))*kind(ibbzzx_c)
  deallocate(ibbzzx_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbzzx_c',subname)

  i_all=-product(shape(ibbyyzz_c))*kind(ibbyyzz_c)
  deallocate(ibbyyzz_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbyyzz_c',subname)

  i_all=-product(shape(ibbxy_ff))*kind(ibbxy_ff)
  deallocate(ibbxy_ff,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxy_ff',subname)

  i_all=-product(shape(ibbzzx_f))*kind(ibbzzx_f)
  deallocate(ibbzzx_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbzzx_f',subname)

  i_all=-product(shape(ibbyyzz_f))*kind(ibbyyzz_f)
  deallocate(ibbyyzz_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbyyzz_f',subname)

  i_all=-product(shape(ibbzxx_c))*kind(ibbzxx_c)
  deallocate(ibbzxx_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbzxx_c',subname)

  i_all=-product(shape(ibbxxyy_c))*kind(ibbxxyy_c)
  deallocate(ibbxxyy_c,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxxyy_c',subname)

  i_all=-product(shape(ibbyz_ff))*kind(ibbyz_ff)
  deallocate(ibbyz_ff,stat=i_stat)
  call memocc(i_stat,i_all,'ibbyz_ff',subname)

  i_all=-product(shape(ibbzxx_f))*kind(ibbzxx_f)
  deallocate(ibbzxx_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbzxx_f',subname)

  i_all=-product(shape(ibbxxyy_f))*kind(ibbxxyy_f)
  deallocate(ibbxxyy_f,stat=i_stat)
  call memocc(i_stat,i_all,'ibbxxyy_f',subname)

  i_all=-product(shape(ibbyyzz_r))*kind(ibbyyzz_r)
  deallocate(ibbyyzz_r,stat=i_stat)
  call memocc(i_stat,i_all,'ibbyyzz_r',subname)
  !!***********************************************************************************************

  if (nproc > 1) then
     !if (iproc.eq.0) then
     !   write(*,'(1x,a,f27.14)')'Tail calculation ended'
     !endif
     allocate(wrkallred(3,2+ndebug),stat=i_stat)
     call memocc(i_stat,wrkallred,'wrkallred',subname)
     wrkallred(1,2)=ekin_sum
     wrkallred(2,2)=epot_sum 
     wrkallred(3,2)=eproj_sum 
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     ekin_sum=wrkallred(1,1) 
     epot_sum=wrkallred(2,1) 
     eproj_sum=wrkallred(3,1)
     i_all=-product(shape(wrkallred))*kind(wrkallred)
     deallocate(wrkallred,stat=i_stat)
     call memocc(i_stat,i_all,'wrkallred',subname)
  endif

end subroutine CalculateTailCorrection
!!***


subroutine transform_fortail(n1,n2,n3,nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,& 
     mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     msegb_c,mvctrb_c,keybg_c,keybv_c,msegb_f,mvctrb_f,keybg_f,keybv_f,  & 
     nbuf,psi_c,psi_f,psig_c,psig_f,psib_c,psib_f)
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension keybg_c(2,msegb_c),keybv_c(msegb_c),keybg_f(2,msegb_f),keybv_f(msegb_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
  dimension psib_c(mvctrb_c),psib_f(7,mvctrb_f)
  dimension psig_c(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf)
  dimension psig_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3)

  call razero((n1+1+2*nbuf)*(n2+1+2*nbuf)*(n3+1+2*nbuf),psig_c)
  call razero(7*(nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),psig_f)

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i+nbuf,i2+nbuf,i3+nbuf)=psi_c(i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_c
     jj=keybv_c(iseg)
     j0=keybg_c(1,iseg)
     j1=keybg_c(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_c(i-i0+jj)=psig_c(i,i2,i3)
     enddo
  enddo


  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(1,i-i0+jj)
        psig_f(2,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(2,i-i0+jj)
        psig_f(3,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(3,i-i0+jj)
        psig_f(4,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(4,i-i0+jj)
        psig_f(5,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(5,i-i0+jj)
        psig_f(6,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(6,i-i0+jj)
        psig_f(7,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(7,i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_f
     jj=keybv_f(iseg)
     j0=keybg_f(1,iseg)
     j1=keybg_f(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
        psib_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
        psib_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
        psib_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
        psib_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
        psib_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
        psib_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
     enddo
  enddo

end subroutine transform_fortail

subroutine transform_fortail_prev(n1,n2,n3,nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,& 
     mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     msegb_c,mvctrb_c,keybg_c,keybv_c,msegb_f,mvctrb_f,keybg_f,keybv_f,  & 
     nbuf,psi_c,psi_f,psig_c,psig_fc,psig_f,psib_c,psib_f)
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension keybg_c(2,msegb_c),keybv_c(msegb_c),keybg_f(2,msegb_f),keybv_f(msegb_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
  dimension psib_c(mvctrb_c),psib_f(7,mvctrb_f)
  dimension psig_c(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf)
  dimension psig_fc(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf,3),  & 
       psig_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3)

  call razero((n1+1+2*nbuf)*(n2+1+2*nbuf)*(n3+1+2*nbuf),psig_c)
  call razero(3*(n1+1+2*nbuf)*(n2+1+2*nbuf)*(n3+1+2*nbuf),psig_fc)
  call razero(7*(nbfu1-nbfl1+1)*(nbfu2-nbfl2+1)*(nbfu3-nbfl3+1),psig_f)

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i+nbuf,i2+nbuf,i3+nbuf)=psi_c(i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_c
     jj=keybv_c(iseg)
     j0=keybg_c(1,iseg)
     j1=keybg_c(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_c(i-i0+jj)=psig_c(i,i2,i3)
     enddo
  enddo


  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(1,i-i0+jj)
        psig_fc(i+nbuf,i2+nbuf,i3+nbuf,1)=psig_f(1,i+nbuf,i2+nbuf,i3+nbuf)
        psig_f(2,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(2,i-i0+jj)
        psig_fc(i+nbuf,i2+nbuf,i3+nbuf,2)=psig_f(2,i+nbuf,i2+nbuf,i3+nbuf)
        psig_f(3,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(3,i-i0+jj)
        psig_f(4,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(4,i-i0+jj)
        psig_fc(i+nbuf,i2+nbuf,i3+nbuf,3)=psig_f(4,i+nbuf,i2+nbuf,i3+nbuf)
        psig_f(5,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(5,i-i0+jj)
        psig_f(6,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(6,i-i0+jj)
        psig_f(7,i+nbuf,i2+nbuf,i3+nbuf)=psi_f(7,i-i0+jj)
     enddo
  enddo

  do iseg=1,msegb_f
     jj=keybv_f(iseg)
     j0=keybg_f(1,iseg)
     j1=keybg_f(2,iseg)
     ii=j0-1
     i3=ii/((nb1+1)*(nb2+1))
     ii=ii-i3*(nb1+1)*(nb2+1)
     i2=ii/(nb1+1)
     i0=ii-i2*(nb1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psib_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
        psib_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
        psib_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
        psib_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
        psib_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
        psib_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
        psib_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
     enddo
  enddo

end subroutine transform_fortail_prev
