!> @file
!!  Routines to do tail calculation (correct effect of finite size)
!! @author
!!    Copyright (C) 2007-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate the finite size corrections over wavefunctions
!! Conceived only for isolated Boundary Conditions, no SIC correction
subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
     Glr,nlpspd,ncongt,pot,hgrid,rxyz,radii_cf,crmult,frmult,nspin,&
     proj,psi,output_denspot,ekin_sum,epot_sum,eproj_sum)
  use module_base
  use module_types
  use module_interfaces, except_this_one => CalculateTailCorrection
  implicit none
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
  integer, intent(in) :: iproc,nproc,ncongt,nspin
  logical, intent(in) :: output_denspot
  real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
  real(kind=8), dimension(at%ntypes,3), intent(in) :: radii_cf
  real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
  real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nspin), intent(in) :: pot
  real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
  real(kind=8), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
  real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
  !local variables
  type(locreg_descriptors) :: lr
  character(len=*), parameter :: subname='CalculateTailCorrection'
  integer :: iseg,i0,j0,i1,j1,i2,i3,ii,iat,iorb,npt,ipt,i,ierr,i_all,i_stat,nbuf,ispin
  integer :: nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3
  integer :: n1,n2,n3,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f
  real(kind=8) :: alatb1,alatb2,alatb3,ekin,epot,eproj,tt,cprecr,sum_tail !n(c) eproj1 epot1,ekin1
  type(orbitals_data) :: orbsb
  type(wavefunctions_descriptors) :: wfdb
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f
  integer, dimension(:,:,:), allocatable :: ibbyz_c,ibbyz_f,ibbxz_c,ibbxz_f,ibbxy_c,ibbxy_f
  real(kind=8), dimension(:,:), allocatable :: txyz,wrkallred
  real(kind=8), dimension(:), allocatable :: psib,hpsib,psir
  integer, dimension(:), pointer :: keyv
  integer, dimension(:,:), pointer :: keyg

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

  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3

  nbuf=nint(rbuf/hgrid)
  !    --- new grid sizes n1,n2,n3
  nb1=n1+2*nbuf
  nb2=n2+2*nbuf
  nb3=n3+2*nbuf

  ! Create new structure with modified grid sizes
  call create_Glr(Glr%geocode,nb1,nb2,nb3,Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3,Glr%d%nfu1,&
             Glr%d%nfu2,Glr%d%nfu3,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,Glr%wfd,Glr%bounds,lr)
 
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

  do iat=1,at%nat
     do iseg=1,nlpspd%plr(iat)%wfd%nseg_c+nlpspd%plr(iat)%wfd%nseg_f
        j0=nlpspd%plr(iat)%wfd%keyglob(1,iseg)
        j1=nlpspd%plr(iat)%wfd%keyglob(2,iseg)
        !do iseg=1,nlpspd%nseg_p(2*at%nat)
        !j0=nlpspd%keyg_p(1,iseg)
        !j1=nlpspd%keyg_p(2,iseg)
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
        nlpspd%plr(iat)%wfd%keyglob(1,iseg)=j0
        nlpspd%plr(iat)%wfd%keyglob(2,iseg)=j1
!!$        nlpspd%keyg_p(1,iseg)=j0
!!$        nlpspd%keyg_p(2,iseg)=j1
     end do
  end do
!end do

  !---reformat wavefunctions

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nbfl1=Glr%d%nfl1+nbuf ; nbfl2=Glr%d%nfl2+nbuf ; nbfl3=Glr%d%nfl3+nbuf
  nbfu1=Glr%d%nfu1+nbuf ; nbfu2=Glr%d%nfu2+nbuf ; nbfu3=Glr%d%nfu3+nbuf
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
  if (iproc ==0 .and. output_denspot) then
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
  allocate(keyg(2,nsegb_c+nsegb_f+ndebug),stat=i_stat)
  call memocc(i_stat,keyg,'keyg',subname)
  allocate(keyv(nsegb_c+nsegb_f+ndebug),stat=i_stat)
  call memocc(i_stat,keyv,'keyv',subname)
  ! coarse grid quantities
  call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_c,nsegb_c,keyg(1,1),keyv(1))

  ! fine grid quantities
  call segkeys(nb1,nb2,nb3,0,nb1,0,nb2,0,nb3,logrid_f,nsegb_f,keyg(1,nsegb_c+1),keyv(nsegb_c+1))

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
     wfdb%keyvloc => keyv
     wfdb%keyvglob => keyv
     wfdb%keygloc => keyg
     wfdb%keyglob => keyg
!!$     allocate(wfdb%keygloc(2,nsegb_c+nsegb_f+ndebug),stat=i_stat)
!!$     call memocc(i_stat,wfdb%keygloc,'wfdb%keygloc',subname)
!!$     allocate(wfdb%keyglob(2,nsegb_c+nsegb_f+ndebug),stat=i_stat)
!!$     call memocc(i_stat,wfdb%keyglob,'wfdb%keyglob',subname)
!!$     do i = 1, 2
!!$        do j = 1, nsegb_c+nsegb_f
!!$           wfdb%keygloc(i,j) = keyg(i,j)
!!$           wfdb%keyglob(i,j) = keyg(i,j)
!!$         end do
!!$     end do
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
          'Calculating tail corrections on ',orbs%norbp,' orbitals per processor.'
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

  !allocate the fake orbital structure for the application of projectors
  call orbitals_descriptors(0,1,1,1,0,1,1,1, &
       reshape((/0._gp,0._gp,0._gp/),(/3,1/)),(/1._gp /),orbsb,.false.)

  do iorb=1,orbs%norbp

     !build the compressed wavefunction in the enlarged box
     call transform_fortail(n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,&
          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),&
          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,Glr%wfd%keygloc(1,Glr%wfd%nseg_c+1),Glr%wfd%keyvloc(Glr%wfd%nseg_c+1),  &
          nsegb_c,nvctrb_c,keyg(1,1),keyv(1),nsegb_f,nvctrb_f,&
          keyg(1,nsegb_c+1),keyv(nsegb_c+1),&
          nbuf,psi(1,iorb),psi(Glr%wfd%nvctr_c+1,iorb),  & 
          x_c,x_f,psib(1),psib(nvctrb_c+1))

     !write(*,*) 'transform_fortail finished',iproc,iorb

     npt=2
     tail_adding: do ipt=1,npt

        if(orbs%spinsgn(iorb+orbs%isorb)>0.0d0) then
           ispin=1
        else
           ispin=2
        end if
       
        !for the tail application leave the old-fashioned hamiltonian
        !since it only deals with Free BC and thus no k-points or so

        !calculate gradient
        call applylocpotkinone(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,nbuf, &
             hgrid,nsegb_c,nsegb_f,nvctrb_c,nvctrb_f,keyg,keyv,  &
             ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,y_c,y_f,psir,  &
             psib,pot(1,1,1,ispin),hpsib,epot,ekin, &
             x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
             ibbzzx_c,ibbyyzz_c,ibbxy_ff,ibbzzx_f,ibbyyzz_f,&
             ibbzxx_c,ibbxxyy_c,ibbyz_ff,ibbzxx_f,ibbxxyy_f,nw1,nw2,ibbyyzz_r,1,1)
        !write(*,'(a,3i3,2f12.8)') 'applylocpotkinone finished',iproc,iorb,ipt,epot,ekin

        if (DistProjApply) then
           call applyprojectorsonthefly(0,orbsb,at,lr,&
                txyz,hgrid,hgrid,hgrid,wfdb,nlpspd,proj,psib,hpsib,eproj)
           !only the wavefunction descriptors must change
        else
           call applyprojectorsone(at%ntypes,at%nat,at%iatype,&
                at%psppar,at%npspcode, &
                nlpspd%nprojel,nlpspd%nproj,proj,nlpspd,&
                !nlpspd%nseg_p,nlpspd%keyg_p,nlpspd%keyv_p,nlpspd%nvctr_p,&
                !proj,&
                nsegb_c,nsegb_f,keyg,keyv,nvctrb_c,nvctrb_f,  & 
                psib,hpsib,eproj)
           !write(*,'(a,2i3,2f12.8)') 'applyprojectorsone finished',iproc,iorb,eproj,sum_tail

        end if

        !calculate residue for the single orbital
        tt=0.d0
        do i=1,nvctrb_c+7*nvctrb_f
           hpsib(i)=hpsib(i)-orbs%eval(iorb+orbs%isorb)*psib(i)
           tt=tt+hpsib(i)**2
        enddo
        tt=sqrt(tt)

        if (ipt == npt) exit tail_adding

        !calculate tail using the preconditioner as solver for the green function application
        cprecr=-orbs%eval(iorb+orbs%isorb)
        call precong(nb1,nb2,nb3,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3, &
             nsegb_c,nvctrb_c,nsegb_f,nvctrb_f,keyg,keyv, &
             ncongt,cprecr,hgrid,ibbyz_c,ibbxz_c,ibbxy_c,ibbyz_f,ibbxz_f,ibbxy_f,hpsib)
        !call plot_wf(10,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keyg,keyv,nsegb_f,nvctrb_f,  & 
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
        !n(c) ekin1=ekin
        !n(c) epot1=epot
        !n(c) eproj1=eproj
        !write(*,'(1x,a,1x,i0,f18.14)') 'norm orbital + tail',iorb,sum_tail
        !call plot_wf(20,nb1,nb2,nb3,hgrid,nsegb_c,nvctrb_c,keyg,keyv,nsegb_f,nvctrb_f,  & 
        !      txyz(1,1),txyz(2,1),txyz(3,1),psib)

        sum_tail=1.d0/sum_tail
        do i=1,nvctrb_c+7*nvctrb_f
           psib(i)=psib(i)*sum_tail
        enddo

     end do tail_adding

!!!     write(*,'(1x,a,i3,3(1x,1pe13.6),2(1x,1pe9.2))') &
!!!          'BIG: iorb,denergies,gnrm,dnorm',&
!!!          iorb,ekin-ekin1,epot-epot1,eproj-eproj1,tt,sum_tail-1.d0

     if (iproc == 0) then
        write(*,'(a)',advance='no') &
             repeat('.',((iorb+orbs%isorb)*40)/orbs%norbp-((iorb-1)*40)/orbs%norbp)
     end if
     ekin_sum=ekin_sum+ekin*orbs%occup(iorb+orbs%isorb)
     epot_sum=epot_sum+epot*orbs%occup(iorb+orbs%isorb)
     eproj_sum=eproj_sum+eproj*orbs%occup(iorb+orbs%isorb)

  end do

  if (iproc == 0) then
     write(*,'(1x,a)')'done.'
  end if
  call deallocate_orbs(orbsb,subname)

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
     i_all=-product(shape(keyg))*kind(keyg)
     deallocate(keyg,stat=i_stat)
     call memocc(i_stat,i_all,'keyg',subname)
     i_all=-product(shape(keyv))*kind(keyv)
     deallocate(keyv,stat=i_stat)
     call memocc(i_stat,i_all,'keyv',subname)
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

END SUBROUTINE CalculateTailCorrection


subroutine transform_fortail(n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,& 
     mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     msegb_c,mvctrb_c,keybg_c,keybv_c,msegb_f,mvctrb_f,keybg_f,keybv_f,  & 
     nbuf,psi_c,psi_f,psig_c,psig_f,psib_c,psib_f)
  implicit none
  integer, intent(in) :: n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f,msegb_c,mvctrb_c,msegb_f,mvctrb_f,nbuf
  integer, intent(in) :: keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  integer, intent(in) :: keybg_c(2,msegb_c),keybv_c(msegb_c),keybg_f(2,msegb_f),keybv_f(msegb_f)
  real(kind=8) :: psi_c(mvctr_c),psi_f(7,mvctr_f)
  real(kind=8) :: psib_c(mvctrb_c),psib_f(7,mvctrb_f)
  real(kind=8) :: psig_c(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf)
  real(kind=8) :: psig_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3)
  !Local variables
  integer :: iseg,jj,j0,j1,i0,i1,i2,i3,ii,i

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

END SUBROUTINE transform_fortail


subroutine transform_fortail_prev(n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3,& 
     mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     msegb_c,mvctrb_c,keybg_c,keybv_c,msegb_f,mvctrb_f,keybg_f,keybv_f,  & 
     nbuf,psi_c,psi_f,psig_c,psig_fc,psig_f,psib_c,psib_f)
  implicit none
  integer, intent(in) :: n1,n2,n3,nb1,nb2,nbfl1,nbfu1,nbfl2,nbfu2,nbfl3,nbfu3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, intent(in) :: msegb_c,mvctrb_c,msegb_f,mvctrb_f,nbuf
  integer, intent(in) :: keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  integer, intent(in) :: keybg_c(2,msegb_c),keybv_c(msegb_c),keybg_f(2,msegb_f),keybv_f(msegb_f)
  real(kind=8) :: psi_c(mvctr_c),psi_f(7,mvctr_f)
  real(kind=8) :: psib_c(mvctrb_c),psib_f(7,mvctrb_f)
  real(kind=8) :: psig_c(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf)
  real(kind=8) :: psig_fc(0:n1+2*nbuf,0:n2+2*nbuf,0:n3+2*nbuf,3),  & 
        & psig_f(7,nbfl1:nbfu1,nbfl2:nbfu2,nbfl3:nbfu3)
  !Local variables
  integer :: iseg,j0,jj,j1,i0,i1,i2,i3,ii,i

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

END SUBROUTINE transform_fortail_prev


subroutine applylocpotkinone(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     y_c,y_f,psir,  &
     psi,pot,hpsi,epot,ekin,x_c,x_f1,x_f2,x_f3,x_f,w1,w2,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r,nspinor,npot)!
  !  Applies the local potential and kinetic energy operator to one wavefunction 
  ! Input: pot,psi
  ! Output: hpsi,epot,ekin
  use module_base
  use module_interfaces
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf,nw1,nw2
  integer, intent(in) :: nseg_c,nseg_f,nvctr_c,nvctr_f,nspinor,npot
  real(gp), intent(in) :: hgrid
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  integer, dimension(2,-14:2*n3+16,0:n1), intent(in) :: ibzzx_c
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_c
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy_ff
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx_f
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz_f
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ibzxx_c
  integer, dimension(2,-14:2*n1+16,-14:2*n2+16), intent(in) :: ibxxyy_c
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ibyz_ff
  integer, dimension(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16), intent(in) :: ibzxx_f
  integer, dimension(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16), intent(in) :: ibxxyy_f
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_r
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(in) :: psi
  real(wp), dimension((2*n1+31)*(2*n2+31)*(2*n3+31),npot), intent(in) :: pot
  real(wp), dimension(nw1), intent(inout) :: w1
  real(wp), dimension(nw2), intent(inout) :: w2
  real(wp), dimension(0:n1,0:n2,0:n3,nspinor), intent(inout) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor), intent(inout) :: y_f
  real(wp), dimension((2*n1+31)*(2*n2+31)*(2*n3+31),nspinor), intent(inout) :: psir
  real(wp), dimension(0:n1,0:n2,0:n3,nspinor), intent(inout) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,nspinor), intent(inout) :: x_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3,NSPINOR),intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3,NSPINOR),intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2,NSPINOR),intent(inout) :: x_f3
  real(gp), intent(out) :: epot,ekin
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(out) :: hpsi
  !local variables
  integer :: i,idx
  real(gp) :: ekino
  real(wp), dimension(0:3) :: scal

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !call razero((2*n1+31)*(2*n2+31)*(2*n3+31)*nspinor,psir)

  do idx=1,nspinor  
     call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
          nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          scal,psi(1,IDX),psi(nvctr_c+1,IDX),  &
          x_c(0,0,0,idx),x_f(1,nfl1,nfl2,nfl3,idx),&
          x_f1(nfl1,nfl2,nfl3,idx),x_f2(nfl2,nfl1,nfl3,idx),x_f3(nfl3,nfl1,nfl2,idx))
     
     call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          w1,w2,x_c(0,0,0,idx),x_f(1,nfl1,nfl2,nfl3,idx), & 
          psir(1,IDX),ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
     
  end do

  call apply_potential(n1,n2,n3,1,1,1,nbuf,nspinor,npot,psir,pot,epot,&
       ibyyzz_r) !optional

!!  epot=0.0_gp
!!  if (nspinor==1 .or. nspinor == 2) then
!!     do ispinor=1,nspinor
!!        if (nbuf == 0) then
!!           call realspace(ibyyzz_r,pot,psir(1,ispinor),epots,n1,n2,n3)
!!        else
!!           !this is for the tails. In principle it should work only for 
!!           call realspace_nbuf(ibyyzz_r,pot,psir(1,ispinor),epot,n1,n2,n3,nbuf)
!!        endif
!!           epot=epot+epots
!!        end do
!!  else
!!     call realspaceINPLACE(ibyyzz_r,pot,psir,epot,n1,n2,n3)
!!  end if
  
  ekin=0.0_gp
  do idx=1,nspinor
     call comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          w1,w2,psir(1,IDX),&
          ibxy_c,ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
          y_c(0,0,0,IDX),y_f(1,nfl1,nfl2,nfl3,IDX))!,ibyz_c,ibyz_f)
     
     call ConvolkineticT(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
          hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, &
          x_c(0,0,0,IDX),x_f(1,nfl1,nfl2,nfl3,IDX),&
          y_c(0,0,0,IDX),y_f(1,nfl1,nfl2,nfl3,IDX),EKINO, &
          x_f1(nfl1,nfl2,nfl3,IDX),x_f2(nfl2,nfl1,nfl3,IDX),x_f3(nfl3,nfl1,nfl2,IDX))
     ekin=ekin+ekino
     
     call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
          nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          scal,y_c(0,0,0,IDX),y_f(1,nfl1,nfl2,nfl3,IDX),hpsi(1,IDX),hpsi(nvctr_c+1,IDX))
  end do
  
END SUBROUTINE applylocpotkinone


!> Applies all the projectors onto a single wavefunction
!! Input: psi_c,psi_f
!! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initialized to zero at the beginning)
subroutine applyprojectorsone(ntypes,nat,iatype,psppar,npspcode, &
     nprojel,nproj,&
     !nseg_p,keyg_p,keyv_p,nvctr_p,&
     proj,nlpspd,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,&
     psi,hpsi,eproj)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ntypes,nat,nprojel,nproj,nseg_c,nseg_f,nvctr_c,nvctr_f
  integer, dimension(ntypes), intent(in) :: npspcode
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
!!$  integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
!!$  integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
!!$  integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
  real(gp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(in) :: psi
  real(wp), dimension(nprojel), intent(in) :: proj
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(inout) :: hpsi
  real(gp), intent(out) :: eproj
  !local variables
  integer :: i,l,iat,iproj,istart_c,mbseg_c,mbseg_f,jseg_c,mbvctr_c,mbvctr_f,ityp !n(c) jseg_f

  ! loop over all projectors
  iproj=0
  eproj=0.0_gp
  istart_c=1
  do iat=1,nat
     call plr_segs_and_vctrs(nlpspd%plr(iat),&
          mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
     jseg_c=1

!!$     mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
!!$     mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
!!$     jseg_c=nseg_p(2*iat-2)+1
!!$     !n(c) jseg_f=nseg_p(2*iat-1)+1
!!$     mbvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
!!$     mbvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)

     ityp=iatype(iat)
     !GTH and HGH pseudopotentials
     do l=1,4
        do i=1,3
           if (psppar(l,i,ityp) /= 0.0_gp) then
              !in this case the ncplx value is 1 mandatory 
              call applyprojector(1,l,i,psppar(0,0,ityp),npspcode(ityp),&
                   nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,&
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                   nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                   nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
!!$                   keyv_p(jseg_c),keyg_p(1,jseg_c),&
                   proj(istart_c),psi,hpsi,eproj)
              iproj=iproj+2*l-1
              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)
           end if
        enddo
     enddo
  enddo
  if (iproj /= nproj) stop '1:applyprojectorsone'
  if (istart_c-1 /= nprojel) stop '2:applyprojectorsone'

END SUBROUTINE applyprojectorsone
