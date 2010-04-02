!!****f* BigDFT/local_analysis
!! FUNCTION
!!    Perform all the projection associated to local variables
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine local_analysis(iproc,nproc,hx,hy,hz,shift,lr,orbs,orbsv,psi,psivirt)
  use module_base
  use module_types
  use module_interfaces, except_this_one => local_analysis
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs,orbsv
  real(gp), dimension(3),intent(in) :: shift
  real(wp), dimension(:), pointer :: psi,psivirt
  !local variables
  character(len=*), parameter :: subname='local_analysis'
  integer :: nelec,norb,norbu,norbd,iunit,i_all,i_stat,norbpv,iat
  type(input_variables) :: inc
  type(atoms_data) :: atc
  type(gaussian_basis) :: G
  real(wp), dimension(:,:), allocatable :: allpsigau,dualcoeffs
  real(gp), dimension(:,:), allocatable :: radii_cf_fake,thetaphi
  real(gp), dimension(:), pointer :: Gocc
  real(gp), dimension(:,:), pointer :: cxyz

  !the number of virtual orbitals in parallel is known only if orbsv%norb>0
  if (orbsv%norb >0) then
     norbpv=orbsv%norbp
  else
     norbpv=0
  end if

  !define the local basis starting from the input files
  call read_input_variables(iproc,'posinp','input.dft','','','',inc,atc,cxyz)

  allocate(radii_cf_fake(atc%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf_fake,'radii_cf_fake',subname)

  call read_system_variables('input.occup',iproc,inc,atc,radii_cf_fake,nelec,&
       norb,norbu,norbd,iunit)

  !shift the positions with the same value of the original positions
  do iat=1,atc%nat
     cxyz(1,iat)=cxyz(1,iat)-shift(1)
     cxyz(2,iat)=cxyz(2,iat)-shift(2)
     cxyz(3,iat)=cxyz(3,iat)-shift(3)
  end do
  
  nullify(G%rxyz)

  !extract the gaussian basis from the pseudowavefunctions
  call gaussian_pswf_basis(21,iproc,inc%nspin,atc,cxyz,G,Gocc)

  allocate(thetaphi(2,G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,thetaphi,'thetaphi',subname)
  call razero(2*G%nat,thetaphi)
  allocate(allpsigau(G%ncoeff,orbs%norbp+norbpv+ndebug),stat=i_stat)
  call memocc(i_stat,allpsigau,'allpsigau',subname)

  !this routine should be simplified like gaussians_to_wavelets
  call wavelets_to_gaussians(lr%geocode,orbs%norbp,orbs%nspinor,&
       lr%d%n1,lr%d%n2,lr%d%n3,G,thetaphi,hx,hy,hz,lr%wfd,psi,allpsigau)

  !the same can be done for virtual orbitals if orbsv%norb > 0
  if (orbsv%norb > 0) then
     call wavelets_to_gaussians(lr%geocode,norbpv,orbsv%nspinor,&
          lr%d%n1,lr%d%n2,lr%d%n3,G,thetaphi,hx,hy,hz,lr%wfd,psivirt,&
          allpsigau(1,orbs%norbp+1))
  end if
  !calculate dual coefficients
  allocate(dualcoeffs(G%ncoeff,orbs%norbp+norbpv+ndebug),stat=i_stat)
  call memocc(i_stat,dualcoeffs,'dualcoeffs',subname)
  call dcopy(G%ncoeff*(orbs%norbp+norbpv),allpsigau,1,dualcoeffs,1)
  !build dual coefficients
  call dual_gaussian_coefficients(orbs%norbp+norbpv,G,dualcoeffs)

  !here we can calculate the Mulliken charge population
  !for any of the elements of the basis, ordered by angular momentum
  !do that only for the occupied orbitals
  call mulliken_charge_population(iproc,nproc,orbs,Gocc,G,allpsigau,dualcoeffs)

  !also partial density of states can be analysed here

  call deallocate_gwf(G,subname)
  nullify(G%rxyz)

  i_all=-product(shape(allpsigau))*kind(allpsigau)
  deallocate(allpsigau,stat=i_stat)
  call memocc(i_stat,i_all,'allpsigau',subname)
  i_all=-product(shape(dualcoeffs))*kind(dualcoeffs)
  deallocate(dualcoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'dualcoeffs',subname)


  i_all=-product(shape(thetaphi))*kind(thetaphi)
  deallocate(thetaphi,stat=i_stat)
  call memocc(i_stat,i_all,'thetaphi',subname)

  !deallocate the auxiliary structures for the calculations
  call deallocate_atoms(atc,subname) 
  call deallocate_atoms_scf(atc,subname) 
  call free_input_variables(inc)
  i_all=-product(shape(radii_cf_fake))*kind(radii_cf_fake)
  deallocate(radii_cf_fake,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf_fake',subname)
  i_all=-product(shape(cxyz))*kind(cxyz)
  deallocate(cxyz,stat=i_stat)
  call memocc(i_stat,i_all,'cxyz',subname)
  i_all=-product(shape(Gocc))*kind(Gocc)
  deallocate(Gocc,stat=i_stat)
  call memocc(i_stat,i_all,'Gocc',subname)
  
END SUBROUTINE local_analysis
!!***


!!****f* BigDFT/mulliken_charge_population
!! SOURCE
!! 
subroutine mulliken_charge_population(iproc,nproc,orbs,Gocc,G,coeff,duals)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(gaussian_basis), intent(in) :: G
  real(gp), dimension(G%ncoeff), intent(in) :: Gocc
  real(wp), dimension(G%ncoeff,orbs%norbp), intent(in) :: coeff,duals
  !local variables
  character(len=*), parameter :: subname='mulliken_charge_population'
  character(len=11) :: shname
  integer :: icoeff,i_all,i_stat,ierr,ishell,iexpo,iat,l,ng,iorb,isat,m,ispin,ig,nchannels
  real(wp) :: msum,rad,radnorm,r,sumch
  real(wp), dimension(2) :: msumiat
  real(wp), dimension(:,:), allocatable :: mchg
  
  !allocate both for spins up and down
  allocate(mchg(G%ncoeff,2+ndebug),stat=i_stat)
  call memocc(i_stat,mchg,'mchg',subname)

  !for any of the orbitals calculate the Mulliken charge
  do icoeff=1,G%ncoeff
     mchg(icoeff,1)=0.0_wp
     mchg(icoeff,2)=0.0_wp
     !print '(a,100(1pe12.5))','icoeff,iorb',coeff(icoeff,:)
     !print '(a,100(1pe12.5))','idualc,iorb',duals(icoeff,:)
     do iorb=1,orbs%norbp
        if (orbs%spinsgn(orbs%isorb+iorb) == 1.0_gp) then
           ispin=1
        else
           ispin=2
        end if
        mchg(icoeff,ispin)=mchg(icoeff,ispin)+&
             orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)*&
             coeff(icoeff,iorb)*duals(icoeff,iorb)
             !duals(icoeff,iorb)**2

     end do
  end do

  !reduce the results
  if (nproc > 1) then
     call mpiallred(mchg(1,1),2*G%ncoeff,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  if (iproc == 0) then
     write(*,'(1x,a)')repeat('-',48)//' Mulliken Charge Population Analysis'
     write(*,'(1x,a)')'Center No. |    Shell    | Rad (AU) | Chg (up) | Chg (down) | Net Pol  | Gross Chg'
  end if

!  do iorb=1,orbs%norbp  
!     msum=0.0_wp
!     do icoeff=1,G%ncoeff
!        msum=msum+coeff(icoeff,iorb)*duals(icoeff,iorb)
!     end do
!     print *,'total sum,iorb',iorb,msum,&
!          orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)
!  end do

  !print the results as a function of the shell
  ishell=0
  iexpo=1
  icoeff=1
  msum=0.0_wp
  do iat=1,G%nat
     msumiat(1)=0.0_wp
     msumiat(2)=0.0_wp
     nchannels=0
     sumch=0.0_gp
     do isat=1,G%nshell(iat)
        ishell=ishell+1
        ng=G%ndoc(ishell)
        l=G%nam(ishell)
        !calculate mean radius (a.u.)
        rad=0.0_wp
        radnorm=0.0_wp
        do ig=1,ng
           r=G%xp(iexpo)
           rad=rad+(G%psiat(iexpo))**2*r
           radnorm=radnorm+(G%psiat(iexpo))**2
           iexpo=iexpo+1
        end do
        rad=rad/radnorm
        do m=1,2*l-1
           call shell_name(l,m,shname)
           msumiat(1)=msumiat(1)+mchg(icoeff,1)
           msumiat(2)=msumiat(2)+mchg(icoeff,2)
           if (iproc == 0) then
              write(*,'(1x,(i6),5x,a,2x,a,a,1x,f7.2,2x,2("|",1x,f8.5,1x),2(a,f8.5))')&
                   iat,'|',shname,'|',rad,(mchg(icoeff,ispin),ispin=1,2),'  | ',&
                   mchg(icoeff,1)-mchg(icoeff,2),' | ',Gocc(icoeff)-(mchg(icoeff,1)+mchg(icoeff,2))
           end if
           sumch=sumch+Gocc(icoeff)
           icoeff=icoeff+1
           nchannels=nchannels+1
        end do
     end do
     if (iproc == 0) write(*,'(15x,a,2("|",1x,f8.5,1x),2(a,f8.5))')&
          '  Center Quantities : ',&
          (msumiat(ispin),ispin=1,2),'  | ',msumiat(1)-msumiat(2),' | ',&
          sumch-(msumiat(1)+msumiat(2))
     msum=msum+msumiat(1)+msumiat(2)
     if (iproc == 0) write(*,'(1x,a)')repeat('-',82)
  end do

  if (iproc == 0) write(*,'(24x,a,f21.12)')'Total Charge considered on the centers: ',msum
   
  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  i_all=-product(shape(mchg))*kind(mchg)
  deallocate(mchg,stat=i_stat)
  call memocc(i_stat,i_all,'mchg',subname)
  
END SUBROUTINE mulliken_charge_population
!!***


!!****f* BigDFT/shell_name
!! SOURCE
!! 
subroutine shell_name(l,m,name)
  implicit none
  integer, intent(in) :: l,m
  character(len=11), intent(out) :: name
  
  select case(l)
  case(1)
     name(1:1)='s'
     select case(m)
     case(1)
        name(2:11)='          '
     case default
        stop 'wrong m'
     end select
  case(2)
     name(1:1)='p'
     select case(m)
     case(1)
        name(2:11)='x         '
     case(2)
        name(2:11)='y         '
     case(3)
        name(2:11)='z         '
     case default
        stop 'wrong m'
     end select
  case(3)
     name(1:1)='d'        
     select case(m)
     case(1)
        name(2:11)='yz        '
     case(2)
        name(2:11)='xz        '
     case(3)
        name(2:11)='xy        '
     case(4)
        name(2:11)='x2-y2     '
     case(5)
        name(2:11)='2z2-r2    '
     case default
        stop 'wrong m'
     end select
  case(4)
     name(1:1)='f'        
     select case(m)
     case(1)
        name(2:11)='x(r2-5z2) '
     case(2)
        name(2:11)='y(r2-5z2) '
     case(3)
        name(2:11)='z(3r2-5z2)'
     case(4)
        name(2:11)='x(x2-3y2) '
     case(5)
        name(2:11)='y(y2-3x2) '
     case(6)
        name(2:11)='z(x2-y2)  '
     case(7)
        name(2:11)='xyz       '
     case default
        stop 'wrong m'
     end select
  case default
     stop 'l not recognized'
  end select

END SUBROUTINE shell_name
!!***
