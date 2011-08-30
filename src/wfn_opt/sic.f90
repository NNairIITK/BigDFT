!> Construct a Self-Interaction-Corrected potential based on the 
!! Perdew-Zunger prescription (Phys. Rev. B 23, 10, 5048 (1981))
subroutine PZ_SIC_potential(iorb,lr,orbs,ixc,hxh,hyh,hzh,pkernel,psir,vpsir,eSICi)
  use module_base
  use module_types
  use module_interfaces
  use Poisson_Solver
  use module_xc
  implicit none
  integer, intent(in) :: iorb,ixc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor), intent(in) :: psir
  real(dp), dimension(*), intent(in) :: pkernel
  real(gp), intent(out) :: eSICi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor), intent(out) :: vpsir
  !local variables
  character(len=*), parameter :: subname='PZ_SIC_potential' 
  integer :: npsir,nspinn,ncomplex,i_all,i_stat,ispin,icomplex,jproc,nproc
  real(gp) :: vexi,eexi,ehi,spinval,hfac
  integer, dimension(:,:), allocatable :: nscarr_fake
  real(dp), dimension(:,:), allocatable :: rhopoti,vSICi
  real(dp), dimension(:), pointer :: rhocore_fake


  !fake nscatterarr array for compatibility with partial_density interface
  !monoprocessor calculation
  nproc=1
  allocate(nscarr_fake(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscarr_fake,'nscarr_fake',subname)
  !fill it just for completeness, even if not used
  do jproc=0,nproc-1
     nscarr_fake(jproc,1)=lr%d%n3i
     nscarr_fake(jproc,2)=lr%d%n3i
     nscarr_fake(jproc,3)=0
     nscarr_fake(jproc,4)=0
  end do
  !fake rhocore array
  nullify(rhocore_fake)
  
  !components of wavefunction in real space which must be considered simultaneously
  !and components of the charge density
  if (orbs%nspinor ==4) then
     npsir=4
     nspinn=4
     ncomplex=0
  else
     npsir=1
     nspinn=orbs%nspin
     ncomplex=orbs%nspinor-1
  end if
  !with this definition (ncomplex+1)*npsir is always equal to nspinor

  !npsir is equal to npot in local_hamiltonian

  !wavefunction squares (spin dependent)
  allocate(rhopoti(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinn+ndebug),stat=i_stat)
  call memocc(i_stat,rhopoti,'rhopoti',subname)
  !initialize the charge density
  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  !otherwise use libXC routine
  call xc_init_rho(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspinn,rhopoti,nproc)

  allocate(vSICi(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinn+ndebug),stat=i_stat)
  call memocc(i_stat,vSICi,'vSICi',subname)

  !occupation number and k-point weight
  hfac=orbs%kwgts(orbs%iokpt(iorb))*(orbs%occup(orbs%isorb+iorb)/(hxh*hyh*hzh))
  !value of the spin state
  spinval=orbs%spinsgn(orbs%isorb+iorb)
  !spin up or down depending of spinval
  if (spinval >0.0_gp) then
     ispin=1
  else
     ispin=2
  end if
  eSICi=0.0_gp
  if (hfac /= 0.d0) then

     !copy the psir function in the vpsir array
     call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir(1,1),1,vpsir(1,1),1)

     !accumulate the density in the rhopoti array, in the complex case
     do icomplex=0,ncomplex

        !print *,'iorb,nrm',iorb,&
        !     nrm2(lr%d%n1i*lr%d%n2i*lr%d%n3i*npsir,psir(1,1),1)

        select case(lr%geocode)
        case('F')
           call partial_density_free(.false.,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                npsir,nspinn,lr%d%n3i,&
                hfac,nscarr_fake,spinval,psir(1,icomplex*npsir+1),rhopoti,lr%bounds%ibyyzz_r)
        case('P')
           call partial_density(.false.,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                npsir,nspinn,lr%d%n3i,&
                hfac,nscarr_fake,spinval,psir(1,icomplex*npsir+1),rhopoti)
        case('S')
           call partial_density(.false.,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                npsir,nspinn,lr%d%n3i,&
                hfac,nscarr_fake,spinval,psir(1,icomplex*npsir+1),rhopoti)
        end select
     end do

     !here the density should be transformed into a potential which is associated to the orbital
     if(orbs%nspinor==4) then
        !this wrapper can be inserted inside the XC_potential routine
        call PSolverNC(lr%geocode,'D',0,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,lr%d%n3i,&
             ixc,hxh,hyh,hzh,&
             rhopoti,pkernel,rhopoti,ehi,eexi,vexi,0.d0,.false.,4)
        !the potential is here ready to be applied to psir
     else
        call XC_potential(lr%geocode,'D',0,1,&
             lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,hxh,hyh,hzh,&
             rhopoti,eexi,vexi,orbs%nspin,rhocore_fake,vSICi) 

        call H_potential(lr%geocode,'D',0,1,&
             lr%d%n1i,lr%d%n2i,lr%d%n3i,hxh,hyh,hzh,&
             rhopoti,pkernel,rhopoti,ehi,0.0_dp,.false.,&
             quiet='YES') !optional argument

        !start to fill the potential with the hartree potential

        !sum the two potentials in the potential array
        !fill the other part, for spin polarised
        if (orbs%nspin == 2) then
           call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i,rhopoti(1,1),1,&
                rhopoti(1,2),1)
        end if
        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,1.0_dp,vSICi(1,1),1,&
             rhopoti(1,1),1)
        !note: this filling is redundant since the orbital has one one of the spins by definition

     end if

     !apply the potential to the psir wavefunction and calculate potential energy
     select case(lr%geocode)
     case('F')
        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npsir,vpsir,&
             rhopoti(1,ispin),eSICi,&
             lr%bounds%ibyyzz_r) !optional

     case('P') 
        !here the hybrid BC act the same way
        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,0,0,0,orbs%nspinor,npsir,vpsir,&
             rhopoti(1,ispin),eSICi)

     case('S')

        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,1,0,0,orbs%nspinor,npsir,vpsir,&
             rhopoti(1,ispin),eSICi)
     end select

  else
     !put to zero the corresponding potential
     call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,vpsir(1,1))
  end if


  i_all=-product(shape(rhopoti))*kind(rhopoti)
  deallocate(rhopoti,stat=i_stat)
  call memocc(i_stat,i_all,'rhopoti',subname)
  i_all=-product(shape(vSICi))*kind(vSICi)
  deallocate(vSICi,stat=i_stat)
  call memocc(i_stat,i_all,'vSICi',subname)
  
  i_all=-product(shape(nscarr_fake))*kind(nscarr_fake)
  deallocate(nscarr_fake,stat=i_stat)
  call memocc(i_stat,i_all,'nscarr_fake',subname)

end subroutine PZ_SIC_potential
