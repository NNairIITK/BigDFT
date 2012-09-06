!> Construct a Self-Interaction-Corrected potential based on the 
!! Perdew-Zunger prescription (Phys. Rev. B 23, 10, 5048 (1981))
subroutine PZ_SIC_potential(iorb,lr,orbs,ixc,hxh,hyh,hzh,pkernel,psir,vpsir,eSICi,eSIC_DCi)
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
  type(coulomb_operator), intent(in) :: pkernel
  real(gp), intent(out) :: eSICi,eSIC_DCi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor), intent(out) :: vpsir
  !local variables
  character(len=*), parameter :: subname='PZ_SIC_potential' 
  integer :: npsir,nspinn,ncomplex,i_all,i_stat,ispin,icomplex,jproc,nproc
  real(gp) :: spinval,hfac,fi,vexi,eexi,ehi
  integer, dimension(:,:), allocatable :: nscarr_fake
  real(dp), dimension(:,:), allocatable :: rhopoti,vSICi
  real(dp), dimension(:,:,:,:), pointer :: rhocore_fake
  real(dp), dimension(6) :: xcstr

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

  !occupation number and k-point weight of the given orbital
  fi=orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)

  hfac=fi/(hxh*hyh*hzh)
  !value of the spin state
  spinval=orbs%spinsgn(orbs%isorb+iorb)
  !spin up or down depending of spinval
  if (spinval >0.0_gp) then
     ispin=1
  else
     ispin=2
  end if
  eSICi=0.0_gp
  ehi=0.0_gp
  eexi=0.0_gp
  vexi=0.0_gp
  eSIC_DCi=0.0_gp
  if (fi /= 0.d0) then

     !copy the psir function in the vpsir array
     call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir(1,1),1,vpsir(1,1),1)

     !here the psir_to_rhoi routine can be used
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
             rhopoti,pkernel%kernel,rhopoti,ehi,eexi,vexi,0.d0,.false.,4)
        !the potential is here ready to be applied to psir
     else
        call XC_potential(lr%geocode,'D',0,1,MPI_COMM_WORLD,&
             lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,hxh,hyh,hzh,&
             rhopoti,eexi,vexi,orbs%nspin,rhocore_fake,vSICi,xcstr) 

        call H_potential('D',pkernel,rhopoti,rhopoti,ehi,0.0_dp,.false.,&
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

     !calculate the contribution to the double-counting SIC energy (to be multiplied by alphaSIC)
     
     eSIC_DCi=-ehi-vexi+eexi

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

!> Construct a Self-Interaction-Corrected potential based on the 
!! Koopmans' correction for DFT (Phys. Rev. B 82 115121 (2010))
!! The input potential is contains the full potential and the charge density (if norbp > 0)
!! the output potential is orbital-dependent and is given by the NK Hamiltonian
!! @ param poti the density in the input, the orbital-wise potential in the output
!!              unless wxdsave is present. In this case poti in unchanged on exit
subroutine NK_SIC_potential(lr,orbs,ixc,fref,hxh,hyh,hzh,pkernel,psi,poti,eSIC_DC,potandrho,wxdsave)
  use module_base
  use module_types
  use module_xc
  use module_interfaces, except_this_one => NK_SIC_potential
  use Poisson_Solver
  implicit none
  integer, intent(in) :: ixc
  real(gp), intent(in) :: hxh,hyh,hzh,fref
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  type(coulomb_operator), intent(in) :: pkernel
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(wp), dimension((lr%d%n1i*lr%d%n2i*lr%d%n3i*((orbs%nspinor/3)*3+1)),max(orbs%norbp,orbs%nspin)), intent(inout) :: poti
  real(gp), intent(out) :: eSIC_DC
  real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,2*orbs%nspin), intent(in), optional :: potandrho !< array which should be used in the virtual case
  real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin), intent(out), optional :: wxdsave !< memory space save the wxd term (preparation of the virtual case)
  !local variables
  character(len=*), parameter :: subname='NK_SIC_potential' 
  logical :: virtual,savewxd
  integer :: npot,i_all,i_stat,ispin,i,jspin,ierr,iorb,ispinor,i1,i2,i3,nproc
  real(gp) :: spinval,fi,oneoh,vexi,eexi,ehi,eexu,vexu,eSIC_DCi,fac1,fac2,rnorboccp,constadd
  type(workarr_sumrho) :: w
  real(dp), dimension(:,:), allocatable :: ni,deltarho,vxci,rho,psir,wxd
  real(dp), dimension(:,:,:,:), allocatable :: fxci
  real(dp), dimension(:,:,:,:), pointer :: rhocore_fake
  real(dp), dimension(6) :: xcstr
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
 
  !stop for non-collinear case (lacking of unified XC routine)
  if (orbs%nspinor==4) stop 'SIC not available for non-collinear case'

  virtual=present(potandrho)

  savewxd=present(wxdsave)

  if (virtual .and. savewxd) stop 'NKpot: options are mutually exclusive'

  !XC potential and work array for the cross-derivative summation (to be used in the occupied case and for saving)
  allocate(wxd(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin+ndebug),stat=i_stat)
  call memocc(i_stat,wxd,'wxd',subname)
  if (orbs%norbp==0) call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,wxd(1,1))

  !quick return if norbp=0 (not possible due to the MPI_ALLREDUCE)
  if (orbs%norbp /= 0) then

     !npot=1
     !if (orbs%nspinor==4) npot=4
     npot=(orbs%nspinor/3)*3+1 !alternative form for the definition above

     !fake rhocore array
     nullify(rhocore_fake)

     !reals-space wavefunction
     allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)

     !total density of the system (for the virtual wfns it is already in density)
     if (.not. virtual) then
        allocate(rho(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin+ndebug),stat=i_stat)
        call memocc(i_stat,rho,'rho',subname)
     end if

     !wavefunction squares (spin dependent)
     allocate(ni(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,ni,'ni',subname)
     !work array deltarho
     allocate(deltarho(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,deltarho,'deltarho',subname)
     !orbital-wise XC potential
     allocate(vxci(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,vxci,'vxci',subname)
     allocate(fxci(lr%d%n1i,lr%d%n2i,lr%d%n3i,orbs%nspin+1+ndebug),stat=i_stat)
     call memocc(i_stat,fxci,'fxci',subname)

     !take total density from the first poti component or the second potand rho component
     if (virtual) then
        call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,potandrho(1,orbs%nspin+1),1,&
             deltarho(1,1),1)
     else
        call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,poti(1,1),1,&
             deltarho(1,1),1)
     end if
     !print *,'here',poti(1,1),deltarho(1,1)

     !put the XC potential in the wxd term, which is the same for all the orbitals
     call XC_potential(lr%geocode,'D',0,1,MPI_COMM_WORLD,&
          lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,hxh,hyh,hzh,&
          deltarho,eexu,vexu,orbs%nspin,rhocore_fake,wxd,xcstr)

     if (.not. virtual) then
        !rescale wxd, pay attention to the occupation of the orbitals
        rnorboccp=0
        do iorb=1,orbs%norbp
           rnorboccp=rnorboccp+&
                orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)
        end do
     else
        rnorboccp=-1.0_gp !Vxc should be positive in wxd for virtual case
     end if

     call dscal(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,-rnorboccp,wxd(1,1),1)

     !put the density in the rho array
     if (.not. virtual) then
        call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,poti(1,1),1,&
             rho(1,1),1)
     end if

     !initalize the double-counting SIC corrections
     eSIC_DC=0.0_gp
     call initialize_work_arrays_sumrho(lr,w)

     do iorb=1,orbs%norbp

        !initialize the charge density
        !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
        !otherwise use libXC routine
        call xc_init_rho(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,ni,1)

        !occupation number and k-point weight of the given orbital
        fi=orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(orbs%isorb+iorb)
        if (virtual) fi=0.0_gp
        oneoh=1.0_gp/(hxh*hyh*hzh)

        !value of the spin state
        spinval=orbs%spinsgn(orbs%isorb+iorb)

        !spin up or down depending of spinval
        if (spinval >0.0_gp) then
           ispin=1
        else
           ispin=2
        end if

        !do something also if the occupation number is  zero
        do ispinor=1,orbs%nspinor
           call daub_to_isf(lr,w,psi(1,ispinor,iorb),psir(1,ispinor))
        end do

        call psir_to_rhoi(oneoh,spinval,orbs%nspin,orbs%nspinor,lr,psir,ni)

        !---STEP 1
        !print *,'step1',iorb+orbs%isorb
        !build the reference density
        !rhoref=rho+(fref-fi)*ni
        if (virtual) then
           call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,potandrho(1,orbs%nspin+1),1,&
                deltarho(1,1),1)
        else
           call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,rho(1,1),1,&
                deltarho(1,1),1)
           end if
        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,fref-fi,ni(1,1),1,deltarho(1,1),1)

        !call xc_clean_rho(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,deltarho,1)

        !if (savewxd) call xc_clean_rho(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,deltarho,1)

        !calculate its vXC and fXC
        call XC_potential(lr%geocode,'D',0,1,MPI_COMM_WORLD,&
             lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,hxh,hyh,hzh,&
             deltarho,eexi,vexi,orbs%nspin,rhocore_fake,vxci,xcstr,fxci)

        !copy the relevant part of Vxc[rhoref] in the potential
        if (.not. savewxd) &
             call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*npot,vxci(1,ispin),1,poti(1,iorb),1)

        !start the definition of the wxd term, in the diagonal part
        !calculate the contribution to the Double Counting energy
        !put it in deltarho array, then put it in fxci
        if (fi /= 0.0_gp) then
           fac1=fi
           fac2=fref/fi
           eSIC_DCi=dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,vxci(1,ispin),1,ni(1,ispin),1)
        else
           fac1=1.0_gp
           fac2=fref
           eSIC_DCi=0.0_gp
        end if

        constadd=0.0_gp
        do jspin=1,orbs%nspin
           do i3=1,lr%d%n3i
              do i2=1,lr%d%n2i
                 do i1=1,lr%d%n1i
                    i=i1+(i2-1)*lr%d%n1i+(i3-1)*lr%d%n1i*lr%d%n2i
                    !use the convention fxc(ispin,jspin)=fxc(:,ispin+jspin-1)
                    deltarho(i,jspin)=fac1*fxci(i1,i2,i3,ispin+jspin-1)*ni(i,ispin)
                    if (ispin==jspin) then 
                       constadd=constadd+deltarho(i,jspin)*ni(i,ispin)
                    end if
                 end do
              end do
           end do
        end do
        eSIC_DCi=eSIC_DCi*fac1*hxh*hyh*hzh
        constadd=constadd*fac2*hxh*hyh*hzh

        !put the term in the potential
        !poti=Vxc[rhoref]+fref ni fxc[rhoref]
        if (.not. savewxd) &
             call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*npot,fac2,deltarho(1,ispin),1,poti(1,iorb),1)

        !start accumulating the contribution to the wxd array in the fxci array
        if (fi /= 0.0_gp) then
           call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,deltarho(1,1),1,fxci(1,1,1,1),1)
        else
           !put the total Vxc in the fxc, which should be copied into pot afterwards
           call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,wxd(1,1),1,fxci(1,1,1,1),1)
        end if

        !---STEP 2
        !print *,'step2',iorb+orbs%isorb
        !build the density variation (only if fi/=0)
        !deltarho=rho-rhoi
        !pay attention here to the XC routines for zero values of the density difference
        if (fi /= 0.0_gp) then
           call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,rho(1,1),1,&
                deltarho(1,1),1)
           call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,-fi,ni(1,1),1,deltarho(1,1),1)

           call xc_clean_rho(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,deltarho,1)

           !calculate its XC potential
           call XC_potential(lr%geocode,'D',0,1,MPI_COMM_WORLD,&
                lr%d%n1i,lr%d%n2i,lr%d%n3i,ixc,hxh,hyh,hzh,&
                deltarho,eexi,vexi,orbs%nspin,rhocore_fake,vxci,xcstr) 
           !saves the values for the double-counting term
           eSIC_DC=eSIC_DC+eexi-eexu+eSIC_DCi

           !and accumulate it in the local cross-diagonal value
           !wxd=vxc[rho-rhoi]+fxci[rhoref]*rhoi
           call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,1.0_wp,vxci(1,1),1,fxci(1,1,1,1),1)

           !wxd=vxc[rho-rhoi]-vxc[rho]+fxci[rhoref]*rhoi
           !subtract again vxc to make the correct wxd (the vxc  should have been in poti also, so it has not been subtracted)
           !to do this update the wxd array which has vxc inside
           call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,1.0_wp,fxci(1,1,1,1),1,wxd(1,1),1)
        end if

        !subtract the diagonal part of the off-diagonal term (in order to avoid the i=j term after summation)
        !in the virtual case this is vxc, which is needed instaed of Vxc[rho-rhoi]
        if (.not. savewxd) call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,-1.0_wp,fxci(1,1,1,ispin),1,poti(1,iorb),1)


        !--STEP 3
        !print *,'step3',iorb+orbs%isorb
        !accumulate the diagonal term in the local array
        !vi=(2*fref-fi)Vh[ni]+fref*ni*fxc[rhoref]+vxc[rhoref]-vxc[rho]
        !create the hartree potential of rhoi, if necessary
        !this should be done in any case due to the constant term which should be added to the potential
        !however, such constant term is ininfluent so we might remove it from here
        !this step is useless if savewxd is activated
        if (.not. savewxd) then
           call H_potential('D',pkernel,ni(1,ispin),ni(1,ispin),ehi,0.0_dp,.false.,&
                quiet='YES') !optional argument
           !saves eexi for the double-counting term
           eSIC_DC=eSIC_DC+fi*(2.0_wp*fref-fi)*ehi
           
           constadd=constadd+2.0_gp*fref*ehi

           !print *,'constadd,iorb',constadd,iorb+orbs%isorb
           !add this result to the potential, and subtract also the constant
           !if (.not. virtual) then
              do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                 poti(i,iorb)=poti(i,iorb)-constadd
              end do
           !end if
           call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,(2.0_wp*fref-fi),ni(1,ispin),1,poti(1,iorb),1)
           
           !add the occupied-space potential to the poti
           if (virtual) call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*npot,1.0_wp,potandrho(1,ispin),1,poti(1,iorb),1)
        end if
        

     end do

     call deallocate_work_arrays_sumrho(w)

  end if

  if (.not. virtual) then
     !sum up the results of the off diagonal term
     if (nproc >1) call mpiallred(wxd(1,1),lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,MPI_SUM,MPI_COMM_WORLD,ierr)

     if (.not. savewxd) then
        !add to the potential orbital per orbital
        do iorb=1,orbs%norbp
           !value of the spin state
           spinval=orbs%spinsgn(orbs%isorb+iorb)
           !spin up or down depending of spinval
           if (spinval >0.0_gp) then
              ispin=1
           else
              ispin=2
           end if
           
           !reduce the result in the proper part of the array
           call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,1.0_wp,wxd(1,ispin),1,poti(1,iorb),1)
           
        end do
     end if
  end if

  !copy the wxd array in the wxdsave memory space
  if (savewxd) call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin,wxd(1,1),1,wxdsave(1,1),1)

  i_all=-product(shape(ni))*kind(ni)
  deallocate(ni,stat=i_stat)
  call memocc(i_stat,i_all,'ni',subname)
  i_all=-product(shape(vxci))*kind(vxci)
  deallocate(vxci,stat=i_stat)
  call memocc(i_stat,i_all,'vxci',subname)
  i_all=-product(shape(wxd))*kind(wxd)
  deallocate(wxd,stat=i_stat)
  call memocc(i_stat,i_all,'wxd',subname)
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  i_all=-product(shape(fxci))*kind(fxci)
  deallocate(fxci,stat=i_stat)
  call memocc(i_stat,i_all,'fxci',subname)
  if (.not. virtual) then
     i_all=-product(shape(rho))*kind(rho)
     deallocate(rho,stat=i_stat)
     call memocc(i_stat,i_all,'rho',subname)
  end if
  i_all=-product(shape(deltarho))*kind(deltarho)
  deallocate(deltarho,stat=i_stat)
  call memocc(i_stat,i_all,'deltarho',subname)

end subroutine NK_SIC_potential

!> Squares the wavefunction in the real space for building the density.
!! Accumulate the result in rhoi, which has to be previously initialized
!! It squares complex and spinorial wavefunctions in the proper way
!! @ param fi occupation number times k-point weigth divided by the volume unit
!! @ param spinval value of the spin of the psir 
subroutine psir_to_rhoi(fi,spinval,nspinrho,nspinor,lr,psir,rhoi)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: nspinrho,nspinor
  real(gp), intent(in) :: fi,spinval
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(in) :: psir
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinrho), intent(inout) :: rhoi
  !local variables
  integer :: ncomplex,icomplex,npsir
  integer, dimension(4) :: nscarr_fake
  
  !quick return if fi=0
  if (fi==0.0_gp) return

  !scatterarr fake array
  nscarr_fake(1)=lr%d%n3i
  nscarr_fake(2)=lr%d%n3i
  nscarr_fake(3)=0
  nscarr_fake(4)=0

  if (nspinor ==4) then
     npsir=4
     ncomplex=0
  else
     npsir=1
     ncomplex=nspinor-1
  end if
  
  do icomplex=0,ncomplex
     select case(lr%geocode)
     case('F')
        call partial_density_free(.false.,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             npsir,nspinrho,lr%d%n3i,&
             fi,nscarr_fake,spinval,psir(1,icomplex*npsir+1),rhoi,lr%bounds%ibyyzz_r)
     case('P')
        call partial_density(.false.,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             npsir,nspinrho,lr%d%n3i,&
             fi,nscarr_fake,spinval,psir(1,icomplex*npsir+1),rhoi)
     case('S')
        call partial_density(.false.,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             npsir,nspinrho,lr%d%n3i,&
             fi,nscarr_fake,spinval,psir(1,icomplex*npsir+1),rhoi)
     end select
  end do

end subroutine psir_to_rhoi

