!!****p* BigDFT/rism
!!
!! COPYRIGHT
!!    Copyright (C) 2010 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
program rism
  use BigDFT_API
  implicit none
  character(len=*), parameter :: subname='rism'
  integer :: n1i,n2i,n3i,iproc,nproc,i_stat,i_all,nelec
  integer :: n3d,n3p,n3pi,i3xcsh,i3s
  real(gp) :: hxh,hyh,hzh
  type(atoms_data) :: atoms
  type(input_variables) :: in
  type(orbitals_data) :: orbs
  type(locreg_descriptors) :: Glr
  real(gp), dimension(3) :: shift
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(gp), dimension(:), allocatable :: atchgs,radii
  real(gp), dimension(:,:), allocatable :: radii_cf
  real(gp), dimension(:,:), pointer :: rxyz
  real(dp), dimension(:,:), pointer :: rho,pot,pot_ion

  !for the moment no need to have parallelism
  iproc=0
  nproc=1

  !initalise the varaibles for the calculation
  !standard names
  call standard_inputfile_names(in)
  call read_input_variables(iproc,'posinp',in,atoms,rxyz)

  if (iproc == 0) then
     call print_general_parameters(in,atoms)
  end if
       
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  call system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,atoms,rxyz,radii_cf,in%crmult,in%frmult,in%hx,in%hy,in%hz,&
       Glr,shift)

  call createWavefunctionsDescriptors(iproc,in%hx,in%hy,in%hz,&
       atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr)

  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)

  call createDensPotDescriptors(iproc,nproc,atoms%geocode,'D',&
       Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  !read the files from the .cubes on the disk
  call read_cube('electronic_density',atoms%geocode,&
       n1i,n2i,n3i,1,hxh,hyh,hzh,rho)

  call read_cube('hartree_potential',atoms%geocode,&
       n1i,n2i,n3i,1,hxh,hyh,hzh,pot)

  !not needed for the moment
  call read_cube('ionic_potential',atoms%geocode,&
       n1i,n2i,n3i,1,hxh,hyh,hzh,pot_ion)

  !also extract the number of atoms and the positions
  !for the atomic densities, if they are present
  

  !start the atomic charges calculation
  allocate(atchgs(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atchgs,'atchgs',subname)

  allocate(radii(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,radii,'radii',subname)

  call assign_atomic_radii(atoms,radii)

  !calculation of the atomic charges to compare wrt Mulliken
  !radii to be defined
  !call atomic_charges(iproc,nproc,rxyz,radii,atoms,0,Glr,ngatherarr,&
  call atomic_charges(iproc,nproc,rxyz,radii,atoms,nelec,Glr,ngatherarr,&
       !hxh,hyh,hzh,n3p,i3s+i3xcsh,rho,pot_ion,atchgs)
       hxh,hyh,hzh,n3p,i3s+i3xcsh,rho,pot,atchgs)


  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)
  i_all=-product(shape(ngatherarr))*kind(ngatherarr)
  deallocate(ngatherarr,stat=i_stat)
  call memocc(i_stat,i_all,'ngatherarr',subname)

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)


  call deallocate_lr(Glr,subname)
  
  call deallocate_orbs(orbs,subname)
  call deallocate_atoms_scf(atoms,subname) 

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)
  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
  i_all=-product(shape(radii))*kind(radii)
  deallocate(radii,stat=i_stat)
  call memocc(i_stat,i_all,'radii',subname)
  i_all=-product(shape(atchgs))*kind(atchgs)
  deallocate(atchgs,stat=i_stat)
  call memocc(i_stat,i_all,'atchgs',subname)
  
  call deallocate_atoms(atoms,subname) 
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  call free_input_variables(in)

  !finalize memory counting
  call memocc(0,0,'count','stop')


end program rism
!!***


!!****f* BigDFT/two_center_two_electrons
!! FUNCTION
!! COPYRIGHT
!!   Copyright (C) 2010 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine two_center_two_electrons(nat,a1,a2,a3,rxyz,radii,H)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), intent(in) :: a1,a2,a3
  real(gp), dimension(nat), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(gp), dimension(nat,nat), intent(out) :: H
  !local variables
  integer, parameter :: n_gauss = 89
  real(gp), parameter :: oneosqrtpi=0.564189583547756286948079451_gp
  integer :: iat,jat,i_gauss
  real(dp) :: ur_gauss,dr_gauss,acc_gauss,factor,factor2
  real(gp) :: a_range,ra2,ra2pb2,rab2,oneopk,oneoexpo,expo,oneofac,fac,ra
  real(gp), dimension(3) :: A
  real(dp), dimension(n_gauss) :: p_gauss,w_gauss

  !first calculate the diagonal elements of the matrix (analytic)
  do iat=1,nat
     ra=radii(iat)
     H(iat,iat)=-oneosqrtpi/ra
  end do

  !then the off-diagonal elements (upper part)

  !Initialization of the gaussian (Beylkin)
  call gequad(p_gauss,w_gauss,ur_gauss,dr_gauss,acc_gauss)
  !In order to have a range from a_range=sqrt(a1*a1+a2*a2+a3*a3)
  !(biggest length in the cube)
  !We divide the p_gauss by a_range**2 and a_gauss by a_range
  a_range = sqrt(a1*a1+a2*a2+a3*a3)
  factor = 1._dp/a_range
  !factor2 = factor*factor
  factor2 = 1._dp/(a1*a1+a2*a2+a3*a3)
  do i_gauss=1,n_gauss
     p_gauss(i_gauss) = factor2*p_gauss(i_gauss)
  end do
  do i_gauss=1,n_gauss
     w_gauss(i_gauss) = factor*w_gauss(i_gauss)
  end do

  do iat=1,nat
     ra2=radii(iat)**2
     A(1)=rxyz(1,iat)
     A(2)=rxyz(2,iat)
     A(3)=rxyz(3,iat)
     do jat=iat,nat
        ra2pb2=ra2+radii(jat)**2
        rab2=(rxyz(1,jat)-A(1))**2
        rab2=rab2+(rxyz(2,jat)-A(2))**2
        rab2=rab2+(rxyz(3,jat)-A(3))**2
        H(iat,jat)=0.0_gp
        !the correct order in which to do the sum seems to be the normal one
        do i_gauss=1,n_gauss
           oneopk=1.0_gp/p_gauss(i_gauss)
           oneoexpo=2.0_gp*ra2pb2+oneopk
           expo=rab2/oneoexpo
           expo=exp(-expo)
           oneofac=2.0_gp*ra2pb2*p_gauss(i_gauss)+1.0_gp
           oneofac=sqrt(oneofac)
           oneofac=oneofac**3
           fac=w_gauss(i_gauss)/oneofac
           H(iat,jat)=H(iat,jat)-fac*expo
        end do
        !calculate the analytic results via the Boys function
        !print *,'Analytic-approx:',iat,jat,H(iat,jat)+erfor(sqrt(rab2),sqrt(ra2pb2))
     end do
  end do

  !copy the same values on the lower triangular part
  do iat=1,nat
     do jat=iat+1,nat
        H(jat,iat)=H(iat,jat)
     end do
  end do

END SUBROUTINE two_center_two_electrons
!!***


!!****f* BigDFT/calculate_rho
!! FUNCTION
!!
!! SOURCE
!!
subroutine calculate_rho(iproc,nproc,geocode,nat,radii,rxyz,hxh,hyh,hzh,&
     n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pot,rhoarr)
  use module_base
  use module_types
  implicit none
  !Arguments---------
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,nat
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(nat), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(dp), dimension(n1i*n2i*n3pi), intent(in) :: pot
  real(gp), dimension(nat), intent(out) :: rhoarr
  !Local variables---------
  logical :: perx,pery,perz,gox,goy,goz
  real(gp) :: pi,prefactor,cutoff,rloc,Vel,cutofrac
  real(gp) :: rx,ry,rz,x,y,z,arg,r2,xp,charge
  integer :: i1,i2,i3,ind,iat,ierr
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,j1,j2,j3,isx,isy,isz,iex,iey,iez
  
  !experimental, just to verify the effect of the cutoff
  open(11)
  read(11,*)prefactor,cutofrac
  close(11)

  
  pi=4.d0*atan(1.d0)

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  charge=0.0_gp
  do iat=1,nat

     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)

     !print '(a,i0,3(1pe15.7))','coords',iat,rx,ry,rz
     !charge=0.0_gp

     !local part
     rloc=radii(iat)
     prefactor=1.0_gp/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
     !maximum extension of the gaussian
     !cutoff=10.d0*rloc
     cutoff=cutofrac!*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)
     
     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     !calculate the forces near the atom due to the error function part of the potential
     !calculate forces for all atoms only in the distributed part of the simulation box
     rhoarr(iat)=0.0_gp
     if (n3pi >0 ) then
        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    Vel=pot(ind)
                    rhoarr(iat)=rhoarr(iat)+Vel*xp
                    charge=charge+xp*prefactor
                    !if (xp > 1e-3_gp) write(16+iat-1,'(3(1x,i6),5(1x,1pe15.7))')j1+nbl1+1,j2+nbl2+1,j3,Vel,xp
                 end if
              end do
           end do
        end do
     end if

     !final results
     rhoarr(iat)=(hxh*hyh*hzh*prefactor)*rhoarr(iat)

  end do
  charge=hxh*hyh*hzh*charge
  !reduce the results
  if (nproc > 1) then
     call mpiallred(rhoarr(1),nat,MPI_SUM,MPI_COMM_WORLD,ierr)
     !the same should be done for the charge
     call mpiallred(charge,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  if (iproc == 0) write(*,*)'Charge:',charge

END SUBROUTINE calculate_rho
!!***

!!****f* BigDFT/calculate_rho_longrange
!! FUNCTION
!!
!! SOURCE
!!
subroutine calculate_rho_longrange(iproc,nproc,at,radii,rxyz,hxh,hyh,hzh,&
     n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,rho,rhoarr)
  use module_base
  use module_types
  implicit none
  !Arguments---------
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  real(gp), dimension(at%nat), intent(in) :: radii
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(n1i*n2i*n3pi), intent(inout) :: rho
  real(gp), dimension(at%nat), intent(out) :: rhoarr
  !Local variables---------
  logical :: perx,pery,perz,gox,goy,goz
  integer :: i1,i2,i3,ind,iat,ierr,ityp
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,j1,j2,j3,isx,isy,isz,iex,iey,iez
  real(gp) :: pi,prefactor,cutoff,rloc,Rel,cutofrac
  real(gp) :: rx,ry,rz,x,y,z,arg,r2,xp,charge,erfor

  !experimental, just to verify the effect of the cutoff
  open(11)
  read(11,*)prefactor,cutofrac
  close(11)
  
  pi=4.d0*atan(1.d0)

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  !substract the local density to the input density
  if (n3pi >0 .and. .false.) then
     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        cutoff=10.0_gp*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        !always calculate on the whole box
        do i3=isz,iez!1,n3pi
           !j3=i3+i3s-nbl3-1
           z=real(i3,gp)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey!1,n2i
              !j2=i2-nbl2-1
              y=real(i2,gp)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=isx,iex!1,n1i
                 !j1=i1-nbl1-1
                 x=real(i1,gp)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5_gp*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                 !if (r2 <= cutoff**2 ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    !ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    rho(ind)=rho(ind)-xp*charge
                 endif
              enddo
           enddo
        enddo
     end do
  end if

  charge=0.0_gp
  do iat=1,at%nat

     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)


     rloc=radii(iat)

     !cutoff=10.d0*rloc
     cutoff=cutofrac!*rloc

     rhoarr(iat)=0.0_gp
     if (n3pi >0 ) then
        do i3=-nbl3,2*n3+1+nbr3
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=-nbl2,2*n2+1+nbr2
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=-nbl1,2*n1+1+nbr1
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox .and.&
                     r2 > cutoff**2 ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    Rel=rho(ind)
                    rhoarr(iat)=rhoarr(iat)+Rel*erfor(sqrt(r2),rloc)
                    charge=charge+erfor(sqrt(r2),rloc)
                 end if
              end do
           end do
        end do
     end if

     !final results
     rhoarr(iat)=(hxh*hyh*hzh)*rhoarr(iat)

  end do
  charge=hxh*hyh*hzh*charge
  !reduce the results
  if (nproc > 1) then
     call mpiallred(rhoarr(1),at%nat,MPI_SUM,MPI_COMM_WORLD,ierr)
     !the same should be done for the charge
     call mpiallred(charge,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  if (iproc == 0) write(*,*)'Charge:',charge

END SUBROUTINE calculate_rho_longrange
!!***

!erf(r/(sqrt(2)rl)/r
function erfor(r,rl)
  use module_base
  implicit none
  real(gp), intent(in) :: r,rl
  real(gp) :: erfor
  !local variables
  real(gp) :: sq2rl,pi,derf_val

  pi=4.0_gp*atan(1.0_gp)

  sq2rl=sqrt(2.0_gp)*rl
  
  if (r == 0.0_gp) then
     erfor=2.0_gp/(sqrt(pi)*sq2rl)
  else
     call derf(derf_val,r/sq2rl)
     erfor=derf_val/r
  end if

end function erfor


!!****f* BigDFT/atomic_charges
!! FUNCTION
!!   calculate atomic charges using Lee, York and Yang method 
!!   Ref: J.Chem.Phys. 102(19),7549 (1995)
!!   use a basis of error functions centered on the atoms, with atom-defined radii
!!
!! SOURCE
!!
subroutine atomic_charges(iproc,nproc,rxyz,radii,atoms,nelec,lr,ngatherarr,&
     hxh,hyh,hzh,n3p,i3s,rho,pot,C)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,n3p,i3s,nelec
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(atoms_data), intent(in) :: atoms
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p), intent(in) :: pot
  real(dp), dimension(lr%d%n1i*lr%d%n2i*n3p), intent(inout) :: rho
  real(gp), dimension(atoms%nat) :: radii
  real(gp), dimension(3,atoms%nat) :: rxyz
  real(gp), dimension(atoms%nat) :: C
  !local variables
  character(len=*), parameter :: subname='atomic_charges'
  logical, parameter :: higherorder=.true.
  integer :: iat,info,nwork,i_all,i_stat,nbasis,i,j
  real(gp) :: ddotu,ddotv,gammafac
  type(gaussian_basis) :: Gpswf,Glongr
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: D,rhoarr,u,v,work,Caux
  real(gp), dimension(:,:), allocatable :: Hlrlr,Hwork,H,ovrlp

  if (atoms%geocode /= 'F') then
     write(*,*)'ERROR: the atomic charges can be calculated only in isolcated BC!'
     stop
  end if

  if (higherorder) then
     !let us first calculate the structure for the basis functions
     !extract the gaussian basis from the pseudowavefunctions
     nullify(Gpswf%rxyz)

     !here we can put the gaussian with herimte polynomials in place of that
     !call gaussian_pswf_basis(31,.false.,iproc,1,atoms,rxyz,Gpswf,Gocc)
     !if (associated(Gocc)) then
     !   i_all=-product(shape(Gocc))*kind(Gocc)
     !   deallocate(Gocc,stat=i_stat)
     !   call memocc(i_stat,i_all,'Gocc',subname)
     !   nullify(Gocc)
     !end if

     call gaussian_hermite_basis(1,atoms%nat,radii,rxyz,Gpswf)     

     call gaussian_rism_basis(atoms%nat,radii,rxyz,Glongr)

  !print *,'nat',atoms%nat,Glongr%ncoeff

     !after having determined the atomic basis funcitons, calculate the number of basis elements
     !taking also into account the long-range part of the basis
     nbasis=atoms%nat+Gpswf%ncoeff
  else
     nbasis=atoms%nat
  end if

  !allocate all the arrays
  allocate(H(nbasis,nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,H,'H',subname)
  allocate(rhoarr(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,rhoarr,'rhoarr',subname)
  allocate(D(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,D,'D',subname)
  allocate(u(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,u,'u',subname)
  allocate(v(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,v,'v',subname)
  allocate(Caux(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,Caux,'Caux',subname)

  if (higherorder) then
     !calculate the overlap matrix as well as the kinetic overlap
     !in view of complete gaussian calculation
     allocate(ovrlp(Gpswf%ncoeff,Gpswf%ncoeff),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     
     !overlap calculation of the kinetic operator
     call kinetic_overlap_h(Gpswf,Gpswf,ovrlp)
  
     !fill the last part of the H matrix
     !for the same structure the kinetic overlap is symmetric
     do j=1,Gpswf%ncoeff
        do i=1,j-1
           H(atoms%nat+i,atoms%nat+j)=-1.0_gp/(8.0_gp*atan(1.0_gp))*ovrlp(i,j)
           H(atoms%nat+j,atoms%nat+i)=-1.0_gp/(8.0_gp*atan(1.0_gp))*ovrlp(i,j)
        end do
        !diagonal elements
        H(atoms%nat+j,atoms%nat+j)=-1.0_gp/(8.0_gp*atan(1.0_gp))*ovrlp(j,j)
     end do


     !test the overlap matrices
     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     !calculate the overlap matrix as well as the kinetic overlap
     !in view of complete gaussian calculation
!!$     allocate(ovrlp(Gpswf%ncoeff,Gpswf%ncoeff),stat=i_stat)
!!$     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!$     call gaussian_overlap(Gpswf,Gpswf,ovrlp)
!!$     if (iproc == 0) then
!!$        do iat=1,Gpswf%ncoeff
!!$           write(*,'(a,i0,10(1pe15.7))')'Gpswf',iat,ovrlp(1:iat,iat)
!!$        end do
!!$     end if
!!$     
!!$     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!$     deallocate(ovrlp,stat=i_stat)
!!$     call memocc(i_stat,i_all,'ovrlp',subname)
     !calculate the overlap matrix as well as the kinetic overlap
     !in view of complete gaussian calculation
     allocate(ovrlp(Glongr%ncoeff,Glongr%ncoeff),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     call gaussian_overlap(Glongr,Glongr,ovrlp)
     if (iproc == 0) then
        do iat=1,Glongr%ncoeff
           write(*,'(a,i0,10(1pe15.7))')'Glongr',iat,ovrlp(1:iat,iat)
        end do
     end if

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     allocate(ovrlp(Gpswf%ncoeff,Glongr%ncoeff),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     
     !overlap between longrange basis and short-range basis
     call gaussian_overlap_h(Gpswf,Glongr,ovrlp)
     
     !fill the block off-diagonal part of the H matrix
     print *,'coeffs',Gpswf%ncoeff,Glongr%ncoeff
     do j=1,Glongr%ncoeff
        do i=1,Gpswf%ncoeff
           if (Gpswf%ncoeff == Glongr%ncoeff .and. i > j) then
              H(atoms%nat+i,j)=ovrlp(j,i)
              H(j,atoms%nat+i)=ovrlp(j,i)
           else
              H(atoms%nat+i,j)=ovrlp(i,j)
              H(j,atoms%nat+i)=ovrlp(i,j)
           end if
        end do
     end do
     
     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
  end if

  !here nat coincides with the number of long range basis functions
  !calculate H matrix
  allocate(Hlrlr(atoms%nat,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,Hlrlr,'Hlrlr',subname)

  call two_center_two_electrons(atoms%nat,atoms%alat1,atoms%alat2,atoms%alat3,&
       rxyz,radii,Hlrlr)
 
  !fill the first part of the H matrix
  do j=1,atoms%nat
     do i=1,atoms%nat
        H(i,j)=Hlrlr(i,j)
     end do
  end do

  i_all=-product(shape(Hlrlr))*kind(Hlrlr)
  deallocate(Hlrlr,stat=i_stat)
  call memocc(i_stat,i_all,'Hlrlr',subname)

  if (iproc == 0) then
     do iat=1,nbasis
        write(*,'(a,i0,10(1pe15.7))')'H',iat,H(1:iat,iat)
     end do
  end if
 
  !calculate the first part of rho array
  call calculate_rho(iproc,nproc,atoms%geocode,atoms%nat,radii,rxyz,hxh,hyh,hzh,&
       lr%d%n1,lr%d%n2,lr%d%n3,n3p,i3s,lr%d%n1i,lr%d%n2i,lr%d%n3i,pot,rhoarr)

  if (iproc == 0) then
     do iat=1,atoms%nat
        write(*,'(a,i0,10(1pe15.7))')'rhoarrV',iat,rhoarr(iat)
     end do
  end if

  call calculate_rho_longrange(iproc,nproc,atoms,radii,rxyz,hxh,hyh,hzh,&
       lr%d%n1,lr%d%n2,lr%d%n3,n3p,i3s,lr%d%n1i,lr%d%n2i,lr%d%n3i,rho,rhoarr)

  if (iproc == 0) then
     do iat=1,atoms%nat
        write(*,'(a,i0,10(1pe15.7))')'rhoarrrho',iat,rhoarr(iat)
     end do
  end if

!!$  call MPI_FINALIZE(i_stat)
!!$  stop



!!$  if (iproc == 0) then
!!$        write(*,'(a,4(1pe15.7))')'rho',rho(:)
!!$  end if
  
  if (higherorder) then
     !calculate the shortrange part of the rho array
     call calculate_rho_shortrange(iproc,nproc,atoms,lr,Gpswf,hxh,hyh,hzh,rxyz,ngatherarr,&
          rho,rhoarr(atoms%nat+1))
  end if


  if (iproc == 0) then
     do iat=1,nbasis
        write(*,'(a,i0,10(1pe15.7))')'rhoarr',iat,rhoarr(iat)
     end do
  end if

  if (higherorder) then
     nullify(Gpswf%rxyz)
     call deallocate_gwf(Gpswf,subname)
     nullify(Glongr%rxyz)
     call deallocate_gwf(Glongr,subname)
  end if


  !initalise D array
  D=0.0_gp
  do iat=1,atoms%nat
     D(iat)=1.0_gp
  end do


  !calculate workspace array for dsysv
  !temporary allocation of the work array, workspace query in dsysv
  allocate(Hwork(nbasis,nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,Hwork,'Hwork',subname)
  allocate(iwork(nbasis+ndebug),stat=i_stat)
  call memocc(i_stat,iwork,'iwork',subname)
  allocate(work(100+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  call dcopy(nbasis*nbasis,H,1,Hwork,1)
  call dcopy(nbasis,rhoarr,1,u,1)
  !workspace query
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),u(1),nbasis,work(1),-1,info)
  nwork=work(1)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  !here we start with the linear algebra
  call dcopy(nbasis,rhoarr,1,u,1)
  call dcopy(nbasis*nbasis,H,1,Hwork,1)
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),u(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of u',info
  end if

  !here we start with the linear algebra
  call dcopy(nbasis,D,1,v,1)
  call dcopy(nbasis*nbasis,H,1,Hwork,1)
!!$  if (iproc == 0) print *,'here',v(:)
!!$  if (iproc == 0) print '(a,4(1pe15.7))','000',Hwork(1,1),Hwork(1,2),Hwork(2,1),Hwork(2,2)
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),v(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of v',info
  end if
!!$  if (iproc == 0) print '(a,4(1pe15.7))','there',v(:)
!!$  if (iproc == 0) print '(a,4(1pe15.7))','AAA',H(1,1),H(1,2),H(2,1),H(2,2)
!!$  !determinant of the matrix, temporary
!!$  gammafac=H(1,1)*H(2,2)-H(1,2)*H(2,1)
!!$  if (iproc == 0) print '(a,4(1pe15.7))','res',(H(2,2)-H(1,2))/gammafac,&
!!$       (H(1,1)-H(2,1))/gammafac,real(nelec,gp)
  ddotu=dot(nbasis,D(1),1,u(1),1)
  ddotv=dot(nbasis,D(1),1,v(1),1)
  

  gammafac=(real(nelec,gp)+ddotu)/ddotv

  !!zero has to be put when the potential is the deformation potential
  !gammafac=ddotu/ddotv

  !calculate the eigenvalues of H
  call dcopy(nbasis*nbasis,H,1,Hwork,1)
  !print *,'nwork',nwork,3*nbasis-1
  call dsyev('N','U',nbasis,Hwork(1,1),nbasis,v(1),work(1),nwork,info)

  if (iproc == 0) then
     do iat=1,nbasis
        write(*,'(a,i0,10(1pe15.7))')'eigenvalues',iat,v(iat)
     end do
     write(*,*)'Condition Number (valid only if all eigenvalues are negative):',v(1)/v(nbasis)
  end if



!!$  if (iproc == 0 )print *,'gamma',gammafac,nelec,v(:)

  !rescale D
  call vscal(nbasis,gammafac,D(1),1)

  !initalise C
  call dcopy(nbasis,D,1,Caux,1)
  
  !fill it 
  call axpy(nbasis,-1.0_gp,rhoarr(1),1,Caux(1),1)

  !solve the system
  call dsysv('U',nbasis,1,H(1,1),nbasis,iwork(1),Caux(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of charges',info
  end if

  !print the charges
  if (iproc == 0) then
     do iat=1,atoms%nat
        write(*,'(1x,a,i4,3(1x,f12.7))')'atom, charge',iat,Caux(iat),Caux(iat)-real(atoms%nelpsp(atoms%iatype(iat)),gp),radii(iat)
        C(iat)=Caux(iat)
     end do
  end if

  i_all=-product(shape(iwork))*kind(iwork)
  deallocate(iwork,stat=i_stat)
  call memocc(i_stat,i_all,'iwork',subname)
  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  i_all=-product(shape(Hwork))*kind(Hwork)
  deallocate(Hwork,stat=i_stat)
  call memocc(i_stat,i_all,'Hwork',subname)

  i_all=-product(shape(H))*kind(H)
  deallocate(H,stat=i_stat)
  call memocc(i_stat,i_all,'H',subname)
  i_all=-product(shape(rhoarr))*kind(rhoarr)
  deallocate(rhoarr,stat=i_stat)
  call memocc(i_stat,i_all,'rhoarr',subname)

  i_all=-product(shape(D))*kind(D)
  deallocate(D,stat=i_stat)
  call memocc(i_stat,i_all,'D',subname)

  i_all=-product(shape(u))*kind(u)
  deallocate(u,stat=i_stat)
  call memocc(i_stat,i_all,'u',subname)

  i_all=-product(shape(v))*kind(v)
  deallocate(v,stat=i_stat)
  call memocc(i_stat,i_all,'v',subname)

  i_all=-product(shape(Caux))*kind(Caux)
  deallocate(Caux,stat=i_stat)
  call memocc(i_stat,i_all,'Caux',subname)


END SUBROUTINE atomic_charges
!!***


subroutine assign_atomic_radii(at,radii)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(at%nat), intent(out) :: radii
  !local variables
  real(gp), parameter :: xi=1.1839527_gp
  integer :: iat,ityp
  real(gp) :: lambda,lambdafrac,cutoff

  open(11)
  read(11,*)lambdafrac,cutoff
  close(11)

  !take York's paper radii and convert them in the new basis
  do iat=1,at%nat
     ityp=at%iatype(iat)
     select case (at%atomnames(ityp))
     case ('H')
        lambda=2.66_gp
     case('Li')
        lambda=2.07_gp
     case('B')
        lambda=3.58_gp
     case('C')
        lambda=4.27_gp
     case('N')
        lambda=4.91_gp
     case('O')
        lambda=5.64_gp
     case('F')
        lambda=6.29_gp
     case('Cl')
        lambda=1.71_gp
     case default
        write(*,*)'ERROR: the radius is not yet defined for this atom'
        stop
     end select
     !calculate the radius by minimizing the difference between 
     !York's function and our basis
     !radii(iat)=0.5_gp*xi/lambda

     !otherwise use the local psp radii
     radii(iat)=lambdafrac*lambda!at%psppar(0,0,ityp)
     !radii(iat)=lambdafrac*at%psppar(0,0,ityp)
     !if (at%atomnames(ityp) == 'O') then
     !   radii(iat)=lambdafrac*at%psppar(0,0,ityp)
     !else
     !   radii(iat)=at%psppar(0,0,ityp)
     !end if
  end do
end subroutine assign_atomic_radii
  
!!****f* BigDFT/gaussian_rism_basis
!! FUNCTION
!!   Gaussian basis associated to the long-range term of rism calculation
!! SOURCE
!!
subroutine gaussian_rism_basis(nat,radii,rxyz,G)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(nat), intent(in) :: radii
  real(gp), dimension(3,nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: G  
  !local variables
  character(len=*), parameter :: subname='gaussian_psp_basis'
  real(gp), parameter :: oneo2pi3halves=0.0634936359342409697857633_gp
  integer :: iat,nshell,iexpo,l,ishell,i_stat

  G%nat=nat
  G%rxyz => rxyz
  allocate(G%nshell(G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)
 
  G%nshltot=0
  do iat=1,G%nat
     nshell=1
     G%nshell(iat)=nshell
     G%nshltot=G%nshltot+nshell
  end do

  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)

  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  do iat=1,G%nat
     do l=1,4 
        if (l==1) then
           ishell=ishell+1
           G%ndoc(ishell)=1
           G%nam(ishell)=l
           G%nexpo=G%nexpo+1
           G%ncoeff=G%ncoeff+2*l-1
        end if
     enddo
  end do

  !allocate and assign the exponents and the coefficients
  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  ishell=0
  iexpo=0
  do iat=1,G%nat
     do l=1,4 
        if (l==1) then
           ishell=ishell+1
           iexpo=iexpo+1
           G%psiat(iexpo)=-oneo2pi3halves/radii(iat)**3
           G%xp(iexpo)=radii(iat)
        end if
     end do
  end do

END SUBROUTINE gaussian_rism_basis
!!***

!!****f* BigDFT/gaussian_hermite_basis
!! FUNCTION
!!   Gaussian basis associated to Hermite Polyomials multiplied by s-terms
!!   given by the radii
!! SOURCE
!!
subroutine gaussian_hermite_basis(nhermitemax,nat,radii,rxyz,G)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nat,nhermitemax
  real(gp), dimension(nat), intent(in) :: radii
  real(gp), dimension(3,nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: G  
  !local variables
  character(len=*), parameter :: subname='gaussian_psp_basis'
  real(gp), parameter :: oneo2pi3halves=0.0634936359342409697857633_gp
  integer :: iat,nshell,iexpo,l,ishell,i_stat

  G%nat=nat
  G%rxyz => rxyz
  allocate(G%nshell(G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)
 
  G%nshltot=0
  do iat=1,G%nat
     nshell=nhermitemax
     G%nshell(iat)=nshell
     G%nshltot=G%nshltot+nshell
  end do

  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)

  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  do iat=1,G%nat
     do l=1,G%nshell(iat) 
        ishell=ishell+1
        G%ndoc(ishell)=1
        G%nam(ishell)=1
        G%nexpo=G%nexpo+1
        G%ncoeff=G%ncoeff+1
     enddo
  end do

  !allocate and assign the exponents and the coefficients
  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  ishell=0
  iexpo=0
  do iat=1,G%nat
     do l=1,G%nshell(iat) 
        ishell=ishell+1
        iexpo=iexpo+1
        G%psiat(iexpo)=1.0_gp
        G%xp(iexpo)=radii(iat)
     end do
  end do

end subroutine gaussian_hermite_basis
!!***


!calculate the second part, by expressing the atomic wavefunctions on a real grid
subroutine calculate_rho_shortrange(iproc,nproc,at,lr,Gpswf,hxh,hyh,hzh,rxyz,ngatherarr,&
     rho,rhoarr)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  type(gaussian_basis), intent(in) :: Gpswf
  type(locreg_descriptors), intent(in) :: lr
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(ngatherarr(iproc,1)), intent(inout) :: rho
  real(dp), dimension(Gpswf%ncoeff), intent(out) :: rhoarr
  !local variables
  character(len=*), parameter :: subname='calculate_rho_shortrange'
  logical, parameter :: deltarho=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,iat,isx,isy,isz,iex,iey,iez,n3pi,i3s
  integer :: i1,i2,i3,j1,j2,j3,ityp,ind,i_stat,i_all,ierr,jorb,isorb,jproc,kproc
  real(gp) :: pi,rx,ry,rz,charge,rloc,cutoff,x,y,z,r2,arg,xp,hfac
  type(orbitals_data) :: orbspswf
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncoeff_par
  real(wp), dimension(:,:), allocatable :: psigaup,psi
  real(wp), dimension(:,:,:), allocatable :: psir
  real(dp), dimension(:,:,:), allocatable :: rhotot

  !here we have to add the density the compensation value and gather the complete result
  !----extracted from createIonicPotential
  pi=4.0_gp*atan(1.0_gp)
  ! Ionic energy (can be calculated for all the processors)

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  n3pi=ngatherarr(iproc,1)/(lr%d%n1i*lr%d%n2i)
  i3s=ngatherarr(iproc,2)/(lr%d%n1i*lr%d%n2i)

  if (n3pi >0 .and. deltarho) then

     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        cutoff=10.0_gp*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        do i3=isz,iez
           z=real(i3,gp)*hzh-rz
           call ind_positions(perz,i3,lr%d%n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,gp)*hyh-ry
              call ind_positions(pery,i2,lr%d%n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,gp)*hxh-rx
                 call ind_positions(perx,i1,lr%d%n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5_gp*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*lr%d%n1i+(j3-i3s+1-1)*lr%d%n1i*lr%d%n2i
                    rho(ind)=rho(ind)-xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if
  !-----end of extraction

  !the gathering should be done here
  allocate(rhotot(lr%d%n1i,lr%d%n2i,lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,rhotot,'rhotot',subname)

  if (nproc > 1) then
     call MPI_ALLGATHERV(rho,ngatherarr(iproc,1),&
          mpidtypw,rhotot,ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
  else
     !not efficient, the copy can be saved
     call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i,rho,1,rhotot,1)
  end if

  allocate(ncoeff_par(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ncoeff_par,'ncoeff_par',subname)


  !determine the number of basis functions treated by each processor
  ncoeff_par(0:nproc-1,1)=0
  do jorb=1,Gpswf%ncoeff
     jproc=mod(jorb-1,nproc)
     ncoeff_par(jproc,1)=ncoeff_par(jproc,1)+1
  end do
  !starting point
  do jproc=0,nproc-1
     isorb=0
     do kproc=0,jproc-1
        isorb=isorb+ncoeff_par(kproc,1)
     end do
     ncoeff_par(jproc,2)=isorb
  end do
  isorb=ncoeff_par(iproc,2)

  !fake orbitals descriptor to calculate the wavelet expansion
  call orbitals_descriptors(0,1,ncoeff_par(iproc,1),ncoeff_par(iproc,1),0,1,1,&
       (/0.0_gp,0.0_gp,0.0_gp/),(/1.0_gp /),orbspswf)

  !allocate the wavefunctions in wavelets and in gaussians
  allocate(psi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncoeff_par(iproc,1)+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  allocate(psigaup(Gpswf%ncoeff,ncoeff_par(iproc,1)+ndebug),stat=i_stat)
  call memocc(i_stat,psigaup,'psigaup',subname)

  !the coefficients are elements of the identity matrix
  call razero(Gpswf%ncoeff*ncoeff_par(iproc,1),psigaup)
  do jorb=1,ncoeff_par(iproc,1)
     psigaup(isorb+jorb,jorb)=1.0_gp
  end do

  !calculate the atomic wavefunctions in Daubechies basis
  !put nproc=1 since the orbitals descriptors is not parallelised by choice
  call gaussians_to_wavelets_new_h(0,1,lr,orbspswf,2.0_gp*hxh,2.0_gp*hyh,2.0_gp*hzh,Gpswf,&
       psigaup,psi)

  call deallocate_orbs(orbspswf,subname)

  !initalize work arrays for decompressing the wavefunction
  call initialize_work_arrays_sumrho(lr,w)

  allocate(psir(lr%d%n1i,lr%d%n2i,lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if

  do jorb=1,ncoeff_par(iproc,1)

     hfac=(hxh*hyh*hzh)
     call daub_to_isf(lr,w,psi(1,jorb),psir)

     !now the integral with the density
     !the square root of the volume units should be put because the psi is normalised to 
     !have unit norm without volume unit
     rhoarr(jorb+isorb)=sqrt(hfac)*dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,rhotot(1,1,1),1,psir(1,1,1),1)
     !check the norm of the wavefunctions
!!$     print *,'norm',jorb+isorb,&
!!$          dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir(1,1,1),1,psir(1,1,1),1),&
!!$          rhoarr(jorb+isorb),&
!!$          sum(rhotot),hfac

     
  end do

  !gather the results in the same array
  if (nproc > 1) then
     call MPI_ALLGATHERV(rhoarr(isorb+min(1,ncoeff_par(iproc,1))),ncoeff_par(iproc,1),&
          mpidtypw,rhoarr(1),ncoeff_par(0,1),&
          ncoeff_par(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
  end if


  !deallocate the deallocatable
  call deallocate_work_arrays_sumrho(w)

  i_all=-product(shape(rhotot))*kind(rhotot)
  deallocate(rhotot,stat=i_stat)
  call memocc(i_stat,i_all,'rhotot',subname)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

  i_all=-product(shape(psigaup))*kind(psigaup)
  deallocate(psigaup,stat=i_stat)
  call memocc(i_stat,i_all,'psigaup',subname)

  i_all=-product(shape(ncoeff_par))*kind(ncoeff_par)
  deallocate(ncoeff_par,stat=i_stat)
  call memocc(i_stat,i_all,'ncoeff_par',subname)
  
end subroutine calculate_rho_shortrange
