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
  real(dp) :: ur_gauss,dr_gauss,acc_gauss,pgauss,factor,factor2
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
     do jat=iat+1,nat
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
     end do
  end do

  !copy the same values on the lower triangular part
  do iat=1,nat
     do jat=iat+1,nat
        H(jat,iat)=H(iat,jat)
     end do
  end do

end subroutine two_center_two_electrons

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
  real(gp) :: pi,prefactor,cutoff,rloc,Vel,rhoel
  real(gp) :: rx,ry,rz,x,y,z,arg,r2,xp,tt,charge
  integer :: i1,i2,i3,ind,iat,ityp,nloc,iloc,ierr
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,j1,j2,j3,isx,isy,isz,iex,iey,iez
  
  pi=4.d0*atan(1.d0)


  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  do iat=1,nat
     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)

     !local part
     rloc=radii(iat)
     prefactor=1.0_gp/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
     !maximum extension of the gaussian
     cutoff=10.d0*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)
     
     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     !calculate the forces near the atom due to the error function part of the potential
     !calculate forces for all atoms only in the distributed part of the simulation box
     charge=0.0_gp
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
                    charge=charge*xp
                 end if
              end do
           end do
        end do
     end if
     
     !final results
     rhoarr(iat)=(hxh*hyh*hzh*prefactor)*rhoarr(iat)
     charge=(hxh*hyh*hzh*prefactor)*charge

  end do

  !reduce the results
  if (nproc > 1) then
     call mpiallred(rhoarr(1),nat,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  !the same should be done for the charge

end subroutine calculate_rho
!!***

!calculate atomic charges using Lee, York and Yang method 
!Ref: J.Chem.Phys. 102(19),7549 (1995)
!use a basis of error functions centered on the atoms, with atom-defined radii
subroutine atomic_charges(iproc,nproc,geocode,rxyz,radii,alat1,alat2,alat3,nelec,nat,grid,&
     hxh,hyh,hzh,n3p,i3s,pot,C)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,nat,n3p,i3s,nelec
  real(gp), intent(in) :: hxh,hyh,hzh,alat1,alat2,alat3
  type(grid_dimensions), intent(in) :: grid
  real(wp), dimension(grid%n1i*grid%n2i*n3p), intent(in) :: pot
  real(gp), dimension(nat) :: radii
  real(gp), dimension(3,nat) :: rxyz
  real(gp), dimension(nat) :: C
  !local variables
  character(len=*), parameter :: subname='atomic_charges'
  integer :: iat,info,nwork,i_all,i_stat
  real(gp) :: ddotu,ddotv,gammafac
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: D,rho,u,v,work
  real(gp), dimension(:,:), allocatable :: H,Hwork

  if (geocode /= 'F') then
     write(*,*)'ERROR: the atomic charges can be calculated only in isolcated BC!'
     stop
  end if

  !here nat coincides with the number of basis functions

  !allocate all the arrays
  allocate(H(nat,nat+ndebug),stat=i_stat)
  call memocc(i_stat,H,'H',subname)
  allocate(rho(nat+ndebug),stat=i_stat)
  call memocc(i_stat,rho,'rho',subname)
  allocate(D(nat+ndebug),stat=i_stat)
  call memocc(i_stat,D,'D',subname)
  allocate(u(nat+ndebug),stat=i_stat)
  call memocc(i_stat,u,'u',subname)
  allocate(v(nat+ndebug),stat=i_stat)
  call memocc(i_stat,v,'v',subname)


  !calculate H matrix
  call two_center_two_electrons(nat,alat1,alat2,alat3,rxyz,radii,H)

  !calculate rho array
  call calculate_rho(iproc,nproc,geocode,nat,radii,rxyz,hxh,hyh,hzh,&
       grid%n1,grid%n2,grid%n3,n3p,i3s,grid%n1i,grid%n2i,grid%n3i,pot,rho)

  !initalise D array
  do iat=1,nat
     D(iat)=1.0_gp
  end do

  !calculate workspace array for dsysv
  !temporary allocation of the work array, workspace query in dsysv
  allocate(Hwork(nat,nat+ndebug),stat=i_stat)
  call memocc(i_stat,Hwork,'Hwork',subname)
  allocate(iwork(nat+ndebug),stat=i_stat)
  call memocc(i_stat,iwork,'iwork',subname)
  allocate(work(100+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  call dcopy(nat*nat,H,1,Hwork,1)
  call dcopy(nat,rho,1,u,1)
  call dsysv('U',nat,1,Hwork(1,1),nat,iwork(1),u(1),nat,work(1),-1,info)
  nwork=work(1)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  !here we start with the linear algebra
  call dcopy(nat,rho,1,u,1)
  call dcopy(nat*nat,H,1,Hwork,1)
  call dsysv('U',nat,1,Hwork(1,1),nat,iwork(1),u(1),nat,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of u',info
  end if

  !here we start with the linear algebra
  call dcopy(nat,D,1,v,1)
  call dcopy(nat*nat,H,1,Hwork,1)
  call dsysv('U',nat,1,Hwork(1,1),nat,iwork(1),v(1),nat,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of v',info
  end if

  ddotu=dot(nat,D(1),1,u(1),1)
  ddotv=dot(nat,D(1),1,v(1),1)

  gammafac=(real(nelec,gp)+ddotu)/ddotv

  !rescale D
  call vscal(nat,gammafac,D(1),1)

  !initalise C
  call dcopy(nat,D,1,C,1)
  
  !fill it 
  call axpy(nat,1.0_gp,rho(1),1,C(1),1)

  !solve the system
  call dsysv('U',nat,1,H(1,1),nat,iwork(1),C(1),nat,work(1),nwork,info)


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
  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)
  i_all=-product(shape(D))*kind(D)
  deallocate(D,stat=i_stat)
  call memocc(i_stat,i_all,'D',subname)
  i_all=-product(shape(u))*kind(u)
  deallocate(u,stat=i_stat)
  call memocc(i_stat,i_all,'u',subname)
  i_all=-product(shape(v))*kind(v)
  deallocate(v,stat=i_stat)
  call memocc(i_stat,i_all,'v',subname)


end subroutine atomic_charges
