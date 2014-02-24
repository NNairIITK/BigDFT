!> @file
!!   In this file, we have the analytic routines for the calculation of the overlap of short-range functions
!! @author
!!   Copyright (C) 2010-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


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


!>   calculate atomic charges using Lee, York and Yang method 
!!   Ref: J.Chem.Phys. 102(19),7549 (1995)
!!   use a basis of error functions centered on the atoms, with atom-defined radii
subroutine atomic_charges_york(iproc,nproc,rxyz,radii,atoms,nelec,lr,ngatherarr,&
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
     write(*,*)'ERROR: the atomic charges can be calculated only in isolated BC!'
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
       lr%d%n1,lr%d%n2,lr%d%n3,n3p,i3s,lr%d%n1i,lr%d%n2i,pot,rhoarr)

  if (iproc == 0) then
     do iat=1,atoms%nat
        write(*,'(a,i0,10(1pe15.7))')'rhoarrV',iat,rhoarr(iat)
     end do
  end if

  call calculate_rho_longrange(iproc,nproc,atoms,radii,rxyz,hxh,hyh,hzh,&
       lr%d%n1,lr%d%n2,lr%d%n3,n3p,i3s,lr%d%n1i,lr%d%n2i,rho,rhoarr)

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

  call vcopy(nbasis*nbasis,H,1,Hwork,1)
  call vcopy(nbasis,rhoarr,1,u,1)
  !workspace query
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),u(1),nbasis,work(1),-1,info)
  nwork=work(1)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  !here we start with the linear algebra
  call vcopy(nbasis,rhoarr,1,u,1)
  call vcopy(nbasis*nbasis,H,1,Hwork,1)
  call dsysv('U',nbasis,1,Hwork(1,1),nbasis,iwork(1),u(1),nbasis,work(1),nwork,info)
  if (info /= 0) then
     write(*,*) 'info calculation of u',info
  end if

  !here we start with the linear algebra
  call vcopy(nbasis,D,1,v,1)
  call vcopy(nbasis*nbasis,H,1,Hwork,1)
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
  call vcopy(nbasis*nbasis,H,1,Hwork,1)
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
  call vcopy(nbasis,D,1,Caux,1)
  
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


END SUBROUTINE atomic_charges_york


subroutine assign_atomic_radii_york(at,radii)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(at%nat), intent(out) :: radii
  !local variables
  !n(c) real(gp), parameter :: xi=1.1839527_gp
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
END SUBROUTINE assign_atomic_radii_york


!>   Gaussian basis associated to Hermite Polyomials multiplied by s-terms
!!   given by the radii
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
  !n(c) real(gp), parameter :: oneo2pi3halves=0.0634936359342409697857633_gp
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
  G%ncplx=1
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
  allocate(G%xp(G%ncplx,G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)
  allocate(G%psiat(G%ncplx,G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  ishell=0
  iexpo=0
  do iat=1,G%nat
     do l=1,G%nshell(iat) 
        ishell=ishell+1
        iexpo=iexpo+1
        G%psiat(1,iexpo)=1.0_gp
        G%xp(1,iexpo)=radii(iat)
     end do
  end do

END SUBROUTINE gaussian_hermite_basis


subroutine calculate_rho(iproc,nproc,geocode,nat,radii,rxyz,hxh,hyh,hzh,&
     n1,n2,n3,n3pi,i3s,n1i,n2i,pot,rhoarr)
  use module_base
  use module_types
  implicit none
  !Arguments---------
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,nat
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


!> Overlap kinetic matrix between two different basis structures
!! the kinetic operator is applicated on the A basis structure
!! The basis structure is supposed to be based on s-functions times Hermite polynomials
subroutine kinetic_overlap_h(A,B,ovrlp)
  use module_base
  use module_types
  implicit none
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables
  integer, parameter :: niw=304,nrw=304
  integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
  integer :: iovrlp,jovrlp,jcoeff,jexpo
  integer :: ngA,ngB,lA,lB,mA,mB
  real(gp) :: dx,dy,dz
  integer, dimension(niw) :: iw
  real(gp), dimension(nrw) :: rw

  iovrlp=0
  ishell=0
  iexpo=1
  icoeff=1
  !loop on each shell (intensive calculation)
  do iat=1,A%nat
     do isat=1,A%nshell(iat)
        ishell=ishell+1
        ngA=A%ndoc(ishell)
        lA=A%nam(ishell)
        if (lA /= 1 ) then
           stop 'only s case is allowed for Hermite polynomials, A' 
        end if
        do mA=1,2*lA-1
           iovrlp=iovrlp+1

           jovrlp=0
           jshell=0
           jexpo=1
           jcoeff=1

           do jat=1,B%nat
              dx=B%rxyz(1,jat)-A%rxyz(1,iat)
              dy=B%rxyz(2,jat)-A%rxyz(2,iat)
              dz=B%rxyz(3,jat)-A%rxyz(3,iat)
              do jsat=1,B%nshell(jat)
                 jshell=jshell+1
                 ngB=B%ndoc(jshell)
                 lB=B%nam(jshell)
                 if (lA /= 1 ) then
                    stop 'only s case is allowed for Hermite polynomials, B' 
                 end if
                 do mB=1,2*lB-1
                    jovrlp=jovrlp+1
                    if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff .or. &
                         A%ncoeff /= B%ncoeff ) then
                       call kineticovrlp_h(A%xp(1,iexpo),A%psiat(1,iexpo),&
                            B%xp(1,jexpo),B%psiat(1,jexpo),&
                            ngA,ngB,lA,isat,lB,jsat,dx,dy,dz,&
                            niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                    end if
                 end do
                 jexpo=jexpo+ngB
                 jcoeff=jcoeff+2*lB-1
              end do
           end do
        end do
        iexpo=iexpo+ngA
        icoeff=icoeff+2*lA-1
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
  call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)
  
END SUBROUTINE kinetic_overlap_h


!>   Calculates the scalar product between two shells
!!   by considering only the nonzero coefficients
!!   actual building block for calculating overlap matrix
!!   inserted work arrays for calculation
!!   Only Hermite polynomials of r^2 are to be considered
subroutine kineticovrlp_h(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,ih1,l2,ih2,dx,dy,dz,&
     niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: ng1,ng2,l1,ih1,l2,ih2,niw,nrw
  real(gp), intent(in) :: dx,dy,dz
  integer, dimension(niw) :: iw
  real(gp), dimension(nrw) :: rw
  real(gp), dimension(ng1), intent(in) :: expo1,coeff1
  real(gp), dimension(ng2), intent(in) :: expo2,coeff2
  real(gp), intent(out) :: ovrlp
  !local variables
  integer :: i1,i2
  real(gp) :: a1,a2,c1,c2,govrlpr

  ovrlp=0.d0
  do i1=1,ng1
     a1=expo1(i1)
     a1=0.5_gp/a1**2
     c1=coeff1(i1)
     do i2=1,ng2
        a2=expo2(i2)
        a2=0.5_gp/a2**2
        c2=coeff2(i2)
        call kinprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,ih2,niw,nrw,iw,rw,govrlpr)
        govrlpr=c1*govrlpr*c2
        !print *,c1,c2,govrlpr
        ovrlp=ovrlp+govrlpr
     end do
  end do
  
END SUBROUTINE kineticovrlp_h


!>  Kinetic overlap between gaussians, based on cartesian coordinates
!!  calculates a dot product between two differents gaussians times spherical harmonics
!!  only hermite polynomials
subroutine kinprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,ih2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,ih1,ih2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=48
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
  real(gp) :: fx,fy,fz,fa,fb,govrlp,kinovrlp,d2fx,d2fy,d2fz

  !calculate the coefficients of the hermite polynomials
  !calculates the number of different couples
  call calc_coeff_hermite_r2(l1,ih1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
  call calc_coeff_hermite_r2(l2,ih2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))

  ovrlp=0.0_gp
  do i2=1,n2
     qx=iw(3*nx+i2)
     qy=iw(4*nx+i2)
     qz=iw(5*nx+i2)
     fb=rw(n1+i2)
     do i1=1,n1
        px=iw(i1)
        py=iw(nx+i1)
        pz=iw(2*nx+i1)
        fa=rw(i1)

        fx=govrlp(a1,a2,dx,px,qx)
        fy=govrlp(a1,a2,dy,py,qy)
        fz=govrlp(a1,a2,dz,pz,qz)

        d2fx=kinovrlp(a1,a2,dx,px,qx)
        d2fy=kinovrlp(a1,a2,dy,py,qy)
        d2fz=kinovrlp(a1,a2,dz,pz,qz)

        ovrlp=ovrlp-0.5_gp*fa*fb*(d2fx*fy*fz+fx*d2fy*fz+fx*fy*d2fz)
        !print *,i1,i2,fx,fy,fz,fa,fb
     end do
  end do
 
END SUBROUTINE kinprod_h


subroutine calc_coeff_hermite_r2(l,ih,nterm_max,nterm,lx,ly,lz,fac_arr)
  use module_base
  implicit none
  integer, intent(in) :: l,ih,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(gp), dimension(nterm_max), intent(out) :: fac_arr
  !local variables
  integer, parameter :: norder_max=16
  integer :: ntermr,iterm
  integer, dimension(norder_max) :: lr
  integer, dimension(norder_max) :: nfac_arr

  if (l /=1) then
     stop 'error coeff_hermite'
  end if

  !take the coefficients for the herimt polynomials
  select case (ih)
  case (1)
     ntermr=1
     lr(1)=0
     nfac_arr(1)=1
  case (2)
     ntermr=1
     lr(1)=1
     nfac_arr(1)=1
  case (3)
     ntermr=2
     lr(1)=0
     lr(1)=2
     nfac_arr(1)=1
     nfac_arr(2)=-1
     !to be continued with all the coefficients
  case default
     stop 'order of hermite polynomial not provided'
  end select

  !and write them on the output arrays 
  nterm=1
  do iterm=1,ntermr
     if (lr(iterm) == 0) then
        lx(nterm)=0
        ly(nterm)=0
        lz(nterm)=0
        fac_arr(nterm)=real(nfac_arr(iterm),gp)
        nterm=nterm+1
     else
        lx(nterm)=2*lr(iterm) ;lx(nterm+1)=0          ;lx(nterm+2)=0
        ly(nterm)=0           ;ly(nterm+1)=2*lr(iterm);ly(nterm+2)=0          
        lz(nterm)=0           ;lz(nterm+1)=0          ;lz(nterm+2)=2*lr(iterm)          
        fac_arr(nterm)=real(nfac_arr(iterm),gp)
        fac_arr(nterm+1)=fac_arr(nterm)
        fac_arr(nterm+2)=fac_arr(nterm)
        nterm=nterm+3
     end if
  end do
  !last term
  nterm=nterm-1

END SUBROUTINE calc_coeff_hermite_r2


subroutine gaussians_to_wavelets_new_h(iproc,nproc,lr,orbs,hx,hy,hz,G,wfn_gau,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,orbs%nspinor,orbs%norbp), intent(in) :: wfn_gau
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
  !local variables
  integer :: iorb,ierr,ispinor,ncplx
  real(dp) :: normdev,tt,scpr,totnorm
  real(gp) :: kx,ky,kz

  if(iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
       'Writing wavefunctions in wavelet form...'

  normdev=0.0_dp
  tt=0.0_dp
  do iorb=1,orbs%norbp
     !features of the k-point ikpt
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))

     !evaluate the complexity of the k-point
     if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
        ncplx=1
     else
        ncplx=2
     end if
     totnorm=0.0_dp
     do ispinor=1,orbs%nspinor,ncplx
        !print *,'start',ispinor,ncplx,iorb,orbs%nspinor
        !the Block wavefunctions are exp(-Ikr) psi(r) (with MINUS k)
        call gaussians_to_wavelets_orb(ncplx,lr,hx,hy,hz,kx,ky,kz,G,&
             wfn_gau(1,ispinor,iorb),psi(1,ispinor,iorb))
        !print *,'end',ispinor,ncplx,iorb,orbs%nspinor
        call wnrm_wrap(ncplx,lr%wfd%nvctr_c,lr%wfd%nvctr_f,psi(1,ispinor,iorb),scpr) 
        totnorm=totnorm+scpr
     end do
     !write(*,'(1x,a,i5,1pe14.7,i3)')'norm of orbital ',iorb,totnorm,ncplx
     do ispinor=1,orbs%nspinor
        call wscal_wrap(lr%wfd%nvctr_c,lr%wfd%nvctr_f,real(1.0_dp/sqrt(totnorm),wp),&
             psi(1,ispinor,iorb))
     end do
     tt=max(tt,abs(1.0_dp-totnorm))
     !print *,'iorb,norm',totnorm
  end do

  if (iproc ==0  .and. verbose > 1) write(*,'(1x,a)')'done.'
  !renormalize the orbitals
  !calculate the deviation from 1 of the orbital norm
  if (nproc > 1) then
     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     normdev=tt
  end if
  if (iproc ==0) write(*,'(1x,a,1pe12.2)')&
       'Deviation from normalization of the imported orbitals',normdev

END SUBROUTINE gaussians_to_wavelets_new_h


subroutine gaussians_to_wavelets_orb_h(ncplx,lr,hx,hy,hz,kx,ky,kz,G,wfn_gau,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff), intent(in) :: wfn_gau
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='gaussians_to_wavelets_orb'
  integer, parameter :: nterm_max=48,maxsizeKB=2048,nw=65536
  integer,parameter:: ncplx_g=1 !true for NC pseudos
  logical :: perx,pery,perz
  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,i,nterm,ig
  integer :: nterms_max,nterms,iterm,n_gau,ml1,mu1,ml2,mu2,ml3,mu3 !n(c) iscoeff
  real(gp) :: rx,ry,rz,gau_a(ncplx_g)
  real(gp),parameter:: gau_cut=1.0_d0 !only useful for PAW
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(wp), allocatable, dimension(:,:,:) :: work
  real(wp), allocatable, dimension(:,:,:,:) :: wx,wy,wz

  !calculate nterms_max:
  !allows only maxsizeKB per one-dimensional array
  !(for a grid of dimension 100 nterms_max=655)
  !bu with at least ngx*nterm_max ~= 100 elements
  nterms_max=max(maxsizeKB*1024/(2*ncplx*max(lr%d%n1,lr%d%n2,lr%d%n3)),100)

  allocate(work(0:nw,2,2+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  allocate(wx(ncplx,0:lr%d%n1,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wx,'wx',subname)
  allocate(wy(ncplx,0:lr%d%n2,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wy,'wy',subname)
  allocate(wz(ncplx,0:lr%d%n3,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wz,'wz',subname)

  !conditions for periodicity in the three directions
  perx=(lr%geocode /= 'F')
  pery=(lr%geocode == 'P')
  perz=(lr%geocode /= 'F')

  !initialize the wavefunction
  call razero((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx,psi)

  !calculate the number of terms for this orbital
  nterms=0
  !loop over the atoms
  ishell=0
  iexpo=1
  icoeff=1
  !n(c) iscoeff=1
  iterm=1
  do iat=1,G%nat
     rx=G%rxyz(1,iat)
     ry=G%rxyz(2,iat)
     rz=G%rxyz(3,iat)
     !loop over the number of shells of the atom type
     do isat=1,G%nshell(iat)
        ishell=ishell+1
        !the degree of contraction of the basis function
        !is the same as the ng value of the createAtomicOrbitals routine
        ng=G%ndoc(ishell)
        !angular momentum of the basis set(shifted for compatibility with BigDFT routines
        l=G%nam(ishell)
        if (l/=1) then
           stop 'error l gaussians_to_wavelets_orb_h'
        end if
        !print *,iproc,iat,ishell,G%nam(ishell),G%nshell(iat)
        !multiply the values of the gaussian contraction times the orbital coefficient

        do m=1,2*l-1
           call calc_coeff_inguess(l,isat,nterm_max,nterm,lx,ly,lz,fac_arr)
           !control whether the basis element may be
           !contribute to some of the orbital of the processor
           if (wfn_gau(icoeff) /= 0.0_wp) then
              if (nterms + nterm*ng > nterms_max) then
                 !accumulate wavefuncton
                 call wfn_from_tensprod(lr,ncplx,nterms,wx,wy,wz,psi)
                 iterm=1
                 nterms=0
              end if
              !assign the arrays
              !make sure that the coefficients returned by 
              !gauss_to_daub are zero outside [ml:mr] 
              do ig=1,ng
                 do i=1,nterm
                    !print *,iat,ig,i,fac_arr(i),wfn_gau(icoeff),G%xp(1,iexpo+ig-1)
                    gau_a(:)=G%xp(1,iexpo+ig-1)
                    n_gau=lx(i)
                    !print *,'x',gau_a,nterm,ncplx,kx,ky,kz,ml1,mu1,lr%d%n1
                    call gauss_to_daub_k(hx,kx*hx,ncplx,ncplx_g,ncplx,fac_arr(i),rx,gau_a,n_gau,&
                         lr%ns1,lr%d%n1,ml1,mu1,&
                         wx(1,0,1,iterm),work,nw,perx,gau_cut) 
                    n_gau=ly(i)
                    !print *,'y',ml2,mu2,lr%d%n2
                    call gauss_to_daub_k(hy,ky*hy,ncplx,ncplx_g,ncplx,wfn_gau(icoeff),ry,gau_a,n_gau,&
                         lr%ns2,lr%d%n2,ml2,mu2,&
                         wy(1,0,1,iterm),work,nw,pery,gau_cut) 
                    n_gau=lz(i) 
                    !print *,'z',ml3,mu3,lr%d%n3
                    call gauss_to_daub_k(hz,kz*hz,ncplx,ncplx_g,ncplx,G%psiat(1,iexpo+ig-1),rz,gau_a,n_gau,&
                         lr%ns3,lr%d%n3,ml3,mu3,&
                         wz(1,0,1,iterm),work,nw,perz,gau_cut)
                    iterm=iterm+1
                 end do
              end do
              nterms=nterms+nterm*ng
           end if
           icoeff=icoeff+1
        end do
        iexpo=iexpo+ng
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  !accumulate wavefuncton
  call wfn_from_tensprod(lr,ncplx,nterms,wx,wy,wz,psi)
!psi=1.d0
  i_all=-product(shape(wx))*kind(wx)
  deallocate(wx,stat=i_stat)
  call memocc(i_stat,i_all,'wx',subname)
  i_all=-product(shape(wy))*kind(wy)
  deallocate(wy,stat=i_stat)
  call memocc(i_stat,i_all,'wy',subname)
  i_all=-product(shape(wz))*kind(wz)
  deallocate(wz,stat=i_stat)
  call memocc(i_stat,i_all,'wz',subname)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

END SUBROUTINE gaussians_to_wavelets_orb_h


!>   Overlap matrix between two different basis structures
!!   The first one is a gaussian hermite basis
subroutine gaussian_overlap_h(A,B,ovrlp)
  use module_base
  use module_types
  implicit none
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables
  integer, parameter :: niw=302,nrw=302
  integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
  integer :: iovrlp,jovrlp,jcoeff,jexpo
  integer :: ngA,ngB,lA,lB,mA,mB
  real(gp) :: dx,dy,dz
  integer, dimension(niw) :: iw
  real(gp), dimension(nrw) :: rw

  iovrlp=0
  ishell=0
  iexpo=1
  icoeff=1

  !loop on each shell (intensive calculation)
  do iat=1,A%nat
     do isat=1,A%nshell(iat)
        ishell=ishell+1
        ngA=A%ndoc(ishell)
        lA=A%nam(ishell)
        if (lA /= 1) then
           stop 'only s terms supported, gauovrlp'
        end if
        do mA=1,2*lA-1
           iovrlp=iovrlp+1

           jovrlp=0
           jshell=0
           jexpo=1
           jcoeff=1

           do jat=1,B%nat
              dx=B%rxyz(1,jat)-A%rxyz(1,iat)
              dy=B%rxyz(2,jat)-A%rxyz(2,iat)
              dz=B%rxyz(3,jat)-A%rxyz(3,iat)
              do jsat=1,B%nshell(jat)
                 jshell=jshell+1
                 ngB=B%ndoc(jshell)
                 lB=B%nam(jshell)
                 do mB=1,2*lB-1
                    jovrlp=jovrlp+1
                    if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                         A%ncoeff /= B%ncoeff ) then
                       call gbasovrlp_h(A%xp(1,iexpo),A%psiat(1,iexpo),&
                            B%xp(1,jexpo),B%psiat(1,jexpo),&
                            ngA,ngB,lA,isat,lB,mB,dx,dy,dz,&
                            niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                    end if
                 end do
                 jexpo=jexpo+ngB
                 jcoeff=jcoeff+2*lB-1
              end do
           end do
        end do
        iexpo=iexpo+ngA
        icoeff=icoeff+2*lA-1
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
  call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)
  
END SUBROUTINE gaussian_overlap_h


!>   Calculates the scalar product between two shells
!!   by considering only the nonzero coefficients
!!   actual building block for calculating overlap matrix
!!   inserted work arrays for calculation
!!   The first are Hermite polynomial basis
subroutine gbasovrlp_h(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,ih1,l2,m2,dx,dy,dz,&
     niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: ng1,ng2,l1,ih1,l2,m2,niw,nrw
  real(gp), intent(in) :: dx,dy,dz
  integer, dimension(niw) :: iw
  real(gp), dimension(nrw) :: rw
  real(gp), dimension(ng1), intent(in) :: expo1,coeff1
  real(gp), dimension(ng2), intent(in) :: expo2,coeff2
  real(gp), intent(out) :: ovrlp
  !local variables
  integer :: i1,i2
  real(gp) :: a1,a2,c1,c2,govrlpr

  ovrlp=0.d0
  do i1=1,ng1
     a1=expo1(i1)
     a1=0.5_gp/a1**2
     c1=coeff1(i1)
     do i2=1,ng2
        a2=expo2(i2)
        a2=0.5_gp/a2**2
        c2=coeff2(i2)
        call gprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,m2,niw,nrw,iw,rw,govrlpr)
        govrlpr=c1*govrlpr*c2
        !print *,c1,c2,govrlpr
        ovrlp=ovrlp+govrlpr
     end do
  end do
  
END SUBROUTINE gbasovrlp_h


!> calculates a dot product between two differents gaussians times spherical harmonics
!! the first one is supposed to be hermite polynomial matrix
subroutine gprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,m2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,ih1,m2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=48
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
  real(gp) :: fx,fy,fz,fa,fb,govrlp

  !calculates the number of different couples
  call calc_coeff_hermite_r2(l1,ih1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
  call calc_coeff_inguess(l2,m2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
  ovrlp=0.d0
  do i2=1,n2
     qx=iw(3*nx+i2)
     qy=iw(4*nx+i2)
     qz=iw(5*nx+i2)
     fb=rw(n1+i2)
     do i1=1,n1
        px=iw(i1)
        py=iw(nx+i1)
        pz=iw(2*nx+i1)
        fa=rw(i1)

        fx=govrlp(a1,a2,dx,px,qx)
        fy=govrlp(a1,a2,dy,py,qy)
        fz=govrlp(a1,a2,dz,pz,qz)

        ovrlp=ovrlp+fa*fb*fx*fy*fz
        !print *,i1,i2,fx,fy,fz,fa,fb
     end do
  end do
 
END SUBROUTINE gprod_h
