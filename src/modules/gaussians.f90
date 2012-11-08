!> @file
!!  Define the fortran types
!! @author
!!    Copyright (C) 2008-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Modules which contains the handling of operations in Gaussian basis
!! Spherical harmonics are used in the cartesian form
module gaussians

  use module_base, only:gp,memocc,ndebug

  private

  integer, parameter :: NSD_=2,EXPO_=1,COEFF_=2  !< positions of exponents and coefficients in the storage space
  integer, parameter :: NSHID_=3,DOC_=1,L_=2,N_=3  !< positions of shell identification numbers in shell id space
  integer, parameter :: NTERM_MAX_OVERLAP=62 !<maximum number of terms for the considered shells
  integer, parameter :: NTERM_MAX_KINETIC=190 !<maximum number of terms for the considered shells in the case of laplacian
  integer, parameter :: L_MAX=3 !< maximum number of angular momentum considered

  !> Structures of basis of gaussian functions
  type, public :: gaussian_basis
     integer :: nat  !< number of centers
     integer :: ncoeff !< number of total basis elements
     integer :: nshltot !< total number of shells (m quantum number ignored) 
     integer :: nexpo !< number of exponents (sum of the contractions)
     !storage units
     integer, dimension(:), pointer :: nshell !< number of shells for any of the centers
     integer, dimension(:), pointer :: ndoc,nam !< degree of contraction, angular momentum of any shell
     real(gp), dimension(:), pointer :: xp,psiat !<factors and values of the exponents (complex numbers are allowed)
     real(gp), dimension(:,:), pointer :: rxyz !<positions of the centers
  end type gaussian_basis

  !> Structures of basis of gaussian functions
  type :: gaussian_basis_new
     integer :: nat  !< number of centers
     integer :: ncoeff !< number of total basis elements
     integer :: nshltot !< total number of shells (m quantum number ignored) 
     integer :: nexpo !< number of exponents (sum of the contractions)
     integer :: ncplx !< =1 if traditional, =2 if complex gaussians
     !storage units
     integer, dimension(:), pointer :: nshell !< number of shells for any of the centers
     integer, dimension(:,:), pointer :: shid !< degree of contraction, angular momentum and principal quantum number
     real(gp), dimension(:,:), pointer :: sd !<sigma and contraction coefficients the exponents (complex numbers are allowed)
     real(gp), dimension(:,:), pointer :: rxyz !<positions of the centers
  end type gaussian_basis_new

  public :: gaudim_check,normalize_shell,gaussian_overlap,kinetic_overlap

contains

  function gaussian_basis_null() result(G)
    implicit none
    type(gaussian_basis_new) :: G
    G%nat=0
    G%ncoeff=0
    G%nshltot=0
    G%nexpo=0
    G%ncplx=1
    nullify(G%nshell)
    nullify(G%shid)
    nullify(G%sd)
    nullify(G%rxyz)
  end function gaussian_basis_null

  function gaussian_basis_init(nat,nshell,rxyz) result(G)
    implicit none
    integer, intent(in) :: nat
    integer, dimension(nat), intent(in) :: nshell
    real(gp), dimension(3,nat), intent(in), target :: rxyz
    type(gaussian_basis_new) :: G
    !local variables
    character(len=*), parameter :: subname='gaussian_basis_init'
    integer :: i_stat,iat

    G=gaussian_basis_null()

    G%nat=nat
    G%rxyz => rxyz

    !number of shells per atoms
    allocate(G%nshell(G%nat+ndebug),stat=i_stat)
    call memocc(i_stat,G%nat,'G%nshell',subname)

    G%nshltot=0
    do iat=1,nat
       G%nshell(iat)=nshell(iat)
       G%nshltot=G%nshltot+nshell(iat)
    end do

    allocate(G%shid(NSHID_,G%nshltot+ndebug),stat=i_stat)
    call memocc(i_stat,G%shid,'G%shid',subname)

  end function gaussian_basis_init

  subroutine gaussian_basis_convert(G,Gold)
    implicit none
    type(gaussian_basis_new), intent(out) :: G
    type(gaussian_basis), intent(in) :: Gold
    !local variables
    character(len=*), parameter :: subname='gaussian_basis_convert'
    integer :: ishell,i_stat,iexpo

    G=gaussian_basis_init(Gold%nat,Gold%nshell,Gold%rxyz)
    G%ncplx=1
    G%nexpo=0
    do ishell=1,G%nshltot
       G%shid(DOC_,ishell)=Gold%ndoc(ishell)
       G%shid(L_,ishell)=Gold%nam(ishell)-1 !<traditional convention for l restored
       G%shid(N_,ishell)=1 !<old structure only has N=1
       G%nexpo=G%nexpo+Gold%ndoc(ishell)
       G%ncoeff=G%ncoeff+2*Gold%nam(ishell)-1
    end do
    !allocate storage space (real exponents and coeffs for the moment)
    allocate(G%sd(G%ncplx*NSD_,G%nexpo+ndebug),stat=i_stat)
    call memocc(i_stat,G%sd,'G%sd',subname)

    do iexpo=1,G%nexpo
       G%sd(EXPO_,iexpo)=0.5_gp/Gold%xp(iexpo)**2
       G%sd(COEFF_,iexpo)=Gold%psiat(iexpo)
    end do

  end subroutine gaussian_basis_convert

  subroutine gaussian_basis_free(G,subname)
    implicit none
    character(len=*), intent(in) :: subname
    type(gaussian_basis_new), intent(inout) :: G
    !local variables
    integer :: i_all,i_stat

    !do not deallocate the atomic centers
    if (associated(G%rxyz)) nullify(G%rxyz)

    if (associated(G%sd)) then
       i_all=-product(shape(G%sd))*kind(G%sd)
       deallocate(G%sd,stat=i_stat)
       call memocc(i_stat,i_all,'G%sd',subname)
    end if
    if (associated(G%shid)) then
       i_all=-product(shape(G%shid))*kind(G%shid)
       deallocate(G%shid,stat=i_stat)
       call memocc(i_stat,i_all,'G%shid',subname)
    end if
    if (associated(G%nshell)) then
       i_all=-product(shape(G%nshell))*kind(G%nshell)
       deallocate(G%nshell,stat=i_stat)
       call memocc(i_stat,i_all,'G%nshell',subname)
    end if

    G=gaussian_basis_null()

  end subroutine gaussian_basis_free

  !>   Overlap matrix between two different basis structures
  subroutine gaussian_overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp
    !only lower triangular part for A%ncoeff=B%ncoeff
    !local variables
    integer, parameter :: niw=18,nrw=6
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo
    integer :: ngA,ngB,lA,lB,mA,mB
    real(gp) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    type(gaussian_basis_new) :: G,H

    call gaussian_basis_convert(G,A)
    call gaussian_basis_convert(H,B)

    call overlap(G,H,ovrlp)

    call gaussian_basis_free(G,'gaussian_overlap')
    call gaussian_basis_free(H,'gaussian_overlap')

    return

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
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
                      call gbasovrlp(A%xp(iexpo),A%psiat(iexpo),&
                           B%xp(jexpo),B%psiat(jexpo),&
                           ngA,ngB,lA,mA,lB,mB,dx,dy,dz,&
                           niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                      !end if
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

  END SUBROUTINE gaussian_overlap

  !>   Overlap kinetic matrix between two different basis structures
  !!   the kinetic operator is applicated on the A basis structure
  subroutine kinetic_overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
    !only lower triangular part for A%ncoeff=B%ncoeff
    !local variables
    integer, parameter :: niw=18,nrw=6
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo
    integer :: ngA,ngB,lA,lB,mA,mB
    real(gp) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    type(gaussian_basis_new) :: G,H

    call gaussian_basis_convert(G,A)
    call gaussian_basis_convert(H,B)

    call kinetic(G,H,ovrlp)

    call gaussian_basis_free(G,'gaussian_overlap')
    call gaussian_basis_free(H,'gaussian_overlap')

    return


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
                      if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff .or. &
                           A%ncoeff /= B%ncoeff ) then
                         call kineticovrlp(A%xp(iexpo),A%psiat(iexpo),&
                              B%xp(jexpo),B%psiat(jexpo),&
                              ngA,ngB,lA,mA,lB,mB,dx,dy,dz,&
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

  END SUBROUTINE kinetic_overlap

  !>   Calculates the scalar product between two shells
  !!   by considering only the nonzero coefficients
  !!   actual building block for calculating overlap matrix
  !!   inserted work arrays for calculation
  subroutine gbasovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
       niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
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
          call gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
          govrlpr=c1*govrlpr*c2
          !print *,c1,c2,govrlpr
          ovrlp=ovrlp+govrlpr
       end do
    end do

  END SUBROUTINE gbasovrlp

  !>   Calculates the scalar product between two shells
  !!   by considering only the nonzero coefficients
  !!   actual building block for calculating overlap matrix
  !!   inserted work arrays for calculation
  subroutine kineticovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
       niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
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
          call kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
          govrlpr=c1*govrlpr*c2
          !print *,c1,c2,govrlpr
          ovrlp=ovrlp+govrlpr
       end do
    end do

  END SUBROUTINE kineticovrlp

  !>   Calculates a dot product between two differents gaussians times spherical harmonics
  !!   valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
  !!   to be rearranged when only some of them is zero
  subroutine gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
    real(gp), intent(in) :: a1,a2,dx,dy,dz
    integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
    real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
    real(gp), intent(out) :: ovrlp
    !local variables
    integer, parameter :: nx=3
    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
    real(gp) :: fx,fy,fz,fa,fb!,govrlp

    !calculates the number of different couples
    call calc_coeff_inguess(l1,m1,nx,n1,&
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

  END SUBROUTINE gprod


  !>   Overlap matrix between two different basis structures
  subroutine overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis_new), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp 
    !local variables
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo,itpdA,itpdB
    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB
    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
    integer, dimension(3,NTERM_MAX_OVERLAP) :: powA,powB
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdA,ftpdB
    real(gp), dimension(3) :: dr,rA

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       rA(1)=A%rxyz(1,iat)
       rA(2)=A%rxyz(2,iat)
       rA(3)=A%rxyz(3,iat)
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%shid(DOC_,ishell)
          lA=A%shid(L_,ishell)
          nA=A%shid(N_,ishell)
          call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
          itpdA=1
          do mA=1,2*lA+1
             !here the entire array should be extracted
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                !here boundary conditions should be considered
                dr(1)=B%rxyz(1,jat)-rA(1)
                dr(2)=B%rxyz(2,jat)-rA(2)
                dr(3)=B%rxyz(3,jat)-rA(3)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%shid(DOC_,jshell)
                   lB=B%shid(L_,jshell)
                   nB=B%shid(N_,jshell)
                   call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)
                   itpdB=1
                   do mB=1,2*lB+1
                      jovrlp=jovrlp+1
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
                      ovrlp(iovrlp,jovrlp)=&
                           gdot(ngA,A%sd(1,iexpo),ntpdA(mA),powA(1,itpdA),ftpdA(itpdA),&
                           ngB,B%sd(1,jexpo),ntpdB(mB),powB(1,itpdB),ftpdB(itpdB),dr)
                      !end if
                      itpdB=itpdB+ntpdB(mB)
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB+1
                end do
             end do
             itpdA=itpdA+ntpdA(mA)
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA+1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE overlap


  !> Calculates a dot product between two basis functions
  !! Basis function is identified by its tensor product decompositions and sigmas
  !! the contraction coefficients are also given
  pure function gdot(ng1,sd1,ntpd1,pws1,ftpd1,ng2,sd2,ntpd2,pws2,ftpd2,dr) result(ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2 !< degrees of contractions of the basis
    integer, intent(in) :: ntpd1,ntpd2 !<number of tensor product decopositions
    integer, dimension(3,ntpd1), intent(in) :: pws1 !<power coefficients for each term and direction
    integer, dimension(3,ntpd2), intent(in) :: pws2 !<power coefficients for each term and directio
    real(gp), dimension(ntpd1), intent(in) :: ftpd1 !<factors of tensor product decompositions
    real(gp), dimension(ntpd2), intent(in) :: ftpd2 !<factors of tensor product decompositions
    real(gp), dimension(3), intent(in) :: dr !<separations between basis functions
    real(gp), dimension(NSD_,ng1), intent(in) :: sd1 !<exponents and coefficient
    real(gp), dimension(NSD_,ng2), intent(in) :: sd2 !<exponents and coefficient
    real(gp):: ovrlp !<scalar product
    !local variables
    integer :: ig1,ig2,i2,i1
    real(gp) ::  f,s1,s2,d1,d2,integral

    ovrlp=0.0_gp
    do ig2=1,ng2
       s2=sd2(EXPO_,ig2)
       d2=sd2(COEFF_,ig2)
       do ig1=1,ng1
          s1=sd1(EXPO_,ig1)
          d1=sd1(COEFF_,ig1)
          integral=0.0_gp
          do i2=1,ntpd2
             do i1=1,ntpd1
                f=  govrlp(s1,s2,dr(1),pws1(1,i1),pws2(1,i2))
                f=f*govrlp(s1,s2,dr(2),pws1(2,i1),pws2(2,i2))
                f=f*govrlp(s1,s2,dr(3),pws1(3,i1),pws2(3,i2))
                integral=integral+ftpd1(i1)*ftpd2(i2)*f
             end do
          end do
          ovrlp=ovrlp+d1*d2*integral
       end do
    end do
  end function gdot

  !perform the gaussian product for all the terms in the shell
  pure subroutine gdot_shell(sd1,l1,ntpdsh1,ntpd1,pws1,ftpd1,&
       sd2,l2,ntpdsh2,ntpd2,pws2,ftpd2,dr,overlap)
    implicit none
    integer, intent(in) :: l1,l2 !<angular momenta of the shell
    integer, intent(in) :: ntpdsh1,ntpdsh2 !<total number of shell tpd
    real(gp), dimension(NSD_), intent(in) :: sd1,sd2 !<exponents and coefficient
    integer, dimension(2*l1+1), intent(in) :: ntpd1 !<number of terms
    integer, dimension(2*l2+1), intent(in) :: ntpd2 !<number of terms
    integer, dimension(3,ntpdsh1), intent(in) :: pws1 !<exponents
    integer, dimension(3,ntpdsh2), intent(in) :: pws2 !<exponents
    real(gp), dimension(ntpdsh1), intent(in) :: ftpd1
    real(gp), dimension(ntpdsh2), intent(in) :: ftpd2
    real(gp), dimension(3), intent(in) :: dr !<separations between basis functions
    real(gp), dimension(2*l1+1,2*l2+1), intent(inout) :: overlap !<overlap of the shell
    !local variables
    integer :: m1,m2,i1,i2,itpd1,itpd2
    real(gp) :: f,integral

    itpd2=0
    do m2=1,2*l2+1
       itpd1=0
       do m1=1,2*l1+1
          integral=0.0_gp
          do i2=1,ntpd2(m2)
             do i1=1,ntpd1(m1)
                f=  govrlp(sd1(EXPO_),sd2(EXPO_),&
                     dr(1),pws1(1,i1+itpd1),pws2(1,i2+itpd2))
                f=f*govrlp(sd1(EXPO_),sd2(EXPO_),&
                     dr(2),pws1(2,i1+itpd1),pws2(2,i2+itpd2))
                f=f*govrlp(sd1(EXPO_),sd2(EXPO_),&
                     dr(3),pws1(3,i1+itpd1),pws2(3,i2+itpd2))
                integral=integral+ftpd1(i1+itpd1)*ftpd2(i2+itpd2)*f
             end do
          end do
          overlap(m1,m2)=overlap(m1,m2)+sd1(COEFF_)*sd2(COEFF_)*integral
          itpd1=itpd1+ntpd1(m1)
       end do
       itpd2=itpd2+ntpd2(m2)
    end do

  end subroutine gdot_shell

  !>   Overlap matrix between two different basis structures
  !!   laplacian is applied to the first one
  subroutine kinetic(A,B,ovrlp)
    implicit none
    type(gaussian_basis_new), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp
    !local variables
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo,igA,igB,ioverlap
    integer :: ngA,ngB,lA,lB,mA,mB,nA,nB,ntpdshA,ntpdshB
    integer, dimension(2*L_MAX+1) :: ntpdA,ntpdB
    integer, dimension(3,NTERM_MAX_KINETIC) :: powA
    real(gp), dimension(NTERM_MAX_KINETIC) :: ftpdA
    integer, dimension(3,NTERM_MAX_OVERLAP) :: powB
    real(gp), dimension(NTERM_MAX_OVERLAP) :: ftpdB
    real(gp), dimension((2*L_MAX+1)*(2*L_MAX+1)) :: shell_overlap
    real(gp), dimension(3) :: dr,rA

    iovrlp=0
    ishell=0
    iexpo=0
    icoeff=1
    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       rA(1)=A%rxyz(1,iat)
       rA(2)=A%rxyz(2,iat)
       rA(3)=A%rxyz(3,iat)
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%shid(DOC_,ishell)
          lA=A%shid(L_,ishell)
          nA=A%shid(N_,ishell)

          jovrlp=0
          jshell=0
          jexpo=0
          jcoeff=1
          do jat=1,B%nat
             !here boundary conditions should be considered
             dr(1)=B%rxyz(1,jat)-rA(1)
             dr(2)=B%rxyz(2,jat)-rA(2)
             dr(3)=B%rxyz(3,jat)-rA(3)
             do jsat=1,B%nshell(jat)
                jshell=jshell+1
                ngB=B%shid(DOC_,jshell)
                lB=B%shid(L_,jshell)
                nB=B%shid(N_,jshell)
                !calculation of shell decomposition is independent of exponent
                call tensor_product_decomposition(nB,lB,ntpdshB,ntpdB,powB,ftpdB)

                shell_overlap=0.0_gp !overlap of the different terms
                do igA=1,ngA
                   call tensor_product_decomposition_laplacian(A%sd(EXPO_,igA+iexpo),nA,lA,&
                        ntpdshA,ntpdA,powA,ftpdA)
!                   call tensor_product_decomposition(nA,lA,ntpdshA,ntpdA,powA,ftpdA)
                   do igB=1,ngB
                      call gdot_shell(A%sd(1,igA+iexpo),lA,ntpdshA,ntpdA,powA,ftpdA,&
                           B%sd(1,igB+jexpo),lB,ntpdshB,ntpdB,powB,ftpdB,dr,&
                           shell_overlap)
                   end do
                end do

                !here the entire array should be copied in the right place
                do mB=1,2*lB+1
                   do mA=1,2*lA+1
                      ovrlp(iovrlp+mA,jovrlp+mB)=-0.5_gp*&
                           shell_overlap(mA+(mB-1)*(2*lA+1))
                   end do
                end do
                jexpo=jexpo+ngB
                jcoeff=jcoeff+2*lB+1
                jovrlp=jovrlp+2*lB+1
             end do
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA+1
          iovrlp=iovrlp+2*lA+1
       end do
    end do

    call gaudim_check(iexpo+1,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo+1,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE kinetic

  !>   Calculates @f$\int \exp^{-a1*x^2} x^l1 \exp^{-a2*(x-d)^2} (x-d)^l2 dx@f$
  !!   Uses gauint0 if d==0
  !!
  pure function govrlp(a1,a2,d,l1,l2)
    implicit none
    integer, intent(in) :: l1,l2
    real(gp), intent(in) :: a1,a2,d
    real(gp) :: govrlp
    !local variables
    integer :: p
    real(gp) :: prefac,stot,aeff,ceff,tt,fsum!,gauint,gauint0

    !quick check
    if (d==0.0_gp) then
       govrlp=gauint0(a1+a2,l1+l2)
       return
    end if

    !build the prefactor
    prefac=a1+a2
    prefac=a2/prefac
    prefac=a1*prefac
    prefac=-d**2*prefac
    prefac=dexp(prefac)

    !build the effective exponent and coefficients
    aeff=a1+a2
    ceff=a2*d
    ceff=ceff/aeff

    !build the first term in the sum
    stot=(-d)**l2
    stot=gauint(aeff,ceff,l1)*stot

    !perform the sum
    do p=1,l2/2
       tt=rfac(1,p)
       fsum=rfac(l2-p+1,l2)
       fsum=fsum/tt
       tt=(-d)**(l2-p)
       fsum=fsum*tt
       tt=gauint(aeff,ceff,l1+p)
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l2/2+1,l2
       tt=rfac(1,l2-p)
       fsum=rfac(p+1,l2)
       fsum=fsum/tt
       tt=(-d)**(l2-p)
       fsum=fsum*tt
       tt=gauint(aeff,ceff,l1+p)
       fsum=fsum*tt
       stot=stot+fsum
    end do

    !final result
    govrlp=prefac*stot
  END FUNCTION govrlp
  
  !>   Kinetic overlap between gaussians, based on cartesian coordinates
  !!   calculates a dot product between two differents gaussians times spherical harmonics
  !!   valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
  !!   to be rearranged when only some of them is zero
  !!
  subroutine kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
    real(gp), intent(in) :: a1,a2,dx,dy,dz
    integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
    real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
    real(gp), intent(out) :: ovrlp
    !local variables
    integer, parameter :: nx=3
    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
    real(gp) :: fx,fy,fz,fa,fb,d2fx,d2fy,d2fz!,govrlp,kinovrlp

    !calculates the number of different couples
    call calc_coeff_inguess(l1,m1,nx,n1,&
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

          d2fx=kinovrlp(a1,a2,dx,px,qx)
          d2fy=kinovrlp(a1,a2,dy,py,qy)
          d2fz=kinovrlp(a1,a2,dz,pz,qz)

          ovrlp=ovrlp-0.5_gp*fa*fb*(d2fx*fy*fz+fx*d2fy*fz+fx*fy*d2fz)
          !print *,i1,i2,fx,fy,fz,fa,fb
       end do
    end do

  END SUBROUTINE kinprod

  !>   Calculates @f$\int d^2/dx^2(\exp^{-a1*x^2} x^l1) \exp^{-a2*(x-d)^2} (x-d)^l2 dx@f$
  !!   in terms of the govrlp function below
  pure function kinovrlp(a1,a2,d,l1,l2)
    implicit none
    integer, intent(in) :: l1,l2
    real(gp), intent(in) :: a1,a2,d
    real(gp) :: kinovrlp
    !local variables
    real(gp) :: fac,ovrlp!govrlp

    !case l1+2
    fac=4._gp*a1**2
    ovrlp=govrlp(a1,a2,d,l1+2,l2)
    kinovrlp=fac*ovrlp
    !case l1
    fac=2._gp*a1*real(2*l1+1,gp)
    ovrlp=govrlp(a1,a2,d,l1,l2)
    kinovrlp=kinovrlp-fac*ovrlp
    !case l1-2 (if applicable)
    if (l1 >=2) then
       fac=real(l1*(l1-1),gp)
       ovrlp=govrlp(a1,a2,d,l1-2,l2)
       kinovrlp=kinovrlp+fac*ovrlp
    end if
  END FUNCTION kinovrlp

  !>   Calculates @f$\int \exp^{-a*x^2} x^l dx@f$
  !!   this works for all l
  pure function gauint0(a,l)
    implicit none
    integer, intent(in) :: l
    real(gp), intent(in) :: a
    real(gp) :: gauint0
    !local variables
    real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
    integer :: p
    real(gp) :: prefac,tt
    !build the prefactor
    prefac=sqrt(a)
    prefac=1.d0/prefac
    prefac=gammaonehalf*prefac**(l+1)

    p=l/2
    if (2*p < l) then
       gauint0=0.0_gp
       return
    end if
    tt=xfac(1,p,-0.5d0)
    !final result
    gauint0=prefac*tt

  END FUNCTION gauint0



  !>   Calculates @f$\int \exp^{-a*(x-c)^2} x^l dx@f$
  !!   this works ONLY when c/=0.d0
  !!
  !!
  pure function gauint(a,c,l)
    implicit none
    integer, intent(in) :: l
    real(gp), intent(in) :: a,c
    real(gp) :: gauint
    !local variables
    real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
    integer :: p
    real(gp) :: prefac,stot,fsum,tt!,firstprod

    !quick check
    !if (c==0.0_gp) then
    !   stop 'gauint0 should be called'
    !end if

    !build the prefactor
    prefac=sqrt(a)
    prefac=1.d0/prefac
    prefac=gammaonehalf*prefac

    !the first term of the sum is one
    !but we have to multiply for the prefactor
    stot=c**l

    !if (c==0.0_gp .and. l==0) then
    !   do p=0,20
    !      print *,'stot,p',stot,a,p,gauint0(a,p)
    !   end do
    !end if

    !calculate the sum
    do p=1,l/4
       tt=rfac(p+1,2*p)
       fsum=rfac(l-2*p+1,l)
       fsum=fsum/tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l/4+1,l/3
       tt=rfac(p+1,l-2*p)
       fsum=rfac(2*p+1,l)
       fsum=fsum/tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l/3+1,l/2
       tt=rfac(l-2*p+1,p)
       fsum=rfac(2*p+1,l)
       fsum=fsum*tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do

    !final result
    gauint=stot*prefac

  END FUNCTION gauint

  !>
  pure function firstprod(p)
    implicit none
    integer, intent(in) :: p
    real(gp) :: firstprod
    !local variables
    integer :: i
    real(gp) :: tt
    firstprod=1.0_gp
    do i=1,p
       tt=real(2*i,gp)
       tt=1.0_gp/tt
       tt=1.0_gp-tt
       firstprod=firstprod*tt
    end do
  END FUNCTION firstprod

  !>
  subroutine gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)
    implicit none
    integer, intent(in) :: iexpo,icoeff,ishell,nexpo,ncoeff,nshltot
    !check of the dimensions
    if (iexpo /= nexpo+1) then
       write(*,*)' ERROR: nexpo+1 <> iexpo',nexpo,iexpo
       stop
    else if (icoeff /= ncoeff+1) then
       write(*,*)' ERROR: ncoeff+1 <> icoeff',ncoeff,icoeff
       stop
    else if (ishell /= nshltot) then
       write(*,*)' ERROR: nshltot <> ishell',nshltot,ishell
       stop
    end if
  END SUBROUTINE gaudim_check

  !>   Normalize a given atomic shell following the angular momentum
  !!
  !!
  pure subroutine normalize_shell(ng,l,expo,coeff)
    implicit none
    integer, intent(in) :: ng,l
    real(gp), dimension(ng), intent(in) :: expo
    real(gp), dimension(ng), intent(inout) :: coeff
    !local variables
    integer :: i,j
    real(gp) :: norm,tt,e1,ex,c1,c2!,gauint0

    norm=0.d0
    do i=1,ng
       e1=expo(i)
       c1=coeff(i)
       do j=1,ng
          ex=expo(j)+e1
          c2=coeff(j)
          tt=gauint0(ex,2*l+2)
          norm=norm+c1*tt*c2
       end do
    end do
    norm=sqrt(0.5_gp*norm)
    norm=1.0_gp/norm
    do i=1,ng
       coeff(i)=coeff(i)*norm
    end do

    !print *,'l=',l,'norm=',norm

  END SUBROUTINE normalize_shell

  !> Factorial (float)
  pure function rfac(is,ie)
    implicit none
    integer, intent(in) :: is,ie
    real(gp) :: rfac
    !local variables
    integer :: i
    real(gp) :: tt
    rfac=1.d0
    do i=is,ie
       tt=real(i,gp)
       rfac=rfac*tt
    end do
  END FUNCTION rfac

  !> Partial factorial, with real shift
  !!With this function n!=xfac(1,n,0.d0)
  pure function xfac(is,ie,sh)
    implicit none
    integer, intent(in) :: is,ie
    real(gp), intent(in) :: sh
    real(gp) :: xfac
    !local variables
    integer :: i
    real(gp) :: tt
    xfac=1.d0
    do i=is,ie
       tt=real(i,gp)+sh
       xfac=xfac*tt
    end do
  END FUNCTION xfac

  !> Routine to extract the coefficients from the quantum numbers and the operation
  pure subroutine tensor_product_decomposition(n,l,ntpd_shell,ntpd,pow,ftpd)
    implicit none
    integer, intent(in) :: n,l
    integer, intent(out) :: ntpd_shell !< No. of terms for the whole shell
    integer, dimension(2*l+1), intent(out) :: ntpd !< number of terms per shell element
    integer, dimension(3,NTERM_MAX_OVERLAP), intent(out) :: pow !< tensor product decompositions
    real(gp), dimension(NTERM_MAX_OVERLAP), intent(out) :: ftpd !<factors
    !All data alltogether
    real(gp), parameter :: rsp=0.564189583547756286948079451560772585844050629328998856844086_gp !1/sqrt(pi)
    select case(n)
    case(1)
       select case(l)
       case(0)
          ntpd_shell=1 ! l=0,n=1
          ntpd(1:1)=(/1/) ! l=0,n=1
          pow(1:3,1:1)=reshape((/0,0,0/),(/3,1/)) ! l=0,n=1
          ftpd(1:1)=rsp*(/0.5_gp/) ! m=1, l=0, n=1
       case(1)
          ntpd_shell=3 ! l=1,n=1
          ntpd(1:3)=(/1,1,1/) ! l=1,n=1
          pow(1:3,1:3)=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/)) ! l=1,n=1
          ftpd(1:1)=rsp*(/sqrt(0.75_gp)/) ! m=1, l=1, n=1
          ftpd(2:2)=rsp*(/sqrt(0.75_gp)/) ! m=2, l=1, n=1
          ftpd(3:3)=rsp*(/sqrt(0.75_gp)/) ! m=3, l=1, n=1
       case(2)
          ntpd_shell=8 ! l=2,n=1
          ntpd(1:5)=(/1,1,1,2,3/) ! l=2,n=1
          pow(1:3,1:8)=reshape((/0,1,1,1,0,1,1,1,0,2,0,0,0,2,0,2,0,0,0,2,0,0,0,2/),(/3,8/)) ! l=2,n=1
          ftpd(1:1)=rsp*(/sqrt(3.75_gp)/) ! m=1, l=2, n=1
          ftpd(2:2)=rsp*(/sqrt(3.75_gp)/) ! m=2, l=2, n=1
          ftpd(3:3)=rsp*(/sqrt(3.75_gp)/) ! m=3, l=2, n=1
          ftpd(4:5)=rsp*(/sqrt(0.9375_gp),-sqrt(0.9375_gp)/) ! m=4, l=2, n=1
          ftpd(6:8)=rsp*(/-sqrt(0.3125_gp),-sqrt(0.3125_gp),sqrt(1.25_gp)/) ! m=5, l=2, n=1
       case(3)
          ntpd_shell=16 ! l=3,n=1
          ntpd(1:7)=(/3,3,3,2,2,2,1/) ! l=3,n=1
          pow(1:3,1:16)=reshape((/3,0,0,1,2,0,1,0,2,2,1,0,0,3,0,0,1,2,2,0,1,0,2,1,0,0,3,3,0,0,1,&
               2,0,2,1,0,0,3,0,2,0,1,0,2,1,1,1,1/),(/3,16/)) ! l=3,n=1
          ftpd(1:3)=rsp*(/sqrt(0.65625_gp),sqrt(0.65625_gp),-sqrt(10.5_gp)/) ! m=1, l=3, n=1
          ftpd(4:6)=rsp*(/sqrt(0.65625_gp),sqrt(0.65625_gp),-sqrt(10.5_gp)/) ! m=2, l=3, n=1
          ftpd(7:9)=rsp*(/sqrt(3.9375_gp),sqrt(3.9375_gp),-sqrt(1.75_gp)/) ! m=3, l=3, n=1
          ftpd(10:11)=rsp*(/sqrt(1.09375_gp),-sqrt(9.84375_gp)/) ! m=4, l=3, n=1
          ftpd(12:13)=rsp*(/-sqrt(9.84375_gp),sqrt(1.09375_gp)/) ! m=5, l=3, n=1
          ftpd(14:15)=rsp*(/sqrt(6.5625_gp),-sqrt(6.5625_gp)/) ! m=6, l=3, n=1
          ftpd(16:16)=rsp*(/sqrt(26.25_gp)/) ! m=7, l=3, n=1
       end select
    case(2)
       select case(l)
       case(0)
          ntpd_shell=3 ! l=0,n=2
          ntpd(1:1)=(/3/) ! l=0,n=2
          pow(1:3,1:3)=reshape((/2,0,0,0,2,0,0,0,2/),(/3,3/)) ! l=0,n=2
          ftpd(1:3)=rsp*(/0.5_gp,0.5_gp,0.5_gp/) ! m=1, l=0, n=2
       case(1)
          ntpd_shell=9 ! l=1,n=2
          ntpd(1:3)=(/3,3,3/) ! l=1,n=2
          pow(1:3,1:9)=reshape((/3,0,0,1,2,0,1,0,2,2,1,0,0,3,0,0,1,2,2,0,1,0,2,1,0,0,3/),(/3,9/)) ! l=1,n=2
          ftpd(1:3)=rsp*(/sqrt(0.75_gp),sqrt(0.75_gp),sqrt(0.75_gp)/) ! m=1, l=1, n=2
          ftpd(4:6)=rsp*(/sqrt(0.75_gp),sqrt(0.75_gp),sqrt(0.75_gp)/) ! m=2, l=1, n=2
          ftpd(7:9)=rsp*(/sqrt(0.75_gp),sqrt(0.75_gp),sqrt(0.75_gp)/) ! m=3, l=1, n=2
       case(2)
          ntpd_shell=19 ! l=2,n=2
          ntpd(1:5)=(/3,3,3,4,6/) ! l=2,n=2
          pow(1:3,1:19)=reshape((/2,1,1,0,3,1,0,1,3,3,0,1,1,2,1,1,0,3,3,1,0,1,3,0,1,1,2,4,0,0,0,4,&
               0,2,0,2,0,2,2,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,19/)) ! l=2,n=2
          ftpd(1:3)=rsp*(/sqrt(3.75_gp),sqrt(3.75_gp),sqrt(3.75_gp)/) ! m=1, l=2, n=2
          ftpd(4:6)=rsp*(/sqrt(3.75_gp),sqrt(3.75_gp),sqrt(3.75_gp)/) ! m=2, l=2, n=2
          ftpd(7:9)=rsp*(/sqrt(3.75_gp),sqrt(3.75_gp),sqrt(3.75_gp)/) ! m=3, l=2, n=2
          ftpd(10:13)=rsp*(/sqrt(0.9375_gp),-sqrt(0.9375_gp),sqrt(0.9375_gp),-sqrt(0.9375_gp)/) ! m=4, l=2, n=2
          ftpd(14:19)=rsp*(/-sqrt(0.3125_gp),-sqrt(1.25_gp),-sqrt(0.3125_gp),sqrt(0.3125_gp),sqrt(0.3125_gp),&
               sqrt(1.25_gp)/) ! m=5, l=2, n=2
       case(3)
          ntpd_shell=35 ! l=3,n=2
          ntpd(1:7)=(/6,6,6,5,5,4,3/) ! l=3,n=2
          pow(1:3,1:35)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,4,&
               0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,&
               4,0,1,0,4,1,2,0,3,0,2,3,3,1,1,1,3,1,1,1,3/),(/3,35/)) ! l=3,n=2
          ftpd(1:6)=rsp*(/sqrt(0.65625_gp),sqrt(2.625_gp),sqrt(0.65625_gp),-sqrt(5.90625_gp),-sqrt(5.90625_gp),&
               -sqrt(10.5_gp)/) ! m=1, l=3, n=2
          ftpd(7:12)=rsp*(/sqrt(0.65625_gp),sqrt(2.625_gp),sqrt(0.65625_gp),-sqrt(5.90625_gp),-sqrt(5.90625_gp),&
               -sqrt(10.5_gp)/) ! m=2, l=3, n=2
          ftpd(13:18)=rsp*(/sqrt(3.9375_gp),sqrt(15.75_gp),sqrt(3.9375_gp),sqrt(0.4375_gp),sqrt(0.4375_gp),&
               -sqrt(1.75_gp)/) ! m=3, l=3, n=2
          ftpd(19:23)=rsp*(/sqrt(1.09375_gp),-sqrt(4.375_gp),-sqrt(9.84375_gp),sqrt(1.09375_gp),-sqrt(9.84375_gp)/) ! m=4, l=3, n=2
          ftpd(24:28)=rsp*(/-sqrt(9.84375_gp),-sqrt(4.375_gp),sqrt(1.09375_gp),-sqrt(9.84375_gp),sqrt(1.09375_gp)/) ! m=5, l=3, n=2
          ftpd(29:32)=rsp*(/sqrt(6.5625_gp),-sqrt(6.5625_gp),sqrt(6.5625_gp),-sqrt(6.5625_gp)/) ! m=6, l=3, n=2
          ftpd(33:35)=rsp*(/sqrt(26.25_gp),sqrt(26.25_gp),sqrt(26.25_gp)/) ! m=7, l=3, n=2
       end select
    case(3)
       select case(l)
       case(0)
          ntpd_shell=6 ! l=0,n=3
          ntpd(1:1)=(/6/) ! l=0,n=3
          pow(1:3,1:6)=reshape((/4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,6/)) ! l=0,n=3
          ftpd(1:6)=rsp*(/0.5_gp,1._gp,0.5_gp,1._gp,1._gp,0.5_gp/) ! m=1, l=0, n=3
       case(1)
          ntpd_shell=18 ! l=1,n=3
          ntpd(1:3)=(/6,6,6/) ! l=1,n=3
          pow(1:3,1:18)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,4,&
               0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5/),(/3,18/)) ! l=1,n=3
          ftpd(1:6)=rsp*(/sqrt(0.75_gp),sqrt(3._gp),sqrt(0.75_gp),sqrt(3._gp),sqrt(3._gp),sqrt(0.75_gp)/) ! m=1, l=1, n=3
          ftpd(7:12)=rsp*(/sqrt(0.75_gp),sqrt(3._gp),sqrt(0.75_gp),sqrt(3._gp),sqrt(3._gp),sqrt(0.75_gp)/) ! m=2, l=1, n=3
          ftpd(13:18)=rsp*(/sqrt(0.75_gp),sqrt(3._gp),sqrt(0.75_gp),sqrt(3._gp),sqrt(3._gp),sqrt(0.75_gp)/) ! m=3, l=1, n=3
       case(2)
          ntpd_shell=33 ! l=2,n=3
          ntpd(1:5)=(/6,6,6,8,7/) ! l=2,n=3
          pow(1:3,1:33)=reshape((/4,1,1,2,3,1,0,5,1,2,1,3,0,3,3,0,1,5,5,0,1,3,2,1,1,4,1,3,0,3,1,2,3,1,0,5,5,1,&
               0,3,3,0,1,5,0,3,1,2,1,3,2,1,1,4,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,0,4,2,2,0,4,0,2,4,6,0,0,4,2,0,2,4,&
               0,0,6,0,2,0,4,0,2,4,0,0,6/),(/3,33/)) ! l=2,n=3
          ftpd(1:6)=rsp*(/sqrt(3.75_gp),sqrt(15._gp),sqrt(3.75_gp),sqrt(15._gp),sqrt(15._gp),sqrt(3.75_gp)/) ! m=1, l=2, n=3
          ftpd(7:12)=rsp*(/sqrt(3.75_gp),sqrt(15._gp),sqrt(3.75_gp),sqrt(15._gp),sqrt(15._gp),sqrt(3.75_gp)/) ! m=2, l=2, n=3
          ftpd(13:18)=rsp*(/sqrt(3.75_gp),sqrt(15._gp),sqrt(3.75_gp),sqrt(15._gp),sqrt(15._gp),sqrt(3.75_gp)/) ! m=3, l=2, n=3
          ftpd(19:26)=rsp*(/sqrt(0.9375_gp),sqrt(0.9375_gp),-sqrt(0.9375_gp),-sqrt(0.9375_gp),sqrt(3.75_gp),&
               -sqrt(3.75_gp),sqrt(0.9375_gp),-sqrt(0.9375_gp)/) ! m=4, l=2, n=3
          ftpd(27:33)=rsp*(/-sqrt(0.3125_gp),-sqrt(2.8125_gp),-sqrt(2.8125_gp),-sqrt(0.3125_gp),sqrt(2.8125_gp),&
               sqrt(2.8125_gp),sqrt(1.25_gp)/) ! m=5, l=2, n=3
       case(3)
          ntpd_shell=62 ! l=3,n=3
          ntpd(1:7)=(/10,10,10,9,9,8,6/) ! l=3,n=3
          pow(1:3,1:62)=reshape((/7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,1,0,6,6,1,0,4,3,0,2,5,0,0,&
               7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,0,1,6,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,2,2,3,0,4,3,2,0,5,0,2,5,0,0,7,7,&
               0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,&
               3,4,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,0,4,3,2,0,5,0,2,5,5,1,1,3,3,1,1,5,1,3,1,3,1,3,3,1,1,5/),(/3,62/)) ! l=3,n=3
          ftpd(1:10)=rsp*(/sqrt(0.65625_gp),sqrt(5.90625_gp),sqrt(5.90625_gp),sqrt(0.65625_gp),-sqrt(2.625_gp),&
               -sqrt(10.5_gp),-sqrt(2.625_gp),-sqrt(32.15625_gp),-sqrt(32.15625_gp),-sqrt(10.5_gp)/) ! m=1, l=3, n=3
          ftpd(11:20)=rsp*(/sqrt(0.65625_gp),sqrt(5.90625_gp),sqrt(5.90625_gp),sqrt(0.65625_gp),-sqrt(2.625_gp),&
               -sqrt(10.5_gp),-sqrt(2.625_gp),-sqrt(32.15625_gp),-sqrt(32.15625_gp),-sqrt(10.5_gp)/) ! m=2, l=3, n=3
          ftpd(21:30)=rsp*(/sqrt(3.9375_gp),sqrt(35.4375_gp),sqrt(35.4375_gp),sqrt(3.9375_gp),sqrt(7._gp),&
               sqrt(28._gp),sqrt(7._gp),-sqrt(0.4375_gp),-sqrt(0.4375_gp),-sqrt(1.75_gp)/) ! m=3, l=3, n=3
          ftpd(31:39)=rsp*(/sqrt(1.09375_gp),-sqrt(1.09375_gp),-sqrt(27.34375_gp),-sqrt(9.84375_gp),sqrt(4.375_gp),&
               -sqrt(17.5_gp),-sqrt(39.375_gp),sqrt(1.09375_gp),-sqrt(9.84375_gp)/) ! m=4, l=3, n=3
          ftpd(40:48)=rsp*(/-sqrt(9.84375_gp),-sqrt(27.34375_gp),-sqrt(1.09375_gp),sqrt(1.09375_gp),-sqrt(39.375_gp),&
               -sqrt(17.5_gp),sqrt(4.375_gp),-sqrt(9.84375_gp),sqrt(1.09375_gp)/) ! m=5, l=3, n=3
          ftpd(49:56)=rsp*(/sqrt(6.5625_gp),sqrt(6.5625_gp),-sqrt(6.5625_gp),-sqrt(6.5625_gp),sqrt(26.25_gp),&
               -sqrt(26.25_gp),sqrt(6.5625_gp),-sqrt(6.5625_gp)/) ! m=6, l=3, n=3
          ftpd(57:62)=rsp*(/sqrt(26.25_gp),sqrt(105._gp),sqrt(26.25_gp),sqrt(105._gp),sqrt(105._gp),&
               sqrt(26.25_gp)/) ! m=7, l=3, n=3
       end select
    end select
  end subroutine tensor_product_decomposition


  !> Routine to extract the coefficients from the quantum numbers and the operation
  !> it provides the tensor product decomposition of the laplacian of a given shell
  pure subroutine tensor_product_decomposition_laplacian(a,n,l,ntpd_shell,ntpd,pow,ftpd)
    implicit none
    integer, intent(in) :: n,l
    real(gp), intent(in) :: a !<exponent of the gaussian
    integer, intent(out) :: ntpd_shell !< No. of terms for the whole shell
    integer, dimension(2*l+1), intent(out) :: ntpd !< number of terms per shell element
    integer, dimension(3,NTERM_MAX_KINETIC), intent(out) :: pow !< tensor product decompositions
    real(gp), dimension(NTERM_MAX_KINETIC), intent(out) :: ftpd !<factors
    !All data alltogether
    real(gp), parameter :: rsp=0.564189583547756286948079451560772585844050629328998856844086_gp !1/sqrt(pi)
    select case(n)
    case(1)
       select case(l)
       case(0)
          ntpd_shell=4 ! l=0,n=1
          ntpd(1:1)=(/4/) ! l=0,n=1
          pow(1:3,1:4)=reshape((/2,0,0,0,2,0,0,0,2,0,0,0/),(/3,4/)) ! l=0,n=1
          ftpd(1:4)=rsp*(/2._gp*a**2,2._gp*a**2,2._gp*a**2,-3._gp*a/) ! m=1, l=0, n=1
       case(1)
          ntpd_shell=12 ! l=1,n=1
          ntpd(1:3)=(/4,4,4/) ! l=1,n=1
          pow(1:3,1:12)=reshape((/3,0,0,1,2,0,1,0,2,1,0,0,2,1,0,0,3,0,0,1,2,0,1,0,2,0,1,0,2,1,0,0,3,0,0,1/),(/3,12/)) ! l=1,n=1
          ftpd(1:4)=rsp*(/sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(75._gp)*a/) ! m=1, l=1, n=1
          ftpd(5:8)=rsp*(/sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(75._gp)*a/) ! m=2, l=1, n=1
          ftpd(9:12)=rsp*(/sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(75._gp)*a/) ! m=3, l=1, n=1
       case(2)
          ntpd_shell=27 ! l=2,n=1
          ntpd(1:5)=(/4,4,4,6,9/) ! l=2,n=1
          pow(1:3,1:27)=reshape((/2,1,1,0,3,1,0,1,3,0,1,1,3,0,1,1,2,1,1,0,3,1,0,1,3,1,0,1,3,0,1,1,2,1,1,0,&
               4,0,0,0,4,0,2,0,2,0,2,2,2,0,0,0,2,0,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4,2,0,0,0,2,0,0,0,2/),(/3,27/)) ! l=2,n=1
          ftpd(1:4)=rsp*(/sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(735._gp)*a/) ! m=1, l=2, n=1
          ftpd(5:8)=rsp*(/sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(735._gp)*a/) ! m=2, l=2, n=1
          ftpd(9:12)=rsp*(/sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(735._gp)*a/) ! m=3, l=2, n=1
          ftpd(13:18)=rsp*(/sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(183.75_gp)*a,sqrt(183.75_gp)*a/) ! m=4, l=2, n=1
          ftpd(19:27)=rsp*(/-sqrt(5._gp)*a**2,-sqrt(20._gp)*a**2,-sqrt(5._gp)*a**2,sqrt(5._gp)*a**2,&
               sqrt(5._gp)*a**2,sqrt(20._gp)*a**2,sqrt(61.25_gp)*a,sqrt(61.25_gp)*a,-sqrt(245._gp)*a/) ! m=5, l=2, n=1
       case(3)
          ntpd_shell=51 ! l=3,n=1
          ntpd(1:7)=(/9,9,9,7,7,6,4/) ! l=3,n=1
          pow(1:3,1:51)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,3,0,0,1,2,0,1,0,2,4,1,0,2,3,0,0,5,0,&
               2,1,2,0,3,2,0,1,4,2,1,0,0,3,0,0,1,2,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,2,0,1,0,2,1,0,0,3,5,&
               0,0,3,2,0,1,4,0,3,0,2,1,2,2,3,0,0,1,2,0,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,2,1,0,0,3,0,4,0,1,0,4,1,2,&
               0,3,0,2,3,2,0,1,0,2,1,3,1,1,1,3,1,1,1,3,1,1,1/),(/3,51/)) ! l=3,n=1
          ftpd(1:9)=rsp*(/sqrt(10.5_gp)*a**2,sqrt(42._gp)*a**2,sqrt(10.5_gp)*a**2,-sqrt(94.5_gp)*a**2,&
               -sqrt(94.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(212.625_gp)*a,-sqrt(212.625_gp)*a,&
               sqrt(3402._gp)*a/) ! m=1, l=3, n=1
          ftpd(10:18)=rsp*(/sqrt(10.5_gp)*a**2,sqrt(42._gp)*a**2,sqrt(10.5_gp)*a**2,&
               -sqrt(94.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(212.625_gp)*a,&
               -sqrt(212.625_gp)*a,sqrt(3402._gp)*a/) ! m=2, l=3, n=1
          ftpd(19:27)=rsp*(/sqrt(63._gp)*a**2,sqrt(252._gp)*a**2,sqrt(63._gp)*a**2,sqrt(7._gp)*a**2,&
               sqrt(7._gp)*a**2,-sqrt(28._gp)*a**2,-sqrt(1275.75_gp)*a,-sqrt(1275.75_gp)*a,sqrt(567._gp)*a/) ! m=3, l=3, n=1
          ftpd(28:34)=rsp*(/sqrt(17.5_gp)*a**2,-sqrt(70._gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(17.5_gp)*a**2,&
               -sqrt(157.5_gp)*a**2,-sqrt(354.375_gp)*a,sqrt(3189.375_gp)*a/) ! m=4, l=3, n=1
          ftpd(35:41)=rsp*(/-sqrt(157.5_gp)*a**2,-sqrt(70._gp)*a**2,sqrt(17.5_gp)*a**2,-sqrt(157.5_gp)*a**2,&
               sqrt(17.5_gp)*a**2,sqrt(3189.375_gp)*a,-sqrt(354.375_gp)*a/) ! m=5, l=3, n=1
          ftpd(42:47)=rsp*(/sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,&
               -sqrt(2126.25_gp)*a,sqrt(2126.25_gp)*a/) ! m=6, l=3, n=1
          ftpd(48:51)=rsp*(/sqrt(420._gp)*a**2,sqrt(420._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(8505._gp)*a/) ! m=7, l=3, n=1
       end select
    case(2)
       select case(l)
       case(0)
          ntpd_shell=10 ! l=0,n=2
          ntpd(1:1)=(/10/) ! l=0,n=2
          pow(1:3,1:10)=reshape((/0,0,0,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4,2,0,0,0,2,0,0,0,2/),(/3,10/)) ! l=0,n=2
          ftpd(1:10)=rsp*(/3._gp,2._gp*a**2,4._gp*a**2,2._gp*a**2,4._gp*a**2,4._gp*a**2,2._gp*a**2,&
               -7._gp*a,-7._gp*a,-7._gp*a/) ! m=1, l=0, n=2
       case(1)
          ntpd_shell=30 ! l=1,n=2
          ntpd(1:3)=(/10,10,10/) ! l=1,n=2
          pow(1:3,1:30)=reshape((/1,0,0,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,3,0,0,1,2,0,1,0,2,0,1,0,4,1,0,2,3,&
               0,0,5,0,2,1,2,0,3,2,0,1,4,2,1,0,0,3,0,0,1,2,0,0,1,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,2,0,1,0,2,&
               1,0,0,3/),(/3,30/)) ! l=1,n=2
          ftpd(1:10)=rsp*(/sqrt(75._gp),sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,&
               sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(243._gp)*a,-sqrt(243._gp)*a,-sqrt(243._gp)*a/) ! m=1, l=1, n=2
          ftpd(11:20)=rsp*(/sqrt(75._gp),sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,&
               sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(243._gp)*a,-sqrt(243._gp)*a,-sqrt(243._gp)*a/) ! m=2, l=1, n=2
          ftpd(21:30)=rsp*(/sqrt(75._gp),sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,sqrt(48._gp)*a**2,&
               sqrt(48._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(243._gp)*a,-sqrt(243._gp)*a,-sqrt(243._gp)*a/) ! m=3, l=1, n=2
       case(2)
          ntpd_shell=60 ! l=2,n=2
          ntpd(1:5)=(/10,10,10,14,16/) ! l=2,n=2
          pow(1:3,1:60)=reshape((/0,1,1,4,1,1,2,3,1,0,5,1,2,1,3,0,3,3,0,1,5,2,1,1,0,3,1,0,1,3,1,0,1,5,0,1,3,2,1,1,4,1,3,&
               0,3,1,2,3,1,0,5,3,0,1,1,2,1,1,0,3,1,1,0,5,1,0,3,3,0,1,5,0,3,1,2,1,3,2,1,1,4,3,1,0,1,3,0,1,1,2,2,0,0,0,2,&
               0,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,0,4,2,2,0,4,0,2,4,4,0,0,0,4,0,2,0,2,0,2,2,2,0,0,0,2,0,0,0,2,6,0,0,4,2,0,2,&
               4,0,0,6,0,2,0,4,0,2,4,0,0,6,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,60/)) ! l=2,n=2
          ftpd(1:10)=rsp*(/sqrt(735._gp),sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,&
               sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a/) ! m=1, l=2, n=2
          ftpd(11:20)=rsp*(/sqrt(735._gp),sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,&
               sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a/) ! m=2, l=2, n=2
          ftpd(21:30)=rsp*(/sqrt(735._gp),sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,sqrt(240._gp)*a**2,&
               sqrt(240._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a,-sqrt(1815._gp)*a/) ! m=3, l=2, n=2
          ftpd(31:44)=rsp*(/sqrt(183.75_gp),-sqrt(183.75_gp),sqrt(15._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(15._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(60._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(453.75_gp)*a,sqrt(453.75_gp)*a,-sqrt(453.75_gp)*a,sqrt(453.75_gp)*a/) ! m=4, l=2, n=2
          ftpd(45:60)=rsp*(/-sqrt(61.25_gp),-sqrt(61.25_gp),sqrt(245._gp),-sqrt(5._gp)*a**2,-sqrt(45._gp)*a**2,&
               -sqrt(45._gp)*a**2,-sqrt(5._gp)*a**2,sqrt(45._gp)*a**2,sqrt(45._gp)*a**2,sqrt(20._gp)*a**2,&
               sqrt(151.25_gp)*a,sqrt(605._gp)*a,sqrt(151.25_gp)*a,-sqrt(151.25_gp)*a,-sqrt(151.25_gp)*a,&
               -sqrt(605._gp)*a/) ! m=5, l=2, n=2
       case(3)
          ntpd_shell=113 ! l=3,n=2
          ntpd(1:7)=(/19,19,19,16,16,14,10/) ! l=3,n=2
          pow(1:3,1:113)=reshape((/3,0,0,1,2,0,1,0,2,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,1,0,6,&
               5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,2,1,0,0,3,0,0,1,2,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,&
               1,4,0,3,4,0,1,6,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,2,0,1,0,2,1,0,0,3,6,0,1,4,2,1,2,4,1,0,6,1,4,0,&
               3,2,2,3,0,4,3,2,0,5,0,2,5,0,0,7,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,3,0,0,1,2,0,7,0,0,5,2,0,3,4,0,&
               1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,2,1,0,0,3,0,6,1,0,4,3,0,2,5,0,0,&
               7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,2,0,1,0,2,1,6,0,1,4,2,1,2,4,1,0,6,&
               1,4,0,3,0,4,3,2,0,5,0,2,5,4,0,1,0,4,1,2,0,3,0,2,3,1,1,1,5,1,1,3,3,1,1,5,1,3,1,3,1,3,3,1,1,5,3,1,1,&
               1,3,1,1,1,3/),(/3,113/)) ! l=3,n=2
          ftpd(1:19)=rsp*(/sqrt(212.625_gp),sqrt(212.625_gp),-sqrt(3402._gp),sqrt(10.5_gp)*a**2,sqrt(94.5_gp)*a**2,&
               sqrt(94.5_gp)*a**2,sqrt(10.5_gp)*a**2,-sqrt(42._gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(42._gp)*a**2,&
               -sqrt(514.5_gp)*a**2,-sqrt(514.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(443.625_gp)*a,-sqrt(1774.5_gp)*a,&
               -sqrt(443.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(7098._gp)*a/) ! m=1, l=3, n=2
          ftpd(20:38)=rsp*(/sqrt(212.625_gp),sqrt(212.625_gp),-sqrt(3402._gp),sqrt(10.5_gp)*a**2,sqrt(94.5_gp)*a**2,&
               sqrt(94.5_gp)*a**2,sqrt(10.5_gp)*a**2,-sqrt(42._gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(42._gp)*a**2,&
               -sqrt(514.5_gp)*a**2,-sqrt(514.5_gp)*a**2,-sqrt(168._gp)*a**2,-sqrt(443.625_gp)*a,-sqrt(1774.5_gp)*a,&
               -sqrt(443.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(3992.625_gp)*a,sqrt(7098._gp)*a/) ! m=2, l=3, n=2
          ftpd(39:57)=rsp*(/sqrt(1275.75_gp),sqrt(1275.75_gp),-sqrt(567._gp),sqrt(63._gp)*a**2,sqrt(567._gp)*a**2,&
               sqrt(567._gp)*a**2,sqrt(63._gp)*a**2,sqrt(112._gp)*a**2,sqrt(448._gp)*a**2,sqrt(112._gp)*a**2,&
               -sqrt(7._gp)*a**2,-sqrt(7._gp)*a**2,-sqrt(28._gp)*a**2,-sqrt(2661.75_gp)*a,-sqrt(10647._gp)*a,&
               -sqrt(2661.75_gp)*a,-sqrt(295.75_gp)*a,-sqrt(295.75_gp)*a,sqrt(1183._gp)*a/) ! m=3, l=3, n=2
          ftpd(58:73)=rsp*(/sqrt(354.375_gp),-sqrt(3189.375_gp),sqrt(17.5_gp)*a**2,-sqrt(17.5_gp)*a**2,&
               -sqrt(437.5_gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(70._gp)*a**2,-sqrt(280._gp)*a**2,-sqrt(630._gp)*a**2,&
               sqrt(17.5_gp)*a**2,-sqrt(157.5_gp)*a**2,-sqrt(739.375_gp)*a,sqrt(2957.5_gp)*a,sqrt(6654.375_gp)*a,&
               -sqrt(739.375_gp)*a,sqrt(6654.375_gp)*a/) ! m=4, l=3, n=2
          ftpd(74:89)=rsp*(/-sqrt(3189.375_gp),sqrt(354.375_gp),-sqrt(157.5_gp)*a**2,-sqrt(437.5_gp)*a**2,-sqrt(17.5_gp)*a**2,&
               sqrt(17.5_gp)*a**2,-sqrt(630._gp)*a**2,-sqrt(280._gp)*a**2,sqrt(70._gp)*a**2,-sqrt(157.5_gp)*a**2,&
               sqrt(17.5_gp)*a**2,sqrt(6654.375_gp)*a,sqrt(2957.5_gp)*a,-sqrt(739.375_gp)*a,sqrt(6654.375_gp)*a,&
               -sqrt(739.375_gp)*a/) ! m=5, l=3, n=2
          ftpd(90:103)=rsp*(/sqrt(2126.25_gp),-sqrt(2126.25_gp),sqrt(105._gp)*a**2,sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,&
               -sqrt(105._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(420._gp)*a**2,sqrt(105._gp)*a**2,-sqrt(105._gp)*a**2,&
               -sqrt(4436.25_gp)*a,sqrt(4436.25_gp)*a,-sqrt(4436.25_gp)*a,sqrt(4436.25_gp)*a/) ! m=6, l=3, n=2
          ftpd(104:113)=rsp*(/sqrt(8505._gp),sqrt(420._gp)*a**2,sqrt(1680._gp)*a**2,sqrt(420._gp)*a**2,sqrt(1680._gp)*a**2,&
               sqrt(1680._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(17745._gp)*a,-sqrt(17745._gp)*a,-sqrt(17745._gp)*a/) ! m=7, l=3, n=2
       end select
    case(3)
       select case(l)
       case(0)
          ntpd_shell=19 ! l=0,n=3
          ntpd(1:1)=(/19/) ! l=0,n=3
          pow(1:3,1:19)=reshape((/2,0,0,0,2,0,0,0,2,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,2,2,2,0,4,2,2,0,4,0,2,4,0,0,6,4,&
               0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4/),(/3,19/)) ! l=0,n=3
          ftpd(1:19)=rsp*(/10._gp,10._gp,10._gp,2._gp*a**2,6._gp*a**2,6._gp*a**2,2._gp*a**2,6._gp*a**2,12._gp*a**2,&
               6._gp*a**2,6._gp*a**2,6._gp*a**2,2._gp*a**2,-11._gp*a,-22._gp*a,-11._gp*a,-22._gp*a,-22._gp*a,&
               -11._gp*a/) ! m=1, l=0, n=3
       case(1)
          ntpd_shell=57 ! l=1,n=3
          ntpd(1:3)=(/19,19,19/) ! l=1,n=3
          pow(1:3,1:57)=reshape((/3,0,0,1,2,0,1,0,2,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,0,4,1,2,4,1,0,6,5,0,0,&
               3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,2,1,0,0,3,0,0,1,2,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,&
               0,1,6,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,2,0,1,0,2,1,0,0,3,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,2,2,3,0,4,3,2,&
               0,5,0,2,5,0,0,7,4,0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5/),(/3,57/)) ! l=1,n=3
          ftpd(1:19)=rsp*(/sqrt(588._gp),sqrt(588._gp),sqrt(588._gp),sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(432._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,&
               -sqrt(2028._gp)*a,-sqrt(507._gp)*a/) ! m=1, l=1, n=3
          ftpd(20:38)=rsp*(/sqrt(588._gp),sqrt(588._gp),sqrt(588._gp),sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(432._gp)*a**2,sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,&
               -sqrt(2028._gp)*a,-sqrt(507._gp)*a/) ! m=2, l=1, n=3
          ftpd(39:57)=rsp*(/sqrt(588._gp),sqrt(588._gp),sqrt(588._gp),sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,sqrt(108._gp)*a**2,sqrt(432._gp)*a**2,sqrt(108._gp)*a**2,&
               sqrt(108._gp)*a**2,sqrt(108._gp)*a**2,sqrt(12._gp)*a**2,-sqrt(507._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a,&
               -sqrt(2028._gp)*a,-sqrt(2028._gp)*a,-sqrt(507._gp)*a/) ! m=3, l=1, n=3
       case(2)
          ntpd_shell=109 ! l=2,n=3
          ntpd(1:5)=(/19,19,19,24,28/) ! l=2,n=3
          pow(1:3,1:109)=reshape((/2,1,1,0,3,1,0,1,3,6,1,1,4,3,1,2,5,1,0,7,1,4,1,3,2,3,3,0,5,3,2,1,5,0,3,5,0,1,7,&
               4,1,1,2,3,1,0,5,1,2,1,3,0,3,3,0,1,5,3,0,1,1,2,1,1,0,3,7,0,1,5,2,1,3,4,1,1,6,1,5,0,3,3,2,3,1,4,3,3,&
               0,5,1,2,5,1,0,7,5,0,1,3,2,1,1,4,1,3,0,3,1,2,3,1,0,5,3,1,0,1,3,0,1,1,2,7,1,0,5,3,0,3,5,0,1,7,0,5,1,&
               2,3,3,2,1,5,2,3,1,4,1,3,4,1,1,6,5,1,0,3,3,0,1,5,0,3,1,2,1,3,2,1,1,4,4,0,0,0,4,0,2,0,2,0,2,2,8,0,0,&
               6,2,0,2,6,0,0,8,0,6,0,2,4,2,2,2,4,2,0,6,2,4,0,4,0,4,4,2,0,6,0,2,6,6,0,0,4,2,0,2,4,0,0,6,0,4,0,2,0,&
               4,2,2,0,4,0,2,4,4,0,0,2,2,0,0,4,0,2,0,2,0,2,2,0,0,4,8,0,0,6,2,0,4,4,0,2,6,0,0,8,0,6,0,2,4,2,2,2,4,&
               2,0,6,2,4,0,4,2,2,4,0,4,4,2,0,6,0,2,6,0,0,8,6,0,0,4,2,0,2,4,0,0,6,0,2,0,4,0,2,4,0,0,6/),(/3,109/)) ! l=2,n=3
          ftpd(1:19)=rsp*(/sqrt(4860._gp),sqrt(4860._gp),sqrt(4860._gp),sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(2160._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,&
               -sqrt(3375._gp)*a,-sqrt(13500._gp)*a,-sqrt(13500._gp)*a,-sqrt(3375._gp)*a/) ! m=1, l=2, n=3
          ftpd(20:38)=rsp*(/sqrt(4860._gp),sqrt(4860._gp),sqrt(4860._gp),sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(2160._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,&
               -sqrt(3375._gp)*a,-sqrt(13500._gp)*a,-sqrt(13500._gp)*a,-sqrt(3375._gp)*a/) ! m=2, l=2, n=3
          ftpd(39:57)=rsp*(/sqrt(4860._gp),sqrt(4860._gp),sqrt(4860._gp),sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(60._gp)*a**2,sqrt(540._gp)*a**2,sqrt(2160._gp)*a**2,sqrt(540._gp)*a**2,sqrt(540._gp)*a**2,&
               sqrt(540._gp)*a**2,sqrt(60._gp)*a**2,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,-sqrt(3375._gp)*a,-sqrt(13500._gp)*a,&
               -sqrt(13500._gp)*a,-sqrt(3375._gp)*a/) ! m=3, l=2, n=3
          ftpd(58:81)=rsp*(/sqrt(1215._gp),-sqrt(1215._gp),sqrt(1215._gp),-sqrt(1215._gp),sqrt(15._gp)*a**2,sqrt(60._gp)*a**2,&
               -sqrt(60._gp)*a**2,-sqrt(15._gp)*a**2,sqrt(135._gp)*a**2,sqrt(135._gp)*a**2,-sqrt(135._gp)*a**2,&
               -sqrt(135._gp)*a**2,sqrt(135._gp)*a**2,-sqrt(135._gp)*a**2,sqrt(15._gp)*a**2,-sqrt(15._gp)*a**2,&
               -sqrt(843.75_gp)*a,-sqrt(843.75_gp)*a,sqrt(843.75_gp)*a,sqrt(843.75_gp)*a,-sqrt(3375._gp)*a,sqrt(3375._gp)*a,&
               -sqrt(843.75_gp)*a,sqrt(843.75_gp)*a/) ! m=4, l=2, n=3
          ftpd(82:109)=rsp*(/-sqrt(405._gp),-sqrt(1620._gp),-sqrt(405._gp),sqrt(405._gp),sqrt(405._gp),sqrt(1620._gp),&
               -sqrt(5._gp)*a**2,-sqrt(80._gp)*a**2,-sqrt(180._gp)*a**2,-sqrt(80._gp)*a**2,-sqrt(5._gp)*a**2,&
               -sqrt(5._gp)*a**2,-sqrt(45._gp)*a**2,-sqrt(45._gp)*a**2,-sqrt(5._gp)*a**2,sqrt(45._gp)*a**2,&
               sqrt(180._gp)*a**2,sqrt(45._gp)*a**2,sqrt(125._gp)*a**2,sqrt(125._gp)*a**2,sqrt(20._gp)*a**2,&
               sqrt(281.25_gp)*a,sqrt(2531.25_gp)*a,sqrt(2531.25_gp)*a,sqrt(281.25_gp)*a,-sqrt(2531.25_gp)*a,&
               -sqrt(2531.25_gp)*a,-sqrt(1125._gp)*a/) ! m=5, l=2, n=3
       case(3)
          ntpd_shell=190 ! l=3,n=3
          ntpd(1:7)=(/31,31,31,27,27,24,19/) ! l=3,n=3
          pow(1:3,1:190)=reshape((/5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,1,0,4,9,0,0,7,2,0,5,4,0,3,6,0,1,8,0,7,0,2,5,&
               2,2,3,4,2,1,6,2,5,0,4,3,2,4,1,4,4,3,0,6,1,2,6,1,0,8,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,4,2,3,&
               0,4,1,2,4,1,0,6,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,0,1,4,8,1,0,6,3,0,4,5,0,2,7,0,0,9,0,6,1,2,4,3,2,2,5,2,&
               0,7,2,4,1,4,2,3,4,0,5,4,2,1,6,0,3,6,0,1,8,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,0,1,6,4,&
               0,1,2,2,1,0,4,1,2,0,3,0,2,3,0,0,5,8,0,1,6,2,1,4,4,1,2,6,1,0,8,1,6,0,3,4,2,3,2,4,3,0,6,3,4,0,5,2,2,5,0,4,5,2,&
               0,7,0,2,7,0,0,9,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,2,2,3,0,4,3,2,0,5,0,2,5,0,0,7,5,0,0,3,2,0,1,4,0,3,0,2,1,2,2,9,0,&
               0,5,4,0,3,6,0,1,8,0,7,0,2,5,2,2,3,4,2,1,6,2,5,0,4,3,2,4,1,4,4,3,0,6,1,2,6,7,0,0,5,2,0,3,4,0,1,6,0,5,0,2,3,2,2,1,&
               4,2,3,0,4,1,2,4,4,1,0,2,3,0,0,5,0,2,1,2,0,3,2,8,1,0,6,3,0,4,5,0,0,9,0,6,1,2,4,3,2,2,5,2,0,7,2,4,1,4,2,3,4,0,5,4,2,&
               1,6,0,3,6,6,1,0,4,3,0,2,5,0,0,7,0,4,1,2,2,3,2,0,5,2,2,1,4,0,3,4,4,0,1,0,4,1,2,0,3,0,2,3,8,0,1,6,2,1,2,6,1,0,8,1,6,0,&
               3,4,2,3,2,4,3,0,6,3,4,0,5,0,4,5,2,0,7,0,2,7,6,0,1,4,2,1,2,4,1,0,6,1,4,0,3,0,4,3,2,0,5,0,2,5,3,1,1,1,3,1,1,1,3,7,1,1,&
               5,3,1,3,5,1,1,7,1,5,1,3,3,3,3,1,5,3,3,1,5,1,3,5,1,1,7,5,1,1,3,3,1,1,5,1,3,1,3,1,3,3,1,1,5/),(/3,190/)) ! l=3,n=3
          ftpd(1:31)=rsp*(/sqrt(1270.5_gp),sqrt(5082._gp),sqrt(1270.5_gp),-sqrt(11434.5_gp),-sqrt(11434.5_gp),&
               -sqrt(20328._gp),sqrt(10.5_gp)*a**2,sqrt(168._gp)*a**2,sqrt(378._gp)*a**2,sqrt(168._gp)*a**2,&
               sqrt(10.5_gp)*a**2,-sqrt(10.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(10.5_gp)*a**2,&
               -sqrt(850.5_gp)*a**2,-sqrt(3402._gp)*a**2,-sqrt(850.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,&
               -sqrt(168._gp)*a**2,-sqrt(758.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(758.625_gp)*a,&
               sqrt(3034.5_gp)*a,sqrt(12138._gp)*a,sqrt(3034.5_gp)*a,sqrt(37172.625_gp)*a,sqrt(37172.625_gp)*a,&
               sqrt(12138._gp)*a/) ! m=1, l=3, n=3
          ftpd(32:62)=rsp*(/sqrt(1270.5_gp),sqrt(5082._gp),sqrt(1270.5_gp),-sqrt(11434.5_gp),-sqrt(11434.5_gp),&
               -sqrt(20328._gp),sqrt(10.5_gp)*a**2,sqrt(168._gp)*a**2,sqrt(378._gp)*a**2,sqrt(168._gp)*a**2,sqrt(10.5_gp)*a**2,&
               -sqrt(10.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(94.5_gp)*a**2,-sqrt(10.5_gp)*a**2,-sqrt(850.5_gp)*a**2,&
               -sqrt(3402._gp)*a**2,-sqrt(850.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,-sqrt(1270.5_gp)*a**2,-sqrt(168._gp)*a**2,&
               -sqrt(758.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(6827.625_gp)*a,-sqrt(758.625_gp)*a,sqrt(3034.5_gp)*a,&
               sqrt(12138._gp)*a,sqrt(3034.5_gp)*a,sqrt(37172.625_gp)*a,sqrt(37172.625_gp)*a,sqrt(12138._gp)*a/) ! m=2, l=3, n=3
          ftpd(63:93)=rsp*(/sqrt(7623._gp),sqrt(30492._gp),sqrt(7623._gp),sqrt(847._gp),sqrt(847._gp),-sqrt(3388._gp),&
               sqrt(63._gp)*a**2,sqrt(1008._gp)*a**2,sqrt(2268._gp)*a**2,sqrt(1008._gp)*a**2,sqrt(63._gp)*a**2,sqrt(343._gp)*a**2,&
               sqrt(3087._gp)*a**2,sqrt(3087._gp)*a**2,sqrt(343._gp)*a**2,sqrt(63._gp)*a**2,sqrt(252._gp)*a**2,sqrt(63._gp)*a**2,&
               -sqrt(63._gp)*a**2,-sqrt(63._gp)*a**2,-sqrt(28._gp)*a**2,-sqrt(4551.75_gp)*a,-sqrt(40965.75_gp)*a,&
               -sqrt(40965.75_gp)*a,-sqrt(4551.75_gp)*a,-sqrt(8092._gp)*a,-sqrt(32368._gp)*a,-sqrt(8092._gp)*a,&
               sqrt(505.75_gp)*a,sqrt(505.75_gp)*a,&
               sqrt(2023._gp)*a/) ! m=3, l=3, n=3
          ftpd(94:120)=rsp*(/sqrt(2117.5_gp),-sqrt(8470._gp),-sqrt(19057.5_gp),sqrt(2117.5_gp),-sqrt(19057.5_gp),&
               sqrt(17.5_gp)*a**2,-sqrt(630._gp)*a**2,-sqrt(1120._gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(157.5_gp)*a**2,&
               -sqrt(157.5_gp)*a**2,-sqrt(3937.5_gp)*a**2,-sqrt(1417.5_gp)*a**2,sqrt(157.5_gp)*a**2,&
               -sqrt(630._gp)*a**2,-sqrt(1417.5_gp)*a**2,sqrt(17.5_gp)*a**2,-sqrt(157.5_gp)*a**2,&
               -sqrt(1264.375_gp)*a,sqrt(1264.375_gp)*a,sqrt(31609.375_gp)*a,sqrt(11379.375_gp)*a,&
               -sqrt(5057.5_gp)*a,sqrt(20230._gp)*a,sqrt(45517.5_gp)*a,-sqrt(1264.375_gp)*a,sqrt(11379.375_gp)*a/) ! m=4, l=3, n=3
          ftpd(121:147)=rsp*(/-sqrt(19057.5_gp),-sqrt(8470._gp),sqrt(2117.5_gp),-sqrt(19057.5_gp),sqrt(2117.5_gp),&
               -sqrt(157.5_gp)*a**2,-sqrt(1120._gp)*a**2,-sqrt(630._gp)*a**2,sqrt(17.5_gp)*a**2,-sqrt(1417.5_gp)*a**2,&
               -sqrt(3937.5_gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(157.5_gp)*a**2,-sqrt(1417.5_gp)*a**2,-sqrt(630._gp)*a**2,&
               sqrt(157.5_gp)*a**2,-sqrt(157.5_gp)*a**2,sqrt(17.5_gp)*a**2,sqrt(11379.375_gp)*a,sqrt(31609.375_gp)*a,&
               sqrt(1264.375_gp)*a,-sqrt(1264.375_gp)*a,sqrt(45517.5_gp)*a,sqrt(20230._gp)*a,-sqrt(5057.5_gp)*a,&
               sqrt(11379.375_gp)*a,-sqrt(1264.375_gp)*a/) ! m=5, l=3, n=3
          ftpd(148:171)=rsp*(/sqrt(12705._gp),-sqrt(12705._gp),sqrt(12705._gp),-sqrt(12705._gp),sqrt(105._gp)*a**2,&
               sqrt(420._gp)*a**2,-sqrt(420._gp)*a**2,-sqrt(105._gp)*a**2,sqrt(945._gp)*a**2,sqrt(945._gp)*a**2,&
               -sqrt(945._gp)*a**2,-sqrt(945._gp)*a**2,sqrt(945._gp)*a**2,-sqrt(945._gp)*a**2,sqrt(105._gp)*a**2,&
               -sqrt(105._gp)*a**2,-sqrt(7586.25_gp)*a,-sqrt(7586.25_gp)*a,sqrt(7586.25_gp)*a,sqrt(7586.25_gp)*a,&
               -sqrt(30345._gp)*a,sqrt(30345._gp)*a,-sqrt(7586.25_gp)*a,sqrt(7586.25_gp)*a/) ! m=6, l=3, n=3
          ftpd(172:190)=rsp*(/sqrt(50820._gp),sqrt(50820._gp),sqrt(50820._gp),sqrt(420._gp)*a**2,sqrt(3780._gp)*a**2,&
               sqrt(3780._gp)*a**2,sqrt(420._gp)*a**2,sqrt(3780._gp)*a**2,sqrt(15120._gp)*a**2,sqrt(3780._gp)*a**2,&
               sqrt(3780._gp)*a**2,sqrt(3780._gp)*a**2,sqrt(420._gp)*a**2,-sqrt(30345._gp)*a,-sqrt(121380._gp)*a,&
               -sqrt(30345._gp)*a,-sqrt(121380._gp)*a,-sqrt(121380._gp)*a,-sqrt(30345._gp)*a/) ! m=7, l=3, n=3
       end select
    end select
  end subroutine tensor_product_decomposition_laplacian
end module gaussians