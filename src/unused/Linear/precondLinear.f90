!!!subroutine choosePreconditioner(iproc,nproc,orbs,lin,lr,hx,hy,hz,ncong,hpsi, nat, rxyz, at, it)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Preconditions all orbitals belonging to iproc.
!!!!  
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc     process ID
!!!!     nproc     total number of processes
!!!!     orbs      type describing the physical orbitals psi
!!!!     lin       type containing parameters for the linear version
!!!!     lr        type describing the localization region
!!!!     hx        grid spacing in x direction
!!!!     hy        grid spacing in y direction
!!!!     hz        grid spacing in z direction
!!!!     ncong     number of CG iterations 
!!!!     rxyz      the center of the confinement potential
!!!!     at        type containing the paraneters for the atoms
!!!!     it        iteration -- delete maybe??
!!!!  Input/Output arguments
!!!!  ---------------------
!!!!     hpsi      the gradient to be preconditioned
!!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!integer, intent(in) :: iproc,nproc,ncong
!!!real(gp), intent(in) :: hx,hy,hz
!!!type(locreg_descriptors), intent(in) :: lr
!!!type(orbitals_data), intent(in) :: orbs
!!!type(linearParameters):: lin
!!!real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: hpsi
!!!integer,intent(in):: nat, it
!!!real(8),dimension(3,nat),intent(in):: rxyz
!!!type(atoms_data), intent(in) :: at
!!!!local variables
!!!integer :: iorb, inds, ncplx, ikpt, ierr, iiAt
!!!real(wp) :: cprecr,scpr,eval_zero,evalmax 
!!!real(gp) :: kx,ky,kz
!!!
!!!
!!!
!!!   evalmax=orbs%eval(orbs%isorb+1)
!!!   do iorb=1,orbs%norbp
!!!     evalmax=max(orbs%eval(orbs%isorb+iorb),evalmax)
!!!   enddo
!!!   call MPI_ALLREDUCE(evalmax,eval_zero,1,mpidtypd,&
!!!        MPI_MAX,MPI_COMM_WORLD,ierr)
!!!
!!!
!!!  do iorb=1,orbs%norbp
!!!!     ! define zero energy for preconditioning 
!!!!     eval_zero=max(orbs%eval(orbs%norb),0.d0)  !  Non-spin pol
!!!!     if (orbs%spinsgn(orbs%isorb+iorb) > 0.0_gp) then    !spin-pol
!!!!        eval_zero=max(orbs%eval(orbs%norbu),0.d0)  !up orbital
!!!!     else if (orbs%spinsgn(orbs%isorb+iorb) < 0.0_gp) then
!!!!        eval_zero=max(orbs%eval(orbs%norbu+orbs%norbd),0.d0)  !down orbital
!!!!     end if
!!!     !indo=(iorb-1)*nspinor+1
!!!     !loop over the spinorial components
!!!     !k-point values, if present
!!!     kx=orbs%kpts(1,orbs%iokpt(iorb))
!!!     ky=orbs%kpts(2,orbs%iokpt(iorb))
!!!     kz=orbs%kpts(3,orbs%iokpt(iorb))
!!!
!!!     !real k-point different from Gamma still not implemented
!!!     if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2 ) then
!!!        ncplx=2
!!!     else
!!!        ncplx=1
!!!     end if
!!!
!!!     do inds=1,orbs%nspinor,ncplx
!!!
!!!
!!!       if (.true.) then
!!!           select case(lr%geocode)
!!!           case('F')
!!!              cprecr=sqrt(.2d0**2+min(0.d0,orbs%eval(orbs%isorb+iorb))**2)
!!!           case('S')
!!!              cprecr=sqrt(0.2d0**2+(orbs%eval(orbs%isorb+iorb)-eval_zero)**2)
!!!           case('P')
!!!              cprecr=sqrt(0.2d0**2+(orbs%eval(orbs%isorb+iorb)-eval_zero)**2)
!!!           end select
!!!
!!!           !cases with no CG iterations, diagonal preconditioning
!!!           !for Free BC it is incorporated in the standard procedure
!!!           if (ncong == 0 .and. lr%geocode /= 'F') then
!!!              select case(lr%geocode)
!!!              case('F')
!!!              case('S')
!!!                 call prec_fft_slab(lr%d%n1,lr%d%n2,lr%d%n3, &
!!!                      lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
!!!                      lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyv, &
!!!                      cprecr,hx,hy,hz,hpsi(1,inds,iorb))
!!!              case('P')
!!!                 call prec_fft(lr%d%n1,lr%d%n2,lr%d%n3, &
!!!                      lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!!                      lr%wfd%keygloc,lr%wfd%keyv, &
!!!                      cprecr,hx,hy,hz,hpsi(1,inds,iorb))
!!!              end select
!!!
!!!           else !normal preconditioner
!!!              
!!!              ! iiAt indicates on which atom orbital iorb is centered.
!!!              !iiAt=lin%onWhichAtom(iorb)
!!!              !iiAt=lin%orbs%inWhichLocregp(iorb)
!!!              iiAt=lin%orbs%inWhichLocreg(lin%orbs%isorb+iorb)
!!!              call solvePrecondEquation(lr,ncplx,ncong,cprecr,&
!!!                   hx,hy,hz,kx,ky,kz,hpsi(1,inds,iorb), rxyz(1,iiAt), orbs,&
!!!                   lin%potentialPrefac(at%iatype(iiAt)), lin%confPotOrder, it)
!!!
!!!           end if
!!!
!!!        end if
!!!
!!!     end do
!!!  enddo
!!!
!!!END SUBROUTINE choosePreconditioner

!!!subroutine solvePrecondEquation2(lr,ncplx,ncong,cprecr,&
!!!     hx,hy,hz,kx,ky,kz,x,  rxyzParab, orbs, potentialPrefac, confPotOrder, it)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Solves the preconditioning equation by conjugate gradient iterations.
!!!!   The equation reads ( kin.energy + cprecr*Id + potentialPrefac*(r-r0)^4 )x=y
!!!!   Solves (KE+cprecr*I)*xx=yy by conjugate gradient method.
!!!! 
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------   
!!!!     lr               type describing the localization region
!!!!     ncplx            real or complex??
!!!!     ncong            number of CG iterations
!!!!     cprecr           preconditioning constant
!!!!     hx               hgrid in x direction
!!!!     hy               hgrid in y direction
!!!!     hz               hgrid in z direction
!!!!     kx               kpoints in x direction?
!!!!     ky               kpoints in y direction?
!!!!     kz               kpoints in z direction?
!!!!     rxyzParab        the center of the confinement potential
!!!!     orbs             type describing the orbitals
!!!!     potentialPrefac  prefactor for the confinement potential
!!!!     it               delete later??
!!!!   Input / Output arguments:
!!!!   -------------------------
!!!!     x                on input: the right hand side of the equation (i.e. y)
!!!!                      on output: the solution of the equation (i.e. x)
!!!!
!!!use module_base
!!!use module_types
!!!! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!!!! x is the right hand side on input and the solution on output
!!!! Calling arguments
!!!implicit none
!!!integer, intent(in) :: ncong,ncplx
!!!real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
!!!type(locreg_descriptors), intent(in) :: lr
!!!real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(inout) :: x
!!!real(8),dimension(3),intent(in):: rxyzParab
!!!type(orbitals_data), intent(in):: orbs
!!!real(8):: potentialPrefac
!!!integer:: confPotOrder, it
!!!
!!!! Local variables
!!!character(len=*), parameter :: subname='precondition_residue'
!!!real(gp), dimension(0:7) :: scal
!!!real(wp) :: rmr_old,rmr_new,alpha,beta
!!!integer :: i_stat,i_all,icong
!!!type(workarr_precond) :: w
!!!real(wp), dimension(:), allocatable :: b,r,d
!!!
!!!  !arrays for the CG procedure
!!!  allocate(b(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)+ndebug),stat=i_stat)
!!!  call memocc(i_stat,b,'b',subname)
!!!  allocate(r(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)+ndebug),stat=i_stat)
!!!  call memocc(i_stat,r,'r',subname)
!!!  allocate(d(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)+ndebug),stat=i_stat)
!!!  call memocc(i_stat,d,'d',subname)
!!!
!!!  call allocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,lr%d,w)
!!!
!!!  call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,x,b)
!!!
!!!  call differentiateBetweenBoundaryConditions(ncplx,lr,hx,hy,hz,kx,ky,kz,cprecr,x,d,w,scal,&
!!!       rxyzParab, orbs, potentialPrefac, confPotOrder, it)
!!!
!!!!!  rmr_new=dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),d(1),1,d(1),1)
!!!!!  write(*,*)'debug1',rmr_new
!!!
!!!  !this operation should be rewritten in a better way
!!!  r=b-d ! r=b-Ax
!!!
!!!  call calculate_rmr_new(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,r,d,rmr_new)
!!!  !stands for
!!!  !d=r
!!!  !rmr_new=dot_product(r,r)
!!!
!!!
!!!  do icong=1,ncong 
!!!     !write(*,*)icong,rmr_new
!!!
!!!     call differentiateBetweenBoundaryConditions(ncplx,lr,hx,hy,hz,kx,ky,kz,cprecr,d,b,w,scal,&
!!!          rxyzParab, orbs, potentialPrefac, confPotOrder, it)
!!!
!!!     !in the complex case these objects are to be supposed real
!!!     alpha=rmr_new/dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),d(1),1,b(1),1)
!!!
!!!     call axpy(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),alpha,d(1),1,x(1),1)
!!!     call axpy(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),-alpha,b(1),1,r(1),1)
!!!
!!!     if (icong==ncong) exit
!!!
!!!     rmr_old=rmr_new
!!!
!!!     call calculate_rmr_new(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,r,b,rmr_new)
!!!
!!!     beta=rmr_new/rmr_old
!!!
!!!     d=b+beta*d
!!!    
!!!  enddo
!!!
!!!  call finalise_precond_residue(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,x)
!!!
!!!  i_all=-product(shape(b))*kind(b)
!!!  deallocate(b,stat=i_stat)
!!!  call memocc(i_stat,i_all,'b',subname)
!!!  i_all=-product(shape(r))*kind(r)
!!!  deallocate(r,stat=i_stat)
!!!  call memocc(i_stat,i_all,'r',subname)
!!!  i_all=-product(shape(d))*kind(d)
!!!  deallocate(d,stat=i_stat)
!!!  call memocc(i_stat,i_all,'d',subname)
!!!
!!!  call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)
!!!
!!!END SUBROUTINE solvePrecondEquation2

!!!subroutine differentiateBetweenBoundaryConditions2(ncplx,lr,hx,hy,hz,kx,ky,kz,&
!!!     cprecr,x,y,w,scal, rxyzParab, orbs, parabPrefac, confPotOrder, it)! y:=Ax
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  integer, intent(in) :: ncplx
!!!  real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
!!!  type(locreg_descriptors), intent(in) :: lr
!!!  real(gp), dimension(0:7), intent(in) :: scal
!!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(in) ::  x
!!!  type(workarr_precond), intent(inout) :: w
!!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(out) ::  y
!!!real(8),dimension(3),intent(in):: rxyzParab
!!!type(orbitals_data), intent(in) :: orbs
!!!real(8):: parabPrefac
!!!integer:: confPotOrder, it
!!!  !local variables
!!!  integer :: idx,nf
!!!
!!!  if (lr%geocode == 'F') then
!!!     do idx=1,ncplx
!!!        call applyOperator(lr%d%n1,lr%d%n2,lr%d%n3,&
!!!             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3, &
!!!             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyg,lr%wfd%keyv,&
!!!             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!!             lr%wfd%keyg(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),&
!!!             lr%wfd%keyv(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!!             scal,cprecr,hx,&
!!!             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
!!!             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f,&
!!!             x(1,idx),x(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
!!!             y(1,idx),y(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
!!!             w%xpsig_c,w%xpsig_f,w%ypsig_c,w%ypsig_f,&
!!!             w%x_f1,w%x_f2,w%x_f3, rxyzParab, orbs, lr, parabPrefac, confPotOrder, it)
!!!     end do
!!!  else if (lr%geocode == 'P') then
!!!     if (lr%hybrid_on) then
!!!
!!!        nf=(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1)
!!!        do idx=1,ncplx
!!!           call apply_hp_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
!!!                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!!                lr%wfd%keyg,lr%wfd%keyv, &
!!!                cprecr,hx,hy,hz,x(1,idx),y(1,idx),&
!!!                w%x_f,w%x_c,w%x_f1,w%x_f2,w%x_f3,w%y_f,w%z1,&
!!!                lr%d%nfl1,lr%d%nfl2,lr%d%nfl3,lr%d%nfu1,lr%d%nfu2,lr%d%nfu3,nf,&
!!!                lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
!!!        end do
!!!     else
!!!        if (ncplx == 1) then
!!!           call apply_hp_scal(lr%d%n1,lr%d%n2,lr%d%n3,&
!!!                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
!!!                lr%wfd%nvctr_f,lr%wfd%keyg,lr%wfd%keyv, &
!!!                cprecr,hx,hy,hz,x,y,w%psifscf,w%ww,w%modul1,w%modul2,w%modul3,&
!!!                w%af,w%bf,w%cf,w%ef,scal) 
!!!        else
!!!           call apply_hp_per_k(lr%d%n1,lr%d%n2,lr%d%n3,&
!!!                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
!!!                lr%wfd%nvctr_f,lr%wfd%keyg,lr%wfd%keyv, &
!!!                !cprecr,hx,hy,hz,0.0_gp,0.0_gp,0.0_gp,x,y,w%psifscf,w%ww,scal) 
!!!                cprecr,hx,hy,hz,kx,ky,kz,x,y,w%psifscf,w%ww,scal) 
!!!        end if
!!!     end if
!!!  else if (lr%geocode == 'S') then
!!!     if (ncplx == 1) then
!!!        call apply_hp_slab_sd(lr%d%n1,lr%d%n2,lr%d%n3,&
!!!             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
!!!             lr%wfd%nvctr_f,lr%wfd%keyg,lr%wfd%keyv, &
!!!             cprecr,hx,hy,hz,x,y,w%psifscf,w%ww,w%modul1,w%modul3,&
!!!             w%af,w%bf,w%cf,w%ef)
!!!     else
!!!        call apply_hp_slab_k(lr%d%n1,lr%d%n2,lr%d%n3,&
!!!             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
!!!             lr%wfd%nvctr_f,lr%wfd%keyg,lr%wfd%keyv, &
!!!             cprecr,hx,hy,hz,kx,ky,kz,x,y,w%psifscf,w%ww) 
!!!
!!!     end if
!!!   end if
!!!END SUBROUTINE differentiateBetweenBoundaryConditions2
!!!
!!!
!!!
!!!subroutine applyOperator2(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
!!!     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
!!!     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
!!!     xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
!!!     xpsig_c,xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3, rxyzParab, orbs, lr, parabPrefac, confPotOrder, it)
!!!!
!!!! Purpose:
!!!! ========
!!!!   This subroutine uncompresses the wave function, applies the operators on it, 
!!!!   and compresses it again. The operators are: kinetic energy + cprec*Id + r^4.
!!!!   Here cprecr is a consatnt and r^4 is the confinement potential of the form
!!!!   lin%potentialPrefac*[(x-x0)^4 + (y-y0)^4 + (z-z0)^4]
!!!!
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
!!!  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
!!!  real(wp), intent(in) :: cprecr
!!!  real(gp), intent(in) :: hgrid
!!!  integer, dimension(nseg_c), intent(in) :: keyv_c
!!!  integer, dimension(nseg_f), intent(in) :: keyv_f
!!!  integer, dimension(2,nseg_c), intent(in) :: keyg_c
!!!  integer, dimension(2,nseg_f), intent(in) :: keyg_f
!!!  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
!!!  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
!!!  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
!!!  real(wp), dimension(0:3), intent(in) :: scal
!!!  real(wp), dimension(nvctr_c), intent(in) :: xpsi_c
!!!  real(wp), dimension(7,nvctr_f), intent(in) :: xpsi_f
!!!  real(wp), dimension(nvctr_c), intent(out) :: ypsi_c
!!!  real(wp), dimension(7,nvctr_f), intent(out) :: ypsi_f
!!!  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xpsig_c,ypsig_c
!!!  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xpsig_f,ypsig_f
!!!  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x_f1
!!!  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(inout) :: x_f2
!!!  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(inout) :: x_f3
!!!  real(8),dimension(3),intent(in):: rxyzParab
!!!  type(orbitals_data), intent(in) :: orbs
!!!  type(locreg_descriptors), intent(in) :: lr
!!!  real(8):: parabPrefac
!!!  integer:: it
!!!
!!!
!!!  ! Uncompress the wavefunction.
!!!  call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
!!!       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
!!!       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
!!!       scal,xpsi_c,xpsi_f,xpsig_c,xpsig_f,x_f1,x_f2,x_f3)
!!!
!!!  ! Apply the  following operators to the wavefunctions: kinetic energy + cprec*Id + r^4.
!!!  if(confPotOrder==4) then
!!!      call ConvolkineticQuartic(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
!!!           cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,&
!!!           xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3, rxyzParab(1), parabPrefac, it)
!!!  else if(confPotOrder==6) then
!!!      call ConvolkineticSextic(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
!!!           cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,&
!!!           xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3, rxyzParab(1), parabPrefac, it)
!!!  end if
!!!
!!!  ! Compress the wavefunctions.
!!!  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
!!!       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
!!!       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
!!!       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)
!!!
!!!
!!!END SUBROUTINE applyOperator2

