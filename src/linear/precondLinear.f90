!> @file
!! Solve the preconditioning equation for the linear version
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Solves the preconditioning equation by conjugate gradient iterations.
!! The equation reads ( kin.energy + cprecr*Id + potentialPrefac*(r-r0)^4 )x=y
!! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method.
!! x is the right hand side on input and the solution on output
subroutine solvePrecondEquation(iproc,nproc,lr,ncplx,ncong,cprecr,&
     hx,hy,hz,kx,ky,kz,x,  rxyzParab, orbs, potentialPrefac, confPotOrder, &
     work_conv, w)

  use module_base
  use module_types
  use locregs
  use locreg_operations

  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: ncong                !> number of CG iterations
  integer, intent(in) :: ncplx                !> real or complex??
  integer, intent(in) :: confPotOrder         
  real(gp), intent(in) :: hx,hy,hz            !> hgrid in x, y and z direction
  real(gp), intent(in) :: cprecr              !> preconditioning constant
  real(gp), intent(in) :: kx,ky,kz            !> kpoints in x, y and z direction
  type(locreg_descriptors), intent(in) :: lr  !> Type describing the localization region
  !> on input: the right hand side of the equation (i.e. y)
  !! on output: the solution of the equation (i.e. x)
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(inout) :: x
  real(kind=8), dimension(3), intent(in) :: rxyzParab !> the center of the confinement potential
  type(orbitals_data), intent(in) :: orbs     !> type describing the orbitals
  real(kind=8) :: potentialPrefac             !> prefactor for the confinement potential
  type(workarrays_quartic_convolutions),intent(inout):: work_conv !< workarrays for the convolutions
  type(workarr_precond),intent(inout) :: w !< workarrays

  ! Local variables
  character(len=*), parameter :: subname='precondition_residue'
  real(gp), dimension(0:7) :: scal
  real(wp) :: rmr_old,rmr_new,alpha,beta
  integer :: icong
  real(wp), dimension(:), allocatable :: b,r,d
  logical:: with_confpot

  call f_routine(id='solvePrecondEquation')

  !!type(workarrays_quartic_convolutions):: work_conv

  !arrays for the CG procedure
  b = f_malloc(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),id='b')
  r = f_malloc(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),id='r')
  d = f_malloc(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),id='d')

  !call allocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,lr%d,w)

  !call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,x,b)
  call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,x,r)

  with_confpot=(potentialPrefac/=0.d0)
  !!call init_local_work_arrays(lr%d%n1, lr%d%n2, lr%d%n3, &
  !!     lr%d%nfl1, lr%d%nfu1, lr%d%nfl2, lr%d%nfu2, lr%d%nfl3, lr%d%nfu3, &
  !!     with_confpot, work_conv)
  call zero_local_work_arrays(lr%d%n1, lr%d%n2, lr%d%n3, &
       lr%d%nfl1, lr%d%nfu1, lr%d%nfl2, lr%d%nfu2, lr%d%nfl3, lr%d%nfu3, &
       with_confpot, work_conv, subname)
  !!call allocate_workarrays_quartic_convolutions(lr, subname, work_conv)
  call differentiateBetweenBoundaryConditions(iproc,nproc,ncplx,lr,hx,hy,hz,kx,ky,kz,cprecr,x,d,w,scal,&
       rxyzParab, potentialPrefac, confPotOrder, work_conv)

!!  rmr_new=dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),d(1),1,d(1),1)
!!  write(*,*)'debug1',rmr_new

  !this operation should be rewritten in a better way
  !r=b-d ! r=b-Ax
  ! Rewritten using axpy since precondition_preconditioner is now called with r instead of b
  call axpy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, -1.d0, d(1), 1, r(1), 1)

  call calculate_rmr_new(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,r,d,rmr_new)
  !stands for
  !d=r
  !rmr_new=dot_product(r,r)


  do icong=1,ncong 
     !write(*,*)icong,rmr_new
     call differentiateBetweenBoundaryConditions(iproc,nproc,ncplx,lr,hx,hy,hz,kx,ky,kz,cprecr,d,b,w,scal,&
          rxyzParab, potentialPrefac, confPotOrder, work_conv)

     !in the complex case these objects are to be supposed real
     ! 0/0 here!
     alpha=rmr_new/dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),d(1),1,b(1),1)

     call axpy(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),alpha,d(1),1,x(1),1)
     call axpy(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),-alpha,b(1),1,r(1),1)

     if (icong==ncong) exit

     rmr_old=rmr_new

     call calculate_rmr_new(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,r,b,rmr_new)

     beta=rmr_new/rmr_old

     !d=b+beta*d
     call swap(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, b(1), 1, d(1), 1)
     call axpy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, beta, b(1), 1, d(1), 1)
    
  enddo

  call finalise_precond_residue(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,x)

  !!call deallocate_workarrays_quartic_convolutions(work_conv)

  call f_free(b)
  call f_free(r)
  call f_free(d)

  !call timing(iproc,'deallocprec','ON') ! lr408t
  !call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)
  !call timing(iproc,'deallocprec','OF') ! lr408t

  call f_release_routine()

END SUBROUTINE solvePrecondEquation


subroutine differentiateBetweenBoundaryConditions(iproc,nproc,ncplx,lr,hx,hy,hz,kx,ky,kz,&
     cprecr,x,y,w,scal, rxyzParab, parabPrefac, confPotOrder, work_conv)! y:=Ax
  use module_base
  use locregs
  use locreg_operations
  implicit none
  integer, intent(in) :: iproc,nproc,ncplx
  real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(in) ::  x
  type(workarr_precond), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(out) ::  y
  real(8),dimension(3),intent(in):: rxyzParab
  real(8):: parabPrefac
  integer:: confPotOrder
  type(workarrays_quartic_convolutions),intent(inout):: work_conv
  !local variables
  integer :: idx,nf

  call f_routine(id='differentiateBetweenBoundaryConditions')

  if (lr%geocode == 'F') then
     do idx=1,ncplx
        call applyOperator(iproc,nproc,lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3, lr%ns1, lr%ns2, lr%ns3, &
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keygloc,lr%wfd%keyvloc,&
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),&
             lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
             scal,cprecr,hx,hy,hz,&
             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f,&
             x(1,idx),x(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
             y(1,idx),y(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
             rxyzParab, parabPrefac, confPotOrder, &
             w%xpsig_c,w%xpsig_f,w%ypsig_c,w%ypsig_f,&
             w%x_f1,w%x_f2,w%x_f3, work_conv)
     end do
  else if (lr%geocode == 'P') then
     if (lr%hybrid_on) then

        nf=(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1)
        do idx=1,ncplx
           call apply_hp_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc,lr%wfd%keyvloc, &
                cprecr,hx,hy,hz,x(1,idx),y(1,idx),&
                w%x_f,w%x_c,w%x_f1,w%x_f2,w%x_f3,w%y_f,w%z1,&
                lr%d%nfl1,lr%d%nfl2,lr%d%nfl3,lr%d%nfu1,lr%d%nfu2,lr%d%nfu3,nf,&
                lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
        end do
     else
        if (ncplx == 1) then
           call apply_hp_scal(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
                lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
                cprecr,hx,hy,hz,x,y,w%psifscf,w%ww,w%modul1,w%modul2,w%modul3,&
                w%af,w%bf,w%cf,w%ef,scal) 
        else
           call apply_hp_per_k(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
                lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
                !cprecr,hx,hy,hz,0.0_gp,0.0_gp,0.0_gp,x,y,w%psifscf,w%ww,scal) 
                cprecr,hx,hy,hz,kx,ky,kz,x,y,w%psifscf,w%ww,scal) 
        end if
     end if
  else if (lr%geocode == 'S') then
     if (ncplx == 1) then
        call apply_hp_slab_sd(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
             lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
             cprecr,x,y,w%psifscf,w%ww,w%modul1,w%modul3,&
             w%af,w%bf,w%cf,w%ef)
     else
        call apply_hp_slab_k(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
             lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
             cprecr,hx,hy,hz,kx,ky,kz,x,y,w%psifscf,w%ww) 

     end if
   end if

  call f_release_routine()

END SUBROUTINE differentiateBetweenBoundaryConditions



!> WRONG DESCRIPTION
!! This subroutine uncompresses the wave function, applies the operators on it, 
!! and compresses it again. The operators are: kinetic energy + cprec*Id + r^4.
!! Here cprecr is a constant and r^4 is the confinement potential of the form
!! lin%potentialPrefac*[(x-x0)^4 + (y-y0)^4 + (z-z0)^4]
subroutine applyOperator(iproc,nproc,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, ns1, ns2, ns3, &
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hx,hy,hz,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
     rxyzParab, parabPrefac, confPotOrder, &
     xpsig_c,xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3, work_conv)

  use module_base
  use module_types
  use module_interfaces
  use locreg_operations

  implicit none
  integer, intent(in) :: iproc, nproc, n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, ns1, ns2, ns3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
  real(wp), intent(inout) :: cprecr
  real(gp), intent(in) :: hx, hy, hz
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: xpsi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: xpsi_f
  real(wp), dimension(nvctr_c), intent(out) :: ypsi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: ypsi_f
  real(8),dimension(3),intent(in):: rxyzParab
  real(8):: parabPrefac
  type(workarrays_quartic_convolutions),intent(inout):: work_conv

  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xpsig_c,ypsig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xpsig_f,ypsig_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(inout) :: x_f3


  ! Local variables
  character(len=*),parameter:: subname='applyOperator'


  call f_routine(id='applyOperator')

  ! Uncompress the wavefunction.
  call f_routine(id='call_to_uncompress_for_quartic_convolutions')
  call uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, &
       nseg_c, nvctr_c, keyg_c, keyv_c, nseg_f, nvctr_f,  keyg_f, keyv_f, &
       scal, xpsi_c, xpsi_f, &
       work_conv)
  call f_release_routine()

  ! Apply the  following operators to the wavefunctions: kinetic energy + cprec*Id + r^4.
  if(confPotOrder==4) then
     call timing(iproc,'convolQuartic ','ON')
      call f_routine(id='call_to_ConvolQuartic4')
      call ConvolQuartic4(iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, &
           hx, hy, hz, ns1, ns2, ns3, ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
           rxyzParab, parabPrefac, .true., cprecr, max(n1,n2,n3), &
           work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
           work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
           work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
           work_conv%aeff0array, work_conv%beff0array, work_conv%ceff0array, work_conv%eeff0array, &
           work_conv%aeff0_2array, work_conv%beff0_2array, work_conv%ceff0_2array, work_conv%eeff0_2array, &
           work_conv%aeff0_2auxarray, work_conv%beff0_2auxarray, work_conv%ceff0_2auxarray, work_conv%eeff0_2auxarray, &
           work_conv%xya_c, work_conv%xyc_c, &
           work_conv%xza_c, work_conv%xzc_c, &
           work_conv%yza_c, work_conv%yzc_c, &
           work_conv%xya_f, work_conv%xyb_f, work_conv%xyc_f, work_conv%xye_f, &
           work_conv%xza_f, work_conv%xzb_f, work_conv%xzc_f, work_conv%xze_f, &
           work_conv%yza_f, work_conv%yzb_f, work_conv%yzc_f, work_conv%yze_f, &
!           work_conv%aeff0, work_conv%aeff1, work_conv%aeff2, work_conv%aeff3, &
!           work_conv%beff0, work_conv%beff1, work_conv%beff2, work_conv%beff3, &
!           work_conv%ceff0, work_conv%ceff1, work_conv%ceff2, work_conv%ceff3, &
!           work_conv%eeff0, work_conv%eeff1, work_conv%eeff2, work_conv%eeff3, &
!           work_conv%aeff0_2, work_conv%aeff1_2, work_conv%aeff2_2, work_conv%aeff3_2, &
!           work_conv%beff0_2, work_conv%beff1_2, work_conv%beff2_2, work_conv%beff3_2, &
!           work_conv%ceff0_2, work_conv%ceff1_2, work_conv%ceff2_2, work_conv%ceff3_2, &
!           work_conv%eeff0_2, work_conv%eeff1_2, work_conv%eeff2_2, work_conv%eeff3_2, & 
           work_conv%y_c, work_conv%y_f)
      call f_release_routine()
      call timing(iproc,'convolQuartic ','OF')
  else if(confPotOrder==6) then

      stop 'sextic potential deprecated'

  end if

  ! Compress the wavefunctions.
  call f_routine(id='compress_forstandard')
  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,work_conv%y_c,work_conv%y_f,ypsi_c,ypsi_f)
  call f_release_routine()

  call f_release_routine()

END SUBROUTINE applyOperator


