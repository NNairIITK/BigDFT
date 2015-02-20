!! module forinterpolation between structures
!! techniques implemented:
!! linear synchronous transit:
!!     subroutine lstpthpnt(nat,rxyzR,rxyzP,lambda,rxyz)
!!     
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> @file
!> minor module for wrapping the lst_penalty routine
!! intended as a temporary patch to conceptual problem in using transforce
!! routine
module lst_penalty_wrapper
  use module_defs, only: gp
  implicit none  

  !usage of pointers, just for fun
!  real(gp), pointer, private :: lambda_ptr=>null()
!  real(gp), dimension(:,:), pointer, private :: rxyzR_ptr=>null(),rxyzP_ptr=>null()
  real(gp), private :: lambda_int
  real(gp), dimension(:,:), allocatable, private :: rxyzR_int,rxyzP_int
 
contains

subroutine init_lst_wrapper(nat,rxyzR,rxyzP,lambda)
    use module_base
    implicit none
    integer, intent(in) :: nat
!    real(gp), intent(in), target :: rxyzR(3,nat) !reactant
!    real(gp), intent(in), target :: rxyzP(3,nat) !product
!    real(gp), intent(in), target :: lambda !interpolation parameter
    real(gp), intent(in), target :: rxyzR(3,nat) !reactant
    real(gp), intent(in), target :: rxyzP(3,nat) !product
    real(gp), intent(in), target :: lambda !interpolation parameter

!pointer lead to segault after some hours
!(intel ifort 14)
!    rxyzR_ptr => rxyzR
!    rxyzP_ptr => rxyzP
!    lambda_ptr=> lambda
    rxyzR_int = f_malloc((/ 1.to.3, 1.to.nat/),id='rxyzR_int')
    rxyzP_int = f_malloc((/ 1.to.3, 1.to.nat/),id='rxyzP_int')

    rxyzR_int = rxyzR
    rxyzP_int = rxyzP
    lambda_int = lambda

end subroutine init_lst_wrapper
!=====================================================================
subroutine finalize_lst_wrapper()
    use module_base
    implicit none

    call f_free(rxyzR_int)
    call f_free(rxyzP_int)
end subroutine finalize_lst_wrapper
!=====================================================================
subroutine valforce(nat,rat,fat,epot)
     !wrapper function for lst_penalty
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp), intent(in) :: rat(3,nat)
    real(gp), intent(out) :: fat(3,nat)
    real(gp), intent(out) :: epot
    !local variables
!    call lst_penalty(nat,rxyzR_ptr,rxyzP_ptr,rat,lambda_ptr,epot,fat)
    call lst_penalty(nat,rxyzR_int,rxyzP_int,rat,lambda_int,epot,fat)
end subroutine valforce

!=====================================================================
subroutine lst_penalty(nat,rxyzR,rxyzP,rxyz,lambda,val,force)
    !computes the linear synchronous penalty function
    !and gradient
    !
    !see:
    !
    !Halgren, T. A., & Lipscomb, W. N. (1977). The synchronous-transit
    !method for determining reaction pathways and locating molecular
    !transition states. Chemical Physics Letters, 49(2), 225–232.
    !doi:10.1016/0009-2614(77)80574-5
    !
    !and
    !
    !Behn, A., Zimmerman, P. M., Bell, A. T., & Head-Gordon, M. (2011)
    !Incorporating Linear Synchronous Transit Interpolation into
    !the Growing String Method: Algorithm and Applications.
    !Journal of Chemical Theory and Computation, 7(12), 4019–4025.
    !doi:10.1021/ct200654u
    !
    !and
    !
    !BS: For computation of gradient see also personal notes in intel
    !notebook (paper version) from June 18th, 2014
    !
    use module_base
    implicit none
    !parameters
    integer, intent(in)  :: nat
    real(gp), intent(in) :: rxyzR(3,nat) !reactant
    real(gp), intent(in) :: rxyzP(3,nat) !product
    real(gp), intent(in) :: rxyz(3,nat)  !the positon at which
    !the penalty fct. is to
    !be evaluated
    real(gp), intent(in)  :: lambda !interpolation parameter
    !(0<=lambda<=1)
    real(gp), intent(out) :: val  !the value of the penalty function
    real(gp), intent(out) :: force(3,nat) !the negative gradient of
    !the penal. fct. at
    !rxyz
    !internal
    integer :: b !outer loop
    integer :: a !inner loop
    real(gp) :: rabi !interpolated interatomic distances
    real(gp) :: rabim4 !rim4=ri**(-4)
    real(gp) :: rabC !computed interatomic distances
    real(gp) :: rabR !interatomic dist. of reactant
    real(gp) :: rabP !interatomic dist. of product
    real(gp) :: oml
    real(gp) :: rxRb, ryRb, rzRb
    real(gp) :: rxPb, ryPb, rzPb
    real(gp) :: rxb, ryb, rzb
    real(gp) :: rxd, ryd, rzd
    real(gp) :: rabiMrabC
    real(gp) :: tt, ttx, tty, ttz
    real(gp) :: waxi, wayi, wazi 

    oml=1.0_gp-lambda!one minus lambda
    val=0.0_gp!function value
    force=0.0_gp!negative gradient
    !first sum
    do b = 1, nat-1
       rxRb = rxyzR(1,b)
       ryRb = rxyzR(2,b)
       rzRb = rxyzR(3,b)
       rxPb = rxyzP(1,b)
       ryPb = rxyzP(2,b)
       rzPb = rxyzP(3,b)
       rxb = rxyz(1,b)
       ryb = rxyz(2,b)
       rzb = rxyz(3,b)
       do a = b+1, nat
          !compute interatomic distances of reactant
          rabR = (rxyzR(1,a)-rxRb)**2+&
               (rxyzR(2,a)-ryRb)**2+&
               (rxyzR(3,a)-rzRb)**2
          rabR = sqrt(rabR)
          !compute interatomic distances of product
          rabP = (rxyzP(1,a)-rxPb)**2+&
               (rxyzP(2,a)-ryPb)**2+&
               (rxyzP(3,a)-rzPb)**2
          rabP = sqrt(rabP)
          !compute interpolated interatomic distances
          rabi = oml*rabR+lambda*rabP
          rabim4 = 1.0_gp / (rabi**4)
          !compute interatomic distances at rxyz
          rxd = rxyz(1,a)-rxb
          ryd = rxyz(2,a)-ryb
          rzd = rxyz(3,a)-rzb
          rabC = rxd**2+ryd**2+rzd**2
          rabC = sqrt(rabC)
          !compute function value
          rabiMrabC = rabi - rabC
          val = val + rabiMrabC**2 * rabim4 
          !compute negative gradient
          rabiMrabC = -2.0_gp*rabiMrabC * rabim4 / rabC
          ttx = rabiMrabC * rxd
          tty = rabiMrabC * ryd
          ttz = rabiMrabC * rzd
          force(1,b) = force(1,b) + ttx
          force(2,b) = force(2,b) + tty
          force(3,b) = force(3,b) + ttz
          force(1,a) = force(1,a) - ttx
          force(2,a) = force(2,a) - tty
          force(3,a) = force(3,a) - ttz
       enddo
    enddo

    !second sum
    do a = 1, nat
       waxi = oml*rxyzR(1,a) + lambda*rxyzP(1,a)
       wayi = oml*rxyzR(2,a) + lambda*rxyzP(2,a)
       wazi = oml*rxyzR(3,a) + lambda*rxyzP(3,a)
       ttx = waxi-rxyz(1,a)
       tty = wayi-rxyz(2,a)
       ttz = wazi-rxyz(3,a)
       tt   =  ttx**2 + tty**2 + ttz**2
       val = val + 1.e-6_gp * tt
       force(1,a) = force(1,a) + 2.e-6_gp * ttx
       force(2,a) = force(2,a) + 2.e-6_gp * tty
       force(3,a) = force(3,a) + 2.e-6_gp * ttz
    enddo
  end subroutine lst_penalty

end module lst_penalty_wrapper
module module_interpol
    implicit none

    private

    public :: lstpthpnt


contains
!=====================================================================
subroutine lstpthpnt(nat,rxyzR,rxyzP,lambda,rxyz)
!returns rxyz which is a point on the
!linear synchronous transit path
!corresponding to the interpolation parameter
!lambda
    use module_base
    use module_global_variables, only: fmax_tol => lst_fmax_tol
    use lst_penalty_wrapper
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp), intent(in) :: rxyzR(3,nat) !reactant
    real(gp), intent(in) :: rxyzP(3,nat) !product
    real(gp), intent(in) :: lambda       !interpol. parameter
    real(gp), intent(out) :: rxyz(3,nat)  !the positon on the
                                          !linear synchr. transit
                                          !path that corresponds
                                          !to lambda
    !internal
    real(gp) :: oml
!    real(gp), parameter :: fmax_tol=1.e-4_gp
!    real(gp), parameter :: fmax_tol=5.e-3_gp
    real(gp) :: fxyz(3,nat), epot

!<-DEBUG START------------------------------------------------------>
!character(len=5) :: fc5
!character(len=200) :: filename,line
!integer :: iat,istat
!integer, save :: ic
!real(gp) :: dmy
!character(len=5):: xat(22)
!<-DEBUG END-------------------------------------------------------->

    oml = 1._gp-lambda

!<-DEBUG START------------------------------------------------------>
!ic=ic+1
!close(33)
!open(unit=33,file='input001/pos001.ascii')
!read(33,*)
!read(33,*)
!read(33,*)
!read(33,*)
!do iat=1,22
!    read(33,'(a)',iostat=istat)line
!    if(istat/=0)exit
!    read(line,*)dmy,dmy,dmy,xat(iat)
!enddo
!close(33)
!<-DEBUG END-------------------------------------------------------->


    !input guess
    rxyz = oml*rxyzR+lambda*rxyzP
!<-DEBUG START------------------------------------------------------>
!if(mod(ic,10)==0)then
!write(fc5,'(i5.5)')ic
!write(filename,*)'posn_'//fc5//'.ascii'
!open(99,file=trim(adjustl((filename))))
!write(99,*)
!write(99,*)10.0 ,0, 10.0 
!write(99,*)0, 0, 10.0 
!do iat=1,nat
!write(99,'(3(1xes24.17),a)')rxyz(1,iat),rxyz(2,iat),&
!          rxyz(3,iat),xat(iat)
!enddo
!close(99)
!endif
!<-DEBUG END-------------------------------------------------------->
    call init_lst_wrapper(nat,rxyzR,rxyzP,lambda)
    call fire(nat,valforce,fmax_tol,rxyz,fxyz,epot)
    call finalize_lst_wrapper()
!<-DEBUG START------------------------------------------------------>
!if(mod(ic,10)==0)then
!write(fc5,'(i5.5)')ic
!write(filename,*)'posl_'//fc5//'.ascii'
!open(99,file=trim(adjustl((filename))))
!write(99,*)
!write(99,*)10.0 ,0, 10.0 
!write(99,*)0, 0, 10.0 
!do iat=1,nat
!write(99,'(3(1xes24.17),1x,a)')rxyz(1,iat)*0.529_gp,rxyz(2,iat)&
!          *0.529_gp,rxyz(3,iat)*0.529_gp,xat(iat)
!enddo
!close(99)
!endif
!<-DEBUG END-------------------------------------------------------->

!!$contains
!!$
!!$  !Fortran specifications forbids to pass internal routines as arguments
!!$  !otherwise, better to leave it included in module but require interface in fire
!!$  subroutine valforce(nat,rat,fat,epot)
!!$    !wrapper function for lst_penalty
!!$    use module_base
!!$    implicit none
!!$    !parameters
!!$    integer, intent(in) :: nat
!!$    real(gp), intent(in) :: rat(3,nat)
!!$    real(gp), intent(out) :: fat(3,nat)
!!$    real(gp), intent(out) :: epot
!!$    !local variables
!!$    call lst_penalty(nat,rxyzR,rxyzP,rat,lambda,epot,fat)
!!$  end subroutine valforce
  end subroutine lstpthpnt


!=====================================================================
subroutine fire(nat,valforce,fmax_tol,rxyz,fxyz,epot)
    use module_base
    use module_global_variables, only: iproc,lst_dt_max
    use yaml_output
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp)  :: fmax_tol
    real(gp), intent(inout) :: rxyz(3,nat),fxyz(3,nat)
    real(gp)  :: epot,fmax
    logical :: success
    !external :: valforce
    interface
       subroutine valforce(nat,rat,fat,epot)
         use module_defs, only: gp
         implicit none
         !parameters
         integer, intent(in) :: nat
         real(gp), intent(in) :: rat(3,nat)
         real(gp), intent(out) :: fat(3,nat)
         real(gp), intent(out) :: epot
       end subroutine valforce
    end interface
    !internal
    integer :: maxit
    real(gp) :: count_fr,fnrm
    real(gp) :: vxyz(3,nat),ff(3,nat),power,dt_max
    real(gp) :: at1,at2,at3,vxyz_norm,fxyz_norm,alpha
    integer :: iat,iter,check,cut
    integer, parameter :: n_min=5
    real(gp), parameter :: f_inc=1.1_gp, f_dec=0.5_gp
    real(gp), parameter :: alpha_start=0.1_gp,f_alpha=0.99_gp
    real(gp) :: dt
    real(gp) :: ddot,dnrm2
    real(gp) :: epotold
!!<-DEBUG START------------------------------------------------------>
!character(len=5), allocatable :: xat(:)
!character(len=5) :: fc5
!character(len=50) :: filename
!!<-DEBUG END-------------------------------------------------------->

    count_fr=0._gp
   ! dt_max = 0.5_gp !for lj
   ! dt_max = 1.8_gp
   ! dt_max = 5.0_gp !for bigdft
    dt_max=lst_dt_max
    maxit=15000

    success=.false.
    dt=dt_max
    alpha=alpha_start
    vxyz=0.0_gp
    check=0
    cut=1
    call valforce(nat,rxyz,ff,epot)
    epotold=epot
    count_fr=count_fr+1.0_gp
    do iter=1,maxit
!!!! 1000   continue !!!to be avoided 
        call daxpy(3*nat,dt,vxyz(1,1),1,rxyz(1,1),1)
        call daxpy(3*nat,0.5_gp*dt*dt,ff(1,1),1,rxyz(1,1),1)
!!<-DEBUG START------------------------------------------------------>
!write(fc5,'(i5.5)')iter
!write(filename,*)'pos_'//fc5//'.ascii'
!open(99,file=trim(adjustl((filename))))
!write(99,*)
!write(99,*)10.0 ,0, 10.0 
!write(99,*)0, 0, 10.0 
!do iat=1,nat
!write(99,'(3(1xes24.17),a)')rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),&
!xat(iat)
!enddo
!close(99)
!!<-DEBUG END-------------------------------------------------------->


        call valforce(nat,rxyz,fxyz,epot)
        count_fr=count_fr+1.0_gp
        do iat=1,nat
           at1=fxyz(1,iat)
           at2=fxyz(2,iat)
           at3=fxyz(3,iat)
           !Evolution of the velocities of the system
           vxyz(1,iat)=vxyz(1,iat) + (.5_gp*dt) * (at1 + ff(1,iat))
           vxyz(2,iat)=vxyz(2,iat) + (.5_gp*dt) * (at2 + ff(2,iat))
           vxyz(3,iat)=vxyz(3,iat) + (.5_gp*dt) * (at3 + ff(3,iat))
           !Memorization of old forces
           ff(1,iat) = at1
           ff(2,iat) = at2
           ff(3,iat) = at3
        end do
        call fnrmandforcemax(fxyz,fnrm,fmax,nat)
        fnrm=sqrt(fnrm) 
        call convcheck(fmax,fmax_tol,check)
        if(check > 5)then
!!<-DEBUG START------------------------------------------------------>
!            write(*,'(a,x,i0,5(1x,es14.7))')&
!            'FIRE converged # e evals, epot, &
!            fmax, fnrm, dt, alpha: ',&
!            int(count_fr),epot,fmax,sqrt(fnrm),dt,alpha
!!<-DEBUG END-------------------------------------------------------->
            success=.true.
            return
!!<-DEBUG START------------------------------------------------------>
!        else
!            write(*,'(a,x,i0,5(1x,es14.7))')&
!            'FIRE # e evals, epot, &
!            fmax, fnrm, dt, alpha: ',&
!            int(count_fr),epot,fmax,sqrt(fnrm),dt,alpha
!!<-DEBUG END-------------------------------------------------------->
        endif
        power = ddot(3*nat,fxyz,1,vxyz,1)
        vxyz_norm = dnrm2(3*nat,vxyz,1)
        fxyz_norm = dnrm2(3*nat,fxyz,1)
        vxyz = (1.0_gp-alpha)*vxyz + &
               alpha * fxyz * vxyz_norm / fxyz_norm
        if(power<=0 .or. epot>epotold)then
            vxyz=0.0_gp
            cut=iter
            dt=dt*f_dec
            alpha=alpha_start
        else if(power > 0 .and. iter-cut >n_min)then
            dt= min(dt*f_inc,dt_max)
            alpha = alpha*f_alpha
        endif
        epotold=epot
    enddo
    if(fmax > fmax_tol .and. iproc==0)then
        call yaml_warning('(MHGPS) Minimization of Linear '//&
                          'Synchronous Transit Path not converged:')
        call yaml_map('(MHGPS) epot, fmax, fnrm, dt, alpha',&
             (/epot,fmax,fnrm,dt,alpha/))
    endif
end subroutine
!=====================================================================
subroutine convcheck(fmax,forcemax,check)
  use module_base
  implicit none
  real(gp), intent(in):: fmax, forcemax
  integer, intent(inout)::check

  if ( fmax < forcemax) then
    check=check+1
  else
    check=0
  endif

end subroutine convcheck
!=====================================================================
subroutine fnrmandforcemax(ff,fnrm,fmax,nat)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), intent(in):: ff(3,nat)
  real(gp), intent(out):: fnrm, fmax
  real(gp):: t1,t2,t3
  integer:: iat
  real(gp) :: ddot

  fmax=0._gp
  do iat=1,nat
     t1=ff(1,iat)**2
     t2=ff(2,iat)**2
     t3=ff(3,iat)**2
     fmax=max(fmax,sqrt(t1+t2+t3))
  enddo

  !this is the norm of the forces of non-blocked atoms
  !one has to discuss whether also the center
  !mass shift should be added
  fnrm=ddot(3*nat,ff(1,1),1,ff(1,1),1)
end subroutine fnrmandforcemax
!=====================================================================

end module module_interpol

