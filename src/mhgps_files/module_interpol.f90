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

module module_interpol
    implicit none

    private

    public :: lstpthpnt


contains
!=====================================================================
subroutine lstpthpnt(runObj,mhgpsst,uinp,rxyzR,rxyzP,lambda,rxyz)
!returns rxyz which is a point on the
!linear synchronous transit path
!corresponding to the interpolation parameter
!lambda
    use module_base
    use bigdft_run
    use module_forces
    use module_mhgps_state
    use module_userinput
    implicit none
    !parameters
    type(run_objects), intent(in) :: runObj
    type(mhgps_state), intent(in) :: mhgpsst
    type(userinput), intent(in) :: uinp
    real(gp), intent(in) :: rxyzR(3,runObj%atoms%astruct%nat) !reactant
    real(gp), intent(in) :: rxyzP(3,runObj%atoms%astruct%nat) !product
    real(gp), intent(in) :: lambda       !interpol. parameter
    real(gp), intent(out) :: rxyz(3,runObj%atoms%astruct%nat)  !the positon on the
                                          !linear synchr. transit
                                          !path that corresponds
                                          !to lambda
    !internal
    real(gp) :: oml
    real(gp) :: fxyz(3,runObj%atoms%astruct%nat), epot
    real(gp) :: drxyz(3,runObj%atoms%astruct%nat)

    drxyz = rxyzP-rxyzR
    !clean forces is just a precaution. In principle
    !the coressponding frozen coordinates should already be
    !identical in rxyzP and rxyzR
    call clean_forces_base(runObj%atoms,drxyz)

    rxyz = rxyzR + lambda * drxyz
    call fire(mhgpsst,uinp,runObj,rxyzR,rxyzP,lambda,rxyz,fxyz,epot)
  end subroutine lstpthpnt


!=====================================================================
subroutine lst_penalty(runObj,rxyzR,rxyzP,rxyz,lambda,val,force)
    !computes the linear synchronous penalty function
    !and gradient.
    !OpenMP parallelized.
    !Attention: Cells of rxyzR and rxyzP and rxyz must be identical
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
    !HISTORY:
    !June 2014   : finished non-periodic version
    !January 2015: OpenMP parallelized
    !March 2015  : generalized to full periodicity
    use module_base
    use bigdft_run
    use module_forces
    implicit none
    !parameters
    type(run_objects), intent(in) :: runObj
    real(gp), intent(in) :: rxyzR(3,runObj%atoms%astruct%nat) !reactant
    real(gp), intent(in) :: rxyzP(3,runObj%atoms%astruct%nat) !product
    real(gp), intent(in) :: rxyz(3,runObj%atoms%astruct%nat)  !the positon at which
    !the penalty fct. is to
    !be evaluated
    real(gp), intent(in)  :: lambda !interpolation parameter
    !(0<=lambda<=1)
    real(gp), intent(out) :: val  !the value of the penalty function
    real(gp), intent(out) :: force(3,runObj%atoms%astruct%nat) !the negative gradient of
                                                               !the penal. fct. at
                                                               !rxyz
    !internal
    integer, parameter :: necmin=1
    integer :: b !outer loop
    integer :: a !inner loop
    integer :: i,j,k !loops over periodic images
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
    real(gp) :: rxyzRint(3,runObj%atoms%astruct%nat) !Reactant in internal coord.
    real(gp) :: rxyzPint(3,runObj%atoms%astruct%nat) !Product in internal coord.
    real(gp) :: rxyzint(3,runObj%atoms%astruct%nat) !the positon at which the 
                                                 !penalty fct. is to be
                                                 !evaluated in internal coord.
    real(gp) :: latvec(3,3), cell(3), weight, trans(3)
    integer :: nec1, nec2, nec3
    real(gp) :: rcut

    !BigDFT only supports orthorombic unit cells (as of Mar. 2015).
    !Therefore, the corresponding cell(:) data structure can only hold three numbers.
    !We therefore build here the full lattice vectors.
    !If the BigDFT data structure is changed one day, or if the code
    !is used outside of the bigdft suite, the lattice vectors should
    !be passed and not be constructed here.
    !In principle, this subroutine should also support
    !non-orthorombic unit cells, if latvec is passed appropriately.
    cell = bigdft_get_cell(runObj)
    rcut=huge(1.0_gp)
    if(bigdft_get_geocode(runObj)=='F')then
        latvec(1,1)=0.0_gp
        latvec(2,1)=0.0_gp
        latvec(3,1)=0.0_gp
        latvec(1,2)=0.0_gp
        latvec(2,2)=0.0_gp
        latvec(3,2)=0.0_gp
        latvec(1,3)=0.0_gp
        latvec(2,3)=0.0_gp
        latvec(3,3)=0.0_gp
        nec1=1
        nec2=1
        nec3=1
    else if(bigdft_get_geocode(runObj)=='P')then
        latvec(1,1)=cell(1)
        latvec(2,1)=0.0_gp
        latvec(3,1)=0.0_gp
        latvec(1,2)=0.0_gp
        latvec(2,2)=cell(2)
        latvec(3,2)=0.0_gp
        latvec(1,3)=0.0_gp
        latvec(2,3)=0.0_gp
        latvec(3,3)=cell(3)
        call give_rcut(latvec,necmin,rcut)
        call n_rep_dim(latvec,rcut,nec1,nec2,nec3)
    else if(bigdft_get_geocode(runObj)=='S')then
        !in BigDFT, the non-periodic direction is the y-direction
        !=> cell(2)=0.0_gp on input
        latvec(1,1)=cell(1)
        latvec(2,1)=0.0_gp
        latvec(3,1)=0.0_gp
        latvec(1,2)=0.0_gp
        latvec(2,2)=1.0_gp
        latvec(3,2)=0.0_gp
        latvec(1,3)=0.0_gp
        latvec(2,3)=0.0_gp
        latvec(3,3)=cell(3)
        call give_rcut(latvec,necmin,rcut)
        call n_rep_dim(latvec,rcut,nec1,nec2,nec3)
        nec2=1
    else
        call f_err_throw('lst_penalty: specified boundary conditions'//&
             ' not supported.')
    endif
!write(199,*)rcut
!write(199,*)nec1,nec2,nec3
!!stop
!nec1=1
!nec2=1
!nec3=1



    oml=1.0_gp-lambda!one minus lambda
    val=0.0_gp!function value
    force=0.0_gp!negative gradient
    !$omp parallel default(private) shared(runObj, rxyzR, rxyzP, rxyz, &
    !$omp val, force, oml, lambda, latvec, nec1, nec2, nec3, rcut)
    !$omp do schedule(dynamic) reduction(+:force,val)
    !first sum
    do b = 1, runObj%atoms%astruct%nat-1
       rxRb = rxyzR(1,b)
       ryRb = rxyzR(2,b)
       rzRb = rxyzR(3,b)
       rxPb = rxyzP(1,b)
       ryPb = rxyzP(2,b)
       rzPb = rxyzP(3,b)
       rxb = rxyz(1,b)
       ryb = rxyz(2,b)
       rzb = rxyz(3,b)
       do i=0,nec1-1
        do j=0,nec2-1
          do k=0,nec3-1
           do a = b+1, runObj%atoms%astruct%nat
              trans(1) = latvec(1,1)*real(i,gp)+latvec(1,2)*real(j,gp)+latvec(1,3)*real(k,gp)
              trans(2) = latvec(2,1)*real(i,gp)+latvec(2,2)*real(j,gp)+latvec(2,3)*real(k,gp)
              trans(3) = latvec(3,1)*real(i,gp)+latvec(3,2)*real(j,gp)+latvec(3,3)*real(k,gp)
              !compute interatomic distances of reactant
              rabR = (rxyzR(1,a)+trans(1)-rxRb)**2+&
                     (rxyzR(2,a)+trans(2)-ryRb)**2+&
                     (rxyzR(3,a)+trans(3)-rzRb)**2
              rabR = sqrt(rabR)
              !compute interatomic distances of product
              rabP = (rxyzP(1,a)+trans(1)-rxPb)**2+&
                     (rxyzP(2,a)+trans(2)-ryPb)**2+&
                     (rxyzP(3,a)+trans(3)-rzPb)**2
              rabP = sqrt(rabP)
              !compute interatomic distances at rxyz
              rxd = rxyz(1,a)+trans(1)-rxb
              ryd = rxyz(2,a)+trans(2)-ryb
              rzd = rxyz(3,a)+trans(3)-rzb
              rabC = rxd**2+ryd**2+rzd**2
              rabC = sqrt(rabC)
              !compute interpolated interatomic distances
              rabi = oml*rabR+lambda*rabP
              rabim4 = 1.0_gp / (rabi**4)
!              if(rabC>rcut)cycle
              if(rabi>rcut)cycle
              if(i==0 .and. j==0 .and. k==0)then
                  weight= rabim4
              else
                  weight = rabim4 - 1.0_gp/ (rcut**4)
              endif
              !compute function value
              rabiMrabC = rabi - rabC
              val = val + rabiMrabC**2 * weight 
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
         enddo
        enddo
    enddo
    !$omp end do

    !$omp do schedule(dynamic) reduction(+:force,val)
    !second sum
    do a = 1, runObj%atoms%astruct%nat
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
    !$omp end do
    !$omp end parallel

    call clean_forces_base(runObj%atoms,force)
end subroutine lst_penalty
!=====================================================================
subroutine fire(mhgpsst,uinp,runObj,rxyzR,rxyzP,lambda,rxyz,fxyz,epot)
    use module_base
    use bigdft_run
    use module_userinput
    use yaml_output
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(in) :: mhgpsst
    type(userinput), intent(in) :: uinp
    type(run_objects), intent(in) :: runObj
    real(gp), intent(in) :: rxyzR(3,runObj%atoms%astruct%nat) !reactant
    real(gp), intent(in) :: rxyzP(3,runObj%atoms%astruct%nat) !product
    real(gp), intent(in) :: lambda !interpolation parameter
    real(gp), intent(inout) :: rxyz(3,runObj%atoms%astruct%nat),fxyz(3,runObj%atoms%astruct%nat)
    real(gp)  :: epot,fmax
    logical :: success
    !external :: valforce
!    interface
!       subroutine valforce(runObj%atoms%astruct%nat,rat,fat,epot)
!         use module_defs, only: gp
!         implicit none
!         !parameters
!         integer, intent(in) :: runObj%atoms%astruct%nat
!         real(gp), intent(in) :: rat(3,runObj%atoms%astruct%nat)
!         real(gp), intent(out) :: fat(3,runObj%atoms%astruct%nat)
!         real(gp), intent(out) :: epot
!       end subroutine valforce
!    end interface
    !internal
    integer :: maxit
    real(gp) :: count_fr,fnrm
    real(gp) :: vxyz(3,runObj%atoms%astruct%nat),ff(3,runObj%atoms%astruct%nat),power,dt_max
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
    dt_max=uinp%lst_dt_max
    maxit=15000

    success=.false.
    dt=dt_max
    alpha=alpha_start
    vxyz=0.0_gp
    check=0
    cut=1
!    call valforce(nat,rxyz,ff,epot)
    call lst_penalty(runObj,rxyzR,rxyzP,rxyz,lambda,epot,ff)
    epotold=epot
    count_fr=count_fr+1.0_gp
    do iter=1,maxit
!!!! 1000   continue !!!to be avoided 
        call daxpy(3*runObj%atoms%astruct%nat,dt,vxyz(1,1),1,rxyz(1,1),1)
        call daxpy(3*runObj%atoms%astruct%nat,0.5_gp*dt*dt,ff(1,1),1,rxyz(1,1),1)
!!<-DEBUG START------------------------------------------------------>
!write(fc5,'(i5.5)')iter
!write(filename,*)'pos_'//fc5//'.ascii'
!open(99,file=trim(adjustl((filename))))
!write(99,*)
!write(99,*)10.0 ,0, 10.0 
!write(99,*)0, 0, 10.0 
!do iat=1,runObj%atoms%astruct%nat
!write(99,'(3(1xes24.17),a)')rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),&
!'  Pt'!xat(iat)
!enddo
!close(99)
!!<-DEBUG END-------------------------------------------------------->


!        call valforce(nat,rxyz,fxyz,epot)
        call lst_penalty(runObj,rxyzR,rxyzP,rxyz,lambda,epot,fxyz)
        count_fr=count_fr+1.0_gp
        do iat=1,runObj%atoms%astruct%nat
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
        call fnrmandforcemax(fxyz,fnrm,fmax,runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        call convcheck(fmax,uinp%lst_fmax_tol,check)
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
        power = ddot(3*runObj%atoms%astruct%nat,fxyz,1,vxyz,1)
        vxyz_norm = dnrm2(3*runObj%atoms%astruct%nat,vxyz,1)
        fxyz_norm = dnrm2(3*runObj%atoms%astruct%nat,fxyz,1)
        if(abs(fxyz_norm)>1.e-14_gp)then
            vxyz = (1.0_gp-alpha)*vxyz + &
                   alpha * fxyz * vxyz_norm / fxyz_norm
        else
            vxyz = (1.0_gp-alpha)*vxyz
        endif
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
    if(fmax > uinp%lst_fmax_tol .and. mhgpsst%iproc==0)then
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
subroutine dist2plane(point,nvec,ppoint,dist)
!This subroutine will calculate  the distance between a plane and a point in
!space. The point is 'point', the normalized normal vector of the plane is
!'nvec', 'ppoint' is an arbitrary point on the plane !and the output is the
!distance 'dist' 
    use module_base
    implicit none
    !parameters 
    real(gp), intent(in) :: point(3),nvec(3),ppoint(3)
    real(gp), intent(out):: dist
    !local
    real(gp) :: p
    !functions
    real(gp) :: ddot
    
    p=ddot(3,nvec(1),1,ppoint(1),1)
    dist=ddot(3,nvec(1),1,point(1),1)-p
end subroutine
!=====================================================================
subroutine nveclatvec(latvec,nvec)
!Calculates the normalized normal vectors to the 3 planes of the cell
    use module_base
    implicit none
    !parameters
    real(gp), intent(in) :: latvec(3,3)
    real(gp), intent(out):: nvec(3,3)
    !internal
    real(gp) :: a(3),b(3),crossp(3),norm
    integer  :: i
    !functions
    real(gp) :: dnrm2
    do i=1,3
        a=latvec(:,i)
        b=latvec(:,mod(i,3)+1)
        crossp(1)=a(2)*b(3)-a(3)*b(2)
        crossp(2)=a(3)*b(1)-a(1)*b(3)
        crossp(3)=a(1)*b(2)-a(2)*b(1)
        norm=dnrm2(3,crossp(1),1)
        nvec(:,i)=crossp(:)/norm
    enddo
end subroutine
!=====================================================================
subroutine give_rcut(latvec,necmin,rcut)
    use module_base
    implicit none
    !parameters
    real(gp), intent(in) :: latvec(3,3)
    integer, intent(in) :: necmin
    real(gp), intent(out) :: rcut
    !internal
    integer :: i
    real(gp) :: nvec(3,3), point0(3),dist(3)
    !calculate normal vectors
    call nveclatvec(latvec,nvec)
    point0=(/0.0_gp,0.0_gp,0.0_gp/)
    do i=1,3
        call dist2plane(latvec(1,mod(i+1,3)+1),nvec(1,i),point0(1),dist(i))
    enddo
    dist=abs(dist)
    rcut=real(necmin-1,gp)*maxval(dist)
end subroutine give_rcut
!=====================================================================
subroutine n_rep_dim(latvec,cut,nec1,nec2,nec3)
!This subroutine will return how many periodic expansions for each lattice
!vector direction are necessary for the periodic boundary conditions
!for the given cut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
    use module_base
    implicit none
    !parameters
    real(gp), intent(in) :: latvec(3,3)
    real(gp), intent(in) :: cut
    integer, intent(out)  :: nec1,nec2,nec3
    !internal
    real(gp) :: nvec(3,3)
    real(gp) :: point0(3)
    real(gp) :: dist(3)
    integer:: i
    nec1=0
    nec2=0
    nec3=0
    !calculate normal vectors
    call nveclatvec(latvec,nvec)
    point0=(/0.0_gp,0.0_gp,0.0_gp/)
    do i=1,3
        call dist2plane(latvec(1,mod(i+1,3)+1),nvec(1,i),point0(1),dist(i))
    enddo
    dist=abs(dist)
    nec1=int(cut/dist(2))+1
    nec2=int(cut/dist(3))+1
    nec3=int(cut/dist(1))+1
end subroutine
!=====================================================================
!subroutine rxyz_cart2int(nat,latvec,rxyzint,rxyzcart)
!    use module_base
!!This subrouine will convert cartesian to internal coordinates
!    implicit none
!    !parameters
!    integer, intent(in) :: nat
!    real(gp), intent(in) :: latvec(3,3)
!    real(gp), intent(in) :: rxyzint(3,nat)
!    real(gp), intent(out) :: rxyzcart(3,nat)
!    !internal
!    real(gp) :: latvecinv(3,3)
!    real(gp) :: div
!    integer:: iat
!
!    !compute latvecinv, the inverse of matrix latvec
!    div=(latvec(1,1)*latvec(2,2)*latvec(3,3)-latvec(1,1)*latvec(2,3)*latvec(3,2)-&
!        latvec(1,2)*latvec(2,1)*latvec(3,3)+latvec(1,2)*latvec(2,3)*latvec(3,1)+&
!        latvec(1,3)*latvec(2,1)*latvec(3,2)-latvec(1,3)*latvec(2,2)*latvec(3,1))
!    div=1.0_gp/div
!    latvecinv(1,1) = (latvec(2,2)*latvec(3,3)-latvec(2,3)*latvec(3,2))*div
!    latvecinv(1,2) =-(latvec(1,2)*latvec(3,3)-latvec(1,3)*latvec(3,2))*div
!    latvecinv(1,3) = (latvec(1,2)*latvec(2,3)-latvec(1,3)*latvec(2,2))*div
!    latvecinv(2,1) =-(latvec(2,1)*latvec(3,3)-latvec(2,3)*latvec(3,1))*div
!    latvecinv(2,2) = (latvec(1,1)*latvec(3,3)-latvec(1,3)*latvec(3,1))*div
!    latvecinv(2,3) =-(latvec(1,1)*latvec(2,3)-latvec(1,3)*latvec(2,1))*div
!    latvecinv(3,1) = (latvec(2,1)*latvec(3,2)-latvec(2,2)*latvec(3,1))*div
!    latvecinv(3,2) =-(latvec(1,1)*latvec(3,2)-latvec(1,2)*latvec(3,1))*div
!    latvecinv(3,3) = (latvec(1,1)*latvec(2,2)-latvec(1,2)*latvec(2,1))*div
!
!    do iat=1,nat
!        call dgemv('n',3,3,1.0_gp,latvecinv(1,1),3,rxyzcart(1,iat),1,&
!             0.0_gp, rxyzint(1,iat),1)
!!        rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
!    enddo
!end subroutine rxyz_cart2int


end module module_interpol

