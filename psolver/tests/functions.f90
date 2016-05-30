!> @file
!!  Provide the test functions for the different examples of Poisson Solver
!! @author
!!    Copyright (C) 2002-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine PS_Check_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse), intent(inout) :: parser

  call yaml_cl_parse_option(parser,'ndim','None',&
       'Domain Sizes','n',&
       dict_new('Usage' .is. &
       'Sizes of the simulation domain of the check',&
       'Allowed values' .is. &
       'Yaml list of integers. If a scalar integer is given, all the dimensions will have this size.'),first_option=.true.)

  call yaml_cl_parse_option(parser,'geocode','F',&
       'Boundary conditions','g',&
       dict_new('Usage' .is. &
       'Set the boundary conditions of the run',&
       'Allowed values' .is. &
       'String scalar. "F","S","W","P" boundary conditions are allowed'))

  call yaml_cl_parse_option(parser,'method','None',&
       'Embedding method','m',&
       dict_new('Usage' .is. &
       'Set the embedding method used. A non present value implies vacuum treatment.',&
       'Allowed values' .is. &
       dict_new("PI" .is. 'Polarization iteration Method',&
       "PCG" .is. 'Preconditioned Conjugate Gradient')))

  call yaml_cl_parse_option(parser,'seteps','4',&
       'Epsilon determination method','e',&
       dict_new('Usage' .is. &
       'Set the dielectric constant determination method.',&
       'Allowed values' .is. &
       dict_new('1' .is. 'Analytical epsilon' ,&
       '2' .is. 'analytical electron dependence',&
       '3' .is. 'real electron density from cube file (need electroninc_density.cube)',&
       '4' .is. 'calculate the cavity and dump it on disk',&
       '5' .is. 'Solves GPe with PCG customized (should be identical to 4 + PCG)',&
       '6' .is. 'Modified Poisson Botzmann Equation solver',&
       '7' .is. 'Solves GPe with PCG customized and a restart is implemented',&
       '8' .is. 'Solves GPe with PI customized (should be identical to 4 + PI)',&
       '9' .is. 'Solves GPe with PSD',&
       '10' .is. 'Solves GPe with PCG customized and is coupled with PSD')))

  call yaml_cl_parse_option(parser,'accel','No',&
       'GPU Acceleration','a',&
       dict_new('Usage' .is. &
       'Boolean, set the GPU acceleration'))

  call yaml_cl_parse_option(parser,'logfile','Yes',&
       'Write logfile','l',&
       dict_new('Usage' .is. &
       'Boolean, set the logfile as log.yaml'))

  call yaml_cl_parse_option(parser,'deltacav','None',&
       'Delta rigid cavity','d',&
       dict_new('Usage' .is. &
       'Sizes of the delta for error function',&
       'Allowed values' .is. &
       'Real value'))

  call yaml_cl_parse_option(parser,'input','None',&
       'Inputfile of Poisson solver','i',&
       dict_new('Usage' .is. &
       'Put the dictionary of PS inputfile to preload some parameters',&
       'Allowed values' .is. &
       'yaml Dictionary'))


end subroutine PS_Check_options

!> Error function in double precision
subroutine derf_local(derf_yy,yy)
  use PSbase, only: dp
  implicit none
  real(dp),intent(in) :: yy
  real(dp),intent(out) :: derf_yy
  integer          ::  done,ii,isw
  real(dp), parameter :: &
                                ! coefficients for 0.0 <= yy < .477
       &  pp(5)=(/ 113.8641541510502e0_dp, 377.4852376853020e0_dp,  &
       &           3209.377589138469e0_dp, .1857777061846032e0_dp,  &
       &           3.161123743870566e0_dp /)
  real(dp), parameter :: &
       &  qq(4)=(/ 244.0246379344442e0_dp, 1282.616526077372e0_dp,  &
       &           2844.236833439171e0_dp, 23.60129095234412e0_dp/)
  ! coefficients for .477 <= yy <= 4.0
  real(dp), parameter :: &
       &  p1(9)=(/ 8.883149794388376e0_dp, 66.11919063714163e0_dp,  &
       &           298.6351381974001e0_dp, 881.9522212417691e0_dp,  &
       &           1712.047612634071e0_dp, 2051.078377826071e0_dp,  &
       &           1230.339354797997e0_dp, 2.153115354744038e-8_dp, &
       &           .5641884969886701e0_dp /)
  real(dp), parameter :: &
       &  q1(8)=(/ 117.6939508913125e0_dp, 537.1811018620099e0_dp,  &
       &           1621.389574566690e0_dp, 3290.799235733460e0_dp,  &
       &           4362.619090143247e0_dp, 3439.367674143722e0_dp,  &
       &           1230.339354803749e0_dp, 15.74492611070983e0_dp/)
  ! coefficients for 4.0 < y,
  real(dp), parameter :: &
       &  p2(6)=(/ -3.603448999498044e-01_dp, -1.257817261112292e-01_dp,   &
       &           -1.608378514874228e-02_dp, -6.587491615298378e-04_dp,   &
       &           -1.631538713730210e-02_dp, -3.053266349612323e-01_dp/)
  real(dp), parameter :: &
       &  q2(5)=(/ 1.872952849923460e0_dp   , 5.279051029514284e-01_dp,    &
       &           6.051834131244132e-02_dp , 2.335204976268692e-03_dp,    &
       &           2.568520192289822e0_dp /)
  real(dp), parameter :: &
       &  sqrpi=.5641895835477563e0_dp, xbig=13.3e0_dp, xlarge=6.375e0_dp, xmin=1.0e-10_dp
  real(dp) ::  res,xden,xi,xnum,xsq,xx

  xx = yy
  isw = 1
  !Here change the sign of xx, and keep track of it thanks to isw
  if (xx<0.0e0_dp) then
     isw = -1
     xx = -xx
  end if

  done=0

  !Residual value, if yy < -6.375e0_dp
  res=-1.0e0_dp

  !abs(yy) < .477, evaluate approximation for erfc
  if (xx<0.477e0_dp) then
     ! xmin is a very small number
     if (xx<xmin) then
        res = xx*pp(3)/qq(3)
     else
        xsq = xx*xx
        xnum = pp(4)*xsq+pp(5)
        xden = xsq+qq(4)
        do ii = 1,3
           xnum = xnum*xsq+pp(ii)
           xden = xden*xsq+qq(ii)
        end do
        res = xx*xnum/xden
     end if
     if (isw==-1) res = -res
     done=1
  end if

  !.477 < abs(yy) < 4.0 , evaluate approximation for erfc
  if (xx<=4.0e0_dp .and. done==0 ) then
     xsq = xx*xx
     xnum = p1(8)*xx+p1(9)
     xden = xx+q1(8)
     do ii=1,7
        xnum = xnum*xx+p1(ii)
        xden = xden*xx+q1(ii)
     end do
     res = xnum/xden
     res = res* exp(-xsq)
     if (isw.eq.-1) then
        res = res-1.0e0_dp
     else
        res=1.0e0_dp-res
     end if
     done=1
  end if

  !y > 13.3e0_dp
  if (isw > 0 .and. xx > xbig .and. done==0 ) then
     res = 1.0e0_dp
     done=1
  end if

  !4.0 < yy < 13.3e0_dp  .or. -6.375e0_dp < yy < -4.0
  !evaluate minimax approximation for erfc
  if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
     xsq = xx*xx
     xi = 1.0e0_dp/xsq
     xnum= p2(5)*xi+p2(6)
     xden = xi+q2(5)
     do ii = 1,4
        xnum = xnum*xi+p2(ii)
        xden = xden*xi+q2(ii)
     end do
     res = (sqrpi+xi*xnum/xden)/xx
     res = res* exp(-xsq)
     if (isw.eq.-1) then
        res = res-1.0e0_dp
     else
        res=1.0e0_dp-res
     end if
  end if

  !All cases have been investigated
  derf_yy = res

end subroutine derf_local

!> This subroutine builds some analytic functions that can be used for 
!! testing the poisson solver.
!! The default choice is already well-tuned for comparison.
!! WARNING: not all the test functions can be used for all the boundary conditions of
!! the poisson solver, in order to have a reliable analytic comparison.
!! The parameters of the functions must be adjusted in order to have a sufficiently localized
!! function in the isolated direction and an explicitly periodic function in the periodic ones.
!! Beware of the high-frequency components that may false the results when hgrid is too high.
recursive subroutine test_functions_box(mesh,nspden,a_gauss,&
     density,potential,rhopot,pot_ion,offset)
  use box
  use f_utils
  use f_precisions
  use PSbase
  use f_jmp
  use time_profiling
  use numerics
  use dynamic_memory
  implicit none
  type(cell), intent(in) :: mesh !<definition of the cell
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: a_gauss
  real(kind=8), intent(out) :: offset
  real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(out) :: pot_ion,potential
  real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3),nspden), intent(out) :: density,rhopot
  !local variables
  integer, parameter :: nrep=10
  logical ::  separable=.true.,old=.true. !save attribute for the profiling
  integer :: i1,i2,i3,ifx,ify,ifz,i,irep,n01,n02,n03
  
  real(kind=8) :: x1,x2,x3,length,denval,a2,derf_tt,factor,r,r2,hx,hy,hz
  real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,acell
  integer(f_long) :: t1,t0
  type(box_iterator) :: bit
  type(f_jmpbuf), save :: jmpbuf

  !backward compatibility
  acell=mesh%ndims(1)*mesh%hgrids(1)
  hx=mesh%hgrids(1)
  hy=mesh%hgrids(2)
  hz=mesh%hgrids(3)
  n01=mesh%ndims(1)
  n02=mesh%ndims(2)
  n03=mesh%ndims(3)

  !select the specifications for the loop
  select case (cell_geocode(mesh))
  case('P')
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=5
     ify=5
     ifz=5
     !parameters of the test functions
     ax=length
     ay=length
     az=length
     bx=2.d0!real(nu,kind=8)
     by=2.d0!real(nu,kind=8)
     bz=2.d0
     !other factors
     factor = 1.0_dp
  case('S')
     !parameters for the test functions
     length=acell
     a=0.5d0/a_gauss**2
     !test functions in the three directions
     ifx=5
     ifz=5
     !non-periodic dimension
     ify=6
     !parameters of the test functions
     ax=length
     az=length
     bx=2.d0!real(nu,kind=8)
     bz=2.d0!real(nu,kind=8)
     !non-periodic dimension
     ay=length
     by=a
     factor = oneofourpi
  case('F')
     a2 = a_gauss**2
     !Normalization
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     separable=.false.
  case default
     print *,'geometry code not admitted',cell_geocode(mesh)
     stop
  end select


  t0=f_time()
  !     do irep=1,nrep 
  entry test_new()
  if (separable .and. .not. old) then
     denval=0.d0 !value for keeping the density positive
     bit=box_iter(mesh,centered=.true.)
     !here separability is restored
     do while(box_next_z(bit))
        call functions(bit%rxyz(3),az,bz,fz,fz2,ifz)
        do while(box_next_y(bit))
           call functions(bit%rxyz(2),ay,by,fy,fy2,ify)
           do while(box_next_x(bit))
              call functions(bit%rxyz(1),ax,bx,fx,fx2,ifx)
              do i=1,nspden
                 density(bit%i,bit%j,bit%k,i) = factor/real(nspden,kind=8)*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
              end do
              potential(bit%i,bit%j,bit%k) = -factor*fourpi*fx*fy*fz
              denval=max(denval,-density(bit%i,bit%j,bit%k,1))
           end do
        end do
     end do
  end if
  call f_profile_end(test_new,jmpbuf)

  !     end do
  t1=f_time()
!!$  print *,f_humantime(t1-t0)
!!$  print *,'there2',denval
  t0=f_time()
  !     do irep=1,nrep 
  entry test_old()

!print *,'here',separable,old,nrep,associated(jmpbuf%jmp_buf)
  if (separable .and. old) then
     denval=0.d0 !value for keeping the density positive
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n01/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx2,ifx)
              print *,'i1,i2,i3',x1,x2,x3
              do i=1,nspden
                 density(i1,i2,i3,i) = factor/real(nspden,kind=8)*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
              end do
              potential(i1,i2,i3) = -fourpi*factor*fx*fy*fz
              denval=max(denval,-density(i1,i2,i3,1))
           end do
        end do
     end do
  end if
  call f_profile_end(test_old,jmpbuf)
  !     end do
  t1=f_time()

!!$  print *,f_humantime(t1-t0)
!!$  print *,'there3',denval

  if (.not. separable .and. .not. old) then
     bit=box_iter(mesh,centered=.true.)
     do while(box_next_point(bit))
        r2=square(mesh,bit%rxyz)
        do i=1,nspden
           density(bit%i,bit%j,bit%k,i) = &
                1.d0/real(nspden,kind=8)*max(factor*safe_exp(-r2/a2),1d-24)
        end do
        r = sqrt(r2)
        !Potential from a gaussian
        if (r == 0.0_dp) then
           potential(bit%i,bit%j,bit%k) = 2.d0/(sqrt(pi)*a_gauss)
        else
           call derf_local(derf_tt,r/a_gauss)
           potential(bit%i,bit%j,bit%k) = derf_tt/r
        end if
     end do
  end if
  if (.not. separable .and. old) then
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              do i=1,nspden
                 density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*max(factor*exp(-r2/a2),1d-24)
              end do
              r = sqrt(r2)
              !Potential from a gaussian
              if (r == 0.d0) then
                 potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
              else
                 call derf_local(derf_tt,r/a_gauss)
                 potential(i1,i2,i3) = derf_tt/r
              end if
           end do
        end do
     end do
  end if
  !here we can profile both cases
  !call f_profile(test_old,'Test separable iterator',&
  !repeat=nrep,jmpbuf=jmpbuf,dump_results=.true.)
  !call f_profile(test_new,'Original loop',&
  !          repeat=nrep,jmpbuf=jmpbuf,dump_results=.true.)


  denval=0.d0
     
  ! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
  ! possible. In that case the only possible comparison is between the serial and the parallel case
  ! To ease the comparison between the serial and the parallel case we add a random pot_ion
  ! to the potential.
  call f_memcpy(src=density,dest=rhopot)

  offset=0.d0
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
           pot_ion(i1,i2,i3)=tt
           offset=offset+potential(i1,i2,i3)
           !add the case for offset in the surfaces case 
           !(for periodic case it is absorbed in offset)
           if (cell_geocode(mesh) == 'S' .and. denval /= 0.d0) then
              x2 = hy*real(i2-1,kind=8)-0.5d0*acell+0.5d0*hy
              potential(i1,i2,i3)=potential(i1,i2,i3)&
                   -8.d0*datan(1.d0)*denval*real(nspden,kind=8)*(x2**2+0.25d0*acell**2)
              !this stands for
              !denval*2pi*Lx*Lz/Ly^2(y^2-Ly^2/4), less accurate in hgrid
           end if
        end do
     end do
  end do
  if (denval /= 0.d0) density=rhopot
  offset=offset*hx*hy*hz

  !print *,'offset',offset

END SUBROUTINE test_functions_box

subroutine functions(x,a,b,f,f2,whichone)
  implicit none
  integer, intent(in) :: whichone
  real(kind=8), intent(in) :: x,a,b
  real(kind=8), intent(out) :: f,f2
  !local variables
  real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
  real(kind=8) :: length,frequency,nu,sigma,agauss

  pi = 4.d0*datan(1.d0)
  select case(whichone)
  case(1)
     !constant
     f=1.d0
     f2=0.d0
  case(2)
     !gaussian of sigma s.t. a=1/(2*sigma^2)
     r2=a*x**2
     f=dexp(-r2)
     f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
  case(3)
     !gaussian "shrinked" with a=length of the system
     length=a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     f2=factor*dexp(-y**2)
     f=dexp(-y**2)
  case(4)
     !cosine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dcos(r)
     f2=-(frequency*pi/length)**2*dcos(r)
  case(5)
     !exp of a cosine, a=length
     nu=2.d0
     r=pi*nu/a*x
     y=dcos(r)
     yp=dsin(r)
     f=dexp(y)
     factor=(pi*nu/a)**2*(-y+yp**2)
     f2=factor*f
  case(6)
     !gaussian times "shrinked" gaussian, sigma=length/10
     length=a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     g=dexp(-y**2)
     g1=-2.d0*y*yp*g
     g2=factor*dexp(-y**2)

     sigma=length/10
     agauss=0.5d0/sigma**2
     r2=agauss*x**2
     h=dexp(-r2)
     h1=-2.d0*agauss*x*h
     h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
     f=g*h
     f2=g2*h+g*h2+2.d0*g1*h1
  case(7)
     !sine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dsin(r)
     f2=-(frequency*pi/length)**2*dsin(r)
  end select

END SUBROUTINE functions
  
