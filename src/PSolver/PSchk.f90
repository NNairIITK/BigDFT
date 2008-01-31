
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!!****p* BigDFT/PSchk
!! NAME
!!   PSchk
!!
!! FUNCTION
!!    Performs a check of the Poisson Solver suite by running with different regimes
!!    and for different choices of the XC functionals
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2007 CEA
!! CREATION DATE
!!    February 2007
!!
!! SOURCE
!!
program PSchk

  use Poisson_Solver

  implicit none
  include 'mpif.h'
  !Length of the box
  real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
  real(kind=8), parameter :: acell = 10.d0
  character(len=50) :: chain
  character(len=1) :: geocode
  character(len=1) :: datacode
  real(kind=8), dimension(:), allocatable :: density,rhopot,potential,pot_ion,xc_pot
  real(kind=8), pointer :: pkernel(:)
  real(kind=8) :: hx,hy,hz,max_diff,length,eh,exc,vxc,hgrid,diff_parser,offset
  real(kind=8) :: ehartree,eexcu,vexcu,diff_par,diff_ser
  integer :: n01,n02,n03,itype_scf,i_all,i_stat
  integer :: i1_max,i2_max,i3_max,iproc,nproc,ierr,i3sd
  integer :: n_cell,ixc,n3d,n3p,n3pi,i3xcsh,i3s
  integer, dimension(3) :: nxyz

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  !initialize memory counting and timings
  call memocc(0,iproc,'count','start')
  call timing(iproc,'parallel      ','IN')

  !the first proc read the data and then send them to the others
  if (iproc==0) then
     !Use arguments
     call getarg(1,chain)
     read(unit=chain,fmt=*) n01
     call getarg(2,chain)
     read(unit=chain,fmt=*) n02
     call getarg(3,chain)
     read(unit=chain,fmt=*) n03
  end if

  nxyz(1)=n01
  nxyz(2)=n02
  nxyz(3)=n03

  call MPI_BCAST(nxyz,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  n01=nxyz(1)
  n02=nxyz(2)
  n03=nxyz(3)

  print *,iproc,n01,n02,n03

  !Step size
  n_cell = max(n01,n02,n03)
  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)

  !grid for the free BC case
  hgrid=max(hx,hy,hz)
  !hgrid=hx

  !order of the scaling functions choosed
  itype_scf=16

  ixc=11
  geocode='F'


  !calculate the kernel in serial for each processor
  call createKernel(geocode,n01,n02,n03,hx,hy,hz,itype_scf,0,1,pkernel)

  !Allocations
  !Density
  allocate(density(n01*n02*n03),stat=i_stat)
  call memocc(i_stat,product(shape(density))*kind(density),'density','poisson_solver')
  !Density then potential
  allocate(potential(n01*n02*n03),stat=i_stat)
  call memocc(i_stat,product(shape(potential))*kind(potential),'potential','poisson_solver')
  !ionic potential
  allocate(pot_ion(n01*n02*n03),stat=i_stat)
  call memocc(i_stat,product(shape(pot_ion))*kind(pot_ion),'pot_ion','poisson_solver')
  !XC potential
  allocate(xc_pot(n01*n02*n03),stat=i_stat)
  call memocc(i_stat,product(shape(xc_pot))*kind(xc_pot),'xc_pot','poisson_solver')
  allocate(rhopot(n01*n02*n03),stat=i_stat)
  call memocc(i_stat,product(shape(rhopot))*kind(rhopot),'rhopot','poisson_solver')

  !then assign the value of the analytic density and the potential
  call test_functions(geocode,0,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
       density,potential,rhopot,pot_ion)

  !now calculate the Poisson potential in serial for all the processors
  call PSolver(geocode,'G',0,1,n01,n02,n03,ixc,hx,hy,hz,&
       rhopot,pkernel,xc_pot,ehartree,eexcu,vexcu,offset,.false.,1)


  if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies:',ehartree,eexcu,vexcu
  if (iproc == 0) then
     !compare the values of the analytic results
     call compare(0,1,n01,n02,n03,potential,rhopot,&
          i1_max,i2_max,i3_max,max_diff,'ANALYTIC  ')
  end if
  !if the latter test pass, we have a reference for all the other calculations
  !build the reference quantities (based on the numerical result, not the analytic)
  potential=rhopot

  !test for the serial solver
  if (iproc == 0) then
     call compare_with_reference(0,1,geocode,n01,n02,n03,ixc,hx,hy,hz,ehartree,eexcu,vexcu,&
          density,potential,pot_ion,xc_pot,pkernel,rhopot)
  end if

  !now the parallel calculation part
  if (nproc > 1) then
     i_all=-product(shape(pkernel))*kind(pkernel)
     deallocate(pkernel,stat=i_stat)
     call memocc(i_stat,i_all,'pkernel','poisson_solver')

     !calculate the kernel 
     call createKernel(geocode,n01,n02,n03,hx,hy,hz,itype_scf,iproc,nproc,pkernel)

     call compare_with_reference(iproc,nproc,geocode,n01,n02,n03,ixc,hx,hy,hz,ehartree,eexcu,vexcu,&
          density,potential,pot_ion,xc_pot,pkernel,rhopot)
 
  end if
  i_all=-product(shape(pkernel))*kind(pkernel)
  deallocate(pkernel,stat=i_stat)
  call memocc(i_stat,i_all,'pkernel','poisson_solver')

  i_all=-product(shape(rhopot))*kind(rhopot)
  deallocate(rhopot,stat=i_stat)
  call memocc(i_stat,i_all,'rhopot','poisson_solver')
  i_all=-product(shape(density))*kind(density)
  deallocate(density,stat=i_stat)
  call memocc(i_stat,i_all,'density','poisson_solver')
  i_all=-product(shape(potential))*kind(potential)
  deallocate(potential,stat=i_stat)
  call memocc(i_stat,i_all,'potential','poisson_solver')
  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion','poisson_solver')
  i_all=-product(shape(xc_pot))*kind(xc_pot)
  deallocate(xc_pot,stat=i_stat)
  call memocc(i_stat,i_all,'xc_pot','poisson_solver')

  call timing(iproc,'              ','RE')
  !finalize memory counting
  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)  

contains

  subroutine compare_with_reference(iproc,nproc,geocode,n01,n02,n03,ixc,hx,hy,hz,ehref,excref,vxcref,&
       density,potential,pot_ion,xc_pot,pkernel,rhopot)
    use Poisson_Solver
    implicit none
    character(len=1), intent(in) :: geocode
    integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc
    real(kind=8), intent(in) :: hx,hy,hz,ehref,excref,vxcref
    real(kind=8), dimension(n01*n02*n03), intent(in) :: density,potential
    real(kind=8), dimension(n01*n02*n03), intent(inout) :: rhopot,pot_ion,xc_pot
    real(kind=8), dimension(:), pointer :: pkernel
    !local variables
    integer :: n3d,n3p,n3pi,i3xcsh,i3s,istden,istpot,i1_max,i2_max,i3_max,i_all,i_stat
    real(kind=8) :: eexcu,vexcu,offset,max_diff,ehartree
    real(kind=8), dimension(:), allocatable :: test,test_xc


    call PS_dim4allocation(geocode,'D',iproc,nproc,n01,n02,n03,ixc,&
         n3d,n3p,n3pi,i3xcsh,i3s)

    !starting point of the three-dimensional arrays
    istden=n01*n02*(i3s-1)+1
    istpot=n01*n02*(i3s+i3xcsh-1)+1

    !test arrays for comparison
    allocate(test(n01*n02*n03),stat=i_stat)
    call memocc(i_stat,product(shape(test))*kind(test),'test','poisson_solver')
    !XC potential
    allocate(test_xc(n01*n02*n03),stat=i_stat)
    call memocc(i_stat,product(shape(xc_pot))*kind(xc_pot),'xc_pot','poisson_solver')

    if (ixc /= 0) then
       test=potential+pot_ion+xc_pot
    else
       test=potential!+pot_ion
    end if

    rhopot=density
    call PSolver(geocode,'G',iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
         rhopot,pkernel,test_xc,ehartree,eexcu,vexcu,offset,.false.,1)

    !compare the values of the analytic results
    call compare(iproc,nproc,n01,n02,n03,potential,rhopot,&
         i1_max,i2_max,i3_max,max_diff,'ANACOMPLET')

    !compare also the xc_potential
    if (ixc/=0) call compare(iproc,nproc,n01,n02,n03,xc_pot,test_xc,&
         i1_max,i2_max,i3_max,max_diff,'XCCOMPLETE')
    if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

    rhopot=density
    !now we can try with the sumpotion=.true. variable
    call PSolver(geocode,'G',iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
         rhopot,pkernel,pot_ion(istpot),ehartree,eexcu,vexcu,offset,.true.,1)

    !then compare again, but the complete result
    call compare(iproc,nproc,n01,n02,n03,test,rhopot,&
         i1_max,i2_max,i3_max,max_diff,'COMPLETE  ')
    if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

    !now the same thing with the 'D' flag
    rhopot=density
    call PSolver(geocode,'D',iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
         rhopot(istden),pkernel,test_xc,ehartree,eexcu,vexcu,offset,.false.,1)

    !compare the values of the analytic results
    call compare(iproc,nproc,n01,n02,n3p,potential(istpot),rhopot(istpot),&
         i1_max,i2_max,i3_max,max_diff,'ANADISTRIB')

    !compare also the xc_potential
    if (ixc/=0) call compare(iproc,nproc,n01,n02,n3p,xc_pot(istpot),test_xc,&
         i1_max,i2_max,i3_max,max_diff,'XCDISTRIBU')
    if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

    rhopot=density
    !now we can try with the sumpotion=.true. variable
    call PSolver(geocode,'D',iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
         rhopot(istden),pkernel,pot_ion(istpot),ehartree,eexcu,vexcu,offset,.true.,1)

    !then compare again, but the complete result
    call compare(iproc,nproc,n01,n02,n3p,test(istpot),rhopot(istpot),&
         i1_max,i2_max,i3_max,max_diff,'COMPLETEDI')
    if (iproc==0) write(unit=*,fmt="(1x,a,3(1pe20.12))") 'Energies diff:',ehref-ehartree,excref-eexcu,vxcref-vexcu

    i_all=-product(shape(test))*kind(test)
    deallocate(test,stat=i_stat)
    call memocc(i_stat,i_all,'test','poisson_solver')
    i_all=-product(shape(test_xc))*kind(test_xc)
    deallocate(test_xc,stat=i_stat)
    call memocc(i_stat,i_all,'test_xc','poisson_solver')

  end subroutine compare_with_reference

end program PSchk

subroutine compare(iproc,nproc,n01,n02,n03,potential,density,i1_max,i2_max,i3_max,max_diff,description)
  implicit none
  include 'mpif.h'
  character(len=10), intent(in) :: description
  integer, intent(in) :: iproc,nproc,n01,n02,n03
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
  integer, intent(out) :: i1_max,i2_max,i3_max
  real(kind=8), intent(out) :: max_diff
  !local variables
  integer :: i1,i2,i3,ierr
  real(kind=8) :: factor,diff_par
  max_diff = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02 
        do i1=1,n01
           factor=abs(potential(i1,i2,i3)-density(i1,i2,i3))
           if (max_diff < factor) then
              max_diff = factor
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do

!!$  print *,'iproc,i3xcsh,i3s,max_diff',iproc,i3xcsh,i3s,max_diff
  
  if (nproc > 1) then
     !extract the max
     call MPI_ALLREDUCE(max_diff,diff_par,1,MPI_double_precision,  &
          MPI_MAX,MPI_COMM_WORLD,ierr)
  else
     diff_par=max_diff
  end if

  if (iproc == 0) then
     if (nproc == 1) then
        write(unit=*,fmt="(1x,a,1pe20.12)") trim(description)// '    Max diff:',diff_par
        !write(unit=*,fmt="(1x,a,1pe20.12)")'      result:',density(i1_max,i2_max,i3_max),&
        !     '    original:',potential(i1_max,i2_max,i3_max)
        !write(*,'(a,3(i0,1x))')'  Max diff at: ',i1_max,i2_max,i3_max
     else
        write(unit=*,fmt="(1x,a,1pe20.12)") trim(description)//'    Max diff:',diff_par
     end if
  end if

  max_diff=diff_par

end subroutine compare


! this subroutine builds some analytic functions that can be used for 
! testing the poisson solver.
! The default choice is already well-tuned for comparison.
! WARNING: not all the test functions can be used for all the boundary conditions of
! the poisson solver, in order to have a reliable analytic comparison.
! The parameters of the functions must be adjusted in order to have a sufficiently localized
! function in the isolated direction and an explicitly periodic function in the periodic ones.
! Beware of the high-frequency components that may falsify the results when hgrid is too high.
subroutine test_functions(geocode,ixc,n01,n02,n03,acell,a_gauss,hx,hy,hz,&
     density,potential,rhopot,pot_ion)
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,ixc
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03), intent(out) :: density,potential,rhopot,pot_ion

  !local variables
  integer :: i1,i2,i3,nu,ifx,ify,ifz
  real(kind=8) :: x,x1,x2,x3,y,length,denval,pi,a2,derf,hgrid,factor,r,r2
  real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,potion_fac

  if (trim(geocode) == 'P') then
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

!!$     !plot of the functions used
!!$     do i1=1,n03
!!$        x = hx*real(i1,kind=8)!valid if hy=hz
!!$        y = hz*real(i1,kind=8) 
!!$        call functions(x,ax,bx,fx,fx2,ifx)
!!$        call functions(y,az,bz,fz,fz2,ifz)
!!$        write(20,*)i1,fx,fx2,fz,fz2
!!$     end do

     !Initialization of density and potential
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
              density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2
              denval=max(denval,-density(i1,i2,i3))
              potential(i1,i2,i3) = fx*fy*fz!density(i1,i2,i3)
           end do
        end do
     end do

     if (ixc==0) denval=0.d0

  else if (trim(geocode) == 'S') then
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
     bx=real(nu,kind=8)
     bz=real(nu,kind=8)
     !non-periodic dimension
     ay=length!4.d0*a
     by=a
     density(:,:,:) = 0.d0!1d-20 !added

!!$     !plot of the functions used
!!$     do i1=1,n02
!!$        x = hx*real(i1-n02/2-1,kind=8)!valid if hy=hz
!!$        y = hy*real(i1-n02/2-1,kind=8) 
!!$        call functions(x,ax,bx,fx,fx2,ifx)
!!$        call functions(y,ay,by,fy,fy2,ify)
!!$        write(20,*)i1,fx,fx2,fy,fy2
!!$     end do

     !Initialisation of density and potential
     !Normalisation
     do i3=1,n03
        x3 = hz*real(i3-n03/2-1,kind=8)
        call functions(x3,az,bz,fz,fz2,ifz)
        do i2=1,n02
           x2 = hy*real(i2-n02/2-1,kind=8)
           call functions(x2,ay,by,fy,fy2,ify)
           do i1=1,n01
              x1 = hx*real(i1-n02/2-1,kind=8)
              call functions(x1,ax,bx,fx,fx2,ifx)
              density(i1,i2,i3) = fx2*fy*fz+fx*fy2*fz+fx*fy*fz2
              !density(i1,i2,i3) = max(abs(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2),1.d-24)
              !if (abs(density(i1,i2,i3)) < 1.d-20) density(i1,i2,i3)=1.d-20
              !denval=max(denval,-density(i1,i2,i3))
              potential(i1,i2,i3) = -fx*fy*fz*16.d0*datan(1.d0)
           end do
        end do
     end do

     !if (ixc==0) denval=0.d0

  else if (trim(geocode) == 'F') then

     !grid for the free BC case
     !hgrid=max(hx,hy,hz)

     pi = 4.d0*atan(1.d0)
     a2 = a_gauss**2

     !Normalization
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     !gaussian function
     do i3=1,n03
        x3 = hz*real(i3-n03/2,kind=8)
        do i2=1,n02
           x2 = hy*real(i2-n02/2,kind=8)
           do i1=1,n01
              x1 = hx*real(i1-n01/2,kind=8)
              r2 = x1*x1+x2*x2+x3*x3
              density(i1,i2,i3) = max(factor*exp(-r2/a2),1d-24)
              r = sqrt(r2)
              !Potential from a gaussian
              if (r == 0.d0) then
                 potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
              else
                 potential(i1,i2,i3) = derf(r/a_gauss)/r
              end if
           end do
        end do
     end do
     
     denval=0.d0

  else

     print *,'geometry code not admitted',geocode
     stop

  end if

! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
! possible. In that case the only possible comparison is between the serial and the parallel case
! To ease the comparison between the serial and the parallel case we add a random pot_ion
! to the potential.

  if (ixc==0) then
     potion_fac=0.d0
  else
     potion_fac=1.d0
  end if

     rhopot(:,:,:) = density(:,:,:) + denval
     do i3=1,n03
        do i2=1,n02
           do i1=1,n01
              call random_number(tt)
              !tt=0.d0!1.d0
              tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
              pot_ion(i1,i2,i3)=tt
              !potential(i1,i2,i3)=potential(i1,i2,i3)+potion_fac*tt
!!$              !for the ixc/=0 case
!!$              call random_number(tt)
!!$              rhopot(i1,i2,i3)=abs(tt)
           end do
        end do
     end do
     if (denval /= 0.d0) density=rhopot

end subroutine test_functions

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

end subroutine functions

!!***

!!$!fake ABINIT subroutines
!!$subroutine wrtout(unit,message,mode_paral)
!!$  implicit none
!!$
!!$  !Arguments ------------------------------------
!!$  integer,intent(in) :: unit
!!$  character(len=4),intent(in) :: mode_paral
!!$  character(len=500),intent(inout) :: message
!!$
!!$  print *,message
!!$end subroutine wrtout
!!$
!!$subroutine leave_new(mode_paral)
!!$  implicit none
!!$
!!$  !Arguments ------------------------------------
!!$  character(len=4),intent(in) :: mode_paral
!!$
!!$  print *,'exiting...'
!!$  stop
!!$end subroutine leave_new


