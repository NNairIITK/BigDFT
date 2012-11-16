!> @file
!!  Use integral form for Poisson solver
!! @author
!!    Copyright (c) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
program PS_Integral

  use module_base
  use module_types
   use module_interfaces
   use module_xc
   use Poisson_Solver

  implicit none
  integer :: n_points, n_range, n_scf, itype_scf
  integer, dimension(1:7) :: n_points_list = 0
  real(dp) :: hgrid, dx, t0, t1
  real(dp), dimension(7) :: pgauss
  real(dp), dimension(:), allocatable :: x_scf
  real(dp), dimension(:), allocatable :: y_scf
  !real(dp), dimension(-n_range:n_range) :: work
  !real(dp), dimension(-n_range:n_range) :: kernel_scf
  !real(dp), dimension(-n_range:n_range) :: work
  real(dp), dimension(:), allocatable :: work
  !real(dp), dimension(-n_range:n_range) :: kernel_scf
  real(dp), dimension(:), allocatable :: kernel_scf
  real(dp), dimension(:,:), allocatable :: analytic_vs_naive, gaussian
  real(dp), dimension(:,:,:), allocatable :: multiple_naive, analytic_integral_result, timings
  !local variables
  real(dp), parameter :: p0_ref = 1.0_dp
  integer :: n_iter,i_kern,i,i_stat,j,k
  real(dp) :: p0_cell,p0gauss,absci,kern
  character(len=64) :: chain
  !character(len=*) :: chain
  logical :: timings_switch = .false. 
  real(dp), dimension(:), pointer :: fISF
  
  include 'lazy_ISF_8_2048.inc'
  include 'lazy_ISF_14_2048.inc'
  include 'lazy_ISF_16_2048.inc'
  include 'lazy_ISF_20_2048.inc'
  include 'lazy_ISF_24_2048.inc'
  include 'lazy_ISF_30_2048.inc'
  include 'lazy_ISF_40_2048.inc'
  include 'lazy_ISF_50_2048.inc'
  include 'lazy_ISF_60_2048.inc'
  include 'lazy_ISF_100_2048.inc'

  interface 
     subroutine my_analytic_integral(alpha,ntot,m,fwork,fISF,argument_nf)
       use module_base
       implicit none
       integer, intent(in) :: ntot,m
       real(dp), intent(in) :: alpha
       real(dp), dimension(0:ntot), intent(inout) :: fwork
       real(dp), dimension(0:2048), intent(in) :: fISF
       integer, optional, intent(in) :: argument_nf  
     end subroutine my_analytic_integral
  end interface

  !pgauss = 1.0e21_dp
  hgrid = 1.0_dp


  call get_command_argument(1,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Integral itype_scf [timings]'
     stop
  end if
  read(unit=chain,fmt=*) itype_scf
  
  call get_command_argument(2,chain)
  if(trim(chain)=='') then
     timings_switch = .false.
  else if (trim(chain)=='timings') then
     timings_switch = .true.
  else 
     write(*,'(1x,a)')&
       'Usage: ./PS_Integral itype_scf [timings]'
     stop
  end if

  !itype_scf = 14

  !Only itype=8,14,16,20,24,30,40,50,60,100
  select case(itype_scf)
  case(8)
     fISF => fISF8   
  case(14)
     fISF => fISF14
  case(16)
     fISF => fISF16
  case(20)
     fISF => fISF20    
  case(24)
     fISF => fISF24
  case(30)
     fISF => fISF30
  case(40)
     fISF => fISF40
  case(50)
     fISF => fISF50
  case(60)
     fISF => fISF60
  case(100)
     fISF => fISF100
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select

  n_range=2*itype_scf  



  !Allocations
  allocate(work(-n_range:n_range), stat=i_stat)
  allocate(kernel_scf(-n_range:n_range), stat=i_stat)


  !Convergence test for the subroutine gauss_conv_scf 
  !with 6 different values of n_points and 7 different pgauss
  allocate(multiple_naive(1:6,1:7,-n_range:n_range), stat = i_stat) 
  allocate(timings(1:2,1:7,1:7), stat = i_stat)
  timings = 0.0_dp


  
  n_points_list = (/ 2, 8, 32, 64, 256, 512, 0 /) 
  


  do i=1,6
     
     !n_points_list(i) = 8*2**i
     n_points = n_points_list(i)

     !Number of integration points = 2*itype_scf*n_points
     n_scf=2*itype_scf*n_points
     
     !Other allocations
     allocate(x_scf(0:n_scf),stat=i_stat)
     allocate(y_scf(0:n_scf),stat=i_stat)
     allocate(gaussian(1:7,0:n_scf),stat=i_stat)

     !Build the scaling function
     call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)

     !Step grid for the integration
     dx = real(n_range,dp)/real(n_scf,dp)

     do j = 1,7
        !pgauss(j) = real(2.0,dp)**(2*jx-3)*(8.0/real(itype_scf,dp))**2
        pgauss(j) = 1.0e-8_dp*2.0_dp**(6*(j-1))
     end do

     
     if (i == 6) then
        do j = 1,7
           !Test: is the width of the Gaussian ok?
           kern = 0.0_dp
           do k=0,n_scf
              gaussian(j,k) = dexp(-pgauss(j)*(hgrid*x_scf(k))**2) ! <<< For subsequent plotting purposes...
              absci = x_scf(k)
              absci = absci*absci*hgrid**2
              kern = kern + dexp(-pgauss(j)*absci)
           end do
           kern = kern*dx
           write(*,'(a28,1pe9.3,a7,1pe23.16)') 'Gaussian test ---> pgauss = ', pgauss(j), ' err = ', &
                abs(kern-sqrt(4*datan(1.d0)/pgauss(j)))
        end do
        ! plot to file the Gaussian & the scaling function 
        ! >>> only in the case of the highest n_points <<<
        write(chain,'(a18)') 'gauss_scf_plot.out'
        !write(chain,j)
        open(unit=51, file = chain, status = 'unknown')
        write(51, *) '# The columns correspond to x_scf, y_scf, gaussian(pgauss),'
        write(51, '(a18,7(1pe11.3))') '# where pgauss =', pgauss(:)
        write(51, *) ''
        do k = 0, n_scf
           write(51, '(9(1pe21.12e3))') x_scf(k), y_scf(k), &
                gaussian(1,k), gaussian(2,k), gaussian(3,k), &
                gaussian(4,k), gaussian(5,k), gaussian(6,k), gaussian(7,k)
        end do
        close(51) 
     end if


        ! !Test: is the norm of the ISF correct?
        ! kern = 0.0_dp
        ! do k=0,n_scf
        !    kern = kern + y_scf(k)
        ! end do
        ! kern = kern*dx
        ! write(*,*) 'ISF test      --->', pgauss, kern

     do j = 1,7
        !template: gauss_conv_scf(itype_scf,pgauss,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
        call cpu_time(t0)
        if (i == 1 .and. timings_switch) then
           do k = 1,1000
              call gauss_conv_scf(itype_scf, pgauss(j), hgrid, dx, n_range, n_scf, x_scf, y_scf, kernel_scf, work)
           end do
        else
           call gauss_conv_scf(itype_scf, pgauss(j), hgrid, dx, n_range, n_scf, x_scf, y_scf, kernel_scf, work)
        end if
        call cpu_time(t1)
        timings(1,i,j) = t1-t0
        !stores the result in multiple_naive
        multiple_naive(i,j,:) = kernel_scf(:)
     end do

     deallocate(gaussian, stat = i_stat)
     deallocate(x_scf, stat = i_stat)
     deallocate(y_scf, stat = i_stat)

  end do
  
  do j = 1,7
     write(chain,'(a28,i1,a4)') 'gauss_conv_scf_integral_plot', j, '.out'
     !write(chain,j)
     open(unit=51, file = chain, status = 'unknown')
     write(51, *) '# pgauss =', pgauss(j)
     write(51, *) '# 1st column    >>> k in [-n_range:n_range]'
     write(51, *) '# other columns >>> ', &
       'naive integrals with no. of integration points = 2*itype_scf*N, where N =', n_points_list(1:6)
     write(51, *) ' '
     do k = -n_range, n_range
        write(51, '(i4, 6(1pe21.12e3))') k, &
             multiple_naive(1,j,k), & 
             multiple_naive(2,j,k), &
             multiple_naive(3,j,k), &
             multiple_naive(4,j,k), &
             multiple_naive(5,j,k), &
             multiple_naive(6,j,k)
     end do
     close(51)
  end do

 
  ! Evaluate the absolute error with respect to the integral computed with the highest n_points
  do i = 1,5
     do j = 1,7
        do k = -n_range, n_range
           multiple_naive(i,j,k) = abs(multiple_naive(i,j,k)-multiple_naive(6,j,k))
        end do
     end do
  end do



  ! write(chain,'(a26)') 'naive_integral_summary.out'
  ! open(unit=52, file = chain, status = 'replace', position = 'append')
  ! write(52, *) '# 1st column    >>> log10(pgauss)'
  ! write(52, *) '# other columns >>> ', &
  !      'maxval(abs("integral with 2*itype_scf*512 points" - "integral with 2*itype_scf*N points")), where N =', n_points_list(1:5)
  ! write(52, *) ' '
  ! do j = 1,7
  !    write(52, '(1pe10.3,5(1pe20.12e3))') log10(pgauss(j)), &
  !         maxval(multiple_naive(1,j,:)), &
  !         maxval(multiple_naive(2,j,:)), &
  !         maxval(multiple_naive(3,j,:)), &
  !         maxval(multiple_naive(4,j,:)), &
  !         maxval(multiple_naive(5,j,:))
  !    !maxval(multiple_naive(6,j,:))
  ! end do
  ! close(52)

  
  write(chain,'(a26)') 'naive_integral_summary.out'
  open(unit=52, file = chain, status = 'replace', position = 'append')
  write(52, *) '# 1st column    >>> N (N.B.: the no. of integration points is = 2*itype_scf*N)'
  write(52, *) '# other columns >>> ', &
       'maxval(abs("integral with 2*itype_scf*512 points" - "integral with 2*itype_scf*N points"))'
  write(52, *) '# for pgauss =', pgauss(1:7)
  write(52, *) ' '
  do i = 1,5
     !(i3,7(1pe20.12e3)
     write(52, '(i3, 7(1pe21.12e3))') n_points_list(i), &
          maxval(multiple_naive(i,1,:)), &
          maxval(multiple_naive(i,2,:)), &
          maxval(multiple_naive(i,3,:)), &
          maxval(multiple_naive(i,4,:)), &
          maxval(multiple_naive(i,5,:)), &
          maxval(multiple_naive(i,6,:)), &
          maxval(multiple_naive(i,7,:))
          !maxval(multiple_naive(6,j,:))
  end do
  close(52)

 
  ! Convergence test for my_analytic_integral routine


  allocate(analytic_integral_result(1:7,1:7,0:n_range), stat=i_stat)
  allocate(analytic_vs_naive(0:n_range,0:n_range), stat=i_stat)

  
  do i = 1,7
     n_points_list(i) = 2**(i+4)
  end do

  
  open(unit=54, file = 'analytic_vs_naive_summary.out', status = 'replace', position = 'append')
  write(54, *) '# 1st column    >>> N, number of Fourier terms taken into account'
  write(54, *) '# other columns >>> ', &
       'maxval(abs("naive integral with 2*itype_scf*512 integration points" - ', &
       '"analytic_integral with N Fourier terms")) for pgauss =', pgauss(:)
  write(54, *) ' '


  do i = 1,7 ! loop over integrals computed with different numbers of Fourier components
     do j = 1,7 ! loop over different Gaussians
        call cpu_time(t0)
        ! template: my_analytic_integral(alpha,ntot,m,fwork,nf)
        if (i == 4 .and. timings_switch) then
           do k = 1, 1000
              call my_analytic_integral(hgrid*sqrt(pgauss(j)),&
                   n_range, itype_scf, analytic_integral_result(i,j,0:n_range), fISF, n_points_list(i))
           end do
        else 
           call my_analytic_integral(hgrid*sqrt(pgauss(j)),&
                n_range, itype_scf, analytic_integral_result(i,j,0:n_range), fISF, n_points_list(i))
        end if
        call cpu_time(t1)
        timings(2,i,j) = t1-t0
        ! norm infinity with respect to the naive integral computed with the highest n_points
        analytic_vs_naive(i,j) = maxval(abs(analytic_integral_result(i,j,0:n_range) - & 
             multiple_naive(6,j,0:n_range)))
     end do

     write(54,'(i5, 7(1pe21.12e3))') n_points_list(i), analytic_vs_naive(i,1),&
          analytic_vs_naive(i,2),&
          analytic_vs_naive(i,3),&
          analytic_vs_naive(i,4),&
          analytic_vs_naive(i,5),&
          analytic_vs_naive(i,6),&
          analytic_vs_naive(i,7)
  end do

  close(54)


  do j = 1,7
     write(chain,'(a25,i1,a4)') 'my_analytic_integral_plot', j, '.out'
     open(unit=61, file = chain, status = 'unknown')
     write(61, *) '# pgauss =', pgauss(j)
     write(61, *) '# 1st column    >>> k in [0:n_range]'
     write(61, *) '# other columns >>> ', &
       'analytic_integral with N Fourier terms, where N =', n_points_list(:)
     write(61, *) '# The last column contains the result from the naive method, to be used as benchmark.'
     write(61, *) ' '
     do k = 0, n_range
        write(61, '(i3, 8(1pe21.12e3))') k, &
             analytic_integral_result(1,j,k), & 
             analytic_integral_result(2,j,k), &
             analytic_integral_result(3,j,k), &
             analytic_integral_result(4,j,k), &
             analytic_integral_result(5,j,k), &
             analytic_integral_result(6,j,k), &
             analytic_integral_result(7,j,k), & 
             multiple_naive(6,j,k) ! <<< reference
     end do
     close(61)
  end do


  !Timings
  if (timings_switch) then
     open(unit=62, file = 'timings.out', status = 'unknown')
     ! write(62,*) '# 1st column >>> index running over different n_points or nf'
     ! write(62,*) '# 2nd column >>> index running over Gaussian widths'
     ! write(62,*) '# 3rd column >>> naive integral'
     ! write(62,*) '# 4th column >>> analytic integral'
     write(62,*) '# >>> Elapsed times for a bunch of 1,000 integrations <<<'
     write(62,*) '# 1st column >>> index running over Gaussian widths'
     write(62,*) '# 2nd column >>> naive integral with n_points = 16'
     write(62,*) '# 3rd column >>> analytic integral with nf = 256'

     do j = 1,7
        !do i = 1,7
        !   write(62, '(2(i2), 2(1pe15.7))') i, j, timings(1,i,j), timings(2,i,j)
        !end do
        ! timings(1,1,:) >>> n_points = 16 (already enough for achieving an acceptable accuracy)
        ! timings(2,4,:) >>> nf = 256 (enough for achieving an acceptable accuracy)
        write(62, '(1(i2), 2(1pe15.7))') j, timings(1,1,j), timings(2,4,j)
     end do
     close(62)
  end if


  write(*,*) ' '
  write(*,*) 'Please give a look at the output files!'
  write(*,*) "Here's the list:"
  write(*,*) 'naive_integral_summary.out'
  write(*,*) 'analytic_vs_naive_summary.out'
  if (timings_switch) write(*,*) 'timings.out'
  
  write(*,*) ' '
  write(*,*) 'Moreover:'
  write(*,*) 'gauss_scf_plot.out'
  write(*,*) 'gauss_conv_scf_integral_plot(1-7).out'
  write(*,*) 'my_analytic_integral_plot(1-7).out'




  deallocate(work, stat = i_stat)
  deallocate(kernel_scf, stat = i_stat)
  deallocate(x_scf, stat = i_stat)
  deallocate(y_scf, stat = i_stat)
  deallocate(analytic_integral_result, stat=i_stat)
  deallocate(analytic_vs_naive, stat=i_stat)
  deallocate(timings, stat=i_stat)


end program PS_Integral



subroutine gauss_conv_scf(itype_scf,pgauss,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
  use module_base
  implicit none
  integer, intent(in) :: n_range,n_scf,itype_scf
  real(dp), intent(in) :: pgauss,hgrid,dx
  real(dp), dimension(0:n_scf), intent(in) :: x_scf
  real(dp), dimension(0:n_scf), intent(in) :: y_scf
  real(dp), dimension(-n_range:n_range), intent(inout) :: work
  real(dp), dimension(-n_range:n_range), intent(inout) :: kernel_scf
  !local variables
  real(dp), parameter :: p0_ref = 1.0_dp
  integer :: n_iter,i_kern,i
  real(dp) :: p0_cell,p0gauss,absci,kern

  !Step grid for the integration
  !dx = real(n_range,dp)/real(n_scf,dp)

  !To have a correct integration
  p0_cell = p0_ref/(hgrid*hgrid)

  !write(*,*) 'p0_cell = ', p0_cell 

  !We calculate the number of iterations to go from pgauss to p0_ref
  n_iter = nint((log(pgauss) - log(p0_cell))/log(4.0_dp))

  !write(*,*) 'n_iter = ', n_iter

  if (n_iter <= 0)then
     n_iter = 0
     p0gauss = pgauss
  else
     p0gauss = pgauss/4._dp**n_iter
  end if

  !write(*,*) 'p0gauss =',  p0gauss

  !Stupid integration
  !Do the integration with the exponential centered in i_kern
  kernel_scf(:) = 0.0_dp
  do i_kern=0,n_range
     kern = 0.0_dp
     do i=0,n_scf
        absci = x_scf(i) - real(i_kern,dp)
        absci = absci*absci*hgrid**2
        kern = kern + y_scf(i)*dexp(-p0gauss*absci)
     end do
     kernel_scf(i_kern) = kern*dx
     kernel_scf(-i_kern) = kern*dx
     !if (abs(kern) < 1.d-18) then
     !   !Too small not useful to calculate
     !   exit
     !end if
  end do

  !Start the iteration to go from p0gauss to pgauss
  call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,work)

END SUBROUTINE gauss_conv_scf






!! ********************************************* !!


!> Here alpha correspondes to sqrt(alpha) in mathematica
!! the final result is fwork(j+m)-fwork(j-m)
subroutine my_analytic_integral(alpha,ntot,m,fwork,fISF,argument_nf)
  use module_base
  implicit none
  integer, intent(in) :: ntot,m
  real(dp), intent(in) :: alpha
  real(dp), dimension(0:ntot), intent(inout) :: fwork
  real(dp), dimension(0:2048), intent(in) :: fISF
  integer, optional, intent(in) :: argument_nf 
 
  !local variables
  integer :: nf
  real(dp), parameter :: pi=3.1415926535897932384_dp
  logical :: flag,flag1,flag2
  integer :: j,q,jz
  real(dp) :: if,r1,r2,res,ypm,ymm,erfcpm,erfcmm,factor,re,ro,factorend
  

  !write(*,*) fISF(1000), fISF16(1000)


  if(present(argument_nf)) then 
     nf=argument_nf
  else
     nf = 64 ! "default value"
  endif


  flag=.false.
  factor=pi/real(2*m,dp)/alpha
  factorend=sqrt(pi)/alpha/real(4*m,dp)
  !fill work array
  !the calculation for j=0 can be separated from the rest
  !since it only contains error functions
  loop_nonzero: do j=0,ntot
     ypm=alpha*real(j+m,dp)
     ymm=alpha*real(j-m,dp)
     call derfcf(erfcpm,ypm) ! erfcpm = erfc(ypm)
     call derfcf(erfcmm,ymm)
     !assume nf even
     res=0._dp
     !reso=0._dp
     do q=nf,2,-2
        !the sign of q only influences the imaginary part
        !so we can multiply by a factor of two

        !call wofz_mod(alpha,m,-q,-j,r1,if,flag1)
        call wofz_mod(alpha,m,q,-j-m,r1,if,flag1)
        call wofz_mod(alpha,m,q,-j+m,r2,if,flag2)
!!$        call wofz_mod(-xe,-y,r1,if,flag1)
!!$        call wofz_mod(xe,-y,r2,if,flag2)
        re=r1-r2
        flag=flag1 .or. flag2 .or. flag
        !if (flag) then
        !   print *,'here',xe,y,q,j
        !   stop 
        !end if
        !call wofz_mod(alpha,m,-q+1,-j,r1,if,flag1)
        call wofz_mod(alpha,m,q-1,-j-m,r1,if,flag1)
        call wofz_mod(alpha,m,q-1,-j+m,r2,if,flag2)
!!$        call wofz_mod(-xo,-y,r1,if,flag1)
!!$        call wofz_mod(xo,-y,r2,if,flag2)
        ro=r1-r2
        flag=flag1 .or. flag2 .or. flag
        !if (flag) then
        !   print *,'there',xo,y
        !   stop 
        !end if
        !write(16,'(2(i4),6(1pe15.7))')j,q,re,ro,erfcmm-erfcpm
        re=re*fISF(q)
        ro=ro*fISF(q-1)
        res=res+re-ro
     end do
     !q=0 
     !fwork(j)=derf(y)+rese-reso
     fwork(j)=erfcmm-erfcpm+2.0_dp*res!e-reso
     fwork(j)=factorend*fwork(j)
     !exit form the loop if it is close to zero
     if (abs(fwork(j)) < 1.e-25_dp) exit loop_nonzero
     !write(17,'(i4,8(1pe15.7))')j,derf(y),erfcsgn,rese,reso,derf(y)+rese-reso,&
     !     -erfcsgn+rese-reso,erfcsgn+rese-reso
  end do loop_nonzero

  !check flag
  if (flag) then
     write (*,*)'value of alpha',alpha
     stop 'problem occurred in wofz'
  end if

  !put to zero the rest
  do jz=j+1,ntot
     fwork(jz)=0.0_dp
  end do

END SUBROUTINE my_analytic_integral



! >>> GARBAGE <<<


 ! open(unit=54, file = 'analytic_vs_naive_summary.out', status = 'replace', position = 'append')
  ! write(54, *) '# 1st column    >>> log10(pgauss)'
  ! write(54, *) '# other columns >>> ', &
  !      'maxval(abs("naive integral with 512 points" - "analytic_integral with N Fourier terms")), where N =', n_points_list(1:7)
  ! write(54, *) '# Note that n_scf = 2*2*itype_scf*N'
  ! write(54, *) ' '

 

  ! do j = 1,7
  !    ! Convergence test for my_analytic_integral
  !    do i = 1,7
  !       call my_analytic_integral(hgrid*sqrt(pgauss(j)),&
  !            n_range, itype_scf, analytic_integral_result(i, 0:n_range), n_points_list(i))
  !       !write(*,*) 'analytic_integral with nf =', 2**i, ' ---> done!'
  !       analytic_vs_naive(i,j) = maxval(abs(analytic_integral_result(i,0:n_range) - & 
  !            multiple_naive(6,j,0:n_range))) ! NB: norm infinity
  !       !write(*,*)  analytic_integral_result(0), multiple_naive(6,j,0)
  !       !analytic_vs_naive(i-5,j) = analytic_integral_result(0)
  !    end do

  !    write(54,'(1pe8.2e02,7(1pe20.12e3))') log10(pgauss(j)),&
  !         analytic_vs_naive(1,j),&
  !         analytic_vs_naive(2,j),&
  !         analytic_vs_naive(3,j),&
  !         analytic_vs_naive(4,j),&
  !         analytic_vs_naive(5,j),&
  !         analytic_vs_naive(6,j),&
  !         analytic_vs_naive(7,j)

  !    ! write(54,'(i1,7(1pe20.12e3))') i,&
  !    !      analytic_vs_naive(i,1),&
  !    !      analytic_vs_naive(i,2),&
  !    !      analytic_vs_naive(i,3),&
  !    !      analytic_vs_naive(i,4),&
  !    !      analytic_vs_naive(i,5),&
  !    !      analytic_vs_naive(i,6),&
  !    !      analytic_vs_naive(i,7)

 

  ! end do

  ! close(54)
