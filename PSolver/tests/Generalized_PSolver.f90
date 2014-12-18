!!
!! Solver for the Generalized Poisson Equation. 
!! @author Giuseppe Fisicaro

program GPS_3D

   use wrapper_mpi
   use Poisson_Solver
   use yaml_output
   use dynamic_memory
   use dictionaries
   use time_profiling
   use f_utils
   implicit none
   
   integer, parameter :: n01 = 200
   integer, parameter :: n02 = 200
   integer, parameter :: n03 = 200
   integer, parameter :: nspden = 1
   character(len=40) :: PSol !, parameter :: PSol='PI'!'PCG' !    PCG = Preconditioned Conjugate Gradient Method.
                                             !    PI  = Polarization Iterative Method.
   !> To set 1 for Vacuum, 2 for error function in y direction, 3 for gaussian, 4 for electron dependence.
   integer, parameter :: SetEps = 3 
   real(kind=8), parameter :: acell = 10.d0
   real(kind=8), parameter :: erfL = 80.0d0 ! To set 1 for Vacuum and correct analitic comparison with gaussian potential.
   real(kind=8), parameter :: erfR = 1.0d0
   real(kind=8), parameter :: sigmaeps = 0.05d0*acell
   !> To set 1 for gaussian rho, 2 for gaussian potential. Beaware to use 2 or 3 when SetEps is setted to 3!!! 
   integer, parameter :: Setrho = 2 
   character(len=2), parameter :: geocode = 'F'
   character(len=2), parameter :: datacode = 'G'
   !> Order of accuracy for derivatives into ApplyLaplace subroutine = Total number of points at left and right of the x0 where we want to calculate the derivative.
   integer, parameter :: nord = 16 
   integer, dimension(3) :: ndims
   real(8), dimension(3) :: hgrids
   real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
   !integer :: m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,n1,n2,n3,
   integer :: itype_scf,i_all,i_stat,n_cell,iproc,nproc,ixc
   real(kind=8) :: hx,hy,hz,freq,fz,fz1,fz2,pi,curr,average,CondNum,wcurr,ave1,ave2,rhores2,En1,En2,dVnorm,hgrid
   real(kind=8) :: Adiag,ersqrt,ercurr,factor,r,r2,max_diff,max_diffpot,fact,x1,x2,x3,derf_tt,diffcurr,diffcurrS,divprod,sum
   real(kind=8) :: ehartree,offset
   real(kind=8), dimension(:,:,:,:), allocatable :: density,rhopot,rvApp
   real(kind=8), pointer :: kernel(:)
   type(coulomb_operator) :: pkernel
!   type(mpi_environment), intent(in), optional :: mpi_env
   real(kind=8), dimension(:,:,:), allocatable :: eps,potential,pot_ion
   integer :: i1,i2,i3,whichone,i,ii,j,info,icurr,ip,isd,i1_max,i2_max,i3_max,n3d,n3p,n3pi,i3xcsh,i3s,n3pr2,n3pr1,ierr
!   type(mpi_environment) :: bigdft_mpi
  type(dictionary), pointer :: options

 

   call f_lib_initialize()

   !read command line
   call PS_Check_command_line_options(options)

   call f_zero(PSol)
   PSol=options .get. 'method'
   call dict_free(options)

   call yaml_set_stream(record_length=92,tabbing=30)!unit=70,filename='log.yaml')
   call yaml_new_document()


!   call MPI_INIT(ierr)
!   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
!   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!   bigdft_mpi%mpi_comm=MPI_COMM_WORLD !workaround to be removed

   !allocate arrays to avoid stack overflow
   density=f_malloc([n01,n02,n03,nspden],id='density')
   rhopot =f_malloc([n01,n02,n03,nspden],id='rhopot')
   rvApp  =f_malloc([n01,n02,n03,nspden],id='rvApp')
   
   eps=f_malloc([n01,n02,n03],id='eps')
   potential=f_malloc([n01,n02,n03],id='potential')
   pot_ion=f_malloc([n01,n02,n03],id='pot_ion')

   n_cell = max(n01,n02,n03)
   hx=acell/real(n01,kind=8)
   hy=acell/real(n02,kind=8)
   hz=acell/real(n03,kind=8)
   hgrid=max(hx,hy,hz)

   ixc=0
   itype_scf=16
   iproc=0
   nproc=1
!-------------------------------------------------------------------------

   ! Create the Kernel.
   ! Calculate the kernel in parallel for each processor.

   ndims=(/n01,n02,n03/)
   hgrids=(/hx,hy,hz/)

   pkernel=pkernel_init(.true.,iproc,nproc,0,geocode,ndims,hgrids,itype_scf)

   call pkernel_set(pkernel,.true.)
   !call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)

!------------------------------------------------------------------------

!------------------------------------------------------------------------
   ! Set environment, namely permittivity epsilon.

    call SetEpsilon(n01,n02,n03,eps,acell,a_gauss,hx,hy,hz,erfL,erfR,sigmaeps,SetEps)

!------------------------------------------------------------------------

   ! Set initial density, and the associated analitical potential for the Standard Poisson Equation.
   call SetInitDensPot(n01,n02,n03,nspden,eps,sigmaeps,erfR,acell,a_gauss,a2,hx,hy,hz,Setrho,density,potential)
!  call SetRhoSoluto(n03,rhosol,acell)

!------------------------------------------------------------------------

   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(n01,n02,n03,nspden,potential,rvApp,acell,eps,nord)

   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)

!!$   max_diffpot = 0.d0
!!$   i1_max = 1
!!$   i2_max = 1
!!$   i3_max = 1
!!$   do i3=1,n03
!!$    do i2=1,n02
!!$     do i1=1,n01
!!$     fact=abs(density(i1,i2,i3,1) - rvApp(i1,i2,i3,1))
!!$      if (max_diffpot < fact) then
!!$       max_diffpot = fact
!!$       i1_max = i1
!!$       i2_max = i2
!!$       i3_max = i3
!!$      end if
!!$     end do
!!$    end do
!!$   end do
!!$   write(*,*)'Max abs difference between analytic density - ApplyLaplace to potential'
!!$   write(*,'(3(1x,I4),2(1x,e14.7))')i1_max,i2_max,i3_max,max_diffpot,abs(density(n01/2,n02/2,n03/2,1)-rvApp(n01/2,n02/2,n03/2,1))

  rhopot(:,:,:,:) = density(:,:,:,:)

  if ( trim(PSol)=='PCG') then

   call Prec_conjugate_gradient(n01,n02,n03,nspden,rhopot,acell,eps,nord,pkernel,potential)

  else if (trim(PSol)=='PI') then

   call PolarizationIteration(n01,n02,n03,nspden,rhopot,acell,eps,nord,pkernel,potential)

  end if

  call pkernel_free(pkernel)
  call f_free(density)
  call f_free(rhopot)
  call f_free(rvApp)
  call f_free(eps)
  call f_free(potential)
  call f_free(pot_ion)

  call yaml_release_document()
  call yaml_close_all_streams()

  
  call f_lib_finalize()
  
contains

  !>identify the options from command line
  !! and write the result in options dict
  subroutine PS_Check_command_line_options(options)
    use yaml_parse
    use dictionaries
    implicit none
    !> dictionary of the options of the run
    !! on entry, it contains the options for initializing
    !! on exit, it contains in the key "BigDFT", a list of the 
    !! dictionaries of each of the run that the local instance of BigDFT
    !! code has to execute
    type(dictionary), pointer :: options
    !local variables
    type(yaml_cl_parse) :: parser !< command line parser

    !define command-line options
    parser=yaml_cl_parse_null()
    !between these lines, for another executable using BigDFT as a blackbox,
    !other command line options can be specified
    !then the bigdft options can be specified
    call PS_check_options(parser)
    !parse command line, and retrieve arguments
    call yaml_cl_parse_cmd_line(parser,args=options)
    !free command line parser information
    call yaml_cl_parse_free(parser)

  end subroutine PS_Check_command_line_options

end program GPS_3D

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

end subroutine PS_Check_options


subroutine PolarizationIteration(n01,n02,n03,nspden,b,acell,eps,nord,pkernel,potential)
  use yaml_output
  use Poisson_Solver

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  type(coulomb_operator), intent(in) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps,potential
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b

  real(kind=8), dimension(n01,n02,n03)  :: pot_ion
  real(kind=8), parameter :: eta = 1.0d0 ! Polarization Iterative Method parameter.
  real(kind=8), parameter :: taupol = 1.0d-20 ! Polarization Iterative Method parameter.
  integer, parameter :: maxiterpol=50
  real(kind=8), dimension(n01,n02,n03,nspden,3) :: dlv
  real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(n01,n02,n03,nspden) :: rhosol,rhopol,rhotot,rhopolnew,rhopolold,rhores,lv
  integer :: i1,i2,i3,i,j,ip,isp
  real(kind=8) :: divprod,rhores2,hx,hy,hz,offset,diffcurr,pi,ehartree

  pi = 4.d0*datan(1.d0)
  hx = acell/real(n01,kind=8)
  hy = acell/real(n02,kind=8)
  hz = acell/real(n03,kind=8)


  !write(*,'(a)')'--------------------------------------------------------------------------------------------'
  !write(*,'(a)')'Starting Polarization Iteration '

  open(unit=18,file='PolConvergence.out',status='unknown')
  open(unit=38,file='MaxAnalysisPI.out',status='unknown')


  call fssnord3DmatNabla3var(n01,n02,n03,nspden,eps,deps,nord,acell)

  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     rhopol(i1,i2,i3,isp) = 0.d0
     rhosol(i1,i2,i3,isp) = b(i1,i2,i3,isp)
    end do
   end do
  end do


    call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
!    call writeroutine(n01,n02,n03,nspden,rhosol,0)

    call yaml_sequence_open('Embedded PSolver, Polarization Iteration Method')

    do ip=1,maxiterpol

       !write(*,'(a)')'--------------------------------------------------------------------------------------------'
       !write(*,*)'Starting PI iteration ',ip

     isp=1
     do i3=1,n03
      do i2=1,n02
       do i1=1,n01
        rhotot(i1,i2,i3,isp)=(1/eps(i1,i2,i3))*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
        lv(i1,i2,i3,isp) = rhotot(i1,i2,i3,isp)
       end do
      end do
     end do
     
     call yaml_sequence(advance='no')
     call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

     call writeroutinePot(n01,n02,n03,nspden,lv,ip,potential)
     call fssnord3DmatNabla(n01,n02,n03,nspden,lv,dlv,nord,acell)

     isp=1
     do i3=1,n03
      do i2=1,n02
       do i1=1,n01
        divprod = 0.d0
        do j=1,3
           divprod = divprod + deps(i1,i2,i3,j)*dlv(i1,i2,i3,isp,j)
        end do
        rhopolnew(i1,i2,i3,isp)=(1/(4.d0*pi))*(1/eps(i1,i2,i3))*divprod
        rhopolold(i1,i2,i3,isp)=rhopol(i1,i2,i3,isp)
        rhopol(i1,i2,i3,isp)=eta*rhopolnew(i1,i2,i3,isp) + (1.d0-eta)*rhopolold(i1,i2,i3,isp)
        rhores(i1,i2,i3,isp) = rhopol(i1,i2,i3,isp) - rhopolold(i1,i2,i3,isp)
       end do
      end do
     end do

     rhores2 = 0.d0
     isp=1
     do i3=1,n03
      do i2=1,n02
       do i1=1,n01
        rhores2 = rhores2 + rhores(i1,i2,i3,isp)*rhores(i1,i2,i3,isp)
       end do
      end do
     end do

     write(18,'(1x,I8,1x,e14.7)')ip,rhores2
     !write(*,'(1x,I8,1x,e14.7)')ip,rhores2

     call EPS_iter_output(ip,0.0_dp,rhores2,0.0_dp,0.0_dp,0.0_dp)
     if (rhores2.lt.taupol) exit

!     call writeroutine(n01,n02,n03,nspden,rhores,ip)

  end do

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       b(i1,i2,i3,isp) = lv(i1,i2,i3,isp)
      end do
     end do
    end do

    call yaml_sequence_close()

    close(unit=18)
    close(unit=38)
    !write(*,*)
    !write(*,'(1x,a,1x,i8)')'Polarization iterations =',ip
    !write(*,'(1x,a,1x,e14.7)')'rhores polarization square =',rhores2
    !write(*,*)
    !write(*,'(a)')'Max abs difference between analytic potential and the computed one'
    call writeroutinePot(n01,n02,n03,nspden,b,ip,potential)
    !write(*,*)
    !write(*,'(a)')'Termination of Polarization Iteration'
    !write(*,'(a)')'--------------------------------------------------------------------------------------------'



end subroutine PolarizationIteration

subroutine Prec_conjugate_gradient(n01,n02,n03,nspden,b,acell,eps,nord,pkernel,potential)

  use Poisson_Solver
  use yaml_output

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  type(coulomb_operator), intent(in) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b

  real(kind=8), dimension(n01,n02,n03,nspden) :: x,r,z,p,q,qold,lv,corr
  real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(n01,n02,n03) :: de2,ddeps
  integer, parameter :: max_iter = 50
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k
  integer :: i,ii,j,i1,i2,i3,isp
  real(kind=8), parameter :: error = 1.0d-20
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,offset,pi

  pi = 4.d0*datan(1.d0)   

  open(unit=18,file='PCGConvergence.out',status='unknown')
  open(unit=38,file='MaxAnalysisPCG.out',status='unknown')

  call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')

  !write(*,'(a)')'--------------------------------------------------------------------------------------------'
  !write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  !write(*,'(a)')'Starting PCG iteration 1'

!------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnordEpsilonDerivative(n01,n02,n03,nspden,eps,de2,ddeps,nord,acell)
  call fssnord3DmatNabla3varde2(n01,n02,n03,nspden,eps,deps,de2,nord,acell)
  call fssnord3DmatDiv3var(n01,n02,n03,nspden,deps,ddeps,nord,acell)

  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     corr(i1,i2,i3,isp)=(-0.125d0/pi)*(0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
    end do
   end do
  end do

!------------------------------------------------------------------------------------
! Apply the Preconditioner

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
  normb=dsqrt(normb)

  call yaml_sequence(advance='no')
  call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     p(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do

!------------------------------------------------------------------------------------
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction

  beta=0.d0
  k=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     q(i1,i2,i3,isp)=b(i1,i2,i3,isp)+p(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
     beta=beta+b(i1,i2,i3,isp)*p(i1,i2,i3,isp)
     k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
    end do
   end do
  end do

!------------------------------------------------------------------------------------

  alpha = beta/k
  !write(*,*)alpha
  normr=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     x(i1,i2,i3,isp) = alpha*p(i1,i2,i3,isp)
     r(i1,i2,i3,isp) = b(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
     normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
     lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
  normr=dsqrt(normr)

  ratio=normr/normb

  call EPS_iter_output(1,normb,normr,ratio,alpha,beta)
!  call writeroutine(n01,n02,n03,nspden,r,1)
  call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
  call writeroutinePot(n01,n02,n03,nspden,x,1,potential)

  write(18,'(1x,I8,2(1x,e14.7))')1,ratio,beta
  !write(*,'(1x,I8,2(1x,e14.7))')1,ratio,beta

  do i=2,max_iter

   if (ratio.lt.error) exit
   if (ratio.gt.max_ratioex) exit

   !write(*,'(a)')'--------------------------------------------------------------------------------------------!'
   !write(*,*)'Starting PCG iteration ',i

!  Apply the Preconditioner

   call yaml_sequence(advance='no')
   call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do

   k=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
      q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
      qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
      k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha

   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
     end do
    end do
   end do
   normr=dsqrt(normr)

   ratio=normr/normb
   write(18,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)

  end do

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
  call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
  write(*,*)

  close(unit=18)
  close(unit=38)

  !write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
  !write(*,'(a)')'--------------------------------------------------------------------------------------------'

end subroutine  Prec_conjugate_gradient

subroutine EPS_iter_output(iter,normb,normr,ratio,alpha,beta)
  use module_defs, only: dp
  use yaml_output
  implicit none
  integer, intent(in) :: iter
  real(dp), intent(in) :: normb,normr,ratio,beta,alpha

  call yaml_mapping_open('Iteration quality',flow=.true.)
  call yaml_comment('Iteration '//trim(yaml_toa(iter)),hfill='_')
  !write the PCG iteration
  call yaml_map('iter',iter,fmt='(i4)')
  !call yaml_map('rho_norm',normb)
  if (normr/=0.0_dp) call yaml_map('res',normr,fmt='(1pe16.4)')
  if (ratio /= 0.0_dp) call yaml_map('ratio',ratio,fmt='(1pe16.4)')
  if (alpha /= 0.0_dp) call yaml_map('alpha',alpha,fmt='(1pe16.4)')
  if (beta /= 0.0_dp) call yaml_map('beta',beta,fmt='(1pe16.4)')

  call yaml_mapping_close()
end subroutine EPS_iter_output


subroutine writeroutine(n01,n02,n03,nspden,r,i)

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: i
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: r
  integer :: i1,i2,i3,j,i1_max,i2_max,i3_max,jj
  real(kind=8) :: max_val,fact

     j=i+40
     jj=i+100

    if (i.le.50) then

     i3=n03/2
     do i2=1,n02
      do i1=1,n01
       write(j,'(2(1x,I4),1x,e14.7)') i1,i2,r(i1,i2,i3,1)
      end do
      write(j,*) 
     end do

      do i1=1,n01
!       write(jj,'(1x,I4,2(1x,e22.15))') i1,r(i1,n02/2,i3,1),r(i1,1,i3,1)
       write(jj,'(1x,I4,2(1x,e22.15))') i1,r(i1,i1,i1,1),r(i1,n02/2,i3,1)
      end do

    end if

      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               fact=abs(r(i1,i2,i3,1))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do
      write(39,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
           r(n01/2,n02/2,n03/2,1),r(2,n02/2,n03/2,1),r(10,n02/2,n03/2,1)

end subroutine writeroutine

subroutine writeroutinePot(n01,n02,n03,nspden,ri,i,potential)
  use yaml_output
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: i
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: ri
  real(kind=8), dimension(n01,n02,n03),intent(in) :: potential
  !automatic array, to be check is stack poses problem
  real(kind=8), dimension(n01,n02,n03,nspden) :: re
  integer :: i1,i2,i3,j,i1_max,i2_max,i3_max,jj
  real(kind=8) :: max_val,fact

      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               re(i1,i2,i3,1) = ri(i1,i2,i3,1) - potential(i1,i2,i3)
               fact=abs(re(i1,i2,i3,1))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do
      write(38,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
           re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)
      !write(*,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
      !     re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)
      if (max_val == 0.d0) then
         call yaml_map('Inf. Norm difference with reference',0.d0)
      else
         call yaml_mapping_open('Inf. Norm difference with reference')
         call yaml_map('Value',max_val,fmt='(1pe22.15)')
         call yaml_map('Point',[i1_max,i2_max,i3_max],fmt='(i4)')
         call yaml_map('Some values',[re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)],&
              fmt='(1pe22.15)')
         call yaml_mapping_close()
      end if

end subroutine writeroutinePot

subroutine FluxSurface(n01,n02,n03,nspden,x,acell,eps,nord)
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), dimension(n01,n02,n03), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03,nspden,3) :: dx
  real(kind=8) :: hx,hy,hz,pi,flux
  integer :: i1,i2,i3,isp,i
   pi = 4.d0*datan(1.d0)

  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)

   call fssnord3DmatNabla(n01,n02,n03,nspden,x,dx,nord,acell)

     flux=0.d0
      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01

         if (i1.eq.1) then
          flux=flux-eps(i1,i2,i3)*dx(i1,i2,i3,isp,1)*hy*hz
         end if
         if (i1.eq.n01) then
          flux=flux+eps(i1,i2,i3)*dx(i1,i2,i3,isp,1)*hy*hz
         end if

         if (i2.eq.1) then
          flux=flux-eps(i1,i2,i3)*dx(i1,i2,i3,isp,2)*hx*hz
         end if
         if (i2.eq.n02) then
          flux=flux+eps(i1,i2,i3)*dx(i1,i2,i3,isp,2)*hx*hz
         end if

         if (i3.eq.1) then
          flux=flux-eps(i1,i2,i3)*dx(i1,i2,i3,isp,3)*hx*hy
         end if
         if (i3.eq.n03) then
          flux=flux+eps(i1,i2,i3)*dx(i1,i2,i3,isp,3)*hx*hy
         end if

        end do
       end do
      end do

    write(*,'(1x,a,1x,e14.7)')'Surface flux is',flux

end subroutine FluxSurface 

subroutine ApplyLaplace(n01,n02,n03,nspden,x,y,acell,eps,nord)
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  real(kind=8), dimension(n01,n02,n03,nspden) :: ddx
  real(kind=8), dimension(n01,n02,n03,nspden,3) :: dx
  real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8) :: hx,hy,hz,pi
  integer :: i1,i2,i3,isp,i
   pi = 4.d0*datan(1.d0)   

  hx=acell/real(n01,kind=8)
  hy=acell/real(n02,kind=8)
  hz=acell/real(n03,kind=8)

   call fssnord3DmatNabla(n01,n02,n03,nspden,x,dx,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
         do i=1,3
          dx(i1,i2,i3,isp,i)=eps(i1,i2,i3)*dx(i1,i2,i3,isp,i)
         end do
        end do
       end do
      end do

   call fssnord3DmatDiv(n01,n02,n03,nspden,dx,y,nord,acell)

   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)

end subroutine ApplyLaplace

subroutine fssnord3DmatNabla(n01,n02,n03,nspden,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03,nspden) :: u
      real(kind=8), dimension(n01,n02,n03,nspden,3) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,isp,i1_max,i2_max
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
      real(kind=8) :: hx,hy,hz,max_diff,fact

      n = nord+1
      m = nord/2
      hx = acell/real(n01,kind=8)
      hy = acell/real(n02,kind=8)
      hz = acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c1DF(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,isp,1) = 0.0d0

             if (i1.le.m) then
              do j=-m,m
               du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3,isp)!/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,isp)!/hx
              end do
             else
              do j=-m,m
               du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(i1 + j,i2,i3,isp)!/hx
              end do
             end if
             du(i1,i2,i3,isp,1)=du(i1,i2,i3,isp,1)/hx

             du(i1,i2,i3,isp,2) = 0.0d0

             if (i2.le.m) then
             do j=-m,m
              du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3,isp)!/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,isp)!/hy
              end do
             else
              do j=-m,m
               du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,i2 + j,i3,isp)!/hy
              end do
             end if
              du(i1,i2,i3,isp,2)=du(i1,i2,i3,isp,2)/hy

             du(i1,i2,i3,isp,3) = 0.0d0

             if (i3.le.m) then
             do j=-m,m
              du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1,isp)!/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,isp)!/hz
              end do
             else
              do j=-m,m
               du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,i3 + j,isp)!/hz
              end do
             end if
              du(i1,i2,i3,isp,3)=du(i1,i2,i3,isp,3)/hz

            end do
         end do
      end do

end subroutine fssnord3DmatNabla

subroutine fssnord3DmatNabla3var(n01,n02,n03,nspden,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03) :: u
      real(kind=8), dimension(n01,n02,n03,3) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: hx,hy,hz

      n = nord+1
      m = nord/2
      hx = acell/real(n01,kind=8)
      hy = acell/real(n02,kind=8)
      hz = acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,1) = 0.0d0

             if (i1.le.m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3)/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)/hx
              end do
             else
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(i1 + j,i2,i3)/hx
              end do
             end if

             du(i1,i2,i3,2) = 0.0d0

             if (i2.le.m) then
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3)/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)/hy
              end do
             else
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,i2 + j,i3)/hy
              end do
             end if

             du(i1,i2,i3,3) = 0.0d0

             if (i3.le.m) then
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1)/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)/hz
              end do
             else
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,i3 + j)/hz
              end do
             end if

            end do
         end do
      end do

end subroutine fssnord3DmatNabla3var

subroutine fssnord3DmatNabla3varde2(n01,n02,n03,nspden,u,du,du2,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03) :: u
      real(kind=8), dimension(n01,n02,n03,3) :: du
      real(kind=8), dimension(n01,n02,n03) :: du2

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: hx,hy,hz

      n = nord+1
      m = nord/2
      hx = acell/real(n01,kind=8)
      hy = acell/real(n02,kind=8)
      hz = acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,1) = 0.0d0
             du2(i1,i2,i3) = 0.0d0

             if (i1.le.m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3)/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)/hx
              end do
             else
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(i1 + j,i2,i3)/hx
              end do
             end if

             du2(i1,i2,i3) = du(i1,i2,i3,1)*du(i1,i2,i3,1)
             du(i1,i2,i3,2) = 0.0d0

             if (i2.le.m) then
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3)/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)/hy
              end do
             else
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,i2 + j,i3)/hy
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + du(i1,i2,i3,2)*du(i1,i2,i3,2)

             du(i1,i2,i3,3) = 0.0d0

             if (i3.le.m) then
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1)/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)/hz
              end do
             else
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,i3 + j)/hz
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + du(i1,i2,i3,3)*du(i1,i2,i3,3)

            end do
         end do
      end do

end subroutine fssnord3DmatNabla3varde2

subroutine fssnordEpsilonDerivative(n01,n02,n03,nspden,u,du2,ddu,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first and second derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.


      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u
      real(kind=8), dimension(n01,n02,n03), intent(out) :: du2,ddu

      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,i1_max,i2_max
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c2D
      real(kind=8) :: hx,hy,hz,d,dd

      n = nord+1
      m = nord/2
      hx = acell/real(n01,kind=8)
      hy = acell/real(n02,kind=8)
      hz = acell/real(n03,kind=8)
      n_cell = min(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first and second derivative coefficients from 'Matematica'.
      !Only nord=2,4,6,8,16 admitted.

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c2D(i,j)=0.d0
       end do
      end do


       include 'FiniteDiffCorff.inc'
       include 'FiniteDiffCorff_2der.inc'

      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du2(i1,i2,i3) = 0.0d0
             ddu(i1,i2,i3) = 0.0d0
             d=0.d0
             dd=0.d0

             if (i1.le.m) then
              do j=-m,m
               d = d + c1D(j,i1-m-1)*u(j+m+1,i2,i3)
               dd = dd + c2D(j,i1-m-1)*u(j+m+1,i2,i3)
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               d = d + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
               dd = dd + c2D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
              end do
             else
              do j=-m,m
               d = d + c1D(j,0)*u(i1 + j,i2,i3)
               dd = dd + c2D(j,0)*u(i1 + j,i2,i3)
              end do
             end if

             du2(i1,i2,i3) = d*d*(real(n01,kind=8)/acell)**2
             ddu(i1,i2,i3) = dd*(real(n01,kind=8)/acell)**2

             d=0.d0
             dd=0.d0

             if (i2.le.m) then
             do j=-m,m
              d = d + c1D(j,i2-m-1)*u(i1,j+m+1,i3)
              dd = dd + c2D(j,i2-m-1)*u(i1,j+m+1,i3)
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               d = d + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
               dd = dd + c2D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
              end do
             else
              do j=-m,m
               d = d + c1D(j,0)*u(i1,i2 + j,i3)
               dd = dd + c2D(j,0)*u(i1,i2 + j,i3)
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d*(real(n02,kind=8)/acell)**2
             ddu(i1,i2,i3) = ddu(i1,i2,i3) + dd*(real(n02,kind=8)/acell)**2

             d=0.d0
             dd=0.d0

             if (i3.le.m) then
             do j=-m,m
              d = d + c1D(j,i3-m-1)*u(i1,i2,j+m+1)
              dd = dd + c2D(j,i3-m-1)*u(i1,i2,j+m+1)
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               d = d + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
               dd = dd + c2D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
              end do
             else
              do j=-m,m
               d = d + c1D(j,0)*u(i1,i2,i3 + j)
               dd = dd + c2D(j,0)*u(i1,i2,i3 + j)
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d*(real(n03,kind=8)/acell)**2
             ddu(i1,i2,i3) = ddu(i1,i2,i3) + dd*(real(n03,kind=8)/acell)**2

            end do
         end do
      end do

end subroutine fssnordEpsilonDerivative

subroutine fssnord3DmatDiv(n01,n02,n03,nspden,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03,nspden,3) :: u
      real(kind=8), dimension(n01,n02,n03,nspden) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,isp
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: hx,hy,hz,d1,d2,d3
      real(kind=8), parameter :: zero = 0.d0! 1.0d-11

      n = nord+1
      m = nord/2
      hx = acell/real(n01,kind=8)
      hy = acell/real(n02,kind=8)
      hz = acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,isp) = 0.0d0

             d1 = 0.d0
             if (i1.le.m) then
              do j=-m,m
               d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,isp,1)!/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,isp,1)!/hx
              end do
             else
              do j=-m,m
               d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,isp,1)!/hx
              end do
             end if
              d1=d1/hx

             d2 = 0.d0
             if (i2.le.m) then
             do j=-m,m
              d2 = d2 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,isp,2)!/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               d2 = d2 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,isp,2)!/hy
              end do
             else
              do j=-m,m
               d2 = d2 + c1D(j,0)*u(i1,i2 + j,i3,isp,2)!/hy
              end do
             end if
              d2=d2/hy

             d3 = 0.d0
             if (i3.le.m) then
             do j=-m,m
              d3 = d3 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,isp,3)!/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               d3 = d3 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,isp,3)!/hz
              end do
             else
              do j=-m,m
               d3 = d3 + c1D(j,0)*u(i1,i2,i3 + j,isp,3)!/hz
              end do
             end if
              d3=d3/hz

             du(i1,i2,i3,isp) = d1+d2+d3

            end do
         end do
      end do

end subroutine fssnord3DmatDiv

subroutine fssnord3DmatDiv3var(n01,n02,n03,nspden,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03,3) :: u
      real(kind=8), dimension(n01,n02,n03) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,isp
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: hx,hy,hz,d1,d2,d3
      real(kind=8), parameter :: zero = 0.d0! 1.0d-11

      n = nord+1
      m = nord/2
      hx = acell/real(n01,kind=8)
      hy = acell/real(n02,kind=8)
      hz = acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3) = 0.0d0

             d1 = 0.d0
             if (i1.le.m) then
              do j=-m,m
               d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,1)!/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)!/hx
              end do
             else
              do j=-m,m
               d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,1)!/hx
              end do
             end if
              d1=d1/hx

             d2 = 0.d0
             if (i2.le.m) then
             do j=-m,m
              d2 = d2 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,2)!/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               d2 = d2 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)!/hy
              end do
             else
              do j=-m,m
               d2 = d2 + c1D(j,0)*u(i1,i2 + j,i3,2)!/hy
              end do
             end if
              d2=d2/hy

             d3 = 0.d0
             if (i3.le.m) then
             do j=-m,m
              d3 = d3 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,3)!/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               d3 = d3 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)!/hz
              end do
             else
              do j=-m,m
               d3 = d3 + c1D(j,0)*u(i1,i2,i3 + j,3)!/hz
              end do
             end if
              d3=d3/hz

             du(i1,i2,i3) = d1+d2+d3

            end do
         end do
      end do

end subroutine fssnord3DmatDiv3var

subroutine SetInitDensPot(n01,n02,n03,nspden,eps,sigmaeps,erfR,acell,a_gauss,a2,hx,hy,hz,Setrho,density,potential)
  use yaml_output
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: Setrho
  real(kind=8), intent(in) :: acell,a_gauss,a2,hx,hy,hz,sigmaeps,erfR
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: density
  real(kind=8), dimension(n01,n02,n03), intent(out) :: potential
  real(kind=8), dimension(n01,n02,n03,nspden) :: density1,density2
  real(kind=8), dimension(n01,n02,n03) :: potential1,potential2
  integer :: i,i1,i2,i3
  real(kind=8) :: sigma,sigma1,sigma2,pi,sum,sump,tt1,tt2,x0,r12
  real(kind=8) :: x1,x2,x3,r,r2,r1,r22,derf_tt1,derf_tt2,factor,factor1,factor2

  pi = 4.d0*datan(1.d0)

  if (Setrho.eq.1) then
! Set initial density as gaussian (or double gaussian with zero total charge) and potential as error function. It works only
! in a vacuum environment.

   sigma1 = 0.03d0*acell
   sigma2 = 2.d0*sigma1
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor1 = 1.d0/((sigma1**3)*sqrt((2.d0*pi)**3))
!         factor2 = 1.d0/((sigma2**3)*sqrt((2.d0*pi)**3))
         factor2 = 0.d0! 1.d0/((sigma2**3)*sqrt((2.d0*pi)**3))
         !gaussian function for the density.
         sum=0.d0
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
                  r12=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r22 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                     density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*max(factor1*exp(-0.5d0*r12/(sigma1**2)),1d-24)&
                                         - 1.d0/real(nspden,kind=8)*max(factor2*exp(-0.5d0*r22/(sigma2**2)),1d-24)
                  sum=sum+density(i1,i2,i3,i)
                  end do
                  r1 = sqrt(r12)
                  r2 = sqrt(r22)
                  !Potential from a gaussian
                  if (r1 == 0.d0) then
                     potential(i1,i2,i3) = 2.d0/(sqrt(2.d0*pi)*sigma1)
                  else
                     call derf_local(derf_tt1,r1/(sqrt(2.d0)*sigma1))
                     potential(i1,i2,i3) = derf_tt1/r1
                  end if
                  sump=sump+potential(i1,i2,i3)
!                  if (r2 == 0.d0) then
!                     potential(i1,i2,i3) = potential(i1,i2,i3) - 2.d0/(sqrt(2.d0*pi)*sigma2)
!                  else
!                     call derf_local(derf_tt2,r2/(sqrt(2.d0)*sigma2))
!                     potential(i1,i2,i3) = potential(i1,i2,i3) - derf_tt2/r2
!                  end if
               end do
            end do
         end do

  else if (Setrho.eq.2) then

! Set initial potential as gaussian and density as the correct Generalized Laplace operator. It works with a gaussian epsilon.

   sigma = 0.03d0*acell
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))
         !gaussian function for the potential.
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  potential(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24)
                   sump=sump+potential(i1,i2,i3)
               end do
            end do
          end do

         sum=0.d0
         !analitic density calculation.
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                   density(i1,i2,i3,i) =(-1.d0/(4.d0*pi))*potential(i1,i2,i3)*(1.d0/(sigma**2))*&
                                        (r2*(1.d0/(sigmaeps**2))*(eps(i1,i2,i3)-erfR)+eps(i1,i2,i3)*(r2*(1.d0/(sigma**2))-3.d0))
                   sum=sum+density(i1,i2,i3,i)
                  end do
               end do
            end do
          end do

  else if (Setrho.eq.3) then

! Set initial potential as double gaussian and density as the correct 
!!Generalized Laplace operator. It works with a gaussian epsilon.

!   sigma1 = 0.033d0*acell
!   sigma2 = 2.d0*sigma1
   sigma2 = 0.033d0*acell
   sigma1 = 2.d0*sigma2
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor1 = 1.d0/((sigma1**3)*sqrt((2.d0*pi)**3))
         factor2 = 1.d0/((sigma2**3)*sqrt((2.d0*pi)**3))
         !gaussian function for the potential.
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  potential1(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor1*exp(-0.5d0*r2/(sigma1**2)),1d-24)
                  potential2(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor2*exp(-0.5d0*r2/(sigma2**2)),1d-24)
                  potential(i1,i2,i3) = potential1(i1,i2,i3) - potential2(i1,i2,i3)
                  sump=sump+potential(i1,i2,i3)
               end do
            end do
          end do

         sum=0.d0
         !analitic density calculation.
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                   density1(i1,i2,i3,i) = (-1.d0/(4.d0*pi))*potential1(i1,i2,i3)*(1.d0/(sigma1**2))*&
                                        (r2*(1.d0/(sigmaeps**2))*(eps(i1,i2,i3)-erfR)+eps(i1,i2,i3)*(r2*(1.d0/(sigma1**2))-3.d0))
                   density2(i1,i2,i3,i) = (-1.d0/(4.d0*pi))*potential2(i1,i2,i3)*(1.d0/(sigma2**2))*&
                                        (r2*(1.d0/(sigmaeps**2))*(eps(i1,i2,i3)-erfR)+eps(i1,i2,i3)*(r2*(1.d0/(sigma2**2))-3.d0))
                   density(i1,i2,i3,i) = density1(i1,i2,i3,i) - density2(i1,i2,i3,i)

                   sum=sum+density(i1,i2,i3,i)
                  end do
               end do
            end do
          end do


  end if

  call yaml_map('Total Charge',sum*hx*hy*hz)
  call yaml_map('Potential monopole',sump*hx*hy*hz)
  call yaml_map('Potential at the boundary 1 n02/2 1',&
       potential(1,n02/2,1))
  call yaml_map('Density at the boundary 1 n02/2 1',density(1,n02/2,1,1))
  !write(*,*) 'charge sum',sum*hx*hy*hz,'potential sum',sump*hx*hy*hz
  !write(*,'(1x,a,1x,e14.7)')'Potential at the boundary 1 n02/2 1',poteantial(1,n02/2,1)
  !write(*,'(1x,a,1x,e14.7)')'Density at the boundary 1 n02/2 1',density(1,n02/2,1,1)

end subroutine SetInitDensPot

subroutine SetEpsilon(n01,n02,n03,eps,acell,a_gauss,hx,hy,hz,erfL,erfR,sigmaeps,SetEps)
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz,erfL,erfR,sigmaeps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(out) :: eps
  real(kind=8), dimension(n01,n02,n03) :: edens
  integer :: i,i1,i2,i3
  real(kind=8), parameter :: edensmax = 1.5d0
  real(kind=8), parameter :: edensmin = 0.01d0
  real(kind=8), parameter :: eps0 = 80.0d0
  real(kind=8) :: x1,x2,x3,r,t,pi,r2,sigma,x0,factor
    
!  open(unit=21,file='Epsilon.out',status='unknown')

  if (SetEps.eq.1) then

   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      eps(i1,i2,i3) = 1.d0
     end do
    end do
   end do

  else if (SetEps.le.2) then

   do i2=1,n02
    x2 = hy*real(i2-n02/2,kind=8)
    do i3=1,n03
     do i1=1,n01

      eps(i1,i2,i3) = (erfL-erfR)*(1.d0-erf((x2)/(acell*0.1d0)))/2.d0 + erfR
     
     end do
    end do
   end do

!     i3=1!n03/2
!     do i2=1,n02
!      do i1=1,n01
!       write(21,'(2(1x,I4),2(1x,e14.7))')i1,i2,eps(i1,i2,i3),eps(i1,i2,n03/2)
!      end do
!      write(21,*)
!     end do

  else if (SetEps.le.3) then

   sigma = sigmaeps
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor = erfL - erfR
         !factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))
         !gaussian function
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  !r2 = x1*x1+x2*x2+x3*x3
                     eps(i1,i2,i3) = max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24) + erfR
               end do
            end do
         end do
  
!     i3=1!n03/2
!     do i2=1,n02
!      do i1=1,n01
!       write(21,'(2(1x,I4),2(1x,e14.7))')i1,i2,eps(i1,i2,i3),eps(i1,i2,n03/2)
!      end do
!      write(21,*)
!     end do

!     do i1=1,n01
!      write(35,'(1x,I8,2(1x,e22.15))') i1,eps(i1,n02/2,1),eps(i1,1,1)
!     end do

   else if (SetEps.eq.4) then

   call SetEledens(n01,n02,n03,edens,acell,a_gauss,hx,hy,hz)

   pi = 4.d0*datan(1.d0)
   r=0.d0
   t=0.d0

   do i3=1,n03
    do i2=1,n02
     do i1=1,n01

      if (dabs(edens(i1,i2,i3)).gt.edensmax) then
       eps(i1,i2,i3)=1.d0
      else if (dabs(edens(i1,i2,i3)).lt.edensmin) then
       eps(i1,i2,i3)=eps0
      else
       r=2.d0*pi*(dlog(edensmax)-dlog(dabs(edens(i1,i2,i3))))/(dlog(edensmax)-dlog(edensmin))
       t=((dlog(eps0))/(2.d0*pi))*(r-dsin(r))
       eps(i1,i2,i3)=dexp(t)
      end if

     end do
    end do
   end do
    
  end if

  close(unit=21)

end subroutine SetEpsilon

subroutine SetEledens(n01,n02,n03,edens,acell,a_gauss,hx,hy,hz)
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03), intent(out) :: edens

  real(kind=8) :: x1,x2,x3,freq,fz,fz1,fz2,average,pi,a2,factor,r2
  integer :: i,i1,i2,i3

!  open(unit=18,file='ElectronDensity.out',status='unknown')

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
     edens(i1,i2,i3) = max(factor*exp(-r2/a2),1d-24)
    end do
   end do
  end do

! Here averaged to zero in the old one dimensional code!

!  close(unit=18)

end subroutine SetEledens

subroutine derf_local(derf_yy,yy)

 implicit none
 integer, parameter :: dp=8
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
