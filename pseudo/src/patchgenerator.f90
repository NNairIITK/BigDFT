
subroutine paw_generator(izatom,zion, lmx,  lpmx, lmax,  hsep, gpot, &
     alpz, alps, &
     ng, noccmax,noccmx, expo,&
     psi, aeval, occup, &
     Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid ,rw,rd,  psigrid, Npaw,&
     PAWpatch , psipsigrid, rcov, rprb, rcore, zcore , Ngrid_box_larger)

  implicit none
  integer, intent(in) :: izatom, ng,noccmax,noccmx,Nsol, labs, Ngrid,  Ngrid_box, Ngrid_box_larger

  integer, intent(in) :: Npaw
  
  real(8), dimension(ng+1), intent(out) :: expo
  integer , intent(in)::lpmx, lmx, lmax

  integer, parameter :: n_int=1000

  real(8), intent(in) :: rgrid(Ngrid), rd(Ngrid),rw(Ngrid), rcore, zcore

  real(8), dimension(0:ng,noccmax, lmx), intent(out) :: psi, Egrid(Nsol),&
       psigrid(Ngrid,Nsol  )
  
  real(8), dimension(4), intent(in) :: gpot !! gpot dim e diventata 4!!!
  real(8),   intent(out)  :: psipsigrid(Ngrid,Nsol  )
  real(8), intent(in) :: rcov, rprb, zion
  


  real(8), dimension(noccmx,lmx  ), intent(in) ::  occup
  real(8), dimension(noccmx,lmx  ), intent(out) ::  aeval
  real(8), intent(out):: PAWpatch(Npaw,Npaw)
  real(8), intent(in):: hsep(6,lpmx)
  real(8), intent(in) :: alpz, alps(*)
 
  !local variables
  real(8) alpl
  real(8), parameter :: fact=4.0_8
  real(8), dimension(noccmx,lmx) ::chrg,res
  real(8), dimension(:), allocatable :: xp
  real(8), dimension(:,:), allocatable :: vh

  real(8), dimension(:,:,:,:), allocatable :: rmt
  integer :: l,i,iocc,i_all,i_stat,  j 
  real(8) :: alrcov, rij,a,a0,a0in,tt,ehomo
  real(8) :: value_at_r
  integer :: igrid, isol, lpx
  logical :: pawisactive

  !filename = 'psppar.'//trim(atomname)

  lpx=0

  lpx_determination: do i=1,4
     if (alps(i) == 0.0_8) then
        exit lpx_determination
     else
        lpx=i-1
     end if
  end do lpx_determination

  alpl=alpz
  
  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2))
  
  allocate(xp(0:ng))
  allocate(rmt(n_int,0:ng,0:ng,lmax+1))
  
  !can be switched on for debugging
  !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
  !     'Input Guess Generation for atom',trim(atomname),&
  !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb
  
  rij=3._8
  ! exponents of gaussians
  ! print *, " ESPONENTI " 
  ! print *, " alpz " , alpz
  a0in=alpz
  a0=a0in/rij
  !       tt=sqrt(sqrt(2._8))
  tt=2._8**(.3_8)
  do i=0,ng
     a=a0*tt**i
     xp(i)=.5_8/a**2
     ! print *, " xp(", i,")", xp(i)
  end do
  
  ! initial guess
  do l=0,lmx-1
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.0_8
        end do
     end do
  end do
  
  call crtvh_paw(ng,lmax,xp,vh,rprb,fact,n_int,rmt)
  
!!!  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
!!!       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
!!!       aeval,ng,psi,res,chrg)
  
  PAWpatch=0.0_8
  call gatom_modified(rcov,rprb,lmax,lpx,lpmx, noccmax,noccmx, occup,&
       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
       aeval,ng,psi,res,chrg,&
       Nsol, Labs, Ngrid,Ngrid_box,Egrid,  rgrid,rw, rd,  psigrid,Npaw,  PAWpatch,&
       psipsigrid,rcore,zcore , Ngrid_box_larger   )              
 
  
  do i=1,ng+1
     expo(i)=sqrt(0.5_8/xp(i-1))
  end do
  
  do l=0,lmx-1
     do iocc=1,noccmax
        if( value_at_r(rprb, ng , expo,psi(0,iocc,l+1)).lt.0.0     ) then
           do i=0,ng
              psi(i,iocc,l+1)=-psi(i,iocc,l+1)
           enddo
        endif
     enddo
  enddo
  
   
END SUBROUTINE paw_generator

function value_at_r(r, ng , expo,psi     )

  implicit none
  integer, parameter :: gp=kind(1.0d0) 
 
  real(gp) , intent(in) ::  r
  integer, intent(in) :: ng
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=1000

  integer ig
  real(gp) sum
  real(gp) :: value_at_r

  sum=0.0
  do ig = 0,ng
     sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
  enddo
  value_at_r=sum

end function value_at_r
 
