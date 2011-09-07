
subroutine paw_generator(izatom,ielpsp, lpmx, hsep, gpot, &
     alpz, alps, &
     ng, noccmax,expo,&
     psi, aeval, occup, &
     Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw,  PAWpatch , psipsigrid )

  use module_base, only: gp, memocc,ndebug
  use module_interfaces
  implicit none
  integer, intent(in) :: izatom,ielpsp,ng,noccmax, Nsol, labs, Ngrid,  Ngrid_box

  integer, intent(in) :: Npaw
  
  real(gp), dimension(ng+1), intent(out) :: expo
  integer , intent(in)::lpmx

  integer, parameter :: n_int=1000

  real(gp), intent(in) :: rgrid(Ngrid)

  real(gp), dimension(0:ng,noccmax, lpmx), intent(out) :: psi, Egrid(Nsol),&
       psigrid(Ngrid,Nsol  )


  real(gp),   intent(out), optional  :: psipsigrid(Ngrid,Nsol  )
  real(gp), dimension(noccmax,lpmx  ), intent(in) ::  occup
  real(gp), dimension(noccmax,lpmx  ), intent(out) ::  aeval
  real(gp), intent(out):: PAWpatch(Npaw,Npaw)
  real(gp), intent(in):: hsep(6,lpmx)
  real(gp), intent(in) :: alpz, alps(:)
 
  !local variables
  real(gp), parameter :: fact=4.0_gp
 real(gp), dimension(4) :: gpot !! gpot dim e diventata 4!!!
  real(gp), dimension(noccmax,lpmx) ::chrg,res
  real(gp), dimension(:), allocatable :: xp
  real(gp), dimension(:,:), allocatable :: vh

  real(gp), dimension(:,:,:,:), allocatable :: rmt
  integer :: nsccode,mxpl,mxchg
  integer :: l,i,iocc,i_all,i_stat,  j 
  real(gp) :: alrcov,rprb,zion,rij,a,a0,a0in,tt,ehomo
  real(gp) :: value_at_r,amu
  integer :: igrid, isol
  logical :: pawisactive
  !filename = 'psppar.'//trim(atomname)

  lpx=0
  lpx_determination: do i=1,4
     if (psppar(i,0) == 0.0_gp) then
        exit lpx_determination
     else
        lpx=i-1
     end if
  end do lpx_determination


  alpl=alpz
  
  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2+ndebug),stat=i_stat)
  call memocc(i_stat,vh,'vh',subname)
  
  allocate(xp(0:ng+ndebug),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(rmt(n_int,0:ng,0:ng,lpmx),stat=i_stat)
  call memocc(i_stat,rmt,'rmt',subname)
  
  !can be switched on for debugging
  !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
  !     'Input Guess Generation for atom',trim(atomname),&
  !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb
  
  rij=3._gp
  ! exponents of gaussians
  ! print *, " ESPONENTI " 
  ! print *, " alpz " , alpz
  a0in=alpz
  a0=a0in/rij
  !       tt=sqrt(sqrt(2._gp))
  tt=2._gp**(.3_gp)
  do i=0,ng
     a=a0*tt**i
     xp(i)=.5_gp/a**2
     ! print *, " xp(", i,")", xp(i)
  end do
  
  ! initial guess
  do l=0,lpmx-1
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.0_gp
        end do
     end do
  end do
  
  call crtvh(ng,lpmx-1,xp,vh,rprb,fact,n_int,rmt)
  
!!!  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
!!!       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
!!!       aeval,ng,psi,res,chrg)
  
  
  if(.not. present(psipsigrid) ) then
     call gatom_modified(rcov,rprb,lpmx-1,lpx,noccmax,occup,&
          zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
          aeval,ng,psi,res,chrg,&
          Nsol, Labs, Ngrid,Ngrid_box,Egrid,  rgrid , psigrid,Npaw,  PAWpatch )
  else
     print *, "chiamo gatom_modified con psipsigrid "
     PAWpatch=0.0_gp
     call gatom_modified(rcov,rprb,lpmx-1,lpx,noccmax,occup,&
          zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
          aeval,ng,psi,res,chrg,&
          Nsol, Labs, Ngrid,Ngrid_box,Egrid,  rgrid , psigrid,Npaw,  PAWpatch,&
          psipsigrid)           
  endif
  
  !! the operation below ensure that AE and pseudo wf 
  !! coincides for big r
  !! In the case of paw fitting, however, psigrid is the dual
  !! ptilde, while psitilde has been fitted to psigrid from input
  !! In this latter case the operation below is no more necessary
  !! ( Moreover ptilde might have another  sign than psitilde due the duality operation ?)
  pawisactive=.false.
  do i=1, Npaw
     do j=1, Npaw
        if( PAWpatch(i,j) /= 0.0_gp ) then
           pawisactive= .true.
        endif
     end do
  end do
  if( .not. pawisactive) then
     do isol=1,nsol
        if( psigrid(Ngrid_box-1, isol)<0) then
           do igrid=1, ngrid
              psigrid(igrid, isol)=-psigrid(igrid, isol)
           end do
        endif
     enddo
  endif
  
  
  if (psp_modifier/=1 ) then
     !post-treatment of the inguess data
     do i=1,ng+1
        expo(i)=sqrt(0.5_gp/xp(i-1))
     end do
     
     do l=0,lpmx-1
        do iocc=1,noccmax
           if( value_at_r(rprb, ng , expo,psi(0,iocc,l+1)).lt.0.0     ) then
              do i=0,ng
                 psi(i,iocc,l+1)=-psi(i,iocc,l+1)
              enddo
           endif
        enddo
     enddo
     
  endif
  
  
END SUBROUTINE paw_generator

