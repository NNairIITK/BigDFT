!create the descriptors for the density and the potential
subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1i,n2i,n3i,ixc,&
     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)

  use Poisson_Solver

  implicit none
  character(len=1), intent(in) :: geocode,datacode
  integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,ixc
  integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
  integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
  !local variables
  integer :: jproc

  if (datacode == 'D') then
     do jproc=0,iproc-1
        call PS_dim4allocation(geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d            !number of planes for the density
        nscatterarr(jproc,2)=n3p            !number of planes for the potential
        nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
        nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
     end do
     do jproc=iproc+1,nproc-1
        call PS_dim4allocation(geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d
        nscatterarr(jproc,2)=n3p
        nscatterarr(jproc,3)=i3s+i3xcsh-1
        nscatterarr(jproc,4)=i3xcsh
     end do
  end if

  call PS_dim4allocation(geocode,datacode,iproc,nproc,n1i,n2i,n3i,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s)
  nscatterarr(iproc,1)=n3d
  nscatterarr(iproc,2)=n3p
  nscatterarr(iproc,3)=i3s+i3xcsh-1
  nscatterarr(iproc,4)=i3xcsh

  ngatherarr(:,1)=n1i*n2i*nscatterarr(:,2)
  ngatherarr(:,2)=n1i*n2i*nscatterarr(:,3)

end subroutine createDensPotDescriptors


!this routine determine the limits of the coarse grid third dimension which must be
!given to the convolution routines such as to generate the correct portion of fine grid
!compatible with the Poisson Solver data distribution.
!for the moment this is based only on isolated Boundary Conditions
subroutine CoarseForSolverDescriptors(i3s,i3e,n3,i3sc,i3ec)
  implicit none
  integer, intent(in) :: i3s,i3e,n3
  integer, intent(out) :: i3sc,i3ec  
  !local variables
  integer, parameter :: nbl=14,nbr=15 ! the fine dimensions of 0:n3 are -nbl:2*n3+1+nbr
  
  !i3s-1-nbl:i3e-1-nbl should be contained inside -nbl+2*i3sc:2*i3ec+1+nbr
  !thus
  !-nbl+2*i3sc <= i3s-1-nbl  <==> 2*i3sc <= i3s-1
  !2*i3c+1+nbr >= i3e-1-nbl  <==> 2*i3ec >= i3e-2-nbl-nbr
  !and finally
  ! i3sc=(i3s-1)/2
  ! i3ec=(i3e-nbl-nbr-1)/2

  !then we should also calculate the displacement for recuperating the good value from the
  !psir
  !this is
  ! i3s-1-2*i3sc
  !all such values should be inserted into arrays, and rescaled sich that i3sc=0 and i3ec=n3
  !the values of nfl3 and nfu3 should also be shifted and truncated following the determined
  !interval

  !also the bounds arrays should be modified accordingly
  !the kinetic bounds can be left as they are, simply by passing the correct starting address
  !the grow bounds arrays should be modified pointwise
  

end subroutine CoarseForSolverDescriptors

