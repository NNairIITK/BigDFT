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


!this subroutine define other wavefunctions descriptors starting from the original descriptors 
!and the limits of a given localisation region
!it also returns an array which is used to mask the compressed wavefunction into the new one
subroutine loc_wfd(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,wfdg,wfdl,maskarr)
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfdg
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec
  logical, dimension(wfdg%nvctr_c+7*wfdg%nvctr_f), intent(out) :: maskarr
  type(wavefunctions_descriptors), intent(out) :: wfdl
  !local variables
  character(len=*), parameter :: subname='loc_wfd'

  !calculate the number of segments of the new descriptors
  !and define the mask array for the wavefunction

  !coarse part
  call num_segkeys_loc(1,n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
       wfdl%nseg_c,wfdl%nvctr_c,maskarr(1))

  !fine part
  call num_segkeys_loc(7,n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+1),wfdg%keyv(wfdg%nseg_c+1),&
       wfdl%nseg_f,wfdl%nvctr_f,maskarr(wfdg%nvctr_c+1))

  call allocate_wfd(wfdl,subname)

  !a little check
  if (count(maskarr) /= wfdl%nvctr_c+7*wfdl%nvctr_f) then
     write(*,'(1x,a)')'ERROR : Masking problem, check maskarr'
     stop
  end if

  !now fill the local wavefunction descriptors
  !coarse part
  call segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
       wfdl%nseg_c,wfdl%nvctr_c,wfdl%keyg(1,1),wfdl%keyv(1))


  !fine part
  call segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+1),wfdg%keyv(wfdg%nseg_c+1),&
       wfdl%nseg_f,wfdl%nvctr_f,wfdl%keyg(1,wfdl%nseg_c+1),wfdl%keyv(wfdl%nseg_f+1))


end subroutine loc_wfd

subroutine segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(nseg_loc), intent(out) :: keyv_loc
  integer, dimension(2,nseg_loc), intent(out) :: keyg_loc
  !local variables
  logical :: go,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,j,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l
  integer :: ngridp

  !dimensions of the localisation region (O:nIl)
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     go=(i3sc <= i3 .and. i3 <= i3ec) .and. (i2sc <= i2 .and. i2 <= i2ec)
     lseg=.false.
     do i=i0,i1
        i1l=i-i1sc
        i2l=i2-i2sc
        i3l=i3-i3sc
        ngridp=i3l*((n1l+1)*(n2l+1)) + i2l*(n1l+1) + i1l+1
        if (go .and. (i1sc <= i .and. i <= i1ec)) then
           nvctr_check=nvctr_check+1
           if (.not. lseg) then
              nsrt=nsrt+1
              keyg_loc(1,nsrt)=ngridp
              keyv_loc(nsrt)=nvctr_check
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              keyg_loc(2,nend)=ngridp-1
              lseg=.false. 
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
        keyg_loc(2,nend)=ngridp
     end if
  end do

  !check
  if (nvctr_check /= nvctr_loc .or. nend /= nsrt .or. nend /= nseg_loc) then
     write(*,'(1x,a,2(i6))')&
          'ERROR: problem in segkeys_loc',&
          nvctr_check,nvctr_loc,nend,nsrt,nseg_loc
     stop
  end if

end subroutine segkeys_loc

subroutine num_segkeys_loc(num,n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,mask)
  implicit none
  integer, intent(in) :: num,n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
  logical, dimension(num,nvctr), intent(out) :: mask
  !local variables
  logical :: go,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,j,nsrt,nend,nvctr_check

  nvctr_loc=0
  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     go=(i3sc <= i3 .and. i3 <= i3ec) .and. (i2sc <= i2 .and. i2 <= i2ec)
     lseg=.false.
     do i=i0,i1
        !index of the compressed function
        ind=i-i0+jj
        nvctr_check=nvctr_check+1
        if (go .and. (i1sc <= i .and. i <= i1ec)) then
           do j=1,num
              mask(j,ind)=.true.
           end do
           nvctr_loc=nvctr_loc+1
           if (.not. lseg) then
              nsrt=nsrt+1
           end if
           lseg=.true.
        else
           do j=1,num
              mask(j,ind)=.false.
           end do
           if (lseg) then
              nend=nend+1
              lseg=.false. 
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
     end if
  end do
  nseg_loc=nend  

  !check
  if (nend /= nsrt) then 
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif

  if (nvctr_check /= nvctr) then
     write(*,'(1x,a,2(i6))')&
          'ERROR: incorrect number of coarse points examined for reducing the localisation region',&
          nvctr_check,nvctr
     stop
  end if
    
end subroutine num_segkeys_loc


!with this subroutine we should convert a bound array from its global version to 
!the version which
subroutine bound_conversion(nl1,nu1,nl2,nu2,nl3,nu3)
end subroutine bound_conversion
