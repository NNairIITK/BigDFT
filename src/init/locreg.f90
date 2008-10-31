
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


!create a routine which transpose the wavefunctions into a form ready for the density construction
!starting from the mask arrays
!this routine should be rethought, the maskarray is much too big
subroutine rhotranspose(iproc,nproc,norbp,norb,nspinor,wfd,wfd_loc,psi_loc,&
     maskrho,ncountrhoreg,psit)
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd,wfd_loc
  integer , intent(in) :: iproc,nproc,norbp,norb,nspinor
  logical, dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: maskrho
  integer, dimension(nproc), intent(in) :: ncountrhoreg
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp*nspinor), intent(in) :: psi_loc
  !this array may also have dimension zero
  real(wp), dimension(wfd_loc%nvctr_c+7*wfd_loc%nvctr_f,norb*nspinor), intent(out) :: psit
  !local variables
  character(len=*), parameter :: subname='rhotranspose'
  integer :: i_stat,i_all,iorb

  !first rearrange the wavefunction in each of thel ocalisation regions which belong to a
  !processor
!!$  call rhoswitch_waves(nlr,norbp,nvctr_c,nvctr_f,nseg_tot,nvctr_tot,nseglr,keymask,&
!!$     ncount,psi,psiloc)

  

end subroutine rhotranspose

!pass from the global wavefunction for a set of orbitals to the local wavefunctions
!for all the orbitals
!INPUTS
! nlr           number of localisation regions contained by the compressed form
! nvctr_c,f     compressed elements of the global function
! nseglr     array of compressed elements for each localisation region
!               number of segments of keymask array
! nvctr_tot   number of components of the distributed wavefunction
! nseg_tot      total number of segments for the mask array
! keymask       array of masking element for each localisation region           
! psi           global wavefunctions with distributed orbitals
!OUTPUTS
! psiloc        wavefunctions of all the localisation regions 
! ncount        array of elements to be sent to each processor for MPI_ALLTOALLV procedure
subroutine rhoswitch_waves(nlr,norbp,nvctr_c,nvctr_f,nseg_tot,nvctr_tot,nseglr,keymask,&
     ncount,psi,psiloc)
  use module_base
  implicit none
  integer, intent(in) :: nlr,nvctr_c,nvctr_f,nseg_tot,nvctr_tot,norbp
  integer, dimension(nlr,2) , intent(in) :: nseglr
  integer, dimension(2,nseg_tot), intent(in) :: keymask
  real(wp), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
  integer, dimension(nlr), intent(out) :: ncount
  real(wp), dimension(nvctr_tot), intent(out) :: psiloc
  !local variables
  integer :: jsh,jsegsh,ilr,iseg,iorb,i,j,k,jseg
  
  jsh=0
  jsegsh=0
  do ilr=1,nlr
     j=0
     do iorb=1,norbp
        jseg=0 !each orbital has the same descriptors
        do iseg=1,nseglr(ilr,1) !coarse segments
           jseg=jseg+1
           do i=keymask(1,jseg+jsegsh),keymask(2,jseg+jsegsh)
              j=j+1
              psiloc(j+jsh)=psi(i,iorb)
           end do
        end do
        do iseg=1,nseglr(ilr,2) !fine segments
           jseg=jseg+1
           do i=keymask(1,jseg+jsegsh),keymask(2,jseg+jsegsh)
              do k=1,7
                 j=j+1
                 psiloc(j+jsh)=psi(k+nvctr_c+7*(i-1),iorb)
              end do
           end do
        end do
     end do
     ncount(ilr)=j
     jsh=jsh+ncount(ilr)
     jsegsh=jsegsh+j
  end do
end subroutine rhoswitch_waves


!build the wavefunction in real space starting from its compressed form
!one of the building blocks for a localised region call
subroutine build_psir(d,wfd,gb,ibyz_c,ibyyzz_r,w1,w2,x_c,x_f,scal,psi,psir)
  use module_base
  use module_types
  implicit none
  type(grid_dimensions), intent(in) :: d
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(grow_bounds), intent(in) :: gb
  integer, dimension(2,0:d%n2,0:d%n3), intent(in) :: ibyz_c
  integer, dimension(2,-14:2*d%n2+16,-14:2*d%n3+16), intent(in):: ibyyzz_r
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(max(4*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)*(2*(d%nfu1-d%nfl1)+31),&
       (2*d%n1+31)*(d%n2+1)*(d%n3+1))), &
       intent(inout) :: w1 !work
  real(wp), dimension(max((d%n3+1)*(2*d%n1+31)*(2*d%n2+31),&
       2*(d%nfu3-d%nfl3+1)*(2*(d%nfu1-d%nfl1)+31)*(2*(d%nfu2-d%nfl2)+31))), &
       intent(inout) :: w2 ! work
  real(wp), dimension(0:d%n1,0:d%n2,0:d%n3), intent(inout) :: x_c
  real(wp), dimension(7,d%nfl1:d%nfu1,d%nfl2:d%nfu2,d%nfl3:d%nfu3), intent(inout) :: x_f
  real(wp), dimension(d%n1i,d%n2i,d%n3i), intent(out) :: psir
  
  call uncompress_forstandard_short(d%n1,d%n2,d%n3,d%nfl1,d%nfu1,d%nfl2,d%nfu2,d%nfl3,d%nfu3, & 
       wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
       wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
       scal,psi(1),psi(wfd%nvctr_c+1),x_c,x_f)
  
  call comb_grow_all(d%n1,d%n2,d%n3,d%nfl1,d%nfu1,d%nfl2,d%nfu2,d%nfl3,d%nfu3,w1,w2,x_c,x_f,  & 
       psir,ibyz_c,gb%ibzxx_c,gb%ibxxyy_c,gb%ibyz_ff,gb%ibzxx_f,gb%ibxxyy_f,ibyyzz_r)
  
end subroutine build_psir


!this subroutine define other wavefunctions descriptors starting from the original descriptors 
!and the limits of a given localisation region
!it also returns an array which is used to mask the compressed wavefunction into the new one
! INPUTS
! ilocreg               localisation region to be considered
! nlocreg               total number of localisation regions
! n1,n2,n3              original dimensions of the global box
! lrlims                array of limits of the localisation regions (global system coordinates)
! wfdg                  global wavefunction descriptors structure
! OUTPUT
! wfdl                  local wavefunction descriptors structure in local system coordinates
! keymask               mask array for traducing the wavefunction in compressed form
!                       to the wavefunction in compressed form for the local system
! ncountlocreg          array of elements for each localisation region
subroutine loc_wfd(ilocreg,nlocreg,n1,n2,n3,lrlims,wfdg,wfdl,keymask,ncountlocreg)
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfdg
  integer, intent(in) :: ilocreg,nlocreg,n1,n2,n3
  integer, dimension(2,3,nlocreg), intent(in) :: lrlims
  type(wavefunctions_descriptors), intent(out) :: wfdl
  integer, dimension(nlocreg), intent(out) :: ncountlocreg
  integer, dimension(:), pointer :: keymask
  !local variables
  character(len=*), parameter :: subname='loc_wfd'
  integer :: i_stat,i_all
  integer :: iloc,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nvctr_c,nseg_c,nseg_f,nvctr_f,ndimkey

  !calculate the number of segments of the new descriptors for each localisation region
  !and the dimension of the array for the translation of the localisation regions
  ndimkey=0
  do iloc=1,nlocreg
     if (iloc /= ilocreg) then
        !coarse part
        call num_segkeys_loc(n1,n2,n3,lrlims(1,1,iloc),lrlims(2,1,iloc),&
             lrlims(1,2,iloc),lrlims(2,2,iloc),lrlims(1,3,iloc),lrlims(2,3,iloc),&
             wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
             nseg_c,nvctr_c)
        !fine part
        call num_segkeys_loc(n1,n2,n3,lrlims(1,1,iloc),lrlims(2,1,iloc),&
             lrlims(1,2,iloc),lrlims(2,2,iloc),lrlims(1,3,iloc),lrlims(2,3,iloc),&
             wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+1),wfdg%keyv(wfdg%nseg_c+1),&
             nseg_f,nvctr_f)
        ncountlocreg(iloc)=nvctr_c+7*nvctr_f
        ndimkey=ndimkey+nseg_c+nseg_f
     end if
  end do

  i1sc=lrlims(1,1,ilocreg)
  i1ec=lrlims(2,1,ilocreg)
  i2sc=lrlims(1,2,ilocreg)
  i2ec=lrlims(2,2,ilocreg)
  i3sc=lrlims(1,3,ilocreg)
  i3ec=lrlims(2,3,ilocreg)

  !coarse part
  call num_segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
       wfdl%nseg_c,wfdl%nvctr_c)

  !fine part
  call num_segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+1),wfdg%keyv(wfdg%nseg_c+1),&
       wfdl%nseg_f,wfdl%nvctr_f)

  ncountlocreg(ilocreg)=wfdl%nvctr_c+7*wfdl%nvctr_f
  ndimkey=ndimkey+wfdl%nseg_c+wfdl%nseg_f

  call allocate_wfd(wfdl,subname)

  allocate(keymask(ndimkey),stat=i_stat)
  call memocc(i_stat,keymask,'keymask',subname)

  !now fill the local wavefunction descriptors
  !and define the mask array for the wavefunction
  !coarse part
  call segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
       wfdl%nseg_c,wfdl%nvctr_c,wfdl%keyg(1,1),wfdl%keyv(1),keymask(1))

  !fine part
  call segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+1),wfdg%keyv(wfdg%nseg_c+1),&
       wfdl%nseg_f,wfdl%nvctr_f,wfdl%keyg(1,wfdl%nseg_c+1),wfdl%keyv(wfdl%nseg_c+1),&
       keymask(wfdg%nseg_c+1))

  !a little check on the masking array
!!$  if (count(maskarr) /= wfdl%nvctr_c+7*wfdl%nvctr_f) then
!!$     write(*,'(1x,a)')'ERROR : Masking problem, check maskarr'
!!$     stop
!!$  end if

end subroutine loc_wfd

subroutine build_keymask(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg_tot,keyg,keyv,&
     nseg_loc,keymask)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg_tot,nseg_loc
  integer, dimension(nseg_tot), intent(in) :: keyv
  integer, dimension(2,nseg_tot), intent(in) :: keyg
  integer, dimension(2,nseg_loc), intent(out) :: keymask
  !local variables
  logical :: go,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,j,nsrt,nend

  !start and end points
  nsrt=0
  nend=0
  do iseg=1,nseg_tot
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
        if (go .and. (i1sc <= i .and. i <= i1ec)) then
           if (.not. lseg) then
              nsrt=nsrt+1
              keymask(1,nsrt)=ind
           end if
           lseg=.true.
        else
           if (lseg) then
              keymask(2,nend)=ind-1
              nend=nend+1
              lseg=.false. 
           end if
        end if
     end do
     if (lseg) then
        keymask(2,nend)=ind
        nend=nend+1
     end if
  end do

  !check
  if (nend /= nsrt .or. nend /= nseg_loc) then
     write(*,'(1x,a,2(i6))')&
          'ERROR: problem in build_keymask',&
          nend,nsrt,nseg_loc
     stop
  end if

end subroutine build_keymask

subroutine segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc,keymask)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(nseg_loc), intent(out) :: keyv_loc
  integer, dimension(2,nseg_loc), intent(out) :: keyg_loc,keymask
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
        !index of the compressed function
        ind=i-i0+jj
        i1l=i-i1sc
        i2l=i2-i2sc
        i3l=i3-i3sc
        ngridp=i3l*((n1l+1)*(n2l+1)) + i2l*(n1l+1) + i1l+1
        if (go .and. (i1sc <= i .and. i <= i1ec)) then
           nvctr_check=nvctr_check+1
           if (.not. lseg) then
              nsrt=nsrt+1
              keymask(1,nsrt)=ind
              keyg_loc(1,nsrt)=ngridp
              keyv_loc(nsrt)=nvctr_check
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              keymask(2,nend)=ind-1
              keyg_loc(2,nend)=ngridp-1
              lseg=.false. 
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
        keymask(2,nend)=ind
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

subroutine num_segkeys_loc(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
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
        nvctr_check=nvctr_check+1
        if (go .and. (i1sc <= i .and. i <= i1ec)) then
           nvctr_loc=nvctr_loc+1
           if (.not. lseg) then
              nsrt=nsrt+1
           end if
           lseg=.true.
        else
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

!conversion of the global bounds into a given localisation region
subroutine make_bounds_loc(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,bounds,bounds_loc)
  use module_base
  use module_types
  implicit none
  type(convolutions_bounds), intent(in) :: bounds
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec
  type(convolutions_bounds), intent(out) :: bounds_loc
  !local variables
  character(len=*), parameter :: subname='kb_conversion'
  integer :: n1l,n2l,n3l,i_stat
  
  !dimensions of the localisation region (O:nIl)
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  call kb_conversion(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,bounds%kb,bounds_loc%kb)

  call sb_conversion(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,bounds%sb,bounds_loc%sb)

  call gb_conversion(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,bounds%gb,bounds_loc%gb)

  !convert the real bounds array
  allocate(bounds_loc%ibyyzz_r(2,-14:2*n2l+16,-14:2*n3l+16+ndebug),stat=i_stat)
  call memocc(i_stat,bounds_loc%ibyyzz_r,'bounds_loc%ibyyzz_r',subname)
  call bound_conversion(-14,2*n2+16,-14,2*n3+16,&
       0,2*i1sc,2*i1ec+30,-14+2*i2sc,2*i2ec+16,-14+2*i3sc,2*i3ec+16,&
       bounds%ibyyzz_r,bounds_loc%ibyyzz_r)

end subroutine make_bounds_loc

!convert the kinetic bounds
subroutine kb_conversion(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,kb,kb_loc)
  use module_base
  use module_types
  implicit none
  type(kinetic_bounds), intent(in) :: kb
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec
  type(kinetic_bounds), intent(out) :: kb_loc
  !local variables
  character(len=*), parameter :: subname='kb_conversion'
  integer :: n1l,n2l,n3l,i_stat
  
  !dimensions of the localisation region (O:nIl)
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  allocate(kb_loc%ibyz_c(2,0:n2l,0:n3l+ndebug),stat=i_stat)
  call memocc(i_stat,kb_loc%ibyz_c,'kb_loc%ibyz_c',subname)
  allocate(kb_loc%ibyz_f(2,0:n2l,0:n3l+ndebug),stat=i_stat)
  call memocc(i_stat,kb_loc%ibyz_f,'kb_loc%ibyz_f',subname)
  call bound_conversion(0,n2,0,n3,&
       0,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       kb%ibyz_c,kb_loc%ibyz_c)
  call bound_conversion(0,n2,0,n3,&
       0,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
       kb%ibyz_f,kb_loc%ibyz_f)

  allocate(kb_loc%ibxz_c(2,0:n1l,0:n3l+ndebug),stat=i_stat)
  call memocc(i_stat,kb_loc%ibxz_c,'kb_loc%ibxz_c',subname)
  allocate(kb_loc%ibxz_f(2,0:n1l,0:n3l+ndebug),stat=i_stat)
  call memocc(i_stat,kb_loc%ibxz_f,'kb_loc%ibxz_f',subname)
  call bound_conversion(0,n1,0,n3,&
       0,i2sc,i2ec,i1sc,i1ec,i3sc,i3ec,&
       kb%ibxz_c,kb_loc%ibxz_c)
  call bound_conversion(0,n1,0,n3,&
       0,i2sc,i2ec,i1sc,i1ec,i3sc,i3ec,&
       kb%ibxz_f,kb_loc%ibxz_f)

  allocate(kb_loc%ibxy_c(2,0:n1l,0:n2l+ndebug),stat=i_stat)
  call memocc(i_stat,kb_loc%ibxy_c,'kb_loc%ibxy_c',subname)
  allocate(kb_loc%ibxy_f(2,0:n1l,0:n2l+ndebug),stat=i_stat)
  call memocc(i_stat,kb_loc%ibxy_f,'kb_loc%ibxy_f',subname)
  call bound_conversion(0,n1,0,n2,&
       0,i3sc,i3ec,i1sc,i1ec,i2sc,i2ec,&
       kb%ibxy_c,kb_loc%ibxy_c)
  call bound_conversion(0,n1,0,n2,&
       0,i3sc,i3ec,i1sc,i1ec,i2sc,i2ec,&
       kb%ibxy_f,kb_loc%ibxy_f)
  
end subroutine kb_conversion

!convert the shrink bounds
subroutine sb_conversion(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,sb,sb_loc)
  use module_base
  use module_types
  implicit none
  type(shrink_bounds), intent(in) :: sb
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec
  type(shrink_bounds), intent(out) :: sb_loc
  !local variables
  character(len=*), parameter :: subname='sb_conversion'
  integer :: n1l,n2l,n3l,i_stat,nfl1l,nfu1l,nfl2l,nfu2l,nfl3l,nfu3l
  
  !dimensions of the localisation region (O:nIl)
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  !dimensions for the fine localisation region
  nfl1l=max(nfl1,i1sc)
  nfl2l=max(nfl2,i2sc)
  nfl3l=max(nfl3,i3sc)

  nfu1l=min(nfu1,i1ec)
  nfu2l=min(nfu2,i2ec)
  nfu3l=min(nfu3,i3ec)


  allocate(sb_loc%ibzzx_c(2,-14:2*n3l+16,0:n1l+ndebug),stat=i_stat)
  call memocc(i_stat,sb_loc%ibzzx_c,'sb_loc%ibzzx_c',subname)
  call bound_conversion(-14,2*n3+16,0,n1,&
       0,i2sc,i2ec,-14+2*i3sc,2*i3ec+16,i1sc,i1ec,&
       sb%ibzzx_c,sb_loc%ibzzx_c)

  allocate(sb_loc%ibyyzz_c(2,-14:2*n2l+16,-14:2*n3l+16+ndebug),stat=i_stat)
  call memocc(i_stat,sb_loc%ibyyzz_c,'sb_loc%ibyyzz_c',subname)
  call bound_conversion(-14,2*n2+16,-14,2*n3+16,&
       0,i1sc,i1ec,-14+2*i2sc,2*i2ec+16,&
       -14+2*i3sc,2*i3ec+16,sb%ibyyzz_c,sb_loc%ibyyzz_c)

  allocate(sb_loc%ibzzx_f(2,-14+2*nfl3l:2*nfu3l+16,nfl1l:nfu1l+ndebug),stat=i_stat)
  call memocc(i_stat,sb_loc%ibzzx_f,'sb_loc%ibzzx_f',subname)
  call bound_conversion(-14+2*nfl3,2*nfu3+16,nfl1,nfu1,&
       nfl2l,nfl2l,nfu2l,-14+2*nfl3l,2*nfu3l+16,&
       nfl1l,nfu1l,sb%ibzzx_f,sb_loc%ibzzx_f)

  allocate(sb_loc%ibyyzz_f(2,-14+2*nfl2l:2*nfu2l+16,-14+2*nfl3l:2*nfu3l+16+ndebug),stat=i_stat)
  call memocc(i_stat,sb_loc%ibyyzz_f,'sb_loc%ibyyzz_f',subname)
  call bound_conversion(-14+2*nfl2,2*nfu2+16,-14+2*nfl3,2*nfu3+16,&
       nfl1l,nfl1l,nfu1l,-14+2*nfl2l,2*nfu2l+16,-14+2*nfl3l,2*nfu3l+16,&
       sb%ibyyzz_c,sb_loc%ibyyzz_c)

  allocate(sb_loc%ibxy_ff(2,nfl1l:nfu1l,nfl2l:nfu2l+ndebug),stat=i_stat)
  call memocc(i_stat,sb_loc%ibxy_ff,'sb_loc%ibxy_ff',subname)
  call bound_conversion(nfl1,nfu1,nfl2,nfu2,&
       nfl3l,nfl3l,nfu3l,nfl1l,nfu1l,nfl2l,nfu2l,&
       sb%ibyyzz_c,sb_loc%ibyyzz_c)
 
end subroutine sb_conversion


!convert the grow bounds
subroutine gb_conversion(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,gb,gb_loc)
  use module_base
  use module_types
  implicit none
  type(grow_bounds), intent(in) :: gb
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec
  type(grow_bounds), intent(out) :: gb_loc
  !local variables
  character(len=*), parameter :: subname='gb_conversion'
  integer :: n1l,n2l,n3l,i_stat,nfl1l,nfu1l,nfl2l,nfu2l,nfl3l,nfu3l
  
  !dimensions of the localisation region (O:nIl)
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  !dimensions for the fine localisation region
  nfl1l=max(nfl1,i1sc)
  nfl2l=max(nfl2,i2sc)
  nfl3l=max(nfl3,i3sc)

  nfu1l=min(nfu1,i1ec)
  nfu2l=min(nfu2,i2ec)
  nfu3l=min(nfu3,i3ec)

  allocate(gb_loc%ibzxx_c(2,0:n3l,-14:2*n1l+16+ndebug),stat=i_stat)
  call memocc(i_stat,gb_loc%ibzxx_c,'gb_loc%ibzxx_c',subname)
  call bound_conversion(0,n3,-14,2*n1+16,&
       0,i2sc,i2ec,i3sc,i3ec,-14+2*i1sc,2*i1ec+16,&
       gb%ibzxx_c,gb_loc%ibzxx_c)

  allocate(gb_loc%ibxxyy_c(2,-14:2*n1l+16,-14:2*n2l+16+ndebug),stat=i_stat)
  call memocc(i_stat,gb_loc%ibxxyy_c,'gb_loc%ibxxyy_c',subname)
  call bound_conversion(-14,2*n1+16,-14,2*n2+16,&
       0,i3sc,i3ec,-14+2*i1sc,2*i1ec+16,-14+2*i2sc,2*i2ec+16,&
       gb%ibxxyy_c,gb_loc%ibxxyy_c)

  allocate(gb_loc%ibyz_ff(2,nfl2l:nfu2l,nfl3l:nfu3l+ndebug),stat=i_stat)
  call memocc(i_stat,gb_loc%ibyz_ff,'gb_loc%ibyz_ff',subname)
  call bound_conversion(nfl2,nfu2,nfl3,nfu3,&
       nfl1l,nfl1l,nfu1l,nfl2l,nfu2l,nfl3l,nfu3l,&
       gb%ibyz_ff,gb_loc%ibyz_ff)

  allocate(gb_loc%ibzxx_f(2,nfl3l:nfu3l,2*nfl1l-14:2*nfu1l+16+ndebug),stat=i_stat)
  call memocc(i_stat,gb_loc%ibzxx_f,'gb_loc%ibzxx_f',subname)
  call bound_conversion(nfl3,nfu3,2*nfl1-14,2*nfu1+16,&
       nfl2l,nfl2l,nfu2l,nfl3l,nfu3l,2*nfl1l-14,2*nfu1l+16,&
       gb%ibzxx_f,gb_loc%ibzxx_f)

  allocate(gb_loc%ibxxyy_f(2,2*nfl1l-14:2*nfu1l+16,2*nfl2l-14:2*nfu2l+16+ndebug),stat=i_stat)
  call memocc(i_stat,gb_loc%ibxxyy_f,'gb_loc%ibxxyy_f',subname)
  call bound_conversion(2*nfl1-14,2*nfu1+16,2*nfl2-14,2*nfu2+16,&
       nfl3l,nfl3l,nfu3l,2*nfl1l-14,2*nfu1l+16,2*nfl2l-14,2*nfu2l+16,&
       gb%ibxxyy_f,gb_loc%ibxxyy_f)
  
end subroutine gb_conversion


!with this subroutine we should convert a bound array from its global version to 
!the version which is compatible for a given localisation region
!! INPUTS
!!   nl2,nu2,nl3,nu3  limits of the global bounds array in the orthogonal directions
!!   nl1_loc          lower limit of the local bounds array in the interested direction
!!   i1s,i1e,i2s,i2e
!!   i3s,i3e          limits of the localisation region in the global system coordinates
!!   ib               global bounds array
!! OUTPUT
!!   ib_loc           local bounds array
subroutine bound_conversion(nl2,nu2,nl3,nu3,nl1_loc,i1s,i1e,i2s,i2e,i3s,i3e,ib,ib_loc)
  implicit none
  integer, intent(in) :: nl1_loc,nl2,nu2,nl3,nu3,i1s,i1e,i2s,i2e,i3s,i3e
  integer, dimension(2,nl2:nu2,nl3:nu3), intent(in) :: ib
  integer, dimension(2,i2s:i2e,i3s:i3e), intent(out) :: ib_loc
  !local variables
  integer :: j2,j3,limlow,limup

  !control whether the localisation region interval is inside the allowed dimensions
  if (i2s>nu2 .or. i2e < nl2 .or. i3s>nu3 .or. i3e < nl3) then
     do j3=i3s,i3e
        do j2=i2s,i2e
           ib_loc(1,j2,j3)=1000
           ib_loc(2,j2,j3)=-1000
        end do
     end do
  else
     do j3=i3s,i3e
        do j2=i2s,i2e
           limlow=max(ib(1,j2,j3),i1s)
           limup=min(ib(2,j2,j3),i1e)
           if (limlow < limup) then
              ib_loc(1,j2,j3)=limlow-i1s+nl1_loc
              ib_loc(2,j2,j3)=limup-i1s+nl1_loc
           else
              ib_loc(1,j2,j3)=1000
              ib_loc(2,j2,j3)=-1000
           end if
        end do
     end do
  end if

end subroutine bound_conversion
