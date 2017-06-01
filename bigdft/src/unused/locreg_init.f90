!> Pass from the global wavefunction for a set of orbitals to the local wavefunctions
!! for all the orbitals
!! INPUTS
!! @param nlr           number of localisation regions contained by the compressed form
!! @param nvctr_c,f     compressed elements of the global function
!! @param nseglr        array of compressed elements for each localisation region
!!                      number of segments of keymask array
!! @param nvctr_tot     number of components of the distributed wavefunction
!! @param nseg_tot      total number of segments for the mask array
!! @param keymask       array of masking element for each localisation region           
!! @param psi           global wavefunctions with distributed orbitals
!! OUTPUTS
!! @param psiloc        wavefunctions of all the localisation regions 
!! @param ncount        array of elements to be sent to each processor for MPI_ALLTOALLV procedure
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
END SUBROUTINE rhoswitch_waves



subroutine segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,keyg,keyv,& !n(c) nvctr (arg:11)
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc)!,keymask)
  implicit none
  integer, intent(in) :: n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nseg_loc,nvctr_loc !n(c) nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(nseg_loc), intent(out) :: keyv_loc
  integer, dimension(2,nseg_loc), intent(out) :: keyg_loc!,keymask
  !local variables
  logical :: go,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend,nvctr_check,n1l,n2l,i1l,i2l,i3l !n(c) n3l
  integer :: ngridp

  !dimensions of the localisation region (O:nIl)
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  !n(c) n3l=i3ec-i3sc

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
              !keymask(1,nsrt)=ind
              keyg_loc(1,nsrt)=ngridp
              keyv_loc(nsrt)=nvctr_check
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              !keymask(2,nend)=ind-1
              keyg_loc(2,nend)=ngridp-1
              lseg=.false. 
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
        !keymask(2,nend)=ind
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

END SUBROUTINE segkeys_loc



subroutine num_segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc)
  implicit none
  integer, intent(in) :: n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
  !local variables
  logical :: go,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check

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
    
END SUBROUTINE num_segkeys_loc


!> Conversion of the global bounds into a given localisation region
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
  integer :: n2l,n3l,i_stat !n(c) n1l
  
  !dimensions of the localisation region (O:nIl)
  !n(c) n1l=i1ec-i1sc
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

END SUBROUTINE make_bounds_loc


!> Convert the kinetic bounds
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
  
END SUBROUTINE kb_conversion



!> Convert the shrink bounds
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
 
END SUBROUTINE sb_conversion




!> Convert the grow bounds
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
  
END SUBROUTINE gb_conversion




!> With this subroutine we should convert a bound array from its global version to 
!! the version which is compatible for a given localisation region
!! INPUTS
!!   @param nl2,nu2,nl3,nu3  limits of the global bounds array in the orthogonal directions
!!   @param nl1_loc          lower limit of the local bounds array in the interested direction
!!   @param i1s,i1e,i2s,i2e
!!   @param i3s,i3e          limits of the localisation region in the global system coordinates
!!   @param ib               global bounds array
!! OUTPUT
!!   @param ib_loc           local bounds array
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

END SUBROUTINE bound_conversion

!>   This subroutine define other wavefunctions descriptors starting from the original descriptors 
!!   and the limits of a given localisation region
!!   it also returns an array which is used to mask the compressed wavefunction into the new one
!! INPUTS
!!   @param ilocreg        localisation region to be considered
!!   @param nlocreg        total number of localisation regions
!!   @param n1,n2          original dimensions of the global box
!!   @param lrlims         array of limits of the localisation regions (global system coordinates)
!!   @param wfdg           global wavefunction descriptors structure
!! OUTPUT
!!   @param wfdl           local wavefunction descriptors structure in local system coordinates
!!   @param keymask        mask array for traducing the wavefunction in compressed form
!!                         to the wavefunction in compressed form for the local system
!!   @param ncountlocreg   array of elements for each localisation region
!!subroutine loc_wfd(ilocreg,nlocreg,n1,n2,lrlims,wfdg,wfdl,keymask,ncountlocreg)
!!  use module_base
!!  use module_types
!!  implicit none
!!  type(wavefunctions_descriptors), intent(in) :: wfdg
!!  integer, intent(in) :: ilocreg,nlocreg,n1,n2
!!  integer, dimension(2,3,nlocreg), intent(in) :: lrlims
!!  type(wavefunctions_descriptors), intent(out) :: wfdl
!!  integer, dimension(nlocreg), intent(out) :: ncountlocreg
!!  integer, dimension(:), pointer :: keymask
!!  !local variables
!!  character(len=*), parameter :: subname='loc_wfd'
!!  integer :: i_stat
!!  integer :: iloc,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nvctr_c,nseg_c,nseg_f,nvctr_f,ndimkey
!!
!!  !calculate the number of segments of the new descriptors for each localisation region
!!  !and the dimension of the array for the translation of the localisation regions
!!  ndimkey=0
!!  do iloc=1,nlocreg
!!     if (iloc /= ilocreg) then
!!        !coarse part
!!        call num_segkeys_loc(n1,n2,lrlims(1,1,iloc),lrlims(2,1,iloc),&
!!             lrlims(1,2,iloc),lrlims(2,2,iloc),lrlims(1,3,iloc),lrlims(2,3,iloc),&
!!             wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
!!             nseg_c,nvctr_c)
!!        !fine part
!!        call num_segkeys_loc(n1,n2,lrlims(1,1,iloc),lrlims(2,1,iloc),&
!!             lrlims(1,2,iloc),lrlims(2,2,iloc),lrlims(1,3,iloc),lrlims(2,3,iloc),&
!!             wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!             wfdg%keyv(wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!             nseg_f,nvctr_f)
!!        ncountlocreg(iloc)=nvctr_c+7*nvctr_f
!!        ndimkey=ndimkey+nseg_c+nseg_f
!!     end if
!!  end do
!!
!!  i1sc=lrlims(1,1,ilocreg)
!!  i1ec=lrlims(2,1,ilocreg)
!!  i2sc=lrlims(1,2,ilocreg)
!!  i2ec=lrlims(2,2,ilocreg)
!!  i3sc=lrlims(1,3,ilocreg)
!!  i3ec=lrlims(2,3,ilocreg)
!!
!!  !coarse part
!!  call num_segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
!!       wfdg%nseg_c,wfdg%nvctr_c,wfdg%keyg(1,1),wfdg%keyv(1),&
!!       wfdl%nseg_c,wfdl%nvctr_c)
!!
!!  !fine part
!!  call num_segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,&
!!       wfdg%nseg_f,wfdg%nvctr_f,wfdg%keyg(1,wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdg%keyv(wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdl%nseg_f,wfdl%nvctr_f)
!!
!!  ncountlocreg(ilocreg)=wfdl%nvctr_c+7*wfdl%nvctr_f
!!  ndimkey=ndimkey+wfdl%nseg_c+wfdl%nseg_f
!!
!!  call allocate_wfd(wfdl,subname)
!!
!!  allocate(keymask(ndimkey+ndebug),stat=i_stat)
!!  call memocc(i_stat,keymask,'keymask',subname)
!!
!!  !now fill the local wavefunction descriptors
!!  !and define the mask array for the wavefunction
!!  !coarse part
!!  call segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,& !n(m)
!!       wfdg%nseg_c,wfdg%keyg(1,1),wfdg%keyv(1),&
!!       wfdl%nseg_c,wfdl%nvctr_c,wfdl%keyg(1,1),wfdl%keyv(1))!,keymask(1))
!!
!!  !fine part
!!  call segkeys_loc(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,& !n(m)
!!       wfdg%nseg_f,wfdg%keyg(1,wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdg%keyv(wfdg%nseg_c+min(1,wfdg%nseg_f)),&
!!       wfdl%nseg_f,wfdl%nvctr_f,wfdl%keyg(1,wfdl%nseg_c+min(1,wfdl%nseg_f)),&
!!       wfdl%keyv(wfdl%nseg_c+min(1,wfdl%nseg_f)))!,&
!!       !keymask(wfdg%nseg_c+1))
!!
!!  !a little check on the masking array
!!!!!  if (count(maskarr) /= wfdl%nvctr_c+7*wfdl%nvctr_f) then
!!!!!     write(*,'(1x,a)')'ERROR : Masking problem, check maskarr'
!!!!!     stop
!!!!!  end if
!!
!!END SUBROUTINE loc_wfd


!!subroutine build_keymask(n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg_tot,keyg,keyv,&
!!     nseg_loc,keymask)
!!  implicit none
!!  integer, intent(in) :: n1,n2,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg_tot,nseg_loc
!!  integer, dimension(nseg_tot), intent(in) :: keyv
!!  integer, dimension(2,nseg_tot), intent(in) :: keyg
!!  integer, dimension(2,nseg_loc), intent(out) :: keymask
!!  !local variables
!!  logical :: go,lseg
!!  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend
!!
!!  !start and end points
!!  nsrt=0
!!  nend=0
!!  do iseg=1,nseg_tot
!!     jj=keyv(iseg)
!!     j0=keyg(1,iseg)
!!     j1=keyg(2,iseg)
!!     ii=j0-1
!!     i3=ii/((n1+1)*(n2+1))
!!     ii=ii-i3*(n1+1)*(n2+1)
!!     i2=ii/(n1+1)
!!     i0=ii-i2*(n1+1)
!!     i1=i0+j1-j0
!!     go=(i3sc <= i3 .and. i3 <= i3ec) .and. (i2sc <= i2 .and. i2 <= i2ec)
!!     lseg=.false.
!!     do i=i0,i1
!!        !index of the compressed function
!!        ind=i-i0+jj
!!        if (go .and. (i1sc <= i .and. i <= i1ec)) then
!!           if (.not. lseg) then
!!              nsrt=nsrt+1
!!              keymask(1,nsrt)=ind
!!           end if
!!           lseg=.true.
!!        else
!!           if (lseg) then
!!              keymask(2,nend)=ind-1
!!              nend=nend+1
!!              lseg=.false. 
!!           end if
!!        end if
!!     end do
!!     if (lseg) then
!!        keymask(2,nend)=ind
!!        nend=nend+1
!!     end if
!!  end do
!!
!!  !check
!!  if (nend /= nsrt .or. nend /= nseg_loc) then
!!     write(*,'(1x,a,2(i6))')&
!!          'ERROR: problem in build_keymask',&
!!          nend,nsrt,nseg_loc
!!     stop
!!  end if
!!
!!END SUBROUTINE build_keymask
