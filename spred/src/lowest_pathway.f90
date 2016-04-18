module previous_graph
implicit none
type pvari_arr
    integer, allocatable :: previous(:) 
end type
type padjacency_list
    integer :: nnodes
    integer, allocatable :: nprevious(:) !stores the number of neighbours for each minimum
    type(pvari_arr), allocatable :: mini(:)
end type
contains
subroutine init_padjacency_list(ajl,nmin)
    implicit none
    !params
    type(padjacency_list), intent(inout) :: ajl
    integer, intent(in) :: nmin
    !internal
    integer :: nPrevPerMin
    integer :: i

    nPrevPerMin=5
    allocate(ajl%mini(nmin),ajl%nprevious(nmin))
    ajl%nprevious=0

    do i=1,nmin
        allocate(ajl%mini(i)%previous(nPrevPerMin))
        ajl%mini(i)%previous=-1
    enddo

end subroutine
integer function get_prev(ajl,i,iprev)
    implicit none
    !params
    type(padjacency_list), intent(in) :: ajl
    integer, intent(in) :: i,iprev
    !internal
    integer :: nprevOfI

    nprevOfI=ajl%nprevious(i)
    if(iprev>nprevOfI)stop 'j>nprevOfI in function get_prev'

    get_prev=ajl%mini(i)%previous(iprev)

end function
subroutine add_prev(ajl,i,iprev)
    implicit none
   !params
    type(padjacency_list), intent(inout) :: ajl 
    integer, intent(in) :: i,iprev
    !internal

    call pinsert(ajl%nprevious(i),ajl%mini(i)%previous,iprev)

    
end subroutine
subroutine pinsert(nprevious,previous,iprev)
    implicit none
    !params
    integer, intent(inout) :: nprevious
    integer, intent(inout),allocatable :: previous(:)
    integer, intent(in) :: iprev
    !internal
    integer, parameter :: resize=5
    integer :: k,nprevmax
    integer, allocatable :: previousTmp(:)

    nprevious=nprevious+1
!write(*,*)nneighbors
    nprevmax=size(previous)
    if(nprevious>nprevmax)then!resize
        allocate(previousTmp(nprevmax))
!write(*,*)'resize',nneighbors
        previousTmp=previous
        deallocate(previous)
        allocate(previous(nprevmax+resize))
        do k=1,nprevmax
            previous(k)=previousTmp(k)
        enddo
        deallocate(previousTmp)
    endif
    previous(nprevious)=iprev
    
end subroutine
end module

module hoppcount_adjacency_list
implicit none
type ivari_arr
    integer, allocatable :: neighbors(:) !stores the id of the neighbor
end type
type iadjacency_list
    integer :: nnodes
    integer, allocatable :: nneighbors(:) !stores the number of neighbours for each minimum
    type(ivari_arr), allocatable :: conn(:)
end type
contains
subroutine init_iadjacency_list(ajl,nmin)
    implicit none
    !params
    type(iadjacency_list), intent(inout) :: ajl
    integer, intent(in) :: nmin
    !internal
    integer :: nneighPerMin
    integer :: i

    nneighPerMin=10
    allocate(ajl%conn(nmin),ajl%nneighbors(nmin))
    ajl%nneighbors=0

    do i=1,nmin
        allocate(ajl%conn(i)%neighbors(nneighPerMin))
        ajl%conn(i)%neighbors=-1
    enddo

end subroutine
subroutine dealloc_iadjacency_list(ajl,nmin)
    implicit none
    !params
    type(iadjacency_list), intent(inout) :: ajl
    integer, intent(in) :: nmin
    !internal
    integer :: i

    do i=1,nmin
        deallocate(ajl%conn(i)%neighbors)
    enddo
    deallocate(ajl%conn,ajl%nneighbors)


end subroutine
integer function get_dist(ajl,i,j)
    implicit none
    !params
    type(iadjacency_list), intent(in) :: ajl
    integer, intent(in) :: i,j
    !internal
    integer :: nneighbOfI,k

    nneighbOfI=ajl%nneighbors(i)
    do k=1,nneighbOfI
        if(ajl%conn(i)%neighbors(k)==j)then
            get_dist=1
            return
        endif
    enddo   
    get_dist=huge(1)
end function
subroutine add_connection(ajl,i,j)
    implicit none
   !params
    type(iadjacency_list), intent(inout) :: ajl 
    integer, intent(in) :: i,j 
    !internal

    call insrt(ajl%nneighbors(i),ajl%conn(i)%neighbors,j)
    call insrt(ajl%nneighbors(j),ajl%conn(j)%neighbors,i)

    
end subroutine
subroutine insrt(nneighbors,neighbors,neighbor)
    implicit none
    !params
    integer, intent(inout) :: nneighbors
    integer, intent(inout),allocatable :: neighbors(:)
    integer, intent(in) :: neighbor
    !internal
    integer, parameter :: resize=10
    integer :: nneighmax,k
    integer, allocatable :: neighborsTmp(:)

    do k=1,nneighbors
        if(neighbors(k)==neighbor)then
            return
        endif
    enddo

    nneighbors=nneighbors+1
!write(*,*)nneighbors
    nneighmax=size(neighbors)
    if(nneighbors>nneighmax)then!resize
        allocate(neighborsTmp(nneighmax))
!write(*,*)'resize',nneighbors
        neighborsTmp=neighbors
        deallocate(neighbors)
        allocate(neighbors(nneighmax+resize))
        do k=1,nneighmax
            neighbors(k)=neighborsTmp(k)
        enddo
        deallocate(neighborsTmp)
    endif
    neighbors(nneighbors)=neighbor
    
end subroutine
end module

module barrier_adjacency_list
implicit none
type vari_arr
    integer, allocatable :: neighbors(:) !stores the id of the neighbor
    integer, allocatable :: barrier_id(:)
    real(8), allocatable :: barriers(:)  !stores the energy of minimum id
end type
type adjacency_list
    integer :: nnodes
    integer, allocatable :: nneighbors(:) !stores the number of neighbours for each minimum
    type(vari_arr), allocatable :: conn(:)
end type
contains
subroutine init_adjacency_list(ajl,nmin)
    implicit none
    !params
    type(adjacency_list), intent(inout) :: ajl
    integer, intent(in) :: nmin
    !internal
    integer :: nneighPerMin
    integer :: i

    nneighPerMin=10
    allocate(ajl%conn(nmin),ajl%nneighbors(nmin))
    ajl%nneighbors=0

    do i=1,nmin
        allocate(ajl%conn(i)%neighbors(nneighPerMin),ajl%conn(i)%barriers(nneighPerMin),&
                 ajl%conn(i)%barrier_id(nneighPerMin))
        ajl%conn(i)%neighbors=-1
        ajl%conn(i)%barriers=huge(1.d0)
        ajl%conn(i)%barrier_id=-1
    enddo

end subroutine
subroutine dealloc_adjacency_list(ajl,nmin)
    implicit none
    !params
    type(adjacency_list), intent(inout) :: ajl
    integer, intent(in) :: nmin
    !internal
    integer :: i

    do i=1,nmin
        deallocate(ajl%conn(i)%neighbors,ajl%conn(i)%barriers)
    enddo
    deallocate(ajl%conn,ajl%nneighbors)


end subroutine

subroutine get_barrier(ajl,i,j,barrener,barrid)
    implicit none
    !params
    type(adjacency_list), intent(in) :: ajl
    integer, intent(in) :: i,j
    real(8), intent(out) :: barrener
    integer, intent(out) :: barrid
    !internal
    integer :: nneighbOfI,k

    nneighbOfI=ajl%nneighbors(i)
    do k=1,nneighbOfI
        if(ajl%conn(i)%neighbors(k)==j)then
            barrener=ajl%conn(i)%barriers(k)
            barrid=ajl%conn(i)%barrier_id(k)
            return
        endif
    enddo   
    barrener=huge(1.d0)
    barrid=-1
end subroutine
subroutine add_barrier(ajl,i,j,barr_id,ebarr)
    implicit none
   !params
    type(adjacency_list), intent(inout) :: ajl 
    integer, intent(in) :: i,j,barr_id 
    real(8), intent(in) :: ebarr
    !internal

    call insert(ajl%nneighbors(i),ajl%conn(i)%neighbors,&
         ajl%conn(i)%barriers,ajl%conn(i)%barrier_id,barr_id,j,ebarr)
    call insert(ajl%nneighbors(j),ajl%conn(j)%neighbors,&
         ajl%conn(j)%barriers,ajl%conn(j)%barrier_id,barr_id,i,ebarr)

    
end subroutine
subroutine insert(nneighbors,neighbors,barriers,barrier_id,barr_id,neighbor,ebarr)
    implicit none
    !params
    integer, intent(inout) :: nneighbors
    integer, intent(inout),allocatable :: neighbors(:)
    integer, intent(inout),allocatable :: barrier_id(:)
    real(8), intent(inout),allocatable :: barriers(:)
    integer, intent(in) :: neighbor, barr_id
    real(8), intent(in) :: ebarr
    !internal
    integer, parameter :: resize=10
    integer :: nneighmax,k
    integer, allocatable :: neighborsTmp(:)
    integer, allocatable :: barrier_idTmp(:)
    real(8), allocatable :: barriersTmp(:)

    do k=1,nneighbors
        if(neighbors(k)==neighbor)then
   !         barriers(k)=min(ebarr,barriers(k))
             if(ebarr<barriers(k))then
                barriers(k)=ebarr
                barrier_id(k) = barr_id
             endif
            return
        endif
    enddo

    nneighbors=nneighbors+1
!write(*,*)nneighbors
    nneighmax=size(neighbors)
    if(nneighbors>nneighmax)then!resize
        allocate(barriersTmp(nneighmax),neighborsTmp(nneighmax),&
                 barrier_idTmp(nneighmax))
!write(*,*)'resize',nneighbors
        barriersTmp=barriers
        barrier_idTmp=barrier_id
        neighborsTmp=neighbors
        deallocate(barriers,neighbors,barrier_id)
        allocate(barriers(nneighmax+resize),neighbors(nneighmax+resize),&
                 barrier_id(nneighmax+resize))
        do k=1,nneighmax
            barriers(k)=barriersTmp(k)
            barrier_id(k)=barrier_idTmp(k)
            neighbors(k)=neighborsTmp(k)
        enddo
        deallocate(barriersTmp,neighborsTmp,barrier_idTmp)
    endif
    barriers(nneighbors)=ebarr
    barrier_id(nneighbors)=barr_id
    neighbors(nneighbors)=neighbor
    
end subroutine
end module


module dij
contains
subroutine dijkstra_pruned(nminin,nts,ts,neig,threshIn,istartIn,idestIn,najl,ajl,idx,idx2,mbener,istat,nspaths,pajl,ntsinpath)
use barrier_adjacency_list
use hoppcount_adjacency_list
use previous_graph
    implicit none
    !parameter
    integer, intent(in) :: nts,nminin,neig(2,nts),istartIn,idestIn
    real(8), intent(in) :: ts(nts),threshIn
    real(8), intent(out) :: mbener
    type(adjacency_list), intent(out) :: ajl
    type(padjacency_list), intent(out) :: pajl
    integer, intent(out) :: istat,nspaths,ntsinpath,najl
    integer, allocatable, intent(out) :: idx(:)
    !internal
    integer :: nmin,nc,i,left,right,istart,idest
    real(8), allocatable :: maxbarr(:)
    integer, allocatable :: mind(:),count(:),idx2(:)
    type(iadjacency_list) :: iajl
    real(8) :: thresh
    integer :: icurr,itmp
    !functions
!    integer :: indexof


    istat=0
    nmin=nminin
    istart=istartIn
    idest=idestIn
    thresh=threshIn

    !trim the tree according to user defined threshold
    !this threshold should be made by an educated guess
    !can significantly speed up the computation (several order of magnitudes!)
    allocate(idx(nmin))
    idx=huge(1)
write(*,*)'nts',nts
    do i=1,nts
        if(ts(i)<=thresh)then
            idx(neig(1,i)) = neig(1,i)
            idx(neig(2,i)) = neig(2,i)
        endif
    enddo
    nc=0
    do i=1,nmin
        if(idx(i)<huge(1))then
            nc=nc+1
            idx(i)=nc
        endif
    enddo
    nmin=nc
    najl=nmin
write(*,*)'nmin',nmin
    call init_adjacency_list(ajl,najl)
    allocate(maxbarr(nmin))
    do i=1,nts
        if(ts(i)<=thresh)then
        left=neig(1,i)
        right=neig(2,i)
        call add_barrier(ajl,idx(left),idx(right),i,ts(i))
        endif
    enddo
    istart=idx(istart)
    idest=idx(idest)
    if(istart==huge(1).or.idest==huge(1))then
        istat = -1
!        write(*,*)trim(adjustl(ci1))//' and '//trim(adjustl(ci2))//' are not connected by a path with barriers lower than thresh'
!        stop
        mbener=huge(1.d0)
        return
    endif

  call dijkstra_barr( najl, ajl, istart,idest,maxbarr,istat )
  mbener = maxbarr(idest)
!  call dealloc_adjacency_list(ajl,najl)

    !ok, we have found the lowest energy barrier between istart and idest
    !now trim the tree again, using mbener as threshold in order
    !to find the shortest (least number of transition states) path
    nmin=nminin
    istart=istartIn
    idest=idestIn
    thresh=mbener
    allocate(idx2(nmin))
    idx2=huge(1)
    do i=1,nts
        if(ts(i)<=thresh)then
            idx2(neig(1,i)) = neig(1,i)
            idx2(neig(2,i)) = neig(2,i)
        endif
    enddo
    nc=0
    do i=1,nmin
        if(idx2(i)<huge(1))then
            nc=nc+1
            idx2(i)=nc
        endif
    enddo
    nmin=nc
    deallocate(maxbarr)
    allocate(mind(nmin),count(nmin))
    call init_iadjacency_list(iajl,nmin)
    call init_padjacency_list(pajl,nmin)
    do i=1,nts
        if(ts(i)<=thresh)then
        left=neig(1,i)
        right=neig(2,i)
        call add_connection(iajl,idx2(left),idx2(right))
        endif
    enddo
    istart=idx2(istart)
    idest=idx2(idest)
    if(istart==huge(1).or.idest==huge(1))then
        istat = -1
!        write(*,*)trim(adjustl(ci1))//' and '//trim(adjustl(ci2))//' are not connected by a path with barriers lower than thresh'
!        stop
        return
    endif


  call dijkstra_shortest( nmin, iajl, istart,idest,mind,pajl,count,istat )
nspaths=count(idest)
ntsinpath=mind(idest)
!write(*,*)'Number of shortest paths: ',count(idest)
!write(*,*)'Number of transition states in shortest path: ',mind(idest)
    
!!    icurr=idest
!!    i=0
!!    do
!!        itmp=indexof(nminin,idx2,icurr)
!!        i=i+1
!!        path(i)=itmp
!!!        write(99,*)icurr,itmp
!!!        write(99,*)icurr
!!        if(itmp==istartIn)return
!!            icurr=previous(icurr)
!!    enddo
end subroutine
subroutine dijkstra_barr( nv, ajl,istart, idest,  maxbarr ,istat)
  use barrier_adjacency_list
  implicit none
    !parameters
  integer, intent(in) :: nv,istart,idest
  type(adjacency_list), intent(in) :: ajl
  real(8), intent(out) :: maxbarr(nv)
  integer, intent(out) :: istat
    !internal
  logical :: connected(nv)
  integer :: i,mv,step
  real(8) :: md
!
!  Start out with only node 1 connected to the tree.
!
  connected = .false.
!  connected(istart) = .true.
!
!  Initialize the minimum distance to the one-step distance.
!
!  maxbarr(1:nv) = ohd(1:nv,istart)
!  maxbarr(1:nv) = ohd(istart,1:nv)
    maxbarr=huge(1.d0)
    maxbarr(istart)=-huge(1.d0)
!
!  Attach one more node on each iteration.
!
  do step = 2,nv
if(mod(step,1000)==0)then
write(*,FMT='(A1,A,F6.2,A,$)')achar(13),' Percent completed: ', (real(step)/real(nv))*100.0,'%'
endif
    call find_lowest(nv,maxbarr,connected,md,mv)
    if( mv == - 1 ) then
!      write( *, '(a)' )  ' '
!      write( *, '(a)' )  'DIJKSTRA_DISTANCE - Warning!'
    write(*,*)
      write( *, '(a,i0.0)' )  '  Search terminated early.',step
!      write( *, '(a)' )  '  Graph might not be connected.'
      istat=-1
      return
    end if
    connected(mv)=.true.
!    if(mv==idest)return
    call update_maxbarr(nv,connected,ajl,mv,maxbarr)
  enddo
write(*,*)
end subroutine
subroutine find_lowest( nv, maxbarr, connected, d, v)
  implicit none
    !parameter
    integer, intent(in) :: nv
    real(8), intent(in) :: maxbarr(nv)
    logical, intent(in) :: connected(nv)
    real(8), intent(out) :: d
    integer, intent(out) :: v
    !internal
    integer :: i
    real(8) :: tmp

  d = huge(1.d0)
  v = -1
  do i = 1, nv
    if( .not. connected(i) .and. maxbarr(i) < d ) then
      d = maxbarr(i)
      v = i
    end if
  end do
end subroutine
subroutine update_maxbarr( nv, connected, ajl, mv, maxbarr )
  use barrier_adjacency_list
  implicit none
    !parameters
    integer, intent(in) :: nv,mv
    logical, intent(in) :: connected(nv)
    type(adjacency_list),intent(in) :: ajl
    real(8), intent(inout) :: maxbarr(nv)
    !internal
    integer :: i,neighb,barrid
    real(8) :: hg,tmp
    hg=huge(1.d0)

  !do i = 1, nv
  do i = 1, ajl%nneighbors(mv)
    neighb=ajl%conn(mv)%neighbors(i)
    if( .not. connected(neighb) ) then
      call get_barrier(ajl,mv,neighb,tmp,barrid)!get_barrier(ajl,i,mv)!ohd(i,mv)
      if( tmp < hg ) then
          !maxbarr(i) = min(maxbarr(i),max(maxbarr(mv),tmp))
          if(maxbarr(mv)>tmp)then
            if(maxbarr(neighb)>maxbarr(mv))then
                maxbarr(neighb)=maxbarr(mv)
            endif
          else
            if(maxbarr(neighb)>tmp)then
                maxbarr(neighb)=tmp
            endif
          endif
      else
stop'Serious problem'
      end if
    end if
  end do

  return
end subroutine

integer function indexof(n,iarr,i)
    implicit none
    !params
    integer :: iarr(n),i,n
    !internal
    integer :: j

    indexof=-1
    do j=1,n
        if(iarr(j)==i)then
            indexof =j
            return
        endif
    enddo
end function
subroutine dijkstra_shortest( nv, iajl,istart, idest,  mind ,pajl,count,istat)
use hoppcount_adjacency_list
use previous_graph
  implicit none
    !parameters
  integer, intent(in) :: nv,istart,idest
  type(iadjacency_list), intent(in) :: iajl
  type(padjacency_list), intent(inout) :: pajl
  integer, intent(out) :: istat,mind(nv),count(nv)
    !internal
  logical :: connected(nv),l
  integer :: mv,step,d
  real(8) :: md
!
!  Start out with only node 1 connected to the tree.
!
  connected = .false.
!  connected(istart) = .true.
!
!  Initialize the minimum distance to the one-step distance.
!
   mind=huge(1)
   mind(istart)=0
   count=0
   count(istart)=1

!
!  Attach one more node on each iteration.
!
  do step = 2,nv
if(mod(step,1000)==0)then
write(*,FMT='(A1,A,F6.2,A,$)')achar(13),' Percent completed: ', (real(step)/real(nv))*100.0,'%'
endif
    call find_closest(nv,mind,connected,d,mv)
    if( mv == - 1 ) then
!      write( *, '(a)' )  ' '
!      write( *, '(a)' )  'DIJKSTRA_DISTANCE - Warning!'
    write(*,*)
      write( *, '(a,i0.0)' )  '  Search terminated early.',step
!      write( *, '(a)' )  '  Graph might not be connected.'
      istat=-1
      return
    end if
    connected(mv)=.true.
!    if(mv==idest)return
    call update_mind(nv,connected,iajl,mv,mind,pajl,count)
  enddo
end subroutine
subroutine find_closest( nv, mind, connected, d, v )
  implicit none
    !parameter
    integer, intent(in) :: nv
    integer, intent(in) :: mind(nv)
    logical, intent(in) :: connected(nv)
    integer, intent(out) :: d
    integer, intent(out) :: v
    !internal
    integer :: i

  d = huge(1)
  v = -1
  do i = 1, nv
    if( .not. connected(i) .and. mind(i) < d ) then
      d = mind(i)
      v = i
    end if
  end do
end subroutine
subroutine update_mind( nv, connected, iajl, mv, mind ,pajl,count)
use hoppcount_adjacency_list
use previous_graph
  implicit none
    !parameters
    integer, intent(in) :: nv,mv
    logical, intent(in) :: connected(nv)
    type(iadjacency_list), intent(in) :: iajl
    type(padjacency_list), intent(inout) :: pajl
    integer, intent(inout) :: mind(nv),count(nv)
    !internal
    integer :: i,hg,tmp,neighb
    hg=huge(1)

! do i = 1, nv
  do i = 1, iajl%nneighbors(mv)
      neighb=iajl%conn(mv)%neighbors(i)
    if ( .not. connected(neighb) ) then
      tmp=get_dist(iajl,mv,neighb)!get_dist(iajl,mv,i)!ohd(i,mv)
      if ( tmp < hg ) then
        if ( mind(mv) + tmp < mind(neighb) ) then
          mind(neighb) = mind(mv) + tmp
          !new shortest path, reset saved paths
          pajl%nprevious(neighb)=0
          call add_prev(pajl,neighb,mv)
!          previous(neighb,0) = 1
!          previous(neighb,1)=mv
          count(neighb)=count(mv)
        else if(mind(mv) + tmp == mind(neighb) )then
            count(neighb)=count(neighb)+count(mv)
            call add_prev(pajl,neighb,mv)
!            previous(neighb,0)=previous(neighb,0)+1
!            previous(neighb,previous(neighb,0))=mv
        end if
      else
stop 'serious problem'
      end if
    end if
  end do

  return

end subroutine
end module
program main
use dij
use barrier_adjacency_list
use previous_graph
    implicit none
    character(len=10) :: ci1,ci2
    integer :: istart, idest, nmin,err,idmy,left,right,nc,i,nts,istat,j
    real(8) :: rdmy,tsener,mbener
    real(8) :: time1,time2
    real(8), allocatable :: ts(:),conn(:,:),minener(:)
    type(adjacency_list) :: ajl
    type(padjacency_list) :: pajl
    integer :: najl
    integer, allocatable :: neig(:,:),idx(:),idx2(:),paths(:,:)
    real(8), parameter :: thresh=huge(1)
!    real(8), parameter :: thresh=-169.5d0
!    real(8), parameter :: thresh=-389.d0
    integer :: nspaths,nconn,nmindat,nprev,ntsinpath
    logical :: mindatexists
    real(8) :: etmp,barrierener
    integer :: dest, depth, pathcounter,icurr,itmp,inext,iprev,barrierid
    character(len=100) :: ci
    character(len=500) :: outdir
    character(len=1500) :: line
    character(len=1000),allocatable :: tsfiledir(:)
    integer, allocatable :: tsfileid(:)
    logical :: tsfileexists
    

    nmindat=0
    mindatexists=.false.
    istart=1
    idest=2
    
    write(*,*)'output directory?'
    read(*,*)outdir
    write(*,*)'lowest barrier between which minima should be searched?'
    read(*,*)istart,idest

    !read minima energy
    inquire(file='mindat',exist=mindatexists)
    if(mindatexists)then
        open(33,file='mindat',status='old',action='read')
        nmindat=0
        do
            read(33,*,iostat=err)etmp
            if(err/=0)exit
            nmindat=nmindat+1
        enddo
        close(33)
        allocate(minener(nmindat))
        open(33,file='mindat',status='old',action='read')
        do i=1,nmindat
            read(33,*)minener(i)
        enddo
        close(33)
    endif

    !counter number of minima
    nmin=0
    nts=0
    open(33,file='tsdat',status='old',action='read',iostat=err)
        if(err/=0)then
            write(*,*)'cannot read tsdat.'
            stop
        endif
        do
            read(33,*,iostat=err)tsener,idmy,idmy,left,right
            if(err/=0)exit
            nts=nts+1
        enddo 
    close(33)
write(*,*)'nts',nts
    allocate(ts(nts),neig(2,nts))
    open(33,file='tsdat',status='old',action='read',iostat=err)
        if(err/=0)then
            write(*,*)'cannot read tsdat.'
            stop
        endif
        do i=1,nts
            read(33,*,iostat=err)ts(i),idmy,idmy,neig(1,i),neig(2,i)
            nmin=max(nmin,neig(1,i))
            nmin=max(nmin,neig(2,i))
        enddo
    close(33)

    inquire(file='tsfiles',exist=tsfileexists)
    if(tsfileexists)then
        open(34,file='tsfiles')
        i=0
        do
            read(34,'(a)',iostat=istat)line
            if(istat/=0)exit
            i=i+1
        enddo
        if(i/=nts)stop'Numer of files in tsfiles does not match number of TS in tsdat'
        close(34)
    endif

    if(nmin==huge(1))stop'too many minima'

!
  call cpu_time(time1)
    call dijkstra_pruned(nmin,nts,ts,neig,thresh,istart,idest,najl,ajl,idx,idx2,mbener,istat,nspaths,pajl,ntsinpath)
  call cpu_time(time2)
!
        write(ci1,'(i10)')istart
        write(ci2,'(i10)')idest
!    if(istat==0)then
!        write(*,*)'min. barr. between '//trim(adjustl(ci1))//' and '//trim(adjustl(ci2))//': '
!        write(*,*)mbener
!        write(*,*)
!    else
!        write(*,*)trim(adjustl(ci1))//' and '//trim(adjustl(ci2))//' not connected (below threshold energy)'
!        write(*,*)mbener
!        write(*,*)
!    endif
    write(*,*)'Number of shortest paths: ',nspaths
    if(nspaths>0)then
        write(*,*)'Number of transition states in shortest path: ',ntsinpath
        write(*,*)'highest energy along lowest energy path: ',mbener
        write(*,*)'time for finding lowest barrier path (s): ',time2-time1
    
        dest=idx2(idest)
        depth=0
        pathcounter=0
        call checkNumberpath(pajl,nmin,idx2,idx2(istart),dest,depth,pathcounter)
        if(pathcounter/=nspaths)then
            write(*,*)'STOP ERROR: pathcounter/=nspaths STOP'
            stop
        endif
        write(*,*)'pathcounter: ',pathcounter
        write(*,*)'Shortest paths as tree: '
        dest=idx2(idest)
        depth=0
        pathcounter=0
        allocate(paths(nspaths,ntsinpath+1))
        paths=huge(1)
        !print all pathes to fort.50, this does NOT destroy the 'previous' tree
        call printpath(pajl,nmin,idx2,idx2(istart),dest,depth,nspaths,ntsinpath,paths,pathcounter)
        write(*,*)
        write(*,*)'Shortests paths as minima in consecutive order: '
        do i=1,nspaths
            write(*,*)'Path no.',i
            do j=1,ntsinpath+1
                if(paths(i,j)==huge(1))then
                   paths(i,j)=paths(i-1,j) 
                endif
                write(ci,'(i0.0)')paths(i,j)
                write(*,'(xa)',advance='no')trim(adjustl(ci))
            enddo
            write(*,*)
        enddo

        if(tsfileexists)then
            allocate(tsfiledir(nts),tsfileid(nts))
            open(34,file='tsfiles')
            do i=1,nts
                read(34,*,iostat=istat)tsfiledir(i),tsfileid(i)
            enddo
            close(34)
        endif
        write(*,*)
        write(*,*)'Shortests paths as saddles in consecutive order: '
        if(mindatexists)then
        do i=1,nspaths
            write(ci,'(i0.0)')i
            open(173,file=trim(adjustl(outdir))//'/path_'//trim(adjustl(ci)))
            if(tsfileexists)then
            open(174,file=trim(adjustl(outdir))//'/saddles_'//trim(adjustl(ci))//'.inp')
            endif
        do j=1,ntsinpath
             write(173,*)minener(paths(i,j)),0
             call get_barrier(ajl,idx(paths(i,j)),idx(paths(i,j+1)),barrierener,barrierid)
             write(173,*)barrierener,3
             write(*,'(x,i0.0)',advance='no')barrierid
            if(tsfileexists)then
                write(174,'(a,x,i9.9)')'"'//trim(adjustl(tsfiledir(barrierid)))//'"',tsfileid(barrierid)
            endif
        enddo
            write(*,*)
             write(173,*)minener(paths(i,ntsinpath+1)),0
            close(173)
            if(tsfileexists)then
            close(174)
            endif
        enddo
            write(*,*)
        else
        write(*,*)'Mindat does not exist. Will not write pathfiles.'
        endif
            write(*,*)
    
        write(*,'(a)')'command for gnuplot:'
        write(*,'(a)')"plot './path_i' u 0:1 with lines, './path_i' u 0:1:2 with points linecolor variable"
        write(*,*)
    endif
end
recursive subroutine checkNumberpath(pajl,nmin,idx,start,dest,depth,pathcounter)
    use dij
    use previous_graph
    implicit none
    !parameters
    type(padjacency_list) :: pajl 
    integer :: nmin, idx(nmin)
    integer :: start, dest, depth
    integer :: pathcounter
    !internal
    integer :: i,j,idest
    character(len=100) :: ci
    idest=indexof(nmin,idx,dest)
    if(dest==start)then
        pathcounter=pathcounter+1
    endif
    write(ci,'(i0.0)')idest
    do i=1,pajl%nprevious(dest)!,prev(dest,0)
        call checkNumberpath(pajl,nmin,idx,start,get_prev(pajl,dest,i),depth+1,pathcounter)
    enddo
end subroutine

recursive subroutine printpath(pajl,nmin,idx,start,dest,depth,nspaths,ntsinpath,paths,pathcounter)
    use dij
    use previous_graph
    implicit none
    !parameters
    type(padjacency_list) :: pajl 
    integer :: nmin, idx(nmin)
    integer :: start, dest, depth
    integer :: nspaths,ntsinpath, paths(nspaths,ntsinpath+1)
    integer :: pathcounter
    !internal
    integer :: i,j,idest
    character(len=100) :: ci
    idest=indexof(nmin,idx,dest)
    write(ci,'(i0.0)')idest
    write(*,'(a)',advance='no')'-'//trim(adjustl(ci))
    write(*,*)
    paths(pathcounter+1,depth+1)=idest
    if(dest==start)then
        pathcounter=pathcounter+1
    endif
    do i=1,pajl%nprevious(dest)!prev(dest,0)
        do j=0, depth
         write(*,'(a)',advance='no')' |'
        enddo
        call printpath(pajl,nmin,idx,start,get_prev(pajl,dest,i),depth+1,nspaths,ntsinpath,paths,pathcounter)
    enddo
end subroutine
