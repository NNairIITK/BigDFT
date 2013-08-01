program kpts_check
  use BigDFT_API
  implicit none
  integer :: norb,nvctr,nproc,nkpts,norbv
  integer, dimension(:,:), allocatable :: norb_par,nvctr_par,norb_parv,nvctr_parv
  !local variables
  character(len=*), parameter :: subname='kpts_check'
  character(len=20) :: argstring
  integer :: i_all,i_stat,info,lubo,lubc

  !@todo : load balancing verifications of particular k-points distributions

  call getarg(1,argstring)
  read(argstring,*)nproc
  call getarg(2,argstring)
  read(argstring,*)nkpts
  call getarg(3,argstring)
  read(argstring,*)norb
  call getarg(4,argstring)
  read(argstring,*)nvctr
  call getarg(5,argstring)
  read(argstring,*)norbv
!
!!$  !ready for the test, should print the maximum load unbalancing
!!$  do norb=10,10000,165
!!$     do nvctr=32456,32456
!!$        do nkpts=1,200,14
!!$           do nproc=1,10
              write(*,*)'TEST nproc,nkpts,norb,nvctr',nproc,nkpts,norb,nvctr
              allocate(norb_par(0:nproc-1,nkpts+ndebug),stat=i_stat)
              call memocc(i_stat,norb_par,'norb_par',subname)
              allocate(nvctr_par(0:nproc-1,nkpts+ndebug),stat=i_stat)
              call memocc(i_stat,nvctr_par,'nvctr_par',subname)
              
              call kpts_to_procs_via_obj(nproc,nkpts,norb,norb_par)
              call kpts_to_procs_via_obj(nproc,nkpts,nvctr,nvctr_par)
              info=0
              call check_kpt_distributions(nproc,nkpts,norb,nvctr,norb_par,nvctr_par,info,lubo,lubc)
              if (info/=0) then !redo the distribution based on the orbitals scheme
                 info=0
                 call components_kpt_distribution(nproc,nkpts,norb,nvctr,norb_par,nvctr_par)
                 call check_kpt_distributions(nproc,nkpts,norb,nvctr,norb_par,nvctr_par,info,lubo,lubc)
              end if
              if (info /=0) then
                 write(*,*)'ERROR for nproc,nkpts,norb,nvctr',nproc,nkpts,norb,nvctr
                 stop 'info'
              end if
!once the data for normal optimization have been decided, check for Davidson distribution

              write(*,*)'TESTD nproc,nkpts,norbv,nvctrv',nproc,nkpts,norbv,nvctr
              allocate(norb_parv(0:nproc-1,nkpts+ndebug),stat=i_stat)
              call memocc(i_stat,norb_parv,'norb_parv',subname)
              allocate(nvctr_parv(0:nproc-1,nkpts+ndebug),stat=i_stat)
              call memocc(i_stat,nvctr_parv,'nvctr_parv',subname)

              !do the orbital repartition for the virtual orbitals
              call components_kpt_distribution(nproc,nkpts,norb,norbv,norb_par,norb_parv)
              !associate the same distribution for the components
              nvctr_parv=nvctr_par
              !then check
              info=0
              call check_kpt_distributions(nproc,nkpts,norbv,nvctr,norb_parv,nvctr_parv,info,lubo,lubc)
              if (info /=0) then
                 write(*,*)'ERROR, VIRTUAL for nproc,nkpts,norb,nvctr,norbv',nproc,nkpts,norb,nvctr,norbv
                 stop 'info'
              end if
              
 
               i_all=-product(shape(nvctr_parv))*kind(nvctr_parv)
              deallocate(nvctr_parv,stat=i_stat)
              call memocc(i_stat,i_all,'nvctr_parv',subname)
              i_all=-product(shape(norb_parv))*kind(norb_parv)
              deallocate(norb_parv,stat=i_stat)
              call memocc(i_stat,i_all,'norb_parv',subname)
!             
              i_all=-product(shape(nvctr_par))*kind(nvctr_par)
              deallocate(nvctr_par,stat=i_stat)
              call memocc(i_stat,i_all,'nvctr_par',subname)
              i_all=-product(shape(norb_par))*kind(norb_par)
              deallocate(norb_par,stat=i_stat)
              call memocc(i_stat,i_all,'norb_par',subname)
!!$           end do
!!$        end do
!!$     end do
!!$  end do



end program kpts_check
