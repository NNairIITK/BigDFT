!> template for the wpdot routine

  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1
  real(wp) :: scpr1,scpr0
  integer :: iaseg0,ibsegs,ibsege
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_wp

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,scpr0,scpr1) &
  !$omp private(jbj,ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  scpr0=0.0_wp
  scpr1=0.0_wp

!!!!start of general region

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr0)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr1)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 

 !!!!end of general region

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel
