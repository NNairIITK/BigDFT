subroutine copy_old_wavefunctions(iproc,nproc,norb,norbp,nspinor,hgrid,n1,n2,n3,eval,wfd,psi,&
     hgrid_old,n1_old,n2_old,n3_old,eval_old,wfd_old,psi_old)
  use module_base
  use module_types

  implicit none
  type(wavefunctions_descriptors) :: wfd,wfd_old
  integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,nspinor
  real(kind=8), intent(in) :: hgrid
  integer, intent(out) :: n1_old,n2_old,n3_old
  real(kind=8), intent(out) :: hgrid_old
  real(kind=8), dimension(:), pointer :: eval,eval_old
  real(kind=8), dimension(:), pointer :: psi,psi_old
  !local variables
  character(len=*), parameter :: subname='copy_old_wavefunctions'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: iseg,nvctrp_old,i1,i2,j,ind1,ind2,iorb,i_all,i_stat,oidx,sidx
  real(kind=8) :: tt

  wfd_old%nvctr_c = wfd%nvctr_c
  wfd_old%nvctr_f = wfd%nvctr_f
  wfd_old%nseg_c  = wfd%nseg_c
  wfd_old%nseg_f  = wfd%nseg_f

  !allocations
  call allocate_wfd(wfd_old,'copy_old_wavefunctions')

  do iseg=1,wfd_old%nseg_c+wfd_old%nseg_f
     wfd_old%keyg(1,iseg)    = wfd%keyg(1,iseg)
     wfd_old%keyg(2,iseg)    = wfd%keyg(2,iseg)
     wfd_old%keyv(iseg)      = wfd%keyv(iseg)
  enddo
  !deallocation
  call deallocate_wfd(wfd,'copy_old_wavefunctions')

  hgrid_old   = hgrid
  n1_old      = n1
  n2_old      = n2
  n3_old      = n3
  !add the number of distributed point for the compressed wavefunction
  tt=dble(wfd_old%nvctr_c+7*wfd_old%nvctr_f)/dble(nproc)
  nvctrp_old=int((1.d0-eps_mach*tt) + tt)

  allocate(psi_old((wfd_old%nvctr_c+7*wfd_old%nvctr_f)*norbp*nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,psi_old,'psi_old',subname)
  allocate(eval_old(norb+ndebug),stat=i_stat)
  call memocc(i_stat,eval_old,'eval_old',subname)

!  do iorb=iproc*norbp*nspinor+1,min((iproc+1)*norbp,norb)*nspinor,nspinor
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     tt=0.d0
     oidx=(iorb-1)*nspinor+1-iproc*norbp*nspinor
     do sidx=oidx,oidx+nspinor-1
        do j=1,wfd_old%nvctr_c+7*wfd_old%nvctr_f
           ind1=j+(wfd_old%nvctr_c+7*wfd_old%nvctr_f)*(sidx-1)
           psi_old(ind1)= psi(ind1)
           tt=tt+psi(ind1)**2
        enddo
     end do

     
!!$     do sidx=oidx,oidx+nspinor-1
!!$     !the control between the different psi
!!$     !must be such that the allocation dimensions are respected
!!$     !(remember that the allocations are slightly enlarged to avoid
!!$     !allocation of work arrays for transposition in the parallel case)
!!$!     do sidx=0,nspinor-1
!!$        do j=1,wfd_old%nvctr_c+7*wfd_old%nvctr_f
!!$           !starting address of the direct array in transposed form
!!$           !changed for obrtaining the good one
!!$           call trans_address(nvctrp_old,wfd_old%nvctr_c+7*wfd_old%nvctr_f,&
!!$                j,sidx, i1,i2)
!!$        !call trans_address(nvctrp_old*nspinor,(wfd_old%nvctr_c+7*wfd_old%nvctr_f)*nspinor,&
!!$        !     j,sidx, i1,i2)
!!$           psi_old(j,sidx)     = psi(i1,i2)
!!$           tt=tt+psi(i1,i2)**2
!!$        enddo
!!$     end do
     tt=sqrt(tt)
     if (abs(tt-1.d0).gt.1.d-8) then
        write(*,*)'wrong psi_old',iorb,tt
        stop 
     end if
     eval_old(iorb) = eval(iorb)
  enddo
  !deallocation
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(eval))*kind(eval)
  deallocate(eval,stat=i_stat)
  call memocc(i_stat,i_all,'eval',subname)

end subroutine copy_old_wavefunctions

subroutine reformatmywaves(iproc,norb,norbp,nat,&
     & hgrid_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
     & hgrid,n1,n2,n3,rxyz,wfd,psi)
  use module_base
  use module_types

  implicit real(kind=8) (a-h,o-z)
  type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
  dimension :: rxyz(3,nat), rxyz_old(3,nat)
  dimension :: psi_old(wfd_old%nvctr_c + 7 * wfd_old%nvctr_f, norbp), &
       psi(wfd%nvctr_c + 7 * wfd%nvctr_f, norbp)
  character(len=*), parameter :: subname='reformatmywaves'
  logical :: reformat

  allocatable :: psifscf(:,:,:), psigold(:,:,:,:,:,:)

  allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8+ndebug),stat=i_stat)
  call memocc(i_stat,psifscf,'psifscf',subname)

  tx=0.d0 ; ty=0.d0 ; tz=0.d0
  do iat=1,nat
  tx=tx+(rxyz(1,iat)-rxyz_old(1,iat))**2
  ty=ty+(rxyz(2,iat)-rxyz_old(2,iat))**2
  tz=tz+(rxyz(3,iat)-rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)
!  write(100+iproc,*) 'displacement',dis
!  write(100+iproc,*) 'rxyz ',rxyz
!  write(100+iproc,*) 'rxyz_old ',rxyz_old

  !reformatting criterion
  if (hgrid_old.eq. hgrid .and. wfd_old%nvctr_c.eq.wfd%nvctr_c .and. &
       wfd_old%nvctr_f .eq. wfd%nvctr_f .and.&
       n1_old.eq.n1  .and. n2_old.eq.n2 .and. n3_old.eq.n3  .and.  displ.lt.1.d-3  ) then
     reformat=.false.
     if (iproc==0) then
        write(*,'(1x,a)',advance='NO')&
         'The wavefunctions do not need reformatting and can be imported directly...   '
       !  '-------------------------------------------------------------- Wavefunctions Restart'
     end if
  else
     reformat=.true.
     if (iproc==0) then
        write(*,'(1x,a)')&
         'The wavefunctions need reformatting because:                                 '
        if (hgrid_old.ne.hgrid) then 
           write(*,"(4x,a,1pe20.12)") &
                '  hgrid_old >< hgrid  ',hgrid_old, hgrid
        else if (wfd_old%nvctr_c.ne.wfd%nvctr_c) then
           write(*,"(4x,a,2i8)") &
                'nvctr_c_old >< nvctr_c',wfd_old%nvctr_c,wfd%nvctr_c
        else if (wfd_old%nvctr_f.ne.wfd%nvctr_f)  then
           write(*,"(4x,a,2i8)") &
                'nvctr_f_old >< nvctr_f',wfd_old%nvctr_f,wfd%nvctr_f
        else if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 )  then  
           write(*,"(4x,a,6i5)") &
                'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
        else
           write(*,"(4x,a,3(1pe19.12))") &
                'molecule was shifted  ' , tx,ty,tz
        endif
           write(*,"(1x,a)",advance='NO')& 
                'Reformatting...'
     end if
!check
!        write(100+iproc,'(1x,a)')&
!         'The wavefunctions need reformatting because:                                 '
!        if (hgrid_old.ne.hgrid) then 
!           write(100+iproc,"(4x,a,1pe20.12)") &
!                '  hgrid_old >< hgrid  ',hgrid_old, hgrid
!        else if (wfd_old%nvctr_c.ne.wfd%nvctr_c) then
!           write(100+iproc,"(4x,a,2i8)") &
!                'nvctr_c_old >< nvctr_c',wfd_old%nvctr_c,wfd%nvctr_c
!        else if (wfd_old%nvctr_f.ne.wfd%nvctr_f)  then
!           write(100+iproc,"(4x,a,2i8)") &
!                'nvctr_f_old >< nvctr_f',wfd_old%nvctr_f,wfd%nvctr_f
!        else if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 )  then  
!           write(100+iproc,"(4x,a,6i5)") &
!                'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
!        else
!           write(100+iproc,"(4x,a,3(1pe19.12))") &
!                'molecule was shifted  ' , tx,ty,tz
!        endif
!checkend
  end if


  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

     if (.not. reformat) then
!write(100+iproc,*) 'no reformatting' 

        do j=1,wfd_old%nvctr_c
           psi(j,iorb-iproc*norbp)=psi_old(j, iorb - iproc * norbp)
        enddo
        do j=1,7*wfd_old%nvctr_f-6,7
           psi(wfd%nvctr_c+j+0,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+0,iorb-iproc*norbp)
           psi(wfd%nvctr_c+j+1,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+1,iorb-iproc*norbp)
           psi(wfd%nvctr_c+j+2,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+2,iorb-iproc*norbp)
           psi(wfd%nvctr_c+j+3,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+3,iorb-iproc*norbp)
           psi(wfd%nvctr_c+j+4,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+4,iorb-iproc*norbp)
           psi(wfd%nvctr_c+j+5,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+5,iorb-iproc*norbp)
           psi(wfd%nvctr_c+j+6,iorb-iproc*norbp)=psi_old(wfd%nvctr_c+j+6,iorb-iproc*norbp)
        enddo

     else

        allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2+ndebug),stat=i_stat)
        call memocc(i_stat,psigold,'psigold',subname)

        call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)


        ! coarse part
        do iseg=1,wfd_old%nseg_c
           jj=wfd_old%keyv(iseg)
           j0=wfd_old%keyg(1,iseg)
           j1=wfd_old%keyg(2,iseg)
           ii=j0-1
           i3=ii/((n1_old+1)*(n2_old+1))
           ii=ii-i3*(n1_old+1)*(n2_old+1)
           i2=ii/(n1_old+1)
           i0=ii-i2*(n1_old+1)
           i1=i0+j1-j0
           do i=i0,i1
              psigold(i,1,i2,1,i3,1) = psi_old(i-i0+jj,iorb-iproc*norbp)
           enddo
        enddo

        ! fine part
        do iseg=1,wfd_old%nseg_f
           jj=wfd_old%keyv(wfd_old%nseg_c + iseg)
           j0=wfd_old%keyg(1,wfd_old%nseg_c + iseg)
           j1=wfd_old%keyg(2,wfd_old%nseg_c + iseg)
           ii=j0-1
           i3=ii/((n1_old+1)*(n2_old+1))
           ii=ii-i3*(n1_old+1)*(n2_old+1)
           i2=ii/(n1_old+1)
           i0=ii-i2*(n1_old+1)
           i1=i0+j1-j0
           do i=i0,i1
              psigold(i,2,i2,1,i3,1)=psi_old(wfd_old%nvctr_c+1+7*(i-i0+jj-1), iorb-iproc*norbp)
              psigold(i,1,i2,2,i3,1)=psi_old(wfd_old%nvctr_c+2+7*(i-i0+jj-1), iorb-iproc*norbp)
              psigold(i,2,i2,2,i3,1)=psi_old(wfd_old%nvctr_c+3+7*(i-i0+jj-1), iorb-iproc*norbp)
              psigold(i,1,i2,1,i3,2)=psi_old(wfd_old%nvctr_c+4+7*(i-i0+jj-1), iorb-iproc*norbp)
              psigold(i,2,i2,1,i3,2)=psi_old(wfd_old%nvctr_c+5+7*(i-i0+jj-1), iorb-iproc*norbp)
              psigold(i,1,i2,2,i3,2)=psi_old(wfd_old%nvctr_c+6+7*(i-i0+jj-1), iorb-iproc*norbp)
              psigold(i,2,i2,2,i3,2)=psi_old(wfd_old%nvctr_c+7+7*(i-i0+jj-1), iorb-iproc*norbp)
           enddo
        enddo

!write(100+iproc,*) 'norm psigold ',dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold,1)

        call reformatonewave(iproc, displ, hgrid_old, &
             & n1_old, n2_old, n3_old, nat, rxyz_old, psigold, hgrid, &
             & wfd%nvctr_c, wfd%nvctr_f, n1, n2, n3, rxyz, &
             & wfd%nseg_c, wfd%nseg_f, wfd%keyg, wfd%keyv, psifscf, & 
             & psi(1,iorb - iproc * norbp))

        i_all=-product(shape(psigold))*kind(psigold)
        deallocate(psigold,stat=i_stat)
        call memocc(i_stat,i_all,'psigold',subname)
     end if
  end do

  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf',subname)

  if (iproc==0) write(*,"(1x,a)")'done.'

END SUBROUTINE reformatmywaves

subroutine readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz_old,rxyz,  & 
     wfd,psi,eval)
  ! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
  use module_base
  use module_types

  implicit real(kind=8) (a-h,o-z)
  type(wavefunctions_descriptors), intent(in) :: wfd
  character(len=*), parameter :: subname='readmywaves'
  character(len=50) filename
  character(len=4) f4
  dimension psi(wfd%nvctr_c+7*wfd%nvctr_f,norbp)
  dimension rxyz_old(3,nat),rxyz(3,nat),eval(norb)
  allocatable :: psifscf(:,:,:)

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8+ndebug),stat=i_stat)
  call memocc(i_stat,psifscf,'psifscf',subname)

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

     write(f4,'(i4.4)') iorb
     filename = 'wavefunction.'//f4
     open(unit=99,file=filename,status='unknown')

     call readonewave(99, .true., iorb,iproc,n1,n2,n3, &
          & hgrid,nat,rxyz_old,rxyz,wfd%nseg_c,wfd%nseg_f,wfd%nvctr_c,wfd%nvctr_f,wfd%keyg,wfd%keyv,&
          psi(1,iorb-iproc*norbp),eval(iorb),psifscf)
     close(99)
  end do

  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf',subname)

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  write(*,'(a,i4,2(1x,e10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
END SUBROUTINE readmywaves

subroutine writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,wfd,psi,eval)
  ! write all my wavefunctions in files by calling writeonewave
  
  use module_types

  implicit real(kind=8) (a-h,o-z)
  type(wavefunctions_descriptors), intent(in) :: wfd
  character(len=4) f4
  character(len=50) filename
  dimension rxyz(3,nat),eval(norb)
  dimension psi(wfd%nvctr_c+7*wfd%nvctr_f,norbp)

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

     write(f4,'(i4.4)') iorb
     filename = 'wavefunction.'//f4
     write(*,*) 'opening ',filename
     open(unit=99,file=filename,status='unknown')

     call writeonewave(99, .true., iorb,n1,n2,n3,hgrid,nat,rxyz,  & 
          wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1)  & 
          ,wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1), & 
          psi(1,iorb-iproc*norbp),psi(wfd%nvctr_c+1,iorb-iproc*norbp),norb,eval)
     close(99)

  enddo

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  write(*,'(a,i4,2(1x,e10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel


  return
END SUBROUTINE writemywaves
