!wrapper for simplifying the call
subroutine wnrm_wrap(mvctr_c,mvctr_f,psi,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(mvctr_c+7*mvctr_f), intent(in) :: psi
  real(dp), intent(out) :: scpr
  !local variables
  integer :: i_f

  i_f=min(mvctr_f,1)
 
  call wnrm(mvctr_c,mvctr_f,psi,psi(mvctr_c+i_f),scpr)
  
end subroutine wnrm_wrap

! calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(dp), intent(out) :: scpr
  !local variables
  integer :: i,ithread,nthread
  real(dp) :: pc,pf1,pf2,pf3,pf4,pf5,pf6,pf7
  real(dp) :: scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7

!!$omp parallel default(private) shared(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
!!$    ithread=omp_get_thread_num()
!!$    nthread=omp_get_num_threads()
    scpr0=0.0_dp
!!$  if (ithread .eq. 0) then
  do i=1,mvctr_c
     !scpr0=scpr0+psi_c(i)**2
     pc=real(psi_c(i),dp)
     scpr0=scpr0+pc**2
  enddo
!!$  endif
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!!$  if (ithread .eq. 1  .or. nthread .eq. 1) then
  do i=1,mvctr_f
!!$     scpr1=scpr1+psi_f(1,i)**2
!!$     scpr2=scpr2+psi_f(2,i)**2
!!$     scpr3=scpr3+psi_f(3,i)**2
!!$     scpr4=scpr4+psi_f(4,i)**2
!!$     scpr5=scpr5+psi_f(5,i)**2
!!$     scpr6=scpr6+psi_f(6,i)**2
!!$     scpr7=scpr7+psi_f(7,i)**2
     pf1=real(psi_f(1,i),dp)
     pf2=real(psi_f(2,i),dp)
     pf3=real(psi_f(3,i),dp)
     pf4=real(psi_f(4,i),dp)
     pf5=real(psi_f(5,i),dp)
     pf6=real(psi_f(6,i),dp)
     pf7=real(psi_f(7,i),dp)
     scpr1=scpr1+pf1**2
     scpr2=scpr2+pf2**2
     scpr3=scpr3+pf3**2
     scpr4=scpr4+pf4**2
     scpr5=scpr5+pf5**2
     scpr6=scpr6+pf6**2
     scpr7=scpr7+pf7**2
  enddo
!!$  endif
!!$omp critical
  scpr=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
!!$omp end critical
!!$omp end parallel

end subroutine wnrm


!wrapper for simplifying the call
subroutine wscal_wrap(mvctr_c,mvctr_f,scal,psi)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), intent(in) :: scal
  real(wp), dimension(mvctr_c+7*mvctr_f), intent(in) :: psi
  !local variables
  integer :: i_f

  i_f=min(mvctr_f,1)
 
  call wscal(mvctr_c,mvctr_f,scal,psi,psi(mvctr_c+i_f))
  
end subroutine wscal_wrap

! multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
subroutine wscal(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: i
!$omp parallel default(private) shared(mvctr_c,mvctr_f,scal,psi_c,psi_f)
!$omp do
  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal
  enddo
!$omp enddo
!$omp do
  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal
     psi_f(2,i)=psi_f(2,i)*scal
     psi_f(3,i)=psi_f(3,i)*scal
     psi_f(4,i)=psi_f(4,i)*scal
     psi_f(5,i)=psi_f(5,i)*scal
     psi_f(6,i)=psi_f(6,i)*scal
     psi_f(7,i)=psi_f(7,i)*scal
  enddo
!$omp enddo
!$omp end parallel

end subroutine wscal

!wrapper for simplifying the call
subroutine wscalv_wrap(mvctr_c,mvctr_f,scal,psi)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c+7*mvctr_f), intent(in) :: psi
  !local variables
  integer :: i_f

  i_f=min(mvctr_f,1)
 
  call wscalv(mvctr_c,mvctr_f,scal,psi,psi(mvctr_c+i_f))
  
end subroutine wscalv_wrap


! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
subroutine wscalv(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: i
!$omp parallel default(private) shared(mvctr_c,mvctr_f,scal,psi_c,psi_f)
!$omp do
  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal(0)           !  1 1 1
  enddo
!$omp enddo
!$omp do
  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
     psi_f(2,i)=psi_f(2,i)*scal(1)       !  1 2 1
     psi_f(3,i)=psi_f(3,i)*scal(2)       !  2 2 1
     psi_f(4,i)=psi_f(4,i)*scal(1)       !  1 1 2
     psi_f(5,i)=psi_f(5,i)*scal(2)       !  2 1 2
     psi_f(6,i)=psi_f(6,i)*scal(2)       !  1 2 2
     psi_f(7,i)=psi_f(7,i)*scal(3)       !  2 2 2
  enddo
!$omp enddo
!$omp end parallel
end subroutine wscalv

! initializes a wavefunction to zero
subroutine wzero(mvctr_c,mvctr_f,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: i

! seems to be never called
!write(*,*) ' i am in wzero'
!$omp parallel default(private) shared(mvctr_c,mvctr_f,psi_c,psi_f)
!$omp do
  do i=1,mvctr_c
     psi_c(i)=0.0_wp
  enddo
!$omp enddo
!$omp do
  do i=1,mvctr_f
     psi_f(1,i)=0.0_wp
     psi_f(2,i)=0.0_wp
     psi_f(3,i)=0.0_wp
     psi_f(4,i)=0.0_wp
     psi_f(5,i)=0.0_wp
     psi_f(6,i)=0.0_wp
     psi_f(7,i)=0.0_wp
  enddo
!$omp enddo
!$omp end parallel
end subroutine wzero


!wrapper of wpdot to avoid boundary problems in absence of wavelets
subroutine wpdot_wrap(mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
  real(wp), dimension(mavctr_c+7*mavctr_f), intent(in) :: apsi
  real(wp), dimension(mbvctr_c+7*mbvctr_f), intent(in) :: bpsi
  real(dp), intent(out) :: scpr
  !local variables
  integer :: ia_f,ib_f,iaseg_f,ibseg_f

  ia_f=min(mavctr_f,1)
  ib_f=min(mbvctr_f,1)

  iaseg_f=min(maseg_f,1)
  ibseg_f=min(mbseg_f,1)


  call wpdot(mavctr_c,mavctr_f,maseg_c,maseg_f,&
       keyav,keyav(maseg_c+iaseg_f),&
       keyag,keyag(1,maseg_c+iaseg_f),&
       apsi,apsi(mavctr_c+ia_f),  &
       mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
       keybv,keybv(mbseg_c+ibseg_f),&
       keybg,keybg(1,mbseg_c+ibseg_f),&
       bpsi,bpsi(mbvctr_c+ib_f),scpr)

end subroutine wpdot_wrap

!this function must be generalized for the linear scaling code
! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
! Warning: the subroutine assumes that bpsi has only one segment along each line,
! whereas apsi can have several segments. This assumption is true if bpsi is a projector 
! To be more precise, it is assumed that the segments of bpsi are always contained inside
! the segments of apsi, no matter whether they are in the same line or not.
subroutine wpdot(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(dp), intent(out) :: scpr
  !local variables
  integer :: iaseg,ibseg,llc,jaj,ja0,ja1,jb1,jb0,jbj,iaoff,iboff,length,llf,i,ithread,nthread
  real(dp) :: pac,paf1,paf2,paf3,paf4,paf5,paf6,paf7,pbc,pbf1,pbf2,pbf3,pbf4,pbf5,pbf6,pbf7
  real(dp) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  !  integer :: ncount0,ncount2,ncount_rate,ncount_max
  !  real(gp) :: tel


  !dee
  !  open(unit=97,file='time_wpdot',status='unknown')
  !  call system_clock(ncount0,ncount_rate,ncount_max)

  scpr=0.0_dp

  !dee
!$omp parallel default (private) &
!$omp shared (maseg_c,keyav_c,keyag_c,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
!$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,mbseg_f,keybg_f,keyag_f,keyav_f)&
!$omp shared (apsi_f,scpr)
!$    ithread=omp_get_thread_num()
!$    nthread=omp_get_num_threads()
    scpr0=0.0_dp
!$  if (ithread .eq. 0) then
  llc=0
  !coarse part
  ibseg=1
  !for each segment of the first function
  loop_jac: do iaseg=1,maseg_c
     jaj=keyav_c(iaseg)
     ja0=keyag_c(1,iaseg)
     ja1=keyag_c(2,iaseg)

     !control if it intersects a segment of the second function
     loop_jbc: do
        jb1=keybg_c(2,ibseg)
        jb0=keybg_c(1,ibseg)
        if (jb1 < ja0) then !not yet there increase the count
           ibseg=ibseg+1
           if (ibseg > mbseg_c) exit loop_jac !function b ended
           cycle loop_jbc ! next segment
        end if
        if (jb0 > ja1) exit loop_jbc  !we went through, don't increase count, exit inner loop

        jbj=keybv_c(ibseg)
        if (ja0 .gt. jb0) then 
           iaoff=0
           iboff=ja0-jb0
           length=min(ja1,jb1)-ja0
        else
           iaoff=jb0-ja0
           iboff=0
           length=min(ja1,jb1)-jb0
        endif
        !write(*,*) 'ja0,ja1,jb0,jb1',ja0,ja1,jb0,jb1,length
        !write(*,'(5(a,i5))') 'C:from ',jaj+iaoff,' to ',jaj+iaoff+length,' and from ',jbj+iboff,' to ',jbj+iboff+length
        do i=0,length
           llc=llc+1
           pac=real(apsi_c(jaj+iaoff+i),dp)
           pbc=real(bpsi_c(jbj+iboff+i),dp)
           scpr0=scpr0+pac*pbc
        enddo
        ibseg=ibseg+1
        if (ibseg > mbseg_c) exit loop_jac !function b ended 
     end do loop_jbc
  enddo loop_jac

!$  endif
  !print *,'nvctr_c',llc,mavctr_c,mbvctr_c


  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!$  if (ithread .eq. 1  .or. nthread .eq. 1) then
  llf=0
  ! fine part
  !add possibility of zero fine segments for the projectors
  ibseg=1
  if (mbseg_f /= 0) then
     loop_jaf: do iaseg=1,maseg_f
        jaj=keyav_f(iaseg)
        ja0=keyag_f(1,iaseg)
        ja1=keyag_f(2,iaseg)

        loop_jbf: do
           jb1=keybg_f(2,ibseg)
           jb0=keybg_f(1,ibseg)
           if (jb1 < ja0) then
              ibseg=ibseg+1
              if (ibseg > mbseg_f) exit loop_jbf
              cycle loop_jbf
           end if
           if (jb0 > ja1) exit loop_jbf

           jbj=keybv_f(ibseg)
           if (ja0 .gt. jb0) then 
              iaoff=0
              iboff=ja0-jb0
              length=min(ja1,jb1)-ja0
           else
              iaoff=jb0-ja0
              iboff=0
              length=min(ja1,jb1)-jb0
           endif
           do i=0,length
              llf=llf+1
              paf1=real(apsi_f(1,jaj+iaoff+i),dp)
              pbf1=real(bpsi_f(1,jbj+iboff+i),dp)
              paf2=real(apsi_f(2,jaj+iaoff+i),dp)
              pbf2=real(bpsi_f(2,jbj+iboff+i),dp)
              paf3=real(apsi_f(3,jaj+iaoff+i),dp)
              pbf3=real(bpsi_f(3,jbj+iboff+i),dp)
              paf4=real(apsi_f(4,jaj+iaoff+i),dp)
              pbf4=real(bpsi_f(4,jbj+iboff+i),dp)
              paf5=real(apsi_f(5,jaj+iaoff+i),dp)
              pbf5=real(bpsi_f(5,jbj+iboff+i),dp)
              paf6=real(apsi_f(6,jaj+iaoff+i),dp)
              pbf6=real(bpsi_f(6,jbj+iboff+i),dp)
              paf7=real(apsi_f(7,jaj+iaoff+i),dp)
              pbf7=real(bpsi_f(7,jbj+iboff+i),dp)

              scpr1=scpr1+paf1*pbf1
              scpr2=scpr2+paf2*pbf2
              scpr3=scpr3+paf3*pbf3
              scpr4=scpr4+paf4*pbf4
              scpr5=scpr5+paf5*pbf5
              scpr6=scpr6+paf6*pbf6
              scpr7=scpr7+paf7*pbf7
           enddo
           ibseg=ibseg+1
           if (ibseg > mbseg_f) exit loop_jaf
        end do loop_jbf
     enddo loop_jaf
  end if

!$  endif
  !print *,'nvctr_f',llf,mavctr_f,mbvctr_f
!$omp critical 
  scpr=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
!$omp end critical
!$omp end parallel
  !        write(*,*) 'llc,llf',llc,llf
  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount0)/dble(ncount_rate)
  !  write(97,'(a40,1x,e10.3,1x,f6.1)') 'wpdot:',tel
  !  close(97)

end subroutine wpdot

subroutine waxpy_wrap(scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,&
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  real(dp), intent(in) :: scpr
  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
  real(wp), dimension(mbvctr_c+7*mbvctr_f), intent(in) :: bpsi
  real(wp), dimension(mavctr_c+7*mavctr_f), intent(inout) :: apsi

  !local variables
  integer :: ia_f,ib_f,iaseg_f,ibseg_f

  ia_f=min(mavctr_f,1)
  ib_f=min(mbvctr_f,1)

  iaseg_f=min(maseg_f,1)
  ibseg_f=min(mbseg_f,1)

  call waxpy(scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
       keybv,keybv(mbseg_c+ibseg_f),&
       keybg,keybg(1,mbseg_c+ibseg_f),&
       bpsi,bpsi(mbvctr_c+ib_f), &
       mavctr_c,mavctr_f,maseg_c,maseg_f,&
       keyav,keyav(maseg_c+iaseg_f),&
       keyag,keyag(1,maseg_c+iaseg_f),&
       apsi,apsi(mavctr_c+ia_f))

end subroutine waxpy_wrap


! rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
! The update is only done in the localization region of apsi
subroutine waxpy(  & 
     scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f, & 
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  real(dp), intent(in) :: scpr
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f
  !local variables
  integer :: iaseg,ibseg,jaj,ja0,ja1,jb1,jb0,jbj,iaoff,iboff,length,i,ithread,nthread
  !  integer :: ncount0,ncount2,ncount_rate,ncount_max
  !  real(gp) :: tel 
  real(wp) :: scprwp
  !dee
  !  open(unit=97,file='time_waxpy',status='unknown')
  !  call system_clock(ncount0,ncount_rate,ncount_max)

  scprwp=real(scpr,wp)
  !dee
!$omp parallel default (private) &
!$omp shared (maseg_c,keyav_c,keyag_c,keybg_c,mbseg_c,mbseg_f,maseg_f)&
!$omp shared (keyav_f,keyag_f,keybg_f,keybv_f,scprwp,bpsi_c,bpsi_f)&
!$omp shared (apsi_f,apsi_c,keybv_c)
!$   ithread=omp_get_thread_num()
!$   nthread=omp_get_num_threads()
  !        llc=0
  ! coarse part
  ibseg=1

!$  if (ithread .eq. 0) then
     loop_jac: do iaseg=1,maseg_c
        jaj=keyav_c(iaseg)
        ja0=keyag_c(1,iaseg)
        ja1=keyag_c(2,iaseg)

        loop_jbc: do
           jb1=keybg_c(2,ibseg)
           jb0=keybg_c(1,ibseg)
           if (jb1 < ja0) then
              ibseg=ibseg+1
              if (ibseg > mbseg_c) exit loop_jac
              cycle loop_jbc
           end if
           if (jb0 > ja1) exit loop_jbc

           jbj=keybv_c(ibseg)
           if (ja0 .gt. jb0) then 
              iaoff=0
              iboff=ja0-jb0
              length=min(ja1,jb1)-ja0
           else
              iaoff=jb0-ja0
              iboff=0
              length=min(ja1,jb1)-jb0
           endif
           do i=0,length
              !          llc=llc+1
              apsi_c(jaj+iaoff+i)=apsi_c(jaj+iaoff+i)+scprwp*bpsi_c(jbj+iboff+i) 
           enddo
           ibseg=ibseg+1
           if (ibseg > mbseg_c) exit loop_jac
        end do loop_jbc
     enddo loop_jac
!$  endif

  !        llf=0
  ! fine part
  ibseg=1
!$  if (ithread .eq. 1 .or. nthread .eq. 1) then
     if (mbseg_f /= 0) then
        loop_jaf: do iaseg=1,maseg_f
           jaj=keyav_f(iaseg)
           ja0=keyag_f(1,iaseg)
           ja1=keyag_f(2,iaseg)

           loop_jbf: do
              jb1=keybg_f(2,ibseg)
              jb0=keybg_f(1,ibseg)
              if (jb1 < ja0) then
                 ibseg=ibseg+1
                 if (ibseg > mbseg_f) exit loop_jbf
                 cycle loop_jbf
              end if
              if (jb0 > ja1) exit loop_jbf

              jbj=keybv_f(ibseg)
              if (ja0 .gt. jb0) then 
                 iaoff=0
                 iboff=ja0-jb0
                 length=min(ja1,jb1)-ja0
              else
                 iaoff=jb0-ja0
                 iboff=0
                 length=min(ja1,jb1)-jb0
              endif
              do i=0,length
                 !          llf=llf+1
                 apsi_f(1,jaj+iaoff+i)=apsi_f(1,jaj+iaoff+i)+scprwp*bpsi_f(1,jbj+iboff+i)
                 apsi_f(2,jaj+iaoff+i)=apsi_f(2,jaj+iaoff+i)+scprwp*bpsi_f(2,jbj+iboff+i)
                 apsi_f(3,jaj+iaoff+i)=apsi_f(3,jaj+iaoff+i)+scprwp*bpsi_f(3,jbj+iboff+i)
                 apsi_f(4,jaj+iaoff+i)=apsi_f(4,jaj+iaoff+i)+scprwp*bpsi_f(4,jbj+iboff+i)
                 apsi_f(5,jaj+iaoff+i)=apsi_f(5,jaj+iaoff+i)+scprwp*bpsi_f(5,jbj+iboff+i)
                 apsi_f(6,jaj+iaoff+i)=apsi_f(6,jaj+iaoff+i)+scprwp*bpsi_f(6,jbj+iboff+i)
                 apsi_f(7,jaj+iaoff+i)=apsi_f(7,jaj+iaoff+i)+scprwp*bpsi_f(7,jbj+iboff+i)
              enddo
              ibseg=ibseg+1
              if (ibseg > mbseg_f) exit loop_jaf
           end do loop_jbf
        enddo loop_jaf
     end if
!$  endif
  !        write(*,*) 'waxpy,llc,llf',llc,llf
!$omp end parallel

  !  call system_clock(ncount2,ncount_rate,ncount_max)
  !  tel=dble(ncount2-ncount0)/dble(ncount_rate)
  !  write(97,'(a40,1x,e10.3,1x,f6.1)') 'waxpy:',tel
  !  close(97)


end subroutine waxpy





