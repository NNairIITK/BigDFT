subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
  ! calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
  implicit real(kind=8) (a-h,o-z)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

  scpr=0.d0
  do i=1,mvctr_c
     scpr=scpr+psi_c(i)**2
  enddo
  scpr1=0.d0
  scpr2=0.d0
  scpr3=0.d0
  scpr4=0.d0
  scpr5=0.d0
  scpr6=0.d0
  scpr7=0.d0
  do i=1,mvctr_f
     scpr1=scpr1+psi_f(1,i)**2
     scpr2=scpr2+psi_f(2,i)**2
     scpr3=scpr3+psi_f(3,i)**2
     scpr4=scpr4+psi_f(4,i)**2
     scpr5=scpr5+psi_f(5,i)**2
     scpr6=scpr6+psi_f(6,i)**2
     scpr7=scpr7+psi_f(7,i)**2
  enddo
  scpr=scpr+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

  return
END SUBROUTINE wnrm

subroutine wscal(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  ! multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
  implicit real(kind=8) (a-h,o-z)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal
  enddo
  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal
     psi_f(2,i)=psi_f(2,i)*scal
     psi_f(3,i)=psi_f(3,i)*scal
     psi_f(4,i)=psi_f(4,i)*scal
     psi_f(5,i)=psi_f(5,i)*scal
     psi_f(6,i)=psi_f(6,i)*scal
     psi_f(7,i)=psi_f(7,i)*scal
  enddo

  return
END SUBROUTINE wscal

subroutine wscalv(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  ! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
  implicit real(kind=8) (a-h,o-z)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f),scal(0:3)

  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal(0)           !  1 1 1
  enddo
  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
     psi_f(2,i)=psi_f(2,i)*scal(1)       !  1 2 1
     psi_f(3,i)=psi_f(3,i)*scal(2)       !  2 2 1
     psi_f(4,i)=psi_f(4,i)*scal(1)       !  1 1 2
     psi_f(5,i)=psi_f(5,i)*scal(2)       !  2 1 2
     psi_f(6,i)=psi_f(6,i)*scal(2)       !  1 2 2
     psi_f(7,i)=psi_f(7,i)*scal(3)       !  2 2 2
  enddo

  return
end subroutine wscalv

subroutine wzero(mvctr_c,mvctr_f,psi_c,psi_f)
  ! initializes a wavefunction to zero
  implicit real(kind=8) (a-h,o-z)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

  do i=1,mvctr_c
     psi_c(i)=0.d0
  enddo
  do i=1,mvctr_f
     psi_f(1,i)=0.d0
     psi_f(2,i)=0.d0
     psi_f(3,i)=0.d0
     psi_f(4,i)=0.d0
     psi_f(5,i)=0.d0
     psi_f(6,i)=0.d0
     psi_f(7,i)=0.d0
  enddo

  return
END SUBROUTINE wzero

!this function must be generalized for the linear scaling code
subroutine wpdot(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
  ! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
  ! Warning: the subroutine assumes that bpsi has only one segment along each line,
  ! whereas apsi can have several segments. This assumption is true if bpsi is a projector 
  ! To be more precise, it is assumed that the segments of bpsi are always contained inside
  ! the segments of apsi, no matter whether they are in the same line or not.
  implicit real(kind=8) (a-h,o-z)
  dimension keyav_c(maseg_c),keyag_c(2,maseg_c),keyav_f(maseg_f),keyag_f(2,maseg_f)
  dimension keybv_c(mbseg_c),keybg_c(2,mbseg_c),keybv_f(mbseg_f),keybg_f(2,mbseg_f)
  dimension apsi_c(mavctr_c),apsi_f(7,mavctr_f),bpsi_c(mbvctr_c),bpsi_f(7,mbvctr_f)

  !        llc=0
  scpr=0.d0
  ! coarse part
  ibseg=1
  do iaseg=1,maseg_c
     jaj=keyav_c(iaseg)
     ja0=keyag_c(1,iaseg)
     ja1=keyag_c(2,iaseg)

100  jb1=keybg_c(2,ibseg)
     if (jb1.lt.ja0) then
        ibseg=ibseg+1
        if (ibseg.gt.mbseg_c) goto 111
        goto 100
     endif
     jb0=keybg_c(1,ibseg)
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
     !           write(*,*) 'ja0,ja1,jb0,jb1',ja0,ja1,jb0,jb1,length
     !          write(*,'(5(a,i5))') 'C:from ',jaj+iaoff,' to ',jaj+iaoff+length,' and from ',jbj+iboff,' to ',jbj+iboff+length
     do i=0,length
        !          llc=llc+1
        scpr=scpr+apsi_c(jaj+iaoff+i)*bpsi_c(jbj+iboff+i) 
     enddo
  enddo
111 continue


  !        llf=0
  scpr1=0.d0
  scpr2=0.d0
  scpr3=0.d0
  scpr4=0.d0
  scpr5=0.d0
  scpr6=0.d0
  scpr7=0.d0
  ! fine part
  ibseg=1
  do iaseg=1,maseg_f
     jaj=keyav_f(iaseg)
     ja0=keyag_f(1,iaseg)
     ja1=keyag_f(2,iaseg)

200  jb1=keybg_f(2,ibseg)
     if (jb1.lt.ja0) then
        ibseg=ibseg+1
        if (ibseg.gt.mbseg_f) goto 222
        goto 200
     endif
     jb0=keybg_f(1,ibseg)
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
        scpr1=scpr1+apsi_f(1,jaj+iaoff+i)*bpsi_f(1,jbj+iboff+i) 
        scpr2=scpr2+apsi_f(2,jaj+iaoff+i)*bpsi_f(2,jbj+iboff+i) 
        scpr3=scpr3+apsi_f(3,jaj+iaoff+i)*bpsi_f(3,jbj+iboff+i) 
        scpr4=scpr4+apsi_f(4,jaj+iaoff+i)*bpsi_f(4,jbj+iboff+i) 
        scpr5=scpr5+apsi_f(5,jaj+iaoff+i)*bpsi_f(5,jbj+iboff+i) 
        scpr6=scpr6+apsi_f(6,jaj+iaoff+i)*bpsi_f(6,jbj+iboff+i) 
        scpr7=scpr7+apsi_f(7,jaj+iaoff+i)*bpsi_f(7,jbj+iboff+i) 
     enddo
  enddo
222 continue

  scpr=scpr+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
  !        write(*,*) 'llc,llf',llc,llf

  return
END SUBROUTINE wpdot

subroutine waxpy(  & 
     scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f, & 
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f)
  ! rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
  ! The update is only done in the localization region of apsi
  implicit real(kind=8) (a-h,o-z)
  dimension keyav_c(maseg_c),keyag_c(2,maseg_c),keyav_f(maseg_f),keyag_f(2,maseg_f)
  dimension keybv_c(mbseg_c),keybg_c(2,mbseg_c),keybv_f(mbseg_f),keybg_f(2,mbseg_f)
  dimension apsi_c(mavctr_c),apsi_f(7,mavctr_f),bpsi_c(mbvctr_c),bpsi_f(7,mbvctr_f)

  !        llc=0
  ! coarse part
  ibseg=1
  do iaseg=1,maseg_c
     jaj=keyav_c(iaseg)
     ja0=keyag_c(1,iaseg)
     ja1=keyag_c(2,iaseg)

100  jb1=keybg_c(2,ibseg)
     if (jb1.lt.ja0) then
        ibseg=ibseg+1
        if (ibseg.gt.mbseg_c) goto 111
        goto 100
     endif
     jb0=keybg_c(1,ibseg)
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
        apsi_c(jaj+iaoff+i)=apsi_c(jaj+iaoff+i)+scpr*bpsi_c(jbj+iboff+i) 
     enddo
  enddo
111 continue

  !        llf=0
  ! fine part
  ibseg=1
  do iaseg=1,maseg_f
     jaj=keyav_f(iaseg)
     ja0=keyag_f(1,iaseg)
     ja1=keyag_f(2,iaseg)

200  jb1=keybg_f(2,ibseg)
     if (jb1.lt.ja0) then
        ibseg=ibseg+1
        if (ibseg.gt.mbseg_f) goto 222
        goto 200
     endif
     jb0=keybg_f(1,ibseg)
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
        apsi_f(1,jaj+iaoff+i)=apsi_f(1,jaj+iaoff+i)+scpr*bpsi_f(1,jbj+iboff+i)
        apsi_f(2,jaj+iaoff+i)=apsi_f(2,jaj+iaoff+i)+scpr*bpsi_f(2,jbj+iboff+i)
        apsi_f(3,jaj+iaoff+i)=apsi_f(3,jaj+iaoff+i)+scpr*bpsi_f(3,jbj+iboff+i)
        apsi_f(4,jaj+iaoff+i)=apsi_f(4,jaj+iaoff+i)+scpr*bpsi_f(4,jbj+iboff+i)
        apsi_f(5,jaj+iaoff+i)=apsi_f(5,jaj+iaoff+i)+scpr*bpsi_f(5,jbj+iboff+i)
        apsi_f(6,jaj+iaoff+i)=apsi_f(6,jaj+iaoff+i)+scpr*bpsi_f(6,jbj+iboff+i)
        apsi_f(7,jaj+iaoff+i)=apsi_f(7,jaj+iaoff+i)+scpr*bpsi_f(7,jbj+iboff+i)
     enddo
  enddo
222 continue
  !        write(*,*) 'waxpy,llc,llf',llc,llf

  return
END SUBROUTINE waxpy





