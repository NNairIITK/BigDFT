subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
  ! calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
  implicit real*8 (a-h,o-z)
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
  implicit real*8 (a-h,o-z)
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
  implicit real*8 (a-h,o-z)
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
  implicit real*8 (a-h,o-z)
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

subroutine wpdot(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
  ! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
  ! Warning: the subroutine assumes that bpsi has only one segment along each line,
  ! whereas apsi can have several segments. This assumption is true if bpsi is a projector 
  implicit real*8 (a-h,o-z)
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
  implicit real*8 (a-h,o-z)
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



!subroutine compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
!     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
!     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
!     psig,psi_c,psi_f)
!  ! Compresses a psig wavefunction into psi_c,psi_f form
!  implicit real*8 (a-h,o-z)
!  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
!  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
!  dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)
!
!  ! coarse part
!  do iseg=1,mseg_c
!     jj=keyv_c(iseg)
!     j0=keyg_c(1,iseg)
!     j1=keyg_c(2,iseg)
!     ii=j0-1
!     i3=ii/((n1+1)*(n2+1))
!     ii=ii-i3*(n1+1)*(n2+1)
!     i2=ii/(n1+1)
!     i0=ii-i2*(n1+1)
!     i1=i0+j1-j0
!     do i=i0,i1
!        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
!     enddo
!  enddo
!
!  ! fine part
!  do iseg=1,mseg_f
!     jj=keyv_f(iseg)
!     j0=keyg_f(1,iseg)
!     j1=keyg_f(2,iseg)
!     ii=j0-1
!     i3=ii/((n1+1)*(n2+1))
!     ii=ii-i3*(n1+1)*(n2+1)
!     i2=ii/(n1+1)
!     i0=ii-i2*(n1+1)
!     i1=i0+j1-j0
!     do i=i0,i1
!        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
!        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
!        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
!        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
!        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
!        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
!        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
!     enddo
!  enddo
!
!END SUBROUTINE compress

subroutine uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_f)
  ! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
  implicit real*8 (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f),scal(0:3)
  dimension psig_c(0:n1,0:n2,0:n3)
  dimension psig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
		
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
		
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
		
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo


end subroutine uncompress_forstandard_short

subroutine uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_f,&
	 x_f1,x_f2,x_f3)
  ! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
  implicit real*8 (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f),scal(0:3)
  dimension psig_c(0:n1,0:n2,0:n3)
  dimension psig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
	real*8::x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
	real*8::x_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3)
	real*8::x_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2)

  

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
		x_f1(i,i2,i3)=psig_f(1,i,i2,i3)
		
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
		x_f2(i2,i,i3)=psig_f(2,i,i2,i3)
		
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
		x_f3(i3,i,i2)=psig_f(4,i,i2,i3)
		
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo


end subroutine uncompress_forstandard


subroutine compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psig_c,psig_f,psi_c,psi_f)
  ! Compresses a psig wavefunction into psi_c,psi_f form
  implicit real*8 (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f),scal(0:3)
  dimension psig_c(0:n1,0:n2,0:n3)
  dimension psig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)

  ! coarse part
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_c(i-i0+jj)=psig_c(i,i2,i3)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)*scal(1)
        psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)*scal(1)
        psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)*scal(2)
        psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)*scal(1)
        psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)*scal(2)
        psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)*scal(2)
        psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)*scal(3)
     enddo
  enddo

end subroutine compress_forstandard


subroutine plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,rx,ry,rz,psi,&
 ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
 nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  implicit real*8 (a-h,o-z)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f)
  !	for grow:
  integer ibyz_c(2,0:n2,0:n3)
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
  !	for real space:
  integer,intent(in):: ibyyzz_r(2,2*n2+31,2*n3+31)
  real*8 scal(0:3)
  real*8, allocatable :: psir(:)

  real*8,allocatable,dimension(:,:,:)::x_c!input 
  real*8,allocatable::x_f(:,:,:,:)
  real*8,allocatable,dimension(:):: w1,w2

  nw1=max((n3+1)*(2*n1+31)*(2*n2+31),&   		
       (n1+1)*(2*n2+31)*(2*n3+31),&
       2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))

  nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
       (n1+1)*(n2+1)*(2*n3+31))

  do i=0,3
     scal(i)=1.d0
  enddo

  allocate(x_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(x_c))*kind(x_c),'x_c','plot_wf')
  
  allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)! work
  call memocc(i_stat,product(shape(x_f))*kind(x_f),'x_f','plot_wf')
  allocate(w1(nw1),stat=i_stat)
  call memocc(i_stat,product(shape(w1))*kind(w1),'w1','plot_wf')
  allocate(w2(nw2),stat=i_stat) ! work
  call memocc(i_stat,product(shape(w2))*kind(w2),'w2','plot_wf')
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','plot_wf')

  call razero((n1+1)*(n2+1)*(n3+1),x_c)
  call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f)

        call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
             nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
             scal,psi(1),psi(nvctr_c+1),x_c,x_f)

        call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
             psir,ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)

  call plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,psir)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','plot_wf')

  i_all=-product(shape(x_c))*kind(x_c)
  deallocate(x_c,stat=i_stat)
  call memocc(i_stat,i_all,'x_c','plot_wf')

  i_all=-product(shape(x_f))*kind(x_f)
  deallocate(x_f,stat=i_stat)
  call memocc(i_stat,i_all,'x_f','plot_wf')

  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1','plot_wf')

  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2','plot_wf')
  return
END SUBROUTINE plot_wf

subroutine plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,pot)
  implicit real*8 (a-h,o-z)
  dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

  hgridh=.5d0*hgrid
  open(iounit) 
  open(iounit+1) 
  open(iounit+2) 

  i3=nint(rz/hgridh)
  i2=nint(ry/hgridh)
  write(*,*) 'plot_p, i2,i3,n2,n3 ',i2,i3,n2,n3
  do i1=-14,2*n1+16
     write(iounit,*) i1*hgridh,pot(i1,i2,i3)
  enddo

  i1=nint(rx/hgridh)
  i2=nint(ry/hgridh)
  write(*,*) 'plot_p, i1,i2 ',i1,i2
  do i3=-14,2*n3+16
     write(iounit+1,*) i3*hgridh,pot(i1,i2,i3)
  enddo

  i1=nint(rx/hgridh)
  i3=nint(rz/hgridh)
  write(*,*) 'plot_p, i1,i3 ',i1,i3
  do i2=-14,2*n2+16
     write(iounit+2,*) i2*hgridh,pot(i1,i2,i3)
  enddo

  close(iounit) 
  close(iounit+1) 
  close(iounit+2) 

  return
END SUBROUTINE plot_pot



subroutine plot_psifscf(iunit,hgrid,n1,n2,n3,psifscf)
  implicit real*8 (a-h,o-z)
  dimension psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

  hgridh=.5d0*hgrid

  ! along x-axis
  i3=n3
  i2=n2
  do i1=-7,2*n1+8
     write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
  enddo

  ! 111 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=i ; i2=i ; i3=i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
  enddo

  ! 1-1-1 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=i ; i2=-i ; i3=-i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
  enddo

  ! -11-1 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=-i ; i2=i ; i3=-i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
  enddo

  ! -1-11 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=-i ; i2=-i ; i3=i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
  enddo

  return
END SUBROUTINE plot_psifscf
