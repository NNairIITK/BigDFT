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

subroutine wpdot(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
  ! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
  ! Warning: the subroutine assumes that bpsi has only one segment along each line,
  ! whereas apsi can have several segments. This assumption is true if bpsi is a projector 
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

subroutine uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psig)
  ! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
  dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)

  call razero(8*(nu1-nl1+1)*(nu2-nl2+1)*(nu3-nl3+1),psig)

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
        ii=ii+1
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
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
        psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
        psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
        psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
        psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
        psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
        psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
        psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
     enddo
  enddo

END SUBROUTINE uncompress


subroutine compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     psig,psi_c,psi_f)
  ! Compresses a psig wavefunction into psi_c,psi_f form
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
  dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)

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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
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
        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
     enddo
  enddo

END SUBROUTINE compress


subroutine uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_fc,psig_f)
  ! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f),scal(0:3)
  dimension psig_c(0:n1,0:n2,0:n3)
  dimension psig_fc(0:n1,0:n2,0:n3,3),psig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)

  call razero((n1+1)*(n2+1)*(n3+1),psig_c)
  call razero(3*(n1+1)*(n2+1)*(n3+1),psig_fc)
  call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),psig_f)

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
        psig_fc(i,i2,i3,1)=psig_f(1,i,i2,i3)
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        psig_fc(i,i2,i3,2)=psig_f(2,i,i2,i3)
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        psig_fc(i,i2,i3,3)=psig_f(4,i,i2,i3)
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
  implicit real(kind=8) (a-h,o-z)
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

subroutine convolut_magic_n(n1,n2,n3,x,y)
  ! Applies the magic filter matrix ( no transposition) ; data set grows
  ! The input array x is not overwritten
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
  dimension x(*)
  dimension y(*)

  !  (i1,i2*i3) -> (i2*i3,I1)
  ndat=(n2+1)*(n3+1)
  call convrot_grow(n1,ndat,x,y)
  !  (i2,i3*I1) -> (i3*I1,I2)
  ndat=(n3+1)*(n1+1+lupfil-lowfil)
  call convrot_grow(n2,ndat,y,x)
  !  (i3,I1*I2) -> (iI*I2,I3)
  ndat=(n1+1+lupfil-lowfil)*(n2+1+lupfil-lowfil)
  call convrot_grow(n3,ndat,x,y)


END SUBROUTINE convolut_magic_n


subroutine convolut_magic_t(n1,n2,n3,x,y)
  ! Applies the magic filter matrix transposed ; data set shrinks
  ! The input array x is overwritten
  implicit real(kind=8) (a-h,o-z)
  parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
  dimension x(*),y(*)

  !  (I1,I2*I3) -> (I2*I3,i1)
  ndat=(n2+1+lupfil-lowfil)*(n3+1+lupfil-lowfil)
  call convrot_shrink(n1,ndat,x,y)
  !  (I2,I3*i1) -> (I3*i1,i2)
  ndat=(n3+1+lupfil-lowfil)*(n1+1)
  call convrot_shrink(n2,ndat,y,x)
  !  (I3,i1*i2) -> (i1*i2,i3)
  ndat=(n1+1)*(n2+1)
  call convrot_shrink(n3,ndat,x,y)

END SUBROUTINE convolut_magic_t


subroutine synthese_grow(n1,n2,n3,ww,x,y)
  ! A synthesis wavelet transformation where the size of the data is allowed to grow
  ! The input array x is not overwritten
  implicit real(kind=8) (a-h,o-z)
  dimension x(0:n1,2,0:n2,2,0:n3,2)
  dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
  dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

  ! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,ww)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_grow(n3,nt,ww,y)

END SUBROUTINE synthese_grow



subroutine analyse_shrink(n1,n2,n3,ww,y,x)
  ! A analysis wavelet transformation where the size of the data is forced to shrink
  ! The input array y is overwritten
  implicit real(kind=8) (a-h,o-z)
  dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
  dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
  dimension x(0:n1,2,0:n2,2,0:n3,2)

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+16)*(2*n3+16)
  call  ana_rot_shrink(n1,nt,y,ww)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+16)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_shrink(n3,nt,y,x)

  return
END SUBROUTINE analyse_shrink

subroutine plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,rx,ry,rz,psi)
  implicit real(kind=8) (a-h,o-z)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f)
  real(kind=8), allocatable :: psifscf(:),psir(:),psig(:,:,:,:,:,:)

  allocate(psig(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
  call memocc(i_stat,product(shape(psig))*kind(psig),'psig','plot_wf')
  allocate(psifscf((2*n1+31)*(2*n2+31)*(2*n3+16)),stat=i_stat)
  call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','plot_wf')
  allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','plot_wf')

  call uncompress(n1,n2,n3,0,n1,0,n2,0,n3, & 
       nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       psi(1),psi(nvctr_c+1),psig)
  call synthese_grow(n1,n2,n3,psir,psig,psifscf)  !psir=ww(((2*n1+16)*(2*n2+16)*(2*n1+2))
  call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

  call plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,psir)

  i_all=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=i_stat)
  call memocc(i_stat,i_all,'psig','plot_wf')
  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf','plot_wf')
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','plot_wf')
  return
END SUBROUTINE plot_wf

subroutine plot_pot(rx,ry,rz,hgrid,n1,n2,n3,iounit,pot)
  implicit real(kind=8) (a-h,o-z)
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
  implicit real(kind=8) (a-h,o-z)
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
