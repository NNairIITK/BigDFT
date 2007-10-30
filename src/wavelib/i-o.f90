subroutine reformatonewave(iproc, hgrid_old, n1_old, n2_old, n3_old, &
     & center_old, psigold, hgrid, nvctr_c, nvctr_f, n1, n2, n3, center, nseg_c, nseg_f, &
     & keyg, keyv, psifscf, psi)
  implicit real(kind=8) (a-h,o-z)
  logical cif1,cif2,cif3
  dimension xya(-1:1,-1:1),xa(-1:1)
  dimension :: center(3), center_old(3)
  dimension :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
  dimension :: psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2), psi(nvctr_c + 7 * nvctr_f)
  dimension :: psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

  allocatable :: psifscfold(:,:,:),psifscfoex(:,:,:),psig(:,:,:,:,:,:),ww(:)

  allocate(psifscfold(-7:2*n1_old+8,-7:2*n2_old+8,-7:2*n3_old+8),stat=i_stat)
  call memocc(i_stat,product(shape(psifscfold))*kind(psifscfold),'psifscfold','reformatonewave')
  allocate(psifscfoex(-8:2*n1_old+9,-8:2*n2_old+9,-8:2*n3_old+9),stat=i_stat)
  call memocc(i_stat,product(shape(psifscfoex))*kind(psifscfoex),'psifscfoex','reformatonewave')

  ! calculate fine scaling functions, psifscfoex=wwold((2*n1_old+16)*(2*n2_old+16)*(2*n3_old+16))
  call synthese_grow(n1_old,n2_old,n3_old,psifscfoex,psigold,psifscfold) 

  do i3=-7,2*n3_old+8
     do i2=-7,2*n2_old+8
        i1=-8
        psifscfoex(i1,i2,i3)=0.d0
        do i1=-7,2*n1_old+8
           psifscfoex(i1,i2,i3)=psifscfold(i1,i2,i3)
        enddo
        i1=2*n1_old+9
        psifscfoex(i1,i2,i3)=0.d0
     enddo
  enddo

  i3=-8
  do i2=-8,2*n2_old+9
     do i1=-8,2*n1_old+9
        psifscfoex(i1,i2,i3)=0.d0
     enddo
  enddo
  i3=2*n3_old+9
  do i2=-8,2*n2_old+9
     do i1=-8,2*n1_old+9
        psifscfoex(i1,i2,i3)=0.d0
     enddo
  enddo

  i2=-8
  do i3=-8,2*n3_old+9
     do i1=-8,2*n1_old+9
        psifscfoex(i1,i2,i3)=0.d0
     enddo
  enddo
  i2=2*n2_old+9
  do i3=-8,2*n3_old+9
     do i1=-8,2*n1_old+9
        psifscfoex(i1,i2,i3)=0.d0
     enddo
  enddo

  ! transform to new structure    
  dx=center(1)-center_old(1)
  dy=center(2)-center_old(2)
  dz=center(3)-center_old(3)
  !   write(*,*) 'dxyz',dx,dy,dz
  hgridh=.5d0*hgrid
  hgridh_old=.5d0*hgrid_old
  call razero((2*n1+16)*(2*n2+16)*(2*n3+16),psifscf)
  do i3=-7,2*n3+8
     z=real(i3,kind=8)*hgridh
     j3=nint((z-dz)/hgridh_old)
     cif3=(j3.ge.-7 .and. j3.le.2*n3_old+8)
     do i2=-7,2*n2+8
        y=real(i2,kind=8)*hgridh
        j2=nint((y-dy)/hgridh_old)
        cif2=(j2.ge.-7 .and. j2.le.2*n2_old+8)
        do i1=-7,2*n1+8
           x=real(i1,kind=8)*hgridh
           j1=nint((x-dx)/hgridh_old)
           cif1=(j1.ge.-7 .and. j1.le.2*n1_old+8)

           !        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfold(j1,j2,j3)
           !        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfoex(j1,j2,j3)

           if (cif1 .and. cif2 .and. cif3) then 
              zr = ((z-dz)-real(j3,kind=8)*hgridh_old)/hgridh_old
              do l2=-1,1
                 do l1=-1,1
                    ym1=psifscfoex(j1+l1,j2+l2,j3-1)
                    y00=psifscfoex(j1+l1,j2+l2,j3  )
                    yp1=psifscfoex(j1+l1,j2+l2,j3+1)
                    xya(l1,l2)=ym1 + (1.d0 + zr)*(y00 - ym1 + zr*(.5d0*ym1 - y00  + .5d0*yp1))
                 enddo
              enddo

              yr = ((y-dy)-real(j2,kind=8)*hgridh_old)/hgridh_old
              do l1=-1,1
                 ym1=xya(l1,-1)
                 y00=xya(l1,0)
                 yp1=xya(l1,1)
                 xa(l1)=ym1 + (1.d0 + yr)*(y00 - ym1 + yr*(.5d0*ym1 - y00  + .5d0*yp1))
              enddo

              xr = ((x-dx)-real(j1,kind=8)*hgridh_old)/hgridh_old
              ym1=xa(-1)
              y00=xa(0)
              yp1=xa(1)
              psifscf(i1,i2,i3)=ym1 + (1.d0 + xr)*(y00 - ym1 + xr*(.5d0*ym1 - y00  + .5d0*yp1))

           endif

        enddo
     enddo
  enddo

  i_all=-product(shape(psifscfold))*kind(psifscfold)
  deallocate(psifscfold,stat=i_stat)
  call memocc(i_stat,i_all,'psifscfold','reformatonewave')
  i_all=-product(shape(psifscfoex))*kind(psifscfoex)
  deallocate(psifscfoex,stat=i_stat)
  call memocc(i_stat,i_all,'psifscfoex','reformatonewave')
  allocate(psig(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
  call memocc(i_stat,product(shape(psig))*kind(psig),'psig','reformatonewave')
  allocate(ww((2*n1+16)*(2*n2+16)*(2*n3+16)),stat=i_stat)
  call memocc(i_stat,product(shape(ww))*kind(ww),'ww','reformatonewave')

  call analyse_shrink(n1,n2,n3,ww,psifscf,psig)
  call compress(n1,n2,n3,0,n1,0,n2,0,n3,  &
       nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       psig,psi(1),psi(nvctr_c+1))

  i_all=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=i_stat)
  call memocc(i_stat,i_all,'psig','reformatonewave')
  i_all=-product(shape(ww))*kind(ww)
  deallocate(ww,stat=i_stat)
  call memocc(i_stat,i_all,'ww','reformatonewave')
END SUBROUTINE reformatonewave

subroutine readonewave(unitwf, useFormattedInput, iorb,iproc,n1,n2,n3, &
     & hgrid,center,nseg_c,nseg_f, nvctr_c,nvctr_f,keyg,keyv,psi,eval,psifscf)
  implicit real(kind=8) (a-h,o-z)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f)
  dimension center(3),center_old(3)
  dimension :: psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
  allocatable :: psigold(:,:,:,:,:,:)
  integer :: unitwf
  logical :: useFormattedInput

  if (useFormattedInput) then
     read(unitwf,*) iorb_old,eval
  else
     read(unitwf) iorb_old,eval
  end if
  if (iorb_old.ne.iorb) stop 'readonewave'
  if (useFormattedInput) then
     read(unitwf,*) hgrid_old
     read(unitwf,*) n1_old,n2_old,n3_old
     read(unitwf,*) (center_old(j),j=1,3)
  else
     read(unitwf) hgrid_old
     read(unitwf) n1_old,n2_old,n3_old
     read(unitwf) (center_old(j),j=1,3)
  end if
  write(*,'(1x,i2,6(1x,e14.7))') iproc,(center(j),j=1,3),(center_old(j),j=1,3)
  if (useFormattedInput) then
     read(unitwf,*) nvctr_c_old, nvctr_f_old
  else
     read(unitwf) nvctr_c_old, nvctr_f_old
  end if

  !           write(*,*) iorb,' hgrid_old,hgrid ',hgrid_old,hgrid
  !           write(*,*) iorb,' nvctr_c_old,nvctr_c ',nvctr_c_old,nvctr_c
  !           write(*,*) iorb,' nvctr_f_old,nvctr_f ',nvctr_f_old,nvctr_f
  !           write(*,*) iorb,' n1_old,n1 ',n1_old,n1
  !           write(*,*) iorb,' n2_old,n2 ',n2_old,n2
  !           write(*,*) iorb,' n3_old,n3 ',n3_old,n3

  if (hgrid_old.eq. hgrid .and. nvctr_c_old.eq.nvctr_c .and. nvctr_f_old.eq.nvctr_f  & 
       .and. n1_old.eq.n1  .and. n2_old.eq.n2 .and. n3_old.eq.n3 ) then

     write(*,*) 'wavefunction ',iorb,' needs NO reformatting on processor',iproc
     do j=1,nvctr_c_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,tt
        else
           read(unitwf) i1,i2,i3,tt
        end if
        psi(j)=tt
     enddo
     do j=1,7*nvctr_f_old-6,7
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        else
           read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        end if
        psi(nvctr_c+j+0)=t1
        psi(nvctr_c+j+1)=t2
        psi(nvctr_c+j+2)=t3
        psi(nvctr_c+j+3)=t4
        psi(nvctr_c+j+4)=t5
        psi(nvctr_c+j+5)=t6
        psi(nvctr_c+j+6)=t7
     enddo

  else
     write(*,*) 'wavefunction ',iorb,' needs reformatting on processor',iproc
     if (hgrid_old.ne.hgrid) write(*,*) 'because hgrid_old >< hgrid',hgrid_old,hgrid
     if (nvctr_c_old.ne.nvctr_c) write(*,*) 'because nvctr_c_old >< nvctr_c',nvctr_c_old,nvctr_c
     if (nvctr_f_old.ne.nvctr_f) write(*,*) 'because nvctr_f_old >< nvctr_f',nvctr_f_old,nvctr_f
     if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 ) &
          write(*,*) 'because cell size has changed',n1_old,n1  , n2_old,n2 , n3_old,n3

     allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2),stat=i_stat)
     call memocc(i_stat,product(shape(psigold))*kind(psigold),'psigold','readonewave')

     call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)
     do iel=1,nvctr_c_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,tt
        else
           read(unitwf) i1,i2,i3,tt
        end if
        psigold(i1,1,i2,1,i3,1)=tt
     enddo
     do iel=1,nvctr_f_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        else
           read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        end if
        psigold(i1,2,i2,1,i3,1)=t1
        psigold(i1,1,i2,2,i3,1)=t2
        psigold(i1,2,i2,2,i3,1)=t3
        psigold(i1,1,i2,1,i3,2)=t4
        psigold(i1,2,i2,1,i3,2)=t5
        psigold(i1,1,i2,2,i3,2)=t6
        psigold(i1,2,i2,2,i3,2)=t7
     enddo

     ! I put nat = 1 here, since only one position is saved in wavefunction files.
     call reformatonewave(iproc, hgrid_old, n1_old, n2_old, n3_old, &
          & center_old, psigold, hgrid, nvctr_c, nvctr_f, n1, n2, n3, center, nseg_c, nseg_f, &
          & keyg, keyv, psifscf, psi)

     i_all=-product(shape(psigold))*kind(psigold)
     deallocate(psigold,stat=i_stat)
     call memocc(i_stat,i_all,'psigold','readonewave')

  endif
END SUBROUTINE readonewave


subroutine writeonewave(unitwf, useFormattedOutput, iorb,n1,n2,n3,hgrid,center,  & 
     nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f, & 
     psi_c,psi_f,norb,eval)
  implicit real(kind=8) (a-h,o-z)
  logical :: useFormattedOutput
  integer :: unitwf
  dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
  dimension psi_c(nvctr_c),psi_f(7,nvctr_f),center(3),eval(norb)


  if (useFormattedOutput) then
     write(unitwf,*) iorb,eval(iorb)
     write(unitwf,*) hgrid
     write(unitwf,*) n1,n2,n3
     write(unitwf,'(3(1x,e24.17))') (center(j),j=1,3)
     write(unitwf,*) nvctr_c, nvctr_f
  else
     write(unitwf) iorb,eval(iorb)
     write(unitwf) hgrid
     write(unitwf) n1,n2,n3
     write(unitwf) (center(j),j=1,3)
     write(unitwf) nvctr_c, nvctr_f
  end if

  ! coarse part
  do iseg=1,nseg_c
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
        tt=psi_c(i-i0+jj) 
        if (useFormattedOutput) then
           write(unitwf,'(3(i4),1x,e19.12)') i,i2,i3,tt
        else
           write(unitwf) i,i2,i3,tt
        end if
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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
        t1=psi_f(1,i-i0+jj)
        t2=psi_f(2,i-i0+jj)
        t3=psi_f(3,i-i0+jj)
        t4=psi_f(4,i-i0+jj)
        t5=psi_f(5,i-i0+jj)
        t6=psi_f(6,i-i0+jj)
        t7=psi_f(7,i-i0+jj)
        if (useFormattedOutput) then
           write(unitwf,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
        else
           write(unitwf) i,i2,i3,t1,t2,t3,t4,t5,t6,t7
        end if
     enddo
  enddo

  write(*,'(1x,i0,a)') iorb,'th wavefunction written'


END SUBROUTINE writeonewave
