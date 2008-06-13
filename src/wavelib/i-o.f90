subroutine reformatonewave(iproc,displ,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,nat,&
     & rxyz_old,psigold,hx,hy,hz,nvctr_c,nvctr_f,n1,n2,n3,rxyz,nseg_c,nseg_f,&
     & keyg,keyv,psifscf,psi)
  use module_base
  implicit none
  integer, intent(in) :: iproc,n1_old,n2_old,n3_old,nat,nvctr_c,nvctr_f,n1,n2,n3,nseg_c,nseg_f
  real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(gp), dimension(3,nat), intent(in) :: rxyz_old,rxyz
  real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(out) :: psifscf
  real(wp), dimension(nvctr_c + 7 * nvctr_f), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='reformatonewave'
  logical :: cif1,cif2,cif3
  integer :: i_stat,i_all,i1,i2,i3,j1,j2,j3,l1,l2,iat
  real(gp) :: hxh,hyh,hzh,hxh_old,hyh_old,hzh_old,x,y,z,dx,dy,dz,xold,yold,zold
  real(wp) :: zr,yr,xr,y01,ym1,y00,yp1
  real(wp), dimension(-1:1,-1:1) :: xya
  real(wp), dimension(-1:1) :: xa
  real(wp), dimension(:), allocatable :: ww,wwold
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psig
  real(wp), dimension(:,:,:), allocatable :: psifscfold

  allocate(psifscfold(-7:2*n1_old+8,-7:2*n2_old+8,-7:2*n3_old+8+ndebug),stat=i_stat)
  call memocc(i_stat,psifscfold,'psifscfold',subname)
  allocate(wwold((2*n1_old+16)*(2*n2_old+16)*(2*n3_old+16)+ndebug),stat=i_stat)
  call memocc(i_stat,wwold,'wwold',subname)

  call synthese_grow(n1_old,n2_old,n3_old,wwold,psigold,psifscfold) 

  i_all=-product(shape(wwold))*kind(wwold)
  deallocate(wwold,stat=i_stat)
  call memocc(i_stat,i_all,'wwold',subname)

  !write(100+iproc,*) 'norm psifscfold ',dnrm2((2*n1_old+16)*(2*n2_old+16)*(2*n3_old+16),&
  !     psifscfold,1)
  
  !write(*,*) iproc,' displ ',displ
  if (hx == hx_old .and. hy == hy_old .and. hz == hz_old .and. &
       n1_old==n1 .and. n2_old==n2 .and. n3_old==n3 .and. &
       displ<= 1.d-2) then
     !if (iproc==0) write(*,*) iproc,' orbital just copied'
     do i3=-7,2*n3+8
        do i2=-7,2*n2+8
           do i1=-7,2*n1+8
              psifscf(i1,i2,i3)=psifscfold(i1,i2,i3)
           enddo
        enddo
     enddo
     
  else

     dx=0.0_gp
     dy=0.0_gp 
     dz=0.0_gp
     do iat=1,nat ! Calculate average shift
        dx=dx+(rxyz(1,iat)-rxyz_old(1,iat))
        dy=dy+(rxyz(2,iat)-rxyz_old(2,iat))
        dz=dz+(rxyz(3,iat)-rxyz_old(3,iat))
     enddo
     dx=dx/real(nat,gp)
     dy=dy/real(nat,gp)
     dz=dz/real(nat,gp)
     
     ! transform to new structure    
     !if (iproc==0) write(*,*) iproc,' orbital fully transformed'
     hxh=.5_gp*hx
     hxh_old=.5_gp*hx_old
     hyh=.5_gp*hy
     hyh_old=.5_gp*hy_old
     hzh=.5_gp*hz
     hzh_old=.5_gp*hz_old

     call razero((2*n1+16)*(2*n2+16)*(2*n3+16),psifscf)

     do i3=-7,2*n3+8
        z=real(i3,gp)*hzh
        do i2=-7,2*n2+8
           y=real(i2,gp)*hyh
           do i1=-7,2*n1+8
              x=real(i1,gp)*hxh
!! The lines below might be reactivated in connection with a better interpolation scheme
!           if (nat.le.10) then ! calculate individual shifts for each atom
!             dx=0.d0 ; dy=0.d0 ; dz=0.d0 ; w=0.d0
!             do iat=1,nat ! Warning: this unoptimized loop can become very expensive for a large NAT
!               rbond=1.d0
!               arg=(x-rxyz(1,iat))**2+(x-rxyz(2,iat))**2+(x-rxyz(3,iat))**2
!               arg=-.5d0*(arg/rbond**2)
!               xp=exp(arg)
!               w=w+xp
!               dx=dx+xp*(rxyz(1,iat)-rxyz_old(1,iat))
!               dy=dy+xp*(rxyz(2,iat)-rxyz_old(2,iat))
!               dz=dz+xp*(rxyz(3,iat)-rxyz_old(3,iat))
!             enddo
!             w=w+1.d-10
!             dx=dx*(1.d0/w) ; dy=dy*(1.d0/w) ; dz=dz*(1.d0/w)
!           endif
              xold=x-dx 
              yold=y-dy
              zold=z-dz
!if (i1.eq.n1 .and. i2.eq.n2 .and. i3.eq.n3  .or.  & 
!    i1.eq.-7 .and. i2.eq.-7 .and. i3.eq.-7  .or.  & 
!    i1.eq.0 .and. i2.eq.0 .and. i3.eq.0  .or.  & 
!    i1.eq.n1/2 .and. i2.eq.n2/2 .and. i3.eq.n3/2 .or.  &
!    i1.eq.3*n1/2 .and. i2.eq.3*n2/2 .and. i3.eq.3*n3/2 .or.  &
!    i1.eq.n1/4 .and. i2.eq.n2/4 .and. i3.eq.n3/4  ) then
!write(iproc+200,'(a,6(e12.4))') 'new ',x,y,z,dx,dy,dz
!write(iproc+200,'(a,3(e12.4))') 'old ',xold,yold,zold
!endif
              j1=nint((xold)/hxh_old)
              cif1=(j1 >= -6 .and. j1 <= 2*n1_old+7)
              j2=nint((yold)/hyh_old)
              cif2=(j2 >= -6 .and. j2 <= 2*n2_old+7)
              j3=nint((zold)/hzh_old)
              cif3=(j3 >= -6 .and. j3 <= 2*n3_old+7)
              
           !        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfold(j1,j2,j3)
           !        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfoex(j1,j2,j3)

              if (cif1 .and. cif2 .and. cif3) then 
                 zr =real(((z-dz)-real(j3,gp)*hzh_old)/hzh_old,wp)
                 do l2=-1,1
                    do l1=-1,1
                       ym1=psifscfold(j1+l1,j2+l2,j3-1)
                       y00=psifscfold(j1+l1,j2+l2,j3  )
                       yp1=psifscfold(j1+l1,j2+l2,j3+1)
                       xya(l1,l2)=ym1 + &
                            (1.0_wp + zr)*(y00 - ym1 + zr*(.5_wp*ym1 - y00  + .5_wp*yp1))
                    enddo
                 enddo

                 yr = real(((y-dy)-real(j2,gp)*hyh_old)/hyh_old,wp)
                 do l1=-1,1
                    ym1=xya(l1,-1)
                    y00=xya(l1,0)
                    yp1=xya(l1,1)
                    xa(l1)=ym1 + &
                         (1.0_wp + yr)*(y00 - ym1 + yr*(.5_wp*ym1 - y00  + .5_wp*yp1))
                 enddo

                 xr = real(((x-dx)-real(j1,gp)*hxh_old)/hxh_old,wp)
                 ym1=xa(-1)
                 y00=xa(0)
                 yp1=xa(1)
                 psifscf(i1,i2,i3)=ym1 + &
                      (1.0_wp + xr)*(y00 - ym1 + xr*(.5_wp*ym1 - y00  + .5_wp*yp1))

              endif

           enddo
        enddo
     enddo
  endif

!write(100+iproc,*) 'norm of psifscf ',dnrm2((2*n1+16)*(2*n2+16)*(2*n3+16),psifscf,1)

  i_all=-product(shape(psifscfold))*kind(psifscfold)
  deallocate(psifscfold,stat=i_stat)
  call memocc(i_stat,i_all,'psifscfold',subname)
  allocate(psig(0:n1,2,0:n2,2,0:n3,2+ndebug),stat=i_stat)
  call memocc(i_stat,psig,'psig',subname)
  allocate(ww((2*n1+16)*(2*n2+16)*(2*n3+16)+ndebug),stat=i_stat)
  call memocc(i_stat,ww,'ww',subname)

  call analyse_shrink(n1,n2,n3,ww,psifscf,psig)

  !write(100+iproc,*) 'norm new psig ',dnrm2(8*(n1+1)*(n2+1)*(n3+1),psig,1)
  call compress(n1,n2,n3,0,n1,0,n2,0,n3,  &
       nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       psig,psi(1),psi(nvctr_c+1))

  !write(100+iproc,*) 'norm of reformatted psi ',dnrm2(nvctr_c+7*nvctr_f,psi,1)

  i_all=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=i_stat)
  call memocc(i_stat,i_all,'psig',subname)
  i_all=-product(shape(ww))*kind(ww)
  deallocate(ww,stat=i_stat)
  call memocc(i_stat,i_all,'ww',subname)
end subroutine reformatonewave

subroutine readonewave(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
     & hx,hy,hz,nat,rxyz_old,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval,psifscf)
  use module_base
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3,nat,nseg_c,nseg_f,nvctr_c,nvctr_f
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), intent(out) :: eval
  real(gp), dimension(3,nat), intent(out) :: rxyz_old
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(out) :: psi
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(out) :: psifscf
  !local variables
  character(len=*), parameter :: subname='readonewave'
  integer :: iorb_old,n1_old,n2_old,n3_old,iat,j,iel,nvctr_c_old,nvctr_f_old,i_stat,i_all,i1,i2,i3
  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold

  !write(*,*) 'INSIDE readonewave'

  if (useFormattedInput) then
     read(unitwf,*) iorb_old,eval
  else
     read(unitwf) iorb_old,eval
  end if
  if (iorb_old.ne.iorb) stop 'readonewave'
  if (useFormattedInput) then
     read(unitwf,*) hx_old,hy_old,hz_old
     read(unitwf,*) n1_old,n2_old,n3_old
     !write(*,*) 'reading ',nat,' atomic positions'
     do iat=1,nat
        read(unitwf,*)(rxyz_old(j,iat),j=1,3)
     enddo
     read(unitwf,*) nvctr_c_old, nvctr_f_old
  else
     read(unitwf) hx_old,hy_old,hz_old
     read(unitwf) n1_old,n2_old,n3_old
     do iat=1,nat
        read(unitwf)(rxyz_old(j,iat),j=1,3)
     enddo
     read(unitwf) nvctr_c_old, nvctr_f_old
  end if

  tx=0.0_gp 
  ty=0.0_gp  
  tz=0.0_gp 
  do iat=1,nat
     tx=tx+(rxyz(1,iat)-rxyz_old(1,iat))**2
     ty=ty+(rxyz(2,iat)-rxyz_old(2,iat))**2
     tz=tz+(rxyz(3,iat)-rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)

  if (hx_old == hx .and. hy_old == hy .and. hz_old == hz .and.&
       nvctr_c_old == nvctr_c .and. nvctr_f_old == nvctr_f .and. & 
       n1_old == n1  .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3) then

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
     if (hx_old /= hx .or. hy_old /= hy .or. hz_old /= hz) write(*,*) &
          'because hgrid_old >< hgrid',hx_old,hy_old,hz_old,hx,hy,hz
     if (nvctr_c_old /= nvctr_c) write(*,*) 'because nvctr_c_old >< nvctr_c',nvctr_c_old,nvctr_c
     if (nvctr_f_old /= nvctr_f) write(*,*) 'because nvctr_f_old >< nvctr_f',nvctr_f_old,nvctr_f
     if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
          write(*,*) 'because cell size has changed',n1_old,n1,n2_old,n2,n3_old,n3
     if (displ > 1.d-3 ) write(*,*) 'large displacement of molecule'

     allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2+ndebug),stat=i_stat)
     call memocc(i_stat,psigold,'psigold',subname)

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
     call reformatonewave(iproc,displ,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,nat,&
          rxyz_old,psigold,hx,hy,hz,nvctr_c,nvctr_f,n1,n2,n3,rxyz,nseg_c,nseg_f,&
          keyg,keyv,psifscf,psi)

     i_all=-product(shape(psigold))*kind(psigold)
     deallocate(psigold,stat=i_stat)
     call memocc(i_stat,i_all,'psigold',subname)

  endif

end subroutine readonewave

subroutine writeonewave(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,nat,rxyz,  & 
     nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f, & 
     psi_c,psi_f,norb,eval)
  use module_base
  implicit none
  logical, intent(in) :: useFormattedOutput
  integer, intent(in) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,norb
  real(gp), intent(in) :: hx,hy,hz
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(norb), intent(in) :: eval
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(gp), dimension(3,nat), intent(in) :: rxyz
  !local variables
  integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7

  if (useFormattedOutput) then
     write(unitwf,*) iorb,eval(iorb)
     write(unitwf,*) hx,hy,hz
     write(unitwf,*) n1,n2,n3
     do iat=1,nat
     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
     enddo
     write(unitwf,*) nvctr_c, nvctr_f
  else
     write(unitwf) iorb,eval(iorb)
     write(unitwf) hx,hy,hz
     write(unitwf) n1,n2,n3
     do iat=1,nat
     write(unitwf) (rxyz(j,iat),j=1,3)
     enddo
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


end subroutine writeonewave
