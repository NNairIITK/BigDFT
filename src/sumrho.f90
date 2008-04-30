subroutine sumrho(geocode,iproc,nproc,norb,norbp,n1,n2,n3,hxh,hyh,hzh,occup,  & 
     wfd,psi,rho,nrho,nscatterarr,nspin,nspinor,spinsgn,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(in) :: bounds
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin,nspinor
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(kind=8), intent(in) :: hxh,hyh,hzh
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8), dimension(norb), intent(in) :: occup,spinsgn
  real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp*nspinor), intent(in) :: psi
  real(kind=8), dimension(max(nrho,1),nspin), intent(out), target :: rho
  !local variables
  include 'mpif.h'
  integer :: nw1,nw2,nrhotot,n3d,n1i,n2i,n3i,nxc,nxf
  integer :: ind1,ind2,ind3,ind1s,ind2s,ind3s,oidx,sidx,nspinn
  integer :: i00,i0,i1,i2,i3,i3off,i3s,isjmp,i,ispin,iorb,jproc,i_all,i_stat,ierr
  real(kind=8) :: hfac,hgridh,tt,charge,hfac2
  real(kind=8), dimension(0:3) :: scal
  real(kind=8), dimension(:,:), allocatable :: psir
  real(kind=8), dimension(:), allocatable :: x_c_psifscf,x_f_psig,w1,w2
  real(kind=8), dimension(:,:), pointer :: rho_p


  call timing(iproc,'Rho_comput    ','ON')

  if (iproc==0) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if

  do i=0,3
     scal(i)=1.d0
  enddo

  select case(geocode)
     case('F')
        n1i=2*n1+31
        n2i=2*n2+31
        n3i=2*n3+31
        
        !dimension of the work arrays
        ! shrink convention: nw1>nw2
        nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
             (n1+1)*(2*n2+31)*(2*n3+31),&
             2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
             2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
        nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
             4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
             (n1+1)*(n2+1)*(2*n3+31),&
             (2*n1+31)*(n2+1)*(n3+1))
        nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
        nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

     case('S')
        n1i=2*n1+2
        n2i=2*n2+31
        n3i=2*n3+2

        !dimension of the work arrays (tentative values)
        nw1=(2*n1+2)*(2*n2+31)*(2*n3+2)
        nw2=1
        nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
        nxf=(n1+1)*(n2+1)*(n3+1)*8

     case('P')
        n1i=2*n1+2
        n2i=2*n2+2
        n3i=2*n3+2
        
        !dimension of the work arrays
        nw1=(2*n1+2)*(2*n2+2)*(2*n3+2)
        nw2=1
        nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
        nxf=(n1+1)*(n2+1)*(n3+1)*8

  end select

  !work arrays
  allocate(x_c_psifscf(nxc),stat=i_stat)
  call memocc(i_stat,product(shape(x_c_psifscf))*kind(x_c_psifscf),'x_c_psifscf','sumrho')
  allocate(x_f_psig(nxf),stat=i_stat)! work
  call memocc(i_stat,product(shape(x_f_psig))*kind(x_f_psig),'x_f_psig','sumrho')
  allocate(w1(nw1),stat=i_stat)
  call memocc(i_stat,product(shape(w1))*kind(w1),'w1','sumrho')
  allocate(w2(nw2),stat=i_stat) ! work
  call memocc(i_stat,product(shape(w2))*kind(w2),'w2','sumrho')

  ! Wavefunction in real space
  nspinn=max(nspin,nspinor)
  allocate(psir(n1i*n2i*n3i,nspinn),stat=i_stat)
  call memocc(i_stat,product(shape(psir))*kind(psir),'psir','sumrho')

  !initialisation
  if (geocode == 'F') then
     call razero(nxc,x_c_psifscf)
     call razero(nxf,x_f_psig)
     
     call razero(n1i*n2i*n3i*nspinn,psir)
  end if


  !calculate dimensions of the complete array to be allocated before the reduction procedure
  nrhotot=0
  do jproc=0,nproc-1
     nrhotot=nrhotot+nscatterarr(jproc,1)
  end do

  if (nproc > 1) then
     allocate(rho_p(n1i*n2i*nrhotot,nspinn),stat=i_stat)
     call memocc(i_stat,product(shape(rho_p))*kind(rho_p),'rho_p','sumrho')
  else
     rho_p => rho
  end if

  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  if(nspinor==4) then 
     call razero(n1i*n2i*nrhotot*nspinor,rho_p)
     call tenminustwenty(n1i*n2i*nrhotot,rho_p,nproc)
  else
     call tenminustwenty(n1i*n2i*nrhotot*nspinn,rho_p,nproc)
  end if

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     hfac=(occup(iorb)/(hxh*hyh*hzh))
     hfac2=hfac*2.0d0
     
     oidx=(iorb-iproc*norbp-1)*nspinor

     if (hfac /= 0.d0) then

        select case(geocode)
           case('F')
              
              do sidx=1,nspinor
                 
                 call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                      wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
                      wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
                      scal,psi(1,oidx+sidx),psi(wfd%nvctr_c+1,oidx+sidx),x_c_psifscf,x_f_psig)
               
                 
                 call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c_psifscf,x_f_psig,  & 
                      psir(1,sidx),bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
                      bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,bounds%ibyyzz_r)
                 
              end do
              !sum different slices by taking into account the overlap
              i3s=0
              if(nspinor==1) then
                 loop_xc_overlap_F: do jproc=0,nproc-1
                    i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
                    n3d=nscatterarr(jproc,1)
                    if (n3d==0) exit loop_xc_overlap_F
                    if(spinsgn(iorb)>0.0d0) then
                       isjmp=1
                    else
                       isjmp=2
                    end if
                    do i3=i3off+1,i3off+n3d
                       i3s=i3s+1
                       ind3=(i3-1)*n1i*n2i
                       ind3s=(i3s-1)*n1i*n2i
                       do i2=1,n2i
                          ind2=(i2-1)*n1i+ind3
                          ind2s=(i2-1)*n1i+ind3s
                          !                 do i1=1,2*n1+31
                          do i1=bounds%ibyyzz_r(1,i2-15,i3-15)+1,bounds%ibyyzz_r(2,i2-15,i3-15)+1
                             ind1=i1+ind2
                             ind1s=i1+ind2s
                             !do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
                             rho_p(ind1s,isjmp)=rho_p(ind1s,isjmp)+hfac*psir(ind1,nspinor)**2
                          end do
                       end do
                    end do
                 end do loop_xc_overlap_F
              else  !similar loop for nspinor=4
                 loop_xc_overlap_FS: do jproc=0,nproc-1
                    i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
                    n3d=nscatterarr(jproc,1)
                    if (n3d==0) exit loop_xc_overlap_FS
                    do i3=i3off+1,i3off+n3d
                       i3s=i3s+1
                       ind3=(i3-1)*n1i*n2i
                       ind3s=(i3s-1)*n1i*n2i
                       do i2=1,n2i
                          ind2=(i2-1)*n1i+ind3
                          ind2s=(i2-1)*n1i+ind3s
                          !                 do i1=1,2*n1+31
                          do i1=bounds%ibyyzz_r(1,i2-15,i3-15)+1,bounds%ibyyzz_r(2,i2-15,i3-15)+1
                             ind1=i1+ind2
                             ind1s=i1+ind2s
                             !rho
                             rho_p(ind1s,1)=rho_p(ind1s,1)+hfac*psir(ind1,1)**2 &
                                  +hfac*psir(ind1,2)**2 &
                                  +hfac*psir(ind1,3)**2 &
                                  +hfac*psir(ind1,4)**2 
                             !m_x
                             rho_p(ind1s,2)= &
                             rho_p(ind1s,2)+hfac2*psir(ind1,1)*psir(ind1,3) &
                                  +hfac2*psir(ind1,2)*psir(ind1,4)
                             !m_y
                             rho_p(ind1s,3)= &
                             rho_p(ind1s,3)+hfac2*psir(ind1,1)*psir(ind1,4) &
                                  -hfac2*psir(ind1,2)*psir(ind1,3)
                             !m_z
                             rho_p(ind1s,4)=rho_p(ind1s,4)+hfac*psir(ind1,1)**2 &
                                  +hfac*psir(ind1,2)**2 &
                                  -hfac*psir(ind1,3)**2 &
                                  -hfac*psir(ind1,4)**2 
                          end do
                       end do
                    end do
                 end do loop_xc_overlap_FS
              end if

           case('P')

              do sidx=1,nspinor
                 
                 call uncompress_per(n1,n2,n3,wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),   &
                      wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),   &
                      psi(1,oidx+sidx),psi(wfd%nvctr_c+1,oidx+sidx),x_c_psifscf,x_f_psig,w1)
                 
                 call convolut_magic_n_per(2*n1+1,2*n2+1,2*n3+1,x_c_psifscf,psir(1,sidx)) 
                 
              end do
              !sum different slices by taking into account the overlap
              i3s=0
              if(nspinor==1) then
                 loop_xc_overlap_P: do jproc=0,nproc-1
                    i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
                    n3d=nscatterarr(jproc,1)
                    if (n3d==0) exit loop_xc_overlap_P
                    if(spinsgn(iorb)>0.0d0) then
                       isjmp=1
                    else
                       isjmp=2
                    end if
                    do i3=i3off+1,i3off+n3d
                       i3s=i3s+1
                       ind3=(i3-1)*n1i*n2i
                       ind3s=(i3s-1)*n1i*n2i
                       do i2=1,n2i
                          ind2=(i2-1)*n1i+ind3
                          ind2s=(i2-1)*n1i+ind3s
                          do i1=1,n1i
                             ind1=i1+ind2
                             ind1s=i1+ind2s
                             rho_p(ind1s,isjmp)=rho_p(ind1s,isjmp)+hfac*psir(ind1,nspinor)**2
                          end do
                       end do
                    end do
                 end do loop_xc_overlap_P
              else
                 loop_xc_overlap_PS: do jproc=0,nproc-1
                    i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
                    n3d=nscatterarr(jproc,1)
                    if (n3d==0) exit loop_xc_overlap_PS
                    do i3=i3off+1,i3off+n3d
                       i3s=i3s+1
                       ind3=(i3-1)*n1i*n2i
                       ind3s=(i3s-1)*n1i*n2i
                       do i2=1,n2i
                          ind2=(i2-1)*n1i+ind3
                          ind2s=(i2-1)*n1i+ind3s
                          do i1=1,n1i
                             ind1=i1+ind2
                             ind1s=i1+ind2s
                             !rho
                             rho_p(ind1s,1)=rho_p(ind1s,1)+hfac*psir(ind1,1)**2 &
                                  +hfac*psir(ind1,2)**2 &
                                  +hfac*psir(ind1,3)**2 &
                                  +hfac*psir(ind1,4)**2 
                             !m_x
                             rho_p(ind1s,2)=rho_p(ind1s,2)+hfac2*psir(ind1,1)*psir(ind1,3) &
                                  +hfac2*psir(ind1,2)*psir(ind1,4)
                             !m_y
                             rho_p(ind1s,3)=rho_p(ind1s,3)+hfac2*psir(ind1,1)*psir(ind1,4) &
                                  -hfac2*psir(ind1,2)*psir(ind1,3)
                             !m_z
                             rho_p(ind1s,4)=rho_p(ind1s,4)+hfac*psir(ind1,1)**2 &
                                  +hfac*psir(ind1,2)**2 &
                                  -hfac*psir(ind1,3)**2 &
                                  -hfac*psir(ind1,4)**2 
                          end do
                       end do
                    end do
                 end do loop_xc_overlap_PS
              end if

           end select
        if (i3s /= nrhotot) then
           print *,'problem with rhopot array in sumrho,i3s,nrhotot,',i3s,nrhotot
           stop
        end if
     end if

  enddo

  if (nproc > 1) then
     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     do ispin=1,nspin
        call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),n1i*n2i*nscatterarr(:,1),&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end do
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
  end if

  ! Check
  tt=0.d0
  i3off=n1i*n2i*nscatterarr(iproc,4)
  if(nspinor==4) then
     nspinn=1
  else
     nspinn=nspin
  end if
  do ispin=1,nspinn
     do i=1,n1i*n2i*nscatterarr(iproc,2)
        tt=tt+rho(i+i3off,ispin)
!!$        !temporary check for debugging purposes
!!$        if (rho(i+i3off) < 9.d-21) then
!!$           print *,iproc,'error in density construction',rho(i+i3off)
!!$        end if
     enddo
  end do

  if (nproc > 1) then
     i_all=-product(shape(rho_p))*kind(rho_p)
     deallocate(rho_p,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p','sumrho')

     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
  else
     !useless, only for completeness
     nullify(rho_p)

     charge=tt
  end if

  if (iproc.eq.0) write(*,'(1x,a,f21.12)')&
       'done. Total electronic charge=',charge*hxh*hyh*hzh

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir','sumrho')
  i_all=-product(shape(x_c_psifscf))*kind(x_c_psifscf)
  deallocate(x_c_psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'x_c_psifscf','sumrho')
  i_all=-product(shape(x_f_psig))*kind(x_f_psig)
  deallocate(x_f_psig,stat=i_stat)
  call memocc(i_stat,i_all,'x_f_psig','sumrho')
  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1','sumrho')
  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2','sumrho')

  call timing(iproc,'Rho_comput    ','OF')

END SUBROUTINE sumrho





