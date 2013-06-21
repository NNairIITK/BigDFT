!> @file
!!  Routines to reformat wavefunctions
!! @author
!!    Copyright (C) 2010-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 
subroutine reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,& !n(c) iproc (arg:1)
     rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1_old,n2_old,n3_old,n1,n2,n3  !n(c) iproc
  real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz_old,rxyz
  real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
  real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
  !local variables
  character(len=*), parameter :: subname='reformatonewave'
  logical :: cif1,cif2,cif3,perx,pery,perz
  integer :: i_stat,i_all,i1,i2,i3,j1,j2,j3,l1,l2,iat,nb1,nb2,nb3,ind,jj1,jj2,jj3a,jj3b,jj3c
  real(gp) :: hxh,hyh,hzh,hxh_old,hyh_old,hzh_old,x,y,z,dx,dy,dz,xold,yold,zold,mindist
  real(wp) :: zr,yr,xr,ym1,y00,yp1
  real(wp), dimension(-1:1,-1:1) :: xya
  real(wp), dimension(-1:1) :: xa
  real(wp), dimension(:), allocatable :: ww,wwold
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psig
  real(wp), dimension(:,:,:), allocatable :: psifscfold

  !conditions for periodicity in the three directions
  perx=(at%astruct%geocode /= 'F')
  pery=(at%astruct%geocode == 'P')
  perz=(at%astruct%geocode /= 'F')

  !buffers related to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(perx,nb1)
  call ext_buffers_coarse(pery,nb2)
  call ext_buffers_coarse(perz,nb3)

  allocate(psifscfold(-nb1:2*n1_old+1+nb1,-nb2:2*n2_old+1+nb2,-nb3:2*n3_old+1+nb3+ndebug),stat=i_stat)
  call memocc(i_stat,psifscfold,'psifscfold',subname)
  allocate(wwold((2*n1_old+2+2*nb1)*(2*n2_old+2+2*nb2)*(2*n3_old+2+2*nb3)+ndebug),stat=i_stat)
  call memocc(i_stat,wwold,'wwold',subname)

  if (at%astruct%geocode=='F') then
     call synthese_grow(n1_old,n2_old,n3_old,wwold,psigold,psifscfold) 
  else if (at%astruct%geocode=='S') then     
     call synthese_slab(n1_old,n2_old,n3_old,wwold,psigold,psifscfold) 
  else if (at%astruct%geocode=='P') then     
     call synthese_per(n1_old,n2_old,n3_old,wwold,psigold,psifscfold) 
  end if

  i_all=-product(shape(wwold))*kind(wwold)
  deallocate(wwold,stat=i_stat)
  call memocc(i_stat,i_all,'wwold',subname)
  
  !write(*,*) iproc,' displ ',displ
  if (hx == hx_old .and. hy == hy_old .and. hz == hz_old .and. &
       n1_old==n1 .and. n2_old==n2 .and. n3_old==n3 .and. &
       displ<= 1.d-2) then
     !if (iproc==0) write(*,*) iproc,' orbital just copied'
     call dcopy((2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(2*n3+2+2*nb3),psifscfold(-nb1,-nb2,-nb3),1,&
          psifscf(1),1)
!!$     do i3=-nb3,2*n3+1+nb3
!!$        do i2=-nb2,2*n2+1+nb2
!!$           do i1=-nb1,2*n1+1+nb1
!!$              ind=i1+nb1+1+(2*n1+2+2*nb1)*(i2+nb2)+&
!!$                   (2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(i3+nb3)
!!$              psifscf(ind)=psifscfold(i1,i2,i3)
!!$           enddo
!!$        enddo
!!$     enddo
     
  else

     dx=0.0_gp
     dy=0.0_gp 
     dz=0.0_gp
     !Calculate average shift
     !Take into account the modulo operation which should be done for non-isolated BC
     do iat=1,at%astruct%nat 
        dx=dx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))
        dy=dy+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))
        dz=dz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))
     enddo
     dx=dx/real(at%astruct%nat,gp)
     dy=dy/real(at%astruct%nat,gp)
     dz=dz/real(at%astruct%nat,gp)
     
     ! transform to new structure    
     !if (iproc==0) write(*,*) iproc,' orbital fully transformed'
     hxh=.5_gp*hx
     hxh_old=.5_gp*hx_old
     hyh=.5_gp*hy
     hyh_old=.5_gp*hy_old
     hzh=.5_gp*hz
     hzh_old=.5_gp*hz_old

     call razero((2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(2*n3+2+2*nb3),psifscf)

     do i3=-nb3,2*n3+1+nb3
        z=real(i3,gp)*hzh
        do i2=-nb2,2*n2+1+nb2
           y=real(i2,gp)*hyh
           do i1=-nb1,2*n1+1+nb1
              x=real(i1,gp)*hxh

              xold=x-dx 
              yold=y-dy
              zold=z-dz

              j1=nint((xold)/hxh_old)
              cif1=(j1 >= -6 .and. j1 <= 2*n1_old+7) .or. perx
              j2=nint((yold)/hyh_old)
              cif2=(j2 >= -6 .and. j2 <= 2*n2_old+7) .or. pery
              j3=nint((zold)/hzh_old)
              cif3=(j3 >= -6 .and. j3 <= 2*n3_old+7) .or. perz

              ind=i1+nb1+1+(2*n1+2+2*nb1)*(i2+nb2)+&
                   (2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(i3+nb3)

              
              if (cif1 .and. cif2 .and. cif3) then 
                 zr =real(((z-dz)-real(j3,gp)*hzh_old)/hzh_old,wp)
                 do l2=-1,1
                    do l1=-1,1
                       !the modulo has no effect on free BC thanks to the
                       !if statement above
                       jj1=modulo(j1+l1+nb1,2*n1_old+1+2*nb1+1)-nb1
                       jj2=modulo(j2+l2+nb2,2*n2_old+1+2*nb2+1)-nb2
                       jj3a=modulo(j3-1+nb3,2*n3_old+1+2*nb3+1)-nb3
                       jj3b=modulo(j3  +nb3,2*n3_old+1+2*nb3+1)-nb3
                       jj3c=modulo(j3+1+nb3,2*n3_old+1+2*nb3+1)-nb3
                       

                       ym1=psifscfold(jj1,jj2,jj3a)
                       y00=psifscfold(jj1,jj2,jj3b)
                       yp1=psifscfold(jj1,jj2,jj3c)

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
                 psifscf(ind)=ym1 + &
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
  allocate(ww((2*n1+2+2*nb1)*(2*n2+2+2*nb2)*(2*n3+2+2*nb3)+ndebug),stat=i_stat)
  call memocc(i_stat,ww,'ww',subname)

  if (at%astruct%geocode=='F') then
     call analyse_shrink(n1,n2,n3,ww,psifscf,psig)
  else if (at%astruct%geocode == 'S') then
     call analyse_slab(n1,n2,n3,ww,psifscf,psig)
  else if (at%astruct%geocode == 'P') then
     call analyse_per(n1,n2,n3,ww,psifscf,psig)
  end if

  !write(100+iproc,*) 'norm new psig ',dnrm2(8*(n1+1)*(n2+1)*(n3+1),psig,1)
  call compress_plain(n1,n2,0,n1,0,n2,0,n3,  &
       wfd%nseg_c,wfd%nvctr_c,wfd%keygloc(1,1),wfd%keyvloc(1),   &
       wfd%nseg_f,wfd%nvctr_f,&
       wfd%keygloc(1,wfd%nseg_c+min(1,wfd%nseg_f)),&
       wfd%keyvloc(wfd%nseg_c+min(1,wfd%nseg_f)),   &
       psig,psi(1),psi(wfd%nvctr_c+min(1,wfd%nvctr_f)))

  !write(100+iproc,*) 'norm of reformatted psi ',dnrm2(nvctr_c+7*nvctr_f,psi,1)

  i_all=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=i_stat)
  call memocc(i_stat,i_all,'psig',subname)
  i_all=-product(shape(ww))*kind(ww)
  deallocate(ww,stat=i_stat)
  call memocc(i_stat,i_all,'ww',subname)

END SUBROUTINE reformatonewave


!> Calculates the minimum difference between two coordinates
!! knowing that there could have been a modulo operation
function mindist(periodic,alat,r,r_old)
  use module_base
  implicit none
  logical, intent(in) :: periodic
  real(gp), intent(in) :: r,r_old,alat
  real(gp) :: mindist

  !for periodic BC calculate mindist only if the center of mass can be defined without the modulo
  if (periodic) then
     if (r_old > 0.5_gp*alat) then
        if (r < 0.5_gp*alat) then
           !mindist=r+alat-r_old
           mindist=0.0_gp
        else
           mindist=r-r_old
        end if
     else
        if (r > 0.5_gp*alat) then
           !mindist=r-alat-r_old
           mindist=0.0_gp
        else
           mindist=r-r_old
        end if
     end if
  else
     mindist=r-r_old
  end if

end function mindist


subroutine ext_buffers_coarse(periodic,nb)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(out) :: nb
  if (periodic) then
     nb=0
  else
     nb=7
  end if
END SUBROUTINE ext_buffers_coarse


module internal_io
  implicit none

contains
  subroutine io_error(error)
    use module_defs

    implicit none

    character(len = *), intent(in) :: error
    integer :: ierr

    call io_warning(error)
    call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
  END SUBROUTINE io_error

  subroutine io_warning(error)
    use module_defs

    implicit none

    character(len = *), intent(in) :: error

    write(0,"(2A)") "WARNING! ", trim(error)
  END SUBROUTINE io_warning

  subroutine io_read_descr(unitwf, formatted, iorb_old, eval, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
    use module_base
    use module_types

    implicit none

    integer, intent(in) :: unitwf
    logical, intent(in) :: formatted
    integer, intent(out) :: iorb_old
    integer, intent(out) :: n1_old, n2_old, n3_old
    real(gp), intent(out) :: hx_old, hy_old, hz_old
    logical, intent(out) :: lstat
    real(wp), intent(out) :: eval
    character(len =256), intent(out) :: error
    ! Optional arguments
    integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
    integer, intent(in), optional :: nat
    real(gp), dimension(:,:), intent(out), optional :: rxyz_old

    integer :: i, iat, i_stat, nat_
    real(gp) :: rxyz(3)

    lstat = .false.
    write(error, "(A)") "cannot read psi description."
    if (formatted) then
       read(unitwf,*,iostat=i_stat) iorb_old,eval
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
       if (i_stat /= 0) return
       !write(*,*) 'reading ',nat,' atomic positions'
       if (present(nat) .And. present(rxyz_old)) then
          read(unitwf,*,iostat=i_stat) nat_
          if (i_stat /= 0) return
          ! Sanity check
          if (size(rxyz_old, 2) /= nat) call io_error("Mismatch in coordinate array size.")
          if (nat_ /= nat) call io_error("Mismatch in coordinate array size.")
          do iat=1,nat
             read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
             if (i_stat /= 0) return
          enddo
       else
          read(unitwf,*,iostat=i_stat) nat_
          if (i_stat /= 0) return
          do iat=1,nat_
             read(unitwf,*,iostat=i_stat)
             if (i_stat /= 0) return
          enddo
       end if
       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
          read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
          if (i_stat /= 0) return
       else
          read(unitwf,*,iostat=i_stat) i, iat
          if (i_stat /= 0) return
       end if
    else
       read(unitwf,iostat=i_stat) iorb_old,eval
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
       if (i_stat /= 0) return
       if (present(nat) .And. present(rxyz_old)) then
          read(unitwf,iostat=i_stat) nat_
          if (i_stat /= 0) return
          ! Sanity check
          if (size(rxyz_old, 2) /= nat) call io_error("Mismatch in coordinate array size.")
          if (nat_ /= nat) call io_error("Mismatch in coordinate array size.")
          do iat=1,nat
             read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
             if (i_stat /= 0) return
          enddo
       else
          read(unitwf,iostat=i_stat) nat_
          if (i_stat /= 0) return
          do iat=1,nat_
             read(unitwf,iostat=i_stat) rxyz
             if (i_stat /= 0) return
          enddo
       end if
       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
          read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
          if (i_stat /= 0) return
       else
          read(unitwf,iostat=i_stat) i, iat
          if (i_stat /= 0) return
       end if
    end if
    lstat = .true.
  END SUBROUTINE io_read_descr

  subroutine io_gcoordToLocreg(n1, n2, n3, nvctr_c, nvctr_f, gcoord_c, gcoord_f, lr)
    use module_defs
    use module_types

    implicit none

    integer, intent(in) :: n1, n2, n3, nvctr_c, nvctr_f
    integer, dimension(3, nvctr_c), intent(in) :: gcoord_c
    integer, dimension(3, nvctr_f), intent(in) :: gcoord_f
    type(locreg_descriptors), intent(out) :: lr

    character(len = *), parameter :: subname = "io_gcoordToLocreg"
    integer :: i, i_stat, i_all
    logical, dimension(:,:,:), allocatable :: logrid_c, logrid_f

    lr%geocode = "P"
    lr%hybrid_on = .false.

    lr%ns1 = 0
    lr%ns2 = 0
    lr%ns3 = 0

    lr%d%n1 = n1
    lr%d%n2 = n2
    lr%d%n3 = n3

    lr%d%n1i = 2 * n1 + 2
    lr%d%n2i = 2 * n2 + 2
    lr%d%n3i = 2 * n3 + 2

    allocate(logrid_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,logrid_c,'logrid_c',subname)
    allocate(logrid_f(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,logrid_f,'logrid_f',subname)

    lr%d%nfl1 = n1
    lr%d%nfl2 = n2
    lr%d%nfl3 = n3
    lr%d%nfu1 = 0
    lr%d%nfu2 = 0
    lr%d%nfu3 = 0

    logrid_c(:,:,:) = .false.
    do i = 1, nvctr_c, 1
       logrid_c(gcoord_c(1, i), gcoord_c(2, i), gcoord_c(3, i)) = .true.
    end do
    logrid_f(:,:,:) = .false.
    do i = 1, nvctr_f, 1
       logrid_f(gcoord_f(1, i), gcoord_f(2, i), gcoord_f(3, i)) = .true.
       lr%d%nfl1 = min(lr%d%nfl1, gcoord_f(1, i))
       lr%d%nfl2 = min(lr%d%nfl2, gcoord_f(2, i))
       lr%d%nfl3 = min(lr%d%nfl3, gcoord_f(3, i))
       lr%d%nfu1 = max(lr%d%nfu1, gcoord_f(1, i))
       lr%d%nfu2 = max(lr%d%nfu2, gcoord_f(2, i))
       lr%d%nfu3 = max(lr%d%nfu3, gcoord_f(3, i))
    end do

    !correct the values of the delimiter if there are no wavelets
    if (lr%d%nfl1 == n1 .and. lr%d%nfu1 == 0) then
       lr%d%nfl1 = n1 / 2
       lr%d%nfu1 = n1 / 2
    end if
    if (lr%d%nfl2 == n2 .and. lr%d%nfu2 == 0) then
       lr%d%nfl2 = n2 / 2
       lr%d%nfu2 = n2 / 2
    end if
    if (lr%d%nfl3 == n3 .and. lr%d%nfu3 == 0) then
       lr%d%nfl3 = n3 / 2
       lr%d%nfu3 = n3 / 2
    end if

    call wfd_from_grids(logrid_c, logrid_f, lr)

    i_all=-product(shape(logrid_c))*kind(logrid_c)
    deallocate(logrid_c,stat=i_stat)
    call memocc(i_stat,i_all,'logrid_c',subname)
    i_all=-product(shape(logrid_f))*kind(logrid_f)
    deallocate(logrid_f,stat=i_stat)
    call memocc(i_stat,i_all,'logrid_f',subname)
  END SUBROUTINE io_gcoordToLocreg

  subroutine read_psi_compress(unitwf, formatted, nvctr_c, nvctr_f, psi, lstat, error, gcoord_c, gcoord_f)
    use module_base
    use module_types

    implicit none

    integer, intent(in) :: unitwf, nvctr_c, nvctr_f
    logical, intent(in) :: formatted
    real(wp), dimension(nvctr_c+7*nvctr_f), intent(out) :: psi
    logical, intent(out) :: lstat
    character(len =256), intent(out) :: error
    integer, dimension(3, nvctr_c), optional, intent(out) :: gcoord_c
    integer, dimension(3, nvctr_f), optional, intent(out) :: gcoord_f

    integer :: j, i1, i2, i3, i_stat
    real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7

    lstat = .false.
    write(error, "(A)") "cannot read psi values."
    if (present(gcoord_c)) then
       do j=1,nvctr_c
          if (formatted) then
             read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
          else
             read(unitwf,iostat=i_stat) i1,i2,i3,tt
          end if
          if (i_stat /= 0) return
          psi(j)=tt
          gcoord_c(:, j) = (/ i1, i2, i3 /)
       enddo
    else
       do j=1,nvctr_c
          if (formatted) then
             read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
          else
             read(unitwf,iostat=i_stat) i1,i2,i3,tt
          end if
          if (i_stat /= 0) return
          psi(j)=tt
       enddo
    end if
    if (present(gcoord_f)) then
       do j=1,7*nvctr_f-6,7
          if (formatted) then
             read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
          else
             read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
          end if
          if (i_stat /= 0) return
          psi(nvctr_c+j+0)=t1
          psi(nvctr_c+j+1)=t2
          psi(nvctr_c+j+2)=t3
          psi(nvctr_c+j+3)=t4
          psi(nvctr_c+j+4)=t5
          psi(nvctr_c+j+5)=t6
          psi(nvctr_c+j+6)=t7
          gcoord_f(:, (j - 1) / 7 + 1) = (/ i1, i2, i3 /)
       enddo
    else
       do j=1,7*nvctr_f-6,7
          if (formatted) then
             read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
          else
             read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
          end if
          if (i_stat /= 0) return
          psi(nvctr_c+j+0)=t1
          psi(nvctr_c+j+1)=t2
          psi(nvctr_c+j+2)=t3
          psi(nvctr_c+j+3)=t4
          psi(nvctr_c+j+4)=t5
          psi(nvctr_c+j+5)=t6
          psi(nvctr_c+j+6)=t7
       enddo
    end if
    lstat = .true.
  END SUBROUTINE read_psi_compress

  subroutine read_psig(unitwf, formatted, nvctr_c, nvctr_f, n1, n2, n3, psig, lstat, error)
    use module_base
    use module_types

    implicit none

    integer, intent(in) :: unitwf, nvctr_c, nvctr_f, n1, n2, n3
    logical, intent(in) :: formatted
    real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psig
    logical, intent(out) :: lstat
    character(len =256), intent(out) :: error

    integer :: i1, i2, i3, i_stat, iel
    real(wp) :: tt, t1, t2, t3, t4, t5, t6, t7

    lstat = .false.
    write(error, "(A)") "cannot read psig values."

    call razero(8*(n1+1)*(n2+1)*(n3+1),psig)
    do iel=1,nvctr_c
       if (formatted) then
          read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
       else
          read(unitwf,iostat=i_stat) i1,i2,i3,tt
       end if
       if (i_stat /= 0) return
       psig(i1,1,i2,1,i3,1)=tt
    enddo
    do iel=1,nvctr_f
       if (formatted) then
          read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
       else
          read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
       end if
       if (i_stat /= 0) return
       psig(i1,2,i2,1,i3,1)=t1
       psig(i1,1,i2,2,i3,1)=t2
       psig(i1,2,i2,2,i3,1)=t3
       psig(i1,1,i2,1,i3,2)=t4
       psig(i1,2,i2,1,i3,2)=t5
       psig(i1,1,i2,2,i3,2)=t6
       psig(i1,2,i2,2,i3,2)=t7
    enddo
    lstat = .true.
  END SUBROUTINE read_psig

  subroutine io_open(unitwf, filename, formatted)
    character(len = *), intent(in) :: filename
    logical, intent(in) :: formatted
    integer, intent(out) :: unitwf

    integer :: i_stat

    ! We open the Fortran file
    unitwf = 99
    if (.not. formatted) then
       open(unit=unitwf,file=trim(filename),status='unknown',form="unformatted", iostat=i_stat)
    else
       open(unit=unitwf,file=trim(filename),status='unknown', iostat=i_stat)
    end if
    if (i_stat /= 0) then
       call io_warning("Cannot open file '" // trim(filename) // "'.")
       unitwf = -1
       return
    end if
  END SUBROUTINE io_open

END MODULE internal_io

subroutine readonewave(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
     & hx,hy,hz,at,wfd,rxyz_old,rxyz,psi,eval,psifscf)
  use module_base
  use module_types
  use internal_io
  use yaml_output
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), intent(out) :: eval
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
  real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
  !local variables
  character(len=*), parameter :: subname='readonewave'
  character(len = 256) :: error
  logical :: perx,pery,perz,lstat
  integer :: iorb_old,n1_old,n2_old,n3_old,iat,iel,nvctr_c_old,nvctr_f_old,i_stat,i_all,i1,i2,i3
  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old,mindist
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold

  !write(*,*) 'INSIDE readonewave'
  call io_read_descr(unitwf, useFormattedInput, iorb_old, eval, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, at%astruct%nat)
  if (.not. lstat) call io_error(trim(error))
  if (iorb_old /= iorb) stop 'readonewave'

  !conditions for periodicity in the three directions
  perx=(at%astruct%geocode /= 'F')
  pery=(at%astruct%geocode == 'P')
  perz=(at%astruct%geocode /= 'F')

  tx=0.0_gp 
  ty=0.0_gp
  tz=0.0_gp
  do iat=1,at%astruct%nat
     tx=tx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
     ty=ty+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
     tz=tz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)

  if (hx_old == hx .and. hy_old == hy .and. hz_old == hz .and.&
       nvctr_c_old == wfd%nvctr_c .and. nvctr_f_old == wfd%nvctr_f .and. & 
       n1_old == n1  .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3) then

     !if (iproc == 0) write(*,*) 'wavefunction ',iorb,' needs NO reformatting on processor',iproc
     if (iproc == 0) call yaml_map('Need to reformat wavefunctions',.false.)
     !if (iproc == 0) write(*,*) 'wavefunctions need NO reformatting'
     call read_psi_compress(unitwf, useFormattedInput, wfd%nvctr_c, wfd%nvctr_f, psi, lstat, error)
     if (.not. lstat) call io_error(trim(error))

  else

     if (iproc == 0 .and. iorb == 1) then
        call yaml_map('Need to reformat wavefunctions',.false.)
        !write(*,*) 'wavefunctions need reformatting'
        if (hx_old /= hx .or. hy_old /= hy .or. hz_old /= hz) &
           & call yaml_comment('because hgrid_old /= hgrid' // &
             & trim(yaml_toa((/ hx_old,hy_old,hz_old,hx,hy,hz /), fmt='(f14.10)')))
           ! & write(*,"(1x,A,6F14.10)") 'because hgrid_old /= hgrid',hx_old,hy_old,hz_old,hx,hy,hz
        if (nvctr_c_old /= wfd%nvctr_c) &
           & call yaml_comment('because nvctr_c_old /= nvctr_c' // trim(yaml_toa((/ nvctr_c_old,wfd%nvctr_c /))))
           ! & write(*,*) 'because nvctr_c_old /= nvctr_c', nvctr_c_old,wfd%nvctr_c
        if (nvctr_f_old /= wfd%nvctr_f) &
           & call yaml_comment('because nvctr_f_old /= nvctr_f' // trim(yaml_toa((/ nvctr_f_old,wfd%nvctr_f /))))
           ! & write(*,*) 'because nvctr_f_old /= nvctr_f', nvctr_f_old,wfd%nvctr_f
        if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
           call yaml_comment('because cell size has changed' // trim(yaml_toa((/ n1_old,n1,n2_old,n2,n3_old,n3 /))))
           ! & write(*,*) 'because cell size has changed',n1_old,n1,n2_old,n2,n3_old,n3
        if (displ > 1.d-3 ) call yaml_comment('large displacement of molecule' // trim(yaml_toa(displ)))
        !if (displ > 1.d-3 ) write(*,*) 'large displacement of molecule',displ
     end if

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
     call reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,&
          rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)

     i_all=-product(shape(psigold))*kind(psigold)
     deallocate(psigold,stat=i_stat)
     call memocc(i_stat,i_all,'psigold',subname)

  endif

END SUBROUTINE readonewave

subroutine readwavetoisf(lstat, filename, formatted, hx, hy, hz, &
     & n1, n2, n3, nspinor, psiscf)
  use module_base
  use module_types
  use internal_io

  implicit none

  character(len = *), intent(in) :: filename
  logical, intent(in) :: formatted
  integer, intent(out) :: n1, n2, n3, nspinor
  real(gp), intent(out) :: hx, hy, hz
  real(wp), dimension(:,:,:,:), pointer :: psiscf
  logical, intent(out) :: lstat

  character(len = *), parameter :: subname = "readwavetoisf"
  integer :: unitwf, iorb, i_all, i_stat, ispinor, ispin, ikpt
  integer, dimension(:,:), allocatable :: gcoord_c, gcoord_f
  real(wp) :: eval
  real(wp), dimension(:), allocatable :: psi
  type(locreg_descriptors) :: lr
  character(len = 256) :: error
  type(workarr_sumrho) :: w
  character(len = 1024) :: fileRI
  integer :: n1_old, n2_old, n3_old, nvctr_c_old, nvctr_f_old
  real(gp) :: hx_old, hy_old, hz_old

  ! We open the Fortran file
  call io_open(unitwf, filename, formatted)
  if (unitwf < 0) then
     return
  end if

  ! We read the basis set description and the atomic definition.
  call io_read_descr(unitwf, formatted, iorb, eval, n1, n2, n3, &
       & hx, hy, hz, lstat, error, lr%wfd%nvctr_c, lr%wfd%nvctr_f)
  if (.not. lstat) then
     call io_warning(trim(error))
     return
  end if
  ! Do a magic here with the filenames...
  call readwavedescr(lstat, filename, iorb, ispin, ikpt, ispinor, nspinor, fileRI)
  if (.not. lstat) then
     call io_warning("cannot read wave ids from filename.")
     return
  end if

  ! Initial allocations.
  allocate(gcoord_c(3,lr%wfd%nvctr_c + ndebug),stat=i_stat)
  call memocc(i_stat,gcoord_c,'gcoord_c',subname)
  allocate(gcoord_f(3,lr%wfd%nvctr_f + ndebug),stat=i_stat)
  call memocc(i_stat,gcoord_f,'gcoord_f',subname)
  allocate(psi(lr%wfd%nvctr_c + 7 * lr%wfd%nvctr_f + ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  ! Read psi and the basis-set
  call read_psi_compress(unitwf, formatted, lr%wfd%nvctr_c, lr%wfd%nvctr_f, psi, lstat, error, gcoord_c, gcoord_f)
  if (.not. lstat) then
     call io_warning(trim(error))
     call deallocate_local()
     return
  end if
  call io_gcoordToLocreg(n1, n2, n3, lr%wfd%nvctr_c, lr%wfd%nvctr_f, &
       & gcoord_c, gcoord_f, lr)

  allocate(psiscf(lr%d%n1i, lr%d%n2i, lr%d%n3i, nspinor + ndebug),stat=i_stat)
  call memocc(i_stat,psiscf,'psiscf',subname)

  call initialize_work_arrays_sumrho(lr,w)

  ! Magic-filter to isf
  call daub_to_isf(lr, w, psi, psiscf(1,1,1,ispinor))

  ! Read the other psi part, if any
  if (nspinor > 1) then
     close(unitwf)
     n1_old = n1
     n2_old = n2
     n3_old = n3
     hx_old = hx
     hy_old = hy
     hz_old = hz
     nvctr_c_old = lr%wfd%nvctr_c
     nvctr_f_old = lr%wfd%nvctr_f
     
     ispinor = modulo(ispinor, 2) + 1
     call io_open(unitwf, trim(fileRI), formatted)
     if (unitwf < 0) then
        call io_warning("cannot read other spinor part from '" // trim(fileRI) // "'.")
        call deallocate_local()
        return
     end if
     ! We read the basis set description and the atomic definition.
     call io_read_descr(unitwf, formatted, iorb, eval, n1, n2, n3, &
          & hx, hy, hz, lstat, error, lr%wfd%nvctr_c, lr%wfd%nvctr_f)
     if (.not. lstat) then
        call io_warning(trim(error))
        call deallocate_local()
        return
     end if
     ! Check consistency of the basis-set.
     if (n1_old == n1 .and. n2_old == n2 .and. n3_old == n3 .and. &
          & hx_old == hx .and. hy_old == hy .and. hz_old == hz .and. &
          & nvctr_c_old == lr%wfd%nvctr_c .and. nvctr_f_old == lr%wfd%nvctr_f) then
        call read_psi_compress(unitwf, formatted, lr%wfd%nvctr_c, lr%wfd%nvctr_f, psi, lstat, error)
        if (.not. lstat) then
           call io_warning(trim(error))
           call deallocate_local()
           return
        end if
        call daub_to_isf(lr, w, psi, psiscf(1,1,1,ispinor))
     else
        call io_warning("It exists a file with the same naming convention" // &
             & " but with a different basis-set.")
        hx = hx_old
        hy = hy_old
        hz = hz_old
        psiscf(:,:,:,ispinor) = real(0, wp)
     end if
  end if

  ! We update the size values to match the allocation of psiscf.
  n1 = lr%d%n1i
  n2 = lr%d%n2i
  n3 = lr%d%n3i
  hx = hx * 0.5d0
  hy = hy * 0.5d0
  hz = hz * 0.5d0

  call deallocate_local()
  lstat = .true.

contains

  subroutine deallocate_local()
    character(len = *), parameter :: subname = "readwavetoisf"

    ! We close the file.
    close(unit=unitwf)

    if (allocated(psi)) then
       i_all=-product(shape(psi))*kind(psi)
       deallocate(psi,stat=i_stat)
       call memocc(i_stat,i_all,'psi',subname)
    end if

    if (allocated(gcoord_c)) then
       i_all=-product(shape(gcoord_c))*kind(gcoord_c)
       deallocate(gcoord_c,stat=i_stat)
       call memocc(i_stat,i_all,'gcoord_c',subname)
    end if
    if (allocated(gcoord_f)) then
       i_all=-product(shape(gcoord_f))*kind(gcoord_f)
       deallocate(gcoord_f,stat=i_stat)
       call memocc(i_stat,i_all,'gcoord_f',subname)
    end if

    if (associated(w%x_c)) then
       call deallocate_work_arrays_sumrho(w)
    end if
    if (associated(lr%bounds%kb%ibyz_f)) then
       call deallocate_bounds(lr%geocode, lr%hybrid_on, lr%bounds, subname)
    end if
    call deallocate_wfd(lr%wfd, subname)
  END SUBROUTINE deallocate_local
END SUBROUTINE readwavetoisf

subroutine readwavedescr(lstat, filename, iorb, ispin, ikpt, ispinor, nspinor, fileRI)
  use module_base
  use module_types

  implicit none

  character(len = *), intent(in) :: filename
  integer, intent(out) :: iorb, ispin, ikpt, nspinor, ispinor
  logical, intent(out) :: lstat
  character(len = 1024), intent(out) :: fileRI

  logical :: exists
  integer :: i, i_stat
  character(len = 1) :: code

  lstat = .false.

  !find the value of iorb
  read(filename(index(filename, ".", back = .true.)+2:len(filename)),*,iostat = i_stat) iorb
  if (i_stat /= 0) return
  i = index(filename, "-k", back = .true.)+2
  read(filename(i:i+2),*,iostat = i_stat) ikpt
  if (i_stat /= 0) return
  i = index(filename, "-", back = .true.)+1
  read(filename(i:i),*,iostat = i_stat) code
  if (i_stat /= 0) return
  if (code == "U" .or. code == "N") ispin = 1
  if (code == "D") ispin = 2
  ! Test file for other spinor part.
  nspinor = 1
  ispinor = 1
  read(filename(i+1:i+1),*,iostat = i_stat) code
  if (i_stat /= 0) return
  if (code == "R") then
     write(fileRI, "(A)") filename
     fileRI(i+1:i+1) = "I"
     inquire(file=trim(fileRI), exist=exists)
     if (exists) then
        ispinor = 1
        nspinor = 2
     end if
  end if
  if (code == "I") then
     write(fileRI, "(A)") filename
     fileRI(i+1:i+1) = "R"
     inquire(file=trim(fileRI), exist=exists)
     if (exists) then
        ispinor = 2
        nspinor = 2
     end if
  end if

  lstat = .true.
END SUBROUTINE readwavedescr


subroutine writeonewave(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,nat,rxyz,  & 
     nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f, & 
     psi_c,psi_f,eval)
  use module_base
  use yaml_output
  implicit none
  logical, intent(in) :: useFormattedOutput
  integer, intent(inout) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz
  real(wp), intent(in) :: eval
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(gp), dimension(3,nat), intent(in) :: rxyz
  !local variables
  integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7

  if (useFormattedOutput) then
     write(unitwf,*) iorb,eval
     write(unitwf,*) hx,hy,hz
     write(unitwf,*) n1,n2,n3
     write(unitwf,*) nat
     do iat=1,nat
     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
     enddo
     write(unitwf,*) nvctr_c, nvctr_f
  else
     write(unitwf) iorb,eval
     write(unitwf) hx,hy,hz
     write(unitwf) n1,n2,n3
     write(unitwf) nat
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

  if (verbose >= 2) call yaml_comment(trim(yaml_toa(iorb)) //'th wavefunction written')
  !if (verbose >= 2) write(*,'(1x,i0,a)') iorb,'th wavefunction written'

END SUBROUTINE writeonewave



!> @file
!!  Routine to reformat one support function (linear scaling). Adapted version of reformatonewave.
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 
! make frag_trans the argument so can eliminate need for interface
subroutine reformat_one_supportfunction(wfd,geocode,hgrids_old,n_old,psigold,& 
     hgrids,n,centre_old,centre_new,da,frag_trans,psi)
  use module_base
  use module_types
  use module_fragments
  implicit none
  integer, dimension(3), intent(in) :: n,n_old
  real(gp), dimension(3), intent(in) :: hgrids,hgrids_old
  type(wavefunctions_descriptors), intent(in) :: wfd
  character(len=1), intent(in) :: geocode
  real(gp), dimension(3), intent(inout) :: centre_old,centre_new,da
  type(fragment_transformation), intent(in) :: frag_trans
  real(wp), dimension(0:n_old(1),2,0:n_old(2),2,0:n_old(3),2), intent(in) :: psigold
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi

  !local variables
  character(len=*), parameter :: subname='reformatonesupportfunction'
  logical, dimension(3) :: per
  integer :: i_stat,i_all
  integer, dimension(3) :: nb, ndims_tmp
  real(gp), dimension(3) :: hgridsh,hgridsh_old
  !!real(wp) :: dnrm2
  real(wp), dimension(:), allocatable :: ww,wwold
  real(wp), dimension(:), allocatable :: x_phi, y_phi
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psig
  real(wp), dimension(:,:,:), pointer :: psifscfold, psifscf, psifscf_tmp
  integer :: itype, nd, nrange

  !if (size(frag_trans%discrete_operations)>0) then
  !   print*,'Error discrete operations not yet implemented',size(frag_trans%discrete_operations),&
  !        frag_trans%discrete_operations(1)
  !   stop
  !end if

  !conditions for periodicity in the three directions
  per(1)=(geocode /= 'F')
  per(2)=(geocode == 'P')
  per(3)=(geocode /= 'F')

  !buffers related to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(per(1),nb(1))
  call ext_buffers_coarse(per(2),nb(2))
  call ext_buffers_coarse(per(3),nb(3))

  allocate(psifscf(-nb(1):2*n(1)+1+nb(1),-nb(2):2*n(2)+1+nb(2),-nb(3):2*n(3)+1+nb(3)+ndebug),stat=i_stat)
  call memocc(i_stat,psifscf,'psifscf',subname)
  allocate(psifscfold(-nb(1):2*n_old(1)+1+nb(1),-nb(2):2*n_old(2)+1+nb(2),-nb(3):2*n_old(3)+1+nb(3)+ndebug),stat=i_stat)
  call memocc(i_stat,psifscfold,'psifscfold',subname)
  allocate(wwold((2*n_old(1)+2+2*nb(1))*(2*n_old(2)+2+2*nb(2))*(2*n_old(3)+2+2*nb(3))+ndebug),stat=i_stat)
  call memocc(i_stat,wwold,'wwold',subname)

  if (geocode=='F') then
     call synthese_grow(n_old(1),n_old(2),n_old(3),wwold,psigold,psifscfold) 
  else if (geocode=='S') then     
     call synthese_slab(n_old(1),n_old(2),n_old(3),wwold,psigold,psifscfold) 
  else if (geocode=='P') then     
     call synthese_per(n_old(1),n_old(2),n_old(3),wwold,psigold,psifscfold) 
  end if

  i_all=-product(shape(wwold))*kind(wwold)
  deallocate(wwold,stat=i_stat)
  call memocc(i_stat,i_all,'wwold',subname)

  ! transform to new structure    
  hgridsh=.5_gp*hgrids
  hgridsh_old=.5_gp*hgrids_old

  !create the scaling function array
  !use lots of points (to optimize one can determine how many points are needed at max)
  itype=16
  nd=2**20

  allocate(x_phi(0:nd+ndebug),stat=i_stat )
  call memocc(i_stat,x_phi,'x_phi',subname)
  allocate(y_phi(0:nd+ndebug) ,stat=i_stat )
  call memocc(i_stat,y_phi,'y_phi',subname)

  call my_scaling_function4b2B(itype,nd,nrange,x_phi,y_phi) 
  if( abs(y_phi(nd/2)-1)>1.0e-10 ) then
     stop " wrong scaling function 4b2B: not a centered one "
  endif

  i_all=-product(shape(x_phi))*kind(x_phi)
  deallocate(x_phi,stat=i_stat)
  call memocc(i_stat,i_all,'x_phi',subname)

  nullify(psifscf_tmp)
  if (size(frag_trans%discrete_operations)==0) then
     psifscf_tmp => psifscfold
     ndims_tmp=(2*n_old+2+2*nb)
  else if (size(frag_trans%discrete_operations)==1) then
     allocate(psifscf_tmp(-nb(1):2*n_old(1)+1+nb(1),-nb(2):2*n_old(2)+1+nb(2),-nb(3):2*n_old(3)+1+nb(3)+ndebug),stat=i_stat)
     call memocc(i_stat,psifscf_tmp,'psifscf_tmp',subname)
     call switch_axes(nd,nrange,y_phi,centre_old,hgridsh_old,(2*n_old+2+2*nb),psifscfold,&
          hgridsh,ndims_tmp,psifscf_tmp,frag_trans%discrete_operations(1),da)
  else if (size(frag_trans%discrete_operations)==2) then
     stop 'only 1 discrete operation allowed right now'
  else if (size(frag_trans%discrete_operations)==3) then
     stop 'only 1 discrete operation allowed right now'
  end if

  !call field_rototranslation(nd,nrange,y_phi,da,frag_trans%rot_axis,centre_old,centre_new,frag_trans%theta,&
  !     hgridsh_old,ndims_tmp,psifscf_tmp,hgridsh,(2*n+2+2*nb),psifscf)
  call field_rototranslation3D(nd+1,nrange,y_phi,da,frag_trans%rot_axis,centre_old,centre_new,frag_trans%theta,&
       hgridsh_old,ndims_tmp,psifscf_tmp,hgridsh,(2*n+2+2*nb),psifscf)

  if (size(frag_trans%discrete_operations)>0) then
     i_all=-product(shape(psifscf_tmp))*kind(psifscf_tmp)
     deallocate(psifscf_tmp,stat=i_stat)
     call memocc(i_stat,i_all,'psifscf_tmp',subname)
  else
     nullify(psifscf_tmp)    
  end if

  i_all=-product(shape(y_phi))*kind(y_phi)
  deallocate(y_phi,stat=i_stat)
  call memocc(i_stat,i_all,'y_phi',subname)

  !!print*, 'norm of psifscf ',dnrm2((2*n(1)+16)*(2*n(2)+16)*(2*n(3)+16),psifscf,1)
  i_all=-product(shape(psifscfold))*kind(psifscfold)
  deallocate(psifscfold,stat=i_stat)
  call memocc(i_stat,i_all,'psifscfold',subname)
  allocate(psig(0:n(1),2,0:n(2),2,0:n(3),2+ndebug),stat=i_stat)
  call memocc(i_stat,psig,'psig',subname)
  allocate(ww((2*n(1)+2+2*nb(1))*(2*n(2)+2+2*nb(2))*(2*n(3)+2+2*nb(3))+ndebug),stat=i_stat)
  call memocc(i_stat,ww,'ww',subname)

  if (geocode=='F') then
     call analyse_shrink(n(1),n(2),n(3),ww,psifscf,psig)
  else if (geocode == 'S') then
     call analyse_slab(n(1),n(2),n(3),ww,psifscf,psig)
  else if (geocode == 'P') then
     call analyse_per(n(1),n(2),n(3),ww,psifscf,psig)
  end if

  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf',subname)

  !!print*, 'norm new psig ',dnrm2(8*(n(1)+1)*(n(2)+1)*(n(3)+1),psig,1),n(1),n(2),n(3)
  call compress_plain(n(1),n(2),0,n(1),0,n(2),0,n(3),  &
       wfd%nseg_c,wfd%nvctr_c,wfd%keygloc(1,1),wfd%keyvloc(1),   &
       wfd%nseg_f,wfd%nvctr_f,&
       wfd%keygloc(1,wfd%nseg_c+min(1,wfd%nseg_f)),&
       wfd%keyvloc(wfd%nseg_c+min(1,wfd%nseg_f)),   &
       psig,psi(1),psi(wfd%nvctr_c+min(1,wfd%nvctr_f)))
  !!print*, 'norm of reformatted psi ',dnrm2(wfd%nvctr_c+7*wfd%nvctr_f,psi,1),wfd%nvctr_c,wfd%nvctr_f
  !!print*, 'norm of reformatted psic ',dnrm2(wfd%nvctr_c,psi,1)
  !!print*, 'norm of reformatted psif ',dnrm2(wfd%nvctr_f*7,psi(wfd%nvctr_c+min(1,wfd%nvctr_f)),1)

  i_all=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=i_stat)
  call memocc(i_stat,i_all,'psig',subname)
  i_all=-product(shape(ww))*kind(ww)
  deallocate(ww,stat=i_stat)
  call memocc(i_stat,i_all,'ww',subname)

END SUBROUTINE reformat_one_supportfunction


!call the routine which performs the interpolation in each direction
subroutine morph_and_transpose(t0_field,nphi,nrange,phi,ndat,nin,psi_in,nout,psi_out)
 use module_base
 implicit none
 integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
 integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
 integer, intent(in) :: nin,nout !< sizes of the input and output array in interpolating direction
 integer, intent(in) :: ndat !< size of the array in orthogonal directions
! real(gp), intent(in) :: h !< grid spacing in the interpolating direction
 real(gp), dimension(nin,ndat), intent(in) :: t0_field !< field of shifts to be applied for each point in grid spacing units
 real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
 real(gp), dimension(nin,ndat), intent(in) :: psi_in !< input wavefunction psifscf
 real(gp), dimension(ndat,nout), intent(out) :: psi_out !< input wavefunction psifscf
 !local variables
 character(len=*), parameter :: subname='morph_and_transpose'
 real(gp), parameter  :: tol=1.e-14_gp
 integer :: i_all,i_stat,nunit,m_isf,ish,ipos,i,j,l,ms,me,k2,k1
 real(gp) :: dt,tt,t0_l,ksh1,ksh2,k,kold,alpha,diff
 real(gp), dimension(:), allocatable :: shf !< shift filter

 !assume for the moment that the grid spacing is constant
 !call f_malloc_routine_id(subname)
 m_isf=nrange/2

 !shf=f_malloc(bounds=(/-m_isf .to. m_isf/),id='shf')

 !calculate the shift filter for the given t0
 allocate(shf(-m_isf:m_isf+ndebug),stat=i_stat )
 call memocc(i_stat,shf,'shf',subname)

 !number of points for a unit displacement
 nunit=nphi/nrange 

 !apply the interpolating filter to the output
 do j=1,ndat
    psi_out(j,:)=0.0_gp
    do i=1,nout

       !find inverse
       call find_inverse(nin,i,t0_field(1,j),t0_l,k1)
!!$       kold=-1000.0_gp
!!$       find_trans: do l=1,nin
!!$          k=real(l,gp)+t0_field(l,j)
!!$          if (k-real(i,gp) > tol) exit find_trans
!!$          kold=k
!!$          !print *,l,k,t0_field(l,j),i
!!$       end do find_trans
!!$       ! want to use either l or l-1 to give us point i - pick closest
!!$       if (k-real(i,gp) < -kold+real(i,gp)) then
!!$          ksh1=k-real(i,gp)
!!$          ksh2=-kold+real(i,gp)
!!$          k1=l
!!$          k2=l-1
!!$          if (k2==0) then
!!$             k2=1
!!$             ksh2=ksh1
!!$          end if
!!$          if (k1==nin+1) then
!!$             k1=nin
!!$             ksh1=ksh2
!!$          end if
!!$       else
!!$          ksh1=-kold+real(i,gp)
!!$          ksh2=k-real(i,gp)
!!$          k1=l-1
!!$          k2=l
!!$          if (k1==0) then
!!$             k1=1
!!$             ksh1=ksh2
!!$          end if
!!$          if (k2==nin+1) then
!!$             k2=nin
!!$             ksh2=ksh1
!!$          end if
!!$       end if
!!$
!!$       if (ksh1==0.0_gp .or. k1==k2) then !otherwise already have exactly on point
!!$          ksh2=1.0_gp
!!$          ksh1=0.0_gp
!!$       end if 
!!$
!!$       alpha=ksh2/(ksh1+ksh2)
!!$
!!$       t0_l=alpha*t0_field(k1,j)+(1.0_gp-alpha)*t0_field(k2,j)
!!$       !end of find_inverse


       dt=t0_l-nint(t0_l)
   
       diff=real(i,gp)-(k1+t0_l)   

       if (abs(diff - dt) < abs(diff+dt)) dt=-dt

       !define filter for the interpolation starting from a constant shift
       call define_filter(dt,nrange,nphi,phi,shf)

       !here the boundary conditions have to be considered
       tt=0.0_gp
       ms=-min(m_isf,k1-1)
       me=min(m_isf,nin-k1)
       do l=ms,me
          tt=tt+shf(l)*psi_in(k1+l,j)
       end do
       !end of interpolate coefficient


       if (i > 0 .and. i < nout) psi_out(j,i)=tt

    end do
 end do

! call f_free(shf)
! call f_malloc_free_routine()

 i_all=-product(shape(shf))*kind(shf)
 deallocate(shf,stat=i_stat)
 call memocc(i_stat,i_all,'shf',subname)

end subroutine morph_and_transpose

subroutine define_filter(dt,nrange,nphi,phi,shf)
  use module_base
  implicit none
  integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
  integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
  real(gp), intent(in) :: dt
  real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
  real(gp), dimension(-nrange/2:nrange/2), intent(out) :: shf !< interpolating filter to be applied
  !local variables
  integer :: nunit,ish,ipos,m_isf,l,jisf
  
  m_isf=nrange/2
  !number of points for a unit displacement
  nunit=nphi/nrange 
  
  !evaluate the shift
  ish=nint(real(nunit,gp)*dt)

  !starting point in the filter definition
  ipos=ish!+1
  if (ish<= 0) then
     jisf=-(abs(ish))/nunit-1
  else if (ish > 0) then
     jisf=ish/nunit+1
  else
     jisf=0
  end if
  jisf=jisf-m_isf

  !fill the filters in its nonzero coefficients
  do l=-m_isf,m_isf
     if (jisf >= -m_isf .and. jisf <= m_isf) then
        shf(l)=phi(ipos)
     else
        shf(l)=0.0_gp
     end if
     jisf=jisf+1
     ipos=ipos+nunit
  end do

!!$  !old method
!!$  if (ish<=0) then
!!$     shf(-m_isf)=0.0_gp
!!$  else
!!$     shf(-m_isf)=phi(ish)  
!!$  end if
!!$  ipos=ish
!!$  
!!$  do l=-m_isf+1,m_isf-1 !extremes excluded
!!$     !position of the shifted argument in the phi array
!!$     ipos=ipos+nunit
!!$     shf(l)=phi(ipos)  
!!$  end do
!!$  
!!$  if (ish<=0) then
!!$     shf(m_isf)=phi(ipos+nunit)
!!$  else
!!$     shf(m_isf)=0.0_gp
!!$  end if

end subroutine define_filter

!> given a translation vector, find the inverse one
subroutine find_inverse(nin,iout,t0_field,t0_l,k1)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: iout !< point of the new grid from which the inverse has to be found
  integer, intent(in) :: nin !< number of points of the input grid
  real(gp), dimension(nin), intent(in) :: t0_field !<array of displacements of the input grid
  integer, intent(out) :: k1 !< starting point of the input grid from which the interplation should be calculated
  real(gp), intent(out) :: t0_l !< resulting shift from the starting point, from which the filter has to be calculated
  !local variables
  real(gp), parameter  :: tol=1.e-14_gp
  integer :: j,l,ms,me,k2
  real(gp) :: dt,tt,ksh1,ksh2,k,kold,alpha,diff

  kold=-1000.0_gp
  find_trans: do l=1,nin
     k=real(l,gp)+t0_field(l)
     if (k-real(iout,gp) > tol) exit find_trans
     kold=k
  end do find_trans
  ! want to use either l or l-1 to give us point i - pick closest
  if (k-real(iout,gp) < -kold+real(iout,gp)) then
     ksh1=k-real(iout,gp)
     ksh2=-kold+real(iout,gp)
     k1=l
     k2=l-1
     if (k2==0) then
        k2=1
        ksh2=ksh1
     end if
     if (k1==nin+1) then
        k1=nin
        ksh1=ksh2
     end if
  else
     ksh1=-kold+real(iout,gp)
     ksh2=k-real(iout,gp)
     k1=l-1
     k2=l
     if (k1==0) then
        k1=1
        ksh1=ksh2
     end if
     if (k2==nin+1) then
        k2=nin
        ksh2=ksh1
     end if
  end if

  if (ksh1==0.0_gp .or. k1==k2) then !otherwise already have exactly on point
     ksh2=1.0_gp
     ksh1=0.0_gp
  end if

  alpha=ksh2/(ksh1+ksh2)

  t0_l=alpha*t0_field(k1)+(1.0_gp-alpha)*t0_field(k2)

end subroutine find_inverse

subroutine my_scaling_function4b2B(itype,nd,nrange,a,x)
   use module_base
   implicit none
   !Arguments
   !Type of interpolating functions
   integer, intent(in) :: itype
   !Number of points: must be 2**nex
   integer, intent(in) :: nd
   integer, intent(out) :: nrange
   real(kind=8), dimension(0:nd), intent(out) :: a,x
   !Local variables
   character(len=*), parameter :: subname='scaling_function4b2B'
   real(kind=8), dimension(:), allocatable :: y
   integer :: i,nt,ni,i_all,i_stat  

   !Only itype=8,14,16,20,24,30,40,50,60,100
   select case(itype)
   case(8,14,16,20,24,30,40,50,60,100)
      !O.K.
   case default
      print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
      stop
   end select
   !!$  write(unit=*,fmt="(1x,a,i0,a)") &
   !!$       "Use interpolating scaling functions of ",itype," order"

   !Give the range of the scaling function
   !from -itype to itype
   ni=2*itype
   nrange = ni
   allocate(y(0:nd+ndebug),stat=i_stat)
   call memocc(i_stat,y,'y',subname)

   ! plot scaling function
   call zero(nd+1,x)
   call zero(nd+1,y)
   nt=ni
   x(nt/2)=1.d0
   loop1: do
      nt=2*nt
      ! write(6,*) 'nd,nt',nd,nt
      select case(itype)
      case(8)
         stop
      case(14)
         stop
      case(16)
         call back_trans_16(nd,nt,x,y)
      case(20)
         stop
      case(24)
         stop
      case(30)
         stop
      case(40)
         stop
      case(50)
         stop
      case(60)
         stop
      case(100)
         stop
      end select

      do i=0,nt-1
         x(i)=y(i)
      end do
      if (nt.eq.nd) then
         exit loop1
      end if
   end do loop1

   !open (unit=1,file='scfunction',status='unknown')
   do i=0,nd
      a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
      !write(1,*) a(i),x(i)
   end do
   !close(1)

   i_all=-product(shape(y))*kind(y)
   deallocate(y,stat=i_stat)
   call memocc(i_stat,i_all,'y',subname)
END SUBROUTINE my_scaling_function4b2B

subroutine define_rotations(da,newz,boxc_old,boxc_new,theta,hgrids_old,ndims_old,&
     hgrids_new,ndims_new,dx,dy,dz)
  use module_base
  implicit none
  real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
  real(gp), intent(in) :: theta !< rotation wrt newzeta vector
  real(gp), dimension(3), intent(in) :: newz !<coordinates of new z vector (should be of norm one)
  real(gp), dimension(3), intent(in) :: boxc_old,boxc_new !<centre of rotation (the coordinates are zero at these points)
  real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
  integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
  real(gp), dimension(ndims_old(1),ndims_old(2)*ndims_old(3)), intent(out) :: dx !<first deformation
  real(gp), dimension(ndims_old(2),ndims_new(1)*ndims_old(3)), intent(out) :: dy !<second deformation
  real(gp), dimension(ndims_old(3),ndims_new(1)*ndims_new(2)), intent(out) :: dz !<third deformation
  !local variables
  integer :: i,j,k
  real(gp) :: x,y,z,cost,sint,onemc,dotp

  !save recalculation
  sint=sin(theta)
  cost=cos(theta)
  onemc=1.0_gp-cost

  z=-boxc_old(3)
  do k=1,ndims_old(3)
     y=-boxc_old(2)
     do j=1,ndims_old(2)
        x=-boxc_old(1)
        do i=1,ndims_old(1)
           dotp=newz(1)*x+newz(2)*y+newz(3)*z
!!$           dx(i,j+(k-1)*ndims_old(2)) = ( da(1) - x + dotp*newz(1) - cost*(-x + dotp*newz(1)) &
!!$                                      + sint*(z*newz(2) - y*newz(3)) ) / hgrids_old(1)
           dx(i,j+(k-1)*ndims_old(2)) = xp_xyz(da(1),theta,newz,x,y,z) / hgrids_old(1)
           x=x+hgrids_old(1)

        end do
        !if (j==1 .and. k==1) print *,'test vector',dx(:,j+(k-1)*ndims_old(2)),'x,y,z',x,y,z,theta,sint,cost,'newz',newz
        !fill the second vector
        x=-boxc_new(1)
        do i=1,ndims_new(1)
           dy(j,k+(i-1)*ndims_old(3)) = yp_xpyz(da(2),theta,newz,x,y,z) / hgrids_old(2)
!!$           dy(j,k+(i-1)*ndims_old(3)) = ( da(2) - y &
!!$                                      + y*(cost + onemc*newz(2)**2) + z*(-(sint*newz(1)) &
!!$                                      + onemc*newz(2)*newz(3)) - ((onemc*newz(1)*newz(2) &
!!$                                      + sint*newz(3))*(-x  + sint*(z*newz(2) - y*newz(3)) &
!!$                                      + onemc*newz(1)*(y*newz(2) + z*newz(3)))) &
!!$                                      / (cost + onemc*newz(1)**2) ) / hgrids_old(2)

           x=x+hgrids_new(1)
        end do
        y=y+hgrids_old(2)
     end do
     !fill the third vector
     y=-boxc_new(2)
     do j=1,ndims_new(2)
        x=-boxc_new(1)
        do i=1,ndims_new(1)
           dz(k,i+(j-1)*ndims_new(1)) = zp_xpypz(da(3),theta,newz,x,y,z) / hgrids_old(3)
!!$           dz(k,i+(j-1)*ndims_new(1)) = ( da(3) - z &
!!$                                      + (sint*(newz(1)*y - newz(2)*x) - onemc*newz(3)*(newz(1)*x + newz(2)*y) + z) &
!!$                                      / (cost + onemc*newz(3)**2) ) / hgrids_old(3)

           x=x+hgrids_new(1)
        end do
        y=y+hgrids_new(2)
     end do
     z=z+hgrids_old(3)
  end do

contains

  pure function xp_xyz(dx,theta,newz,x,y,z)
    use module_base
    implicit none
    real(gp), intent(in) :: dx,theta,x,y,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: xp_xyz
    !local variables
    real(gp) :: cost,sint,onemc,dotp

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost
    dotp=newz(1)*x+newz(2)*y+newz(3)*z

    xp_xyz=dx-x+dotp*newz(1) - cost*(-x + dotp*newz(1)) + sint*(z*newz(2)-y*newz(3))

  end function xp_xyz

  pure function yp_xpyz(dy,theta,newz,xp,y,z)
    use module_base
    implicit none
    real(gp), intent(in) :: dy,theta,xp,y,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: yp_xpyz
    !local variables
    real(gp) :: cost,sint,onemc

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost

    yp_xpyz=dy - y  + y*(cost + onemc*newz(2)**2) + z*(-(sint*newz(1)) &
         + onemc*newz(2)*newz(3)) - &
         ((onemc*newz(1)*newz(2) + sint*newz(3))*(-xp  + sint*(z*newz(2) - y*newz(3)) &
         + onemc*newz(1)*(y*newz(2) + z*newz(3))))  / (cost + onemc*newz(1)**2)

  end function yp_xpyz

  pure function zp_xpypz(dz,theta,newz,xp,yp,z)
    use module_base
    implicit none
    real(gp), intent(in) :: dz,theta,xp,yp,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: zp_xpypz
    !local variables
    real(gp) :: cost,sint,onemc

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost

    zp_xpypz=dz - z +&
         (sint*(newz(1)*yp - newz(2)*xp) - onemc*newz(3)*(newz(1)*xp + newz(2)*yp) + z) / (cost + onemc*newz(3)**2)

  end function zp_xpypz

end subroutine define_rotations


subroutine define_rotations_switch(da,newz,boxc_old,boxc_new,theta,hgrids_old,ndims_old,&
     hgrids_new,ndims_new,dx,dy,dz)
  use module_base
  implicit none
  real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
  real(gp), intent(in) :: theta !< rotation wrt newzeta vector
  real(gp), dimension(3), intent(in) :: newz !<coordinates of new z vector (should be of norm one)
  real(gp), dimension(3), intent(in) :: boxc_old,boxc_new !<centre of rotation (the coordinates are zero at these points)
  real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
  integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
  real(gp), dimension(ndims_old(1),ndims_old(2)*ndims_old(3)), intent(out) :: dx !<first deformation
  real(gp), dimension(ndims_old(2),ndims_new(1)*ndims_old(3)), intent(out) :: dy !<second deformation
  real(gp), dimension(ndims_old(3),ndims_new(1)*ndims_new(2)), intent(out) :: dz !<third deformation
  !local variables
  integer :: i,j,k
  real(gp) :: x,y,z,cost,sint,onemc,dotp

  !save recalculation
  sint=sin(theta)
  cost=cos(theta)
  onemc=1.0_gp-cost

  z=-boxc_old(3)
  do k=1,ndims_old(3)
     y=-boxc_old(2)
     do j=1,ndims_old(2)
        x=-boxc_old(1)
        do i=1,ndims_old(1)
           dotp=newz(1)*x+newz(2)*y+newz(3)*z
           dx(i,j+(k-1)*ndims_old(2)) = ( da(1) - x + x*newz(1)**2 - y*newz(1)*newz(2) + sint*y*newz(3) &
                                      - z*(sint*newz(2) + newz(1)*newz(3)) + cost*(x - x*newz(1)**2 &
                                      + y*newz(1)*newz(2) + z*newz(1)*newz(3)) ) / hgrids_old(1)
           x=x+hgrids_old(1)
        end do
        !fill the second vector
        x=-boxc_new(1)
        do i=1,ndims_new(1)
           dy(j,k+(i-1)*ndims_old(3)) = ( -da(2) - y + (cost*y - sint*z*newz(1) - onemc*x*newz(1)*newz(2) &
                                      - newz(3)*(sint*x - onemc*(-(z*newz(2)) + y*newz(3)))) &
                                      /(-cost - onemc*newz(1)**2) ) / hgrids_old(2)
           x=x+hgrids_new(1)
        end do
        y=y+hgrids_old(2)
     end do
     !fill the third vector
     y=-boxc_new(2)
     do j=1,ndims_new(2)
        x=-boxc_new(1)
        do i=1,ndims_new(1)
           dz(k,i+(j-1)*ndims_new(1)) = ( -da(3) - z + (z + sint*(-(y*newz(1)) + x*newz(2)) + onemc*(x*newz(1) &
                                      + y*newz(2))*newz(3))/(-cost - onemc*newz(3)**2) ) / hgrids_old(3)
           x=x+hgrids_new(1)
        end do
        y=y+hgrids_new(2)
     end do
     z=z+hgrids_old(3)
  end do

end subroutine define_rotations_switch

subroutine field_rototranslation(n_phi,nrange_phi,phi_ISF,da,newz,centre_old,centre_new,theta,hgrids_old,ndims_old,f_old,&
     hgrids_new,ndims_new,f_new)
  use module_base
  implicit none
  integer, intent(in) :: n_phi,nrange_phi !< number of points of ISF array and real-space range
  real(gp), intent(in) :: theta !< rotation wrt newzeta vector
  real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
  real(gp), dimension(3), intent(in) :: newz !<coordinates of new z vector (should be of norm one)
  real(gp), dimension(3), intent(in) :: centre_old,centre_new !<centre of rotation
  real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
  integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
  real(gp), dimension(n_phi), intent(in) :: phi_ISF
  real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(in) :: f_old
  real(gp), dimension(ndims_new(1),ndims_new(2),ndims_new(3)), intent(out) :: f_new
  !local variables
  real(gp), dimension(:,:), allocatable :: dx,dy,dz
  real(gp), dimension(:,:), allocatable :: work,work2
  
  call f_routine(id='field_rototranslation')
  
  dx=f_malloc((/ndims_old(1),ndims_old(2)*ndims_old(3)/),id='dx')
  dy=f_malloc((/ndims_old(2),ndims_new(1)*ndims_old(3)/),id='dy')
  dz=f_malloc((/ndims_old(3),ndims_new(1)*ndims_new(2)/),id='dz')
  work =f_malloc(shape(dy),id='work')
  work2=f_malloc(shape(dz),id='work2')

  call define_rotations(da,newz,centre_old,centre_new,theta,hgrids_old,ndims_old,&
       hgrids_new,ndims_new,dx,dy,dz)

  !call define_rotations(da,newz,centre_old,centre_new,0.498_gp*pi_param,hgrids_old,ndims_old,&
  !     hgrids_new,ndims_new,dx,dy,dz)


  !call define_rotations_switch(da,newz,centre_old,centre_new,theta,hgrids_old,ndims_old,&
  !     hgrids_new,ndims_new,dx,dy,dz)
  
  !perform interpolation
  call morph_and_transpose(dx,n_phi,nrange_phi,phi_ISF,ndims_old(2)*ndims_old(3),&
         ndims_old(1),f_old,ndims_new(1),work)

  call morph_and_transpose(dy,n_phi,nrange_phi,phi_ISF,ndims_new(1)*ndims_old(3),&
         ndims_old(2),work,ndims_new(2),work2)

  call morph_and_transpose(dz,n_phi,nrange_phi,phi_ISF,ndims_new(1)*ndims_new(2),&
         ndims_old(3),work2,ndims_new(3),f_new)
  
  call f_free(dx)
  call f_free(dy)
  call f_free(dz)
  call f_free(work)
  call f_free(work2)

  call f_release_routine()
end subroutine field_rototranslation

!> routine which directly applies the 3D transformation of the rototranslation
subroutine field_rototranslation3D(n_phi,nrange_phi,phi_ISF,da,newz,centre_old,centre_new,theta,hgrids_old,ndims_old,f_old,&
     hgrids_new,ndims_new,f_new)
  use module_base
  implicit none
  integer, intent(in) :: n_phi,nrange_phi !< number of points of ISF array and real-space range
  real(gp), intent(in) :: theta !< rotation wrt newzeta vector
  real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
  real(gp), dimension(3), intent(in) :: newz !<coordinates of new z vector (should be of norm one)
  real(gp), dimension(3), intent(in) :: centre_old,centre_new !<centre of rotation
  real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
  integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
  real(gp), dimension(n_phi), intent(in) :: phi_ISF
  real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(in) :: f_old
  real(gp), dimension(ndims_new(1),ndims_new(2),ndims_new(3)), intent(out) :: f_new
  !local variables
  integer :: m_isf
  real(gp), dimension(:), allocatable :: shf
!  real(gp), dimension(:,:), allocatable :: dx,dy,dz
  real(gp), dimension(:,:), allocatable :: work,work2
  
!  print *,'3d'

  call f_routine(id='field_rototranslation3D')
  
!!$  dx=f_malloc((/ndims_old(1),ndims_old(2)*ndims_old(3)/),id='dx')
!!$  dy=f_malloc((/ndims_old(2),ndims_new(1)*ndims_old(3)/),id='dy')
!!$  dz=f_malloc((/ndims_old(3),ndims_new(1)*ndims_new(2)/),id='dz')
  work =f_malloc((/ndims_old(2),ndims_new(1)*ndims_old(3)/),id='work')
  work2=f_malloc((/ndims_old(3),ndims_new(1)*ndims_new(2)/),id='work2')

  m_isf=nrange_phi/2
  shf=f_malloc(-m_isf .to. m_isf,id='shf')
  !for each of the dimensions build the interpolating vector which is needed
  

!!$  call define_rotations(da,newz,centre_old,centre_new,theta,hgrids_old,ndims_old,&
!!$       hgrids_new,ndims_new,dx,dy,dz)

  !call define_rotations(da,newz,centre_old,centre_new,0.498_gp*pi_param,hgrids_old,ndims_old,&
  !     hgrids_new,ndims_new,dx,dy,dz)


  !call define_rotations_switch(da,newz,centre_old,centre_new,theta,hgrids_old,ndims_old,&
  !     hgrids_new,ndims_new,dx,dy,dz)

  !the condition for activating the switch can be generalized further
  if (theta==0.0_gp .and. .false.) then
     call shift_only(ndims_old(2)*ndims_old(3),ndims_old(1),&
          ndims_new(1),da(1)/hgrids_old(1),nrange_phi,n_phi,phi_ISF,shf,&
          f_old,work)
     call shift_only(ndims_new(1)*ndims_old(3),ndims_old(2),&
          ndims_new(2),da(2)/hgrids_old(2),nrange_phi,n_phi,phi_ISF,shf,&
          work,work2)
     call shift_only(ndims_new(1)*ndims_new(2),ndims_old(3),&
          ndims_new(3),da(3)/hgrids_old(3),nrange_phi,n_phi,phi_ISF,shf,&
          work2,f_new)
  else
!print *,'x'
     !perform interpolation
     call interpolate_xp_from_x(ndims_new(1),ndims_old(2),ndims_old(3),&
          centre_old(1),centre_old(2),centre_old(3),centre_new(1),&
          hgrids_old(1),hgrids_old(2),hgrids_old(3),hgrids_new(1),&
          ndims_old(1),da(1),theta,newz,nrange_phi,n_phi,phi_ISF,shf,&
          f_old,work)

!!$  call morph_and_transpose(dx,n_phi,nrange_phi,phi_ISF,ndims_old(2)*ndims_old(3),&
!!$         ndims_old(1),f_old,ndims_new(1),work)
!print *,'y'
     call interpolate_yp_from_y(ndims_new(1),ndims_new(2),ndims_old(3),&
          centre_new(1),centre_old(2),centre_old(3),centre_new(2),&
          hgrids_new(1),hgrids_old(2),hgrids_old(3),hgrids_new(2),&
          ndims_old(2),da(2),theta,newz,nrange_phi,n_phi,phi_ISF,shf,&
          work,work2)

!!$  call morph_and_transpose(dy,n_phi,nrange_phi,phi_ISF,ndims_new(1)*ndims_old(3),&
!!$         ndims_old(2),work,ndims_new(2),work2)
!print *,'z'
     call interpolate_zp_from_z(ndims_new(1),ndims_new(2),ndims_new(3),&
          centre_new(1),centre_new(2),centre_old(3),centre_new(3),&
          hgrids_new(1),hgrids_new(2),hgrids_old(3),hgrids_new(3),&
          ndims_old(3),da(3),theta,newz,nrange_phi,n_phi,phi_ISF,shf,&
          work2,f_new)
  end if
!!$  call morph_and_transpose(dz,n_phi,nrange_phi,phi_ISF,ndims_new(1)*ndims_new(2),&
!!$         ndims_old(3),work2,ndims_new(3),f_new)
  
!!$  call f_free(dx)
!!$  call f_free(dy)
!!$  call f_free(dz)
  call f_free(work)
  call f_free(work2)
  call f_free(shf)
  call f_release_routine()

end subroutine field_rototranslation3D

subroutine shift_only(ndat,nin,nout,t0,nrange,nphi,phi,shf,&
     psi_in,psi_out)
  use module_base
  implicit none
  integer, intent(in) :: ndat,nin,nout
  real(gp), intent(in) :: t0
  integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
  integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
  real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
  real(gp), dimension(-nrange/2:nrange/2), intent(inout) :: shf !< interpolating
  real(gp), dimension(nin,ndat), intent(in) :: psi_in
  real(gp), dimension(ndat,nout), intent(out) :: psi_out
  !local variables
  integer :: i,j,ms,me,l,ish,m_isf
  real(gp) :: dt,tt

  m_isf=nrange/2  
  dt=t0-nint(t0)
  !define the shift for output results
  ish=nint(t0)

  !define filter for the interpolation starting from a constant shift
  call define_filter(dt,nrange,nphi,phi,shf)
  !apply the interpolating filter to the output
  do j=1,ndat
     psi_out(j,:)=0.0_gp
     do i=1,nin
        !here the boundary conditions have to be considered
        tt=0.0_gp
        ms=-min(m_isf,i-1)
        me=min(m_isf,nin-i)
        do l=ms,me
           tt=tt+shf(l)*psi_in(i+l,j)
        end do
        !      tt=h*tt

        if (i+ish > 0 .and. i+ish < nout) psi_out(j,i+ish)=tt
        !	      write(102,*)i+ish,psi_out(j,i+ish)
     end do
  end do

end subroutine shift_only

subroutine interpolate_xp_from_x(nx,ny,nz,cx,cy,cz,cx_new,hx,hy,hz,hx_new,&
     nx_old,dx,theta,newz,nrange,nphi,phi,shf,&
     psi_in,psi_out)
  use module_base
  implicit none
  integer, intent(in) :: nx,ny,nz,nx_old
  real(gp), intent(in) :: cx,cy,cz,hx,hy,hz,dx,theta,cx_new,hx_new
  real(gp), dimension(3), intent(in) :: newz
  integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
  integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
  real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
  real(gp), dimension(-nrange/2:nrange/2), intent(inout) :: shf !< interpolating filter to be applied
  real(gp), dimension(nx_old,ny,nz), intent(in) :: psi_in 
  real(gp), dimension(ny,nz,nx), intent(out) :: psi_out 
  !local variables
  integer :: i,j,k,ms,me,l,m_isf,k1
  real(gp) :: x,y,z,tt,dt,t0_l,diff
  !to be removed
  !real(gp), dimension(nx_old) :: tx

  m_isf=nrange/2

  z=-cz
  do k=1,nz
     y=-cy
     do j=1,ny
!!$        !fill the translation vector for this slice
!!$        x=-cx
!!$        do i=1,nx_old
!!$           tx(i)=xp_xyz(dx,theta,newz,x,y,z)/hx
!!$           x=x+hx
!!$        end do
        !start interpolation
        !restart x over new coordinates
        x=-cx_new
        do i=1,nx
           !find x of xp,y,z, put it in k1 and give the shift from it
           !this routine can integrate the vector definition above
!!$           call find_inverse(nx_old,i,tx,t0_l,k1)
!!$           print '(a,1pg25.17,5(i4),6(1pg25.17))','final values',t0_l,k1,&
!!$                i,j,k,&
!!$                !value to be rounded
!!$                min(max(1,nint((x_xpyz(theta,newz,x,y,z)+cx+hx)/hx)),nx_old),& 
!!$                !vector to be applied
!!$                (x-x_xpyz(theta,newz,x,y,z)+dx)/hx 
           t0_l=(x-x_xpyz(theta,newz,x,y,z)+dx)/hx 
           k1=min(max(1,nint((x_xpyz(theta,newz,x,y,z)+cx+hx)/hx)),nx_old)

!!$           dt=t0_l-nint(t0_l)
           diff=real(i,gp)-(k1+t0_l)   
!!$           if (abs(diff) < 1.d0 .or. .true.) then
              dt=-diff
!!$           else if (abs(diff - dt) < abs(diff+dt)) then
!!$              !if (abs(diff) < 1.d0) print *,'herex',dt,diff
!!$              dt=-dt
!!$           end if

           !if (abs(abs(diff)-abs(dt)) > 1.d-12) print '(a,1pg25.17,4(i4),6(1pg25.17))','final values',t0_l,k1,&
           !     i,j,k,dt,diff

           !define filter for the interpolation starting from a constant shift
           call define_filter(dt,nrange,nphi,phi,shf)
           !interpolate the result
           tt=0.0_gp
           ms=-min(m_isf,k1-1)
           me=min(m_isf,nx_old-k1)
           do l=ms,me
              tt=tt+shf(l)*psi_in(k1+l,j,k)
           end do
           !end of interpolate coefficient
           psi_out(j,k,i)=tt
           x=x+hx_new
        end do
!stop
        y=y+hy
     end do
     z=z+hz
  end do

contains

  pure function xp_xyz(dx,theta,newz,x,y,z)
    use module_base
    implicit none
    real(gp), intent(in) :: dx,theta,x,y,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: xp_xyz
    !local variables
    real(gp) :: cost,sint,onemc,dotp

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost
    dotp=newz(1)*x+newz(2)*y+newz(3)*z

    xp_xyz=dx-x+dotp*newz(1) - cost*(-x + dotp*newz(1)) + sint*(z*newz(2)-y*newz(3))

  end function xp_xyz

  !> inverse function to be put in the scheme
  pure function x_xpyz(theta,newz,xp,y,z)
    use module_base
    implicit none
    real(gp), intent(in) :: theta,xp,y,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: x_xpyz
    !local variables
    real(gp) :: cost,sint,onemc

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost
    x_xpyz= xp+sint*newz(3)*y-sint*newz(2)*z-&
         onemc*newz(1)*(newz(2)*y+newz(3)*z)
    x_xpyz=x_xpyz/(cost+onemc*(newz(1)**2))
!!$    tx_xpyz= dx - xp  - sint*newz(3)*y+sint*newz(2)*z+&
!!$         onemc*newz(1)*(newz(2)*y+newz(3)*z)
!!$    tx_xpyz=tx_xpyz/(onemc*(newz(1)**2-1.0_gp))

  end function x_xpyz


end subroutine interpolate_xp_from_x

 
subroutine interpolate_yp_from_y(nx,ny,nz,cx,cy,cz,cy_new,hx,hy,hz,hy_new,&
     ny_old,dy,theta,newz,nrange,nphi,phi,shf,&
     psi_in,psi_out)
  use module_base
  implicit none
  integer, intent(in) :: nx,ny,nz,ny_old
  real(gp), intent(in) :: cx,cy,cz,hx,hy,hz,dy,theta,cy_new,hy_new
  real(gp), dimension(3), intent(in) :: newz
  integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
  integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
  real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
  real(gp), dimension(-nrange/2:nrange/2), intent(inout) :: shf !< interpolating filter to be applied
  real(gp), dimension(ny_old,nz,nx), intent(in) :: psi_in 
  real(gp), dimension(nz,nx,ny), intent(out) :: psi_out 
  !local variables
  integer :: i,j,k,ms,me,l,m_isf,k1
  real(gp) :: x,y,z,tt,dt,t0_l,diff
  !to be removed
  !real(gp), dimension(ny_old) :: ty

  m_isf=nrange/2
 
  x=-cx
  do i=1,nx
     z=-cz
     do k=1,nz
!!$        !fill the translation vector for this slice
!!$        y=-cy
!!$        do j=1,ny_old
!!$           ty(j)=yp_xpyz(dy,theta,newz,x,y,z)/hy
!!$           y=y+hy
!!$        end do
        !start interpolation
        y=-cy_new
        do j=1,ny
           !find x of xp,y,z, put it in k1 and give the shift from it
           !this routine can integrate the vector definition above
!!$           call find_inverse(ny_old,j,ty,t0_l,k1)
!!$           print '(a,1pg25.17,5(i4),6(1pg25.17))','final values',t0_l,k1,&
!!$                i,j,k,&
!!$                !value to be rounded
!!$                min(max(1,nint((y_xpypz(theta,newz,x,y,z)+cy+hy)/hy)),ny_old),& 
!!$                !vector to be applied
!!$                (y-y_xpypz(theta,newz,x,y,z)+dy)/hy 
           t0_l=(y-y_xpypz(theta,newz,x,y,z)+dy)/hy 
           k1=min(max(1,nint((y_xpypz(theta,newz,x,y,z)+cy+hy)/hy)),ny_old)

!!$           dt=t0_l-nint(t0_l)
           diff=real(j,gp)-(k1+t0_l)   
!!$!           if (abs(diff - dt) < abs(diff+dt)) dt=-dt
!!$           if (abs(diff) < 1.d0 .or. .true.) then
              dt=-diff
!!$           else if (abs(diff - dt) < abs(diff+dt)) then
!!$              !if (abs(diff) < 1.d0) print *,'herey',dt,diff
!!$              dt=-dt
!!$           end if

           !define filter for the interpolation starting from a constant shift
           call define_filter(dt,nrange,nphi,phi,shf)
           !interpolate the result
           tt=0.0_gp
           ms=-min(m_isf,k1-1)
           me=min(m_isf,ny_old-k1)
           do l=ms,me
              tt=tt+shf(l)*psi_in(k1+l,k,i)
           end do
           !end of interpolate coefficient
           psi_out(k,i,j)=tt
           y=y+hy_new
        end do
        z=z+hz
     end do
     x=x+hx
  end do

contains

  pure function yp_xpyz(dy,theta,newz,xp,y,z)
    use module_base
    implicit none
    real(gp), intent(in) :: dy,theta,xp,y,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: yp_xpyz
    !local variables
    real(gp) :: cost,sint,onemc

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost

    yp_xpyz=dy - y  + y*(cost + onemc*newz(2)**2) + z*(-(sint*newz(1)) &
         + onemc*newz(2)*newz(3)) - &
         ((onemc*newz(1)*newz(2) + sint*newz(3))*(-xp  + sint*(z*newz(2) - y*newz(3)) &
         + onemc*newz(1)*(y*newz(2) + z*newz(3))))  / (cost + onemc*newz(1)**2)

  end function yp_xpyz

    pure function y_xpypz(theta,newz,xp,yp,z)
    use module_base
    implicit none
    real(gp), intent(in) :: theta,xp,yp,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: y_xpypz
    !local variables
    real(gp) :: cost,sint,onemc
    
    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost

    y_xpypz=yp-onemc*newz(1)*newz(2)*xp-sint*newz(3)*xp+sint*newz(1)*z-&
         onemc*((1.d0-newz(1)**2)*yp-newz(2)*newz(3)*z)
    y_xpypz=y_xpypz/(cost+onemc*(newz(3)**2))

  end function y_xpypz

end subroutine interpolate_yp_from_y

subroutine interpolate_zp_from_z(nx,ny,nz,cx,cy,cz,cz_new,hx,hy,hz,hz_new,&
     nz_old,dz,theta,newz,nrange,nphi,phi,shf,&
     psi_in,psi_out)
  use module_base
  implicit none
  integer, intent(in) :: nx,ny,nz,nz_old
  real(gp), intent(in) :: cx,cy,cz,hx,hy,hz,dz,theta,cz_new,hz_new
  real(gp), dimension(3), intent(in) :: newz
  integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
  integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
  real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
  real(gp), dimension(-nrange/2:nrange/2), intent(inout) :: shf !< interpolating filter to be applied
  real(gp), dimension(nz_old,nx,ny), intent(in) :: psi_in 
  real(gp), dimension(nx,ny,nz), intent(out) :: psi_out 
  !local variables
  integer :: i,j,k,ms,me,l,m_isf,k1
  real(gp) :: x,y,z,tt,dt,t0_l,diff
  !to be removed
  !real(gp), dimension(nz_old) :: tz

  m_isf=nrange/2

  y=-cy
  do j=1,ny
     x=-cx
     do i=1,nx
!!$        z=-cz
!!$        do k=1,nz_old
!!$        !fill the translation vector for this slice
!!$           tz(k)=zp_xpypz(dz,theta,newz,x,y,z)/hz
!!$           z=z+hz
!!$        end do
        !start interpolation
        z=-cz_new
        do k=1,nz
           !find z of xp,y,z, put it in k1 and give the shift from it
           !this routine can integrate the vector definition above
!!$           call find_inverse(nz_old,k,tz,t0_l,k1)
!!$           print '(a,1pg25.17,5(i4),6(1pg25.17))','final values',t0_l,k1,&
!!$                i,j,k,&
!!$                !value to be rounded
!!$                min(max(1,nint((z_xpypzp(theta,newz,x,y,z)+cz+hz)/hz)),nz_old),& 
!!$                !vector to be applied
!!$                (z-z_xpypzp(theta,newz,x,y,z)+dz)/hz 

           t0_l=(z-z_xpypzp(theta,newz,x,y,z)+dz)/hz 
           k1=min(max(1,nint((z_xpypzp(theta,newz,x,y,z)+cz+hz)/hz)),nz_old)

!!$           dt=t0_l-nint(t0_l)
           diff=real(k,gp)-(k1+t0_l)   
!!$           !if (abs(diff - dt) < abs(diff+dt)) dt=-dt
!!$           if (abs(diff) < 1.d0 .or. .true.) then
              dt=-diff
!!$           else if (abs(diff - dt) < abs(diff+dt)) then
!!$              !if (abs(diff) < 1.d0) print *,'herez',dt,diff
!!$              dt=-dt
!!$           end if

           !define filter for the interpolation starting from a constant shift
           call define_filter(dt,nrange,nphi,phi,shf)
           !interpolate the result
           tt=0.0_gp
           ms=-min(m_isf,k1-1)
           me=min(m_isf,nz_old-k1)
           do l=ms,me
              tt=tt+shf(l)*psi_in(k1+l,i,j)
           end do
           !end of interpolate coefficient
           psi_out(i,j,k)=tt
           z=z+hz_new
        end do
        x=x+hx
     end do
     y=y+hy
  end do

contains

  pure function zp_xpypz(dz,theta,newz,xp,yp,z)
    use module_base
    implicit none
    real(gp), intent(in) :: dz,theta,xp,yp,z
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: zp_xpypz
    !local variables
    real(gp) :: cost,sint,onemc

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost

    zp_xpypz=dz - z +&
         (sint*(newz(1)*yp - newz(2)*xp) - onemc*newz(3)*(newz(1)*xp + newz(2)*yp) + z) / (cost + onemc*newz(3)**2)

  end function zp_xpypz

  pure function z_xpypzp(theta,newz,xp,yp,zp)
    use module_base
    implicit none
    real(gp), intent(in) :: theta,xp,yp,zp
    real(gp), dimension(3), intent(in) :: newz
    real(gp) :: z_xpypzp
    !local variables
    real(gp) :: cost,sint,onemc

    sint=sin(theta)
    cost=cos(theta)
    onemc=1.0_gp-cost

    z_xpypzp=sint*(newz(2)*xp-newz(1)*yp)+onemc*newz(3)*(newz(1)*xp+newz(2)*yp)+&
         (cost+onemc*newz(3)**2)*zp

  end function z_xpypzp

end subroutine interpolate_zp_from_z




subroutine switch_axes(n_phi,nrange_phi,phi_ISF,centre_old,hgrids_old,ndims_old,f_old,&
     hgrids_new,ndims_new,f_new,discrete_op,da_global)
  use module_base
  implicit none
  integer, intent(in) :: n_phi,nrange_phi !< number of points of ISF array and real-space range
  character(len=2) :: discrete_op
  real(gp), dimension(3), intent(inout) :: centre_old
  !real(gp), dimension(3), intent(out) :: centre_new !<centre of rotation
  real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
  integer, dimension(3), intent(in) :: ndims_old !<dimension of old and new box
  integer, dimension(3), intent(out) :: ndims_new !<dimension of old and new box
  real(gp), dimension(3), intent(inout) :: da_global ! shift to be used in final interpolation with non-discrete rotation
  real(gp), dimension(n_phi), intent(in) :: phi_ISF
  real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(in) :: f_old
  real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(out) :: f_new ! in general allocate here to ndims_new, but size doesn't change for now so leave like this

  ! local variables
  integer :: i, ipos
  integer, dimension(3) :: icentre_old, icentre_new
  real(gp), dimension(3) :: da, centre_new
  real(gp), dimension(:,:,:), allocatable :: f_tmp, f_tmp2

  ! For now do in the least efficient, easiest way
  call f_routine(id='switch_axes')

  f_tmp=f_malloc((/ndims_old(1),ndims_old(2),ndims_old(3)/),id='f_tmp')

  !do i=1,size(discrete_ops)
     !theta = 90 => y -> -z, z -> y
     if (discrete_op=='x1') then
        stop '90 degrees not yet treated'

     !theta = 180 => y -> -y, z -> -z
     else if (discrete_op=='x2') then
        !shift so that centre is on a grid point: +1 to account for the fact that point 1 corresponds to x=0.0
        icentre_old(1)=0
        da(1)=0.0_gp
        icentre_old(2)=nint(centre_old(2)/hgrids_old(2))+1
        da(2)=real((icentre_old(2)-1)*hgrids_old(2),gp)-centre_old(2)
        icentre_old(3)=nint(centre_old(3)/hgrids_old(3))+1
        da(3)=real((icentre_old(3)-1)*hgrids_old(3),gp)-centre_old(3)

        icentre_new(1)=0
        icentre_new(2)=ndims_old(2)-icentre_old(2)+1
        icentre_new(3)=ndims_old(3)-icentre_old(3)+1

print*,'da',da

        centre_new = centre_old+da
        call field_rototranslation(n_phi,nrange_phi,phi_ISF,da,(/0.0_gp,0.0_gp,0.0_gp/),centre_old,centre_new,0.0_gp,&
            hgrids_old,ndims_old,f_old,hgrids_old,ndims_old,f_tmp)

        f_tmp2=f_malloc((/ndims_old(1),ndims_old(2),ndims_old(3)/),id='f_tmp2')

        ! only switching within axes, so no change in dimensions here
        ndims_new = ndims_old

        ! switch axes so that y -> -y and z -> -z
        do i=1,ndims_old(2)
           ipos=i-icentre_new(2)
           if (-ipos+icentre_old(2)>=1.and.-ipos+icentre_old(2)<=ndims_old(2)) then
              f_tmp2(:,i,:)=f_tmp(:,-ipos+icentre_old(2),:)
           else
stop 'out of range y'
              f_tmp2(:,i,:)=0.0_gp
           end if
        end do

        do i=1,ndims_old(3)
           ipos=i-icentre_new(3)
           if (-ipos+icentre_old(3)>=1.and.-ipos+icentre_old(3)<=ndims_old(3)) then
              f_new(:,:,i)=f_tmp2(:,:,-ipos+icentre_old(3))
           else
stop 'out of range z'
              f_new(:,:,i)=0.0_gp
           end if
        end do

        ! correct the shift to be done afterwards
        da_global=da_global+da-(icentre_new-icentre_old)*hgrids_old
print*,'n',ndims_old
print*,'centres',da,icentre_new,icentre_old
        ! correct centre_old as well ?!  - check when you have a non-zero rotation after!
        centre_old=centre_old+da+(icentre_new-icentre_old)*hgrids_old

     !theta = -90 => y -> z, z -> -y
     else if (discrete_op=='x3') then
        stop '90 degrees not yet treated'

     !theta = 90 => x -> z, z -> -x
     else if (discrete_op=='y1') then
        stop '90 degrees not yet treated'

     !theta = 180 => x -> -x, z -> -z
     else if (discrete_op=='y2') then

     !theta = -90 => x -> -z, z -> x
     else if (discrete_op=='y3') then
        stop '90 degrees not yet treated'

     !theta = 90 => x -> -y, y -> x
     else if (discrete_op=='z1') then
        stop '90 degrees not yet treated'

     !theta = 180 => x -> -x, y -> -y
     else if (discrete_op=='z2') then

     !theta = -90 => x -> y, y -> -x
     else if (discrete_op=='z3') then
        stop '90 degrees not yet treated'
     end if
  !end do

  call f_free(f_tmp)
  call f_free(f_tmp2)

  call f_release_routine()

end subroutine switch_axes





