!> @file
!!    Routines to create the kernel for Poisson solver
!! @author
!!    Copyright (C) 2002-2011 BigDFT group  (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Allocate a pointer which corresponds to the zero-padded FFT slice needed for
!! calculating the convolution with the kernel expressed in the interpolating scaling
!! function basis. The kernel pointer is unallocated on input, allocated on output.
!! SYNOPSIS
!!    @param geocode  Indicates the boundary conditions (BC) of the problem:
!!              - 'F' free BC, isolated systems.
!!                    The program calculates the solution as if the given density is
!!                    "alone" in R^3 space.
!!              - 'S' surface BC, isolated in y direction, periodic in xz plane                
!!                    The given density is supposed to be periodic in the xz plane,
!!                    so the dimensions in these direction mus be compatible with the FFT
!!                    Beware of the fact that the isolated direction is y!
!!              - 'P' periodic BC.
!!                    The density is supposed to be periodic in all the three directions,
!!                    then all the dimensions must be compatible with the FFT.
!!                    No need for setting up the kernel.
!!              - 'W' Wires BC.
!!                    The density is supposed to be periodic in z direction, 
!!                    which has to be compatible with the FFT.
!!    @param iproc,nproc number of process, number of processes
!!    @param n01,n02,n03 dimensions of the real space grid to be hit with the Poisson Solver
!!    @param itype_scf   order of the interpolating scaling functions used in the decomposition
!!    @param hx,hy,hz grid spacings. For the isolated BC case for the moment they are supposed to 
!!                    be equal in the three directions
!!    @param kernel   pointer for the kernel FFT. Unallocated on input, allocated on output.
!!                    Its dimensions are equivalent to the region of the FFT space for which the
!!                    kernel is injective. This will divide by two each direction, 
!!                    since the kernel for the zero-padded convolution is real and symmetric.
!!
!! @warning
!!    Due to the fact that the kernel dimensions are unknown before the calling, the kernel
!!    must be declared as pointer in input of this routine.
!!    To avoid that, one can properly define the kernel dimensions by adding 
!!    the nd1,nd2,nd3 arguments to the PS_dim4allocation routine, then eliminating the pointer
!!    declaration.
subroutine createKernel(iproc,nproc,geocode,n01,n02,n03,hx,hy,hz,itype_scf,kernel,wrtmsg)
  use module_base, only: ndebug
  use yaml_output
  implicit none
 ! include 'mpif.h'
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,itype_scf,iproc,nproc
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), pointer :: kernel(:)
  logical, intent(in) :: wrtmsg
  !local variables
  character(len=*), parameter :: subname='createKernel'
  integer :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,i_stat
  integer :: jproc,nlimd,nlimk,jfd,jhd,jzd,jfk,jhk,jzk,npd,npk

  call timing(iproc,'PSolvKernel   ','ON')

  if (iproc==0 .and. wrtmsg) write(*,'(1x,a)')&
          '------------------------------------------------------------ Poisson Kernel Creation'

  if (geocode == 'P') then
     
     if (iproc==0 .and. wrtmsg) write(*,'(1x,a)',advance='no')&
          'Poisson solver for periodic BC, no kernel calculation...'
     
     call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

     allocate(kernel(nd1*nd2*nd3/nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel,'kernel',subname)

     call Periodic_Kernel(n1,n2,n3,nd1,nd2,nd3,hx,hy,hz,itype_scf,kernel,iproc,nproc)

     nlimd=n2
     nlimk=n3/2+1

  else if (geocode == 'S') then
     
     if (iproc==0 .and. wrtmsg) write(*,'(1x,a)',advance='no')&
          'Calculating Poisson solver kernel, surfaces BC...'

     !Build the Kernel
     call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

     allocate(kernel(nd1*nd2*nd3/nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel,'kernel',subname)

     !the kernel must be built and scattered to all the processes
     call Surfaces_Kernel(n1,n2,n3,m3,nd1,nd2,nd3,hx,hz,hy,itype_scf,kernel,iproc,nproc)

     !last plane calculated for the density and the kernel
     nlimd=n2
     nlimk=n3/2+1

  else if (geocode == 'F') then

     if (iproc==0 .and. wrtmsg) write(*,'(1x,a)',advance='no')&
          'Calculating Poisson solver kernel, free BC...'

     !Build the Kernel
     call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

     allocate(kernel(nd1*nd2*nd3/nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel,'kernel',subname)

     !the kernel must be built and scattered to all the processes
     call Free_Kernel(n01,n02,n03,n1,n2,n3,nd1,nd2,nd3,hx,hy,hz,itype_scf,iproc,nproc,kernel)

     !last plane calculated for the density and the kernel
     nlimd=n2/2
     nlimk=n3/2+1
     
  else if (geocode == 'W') then

     if (iproc==0 .and. wrtmsg) write(*,'(1x,a)',advance='no')&
          'Calculating Poisson solver kernel, wires BC...'

     call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc)

     allocate(kernel(nd1*nd2*nd3/nproc+ndebug),stat=i_stat)
     call memocc(i_stat,kernel,'kernel',subname)

     call Wires_Kernel(iproc,nproc,n1,n2,n3,nd1,nd2,nd3,hx,hy,hz,itype_scf,kernel)

     nlimd=n2
     nlimk=n3/2+1

  else
     
     !if (iproc==0) 
     write(*,'(1x,a,3a)')'createKernel, geocode not admitted',geocode

     stop
  end if

  if (iproc==0 .and. wrtmsg) then
     write(*,'(a)')'done.'
     !if (geocode /= 'P') then 
        write(*,'(1x,2(a,i0))')&
             'Memory occ. per proc. (Bytes):  Density=',md1*md3*md2/nproc*8,&
             '  Kernel=',nd1*nd2*nd3/nproc*8
     !else
     !   write(*,'(1x,2(a,i0))')&
     !        'Memory occ. per proc. (Bytes):  Density=',md1*md3*md2/nproc*8,&
     !        '  Kernel=',8
     !end if
     write(*,'(1x,a,i0)')&
          '                                Full Grid Arrays=',n01*n02*n03*8
     !print the load balancing of the different dimensions on screen
     if (nproc > 1) then
        write(*,'(1x,a)')&
             'Load Balancing for Poisson Solver related operations:'
        jhd=10000
        jzd=10000
        npd=0
        load_balancing: do jproc=0,nproc-1
           !print *,'jproc,jfull=',jproc,jproc*md2/nproc,(jproc+1)*md2/nproc
           if ((jproc+1)*md2/nproc <= nlimd) then
              jfd=jproc
           else if (jproc*md2/nproc <= nlimd) then
              jhd=jproc
              npd=nint(real(nlimd-(jproc)*md2/nproc,kind=8)/real(md2/nproc,kind=8)*100.d0)
           else
              jzd=jproc
              exit load_balancing
           end if
        end do load_balancing
        write(*,'(1x,a,i3,a)')&
             'LB_density        : processors   0  -',jfd,' work at 100%'
        if (jfd < nproc-1) write(*,'(1x,a,i5,a,i5,1a)')&
             '                    processor     ',jhd,&
             '   works at ',npd,'%'
        if (jhd < nproc-1) write(*,'(1x,a,i5,1a,i5,a)')&
             '                    processors ',&
             jzd,'  -',nproc-1,' work at   0%'
        jhk=10000
        jzk=10000
        npk=0
        if (geocode /= 'P') then
           load_balancingk: do jproc=0,nproc-1
              !print *,'jproc,jfull=',jproc,jproc*nd3/nproc,(jproc+1)*nd3/nproc
              if ((jproc+1)*nd3/nproc <= nlimk) then
                 jfk=jproc
              else if (jproc*nd3/nproc <= nlimk) then
                 jhk=jproc
                 npk=nint(real(nlimk-(jproc)*nd3/nproc,kind=8)/real(nd3/nproc,kind=8)*100.d0)
              else
                 jzk=jproc
                 exit load_balancingk
              end if
           end do load_balancingk
           write(*,'(1x,a,i3,a)')&
                ' LB_kernel        : processors   0  -',jfk,' work at 100%'
           if (jfk < nproc-1) write(*,'(1x,a,i5,a,i5,1a)')&
                '                    processor     ',jhk,&
                '   works at ',npk,'%'
           if (jhk < nproc-1) write(*,'(1x,a,i5,1a,i5,a)')&
                '                    processors ',jzk,'  -',nproc-1,&
                ' work at   0%'
        end if
        write(*,'(1x,a)')&
             'Complete LB per proc.= 1/3 LB_density + 2/3 LB_kernel'
     end if

  end if

  call timing(iproc,'PSolvKernel   ','OF')

END SUBROUTINE createKernel
