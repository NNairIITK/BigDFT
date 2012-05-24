!> @file
!!  Routines to handle posinp files
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Match rotation
program matchrotation
   implicit none
   integer::nat,iat,i,iproc,icycle,ibest
   real(kind=8) :: cell(3),rcm(3),rcminit(3),x(6),f(6),evec(6,6),eval(6),direction(6)
   real(kind=8) :: tt1,tt2,tt3,tt4,tt5,tt6
   real(kind=8) :: dis,dis_sq,dx,dy,dz,dis_sq_min,fnrm,alpha,dis_sq_old,tcpuinit
   real(8), allocatable::rat1(:,:),rat2(:,:),rat1min(:,:),rat2min(:,:)
   character(5), allocatable::an(:),anmin(:)
   character(50)::filename
   logical::cputimeexceeded
   filename='posinp.xyz'
   call readnumberatoms(filename,nat)
   allocate(rat1(3,nat),rat2(3,nat),an(nat),rat1min(3,nat),rat2min(3,nat),anmin(nat))
   call readpositionatoms(nat,filename,cell,rat1,an)
   filename='posinp2.xyz'
   call readpositionatoms(nat,filename,cell,rat2,an)
   !---------------------------------------------------------------------------
   call writepathway(nat,rat1,rat2,an,cell,'pathinp.xyz')
   !---------------------------------------------------------------------------
   call calcenterofmass(nat,rat1,rcm)
   rcminit(1:3)=rcm(1:3) !save it for printing output configuration
   do iat=1,nat
      rat1(1:3,iat)=rat1(1:3,iat)-rcm(1:3)
   enddo
   call calcenterofmass(nat,rat1,rcm)
   !---------------------------------------------------------------------------
   call cpu_time(tcpuinit)
   ibest=0
   !alpha=5.d-3
   alpha=1.d0
   dis_sq_min=1.d10
   do icycle=1,10**6
      call calcenterofmass(nat,rat2,rcm)
      do iat=1,nat
         rat2(1:3,iat)=rat2(1:3,iat)-rcm(1:3)
      enddo
      call calcenterofmass(nat,rat2,rcm)
      !---------------------------------------------------------------------------
      dis=0.d0
      do iat=1,nat
         dx=rat1(1,iat)-rat2(1,iat)
         dy=rat1(2,iat)-rat2(2,iat)
         dz=rat1(3,iat)-rat2(3,iat)
         dis=dis+dx**2+dy**2+dz**2
      enddo
      dis=sqrt(dis)
      !---------------------------------------------------------------------------
      if(icycle>1) call reindex(nat,rat1,rat2)
      !---------------------------------------------------------------------------
      iproc=0
      !x(1:6)=0.d0
      call random_number(x)
      x(1:3)=x(1:3)*3.14d0/2.d0
      !x(4:6)=(x(4:6)-0.5d0)*2.d0
      x(4:6)=0.d0
      direction(1:6)=0.d0
      do i=0,500
         x(1:6)=x(1:6)+direction(1:6)
         call calenergyforces(nat,rat1,rat2,x,dis_sq,f,evec,eval)
         if(i==0) then
            dis_sq_old=dis_sq
            alpha=2.d0*eval(6)
         endif
         fnrm=sqrt(f(1)**2+f(2)**2+f(3)**2+f(4)**2+f(5)**2+f(6)**2)
         !write(*,'(1a,1i5,1es24.15,9es10.2)') 'ROTATION ',i,dis_sq,dis_sq-dis_sq_old,fnrm,alpha,&
         !    eval(1),eval(2),eval(3),eval(4),eval(5),eval(6)
         !write(23,'(i7,i4,6es15.5)') icycle,i,x(1:6)
         tt1=max(abs(eval(1)),abs(eval(2)))*1.d-2
         if(i==10) alpha=max((2.d0**(-2 ))*eval(6),tt1)
         if(i==30) alpha=max((2.d0**(-4 ))*eval(6),tt1)
         if(i==50) alpha=max((2.d0**(-6 ))*eval(6),tt1)
         if(i==70) alpha=max((2.d0**(-8 ))*eval(6),tt1)
         if(i==90) alpha=max((2.d0**(-10))*eval(6),tt1)
         if(fnrm<1.d-2) exit
         !if(i>100) then
         tt1=f(1)*evec(1,1)+f(2)*evec(2,1)+f(3)*evec(3,1)+f(4)*evec(4,1)+f(5)*evec(5,1)+f(6)*evec(6,1)
         tt2=f(1)*evec(1,2)+f(2)*evec(2,2)+f(3)*evec(3,2)+f(4)*evec(4,2)+f(5)*evec(5,2)+f(6)*evec(6,2)
         tt3=f(1)*evec(1,3)+f(2)*evec(2,3)+f(3)*evec(3,3)+f(4)*evec(4,3)+f(5)*evec(5,3)+f(6)*evec(6,3)
         tt4=f(1)*evec(1,4)+f(2)*evec(2,4)+f(3)*evec(3,4)+f(4)*evec(4,4)+f(5)*evec(5,4)+f(6)*evec(6,4)
         tt5=f(1)*evec(1,5)+f(2)*evec(2,5)+f(3)*evec(3,5)+f(4)*evec(4,5)+f(5)*evec(5,5)+f(6)*evec(6,5)
         tt6=f(1)*evec(1,6)+f(2)*evec(2,6)+f(3)*evec(3,6)+f(4)*evec(4,6)+f(5)*evec(5,6)+f(6)*evec(6,6)
         direction(1:6)=               tt1*evec(1:6,1)/sqrt(alpha**2+eval(1)**2)
         direction(1:6)=direction(1:6)+tt2*evec(1:6,2)/sqrt(alpha**2+eval(2)**2)
         direction(1:6)=direction(1:6)+tt3*evec(1:6,3)/sqrt(alpha**2+eval(3)**2)
         direction(1:6)=direction(1:6)+tt4*evec(1:6,4)/sqrt(alpha**2+eval(4)**2)
         direction(1:6)=direction(1:6)+tt5*evec(1:6,5)/sqrt(alpha**2+eval(5)**2)
         direction(1:6)=direction(1:6)+tt6*evec(1:6,6)/sqrt(alpha**2+eval(6)**2)
         !else
         !write(*,*) 'ALIREZA ',i
         !direction(1:6)=f(1:6)/alpha
         !endif
         dis_sq_old=dis_sq
      enddo
      !stop
      call transform(x,nat,rat2)
      !---------------------------------------------------------------------------
      dis=0.d0
      do iat=1,nat
         dx=rat1(1,iat)-rat2(1,iat)
         dy=rat1(2,iat)-rat2(2,iat)
         dz=rat1(3,iat)-rat2(3,iat)
         dis=dis+dx**2+dy**2+dz**2
      enddo
      dis=sqrt(dis)
      !---------------------------------------------------------------------------
      write(*,'(a,i6,es24.15)') 'icycle,current RMSD ',icycle,sqrt(dis_sq/nat)
      !call finalize_energyandforces
      if(dis_sq<dis_sq_min) then
         rat1min(1:3,1:nat)=rat1(1:3,1:nat)
         rat2min(1:3,1:nat)=rat2(1:3,1:nat)
         anmin(1:nat)=an(1:nat)
         dis_sq_min=dis_sq
         ibest=ibest+1
         write(*,'(a,i3,i7,es20.10)') 'new lowest RMSD ',ibest,icycle,sqrt(dis_sq/nat)
         write(filename,'(a8,i3.3,a4)') 'posbest_',ibest,'.xyz'
         call writepositionatoms(nat,filename,cell,rat1,an)
         write(filename,'(a9,i3.3,a4)') 'posbest2_',ibest,'.xyz'
         call writepositionatoms(nat,filename,cell,rat2,an)
      endif
      if(mod(icycle,100)==0) then
         call checkwhetherCPUtimeexceeded(tcpuinit,cputimeexceeded)
         if(cputimeexceeded) exit
      endif
      !exit
   enddo !end of loop over icycle
   !---------------------------------------------------------------------------
   write(*,'(a,1es24.15)') 'FINAL RMSD ',sqrt(dis_sq_min/nat)
   rat1(1:3,1:nat)=rat1min(1:3,1:nat)
   rat2(1:3,1:nat)=rat2min(1:3,1:nat)
   an(1:nat)=anmin(1:nat)
   do iat=1,nat
      rat1(1:3,iat)=rat1(1:3,iat)+rcminit(1:3)
      rat2(1:3,iat)=rat2(1:3,iat)+rcminit(1:3)
   enddo
   filename='posout.xyz'
   call writepositionatoms(nat,filename,cell,rat1,an)
   filename='posout2.xyz'
   call writepositionatoms(nat,filename,cell,rat2,an)
   call writepathway(nat,rat1,rat2,an,cell,'pathout.xyz')
   deallocate(rat1,rat2,an,rat1min,rat2min,anmin)
END PROGRAM matchrotation


subroutine reindex(nat,rat1,rat2)
   implicit none
   integer::nat,iat,itry,mat1,mat2
   real(kind=8) :: rat1(3,nat),rat2(3,nat),tt,dx,dy,dz,dis,disold,dis_init,xyz(1:3)
   real(8), save::ediff=1.d0
   dis=0.d0
   do iat=1,nat
      dx=rat1(1,iat)-rat2(1,iat)
      dy=rat1(2,iat)-rat2(2,iat)
      dz=rat1(3,iat)-rat2(3,iat)
      dis=dis+dx**2+dy**2+dz**2
   enddo
   !write(*,'(a,1es24.15)') 'distance square',dis
   disold=dis
   dis_init=dis
   ediff=ediff*0.5d0
   do itry=1,2000
      call random_number(tt) ; mat1=int(tt*nat)+1
      call random_number(tt) ; mat2=int(tt*nat)+1
      if(mat1<1 .or. mat1>nat) write(*,*) 'ERROR: mat1<1 .or. mat1>nat'
      if(mat2<1 .or. mat2>nat) write(*,*) 'ERROR: mat2<1 .or. mat2>nat'
      if(mat1==mat2) cycle
      dx=rat1(1,mat1)-rat2(1,mat1)
      dy=rat1(2,mat1)-rat2(2,mat1)
      dz=rat1(3,mat1)-rat2(3,mat1)
      dis=dis-(dx**2+dy**2+dz**2)
      dx=rat1(1,mat2)-rat2(1,mat2)
      dy=rat1(2,mat2)-rat2(2,mat2)
      dz=rat1(3,mat2)-rat2(3,mat2)
      dis=dis-(dx**2+dy**2+dz**2)
      xyz(1:3)=rat2(1:3,mat1)
      rat2(1:3,mat1)=rat2(1:3,mat2)
      rat2(1:3,mat2)=xyz(1:3)
      dx=rat1(1,mat1)-rat2(1,mat1)
      dy=rat1(2,mat1)-rat2(2,mat1)
      dz=rat1(3,mat1)-rat2(3,mat1)
      dis=dis+(dx**2+dy**2+dz**2)
      dx=rat1(1,mat2)-rat2(1,mat2)
      dy=rat1(2,mat2)-rat2(2,mat2)
      dz=rat1(3,mat2)-rat2(3,mat2)
      dis=dis+(dx**2+dy**2+dz**2)
      if(.not. dis<disold+ediff) then
         xyz(1:3)=rat2(1:3,mat1)
         rat2(1:3,mat1)=rat2(1:3,mat2)
         rat2(1:3,mat2)=xyz(1:3)
         dis=disold
         ediff=ediff*1.01d0
      else
         ediff=ediff/1.01d0
         !write(*,'(a,i6,1es24.15,2es14.5)') 'MC ',itry,dis,dis-disold,ediff
      endif
      if(dis<dis_init*0.9999d0) exit
      disold=dis
   enddo
END SUBROUTINE reindex


subroutine writepathway(nat,rat1,rat2,an,cell,filename)
   implicit none
   integer::nat,iat,ip
   real(kind=8) :: rat1(3,nat),rat2(3,nat),cell(3),dt,t
   character(5)::an(nat)
   character(*)::filename
   real(8), allocatable::rat_t(:,:)
   allocate(rat_t(3,nat))
   dt=1.d0/100.d0
   open(unit=9,file=filename,status='replace')
   do ip=0,100
      t=ip*dt
      rat_t(1:3,1:nat)=(1.d0-t)*rat1(1:3,1:nat)+t*rat2(1:3,1:nat)
      write(9,*) nat,' angstroemd0'
      !write(9,*) nat,' atomicd0'
      write(9,*) cell(1),cell(2),cell(3)
      do iat=1,nat
         write(9,'(a,1x,3es24.15)') trim(an(iat)),rat_t(1,iat),rat_t(2,iat),rat_t(3,iat)
      enddo
   enddo
   close(9)
   deallocate(rat_t)
END SUBROUTINE writepathway


subroutine calcenterofmass(nat,rat,rcm)
   implicit none
   integer::nat,iat
   real(kind=8) :: rat(3,nat),rcm(3)
   rcm(1:3)=0.d0
   do iat=1,nat
      rcm(1:3)=rcm(1:3)+rat(1:3,iat)
   enddo
   rcm(1:3)=rcm(1:3)/nat
   !write(*,'(a,3es24.15)') 'center of mass ',rcm(1:3)
END SUBROUTINE calcenterofmass


subroutine readnumberatoms(filename,nat)
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(out) :: nat
   open(unit=9,file=filename,status='old')
   read(9,*) nat
   close(9)
   write(*,*) 'nat ',nat
END SUBROUTINE readnumberatoms


subroutine readpositionatoms(nat,filename,cell,rat,an)
   implicit none
   integer :: nat,iat,nat_tmp
   character(len=*) :: filename
   real(kind=8) :: cell(3),rat(3,nat)
   character(len=5) :: an(nat)
   character(len=40) :: units
   open(unit=9,file=filename,status='old')
   read(9,*) nat_tmp,units
   !if(trim(units)/='angstroem' .or. trim(units)/='angstroemd0') &
      !    write(*,*) 'WARNING: length units unknown or different from angstroem'
   read(9,*) 
   cell(1:3)=20.d0
   !read(9,*) cell(1),cell(2),cell(3)
   do iat=1,nat
      read(9,*) an(iat),rat(1,iat),rat(2,iat),rat(3,iat)
   enddo
   close(9)
   write(*,'(a,3es24.15)') 'cell ',cell(1:3)
END SUBROUTINE readpositionatoms


subroutine writepositionatoms(nat,filename,cell,rat,an)
   implicit none
   integer::nat,iat
   character(*)::filename
   real(kind=8) :: cell(3),rat(3,nat)
   character(5)::an(nat)
   open(unit=9,file=filename,status='replace')
   write(9,*) nat,'angstroemd0'
   write(9,*) cell(1),cell(2),cell(3)
   do iat=1,nat
      write(9,'(a,1x,3es24.15)') an(iat),rat(1,iat),rat(2,iat),rat(3,iat)
   enddo
   close(9)
END SUBROUTINE writepositionatoms


subroutine transform(x,nat,rat)
   implicit none
   integer::nat,iat
   real(kind=8) :: x(6),rat(3,nat),xat1,yat1,zat1
   real(kind=8) :: u(3,3),ud(3,3,3),udd(3,3,3,3)
   call buildorthogonal(x,u,ud,udd)
   do iat=1,nat
      xat1=rat(1,iat) ; yat1=rat(2,iat) ; zat1=rat(3,iat)
      rat(1,iat)=u(1,1)*xat1+u(1,2)*yat1+u(1,3)*zat1+x(4)
      rat(2,iat)=u(2,1)*xat1+u(2,2)*yat1+u(2,3)*zat1+x(5)
      rat(3,iat)=u(3,1)*xat1+u(3,2)*yat1+u(3,3)*zat1+x(6)
   enddo
END SUBROUTINE transform


subroutine calenergyforces(nat,rat1,rat2,x,chi,f,evec,eval)
   implicit none
   integer::nat,iat,info
   real(kind=8) :: rat1(3,nat),rat2(3,nat),chi,x(6),f(6),evec(6,6),eval(6)
   real(kind=8) :: u(3,3),ud(3,3,3),udd(3,3,3,3),work(500)
   real(kind=8) :: tt1,tt2,tt3,tt4,xat1,yat1,zat1,xat2,yat2,zat2
   real(kind=8) :: a(3),b1(3),b2(3),b3(3),c11(3),c22(3),c33(3),c12(3),c13(3),c23(3)
   call buildorthogonal(x,u,ud,udd)
   chi=0.d0
   f(1:6)=0.d0
   evec(1:6,1:6)=0.d0
   do iat=1,nat
      xat1=rat1(1,iat) ; yat1=rat1(2,iat) ; zat1=rat1(3,iat)
      xat2=rat2(1,iat) ; yat2=rat2(2,iat) ; zat2=rat2(3,iat)
      !calculation of the function
      a(1)=u(1,1)*xat2+u(1,2)*yat2+u(1,3)*zat2
      a(2)=u(2,1)*xat2+u(2,2)*yat2+u(2,3)*zat2
      a(3)=u(3,1)*xat2+u(3,2)*yat2+u(3,3)*zat2
      tt1=xat1**2+yat1**2+zat1**2+xat2**2+yat2**2+zat2**2
      tt2=x(4)**2+x(5)**2+x(6)**2
      tt3=-2.d0*(xat1*x(4)+yat1*x(5)+zat1*x(6))
      tt4=-2.d0*((xat1-x(4))*a(1)+(yat1-x(5))*a(2)+(zat1-x(6))*a(3))
      chi=chi+tt1+tt2+tt3+tt4
      !calculation of the first derivative
      b1(1)=ud(1,1,1)*xat2+ud(1,2,1)*yat2+ud(1,3,1)*zat2
      b1(2)=ud(2,1,1)*xat2+ud(2,2,1)*yat2+ud(2,3,1)*zat2
      b1(3)=ud(3,1,1)*xat2+ud(3,2,1)*yat2+ud(3,3,1)*zat2
      b2(1)=ud(1,1,2)*xat2+ud(1,2,2)*yat2+ud(1,3,2)*zat2
      b2(2)=ud(2,1,2)*xat2+ud(2,2,2)*yat2+ud(2,3,2)*zat2
      b2(3)=ud(3,1,2)*xat2+ud(3,2,2)*yat2+ud(3,3,2)*zat2
      b3(1)=ud(1,1,3)*xat2+ud(1,2,3)*yat2+ud(1,3,3)*zat2
      b3(2)=ud(2,1,3)*xat2+ud(2,2,3)*yat2+ud(2,3,3)*zat2
      b3(3)=ud(3,1,3)*xat2+ud(3,2,3)*yat2+ud(3,3,3)*zat2
      f(1)=f(1)+2.d0*((xat1-x(4))*b1(1)+(yat1-x(5))*b1(2)+(zat1-x(6))*b1(3))
      f(2)=f(2)+2.d0*((xat1-x(4))*b2(1)+(yat1-x(5))*b2(2)+(zat1-x(6))*b2(3))
      f(3)=f(3)+2.d0*((xat1-x(4))*b3(1)+(yat1-x(5))*b3(2)+(zat1-x(6))*b3(3))
      f(4)=f(4)+2.d0*(-x(4)+xat1-a(1))
      f(5)=f(5)+2.d0*(-x(5)+yat1-a(2))
      f(6)=f(6)+2.d0*(-x(6)+zat1-a(3))
      !calculation of the second derivative
      c11(1)=udd(1,1,1,1)*xat2+udd(1,2,1,1)*yat2+udd(1,3,1,1)*zat2
      c11(2)=udd(2,1,1,1)*xat2+udd(2,2,1,1)*yat2+udd(2,3,1,1)*zat2
      c11(3)=udd(3,1,1,1)*xat2+udd(3,2,1,1)*yat2+udd(3,3,1,1)*zat2
      c22(1)=udd(1,1,2,2)*xat2+udd(1,2,2,2)*yat2+udd(1,3,2,2)*zat2
      c22(2)=udd(2,1,2,2)*xat2+udd(2,2,2,2)*yat2+udd(2,3,2,2)*zat2
      c22(3)=udd(3,1,2,2)*xat2+udd(3,2,2,2)*yat2+udd(3,3,2,2)*zat2
      c33(1)=udd(1,1,3,3)*xat2+udd(1,2,3,3)*yat2+udd(1,3,3,3)*zat2
      c33(2)=udd(2,1,3,3)*xat2+udd(2,2,3,3)*yat2+udd(2,3,3,3)*zat2
      c33(3)=udd(3,1,3,3)*xat2+udd(3,2,3,3)*yat2+udd(3,3,3,3)*zat2
      c12(1)=udd(1,1,1,2)*xat2+udd(1,2,1,2)*yat2+udd(1,3,1,2)*zat2
      c12(2)=udd(2,1,1,2)*xat2+udd(2,2,1,2)*yat2+udd(2,3,1,2)*zat2
      c12(3)=udd(3,1,1,2)*xat2+udd(3,2,1,2)*yat2+udd(3,3,1,2)*zat2
      c13(1)=udd(1,1,1,3)*xat2+udd(1,2,1,3)*yat2+udd(1,3,1,3)*zat2
      c13(2)=udd(2,1,1,3)*xat2+udd(2,2,1,3)*yat2+udd(2,3,1,3)*zat2
      c13(3)=udd(3,1,1,3)*xat2+udd(3,2,1,3)*yat2+udd(3,3,1,3)*zat2
      c23(1)=udd(1,1,2,3)*xat2+udd(1,2,2,3)*yat2+udd(1,3,2,3)*zat2
      c23(2)=udd(2,1,2,3)*xat2+udd(2,2,2,3)*yat2+udd(2,3,2,3)*zat2
      c23(3)=udd(3,1,2,3)*xat2+udd(3,2,2,3)*yat2+udd(3,3,2,3)*zat2
      evec(1,1)=evec(1,1)-2.d0*((xat1-x(4))*c11(1)+(yat1-x(5))*c11(2)+(zat1-x(6))*c11(3))
      evec(2,2)=evec(2,2)-2.d0*((xat1-x(4))*c22(1)+(yat1-x(5))*c22(2)+(zat1-x(6))*c22(3))
      evec(3,3)=evec(3,3)-2.d0*((xat1-x(4))*c33(1)+(yat1-x(5))*c33(2)+(zat1-x(6))*c33(3))
      evec(1,2)=evec(1,2)-2.d0*((xat1-x(4))*c12(1)+(yat1-x(5))*c12(2)+(zat1-x(6))*c12(3))
      evec(1,3)=evec(1,3)-2.d0*((xat1-x(4))*c13(1)+(yat1-x(5))*c13(2)+(zat1-x(6))*c13(3))
      evec(2,3)=evec(2,3)-2.d0*((xat1-x(4))*c23(1)+(yat1-x(5))*c23(2)+(zat1-x(6))*c23(3))
      evec(1,4)=evec(1,4)+2.d0*b1(1)
      evec(1,5)=evec(1,5)+2.d0*b1(2)
      evec(1,6)=evec(1,6)+2.d0*b1(3)
      evec(2,4)=evec(2,4)+2.d0*b2(1)
      evec(2,5)=evec(2,5)+2.d0*b2(2)
      evec(2,6)=evec(2,6)+2.d0*b2(3)
      evec(3,4)=evec(3,4)+2.d0*b3(1)
      evec(3,5)=evec(3,5)+2.d0*b3(2)
      evec(3,6)=evec(3,6)+2.d0*b3(3)
   enddo
   evec(4,4)=2.d0*real(nat,8)
   evec(5,5)=2.d0*real(nat,8)
   evec(6,6)=2.d0*real(nat,8)
   evec(4,5)=0.d0
   evec(4,6)=0.d0
   evec(5,6)=0.d0
   !-------------------------------------------------------------
   evec(2,1)=evec(1,2) ; evec(3,1)=evec(1,3) ; evec(4,1)=evec(1,4); evec(5,1)=evec(1,5) ; evec(6,1)=evec(1,6)
   evec(3,2)=evec(2,3) ; evec(4,2)=evec(2,4) ; evec(5,2)=evec(2,5) ; evec(6,2)=evec(2,6)
   evec(4,3)=evec(3,4) ; evec(5,3)=evec(3,5) ; evec(6,3)=evec(3,6)
   evec(5,4)=evec(4,5) ; evec(6,4)=evec(4,6)
   evec(6,5)=evec(5,6)
   !-------------------------------------------------------------
   !write(21,'(6es15.5)') evec(1,1:6)
   !write(21,'(6es15.5)') evec(2,1:6)
   !write(21,'(6es15.5)') evec(3,1:6)
   !write(21,'(6es15.5)') evec(4,1:6)
   !write(21,'(6es15.5)') evec(5,1:6)
   !write(21,'(6es15.5)') evec(6,1:6)
   !write(21,'(a)') '---------------------------------'
   !stop
   !-------------------------------------------------------------
   call DSYEV('V','L',6,evec,6,eval,work,500,info)
   if(info/=0) stop 'ERROR: DSYEV failed.'
   !write(22,'(6es15.5)') f(1),f(2),f(3),f(4),f(5),f(6)
   !stop
   !-------------------------------------------------------------
   return
END SUBROUTINE calenergyforces


subroutine buildorthogonal(x,u,ud,udd)
   real(kind=8) :: x(3),u(3,3),ud(3,3,3),udd(3,3,3,3)
   real(kind=8) :: cx1,cx2,cx3,sx1,sx2,sx3
   !-------------------------------------------------------------
   cx1=cos(x(1)) ; cx2=cos(x(2)) ; cx3=cos(x(3))
   sx1=sin(x(1)) ; sx2=sin(x(2)) ; sx3=sin(x(3))
   !-------------------------------------------------------------
   u(1,1)=cx2
   u(2,1)=cx3*sx2
   u(3,1)=sx3*sx2
   u(1,2)=-sx2*cx1
   u(2,2)=cx3*cx2*cx1-sx3*sx1
   u(3,2)=sx3*cx2*cx1+cx3*sx1
   u(1,3)=sx2*sx1
   u(2,3)=-cx3*cx2*sx1-sx3*cx1
   u(3,3)=-sx3*cx2*sx1+cx3*cx1
   !-------------------------------------------------------------
   ud(1,1,1)=0.d0
   ud(2,1,1)=0.d0
   ud(3,1,1)=0.d0
   ud(1,2,1)=sx2*sx1
   ud(2,2,1)=-cx3*cx2*sx1-sx3*cx1
   ud(3,2,1)=-sx3*cx2*sx1+cx3*cx1
   ud(1,3,1)=sx2*cx1
   ud(2,3,1)=-cx3*cx2*cx1+sx3*sx1
   ud(3,3,1)=-sx3*cx2*cx1-cx3*sx1
   !-------------------------------------------------------------
   ud(1,1,2)=-sx2
   ud(2,1,2)=cx3*cx2
   ud(3,1,2)=sx3*cx2
   ud(1,2,2)=-cx2*cx1
   ud(2,2,2)=-cx3*sx2*cx1
   ud(3,2,2)=-sx3*sx2*cx1
   ud(1,3,2)=cx2*sx1
   ud(2,3,2)=cx3*sx2*sx1
   ud(3,3,2)=sx3*sx2*sx1
   !-------------------------------------------------------------
   ud(1,1,3)=0.d0
   ud(2,1,3)=-sx3*sx2
   ud(3,1,3)=cx3*sx2
   ud(1,2,3)=0.d0
   ud(2,2,3)=-sx3*cx2*cx1-cx3*sx1
   ud(3,2,3)=cx3*cx2*cx1-sx3*sx1
   ud(1,3,3)=0.d0
   ud(2,3,3)=sx3*cx2*sx1-cx3*cx1
   ud(3,3,3)=-cx3*cx2*sx1-sx3*cx1
   !-------------------------------------------------------------
   udd(1,1,1,1)=0.d0
   udd(2,1,1,1)=0.d0
   udd(3,1,1,1)=0.d0
   udd(1,2,1,1)=sx2*cx1
   udd(2,2,1,1)=-cx3*cx2*cx1+sx3*sx1
   udd(3,2,1,1)=-sx3*cx2*cx1-cx3*sx1
   udd(1,3,1,1)=-sx2*sx1
   udd(2,3,1,1)=cx3*cx2*sx1+sx3*cx1
   udd(3,3,1,1)=sx3*cx2*sx1-cx3*cx1
   !-------------------------------------------------------------
   udd(1,1,2,2)=-cx2
   udd(2,1,2,2)=-cx3*sx2
   udd(3,1,2,2)=-sx3*sx2
   udd(1,2,2,2)=sx2*cx1
   udd(2,2,2,2)=-cx3*cx2*cx1
   udd(3,2,2,2)=-sx3*cx2*cx1
   udd(1,3,2,2)=-sx2*sx1
   udd(2,3,2,2)=cx3*cx2*sx1
   udd(3,3,2,2)=sx3*cx2*sx1
   !-------------------------------------------------------------
   udd(1,1,3,3)=0.d0
   udd(2,1,3,3)=-cx3*sx2
   udd(3,1,3,3)=-sx3*sx2
   udd(1,2,3,3)=0.d0
   udd(2,2,3,3)=-cx3*cx2*cx1+sx3*sx1
   udd(3,2,3,3)=-sx3*cx2*cx1-cx3*sx1
   udd(1,3,3,3)=0.d0
   udd(2,3,3,3)=cx3*cx2*sx1+sx3*cx1
   udd(3,3,3,3)=sx3*cx2*sx1-cx3*cx1
   !-------------------------------------------------------------
   udd(1,1,1,2)=0.d0
   udd(2,1,1,2)=0.d0
   udd(3,1,1,2)=0.d0
   udd(1,2,1,2)=cx2*sx1
   udd(2,2,1,2)=cx3*sx2*sx1
   udd(3,2,1,2)=sx3*sx2*sx1
   udd(1,3,1,2)=cx2*cx1
   udd(2,3,1,2)=cx3*sx2*cx1
   udd(3,3,1,2)=sx3*sx2*cx1
   !-------------------------------------------------------------
   udd(1,1,1,3)=0.d0
   udd(2,1,1,3)=0.d0
   udd(3,1,1,3)=0.d0
   udd(1,2,1,3)=0.d0
   udd(2,2,1,3)=sx3*cx2*sx1-cx3*cx1
   udd(3,2,1,3)=-cx3*cx2*sx1-sx3*cx1
   udd(1,3,1,3)=0.d0
   udd(2,3,1,3)=sx3*cx2*cx1+cx3*sx1
   udd(3,3,1,3)=-cx3*cx2*cx1+sx3*sx1
   !-------------------------------------------------------------
   udd(1,1,2,3)=0.d0
   udd(2,1,2,3)=-sx3*cx2
   udd(3,1,2,3)=cx3*cx2
   udd(1,2,2,3)=0.d0
   udd(2,2,2,3)=-sx3*sx2*cx1
   udd(3,2,2,3)=-cx3*sx2*cx1
   udd(1,3,2,3)=0.d0
   udd(2,3,2,3)=-sx3*sx2*sx1
   udd(3,3,2,3)=cx3*sx2*sx1
   !-------------------------------------------------------------
   udd(1:3,1:3,2,1)=udd(1:3,1:3,1,2)
   udd(1:3,1:3,3,1)=udd(1:3,1:3,1,3)
   udd(1:3,1:3,3,2)=udd(1:3,1:3,2,3)
   !-------------------------------------------------------------
END SUBROUTINE buildorthogonal


subroutine checkwhetherCPUtimeexceeded(tcpuinit,cputimeexceeded)
   implicit none
   logical::cputimeexceeded
   real(kind=8) :: cpulimit,tcpuinit,tcpucurr
   integer::k
   cputimeexceeded=.false.
   open(unit=55,file='CPUlimit',status='unknown')
   read(55,*,iostat=k) cpulimit
   if(k==0) then !k=0 means there was no error, nor was EOF encountered.
      call cpu_time(tcpucurr)
      if(tcpucurr-tcpuinit>cpulimit) then
         write(*,'(a)') 'CPU time exceeded.'
         cputimeexceeded=.true.
      endif
   endif
   close(55)
END SUBROUTINE checkwhetherCPUtimeexceeded
