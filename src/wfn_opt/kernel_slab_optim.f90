
subroutine hit_with_kernel_slab(x,zx,kern_k1,kern_k3,n1,n2,n3,c,hgrid) 
! hits the input array x with the kernel
! ((-1/2\Delta+C)_{ij})^{-1}
 use module_base
 implicit none
 integer,intent(in) :: n1,n2,n3
 real(gp),intent(in) :: kern_k1(0:n1)
 real(gp),intent(in) :: kern_k3(0:n3)
 real(gp),intent(in) :: c,hgrid

 real(wp),intent(inout) :: x(0:n1,0:n2,0:n3)! input

 real(wp) :: zx(2,0:(n1+1)/2,0:n2,0:n3)! work array
 integer :: i1,i2,i3,isign,inzee

 !*******************************************************************************************
 ! for 2-dimensional FFT:
 real(wp),allocatable :: z(:,:,:,:),zw(:)
 integer,parameter :: ncache=4*1024

 real(wp) :: t1,t2,tela
 integer :: count1,count2,count_rate,count_max
 integer :: nd1,nd2,nd3

nd1=n1+2
nd2=n2+2
nd3=n3+2

allocate(z(2,nd1,nd3,2))
allocate(zw(ncache+1))

call cpu_time(t1)
call system_clock(count1,count_rate,count_max)      

call forward_fft

call segment_invert(n1,n2,n3,kern_k1,kern_k3,c,zx,hgrid)

call backward_fft

call cpu_time(t2)
call system_clock(count2,count_rate,count_max)      
tela=(count2-count1)/real(count_rate,kind(wp))
!write(*,*) 'Time (CPU,ELA)  (sec):' ,t2-t1,tela

deallocate(z,zw)

contains

 subroutine forward_fft
 implicit none
 ! In order to use complex Fourier transforms, we make the ansatz: 
 ! for i2=0..(n2-1)/2,
 ! z(:,:)=x(:,2*i2,:)+i x(:,2*i2+1,:)

 ! We use the following formulas: if z(j1,j2)=x(j1,j2)+i y(j1,j2)
 ! and x,y are real, then for the Fourier tramsforms X(k1,k2),Y(k1,k2) of x
 ! and y and Z=X+iY, one can prove:
 ! Re X(k1,k2)=( Re Z(-k1,-k2)+Re Z(k1,k2))/2
 ! Im X(k1,k2)=(-Im Z(-k1,-k2)+Im Z(k1,k2))/2
 ! Re Y(k1,k2)=( Im Z(-k1,-k2)+Im Z(k1,k2))/2
 ! Im Y(k1,k2)=( Re Z(-k1,-k2)-Re Z(k1,k2))/2
 isign=1 ; inzee=1
 do i2=0,(n2-1)/2
  do i3=0,n3
   do i1=0,n1
    z(1,i1+1,i3+1,inzee)=x(i1,2*i2  ,i3)
    z(2,i1+1,i3+1,inzee)=x(i1,2*i2+1,i3)
   enddo
  enddo
  
  call fft2d(n1+1,n3+1,nd1,nd3,z,isign,inzee,zw,ncache)

  ! i3=0
  zx(1,0,2*i2  ,0)=2._wp*z(1,1,1,inzee)
  zx(2,0,2*i2  ,0)=0._wp ! the zeroth mode coefficient is real for real data
  zx(1,0,2*i2+1,0)=2._wp*z(2,1,1,inzee)
  zx(2,0,2*i2+1,0)=0._wp ! the zeroth mode coefficient is real for real data 

  do i1=1,(n1+1)/2
   zx(1,i1,2*i2  ,0)= z(1,i1+1,1,inzee)+z(1,n1+2-i1,1,inzee)
   zx(2,i1,2*i2  ,0)= z(2,i1+1,1,inzee)-z(2,n1+2-i1,1,inzee)
   zx(1,i1,2*i2+1,0)= z(2,i1+1,1,inzee)+z(2,n1+2-i1,1,inzee)
   zx(2,i1,2*i2+1,0)=-z(1,i1+1,1,inzee)+z(1,n1+2-i1,1,inzee)
  enddo
  ! nonzero i3
  do i3=1,n3
   zx(1,0,2*i2  ,i3)= z(1,1,i3+1,inzee)+z(1,1,n3+2-i3,inzee)
   zx(2,0,2*i2  ,i3)= z(2,1,i3+1,inzee)-z(2,1,n3+2-i3,inzee)
   zx(1,0,2*i2+1,i3)= z(2,1,i3+1,inzee)+z(2,1,n3+2-i3,inzee)
   zx(2,0,2*i2+1,i3)=-z(1,1,i3+1,inzee)+z(1,1,n3+2-i3,inzee)
   do i1=1,(n1+1)/2
    ! real data leads to symmetric real, antisymmetric imaginary
    zx(1,i1,2*i2  ,i3)= z(1,i1+1,i3+1,inzee)+z(1,n1+2-i1,n3+2-i3,inzee)
    zx(2,i1,2*i2  ,i3)= z(2,i1+1,i3+1,inzee)-z(2,n1+2-i1,n3+2-i3,inzee)
    zx(1,i1,2*i2+1,i3)= z(2,i1+1,i3+1,inzee)+z(2,n1+2-i1,n3+2-i3,inzee)
    zx(2,i1,2*i2+1,i3)=-z(1,i1+1,i3+1,inzee)+z(1,n1+2-i1,n3+2-i3,inzee)
   enddo
  enddo
 enddo

 ! if n2 is even, 2((n2-1)/2)+1=n2-1, so i2=n2 is computed separately
 if (2*(n2/2).eq.n2) then

  do i3=0,n3
   do i1=0,n1
    z(1,i1+1,i3+1,inzee)=x(i1,n2  ,i3)*2._wp
    z(2,i1+1,i3+1,inzee)=0._wp
   enddo
  enddo
  
  call fft2d(n1+1,n3+1,nd1,nd3,z,isign,inzee,zw,ncache)

  do i3=0,n3
   do i1=0,(n1+1)/2
    zx(1,i1,n2,i3)=z(1,i1+1,i3+1,inzee)
    zx(2,i1,n2,i3)=z(2,i1+1,i3+1,inzee)
   enddo
  enddo
 endif

 end subroutine forward_fft

 subroutine backward_fft
 implicit none
 ! In order to use complex Fourier transforms, we make the ansatz: 
 ! for i2=0..(n2-1)/2,
 ! z(:,:,)=x(:,2*i2,:)+i x(:,2*i2+1,:)

 ! We use the following formulas: if z(j1,j2)=x(j1,j2)+i y(j1,j2)
 ! and x,y are real, then for the Fourier tramsforms X(k1,k2),Y(k1,k2) of x
 ! and y and Z=X+iY, we have:
 ! Re Z( k1, k2) =  Re X(k1,k2) - Im Y(k1,k2)
 ! Im Z( k1, k2) =  Im X(k1,k2) + Re Y(k1,k2)
 ! Re Z(-k1,-k2) =  Re X(k1,k2) + Im Y(k1,k2)
 ! Im Z(-k1,-k2) = -Im X(k1,k2) + Re Y(k1,k2)
 
 isign=-1 ; inzee=1
 do i2=0,(n2-1)/2
  ! i3=0
  ! zx(i1,i2,0)=zx*(n1+1-i1,i2,0) for i1 != 0
  z(1,1,1,inzee)=zx(1,0,2*i2,0)-zx(2,0,2*i2+1,0)
  z(2,1,1,inzee)=zx(2,0,2*i2,0)+zx(1,0,2*i2+1,0)
  do i1=1,(n1+1)/2
   z(1,i1+1,1,inzee)=zx(1,i1,2*i2,0)-zx(2,i1,2*i2+1,0)
   z(2,i1+1,1,inzee)=zx(2,i1,2*i2,0)+zx(1,i1,2*i2+1,0)
  
   z(1,n1+2-i1,1,inzee)= zx(1,i1,2*i2,0)+zx(2,i1,2*i2+1,0)
   z(2,n1+2-i1,1,inzee)=-zx(2,i1,2*i2,0)+zx(1,i1,2*i2+1,0)
  enddo

  ! nonzero i3
  ! zx(i1,i2,i3)=zx*(n1+1-i1,i2,n3+1-i3) for i1 != 0
  do i3=1,n3

   z(1,1,i3+1,inzee)=zx(1,0,2*i2,i3)-zx(2,0,2*i2+1,i3)
   z(2,1,i3+1,inzee)=zx(2,0,2*i2,i3)+zx(1,0,2*i2+1,i3)
   do i1=1,(n1+1)/2
    z(1,i1+1,i3+1,inzee)=zx(1,i1,2*i2,i3)-zx(2,i1,2*i2+1,i3)
    z(2,i1+1,i3+1,inzee)=zx(2,i1,2*i2,i3)+zx(1,i1,2*i2+1,i3)
   
    z(1,n1+2-i1,n3+2-i3,inzee)= zx(1,i1,2*i2,i3)+zx(2,i1,2*i2+1,i3)
    z(2,n1+2-i1,n3+2-i3,inzee)=-zx(2,i1,2*i2,i3)+zx(1,i1,2*i2+1,i3)
   enddo

  enddo
  

  call fft2d(n1+1,n3+1,nd1,nd3,z,isign,inzee,zw,ncache)

  do i3=0,n3
   do i1=0,n1
    x(i1,2*i2  ,i3)=z(1,i1+1,i3+1,inzee)*.5d0
    x(i1,2*i2+1,i3)=z(2,i1+1,i3+1,inzee)*.5d0
   enddo
  enddo
 enddo

 ! if n2 is even, 2((n2-1)/2)+1=n2-1, so i2=n2 is computed separately
 if (2*(n2/2).eq.n2) then
  ! i3=0
  ! zx(i1,n2,0)=zx*(n1+1-i1,n2,0) for i1 != 0
  z(1,1,1,inzee)=zx(1,0,n2,0)
  z(2,1,1,inzee)=zx(2,0,n2,0)
  do i1=1,(n1+1)/2
   z(1,i1+1,1,inzee)=zx(1,i1,n2,0)
   z(2,i1+1,1,inzee)=zx(2,i1,n2,0)
  
   z(1,n1+2-i1,1,inzee)=zx(1,i1,n2,0)
   z(2,n1+2-i1,1,inzee)=-zx(2,i1,n2,0)
  enddo

  ! nonzero i3
  ! zx(i1,n2,i3)=zx*(n1+1-i1,n2,n3+1-i3) for i1 != 0
  do i3=1,n3

   z(1,1,i3+1,inzee)=zx(1,0,n2,i3)
   z(2,1,i3+1,inzee)=zx(2,0,n2,i3)
   do i1=1,(n1+1)/2
    z(1,i1+1,i3+1,inzee)=zx(1,i1,n2,i3)
    z(2,i1+1,i3+1,inzee)=zx(2,i1,n2,i3)
   
    z(1,n1+2-i1,n3+2-i3,inzee)=zx(1,i1,n2,i3)
    z(2,n1+2-i1,n3+2-i3,inzee)=-zx(2,i1,n2,i3)
   enddo

  enddo
  

  call fft2d(n1+1,n3+1,nd1,nd3,z,isign,inzee,zw,ncache)

  do i3=0,n3
   do i1=0,n1
    x(i1,n2,i3)=z(1,i1+1,i3+1,inzee)*.5d0
   enddo
  enddo
 endif
end subroutine backward_fft


end subroutine hit_with_kernel_slab


