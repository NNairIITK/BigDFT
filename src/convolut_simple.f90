! Simple non-optimized version of the major convolution routines

        subroutine convolut_kinetic(n1,n2,n3,hgrid,x,y)
!   applies the kinetic energy operator onto x to get y
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-14,lupfil=14)
        dimension fil(lowfil:lupfil),x(0:n1,0:n2,0:n3),y(0:n1,0:n2,0:n3)
        scale=-.5d0/hgrid**2

! second derivative filters for Daubechies 16
         fil(0)=   -3.5536922899131901941296809374d0*scale
         fil(1)=    2.2191465938911163898794546405d0*scale
         fil(2)=   -0.6156141465570069496314853949d0*scale
         fil(3)=    0.2371780582153805636239247476d0*scale
         fil(4)=   -0.0822663999742123340987663521d0*scale
         fil(5)=    0.02207029188482255523789911295638968409d0*scale
         fil(6)=   -0.409765689342633823899327051188315485d-2*scale
         fil(7)=    0.45167920287502235349480037639758496d-3*scale
         fil(8)=   -0.2398228524507599670405555359023135d-4*scale
         fil(9)=    2.0904234952920365957922889447361d-6*scale
         fil(10)=  -3.7230763047369275848791496973044d-7*scale
         fil(11)=  -1.05857055496741470373494132287d-8*scale
         fil(12)=  -5.813879830282540547959250667d-11*scale
         fil(13)=   2.70800493626319438269856689037647576d-13*scale
         fil(14)=  -6.924474940639200152025730585882d-18*scale

	do i=1,14
        fil(-i)=fil(i)
        enddo

	do i3=0,n3

! (1/2) d^2/dx^2
	do i2=0,n2
	do i1=0,n1
        tt=0.d0
        do l=max(-i1,lowfil),min(lupfil,n1-i1)
	tt=tt+x(i1+l,i2,i3)*fil(l)
        enddo
	y(i1,i2,i3)=tt
        enddo
        enddo

! + (1/2) d^2/dy^2
	do i1=0,n1
	do i2=0,n2
        tt=0.d0
        do l=max(-i2,lowfil),min(lupfil,n2-i2)
	tt=tt+x(i1,i2+l,i3)*fil(l)
        enddo
	y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
        enddo

        enddo

! + (1/2) d^2/dz^2
	do i2=0,n2
	do i1=0,n1
	do i3=0,n3
        tt=0.d0
        do l=max(-i3,lowfil),min(lupfil,n3-i3)
	tt=tt+x(i1,i2,i3+l)*fil(l)
        enddo
	y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
        enddo
        enddo

	return
	end


        subroutine convrot_grow(n1,ndat,x,y)
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7)
        dimension x(0:n1,ndat),y(ndat,-lupfil:n1-lowfil)
! the filtered output data structure has grown by the filter length

!          THE MAGIC FILTER FOR DAUBECHIES-16
        REAL*8 fil(lowfil:lupfil)
        DATA fil / &
        8.4334247333529341094733325815816D-7,&
       -0.1290557201342060969516786758559028D-4,&
        0.8762984476210559564689161894116397D-4,&
       -0.30158038132690463167163703826169879D-3,&
        0.174723713672993903449447812749852942D-2,&
       -0.942047030201080385922711540948195075D-2,&
        0.2373821463724942397566389712597274535D-1,&
        0.612625895831207982195380597D-1,&
        0.9940415697834003993178616713D0,&
       -0.604895289196983516002834636D-1, &
       -0.2103025160930381434955489412839065067D-1,&
        0.1337263414854794752733423467013220997D-1,&
       -0.344128144493493857280881509686821861D-2,&
        0.49443227688689919192282259476750972D-3,&
       -0.5185986881173432922848639136911487D-4,&
         2.72734492911979659657715313017228D-6 /



	do j=1,ndat
	do i=-lupfil,n1-lowfil
        
        tt=0.d0
        do l=max(-i,lowfil),min(lupfil,n1-i)
	tt=tt+x(i+l,j)*fil(l)
        enddo
	y(j,i)=tt

        enddo
        enddo

	return
	end


        subroutine convrot_shrink(n1,ndat,x,y)
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7)
        dimension x(lowfil:n1+lupfil,ndat),y(ndat,0:n1)
! the filtered output data structure has shrunk by the filter length

!          THE MAGIC FILTER FOR DAUBECHIES-16
        REAL*8 fil(lowfil:lupfil)
        DATA fil / &
         2.72734492911979659657715313017228D-6,&
       -0.5185986881173432922848639136911487D-4,&
        0.49443227688689919192282259476750972D-3,&
       -0.344128144493493857280881509686821861D-2,&
        0.1337263414854794752733423467013220997D-1,&
       -0.2103025160930381434955489412839065067D-1,&
       -0.604895289196983516002834636D-1,&
        0.9940415697834003993178616713D0,&
        0.612625895831207982195380597D-1,&
        0.2373821463724942397566389712597274535D-1,&
       -0.942047030201080385922711540948195075D-2,&
        0.174723713672993903449447812749852942D-2,&
       -0.30158038132690463167163703826169879D-3,&
        0.8762984476210559564689161894116397D-4,&
       -0.1290557201342060969516786758559028D-4,&
        8.4334247333529341094733325815816D-7 /


	do j=1,ndat
	do i=0,n1
        
        tt=0.d0
        do l=lowfil,lupfil
	tt=tt+x(i+l,j)*fil(l)
        enddo
	y(j,i)=tt

        enddo
        enddo

	return
	end

