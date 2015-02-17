module module_BornHugginsMayerTosiFumi
    use module_base
    implicit none
    integer, save ::ntypinter
    real(gp), save ::hsp
    real(gp), allocatable, save ::aaa(:),bbb(:),ccc(:),ddd(:),eee(:),fff(:),fdsp(:,:,:),fsp(:,:,:)
    integer, allocatable, save ::typinteraction(:,:)
    integer, allocatable, save :: typat(:)
    logical, save :: initialized = .false.

    private

    public :: init_bhmtf
    public :: energyandforces_bhmtf

contains
subroutine allocateshrtrngpotentialcoeffarrays(ntypat,nsp)
! parameters are read in from file 'parameters.dat'
    implicit none
    integer, intent(in) :: ntypat
    integer::istat,nsp
    ntypinter=(ntypat**2 + ntypat)/2
    aaa = f_malloc((/1 .to. ntypinter/),id='aaa')
    bbb = f_malloc((/1 .to. ntypinter/),id='bbb')
    ccc = f_malloc((/1 .to. ntypinter/),id='ccc')
    ddd = f_malloc((/1 .to. ntypinter/),id='ddd')
    eee = f_malloc((/1 .to. ntypinter/),id='eee')
    fff = f_malloc((/1 .to. ntypinter/),id='fff')
    fdsp = f_malloc((/0 .to. 3, 0 .to. nsp, 1 .to. ntypinter/),id='fdsp')
    fsp = f_malloc((/0 .to. 4, 0 .to. nsp-1, 1 .to. ntypinter/),id='fsp')
    typinteraction = f_malloc((/1 .to. ntypat, 1 .to. ntypat/),id='typinteraction')
end subroutine allocateshrtrngpotentialcoeffarrays
subroutine deallocateshrtrngpotentialcoeffarrays
    implicit none
    integer::istat
    call f_free(aaa) 
    call f_free(bbb)
    call f_free(ccc)
    call f_free(ddd)
    call f_free(eee)
    call f_free(fff)
    call f_free(fdsp)
    call f_free(fsp)
    call f_free(typinteraction)

end subroutine deallocateshrtrngpotentialcoeffarrays

subroutine finalize_bhmtf()
    use module_base
    use yaml_output
    implicit none
    call f_free(typat)
    call deallocateshrtrngpotentialcoeffarrays
end subroutine

subroutine init_bhmtf(nat,astruct,paramset,paramfile,geocode)
    use module_base
    use yaml_output
    use module_atoms
    implicit none
    !parameter
    integer, intent(in) :: nat
    type(atomic_structure), intent(in) :: astruct
    character(len=*), intent(in) :: paramset
    character(len=*), intent(in) :: paramfile
    character(len=*), intent(in) :: geocode
    !internal
    integer, parameter :: nall=3
    integer, parameter :: nsp=5*10**5
    real(gp), parameter :: rcut = sqrt(3.0_gp)*130.0_gp
    integer :: iat
    integer :: iall
    integer :: ityp
    integer :: itypat
    integer :: jtypat
    integer :: itypinter
    integer :: ntypat
    integer :: leni, lenj
    character(len=2) :: namat(100)
    character(20)::strtmpi,strtmpj
    character(2):: namatnamat1,namatnamat2
    real(gp)::parameters(6,nall)
    character(len=4) :: nameinteraction(nall)
    logical :: condition

    initialized=.false.

    if(trim(geocode)/='F')then
        initialized=.false.
        call f_err_throw('BHMTF potential only works with free boundary '//&
             'conditions. Specified boundary conditions are: '//&
             trim(adjustl(geocode)))
    endif

    typat = f_malloc((/ 1 .to. nat/),id='typat')


    ntypat=0
    namat=''
    do iat=1,nat
        condition=.true.
        do ityp=1,ntypat
            if(trim(astruct%atomnames(astruct%iatype(iat)))==namat(ityp))then
                typat(iat)=ityp
                condition=.false.
                exit
            endif
        enddo
        if(condition) then
            ntypat=ntypat+1
            if (ntypat > 100) stop 'more than 100 atomnames not permitted'
            namat(ityp)=trim(astruct%atomnames(astruct%iatype(iat)))
            typat(iat)=ntypat
        endif
    enddo

    call allocateshrtrngpotentialcoeffarrays(ntypat,nsp)
    call preparespline(nsp,rcut)

    if(trim(paramfile)/='none')then                                    
         call f_err_throw('Reading parameters from file not '//&        
              'implemented for BHMTF potential')                                     
    else  

        select case(trim(paramset))
        case('NaCl')
            call yaml_mapping_open('Using NaCL parameters from'//&           
                 ' Adams, D. J., McDonald, I. R. J. Phys. C: Solid State Phys., Vol. 7, 1974')
            nameinteraction(1) = 'NaNa' 
            parameters(1,1) = 15.5672851674_gp !A(Na<->Na) 
            parameters(2,1) = 1.6693287451_gp !1/rho (Na<->Na)
            parameters(3,1) = 1.75485576_gp !C(Na<->Na) 
            parameters(4,1) = 2.9841448_gp !D(Na<->Na) 
            parameters(5,1) = 1.0_gp !q_Na*q_Na
            parameters(6,1) = 0.0_gp !???
            nameinteraction(2) = 'NaCL' 
            parameters(1,2) =  46.1167544593_gp!A(Na<->Cl) 
            parameters(2,2) =  1.6693287451_gp !1/rho (Na<->Cl)
            parameters(3,2) =  11.6990384_gp !C(Na<->Cl) 
            parameters(4,2) =  51.8495159_gp!D(Na<->Cl) 
            parameters(5,2) = -1.0_gp !q_Na*q_Cl
            parameters(6,2) = 0.0_gp !???
            nameinteraction(3) = 'ClCl' 
            parameters(1,3) =  128.0783919804_gp!A(Cl<->Cl) 
            parameters(2,3) =  1.6693287451_gp !1/rho (Cl<->Cl)
            parameters(3,3) =  121.168612_gp!C(Cl<->Cl) 
            parameters(4,3) =  869.132173_gp!D(Cl<->Cl) 
            parameters(5,3) = 1.0_gp !q_Cl*q_Cl
            parameters(6,3) = 0.0_gp !???
            call yaml_map('A(Na<->Na) [Hartree]', parameters(1,1) ,  fmt='(1pe10.4)')           
            call yaml_map('A(Na<->Cl) [Hartree]', parameters(1,2),  fmt='(1pe10.4)')           
            call yaml_map('A(Cl<->Cl) [Hartree]', parameters(1,3),  fmt='(1pe10.4)')           
            call yaml_map('C(Na<->Na) [Bohr^6 Hartree]', parameters(3,1),  fmt='(1pe10.4)')           
            call yaml_map('C(Na<->Cl) [Bohr^6 Hartree]', parameters(3,2),  fmt='(1pe10.4)')           
            call yaml_map('C(Cl<->Cl) [Bohr^6 Hartree]', parameters(3,3),  fmt='(1pe10.4)')           
            call yaml_map('D(Na<->Na) [Bohr^8 Hartree]', parameters(4,1),  fmt='(1pe10.4)')           
            call yaml_map('D(Na<->Cl) [Bohr^8 Hartree]', parameters(4,2),  fmt='(1pe10.4)')           
            call yaml_map('D(Cl<->Cl) [Bohr^8 Hartree]', parameters(4,3),  fmt='(1pe10.4)')           
            call yaml_map('rho^-1 (Na<->Na)[Bohr^-1]',    parameters(2,1),   fmt='(1pe10.4)')           
            call yaml_map('rho^-1 (Na<->Cl)[Bohr^-1]',    parameters(2,2),   fmt='(1pe10.4)')           
            call yaml_map('rho^-1 (Cl<->Cl)[Bohr^-1]',    parameters(2,3),   fmt='(1pe10.4)')           
            call yaml_mapping_close() 
        case('NaF')
            call yaml_mapping_open('Using NaF parameters from'//&           
                 ' Adams, D. J., McDonald, I. R. J. Phys. C: Solid State Phys., Vol. 7, 1974')
            nameinteraction(1) = 'NaNa' 
            parameters(1,1) = 11.6391827398_gp !A(Na<->Na) 
            parameters(2,1) = 1.6035673097_gp !1/rho (Na<->Na)
            parameters(3,1) = 1.75485576_gp !C(Na<->Na) 
            parameters(4,1) = 2.9841448_gp !D(Na<->Na) 
            parameters(5,1) = 1.0_gp !q_Na*q_Na
            parameters(6,1) = 0.0_gp !???
            nameinteraction(2) = 'NaF' 
            parameters(1,2) = 9.5687865911_gp !A(Na<->F) 
            parameters(2,2) = 1.6035673097_gp  !1/rho (Na<->F)
            parameters(3,2) = 4.7005065_gp  !C(Na<->F) 
            parameters(4,2) = 14.1746878_gp !D(Na<->F) 
            parameters(5,2) = -1.0_gp !q_Na*q_F
            parameters(6,2) = 0.0_gp !???
            nameinteraction(3) = 'FF' 
            parameters(1,3) = 7.3750085331_gp !A(F<->F) 
            parameters(2,3) = 1.6035673097_gp  !1/rho (F<->F)
            parameters(3,3) = 17.2351905_gp !C(F<->F) 
            parameters(4,3) = 74.60362_gp !D(F<->F) 
            parameters(5,3) = 1.0_gp !q_F*q_F
            parameters(6,3) = 0.0_gp !???
            call yaml_map('A(Na<->Na) [Hartree]', parameters(1,1) ,  fmt='(1pe10.4)')           
            call yaml_map('A(Na<->F) [Hartree]', parameters(1,2),  fmt='(1pe10.4)')           
            call yaml_map('A(F<->F) [Hartree]', parameters(1,3),  fmt='(1pe10.4)')           
            call yaml_map('C(Na<->Na) [Bohr^6 Hartree]', parameters(3,1),  fmt='(1pe10.4)')           
            call yaml_map('C(Na<->F) [Bohr^6 Hartree]', parameters(3,2),  fmt='(1pe10.4)')           
            call yaml_map('C(F<->F) [Bohr^6 Hartree]', parameters(3,3),  fmt='(1pe10.4)')           
            call yaml_map('D(Na<->Na) [Bohr^8 Hartree]', parameters(4,1),  fmt='(1pe10.4)')           
            call yaml_map('D(Na<->F) [Bohr^8 Hartree]', parameters(4,2),  fmt='(1pe10.4)')           
            call yaml_map('D(F<->F) [Bohr^8 Hartree]', parameters(4,3),  fmt='(1pe10.4)')           
            call yaml_map('rho^-1 (Na<->Na)[Bohr^-1]',    parameters(2,1),   fmt='(1pe10.4)')           
            call yaml_map('rho^-1 (Na<->F)[Bohr^-1]',    parameters(2,2),   fmt='(1pe10.4)')           
            call yaml_map('rho^-1 (F<->F)[Bohr^-1]',    parameters(2,3),   fmt='(1pe10.4)')           
            call yaml_mapping_close() 
        case('default')
            initialized=.false.
            call f_err_throw('No "default" parameter set for BHMTF defined.')
        case default
            initialized=.false.
            call f_err_throw('Following parameter set for BHMTF force field '//&
                'is unknown: '//trim(paramset))
        end select
    endif

    itypinter=0
    do itypat=1,ntypat
        strtmpi=namat(itypat)
        leni=len_trim(strtmpi)
        do jtypat=itypat,ntypat
            itypinter=itypinter+1
            typinteraction(itypat,jtypat)=itypinter
            typinteraction(jtypat,itypat)=itypinter
            strtmpj=namat(jtypat)
            lenj=len_trim(strtmpj)
            namatnamat1=strtmpi(1:leni)//strtmpj(1:lenj)
            namatnamat2=strtmpj(1:lenj)//strtmpi(1:leni)
            do iall=1,nall
                if(namatnamat1==trim(adjustl(nameinteraction(iall))) .or. namatnamat2==trim(adjustl(nameinteraction(iall)))) then
                aaa(itypinter)=parameters(1,iall)
                bbb(itypinter)=parameters(2,iall)
                ccc(itypinter)=parameters(3,iall)
                ddd(itypinter)=parameters(4,iall)
                eee(itypinter)=parameters(5,iall)
                fff(itypinter)=parameters(6,iall)
                endif
            enddo
            !write(*,*) namatnamat1,namatnamat2
        enddo
    enddo
    initialized=.true.
end subroutine init_bhmtf

subroutine energyandforces_bhmtf(nat,rat,fat,epot)
! calculates the  Fumi-Tosi-Born-Meyer potential and forces in a spline representation
! THe array fps holds the spline coefficeints for the energy and fdsp for the forces
    use module_base
!    use shrtrngpotentialcoeff, &
!    only:aaa,bbb,ccc,ddd,eee,fff,fsp,fdsp,typinteraction,hsp
!    use atoms, only:typat
    implicit none
    integer::itypinter
    real(gp)::a2,a3,a4,a5,rinv2,rinv6,rinv8,rinv10,rinvqq,rinv28,rinv30
    integer::nat !number of atoms.
    real(gp)::rat(3,nat) !positions x,y,z of atoms.
    real(gp)::epot !short range electrostatic energy
    real(gp)::fat(3,nat) !forces exerted on atoms.
    !dummy variables
    real(gp)::dx,dy,dz,sclinv,r,rsq,xiat,yiat,ziat,alphainv
    real(gp)::t,tt,tt1,t2,tt3,ttt
    real(gp)::fx,fy,fz,pi,hspinv,rhspinv,rinv,spf,spfd
    integer::iat,jat,isp

    if(initialized==.false.)then
    endif

    pi=4.0_gp*atan(1.0_gp)
    hspinv=1.0_gp/hsp
    !write(*,*) 'inside shortenergy  hsp=',hsp
    epot=0.0_gp;fat=0.0_gp
    do iat=1,nat
        xiat=rat(1,iat);yiat=rat(2,iat);ziat=rat(3,iat)
        do jat=iat+1,nat
            dx=xiat-rat(1,jat)
            dy=yiat-rat(2,jat)
            dz=ziat-rat(3,jat)
            rsq= dx*dx+dy*dy+dz*dz
            r=sqrt(rsq)
            rinv=1.0_gp/r
            itypinter=typinteraction(typat(iat),typat(jat))
            !write(*,*) 'itypinter=',itypinter
            !-------------------------------------------------------------------
!            a2=ccc(itypinter);a3=ddd(itypinter)  !*eee(itypinter)
!            a4=aaa(itypinter);a5=bbb(itypinter)  !;a6=fff(itypinter)
!            !rinv2=rinv**2;rinv4=rinv2**2;rinv6=rinv4*rinv2;rinv7=rinv6*rinv
!            !rinv8=rinv6*rinv2;rinv9=rinv8*rinv
!            rinv2=rinv**2;rinv6=rinv2**3;rinv8=rinv6*rinv2;rinv10=rinv8*rinv2
!            rinv28=rinv10*rinv10*rinv8;rinv30=rinv28*rinv2
!            rinvqq=rinv*eee(itypinter)
!            t2=a4*exp(-a5*r)
!            epot=epot+rinvqq - a2*rinv6 -a3*rinv8 + t2 !+ 1.d8*rinv28
!            !epot=(((epot+ t2) -a3*rinv8) - a2*rinv6) + rinvqq
!            ttt=rinv2*rinvqq -6.0_gp*a2*rinv8 -8.0_gp*a3*rinv10 + a5*t2*rinv !+ 3.e9_gp*rinv30
!            fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
            !-------------------------------------------------------------------
            !for using splines.
            rhspinv=r*hspinv
            isp=floor(rhspinv)
            t=rhspinv-isp
            spf=fsp(0,isp,itypinter)+(fsp(1,isp,itypinter)+(fsp(2,isp,itypinter)+(fsp(3,isp,itypinter)+ &
                fsp(4,isp,itypinter)*t)*t)*t)*t
            spfd=fdsp(0,isp,itypinter)+(fdsp(1,isp,itypinter)+(fdsp(2,isp,itypinter)+fdsp(3,isp,itypinter)*t)*t)*t
            epot=epot+spf
            ttt=-spfd*rinv
            fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
            !-------------------------------------------------------------------
            fat(1,iat)=fat(1,iat)+fx
            fat(2,iat)=fat(2,iat)+fy
            fat(3,iat)=fat(3,iat)+fz
            fat(1,jat)=fat(1,jat)-fx
            fat(2,jat)=fat(2,jat)-fy
            fat(3,jat)=fat(3,jat)-fz
        enddo
    enddo
end subroutine energyandforces_bhmtf

!*******************************************************************************

subroutine preparespline(nsp,rcut)
!nsp=5*10**5
!rcut=sqrt(3)*130.0_gp
! nsp is the user determined value of the spline intervalls
! rcut the largest inter-particle distance possible within the simulations cell
    implicit none
    real(gp)::rcut,rco1,rco2,a2,a3,a4,a5,a6
    integer::nsp,itypinter
    hsp=rcut/real(nsp,8)
    write(*,*) 'nsp,hsp',nsp,hsp
    write(*,*) 'reasonable values of hsp are of the order of 1.e-3'
    rco2=rcut
    rco1=1.5_gp/0.529177_gp
    do itypinter=1,ntypinter
        a2=ccc(itypinter);a3=ddd(itypinter);a6=eee(itypinter)
        a4=aaa(itypinter);a5=bbb(itypinter)  !;a6=fff(itypinter)
        call spline_bm(fdsp(0,0,itypinter),fsp(0,0,itypinter),nsp,rco1,rco2,a2,a3,a4,a5,a6,hsp)
    enddo
end subroutine preparespline
!*******************************************************************************
real*16 function func(r,rco1,a2,a3,a4,a5,a6)
    implicit none
    real(16)::r,rco1,a2,a3,a4,a5,a6,alpha,beta,frco1,rinv,rinv2,rinv6,rinv8,rinv10
    if(r<rco1) then
        rinv=1.q0/rco1
        rinv2=rinv**2;rinv6=rinv2**3;rinv8=rinv6*rinv2;rinv10=rinv8*rinv2
        frco1=rinv*a6 - a2*rinv6 -a3*rinv8 + a4*exp(-a5*rco1)
        beta=10.q0*abs(frco1)
        alpha=(frco1-beta)/rco1
        func=alpha*r+beta
    else
        rinv=1.q0/r
        rinv2=rinv**2;rinv6=rinv2**3;rinv8=rinv6*rinv2;rinv10=rinv8*rinv2
        func=rinv*a6 - a2*rinv6 -a3*rinv8 + a4*exp(-a5*r)
    endif
end function func
!*******************************************************************************
real*16 function funcder(r,rco1,a2,a3,a4,a5,a6)
    implicit none
    real(16)::r,rco1,a2,a3,a4,a5,a6,alpha,beta,frco1,rinv,rinv2,rinv6,rinv8,rinv10
    if(r<rco1) then
        rinv=1.q0/rco1
        rinv2=rinv**2;rinv6=rinv2**3;rinv8=rinv6*rinv2;rinv10=rinv8*rinv2
        frco1=rinv*a6 - a2*rinv6 -a3*rinv8 + a4*exp(-a5*rco1)
        beta=10.q0*abs(frco1)
        alpha=(frco1-beta)/rco1
        funcder=alpha
    else
        rinv=1.q0/r
        rinv2=rinv**2;rinv6=rinv2**3;rinv8=rinv6*rinv2;rinv10=rinv8*rinv2
        funcder=-rinv2*a6 + 6.q0*a2*rinv6*rinv + 8.q0*a3*rinv8*rinv - a5*a4*exp(-a5*r)
    endif
end function funcder
!*******************************************************************************
real*16 function funcsecder(r,hsp,rco1,a2,a3,a4,a5,a6)
    implicit none
    real(16)::r,hsp,rco1,a2,a3,a4,a5,a6,rinv,rinv2,rinv6,rinv8,rinv10
    if(r<rco1) then
        funcsecder=0.q0
    else
        rinv=1.q0/r
        rinv2=rinv**2;rinv6=rinv2**3;rinv8=rinv6*rinv2;rinv10=rinv8*rinv2
        funcsecder=hsp*(2.q0*rinv2*a6*rinv - 42.q0*a2*rinv8 - 72.q0*a3*rinv10 + a5*a5*a4*exp(-a5*r))
    endif
end function funcsecder
!*******************************************************************************
!*******************************************************************************
subroutine spline_bm(fdsp,fsp,nsp,rco1,rco2,a2,a3,a4,a5,a6,hsp)
    use module_base
    !computes the spline coefficients for forces and energy,
    !         intermediate results are calculated in quadruple precision
    !nsp                ---->  number of spline nodes
    !fdsp(0:3,0:nsp) ---->  spline coefficients for short-range force
    !fsp(0:4,0:nsp)  ---->  spline coefficients for short range potential
    !hsp            ---->  length of interval between two spline nodes
    implicit none
    real(gp)::fdsp(0:3,0:nsp),fsp(0:4,0:nsp-1),rco1,rco2,a2,a3,a4,a5,a6,hsp
    integer::nsp,i,j
    real(16)::qdel,qfsp_int,qhsp,qrco1,qrco2,qa2,qa3,qa4,qa5,qa6,qr
    real(16), allocatable::qfdsp(:,:),qfsp(:,:)
    allocate(qfdsp(0:3,0:nsp),qfsp(0:4,0:nsp-1))
    qhsp=real(hsp,16);qa2=real(a2,16);qa3=real(a3,16);qa4=real(a4,16)
    qa5=real(a5,16);qa6=real(a6,16);qrco1=real(rco1,16);qrco2=real(rco2,16)
    qfdsp(0,0)=funcder(0.q0,qrco1,qa2,qa3,qa4,qa5,qa6)
    qfdsp(1,0)=0.q0
    !f and f' outside the origin.
    do i=1,nsp
        qr=qhsp*i
        qfdsp(0,i)=funcder(qr,qrco1,qa2,qa3,qa4,qa5,qa6)
        qfdsp(1,i)=funcsecder(qr,qhsp,qrco1,qa2,qa3,qa4,qa5,qa6)
    enddo
    !compute the 3-d and 4-th order coefficients of the local polynomials
    do i=0,nsp-1
        qfdsp(2,i)=3.q0*(-qfdsp(0,i)+qfdsp(0,i+1))-2.q0*qfdsp(1,i)-qfdsp(1,i+1)
        qfdsp(3,i)=2.q0*( qfdsp(0,i)-qfdsp(0,i+1))+     qfdsp(1,i)+qfdsp(1,i+1)
    enddo
    qfsp_int=1.q0 !for the time being assume that the potential starts from this value.
    do j=0,nsp-1
        !the zeroth coefficient is the value of the integral
        qfsp(0,j)=qfsp_int
        !the coefficients of the energy spline are obtained by integration
        qfsp(1,j)=qhsp*qfdsp(0,j)
        qfsp(2,j)=qhsp*qfdsp(1,j)*0.5q0
        qfsp(3,j)=qhsp*qfdsp(2,j)/3.q0
        qfsp(4,j)=qhsp*qfdsp(3,j)*0.25q0
        ! calculate the value of integral over the subinterval
        qfsp_int=qfsp_int+qfsp(1,j)+qfsp(2,j)+qfsp(3,j)+qfsp(4,j)
    enddo
    !correct the additive constant:
    !the value of the potential at rco2 should be correct
    qr=qrco2
    qdel=func(qr,qrco1,qa2,qa3,qa4,qa5,qa6)-qfsp_int
    do j=0,nsp-1
        qfsp(0,j)=qfsp(0,j)+qdel
    enddo
    fdsp(0:3,0:nsp)=real(qfdsp(0:3,0:nsp),8)
    fsp(0:4,0:nsp-1)=real(qfsp(0:4,0:nsp-1),8)
    deallocate(qfdsp,qfsp)
end subroutine spline_bm
end module module_BornHugginsMayerTosiFumi
