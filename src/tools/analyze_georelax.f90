    implicit real*8 (a-h,o-z)
    parameter(npmax=1000,natx=1000,bohr=0.5291772108d0)
    character*4 fn4
    character*15 filename
    character*20 atmn,units,label
    dimension epos(npmax),grad(npmax),displ(npmax)
    dimension rxyz(3,natx,npmax)

    np=0
     do ip=1,npmax
        write(fn4,'(i4.4)') ip
        filename = 'posout_'//fn4//'.xyz'
        open(unit=9,file=filename,status='old',iostat=ierror)
     if (ierror == 0) then
        np=np+1
        read(9,*) nat,units,epos(np),label,grad(np)
        read(9,*)
        do iat=1,nat
          read(9,*) atmn,t1,t2,t3
          if (units=='angstroem' .or. units=='angstroemd0') then ! if Angstroem convert to Bohr
            write(*,*) 'converted'
              rxyz(1,iat,np)=t1/bohr
              rxyz(2,iat,np)=t2/bohr
              rxyz(3,iat,np)=t3/bohr
          else
              rxyz(1,iat,np)=t1
              rxyz(2,iat,np)=t2
              rxyz(3,iat,np)=t3
          endif
        enddo
        close(9)
     endif
     end do

     do ip=1,np-1
        tt=0.d0
        do iat=1,nat
        do l=1,3
        tt=tt+(rxyz(l,iat,ip)-rxyz(l,iat,ip+1))**2
        enddo
        enddo
     displ(ip)=tt
     enddo

     displ(np)=0.d0
     write(*,*) 'energy diff to last config,    gradient norm   , displacement from previous config'
     do ip=1,np
     write(*,*) epos(ip)-epos(np),grad(ip),displ(ip)
     enddo

     end
