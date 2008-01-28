subroutine createIonicPotential(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,&
     hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,eion)
  use Poisson_Solver
  implicit none
  include 'mpif.h'
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  real(kind=8), intent(in) :: hxh,hyh,hzh,elecfield
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: pot_ion
  real(kind=8), intent(out) :: eion
  !local variables
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,jat,i1,i2,i3,j1,j2,j3,ii1,ii2,ii3,ix,iy,iz,ierr,ityp,jtyp,nspin,i3end,i3start
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc
  real(kind=8) :: hgridh,pi,rholeaked,dist,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot
  real(kind=8) :: ehart,eexcu,vexcu
  real(kind=8), dimension(4) :: charges_mpi

  call timing(iproc,'CrtLocPot     ','ON')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------- Ionic Potential Creation'
  end if

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)
  eion=0.d0
  do iat=1,nat
     ityp=iatype(iat)
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     !    ion-ion interaction
     do jat=1,iat-1
        dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
        jtyp=iatype(jat)
        eion=eion+real(nelpsp(jtyp)*nelpsp(ityp),kind=8)/dist
     enddo
  end do
  if (iproc.eq.0) write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion

  !Creates charge density arising from the ionic PSP cores
  if (n3pi >0 ) then

     !conditions for periodicity in the three directions
     perx=(geocode /= 'F')
     pery=(geocode == 'P')
     perz=(geocode /= 'F')

     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)

     call razero(n1i*n2i*n3pi,pot_ion)

     do iat=1,nat
        ityp=iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        ix=nint(rx/hxh) 
        iy=nint(ry/hyh) 
        iz=nint(rz/hzh)

        rloc=psppar(0,0,ityp)
        charge=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc
        ii1=nint(cutoff/hxh)
        ii2=nint(cutoff/hyh)
        ii3=nint(cutoff/hzh)

        do i3=iz-ii3,iz+ii3
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=iy-ii2,iy+ii2
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=ix-ii1,ix+ii1
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3.ge.i3s .and. j3.le. i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)-xp*charge
                 else if (.not. goz ) then
                    rholeaked=rholeaked+xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if
  ! Check
  tt=0.d0
  do j3= 1,n3pi!i3start,i3end
     !j3=i3+15-i3s+1
     do i2= -nbl2,2*n2+1+nbr2
        do i1= -nbl1,2*n1+1+nbr1
           ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked**hxh*hyh*hzh

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt

  if (nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked
     call MPI_ALLREDUCE(charges_mpi(1),charges_mpi(3),2,MPI_double_precision,  &
          MPI_SUM,MPI_COMM_WORLD,ierr)
     tt_tot=charges_mpi(3)
     rholeaked_tot=charges_mpi(4)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc.eq.0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

  call timing(iproc,'CrtLocPot     ','OF')


  !here the value of the datacode must be kept fixed
  nspin=1
  call PSolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
       pot_ion,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.false.,nspin)

  call timing(iproc,'CrtLocPot     ','ON')

  !print *,'ehartree',ehart
  if (n3pi > 0) then
     do iat=1,nat
        ityp=iatype(iat)

        rx=rxyz(1,iat)
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        ix=nint(rx/hxh)
        iy=nint(ry/hyh)
        iz=nint(rz/hzh)

        ! determine number of local terms
        nloc=0
        do iloc=1,4
           if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
        enddo
        rloc=psppar(0,0,ityp)
        cutoff=10.d0*rloc
        ii1=nint(cutoff/hxh)
        ii2=nint(cutoff/hyh)
        ii3=nint(cutoff/hzh)
        
        if (nloc /= 0) then

           !this part should be changed with a modulo, in order to preserve the periodicity

           !calculate start and end of the distributed pot
           i3start=max(max(-14,iz-ii3),i3s-15)
           i3end=min(min(2*n3+16,iz+ii3),i3s-1+n3pi-1-14)

           do i3=iz-ii3,iz+ii3
              z=real(i3,kind=8)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                 do i2=iy-ii2,iy+ii2
                    y=real(i2,kind=8)*hyh-ry
                    call ind_positions(pery,i2,n2,j2,goy)
                    if (goy) then
                       do i1=ix-ii1,ix+ii1
                          x=real(i1,kind=8)*hxh-rx
                          call ind_positions(perx,i1,n1,j1,gox)
                          if (gox) then
                             r2=x**2+y**2+z**2
                             arg=r2/rloc**2
                             xp=exp(-.5d0*arg)
                             tt=psppar(0,nloc,ityp)
                             do iloc=nloc-1,1,-1
                                tt=arg*tt+psppar(0,iloc,ityp)
                             enddo
                             ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                             pot_ion(ind)=pot_ion(ind)+xp*tt
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do

!!$
!!$           do i3=i3start,i3end
!!$              z=real(i3,kind=8)*hzh-rz
!!$              j3=i3+15-i3s+1
!!$              do i2=max(-nbl2,iy-ii2),min(2*n2+1+nbr2,iy+ii2)
!!$                 y=real(i2,kind=8)*hyh-ry
!!$                 do i1=max(-nbl1,ix-ii1),min(2*n1+1+nbr1,ix+ii1)
!!$                    x=real(i1,kind=8)*hxh-rx
!!$                    r2=x**2+y**2+z**2
!!$                    arg=r2/rloc**2
!!$                    xp=exp(-.5d0*arg)
!!$                    tt=psppar(0,nloc,ityp)
!!$                    do iloc=nloc-1,1,-1
!!$                       tt=arg*tt+psppar(0,iloc,ityp)
!!$                    enddo
!!$                    ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
!!$                    pot_ion(ind)=pot_ion(ind)+xp*tt
!!$                 enddo
!!$              enddo
!!$           enddo

        end if

     enddo

  end if

  !use rhopot to calculate the potential from a constant electric field along x direction
  if (elecfield /= 0.d0) then
     !constant electric field allowed only for free BC
     if (geocode == 'F') then
     if (iproc.eq.0) write(*,'(1x,a)') &
          'The constant electric field is allowed only for Free BC'
     stop
     end if
     if (iproc.eq.0) write(*,'(1x,a,1pe10.2)') &
          'Adding constant electric field of intensity',elecfield,&
          'Ha*Bohr'

     if (n3pi > 0) then
        do i3=1,n3pi
           do i2= -14,2*n2+16
              do i1= -14,2*n1+16
                 ind=i1+15+(i2+14)*(2*n1+31)+(i3-1)*(2*n1+31)*(2*n2+31)
                 pot_ion(ind)=pot_ion(ind)+0.5d0*elecfield*hxh*real(i1-n1,kind=8)
              enddo
           enddo
        enddo
     end if

  end if

  call timing(iproc,'CrtLocPot     ','OF')



end subroutine createIonicPotential


subroutine createIonicPotential_old(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,&
     hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,eion)

  use Poisson_Solver
  
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,nat,ntypes,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  real(kind=8), intent(in) :: hxh,hyh,hzh,elecfield
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), intent(out) :: eion
  real(kind=8), dimension(*), intent(out) :: pot_ion
  !local variables
  real(kind=8) :: ehart,eexcu,vexcu
  integer :: nspin
  call timing(iproc,'CrtLocPot     ','ON')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------- Ionic Potential Creation'
  end if


  ! Precalculate ionic potential from PSP charge densities and local Gaussian terms
  call input_rho_ion(geocode,iproc,nproc,ntypes,nat,iatype,rxyz,psppar, &
       & nelpsp,n1,n2,n3,n3pi,i3s,hxh,hyh,hzh,pot_ion,eion)
  if (iproc.eq.0) write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion

  call timing(iproc,'CrtLocPot     ','OF')


  !here the value of the datacode must be kept fixed
  nspin=1
  call PSolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
       pot_ion,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.false.,nspin)

  call timing(iproc,'CrtLocPot     ','ON')


  !print *,'ehartree',ehart
  if (n3pi > 0) then
     call addlocgauspsp(iproc,ntypes,nat,iatype,rxyz,psppar,&
          n1,n2,n3,n3pi,i3s,hxh,hyh,hzh,pot_ion)
  end if

  !use rhopot to calculate the potential from a constant electric field along x direction
  if (elecfield /= 0.d0) then
     !constant electric field allowed only for free BC
     if (geocode == 'F') then
     if (iproc.eq.0) write(*,'(1x,a)') &
          'The constant electric field is allowed only for Free BC'
     stop
     end if
     if (iproc.eq.0) write(*,'(1x,a,1pe10.2)') &
          'Adding constant electric field of intensity',elecfield,&
          'Ha*Bohr'

     if (n3pi > 0) call pot_constantfield(iproc,n1,n2,n3,n3pi,pot_ion,hxh,elecfield)

  end if

  call timing(iproc,'CrtLocPot     ','OF')

end subroutine createIonicPotential_old

subroutine input_rho_ion(geocode,iproc,nproc,ntypes,nat,iatype,rxyz,psppar, &
     & nelpsp,n1,n2,n3,n3pi,i3s,hxh,hyh,hzh,rho,eion)
  !Creates charge density arising from the ionic PSP cores
  implicit none
  include 'mpif.h'
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3,n3pi,i3s
  real(kind=8), intent(in) :: hxh,hyh,hzh
  real(kind=8), intent(out) :: eion
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(inout) :: rho
  !local variables
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,jat,i1,i2,i3,j1,j2,j3,ii1,ii2,ii3,ix,iy,iz,ierr,ityp,jtyp
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,n1i,n2i,n3i
  real(kind=8) :: hgridh,pi,rholeaked,dist,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot
  real(kind=8), dimension(:), allocatable :: charges_mpi
  
  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)
  eion=0.d0
  do iat=1,nat
     ityp=iatype(iat)
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     !    ion-ion interaction
     do jat=1,iat-1
        dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
        jtyp=iatype(jat)
        eion=eion+real(nelpsp(jtyp)*nelpsp(ityp),kind=8)/dist
     enddo
  end do

  if (n3pi >0 ) then

     !conditions for periodicity in the three directions
     perx=(geocode /= 'F')
     pery=(geocode == 'P')
     perz=(geocode /= 'F')

     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)

     n1i=2*n1+2+nbl1+nbr1
     n2i=2*n2+2+nbl2+nbr2
     n3i=2*n3+2+nbl3+nbr3

     call razero(n1i*n2i*n3pi,rho)

     do iat=1,nat
        ityp=iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        ix=nint(rx/hxh) 
        iy=nint(ry/hyh) 
        iz=nint(rz/hzh)

        rloc=psppar(0,0,ityp)
        charge=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc
        ii1=nint(cutoff/hxh)
        ii2=nint(cutoff/hyh)
        ii3=nint(cutoff/hzh)

        do i3=iz-ii3,iz+ii3
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=iy-ii2,iy+ii2
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=ix-ii1,ix+ii1
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3.ge.i3s .and. j3.le. i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    rho(ind)=rho(ind)-xp*charge
                 else if (.not. goz ) then
                    rholeaked=rholeaked+xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if
  ! Check
  tt=0.d0
  do j3= 1,n3pi!i3start,i3end
     !j3=i3+15-i3s+1
     do i2= -nbl2,2*n2+1+nbr2
        do i1= -nbl1,2*n1+1+nbr1
           ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
           tt=tt+rho(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked**hxh*hyh*hzh

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt

  if (nproc > 1) then
     allocate(charges_mpi(4),stat=i_stat)
     call memocc(i_stat,product(shape(charges_mpi))*kind(charges_mpi),'charges_mpi','input_rho_ion')
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked
     call MPI_ALLREDUCE(charges_mpi(1),charges_mpi(3),2,MPI_double_precision,  &
          MPI_SUM,MPI_COMM_WORLD,ierr)
     tt_tot=charges_mpi(3)
     rholeaked_tot=charges_mpi(4)
     i_all=-product(shape(charges_mpi))*kind(charges_mpi)
     deallocate(charges_mpi,stat=i_stat)
     call memocc(i_stat,i_all,'charges_mpi','input_rho_ion')
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc.eq.0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

end subroutine input_rho_ion

!determine the index in which the potential must be inserted, following the BC
!determine also whether the index is inside or outside the box for free BC
subroutine ind_positions(periodic,i,n,j,go)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: i,n
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic) then
     go=.true.
     j=modulo(i,n+1)
  else
     j=i
     if (i >= -14 .and. i <= 2*n+16) then
        go=.true.
     else
        go=.false.
     end if
  end if

end subroutine ind_positions

subroutine ext_buffers(periodic,nl,nr)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(out) :: nl,nr

  if (periodic) then
     nl=0
     nr=0
  else
     nl=14
     nr=15
  end if
end subroutine ext_buffers

subroutine pot_constantfield(iproc,n1,n2,n3,n3pi,pot,hgridh,elecfield)
  !Creates charge density arising from the ionic PSP cores
  implicit none
  include 'mpif.h'
  integer, intent(in) :: iproc,n1,n2,n3,n3pi
  real(kind=8), intent(in) :: hgridh,elecfield
  real(kind=8), dimension(*), intent(inout) :: pot
  !local variables
  integer :: i1,i2,i3,ind
  

  do i3=1,n3pi
     do i2= -14,2*n2+16
        do i1= -14,2*n1+16
           ind=i1+15+(i2+14)*(2*n1+31)+(i3-1)*(2*n1+31)*(2*n2+31)
           pot(ind)=pot(ind)+0.5d0*elecfield*hgridh*real(i1-n1,kind=8)
        enddo
     enddo
  enddo

end subroutine pot_constantfield

subroutine addlocgauspsp(iproc,ntypes,nat,iatype,rxyz,psppar,&
     n1,n2,n3,n3pi,i3s,hxh,hyh,hzh,pot)
  ! Add local Gaussian terms of the PSP to pot, where pot is distributed
  implicit none
  integer, intent(in) :: ntypes,nat,n1,n2,n3,n3pi,iproc,i3s
  real(kind=8), intent(in) :: hxh,hyh,hzh
  integer, dimension(nat), intent(in) :: iatype
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(-14:2*n1+16,-14:2*n2+16,n3pi), intent(inout) :: pot
  !local variables
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,i1,i2,i3,ii1,ii2,ii3,ix,iy,iz,ityp,iloc,nloc,i3start,i3end,j1,j2,j3
  real(kind=8) :: hgridh,rloc,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
 
  do iat=1,nat
     ityp=iatype(iat)

     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     ix=nint(rx/hxh)
     iy=nint(ry/hyh)
     iz=nint(rz/hzh)

     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo
     rloc=psppar(0,0,ityp)
     cutoff=10.d0*rloc
     ii1=nint(cutoff/hxh)
     ii2=nint(cutoff/hyh)
     ii3=nint(cutoff/hzh)

     if (nloc /= 0) then

        !this part should be changed with a modulo, in order to preserve the periodicity

        !calculate start and end of the distributed pot
        i3start=max(max(-14,iz-ii3),i3s-15)
        i3end=min(min(2*n3+16,iz+ii3),i3s-1+n3pi-1-14)

        do i3=i3start,i3end
           j3=i3+15-i3s+1
           do i2=max(-14,iy-ii2),min(2*n2+16,iy+ii2)
              do i1=max(-14,ix-ii1),min(2*n1+16,ix+ii1)
                 x=real(i1,kind=8)*hxh-rx
                 y=real(i2,kind=8)*hyh-ry
                 z=real(i3,kind=8)*hzh-rz
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 tt=psppar(0,nloc,ityp)
                 do iloc=nloc-1,1,-1
                    tt=arg*tt+psppar(0,iloc,ityp)
                 enddo
                 pot(i1,i2,j3)=pot(i1,i2,j3)+xp*tt
              enddo
           enddo
        enddo

     end if

  enddo
 end subroutine addlocgauspsp

