!> @file
!!  Contains a part of code which is included
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

  if(  iand( in%potshortcut, 16)  > 0 ) then
     b2BN=2
  else
     b2BN=1
  endif


  do b2Bcounter=1,b2BN

     if(b2Bcounter==1) then
        write(filename,'(A)' ) 'b2B_xanes'
        rhotarget=>rhopot
     else
        write(filename,'(A)') 'b2B_rho'
        rhotarget=>rhoXanes
     endif


     inquire(file=trim(trim(filename)//'.cube'),exist=exists)
     print *, "controllo ",  trim(filename)//'.cube', exists
     if(exists) then

        call read_density_cube_old(trim(filename), n1i_bB,n2i_bB,n3i_bB, 1 , hx_old ,hy_old ,hz_old , nat_b2B, rxyz_b2B, pot_bB )
        hx_old=hx_old*2
        hy_old=hy_old*2
        hz_old=hz_old*2


        if( (atoms%nat/nat_b2B)*nat_b2B /=  atoms%nat ) then
           if(iproc==0) write(*,*)  "   b2B_xanes cube  is not compatible with actual positions" 
           if(nproc>1) call MPI_Finalize(ierr)
           stop '      b2B_xanes cube  is not compatible with actual positions          '
        end if

     else

        if(b2BN>1) then
           if(iproc==0) write(*,*)  " b2B must be read only from *.cube when potential is energy dependent  " 
           if(nproc>1) call MPI_Finalize(ierr)
           stop '   b2B must be read only from *.cube when potential is energy dependent         '
        endif
        print  *, " reading atomic positions from file ","b2B_xanes.xyz"
        call read_atomic_file("b2B_xanes.xyz",iproc, atoms_b2B, rxyz_b2B )
        print *, "OK ", shape( rxyz_b2B )

        nat_b2B= (   Ubound(rxyz_b2B,2)  - Lbound(rxyz_b2B,2)  +  1 ) - ndebug


        if( (atoms%nat/nat_b2B)*nat_b2B /= atoms%nat ) then
           if(iproc==0) write(*,*)  "   b2B_xanes.xyz  is not compatible with actual positions" 
           if(nproc>1) call MPI_Finalize(ierr)
           stop '      b2B_xanes.xyz  is not compatible with actual positions          '
        end if

        print  *, " reading potential from file ","b2B_xanes.pot"
        call  read_potfile4b2B("b2B_xanes.pot",n1i_bB,n2i_bB,n3i_bB, pot_bB, alat1_bB, alat2_bB, alat3_bB)
        print  *, " reading OK "


        if( atoms_b2B%geocode/='F') then
           hx_old = 2*alat1_bB / (n1i_bB)
           hy_old = 2*alat2_bB / (n2i_bB)
           hz_old = 2*alat3_bB / (n3i_bB)
        else
           hx_old = 2*alat1_bB / (n1i_bB-2)
           hy_old = 2*alat2_bB / (n2i_bB-2)
           hz_old = 2*alat3_bB / (n3i_bB-2)
        endif

        call deallocate_atoms(atoms_b2B,subname) 

     endif


     allocate(rhopottmp( max(n1i_bB,n1i),max(n2i_bB,n2i),max(n3i_bB,n3d),in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rhopottmp,'rhopottmp',subname)


     rhotarget=0.0_gp



     itype=16
     nd=2**20

     allocate(  intfunc_x(0:nd+ndebug),stat=i_stat )
     call memocc(i_stat,intfunc_x,'intfunc_x',subname)
     allocate( intfunc_y(0:nd+ndebug) ,stat=i_stat )
     call memocc(i_stat,intfunc_y,'intfunc_y',subname)

     print *, " scaling function for interpolation "

     call scaling_function4b2B(itype,nd,nrange,intfunc_x,intfunc_y)  ! intervallo di 32 con 2**20 punti
     if( abs(intfunc_y(nd/2)-1)>1.0e-10 ) then
        stop " wrong scaling function 4b2B: not a centered one "
     endif

     i_all=-product(shape(intfunc_x))*kind(intfunc_x)
     deallocate(intfunc_x,stat=i_stat)
     call memocc(i_stat,i_all,'intfunc_x',subname)

     allocate(auxint(n1i+n2i+n3i+ndebug),stat=i_stat)
     call memocc(i_stat,auxint,'auxint',subname)

     Nreplicas = atoms%nat / nat_b2B
     dumvect3d(1)=atoms%alat1
     dumvect3d(2)=atoms%alat2
     dumvect3d(3)=atoms%alat3

     do ireplica=0, Nreplicas-1

        replicaoffset = (ireplica)*(   nat_b2B      )

        do j=1,3
           shift_b2B(j) = rxyz(j,1+replicaoffset) - rxyz_b2B(j,1) 

           do iat=2+ireplica*nat_b2B, (ireplica+1)*nat_b2B

              shiftdiff = shift_b2B(j) - (rxyz(j,iat) -rxyz_b2B(j,iat-replicaoffset)   )

              if( abs( shiftdiff )>1.0e-4 .and.  abs(abs(shiftdiff)-dumvect3d(j))  >1.0e-4) then
                 if(iproc==0) write(*,*)  "   b2B_xanes  positions are not compatible with actual positions" 
                 if(nproc>1) call MPI_Finalize(ierr)
                 stop '      b2B_xanes positions are not compatible with actual positions          '
              end if
           enddo
        enddo

        print *,"for replica ", ireplica,  "SHIFT " , shift_b2B



        rhopottmp=0.0_gp
        do iz_bB = 1,n3i_bB
           do iy_bB=1,n2i_bB
              do ix_bB=1,n1i_bB
                 rhopottmp(ix_bB,iy_bB,iz_bB +i3xcsh,1) =  pot_bB(ix_bB  + (iy_bB-1)*n1i_bB  + (iz_bB-1)*n1i_bB*n2i_bB)
              enddo
           enddo
        enddo

        i_all=-product(shape(pot_bB))*kind(pot_bB)
        deallocate(pot_bB,stat=i_stat)
        call memocc(i_stat,i_all,'rho',subname)



        do iz_bB = 1,n3i_bB-1
           do iy_bB=1,n2i_bB
              auxint = 0.0_gp
              do ix_bB=1,n1i_bB
                 rx_bB = hx_old*(ix_bB-1)           /2.0   +  shift_b2B(1)
                 minX_B  =  max(1,NINT((rx_bB -8*hx_old/2)/(hx/2.0)))
                 maxX_B  =  min(n1i,NINT((rx_bB +8*hx_old/2)/(hx/2.0)))

                 minX_B  =  NINT((rx_bB -8*hx_old/2)/(hx/2.0))
                 maxX_B  =  NINT((rx_bB +8*hx_old/2)/(hx/2.0))

                 do ixnl= minX_B , maxX_B 
                    ix = mod(ixnl-1 + n1i , n1i) +1

                    rx = hx*(ix-1  )/2.0  

                    shiftdiff = (rx-rx_bB)
                    if ( abs(shiftdiff -atoms%alat1) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat1
                    if ( abs(shiftdiff +atoms%alat1) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat1

                    idelta = NINT( shiftdiff *2**15/(hx_old/2))  
                    factx = intfunc_y(nd/2+idelta)
!!$                    print *, rx, rx_bB, ix, ix_bB      , factx
                    auxint(ix) = auxint(ix) + &
                         factx * rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
!!$                 print *, auxint(ix_bB) ,  rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
              enddo
              rhopottmp(:,iy_bB,iz_bB+i3xcsh,1)=auxint(1:n1i)
           enddo
        enddo

        do iz_bB = 1,n3i_bB-1
           do ix_bB=1,n1i
              auxint = 0.0_gp
              do iy_bB=1,n2i_bB
                 ry_bB = hy_old*(iy_bB-1)/2.0   +  shift_b2B(2)
                 minY_B  =  max(1  ,NINT((ry_bB -8*hy_old/2)/(hy/2.0)))
                 maxY_B  =  min(n2i,NINT((ry_bB +8*hy_old/2)/(hy/2.0)))

                 minY_B  =  NINT((ry_bB -8*hy_old/2)/(hy/2.0))
                 maxY_B  =  NINT((ry_bB +8*hy_old/2)/(hy/2.0))

                 do iynl= minY_B , maxY_B 
                    iy = mod(iynl-1 + n2i , n2i) +1

                    ry = hy*(iy-1  )/2.0  

                    shiftdiff = (ry-ry_bB)
                    if ( abs(shiftdiff -atoms%alat2) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat2
                    if ( abs(shiftdiff +atoms%alat2) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat2


                    idelta = NINT(shiftdiff *2**15/(hy_old/2))
                    facty = intfunc_y(nd/2+idelta)
                    auxint(iy) = auxint(iy) + &
                         facty * rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
              enddo
              rhopottmp(ix_bB ,:,iz_bB+i3xcsh,1)=auxint(1:n2i)
           enddo
        enddo

        do ix_bB=1,n1i
           do iy_bB=1,n2i
              auxint = 0.0_gp
              do iz_bB = 1,n3i_bB-1
                 rz_bB = hz_old*(iz_bB-1)           /2.0   +  shift_b2B(3)

                 minZ_B  =  max(1  ,  NINT((rz_bB -8*hz_old/2)/(hz/2.0)))
                 maxZ_B  =  min(n3i-i3xcsh , NINT((rz_bB +8*hz_old/2)/(hz/2.0)))

                 minZ_B  =    NINT((rz_bB -8*hz_old/2)/(hz/2.0))
                 maxZ_B  =    NINT((rz_bB +8*hz_old/2)/(hz/2.0))

                 do iznl= minZ_B , maxZ_B 

                    iz = mod(iznl-1 + n3i , n3i) +1

                    rz = hz*(iz-1  )/2.0  

                    shiftdiff = (rz-rz_bB)
                    if ( abs(shiftdiff -atoms%alat3) < abs(shiftdiff)) shiftdiff=shiftdiff -atoms%alat3
                    if ( abs(shiftdiff +atoms%alat3) < abs(shiftdiff)) shiftdiff=shiftdiff +atoms%alat3

                    idelta = NINT( shiftdiff *2**15/(hz_old/2.0))     
                    factz = intfunc_y(nd/2+idelta)
                    auxint(iz+i3xcsh) = auxint(iz+i3xcsh) + &
                         factz * rhopottmp(ix_bB,iy_bB,iz_bB+i3xcsh,1)
                 enddo
              enddo
              rhotarget(ix_bB ,iy_bB, : ,1)= rhotarget(ix_bB ,iy_bB, : ,1)+auxint(1:n3i)
           enddo
        enddo
     enddo
     i_all=-product(shape(rxyz_b2B))*kind(rxyz_b2B)
     deallocate(rxyz_b2B,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz',subname)
     i_all=-product(shape(rhopottmp))*kind(rhopottmp)
     deallocate(rhopottmp,stat=i_stat)
     call memocc(i_stat,i_all,'rhopottmp',subname)
     i_all=-product(shape(auxint))*kind(auxint)
     deallocate(auxint,stat=i_stat)
     call memocc(i_stat,i_all,'auxint',subname)
     i_all=-product(shape(intfunc_y))*kind(intfunc_y)
     deallocate(intfunc_y,stat=i_stat)
     call memocc(i_stat,i_all,'intfunc_y',subname)
  enddo


  if (iproc == 0) write(*,*) 'writing NEW local_potential.pot'

  call plot_density_old(atoms%geocode,'local_potentialb2BNEW.pot',iproc,nproc,&
       n1,n2,n3,n1i,n2i,n3i,n3p,&
       atoms%alat1,atoms%alat2,atoms%alat3,ngatherarr,rhopot(1,1,1+i3xcsh,1))

  print *," exiting b2B"
