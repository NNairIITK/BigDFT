  use libBigDFT

        implicit real*8 (a-h,o-z)

! atomic coordinates, forces
        real*8, allocatable, dimension(:,:) :: rxyz, fxyz, rxyz_old
        logical output_wf,output_grid
        character*20 tatonam
! atomic types
        integer, allocatable, dimension(:) :: iatype
        character*20 :: atomnames(100), units
        real*8, pointer :: psi(:,:), eval(:)
        integer, pointer :: keyv(:), keyg(:,:)
!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface
        include 'mpif.h'
        include 'parameters.h'

! Start MPI in parallel version
        if (parallel) then
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        write(6,*) 'mpi started',iproc,nproc
        call system("hostname")
        else
        nproc=1
        iproc=0
        endif

!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$       write(*,*) 'iproc,iam,npr',iproc,iam,npr
!$omp end parallel


! read atomic positions
        open(unit=9,file='posinp',status='old')
        read(9,*) nat,units
        if (iproc.eq.0) write(6,*) 'nat=',nat
        allocate(rxyz_old(3,nat),rxyz(3,nat),iatype(nat),fxyz(3,nat))
        ntypes=0
        do iat=1,nat
        read(9,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),tatonam
         do ityp=1,ntypes
           if (tatonam.eq.atomnames(ityp)) then
              iatype(iat)=ityp
              goto 200
           endif
         enddo
         ntypes=ntypes+1
         if (ntypes.gt.100) stop 'more than 100 atomnames not permitted'
         atomnames(ityp)=tatonam
         iatype(iat)=ntypes
200        continue
        if (units.eq.'angstroem') then
! if Angstroem convert to Bohr
        do i=1,3 ;  rxyz(i,iat)=rxyz(i,iat)/.529177d0  ; enddo
        else if  (units.eq.'atomic' .or. units.eq.'bohr') then
        else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
        endif
        enddo
        close(9)
        do ityp=1,ntypes
        if (iproc.eq.0) write(*,*) 'atoms of type ',ityp,' are ',atomnames(ityp)
        enddo

           ampl=1.d-1  ! amplitude for random displacement away from input file geometry (usually equilibrium geom.)
           if (iproc.eq.0) write(*,*) 'random displacemnt amplitude',ampl
           do iat=1,nat
              call random_number(tt)
              rxyz(1,iat)=rxyz(1,iat)+ampl*tt
              call random_number(tt)
              rxyz(2,iat)=rxyz(2,iat)+ampl*tt
              call random_number(tt)
              rxyz(3,iat)=rxyz(3,iat)+ampl*tt
           enddo
! geometry optimization
!        betax=2.d0   ! Cincodinine
!        betax=4.d0  ! Si H_4
        betax=7.5d0  ! silicon systems
!         betax=10.d0  !  Na_Cl clusters
        beta=.75d0*betax
        energyold=1.d100
       fluct=0.d0
       flucto=0.d0
       ngeostep=500
       do 500, igeostep=1,ngeostep

        output_grid=.false. 
        if (igeostep.eq.1) then 
            inputPsiId=0
        else
            inputPsiId=1
        endif
        output_wf=.true. 
        call cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz, &
                   & psi, keyg, keyv, nvctr_c, nvctr_f, nseg_c, nseg_f, norbp, norb, eval, &
                   & inputPsiId, output_grid, output_wf, n1, n2, n3, hgrid, rxyz_old)
        rxyz_old = rxyz

        if (iproc.eq.0) call wtposout(igeostep-1,nat,rxyz,atomnames,iatype)

        
        if (energy.gt.energyold) then 
           beta=.5d0*beta
        else
           beta=min(1.05d0*beta,betax)
        endif

           sumx=0.d0
           sumy=0.d0
           sumz=0.d0
           sum=0.d0
           do iat=1,nat
              rxyz(1,iat)=rxyz(1,iat)+beta*fxyz(1,iat)
              rxyz(2,iat)=rxyz(2,iat)+beta*fxyz(2,iat)
              rxyz(3,iat)=rxyz(3,iat)+beta*fxyz(3,iat)
              sum=sum+fxyz(1,iat)**2+fxyz(2,iat)**2+fxyz(3,iat)**2
              sumx=sumx+fxyz(1,iat)
              sumy=sumy+fxyz(2,iat)
              sumz=sumz+fxyz(3,iat)
              if (iproc.eq.0) write(*,'(a,i3,3(1x,e14.7))') 'fxyz ',iat,(fxyz(j,iat),j=1,3)
           end do
           fluctoo=flucto
           flucto=fluct
           fluct=sumx**2+sumy**2+sumz**2
        if (iproc.eq.0) then
           write(*,'(a,1x,e21.14,1x,e10.3)')'ANALYSIS OF FORCES energy, beta',energy,beta
           write(*,'(a,3(1x,e11.4))')'the norm of the forces is', sqrt(sum),sqrt(sum/nat),sqrt(sum/(3*nat))
           write(*,*) 'fluct',fluct
           write(*,*) 'stop comparison',sum,sqrt(1.d0*nat)*(fluct+flucto+fluctoo)/3.d0
           write(*,*)'the sum of the forces is'
           write(*,'(a16,3x,e16.8)')'x direction',sumx
           write(*,'(a16,3x,e16.8)')'y direction',sumy
           write(*,'(a16,3x,e16.8)')'z direction',sumz
        endif

        if (sum.lt.sqrt(1.d0*nat)*(fluct+flucto+fluctoo)/3.d0) then   ! assume that fluct increases as sqrt(nat)
        if (iproc.eq.0) then
           write(*,*) 'Final positions'
           do iat=1,nat
              write(*,'(3(1x,e14.7),2x,a20)') (rxyz(j,iat),j=1,3),atomnames(iatype(iat))
           enddo
           call wtposout(igeostep,nat,rxyz,atomnames,iatype)
        endif
           write(*,*) 'No better convergence possible'
           goto 501
        endif
        energyold=energy

500 continue
501 continue

!!  write all the wavefunctions into files
  call  writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
       nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
  write(*,*) iproc,' finished writing waves of relaxed geometry'

  deallocate(psi, eval, keyg, keyv)
  deallocate(rxyz,rxyz_old,iatype,fxyz)

        if (parallel) call MPI_FINALIZE(ierr)

	end


         subroutine wtposout(igeostep,nat,rxyz,atomnames,iatype)
         implicit real*8 (a-h,o-z)
         character*20 :: atomnames(100), filename
         character*3  fn
         dimension rxyz(3,nat),iatype(nat)

         write(fn,'(i3.3)') igeostep
         filename = 'posout_'//fn//'.ascii'
        open(unit=9,file=filename)
         xmax=0.d0 ; ymax=0.d0 ; zmax=0.d0
         do iat=1,nat
            xmax=max(rxyz(1,iat),xmax)
            ymax=max(rxyz(2,iat),ymax)
            zmax=max(rxyz(3,iat),zmax)
         enddo
         write(9,*) nat,' atomic ', igeostep
         write(9,*) xmax+5.d0, 0.d0, ymax+5.d0
         write(9,*) 0.d0, 0.d0, zmax+5.d0
         do iat=1,nat
            write(9,'(3(1x,e14.7),2x,a20)') (rxyz(j,iat),j=1,3),atomnames(iatype(iat))
         enddo

         return
         end

