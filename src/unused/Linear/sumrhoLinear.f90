!!!!subroutine sumrholinear_auxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, phi, at, nscatterarr)
!!!!!
!!!!use module_base
!!!!use module_types
!!!!use libxc_functionals
!!!!use module_interfaces, exceptThisOne => sumrholinear_auxiliary
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc
!!!!type(orbitals_data),intent(in):: orbs
!!!!type(locreg_descriptors),intent(in):: Glr
!!!!type(input_variables),intent(in):: input
!!!!type(linearParameters),intent(inout):: lin
!!!!real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
!!!!real(8),dimension(max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp)),intent(in):: phi
!!!!type(atoms_data),intent(in):: at
!!!!integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!!
!!!!! Local variables
!!!!integer:: iorb, jorb, korb, istat, indLarge, i1, i2, i3, ilr, jlr, ind
!!!!integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, iall, istr, istri, istrj
!!!!integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, iiorb, jjorb
!!!!integer:: ierr, jproc, is, ie, nreceives
!!!!integer:: nfast, nslow, nsameproc, m, i1d0, j1d0, indri0, indrj0, indLarge0
!!!!real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern, t1, t2, time
!!!!real(8),dimension(:,:),allocatable:: densKern
!!!!character(len=*),parameter:: subname='sumrhoForLocalizedBasis2'
!!!!integer,dimension(mpi_status_size):: stat
!!!!logical:: sendComplete, receiveComplete
!!!!
!!!!
!!!!if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating auxiliary arrays for charge density...'
!!!!
!!!!!lin%comsr%communComplete=.false.
!!!!!lin%comsr%computComplete=.false.
!!!!
!!!!
!!!!!!! Allocate the density kernel.
!!!!!!allocate(densKern(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
!!!!!!call memocc(istat, densKern, 'densKern', subname)
!!!!!!
!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!call cpu_time(t1)
!!!!!!! Calculate the density kernel.
!!!!!!call dgemm('n', 't', lin%lb%orbs%norb, lin%lb%orbs%norb, orbs%norb, 1.d0, coeff(1,1), lin%lb%orbs%norb, &
!!!!!!     coeff(1,1), lin%lb%orbs%norb, 0.d0, densKern(1,1), lin%lb%orbs%norb)
!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!call cpu_time(t2)
!!!!!!time=t2-t1
!!!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time
!!!!!!
!!!!!!
!!!!!!! Define some constant factors.
!!!!!!hxh=.5d0*input%hx
!!!!!!hyh=.5d0*input%hy
!!!!!!hzh=.5d0*input%hz
!!!!!!if(input%nspin==1) then
!!!!!!    factor=2.d0/(hxh*hyh*hzh)
!!!!!!else
!!!!!!    factor=1.d0/(hxh*hyh*hzh)
!!!!!!end if
!!!!!!
!!!!!!! Initialize rho.
!!!!!!if (libxc_functionals_isgga()) then
!!!!!!    call razero(nrho, rho)
!!!!!!else
!!!!!!    ! There is no mpi_allreduce, therefore directly initialize to
!!!!!!    ! 10^-20 and not 10^-20/nproc.
!!!!!!    rho=1.d-20
!!!!!!    !call tenminustwenty(nrho, rho, nproc)
!!!!!!end if
!!!!
!!!!
!!!!! Check whether the communication has completed.
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t1)
!!!!nfast=0
!!!!nsameproc=0
!!!!testLoop: do
!!!!    do jproc=0,nproc-1
!!!!        do korb=1,lin%comsr%noverlaps(jproc)
!!!!            if(lin%comsr%communComplete(korb,jproc)) cycle
!!!!            call mpi_test(lin%comsr%comarr(8,korb,jproc), sendComplete, stat, ierr)
!!!!            call mpi_test(lin%comsr%comarr(9,korb,jproc), receiveComplete, stat, ierr)
!!!!            ! Attention: mpi_test is a local function.
!!!!            if(sendComplete .and. receiveComplete) lin%comsr%communComplete(korb,jproc)=.true.
!!!!        end do
!!!!    end do
!!!!    ! If we made it until here, either all all the communication is
!!!!    ! complete or we better wait for each single orbital.
!!!!    exit testLoop
!!!!end do testLoop
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t2)
!!!!time=t2-t1
!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for test:',time
!!!!
!!!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t1)
!!!!call mpiallred(lin%comsr%communComplete(1,0), nproc*maxval(lin%comsr%noverlaps), mpi_land, mpi_comm_world, ierr)
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t2)
!!!!time=t2-t1
!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for allreduce:',time
!!!!
!!!!
!!!!
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t1)
!!!!! Wait for the communications that have not completed yet
!!!!nslow=0
!!!!do jproc=0,nproc-1
!!!!    do korb=1,lin%comsr%noverlaps(jproc)
!!!!        if(lin%comsr%communComplete(korb,jproc)) then
!!!!            mpisource=lin%comsr%comarr(1,korb,jproc)
!!!!            mpidest=lin%comsr%comarr(5,korb,jproc)
!!!!            if(mpisource==mpidest) then
!!!!                nsameproc=nsameproc+1
!!!!            else
!!!!                nfast=nfast+1
!!!!            end if
!!!!            cycle
!!!!        end if
!!!!        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
!!!!        nslow=nslow+1
!!!!        call mpi_wait(lin%comsr%comarr(8,korb,jproc), stat, ierr)
!!!!        call mpi_wait(lin%comsr%comarr(9,korb,jproc), stat, ierr)
!!!!        lin%comsr%communComplete(korb,jproc)=.true.
!!!!        lin%comsr%computComplete(korb,jproc)=.true.
!!!!    end do
!!!!end do
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t2)
!!!!time=t2-t1
!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for wait:',time
!!!!!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!!                       nfast, ' could be overlapped with computation.'
!!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!!
!!!!
!!!!do iorb=1,lin%comsr%noverlaps(iproc)
!!!!    if(.not. lin%comsr%communComplete(iorb,iproc)) then
!!!!        write(*,'(a,i0,a,i0,a)') 'ERROR: communication of orbital ', iorb, ' to process ', iproc, ' failed!'
!!!!        stop
!!!!    end if
!!!!    !!if(.not. lin%comsr%computComplete(iorb,iproc)) then
!!!!    !!    write(*,'(a,i0,a,i0,a)') 'ERROR: computation of orbital ', iorb, ' on process ', iproc, ' failed!'
!!!!    !!    stop
!!!!    !!end if
!!!!end do
!!!!
!!!!
!!!!
!!!!! Now calculate the charge density. Each process calculates only one slice of the total charge density.
!!!!! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
!!!!! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
!!!!! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t1)
!!!!
!!!!
!!!!! First determine the bounds of the auxiliary array
!!!!
!!!!! Bounds of the slice in global coordinates.
!!!!is=nscatterarr(iproc,3)-14
!!!!ie=is+nscatterarr(iproc,1)-1
!!!!
!!!!
!!!!
!!!!!allocate(comsr%auxarray(comsr%maxsize_auxarray,lin%comsr%noverlaps(iproc),lin%comsr%noverlaps(iproc), stat=istat)
!!!!!call memocc(istat, comsr%auxarray, 'comsr%auxarray', subname)
!!!!
!!!!
!!!!
!!!!
!!!!! Bounds of the slice in global coordinates.
!!!!is=nscatterarr(iproc,3)-14
!!!!ie=is+nscatterarr(iproc,1)-1
!!!!
!!!!totalCharge=0.d0
!!!!ind=1
!!!!do iorb=1,lin%comsr%noverlaps(iproc)
!!!!    iiorb=lin%comsr%overlaps(iorb) !global index of orbital iorb
!!!!    ilr=lin%comsr%comarr(4,iorb,iproc) !localization region of orbital iorb
!!!!    istri=lin%comsr%comarr(6,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
!!!!    !do jorb=1,lin%comsr%noverlaps(iproc)
!!!!    do jorb=iorb,lin%comsr%noverlaps(iproc)
!!!!        jjorb=lin%comsr%overlaps(jorb) !global indes of orbital jorb
!!!!        jlr=lin%comsr%comarr(4,jorb,iproc) !localization region of orbital jorb
!!!!        istrj=lin%comsr%comarr(6,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
!!!!        ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
!!!!        i1s=max(2*lin%lzd%llr(ilr)%ns1-14,2*lin%lzd%llr(jlr)%ns1-14)
!!!!        i1e=min(2*lin%lzd%llr(ilr)%ns1-14+lin%lzd%llr(ilr)%d%n1i-1,2*lin%lzd%llr(jlr)%ns1-14+lin%lzd%llr(jlr)%d%n1i-1)
!!!!        i2s=max(2*lin%lzd%llr(ilr)%ns2-14,2*lin%lzd%llr(jlr)%ns2-14)
!!!!        i2e=min(2*lin%lzd%llr(ilr)%ns2-14+lin%lzd%llr(ilr)%d%n2i-1,2*lin%lzd%llr(jlr)%ns2-14+lin%lzd%llr(jlr)%d%n2i-1)
!!!!        i3s=max(2*lin%lzd%llr(ilr)%ns3-14,2*lin%lzd%llr(jlr)%ns3-14,is)
!!!!        i3e=min(2*lin%lzd%llr(ilr)%ns3-14+lin%lzd%llr(ilr)%d%n3i-1,2*lin%lzd%llr(jlr)%ns3-14+lin%lzd%llr(jlr)%d%n3i-1,ie)
!!!!        !factorTimesDensKern = factor*densKern(iiorb,jjorb)
!!!!        ! Now loop over all points in the box in which the orbitals overlap.
!!!!        do i3=i3s,i3e !bounds in z direction
!!!!            !!i3d=i3-i3s+1 !z coordinate of orbital iorb with respect to the overlap box
!!!!            !!j3d=i3-i3s+1 !z coordinate of orbital jorb with respect to the overlap box
!!!!            i3d=i3-max(is,2*lin%lzd%llr(ilr)%ns3-14)+1 !z coordinate of orbital iorb with respect to the overlap box
!!!!            j3d=i3-max(is,2*lin%lzd%llr(jlr)%ns3-14)+1 !z coordinate of orbital jorb with respect to the overlap box
!!!!            indi3=(i3d-1)*lin%lzd%llr(ilr)%d%n2i*lin%lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
!!!!            indj3=(j3d-1)*lin%lzd%llr(jlr)%d%n2i*lin%lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
!!!!            indl3=(i3-is)*Glr%d%n2i*Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!!!            do i2=i2s,i2e !bounds in y direction
!!!!                i2d=i2-2*lin%lzd%llr(ilr)%ns2 !y coordinate of orbital iorb with respect to the overlap box
!!!!                j2d=i2-2*lin%lzd%llr(jlr)%ns2 !y coordinate of orbital jorb with respect to the overlap box
!!!!                indi2=(i2d+15-1)*lin%lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
!!!!                indj2=(j2d+15-1)*lin%lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
!!!!                indl2=(i2+15-1)*Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!!!                !!!! This is the old version.
!!!!                !!do i1=i1s,i1e
!!!!                !!    i1d=i1-2*lin%lzd%llr(ilr)%ns1
!!!!                !!    j1d=i1-2*lin%lzd%llr(jlr)%ns1
!!!!                !!    ! Now calculate the index in the boxes.
!!!!                !!    indri = indi3 + indi2 + i1d+15 + istri
!!!!                !!    indrj = indj3 + indj2 + j1d+15 + istrj
!!!!                !!    indLarge = indl3 + indl2 + i1+15
!!!!                !!    tt = factor*densKern(iiorb,jjorb)*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
!!!!                !!    rho(indLarge) = rho(indLarge) + tt
!!!!                !!    totalCharge = totalCharge + tt
!!!!                !!end do
!!!!                ! #####################################################################
!!!!                ! This is the new version.
!!!!                m=mod(i1e-i1s+1,4)
!!!!                if(m/=0) then
!!!!                    ! The following five variables hold some intermediate results to speed up the code.
!!!!                    i1d0=-2*lin%lzd%llr(ilr)%ns1 
!!!!                    j1d0=-2*lin%lzd%llr(jlr)%ns1
!!!!                    indri0 = indi3 + indi2 + 15 + istri
!!!!                    indrj0 = indj3 + indj2 + 15 + istrj
!!!!                    indLarge0 = indl3 + indl2 + 15
!!!!                    do i1=i1s,i1s+m-1
!!!!                        i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
!!!!                        j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
!!!!                        indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
!!!!                        indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
!!!!                        indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
!!!!                        !tt = factorTimesDensKern*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
!!!!                        lin%comsr%auxarray(ind)=lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
!!!!                        !rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
!!!!                        !totalCharge = totalCharge + tt !add the contribution to the total charge
!!!!                        ind=ind+1
!!!!                        !write(5000+iproc,*) indri, indrj
!!!!                    end do
!!!!                end if
!!!!                ! This is the same again, this time with unrolled loops.
!!!!                if(i1e-i1s+1>4) then
!!!!                    i1d0=-2*lin%lzd%llr(ilr)%ns1
!!!!                    j1d0=-2*lin%lzd%llr(jlr)%ns1
!!!!                    indri0 = indi3 + indi2 + 15 + istri
!!!!                    indrj0 = indj3 + indj2 + 15 + istrj
!!!!                    indLarge0 = indl3 + indl2 + 15
!!!!                    do i1=i1s+m,i1e,4
!!!!                        i1d=i1d0+i1
!!!!                        j1d=j1d0+i1
!!!!                        indri = indri0 + i1d
!!!!                        indrj = indrj0 + j1d
!!!!                        indLarge = indLarge0 + i1
!!!!                        !tt0 = factorTimesDensKern*lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
!!!!                        !tt1 = factorTimesDensKern*lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
!!!!                        !tt2 = factorTimesDensKern*lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
!!!!                        !tt3 = factorTimesDensKern*lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
!!!!                        lin%comsr%auxarray(ind  )=lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
!!!!                        lin%comsr%auxarray(ind+1)=lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
!!!!                        lin%comsr%auxarray(ind+2)=lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
!!!!                        lin%comsr%auxarray(ind+3)=lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
!!!!                        !rho(indLarge  ) = rho(indLarge  ) + tt0
!!!!                        !rho(indLarge+1) = rho(indLarge+1) + tt1
!!!!                        !rho(indLarge+2) = rho(indLarge+2) + tt2
!!!!                        !rho(indLarge+3) = rho(indLarge+3) + tt3
!!!!                        !totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
!!!!                        ind=ind+4
!!!!                        !write(5000+iproc,*) indri, indrj
!!!!                        !write(5000+iproc,*) indri+1, indrj+1
!!!!                        !write(5000+iproc,*) indri+2, indrj+2
!!!!                        !write(5000+iproc,*) indri+3, indrj+3
!!!!                    end do
!!!!                end if
!!!!            end do
!!!!        end do
!!!!    end do
!!!!end do
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t2)
!!!!time=t2-t1
!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for large loop:',time
!!!!
!!!!!!call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!if(iproc==0) write(*,'(1x,a,es20.12)') 'done. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh
!!!!!!
!!!!!!iall=-product(shape(densKern))*kind(densKern)
!!!!!!deallocate(densKern, stat=istat)
!!!!!!call memocc(istat, iall, 'densKern', subname)
!!!!
!!!!
!!!!end subroutine sumrholinear_auxiliary
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!subroutine sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, nrho, rho, at, nscatterarr)
!!!!!
!!!!use module_base
!!!!use module_types
!!!!use libxc_functionals
!!!!use module_interfaces, exceptThisOne => sumrholinear_withauxiliary
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, nrho
!!!!type(orbitals_data),intent(in):: orbs
!!!!type(locreg_descriptors),intent(in):: Glr
!!!!type(input_variables),intent(in):: input
!!!!type(linearParameters),intent(inout):: lin
!!!!real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
!!!!real(8),dimension(nrho),intent(out),target:: rho
!!!!type(atoms_data),intent(in):: at
!!!!integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!!
!!!!! Local variables
!!!!integer:: iorb, jorb, korb, istat, indLarge, i1, i2, i3, ilr, jlr
!!!!integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, iall, istr, istri, istrj
!!!!integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, iiorb, jjorb
!!!!integer:: ierr, jproc, is, ie, nreceives, ind
!!!!integer:: nfast, nslow, nsameproc, m, i1d0, j1d0, indri0, indrj0, indLarge0
!!!!real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern, t1, t2, time
!!!!real(8),dimension(:,:),allocatable:: densKern
!!!!character(len=*),parameter:: subname='sumrholinear_withauxiliary'
!!!!integer,dimension(mpi_status_size):: stat
!!!!logical:: sendComplete, receiveComplete
!!!!
!!!!
!!!!if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating charge density...'
!!!!
!!!!!lin%comsr%communComplete=.false.
!!!!!lin%comsr%computComplete=.false.
!!!!
!!!!
!!!!! Allocate the density kernel.
!!!!allocate(densKern(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
!!!!call memocc(istat, densKern, 'densKern', subname)
!!!!
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t1)
!!!!! Calculate the density kernel.
!!!!call dgemm('n', 't', lin%lb%orbs%norb, lin%lb%orbs%norb, orbs%norb, 1.d0, coeff(1,1), lin%lb%orbs%norb, &
!!!!     coeff(1,1), lin%lb%orbs%norb, 0.d0, densKern(1,1), lin%lb%orbs%norb)
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t2)
!!!!time=t2-t1
!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time
!!!!
!!!!
!!!!! Define some constant factors.
!!!!hxh=.5d0*input%hx
!!!!hyh=.5d0*input%hy
!!!!hzh=.5d0*input%hz
!!!!if(input%nspin==1) then
!!!!    factor=2.d0/(hxh*hyh*hzh)
!!!!else
!!!!    factor=1.d0/(hxh*hyh*hzh)
!!!!end if
!!!!
!!!!! Initialize rho.
!!!!if (libxc_functionals_isgga()) then
!!!!    call razero(nrho, rho)
!!!!else
!!!!    ! There is no mpi_allreduce, therefore directly initialize to
!!!!    ! 10^-20 and not 10^-20/nproc.
!!!!    rho=1.d-20
!!!!    !call tenminustwenty(nrho, rho, nproc)
!!!!end if
!!!!
!!!!
!!!!!!!! Check whether the communication has completed.
!!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!!call cpu_time(t1)
!!!!!!!nfast=0
!!!!!!!nsameproc=0
!!!!!!!testLoop: do
!!!!!!!    do jproc=0,nproc-1
!!!!!!!        do korb=1,lin%comsr%noverlaps(jproc)
!!!!!!!            if(lin%comsr%communComplete(korb,jproc)) cycle
!!!!!!!            call mpi_test(lin%comsr%comarr(8,korb,jproc), sendComplete, stat, ierr)
!!!!!!!            call mpi_test(lin%comsr%comarr(9,korb,jproc), receiveComplete, stat, ierr)
!!!!!!!            ! Attention: mpi_test is a local function.
!!!!!!!            if(sendComplete .and. receiveComplete) lin%comsr%communComplete(korb,jproc)=.true.
!!!!!!!        end do
!!!!!!!    end do
!!!!!!!    ! If we made it until here, either all all the communication is
!!!!!!!    ! complete or we better wait for each single orbital.
!!!!!!!    exit testLoop
!!!!!!!end do testLoop
!!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!!call cpu_time(t2)
!!!!!!!time=t2-t1
!!!!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for test:',time
!!!!!!!
!!!!!!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!!call cpu_time(t1)
!!!!!!!call mpiallred(lin%comsr%communComplete(1,0), nproc*maxval(lin%comsr%noverlaps), mpi_land, mpi_comm_world, ierr)
!!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!!call cpu_time(t2)
!!!!!!!time=t2-t1
!!!!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for allreduce:',time
!!!!!!!
!!!!!!!
!!!!!!!
!!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!!call cpu_time(t1)
!!!!!!!! Wait for the communications that have not completed yet
!!!!!!!nslow=0
!!!!!!!do jproc=0,nproc-1
!!!!!!!    do korb=1,lin%comsr%noverlaps(jproc)
!!!!!!!        if(lin%comsr%communComplete(korb,jproc)) then
!!!!!!!            mpisource=lin%comsr%comarr(1,korb,jproc)
!!!!!!!            mpidest=lin%comsr%comarr(5,korb,jproc)
!!!!!!!            if(mpisource==mpidest) then
!!!!!!!                nsameproc=nsameproc+1
!!!!!!!            else
!!!!!!!                nfast=nfast+1
!!!!!!!            end if
!!!!!!!            cycle
!!!!!!!        end if
!!!!!!!        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
!!!!!!!        nslow=nslow+1
!!!!!!!        call mpi_wait(lin%comsr%comarr(8,korb,jproc), stat, ierr)
!!!!!!!        call mpi_wait(lin%comsr%comarr(9,korb,jproc), stat, ierr)
!!!!!!!        lin%comsr%communComplete(korb,jproc)=.true.
!!!!!!!        lin%comsr%computComplete(korb,jproc)=.true.
!!!!!!!    end do
!!!!!!!end do
!!!!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!!!!call cpu_time(t2)
!!!!!!!time=t2-t1
!!!!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for wait:',time
!!!!!!!!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!!!!!                       nfast, ' could be overlapped with computation.'
!!!!!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!!!!!
!!!!!!!
!!!!!!!do iorb=1,lin%comsr%noverlaps(iproc)
!!!!!!!    if(.not. lin%comsr%communComplete(iorb,iproc)) then
!!!!!!!        write(*,'(a,i0,a,i0,a)') 'ERROR: communication of orbital ', iorb, ' to process ', iproc, ' failed!'
!!!!!!!        stop
!!!!!!!    end if
!!!!!!!    !!if(.not. lin%comsr%computComplete(iorb,iproc)) then
!!!!!!!    !!    write(*,'(a,i0,a,i0,a)') 'ERROR: computation of orbital ', iorb, ' on process ', iproc, ' failed!'
!!!!!!!    !!    stop
!!!!!!!    !!end if
!!!!!!!end do
!!!!
!!!!
!!!!
!!!!! Now calculate the charge density. Each process calculates only one slice of the total charge density.
!!!!! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
!!!!! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
!!!!! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t1)
!!!!
!!!!
!!!!! Bounds of the slice in global coordinates.
!!!!is=nscatterarr(iproc,3)-14
!!!!ie=is+nscatterarr(iproc,1)-1
!!!!
!!!!totalCharge=0.d0
!!!!do iorb=1,lin%comsr%noverlaps(iproc)
!!!!    iiorb=lin%comsr%overlaps(iorb) !global index of orbital iorb
!!!!    ilr=lin%comsr%comarr(4,iorb,iproc) !localization region of orbital iorb
!!!!    !istri=lin%comsr%comarr(6,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
!!!!    do jorb=1,lin%comsr%noverlaps(iproc)
!!!!        jjorb=lin%comsr%overlaps(jorb) !global indes of orbital jorb
!!!!        jlr=lin%comsr%comarr(4,jorb,iproc) !localization region of orbital jorb
!!!!        !istrj=lin%comsr%comarr(6,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
!!!!        ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
!!!!        i1s=max(2*lin%lzd%llr(ilr)%ns1-14,2*lin%lzd%llr(jlr)%ns1-14)
!!!!        i1e=min(2*lin%lzd%llr(ilr)%ns1-14+lin%lzd%llr(ilr)%d%n1i-1,2*lin%lzd%llr(jlr)%ns1-14+lin%lzd%llr(jlr)%d%n1i-1)
!!!!        i2s=max(2*lin%lzd%llr(ilr)%ns2-14,2*lin%lzd%llr(jlr)%ns2-14)
!!!!        i2e=min(2*lin%lzd%llr(ilr)%ns2-14+lin%lzd%llr(ilr)%d%n2i-1,2*lin%lzd%llr(jlr)%ns2-14+lin%lzd%llr(jlr)%d%n2i-1)
!!!!        i3s=max(2*lin%lzd%llr(ilr)%ns3-14,2*lin%lzd%llr(jlr)%ns3-14,is)
!!!!        i3e=min(2*lin%lzd%llr(ilr)%ns3-14+lin%lzd%llr(ilr)%d%n3i-1,2*lin%lzd%llr(jlr)%ns3-14+lin%lzd%llr(jlr)%d%n3i-1,ie)
!!!!        factorTimesDensKern = factor*densKern(iiorb,jjorb)
!!!!        ! Now loop over all points in the box in which the orbitals overlap.
!!!!        if(jorb>=iorb) then
!!!!            ind=lin%comsr%startingindex(jorb,iorb)
!!!!        else
!!!!            ind=lin%comsr%startingindex(iorb,jorb)
!!!!        end if
!!!!        do i3=i3s,i3e !bounds in z direction
!!!!            !!i3d=i3-i3s+1 !z coordinate of orbital iorb with respect to the overlap box
!!!!            !!j3d=i3-i3s+1 !z coordinate of orbital jorb with respect to the overlap box
!!!!            !i3d=i3-max(is,2*lin%lzd%llr(ilr)%ns3-14)+1 !z coordinate of orbital iorb with respect to the overlap box
!!!!            !j3d=i3-max(is,2*lin%lzd%llr(jlr)%ns3-14)+1 !z coordinate of orbital jorb with respect to the overlap box
!!!!            !indi3=(i3d-1)*lin%lzd%llr(ilr)%d%n2i*lin%lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
!!!!            !indj3=(j3d-1)*lin%lzd%llr(jlr)%d%n2i*lin%lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
!!!!            indl3=(i3-is)*Glr%d%n2i*Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!!!            do i2=i2s,i2e !bounds in y direction
!!!!                !i2d=i2-2*lin%lzd%llr(ilr)%ns2 !y coordinate of orbital iorb with respect to the overlap box
!!!!                !j2d=i2-2*lin%lzd%llr(jlr)%ns2 !y coordinate of orbital jorb with respect to the overlap box
!!!!                !indi2=(i2d+15-1)*lin%lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
!!!!                !indj2=(j2d+15-1)*lin%lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
!!!!                indl2=(i2+15-1)*Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!!!                !!!! This is the old version.
!!!!                !!do i1=i1s,i1e
!!!!                !!    i1d=i1-2*lin%lzd%llr(ilr)%ns1
!!!!                !!    j1d=i1-2*lin%lzd%llr(jlr)%ns1
!!!!                !!    ! Now calculate the index in the boxes.
!!!!                !!    indri = indi3 + indi2 + i1d+15 + istri
!!!!                !!    indrj = indj3 + indj2 + j1d+15 + istrj
!!!!                !!    indLarge = indl3 + indl2 + i1+15
!!!!                !!    tt = factor*densKern(iiorb,jjorb)*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
!!!!                !!    rho(indLarge) = rho(indLarge) + tt
!!!!                !!    totalCharge = totalCharge + tt
!!!!                !!end do
!!!!                ! #####################################################################
!!!!                ! This is the new version.
!!!!                m=mod(i1e-i1s+1,4)
!!!!                if(m/=0) then
!!!!                    ! The following five variables hold some intermediate results to speed up the code.
!!!!                    !i1d0=-2*lin%lzd%llr(ilr)%ns1 
!!!!                    !j1d0=-2*lin%lzd%llr(jlr)%ns1
!!!!                    !indri0 = indi3 + indi2 + 15 + istri
!!!!                    !indrj0 = indj3 + indj2 + 15 + istrj
!!!!                    indLarge0 = indl3 + indl2 + 15
!!!!                    do i1=i1s,i1s+m-1
!!!!                        !i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
!!!!                        !j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
!!!!                        !indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
!!!!                        !indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
!!!!                        indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
!!!!                        !tt = factorTimesDensKern*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
!!!!                        tt = factorTimesDensKern*lin%comsr%auxarray(ind)
!!!!                        rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
!!!!                        totalCharge = totalCharge + tt !add the contribution to the total charge
!!!!                        ind=ind+1
!!!!                    end do
!!!!                end if
!!!!                ! This is the same again, this time with unrolled loops.
!!!!                if(i1e-i1s+1>4) then
!!!!                    !i1d0=-2*lin%lzd%llr(ilr)%ns1
!!!!                    !j1d0=-2*lin%lzd%llr(jlr)%ns1
!!!!                    !indri0 = indi3 + indi2 + 15 + istri
!!!!                    !indrj0 = indj3 + indj2 + 15 + istrj
!!!!                    indLarge0 = indl3 + indl2 + 15
!!!!                    do i1=i1s+m,i1e,4
!!!!                        !i1d=i1d0+i1
!!!!                        !j1d=j1d0+i1
!!!!                        !indri = indri0 + i1d
!!!!                        !indrj = indrj0 + j1d
!!!!                        indLarge = indLarge0 + i1
!!!!                        !!tt0 = factorTimesDensKern*lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
!!!!                        !!tt1 = factorTimesDensKern*lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
!!!!                        !!tt2 = factorTimesDensKern*lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
!!!!                        !!tt3 = factorTimesDensKern*lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
!!!!                        tt0 = factorTimesDensKern*lin%comsr%auxarray(ind  )
!!!!                        tt1 = factorTimesDensKern*lin%comsr%auxarray(ind+1)
!!!!                        tt2 = factorTimesDensKern*lin%comsr%auxarray(ind+2)
!!!!                        tt3 = factorTimesDensKern*lin%comsr%auxarray(ind+3)
!!!!                        rho(indLarge  ) = rho(indLarge  ) + tt0
!!!!                        rho(indLarge+1) = rho(indLarge+1) + tt1
!!!!                        rho(indLarge+2) = rho(indLarge+2) + tt2
!!!!                        rho(indLarge+3) = rho(indLarge+3) + tt3
!!!!                        totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
!!!!                        ind=ind+4
!!!!                    end do
!!!!                end if
!!!!            end do
!!!!        end do
!!!!    end do
!!!!end do
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!call cpu_time(t2)
!!!!time=t2-t1
!!!!if(iproc==0) write(*,'(a,es12.4)') 'time for large loop:',time
!!!!
!!!!call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
!!!!if(iproc==0) write(*,'(1x,a,es20.12)') 'done. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh
!!!!
!!!!iall=-product(shape(densKern))*kind(densKern)
!!!!deallocate(densKern, stat=istat)
!!!!call memocc(istat, iall, 'densKern', subname)
!!!!
!!!!
!!!!end subroutine sumrholinear_withauxiliary

!!subroutine setNscatterarr(iproc,nproc,datacode,atoms,n1i,n2i,n3i,ixc,nscatterarr)
!!  use module_base
!!  use module_types
!!  use Poisson_Solver
!!  use module_xc
!!  implicit none
!!  !Arguments
!!  character(len=1), intent(in) :: datacode
!!  integer, intent(in) :: iproc,nproc,ixc,n1i,n2i,n3i
!!  real(gp), intent(in) :: hxh,hyh,hzh
!!  type(atoms_data), intent(in) :: atoms
!!  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
!!  !Local Variables
!!   integer :: jproc
!!
!!  if (datacode == 'D') then
!!     do jproc=0,iproc-1
!!        call PS_dim4allocation(atoms%geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,n3d,n3p,n3pi,i3xcsh,i3s)
!!        nscatterarr(jproc,1)=n3d            !number of planes for the density
!!        nscatterarr(jproc,2)=n3p            !number of planes for the potential
!!        nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
!!        nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
!!     end do
!!     do jproc=iproc+1,nproc-1
!!        call PS_dim4allocation(atoms%geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,n3d,n3p,n3pi,i3xcsh,i3s)
!!        nscatterarr(jproc,1)=n3d
!!        nscatterarr(jproc,2)=n3p
!!        nscatterarr(jproc,3)=i3s+i3xcsh-1
!!        nscatterarr(jproc,4)=i3xcsh
!!     end do
!!  end if
!!
!!  call PS_dim4allocation(atoms%geocode,datacode,iproc,nproc,n1i,n2i,n3i,ixc,n3d,n3p,n3pi,i3xcsh,i3s)
!!  nscatterarr(iproc,1)=n3d
!!  nscatterarr(iproc,2)=n3p
!!  nscatterarr(iproc,3)=i3s+i3xcsh-1
!!  nscatterarr(iproc,4)=i3xcsh
!!
!!end subroutine setNscatterarr

!!!!subroutine postCommunicationSumrho2(iproc, nproc, comsr, sendBuf, recvBuf)
!!!!use module_base
!!!!use module_types
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc
!!!!!type(p2pCommsSumrho),intent(inout):: comsr
!!!!type(p2pComms),intent(inout):: comsr
!!!!real(8),dimension(comsr%nsendBuf),intent(inout):: sendBuf
!!!!real(8),dimension(comsr%nrecvBuf),intent(out):: recvBuf
!!!!
!!!!! Local variables
!!!!integer:: jproc, nreceives, nsends, iorb, mpisource, istsource, ncount, lrsource, mpidest, istdest, tag, ierr
!!!!integer:: ist, istr, ilr
!!!!
!!!!
!!!!! Communicate the orbitals for the calculation of the charge density.
!!!!! Since we use non blocking point to point communication, only post the message
!!!!! and continues with other calculations.
!!!!! Be aware that you must not modify the send buffer without checking whether
!!!!! the communications has completed.
!!!!if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for the calculation of the charge density... '
!!!!nreceives=0
!!!!nsends=0
!!!!comsr%communComplete=.false.
!!!!procLoop1: do jproc=0,nproc-1
!!!!    orbsLoop1: do iorb=1,comsr%noverlaps(jproc)
!!!!        mpisource=comsr%comarr(1,iorb,jproc)
!!!!        istsource=comsr%comarr(2,iorb,jproc)
!!!!        ncount=comsr%comarr(3,iorb,jproc)
!!!!        lrsource=comsr%comarr(4,iorb,jproc)
!!!!        mpidest=comsr%comarr(5,iorb,jproc)
!!!!        istdest=comsr%comarr(6,iorb,jproc)
!!!!        tag=comsr%comarr(7,iorb,jproc)
!!!!        if(ncount==0) then
!!!!            ! No communication is needed. This should be improved in the initialization, i.e. this communication
!!!!            ! with 0 elements should be removed from comgp%noverlaps etc.
!!!!            comsr%comarr(8,iorb,jproc)=mpi_request_null
!!!!            comsr%comarr(9,iorb,jproc)=mpi_request_null
!!!!            comsr%communComplete(iorb,jproc)=.true.
!!!!            if(iproc==mpidest) then
!!!!                ! This is just to make the check at the end happy.
!!!!                nreceives=nreceives+1
!!!!            end if
!!!!        else
!!!!            if(mpisource/=mpidest) then
!!!!                ! The orbitals are on different processes, so we need a point to point communication.
!!!!                if(iproc==mpisource) then
!!!!                    !write(*,'(6(a,i0))') 'sumrho: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
!!!!                    !call mpi_isend(lphi(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, comsr%comarr(8,iorb,jproc), ierr)
!!!!                    call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
!!!!                         comsr%comarr(8,iorb,jproc), ierr)
!!!!                    comsr%comarr(9,iorb,jproc)=mpi_request_null !is this correct?
!!!!                    nsends=nsends+1
!!!!                else if(iproc==mpidest) then
!!!!                   !write(*,'(6(a,i0))') 'sumrho: process ', mpidest, ' receives ', ncount, &
!!!!                   !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!!                    call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
!!!!                         comsr%comarr(9,iorb,jproc), ierr)
!!!!                    comsr%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
!!!!                    nreceives=nreceives+1
!!!!                else
!!!!                    comsr%comarr(8,iorb,jproc)=mpi_request_null
!!!!                    comsr%comarr(9,iorb,jproc)=mpi_request_null
!!!!                end if
!!!!            else
!!!!                ! The orbitals are on the same process, so simply copy them.
!!!!                if(iproc==mpisource) then
!!!!                    !write(*,'(6(a,i0))') 'sumrho: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
!!!!                    call dcopy(ncount, sendBuf(istsource), 1, recvBuf(istdest), 1)
!!!!                    comsr%comarr(8,iorb,jproc)=mpi_request_null
!!!!                    comsr%comarr(9,iorb,jproc)=mpi_request_null
!!!!                    nsends=nsends+1
!!!!                    nreceives=nreceives+1
!!!!                    comsr%communComplete(iorb,mpisource)=.true.
!!!!                else
!!!!                    comsr%comarr(8,iorb,jproc)=mpi_request_null
!!!!                    comsr%comarr(9,iorb,jproc)=mpi_request_null
!!!!                    comsr%communComplete(iorb,mpisource)=.true.
!!!!                end if
!!!!            end if
!!!!        end if
!!!!    end do orbsLoop1
!!!!end do procLoop1
!!!!if(iproc==0) write(*,'(a)') 'done.'
!!!!
!!!!if(nreceives/=comsr%noverlaps(iproc)) then
!!!!    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comsr%noverlaps(iproc)', nreceives,&
!!!!         comsr%noverlaps(iproc)
!!!!    stop
!!!!end if
!!!!call mpi_barrier(mpi_comm_world, ierr)
!!!!
!!!!end subroutine postCommunicationSumrho2

!!!!> Initializes the parameters needed for the communication of the orbitals
!!!!! when calculating the charge density.
!!!!!
!!!!! input arguments
!!!!!  @param jproc        process to which the orbital shall be sent
!!!!!  @param iorb         orbital that is to be sent
!!!!!  @param istDest      the position on the MPI process to which it should be sent
!!!!!  @param tag          communication tag
!!!!!  @param lin          type containing the parameters for the linear scaling version
!!!!! output arguments
!!!!!  @param commsSumrho  contains the parameters
!!!subroutine setCommunicationInformation(jproc, iorb, istDest, tag, lin, commsSumrho)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: jproc, iorb, istDest, tag
!!!type(linearParameters),intent(in):: lin
!!!integer,dimension(9),intent(out):: commsSumrho
!!!
!!!! Local variables
!!!integer:: mpisource, ist, jorb, jlr
!!!
!!!! on which MPI process is the orbital that has to be sent to jproc
!!!mpisource=lin%orbs%onWhichMPI(iorb)
!!!commsSumrho(1)=mpisource
!!!
!!!! starting index of the orbital on that MPI process
!!!ist=1
!!!do jorb=lin%orbs%isorb_par(mpisource)+1,iorb-1
!!!    !jlr=lin%onWhichAtomAll(jorb)
!!!    jlr=lin%orbs%inWhichLocreg(jorb)
!!!    ist=ist+lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
!!!end do
!!!commsSumrho(2)=ist
!!!
!!!! amount of data to be sent
!!!!jlr=lin%onWhichAtomAll(iorb)
!!!jlr=lin%orbs%inWhichLocreg(iorb)
!!!commsSumrho(3)=lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
!!!
!!!! localization region to which this orbital belongs to
!!!!commsSumrho(4)=lin%onWhichAtomAll(iorb)
!!!commsSumrho(4)=lin%orbs%inWhichLocreg(iorb)
!!!
!!!! to which MPI process should this orbital be sent
!!!commsSumrho(5)=jproc
!!!
!!!! the position on the MPI process to which it should be sent
!!!commsSumrho(6)=istDest
!!!
!!!! the tag for this communication
!!!commsSumrho(7)=tag
!!!
!!!! commsSumrho(8): this entry is used a request for the mpi_isend.
!!!
!!!! commsSumrho(9): this entry is used a request for the mpi_irecv.
!!!
!!!
!!!end subroutine setCommunicationInformation


!!!!!> Initializes the parameters needed for the communication of the orbitals
!!!!!! when calculating the charge density.
!!!!!!
!!!!!! input arguments
!!!!!!  @param jproc        process to which the orbital shall be sent
!!!!!!  @param iorb         orbital that is to be sent
!!!!!!  @param istDest      the position on the MPI process to which it should be sent
!!!!!!  @param tag          communication tag
!!!!!!  @param lin          type containing the parameters for the linear scaling version
!!!!!! output arguments
!!!!!!  @param commsSumrho  contains the parameters
!!!!subroutine setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, istDest, tag, nlr, Llr, &
!!!!           onWhichAtomAll, orbs, commsSumrho)
!!!!use module_base
!!!!use module_types
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: jproc, iorb, is3ovrlp, n3ovrlp, istDest, tag, nlr
!!!!type(locreg_descriptors),dimension(nlr),intent(in):: Llr
!!!!type(orbitals_data):: orbs
!!!!integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
!!!!integer,dimension(6),intent(out):: commsSumrho
!!!!
!!!!! Local variables
!!!!integer:: mpisource, ist, jorb, jlr
!!!!
!!!!! on which MPI process is the orbital that has to be sent to jproc
!!!!mpisource=orbs%onWhichMPI(iorb)
!!!!commsSumrho(1)=mpisource
!!!!
!!!!! starting index of the orbital on that MPI process
!!!!ist=1
!!!!do jorb=orbs%isorb_par(mpisource)+1,iorb-1
!!!!    jlr=onWhichAtomAll(jorb)
!!!!    !ist=ist+lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
!!!!    ist = ist + Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*Llr(jlr)%d%n3i
!!!!end do
!!!!jlr=onWhichAtomAll(iorb)
!!!!ist = ist + Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*(is3ovrlp-1)
!!!!commsSumrho(2)=ist
!!!!
!!!!! amount of data to be sent
!!!!jlr=onWhichAtomAll(iorb)
!!!!!commsSumrho(3)=lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
!!!!commsSumrho(3)=Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*n3ovrlp
!!!!
!!!!!!! localization region to which this orbital belongs to
!!!!!!commsSumrho(4)=onWhichAtomAll(iorb)
!!!!
!!!!! to which MPI process should this orbital be sent
!!!!!commsSumrho(5)=jproc
!!!!commsSumrho(4)=jproc
!!!!
!!!!! the position on the MPI process to which it should be sent
!!!!!commsSumrho(6)=istDest
!!!!commsSumrho(5)=istDest
!!!!
!!!!! the tag for this communication
!!!!!commsSumrho(7)=tag
!!!!commsSumrho(6)=tag
!!!!
!!!!! commsSumrho(8): this entry is used as request for the mpi_isend.
!!!!
!!!!! commsSumrho(9): this entry is used as request for the mpi_irecv.
!!!!
!!!!
!!!!end subroutine setCommunicationInformation2

subroutine calculate_density_kernel(iproc, nproc, ld_coeff, orbs, orbs_tmb, coeff, kernel,  ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, ld_coeff
  type(orbitals_data),intent(in):: orbs, orbs_tmb
  real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
  real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out):: kernel
  real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),optional,intent(in):: ovrlp

  ! Local variables
  integer:: istat, iall, ierr
  real(8),dimension(:,:),allocatable:: density_kernel_partial
  character(len=*),parameter:: subname='calculate_density_kernel'
!!  integer :: i, j, lwork, k
!!  real(8),dimension(:,:), allocatable :: rhoprime, ks, ksk, ksksk,vr,vl
!!  real(8),dimension(:),allocatable:: eval,eval1,beta
!!  real(8),dimension(:),allocatable:: work
!!  real(8),dimension(norb,norb) :: kernelij
!!  real(8),dimension(ld_coeff,norb):: ocoeff
!!  real(dp) :: av_idem, av_k_sym_diff, tr, tr2
!!  integer :: num_counted
  call timing(iproc,'calc_kernel','ON') !lr408t

  if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
  allocate(density_kernel_partial(orbs_tmb%norb,max(orbs_tmb%norbp,1)), stat=istat)
  call memocc(istat, density_kernel_partial, 'density_kernel_partial', subname)
  if(orbs_tmb%norbp>0) then
      call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff(1,1), ld_coeff, &
           coeff(orbs_tmb%isorb+1,1), ld_coeff, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
  end if
  call timing(iproc,'calc_kernel','OF') !lr408t

  if (nproc > 1) then
     call timing(iproc,'commun_kernel','ON') !lr408t
     call mpi_allgatherv(density_kernel_partial(1,1), orbs_tmb%norb*orbs_tmb%norbp, mpi_double_precision, &
          kernel(1,1), orbs_tmb%norb*orbs_tmb%norb_par(:,0), orbs_tmb%norb*orbs_tmb%isorb_par, mpi_double_precision, &
          mpi_comm_world, ierr)
     call timing(iproc,'commun_kernel','OF') !lr408t
  else
     call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,kernel(1,1),1)
  end if

  iall=-product(shape(density_kernel_partial))*kind(density_kernel_partial)
  deallocate(density_kernel_partial,stat=istat)
  call memocc(istat,iall,'density_kernel_partial',subname)

!!$  ! calculate kernelij and print
!!!  if (present(ovrlp)) then
!!$     call gemm('t', 'n', norb_tmb, norb, norb_tmb, 1.d0, ovrlp(1,1), ld_coeff, &
!!$          coeff(1,1), ld_coeff, 0.d0, ocoeff(1,1), norb_tmb)
!!$  else
!!$     call dcopy(ld_coeff*norb,coeff(1,1),1,ocoeff(1,1),1)
!!$  end if
!!$
!!$  call gemm('t', 'n', norb, norb, norb_tmb, 1.d0, coeff(1,1), ld_coeff, &
!!$       ocoeff(1,1), ld_coeff, 0.d0, kernelij(1,1), norb)
!!$
!!$  !!open(33,file='kernelij.dat',status='replace')
!!$  !!do i=1,norb
!!$  !!do j=1,norb
!!$  !!  write(33,*) i,j,kernelij(i,j)
!!$  !!end do
!!$  !!end do
!!$  !!close(33)
!!$
!!$  ! calculate eigenvalues
!!$if (present(ovrlp)) then
!!$
!!$if (.false.) then
!!$  allocate(rhoprime(norb_tmb,norb_tmb))
!!$
!!$  allocate(vl(1:norb_tmb,1:norb_tmb))
!!$  allocate(vr(1:norb_tmb,1:norb_tmb))
!!$  allocate(eval1(1:norb_tmb))
!!$  allocate(beta(1:norb_tmb))
!!$
!!$  allocate(eval(1:norb_tmb))
!!$  allocate(work(1:1))
!!$  call dcopy(norb_tmb*norb_tmb,kernel(1,1),1,rhoprime(1,1),1)
!!$
!!$  lwork = -1
!!$  !call dsygv(1, 'v', 'l',norb_tmb,&
!!$  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, work, lwork, ierr)
!!$  !call dggev('n', 'n',norb_tmb,&
!!$  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, eval1, beta, &
!!$  !      vl,1,vr,1,work, lwork, ierr)
!!$      call DGEEV( 'v','v', norb_tmb, rhoprime(1,1), norb_tmb, eval, eval1, VL, norb_tmb, VR,&
!!$                  norb_tmb, WORK, LWORK, ierr )
!!$
!!$  lwork = work(1)
!!$  deallocate(work)
!!$  allocate(work(1:lwork))
!!$
!!$  !call dsygv(1, 'v', 'l',norb_tmb,&
!!$  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, work, lwork, ierr)
!!$  !call dggev('n', 'n',norb_tmb,&
!!$  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, eval1, beta, &
!!$  !      vl,1,vr,1,work, lwork, ierr)
!!$      call DGEEV( 'v','v', norb_tmb, rhoprime(1,1), norb_tmb, eval, eval1, VL, norb_tmb, VR,&
!!$                  norb_tmb, WORK, LWORK, ierr )
!!$
!!$  if (ierr /= 0) print*,'Error in kernel diag, info=',ierr
!!$  write(32,*) 'Before',eval
!!$  write(32,*) 'eval1',eval1
!!$  !write(32,*) 'beta',beta
!!$  write(32,*) sum(eval)  
!!$  write(32,*) ''
!!$  !call dcopy(norb_tmb*norb_tmb,rhoprime(1,1),1,kernel(1,1),1)
!!$
!!$  deallocate(work)
!!$  deallocate(eval)
!!$end if
!!$
!!$  !!!! Check symmetry and idempotency of kernel
!!!  allocate(ks(1:norb_tmb,1:norb_tmb))
!!!  allocate(ksk(1:norb_tmb,1:norb_tmb))
!!!  allocate(ksksk(1:norb_tmb,1:norb_tmb))
!!$
!!$  !!av_k_sym_diff = 0.0_dp
!!$  !!num_counted = 0
!!$  !!do i=1,norb_tmb
!!$  !!do j=i+1,norb_tmb
!!$  !!   av_k_sym_diff = av_k_sym_diff + abs(kernel(j,i)-kernel(i,j))
!!$  !!   kernel(i,j) = 0.5_dp * (kernel(j,i) + kernel(i,j))
!!$  !!   kernel(j,i) = 0.5_dp * (kernel(j,i) + kernel(i,j))
!!$  !!   num_counted = num_counted + 1
!!$  !!   if (abs(kernel(j,i)-kernel(i,j)) > 1.0d-3) then !lr408 debug
!!$  !!      if (iproc == 0) write(47,*) 'K not symmetric',i,j,kernel(j,i),kernel(i,j),&
!!$  !!           kernel(j,i)-kernel(i,j)
!!$  !!   end if
!!$  !!end do
!!$  !!end do
!!$  !!!!if (iproc == 0) write (48,*) 'av_k_sym_diff',av_k_sym_diff / num_counted
!!$
!!$  ! purification
!!!  do k=1,1
!!!    call dgemm('n','t', norb_tmb,norb_tmb,norb_tmb,1.d0,kernel(1,1),norb_tmb,&
!!!         ovrlp(1,1),norb_tmb,0.d0,ks(1,1),norb_tmb)
!!!     !call dcopy(norb_tmb*norb_tmb,kernel(1,1),1,ks(1,1),1)
!!!
!!!     call dgemm('n','t', norb_tmb,norb_tmb,norb_tmb,1.d0,ks(1,1),norb_tmb,&
!!!          kernel(1,1),norb_tmb,0.d0,ksk(1,1),norb_tmb)
!!!
!!!     !call dgemm('n','t', norb_tmb,norb_tmb,norb_tmb,1.d0,ks(1,1),norb_tmb,&
!!!     !     ksk(1,1),norb_tmb,0.d0,ksksk(1,1),norb_tmb)
!!!
!!!     !kernel = 3.0_dp * ksk - 2.0_dp * ksksk
!!!
!!!     !do i=1,norb_tmb
!!!     !do j=i,norb_tmb
!!!     !  kernel(i,j) = 0.5_dp * (kernel(j,i) + kernel(i,j))
!!!     !  kernel(j,i) = 0.5_dp * (kernel(j,i) + kernel(i,j))
!!!     !end do
!!!     !end do
!!!
!!!     !DEBUG test idempotency of density matrix
!!!     !call dgemm('n','t', norb_tmb,norb_tmb,norb_tmb,1.d0,kernel(1,1),norb_tmb,&
!!!     !     kernel(1,1),norb_tmb,0.d0,rhoprime(1,1),norb_tmb)
!!!     av_idem = 0.0_dp
!!!     open(31)
!!!     do i = 1, norb_tmb
!!!       do j = 1, norb_tmb
!!!         av_idem = av_idem + abs(ksk(i,j) - kernel(i,j))
!!!         if(abs(ksk(i,j) - kernel(i,j))>1.0d-3) then
!!!           if (iproc==0) write(31,*) 'Not indempotent',i,j,ksk(i,j), kernel(i,j),ksk(i,j) - kernel(i,j)
!!!         end if
!!!       end do
!!!       tr = tr + ks(i,i)
!!!       tr2 = tr2 + kernel(i,i)
!!!     end do
!!!     close(31)
!!!     if (iproc == 0) write(49,*) 'av_idem',av_idem / (norb_tmb **2)
!!!     if (iproc == 0) write(49,*) 'trace',tr,tr2
!!!
!!!  end do

!!$if (.false.) then
!!$  allocate(eval(1:norb_tmb))
!!$  allocate(work(1:1))
!!$
!!$  deallocate(beta)
!!$  deallocate(eval1)
!!$  deallocate(vr)
!!$  deallocate(vl)
!!$
!!$!  call dcopy(norb_tmb*norb_tmb,kernel(1,1),1,rhoprime(1,1),1)
!!$
!!$!  lwork = -1
!!$
!!$!  allocate(vl(1:norb_tmb,1:norb_tmb))
!!$!  allocate(vr(1:norb_tmb,1:norb_tmb))
!!$!  allocate(eval1(1:norb_tmb))
!!$!  allocate(beta(1:norb_tmb))
!!$
!!$!  !call dsygv(1, 'v', 'l',norb_tmb,&
!!$!  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, work, lwork, ierr)
!!$!  !call dggev('n', 'n',norb_tmb,&
!!$!  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, eval1, beta, &
!!$!   !     vl,1,vr,1,work, lwork, ierr)
!!$!      call DGEEV( 'v','v', norb_tmb, rhoprime(1,1), norb_tmb, eval, eval1, VL, norb_tmb, VR,&
!!$!                  norb_tmb, WORK, LWORK, ierr )
!!$
!!$!  lwork = work(1)
!!$!  deallocate(work)
!!$!  allocate(work(1:lwork))
!!$
!!$!  !call dsygv(1, 'v', 'l',norb_tmb,&
!!$!  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, work, lwork, ierr)
!!$!  !call dggev('n', 'n',norb_tmb,&
!!$!  !      rhoprime(1,1), norb_tmb, ovrlp, norb_tmb, eval, eval1, beta, &
!!$!  !      vl,1,vr,1,work, lwork, ierr)
!!$!      call DGEEV( 'v','v', norb_tmb, rhoprime(1,1), norb_tmb, eval, eval1, VL, norb_tmb, VR,&
!!$!                  norb_tmb, WORK, LWORK, ierr )
!!$
!!$!  if (ierr /= 0) print*,'Error in kernel diag, info=',ierr
!!$!  write(32,*) 'After',eval
!!$!  write(32,*) 'eval1',eval1
!!$!  !write(32,*) 'beta',beta
!!$!  write(32,*) sum(eval)  
!!$!  write(32,*) ''
!!$!  !call dcopy(norb_tmb*norb_tmb,rhoprime(1,1),1,kernel(1,1),1)
!!$
!!$!  deallocate(work)
!!$!  deallocate(eval)
!!$
!!$!  deallocate(beta)
!!$!  deallocate(eval1)
!!$!  deallocate(vr)
!!$!  deallocate(vl)
!!$
!!$  deallocate(rhoprime)
!!$end if
!!$
!!$  !!!deallocate(ksksk)
!!$  !!!deallocate(ksk)
!!$  !!!deallocate(ks)
!!$  !END DEBUG
!!$
!!$  !write(*,*) 'DEBUG KERNEL'
!!$  !do i=1,norb_tmb
!!$  !    do j=1,norb_tmb
!!$  !        if(j==i) then
!!$  !            kernel(j,i)=1.d0
!!$  !        else
!!$  !            kernel(j,i)=0.d0
!!$  !        end if
!!$  !    end do
!!$  !end do
!!$
!!!end if

!!!if (iproc==0) then ! lr408
!!!   open(127,file='kernel.dat',status='replace')
!!!   do i=1,norb_tmb
!!!   do j=1,norb_tmb
!!!      write(127,*) i,j,kernel(i,j)
!!!   end do
!!!   write(127,*) ''
!!!   end do
!!!   write(127,*) ''
!!!   close(127)
!!!end if

  !call timing(iproc,'sumrho_TMB    ','OF')

end subroutine calculate_density_kernel



!!subroutine calculate_charge_density(iproc, nproc, lzd, hxh, hyh, hzh, orbs, densKern, factor, &
!!           nscatterarr, maxval_noverlaps, noverlaps, overlaps, comarr, ise3, nrecvBuf, recvBuf, nrho, rho)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, nrho, maxval_noverlaps, nrecvBuf
!!  real(8),intent(in):: hxh, hyh, hzh, factor
!!  type(local_zone_descriptors),intent(in) :: lzd
!!  type(orbitals_data),intent(in) :: orbs
!!  real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: densKern
!!  integer, dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!  integer,dimension(0:nproc-1),intent(in):: noverlaps
!!  integer,dimension(noverlaps(iproc)),intent(in):: overlaps
!!  integer,dimension(6,maxval_noverlaps,0:nproc-1),intent(in):: comarr
!!  integer,dimension(noverlaps(iproc),2),intent(in):: ise3
!!  real(8),dimension(nrecvBuf),intent(in):: recvBuf
!!  real(kind=8),dimension(nrho),intent(out),target :: rho
!!
!!  ! Local variables
!!  integer:: is, ie, iorb, iiorb, ilr, istri, jorb, jjorb, jlr, istrj, azones, bzones, ii, izones, jzones
!!  integer:: i1s, i1e, i2s, i2e, i3s, i3e, i3, i3d, j3d, indi3, indj3, indl3, i2, i2d, j2d, indi2, indj2, indl2
!!  integer:: m, i1d0, j1d0, indri0, indrj0, indLarge0, i1, i1d, j1d, indLarge
!!  integer:: x,y,z,ishift1,ishift2,ishift3,jshift1,jshift2,jshift3, indri, indrj, ierr
!!  real(kind=8):: totalCharge, tt, tt0, tt1, tt2, tt3, factorTimesDensKern
!!  integer,dimension(3,4) :: astart, aend, bstart,bend
!!
!!  call timing(iproc,'sumrho_TMB    ','ON')
!!  
!!  ! Bounds of the slice in global coordinates.
!!  is=nscatterarr(iproc,3) 
!!  ie=is+nscatterarr(iproc,1)-1
!!  
!!  ! This sum is "symmetric", so only do the second loop (jorb) only up to iorb and multiply by two if iorb/=jorb.
!!  totalCharge=0.d0
!!  !print*,'noverlaps(iproc)',iproc,noverlaps(iproc)
!!  do iorb=1,noverlaps(iproc)
!!      iiorb=overlaps(iorb) !global index of orbital iorb
!!      ilr=orbs%inwhichlocreg(iiorb) !localization region of orbital iorb
!!      istri=comarr(5,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
!!      do jorb=1,iorb
!!          jjorb=overlaps(jorb) !global indes of orbital jorb
!!          jlr=orbs%inwhichlocreg(jjorb) !localization region of orbital jorb
!!          istrj=comarr(5,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
!!          !!tt = (lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2 &
!!          !!    +(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2 &
!!          !!    +(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
!!          !!tt=sqrt(tt)
!!          !!write(200,*) tt, densKern(iiorb,jjorb)
!!          !!if(tt>6.d0) cycle
!!  
!!          azones = 1
!!          bzones = 1
!!          !Calculate the number of regions to cut alr and blr
!!          do ii=1,2
!!             if(lzd%llr(ilr)%outofzone(ii) > 0) azones = azones * 2
!!             if(lzd%llr(jlr)%outofzone(ii) > 0) bzones = bzones * 2
!!          end do
!!  
!!          !FRACTURE THE FIRST LOCALIZATION REGION
!!          call fracture_periodic_zone_ISF(azones,lzd%Glr,lzd%Llr(ilr),lzd%Llr(ilr)%outofzone(:),astart,aend)
!!         
!!          !FRACTURE SECOND LOCREG
!!          call fracture_periodic_zone_ISF(bzones,lzd%Glr,lzd%Llr(jlr),lzd%Llr(jlr)%outofzone(:),bstart,bend)
!!          do izones=1,azones
!!             do jzones=1,bzones
!!                ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
!!                i1s=max(astart(1,izones),bstart(1,jzones))
!!                i1e=min(aend(1,izones)-1,bend(1,jzones)-1)
!!                i2s=max(astart(2,izones),bstart(2,jzones))
!!                i2e=min(aend(2,izones)-1,bend(2,jzones)-1)
!!                i3s=max(ise3(iorb,1),ise3(jorb,1))
!!                i3e=min(ise3(iorb,2),ise3(jorb,2))
!!                call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(ilr),x,y,z,ishift1, ishift2, ishift3)
!!                call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(jlr),x,y,z,jshift1, jshift2, jshift3)
!!                factorTimesDensKern = factor*densKern(iiorb,jjorb)
!!                if(iorb/=jorb) then
!!                    ! Multiply by two since these elements appear twice (but are calculated only once).
!!                    factorTimesDensKern = 2.d0*factorTimesDensKern
!!                end if
!!                ! Now loop over all points in the box in which the orbitals overlap.
!!                !if(i3s>i3e) write(*,*) 'no calculation done'
!!  
!!                do i3=i3s,i3e !bounds in z direction
!!                    i3d=i3 -max(is,-ishift3) !z coordinate of orbital iorb with respect to the overlap box
!!                    j3d=i3 -max(is,-jshift3) !z coordinate of orbital jorb with respect to the overlap box
!!                    indi3=i3d*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
!!                    indj3=j3d*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
!!                    !indl3=(i3-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!                    !indl3=(modulo(i3-1,Lzd%Glr%d%n3i)-is+1)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!                    indl3=(modulo(i3,Lzd%Glr%d%n3i+1)-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!                    do i2=i2s,i2e !bounds in y direction
!!                        i2d=i2 + ishift2 !y coordinate of orbital iorb with respect to the overlap box
!!                        j2d=i2 + jshift2 !y coordinate of orbital jorb with respect to the overlap box
!!                        indi2=i2d*lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
!!                        indj2=j2d*lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
!!                        !indl2=i2*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!                        !indl2=(modulo(i2-1,Lzd%Glr%d%n2i)+1)*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!                        indl2=(modulo(i2,Lzd%Glr%d%n2i+1))*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!                        ! For all other than free BC, choose m such that the unrolled part is never used.
!!                        if(Lzd%Glr%geocode=='F') then
!!                            m=mod(i1e-i1s+1,4)
!!                        else
!!                            m=i1e-i1s+1
!!                        end if
!!                        if(m/=0) then
!!                            ! The following five variables hold some intermediate results to speed up the code.
!!                            i1d0= ishift1 
!!                            j1d0= jshift1
!!                            indri0 = indi3 + indi2 + istri + 1
!!                            indrj0 = indj3 + indj2 + istrj + 1
!!                            indLarge0 = indl3 + indl2 + 1   
!!                            do i1=i1s,i1s+m-1
!!                                i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
!!                                j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
!!                                indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
!!                                indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
!!                                !indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
!!                                !indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1 !index for which the charge density is beeing calculated
!!                                indLarge = indLarge0 + modulo(i1,Lzd%Glr%d%n1i+1) !index for which the charge density is beeing calculated
!!                                tt = factorTimesDensKern*recvBuf(indri)*recvBuf(indrj)
!!                                rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
!!                                totalCharge = totalCharge + tt !add the contribution to the total charge
!!                            end do
!!                        end if
!!                        ! This is the same again, this time with unrolled loops.
!!                        if(i1e-i1s+1>4) then
!!                            i1d0= ishift1 
!!                            j1d0= jshift1
!!                            indri0 = indi3 + indi2 + istri + 1
!!                            indrj0 = indj3 + indj2 + istrj + 1
!!                            indLarge0 = indl3 + indl2 + 1
!!                            do i1=i1s+m,i1e,4
!!                                i1d=i1d0+i1
!!                                j1d=j1d0+i1
!!                                indri = indri0 + i1d
!!                                indrj = indrj0 + j1d
!!                                !indLarge = indLarge0 + i1
!!                                tt0 = factorTimesDensKern*recvBuf(indri  )*recvBuf(indrj  )
!!                                tt1 = factorTimesDensKern*recvBuf(indri+1)*recvBuf(indrj+1)
!!                                tt2 = factorTimesDensKern*recvBuf(indri+2)*recvBuf(indrj+2)
!!                                tt3 = factorTimesDensKern*recvBuf(indri+3)*recvBuf(indrj+3)
!!                                !indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1
!!                                indLarge = indLarge0 + i1
!!                                rho(indLarge  ) = rho(indLarge  ) + tt0
!!                                !indLarge = indLarge0 +modulo(i1,Lzd%Glr%d%n1i)+1
!!                                indLarge = indLarge0 + i1+1
!!                                rho(indLarge) = rho(indLarge) + tt1
!!                                !rho(indLarge+1) = rho(indLarge+1) + tt1
!!                                !indLarge = indLarge0 + modulo(i1+1,Lzd%Glr%d%n1i)+1
!!                                indLarge = indLarge0 + i1+2
!!                                rho(indLarge) = rho(indLarge) + tt2
!!                                !rho(indLarge+2) = rho(indLarge+2) + tt1
!!                                !indLarge = indLarge0 + modulo(i1+2,Lzd%Glr%d%n1i)+1
!!                                indLarge = indLarge0 + i1+3
!!                                rho(indLarge) = rho(indLarge) + tt3
!!                                !rho(indLarge+3) = rho(indLarge+3) + tt1
!!                                totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
!!                            end do
!!                        end if
!!                    end do
!!                end do
!!            end do !jzones
!!         end do !izones
!!      end do
!!  end do
!!  !!call mpi_barrier(mpi_comm_world, ierr)
!!  !!call cpu_time(t2)
!!  !!time=t2-t1
!!  
!!  call timing(iproc,'sumrho_TMB    ','OF')
!!  
!!  call timing(iproc,'sumrho_allred','ON')
!!  
!!  call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
!!  if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh
!!  
!!  call timing(iproc,'sumrho_allred','OF')
!!
!!
!!end subroutine calculate_charge_density



!!!subroutine sumrhoForLocalizedBasis2(iproc,nproc,lzd,orbs,&
!!!     comsr,densKern,nrho,rho,at,nscatterarr)
!!!!
!!!use module_base
!!!use module_types
!!!use libxc_functionals
!!!use module_interfaces, exceptThisOne => sumrhoForLocalizedBasis2
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in) :: iproc, nproc, nrho
!!!type(local_zone_descriptors),intent(in) :: lzd
!!!type(orbitals_data),intent(in) :: orbs
!!!!type(p2pCommsSumrho),intent(inout) :: comsr
!!!type(p2pComms),intent(inout) :: comsr
!!!!real(kind=8),dimension(orbs%norb,norb),intent(in) :: coeff
!!!!real(kind=8),dimension(ld_coeff,norb),intent(in) :: coeff
!!!real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: densKern
!!!real(8),dimension(nrho),intent(out),target:: rho
!!!type(atoms_data),intent(in) :: at
!!!integer, dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!
!!!! Local variables
!!!integer :: iorb, jorb, indLarge, i1, i2, i3, ilr, jlr
!!!integer :: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, istri, istrj
!!!integer :: indi2, indi3, indj2, indj3, indl2, indl3, iiorb, jjorb
!!!integer :: ierr, is, ie
!!!integer :: m, i1d0, j1d0, indri0, indrj0, indLarge0
!!!integer :: azones,bzones,ii,izones,jzones,x,y,z,ishift1,ishift2,ishift3,jshift1,jshift2,jshift3
!!!integer,dimension(3,4) :: astart, aend, bstart,bend
!!!real(kind=8) :: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern
!!!!real(kind=8),dimension(:,:),allocatable :: densKern
!!!!character(len=*),parameter :: subname='sumrhoForLocalizedBasis2'
!!!
!!!if(iproc==0) write(*,'(a)',advance='no') 'Calculating charge density... '
!!!
!!!
!!!
!!!
!!!! Define some constant factors.
!!!hxh=.5d0*lzd%hgrids(1)
!!!hyh=.5d0*lzd%hgrids(2)
!!!hzh=.5d0*lzd%hgrids(3)
!!!factor=1.d0/(hxh*hyh*hzh)
!!!
!!!! Initialize rho.
!!!if (libxc_functionals_isgga()) then
!!!    call razero(nrho, rho)
!!!else
!!!    ! There is no mpi_allreduce, therefore directly initialize to
!!!    ! 10^-20 and not 10^-20/nproc.
!!!    !rho=1.d-20
!!!    call tenminustwenty(nrho, rho, nproc)
!!!end if
!!!call timing(iproc,'p2pSumrho_wait','ON')
!!!
!!!
!!!call wait_p2p_communication(iproc, nproc, comsr)
!!!
!!!
!!!
!!!! Now calculate the charge density. Each process calculates only one slice of the total charge density.
!!!! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
!!!! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
!!!! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
!!!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!!!call cpu_time(t1)
!!!
!!!call timing(iproc,'p2pSumrho_wait','OF')
!!!
!!!
!!!!!call calculate_charge_density(iproc, nproc, lzd, hxh, hyh, hzh, orbs, densKern, factor, &
!!!!!           nscatterarr, maxval(comsr%noverlaps), comsr%noverlaps, comsr%overlaps, comsr%comarr, comsr%ise3, comsr%nrecvBuf, comsr%recvBuf, nrho, rho)
!!!
!!!call timing(iproc,'sumrho_TMB    ','ON')
!!!
!!!! Bounds of the slice in global coordinates.
!!!is=nscatterarr(iproc,3) 
!!!ie=is+nscatterarr(iproc,1)-1
!!!
!!!! This sum is "symmetric", so only do the second loop (jorb) only up to iorb and multiply by two if iorb/=jorb.
!!!totalCharge=0.d0
!!!!print*,'comsr%noverlaps(iproc)',iproc,comsr%noverlaps(iproc)
!!!do iorb=1,comsr%noverlaps(iproc)
!!!    iiorb=comsr%overlaps(iorb) !global index of orbital iorb
!!!    ilr=orbs%inwhichlocreg(iiorb) !localization region of orbital iorb
!!!    istri=comsr%comarr(5,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
!!!    do jorb=1,iorb
!!!        jjorb=comsr%overlaps(jorb) !global indes of orbital jorb
!!!        jlr=orbs%inwhichlocreg(jjorb) !localization region of orbital jorb
!!!        istrj=comsr%comarr(5,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
!!!        !!tt = (lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2 &
!!!        !!    +(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2 &
!!!        !!    +(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
!!!        !!tt=sqrt(tt)
!!!        !!write(200,*) tt, densKern(iiorb,jjorb)
!!!        !!if(tt>6.d0) cycle
!!!
!!!        azones = 1
!!!        bzones = 1
!!!        !Calculate the number of regions to cut alr and blr
!!!        do ii=1,2
!!!           if(lzd%llr(ilr)%outofzone(ii) > 0) azones = azones * 2
!!!           if(lzd%llr(jlr)%outofzone(ii) > 0) bzones = bzones * 2
!!!        end do
!!!
!!!        !FRACTURE THE FIRST LOCALIZATION REGION
!!!        call fracture_periodic_zone_ISF(azones,lzd%Glr,lzd%Llr(ilr),lzd%Llr(ilr)%outofzone(:),astart,aend)
!!!       
!!!        !FRACTURE SECOND LOCREG
!!!        call fracture_periodic_zone_ISF(bzones,lzd%Glr,lzd%Llr(jlr),lzd%Llr(jlr)%outofzone(:),bstart,bend)
!!!        do izones=1,azones
!!!           do jzones=1,bzones
!!!              ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
!!!              i1s=max(astart(1,izones),bstart(1,jzones))
!!!              i1e=min(aend(1,izones)-1,bend(1,jzones)-1)
!!!              i2s=max(astart(2,izones),bstart(2,jzones))
!!!              i2e=min(aend(2,izones)-1,bend(2,jzones)-1)
!!!              i3s=max(comsr%ise3(iorb,1),comsr%ise3(jorb,1))
!!!              i3e=min(comsr%ise3(iorb,2),comsr%ise3(jorb,2))
!!!              call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(ilr),x,y,z,ishift1, ishift2, ishift3)
!!!              call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(jlr),x,y,z,jshift1, jshift2, jshift3)
!!!              factorTimesDensKern = factor*densKern(iiorb,jjorb)
!!!              if(iorb/=jorb) then
!!!                  ! Multiply by two since these elements appear twice (but are calculated only once).
!!!                  factorTimesDensKern = 2.d0*factorTimesDensKern
!!!              end if
!!!              ! Now loop over all points in the box in which the orbitals overlap.
!!!              !if(i3s>i3e) write(*,*) 'no calculation done'
!!!
!!!              do i3=i3s,i3e !bounds in z direction
!!!                  i3d=i3 -max(is,-ishift3) !z coordinate of orbital iorb with respect to the overlap box
!!!                  j3d=i3 -max(is,-jshift3) !z coordinate of orbital jorb with respect to the overlap box
!!!                  indi3=i3d*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
!!!                  indj3=j3d*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
!!!                  !indl3=(i3-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!!                  !indl3=(modulo(i3-1,Lzd%Glr%d%n3i)-is+1)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!!                  indl3=(modulo(i3,Lzd%Glr%d%n3i+1)-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
!!!                  do i2=i2s,i2e !bounds in y direction
!!!                      i2d=i2 + ishift2 !y coordinate of orbital iorb with respect to the overlap box
!!!                      j2d=i2 + jshift2 !y coordinate of orbital jorb with respect to the overlap box
!!!                      indi2=i2d*lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
!!!                      indj2=j2d*lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
!!!                      !indl2=i2*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!!                      !indl2=(modulo(i2-1,Lzd%Glr%d%n2i)+1)*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!!                      indl2=(modulo(i2,Lzd%Glr%d%n2i+1))*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
!!!                      ! For all other than free BC, choose m such that the unrolled part is never used.
!!!                      if(Lzd%Glr%geocode=='F') then
!!!                          m=mod(i1e-i1s+1,4)
!!!                      else
!!!                          m=i1e-i1s+1
!!!                      end if
!!!                      if(m/=0) then
!!!                          ! The following five variables hold some intermediate results to speed up the code.
!!!                          i1d0= ishift1 
!!!                          j1d0= jshift1
!!!                          indri0 = indi3 + indi2 + istri + 1
!!!                          indrj0 = indj3 + indj2 + istrj + 1
!!!                          indLarge0 = indl3 + indl2 + 1   
!!!                          do i1=i1s,i1s+m-1
!!!                              i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
!!!                              j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
!!!                              indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
!!!                              indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
!!!                              !indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
!!!                              !indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1 !index for which the charge density is beeing calculated
!!!                              indLarge = indLarge0 + modulo(i1,Lzd%Glr%d%n1i+1) !index for which the charge density is beeing calculated
!!!                              tt = factorTimesDensKern*comsr%recvBuf(indri)*comsr%recvBuf(indrj)
!!!                              rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
!!!                              totalCharge = totalCharge + tt !add the contribution to the total charge
!!!                          end do
!!!                      end if
!!!                      ! This is the same again, this time with unrolled loops.
!!!                      if(i1e-i1s+1>4) then
!!!                          i1d0= ishift1 
!!!                          j1d0= jshift1
!!!                          indri0 = indi3 + indi2 + istri + 1
!!!                          indrj0 = indj3 + indj2 + istrj + 1
!!!                          indLarge0 = indl3 + indl2 + 1
!!!                          do i1=i1s+m,i1e,4
!!!                              i1d=i1d0+i1
!!!                              j1d=j1d0+i1
!!!                              indri = indri0 + i1d
!!!                              indrj = indrj0 + j1d
!!!                              !indLarge = indLarge0 + i1
!!!                              tt0 = factorTimesDensKern*comsr%recvBuf(indri  )*comsr%recvBuf(indrj  )
!!!                              tt1 = factorTimesDensKern*comsr%recvBuf(indri+1)*comsr%recvBuf(indrj+1)
!!!                              tt2 = factorTimesDensKern*comsr%recvBuf(indri+2)*comsr%recvBuf(indrj+2)
!!!                              tt3 = factorTimesDensKern*comsr%recvBuf(indri+3)*comsr%recvBuf(indrj+3)
!!!                              !indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1
!!!                              indLarge = indLarge0 + i1
!!!                              rho(indLarge  ) = rho(indLarge  ) + tt0
!!!                              !indLarge = indLarge0 +modulo(i1,Lzd%Glr%d%n1i)+1
!!!                              indLarge = indLarge0 + i1+1
!!!                              rho(indLarge) = rho(indLarge) + tt1
!!!                              !rho(indLarge+1) = rho(indLarge+1) + tt1
!!!                              !indLarge = indLarge0 + modulo(i1+1,Lzd%Glr%d%n1i)+1
!!!                              indLarge = indLarge0 + i1+2
!!!                              rho(indLarge) = rho(indLarge) + tt2
!!!                              !rho(indLarge+2) = rho(indLarge+2) + tt1
!!!                              !indLarge = indLarge0 + modulo(i1+2,Lzd%Glr%d%n1i)+1
!!!                              indLarge = indLarge0 + i1+3
!!!                              rho(indLarge) = rho(indLarge) + tt3
!!!                              !rho(indLarge+3) = rho(indLarge+3) + tt1
!!!                              totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
!!!                          end do
!!!                      end if
!!!                  end do
!!!              end do
!!!          end do !jzones
!!!       end do !izones
!!!    end do
!!!end do
!!!
!!!if(iproc==0) write(*,'(a)') 'done.'
!!!
!!!call timing(iproc,'sumrho_TMB    ','OF')
!!!
!!!call timing(iproc,'sumrho_allred','ON')
!!!
!!!call mpiallred(totalCharge, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!!if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh
!!!
!!!call timing(iproc,'sumrho_allred','OF')
!!!
!!!
!!!
!!!end subroutine sumrhoForLocalizedBasis2



subroutine purify_kernel(orbs_tmb, overlap, kernel)
  !! WARNING: NOT GUARANTIED TO WORK
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs_tmb
  real(kind=8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(in) :: overlap
  real(kind=8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(inout) :: kernel

  ! Local variables
  integer :: istat, iall
  real(kind=8),dimension(:,:),allocatable :: ks, ksk, ksksk
  character(len=*),parameter :: subname='purify_kernel'

   allocate(ks(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
   call memocc(istat, ks, 'ks', subname) 
   allocate(ksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
   call memocc(istat, ksk, 'ksk', subname) 
   allocate(ksksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
   call memocc(istat, ksksk, 'ksksk', subname) 

   call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, kernel(1,1), orbs_tmb%norb, &
              overlap(1,1), orbs_tmb%norb, 0.d0, ks(1,1), orbs_tmb%norb) 
   call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), orbs_tmb%norb, &
              kernel(1,1), orbs_tmb%norb, 0.d0, ksk(1,1), orbs_tmb%norb)
   call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), orbs_tmb%norb, &
              ksk(1,1), orbs_tmb%norb, 0.d0, ksksk(1,1), orbs_tmb%norb)
   print *,'PURIFYING THE KERNEL'
   kernel = 3*ksk-2*ksksk
   
   iall = -product(shape(ks))*kind(ks)
   deallocate(ks,stat=istat)
   call memocc(istat, iall, 'ks', subname)
   iall = -product(shape(ksk))*kind(ksk)
   deallocate(ksk,stat=istat)
   call memocc(istat, iall, 'ksk', subname)
   iall = -product(shape(ksksk))*kind(ksksk)
   deallocate(ksksk,stat=istat)
   call memocc(istat, iall, 'ksksk', subname)
end subroutine purify_kernel
