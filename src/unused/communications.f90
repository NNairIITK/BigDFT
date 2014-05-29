
    subroutine transpose_switch_psirt(collcom_sr, psirt, psirtwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirt
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork
    
      ! Local variables
      integer :: i, ind, sum_c, m
    
      sum_c = sum(collcom_sr%nrecvcounts_c)
    
      !$omp parallel default(private) &
      !$omp shared(collcom_sr, psirt, psirtwork, sum_c, m)
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i=1,m
           ind = collcom_sr%iexpand_c(i)
           psirtwork(ind) = psirt(i)
        end do
      end if
    
    
      !$omp do
      do i=m+1,sum_c,7
          psirtwork(collcom_sr%iexpand_c(i+0))=psirt(i+0)
          psirtwork(collcom_sr%iexpand_c(i+1))=psirt(i+1)
          psirtwork(collcom_sr%iexpand_c(i+2))=psirt(i+2)
          psirtwork(collcom_sr%iexpand_c(i+3))=psirt(i+3)
          psirtwork(collcom_sr%iexpand_c(i+4))=psirt(i+4)
          psirtwork(collcom_sr%iexpand_c(i+5))=psirt(i+5)
          psirtwork(collcom_sr%iexpand_c(i+6))=psirt(i+6)
      end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_switch_psirt



    subroutine transpose_communicate_psirt(iproc, nproc, collcom_sr, psirtwork, psirwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork
    
      ! Local variables
      integer :: ierr
    
      if (nproc>1) then
      call mpi_alltoallv(psirtwork, collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, psirwork, &
           collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
          call vcopy(collcom_sr%ndimpsi_c, psirtwork(1), 1, psirwork(1), 1)
      end if
    
    end subroutine transpose_communicate_psirt



    subroutine transpose_unswitch_psir(collcom_sr, psirwork, psir)
      use module_base
      use module_types
      implicit none
    
      ! Caling arguments
      type(comms_linear),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psir
    
      ! Local variables
      integer :: i, ind, m
    
    
      !$omp parallel default(private) &
      !$omp shared(collcom_sr, psirwork, psir, m)
    
      m = mod(collcom_sr%ndimpsi_c,7)
    
      if(m/=0) then
        do i = 1,m
         ind=collcom_sr%irecvbuf_c(i)
         psir(ind)=psirwork(i)
        end do
      end if
    
      ! coarse part
    
      !$omp do
        do i=m+1,collcom_sr%ndimpsi_c,7
            psir(collcom_sr%irecvbuf_c(i+0))=psirwork(i+0)
            psir(collcom_sr%irecvbuf_c(i+1))=psirwork(i+1)
            psir(collcom_sr%irecvbuf_c(i+2))=psirwork(i+2)
            psir(collcom_sr%irecvbuf_c(i+3))=psirwork(i+3)
            psir(collcom_sr%irecvbuf_c(i+4))=psirwork(i+4)
            psir(collcom_sr%irecvbuf_c(i+5))=psirwork(i+5)
            psir(collcom_sr%irecvbuf_c(i+6))=psirwork(i+6)
        end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_unswitch_psir



    subroutine untoglobal_and_transpose(iproc,nproc,orbs,Lzd,comms,psi,&
         work,outadd) !optional
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc
      type(orbitals_data), intent(in) :: orbs
      type(local_zone_descriptors), intent(in) :: Lzd
      type(comms_cubic), intent(in) :: comms
      real(wp), dimension(:), pointer :: psi !< Input psi should always be in global region, while output psi is in locregs
      real(wp), dimension(:), pointer, optional :: work
      real(wp), dimension(*), intent(out), optional :: outadd !< Optional argument
      !local variables
      character(len=*), parameter :: subname='untoglobal_and_transpose'
      integer :: ierr,i_all,i_stat
      integer :: psishift1,totshift,iorb,ilr,ldim,Gdim
      real(wp), dimension(:), pointer :: workarr
    
      ! Input psi should always be in global region !
      call timing(iproc,'Un-TransSwitch','ON')
    
      if (nproc > 1) then
         !control check
         if (.not. present(work) .or. .not. associated(work)) then
            !if(iproc == 0) 
                 write(*,'(1x,a)')&
                 "ERROR: Unproper work array for untransposing in parallel"
            stop
         end if
         call timing(iproc,'Un-TransSwitch','OF')
         call timing(iproc,'Un-TransComm  ','ON')
         call MPI_ALLTOALLV(psi,comms%ncntt,comms%ndsplt,mpidtypw,  &
              work,comms%ncntd,comms%ndspld,mpidtypw,bigdft_mpi%mpi_comm,ierr)
         call timing(iproc,'Un-TransComm  ','OF')
         call timing(iproc,'Un-TransSwitch','ON')
         if (present(outadd)) then
            !!call unswitch_waves_v(nproc,orbs,&
            !!     Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,comms%nvctr_par(0,1),work,outadd)
            call unswitch_waves_v(nproc,orbs,&
                 Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,comms%nvctr_par,work,outadd)
         else
            !!call unswitch_waves_v(nproc,orbs,&
            !!     Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,comms%nvctr_par(0,1),work,psi)
            call unswitch_waves_v(nproc,orbs,&
                 Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,comms%nvctr_par,work,psi)
         end if
      else
         if(orbs%nspinor /= 1) then
            call psitransspi(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs,psi,.false.)
         end if
      end if
    
      !for linear scaling must project the wavefunctions back into the locregs
      if(Lzd%linear) then
         psishift1 = 1 
         totshift = 0
         Gdim = max((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor,&
               sum(comms%ncntt(0:nproc-1)))
         allocate(workarr(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
         call memocc(i_stat,workarr,'workarr',subname)
         call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp),workarr)
         do iorb=1,orbs%norbp
            ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
            ldim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
            if(present(outadd)) then
                call psi_to_locreg2(iproc, ldim, Gdim, Lzd%Llr(ilr), Lzd%Glr, psi(totshift), outadd(psishift1))
            else
                call psi_to_locreg2(iproc, ldim, Gdim, Lzd%Llr(ilr), Lzd%Glr, psi(totshift), workarr(psishift1))
            end if
            psishift1 = psishift1 + ldim
            totshift = totshift + (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor
         end do
         !reallocate psi to the locreg dimensions
         i_all=-product(shape(psi))*kind(psi)
         deallocate(psi,stat=i_stat)
         call memocc(i_stat,i_all,'psi',subname)
         allocate(psi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
         call memocc(i_stat,psi,'psi',subname)
         call vcopy(max(orbs%npsidim_orbs,orbs%npsidim_comp),workarr(1),1,psi(1),1) !psi=work
         i_all=-product(shape(workarr))*kind(workarr)
         deallocate(workarr,stat=i_stat)
         call memocc(i_stat,i_all,'workarr',subname)
      end if
    
      call timing(iproc,'Un-TransSwitch','OF')
    END SUBROUTINE untoglobal_and_transpose
