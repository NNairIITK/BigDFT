module io
  implicit none

  private

  !> Public routines
  public :: read_linear_matrix_dense
  public :: writeonewave_linear
  public :: writemywaves_linear
  public :: writemywaves_linear_fragments
  public :: read_coeff_minbasis
  public :: io_read_descr_linear

  public :: io_error, io_warning, io_open
  public :: io_read_descr, read_psi_compress
  public :: io_gcoordToLocreg
  public :: read_psig

  contains


    !> Write all my wavefunctions in files by calling writeonewave
    subroutine writemywaves_linear(iproc,filename,iformat,npsidim,Lzd,orbs,nelec,at,rxyz,psi,nfvctr,coeff)
      use module_types
      use module_base
      use yaml_output
      use module_interfaces, except_this_one => writeonewave
      implicit none
      integer, intent(in) :: iproc,iformat,npsidim,nelec,nfvctr
      !integer, intent(in) :: norb   !< number of orbitals, not basis functions
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(in) :: orbs         !< orbs describing the basis functions
      type(local_zone_descriptors), intent(in) :: Lzd
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      real(wp), dimension(npsidim), intent(in) :: psi  ! Should be the real linear dimension and not the global
      real(wp), dimension(nfvctr,orbs%norb), intent(in) :: coeff
      character(len=*), intent(in) :: filename
      !Local variables
      logical :: binary, is_etsf
      integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr,shift,ii,iat,unitwf
      integer :: jorb,jlr
      real(kind=4) :: tr0,tr1
      real(kind=8) :: tel
    
      unitwf=99
      binary=(iformat/=WF_FORMAT_PLAIN)
      is_etsf=(iformat==WF_FORMAT_ETSF)
    
      if (iproc == 0) call yaml_map('Write wavefunctions to file', trim(filename)//'.*')
      !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
    
      !if (binary) then
      if (is_etsf) then
         call f_err_throw('Linear scaling with ETSF writing not implemented yet')
    !     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
      else
         call cpu_time(tr0)
         call system_clock(ncount1,ncount_rate,ncount_max)
    
         ! Write the TMBs in the Plain BigDFT files.
         ! Use same ordering as posinp and llr generation
         ii = 0
         do iat = 1, at%astruct%nat
            do iorb=1,orbs%norbp
               if(iat == orbs%onwhichatom(iorb+orbs%isorb)) then
                  shift = 1
                  do jorb = 1, iorb-1 
                     jlr = orbs%inwhichlocreg(jorb+orbs%isorb)
                     shift = shift + Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f
                  end do
                  ii = ii + 1
                  ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
                  do ispinor=1,orbs%nspinor
                     call open_filename_of_iorb(unitwf,binary,filename, &
                        & orbs,iorb,ispinor,iorb_out)
                     call writeonewave_linear(unitwf,.not. binary,iorb_out,&
                        & Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,&
                        & Lzd%Llr(ilr)%ns1,Lzd%Llr(ilr)%ns2,Lzd%Llr(ilr)%ns3,& 
                        & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                        & Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
                        & at%astruct%nat,rxyz,Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                        & Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                        & Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                        & Lzd%Llr(ilr)%wfd%keygloc(1:,Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & Lzd%Llr(ilr)%wfd%keyvloc(Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & psi(shift),psi(Lzd%Llr(ilr)%wfd%nvctr_c+shift),orbs%eval(iorb+orbs%isorb),&
                        & orbs%onwhichatom(iorb+orbs%isorb))
                     call f_close(unitwf)
                  end do
               end if
            enddo
         end do
    
        ! Now write the coefficients to file
        ! Must be careful, the orbs%norb is the number of basis functions
        ! while the norb is the number of orbitals.
        if(iproc == 0) then
           call f_open_file(unitwf,file=filename//'_coeff.bin',&
                binary=binary)
           !if(iformat == WF_FORMAT_PLAIN) then
           !  open(99, file=filename//'_coeff.bin', status='unknown',form='formatted')
           !else
           !open(99, file=filename//'_coeff.bin', status='unknown',form='unformatted')
           !end if
          call writeLinearCoefficients(unitwf,.not. binary,at%astruct%nat,rxyz,orbs%norb,&
               nelec,nfvctr,coeff,orbs%eval)
          call f_close(unitwf)
       end if
         call cpu_time(tr1)
         call system_clock(ncount2,ncount_rate,ncount_max)
         tel=dble(ncount2-ncount1)/dble(ncount_rate)
         if (iproc == 0) then
            call yaml_sequence_open('Write Waves Time')
            call yaml_sequence(advance='no')
            call yaml_mapping_open(flow=.true.)
            call yaml_map('Process',iproc)
            call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
            call yaml_mapping_close()
            call yaml_sequence_close()
         end if
         !write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
         !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
      end if
    
    END SUBROUTINE writemywaves_linear


    !> Write all my wavefunctions in files by calling writeonewave
    subroutine writemywaves_linear_fragments(iproc,filename,iformat,npsidim,Lzd,orbs,nelec,at,rxyz,psi,coeff,&
         dir_output,input_frag,ref_frags)
      use module_types
      use module_base
      use module_fragments
      use yaml_output
      use module_interfaces
      implicit none
      integer, intent(in) :: iproc,iformat,npsidim,nelec
      !integer, intent(in) :: norb   !< number of orbitals, not basis functions
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(in) :: orbs         !< orbs describing the basis functions
      type(local_zone_descriptors), intent(in) :: Lzd
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      real(wp), dimension(npsidim), intent(in) :: psi  ! Should be the real linear dimension and not the global
      real(wp), dimension(orbs%norb,orbs%norb), intent(in) :: coeff !SM: IS this correcy even with spin?
      character(len=*), intent(in) :: dir_output, filename
      type(fragmentInputParameters), intent(in) :: input_frag
      type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
      !Local variables
      integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr,shift,ii,iat
      integer :: jorb,jlr,isforb,isfat,ifrag,ifrag_ref,iforb,iiorb,iorbp,iiat,unitwf
      real(kind=4) :: tr0,tr1
      real(kind=8) :: tel
      character(len=256) :: full_filename
      logical, allocatable, dimension(:) :: fragment_written
    
      if (iproc == 0) call yaml_map('Write wavefunctions to file', trim(filename)//'.*')
      !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
    
      if (iformat == WF_FORMAT_ETSF) then
          stop 'Linear scaling with ETSF writing not implemented yet'
    !     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
      else
         call cpu_time(tr0)
         call system_clock(ncount1,ncount_rate,ncount_max)
    
         ! Write the TMBs in the Plain BigDFT files.
         ! For now only output one (first) set per reference fragment
        
         ! array to check if we already outputted tmbs for this fragment type
         fragment_written=f_malloc((/input_frag%nfrag_ref/),id='fragment_written')
         fragment_written=.false.
    
         unitwf=99
         isforb=0
         isfat=0
         loop_ifrag: do ifrag=1,input_frag%nfrag
            ! find reference fragment this corresponds to and check if we already outputted tmbs for this reference fragment
            ifrag_ref=input_frag%frag_index(ifrag)
            if (fragment_written(ifrag_ref)) then
               isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
               isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat   
               cycle
            end if
            fragment_written(ifrag_ref)=.true.
    
            ! loop over orbitals of this fragment
            loop_iforb: do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
               loop_iorb: do iorbp=1,orbs%norbp
                  iiorb=iorbp+orbs%isorb
    
                  ! check if this ref frag orbital corresponds to the orbital we want
                  if (iiorb/=iforb+isforb) cycle
    
                  ilr=orbs%inwhichlocreg(iiorb)
                  iiat=orbs%onwhichatom(iiorb)
    
                  shift = 1
                  do jorb = 1, iorbp-1 
                     jlr = orbs%inwhichlocreg(jorb+orbs%isorb)
                     shift = shift + Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f
                  end do
    
                  loop_ispinor: do ispinor=1,orbs%nspinor
                     ! as this is a fragment calculation frag%dirname should contain fragment directory (otherwise it would be empty - should add check)
                     ! bit of a hack to use orbs here not forbs, but different structures so this is necessary - to clean somehow
                     full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)
    
                     call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),full_filename, &
                          & orbs,iorbp,ispinor,iorb_out,iforb)
    
                     !also what to do with eval? - at the moment completely arbitrary
                     call writeonewave_linear(unitwf,(iformat == WF_FORMAT_PLAIN),iorb_out,&
                        & Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,&
                        & Lzd%Llr(ilr)%ns1,Lzd%Llr(ilr)%ns2,Lzd%Llr(ilr)%ns3,& 
                        & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                        & Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
                        & ref_frags(ifrag_ref)%astruct_frg%nat,rxyz(:,isfat+1:isfat+ref_frags(ifrag_ref)%astruct_frg%nat),&
                        & Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                        & Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                        & Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                        & Lzd%Llr(ilr)%wfd%keygloc(1:,Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & Lzd%Llr(ilr)%wfd%keyvloc(Lzd%Llr(ilr)%wfd%nseg_c+1:), &
                        & psi(shift),psi(Lzd%Llr(ilr)%wfd%nvctr_c+shift),-0.5d0, & !orbs%eval(iiorb),&
                        & orbs%onwhichatom(iiorb)-isfat)
    
                     close(unitwf)
    
                  end do loop_ispinor
               end do loop_iorb
            end do loop_iforb
    
    
            ! NEED to think about this - just make it diagonal for now? or random?  or truncate so they're not normalized?  or normalize after truncating?
            ! Or maybe don't write coeffs at all but assume we're always doing frag to frag and can use isolated frag coeffs?
    
            ! Now write the coefficients to file
            ! Must be careful, the orbs%norb is the number of basis functions
            ! while the norb is the number of orbitals.
            if(iproc == 0) then
               full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)
     
               call f_open_file(unitwf,file=trim(full_filename)//'_coeff.bin',&
                    binary=(iformat /= WF_FORMAT_PLAIN))
               !if(iformat == WF_FORMAT_PLAIN) then
               !   open(unitwf, file=trim(full_filename)//'_coeff.bin', status='unknown',form='formatted')
               !else
               !   open(unitwf, file=trim(full_filename)//'_coeff.bin', status='unknown',form='unformatted')
               !end if
               
               ! Not sure whether this is correct for nspin=2...
               call writeLinearCoefficients(unitwf,(iformat == WF_FORMAT_PLAIN),ref_frags(ifrag_ref)%astruct_frg%nat,&
                    rxyz(:,isfat+1:isfat+ref_frags(ifrag_ref)%astruct_frg%nat),ref_frags(ifrag_ref)%fbasis%forbs%norb,&
                    ref_frags(ifrag_ref)%nelec,&
                    ref_frags(ifrag_ref)%fbasis%forbs%norb, &
                    coeff(isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb,&
                    isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb),&
                    orbs%eval(isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb)) !-0.5d0
               call f_close(unitwf)
            end if
            call cpu_time(tr1)
            call system_clock(ncount2,ncount_rate,ncount_max)
            tel=dble(ncount2-ncount1)/dble(ncount_rate)
            if (iproc == 0) then
               call yaml_sequence_open('Write Waves Time')
               call yaml_sequence(advance='no')
               call yaml_mapping_open(flow=.true.)
               call yaml_map('Process',iproc)
               call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
               call yaml_mapping_close()
               call yaml_sequence_close()
            end if
    
            isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
            isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat    
         end do loop_ifrag
      end if
    
      call f_free(fragment_written)
    
    
    END SUBROUTINE writemywaves_linear_fragments




    subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,ns1,ns2,ns3,hx,hy,hz,locregCenter,&
         locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  &
         nseg_f,nvctr_f,keyg_f,keyv_f, &
         psi_c,psi_f,eval,onwhichatom)
      use module_base
      use yaml_output
      implicit none
      logical, intent(in) :: useFormattedOutput
      integer, intent(in) :: unitwf,iorb,n1,n2,n3,ns1,ns2,ns3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
      real(gp), intent(in) :: hx,hy,hz,locrad,confPotprefac
      real(wp), intent(in) :: eval
      integer, dimension(nseg_c), intent(in) :: keyv_c
      integer, dimension(nseg_f), intent(in) :: keyv_f
      integer, dimension(2,nseg_c), intent(in) :: keyg_c
      integer, dimension(2,nseg_f), intent(in) :: keyg_f
      real(wp), dimension(nvctr_c), intent(in) :: psi_c
      real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
      real(gp), dimension(3,nat), intent(in) :: rxyz
      real(gp), dimension(3), intent(in) :: locregCenter
      integer, intent(in) :: onwhichatom
      !local variables
      integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j,np,n1p1
      real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
    
      if (useFormattedOutput) then
         write(unitwf,*) iorb,eval
         write(unitwf,*) hx,hy,hz
         write(unitwf,*) n1,n2,n3
         write(unitwf,*) ns1,ns2,ns3
         write(unitwf,*) locregCenter(1),locregCenter(2),locregCenter(3),onwhichatom,locrad,&
              confPotOrder,confPotprefac
         write(unitwf,*) nat
         do iat=1,nat
            write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
         enddo
         write(unitwf,*) nvctr_c, nvctr_f
      else
         write(unitwf) iorb,eval
         write(unitwf) hx,hy,hz
         write(unitwf) n1,n2,n3
         write(unitwf) ns1,ns2,ns3
         write(unitwf) locregCenter(1),locregCenter(2),locregCenter(3),onwhichatom,locrad,&
              confPotOrder,confPotprefac
         write(unitwf) nat
         do iat=1,nat
            write(unitwf) (rxyz(j,iat),j=1,3)
         enddo
         write(unitwf) nvctr_c, nvctr_f
      end if
    
      n1p1=n1+1
      np=n1p1*(n2+1)
    
      ! coarse part
      do iseg=1,nseg_c
         jj=keyv_c(iseg)
         j0=keyg_c(1,iseg)
         j1=keyg_c(2,iseg)
         ii=j0-1
         i3=ii/np
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         do i=i0,i1
            tt=psi_c(i-i0+jj)
            if (useFormattedOutput) then
               write(unitwf,'(3(i4),1x,e19.12)') i,i2,i3,tt
            else
               write(unitwf) i,i2,i3,tt
            end if
         enddo
      enddo
    
      ! fine part
      do iseg=1,nseg_f
         jj=keyv_f(iseg)
         j0=keyg_f(1,iseg)
         j1=keyg_f(2,iseg)
         ii=j0-1
         i3=ii/np
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         do i=i0,i1
            t1=psi_f(1,i-i0+jj)
            t2=psi_f(2,i-i0+jj)
            t3=psi_f(3,i-i0+jj)
            t4=psi_f(4,i-i0+jj)
            t5=psi_f(5,i-i0+jj)
            t6=psi_f(6,i-i0+jj)
            t7=psi_f(7,i-i0+jj)
            if (useFormattedOutput) then
               write(unitwf,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
               write(unitwf) i,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
         enddo
      enddo
    
      if (verbose >= 2 .and. bigdft_mpi%iproc==0) call yaml_map('Wavefunction written No.',iorb)
      !if (verbose >= 2) write(*,'(1x,i0,a)') iorb,'th wavefunction written'
    
    END SUBROUTINE writeonewave_linear


    subroutine read_linear_matrix_dense(iunit, ntmb, nat, matrix, rxyz, on_which_atom)
      use module_base
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iunit, ntmb, nat
      real(kind=8),dimension(ntmb,ntmb),intent(out) :: matrix
      real(kind=8),dimension(nat,nat),intent(out),optional :: rxyz
      integer,dimension(ntmb),intent(out),optional :: on_which_atom
    
      ! Local variables
      integer :: itmb, jtmb, ii, jj, iat
      logical :: read_rxyz, read_on_which_atom
      real(kind=8),dimension(3) :: dummy
    
      read_on_which_atom = present(on_which_atom)
      read_rxyz = present(rxyz)

      do iat=1,nat
          if (read_rxyz) then
              read(iunit,*) rxyz(1:3,iat)
          else
              read(iunit,*) dummy(1:3)
          end if
      end do
    
      do itmb=1,ntmb
          do jtmb=1,ntmb
              if(read_on_which_atom .and. jtmb==1) then
                  read(iunit,*) ii, jj, matrix(ii,jj), on_which_atom(itmb)
              else
                  read(iunit,*) ii, jj, matrix(ii,jj)
              end if
              if (ii/=itmb) call f_err_throw('ii/=itmb',err_name='BIGDFT_RUNTIME_ERROR')
              if (jj/=jtmb) call f_err_throw('jj/=jtmb',err_name='BIGDFT_RUNTIME_ERROR')
          end do
      end do
    
    end subroutine read_linear_matrix_dense


    subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, &
           & lstat, error, nat, rxyz_old)
        use module_base
        use module_types
        !use internal_io
        implicit none
        integer, intent(in) :: unitwf
        logical, intent(in) :: formatted
        integer, intent(out) :: norb_old, ntmb_old
        logical, intent(out) :: lstat
        character(len =256), intent(out) :: error
        ! Optional arguments
        integer, intent(in), optional :: nat
        real(gp), dimension(:,:), intent(out), optional :: rxyz_old
    
        integer :: i, iat, i_stat, nat_
        real(gp) :: rxyz(3)
    
        lstat = .false.
        write(error, "(A)") "cannot read coeff description."
        if (formatted) then
           read(unitwf,*,iostat=i_stat) ntmb_old, norb_old
           if (i_stat /= 0) return
           !write(*,*) 'reading ',nat,' atomic positions'
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
              enddo
           else
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,*,iostat=i_stat)
                 if (i_stat /= 0) return
              enddo
           end if
           !read(unitwf,*,iostat=i_stat) i, iat
           !if (i_stat /= 0) return
        else
           read(unitwf,iostat=i_stat) ntmb_old, norb_old
           if (i_stat /= 0) return
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
              enddo
           else
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,iostat=i_stat) rxyz
                 if (i_stat /= 0) return
              enddo
           end if
           !read(unitwf,iostat=i_stat) i, iat
           !if (i_stat /= 0) return
        end if
        lstat = .true.
    END SUBROUTINE io_read_descr_coeff


    subroutine read_coeff_minbasis(unitwf,useFormattedInput,iproc,ntmb,norb_old,nfvctr,coeff,eval,nat,rxyz_old)
      use module_base
      use module_types
      !use internal_io
      use yaml_output
      implicit none
      logical, intent(in) :: useFormattedInput
      integer, intent(in) :: unitwf,iproc,ntmb,nfvctr
      integer, intent(out) :: norb_old
      real(wp), dimension(nfvctr,ntmb), intent(out) :: coeff
      real(wp), dimension(ntmb), intent(out) :: eval
      integer, optional, intent(in) :: nat
      real(gp), dimension(:,:), optional, intent(out) :: rxyz_old
    
      !local variables
      character(len = 256) :: error
      logical :: lstat
      integer :: i_stat
      integer :: ntmb_old, i1, i2,i,j,iorb,iorb_old
      real(wp) :: tt
    
      call io_read_descr_coeff(unitwf, useFormattedInput, norb_old, ntmb_old, &
           & lstat, error, nat, rxyz_old)
      if (.not. lstat) call io_error(trim(error))
    
      if (ntmb_old /= ntmb_old) then
         if (iproc == 0) write(error,"(A)") 'error in read coeffs, ntmb_old/=ntmb'
         call io_error(trim(error))
      end if
    
      ! read the eigenvalues
      if (useFormattedInput) then
         do iorb=1,ntmb
            read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
            if (iorb_old /= iorb) stop 'read_coeff_minbasis'
         enddo
      else 
         do iorb=1,ntmb
            read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
            if (iorb_old /= iorb) stop 'read_coeff_minbasis'
         enddo
         if (i_stat /= 0) stop 'Problem reading the eigenvalues'
      end if
    
      !if (iproc == 0) write(*,*) 'coefficients need NO reformatting'
    
      ! Now read the coefficients
      do i = 1, ntmb
         do j = 1,nfvctr
            if (useFormattedInput) then
               read(unitwf,*,iostat=i_stat) i1,i2,tt
            else
               read(unitwf,iostat=i_stat) i1,i2,tt
            end if
            if (i_stat /= 0) stop 'Problem reading the coefficients'
            coeff(j,i) = tt  
         end do
      end do
    
      ! rescale so first significant element is +ve
      do i = 1, ntmb
         do j = 1,nfvctr
            if (abs(coeff(j,i))>1.0e-1) then
               if (coeff(j,i)<0.0_gp) call dscal(ntmb,-1.0_gp,coeff(1,i),1)
               exit
            end if
         end do
         if (j==ntmb+1) print*,'Error finding significant coefficient, coefficients not scaled to have +ve first element'
      end do
    
    END SUBROUTINE read_coeff_minbasis


    subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n_old1, n_old2, n_old3, &
           & ns_old1, ns_old2, ns_old3, hgrids_old, lstat, error, onwhichatom, locrad, locregCenter, &
           & confPotOrder, confPotprefac, nvctr_c_old, nvctr_f_old, nat, rxyz_old)
        use module_base
        use module_types
        !use internal_io
        use yaml_output
        implicit none
    
        integer, intent(in) :: unitwf
        logical, intent(in) :: formatted
        integer, intent(out) :: iorb_old
        integer, intent(out) :: n_old1, n_old2, n_old3, ns_old1, ns_old2, ns_old3
        real(gp), dimension(3), intent(out) :: hgrids_old
        logical, intent(out) :: lstat
        real(wp), intent(out) :: eval
        real(gp), intent(out) :: locrad
        real(gp), dimension(3), intent(out) :: locregCenter
        character(len =256), intent(out) :: error
        integer, intent(out) :: onwhichatom
        integer, intent(out) :: confPotOrder
        real(gp), intent(out) :: confPotprefac
        ! Optional arguments
        integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
        integer, intent(in), optional :: nat
        real(gp), dimension(:,:), intent(out), optional :: rxyz_old
    
        integer :: i, iat, i_stat, nat_
        real(gp) :: rxyz(3)
    
        lstat = .false.
        write(error, "(A)") "cannot read psi description."
        if (formatted) then
           read(unitwf,*,iostat=i_stat) iorb_old,eval
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) hgrids_old(1),hgrids_old(2),hgrids_old(3)
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) n_old1,n_old2,n_old3
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) ns_old1,ns_old2,ns_old3
           if (i_stat /= 0) return
           read(unitwf,*,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
                locrad,confPotOrder, confPotprefac
           if (i_stat /= 0) return
           !call yaml_map('Reading atomic positions',nat)
           !write(*,*) 'reading ',nat,' atomic positions' !*
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
    
              enddo
           else
              read(unitwf,*,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,*,iostat=i_stat)
                 if (i_stat /= 0) return
              enddo
           end if
           if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
              read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
              if (i_stat /= 0) return
           else
              read(unitwf,*,iostat=i_stat) i, iat
              if (i_stat /= 0) return
           end if
        else
           read(unitwf,iostat=i_stat) iorb_old,eval
           if (i_stat /= 0) return
    
           read(unitwf,iostat=i_stat) hgrids_old(1),hgrids_old(2),hgrids_old(3)
           if (i_stat /= 0) return
           read(unitwf,iostat=i_stat) n_old1,n_old2,n_old3
           if (i_stat /= 0) return
           read(unitwf,iostat=i_stat) ns_old1,ns_old2,ns_old3
           if (i_stat /= 0) return
           read(unitwf,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
                locrad,confPotOrder, confPotprefac
           if (i_stat /= 0) return
           if (present(nat) .And. present(rxyz_old)) then
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              ! Sanity check
              if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
              if (nat_ /= nat) stop "Mismatch in coordinate array size."
              do iat=1,nat
                 read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
                 if (i_stat /= 0) return
              enddo
           else
              read(unitwf,iostat=i_stat) nat_
              if (i_stat /= 0) return
              do iat=1,nat_
                 read(unitwf,iostat=i_stat) rxyz
                 if (i_stat /= 0) return
              enddo
           end if
           if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
              read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
              if (i_stat /= 0) return
           else
              read(unitwf,iostat=i_stat) i, iat
              if (i_stat /= 0) return
           end if
        end if
        lstat = .true.
    
    END SUBROUTINE io_read_descr_linear


    subroutine io_error(error)
      use module_defs
  
      implicit none
  
      character(len = *), intent(in) :: error
      integer :: ierr
  
      call io_warning(error)
      call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
    END SUBROUTINE io_error
  
  
    subroutine io_warning(error)
      use module_defs
  
      implicit none
  
      character(len = *), intent(in) :: error
  
      write(0,"(2A)") "WARNING! ", trim(error)
    END SUBROUTINE io_warning
  
  
    !> Read the input/output descriptors (for a wavefunction for instance)
    subroutine io_read_descr(unitwf, formatted, iorb_old, eval, n1_old, n2_old, n3_old, &
         & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
      use module_base
      use module_types
  
      implicit none
  
      integer, intent(in) :: unitwf
      logical, intent(in) :: formatted
      integer, intent(out) :: iorb_old
      integer, intent(out) :: n1_old, n2_old, n3_old
      real(gp), intent(out) :: hx_old, hy_old, hz_old
      logical, intent(out) :: lstat
      real(wp), intent(out) :: eval
      character(len =256), intent(out) :: error
      ! Optional arguments
      integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
      integer, intent(in), optional :: nat
      real(gp), dimension(:,:), intent(out), optional :: rxyz_old
  
      integer :: i, iat, i_stat, nat_
      real(gp) :: rxyz(3)
  
      lstat = .false.
      write(error, "(A)") "cannot read psi description."
      if (formatted) then
         read(unitwf,*,iostat=i_stat) iorb_old,eval
         if (i_stat /= 0) return
         read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
         if (i_stat /= 0) return
         read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
         if (i_stat /= 0) return
         !write(*,*) 'reading ',nat,' atomic positions'
         if (present(nat) .And. present(rxyz_old)) then
            read(unitwf,*,iostat=i_stat) nat_
            if (i_stat /= 0) return
            ! Sanity check
            if (size(rxyz_old, 2) /= nat) call io_error("Mismatch in coordinate array size.")
            if (nat_ /= nat) call io_error("Mismatch in coordinate array size.")
            do iat=1,nat
               read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
               if (i_stat /= 0) return
            enddo
         else
            read(unitwf,*,iostat=i_stat) nat_
            if (i_stat /= 0) return
            do iat=1,nat_
               read(unitwf,*,iostat=i_stat)
               if (i_stat /= 0) return
            enddo
         end if
         if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
            read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
            if (i_stat /= 0) return
         else
            read(unitwf,*,iostat=i_stat) i, iat
            if (i_stat /= 0) return
         end if
      else
         read(unitwf,iostat=i_stat) iorb_old,eval
         if (i_stat /= 0) return
         read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
         if (i_stat /= 0) return
         read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
         if (i_stat /= 0) return
         if (present(nat) .And. present(rxyz_old)) then
            read(unitwf,iostat=i_stat) nat_
            if (i_stat /= 0) return
            ! Sanity check
            if (size(rxyz_old, 2) /= nat) call io_error("Mismatch in coordinate array size.")
            if (nat_ /= nat) call io_error("Mismatch in coordinate array size.")
            do iat=1,nat
               read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
               if (i_stat /= 0) return
            enddo
         else
            read(unitwf,iostat=i_stat) nat_
            if (i_stat /= 0) return
            do iat=1,nat_
               read(unitwf,iostat=i_stat) rxyz
               if (i_stat /= 0) return
            enddo
         end if
         if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
            read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
            if (i_stat /= 0) return
         else
            read(unitwf,iostat=i_stat) i, iat
            if (i_stat /= 0) return
         end if
      end if
      lstat = .true.
    END SUBROUTINE io_read_descr
  
  
    subroutine io_gcoordToLocreg(n1, n2, n3, nvctr_c, nvctr_f, gcoord_c, gcoord_f, lr)
      use module_base
      use locregs
  
      implicit none
      !Arguments
      integer, intent(in) :: n1, n2, n3, nvctr_c, nvctr_f
      integer, dimension(3, nvctr_c), intent(in) :: gcoord_c
      integer, dimension(3, nvctr_f), intent(in) :: gcoord_f
      type(locreg_descriptors), intent(out) :: lr
      !Local variables
      character(len = *), parameter :: subname = "io_gcoordToLocreg"
      integer :: i
      logical, dimension(:,:,:), allocatable :: logrid_c, logrid_f
  
      call f_routine(id=subname)
  
      call nullify_locreg_descriptors(lr)
  
      lr%geocode = "P"
      lr%hybrid_on = .false.
  
      lr%ns1 = 0
      lr%ns2 = 0
      lr%ns3 = 0
  
      lr%d%n1 = n1
      lr%d%n2 = n2
      lr%d%n3 = n3
  
      lr%d%n1i = 2 * n1 + 2
      lr%d%n2i = 2 * n2 + 2
      lr%d%n3i = 2 * n3 + 2
  
      logrid_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_c')
      logrid_f = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_f')
  
      lr%d%nfl1 = n1
      lr%d%nfl2 = n2
      lr%d%nfl3 = n3
      lr%d%nfu1 = 0
      lr%d%nfu2 = 0
      lr%d%nfu3 = 0
  
      logrid_c(:,:,:) = .false.
      do i = 1, nvctr_c, 1
         logrid_c(gcoord_c(1, i), gcoord_c(2, i), gcoord_c(3, i)) = .true.
      end do
      logrid_f(:,:,:) = .false.
      do i = 1, nvctr_f, 1
         logrid_f(gcoord_f(1, i), gcoord_f(2, i), gcoord_f(3, i)) = .true.
         lr%d%nfl1 = min(lr%d%nfl1, gcoord_f(1, i))
         lr%d%nfl2 = min(lr%d%nfl2, gcoord_f(2, i))
         lr%d%nfl3 = min(lr%d%nfl3, gcoord_f(3, i))
         lr%d%nfu1 = max(lr%d%nfu1, gcoord_f(1, i))
         lr%d%nfu2 = max(lr%d%nfu2, gcoord_f(2, i))
         lr%d%nfu3 = max(lr%d%nfu3, gcoord_f(3, i))
      end do
  
      !correct the values of the delimiter if there are no wavelets
      if (lr%d%nfl1 == n1 .and. lr%d%nfu1 == 0) then
         lr%d%nfl1 = n1 / 2
         lr%d%nfu1 = n1 / 2
      end if
      if (lr%d%nfl2 == n2 .and. lr%d%nfu2 == 0) then
         lr%d%nfl2 = n2 / 2
         lr%d%nfu2 = n2 / 2
      end if
      if (lr%d%nfl3 == n3 .and. lr%d%nfu3 == 0) then
         lr%d%nfl3 = n3 / 2
         lr%d%nfu3 = n3 / 2
      end if
  
      call wfd_from_grids(logrid_c, logrid_f, .true., lr)
  
      call f_free(logrid_c)
      call f_free(logrid_f)
  
      call f_release_routine()
  
    END SUBROUTINE io_gcoordToLocreg
  
    subroutine read_psi_compress(unitwf, formatted, nvctr_c, nvctr_f, psi, lstat, error, gcoord_c, gcoord_f)
      use module_base
      use module_types
  
      implicit none
  
      integer, intent(in) :: unitwf, nvctr_c, nvctr_f
      logical, intent(in) :: formatted
      real(wp), dimension(nvctr_c+7*nvctr_f), intent(out) :: psi
      logical, intent(out) :: lstat
      character(len =256), intent(out) :: error
      integer, dimension(3, nvctr_c), optional, intent(out) :: gcoord_c
      integer, dimension(3, nvctr_f), optional, intent(out) :: gcoord_f
  
      integer :: j, i1, i2, i3, i_stat
      real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
  
      lstat = .false.
      write(error, "(A)") "cannot read psi values."
      if (present(gcoord_c)) then
         do j=1,nvctr_c
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,tt
            end if
            if (i_stat /= 0) return
            psi(j)=tt
            gcoord_c(:, j) = (/ i1, i2, i3 /)
         enddo
      else
         do j=1,nvctr_c
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,tt
            end if
            if (i_stat /= 0) return
            psi(j)=tt
         enddo
      end if
      if (present(gcoord_f)) then
         do j=1,7*nvctr_f-6,7
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
            if (i_stat /= 0) return
            psi(nvctr_c+j+0)=t1
            psi(nvctr_c+j+1)=t2
            psi(nvctr_c+j+2)=t3
            psi(nvctr_c+j+3)=t4
            psi(nvctr_c+j+4)=t5
            psi(nvctr_c+j+5)=t6
            psi(nvctr_c+j+6)=t7
            gcoord_f(:, (j - 1) / 7 + 1) = (/ i1, i2, i3 /)
         enddo
      else
         do j=1,7*nvctr_f-6,7
            if (formatted) then
               read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            else
               read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            end if
            if (i_stat /= 0) return
            psi(nvctr_c+j+0)=t1
            psi(nvctr_c+j+1)=t2
            psi(nvctr_c+j+2)=t3
            psi(nvctr_c+j+3)=t4
            psi(nvctr_c+j+4)=t5
            psi(nvctr_c+j+5)=t6
            psi(nvctr_c+j+6)=t7
         enddo
      end if
      lstat = .true.
    END SUBROUTINE read_psi_compress
  
  
    subroutine read_psig(unitwf, formatted, nvctr_c, nvctr_f, n1, n2, n3, psig, lstat, error)
      use module_base
      use module_types
  
      implicit none
  
      integer, intent(in) :: unitwf, nvctr_c, nvctr_f, n1, n2, n3
      logical, intent(in) :: formatted
      real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psig
      logical, intent(out) :: lstat
      character(len =256), intent(out) :: error
  
      integer :: i1, i2, i3, i_stat, iel
      real(wp) :: tt, t1, t2, t3, t4, t5, t6, t7
  
      lstat = .false.
      write(error, "(A)") "cannot read psig values."
  
      call f_zero(psig)
      do iel=1,nvctr_c
         if (formatted) then
            read(unitwf,*,iostat=i_stat) i1,i2,i3,tt
         else
            read(unitwf,iostat=i_stat) i1,i2,i3,tt
         end if
         if (i_stat /= 0) return
         psig(i1,1,i2,1,i3,1)=tt
      enddo
      do iel=1,nvctr_f
         if (formatted) then
            read(unitwf,*,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         else
            read(unitwf,iostat=i_stat) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
         end if
         if (i_stat /= 0) return
         psig(i1,2,i2,1,i3,1)=t1
         psig(i1,1,i2,2,i3,1)=t2
         psig(i1,2,i2,2,i3,1)=t3
         psig(i1,1,i2,1,i3,2)=t4
         psig(i1,2,i2,1,i3,2)=t5
         psig(i1,1,i2,2,i3,2)=t6
         psig(i1,2,i2,2,i3,2)=t7
      enddo
      lstat = .true.
    END SUBROUTINE read_psig
  
    subroutine io_open(unitwf, filename, formatted)
      use f_utils, only: f_open_file
      implicit none
      character(len = *), intent(in) :: filename
      logical, intent(in) :: formatted
      integer, intent(out) :: unitwf
  
      integer :: i_stat
  
      ! We open the Fortran file
      unitwf = 99
      call f_open_file(unitwf,file=trim(filename),binary=.not. formatted)
  !!$    if (.not. formatted) then
  !!$       open(unit=unitwf,file=trim(filename),status='unknown',form="unformatted", iostat=i_stat)
  !!$    else
  !!$       open(unit=unitwf,file=trim(filename),status='unknown', iostat=i_stat)
  !!$    end if
  !!$    if (i_stat /= 0) then
  !!$       call io_warning("Cannot open file '" // trim(filename) // "'.")
  !!$       unitwf = -1
  !!$       return
  !!$    end if
    END SUBROUTINE io_open


end module io
