
      !######################################################################################################################################################################
      ! CUBE TREATMENT BEGINS HERE
      ! TO DO: Redefine the timing
      !        Do the Parallelization
      !######################################################################################################################################################################

!!$   else if ( (filetype == 'cube' .or. filetype == 'CUBE') .and. nproc==1 ) then
!!$
!!$
!!$      ! Read integers in order to allocate tables used to store
!!$      ! wavefunctions, l, mr, rvalue, zona, ...
!!$
!!$      call read_cube_header_1(nx, ny, nz, n_at, bx, by, bz)
!!$      allocate(Z(n_at)) 
!!$      allocate(at_pos(n_at,3))
!!$      call read_cube_header_2(n_at, Z, at_pos)
!!$      call read_nnkp_int_alloc(iproc,seedname, n_kpts, n_proj, n_nnkpts, n_excb)
!!$
!!$
!!$      !! Find if the values read do not reach limitations of the program
!!$      !! Stops the program and send an error message in case there are problems
!!$      !
!!$      !   call limitations(seedname, n_proj, n_occ, n_virt_tot, n_bands, n_kpts)
!!$
!!$
!!$      ! Allocations
!!$
!!$      allocate(kpts(n_kpts,3))
!!$      allocate(ctr_proj(n_proj,3))
!!$      allocate(x_proj(n_proj,3))
!!$      allocate(y_proj(n_proj,3))
!!$      allocate(z_proj(n_proj,3))
!!$      allocate(l(n_proj))
!!$      allocate(mr(n_proj))
!!$      allocate(rvalue(n_proj))
!!$      allocate(zona(n_proj))
!!$      allocate(k_plus_b(n_kpts*n_nnkpts,2))
!!$      allocate(G_vec(n_kpts*n_nnkpts,3))
!!$      allocate(excb(n_excb))
!!$
!!$
!!$      ! Read Wannier90 .nnkp file.
!!$      ! The most important informations to be read are : 
!!$      !  - ctr_proj : the coordinates of the center of the projections
!!$      !  - l, mr, rvalue and zona = Z/a : the parameters used to build spherical harmonics
!!$      !  - k_plus_b and G_vec : the parameters used to build the nearest neighbours k-points
!!$
!!$      call read_nnkp(iproc,seedname, calc_only_A, real_latt, recip_latt, n_kpts, n_proj, n_nnkpts, &
!!$         &   n_excb, kpts, ctr_proj, x_proj, z_proj, l, mr, rvalue, &
!!$         &   zona, k_plus_b, G_vec, excb)
!!$
!!$      ! Check that z_proj and x_proj are orthonormal 
!!$      do np = 1, n_proj
!!$         znorm = z_proj(np,1)**2+z_proj(np,2)**2+z_proj(np,3)**2
!!$         xnorm = x_proj(np,1)**2+x_proj(np,2)**2+x_proj(np,3)**2
!!$         ortho = z_proj(np,1)*x_proj(np,1) + z_proj(np,2)*x_proj(np,2) + z_proj(np,3)*x_proj(np,3)
!!$         if(abs(znorm - 1.d0) > eps8 .or. abs(znorm - 1.d0) > eps8 .or. abs(ortho) > eps8) then
!!$            if(iproc == 0) then
!!$               write(*,'(A)') 'Checkorthonormality of z_proj and x_proj:'
!!$               write(*,'(A,e3.2)') 'z norm: ',znorm
!!$               write(*,'(A,e3.2)') 'x norm: ',xnorm
!!$               write(*,'(A,e3.2)') 'x dot z: ',ortho
!!$            end if
!!$         end if
!!$         !Now we need to calculate the y direction
!!$         y_proj(np,1) = x_proj(np,3)*z_proj(np,2) - x_proj(np,2)*z_proj(np,3)
!!$         y_proj(np,2) = x_proj(np,1)*z_proj(np,3) - x_proj(np,3)*z_proj(np,1)
!!$         y_proj(np,3) = x_proj(np,2)*z_proj(np,1) - x_proj(np,1)*z_proj(np,2)
!!$      end do
!!$
!!$      if (pre_check .eqv. .true.) then 
!!$         ! If pre-check mode is chosen, calculation of the the scalar product of all unoccupied wavefunctions calculated by BigDFT and spherical harmonics.
!!$
!!$         ! Define which unoccupied orbitals have to be read in pre-check mode.
!!$         ! virt_list nominates the number of all the virtual orbitals up to n_virt_tot (integer defined in input.inter)
!!$         n_bands=n_occ+n_virt_tot
!!$         if (n_virt_tot .ne. 0) then
!!$            allocate(virt_list(n_virt_tot))
!!$            do i=1,n_virt_tot
!!$               virt_list(i)=i
!!$            end do
!!$            if (n_virt == 0) then
!!$               print *, 'There is no need to do the pre-check if there are no unoccupied states to add to the localization process'
!!$               STOP
!!$            end if
!!$         end if
!!$
!!$         ! Algorithm to compute the scalar product :
!!$         ! The term 'sqrt(bx(1)*by(2)*bz(3))' is there to normalize spherical harmonics.
!!$         ! Wavefunctions calculated by BigDFT already are normalized.
!!$         write(*,*) '!==================================!'
!!$         write(*,*) '! Calculating amnk=<virt|sph_har>  !'
!!$         write(*,*) '!       in pre-check mode :        !'
!!$         write(*,*) '!==================================!'
!!$         write(*,'(A12,4x,A15)') 'Virtual band', 'amnk_guess(nb)='
!!$         ! b1, b2 and b3 are the norm of the lattice parameters
!!$         ! bx(:), by(:) and bz(:) define the basis used in BigDFT calculation (in Cartesian)
!!$         ! nx, ny, and nz define the number of points where wavefunctions are calculated, following each cartesian axis
!!$         allocate(amnk(n_virt_tot,n_proj))
!!$         allocate(amnk_guess(n_virt_tot))
!!$         allocate(amnk_guess_sorted(n_virt))
!!$         allocate(amnk_bands_sorted(n_virt))
!!$         b1=sqrt(nx*bx(1)*nx*bx(1)+ny*bx(2)*ny*bx(2)+nz*bx(3)*nz*bx(3))
!!$         b2=sqrt(nx*by(1)*nx*by(1)+ny*by(2)*ny*by(2)+nz*by(3)*nz*by(3))
!!$         b3=sqrt(nx*bz(1)*nx*bz(1)+ny*bz(2)*ny*bz(2)+nz*bz(3)*nz*bz(3))
!!$         do nb=1,n_virt_tot
!!$            amnk_guess(nb)=0.d0
!!$            allocate(virt(nx,ny,nz))
!!$            call read_virtcube_1(nx, ny, nz, n_at, virt, nb+n_occ, virt_list, n_virt_tot, n_occ)
!!$            do np=1, n_proj
!!$               amnk(nb,np)=0.d0
!!$               r0x=ctr_proj(np,1)*b1
!!$               r0y=ctr_proj(np,2)*b2
!!$               r0z=ctr_proj(np,3)*b3
!!$               do k=1,nz
!!$                  if ( nb == 1 ) then
!!$                     zz=k*bz(3)-r0z
!!$                  end if
!!$                  do j=1,ny
!!$                     if ( nb == 1 ) then
!!$                        yy=j*by(2)-r0y
!!$                     end if
!!$                     do i=1,nx
!!$                        if ( nb == 1 ) then
!!$                           if ( (np==1) .and. (i==1) .and. (j==1) .and. (k==1) ) then
!!$                              allocate(ylm(nx,ny,nz))
!!$                              allocate(func_r(nx,ny,nz))
!!$                              allocate(sph_har(nx,ny,nz,n_proj))
!!$                           end if
!!$                           xx=i*bx(1)-r0x
!!$                           call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
!!$                              &   xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
!!$                           call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
!!$                              &   xx, yy, zz, n_proj, func_r)
!!$                           sph_har(i,j,k,np)=func_r(i,j,k)*ylm(i,j,k)
!!$                           if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then
!!$                              call write_cube(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, n_proj, &
!!$                                 &   nx, ny, nz, n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
!!$                              if ( np==n_proj) then
!!$                                 deallocate(func_r)
!!$                                 deallocate(ylm)
!!$                              end if
!!$                           end if
!!$                        end if
!!$                        amnk(nb,np)=amnk(nb,np)+virt(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
!!$                     end do
!!$                  end do
!!$               end do
!!$               ! sqrt(amnk_tot) verifies the normalization in each band (it must tend to 1)
!!$               amnk_guess(nb)=amnk_guess(nb)+(amnk(nb,np))**2
!!$            end do
!!$            deallocate(virt)
!!$            write(*,'(I4,11x,F12.6)') nb, sqrt(amnk_guess(nb))
!!$         end do
!!$         deallocate(virt_list)
!!$         write(*,*)
!!$         write(*,'(1a)') 'These are the virtual bands to use to construct the actual Amn and Mmn matrices :'
!!$         write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_guess)='
!!$         do nb1=1,n_virt
!!$            amnk_guess_sorted(nb1)=maxval(amnk_guess,n_virt_tot)
!!$            amnk_bands_sorted(nb1)=maxloc(amnk_guess,n_virt_tot)
!!$            amnk_guess(amnk_bands_sorted(nb1))=0.d0
!!$            write(*,'(I4,3x,F12.6)') amnk_bands_sorted(nb1), sqrt(amnk_guess_sorted(nb1))
!!$         end do
!!$         deallocate(amnk)
!!$         deallocate(amnk_guess)
!!$         deallocate(amnk_guess_sorted)
!!$         write(*,*) '!==================================!'
!!$         write(*,*) '! Calculating amnk=<virt|sph_har>  !'
!!$         write(*,*) '!     in pre-check mode done       !'
!!$         write(*,*) '!==================================!'
!!$         write(*,*)
!!$         write(*,*)
!!$
!!$         ! Rewrite the input.inter file to add the unoccupied states at the end and put w_rad, w_sph and w_ang to false.
!!$         call write_inter(n_virt, amnk_bands_sorted)
!!$
!!$      end if
!!$
!!$      ! In this part of the program, the pre-check has already been done.
!!$      ! Calculate the scalar product of chosen wavefunctions calculated by BigDFT and spherical harmonics computed above.
!!$      ! This values are written in a .amn file that will be then used by Wannier90.
!!$
!!$      ! Define which unoccupied orbitals have to be read.
!!$      ! Here, virt_list nominates the virtual orbitals chosen in pre-check mode.
!!$      n_bands=n_occ+n_virt
!!$      if (n_virt .ne. 0) then
!!$         allocate(virt_list(n_virt))
!!$         if (pre_check .eqv. .true.) virt_list(:)=amnk_bands_sorted(:)
!!$         if (pre_check .eqv. .false.) call read_inter_list(iproc, n_virt, virt_list)
!!$         !   else 
!!$         !      STOP
!!$      end if
!!$
!!$      call timing(iproc,'CrtProjectors ','ON')
!!$
!!$      ! Quick algorithm to compute the amnk matrix :
!!$      ! The term 'sqrt(bx(1)*by(2)*bz(3))' is there to normalize spherical harmonics.
!!$      ! Wavefunctions already are normalized.
!!$      write(*,*) '!==================================!'
!!$      write(*,*) '! Calculating amnk=<psi|sph_har> : !'
!!$      write(*,*) '!==================================!'
!!$      write(*,*) 'The values of sqrt(amnk_tot) check the normalization in each band.'
!!$      write(*,*) 'They must tend to their lower limit value, which is equal to 1 :'
!!$      write(*,'(A4,4x,A15)') 'Band', 'sqrt(amnk_tot)='
!!$      ! b1, b2 and b3 are the norm of the lattice parameters
!!$      ! bx(:), by(:) and bz(:) define the basis used in BigDFT calculation (in Cartesian)
!!$      ! nx, ny, and nz define the number of points where wavefunctions are calculated, following each cartesian axis
!!$      allocate(amnk(n_bands,n_proj))
!!$      if (pre_check .eqv. .false.) allocate(amnk_bands_sorted(n_virt))
!!$      allocate(amnk_tot(n_bands))
!!$      b1=sqrt(nx*bx(1)*nx*bx(1)+ny*bx(2)*ny*bx(2)+nz*bx(3)*nz*bx(3))
!!$      b2=sqrt(nx*by(1)*nx*by(1)+ny*by(2)*ny*by(2)+nz*by(3)*nz*by(3))
!!$      b3=sqrt(nx*bz(1)*nx*bz(1)+ny*bz(2)*ny*bz(2)+nz*bz(3)*nz*bz(3))
!!$      if (n_virt == n_occ) then
!!$         do nb=1,n_occ
!!$            amnk_tot(nb)=0.d0
!!$            amnk_tot(nb+n_occ)=0.d0
!!$            allocate(orb(nx,ny,nz))
!!$            call read_orbcube_1(nx, ny, nz, n_at, orb, nb) 
!!$            allocate(virt(nx,ny,nz))
!!$            call read_virtcube_1(nx, ny, nz, n_at, virt, nb+n_occ, virt_list, n_virt, n_occ)
!!$            do np=1,n_proj
!!$               amnk(nb,np)=0.d0
!!$               amnk(nb+n_occ,np)=0.d0
!!$               r0x=ctr_proj(np,1)*b1
!!$               r0y=ctr_proj(np,2)*b2
!!$               r0z=ctr_proj(np,3)*b3
!!$               do k=1,nz
!!$                  if ( nb == 1 ) zz=k*bz(3)-r0z
!!$                  do j=1,ny
!!$                     if ( nb == 1 ) yy=j*by(2)-r0y
!!$                     do i=1,nx
!!$                        if ( nb == 1 .and. (pre_check .eqv. .false.) ) then
!!$                           if ( (np==1) .and. (i==1) .and. (j==1) .and. (k==1) ) then
!!$                              allocate(ylm(nx,ny,nz))
!!$                              allocate(func_r(nx,ny,nz))
!!$                              allocate(sph_har(nx,ny,nz,n_proj))
!!$                           end if
!!$                           xx=i*bx(1)-r0x
!!$                           call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
!!$                              &   xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
!!$                           call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
!!$                              &   xx, yy, zz, n_proj, func_r)
!!$                           sph_har(i,j,k,np)=func_r(i,j,k)*ylm(i,j,k)
!!$                           if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then
!!$                              call write_cube(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, n_proj, &
!!$                                 &   nx, ny, nz, n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
!!$                              if ( np==n_proj ) then
!!$                                 deallocate(func_r)
!!$                                 deallocate(ylm)
!!$                              end if
!!$                           end if
!!$                        end if
!!$                        amnk(nb,np)=amnk(nb,np)+orb(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
!!$                        amnk(nb+n_occ,np)=amnk(nb+n_occ,np)+virt(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
!!$                     end do
!!$                  end do
!!$               end do
!!$               ! sqrt(amnk_tot) verifies the normalization in each band (it must tend to 1)
!!$               amnk_tot(nb)=amnk_tot(nb)+(amnk(nb,np))**2
!!$               amnk_tot(nb+n_occ)=amnk_tot(nb+n_occ)+(amnk(nb+n_occ,np))**2
!!$            end do
!!$            if (allocated(orb)) deallocate(orb)
!!$            if (allocated(virt)) deallocate(virt)
!!$         end do
!!$         do nb=1,n_bands
!!$            write(*,'(I4,3x,F12.6)') nb, sqrt(amnk_tot(nb))
!!$         end do
!!$         deallocate(amnk_tot)
!!$         deallocate(sph_har)
!!$      else
!!$         do nb=1,n_bands
!!$            amnk_tot(nb)=0.d0
!!$            if ( (nb .gt. 0) .and. (nb .le. n_occ) ) then
!!$               allocate(orb(nx,ny,nz))
!!$               call read_orbcube_1(nx, ny, nz, n_at, orb, nb)        
!!$            else
!!$               allocate(virt(nx,ny,nz))
!!$               call read_virtcube_1(nx, ny, nz, n_at, virt, nb, virt_list, n_virt, n_occ)
!!$            end if
!!$            do np=1, n_proj
!!$               amnk(nb,np)=0.d0
!!$               r0x=ctr_proj(np,1)*b1
!!$               r0y=ctr_proj(np,2)*b2
!!$               r0z=ctr_proj(np,3)*b3
!!$               do k=1,nz
!!$                  if ( nb == 1 ) then
!!$                     zz=k*bz(3)-r0z
!!$                  end if
!!$                  do j=1,ny
!!$                     if ( nb == 1 ) then
!!$                        yy=j*by(2)-r0y
!!$                     end if
!!$                     do i=1,nx
!!$                        if ( nb == 1 .and. (pre_check .eqv. .false.) ) then
!!$                           if ( (np==1) .and. (i==1) .and. (j==1) .and. (k==1) ) then
!!$                              allocate(ylm(nx,ny,nz))
!!$                              allocate(func_r(nx,ny,nz))
!!$                              allocate(sph_har(nx,ny,nz,n_proj))
!!$                           end if
!!$                           xx=i*bx(1)-r0x
!!$                           call angularpart(l, mr, np, nx, ny, nz, i, j, k, &
!!$                              &   xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
!!$                           call radialpart(rvalue, zona, np, nx, ny, nz, i, j, k, &
!!$                              &   xx, yy, zz, n_proj, func_r)
!!$                           sph_har(i,j,k,np)=func_r(i,j,k)*ylm(i,j,k)
!!$                           if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then
!!$                              call write_cube(w_sph, w_ang, w_rad, 'sph_har', 'func_r', 'ylm', np, n_proj, &
!!$                                 &   nx, ny, nz, n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
!!$                              if ( np==n_proj) then
!!$                                 deallocate(func_r)
!!$                                 deallocate(ylm)
!!$                              end if
!!$                           end if
!!$                        end if
!!$                        if ( (nb .gt. 0) .and. (nb .le. n_occ) ) then
!!$                           amnk(nb,np)=amnk(nb,np)+orb(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
!!$                        else
!!$                           amnk(nb,np)=amnk(nb,np)+virt(i,j,k)*sph_har(i,j,k,np)*sqrt(bx(1)*by(2)*bz(3))
!!$                        end if
!!$                     end do
!!$                  end do
!!$               end do
!!$               ! sqrt(amnk_tot) verifies the normalization in each band (it must tend to 1)
!!$               amnk_tot(nb)=amnk_tot(nb)+(amnk(nb,np))**2
!!$            end do
!!$            if (allocated(orb)) deallocate(orb)
!!$            if (allocated(virt)) deallocate(virt)
!!$            write(*,'(I4,3x,F12.6)') nb, sqrt(amnk_tot(nb))
!!$         end do
!!$         deallocate(amnk_tot)
!!$         deallocate(sph_har)
!!$      end if
!!$      write(*,*) '!==================================!'
!!$      write(*,*) '!  Calculating amnk=<psi|sph_har>  !'
!!$      write(*,*) '!               done               !'
!!$      write(*,*) '!==================================!'
!!$      write(*,*)
!!$      write(*,*)
!!$
!!$
!!$      ! Write the .amn file
!!$      call write_amn(seedname, n_bands, n_kpts, n_proj, amnk)
!!$      deallocate(amnk)
!!$
!!$      call timing(iproc,'CrtProjectors ','OF')
!!$      call timing(iproc,'CrtDescriptors','ON')
!!$
!!$      ! Calculate the scalar product of wavefunctions psi calculated by BigDFT and themselves.
!!$      ! This values are written in a .mmn file that will be then used by Wannier90.
!!$      ! Quick algorithm to calculate mmnk matrix
!!$      write(*,*) '!==================================!'
!!$      write(*,*) '!   Calculating mmnk=<psi|psi> :   !'
!!$      write(*,*) '!==================================!'
!!$      write(*,*) 'The values of sqrt(mmnk_tot) check the normalization in each band.'
!!$      write(*,*) 'They must tend to their lower limit value, which is equal to 1 :'
!!$      allocate(mmnk_re(n_bands,n_bands,n_kpts*n_nnkpts))
!!$      allocate(mmnk_im(n_bands,n_bands,n_kpts*n_nnkpts))
!!$      allocate(mmnk_tot(n_bands,n_kpts*n_nnkpts))
!!$      ! Algorithm to compute the scalar product :
!!$      do inn=1,n_kpts*n_nnkpts
!!$         write(*,*)
!!$         write(*,'(A21,3(I4,1x))') 'k-point coordinates :', (G_vec(inn,np), np=1,3)!G_vec(inn,1), G_vec(inn,2), G_vec(inn,3)
!!$         write(*,'(A4,4x,A15)') 'Band', 'sqrt(mmnk_tot)='
!!$         if (n_occ == n_virt) then 
!!$            do nb1=1,n_occ
!!$               mmnk_tot(nb1,inn)=0.d0
!!$               mmnk_tot(nb1+n_occ,inn)=0.d0
!!$               if ( (nb1 == 1) .and. (inn == 1) ) allocate(psi(nx,ny,nz,n_bands))
!!$               do nb2=1,n_occ
!!$                  mmnk_re(nb1,nb2,inn)=0.d0
!!$                  mmnk_im(nb1,nb2,inn)=0.d0
!!$                  mmnk_re(nb1+n_occ,nb2,inn)=0.d0
!!$                  mmnk_im(nb1+n_occ,nb2,inn)=0.d0
!!$                  mmnk_re(nb1,nb2+n_occ,inn)=0.d0
!!$                  mmnk_im(nb1,nb2+n_occ,inn)=0.d0
!!$                  mmnk_re(nb1+n_occ,nb2+n_occ,inn)=0.d0
!!$                  mmnk_im(nb1+n_occ,nb2+n_occ,inn)=0.d0
!!$                  do k=1,nz
!!$                     zz=k*bz(3)
!!$                     do j=1,ny
!!$                        yy=j*by(2)
!!$                        do i=1,nx
!!$                           xx=i*bx(1)
!!$                           if ( (inn==1) .and. (nb1 == 1) &
!!$                              &   .and. (i==1) .and. (j==1) .and. (k==1)) then
!!$                           call read_orbcube(nx, ny, nz, n_at, psi, n_bands, nb2)
!!$                           call read_virtcube(nx, ny, nz, n_at, psi, nb2+n_occ, virt_list, n_virt, n_bands, n_occ)
!!$                        end if
!!$                        mmnk_re(nb1,nb2,inn)=mmnk_re(nb1,nb2,inn)+psi(i,j,k,nb1)*psi(i,j,k,nb2)*cos( 2*pi* &
!!$                           &   (xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                        mmnk_im(nb1,nb2,inn)=mmnk_im(nb1,nb2,inn)-psi(i,j,k,nb1)*psi(i,j,k,nb2)*sin( 2*pi* &
!!$                           &   (xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                        mmnk_re(nb1+n_occ,nb2,inn)=mmnk_re(nb1+n_occ,nb2,inn)+psi(i,j,k,nb1+n_occ)* &
!!$                           &   psi(i,j,k,nb2)*cos( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+ &
!!$                           &   zz*G_vec(inn,3)/b3) )
!!$                        mmnk_im(nb1+n_occ,nb2,inn)=mmnk_im(nb1+n_occ,nb2,inn)-psi(i,j,k,nb1+n_occ)* &
!!$                           &   psi(i,j,k,nb2)*sin( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+ &
!!$                           &   zz*G_vec(inn,3)/b3) )
!!$                        mmnk_re(nb1,nb2+n_occ,inn)=mmnk_re(nb1,nb2+n_occ,inn)+psi(i,j,k,nb1)* &
!!$                           &   psi(i,j,k,nb2+n_occ)*cos( 2*pi*(xx*G_vec(inn,1)/b1+ &
!!$                           &   yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                        mmnk_im(nb1,nb2+n_occ,inn)=mmnk_im(nb1,nb2+n_occ,inn)-psi(i,j,k,nb1)* &
!!$                           &   psi(i,j,k,nb2+n_occ)*sin( 2*pi*(xx*G_vec(inn,1)/b1+ &
!!$                           &   yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                        mmnk_re(nb1+n_occ,nb2+n_occ,inn)=mmnk_re(nb1+n_occ,nb2+n_occ,inn)+psi(i,j,k,nb1+n_occ)* &
!!$                           &   psi(i,j,k,nb2+n_occ)*cos( 2*pi*(xx*G_vec(inn,1)/b1+ &
!!$                           &   yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                        mmnk_im(nb1+n_occ,nb2+n_occ,inn)=mmnk_im(nb1+n_occ,nb2+n_occ,inn)-psi(i,j,k,nb1+n_occ)* &
!!$                           &   psi(i,j,k,nb2+n_occ)* sin( 2*pi*(xx*G_vec(inn,1)/b1+ &
!!$                           &   yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                     end do
!!$                  end do
!!$               end do
!!$               ! sqrt(mmnk_tot) verifies the normalization in each band (it must tend to 1)
!!$               mmnk_tot(nb1,inn)=mmnk_tot(nb1,inn)+mmnk_re(nb1,nb2,inn)**2+mmnk_im(nb1,nb2,inn)**2
!!$               mmnk_tot(nb1+n_occ,inn)=mmnk_tot(nb1+n_occ,inn)+mmnk_re(nb1+n_occ,nb2,inn)**2+ &
!!$                  &   mmnk_im(nb1+n_occ,nb2,inn)**2
!!$               mmnk_tot(nb1,inn)=mmnk_tot(nb1,inn)+mmnk_re(nb1,nb2+n_occ,inn)**2+ &
!!$                  &   mmnk_im(nb1,nb2+n_occ,inn)**2
!!$               mmnk_tot(nb1+n_occ,inn)=mmnk_tot(nb1+n_occ,inn)+mmnk_re(nb1+n_occ,nb2+n_occ,inn)**2 &
!!$                  &   +mmnk_im(nb1+n_occ,nb2+n_occ,inn)**2
!!$            end do
!!$         end do
!!$         do nb1=1,n_bands
!!$            write(*,'(I4,3x,F12.6)') nb1, sqrt(mmnk_tot(nb1,inn))       
!!$         end do
!!$      else
!!$         do nb1=1,n_bands
!!$            mmnk_tot(nb1,inn)=0.d0
!!$            if ( (nb1 == 1) .and. (inn == 1) ) then
!!$               allocate(psi(nx,ny,nz,n_bands))
!!$               call read_orbcube(nx, ny, nz, n_at, psi, n_bands, nb1)
!!$            end if
!!$            do nb2=1,n_bands
!!$               mmnk_re(nb1,nb2,inn)=0.d0
!!$               mmnk_im(nb1,nb2,inn)=0.d0
!!$               do k=1,nz
!!$                  zz=k*bz(3)
!!$                  do j=1,ny
!!$                     yy=j*by(2)
!!$                     do i=1,nx
!!$                        xx=i*bx(1)
!!$                        if ( (inn==1) .and. (nb1 == 1) .and. (nb2 .ne. 1) .and. (nb2 .le. n_occ) &
!!$                           &   .and. (i==1) .and. (j==1) .and. (k==1)) then
!!$                           call read_orbcube(nx, ny, nz, n_at, psi, n_bands, nb2)
!!$                        end if
!!$                        if ( (inn==1) .and. (nb1 == 1) .and. (nb2 .gt. n_occ) &
!!$                        &   .and. (i==1) .and. (j==1) .and. (k==1)) then
!!$                           call read_virtcube(nx, ny, nz, n_at, psi, nb2, virt_list, n_virt, n_bands, n_occ)
!!$                        end if
!!$                        mmnk_re(nb1,nb2,inn)=mmnk_re(nb1,nb2,inn)+psi(i,j,k,nb1)*psi(i,j,k,nb2)* &
!!$                        &   cos( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                        mmnk_im(nb1,nb2,inn)=mmnk_im(nb1,nb2,inn)-psi(i,j,k,nb1)*psi(i,j,k,nb2)* &
!!$                        &   sin( 2*pi*(xx*G_vec(inn,1)/b1+yy*G_vec(inn,2)/b2+zz*G_vec(inn,3)/b3) )
!!$                     end do
!!$                  end do
!!$               end do
!!$               ! sqrt(mmnk_tot) verifies the normalization in each band (it must tend to 1)
!!$               mmnk_tot(nb1,inn)=mmnk_tot(nb1,inn)+mmnk_re(nb1,nb2,inn)**2+mmnk_im(nb1,nb2,inn)**2
!!$            end do
!!$            write(*,'(I4,3x,F12.6)') nb1, sqrt(mmnk_tot(nb1,inn))
!!$         end do
!!$      end if
!!$   end do
!!$   deallocate(mmnk_tot)
!!$   write(*,*) '!==================================!'
!!$   write(*,*) '! Calculating mmnk=<psi|psi> done  !'
!!$   write(*,*) '!==================================!'
!!$   write(*,*)
!!$   write(*,*)
!!$
!!$
!!$   ! Write the .mmn file
!!$   call write_mmn(seedname, n_bands, n_kpts, n_nnkpts, k_plus_b, G_vec, mmnk_re, mmnk_im)
!!$   deallocate(mmnk_re)
!!$   deallocate(mmnk_im)
!!$
!!$   call timing(iproc,'CrtDescriptors','OF')
!!$
!!$   ! Write UNKnk.s files containing all the Bloch states calculated by BigDFT. 
!!$   ! These files are used for plotting.
!!$   ! Beware : it is not implemented in the case where there are more than 999 k-points
!!$   do nk=1, n_kpts
!!$      ! s is the spin, set by default to 1
!!$      s=1
!!$      if (w_unk .eqv. .true. ) then 
!!$         call write_unk(n_bands, nx, ny, nz, nk, s, psi)
!!$      end if
!!$   end do
!!$
!!$
!!$   ! Deallocations
!!$   if(allocated(psi))                   deallocate(psi)
!!$   if(allocated(Z))                     deallocate(Z)
!!$   if(allocated(at_pos))                deallocate(at_pos)
!!$   if(allocated(kpts))                  deallocate(kpts)
!!$   if(allocated(ctr_proj))              deallocate(ctr_proj)
!!$   if(allocated(x_proj))                deallocate(x_proj)
!!$   if(allocated(y_proj))                deallocate(y_proj)
!!$   if(allocated(z_proj))                deallocate(z_proj)
!!$   if(allocated(l))                     deallocate(l)
!!$   if(allocated(mr))                    deallocate(mr)
!!$   if(allocated(rvalue))                deallocate(rvalue)
!!$   if(allocated(zona))                  deallocate(zona)
!!$   if(allocated(k_plus_b))              deallocate(k_plus_b)
!!$   if(allocated(G_vec))                 deallocate(G_vec)
!!$   if(allocated(excb))                  deallocate(excb)
!!$   if(allocated(amnk_bands_sorted))     deallocate(amnk_bands_sorted)
!!$
!!$else
!!$   if (iproc==0) write(*,*) 'Cubic code not parallelized'
!!$end if


!subroutine angularpart(l, mr, np, nx, ny, nz, ix, iy, iz, &
!                    xx, yy, zz, n_proj, ylm)
!
!   ! This routine returns the angular part of the spherical harmonic identified by indices (l,mr)
!   ! Calcutations are made in Cartesian coordinates
!
!   implicit none
!
!   ! I/O variables
!   integer, intent(in) :: l(n_proj), mr(n_proj)
!   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
!   real(kind=8), intent(in) :: xx, yy, zz
!   real(kind=8), dimension(nx,ny,nz), intent(out) :: ylm
!
!   ! local variables
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8), parameter :: eps8  = 1.0e-8
!   real(kind=8), external :: s, pz, px, py, dz2, dxz, dyz, dx2my2, dxy, &
!                          fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
!   real(kind=8) :: rr
!   real(kind=8) :: bs2, bs3, bs6, bs12
!
!   bs2 = 1.d0/sqrt(2.d0)
!   bs3 = 1.d0/sqrt(3.d0)
!   bs6 = 1.d0/sqrt(6.d0)
!   bs12 = 1.d0/sqrt(12.d0)
!
!
!   if (l(np) > 3 .OR. l(np) < -5 ) then 
!      write(*,*) 'error, l out of range '
!   else
!      if (l(np)>=0) then
!         if (mr(np) < 1 .OR. mr(np) > 2*l(np)+1) then
!	    write(*,*) 'error, mr out of range'
!	 end if
!      else
!         if (mr(np) < 1 .OR. mr(np) > abs(l(np))+1 ) then 
!	    write(*,*) 'error, mr out of range'
!         end if
!      end if
!   end if
!
!   rr = sqrt( xx*xx + yy*yy + zz*zz )
!    
!   if (l(np)==0) then   ! s orbital
!      ylm(ix,iy,iz) = s(xx,yy,zz,rr)  
!   end if
!
!   if (l(np)==1) then   ! p orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = pz(xx,yy,zz,rr) 
!      if (mr(np)==2) ylm(ix,iy,iz) = px(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = py(xx,yy,zz,rr)
!   end if
!
!   if (l(np)==2) then   ! d orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = dz2(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = dxz(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = dyz(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = dx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = dxy(xx,yy,zz,rr)
!   endif
!
!   if (l(np)==3) then   ! f orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = fz3(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = fxz2(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = fyz2(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = fzx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = fxyz(xx,yy,zz,rr)
!      if (mr(np)==6) ylm(ix,iy,iz) = fxx2m3y2(xx,yy,zz,rr)
!      if (mr(np)==7) ylm(ix,iy,iz) = fy3x2my2(xx,yy,zz,rr)
!   endif
!
!   if (l(np)==-1) then  !  sp hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs2 * ( s(xx,yy,zz,rr) + px(xx,yy,zz,rr) ) 
!      if (mr(np)==2) ylm(ix,iy,iz) = bs2 * ( s(xx,yy,zz,rr) - px(xx,yy,zz,rr) ) 
!   end if
!
!   if (l(np)==-2) then  !  sp2 hybrids 
!      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr) +2.d0*bs6*px(xx,yy,zz,rr) 
!   end if
!
!   if (l(np)==-3) then  !  sp3 hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)+px(xx,yy,zz,rr)+py(xx,yy,zz,rr)+pz(xx,yy,zz,rr))
!      if (mr(np)==2) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)+px(xx,yy,zz,rr)-py(xx,yy,zz,rr)-pz(xx,yy,zz,rr))
!      if (mr(np)==3) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)-px(xx,yy,zz,rr)+py(xx,yy,zz,rr)-pz(xx,yy,zz,rr))
!      if (mr(np)==4) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)-px(xx,yy,zz,rr)-py(xx,yy,zz,rr)+pz(xx,yy,zz,rr))
!   end if
!
!   if (l(np)==-4) then  !  sp3d hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr) +2.d0*bs6*px(xx,yy,zz,rr) 
!      if (mr(np)==4) ylm(ix,iy,iz) = bs2*pz(xx,yy,zz,rr)+bs2*dz2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) =-bs2*pz(xx,yy,zz,rr)+bs2*dz2(xx,yy,zz,rr)
!   end if
!
!   if (l(np)==-5) then  ! sp3d2 hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*px(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)+.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*px(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)+.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)-.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)-.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*pz(xx,yy,zz,rr)+bs3*dz2(xx,yy,zz,rr)
!      if (mr(np)==6) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*pz(xx,yy,zz,rr)+bs3*dz2(xx,yy,zz,rr)
!   end if
!
!END SUBROUTINE angularpart
!
!
!! The following functions are used to calculate angular parts of the spherical harmonics
!!======== l = 0 =====================================================================
!function s(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) :: s, xx, yy, zz, rr
!   s = 1.d0/ sqrt(4*pi)
!END FUNCTION s
!
!
!!======== l = 1 =====================================================================
!function pz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::pz, xx, yy, zz, rr
!   pz =  sqrt(3.d0/(4*pi)) * (zz/rr)
!END FUNCTION pz
!
!function px(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::px, xx, yy, zz, rr
!   px =  sqrt(3.d0/(4*pi)) * (xx/rr)
!END FUNCTION px
!
!function py(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::py, xx, yy, zz, rr
!   py =  sqrt(3.d0/(4*pi)) * (yy/rr)
!END FUNCTION py
!
!
!!======== l = 2 =====================================================================
!function dz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dz2, xx, yy, zz, rr
!   dz2 =  sqrt(1.25d0/(4*pi)) * (-xx*xx-yy*yy+2.d0*zz*zz)/(rr*rr)
!END FUNCTION dz2
!
!function dxz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dxz, xx, yy, zz, rr
!   dxz =  sqrt(15.d0/(4*pi)) * (xx*zz)/(rr*rr)
!END FUNCTION dxz
!
!function dyz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dyz, xx, yy, zz, rr
!   dyz =  sqrt(15.d0/(4*pi)) * (yy*zz)/(rr*rr)
!END FUNCTION dyz
!
!function dx2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dx2my2, xx, yy, zz, rr
!   dx2my2 =  sqrt(3.75d0/(4*pi)) * (xx*xx-yy*yy)/(rr*rr)
!END FUNCTION dx2my2
!
!function dxy(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dxy, xx, yy, zz, rr
!   dxy =  sqrt(3.75d0/(4*pi)) * (xx*yy)/(rr*rr)
!END FUNCTION dxy
!
!
!!======== l = 3 =====================================================================
!function fz3(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fz3, xx, yy, zz, rr
!   fz3 =  0.25d0*sqrt(7.d0/pi) * (zz*(2.d0*zz*zz-3.d0*xx*xx-3.d0*yy*yy)) / (rr*rr*rr)
!END FUNCTION fz3
!
!function fxz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxz2, xx, yy, zz, rr
!   fxz2 =  0.25d0*sqrt(10.5d0/pi) * (xx*(4.d0*zz*zz-xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fxz2
!
!function fyz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fyz2, xx, yy, zz, rr
!   fyz2 =  0.25d0*sqrt(10.5d0/pi) * (yy*(4.d0*zz*zz-xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fyz2
!
!function fzx2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fzx2my2, xx, yy, zz, rr
!   fzx2my2 =  0.25d0*sqrt(105d0/pi) * (zz*(xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fzx2my2
!
!function fxyz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxyz, xx, yy, zz, rr
!   fxyz =  0.25d0*sqrt(105d0/pi) * (xx*yy*zz) / (rr*rr*rr)
!END FUNCTION fxyz
!
!function fxx2m3y2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxx2m3y2, xx, yy, zz, rr
!   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * (xx*(xx*xx-3.d0*yy*yy)) / (rr*rr*rr)
!END FUNCTION fxx2m3y2
!
!function fy3x2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fy3x2my2, xx, yy, zz, rr
!   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * (yy*(3.d0*xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fy3x2my2


!!subroutine make_precheck(iproc,nproc,input,Glr,orbsv,commsv,orbsp,commsp,atoms,w,rxyz,n_proj,ctr_proj,&
!!           x_proj,y_proj,z_proj,l,mr,rvalue,zona,amnk_bands_sorted,sph_daub)
!!   use BigDFT_API
!!   use Poisson_Solver
!!   implicit none
!!   integer, intent(in) :: iproc, nproc, n_proj
!!   type(input_variables),intent(in) :: input
!!   type(locreg_descriptors), intent(in) :: Glr
!!   type(orbitals_data), intent(inout) :: orbsv,orbsp
!!   type(communications_arrays), target :: commsv,commsp
!!   type(atoms_data), intent(in) :: atoms
!!   type(workarr_sumrho), intent(in) :: w
!!   real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!   real(kind=8), dimension (n_proj,3), intent(in) :: ctr_proj, x_proj, y_proj, z_proj
!!   integer, dimension (n_proj),intent(in) :: l, mr, rvalue
!!   real, dimension (n_proj), intent(in) :: zona
!!   integer, dimension (:), pointer :: amnk_bands_sorted
!!   real(wp),dimension(:),pointer :: sph_daub
!!   !local variables
!!   character(len=*), parameter :: subname='make_precheck'
!!   integer :: i_stat, i_all, npsidim, i, j, k, np, npp, pshft
!!   integer :: ind, nb, ierr, npsidim2,nvctrp
!!   real(kind=8) :: b1, b2, b3, r0x, r0y, r0z, zz, yy, xx
!!   real(wp), allocatable :: psi_etsfv(:,:),sph_har_etsf(:)!,sph_daub(:)
!!   real(wp), allocatable :: psi_etsf2(:)
!!   real(wp), pointer :: pwork(:)
!!   character(len=60) :: filename
!!   real(kind=8), allocatable :: ylm(:,:,:), func_r(:,:,:)
!!   real(kind=8), allocatable :: amnk(:,:), amnk_guess(:)
!!   integer :: n_virt, n_virt_tot
!!   real(kind=8), allocatable, dimension(:) :: amnk_guess_sorted
!!   real(gp), dimension(3,atoms%nat) :: rxyz_old
!!
!!   call timing(iproc,'CrtProjectors ','ON')
!!
!!   ! Read wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
!!   allocate(psi_etsfv(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,max(orbsv%norbp*orbsv%nspinor,1)),stat=i_stat)
!!   call memocc(i_stat,psi_etsfv,'psi_etsfv',subname)
!!   if(associated(orbsv%eval)) nullify(orbsv%eval)
!!   allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
!!   call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)
!!
!!   filename= trim(input%dir_output) // 'virtuals'// trim(wfformat_read)
!!   call readmywaves(iproc,filename,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
!!      Glr%wfd,psi_etsfv)
!!   i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
!!   deallocate(orbsv%eval,stat=i_stat)
!!   nullify(orbsv%eval)
!!   call memocc(i_stat,i_all,'orbsv%eval',subname)
!!
!!   ! Tranposition of the distribution of the BigDFT wavefunctions : orbitals -> components.
!!   npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsv%norbp*orbsv%nspinor,sum(commsv%ncntt(0:nproc-1)))
!!   npsidim2=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsp%norbp,sum(commsp%ncntt(0:nproc-1)))
!!   allocate(psi_etsf2(npsidim),stat=i_stat) !!doing this because psi_etsfv does not incorporate enough space for transpose
!!   call memocc(i_stat,psi_etsf2,'psi_etsf2',subname)
!!
!!   call to_zero(npsidim,psi_etsf2)
!!   if(nproc > 1) then
!!     allocate(pwork(npsidim),stat=i_stat)
!!     call memocc(i_stat,pwork,'pwork',subname)
!!     call transpose_v(iproc,nproc,orbsv,Glr%wfd,commsv,psi_etsfv(1,1),work=pwork,outadd=psi_etsf2(1))
!!     i_all = -product(shape(pwork))*kind(pwork)
!!     deallocate(pwork,stat=i_stat)
!!     call memocc(i_stat,i_all,'pwork',subname)
!!   else
!!      ! just copy the wavefunctions 
!!      k=0
!!      do j=1,orbsv%norbp*orbsv%nspinor
!!      do i=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!         k=k+1
!!         psi_etsf2(k) = psi_etsfv(i,j)
!!      end do
!!      end do
!!   end if
!!   i_all=-product(shape(psi_etsfv))*kind(psi_etsfv)
!!   deallocate(psi_etsfv,stat=i_stat)
!!   call memocc(i_stat,i_all,'psi_etsfv',subname)
!!
!!   ! - b1, b2 and b3 are the norm of the lattice parameters.
!!   b1=atoms%alat1
!!   b2=atoms%alat2
!!   b3=atoms%alat3
!!   ! - Allocations
!!   allocate(amnk(orbsv%norb,orbsp%norb),stat=i_stat)
!!   call memocc(i_stat,amnk,'amnk',subname)
!!   allocate(amnk_guess(orbsv%norb),stat=i_stat)
!!   call memocc(i_stat,amnk_guess,'amnk_guess',subname)
!!   allocate(sph_daub(npsidim2), stat=i_stat)
!!   call memocc(i_stat,sph_daub,'sph_daub',subname)
!!
!!   ! Begining of the algorithm to compute the scalar product in order to find the best unoccupied orbitals to use to compute the actual Amnk matrix :
!!   if (iproc==0) then
!!      write(*,*) '!==================================!'
!!      write(*,*) '! Calculating amnk=<virt|sph_har>  !'
!!      write(*,*) '!       in pre-check mode :        !'
!!      write(*,*) '!==================================!'
!!      write(*,'(A12,4x,A15)') 'Virtual band', 'amnk_guess(nb)='
!!   end if
!!
!!   ! Calculation of the spherical harmonics in parallel.
!!   ! It is done in the real space and then converted in the Daubechies representation.
!!   pshft = 0
!!   do npp=1, orbsp%norbp
!!      np = npp + orbsp%isorb
!!      ! Convolution buffer : n1i=2*n1+31 -> explains the '13*input%hx*0.5' term
!!      r0x=ctr_proj(np,1)*b1+13*input%hx*0.5
!!      r0y=ctr_proj(np,2)*b2+13*input%hy*0.5
!!      r0z=ctr_proj(np,3)*b3+13*input%hz*0.5
!!      do k=1,Glr%d%n3i
!!         zz=(k-1)*input%hz*0.5-r0z
!!         do j=1,Glr%d%n2i
!!            yy=(j-1)*input%hy*0.5-r0y
!!            do i=1,Glr%d%n1i
!!               ind=(k-1)*Glr%d%n2i*Glr%d%n1i+(j-1)*Glr%d%n1i+i
!!               xx=(i-1)*input%hx*0.5-r0x
!!               call angularpart(l, mr, np, Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, i, j, k, &
!!                     xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
!!               call radialpart(rvalue, zona, np, Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, i, j, k, &
!!                     xx, yy, zz, n_proj, func_r)
!!               ! The 'sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)' term is here to normalize spherical harmonics
!!               sph_har_etsf(ind)=func_r(i,j,k)*ylm(i,j,k)*sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)
!!            end do
!!         end do
!!      end do
!!      call isf_to_daub(Glr,w,sph_har_etsf(1),sph_daub(1+pshft))
!!      pshft=pshft + max(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,commsp%ncntt(iproc)/orbsp%norbp)
!!   end do
!!
!!   call timing(iproc,'CrtProjectors ','OF')
!!
!!   ! Tranposition of the distribution of the spherical harmonics: orbitals -> components.
!!   allocate(pwork(npsidim2),stat=i_stat)
!!   call memocc(i_stat,pwork,'pwork',subname)
!!   call transpose_v(iproc,nproc,orbsp,Glr%wfd,commsp,sph_daub,work=pwork)
!!   i_all = -product(shape(pwork))*kind(pwork)
!!   deallocate(pwork,stat=i_stat)
!!   call memocc(i_stat,i_all,'pwork',subname)
!!   call timing(iproc,'ApplyProj     ','ON')
!!
!!   ! Scalar product of amnk=<sph_daub|psi> in parallel.
!!   call to_zero(orbsp%norb*orbsv%norb,amnk)
!!   nvctrp=commsv%nvctr_par(iproc,1)
!!   call gemm('T','N',orbsv%norb,orbsp%norb,nvctrp,1.0_wp,psi_etsf2(1),max(1,nvctrp),&
!!        sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),orbsv%norb)
!!      
!!   ! Construction of the whole Amnk_guess matrix.
!!   call mpiallred(amnk(1,1),orbsv%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
!!
!!   ! For each unoccupied orbitals, check how they project on spherical harmonics.
!!   ! The greater amnk_guess(nb) is, the more they project on spherical harmonics.
!!   do nb=1,orbsv%norb
!!      amnk_guess(nb)=0.0
!!      do np=1,orbsp%norb
!!         amnk_guess(nb)=amnk_guess(nb)+(amnk(nb,np))**2
!!      end do
!!      if (iproc==0) write(*,'(I4,11x,F12.6)') nb, sqrt(amnk_guess(nb))
!!   end do
!!
!!   ! Choice of the unoccupied orbitals to calculate the Amnk matrix
!!   if (iproc==0) then
!!      write(*,*) 
!!      write(*,'(1a)') 'These are the virtual bands to use to construct the actual Amn and Mmn matrices :'
!!      write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_guess)='
!!   end if
!!   allocate(amnk_guess_sorted(n_virt),stat=i_stat)
!!   call memocc(i_stat,amnk_guess_sorted,'amnk_guess_sorted',subname)
!!   do nb=1,n_virt
!!      amnk_guess_sorted(nb)=maxval(amnk_guess,1)
!!      amnk_bands_sorted(nb)=maxloc(amnk_guess,1)
!!      amnk_guess(amnk_bands_sorted(nb))=0.d0
!!   if (iproc==0) write(*,'(I4,3x,F12.6)') amnk_bands_sorted(nb), sqrt(amnk_guess_sorted(nb))
!!   end do
!!
!!   ! End of the pre-check mode
!!   i_all = -product(shape(psi_etsf2))*kind(psi_etsf2)
!!   deallocate(psi_etsf2,stat=i_stat)
!!   call memocc(i_stat,i_all,'psi_etsf2',subname)
!!   i_all = -product(shape(amnk))*kind(amnk)
!!   deallocate(amnk,stat=i_stat)
!!   call memocc(i_stat,i_all,'amnk',subname)
!!   i_all = -product(shape(amnk_guess_sorted))*kind(amnk_guess_sorted)
!!   deallocate(amnk_guess_sorted,stat=i_stat)
!!   call memocc(i_stat,i_all,'amnk_guess_sorted',subname)
!!   i_all = -product(shape(amnk_guess))*kind(amnk_guess)
!!   deallocate(amnk_guess,stat=i_stat)
!!   call memocc(i_stat,i_all,'amnk_guess',subname)
!!
!!   if (iproc==0) then
!!      write(*,*) '!==================================!'
!!      write(*,*) '! Calculating amnk=<virt|sph_har>  !'
!!      write(*,*) '!     in pre-check mode done       !'
!!      write(*,*) '!==================================!'
!!      write(*,*)
!!      write(*,*)
!!   end if
!!
!!   ! Rewrite the input.inter file to add the chosen unoccupied states.
!!   if (iproc==0) call write_inter(n_virt, amnk_bands_sorted)
!!
!!   call timing(iproc,'ApplyProj     ','OF')
!!
!!END SUBROUTINE make_precheck

!>
!!subroutine limitations(seedname, n_proj, n_occ, n_virt, n_bands, n_kpts)
!!
!!   ! This routine sends a message in case there are problems with some values
!!
!!   implicit none
!!
!!   ! I/O variables
!!   character, intent(in) :: seedname*16
!!   integer, intent(in) :: n_proj, n_occ, n_virt, n_bands, n_kpts
!!
!!
!!   ! Beware : it is not implemented in the case where there are more than 9999 occupied bands 
!!   if (n_occ>9999) then
!!      write(*,*) 'There are too many occupied bands'
!!      STOP
!!   else
!!      if (n_occ<=0) then
!!         write(*,*) 'Wrong number of occupied bands (lower than 0)'
!!         STOP
!!      end if
!!   end if
!!
!!   ! Beware : it is not implemented in the case where there are more than 9999 unoccupied bands 
!!   if (n_virt>9999) then
!!      write(*,*) 'There are too many virtual bands'
!!      STOP
!!   else
!!      if (n_virt<0) then
!!         write(*,*) 'Wrong number of virtual bands (lower than 0)'
!!         STOP
!!      end if
!!   end if
!!
!!   ! Beware : it is not implemented in the case where there are more than 9999 projections 
!!   if (n_proj>9999) then
!!      write(*,*) 'There are too many projections'
!!      STOP
!!   else
!!      if (n_proj<0) then
!!         write(*,*) 'Wrong number of projections (lower than 0)'
!!         STOP
!!      end if
!!   end if
!!
!!   ! Beware : it is not implemented in the case where there are more than 9999 k-points 
!!   if (n_kpts>9999) then
!!      write(*,*) 'There are too many k-points '
!!      STOP
!!   else
!!      if (n_kpts<0) then
!!         write(*,*) 'Wrong number of k-points (lower than 0)'
!!         STOP
!!      end if
!!   end if
!!
!!   ! Beware : there cannot be more projections than bands
!!   if (n_proj > n_bands) then
!!      write(*,*) 'There cannot be more projections than bands'
!!      write(*,*) 'Possible reasons for problems :'
!!      write(*,*) '- Wrong specifications for the number of occupied and unoccupied orbitals in input.inter'
!!      write(*,*) '- Wrong specifications for the projections in ',  trim(seedname)//'.win'
!!      write(*,*) 'In this second case, do not forget to restart Wannier90 in the Post-processing mode before restarting&
!!         &   this interface program'
!!      STOP
!!   end if
!!
!!END SUBROUTINE limitations

!>
!!subroutine read_cube_header_1(nx, ny, nz, n_at, bx, by, bz)
!!
!!   ! This routine reads the first lines of a .cube file
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(out) :: nx, ny, nz, n_at
!!   real(kind=8), intent(out) :: bx(3), by(3), bz(3)
!!
!!   ! Local variables
!!   integer :: i
!!
!!
!!   OPEN(11, FILE='orbital0001.cube', STATUS='OLD')
!!
!!   write(*,*) '!==================================!'
!!   write(*,*) '!       Reading .cube header :     !'
!!   write(*,*) '!==================================!'
!!
!!   read(11,*) ! skip first line
!!   read(11,*) ! skip second line
!!   read(11,*) n_at
!!   read(11,*) nx, (bx(i), i=1,3)
!!   read(11,*) ny, (by(i), i=1,3)
!!   read(11,*) nz, (bz(i), i=1,3)
!!   print *, 'Values read :'
!!   write(*,'(A5,4(I4,A5))') 'n_at=',  n_at, '; nx=', nx, '; ny=', ny, '; nz=', nz
!!   CLOSE(11)
!!
!!END SUBROUTINE read_cube_header_1
!!
!!
!!
!!subroutine read_cube_header_2(n_at, Z, at_pos)
!!
!!   ! This routine reads the first lines of a .cube file
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(in) :: n_at
!!   integer, dimension(n_at), intent(out) :: Z
!!   real(kind=8), dimension(n_at,3), intent(out) :: at_pos
!!
!!   ! Local variables
!!   integer :: i,j
!!   real :: zero
!!
!!
!!   OPEN(11, FILE='orbital0001.cube', STATUS='OLD')
!!
!!   do i=1,6
!!      read(11,*) ! skip first lines
!!   end do
!!   read(11,*)   (Z(i), zero, (at_pos(i,j), j=1,3), i=1,n_at)
!!   CLOSE(11)
!!
!!   write(*,*) '!==================================!'
!!   write(*,*) '!    Reading .cube header done     !'
!!   write(*,*) '!==================================!'
!!   print *
!!   print *
!!
!!END SUBROUTINE read_cube_header_2


!>
!!subroutine read_orbcube_1(nx, ny, nz, n_at, psi, io)
!!
!!   ! This routine reads a .cube file
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(in) :: nx, ny, nz, n_at, io
!!   real(kind=8), dimension(nx, ny, nz), intent(out) :: psi
!!
!!   ! Local variables
!!   integer :: i, j, k
!!   character(len=13) :: subname
!!   character(len=3) :: io_c
!!
!!
!!   ! Finds the name of the file to be read
!!   if (io>0 .and. io<10) then 
!!      write(io_c, '(i1)') io
!!      subname='orbital000'//io_c
!!   else  
!!      if (io>9 .and. io<100) then
!!         write(io_c, '(i2)') io
!!         subname='orbital00'//io_c
!!      else
!!         if (io>99 .and. io<1000) then
!!            write(io_c, '(i3)') io
!!            subname='orbital0'//io_c
!!         else
!!            if (io>999 .and. io<10000) then
!!               write(io_c, '(i4)') io
!!               subname='orbital'//io_c
!!            end if
!!         end if
!!      end if
!!   end if
!!
!!   ! Reading of the file
!!   OPEN(11, FILE=subname//'.cube', STATUS='OLD')
!!   do i=1,6+n_at
!!      read(11,*)   ! skip header
!!   end do
!!   read(11,*)   (((psi(i,j,k), k=1,nz), j=1,ny), i=1,nx)
!!   CLOSE(11)
!!
!!END SUBROUTINE read_orbcube_1


!>
!!subroutine read_orbcube(nx, ny, nz, n_at, psi, n_bands, io)
!!
!!   ! This routine reads a .cube file
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(in) :: nx, ny, nz, n_at, io, n_bands
!!   real(kind=8), dimension(nx, ny, nz, n_bands), intent(out) :: psi
!!
!!   ! Local variables
!!   integer :: i, j, k
!!   character(len=13) :: subname
!!   character(len=3) :: io_c
!!
!!
!!   ! Finds the name of the file to be read
!!   if (io>0 .and. io<10) then 
!!      write(io_c, '(i1)') io
!!      subname='orbital000'//io_c
!!   else  
!!      if (io>9 .and. io<100) then
!!         write(io_c, '(i2)') io
!!         subname='orbital00'//io_c
!!      else
!!         if (io>99 .and. io<1000) then
!!            write(io_c, '(i3)') io
!!            subname='orbital0'//io_c
!!         else
!!            if (io>999 .and. io<10000) then
!!               write(io_c, '(i4)') io
!!               subname='orbital'//io_c
!!            end if
!!         end if
!!      end if
!!   end if
!!
!!   ! Reading of the file
!!   OPEN(11, FILE=subname//'.cube', STATUS='OLD')
!!
!!   do i=1,6+n_at
!!      read(11,*)   ! skip header
!!   end do
!!   read(11,*)   (((psi(i,j,k,io), k=1,nz), j=1,ny), i=1,nx)
!!
!!   CLOSE(11)
!!
!!END SUBROUTINE read_orbcube


!>
!!subroutine read_virtcube_1(nx, ny, nz, n_at, psi, nb, virt_list, n_virt_tot, n_occ)
!!
!!   ! This routine reads a .cube file
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(in) :: nx, ny, nz, n_at, nb, n_virt_tot, n_occ
!!   integer, dimension(n_virt_tot), intent (in) :: virt_list
!!   real(kind=8), dimension(nx, ny, nz), intent(out) :: psi
!!
!!   ! Local variables
!!   integer :: i, j, k, iv
!!   character(len=13) :: subname
!!   character(len=3) :: vl_c
!!
!!   iv=nb-n_occ
!!
!!
!!   ! Finds the name of the file to be read
!!   if (virt_list(iv)>0 .and. virt_list(iv)<10) then 
!!      write(vl_c, '(i1)') virt_list(iv)
!!      subname='virtual000'//vl_c
!!   else  
!!      if (virt_list(iv)>9 .and. virt_list(iv)<100) then
!!         write(vl_c, '(i2)') virt_list(iv)
!!         subname='virtual00'//vl_c
!!      else
!!         if (virt_list(iv)>99 .and. virt_list(iv)<1000) then
!!            write(vl_c, '(i3)') virt_list(iv)
!!            subname='virtual0'//vl_c
!!         else
!!            if (virt_list(iv)>999 .and. virt_list(iv)<10000) then
!!               write(vl_c, '(i4)') virt_list(iv)
!!               subname='virtual'//vl_c
!!            end if
!!         end if
!!      end if
!!   end if
!!
!!
!!   ! Reading of the file
!!   OPEN(11, FILE=subname//'.cube', STATUS='OLD')
!!
!!   do i=1,6+n_at
!!      read(11,*)   ! skip header
!!   end do
!!   read(11,*)   (((psi(i,j,k), k=1,nz), j=1,ny), i=1,nx)
!!
!!   CLOSE(11)
!!
!!END SUBROUTINE read_virtcube_1
!!
!!
!!!> This routine reads a .cube file
!!subroutine read_virtcube(nx, ny, nz, n_at, psi, nb, virt_list, n_virt, n_bands, n_occ)
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(in) :: nx, ny, nz, n_at, nb, n_virt, n_bands, n_occ
!!   integer, dimension(n_virt), intent (in) :: virt_list
!!   real(kind=8), dimension(nx, ny, nz, n_bands), intent(out) :: psi
!!
!!   ! Local variables
!!   integer :: i, j, k, iv
!!   character(len=13) :: subname
!!   character(len=3) :: vl_c
!!
!!   iv=nb-n_occ
!!
!!   ! Finds the name of the file to be read
!!   if (virt_list(iv)>0 .and. virt_list(iv)<10) then 
!!      write(vl_c, '(i1)') virt_list(iv)
!!      subname='virtual000'//vl_c
!!   else  
!!      if (virt_list(iv)>9 .and. virt_list(iv)<100) then
!!         write(vl_c, '(i2)') virt_list(iv)
!!         subname='virtual00'//vl_c
!!      else
!!         if (virt_list(iv)>99 .and. virt_list(iv)<1000) then
!!            write(vl_c, '(i3)') virt_list(iv)
!!            subname='virtual0'//vl_c
!!         else
!!            if (virt_list(iv)>999 .and. virt_list(iv)<10000) then
!!               write(vl_c, '(i4)') virt_list(iv)
!!               subname='virtual'//vl_c
!!            end if
!!         end if
!!      end if
!!   end if
!!
!!
!!   ! Reading of the file
!!   OPEN(11, FILE=subname//'.cube', STATUS='OLD')
!!
!!   do i=1,6+n_at
!!      read(11,*)   ! skip header
!!   end do
!!   read(11,*)   (((psi(i,j,k,nb), k=1,nz), j=1,ny), i=1,nx)
!!
!!   CLOSE(11)
!!
!!END SUBROUTINE read_virtcube

!!subroutine write_cube(w_sph, w_ang, w_rad, fn1, fn2, fn3, np, n_proj, nx, ny, nz, &
!!      &   n_at, bx, by, bz, Z, at_pos, sph_har, func_r, ylm)
!!
!!   ! This routine writes .cube files for radial part, angular part and
!!   ! the spherical harmonic given in argument
!!
!!   implicit none
!!
!!   ! I/O variables
!!   logical, intent(in) :: w_sph, w_ang, w_rad
!!   character(len=7), intent(in) :: fn1
!!   character(len=6), intent(in) :: fn2
!!   character(len=3), intent(in) :: fn3
!!   integer, intent(in) :: np, nx, ny, nz, n_at, n_proj
!!   real(kind=8), dimension(3), intent(in) :: bx, by, bz
!!   integer, dimension(n_at), intent(in) :: Z
!!   real(kind=8), dimension(n_at,3), intent(in) :: at_pos
!!   real(kind=8), dimension(nx,ny,nz,n_proj), intent(in) :: sph_har
!!   real(kind=8), dimension(nx,ny,nz), intent(in) :: func_r, ylm
!!
!!   ! Local variables
!!   character(len=13) :: subname1
!!   character(len=12) :: subname2
!!   character(len=9) :: subname3
!!   character(len=3) :: np_c
!!   integer :: i, j, rem, ix, iy, iz
!!
!!
!!   ! Concatenations to write the names of .cube files
!!   if (np>0 .and. np<10) then 
!!      write(np_c, '(i1)') np
!!      subname1=fn1//'000'//np_c
!!      subname2=fn2//'000'//np_c
!!      subname3=fn3//'000'//np_c
!!   else  
!!      if (np>9 .and. np<100) then
!!         write(np_c, '(i2)') np
!!         subname1=fn1//'00'//np_c
!!         subname2=fn2//'00'//np_c
!!         subname3=fn3//'00'//np_c
!!      else
!!         if (np>99 .and. np<1000) then
!!            write(np_c, '(i3)') np
!!            subname1=fn1//'0'//np_c
!!            subname2=fn2//'0'//np_c
!!            subname3=fn3//'0'//np_c
!!         else
!!            if (np>999 .and. np<10000) then
!!               write(np_c, '(i4)') np
!!               subname1=fn1//np_c
!!               subname2=fn2//np_c
!!               subname3=fn3//np_c
!!            end if
!!         end if
!!      end if
!!   end if
!!
!!   rem=nz-floor(nz/6.d0)*6
!!
!!   ! Write the sph_harxxx.cube files
!!   if (w_sph) then
!!      OPEN(12, FILE=subname1//'.cube', STATUS='unknown')
!!      write(12,*) ' CUBE file for ISF field'
!!      write(12,*) ' Case for'
!!      write(12,'(I4,1X,F12.6,2(1X,F12.6))') n_at, real(0.d0), real(0.d0), real(0.d0)
!!      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nx, bx(1), bx(2), bx(3)
!!      write(12,'(I4,1X,F12.6,2(1X,F12.6))') ny, by(1), by(2), by(3)
!!      write(12,'(I4,1X,F12.6,2(1X,F12.6))') nz, bz(1), bz(2), bz(3)
!!      do i=1, n_at
!!         write(12,'(I4,1X,F12.6,3(1X,F12.6))') Z(i), real(0.d0), (real(at_pos(i,j)), j=1,3)
!!      end do
!!      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
!!      do ix=nx,1,-1
!!         do iy=ny,1,-1
!!            do iz=nz,1,-1
!!               write(12,'(E14.6)',advance='no') real(sph_har(ix,iy,iz,np))
!!               if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. nz) ) .or. (iz .eq. 1) ) then
!!                  write(12,'(a)') ''
!!               end if
!!            end do
!!         end do
!!      end do
!!      CLOSE(12)
!!   end if
!!
!!   ! Write the func_rxxx.cube file
!!   if (w_rad) then
!!      OPEN(13, FILE=subname2//'.cube', STATUS='unknown')
!!      write(13,*) ' CUBE file for ISF field'
!!      write(13,*) ' Case for'
!!      write(13,'(I4,1X,F12.6,2(1X,F12.6))') n_at, real(0.d0), real(0.d0), real(0.d0)
!!      write(13,'(I4,1X,F12.6,2(1X,F12.6))') nx, bx(1), bx(2), bx(3)
!!      write(13,'(I4,1X,F12.6,2(1X,F12.6))') ny, by(1), by(2), by(3)
!!      write(13,'(I4,1X,F12.6,2(1X,F12.6))') nz, bz(1), bz(2), bz(3)
!!      do i=1, n_at
!!         write(13,'(I4,1X,F12.6,3(1X,F12.6))') Z(i), real(0.d0), (real(at_pos(i,j)), j=1,3)
!!      end do
!!      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
!!      do ix=nx,1,-1
!!         do iy=ny,1,-1
!!            do iz=nz,1,-1
!!               write(13,'(E14.6)',advance='no') real(func_r(ix,iy,iz))
!!               if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. nz) ) .or. (iz .eq. 1) ) then
!!                  write(13,'(a)') ''
!!               end if
!!            end do
!!         end do
!!      end do
!!      CLOSE(13)
!!   end if
!!
!!   ! Write the ylmxxx.cube file
!!   if (w_ang) then
!!      OPEN(14, FILE=subname3//'.cube', STATUS='unknown')
!!      write(14,*) ' CUBE file for ISF field'
!!      write(14,*) ' Case for'
!!      write(14,'(I4,1X,F12.6,2(1X,F12.6))') n_at, real(0.d0), real(0.d0), real(0.d0)
!!      write(14,'(I4,1X,F12.6,2(1X,F12.6))') nx, bx(1), bx(2), bx(3)
!!      write(14,'(I4,1X,F12.6,2(1X,F12.6))') ny, by(1), by(2), by(3)
!!      write(14,'(I4,1X,F12.6,2(1X,F12.6))') nz, bz(1), bz(2), bz(3)
!!      do i=1, n_at
!!         write(14,'(I4,1X,F12.6,3(1X,F12.6))') Z(i), real(0.d0), (real(at_pos(i,j)), j=1,3)
!!      end do
!!      ! Volumetric data in batches of 6 values per line, 'z'-direction first.
!!      do ix=nx,1,-1
!!         do iy=ny,1,-1
!!            do iz=nz,1,-1
!!               write(14,'(E14.6)',advance='no') real(ylm(ix,iy,iz))
!!               if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. nz) ) .or. (iz .eq. 1) ) then
!!                  write(14,'(a)') ''
!!               end if
!!            end do
!!         end do
!!      end do
!!      CLOSE(14)
!!   end if
!!
!!END SUBROUTINE write_cube

!>
!!subroutine write_unk(n_bands, nx, ny, nz, nk, s, psi)
!!
!!   ! This routine writes a UNKnk.s used for plotting
!!
!!   implicit none
!!
!!   ! I/O variables
!!   integer, intent(in) :: n_bands, nx, ny, nz, nk, s
!!   real(kind=8), dimension(nx,ny,nz,n_bands), intent(in) :: psi
!!
!!   ! Local variables
!!   integer :: nb, i, j, k
!!   character :: s_c*1, nk_c*3, seedname*10
!!
!!
!!   ! Concatenations to write the names of UNKnk.s files.
!!   if (nk>0 .and. nk<10) then 
!!      write(nk_c, '(i1)') nk
!!      write(s_c, '(i1)') s
!!      seedname=(('UNK0000'//trim(nk_c))//'.')//s_c
!!   else  
!!      if (nk>9 .and. nk<100) then
!!         write(nk_c, '(i2)') nk
!!         write(s_c, '(i1)') s
!!         seedname=(('UNK000'//trim(nk_c))//'.')//s_c
!!      else
!!         if (nk>99 .and. nk<1000) then
!!            write(nk_c, '(i3)') nk
!!            write(s_c, '(i1)') s
!!            seedname=(('UNK00'//trim(nk_c))//'.')//s_c
!!         else
!!            if (nk>999 .and. nk<10000) then
!!               write(nk_c, '(i4)') nk
!!               write(s_c, '(i1)') s
!!               seedname=(('UNK0'//trim(nk_c))//'.')//s_c
!!            end if
!!         end if
!!      end if
!!   end if
!!
!!   ! Writing the UNKnk.s file
!!   OPEN(12, FILE=seedname, STATUS='unknown')
!!   write(*,*) '!==================================!'
!!   write(*,*) '!     Writing a UNKnk.s file :     !'
!!   write(*,*) '!==================================!'
!!   write(12,'(I4,4(1X,I4))') nx, ny, nz, nk, n_bands
!!   do nb=1, n_bands
!!      do k=1, nz
!!         do j=1, ny
!!            do i=1, nx
!!               write(12,'(E13.6, 1X, E13.6)') psi(i,j,k,nb), 0.d0
!!            end do
!!         end do
!!      end do
!!   end do
!!   CLOSE(12)
!!   write(*,*) '!==================================!'
!!   write(*,*) '!    Writing a UNKnk.s file done   !'
!!   write(*,*) '!==================================!'
!!   write(*,*)
!!   write(*,*)
!!
!!END SUBROUTINE write_unk
