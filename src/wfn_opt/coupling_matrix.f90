!> @file
!!  Routines related to coupling matrix (TD-DFT Casida's formalism) modified by MM
!! @author
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine center_of_charge(at,rxyz,cc)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3), intent(out) :: cc
  !local variables
  integer :: iat,ityp
  real(gp) :: zatom,rx,ry,rz,qtot

  cc(1)=0.0_gp
  cc(2)=0.0_gp
  cc(3)=0.0_gp
  qtot=0.0_gp
  do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     zatom=real(at%nelpsp(ityp),gp)
     qtot=qtot+zatom
     !coordinates of the center
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     cc(1)=cc(1)+rx*zatom
     cc(2)=cc(2)+ry*zatom
     cc(3)=cc(3)+rz*zatom
  end do
  if (qtot /= 0.0_gp) then
     cc(1)=cc(1)/qtot
     cc(2)=cc(2)/qtot
     cc(3)=cc(3)/qtot
  end if
END SUBROUTINE center_of_charge


!> Calculate the coupling matrix needed for Casida's TDDFT approach
subroutine coupling_matrix_prelim(iproc,nproc,geocode,tddft_approach,nspin,lr,orbsocc,orbsvirt,i3s,n3p,&
     hxh,hyh,hzh,chargec,pkernel,dvxcdrho,psirocc,psivirtr)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  use bounds, only: ext_buffers
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=4), intent(in) :: tddft_approach
  integer, intent(in) :: iproc,nproc,n3p,nspin,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3) :: chargec
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p*orbsocc%norb), intent(in) :: psirocc
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb), intent(in) :: psivirtr
  real(wp), dimension(lr%d%n1i,lr%d%n2i,n3p,max((nspin*(nspin+1))/2,2)), intent(in) :: dvxcdrho
  type(coulomb_operator) :: pkernel
  !local variables
  character(len=*), parameter :: subname='coupling_matrix_prelim'
  !logical :: onlyfxc=.false.,dofxc=.true.,perx,pery,perz
  logical :: dofH, dofxc, perx, pery, perz
  !logical :: tda=.true.
  integer :: imulti,jmulti,jorba,jorbi,spinindex
  integer :: i1,i2,i3p,iorbi,iorba,indi,inda,ind2,ind3,ntda,ispin,jspin
  integer :: ik,jk,nmulti,lwork,info,nbl1,nbl2,nbl3,nbr3,nbr2,nbr1,ndipoles
  real(gp) :: ehart,hfac,x,y,z,fsumrule_test
  real(wp), dimension(:), allocatable :: omega,work
  real(wp), dimension(:,:), allocatable :: K,Kbig,Kaux,dipoles,fi
  real(wp), dimension(:,:,:), allocatable :: v_ias
  real(wp), dimension(:,:,:,:), allocatable :: rho_ias

  if(iproc==0) call yaml_comment('Linear-Response TDDFT calculations',hfill='-')
  !if(iproc==0) write(*,'(1x,a)')"=========================================================="
  !if(iproc==0) write(*,'(1x,a)')" Linear-Response TDDFT calculations"


  !write(*,*) "dofxc=",dofxc,";  dofH=",dofH
  dofxc=.true.!.true.
  dofH=.true.!onlyfxc=.false.
  !write(*,*) "dofxc=",dofxc,";  dofH=",dofH

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)


  !Calculate nmulti, the number of allowed transitions between an occupied KS state and a virtual KS state.
  !Both states must have the same spin for the transition to occur.
  !For that purpose, we use a multi-index imulti: imulti=a+(i-1)*nvirt, where:
  !a is the virtual state index, i is the occupied state index, and nvirt is the total number of virtual states.
  !The length of this multi-index is the number of occupied states times the number of virtual states (nocc*nvirt)
  nmulti=0 !initialize the counter of allowed transitions
  do imulti=1,orbsvirt%norb*orbsocc%norb
     !calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1 !index of the occupied state considered
     iorba=imulti-(iorbi-1)*orbsvirt%norb !index of the virtual state considered
     !check the spin of both states is the same
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        !if same spin, then increment the counter of allowed transitions
        nmulti=nmulti+1
     end if
  end do

  !Allocate partial densities and potentials
  rho_ias = f_malloc((/ lr%d%n1i, lr%d%n2i, n3p, nmulti /),id='rho_ias')
  v_ias = f_malloc((/ lr%d%n1i, lr%d%n2i, n3p /),id='v_ias')

  !factor used during 3D integration
  hfac=1.0_gp/(hxh*hyh*hzh)

  !Choice between doing the Full TDDFT or using the Tamm-Damcoff Approximation (TDA)
  if (tddft_approach == 'full') then !if full-TDDFT
     !We want to solve \Omega F = \omega F with \Omega being the following matrix:
     !\Omega(i,j) = {\omega_i}^2 \delta_{i,j} + 2*sqrt(\omega_i \omega_j)*K(i,j)
     !where K(i,j) is the coupling matrix, coupling two transitions i and j.
     !Each transition is between an occupied and a virtual state of same spin (their energy is noted E_occ and E_virt, respectively).
     !The energy of the transition is then \omega_i = E_virt - E_occ.
     !\delta_{i,j} is the Kronecker symbol.
     !We are interested in the vector \omega, which gives the excitation energies of the system.
 
     !Is the calculation done spin averaged (nspin==1) or spin-polarized (nspin==2) ?
     !This plays a role when we compute the fxc part.
     if (nspin == 1) then !spin-averaged
        !In the spin-averaged case, the coupling matrix K is made of 4 symmetric sub-matrices of size nmulti*nmulti.
        !The K_H terms are the same in each of these 4 sub-matrices.
        !However, the exchange correlation part depends on the spin,
        !so that for the diagonal sub-matrices, one definition of f_xc is used, 
        !while another one is for the off-diagonal sub-matrices.

        !We have to take care of the fact that orbitals have a number of occupation, 
        !and therefore double the size of the coupling matrix.
        ndipoles = 2*nmulti

        !allocate the dipoles (computed in order to get the oscillator strength)
        dipoles = f_malloc0((/ 3, ndipoles /),id='dipoles')
        !allocate coupling matrix
        K = f_malloc0((/ ndipoles, ndipoles /),id='K')

        !Now we can start to build the partial densities and the corresponding partial potentials.
        !loop_i is a loop over a multiindex, the same as the one used above to find the number of allowed transitions.
        ik=0
        loop_i: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !occupied orbital index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !virtual orbital index

           !check the spin of the occupied and virtual orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              !if they have the same spin, then increment the counter of transitions
              ik=ik+1
              !if (orbsocc%spinsgn(iorbi) == 1.0_gp) then
              !   ispin=1
              !else if (orbsocc%spinsgn(iorbi) == -1.0_gp) then
              !   ispin=2
              !end if
           else
              !if the orbitals do not have the same spin, then cycle
              cycle loop_i
           end if


           !Partial density definition (psi_occ*psi_virt) and dipole moments calculation.
           !We use three loops, one for each direction of space.
           !Another multi-index is used here for the position of the grid point.
           do i3p=1,n3p !loop over z
              ind3=(i3p-1)*lr%d%n1i*lr%d%n2i !z-index
              z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3) !z (origin at the center of charge)

              do i2=1,lr%d%n2i !loop over y
                 ind2=(i2-1)*lr%d%n1i+ind3 !y-index
                 y=real(i2-nbl2-1,gp)*hyh-chargec(2) !y (origin at the center of charge)

                 do i1=1,lr%d%n1i !loop over x
                    x=real(i1-nbl1-1,gp)*hxh-chargec(1) !x (origin at the center of charge)
                    
                    indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the occupied state wavefunction
                    inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the virtual state wavefunction
                    rho_ias(i1,i2,i3p,ik)=hfac*psivirtr(inda)*psirocc(indi) !partial density, multiplied by the integration factor hfac

                    !Calculate dipole moments
                    dipoles(1,ik)=dipoles(1,ik)+x*rho_ias(i1,i2,i3p,ik) !dipole moment along x-axis: \int dx dy dz \psi_occ(x,y,z) x \psi_virt(x,y,z)
                    dipoles(2,ik)=dipoles(2,ik)+y*rho_ias(i1,i2,i3p,ik) !dipole moment along y-axis: \int dx dy dz \psi_occ(x,y,z) y \psi_virt(x,y,z)
                    dipoles(3,ik)=dipoles(3,ik)+z*rho_ias(i1,i2,i3p,ik) !dipole moment along z-axis: \int dx dy dz \psi_occ(x,y,z) z \psi_virt(x,y,z)
                 end do
              end do
           end do
 
           !Calculate the Hartree potential corresponding to the partial density
           if (dofH) then
              !Copy the partial density onto the partial potential space to pass it to PSolver
              call vcopy(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,ik),1,v_ias(1,1,1),1)
              !Partial potential term for each partial density
              call H_potential('D',pkernel,v_ias(1,1,1),rho_ias,ehart,0.0_dp,.false.,&
                   quiet='YES')
           end if

           !After the Poisson Solver we can calculate the lower triangular part of two sub-matrices of the coupling matrix
           jk=0 !initialize another counter of transition
           loop_j: do jmulti=1,imulti 
              !Calculate the orbital index for the second transition of the coupling matrix element
              !We use the same multi index as above for the first transtion
              jorbi=(jmulti-1)/orbsvirt%norb+1 !index of the occupied state
              jorba=jmulti-(jorbi-1)*orbsvirt%norb !index of the virtual state

              !Check the spin of the occupied and virtual orbitals considered (needed for the exchange correlation part)
              !This is the same convention as above. 
              if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
                 jk=jk+1
                 !if (orbsocc%spinsgn(jorbi) == 1.0_gp) then
                 !   jspin=1
                 !else if (orbsocc%spinsgn(jorbi) == -1.0_gp) then
                 !   jspin=2
                 !end if
              else
                 cycle loop_j
              end if

              !Calculation of the Hartree part of the coupling matrix K_H(ik,jk)
              !K_H(ik,jk) = \int V_Hartree(x,y,z) \rho_bj(x,y,z) dx dy dz
              if (dofH) then
                 K(ik,jk)=hxh*hyh*hzh*&
                      dot(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,jk),1,v_ias(1,1,1),1)
                 !Copy the Hartree term in one of the off-diagonal sub-matrix.
                 K(ik+nmulti,jk)=K(ik,jk)
              end if

              !Add the XC contribution
              !Map the spin couples to the index of dvxcdrho, in the abinit convention
              if (dofxc) then
                 !Add the XC contribution
                 do i3p=1,n3p
                    do i2=1,lr%d%n2i
                       do i1=1,lr%d%n1i
                          !Add the XC contribution to one of the diagonal matrix
                          K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                               rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                               dvxcdrho(i1,i2,i3p,1) !index=1
                          !Add the XC contribution to one of the off-diagonal matrix
                          K(ik+nmulti,jk)=K(ik+nmulti,jk)+hxh*hyh*hzh*&
                               rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                               dvxcdrho(i1,i2,i3p,2) !index=2
                       end do
                    end do
                 end do

                 !!Add the XC contribution to one of the off-diagonal matrix
                 !do i3p=1,n3p
                 !   do i2=1,lr%d%n2i
                 !      do i1=1,lr%d%n1i
                 !         K(ik+nmulti,jk)=K(ik+nmulti,jk)+hxh*hyh*hzh*&
                 !              rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !              dvxcdrho(i1,i2,i3p,2) !index=2
                 !      end do
                 !   end do
                 !end do
              end if

              !Multiply K by the factor 2*\sqrt(\omega_i \omega_j)
              K(ik       ,jk)=K(ik       ,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
                   sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))
              K(ik+nmulti,jk)=K(ik+nmulti,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
                   sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))

           end do loop_j
        end do loop_i

        !We now build the lower triangular parts of the two sub-matrices.
        do imulti=1,nmulti
           do jmulti=imulti+1,nmulti
              K(imulti,jmulti)=K(jmulti,imulti)
              K(imulti+nmulti,jmulti)=K(jmulti+nmulti,imulti)
           end do
        end do

        !Add the diagonal part: {\omega_i}^2 \delta_{i,j}
        !loop over the transitions
        ik=0 !initialize the transition counter
        loop_i2: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !occupied orbital index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !virtual orbital index
           !check the spin of the orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              ik=ik+1
           else
              cycle loop_i2
           end if
           !Add the square of the energy difference of the eigenvalues
           K(ik,ik)=K(ik,ik)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**2
           !K(ik+nmulti,ik+nmulti)=K(ik+nmulti,ik+nmulti)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**2
        end do loop_i2

        !We finally build the two other submatrices
        do imulti=1,nmulti
           do jmulti=1,nmulti
              K(imulti,jmulti+nmulti)=K(imulti+nmulti,jmulti)
              K(imulti+nmulti,jmulti+nmulti)=K(imulti,jmulti)
           end do
        end do

        !If more than one processor, then perform two MPI_all_reduce
        if (nproc > 1) then
           call mpiallred(K,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the coupling matrix
           call mpiallred(dipoles(1,1),3*nmulti,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the dipoles
        end if
      
        !Copy the values of the dipoles in the second part of the array
        call vcopy(3*nmulti,dipoles(1,1),1,dipoles(1,nmulti+1),1)
        call dscal(3*ndipoles,hxh*hyh*hzh,dipoles(1,1),1)

        call f_free(v_ias)

        !Allocate the oscillator strength
        fi = f_malloc((/ 3, ndipoles /),id='fi')
        !Allocate the excitation energy vector (actually \omega^2 in this case)
        omega = f_malloc(ndipoles,id='omega')

        lwork = 3*ndipoles !safe value
        work = f_malloc(lwork,id='work')

        !!test: print out the matrix elements of K
        !do jmulti = 1, ndipoles
        !   fsumrule_test=0.0
        !   do imulti = 1, ndipoles
        !      write(*,*) jmulti, imulti, K(jmulti, imulti)
        !      fsumrule_test=fsumrule_test+K(jmulti, imulti)**2
        !   end do
        !   write(*,*) jmulti, fsumrule_test
        !end do

        call DSYEV('V','U',ndipoles,K,ndipoles,omega,work,lwork,info)
        if (info /= 0) then
           call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
           !write(*,*) 'Error, DSYEV',info
        end if

        !Calculation of the transition dipoles (K now represents the coefficients of the KS transition expansion of the true excitations)
        call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
             K(1,1),ndipoles,0.0_wp,fi(1,1),3)

        !Summary of the results and pretty printing
        if (iproc == 0) then

           !if (tddft_approach=='TDA') call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')
           !if (tddft_approach=='full') call yaml_comment('FULL TDDFT',hfill='-')
           call yaml_comment('FULL TDDFT',hfill='-')

           call yaml_sequence_open('Excitation Energy and Oscillator Strength')

           do imulti = 1, ndipoles
              call yaml_sequence(trim(yaml_toa((/ Ha_eV*sqrt(omega(imulti)),&
                 sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2) /),&
                 fmt='(1pe10.3)')),advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
           end do

           call yaml_sequence_close()


           !Test of the Thomas-Reiche-Kuhn sum rule
           !This is not what is actually coded
           !fsumrule_test=0_wp
           !do imulti = 1, ndipoles
           !   fsumrule_test = fsumrule_test + &
           !                 sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           !end do
           !call yaml_map('Test of the f-sum rule',fsumrule_test,fmt='(1pe10.3)')

           !Write an output file containing excitation energies and oscillator strength.
           !It will be used to plot absorption spectr.a
           open(unit=9, file='td_spectra.txt')
           !write(9,'(i4,5x,a19)') ndipoles, '#(results in eV)' 
           write(9,'(i4)') ndipoles
           !do imulti = 1, min(100, ndipoles) 
           do imulti = 1, ndipoles 
              write(9,'(f9.4,5x,1pe10.3)') Ha_eV*sqrt(omega(imulti)),&
                   sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           close(unit=9)

           !Count the number of KS transitions needed to account for all true excitations
           ik=0
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then !We chose a minimal value for the transition to be taken into account
                       ik=ik+1
                    end if
                    if (abs(K(jmulti+nmulti,imulti)) > 5.d-02) then
                       ik=ik+1
                    end if
                 end do
              end do
           end do

           !Write a file containing all data concerning the transitions.
           !This will be used to plot further outputs.
           open(unit=10, file='transitions.txt')
           !The first line contains the number of transitions that will be written
           write(10,*) ik
           !Then we loop ever the matrix elements of K 
           !(K is now containing the coefficients of the KS transitions reproducing the true excitations.)
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then
                       !iorbi=jmulti/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb
                       write(10,*) Ha_eV*sqrt(omega(imulti)), iorbi, orbsocc%eval(iorbi),&
                              &iorba, orbsvirt%eval(iorba), abs(K(jmulti,imulti)),&
                              &sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                    end if
                    if (abs(K(jmulti+nmulti,imulti)) > 5.d-02) then
                    write(10,*) Ha_eV*sqrt(omega(imulti)), iorbi, orbsocc%eval(iorbi),&
                              &iorba, orbsvirt%eval(iorba), abs(K(jmulti+nmulti,imulti)),&
                              &sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                    end if
                 end do
              end do
           end do
           close(unit=10)

           !Pretty printing of the transition energies in the output
           call yaml_sequence_open('Transition energies (eV)')
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti

              call yaml_sequence(advance='no')
              call yaml_sequence_open(advance='no',flow=.true.)
              !if (tddft_approach=='TDA') call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
              !if (tddft_approach=='full') call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))
              call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))

              ik=0
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then
                       !iorbi=jmulti/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_mapping_open(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(K(jmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_mapping_close()   
                    end if
                    if (abs(K(jmulti+nmulti,imulti)) > 5.d-02) then
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_mapping_open(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(K(jmulti+nmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_mapping_close()
                    end if
                 end do
              end do

              call yaml_sequence_close(advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))

           end do
           call yaml_sequence_close()

        end if !iproc == 0

        call f_free(omega)
        call f_free(work)
        call f_free(fi)
        call f_free(rho_ias)
        call f_free(K)
        call f_free(dipoles)


     else if (nspin == 2) then !spin-polarized
        !In the spin-polarized case, the coupling matrix K is made is a matrix of size nmulti*nmulti.
        !There is no spin dependence on Hartree part of the coupling matrix K_H
        !However, the exchange correlation part depends on the spin,
        !so that three definitions of the exchange correlation are used :
        !one for f_xc^{up,up}, one for f_xc^{down,down} and one for f_xc^{down,up}==f_xc^{up,down}.

        !The number of dipoles in this case is the same as the number of allowed transitions in this case
        !as the occupation number of each KS orbital is one.
        ndipoles = nmulti

        !allocate the dipoles (computed in order to get the oscillator strength)
        dipoles = f_malloc0((/ 3, ndipoles /),id='dipoles')
        !allocate coupling matrix
        K = f_malloc0((/ ndipoles, ndipoles /),id='K')

        !Now we can start to build the partial densities and the corresponding partial potentials.
        !loop_i is a loop over a multiindex, the same as the one used above to find the number of allowed transitions.
        ik=0
        loop_i3: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !occupied state index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !virtual state index

           !check the spin of the occupied and virtual orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              !if they have the same spin, then increment the counter of transitions
              ik=ik+1
              !if (orbsocc%spinsgn(iorbi) == 1.0_gp) then
              !   ispin=1
              !else if (orbsocc%spinsgn(iorbi) == -1.0_gp) then
              !   ispin=2
              !end if
           else
              !if the orbitals do not have the same spin, then cycle
              cycle loop_i3
           end if

           !Partial density definition (psi_occ*psi_virt) and dipole moments calculation.
           !We use three loops, one for each direction of space.
           !Another multi-index is used here for the position of the grid point.
           do i3p=1,n3p !loop over z
              ind3=(i3p-1)*lr%d%n1i*lr%d%n2i !z-index
              z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3) !z (origin at the center of charge)

              do i2=1,lr%d%n2i !loop over y
                 ind2=(i2-1)*lr%d%n1i+ind3 !y-index
                 y=real(i2-nbl2-1,gp)*hyh-chargec(2) !y (origin at the center of charge)

                 do i1=1,lr%d%n1i !loop over x
                    x=real(i1-nbl1-1,gp)*hxh-chargec(1) !x (origin at the center of charge)
                    
                    indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the occupied state wavefunction
                    inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the virtual state wavefunction
                    rho_ias(i1,i2,i3p,ik)=hfac*psivirtr(inda)*psirocc(indi) !partial density, multiplied by the integration factor hfac

                    !Calculate dipole moments
                    dipoles(1,ik)=dipoles(1,ik)+x*rho_ias(i1,i2,i3p,ik) !dipole moment along x-axis: \int dx dy dz \psi_occ(x,y,z) x \psi_virt(x,y,z)
                    dipoles(2,ik)=dipoles(2,ik)+y*rho_ias(i1,i2,i3p,ik) !dipole moment along y-axis: \int dx dy dz \psi_occ(x,y,z) y \psi_virt(x,y,z)
                    dipoles(3,ik)=dipoles(3,ik)+z*rho_ias(i1,i2,i3p,ik) !dipole moment along z-axis: \int dx dy dz \psi_occ(x,y,z) z \psi_virt(x,y,z)
                 end do
              end do
           end do

           !Calculate the Hartree potential corresponding to the partial density
           if (dofH) then
              !Copy the partial density onto the partial potential space to pass it to PSolver
              call vcopy(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,ik),1,v_ias(1,1,1),1)
              !Partial potential term for each partial density
              call H_potential('D',pkernel,v_ias(1,1,1),rho_ias,ehart,0.0_dp,.false.,&
                   quiet='YES')
           end if

           !After the Poisson Solver we can calculate the lower triangular of the coupling matrix
           jk=0 !initialize another counter of transition
           loop_j2: do jmulti=1,imulti 
              !Calculate the orbital index for the second transition of the coupling matrix element
              !We use the same multi index as above for the first transtion
              jorbi=(jmulti-1)/orbsvirt%norb+1 !index of the occupied state
              jorba=jmulti-(jorbi-1)*orbsvirt%norb !index of the virtual state

              !Check the spin of the occupied and virtual orbitals considered (needed for the exchange correlation part)
              !This is the same convention as above. 
              if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
                 jk=jk+1
                 !if (orbsocc%spinsgn(jorbi) == 1.0_gp) then
                 !   jspin=1
                 !else if (orbsocc%spinsgn(jorbi) == -1.0_gp) then
                 !   jspin=2
                 !end if
              else
                 cycle loop_j2
              end if

              !Calculation of the Hartree part of the coupling matrix K_H(ik,jk)
              !K_H(ik,jk) = \int V_Hartree(x,y,z) \rho_bj(x,y,z) dx dy dz
              if (dofH) then
                 K(ik,jk)=hxh*hyh*hzh*&
                      dot(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,jk),1,v_ias(1,1,1),1)
              end if

              !Add the XC contribution
              !Map the spin couples to the index of dvxcdrho, in the abinit convention
              if (dofxc) then

                 !fxc depends on the spin sign of the two transitions.
                 !We define the right spin contribution using the variable spin-index.
                 !spinindex=1: up-up
                 !spinindex=2: up-down and down-up
                 !spinindex=3: down-down
                 spinindex=ispin+jspin-1

                 !Add the right XC contribution
                 do i3p=1,n3p
                    do i2=1,lr%d%n2i
                       do i1=1,lr%d%n1i
                          K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                               rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                               dvxcdrho(i1,i2,i3p,spinindex) !index=1, 2 or 3
                       end do
                    end do
                 end do

                 !!Add the right XC contribution (according to the the spins of the transitions)
                 !if ( orbsocc%spinsgn(iorbi) == 1.0_gp .and. orbsocc%spinsgn(jorbi) == 1.0_gp ) then

                 !   do i3p=1,n3p
                 !      do i2=1,lr%d%n2i
                 !         do i1=1,lr%d%n1i
                 !            K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                 !                 rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !                 dvxcdrho(i1,i2,i3p,1) !index=1
                 !         end do
                 !      end do
                 !   end do

                 !else if ( (orbsocc%spinsgn(iorbi) == 1.0_gp .and. orbsocc%spinsgn(jorbi) == -1.0_gp) &
                 !     .or. (orbsocc%spinsgn(iorbi) == -1.0_gp .and. orbsocc%spinsgn(jorbi) == 1.0_gp) ) then

                 !   do i3p=1,n3p
                 !      do i2=1,lr%d%n2i
                 !         do i1=1,lr%d%n1i
                 !            K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                 !                 rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !                 dvxcdrho(i1,i2,i3p,2) !index=2
                 !         end do
                 !      end do
                 !   end do

                 !else

                 !   do i3p=1,n3p
                 !      do i2=1,lr%d%n2i
                 !         do i1=1,lr%d%n1i
                 !            K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                 !                 rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !                 dvxcdrho(i1,i2,i3p,3) !index=3
                 !         end do
                 !      end do
                 !   end do

                 !end if

              end if

              !Multiply K by the factor 2*\sqrt(\omega_i \omega_j)
              K(ik,jk)=K(ik,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
                   sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))

           end do loop_j2
        end do loop_i3

        !We now build the lower triangular parts of the two sub-matrices.
        do imulti=1,nmulti
           do jmulti=imulti+1,nmulti
              K(imulti,jmulti)=K(jmulti,imulti)
           end do
        end do

        !Add the diagonal part: {\omega_i}^2 \delta_{i,j}
        !loop over the transitions
        ik=0 !initialize the transition counter
        loop_i4: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !occupied orbital index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !virtual orbital index
           !check the spin of the orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              ik=ik+1
           else
              cycle loop_i4
           end if
           !Add the square of the energy difference of the eigenvalues
           K(ik,ik)=K(ik,ik)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**2
           !K(ik+nmulti,ik+nmulti)=K(ik+nmulti,ik+nmulti)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**2
        end do loop_i4

        !If more than one processor, then perform two MPI_all_reduce
        if (nproc > 1) then
           call mpiallred(K,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the coupling matrix
           call mpiallred(dipoles(1,1),3*nmulti,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the dipoles
        end if
      
        !Copy the values of the dipoles in the second part of the array
        !call vcopy(3*nmulti,dipoles(1,1),1,dipoles(1,nmulti+1),1)
        call dscal(3*ndipoles,hxh*hyh*hzh,dipoles(1,1),1)

        call f_free(v_ias)

        !Allocate the oscillator strength
        fi = f_malloc((/ 3, ndipoles /),id='fi')
        !Allocate the excitation energy vector (actually \omega^2 in this case)
        omega = f_malloc(ndipoles,id='omega')

        lwork = 3*ndipoles !safe value
        work = f_malloc(lwork,id='work')

        !!test: print out the matrix elements of K
        !do jmulti = 1, ndipoles
        !   fsumrule_test=0.0
        !   do imulti = 1, ndipoles
        !      write(*,*) jmulti, imulti, K(jmulti, imulti)
        !      fsumrule_test=fsumrule_test+K(jmulti, imulti)**2
        !   end do
        !   write(*,*) jmulti, fsumrule_test
        !end do

        call DSYEV('V','U',ndipoles,K,ndipoles,omega,work,lwork,info)
        if (info /= 0) then
           call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
           !write(*,*) 'Error, DSYEV',info
        end if

        !Calculation of the transition dipoles (K now represents the coefficients of the KS transition expansion of the true excitations)
        call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
             K(1,1),ndipoles,0.0_wp,fi(1,1),3)

        !Summary of the results and pretty printing
        if (iproc == 0) then

           !if (tddft_approach=='TDA') call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')
           !if (tddft_approach=='full') call yaml_comment('FULL TDDFT',hfill='-')
           call yaml_comment('FULL TDDFT',hfill='-')

           call yaml_sequence_open('Excitation Energy and Oscillator Strength')

           do imulti = 1, ndipoles
              call yaml_sequence(trim(yaml_toa((/ Ha_eV*sqrt(omega(imulti)),&
                 sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2) /),&
                 fmt='(1pe10.3)')),advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
           end do

           call yaml_sequence_close()


           !Test of the Thomas-Reiche-Kuhn sum rule
           !This is not what is actually coded
           !fsumrule_test=0_wp
           !do imulti = 1, ndipoles
           !   fsumrule_test = fsumrule_test + &
           !                 sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           !end do
           !call yaml_map('Test of the f-sum rule',fsumrule_test,fmt='(1pe10.3)')

           !Write an output file containing excitation energies and oscillator strength.
           !It will be used to plot absorption spectr.a
           open(unit=9, file='td_spectra.txt')
           !write(9,'(i4,5x,a19)') ndipoles, '#(results in eV)' 
           write(9,'(i4)') ndipoles
           !do imulti = 1, min(100, ndipoles) 
           do imulti = 1, ndipoles 
              write(9,'(f9.4,5x,1pe10.3)') Ha_eV*sqrt(omega(imulti)),&
                   sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           close(unit=9)

           !Count the number of KS transitions needed to account for all true excitations
           ik=0
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then !We chose a minimal value for the transition to be taken into account
                       ik=ik+1
                    end if
                 end do
              end do
           end do

           !Write a file containing all data concerning the transitions.
           !This will be used to plot further outputs.
           open(unit=10, file='transitions.txt')
           !The first line contains the number of transitions that will be written
           write(10,*) ik
           !Then we loop ever the matrix elements of K 
           !(K is now containing the coefficients of the KS transitions reproducing the true excitations.)
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then
                       !iorbi=jmulti/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb
                       write(10,*) Ha_eV*sqrt(omega(imulti)), iorbi, orbsocc%eval(iorbi),&
                              &iorba, orbsvirt%eval(iorba), abs(K(jmulti,imulti)),&
                              &sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                    end if
                 end do
              end do
           end do
           close(unit=10)

           !Pretty printing of the transition energies in the output
           call yaml_sequence_open('Transition energies (eV)')
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti

              call yaml_sequence(advance='no')
              call yaml_sequence_open(advance='no',flow=.true.)
              !if (tddft_approach=='TDA') call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
              !if (tddft_approach=='full') call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))
              call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))

              ik=0
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then
                       !iorbi=jmulti/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_mapping_open(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(K(jmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_mapping_close()   
                    end if
                 end do
              end do

              call yaml_sequence_close(advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))

           end do
           call yaml_sequence_close()

        end if !iproc == 0

        call f_free(omega)
        call f_free(work)
        call f_free(fi)
        call f_free(rho_ias)
        call f_free(K)
        call f_free(dipoles)

     end if !nspin==2

  else if (tddft_approach == 'TDA') then !if Tamm-Dancoff Approximation

     !!exponent: 0 for Tamm-Dancoff Approximation
     !ntda=0
 
     !Is the calculation done spin averaged (nspin==1) or spin-polarized (nspin==2) ?
     !This plays a role when we compute the fxc part.
     if (nspin == 1) then !spin-averaged
        !In the spin-averaged case, the coupling matrix K is made of 4 symmetric sub-matrices of size nmulti*nmulti.
        !The K_H terms are the same in each of these 4 sub-matrices.
        !However, the exchange correlation part depends on the spin,
        !so that for the diagonal sub-matrices, one definition of f_xc is used, 
        !while another one is for the off-diagonal sub-matrices.

        !We have to take care of the fact that orbitals have a number of occupation, 
        !and therefore double the size of the coupling matrix.
        ndipoles = 2*nmulti

        !allocate the dipoles (computed in order to get the oscillator strength)
        dipoles = f_malloc0((/ 3, ndipoles /),id='dipoles')
        !allocate coupling matrix
        K = f_malloc0((/ ndipoles, ndipoles /),id='K')

        !Now we can start to build the partial densities and the corresponding partial potentials.
        !loop_i is a loop over a multiindex, the same as the one used above to find the number of allowed transitions.
        ik=0
        loop_i5: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !occupied orbital index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !virtual orbital index

           !check the spin of the occupied and virtual orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              !if they have the same spin, then increment the counter of transitions
              ik=ik+1
           else
              !if the orbitals do not have the same spin, then cycle
              cycle loop_i5
           end if


           !Partial density definition (psi_occ*psi_virt) and dipole moments calculation.
           !We use three loops, one for each direction of space.
           !Another multi-index is used here for the position of the grid point.
           do i3p=1,n3p !loop over z
              ind3=(i3p-1)*lr%d%n1i*lr%d%n2i !z-index
              z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3) !z (origin at the center of charge)

              do i2=1,lr%d%n2i !loop over y
                 ind2=(i2-1)*lr%d%n1i+ind3 !y-index
                 y=real(i2-nbl2-1,gp)*hyh-chargec(2) !y (origin at the center of charge)

                 do i1=1,lr%d%n1i !loop over x
                    x=real(i1-nbl1-1,gp)*hxh-chargec(1) !x (origin at the center of charge)
                    
                    indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the occupied state wavefunction
                    inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the virtual state wavefunction
                    rho_ias(i1,i2,i3p,ik)=hfac*psivirtr(inda)*psirocc(indi) !partial density, multiplied by the integration factor hfac

                    !Calculate dipole moments
                    dipoles(1,ik)=dipoles(1,ik)+x*rho_ias(i1,i2,i3p,ik) !dipole moment along x-axis: \int dx dy dz \psi_occ(x,y,z) x \psi_virt(x,y,z)
                    dipoles(2,ik)=dipoles(2,ik)+y*rho_ias(i1,i2,i3p,ik) !dipole moment along y-axis: \int dx dy dz \psi_occ(x,y,z) y \psi_virt(x,y,z)
                    dipoles(3,ik)=dipoles(3,ik)+z*rho_ias(i1,i2,i3p,ik) !dipole moment along z-axis: \int dx dy dz \psi_occ(x,y,z) z \psi_virt(x,y,z)
                 end do
              end do
           end do
 
           !Calculate the Hartree potential corresponding to the partial density
           if (dofH) then
              !Copy the partial density onto the partial potential space to pass it to PSolver
              call vcopy(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,ik),1,v_ias(1,1,1),1)
              !Partial potential term for each partial density
              call H_potential('D',pkernel,v_ias(1,1,1),rho_ias,ehart,0.0_dp,.false.,&
                   quiet='YES')
           end if

           !After the Poisson Solver we can calculate the lower triangular part of two sub-matrices of the coupling matrix
           jk=0 !initialize another counter of transition
           loop_j3: do jmulti=1,imulti 
              !Calculate the orbital index for the second transition of the coupling matrix element
              !We use the same multi index as above for the first transtion
              jorbi=(jmulti-1)/orbsvirt%norb+1 !index of the occupied state
              jorba=jmulti-(jorbi-1)*orbsvirt%norb !index of the virtual state

              !Check the spin of the occupied and virtual orbitals considered (needed for the exchange correlation part)
              !This is the same convention as above. 
              if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
                 jk=jk+1
                 !if (orbsocc%spinsgn(jorbi) == 1.0_gp) then
                 !   jspin=1
                 !else if (orbsocc%spinsgn(jorbi) == -1.0_gp) then
                 !   jspin=2
                 !end if
              else
                 cycle loop_j3
              end if

              !Calculation of the Hartree part of the coupling matrix K_H(ik,jk)
              !K_H(ik,jk) = \int V_Hartree(x,y,z) \rho_bj(x,y,z) dx dy dz
              if (dofH) then
                 K(ik,jk)=hxh*hyh*hzh*&
                      dot(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,jk),1,v_ias(1,1,1),1)
                 !Copy the Hartree term in one of the off-diagonal sub-matrix.
                 K(ik+nmulti,jk)=K(ik,jk)
              end if

              !Add the XC contribution
              !Map the spin couples to the index of dvxcdrho, in the abinit convention
              if (dofxc) then
                 do i3p=1,n3p
                    do i2=1,lr%d%n2i
                       do i1=1,lr%d%n1i
                          !Add the XC contribution to one of the diagonal matrix
                          K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                               rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                               dvxcdrho(i1,i2,i3p,1) !index=1
                          !Add the XC contribution to one of the off-diagonal matrix
                          K(ik+nmulti,jk)=K(ik+nmulti,jk)+hxh*hyh*hzh*&
                               rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                               dvxcdrho(i1,i2,i3p,2) !index=2
                       end do
                    end do
                 end do

                 !!Add the XC contribution to one of the off-diagonal matrix
                 !do i3p=1,n3p
                 !   do i2=1,lr%d%n2i
                 !      do i1=1,lr%d%n1i
                 !         K(ik+nmulti,jk)=K(ik+nmulti,jk)+hxh*hyh*hzh*&
                 !              rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !              dvxcdrho(i1,i2,i3p,2) !index=2
                 !      end do
                 !   end do
                 !end do
              end if

              !!Multiply K by the factor 2*\sqrt(\omega_i \omega_j)
              !K(ik       ,jk)=K(ik       ,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
              !     sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))
              !K(ik+nmulti,jk)=K(ik+nmulti,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
              !     sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))

           end do loop_j3
        end do loop_i5

        !If more than one processor, then perform two MPI_all_reduce
        if (nproc > 1) then
           call mpiallred(K,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the coupling matrix
           call mpiallred(dipoles(1,1),3*nmulti,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the dipoles
        end if
      
        !Add the diagonal part: \omega_i \delta_{i,j}
        !loop over the transitions
        ik=0 !initialize the transition counter
        loop_i6: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !virtual orbital index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !ocupied orbital index
           !check the spin of the orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              ik=ik+1
           else
              cycle loop_i6
           end if
           !Add the energy difference of the eigenvalues
           K(ik,ik)=K(ik,ik)+orbsvirt%eval(iorba)-orbsocc%eval(iorbi)
           !K(ik+nmulti,ik+nmulti)=K(ik,ik)
        end do loop_i6

        !We now build the lower triangular parts of the two sub-matrices.
        do imulti=1,nmulti
           do jmulti=imulti+1,nmulti
              K(imulti,jmulti)=K(jmulti,imulti)
              K(imulti+nmulti,jmulti)=K(jmulti+nmulti,imulti)
           end do
        end do

        !We finally build the two other submatrices
        do imulti=1,nmulti
           do jmulti=1,nmulti
              K(imulti,jmulti+nmulti)=K(imulti+nmulti,jmulti)
              K(imulti+nmulti,jmulti+nmulti)=K(imulti,jmulti)
           end do
        end do

        !Copy the values of the dipoles in the second part of the array
        call vcopy(3*nmulti,dipoles(1,1),1,dipoles(1,nmulti+1),1)
        call dscal(3*ndipoles,hxh*hyh*hzh,dipoles(1,1),1)

        call f_free(v_ias)

        !Allocate the oscillator strength
        fi = f_malloc((/ 3, ndipoles /),id='fi')
        !Allocate the excitation energy vector (actually \omega^2 in this case)
        omega = f_malloc(ndipoles,id='omega')

        lwork = 3*ndipoles !safe value
        work = f_malloc(lwork,id='work')

        !test: print out the matrix elements of K
        do jmulti = 1, ndipoles
           fsumrule_test=0.0
           do imulti = 1, ndipoles
              write(*,*) jmulti, imulti, K(jmulti, imulti)
              fsumrule_test=fsumrule_test+K(jmulti, imulti)**2
           end do
           write(*,*) jmulti, fsumrule_test
        end do

        call DSYEV('V','U',ndipoles,K,ndipoles,omega,work,lwork,info)
        if (info /= 0) then
           call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
           !write(*,*) 'Error, DSYEV',info
        end if

        !Calculation of the transition dipoles (K now represents the coefficients of the KS transition expansion of the true excitations)
        call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
             K(1,1),ndipoles,0.0_wp,fi(1,1),3)

        !Summary of the results and pretty printing
        if (iproc == 0) then

           !if (tddft_approach=='TDA') call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')
           !if (tddft_approach=='full') call yaml_comment('FULL TDDFT',hfill='-')
           call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')

           call yaml_sequence_open('Excitation Energy and Oscillator Strength')

           do imulti = 1, ndipoles
              call yaml_sequence(trim(yaml_toa((/ Ha_eV*omega(imulti),&
                 omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2) /),&
                 fmt='(1pe10.3)')),advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
           end do

           call yaml_sequence_close()


           !Test of the Thomas-Reiche-Kuhn sum rule
           !This is not what is actually coded ???
           fsumrule_test=0_wp
           do imulti = 1, ndipoles
              fsumrule_test = fsumrule_test + &
                            omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           call yaml_map('Test of the f-sum rule',fsumrule_test,fmt='(1pe10.3)')

           !Write an output file containing excitation energies and oscillator strength.
           !It will be used to plot absorption spectr.a
           open(unit=9, file='td_spectra.txt')
           !write(9,'(i4,5x,a19)') ndipoles, '#(results in eV)' 
           write(9,'(i4)') ndipoles
           !do imulti = 1, min(100, ndipoles) 
           do imulti = 1, ndipoles 
              write(9,'(f9.4,5x,1pe10.3)') Ha_eV*omega(imulti),&
                   omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           close(unit=9)

           !Count the number of KS transitions needed to account for all true excitations
           ik=0
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then !We chose a minimal value for the transition to be taken into account
                       ik=ik+1
                    end if
                    if (abs(K(jmulti+nmulti,imulti)) > 5.d-02) then
                       ik=ik+1
                    end if
                 end do
              end do
           end do

           !Write a file containing all data concerning the transitions.
           !This will be used to plot further outputs.
           open(unit=10, file='transitions.txt')
           !The first line contains the number of transitions that will be written
           write(10,*) ik
           !Then we loop ever the matrix elements of K 
           !(K is now containing the coefficients of the KS transitions reproducing the true excitations.)
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti = (iorbi-1)*orbsvirt%norb + iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then
                       !iorbi=jmulti/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb
                       write(10,*) Ha_eV*omega(imulti), iorbi, orbsocc%eval(iorbi),&
                              &iorba, orbsvirt%eval(iorba), abs(K(jmulti,imulti)),&
                              &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                    end if
                    if (abs(K(jmulti+nmulti,imulti)) > 5.d-02) then
                    write(10,*) Ha_eV*omega(imulti), iorbi, orbsocc%eval(iorbi),&
                              &iorba, orbsvirt%eval(iorba), abs(K(jmulti+nmulti,imulti)),&
                              &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                    end if
                 end do
              end do
           end do
           close(unit=10)

           !Pretty printing of the transition energies in the output
           call yaml_sequence_open('Transition energies (eV)')
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti

              call yaml_sequence(advance='no')
              call yaml_sequence_open(advance='no',flow=.true.)
              !if (tddft_approach=='TDA') call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
              !if (tddft_approach=='full') call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))
              call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))

              ik=0
              !do jmulti = 1, ndipoles
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(K(jmulti,imulti)) > 5.d-02) then
                       !iorbi=jmulti/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_mapping_open(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(K(jmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_mapping_close()   
                    end if
                    if (abs(K(jmulti+nmulti,imulti)) > 5.d-02) then
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_mapping_open(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(K(jmulti+nmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_mapping_close()
                    end if
                 end do
              end do

              call yaml_sequence_close(advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))

           end do
           call yaml_sequence_close()

        end if !iproc == 0

        call f_free(omega)
        call f_free(work)
        call f_free(fi)
        call f_free(rho_ias)
        call f_free(K)
        call f_free(dipoles)


     else if (nspin == 2) then !spin-polarized
        !In the spin-polarized case, the coupling matrix K is made is a matrix of size nmulti*nmulti.
        !There is no spin dependence on Hartree part of the coupling matrix K_H
        !However, the exchange correlation part depends on the spin,
        !so that three definitions of the exchange correlation are used :
        !one for f_xc^{up,up}, one for f_xc^{down,down} and one for f_xc^{down,up}==f_xc^{up,down}.

        !The number of dipoles in this case is the same as the number of allowed transitions in this case
        !as the occupation number of each KS orbital is one.
        ndipoles = nmulti

        !allocate the dipoles (computed in order to get the oscillator strength)
        dipoles = f_malloc0((/ 3, ndipoles /),id='dipoles')
        !allocate coupling matrix
        K = f_malloc0((/ ndipoles, ndipoles /),id='K')

        !Now we can start to build the partial densities and the corresponding partial potentials.
        !loop_i is a loop over a multiindex, the same as the one used above to find the number of allowed transitions.
        ik=0
        loop_i7: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !index of the occupied state involved in the first transition of the coupling matrix element
           iorba=imulti-(iorbi-1)*orbsvirt%norb !index of the virtual state involved in the first transition of the coupling matrix element

           !check the spin of the occupied and virtual orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              !if they have the same spin, then increment the counter of transitions
              ik=ik+1
              !We also memorize the spin sign of the KS states involved in the first transition of the coupling matrix element
              if (orbsocc%spinsgn(iorbi) == 1.0_gp) then
                 ispin=1
              else if (orbsocc%spinsgn(iorbi) == -1.0_gp) then
                 ispin=2
              end if
           else
              !if the orbitals do not have the same spin, then cycle
              cycle loop_i7
           end if

           !Partial density definition (psi_occ*psi_virt) and dipole moments calculation.
           !We use three loops, one for each direction of space.
           !Another multi-index is used here for the position of the grid point.
           do i3p=1,n3p !loop over z
              ind3=(i3p-1)*lr%d%n1i*lr%d%n2i !z-index
              z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3) !z (origin at the center of charge)

              do i2=1,lr%d%n2i !loop over y
                 ind2=(i2-1)*lr%d%n1i+ind3 !y-index
                 y=real(i2-nbl2-1,gp)*hyh-chargec(2) !y (origin at the center of charge)

                 do i1=1,lr%d%n1i !loop over x
                    x=real(i1-nbl1-1,gp)*hxh-chargec(1) !x (origin at the center of charge)
                    
                    indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the occupied state wavefunction
                    inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the virtual state wavefunction
                    rho_ias(i1,i2,i3p,ik)=hfac*psivirtr(inda)*psirocc(indi) !partial density, multiplied by the integration factor hfac

                    !Calculate dipole moments
                    dipoles(1,ik)=dipoles(1,ik)+x*rho_ias(i1,i2,i3p,ik) !dipole moment along x-axis: \int dx dy dz \psi_occ(x,y,z) x \psi_virt(x,y,z)
                    dipoles(2,ik)=dipoles(2,ik)+y*rho_ias(i1,i2,i3p,ik) !dipole moment along y-axis: \int dx dy dz \psi_occ(x,y,z) y \psi_virt(x,y,z)
                    dipoles(3,ik)=dipoles(3,ik)+z*rho_ias(i1,i2,i3p,ik) !dipole moment along z-axis: \int dx dy dz \psi_occ(x,y,z) z \psi_virt(x,y,z)
                 end do
              end do
           end do

           !Calculate the Hartree potential corresponding to the partial density
           if (dofH) then
              !Copy the partial density onto the partial potential space to pass it to PSolver
              call vcopy(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,ik),1,v_ias(1,1,1),1)
              !Partial potential term for each partial density
              call H_potential('D',pkernel,v_ias(1,1,1),rho_ias,ehart,0.0_dp,.false.,&
                   quiet='YES')
           end if

           !After the Poisson Solver we can calculate the lower triangular of the coupling matrix
           jk=0 !initialize another counter of transition
           loop_j4: do jmulti=1,imulti 
              !Calculate the orbital index for the second transition of the coupling matrix element
              !We use the same multi index as above for the first transtion
              jorbi=(jmulti-1)/orbsvirt%norb+1 !index of the occupied state involved in the second transition of the coupling matrix element
              jorba=jmulti-(jorbi-1)*orbsvirt%norb !index of the virtual state involved in the second transition of the coupling matrix element

              !Check the spin of the occupied and virtual orbitals considered (needed for the exchange correlation part)
              !This is the same convention as above. 
              if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
                 jk=jk+1
                 !We also memorize the spin sign of the KS states involved in the second transition of the coupling matrix element
                 if (orbsocc%spinsgn(jorbi) == 1.0_gp) then
                    jspin=1
                 else if (orbsocc%spinsgn(jorbi) == -1.0_gp) then
                    jspin=2
                 end if
              else
                 cycle loop_j4
              end if

              !Calculation of the Hartree part of the coupling matrix K_H(ik,jk)
              !K_H(ik,jk) = \int V_Hartree(x,y,z) \rho_bj(x,y,z) dx dy dz
              if (dofH) then
                 K(ik,jk)=hxh*hyh*hzh*&
                      dot(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,jk),1,v_ias(1,1,1),1)
              end if

              !Add the XC contribution
              !Map the spin couples to the index of dvxcdrho, in the abinit convention
              if (dofxc) then
                 
                 !fxc depends on the spin sign of the two transitions.
                 !We define the right spin contribution using the variable spin-index.
                 !spinindex=1: up-up
                 !spinindex=2: up-down and down-up
                 !spinindex=3: down-down
                 spinindex=ispin+jspin-1

                 !Add the right XC contribution
                 do i3p=1,n3p
                    do i2=1,lr%d%n2i
                       do i1=1,lr%d%n1i
                          K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                               rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                               dvxcdrho(i1,i2,i3p,spinindex) !index=1, 2 or 3
                       end do
                    end do
                 end do

                 !!Add the right XC contribution (according to the the spins of the transitions)
                 !if ( orbsocc%spinsgn(iorbi) == 1.0_gp .and. orbsocc%spinsgn(jorbi) == 1.0_gp ) then

                 !   do i3p=1,n3p
                 !      do i2=1,lr%d%n2i
                 !         do i1=1,lr%d%n1i
                 !            K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                 !                 rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !                 dvxcdrho(i1,i2,i3p,1) !index=1
                 !         end do
                 !      end do
                 !   end do

                 !else if ( (orbsocc%spinsgn(iorbi) == 1.0_gp .and. orbsocc%spinsgn(jorbi) == -1.0_gp) &
                 !     .or. (orbsocc%spinsgn(iorbi) == -1.0_gp .and. orbsocc%spinsgn(jorbi) == 1.0_gp) ) then

                 !   do i3p=1,n3p
                 !      do i2=1,lr%d%n2i
                 !         do i1=1,lr%d%n1i
                 !            K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                 !                 rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !                 dvxcdrho(i1,i2,i3p,2) !index=2
                 !         end do
                 !      end do
                 !   end do

                 !else

                 !   do i3p=1,n3p
                 !      do i2=1,lr%d%n2i
                 !         do i1=1,lr%d%n1i
                 !            K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                 !                 rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                 !                 dvxcdrho(i1,i2,i3p,3) !index=3
                 !         end do
                 !      end do
                 !   end do

                 !end if

              end if

           end do loop_j4
        end do loop_i7

        !If more than one processor, then perform two MPI_all_reduce
        if (nproc > 1) then
           call mpiallred(K,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the coupling matrix
           call mpiallred(dipoles(1,1),3*nmulti,MPI_SUM,comm=bigdft_mpi%mpi_comm) !MPI_all_reduce of the dipoles
        end if
      
        !Add the diagonal part: \omega_i \delta_{i,j}
        !loop over the transitions
        ik=0 !initialize the transition counter
        loop_i8: do imulti=1,orbsvirt%norb*orbsocc%norb
           !calculate the orbital index
           iorbi=(imulti-1)/orbsvirt%norb+1 !occupied orbital index
           iorba=imulti-(iorbi-1)*orbsvirt%norb !virtual orbital index
           !check the spin of the orbitals considered
           if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
              ik=ik+1
           else
              cycle loop_i8
           end if
           !Add the energy difference of the eigenvalues
           K(ik,ik)=K(ik,ik)+orbsvirt%eval(iorba)-orbsocc%eval(iorbi)
           !K(ik+nmulti,ik+nmulti)=K(ik+nmulti,ik+nmulti)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**2
        end do loop_i8

        !We now build the upper triangular part of the coupling matrix.
        do imulti=1,nmulti
           do jmulti=imulti+1,nmulti
              K(imulti,jmulti)=K(jmulti,imulti)
              !K(imulti+nmulti,jmulti)=K(jmulti+nmulti,imulti)
           end do
        end do

        !Copy the values of the dipoles in the second part of the array
        !call vcopy(3*nmulti,dipoles(1,1),1,dipoles(1,nmulti+1),1)
        call dscal(3*ndipoles,hxh*hyh*hzh,dipoles(1,1),1)

        call f_free(v_ias)

        !Allocate the oscillator strength
        fi = f_malloc((/ 3, ndipoles /),id='fi')
        !Allocate the excitation energy vector (actually \omega^2 in this case)
        omega = f_malloc(ndipoles,id='omega')

        lwork = 3*ndipoles !safe value
        work = f_malloc(lwork,id='work')

        !test: print out the matrix elements of K
        do jmulti = 1, ndipoles
           fsumrule_test=0.0
           do imulti = 1, ndipoles
              write(*,*) jmulti, imulti, K(jmulti, imulti)
              fsumrule_test=fsumrule_test+K(jmulti, imulti)**2
           end do
           write(*,*) jmulti, fsumrule_test
        end do

        call DSYEV('V','U',ndipoles,K,ndipoles,omega,work,lwork,info)
        if (info /= 0) then
           call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
           !write(*,*) 'Error, DSYEV',info
        end if

        !Calculation of the transition dipoles (K now represents the coefficients of the KS transition expansion of the true excitations)
        call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
             K(1,1),ndipoles,0.0_wp,fi(1,1),3)

        !Summary of the results and pretty printing
        if (iproc == 0) then

           !if (tddft_approach=='TDA') call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')
           !if (tddft_approach=='full') call yaml_comment('FULL TDDFT',hfill='-')
           call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')

           call yaml_sequence_open('Excitation Energy and Oscillator Strength')

           do imulti = 1, ndipoles
              call yaml_sequence(trim(yaml_toa((/ Ha_eV*omega(imulti),&
                 omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2) /),&
                 fmt='(1pe10.3)')),advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
           end do

           call yaml_sequence_close()


           !Test of the Thomas-Reiche-Kuhn sum rule
           !This is not what is actually coded ???
           fsumrule_test=0_wp
           do imulti = 1, ndipoles
              fsumrule_test = fsumrule_test + &
                            sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           call yaml_map('Test of the f-sum rule',fsumrule_test,fmt='(1pe10.3)')

           !Write an output file containing excitation energies and oscillator strength.
           !It will be used to plot absorption spectr.a
           open(unit=9, file='td_spectra.txt')
           !write(9,'(i4,5x,a19)') ndipoles, '#(results in eV)' 
           write(9,'(i4)') ndipoles
           !do imulti = 1, min(100, ndipoles) 
           do imulti = 1, ndipoles 
              write(9,'(f9.4,5x,1pe10.3)') Ha_eV*omega(imulti),&
                   omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           close(unit=9)

           !Count the number of KS transitions needed to account for all true excitations
           ik=0 !counter of transitions
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              jk=0 !counter of dipoles
              do jmulti = 1, orbsocc%norb*orbsvirt%norb
                 jorbi=(jmulti-1)/orbsvirt%norb+1
                 jorba=jmulti-(jorbi-1)*orbsvirt%norb  
                 if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
                    jk=jk+1
              !do jmulti = 1, ndipoles
              !do jorbi = 1, orbsocc%norb
              !   do jorba = 1, orbsvirt%norb
              !      jmulti =  (jorbi-1)*orbsvirt%norb+ jorba
              !      write (*,*) "iorbi", iorbi, "iorba", iorba, "jmulti", jmulti
                    if (abs(K(jk,imulti)) > 5.d-02) then !We chose a minimal value for the transition to be taken into account
                       ik=ik+1
                    end if
                 else
                    cycle
                 end if
              !   end do
              end do
           end do

           !Write a file containing all data concerning the transitions.
           !This will be used to plot further outputs.
           open(unit=10, file='transitions.txt')
           !The first line contains the number of transitions that will be written
           write(10,*) ik
           !Then we loop ever the matrix elements of K 
           !(K is now containing the coefficients of the KS transitions reproducing the true excitations.)
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti
              ik=0
              do jmulti = 1, orbsocc%norb*orbsvirt%norb
                 jorbi=(jmulti-1)/orbsvirt%norb+1
                 jorba=jmulti-(jorbi-1)*orbsvirt%norb  
                 if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
              !do jmulti = 1, ndipoles
              !do jorbi = 1, orbsocc%norb
              !   do jorba = 1, orbsvirt%norb
              !      jmulti =  (jorbi-1)*orbsvirt%norb + jorba
                    ik=ik+1
                    !write (*,*) "jorbi", jorbi, "jorba", jorba, "jmulti", jmulti, "ik", ik
                    if (abs(K(ik,imulti)) > 5.d-02) then
                       !iorbi=(jmulti-1)/orbsvirt%norb+1
                       !iorba=jmulti-(iorbi-1)*orbsvirt%norb  
                       !write (*,*) "jorbi", jorbi, "jorba", jorba, "jmulti", jmulti, "ik", ik
                       write(10,*) Ha_eV*omega(imulti), jorbi, orbsocc%eval(jorbi),&
                              &jorba, orbsvirt%eval(jorba), abs(K(ik,imulti)),&
                              &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                    end if
              !   end do
                 else 
                    cycle
                 end if
              end do
           end do
           close(unit=10)

           !Pretty printing of the transition energies in the output
           call yaml_sequence_open('Transition energies (eV)')
           do imulti = 1, ndipoles
           !do imulti = 1, nmulti

              call yaml_sequence(advance='no')
              call yaml_sequence_open(advance='no',flow=.true.)
              !if (tddft_approach=='TDA') call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
              !if (tddft_approach=='full') call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))
              call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))

              ik=0 !counter for the yaml_newline
              jk=0
              do jmulti = 1, orbsocc%norb*orbsvirt%norb
                 jorbi=(jmulti-1)/orbsvirt%norb+1
                 jorba=jmulti-(jorbi-1)*orbsvirt%norb
                 if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
              !do iorbi = 1, orbsocc%norb
              !   do iorba = 1, orbsvirt%norb
              !      jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    jk=jk+1 
                    if (abs(K(jk,imulti)) > 5.d-02) then
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_mapping_open(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ jorbi, jorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(K(jk,imulti)),fmt='(1pe10.3)')))
                       call yaml_mapping_close()   
                    end if
              !   end do
                 else
                    cycle
                 end if
              end do

              call yaml_sequence_close(advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))

           end do
           call yaml_sequence_close()

        end if !iproc == 0

        call f_free(omega)
        call f_free(work)
        call f_free(fi)
        call f_free(rho_ias)
        call f_free(K)
        call f_free(dipoles)

     end if !nspin==2

  end if


  !write(*,*) "dofxc=",dofxc,";  dofH=", dofH


END SUBROUTINE coupling_matrix_prelim
