interface
   subroutine read_atomic_positions(iproc,ifile,units,nat,ntypes,iatype,atomnames,rxyz)
     implicit none
     character(len=20), intent(in) :: units
     integer, intent(in) :: iproc,ifile,nat
     integer, intent(out) :: ntypes
     character(len=20), dimension(100), intent(out) :: atomnames
     integer, dimension(nat), intent(out) :: iatype
     real(kind=8), dimension(3,nat), intent(out) :: rxyz
   end subroutine read_atomic_positions
end interface
public :: read_atomic_positions

interface
   subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,occup,spinar)
     implicit none
     ! Arguments
     integer, intent(in) :: nelec,nspin,iproc,norb,norbu,norbd,iunit
     real(kind=8), intent(out) :: occup(norb),spinar(norb)
   end subroutine input_occup
end interface
public :: input_occup

interface
   subroutine read_system_variables(iproc,nproc,nat,ntypes,nspin,ncharge,mpol,atomnames,iatype,&
        psppar,radii_cf,npspcode,iasctype,nelpsp,nzatom,nelec,natsc,norb,norbu,norbd,norbp,iunit)
     implicit none
     integer, intent(in) :: iproc,nproc,nat,ntypes,nspin,ncharge,mpol
     integer, intent(out) :: nelec,natsc,norb,norbu,norbd,norbp,iunit
     character(len=20), dimension(ntypes), intent(in) :: atomnames
     integer, dimension(ntypes), intent(in) :: iatype
     integer, dimension(ntypes), intent(out) :: npspcode,iasctype,nelpsp,nzatom
     real(kind=8), dimension(ntypes,2), intent(out) :: radii_cf
     real(kind=8), dimension(0:4,0:6,ntypes), intent(out) :: psppar
   end subroutine read_system_variables
end interface
public :: read_system_variables

interface
   subroutine system_size(iproc,nat,ntypes,rxyz,radii_cf,crmult,frmult,hgrid,iatype,atomnames, &
        alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)
     ! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
     !and shifts the atoms such that their position is the most symmetric possible
     implicit none
     integer, intent(in) :: iproc,nat,ntypes
     real(kind=8), intent(in) :: hgrid,crmult,frmult
     character(len=20), dimension(ntypes), intent(in) :: atomnames
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(3,nat), intent(inout) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     integer, intent(out) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
     real(kind=8), intent(out) :: alat1,alat2,alat3
   end subroutine system_size
end interface
public :: system_size

interface
   subroutine reformatmywaves(iproc, norb, norbp, nat, &
        & hgrid_old, nvctr_c_old, nvctr_f_old, n1_old, n2_old, n3_old, rxyz_old, &
        & nseg_c_old, nseg_f_old, keyg_old, keyv_old, psi_old, &
        & hgrid, nvctr_c, nvctr_f, n1, n2, n3, rxyz, &
        & nseg_c, nseg_f, keyg, keyv, psi)
     implicit real(kind=8) (a-h,o-z)
     dimension :: rxyz(3,nat), rxyz_old(3,nat), center(3), center_old(3)
     dimension :: keyg_old(2, nseg_c_old + nseg_f_old), keyv_old(nseg_c_old + nseg_f_old)
     dimension :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
     dimension :: psi_old(nvctr_c_old + 7 * nvctr_f_old, norbp), psi(nvctr_c + 7 * nvctr_f, norbp)
   end subroutine reformatmywaves
end interface
public :: reformatmywaves

interface
   subroutine readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,nat,rxyz,  & 
        nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
     ! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
     implicit real(kind=8) (a-h,o-z)
     dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
     dimension psi(nvctr_c+7*nvctr_f,norbp)
     dimension rxyz(3,nat),eval(norb),center(3)
   end subroutine readmywaves
end interface
public :: readmywaves

interface
   subroutine writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
        nat,rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,eval)
     ! write all my wavefunctions in files by calling writeonewave
     implicit real(kind=8) (a-h,o-z)
     character(len=4) f4
     character(len=50) filename
     dimension rxyz(3,nat),eval(norb),center(3)
     dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
     dimension psi(nvctr_c+7*nvctr_f,norbp)
   end subroutine writemywaves
end interface
public :: writemywaves


