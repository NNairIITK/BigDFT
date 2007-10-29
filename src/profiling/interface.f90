interface
   subroutine MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
        rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid,nspin)

     use Poisson_Solver

     implicit none
     !Arguments
     logical, intent(in) :: output_grid
     integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin
     integer, dimension(nat), intent(in) :: iatype
     character(len=20), dimension(100), intent(in) :: atomnames
     real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
   end subroutine MemoryEstimator
end interface
public :: MemoryEstimator
