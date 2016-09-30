module module_alborz
  implicit none

  private
  character(len=5), allocatable, save :: sat(:)
  character(len=4), save :: boundcond

  public initialize_alborz, finalize_alborz, call_to_alborz_get

  contains

  subroutine initialize_alborz(nat,astruct)
    use module_base
    use yaml_output
    use module_atoms
    !parameter
    implicit none
    integer, intent(in) :: nat
    type(atomic_structure), intent(in) :: astruct   
    !internal
    integer :: iat

     !!!sat = f_malloc((/nat/),id='sat')
     allocate(sat(nat))
     do iat=1,nat
         sat(iat)=astruct%atomnames(astruct%iatype(iat))
     enddo

    !Boundary Conditions
    select case(astruct%geocode)
    case('P')
      boundcond='bulk'
    case('S')
      boundcond='slab'
    case('F')
      boundcond='free'
    case('W')
      boundcond='wire'
    end select
     call alborz_as_potential_init(nat,sat)

  end subroutine initialize_alborz

  subroutine finalize_alborz
    implicit none
    call alborz_as_potential_final
    boundcond=''
    !!call f_free(sat)
    deallocate(sat)
  end subroutine finalize_alborz

  subroutine call_to_alborz_get(nat,cell,xcart,fcart,energy,strten)
    use module_base
      implicit none
      integer, intent(in):: nat
      real(8), intent(in) :: cell(3)
      real(8), intent(in):: xcart(3,nat)
      real(8), intent(out):: fcart(3,nat), energy, strten(6)
      integer:: iat
      real(8):: stress(3,3),latvec(3,3), vol

      !hack for cell usage for free BC in alborz
      real(8) :: cmx, cmy, cmz
      real(8) :: xcarts(3,nat)
      real(8) :: cellin(3)

      cellin=cell

    if(boundcond=='free')then
      ! center of mass
      cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
      do iat=1,nat
         cmx=cmx+xcart(1,iat)
         cmy=cmy+xcart(2,iat)
         cmz=cmz+xcart(3,iat)
      enddo
      cmx=cmx/real(nat,gp)
      cmy=cmy/real(nat,gp)
      cmz=cmz/real(nat,gp)
      cellin(1) = 2.d0*(maxval(xcart(1,:))-minval(xcart(1,:)))+3.d0
      cellin(2) = 2.d0*(maxval(xcart(2,:))-minval(xcart(2,:)))+3.d0
      cellin(3) = 2.d0*(maxval(xcart(3,:))-minval(xcart(3,:)))+3.d0
      do iat=1,nat
        xcarts(1,iat)=xcart(1,iat)-cmx + 0.5d0*cellin(1)
        xcarts(2,iat)=xcart(2,iat)-cmy + 0.5d0*cellin(2)
        xcarts(3,iat)=xcart(3,iat)-cmz + 0.5d0*cellin(3)
      enddo
    endif
  !hack end

      latvec(1,1)=cellin(1)
      latvec(2,1)=0.0_gp
      latvec(3,1)=0.0_gp
      latvec(1,2)=0.0_gp
      latvec(2,2)=cellin(2)
      latvec(3,2)=0.0_gp
      latvec(1,3)=0.0_gp
      latvec(2,3)=0.0_gp
      latvec(3,3)=cellin(3)
      call alborz_as_potential_get(boundcond,nat,latvec,xcarts,sat,fcart,energy,stress)
      strten(1) = -stress(1,1)
      strten(2) = -stress(2,2)
      strten(3) = -stress(3,3)
      strten(6) = -stress(2,1)
      strten(5) = -stress(3,1)
      strten(4) = -stress(3,2)
      call getvol_alborz(latvec,vol)
      strten=strten/vol
  end subroutine call_to_alborz_get
end module module_alborz

