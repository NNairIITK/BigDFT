!> this routine gives the covalent radius of a atomi provided its symbol in
!! standard format (capitals at the beginning)
subroutine covalent_radius(symbol,rcov)
  implicit none
  character(len=*), intent(in) :: symbol  !< Atomic symbol
  double precision, intent(out) :: rcov   !< Covalent radius in AU

  select case(trim(adjustl(symbol)))
  case( "LJ") ; rcov=0.56d0
  case( "H") ; rcov=0.75d0
  case( "He"); rcov=0.75d0
  case( "Li"); rcov=3.40d0
  case( "Be"); rcov=2.30d0
  case( "B") ; rcov=1.55d0
  case( "C") ; rcov=1.45d0
  case( "N") ; rcov=1.42d0
  case( "O") ; rcov=1.38d0
  case( "F") ; rcov=1.35d0
  case( "Ne"); rcov=1.35d0
  case( "Na"); rcov=3.40d0
  case( "Mg"); rcov=2.65d0
  case( "Al"); rcov=2.23d0
  case( "Si"); rcov=2.09d0
  case( "P") ; rcov=2.00d0
  case( "S") ; rcov=1.92d0
  case( "Cl"); rcov=1.87d0
  case( "Ar"); rcov=1.80d0
  case( "K") ; rcov=4.00d0
  case( "Ca"); rcov=3.80d0
  case( "Sc"); rcov=2.70d0
  case( "Ti"); rcov=2.70d0
  case( "V") ; rcov=2.60d0
  case( "Cr"); rcov=2.60d0
  case( "Mn" ); rcov=2.50d0
  case( "Fe" ); rcov=2.50d0
  case( "Co" ); rcov=2.40d0
  case( "Ni" ); rcov=2.30d0
  case( "Cu" ); rcov=2.80d0
  case( "Zn" ); rcov=2.70d0
  case( "Ga" ); rcov=2.40d0
  case( "Ge" ); rcov=2.40d0
  case( "As" ); rcov=2.30d0
  case( "Se" ); rcov=2.30d0
  case( "Br" ); rcov=2.20d0
  case( "Kr" ); rcov=2.20d0
  case( "Rb" ); rcov=4.50d0
  case( "Sr" ); rcov=4.00d0
  case( "Y" ); rcov=3.50d0
  case( "Zr" ); rcov=3.00d0
  case( "Nb" ); rcov=2.92d0
  case( "Mo" ); rcov=2.83d0
  case( "Tc" ); rcov=2.75d0
  case( "Ru" ); rcov=2.67d0
  case( "Rh" ); rcov=2.58d0
  case( "Pd" ); rcov=2.50d0
  case( "Ag" ); rcov=2.90d0
  case( "Cd" ); rcov=2.80d0
  case( "In" ); rcov=2.70d0
  case( "Sn" ); rcov=2.66d0
  case( "Sb" ); rcov=2.66d0
  case( "Te" ); rcov=2.53d0
  case( "I" ); rcov=2.50d0
  case( "Xe" ); rcov=2.50d0
  case( "Cs" ); rcov=4.50d0
  case( "Ba" ); rcov=4.00d0
  case( "La" ); rcov=3.50d0
  case( "Ce" ); rcov=3.50d0
  case( "Pr" ); rcov=3.44d0
  case( "Nd" ); rcov=3.38d0
  case( "Pm" ); rcov=3.33d0
  case( "Sm" ); rcov=3.27d0
  case( "Eu" ); rcov=3.21d0
  case( "Gd" ); rcov=3.15d0
  case( "Tb" ); rcov=3.09d0
  case( "Dy" ); rcov=3.03d0
  case( "Ho" ); rcov=2.97d0
  case( "Er" ); rcov=2.92d0
  case( "Tm" ); rcov=2.92d0
  case( "Yb" ); rcov=2.80d0
  case( "Lu" ); rcov=2.80d0
  case( "Hf" ); rcov=2.90d0
  case( "Ta" ); rcov=3.10d0
  case( "W" ); rcov=2.60d0
  case( "Re" ); rcov=2.60d0
  case( "Os" ); rcov=2.50d0
  case( "Ir" ); rcov=2.60d0
  case( "Pt" ); rcov=2.60d0
  case( "Au" ); rcov=4.00d0
  case( "Hg" ); rcov=3.20d0
  case( "Tl" ); rcov=3.20d0
  case( "Pb" ); rcov=3.30d0
  case( "Bi" ); rcov=2.90d0
  case( "Po" ); rcov=2.80d0
  case( "At" ); rcov=2.60d0
  case( "Rn" ); rcov=2.60d0
  case( "U" ); rcov=3.38d0
  case default
     stop 'covalent radius cannot be determined, invalid symbol'
  end select
end subroutine covalent_radius
