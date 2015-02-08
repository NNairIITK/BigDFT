!> @file 
!!   Routines to do special convolution of linear toolbox
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS

subroutine ConvolQuartic4(iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
           hgrid, offsetx, offsety, offsetz, ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
           rxyzConf, potentialPrefac, with_kinetic, cprecr, maxdim, &
           xx_c, xx_f1, xx_f, xy_c, xy_f2, xy_f,  xz_c, xz_f4, xz_f, &
           aeff0array, beff0array, ceff0array, eeff0array, &
           aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array, &
           aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray, &
           xya_c, xyc_c, xza_c, xzc_c, &
           yza_c, yzc_c, xya_f, xyb_f, xyc_f, xye_f, &
           xza_f, xzb_f, xzc_f, xze_f, yza_f, yzb_f, yzc_f, yze_f, &
!           aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, &
!           ceff0, ceff1, ceff2, ceff3, eeff0, eeff1, eeff2, eeff3, &
!           aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2, &
!           ceff0_2, ceff1_2, ceff2_2, ceff3_2, eeff0_2, eeff1_2, eeff2_2, eeff3_2, & 
           y_c, y_f)

  use module_base, only: wp,gp,f_zero
  use dynamic_memory
  implicit none
  
  integer,parameter:: lb=-14 !< The lower bound for the filters.
  integer,parameter:: ub=14  !< The upper bound for the filters.

  real(8),dimension(lb:ub),parameter:: a = [ -6.92447490505951028E-18, &
       2.70800498995346639E-13, -5.81387993303650319E-11, -1.05857056453828591E-08, &
       -3.72307624729728559E-07, 2.09042354981647804E-06, -2.39822857110993937E-05, &
       4.51679195975884795E-04, -4.09765681251883507E-03,  2.20702923834323883E-02, &
       -8.22663977742195129E-02, 2.37178057432174683E-01, -6.15614175796508789E-01, &
       2.21914649009704590E+00, -3.55369234085083008E+00,  2.21914649009704590E+00, &
       -6.15614175796508789E-01, 2.37178057432174683E-01, -8.22663977742195129E-02, &
       2.20702923834323883E-02, -4.09765681251883507E-03,  4.51679195975884795E-04, &
       -2.39822857110993937E-05, 2.09042354981647804E-06, -3.72307624729728559E-07, &
       -1.05857056453828591E-08,-5.81387993303650319E-11,  2.70800498995346639E-13, &
       -6.92447490505951028E-18  ]

  real(8),dimension(lb:ub),parameter:: b = [  1.23926298748915144E-17, &
       -4.84665694980683149E-13, &
       1.04757348540220618E-10, &
       1.87909231522041173E-08, &
       6.39191314318015945E-07, &
       -4.62561549613067677E-06, &
       5.07580797039626444E-05, &
       -9.34236181144083580E-04, &
       8.08720235045630956E-03, &
       -5.65773658045679495E-02, &
       2.33704891436108609E-01, &
       -5.32982207596253462E-01, &
       6.34570315769713792E-01, &
       -2.17161162613989146E-01, &
       -4.20816528718038670E-01, &
       6.57030701554058305E-01, &
       -4.24980483446119562E-01, &
       1.43273284054659191E-01, &
       -2.84696145431842073E-02, &
       7.13761173990859820E-03, &
       -2.13562905732165138E-03, &
       2.17107937709352954E-04, &
       -1.13345668089136977E-05, &
       8.79697176379788384E-07, &
       -2.16568315068759430E-07, &
       -5.96264927673747788E-09, &
       -3.22647031147176524E-11, &
       1.51306166977329773E-13, &
       -3.86910248985843272E-18  ]

  real(8),dimension(lb:ub),parameter:: c = [ -3.86910248985843272E-18, &
       1.51306166977329773E-13, &
       -3.22647031147176524E-11, &
       -5.96264927673747788E-09, &
       -2.16568315068759430E-07, &
       8.79697176379788384E-07, &
       -1.13345668089136977E-05, &
       2.17107937709352954E-04, &
       -2.13562905732165138E-03, &
       7.13761173990859820E-03, &
       -2.84696145431842073E-02, &
       1.43273284054659191E-01, &
       -4.24980483446119562E-01, &
       6.57030701554058305E-01, &
       -4.20816528718038670E-01, &
       -2.17161162613989146E-01, &
       6.34570315769713792E-01, &
       -5.32982207596253462E-01, &
       2.33704891436108609E-01, &
       -5.65773658045679495E-02, &
       8.08720235045630956E-03, &
       -9.34236181144083580E-04, &
       5.07580797039626444E-05, &
       -4.62561549613067677E-06, &
       6.39191314318015945E-07, &
       1.87909231522041173E-08, &
       1.04757348540220618E-10, &
       -4.84665694980683149E-13, &
       1.23926298748915144E-17  ]

  real(8),dimension(lb:ub),parameter:: e = [  6.92447507792733360E-18, &
       -2.70800496076004922E-13, &
       5.81387997394137139E-11, &
       1.05857053487793953E-08, &
       3.72307646478137234E-07, &
       -2.09042376230958760E-06, &
       2.39822849010384008E-05, &
       -4.51679172658164203E-04, &
       4.09765614832326392E-03, &
       -2.20732688908272731E-02, &
       8.20745364532223690E-02, &
       -2.69959298589643554E-01, &
       -4.25170863143320842E-02, &
       -7.14405969106863975E+00, &
       -2.48758457220106488E+01, &
       -7.14405969106863975E+00, &
       -4.25170863143320842E-02, &
       -2.69959298589643554E-01, &
       8.20745364532223690E-02, &
       -2.20732688908272731E-02, &
       4.09765614832326392E-03, &
       -4.51679172658164203E-04, &
       2.39822849010384008E-05, &
       -2.09042376230958760E-06, &
       3.72307646478137234E-07, &
       1.05857053487793953E-08, &
       5.81387997394137139E-11, &
       -2.70800496076004922E-13, &
       6.92447507792733360E-18  ]

  real(8),dimension(lb:ub),parameter:: a1 = [ -4.27782380967095731E-22, &
       1.33836790789675412E-16, &
       -2.23077388809486687E-13, &
       -4.18701080751038290E-11, &
       3.84061421554449112E-10, &
       6.90677026682351425E-08, &
       -1.50667722209618660E-06, &
       1.43009565363172442E-05, &
       -7.38484450266696513E-05, &
       1.17437397420872003E-04, &
       4.36133099719882011E-03, &
       -2.15712171047925949E-02, &
       2.26678475737571716E-02, &
       9.18613970279693604E-02, &
       -3.43336127698421478E-02, &
       9.18613970279693604E-02, &
       2.26678475737571716E-02, &
       -2.15712171047925949E-02, &
       4.36133099719882011E-03, &
       1.17437397420872003E-04, &
       -7.38484450266696513E-05, &
       1.43009565363172442E-05, &
       -1.50667722209618660E-06, &
       6.90677026682351425E-08, &
       3.84061421554449112E-10, &
       -4.18701080751038290E-11, &
       -2.23077388809486687E-13, &
       1.33836790789675412E-16, &
       -4.27782380967095731E-22  ]

  real(8),dimension(lb:ub),parameter:: b1 = [  7.65595806716108416E-22, &
       -2.39526829029937495E-16, &
       3.99587744885038026E-13, &
       7.43506281339207227E-11, &
       -7.94750732307692066E-10, &
       -1.22260456570109001E-07, &
       2.87022347262966861E-06, &
       -3.01746671145040958E-05, &
       1.87517678388117517E-04, &
       -5.06159000963081146E-04, &
       -7.78069985307987407E-03, &
       5.60208708828413971E-02, &
       -1.59595684240634494E-01, &
       2.45149541954812523E-01, &
       -2.14819128613303770E-01, &
       9.90323037445533675E-02, &
       -1.30155486978821159E-02, &
       -6.91564874276195792E-03, &
       2.28647192633767860E-03, &
       3.94195663735359520E-06, &
       -2.62413561954799163E-05, &
       6.63758559364207982E-06, &
       -7.86922176158694678E-07, &
       3.89648346451391969E-08, &
       1.80799149218002746E-10, &
       -2.35773089429422988E-11, &
       -1.24537043995966775E-13, &
       7.47819644251633595E-17, &
       -2.39026636952344337E-22  ]

  real(8),dimension(lb:ub),parameter:: c1 = [ -2.39026636952344337E-22, &
       7.47819644251633595E-17, &
       -1.24537043995966775E-13, &
       -2.35773089429422988E-11, &
       1.80799149218002746E-10, &
       3.89648346451391969E-08, &
       -7.86922176158694678E-07, &
       6.63758559364207982E-06, &
       -2.62413561954799146E-05, &
       3.94195663735364602E-06, &
       2.28647192633767859E-03, &
       -6.91564874276195653E-03, &
       -1.30155486978821149E-02, &
       9.90323037445533658E-02, &
       -2.14819128613303770E-01, &
       2.45149541954812528E-01, &
       -1.59595684240634471E-01, &
       5.60208708828414089E-02, &
       -7.78069985307987367E-03, &
       -5.06159000963081133E-04, &
       1.87517678388117524E-04, &
       -3.01746671145040958E-05, &
       2.87022347262966861E-06, &
       -1.22260456570109001E-07, &
       -7.94750732307692066E-10, &
       7.43506281339207227E-11, &
       3.99587744885038026E-13, &
       -2.39526829029937495E-16, &
       7.65595806716108416E-22  ]

  real(8),dimension(lb:ub),parameter:: e1 = [  4.27782410746594903E-22, &
       -1.33836791682177454E-16, &
       2.23077373341138633E-13, &
       4.18701078507265792E-11, &
       -3.84061437780811023E-10, &
       -6.90677057885295201E-08, &
       1.50667719157836429E-06, &
       -1.43009565685797134E-05, &
       7.38484438588798723E-05, &
       -1.17437022042463083E-04, &
       -4.36283743232178177E-03, &
       2.14973685507261842E-02, &
       -1.83065170245152523E-02, &
       -6.91935470160085394E-02, &
       -1.11786975230691315E-09, &
       -6.91935470160085264E-02, &
       -1.83065170245152626E-02, &
       2.14973685507261885E-02, &
       -4.36283743232178238E-03, &
       -1.17437022042463063E-04, &
       7.38484438588798689E-05, &
       -1.43009565685797134E-05, &
       1.50667719157836429E-06, &
       -6.90677057885295201E-08, &
       -3.84061437780811023E-10, &
       4.18701078507265792E-11, &
       2.23077373341138633E-13, &
       -1.33836791682177454E-16, &
       4.27782410746594903E-22  ]

  real(8),dimension(lb:ub),parameter:: a2 = [  5.57613113797079647E-21, &
       -1.61540205057621221E-15, &
       2.44941827237665777E-12, &
       4.24629970074974494E-10, &
       -3.50622375577813727E-09, &
       -5.49912783753825352E-07, &
       1.07758069134433754E-05, &
       -8.87797723407857120E-05, &
       3.82019672542810440E-04, &
       -4.44838660769164562E-04, &
       -1.42199024558067322E-02, &
       3.99360917508602142E-02, &
       3.36528867483139038E-02, &
       -2.19314932823181152E-01, &
       1.65585339069366455E-01, &
       -3.55921462178230286E-02, &
       1.24324277043342590E-01, &
       -8.94912108778953552E-02, &
       2.06707436591386795E-02, &
       7.29535357095301151E-04, &
       -5.04161638673394918E-04, &
       1.11433619167655706E-04, &
       -1.33310277306009084E-05, &
       6.93305878485261928E-07, &
       4.17500478633314742E-09, &
       -4.96512386760628033E-10, &
       -2.90443884221058823E-12, &
       1.86435445701578929E-15, &
       -6.40177613495305920E-21  ]

  real(8),dimension(lb:ub),parameter:: b2 = [ -1.07451146172287577E-20, &
       3.13060246628943970E-15, &
       -4.78749340108647095E-12, &
       -8.27894411800580670E-10, &
       8.15738084296199976E-09, &
       1.09377804413390625E-06, &
       -2.35274921407857363E-05, &
       2.22483374802824865E-04, &
       -1.23566737058074802E-03, &
       2.81931962789967665E-03, &
       3.56570274706771702E-02, &
       -1.99870719176705833E-01, &
       4.09070789978144124E-01, &
       -3.83487854930014397E-01, &
       1.22706132164835730E-01, &
       4.07716120232467710E-02, &
       -1.73643908349863515E-02, &
       -1.68596657049951048E-02, &
       7.63552926096439174E-03, &
       5.38805252160585462E-05, &
       -1.33751751720543844E-04, &
       4.34601260064193826E-05, &
       -6.11241925587226617E-06, &
       3.52531960342273632E-07, &
       1.75161930985329523E-09, &
       -2.56223460315090577E-10, &
       -1.49681996733928602E-12, &
       9.66934874772265558E-16, &
       -3.33801397186598040E-21  ]

  real(8),dimension(lb:ub),parameter:: c2 = [  3.11570517857345129E-21, &
       -9.02614217267229018E-16, &
       1.36753205696331864E-12, &
       2.38900027287721328E-10, &
       -1.68356448485162769E-09, &
       -3.09870229248918553E-07, &
       5.69141338437313638E-06, &
       -4.28284866174613940E-05, &
       1.54903165637988905E-04, &
       1.84029190592224884E-05, &
       -8.36977423915950192E-03, &
       1.77185780796068396E-02, &
       2.16822550214745961E-02, &
       -5.82606911845323180E-02, &
       -9.21129972751010266E-02, &
       3.51960771815418264E-01, &
       -3.88907631878252118E-01, &
       1.92275377336910452E-01, &
       -3.43692713216321939E-02, &
       -2.74842935720907677E-03, &
       1.20206244492311618E-03, &
       -2.30136631629999615E-04, &
       2.52663068944202945E-05, &
       -1.22917063308546998E-06, &
       -8.53238440365876605E-09, &
       8.82170034527123439E-10, &
       5.20220026056775308E-12, &
       -3.33662185797648242E-15, &
       1.14571636261214844E-20  ]

  real(8),dimension(lb:ub),parameter:: e2 = [ -6.00391354600385591E-21, &
       1.74924014379757951E-15, &
       -2.67288615604453156E-12, &
       -4.65847082273536233E-10, &
       4.00909115868514137E-09, &
       6.17273186723220065E-07, &
       -1.24681741129386372E-05, &
       1.08808259548054617E-04, &
       -5.33848387091264561E-04, &
       6.96357808776253802E-04, &
       2.17704467254226292E-02, &
       -9.27586349639902768E-02, &
       1.17357770175067983E-01, &
       -1.44609721723709135E-02, &
       2.28332297558082387E-01, &
       -1.52848064755532954E-01, &
       4.41317011994163056E-02, &
       3.62255764829988165E-02, &
       -1.31322527750372515E-02, &
       -4.78012401544689226E-04, &
       3.52332937508960380E-04, &
       -9.14051322471805596E-05, &
       1.16386609504637462E-05, &
       -6.25945518694966712E-07, &
       -3.67213752393608442E-09, &
       4.55295290054135748E-10, &
       2.68097082627805916E-12, &
       -1.73051640666939972E-15, &
       5.97399387029546928E-21  ]

  real(8),dimension(lb:ub),parameter:: a3 = [ -5.45235389905989279E-20, &
       1.46358789074752665E-14, &
       -2.01811241329341584E-11, &
       -3.24547011487652526E-09, &
       2.42857414178843101E-08, &
       3.29380782204680145E-06, &
       -5.91745010751765221E-05, &
       4.36957285273820162E-04, &
       -1.67352124117314816E-03, &
       1.98640045709908009E-03, &
       4.50586304068565369E-02, &
       -6.27721697092056274E-02, &
       -1.87771692872047424E-01, &
       3.47782939672470093E-01, &
       -7.22572430968284607E-02, &
       -3.45776788890361786E-02, &
       2.86159813404083252E-01, &
       -2.85770207643508911E-01, &
       8.37636739015579224E-02, &
       4.12162579596042633E-03, &
       -2.77279899455606937E-03, &
       6.74822716973721981E-04, &
       -8.98371581570245326E-05, &
       5.22961454407777637E-06, &
       3.43174519912281539E-08, &
       -4.43153069795698684E-09, &
       -2.83714943899449068E-11, &
       1.94904500389536314E-14, &
       -7.18620847350200121E-20  ]

  real(8),dimension(lb:ub),parameter:: b3 = [  1.13123436610365484E-19, &
       -3.07100166700629325E-14, &
       4.30376242080603376E-11, &
       6.94162274675861743E-09, &
       -6.34036165941134668E-08, &
       -7.35379995218091881E-06, &
       1.47173214489243831E-04, &
       -1.28034683673851311E-03, &
       6.59497106305552553E-03, &
       -1.36921541340832689E-02, &
       -1.46897049595819335E-01, &
       6.37478980242085363E-01, &
       -9.30556215252272467E-01, &
       5.59458823389911942E-01, &
       -1.72942680845484709E-01, &
       1.62694799746086264E-01, &
       -1.03342326611766857E-01, &
       -2.08251898321085647E-02, &
       2.27089101547702390E-02, &
       8.34163098143173926E-04, &
       -5.63993287410049933E-04, &
       2.23449248084776174E-04, &
       -3.63125559271239384E-05, &
       2.39828253857973364E-06, &
       1.28421664625344509E-08, &
       -2.09729277533997787E-09, &
       -1.34982620213705742E-11, &
       9.38388852480617887E-15, &
       -3.49671886957245740E-20  ]

  real(8),dimension(lb:ub),parameter:: c3 = [ -3.04654359861049767E-20, &
       8.17787588062863309E-15, &
       -1.12680455877343732E-11, &
       -1.82444868673102661E-09, &
       1.18723849727281025E-08, &
       1.85434546705993772E-06, &
       -3.15762398492331581E-05, &
       2.17290763713708732E-04, &
       -7.38492447647641646E-04, &
       3.46249835019798076E-04, &
       2.65636962851717066E-02, &
       -2.40461111279288945E-02, &
       -1.13057521285956212E-01, &
       1.75811609409945607E-01, &
       -1.49997829824147128E-01, &
       4.88522886262723873E-01, &
       -8.54944372023607445E-01, &
       5.97603435176744222E-01, &
       -1.38204695319924360E-01, &
       -1.31073095034589249E-02, &
       6.26732305736694859E-03, &
       -1.36644597735980215E-03, &
       1.69343102660099443E-04, &
       -9.28314434210991049E-06, &
       -6.93099235960035602E-08, &
       7.87787729799378351E-09, &
       5.08133777611095321E-11, &
       -3.48819103637274518E-14, &
       1.28610501418155119E-19  ]

  real(8),dimension(lb:ub),parameter:: e3 = [  6.32085703717368256E-20, &
       -1.71593847930409723E-14, &
       2.40296297678319361E-11, &
       3.90306169567011521E-09, &
       -3.16398262659731984E-08, &
       -4.14682388687797182E-06, &
       7.86772014015926421E-05, &
       -6.42360316804161536E-04, &
       3.03629694962606823E-03, &
       -4.23534856061406406E-03, &
       -9.01621008293589939E-02, &
       3.08186486078434389E-01, &
       -3.02648997344974718E-01, &
       -1.63291833641518359E-03, &
       2.40612326154245674E-02, &
       -2.52596470969653632E-01, &
       1.81819415713290495E-01, &
       5.37877230562302454E-02, &
       -3.83329370344230677E-02, &
       -2.59775806368091250E-03, &
       1.40265791479232433E-03, &
       -4.59627480916031660E-04, &
       6.87230434889209774E-05, &
       -4.26390036678194569E-06, &
       -2.65855222668025459E-08, &
       3.72895715741983696E-09, &
       2.41751537799131406E-11, &
       -1.67942725094584736E-14, &
       6.25802566107744820E-20  ]

  real(8),dimension(lb:ub),parameter:: a4 = [  4.73982297509915827E-19, &
       -1.17970778003122223E-13, &
       1.47867218469599493E-10, &
       2.21538325462233843E-08, &
       -1.50909258422871062E-07, &
       -1.75753048097249120E-05, &
       2.94731522444635630E-04, &
       -2.00139172375202179E-03, &
       7.15284887701272964E-03, &
       -9.75563097745180130E-03, &
       -1.39017373323440552E-01, &
       3.75488102436065674E-02, &
       6.29527032375335693E-01, &
       -6.22455954551696777E-01, &
       2.39198490977287292E-01, &
       -1.79768219590187073E-01, &
       6.60393893718719482E-01, &
       -8.88859748840332031E-01, &
       3.33310693502426147E-01, &
       2.19652801752090454E-02, &
       -1.43004665151238441E-02, &
       3.75307234935462475E-03, &
       -5.46617549844086170E-04, &
       3.51455855707172304E-05, &
       2.53031799957170733E-07, &
       -3.52819533588899503E-08, &
       -2.46440173823359032E-10, &
       1.81234964102133800E-13, &
       -7.17145273459679160E-19  ]

  real(8),dimension(lb:ub),parameter:: b4 = [ -1.05879107559698705E-18, &
       2.67975668991362477E-13, &
       -3.44037604990392807E-10, &
       -5.19398195117999355E-08, &
       4.41615034358723889E-07, &
       4.40071230383561964E-05, &
       -8.30762190462945196E-04, &
       6.77564167923146284E-03, &
       -3.32405637216197889E-02, &
       6.54637771233801294E-02, &
       5.95097451648797403E-01, &
       -2.00500735905367301E+00, &
       2.09342664719743266E+00, &
       -8.21319709469621505E-01, &
       1.49823989402031084E-01, &
       1.42193968646241480E-01, &
       -2.58174510470954888E-01, &
       -6.34201291534517743E-03, &
       6.66580611749571090E-02, &
       6.79364923822744610E-03, &
       -2.24574587498619634E-03, &
       1.06412843446368771E-03, &
       -1.95665138116746604E-04, &
       1.45486779922684776E-05, &
       8.45902475713471526E-08, &
       -1.53245377007307429E-08, &
       -1.08248520437343012E-10, &
       8.10102415572423511E-14, &
       -3.25649385965431673E-19  ]

  real(8),dimension(lb:ub),parameter:: c4 = [  2.64841177354852000E-19, &
       -6.59168226006457738E-14, &
       8.25659985103387980E-11, &
       1.24446685547390328E-08, &
       -7.49608949097184224E-08, &
       -9.88743988180540206E-06, &
       1.58701212643154945E-04, &
       -1.01979783801389721E-03, &
       3.34978597853062363E-03, &
       -3.11164570844213000E-03, &
       -8.21852160064285005E-02, &
       1.90046894767111568E-03, &
       3.03170077824383044E-01, &
       -1.71554362097906985E-01, &
       -1.19411738947435464E-01, &
       6.67866005617359376E-01, &
       -1.84671114749237439E+00, &
       1.83677986289889223E+00, &
       -5.52785700365708630E-01, &
       -6.09059165980459183E-02, &
       3.09751752376641818E-02, &
       -7.46637518078680060E-03, &
       1.02466322093907858E-03, &
       -6.24487003129975085E-05, &
       -5.05322701667792061E-07, &
       6.27526586982342632E-08, &
       4.41347806410255367E-10, &
       -3.24354727486844695E-13, &
       1.28346416737257737E-18  ]

  real(8),dimension(lb:ub),parameter:: e4 = [ -5.91607470707858413E-19, &
       1.49732874622216788E-13, &
       -1.92100738129978935E-10, &
       -2.91841511215122459E-08, &
       2.23394316572768974E-07, &
       2.48001476795504746E-05, &
       -4.47594374637332827E-04, &
       3.47128223750824498E-03, &
       -1.60161771714504936E-02, &
       2.38282684039843004E-02, &
       3.61599000680735382E-01, &
       -9.82482562143621237E-01, &
       7.52753810635772963E-01, &
       -1.14204696883252737E-01, &
       2.51501726273655956E-01, &
       -4.84276384425144080E-01, &
       5.62339757401709483E-01, &
       2.85047897727701947E-02, &
       -1.07918110020031040E-01, &
       -1.51435426900326290E-02, &
       5.34875354394726576E-03, &
       -2.14609073676606404E-03, &
       3.67972099703273657E-04, &
       -2.58921736961675653E-05, &
       -1.72989790542752860E-07, &
       2.72620366731214978E-08, &
       1.93858659954877273E-10, &
       -1.44983357227351254E-13, &
       5.82809854008927484E-19  ]


  ! Calling arguments
  integer,intent(in) :: iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz, maxdim
  real(gp),intent(in) :: hgrid, potentialPrefac, cprecr
  logical,intent(in) :: with_kinetic
  real(8),dimension(3) :: rxyzConf
  integer,dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer,dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer,dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp),dimension(0:n1,0:n2,0:n3),intent(in) :: xx_c
  real(wp),dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in) :: xx_f1
  real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in) :: xx_f
  real(wp),dimension(0:n2,0:n1,0:n3),intent(in) :: xy_c
  real(wp),dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in) :: xy_f2
  real(wp),dimension(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in) :: xy_f
  real(wp),dimension(0:n3,0:n1,0:n2),intent(in) :: xz_c
  real(wp),dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in) :: xz_f4
  real(wp),dimension(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in) :: xz_f
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: aeff0array
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: beff0array
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: ceff0array
  real(wp),dimension(-14:14,0:maxdim),intent(inout):: eeff0array
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: aeff0_2array
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: beff0_2array
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: ceff0_2array
  real(wp),dimension(-14:14,0:maxdim),intent(inout):: eeff0_2array
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: aeff0_2auxarray
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: beff0_2auxarray
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: ceff0_2auxarray
  real(wp),dimension(-17:17,0:maxdim),intent(inout):: eeff0_2auxarray
  real(wp),dimension(0:n2,0:n1,0:n3):: xya_c, xyc_c
  real(wp),dimension(0:n3,0:n1,0:n2):: xza_c, xzc_c, yza_c, yzc_c
  real(wp),dimension(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xya_f
  real(wp),dimension(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xyb_f
  real(wp),dimension(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xyc_f
  real(wp),dimension(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xye_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xza_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xzb_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xzc_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xze_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yza_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yzb_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yzc_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yze_f
!  real(wp),dimension(35):: aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, ceff0, ceff1, ceff2, ceff3
!  real(wp),dimension(29):: eeff0, eeff1, eeff2, eeff3
!  real(wp),dimension(35):: aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2
!  real(wp),dimension(35):: ceff0_2, ceff1_2, ceff2_2, ceff3_2
!  real(wp),dimension(29):: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f

  ! Local variables
  !character(len=*),parameter :: subname='ConvolQuartic4'
  integer,parameter :: lowfil=-14,lupfil=14
  integer :: t,i1,i2,i3, icur,istart,iend,l
  real(wp) :: dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp) :: tt112, tt121, tt122, tt212, tt221, tt222, tt211
  real(kind=8) :: x0, y0, z0
  real(kind=8) :: tt10, tt11, tt12, tt13
  real(kind=8) :: tt20, tt21, tt22, tt23
  real(kind=8) :: tt30, tt31, tt32, tt33
  real(kind=8) :: tt40, tt41, tt42, tt43
  real(kind=8) :: tt50, tt51, tt52, tt53
  real(kind=8) :: tt60, tt61, tt62, tt63
  real(kind=8) :: tt70
  real(kind=8) :: tt0a0, tt0a1, tt0a2, tt0a3
  real(kind=8) :: tt0b0, tt0b1, tt0b2, tt0b3
  real(kind=8) :: tt0c0, tt0c1, tt0c2, tt0c3
  real(kind=8) :: tt0e0, tt0e1, tt0e2, tt0e3
  real(kind=8) :: tt1a0, tt1b0, tt1c0, tt1e0                     
  real(kind=8) :: tt2a0, tt2b0, tt2c0, tt2e0                     
  real(kind=8) :: tt3a0, tt3b0, tt3c0, tt3e0                     
  real(kind=8) :: tt4a0, tt4b0, tt4c0, tt4e0                     
  real(kind=8) :: tt5a0, tt5b0, tt5c0, tt5e0                     
  real(kind=8) :: tt6a0, tt6b0, tt6c0, tt6e0                     
  real(kind=8) :: tt7a0, tt7b0, tt7c0, tt7e0                     
  logical:: with_confpot
  real(kind=8) :: ddot,prefac1,hgrid2,hgrid3


  call f_routine(id='ConvolQuartic4')

  ! Flag indicating whether a confining quartic potential is used
  with_confpot=(potentialPrefac/=0.d0)


  prefac1=-.5d0/hgrid**2
  hgrid2=hgrid**2
  hgrid3=hgrid**3


  !initialize the arrays to zero.
  !This is important since the  bounding region can be concave
  call f_zero(y_c)
  call f_zero(y_f)

  !write(*,*) 'before: ddot',ddot((n1+1)*(n2+1)*(n3+1), y_c, 1, y_c, 1)

  !$omp parallel default(shared) &
  !$omp private(i1,i2,i3,x0,y0,z0,dyi,dyi0,dyi1,dyi2,dyi3,icur,istart,iend) &
  !$omp private(tt0a0,tt0a1,tt0a2,tt0a3,tt0b0,tt0b1,tt0b2,tt0b3) &
  !$omp private(tt0c0,tt0c1,tt0c2,tt0c3,tt0e0,tt0e1,tt0e2,tt0e3) &
  !$omp private(tt1a0,tt1b0,tt1c0,tt1e0,tt2a0,tt2b0,tt2c0,tt2e0,tt7a0,tt7b0) &                   
  !$omp private(tt3a0,tt3b0,tt3c0,tt3e0,tt4a0,tt4b0,tt4c0,tt4e0,tt7c0,tt7e0) & 
  !$omp private(tt5a0,tt5b0,tt5c0,tt5e0,tt6a0,tt6b0,tt6c0,tt6e0) &
  !$omp private(t112,t121,t122,t212,t221,t222,t211,tt112,tt121,tt122,tt212,tt221,tt222,tt211) &
  !$omp private(tt10,tt11,tt12,tt13,tt20,tt21,tt22,tt23,tt30,tt31,tt32,tt33,tt40,tt41,tt42,tt43) &
  !$omp private(tt50,tt51,tt52,tt53,tt60,tt61,tt62,tt63,tt70,t,l) !&
  !!$omp shared(hgrid,prefac1,hgrid2,hgrid3,offsetx,offsety,offsetz,rxyzConf) &
  !!$omp shared(with_kinetic,potentialPrefac,with_confpot,cprecr) &
  !!$omp shared(nfu1,nfu2,nfu3,n1,n2,n3,nfl1,nfl2,nfl3)&
  !!$omp shared(ibxy_c,ibxy_f,ibxz_c,ibyz_c,ibxz_f,ibyz_f,xx_c,xx_f,xx_f1,xy_c,xy_f,xz_f,xy_f2,xz_f4,xz_c)&
  !!$omp shared(y_c,y_f)&
  !!$omp shared(aeff0_2auxarray,beff0_2auxarray,ceff0_2auxarray,eeff0_2auxarray,aeff0array,beff0array,ceff0array,eeff0array)&
  !!$omp shared(aeff0_2array,beff0_2array,ceff0_2array,eeff0_2array)&
  !!$omp shared(xya_c,xyc_c,xza_c,xzc_c,yza_c,yzc_c)&
  !!$omp shared(xya_f,xyb_f,xyc_f,xye_f,xza_f,xzb_f,xzc_f,xze_f,yza_f,yzb_f,yzc_f,yze_f)
 
  !$omp do
  do i1=0,n1
     x0=hgrid*(i1+offsetx)-rxyzConf(1)
     if(.not. with_kinetic) then
        call position_dependent_filters(potentialPrefac, x0, aeff0array(lowfil,i1), 'a')
        call position_dependent_filters(potentialPrefac, x0, beff0array(lowfil,i1), 'b')
        call position_dependent_filters(potentialPrefac, x0, ceff0array(lowfil,i1), 'c')
        call position_dependent_filters(potentialPrefac, x0, eeff0array(lowfil,i1), 'e')
     else
        call position_dependent_filters(potentialPrefac,x0, aeff0array(lowfil,i1), 'aeff')
        call position_dependent_filters(potentialPrefac,x0, beff0array(lowfil,i1), 'beff')
        call position_dependent_filters(potentialPrefac,x0, ceff0array(lowfil,i1), 'ceff')
        call position_dependent_filters(potentialPrefac,x0, eeff0array(lowfil,i1), 'eeff')
     end if
     if(with_confpot) then
        call position_dependent_filters(1.d0,x0, aeff0_2auxarray(lowfil,i1), 'a2')
        call position_dependent_filters(1.d0,x0, beff0_2auxarray(lowfil,i1), 'b2')
        call position_dependent_filters(1.d0,x0, ceff0_2auxarray(lowfil,i1), 'c2')
        call position_dependent_filters(1.d0,x0, eeff0_2auxarray(lowfil,i1), 'e2')
     end if
  end do
  !$omp end do
  
  !$omp do schedule(static,1)
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 

              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + xx_c(t,i2,i3)*aeff0array(t-i1-0,i1+0)
                 dyi1=dyi1 + xx_c(t,i2,i3)*aeff0array(t-i1-1,i1+1)
                 dyi2=dyi2 + xx_c(t,i2,i3)*aeff0array(t-i1-2,i1+2)
                 dyi3=dyi3 + xx_c(t,i2,i3)*aeff0array(t-i1-3,i1+3)
              end do

              y_c(i1+0,i2,i3)=dyi0+cprecr*xx_c(i1+0,i2,i3)
              y_c(i1+1,i2,i3)=dyi1+cprecr*xx_c(i1+1,i2,i3)
              y_c(i1+2,i2,i3)=dyi2+cprecr*xx_c(i1+2,i2,i3)
              y_c(i1+3,i2,i3)=dyi3+cprecr*xx_c(i1+3,i2,i3)

              ! sss coefficients
              if(with_confpot) then

                 tt0a0=0.d0 ; tt0a1=0.d0 ; tt0a2=0.d0 ; tt0a3=0.d0
                 tt0b0=0.d0 ; tt0b1=0.d0 ; tt0b2=0.d0 ; tt0b3=0.d0
                 tt0c0=0.d0 ; tt0c1=0.d0 ; tt0c2=0.d0 ; tt0c3=0.d0
                 tt0e0=0.d0 ; tt0e1=0.d0 ; tt0e2=0.d0 ; tt0e3=0.d0

                 do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                    tt0a0=tt0a0 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-0,i1+0)
                    tt0a1=tt0a1 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-1,i1+1)
                    tt0a2=tt0a2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-2,i1+2)
                    tt0a3=tt0a3 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-3,i1+3)

                    tt0b0=tt0b0 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-0,i1+0)
                    tt0b1=tt0b1 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-1,i1+1)
                    tt0b2=tt0b2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-2,i1+2)
                    tt0b3=tt0b3 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-3,i1+3)

                    tt0c0=tt0c0 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-0,i1+0)
                    tt0c1=tt0c1 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-1,i1+1)
                    tt0c2=tt0c2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-2,i1+2)
                    tt0c3=tt0c3 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-3,i1+3)

                    tt0e0=tt0e0 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-0,i1+0)
                    tt0e1=tt0e1 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-1,i1+1)
                    tt0e2=tt0e2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-2,i1+2)
                    tt0e3=tt0e3 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-3,i1+3)
                 end do


                 xya_c(i2,i1+0,i3)=tt0a0
                 xya_c(i2,i1+1,i3)=tt0a1
                 xya_c(i2,i1+2,i3)=tt0a2
                 xya_c(i2,i1+3,i3)=tt0a3
                 xza_c(i3,i1+0,i2)=tt0a0
                 xza_c(i3,i1+1,i2)=tt0a1
                 xza_c(i3,i1+2,i2)=tt0a2
                 xza_c(i3,i1+3,i2)=tt0a3

                 xyc_c(i2,i1+0,i3)=tt0c0
                 xyc_c(i2,i1+1,i3)=tt0c1
                 xyc_c(i2,i1+2,i3)=tt0c2
                 xyc_c(i2,i1+3,i3)=tt0c3
                 xzc_c(i3,i1+0,i2)=tt0c0
                 xzc_c(i3,i1+1,i2)=tt0c1
                 xzc_c(i3,i1+2,i2)=tt0c2
                 xzc_c(i3,i1+3,i2)=tt0c3
              end if

           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp ; tt0a0=0.d0 ; tt0b0=0.d0 ; tt0c0=0.d0 ; tt0e0=0.d0

           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + xx_c(t,i2,i3)*aeff0array(t-i1,i1)
           end do

           y_c(i1,i2,i3)=dyi+cprecr*xx_c(i1,i2,i3)


           if(with_confpot) then
              ! sss coefficients
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                 tt0a0=tt0a0 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1,i1)
                 tt0b0=tt0b0 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1,i1)
                 tt0c0=tt0c0 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1,i1)
                 tt0e0=tt0e0 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1,i1)

                 xya_c(i2,i1,i3)=tt0a0
                 xza_c(i3,i1,i2)=tt0a0

                 xyc_c(i2,i1,i3)=tt0c0
                 xzc_c(i3,i1,i2)=tt0c0

              enddo
           end if

        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (iend-istart.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp

              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + xx_f1(t,i2,i3)*beff0array(t-i1-0,i1+0)
                 dyi1=dyi1 + xx_f1(t,i2,i3)*beff0array(t-i1-1,i1+1)
                 dyi2=dyi2 + xx_f1(t,i2,i3)*beff0array(t-i1-2,i1+2)
                 dyi3=dyi3 + xx_f1(t,i2,i3)*beff0array(t-i1-3,i1+3)
              enddo

              y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
              y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
              y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
              y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3

           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0.0_wp

           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + xx_f1(t,i2,i3)*beff0array(t-i1,i1)
           enddo

           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi

        enddo

        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 

              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + xx_c(t,i2,i3)*ceff0array(t-i1-0,i1+0)
                 dyi1=dyi1 + xx_c(t,i2,i3)*ceff0array(t-i1-1,i1+1)
                 dyi2=dyi2 + xx_c(t,i2,i3)*ceff0array(t-i1-2,i1+2)
                 dyi3=dyi3 + xx_c(t,i2,i3)*ceff0array(t-i1-3,i1+3)
              enddo

              y_f(1,i1+0,i2,i3)=dyi0
              y_f(1,i1+1,i2,i3)=dyi1
              y_f(1,i1+2,i2,i3)=dyi2
              y_f(1,i1+3,i2,i3)=dyi3

           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + xx_c(t,i2,i3)*ceff0array(t-i1,i1)
           enddo

           y_f(1,i1,i2,i3)=dyi

        enddo
     enddo
  enddo
  !$omp end do

  ! wavelet part
  !$omp do schedule(static,1) 
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 

           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + xx_f(4,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(5,i1+l,i2,i3)*beff0array(l,i1)
              t121=t121 + xx_f(2,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(3,i1+l,i2,i3)*beff0array(l,i1)
              t122=t122 + xx_f(6,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(7,i1+l,i2,i3)*beff0array(l,i1)
              t212=t212 + xx_f(4,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(5,i1+l,i2,i3)*eeff0array(l,i1)
              t221=t221 + xx_f(2,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(3,i1+l,i2,i3)*eeff0array(l,i1)
              t222=t222 + xx_f(6,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(7,i1+l,i2,i3)*eeff0array(l,i1)
              t211=t211 + xx_f(1,i1+l,i2,i3)*eeff0array(l,i1)
           end do


           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*xx_f(1,i1,i2,i3)
           y_f(2,i1,i2,i3)=t121+cprecr*xx_f(2,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221+cprecr*xx_f(3,i1,i2,i3)
           y_f(4,i1,i2,i3)=t112+cprecr*xx_f(4,i1,i2,i3)   
           y_f(5,i1,i2,i3)=t212+cprecr*xx_f(5,i1,i2,i3)   
           y_f(6,i1,i2,i3)=t122+cprecr*xx_f(6,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222+cprecr*xx_f(7,i1,i2,i3)

           if(with_confpot) then


              tt1a0=0.d0 ; tt1b0=0.d0 ; tt1c0=0.d0 ; tt1e0=0.d0
              tt2a0=0.d0 ; tt2b0=0.d0 ; tt2c0=0.d0 ; tt2e0=0.d0
              tt3a0=0.d0 ; tt3b0=0.d0 ; tt3c0=0.d0 ; tt3e0=0.d0
              tt4a0=0.d0 ; tt4b0=0.d0 ; tt4c0=0.d0 ; tt4e0=0.d0
              tt5a0=0.d0 ; tt5b0=0.d0 ; tt5c0=0.d0 ; tt5e0=0.d0
              tt6a0=0.d0 ; tt6b0=0.d0 ; tt6c0=0.d0 ; tt6e0=0.d0
              tt7a0=0.d0 ; tt7b0=0.d0 ; tt7c0=0.d0 ; tt7e0=0.d0

              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                 ! dss coefficients
                 tt1b0=tt1b0 + xx_f(1,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                 tt1e0=tt1e0 + xx_f(1,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                 ! sds coefficients
                 tt2a0=tt2a0 + xx_f(2,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                 tt2c0=tt2c0 + xx_f(2,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                 ! dds coefficients
                 tt3b0=tt3b0 + xx_f(3,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                 tt3e0=tt3e0 + xx_f(3,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                 ! ssd coefficients
                 tt4a0=tt4a0 + xx_f(4,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                 tt4c0=tt4c0 + xx_f(4,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                 ! dsd coefficients
                 tt5b0=tt5b0 + xx_f(5,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                 tt5e0=tt5e0 + xx_f(5,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                 ! sdd coefficients
                 tt6a0=tt6a0 + xx_f(6,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                 tt6c0=tt6c0 + xx_f(6,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                 ! ddd coefficients
                 tt7b0=tt7b0 + xx_f(7,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                 tt7e0=tt7e0 + xx_f(7,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
              enddo

              ! dss coefficients
              xyb_f(1,i2,i1,i3)=tt1b0
              xye_f(1,i2,i1,i3)=tt1e0
              xzb_f(1,i3,i1,i2)=tt1b0
              xze_f(1,i3,i1,i2)=tt1e0
              ! sds coefficients
              xya_f(1,i2,i1,i3)=tt2a0
              xyc_f(1,i2,i1,i3)=tt2c0
              xza_f(1,i3,i1,i2)=tt2a0
              xzc_f(1,i3,i1,i2)=tt2c0
              ! dds coefficients
              !xyb_f3(i2,i1,i3)=tt3b0
              xyb_f(2,i2,i1,i3)=tt3b0
              xye_f(2,i2,i1,i3)=tt3e0
              xzb_f(2,i3,i1,i2)=tt3b0
              xze_f(2,i3,i1,i2)=tt3e0
              ! ssd coefficients
              xya_f(2,i2,i1,i3)=tt4a0
              xyc_f(2,i2,i1,i3)=tt4c0
              xza_f(2,i3,i1,i2)=tt4a0
              xzc_f(2,i3,i1,i2)=tt4c0
              ! dsd coefficients
              xyb_f(3,i2,i1,i3)=tt5b0
              xye_f(3,i2,i1,i3)=tt5e0
              xzb_f(3,i3,i1,i2)=tt5b0
              xze_f(3,i3,i1,i2)=tt5e0
              ! sdd coefficients
              xya_f(3,i2,i1,i3)=tt6a0
              xyc_f(3,i2,i1,i3)=tt6c0
              xza_f(3,i3,i1,i2)=tt6a0
              xzc_f(3,i3,i1,i2)=tt6c0
              ! sdd coefficients
              xyb_f(4,i2,i1,i3)=tt7b0
              xye_f(4,i2,i1,i3)=tt7e0
              xzb_f(4,i3,i1,i2)=tt7b0
              xze_f(4,i3,i1,i2)=tt7e0
           end if
        enddo
     enddo
  enddo
  !$omp enddo

 
  !$omp do 
  do i2=0,n2

     y0=hgrid*(i2+offsety)-rxyzConf(2)
     if(.not. with_kinetic) then
        call position_dependent_filters(potentialPrefac,y0, aeff0array(lowfil,i2), 'a')
        call position_dependent_filters(potentialPrefac,y0, beff0array(lowfil,i2), 'b')
        call position_dependent_filters(potentialPrefac,y0, ceff0array(lowfil,i2), 'c')
        call position_dependent_filters(potentialPrefac,y0, eeff0array(lowfil,i2), 'e')
     else
        call position_dependent_filters(potentialPrefac,y0, aeff0array(lowfil,i2), 'aeff')
        call position_dependent_filters(potentialPrefac,y0, beff0array(lowfil,i2), 'beff')
        call position_dependent_filters(potentialPrefac,y0, ceff0array(lowfil,i2), 'ceff')
        call position_dependent_filters(potentialPrefac,y0, eeff0array(lowfil,i2), 'eeff')
     end if

     if(with_confpot) then
        call position_dependent_filters(potentialPrefac,y0, aeff0_2array(lowfil,i2), 'a2')
        call position_dependent_filters(potentialPrefac,y0, beff0_2array(lowfil,i2), 'b2')
        call position_dependent_filters(potentialPrefac,y0, ceff0_2array(lowfil,i2), 'c2')
        call position_dependent_filters(potentialPrefac,y0, eeff0_2array(lowfil,i2), 'e2')

        call position_dependent_filters(1.d0,y0, aeff0_2auxarray(lowfil,i2), 'a2')
        call position_dependent_filters(1.d0,y0, beff0_2auxarray(lowfil,i2), 'b2')
        call position_dependent_filters(1.d0,y0, ceff0_2auxarray(lowfil,i2), 'c2')
        call position_dependent_filters(1.d0,y0, eeff0_2auxarray(lowfil,i2), 'e2')
     end if

  end do
  !$omp end do

  !$omp do schedule(static,1) 
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp ;  dyi1=0.0_wp ;  dyi2=0.0_wp ;  dyi3=0.0_wp 

              if(with_confpot) then

                 tt0a0=0.d0 ; tt0a1=0.d0 ; tt0a2=0.d0 ; tt0a3=0.d0
                 tt0b0=0.d0 ; tt0b1=0.d0 ; tt0b2=0.d0 ; tt0b3=0.d0
                 tt0c0=0.d0 ; tt0c1=0.d0 ; tt0c2=0.d0 ; tt0c3=0.d0
                 tt0e0=0.d0 ; tt0e1=0.d0 ; tt0e2=0.d0 ; tt0e3=0.d0

                 do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                    dyi0=dyi0 + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0)
                    dyi1=dyi1 + xy_c(t,i1,i3)*aeff0array(t-i2-1,i2+1) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1)
                    dyi2=dyi2 + xy_c(t,i1,i3)*aeff0array(t-i2-2,i2+2) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2)
                    dyi3=dyi3 + xy_c(t,i1,i3)*aeff0array(t-i2-3,i2+3) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3)

                    tt0a0=tt0a0 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                    tt0a1=tt0a1 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                    tt0a2=tt0a2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                    tt0a3=tt0a3 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                    tt0b0=tt0b0 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                    tt0b1=tt0b1 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                    tt0b2=tt0b2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                    tt0b3=tt0b3 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                    tt0c0=tt0c0 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                    tt0c1=tt0c1 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                    tt0c2=tt0c2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                    tt0c3=tt0c3 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                    tt0e0=tt0e0 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                    tt0e1=tt0e1 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                    tt0e2=tt0e2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                    tt0e3=tt0e3 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)

                 end do

                 yza_c(i3,i1,i2+0)=tt0a0
                 yza_c(i3,i1,i2+1)=tt0a1
                 yza_c(i3,i1,i2+2)=tt0a2
                 yza_c(i3,i1,i2+3)=tt0a3

                 !!yzb_c(i3,i1,i2+0)=tt0b0
                 !!yzb_c(i3,i1,i2+1)=tt0b1
                 !!yzb_c(i3,i1,i2+2)=tt0b2
                 !!yzb_c(i3,i1,i2+3)=tt0b3

                 yzc_c(i3,i1,i2+0)=tt0c0
                 yzc_c(i3,i1,i2+1)=tt0c1
                 yzc_c(i3,i1,i2+2)=tt0c2
                 yzc_c(i3,i1,i2+3)=tt0c3

                 !!yze_c(i3,i1,i2+0)=tt0e0
                 !!yze_c(i3,i1,i2+1)=tt0e1
                 !!yze_c(i3,i1,i2+2)=tt0e2
                 !!yze_c(i3,i1,i2+3)=tt0e3
              else
                 do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                    dyi0=dyi0 + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0)
                    dyi1=dyi1 + xy_c(t,i1,i3)*aeff0array(t-i2-1,i2+1)
                    dyi2=dyi2 + xy_c(t,i1,i3)*aeff0array(t-i2-2,i2+2)
                    dyi3=dyi3 + xy_c(t,i1,i3)*aeff0array(t-i2-3,i2+3)
                 end do
              end if

              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3

           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0.0_wp
           if(with_confpot) then
              tt0a0=0.d0 ; tt0b0=0.d0 ; tt0c0=0.d0 ; tt0e0=0.d0
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                 dyi=dyi + xy_c(t,i1,i3)*aeff0array(t-i2,i2) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2,i2)

                 tt0a0=tt0a0 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2)

                 tt0b0=tt0b0 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2)

                 tt0c0=tt0c0 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2)

                 tt0e0=tt0e0 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2)
              end do
              yza_c(i3,i1,i2)=tt0a0

              !!yzb_c(i3,i1,i2)=tt0b0

              yzc_c(i3,i1,i2)=tt0c0

              !!yze_c(i3,i1,i2)=tt0e0
           else
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                 dyi=dyi + xy_c(t,i1,i3)*aeff0array(t-i2,i2)
              end do
           end if
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi

        enddo
        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (iend-istart.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp
              if(with_confpot) then
                 do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                    dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                         2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2+0)
                    dyi1=dyi1 + xy_f2(t,i1,i3)*beff0array(t-i2-1,i2+1) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-1,i2+1) + &
                         2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-1,i2+1)
                    dyi2=dyi2 + xy_f2(t,i1,i3)*beff0array(t-i2-2,i2+2) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-2,i2+2) + &
                         2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-2,i2+2)
                    dyi3=dyi3 + xy_f2(t,i1,i3)*beff0array(t-i2-3,i2+3) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-3,i2+3) + &
                         2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-3,i2+3)
                 enddo
              else
                 do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                    dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0)
                    dyi1=dyi1 + xy_f2(t,i1,i3)*beff0array(t-i2-1,i2+1)
                    dyi2=dyi2 + xy_f2(t,i1,i3)*beff0array(t-i2-2,i2+2)
                    dyi3=dyi3 + xy_f2(t,i1,i3)*beff0array(t-i2-3,i2+3)
                 enddo
              end if
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
           enddo
           istart=i2
        endif

        do i2=istart,iend
           dyi0=0.0_wp
           if(with_confpot) then
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2) + &
                      2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2)
              enddo
           else
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2)
              enddo
           end if
           y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
        enddo

        if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
              tt10=0.0_wp ; tt11=0.0_wp ; tt12=0.0_wp ; tt13=0.0_wp 
              tt20=0.0_wp ; tt21=0.0_wp ; tt22=0.0_wp ; tt23=0.0_wp 
              tt30=0.0_wp ; tt31=0.0_wp ; tt32=0.0_wp ; tt33=0.0_wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2+0)
                 dyi1=dyi1 + xy_c(t,i1,i3)*ceff0array(t-i2-1,i2+1)
                 dyi2=dyi2 + xy_c(t,i1,i3)*ceff0array(t-i2-2,i2+2)
                 dyi3=dyi3 + xy_c(t,i1,i3)*ceff0array(t-i2-3,i2+3)
              end do
              y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
              y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
              y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
              y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3

              if(with_confpot) then
                 do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                    tt10=tt10 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0)
                    tt11=tt11 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1)
                    tt12=tt12 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2)
                    tt13=tt13 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3)

                    tt20=tt20 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0)
                    tt21=tt21 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-1,i2+1)
                    tt22=tt22 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-2,i2+2)
                    tt23=tt23 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-3,i2+3)

                    tt30=tt30 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0)
                    tt31=tt31 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-1,i2+1)
                    tt32=tt32 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-2,i2+2)
                    tt33=tt33 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-3,i2+3)
                 enddo

                 y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
                 y_f(1,i1,i2+1,i3)=y_f(1,i1,i2+1,i3)+tt11
                 y_f(1,i1,i2+2,i3)=y_f(1,i1,i2+2,i3)+tt12
                 y_f(1,i1,i2+3,i3)=y_f(1,i1,i2+3,i3)+tt13

                 y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
                 y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+tt21
                 y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+tt22
                 y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+tt23

                 y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
                 y_f(3,i1,i2+1,i3)=y_f(3,i1,i2+1,i3)+tt31
                 y_f(3,i1,i2+2,i3)=y_f(3,i1,i2+2,i3)+tt32
                 y_f(3,i1,i2+3,i3)=y_f(3,i1,i2+3,i3)+tt33
              end if
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi0=0.0_wp ; tt10=0.0_wp ; tt20=0.0_wp ; tt30=0.0_wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2)
           end do
           y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0

           if(with_confpot) then
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                 tt10=tt10 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2)

                 tt20=tt20 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2)

                 tt30=tt30 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2)
              enddo
              y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
              y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
              y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
           end if
        enddo
     enddo
  enddo
  !$omp end do 



  ! wavelet part

  !$omp do schedule(static,1) 
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           ! Get the effective filters for the y dimension
           tt10 = 0.d0 ; tt20 = 0.d0 ; tt30 = 0.d0 ; tt40 = 0.d0 ; tt50 = 0.d0 ; tt60 = 0.d0 ; tt70 = 0.d0

           if(with_confpot) then

              tt1a0=0.d0 ; tt1b0=0.d0 ; tt1c0=0.d0 ; tt1e0=0.d0
              tt2a0=0.d0 ; tt2b0=0.d0 ; tt2c0=0.d0 ; tt2e0=0.d0
              tt3a0=0.d0 ; tt3b0=0.d0 ; tt3c0=0.d0 ; tt3e0=0.d0
              tt4a0=0.d0 ; tt4b0=0.d0 ; tt4c0=0.d0 ; tt4e0=0.d0
              tt5a0=0.d0 ; tt5b0=0.d0 ; tt5c0=0.d0 ; tt5e0=0.d0
              tt6a0=0.d0 ; tt6b0=0.d0 ; tt6c0=0.d0 ; tt6e0=0.d0
              tt7a0=0.d0 ; tt7b0=0.d0 ; tt7c0=0.d0 ; tt7e0=0.d0

              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 tt10 = tt10 + xy_f(1,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(3,i2+l,i1,i3)*beff0array(l,i2) + &
                      2.d0*xye_f(1,i2+l,i1,i3)*aeff0_2array(l,i2) + &
                      2.d0*(xyc_f(1,i2+l,i1,i3)+xye_f(2,i2+l,i1,i3))*beff0_2array(l,i2)

                 tt20 = tt20 + xy_f(2,i2+l,i1,i3)*eeff0array(l,i2) +                                      &
                      2.d0*xyb_f(1,i2+l,i1,i3)*ceff0_2array(l,i2) + &
                      2.d0*(xya_f(1,i2+l,i1,i3)+xyb_f(2,i2+l,i1,i3))*eeff0_2array(l,i2)

                 tt30 = tt30 + xy_f(1,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(3,i2+l,i1,i3)*eeff0array(l,i2) + &
                      2.d0*xye_f(1,i2+l,i1,i3)*ceff0_2array(l,i2) + &
                      2.d0*(xyc_f(1,i2+l,i1,i3)+xye_f(2,i2+l,i1,i3))*eeff0_2array(l,i2)

                 tt40 = tt40 + xy_f(4,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(6,i2+l,i1,i3)*beff0array(l,i2) + &
                      2.d0*(xya_f(2,i2+l,i1,i3)+xyb_f(3,i2+l,i1,i3))*aeff0_2array(l,i2) + &
                      2.d0*(xya_f(3,i2+l,i1,i3)+xyb_f(4,i2+l,i1,i3))*beff0_2array(l,i2)

                 tt50 = tt50 + xy_f(5,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(7,i2+l,i1,i3)*beff0array(l,i2) + &
                      2.d0*(xyc_f(2,i2+l,i1,i3)+xye_f(3,i2+l,i1,i3))*aeff0_2array(l,i2) + &
                      2.d0*(xyc_f(3,i2+l,i1,i3)+xye_f(4,i2+l,i1,i3))*beff0_2array(l,i2)

                 tt60 = tt60 + xy_f(4,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(6,i2+l,i1,i3)*eeff0array(l,i2) + &
                      2.d0*(xya_f(2,i2+l,i1,i3)+xyb_f(3,i2+l,i1,i3))*ceff0_2array(l,i2) + &
                      2.d0*(xya_f(3,i2+l,i1,i3)+xyb_f(4,i2+l,i1,i3))*eeff0_2array(l,i2)

                 tt70 = tt70 + xy_f(5,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(7,i2+l,i1,i3)*eeff0array(l,i2) + &
                      2.d0*(xyc_f(2,i2+l,i1,i3)+xye_f(3,i2+l,i1,i3))*ceff0_2array(l,i2) + &
                      2.d0*(xyc_f(3,i2+l,i1,i3)+xye_f(4,i2+l,i1,i3))*eeff0_2array(l,i2)

                 ! dss coefficients
                 tt1a0=tt1a0 + xy_f(1,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                 tt1c0=tt1c0 + xy_f(1,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                 ! sds coefficients
                 tt2b0=tt2b0 + xy_f(2,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                 tt2e0=tt2e0 + xy_f(2,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                 ! dds coefficients
                 tt3b0=tt3b0 + xy_f(3,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                 tt3e0=tt3e0 + xy_f(3,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                 ! ssd coefficients
                 tt4a0=tt4a0 + xy_f(4,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                 tt4c0=tt4c0 + xy_f(4,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                 ! dsd coefficients
                 tt5a0=tt5a0 + xy_f(5,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                 tt5c0=tt5c0 + xy_f(5,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                 ! sdd coefficients
                 tt6b0=tt6b0 + xy_f(6,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                 tt6e0=tt6e0 + xy_f(6,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                 ! ddd coefficients
                 tt7b0=tt7b0 + xy_f(7,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                 tt7e0=tt7e0 + xy_f(7,i2+l,i1,i3)*eeff0_2auxarray(l,i2)

              end do
              ! dss coefficients
              yza_f(1,i3,i1,i2)=tt1a0
              yzc_f(1,i3,i1,i2)=tt1c0
              ! sds coefficients
              yzb_f(1,i3,i1,i2)=tt2b0
              yze_f(1,i3,i1,i2)=tt2e0
              ! dds coefficients
              yzb_f(2,i3,i1,i2)=tt3b0
              yze_f(2,i3,i1,i2)=tt3e0
              ! ssd coefficients
              yza_f(2,i3,i1,i2)=tt4a0
              yzc_f(2,i3,i1,i2)=tt4c0
              ! dsd coefficients
              yza_f(3,i3,i1,i2)=tt5a0
              yzc_f(3,i3,i1,i2)=tt5c0
              ! sdd coefficients
              yzb_f(3,i3,i1,i2)=tt6b0
              yze_f(3,i3,i1,i2)=tt6e0
              ! sdd coefficients
              yzb_f(4,i3,i1,i2)=tt7b0
              yze_f(4,i3,i1,i2)=tt7e0

           else
              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 tt10 = tt10 + xy_f(1,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(3,i2+l,i1,i3)*beff0array(l,i2)

                 tt20 = tt20 + xy_f(2,i2+l,i1,i3)*eeff0array(l,i2)

                 tt30 = tt30 + xy_f(1,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(3,i2+l,i1,i3)*eeff0array(l,i2)

                 tt40 = tt40 + xy_f(4,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(6,i2+l,i1,i3)*beff0array(l,i2)

                 tt50 = tt50 + xy_f(5,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(7,i2+l,i1,i3)*beff0array(l,i2)

                 tt60 = tt60 + xy_f(4,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(6,i2+l,i1,i3)*eeff0array(l,i2)

                 tt70 = tt70 + xy_f(5,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(7,i2+l,i1,i3)*eeff0array(l,i2)
              end do
           end if
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70

        enddo
     enddo
  enddo
  !$omp enddo 


  !$omp do 
  do i3=0,n3
     z0=hgrid*(i3+offsetz)-rxyzConf(3)
     if(.not. with_kinetic) then
        call position_dependent_filters(potentialPrefac,z0, aeff0array(lowfil,i3), 'a')
        call position_dependent_filters(potentialPrefac,z0, beff0array(lowfil,i3), 'b')
        call position_dependent_filters(potentialPrefac,z0, ceff0array(lowfil,i3), 'c')
        call position_dependent_filters(potentialPrefac,z0, eeff0array(lowfil,i3), 'e')
     else
        call position_dependent_filters(potentialPrefac, z0, aeff0array(lowfil,i3), 'aeff')
        call position_dependent_filters(potentialPrefac, z0, beff0array(lowfil,i3), 'beff')
        call position_dependent_filters(potentialPrefac, z0, ceff0array(lowfil,i3), 'ceff')
        call position_dependent_filters(potentialPrefac, z0, eeff0array(lowfil,i3), 'eeff')
     end if
     if(with_confpot) then
        call position_dependent_filters(potentialPrefac,z0, aeff0_2array(lowfil,i3), 'a2')
        call position_dependent_filters(potentialPrefac,z0, beff0_2array(lowfil,i3), 'b2')
        call position_dependent_filters(potentialPrefac,z0, ceff0_2array(lowfil,i3), 'c2')
        call position_dependent_filters(potentialPrefac,z0, eeff0_2array(lowfil,i3), 'e2')
     end if
  end do
  !$omp end do

  !!!in this last section a bug of ifort with fpe0 is found with OpenMP
  
  ! + (1/2) d^2/dz^2
  !$omp do schedule(static,1) 
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
              if(with_confpot) then
                 do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                    dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0) + &
                         2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0)
                    dyi1=dyi1 + xz_c(t,i1,i2)*aeff0array(t-i3-1,i3+1) + &
                         2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-1,i3+1)
                    dyi2=dyi2 + xz_c(t,i1,i2)*aeff0array(t-i3-2,i3+2) + &
                         2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-2,i3+2)
                    dyi3=dyi3 + xz_c(t,i1,i2)*aeff0array(t-i3-3,i3+3) + &
                         2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-3,i3+3)
                 enddo
              else
                 do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                    dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0)
                    dyi1=dyi1 + xz_c(t,i1,i2)*aeff0array(t-i3-1,i3+1)
                    dyi2=dyi2 + xz_c(t,i1,i2)*aeff0array(t-i3-2,i3+2)
                    dyi3=dyi3 + xz_c(t,i1,i2)*aeff0array(t-i3-3,i3+3)
                 enddo
              end if
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0.0_wp
           if(with_confpot) then
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                 dyi=dyi + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3) + &
                      2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0)
              enddo
           else
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                 dyi=dyi + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3)
              enddo
           end if
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        !the openMP problem is between here ....
        if (iend-istart.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp
              if(with_confpot) then
                 do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                    dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0) + &
                         2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-0,i3+0) + &
                         2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                         *beff0_2array(t-i3-0,i3+0)
                    dyi1 = dyi1 + xz_f4(t,i1,i2)*beff0array(t-i3-1,i3+1) + &
                         2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-1,i3+1) + &
                         2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                         *beff0_2array(t-i3-1,i3+1)
                    dyi2 = dyi2 + xz_f4(t,i1,i2)*beff0array(t-i3-2,i3+2) + &
                         2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-2,i3+2) + &
                         2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                         *beff0_2array(t-i3-2,i3+2)
                    dyi3 = dyi3 + xz_f4(t,i1,i2)*beff0array(t-i3-3,i3+3) + &
                         2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-3,i3+3) + &
                         2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                         *beff0_2array(t-i3-3,i3+3)
                 enddo
              else
                 do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                    dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0)
                    dyi1 = dyi1 + xz_f4(t,i1,i2)*beff0array(t-i3-1,i3+1)
                    dyi2 = dyi2 + xz_f4(t,i1,i2)*beff0array(t-i3-2,i3+2)
                    dyi3 = dyi3 + xz_f4(t,i1,i2)*beff0array(t-i3-3,i3+3)
                 enddo
              end if
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i3
        endif

        do i3=istart,iend
           dyi=0.0_wp
           if(with_confpot) then
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
                 dyi=dyi + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3) + &
                      2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-0,i3) + &
                      2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                      *beff0_2array(t-i3-0,i3)
              enddo
           else
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
                 dyi=dyi + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3)
              enddo
           end if
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo
        !...and here


        if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 

              tt10 = 0.d0 ; tt11 = 0.d0 ; tt12 = 0.d0 ; tt13 = 0.d0
              tt40 = 0.d0 ; tt41 = 0.d0 ; tt42 = 0.d0 ; tt43 = 0.d0
              tt50 = 0.d0 ; tt51 = 0.d0 ; tt52 = 0.d0 ; tt53 = 0.d0
              tt20 = 0.d0 ; tt21 = 0.d0 ; tt22 = 0.d0 ; tt23 = 0.d0
              tt60 = 0.d0 ; tt61 = 0.d0 ; tt62 = 0.d0 ; tt63 = 0.d0
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3+0)
                 dyi1=dyi1 + xz_c(t,i1,i2)*ceff0array(t-i3-1,i3+1)
                 dyi2=dyi2 + xz_c(t,i1,i2)*ceff0array(t-i3-2,i3+2)
                 dyi3=dyi3 + xz_c(t,i1,i2)*ceff0array(t-i3-3,i3+3)
              end do

              if(with_confpot) then


                 do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                    tt10 = tt10 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                    tt11 = tt11 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                    tt12 = tt12 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                    tt13 = tt13 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                    tt40 = tt40 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt41 = tt41 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt42 = tt42 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt43 = tt43 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                    tt50 = tt50 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt51 = tt51 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt52 = tt52 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt53 = tt53 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                    tt20 = tt20 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                    tt21 = tt21 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                    tt22 = tt22 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                    tt23 = tt23 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                    tt40 = tt40 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt41 = tt41 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt42 = tt42 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt43 = tt43 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                    tt60 = tt60 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt61 = tt61 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt62 = tt62 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt63 = tt63 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)
                 enddo
                 y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
                 y_f(1,i1,i2,i3+1) = y_f(1,i1,i2,i3+1) + tt11
                 y_f(1,i1,i2,i3+2) = y_f(1,i1,i2,i3+2) + tt12
                 y_f(1,i1,i2,i3+3) = y_f(1,i1,i2,i3+3) + tt13

                 y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
                 y_f(5,i1,i2,i3+1) = y_f(5,i1,i2,i3+1) + tt51
                 y_f(5,i1,i2,i3+2) = y_f(5,i1,i2,i3+2) + tt52
                 y_f(5,i1,i2,i3+3) = y_f(5,i1,i2,i3+3) + tt53

                 y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
                 y_f(2,i1,i2,i3+1) = y_f(2,i1,i2,i3+1) + tt21
                 y_f(2,i1,i2,i3+2) = y_f(2,i1,i2,i3+2) + tt22
                 y_f(2,i1,i2,i3+3) = y_f(2,i1,i2,i3+3) + tt23

                 y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
                 y_f(6,i1,i2,i3+1) = y_f(6,i1,i2,i3+1) + tt61
                 y_f(6,i1,i2,i3+2) = y_f(6,i1,i2,i3+2) + tt62
                 y_f(6,i1,i2,i3+3) = y_f(6,i1,i2,i3+3) + tt63
              end if
              y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
              y_f(4,i1,i2,i3+1) = y_f(4,i1,i2,i3+1) + dyi1 + tt41
              y_f(4,i1,i2,i3+2) = y_f(4,i1,i2,i3+2) + dyi2 + tt42
              y_f(4,i1,i2,i3+3) = y_f(4,i1,i2,i3+3) + dyi3 + tt43
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi0=0.0_wp ;  tt10 = 0.d0 ; tt40 = 0.d0 ; tt50 = 0.d0 ; tt20 = 0.d0 ; tt60 = 0.d0
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3)
           end do

           if(with_confpot) then

              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                 tt10 = tt10 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3)
                 tt40 = tt40 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)
                 tt50 = tt50 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)
                 tt20 = tt20 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3)
                 tt40 = tt40 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)
                 tt60 = tt60 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)
              enddo
           end if
           y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
           y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
           y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
           y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
           y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
        enddo
     enddo
  enddo
  !$omp enddo


  ! wavelet part

  !$omp do schedule(static,1) 
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           tt10 = 0.d0 ; tt20 = 0.d0 ; tt30 = 0.d0 ; tt40 = 0.d0 ; tt50 = 0.d0 ; tt60 = 0.d0 ; tt70 = 0.d0

           if(with_confpot) then
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 tt10 = tt10 + xz_f(1,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(5,i3+l,i1,i2)*beff0array(l,i3) + &
                      2.d0*                    xze_f(1,i3+l,i1,i2) *aeff0_2array(l,i3) + &
                      2.d0*(xzc_f(2,i3+l,i1,i2)+xze_f(3,i3+l,i1,i2))*beff0_2array(l,i3) + &
                      2.d0*(yza_f(1,i3+l,i1,i2)+yzb_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                      2.d0*(yza_f(3,i3+l,i1,i2)+yzb_f(4,i3+l,i1,i2))*beff0_2array(l,i3)

                 tt20 = tt20 + xz_f(2,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(6,i3+l,i1,i2)*beff0array(l,i3) + &
                      2.d0*(xza_f(1,i3+l,i1,i2)+xzb_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                      2.d0*(xza_f(3,i3+l,i1,i2)+xzb_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                      2.d0*                    yze_f(1,i3+l,i1,i2) *aeff0_2array(l,i3) + &
                      2.d0*(yzc_f(2,i3+l,i1,i2)+yze_f(3,i3+l,i1,i2))*beff0_2array(l,i3)

                 tt30 = tt30 + xz_f(3,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(7,i3+l,i1,i2)*beff0array(l,i3) + &
                      2.d0*(xzc_f(1,i3+l,i1,i2)+xze_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                      2.d0*(xzc_f(3,i3+l,i1,i2)+xze_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                      2.d0*(yzc_f(1,i3+l,i1,i2)+yze_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                      2.d0*(yzc_f(3,i3+l,i1,i2)+yze_f(4,i3+l,i1,i2))*beff0_2array(l,i3)

                 tt40 = tt40 + xz_f(4,i3+l,i1,i2)*eeff0array(l,i3)                                      + &
                      2.d0*                    xzb_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                      2.d0*(xza_f(2,i3+l,i1,i2)+xzb_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                      2.d0*                    yzb_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                      2.d0*(yza_f(2,i3+l,i1,i2)+yzb_f(3,i3+l,i1,i2))*eeff0_2array(l,i3)

                 tt50 = tt50 + xz_f(1,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(5,i3+l,i1,i2)*eeff0array(l,i3) + &
                      2.d0*                    xze_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                      2.d0*(xzc_f(2,i3+l,i1,i2)+xze_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                      2.d0*(yza_f(1,i3+l,i1,i2)+yzb_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                      2.d0*(yza_f(3,i3+l,i1,i2)+yzb_f(4,i3+l,i1,i2))*eeff0_2array(l,i3)

                 tt60 = tt60 + xz_f(2,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(6,i3+l,i1,i2)*eeff0array(l,i3) + &
                      2.d0*(xza_f(1,i3+l,i1,i2)+xzb_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                      2.d0*(xza_f(3,i3+l,i1,i2)+xzb_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                      2.d0*                    yze_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                      2.d0*(yzc_f(2,i3+l,i1,i2)+yze_f(3,i3+l,i1,i2))*eeff0_2array(l,i3)

                 tt70 = tt70 + xz_f(3,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(7,i3+l,i1,i2)*eeff0array(l,i3) + &
                      2.d0*(xzc_f(1,i3+l,i1,i2)+xze_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                      2.d0*(xzc_f(3,i3+l,i1,i2)+xze_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                      2.d0*(yzc_f(1,i3+l,i1,i2)+yze_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                      2.d0*(yzc_f(3,i3+l,i1,i2)+yze_f(4,i3+l,i1,i2))*eeff0_2array(l,i3)
              enddo
           else
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 tt10 = tt10 + xz_f(1,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(5,i3+l,i1,i2)*beff0array(l,i3)
                 tt20 = tt20 + xz_f(2,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(6,i3+l,i1,i2)*beff0array(l,i3)
                 tt30 = tt30 + xz_f(3,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(7,i3+l,i1,i2)*beff0array(l,i3)
                 tt40 = tt40 + xz_f(4,i3+l,i1,i2)*eeff0array(l,i3) 
                 tt50 = tt50 + xz_f(1,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(5,i3+l,i1,i2)*eeff0array(l,i3)
                 tt60 = tt60 + xz_f(2,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(6,i3+l,i1,i2)*eeff0array(l,i3)
                 tt70 = tt70 + xz_f(3,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(7,i3+l,i1,i2)*eeff0array(l,i3)
              enddo
           end if
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70
        enddo
     enddo
  enddo
  !$omp enddo

  !$omp end parallel

!  write(*,*) 'after: ddot',ddot((n1+1)*(n2+1)*(n3+1), y_c, 1, y_c, 1),&
!       ddot(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),y_f,1,y_f,1)

  call f_release_routine()

  contains
    !> identify and evaluate filters associated to the position
    pure subroutine position_dependent_filters(parabPrefac, x0, eff, filterCode)
!      use filterModule, only:lb,ub,a,b,c,e,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,e1,e2,e3,e4
      implicit none
      ! Calling arguments
      real(kind=8),intent(in) :: parabPrefac !<prefactor of the potential
      real(kind=8),intent(in) :: x0          !< Center of the parabolic potential (x-x0)^2
      real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
      character(len=*), intent(in) :: filterCode
      ! Local variables
      integer :: i
      real(kind=8) :: fac, fac2, x02, x03

      fac=parabPrefac
      fac2=parabPrefac*hgrid
      x02=x0**2
      x03=x0**3

      ! Determine which filter we have to calculate
      select case(filterCode)
      case('aeff')
         do i=lb,ub
            eff(i)=prefac1*a(i) + fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
         end do
         eff(0)=eff(0)+fac*x0**4
      case('beff')
         do i=lb,ub
            eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
         end do
      case('ceff')
         do i=lb,ub
            eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
         end do
      case('eeff')
         do i=lb,ub
            eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
         end do
         eff(0)=eff(0)+fac*x0**4
      case('a')
         do i=lb,ub
            eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
         end do
         eff(0)=eff(0)+fac*x0**4
      case('b')
         do i=lb,ub
            eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
         end do
      case('c')
         do i=lb,ub
            eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
         end do
      case('e')
         do i=lb,ub
            eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
         end do
         eff(0)=eff(0)+fac*x0**4
      case('a2')
         do i=lb,ub
            eff(i) = fac2*( hgrid*a2(i) + 2.d0*x0*a1(i) )
         end do
         eff(0)=eff(0)+fac*x0**2
      case('b2')
         do i=lb,ub
            eff(i) = fac2*( hgrid*b2(i) + 2.d0*x0*b1(i) )
         end do
      case('c2')
         do i=lb,ub
            eff(i) = fac2*( hgrid*c2(i) + 2.d0*x0*c1(i) )
         end do
      case('e2')
         do i=lb,ub
            eff(i) = fac2*( hgrid*e2(i) + 2.d0*x0*e1(i) )
         end do
         eff(0)=eff(0)+fac*x0**2
      end select

    end subroutine position_dependent_filters



END SUBROUTINE ConvolQuartic4

