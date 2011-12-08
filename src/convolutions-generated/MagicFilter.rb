
FILTER = [ "8.4334247333529341094733325815816e-7",
       "-0.1290557201342060969516786758559028e-4",
       "0.8762984476210559564689161894116397e-4",
       "-0.30158038132690463167163703826169879e-3",
       "0.174723713672993903449447812749852942e-2",
       "-0.942047030201080385922711540948195075e-2",
       "0.2373821463724942397566389712597274535e-1",
       "0.612625895831207982195380597e-1",
       "0.9940415697834003993178616713",
       "-0.604895289196983516002834636e-1",
       "-0.2103025160930381434955489412839065067e-1",
       "0.1337263414854794752733423467013220997e-1",
       "-0.344128144493493857280881509686821861e-2",
       "0.49443227688689919192282259476750972e-3",
       "-0.5185986881173432922848639136911487e-4",
       "2.72734492911979659657715313017228e-6" ]



def MagicFilter(filt, center, unroll, invert, free=false )
  function_name = "magicfilter"
  if free then
    function_name += "_free"
  else 
    function_name += "_per"
  end
  if invert then
    function_name += "_inv"
  end
  dim_in_min = "0"
  dim_in_max = "n-1"
  dim_out_min = "0"
  dim_out_max = "n-1"
  if free then
    if invert then
      dim_in_min = "lowfil"
      dim_in_max = "n-1+upfil"
    else
      dim_out_min = "-upfil"
      dim_out_max = "n-1-lowfil"
    end
  end
 
  puts "\
subroutine #{function_name}(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat"
  if invert then
  puts "\
  integer, parameter :: lowfil=-#{filt.length-center-1},upfil=#{center}"
  else
  puts "\
  integer, parameter :: lowfil=-#{center},upfil=#{filt.length-center-1}"
  end
  puts "\ 
  real(wp), dimension(#{dim_in_min}:#{dim_in_max},ndat), intent(in) :: x
  real(wp), dimension(ndat,#{dim_out_min}:#{dim_out_max}), intent(out) :: y
"
 
  print"\
  integer :: i,j,k,l
  real(wp) :: tt1"
  2.upto(unroll) { |i|
    print ",tt#{i}"
  }
  puts

  puts"\
  real(wp) fil(lowfil:upfil)
  DATA fil / &"
  if invert then sub = filt.reverse else sub =filt end
  sub[0..-2].each { |v|
    puts v+"_wp,&"
  }
  puts sub.last+"_wp /"
  
  puts "
  do j=1,ndat#{unroll>1?"#{-unroll+1},#{unroll}":",1"}
    do i=#{dim_out_min},#{dim_out_max}
      tt1=0.e0_wp"
  2.upto(unroll){ |i| puts"\
      tt#{i}=0.e0_wp"
  }
  if free and not invert then
    puts"\
      do l=max(-i,lowfil),min(lupfil,n1-i)"
  else
    puts "\
      do l=lowfil,upfil"
  end
  if not free then
    puts "\
        k=modulo(i+l,n)"
  else
    puts "\
        k=i+l"
  end
  puts "\
        tt1=tt1+x( k,j  )*fil(l)"
  2.upto(unroll){ |i| puts "\
        tt#{i}=tt#{i}+x( k,j+#{i-1})*fil(l)"
  }
  puts "\
      enddo
      y(j  ,i)=tt1"
  2.upto(unroll){ |i| puts "\
      y(j+#{i-1},i)=tt#{i}"
  }
  puts "\
    enddo
  enddo"

  if unroll>2 then puts "\
  do j=ndat-modulo(ndat,#{unroll})+1,ndat
    do i=#{dim_out_min},#{dim_out_max}
      tt1=0.e0_wp"
  if free and not invert then
    puts "\
      do l=max(-i,lowfil),min(lupfil,n1-i)"
  else
    puts "\
      do l=lowfil,upfil"
  end
  if not free then
    puts "\
        k=modulo(i+l,n)"
  else
    puts "\
        k=i+l"
  end
  puts "\
        tt1=tt1+x( k,j  )*fil(l)
      enddo
      y(j  ,i)=tt1
    enddo
  end"
  end
  puts "\
END SUBROUTINE #{function_name}"
end

MagicFilter(FILTER,8,0,false)
MagicFilter(FILTER,8,5,true)
MagicFilter(FILTER,8,3,false,true)
MagicFilter(FILTER,8,4,true,true)
  

