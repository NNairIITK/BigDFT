require 'yaml'


MAX_N1=256
MAX_N2=256
MAX_N3=256
N1_STEP=64
N2_STEP=64
N3_STEP=64
N1_URANGE=16 
N2_URANGE=16 
N3_URANGE=0
N1_USTEP=8
N2_USTEP=1
N3_USTEP=1

sizes = []

for n1 in (N1_STEP..MAX_N1).step(N1_STEP)
  for un1 in ((n1-N1_URANGE)..(n1+N1_URANGE)).step(N1_USTEP)
    for n2 in (N2_STEP..MAX_N2).step(N2_STEP)
      for un2 in ((n2-N2_URANGE)..(n2+N2_URANGE)).step(N2_USTEP)
        for n3 in (n2..MAX_N3).step(N3_STEP)
          for un3 in ((n3-N3_URANGE)..(n3+N3_URANGE)).step(N3_USTEP)
            sizes.push( [un1, un2, un3 , un2*un3, un1*un2*un3] )
#            puts [un1, un2, un3].inspect
          end
        end
      end
    end
  end
end

results = []

f = File::new("opencl_profile_0.log","r")
4.times {f.gets}
while line = f.gets
  res = line.scan(/\s*(.*?)=\[\s*(.*?)\s*?\]/)
  res.each { |elem|
    sub = elem[1].scan(/([\d\.]+)/)
    if sub.length > 1 then
      elem[1] = sub.flatten!
    end
  }
  h = {}
  res.each{ |r|
    h[r[0]]=r[1]
  }
#  puts res.inspect
  results.push(h) if not (h["method"].match(/memcpy/) or h["method"].match(/init/))
end
times = Hash::new { |hash,key| hash[key] = [] }
f.close

results.each { |res|
  times[res["method"]].push([ res["gputime"].to_f ])
}

puts sizes.length

times.each{ |key,value|
  puts key
  puts value.length
#  puts value.inspect
  for i in 0...sizes.length
    value[i].concat(sizes[i])
  end
}

times.each{ |key,value|
  f = File::new(key+".dat","w")
  value.sort! { |x,y| 
#    puts x.inspect
#    puts y.inspect
    [ x[1], x[4] ] <=> [ y[1], y[4] ]
  }
  value.each { |val|
    f.puts "#{val[1]} #{val[4]} #{val[0]} #{val[5]/val[0]}"
  }
  f.close
}

