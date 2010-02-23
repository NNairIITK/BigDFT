require 'yaml'


MAX_N1=256
MAX_N2=256
MAX_N3=256
N1_STEP=64
N2_STEP=64
N3_STEP=64
N1_URANGE=8 
N2_URANGE=8 
N3_URANGE=0

sizes = []

for n1 in (N1_STEP..MAX_N1).step(N1_STEP)
  for un1 in (n1-N1_URANGE)..(n1+N1_URANGE)
    for n2 in (N2_STEP..MAX_N2).step(N2_STEP)
      for un2 in (n2-N2_URANGE)..(n2+N2_URANGE)
        for n3 in (n2..MAX_N3).step(N3_STEP)
          for un3 in (n3-N3_URANGE)..(n3+N3_URANGE)
            sizes.push( [un1, un2, un3 , un2*un3, un1*un2*un3] )
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
  results.push(h) if not h["method"].match(/memcpy/)
end
times = Hash::new { |hash,key| hash[key] = [] }
f.close

results.each { |res|
  times[res["method"]].push([ res["gputime"].to_f ])
}

times.each{ |key,value|
  for i in 0...sizes.length
    value[i].concat(sizes[i])
  end
}

times.each{ |key,value|
  f = File::new(key+".dat","w")
  value.each { |val|
    f.puts "#{val[1]} #{val[4]} #{val[0]}"
  }
  f.close
}

