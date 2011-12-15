require 'yaml'

#Classes which have to be analyzed
CLASSES = [ "Communications",
    "Convolutions",   
    "Linear Algebra", 
    "Other",          
    "Potential",      
#    "Initialization", 
#    "Finalization",   
    "Total"] 

# Subcategories 
CATEGORIES = []#["ApplyProj"]

#print the legend of the results
print "#ncores nproc nthds" 
CLASSES.each { |cls|
  printf(" %#{cls.length<9 ? 9 : cls.length}s", cls)
}
CATEGORIES.each { |cat|
  printf(" %#{cat.length<9 ? 9 : cat.length}s", cat)
}
print ' | References: Time, Procs'
puts


#parse the first argument
docs = YAML::load_stream( File::read(ARGV.first) ).documents

puts "# Reference file (#{ARGV.first}):"

#find the reference number of cores and the total time (First run of first file)
nprocref=docs.first["SUMMARY"]["CPU Parallelism"]["MPI procs"]
nthdsref=docs.first["SUMMARY"]["CPU Parallelism"]["OMP thrds"]
timeref= docs.first["WFN_OPT"]["Classes"]["Total"][1]



nproc = docs.first["SUMMARY"]["CPU Parallelism"]["MPI procs"]
nthds = docs.first["SUMMARY"]["CPU Parallelism"]["OMP thrds"]
printf( " %#{'ncores'.length}d %#{'nproc'.length}d %#{'nthds'.length}d",(nthds==0 ? nproc : nthds*nproc),nproc, nthds)
CLASSES.each { |cls|
  printf( " %#{cls.length<9 ? 9 : cls.length}.2e",docs.first["WFN_OPT"]["Classes"][cls][1])
}
CATEGORIES.each { |cls|
  printf( " %#{cls.length<9 ? 9 : cls.length}.2e",docs.first["WFN_OPT"]["Categories"][cls]["Data"][1])
}
printf( " %#{'| Ref.: Time '.length}.2e %5d %5d",timeref, nprocref,nthdsref)
puts
##end of reference file writing

##start with writing of all the files
ARGV.each { |arg|
  puts "# #{arg}:"

  #read each file
  docs = YAML::load_stream( File::read(arg) ).documents
  if docs.length>1
    puts "#Read #{docs.length} YAML documents."
  end
 
  docs.each { |doc|
    nproc = doc["SUMMARY"]["CPU Parallelism"]["MPI procs"]
    nthds = doc["SUMMARY"]["CPU Parallelism"]["OMP thrds"]
    printf( " %#{'ncores'.length}d %#{'nproc'.length}d %#{'nthds'.length}d",(nthds==0 ? nproc : nthds*nproc),nproc, nthds)
    CLASSES.each { |cls|
      printf( " %#{cls.length<9 ? 9 : cls.length}.2e",doc["WFN_OPT"]["Classes"][cls][1])
    }
    CATEGORIES.each { |cls|
      printf( " %#{cls.length<9 ? 9 : cls.length}.2e",doc["WFN_OPT"]["Categories"][cls]["Data"][1])
    }
    printf( " %#{'| Ref.: Time '.length}.2e %5d %5d",timeref, nprocref,nthdsref)
    puts


    #write the hostname associated to the iproc=0
#    iproc=0
#    hostnames=doc["SUMMARY"]["CPU Parallelism"]["Hostnames"]
#    puts "# Hostame for processor #{iproc}: #{hostnames[iproc]}" if hostnames


  }
}

