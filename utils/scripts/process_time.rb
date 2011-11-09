require 'yaml'

CLASSES = [ "Communications",
    "Convolutions",   
    "Linear Algebra", 
    "Other",          
    "Potential",      
#    "Initialization", 
#    "Finalization",   
    "Total"] 

CATEGORIES = ["ApplyProj"]

#find the reference number of cores and the total time (First run of first file)
docs = YAML::load_stream( File::read(ARGV.first) ).documents
nprocref=docs.first["SUMMARY"]["CPU Parallelism"]["MPI procs"]
nthdsref=docs.first["SUMMARY"]["CPU Parallelism"]["OMP thrds"]
timeref= docs.first["WFN_OPT"]["Classes"]["Total"][1]

print "#ncores nproc nthds" 
CLASSES.each { |cls|
  printf(" %#{cls.length<9 ? 9 : cls.length}s", cls)
}
print ' | References: Time, Procs'
puts
puts "# Reference file (#{ARGV.first}):"
nproc = docs.first["SUMMARY"]["CPU Parallelism"]["MPI procs"]
nthds = docs.first["SUMMARY"]["CPU Parallelism"]["OMP thrds"]
printf( " %#{'ncores'.length}d %#{'nproc'.length}d %#{'nthds'.length}d",(nthds==0 ? nproc : nthds*nproc),nproc, nthds)
CLASSES.each { |cls|
  printf( " %#{cls.length<9 ? 9 : cls.length}.2e",docs.first["WFN_OPT"]["Classes"][cls][1])
}
printf( " %#{'| Ref.: Time '.length}.2e %5d %5d",timeref, nprocref,nthdsref)
puts

ARGV.each { |arg|
  puts "# #{arg}:"
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
# To extract detailed categories
#    CATEGORIES.each { |cls|
#      printf( " %#{cls.length<9 ? 9 : cls.length}.2e",doc["WFN_OPT"]["Categories"][cls]["Data"][1])
#    }
    printf( " %#{'| Ref.: Time '.length}.2e %5d %5d",timeref, nprocref,nthdsref)
    puts
  }
}

