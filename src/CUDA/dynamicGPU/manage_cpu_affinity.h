#ifndef __managecpuaff__
#define __managecpuaff__


#include <vector>
#include <list>

#include <string>

//define the ID of CPUs connected to one GPU
class gpu_cpus_connexion
{
public:

 

  //  gpu_cpus_connexion():num_connected(0){}
  gpu_cpus_connexion(const std::string& st);

  void add_cpu_from_str(const std::string&); // string syntax : "cpu1,cpu2,....,cpun

  void add_cpu(int cpuNum);

 


  int set_affinity(int cpuID) const; //set current tread/process affinity with an avaible cpu number 

  int get_num_cpus() const {return cpus_connected.size();}
private:
  std::vector<int> cpus_connected;

 

 
};


//======================================

class manage_cpu_affinity
{
public:

  manage_cpu_affinity(int _iproc):iproc(_iproc){}

  void add_connexion(const gpu_cpus_connexion&);
  

  //set the calling process affinity to a processor that we have to connect to gpu_to_atach GPU
  void set_affinity(int cpu_aff) const;

private:

  std::list<gpu_cpus_connexion> connexions;

  int iproc;
};

#endif
