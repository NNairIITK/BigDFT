#ifndef __traceexech__
#define __traceexech__


#include <string>
#include <iostream>
#include <fstream>

#include <sys/time.h>
#include <time.h>

typedef  long long time_mesure;

class trace_exec;

enum trace_t {GPU_TRSF, GPU_CALC, CPU_MEMCPY};
class write_trace
{
public:
  write_trace(trace_t trace_,int traceNum_,const std::string& name_, trace_exec *te_)
    :trace(trace_),
     traceNum(traceNum_),
     name(name_),
     te(te_){}

  void begTimeCount();
  void endTimeCount();

  
private:
  void writeToFile(time_mesure);
  
  
  const trace_t trace;
  const int traceNum;
  const std::string name;

  time_mesure getTime() const;
  trace_exec *te;
};


class trace_exec
{
public:
  trace_exec(const std::string&,bool enableTrace);
  ~trace_exec();

  write_trace getNewWriteTrace(trace_t,const std::string& name );
  
 
private:
  void writeToFile(trace_t,int traceNumber,time_mesure time, const std::string& name);
  std::ofstream *traceFile;
  
  bool enableTrace;
  
  int nextTraceNumber;
  friend class write_trace;
};

#endif
