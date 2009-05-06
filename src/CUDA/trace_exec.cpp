
#include "trace_exec.h"

trace_exec::trace_exec(const std::string& traceFile_)
{
  traceFile = new std::ofstream(traceFile_.c_str());
  nextTraceNumber = 0;
}
trace_exec::~trace_exec()
{
  delete traceFile;
}



write_trace trace_exec::getNewWriteTrace(trace_t tr,const std::string& name)
{
  return write_trace(tr,nextTraceNumber++,name,this);
 
}


void trace_exec::writeToFile(trace_t tr,
			     int traceNumber,
			     time_mesure time, 
			     const std::string& name)
{
  std::string type;
  switch(tr)
    {
    case GPU_TRSF:
      type = "GPU_TRSF";
      break;
    case GPU_CALC:
      type = "GPU_CALC";
      break;
    case CPU_MEMCPY:
      type = "CPU_MEMCPY";
      break;
    default:
      type = "DEFAULT";
      
    }

  (*traceFile) << type << "," 
	       << traceNumber << "," 
	       << time << "," 
	       << name << std::endl;
}


// ============== write_trace ==========

time_mesure write_trace::getTime() const
{
  timeval tim;
  gettimeofday(&tim,NULL);
  

  return  tim.tv_usec + tim.tv_sec * 1000000L;
}


void write_trace::writeToFile(time_mesure usec)
{
  te->writeToFile(trace,traceNum,usec,name);
}


void write_trace::begTimeCount()
{
  writeToFile(getTime());
}


void write_trace::endTimeCount()
{
writeToFile(getTime());

}


