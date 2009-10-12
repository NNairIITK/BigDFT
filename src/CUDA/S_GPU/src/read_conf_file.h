/****c* CUDA/read_conf_file.h
** 
** AUTHOR
**  Matthieu Ospici
** 
** CHANGELOG
**  Started on  Wed May 21 14:22:40 2008 Matthieu Ospici
**
**  Last update Wed May 21 14:22:40 2008 Matthieu Ospici
**
** SOURCE
*/

#ifndef   	READ_CONF_FILE_H_
#define   	READ_CONF_FILE_H_


#include <map>
#include <string>

#include "read_conf_exception.h"
#include "class_utils.h"

typedef std::map<std::string , std::string> mapFile_t;








class readConfFile
{
public:
  readConfFile(const std::string& filename) throw(file_not_found);
  void get(const std::string& key, std::string& value) const throw(read_not_found);
  void get(const std::string& key, int *value) const throw(read_not_found);
  
private:
  mapFile_t mFile;
};






class readConfFileGPU_CPU : public readConfFile
{
public:
  readConfFileGPU_CPU( const std::string& filename) : readConfFile(filename){};

  int getGPU(int MPI_ID) const throw(read_not_found_GPU);
  int getCPU(int MPI_ID) const throw(read_not_found_CPU);
  int getFlag(int MPI_ID) const throw(); //0 CUDA, 1 BLAS
};

#endif 	    /* !READ_CONF_FILE_H_ */

/****/
