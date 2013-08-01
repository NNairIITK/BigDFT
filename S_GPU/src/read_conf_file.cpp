#include <fstream>
#include <sstream>

#include "read_conf_file.h"

readConfFile::readConfFile ( const std::string& filename ) throw ( file_not_found )
{
	std::ifstream configfile ( filename.c_str() ) ;
	if ( !configfile )
		throw file_not_found ( filename );
	std::string s ;
	while ( configfile )
	{
		getline ( configfile, s ) ;
		if ( s.size() == 0 || s[0] == '#' )
			continue ; // Skip empty lines and comment
		int offset = s.find ( '=' );
		std::string key ( s.substr ( 0, offset ) ) ;
		std::string value ( s.substr ( offset + 1, s.size() ) ) ;
		mFile[key] = value ;
	}
}

void readConfFile::get ( const std::string& key, std::string& value ) const throw ( read_not_found )
{
	mapFile_t::const_iterator found = mFile.find ( key ) ;
	if ( found == mFile.end() ) // not found
		throw read_not_found ( key ) ;
	else
		value = ( *found ).second ;
}

void readConfFile::get ( const std::string& key, int *value ) const throw ( read_not_found )
{
	std::string tmp;
	get ( key, tmp );
	*value =  strTo<int> ( tmp );
}

//========================================================================

int readConfFileGPU_CPU::getGPU ( int MPI_ID ) const throw ( read_not_found_GPU )
{
	std::ostringstream oss;
	oss << "ID_" <<  MPI_ID << "_GPU";
	std::string strSearch = oss.str();
	std::string searchRes;
	try
	{
		get ( strSearch, searchRes );
		return  strTo<int> ( searchRes );
	}
	catch ( read_not_found e )
	{
		throw read_not_found_GPU ( e.what() );
	}
}

int readConfFileGPU_CPU::getCPU ( int MPI_ID ) const throw ( read_not_found_CPU )
{
	std::ostringstream oss;
	oss << "ID_" <<  MPI_ID << "_CPU";
	std::string strSearch = oss.str();
	std::string searchRes;
	try
	{
		get ( strSearch, searchRes );
		return  strTo<int> ( searchRes );
	}
	catch ( read_not_found e )
	{
		throw read_not_found_CPU ( e.what() );
	}
}

//1 CUDA, 0 CUDABLAS
int readConfFileGPU_CPU::getFlag ( int MPI_ID ) const throw()
{
	std::ostringstream oss;
	oss << "ID_" <<  MPI_ID << "_FLAG";
	std::string strSearch = oss.str();
	std::string searchRes;
	try
	{
		get ( strSearch, searchRes );
		return  strTo<int> ( searchRes );
	}
	catch ( read_not_found e )
	{
		return 1; //default value : return 1 for CUDA
	}
}
