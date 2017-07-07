
#include <sample.h>


// STL imports 
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace nimbus {


	void parseSamplesFromHeader( 
		char* text, 
		std::size_t nchar, 
		std::vector<std::pair<std::string,std::string> > &samples
		) {
		
		//
		std::string line = "" ;
		std::stringstream stream( text ) ;
		std::pair<std::string, std::string> tmp = std::pair<std::string, std::string>() ; 
	
		// while we can pull lines from the stream
		while( std::getline(stream, line, '\n') ) {

			// parse @RG lines
			if( line.substr(0, 3) == "@RG" ) {

				// add a new Sample				
				bool rval = fromRGline( line, tmp ) ;
				if( rval ) {
					samples.push_back( tmp ) ;
				}
			}
		}
	}

	bool fromRGline( 
		std::string text, 
		std::pair< std::string, std::string > &sample ) {

		//
		std::size_t it   = 0 ;
		std::string line = std::string( text ) ;

		// Set the sample name
		it = line.find( "SM:" ) ;
		if( it == std::string::npos ) {
			// throw std::runtime_error( "BAM header line does not contain a sample ID" ) ;
			return false ;
		}

		std::string sname = line.substr( it + 3 ) ;
		it = sname.find( "\t" ) ;
		if( it != std::string::npos ) {
			sname = sname.substr( 0, it ) ;
		}
		
		sample.first = sname ;

		// Set the readgroup id
		it = line.find( "ID:" ) ;
		if( it == std::string::npos ) {
			//throw std::runtime_error( "BAM header line does not contain a readgroup identifier" ) ;
			return false ;
		}
		std::string rgid = line.substr( it + 3 ) ;
		it = rgid.find( "\t" ) ;
		if( it != std::string::npos ) {
			rgid = rgid.substr( 0, it ) ;
		}

		sample.second = rgid ;
		return true ;
	}

	bool eq( 
		std::pair< std::string, std::string > a, 
		std::pair< std::string, std::string > b 
		) {		
		if( a.first.compare( b.first )  != 0 ) {
			return false ;
		}

		if( a.second.compare( b.second )  != 0 ) {
				return false ;
		}
		return true ;
	}


}
