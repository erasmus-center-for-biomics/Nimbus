#pragma once 

// STL imports 
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

// SAM tools
#include <sam.h>

namespace nimbus {
	
	/**
	 * Parses the sample information from the BAM header
	 *
	 */ 
	void parseSamplesFromHeader( 
		char* text, 
		std::size_t nchar, 
		std::vector<std::pair<std::string,std::string> > &samples
		) ;

	bool fromRGline( 
		std::string text, 
		std::pair< std::string, std::string > &sample 
		) ;
	/*
	 *
	 *
	 */
	bool eq( 
		std::pair< std::string, std::string > a, 
		std::pair< std::string, std::string > b 
		) ;


}
