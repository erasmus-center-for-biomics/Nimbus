#ifndef sequtils_fastq_hpp
#define sequtils_fastq_hpp

// STL
#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

// own headers
#include <rwwb/sequtils/types.hpp> 

namespace rwwb {

namespace sequtils {
	
	//
	// A continuous reader with minimal copies
	//
	//
	//

	typedef struct __read__ {
		std::size_t uid ;
		std::string name ;
		std::vector<base_t> sequence ;
		std::string quality ;
	} read ;

	//
	// Get multiple reads from a FastQ file. Rhe number of reads 
	//   to obtain is determined  by the size of reads.
	//
	// 

	std::size_t reads_from_fastq(std::istream& handle, std::size_t& counter, std::vector<read>& reads) {
		
		// the number of reads obtained
		std::size_t rval = counter ;
		
		// the l variable
		std::string l = "" ;

		for(std::size_t i=0; i<reads.size(); ++i) {
			
			// read header
			if(!std::getline(handle, l).good())
				break ;
			
			reads[i].name = std::string(l.substr(1)) ;
			
			// sequence
			if(!std::getline(handle, l).good())
				break ;

			reads[i].sequence.resize(l.size()) ;
			std::transform(l.begin(), l.end(), reads[i].sequence.begin(), char_to_base ) ;

			// +
			if(!std::getline(handle, l).good())
				break ;

			if(!std::getline(handle, l).good())
				break ;
				
			reads[i].quality = std::string(l) ;

			// set the uid
			counter += 1 ;
			reads[i].uid = counter ;
		}

		// return the number of reads from a fastq.
		return counter - rval ;
	}

} ;
} ;
#endif