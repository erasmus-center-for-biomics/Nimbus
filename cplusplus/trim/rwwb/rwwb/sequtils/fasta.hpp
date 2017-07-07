#ifndef sequtils_fasta_hpp
#define sequtils_fasta_hpp

// STL
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

// own headers
#include <rwwb/sequtils/types.hpp> 

namespace rwwb {

namespace sequtils {
	//
	// A functor to parse FastA files. 	
	//
	// a functor allows us to keep the last line in memmory without 
	//  the overhead of making a full fledged class. This seems a bit 
	//  wacky, but seems to be a good use case for a functor.
	//
	struct fasta {
		// initializer
		fasta(): last_line(""){}
		
		//  
		bool operator()( std::istream &h, std::string &nm, std::vector<base_t>& seq ) {
			
			//default states of theobjects 
			nm.clear() ;
			seq.clear() ;

			// if we can get data from the handle
			while( h.good() ) {
			
				// if we have a line left over from the last 
				// round, keep that because it should be the 
				// header  
				if( last_line == "" ) {
					std::getline( h, last_line ) ;
				}

				// we found a new header and we did not encounter one yet
				if( nm == "" && last_line != "" && last_line[0] == '>' ) {
					nm = last_line.substr(1) ;
					last_line = "" ;

				// we found a header after we encountered one already and should return the data
				} else if(nm != "" && last_line != "" && last_line[0] == '>') {
					return seq.size() > 0 ? true : false ;
				
				// if we already found an sequence identifier, we can add 
				// bases to seq
				} else if( nm != "" ) {
					
					// make sure we enough space to put our extra base
					auto offset = seq.size() ;
					seq.resize( seq.size() + last_line.size() ) ; 

					// add the bases to the sequence
					std::transform( last_line.begin(), last_line.end(), seq.begin() + offset, char_to_base ) ;
					last_line = "" ;
				}
			}

			// if we could not get further 
			// sequences return false 
			return seq.size() > 0 ? true : false ; ;
		}
	private:
		std::string last_line ;
	} ;
} ;
} ;

#endif