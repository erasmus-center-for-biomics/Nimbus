
// header to implement
#include <aspects.h>

// STL
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <limits>

namespace nimbus {

	//
	//
	// STL type based 
	//
	//

	int compare( 
		std::vector<std::string> a,
		std::vector<std::string> b  ) {

		// get the comparator size
		std::size_t s = a.size() > b.size() ? b.size() : a.size() ;

		// check that the values are not larger in the 
		// current than in the oth object 
		for( std::size_t i=0; i<s; ++i ) { 
			int c = a[i].compare(0, std::string::npos, b[i] ) ;
			// int c = cmp( a[i], b[i] ) ;
			if( c != 0 ) {
				return c ;
			}
		}
		
		//
		// we reached this far, so 
		// all values up to msize are equal
		// 
		// now the shortest array is defined as 
		// less, the longer is more. If both are 
		// the same length define them as equal
		if( a.size() == b.size() ) {
			return 0 ;
		} else if( a.size() < b.size() ) {
			return -1 ;
		} else {
			return 1 ;
		}
	}

	int compare( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) {
		
		//
		return compare( a.first, b.first ) ;		
	}

	bool lessthan( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) {
		int c = compare( a.first, b.first ) ;
		return c < 0 ? true : false ;
	}

	bool morethan( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) {
		int c = compare( a.first, b.first ) ;
		return c > 0 ? true : false ;
	}

	bool equals( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) {
		int c = compare( a.first, b.first ) ;
		return c == 0 ? true : false ;
	}

}
