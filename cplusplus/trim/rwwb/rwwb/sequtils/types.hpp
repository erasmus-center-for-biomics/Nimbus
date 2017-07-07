#ifndef sequtils_types_hpp
#define sequtils_types_hpp

//
// LICENSE BIOLERPLATE
//
// In this file, the basic types and converters are 
//   defined related to DNA sequences. 
//
// 

#include <cstddef> 
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>

namespace rwwb {

namespace sequtils {

	// bases are expressed as uint8_t types
	//  where ACGT is 0,1,2,3 and N is -1 
	typedef int8_t base_t ; 
	
	// Copy to kmer
	//
	//
	inline bool copy_to_kmer( std::vector<base_t>& target, const std::vector<base_t>& source, std::size_t offset, std::size_t size, bool stop_on_N){
		bool rval = false ;
		std::size_t i = offset ; 
		for(size_t j=0; j<size; ++j) {
			
			i = offset + j ;
			if(i>source.size()) {
				throw std::out_of_range("") ;
			}
			
			target[j] = source[i] ;

			if(target[j] < 0){
				rval = true ;

				// direct return on Ns
				if(stop_on_N){
					return rval ;
				}
			}
		}
		return rval ;
	}

	/*
IUPAC table

A	Adenine
C	Cytosine
G	Guanine
T (or U)	Thymine (or Uracil)
R	A or G
Y	C or T
S	G or C
W	A or T
K	G or T
M	A or C
B	C or G or T
D	A or G or T
H	A or C or T
V	A or C or G
N	any base
. or -	gap

	 */

	//
	// A converter function to convert characters to bases. 
	//
	//
	//
	inline base_t char_to_base( char x ) {
		base_t r = -1 ;
		switch(x) {
		case 'A':
		case 'a':
			r = 0 ;
			break ;
		case 'C':
		case 'c':
			r = 1 ;
			break ;
		case 'G':
		case 'g':
			r = 2 ;
			break ;
		case 'T':
		case 't':
			r = 3 ;
			break ;
		// extended IUPAC codes
		case 'R':
		case 'r':
			r = 4 ;
			break ;
		case 'Y':
		case 'y':
			r = 5 ;
			break ;
		case 'S':
		case 's':
			r = 6 ;
			break ;
		case 'W':
		case 'w':
			r = 7 ;
			break ;
		case 'K':
		case 'k':
			r = 8 ;
			break ;
		case 'M':
		case 'm':
			r = 9 ;
			break ;
		case 'B':
		case 'b':
			r = 10 ;
			break ;
		case 'D':
		case 'd':
			r = 11 ;
			break ;
		case 'H':
		case 'h':
			r = 12 ;
			break ;		
		case 'V':
		case 'v':
			r = 13 ;
			break ;		
		default: {
				r = -1 ;
				break ;
				 }
		} ;
		return r ;
	}

	
	//
	// A converter function to convert bases to characters. 
	//
	//
	//
	inline char base_to_char( base_t x ) {
		char r ;
		switch(x) {
		case 0:
			r = 'A' ;
			break ;
		case 1:
			r = 'C' ;
			break ;
		case 2:
			r = 'G' ;
			break ;
		case 3:
			r = 'T' ;
			break ;
		case 4:
			r = 'R' ;
			break ;
		case 5:
			r = 'Y' ;
			break ;
		case 6:
			r = 'S' ;
			break ;
		case 7:
			r = 'W' ;
			break ;
		case 8:
			r = 'K' ;
			break ;
		case 9:
			r = 'M' ;
			break ;
		case 10:
			r = 'B' ;
			break ;
		case 11:
			r = 'D' ;
			break ;
		case 12:
			r = 'H' ;
			break ;
		case 13:
			r = 'V' ;
			break ;
		default:
			r = 'N' ;
			break ;
		} ;
		return r ;
	}

	//
	// Generates a vector of bases from a string
	//
	inline std::vector<base_t> string_to_base( const std::string seq ) {
		std::vector<base_t> r = std::vector<base_t>( seq.size(), -1  ) ;
		std::transform( seq.begin(), seq.end(), r.begin(), char_to_base ) ;
		return r ;
	}

	//
	//
	//
	struct char_to_qual {
		char_to_qual( uint8_t offset ): offset(offset) {} 
		uint8_t operator()( char v ) const { 
			uint8_t dv = (uint8_t) v ;
			return dv > offset ? dv - offset : 0 ;
		} 
	private:
		uint8_t offset ;
	} ;

} ;
} ;

#endif