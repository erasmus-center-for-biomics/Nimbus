
// STL headers
#include <cstdlib>
#include <string>
#include <iostream>

// SAM headers 
#include <faidx.h>

// Own headers
#include <refsequence.h>

namespace nimbus {

	GenomeSequence::GenomeSequence( ) {
		fai = NULL ;
		filename = "" ;
	}

	void GenomeSequence::set( std::string fn ) {
		fai = fai_load( fn.c_str() ) ;
		filename = fn ;
	}

	GenomeSequence::~GenomeSequence() {
		if( fai != NULL ) {			
			fai_destroy( fai ) ;			
		}
	}

	std::string GenomeSequence::get( const char* rname, int start, int end ) {
		
		// if the index is NULL return a single 
		if( fai == NULL )
			return std::string("N") ;

		// default return value is ""
		std::string rval = "" ;

		// use htslib to get the sequence
		int l   = 0 ;	
		char* s = faidx_fetch_seq( fai, rname, start, end, &l ) ;			

		// save the sequence in rval and free the c-string
		if( l > 0 ) {
			rval = std::string( s ) ;
			free(s) ;
		}

		// return the sequence
		return rval ;
	}

	std::string GenomeSequence::get( std::string rname, int start, int end ) {
		return get( rname.c_str(), start, end ) ;
	}

}