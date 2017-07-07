#pragma once

// STL headers
#include <cstdlib>
#include <string>

// SAM headers 
#include <faidx.h>

namespace nimbus {


	class GenomeSequence {

		//
		faidx_t *fai ;
		
	public:
		
		std::string filename ;

		/**
		 * Constructs a new GenomeSequence object that 
		 *  will only yield N as sequence
		 */
		GenomeSequence( ) ;

		/**
		 * Sets the FastA file from which to obtain
		 *  sequences
		 *
		 * @param fn - the name of the samtools indexed FastA file
		 */
		void set( std::string fn ) ;

		/**
		 * Destroys the GenomeSequence object
		 */
		~GenomeSequence() ;

		/**
		 * Get the sequence from sequence rname from start to end
		 *
		 * @param rname - the sequence name
		 * @param start - the start of the region
		 * @param end   - the end of the region 
		 *
		 * @returns: the sequence obtained
		 */
		std::string get( const char* rname, int start, int end ) ;

		std::string get( std::string rname, int start, int end ) ;
	} ;
	 
}