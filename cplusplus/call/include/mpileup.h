#pragma once 

// Standard libraries
#include <cstdlib>

// SAM libraries
#include <sam.h>

// Own headers
#include <read_provider.h>


namespace nimbus {

	/**
	 * A class to wrap the samtools bam_mplp functions 
	 *  in a nice C++ style object.
	 */
	class Mpileup {
				
		// the mpileup iterator
		bam_mplp_t iter ;

	public:
		// refreshed with each iteration
		int tid ;
		int pos ;
		const bam_pileup1_t **plp ;
		int *n_plp ;

	// set before performing init
	public: 
		std::size_t maximum_depth ;
		
	public:

		Mpileup(  ) ;

		~Mpileup() ;

		void initialize( Provider** P, int n ) ;

		/**
		 * Iterates to the next base
		 *
		 */
		bool next() ;

	} ;
}