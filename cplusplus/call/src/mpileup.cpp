// Standard libraries
#include <cstdlib>
#include <iostream>

// SAM libraries
#include <sam.h>

// Own headers
#include <mpileup.h>
#include <read_provider.h>

namespace nimbus {

	Mpileup::Mpileup( ) {
		iter  = NULL ;
		tid   = 0 ;
		pos   = 0 ;
		plp   = NULL ;
		n_plp = NULL ;

		maximum_depth   = 10000 ;		
	}

	Mpileup::~Mpileup( ) {		
		if( iter != NULL ) {
			std::cerr << "iterator" << std::endl ;
			free(iter) ;
			iter = NULL ;
		}		
		if( plp != NULL ) {
			std::cerr << "pileup" << std::endl ;
			free( plp ) ;	
			plp = NULL ;
		}		
		if( n_plp != NULL ) {
			std::cerr << "count" << std::endl ;
			free(n_plp) ;
			n_plp = NULL ;
		}						
	}

	void Mpileup::initialize( Provider** P, int n ) {
		
		// Prepare the iterators over the 
		iter = bam_mplp_init( n, ReadProvider, (void**) P ) ;		
		bam_mplp_set_maxcnt( iter, (int) maximum_depth ) ;

		// Allocate memory for the results
		plp   = (const bam_pileup1_t**) calloc( n, sizeof(bam_pileup1_t*) ) ;
		n_plp = (int*) calloc(  n, sizeof(int) ) ;

	}

	bool Mpileup::next() {

		//
		int ret = bam_mplp_auto(iter, &tid, &pos, n_plp, plp ) ;
		
		// 
		return ret > 0 ? true : false ;
	}

} 