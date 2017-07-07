
// C standard library
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>

// Own headers
#include <cigar.h>
#include <sam.h>
#include <alignment_functions.h>

//
namespace nimbus {

	//
	// Constructors and destructors
	//

	Cigar::Cigar() {
		_segments = NULL ;
		_n = 0 ;
	}

	Cigar::~Cigar(){
		if( _segments != NULL ) {
			delete[] _segments ;
		}
	}

	
	void Cigar::reset( ) {
		_n = 0 ;
		if( _segments != NULL ) {
			delete[] _segments ;
		}
	}

	void Cigar::from_hts_cigar( uint32_t* cig, int clen ) {
		
		reset() ;

		// expand the number of operations
		expand( (std::size_t) clen )  ;

		//  add the operations
		for( std::size_t i=0; i<(std::size_t)clen; ++i ) {

			// interpret the cigar
			int v = bam_cigar_op( cig[i] ) ;
			std::size_t l = (std::size_t) bam_cigar_oplen( cig[i] ) ;

			//
			_segments[i].set( int2op(v), l ) ;
		}
	}
	

	//
	// Utility
	//


	std::size_t Cigar::query_length( ) const {
		std::size_t rval = 0 ;
		for( std::size_t i=0; i<_n; ++i ) {
			char op = _segments[i].operation() ;
			if( on_query(op) ) {
				rval += _segments[i].runLength() ;
			}						
		}
		return rval ;
	}

	std::size_t Cigar::reference_length( ) const {
		std::size_t rval = 0 ;
		for( std::size_t i=0; i<_n; ++i ) {
			char op = _segments[i].operation() ;
			if( on_reference(op) ) {
				rval += _segments[i].runLength() ;
			}			
		}
		return rval ;
	}

	bool Cigar::only_on_query( char op ) const {
		bool rval = false ;
		switch( op ) {		
		case 'I':
			rval = true ;
			break ;
		case 'S':
			rval = true ;
			break ;		
		}
		return rval ;
	}


	bool Cigar::on_query( char op ) const {
		bool rval = false ;
		switch( op ) {
		case 'M':
			rval = true ;
			break ;
		case 'I':
			rval = true ;
			break ;
		case 'S':
			rval = true ;
			break ;
		case '=':
			rval = true ;
			break ;
		case 'X':
			rval = true ;
			break ;
		}
		return rval ;
	}

	bool Cigar::only_on_reference( char op ) const {
		bool rval = false ;
		switch( op ) {		
		case 'D':
			rval = true ;
			break ;			
		case 'N':
			rval = true ;
			break ;					
		}
		return rval ;
	}


	bool Cigar::on_reference( char op ) const {
		bool rval = false ;
		switch( op ) {
		case 'M':
			rval = true ;
			break ;
		case 'D':
			rval = true ;
			break ;			
		case 'N':
			rval = true ;
			break ;			
		case '=':
			rval = true ;
			break ;
		case 'X':
			rval = true ;
			break ;
		}
		return rval ;
	}

	std::size_t Cigar::length() const {
		return _n ;
	}

	/*
	void Cigar::query_at( std::size_t query, std::size_t &i_operation, std::size_t &pos_in_op ) const {
		
		//
		std::size_t crd = 0 ;
		i_operation = _n + 1 ;
		pos_in_op   = 0 ;

		for( std::size_t i=0; i<_n; ++i ) {

			// if the data is on the query
			if( on_query( _segments[i].operation() ) ) {

				// increase the coordinate	
				crd += _segments[i].runLength() ;

				//
				if( crd >= query ) {
					i_operation = i ;
					pos_in_op   = _segments[i].runLength() - (crd - query) ;
					break ;
				}
			}
		}
	}
	*/
	
	/*
	void Cigar::reference_at( std::size_t r, std::size_t &i_operation, std::size_t &pos_in_op ) const {

		//
		std::size_t crd  = 0 ;
		i_operation      = _n + 1 ;
		pos_in_op        = 0 ;

		for( std::size_t i=0; i<_n; ++i ) {

			// if the data is on the query
			if( on_query( _segments[i].operation() ) ) {

				// increase the coordinate	
				crd += _segments[i].runLength() ;

				// 
				if( crd >= r ) {
					i_operation = i ;
					pos_in_op   = _segments[i].runLength() - (crd - r) ;
					break ;
				}
			}
		}		
	}
	*/

	std::size_t Cigar::reference_offset_to_query_offset( std::size_t relpos, std::size_t &i_operation, std::size_t &pos_in_op  ) const {

		// set the coordinates to account for
		std::size_t crd  = relpos ;
		
		// set default values for the returns
		std::size_t rval = 0 ;
		pos_in_op   = 0 ;
		i_operation = _n + 1 ;

		/* std::cout << "New" << std::endl ; */

		// iterate over the cigar instances
		for( std::size_t i=0; i<_n; ++i ) {
			
			// if the operation should account for the relative position
			if( on_reference( _segments[i].operation() ) ) {

				// if we can fully account for the coordinates
				if( _segments[i].runLength() - 1 < crd ){ 			
					crd -= _segments[i].runLength() ;			
				} else {

					// we will only exit on a reference base
					pos_in_op   = crd ;
					i_operation = i ;

					// only increase the return value if the base is on the query
					if( on_query( _segments[i].operation() ) ){
						rval += crd  ;
					}

					/* // debug message
					std::cout << "cycle: " << i 
						<< ", left: " << crd 
						<< ", position: " << rval 
						<< ", operation: " << _segments[i].operation() 
						<< ", length: " << _segments[i].runLength() 
						<< std::endl ;
					*/

					// break the loop
					break ;								
				}
			}

			// increase the query position
			if( on_query( _segments[i].operation() ) ) {
				rval += _segments[i].runLength() ;
			}

			/* // debug message
			std::cout << "cycle: " << i 
				<< ", left: " << crd 
				<< ", position: " << rval 
				<< ", operation: " << _segments[i].operation() 
				<< ", length: " << _segments[i].runLength() 
				<< std::endl ;
			*/
		}		
		return rval ;
	}

	std::size_t Cigar::reference_offset_to_query_offset( std::size_t relpos ) const {		
		std::size_t i_operation ; 
		std::size_t pos_in_op ;
		
		return reference_offset_to_query_offset( relpos, i_operation, pos_in_op ) ;		
	}

	CigarOperation* Cigar::at( std::size_t i ) const {
		if( i < _n ) {
			return _segments + i ;
		}
		return NULL ;
	}

	//
	// Protected functions
	//

	void Cigar::expand( std::size_t extra ) {

		// 
		CigarOperation* tmp = new CigarOperation[ _n + extra ] ;
		
		for(std::size_t i=0; i<_n; ++i ) {
			tmp[i].set( _segments[i].operation(), _segments[i].runLength() ) ;
		}

		// delete the segments if present
		if( _segments != NULL ) 
			delete[] _segments ;
		
		// 
		_segments = tmp ;
		_n       += extra ;
	}

	char Cigar::int2op( int v ) {
		char op = '?' ;

		switch(v) {
		case 0:
			op = 'M' ;
			break ;
		case 1:
			op = 'I' ;
			break ;
		case 2:
			op = 'D' ;
			break ;
		case 3:
			op = 'N' ;
			break ;
		case 4:
			op = 'S' ;
			break ;
		case 5:
			op = 'H' ;
			break ;
		case 6:
			op = 'P' ;
			break ;
		case 7:
			op = '=' ;
			break ;
		case 8:
			op = 'X' ;
			break ;
		case 9:
			op = 'M' ;
			break ;
		default:
			op = '?' ;
			break ;
		}
		return op ;
	}

	std::string Cigar::str() const {
		std::stringstream s ;
		for( std::size_t i=0; i<length(); ++i ) {
			s << _segments[i].str() ;
		}
		return s.str() ;
	}

} 