#pragma once

// STL
#include <cstdlib>
#include <cstddef>
#include <string>
#include <sstream>
#include <cstdint>

namespace nimbus {

	/**
	 * A CIGAR operation
	 *
	 */
	class CigarOperation {
		char _operation ;
		std::size_t _runlength ;

	public:
		CigarOperation() {
			_operation = '?' ;
			_runlength = 0 ;
		}

		CigarOperation( char o, std::size_t n ) {
			set( o, n ) ;
		}

		void set( char o, std::size_t n ) {
			_operation = o ;
			_runlength = n ;
		}

		std::size_t runLength() const {
			return _runlength ;
		}

		char operation() const {
			return _operation ;
		}

		std::string str() const {
			std::stringstream s ;
			s << operation() << runLength() ;
			return s.str() ;
		}
	} ;

	/**
	 * A CIGAR 
	 *
	 */
	class Cigar {

		CigarOperation* _segments ;
		std::size_t _n ;

	public:
		Cigar() ;		

		~Cigar() ;

		void reset( ) ;

		void from_hts_cigar( uint32_t* content, int clen ) ;

		std::size_t query_length( ) const ;

		std::size_t reference_length( ) const ;

		std::size_t length( ) const ;

		std::string str() const ;

		std::size_t reference_offset_to_query_offset( std::size_t r ) const ;
		
		std::size_t reference_offset_to_query_offset( std::size_t relpos, std::size_t &i_operation, std::size_t &pos_in_op  ) const ;

		/**
		 * get the operation at query coordinate 
		 *
		 */
		// void query_at( std::size_t q, std::size_t &i_operation, std::size_t &pos_in_op ) const ;

		// void reference_at( std::size_t r, std::size_t &i_operation, std::size_t &pos_in_op ) const ;
		
		/**
		 * Get a pointer to the cigar operation at i or NULL
		 *
		 */
		CigarOperation* at( std::size_t i ) const ;

	protected:

		void expand( std::size_t extra ) ;

		char int2op( int v ) ;

		bool on_reference( char op ) const ; 

		bool only_on_reference( char op ) const ; 		

		bool on_query( char op ) const ;

		bool only_on_query( char op ) const ; 
	} ;

}