#include "stdafx.h"

#include "Utils.h"
#include "Read.h"

namespace Nimbus {

	namespace basic {


		//
		//
		//
		//
		//

		/* 
		 Default constructor
		 */
		Read::Read( std::string name, std::string sequence, std::string quality ) {
			_name = std::string(name)  ;
			_seq  = std::string(sequence) ;
			_qual = std::string(quality) ;

			_r_qual = "" ;
			_rc_seq = "" ;
		}

		Read::Read() {
			_name = "" ;
			_seq  = "" ;
			_qual = "" ;
			_r_qual = "" ;
			_rc_seq = "" ;
		}

		Read::Read( const Read& other ) {
			_name = other.name() ;
			_seq  = other.sequence() ;
			_qual = other.quality() ;
		}

		Read::~Read(void)
		{
		}

		// get the reverse complement of the sequence
		std::string Read::rc_sequence() {
			std::string rval ;

			// set the reverse complement read on the first call of this function
			if( _rc_seq.size() == 0 ) {

				// complement the read via a stringstream
				std::stringstream tmp ;
				for( std::string::reverse_iterator it=_seq.rbegin(); it!=_seq.rend(); ++it ) {
					tmp << utils::complement_base( *it ) ; 
				}
				_rc_seq = tmp.str() ;
			}
			rval = _rc_seq ;
			return rval ;
		}

		// get the reverse complement of the quality
		std::string Read::r_quality() {
			std::string rval ;

			// set the reverse complement read on the first call of this function
			if( _r_qual.size() == 0 ) {
				_r_qual = std::string( _qual.rbegin(), _qual.rend() ) ;
			} 
			rval = _r_qual ;
			return rval ;
		}

		// get the reverse complement for this read
		std::string Read::name() const {
			return _name ;
		}

		std::string Read::sequence() const {
			return _seq ;
		}

		std::string Read::quality() const {
			return _qual ;
		}

		std::string Read::fastq() const {
			std::stringstream rval ;
			rval << "@" << name() << std::endl 
				<< sequence() << std::endl 
				<< "+" << std::endl 
				<< quality() ; 
			return rval.str() ;
		}
	
		unsigned int Read::size() const {
			return (unsigned int)_seq.size() ;
		}

		std::string Read::str() const {
			std::stringstream s ;
			s << "" << name() << "-" << sequence () << "-" << quality() ;
			return s.str() ;
		}

		
		//
		// protected functions
		//
		void Read::sequence( std::string s ) {
			_seq    = std::string(s) ;
			_rc_seq = "" ;
		}

		void Read::quality( std::string q ) {
			_qual   = std::string(q) ;
			_r_qual = "" ;
		}

	}
}