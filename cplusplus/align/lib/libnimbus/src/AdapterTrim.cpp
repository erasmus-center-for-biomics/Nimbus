#include "stdafx.h"
#include "AdapterTrim.h"


namespace Nimbus {

	namespace utils {
		AdapterTrim::AdapterTrim( std::string sequence, int seedsize ) {
			_seq  = sequence ;
			_seed = sequence.substr( 0, seedsize ) ;
		}


		AdapterTrim::~AdapterTrim(void) {
		}

		basic::Read* AdapterTrim::trim( basic::Read* r ) {
			basic::Read* rval = NULL ;

			//
			std::string oseq = r->sequence() ;			
			std::string oqua = r->quality() ;
			
			// trim the sequence 
			std::string tseq = trim( oseq ) ;
			std::string tqua = oqua ; 

			// std::cout << tseq << std::endl ;

			// if we did not trim away the entire adapter
			if( tseq.size() > 0 ) {
				if( tseq.size() < oseq.size() ) tqua = oqua.substr( 0, tseq.size() ) ;
			} else {
				// if we did, fill the read with N characters (some 
				// aligners do not allow empty sequences)
				tseq = std::string( oseq.size(), 'N' ) ;
			}

			// construct a new read
			rval = new basic::Read( r->name(), tseq, tqua ) ; 

			return rval ;
		}

	
		std::string AdapterTrim::trim( std::string seq ) {
			return trim( seq, std::string("") ) ;
		}

		std::string AdapterTrim::trim( std::string seq, std::string pre ) {
			// declare the return value, we will modify this only if required
			std::string rval = pre + seq ;
			size_t pos = -1 ;
		
			// find the seed in the sequence 
			pos = seq.find( _seed ) ;

			// if we found a seed and the string is not empty
			if( pos != std::string::npos ) {

				// get the possible match
				std::string pseq  = seq.substr( pos ) ;

				// get the matching portion of the query
				std::string query = _seq ;
				if( pseq.size() > query.size() ){
					pseq = pseq.substr( 0, query.size() ) ;
				} else if (pseq.size() < query.size()) {
					query = query.substr( 0, pseq.size() ) ;
				}

				// test whether the query and the sequence are the same
				if( pseq == query ) {
					rval = pre + seq.substr(0, pos ) ;
				} else {
					// if not move the seed-length positions to the left
					pos += _seed.size() ;
					// recursively search the rest of the string
					if( pos < seq.size() ) {
						rval = trim( seq.substr( pos ), pre + seq.substr(0, pos ) ) ;  
					}
				}
			}
			// return the trimmed sequence
			return rval ;
		}
	}
}