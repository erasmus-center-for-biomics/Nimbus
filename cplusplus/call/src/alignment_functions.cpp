
// Standard library functions
#include <cstdlib>
#include <string>
#include <sstream> 
#include <iostream>
#include <vector>
#include <algorithm>

// HTS lib
#include <sam.h>

// own headers
#include <alignment_functions.h>
#include <cigar.h> 
#include <sample.h>

// all our functions live in the nimbus namespace
namespace nimbus {


	/**
	 * Get the sequence from an alignment 
	 *
	 */
	std::string alignment_sequence( const bam1_t* alignment, int start, int end ) {		
		//
		std::stringstream rval ;		
		
		// 
		uint8_t* s = bam_get_seq( alignment ) ;

		for( int i=start; i<=end; ++i ) {
			//
			int b = bam_seqi( s, i ) ;
			switch(b) {
			case 1:
				rval << "A" ;
				break ;
			case 2:
				rval << "C" ; 
				break ;
			case 4:
				rval << "G" ;
				break ;
			case 8:
				rval << "T" ;
				break ;
			default:
				rval << "N" ;
				break ;
			}
		}
		//
		return rval.str() ;
	}

	int alignment_sequence_quality( const bam1_t* alignment, int start, int end ) {	
		int rval = 0 ;
		uint8_t* q = bam_get_qual( alignment ) ;
		for( int i=start; i<=end; ++i ) {
			rval += (int) q[i] ;
		}
		rval /= end - start + 1 ;
		return rval ;
	}
	
	std::string Sequence( int tid, int pos, const bam1_t* alignment, void* results ) {
		
		// set the default return value to nothing
		std::string rval = "" ;
		int variantqual  = 0 ;

		// get the relative reference position in the read
		std::size_t relpos = 0 ;
		if( alignment->core.pos < pos ) {
			relpos = (std::size_t) pos - alignment->core.pos ;
		}

		// get the cigar 
		Cigar cigar = Cigar( ) ;
		cigar.from_hts_cigar( bam_get_cigar(alignment), alignment->core.n_cigar ) ;
		
		// get the relative position in the CIGAR string and the query position
		std::size_t opbin = 0 ;
		std::size_t oppos = 0 ;
		std::size_t q_pos = cigar.reference_offset_to_query_offset( relpos, opbin, oppos ) ;

		// Figure out the read-sequence
		
		// A reference call or SNP will be annotated 
		// as either an M, X or =
		//
		// An insertion is present in the read but not 
		// in the reference. We will call the insertion 
		// if the opbin + 1 operation equals I and oppos equals 
		// the runLength. In this manner, calling insertions 
		// is normalized with 
		//
		// A deletion is present in the reference but not 
		// in the read. Deletions are called if oppos equals
		// the runLength of the operation and opbin + 1 equals
		// D.
		//
		// With the following rules, InDel variants are called 
		// on the base before they occur.
		//
		// These rules need to applied in the following order, 
		// deletion, insertion, and finally SNPs. 
		//

		// did we call an indel previously
		bool calledindel = false ;

		// check that the current base is a regular base
		if( cigar.at( opbin )->operation() == 'M' || cigar.at( opbin )->operation() == '=' || cigar.at( opbin )->operation() == 'X') {

			// check that we are on the correct position
			if( (oppos) == cigar.at( opbin )->runLength() - 1) {

				// we can't call an indel on the last base of the alignment
				if ( opbin < cigar.length() - 1 ) {

					// check for indels
					if( cigar.at( opbin + 1 )->operation() == 'D' ) {
						//
						rval  = alignment_sequence( alignment, q_pos, q_pos ) ;
						rval += std::string( cigar.at( opbin + 1 )->runLength(), '-' ) ;

						variantqual = alignment_sequence_quality( alignment, q_pos, q_pos ) ;
						//
						calledindel = true ;
					} else if( cigar.at( opbin + 1)->operation() == 'I' ) {
						calledindel = true ;						
						rval = alignment_sequence( alignment, q_pos, q_pos + cigar.at( opbin + 1 )->runLength()  ) ;						

						// get the variant quality
						variantqual = alignment_sequence_quality( alignment, q_pos, q_pos + cigar.at( opbin + 1 )->runLength() ) ;						
					}
				}
			}

			// we did not call any indels, so we are 
			// free to call the reference sequence or 
			// a SNP
			if( ! calledindel ) {
				rval = alignment_sequence( alignment, q_pos, q_pos ) ;
				variantqual = alignment_sequence_quality( alignment, q_pos, q_pos ) ;
			}
		}

		if( results != NULL ) {
			SequenceInformation* ret = (SequenceInformation*) results ;
			ret->query_position      = q_pos ;			
			ret->quality             = variantqual ;
		}

		// print the currently inferred data
		/*
		std::cerr << "tid: " << tid 
			<< ", position: " << pos 
			<< ", sequence: " << rval 
			<< ", indel: " << calledindel
			<< ", alignment start: " << alignment->core.pos
			<< ", relative position: " << relpos 
			<< ", query position: " << q_pos
			<< ", cigar length:"<< alignment->core.n_cigar 
			<< ", is unmapped: " << (alignment->core.flag & BAM_FUNMAP) 
			<< ", cigar-operation: " << cigar.at(opbin)->operation() 
			<< ", length: " << cigar.at(opbin)->runLength()  
			<< ", at " << oppos 
			<< std::endl ;
		*/

		// the sequence
		return rval ;
	}


	std::string ReadGroup( const bam1_t* alignment ) {
		// initialize the return value to unknown
		std::string rval = "unknown" ;

		// get the readgroup tag
		uint8_t* p = bam_aux_get( alignment, "RG" ) ; 
		if( p ) {

			// set the readgroup
			rval = std::string( bam_aux2Z(p) ) ;
		}
		return rval ;
	}


	std::string GetLabel( const bam1_t* alignment, std::string label )  {
		// initialize the return value to unknown
		std::string rval = "unknown" ;

		if( label.size() != 2 )
			return rval ;

		// get the readgroup tag
		uint8_t* p = bam_aux_get( alignment, label.c_str() ) ; 
		if( p ) {

			// set the readgroup
			rval = std::string( bam_aux2Z(p) ) ;
		}
		return rval ;
	}
	
	std::string Strand( const bam1_t* alignment ) {
		return bam_is_rev(alignment) ? "reverse" : "forward" ;
	}
	
	std::string Amplicon( const bam1_t* alignment ) {
		uint8_t* p = bam_aux_get( alignment, "am" ) ; 
		if( p ) {
			return std::string( bam_aux2Z(p) ) ;
		} 
		return "unknown" ;		
	}


}
