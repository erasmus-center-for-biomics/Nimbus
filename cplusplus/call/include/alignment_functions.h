#pragma once 

// Standard libraries
#include <cstdlib>
#include <string>

// HTS lib
#include <sam.h>

/**
 *
 *
 */
namespace nimbus {

	// the function definition for an aspect provider function
	//		this definition allows us to store functions in a 
	//		vector or array
	typedef std::string (*AspectProvider)( int tid, int pos, const bam1_t* alignment, void* options, void* results ) ;


	/**
	 * a sequence result 
	 */
	typedef struct __sequence_information_ {
		std::size_t query_position ;
		int quality ;
	} SequenceInformation ;


	/**
	 * Gets the sequence at the current position for alignment bam1_t
	 *
	 */
	std::string Sequence( int tid, int pos, const bam1_t* alignment, void* results ) ;

	/**
	 * Gets the amplicon for the current read
	 *
	 */
	std::string Amplicon( const bam1_t* alignment ) ;

	std::string Strand( const bam1_t* alignment ) ;

	std::string ReadGroup( const bam1_t* alignment ) ;

	std::string GetLabel( const bam1_t* alignment, std::string label ) ;

}