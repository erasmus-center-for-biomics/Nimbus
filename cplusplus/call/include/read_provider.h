#pragma once 

// STL imports 
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>	// std::pair

// HTS lib imports
#include <sam.h>

namespace nimbus {


	typedef struct __provider_options__ {
		uint8_t minimum_mapping_quality ;
	} ProviderOptions ;	

	typedef struct __provider__ {
		samFile* sam ;
		bam_hdr_t* header ; 
		ProviderOptions* options ;
	} Provider ;

	/**
	 * A struct to hold the mpileup results
	 *
	 */
	typedef struct __mpileup_result__ {
		int tid ;
		int pos ;
		const bam_pileup1_t **plp ;
		int *n_plp ;		
	} MpileupResult ;

	class SequenceProvider {

		// the input files
		std::vector<samFile*> samfiles ;
		std::vector<bam_hdr_t*> headers ;
		std::vector<std::string> filenames ;

		// the provider for the data
		Provider** data ; 
		ProviderOptions* p ;

		// the mpileup iterator
		bam_mplp_t iter ;

	public:
		
		// the pileup
		MpileupResult pileup ;

		// the maximum depth
		std::size_t maximum_depth ; 
		std::size_t minimum_mapping_quality ; 

		// information obtained from the sam headers
		std::vector< std::size_t > lengths ;
		std::vector< std::string > names ;
			
		// a vector to hold the sample vs readgroup mapping
		std::vector< std::pair<std::string, std::string> > samples ;

	public:

		/**
		 * Construct a new SequenceProvider object
		 */
		SequenceProvider() ;

		/**
		 * Destroys a SequenceProvider object
		 */
		~SequenceProvider() ;

		/**
		 * Add a new alignment file to the sequence provider object
		 *
		 * @param fn - the input BAM file
		 *
		 * @returns: success or failure
		 */
		void addSamFile( std::string fn ) ;

		/**
		 * Retrieves the header information from each of the 
		 *  BAM files and fills the relevant fields		 
		 */
		void getInformation() ;

		/**
		 * Converts the stored information in such 
		 *  a manner that it is useable by bam_mplp_init
		 *
		 * @param po - provider options that are constant 
		 *			   for all the alignment files
		 *
		 * @returns: void** data
		 */
		void initializePileup( ) ;

		/**
		 * closes the mpileup 
		 *
		 */
		void closePileup( ) ;

		/**
		 *
		 *
		 */
		bool next() ;

		/**
		 *
		 *
		 */
		std::size_t n_entries() const ;

		std::vector<std::string> getFileNames() ;

	protected:

		/**
		 * Parses the sample information from the SAM file header 
		 *  and registers its contents in the samples vector. 
		 *
		 * @param text - the text of the header to parse 
		 */
		void parseSamples( std::string text ) ;

	} ;


	/**
	 * Gets reads from a BAM file
	 *
	 */
	int ReadProvider( void* data, bam1_t *b ) ;

}