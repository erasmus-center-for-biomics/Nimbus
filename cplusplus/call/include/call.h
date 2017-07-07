#pragma once 

// STL 
#include <cstdlib> 
#include <string>
#include <vector>
#include <algorithm>

// HTSlib
#include <sam.h>

// own headers
#include <aspects.h>
#include <read_provider.h>
#include <alignment_functions.h>
#include <mpileup.h>
#include <sample.h>
#include <refsequence.h>

// all code goes in the nimbus namespace
namespace nimbus {


	typedef struct __calling_options__ {
		
		std::size_t maximum_reads_in_pileup ;
		int minimum_mapping_quality ;
		std::size_t maximum_alleles ;

		int minimum_allelic_depth ;		
		int minimum_allelic_quality ;
		double minimum_allelic_depth_f ;
		double minimum_allelic_quality_f ;

		bool non_reference ;		
	} CallingOptions ;

	/**
	 * an object that performs the variant calling
	 */
	class CallingVariants {

			
		// reusable variables for run
		std::vector< aspect > AB ;
		std::vector<SequenceInformation> SI ;
		std::vector< allele_obj > alleles ; 	

		//
		Mpileup mp ;
	public:
		// the genome sequence 
		GenomeSequence* genome ;

		// the sequence provider
		SequenceProvider* provider ;		

		// the calling options
		CallingOptions options ;
						
		// The optional info fields on which to split alleles
		std::vector<std::string> infofields ;

	
	public:
		/**
		 * the object constructor
		 * 
		 */
		CallingVariants() ;
	
		/**
		 *
		 *
		 *
		 */
		~CallingVariants() ;

		/**
		 * Runs the variant calling
		 *
		 * @param out - the output stream to which to print the results
		 */
		void run( std::ostream& out ) ;
	
	protected:
				
		/**
		 *
		 *
		 */
		void initialize_vectors( std::size_t n_columns ) ;

		/**
		 *
		 *
		 */
		std::size_t processAlignments( ) ;
		
		/**
		 *
		 *
		 *
		 */
		std::size_t aggregateAlignments( std::size_t totaldepth ) ;

		/**
		 *
		 *
		 */
		bool shouldReport( std::string gseq, std::size_t totalvar ) ;

		/**
		 *
		 *
		 *
		 */
		void reportAlleles( std::ostream& out, std::string gseq, std::size_t totalvar ) ;


		void writeHeader( std::ostream& out ) ;
	} ;



}