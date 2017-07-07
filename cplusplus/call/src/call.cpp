
// to implement
#include <call.h>

// STL library
#include <cstdlib>
#include <string>
#include <sstream>
#include <limits>
#include <iostream>
#include <algorithm>
#include <tuple>		// std::tuple
#include <vector>		// std::vector

// HTSlib
#include <sam.h>
#include <faidx.h>

// own headers
#include <read_provider.h>
#include <sample.h>
#include <aspects.h>
#include <alignment_functions.h>
#include <refsequence.h>

namespace nimbus {

	//
	std::string UNKNOWN_SAMPLE = "UNKNOWN"  ;


	class CompareFields {
		std::size_t idx ; 
	public:
		CompareFields( std::size_t i ){ 
			idx = i ;
		}

		bool operator()( const aspect& a, const aspect& b ) {
			// return a.first[idx].compare(b.first[idx]) < 0 ? true : false ;
			return a.first[idx] < b.first[idx] ;
		}
	} ;


	//
	// Constructors
	//

	CallingVariants::CallingVariants() {

		// set empty info fields 
		infofields = std::vector<std::string>() ;

		// set the accessible objects
		genome   = new GenomeSequence() ;
		provider = new SequenceProvider() ;		

		// calling options
		options = CallingOptions() ;
		options.minimum_mapping_quality   = 0 ;
		options.maximum_reads_in_pileup   = 10000 ;
		options.maximum_alleles           = 4096 ;
		options.minimum_allelic_depth     = 0 ;
		options.minimum_allelic_quality   = 0 ;
		options.minimum_allelic_depth_f   = 0 ;
		options.minimum_allelic_quality_f = 0 ; 		
		options.non_reference             = false ;
		
		//
		AB = std::vector< aspect >() ;
		SI = std::vector<SequenceInformation>() ;
		alleles = std::vector< allele_obj >( ) ;
	}

	/**
	 *
	 *
	 *
	 */
	CallingVariants::~CallingVariants() {
		if( genome != NULL)
			delete genome ;
		if( provider != NULL)
			delete provider ;
	}


	//
	// Run functions
	//

	void CallingVariants::run( std::ostream& out ) {

		//
		std::size_t n_columns = 3 + infofields.size() ;
				
		provider->maximum_depth = options.maximum_reads_in_pileup ;		
		provider->minimum_mapping_quality = options.minimum_mapping_quality ;	

		// initialize the PileUp engine
		provider->initializePileup() ;
		
		//
		initialize_vectors( n_columns ) ;
		
		// write the header
		writeHeader(out) ;

		// Foreach mpileup generated from the BAM file
		while( provider->next() ) {
			
			// Proces the alignments 
			std::size_t totaldepth = processAlignments() ;
			
			// Sort the aspects			
			//std::sort( AB.begin(), AB.begin() + totaldepth, lessthan ) ;			
			for( std::size_t idx=n_columns; idx>0; --idx) {				
				CompareFields cmp( idx-1 ) ;
				std::stable_sort( AB.begin(), AB.begin() + totaldepth, cmp ) ;
			}
			
			// Get the genome sequence at the current position			
			std::string gseq = genome->get( provider->names[provider->pileup.tid], provider->pileup.pos, provider->pileup.pos ) ;			
			if( gseq.empty() ) 
				continue ;

			// Aggregate alleles			
			std::size_t totalvar = aggregateAlignments( totaldepth ) ;
			
			// Report alleles
			bool report = shouldReport( std::string( 1, toupper(gseq[0])), totalvar ) ;

			// Report a variant or not
			if( report ) 
				reportAlleles( out, gseq, totalvar ) ;
		}

		// close the mpileup engine
		provider->closePileup() ;
	}


	//
	//
	//
	//
	//

	void CallingVariants::writeHeader( std::ostream& out ) {
		out << "## " << std::endl ;

		// write the options
		out << "## Options:"  << std::endl ;
		out << "##    maximum-number-of-reads-in-pileup: " <<  options.maximum_reads_in_pileup << std::endl ;
		out << "##    minimum-mapping-quality: " <<  options.minimum_mapping_quality << std::endl ;
		out << "##    maximum-number-of-alleles: " <<  options.maximum_alleles << std::endl ;
		out << "##    minimum-allelic-depth: " <<  options.minimum_allelic_depth << std::endl ;
		out << "##    minimum-allelic-quality: " <<  options.minimum_allelic_quality << std::endl ;
		out << "##    minimum-allelic-depth-frequency: " <<  options.minimum_allelic_depth_f << std::endl ;
		out << "##    minimum-allelic-quality-frequency: " <<  options.minimum_allelic_quality_f << std::endl ;
		out << "##    report-only-non-reference-alleles: " <<  options.non_reference << std::endl ;
		out << "## " << std::endl ;

		// write the info field
		out << "## Info: sample, strand" ;
		for( std::size_t i=0; i<infofields.size(); ++i ) {
			out << ", " << infofields[i]  ;
		}
		out << std::endl ;
		out << "## " << std::endl ;

		// write the FastA file
		out << "## FastA file: " << genome->filename << std::endl ;
		out << "## " << std::endl ;

		// write the BAM files
		std::vector<std::string> filenames = provider->getFileNames() ;

		out << "## BAM files:"  << std::endl ;
		for( std::size_t i=0; i<filenames.size(); ++i ) {
			out << "## " <<  filenames[i] << std::endl ;  ;
		}
	}

	void CallingVariants::initialize_vectors( std::size_t n_columns ) {

		// prepare the aspect list and set it to hold maximumDepth elements
		aspect da = aspect() ;
		da.first.resize( n_columns, "__NOT_INITIALIZED__" ) ;
		AB = std::vector< aspect >( options.maximum_reads_in_pileup, da ) ;
		
		// to store sequence info per alignment
		SI = std::vector<SequenceInformation>( options.maximum_reads_in_pileup ) ;

		// set the default allele
		allele_obj dv = allele_obj() ; 
		dv.sequence = "" ;
		dv.n    = 0 ;
		dv.qual = 0 ;
		dv.info = std::vector<std::string>( n_columns, "__NOT_INITIALIZED__" ) ;
				
		// Initialize the variant vector
		alleles = std::vector< allele_obj >( options.maximum_alleles, dv ) ;								
	}

	std::size_t CallingVariants::processAlignments( ) {
		
		// get the total number of reads overlapping the current position
		std::size_t totaldepth = 0 ;			
			
		// prepare the aspects per read
		std::size_t asit = 0 ; 			
			
		// for each read provider
		for( std::size_t i=0; i<provider->n_entries(); ++i) {

			// iterate over the alignments
			for( std::size_t k=0; k<(std::size_t) provider->pileup.n_plp[i]; ++k ) { 	

				// if the number of alignments is over the maximum 
				//  depth continue without processing them
				if( asit >= options.maximum_reads_in_pileup ) 
					continue ;

				// get the current alignment
				const bam_pileup1_t* p  = provider->pileup.plp[i] + k ;
				const bam1_t* alignment = p->b ;

				// get the sequence 
				std::string sequence = Sequence( provider->pileup.tid, provider->pileup.pos, alignment, &SI[asit] ) ;	// get the sequence													
				std::string strand   = Strand( alignment ) ;						// get the strand 				
				std::string sample   = "unknown" ;									// get the samplename
				std::string readgrp  = ReadGroup( alignment ) ;						// get the readgroup				

				// parse the sample informations
				if( readgrp.compare("unknown") != 0 ) {
					for( std::size_t i=0; i<provider->samples.size(); ++i ) {
						if( provider->samples[i].second.compare(readgrp) == 0 ) {
							sample = provider->samples[i].first ;
							break ;
						}
					}
				}

				// set the quality to the read mapping if this is smaller than the calling quality
				SI[asit].quality = SI[asit].quality < (int) alignment->core.qual ? SI[asit].quality : (int) alignment->core.qual ;
					
				// fill the aspects buffer
				AB[asit].first[0] = sequence ;
				AB[asit].first[1] = sample ;
				AB[asit].first[2] = strand ;

				// set the optional info fields
				for( std::size_t j=0; j<infofields.size(); ++j ) {
					std::string val = GetLabel( alignment, infofields[j] ) ;
					AB[asit].first[ j + 3 ] = val ;
				}

				AB[asit].second   = asit ;

				// increase the aspect iterator
				asit += 1 ;	
			}
		}		

		// save the total number of reads we obtained from the input providers
		totaldepth = asit ;			

		// 
		return totaldepth ;
	}

	std::size_t CallingVariants::aggregateAlignments( std::size_t totaldepth ) {
		std::size_t vi         = 0 ;
		alleles[vi].sequence  = AB[0].first[0] ;
		alleles[vi].n         = 1 ;
		alleles[vi].qual      = SI[ AB[0].second ].quality ;
		for( std::size_t k=0; k<AB[0].first.size(); ++k ) {
			alleles[vi].info[k] = AB[0].first[k] ;
		}

		//
		std::size_t ncol = AB[0].first.size();

		if( totaldepth > 1 ) {				
			for( std::size_t i=1; i<totaldepth; ++i ) {	

				// 
				bool same = true ;
				for( std::size_t idx=ncol; idx>0; --idx ) {
					if( AB[i-1].first[idx-1] != AB[i].first[idx-1] ) {
						same = false ;
						break ;
					}
				}

				// check whether the current read equals the previous read
				//if( compare( AB[i-1], AB[i] ) != 0 ) {
				if( ! same ) { 
					++vi ;

					if( vi >= options.maximum_alleles ) 
						break ;

					// overwrite the current aspects 
					alleles[vi].sequence  = AB[i].first[ 0 ] ;
					alleles[vi].n         = 0 ;
					alleles[vi].qual      = 0 ;
					for( std::size_t k=0; k<AB[i].first.size(); ++k ) {
						alleles[vi].info[k] = AB[i].first[k] ;
					}
				} 

				// update the statistics
				alleles[vi].n    += 1 ;
				alleles[vi].qual += SI[ AB[i].second ].quality ;
			}
		}
		//
		std::size_t totalvar = vi + 1 ;
		return totalvar ;
	}


	bool CallingVariants::shouldReport( std::string gseq, std::size_t totalvar ) {
		
		//	should we report the position	
		bool rval = false ;

		// 
		std::size_t n = 0 ;
		std::size_t q = 0 ;

		// get the totals
		for( std::size_t i=0; i<totalvar; ++i ) {
			n += alleles[i].n ;
			q += alleles[i].qual ;				
		}

		// check each allele
		for( std::size_t i=0; i<totalvar; ++i ) {

			//
			bool report = true ;
				
			// check whether we should not report the allele
			if( alleles[i].n < options.minimum_allelic_depth )
				report = false ;
			if( alleles[i].qual < options.minimum_allelic_quality )
				report = false ;

			double f_n = (double) alleles[i].n / n ;
			if( f_n < options.minimum_allelic_depth_f )
				report = false ;

			double f_q = (double) alleles[i].qual / q ;
			if( f_q < options.minimum_allelic_quality_f )
				report = false ;
			if( options.non_reference && gseq.compare(alleles[i].sequence) == 0 )
				report = false ;

			// 
			if( report )
				rval = true ;
		}
		
		// 
		return rval ;
	}

	void CallingVariants::reportAlleles( std::ostream& out, std::string gseq, std::size_t totalvar ) {
		
		// 
		std::size_t n = 0 ;
		std::size_t q = 0 ;

		// get the totals
		for( std::size_t i=0; i<totalvar; ++i ) {
			n += alleles[i].n ;
			q += alleles[i].qual ;				
		}


		// report the output stream
		out	<< provider->names[ (std::size_t) provider->pileup.tid] 
			<< "\t" << provider->pileup.pos
			<< "\t" << gseq 
			<< "\t" << totalvar 
			<< "\t" << n 
			<< "\t" << q 
			<< std::endl ;
			
		// for each variant report it
		for( std::size_t i=0; i<totalvar; ++i ){

			// aggregate the info fields
			std::stringstream sinfo ;
			for( std::size_t k=1; k<alleles[i].info.size(); ++k ) {
				if(k != 1) 
					sinfo << ", " ;
				sinfo << alleles[i].info[k] ;
			}

			// print some output
			out << "\t" << alleles[i].sequence
				<< "\t" << alleles[i].n
				<< "\t" << alleles[i].qual
				<< "\t" << sinfo.str()
				<< std::endl ;
		}


	}

}