
// Standard library
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Boost libraries
#include <boost/program_options.hpp>

// Own headers
#include <call.h>
#include <refsequence.h>

// set the namespace
using namespace boost ;
namespace po = boost::program_options ;
using namespace nimbus ;

/**
 * The entry point of the program
 *
 */
int main( int argc, char** argv ) {

	// Options to obtain from the commandline
	std::string fnout   = "" ;
	std::string fnfasta = "" ;
	std::vector< std::string > bamfiles   = std::vector< std::string >() ;
	std::vector< std::string > infofields = std::vector< std::string >() ;
	std::size_t maximum_reads_in_pileup = 10000 ;
	int minimum_mapping_quality         = 20 ;
	std::size_t maximum_alleles         = 4096 ;
	int minimum_allelic_depth           = 0 ;		
	int minimum_allelic_quality         = 0 ;
	double minimum_allelic_depth_f      = 0.0 ;
	double minimum_allelic_quality_f    = 0.0 ;
	bool non_reference                  = false ;
	

	// Parse the options with Boost program_options module.
	po::options_description desc( "Allowed options" ) ;
	desc.add_options()
		( "help", "Produce the help message")
		("bam,b", po::value< std::vector<std::string> >(&bamfiles), "The input BAM files." )
		("fasta,f", po::value< std::string >(&fnfasta), "The samtools indexed FastA file containing the reference sequence." )
		("output,o", po::value< std::string >(&fnout ), "The output file." )
		("info", po::value< std::vector<std::string> >(&infofields), "additional BAM fields on which divide alleles." )
		("minimum-mapping-quality", po::value< int >(&minimum_mapping_quality), "The minimum mapping quality of a read for it to be considered, default: 20" )
		("maximum-number-of-reads-in-pileup", po::value< std::size_t >(&maximum_reads_in_pileup), "The maximum number of reads to consider per position in the genome, default: 10000" )
		("maximum-number-of-alleles", po::value< std::size_t >(&maximum_alleles), "The maximum number of alleles to report upon, default: 4096" )
		("minimum-allelic-depth", po::value< int >(&minimum_allelic_depth), "The minimum read-depth of an allele to report upon, default: 0" )
		("minimum-allelic-quality", po::value< int >(&minimum_allelic_quality), "The minimum quality of an allele to report upon, default: 0" )
		("minimum-allelic-depth-frequency", po::value< double >(&minimum_allelic_depth_f), "The minimum read-depth frequency of an allele to report upon, default: 0.0" )
		("minimum-allelic-quality-frequency", po::value< double >(&minimum_allelic_quality_f), "The minimum quality frequency of an allele to report upon, default: 0.0" )
		("report-only-non-reference-alleles", po::value< bool >(&non_reference), "Should we report only positions with non-reference alleles, default: false" )
	;

	po::positional_options_description p ;
	po::variables_map vm ;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm) ;
	po::notify(vm) ;

	if( vm.count("help") ) {
		std::cerr << "Usage" << std::endl ; 
		std::cerr << desc << std::endl; 
		return 0 ;
	}
	if( ! vm.count("bam") ) {
		std::cerr << "No BAM files provided." << std::endl ;
		std::cerr << "Usage" << std::endl ; 
		std::cerr << desc << std::endl ;  
		return 0 ;
	}


	// create a variant caller
	CallingVariants caller = CallingVariants() ;	

	// register the FastA file at the caller	
	if( fnfasta.compare("") != 0 ) {
		caller.genome->set( fnfasta ) ;	
	}

	// register the first BAM file
	for( std::size_t i=0; i<bamfiles.size(); ++i ) {
		caller.provider->addSamFile( bamfiles[i] ) ; 	
	}
	
	// set the options for the caller
	caller.options.maximum_reads_in_pileup   = maximum_reads_in_pileup ;
	caller.options.minimum_mapping_quality   = minimum_mapping_quality ;
	caller.options.maximum_alleles           = maximum_alleles ;
	caller.options.minimum_allelic_depth     = minimum_allelic_depth ;		
	caller.options.minimum_allelic_quality   = minimum_allelic_quality ;
	caller.options.minimum_allelic_depth_f   = minimum_allelic_depth_f ;
	caller.options.minimum_allelic_quality_f = minimum_allelic_quality_f ;
	caller.options.non_reference             = non_reference ;

	// add the info fields to differentiate alleles with
	for( std::size_t i=0; i<infofields.size(); ++i ){
		caller.infofields.push_back( infofields[i] ) ;
	}

	// Get the data from the BAM files
	caller.provider->getInformation() ;
	
	// run the variant caller
	if( !fnout.empty() ) {
		std::ofstream out( fnout.c_str(), std::ofstream::out ) ; 
		caller.run( out ) ;
		out.close() ;
	} else {
		caller.run( std::cout ) ;
	}

	// return 0
	return 0 ;
}

