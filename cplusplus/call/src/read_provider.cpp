
#include <read_provider.h>

// STL
#include <cstdlib>
#include <string>		// std::string
#include <vector>		// std::vector
#include <sstream>		// std::stringstream
#include <iostream>		// std::cerr 

// HTS lib imports
#include <sam.h>

namespace nimbus {

	//
	//
	// Private utility functions
	//
	//

	std::size_t __index_of__( const std::vector<std::string> V, const std::string _ ) {
		std::size_t rval = V.size() ;
		for( std::size_t i=0; i<V.size(); ++i ) {
			if( V[i].compare(_) == 0 ) {
				rval = i ;
				break ;
			}
		}
		return rval ;
	}
	

	std::size_t __index_of__( const std::vector< std::pair<std::string, std::string> > V, 
		const std::pair<std::string, std::string> _ ) {
		std::size_t rval = V.size() ;
		for( std::size_t i=0; i<V.size(); ++i ) {
			if( V[i].first.compare(_.first) == 0 && 
				V[i].second.compare(_.second) == 0) {
				rval = i ;
				break ;
			}
		}
		return rval ;
	}
	
	//
	//
	//
	//
	//

	SequenceProvider::SequenceProvider() {
		
		//
		samfiles = std::vector<samFile*>() ;
		headers  = std::vector<bam_hdr_t*>() ; 
		lengths  = std::vector< std::size_t >() ;
		names    = std::vector< std::string >() ;
		samples  = std::vector< std::pair<std::string, std::string> >() ;

		//
		data  = NULL ;
		p     = NULL ;

		// set the default maximum depth
		maximum_depth = 1000 ;
		minimum_mapping_quality = 20 ;
		
		// initialize the Results
		pileup       = MpileupResult() ;
		pileup.tid   = 0 ;
		pileup.pos   = 0 ;		
		pileup.plp   = NULL ;
		pileup.n_plp = NULL ;

		// 
		iter  = NULL ;

	}

	SequenceProvider::~SequenceProvider() {

		// remove the sam files and their headers
		for( std::size_t i=0; i<samfiles.size(); ++i ) {
			samFile* fp   = samfiles[i] ;
			bam_hdr_t* hd = headers[i] ;
			
			bam_hdr_destroy( hd ) ;
			sam_close( fp ) ;			
		}		
	}

	void SequenceProvider::addSamFile( std::string fn ) {
		samFile* s   = sam_open( fn.c_str(), "rb" ) ;		
		bam_hdr_t* h = sam_hdr_read( s ) ;

		filenames.push_back( fn ) ;
		samfiles.push_back( s ) ;
		headers.push_back( h ) ;
	}

	void SequenceProvider::getInformation( ) {

		// foreach registered samfile
		for( std::size_t i=0; i<samfiles.size(); ++i ) {

			// get the header information
			bam_hdr_t* header = headers[i] ;

			// foreach sequence in the header of the first sam file
			if( i == 0 ) {

				//
				names.resize( (std::size_t) header->n_targets, "?" ) ;
				lengths.resize( (std::size_t) header->n_targets, 0 ) ;

				for( int k=0; k<header->n_targets; ++k ) {
					
					// get the replicon name
					std::string name = std::string( header->target_name[k] ) ;
					names[(std::size_t) k]   = name ;
					lengths[(std::size_t) k] = (std::size_t) header->target_len[k] ;
				}
			}

			// parse the samples from the header
			parseSamples( std::string(header->text) ) ;

		}
	}

	void SequenceProvider::initializePileup( ) {
			
		//
		p = new ProviderOptions() ;
		int minqual = minimum_mapping_quality ;
		minqual = minqual > 255 ? 255 : minqual ;
		minqual = minqual < 0 ? 0 : minqual ;
		p->minimum_mapping_quality = minqual ;

		Provider** data = new Provider*[ samfiles.size() ] ;		//
		for( std::size_t i=0; i<samfiles.size(); ++i ) {
			data[i]          = new Provider() ;
			data[i]->sam     = samfiles[i] ;
			data[i]->header  = headers[i] ;
			data[i]->options = p ;			 
		}

		// Prepare the iterators over the 
		iter = bam_mplp_init( samfiles.size(), ReadProvider, (void**) data ) ;		
		bam_mplp_set_maxcnt( iter, (int) maximum_depth ) ;

		// Allocate memory for the results
		pileup.plp   = new const bam_pileup1_t*[samfiles.size()] ;
		pileup.n_plp = new int[samfiles.size()] ;

	}

	void SequenceProvider::closePileup() {
		
		// close the iterator
		bam_mplp_destroy( iter ) ;

		// remove the pileups
		if( pileup.plp != NULL ) {			
			delete[] pileup.plp ;
		}

		// remove the read depths
		if( pileup.n_plp != NULL ) {
			delete[] pileup.n_plp ;
		}
		
		// remove the read providers
		if( data != NULL ) {
			for( std::size_t i=0; i<samfiles.size(); ++i ){
				delete data[i] ;
			}
			delete[] data ;
		}

		// delete the provider
		if( p != NULL ) {
			delete p ;
		}
	}

	bool SequenceProvider::next( ) {
		bool rval = false ;
		if( iter != NULL ) {			
			int ret = bam_mplp_auto(iter, &pileup.tid, &pileup.pos, pileup.n_plp, pileup.plp ) ;
			rval = ret > 0 ? true : false ;
		}
		return rval ;
	}

	//
	// Protected functions
	//

	void SequenceProvider::parseSamples( std::string text ) {

		// convert the text to a stream 
		//  so can use getline on it
		std::stringstream stream( text ) ;
		std::string line = "" ;

		// parse each line in the text
		while( std::getline(stream, line, '\n').good() ) {

			// if the line represents a header line
			if( line.substr(0,3).compare("@RG") == 0  ) {

				// initialize the variables
				std::size_t b = std::string::npos ;
				std::size_t e = std::string::npos ;
				std::string sample    = "" ;
				std::string readgroup = "" ;

				// parse the sample information
				b = line.find( "SM:" ) ;
				if( b == std::string::npos )
					continue ;			// when no samples are found, move over to the next line				

				b += 3 ;
				e = line.find( "\t", b ) ;
				sample = e != std::string::npos ? line.substr( b, e - b ) : line.substr( b ) ;
				
				// get the readgroup information
				b = line.find( "ID:" ) ;
				if( b == std::string::npos ) 
					continue ;			// when no readgroups are found, move over to the next line
				
				b += 3 ;
				e = line.find( "\t", b ) ;
				readgroup = e != std::string::npos ? line.substr( b, e - b ) : line.substr( b ) ;

				// check whether the pair is already present
				std::pair<std::string, std::string> p = std::pair<std::string, std::string>(sample, readgroup) ;
		
				// add the read-group sample pair if it was not already present
				std::size_t smpidx = __index_of__( samples, p ) ;
				if( smpidx == samples.size() ) 
					samples.push_back( p ) ;
			}
		}

	}

	std::size_t SequenceProvider::n_entries() const {
		return samfiles.size() ;
	}

	std::vector<std::string> SequenceProvider::getFileNames(){
		return filenames ;
	}

	
	int ReadProvider( void* data, bam1_t *b ) {
		// std::cerr << "Called ReadProvider" << std::endl ;
		int rval ;

		// cast the options
		Provider* opt = (Provider*) data ;
		
		// iterate over the reads in the BAM file, until we 
		// have a read that we can process
		while( true ) {

			// get a read from the BAM file
			rval = sam_read1( opt->sam, opt->header, b ) ;
			
			// if we could not get a read, end the loop
			if( rval < 0 ) 
				break ;

			if( b->core.tid < 0 || (b->core.flag & BAM_FUNMAP) ) {
				rval = -1 ;
				continue ;
			}

			// test whether the alignment score of the read is sufficient
			if( b->core.qual < opt->options->minimum_mapping_quality ) {
				rval = -1 ;
				continue ;
			}

			break ;
		}
		// std::cerr << "Finished ReadProvider" << std::endl ;
		return rval ;
	}


}
