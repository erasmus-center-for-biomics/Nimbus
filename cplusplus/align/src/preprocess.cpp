
#include "nimbusheader.h"

#include "opt.h"

using namespace std ;
using namespace Nimbus::basic ;
using namespace Nimbus::IO ;
using namespace Nimbus::utils ;

/*
 * Preprocess the data 
 *
 */
void preprocess( string in, string out, string adapter, int seedsize, bool truncname ) {
	
	// the input and output streams
	ifstream* _hi    = new ifstream( in.c_str(), fstream::in ) ;
	ofstream* _ho    = new ofstream( out.c_str(), fstream::out ) ;	 
	AdapterTrim* _at = new AdapterTrim( adapter, seedsize ) ;

	// the read
	Read* r = FastQReader( *_hi ) ;
	while( r != NULL ) {

		// trim the adapter sequence from the read
		Read* pt = _at->trim( r ) ; 
		
		// do we need to truncate the name?
		if( truncname ) {

			// find the first space character
			Read* tmp = NULL ;
			size_t idx = pt->name().find( " " ) ;

			// if we found a spce
			if( idx != std::string::npos && idx > 1 ) {

				// get the new name of the read which is a substring to the first space
				std::string name = pt->name().substr( 0, idx - 1  ) ;

				// copy and switch the read around
				tmp = pt ;
				pt = new Read( name, tmp->sequence(), tmp->quality() ) ;

				// delete the original
				delete tmp ;
			}
		}

		// print the read to the output stream
		(*_ho) << pt->fastq() << endl ; 

		// clean up 
		delete r ;
		delete pt ;

		// get the next read from the fastq file
		r = FastQReader( *_hi ) ;
	}
	
	// close the input and output handles
	_hi->close() ;
	_ho->close() ;

	//
	delete _at ;
	delete _hi ;
	delete _ho ;
}

//
using namespace commandline ;

/**
 * the main for trimming algorithm
 */
int preprocess_main( int argc, char* argv[] ) {

	std::string fn_in   = "" ;
	std::string fn_out  = "default.fq" ;
	std::string adapter = "AGATCGGAAGAGC" ;
	int seedsize        = 1 ;
	bool truncname      = false ;

	// prepare a new option parser
	OptParser* op = new OptParser( "nimbus trim" ) ;

	//
	// add the required options
	op->add( 'i', "input", true, true, "a FastQ with reads to trim" ) ;
	op->add( 'o', "output", true, true, "the FastQ file where the trimmed reads will be written to" ) ;

	// add the optionals 
	op->add( 'a', "adapter", false, true, "the sequence of the adapter to trim" ) ;
	op->add( 's', "seedsize", false, true, "the minimum length of the adapter to match" ) ;
	op->add( 't', "truncate-names", false, false, "should the readnames be truncated?" ) ;

		// parse the provided options
	op->interpret( argc, argv ) ;

	// check whether all required options are set
	vector<string> miss = op->missingOptions() ;
	if( miss.size() > 0 ) {
		string mess = "Missing options:" ;
		for( unsigned int i=0; i<miss.size(); ++i ) {
			if( i > 0 ) 
				mess += ", " ;
			mess += miss[i] ;
		}

		// quit due to missing arguments
		op->usageInformation( mess, true ) ;
	}

	// check input file presence
	if( ! FileExists(op->getValue( "input")) ) {
		op->usageInformation( "Input FastQ file " +  op->getValue( "input") + " not found", true ) ;
	} else {
		fn_in = op->getValue( "input") ;
	}
	
	// set the output file
	if( op->getValue( "output") != "" ) 
		fn_out = op->getValue( "output") ;

	// set the adapter sequence of required
	if( op->getValue( "adapter") != "" ) 
		adapter = op->getValue( "adapter") ;

	// set the seed size 
	if( op->getValue( "seedsize") != "" ) 
		seedsize = atoi( op->getValue( "seedsize").c_str() ) ;
	
	if( op->getValue( "truncate-names") != "" ) 
		truncname = true ;

	// run the preprocess procedure
	preprocess( fn_in, fn_out, adapter, seedsize, truncname ) ;


	// remove the option parser
	delete op ;

	//
	return 0 ;
}
