// nimbus.cpp : Defines the entry point for the console application.

//
#include "nimbusheader.h"

//
#include "nimbus.h"
#include "opt.h"

//
#include "AmpliconIndex.h"
#include "io.h"
#include "SAMrecord.h"
#include "Amplicon.h" 
#include "Alignment.h"
#include "AmpliconAlignment.h"
#include "Manager.h"

//
using namespace std ;
using namespace Nimbus ;
using namespace Nimbus::alignment ;
using namespace Nimbus::seed ;
using namespace Nimbus::basic ;
using namespace Nimbus::IO ;
using namespace NimApp ;

/**
 * The alignment procedure
 *
 **/
void NimbusAlignment( 
	string fastq_f, 
	string fastq_r, 
	string design, 
	string fasta, 
	string samfile,
	int maxamplicons,
	int keysize, int match, int mismatch, int gapextend, int gapopen, int seedmargin, int threads ) {
	
	cerr << "[Main] Loading index" << endl ;
	
	// load the amplicons and the samheader
	SAMHeader header ;
	vector<Amplicon*> amplicons ;
	AmpliconReader( amplicons, header, design, fasta ) ;

	// make sure to remove duplicate amplicons
	sort( amplicons.begin(), amplicons.end(), cmp_lt_amplicon_p ) ;
	vector<Amplicon*> amp_toadd = vector<Amplicon*>( amplicons.size() );
	vector<Amplicon*>::iterator it = unique_copy( amplicons.begin(), amplicons.end(), amp_toadd.begin(), cmp_eq_amplicon_p ) ;
	amp_toadd.resize( distance( amp_toadd.begin(), it ) );

	// create an amplicon index
	AmpliconIndex* ai = new AmpliconIndex() ;
	for( vector<Amplicon*>::iterator it=amp_toadd.begin(); it!=amp_toadd.end(); ++it ) { 
		ai->add( *it ) ; 				
	}
	ai->build( keysize ) ;

	cerr << "[Main] Loaded " << ai->dbsize() << " bases in " << ai->n_amplicons() << " amplicons" << endl ;
	cerr << "[Main] Preparing alignment" << endl ;

	// create a new score calculator
	AlignmentScore* scores = new AlignmentScore( match, mismatch, gapextend, maxamplicons ) ;
	AmpliconAlignment* aa  = new AmpliconAlignment( ai, scores, seedmargin, gapopen ) ;

	// create the thread manager
	Manager mng = Manager() ;

	// open the input and output
	mng.addForwardInput( fastq_f ) ;
	mng.addReverseInput( fastq_r ) ;
	mng.addOutput( samfile ) ;
	mng.writeToOutput( header.str() ) ;
	mng.writeToOutput( "@PG\tID:nimbus\tPN:nimbus\tVN:beta\n@CO\t\n" ) ;
	mng.finalizeStreams() ;

	// initialize the workers
	mng.addWorkers( aa , threads ) ;
	cerr << "[Main] Performing alignment" << endl ;
	
	// run the tool
	mng.run() ;
	cerr << "[Main] Finished alignment" << endl ;

	// cleanup
	delete ai ;
	delete aa ;
	delete scores ;
	for( vector<Amplicon*>::iterator it=amplicons.begin(); it!=amplicons.end(); ++it ) delete *it ;	
} ;

using namespace commandline ;

int nimbus_main( int argc, char* argv[] ) {

	// define a new option parser
	OptParser* op = new OptParser( "nimbus align" ) ;

	// add the required options
	op->add( '1', "forward", true, true, "the forward read from the sequencing" ) ;
	op->add( '2', "reverse", true, true, "the reverse read from the sequencing" ) ;
	op->add( 'd', "design", true, true, "the BED file with the design" ) ;
	op->add( 'f', "fasta", true, true, "the FastA file with the genome sequence" ) ;
	op->add( 'o', "sam", true, true, "the SAM output file" ) ;

	// add the optionals 
	op->add( 'x', "maximum-amplicons", false, true, "reads that generate more than this number of candidate amplicons are not considered in the alignment (default: 6000)" ) ;
	op->add( 'k', "key-size", false, true, "the key size to use (default: 7)" ) ;
	op->add( 'm', "match", false, true, "the match score (default: 2)" ) ;
	op->add( 'n', "mismatch", false, true, "the mismatch score (default: -1)" ) ;
	op->add( 'e', "gap-extend", false, true, "the gap extend score (default: -1)" ) ;
	op->add( 'g', "gap-open", false, true, "the gap open score (default: -1)" ) ;
	op->add( 's', "seed-margin", false, true, "the seed margin (default: 5)" ) ;
	op->add( 'w', "workers", false, true, "the number of workers (default: 5)" ) ;

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

	// check file presence
	if( ! FileExists(op->getValue( "forward")) ) 
		op->usageInformation( "Forward FastQ file " +  op->getValue( "forward") + " not found", true ) ;
	
	if( ! FileExists(op->getValue( "reverse")) ) 
		op->usageInformation( "Reverse FastQ file " +  op->getValue( "reverse") + " not found", true ) ;

	if( ! FileExists(op->getValue( "design")) ) 
		op->usageInformation( "Design file " +  op->getValue( "design") + " not found", true ) ;

	if( ! FileExists(op->getValue( "fasta")) ) 
		op->usageInformation( "FastA file " +  op->getValue( "fasta") + " not found", true ) ;

	if( op->getValue("sam") == "" ) 
		op->usageInformation( "SAM output file not provided", true ) ;
	

	// set the optional paramters
	int keysize   = 7 ;
	int match     = 2 ;
	int mismatch  = -1 ;
	int gapextend = -1 ;
	int gapopen   = -1 ; 	
	int seedmargin = 5 ;
	int threads   = 5 ;
	int maxamplicons = 6000 ;

	// set the optional data
	if( op->getValue("maximum-amplicons") != "" )
		maxamplicons = atoi( op->getValue("maximum-amplicons").c_str() )  ;

	if( op->getValue("key-size") != "" )
		keysize = atoi( op->getValue("key-size").c_str() )  ;

	if( op->getValue("match") != "" )
		match = atoi( op->getValue("match").c_str() )  ;

	if( op->getValue("mismatch") != "" )
		mismatch = atoi( op->getValue("mismatch").c_str() )  ;

	if( op->getValue("gap-open") != "" )
		gapopen = atoi( op->getValue("gap-open").c_str() )  ;

	if( op->getValue("gap-extend") != "" )
		gapextend = atoi( op->getValue("gap-extend").c_str() )  ;

	if( op->getValue("seed-margin") != "" )
		seedmargin = atoi( op->getValue("seed-margin").c_str() )  ;

	if( op->getValue("workers") != "" )
		threads = atoi( op->getValue("workers").c_str() )  ;

	// report the options
	cerr << "[Align] calling alignment with the following options:" << endl ;
	cerr << "[Align] -1 " << op->getValue( "forward") << endl ;
	cerr << "[Align] -2 " << op->getValue( "reverse" ) << endl ;
	cerr << "[Align] --design " << op->getValue( "design" ) << endl ; 
	cerr << "[Align] --fasta " << op->getValue( "fasta" ) << endl ; 
	cerr << "[Align] --sam " << op->getValue( "sam" ) << endl ;
	cerr << "[Align] --key-size " << keysize << endl ;
	cerr << "[Align] --match " << match << endl ; 
	cerr << "[Align] --mismatch " << mismatch << endl ; 
	cerr << "[Align] --gap-extend " << gapextend << endl ; 
	cerr << "[Align] --gap-open " << gapopen << endl ; 
	cerr << "[Align] --seed-margin " << seedmargin << endl ; 
	cerr << "[Align] --workers " << threads << endl ;
	cerr << "[Align] --maximum-amplicons " << maxamplicons << endl ;


	// call the nimbus function
	NimbusAlignment( 
		op->getValue( "forward"), 
		op->getValue( "reverse" ),
		op->getValue( "design" ), 
		op->getValue( "fasta" ), 
		op->getValue( "sam" ),
		maxamplicons, keysize, match, mismatch, gapextend, gapopen, seedmargin, threads ) ;

	//
	delete op ;

	return 0 ;
}
