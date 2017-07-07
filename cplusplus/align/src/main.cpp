
#include "nimbusheader.h"
#include "nimbus.h"
#include "preprocess.h"
// #include "block.h"

using namespace std ;

/**
 * the usage for the main loop
 *
 */
void main_usage( string message, bool terminate ) {

	// print the error message
	printf( "%s\n\n", message.c_str() ) ;
	printf( "Usage:\n" ) ;
	printf( " nimbus [function] [options]\n" ) ;
	printf( "\n" ) ;
	printf( "Functions:\n" ) ;
	printf( "  trim\ttrims adapter sequences from the reads in a FastQ file\n" ) ;
	printf( "  align\taligns the provided reads to the amplicons\n" ) ;
//	printf( "  count\tcounts the aligned reads per amplicon\n" ) ;
	printf( "\n" ) ;

	// exit if the error code exceeds -1 
	if( terminate ) {
		exit( EXIT_FAILURE ) ;
	}
}

/**
 * The main entry point into this application
 * 
 */
int main(int argc, char* argv[]) {

	// print the arguments when no arguments provided
	if( argc < 2 ) 
		main_usage( "No arguments provided", true ) ;
	
	// check whether the specified function is present
	vector<string> func = vector<string>() ;
	func.push_back( "trim" ) ;
	func.push_back( "align" ) ;
	// func.push_back( "count" ) ;
	int fnum = -1 ;
	for( unsigned int i=0; i<func.size(); i++ ) {
		if( string(argv[1]) == func[i] ) 
			fnum = i ;
	}

	if( fnum == -1 ) 
		main_usage( "Function '" + string(argv[1]) + "' is not recognized", true ) ;

	// run the main for the specified function 
	if( fnum == 0 ) {		
		preprocess_main( argc, argv ) ;
	} else if( fnum == 1 )  {
		nimbus_main( argc, argv ) ;
	} else {
		main_usage( "function out of bounds", true ) ;
	}

	// always return 0
	return 0 ;
}
