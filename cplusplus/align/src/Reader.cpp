#include "nimbusheader.h"
#include "Reader.h"


using namespace std ;
using namespace Nimbus ;
using namespace Nimbus::basic ;
using namespace Nimbus::IO ;
using namespace threadutils ;

namespace NimApp {

	Reader::Reader(  ) {	
		_ha = NULL ;
		_hb = NULL ;

		_limit  = LIMIT ;  // keep a maximum of 10,000 read pairs in the queue
		_queue  = new TQueue< pair<Read*,Read*> >() ;		
		_sigcnt = new Signal<long>( 0 ) ;		
		_stop   = new Signal<bool>( false ) ;
	}

	
	Reader::Reader( istream* xa, istream* xb ) {	
		_ha = xa ;
		_hb = xb ;

		_limit  = LIMIT ;  // keep a maximum of 10,000 read pairs in the queue
		_queue  = new TQueue< pair<Read*,Read*> >() ;
		_sigcnt = new Signal<long>( 0 ) ;
		_stop   = new Signal<bool>( false ) ;
	}

	Reader::Reader( istream* xa, istream* xb, unsigned int l ) {
		_ha = xa ;
		_hb = xb ;

		_limit  = l ;  // keep a maximum of 10,000 read pairs in the queue
		_queue  = new TQueue< pair<Read*,Read*> >() ;
		_sigcnt = new Signal<long>( 0 ) ;
		_stop   = new Signal<bool>( false ) ;
	}

	Reader::Reader( istream* xa, istream* xb, unsigned int l, TQueue<pair<Read*,Read*>>* q, Signal<long>* s ) {
		_ha = xa ;
		_hb = xb ;

		_limit  = l ;  // keep a maximum of 10,000 read pairs in the queue
		_queue  = q ;
		_sigcnt = s ;
		_stop   = new Signal<bool>( false ) ;
	}

	Reader::Reader( istream* xa, istream* xb, unsigned int l, TQueue<pair<Read*,Read*>>* q, Signal<long>* s, Signal<bool>* b ) {
		_ha = xa ;
		_hb = xb ;

		_limit  = l ;  // keep a maximum of 10,000 read pairs in the queue
		_queue  = q ;
		_sigcnt = s ;
		_stop   = b ;
	}

	Reader::~Reader() {
	}



	pair<Read*,Read*> Reader::process( bool& proceed ) {

		// declare the return value
		pair<Read*,Read*> rval ;

		Read* a = NULL ;
		Read* b = NULL ;
		
		if( _ha != NULL ) a = FastQReader( *_ha ) ; 
		if( _hb != NULL ) b = FastQReader( *_hb ) ;

		// if these are ok, add them to the output
		if( a != NULL ) {
			rval = pair<Read*,Read*>( a, b ) ;			
			proceed = true ;
		} else {
			// if not, set the proceed bool to false
			proceed = false ;
		}
		
		// return the value
		return rval ;
	}

	//
	// The processing loop
	//
	void Reader::run() {
		// while we are allowed		
		bool proceed = true ;

		//long inputcnt = 0 ;

		while( proceed ) {
				
			// give the other processes some breathing room if the queue is full		
			while( _queue->size() >= _limit ) {
#if __cplusplus >= 201103L				
				this_thread::sleep_for( chrono::milliseconds(THREADSLEEPTIME) ) ;
#elif _POSIX_VERSION >= 200112L
				usleep( THREADSLEEPTIME * 1000 ) ;
#endif		
			}
			
			// push the results in the queue				
			pair<Read*,Read*> x = process( proceed ) ;			
			
			if( proceed ) { 

				// report the input count
				//if( inputcnt % 1000000 == 0 )
				//	cerr << "[InputReader] Processed " <<  inputcnt << " read (pairs)" << endl ;

				// add a read counter
				//inputcnt += 1 ;

				// add a new input to the stream
				_queue->push( x ) ;
				
				// update the counter signal
				long c = _sigcnt->get() ;
				c += 1 ;
				_sigcnt->set( c ) ;
			}

			// stop if the we get the signal
			if( _stop->get() ) proceed = false ;
		}
	}


}