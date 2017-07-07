#include "nimbusheader.h"
#include "Worker.h"

namespace NimApp {

	using namespace std ;
	using namespace threadutils ;
	using namespace Nimbus ;
	using namespace Nimbus::basic ;

	Worker::Worker( AmpliconAlignment* a, Signal<bool>* s, TQueue< pair<Read*, Read*> >* i, TQueue< AlignmentBuilder* >* o, unsigned int l )  {
		_aa    = a ;
		_stop  = s ;
		_in    = i ;
		_out   = o ;
		_limit = l ;
	}

	Worker::Worker( AmpliconAlignment* a, Signal<bool>* s, TQueue< pair<Read*, Read*> >* i, TQueue< AlignmentBuilder* >* o ) {
		_aa    = a ;
		_stop  = s ;
		_in    = i ;
		_out   = o ;
		_limit = LIMIT ;
	}

	Worker::Worker( AmpliconAlignment* a, TQueue< pair<Read*, Read*> >* i, TQueue< AlignmentBuilder* >* o ) {
		_aa    = a ;
		_stop  = new Signal<bool>( false ) ;
		_in    = i ;
		_out   = o ;
		_limit = LIMIT ;				
	}

	Worker::~Worker() {
	}

	Signal<bool>* Worker::getStopSignal() {
		return _stop ;
	}

	/*
	 	* process the alignments in a paired end manner
		*/ 
	AlignmentBuilder* Worker::process( pair<Read*,Read*> p ) {				
		AlignmentBuilder t = _aa->align( p ) ;
		return new AlignmentBuilder( t ) ;
	}


	/*
		* Runs the processing loop 
		*/
	void Worker::run() {

		// keep running untill the stop signal
		bool proceed = true ;
		while( proceed ) {

			// check whether we should stop processing
			if( _stop->get() ) break ;

			// wait with processing while the output is full			
			while( _out->size() > _limit ) {				
#if __cplusplus >= 201103L				
				this_thread::sleep_for( chrono::milliseconds(THREADSLEEPTIME) ) ;
#elif _POSIX_VERSION >= 200112L
				usleep( THREADSLEEPTIME * 1000 ) ;
#endif
				if( _stop->get() ) {
					proceed = false ;
					break ; 
				}				
			}			

			// sleep if there is no input
			//while( _in->empty() && ! _stop->get() ) {				
				// this_thread::sleep_for( chrono::milliseconds(1) ) ;
			//}

			// get a new read from the stack
			pair<Read*,Read*> p ;
			bool ok = _in->shift( p ) ;
			if( ok ) {
						
				// add the result to the output queue
				_out->push( process( p ) ) ;
			}
		}
	}

}