#include "nimbusheader.h"
#include "Manager.h"

namespace NimApp {

	using namespace std ;
	using namespace threadutils ;


	void closeAndDelete(  ifstream* s ) {
		if( s != NULL ) {
			if( s->is_open() ) {
				s->close() ;
			}
			delete s ;
			s = NULL ;
		}		
	}

	void closeAndDelete(  ofstream* s ) {
		if( s != NULL ) {
			if( s->is_open() ) {
				s->close() ;
			}
			delete s ;
		}
		s = NULL ;
	}

	Manager::Manager( ) {

		// prepare the kill switch
		_stop   = new Signal<bool>(false) ; 
		_oqueue = new TQueue< Nimbus::AlignmentBuilder* >() ;
		
		// set the input to NULL
		_pfa = NULL ;
		_pfb = NULL ;

		_in  = NULL ;
		_out = NULL ;

		// set the output to NULL
		_pfo = NULL ;

		// make an empty worker vector
		_workers = vector<Worker>() ;
	}

	Manager::~Manager(void) {

		// close the first FastQ file
		if( _pfa != NULL ) {
			if( _pfa->is_open() ) _pfa->close() ;
			delete _pfa ;
		}

		// close the second FastQ file
		if( _pfb != NULL ) {
			if( _pfb->is_open() ) _pfb->close() ;
			delete _pfb ;
		}

		// close the output stream
		if( _pfo != NULL ) {
			if( _pfo->is_open() ) _pfo->close() ;
			delete _pfo ;
		}	

		// delete the signals
		if( _stop != NULL ) delete _stop ;
		if( _oqueue != NULL ) delete _oqueue ;
		if( _in != NULL ) delete _in ;
		if( _out != NULL ) delete _out ;
	}

	//
	// 
	//

	void Manager::addForwardInput( string fn ) {		
		_pfa = new ifstream( fn.c_str(), fstream::in ) ;
	}

	void Manager::addReverseInput( string fn ) {		
		_pfb = new ifstream( fn.c_str(), fstream::in ) ;
	}

	void Manager::addInput( string fna, string fnb ) {
		addForwardInput( fna ) ;
		addReverseInput( fnb ) ;
	}

	void Manager::addOutput( string fn ) {
		_pfo = new ofstream( fn.c_str(), fstream::out ) ;

	}
	
	void Manager::finalizeStreams( ) { 
		_out = new Writer( _pfo, _oqueue, _stop ) ;
		_in  = new Reader( _pfa, _pfb ) ;
	}

	void Manager::writeToOutput( std::string s ) {
		if( _pfo != NULL ) {
			(*_pfo) <<  s ;
		}
	}

	void Manager::addWorkers(  Nimbus::AmpliconAlignment* a, int n ) {
		for( int i=0; i<n; i++ ){
			_workers.push_back( Worker( a, _stop, _in->getQueue(), _oqueue ) ) ;
		}
	}

	//
	//
	//

	void Manager::run( ) {

		// get the threads running
		vector<thread> wthreads = vector<thread>() ;
		for( vector<Worker>::iterator it=_workers.begin(); it!=_workers.end(); ++it ) {
			Worker w = *it ;
			wthreads.push_back( thread( &Worker::run, w) );
		}
		thread in  = thread( &Reader::run, *_in ) ;
		thread out = thread( &Writer::run, *_out ) ; ;

		// wait for the input to have finished
		in.join() ;
		cerr << "[Manager] Input has been processed" << endl ;
		
		// wait for the number of input reads to match the output
		while( _in->getCounter()->get() != _out->getCounter()->get() ) {
			// poll the counter every second
#if __cplusplus >= 201103L				
			this_thread::sleep_for( chrono::seconds(1) ) ;
#elif _POSIX_VERSION >= 200112L
			sleep( 1 ) ;
#endif
		}
		cerr << "[Manager] All reads have been written to the output" << endl ;
		cerr << "[Manager] Sending the stop signal to the workers" << endl ;
		_stop->set( true ) ;
		for( vector<thread>::iterator it=wthreads.begin(); it!=wthreads.end(); ++it ) {				
			it->join() ;				
		}
		cerr << "[Manager] processed " << _in->getCounter()->get() << " elements in the input" << endl ;
		cerr << "[Manager] processed " << _out->getCounter()->get() << " elements in the output" << endl ;

	 	// join the writer thread
		out.join() ;
		cerr << "[Manager] Joined all the threads" << endl ;

		if( _in->getQueue() != NULL ) delete _in->getQueue() ;
		if( _in->getStopSignal() != NULL ) delete _in->getStopSignal() ;
		if( _in->getCounter() != NULL ) delete _in->getCounter() ;
		if( _out->getCounter() != NULL ) delete _out->getCounter() ;
	}
	
}