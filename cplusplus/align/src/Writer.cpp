#include "nimbusheader.h"
#include "Writer.h"

namespace NimApp {

		using namespace std ;
		using namespace threadutils ;
		using namespace Nimbus ;
		using namespace Nimbus::basic ;
		using namespace Nimbus::alignment ;

		Writer::Writer( ) {
			_in     = NULL ;
			_out    = NULL ;
			_sigcnt = new Signal<long>( 0 ) ;
			_stop   = new Signal<bool>( false ) ;		
		}

		Writer::Writer( ostream* o ) {
			_in     = NULL ;
			_out    = o ;
			_sigcnt = new Signal<long>( 0 ) ;
			_stop   = new Signal<bool>( false ) ;
		}

		Writer::Writer( ostream* o, TQueue<AlignmentBuilder*>* q, Signal<bool>* b ) {
			_in     = q ;
			_out    = o ;
			_sigcnt = new Signal<long>( 0 ) ;
			_stop   = b ;
		}

		Writer::Writer( ostream* o, TQueue<AlignmentBuilder*>* q, Signal<long>* s, Signal<bool>* b ) {
			_in     = q ;
			_out    = o ;
			_sigcnt = s ;
			_stop   = b ;
		}

		Writer::~Writer() {			
		}

		//
		Signal<bool>* Writer::getStopSignal() {
			return _stop ;
		}

		Signal<long>* Writer::getCounter() {
			return _sigcnt ;
		}

		bool Writer::process( AlignmentBuilder* value ) { 

			// if there is no opened output stream: stop the iteration
			if( _out == NULL ) return false ;

			// should never happen, but if it does don't kill writer
			if( value == NULL ) return true ;

			// write all the generated SAM records
			for( vector<AlnSet>::iterator it=value->entries.begin();  it!=value->entries.end(); ++it ) {
				if( it->f_record != NULL ) (*_out) << it->f_record->str() << endl ; 
				if( it->r_record != NULL ) (*_out) << it->r_record->str() << endl ; 
			}
			
			// if we did not have any SAM entries, write 
			// empty samrecords
			if( ! value->samrecordspresent()  ) {

				// create an empty samrecord
				SAMRecord* f = NULL ;
				SAMRecord* r = NULL ;

				// 
				if( value->forward != NULL ) {
					f = new SAMRecord( *(value->forward) ) ;
					f->setPaired() ;
					f->setFirstSegmentInTemplate() ;
				}				
				if( value->reverse != NULL ) {
					r = new SAMRecord( *(value->reverse) ) ;
					r->setPaired() ;
					r->setLastSegmentInTemplate() ;
				}

				// write and delete samrecords
				if( f != NULL && r != NULL ) {
					f->mate( *r ) ;
					r->mate( *f ) ;
					(*_out) << f->str() << endl ; 
					(*_out) << r->str() << endl ; 
					
				} else if( f != NULL ) {
					(*_out) << f->str() << endl ; 					
				} else if( r != NULL ) {
					(*_out) << r->str() << endl ; 					
				}
				// cleanup 
				if( f != NULL ) delete f ;
				if( r != NULL ) delete r ;
			}


			// clean up the reads			
			if( value->forward != NULL ) delete value->forward ;
			if( value->reverse != NULL ) delete value->reverse ;

			// the alnsets will go out of scope in the run 
			// function and will be cleaned automatically
			for(vector<AlnSet>::iterator it=value->entries.begin();  it!=value->entries.end(); ++it) {
				it->delete_content() ;
			}
			delete value ;

			// returns that we should proceed if we made it this far
			return true ;		
		}

		// run the processor
		void Writer::run( ) {

			// check if we should be active
			bool proceed = true ;

			if( _in == NULL ) proceed = false ;

			// while we need to be active
			while( proceed ) {
		
				// wait for input 
#if __cplusplus >= 201103L
				while( _in->empty() ) {
					this_thread::sleep_for( chrono::milliseconds(1) ) ;
					if( _stop->get() ) {
						proceed = false ;
						break ; 
					}
				}
#endif
				// get a new value from the queue
				AlignmentBuilder* value = NULL ;
				bool ok = _in->shift( value ) ;
				if( ok ) {

					// process the value
					proceed = process( value ) ; 
					
					// update the counter signal
					long c = _sigcnt->get() ;
					c += 1 ;
					_sigcnt->set( c ) ;
				}

				// stop if the we get the signal
				if( _stop->get() ) proceed = false ;
				
			} // end of while loop
		}


}