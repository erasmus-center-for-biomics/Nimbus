#pragma once


#include "nimbusheader.h"

namespace NimApp {

	class Reader {
	protected:
		//
		threadutils::Signal<bool>* _stop ; 
		threadutils::Signal<long>* _sigcnt ; 
		threadutils::TQueue<std::pair<Nimbus::basic::Read*,Nimbus::basic::Read*>>*    _queue ;

		unsigned int _limit ;

		// the input streams for the first 
		// and second data reads
		std::istream* _ha ;
		std::istream* _hb ;

	public:

		//
		// Constructors
		//

		Reader(  ) ;
		
		Reader( std::istream* xa, std::istream* xb ) ;
	
		Reader( std::istream* xa, std::istream* xb, unsigned int l ) ;
	
		Reader( std::istream* xa, std::istream* xb, unsigned int l, threadutils::TQueue<std::pair<Nimbus::basic::Read*,Nimbus::basic::Read*>>* q, threadutils::Signal<long>* s ) ;
	
		Reader( std::istream* xa, std::istream* xb, unsigned int l, threadutils::TQueue<std::pair<Nimbus::basic::Read*,Nimbus::basic::Read*>>* q, threadutils::Signal<long>* s, threadutils::Signal<bool>* b ) ;
	
		~Reader() ;

		//
		// accessors
		//

		threadutils::TQueue< std::pair<Nimbus::basic::Read*,Nimbus::basic::Read*> >* getQueue( ) {
			return _queue ;
		}

		threadutils::Signal<long>* getCounter( ) {
			return _sigcnt ;
		} 

		threadutils::Signal<bool>* getStopSignal( ) {
			return _stop ;
		}

		//
		// the processing function
		//
		std::pair<Nimbus::basic::Read*,Nimbus::basic::Read*> process( bool& proceed ) ;

		//
		// The processing loop
		//
		void run() ;


	} ;

	
}