#pragma once

#include "nimbusheader.h"

namespace NimApp {


	//template< class T> 
	class Writer {
	protected:
		threadutils::TQueue<Nimbus::AlignmentBuilder*>* _in ;

		// control signals
		threadutils::Signal<bool>* _stop ;
		threadutils::Signal<long>* _sigcnt ;

		// the output stream
		std::ostream* _out ;
		
	public:

		//
		// constructors
		//
		Writer(  ) ; 

		Writer( std::ostream* o ) ; 

		Writer( std::ostream* o, threadutils::TQueue<Nimbus::AlignmentBuilder*>* q, threadutils::Signal<bool>* b ) ; 

		Writer( std::ostream* o, threadutils::TQueue<Nimbus::AlignmentBuilder*>* q, threadutils::Signal<long>* s, threadutils::Signal<bool>* b ) ; 

		//
		// destructor
		//
		~Writer(void) ;

		//
		threadutils::Signal<bool>* getStopSignal() ;

		threadutils::Signal<long>* getCounter() ;

		//
		// processors
		//

		bool process( Nimbus::AlignmentBuilder* value ) ;

		// run the processor
		void run( ) ;
	} ;		
	


}