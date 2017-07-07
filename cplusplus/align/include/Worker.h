#pragma once

#include "nimbusheader.h"

namespace NimApp {

		class Worker {
			// the stop signal
			threadutils::Signal<bool>* _stop ;

			// the input Queue
			threadutils::TQueue< std::pair<Nimbus::basic::Read*, Nimbus::basic::Read*> >* _in ;

			// the output Queue
			threadutils::TQueue< Nimbus::AlignmentBuilder* >* _out ;

			unsigned int _limit ;

			Nimbus::AmpliconAlignment* _aa ;

		public:
			Worker( Nimbus::AmpliconAlignment* a, threadutils::Signal<bool>* s, threadutils::TQueue< std::pair<Nimbus::basic::Read*, Nimbus::basic::Read*> >* i, threadutils::TQueue< Nimbus::AlignmentBuilder* >* o, unsigned int l )  ;

			Worker( Nimbus::AmpliconAlignment* a, threadutils::Signal<bool>* s, threadutils::TQueue< std::pair<Nimbus::basic::Read*, Nimbus::basic::Read*> >* i, threadutils::TQueue< Nimbus::AlignmentBuilder* >* o ) ;

			Worker( Nimbus::AmpliconAlignment* a, threadutils::TQueue< std::pair<Nimbus::basic::Read*, Nimbus::basic::Read*> >* i, threadutils::TQueue< Nimbus::AlignmentBuilder* >* o ) ;

			~Worker(void) ;

			threadutils::Signal<bool>* getStopSignal() ;

			/*
		 	 process the alignments in a paired end manner
			 */ 
			Nimbus::AlignmentBuilder* process( std::pair<Nimbus::basic::Read*,Nimbus::basic::Read*> p ) ;


			/**
			 Runs the processing loop 
			 **/
			void run() ;

		};

	
}