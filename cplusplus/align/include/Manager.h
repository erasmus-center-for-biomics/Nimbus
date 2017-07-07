#pragma once

#include "nimbusheader.h"
#include "Reader.h"
#include "Worker.h"
#include "Writer.h"

namespace NimApp {
	
	
	class Manager {
		
		std::ifstream* _pfa ; 
		std::ifstream* _pfb ;
		std::ofstream* _pfo ;

		//
		Reader* _in ;
		Writer* _out ;
		std::vector< Worker > _workers ;

		//
		threadutils::Signal<bool>* _stop ;
		threadutils::TQueue< Nimbus::AlignmentBuilder*>* _oqueue ;
		
	public:
		/*
		 * constructors
		 */
		Manager( ) ;

		// Manager( std::string fnout, std::string fna, std::string fnb ) ;

		// Manager( std::string fnout, std::string fna, std::string fnb, int n_workers ) ;

		/*
		 * add the input files
		 */
		void addForwardInput( std::string fn ) ;

		void addReverseInput( std::string fn ) ;

		void addInput( std::string fna, std::string fnb ) ;

		void finalizeStreams() ;

		/*
		 * opens the output file for writing
		 */ 
		void addOutput( std::string fn ) ;

		/*
		 * worker builder
		 */
		void addWorkers(  Nimbus::AmpliconAlignment* a, int n ) ;
		/*
		 * allows for data to be written to the output stream before 
		 * the manager delegates responsibility to the writer.
		 *
		 * This is specifically meant to write the SAM header
		 */
		void writeToOutput( std::string x ) ;

		/*
		 * destructor
		 */
		~Manager(void);

		/*
		 * runs the threads after everything has been setup
		 */ 
		void run( ) ;

	} ;
	
}