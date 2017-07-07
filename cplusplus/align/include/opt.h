#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <iostream>
#include <fstream>

namespace commandline { 
	/**
		* An option parser class. 
		*/
	class OptParser {

		// the tool name 
		std::string toolname ;

		// the options
		std::vector<std::string> options ;

		// the short option names
		std::vector<std::string> shortoptions ;

		// the long option names
		std::vector<std::string> longoptions ;

		// are the options required
		std::vector<bool> reqoptions ;

		// does the option have a value
		std::vector<bool> hasvalues ;

		// the descriptions of the tools
		std::vector<std::string> descriptions ;

		// a map with the values
		std::map<std::string, std::string> values ;

	public:

		/**
			* the constructor
			*/
		OptParser( std::string tn ) {
			toolname     = tn ;
			options      = std::vector<std::string>() ;
			shortoptions = std::vector<std::string>() ;
			longoptions  = std::vector<std::string>() ;
			reqoptions   = std::vector<bool>() ;
			values       = std::map< std::string, std::string >() ; 
		}

		/**
			* Adds an option to the object
			*
			**/ 
		void add( char shortopt, std::string longopt, bool required, bool hasvalue, std::string desc ) {
			std::string so = std::string( "-" ) ;
			so += shortopt ;
			shortoptions.push_back( so ) ;
			longoptions.push_back( "--" + longopt ) ;
			reqoptions.push_back( required ) ;		
			hasvalues.push_back( hasvalue ) ;
			descriptions.push_back( desc ) ;
		}


		/**
			* set the options
			*/
		void set_options( int argv, char* args[] ) { 
			for( int i=0; i<argv; i++ ){
				options.push_back( std::string( args[i] ) );
			}
		}

		/**
		 * 
		 */
		void interpret( int argv, char* args[] ) {
			set_options( argv, args ) ;
			interpret() ;
		}

		/**
			* Interpret the commandline options
			*
			*/
		void interpret( ) {

			//
			for( unsigned int i=0; i<options.size(); i++ ) {

				// check the short options
				for( unsigned int k=0; k<shortoptions.size(); k++ ) {

					// if the option matches
					if( options[i].compare( shortoptions[k] ) == 0 ) {						

						// add the value to the map (use the long options value)
						std::pair<std::string,std::string> p = std::pair<std::string,std::string>( longoptions[k], "true" ) ;
						if( hasvalues[k] && i < options.size() - 1  ) 
							p.second = options[ i + 1 ] ;
						values.insert( p ) ;
					}				
				}

				// check the long options
				for( unsigned int k=0; k<longoptions.size(); k++ ) {

					//
					if( options[i].compare( longoptions[k] ) == 0 ) {						

						// add the value to the map
						std::pair<std::string,std::string> p = std::pair<std::string,std::string>( longoptions[k], "true" ) ;
						if( hasvalues[k] && i < options.size() - 1) 
							p.second = options[ i + 1 ] ;
						values.insert( p ) ;
					}				
				}

			} // foreach option
		}

		/**
			* usage information
			*
			*/
		std::string usageInformation( std::string mess, bool quit ) const {
			std::stringstream s ;

			// print the Message:
			s << "Message:\n" << mess << "\n\n" ;

			// print the options, if loaded to the object, otherwise don't
			if( options.size() > 0 ) {
				s << "Program called with:\n" ; 
				for( unsigned int i=0; i<options.size(); i++ ) {
					s << " " << options[i] ;
				}
				s << "\n" ;
			}

			// add the usage information
			s << "Usage:\n" ;
			s << toolname  << " ";
			for( unsigned int i=0; i<shortoptions.size(); i++ ) {
				s << shortoptions[i] << " " ;
				if( hasvalues[i] ) s << "[VALUE] " ;
			}
			s << "\n" ;
			s << "\n" ;
			s << "Options\n" ;
			for( unsigned int i=0; i<shortoptions.size(); i++ ) {
				s << shortoptions[i] << "/" << longoptions[i] << " " ;
				
				//
				if( reqoptions[i] ) { 
					s << " required " ;
				} else {
					s << " optional " ;
				}

				if( hasvalues[i] ) { 
					s << " [value] " ;
				} else {
					s << " [] " ;
				}

				s << descriptions[i] << "\n" ;
			}

			// quit if required
			if( quit ){ 
				printf( "%s\n", s.str().c_str()  ) ;
				exit( EXIT_FAILURE ) ;
			}
		
			// 
			return s.str() ;
		}

		/**
			* determine which essential options were missing from the analysis
			*
			*/
		std::vector<std::string> missingOptions( ) {
			std::vector<std::string> rval = std::vector<std::string>() ;

			for( unsigned int i=0; i<longoptions.size(); i++ ) {
				
				//
				std::map<std::string, std::string>::iterator it = values.find( longoptions[i] ) ;

				if( reqoptions[i] && it == values.end() )  {
					rval.push_back( longoptions[i] ) ;
				}
			}
			return rval ;
		}

		/**
			* gets a value after interpreting the commandline options
			*
			*/
		std::string getValue( std::string opt ) {
			std::string rval = "" ;
			std::map<std::string, std::string>::iterator it = values.find( "--" + opt ) ;
			if( it != values.end() ) {
				rval = it->second ;
			}
			return rval ;
		}

	} ;


	bool FileExists( const std::string fn ) ; 	
}