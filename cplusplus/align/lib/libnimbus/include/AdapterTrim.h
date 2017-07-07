#pragma once

#include "stdafx.h"
#include "Read.h"

namespace Nimbus {

	namespace utils {

		
		//
		// A class to trim sequences from reads 
		//
		class AdapterTrim {
			
			std::string _seq ;
			std::string _seed ;

		public:	
			AdapterTrim( std::string sequence, int seedsize );
	
			~AdapterTrim(void);

			/**
			 These functions trim the sequences in read/seq using a seeded and recursive method
			 **/
			basic::Read* trim( basic::Read* read ) ;
			std::string trim( std::string seq, std::string pre ) ;
			std::string trim( std::string seq ) ;
		} ;


	}
}