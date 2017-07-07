#pragma once

#include "Amplicon.h"
#include "_DNANode.h"
#include "Read.h"

namespace Nimbus {

	namespace seed {

	

		//
		// AmpliconIndex
		//

		class AmpliconIndex	{

			//
			_DNANode<basic::Amplicon*>* _idx_f_a ;
			_DNANode<basic::Amplicon*>* _idx_r_a ;
			_DNANode<basic::Amplicon*>* _idx_f_b ;
			_DNANode<basic::Amplicon*>* _idx_r_b ;
			std::vector< basic::Amplicon* > _amplicons ;
			int _keysize ;

		public:
			AmpliconIndex(void) ;
			~AmpliconIndex(void) ;

			/**
			 * the size of the database: the 
			 *  total number of nucleotides covered in the 
			 *  amplicon vector.
			 **/
			long dbsize() const ;

			/**
			 * Return the number of amplicons currently loaded
			 **/
			unsigned int n_amplicons() const ;

			/**
			 * Adds an amplicon to the Index 
			 **/
			bool add( basic::Amplicon* a ) ;

			/**
			 * Builds the index
			 **/
			unsigned int build( int keysize ) ;

			/** 
			 * Gets the amplicons corresponding to f and r
			 **/
			std::vector<basic::Amplicon*> getAmplicons( std::string f, std::string r ) const ; 

			std::vector<basic::Amplicon*> getAmpliconsF( std::string f ) const ;

			std::vector<basic::Amplicon*> getAmpliconsR( std::string r ) const ;

			/** 
			 * Gets the amplicons corresponding to f and r
			 **/
			std::vector<basic::Amplicon*> getAmplicons( basic::Read* f, basic::Read* r ) const ; 

			std::vector<basic::Amplicon*> getAmplicons( std::pair<basic::Read*, basic::Read*> p ) const ; 

		protected:

		} ;

	}
}
