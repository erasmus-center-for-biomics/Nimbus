#include "stdafx.h"
#include "Amplicon.h"

namespace Nimbus {

	namespace basic {
	
		//
		//
		// GenomicRegion
		//
		//

		std::string GenomicRegion::chromosome() const { return _chr ; }

		/*
		 Get the start coordinate of the amplicon
		 */ 
		int GenomicRegion::start() const { return _start ; }

		/*
		 Get the end coordinate of the amplicon
		 */
		int GenomicRegion::end() const { return _end ; }

		/*
		 Return the strand of the amplicon (forward or reverse)
		 */
		bool GenomicRegion::forward() const { return _forward ; }

		/*
		 Get the name of the amplicon
		 */ 
		std::string GenomicRegion::name() const { return _name ; }

		/*
		 Gets the difference between the start and the end of the amplicon
		 */
		int GenomicRegion::width() const { return end() - start() ; }


		std::string GenomicRegion::str( bool nm ) const {
			char strand = '+' ;
			if( !forward() ) strand = '-' ;

			std::stringstream rval ;
			if( name() != "" && nm == true ) {
				rval << chromosome() << ":" << start() << "-" << end() << "(" << strand << ")" << ":" << name() ;
			} else {
				rval << chromosome() << ":" << start() << "-" << end() << "(" << strand << ")"  ;
			}

			return rval.str() ;
		}

		std::string GenomicRegion::str( ) const {
			return str( false ) ;
		}

		std::string GenomicRegion::format() const {
			return str( false ) ;
		}

		bool GenomicRegion::operator<( const GenomicRegion& a ) const {
			bool rval = false ;
			if( chromosome() < a.chromosome() ){ 
				rval = true ;
			} else if( chromosome() == a.chromosome() ) {
				if( start() < a.start() ){ 
					rval = true ;
				} else if( start() == a.start() ) {
					if( end() < a.end() ) {
						rval = true ;
					} else if( end() == a.end() ) {
						if( forward() < a.forward() ) {
							rval = true ;
						}
					}
				}
			}
			return rval ;
		}

		bool GenomicRegion::operator>( const GenomicRegion& a ) const {
			bool rval = false ;
			if( chromosome() > a.chromosome() ){ 
				rval = true ;
			} else if( chromosome() == a.chromosome() ) {
				if( start() > a.start() ){ 
					rval = true ;
				} else if( start() == a.start() ) {
					if( end() > a.end() ) {
							rval = true ;
					} else if( end() == a.end() ) {
						if( forward() > a.forward() ) {
							rval = true ;
						}
					}
				}
			}
			return rval ;
		}

		bool GenomicRegion::operator==( const GenomicRegion& a ) const {
			bool rval = false ;
			if( chromosome() == a.chromosome() && start() == a.start() && end() == a.end() && forward() == a.forward() ) rval = true ;
			return rval ;
		}

		//
		//
		//
		// Amplicon
		//
		//
		//

		/*
		 
		 */
		Amplicon::~Amplicon( ) {
		
		}
		
		/* 
		 Get the sequence of the amplicon
		 */
		std::string Amplicon::sequence() const { return _sequence ; }

		

		/*
		Pretty prints the amplicon
		 */ 
		std::string Amplicon::str() const {			
			std::stringstream rval ;
			rval << GenomicRegion::str() << ":" << sequence() ;
			return rval.str() ;
		}

		//
		//
		// Overloaded operators
		//
		//
		

		//
		//
		// pointer comparison functions for the sort
		//
		//

		bool cmp_lt_amplicon_p( Amplicon* a, Amplicon* b ) {
			return (*a) < (*b) ;
		}


		bool cmp_gt_amplicon_p( Amplicon* a, Amplicon* b ) {
			return (*a) > (*b) ;
		}


		bool cmp_eq_amplicon_p( Amplicon* a, Amplicon* b ) {
			return (*a) == (*b) ;
		}


	}

}