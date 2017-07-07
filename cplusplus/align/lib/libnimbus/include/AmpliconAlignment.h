#pragma once

#include "stdafx.h"
#include "Read.h"
#include "AmpliconIndex.h"
#include "Alignment.h"
#include "AlignmentBuilder.h"

namespace Nimbus {

	//
	// functions outside of the object
	//

	/* 
	 determines whether an AlnSet ought to be filtered based on the 
	 position of the aligments with respect to the ends of the amplicon. 

	 the limit parameter indicates the tolerance for this filter.
	 */
	bool seedPositionFilter( AlnSet a, int limit ) ;
	
	double lambdaCalculator( int _match, int _mismatch, double precission, double nf ) ;

	class MappingQuality {
	public:
		long _dbsize ;
		double _lambda ;
		double _k ;
	public:
		MappingQuality( long dbs, double l ) ;
		 
		double evalue( double bscore, int n ) const ;

		double bscore( int score ) const ;

		int phredEncode( double v ) const ;

	} ;

	//
	// the main class to do the alignment
	//
	class AmpliconAlignment {

		seed::AmpliconIndex* _ai ;
		alignment::AlignmentScore* _scores ;

		int _posd ;
		int _gapopen ;

		bool _reportsecondary ;
		MappingQuality* _mapqual ;

	public:
	
		AmpliconAlignment( seed::AmpliconIndex* ai, alignment::AlignmentScore* scores, int seedpos, int go ) ;
	
		AmpliconAlignment( seed::AmpliconIndex* ai, alignment::AlignmentScore* scores, int seedpos, int go, bool rs ) ;

		~AmpliconAlignment(void);

		/**
		 Aligns the reads provided to the amplicon index
		 **/
		AlignmentBuilder align( std::pair<basic::Read*,basic::Read*> p ) const ;


	protected:

	} ;


}

