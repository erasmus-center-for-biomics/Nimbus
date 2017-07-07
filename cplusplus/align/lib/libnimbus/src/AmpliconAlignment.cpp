 #include "stdafx.h"
#include "AmpliconAlignment.h"
#include "Amplicon.h"

namespace Nimbus {


	//
	//
	//
	
	double lambdaCalculator( int match, int mismatch, double precision, double nf ) {
		// double se   = nf * match + (3 * nf) * mismatch ;
		double rval = 1.0 ;
		double high = 2.0 ;
		double low  = 0.0 ;

		while( high - low > precision ) {
			double sum = nf * nf * exp( rval * match ) * 4 + nf * nf * exp( rval * mismatch ) * 12 ;

			if( sum > 1.0 ) {
				high = rval ;
				rval = (rval + low) / 2 ;
			} else {
				low  = rval ;
				rval = (rval + high) / 2 ;
			}
		}
		
		//
		return rval ;
	}

	MappingQuality::MappingQuality( long dbs, double l ) {
		_dbsize = dbs ;		
		_lambda = l ;
		_k      = 0.1 ;
		// printf("lambda: %f, dbsize: %d, K: %f\n", _lambda, _dbsize, _k ) ;
	}

	double MappingQuality::evalue( double bscore, int n ) const {
		// the not normalized e-value calculation
		double rval = _dbsize * n * pow( 2, -1.0 * bscore ) ;
		return rval ;
	}

	double MappingQuality::bscore( int score ) const {
		return ( _lambda * (double)score - log(_k) ) / log(2) ;
	}

	int MappingQuality::phredEncode( double v ) const {
		return (int) ( -10.0 * log10( v ) ) ;
	}

	//
	//
	//

	bool seedPositionFilter( AlnSet a, int limit ) {
		
		//
		bool rval = false ;
		bool f_ok = true ;
		bool r_ok = true ;

		// if there was a forward alignment
		if( a.f_record != NULL && a.amplicon != NULL ) {
			f_ok  = false ;
			int s = a.f_record->pos() ;
			int e = a.f_record->pos() + a.f_record->rlen() ;

			// printf( "%d\t%d\t%d\n", abs(a.amplicon->start() - s), abs(a.amplicon->end() - e), limit ) ;
			if( abs(a.amplicon->start() - s) <= limit ) f_ok = true ;			
			if( abs(a.amplicon->end() - e) <= limit ) 	f_ok = true ;
			if( abs(a.amplicon->start() - e) <= limit ) f_ok = true ;			
			if( abs(a.amplicon->end() - s) <= limit ) 	f_ok = true ;
		}

		// if we had a reverse alignment
		if( a.r_record != NULL && a.amplicon != NULL ) {
			r_ok  = false ;
			int s = a.r_record->pos() ;
			int e = a.r_record->pos() + a.r_record->rlen() ;

			if( abs(a.amplicon->start() - s) <= limit ) r_ok = true ;			
			if( abs(a.amplicon->end() - e) <= limit ) 	r_ok = true ;
			if( abs(a.amplicon->start() - e) <= limit ) r_ok = true ;			
			if( abs(a.amplicon->end() - s) <= limit ) 	r_ok = true ;
		}

		if( !f_ok || !r_ok ) {
			rval = true ;
		}

		// return the OK (or not)
		return rval ;
	}

	//
	//
	//

	AmpliconAlignment::AmpliconAlignment( seed::AmpliconIndex* ai, alignment::AlignmentScore* scores, int seedpos, int go ) {
		_ai      = ai ;
		_scores  = scores ;
		_posd    = seedpos ;
		_gapopen = go ;
		_reportsecondary = false ;
		// mapping quality score calculator
		_mapqual = new MappingQuality( 
					_ai->dbsize(), 
					lambdaCalculator(_scores->_match, _scores->_mismatch, 0.001, 0.25) 
					) ;
	}

	AmpliconAlignment::AmpliconAlignment( seed::AmpliconIndex* ai, alignment::AlignmentScore* scores, int seedpos, int go, bool rs ) {
		_ai      = ai ;
		_scores  = scores ;
		_posd    = seedpos ;
		_gapopen = go ;
		_reportsecondary = rs ;
		// mapping quality score calculator
		_mapqual = new MappingQuality( 
					_ai->dbsize(), 
					lambdaCalculator(_scores->_match, _scores->_mismatch, 0.001, 0.25) 
					) ;
	}

	AmpliconAlignment::~AmpliconAlignment(void) {
		delete _mapqual ;
	}
		

	AlignmentBuilder AmpliconAlignment::align( std::pair<basic::Read*,basic::Read*> p ) const {

		// declare the output variable
		AlignmentBuilder rval = AlignmentBuilder( p.first, p.second ) ;

		// get the amplicons
		std::vector<basic::Amplicon*> ampset = _ai->getAmplicons( p ) ;

		// add the amplicons
		if((int)ampset.size() < _scores->_maxamp){
			for( std::vector<basic::Amplicon*>::iterator it=ampset.begin(); it!=ampset.end(); ++it ) {
				rval.add( *it ) ;
			}
		}

		// align the reads to the amplicon
		rval.align( _scores, _gapopen ) ; 

		// create the relevant samrecords
		if( _reportsecondary ) {			
			rval.createRecords( ) ;
		} else {
			// only create a record for the best record
			int idx = rval.best() ;
			if(idx != -1 ) {
				rval.createRecord( rval.entries[idx] ) ;
			} 
		}
		
		// calculate the mapping scores via the BLAST like algorithm
		for( std::vector<AlnSet>::iterator it=rval.entries.begin(); it!=rval.entries.end(); ++it ) {
			
			if( it->f_record != NULL && it->f_alignment != NULL ) {
				int s = it->f_alignment->getAlignmentScore()  ;
				int mq = _mapqual->phredEncode( 
					_mapqual->evalue( 
						_mapqual->bscore( s ), 
						it->f_record->qlen() 
					) ) ; 
				it->f_record->mapq( mq ) ;				
			}
			if( it->r_record != NULL && it->r_alignment != NULL ) {
				int s = it->r_alignment->getAlignmentScore()  ;
				int mq = _mapqual->phredEncode( 
					_mapqual->evalue( 
						_mapqual->bscore( s ), 
						it->r_record->qlen() 
					) ) ; 
				it->r_record->mapq( mq ) ;				
			}
		}

		//printf( "%d\n", _posd ) ;

		// check the positions relative to the amplicon end
		if( _posd >= 0 ) {
			for( std::vector<AlnSet>::iterator it=rval.entries.begin(); it!=rval.entries.end(); ++it ) {

				// if we would write this record to a SAM file
				if( it->f_record != NULL || it->r_record != NULL ) {

					// check whether the alignment spans the first x bases of the amplicon
					bool flt = seedPositionFilter( *it, _posd ) ;

					// if not unmap the reads and add a tag indicating why
					if( flt ) {
						if( it->f_record != NULL ) {
							std::stringstream s ;
							s << it->f_record->rname() << ";" << it->f_record->pos() << ";" << it->f_record->cigar() ;
							it->f_record->Unmap() ;						
							it->f_record->setFirstSegmentInTemplate() ;
							it->f_record->add_tag("sp", 'Z', s.str() ) ;
						}
						if( it->r_record != NULL ) { 
							std::stringstream s ;
							s << it->r_record->rname() << ";" << it->r_record->pos() << ";" << it->r_record->cigar() ;
							it->r_record->Unmap() ;
							it->r_record->setLastSegmentInTemplate() ;
							it->r_record->add_tag("sp", 'Z', s.str() ) ;
						}
						if( it->f_record != NULL && it->r_record != NULL ) {
							it->f_record->mate( *it->r_record ) ;
							it->r_record->mate( *it->f_record ) ;
						}
					}
				}
			}
		}

		// returns the results from the alignment
		return rval ;
	}


}