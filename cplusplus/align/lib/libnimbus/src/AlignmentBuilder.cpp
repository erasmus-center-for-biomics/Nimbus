#include "stdafx.h"
#include "AlignmentBuilder.h"
#include "Utils.h"

namespace Nimbus {

	using namespace std ;
	using namespace basic ;
	using namespace alignment ;

	//
	// alignment set
	//

	AlnSet::AlnSet() {
		amplicon    = NULL ;
		f_alignment = NULL ;
		r_alignment = NULL ;
		f_path      = NULL ;
		r_path      = NULL ;
		f_record    = NULL ;
		r_record    = NULL ;
	}

	AlnSet::AlnSet( Amplicon* a ) {
		amplicon    = a ;
		f_alignment = NULL ;
		r_alignment = NULL ;
		f_path      = NULL ;
		r_path      = NULL ;
		f_record    = NULL ;
		r_record    = NULL ;
	}

	AlnSet::~AlnSet() {		

	}

	void AlnSet::delete_content() {		
		if( f_alignment != NULL ) delete f_alignment ;
		if( r_alignment != NULL ) delete r_alignment ;
		if( f_path      != NULL ) delete f_path ;
		if( r_path      != NULL ) delete r_path ;
		if( f_record    != NULL ) delete f_record ;
		if( r_record    != NULL ) delete r_record ;
	}

	void AlnSet::align( AlignmentScore* scores, int gapopen, Read* f ) {
		if( amplicon != NULL ) {
			if( amplicon->forward() ) {
				f_alignment = new SmithWaterman( scores, gapopen, amplicon->sequence(), f->sequence() ) ;
			} else {
				f_alignment = new SmithWaterman( scores, gapopen, amplicon->sequence(), f->rc_sequence() ) ;
			}
		}
	}

	void AlnSet::align( AlignmentScore* scores, int gapopen, Read* f, Read* r ) {
		if( amplicon != NULL ) {
			align( scores, gapopen, f ) ;

			// align the second read
			if( r != NULL && amplicon->forward() ) {
				r_alignment = new SmithWaterman( scores, gapopen, amplicon->sequence(), r->rc_sequence() ) ;
			} else if( r != NULL && !amplicon->forward() ) {
				r_alignment = new SmithWaterman( scores, gapopen, amplicon->sequence(), r->sequence() ) ;
			}
		}
	}

	void AlnSet::SAMrecord( Read* f, Read* r ) {	
		if( amplicon != NULL && f_alignment != NULL && r_alignment != NULL ) {
			SAMrecord_f( f ) ;
			SAMrecord_r( r ) ;
			f_record->mate( *r_record ) ;
			r_record->mate( *f_record ) ;
		}
	}	

	void AlnSet::SAMrecord_f( Read* f ) {
		if( f != NULL && amplicon != NULL && f_alignment != NULL ) {

			// get the position
			f_path  = f_alignment->getPath() ;
			int idx = utils::QueryStart( *f_path ) ;
			int pos = f_path->at(idx).first + amplicon->start() ;

			// set the SAM records
			f_record = new SAMRecord( *f, false, amplicon->chromosome(), pos, amplicon->forward() ) ;
			f_record->cigar( f_alignment->QCIGAR() ) ;
			f_record->mapq( 150 ) ;
			f_record->add_tag( "am", 'Z', amplicon->format() ) ;
			f_record->add_tag( "AS", 'i', f_alignment->getAlignmentScore() ) ;
			f_record->setProperlyAligned() ;
			f_record->setFirstSegmentInTemplate() ;

			// add the nr tag using the Levenshtein distance
			int start     = f_path->at(idx).first - 1 ;
			int len       = f_path->at( utils::QueryEnd(*f_path) ).first  - start + 1 ;
			string rseq   = amplicon->sequence().substr( start, len ) ;
			Levenshtein l = Levenshtein( rseq, f_record->sequence() ) ;
			f_record->add_tag( "NM", 'i', l.distance() ) ;

			// add the reference sequence tag (rs)
			// f_record->add_tag( "rs", 'Z', rseq ) ;
		}
	}

	void AlnSet::SAMrecord_r( Read* r ) {
		if( r != NULL && amplicon != NULL && r_alignment != NULL ) {
			// get the position
			r_path  = r_alignment->getPath() ;
			int idx = utils::QueryStart( *r_path ) ;
			int pos = r_path->at(idx).first + amplicon->start() ;
			
			// set the SAM records
			r_record = new SAMRecord( *r, false, amplicon->chromosome(), pos, !amplicon->forward() ) ;
			r_record->cigar( r_alignment->QCIGAR( *r_path ) ) ;
			r_record->mapq( 150 ) ;
			r_record->add_tag( "am", 'Z', amplicon->format() ) ;
			r_record->add_tag( "AS", 'i', r_alignment->getAlignmentScore() ) ;
			r_record->setProperlyAligned() ;
			r_record->setLastSegmentInTemplate() ;

			// add the nr tag using the Levenshtein distance	
			int start     = r_path->at(idx).first - 1 ;
			int len       = r_path->at( utils::QueryEnd(*r_path) ).first - start + 1  ;
			string rseq   = amplicon->sequence().substr( start, len ) ;				
			Levenshtein l = Levenshtein( rseq, r_record->sequence() ) ;
			r_record->add_tag( "NM", 'i', l.distance() ) ;

			// add the reference sequence tag (rs)
			// r_record->add_tag( "rs", 'Z', rseq ) ;
		}
	}

	//
	// ResultSet functions
	//
	
	AlignmentBuilder::AlignmentBuilder() {
		forward = NULL ;
		reverse = NULL ;
		entries = vector<AlnSet>() ;
	}

	AlignmentBuilder::AlignmentBuilder( Read* f ) {
		forward = f ;
		reverse = NULL ;
		entries = vector<AlnSet>() ;
	}
	
	AlignmentBuilder::AlignmentBuilder( Read* f, Read* r ) {
		forward = f ;
		reverse = r ;
		entries = vector<AlnSet>() ;
	}
	
	AlignmentBuilder::~AlignmentBuilder(void) {
		
	}

	bool AlignmentBuilder::empty() const {
		return entries.size() == 0  ;
	}

	void AlignmentBuilder::add( Amplicon* a ) {
		entries.push_back( AlnSet(a) ) ;
	}

	void AlignmentBuilder::align( AlignmentScore* scores, int gapopen ) {
		for( vector<AlnSet>::iterator it=entries.begin(); it!=entries.end(); ++it ) {
			it->align( scores, gapopen, forward, reverse ) ;
		}
	}


	int AlignmentBuilder::best() {
		int rval = -1 ;
		vector<int> scores = vector<int>() ;
		scores.resize( entries.size() ) ;
		for( unsigned int i=0; i<entries.size(); i++ ) {

			// only consider initialized alignments
			if( entries[i].f_alignment != NULL )  {

				// set rval to the first initialized entry, if not defined 
				if( rval == -1 ) rval = i ; 

				// record the scores
				if(entries[i].r_alignment != NULL)  {
					scores[i] = entries[i].r_alignment->getAlignmentScore() + entries[i].f_alignment->getAlignmentScore() ;
				} else {
					scores[i] = entries[i].f_alignment->getAlignmentScore() ;
				}

				if( scores[i] > scores[rval] ) rval = (int) i ;
			}
		}
		return rval ;
	}

	void AlignmentBuilder::createRecord( AlnSet& a ) {		
		if( forward != NULL && reverse != NULL ) {
			a.SAMrecord( forward, reverse ) ;
		} else if( forward != NULL ) {
			a.SAMrecord_f( forward ) ;
		} else if( reverse != NULL ) {
			a.SAMrecord_r( reverse ) ;
		}
		if( a.f_record != NULL ) a.f_record->add_tag( "nh", 'i', (int) entries.size() ) ;
		if( a.r_record != NULL ) a.r_record->add_tag( "nh", 'i', (int) entries.size() ) ;
	}

	void AlignmentBuilder::createRecords( ) { 
		for( vector<AlnSet>::iterator it=entries.begin(); it!=entries.end(); ++it ) {
			createRecord( *it ) ;
		}
	}

	bool AlignmentBuilder::samrecordspresent( ) const {
		bool rval = false ;
		for( vector<AlnSet>::const_iterator it=entries.cbegin(); it!=entries.cend(); ++it ) {
			if( it->f_record != NULL || it->r_record != NULL ) rval = true ;
		}
		return rval ;
	}

}
