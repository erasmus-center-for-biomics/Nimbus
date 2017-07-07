#include "stdafx.h"
#include "SAMrecord.h"

namespace Nimbus {

	namespace alignment { 

		using namespace std ;

		//
		// SAMFlag functions
		//

		SAMFlag::SAMFlag() { 
			_flag = 0 ;
		}

		SAMFlag::SAMFlag( unsigned int f ) {
			_flag = f ;
		}

		
		// protected setting and unsetting functions
		
		bool SAMFlag::isset( unsigned int f) {
			return (_flag & f) == f ;
		}

		void SAMFlag::set( unsigned int f) {
			//if( !isset(f) ) {
			_flag |= f ;
			//}
		}

		void SAMFlag::unset( unsigned int f) {
			if( isset(f) ) {
				_flag ^= f ;
			}
		}

		// specific field setters
		void SAMFlag::setPaired() { set(0x1) ; }
		void SAMFlag::setUnmapped( ) { set(0x4) ; }
		void SAMFlag::setProperlyAligned( ) { set(0x2) ; }
		void SAMFlag::setMateUnMapped( ) { set(0x8) ; }
		void SAMFlag::setReverseComplemented( ) { set(0x10) ; }
		void SAMFlag::setNexSegmentReverseComplemented( ) { set(0x20) ; }
		void SAMFlag::setFirstSegmentInTemplate( ) { set(0x40) ; }
		void SAMFlag::setLastSegmentInTemplate( ) { set(0x80) ; }
		void SAMFlag::setSecondaryAlignment( ) { set(0x100) ; }
		void SAMFlag::setNotPassingQC( ) { set(0x200) ; }
		void SAMFlag::setOpticalDuplicate( ) { set(0x400) ; }

		// specific field unsetters
		void SAMFlag::unsetPaired() { unset(0x1) ; }
		void SAMFlag::unsetUnmapped( ) { unset(0x4) ; }
		void SAMFlag::unsetProperlyAligned( ) { unset(0x2) ; }
		void SAMFlag::unsetMateUnMapped( ) { unset(0x8) ; }
		void SAMFlag::unsetReverseComplemented( ) { unset(0x10) ; }
		void SAMFlag::unsetNexSegmentReverseComplemented( ) { unset(0x20) ; }
		void SAMFlag::unsetFirstSegmentInTemplate( ) { unset(0x40) ; }
		void SAMFlag::unsetLastSegmentInTemplate( ) { unset(0x80) ; }
		void SAMFlag::unsetSecondaryAlignment( ) { unset(0x100) ; }
		void SAMFlag::unsetNotPassingQC( ) { unset(0x200) ; }
		void SAMFlag::unsetOpticalDuplicate( ) { unset(0x400) ; }

		// specific field checkers
		bool SAMFlag::isPaired() { return isset(0x1) ; }
		bool SAMFlag::isUnmapped( ) { return isset(0x4) ; }
		bool SAMFlag::isProperlyAligned( ) { return isset(0x2) ; }
		bool SAMFlag::isMateUnMapped( ) { return isset(0x8) ; }
		bool SAMFlag::isReverseComplemented( ) { return isset(0x10) ; }
		bool SAMFlag::isNexSegmentReverseComplemented( ) { return isset(0x20) ; }
		bool SAMFlag::isFirstSegmentInTemplate( ) { return isset(0x40) ; }
		bool SAMFlag::isLastSegmentInTemplate( ) { return isset(0x80) ; }
		bool SAMFlag::isSecondaryAlignment( ) { return isset(0x100) ; }
		bool SAMFlag::isNotPassingQC( ) { return isset(0x200) ; }
		bool SAMFlag::isOpticalDuplicate( ) { return isset(0x400) ; }

		unsigned int SAMFlag::flag() const {
			return (unsigned int) int(_flag) ;
		}

		//
		// The SAMCore implementation
		//

		SAMCore::SAMCore() {
			_rname = "*" ;
			_pos   = 0 ;
			_mapq  = 0 ;
			_cigar = vector<char>() ;
			_rlen  = -1 ;
			_qlen  = -1 ;
		}

		void SAMCore::rname( string rn ) {
			_rname = string(rn) ;
		}

		string SAMCore::rname() const {
			return _rname ;
		}

		void SAMCore::cigar( string c ) {
			_cigar = vector<char>( c.begin(), c.end() ) ;
		}

		void SAMCore::cigar( vector<char> c ) {
			_cigar = vector<char>( c.begin(), c.end() ) ;
		}

		string SAMCore::cigar() const {
			if( _cigar.size() > 0 ) {
				return utils::runLengthEncode(_cigar, "", "" ) ;
			} else {
				return "*" ;
			}
		}

		void SAMCore::pos( int p ) {
			_pos = int(p) ;
		}

		int SAMCore::pos( ) const {
			return _pos ;			
		}

		void SAMCore::mapq( int q ) {
			// make sure the value is within the appropriate range
			q = q > 255 ? 255 : q ;
			q = q < 0 ? 0 : q ;
			_mapq = int(q) ;
		}

		int SAMCore::mapq( ) const {			
			return _mapq ;			
		}

		int SAMCore::rlen()  {
			int rval = 0 ;
			if( _rlen == -1 ) {
				for( vector<char>::iterator it=_cigar.begin(); it!=_cigar.end(); ++it ) {
					switch(*it){
					case 'M':
					case '=':
					case 'X':
					case 'D':
					case 'N': 
					case 'H':
						rval++ ;
						break ;					
					} ;
				}
				_rlen = rval ;
			} else {
				rval = _rlen ;
			}
			return rval ;
		} 

		int SAMCore::qlen() {
			int rval = 0 ;
			if( _qlen == -1 ) {
				for( vector<char>::iterator it=_cigar.begin(); it!=_cigar.end(); ++it ) {
					switch(*it){
					case 'M':
					case '=':
					case 'X':
					case 'I':
					case 'S': 				
						rval++ ;
						break ;
					} ;
				}
				_qlen = rval; 
			} else {
				rval = _qlen ;
			}
			return rval ;
		}

		//
		// The SAMRecord implementation
		//

		SAMRecord::~SAMRecord(void) {}

		void SAMRecord::mate( SAMRecord& oth ) {
			setPaired() ;
			oth.setPaired() ;
			//
			if( !isUnmapped() && !oth.isUnmapped() )  {
				if( rname() == oth.rname() ){
					_rnext = "="  ;
					_pnext = oth.pos() ;		
					_tlen  = pos() < oth.pos() ? ( oth.pos() - pos() ) + oth.rlen() : (pos() - oth.pos() + rlen() ) * -1 ;
					_tlen  = abs( _tlen ) ; 
					oth.pnext( pos() ) ;
					oth.rnext( "=" ) ;
					oth.tlen( tlen() ) ; 
				} else {
					_rnext = oth.rname() ;
					_pnext = oth.pos() ;
					_tlen  = 0 ;
					oth.rnext( rname() ) ;
					oth.pnext( pnext() ) ;
					oth.tlen( 0 ) ;
				}
				if( oth.isReverseComplemented() ) setNexSegmentReverseComplemented() ;
				if( isReverseComplemented() ) oth.setNexSegmentReverseComplemented() ;
			} else if( oth.isUnmapped() ) {
				// if the other read is unmapped set the appropriate values
				_rnext = "=" ;
				_pnext = pos() ;
				_tlen  = 0 ; 
				setMateUnMapped( ) ;
				unsetNexSegmentReverseComplemented() ;
				
				oth.pos( pos() ) ;
				oth.rname( "*" ) ;	
				oth.rnext( "=" ) ;
				oth.pnext( pnext() ) ;
				oth.tlen( 0 ) ;
				
			} else if( isUnmapped() ) {
				oth.rnext( "=" ) ;
				oth.pnext( oth.pos() ) ;
				oth.tlen( 0 ) ;				
				oth.setMateUnMapped( ) ;
				oth.unsetNexSegmentReverseComplemented() ;
				
				pos( oth.pos() ) ;
				rname( "*" ) ;		
				rnext( "=" ) ;
				pnext( oth.pnext() ) ;
				tlen( 0 ) ;
			}
		}

		// 
		string SAMRecord::str() const {
			stringstream rval ;
			rval << name() << "\t" 
				<< flag() << "\t" 
				<< rname() << "\t" 
				<< pos() << "\t"
				<< mapq() << "\t"
				<< cigar() << "\t"
				<< rnext() << "\t"
				<< pnext() << "\t"
				<< tlen() << "\t"
				<< sequence() << "\t"
				<< quality() ;

			// add the tags to the format
			for( vector<string>::const_iterator it=_tags.cbegin(); it!=_tags.cend(); ++it) { 
				rval << "\t" << *it ;
			}

			// return the string with the format
			return rval.str() ;
		}

		//
		// getters 
		//
		string SAMRecord::rnext() const { 
			return _rnext ;
		}
		int SAMRecord::pnext() const {
			return _pnext ;
		}
		int SAMRecord::tlen() const {
			return _tlen ;
		}

		//
		// setters
		//
		void SAMRecord::rnext( string r ) { 
			_rnext = r ;
		}
		void SAMRecord::pnext( int p ) {
			_pnext = p ;
		}
		void SAMRecord::tlen( int t ) {
			_tlen = t ;
		}

		//
		//
		//
		vector< string > SAMRecord::tags() const {
			return vector<string>( _tags.begin(), _tags.end() ) ;
		}

		//
		void SAMRecord::empty_tags() {
			_tags.clear() ;
		}
	
		void SAMRecord::Unmap( ) { 
			_flag = 0 ;
			setUnmapped() ;
			/*
			unsetProperlyAligned() ;
			//unsetPaired() ;
			unsetReverseComplemented() ;
			unsetFirstSegmentInTemplate() ;
			unsetLastSegmentInTemplate() ;
			*/

			// if the read was mapped to reverse strand, flip it again
			if( isReverseComplemented() ) { 
				sequence( rc_sequence() ) ;
				quality( r_quality() ) ;				
			}

			set_to_defaults() ;
			_cigar.clear() ;
			mapq(0) ;
			pos(0) ;
			rname("*") ;
		}


		//
		//
		///

		SAMHeader::SAMHeader() {
			_names   = vector<string>() ;
			_lengths = vector<int>() ;
		}

		void SAMHeader::add( string nm, int l ) {
			assert( _names.size() == _lengths.size() ) ;
			_names.push_back( nm ) ;
			_lengths.push_back( l ) ;
		}

		vector<string> SAMHeader::names() const {
			return _names ;
		}

		vector<int> SAMHeader::lengths() const {
			return _lengths ;
		}

		string SAMHeader::str() const {
			assert( _names.size() == _lengths.size() ) ;

			stringstream s ;
			s << "@HD\tVN:1.3\tSO:unsorted" << endl ;
			for( unsigned int i=0; i<_names.size(); i++ ) {
				s << "@SQ\tSN:" << _names[i] << "\tLN:" << _lengths[i] << endl ; 
 			}
			return s.str() ;
		}
	}

}