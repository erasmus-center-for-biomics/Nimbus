#pragma once
#include "stdafx.h"
#include "Read.h"
#include "Utils.h"


namespace Nimbus {

	namespace alignment { 


		/**
		 Overview
		 ========

		 This source file is to define the classes associated with SAM files. 
		 
		 The central class in this source file is the SAMRecord class.

		 SAMRecords are defined from multiple parent classes to 
		 divide the logic for read information, flag information
		 and positional information. One these classes is external 
		 with the Read class. 

		 These classes make up the SAMRecord via multiple inheritance.
		 
		 SAMFlag
		 -------
		 This class contains all the _flag integer variable and all
		 the named functions used to set the different bits in the flag. 

		 SAMCore
		 -------
		 The main determinants for the alignment are defined here, such
		 as the reference name, position and the CIGAR string. 

		 SAMRecord
		 ---------
		 This class is derived from the Read, SAMFlag and SAMRecord 
		 classes. The SAMRecord class itself also defines the mate-pair 
		 related functions.  

		 SAMHeader
		 ---------
		 A class representing the SAM header.  
		**/

		class SAMFlag { 

			//
			// A class to hold all the FLAG related functions
			//

		protected:
			unsigned int _flag ;
		public:

			SAMFlag( ) ;

			SAMFlag( unsigned int _f ) ;

			// set the bits for the respective properties
			void setPaired() ;
			void setUnmapped( ) ;
			void setProperlyAligned( ) ;
			void setMateUnMapped( ) ;
			void setReverseComplemented( ) ;
			void setNexSegmentReverseComplemented( ) ;
			void setFirstSegmentInTemplate( ) ;
			void setLastSegmentInTemplate( ) ;
			void setSecondaryAlignment( ) ;
			void setNotPassingQC( ) ;
			void setOpticalDuplicate( ) ;

			// unset the bits for the respective properties
			void unsetPaired( ) ;
			void unsetUnmapped( ) ;
			void unsetProperlyAligned( ) ;
			void unsetMateUnMapped( ) ;
			void unsetReverseComplemented( ) ;
			void unsetNexSegmentReverseComplemented( ) ;
			void unsetFirstSegmentInTemplate( ) ;
			void unsetLastSegmentInTemplate( ) ;
			void unsetSecondaryAlignment( ) ;
			void unsetNotPassingQC( ) ;
			void unsetOpticalDuplicate( ) ;

			// check whether a bit is set
			bool isPaired( ) ;
			bool isUnmapped( ) ;
			bool isProperlyAligned( ) ;
			bool isMateUnMapped( ) ;
			bool isReverseComplemented( ) ;
			bool isNexSegmentReverseComplemented( ) ;
			bool isFirstSegmentInTemplate( ) ;
			bool isLastSegmentInTemplate( ) ;
			bool isSecondaryAlignment( ) ;
			bool isNotPassingQC( ) ;
			bool isOpticalDuplicate( ) ;

			unsigned int flag() const ;
		protected:
			void set( unsigned int x ) ;

			bool isset( unsigned int x ) ; 

			void unset( unsigned int x ) ;
		} ;


		

		class SAMCore {
			//
			// the alignment core
			//

		protected:
			std::string _rname ;
			int _pos ;
			int _mapq ;
			std::vector<char> _cigar ;
			int _rlen ;
			int _qlen ;

		public:

			//
			// constructors
			//

			SAMCore() ;

			SAMCore(std::string rn, int p, int m, std::string _cig ): _rname(rn), _pos(p), _mapq(m), 
				_cigar(_cig.begin(), _cig.end() ), _rlen(-1), _qlen(-1) { } 

			SAMCore(std::string rn, int p, int m, std::vector<char> _cig ): _rname(rn), _pos(p), _mapq(m), 
				_cigar(_cig.begin(), _cig.end() ), _rlen(-1), _qlen(-1) { } 
			
			~SAMCore() {}
			
			/**
			 get and set reference name
			 **/
			std::string rname() const ;

			void rname(std::string rn ) ;

			/**
			 get and set position
			 **/
			int pos() const ;

			void pos( int p ) ;

			/**
			 get and set mapping quality
			 **/
			int mapq() const ;

			void mapq( int q ) ;

			/**
			 Get the CIGAR string of the read
			 **/ 
			std::string cigar() const ;

			/**
			 set the cigar field
			 **/
			void cigar( std::string c )  ;

			void cigar( std::vector<char> c )  ;

			/**
			 the length of the record on the reference
			 **/ 
			int rlen() ;

			/**
			 the length of the record on the query
			 **/
			int qlen() ;
		} ;

		//
		// create the SAMRecord 
		//

		class SAMRecord: public basic::Read, public SAMFlag, public SAMCore {
			
			std::string _rnext ;
			int _pnext ;
			int _tlen ;

			std::vector< std::string > _tags ;

		public:

			// constructors
			SAMRecord( basic::Read r ): Read(r), SAMFlag(), SAMCore("*", 0, 0, "") {
				setUnmapped() ;
				set_to_defaults() ;
			}


			SAMRecord( basic::Read r, bool unmapped, std::string rname, int pos, bool forward ): Read(r), SAMFlag(), SAMCore(rname, pos, 0, "") {

				// if the record represents an unmapped read
				if( unmapped  ) setUnmapped() ;

				// set the strand to reverse complement and invert 
				// the sequence if not mapped to the forward strand
				if( !forward ) set_to_reverse_strand() ;

				// set our mate fields to default
				set_to_defaults() ;

			}

			SAMRecord( basic::Read r, bool unmapped, std::string rname, int pos, bool forward, std::string cig ): Read(r), SAMFlag(), SAMCore(rname, pos, 0, cig) {

				// if the record represents an unmapped read
				if( unmapped  ) setUnmapped() ;

				// set the strand to reverse complement and invert 
				// the sequence if not mapped to the forward strand
				if( !forward ) set_to_reverse_strand() ;
				
				// set our mate fields to default
				set_to_defaults() ;
			}

			SAMRecord( basic::Read r, bool unmapped, std ::string rname, int pos, bool forward, int mq, std::string cig ): Read(r), SAMFlag(), SAMCore(rname, pos, mq, cig) {

				// if the record represents an unmapped read
				if( unmapped  ) setUnmapped() ;

				// set the strand to reverse complement and invert 
				// the sequence if not mapped to the forward strand
				if( !forward ) set_to_reverse_strand() ;
				
				// set our mate fields to default
				set_to_defaults() ;
			}

			~SAMRecord(void) ;

			// getters						
			std::string rnext() const ;
			int pnext() const ;
			int tlen() const ;
			std::vector< std::string > tags() const ;
			void empty_tags() ;

			void rnext( std::string r ) ;
			void pnext( int ) ;
			void tlen( int ) ;
			/**
			 Set the mate information for this read

			 **/
			void mate( SAMRecord& oth ) ;

			std::string str() const ;

			/*
			 Sets the read unmapped and resets the fields
			 */
			void Unmap() ;

			/** 
			 add any value as a tag
			 **/
			template<class T>
			void add_tag( std::string nm, char type, T content ) {
				std::stringstream s ;
				s << nm << ':' << type << ':' << content ;
				_tags.push_back( s.str() ) ;
			}			

			


		protected:

			void set_to_defaults( ) {
				_rnext = "*" ;
				_pnext = 0 ;
				_tlen  = 0 ;
			} 

			void set_to_reverse_strand() {
				setReverseComplemented() ;
				sequence( rc_sequence() ) ;
				quality( r_quality() ) ;
			}

		} ;

		class SAMHeader {
		
			std::vector<std::string> _names ;
			std::vector<int> _lengths ;

		public:
			SAMHeader() ;

			void add( std::string nm, int l ) ;

			/* 
			 Returns the formatted SAM header  
			 */
			std::string str() const ;

			std::vector<std::string> names() const ;

			std::vector<int> lengths() const ;
		} ;
		
		
	}
}