#pragma once

#include "stdafx.h"
#include "Read.h"
#include "Amplicon.h"
#include "Alignment.h"
#include "SAMrecord.h"

namespace Nimbus {

	class AlnSet {

	public:
		basic::Amplicon* amplicon ;
		alignment::Alignment* f_alignment ;
		alignment::Alignment* r_alignment ;
		std::vector< std::pair<int,int> >* f_path ;
		std::vector< std::pair<int,int> >* r_path ;
		alignment::SAMRecord* f_record ;
		alignment::SAMRecord* r_record ;

	public:
		AlnSet() ;
		AlnSet( basic::Amplicon* a ) ;

		~AlnSet() ;

		void delete_content() ;

		void align( alignment::AlignmentScore* scores, int gapopen, basic::Read* f ) ;

		void align( alignment::AlignmentScore* scores, int gapopen, basic::Read* f, basic::Read* r ) ;

		void SAMrecord( basic::Read* f, basic::Read* r ) ;

		void SAMrecord_f( basic::Read* f ) ;

		void SAMrecord_r( basic::Read* r ) ;
		
	} ;
	 
	class AlignmentBuilder {

	public:
		basic::Read* forward ;
		basic::Read* reverse ;
		std::vector<AlnSet> entries ;

	public:
		AlignmentBuilder( ) ;

		AlignmentBuilder( basic::Read* f ) ;

		AlignmentBuilder( basic::Read* f, basic::Read* r ) ;
		
		~AlignmentBuilder(void);

		bool empty() const ;

		/*
		 adds an amplicon to the result set
		 */
		void add( basic::Amplicon* a ) ;

		/*
		 Aligns the read to the amplicons in the resultset 
		 */
		void align( alignment::AlignmentScore* scores, int gapopen ) ;

		/*
		 gets the alignment with the best combined score
		 */
		int best()  ;

		/*
		 Creates a forward (and if possible  a reverse) SAMrecord
		 */
		void createRecord( AlnSet& a ) ; 

		/* 
		 Creates SAMRecords for each alignment
		 */ 
		void createRecords( ) ; 
		
		bool samrecordspresent( ) const ;

	} ;

}