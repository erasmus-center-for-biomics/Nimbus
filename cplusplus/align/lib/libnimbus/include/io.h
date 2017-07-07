#pragma once

#include "stdafx.h"
#include "Read.h"
#include "Amplicon.h"
#include "AdapterTrim.h"
#include "SAMrecord.h"

namespace Nimbus {

	namespace IO {
		/*
		 * A FastQ file reader that works similar to a python generator. The 
		 * stop signal is a NULL object  
		 */
		basic::Read* FastQReader( std::istream& input ) ;  

		/*
		 * Same as the basicv FastQ file reader, but this one also trims the adapter 
		 * sequences from the read
		 */
		basic::Read* FastQReader( std::istream& input, utils::AdapterTrim* a ) ;  

		/*
		 * A FastA file reader
		 */
		std::pair< std::string, std::string* >* FastAReader( std::istream& input ) ;

		/*
		 * A basic BED file reader that returns a single GenomicRegion instance from an input stream
		 */
		basic::GenomicRegion* BEDReader( std::istream& input ) ;

		/*
		 *  A basic BED file reader that returns all the regions in a BED file
		 */
		std::vector<basic::GenomicRegion*> BEDReader( std::string fname ) ;

		/* 
		 * Reads the amplicons from the bed and fasta file
		 */
		void AmpliconReader( std::vector<basic::Amplicon*>& amplicons, alignment::SAMHeader& header, std::string fname_bed, std::string fname_fasta ) ;

		/*
		 * Reads the amplicons from a combined index file
		 */
		//void AmpliconReader( std::vector<basic::Amplicon*>& amplicons, alignment::SAMHeader& header, std::string fname_index ) ;

	}
}