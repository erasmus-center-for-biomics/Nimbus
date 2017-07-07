#include "stdafx.h"
#include "io.h"
#include "Utils.h"
#include "Amplicon.h"

namespace Nimbus {

	namespace IO {
		
		using namespace std ;
		using namespace basic ;
		using namespace alignment ;

		Read* FastQReader( istream& input ) {

			// declare the return value
			Read* rval = NULL ;
			
			string pname = "" ;
			string seq   = "" ;
			string sname = "" ;
			string qual  = "" ;

			// we can read a ine
			while( input.good() ) { 
				string tmp ;
				getline( input, tmp ) ;

				// cycle through the variables until: 
				// pname.size() > 0 and pname[0] == '@' && sname[0] == '+'  
				pname = seq ;
				seq   = sname ; 
				sname = qual ;
				qual  = tmp ;

				// break as soon as we have something that looks like a FastQ entry 
				if( pname.size() > 0 && pname[0] == '@' && sname.size() > 0 && sname[0] == '+' ) {
					break ;
				}

			}

			// do some more checks
			if( seq.size() > 0 && seq.size() == qual.size() ) {
				rval = new Read( pname.substr(1), seq, qual ) ;
			}

			// return the read pointer
			return rval ;
		}


		Read* FastQReader( istream& input, utils::AdapterTrim* a ) {
			Read* r = FastQReader( input ) ;
			Read* rval = a->trim( r ) ;
			delete r ;
			return rval ;
		}

		pair< string, string* >* FastAReader( istream& input ) {

			// declare the return object
			pair< string, string* >* rval = NULL ; 

			// declare the temporary objects
			string name = "" ;
			string* seq = new string("") ;

			// while the input can be read
			while( input.good() ) {
				// get the previous point in the file
				streamoff p = input.tellg() ;	

				// get a new string from the input file
				string tmp ;
				getline( input, tmp ) ;

				if( tmp.size() > 1 && tmp[0] == '>' && name == "" ) {
					// set the name of the sequence
					name = tmp.substr(1) ;					
				} else if( tmp.size() > 1 && tmp[0] == '>' && name != "" ) {
					// if we find a second sequence name rewind and break
					//
					// we will only rewind once per chromosome
					input.seekg(p) ;
					break ;
				} else if( tmp.size() > 0 ) {
					// convert the sequence to uppercase
					transform( tmp.begin(), tmp.end(), tmp.begin(), ::toupper ) ;

					// append tmp to the sequence
					seq->append( tmp ) ;
				} 				
			}

			//
			if( name == "" && seq->size() == 0 && rval == NULL ) {
				delete seq ;				
			} else if( name != "" && seq->size() > 0 ) {
				rval = new pair< string, string* >(name, seq) ;				
			} 

			// return the chromosome pair
			return rval ;
		}


		GenomicRegion* BEDReader( istream& input ) {

			// prepare the return value
			GenomicRegion* rval = NULL ;

			// declare the fields vector
			vector<string> fields = vector<string>() ; 
				
			// get a line from the outp
			while( input.good() ) {	
				//
				string tmp ;
				getline( input, tmp ) ;

				// get the fields
				fields = utils::split_string( tmp, "\t" ) ;
				if( fields.size() >= 3 ) break ;				
			}

			// interpret the bed region
			if( fields.size() == 3 ) {
				rval = new GenomicRegion( fields[0], atoi( fields[1].c_str() ), atoi( fields[2].c_str() ), true ) ;
			} else if( fields.size() == 4) {
				rval = new GenomicRegion( fields[0], atoi( fields[1].c_str() ), atoi( fields[2].c_str() ), true,  fields[3] ) ;
			} else if( fields.size() >= 6) {
				bool f = fields[5] == "-" ? false : true ;
				rval = new GenomicRegion( fields[0], atoi( fields[1].c_str() ), atoi( fields[2].c_str() ), f,  fields[3] ) ;
			}

			// return the pointer or NULL
			return rval ;
		}

		vector<GenomicRegion*> BEDReader( string fname ) {

			// initialize the return value
			vector<GenomicRegion*> rval = vector<GenomicRegion*>() ;

			// initialize the input handle with the currrent file
			ifstream handle( fname.c_str(), ios::in ) ;

			if( handle.is_open() ) {
				// while we have input 
				GenomicRegion* region = BEDReader( handle ) ;
				while( region != NULL ) {
					rval.push_back( region ) ;
					region = BEDReader( handle ) ;
				}
				
				// close the input handle
				handle.close() ;
			}

			// return the return value
			return rval ; 
		}


		void AmpliconReader( vector<Amplicon*>& amplicons, SAMHeader& header, string fname_bed, string fname_fasta ) {

			// declare the result vector
			amplicons = vector<Amplicon*>() ;
			header    = SAMHeader() ; 

			// load all the bed regions
			vector<GenomicRegion*> bed = BEDReader( fname_bed ) ;
			
			// initialize the input handle with the currrent file
			ifstream fasta( fname_fasta.c_str(), ios::in ) ;

			if( fasta.is_open() ) {

				// iterate over the fasta sequences
				pair< string, string* >* fseq = FastAReader( fasta ) ;
				while( fseq != NULL ) {

					// add the reference to the sam header
					header.add( fseq->first, (int) fseq->second->size() ) ;

					// get the genomic regions for the current chromosome
					vector<GenomicRegion*> cur = vector<GenomicRegion*>( ) ;
					for( vector<GenomicRegion*>::iterator it=bed.begin(); it!=bed.end(); ++it ) {
						GenomicRegion* x = *it ;
						if( x->chromosome() == fseq->first ) {
							cur.push_back( *it ) ;
						}
					}

					// Build the amplicon list
					for( vector<GenomicRegion*>::iterator it=cur.begin(); it!=cur.end(); ++it ) {
						GenomicRegion gr = *(*it) ;
						string ampseq    = fseq->second->substr( gr.start(), gr.width() ) ;
						amplicons.push_back( new Amplicon( gr, ampseq ) ) ; 
					}

					// clean the memory before we move on
					delete fseq->second ;
					delete fseq ;

					// prepare for the next round
					fseq = FastAReader( fasta ) ;
				}
				fasta.close() ;
			}

			// delete the genomic regions
			for( vector<GenomicRegion*>::iterator it=bed.begin(); it!=bed.end(); ++it ) {				
				delete *it ;
			}
		
		}

	}
}
