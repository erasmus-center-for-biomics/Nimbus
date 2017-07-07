#pragma once

namespace Nimbus {

	namespace basic {

		class Read {
			
		protected:
			std::string _name ;
			std::string _seq ;
			std::string _qual ;
			std::string _rc_seq ;
			std::string _r_qual ;

		public:
			Read( std::string name, std::string sequence, std::string qualiy ) ;
			
			Read() ;

			Read( const Read& other ) ;

			~Read(void);
			
			/* get the reverse complement of the sequence */
			std::string rc_sequence() ;

			/* get the reverse complement of the quality */
			std::string r_quality() ;

			/* get the reverse complement for this read */
			std::string name() const ;

			/* get the sequence of the read */
			std::string sequence() const ;

			/* get the quality string of the read */
			std::string quality() const ;

			/* return the read in FastQ format */
			std::string fastq() const ;
	
			unsigned int size() const ;

			std::string str( ) const ;
			

		protected:
			/* set the sequence of the read */
			void sequence( std::string s ) ;

			/* set the quality string of the read */
			void quality( std::string q ) ;
			
		};

	}
}