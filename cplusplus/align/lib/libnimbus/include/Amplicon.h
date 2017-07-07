#pragma once


namespace Nimbus {

	namespace basic {
	

		class GenomicRegion {
		protected: 
			std::string _chr ;
			int _start ;
			int _end ;
			bool _forward ;
			std::string _name ;
		public:
			GenomicRegion( std::string chr, int start, int end, bool forward ): _chr(chr), 
																				_start(start), _end(end), 
																				_forward(forward), _name("") {
			}

			GenomicRegion( std::string chr, int start, int end, bool forward, std::string nm ): _chr(chr), 
																				_start(start), _end(end), 
																				_forward(forward), _name(nm) {
			}
			
			GenomicRegion( const GenomicRegion& other ):_chr(other.chromosome()), 
														_start(other.start()), _end(other.end()), 
														_forward(other.forward()), _name(other.name())  {

			}

			~GenomicRegion(){} 



			// constant getters
			std::string chromosome() const ;
			int start() const ;
			int end() const ; 
			bool forward() const ;
			std::string name() const ;
			int width() const ;
			std::string str() const ;
			std::string str( bool nm ) const ;

			std::string format() const ;

			// overloaded relational operators for sorting
			bool operator<( const GenomicRegion& a ) const ;
			bool operator>( const GenomicRegion& a ) const ;
			bool operator==( const GenomicRegion& a ) const ;
		} ; 


		class Amplicon: public GenomicRegion {

			
			std::string _sequence ;
			
		public:

			// default constructor
			Amplicon( std::string chr, int start, int end, bool forward, std::string sequence ) : GenomicRegion(chr, start, end, forward), _sequence(sequence) { }

			Amplicon( std::string chr, int start, int end, bool forward, std::string sequence, std::string name ) : GenomicRegion(chr, start, end, forward, name), _sequence(sequence) { }

			Amplicon( GenomicRegion g, std::string sequence ) : GenomicRegion(g), _sequence(sequence) { }

			// copy constructor
			Amplicon( const Amplicon& other ): GenomicRegion( (GenomicRegion)other ), _sequence(other.sequence() ){
			}
		
			// standard destroyer
			~Amplicon(void);

			// constant getters
			std::string sequence() const ;
			
			
			std::string str( ) const ;
		} ;

		/**
		 pointer comparators
	 	 **/
		bool cmp_lt_amplicon_p( Amplicon* a, Amplicon* b ) ;
		bool cmp_gt_amplicon_p( Amplicon* a, Amplicon* b ) ;
		bool cmp_eq_amplicon_p( Amplicon* a, Amplicon* b ) ;

	}

	
}