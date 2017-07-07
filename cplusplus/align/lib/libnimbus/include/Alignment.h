#pragma once

namespace Nimbus {

	namespace alignment {

		//
		// A class to 
		//

		

		class AlignmentScore {

		public:
			// the scores
			int _match ;	// match
			int _mismatch ;	// mismatch
			int _gap ;		// gap
			int _maxamp ;   // maximum number of amplicons
		public:
			AlignmentScore( int m, int mm, int g, int maxamp): _match(m), _mismatch(mm), _gap(g), _maxamp(maxamp) {}
			~AlignmentScore( ) {}

			/** 
			 Get the score with reference r and query q
			 **/
			int score( char r, char q ) const ;

			std::string str() const ;
		} ;



		// define the directions
		enum directions_t { d_BOUND, d_DIAG, d_VERTICAL, d_HORIZONTAL } ;
		//int D_BOUND = 0 ;
		//int D_DIAG  = 1 ;
		//int D_DOWN  = 2 ;
		//int D_RIGHT = 3 ;

		class Alignment {
			/*
			The base class of the classes implementing alignment algorithms

			Alignment matrices are represented as follows:

			  - S u b j e c t
			- # # # # # # # #
			Q # # # # # # # #
			u # # # # # # # #
			e # # # # # # # #
			r # # # # # # # #
			y # # # # # # # #

			 */
		protected:

			// the score matrix
			int** scores ;

			// 4 directions:
			//	- 0 boundary
			//	- 1, match or mismatch
			//	- 2, subject -, query base
			//	- 3, subject base, query - 			
			int** directions ;

			// the subject and query functions
			std::string subject ;
			std::string query ;

			// Variables for the trace back
			int _max ;					// maximum score 
			std::pair<int,int> _coord ;	// query and subject positions of the maximum score
			
			AlignmentScore* scorecalc ;

		public:		
			Alignment( AlignmentScore* as ) : scorecalc(as) {
				subject    = "" ;
				query      = "" ;
				scores     = NULL ;
				directions = NULL ;				
				_max       = 0 ;
				_coord     = std::pair<int,int>( 0, 0 ) ;
			}

			~Alignment(void) ; 

			/**
			 Get the path through the matrix for the best alignment
			 **/ 
			std::vector< std::pair<int,int> >* getPath( ) const ;

			/** 
			 Get the path through the matrix for the alignment starting at coordinate x and y in the matrix
			 **/
			std::vector< std::pair<int,int> >* getPath( int x, int y ) const ;
			std::vector< std::pair<int,int> >* getPath( std::pair<int, int> coord ) const ;

			/** 
			 Get the end coordinates of the best alignment
			 **/
			std::pair<int,int> bestAlignment() const {
				return _coord ;
			}

			/**
			 Get the maximum alignment score
			 **/
			int getAlignmentScore( ) const ;			

			/**
			 Get the alignment score at coordinate x and y
			 **/
			int getAlignmentScore( int x, int y ) const ;
			int getAlignmentScore( std::pair<int, int> coord ) const ;

			/** 
			 Get the CIGAR string for the Query against the reference in a character vector
			 **/ 
			std::vector<char> QCIGAR() const ;		
		
			std::vector<char> QCIGAR( std::vector< std::pair<int,int> > path ) const ;

			/**
			 This method is specific foreach alignment method
			 **/
			virtual void fillMatrix( std::string ref, std::string q ) {}

			std::string formatScoreMatrix( ) ;
			std::string formatDirectionMatrix( ) ;

			/**
			 * Returns a presentation of the matrix
			 **/
			std::string str() const ;

			
		protected:
			void _clean_matrix( ) ; 
		} ;


		class SmithWaterman: public Alignment {
			int _gapopen ;

		public:
			SmithWaterman( AlignmentScore* as, int go ): Alignment(as), _gapopen(go) {}

			SmithWaterman( AlignmentScore* as, int go, std::string s, std::string q ): Alignment(as), _gapopen(go) {
				fillMatrix( s, q ) ;
			}

			/*~SmithWaterman() {
				_clean_matrix() ;
			}*/

			void fillMatrix( std::string ref, std::string q ) ;

			void _init_matrix( ) ;
		} ;

		class NeedlemanWunsch: public Alignment {
			int _gapopen ;

		public:
			NeedlemanWunsch( AlignmentScore* as, int go ): Alignment(as), _gapopen(go) {}

			NeedlemanWunsch( AlignmentScore* as, int go, std::string s, std::string q ): Alignment(as), _gapopen(go) {
				fillMatrix( s, q ) ;
			}

			/*~NeedlemanWunsch() {
				_clean_matrix() ;
			}*/

			void fillMatrix( std::string ref, std::string q ) ;

			void _init_matrix( ) ;
		} ;

		class Levenshtein: public Alignment {
		public:
			Levenshtein( ): Alignment(new AlignmentScore( 0, 1, 1, 1)) {
			}
			
			Levenshtein( std::string ref, std::string q ): Alignment(new AlignmentScore( 0, 1, 1, 1)) {
				fillMatrix( ref, q ) ;
			}

			Levenshtein( const Levenshtein& other ): Alignment(new AlignmentScore( 0, 1, 1, 1))  {
				fillMatrix( other.subject, other.query ) ;
			}

			~Levenshtein( ) {
				delete scorecalc ;
			//	_clean_matrix() ;
			}

			int distance() ;

			void fillMatrix( std::string ref, std::string q ) ;

			void _init_matrix( ) ;

		} ;
	}

}
