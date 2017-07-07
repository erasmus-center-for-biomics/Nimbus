#include "stdafx.h"
#include "Alignment.h"


namespace Nimbus {

	namespace alignment {

		//
		//
		// Implementation of the AlignmentScore object
		//
		//

		

		/** 
		 * Get the score with reference r and query q
		 *
		 **/
		int AlignmentScore::score( char r, char q ) const {

			int rval ;
			if( r == '-' || q == '-' ) {
				rval = _gap ;
			} else if( r == 'N' || q == 'N' ) {
				rval = _match ;
			} else if( r == q ) {
				rval = _match ;
			} else {
				rval = _mismatch ;
			}
			return rval ;
		}		

		std::string AlignmentScore::str() const { 
			
			std::stringstream ss ;

			ss << "match: " << _match << " mismatch: " << _mismatch << " gap " << _gap ;
			return ss.str() ;
		}

		//
		//
		// Implementation of the Alignment object
		//
		//
		Alignment::~Alignment(void) {
			_clean_matrix() ;
		}

		int Alignment::getAlignmentScore( ) const {
			return getAlignmentScore( bestAlignment() ) ;
		}

		int Alignment::getAlignmentScore( std::pair<int,int> xy ) const {
			return getAlignmentScore( xy.first, xy.second ) ;
		}

		/*
		 Get the scores at position x and y
		 */
		int Alignment::getAlignmentScore( int x, int y) const {
			int rval = -1 ;
			if( scores != NULL && x < (int)subject.size() + 1 && y < (int)query.size() + 1 ) {
				rval = scores[x][y] ;
			}
			return rval ;
		}


		//
		// pathing algorithms
		//

		/**
		  Get the optimal path through the alignment matrix
		 **/
		std::vector< std::pair<int,int> >* Alignment::getPath( int x, int y ) const {

			// declare the output vector
			std::vector< std::pair<int,int> >* rval = new std::vector< std::pair<int,int> >() ;

			// if we have scores and directions
			if( scores != NULL && directions != NULL ) {

				// get the x and y coordinates 
				int cur_i = x ;
				int cur_j = y ;

				// while we did not hit a bound
				while( directions[cur_i][cur_j] != d_BOUND ) {
					// add the coordinates to the vector
					rval->push_back( std::pair<int,int>( cur_i, cur_j ) ) ;

					// get the directions
					switch( directions[cur_i][cur_j] ) {
					case d_DIAG:
						cur_i-- ;
						cur_j-- ;						
						break ;
					case d_VERTICAL:
						cur_j-- ;						
						break ;
					case d_HORIZONTAL:
						cur_i-- ;						
						break ;
					default:
						break ;
					}										
				}
				// rval.push_back( std::pair<int,int>( cur_i, cur_j ) ) ;
			}

			// reverse the path, so that the 5` entry goes first
			std::reverse( rval->begin(), rval->end() ) ;
			return rval ;
		}

		std::vector< std::pair<int,int> >* Alignment::getPath( ) const {
			return getPath( _coord.first, _coord.second ) ;
		}

		std::vector< std::pair<int,int> >* Alignment::getPath( std::pair<int,int> xy ) const {
			return getPath( xy.first, xy.second ) ;
		}

		/*
		 Convert the path to a CIGAR string
		 */
		std::vector<char> Alignment::QCIGAR( std::vector< std::pair<int,int> > path ) const {

			// declare the return value
			std::vector<char> rval = std::vector<char>() ; 
			
			// return an empty cigar if the directions were not defined
			if( directions == NULL ) {
				return rval ;
			}

			// matrix coordinates
			int cur_i ;
			int cur_j ;

			// foreach pair in the path
			for( unsigned int i=0; i<path.size(); i++ ) {

				// traverse path
				cur_i = path[i].first ;
				cur_j = path[i].second ;

				// make sure the path is compatible with the stored sequences
				if( cur_i > (int) subject.size() + 1 && cur_j > (int) query.size() + 1 ) {

					// otherwise clear the output and break the for loop
					rval.clear() ;
					break ;
				}

				// we did not get a path back to the first base of the query
				if( i == 0 && ( cur_j != 1) ) {
					// add soft clipped bases to the alignment 
					for( int x=1; x<cur_j; x++ ) rval.push_back( 'S' ) ;
				}
				
				// add the CIGAR characters for match, insertion and deletion
				if( directions[cur_i][cur_j] == d_DIAG ) {
					rval.push_back( 'M' ) ; 
				} else if( directions[cur_i][cur_j] == d_HORIZONTAL ) {					
					rval.push_back( 'D' ) ;
				} else if( directions[cur_i][cur_j] == d_VERTICAL ) {
					rval.push_back( 'I' ) ;
				}

				// if the alignment did not span to the end of the read
				if( i == path.size() - 1 && cur_j != (int)query.size() ) {

					// add soft clipped bases
					for( int x=cur_j; x<(int)query.size(); x++ ) rval.push_back( 'S' ) ;
				}
			}

			// return the CIGAR character vector
			return rval ;
		}

		std::vector<char> Alignment::QCIGAR( ) const {
			std::vector< std::pair<int,int> >* path = getPath() ;
			std::vector<char> rval = QCIGAR( *path ) ;
			delete path ;
			return rval ;
		}

		/**
		 * clean the matrix
		 *
		 **/
		void Alignment::_clean_matrix() {

			// clear the scores
			if( scores != NULL ) {
				for( unsigned int i=0; i<subject.size()+1; i++ ) {
					delete[] scores[i] ;
				}
				delete[] scores ;
				scores = NULL ;
			}

			// clear the directions matrix
			if( directions != NULL ) {
				for( unsigned int i=0; i<subject.size()+1; i++ ) {
					delete[] directions[i] ;
				}
				delete[] directions ;
				directions = NULL ;
			}
			subject = "" ;
			query   = "" ;
		}

		std::string Alignment::formatScoreMatrix( ) {
			std::string rval = "" ;
			std::stringstream ss ;

			if( scores != NULL ) {

				// set the counters
				int n_subj = (int) subject.size() + 1 ;
				int n_qry  = (int) query.size() + 1 ;

				// header 
				ss << '\t' << '-' ; 
				for( int j=1; j<n_qry; j++ ) {
					ss << '\t' << j << ":" << query[j-1] ;
				}
				ss << '\n' ;
				// traverse the subject than the query
				for( int i=0; i<n_subj; i++ ) {
					if( i == 0 ) {
						ss << i << ':' << '-';
					} else {
						ss << i << ':' << subject[i-1]  ;
					}
					//
					for( int j=0; j<n_qry; j++ ) {
						 ss << '\t' << scores[i][j] ;
					}
					ss << '\n' ;
				}
				rval = ss.str() ;
			}
			
			return rval ;
		}


		std::string Alignment::formatDirectionMatrix( ) {
			std::string rval = "" ;
			std::stringstream ss ;

			if( directions != NULL ) {

				// set the counters
				int n_subj = (int) subject.size() + 1 ;
				int n_qry  = (int) query.size() + 1 ;

				// header 
				ss << '\t' << '-' ; 
				for( int j=1; j<n_qry; j++ ) {
					ss << '\t' << j << ":" <<  query[j-1] ;
				}
				ss << '\n' ;
				// traverse the subject than the query
				for( int i=0; i<n_subj; i++ ) {
					if( i == 0 ) {
						ss << i << ':' << '-';
					} else {
						ss << i << ':' << subject[i-1]  ;
					}
					//
					for( int j=0; j<n_qry; j++ ) {
						 ss << '\t' << directions[i][j] ;
					}
					ss << '\n' ;
				}
				rval = ss.str() ;
			}
			
			return rval ;
		}

		std::string Alignment::str() const {
			std::stringstream rval ;

			//
			if( scores != NULL && directions != NULL ) {
				rval << " -" ;
				for( unsigned int j=1; j<query.size()+1; j++ ) {
					rval << " " << query[j-1] ;
				}
				rval << std::endl ;

				//
				for( unsigned int i=0; i<subject.size()+1; i++ ) {
					//
					for( unsigned int j=0; j<query.size()+1; j++ ) {

						// before we write any scores
						if( j == 0 ) {
							// print the subject base or -
							if( i == 0 ) {
								rval << "-"  ;
							} else {
								rval << subject[i-1] ;
							}
						}
						// print the score and direction at coordinate i,j
						rval << " " << scores[i][j] << "(" << directions[i][j] << ")" ;						
					}

					// end the line
					rval << std::endl ;
				}
			}
			return rval.str() ;

		}

		//
		//
		// Specific alignments: Smith-Waterman
		//
		//

		/*
		 * fills the score and direction matrix according to the Smith-Waterman method.
		 * 
		 * The 0-indexed column and row will be filled with 0's. The content of the matrix 
		 * will be filled using the maximum of {0, score gap X, score gap Y, or score diagonal}
		 * 
		 */
		void SmithWaterman::fillMatrix( std::string s, std::string q ) {
			// set the class parameters
			subject = s ;
			query   = q ;
			_max    = 0 ;
						
			// initialize the matrix
			_init_matrix() ;

			// set the counters
			int n_subj = (int) subject.size() + 1 ;
			int n_qry  = (int) query.size() + 1 ;

			// traverse the subject than the query
			for( int i=1; i<n_subj; i++ ) {
				for( int j=1; j<n_qry; j++ ) {
					
					//
					//                   |   scores[i-1, j ] + scorecalc->score( '-', query[j-1]) + if( directions[i-1,j] != d_HORIZONTAL, _gapopen, 0 ) |
					//  scores[i,j] = max|   scores[ i ,j-1] + scorecalc->score(subject[i-1], '-') + if( directions[i,j-1] != d_VERTICAL, _gapopen, 0 )  |
					//                   |   scores[i-1,j-1] + scorecalc->score(subject[i-1], query[j-1])                                                |
					//                   |   0                                                                                                           |
					//
					//
					//  directions[i,j] = which_max( score_calc[i,j] )
					//

					// our scoring functions: base       +      per base score                   +    gap open penalty
					int _s_horizontal = scores[i-1][j]   + scorecalc->score( '-', query[j-1] )   + (directions[i-1][j] != d_HORIZONTAL ? _gapopen : 0) ;
					int _s_vertical   = scores[i][j-1]   + scorecalc->score( subject[i-1], '-' ) + (directions[i][j-1] != d_VERTICAL ? _gapopen : 0)  ;
					int _s_diagonal   = scores[i-1][j-1] + scorecalc->score( subject[i-1], query[j-1] ) ; 
				
					//printf("i:%d, j:%d, score: %d,horizontal: %d, vertical: %d, diagonal: %d\n", i, j, scorecalc->score( subject[j-1], query[j-1] ), _s_horizontal, _s_vertical, _s_diagonal) ;

					// Fill the matrix and decide which scores to record
					int lmax = _s_diagonal ;
					if( lmax < _s_horizontal ) lmax = _s_horizontal ;
					if( lmax < _s_vertical ) lmax = _s_vertical ;
					if( lmax < 0 ) lmax = 0 ;

					if( _s_diagonal == lmax ) {
						scores[i][j]     = _s_diagonal ;
						directions[i][j] = d_DIAG ;
					} else if ( _s_vertical == lmax) { 
						scores[i][j]     = _s_vertical ;
						directions[i][j] = d_VERTICAL ;
					} else if ( _s_horizontal == lmax) { 
						scores[i][j]     = _s_horizontal ;
						directions[i][j] = d_HORIZONTAL ;					
					} else if ( 0 == lmax) { 
						scores[i][j] = 0 ;
						directions[i][j] = d_BOUND ;
					}
					/*
					scores[i][j] = 0 ;
					directions[i][j] = d_BOUND ;
					if( _s_diagonal > scores[i][j] ) {
						scores[i][j]     = _s_diagonal ;
						directions[i][j] = d_DIAG ;
					} else if( _s_vertical > scores[i][j] ) {
						scores[i][j]     = _s_vertical ;
						directions[i][j] = d_VERTICAL ;
					} else if( _s_horizontal > scores[i][j] ) {
						scores[i][j]     = _s_horizontal ;
						directions[i][j] = d_HORIZONTAL ;
					}
					*/
					// record the top positions
					if( scores[i][j] > _max ) {
						_max   = scores[i][j] ;
						_coord = std::pair<int,int>( i, j ) ;
					}
				}
			}
			// return void
		}

		/*
		 * Initialize the matrix with 0's
		 */
		void SmithWaterman::_init_matrix( ) {

			// get the number of rows and columns in the matrix
			int n_subj = (int) subject.size() + 1 ;
			int n_qry  = (int) query.size() + 1 ;

			// clean the matrix if required
			if( scores != NULL || directions != NULL ) _clean_matrix() ;

			//initialize the matrix
			scores     = new int*[n_subj] ;
			directions = new int*[n_subj] ;

			// fill the matrix
			for( int i=0; i<n_subj; i++ ) {

				// add a new row
				scores[i]     = new int[n_qry] ;
				directions[i] = new int[n_qry] ;

				// fill the rows with 0's
				for( int j=0; j<n_qry; j++ ) {
					scores[i][j]     = 0 ;
					directions[i][j] = d_BOUND ;
				}
			}
		}

		//
		//
		// Specific alignments: Needleman Wunsch
		//
		//

		
		/*
		 * fills the score and direction matrix according to the Needleman-Wunsch method.
		 * 
		 * The 0-indexed column and row will be filled with 0's. The content of the matrix 
		 * will be filled using the maximum of {0, score gap X, score gap Y, or score diagonal}
		 * 
		 */
		void NeedlemanWunsch::fillMatrix( std::string s, std::string q ) {
			
			// set the class parameters
			subject = s ;
			query   = q ;
			_max    = 0 ;
						
			// initialize the matrix
			_init_matrix() ;

			// set the counters
			int n_subj = (int) subject.size() + 1 ;
			int n_qry  = (int) query.size() + 1 ;

			// traverse the subject than the query
			for( int i=1; i<n_subj; i++ ) {
				for( int j=1; j<n_qry; j++ ) {
					
					//
					//                   |   scores[i-1, j ] + scorecalc->score( '-', query[j-1]) + if( directions[i-1,j] != d_HORIZONTAL, _gapopen, 0 ) |
					//  scores[i,j] = max|   scores[ i ,j-1] + scorecalc->score(subject[i-1], '-') + if( directions[i,j-1] != d_VERTICAL, _gapopen, 0 )  |
					//                   |   scores[i-1,j-1] + scorecalc->score(subject[i-1], query[j-1])                                                |					
					//
					//  directions[i,j] = which_max( score_calc[i,j] )
					//

					// our scoring functions: base     +      per base score                    +    gap open penalty
					int _s_horizontal = scores[i-1][j] + scorecalc->score( '-', query[j-1] )    + directions[i-1][j] != d_HORIZONTAL ? _gapopen : 0 ;
					int _s_vertical   = scores[i][j-1] + scorecalc->score( subject[i-1], '-' )  + directions[i][j-1] != d_VERTICAL ? _gapopen : 0  ;
					int _s_diagonal   = scores[i-1][j-1] + scorecalc->score( subject[i-1], query[j-1] ) ; 
				
					// Fill the matrix and decide which scores to record: default is the diagonal					
					scores[i][j]     = _s_diagonal ;
					directions[i][j] = d_DIAG ;
					if( _s_vertical > scores[i][j] ) {
						scores[i][j]     = _s_vertical ;
						directions[i][j] = d_VERTICAL ;
					} else if( _s_horizontal > scores[i][j] ) {
						scores[i][j]     = _s_horizontal ;
						directions[i][j] = d_HORIZONTAL ;
					}

					// record the top positions
					if( scores[i][j] > _max ) {
						_max   = scores[i][j] ;
						_coord = std::pair<int,int>( i, j ) ;
					}
				}
			}
			// return void
		}

		/*
		 * Initialize the matrix with 0's
		 */
		void NeedlemanWunsch::_init_matrix( ) {

			//
			int n_subj = (int) subject.size() + 1 ;
			int n_qry  = (int) query.size() + 1 ;

			// clean the matrix if required
			if( scores != NULL || directions != NULL ) _clean_matrix() ;

			//initialize the matrix
			scores     = new int*[n_subj] ;
			directions = new int*[n_subj] ;

			// fill the matrix
			for( int i=0; i<n_subj; i++ ) {

				// add a new row
				scores[i]     = new int[n_qry] ;
				directions[i] = new int[n_qry] ;

				// fill the rows with 0's
				for( int j=0; j<n_qry; j++ ) {
					scores[i][j]     = 0 ;
					directions[i][j] = d_BOUND ;
				}
			}

			// set the important 0 row and 0 column bounds
			for( int i=0; i<n_subj; i++ ) { 
				scores[i][0] = i * _gapopen ;
			}
			for( int j=0; j<n_qry; j++ ) { 
				scores[0][j] = j * _gapopen ;
			}
		}

		//
		//
		// Specific alignments: Levenshtein
		//
		//

		void Levenshtein::fillMatrix( std::string s, std::string q ) {
			// set the class parameters
			subject = s ;
			query   = q ;
			_max    = 0 ;
						
			// initialize the matrix
			_init_matrix() ;

			// set the counters
			int n_subj = (int) subject.size() + 1 ;
			int n_qry  = (int) query.size() + 1 ;

			// traverse the subject than the query
			for( int i=1; i<n_subj; i++ ) {
				for( int j=1; j<n_qry; j++ ) {

					//
					//                   |   scores[i-1, j ] + scorecalc->score( '-', query[j-1])         |
					//  scores[i,j] = min|   scores[ i ,j-1] + scorecalc->score(subject[i-1], '-')        |
					//                   |   scores[i-1,j-1] + scorecalc->score(subject[i-1], query[j-1]) |					
					//
					//  directions[i,j] = which_min( score_calc[i,j] )
					//

					//
					int _s_horizontal = scores[i-1][j] + scorecalc->score( '-', query[j-1] ) ;
					int _s_vertical   = scores[i][j-1] + scorecalc->score( subject[i-1], '-' ) ;
					int _s_diagonal   = scores[i-1][j-1] + scorecalc->score( subject[i-1], query[j-1] ) ; 

					// Fill the matrix and decide which scores to record: default is the diagonal					
					scores[i][j]     = _s_diagonal ;
					directions[i][j] = d_DIAG ;
					if( _s_vertical < scores[i][j] ) {
						scores[i][j]     = _s_vertical ;
						directions[i][j] = d_VERTICAL ;
					} else if( _s_horizontal < scores[i][j] ) {
						scores[i][j]     = _s_horizontal ;
						directions[i][j] = d_HORIZONTAL ;
					}
					
				}
			}
		}

		/**

		 **/
		int Levenshtein::distance() {

			// the return value 
			int rval = -1 ;
			
			// if we have defined scores and directions matrices we can proceed
			if( scores != NULL && directions != NULL ) {

				// get the end of the subject and the query
				int n_subj = (int) subject.size() ;
				int n_qry  = (int) query.size() ;

				// set the distance 
				rval = scores[n_subj][n_qry] ;
			}

			// return the distance
			return rval ;
		}

		void Levenshtein::_init_matrix( ) {
			int n_subj = (int) subject.size() + 1 ;
			int n_qry  = (int) query.size() + 1 ;

			// clean the matrix if required
			if( scores != NULL || directions != NULL ) _clean_matrix() ;

			//initialize the matrix
			scores     = new int*[n_subj] ;
			directions = new int*[n_subj] ;

			// fill the matrix
			for( int i=0; i<n_subj; i++ ) {

				// add a new row
				scores[i]     = new int[n_qry] ;
				directions[i] = new int[n_qry] ;

				// fill the rows with 0's
				for( int j=0; j<n_qry; j++ ) {
					scores[i][j]     = 0 ;
					directions[i][j] = d_BOUND ;
				}
			}

			// set the important 0 row and 0 column bounds
			for( int i=0; i<n_subj; i++ ) { 
				scores[i][0] = i ;
			}
			for( int j=0; j<n_qry; j++ ) { 
				scores[0][j] = j ;
			}
		}
	}
}