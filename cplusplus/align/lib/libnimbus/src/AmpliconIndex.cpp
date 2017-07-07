#include "stdafx.h"
#include "Utils.h"
#include "AmpliconIndex.h"

namespace Nimbus {

	namespace seed {

		using namespace std ;
		using namespace basic ;

		//
		//
		// AmpliconIndex implementation
		//
		//

		AmpliconIndex::AmpliconIndex(void) {
			_keysize = 7 ;
			_idx_f_a = NULL ;
			_idx_r_a = NULL ;
			_idx_f_b = NULL ;
			_idx_r_b = NULL ;
			vector<Amplicon*> _amplicons = vector<Amplicon*>() ;
		}

		/*
		 * The index destructor
		 */
		AmpliconIndex::~AmpliconIndex(void) {
			if( _idx_f_a != NULL ) delete _idx_f_a ;
			if( _idx_r_a != NULL ) delete _idx_r_a ;
			if( _idx_f_b != NULL ) delete _idx_f_b ;
			if( _idx_r_b != NULL ) delete _idx_r_b ;
			//for( vector<Amplicon*>::iterator it=_amplicons.begin(); it!=_amplicons.end(); ++it ) {
			//	if( *it != NULL ) delete *it ;
			//}
			//_amplicons.clear() ;
		}


		long AmpliconIndex::dbsize( ) const {
			long rval = 0 ; 
			for( vector<Amplicon*>::const_iterator it=_amplicons.cbegin(); it!=_amplicons.cend(); ++it ) {
				Amplicon* a = *it ;
				rval += a->width() ;
			}
			return rval ;
		}

		unsigned int AmpliconIndex::n_amplicons() const { 
			return (unsigned int) _amplicons.size() ;
		}

		//
		// Add
		//

		/*
		 * Adds an amplicon to the index
		 */
		bool AmpliconIndex::add( Amplicon* a ) {
			bool rval = false ;
			if( a != NULL ) {
				_amplicons.push_back( a ) ;
			}
			return rval ;
		}

		/*
		 * Builds the amplicon index
		 */
		unsigned int AmpliconIndex::build( int ks ) {

			// don't do anything if the index has already been build
			if( _idx_f_a != NULL || _idx_r_a != NULL || _idx_f_b != NULL || _idx_r_b != NULL) {
				return 0 ;
			}

			// set the keysize
			_keysize = ks ;
			
			// build the roots of the index tree
			_idx_f_a = new _DNANode<Amplicon*>( ) ;
			_idx_r_a = new _DNANode<Amplicon*>( ) ;
			_idx_f_b = new _DNANode<Amplicon*>( ) ;
			_idx_r_b = new _DNANode<Amplicon*>( ) ;

			// sort the amplicons prior to assignment 
			sort( _amplicons.begin(), _amplicons.end(), cmp_lt_amplicon_p ) ;

			// iterate over the amplicon pointer vector to add them to the various trees
			for( vector< Amplicon* >::iterator it=_amplicons.begin(); it!=_amplicons.end(); ++it ) {
				
				// get the pointer to the current amplicon
				Amplicon* a = *it ;

				//
				// Key determination procedure
				//
				string k_f_a  = "" ;
				string k_r_a  = "" ;
				string k_f_b  = "" ;
				string k_r_b  = "" ;
				string seq    = a->sequence() ;
				

				// get the key sequences
				if( seq.size() >= (unsigned int) _keysize ) {
					k_f_a  = seq.substr( 0, _keysize ) ;

					// get the reverse complement
					stringstream x ;
					string tmp = seq.substr( seq.size() - _keysize, _keysize ) ;					
					for( string::reverse_iterator rit=tmp.rbegin(); rit!=tmp.rend(); ++rit ) {
						x << utils::complement_base( *rit ) ; 						
					}
					k_r_a = x.str() ;
				}

				// get the second key pair
				if( seq.size() >= (unsigned int) _keysize * 2 ) {
					k_f_b = seq.substr( _keysize, _keysize ) ;

					// reverse complement
					stringstream x ;
					string tmp = seq.substr( seq.size() - (2 *_keysize), _keysize ) ;
					for( string::reverse_iterator rit=tmp.rbegin(); rit!=tmp.rend(); ++rit ) {
						x << utils::complement_base( *rit ) ; 
					}
					k_r_b = x.str() ;
				}
				
				// switch the keys if the amplicon was designed on the opposite strand
				if( !a->forward() ) {
					string tmp = k_f_a ;
					k_f_a = k_r_a ; 
					k_r_a = tmp ;

					tmp   = k_f_b ;
					k_f_b = k_r_b ; 
					k_r_b = tmp ;
				}

				//
				// Add the amplicon pointer to the trees
				//
				// printf("Build keys are %s %s\n", k_f_a.c_str(), k_r_a.c_str() ) ;
				// cout << a->str() << "\t" << k_f_a << "\t" << k_r_a << "\t" << k_f_b << "\t" << k_r_b << endl ;

				// add the amplicon pointer to the various trees
				if( (int)k_f_a.size() == _keysize ) addTreeData( _idx_f_a, k_f_a, a ) ;
				if( (int)k_r_a.size() == _keysize ) addTreeData( _idx_r_a, k_r_a, a ) ;
				if( (int)k_f_b.size() == _keysize ) addTreeData( _idx_f_b, k_f_b, a ) ;
				if( (int)k_r_b.size() == _keysize ) addTreeData( _idx_r_b, k_r_b, a ) ;
			}

			// return the number of amplicons processed
			return (unsigned int)_amplicons.size() ; 
		}

		/*
		 * Gets the amplicons corresponding to read sequences f and r
		 */
		vector<Amplicon*> AmpliconIndex::getAmplicons( string f, string r ) const {

			// declare the return value
			vector<Amplicon*> rval = vector<Amplicon*>() ;
			vector<Amplicon*> vf   = getAmpliconsF( f ) ;
			vector<Amplicon*> vr   = getAmpliconsR( r ) ;			
						
			// get the overlap between the 2 sets 
			rval.resize( vf.size() + vr.size() ) ;
			vector<Amplicon*>::iterator it = set_intersection( vf.begin(), vf.end(), 
																vr.begin(), vr.end(), 
																rval.begin(), cmp_lt_amplicon_p ) ;
			rval.resize( it - rval.begin() ) ;

			// printf( "%d\t%d\t%d\n", vf.size(), vr.size(), rval.size() ) ;
			// return the amplicons
			return rval ;
		}

		vector<Amplicon*> AmpliconIndex::getAmpliconsF( string f ) const {

			// declare the return value
			vector<Amplicon*> rval = vector<Amplicon*>() ;
			vector<Amplicon*> va   = vector<Amplicon*>() ;
			vector<Amplicon*> vb   = vector<Amplicon*>() ;
		
				// process the first key
			if( f.size() >= (unsigned int) _keysize ) {
				string fk = f.substr( 0, _keysize ) ;				
				va = getTreeData( _idx_f_a, fk ) ;	
			}
			if( f.size() >= (unsigned int) 2 * _keysize ) {
				string fk = f.substr( _keysize, _keysize ) ;
				vb = getTreeData( _idx_f_b, fk ) ;
			}	

			// get the union of both target lists
			rval.resize( va.size() + vb.size() ) ;
			vector<Amplicon*>::iterator it = set_union( va.begin(), va.end(), 
														vb.begin(), vb.end(), 
														rval.begin(), cmp_lt_amplicon_p ) ;
			rval.resize( it - rval.begin() ) ;

			//
			return rval ;
		}

		vector<Amplicon*> AmpliconIndex::getAmpliconsR( string r ) const {

			// declare the return value
			vector<Amplicon*> rval = vector<Amplicon*>() ;
			vector<Amplicon*> va   = vector<Amplicon*>() ;
			vector<Amplicon*> vb   = vector<Amplicon*>() ;
		
				// process the first key
			if( r.size() >= (unsigned int) _keysize ) {
				string rk = r.substr( 0, _keysize ) ;				
				va = getTreeData( _idx_r_a, rk ) ;	
			}
			if( r.size() >= (unsigned int) 2 * _keysize ) {
				string rk = r.substr( _keysize, _keysize ) ;	
				vb = getTreeData( _idx_r_b, rk ) ;
			}	

			// get the union of both target lists
			rval.resize( va.size() + vb.size() ) ;
			vector<Amplicon*>::iterator it = set_union( va.begin(), va.end(), 
														vb.begin(), vb.end(), 
														rval.begin(), cmp_lt_amplicon_p ) ;
			rval.resize( it - rval.begin() ) ;

			// return the return value
			return rval ;
		}


		/*
		 * Gets the amplicons corresponding to reads f and r
		 */
		vector<Amplicon*> AmpliconIndex::getAmplicons( Read* f, Read* r ) const {
			if( f != NULL && r != NULL )  {
				return getAmplicons( f->sequence(), r->sequence() ) ;
			} else if( f != NULL &&  r == NULL )  {
				return getAmpliconsF( f->sequence() ) ;
			} else {
				return vector<Amplicon*>() ;
			}
		}

		vector<Amplicon*> AmpliconIndex::getAmplicons( pair<Read*, Read*> p ) const {
			if( p.first != NULL && p.second != NULL )  {
				return getAmplicons( p.first->sequence(), p.second->sequence() ) ;
			} else if( p.first != NULL && p.second == NULL )  {
				return getAmpliconsF( p.first->sequence() ) ;
			} else {
				return vector<Amplicon*>() ;
			}
		}
	}
}