#pragma once

namespace Nimbus {

	namespace seed {

		//
		// _DNANode
		//
		template <class T>
		class _DNANode {
			_DNANode* _A ;
			_DNANode* _T ;
			_DNANode* _C ;
			_DNANode* _G ;
			_DNANode* _N ;
			_DNANode* _parent ;
			std::vector<T>* _data ;

		public:
			_DNANode() {
				_A = NULL ;
				_T = NULL ;
				_C = NULL ;
				_G = NULL ;
				_N = NULL ;
				_parent = NULL ;
				_data   = NULL ;
			}

			_DNANode(_DNANode * p ): _parent(p) {
				_A = NULL ;
				_T = NULL ;
				_C = NULL ;
				_G = NULL ;
				_N = NULL ;
				_data = NULL ;
			}

			~_DNANode() {
				if( _A != NULL ) delete _A ;
				if( _T != NULL ) delete _T ;
				if( _C != NULL ) delete _C ;
				if( _G != NULL ) delete _G ;
				if( _N != NULL ) delete _N ;
				if( _data != NULL ) delete _data ;
			}

			/* Adds data to the DNANode */
			void addData( T l ) {
				if( _data == NULL ) {
					_data = new std::vector<T>() ;
				}
				_data->push_back( l ) ;
			}

			/* Adds a child at base b to the node */
			bool addChild( char b, _DNANode* d ) {
				bool rval = true ;
				switch(b) {
				case 'A':
					_A = d ;
					break ;
				case 'T':
					_T = d ;
					break ;
				case 'C':
					_C = d ;
					break ;
				case 'G':
					_G = d ;
					break ;
				case 'N':
					_N = d ;
					break ;
				default:
					rval = false ;
					break ;
				}
				return rval ;
			}

			/* Get the child at base b */
			_DNANode* getChild( char b ) const {
				_DNANode* rval = NULL ;
				switch(b) {
				case 'A':
					rval = _A ;
					break ;
				case 'T':
					rval = _T ;
					break ;
				case 'C':
					rval = _C ;
					break ;
				case 'G':
					rval = _G ;
					break ;
				case 'N':
					rval = _N ;
					break ;
				}
				return rval ;
			}

			/* get the data loaded to this _DNANode object */
			std::vector<T> data() const {
				if( _data != NULL ) {
					// printf("data is not NULL: contains %d elements\n", _data->size() ) ;
					return std::vector<T>(_data->begin(), _data->end() ) ;
				} else {
					// printf("data is NULL\n") ;
					return std::vector<T>() ;
				}
			}
		} ;

		
		/*
		 Adds a value to the _DNANode tree at branch key
		 */
		template <class T>
		void addTreeData( _DNANode<T>* base, std::string key, T l ) {
			_DNANode<T>* r = base ;

			// printf( "start key: " ) ;
			// foreach character in the string
			for( std::string::iterator it=key.begin(); it!=key.end(); ++it ) {
				_DNANode<T>* n = NULL ;
				n = r->getChild( *it ) ;
				if( n == NULL ) {
					n = new _DNANode<T>( r ) ;
					r->addChild( *it, n ) ;
					// printf( "%c-", *it) ;
				}
				r = n ;
			}

			// after traversing the string add the index to the end
			if( r != NULL ) {
				// printf("adding") ;
				r->addData( l ) ;
			}
			// printf( "\n" ) ;
		}

		/*
		 Gets the value from the _DNANode tree at branch key
		 */
		template <class T>
		std::vector<T> getTreeData( _DNANode<T>* base, std::string key ) {
			std::vector<T> rval = std::vector<T>() ;

			// printf( "retrieval start key: " ) ;
			// traverse the string to get an end-point
			_DNANode<T>* r = base ;
			for( std::string::iterator it=key.begin(); it!=key.end(); ++it ) {
				// printf( "%c-", *it) ;
				_DNANode<T>* n = NULL ;
				n = r->getChild( *it ) ;
				if( n == NULL ) {
					// printf("Break") ;
					break ;
				}
				r = n ;
			}
			// printf( "\n" ) ;

			// if we successfully got to an end-point, copy its contents to the return vector 
			// Note that this should never be too big 
			if( r != NULL ) {
				// printf( "r not NULL\n" ) ;
				rval = r->data() ;				
			}

			// return the std::vector<T>
			return rval ;
		}

	}
}