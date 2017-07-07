#include "stdafx.h"
#include "Utils.h"

namespace Nimbus {

	namespace utils {

		char complement_base( char b ) {
			char rval = 'N' ;

			switch(b) {
			case 'A':
				rval = 'T' ;
				break ;
			case 'T':
				rval = 'A' ;
				break ;
			case 'G':
				rval = 'C' ;
				break ;
			case 'C':
				rval = 'G' ;
				break ;
			case 'a':
				rval = 't' ;
				break ;
			case 't':
				rval = 'a' ;
				break ;
			case 'g':
				rval = 'c' ;
				break ;
			case 'c':
				rval = 'g' ;
				break ;
			default:
				rval = 'N' ;
				break ;
			}

			return rval ;
		}

		std::vector<std::string> split_string( std::string x, std::string d )  {

			//
			std::vector<std::string> rval = std::vector<std::string>() ;

			// fill the rval 
			size_t pos = 0;
			std::string t ;
			while( (pos = x.find(d)) != std::string::npos ) {

				t = x.substr(0, pos) ;
				rval.push_back( t ) ;

				x.erase( 0, pos + d.length() ) ;
			}
			// add the last entry
			rval.push_back( x ) ;

			// return the data
			return rval ;
		}

		int QueryStart( std::vector< std::pair<int,int> > path ) {

			// declare the return value 
			int rval = 0 ;

			// traversing the path, we will see whether the query coorindate exceeds 0. 
			// If it does break as we found the query start. If not try the next position in the path 
			for( std::vector< std::pair<int,int> >::iterator it=path.begin(); it!=path.end(); ++it ) {
				if( it->second > 0 ) {
					break ;
				} else {
					rval++ ;
				}
			}
			return rval ;
		}
		
		
		int QueryEnd( std::vector< std::pair<int,int> > path ) {

			// the return value
			int rval = (int) path.size() - 1 ;

			// the end value of the query
			int endval  = path.at(path.size() - 1).second ;

			// traverse the path from the end to the beginning
			for( std::vector< std::pair<int,int> >::reverse_iterator it=path.rbegin(); it!=path.rend(); ++it ) {

				// if we found that the current query coordinate equals
				// the end coordinate, decrease the counter. If not we have 
				// found the end and should break the loop 
				if( it->second != endval ) {
					break ;
				} else {
					rval-- ;
				}
			}
			return rval ;
		}

	}
}