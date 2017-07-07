#pragma once

#include "stdafx.h"
#include "Read.h"

namespace Nimbus {

	namespace utils {

		/**
		 * Returns the complement of base b. If b did not represent a base return 'N'
		 **/
		char complement_base( char b ) ;

		/**
		 * Splits the string x on the delimiter d
		 **/
		std::vector<std::string> split_string( std::string x, std::string d )  ;

		/** 
		 Determines where the query of the path actually starts  (second > 0)
		 **/
		int QueryStart( std::vector< std::pair<int,int> > path ) ;
		
		/** 
		 * Determines where the query of the path actually ends (second > 0 from the end)
		 **/
		int QueryEnd( std::vector< std::pair<int,int> > path ) ;

		/**
		 * Pack the vector of T in pairs of T and int		 
		 **/
		template <class T>
		std::vector< std::pair<T,int> > runLengthEncodePack( std::vector<T> chars ) {

			// declare the return vector
			std::vector< std::pair<T,int> > rval = std::vector< std::pair<T,int> >() ;

			// values to keep track of current and previous
			bool check = false ;
			T p ;

			// define a pair
			std::pair<T,int> pr = std::pair<T,int>() ;
			pr.second = 0 ;

			// for each element in a vector
			for( typename std::vector<T>::iterator it=chars.begin(); it!=chars.end(); ++it ) {

				// if not the first cycle
				if( check ) {

					// if the same as the previous 
					if( p == *it ) {
						// increase the number of elements
						pr.second += 1 ;
					} else {
						// otherwise add the previous pair to the result and initialize a new path
						rval.push_back( pr ) ;
						pr        = std::pair<T,int>() ;
						pr.first  = *it ;
						pr.second = 1 ;
					}
				// otherwise fill the first pair
				} else {
					pr.first  = *it ;
					pr.second = 1 ;
					check = true ;
				}
				// the previous is set to the current at the end of the loop				
				p = *it ;
			}

			// add the last element if the chars vector was not empty
			if( chars.size() > 0 ) rval.push_back( pr ) ;

			return rval ;
		}

		/**
		 * Converts the packed values to a string
		 * 
		 * fs is the field separator
		 * ns is the number separator
		 **/
		template <class T>
		std::string runLengthEncode( std::vector< std::pair<T,int> > elem, std::string fs, std::string ns ) {
			std::stringstream ss ;		
			bool first = true ;
			for( typename std::vector< std::pair<T,int> >::iterator it=elem.begin(); it!=elem.end(); ++it ) {	
				if( !first ) ss << fs ;
				ss << it->second << ns << it->first ;
				first = false ;
			}
			// return the string value of the stringstream
			return ss.str() ;
		}

		/**
		 * Conveniance function to encode and pack the 
		 * entries in a CIGAR character vector
		 **/
		template <class T>
		std::string runLengthEncode( std::vector<T> chars, std::string fs, std::string ns ) {
			return runLengthEncode( runLengthEncodePack( chars ), fs, ns ) ;
		}


	
	}
}