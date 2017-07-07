#pragma once


// STL headers
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>	// std::pair
#include <tuple>

namespace nimbus {
	

	//
	// A type where the first value points to an array of 
	// strings, and the second is an index
	//
	// The first value is used for sorting and such, the second 
	// to point to data. The vector is fixed with based on runtime 
	// parameters
	//	

	typedef std::pair< 
		std::vector<std::string>, 
		std::size_t > aspect ;

	int compare( 
		std::vector<std::string> a,
		std::vector<std::string> b ) ;

	int compare( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) ;

	bool lessthan( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) ;

	bool morethan( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) ;

	bool equals( 
		std::pair< std::vector<std::string>,std::size_t > a, 
		std::pair< std::vector<std::string>,std::size_t > b ) ;

	//
	// A tuple to hold the variant statistics
	// 
	// #
	// #
	// #
	//

	//
	//
	//
	//
	//
	
	typedef struct __allele_obj__ {		
		std::string sequence ;			// found sequence
		long n ;						// number of reads
		long qual ; 					// the cummulative quality values
		std::vector<std::string> info ;	// a vector of strings describing the variant
	} allele_obj ;
	
	
	
}