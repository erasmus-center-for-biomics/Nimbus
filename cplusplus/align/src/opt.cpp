
#include "nimbusheader.h"
#include "opt.h"


namespace commandline {

	bool FileExists( const std::string fn ) {
		bool rval = false ;
		std::ifstream h( fn.c_str() ) ;
		if( h.good() ) {
			h.close();
			rval = true ;
		} else {
			h.close() ;
		}
		return rval ;
	}

}