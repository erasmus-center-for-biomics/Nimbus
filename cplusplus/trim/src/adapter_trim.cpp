
// STL
#include <cstddef>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

// BOOST libraries
#include <boost/program_options.hpp>

// own libraries
#include <sequence_matcher.hpp>
#include <rwwb/sequtils/types.hpp>
#include <rwwb/sequtils/fastq.hpp>
#include <rwwb/sequtils/fasta.hpp>

//
// Processes the reads in the input stream and writes them to the output stream
//
int process_reads(std::istream& hin, std::ostream& hout, const std::vector< Biomics::SequenceMatcher<rwwb::sequtils::base_t> >& adapters, std::size_t buffersize, std::size_t min_bases_left) {
    
    // prepare the read buffer
    rwwb::sequtils::read default_read = rwwb::sequtils::read() ;    
    default_read.uid = 0 ;
    default_read.name = "" ;
    default_read.sequence = std::vector<rwwb::sequtils::base_t>(0,0) ;
    default_read.quality = "" ;    
    std::vector<rwwb::sequtils::read> buffer = std::vector<rwwb::sequtils::read>(buffersize, default_read) ;     
    
    std::size_t reads_loaded = 0 ;
    std::size_t bases_left = 0 ; 
    
    // read entries from the input stream
    while(true) {
        //
        std::size_t obtained = rwwb::sequtils::reads_from_fastq(hin, reads_loaded, buffer) ;
        
        // process the reads
        for(std::size_t i=0; i<obtained; ++i){
            bases_left = buffer[i].sequence.size() ;
            
            // match the adapters
            for(std::size_t k=0; k<adapters.size(); ++k) {            
                std::size_t tmp = adapters[k].match(buffer[i].sequence) ;
                if(tmp < bases_left) {
                    bases_left = tmp ;
                }
            }
            
            // trim if necessary
            if(bases_left < buffer[i].sequence.size()) {
                
                // If too few bases are left 
                // 
                if(bases_left <= min_bases_left) {
                    // replace the sequence by N's
                    buffer[i].sequence = std::vector<rwwb::sequtils::base_t>(buffer[i].sequence.size(), -1) ;
                    buffer[i].quality = std::string(buffer[i].sequence.size(), '!') ;
                } else {
                    // otherwise truncate the read                    
                    buffer[i].sequence = std::vector<rwwb::sequtils::base_t>(buffer[i].sequence.begin(), buffer[i].sequence.begin() + bases_left) ;
                    buffer[i].quality = std::string(buffer[i].quality.begin(), buffer[i].quality.begin() + bases_left) ;
                }
            }
        }
        
        // write the buffer to the output
        for(std::size_t i=0; i<obtained; ++i){
            std::vector<char> tmp(buffer[i].sequence.size(), 'N') ;
            std::transform(buffer[i].sequence.begin(), buffer[i].sequence.end(), tmp.begin(), rwwb::sequtils::base_to_char) ;            
            hout << "@" << buffer[i].name << std::endl 
                << std::string(tmp.begin(), tmp.end()) << std::endl 
                << "+" << std::endl 
                << buffer[i].quality << std::endl ; 
        }
        
        // if we didn't get the expected number of reads stop the analysis
        if(obtained != buffersize) 
            break ;
    }   
    return 0 ;
}

//
// Add adapters from file_adapters to the adapters vector
//
void adapter_helper(std::vector< std::vector<rwwb::sequtils::base_t> >& sequences, std::string file_adapters ) {
    std::string label ;
    std::vector<rwwb::sequtils::base_t> seq ;
    auto parser = rwwb::sequtils::fasta() ;                
    if(file_adapters == "-") {                        
        while( parser(std::cin, label, seq) ){                
            sequences.push_back(seq) ;
        }
    } else {
        std::ifstream fadapter(file_adapters.c_str(), std::ifstream::in) ;
        while( parser( fadapter, label, seq) ){
            sequences.push_back(seq) ;
        }
        fadapter.close() ;
    }
}


//
// The main program loop
//
int main(int argc, char** argv) {
    
    //
    int return_code = 0 ;
    std::ifstream fin ;
    std::ofstream fout ;
    bool verbose = false ;   
       
    // the input parameters 
    std::string file_input = "-" ;
    std::string file_output = "-" ;
    std::string file_adapters = "" ;
    std::vector<std::string> adapter_sequences(0) ;
    std::size_t maximum_mismatches = 2 ;
    std::size_t minimum_matches = 1 ;
    std::size_t minimum_bases_remaining = 25 ;    
    std::size_t buffer_size = 1000 ;
        
    // parse the command-line parameters
    boost::program_options::options_description desc("Allowed options") ;
    desc.add_options()
        ("help", "Produce the help message")
        ("verbose", "Provide verbose output")
        ("input,i", boost::program_options::value<std::string>(&file_input), "The input FastQ file")
        ("output,o", boost::program_options::value<std::string>(&file_output), "The output FastQ file")
        ("adapter,a", boost::program_options::value< std::vector<std::string> >(&adapter_sequences), "The sequence of an  adapter to trim")
        ("adapter-file,f", boost::program_options::value<std::string>(&file_adapters), "A fasta file with adapter sequences")
        ("maximum-mismatches", boost::program_options::value<std::size_t>(&maximum_mismatches), "The maximum number of mismatches (default 2)")
        ("minimum-matches", boost::program_options::value<std::size_t>(&minimum_matches), "The minimum number of matching bases (default 1)")
        ("minimum-bases-remaining", boost::program_options::value<std::size_t>(&minimum_bases_remaining), "The minimum number of bases in the reads remaining (default 25)")
        ("buffer-size", boost::program_options::value<std::size_t>(&buffer_size), "The size of the read buffer to process (default 1000)") ;    
    boost::program_options::positional_options_description p ;
	boost::program_options::variables_map vm ;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm) ;
	boost::program_options::notify(vm) ;
    
    // check the command-line arguments
    if( vm.count("help")) {
		std::cerr << "Usage" << std::endl ; 
		std::cerr << desc << std::endl; 
		return 0 ;
	}
    if(vm.count("verbose")){
        verbose = true ;
        std::cerr << "Parameters" << std::endl ;
        std::cerr << ".. Input file: " << file_input << std::endl ; 
        std::cerr << ".. Output file: " << file_output << std::endl ;
        std::cerr << ".. Adapter file: " << file_adapters << std::endl ;
        std::cerr << ".. Adapter sequences on commandline: " << std::endl ;
        for(std::size_t i=0; i<adapter_sequences.size(); ++i){
            std::cerr << ".... " << adapter_sequences[i] << std::endl ;
        }
        std::cerr << ".. Maximum number of mismatches: " << maximum_mismatches << std::endl ;
        std::cerr << ".. Minimum number of matches: " << minimum_matches << std::endl ;
        std::cerr << ".. Minimum read size remaining: " << minimum_bases_remaining << std::endl ;    
        std::cerr << ".. Buffer size: " << buffer_size << std::endl ;
        std::cerr << std::endl ;
    }    
    if(file_input == "-" && file_adapters == "-"){
        std::cerr << "Cannot obtain information from 2 input streams" << std::endl << std::endl ;
        std::cerr << "Usage" << std::endl ; 
		std::cerr << desc << std::endl; 
		return 101 ;
    }
        
    // assign the in and output files
    if(file_input != "-"){
        fin.open(file_input.c_str(), std::ifstream::in) ;
    }
    if(file_output != "-"){
        fout.open(file_output, std::ifstream::out) ;
    }
        
    // initialize the adapters
    std::vector<Biomics::SequenceMatcher<rwwb::sequtils::base_t> > adapters ;
    for(std::size_t i=0; i<adapter_sequences.size(); ++i){
        adapters.push_back(Biomics::SequenceMatcher<rwwb::sequtils::base_t>(rwwb::sequtils::string_to_base(adapter_sequences[i]), maximum_mismatches, minimum_matches)) ;    
    }
                
    // adapters should be loaded from a file
    if(file_adapters != ""){        
        std::vector<std::vector<rwwb::sequtils::base_t> > base_t_sequences ; 
        adapter_helper(base_t_sequences, file_adapters) ;
        for(std::size_t i=0; i<base_t_sequences.size(); ++i){
            adapters.push_back(Biomics::SequenceMatcher<rwwb::sequtils::base_t>(base_t_sequences[i], maximum_mismatches, minimum_matches)) ;    
        }
    }  
        
    if(adapters.size() == 0){
        std::cerr << "No adapter sequences added" << std::endl << std::endl ;
        std::cerr << "Usage" << std::endl ; 
		std::cerr << desc << std::endl; 
		return 104 ;
    }
    
    // add the input streams 
    std::istream& hin = fin.is_open() ? fin : std::cin ;
    std::ostream& hout = fout.is_open() ? fout : std::cout ;
        
    // process the reads from standard in 
    return_code = process_reads(hin, hout, adapters, buffer_size, minimum_bases_remaining) ;
          
    // cleanup the opened files
    if(file_input != "-"){
        fin.close() ;
    }
    if(file_output != "-"){
        fout.close() ;
    }
          
    // return the default return code
    return return_code ;
}