#ifndef sequence_matcher_hpp
#define sequence_matcher_hpp

// STL
#include <cstddef>
#include <vector>
#include <functional>

namespace Biomics {

    //
    // An iterative matching algorithm
    //
    // 

    /*
    template<typename T>
    std::size_t vector_matcher(const std::vector<T>& subject, const std::vector<T>& query, std::size_t mismatch_threshold, std::size_t minimum_matches) {
        
        //
        std::size_t match_location = subject.size() ;                                     
        std::size_t mismatches = 0 ;            
        std::size_t matches = 0 ;
        bool ismatch = false ;
        
        for(std::size_t i=0; i<subject.size(); ++i){
            if(subject[i] == query[0]) {
                
                // determine where to match
                ismatch = true ;                    
                mismatches = 0 ;
                matches = 0 ;                     
                for(std::size_t j=i; j<subject.size(); ++j) {
                    
                    if(j-i >= query.size())                         
                        break ;                        
                    
                    if(subject[j] != query[j-i]) {
                        mismatches += 1 ;
                        if(mismatches >= mismatch_threshold){
                            ismatch = false ;
                            break ;
                        }
                    } else {
                        matches += 1 ; 
                    }
                }
                
                // check whether the sequence matches 
                if(ismatch && matches > minimum_matches){
                    match_location = i ;
                    break ;                        
                }                    
            }
        }
        
        //
        return match_location ;
    }
    */
    
    template<typename T>
    inline bool default_comparator(T a, T b){
        return a == b ? true : false ;
    }

    template<typename T>
    std::size_t vector_matcher( 
        bool (*compare)(T, T) , const std::vector<T>& subject, const std::vector<T>& query, 
        std::size_t mismatch_threshold, std::size_t minimum_matches, std::size_t first_base) {
        
        //
        std::size_t match_location = subject.size() ;                                     
        std::size_t mismatches = 0 ;            
        std::size_t matches = 0 ;
        bool ismatch = false ;
        
        // 
        if(subject.size() <= first_base){
            return match_location ;
        }

        for(std::size_t i=first_base; i<subject.size(); ++i){
            if(compare(subject[i], query[0])) {
                
                // determine where to match
                ismatch = true ;                    
                mismatches = 0 ;
                matches = 0 ;                     
                for(std::size_t j=i; j<subject.size(); ++j) {
                    
                    if(j-i >= query.size())
                        break ;                        
                    
                    if(!compare(subject[j], query[j-i])) {
                        mismatches += 1 ;
                        if(mismatches >= mismatch_threshold){
                            ismatch = false ;
                            break ;
                        }
                    } else {
                        matches += 1 ; 
                    }
                }
                
                // check whether the sequence matches 
                if(ismatch && matches >= minimum_matches){
                    match_location = i ;
                    break ;                        
                }                    
            }
        }
        
        //
        return match_location ;
    }

    template<typename T>
    std::size_t vector_matcher(const std::vector<T>& subject, const std::vector<T>& query, std::size_t mismatch_threshold, std::size_t minimum_matches){
        return vector_matcher<T>(&default_comparator, subject, query, mismatch_threshold, minimum_matches) ;
    }
    

    template<typename T>
    class SequenceMatcher {
                    
        std::vector<T> sequence ;
        std::size_t maximum_mismatches ;
        std::size_t minimum_matches ;
        std::size_t first_base ;
        bool (*compare)(T, T) ;
      public:
        // Primary constructor
        //
        // 
        SequenceMatcher(std::vector<T> seq):sequence(seq), maximum_mismatches(1), minimum_matches(1), first_base(0), compare(&default_comparator) {}
        
        // Secondary constructor 
        //
        //
        SequenceMatcher(std::vector<T> seq, std::size_t max_mm, std::size_t min_m, std::size_t fb):sequence(seq), maximum_mismatches(max_mm), minimum_matches(min_m), first_base(fb), compare(&default_comparator){}

        // Tertiary constructor
        //
        // (with function assignment)
        SequenceMatcher(std::vector<T> seq, std::size_t max_mm, std::size_t min_m, std::size_t fb, bool (*c)(T, T)):sequence(seq), maximum_mismatches(max_mm), minimum_matches(min_m), first_base(fb), compare(c){}

        // Quaternary constructor
        //
        // (with function assignment)
        SequenceMatcher(std::vector<T> seq, std::size_t max_mm, std::size_t min_m,  bool (*c)(T, T)):sequence(seq), maximum_mismatches(max_mm), minimum_matches(min_m), first_base(0), compare(c){}

        // Quinternary? constructor 
        //
        //
        SequenceMatcher(std::vector<T> seq, std::size_t max_mm, std::size_t min_m):sequence(seq), maximum_mismatches(max_mm), minimum_matches(min_m), first_base(0), compare(&default_comparator){}

        // Determines the point in the vector after which the sequence ought to be trimmed. 
        //
        //           
        std::size_t match(const std::vector<T>& s) const {
            return vector_matcher<T>(compare, s, sequence, maximum_mismatches, minimum_matches, first_base) ;
        }
                            
    } ;
        
} ;
#endif