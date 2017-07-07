#pragma once

#include "stdafx.h"

namespace threadutils {

	/**
	 A class to send signals between threads 
	 **/
	template <class T>
	class Signal {

	protected:
		T _signal ;
		std::mutex _m ;

	public:
		Signal(void) {
		
		}

		Signal( T s ): _signal(s) { 
		
		}

		Signal( const Signal& other ) {
			std::lock_guard<std::mutex> guard( other.get_mutex() ) ;
			_signal = other._signal ;
		}

		~Signal(void) {}

	private:
		Signal& operator=( const Signal& ) ;

	public:
		/*
		 set the signal value
		 */
		void set( T val ){ 
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			_signal = val ;
		} 

		/*
		 get the signal object
		 */
		T get( ) {
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			return _signal ;
		}

		/*
		 get the mutex
		 */
		std::mutex& get_mutex( ) {
			return _m ;
		}

			
	};

}
