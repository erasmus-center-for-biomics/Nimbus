#pragma once

namespace threadutils {

	/**
	 A queue for used in threads
	 **/
	template<class T>
	class TQueue {

	protected:
		std::queue<T> _b ;
		std::mutex _m ;

	public:

		TQueue(void) {
			_b = std::queue<T>() ;
		}

		~TQueue(void){ }

		TQueue( const TQueue& other ) {
			std::lock_guard<std::mutex> guard( other.get_mutex() ) ;
			_b = other._b ;
		}

	private:
		TQueue& operator=( const TQueue& ) ;

	public:
		void push( T val ){ 
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			_b.push( val ) ;
		} 

		T front( ) {
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			return _b.front() ;
		}

		void pop( ) {
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			_b.pop() ;
		}

		bool shift( T& rval ) {
			bool r = false ;
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			if( ! _b.empty() ) {
				rval = _b.front() ;
				_b.pop() ;			
				r = true ;
			}
			return r ;
		}

		bool empty() {
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			bool r = _b.empty() ;
			return r ;
		}

		size_t size() {
			std::lock_guard<std::mutex> guard( get_mutex() ) ;
			return _b.size() ;
		}

		std::mutex& get_mutex( ) {
			return _m ;
		}
	};

}