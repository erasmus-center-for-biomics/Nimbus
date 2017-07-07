
CC = g++
baseCFLAGS  = -c -g -Wall -O4 -std=c++0x
baseLDFLAGS = -g -L/usr/lib64 -std=c++0x

htsinclude = $${HTS_INCLUDE}
htslib   = $${HTS_LIB}

# find source and targets and set the object files
src = $(wildcard src/*.cpp)
obj = $(patsubst src/%.cpp, build/%.o, $(src))

# what to do 
all: caller

caller: $(obj)
	mkdir -p bin/
	$(CC) $(obj) $(baseLDFLAGS) \
		-pthread \
		$(htslib) -lm -lz \
		-lboost_program_options \
		-o bin/nimbus_call

#
# build recipes
#		
build/%.o: src/%.cpp
	mkdir -p build
	$(CC) $(baseCFLAGS) -Iinclude \
		-I$(htsinclude) \
		src/$*.cpp -o $@
	