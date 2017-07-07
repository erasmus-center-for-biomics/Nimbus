
# is either c++0x or c++11
cversion=c++0x

# basic parameters
CC = g++
baseCFLAGS = -c -g -Wall -O4 -std=$(cversion)
baseLDFLAGS = -g -L/usr/lib64 -std=$(cversion)
threadlib= -pthread

# find source and targets and set the object files
src = $(wildcard src/*.cpp)
obj = $(patsubst src/%.cpp, build/%.o, $(src))

all: adapter_trim

adapter_trim: build/adapter_trim.o
	-mkdir -p bin/
	$(CC) build/adapter_trim.o $(baseLDFLAGS) $(threadlib) \
		-Irwwb \
		-Iinclude \
		-lboost_program_options \
		-o bin/nimbus_trim

clean:
	-rm -rf build/*
	-rm bin/adapter_trim

build/%.o: src/%.cpp
	mkdir -p build
	$(CC) $(baseCFLAGS) \
		-Irwwb \
		-Iinclude \
		$(threadlib) \
		src/$*.cpp -o $@
