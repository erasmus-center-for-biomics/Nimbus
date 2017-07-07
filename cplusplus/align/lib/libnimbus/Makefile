#!/bin/env make 

# set the library name
libname=libnimbus.a

# is either c++0x or c++11
cversion=c++0x

# basic parameters
AR = ar
CC = g++
baseCFLAGS = -c -g -Wall -O4 -std=$(cversion)
baseLDFLAGS = -g -L/usr/lib64 -std=$(cversion)

# find source and targets and set the object files
src = $(wildcard src/*.cpp)
obj = $(patsubst src/%.cpp, build/%.o, $(src))

# 
all: $(obj) 
	$(AR) -cvq $(libname) $(obj)
	rm -r build

clean:
	-rm -rf build/*
	-rm $(libname)
	
build/%.o: src/%.cpp
	mkdir -p build
	$(CC) $(baseCFLAGS) -Iinclude src/$*.cpp -o $@
	