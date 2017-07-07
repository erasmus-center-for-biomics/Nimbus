
# is either c++0x or c++11
cversion=c++0x

#
CC = g++
baseCFLAGS = -c -g -Wall -O4 -std=$(cversion)
baseLDFLAGS = -g -L/usr/lib64 -std=$(cversion)
threadlib= -pthread

# find source and targets and set the object files
src = $(wildcard src/*.cpp)
obj = $(patsubst src/%.cpp, build/%.o, $(src))

# what to do 
all: nimbus

libnimbus:
	$(MAKE) -C lib/libnimbus

libthreadutils:
	$(MAKE) -C lib/libthreadutils

nimbus: libnimbus libthreadutils $(obj) 	
	mkdir -p bin/
	$(CC) $(obj) $(baseLDFLAGS) $(threadlib) \
		-Iinclude \
		-Llib/libnimbus -lnimbus -Ilib/libnimbus/include \
		-Llib/libthreadutils -lthreadutils -Ilib/libthreadutils/include \
		-o bin/nimbus_align
	rm -r build
		
clean:
	$(MAKE) -C lib/libnimbus clean
	$(MAKE) -C lib/libthreadutils clean
	-rm -rf build/*
	-rm -rf bin/nimbus
	
	
build/%.o: src/%.cpp
	mkdir -p build
	$(CC) $(baseCFLAGS) -Iinclude \
		-Ilib/libnimbus/include \
		-Ilib/libthreadutils/include \
		$(threadlib) src/$*.cpp -o $@
		
