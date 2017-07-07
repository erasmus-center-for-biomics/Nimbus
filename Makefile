#
# Compile the C++ parts of Nimbus
#
cversion=c++0x

all: trim align call
	mkdir -p bin
	mv cplusplus/trim/bin/nimbus_trim bin/
	mv cplusplus/align/bin/nimbus_align bin/
	mv cplusplus/call/bin/nimbus_call bin/

trim:
	$(MAKE) -C cplusplus/trim

align:
	$(MAKE) -C cplusplus/align

call:
	$(MAKE) -C cplusplus/call
