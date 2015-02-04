CC = gcc

all: ring srtpb mkbin
	@echo Done.

ring:
	\make -w -C ringserver

srtpb:
	\make -w -C src all
.IGNORE: clean
clean:
	rm ringserver/pcre/Makefile
	rm ringserver/mxml/Makefile	
	\make -w -C ringserver clean
	\make -w -C src clean
	@rm -r bin
.IGNORE: deepcclean
deepclean:
	rm -r ringserver
	\make -w -C src clean
	@rm -r bin
.IGNORE: mkbin
mkbin:
	@mkdir -p bin/
	@mkdir -p bin/ring
	@mkdir -p bin/segments
	@cp src/*.py src/*.so ringserver/ringserver src/ring.conf ringserver/dalitool/dalitool ringserver/slinktool/slinktool bin
