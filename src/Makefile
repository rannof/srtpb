CC = gcc
CCFLAGS = -O2 -Wall
LIBFLAGS = -lm -lrt
LIBS = -L../ringserver/libmseed -L../ringserver/dalitool/libdali
INCLUDS = -I../ringserver/libmseed -I../ringserver/dalitool/libdali
LIBMSEED = -lmseed
LIBDALI  = -ldali
srtpbX = srtpb
sysswig := $(shell which swig)

all: srtpb msrtslice

srtpb:
ifeq (swig,$(findstring swig,$(sysswig)))
	swig -python srtpblib.i
	$(CC) $(CCFLAGS) $(LIBS) $(INCLUDS) -c srtpblib.c $(LIBFLAGS) $(LIBMSEED) $(LIBDALI) srtpblib_wrap.c -I/usr/include/python2.7/ -I/usr/lib64/python2.7/config -fPIC
	$(CC) $(CCFLAGS) -shared $(LIBS) $(INCLUDS) ../ringserver/dalitool/libdali/libdali.a ../ringserver/dalitool/libmseed/libmseed.a srtpblib.o srtpblib_wrap.o -o _srtpblibmodule.so $(LIBFLAGS) $(LIBMSEED) $(LIBDALI)
endif

msrtslice:
	$(CC) $(CCFLAGS) $(LIBS) $(INCLUDS) ../ringserver/libmseed/libmseed.a ../ringserver/dalitool/libdali/libdali.a msrtslice.c -o msrtslice -lm  $(LIBMSEED) $(LIBDALI)


.IGNORE: clean
clean:
	rm -r *.o *_wrap.c *.so *.pyc srtpblib.py msrtslice

  
