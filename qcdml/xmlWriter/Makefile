VPATH= .:../Types

CC=gcc

CFLAGS= -Wall -I . -I ../Types

CXX= g++
CXXFLAGS= $(CFLAGS)

HEADERS= Array.h Types.h cast.h xmlWriter.h outStream.h

LDFLAGS= -L../Types -lArray -lString -L. -lxmlWriter

libs: libxmlWriter.a

libxmlWriter.a: libxmlWriter.a(xmlWriter.o) libxmlWriter.a(outStream.o)

test-xmlWriter: test-xmlWriter.o libxmlWriter.a
	$(CC) -o $@ test-xmlWriter.o $(LDFLAGS)

test-xmlWriter-plusplus: test-xmlWriter-plusplus.o libxmlWriter.a
	$(CXX) -o $@ test-xmlWriter-plusplus.o $(LDFLAGS)

clean:
	/bin/rm -f test-xmlWriter.o test-xmlWriter-plusplus.o \
		xmlWriter.o outStream.o

distclean: clean
	/bin/rm -f test-xmlWriter test-xmlWriter-plusplus libxmlWriter.a *~

test-xmlWriter.o: test-xmlWriter.c $(HEADERS)

test-xmlWriter-plusplus.o: test-xmlWriter-plusplus.cpp $(HEADERS)

libxmlWriter.a(xmlWriter.o): xmlWriter.c $(HEADERS)

libxmlWriter.a(outStream.o): outStream.c outStream.h cast.h
