CC= gcc

CXX= g++

CFLAGS= -Wall -I .

CXXFLAGS= $(CFLAGS)

ALDFLAGS= -L. -lArray

SLDFLAGS= -L. -lString

AHEADERS= Array.h Types.h cast.h

SHEADERS= String.h

libs: libArray.a libString.a

test: libs test-Array test-String

libArray.a: libArray.a(Array.o)

libString.a: libString.a(String.o)

test-Array: test-Array.o libArray.a(Array.o)
	$(CC) -o $@ test-Array.o $(ALDFLAGS)

test-String: test-String.o libString.a(String.o)
	$(CC) -o $@ test-String.o $(SLDFLAGS)

test-Array-plusplus: test-Array-plusplus.o libArray.a(Array.o)
	$(CXX) -o $@ test-Array-plusplus.o $(ALDFLAGS)

clean:
	/bin/rm -f String.o test-String.o Array.o test-Array.o test-Array-plusplus.o

distclean: clean
	/bin/rm -f libArray.a libString.a test-String test-Array test-Array-plusplus *~

libArray.a(Array.o): Array.c $(AHEADERS)

libString.a(String.o): String.c $(SHEADERS)

test-Array.o: test-Array.c $(AHEADERS)

test-String.o: test-String.c $(SHEADERS)

test-Array-plusplus.o: test-Array-plusplus.cpp $(AHEADERS)


