
CC=gcc

VPATH=.:../Types:../xmlWriter

CFLAGS= -Wall -I. -I../Types -I../xmlWriter

LDFLAGS= -L../xmlWriter -lxmlWriter -L../Types -lString

test-xmlOutput: test-xmlOutput.o gaugeConfigMetadata.o
	$(CC) -o $@ test-xmlOutput.o gaugeConfigMetadata.o $(LDFLAGS)

clean:
	/bin/rm -f test-xmlOutput.o gaugeConfigMetadata.o

distclean: clean
	/bin/rm -f test-xmlOutput *~

HEADERS= gaugeConfigMetadata.h xmlWriter.h Types.h

test-xmlOutput.o: test-xmlOutput.c $(HEADERS)

gaugeConfigMetadata.o: gaugeConfigMetadata.c $(HEADERS)


