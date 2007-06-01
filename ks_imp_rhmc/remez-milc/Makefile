CC = gcc
CXX = g++
LDFLAGS = -L$(HOME)/gmp/lib -lmpfr -lgmp -Wl,--rpath,/usr/local/lib64

INCLIST = -I$(HOME)/gmp/include

OBJS = alg_remez.o poly4.o main.o poly4_sqrt.o

CFLAGS=  -O3
CXXFLAGS=  -O3
DFLAGS +=  -DGMP

#  targets

test:	main.o alg_remez.o
	${CXX} alg_remez.o main.o ${LDFLAGS}  -o test

poly4:	poly4.o alg_remez.o
	${CXX} alg_remez.o poly4.o ${LDFLAGS}  -o poly4

poly4_sqrt:	poly4_sqrt.o alg_remez.o
	${CXX} alg_remez.o poly4_sqrt.o ${LDFLAGS}  -o poly4_sqrt

.SUFFIXES:
.SUFFIXES:  .o .C .c

$(BIN):  $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(BIN)

${OBJS}: Makefile

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<

clean:
	rm -f *.dat *.o  *~ *\#* *.stackdump $(BIN)
