# Makefile for inline SSE2 math routines for QCD programs
#

MAKEFILE = Makefile
PERL = perl
N2C = ${PERL} ../bin/nasm2c.pl
LIB = ../nas

all:
	for f in ${LIB}/*.s; do \
	  t=`basename $$f .s`.h; \
	  $(N2C) $$f > $$t; \
	done
# Exceptions
	$(N2C) ${LIB}/sse_m_amv_4vec.s aa,bb,cc0,cc1,cc2,cc3 > sse_m_amv_4vec.h
	$(N2C) ${LIB}/sse_m_mv_s_4dir.s aa,bb0,bb1,bb2,bb3,cc > sse_m_mv_s_4dir.h

# Do not delete inline_sse.h!
clean:
	-/bin/rm -f sse_*.h


