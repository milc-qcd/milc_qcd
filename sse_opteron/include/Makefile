# Makefile for inline SSE2 math routines for QCD programs
#

MAKEFILE = Makefile
PERL = perl
N2C = $(PERL) ../bin/nasm2c.pl
HDR = $(PERL) ../bin/headers.pl
NAS = ../nas

all:
	for f in ${NAS}/*.s; do \
	t=`basename $$f .s`.h; \
	  $(N2C) < $$f > $$t; \
	done
	$(HDR) > inline_headers.h
# Exceptions
	$(N2C) ${NAS}/sse_m_amv_4vec.s aa,bb,cc0,cc1,cc2,cc3 > sse_m_amv_4vec.h
	$(N2C) ${NAS}/sse_m_mv_s_4dir.s aa,bb0,bb1,bb2,bb3,cc > sse_m_mv_s_4dir.h
	$(N2C) ${NAS}/sse_sub4vecs.s aa,bb0,bb1,bb2,bb3 > sse_sub4vecs.h
	$(N2C) ${NAS}/sse_s_m_a_vec.s aa,bb,cc,dd > sse_s_m_a_vec.h
	$(N2C) ${NAS}/sse_s_m_a_mat.s aa,bb,cc,dd > sse_s_m_a_mat.h

# Do not delete inline_sse.h!
clean:
	-/bin/rm -f sse_*.h

