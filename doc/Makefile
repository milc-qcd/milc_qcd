.SUFFIXES: .texi .html .info .dvi .ps

DESTPATH = ${HOME}/public_html/milc
DELIVERABLES = milcv6.html milcv6_toc.html milcv6.ps

.texi.html::
	texi2html $*.texi

.dvi.ps::
	dvips $*.dvi

.texi.dvi::
	env TEX= texi2dvi $*.texi

milcv6_toc.html: milcv6.html

all: ${DELIVERABLES} Makefile
#	/bin/cp -p ${DELIVERABLES} ${DESTPATH}

