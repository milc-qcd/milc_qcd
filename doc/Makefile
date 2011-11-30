.SUFFIXES: .texi .html .info .dvi .ps

DESTPATH = ${HOME}/public_html/milc
DELIVERABLES = milcv7.html milcv7_toc.html milcv7.ps

.texi.html::
	texi2html --l2h --iftex $*.texi

.dvi.ps::
	dvips $*.dvi

.texi.dvi::
	env TEX= texi2dvi $*.texi

milcv7_toc.html: milcv7.html

all: ${DELIVERABLES} Makefile
#	/bin/cp -p ${DELIVERABLES} ${DESTPATH}

