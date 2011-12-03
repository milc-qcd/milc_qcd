.SUFFIXES: .texi .html .info .dvi .ps .pdf

DESTPATH = ${HOME}/public_html/milc
DELIVERABLES = milcv7.html milcv7_toc.html milcv7.ps milcv7.pdf

.texi.html::
	texi2html --l2h --iftex $*.texi

.ps.pdf::
	ps2pdf $*.ps

.dvi.ps::
	dvips $*.dvi

.texi.dvi::
	env TEX= texi2dvi $*.texi

milcv7_toc.html: milcv7.html

all: ${DELIVERABLES} Makefile
#	/bin/cp -p ${DELIVERABLES} ${DESTPATH}

