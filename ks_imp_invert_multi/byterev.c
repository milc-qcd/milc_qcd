/******************* byterev.c *********************/
/* DT  6/17/92 */

/* reverse byte order in file
   Usage:  byterev < in > out
*/
#include <stdio.h>
#define BUFSIZE 1024	/* number of words to read at once */

main(argc,argv)  int argc; char ** argv; {
    int buf[BUFSIZE];
    register int i,words;
    register int old,new;

    while( (words = read(0,buf,BUFSIZE*4)/4) != 0 ){

	for(i=0;i<words;i++){
	    old = buf[i];
	    new = old >> 24 & 0x000000ff;
	    new |= old >> 8 & 0x0000ff00;
	    new |= old << 8 & 0x00ff0000;
	    new |= old << 24 & 0xff000000;
	    buf[i] = new;
	}

	write(1,buf,words*4);
    }
}

