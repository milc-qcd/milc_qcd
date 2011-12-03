/************ myclock.c *************************/
/* MIMD version 7 */

#include "time.h"
#include "stdio.h"
void myclock(int restart, char* msg)
{	
    double secs;
    static double lastsecs;
    static double initsecs;

    
    secs = clock()/(CLOCKS_PER_SEC/100);
    secs/=100;
    if (restart) initsecs=secs; 
    else
    printf("TIMING %s : %e %e\n",msg,secs-lastsecs,secs-initsecs); 

    lastsecs=secs;
} 
