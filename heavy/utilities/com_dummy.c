/* Dummy version of com routines used by su3_sum */

#include <time.h>
#include <stdio.h>

/* Machine initialization */
void initialize_machine(int argc, char **argv){}

/* Tell what kind of machine we are on */
static char name[]="Scalar processor";
char * machine_type(){
    return(name);
}

/* Return my node number */
int mynode(){ return(0); }

/* Return number of nodes */
int numnodes(){ return(1); }

/* Print time stamp */
void time_stamp(char *msg){
  time_t time_stamp;
  
  if(mynode()==0){
    time(&time_stamp);
    printf("%s: %s\n",msg,ctime(&time_stamp));
  }
}

/* Double precision time */
/* This one wraps around after 36 minutes!! It gives the cpu time,
   not the wall clock time */
double dclock(){
long fine;
    fine = clock();
    return( ((double)fine)/1000000.0 );
}

/* version of exit for multinode processes -- kill all nodes */
void terminate(int status) {
  time_stamp("termination");
  printf("Termination: status = %d\n",status);
  fflush(stdout);fflush(stderr);
  exit(status);
}
