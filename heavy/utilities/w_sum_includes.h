/****************** w_heavy_includes.h  ******************************/
/*
 *  Include files for the clover_invert application */

/* Include files */
#include "../../include/config.h"   /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/complex.h"
#include "../../include/su3.h"
#include "../../include/macros.h"
#include "../../include/generic.h"
#include "lattice_k.h"
#include "../../include/io_lat.h"

/*************************************************** 
   prototypes for functions in high level code 
   ***************************************************/

void sum(Real kappa_h, Real kappa_c, int nhop, int file_nhop,
	   FILE * fp_m_in, int fb_m_in, FILE * fp_k_out, int writeflag) ;


void sum_light(int file_nhop, FILE *fp_m_in, int fb_m_in,
	       FILE *fp_k_out) ;


FILE *(r_ascii_m_i(char *filenam,int i1,int *file_hops))  ;

void r_ascii_m(FILE *fp,int spin,int color,int iter,double *prop[]) ;

void r_ascii_m_f(FILE *fp, char *filenam) ;

int r_binary_m_i(char *filenam,int i1,int *file_hops)  ;
void r_binary_m(int fp, int spin, int color, int iter,double *prop[]) ;
void r_binary_m_f(int fp,char *filenam) ;

FILE *(w_ascii_m_i(char *filenam,int i1)) ;
FILE *(a_ascii_m_i(char *filenam,int i1)) ;
void w_ascii_m(FILE *fp,int spin, int color,int iter,double *prop[] ) ;
void w_ascii_m_f(FILE *fp, char *filenam) ;

int w_binary_m_i(char *filenam, int i1) ;
int a_binary_m_i(char *filenam,int i1) ;
void w_binary_m(int fp, int spin, int color, int iter, double *prop[]) ;
void w_binary_m_f(int fp, char *filenam) ;

FILE *(w_ascii_k_i(char *filenam, int i1)) ;
void w_ascii_k(FILE * fp, int iter, double *prop[]) ;
void w_ascii_k_f(FILE * fp, char *filenam) ;


int setup_k() ;
int initial_set() ;
int readin(int prompt) ;

void initialize_machine(int argc, char **argv);
char * machine_type();
int mynode();
int numnodes();
double dclock();
void time_stamp(char *msg);
void terminate(int status);

/* Directions, and a macro to give the opposite direction */
/*  These must go from 0 to 7 because they will be used to index an
    array. */
/* Also define NDIRS = number of directions */
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

/********** end of the function prototypes ***/


/* end of w_sum_includes.h  */




