/********************** io_helpers.c **********************************/
/* MIMD version 6 */
/* DT 8/97 
     General purpose high level routines, to be used by any application
     that wants them.

  cmn 11/24
    I have cut out the routines that I need for the 
    codes that sums up the hopping expansion from the file
    io_helpers.c. This file required a number of other files
    to be included.

*/

#include "../../include/io_lat.h"
#include <string.h>

/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string precede the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */

int get_f( FILE *fp, int prompt, char *variable_name_string, Real *value ){
    int s;
    char checkname[80];

    if(prompt)  {
    	printf("enter %s ",variable_name_string);
#if PRECISION == 1
    	s=fscanf(fp,"%e",value);
#else
    	s=fscanf(fp,"%le",value);
#endif
    	if(s == 1){
	    printf("%s %g\n",variable_name_string,*value);
	    return(0);
	}
    }
    else  {
#if PRECISION == 1
    	s=fscanf(fp,"%s %e",checkname,value);
#else
    	s=fscanf(fp,"%s %le",checkname,value);
#endif
    	if (s == EOF) return(1);
    	if(s == 2 && strcmp(checkname,variable_name_string) == 0){
	    printf("%s %g\n",variable_name_string,*value);
    	    return(0);
	}
    }
    printf("get_f: ERROR IN INPUT: %s\n",variable_name_string);
    return(1);
}

int get_i(FILE *fp, int prompt, char *variable_name_string, int *value ){
    int s;
    char checkname[80];

    if(prompt)  {
    	printf("enter %s ",variable_name_string);
    	s=fscanf(fp,"%d",value);
    	if (s == 1){
	    printf("%s %d\n",variable_name_string,*value);
	    return(0);
	}
    }
    else  {
    	s=fscanf(fp,"%s%d",checkname,value);
    	if (s == EOF) return(1);
    	if(s == 2 && strcmp(checkname,variable_name_string) == 0){
	    printf("%s %d\n",variable_name_string,*value);
    	    return(0);
	}
    }
    printf("get_i: ERROR IN INPUT: %s\n",variable_name_string);
    return(1);
}

/* Read a single word as a string */

int get_s( FILE *fp, int prompt, char *variable_name_string, char *value ){
    int s;
    char checkname[80];

    if(prompt)  {
    	printf("enter %s ",variable_name_string);
    	s=fscanf(fp,"%s",value);
    	if(s == 1){
	    printf("%s %s\n",variable_name_string,value);
	    return(0);
	}
    }
    else  {
    	s=fscanf(fp,"%s %s",checkname,value);
    	if (s == EOF) return(1);
    	if(s == 2 && strcmp(checkname,variable_name_string) == 0){
	    printf("%s %s\n",variable_name_string,value);
    	    return(0);
	}
    }
    printf("get_s: ERROR IN INPUT: %s\n",variable_name_string);
    return(1);
}

/* get_prompt gets the initial value of prompt */
/* 0 for reading from file, 1 prompts for input from terminal */
/* should be called only by node 0 */
/* return 0 if sucessful, 1 if failure */
int get_prompt(FILE *fp, int *prompt ){
    char initial_prompt[80];

    *prompt = -1;
    printf( "type 0 for no prompts  or 1 for prompts\n");
    fscanf(fp, "%s",initial_prompt);
    if(strcmp(initial_prompt,"prompt") == 0)  {
       fscanf(fp, "%d",prompt);
    }
    else if(strcmp(initial_prompt,"0") == 0) *prompt=0;
    else if(strcmp(initial_prompt,"1") == 0) *prompt=1;

    if( *prompt==0 || *prompt==1 )return(0);
    else{
        printf("get_prompt: ERROR IN INPUT: initial prompt\n");
        return(1);
    }
}
