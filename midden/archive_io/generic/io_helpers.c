/********************** io_helpers.c **********************************/
/* MIMD  version 5 */
/* DT 8/97 
     General purpose high level routines, to be used by any application
     that wants them.
*/

#include "generic_includes.h"
#include <io_lat.h>


/* save a lattice in any of the formats:
    SAVE_ASCII, SAVE_SERIAL, SAVE_OLD_BINARY, SAVE_PARALLEL, SAVE_CHECKPOINT
    c1 and c2 are couplings, used in the old (versions 1-4) lattice formats.
*/
gauge_file *save_lattice( int flag, char *filename, Real c1, Real c2 ){
    double dtime;
    gauge_file *gf;
    Real ssplaq, stplaq;
    void plaquette(Real *ss_plaq,Real *st_plaq);

    dtime = -dclock();
    switch( flag ){
        case FORGET:
            gf = NULL;
            break;
	case SAVE_ASCII:
	    gf = save_ascii(filename);
	    break;
	case SAVE_SERIAL:
	    gf = save_serial(filename);
	    break;
	case SAVE_OLD_BINARY:
	    gf = save_old_binary(filename,c1,c2);
	    break;
	case SAVE_PARALLEL:
	    gf = save_parallel(filename);
	    break;
	case SAVE_CHECKPOINT:
	    gf = save_checkpoint(filename);
	    break;
	case SAVE_ARCHIVE:
	    gf = save_archive(filename);
	    break;
	default:
	    printf("save_lattice: ERROR: unknown type for saving lattice\n");
	    terminate(1);
    }
    dtime += dclock();
    if(flag != FORGET)
      node0_printf("Time to save = %e\n",dtime);
    plaquette(&ssplaq,&stplaq);
    node0_printf("CHECK PLAQ: %e %e\n",ssplaq,stplaq);
    return gf;
}

/* reload a lattice in any of the formats, or cold lattice, or keep
     current lattice:
    FRESH, CONTINUE,
    RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL
    c1 and c2 are couplings, used in the old (versions 1-4) binary formats.
*/
gauge_file *reload_lattice( int flag, char *filename){
    double dtime;
    gauge_file *gf;
    Real ssplaq, stplaq;
    Real max_deviation;
    void coldlat();

    dtime = -dclock();
    switch(flag){
	case CONTINUE:	/* do nothing */
            gf = NULL;
	    break;
	case FRESH:	/* cold lattice */
	    coldlat();
            gf = NULL;
	    break;
	case RELOAD_ASCII:	/* read Ascii lattice */
	    gf = restore_ascii(filename);
	    break;
	case RELOAD_SERIAL:	/* read binary lattice serially */
	    gf = restore_serial(filename);
	    break;
	case RELOAD_PARALLEL:	/* read binary lattice in parallel */
	    gf = restore_parallel(filename);
	    break;
	default:
	    if(this_node==0)printf("reload_lattice: Bad startflag %d\n",flag);
	    terminate(1);
    }
    dtime += dclock();
    if(flag != FRESH && flag != CONTINUE)
      node0_printf("Time to reload gauge configuration = %e\n",dtime);
    plaquette(&ssplaq,&stplaq);
    if(this_node==0){
        printf("CHECK PLAQ: %e %e\n",ssplaq,stplaq);fflush(stdout);
    }
    dtime = -dclock();
    max_deviation = check_unitarity();
    g_floatmax(&max_deviation);
    dtime += dclock();
    if(this_node==0)printf("Unitarity checked.  Max deviation %.2e\n",
			   max_deviation); fflush(stdout);
    if(this_node==0)printf("Time to check unitarity = %e\n",dtime);
    return gf;
}

/* find out what kind of starting lattice to use, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_starting_lattice( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) printf(
        "enter 'continue', 'fresh', 'reload_ascii', 'reload_serial', or 'reload_parallel'\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_starting_lattice: ERROR IN INPUT: starting lattice command\n");
        return(1);
    }

    if(strcmp("fresh",savebuf) == 0 ){
       *flag = FRESH;
    }
    else if(strcmp("continue",savebuf) == 0 ) {
        *flag = CONTINUE;
        printf("continuing with current lattice\n");
    }
    else if(strcmp("reload_ascii",savebuf) == 0 ) {
       *flag = RELOAD_ASCII;
    }
    else if(strcmp("reload_serial",savebuf) == 0 ) {
       *flag = RELOAD_SERIAL;
    }
    else if(strcmp("reload_parallel",savebuf) == 0 ) {
       *flag = RELOAD_PARALLEL;
    }
    else{
    	printf("ask_starting_lattice: ERROR IN INPUT: lattice_command is invalid\n"); return(1);
    }

    /*read name of file and load it */
    if( *flag != FRESH && *flag != CONTINUE ){
        if(prompt!=0)printf("enter name of file containing lattice\n");
        status=scanf("%s",filename);
        if(status !=1) {
	    printf("ask_starting_lattice: ERROR IN INPUT: file name read\n"); return(1);
        }
    }
    return(0);
}

/* find out what do to with lattice at end, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_ending_lattice( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) printf(
        "'forget' lattice at end,  'save_ascii', 'save_serial', 'save_old_binary', or 'save_parallel'\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_ending_lattice: ERROR IN INPUT: ending lattice command\n");
        return(1);
    }
    if(strcmp("save_ascii",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_serial",savebuf) == 0 ) {
        *flag=SAVE_SERIAL;
    }
    else if(strcmp("save_old_binary",savebuf) == 0 ) {
        *flag=SAVE_OLD_BINARY;
    }
    else if(strcmp("save_parallel",savebuf) == 0 ) {
      *flag=SAVE_PARALLEL;
    }
    else if(strcmp("save_checkpoint",savebuf) == 0 ) {
        *flag=SAVE_CHECKPOINT;
    }
    else if(strcmp("save_archive",savebuf) == 0 ) {
        *flag=SAVE_ARCHIVE;
    }
    else if(strcmp("forget",savebuf) == 0 ) {
        *flag=FORGET;
    }
    else {
      printf("ask_ending_lattice: ERROR IN INPUT: %s is not a save lattice command\n",savebuf);
      return(1);
    }



    if( *flag != FORGET ){
        if(prompt!=0)printf("enter filename\n");
        status=scanf("%s",filename);
        if(status !=1){
    	    printf("ask_ending_lattice: ERROR IN INPUT: save filename\n");
	    return(1);
        }
/* Don't know how to pass ENSEMBLE_ID and SEQUENCE_NUMBER to save_lattice */ 
	/*
	if(*flag == SAVE_ARCHIVE) {
	*/
        printf("gauge configuration to be saved in %s\n",filename);
    }
    return(0);
}


void coldlat(){
    /* sets link matrices to unit matrices */
    register int i,j,k,dir;
    register site *sit;

    FORALLSITES(i,sit){
	for(dir=XUP;dir<=TUP;dir++){
	    for(j=0; j<3; j++)  {
		for(k=0; k<3; k++)  {
		    if (j != k)  {
		       sit->link[dir].e[j][k] = cmplx(0.0,0.0);
		    }
		    else  {
		       sit->link[dir].e[j][k] = cmplx(1.0,0.0);
		    }
		}
	    }
	}
    }

    node0_printf("unit gauge configuration loaded\n");
}

void funnylat()  {
    /* sets link matrices to funny matrices for debugging */
    register int i,j,k,dir;
    register site *sit;

    FORALLSITES(i,sit){
	for(dir=XUP;dir<=TUP;dir++){
	    for(j=0; j<3; ++j)  {
		for(k=0; k<3; ++k)  {
		    sit->link[dir].e[j][k] = cmplx(0.0,0.0);
		}
	    }
	    sit->link[dir].e[0][0].real = dir;
	    sit->link[dir].e[1][1].real = 10*sit->x;
	    sit->link[dir].e[2][2].real = 100*sit->y;
	    sit->link[dir].e[0][0].imag = dir;
	    sit->link[dir].e[1][1].imag = 10*sit->z;
	    sit->link[dir].e[2][2].imag = 100*sit->t;
	}
    }
}


/* get_f is used to get a Realing point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string precede the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */

int get_f( int prompt, char *variable_name_string, Real *value ){
    int s;
    char checkname[80];

    if(prompt)  {
    	printf("enter %s ",variable_name_string);
    	s=scanf("%leHELP",value);
    	if(s == 1){
	    printf("%s %g\n",variable_name_string,*value);
	    return(0);
	}
    }
    else  {
    	s=scanf("%s%leHELP",checkname,value);
    	if (s == EOF) return(1);
    	if(s == 2 && strcmp(checkname,variable_name_string) == 0){
	    printf("%s %g\n",variable_name_string,*value);
    	    return(0);
	}
    }
    printf("get_f: ERROR IN INPUT: %s\n",variable_name_string);
    return(1);
}

int get_i( int prompt, char *variable_name_string, int *value ){
    int s;
    char checkname[80];

    if(prompt)  {
    	printf("enter %s ",variable_name_string);
    	s=scanf("%d",value);
    	if (s == 1){
	    printf("%s %d\n",variable_name_string,*value);
	    return(0);
	}
    }
    else  {
    	s=scanf("%s%d",checkname,value);
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

int get_s( int prompt, char *variable_name_string, char *value ){
    int s;
    char checkname[80];

    if(prompt)  {
    	printf("enter %s ",variable_name_string);
    	s=scanf("%s",value);
    	if(s == 1){
	    printf("%s %s\n",variable_name_string,value);
	    return(0);
	}
    }
    else  {
    	s=scanf("%s %s",checkname,value);
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
int get_prompt( int *prompt ){
    char initial_prompt[80];

    *prompt = -1;
    printf( "type 0 for no prompts  or 1 for prompts\n");
    scanf("%s",initial_prompt);
    if(strcmp(initial_prompt,"prompt") == 0)  {
       scanf("%d",prompt);
    }
    else if(strcmp(initial_prompt,"0") == 0) *prompt=0;
    else if(strcmp(initial_prompt,"1") == 0) *prompt=1;

    if( *prompt==0 || *prompt==1 )return(0);
    else{
        printf("get_prompt: ERROR IN INPUT: initial prompt\n");
        return(1);
    }
}
