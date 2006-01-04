/********************** io_helpers.c **********************************/
/* MIMD version 7 */
/* DT 8/97 
     General purpose high level routines, to be used by any application
     that wants them.
*/

#include "generic_includes.h"
#include "../include/io_lat.h"
#include "../include/file_types.h"
#ifdef HAVE_QIO
#include <qio.h>
#endif

static file_type gauge_list[N_GAUGE_TYPES] =
  { {FILE_TYPE_GAUGE_V1,      GAUGE_VERSION_NUMBER_V1},
    {FILE_TYPE_GAUGE_V5,      GAUGE_VERSION_NUMBER},
    {FILE_TYPE_GAUGE_1996,    GAUGE_VERSION_NUMBER_1996},
    {FILE_TYPE_GAUGE_FNAL,    GAUGE_VERSION_NUMBER_FNAL},
    {FILE_TYPE_GAUGE_ARCHIVE, GAUGE_VERSION_NUMBER_ARCHIVE},
    {FILE_TYPE_GAUGE_SCIDAC,  LIME_MAGIC_NO}
  };

/* save a lattice in any of the formats:
    SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_CHECKPOINT
*/
gauge_file *save_lattice( int flag, char *filename, char *stringLFN){
    double dtime;
    gauge_file *gf;
    double ssplaq, stplaq;
    int do_scidac = 0;

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
	case SAVE_PARALLEL:
	    gf = save_parallel(filename);
	    break;
	case SAVE_CHECKPOINT:
	    gf = save_checkpoint(filename);
	    break;
        case SAVE_SERIAL_FM:
 	    printf("Save serial FNAL format not implemented\n");
	    break;
        case SAVE_SERIAL_ILDG:
#ifdef HAVE_QIO
 	    gf = save_serial_ildg(filename,stringLFN);
#else
	    node0_printf("save_serial_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARALLEL_ILDG:
#ifdef HAVE_QIO
 	    gf = save_parallel_ildg(filename,stringLFN);
#else
	    node0_printf("save_parallel_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARTITION_ILDG:
#ifdef HAVE_QIO
 	    gf = save_partition_ildg(filename,stringLFN);
#else
	    node0_printf("save_partition_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_MULTIFILE_ILDG:
#ifdef HAVE_QIO
 	    gf = save_multifile_ildg(filename,stringLFN);
#else
	    node0_printf("save_multifile_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_SERIAL_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_serial_scidac(filename);
#else
	    node0_printf("save_serial_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARALLEL_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_parallel_scidac(filename);
#else
	    node0_printf("save_parallel_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
	    gf = save_multifile_scidac(filename);
#else
	    node0_printf("save_multifile_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARTITION_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_partition_scidac(filename);
#else
	    node0_printf("save_partition_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
	case SAVE_SERIAL_ARCHIVE:
	    gf = save_serial_archive(filename);
	    break;
	default:
	    node0_printf("\nsave_lattice: ERROR: unknown type for saving lattice\n");
	    terminate(1);
    }
    dtime += dclock();
    if(flag != FORGET)
      node0_printf("Time to save = %e\n",dtime);
#ifndef NOLINKS
    d_plaquette(&ssplaq,&stplaq);
#endif
#if (PRECISION==1)
    node0_printf("CHECK PLAQ: %e %e\n",ssplaq,stplaq);
#else
    /* Double precision */
    node0_printf("CHECK PLAQ: %.16e %.16e\n",ssplaq,stplaq);
#endif
    return gf;
}

/* reload a lattice in any of the formats, or cold lattice, or keep
     current lattice:
    FRESH, CONTINUE,
    RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL
*/
void coldlat();

gauge_file *reload_lattice( int flag, char *filename){
    double dtime;
    gauge_file *gf;
    double ssplaq, stplaq;
    Real max_deviation;

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
#ifdef SCHROED_FUN
    set_boundary_fields();
#endif
#ifndef NOLINKS
    d_plaquette(&ssplaq,&stplaq);
#endif
    if(this_node==0){
#if (PRECISION==1)
        printf("CHECK PLAQ: %e %e\n",ssplaq,stplaq);fflush(stdout);
#else
	/* Double precision */
	node0_printf("CHECK PLAQ: %.16e %.16e\n",ssplaq,stplaq);fflush(stdout);
#endif
    }
    dtime = -dclock();
    max_deviation = check_unitarity();
    g_floatmax(&max_deviation);
    dtime += dclock();
    if(this_node==0)printf("Unitarity checked.  Max deviation %.2e\n",
			   max_deviation); fflush(stdout);
    if(this_node==0)printf("Time to check unitarity = %e\n",dtime);
#if (PRECISION==2)
    reunitarize();
    dtime = -dclock();
    max_deviation = check_unitarity();
    g_floatmax(&max_deviation);
    dtime += dclock();
    if(this_node==0)printf("Unitarity checked.  Max deviation %.2e\n",
                       max_deviation); fflush(stdout);
    if(this_node==0)printf("Time to check unitarity = %e\n",dtime);
#endif
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
    if (status == EOF){
      printf("ask_starting_lattice: EOF on STDIN.\n");
      return(1);
    }
    if(status !=1) {
        printf("\nask_starting_lattice: ERROR IN INPUT: can't read starting lattice command\n");
        return(1);
    }

    printf("%s ",savebuf);
    if(strcmp("fresh",savebuf) == 0 ){
       *flag = FRESH;
    printf("\n");
    }
    else if(strcmp("continue",savebuf) == 0 ) {
        *flag = CONTINUE;
	printf("\n");
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
    	printf(" is not a valid starting lattice command. INPUT ERROR.\n"); 
	return(1);
    }

    /*read name of file and load it */
    if( *flag != FRESH && *flag != CONTINUE ){
        if(prompt!=0)printf("enter name of file containing lattice\n");
        status=scanf("%s",filename);
        if(status !=1) {
	    printf("\nask_starting_lattice: ERROR IN INPUT: error reading file name\n"); return(1);
        }
	printf("%s\n",filename);
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
        "'forget' lattice at end,  'save_ascii', 'save_serial', 'save_parallel', 'save_checkpoint', 'save_serial_fm', 'save_serial_scidac', 'save_parallel_scidac', 'save_multifile_scidac', 'save_partition_scidac', 'save_serial_archive', 'save_serial_ildg', 'save_parallel_ildg', 'save_partition_ildg', or 'save_multifile_ildg'\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("\nask_ending_lattice: ERROR IN INPUT: error reading ending lattice command\n");
        return(1);
    }
    printf("%s ",savebuf);
    if(strcmp("save_ascii",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_serial",savebuf) == 0 ) {
        *flag=SAVE_SERIAL;
    }
    else if(strcmp("save_parallel",savebuf) == 0 ) {
      *flag=SAVE_PARALLEL;
    }
    else if(strcmp("save_checkpoint",savebuf) == 0 ) {
        *flag=SAVE_CHECKPOINT;
    }
    else if(strcmp("save_serial_fm",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_FM;
    }
    else if(strcmp("save_serial_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_SERIAL_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_parallel_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARALLEL_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_multifile_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_MULTIFILE_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_partition_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARTITION_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_serial_archive",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_ARCHIVE;
    }
    else if(strcmp("save_serial_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_SERIAL_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_parallel_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARALLEL_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_partition_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARTITION_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_multifile_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_MULTIFILE_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("forget",savebuf) == 0 ) {
        *flag=FORGET;
	printf("\n");
    }
    else {
      printf("is not a save lattice command. INPUT ERROR\n");
      return(1);
    }

    if( *flag != FORGET ){
        if(prompt!=0)printf("enter filename\n");
        status=scanf("%s",filename);
        if(status !=1){
    	    printf("\nask_ending_lattice: ERROR IN INPUT: error reading filename\n"); return(1);
        }
	printf("%s\n",filename);

    }
    return 0;
}

int ask_ildg_LFN(int prompt, int flag, char *stringLFN){
  int status = 0;

  /* For ILDG output formats we require a logical file name next */
  if( flag == SAVE_SERIAL_ILDG ||
      flag == SAVE_PARALLEL_ILDG ||
      flag == SAVE_PARTITION_ILDG ||
      flag == SAVE_MULTIFILE_ILDG ){
    status = get_s(prompt, "ILDG_LFN", stringLFN);
  }
  else
    stringLFN[0] = '\0';
  return status;
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


/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string precede the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */

int get_f( int prompt, char *variable_name_string, Real *value ){
    int s;
    char checkname[80];

    if(prompt)  {
	s = 0;
	while(s != 1){
	  printf("enter %s ",variable_name_string);
	  scanf("%s",checkname);
#if PRECISION == 1
	  s=sscanf(checkname,"%e",value);
#else
	  s=sscanf(checkname,"%le",value);
#endif
	  if(s == 1)
	    printf("%s %g\n",variable_name_string,*value);
	  else
	    printf("Data format error.\n");
	}
    }
    else  {
      s = scanf("%s",checkname);
      if (s == EOF){
	printf("get_f: EOF on STDIN.\n");
	return(1);
      }
      else if(s==0){
	printf("get_f: Format error looking for %s\n",variable_name_string);
	return(1);
      }
      else if(strcmp(checkname,variable_name_string) != 0){
	printf("\nget_f: ERROR IN INPUT: expected %s but found %s\n",
	       variable_name_string,checkname);
	return(1);
      }

      printf("%s ",variable_name_string);
	  
#if PRECISION == 1
      s = scanf("%e",value);
#else
      s = scanf("%le",value);
#endif
      if (s == EOF){
	printf("\nget_f: Expecting value for %s but found EOF.\n",
	       variable_name_string);
	return(1);
      }
      else if(s==0){
	printf("\nget_f: Format error reading value for %s\n",
	       variable_name_string);
	return(1);
      }
      printf("%g\n",*value);
    }

    return(0);
}

int get_i( int prompt, char *variable_name_string, int *value ){
    int s;
    char checkname[80];

    if(prompt)  {
      s = 0;
      while(s != 1){
    	printf("enter %s ",variable_name_string);
	scanf("%s",checkname);
    	s=sscanf(checkname,"%d",value);
    	if (s == 1)
	  printf("%s %d\n",variable_name_string,*value);
	else
	  printf("Data format error.\n");
      }
    }
    else  {
      s = scanf("%s",checkname);
      if (s == EOF){
	printf("get_i: EOF on STDIN.\n");
	return(1);
      }
      else if(s==0){
	printf("\nget_i: Format error looking for %s\n",variable_name_string);
	return(1);
      }
      else if(strcmp(checkname,variable_name_string) != 0){
	printf("\nget_i: ERROR IN INPUT: expected %s but found %s\n",
	       variable_name_string,checkname);
	return(1);
      }

      printf("%s ",variable_name_string);
	  
      s = scanf("%d",value);
      if (s == EOF){
	printf("\nget_i: Expecting value for %s but found EOF.\n",
	       variable_name_string);
	return(1);
      }
      else if(s==0){
	printf("\nget_i: Format error reading value for %s\n",
	       variable_name_string);
	return(1);
      }
      printf("%d\n",*value);
    }
    
    return(0);

}

/* Read a single word as a string */

int get_s( int prompt, char *variable_name_string, char *value ){
    int s;
    char checkname[80];

    if(prompt)  {
      s = 0;
      while(s != 1){
    	printf("enter %s ",variable_name_string);
    	s=scanf("%s",value);
    	if(s == 1)
	  printf("%s %s\n",variable_name_string,value);
	else
	  printf("Data format error.\n");
	}
    }
    else  {
      s = scanf("%s",checkname);
      if (s == EOF){
	printf("get_s: EOF on STDIN.\n");
	return(1);
      }
      else if(s==0){
	printf("\nget_s: Format error looking for %s\n",variable_name_string);
	return(1);
      }
      else if(strcmp(checkname,variable_name_string) != 0){
	printf("\nget_s: ERROR IN INPUT: expected %s but found %s\n",
	       variable_name_string,checkname);
	return(1);
      }

      printf("%s ",variable_name_string);
	  
      s = scanf("%s",value);
      if (s == EOF){
	printf("\nget_s: Expecting value for %s but found EOF.\n",
	       variable_name_string);
	return(1);
      }
      else if(s==0){
	printf("\nget_s: Format error reading value for %s\n",
	       variable_name_string);
	return(1);
      }
      printf("%s\n",value);
    }
    return(0);
}

/* Read a vector of integers */
int get_vi( int prompt, char *variable_name_string, int *value, int nvalues ){
    int s;
    int i;
    char checkname[80];

    if(prompt)  {
      s = 0;
      while(s != 1){
    	printf("enter %s with %d values",variable_name_string, nvalues);
	scanf("%s",checkname);
    	s=sscanf(checkname,"%d",value);
    	if (s == 1)
	  printf("%s %d\n",variable_name_string,*value);
	else
	  printf("Data format error.\n");
      }
    }
    else  {
      s = scanf("%s",checkname);
      if (s == EOF){
	printf("get_vi: EOF on STDIN.\n");
	return(1);
      }
      else if(s==0){
	printf("\nget_vi: Format error looking for %s\n",variable_name_string);
	return(1);
      }
      else if(strcmp(checkname,variable_name_string) != 0){
	printf("\nget_vi: ERROR IN INPUT: expected %s but found %s\n",
	       variable_name_string,checkname);
	return(1);
      }

      printf("%s ",variable_name_string);
	  
      for(i = 0; i < nvalues; i++){
	
	s = scanf("%d",value + i);
	if (s == EOF){
	  printf("\nget_vi: Expecting value for %s but found EOF.\n",
		 variable_name_string);
	  return(1);
	}
	else if(s==0){
	  printf("\nget_vi: Format error reading value for %s\n",
		 variable_name_string);
	  return(1);
	}
	printf("%d ",value[i]);
      }
      printf("\n");
    }
    
    return(0);

}

/* get_prompt gets the initial value of prompt */
/* 0 for reading from file, 1 prompts for input from terminal */
/* should be called only by node 0 */
/* return 0 if sucessful, 1 if failure */
int get_prompt( int *prompt ){
    char initial_prompt[80];
    int status;

    *prompt = -1;
    printf( "type 0 for no prompts  or 1 for prompts\n");
    status = scanf("%s",initial_prompt);
    if(status != 1){
      printf("\nget_prompt: Can't read stdin\n");
      terminate(1);
    }
    if(strcmp(initial_prompt,"prompt") == 0)  {
       scanf("%d",prompt);
    }
    else if(strcmp(initial_prompt,"0") == 0) *prompt=0;
    else if(strcmp(initial_prompt,"1") == 0) *prompt=1;

    if( *prompt==0 || *prompt==1 )return(0);
    else{
        printf("\nget_prompt: ERROR IN INPUT: initial prompt\n");
        return(1);
    }
}
