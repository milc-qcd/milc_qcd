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

#ifndef NO_GAUGE_FIELD

/* save a lattice in any of the formats:
    SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_CHECKPOINT
*/
gauge_file *save_lattice( int flag, char *filename, char *stringLFN){
    double dtime;
    gauge_file *gf = NULL;

#ifndef NO_GAUGE_FIELD
    d_plaquette(&g_ssplaq,&g_stplaq);
    d_linktrsum(&linktrsum);
    nersc_checksum = nersc_cksum();
#endif

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
        case SAVE_PARTFILE_ILDG:
#ifdef HAVE_QIO
 	    gf = save_partfile_ildg(filename,stringLFN);
#else
	    node0_printf("save_partfile_ildg requires QIO compilation\n");
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
        case SAVE_PARTFILE_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_partfile_scidac(filename);
#else
	    node0_printf("save_partfile_scidac requires QIO compilation\n");
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
#if (PRECISION==1)
    node0_printf("CHECK PLAQ: %e %e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %e CKSUM: %x\n",
		 linktrsum.real/3.,nersc_checksum);
#else
    /* Double precision */
    node0_printf("CHECK PLAQ: %.16e %.16e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
		 linktrsum.real/3.,nersc_checksum);
#endif
    return gf;
}

/* reload a lattice in any of the formats, or cold lattice, or keep
     current lattice:
    FRESH, CONTINUE,
    RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL
*/
void coldlat(void);

gauge_file *reload_lattice( int flag, char *filename){
    double dtime;
    gauge_file *gf = NULL;
    Real max_deviation;
#if PRECISION == 2
    Real max_deviation2;
#endif

    dtime = -dclock();
    switch(flag){
	case CONTINUE:	/* return NULL.  We lose information if we do this  */
 	    node0_printf("reload_lattice: WARNING: gaugefile struct set to NULL\n");
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
#ifndef NO_GAUGE_FIELD
    d_plaquette(&g_ssplaq,&g_stplaq);
    d_linktrsum(&linktrsum);
    nersc_checksum = nersc_cksum();
#endif
    if(this_node==0){
#if (PRECISION==1)
    node0_printf("CHECK PLAQ: %e %e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %e CKSUM: %x\n",
		 linktrsum.real/3.,nersc_checksum);
    fflush(stdout);
#else
    /* Double precision */
    node0_printf("CHECK PLAQ: %.16e %.16e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
		 linktrsum.real/3.,nersc_checksum);
    fflush(stdout);
#endif
    }
    dtime = -dclock();
    max_deviation = check_unitarity();
    g_floatmax(&max_deviation);
#if (PRECISION==1)
    if(this_node==0)printf("Unitarity checked.  Max deviation %.2e\n",
			   max_deviation); fflush(stdout);
#else
    reunitarize();
    max_deviation2 = check_unitarity();
    g_floatmax(&max_deviation2);
    if(this_node==0)
      printf("Reunitarized for double precision. Max deviation %.2e changed to %.2e\n",
                       max_deviation,max_deviation2); fflush(stdout);
#endif
    dtime += dclock();
    if(this_node==0)printf("Time to check unitarity = %e\n",dtime);
    return gf;
}

#endif /* ifndef NO_GAUGE_FIELD */

/* Get next tag, but skip past end of line if we encounter # for comments */
#define MAX_TAG 512
char *get_next_tag(FILE *fp, char *tag, char *myname){
  static char line[MAX_TAG];
  int s;

  while(1){
    s = fscanf(fp, "%s",line);
    if (s == EOF){
      printf("%s(%d): EOF on input.\n",myname,this_node);
      return NULL;
    }
    if(s == 0){
      printf("%s(%d): Error reading %s\n",myname,this_node,tag);
      return NULL;
    }
    /* Provide for comment lines with # before "prompt" */
    if(strchr(line,'#')==NULL)break;
    else{
      printf("%s",line);  /* Print string with # character*/
      if(fgets(line,MAX_TAG,fp)==NULL){
	printf("%s(%d) EOF on input.\n",myname,this_node);
	return NULL;
      }
      printf("%s",line); /* Print rest of line */
    }
  }
  return line;
}

/* Read and echo the next tag.  Echo any intervening comments */
/* Comments begin with # and apply to the rest of the line */
/* Verify that the input tag agrees with the expected tag */

int get_check_tag(FILE *fp, char *tag, char *myname){
  char *checktag;
  
  checktag = get_next_tag(fp, tag, myname);
  if(checktag == NULL)return 1;
  
  if(strcmp(checktag,tag) != 0){
    printf("\n%s: ERROR IN INPUT: expected %s but found %s\n",
	   myname,tag,checktag);
    return 1;
  }
  printf("%s ",tag);
  return 0;
}

/* Check return value of scanf */
static int check_read(int s, char *myname, char *tag){

  if (s == EOF){
    printf("\n%s: Expecting value for %s but found EOF.\n",
	   myname,tag);
    return 1;
  }
  else if(s == 0){
    printf("\n%s: Format error reading value for %s\n",
	   myname,tag);
    return 1;
  }
  else
    return 0;
}


/* find out what kind of starting lattice to use, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_starting_lattice( FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_starting_lattice";
  
  if (prompt==1) printf(
			"enter 'continue', 'fresh', 'reload_ascii', 'reload_serial', or 'reload_parallel'\n");
  
  savebuf = get_next_tag(fp, "read lattice command", myname);
  if (savebuf == NULL)return 1;
  
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
    return 1;
  }
  
  /*read name of file and load it */
  if( *flag != FRESH && *flag != CONTINUE ){
    if(prompt==1)printf("enter name of file containing lattice\n");
    status=fscanf(fp," %s",filename);
    if(status !=1) {
      printf("\n%s(%d): ERROR IN INPUT: error reading file name\n",
	     myname, this_node); 
      return 1;
    }
    printf("%s\n",filename);
  }
  return 0;
}

/* find out what do to with lattice at end, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_ending_lattice";
  
  if (prompt==1) printf(
			"'forget' lattice at end,  'save_ascii', 'save_serial', 'save_parallel', 'save_checkpoint', 'save_serial_fm', 'save_serial_scidac', 'save_parallel_scidac', 'save_multifile_scidac', 'save_partfile_scidac', 'save_serial_archive', 'save_serial_ildg', 'save_parallel_ildg', 'save_partfile_ildg', or 'save_multifile_ildg'\n");
  
  savebuf = get_next_tag(fp, "save lattice command", myname);
  if (savebuf == NULL)return 1;

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
  else if(strcmp("save_partfile_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
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
  else if(strcmp("save_partfile_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_ILDG;
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
    return 1;
  }
  
  if( *flag != FORGET ){
    if(prompt==1)printf("enter filename\n");
    status=fscanf(fp,"%s",filename);
    if(status !=1){
      printf("\nask_ending_lattice: ERROR IN INPUT: error reading filename\n"); return 1;
    }
    printf("%s\n",filename);
    
  }
  return 0;
}

/*--------------------------------------------------------------------*/

/* For FNAL formatted ASCII correlator files */

int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename){

  char *savebuf;
  int status;
  char myname[] = "ask_corr_file";

  if (prompt==1)
    printf("'forget_corr', 'save_corr_fnal' for correlator file type\n");

  savebuf = get_next_tag(fp, "output correlator file command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("forget_corr",savebuf) == 0 ){
    *flag = FORGET;
  }
  else if(strcmp("save_corr_fnal",savebuf) == 0 ) {
    *flag = SAVE_ASCII;  /* Lazy borrowing of lattice flags! */
  }
  else{
    printf("is not a save correlator command. INPUT ERROR\n");
    return 1;
  }
  
  if( *flag == FORGET )
    printf("\n");
  else{
    if(prompt==1)printf("enter filename\n");
    status=fscanf(fp,"%s",filename);
    if(status !=1){
      printf("\n%s: ERROR reading filename\n",myname); 
      return 1;
    }
    printf("%s\n",filename);
  }

  return 0;

} /* ask_corr_file */

/*--------------------------------------------------------------------*/

int ask_ildg_LFN(FILE *fp, int prompt, int flag, char *stringLFN){
  int status = 0;

  /* For ILDG output formats we require a logical file name next */
  if( flag == SAVE_SERIAL_ILDG ||
      flag == SAVE_PARALLEL_ILDG ||
      flag == SAVE_PARTFILE_ILDG ||
      flag == SAVE_MULTIFILE_ILDG ){
    status = get_s(fp, prompt, "ILDG_LFN", stringLFN);
  }
  else
    stringLFN[0] = '\0';
  return status;
}

#ifndef NO_GAUGE_FIELD

void coldlat(void){
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

void funnylat(void)  {
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

#endif

/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string precede the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */

int get_f( FILE *fp, int prompt, char *tag, Real *value ){
    int s;
    char checkvalue[80];
    char myname[] = "get_f";

    if(prompt==1)  {
	s = 0;
	while(s != 1){
	  printf("enter %s ",tag);
	  fscanf(fp,"%s",checkvalue);
#if PRECISION == 1
	  s=sscanf(checkvalue,"%e",value);
#else
	  s=sscanf(checkvalue,"%le",value);
#endif
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  else printf("%s %g\n",tag,*value);
	}
    }
    else  {
      if(get_check_tag(fp, tag, myname) == 1)return 1;
	  
#if PRECISION == 1
      s = fscanf(fp,"%e",value);
#else
      s = fscanf(fp,"%le",value);
#endif
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%g\n",*value);
    }

    return 0;
}

int get_i( FILE *fp, int prompt, char *tag, int *value ){
    int s;
    char checkvalue[80];
    char myname[] = "get_i";

    if(prompt==1)  {
      s = 0;
      while(s != 1){
    	printf("enter %s ",tag);
	fscanf(fp,"%s",checkvalue);
    	s=sscanf(checkvalue,"%d",value);
	if(s==EOF)return 1;
	if(s==0)printf("Data format error.\n");
	else printf("%s %d\n",tag,*value);
      }
    }
    else  {
      if(get_check_tag(fp, tag, myname) == 1)return 1;
	  
      s = fscanf(fp,"%d",value);
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%d\n",*value);
    }
    
    return 0;

}

/* Read a single word as a string without printing an end-of-line  */

int get_sn( FILE *fp, int prompt, char *tag, char *value ){
    int s;
    char myname[] = "get_sn";

    if(prompt==1)  {
      s = 0;
      while(s != 1){
    	printf("enter %s ",tag);
    	s=fscanf(fp,"%s",value);
	if(s==EOF)return 1;
	if(s==0)printf("Data format error.\n");
	else printf("%s %s",tag,value);
      }
    }
    else  {
      if(get_check_tag(fp, tag, myname) == 1)return 1;

      s = fscanf(fp,"%s",value);
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%s",value);
    }
    return 0;
}

/* Read a single word as a string */

int get_s( FILE *fp, int prompt, char *tag, char *value ){
    int s;

    s = get_sn(fp, prompt, tag, value);
    printf("\n");
    return s;
}

/* Read a vector of integers */
int get_vi( FILE* fp, int prompt, char *tag, 
	    int *value, int nvalues ){
    int s,i;
    char myname[] = "get_vi";

    if(prompt==1)  {
      printf("enter %s with %d values",tag, nvalues);
      for(i = 0; i < nvalues; i++){
	s = 0;
	while(s != 1){
	  printf("\n[%d] ",i);
	  s=fscanf(fp,"%d",value+i);
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  printf("%s %d\n",tag,value[i]);
	}
      }
    }
    else  {
      if(get_check_tag(fp, tag, myname) == 1)return 1;
	  
      for(i = 0; i < nvalues; i++){
	s = fscanf(fp,"%d",value + i);
	if(check_read(s,myname,tag) == 1)return 1;
	printf("%d ",value[i]);
      }
      printf("\n");
    }
    
    return 0;

}

/* Read a vector of reals */
int get_vf( FILE* fp, int prompt, char *tag, 
	    Real *value, int nvalues ){
    int s,i;
    char myname[] = "get_vf";

    if(prompt==1)  {
      printf("enter %s with %d values",tag, nvalues);
      for(i = 0; i < nvalues; i++){
	s = 0;
	while(s != 1){
	  printf("\n[%d] ",i);
#if PRECISION == 1
	  s=scanf("%e",value+i);
#else
	  s=scanf("%le",value+i);
#endif
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  printf("%s %g\n",tag,*(value+i));
	}
      }
    }
    else  {
      if(get_check_tag(fp, tag, myname) == 1)return 1;
	  
      for(i = 0; i < nvalues; i++){
#if PRECISION == 1
	s = fscanf(fp,"%e",value + i);
#else
	s = fscanf(fp,"%le",value + i);
#endif
	if(check_read(s,myname,tag) == 1)return 1;
	printf("%g ",value[i]);
      }
      printf("\n");
    }
    
    return 0;
}

/* Read a vector of strings */

int get_vs( FILE *fp, int prompt, char *tag, char **value, int nvalues ){
  char myname[] = "get_vs";

  int s, i;
  
  if(prompt==1)  {
    s = 0;
    printf("enter %s with %d values", tag, nvalues);
      for(i = 0; i < nvalues; i++){
	while(s != 1){
	  printf("\n[%d] ",i);
	  s=scanf("%s",value[i]);
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  printf("%s %s\n",tag,value[i]);
	}
      }
    } else  {
    if(get_check_tag(fp, tag, myname) == 1)return 1;
	  
    for(i = 0; i < nvalues; i++){
      s = fscanf(fp,"%s",value[i]);
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%s ",value[i]);
    }
    printf("\n");
  }
  
  return 0;
}
    
/* get_prompt gets the initial value of prompt */
/* 0 for reading from file, 1 prompts for input from terminal,
   2 for checking input parameters without execution (for proofreading) */
/* should be called only by node 0 */
/* return 0 if sucessful, 1 if failure */
int get_prompt( FILE *fp, int *prompt ){
    char initial_prompt[512];
    int status;
    char myname[] = "get_prompt";

    *prompt = -1;
    printf( "type 0 for no prompts, 1 for prompts, or 2 for proofreading\n");
    while(1){
      status = fscanf(fp, "%s",initial_prompt);
      if(status != 1){
	printf("\n%s: Can't read input\n",myname);
	terminate(1);
      }
      if(strchr(initial_prompt,'#')==NULL)break;
      /* Provide for comment lines with # before "prompt" */
      else{
	printf("%s",initial_prompt);
	if(fgets(initial_prompt,512,fp)==NULL){
	  printf("%s(%d) EOF on input.\n",myname,this_node);
	  return 1;
	}
	printf("%s",initial_prompt);
      }
    }
    if(strcmp(initial_prompt,"prompt") == 0)  {
      fscanf(fp, "%d",prompt);
    }
    else if(strcmp(initial_prompt,"0") == 0) *prompt=0;
    else if(strcmp(initial_prompt,"1") == 0) *prompt=1;
    else if(strcmp(initial_prompt,"2") == 0) *prompt=2;

    if( *prompt==0 || *prompt==1 || *prompt==2 )return 0;
    else{
        printf("\n%s: ERROR IN INPUT: initial prompt\n",myname);
        return 1;
    }
}
