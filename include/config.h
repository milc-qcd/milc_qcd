#ifndef _CONFIG_H
#define _CONFIG_H
/* config.h.  For now, NOT generated automatically by configure.  */

/* Collects macros for preprocessor tweaks that accommodate
   differences in compilers, architecture and OS */

/********************************************************************/
/* Code version                                                     */
/********************************************************************/

#define MILC_CODE_VERSION "7.7.8"

/********************************************************************/
/* Compiler/Processor-dependent macros */
/********************************************************************/

/* Specify the unsigned 32 bit integer base type for this compiler */
/* Run the script "getint.sh" to find out what to use */
/* One and only one of these should be defined */
#define INT_IS_32BIT 1  /* Most present systems */
#undef SHORT_IS_32BIT   /* Needed on T3E UNICOS, for example */

/* Define if the target processor has native double precision */
/* (For some library routines, gives slightly better performance) */
/* Systems that do: IBM SP.  Most do not. */
#undef NATIVEDOUBLE

/* Define if the cache line is 64 bytes (if not, we assume 32 bytes). */
/* Processors that do: P4 (actually fetches 128), EV67, EV68  */
/* Used only for prefetching, so it only affects performance */
#undef HAVE_64_BYTE_CACHELINE

/********************************************************************/
/* Compiler/OS-dependent macros */
/********************************************************************/

/* Define if you have the <ieeefp.h> header file. */
/* Systems that don't: T3E UNICOS, Exemplar, Linux gcc, SP AIX, HP/Compaq True64 */
#undef HAVE_IEEEFP_H

/* Define if you have the <unistd.h> header file. */
/* Systems that don't: NT */
#define HAVE_UNISTD_H 1

/* Define if you have the <sys/time.h> header file. */
/* Most systems do */
#define HAVE_SYS_TIME_H 1

/* Define if you have ANSI "fseeko" */
/* #undef HAVE_FSEEKO */  
/* Systems that don't: T3E UNICOS */
#define HAVE_FSEEKO 1

#endif /* _CONFIG_H */
