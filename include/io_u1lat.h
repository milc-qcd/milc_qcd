/* ************************************************************ */
/*								*/
/*				IO_U1LAT.H			*/
/*								*/
/* Prototypes for gauge configuration I/O			*/
/*								*/
/* Last Updated on 07.24.07                                     */
/* Last Updated on 01.05.12 by S. Gottlieb                      */
/*								*/
/* ************************************************************ */
#ifndef _IO_U1LAT_H
#define _IO_U1LAT_H

#include "../include/int32type.h"
#include "../include/io_lat.h"

#define U1GAUGE_VERSION_NUMBER 0x7d7
#define U1GAUGE_VERSION_NUMBER_BINARY 0x8d7

/* io_u1lat.c */
int ask_ending_u1_lattice(FILE *fp,
    int prompt,int *flag,char *filename);
gauge_file *save_u1_lattice(int flag,char *filename);
gauge_file *save_u1_ascii(char *filename);
gauge_file *save_u1_serial(char *filename);
gauge_file *save_u1_parallel(char *filename);
gauge_file *save_u1_serial_scidac(char *filename);
gauge_file *save_u1_parallel_scidac(char *filename);
gauge_file *setup_output_u1gauge_file(void);
void write_u1gauge_info_file(gauge_file *gf);

void write_appl_u1gauge_info(FILE *fp, gauge_file *gf);

int ask_starting_u1_lattice(FILE *fp,
    int prompt,int *flag,char *filename);
gauge_file *reload_u1_lattice(int flag, char *filename);
gauge_file *restore_u1_ascii(char *filename);
gauge_file *restore_u1_serial(char *filename);
gauge_file *restore_u1_parallel(char *filename);
gauge_file *setup_input_u1gauge_file(char *filename);

void f2d_4cmplx(fcomplex *a,complex *b);
void d2f_4cmplx(complex *a,fcomplex *b);
void cold_u1lat(void);
gauge_file *setup_u1_output_gauge_file(void);
gauge_file *r_u1_serial_i(char *filename);
int read_u1_gauge_hdr(gauge_file *gf, int parallel);

#endif /* _IO_U1LAT_H */

/* ************************************************************ */

