#ifndef _PVM_H
#define _PVM_H
/* Prototypes for pvm library calls */
/* PVM Version 2.4 */

#ifdef PROTO
int barrier(char *barrier_name, int num);
int enroll(char *component_name);
int getbytes(char *ptr, int num);
int getncplx(Real *ptr, int num);
int getndcplx(double *ptr, int num);
int getndReal(double *ptr, int num);
int getnReal(Real *ptr, int num);
int getnint(int *ptr, int num);
int getstring(char *ptr);
int initiate(char *object_file, char *arch);
void initsend();
void leave();
int probe(int msgtype);
int probemulti(int num, int *msgtypes);
int pstatus(int *ncpu, int *nformat);
int putbytes(char *ptr, int num);
int putncplx(Real *ptr, int num);
int putnReal(Real *ptr, int num);
int putndcplx(double *ptr, int num);
int putndReal(double *ptr, int num);
int putnint(int *ptr, int num);
int putstring(char *ptr);
int rcv(int msgtype);
int rcvinfo(int *bytes, int *msgtype, char *component, int *instance);
int snd(char *component, int instance, int msgtype);
int whoami(char *component, int *instance);

#else	/* sun 4 compiler doesn't speak ANSI */
int barrier();
int enroll();
int getbytes();
int getncplx();
int getndcplx();
int getndReal();
int getnReal();
int getnint();
int getstring();
int initiate();
void initsend();
void leave();
int probe();
int probemulti();
int pstatus();
int putbytes();
int putncplx();
int putnReal();
int putndcplx();
int putndReal();
int putnint();
int putstring();
int rcv();
int rcvinfo();
int snd();
int whoami();

#endif /* end if defined PROTO */
#endif /* _PVM_H */
